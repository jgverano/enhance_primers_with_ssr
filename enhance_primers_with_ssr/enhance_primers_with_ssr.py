#!/usr/bin/env python3
"""
Enhance primer TSV with SSR metadata from MISA files.

Integrates Simple Sequence Repeat (SSR) information from MISA analysis
into primer design output from primers_to_tsv.py.
"""

import argparse
import sys
import re
import pandas as pd
from typing import Dict, List, Optional


def parse_misa_file(misa_file: str) -> Dict[str, List[Dict]]:
    """
    Parse MISA file to extract SSR metadata.
    
    Args:
        misa_file: Path to MISA file (*.misa format)
        
    Returns:
        Dictionary mapping sequence IDs to list of SSR metadata dicts:
        {sequence_id: [{ssr_nr: 1, type: 'p3', motif: '(ATC)6', ...}, ...]}
    """
    ssr_data = {}
    
    with open(misa_file, 'r') as f:
        header = f.readline().strip().split('\t')
        
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 7:
                continue
            
            seq_id = parts[0]
            ssr_nr = int(parts[1])
            ssr_type = parts[2]
            ssr_motif = parts[3]
            ssr_size = int(parts[4])
            ssr_start = int(parts[5])
            ssr_end = int(parts[6])
            
            if seq_id not in ssr_data:
                ssr_data[seq_id] = []
            
            ssr_data[seq_id].append({
                'ssr_nr': ssr_nr,
                'ssr_type': ssr_type,
                'ssr_motif': ssr_motif,
                'ssr_size': ssr_size,
                'ssr_start': ssr_start,
                'ssr_end': ssr_end
            })
    
    return ssr_data


def parse_primer3_output(p3out_file: str) -> Dict[str, Dict]:
    """
    Parse Primer3 output to extract sequence-to-SSR mapping.
    
    Extracts sequence identifiers in format {base_id}_{ssr_nr} and associated
    primer sequences for matching with the primer TSV.
    
    Args:
        p3out_file: Path to Primer3 output file (*.p3out format)
        
    Returns:
        Dictionary mapping Primer3 sequence IDs to primer info:
        {sequence_id: {'base_id': str, 'ssr_nr': int, 'left_primers': [...], 'right_primers': [...]}}
    """
    primer_mapping = {}
    current_sequence = None
    base_id = None
    ssr_nr = None
    
    with open(p3out_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            if line.startswith('SEQUENCE:'):
                current_sequence = line.split(':', 1)[1].strip()
                match = re.match(r'(.+)_(\d+)$', current_sequence)
                if match:
                    base_id = match.group(1)
                    ssr_nr = int(match.group(2))
                    primer_mapping[current_sequence] = {
                        'base_id': base_id,
                        'ssr_nr': ssr_nr,
                        'left_primers': [],
                        'right_primers': []
                    }
                else:
                    base_id = None
                    ssr_nr = None
                    current_sequence = None
            
            elif current_sequence and base_id is not None and ssr_nr is not None:
                if line.startswith('PRIMER_LEFT_') and '_SEQUENCE=' in line:
                    seq = line.split('=', 1)[1]
                    primer_mapping[current_sequence]['left_primers'].append(seq)
                elif line.startswith('PRIMER_RIGHT_') and '_SEQUENCE=' in line:
                    seq = line.split('=', 1)[1]
                    primer_mapping[current_sequence]['right_primers'].append(seq)
    
    return primer_mapping


def find_ssr_for_primer(
    base_id: str,
    primer_seq: str,
    primer_type: str,
    ssr_data: Dict[str, List[Dict]],
    primer_mapping: Dict[str, Dict],
    source_sequence_id: Optional[str] = None
) -> Optional[Dict]:
    """
    Match primer to SSR using Primer3 output-based matching.
    
    Matches primer sequences to left/right primers in Primer3 output
    to extract the SSR number, then looks up SSR metadata in MISA data.
    
    Args:
        base_id: Sequence ID without F/R suffix
        primer_seq: Primer DNA sequence
        primer_type: 'F', 'R', or 'U' (forward, reverse, unknown)
        ssr_data: MISA data parsed by parse_misa_file()
        primer_mapping: Primer3 data parsed by parse_primer3_output()
        source_sequence_id: Exact Primer3 sequence identifier if available
        
    Returns:
        SSR metadata dict if match found, None otherwise
    """
    if not primer_mapping:
        return None

    if source_sequence_id and source_sequence_id in primer_mapping:
        primers = primer_mapping[source_sequence_id]
        primer_list = primers['left_primers'] if primer_type == 'F' else primers['right_primers']
        if primer_seq.upper() in [p.upper() for p in primer_list]:
            mapped_base_id = primers['base_id']
            mapped_ssr_nr = primers['ssr_nr']
            if mapped_base_id in ssr_data:
                for ssr in ssr_data[mapped_base_id]:
                    if ssr['ssr_nr'] == mapped_ssr_nr:
                        return ssr

    for primers in primer_mapping.values():
        if primers['base_id'] == base_id:
            primer_list = primers['left_primers'] if primer_type == 'F' else primers['right_primers']
            if primer_seq.upper() in [p.upper() for p in primer_list]:
                if base_id in ssr_data:
                    for ssr in ssr_data[base_id]:
                        if ssr['ssr_nr'] == primers['ssr_nr']:
                            return ssr
    
    return None


def enhance_primers_with_ssr(
    primer_tsv: str,
    misa_file: str,
    p3out_file: Optional[str] = None,
    output_file: str = None
):
    """
    Main function to enhance primer TSV with SSR metadata.
    
    Augments primer TSV output from primers_to_tsv.py with SSR metadata
    by cross-referencing MISA file and optionally Primer3 output.
    
    Args:
        primer_tsv: Path to input primer TSV file
        misa_file: Path to MISA file (*.misa)
        p3out_file: Optional path to Primer3 output file (*.p3out)
        output_file: Output path (default: input with .ssr_enhanced.tsv suffix)
    """
    print(f"Parsing MISA file: {misa_file}...")
    ssr_data = parse_misa_file(misa_file)
    print(f"  Found SSR data for {len(ssr_data)} sequences")
    
    primer_mapping = None
    if p3out_file:
        print(f"Parsing Primer3 output: {p3out_file}...")
        primer_mapping = parse_primer3_output(p3out_file)
        print(f"  Found primer mappings for {len(primer_mapping)} sequence-SSR pairs")
    
    print(f"Loading primer TSV: {primer_tsv}...")
    df = pd.read_csv(primer_tsv, sep='\t')
    print(f"  Loaded {len(df)} primer records")

    # Preserve the chromosome values exactly as provided by the input TSV.
    # The enhancer should only append SSR metadata, not reinterpret coordinates.
    preserved_chromosome = df['chromosome'].copy() if 'chromosome' in df.columns else None
    
    df['ssr_nr'] = ''
    df['ssr_type'] = ''
    df['ssr_motif'] = ''
    df['ssr_size'] = ''
    df['ssr_start'] = ''
    df['ssr_end'] = ''
    df['ssr_position_rel'] = ''
    df['amplicon_contains_ssr'] = ''
    
    matched_count = 0
    unmatched_count = 0
    
    for idx, row in df.iterrows():
        seq_id_full = row['sequence_id']
        match = re.match(r'(.+)_(F|R|U)$', seq_id_full)
        if match:
            base_id = match.group(1)
            primer_type = match.group(2)
        else:
            base_id = seq_id_full
            primer_type = 'U'
        
        primer_seq = row['primer_sequence']
        source_sequence_id = ''
        if 'source_sequence_id' in df.columns and pd.notna(row['source_sequence_id']):
            source_sequence_id = str(row['source_sequence_id']).strip()
            if source_sequence_id.lower() == 'nan':
                source_sequence_id = ''
        
        ssr = find_ssr_for_primer(
            base_id=base_id,
            primer_seq=primer_seq,
            primer_type=primer_type,
            ssr_data=ssr_data,
            primer_mapping=primer_mapping,
            source_sequence_id=source_sequence_id or None
        )
        
        if ssr:
            df.at[idx, 'ssr_nr'] = ssr['ssr_nr']
            df.at[idx, 'ssr_type'] = ssr['ssr_type']
            df.at[idx, 'ssr_motif'] = ssr['ssr_motif']
            df.at[idx, 'ssr_size'] = ssr['ssr_size']
            df.at[idx, 'ssr_start'] = ssr['ssr_start']
            df.at[idx, 'ssr_end'] = ssr['ssr_end']
            df.at[idx, 'ssr_position_rel'] = 'within_amplicon'
            df.at[idx, 'amplicon_contains_ssr'] = 'TRUE'
            
            matched_count += 1
        else:
            df.at[idx, 'ssr_nr'] = ''
            df.at[idx, 'ssr_type'] = 'NO_MATCH'
            df.at[idx, 'ssr_motif'] = ''
            df.at[idx, 'ssr_size'] = ''
            df.at[idx, 'ssr_start'] = ''
            df.at[idx, 'ssr_end'] = ''
            df.at[idx, 'ssr_position_rel'] = 'not_found'
            df.at[idx, 'amplicon_contains_ssr'] = 'FALSE'
            
            unmatched_count += 1
    
    if output_file is None:
        output_file = primer_tsv.replace('.tsv', '.ssr_enhanced.tsv')

    if preserved_chromosome is not None:
        df['chromosome'] = preserved_chromosome
    
    print(f"\nWriting enhanced output to: {output_file}...")
    df.to_csv(output_file, sep='\t', index=False)
    
    print(f"\nSummary:")
    print(f"  Total primers: {len(df)}")
    print(f"  Matched with SSR: {matched_count} ({matched_count/len(df)*100:.1f}%)")
    print(f"  Unmatched: {unmatched_count} ({unmatched_count/len(df)*100:.1f}%)")
    
    ssr_types = df[df['ssr_type'] != '']['ssr_type'].value_counts()
    if not ssr_types.empty:
        print(f"\nSSR Type Distribution:")
        for ssr_type, count in ssr_types.items():
            if ssr_type != 'NO_MATCH':
                print(f"  {ssr_type}: {count}")


def main():
    parser = argparse.ArgumentParser(
        description='Enhance primer TSV with SSR metadata from MISA files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage (position-based matching)
  python enhance_primers_with_ssr.py -i primers.tsv -m sequences.fa.misa -o primers.ssr.tsv
  
  # With Primer3 output for improved matching
  python enhance_primers_with_ssr.py -i primers.tsv -m sequences.fa.misa -p sequences.p3out -o primers.ssr.tsv
  
  # For defense candidates
  python enhance_primers_with_ssr.py \\
    -i defense_candidates.primers.tsv \\
    -m defense_candidates.fa.misa \\
    -p defense_candidates.p3out \\
    -o defense_candidates.primers.ssr_enhanced.tsv
        """
    )
    
    parser.add_argument('-i', '--input', required=True,
                       help='Input primer TSV file (from primers_to_tsv.py)')
    parser.add_argument('-m', '--misa', required=True,
                       help='MISA file (*.misa) with SSR metadata')
    parser.add_argument('-p', '--primer3-output', default=None,
                       help='Optional: Primer3 output file (*.p3out) for improved matching')
    parser.add_argument('-o', '--output', default=None,
                       help='Output TSV file (default: input with .ssr_enhanced.tsv suffix)')
    
    args = parser.parse_args()
    
    try:
        enhance_primers_with_ssr(
            primer_tsv=args.input,
            misa_file=args.misa,
            p3out_file=args.primer3_output,
            output_file=args.output
        )
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
