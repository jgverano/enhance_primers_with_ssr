# Enhance Primers With SSR

## About

Enhance Primers With SSR is a small bioinformatics command-line utility for
adding Simple Sequence Repeat (SSR) metadata from MISA output to primer TSV
files. It is designed for primer design workflows where candidate primers need
to be connected back to their SSR targets for filtering, validation, or
downstream analysis.

This utility is intended for primer design workflows where primers were exported
to TSV and SSR calls are available from a `.misa` file. If Primer3 output is
also available, the script can use it to match primers back to the exact SSR
target more reliably.

## What It Does

`enhance_primers_with_ssr.py` reads:

- a primer TSV file, usually produced by `primers_to_tsv.py`
- a MISA SSR file, usually named like `sequences.fa.misa`
- optionally, a Primer3 output file

It writes a new TSV with appended SSR fields:

- `ssr_nr`
- `ssr_type`
- `ssr_motif`
- `ssr_size`
- `ssr_start`
- `ssr_end`
- `ssr_position_rel`
- `amplicon_contains_ssr`

## Installation

```bash
git clone <remote-url>
cd enhance_primers_with_ssr
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

For a minimal setup without editable installation:

```bash
python -m pip install -r requirements.txt
```

## Usage

After installation:

```bash
enhance-primers-with-ssr \
  -i primers.tsv \
  -m sequences.fa.misa \
  -p sequences.p3out \
  -o primers.ssr_enhanced.tsv
```

You can also run the script directly:

```bash
python enhance_primers_with_ssr/enhance_primers_with_ssr.py \
  -i primers.tsv \
  -m sequences.fa.misa \
  -p sequences.p3out \
  -o primers.ssr_enhanced.tsv
```

`-p/--primer3-output` is optional. Matching is strongest when Primer3 output is
provided because primer sequences can be mapped back to a specific SSR number.

## Input Requirements

The primer TSV must contain at least:

- `sequence_id`
- `primer_sequence`

The script also uses `source_sequence_id` when present. This is helpful when the
Primer3 sequence ID contains the SSR number, for example:

```text
IHE45_01G034600.1_1
```

The MISA file should contain the standard tab-separated MISA columns:

```text
ID    SSR nr.    SSR type    SSR    size    start    end
```

## Example

Run the included sample files:

```bash
enhance-primers-with-ssr \
  -i examples/primers.tsv \
  -m examples/sequences.fa.misa \
  -p examples/sequences.p3out \
  -o examples/primers.ssr_enhanced.tsv
```

Expected summary:

```text
Total primers: 2
Matched with SSR: 2
Unmatched: 0
```

## Development

```bash
python -m pip install -e .
python -m py_compile enhance_primers_with_ssr/enhance_primers_with_ssr.py
```

## License

This project is released under the MIT License. See `LICENSE`.
