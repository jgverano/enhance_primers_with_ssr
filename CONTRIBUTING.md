# Contributing

Contributions are welcome.

## Local Setup

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

## Checks

Before opening a pull request, run:

```bash
python -m py_compile enhance_primers_with_ssr/enhance_primers_with_ssr.py
```

If you change input or output behavior, add a small example file or update the
README usage notes.
