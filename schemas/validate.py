"""
CLI: validate a TSV against a named schema.

Usage
-----
    uv run python -m schemas.validate <schema_name> <path/to/file.tsv>
    uv run python -m schemas.validate --list

Examples
--------
    python -m schemas.validate AdjustmentValueSchema \\
        ../omics-vEcoli/reconstruction/ecoli/flat/adjustments/rna_expression_adjustments.tsv

Exits non-zero on failure; prints the first few problem rows for quick
diagnosis.
"""

from __future__ import annotations

import argparse
import sys

import pandas as pd

from . import __all__ as _schema_names
from . import __dict__ as _ns


def _schemas():
    return {name: _ns[name] for name in _schema_names}


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("schema", nargs="?", help="Schema name, e.g. 'AdjustmentValueSchema'")
    p.add_argument("path", nargs="?", help="Path to the TSV file to validate")
    p.add_argument("--list", action="store_true", help="List available schemas and exit")
    args = p.parse_args(argv)

    schemas = _schemas()

    if args.list or args.schema is None or args.path is None:
        print("Available schemas:")
        for name in sorted(schemas):
            print(f"  {name}")
        return 0 if args.list else 1

    if args.schema not in schemas:
        print(f"Unknown schema: {args.schema!r}", file=sys.stderr)
        return 2

    schema = schemas[args.schema]
    df = pd.read_csv(args.path, sep="\t", comment="#")
    try:
        schema.validate(df, lazy=True)
    except Exception as e:
        print(f"FAIL {args.schema} on {args.path}", file=sys.stderr)
        print(str(e)[:2000], file=sys.stderr)
        return 1

    print(f"OK {args.schema} on {args.path} ({len(df)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
