"""
Validate a bundle manifest TSV.

Three-stage check:

1. Schema — bundle file matches ``ReferenceBundleSchema`` (column shape +
   primary-key uniqueness + the canonical-key contract from
   ``REQUIRED_CANONICAL_KEYS``).
2. Path resolution — every ``source_path`` resolves to an existing file
   or directory. Relative paths are resolved against the bundle file's
   parent directory.
3. Content — for every row with ``schema_name`` set, the referenced
   file is validated against ``getattr(schemas, schema_name)``.

Usage:

    uv run python scripts/validate_bundle.py path/to/bundle.tsv

Exits non-zero on any failure. Intended for CI gating of campaign-generated
variant bundles before they're handed off to vEcoli's ParCa.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import schemas as _schemas_module  # noqa: E402
from schemas import ReferenceBundleSchema  # noqa: E402
from schemas.reference_bundle import (  # noqa: E402
    REQUIRED_CANONICAL_KEYS,
    _check_required_canonical_keys,
)


def _fail(label: str, err: Exception) -> None:
    print(f"FAIL {label}", file=sys.stderr)
    print(str(err)[:2000], file=sys.stderr)


def validate_bundle(bundle_path: Path) -> list[str]:
    """Run the three-stage validation. Returns a list of failure labels."""
    failures: list[str] = []
    bundle_path = bundle_path.resolve()
    bundle_root = bundle_path.parent

    if not bundle_path.is_file():
        _fail(f"bundle file not found: {bundle_path}", FileNotFoundError(str(bundle_path)))
        return ["<bundle missing>"]

    bundle = pd.read_csv(bundle_path, sep="\t", comment="#")

    # Stage 1: schema (shape + canonical-key contract)
    try:
        ReferenceBundleSchema.validate(bundle, lazy=True)
    except Exception as e:
        _fail(f"ReferenceBundleSchema on {bundle_path}", e)
        # Surface missing canonical keys explicitly when applicable
        missing = getattr(_check_required_canonical_keys, "last_missing", None)
        if missing:
            print(
                f"  Missing required canonical_key{'s' if len(missing) != 1 else ''} "
                f"({len(missing)}):",
                file=sys.stderr,
            )
            for k in missing:
                print(f"    - {k}", file=sys.stderr)
        return ["<bundle schema>"]
    print(f"OK schema ({len(bundle)} rows; {len(REQUIRED_CANONICAL_KEYS)} required keys present)")

    # Stage 2: path resolution
    # Stage 3: content validation (when schema_name set)
    for _, row in bundle.iterrows():
        canonical_key = row["canonical_key"]
        rel = row["source_path"]
        path = (bundle_root / rel).resolve()
        if not path.exists():
            _fail(f"{canonical_key}: missing source_path {path}", FileNotFoundError(str(path)))
            failures.append(canonical_key)
            continue

        schema_name = row.get("schema_name")
        if isinstance(schema_name, str) and schema_name and path.is_file():
            schema = getattr(_schemas_module, schema_name, None)
            if schema is None:
                _fail(
                    f"{canonical_key}: schema_name {schema_name!r} not exported from schemas package",
                    KeyError(schema_name),
                )
                failures.append(canonical_key)
                continue
            try:
                df = pd.read_csv(path, sep="\t", comment="#")
                schema.validate(df, lazy=True)
            except Exception as e:
                _fail(f"{canonical_key}: {schema_name} on {rel}", e)
                failures.append(canonical_key)

    if failures:
        print(
            f"\n{len(failures)} per-row failure(s); first few: {failures[:5]}",
            file=sys.stderr,
        )

    return failures


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate a bundle manifest TSV.")
    parser.add_argument(
        "bundle_path",
        type=Path,
        help="Path to a bundle manifest TSV (e.g., a campaign-generated variant).",
    )
    args = parser.parse_args()

    failures = validate_bundle(args.bundle_path)
    if failures:
        return 1
    print("\nBundle valid.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
