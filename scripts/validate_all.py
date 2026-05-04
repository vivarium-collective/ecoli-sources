"""
Validate the reference bundle manifest and the rnaseq experimental manifest,
plus every file each one references.

Validation order:

1. ``data/reference_bundle.tsv`` — top-level bundle manifest. Validates
   against ``ReferenceBundleSchema``; every ``source_path`` must resolve
   to an existing file or directory under ``data/``.
2. ``data/rnaseq_experimental/manifest.tsv`` — rnaseq sub-registry of
   available TPM datasets. Validates against ``RnaseqSamplesManifestSchema``;
   every referenced TPM file is validated against ``RnaseqTpmTableSchema``.

Exits non-zero on the first failure and prints a summary. Intended for use
in CI; also runnable locally:

    uv run python scripts/validate_all.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import schemas as _schemas_module  # noqa: E402
from schemas import (  # noqa: E402
    ReferenceBundleSchema,
    RnaseqSamplesManifestSchema,
    RnaseqTpmTableSchema,
)

DATA_DIR = REPO_ROOT / "ecoli_sources" / "data"
BUNDLE = DATA_DIR / "reference_bundle.tsv"
RNASEQ_DIR = DATA_DIR / "rnaseq_experimental"
MANIFEST = RNASEQ_DIR / "manifest.tsv"


def _fail(label: str, err: Exception) -> None:
    print(f"FAIL {label}", file=sys.stderr)
    print(str(err)[:2000], file=sys.stderr)


def validate_bundle() -> list[str]:
    """Validate reference_bundle.tsv schema and every source_path. Returns list of failures."""
    failures: list[str] = []

    if not BUNDLE.exists():
        _fail(f"missing bundle: {BUNDLE}", FileNotFoundError(str(BUNDLE)))
        return ["<bundle missing>"]

    bundle = pd.read_csv(BUNDLE, sep="\t", comment="#")
    try:
        ReferenceBundleSchema.validate(bundle, lazy=True)
    except Exception as e:
        _fail(f"ReferenceBundleSchema on {BUNDLE}", e)
        return ["<bundle schema>"]
    print(f"OK reference_bundle ({len(bundle)} rows)")

    for _, row in bundle.iterrows():
        canonical_key = row["canonical_key"]
        rel = row["source_path"]
        path = (DATA_DIR / rel).resolve()
        if not path.exists():
            _fail(f"{canonical_key}: missing source_path {path}", FileNotFoundError(str(path)))
            failures.append(canonical_key)
            continue
        kind = "dir" if path.is_dir() else "file"

        # When schema_name is set on the bundle row, content-validate the
        # file against that schema. Empty/missing schema_name is fine —
        # schemas are accreted incrementally.
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
                continue
            print(f"OK {canonical_key} ({kind}: {rel}) [validated against {schema_name}]")
        else:
            print(f"OK {canonical_key} ({kind}: {rel})")

    return failures


def validate_rnaseq_manifest() -> list[str]:
    """Validate rnaseq manifest schema and every referenced TPM file. Returns list of failures."""
    failures: list[str] = []

    if not MANIFEST.exists():
        _fail(f"missing manifest: {MANIFEST}", FileNotFoundError(str(MANIFEST)))
        return ["<manifest missing>"]

    manifest = pd.read_csv(MANIFEST, sep="\t", comment="#")
    try:
        RnaseqSamplesManifestSchema.validate(manifest, lazy=True)
    except Exception as e:
        _fail(f"RnaseqSamplesManifestSchema on {MANIFEST}", e)
        return ["<manifest schema>"]
    print(f"OK rnaseq manifest ({len(manifest)} rows)")

    for _, row in manifest.iterrows():
        rel = row["file_path"]
        path = (RNASEQ_DIR / rel).resolve()
        dataset_id = row["dataset_id"]
        if not path.exists():
            # Generated perturbation variants live under
            # data/rnaseq_experimental/perturbations/ and are gitignored;
            # skip with a note rather than fail CI.
            if "perturbations/" in str(path) or path.parent.name == "perturbations":
                print(f"SKIP {dataset_id}: generated variant not present ({rel})")
                continue
            _fail(f"{dataset_id}: missing TPM file {path}", FileNotFoundError(str(path)))
            failures.append(dataset_id)
            continue

        df = pd.read_csv(path, sep="\t", comment="#")
        try:
            RnaseqTpmTableSchema.validate(df, lazy=True)
        except Exception as e:
            _fail(f"RnaseqTpmTableSchema on {rel} ({dataset_id})", e)
            failures.append(dataset_id)
            continue
        print(f"OK {dataset_id} ({len(df)} rows)")

    return failures


def main() -> int:
    bundle_failures = validate_bundle()
    rnaseq_failures = validate_rnaseq_manifest()

    failures = bundle_failures + rnaseq_failures
    if failures:
        print(f"\n{len(failures)} validation failure(s):", file=sys.stderr)
        for d in failures:
            print(f"  - {d}", file=sys.stderr)
        return 1

    print("\nAll validations passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
