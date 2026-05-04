"""
Validate ecoli_sources/data/rnaseq_experimental/manifest.tsv and every TPM file it references.

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

from schemas import RnaseqSamplesManifestSchema, RnaseqTpmTableSchema  # noqa: E402

DATA_DIR = REPO_ROOT / "ecoli_sources" / "data" / "rnaseq_experimental"  # noqa: E402
MANIFEST = DATA_DIR / "manifest.tsv"


def _fail(label: str, err: Exception) -> None:
    print(f"FAIL {label}", file=sys.stderr)
    print(str(err)[:2000], file=sys.stderr)


def main() -> int:
    if not MANIFEST.exists():
        print(f"FAIL missing manifest: {MANIFEST}", file=sys.stderr)
        return 1

    manifest = pd.read_csv(MANIFEST, sep="\t", comment="#")
    try:
        RnaseqSamplesManifestSchema.validate(manifest, lazy=True)
    except Exception as e:
        _fail(f"RnaseqSamplesManifestSchema on {MANIFEST}", e)
        return 1
    print(f"OK manifest ({len(manifest)} rows)")

    failures: list[str] = []
    for _, row in manifest.iterrows():
        rel = row["file_path"]
        path = (DATA_DIR / rel).resolve()
        dataset_id = row["dataset_id"]
        if not path.exists():
            # Generated perturbation variants live under
            # ecoli_sources/data/rnaseq_experimental/perturbations/ and are
            # gitignored; skip with a note rather than fail CI.
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

    if failures:
        print(f"\n{len(failures)} dataset(s) failed validation:", file=sys.stderr)
        for d in failures:
            print(f"  - {d}", file=sys.stderr)
        return 1

    print(f"\nAll {len(manifest)} datasets validated.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
