"""
ecoli-sources: curated E. coli omics datasets for vEcoli.

Locator package. Provides ``DATA_DIR`` and ``BUNDLE_PATH`` for resolving
the package's data files at runtime — both in development (when
``ecoli-sources`` is installed editable) and when installed via a
git-pinned dependency.

Typical use from vEcoli or other consumers::

    from ecoli_sources import BUNDLE_PATH
    # BUNDLE_PATH points at the default reference_bundle.tsv

The package's other top-level directories (``schemas/``, ``processing/``,
``analysis/``) remain importable as their own top-level packages — this
module exists only to expose the data root.
"""

from pathlib import Path

DATA_DIR: Path = Path(__file__).resolve().parent / "data"
BUNDLE_PATH: Path = DATA_DIR / "reference_bundle.tsv"
