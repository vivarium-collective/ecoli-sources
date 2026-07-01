"""
Unit test for the validation overlay-union loader (``ecoli_sources.validation``).

Exercises the loader against a synthetic primary+overlay fixture (built in a temp
dir) rather than the committed data, so it checks *loader logic* — the
``canonical_key`` union, both claim-table conventions, the uniqueness contract,
the cultivation views, and overlay-path parsing — independent of any real payload.

Follows the repo's script idiom (``scripts/validate_all.py``): plain asserts,
prints ``OK``/``FAIL`` per check, exits non-zero on the first failure. No pytest
dependency (functions are ``test_*`` so a future pytest run would also discover
them). Runnable locally:

    uv run python scripts/test_validation_loader.py
"""
from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from ecoli_sources import validation as V  # noqa: E402

APPROX = 1e-9


def _write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def build_fixture(root: Path) -> dict:
    """Create a primary (public, per-observable) + overlay (private, tidy) bundle
    fixture. Returns the paths the tests use."""
    primary = root / "primary"
    overlay = root / "overlay"

    # --- primary: per-observable ScalarClaimSchema (the public/literature shape)
    _write_tsv(pd.DataFrame([
        {"canonical_key": "basal__glucose_uptake",
         "source_path": "data/basal/glucose_uptake.tsv",
         "description": "literature glucose uptake",
         "schema_name": "ScalarClaimSchema"},
    ]), primary / "validation_bundle.tsv")
    _write_tsv(pd.DataFrame([
        {"source_id": "long_2017", "value": 8.46, "units": "mmol/gDW/h",
         "kind": "measured", "uncertainty": 0.42},
        {"source_id": "lacroix_2015", "value": 8.59, "units": "mmol/gDW/h",
         "kind": "measured", "uncertainty": 1.42},
        {"source_id": "fba_thermo", "value": 11.0, "units": "mmol/gDW/h",
         "kind": "theoretical_max", "uncertainty": None},
    ]), primary / "data/basal/glucose_uptake.tsv")
    _write_tsv(pd.DataFrame([
        {"cultivation_id": "lit_basal", "mode": "literature", "notes": "aggregate"},
    ]), primary / "cultivations.tsv")

    # --- overlay: tidy ScalarObservationSchema (the private/CD1 shape). One assay
    #     table holds two observables across two cultivations; bundle rows select
    #     slots by (cultivation_id, observable).
    _write_tsv(pd.DataFrame([
        {"canonical_key": "cd1_ginkgo_wt_m9__glucose_uptake",
         "source_path": "fermentation.tsv", "description": "cd1 glc uptake",
         "schema_name": "ScalarObservationSchema",
         "cultivation_id": "cd1_ginkgo_wt_m9", "observable": "glucose_uptake"},
        {"canonical_key": "cd1_ginkgo_wt_m9__acetate_secretion",
         "source_path": "fermentation.tsv", "description": "cd1 acetate",
         "schema_name": "ScalarObservationSchema",
         "cultivation_id": "cd1_ginkgo_wt_m9", "observable": "acetate_secretion"},
    ]), overlay / "validation_bundle.tsv")
    _write_tsv(pd.DataFrame([
        # cd1_ginkgo_wt_m9 × glucose_uptake × 3 reactors
        {"cultivation_id": "cd1_ginkgo_wt_m9", "observable": "glucose_uptake",
         "value": 9.1, "units": "mmol/gDW/h", "kind": "measured", "replicate": "r1",
         "phase": "exponential_batch", "window": "3-4h"},
        {"cultivation_id": "cd1_ginkgo_wt_m9", "observable": "glucose_uptake",
         "value": 9.3, "units": "mmol/gDW/h", "kind": "measured", "replicate": "r2",
         "phase": "exponential_batch", "window": "3-4h"},
        {"cultivation_id": "cd1_ginkgo_wt_m9", "observable": "glucose_uptake",
         "value": 8.9, "units": "mmol/gDW/h", "kind": "measured", "replicate": "r3",
         "phase": "exponential_batch", "window": "3-4h"},
        # cd1_ginkgo_wt_m9 × acetate_secretion × 3 reactors
        {"cultivation_id": "cd1_ginkgo_wt_m9", "observable": "acetate_secretion",
         "value": 1.2, "units": "mmol/gDW/h", "kind": "measured", "replicate": "r1",
         "phase": "exponential_batch", "window": "3-4h"},
        {"cultivation_id": "cd1_ginkgo_wt_m9", "observable": "acetate_secretion",
         "value": 1.4, "units": "mmol/gDW/h", "kind": "measured", "replicate": "r2",
         "phase": "exponential_batch", "window": "3-4h"},
        {"cultivation_id": "cd1_ginkgo_wt_m9", "observable": "acetate_secretion",
         "value": 1.0, "units": "mmol/gDW/h", "kind": "measured", "replicate": "r3",
         "phase": "exponential_batch", "window": "3-4h"},
        # a DIFFERENT cultivation in the same table — must not leak into the above
        {"cultivation_id": "cd1_ginkgo_vio_m9", "observable": "glucose_uptake",
         "value": 7.0, "units": "mmol/gDW/h", "kind": "measured", "replicate": "r1",
         "phase": "exponential_batch", "window": "3-4h"},
    ]), overlay / "fermentation.tsv")
    _write_tsv(pd.DataFrame([
        {"cultivation_id": "cd1_ginkgo_wt_m9", "mode": "EXP", "performer": "Ginkgo"},
        {"cultivation_id": "v2ecoli_baseline", "mode": "SIM", "performer": "Stanford"},
    ]), overlay / "cultivations.tsv")

    # --- a second overlay that duplicates a primary key (collision fixture)
    _write_tsv(pd.DataFrame([
        {"canonical_key": "basal__glucose_uptake", "source_path": "unused.tsv",
         "description": "dup", "schema_name": "ScalarClaimSchema"},
    ]), root / "dup" / "validation_bundle.tsv")

    return {
        "primary_bundle": primary / "validation_bundle.tsv",
        "overlay_bundle": overlay / "validation_bundle.tsv",
        "dup_bundle": root / "dup" / "validation_bundle.tsv",
        "primary_reg": primary / "cultivations.tsv",
        "overlay_reg": overlay / "cultivations.tsv",
    }


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

def test_union_both_conventions(fx):
    m = V.load_scalar_observations(fx["primary_bundle"], [fx["overlay_bundle"]])
    assert set(m) == {
        "basal__glucose_uptake",
        "cd1_ginkgo_wt_m9__glucose_uptake",
        "cd1_ginkgo_wt_m9__acetate_secretion",
    }, f"unexpected union keys: {sorted(m)}"


def test_per_observable_band(fx):
    b = V.load_scalar_observations(fx["primary_bundle"], [])["basal__glucose_uptake"]
    assert sorted(b["measured"]) == [8.46, 8.59], b["measured"]
    assert abs(b["theoretical_max"] - 11.0) < APPROX
    assert b["theoretical_source"] == "fba_thermo"
    assert b["observable"] == "glucose_uptake", b["observable"]  # inferred from key tail
    assert b["cultivation_id"] is None
    assert b["phase"] is None and b["window"] is None  # per-observable carries no temporal tier
    assert b["units"] == "mmol/gDW/h"
    assert set(b["sources"]) == {"long_2017", "lacroix_2015"}


def test_tidy_band(fx):
    b = V.load_scalar_observations(fx["primary_bundle"], [fx["overlay_bundle"]])[
        "cd1_ginkgo_wt_m9__glucose_uptake"]
    assert sorted(b["measured"]) == [8.9, 9.1, 9.3], b["measured"]  # 3 reactors, not the vio row
    assert b["cultivation_id"] == "cd1_ginkgo_wt_m9"
    assert b["observable"] == "glucose_uptake"
    assert b["phase"] == "exponential_batch"
    assert b["window"] == "3-4h"
    assert set(b["sources"]) == {"r1", "r2", "r3"}  # source_col = replicate


def test_observations_for_cultivation(fx):
    obs = V.observations_for_cultivation(
        "cd1_ginkgo_wt_m9", fx["primary_bundle"], [fx["overlay_bundle"]])
    assert set(obs) == {"glucose_uptake", "acetate_secretion"}, sorted(obs)
    assert obs["acetate_secretion"]["cultivation_id"] == "cd1_ginkgo_wt_m9"
    assert sorted(obs["acetate_secretion"]["measured"]) == [1.0, 1.2, 1.4]


def test_prefix_filter(fx):
    m = V.load_scalar_observations(
        fx["primary_bundle"], [fx["overlay_bundle"]], prefix="cd1_ginkgo_wt_m9__")
    assert set(m) == {
        "cd1_ginkgo_wt_m9__glucose_uptake",
        "cd1_ginkgo_wt_m9__acetate_secretion",
    }, f"prefix filter leaked: {sorted(m)}"


def test_canonical_key_collision(fx):
    try:
        V.load_scalar_observations(fx["primary_bundle"], [fx["dup_bundle"]])
    except ValueError as e:
        assert "duplicate validation canonical_key" in str(e), str(e)
        return
    raise AssertionError("expected ValueError on duplicate canonical_key, none raised")


def test_cultivations_union(fx):
    c = V.load_cultivations(fx["primary_reg"], [fx["overlay_reg"]])
    assert set(c) >= {"lit_basal", "cd1_ginkgo_wt_m9", "v2ecoli_baseline"}, sorted(c)
    assert c["cd1_ginkgo_wt_m9"]["mode"] == "EXP"
    assert c["v2ecoli_baseline"]["mode"] == "SIM"  # a SIM is a cultivation too


def test_cultivation_id_collision(fx):
    try:
        V.load_cultivations(fx["primary_reg"], [fx["primary_reg"]])  # same reg twice
    except ValueError as e:
        assert "duplicate cultivation_id" in str(e), str(e)
        return
    raise AssertionError("expected ValueError on duplicate cultivation_id, none raised")


def test_overlay_separator_parsing(fx):
    assert V._split_overlay_separator("a;b;c") == ["a", "b", "c"]
    assert V._split_overlay_separator("a:b") == ["a", "b"]  # legacy separator
    # a cloud URI must not be torn apart by the split
    assert V._split_overlay_separator("s3://bucket/x;./local") == ["s3://bucket/x", "./local"]
    # env-var driven overlay paths
    os.environ[V.VALIDATION_OVERLAYS_ENV] = str(fx["overlay_bundle"])
    try:
        assert V.validation_overlay_paths() == [Path(str(fx["overlay_bundle"]))]
    finally:
        del os.environ[V.VALIDATION_OVERLAYS_ENV]


CHECKS = [
    test_union_both_conventions,
    test_per_observable_band,
    test_tidy_band,
    test_observations_for_cultivation,
    test_prefix_filter,
    test_canonical_key_collision,
    test_cultivations_union,
    test_cultivation_id_collision,
    test_overlay_separator_parsing,
]


def main() -> int:
    with tempfile.TemporaryDirectory() as td:
        fx = build_fixture(Path(td))
        failures = []
        for check in CHECKS:
            try:
                check(fx)
                print(f"OK {check.__name__}")
            except Exception as e:  # noqa: BLE001 — report and continue
                print(f"FAIL {check.__name__}", file=sys.stderr)
                print(f"  {type(e).__name__}: {e}", file=sys.stderr)
                failures.append(check.__name__)
    if failures:
        print(f"\n{len(failures)} loader check(s) failed:", file=sys.stderr)
        for f in failures:
            print(f"  - {f}", file=sys.stderr)
        return 1
    print(f"\nAll {len(CHECKS)} loader checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
