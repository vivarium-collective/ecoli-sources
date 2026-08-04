"""
Unit test for the vector-observation contract and the loader's vector slots.

Exercises ``schemas.VectorObservationSchema`` and
``ecoli_sources.validation.{load_vector_observations,read_vector_observations}``
against synthetic fixtures built in a temp dir. The payloads that motivate this
tier are experimental data living in private overlays, so the mechanism has to
be testable without any of them.

Run with::

    uv run pytest scripts/test_vector_observations.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from ecoli_sources import validation as V  # noqa: E402
from schemas import VectorObservationSchema  # noqa: E402


@pytest.fixture
def root(tmp_path) -> Path:
    """Fresh temp root; the collision check builds its own payload under it."""
    return tmp_path


@pytest.fixture
def bundle(root: Path) -> Path:
    """Path to a synthetic validation bundle manifest."""
    return build_fixture(root)


def _row(**over) -> dict:
    """One valid detected observation; keyword args override any field."""
    row = {
        "cultivation_group_id": "campaign_perf_wt_m9",
        "observable": "proteome",
        "entity_id": "MONOMER-A",
        "symbol": "aaa",
        "units": "normalized_abundance",
        "kind": "measured",
        "detection": "detected",
        "mean_arithmetic": 10.0,
        "mean_geometric": 9.5,
        "sd_log10": 0.2,
        "n": 6, "n_pos": 6, "n_total": 6,
        "phase": "exponential_batch", "window": "3-4h",
        "method": "DIA proteomics", "notes": "",
    }
    row.update(over)
    return row


def _frame(*rows: dict) -> pd.DataFrame:
    return pd.DataFrame(list(rows))


def _rejects(df: pd.DataFrame, fragment: str) -> None:
    try:
        VectorObservationSchema.validate(df, lazy=True)
    except Exception as exc:
        assert fragment in str(exc), f"expected {fragment!r} in error, got: {exc}"
        return
    raise AssertionError(f"expected rejection mentioning {fragment!r}, but it passed")


def _write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def build_fixture(root: Path) -> Path:
    """A bundle with two vector slots (two cultivation groups) + one scalar
    slot, so the loader's schema split is actually exercised."""
    bundle = root / "validation_data"

    _write_tsv(_frame(
        _row(),
        _row(entity_id="MONOMER-B", symbol="bbb", units="mass_fraction",
             mean_arithmetic=0.02, mean_geometric=0.019, sd_log10=0.1),
        # a genuine measured zero — detected, value 0.0, no positives
        _row(entity_id="MONOMER-C", symbol="ccc", mean_arithmetic=0.0,
             mean_geometric=None, sd_log10=None, n=6, n_pos=0),
        # looked for, absent: no value at all
        _row(entity_id="MONOMER-D", symbol="ddd", detection="not_detected",
             mean_arithmetic=None, mean_geometric=None, sd_log10=None,
             n=0, n_pos=0),
    ), bundle / "vectors" / "proteome__wt__exponential_batch.tsv")

    _write_tsv(_frame(
        _row(cultivation_group_id="campaign_perf_mut_m9", n_total=4, n=4, n_pos=4),
    ), bundle / "vectors" / "proteome__mut__exponential_batch.tsv")

    _write_tsv(_frame({
        "cultivation_group_id": ["campaign_perf_wt_m9"], "observable": ["growth_rate"],
        "value": [0.8], "units": ["1/h"], "kind": ["measured"], "replicate": ["r1"],
    }), bundle / "fermentation.tsv")

    _write_tsv(pd.DataFrame([
        {"canonical_key": "wt__proteome", "source_path": "vectors/proteome__wt__exponential_batch.tsv",
         "description": "WT proteome, exponential batch", "schema_name": "VectorObservationSchema",
         "cultivation_group_id": "campaign_perf_wt_m9", "observable": "proteome",
         "units": "normalized_abundance", "phase": "exponential_batch", "window": "3-4h"},
        {"canonical_key": "wt__proteome_mass_fraction",
         "source_path": "vectors/proteome__wt__exponential_batch.tsv",
         "description": "WT proteome as mass fraction", "schema_name": "VectorObservationSchema",
         "cultivation_group_id": "campaign_perf_wt_m9", "observable": "proteome",
         "units": "mass_fraction", "phase": "exponential_batch", "window": "3-4h"},
        {"canonical_key": "mut__proteome", "source_path": "vectors/proteome__mut__exponential_batch.tsv",
         "description": "mutant proteome", "schema_name": "VectorObservationSchema",
         "cultivation_group_id": "campaign_perf_mut_m9", "observable": "proteome",
         "units": "normalized_abundance", "phase": "exponential_batch", "window": "3-4h"},
        {"canonical_key": "wt__growth_rate", "source_path": "fermentation.tsv",
         "description": "WT growth rate", "schema_name": "ScalarObservationSchema",
         "cultivation_group_id": "campaign_perf_wt_m9", "observable": "growth_rate"},
    ]), bundle / "validation_bundle.tsv")

    return bundle / "validation_bundle.tsv"


def test_schema_accepts_the_valid_shapes() -> None:
    df = _frame(
        _row(),
        _row(entity_id="MONOMER-C", mean_arithmetic=0.0, mean_geometric=None,
             sd_log10=None, n=6, n_pos=0),
        _row(entity_id="MONOMER-D", detection="not_detected", mean_arithmetic=None,
             mean_geometric=None, sd_log10=None, n=0, n_pos=0),
        _row(entity_id="MONOMER-E", mean_geometric=9.0, sd_log10=None, n=6, n_pos=1),
    )
    VectorObservationSchema.validate(df, lazy=True)
    print("OK   accepts a detected value, a genuine zero, a non-detection, and n_pos=1")


def test_schema_rejects_incoherent_rows() -> None:
    # a non-detection carrying a number
    _rejects(_frame(_row(detection="not_detected")), "value_matches_detection")
    # a detection missing its number
    _rejects(_frame(_row(mean_arithmetic=None)), "value_matches_detection")
    # counts not nested
    _rejects(_frame(_row(n=6, n_pos=7)), "counts_are_nested")
    _rejects(_frame(_row(n=6, n_pos=6, n_total=4)), "counts_are_nested")
    # log statistics without the positives to support them
    _rejects(_frame(_row(sd_log10=0.3, n_pos=1)), "log_statistics_have_support")
    _rejects(_frame(_row(mean_geometric=1.0, sd_log10=None, n_pos=0)),
             "log_statistics_have_support")
    # closed vocabularies
    _rejects(_frame(_row(detection="missing")), "detection")
    _rejects(_frame(_row(kind="guessed")), "kind")
    # negative abundance
    _rejects(_frame(_row(mean_arithmetic=-1.0)), "mean_arithmetic")
    # duplicate row key
    _rejects(_frame(_row(), _row()), "unique")
    print("OK   rejects incoherent detection/value, counts, log support, vocab, duplicates")


def test_loader_locates_vector_slots(bundle: Path) -> None:
    slots = V.load_vector_observations(bundle, overlay_paths=[])
    assert set(slots) == {"wt__proteome", "wt__proteome_mass_fraction", "mut__proteome"}, slots
    wt = slots["wt__proteome"]
    assert wt["cultivation_group_id"] == "campaign_perf_wt_m9"
    assert wt["units"] == "normalized_abundance"
    assert wt["window"] == "3-4h"
    assert Path(wt["path"]).is_file()

    # the scalar slot in the same bundle must NOT appear here, and vice versa
    scalars = V.load_scalar_observations(bundle, overlay_paths=[])
    assert set(scalars) == {"wt__growth_rate"}, scalars
    print("OK   locates vector slots only; scalar slots stay with the scalar loader")


def test_reading_a_slot_applies_its_filter(bundle: Path) -> None:
    slots = V.load_vector_observations(bundle, overlay_paths=[])

    # Two slots share ONE file, differing only by units — the filter is what
    # separates them, so a units-blind read would silently merge two vectors.
    norm = V.read_vector_observations(slots["wt__proteome"])
    mass = V.read_vector_observations(slots["wt__proteome_mass_fraction"])
    assert len(norm) == 3, len(norm)
    assert len(mass) == 1, len(mass)
    assert set(norm["units"]) == {"normalized_abundance"}
    assert set(mass["units"]) == {"mass_fraction"}

    assert set(norm["detection"]) == {"detected", "not_detected"}
    zero = norm[norm.entity_id == "MONOMER-C"].iloc[0]
    assert zero["mean_arithmetic"] == 0.0 and pd.isna(zero["mean_geometric"])
    print("OK   reading a slot applies its (group, observable, units) filter")


def test_key_collisions_are_caught_across_the_schema_split(root: Path) -> None:
    """A vector row must not shadow a scalar key: they address the same thing."""
    bundle = root / "collide"
    _write_tsv(_frame(_row()), bundle / "v.tsv")
    _write_tsv(pd.DataFrame([
        {"canonical_key": "dup__slot", "source_path": "v.tsv", "description": "vector",
         "schema_name": "VectorObservationSchema"},
        {"canonical_key": "dup__slot", "source_path": "v.tsv", "description": "scalar",
         "schema_name": "ScalarObservationSchema"},
    ]), bundle / "validation_bundle.tsv")
    for fn in (V.load_vector_observations, V.load_scalar_observations):
        try:
            fn(bundle / "validation_bundle.tsv", overlay_paths=[])
        except ValueError as exc:
            assert "duplicate validation canonical_key" in str(exc), exc
        else:
            raise AssertionError(f"{fn.__name__} did not raise on a cross-schema collision")
    print("OK   canonical_key collisions are caught across the scalar/vector split")


