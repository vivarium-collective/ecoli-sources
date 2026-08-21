"""Unit test for the per-replicate vector contract, and for the coupling that
registering a schema does NOT by itself validate it.

Exercises ``schemas.VectorReplicateSchema`` and
``ecoli_sources.validation.read_vector_observations`` against synthetic fixtures
built in a temp dir. The payloads that motivate this tier live in private
overlays, so the mechanism has to be testable without any of them.

⚠ The first version of this schema was landed WITHOUT this file, and a cold
review found the resulting defect by execution: the schema was added to
``_VECTOR_SCHEMAS`` (making a slot resolvable) while the loader's validate
branch still tested one literal schema name, so a caller passing
``validate=True`` had a corrupt table returned to it with no error.
:func:`test_registering_a_schema_also_validates_it` is that regression.

Run with::

    uv run pytest scripts/test_vector_replicates.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ecoli_sources.validation import _VECTOR_SCHEMAS, _schema_for  # noqa: E402
from schemas.vector_replicate import (  # noqa: E402
    DETECTION_STATES,
    TIME_BASES,
    VectorReplicateSchema,
)


def _rows() -> pd.DataFrame:
    """Two coherent per-sample rows: one detected, one never looked for."""
    return pd.DataFrame([
        dict(cultivation_group_id="grp_a", observable="proteome", entity_id="E1",
             units="normalized_abundance", kind="measured", replicate_id="r1_t3",
             replicate="1", replicate_basis="reactor", time_basis="instant",
             timepoint_h=3.0, sample_id="S1", detection="detected", value=0.5,
             phase="exponential_batch", window="3-4h"),
        dict(cultivation_group_id="grp_a", observable="proteome", entity_id="E1",
             units="normalized_abundance", kind="measured", replicate_id="r2_t3",
             replicate="2", replicate_basis="reactor", time_basis="instant",
             timepoint_h=3.0, sample_id="S2", detection="not_detected", value=None,
             phase="exponential_batch", window="3-4h"),
    ])


def test_a_coherent_table_validates():
    """Guards the negatives below: if this fails they prove nothing."""
    VectorReplicateSchema.validate(_rows(), lazy=True)


@pytest.mark.parametrize("mutate,why", [
    (lambda d: d.assign(value=[0.5, 0.7]),
     "a value on a not_detected row invents a measurement"),
    (lambda d: d.assign(value=[None, None]),
     "a detected row with no value drops silently from any re-aggregation"),
    (lambda d: pd.concat([d, d.head(1)]),
     "a duplicate key double-weights one sample: invisible in the mean, wrong in n"),
    (lambda d: d.assign(detection=["detcted", "not_detected"]),
     "a detection token outside the vocabulary"),
    (lambda d: d.assign(replicate_basis=["unicorn", "reactor"]),
     "a basis outside the cultivation registry's vocabulary"),
    (lambda d: d.assign(time_basis=["whenever", "instant"]),
     "a time basis outside the vocabulary"),
    (lambda d: d.assign(detection=["below_limit", "not_detected"]),
     "a below_limit row carrying a number: the aggregate tier refuses one, and "
     "a re-aggregating consumer would average it as if it were a measurement"),
    (lambda d: d.assign(timepoint_h=[None, None]),
     "an 'instant' row with no timepoint — the exact ambiguity time_basis exists "
     "to remove"),
    (lambda d: d.assign(kind=["simulated", "measured"]),
     "a plausible-looking kind outside the vocabulary would validate and then "
     "fall out of every kind == 'measured' filter downstream"),
])
def test_incoherent_rows_are_rejected(mutate, why):
    with pytest.raises(Exception):
        VectorReplicateSchema.validate(mutate(_rows()), lazy=True)


def test_below_limit_carries_no_number_like_the_aggregate_tier():
    """The two tiers must agree on what a token MEANS, not merely on its spelling.

    An earlier revision required a number here while ``VectorObservationSchema``
    refused one — opposite contracts for ``below_limit``, with the producer
    emitting ``0.0``. Pin the agreement in both directions.
    """
    ok = _rows().assign(detection=["below_limit", "not_detected"], value=[None, None])
    VectorReplicateSchema.validate(ok, lazy=True)

    from schemas.vector_observation import VectorObservationSchema
    agg = pd.DataFrame([dict(
        cultivation_group_id="grp_a", observable="proteome", entity_id="E1",
        symbol="e1", units="normalized_abundance", kind="measured",
        detection="below_limit", mean_arithmetic=None, mean_geometric=None,
        sd_log10=None, n=0, n_pos=0, n_total=6)])
    VectorObservationSchema.validate(agg, lazy=True)  # same claim, same shape


def test_interval_mean_rows_need_no_timepoint():
    """The other half of F4: an interval may be untimed, and that is legitimate."""
    ok = _rows().assign(time_basis="interval_mean", timepoint_h=[None, None])
    VectorReplicateSchema.validate(ok, lazy=True)


def test_vocabularies_are_shared_not_restated():
    """Two independent lists are free to drift into disagreeing while both validate."""
    from schemas.vector_observation import _DETECTION
    assert DETECTION_STATES == tuple(_DETECTION)
    from schemas.cultivation import _REPLICATE_BASES
    # Assert the BEHAVIOUR rather than the check's internals: every registry
    # basis must validate here, and a non-basis must not. Reaching into the
    # Check object couples this test to pandera's private shape.
    for basis in _REPLICATE_BASES:
        VectorReplicateSchema.validate(_rows().assign(replicate_basis=basis), lazy=True)
    with pytest.raises(Exception):
        VectorReplicateSchema.validate(_rows().assign(replicate_basis="unicorn"),
                                       lazy=True)
    assert TIME_BASES == ("instant", "interval_mean")


def _slot(tmp_path: Path, frame: pd.DataFrame, schema_name: str) -> dict:
    p = tmp_path / "t.tsv"
    frame.to_csv(p, sep="\t", index=False)
    return {"path": p, "schema_name": schema_name, "canonical_key": "k"}


def test_registering_a_schema_also_validates_it(tmp_path):
    """★ THE REGRESSION. Fails against a loader that tests one literal name.

    A slot whose schema is registered in ``_VECTOR_SCHEMAS`` must be VALIDATED
    when the caller asks for it, not merely located. Before this test the
    replicate tier was resolvable and ungraded, and a corrupt table came back
    clean.
    """
    from ecoli_sources.validation import read_vector_observations

    corrupt = _rows().assign(detection=["detcted", "not_detected"])
    with pytest.raises(Exception):
        read_vector_observations(_slot(tmp_path, corrupt, "VectorReplicateSchema"),
                                 validate=True)


def test_validate_false_still_skips(tmp_path):
    """The escape hatch has to keep working, or renders break on payload defects."""
    from ecoli_sources.validation import read_vector_observations

    corrupt = _rows().assign(detection=["detcted", "not_detected"])
    out = read_vector_observations(_slot(tmp_path, corrupt, "VectorReplicateSchema"),
                                   validate=False)
    assert len(out) == 2


def test_every_registered_schema_this_repo_defines_is_dispatched():
    """Catches the NEXT instance of F1 without anyone knowing where it is.

    ``_VECTOR_SCHEMAS`` legitimately contains names validated elsewhere, so a
    ``None`` from ``_schema_for`` is not automatically wrong. What is wrong is a
    schema module that EXISTS in this repo and is registered, yet is never
    dispatched — registration without validation.
    """
    import importlib

    for name in _VECTOR_SCHEMAS:
        module = "schemas." + "".join(
            "_" + c.lower() if c.isupper() else c
            for c in name.removesuffix("Schema")).lstrip("_")
        try:
            importlib.import_module(module)
        except ModuleNotFoundError:
            continue  # validated elsewhere; not this repo's contract
        assert _schema_for(name) is not None, (
            f"{name} is registered in _VECTOR_SCHEMAS and its schema module "
            f"exists here, but _schema_for does not dispatch it — slots resolve "
            f"and validate=True checks nothing")
