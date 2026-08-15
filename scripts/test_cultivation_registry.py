"""
Unit test for the cultivation-registry contract.

Exercises ``schemas.CultivationRegistrySchema``, with emphasis on the four
columns that carry provenance a consumer cannot reconstruct: ``n_reps`` /
``replicate_basis`` (how many, of what), ``variant`` (which arm), and
``run_commit`` (which code). Synthetic frames only -- the payloads that
motivate this tier live in private overlays, so the contract has to be
testable without them.

Run with::

    uv run pytest scripts/test_cultivation_registry.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from schemas import CultivationRegistrySchema  # noqa: E402


def _row(**over) -> dict:
    """A minimal valid registry row; override any field."""
    base = {
        "cultivation_group_id": "cd1_perf_wt_m9",
        "mode": "EXP",
        "performer": "Some Lab",
        "n_reps": 3,
        "replicate_basis": "reactor",
        "notes": "",
    }
    base.update(over)
    return base


def _frame(*rows) -> pd.DataFrame:
    return pd.DataFrame(list(rows) or [_row()])


# --- the required minimum ---------------------------------------------------

def test_minimal_row_validates():
    out = CultivationRegistrySchema.validate(_frame())
    assert len(out) == 1


def test_mode_is_constrained():
    with pytest.raises(Exception):
        CultivationRegistrySchema.validate(_frame(_row(mode="simulation")))


def test_group_id_is_unique():
    with pytest.raises(Exception):
        CultivationRegistrySchema.validate(
            _frame(_row(), _row(notes="a duplicate key")))


# --- n_reps / replicate_basis ----------------------------------------------

def test_n_reps_is_an_integer_count():
    """The column counts units. It is deliberately no longer a free-text field
    that could hold 'ensemble' or a whole sentence -- three encodings in one
    column is what made it unreadable."""
    for bad in ("ensemble", "6000 cells (1000 seeds x generations 5-10)"):
        with pytest.raises(Exception):
            CultivationRegistrySchema.validate(_frame(_row(n_reps=bad)))


def test_n_reps_accepts_a_numeric_string():
    """Coercion is deliberate: producers writing TSV emit '3', not 3."""
    out = CultivationRegistrySchema.validate(_frame(_row(n_reps="3")))
    assert int(out.n_reps.iloc[0]) == 3


def test_n_reps_may_be_absent_when_unknown():
    """A blank count is honest for a declared-but-unpopulated group; a
    placeholder word is not."""
    out = CultivationRegistrySchema.validate(
        _frame(_row(n_reps=pd.NA, replicate_basis=None,
                    notes="placeholder; nothing curated yet")))
    assert pd.isna(out.n_reps.iloc[0])


def test_replicate_basis_is_constrained():
    with pytest.raises(Exception):
        CultivationRegistrySchema.validate(_frame(_row(replicate_basis="bioreactor")))


def test_simulated_and_measured_bases_coexist():
    """The case the split exists for: 6000 cells and 3 reactors are both 'n',
    and the registry must be able to say they are different units."""
    out = CultivationRegistrySchema.validate(_frame(
        _row(cultivation_group_id="cd1_perf_wt_m9", mode="EXP",
             n_reps=3, replicate_basis="reactor"),
        _row(cultivation_group_id="cd1_sim_wt_m9", mode="SIM",
             n_reps=6000, replicate_basis="cell",
             notes="1000 seeds x generations 5-10"),
    ))
    bases = dict(zip(out.cultivation_group_id, out.replicate_basis))
    assert bases["cd1_perf_wt_m9"] == "reactor"
    assert bases["cd1_sim_wt_m9"] == "cell"


# --- variant / run_commit ---------------------------------------------------

def test_variant_survives_validation():
    """Regression guard for ``strict='filter'``: an UNDECLARED column is dropped
    from validate()'s return, so a producer writing variant into the TSV while
    the schema stayed silent about it would lose the value on any
    validate-then-write round trip."""
    out = CultivationRegistrySchema.validate(
        _frame(_row(mode="SIM", n_reps=4, replicate_basis="seed",
                    variant="1e-5 mM mecillinam, dAcrB")))
    assert "variant" in out.columns
    assert out.variant.iloc[0] == "1e-5 mM mecillinam, dAcrB"


def test_run_commit_survives_validation_and_may_be_blank():
    """Blank is the correct value for an archived run that recorded no git
    identity; the reason belongs in notes. What must NOT happen is the column
    silently disappearing."""
    out = CultivationRegistrySchema.validate(_frame(
        _row(cultivation_group_id="cd1_sim_a", mode="SIM", n_reps=6000,
             replicate_basis="cell", run_commit="",
             notes="run_commit blank: the run recorded no commit"),
        _row(cultivation_group_id="cd1_sim_b", mode="SIM", n_reps=512,
             replicate_basis="cell", run_commit="0b07d7d9db29f61ef17b1"),
    ))
    assert "run_commit" in out.columns
    assert out.run_commit.tolist() == ["", "0b07d7d9db29f61ef17b1"]
