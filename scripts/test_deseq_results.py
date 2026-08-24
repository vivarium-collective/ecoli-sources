"""Unit tests for the differential-expression result contract.

Exercises ``schemas.DeseqResultTableSchema`` and ``schemas.select_significant_genes``
against synthetic fixtures. The payloads that motivate this tier live in private
overlays, so the mechanism has to be testable without any of them — every
identifier below is invented.

⚠ The centre of gravity here is ``select_significant_genes``' ORDER OF
OPERATIONS, not the schema. The function exists because consumers were
re-deriving "the significant genes" independently, and the steps compose in a
way where a reasonable-looking reordering silently returns a different gene set
rather than an error. Two of the tests below fail against the plausible wrong
implementations; without them the others would pass on either.

Run with::

    uv run pytest scripts/test_deseq_results.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from schemas import DeseqResultTableSchema, select_significant_genes  # noqa: E402


def _results(rows) -> pd.DataFrame:
    return pd.DataFrame(
        rows,
        columns=["transcript_id", "gene_id", "log2_fold_change", "padj"],
    )


# --- the schema ---------------------------------------------------------------

def test_a_conforming_table_validates():
    frame = _results([
        ("tx-0001", "GENE00001", 2.5, 0.001),
        ("tx-0002", "GENE00002", -1.8, 0.02),
    ])
    DeseqResultTableSchema.validate(frame)


def test_several_transcripts_may_share_one_gene():
    """gene_id is deliberately NOT unique — that is the table's whole grain, and
    a unique constraint here would reject correct data."""
    frame = _results([
        ("tx-0001", "GENE00001", 2.5, 0.001),
        ("tx-0002", "GENE00001", 1.2, 0.04),
    ])
    DeseqResultTableSchema.validate(frame)


def test_a_transcript_with_no_gene_and_no_estimate_is_representable():
    """A heterologous transcript has no reference gene, and a tested transcript
    can yield no estimate. Both are real rows, not errors: independent filtering
    leaving padj undefined is meaningfully different from 'not significant'."""
    frame = _results([("tx-0003", None, None, None)])
    DeseqResultTableSchema.validate(frame)


def test_a_repeated_transcript_is_rejected():
    """transcript_id IS the grain, so duplicates mean the file is malformed."""
    frame = _results([
        ("tx-0001", "GENE00001", 2.5, 0.001),
        ("tx-0001", "GENE00002", 1.0, 0.01),
    ])
    with pytest.raises(Exception):
        DeseqResultTableSchema.validate(frame)


def test_an_out_of_range_padj_is_rejected():
    frame = _results([("tx-0001", "GENE00001", 2.5, 1.4)])
    with pytest.raises(Exception):
        DeseqResultTableSchema.validate(frame)


# --- the canonical selection --------------------------------------------------

def test_selection_keeps_only_significant_genes():
    frame = _results([
        ("tx-0001", "GENE00001", 2.5, 0.001),   # significant
        ("tx-0002", "GENE00002", 0.4, 0.001),   # |lfc| too small
        ("tx-0003", "GENE00003", 3.0, 0.5),     # padj too high
    ])
    assert list(select_significant_genes(frame).index) == ["GENE00001"]


def test_thresholds_are_strict_on_both_sides():
    """Exactly-at-threshold rows are excluded, matching `padj < t` and `|lfc| > t`."""
    frame = _results([
        ("tx-0001", "GENE00001", 1.0, 0.05),    # |lfc| == thresh -> out
        ("tx-0002", "GENE00002", 2.0, 0.1),     # padj == thresh -> out
    ])
    assert select_significant_genes(frame).empty


def test_transcripts_collapse_to_the_largest_absolute_change():
    frame = _results([
        ("tx-0001", "GENE00001", 1.5, 0.001),
        ("tx-0002", "GENE00001", -4.0, 0.002),  # larger |lfc|, negative
        ("tx-0003", "GENE00001", 2.0, 0.003),
    ])
    out = select_significant_genes(frame)
    assert list(out.index) == ["GENE00001"]
    assert out.loc["GENE00001", "log2_fold_change"] == -4.0
    assert out.loc["GENE00001", "transcript_id"] == "tx-0002"


def test_filtering_happens_BEFORE_collapsing():
    """⚠ The discriminating case for order of operations.

    GENE00001 has a large-|lfc| transcript that is NOT significant, and a
    smaller-|lfc| one that IS. Filter-then-collapse keeps the gene via the
    significant transcript. Collapse-then-filter would pick the large
    non-significant transcript first and then drop the gene entirely.

    Both implementations are one line apart and neither errors — only the gene
    set differs. Without this test the suite passes on either.
    """
    frame = _results([
        ("tx-0001", "GENE00001", 9.0, 0.9),     # biggest |lfc|, NOT significant
        ("tx-0002", "GENE00001", 1.5, 0.001),   # significant
    ])
    out = select_significant_genes(frame)
    assert list(out.index) == ["GENE00001"], "gene lost — collapse ran before the filter"
    assert out.loc["GENE00001", "transcript_id"] == "tx-0002"


def test_rows_without_a_gene_id_are_dropped_not_grouped():
    """⚠ Second discriminating case. Dropping un-mapped transcripts is step 1 of
    the contract, not a tidy-up: a significant transcript with no gene_id must
    not survive as a null-keyed row, which is what a plain groupby would do.
    """
    frame = _results([
        ("tx-0001", "GENE00001", 2.5, 0.001),
        ("tx-0002", None, 5.0, 0.001),          # significant, but maps to no gene
    ])
    out = select_significant_genes(frame)
    assert list(out.index) == ["GENE00001"]
    assert out.index.notna().all()


def test_selection_does_not_mutate_its_input():
    """The implementation adds a working column; the caller must not see it."""
    frame = _results([("tx-0001", "GENE00001", 2.5, 0.001)])
    before = frame.copy(deep=True)
    select_significant_genes(frame)
    pd.testing.assert_frame_equal(frame, before)
