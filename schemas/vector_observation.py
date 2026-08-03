"""
Pandera schema for **tidy vector observations** — a cultivation-keyed assay
table whose observable is a vector rather than a scalar.

Completes a 2×2 that was previously missing a cell:

|            | per-source / literature                    | cultivation-keyed / assay   |
|------------|--------------------------------------------|-----------------------------|
| **scalar** | ``ScalarClaimSchema``                      | ``ScalarObservationSchema`` |
| **vector** | ``ProteinAbundanceSchema``, ``ReactionFluxSchema``, … | ``VectorObservationSchema`` |

``ScalarObservationSchema`` is the cultivation-keyed generalization of
``ScalarClaimSchema``; this is the same move one row down. One row per
``(cultivation_group_id, observable, entity_id, units)`` — e.g. one measured
protein abundance for one cultivation group's exponential-batch aggregate.

**Measurement-only, by design.** A vector observation carries *no* model-side
value and no simulation provenance. Pairing a measurement with a simulation —
resolving the shared entity set, renormalizing, computing concordance, and
rendering two-sided provenance — belongs to whatever grades the comparison, not
to the curated data. Baking a model vector in alongside makes the measurement's
own content model-dependent: the entity set becomes an intersection with that
model's id-space, and any normalization computed over the intersection shifts
when model coverage shifts, though nothing was measured differently.

**Aggregated, and deliberately so.** Unlike the scalar tier — which keeps one row
per replicate and lets the consumer aggregate — a vector observation ships
pre-aggregated over an interval named by ``phase``/``window``. That aggregation
is a *curation* decision: choosing which timepoints represent a condition can
depend on knowledge absent from the data (replicate clustering, a process
deviation that disqualifies an interval). A downstream aggregator lacking that
knowledge would produce a worse summary, not a more neutral one. Multiple
aggregations of one experiment coexist as rows differing in ``phase``/``window``;
``cultivation_group_id`` stays strictly physical.

Because the values are aggregates, each row carries its own dispersion and
counts, so a consumer can weight, filter, or draw a band without the raw tier.
"""

from __future__ import annotations

import pandera.pandas as pa

# Mirrors ``schemas.validation_claim._KINDS``: the epistemic class of the claim.
_KINDS = ["measured", "theoretical_max", "model_predicted"]

# The DETECTION outcome — a different axis from ``kind``. Deliberately worded to
# avoid a collision: a ``kind`` of "measured" and a detection of "measured" would
# be the same word for two unrelated things.
#
# * ``detected``     — a real result. ``0.0`` here is a genuine measured zero
#                      (e.g. an RNA-seq panel gene with no reads mapped).
# * ``below_limit``  — looked for, below the limit of detection/quantification.
#                      A statement about the limit, NOT a measurement of zero.
# * ``not_detected`` — on the panel, absent from this group's measurements.
_DETECTION = ["detected", "below_limit", "not_detected"]


def _value_matches_detection(df) -> bool:
    """``mean_arithmetic`` is present exactly when the entity was detected.
    A non-detection with a number, or a detection without one, is incoherent."""
    return bool(((df["detection"] == "detected") == df["mean_arithmetic"].notna()).all())


def _counts_are_nested(df) -> bool:
    """``n_pos`` ⊆ ``n`` ⊆ ``n_total`` — positives are a subset of the
    measurements contributing, which are a subset of those available."""
    return bool(((df["n_pos"] <= df["n"]) & (df["n"] <= df["n_total"])).all())


def _log_statistics_have_support(df) -> bool:
    """A geometric mean needs ≥1 positive value; a log-scale SD needs ≥2.
    (Both are computed over positives — zeros cannot be logged.)"""
    gm_ok = bool((df["mean_geometric"].notna() <= (df["n_pos"] >= 1)).all())
    sd_ok = bool((df["sd_log10"].notna() <= (df["n_pos"] >= 2)).all())
    return gm_ok and sd_ok


VectorObservationSchema = pa.DataFrameSchema(
    name="vector_observation",
    columns={
        "cultivation_group_id": pa.Column(
            dtype=str, nullable=False,
            description=(
                "FK to the cultivation-group registry (``cultivations.tsv``). The "
                "DERIVED group — strain × medium × run, aggregating replicate "
                "reactors — not a raw, LabKey-atomic per-reactor cultivation. "
                "Strictly physical: an aggregation window is NOT part of this key "
                "(see ``phase``/``window``)."
            ),
        ),
        "observable": pa.Column(
            dtype=str, nullable=False,
            description=(
                "The vector-valued quantity (e.g. proteome, transcriptome). With "
                "``cultivation_group_id`` and the ``phase``/``window`` pair it "
                "addresses one validation slot."
            ),
        ),
        "entity_id": pa.Column(
            dtype=str, nullable=False,
            description=(
                "Within-vector key — the entity as the source identifies it (an "
                "EcoCyc monomer id for a proteome, a gene id for a transcriptome). "
                "Modality-neutral by design, unlike ``ProteinAbundanceSchema.gene`` "
                "or ``ReactionFluxSchema.reaction_id``."
            ),
        ),
        "symbol": pa.Column(
            dtype=str, nullable=True, required=False,
            description=(
                "Human-readable alias for ``entity_id`` (e.g. a gene symbol). "
                "Display and labelling only — never a join key: symbol spaces "
                "differ between sources and are not injective."
            ),
        ),
        "units": pa.Column(
            dtype=str, nullable=False,
            description=(
                "Units of the statistics on this row. Part of the row key: one "
                "measurement reported in two complementary units (e.g. a "
                "detection-normalized abundance and a mass fraction) is two rows, "
                "not two columns — which keeps the statistic columns fixed as "
                "units are added."
            ),
        ),
        "kind": pa.Column(
            dtype=str, nullable=False, checks=pa.Check.isin(_KINDS),
            description=(
                "Epistemic class of the claim: measured | theoretical_max | "
                "model_predicted. Distinct from ``detection``, which is the "
                "assay outcome."
            ),
        ),
        "detection": pa.Column(
            dtype=str, nullable=False, checks=pa.Check.isin(_DETECTION),
            description=(
                "detected | below_limit | not_detected. Separates a genuine "
                "measured zero from an absence — a distinction a bare number "
                "cannot carry, and one that silently disappears if absence is "
                "encoded as row-omission or as a null the reader coerces."
            ),
        ),
        "mean_arithmetic": pa.Column(
            float, nullable=True,
            checks=pa.Check.greater_than_or_equal_to(0),
            description=(
                "Mean over the detected measurements in the window, INCLUDING "
                "genuine zeros. The additive summary: for a unit that is a "
                "fraction of a whole, only this mean preserves the sum. Null "
                "exactly when ``detection != detected``."
            ),
        ),
        "mean_geometric": pa.Column(
            float, nullable=True,
            checks=pa.Check.greater_than_or_equal_to(0),
            description=(
                "Mean over POSITIVE measurements only. The log-space centre — the "
                "right estimate when the quantity is ~log-normal and compared on "
                "log axes, and the only mean coherent with ``sd_log10``. Null when "
                "``n_pos`` is 0."
            ),
        ),
        "sd_log10": pa.Column(
            float, nullable=True,
            checks=pa.Check.greater_than_or_equal_to(0),
            description=(
                "Dispersion of log10 values over positive measurements. On the log "
                "scale rather than linear because a band drawn on a log axis needs "
                "it, and a geometric dispersion cannot be recovered from a linear "
                "one without the raw values. Null when ``n_pos`` < 2 — routine, not "
                "exceptional, for sparsely-detected entities."
            ),
        ),
        "n": pa.Column(
            int, nullable=False, checks=pa.Check.greater_than_or_equal_to(0),
            description="Measurements contributing to ``mean_arithmetic``.",
        ),
        "n_pos": pa.Column(
            int, nullable=False, checks=pa.Check.greater_than_or_equal_to(0),
            description="Positive measurements behind ``mean_geometric`` / ``sd_log10``.",
        ),
        "n_total": pa.Column(
            int, nullable=False, checks=pa.Check.greater_than_or_equal_to(0),
            description=(
                "Measurements AVAILABLE for this cultivation group × window "
                "(replicates × timepoints). The denominator: a mean over 4 of 6 "
                "samples and a mean over 6 of 6 are different claims, and "
                "'detected in most samples' filtering is impossible without it."
            ),
        ),
        "phase": pa.Column(
            dtype=str, nullable=True, required=False,
            description=(
                "The biological interval this aggregate summarizes (e.g. "
                "exponential_batch). Same derived-tier vocabulary as "
                "``ScalarObservationSchema``."
            ),
        ),
        "window": pa.Column(
            dtype=str, nullable=True, required=False,
            description=(
                "Concrete interval for ``phase`` (e.g. 3-4h). The judgment — how "
                "the interval was chosen — lives upstream in the extractor and its "
                "narrative summary; this records the result."
            ),
        ),
        "method": pa.Column(
            dtype=str, nullable=True, required=False,
            description="How the values were obtained (e.g. DIA proteomics, RNA-seq).",
        ),
        "notes": pa.Column(
            dtype=str, nullable=True, required=False,
            description="Free-text caveats (unit assumptions, detection semantics, …).",
        ),
    },
    checks=[
        pa.Check(_value_matches_detection, name="value_matches_detection",
                 error="mean_arithmetic must be non-null exactly when detection == 'detected'"),
        pa.Check(_counts_are_nested, name="counts_are_nested",
                 error="counts must satisfy n_pos <= n <= n_total"),
        pa.Check(_log_statistics_have_support, name="log_statistics_have_support",
                 error="mean_geometric requires n_pos >= 1; sd_log10 requires n_pos >= 2"),
    ],
    unique=["cultivation_group_id", "observable", "entity_id", "units"],
    strict="filter",  # allow extra columns; validate the required ones
    coerce=True,
    description=(
        "Tidy vector-observation assay table: one row per (cultivation_group_id, "
        "observable, entity_id, units), carrying a pre-aggregated measurement with "
        "its dispersion, counts, and detection status. The cultivation-keyed "
        "analogue of the per-source vector claim schemas, and the vector analogue "
        "of ScalarObservationSchema. Measurement-only — no model-side value."
    ),
)
