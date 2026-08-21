"""Per-sample sibling of ``VectorObservationSchema``.

``VectorObservationSchema`` ships one row per entity, pre-aggregated over an
interval. This ships the **individual measurements behind it** — one row per
``(cultivation_group_id, observable, entity_id, units, replicate_id)``.

**Both tiers are curated the same way, by the same producer.** The aggregation's
ownership is about **which samples constitute the condition** — the interval, the
replicates, and the exclusions that reasoning implies — **not** about whether
those samples are averaged before shipping. Averaging was an implementation of
that decision, not the decision itself.

⇒ A per-replicate table holds *exactly* the samples its aggregated sibling
averages: no more and no fewer. Re-aggregating it reproduces that sibling's
``n`` / ``n_pos`` / ``n_total`` exactly and its means to emitter rounding. A
consumer therefore cannot widen the interval or re-admit an excluded replicate
from this file — those samples are not in it. **The guarantee is mechanical
rather than promised**, which is what keeps the curation upstream.

**What this tier is for.** Statistics that need a *distribution* rather than its
summary. A kernel two-sample test compares distributions; a mean and an ``n``
cannot be sampled. ⚠ Synthesizing samples from a mean and a dispersion is not an
adequate substitute when the other side of a comparison carries real ones: drawn
independently per entity, the synthetic side has an identity covariance
structure by construction, so the two sides differ in a way that has nothing to
do with what was measured.

**What it is NOT for.** It is not a raw tier and must not be re-selected from,
pooled with samples outside the interval, or substituted for the aggregate on
axes that grade against the aggregate. The aggregated sibling remains the
default operand.

⚠ **``replicate_id`` is not necessarily a cultivation registry's replicate.** A
registry counts replicates on a declared basis (``replicate_basis``); a row here
is whatever the producer's provenance record says one row is. Where an interval
pools several timepoints, one row is a *sample* — a replicate observed at a
timepoint — not a replicate. ``replicate`` and ``timepoint_h`` are therefore
carried separately and **nothing is pre-collapsed**, so a consumer can reduce to
whichever level it needs. Producers state in prose, in the provenance record,
what one row is.

⚠ **``time_basis`` exists because the two sides of a comparison do not relate to
time the same way.** A measured row is typically a sample drawn at
``timepoint_h``. A simulated per-cell row has no such instant: its value is a
mean over that cell's cycle. Both are legitimate, and a null ``timepoint_h``
alone cannot tell them apart — so the aggregation is *declared* rather than
inferred from an absence.

**A bare ``value`` here does not contradict the aggregated tier's refusal of
one.** That refusal was about a ``value`` sitting *beside* ``mean_arithmetic`` /
``mean_geometric`` alternatives, where a consumer can silently pick the wrong
one. A single sample has one number and no alternative.

**Deliberately narrow.** ``symbol``, ``method`` and ``notes`` are not repeated
per sample — the aggregated sibling carries them over an identical entity set,
and duplicating a display alias once per sample buys nothing at this row count.
"""

from __future__ import annotations

import pandas as pd
import pandera.pandas as pa

from .cultivation import _REPLICATE_BASES
from .vector_observation import _DETECTION, _KINDS

#: Detection vocabulary. ⚠ IMPORTED from ``vector_observation`` rather than
#: restated: a per-sample row and its aggregate must agree on what a token
#: MEANS, and two independent lists are free to drift into disagreeing while
#: both validate. (An earlier revision of this module declared its own copy
#: under a comment claiming it was shared. It was not.)
DETECTION_STATES = tuple(_DETECTION)

#: How a row's value relates to time. ``instant`` is a sample taken at
#: ``timepoint_h``; ``interval_mean`` is an average over an interval, which is
#: what a per-cell simulated value is (a mean over that cell's cycle). Required,
#: because leaving it to a null ``timepoint_h`` would make "aggregated over an
#: interval" indistinguishable from "timepoint unknown" — the same conflation
#: ``detection`` exists to prevent one column over.
TIME_BASES = ("instant", "interval_mean")

VECTOR_REPLICATE_COLUMNS = [
    "cultivation_group_id",
    "observable",
    "entity_id",
    "units",
    "kind",
    "replicate_id",
    "replicate",
    "replicate_basis",
    "time_basis",
    "timepoint_h",
    "sample_id",
    "detection",
    "value",
    "phase",
    "window",
]


def _value_present_iff_detected(df: pd.DataFrame) -> bool:
    """``value`` is present exactly when the sample DETECTED the entity.

    ⚠ Mirrors ``VectorObservationSchema``'s rule for ``mean_arithmetic``
    deliberately, and an earlier revision of this module did not: it required a
    number on ``below_limit`` rows while the aggregate refused one, so the two
    tiers asserted OPPOSITE contracts for the same token.

    That divergence was not cosmetic. A ``below_limit`` row is a statement about
    the limit, not a measurement of zero — and a producer emitting ``0.0`` there
    hands a re-aggregating consumer a number that drags the mean toward zero and
    inflates ``n``. Measured on a worked case: values (10, 0, 20) with the zero
    reported below limit give the aggregate ``n=2, mean=15.0``, while naively
    averaging a shipped ``0.0`` gives ``n=3, mean=10.0``. Silently.

    ⇒ ``below_limit`` and ``not_detected`` both carry a null. What distinguishes
    them is the CLAIM (looked for and under the limit, versus looked for and
    absent), which is what ``detection`` is for.
    """
    return bool(((df["detection"] == "detected") == df["value"].notna()).all())


def _instant_rows_carry_a_timepoint(df: pd.DataFrame) -> bool:
    """An ``instant`` row must say WHEN; an ``interval_mean`` row need not.

    Without this the column does not deliver its own stated purpose: ``instant``
    with a null ``timepoint_h`` is exactly the "is this unknown, or aggregated?"
    ambiguity ``time_basis`` exists to remove, and it validated.

    ``interval_mean`` is left free: an interval may legitimately carry a
    midpoint, and the interval itself is named by ``phase``/``window``.

    ⚠ Known limit, accepted 2026-08-21: this assumes every ``instant``-basis
    modality can state a timepoint. A cultivation modality where a sample is
    genuinely instantaneous but untimed would fail here and should relax this
    rule rather than mislabel itself ``interval_mean``.
    """
    instant = df["time_basis"] == "instant"
    return bool(df.loc[instant, "timepoint_h"].notna().all())


def _replicate_key_unique(df: pd.DataFrame) -> bool:
    """One row per (group, observable, entity, units, replicate).

    A duplicate here silently double-weights a sample in any re-aggregation,
    which is invisible in the mean and wrong in the counts.
    """
    key = ["cultivation_group_id", "observable", "entity_id", "units", "replicate_id"]
    return not bool(df.duplicated(subset=key).any())


VectorReplicateSchema = pa.DataFrameSchema(
    {
        "cultivation_group_id": pa.Column(str, nullable=False),
        "observable": pa.Column(str, nullable=False),
        "entity_id": pa.Column(str, nullable=False),
        "units": pa.Column(str, nullable=False),
        "kind": pa.Column(
            str, nullable=False, checks=pa.Check.isin(_KINDS),
            description=(
                "Same vocabulary as the aggregated tier. Unconstrained, a "
                "plausible-looking ``simulated`` would validate here and then "
                "fall out of every ``kind == 'measured'`` filter downstream."),
        ),
        "replicate_id": pa.Column(
            str, nullable=False,
            description=(
                "Identifier for one sample within the curated interval. NOT "
                "necessarily the cultivation registry's replicate — see the "
                "module docstring."),
        ),
        "replicate": pa.Column(
            str, nullable=True,
            description=(
                "The replicate this sample came from, on the basis named in "
                "``replicate_basis`` — a reactor, a well, a simulated cell. "
                "Deliberately NOT named for one basis: the same shape serves a "
                "measured operand (basis ``reactor``) and a simulated one "
                "(basis ``cell``)."),
        ),
        "replicate_basis": pa.Column(
            str, nullable=True,
            checks=pa.Check.isin(_REPLICATE_BASES),
            description=(
                "What one replicate IS, from the cultivation registry's "
                "vocabulary. Carried per row so the table is self-describing: "
                "``n`` alone cannot tell a consumer whether it is looking at "
                "3 bioreactors or 6000 simulated cells, and a statistic that "
                "treats those alike is wrong."),
        ),
        "time_basis": pa.Column(
            str, nullable=False,
            checks=pa.Check.isin(TIME_BASES),
            description=(
                "Whether the value is a sample at ``timepoint_h`` "
                "(``instant``) or a mean over an interval "
                "(``interval_mean``). Required so that a null ``timepoint_h`` "
                "cannot silently mean two different things."),
        ),
        "timepoint_h": pa.Column(
            float, nullable=True,
            description=(
                "Hours into the cultivation, for an ``instant`` row. Null when "
                "``time_basis`` is ``interval_mean`` — the interval is then "
                "named by ``phase``/``window``. Carried separately from "
                "``replicate`` so a consumer can collapse to either level."),
        ),
        "sample_id": pa.Column(
            str, nullable=True,
            description="Producer's own assay-level sample identifier, for traceback.",
        ),
        "detection": pa.Column(
            str, nullable=False,
            checks=pa.Check.isin(DETECTION_STATES),
            description="Per-sample detection state; same vocabulary as the aggregate.",
        ),
        "value": pa.Column(
            float, nullable=True,
            checks=pa.Check.greater_than_or_equal_to(0),
            description=(
                "The measurement, in ``units``. Null iff ``detection`` is "
                "``not_detected``. Non-negative: this tier carries abundances, "
                "not signed quantities."),
        ),
        # Optional and nullable, matching the aggregated tier. Required here
        # only, a producer could ship an aggregate with no phase/window but not
        # its replicate sibling — an asymmetry with no reason behind it.
        "phase": pa.Column(str, nullable=True, required=False),
        "window": pa.Column(str, nullable=True, required=False),
    },
    checks=[
        pa.Check(_value_present_iff_detected,
                 error="value must be non-null iff detection == 'detected'"),
        pa.Check(_instant_rows_carry_a_timepoint,
                 error="time_basis 'instant' requires a timepoint_h"),
        pa.Check(_replicate_key_unique,
                 error="duplicate (cultivation_group_id, observable, entity_id, units, replicate_id)"),
    ],
    strict="filter",  # allow extra columns; validate the required ones
    coerce=True,
    name="VectorReplicateSchema",
)
