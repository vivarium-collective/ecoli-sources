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

#: Detection vocabulary, shared with ``VectorObservationSchema``. A measured
#: zero and an entity that was never looked for are different claims, and
#: collapsing them is the failure this column exists to prevent.
DETECTION_STATES = ("detected", "below_limit", "not_detected")

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


def _value_present_iff_measured(df: pd.DataFrame) -> bool:
    """``value`` is non-null exactly when the sample yielded a number.

    ``not_detected`` means the entity was looked for in this sample and no
    number came back, so a value there would be an invention. Every other state
    must carry one — a null would silently drop the sample from any
    re-aggregation and change ``n``.
    """
    measured = df["detection"] != "not_detected"
    return bool((df["value"].notna() == measured).all())


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
        "kind": pa.Column(str, nullable=False),
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
        "phase": pa.Column(str, nullable=False),
        "window": pa.Column(str, nullable=False),
    },
    checks=[
        pa.Check(_value_present_iff_measured,
                 error="value must be non-null iff detection != 'not_detected'"),
        pa.Check(_replicate_key_unique,
                 error="duplicate (cultivation_group_id, observable, entity_id, units, replicate_id)"),
    ],
    strict=False,
    coerce=True,
    name="VectorReplicateSchema",
)
