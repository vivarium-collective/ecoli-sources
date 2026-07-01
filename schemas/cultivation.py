"""
Pandera schemas for the cultivation-centric validation layer.

These extend the validation-data subsystem (``schemas/validation_claim.py``)
toward a cultivation-keyed, multi-performer organization:

* ``CultivationRegistrySchema`` — the **cultivation registry**
  (``validation_data/cultivations.tsv``). One row per cultivation = a specific
  strain × medium × run, carrying the who/how/when/why. A SIM is a cultivation
  too (``mode = SIM``), so the registry is the EXP↔SIM join. Assay tables
  reference it by ``cultivation_id``.

* ``ScalarObservationSchema`` — a **tidy** scalar-observation assay table
  (e.g. ``validation_data/fermentation.tsv``). One row per
  ``(cultivation_id, observable, replicate)``. The cultivation-keyed,
  tidy-long generalization of ``ScalarClaimSchema`` (which is per-observable
  file, single implicit cultivation): many observables and cultivations live in
  one assay table, joined to the registry by ``cultivation_id`` and addressed by
  ``(cultivation_id, observable)``. This tidy, long-form layout is a common
  assay-export shape, so one representation serves grading and data delivery.

Both are read by ``ecoli_sources.validation`` (the overlay-union loader), which
unions a primary validation bundle with any overlay bundles on
``canonical_key`` and resolves each scalar key to its measured band — from a
per-observable ``ScalarClaimSchema`` file OR a tidy ``ScalarObservationSchema``
table filtered to ``(cultivation_id, observable)``.
"""

import pandera.pandas as pa

# Mirrors ``schemas.validation_claim._KINDS`` (kept local to avoid a private
# cross-module import): measured = grade target; theoretical_max = first-
# principles ceiling; model_predicted = another model's value (context).
_KINDS = ["measured", "theoretical_max", "model_predicted"]

# A cultivation is experimental, simulated, or a literature value. The registry
# is the join across all three (e.g. grade a SIM cultivation vs an EXP one).
_MODES = ["EXP", "SIM", "literature"]


# ---------------------------------------------------------------------------
# Cultivation registry: one row per cultivation (the who/how/when/why)
# ---------------------------------------------------------------------------

CultivationRegistrySchema = pa.DataFrameSchema(
    name="cultivation_registry",
    columns={
        "cultivation_id": pa.Column(
            dtype=str,
            unique=True,
            nullable=False,
            description=(
                "Primary key (snake_case). Stable token identifying one "
                "cultivation; assay tables and validation-bundle rows reference "
                "it. Convention: ``<cd_track>_<performer>_<strain>_<medium>``; "
                "full identity lives in the columns."
            ),
        ),
        "mode": pa.Column(
            dtype=str,
            nullable=False,
            checks=pa.Check.isin(_MODES),
            description="One of EXP | SIM | literature. Lets the registry join sim vs experiment.",
        ),
        "performer": pa.Column(
            dtype=str, nullable=True, required=False,
            description="Who executed the cultivation (an experimental lab, or the simulator for a SIM).",
        ),
        "platform": pa.Column(
            dtype=str, nullable=True, required=False,
            description="Instrument/platform (e.g. Ambr250, bioreactor, vEcoli).",
        ),
        "experiment_id": pa.Column(
            dtype=str, nullable=True, required=False,
            description="Performer experiment id, or a sim run/commit ref.",
        ),
        "date": pa.Column(
            dtype=str, nullable=True, required=False,
            description="Cultivation date (free-form; e.g. 2025-11). Approximate values flagged in notes.",
        ),
        "cd_track": pa.Column(
            dtype=str, nullable=True, required=False,
            description="Campaign / demonstration-track label, or '—' for none.",
        ),
        "strain": pa.Column(
            dtype=str, nullable=True, required=False,
            description="Strain (e.g. MG1655 rph+).",
        ),
        "genotype": pa.Column(
            dtype=str, nullable=True, required=False,
            description="Detailed genotype, including heterologous content.",
        ),
        "medium": pa.Column(
            dtype=str, nullable=True, required=False,
            description="Growth medium / condition (e.g. Modified M9+N+Fe).",
        ),
        "operation": pa.Column(
            dtype=str, nullable=True, required=False,
            description="Operation mode (e.g. batch, fed-batch). Inferred values flagged in notes.",
        ),
        "n_reps": pa.Column(
            dtype=str, nullable=True, required=False,
            description="Replicate count (e.g. 3) or descriptor (e.g. ensemble). String to allow both.",
        ),
        "notes": pa.Column(
            dtype=str, nullable=True, required=False,
            description="Free-text who/how/when/why + provisional-metadata caveats.",
        ),
    },
    strict="filter",  # allow extra columns; validate the required ones
    coerce=True,
    description=(
        "Cultivation registry: one row per cultivation (strain × medium × run), "
        "carrying provenance. A SIM is a cultivation (mode=SIM), so this is the "
        "EXP↔SIM join. Assay tables reference rows by cultivation_id."
    ),
)


# ---------------------------------------------------------------------------
# Tidy scalar observations: one row per (cultivation, observable, replicate)
# ---------------------------------------------------------------------------

ScalarObservationSchema = pa.DataFrameSchema(
    name="scalar_observation",
    columns={
        "cultivation_id": pa.Column(
            dtype=str,
            nullable=False,
            description="FK to the cultivation registry (cultivations.tsv).",
        ),
        "observable": pa.Column(
            dtype=str,
            nullable=False,
            description=(
                "The measured quantity (e.g. growth_rate, glucose_uptake, "
                "acetate_secretion). With cultivation_id, addresses a validation "
                "slot: the bundle's (cultivation_id, observable) filter selects "
                "the rows for one canonical_key."
            ),
        ),
        "value": pa.Column(
            float, nullable=False,
            description="The observed scalar value, expressed in ``units``.",
        ),
        "units": pa.Column(
            dtype=str, nullable=False,
            description="Units of ``value`` (e.g. mmol/gDW/h, 1/h, gDW/g_glucose).",
        ),
        "kind": pa.Column(
            dtype=str, nullable=False, checks=pa.Check.isin(_KINDS),
            description="One of measured | theoretical_max | model_predicted.",
        ),
        "phase": pa.Column(
            dtype=str, nullable=True, required=False,
            description=(
                "Derived-tier temporal provenance: the biological interval the KPI "
                "summarizes (e.g. exponential_batch, end_of_fermentation). A KPI is "
                "computed OVER an interval (with judgment), not sampled at a point — "
                "so this, not a raw culture age, is the right temporal descriptor for "
                "the derived/summary tier (the tier this table lives in, like the "
                "literature references)."
            ),
        ),
        "window": pa.Column(
            dtype=str, nullable=True, required=False,
            description=(
                "Concrete interval for `phase` (e.g. '3-4h', 'EOF'). The judgment — "
                "how the interval was determined — lives in the upstream extractor; "
                "this records the result. Blank where not yet captured."
            ),
        ),
        "replicate": pa.Column(
            dtype=str, nullable=True, required=False,
            description=(
                "Replicate label within the cultivation (e.g. reactor r1/r2/r3). "
                "Multiple measured replicates form the measured band."
            ),
        ),
        "method": pa.Column(
            dtype=str, nullable=True, required=False,
            description="How the value was obtained (e.g. Ambr250 fermentation off-gas).",
        ),
        "uncertainty": pa.Column(
            float, nullable=True, required=False,
            checks=pa.Check.greater_than_or_equal_to(0),
            description="Optional symmetric ± uncertainty, same units as ``value``.",
        ),
        "notes": pa.Column(
            dtype=str, nullable=True, required=False,
            description="Free-text caveats (units assumptions, sampling phase, etc.).",
        ),
    },
    strict="filter",  # allow extra columns; validate the required ones
    coerce=True,
    description=(
        "Tidy scalar-observation assay table; one row per (cultivation_id, "
        "observable, replicate). The cultivation-keyed, tidy-long generalization "
        "of ScalarClaimSchema, a common tidy assay-export layout. This is the "
        "DERIVED/summary tier (KPIs computed over an interval, like the literature "
        "references) — temporal provenance is phase/window + method, not culture age."
    ),
)
