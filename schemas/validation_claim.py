"""
Pandera schemas for the validation-data subsystem.

The validation-data subsystem (``ecoli_sources/validation_data/``) is a
SEPARATE concern from the ParCa reference bundle (``data/reference_bundle.tsv``).
The reference bundle ships data fed *into* the model; validation data is
curated experimental and theoretical reference values that model *outputs*
are graded *against*. See ``ecoli_sources/validation_data/README.md`` for the
subsystem design and rationale.

Two schemas:

* ``ValidationBundleSchema`` — the validation manifest
  (``validation_data/validation_bundle.tsv``). Same
  ``(canonical_key, source_path, description, schema_name)`` shape as the
  reference bundle, but deliberately WITHOUT a required-canonical-keys
  contract: a missing validation key should make a report-card axis
  ``ungraded``, not raise a fatal load error the way a missing ParCa input
  does.
* ``ScalarClaimSchema`` — a curated table of scalar reference claims for one
  ``(condition, observable)`` slot, e.g. ``data/basal/biomass_yield.tsv``.
  One row per source; multiple sources (measured + theoretical) per observable
  are the norm, distinguished by ``kind``.
* ``ReactionFluxSchema`` — a curated *vector* reference: per-reaction fluxes
  for one ``(condition, observable)`` slot, e.g. ``data/basal/metabolic_fluxes.tsv``.
  One row per ``(source, reaction, timepoint)``; the within-vector key is
  ``reaction_id``. This is the vector analogue of ``ScalarClaimSchema`` for
  ¹³C-MFA flux maps and similar per-reaction datasets, where a single source
  contributes many addressable reactions rather than one scalar.
"""

import pandera.pandas as pa


# ---------------------------------------------------------------------------
# Validation manifest: canonical_key -> claim table (NO required-keys contract)
# ---------------------------------------------------------------------------

ValidationBundleSchema = pa.DataFrameSchema(
    name="validation_bundle",
    columns={
        "canonical_key": pa.Column(
            dtype=str,
            unique=True,
            nullable=False,
            description=(
                "Slot name for a (condition, observable) validation role "
                "(snake_case, primary key). Convention: ``<condition>__<observable>`` "
                "(e.g. ``basal__biomass_yield``)."
            ),
        ),
        "source_path": pa.Column(
            dtype=str,
            nullable=False,
            description=(
                "Path to the curated claim table, relative to the "
                "``validation_data/`` directory (e.g. ``data/basal/biomass_yield.tsv``)."
            ),
        ),
        "description": pa.Column(
            dtype=str,
            nullable=False,
            description="One-liner explaining the validation slot's role.",
        ),
        "schema_name": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description=(
                "Optional Pandera schema for the claim table, resolved from "
                "``ecoli_sources/schemas/`` (e.g. ``ScalarClaimSchema``)."
            ),
        ),
    },
    # Unlike ReferenceBundleSchema, NO ``_check_required_canonical_keys`` check:
    # validation slots are optional by design. A consumer that asks for an
    # absent key degrades that report-card axis to ``ungraded``.
    strict="filter",
    coerce=True,
    description=(
        "Validation manifest mapping (condition, observable) keys to curated "
        "claim tables. Distinct from the ParCa reference bundle: carries NO "
        "required-canonical-keys contract — a missing key degrades a "
        "report-card axis to 'ungraded' rather than failing the load."
    ),
)


# ---------------------------------------------------------------------------
# Scalar claim table: one row per source for a single (condition, observable)
# ---------------------------------------------------------------------------

# ``kind`` separates the three reference classes a card treats differently:
#   - measured        : an experimental observation; the grade target.
#   - theoretical_max  : a first-principles ceiling (thermodynamic or
#                        stoichiometric / FBA). A model OUTPUT that EXCEEDS a
#                        theoretical_max is unambiguously broken — a harder,
#                        differentiated failure than deviating from a measured
#                        value (no adequacy judgment required).
#   - model_predicted  : a value from another model (context / cross-reference).
_KINDS = ["measured", "theoretical_max", "model_predicted"]

ScalarClaimSchema = pa.DataFrameSchema(
    name="scalar_claim",
    columns={
        "source_id": pa.Column(
            dtype=str,
            nullable=False,
            description=(
                "Provenance key. Must match a directory under "
                "``validation_data/references/<source_id>/`` and (ideally) a "
                "BibTeX entry in ``validation_data/references.bib``."
            ),
        ),
        "value": pa.Column(
            float,
            nullable=False,
            description="The reported scalar value, expressed in ``units``.",
        ),
        "units": pa.Column(
            dtype=str,
            nullable=False,
            description="Units of ``value`` (e.g. ``gDW/g_glucose``, ``1/h``).",
        ),
        "kind": pa.Column(
            dtype=str,
            nullable=False,
            checks=pa.Check.isin(_KINDS),
            description=(
                "One of measured | theoretical_max | model_predicted. Cards "
                "grade against ``measured`` and treat exceedance of a "
                "``theoretical_max`` as a harder, differentiated failure."
            ),
        ),
        "method": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description=(
                "How the value was obtained (e.g. batch, chemostat, 13C-MFA, "
                "thermodynamic, FBA)."
            ),
        ),
        "uncertainty": pa.Column(
            float,
            nullable=True,
            required=False,
            checks=pa.Check.greater_than_or_equal_to(0),
            description="Optional symmetric ± uncertainty, same units as ``value``.",
        ),
        "strain": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description="Strain background (e.g. MG1655), if reported.",
        ),
        "notes": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description="Free-text caveats: media, growth mode, extraction notes.",
        ),
    },
    strict="filter",  # allow extra columns; validate the required ones
    coerce=True,
    description=(
        "Curated scalar reference claims for one (condition, observable) slot; "
        "one row per source. Multiple sources (measured + theoretical) per "
        "observable are expected, distinguished by ``kind``."
    ),
)


# ---------------------------------------------------------------------------
# Reaction-flux claim table: a VECTOR reference (per-reaction fluxes)
# ---------------------------------------------------------------------------

# The vector analogue of ScalarClaimSchema. Where a scalar slot holds one value
# per source, a flux slot holds a whole reaction map per source: one row per
# ``(source_id, reaction_id, timepoint)``. The within-vector key is
# ``reaction_id``. Multiple sources accrue into one table (e.g. several
# ¹³C-MFA studies' central-carbon maps), exactly as multiple sources stack in a
# scalar slot.
#
# Columns: ``value`` is the absolute flux (model-comparable units);
# ``value_relative_pct`` is the flux normalized to a reference flux (e.g. glucose
# uptake = 100), provided for scale-invariant comparison of carbon routing.
# ``timepoint`` is carried when a source reports a time course (e.g. nonstationary
# MFA); a time course is a trajectory rather than replicate measurements, so see
# the per-source summary.md before aggregating across timepoints.
ReactionFluxSchema = pa.DataFrameSchema(
    name="reaction_flux",
    columns={
        "source_id": pa.Column(
            dtype=str,
            nullable=False,
            description=(
                "Provenance key. Must match a directory under "
                "``validation_data/references/<source_id>/`` and (ideally) a "
                "BibTeX entry in ``validation_data/references.bib``."
            ),
        ),
        "reaction_id": pa.Column(
            dtype=str,
            nullable=False,
            description=(
                "Within-vector key — the reaction identifier as the source "
                "names it (model-agnostic, e.g. the paper's abbreviation)."
            ),
        ),
        "value": pa.Column(
            float,
            nullable=False,
            description="Absolute reaction flux, expressed in ``units``.",
        ),
        "units": pa.Column(
            dtype=str,
            nullable=False,
            description="Units of ``value`` (e.g. ``mmol/gDW/h``).",
        ),
        "kind": pa.Column(
            dtype=str,
            nullable=False,
            checks=pa.Check.isin(_KINDS),
            description="One of measured | theoretical_max | model_predicted.",
        ),
        "ecocyc_rxn_id": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description=(
                "Optional mapping of ``reaction_id`` to an EcoCyc reaction id, "
                "for model comparison. Blank where no clean single-reaction "
                "mapping exists."
            ),
        ),
        "reaction_name": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description="Optional human-readable reaction name.",
        ),
        "timepoint": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description=(
                "Optional sampling timepoint/phase label when the source "
                "reports a time course (e.g. ``5h``)."
            ),
        ),
        "value_relative_pct": pa.Column(
            float,
            nullable=True,
            required=False,
            description=(
                "Optional flux normalized to a reference flux (e.g. glucose "
                "uptake = 100); a scale-invariant basis for comparing carbon "
                "routing."
            ),
        ),
        "uncertainty": pa.Column(
            float,
            nullable=True,
            required=False,
            checks=pa.Check.greater_than_or_equal_to(0),
            description="Optional symmetric ± uncertainty, same units as ``value``.",
        ),
        "strain": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description="Strain background (e.g. MG1655), if reported.",
        ),
        "method": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description="How the fluxes were obtained (e.g. 13C-MFA, INST 13C-MFA).",
        ),
        "notes": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description="Free-text caveats: media, growth mode, extraction notes.",
        ),
    },
    strict="filter",  # allow extra columns; validate the required ones
    coerce=True,
    description=(
        "Curated per-reaction flux reference (vector) for one (condition, "
        "observable) slot; one row per (source, reaction, timepoint). The "
        "vector analogue of ScalarClaimSchema for ¹³C-MFA flux maps."
    ),
)
