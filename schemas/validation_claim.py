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
