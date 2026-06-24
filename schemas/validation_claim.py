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


# ---------------------------------------------------------------------------
# Protein-abundance claim table: a VECTOR reference (per-gene proteome)
# ---------------------------------------------------------------------------

# The proteome analogue of ReactionFluxSchema: a per-gene protein-abundance
# vector for one (condition, observable) slot, keyed within the vector by
# ``gene``. One row per (source, gene). A consumer grades the model's per-gene
# protein counts against this vector (typically a log-log R²). Abundances are
# absolute copies per cell; ``cv`` carries the reported coefficient of variation
# across replicates when available.
ProteinAbundanceSchema = pa.DataFrameSchema(
    name="protein_abundance",
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
        "gene": pa.Column(
            dtype=str,
            nullable=False,
            description=(
                "Within-vector key — the gene name as the source reports it "
                "(model-agnostic; mapping to a model monomer is a downstream "
                "concern)."
            ),
        ),
        "value": pa.Column(
            float,
            nullable=False,
            description="Protein abundance, expressed in ``units``.",
        ),
        "units": pa.Column(
            dtype=str,
            nullable=False,
            description="Units of ``value`` (e.g. ``copies/cell``).",
        ),
        "kind": pa.Column(
            dtype=str,
            nullable=False,
            checks=pa.Check.isin(_KINDS),
            description="One of measured | theoretical_max | model_predicted.",
        ),
        "uniprot_id": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description="Optional UniProt accession for the protein.",
        ),
        "cv": pa.Column(
            float,
            nullable=True,
            required=False,
            checks=pa.Check.greater_than_or_equal_to(0),
            description=(
                "Optional coefficient of variation across replicates "
                "(relative; the source's reported measurement spread)."
            ),
        ),
        "condition": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description="Growth condition (e.g. ``M9 glucose``), if relevant.",
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
            description="How abundance was measured (e.g. LC-MS/MS iBAQ).",
        ),
        "notes": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description="Free-text caveats.",
        ),
    },
    strict="filter",  # allow extra columns; validate the required ones
    coerce=True,
    description=(
        "Curated per-gene protein-abundance reference (vector) for one "
        "(condition, observable) slot; one row per (source, gene). The proteome "
        "analogue of ReactionFluxSchema."
    ),
)


# ---------------------------------------------------------------------------
# Metabolite-concentration claim table: a VECTOR reference (per-metabolite pools)
# ---------------------------------------------------------------------------

# The metabolome analogue of ReactionFluxSchema / ProteinAbundanceSchema: a
# per-metabolite absolute-concentration vector, keyed within the vector by
# ``metabolite``. One row per (source, metabolite, condition) — a single source
# (e.g. Bennett 2009) reports the same metabolite across several carbon sources,
# so ``condition`` discriminates the rows. ``ci_low`` / ``ci_high`` carry the
# reported (often asymmetric, e.g. 95% CI) concentration bounds; ``significance``
# carries a source-specific flag for cross-condition differences.
MetaboliteConcentrationSchema = pa.DataFrameSchema(
    name="metabolite_concentration",
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
        "metabolite": pa.Column(
            dtype=str,
            nullable=False,
            description=(
                "Within-vector key — the metabolite name as the source reports "
                "it (model-agnostic; mapping to a model metabolite id is a "
                "downstream concern, carried in ``ecocyc_id`` when resolved)."
            ),
        ),
        "value": pa.Column(
            float,
            nullable=False,
            checks=pa.Check.greater_than_or_equal_to(0),
            description="Absolute intracellular concentration, expressed in ``units``.",
        ),
        "units": pa.Column(
            dtype=str,
            nullable=False,
            description="Units of ``value`` (e.g. ``mol/L``).",
        ),
        "kind": pa.Column(
            dtype=str,
            nullable=False,
            checks=pa.Check.isin(_KINDS),
            description="One of measured | theoretical_max | model_predicted.",
        ),
        "condition": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description=(
                "Carbon source / growth condition (e.g. ``glucose``, "
                "``glycerol``, ``acetate``). Discriminates a source's rows when "
                "it reports the same metabolite across conditions."
            ),
        ),
        "ci_low": pa.Column(
            float,
            nullable=True,
            required=False,
            checks=pa.Check.greater_than_or_equal_to(0),
            description="Lower bound of the reported concentration interval (e.g. 95% CI).",
        ),
        "ci_high": pa.Column(
            float,
            nullable=True,
            required=False,
            checks=pa.Check.greater_than_or_equal_to(0),
            description="Upper bound of the reported concentration interval (e.g. 95% CI).",
        ),
        "significance": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description=(
                "Optional source-specific flag for statistically significant "
                "cross-condition differences."
            ),
        ),
        "ecocyc_id": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description=(
                "Optional mapping of ``metabolite`` to a model/EcoCyc metabolite "
                "id, for model comparison. Blank where not yet resolved."
            ),
        ),
        "strain": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description="Strain background (e.g. NCM3722), if reported.",
        ),
        "method": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description="How concentrations were measured (e.g. LC-MS/MS isotope-ratio).",
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
        "Curated per-metabolite absolute-concentration reference (vector) for "
        "one (condition, observable) slot; one row per (source, metabolite, "
        "condition). The metabolome analogue of ReactionFluxSchema."
    ),
)


# ---------------------------------------------------------------------------
# Macromolecular-composition claim table: a VECTOR reference (dry-mass fractions)
# ---------------------------------------------------------------------------

# A per-component dry-mass-fraction vector (protein / RNA / DNA), keyed within
# the vector by ``component``. Distinct from the ParCa input
# ``DryMassCompositionSchema`` (schemas/parameters.py), which is data fed INTO
# the model; this is an independent literature reference model outputs are graded
# AGAINST. Macromolecular composition is growth-rate-dependent, so one row per
# (source, component, doubling_time): a consumer grades against the row (or
# interpolation) matching the model's growth rate.
MacromoleculeCompositionSchema = pa.DataFrameSchema(
    name="macromolecule_composition",
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
        "component": pa.Column(
            dtype=str,
            nullable=False,
            description=(
                "Within-vector key — the dry-mass component (e.g. ``protein``, "
                "``rna``, ``dna``). RNA is total RNA (rRNA+tRNA+mRNA)."
            ),
        ),
        "value": pa.Column(
            float,
            nullable=False,
            checks=pa.Check.in_range(0, 1),
            description="Fraction of total dry weight (0–1), expressed in ``units``.",
        ),
        "units": pa.Column(
            dtype=str,
            nullable=False,
            description="Units of ``value`` (e.g. ``fraction_dry_weight``).",
        ),
        "kind": pa.Column(
            dtype=str,
            nullable=False,
            checks=pa.Check.isin(_KINDS),
            description="One of measured | theoretical_max | model_predicted.",
        ),
        "doubling_time_min": pa.Column(
            float,
            nullable=True,
            required=False,
            checks=pa.Check.greater_than(0),
            description=(
                "Doubling time (min) — the growth-rate index. Composition is "
                "growth-rate-dependent; grade against the matching doubling time."
            ),
        ),
        "growth_rate_doublings_per_h": pa.Column(
            float,
            nullable=True,
            required=False,
            checks=pa.Check.greater_than(0),
            description="Growth rate in doublings/h (= 60 / doubling_time_min).",
        ),
        "strain": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description="Strain background (e.g. B/r), if reported.",
        ),
        "method": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description="How composition was obtained.",
        ),
        "notes": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description="Free-text caveats: media, strain, growth-rate-curve provenance.",
        ),
    },
    strict="filter",  # allow extra columns; validate the required ones
    coerce=True,
    description=(
        "Curated per-component dry-mass-fraction reference (vector) for one "
        "(condition, observable) slot; one row per (source, component, "
        "doubling_time). An independent literature reference (graded AGAINST), "
        "distinct from the ParCa-input DryMassCompositionSchema."
    ),
)


# ---------------------------------------------------------------------------
# Response-curve claim table: a SWEPT-RESPONSE reference (observable vs a driver)
# ---------------------------------------------------------------------------

# A perturbation/response reference: an observable measured along a swept driver
# (e.g. acetate excretion vs growth rate — Basan 2015's overflow curve). Unlike
# the basal (single-condition) schemas, this carries a whole CURVE per source:
# one row per (source, series, level), keyed within by ``series`` (one curve) and
# ordered along ``growth_rate_per_h`` (the x-axis the source plots against, which
# in the model is an EMERGENT output of the swept knob). ``perturbation``/``level``
# record HOW the point was set (e.g. inducer level, carbon source, mutant). A
# single source contributes several series (Basan: titration / carbon-source /
# mutant), discriminated by ``series_type``; a consumer grades the apples-to-apples
# series for its sweep (e.g. the glucose uptake-titration series) with the
# ``curve_response`` criterion (overflow onset + slope). ``value`` may be negative
# (measurement noise around zero below onset), so it carries no non-negativity check.
ResponseCurveSchema = pa.DataFrameSchema(
    name="response_curve",
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
        "series": pa.Column(
            dtype=str,
            nullable=False,
            description=(
                "Within-vector key — the curve this point belongs to (e.g. "
                "``ptsG_glucose``). One source reports several series; a consumer "
                "grades the series matching its sweep."
            ),
        ),
        "growth_rate_per_h": pa.Column(
            float,
            nullable=False,
            checks=pa.Check.greater_than_or_equal_to(0),
            description=(
                "The curve x-axis: growth rate (1/h) the source plots the "
                "observable against. In the model this is an EMERGENT output at "
                "each swept-knob setting, matching the source's causal direction."
            ),
        ),
        "value": pa.Column(
            float,
            nullable=False,
            description=(
                "The observable (curve y-value) at this point, in ``units``. May "
                "be negative (measurement noise around zero below the onset)."
            ),
        ),
        "units": pa.Column(
            dtype=str,
            nullable=False,
            description="Units of ``value`` (e.g. ``mM/OD600/h``).",
        ),
        "kind": pa.Column(
            dtype=str,
            nullable=False,
            checks=pa.Check.isin(_KINDS),
            description="One of measured | theoretical_max | model_predicted.",
        ),
        "series_type": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description=(
                "Optional class of the series — how growth rate was varied (e.g. "
                "``uptake_titration``, ``carbon_source``, ``uptake_mutant``)."
            ),
        ),
        "condition": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description="Media / strain descriptor for this point, as the source reports it.",
        ),
        "perturbation": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description=(
                "Name of the swept variable that sets this point (e.g. "
                "``inducer_3MBA_uM``, ``carbon_source``, ``strain``)."
            ),
        ),
        "level": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description="Value/label of ``perturbation`` at this point (e.g. ``300``, ``glycerol 0.2%``).",
        ),
        "ci_low": pa.Column(
            float,
            nullable=True,
            required=False,
            description="Lower bound of the reported interval on ``value`` (e.g. 95% CI), if reported.",
        ),
        "ci_high": pa.Column(
            float,
            nullable=True,
            required=False,
            description="Upper bound of the reported interval on ``value`` (e.g. 95% CI), if reported.",
        ),
        "strain": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description="Strain background (e.g. NCM3722), if reported.",
        ),
        "method": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description="How the curve was obtained (e.g. batch culture, transporter titration).",
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
        "Curated swept-response reference (curve) for one (perturbation, "
        "observable) slot; one row per (source, series, level), ordered along "
        "``growth_rate_per_h``. The perturbation analogue of the basal vector "
        "schemas, graded with the ``curve_response`` criterion."
    ),
)
