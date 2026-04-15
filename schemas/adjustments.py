"""
Pandera schemas for the ``flat/adjustments/`` family.

These tables encode manual parca-time overrides that are applied on top of
model-derived values. Some adjustments (e.g. 10× RNA expression for specific
genes) are load-bearing; understanding them is central to sensitivity analysis.

Shapes represented here:

* Generic ``name / value / units / _source / _comments`` — used by
  ``rna_expression_adjustments``, ``rna_deg_rates_adjustments``,
  ``translation_efficiencies_adjustments``, ``protein_deg_rates_adjustments``.
* Amino-acid pathway overrides.
* Balanced translation efficiencies (operon-grouped).
* Relative metabolite concentration changes (media-specific).
"""

import pandera.pandas as pa


AdjustmentValueSchema = pa.DataFrameSchema(
    name="adjustment_value",
    columns={
        "name": pa.Column(
            str,
            nullable=False,
            description=(
                "Target identifier. Shape depends on adjusted quantity: "
                "RNA id (e.g. 'EG11493_RNA'), monomer id with compartment "
                "(e.g. 'ADCLY-MONOMER[c]'), or complex id."
            ),
        ),
        "value": pa.Column(
            float,
            nullable=False,
            description="Multiplicative factor (or absolute value, depending on quantity).",
        ),
        "units": pa.Column(
            str,
            nullable=True,
            required=False,
            description="Optional unit string; blank = dimensionless.",
        ),
        "_source": pa.Column(
            str,
            nullable=True,
            required=False,
            description="Origin of the adjustment (script, paper, etc.).",
        ),
        "_comments": pa.Column(
            str,
            nullable=True,
            required=False,
            description="Free-form rationale; often cites gene symbol and purpose.",
        ),
    },
    strict="filter",
    coerce=True,
    description=(
        "Generic name/value adjustment row. Used by "
        "rna_expression_adjustments, rna_deg_rates_adjustments, "
        "translation_efficiencies_adjustments, protein_deg_rates_adjustments."
    ),
)


AminoAcidPathwayAdjustmentSchema = pa.DataFrameSchema(
    name="amino_acid_pathway_adjustment",
    columns={
        "Amino acid": pa.Column(
            str, nullable=False,
            description="Amino acid id with compartment (e.g. 'CYS[c]').",
        ),
        "Parameter": pa.Column(
            str, nullable=False,
            description="Name of the pathway parameter being adjusted (e.g. 'kcat_fwd').",
        ),
        "Factor": pa.Column(
            float, nullable=False,
            description="Multiplicative factor applied to the parameter.",
        ),
        "_comments": pa.Column(
            str, nullable=True, required=False,
            description="Rationale for the adjustment.",
        ),
    },
    strict="filter",
    coerce=True,
    description="Manual overrides on amino-acid biosynthesis pathway parameters.",
)


BalancedTranslationEfficiencyGroupSchema = pa.DataFrameSchema(
    name="balanced_translation_efficiency_group",
    columns={
        "proteins": pa.Column(
            str, nullable=False,
            description=(
                "JSON-encoded list of monomer ids (with compartment) that "
                "share the same operon and are forced to the same translation "
                "efficiency at parca time."
            ),
        ),
        "_comments": pa.Column(
            str, nullable=True, required=False,
        ),
    },
    strict="filter",
    coerce=True,
    description=(
        "Operon-level balance directive: proteins in each row get equal "
        "translation efficiency."
    ),
)


RelativeMetaboliteConcentrationChangeSchema = pa.DataFrameSchema(
    name="relative_metabolite_concentration_change",
    columns={
        "media": pa.Column(
            str, nullable=False,
            description="Media condition label (e.g. 'anaerobic').",
        ),
        "metabolite": pa.Column(
            str, nullable=False,
            description="Metabolite id with compartment.",
        ),
        "fold_change": pa.Column(
            float, nullable=False,
            description="Fold change relative to basal.",
        ),
        "_source": pa.Column(str, nullable=True, required=False),
        "_comments": pa.Column(str, nullable=True, required=False),
    },
    strict="filter",
    coerce=True,
    description=(
        "Manually curated metabolite concentration fold changes applied on "
        "top of relative_metabolite_concentrations.tsv for specific media."
    ),
)
