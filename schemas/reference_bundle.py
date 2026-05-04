"""
Pandera schema for the reference bundle manifest.

The bundle manifest is a flat TSV mapping canonical keys (model-side slot
names) to source paths within the data/ package directory. ParCa-time
consumers in vEcoli resolve each canonical key through this bundle to find
the file (or directory) for that data role. The default reference bundle
ships at ``data/reference_bundle.tsv``; alternative bundles (variants)
follow the same shape.

The schema also enforces the **canonical key contract**: any bundle must
contain every key listed in ``REQUIRED_CANONICAL_KEYS`` below. This
catches variant bundles that drop or rename keys vEcoli's loader expects,
at validation time rather than at consumer call sites deep in ParCa.
"""

import pandera.pandas as pa


# The canonical-key contract — vEcoli expects to be able to resolve each of
# these via ``SourceBundle.get(...)``. Variant bundles are free to add
# extras (and to point any key at a different source file), but they must
# contain every key in this list. Keep alphabetically sorted.
#
# Adding a new key: extend this list AND add the corresponding row to
# ``data/reference_bundle.tsv`` AND make sure vEcoli's loader knows how to
# request it.
REQUIRED_CANONICAL_KEYS: list[str] = [
    "adjustments__amino_acid_pathways",
    "adjustments__balanced_translation_efficiencies",
    "adjustments__protein_deg_rates_adjustments",
    "adjustments__relative_metabolite_concentrations_changes",
    "adjustments__rna_deg_rates_adjustments",
    "adjustments__rna_expression_adjustments",
    "adjustments__translation_efficiencies_adjustments",
    "amino_acid_export_kms",
    "amino_acid_export_kms_removed",
    "amino_acid_pathways",
    "amino_acid_uptake_rates",
    "amino_acid_uptake_rates_removed",
    "base_codes__amino_acids",
    "base_codes__dntp",
    "base_codes__nmp",
    "base_codes__ntp",
    "biomass",
    "cell_wall__murein_strand_length_distribution",
    "compartments",
    "complexation_reactions",
    "complexation_reactions_added",
    "complexation_reactions_modified",
    "complexation_reactions_removed",
    "condition__condition_defs",
    "condition__environment_molecules",
    "condition__media__5X_supplement_EZ",
    "condition__media__MIX0-47",
    "condition__media__MIX0-51",
    "condition__media__MIX0-55",
    "condition__media__MIX0-57",
    "condition__media__MIX0-58",
    "condition__media__MIX0-844",
    "condition__media_recipes",
    "condition__tf_condition",
    "condition__timelines_def",
    "disabled_kinetic_reactions",
    "dna_sites",
    "dna_supercoiling",
    "dry_mass_composition",
    "endoRNases",
    "equilibrium_reaction_rates",
    "equilibrium_reactions",
    "equilibrium_reactions_added",
    "equilibrium_reactions_removed",
    "fold_changes",
    "fold_changes_nca",
    "fold_changes_removed",
    "footprint_sizes",
    "gene_fragments",
    "genes",
    "growth_rate_dependent_parameters",
    "linked_metabolites",
    "mass_fractions__LPS_fractions",
    "mass_fractions__glycogen_fractions",
    "mass_fractions__ion_fractions",
    "mass_fractions__lipid_fractions",
    "mass_fractions__murein_fractions",
    "mass_fractions__soluble_fractions",
    "mass_parameters",
    "metabolic_reactions",
    "metabolic_reactions_added",
    "metabolic_reactions_modified",
    "metabolic_reactions_removed",
    "metabolism_kinetics",
    "metabolite_concentrations",
    "metabolite_concentrations_removed",
    "metabolites",
    "metabolites_added",
    "modified_proteins",
    "molecular_weight_keys",
    "new_gene_data__gfp__gene_sequences",
    "new_gene_data__gfp__genes",
    "new_gene_data__gfp__insertion_location",
    "new_gene_data__gfp__protein_half_lives_measured",
    "new_gene_data__gfp__proteins",
    "new_gene_data__gfp__rna_half_lives",
    "new_gene_data__gfp__rnas",
    "new_gene_data__new_gene_baseline_expression_parameters",
    "new_gene_data__template__gene_sequences",
    "new_gene_data__template__genes",
    "new_gene_data__template__insertion_location",
    "new_gene_data__template__protein_half_lives_measured",
    "new_gene_data__template__proteins",
    "new_gene_data__template__rna_half_lives",
    "new_gene_data__template__rnas",
    "parameters",
    "ppgpp_fc",
    "ppgpp_regulation",
    "ppgpp_regulation_added",
    "ppgpp_regulation_removed",
    "protein_half_lives_measured",
    "protein_half_lives_n_end_rule",
    "protein_half_lives_pulsed_silac",
    "proteins",
    "relative_metabolite_concentrations",
    "rna_half_lives",
    "rna_maturation_enzymes",
    "rna_seq_data__doubling_times",
    "rna_seq_data__rnaseq_rsem_tpm_mean",
    "rna_seq_data__rnaseq_rsem_tpm_std",
    "rna_seq_data__rnaseq_seal_rpkm_mean",
    "rna_seq_data__rnaseq_seal_rpkm_std",
    "rnas",
    "rnaseq_basal_tpms",
    "rnaseq_experimental_tpms",
    "rrna_options__remove_rrff__genes_removed",
    "rrna_options__remove_rrff__mass_parameters_modified",
    "rrna_options__remove_rrff__rnas_removed",
    "rrna_options__remove_rrff__transcription_units_modified",
    "rrna_options__remove_rrna_operons__transcription_units_added",
    "rrna_options__remove_rrna_operons__transcription_units_removed",
    "secretions",
    "sequence",
    "sequence_motifs",
    "tf_one_component_bound",
    "transcription_factors",
    "transcription_units",
    "transcription_units_added",
    "transcription_units_modified",
    "transcription_units_removed",
    "transcriptional_attenuation",
    "transcriptional_attenuation_removed",
    "translation_efficiency",
    "trna_charging_reactions",
    "trna_charging_reactions_added",
    "trna_charging_reactions_removed",
    "trna_data__trna_data",
    "trna_data__trna_growth_rates",
    "trna_data__trna_ratio_to_16SrRNA_0p4",
    "trna_data__trna_ratio_to_16SrRNA_0p7",
    "trna_data__trna_ratio_to_16SrRNA_1p07",
    "trna_data__trna_ratio_to_16SrRNA_1p6",
    "trna_data__trna_ratio_to_16SrRNA_2p5",
    "two_component_system_templates",
    "two_component_systems",
]


def _check_required_canonical_keys(df) -> bool:
    """DataFrame-level check: every REQUIRED_CANONICAL_KEYS entry is present.

    Pandera invokes this with the bundle DataFrame; returns True if every
    required key appears in the ``canonical_key`` column. Pandera's
    failure-cases output names the missing keys when it returns False.
    """
    present = set(df["canonical_key"])
    missing = [k for k in REQUIRED_CANONICAL_KEYS if k not in present]
    if missing:
        # Stash the missing list on the function so downstream tooling
        # (validate_bundle.py) can surface it; pandera's own error path
        # will already report this check failed.
        _check_required_canonical_keys.last_missing = missing  # type: ignore[attr-defined]
        return False
    return True


ReferenceBundleSchema = pa.DataFrameSchema(
    name="reference_bundle",
    columns={
        "canonical_key": pa.Column(
            dtype=str,
            unique=True,
            nullable=False,
            description=(
                "Slot name for an addressable data role in the model "
                "(snake_case, primary key). Stable across file renames."
            ),
        ),
        "source_path": pa.Column(
            dtype=str,
            nullable=False,
            description=(
                "Path to source file or directory, relative to the bundle "
                "root (i.e., the ``data/`` directory of the ecoli-sources "
                "package). Files and directories are both allowed; the "
                "consumer (vEcoli) determines whether to load directly or "
                "traverse based on the canonical key's expected shape."
            ),
        ),
        "description": pa.Column(
            dtype=str,
            nullable=False,
            description="One-liner explaining the slot's role.",
        ),
        "schema_name": pa.Column(
            dtype=str,
            nullable=True,
            required=False,
            description=(
                "Optional: Pandera schema name from ecoli-sources/schemas/ "
                "for content validation (e.g. ``RnaseqTpmTableSchema``). "
                "Empty for canonical keys whose source has no associated "
                "schema yet."
            ),
        ),
    },
    checks=pa.Check(
        _check_required_canonical_keys,
        name="all_required_canonical_keys_present",
        error=(
            "Bundle is missing one or more required canonical keys. See "
            "schemas/reference_bundle.py REQUIRED_CANONICAL_KEYS for the "
            "full contract."
        ),
    ),
    strict="filter",
    coerce=True,
    description=(
        "Bundle manifest mapping canonical keys to source paths. One row "
        "per canonical key. Must contain every entry in "
        "REQUIRED_CANONICAL_KEYS (the canonical-key contract); extras "
        "permitted. The default reference bundle "
        "(``data/reference_bundle.tsv``) is shipped with the package; "
        "alternative bundles (variants) follow the same shape."
    ),
)
