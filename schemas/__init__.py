"""
Pandera schemas for vEcoli input data files.

These schemas define the canonical formats for experimental inputs and
parca-time parameter tables. Validation at ingestion ensures consistent
structure across datasets and enables like-for-like substitution
(e.g. reference vs alternate RNA-seq transcriptomes, measured vs
synthesized adjustment tables).

Organization:

* ``reference_bundle`` — bundle manifest (canonical_key → source_path)
  defining the data package's contract with vEcoli.
* ``rnaseq`` — per-condition TPM tables + samples manifest (Chris).
* ``adjustments`` — parca-time manual overrides (flat/adjustments/*).
* ``parameters`` — growth-rate-dependent physiological parameters.
* ``half_lives`` — RNA and protein half-lives.
* ``translation`` — per-gene translation efficiency.
* ``regulation`` — TF table, fold changes, ppGpp regulation.
"""

from .reference_bundle import ReferenceBundleSchema
from .adjustments import (
    AdjustmentValueSchema,
    AminoAcidPathwayAdjustmentSchema,
    BalancedTranslationEfficiencyGroupSchema,
    RelativeMetaboliteConcentrationChangeSchema,
)
from .half_lives import (
    ProteinHalfLivesMeasuredSchema,
    ProteinHalfLivesNEndRuleSchema,
    ProteinHalfLivesPulsedSilacSchema,
    ProteinHalfLivesSchema,
    RnaHalfLivesSchema,
)
from .parameters import (
    DryMassCompositionSchema,
    GrowthRateDependentParametersSchema,
)
from .regulation import (
    FoldChangesSchema,
    PpgppRegulationSchema,
    TranscriptionFactorsSchema,
)
from .rnaseq import (
    RnaseqSamplesManifestSchema,
    RnaseqTpmTableSchema,
)
from .translation import TranslationEfficiencySchema

__all__ = [
    # reference_bundle
    "ReferenceBundleSchema",
    # rnaseq
    "RnaseqTpmTableSchema",
    "RnaseqSamplesManifestSchema",
    # adjustments
    "AdjustmentValueSchema",
    "AminoAcidPathwayAdjustmentSchema",
    "BalancedTranslationEfficiencyGroupSchema",
    "RelativeMetaboliteConcentrationChangeSchema",
    # parameters
    "GrowthRateDependentParametersSchema",
    "DryMassCompositionSchema",
    # half-lives
    "RnaHalfLivesSchema",
    "ProteinHalfLivesSchema",
    "ProteinHalfLivesMeasuredSchema",
    "ProteinHalfLivesPulsedSilacSchema",
    "ProteinHalfLivesNEndRuleSchema",
    # translation
    "TranslationEfficiencySchema",
    # regulation
    "TranscriptionFactorsSchema",
    "FoldChangesSchema",
    "PpgppRegulationSchema",
]
