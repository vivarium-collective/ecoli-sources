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
* ``validation_claim`` — the validation-data subsystem: a separate manifest
  (no required-keys contract) + scalar reference-claim tables that model
  *outputs* are graded against. See ``ecoli_sources/validation_data/README.md``.
* ``vector_observation`` — cultivation-keyed tidy tables whose observable is a
  VECTOR (a per-gene proteome, a per-gene transcriptome): the vector analogue
  of ``ScalarObservationSchema``. Measurement-only; resolved by
  ``ecoli_sources.validation.load_vector_observations``.
* ``rnaseq`` — per-condition TPM tables + samples manifest (Chris).
* ``adjustments`` — parca-time manual overrides (flat/adjustments/*).
* ``parameters`` — growth-rate-dependent physiological parameters.
* ``half_lives`` — RNA and protein half-lives.
* ``translation`` — per-gene translation efficiency.
* ``regulation`` — TF table, fold changes, ppGpp regulation.
"""

from .reference_bundle import ReferenceBundleSchema
from .validation_claim import (
    MacromoleculeCompositionSchema,
    MetaboliteConcentrationSchema,
    ProteinAbundanceSchema,
    ReactionFluxSchema,
    ScalarClaimSchema,
    ValidationBundleSchema,
)
from .cultivation import (
    CultivationRegistrySchema,
    ScalarObservationSchema,
)
from .vector_observation import VectorObservationSchema
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
    # validation-data subsystem
    "ValidationBundleSchema",
    "ScalarClaimSchema",
    "ReactionFluxSchema",
    "ProteinAbundanceSchema",
    "MetaboliteConcentrationSchema",
    "MacromoleculeCompositionSchema",
    # cultivation-centric validation layer
    "CultivationRegistrySchema",
    "ScalarObservationSchema",
    "VectorObservationSchema",
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
