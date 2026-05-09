"""
Pandera schemas for transcriptional regulation tables:
transcription factors, fold-change targets, and ppGpp regulation.
"""

import pandera.pandas as pa


TranscriptionFactorsSchema = pa.DataFrameSchema(
    name="transcription_factors",
    columns={
        "TF": pa.Column(
            str, nullable=False,
            description="TF gene symbol (e.g. 'nhaR').",
        ),
        "geneId": pa.Column(
            str, nullable=True,
            description=(
                "EcoCyc gene id for the TF (e.g. 'EG11078'). Nullable — "
                "some rows carry a TF symbol but no EcoCyc gene mapping."
            ),
        ),
        "oneComponentId": pa.Column(
            str, nullable=True, required=False,
            description="Monomer/complex id when TF is a one-component system.",
        ),
        "twoComponentId": pa.Column(
            str, nullable=True, required=False,
            description="Monomer/complex id when TF is part of a two-component system.",
        ),
        "nonMetaboliteBindingId": pa.Column(
            str, nullable=True, required=False,
            description="Id of the TF form that binds DNA without a small-molecule ligand.",
        ),
        "activeId": pa.Column(
            str, nullable=True, required=False,
            description="Id of the transcriptionally active TF form.",
        ),
        "_notes": pa.Column(str, nullable=True, required=False),
    },
    strict="filter",
    coerce=True,
    description=(
        "Canonical TF table. Rows without explicit active/bound forms are "
        "inferred by parca from the one-component column."
    ),
)


FoldChangesSchema = pa.DataFrameSchema(
    name="fold_changes",
    columns={
        "TF": pa.Column(
            str, nullable=False,
            description="TF gene symbol (regulator).",
        ),
        "Target": pa.Column(
            str, nullable=False,
            description="Target gene symbol (regulated).",
        ),
        "log2 FC mean": pa.Column(
            float, nullable=False,
            description="Mean log2 fold change of target under TF knockout/overexpression.",
        ),
        "log2 FC std": pa.Column(
            float, nullable=True, required=False,
            checks=pa.Check.greater_than_or_equal_to(0),
            description="Std dev of log2 fold change across replicates.",
        ),
        "Regulation_direct": pa.Column(
            float, nullable=True, required=False,
            description=(
                "Direct-regulation confidence code (typically 1=direct, "
                "2=inferred; meaning defined by the source curation). "
                "Stored as float so all-NaN columns (e.g. NCA-derived "
                "fold-change tables, where this annotation is not "
                "carried) validate alongside curated tables that "
                "populate it with 1/2."
            ),
        ),
    },
    strict="filter",
    coerce=True,
    description="TF → target fold-change measurements used for P-solve fitting.",
)


PpgppRegulationSchema = pa.DataFrameSchema(
    name="ppgpp_regulation",
    columns={
        "Gene": pa.Column(
            str, nullable=False,
            description="Regulated gene symbol.",
        ),
        "Curated Gene": pa.Column(
            str, nullable=True, required=False,
            description="Synonym used to align with wcm gene symbols.",
        ),
        "ppGpp": pa.Column(
            # Cells may be blank; read_value() in transcription.py treats blank as 0.
            pa.Category,
            nullable=True, required=False,
            description="ppGpp regulation direction: '+', '-', or blank.",
        ),
        "DksA-ppGpp": pa.Column(
            pa.Category,
            nullable=True, required=False,
            description="DksA-ppGpp regulation direction: '+', '-', or blank.",
        ),
    },
    strict="filter",
    coerce=True,
    description=(
        "Manually curated ppGpp / DksA-ppGpp regulation targets (EcoCyc, "
        "with duplicate/direction normalization)."
    ),
)
