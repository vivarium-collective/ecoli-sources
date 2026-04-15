"""
Pandera schemas for RNA and protein half-life tables.

These feed first-order degradation rates (k_deg = ln(2)/half_life) used
throughout parca for fitting expression and degradation.
"""

import pandera.pandas as pa


RnaHalfLivesSchema = pa.DataFrameSchema(
    name="rna_half_lives",
    columns={
        "id": pa.Column(
            str, nullable=False, unique=True,
            description="Cistron / gene id (e.g. 'EG10001').",
        ),
        "half_life (units.min)": pa.Column(
            float, nullable=False,
            checks=pa.Check.greater_than_or_equal_to(0),
            description=(
                "RNA half-life in minutes. Zero (or ``-0.0``) is permitted "
                "as a no-data / placeholder marker (e.g. G6492)."
            ),
        ),
    },
    strict="filter",
    coerce=True,
    description="Measured RNA half-lives keyed by gene id.",
)


ProteinHalfLivesMeasuredSchema = pa.DataFrameSchema(
    name="protein_half_lives_measured",
    columns={
        "id": pa.Column(
            str, nullable=False, unique=True,
            description="Monomer id (optionally with compartment tag).",
        ),
        # Note: the "measured" file uses a space; "pulsed_silac" uses an
        # underscore. Kept as-is to match the source.
        "half life (units.min)": pa.Column(
            float, nullable=False,
            checks=pa.Check.greater_than(0),
            description="Protein half-life in minutes.",
        ),
        "_comments": pa.Column(str, nullable=True, required=False),
    },
    strict="filter",
    coerce=True,
    description="Manually measured protein half-lives (small curated set).",
)


ProteinHalfLivesPulsedSilacSchema = pa.DataFrameSchema(
    name="protein_half_lives_pulsed_silac",
    columns={
        "id": pa.Column(
            str, nullable=False, unique=True,
            description="Monomer id.",
        ),
        "half_life (units.min)": pa.Column(
            float, nullable=False,
            checks=pa.Check.greater_than(0),
            description="Protein half-life in minutes (pulsed-SILAC estimate).",
        ),
    },
    strict="filter",
    coerce=True,
    description=(
        "Protein half-lives measured by pulsed SILAC (Li et al.). "
        "Note: the column name here uses an underscore, unlike the "
        "``measured`` table which uses a space. Both are retained as-is."
    ),
)


# Back-compat alias: many call sites just want "protein half-lives".
ProteinHalfLivesSchema = ProteinHalfLivesMeasuredSchema


ProteinHalfLivesNEndRuleSchema = pa.DataFrameSchema(
    name="protein_half_lives_n_end_rule",
    columns={
        "aa_code": pa.Column(
            str, nullable=False, unique=True,
            description="Single-letter amino-acid code for N-terminal residue.",
        ),
        "aa_id": pa.Column(
            str, nullable=False, unique=True,
            description="Full amino-acid id with compartment (e.g. 'L-ALPHA-ALANINE[c]').",
        ),
        "half life (units.min)": pa.Column(
            float, nullable=False,
            checks=pa.Check.greater_than(0),
            description=(
                "Protein half-life (min) assumed for monomers with this "
                "N-terminal residue."
            ),
        ),
    },
    strict="filter",
    coerce=True,
    description=(
        "N-end rule: protein half-life as a function of N-terminal residue "
        "(Tobias et al., 1986/1991)."
    ),
)
