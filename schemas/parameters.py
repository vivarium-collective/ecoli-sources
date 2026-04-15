"""
Pandera schemas for whole-cell parameter tables keyed by doubling time.

These tables bind the physiological state of the cell (mass composition,
elongation rates, ppGpp concentration, RNAP count) to the growth rate at
which each value was measured. They are consumed by parca to parameterize
the model for a given condition via interpolation in doubling time.
"""

import pandera.pandas as pa


GrowthRateDependentParametersSchema = pa.DataFrameSchema(
    name="growth_rate_dependent_parameters",
    columns={
        "doublingTime (units.min)": pa.Column(
            float, nullable=False, unique=True,
            checks=pa.Check.greater_than(0),
            description="Cell doubling time in minutes. Primary key.",
        ),
        "fractionActiveRnap": pa.Column(
            float, nullable=False,
            checks=[pa.Check.in_range(0, 1)],
        ),
        "stableRnaPerTotalRnaSynthesized": pa.Column(
            float, nullable=False,
            checks=[pa.Check.in_range(0, 1)],
        ),
        "fractionActiveRnapSynthesizingStableRna": pa.Column(
            float, nullable=False,
            checks=[pa.Check.in_range(0, 1)],
        ),
        "ratioRProteinToTotalProtein": pa.Column(
            float, nullable=False,
            checks=[pa.Check.in_range(0, 1)],
        ),
        "distanceBetweenRibosomesOnMRna (units.nt)": pa.Column(
            float, nullable=False,
            checks=[pa.Check.greater_than(0)],
        ),
        "ribosomeElongationRate (units.aa/units.s)": pa.Column(
            float, nullable=False,
            checks=[pa.Check.greater_than(0)],
        ),
        "rnaPolymeraseElongationRate (units.nt/units.s)": pa.Column(
            float, nullable=False,
            checks=[pa.Check.greater_than(0)],
        ),
        "fractionActiveRibosome": pa.Column(
            float, nullable=False,
            checks=[pa.Check.in_range(0, 1)],
        ),
        "ppGpp_conc (units.pmol/units.ug)": pa.Column(
            float, nullable=False,
            checks=[pa.Check.greater_than_or_equal_to(0)],
        ),
        "RNAP_per_cell": pa.Column(
            float, nullable=False,
            checks=[pa.Check.greater_than(0)],
        ),
    },
    strict="filter",
    coerce=True,
    description=(
        "Per-doubling-time physiological parameters fit to measurement. "
        "Parca interpolates within this table to parameterize the current "
        "condition."
    ),
)


DryMassCompositionSchema = pa.DataFrameSchema(
    name="dry_mass_composition",
    columns={
        "doublingTime (units.min)": pa.Column(
            float, nullable=False, unique=True,
            checks=pa.Check.greater_than(0),
        ),
        "proteinMassFraction": pa.Column(float, nullable=False, checks=pa.Check.in_range(0, 1)),
        "rnaMassFraction": pa.Column(float, nullable=False, checks=pa.Check.in_range(0, 1)),
        "dnaMassFraction": pa.Column(float, nullable=False, checks=pa.Check.in_range(0, 1)),
        "lipidMassFraction": pa.Column(float, nullable=False, checks=pa.Check.in_range(0, 1)),
        "lpsMassFraction": pa.Column(float, nullable=False, checks=pa.Check.in_range(0, 1)),
        "mureinMassFraction": pa.Column(float, nullable=False, checks=pa.Check.in_range(0, 1)),
        "glycogenMassFraction": pa.Column(float, nullable=False, checks=pa.Check.in_range(0, 1)),
        "solublePoolMassFraction": pa.Column(float, nullable=False, checks=pa.Check.in_range(0, 1)),
        "inorganicIonMassFraction": pa.Column(float, nullable=False, checks=pa.Check.in_range(0, 1)),
        "averageDryMass (units.fg)": pa.Column(
            float, nullable=False,
            checks=pa.Check.greater_than(0),
        ),
    },
    strict="filter",
    coerce=True,
    description=(
        "Dry-mass composition fractions (sum ≈ 1 per row) and average dry "
        "cell mass, as a function of doubling time."
    ),
)
