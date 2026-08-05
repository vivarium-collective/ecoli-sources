"""
Checks for genotype variant generation and bundle composition
(``processing/genotypes.py``).

The properties exercised here are the ones whose failure would be SILENT: a
variant manifest whose unoverridden paths no longer resolve, an identity that
misses an edited input, provenance that validation quietly discards, and two
generators both claiming one canonical key.

Follows the repo's script idiom (``scripts/validate_all.py``,
``scripts/test_validation_loader.py``): plain asserts, prints ``OK``/``FAIL``
per check, exits non-zero if any failed. No pytest dependency (functions are
``test_*`` so a future pytest run would also discover them). Runnable locally:

    uv run python scripts/test_genotypes.py

Writes only into a temp dir — a generator that wrote into the committed data
package would itself be a bug, and one check asserts it does not.
"""

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ecoli_sources import BUNDLE_PATH  # noqa: E402
from processing.genotypes import (  # noqa: E402
    GeneratorResult,
    compose_variant_bundle,
    genotype_id,
    knockdown,
    read_bundle,
)
from schemas import ReferenceBundleSchema  # noqa: E402

KEY = "rnaseq_experimental_tpms"
N_CONTRACT_KEYS = 135


class Fixture:
    """Shared read-only context: the base bundle plus three real gene ids."""

    def __init__(self, tmp: Path):
        self.tmp = tmp
        self.base = read_bundle(BUNDLE_PATH)
        source = self.base.loc[self.base["canonical_key"] == KEY].iloc[0]["source_path"]
        self.source_path = Path(source)
        self.genes = list(pd.read_csv(source, sep="\t")["gene_id"].head(3))

    def dir(self, name: str) -> Path:
        """A fresh subdirectory, so checks cannot contaminate one another."""
        p = self.tmp / name
        p.mkdir(parents=True, exist_ok=True)
        return p


def _expect(exc_type, match: str, fn, *args, **kwargs) -> None:
    """Assert ``fn`` raises ``exc_type`` whose message contains ``match``."""
    try:
        fn(*args, **kwargs)
    except exc_type as e:
        assert match in str(e), f"expected {match!r} in error, got: {e}"
        return
    raise AssertionError(f"expected {exc_type.__name__} containing {match!r}, none raised")


# ---------------------------------------------------------------------------
# read_bundle / path resolution
# ---------------------------------------------------------------------------


def test_read_bundle_resolves_every_path(fx: Fixture) -> None:
    assert len(fx.base) == N_CONTRACT_KEYS, f"expected {N_CONTRACT_KEYS} rows"
    unresolved = [p for p in fx.base["source_path"] if not Path(p).is_file()]
    assert not unresolved, f"{len(unresolved)} unresolvable, e.g. {unresolved[:3]}"


# ---------------------------------------------------------------------------
# The knockdown generator
# ---------------------------------------------------------------------------


def test_knockdown_scales_only_the_named_genes(fx: Fixture) -> None:
    before = pd.read_csv(fx.source_path, sep="\t").set_index("gene_id")["tpm_mean"]

    result = knockdown(fx.genes, 0.1, fx.dir("kd_scale"))
    after = pd.read_csv(result.rows[KEY], sep="\t").set_index("gene_id")["tpm_mean"]

    for gene in fx.genes:
        assert abs(after[gene] - before[gene] * 0.1) < 1e-9, f"{gene} not scaled"
    untouched = [g for g in before.index if g not in set(fx.genes)]
    assert after[untouched].equals(before[untouched]), "non-target genes changed"


def test_knockdown_records_its_provenance(fx: Fixture) -> None:
    result = knockdown(fx.genes, 0.25, fx.dir("kd_prov"))

    assert result.provenance["generator"] == "knockdown"
    assert result.provenance["operator"] == "scale_gene_set"
    assert result.provenance["params"]["factor"] == 0.25
    assert result.provenance["params"]["gene_ids"] == fx.genes
    assert KEY in result.provenance["reads"]


def test_knockdown_rejects_genes_the_table_does_not_measure(fx: Fixture) -> None:
    # Silently doing nothing is the failure mode: a knockdown of an unmeasured
    # gene yields a variant identical to its base, which reads as a null result.
    _expect(KeyError, "absent from", knockdown, ["NOT_A_REAL_GENE"], 0.1, fx.dir("kd_bad"))


def test_knockdown_leaves_the_shared_dataset_manifest_untouched(fx: Fixture) -> None:
    # Genotype variants are evictable derivatives; the rnaseq dataset manifest is
    # a committed registry. A 1000-genotype campaign must not mutate it 1000 times.
    manifest = Path(BUNDLE_PATH).parent / "rnaseq_experimental" / "manifest.tsv"
    before = manifest.read_bytes() if manifest.is_file() else None

    knockdown(fx.genes, 0.5, fx.dir("kd_manifest"))

    if before is not None:
        assert manifest.read_bytes() == before, "dataset manifest was mutated"


# ---------------------------------------------------------------------------
# Composition
# ---------------------------------------------------------------------------


def test_composed_bundle_is_complete_and_valid(fx: Fixture) -> None:
    result = knockdown(fx.genes, 0.1, fx.dir("c1_files"))

    manifest_path = compose_variant_bundle([result], fx.dir("c1_geno"))

    df = pd.read_csv(manifest_path, sep="\t")
    assert len(df) == N_CONTRACT_KEYS
    ReferenceBundleSchema.validate(df)  # every contract key still present


def test_every_path_in_a_composed_bundle_resolves(fx: Fixture) -> None:
    # ★ The silent failure this guards. Manifest paths resolve relative to the
    # MANIFEST'S OWN directory, so copying one elsewhere without rewriting them
    # yields a manifest that validates cleanly and whose every unoverridden key
    # points at nothing.
    result = knockdown(fx.genes, 0.1, fx.dir("c2_files"))
    manifest_path = compose_variant_bundle([result], fx.dir("c2_geno"))

    root = manifest_path.parent
    df = pd.read_csv(manifest_path, sep="\t")
    broken = [
        (r["canonical_key"], r["source_path"])
        for _, r in df.iterrows()
        if not (root / r["source_path"]).resolve().is_file()
    ]
    assert not broken, f"{len(broken)} unresolvable paths, e.g. {broken[:3]}"


def test_overridden_key_points_at_the_variant(fx: Fixture) -> None:
    result = knockdown(fx.genes, 0.1, fx.dir("c3_files"))
    manifest_path = compose_variant_bundle([result], fx.dir("c3_geno"))

    df = pd.read_csv(manifest_path, sep="\t")
    row = df.loc[df["canonical_key"] == KEY].iloc[0]
    resolved = (manifest_path.parent / row["source_path"]).resolve()

    assert resolved == result.rows[KEY].resolve()


def test_overridden_key_keeps_its_schema_name(fx: Fixture) -> None:
    # A variant must stay validatable: the override inherits the key's schema.
    expected = fx.base.loc[fx.base["canonical_key"] == KEY].iloc[0]["schema_name"]
    result = knockdown(fx.genes, 0.1, fx.dir("c4_files"))

    manifest_path = compose_variant_bundle([result], fx.dir("c4_geno"))

    df = pd.read_csv(manifest_path, sep="\t")
    got = df.loc[df["canonical_key"] == KEY].iloc[0]["schema_name"]
    assert got == expected, f"schema_name {got!r} != {expected!r}"


def test_two_generators_claiming_one_key_raises(fx: Fixture) -> None:
    # A knockout and an insertion both touch `genes`. Last-writer-wins would be a
    # correctness bug dressed as a dict update, so the composer refuses.
    a = knockdown(fx.genes[:1], 0.1, fx.dir("c5_a"))
    b = knockdown(fx.genes[1:2], 0.5, fx.dir("c5_b"))

    _expect(
        ValueError, "written by two generator results",
        compose_variant_bundle, [a, b], fx.dir("c5_geno"),
    )


def test_unknown_canonical_key_raises(fx: Fixture) -> None:
    real = knockdown(fx.genes, 0.1, fx.dir("c6_files"))
    bogus = GeneratorResult(rows={"not_a_contract_key": real.rows[KEY]})

    _expect(
        ValueError, "not in the base bundle",
        compose_variant_bundle, [bogus], fx.dir("c6_geno"),
    )


def test_generator_result_rejects_missing_files(fx: Fixture) -> None:
    _expect(
        FileNotFoundError, "do not exist",
        GeneratorResult, {KEY: fx.tmp / "never_written.tsv"},
    )


# ---------------------------------------------------------------------------
# Identity
# ---------------------------------------------------------------------------


def test_identical_genotypes_share_an_id(fx: Fixture) -> None:
    # Reuse depends on this: the same perturbation must not rebuild.
    a = compose_variant_bundle([knockdown(fx.genes, 0.1, fx.dir("i1_a"))], fx.dir("i1_ga"))
    b = compose_variant_bundle([knockdown(fx.genes, 0.1, fx.dir("i1_b"))], fx.dir("i1_gb"))

    assert genotype_id(a) == genotype_id(b)


def test_different_factors_give_different_ids(fx: Fixture) -> None:
    a = compose_variant_bundle([knockdown(fx.genes, 0.1, fx.dir("i2_a"))], fx.dir("i2_ga"))
    b = compose_variant_bundle([knockdown(fx.genes, 0.5, fx.dir("i2_b"))], fx.dir("i2_gb"))

    assert genotype_id(a) != genotype_id(b)


def test_a_variant_differs_from_its_base(fx: Fixture) -> None:
    variant = compose_variant_bundle(
        [knockdown(fx.genes, 0.1, fx.dir("i3_files"))], fx.dir("i3_geno")
    )

    assert genotype_id(variant) != genotype_id(BUNDLE_PATH)


def test_editing_a_referenced_file_changes_the_id(fx: Fixture) -> None:
    # ★ Why identity is a CONTENT hash. Editing a variant file in place leaves
    # the manifest byte-identical, so a manifest-bytes hash would keep pointing a
    # genotype at a stale ParCa build with nothing to indicate it.
    result = knockdown(fx.genes, 0.1, fx.dir("i4_files"))
    manifest_path = compose_variant_bundle([result], fx.dir("i4_geno"))
    before_manifest = manifest_path.read_bytes()
    before_id = genotype_id(manifest_path)

    edited = pd.read_csv(result.rows[KEY], sep="\t")
    edited.loc[0, "tpm_mean"] = float(edited.loc[0, "tpm_mean"]) + 1.0
    edited.to_csv(result.rows[KEY], sep="\t", index=False)

    assert manifest_path.read_bytes() == before_manifest, "manifest changed unexpectedly"
    assert genotype_id(manifest_path) != before_id, "identity missed an edited input"


def test_id_is_independent_of_location(fx: Fixture) -> None:
    # Content-addressed means re-materializing elsewhere reproduces the id, which
    # is what makes eviction safe rather than lossy.
    a = compose_variant_bundle(
        [knockdown(fx.genes, 0.1, fx.dir("i5_a"))], fx.dir("i5/deep/nested/ga")
    )
    b = compose_variant_bundle([knockdown(fx.genes, 0.1, fx.dir("i5_b"))], fx.dir("i5_gb"))

    assert genotype_id(a) == genotype_id(b)


# ---------------------------------------------------------------------------
# Provenance + materialization discipline
# ---------------------------------------------------------------------------


def test_sidecar_records_lineage(fx: Fixture) -> None:
    manifest_path = compose_variant_bundle(
        [knockdown(fx.genes, 0.1, fx.dir("s1_files"))], fx.dir("s1_geno"), name="kd-demo"
    )

    sidecar = json.loads((manifest_path.parent / "genotype.json").read_text())

    assert sidecar["name"] == "kd-demo"
    assert sidecar["genotype_id"] == genotype_id(manifest_path)
    assert sidecar["base_genotype_id"] == genotype_id(BUNDLE_PATH)
    assert sidecar["overridden_keys"] == [KEY]
    assert sidecar["generators"][0]["params"]["factor"] == 0.1


def test_provenance_would_be_lost_in_the_manifest(fx: Fixture) -> None:
    # ★ Why provenance lives in a sidecar. ReferenceBundleSchema is
    # strict="filter", so an extra column is DROPPED by validation rather than
    # rejected — the in-manifest form appears to work and loses the data.
    df = fx.base.copy()
    df["genotype_provenance"] = "kd:factor=0.1"

    validated = ReferenceBundleSchema.validate(df)

    assert "genotype_provenance" not in validated.columns


def test_variants_are_not_written_into_the_package(fx: Fixture) -> None:
    # Materialization is cache-like: nothing a generator writes belongs inside
    # the committed data package.
    package_root = Path(BUNDLE_PATH).parent.resolve()

    result = knockdown(fx.genes, 0.1, fx.dir("s3_files"))
    manifest_path = compose_variant_bundle([result], fx.dir("s3_geno"))

    for written in (result.rows[KEY], manifest_path):
        assert package_root not in Path(written).resolve().parents, f"{written} is inside the package"


CHECKS = [
    test_read_bundle_resolves_every_path,
    test_knockdown_scales_only_the_named_genes,
    test_knockdown_records_its_provenance,
    test_knockdown_rejects_genes_the_table_does_not_measure,
    test_knockdown_leaves_the_shared_dataset_manifest_untouched,
    test_composed_bundle_is_complete_and_valid,
    test_every_path_in_a_composed_bundle_resolves,
    test_overridden_key_points_at_the_variant,
    test_overridden_key_keeps_its_schema_name,
    test_two_generators_claiming_one_key_raises,
    test_unknown_canonical_key_raises,
    test_generator_result_rejects_missing_files,
    test_identical_genotypes_share_an_id,
    test_different_factors_give_different_ids,
    test_a_variant_differs_from_its_base,
    test_editing_a_referenced_file_changes_the_id,
    test_id_is_independent_of_location,
    test_sidecar_records_lineage,
    test_provenance_would_be_lost_in_the_manifest,
    test_variants_are_not_written_into_the_package,
]


def main() -> int:
    failures = []
    with tempfile.TemporaryDirectory() as td:
        fx = Fixture(Path(td))
        for check in CHECKS:
            try:
                check(fx)
                print(f"OK {check.__name__}")
            except Exception as e:  # noqa: BLE001 — report and continue
                print(f"FAIL {check.__name__}", file=sys.stderr)
                print(f"  {type(e).__name__}: {e}", file=sys.stderr)
                failures.append(check.__name__)
    if failures:
        print(f"\n{len(failures)} genotype check(s) failed:", file=sys.stderr)
        for f in failures:
            print(f"  - {f}", file=sys.stderr)
        return 1
    print(f"\nAll {len(CHECKS)} genotype checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
