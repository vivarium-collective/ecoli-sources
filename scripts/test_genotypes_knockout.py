"""
Tests for the chromosome-deletion knockout generator
(``processing.genotypes.knockout``).

Two classes of property are exercised, and the split matters:

* **Coordinate arithmetic** -- a pure case matrix over every relative position
  a feature can hold against a deletion. These are the defects that were live
  in the code this transform was ported from.
* **Plumbing** -- the generator reads five JSON-typed TSVs and a FASTA, edits
  them field-wise, and writes six coupled canonical keys. This is new code and
  is where new bugs live.

COORDINATE CONVENTION under test: left_end_pos / right_end_pos are 1-based and
INCLUSIVE, so a feature spanning [L, R] occupies sequence[L - 1 : R]. The
convention is the reason these tests exist -- feature LENGTHS stay correct
under an off-by-one splice while the SEQUENCE silently shifts, so any test that
only counts bases passes on a genome that is wrong.

Run: ``uv run pytest scripts/test_genotypes_knockout.py``
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ecoli_sources import BUNDLE_PATH  # noqa: E402
from processing.genotypes import (  # noqa: E402
    KNOCKOUT_KEYS,
    _classify_against_deletion,
    _FlatFile,
    _read_fasta,
    _Row,
    _update_coordinates,
    compose_variant_bundle,
    genotype_id,
    knockout,
    read_bundle,
)

N_CONTRACT_KEYS = 135

# The gene the tranche-A knockout study uses, and the one v2ecoli's own
# deletion suite targets. It sits inside the lac operon, so transcription units
# straddle it -- which is what makes the round-trip's second branch fire.
LACY = "EG10526"

# The real SMS double knockout: disjoint, 744 kb apart.
TRPR = "EG11029"
TNAA = "EG11005"


# ---------------------------------------------------------------------------
# Fixtures. The base bundle is read once per module: it is ~5.5 MB of flat
# files and every test wants the same immutable view of it.
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def base_paths() -> dict[str, str]:
    base = read_bundle(BUNDLE_PATH)
    return dict(zip(base["canonical_key"], base["source_path"]))


@pytest.fixture(scope="module")
def wild_type(base_paths) -> dict[str, object]:
    """Base coordinates and genome, keyed the way the tests want to ask."""
    genes = _FlatFile(base_paths["genes"])
    _, sequence, width = _read_fasta(Path(base_paths["sequence"]))
    return {
        "genes": genes,
        "coords": {
            r.get("id"): (r.get("left_end_pos"), r.get("right_end_pos"))
            for r in genes.rows
            if r.get("left_end_pos") is not None
        },
        "sequence": sequence,
        "width": width,
    }


@pytest.fixture(scope="module")
def lacy_ko(tmp_path_factory) -> object:
    """One real single-gene knockout, shared by the genome-scale tests."""
    return knockout([LACY], tmp_path_factory.mktemp("lacy_ko"))


def _feature_seq(sequence: str, left: int, right: int) -> str:
    """1-based inclusive [left, right] -> the bases it occupies."""
    return sequence[left - 1 : right]


def _coord_bearing(files: dict[str, _FlatFile]):
    for kind, flat in files.items():
        for row in flat.rows:
            if row.get("left_end_pos") is not None:
                yield kind, row


def _load_variant(result) -> dict[str, _FlatFile]:
    return {
        k: _FlatFile(result.rows[k]) for k in KNOCKOUT_KEYS if k != "sequence"
    }


# ---------------------------------------------------------------------------
# The case matrix: every relative position of a feature against a deletion.
# ---------------------------------------------------------------------------

DEL_LEFT, DEL_RIGHT = 500, 600
DEL_LEN = DEL_RIGHT - DEL_LEFT + 1

# (label, input coords, expected output coords, or None when it loses its place)
CASES = [
    # Entirely upstream -- untouched, including the feature ending on the base
    # immediately before the cut.
    ("before", (100, 200), (100, 200)),
    ("before_abutting", (400, DEL_LEFT - 1), (400, DEL_LEFT - 1)),
    # Entirely downstream -- shifts left by the deleted length, including the
    # feature starting on the base immediately after.
    ("after", (700, 800), (700 - DEL_LEN, 800 - DEL_LEN)),
    ("after_abutting", (DEL_RIGHT + 1, 700), (DEL_LEFT, 700 - DEL_LEN)),
    # Wholly inside the deletion -- loses its position with it.
    ("contained", (520, 540), None),
    ("contained_exact", (DEL_LEFT, DEL_RIGHT), None),
    ("contained_at_left_edge", (DEL_LEFT, 540), None),
    ("contained_at_right_edge", (560, DEL_RIGHT), None),
    # Starts before, ends inside -- truncated at the cut.
    ("overlaps_left", (450, 550), (450, DEL_LEFT - 1)),
    ("overlaps_left_minimal", (450, DEL_LEFT), (450, DEL_LEFT - 1)),
    ("overlaps_left_to_right_edge", (450, DEL_RIGHT), (450, DEL_LEFT - 1)),
    # Starts inside, ends after -- surviving 3' portion begins at the cut.
    ("overlaps_right", (550, 700), (DEL_LEFT, 700 - DEL_LEN)),
    ("overlaps_right_minimal", (DEL_RIGHT, 700), (DEL_LEFT, 700 - DEL_LEN)),
    ("overlaps_right_from_left_edge", (DEL_LEFT, 700), (DEL_LEFT, 700 - DEL_LEN)),
    # Starts before, ends after -- keeps both flanks, loses del_len.
    ("spans", (400, 700), (400, 700 - DEL_LEN)),
    ("spans_minimal", (DEL_LEFT - 1, DEL_RIGHT + 1), (DEL_LEFT - 1, DEL_LEFT)),
]

_COLUMNS = ["id", "common_name", "left_end_pos", "right_end_pos", "genes"]


def _row(id_: str, left, right, common_name=None, genes=None) -> _Row:
    """Build one in-memory row without touching a file."""
    values = [id_, common_name, left, right, genes]
    raw = "\t".join("" if v is None and i > 1 else json.dumps(v)
                    for i, v in enumerate(values))
    return _Row(raw, _COLUMNS, Path("synthetic.tsv"), 1)


@pytest.mark.parametrize("kind", ["gene", "tu", "dna_site"])
@pytest.mark.parametrize("label,given,expected", CASES, ids=[c[0] for c in CASES])
def test_coordinate_case_matrix(kind, label, given, expected):
    """Every relative position, for every feature kind."""
    rows = [_row(label, *given, common_name=label)]

    if expected is None:
        with pytest.warns(UserWarning, match="contained within the deletion"):
            _update_coordinates(rows, kind, "EG_TEST", DEL_LEFT, DEL_RIGHT)
    else:
        _update_coordinates(rows, kind, "EG_TEST", DEL_LEFT, DEL_RIGHT)

    row = rows[0]
    if expected is None:
        if kind == "gene":
            # D16 extended: a gene never disappears, it is tombstoned.
            assert not row.removed, "a contained GENE must be tombstoned, not removed"
            assert (row.get("left_end_pos"), row.get("right_end_pos")) == (None, None)
        else:
            assert row.removed, f"a contained {kind} should be removed"
    else:
        assert not row.removed, f"{label} should have been retained"
        assert (row.get("left_end_pos"), row.get("right_end_pos")) == expected


def test_classification_is_total():
    """
    Every (left <= right) position classifies into exactly one case -- no
    fallthrough. The sweep must also reach every branch, or it proves nothing.
    """
    seen = set()
    for left in range(495, 606):
        for right in range(left, 610):
            case = _classify_against_deletion(left, right, DEL_LEFT, DEL_RIGHT)
            assert case in {
                "before", "after", "contained",
                "spans", "overlaps_left", "overlaps_right",
            }
            seen.add(case)

    assert seen == {
        "before", "after", "contained",
        "spans", "overlaps_left", "overlaps_right",
    }


def test_removal_does_not_perturb_iteration():
    """
    A contained row must not cause the row after it to be skipped -- the
    classic mutate-while-iterating defect. Flagging rather than deleting makes
    this structural, and this pins it.
    """
    rows = [
        _row("contained_first", 520, 540),
        _row("downstream", 700, 800),
        _row("contained_second", 530, 550),
        _row("downstream_2", 900, 1000),
    ]

    with pytest.warns(UserWarning):
        _update_coordinates(rows, "dna_site", "EG_TEST", DEL_LEFT, DEL_RIGHT)

    survivors = {r.get("id"): (r.get("left_end_pos"), r.get("right_end_pos"))
                 for r in rows if not r.removed}
    assert survivors == {
        "downstream": (700 - DEL_LEN, 800 - DEL_LEN),
        "downstream_2": (900 - DEL_LEN, 1000 - DEL_LEN),
    }


@pytest.mark.parametrize("left,right", [(None, 600), (500, None), (None, None)])
def test_rows_without_usable_coordinates_are_skipped(left, right):
    """
    A half-populated row must be skipped, not compared -- comparing would raise
    on None. Both spellings of absent occur in the shipped data.
    """
    rows = [_row("partial", left, right)]

    _update_coordinates(rows, "gene", "EG_TEST", DEL_LEFT, DEL_RIGHT)

    assert (rows[0].get("left_end_pos"), rows[0].get("right_end_pos")) == (left, right)


def test_only_content_losing_features_are_annotated():
    """
    Truncated features are marked; features that merely shifted are not --
    otherwise every feature downstream of a cut would be tagged.
    """
    shifted = _row("shifted", 700, 800, common_name="shifted")
    truncated = _row("truncated", 450, 550, common_name="truncated")

    _update_coordinates([shifted, truncated], "tu", "EG_X", DEL_LEFT, DEL_RIGHT)

    assert shifted.get("common_name") == "shifted"
    assert truncated.get("common_name") == "truncated_removed_EG_X"


def test_annotation_is_idempotent_and_null_safe():
    already = _row("t", 450, 550, common_name="t_removed_EG_X")
    null_named = _row("u", 450, 550, common_name=None)

    _update_coordinates([already, null_named], "tu", "EG_X", DEL_LEFT, DEL_RIGHT)

    assert already.get("common_name") == "t_removed_EG_X"
    assert null_named.get("common_name") == "_removed_EG_X"


# ---------------------------------------------------------------------------
# Untouched rows are byte-identical -- the property that makes line-oriented
# editing worth doing.
# ---------------------------------------------------------------------------


def test_rows_upstream_of_the_cut_are_byte_identical(lacy_ko, wild_type, base_paths):
    """
    ★ The plumbing invariant. Re-serializing a JSON-typed TSV through a parser
    would rewrite `null` as empty, strip quoting, reformat list columns and
    drop the `#` provenance headers -- for every row, including the ~4,700 that
    the deletion never touches. That output still validates, still hashes
    stably, and only fails at ParCa build time.

    Editing field-wise makes untouched rows byte-identical BY CONSTRUCTION;
    this pins it against the real files rather than trusting the construction.
    """
    del_left, _ = wild_type["coords"][LACY]
    checked = 0

    for key in KNOCKOUT_KEYS:
        if key == "sequence":
            continue
        before = Path(base_paths[key]).read_text().split("\n")
        after = Path(lacy_ko.rows[key]).read_text().split("\n")
        assert len(before) == len(after), f"{key}: line count changed"

        flat = _FlatFile(base_paths[key])
        # Comment/header preamble must survive verbatim.
        assert before[: len(flat.preamble)] == after[: len(flat.preamble)], (
            f"{key}: the comment header was not preserved"
        )

        offset = len(flat.preamble)
        for i, row in enumerate(flat.rows):
            right = row.get("right_end_pos")
            if right is not None and right < del_left and LACY not in (row.get("genes") or []):
                assert before[offset + i] == after[offset + i], (
                    f"{key} row {row.get('id')!r} lies upstream of the deletion "
                    "but its bytes changed"
                )
                checked += 1

    assert checked > 500, f"only {checked} upstream rows checked; too weak to mean anything"


def test_the_fasta_keeps_its_header_and_line_width(lacy_ko, base_paths):
    """
    Biopython's SeqIO.write re-wraps at 60 columns; this file is wrapped at 70.
    Re-wrapping would rewrite all 66,311 lines of a 4.5 MB file for no reason.
    """
    base_header, _, base_width = _read_fasta(Path(base_paths["sequence"]))
    ko_header, ko_seq, ko_width = _read_fasta(Path(lacy_ko.rows["sequence"]))

    assert ko_header == base_header
    assert ko_width == base_width == 70
    # Every line but the last is full width.
    lines = Path(lacy_ko.rows["sequence"]).read_text().rstrip("\n").split("\n")[1:]
    assert all(len(line) == 70 for line in lines[:-1])
    assert len(ko_seq) == sum(len(line) for line in lines)


# ---------------------------------------------------------------------------
# Genome-scale, against the real knowledge base.
# ---------------------------------------------------------------------------


def test_the_genome_shortens_by_exactly_the_gene_length(lacy_ko, wild_type):
    left, right = wild_type["coords"][LACY]
    _, ko_seq, _ = _read_fasta(Path(lacy_ko.rows["sequence"]))

    assert len(ko_seq) == len(wild_type["sequence"]) - (right - left + 1)


def test_the_deleted_gene_is_tombstoned_not_removed(lacy_ko, wild_type):
    genes = _FlatFile(lacy_ko.rows["genes"])
    row = next((r for r in genes.rows if r.get("id") == LACY), None)

    assert row is not None, "the knocked-out gene row must survive as a tombstone"
    assert (row.get("left_end_pos"), row.get("right_end_pos")) == (None, None)
    assert len(genes.rows) == len(wild_type["genes"].rows), "no gene row may vanish"


def test_genome_round_trip_invariant(lacy_ko, wild_type, capsys):
    """
    ★ THE genome-scale invariant, in two branches:

      * a feature whose length is unchanged kept all of its content, so its
        sequence must be byte-identical across the deletion;
      * a feature that shrank straddled the deletion, so its sequence must be
        its original with exactly the deleted span excised -- both flanks
        intact and correctly rejoined.

    An off-by-one splice preserves every LENGTH, so only the second branch can
    catch it. A single-branch form of this test is vacuous for the defect it
    exists to detect.
    """
    del_left, _ = wild_type["coords"][LACY]
    wt_seq = wild_type["sequence"]
    _, ko_seq, _ = _read_fasta(Path(lacy_ko.rows["sequence"]))

    base_files = _load_base_files()
    ko_files = _load_variant(lacy_ko)

    before = {
        (kind, row.get("id")): _feature_seq(
            wt_seq, row.get("left_end_pos"), row.get("right_end_pos")
        )
        for kind, row in _coord_bearing(base_files)
    }

    intact, shrunk, corrupted = 0, 0, []
    for kind, row in _coord_bearing(ko_files):
        if row.get("id") == LACY:
            continue
        key = (kind, row.get("id"))
        if key not in before:
            continue

        old = before[key]
        new = _feature_seq(ko_seq, row.get("left_end_pos"), row.get("right_end_pos"))

        if len(new) == len(old):
            intact += 1
            if new != old:
                corrupted.append((kind, row.get("id"), "content changed"))
            continue

        shrunk += 1
        lost = len(old) - len(new)
        # The kept 5' flank is whatever lay upstream of the cut.
        offset = max(del_left - _base_left(base_files, kind, row.get("id")), 0)
        expected = old[:offset] + old[offset + lost :]
        if new != expected:
            corrupted.append((kind, row.get("id"), "excision misaligned"))

    with capsys.disabled():
        print(f"\n  round-trip: {intact} intact, {shrunk} shrunk, {len(corrupted)} corrupted")

    assert not corrupted, f"{len(corrupted)} corrupted features: {corrupted[:10]}"
    assert intact > 10000, f"expected the full feature set, saw {intact}"
    assert shrunk > 0, "no feature straddled the deletion -- the second branch never fired"


def _load_base_files() -> dict[str, _FlatFile]:
    base = read_bundle(BUNDLE_PATH)
    paths = dict(zip(base["canonical_key"], base["source_path"]))
    return {k: _FlatFile(paths[k]) for k in KNOCKOUT_KEYS if k != "sequence"}


def _base_left(base_files, kind, id_):
    row = next(r for r in base_files[kind].rows if r.get("id") == id_)
    return row.get("left_end_pos")


def test_coordinates_stay_within_the_chromosome(lacy_ko):
    _, ko_seq, _ = _read_fasta(Path(lacy_ko.rows["sequence"]))
    genome_length = len(ko_seq)

    out_of_range = [
        (kind, row.get("id"), row.get("left_end_pos"), row.get("right_end_pos"))
        for kind, row in _coord_bearing(_load_variant(lacy_ko))
        if not 1 <= row.get("left_end_pos") <= row.get("right_end_pos") <= genome_length
    ]

    assert not out_of_range, f"features outside [1, {genome_length}]: {out_of_range[:10]}"


def test_an_abutting_downstream_gene_keeps_its_sequence(tmp_path, wild_type):
    """
    ★ The seam case, and the second independent detector of an off-by-one.

    It has to be a gene starting on the base IMMEDIATELY after the cut. A
    length-preserving off-by-one excises [L+1, R+1] instead of [L, R]: every
    feature separated from the target by even one intergenic base keeps all of
    its own bases and is merely translated, so it survives the defect
    untouched. Only a feature abutting the cut loses a base to it.

    36 exactly-abutting gene pairs exist in the shipped data, so this is real
    rather than synthetic. Testing "the next 25 genes downstream" instead
    catches nothing -- an intergenic gap absorbs the error.

    ★ And the pair has to be CHOSEN, not taken. Under the defect the neighbour
    reads one base too early, so a pair whose boundary bases happen to be equal
    is corrupted identically to its own original and detects nothing. That is
    15 of the 36 pairs -- including the very first one, an A/A collision. A
    test that took `next(abutting_pairs)` would look like a seam test and
    assert nothing.
    """
    wt_seq = wild_type["sequence"]
    by_left = sorted((left, right, gene)
                     for gene, (left, right) in wild_type["coords"].items())
    discriminating = [
        (a, b) for a, b in zip(by_left, by_left[1:])
        if b[0] == a[1] + 1 and wt_seq[a[0] - 1] != wt_seq[b[0] - 1]
    ][:5]

    assert discriminating, "no abutting pair with distinguishable boundary bases"

    for (t_left, _, target_id), (n_left, n_right, neighbour_id) in discriminating:
        result = knockout([target_id], tmp_path / f"abutting_{target_id}")
        _, ko_seq, _ = _read_fasta(Path(result.rows["sequence"]))
        row = next(r for r in _FlatFile(result.rows["genes"]).rows
                   if r.get("id") == neighbour_id)

        after = _feature_seq(ko_seq, row.get("left_end_pos"), row.get("right_end_pos"))
        expected = _feature_seq(wt_seq, n_left, n_right)

        assert after == expected, (
            f"{neighbour_id} abuts the deletion of {target_id} and its sequence "
            f"shifted: {after[:6]!r} != {expected[:6]!r} -- classic off-by-one, "
            "which changes a start codon while preserving every length"
        )


# ---------------------------------------------------------------------------
# Multi-gene: order invariance, disjoint deletions, overlap.
# ---------------------------------------------------------------------------


def test_input_order_does_not_change_the_output(tmp_path):
    """
    ★ The decision-3 test. The genotype id is a content hash over the produced
    files, so if caller order changed the bytes, {trpR, tnaA} and {tnaA, trpR}
    would be two genotypes -- and two ParCa builds -- for one strain.
    """
    forward = knockout([TRPR, TNAA], tmp_path / "forward")
    reverse = knockout([TNAA, TRPR], tmp_path / "reverse")

    for key in KNOCKOUT_KEYS:
        assert Path(forward.rows[key]).read_bytes() == Path(reverse.rows[key]).read_bytes(), (
            f"{key} differs between input orderings"
        )

    a = compose_variant_bundle([forward], tmp_path / "geno_a")
    b = compose_variant_bundle([reverse], tmp_path / "geno_b")
    assert genotype_id(a) == genotype_id(b)

    # The sidecar too, or the artifact is only half order-independent.
    sa = json.loads((a.parent / "genotype.json").read_text())
    sb = json.loads((b.parent / "genotype.json").read_text())
    assert sa["generators"] == sb["generators"]
    assert sa["generators"][0]["params"]["gene_ids"] == [TRPR, TNAA], (
        "provenance should record the canonical applied order (descending left_end_pos)"
    )


def test_two_disjoint_deletions_shift_by_the_right_amounts(tmp_path, wild_type):
    """
    A feature between two deletions shifts by the downstream one only; a
    feature after both shifts by their sum; a feature before both is untouched.
    Getting this wrong is exactly what the descending-coordinate order prevents.
    """
    lo_left, lo_right = wild_type["coords"][TNAA]
    hi_left, hi_right = wild_type["coords"][TRPR]
    assert lo_right < hi_left, "fixture assumption: tnaA lies upstream of trpR"
    lo_len = lo_right - lo_left + 1
    hi_len = hi_right - hi_left + 1

    result = knockout([TRPR, TNAA], tmp_path / "double")
    genes = {r.get("id"): (r.get("left_end_pos"), r.get("right_end_pos"))
             for r in _FlatFile(result.rows["genes"]).rows}

    before_both, between, after_both = 0, 0, 0
    for gene, (left, right) in wild_type["coords"].items():
        if gene in (TRPR, TNAA):
            continue
        got = genes[gene]
        if right < lo_left:
            assert got == (left, right), f"{gene} upstream of both, should not move"
            before_both += 1
        elif left > lo_right and right < hi_left:
            assert got == (left - lo_len, right - lo_len), f"{gene} between the cuts"
            between += 1
        elif left > hi_right:
            assert got == (left - lo_len - hi_len, right - lo_len - hi_len), (
                f"{gene} downstream of both"
            )
            after_both += 1

    assert before_both and between and after_both, (
        f"a region was empty: {before_both}/{between}/{after_both}"
    )


def test_overlapping_targets_raise(tmp_path, wild_type):
    """
    No application order deletes two overlapping genes cleanly: the second
    target's coordinates have already been truncated or nulled by the first.
    """
    coords = wild_type["coords"]
    by_left = sorted(coords.items(), key=lambda kv: kv[1][0])
    pair = next(
        (a, b)
        for (a, (al, ar)), (b, (bl, br)) in zip(by_left, by_left[1:])
        if bl <= ar and al <= br
    )

    with pytest.raises(ValueError, match="overlap"):
        knockout(list(pair), tmp_path / "overlap")


def test_repeated_gene_ids_raise(tmp_path):
    with pytest.raises(ValueError, match="repeated gene ids"):
        knockout([LACY, LACY], tmp_path / "dup")


# ---------------------------------------------------------------------------
# Transcription-unit coupling, including the two overlays.
# ---------------------------------------------------------------------------


def test_transcription_units_lose_the_deleted_gene(lacy_ko):
    for key in ("transcription_units", "transcription_units_added",
                "transcription_units_modified"):
        for row in _FlatFile(lacy_ko.rows[key]).rows:
            members = row.get("genes") or []
            if LACY in members:
                # The only legitimate survivor is a solo TU, which is nulled
                # rather than emptied.
                assert len(members) == 1 and not _has_coords(row), (
                    f"{key} row {row.get('id')!r} still carries the deleted gene"
                )


def _has_coords(row) -> bool:
    return row.get("left_end_pos") is not None and row.get("right_end_pos") is not None


def test_the_modified_overlay_shifts_with_the_base(lacy_ko, wild_type, base_paths):
    """
    ★ Why this is six keys and not four. `transcription_units_modified` carries
    real coordinates that ParCa's _modify_data applies OVER the base table. Omit
    it and those rows stay in pre-deletion coordinate space -- and the stale
    values win, silently, for any deletion upstream of them.
    """
    del_left, del_right = wild_type["coords"][LACY]
    del_len = del_right - del_left + 1

    base = {r.get("id"): (r.get("left_end_pos"), r.get("right_end_pos"))
            for r in _FlatFile(base_paths["transcription_units_modified"]).rows}
    after = {r.get("id"): (r.get("left_end_pos"), r.get("right_end_pos"))
             for r in _FlatFile(lacy_ko.rows["transcription_units_modified"]).rows}

    downstream = {k: v for k, v in base.items() if v[0] is not None and v[0] > del_right}
    assert downstream, "fixture assumption: a modified TU lies downstream of lacY"

    for tu_id, (left, right) in downstream.items():
        assert after[tu_id] == (left - del_len, right - del_len), (
            f"modified-overlay TU {tu_id!r} did not shift with the base table"
        )


# ---------------------------------------------------------------------------
# Composition, identity and materialization discipline, over six coupled keys.
# ---------------------------------------------------------------------------


def test_the_composed_bundle_overrides_exactly_six_keys(lacy_ko, tmp_path):
    manifest_path = compose_variant_bundle([lacy_ko], tmp_path / "geno")

    sidecar = json.loads((manifest_path.parent / "genotype.json").read_text())
    assert sorted(sidecar["overridden_keys"]) == sorted(KNOCKOUT_KEYS)
    assert len(sidecar["overridden_keys"]) == 6


def test_every_path_in_a_composed_bundle_resolves(lacy_ko, tmp_path):
    """
    The silent failure: manifest paths resolve relative to the manifest's own
    directory, so an unrelativized one validates cleanly and points at nothing.
    A knockout overrides six of 135 keys, so 129 must still reach base data.
    """
    manifest_path = compose_variant_bundle([lacy_ko], tmp_path / "geno")
    root = manifest_path.parent
    df = pd.read_csv(manifest_path, sep="\t")

    assert len(df) == N_CONTRACT_KEYS
    broken = [
        (r["canonical_key"], r["source_path"])
        for _, r in df.iterrows()
        if not (root / r["source_path"]).resolve().is_file()
    ]
    assert not broken, f"{len(broken)} unresolvable paths, e.g. {broken[:3]}"

    overridden = {Path(p).resolve() for p in lacy_ko.rows.values()}
    still_base = sum(
        1 for _, r in df.iterrows()
        if (root / r["source_path"]).resolve() not in overridden
    )
    assert still_base == N_CONTRACT_KEYS - 6 == 129


def test_editing_any_of_the_six_changes_the_id(lacy_ko, tmp_path):
    """
    Identity is a content hash over referenced files, not the manifest's bytes:
    editing a variant file in place leaves the manifest byte-identical, so a
    bytes-hash would silently reuse a stale ParCa build.
    """
    manifest_path = compose_variant_bundle([lacy_ko], tmp_path / "geno")
    before_manifest = manifest_path.read_bytes()
    before_id = genotype_id(manifest_path)

    target = Path(lacy_ko.rows["dna_sites"])
    target.write_text(target.read_text() + "\n")

    assert manifest_path.read_bytes() == before_manifest
    assert genotype_id(manifest_path) != before_id, "identity missed an edited input"


def test_variants_are_not_written_into_the_package(lacy_ko, tmp_path):
    """Materialization is cache-like; nothing written belongs in the package."""
    package_root = Path(BUNDLE_PATH).parent.resolve()
    manifest_path = compose_variant_bundle([lacy_ko], tmp_path / "geno")

    for written in list(lacy_ko.rows.values()) + [manifest_path]:
        assert package_root not in Path(written).resolve().parents, (
            f"{written} is inside the committed data package"
        )


def test_provenance_records_what_was_deleted(lacy_ko, wild_type):
    prov = lacy_ko.provenance

    assert prov["generator"] == "knockout"
    assert prov["params"]["gene_ids"] == [LACY]
    assert prov["params"]["coordinates"][LACY] == list(wild_type["coords"][LACY])
    assert set(prov["reads"]) == set(KNOCKOUT_KEYS)


# ---------------------------------------------------------------------------
# Refusals.
# ---------------------------------------------------------------------------


def test_deleting_an_unknown_gene_raises(tmp_path):
    with pytest.raises(KeyError, match="no such gene"):
        knockout(["NOT_A_REAL_GENE"], tmp_path / "unknown")


def test_deleting_a_gene_with_no_coordinates_raises(tmp_path, wild_type):
    """
    21 of 4,747 shipped genes permanently carry null coordinates. Knocking one
    out would silently produce a variant identical to its base.
    """
    phantom = next(
        r.get("id") for r in wild_type["genes"].rows
        if r.get("left_end_pos") is None
    )

    with pytest.raises(KeyError, match="no coordinates"):
        knockout([phantom], tmp_path / "phantom")


def test_an_empty_knockout_raises(tmp_path):
    with pytest.raises(ValueError, match="at least one gene"):
        knockout([], tmp_path / "empty")
