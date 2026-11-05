"""Codon-aware calling from a reference-projected alignment (synthetic reference)."""

from __future__ import annotations

import random

import pytest

from virotrace.align.project import Insertion
from virotrace.drm.caller import (
    call_projected,
    classify_codon,
    insertions_by_gene,
)
from virotrace.drm.reference import CODON_TABLE, translate
from virotrace.settings import load_settings

POL_START, POL_END = 2085, 5096
SENSE_CODONS = sorted(c for c, aa in CODON_TABLE.items() if aa != "*")


@pytest.fixture(scope="module")
def spans(config_path):
    genes = load_settings(config_path).reference.genes
    return {"PR": tuple(genes.pr), "RT": tuple(genes.rt), "IN": tuple(genes.in_)}


@pytest.fixture(scope="module")
def ref_pol():
    """A 3012-nt synthetic pol with no stop codon except the terminal one."""
    rng = random.Random(1)
    codons = [rng.choice(SENSE_CODONS) for _ in range((POL_END - POL_START + 1) // 3 - 1)]
    return "".join(codons) + "TAA"


@pytest.fixture(scope="module")
def ref_protein(ref_pol, spans):
    return {
        gene: translate(ref_pol[start - POL_START : end - POL_START + 1])
        for gene, (start, end) in spans.items()
    }


def codon_index(spans, gene, position):
    """Index in the projected string of the first base of ``gene`` codon ``position``."""
    return spans[gene][0] - POL_START + 3 * (position - 1)


def edit(projected, index, replacement):
    return projected[:index] + replacement + projected[index + len(replacement) :]


def call(projected, spans, ref_protein, **kwargs):
    return call_projected(
        "S1",
        projected,
        pol_start=POL_START,
        gene_spans=spans,
        reference_protein=ref_protein,
        **kwargs,
    )


def alt_codon_for(residue):
    return next(c for c in SENSE_CODONS if CODON_TABLE[c] == residue)


# ---------------------------------------------------------------- classification
@pytest.mark.parametrize(
    ("codon", "kind", "residues"),
    [
        ("ATG", "aa", {"M"}),
        ("RTG", "aa", {"M", "V"}),  # A/G -> M or V
        ("TAA", "stop", set()),
        ("TAR", "stop", set()),  # TAA / TAG
        ("TGR", "aa", {"W"}),  # TGA (stop) / TGG (W): judged on the amino-acid reading
        ("NNN", "too_ambiguous", set()),
        ("AT-", "invalid", set()),
    ],
)
def test_classify_codon(codon, kind, residues):
    got_kind, got_residues = classify_codon(codon, max_residues=4)
    assert got_kind == kind and got_residues == frozenset(residues)


# --------------------------------------------------------------------- full gene
def test_wild_type_gives_full_coverage_and_no_mutations(ref_pol, spans, ref_protein):
    result = call(ref_pol, spans, ref_protein)
    assert result.calls.mutations == []
    assert result.calls.coverage_fraction("PR") == 1.0
    assert result.calls.coverage_fraction("RT") == 1.0
    assert result.calls.coverage_fraction("IN") == 1.0


def test_substitution_is_called_in_gene_local_numbering(ref_pol, spans, ref_protein):
    ref_aa = ref_protein["RT"][183]
    target = "V" if ref_aa != "V" else "I"
    projected = edit(ref_pol, codon_index(spans, "RT", 184), alt_codon_for(target))
    result = call(projected, spans, ref_protein)
    assert result.calls.labels() == {f"RT:{ref_aa}184{target}"}


def test_mixture_keeps_every_residue_under_the_any_policy(ref_pol, spans, ref_protein):
    # RTG is A/G-TG = Met or Val. Whatever the synthetic reference residue is at this
    # position, the call must carry both residues and be flagged as a mixture.
    index = codon_index(spans, "RT", 103)
    result = call(edit(ref_pol, index, "RTG"), spans, ref_protein, mixture_policy="any")
    reference = ref_protein["RT"][102]
    expected = frozenset({"M", "V"})
    if expected == {reference}:  # impossible: a mixture has two residues
        raise AssertionError("test setup error")
    (mutation,) = result.calls.mutations
    assert mutation.position == 103 and mutation.alt_aa == expected and mutation.is_mixture
    assert 103 in result.calls.covered["RT"]


def test_strict_policy_turns_a_mixture_into_a_no_call(ref_pol, spans, ref_protein):
    index = codon_index(spans, "RT", 50)
    projected = edit(ref_pol, index, "RTG")  # M or V
    result = call(projected, spans, ref_protein, mixture_policy="strict")
    assert 50 not in result.calls.covered["RT"]
    assert result.calls.mutations == []


def test_fully_ambiguous_codon_is_a_no_call_not_wild_type(ref_pol, spans, ref_protein):
    projected = edit(ref_pol, codon_index(spans, "RT", 50), "NNN")
    result = call(projected, spans, ref_protein)
    assert 50 not in result.calls.covered["RT"]
    assert result.genes["RT"].n_too_ambiguous == 1
    assert result.calls.mutations == []


# ------------------------------------------------------------------- coverage
def test_uncovered_positions_are_no_calls_never_wild_type(ref_pol, spans, ref_protein):
    """A sequence covering only the first 100 RT codons says nothing about the rest."""
    rt_start = spans["RT"][0] - POL_START
    keep_from, keep_to = rt_start, rt_start + 300
    projected = "-" * keep_from + ref_pol[keep_from:keep_to] + "-" * (len(ref_pol) - keep_to)
    result = call(projected, spans, ref_protein)

    assert result.calls.covered["PR"] == set()
    assert result.calls.covered["IN"] == set()
    assert result.calls.covered["RT"] == set(range(1, 101))
    assert result.calls.status_at("RT", 101) == "no_call"
    assert result.calls.status_at("RT", 100) == "called"
    assert result.calls.coverage_fraction("RT") == pytest.approx(100 / 560)


def test_codon_split_by_the_edge_of_coverage_is_not_called(ref_pol, spans, ref_protein):
    rt_start = spans["RT"][0] - POL_START
    projected = "-" * (rt_start + 1) + ref_pol[rt_start + 1 :]  # first RT codon loses a base
    result = call(projected, spans, ref_protein)
    assert 1 not in result.calls.covered["RT"] and 2 in result.calls.covered["RT"]


# ------------------------------------------------------------- indels and stops
def test_in_frame_codon_deletion_is_covered_and_not_a_substitution(ref_pol, spans, ref_protein):
    projected = edit(ref_pol, codon_index(spans, "RT", 69), "---")
    result = call(projected, spans, ref_protein)
    assert 69 in result.genes["RT"].deleted
    assert 69 in result.calls.covered["RT"]
    assert result.calls.mutations == []
    assert result.genes["RT"].usable


def test_one_base_deletion_is_a_frameshift_and_disqualifies_the_gene(ref_pol, spans, ref_protein):
    projected = edit(ref_pol, codon_index(spans, "RT", 200), "-")
    result = call(projected, spans, ref_protein)
    assert result.genes["RT"].frameshift
    assert result.calls.covered["RT"] == set()
    assert result.calls.coverage_fraction("PR") == 1.0  # other genes are unaffected


def test_deletion_straddling_a_gene_boundary_is_judged_at_full_length(ref_pol, spans, ref_protein):
    """PR and RT are adjacent; a 3-base deletion across the join is in frame overall."""
    join = spans["PR"][1] - POL_START + 1  # first base of RT
    projected = edit(ref_pol, join - 1, "---")  # last PR base + first two RT bases
    result = call(projected, spans, ref_protein)
    assert not result.genes["PR"].frameshift and not result.genes["RT"].frameshift


def test_off_frame_insertion_is_a_frameshift_in_the_gene_it_falls_in(ref_pol, spans, ref_protein):
    rt_mid = spans["RT"][0] - POL_START + 600 + 1  # pol-relative position inside RT
    insertions = [Insertion(after_ref_pos=rt_mid, length=1)]
    by_gene = insertions_by_gene(insertions, spans, POL_START)
    assert by_gene["RT"] == [1] and by_gene["PR"] == [] and by_gene["IN"] == []
    result = call(ref_pol, spans, ref_protein, insertion_lengths=by_gene)
    assert result.genes["RT"].frameshift and not result.genes["PR"].frameshift


def test_in_frame_insertion_is_not_a_frameshift(ref_pol, spans, ref_protein):
    rt_mid = spans["RT"][0] - POL_START + 601
    by_gene = insertions_by_gene([Insertion(rt_mid, 3)], spans, POL_START)
    assert not call(ref_pol, spans, ref_protein, insertion_lengths=by_gene).genes["RT"].frameshift


def test_internal_stop_disqualifies_the_gene(ref_pol, spans, ref_protein):
    projected = edit(ref_pol, codon_index(spans, "IN", 100), "TAG")
    result = call(projected, spans, ref_protein)
    assert result.genes["IN"].internal_stops == [100]
    assert result.calls.covered["IN"] == set()
    assert result.calls.coverage_fraction("RT") == 1.0


def test_empty_sequence_covers_nothing(spans, ref_protein):
    result = call("-" * (POL_END - POL_START + 1), spans, ref_protein)
    assert all(result.calls.covered[g] == set() for g in spans)
