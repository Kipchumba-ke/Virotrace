"""Mutation representation, gene-local numbering, and the no-call guarantee."""

from __future__ import annotations

import pytest
from hypothesis import given
from hypothesis import strategies as st

from virotrace.drm.mutations import (
    GENE_PROTEIN_LENGTH,
    Mutation,
    MutationCalls,
    collapse_mixture,
    parse_mutation,
)
from virotrace.errors import ContractError


def test_label_round_trip():
    assert parse_mutation("RT:M184V").label == "RT:M184V"


def test_mixture_label_is_sorted_and_stable():
    mutation = parse_mutation("RT:M184VI")
    assert mutation.is_mixture
    assert mutation.label == "RT:M184IV"
    assert parse_mutation(mutation.label).label == mutation.label


@given(
    gene=st.sampled_from(sorted(GENE_PROTEIN_LENGTH)),
    ref=st.sampled_from("ACDEFGHIKLMNPQRSTVWY"),
    alt=st.sampled_from("ACDEFGHIKLMNPQRSTVWY"),
)
def test_round_trip_is_a_property(gene, ref, alt):
    position = GENE_PROTEIN_LENGTH[gene] // 2 or 1
    label = f"{gene}:{ref}{position}{alt}"
    assert parse_mutation(label).label == label


def test_hxb2_nucleotide_coordinate_leak_is_rejected():
    """RT 2550 is an HXB2 nucleotide position, not a gene-local codon.

    This confusion is a classic HIV-bioinformatics defect, so the type refuses it
    rather than producing a plausible-looking wrong mutation.
    """
    with pytest.raises(ContractError, match=r"outside 1\.\.560"):
        Mutation(gene="RT", position=2550, ref_aa="M", alt_aa=frozenset("V"))


def test_position_zero_is_rejected():
    with pytest.raises(ContractError):
        Mutation(gene="PR", position=0, ref_aa="M", alt_aa=frozenset("I"))


def test_unparseable_label_raises():
    with pytest.raises(ContractError):
        parse_mutation("POL:M184V")


def test_drug_classes_follow_the_gene():
    assert parse_mutation("RT:M184V").drug_classes == ("NRTI", "NNRTI")
    assert parse_mutation("PR:M46I").drug_classes == ("PI",)
    assert parse_mutation("IN:Q148H").drug_classes == ("INSTI",)


# ------------------------------------------------------------------ no-call guarantee
def test_uncovered_position_is_no_call_not_wild_type():
    """The highest-value test in the repo.

    A caller that reports uncovered positions as wild type produces a reassuringly
    low resistance prevalence and never raises an error.
    """
    calls = MutationCalls(strain="KE_001", covered={"RT": set(range(1, 100))})
    assert calls.status_at("RT", 41) == "called"
    assert calls.status_at("RT", 184) == "no_call"
    assert calls.status_at("IN", 148) == "no_call"


def test_coverage_fraction_and_class_eligibility():
    calls = MutationCalls(
        strain="KE_002",
        covered={"RT": set(range(1, 561)), "PR": set(range(1, 51))},
    )
    assert calls.coverage_fraction("RT") == 1.0
    assert calls.coverage_fraction("PR") == pytest.approx(50 / 99)
    assert calls.coverage_fraction("IN") == 0.0

    assert calls.covers_class("NRTI", 0.80) is True
    assert calls.covers_class("PI", 0.80) is False  # 50/99 is below threshold
    assert calls.covers_class("INSTI", 0.80) is False


def test_unknown_gene_and_class_are_contract_errors():
    calls = MutationCalls(strain="x")
    with pytest.raises(ContractError):
        calls.status_at("POL", 1)
    with pytest.raises(ContractError):
        calls.covers_class("NOPE", 0.8)


# ------------------------------------------------------------------ mixture policy
def test_mixture_policies():
    assert collapse_mixture({"I", "V"}, "any") == frozenset({"I", "V"})
    assert collapse_mixture({"I", "V"}, "strict") == frozenset()
    assert collapse_mixture({"V"}, "strict") == frozenset({"V"})
    # 'majority' is not inferable from sequence alone, so it behaves as strict
    # rather than silently guessing a dominant variant.
    assert collapse_mixture({"I", "V"}, "majority") == frozenset()


def test_mixture_with_no_standard_residue_is_rejected():
    with pytest.raises(ContractError):
        collapse_mixture({"X", "*"}, "any")
