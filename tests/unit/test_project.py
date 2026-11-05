"""Projection of a pairwise alignment onto reference coordinates."""

from __future__ import annotations

import pytest

from virotrace.align.project import Insertion, project_to_reference
from virotrace.errors import ContractError

REF = "ACGTACGTACGTACGTACGT"  # 20 bases


def test_identical_rows_project_unchanged():
    result = project_to_reference(REF, REF, start=1, end=20)
    assert result.projected == REF and result.insertions == ()


def test_span_selects_reference_positions_inclusive():
    result = project_to_reference(REF, REF, start=5, end=10)
    assert result.projected == REF[4:10]


def test_deletion_and_edges_are_gaps_and_case_is_normalised():
    seq = "--cgtacg---ACGTACGTA"
    result = project_to_reference(REF, seq, start=1, end=20)
    assert result.projected == "--CGTACG---ACGTACGTA"


def test_insertion_columns_are_removed_but_recorded():
    # sequence carries 2 extra bases after reference position 4
    ref_row = "ACGT--ACGTACGTACGTACGT"
    seq_row = "ACGTGGACGTACGTACGTACGT"
    result = project_to_reference(ref_row, seq_row, start=1, end=20)
    assert result.projected == REF
    assert result.insertions == (Insertion(after_ref_pos=4, length=2),)


def test_flanking_bases_beyond_the_span_are_not_insertions():
    ref_row = "--ACGTACGTACGTACGTACGT--"
    seq_row = "GGACGTACGTACGTACGTACGTCC"
    result = project_to_reference(ref_row, seq_row, start=1, end=20)
    assert result.insertions == ()


def test_insertion_that_is_a_gap_in_the_sequence_too_is_ignored():
    ref_row = "ACGT--ACGTACGTACGTACGT"
    seq_row = "ACGT--ACGTACGTACGTACGT"
    assert project_to_reference(ref_row, seq_row, start=1, end=20).insertions == ()


def test_reference_too_short_for_the_span_is_a_contract_error():
    with pytest.raises(ContractError, match="only"):
        project_to_reference(REF, REF, start=1, end=30)


def test_rows_of_different_length_are_a_contract_error():
    with pytest.raises(ContractError, match="differ in length"):
        project_to_reference(REF, REF[:-1], start=1, end=10)
