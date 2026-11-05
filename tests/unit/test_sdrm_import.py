"""SDRM transcription helper and the translation it checks against.

The mutation text below is synthetic test data, not the WHO list.
"""

from __future__ import annotations

import pytest

from virotrace.drm.reference import CODON_TABLE, gene_proteins, translate
from virotrace.drm.sdrm import load_sdrm_list
from virotrace.drm.sdrm_import import check_wild_types, parse_sdrm_text, render_tsv
from virotrace.errors import ContractError

SYNTHETIC = """
NRTI: M41L, K65R, L74V/I
NNRTI: K103NS
PI: M46IL; I47V
"""


def test_codon_table_spot_checks():
    assert len(CODON_TABLE) == 64
    assert (CODON_TABLE["ATG"], CODON_TABLE["TGG"], CODON_TABLE["TTT"]) == ("M", "W", "F")
    assert {CODON_TABLE[c] for c in ("TAA", "TAG", "TGA")} == {"*"}
    assert CODON_TABLE["GGG"] == "G" and CODON_TABLE["AAA"] == "K"


def test_translate_and_gene_slicing():
    assert translate("ATGAAATAG") == "MK*"
    with pytest.raises(ContractError, match="multiple of 3"):
        translate("ATGA")
    assert gene_proteins("CCATGAAATAGCC", {"X": (3, 11)}) == {"X": "MK*"}


def test_parse_assigns_gene_by_class_and_merges_residues():
    rows = {(r.drug_class, r.position): r for r in parse_sdrm_text(SYNTHETIC)}
    assert rows[("NRTI", 74)].residues == "IV" and rows[("NRTI", 74)].gene == "RT"
    assert rows[("NNRTI", 103)].gene == "RT" and rows[("NNRTI", 103)].residues == "NS"
    assert rows[("PI", 46)].gene == "PR" and rows[("PI", 46)].residues == "IL"


@pytest.mark.parametrize("bad", ["NRTI: M41", "NRTI: 41L", "NRTI: M41LZ!", "M41L"])
def test_unreadable_tokens_fail_loudly_instead_of_being_dropped(bad):
    with pytest.raises(ContractError):
        parse_sdrm_text(bad)


def test_wild_type_check_reports_mismatch_and_out_of_range():
    rows = parse_sdrm_text("NRTI: M41L, K65R, M999L")
    protein = "M" * 40 + "M" + "A" * 30  # position 41 = M, position 65 = A
    problems = check_wild_types(rows, {"RT": protein})
    assert len(problems) == 2
    assert any("K65R" in p and "HXB2 has A" in p for p in problems)
    assert any("outside" in p for p in problems)


def test_rendered_tsv_round_trips_through_the_real_loader(tmp_path):
    rows = parse_sdrm_text(SYNTHETIC)
    path = tmp_path / "sdrm.tsv"
    path.write_text(render_tsv(rows, source="synthetic", version="test", input_text=SYNTHETIC))
    loaded = load_sdrm_list(path, version="test")
    assert loaded.entries[("RT", 74)] == frozenset("IV")
    assert loaded.entries[("PR", 46)] == frozenset("IL")
    assert "synthetic" in path.read_text() and "sha256" in path.read_text()


def test_insertions_and_deletions_stay_out_of_residues():
    (row,) = [r for r in parse_sdrm_text("NRTI: T69D/Ins/Del") if r.position == 69]
    assert row.residues == "D"  # "Ins" must never be read as the residues I, N, S
    assert row.indels == "Del,Ins"


def test_indels_survive_the_round_trip_and_do_not_affect_substitution_tagging(tmp_path):
    from virotrace.drm.mutations import parse_mutation

    text = "NRTI: K70R/E/Del"
    path = tmp_path / "sdrm.tsv"
    path.write_text(render_tsv(parse_sdrm_text(text), source="s", version="t", input_text=text))
    loaded = load_sdrm_list(path, version="t")
    assert loaded.entries[("RT", 70)] == frozenset("ER")
    assert loaded.indels[("RT", 70)] == frozenset({"Del"})
    assert loaded.is_sdrm(parse_mutation("RT:K70R")) is True
