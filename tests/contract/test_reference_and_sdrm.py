"""Reference coordinate integrity and SDRM list contract.

The docs claim HXB2 gene coordinates; these tests make the claim falsifiable.
Once config/reference/hxb2.gb is in place, test_gene_ranges_match_annotated_record
replaces the arithmetic check with an assertion against the annotated record.
"""

from __future__ import annotations

import itertools
from pathlib import Path

import pytest

from virotrace.drm.mutations import GENE_PROTEIN_LENGTH, parse_mutation
from virotrace.errors import ContractError
from virotrace.settings import load_settings

EXPECTED_PROTEIN_LENGTH = {"pr": 99, "rt": 560, "in_": 288}


def test_gene_spans_are_codon_sized_and_match_protein_lengths(config_path):
    genes = load_settings(config_path).reference.genes
    for key, protein_length in EXPECTED_PROTEIN_LENGTH.items():
        start, end = getattr(genes, key)
        span = end - start + 1
        assert span % 3 == 0, f"{key} span {span} is not a multiple of 3"
        assert span // 3 == protein_length, (
            f"{key} encodes {span // 3} aa, expected {protein_length}"
        )


def test_subgenes_lie_within_pol(config_path):
    genes = load_settings(config_path).reference.genes
    pol_start, pol_end = genes.pol
    for key in ("pr", "rt", "in_"):
        start, end = getattr(genes, key)
        assert pol_start <= start < end <= pol_end, f"{key} falls outside pol"


def test_subgenes_do_not_overlap(config_path):
    genes = load_settings(config_path).reference.genes
    spans = sorted(getattr(genes, key) for key in ("pr", "rt", "in_"))
    for (_, end_a), (start_b, _) in itertools.pairwise(spans):
        assert end_a < start_b, "protease/RT/integrase spans overlap"


def test_protein_length_constants_agree_with_config(config_path):
    genes = load_settings(config_path).reference.genes
    mapping = {"PR": "pr", "RT": "rt", "IN": "in_"}
    for gene, key in mapping.items():
        start, end = getattr(genes, key)
        assert GENE_PROTEIN_LENGTH[gene] == (end - start + 1) // 3


REPO_ROOT_REF = Path(__file__).resolve().parents[2] / "config" / "reference"
HXB2_RECORD = REPO_ROOT_REF / "hxb2.gb"


@pytest.mark.skipif(not HXB2_RECORD.exists(), reason="config/reference/hxb2.gb is not present")
def test_gene_ranges_match_the_pinned_record(config_path):
    """Falsifies the documented coordinates against the real HXB2 sequence."""
    from virotrace.drm.reference import gene_proteins, load_reference_sequence

    genes = load_settings(config_path).reference.genes
    sequence = load_reference_sequence(HXB2_RECORD)
    assert len(sequence) == 9719  # K03455

    proteins = gene_proteins(
        sequence, {"PR": tuple(genes.pr), "RT": tuple(genes.rt), "IN": tuple(genes.in_)}
    )
    pr, rt, integrase = proteins["PR"], proteins["RT"], proteins["IN"]
    assert (len(pr), len(rt), len(integrase)) == (99, 560, 288)
    # Conserved landmarks that only line up if numbering is gene-local and the frame is
    # right: the protease active-site triad (D25-T26-G27), the RT YMDD motif (Y183-M184-
    # D185-D186) and the integrase catalytic DDE triad (D64, D116, E152). HXB2 itself
    # differs from the subtype B consensus at some positions (its PR starts PQVTLW, the
    # consensus PQITLW), so N-terminal residues are not used as the check.
    assert pr[24:27] == "DTG"
    assert rt[182:186] == "YMDD"
    assert (integrase[63], integrase[115], integrase[151]) == ("D", "D", "E")
    for gene, protein in proteins.items():
        assert "*" not in protein, f"{gene} contains an in-frame stop: coordinates are off"

    pol = gene_proteins(sequence, {"POL": tuple(genes.pol)})["POL"]
    assert len(pol) == 1004 and pol.endswith("*") and "*" not in pol[:-1]


# ------------------------------------------------------------------------- SDRM
def test_sdrm_list_loads_and_is_non_empty(sdrm_list):
    assert len(sdrm_list) > 0
    assert sdrm_list.version == "2009"


def test_who_2009_list_has_the_published_93_mutations(sdrm_list):
    """The 2009 list is 34 NRTI + 19 NNRTI + 40 PI mutations (Bennett et al. 2009).

    Substitutions are counted as (position, residue) pairs; the one insertion on the list
    (T69 Ins) is counted separately because it is not a substitution.
    """
    protease = sum(len(v) for (gene, _), v in sdrm_list.entries.items() if gene == "PR")
    rt = sum(len(v) for (gene, _), v in sdrm_list.entries.items() if gene == "RT")
    indels = sum(len(v) for v in sdrm_list.indels.values())
    assert protease == 40
    assert rt + indels == 34 + 19
    assert not any(gene == "IN" for gene, _ in sdrm_list.entries)  # INSTI is not in the 2009 list


def test_stanford_list_is_a_superset_of_the_2009_list(sdrm_list):
    from virotrace.drm.sdrm import load_sdrm_list

    extended = load_sdrm_list(
        REPO_ROOT_REF / "sdrm_stanford_2026-09-25.tsv", version="stanford-2026-09-25"
    )
    for key, residues in sdrm_list.entries.items():
        assert residues <= extended.entries.get(key, frozenset()), key
    assert any(gene == "IN" for gene, _ in extended.entries)


def test_sdrm_positions_are_within_their_gene(sdrm_list):
    for gene, position in sdrm_list.entries:
        assert 1 <= position <= GENE_PROTEIN_LENGTH[gene]


def test_known_sdrm_is_tagged(sdrm_list):
    assert sdrm_list.is_sdrm(parse_mutation("RT:M184V")) is True
    assert sdrm_list.is_sdrm(parse_mutation("RT:K103N")) is True


def test_mixture_is_tagged_if_any_residue_is_an_sdrm(sdrm_list):
    """Mixture policy 'any': resistant if either call is a DRM."""
    mutation = parse_mutation("RT:M184MV")
    assert sdrm_list.is_sdrm(mutation) is True
    assert sdrm_list.matched_residues(mutation) == frozenset({"V"})


def test_non_sdrm_position_is_not_tagged(sdrm_list):
    assert sdrm_list.is_sdrm(parse_mutation("RT:A98G")) is False


def test_missing_list_is_a_contract_error():
    from virotrace.drm.sdrm import load_sdrm_list

    with pytest.raises(ContractError, match="not found"):
        load_sdrm_list("config/reference/does_not_exist.tsv", version="x")
