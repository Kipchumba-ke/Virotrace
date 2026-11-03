"""Reference coordinate integrity and SDRM list contract.

The docs claim HXB2 gene coordinates; these tests make the claim falsifiable.
Once config/reference/hxb2.gb is in place, test_gene_ranges_match_annotated_record
replaces the arithmetic check with an assertion against the annotated record.
"""

from __future__ import annotations

import itertools

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


@pytest.mark.skip(reason="enable once config/reference/hxb2.gb is pinned (Phase 1)")
def test_gene_ranges_match_annotated_record():
    raise NotImplementedError


# ------------------------------------------------------------------------- SDRM
def test_sdrm_list_loads_and_is_non_empty(sdrm_list):
    assert len(sdrm_list) > 0
    assert sdrm_list.version == "placeholder"


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
