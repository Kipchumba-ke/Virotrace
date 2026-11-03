"""Integration: FASTA I/O -> QC -> subsample, on fixtures, no network."""

from __future__ import annotations

import pytest

from virotrace.curate.qc import apply_qc
from virotrace.curate.subsample import stratified_subsample
from virotrace.errors import ContractError
from virotrace.io.fasta import assert_equal_length, read_fasta, write_fasta
from virotrace.provenance import build_provenance, write_provenance

GENE_NT = 99


def _records() -> list[tuple[str, str]]:
    good = "ACG" * 33
    return [
        ("KE_001", good),
        ("KE_002", good[:-3] + "ATG"),
        ("KE_003", good[:-6] + "ATGATG"),
        ("KE_bad", "ACG"),  # fails QC-1
    ]


def test_fasta_round_trip(tmp_path):
    path = tmp_path / "seqs.fasta"
    written = write_fasta(path, _records())
    assert written == 4
    assert list(read_fasta(path)) == _records()


def test_header_policy_rejects_embedded_metadata(tmp_path):
    path = tmp_path / "bad.fasta"
    path.write_text(">KE_001|Kenya|2015\nACGT\n")
    with pytest.raises(ContractError, match="header policy"):
        list(read_fasta(path))


def test_duplicate_strain_on_write_is_rejected(tmp_path):
    with pytest.raises(ContractError, match="duplicate strain"):
        write_fasta(tmp_path / "dup.fasta", [("KE_001", "ACGT"), ("KE_001", "ACGT")])


def test_ragged_alignment_is_rejected():
    with pytest.raises(ContractError, match="ragged"):
        assert_equal_length([("a", "ACGT"), ("b", "ACG")])


def test_qc_then_subsample_preserves_referential_integrity(tmp_path):
    """Every strain leaving subsampling must have survived QC."""
    path = tmp_path / "seqs.fasta"
    write_fasta(path, _records())

    qc_result = apply_qc(
        read_fasta(path),
        gene_length_nt=GENE_NT,
        min_gene_coverage=0.80,
        max_ambiguity_frac=0.05,
    )
    qc_result.assert_conserved()
    assert "KE_bad" not in qc_result.kept

    metadata = [
        {
            "strain": strain,
            "country_code": "KEN",
            "year": 2015 + i,
            "subtype_called": "A1",
            "date_precision": "day",
            "gene_coverage": 0.95,
            "ambiguity_frac": 0.0,
        }
        for i, strain in enumerate(qc_result.kept)
    ]
    retained, summaries = stratified_subsample(
        metadata, strata=["country_code", "year_bin", "subtype_called"], max_per_stratum=10
    )

    assert set(retained) <= set(qc_result.kept)
    assert sum(summary.n_retained for summary in summaries) == len(retained)


def test_provenance_sidecar_is_written_next_to_its_output(tmp_path):
    output = tmp_path / "drm_table.tsv"
    output.write_text("strain\tmutation\n")
    record = build_provenance(
        stage="drm call",
        config_hash="abc123def456",
        snapshot_id="2026-11-05_a3f19c",
        seed=20261104,
        include_tool_versions=False,
    )
    sidecar = write_provenance(output, record)

    assert sidecar.name == "drm_table.tsv.provenance.json"
    assert sidecar.exists()
    assert "config_hash" in sidecar.read_text()
