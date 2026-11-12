"""Integration: fixture adapter -> immutable snapshot -> normalise -> corpus count."""

from __future__ import annotations

import json

import pytest

from virotrace.acquire.fixture import FixtureAdapter
from virotrace.acquire.snapshot import take_snapshot
from virotrace.curate.corpus import count_corpus
from virotrace.curate.normalise import normalise_records
from virotrace.errors import ContractError
from virotrace.io.fasta import read_fasta
from virotrace.io.manifest import make_snapshot_id, verify_manifest


def _snapshot(tmp_path, records_toy_tsv):
    return take_snapshot([FixtureAdapter(records_toy_tsv)], tmp_path, date="2026-11-05")


def test_snapshot_writes_manifest_with_verifiable_checksums(tmp_path, records_toy_tsv):
    snap = _snapshot(tmp_path, records_toy_tsv)
    manifest = json.loads((snap / "manifest.json").read_text())

    assert manifest["snapshot_id"] == snap.name
    assert manifest["sources"][0]["query"] == "fixture:records_toy.tsv"
    assert manifest["sources"][0]["retrieved"] == 9
    assert {f["path"] for f in manifest["files"]} == {"sequences.fasta", "records.tsv"}
    assert verify_manifest(snap / "manifest.json") == []


def test_tampering_with_a_snapshot_is_detected(tmp_path, records_toy_tsv):
    snap = _snapshot(tmp_path, records_toy_tsv)
    (snap / "sequences.fasta").write_text(">KE000001\nACGT\n")
    problems = verify_manifest(snap / "manifest.json")
    assert any("checksum mismatch: sequences.fasta" in p for p in problems)


def test_snapshots_are_immutable(tmp_path, records_toy_tsv):
    _snapshot(tmp_path, records_toy_tsv)
    with pytest.raises(ContractError, match="immutable"):
        _snapshot(tmp_path, records_toy_tsv)


class _Failing:
    name = "failing"
    access_tier = "open"

    def query_description(self):
        return "failing-query"

    def fetch(self):
        raise ContractError("NCBI requires a contact email")
        yield  # pragma: no cover


def test_a_failed_fetch_leaves_nothing_behind_so_the_retry_works(tmp_path, records_toy_tsv):
    with pytest.raises(ContractError, match="contact email"):
        take_snapshot([_Failing()], tmp_path, date="2026-11-05")
    assert list(tmp_path.iterdir()) == []


def test_an_empty_leftover_folder_does_not_block_a_new_snapshot(tmp_path, records_toy_tsv):
    leftover = tmp_path / make_snapshot_id("2026-11-05", ["fixture:records_toy.tsv"])
    leftover.mkdir()
    snap = _snapshot(tmp_path, records_toy_tsv)
    assert snap == leftover and (snap / "manifest.json").exists()


def test_snapshot_id_is_stable_and_order_independent():
    a = make_snapshot_id("2026-11-05", ["q1", "q2"])
    assert a == make_snapshot_id("2026-11-05", ["q2", "q1"])
    assert a != make_snapshot_id("2026-11-05", ["q1", "q3"])
    assert a.startswith("2026-11-05_")


def test_normalise_quarantines_instead_of_dropping(tmp_path, records_toy_tsv):
    snap = _snapshot(tmp_path, records_toy_tsv)
    records = list(FixtureAdapter(records_toy_tsv).fetch())
    result = normalise_records(records, snapshot_id=snap.name, max_record_error_rate=0.5)

    result.assert_conserved()
    assert result.n_in == 9
    quarantined = {q["accession"]: q["rule"] for q in result.quarantine}
    assert quarantined == {"XX000001": "N-2"}  # Atlantis is not a country

    # An undated record stays in (valid for the DRM layer); a pre-1985 date reads as unknown.
    by_strain = {r["strain"]: r for r in result.rows}
    assert by_strain["KE000004"]["date_precision"] == "unknown"
    assert by_strain["KE000005"]["date_precision"] == "unknown"
    assert all(r["treatment_status"] == "unknown" for r in result.rows)


def test_mass_quarantine_aborts():
    from virotrace.acquire.base import RawRecord

    bad = [RawRecord(f"A{i}", "ACGT", "fixture", {"country": "Atlantis"}) for i in range(10)]
    with pytest.raises(ContractError, match="refusing to continue"):
        normalise_records(bad, snapshot_id="s", max_record_error_rate=0.25)


def test_p1_gate_counts_dated_kenyan_sequences(tmp_path, records_toy_tsv):
    snap = _snapshot(tmp_path, records_toy_tsv)
    result = normalise_records(
        FixtureAdapter(records_toy_tsv).fetch(), snapshot_id=snap.name, max_record_error_rate=0.5
    )
    counts = count_corpus(result.rows)

    assert counts["n_ken"] == 5
    assert counts["by_precision"] == {"day": 1, "month": 1, "year": 1, "unknown": 2}
    assert counts["n_time_tree_eligible_strict"] == 2
    assert counts["n_time_tree_eligible_with_year_only"] == 3
    assert counts["year_range"] == (2014, 2016)


def test_fasta_headers_are_strain_only(tmp_path, records_toy_tsv):
    snap = _snapshot(tmp_path, records_toy_tsv)
    strains = [s for s, _ in read_fasta(snap / "sequences.fasta")]
    assert strains == sorted(strains)
    assert "KE000001" in strains
