"""Cross-source deduplication and the latest-snapshot helper."""

from __future__ import annotations

import pytest

from virotrace.acquire.base import RawRecord
from virotrace.acquire.fixture import FixtureAdapter
from virotrace.acquire.snapshot import latest_snapshot, take_snapshot
from virotrace.curate.dedup import base_accession, deduplicate
from virotrace.errors import ContractError


def rec(acc, seq="ACGT" * 10, source="genbank", **fields):
    return RawRecord(acc, seq, source, fields)


def test_base_accession_strips_version():
    assert base_accession("AB123456.2") == "AB123456"
    assert base_accession("AB123456") == "AB123456"


def test_same_accession_from_two_sources_collapses():
    result = deduplicate([rec("AB1.1", source="genbank"), rec("AB1", source="lanl")])
    assert len(result.kept) == 1
    assert [d.level for d in result.log] == ["accession"]
    result.assert_conserved()


def test_identical_sequence_under_different_accessions_collapses_to_the_cleaner_one():
    clean = rec("AB1", "ACGT" * 10)
    noisy = rec("AB2", "ACGT" * 9 + "ACNT")
    # different sequences, so they are NOT collapsed
    assert len(deduplicate([clean, noisy]).kept) == 2

    same_a = rec("AB1", "ACGT" * 10)
    same_b = rec("AB2", "ACGT" * 10)
    result = deduplicate([same_b, same_a])
    assert [r.accession for r in result.kept] == ["AB1"]  # stable tie-break
    assert result.log[0].level == "sequence"
    assert result.log[0].kept_instead == "AB1"


def test_prefers_fewer_ambiguous_bases_when_collapsing():
    # same accession, different content: the unambiguous copy wins
    good = rec("AB1.1", "ACGT" * 10)
    bad = rec("AB1.2", "ACGT" * 9 + "NNNN")
    result = deduplicate([bad, good])
    assert [r.accession for r in result.kept] == ["AB1.1"]


def test_same_isolate_different_sequence_is_flagged_not_collapsed():
    a = rec("AB1", "ACGT" * 10, isolate="K1", country="Kenya", collection_date="2015")
    b = rec("AB2", "ACGA" * 10, isolate="K1", country="Kenya: Nairobi", collection_date="2015")
    result = deduplicate([a, b])
    assert len(result.kept) == 2  # nothing is lost to an unlogged judgement call
    assert [d.level for d in result.log] == ["flag-isolate"]
    result.assert_conserved()


def test_deduplication_is_order_independent():
    records = [rec(f"A{i}", "ACGT" * 10 if i % 2 else "TTTT" * 10) for i in range(6)]
    forward = [r.accession for r in deduplicate(records).kept]
    backward = [r.accession for r in deduplicate(list(reversed(records))).kept]
    assert forward == backward


def test_latest_snapshot_picks_the_newest(tmp_path, records_toy_tsv):
    older = take_snapshot([FixtureAdapter(records_toy_tsv)], tmp_path, date="2026-10-01")
    newer = take_snapshot([FixtureAdapter(records_toy_tsv)], tmp_path, date="2026-11-05")
    assert latest_snapshot(tmp_path) == newer
    assert older != newer


def test_latest_snapshot_with_none_raises(tmp_path):
    with pytest.raises(ContractError, match="no snapshots"):
        latest_snapshot(tmp_path)
