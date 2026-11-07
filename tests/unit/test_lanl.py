"""LANL adapter: headers, country codes, years, reproducibility record, dedup preference."""

from __future__ import annotations

import pytest

from virotrace.acquire.base import RawRecord
from virotrace.acquire.lanl import LanlAdapter, parse_header, parse_year
from virotrace.curate.dedup import deduplicate
from virotrace.curate.normalise import normalise_record
from virotrace.errors import ContractError


def test_header_fields_and_names_that_contain_dots():
    f = parse_header("A1.KE.00.02KE.NDK.12.KC123456")
    assert f == {
        "subtype": "A1", "country": "Kenya", "collection_date": "2000",
        "isolate": "02KE.NDK.12", "accession": "KC123456",
    }  # fmt: skip
    assert parse_header("01_AE.TH.x.name.-")["accession"] == ""  # no usable accession
    assert parse_header("C.ZA.2005.p1.AY772699")["collection_date"] == "2005"


def test_two_digit_years_and_unknown_years():
    assert [parse_year(x) for x in ("00", "98", "05", "x", "-", "60", "2012")] == [
        "2000", "1998", "2005", "", "", "", "2012",
    ]  # fmt: skip


def test_a_header_in_another_format_is_refused_with_the_fix():
    with pytest.raises(ContractError, match="default header"):
        parse_header("sequence_1 some description")


def test_records_are_ungapped_and_normalise_like_genbank_ones(tmp_path):
    (tmp_path / "pol.fasta").write_text(
        ">D.UG.99.isolate_a.AB000001\nAC-GT\nAC..\n>B.KE.10.noacc.-\nACGT\n"
    )
    first, second = list(LanlAdapter(tmp_path).fetch())
    assert (first.accession, first.sequence, first.source) == ("AB000001", "ACGTAC", "lanl")
    assert second.accession == "LANL_noacc"
    row = normalise_record(first, snapshot_id="s1")
    assert (row["country"], row["year"], row["subtype_reported"]) == ("Uganda", 1999, "D")


def test_unmapped_country_codes_survive_to_be_quarantined_and_logged(tmp_path):
    (tmp_path / "x.fa").write_text(">A1.ZZ.00.n.AB000002\nACGT\n")
    (record,) = LanlAdapter(tmp_path).fetch()
    assert record.fields["country"] == "ZZ"


def test_query_description_records_file_names_and_hashes_and_changes_with_content(tmp_path):
    path = tmp_path / "a.fasta"
    path.write_text(">A1.KE.00.n.AB000003\nACGT\n")
    before = LanlAdapter(tmp_path).query_description()
    assert before.startswith("lanl-download:a.fasta:")
    path.write_text(">A1.KE.00.n.AB000003\nACGA\n")
    assert LanlAdapter(tmp_path).query_description() != before


def test_missing_directory_or_files_say_what_to_do(tmp_path):
    with pytest.raises(ContractError, match=r"hiv\.lanl\.gov"):
        list(LanlAdapter(tmp_path / "nope").fetch())
    with pytest.raises(ContractError, match="no FASTA"):
        LanlAdapter(tmp_path).query_description()


def test_dedup_keeps_the_genbank_copy_of_a_record_found_in_both():
    lanl = RawRecord("AB1", "ACGTACGT", "lanl", {"subtype": "A1"})
    genbank = RawRecord("AB1.1", "ACGTACGT", "genbank", {"title": "a study"})
    for order in ([lanl, genbank], [genbank, lanl]):
        (kept,) = deduplicate(order).kept
        assert kept.source == "genbank"
