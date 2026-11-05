"""GenBank flat-file parsing and the Entrez adapter, fully offline."""

from __future__ import annotations

import json

import pytest

from virotrace.acquire.genbank import GenBankAdapter, build_term
from virotrace.acquire.genbank_parse import parse_genbank
from virotrace.curate.normalise import normalise_records
from virotrace.errors import ConfigError, ContractError

FIXTURE = "genbank_toy.gb"


@pytest.fixture
def gb_text(request):
    from pathlib import Path

    return (Path(__file__).parents[1] / "fixtures" / FIXTURE).read_text()


def test_parser_extracts_fields_and_sequence(gb_text):
    records = list(parse_genbank(gb_text))
    assert [r.accession for r in records] == ["KE100001.1", "KE100002.1"]

    first = records[0]
    assert first.sequence.startswith("ACGTACGTAC") and len(first.sequence) == 60
    assert first.fields["country"] == "Kenya: Nairobi"
    assert first.fields["collection_date"] == "12-Mar-2015"
    assert first.fields["isolate"] == "K1"
    # depositor-reported subtype, read from a wrapped /note
    assert first.fields["subtype"] == "A1"


def test_missing_qualifiers_are_empty_not_guessed(gb_text):
    second = list(parse_genbank(gb_text))[1]
    assert second.fields["collection_date"] == ""
    assert second.fields["subtype"] == ""


def test_record_without_sequence_is_skipped():
    text = "LOCUS       X 0 bp\nACCESSION   X1\nVERSION     X1.1\nORIGIN\n//\n"
    assert list(parse_genbank(text)) == []


def test_parsed_records_normalise_end_to_end(gb_text):
    result = normalise_records(parse_genbank(gb_text), snapshot_id="s")
    by_strain = {r["strain"]: r for r in result.rows}
    assert by_strain["KE100001.1"]["date"] == "2015-03-12"
    assert by_strain["KE100001.1"]["admin1"] == "Nairobi"
    assert by_strain["KE100002.1"]["date_precision"] == "unknown"


# ------------------------------------------------------------------- the adapter
class FakeEntrez:
    """Stands in for NCBI; records every call so tests can assert on them."""

    def __init__(self, gb_text, *, count=None, fail_first=0):
        self.gb_text = gb_text
        self.calls: list[tuple[str, dict]] = []
        self.count = count
        self.fail_first = fail_first

    def __call__(self, url, params):
        self.calls.append((url.rsplit("/", 1)[-1], params))
        if self.fail_first > 0:
            self.fail_first -= 1
            raise OSError("simulated network error")
        if url.endswith("esearch.fcgi"):
            ids = ["101", "102"]
            return json.dumps(
                {"esearchresult": {"idlist": ids, "count": str(self.count or len(ids))}}
            )
        return self.gb_text


def _adapter(fake, **kw):
    return GenBankAdapter(email="me@example.org", transport=fake, sleep=lambda s: None, **kw)


def test_adapter_fetches_and_parses(gb_text):
    fake = FakeEntrez(gb_text)
    records = list(_adapter(fake).fetch())
    assert len(records) == 2
    endpoints = [name for name, _ in fake.calls]
    assert endpoints == ["esearch.fcgi", "efetch.fcgi"]


def test_query_is_recorded_verbatim():
    adapter = _adapter(FakeEntrez(""))
    assert adapter.query_description() == (
        "db=nuccore term=txid11676[Organism] AND Kenya[Country] AND 1000:10000[SLEN]"
    )
    assert build_term(country="Uganda", min_length=900) == (
        "txid11676[Organism] AND Uganda[Country] AND 900:10000[SLEN]"
    )


def test_email_is_required(monkeypatch):
    monkeypatch.delenv("VIROTRACE_NCBI_EMAIL", raising=False)
    adapter = GenBankAdapter(transport=FakeEntrez(""), sleep=lambda s: None)
    with pytest.raises(ConfigError, match="VIROTRACE_NCBI_EMAIL"):
        list(adapter.fetch())


def test_truncated_result_set_is_refused_not_silently_accepted(gb_text):
    fake = FakeEntrez(gb_text, count=5000)
    with pytest.raises(ContractError, match="silently truncating"):
        list(_adapter(fake).fetch())


def test_transient_errors_are_retried(gb_text):
    fake = FakeEntrez(gb_text, fail_first=2)
    assert len(list(_adapter(fake).fetch())) == 2


def test_persistent_errors_raise_after_the_retry_budget(gb_text):
    fake = FakeEntrez(gb_text, fail_first=99)
    with pytest.raises(ContractError, match="failed after 3 attempts"):
        list(_adapter(fake).fetch())


def test_batching(gb_text):
    fake = FakeEntrez(gb_text)
    list(_adapter(fake, batch_size=1).fetch())
    assert [n for n, _ in fake.calls].count("efetch.fcgi") == 2


def test_api_key_raises_the_rate_limit():
    with_key = GenBankAdapter(email="a@b.c", api_key="k", transport=FakeEntrez(""))
    without = GenBankAdapter(email="a@b.c", transport=FakeEntrez(""))
    assert with_key.min_interval < without.min_interval


def test_geo_loc_name_is_accepted_when_country_absent(gb_text):
    modern = gb_text.replace("/country=", "/geo_loc_name=")
    records = list(parse_genbank(modern))
    assert records[0].fields["country"] == "Kenya: Nairobi"


def test_dev_cap_spreads_across_the_result_list_not_the_front():
    from virotrace.acquire.genbank import _evenly_spaced

    ids = [str(i) for i in range(1000)]
    picked = _evenly_spaced(ids, 10)
    assert len(picked) == 10 and picked[0] == "0"
    assert int(picked[-1]) >= 900  # reaches the end of the list
    assert picked == _evenly_spaced(ids, 10)  # deterministic
    assert _evenly_spaced(ids, 5000) == ids


_CLONE_RECORD = """LOCUS       AB000001                 60 bp    DNA     linear   VRL 24-JUL-2016
DEFINITION  Human immunodeficiency virus 1 proviral DNA, complete genome,
            clone: pPRD320-16G41.
ACCESSION   AB000001
VERSION     AB000001.1
FEATURES             Location/Qualifiers
     source          1..60
                     /organism="Human immunodeficiency virus 1"
                     /proviral
                     /isolate="HH8793"
                     /clone="pPRD320-16G41"
                     /geo_loc_name="Kenya"
ORIGIN
        1 acgtacgtac acgtacgtac acgtacgtac acgtacgtac acgtacgtac acgtacgtac
//
"""


def test_clone_and_definition_are_captured_from_a_real_style_record():
    (record,) = list(parse_genbank(_CLONE_RECORD))
    assert record.fields["definition"].endswith("clone: pPRD320-16G41.")
    assert "complete genome" in record.fields["definition"]  # wrapped line rejoined
    assert record.fields["clone"] == "pPRD320-16G41"
    assert record.fields["proviral"] == "yes"


def test_clone_derived_flag_is_set_for_clones_and_not_for_ordinary_isolates(gb_text):
    from virotrace.curate.normalise import is_clone_derived

    (clone,) = list(parse_genbank(_CLONE_RECORD))
    assert is_clone_derived(clone.fields)
    assert not any(is_clone_derived(r.fields) for r in parse_genbank(gb_text))
    # a flag on the definition text alone, e.g. when the qualifier is missing
    assert is_clone_derived({"definition": "HIV-1 isolate X infectious molecular clone"})
    assert not is_clone_derived({"definition": "HIV-1 isolate K1 pol gene, partial cds."})
