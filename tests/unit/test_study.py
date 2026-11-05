from virotrace.acquire.genbank_parse import parse_genbank
from virotrace.curate.study import study_key, treatment_from_text

GB = """LOCUS       AB000001                  20 bp    DNA     linear   VRL 01-JAN-2020
DEFINITION  HIV-1 isolate X from Kenya pol protein (pol) gene, partial cds.
ACCESSION   AB000001
VERSION     AB000001.1
DBLINK      BioProject: PRJNA123456
REFERENCE   1  (bases 1 to 20)
  AUTHORS   Smith,J., Doe,A. and Roe,B.
  TITLE     Drug resistance mutations among adults experiencing treatment
            failure in western Kenya
  JOURNAL   J. Acquir. Immune Defic. Syndr. 80 (1), 1-9 (2019)
   PUBMED   30000001
REFERENCE   2  (bases 1 to 20)
  AUTHORS   Smith,J.
  TITLE     Direct Submission
FEATURES             Location/Qualifiers
     source          1..20
                     /organism="Human immunodeficiency virus 1"
                     /country="Kenya"
                     /collection_date="2017"
ORIGIN
        1 acgtacgtac gtacgtacgt
//
"""


def test_reference_and_bioproject_are_read_from_the_first_reference_only():
    (rec,) = list(parse_genbank(GB))
    f = rec.fields
    assert f["pubmed"] == "30000001" and f["bioproject"] == "PRJNA123456"
    assert f["study_title"].endswith("western Kenya") and "Direct" not in f["study_title"]
    assert f["study_authors"].startswith("Smith")
    assert study_key(f) == "PMID:30000001"
    assert treatment_from_text(f) == "experienced"


def test_study_key_falls_back_and_never_invents():
    assert study_key({"bioproject": "PRJNA1"}) == "BP:PRJNA1"
    assert (
        study_key({"study_authors": "Doe,A.", "study_title": "Some Title!"}) == "REF:Doe|some-title"
    )
    assert study_key({}) == "unknown_study"


def test_treatment_hint_is_conservative():
    assert (
        treatment_from_text({"study_title": "Transmitted drug resistance in naive adults"})
        == "naive"
    )
    assert treatment_from_text({"study_title": "HIV diversity in Kenya"}) == "unknown"
    # contradictory statements stay unknown rather than guessing
    assert treatment_from_text({"study_title": "naive and treatment failure cohorts"}) == "unknown"
