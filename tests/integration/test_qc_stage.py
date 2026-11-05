"""The QC stage: ordered rules, gene-level findings, conservation, corpus recount."""

from __future__ import annotations

import random
from pathlib import Path

import pytest

from virotrace.curate.qc_stage import load_lab_strains, qualifying_genes, run_qc
from virotrace.drm.reference import CODON_TABLE
from virotrace.errors import ContractError
from virotrace.io.fasta import write_fasta
from virotrace.io.tables import read_tsv, write_tsv
from virotrace.settings import load_settings

POL_START, POL_END = 2085, 5096
SENSE = sorted(c for c, aa in CODON_TABLE.items() if aa != "*")
COVERAGE_COLUMNS = [
    "strain", "gene", "coverage", "n_called", "n_deleted", "n_too_ambiguous",
    "frameshift", "n_internal_stops", "usable", "identity", "sdrm_position_coverage",
]  # fmt: skip
META_COLUMNS = [
    "strain", "ambiguity_frac", "country_code", "region", "date_precision", "year", "clone_derived",
]  # fmt: skip


def _write_reference(path: Path) -> str:
    rng = random.Random(7)
    pol = "".join(rng.choice(SENSE) for _ in range(1003)) + "TAA"
    seq = "A" * (POL_START - 1) + pol + "C" * (9719 - POL_END)
    lines = ["LOCUS       K03455                9719 bp    DNA     linear   VRL 01-JAN-2000"]
    lines += ["VERSION     K03455.1", "FEATURES             Location/Qualifiers"]
    lines += ["     source          1..9719", "ORIGIN"]
    for i in range(0, len(seq), 60):
        chunk = seq[i : i + 60].lower()
        lines.append(f"{i + 1:>9} " + " ".join(chunk[j : j + 10] for j in range(0, len(chunk), 10)))
    lines.append("//")
    path.write_text("\n".join(lines) + "\n")
    return pol


def _gene_rows(
    strain, *, identity="0.9", pr="1.0", rt="1.0", rt_usable="true", rt_stops=0, rt_fs="false"
):
    def row(gene, cov_pos, usable="true", stops=0, fs="false"):
        return {
            "strain": strain, "gene": gene, "coverage": cov_pos or "0.0", "n_called": 0,
            "n_deleted": 0, "n_too_ambiguous": 0, "frameshift": fs, "n_internal_stops": stops,
            "usable": usable, "identity": identity, "sdrm_position_coverage": cov_pos,
        }  # fmt: skip

    return [
        row("PR", pr),
        row("RT", rt, rt_usable, rt_stops, rt_fs),
        {**row("IN", ""), "coverage": "0.0", "sdrm_position_coverage": ""},
    ]


def _hypermutate(pol: str) -> str:
    """G->A wherever the next reference base is A or G (the APOBEC3G context)."""
    out = list(pol)
    for i in range(len(pol) - 1):
        if pol[i] == "G" and pol[i + 1] in "GA":
            out[i] = "A"
    return "".join(out)


@pytest.fixture
def workspace(config_path, tmp_path):
    pol = _write_reference(tmp_path / "hxb2.gb")
    (tmp_path / "lab_strains.txt").write_text("# lab strains\nLAB1\n")
    settings = load_settings(
        config_path,
        "dev",
        {
            "paths": {
                "raw": str(tmp_path / "raw"),
                "interim": str(tmp_path / "interim"),
                "processed": str(tmp_path / "processed"),
                "reference": str(tmp_path),
                "results": str(tmp_path / "results"),
                "auspice": str(tmp_path / "auspice"),
                "reports": str(tmp_path / "reports"),
            }
        },
    )

    cases = {
        "GOOD": _gene_rows("GOOD"),
        "GARBAGE": _gene_rows(
            "GARBAGE", identity="0.5", rt_stops=40, rt_usable="false", rt_fs="true"
        ),
        "SHORT": _gene_rows("SHORT", pr="0.5", rt="0.5"),
        "LAB1": _gene_rows("LAB1"),
        "AMBIG": _gene_rows("AMBIG"),
        "BADRT": _gene_rows("BADRT", rt_usable="false", rt_stops=2),
        "HYPER": _gene_rows("HYPER"),
        "CLONE": _gene_rows("CLONE"),
    }
    write_tsv(
        tmp_path / "results" / "kenya" / "drm" / "gene_coverage.tsv",
        [row for rows in cases.values() for row in rows],
        COVERAGE_COLUMNS,
    )

    def meta(strain, *, ambiguity="0.0", clone="false", year="2015", precision="day"):
        return {
            "strain": strain, "ambiguity_frac": ambiguity, "country_code": "KEN",
            "region": "eastern_africa", "date_precision": precision, "year": year,
            "clone_derived": clone,
        }  # fmt: skip

    write_tsv(
        tmp_path / "interim" / "kenya" / "metadata.valid.tsv",
        [
            meta(
                s,
                ambiguity="0.2" if s == "AMBIG" else "0.0",
                clone="true" if s == "CLONE" else "false",
            )
            for s in cases
        ],
        META_COLUMNS,
    )
    write_fasta(
        tmp_path / "interim" / "kenya" / "aln.pol.unmasked.fasta",
        [(s, _hypermutate(pol) if s == "HYPER" else pol) for s in cases],
        wrap=0,
    )
    return settings, tmp_path


def test_each_rule_removes_only_what_it_should_and_the_books_balance(workspace):
    settings, tmp = workspace
    summary = run_qc(settings)

    assert summary["n_aligned"] == 8
    assert summary["waterfall"] == {"QC-7": 1, "QC-0": 1, "QC-1": 1, "QC-2": 1, "QC-5": 1}
    passed = {r["strain"]: r for r in read_tsv(tmp / "interim" / "kenya" / "qc.passed.tsv")}
    assert set(passed) == {"GOOD", "BADRT", "CLONE"}
    assert passed["GOOD"]["genes"] == "PR,RT"
    assert passed["BADRT"]["genes"] == "PR"  # a bad RT does not sink the protease
    assert summary["n_passed"] + sum(summary["waterfall"].values()) == summary["n_aligned"]


def test_ledger_names_the_rule_the_gene_and_the_observed_value(workspace):
    settings, tmp = workspace
    run_qc(settings)
    ledger = read_tsv(tmp / "results" / "kenya" / "qc" / "excluded.tsv")
    by = {(r["strain"], r["rule"], r["gene"]): r for r in ledger}

    assert by[("GARBAGE", "QC-0", "")]["observed"] == "0.5"
    assert by[("BADRT", "QC-3", "RT")]["observed"] == "2"
    assert ("HYPER", "QC-5", "") in by
    assert by[("AMBIG", "QC-2", "")]["observed"] == "0.2"
    assert ("LAB1", "QC-7", "") in by
    # a forced mis-alignment has stops and frameshifts everywhere; they must not be logged
    assert not any(r["strain"] == "GARBAGE" and r["rule"] in {"QC-3", "QC-4"} for r in ledger)


def test_corpus_is_recounted_before_and_after_the_clone_policy(workspace):
    settings, _ = workspace
    summary = run_qc(settings)
    assert summary["corpus_after_qc"]["n_total"] == 3
    assert summary["clone_policy"] == "exclude"
    assert summary["corpus_after_qc_and_clone_policy"]["n_total"] == 2  # CLONE leaves


def test_hypermutation_can_be_flag_only(config_path, workspace):
    settings, tmp = workspace
    relaxed = load_settings(
        config_path,
        "dev",
        {**{"paths": settings.paths}, "qc": {"hypermutation": {"action": "flag_only"}}},
    )
    summary = run_qc(relaxed)
    assert summary["waterfall"]["QC-5"] == 0
    passed = {r["strain"]: r for r in read_tsv(tmp / "interim" / "kenya" / "qc.passed.tsv")}
    assert passed["HYPER"]["hypermutated"] == "true"


def test_position_coverage_threshold_is_configurable(config_path, workspace):
    settings, _ = workspace
    loose = load_settings(
        config_path,
        "dev",
        {**{"paths": settings.paths}, "qc": {"min_class_position_coverage": 0.5}},
    )
    summary = run_qc(loose)
    assert summary["waterfall"]["QC-1"] == 0  # SHORT (0.5) now qualifies


def test_missing_coverage_table_points_at_the_previous_stage(config_path, tmp_path):
    settings = load_settings(
        config_path,
        "dev",
        {
            "paths": {
                **load_settings(config_path, "dev").paths,
                "interim": str(tmp_path),
                "results": str(tmp_path),
            }
        },
    )
    with pytest.raises(ContractError, match="drm call"):
        run_qc(settings)


def test_helpers(tmp_path):
    assert load_lab_strains(tmp_path / "absent.txt") == frozenset()
    (tmp_path / "l.txt").write_text("# c\nA\n\nB\n")
    assert load_lab_strains(tmp_path / "l.txt") == {"A", "B"}
    rows = _gene_rows("X", pr="1.0", rt="0.9")
    assert qualifying_genes(rows, min_position_coverage=1.0, min_gene_coverage=0.8) == ["PR"]
    assert qualifying_genes(rows, min_position_coverage=0.9, min_gene_coverage=0.8) == ["PR", "RT"]
