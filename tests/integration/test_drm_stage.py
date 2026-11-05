"""The drm-call stage on a synthetic, reference-projected alignment."""

from __future__ import annotations

import json
import random
import shutil
from pathlib import Path

import pytest

from virotrace.drm.reference import CODON_TABLE, translate
from virotrace.drm.stage import run_drm_call
from virotrace.errors import ContractError
from virotrace.io.fasta import write_fasta
from virotrace.io.tables import read_tsv, write_tsv
from virotrace.settings import load_settings

POL_START, POL_END = 2085, 5096
SENSE = sorted(c for c, aa in CODON_TABLE.items() if aa != "*")
REPO = Path(__file__).resolve().parents[2]


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


@pytest.fixture
def workspace(config_path, tmp_path):
    pol = _write_reference(tmp_path / "hxb2.gb")
    for name in ("sdrm_who.tsv", "sdrm_stanford_2026-09-25.tsv"):
        shutil.copy(REPO / "config" / "reference" / name, tmp_path / name)
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
    return settings, pol, tmp_path


def _rt_codon_start(settings, position):
    return settings.reference.genes.rt[0] - POL_START + 3 * (position - 1)


def _build(settings, pol, tmp):
    rt_ref_aa = translate(pol[_rt_codon_start(settings, 184) : _rt_codon_start(settings, 184) + 3])
    alt = "I" if rt_ref_aa == "V" else "V"  # both V and I are SDRMs at RT 184
    alt_codon = next(c for c in SENSE if CODON_TABLE[c] == alt)

    def edit(seq, index, new):
        return seq[:index] + new + seq[index + len(new) :]

    mutant = edit(pol, _rt_codon_start(settings, 184), alt_codon)
    rt_start = settings.reference.genes.rt[0] - POL_START
    rt_only = "-" * rt_start + pol[rt_start : rt_start + 1680] + "-" * (len(pol) - rt_start - 1680)
    frameshifted = edit(pol, _rt_codon_start(settings, 200), "-")
    rng = random.Random(99)
    garbage = "".join(rng.choice("ACGT") for _ in range(len(pol)))  # forced mis-alignment
    rows = [
        ("GARBAGE", garbage),
        ("WT", pol),
        ("MUT", mutant),
        ("RTONLY", rt_only),
        ("FRAME", frameshifted),
        ("NOPOL", "-" * len(pol)),
    ]
    write_fasta(tmp / "interim" / "kenya" / "aln.pol.unmasked.fasta", rows, wrap=0)
    write_tsv(
        tmp / "interim" / "kenya" / "insertions.tsv", [], ["strain", "after_ref_pos", "length"]
    )
    return alt


def test_stage_reports_coverage_mutations_and_disqualifications(workspace):
    settings, pol, tmp = workspace
    alt = _build(settings, pol, tmp)
    summary = run_drm_call(settings)

    assert summary["n_sequences"] == 6
    assert summary["n_with_almost_no_pol"] == 1  # NOPOL
    assert summary["sdrm_list_version"] == "2009"

    coverage = {
        (r["strain"], r["gene"]): r
        for r in read_tsv(tmp / "results" / "kenya" / "drm" / "gene_coverage.tsv")
    }
    assert coverage[("WT", "RT")]["coverage"] == "1.0"
    assert coverage[("RTONLY", "PR")]["coverage"] == "0.0"  # not sequenced: no call
    assert coverage[("RTONLY", "RT")]["coverage"] == "1.0"
    assert coverage[("FRAME", "RT")]["frameshift"] == "true"
    assert coverage[("FRAME", "RT")]["usable"] == "false"
    assert coverage[("FRAME", "PR")]["usable"] == "true"  # other genes unaffected
    assert coverage[("NOPOL", "RT")]["coverage"] == "0.0"

    table = read_tsv(tmp / "results" / "kenya" / "drm" / "drm_table.tsv")
    assert [(r["strain"], r["gene"], r["position"], r["residues"]) for r in table] == [
        ("MUT", "RT", "184", alt)
    ]
    assert table[0]["is_sdrm"] == "true" and table[0]["sdrm_list_version"] == "2009"

    # one gene row per (strain, gene): 6 sequences x 3 genes
    assert len(coverage) == 18
    # identity separates real pol from a forced mis-alignment, and never reports an
    # uncovered gene as identity zero
    assert coverage[("WT", "RT")]["identity"] == "1.0"
    assert float(coverage[("GARBAGE", "RT")]["identity"]) < 0.5
    assert coverage[("RTONLY", "PR")]["identity"] == ""
    # position coverage is about the SDRM positions, not the whole gene
    assert coverage[("RTONLY", "RT")]["sdrm_position_coverage"] == "1.0"
    assert coverage[("RTONLY", "PR")]["sdrm_position_coverage"] == "0.0"
    assert coverage[("WT", "IN")]["sdrm_position_coverage"] == ""  # no IN positions in 2009
    stored = json.loads((tmp / "results" / "kenya" / "drm" / "drm_summary.json").read_text())
    assert stored["per_gene"]["RT"]["n_covering_all_sdrm_positions"] >= 3  # WT, MUT, RTONLY
    assert stored["per_gene"]["IN"]["n_sdrm_positions_on_list"] == 0
    assert sum(stored["per_gene"]["RT"]["identity_deciles"]) == 5  # all but NOPOL
    assert stored["per_gene"]["RT"]["identity_deciles"][9] >= 3  # WT, MUT, RTONLY

    assert stored["per_gene"]["RT"]["n_meeting_min_coverage"] == 3  # WT, MUT, RTONLY
    assert stored["per_gene"]["PR"]["n_meeting_min_coverage"] == 3  # WT, MUT, FRAME


def test_off_frame_insertion_from_the_side_table_disqualifies_the_gene(workspace):
    settings, pol, tmp = workspace
    _build(settings, pol, tmp)
    rt_mid = settings.reference.genes.rt[0] - POL_START + 601
    write_tsv(
        tmp / "interim" / "kenya" / "insertions.tsv",
        [{"strain": "WT", "after_ref_pos": rt_mid, "length": 2}],
        ["strain", "after_ref_pos", "length"],
    )
    run_drm_call(settings)
    coverage = {
        (r["strain"], r["gene"]): r
        for r in read_tsv(tmp / "results" / "kenya" / "drm" / "gene_coverage.tsv")
    }
    assert coverage[("WT", "RT")]["frameshift"] == "true"
    assert coverage[("WT", "PR")]["frameshift"] == "false"


def test_missing_alignment_points_at_the_previous_stage(workspace):
    settings, _, _ = workspace
    with pytest.raises(ContractError, match="align"):
        run_drm_call(settings)
