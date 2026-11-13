"""The drm-prevalence stage: SDRM prevalence per drug class, with honest denominators.

Population: sequences that passed QC (``qc.passed.tsv``) and survive the clone policy.
Subsampling for the tree is deliberately not applied; it would shape prevalence.

Per class, a sequence is in the denominator only if the class's gene qualified (usable,
and covering the class's SDRM positions); it is in the numerator if it carries at least
one SDRM at a position of that class. Missing data never counts as susceptible.

Two lists are reported, labelled in every row:
    headline     the list in ``drm.sdrm.list`` (the original WHO 2009 list)
    sensitivity  ``drm.sdrm.sensitivity_list`` (Stanford), with its own denominators,
                 because it covers more positions than the 2009 list
    sensitivity_headline_denominator
                 Stanford's hits among the headline's own denominators; comparable to the
                 headline, and a lower bound for the extra positions
A class with no positions on a list (INSTI on the 2009 list) gets a ``no_list_positions``
row instead of an estimate: zero numerators there would be absence of a definition, not
absence of resistance.

Outputs: results/<scope>/drm/prevalence.tsv, trend_tests.tsv, study_summary.tsv and
prevalence_summary.json.
Trend tests share one caveat: sequence deposition changes over time, so a trend can
reflect which studies were published when, not epidemic change.
"""

from __future__ import annotations

import json
from collections import defaultdict
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

from ..curate.clones import apply_population_policies
from ..curate.study import load_study_annotations
from ..errors import ContractError
from ..io.tables import read_tsv, write_tsv
from ..settings import Settings
from .classes import DRUG_CLASSES, drug_class_of, gene_of_class
from .prevalence import cochran_armitage_trend, estimate_prevalence
from .sdrm import load_sdrm_list

PREVALENCE_COLUMNS = [
    "list_role",
    "list_version",
    "stratum_kind",
    "stratum",
    "drug_class",
    "numerator",
    "denominator",
    "proportion",
    "ci_low",
    "ci_high",
    "status",
]
TREND_COLUMNS = [
    "list_role",
    "list_version",
    "drug_class",
    "n_bins",
    "first_bin",
    "last_bin",
    "z",
    "p_two_sided",
    "note",
]
NO_LIST_POSITIONS = "no_list_positions"


def year_bin(year: int, width: int) -> str:
    start = (year // width) * width
    return f"{start}-{start + width - 1}"


def _class_positions(entries: Iterable[tuple[str, int]], drug_class: str) -> set[int]:
    gene = gene_of_class(drug_class)
    return {pos for (g, pos) in entries if g == gene and drug_class_of(g, pos) == drug_class}


def _qualifies(
    rows: list[dict[str, str]], gene: str, sens: bool, threshold: float, min_gene_coverage: float
) -> bool:
    for row in rows:
        if row["gene"] != gene or row["usable"].strip().lower() != "true":
            continue
        column = "sdrm_position_coverage_sensitivity" if sens else "sdrm_position_coverage"
        value = row.get(column, "")
        if value.strip():
            return float(value) >= threshold
        return (
            float(row["coverage"]) >= min_gene_coverage
        )  # no list positions in this gene: plain coverage
    return False


STUDY_COLUMNS = [
    "study_key",
    "study_title",
    "n_sequences",
    "year_min",
    "year_max",
    "treatment_hint",
    "annotation",
    "NRTI",
    "NNRTI",
    "PI",
]


def _study_summary(
    population: Mapping[str, Mapping[str, Any]],
    coverage: Mapping[str, list[dict[str, str]]],
    hits: Mapping[tuple[str, str], set[str]],
    qc: Any,
    headline_list: Any,
) -> list[dict[str, Any]]:
    """Per study: size, years, stated treatment context and headline-list prevalence.

    The diagnostic that shows whether a stratum is one study. ``NRTI`` etc. read
    ``resistant/denominator``; a class with no qualifying sequences reads ``0/0``.
    """
    by_study: dict[str, list[str]] = defaultdict(list)
    for strain, member in population.items():
        by_study[member.get("study_key", "").strip() or "unknown_study"].append(strain)
    rows: list[dict[str, Any]] = []
    for key, strains in by_study.items():
        years = [int(population[s]["year"]) for s in strains if population[s]["year"].strip()]
        statuses: defaultdict[str, int] = defaultdict(int)
        for s in strains:
            statuses[population[s].get("treatment_hint", "").strip() or "unknown"] += 1
        top_status = max(statuses, key=lambda k: statuses[k])
        titles = [
            population[s].get("study_title", "")
            for s in strains
            if population[s].get("study_title")
        ]
        row: dict[str, Any] = {
            "study_key": key,
            "study_title": titles[0] if titles else "",
            "n_sequences": len(strains),
            "year_min": min(years, default=""),
            "year_max": max(years, default=""),
            "treatment_hint": top_status,
            "annotation": population[strains[0]].get("annotation", "") or "none",
        }
        for drug_class in ("NRTI", "NNRTI", "PI"):
            gene = gene_of_class(drug_class)
            den = [
                s
                for s in strains
                if _qualifies(
                    coverage[s], gene, False, qc.min_class_position_coverage, qc.min_gene_coverage
                )
            ]
            row[drug_class] = f"{sum(s in hits[('headline', drug_class)] for s in den)}/{len(den)}"
        rows.append(row)
    return sorted(rows, key=lambda r: -r["n_sequences"])


def run_prevalence(settings: Settings, *, repo_root: Path | str = ".") -> dict[str, Any]:
    root = Path(repo_root)
    scope = settings.scope
    interim = root / settings.paths["interim"] / scope
    results = root / settings.paths["results"] / scope
    cfg = settings.prevalence
    qc = settings.qc

    paths = {
        "drm call": results / "drm" / "drm_table.tsv",
        "drm call (coverage)": results / "drm" / "gene_coverage.tsv",
        "normalise": interim / "metadata.valid.tsv",
        "qc": interim / "qc.passed.tsv",
    }
    for stage, path in paths.items():
        if not path.exists():
            raise ContractError(f"{path} not found; run `virotrace {stage.split(' (')[0]}` first")

    reference_dir = root / settings.paths["reference"]
    headline_list = load_sdrm_list(
        reference_dir / settings.drm.sdrm.list, version=settings.drm.sdrm.list_version
    )
    # (role, list the hits come from, qualify on sensitivity coverage?, hits from sensitivity
    # columns?, list whose positions decide whether a class exists and who qualifies)
    lists: list[tuple[str, Any, bool, bool, Any]] = [
        ("headline", headline_list, False, False, headline_list)
    ]
    if settings.drm.sdrm.sensitivity_list:
        sens_list = load_sdrm_list(
            reference_dir / settings.drm.sdrm.sensitivity_list,
            version=settings.drm.sdrm.sensitivity_list_version or "sensitivity",
        )
        lists.append(("sensitivity", sens_list, True, True, sens_list))
        # Same people as the headline, Stanford's hits: a lower bound, because positions a
        # sequence does not cover cannot be called. Comparable to the headline by construction.
        lists.append(("sensitivity_headline_denominator", sens_list, False, True, headline_list))

    passed = {r["strain"] for r in read_tsv(paths["qc"])}
    meta_rows = [r for r in read_tsv(paths["normalise"]) if r["strain"] in passed]
    annotations = load_study_annotations(reference_dir / cfg.study_annotations)
    calls_path = interim / "subtype_calls.tsv"
    called = (
        {r["strain"]: r["subtype_called"] for r in read_tsv(calls_path)}
        if calls_path.exists()
        else None
    )
    typed: list[dict[str, Any]] = []
    for r in meta_rows:
        item: dict[str, Any] = {
            **r,
            "clone_derived": r.get("clone_derived", "").strip().lower() == "true",
        }
        note = annotations.get(item.get("study_key", ""))
        # The keyword hint is shown for review but never stratifies: it mislabelled a study
        # with 61% NRTI resistance as naive. Only hand-curated annotations do.
        if called is not None:
            item["subtype_called"] = called.get(item["strain"], "unassigned")
        item["treatment_hint"] = item.get("treatment_status", "")
        item["annotation"] = note["treatment_context"] if note else ""
        item["treatment_status"] = (
            note["treatment_context"]
            if note and note["treatment_context"] != "unknown"
            else "unannotated"
        )
        typed.append(item)
    kept, removed = apply_population_policies(
        typed,
        clone_policy=qc.clone_policy,
        individual_policy=qc.individual_policy,
        patterns=qc.individual_id_patterns,
    )
    population = {r["strain"]: r for r in kept}

    coverage: dict[str, list[dict[str, str]]] = defaultdict(list)
    for crow in read_tsv(paths["drm call (coverage)"]):
        if crow["strain"] in population:
            coverage[crow["strain"]].append(crow)

    hits: dict[tuple[str, str], set[str]] = defaultdict(set)  # (role, class) -> strains
    for hrow in read_tsv(paths["drm call"]):
        if hrow["strain"] not in population:
            continue
        for role, _list, _qual, sens, _entries in lists:
            column = "is_sdrm_sensitivity" if sens else "is_sdrm"
            if hrow[column].strip().lower() != "true":
                continue
            hits[(role, drug_class_of(hrow["gene"], int(hrow["position"])))].add(hrow["strain"])

    width = cfg.year_bin_width
    out_rows: list[dict[str, Any]] = []
    trend_rows: list[dict[str, Any]] = []
    summary: dict[str, Any] = {
        "scope": scope,
        "n_population": len(population),
        "n_removed_by_population_policy": removed,
        "lists": {},
    }

    def strata_of(row: Mapping[str, Any]) -> list[tuple[str, str]]:
        found = [("overall", "all")]
        if row["year"].strip():
            found.append(("year_bin", year_bin(int(row["year"]), width)))
        found.append(("treatment_status", row.get("treatment_status", "").strip() or "unannotated"))
        if cfg.by_country:
            found.append(("country", row.get("country", "").strip() or "unknown"))
        if cfg.also_by_subtype and called is not None:
            found.append(("subtype_called", row.get("subtype_called", "") or "unassigned"))
        elif cfg.also_by_subtype:
            found.append(
                ("subtype_reported", row.get("subtype_reported", "").strip() or "unassigned")
            )
        return found

    for role, sdrm_list, qual_sens, _hit_sens, qual_list in lists:
        listed = {}
        for drug_class in DRUG_CLASSES:
            positions = _class_positions(qual_list.entries, drug_class)
            listed[drug_class] = len(positions)
            if not positions:
                out_rows.append(
                    {
                        "list_role": role,
                        "list_version": sdrm_list.version,
                        "stratum_kind": "overall",
                        "stratum": "all",
                        "drug_class": drug_class,
                        "numerator": "",
                        "denominator": "",
                        "proportion": "",
                        "ci_low": "",
                        "ci_high": "",
                        "status": NO_LIST_POSITIONS,
                    }
                )
                continue
            gene = gene_of_class(drug_class)
            counts: dict[tuple[str, str], list[int]] = defaultdict(lambda: [0, 0])
            for strain, member in population.items():
                if not _qualifies(
                    coverage[strain],
                    gene,
                    qual_sens,
                    qc.min_class_position_coverage,
                    qc.min_gene_coverage,
                ):
                    continue
                carries = strain in hits[(role, drug_class)]
                for key in strata_of(member):
                    counts[key][1] += 1
                    counts[key][0] += carries
            for (kind, label), (num, den) in sorted(counts.items()):
                est = estimate_prevalence(
                    numerator=num,
                    denominator=den,
                    drug_class=drug_class,
                    min_stratum_n=cfg.min_stratum_n,
                    ci_level=cfg.ci_level,
                )
                out_rows.append(
                    {
                        "list_role": role,
                        "list_version": sdrm_list.version,
                        "stratum_kind": kind,
                        "stratum": label,
                        "drug_class": drug_class,
                        "numerator": num,
                        "denominator": den,
                        "proportion": "" if est.proportion is None else round(est.proportion, 5),
                        "ci_low": "" if est.ci_low is None else round(est.ci_low, 5),
                        "ci_high": "" if est.ci_high is None else round(est.ci_high, 5),
                        "status": est.status,
                    }
                )
            if cfg.trend_test == "cochran_armitage":
                bins = sorted(
                    (int(label.split("-")[0]), c)
                    for (kind, label), c in counts.items()
                    if kind == "year_bin" and c[1] >= cfg.min_stratum_n
                )
                base = {
                    "list_role": role,
                    "list_version": sdrm_list.version,
                    "drug_class": drug_class,
                }
                if len(bins) < 3:
                    trend_rows.append(
                        {**base, "n_bins": len(bins), "first_bin": "", "last_bin": "", "z": "",
                         "p_two_sided": "", "note": "fewer than 3 year bins with enough sequences"}
                    )  # fmt: skip
                else:
                    z, p = cochran_armitage_trend(
                        [c[0] for _, c in bins],
                        [c[1] for _, c in bins],
                        [float(s) for s, _ in bins],
                    )
                    trend_rows.append(
                        {**base, "n_bins": len(bins), "first_bin": bins[0][0], "last_bin": bins[-1][0],
                         "z": round(z, 4), "p_two_sided": f"{p:.3g}",
                         "note": "composition of sampling changes over time; not proof of epidemic change"}
                    )  # fmt: skip
        summary["lists"][role] = {
            "version": sdrm_list.version,
            "positions_per_class": listed,
            "overall": {
                r["drug_class"]: {
                    k: r[k]
                    for k in (
                        "numerator",
                        "denominator",
                        "proportion",
                        "ci_low",
                        "ci_high",
                        "status",
                    )
                }
                for r in out_rows
                if r["list_role"] == role and r["stratum_kind"] == "overall"
            },
        }

    study_rows = _study_summary(population, coverage, hits, qc, headline_list)
    out = results / "drm"
    write_tsv(out / "study_summary.tsv", study_rows, STUDY_COLUMNS)
    summary["n_studies"] = len(study_rows)
    summary["largest_studies"] = [
        {k: r[k] for k in ("study_key", "n_sequences", "treatment_hint", "annotation")}
        for r in study_rows[:5]
    ]
    write_tsv(out / "prevalence.tsv", out_rows, PREVALENCE_COLUMNS)
    write_tsv(out / "trend_tests.tsv", trend_rows, TREND_COLUMNS)
    (out / "prevalence_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    return summary
