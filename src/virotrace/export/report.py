"""The surveillance report: one self-contained HTML file with its limitations section.

The report computes nothing new. It lays out what the stages wrote (prevalence, trend tests,
QC and tree-set counts, temporal gates, clusters) and always ends with the limitations. The
numbered limitations in ``docs/04-bioinformatics-methodology.md`` section 9 are read from that
file, so the document and the report cannot drift apart; the project-specific limitations
below are added to them.
"""

from __future__ import annotations

import json
import re
from datetime import UTC, datetime
from html import escape
from pathlib import Path
from typing import Any

from ..errors import ContractError
from ..io.tables import read_tsv
from ..settings import Settings

METHODS_DOC = Path("docs/04-bioinformatics-methodology.md")
NOT_RUN = "Not generated for this build."

EXTRA_LIMITATIONS = [
    "Country is the GenBank country qualifier. Multi-site trial cohorts (for example the RV217 "
    "and HVTN 703 isolates) may carry a single country tag for participants enrolled elsewhere; "
    "participant sites have not been verified, and these participants are a small share of "
    "individuals.",
    "Sequences come from GenBank only. LANL was not used, so records without a GenBank "
    "country qualifier, or whose country LANL curated, are not in the corpus.",
    "Prevalence is for published sequences, not transmitted drug resistance: most studies "
    "sampled people on or after treatment, and study treatment context comes from hand-curated "
    "annotations, with unreviewed studies shown as unannotated.",
    "SDRM list insertions and deletions (T69 insertions, K70 deletion) are not detected by the "
    "point-mutation caller, so prevalence is a lower bound for those positions.",
    "The SDRM list was transcribed from the published source and must be verified against it "
    "before a result is quoted.",
    "Subtype placement and the recombination screen are distance-based screens, not RDP or "
    "3SEQ; individual calls and breakpoints are provisional.",
    "Dated trees give rough timelines: estimated substitution rates are low, pooled "
    "multi-subtype trees fail the temporal-signal gate, and branch dates carry wide "
    "uncertainty. Only subtype trees that pass the gate are dated. Where augur's own rate "
    "estimate failed or was implausible, the rate is fixed at the gate's regression slope; the "
    "root dates of such trees mostly reflect the sampling window, are not estimates of when a "
    "subtype arrived, and are fit for relative timing only.",
    "Clusters are a descriptive TN93 statistic at a fixed threshold on the tree set only, and "
    "undercount; they are not transmission chains. Near-identical sequences within one study "
    "can be repeat sequences from one person that the isolate-name patterns did not recognise.",
]


def _read_json(path: Path) -> dict[str, Any] | None:
    return json.loads(path.read_text()) if path.exists() else None


def _table(headers: list[str], rows: list[list[Any]]) -> str:
    head = "".join(f"<th>{escape(str(h))}</th>" for h in headers)
    body = "".join(
        "<tr>" + "".join(f"<td>{escape(str(c))}</td>" for c in row) + "</tr>" for row in rows
    )
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def _pct(value: Any, digits: int = 1) -> str:
    return f"{float(value) * 100:.{digits}f}%" if value not in ("", None) else ""


def _estimate(row: dict[str, str]) -> str:
    if row["status"] == "estimated":
        return (
            f"{_pct(row['proportion'])} ({_pct(row['ci_low'])} to {_pct(row['ci_high'])}), "
            f"{row['numerator']}/{row['denominator']}"
        )
    if row["status"] == "n_too_small":
        return f"suppressed (n={row['denominator']})"
    return "no list positions"


def limitations(root: Path, scope: str = "kenya") -> list[str]:
    path = root / METHODS_DOC
    if not path.exists():
        raise ContractError(f"{path} not found; the report cannot be written without it")
    text = path.read_text()
    match = re.search(r"^## 9\. Stated limitations.*?\n(.*?)(?=^## |\Z)", text, re.M | re.S)
    items = re.findall(r"^\d+\.\s+(.*)$", match.group(1), re.M) if match else []
    if not items:
        raise ContractError(f"{path} has no numbered limitations in section 9")
    if scope != "kenya":
        # the methods document is written for the Kenya build; say what the report is about
        place = "Africa" if scope == "africa" else scope.capitalize()
        items = [
            i.replace("the Kenyan epidemic", f"the {place} epidemic").replace(
                "Kenya's epidemic is", f"the {place} epidemic is"
            )
            for i in items
        ]
    return items + EXTRA_LIMITATIONS


def _prevalence_sections(results: Path) -> str:
    rows = read_tsv(results / "drm" / "prevalence.tsv")
    classes = [
        c for c in ("NRTI", "NNRTI", "PI", "INSTI") if any(r["drug_class"] == c for r in rows)
    ]
    out = []
    for role, title in (
        ("headline", "Headline: WHO 2009 surveillance list"),
        (
            "sensitivity",
            "Sensitivity: current Stanford list (sequences covering all its positions)",
        ),
        (
            "sensitivity_headline_denominator",
            "Sensitivity: current Stanford list, same people as the headline",
        ),
    ):
        sub = [r for r in rows if r["list_role"] == role]
        if not sub:
            continue
        out.append(f"<h3>{escape(title)} (version {escape(sub[0]['list_version'])})</h3>")
        if role == "sensitivity":
            out.append(
                "<p>Only sequences covering every position of this list are counted, so the "
                "denominators are smaller than the headline's and the percentages are not "
                "comparable to it. Use the next table for a like-for-like comparison.</p>"
            )
        if role == "sensitivity_headline_denominator":
            out.append(
                "<p>The same people as the headline, with Stanford's mutations called at the "
                "positions each sequence covers: a lower bound, comparable to the headline.</p>"
            )
        kinds = list(dict.fromkeys(r["stratum_kind"] for r in sub))
        for kind in kinds:
            strata = list(dict.fromkeys(r["stratum"] for r in sub if r["stratum_kind"] == kind))
            table_rows = []
            for stratum in strata:
                cells = []
                for c in classes:
                    match = [
                        r
                        for r in sub
                        if r["stratum_kind"] == kind
                        and r["stratum"] == stratum
                        and r["drug_class"] == c
                    ]
                    cells.append(_estimate(match[0]) if match else "not assessed")
                table_rows.append([stratum, *cells])
            label = "Overall" if kind == "overall" else f"By {kind.replace('_', ' ')}"
            out.append(f"<h4>{escape(label)}</h4>" + _table(["Stratum", *classes], table_rows))
    trends = read_tsv(results / "drm" / "trend_tests.tsv")
    if trends:
        out.append("<h3>Trend tests (headline list)</h3>")
        out.append(
            _table(
                ["Class", "Bins", "First", "Last", "z", "p (two-sided)", "Note"],
                [
                    [
                        t["drug_class"],
                        t["n_bins"],
                        t["first_bin"],
                        t["last_bin"],
                        t["z"],
                        t["p_two_sided"],
                        t["note"],
                    ]
                    for t in trends
                    if t["list_role"] == "headline"
                ],
            )
        )
        out.append(
            "<p>A trend can reflect changing study composition rather than the epidemic "
            "(limitation 7).</p>"
        )
    return "\n".join(out)


def _data_section(results: Path) -> str:
    qc = _read_json(results / "qc" / "qc_summary.json")
    sub = _read_json(results / "qc" / "subsample_summary.json")
    out = []
    if qc:
        waterfall = qc.get("waterfall", {})
        out.append("<h3>Quality control</h3>")
        out.append(
            _table(
                ["Step", "Sequences removed"],
                [[k, v] for k, v in sorted(waterfall.items())]
                + [
                    [f"{k} (population policy)", v]
                    for k, v in qc.get("n_removed_by_population_policy", {}).items()
                ],
            )
        )
        out.append(
            f"<p>{qc.get('n_passed', 'n/a')} sequences passed QC. Population policies: clone "
            f"policy <code>{escape(str(qc.get('clone_policy')))}</code>, individual policy "
            f"<code>{escape(str(qc.get('individual_policy')))}</code>.</p>"
        )
        after = qc.get("corpus_after_qc_and_clone_policy", {})
        if after:
            out.append(
                f"<p>After the policies: {after.get('n_total')} individuals, "
                f"{after.get('n_time_tree_eligible_strict')} with day or month collection dates "
                f"({after.get('n_time_tree_eligible_with_year_only')} counting year-only).</p>"
            )
    if sub:
        out.append("<h3>Tree set</h3>")
        out.append(
            f"<p>{sub['waterfall'].get('retained')} individuals with strictly dated sequences, "
            f"{sub['year_span'][0]} to {sub['year_span'][1]}. Counts by subtype and year:</p>"
        )
        out.append(_table(["Subtype", "n"], [[k, v] for k, v in sub.get("by_subtype", {}).items()]))
        out.append(
            _table(["Year", "n"], [[k, v] for k, v in sorted(sub.get("by_year", {}).items())])
        )
    return "\n".join(out) or f"<p>{NOT_RUN}</p>"


def _gate_row(label: str, gate: dict[str, Any], time_tree: bool, root: str = "") -> list[Any]:
    return [
        label,
        gate.get("n", ""),
        gate.get("slope", ""),
        gate.get("r_squared", ""),
        gate.get("permutation_p", ""),
        "pass" if gate.get("passed") else "fail",
        ("time tree" + (f" ({root})" if root else "")) if time_tree else "divergence tree",
    ]


def _phylo_section(results: Path, subtypes: list[str]) -> str:
    gate = _read_json(results / "phylo" / "temporal_gate.json")
    if gate is None:
        return f"<p>{NOT_RUN}</p>"
    rows = [_gate_row("Pooled (all subtypes)", gate["pooled"], gate.get("time_tree", False))]
    for label, g in sorted(gate.get("by_subtype", {}).items()):
        rows.append(_gate_row(f"{label} within the pooled tree", g, False))
    for s in subtypes:
        sub = _read_json(results / "phylo" / s / "temporal_gate.json")
        if sub:
            rows.append(
                _gate_row(
                    f"{s} tree", sub["pooled"], sub.get("time_tree", False), sub.get("root", "")
                )
            )
    out = [
        "<h3>Temporal signal</h3>",
        _table(["Tree", "Tips", "Slope", "R squared", "Permutation p", "Gate", "Result"], rows),
        "<p>The gate is root-to-tip regression at the best root with a date-randomisation "
        "test. A tree that fails it is shown as a midpoint-rooted divergence tree.</p>",
    ]
    clock_rows = []
    for sub_label in subtypes:
        cs = _read_json(results / "phylo" / sub_label / "clock_sensitivity.json")
        if not cs:
            continue
        est = cs["estimated"]
        clock_rows.append(
            [
                sub_label,
                "fixed by the pipeline" if est.get("rate_fixed") else "estimated",
                f"{est['clock_rate']:.5f}",
                est["root_date"],
                "",
            ]
        )
        for rate, run in sorted(cs["fixed"].items(), key=lambda kv: float(kv[0])):
            clock_rows.append(
                [sub_label, "fixed", rate, run["root_date"], f"{run['root_shift_years']:+.1f}"]
            )
    if clock_rows:
        out += [
            "<h3>Clock-rate sensitivity</h3>",
            _table(
                ["Tree", "Rate", "Substitutions/site/year", "Root date", "Shift (years)"],
                clock_rows,
            ),
            "<p>Root dates move with the assumed rate. Quote a timeline as a range across these "
            "rates, not as the point estimate.</p>",
        ]
    st = _read_json(results / "phylo" / "subtype_tree" / "subtype_tree_summary.json")
    if st:
        c = st["concordance"]
        out += [
            "<h3>Second subtype call (tree-based)</h3>",
            f"<p>Each of {st['n_queries']} tree-set sequences was placed in a maximum-likelihood "
            f"tree of {st['n_references']} reference sequences. A call needs a clade of one "
            f"subtype and {st['min_support']:g}% bootstrap support on the joining branch. "
            f"Agree {c.get('agree', 0)}, disagree {c.get('disagree', 0)}, not comparable "
            f"{c.get('not_comparable', 0)}, unresolved {c.get('unresolved', 0)}. Unresolved is "
            "common because few references represent each subtype, so this check can confirm a "
            "call but not rule one out.</p>",
        ]
    cl = _read_json(results / "phylo" / "clusters_summary.json")
    if cl:
        s = cl["sdrm"]
        adj = cl.get("sdrm_cluster_odds_ratio_adjusted_for_treatment_context", {})
        out += [
            "<h3>Clusters (descriptive)</h3>",
            f"<p>At a TN93 distance of {cl['threshold']}, {cl['n_clustered']} of "
            f"{cl['n_sequences']} sequences ({_pct(cl['proportion_clustered'])}) fall in "
            f"{cl['n_clusters']} clusters, the largest of {cl['largest_cluster']}. These counts "
            "come from the tree set only and are not transmission chains.</p>",
            _table(
                ["", "n", "With SDRM", "Proportion"],
                [
                    [
                        "In a cluster",
                        s["clustered"]["n"],
                        s["clustered"]["with_sdrm"],
                        _pct(s["clustered"]["proportion"]),
                    ],
                    [
                        "Not in a cluster",
                        s["unclustered"]["n"],
                        s["unclustered"]["with_sdrm"],
                        _pct(s["unclustered"]["proportion"]),
                    ],
                ],
            ),
            f"<p>Adjusted for study treatment context (Mantel-Haenszel), the odds ratio for an "
            f"SDRM being in a cluster is {adj.get('odds_ratio')} (95% interval "
            f"{adj.get('ci_low')} to {adj.get('ci_high')}). The crude difference reflects which "
            "cohorts are densely sampled.</p>",
        ]
    return "\n".join(out)


def render(settings: Settings, root: Path) -> str:
    results = root / settings.paths["results"] / settings.scope
    for path, stage in (
        (results / "drm" / "prevalence.tsv", "drm prevalence"),
        (results / "drm" / "trend_tests.tsv", "drm prevalence"),
    ):
        if not path.exists():
            raise ContractError(f"{path} not found; run `virotrace {stage}` first")
    items = "".join(f"<li>{escape(t)}</li>" for t in limitations(root, settings.scope))
    scope = escape(settings.scope.capitalize())
    stamp = datetime.now(UTC).strftime("%Y-%m-%d")
    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>Virotrace: {scope} HIV-1 drug resistance surveillance report</title>
<style>
body{{font:16px/1.5 system-ui,sans-serif;max-width:60rem;margin:2rem auto;padding:0 1rem;color:#222}}
table{{border-collapse:collapse;margin:.5rem 0 1rem;font-size:.9rem}}
th,td{{border:1px solid #ccc;padding:.25rem .5rem;text-align:left;vertical-align:top}}
th{{background:#f3f3f3}} .notice{{background:#fff4e5;border-left:4px solid #e08a00;padding:.5rem 1rem}}
</style></head><body>
<h1>{scope} HIV-1 pol: drug resistance and lineage report</h1>
<p>Generated {stamp} from the {scope} build.</p>
<p class="notice"><strong>Not a clinical tool.</strong> These are descriptive statistics from
publicly deposited sequences. No result should inform the care of an individual patient.</p>
<h2>Summary of what this describes</h2>
<p>Every prevalence figure describes deposited sequences, one per individual, not the
epidemic, and shows its denominator. Subtyping and recombination are screens. Read the
limitations before quoting a number.</p>
<h2>Data and quality control</h2>
{_data_section(results)}
<h2>Drug resistance prevalence</h2>
{_prevalence_sections(results)}
<h2>Lineages and clusters</h2>
{_phylo_section(results, settings.phylo.refine.subtype_trees)}
<h2>Limitations</h2>
<ol>{items}</ol>
<h2>Provenance</h2>
<p>Scope <code>{scope.lower()}</code>. The fully resolved configuration is in
<code>results/{escape(settings.scope)}/resolved_config.yaml</code> and snapshot manifests list the
exact queries, so the corpus can be rebuilt. SDRM lists: WHO 2009 (headline), Stanford
2026-09-25 (sensitivity).</p>
</body></html>
"""


def run_report(settings: Settings, *, repo_root: Path | str = ".") -> dict[str, Any]:
    root = Path(repo_root)
    html = render(settings, root)
    out = root / settings.paths["reports"] / settings.scope / "index.html"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html)
    return {"scope": settings.scope, "report": str(out), "bytes": len(html.encode())}
