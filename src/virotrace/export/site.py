"""The explorer: one self-contained HTML page over the aggregate results.

Like the report, the site computes nothing new. It reads the tables and summaries the stages
wrote, embeds them as JSON in a single page, and lets a reader filter the resistance estimates
by list, stratum and drug class. Nothing leaves the page: no server, no external script, font or
tile. Only aggregate tables are embedded; accessions, strains and per-sequence rows are not.
The limitations are read from the same place as the report's, so the two cannot drift apart.
"""

from __future__ import annotations

import json
import re
import shutil
import unicodedata
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from ..curate.study import load_study_annotations
from ..errors import ContractError
from ..io.tables import read_tsv
from ..settings import Settings
from .africa_map import AFRICA
from .report import REFERENCES, _read_json, limitations
from .site_template import PAGE

# GenBank /country names that differ from the Natural Earth names on the outline.
COUNTRY_ALIASES = {
    "democratic republic of the congo": "COD", "democratic republic of congo": "COD",
    "dr congo": "COD", "drc": "COD", "zaire": "COD",
    "republic of the congo": "COG", "republic of congo": "COG", "congo": "COG",
    "cote d ivoire": "CIV", "ivory coast": "CIV",
    "tanzania": "TZA", "united republic of tanzania": "TZA",
    "swaziland": "SWZ", "eswatini": "SWZ",
    "central african republic": "CAF", "equatorial guinea": "GNQ",
    "south sudan": "SSD", "western sahara": "ESH",
    "sao tome and principe": "STP", "the gambia": "GMB", "gambia": "GMB",
}  # fmt: skip


def _num(value: Any, kind: type = float) -> Any:
    if value in ("", None):
        return None
    try:
        return kind(float(value)) if kind is int else float(value)
    except (TypeError, ValueError):
        return None


def _prevalence(results: Path) -> list[dict[str, Any]]:
    out = []
    for r in read_tsv(results / "drm" / "prevalence.tsv"):
        out.append(
            {
                "list_role": r["list_role"],
                "list_version": r["list_version"],
                "stratum_kind": r["stratum_kind"],
                "stratum": r["stratum"],
                "drug_class": r["drug_class"],
                "numerator": _num(r.get("numerator"), int),
                "denominator": _num(r.get("denominator"), int),
                "proportion": _num(r.get("proportion")),
                "ci_low": _num(r.get("ci_low")),
                "ci_high": _num(r.get("ci_high")),
                "status": r["status"],
            }
        )
    return out


def _trends(results: Path) -> list[dict[str, Any]]:
    return [
        {
            "drug_class": t["drug_class"],
            "n_bins": _num(t.get("n_bins"), int),
            "first_bin": t.get("first_bin", ""),
            "last_bin": t.get("last_bin", ""),
            "z": _num(t.get("z")),
            "p_two_sided": _num(t.get("p_two_sided")),
            "note": t.get("note", ""),
        }
        for t in read_tsv(results / "drm" / "trend_tests.tsv")
        if t["list_role"] == "headline"
    ]


def _gate_row(label: str, gate: dict[str, Any], time_tree: bool) -> dict[str, Any]:
    return {
        "label": label,
        "n": gate.get("n"),
        "slope": gate.get("slope"),
        "r_squared": gate.get("r_squared"),
        "permutation_p": gate.get("permutation_p"),
        "passed": bool(gate.get("passed")),
        "result": "time tree" if time_tree else "divergence tree",
    }


def _gates(results: Path, subtypes: list[str]) -> list[dict[str, Any]]:
    gate = _read_json(results / "phylo" / "temporal_gate.json")
    if gate is None:
        return []
    rows = [_gate_row("Pooled (all subtypes)", gate["pooled"], gate.get("time_tree", False))]
    for label, g in sorted(gate.get("by_subtype", {}).items()):
        rows.append(_gate_row(f"{label} within the pooled tree", g, False))
    for s in subtypes:
        sub = _read_json(results / "phylo" / s / "temporal_gate.json")
        if sub:
            rows.append(_gate_row(f"{s} tree", sub["pooled"], sub.get("time_tree", False)))
    return rows


def _norm(name: str) -> str:
    text = unicodedata.normalize("NFKD", name).encode("ascii", "ignore").decode().lower()
    return re.sub(r"\s+", " ", re.sub(r"[^a-z0-9]+", " ", text.replace("&", " and "))).strip()


def country_ids(strata: list[str]) -> tuple[dict[str, str], list[str]]:
    """Map GenBank country names to outline ids; names with no outline are returned apart."""
    lookup = {_norm(n): f["id"] for f in AFRICA["features"] for n in f["names"]}
    lookup.update(COUNTRY_ALIASES)
    found: dict[str, str] = {}
    missing: list[str] = []
    for name in strata:
        key = lookup.get(_norm(name))
        if key:
            found[name] = key
        else:
            missing.append(name)
    return found, missing


def _annotation_coverage(settings: Settings, root: Path) -> dict[str, Any] | None:
    """Share of QC-passed sequences whose study has a hand-curated treatment context."""
    interim = root / settings.paths["interim"] / settings.scope
    passed_path, meta_path = interim / "qc.passed.tsv", interim / "metadata.valid.tsv"
    if not (passed_path.exists() and meta_path.exists()):
        return None
    passed = {r["strain"] for r in read_tsv(passed_path)}
    ref = root / settings.paths["reference"] / settings.prevalence.study_annotations
    annotations = load_study_annotations(ref) if ref.exists() else {}
    by_context: Counter[str] = Counter()
    studies: set[str] = set()
    annotated_studies: set[str] = set()
    for r in read_tsv(meta_path):
        if r["strain"] not in passed:
            continue
        key = r.get("study_key", "") or "unknown_study"
        studies.add(key)
        context = annotations.get(key, {}).get("treatment_context")
        if context:
            annotated_studies.add(key)
        by_context[context or "unannotated"] += 1
    total = sum(by_context.values())
    if not total:
        return None
    return {
        "n_sequences": total,
        "n_annotated": total - by_context["unannotated"],
        "by_context": dict(by_context),
        "n_studies": len(studies),
        "n_studies_annotated": len(annotated_studies),
    }


def _copy_trees(settings: Settings, root: Path, out: Path) -> list[dict[str, Any]]:
    """Copy the Auspice JSONs beside the page so they can be downloaded or opened elsewhere."""
    src = root / settings.paths["auspice"]
    files = sorted(src.glob(f"virotrace_{settings.scope}*.json")) if src.exists() else []
    if not files:
        return []
    target = out / "trees"
    target.mkdir(parents=True, exist_ok=True)
    listed = []
    for f in files:
        shutil.copy2(f, target / f.name)
        label = f.stem.removeprefix(f"virotrace_{settings.scope}").lstrip("_")
        listed.append(
            {
                "name": f.name,
                "label": label or "pooled",
                "file": f"trees/{f.name}",
                "bytes": f.stat().st_size,
            }
        )
    return listed


def collect(
    settings: Settings, root: Path, *, trees: list[dict[str, Any]] | None = None
) -> dict[str, Any]:
    results = root / settings.paths["results"] / settings.scope
    for path, stage in (
        (results / "drm" / "prevalence.tsv", "drm prevalence"),
        (results / "drm" / "trend_tests.tsv", "drm prevalence"),
    ):
        if not path.exists():
            raise ContractError(f"{path} not found; run `virotrace {stage}` first")
    prevalence = _prevalence(results)
    versions = {r["list_role"]: r["list_version"] for r in prevalence}
    countries = sorted({r["stratum"] for r in prevalence if r["stratum_kind"] == "country"})
    ids, missing = country_ids(countries)
    cl = _read_json(results / "phylo" / "clusters_summary.json")
    st = _read_json(results / "phylo" / "subtype_tree" / "subtype_tree_summary.json")
    return {
        "scope": settings.scope,
        "generated": datetime.now(UTC).strftime("%Y-%m-%d"),
        "min_stratum_n": settings.report.min_stratum_n,
        "list_versions": versions,
        "prevalence": prevalence,
        "trends": _trends(results),
        "map": (
            {
                "width": AFRICA["width"],
                "height": AFRICA["height"],
                "features": [
                    {k: f[k] for k in ("id", "name", "d", "cx", "cy")} for f in AFRICA["features"]
                ],
            }
            if ids
            else None
        ),
        "country_ids": ids,
        "unmapped_countries": missing,
        "annotation": _annotation_coverage(settings, root),
        "qc": _read_json(results / "qc" / "qc_summary.json"),
        "tree_set": _read_json(results / "qc" / "subsample_summary.json"),
        "clusters": cl,
        "gates": _gates(results, settings.phylo.refine.subtype_trees),
        "subtype_tree": (
            {
                "n_queries": st.get("n_queries"),
                "n_references": st.get("n_references"),
                "support_used": bool(st.get("support_used")),
                "concordance": st.get("concordance", {}),
            }
            if st
            else None
        ),
        "trees": trees or [],
        "limitations": limitations(root, settings.scope),
        "references": REFERENCES,
    }


def render(data: dict[str, Any]) -> str:
    # No "<" can appear in the data block, so no "</script>" or "<!--" inside it can end or
    # confuse it; \u003c is the same character to the JSON parser.
    blob = json.dumps(data, ensure_ascii=False, separators=(",", ":")).replace("<", "\\u003c")
    scope = data["scope"].capitalize()
    return PAGE.replace("__TITLE__", f"Virotrace: {scope} HIV-1 resistance explorer").replace(
        "__DATA__", blob
    )


def run_site(
    settings: Settings, *, repo_root: Path | str = ".", base_url: str | None = None
) -> dict[str, Any]:
    root = Path(repo_root)
    out = root / settings.paths["reports"] / settings.scope / "site"
    out.mkdir(parents=True, exist_ok=True)
    trees = _copy_trees(settings, root, out)
    if base_url:
        # nextstrain.org fetches a public JSON by address; it needs the page to be published.
        host = base_url.removeprefix("https://").removeprefix("http://").rstrip("/")
        for t in trees:
            t["nextstrain"] = f"https://nextstrain.org/fetch/{host}/{t['file']}"
    data = collect(settings, root, trees=trees)
    (out / "index.html").write_text(render(data), encoding="utf-8")
    (out / "data.json").write_text(json.dumps(data, indent=1, ensure_ascii=False), "utf-8")
    return {
        "scope": settings.scope,
        "site": str(out / "index.html"),
        "data": str(out / "data.json"),
        "n_trees": len(trees),
        "unmapped_countries": data["unmapped_countries"],
        "annotated_fraction": (
            round(data["annotation"]["n_annotated"] / data["annotation"]["n_sequences"], 4)
            if data["annotation"]
            else None
        ),
        "bytes": (out / "index.html").stat().st_size,
    }
