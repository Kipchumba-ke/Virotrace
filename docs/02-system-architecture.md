# Virotrace — System Architecture

---

## 1. Shape of the system

Virotrace is a **batch DAG, not a service**. There is no API, no database server, and no long-running
process in v1. The whole system is a directed acyclic graph of file transformations orchestrated by
Snakemake, plus a static viewer (Auspice) and a static report.

This is the right shape because:
- the question is retrospective, not real-time;
- the compute budget is one laptop;
- reproducibility is the headline feature, and a DAG over immutable files is the cheapest way to get it;
- it degrades gracefully — any stage can be run, inspected, and cached independently.

Adding a web service in v1 would add deployment surface and subtract nothing from the credibility of
the result. Noted as a v2 option in §7.

## 2. Layers

```
┌──────────────────────────────────────────────────────────────────┐
│ L5  PRESENTATION      Auspice (tree) · HTML report · CSV exports │
├──────────────────────────────────────────────────────────────────┤
│ L4  EXPORT            Auspice v2 JSON assembly · report render   │
├──────────────────────────────────────────────────────────────────┤
│ L3  ANALYSIS                                                     │
│     ┌───────────────────────┐   ┌──────────────────────────────┐ │
│     │ L3a DRM layer         │   │ L3b PHYLO layer              │ │
│     │  mutation calling     │   │  tree · time-tree · traits   │ │
│     │  prevalence + CIs     │   │  clusters                    │ │
│     └───────────────────────┘   └──────────────────────────────┘ │
├──────────────────────────────────────────────────────────────────┤
│ L2  CURATION          QC · align · subtype · recombination       │
│                       subsample · exclusion ledger               │
├──────────────────────────────────────────────────────────────────┤
│ L1  NORMALISATION     metadata harmonisation · schema validation │
├──────────────────────────────────────────────────────────────────┤
│ L0  ACQUISITION       source adapters · immutable snapshot       │
└──────────────────────────────────────────────────────────────────┘
```

The **L1/L2 boundary is the load-bearing interface**: one validated metadata table plus one curated
FASTA. Both analysis layers consume exactly that and nothing upstream of it. This is what lets the
DRM layer ship first without the phylo layer being a rewrite.

## 3. The DAG

```
                 config.yaml  +  sources.yaml
                          │
                 ┌────────▼────────┐
                 │ fetch (per src) │   L0  → data/raw/<snapshot>/
                 └────────┬────────┘        + manifest.json (sha256, query, date)
                          │
                 ┌────────▼────────┐
                 │ normalise       │   L1  → metadata.raw.tsv
                 └────────┬────────┘
                 ┌────────▼────────┐
                 │ validate schema │   L1  → metadata.valid.tsv
                 └────────┬────────┘        quarantine.tsv
                          │
                 ┌────────▼────────┐
                 │ qc sequences    │   L2  → seqs.qc.fasta
                 └────────┬────────┘        excluded.tsv  ◄── ledger
                          │
                 ┌────────▼────────┐
                 │ align to HXB2   │   L2  → aln.pol.fasta (codon-aware)
                 └────────┬────────┘
                   ┌──────┴───────┐
         ┌─────────▼──────┐   ┌───▼──────────────┐
         │ subtype call   │   │ recombination    │   L2
         └─────────┬──────┘   └───┬──────────────┘
                   └──────┬───────┘
                 ┌────────▼────────┐
                 │ subsample       │   L2  → aln.final.fasta  (stratified)
                 └────────┬────────┘        metadata.final.tsv
                          │
            ┌─────────────┴──────────────┐
   ┌────────▼─────────┐        ┌─────────▼──────────┐
   │ L3a  call DRMs   │        │ L3b  tree (IQ-TREE)│
   │      ↓           │        │      ↓             │
   │  prevalence +CI  │        │  refine (TreeTime) │
   │      ↓           │        │      ↓             │
   │  drm_table.tsv   │        │  ancestral/translate
   │  prevalence.tsv  │        │      ↓             │
   └────────┬─────────┘        │  traits (mugration)│
            │                  │  clusters (TN93)   │
            │                  └─────────┬──────────┘
            │          ┌─────────────────┘
            │          │
       ┌────▼──────────▼────┐
       │ L4 export auspice  │  → auspice/virotrace_<scope>_pol.json
       │ L4 render report   │  → reports/<scope>/index.html
       └────────────────────┘
```

DRM output feeds the Auspice export as node colourings — that is the join point between the two
layers, and it is a file contract (`drm_table.tsv` keyed by strain), not a code dependency.

## 4. Components and technology choices

| Component | Choice | Why this, here |
|---|---|---|
| Orchestration | **Snakemake** | File-based DAG, resumable, config-driven, `--dry-run` as a free architecture diagram, conda integration per-rule. Nextflow is equally valid but heavier for a solo laptop project; Make lacks wildcards and provenance |
| Language | **Python 3.11** | Matches existing skills; Biopython, pandas, pydantic ecosystem |
| Metadata contract | **pandera** schema + **pydantic** config models | Validation at the boundary, with failure cases as data not exceptions |
| Alignment | **MAFFT** via `augur align --reference-sequence HXB2` | Reference-anchored so coordinates stay interpretable; codon-aware handling downstream |
| Subtyping | reference-set phylogenetic placement, cross-checked against an external subtyping tool | Avoids a single opaque call; disagreements get flagged rather than averaged |
| Recombination screen | **3SEQ** (or RDP-class equivalent) | Needed because HIV recombination breaks the single-tree assumption |
| Resistance algorithm | local implementation of the **Stanford HIVdb** algorithm, with the public Sierra service as an optional cached cross-check | Offline-first keeps the pipeline deterministic and runnable without network; the cross-check catches algorithm-version drift |
| Tree | **IQ-TREE 2**, GTR+F+R4, UFBoot 1000 (`-fast` at africa scope) | Best accuracy-per-second on a laptop for this alignment size |
| Time tree + traits | **TreeTime** via `augur refine` / `augur traits` | Standard, scriptable, and what Auspice expects |
| Clustering | TN93 pairwise distance with a configured threshold | Transparent and tunable; threshold is a config value, not a magic number in code |
| Viewer | **Auspice v2** | The ask; also means the output format is a published spec rather than bespoke |
| Report | Jinja2 → HTML, figures via matplotlib | No JS build step; one artifact, openable offline |
| Env | conda/mamba + lockfile; optional Docker | Bioinformatics tooling is not pip-installable in practice |

## 5. Data contracts

Every arrow in §3 is a typed file. The three that matter:

**`metadata.valid.tsv`** — one row per sequence. Columns and types are fixed in
`config/schema/metadata.py`; see `03-data-acquisition.md` §4 for the field list. Breaking this schema
is a breaking change and bumps the pipeline minor version.

**`aln.final.fasta`** — reference-anchored *pol* alignment. Headers are `strain` only (no embedded
metadata; that is what the TSV is for). HXB2 retained as the first record for coordinate sanity, and
dropped at export.

**`drm_table.tsv`** — `strain · gene · mutation · drug_class · is_sdrm · hivdb_score · interpretation ·
algorithm_version`. Long format, one mutation per row. The phylo layer joins on `strain` and never
needs to know how the call was made.

Provenance block: every terminal output carries a sidecar `*.provenance.json` with config hash,
snapshot ID, tool versions, seeds, and UTC timestamp.

## 6. Directory layout at runtime

```
data/
  raw/<snapshot_id>/        immutable, never edited, manifest.json alongside
  interim/<scope>/          QC'd, aligned, subtyped — regenerable
  processed/<scope>/        final alignment + metadata, the analysis inputs
  reference/                HXB2, subtype reference set, DRM lists (pinned versions)
results/<scope>/
  drm/                      drm_table.tsv, prevalence.tsv, figures/
  phylo/                    tree.nwk, branch_lengths.json, traits.json, clusters.tsv
  qc/                       excluded.tsv, quarantine.tsv, qc_summary.html
auspice/                    *.json
reports/<scope>/            index.html + assets
logs/<scope>/               per-rule logs
```

`data/raw` being immutable is enforced by convention plus a CI check that no rule declares it as an
output. Snapshots are identified by date + config hash and are the unit of "which data version was
this result built from".

## 7. Named scale-up points

The architecture is laptop-shaped but not laptop-locked. Each of these is a profile or config change,
not a redesign:

| Pressure | Change | Not a redesign because |
|---|---|---|
| More sequences than subsampling caps allow | Raise `subsample.max_per_stratum`; switch IQ-TREE off `-fast` | Caps are config values consumed by one rule |
| Tree too slow | Move the `tree` and `refine` rules to a cluster executor | Snakemake rule-level resources already declared |
| Want Bayesian dating | Add a `beast` rule consuming `aln.final.fasta` | The L2/L3 contract already provides exactly its input |
| Want an interactive query layer | Load `drm_table.tsv` + `prevalence.tsv` into DuckDB behind a thin API | DRM outputs are already long-format tabular |
| Controlled-access data becomes available | New L0 adapter + an access tier flag in metadata | Source adapters are pluggable; nothing downstream knows the source |

## 8. Failure posture

- **Fail loud at the boundary, degrade gracefully inside.** Schema violations halt the run with a
  readable diff. Individual bad sequences are quarantined, counted, and the run continues.
- **No silent fallbacks.** If the resistance cross-check service is unreachable, the run proceeds on
  the local algorithm and writes that fact into provenance — it does not pretend both agreed.
- **Gates, not hopes.** The temporal-signal gate (root-to-tip R² below threshold) stops the time-tree
  and emits the divergence tree with an explanation. A pipeline that always produces a confident
  answer is a pipeline that will eventually produce a wrong one.
