# Virotrace — Technical Specification

**Version:** 0.1 (pre-implementation)
**Owner:** Kipchumba
**Start date:** 2026-11-04
**Status:** Draft, ready to build against

---

## 1. Purpose

Virotrace is a reproducible, configuration-driven pipeline and visualisation layer for
**HIV-1 molecular surveillance in Kenya, situated in a wider African context**. It answers two
classes of question from the same curated sequence corpus:

- **Drug resistance (DRM layer).** Which resistance mutations are present in publicly deposited
  Kenyan HIV-1 *pol* sequences, at what prevalence, in which drug classes, and how has that
  changed over time relative to the rest of East and Southern Africa?
- **Lineage phylodynamics (phylo layer).** What subtypes and circulating recombinant forms are
  represented, how are they structured phylogenetically, and what does a time-resolved tree with
  geographic trait reconstruction suggest about lineage movement?

The output is a public Nextstrain/Auspice build plus a tabular resistance report, both regenerable
end-to-end from a single command and a config file.

## 2. Scope

### In scope (v1)
- Public-tier sequence acquisition from LANL HIV Sequence Database and NCBI Virus/GenBank.
- *pol* gene region (protease, reverse transcriptase, integrase) as the primary analytical target.
- Genotypic resistance interpretation against the Stanford HIVdb algorithm, plus the WHO
  surveillance drug resistance mutation (SDRM) list for surveillance-grade prevalence.
- Subtype assignment and recombination screening.
- Maximum-likelihood tree, time-calibrated refinement, ancestral trait reconstruction over
  country/region, exported as Auspice v2 JSON.
- Static resistance report (tables + figures) and a browsable tree.

### Out of scope (v1) — recorded so it stays out
- **Controlled-access data.** National surveillance repositories, PANGEA, and any dataset requiring
  an institutional data-transfer agreement. v1 is designed to be scientifically defensible on
  public data alone. See `03-data-acquisition.md` §2 for the access-tier policy.
- **Clinical outcome linkage.** No treatment history, viral load, or CD4 data — public deposits
  almost never carry it, and inferring treatment status from mutations is an overreach.
- **Transmission inference as a claim about individuals.** Cluster detection is included as a
  descriptive population-level tool only. Virotrace will not assert who infected whom; see
  §7 Ethical constraints.
- **Whole-genome HIV alignment.** Recombination and indel density make genome-wide HIV alignment a
  research problem in its own right. v1 stays in *pol*.
- Bayesian phylodynamics (BEAST/BEAST2). Documented as a scale-up path, not a v1 deliverable.

### Non-goals
Virotrace is not a clinical decision-support tool. It produces no patient-level resistance
interpretation for care, and every output carries that disclaimer.

## 3. Users

| User | Need | What they touch |
|---|---|---|
| Public health analyst (primary) | Regional DRM prevalence and trend, by class and year | Auspice view, resistance report, CSV exports |
| Bioinformatician / reviewer | Verify the method, rerun with different filters | config.yaml, Snakemake DAG, logs |
| Portfolio reviewer / hiring manager | Judge engineering and scientific judgement in ~10 minutes | README, Auspice link, docs/, test suite, CI badge |

The third user is a real constraint, not a vanity one. It means the repo must be legible and the
build must run on a clean machine without a support call.

## 4. Functional requirements

| ID | Requirement | Acceptance |
|---|---|---|
| FR-1 | Acquire sequences and metadata from configured public sources into an immutable raw snapshot | Snapshot directory with manifest, per-file SHA256, and source query recorded |
| FR-2 | Normalise heterogeneous metadata into one validated schema | Every record passes schema validation or is quarantined with a reason code |
| FR-3 | QC-filter sequences: length, ambiguity fraction, frameshift, internal stop, APOBEC hypermutation | Each exclusion written to `excluded.tsv` with rule and value |
| FR-4 | Assign subtype and flag probable recombinants | Per-sequence subtype call with method and confidence recorded |
| FR-5 | Call DRMs against HIVdb and tag SDRMs | Per-sequence mutation list by drug class; per-drug penalty score and interpretation level |
| FR-6 | Estimate DRM prevalence stratified by scope, year, and drug class with confidence intervals | Wilson CIs; denominators explicit; n<threshold suppressed not silently plotted |
| FR-7 | Build an ML tree and a time-resolved tree over the filtered alignment | Newick + branch-length-annotated node data; bootstrap support retained |
| FR-8 | Reconstruct geographic traits and export Auspice v2 JSON with DRM colourings | Auspice JSON loads in `auspice view`; DRM status selectable as a colour-by |
| FR-9 | Produce a self-contained resistance report | HTML/markdown report with tables, figures, methods, caveats, data version |
| FR-10 | Reproduce any build from a config file and a data snapshot | Rerunning a pinned config regenerates byte-identical outputs for deterministic stages |
| FR-11 | Run two scopes — `kenya` and `africa` — from the same codebase via config | `snakemake --config scope=africa` produces a parallel build without code edits |

## 5. Non-functional requirements

**Compute envelope (binding constraint: single laptop).**

| Budget | Target | Enforcement |
|---|---|---|
| Peak RAM | ≤ 6 GB | Subsampling caps in config; alignment chunking |
| Full `dev` profile wall clock | ≤ 10 min | Toy dataset, ~200 seqs |
| Full `kenya` profile wall clock | ≤ 90 min | Subsample cap ~1,500 seqs, IQ-TREE with `-fast` fallback |
| Full `africa` profile wall clock | ≤ 6 h, resumable | Snakemake checkpointing; documented overnight run |
| Disk | ≤ 20 GB including snapshots | Snapshot pruning policy |

Tree inference is the cost centre. The design response is sequence subsampling **by stratum**
(country × year × subtype) rather than random downsampling, so that reducing n degrades precision
without biasing composition. Scale-up points are named explicitly in `02-system-architecture.md` §7
so the same pipeline can move to a server by changing a profile, not a method.

**Other non-functionals**
- *Reproducibility.* Pinned environment (conda lockfile), pinned reference data versions, recorded
  random seeds, recorded tool versions in every output's provenance block.
- *Transparency of exclusion.* No sequence disappears silently. Every filtered record is
  accounted for in a machine-readable exclusion ledger — this is the single most important
  credibility feature of a surveillance pipeline.
- *Determinism.* Stages that can be deterministic must be: fixed seeds, sorted inputs, stable
  tie-breaking. Stochastic stages (bootstrap, MCMC-free heuristics) record their seed.
- *Portability.* Linux and macOS; conda environment plus an optional container.
- *Legibility.* Any stage runnable standalone from the CLI on its declared inputs.

## 6. Deliverables

1. Public GitHub repository with CI green.
2. Auspice build(s) — `virotrace/kenya/pol` and `virotrace/africa/pol`.
3. Resistance surveillance report (HTML + CSV tables).
4. `docs/` — this spec and the six companion documents.
5. A written limitations section that a reviewer will find before they find a flaw.

## 7. Ethical and interpretive constraints

These are requirements, not boilerplate.

- **Sampling bias is the dominant limitation.** Public HIV sequences from Kenya are a convenience
  sample: clustered by study, site, and era. Prevalence estimates describe *deposited sequences*,
  not the Kenyan epidemic. Every prevalence figure must carry its denominator and this framing.
- **No individual-level inference.** Cluster membership is reported as a descriptive statistic.
  The repository will not frame any result as evidence about transmission between identifiable
  people, and the report says so.
- **Geographic resolution is coarse and must not be overstated.** Most GenBank/LANL records carry
  country, not subnational location. Sub-national mapping will be attempted only where a record
  explicitly supplies it, and the map legend must state the proportion of records resolved only to
  country level. Inventing county assignments to make a prettier map is the failure mode to avoid.
- **Non-B subtype caveat.** The HIVdb algorithm is best validated on subtype B. Kenya's epidemic is
  predominantly A1 with D, C, and recombinant forms. Interpretation levels are reported with that
  caveat attached, and polymorphism-vs-resistance ambiguity is flagged rather than smoothed over.
- **No clinical use.** Prominent disclaimer in README, report, and Auspice build description.

## 8. Milestones

| Phase | Window | Exit criterion |
|---|---|---|
| **P0 — Foundation** | Day 1–3 | Repo scaffold, env builds, CI green, toy dataset runs end-to-end through a stub DAG |
| **P1 — Data layer** | Day 4–9 | FR-1 to FR-3 complete; validated metadata table; exclusion ledger; snapshot manifest |
| **P2 — DRM layer** | Day 10–18 | FR-4 to FR-6 complete; resistance report on real Kenyan data; positive-control tests passing |
| **P3 — Phylo layer** | Day 19–30 | FR-7, FR-8 complete; Auspice build live for `kenya` |
| **P4 — Africa scope + polish** | Day 31–40 | FR-11; `africa` build; report finalised; limitations written; README demoable |

P2 before P3 is deliberate: it produces a defensible result early, and the metadata schema it forces
you to settle is exactly what the phylo layer needs.

## 9. Risks

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Too few Kenyan *pol* sequences with usable collection dates | Medium | High — kills the time-tree | Quantify in P1 **before** building P3. Fallback: pool East Africa for the tree, keep Kenya as a trait |
| Recombinant sequences violate single-tree assumptions | High | High | Screen and either exclude or partition; report the excluded fraction; never silently leave them in |
| Programmatic access to a resistance-calling service is rate-limited or offline | Medium | Medium | Use a local algorithm implementation as the default path; treat any web service as an optional cross-check, cached |
| IQ-TREE runtime exceeds laptop patience at `africa` scope | Medium | Medium | Stratified subsampling caps; `-fast` mode; documented overnight path |
| Weak temporal signal → meaningless time-tree | Medium | High | Root-to-tip regression as a gate. If R² fails the configured threshold, publish the divergence tree and say why, rather than a time-tree you don't believe |
| Scope creep into two half-built projects | High | High | Phase gates above; `africa` and phylo layers only start when their predecessor's exit criteria are met |

## 10. Success criteria

Virotrace v1 succeeds if a competent stranger can, in one sitting:
clone it, build the environment, run the `dev` profile to completion, read the methodology, find the
limitations without being told where to look, and reproduce a published figure from the pinned
config. Scientific novelty is explicitly **not** a success criterion; defensibility and
reproducibility are.
