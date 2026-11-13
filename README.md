# Virotrace

Reproducible HIV-1 molecular surveillance for **Kenya**, in an **Africa-wide** context:
genotypic drug-resistance surveillance and lineage phylogenomics over the same curated
corpus of publicly deposited *pol* sequences, with a Nextstrain/Auspice build and a
resistance report as outputs.

> **Not a clinical tool.** Virotrace produces population-level descriptive statistics from
> public sequence data. No output should be used to inform the care of any individual
> patient.

---

## What it does

| Layer | Question | Output |
|---|---|---|
| **DRM** | Which resistance mutations are present in deposited Kenyan HIV-1 *pol* sequences, at what prevalence, in which drug classes, and how does that compare with the wider region? | `results/<scope>/drm/prevalence.tsv`, resistance report |
| **Phylo** | What subtypes and recombinant forms are represented, how are they structured, and what does a time-resolved tree with geographic traits show? | `auspice/virotrace_<scope>.json`, plus `virotrace_<scope>_<subtype>.json` dated trees |

The two layers meet at a file contract (`drm_table.tsv`), not a code dependency, so
resistance mutations appear as colourings on the tree without either layer importing the
other.

## Quickstart

```bash
git clone https://github.com/Kipchumba-ke/Virotrace && cd virotrace
make setup && conda activate virotrace

virotrace info                     # resolved versions of every external binary
make dag                           # validate the DAG and all three profiles, run nothing
make test                          # fast suite, under a minute
make dev                           # toy build end to end, minutes, no network
```

Then the real builds:

```bash
make kenya                         # primary build
make africa CORES=4                # context build; long, resumable
```

## Scope and honesty about the data

Public HIV sequences from Kenya are a **convenience sample**, clustered by study, site and
era. Every prevalence figure describes *deposited sequences*, not the Kenyan epidemic, and
carries its denominator. Three further limits are built into the design rather than noted
afterwards:

- **Geographic resolution is country-level for most records.** Subnational reconstruction
  runs only when enough records genuinely supply a location, and the map legend states the
  resolved fraction. No invented county assignments.
- **Treatment status is unknown for most records**, and is never inferred from the presence
  of mutations — that would be circular. So pretreatment and acquired resistance cannot be
  cleanly separated in v1.
- **Recombinants leave the tree but stay in the resistance layer.** Recombination
  invalidates a phylogeny, not a genotype. The excluded fraction is reported in both outputs.

The full list is in [`docs/04-bioinformatics-methodology.md`](docs/04-bioinformatics-methodology.md#9-stated-limitations-to-be-carried-verbatim-into-the-report),
and the report cannot be generated without it.

## Documentation

| Doc | Contents |
|---|---|
| [01 Technical specification](docs/01-technical-specification.md) | Scope, requirements, compute envelope, milestones, risks |
| [02 System architecture](docs/02-system-architecture.md) | Layers, DAG, data contracts, technology rationale, scale-up points |
| [03 Data acquisition](docs/03-data-acquisition.md) | Sources, access tiers, metadata schema, cohort and subsampling design |
| [04 Bioinformatics methodology](docs/04-bioinformatics-methodology.md) | QC, alignment, subtyping, recombination, DRM calling, prevalence, phylogenetics |
| [05 Software architecture](docs/05-software-architecture.md) | Repository structure, module boundaries, CLI, error handling |
| [06 Configuration](docs/06-configuration.md) | The config file, profiles, validation, what is deliberately not configurable |
| [07 Testing strategy](docs/07-testing-strategy.md) | Positive and negative controls, invariants, golden files, CI |

## Project status

The Kenya build runs end to end on real data. Each stage is a `virotrace` command, and the
Snakemake rules in `workflow/` run the same stages in the same order:

```bash
virotrace acquire   --profile kenya      # L0 snapshot with manifest
virotrace normalise --profile kenya      # L1 dedup, scope filter, validation
virotrace align     --profile kenya      # L2 MAFFT to HXB2 pol coordinates
virotrace qc        --profile kenya      # L2 QC ledger; clone and individual policies
virotrace drm call  --profile kenya      # L3a mutations at SDRM positions
virotrace subtype   --profile kenya      # L2 subtype placement (depositor label compared)
virotrace recombination --profile kenya  # L2 breakpoint screen
virotrace drm prevalence --profile kenya # L3a per-class prevalence, Wilson CIs, trends
virotrace mask      --profile kenya      # L2 gap-masked alignment for the tree only
virotrace subsample --profile kenya      # L2 tree set: strictly dated, one per individual
virotrace phylo tree     --profile kenya # L3b IQ-TREE
virotrace phylo refine   --profile kenya # L3b temporal gate, then time tree or divergence tree
virotrace phylo subtrees --profile kenya # L3b per-subtype trees and gates (Kenya: A1, D)
virotrace phylo clusters --profile kenya # L3b TN93 clusters (descriptive only)
virotrace phylo clock-sensitivity --profile kenya  # L3b root dates under fixed clock rates
virotrace phylo subtype-tree --profile kenya       # optional: tree-based second subtype call
virotrace export auspice --profile kenya # L4 Auspice datasets in auspice/
virotrace export report  --profile kenya # L4 HTML report in reports/<scope>/index.html
auspice view --datasetDir auspice
```

Implemented: acquisition (GenBank), normalisation, QC, alignment, DRM calling and prevalence,
person identification across studies, subtype placement (distance-based, plus an optional
tree-based second call), a recombination screen, masking, subsampling, tree inference, the
temporal-signal gate with date randomisation, per-subtype dated trees, clock-rate sensitivity,
clusters, the Auspice export, and the HTML report with its limitations and references.
Both the Kenya and the Africa builds have run end to end on real data.
Not yet: ancestral trait reconstruction (`phylo traits`, off in the DAG until it exists), the
LANL adapter, and amino-acid mutation tracks in the export.

### The Africa build

```bash
# VIROTRACE_NCBI_EMAIL (required) and VIROTRACE_NCBI_API_KEY go in .env
virotrace acquire --profile africa     # about 114,000 records; resumable (fetch cache)
virotrace normalise --profile africa   # then the same stages as above, in the same order
```

The Africa query is narrowed to records that look like *pol* (`acquire.extra_term` in
`config/profiles/africa.yaml`; the clause and its recall check are in
[docs/03](docs/03-data-acquisition.md#3-sources-and-how-each-is-queried)), because the unnarrowed query matches about 187,000 records,
most of them not *pol*. Entrez returns at most 100,000 ids per query, so the search is split
by publication date and checked against the total NCBI reports. Each build keeps its own
snapshot (`data/raw/<scope>.current`), so the Kenya and Africa corpora never mix.

The Africa build of November 2026 held 88,882 QC-passed sequences, 40,217 individuals after
the clone and individual policies, and a 4,577-individual tree set (capped at 25 per country,
two-year bin and subtype). Subtype trees are built for A1, C, D and G.

What the current results say, and do not say:

- **Clusters** are a descriptive statistic at a 1.5% TN93 threshold, not transmission chains.
  SDRM prevalence is higher inside clusters than outside, but the difference disappears after
  stratifying by study treatment context; it reflects which cohorts are densely sampled.
- **Pooled trees do not have usable temporal signal**; per-subtype trees do, but weakly
  (root-to-tip R-squared 0.13 to 0.46). Where TreeTime cannot estimate a sensible clock (it
  fails outright, or returns a rate or a root date no HIV-1 could have), the pipeline retries
  at the midpoint root with the rate fixed at the gate's regression slope, and says so in the
  report and in the Auspice title. In the Africa build all four dated trees are of that kind.
  Their root dates mostly reflect the sampling window: use them for relative timing, and quote
  a range across the clock-sensitivity table, never a point estimate. Subtype A1 fits no
  plausible rate across the prior range.
- **Pooled prevalence is not a regional estimate.** In the Africa build South Africa is about
  40% of the pooled sequences, treatment context is annotated for only about 15% of
  individuals, and several countries' figures come from single treatment-failure cohorts.
- **Subtyping and recombination are screens**, not RDP or 3SEQ. The point-mutation caller does
  not detect the T69 insertions or the K70 deletion.
- **Individual identification** relies on isolate-name patterns in
  `qc.individual_id_patterns`. Studies whose naming is not covered are treated as one person
  per record, so duplicate people can remain. Review `near_identical_pairs.tsv` after a run.

Reference data lives in `config/reference/`: the pinned HXB2 record (`hxb2.gb`) and the
WHO surveillance mutation lists. `sdrm_who.tsv` is the original WHO 2009 list (the
headline metric, comparable with published estimates); `sdrm_stanford_2026-09-25.tsv` is
the current Stanford list (2009 plus later additions and INSTI), used as a sensitivity
analysis. Both are generated from the transcriptions in `config/reference/sources/` with
`make sdrm`, which also checks every wild-type residue against HXB2. Verify a
transcription against its source before reporting a result.

## Licence

MIT. Sequence data is used under each source's terms and is not redistributed in bulk —
the snapshot manifests record the exact queries so the corpus can be rebuilt.
