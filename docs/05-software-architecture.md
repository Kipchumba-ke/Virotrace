# Virotrace — Software Architecture and Repository Structure

---

## 1. Design rules

Four rules, applied consistently, are what make this repo readable to a stranger:

1. **Snakemake orchestrates; Python computes.** Rules contain no logic beyond input/output wiring and
   a single call into the package. Any `run:` block longer than three lines becomes a function in
   `src/virotrace/`.
2. **Every stage is independently runnable.** `virotrace <stage> --config config.yaml` works without
   Snakemake, so stages can be debugged, tested, and demonstrated in isolation.
3. **No I/O in logic functions.** Pure functions take and return data structures; thin wrappers
   handle files. This is what makes the test suite cheap.
4. **Config is data, not code.** No module reads a hard-coded threshold. Parameters arrive as a
   validated settings object.

## 2. Repository structure

```
virotrace/
├── README.md                       demo path, 10-minute quickstart, disclaimer
├── LICENSE
├── CITATION.cff
├── pyproject.toml                  package + tool config (ruff, mypy, pytest)
├── environment.yml                 conda env (bio tools)
├── environment.lock.yml            pinned, regenerated deliberately
├── Makefile                        make setup | dev | kenya | africa | test | lint | clean
├── .pre-commit-config.yaml
├── .gitignore                      data/, results/, auspice/, .snakemake/
│
├── config/
│   ├── config.yaml                 the single source of parameters
│   ├── sources.yaml                source queries per scope
│   ├── profiles/
│   │   ├── dev.yaml                toy dataset, ~200 seqs, <10 min
│   │   ├── kenya.yaml              primary build
│   │   └── africa.yaml             context build, overnight
│   ├── schema/
│   │   ├── config.schema.json      JSON Schema for Snakemake's own validation
│   │   └── metadata.py             pandera DataFrameSchema (the L1/L2 contract)
│   └── reference/
│       ├── hxb2.gb                 pinned reference record
│       ├── subtype_references.fasta
│       ├── sdrm_who.tsv            WHO surveillance DRM list, versioned
│       └── lab_strains.txt         exclusion name list
│
├── workflow/
│   ├── Snakefile                   includes rules, defines `rule all`
│   ├── rules/
│   │   ├── acquire.smk             L0
│   │   ├── curate.smk              L1 + L2
│   │   ├── drm.smk                 L3a
│   │   ├── phylo.smk               L3b
│   │   └── export.smk              L4
│   └── envs/                       per-rule conda envs (iqtree.yaml, mafft.yaml, …)
│
├── src/virotrace/
│   ├── __init__.py                 __version__
│   ├── cli.py                      typer app; one subcommand per stage
│   ├── settings.py                 pydantic models for config; single load point
│   ├── provenance.py               provenance sidecar writer, config hashing
│   ├── logging.py                  structured logging setup
│   ├── io/
│   │   ├── fasta.py                read/write, header policy
│   │   ├── tables.py               TSV read/write with dtype enforcement
│   │   └── manifest.py             snapshot manifests, sha256
│   ├── acquire/
│   │   ├── base.py                 SourceAdapter protocol
│   │   ├── genbank.py
│   │   ├── lanl.py
│   │   └── snapshot.py
│   ├── curate/
│   │   ├── normalise.py            source records → schema rows
│   │   ├── geo.py                  country normalisation, region mapping, admin1
│   │   ├── dates.py                date parsing + ambiguity notation + precision
│   │   ├── dedup.py
│   │   ├── qc.py                   QC-1..QC-7, returns (kept, ledger)
│   │   ├── hypermut.py             APOBEC context test
│   │   ├── align.py                MAFFT/augur wrapper, gene slicing
│   │   ├── subtype.py              placement + cross-check + agreement
│   │   ├── recombination.py
│   │   └── subsample.py            stratified, deterministic
│   ├── drm/
│   │   ├── algorithm.py            HIVdb algorithm execution, version pinning
│   │   ├── mutations.py            alignment → gene-local AA mutation list
│   │   ├── mixtures.py             ambiguity → AA set, mixture policy
│   │   ├── sdrm.py                 WHO SDRM tagging
│   │   └── prevalence.py           strata, Wilson CIs, suppression, trend
│   ├── phylo/
│   │   ├── tree.py                 IQ-TREE wrapper
│   │   ├── refine.py               TreeTime/augur refine
│   │   ├── temporal.py             root-to-tip regression, the gate
│   │   ├── traits.py               mugration, admin1 gating
│   │   └── clusters.py             TN93, cluster stats
│   ├── export/
│   │   ├── auspice.py              v2 JSON assembly, DRM colourings
│   │   └── report.py               Jinja2 render, figures
│   └── viz/
│       └── figures.py              prevalence plots, QC plots, composition plots
│
├── tests/
│   ├── conftest.py                 fixtures
│   ├── fixtures/
│   │   ├── toy_pol.fasta           ~20 seqs, planted mutations
│   │   ├── positive_controls.tsv   sequence → expected DRM calls
│   │   ├── hypermut_pair.fasta     edited + unedited parent
│   │   ├── subtype_refs_mini.fasta
│   │   └── metadata_toy.tsv
│   ├── unit/
│   ├── contract/                   schema + data-contract tests
│   ├── integration/                multi-stage on fixtures
│   └── regression/
│       └── golden/                 expected outputs, refreshed deliberately
│
├── docs/
│   ├── 01-technical-specification.md
│   ├── 02-system-architecture.md
│   ├── 03-data-acquisition.md
│   ├── 04-bioinformatics-methodology.md
│   ├── 05-software-architecture.md
│   ├── 06-configuration.md
│   ├── 07-testing-strategy.md
│   └── decisions/                  ADRs: 0001-why-snakemake.md, 0002-exclude-recombinants-from-tree.md, …
│
├── notebooks/                      exploratory only; never imported by src/
│
└── .github/workflows/
    ├── ci.yml                      lint, type, test, dry-run DAG
    └── dev-build.yml               full dev-profile build on a schedule
```

## 3. Module responsibilities and boundaries

| Package | Owns | Must not |
|---|---|---|
| `acquire/` | Talking to the outside world; writing immutable snapshots | Interpret or filter biology |
| `curate/` | Everything between a raw record and an analysis-ready alignment | Know about resistance or trees |
| `drm/` | Mutation extraction and resistance interpretation | Know about phylogenetics |
| `phylo/` | Trees, time, traits, clusters | Know how DRMs were called — it joins a TSV |
| `export/` | Assembling outputs for humans | Compute anything new |
| `io/`, `settings/`, `provenance/` | Cross-cutting plumbing | Contain domain logic |

`drm/` and `phylo/` not importing each other is the structural expression of the phased plan: the DRM
layer can be finished, tested, and shipped while `phylo/` is still empty stubs.

## 4. CLI design

```
virotrace --help
  acquire      --scope kenya --out data/raw
  normalise    --snapshot <id>
  qc           --in ... --out ... --ledger ...
  align
  subtype
  recombination
  subsample
  drm call     --aln ... --out drm_table.tsv
  drm prevalence
  phylo tree
  phylo refine
  phylo traits
  phylo clusters
  export auspice
  export report
  validate config     # schema-check config without running anything
  info                # versions, pinned reference data, env fingerprint
```

Conventions: every command takes `--config`, writes a provenance sidecar, logs to stderr as
structured JSON when `--json-logs` is set, and exits non-zero on contract violation. `validate config`
exists so a bad config fails in a second rather than forty minutes into a tree.

## 5. Error handling

Three exception classes, each with a defined caller response:

| Class | Meaning | Behaviour |
|---|---|---|
| `ConfigError` | Config invalid or internally inconsistent | Abort immediately, print the offending key and the constraint |
| `ContractError` | A data contract was violated (schema, coordinate, frame) | Abort the stage; this is a bug or bad reference data, not bad input |
| `RecordError` | One record is unusable | Quarantine the record, increment a counter, continue |

The split matters: a pipeline that aborts on one malformed GenBank record is unusable, and one that
silently swallows a schema change is untrustworthy. `RecordError` counts are summarised at the end of
every stage, and a configurable `max_record_error_rate` escalates to abort if too much is being
quarantined — silent mass-quarantine being the failure mode that split otherwise invites.

## 6. Conventions

- **Naming:** `snake_case` throughout; files named `<stage>.<scope>.<ext>`; no spaces, no dates in
  filenames except snapshot IDs.
- **Types:** full type hints; `mypy --strict` on `src/virotrace/` in CI.
- **Lint/format:** ruff (lint + format), line length 100.
- **Logging:** `structlog`-style key-value logging; every stage logs `n_in`, `n_out`, `n_excluded`,
  duration. Those four numbers at every stage boundary are the cheapest debugging tool available.
- **Versioning:** semver on the package; the **metadata schema has its own version** and a bump is a
  breaking change.
- **ADRs:** decisions with real trade-offs get a short `docs/decisions/NNNN-*.md`. The two to write on
  day one are *why Snakemake* and *why recombinants leave the tree but stay in the DRM layer* — a
  reviewer will ask both.

## 7. Environment management

- `environment.yml` for the dev environment; per-rule envs in `workflow/envs/` so Snakemake can build
  minimal environments for heavy tools and so a tool upgrade is scoped to one rule.
- `environment.lock.yml` regenerated deliberately and committed; CI builds from the lock.
- Optional `Dockerfile` for a one-command reviewer path.
- `virotrace info` prints the resolved versions of every external binary, which goes into provenance.
