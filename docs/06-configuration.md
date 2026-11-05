# Virotrace — Configuration-Driven Analysis

---

## 1. Principle

**No analytical parameter appears in Python or in a Snakemake rule.** Thresholds, model strings,
caps, seeds, and scope definitions all live in YAML, are validated against a schema before anything
runs, and are hashed into every output's provenance.

Three things follow, and they are the reason to do it this way:
- a reviewer can read the entire analysis configuration on one screen;
- changing the analysis never means editing code, so "what did we change" is a config diff;
- the config hash makes "which parameters produced this figure" answerable after the fact.

## 2. Layering

```
config/config.yaml                 base — every key present with a default
        ↓ overridden by
config/profiles/<profile>.yaml     dev | kenya | africa
        ↓ overridden by
--config key=value                 CLI, for one-off experiments
```

Later layers override earlier ones key-by-key. Deep merge, so a profile can set
`qc.max_ambiguity_frac` without restating the rest of `qc`. The fully resolved config is written to
`results/<scope>/resolved_config.yaml` next to the outputs — the config you can read in the repo is
the input, the one next to the results is the record.

## 3. Base configuration

```yaml
# config/config.yaml
version: 1
schema_version: 1          # metadata contract version; a bump is breaking

scope: kenya               # kenya | africa  (profiles set this)
seed: 20261104
threads: 4                 # laptop default; profiles may raise
snapshot_id: null          # null = use latest; pin a string for reproducible reruns

paths:
  raw: data/raw
  interim: data/interim
  processed: data/processed
  reference: config/reference
  results: results
  auspice: auspice
  reports: reports

scopes:
  kenya:
    countries: [KEN]
    regions: []
  africa:
    countries: []
    regions: [eastern_africa, southern_africa, western_africa, middle_africa, northern_africa]

reference:
  accession: K03455              # HXB2
  record: hxb2.gb
  genes:                         # HXB2 nucleotide coordinates; asserted in tests
    pol: [2085, 5096]
    pr:  [2253, 2549]
    rt:  [2550, 4229]
    in:  [4230, 5096]

acquire:
  sources: [genbank, lanl]
  min_length: 1000
  max_length: 10000
  email: null                    # required by NCBI; set via env VIROTRACE_NCBI_EMAIL
  retries: 3
  cache: true

curate:
  dates:
    min_year: 1985
    max_year: 2026
    allow_year_only: true        # encoded as YYYY-XX-XX
  geo:
    country_lookup: country_map.tsv
    require_admin1: false
  dedup:
    identity_threshold: 1.0      # 1.0 = exact only; lower to collapse near-duplicates
    collapse_strategy: keep_best_coverage

qc:
  min_gene_coverage: 0.80        # QC-1, fraction of the gene
  max_ambiguity_frac: 0.05       # QC-2
  exclude_internal_stop: true    # QC-3
  exclude_frameshift: true       # QC-4
  hypermutation:
    enabled: true                # QC-5
    p_threshold: 0.05
    action: exclude              # exclude | flag_only
  min_identity_to_reference: 0.75  # QC-0: below this the sequence is not pol
  min_class_position_coverage: 1.0 # QC-1: share of a class's SDRM positions covered
  exclude_lab_strains: true      # QC-7
  clone_policy: exclude          # exclude | keep_all: clone-derived records leave the
  individual_policy: one_per_individual   # keep_all | one_per_individual: one sequence per
                                          # identifiable person (earliest, least ambiguous)
  individual_id_patterns: ['^(RV217_\d+)_']  # regexes on the isolate name; group 1 = person id
                                 # prevalence and tree layers, stay in the DRM table
  lab_strain_list: lab_strains.txt
  max_record_error_rate: 0.25    # abort if more than this fraction quarantines

align:
  tool: mafft
  reference_anchored: true
  codon_aware: true
  mask_gap_frac: 0.90            # tree alignment only; never applied to DRM calling

subtype:
  methods: [placement, external]
  reference_panel: subtype_references.fasta
  min_support: 0.70
  on_disagreement: report        # report | prefer_placement | prefer_external

recombination:
  enabled: true
  method: 3seq
  p_threshold: 0.05
  tree_policy: exclude_confirmed # exclude_confirmed | exclude_suspected | keep_all
  drm_policy: keep_all           # recombination does not invalidate a genotype

subsample:
  enabled: true
  strata: [country_code, year_bin, subtype_called]
  year_bin_width: 2
  max_per_stratum: 40
  min_per_stratum: 1
  priority: [date_precision, gene_coverage, ambiguity_frac]

drm:
  algorithm: hivdb
  algorithm_version: null        # null = pinned bundled version; recorded in provenance
  classes: [NRTI, NNRTI, PI, INSTI]
  require_gene_coverage: true    # per-class denominators
  mixtures:
    policy: any                  # any | majority | strict
  sdrm:
    list: sdrm_who.tsv
    list_version: "2009"         # headline prevalence list (original WHO 2009)
    sensitivity_list: sdrm_stanford_2026-09-25.tsv   # optional; reported separately with its
    sensitivity_list_version: "stanford-2026-09-25"  # own denominators, never the headline
  cross_check:
    enabled: false               # optional external re-call
    cache: true
    fail_on_disagreement: true

prevalence:
  strata: [year_bin, drug_class]
  also_by_subtype: true
  ci_method: wilson
  ci_level: 0.95
  min_stratum_n: 10              # below this: "n too small", not a plotted point
  trend_test: cochran_armitage
  study_annotations: study_annotations.tsv   # hand-curated study_key -> treatment context
  year_bin_width: 5              # bins start at multiples of this: 2000-2004, 2005-2009, ...

phylo:
  enabled: true                  # false ships the DRM-only build
  tree:
    tool: iqtree2
    model: GTR+F+R4
    model_finder: false
    bootstrap: ufboot
    replicates: 1000
    fast_mode: false
  refine:
    enabled: true
    clock_rate: null             # null = estimate
    clock_rate_prior: [0.0008, 0.0030]   # subs/site/year, a prior range to validate against
    clock_filter_iqd: 4
    date_confidence: true
  temporal_gate:
    enabled: true
    min_temporal_r2: 0.10
    require_positive_slope: true
    date_randomisation: false    # true at kenya scope once n is known
  traits:
    columns: [country_code, region]
    admin1:
      enabled: auto              # auto = only if resolved fraction clears the threshold
      min_resolved_frac: 0.30
  clusters:
    enabled: true
    distance: tn93
    threshold: 0.015
    min_size: 2

report:
  min_stratum_n: 10
  include_limitations: true       # not switchable to false in profiles; see §6
  figures: [prevalence_by_year, prevalence_by_class, subtype_composition,
            sampling_composition, qc_waterfall, cluster_size_distribution]

logging:
  level: INFO
  json: false
```

## 4. Profiles

```yaml
# config/profiles/dev.yaml  —  CI and first-run smoke; must finish in minutes
scope: kenya
threads: 2
acquire:
  sources: [fixture]            # reads tests/fixtures, no network
subsample:
  max_per_stratum: 5
phylo:
  tree: {model: GTR, bootstrap: none, fast_mode: true}
  refine: {enabled: false}
  clusters: {enabled: false}
prevalence:
  min_stratum_n: 1              # fixtures are tiny; real builds must not do this
report:
  figures: [prevalence_by_class, qc_waterfall]
```

```yaml
# config/profiles/kenya.yaml  —  the primary build
scope: kenya
threads: 4
subsample:
  enabled: false                # expected to fit without capping; revisit after P1 counts
phylo:
  temporal_gate: {date_randomisation: true}
```

```yaml
# config/profiles/africa.yaml  —  context build, overnight
scope: africa
threads: 4
subsample:
  enabled: true
  max_per_stratum: 25
phylo:
  tree: {fast_mode: true, replicates: 200}
report:
  min_stratum_n: 15
```

`dev.yaml` lowering `min_stratum_n` to 1 is a deliberate, visible compromise so fixtures produce
output; the comment exists so nobody copies it into a real profile.

## 5. Validation

Two independent gates:

1. **Snakemake's own** `validate(config, "config/schema/config.schema.json")` at workflow load.
2. **Pydantic settings models** in `src/virotrace/settings.py`, which enforce what JSON Schema cannot:

```python
class SubsampleConfig(BaseModel):
    enabled: bool = True
    max_per_stratum: int = Field(40, ge=1)
    min_per_stratum: int = Field(1, ge=0)

    @model_validator(mode="after")
    def _ordering(self):
        if self.min_per_stratum > self.max_per_stratum:
            raise ConfigError("subsample.min_per_stratum exceeds max_per_stratum")
        return self
```

Cross-field rules worth encoding up front:
- `drm.require_gene_coverage` true ⇒ `qc.min_gene_coverage > 0`;
- `phylo.refine.enabled` true ⇒ `curate.dates.allow_year_only` handled, not silently dropped;
- `traits.admin1.enabled == true` (not `auto`) ⇒ warn loudly, because forcing subnational
  reconstruction past the resolved-fraction gate is exactly the overclaim the spec forbids;
- `report.min_stratum_n < 5` ⇒ warn unless the profile is `dev`.

`virotrace validate config --profile africa` runs both gates in about a second.

## 6. What is *not* configurable

Deliberate omissions, because configurability here would let a future run quietly become
indefensible:

- The exclusion ledger cannot be turned off.
- `report.include_limitations` cannot be set false.
- Provenance writing cannot be disabled.
- `treatment_status` cannot be inferred from DRMs by any config path.
- Wilson CIs cannot be replaced with intervals that misbehave at small n.

## 7. Reproducibility mechanics

| Mechanism | Implementation |
|---|---|
| Config hash | SHA256 of the resolved config, in every provenance sidecar and in output filenames at `africa` scope |
| Data pinning | `snapshot_id` pinned in a profile reproduces the corpus; manifests committed |
| Seeds | One root `seed`; per-stage seeds derived deterministically from `(seed, stage_name)` so adding a stage doesn't shift another stage's stream |
| Tool versions | `virotrace info` output embedded in provenance |
| Determinism | Inputs sorted before iteration; ties broken on `strain`; no reliance on dict or filesystem order |
| Resolved config | Written beside results, so a figure is always traceable to its exact parameters |
