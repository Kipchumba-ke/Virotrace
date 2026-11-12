# Virotrace — Data Acquisition and Dataset Design

---

## 1. The honest starting position

Before any code: **the dataset determines what Virotrace can legitimately claim**, and the public
HIV sequence record for Kenya has three properties that must shape the design rather than be
discovered halfway through.

1. **It is a convenience sample.** Sequences arrive in study-shaped batches — a cohort here, a
   resistance survey there. Composition reflects who got sequenced and published, not the epidemic.
2. **Metadata is thin and uneven.** Collection year is common, collection month less so, subnational
   location rare. Treatment status is usually absent or ambiguous.
3. **Volume is modest.** This is a constraint on statistical power and on how finely results can be
   stratified. Plan for hundreds-to-low-thousands of Kenyan *pol* sequences, not tens of thousands.

Design consequence: **quantify the corpus in Phase 1 before committing to Phase 3.** The milestone
gate in `01-technical-specification.md` §8 exists for this reason. If sequences with usable dates are
too few for a time-tree, that is a finding to report, not a problem to paper over.

## 2. Access-tier policy

Virotrace recognises three tiers and v1 uses only Tier 1.

| Tier | Description | v1 use |
|---|---|---|
| **1 — Open** | LANL HIV Sequence Database, NCBI Virus / GenBank, published supplementary data with open licences | **Yes.** All v1 analysis |
| **2 — Registered** | Resources needing an account or acceptance of terms but no institutional agreement | Only if terms permit redistribution of derived results; recorded per-source |
| **3 — Controlled** | National surveillance data, PANGEA, any data-transfer-agreement dataset | **No.** Out of scope. An adapter interface exists so it *could* be added under a proper agreement, but no v1 result depends on it |

Every record carries an `access_tier` field. Any future mixing of tiers must be visible in the data,
not buried in a changelog. The repo will never contain Tier 3 data, and `.gitignore` plus a CI check
guard `data/` from being committed at all.

## 3. Sources and how each is queried

| Source | Gives | Query axis | Adapter notes |
|---|---|---|---|
| **LANL HIV Sequence Database** (hiv.lanl.gov) | Curated HIV sequences with subtype annotations, country, year; strong search interface | country, subtype, genomic region, year range | Best single source for a curated starting corpus. Record the exact search parameters in the snapshot manifest — this is the reproducibility anchor, since re-querying later returns a superset |
| **NCBI Virus / GenBank** (HIV-1, taxid 11676) | Broadest coverage; full records with `/country`, `/collection_date` qualifiers | Entrez query by taxon + country + region length | Programmatic and stable; use for completeness and for records LANL hasn't curated |
| **Stanford HIVDB published datasets** | Genotype datasets used to validate resistance algorithms | n/a | Not for prevalence estimation — used as **positive controls** for the DRM caller (see `07-testing-strategy.md` §4) |
| **GADM / national admin boundaries** | Kenya county geometries for mapping | n/a | Reference data, pinned version, only applied where a record actually supplies subnational location |
| **WHO HIV drug resistance reports** | Published regional DRM prevalence | n/a | Not ingested. Used in the report's discussion as an external sanity check on magnitude |

**Narrowing a very large query.** The Africa query matches about 187,000 records, and in a
500-record trial only about a third survived QC, mostly because the rest do not cover PR and RT.
`acquire.extra_term` adds an Entrez clause (ANDed, so it can only remove records) to the query,
for example a clause that keeps records naming pol. Compare counts on the NCBI nuccore website
with and without the clause before trusting it, and check a trial's QC funnel: a clause that is
too tight silently drops pol records. A query over 100,000 ids is split by publication year
automatically.

**Running the LANL source.** LANL has no query API, so the adapter reads files you download:
at hiv.lanl.gov run the Sequence Search (HIV-1, pol region, the countries of your scope), download
the result as FASTA with the default header `Subtype.Country.Year.Name.Accession`, put the files in
`data/external/lanl/` (git-ignored: LANL's terms do not permit redistribution), and add `lanl` to
`acquire.sources`. The snapshot manifest records each file's name and SHA256 in place of a query
string. LANL headers carry no title or PubMed id, so a record found only in LANL has no study key
and counts as its own individual; where GenBank holds the same record, deduplication keeps the
GenBank copy. Re-downloading later returns a superset, so keep the files you built from.

Deduplication across LANL and GenBank is mandatory and non-trivial: match on accession first, then on
sequence identity, then on `(strain, country, year)` near-duplicates. The dedup decision for every
collapsed pair is logged.

## 4. Metadata schema

One table, one row per retained sequence. `*` = required.

| Field | Type | Controlled vocabulary / format | Notes |
|---|---|---|---|
| `strain`* | string | `^[A-Za-z0-9_.\-]+$`, unique | Canonical ID used in FASTA headers and all joins |
| `accession`* | string | — | Source accession |
| `source`* | enum | `lanl` \| `genbank` | Which adapter produced it |
| `access_tier`* | enum | `1` \| `2` \| `3` | Per §2 |
| `country`* | string | ISO 3166-1 name, normalised | Free-text country strings are messy; normalise through a lookup, log every unmapped value |
| `country_code`* | string | ISO 3166-1 alpha-3 | For joins and maps |
| `region`* | enum | `eastern_africa` \| `southern_africa` \| `western_africa` \| `middle_africa` \| `northern_africa` \| `non_africa` | UN M49 grouping, derived not supplied |
| `admin1` | string | Kenyan county name, normalised | **Nullable and usually null.** Populated only from an explicit record field |
| `admin1_resolved` | bool | — | Explicit flag so maps can report what fraction is country-only |
| `date`* | string | `YYYY-MM-DD` \| `YYYY-MM-XX` \| `YYYY-XX-XX` | Nextstrain-style ambiguity notation; never invent precision |
| `date_precision`* | enum | `day` \| `month` \| `year` \| `unknown` | Drives whether a record is eligible for the time-tree |
| `year`* | int | 1980–current | Derived |
| `subtype_reported` | string | — | As supplied by the source, verbatim |
| `subtype_called` | string | — | Virotrace's own call |
| `subtype_method` | string | — | Which method produced `subtype_called` |
| `subtype_agreement` | enum | `agree` \| `disagree` \| `reported_missing` | Surfaced in QC, not hidden |
| `recombinant_flag` | enum | `none` \| `suspected` \| `confirmed` \| `not_tested` | From the recombination screen |
| `gene_coverage_pr` | float | 0–1 | Fraction of protease covered |
| `gene_coverage_rt` | float | 0–1 | Reverse transcriptase |
| `gene_coverage_in` | float | 0–1 | Integrase |
| `ambiguity_frac`* | float | 0–1 | Non-ACGT fraction |
| `hypermutation_flag` | bool | — | APOBEC-context G→A screen result |
| `qc_status`* | enum | `pass` \| `excluded` \| `quarantined` | |
| `qc_reason` | string | rule ID | Null iff `pass` |
| `treatment_status` | enum | `naive` \| `experienced` \| `unknown` | **Default `unknown`.** Read only from an explicit statement in the record's own text (title, definition, notes) by keyword (`curate/study.py`); contradictory statements stay `unknown`. Never inferred from mutations, which would be circular. A keyword hint, not clinical metadata |
| `study_key` | string | `PMID:…` \| `BP:…` \| `REF:author\|title` \| `unknown_study` | Groups records of one study; used by `study_summary.tsv` to show when a stratum is one study |
| `snapshot_id`* | string | — | Provenance |

`treatment_status` deserves emphasis: the temptation to infer treatment exposure from the presence of
DRMs, then report DRM prevalence by treatment status, is circular reasoning. Virotrace refuses the
inference and reports `unknown` honestly, which means **pretreatment vs. acquired resistance cannot
be cleanly separated in v1** — stated as a limitation rather than finessed.

## 5. Cohort definition

**Target cohort (scope = `kenya`)**
- Inclusion: HIV-1; `country_code == KEN`; *pol* coverage ≥ configured minimum for at least one of
  PR/RT/IN; `date_precision` in {day, month, year}; QC pass.
- Exclusion: sequences < minimum length; ambiguity fraction above threshold; internal stop codon or
  frameshift in the target gene; hypermutation flagged; laboratory-adapted or clonal strains; exact
  duplicates.

**Context cohort (scope = `africa`)**
- Same inclusion rules, `region` in the African groupings, with stratified subsampling so that
  Kenya is not swamped by better-sequenced countries.

**Time-tree eligible subset**
- Additionally requires `date_precision` in {day, month} **or** year-only with the ambiguity encoded
  for TreeTime. Record the size of this subset prominently — it is the real n for anything temporal.

## 6. Subsampling design

Random downsampling would bias composition toward whichever country or study deposited most. Instead:

```
stratum = (country_code, year_bin, subtype_called)
max_per_stratum  = config.subsample.max_per_stratum
min_per_stratum  = config.subsample.min_per_stratum   # strata below this are kept whole
priority         = date_precision, then gene coverage, then lower ambiguity_frac
seed             = config.seed
```

Within a stratum, records are ranked by `priority` and the top `max_per_stratum` retained. Ties broken
by sorted `strain` so the result is deterministic. The realised stratum table is written to
`results/<scope>/qc/subsample_strata.tsv` and plotted in the report — a reviewer can see exactly what
the tree was built from.

For `scope=kenya`, caps will usually not bind; they exist so `scope=africa` stays inside the laptop
envelope with known composition rather than unknown truncation.

## 7. Provenance and versioning

Each acquisition writes `data/raw/<snapshot_id>/manifest.json`:

```json
{
  "snapshot_id": "2026-11-05_a3f19c",
  "created_utc": "2026-11-05T06:12:44Z",
  "sources": [
    {
      "name": "genbank",
      "query": "txid11676[Organism] AND Kenya[Country] AND 1000:4000[SLEN]",
      "retrieved": 1842,
      "files": [{"path": "genbank.gb", "sha256": "..."}]
    }
  ],
  "tool_versions": {"entrez": "...", "virotrace": "0.1.0"},
  "config_hash": "..."
}
```

`snapshot_id` = date + short hash of the source queries. Results are always labelled with the snapshot
they came from. Because public databases grow, **a rerun months later will not reproduce the same
corpus** — reproducibility is therefore defined as *the same snapshot + the same config → the same
result*, and the snapshot manifest is what makes that checkable. Snapshots are gitignored; manifests
are committed.

## 8. Licensing and attribution

Sequence data is used under each source's terms, with accessions retained throughout so every
downstream row is traceable to its deposit. The report includes a data-sources section with retrieval
dates and query strings. Derived tables are released under the repository's licence; underlying
sequence data is not redistributed in bulk — the manifest plus the query lets anyone rebuild it.
