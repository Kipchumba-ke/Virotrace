# Virotrace — Bioinformatics Methodology

Every parameter named here lives in `config/config.yaml`, not in code. Values below are the
**defaults**, chosen to be defensible and laptop-feasible, and each is justified so that changing it
is an informed decision.

---

## 1. Reference and coordinate system

All coordinates are **HXB2** (GenBank K03455), the standard HIV-1 numbering reference.

| Region | HXB2 nucleotide range | Codons | Protein length |
|---|---|---|---|
| *pol* (full, incl. terminal stop) | 2085–5096 | 1004 | 1003 aa |
| Protease (PR) | 2253–2549 | 99 | 99 aa |
| Reverse transcriptase (RT) | 2550–4229 | 560 | 560 aa |
| Integrase (IN, coding only) | 4230–**5093** | 288 | 288 aa |

> **Integrase end coordinate.** Integrase is widely quoted as 4230–5096, which is 289 codons
> because it includes *pol*'s terminal stop codon. Virotrace uses 4230–5093 so that
> `(end - start + 1) / 3` equals the 288-residue protein length, and gene-local position 288 is
> the last real residue rather than a stop.
>
> **Build-time check:** these ranges are asserted in `tests/contract/test_reference_and_sdrm.py`,
> against the annotated HXB2 record once `config/reference/hxb2.gb` is pinned. The first version of
> this table said 5096 and the test caught it — which is the point of making the claim falsifiable
> rather than trusting the document.

Resistance mutations are reported in **gene-local amino acid numbering** (e.g. RT M184V = RT codon
184), which is the convention every resistance database and clinician uses. Nucleotide positions are
retained internally for alignment work only.

## 2. Quality control

Applied in order; each rule writes to the exclusion ledger with the observed value.

| ID | Rule | Default | Rationale |
|---|---|---|---|
| QC-1 | Minimum ungapped length in target region | ≥ 80% of the gene | Short fragments produce spurious "wild-type" calls at uncovered positions — the most dangerous silent error in DRM work |
| QC-2 | Ambiguity fraction (non-ACGT) | ≤ 5% | High ambiguity usually means mixed populations or poor sequencing; ambiguous codons also inflate false DRM calls |
| QC-3 | Internal stop codon in target gene | exclude | Indicates frameshift, misassembly, or a defective provirus |
| QC-4 | Frameshift (indel not a multiple of 3 relative to reference) | exclude | Breaks codon-aware translation |
| QC-5 | APOBEC hypermutation | flag, exclude by default | G→A in APOBEC3G/F context (GG→AG, GA→AA) generates mutations that look like resistance but are host-editing artefacts. A Hypermut-style context-ratio test with a Fisher exact p-value threshold |
| QC-6 | Exact or near-exact duplicates | collapse, log | Avoids pseudo-replication inflating both prevalence and cluster support |
| QC-7 | Laboratory strains / clones | exclude by name list | HXB2, NL4-3 and relatives are not epidemiological samples |

**Why QC-5 is non-negotiable:** hypermutated sequences are a well-known source of spurious
drug-resistance signal in sequence-database surveys. Omitting this screen would make the headline
prevalence number unreliable, so it runs by default and the excluded count is reported.

**Why coverage is tracked per gene:** a sequence may cover RT well and integrase not at all. Calling
"no INSTI resistance" from a sequence with no integrase is a false negative by construction.
Denominators in §6 are therefore **per drug class**, computed from per-gene coverage — not one global n.

## 3. Alignment

- Reference-anchored alignment with MAFFT via `augur align --reference-sequence HXB2`, keeping
  reference coordinates so gene slicing is exact.
- **Codon-aware handling** for the coding regions: translate, check frame, and treat the alignment as
  codons rather than independent nucleotides when extracting mutations. A nucleotide-level indel that
  shifts frame invalidates every downstream amino-acid call, so frame is validated per sequence, not
  assumed.
- Gene sub-alignments (`pr`, `rt`, `in`) are sliced from the *pol* alignment by HXB2 coordinate.
- Columns with gap fraction above a configured threshold are masked for tree inference only — never
  for DRM calling, where a masked position must become "no call", not "wild type".

That last distinction matters: masking is an acceptable tree-building convenience and an unacceptable
resistance-calling shortcut. The pipeline keeps the two alignments separate for this reason.

## 4. Subtype assignment and recombination

**Subtyping — two independent calls, disagreement surfaced.**
1. *Phylogenetic placement:* align the query set with a curated subtype reference panel and place
   sequences on a reference tree; assign by the supported reference clade.
2. *External/algorithmic call:* an established subtyping tool as a cross-check.

Agreement is recorded in `subtype_agreement`. Disagreements are **reported, not resolved by
fiat** — in a recombinant-rich epidemic, disagreement between subtyping methods is itself
informative, and hiding it would misrepresent confidence.

**As implemented (Phase 2).** Placement is *distance-based*: each QC-passed sequence is compared with
the reference panel over sites where both have an unambiguous base, the call is the nearest label,
and its support is the share of site-bootstrap replicates that agree (`subtype.min_support`, default
0.70; below it the call is `unassigned`). The second, independent call is not yet an external tool:
the depositor's own label is compared (`subtype_agreement`: agree / disagree / unreported /
not_comparable) and disagreement is reported, not resolved. Recombination is screened by a
sliding-window subtype scan (`recombination.method: window_scan`): a sequence is `suspected` when
two or more supported segments carry different labels and `confirmed` when each segment re-placed on
its own columns keeps its label at `confirm_support`. Limits to carry into the report: pol (a PR+RT
genotype is about 1.3 kb) is short, so intra-subtype and near-end breakpoints are invisible; a mosaic
of lineages missing from the panel looks like an unassigned sequence, not a recombinant; and this is
not 3SEQ or RDP and claims none of their statistics. The panel is supplied by the user
(`config/reference/subtype_references.fasta`, LANL headers such as `A1.KE.00.name.ACCESSION`).

**Recombination screening.** Kenya's epidemic is predominantly subtype A1 with D, C, and a
substantial burden of circulating and unique recombinant forms. A single bifurcating tree assumes no
recombination, so:
- screen with a breakpoint-detection method (3SEQ or RDP-class) over the *pol* alignment;
- sequences with supported breakpoints are flagged `suspected`/`confirmed`;
- **default policy: exclude confirmed recombinants from the tree, retain them in the DRM layer.**

That asymmetry is deliberate and is the single most important methodological decision in the project:
recombination invalidates the *phylogeny*, not the *resistance genotype*. Excluding recombinants from
DRM prevalence would bias the resistance estimate; including them in the tree would bias the
phylogeny. The excluded fraction is reported in both outputs.

## 5. Drug resistance calling

- **Algorithm:** the Stanford HIVdb genotypic resistance interpretation algorithm, run from a local
  implementation with the algorithm version pinned and recorded in every output row.
- **Classes covered:** NRTI, NNRTI, PI, INSTI — each conditional on adequate coverage of its gene.
- **Two parallel mutation sets, used for different purposes:**

| Set | Purpose | Why separate |
|---|---|---|
| Full HIVdb DRM list + penalty scores | Per-sequence interpretation (susceptible → high-level resistance) | Clinically oriented, sensitive to polymorphism |
| **WHO SDRM list** (surveillance DRMs) | Population prevalence and trend | Deliberately curated to be comparable across populations and to exclude common polymorphisms. This is the correct denominator-and-numerator set for a *surveillance* claim |

Reporting SDRM prevalence as the headline surveillance metric, with full HIVdb interpretation as a
secondary per-sequence view, is what makes the numbers comparable to published regional estimates.
Using HIVdb penalty scores as a prevalence metric would not be.

- **Non-B caveat, operationalised.** HIVdb is best validated on subtype B. For non-B sequences the
  pipeline: (a) records `subtype_called` on every DRM row, (b) separates known subtype-specific
  polymorphisms from resistance mutations where the algorithm distinguishes them, and (c) stratifies
  the report by subtype so a reader can see whether a signal is concentrated in one lineage. The
  caveat is printed next to the table, not in a footnote.
- **Mixtures.** Ambiguous nucleotides encoding two amino acids at a DRM position are reported as
  mixtures with both calls retained, and counted per a configured policy (`any` = resistant if either
  call is a DRM). The policy is a config value because reasonable analysts differ.
- **Cross-check path.** Optionally re-call a sample through the public Sierra service and assert
  agreement, cached locally. Disagreement fails CI rather than being averaged — it means an algorithm
  version drifted.

## 6. Prevalence estimation

For each (scope, year_bin, drug_class, subtype?) stratum:

- Numerator: sequences with ≥1 SDRM for that class.
- **Denominator: sequences with adequate coverage of that class's gene** — never the full corpus.
- Interval: **Wilson score** 95% CI (correct behaviour at small n and near 0%, unlike Wald).
- Suppression: strata with n below `report.min_stratum_n` are reported as "n too small" rather than
  plotted. An unsuppressed 100% from n=2 is the classic surveillance-dashboard embarrassment.
- Trend: Cochran–Armitage test for trend across ordered year bins, with the caveat that **sampling
  composition changes over time**, so a trend may reflect which studies deposited when. Reported as
  descriptive, with the per-year composition table adjacent so the reader can judge.

No model-based adjustment for sampling bias in v1 — the data do not support a credible adjustment,
and an unadjusted number with an honest caveat beats an adjusted number with a hidden assumption.

## 7. Phylogenetics

**Tree inference.** IQ-TREE 2 on `aln.final.fasta`:
- substitution model: `GTR+F+R4` by default (free-rate heterogeneity handles rate variation across
  *pol* better than a single gamma); `ModelFinder` available via config but off by default to control
  runtime.
- support: ultrafast bootstrap, 1000 replicates. At `scope=africa`, `-fast` mode with support from a
  reduced replicate count, recorded in provenance.
- seed fixed from config.

**Temporal calibration.** `augur refine` (TreeTime):
- clock model with a **prior range** rather than a fixed rate, since *pol* clock estimates vary with
  dataset and sampling window; the estimated rate is reported and sanity-checked against published
  ranges rather than asserted from this document.
- `--date-confidence`, outlier detection (`--clock-filter-iqd`) on, with removed tips logged.
- date ambiguity handled through the `YYYY-XX-XX` notation from the metadata schema.

**Temporal-signal gate.** Before the time-tree is accepted:
1. root-to-tip regression on the divergence tree;
2. require R² ≥ `phylo.min_temporal_r2` and a positive slope;
3. date-randomisation check on a subsample at `kenya` scope.

If the gate fails, the pipeline **emits the divergence tree and a stated reason**, and the report says
no time-calibrated inference was made. A time-tree from data without temporal signal looks exactly as
convincing as one from data with it — which is why the gate is automated rather than left to judgement
at writing time.

**Ancestral traits.** `augur traits` mugration over `country_code` and `region`. Subnational traits
are reconstructed **only if** `admin1_resolved` coverage exceeds a configured fraction; otherwise the
geographic layer stops at country, and the legend states the resolved fraction. Discrete-trait
mugration on a convenience sample yields migration patterns that partly reflect sampling intensity —
stated in the report as a limitation of the geographic reconstruction specifically.

**Clustering.** TN93 pairwise distances with threshold `phylo.cluster_threshold` (default 1.5%,
configurable), minimum cluster size configurable. Reported as a descriptive population statistic:
cluster size distribution, subtype composition, DRM enrichment within clusters. **No individual-level
transmission claims** (`01-technical-specification.md` §7).

## 8. Validation

| Check | Method | Fails the build? |
|---|---|---|
| DRM caller correctness | Positive-control sequences with known mutations and known expected calls | Yes |
| Hypermutation screen | Synthetic hypermutated sequence must be flagged; its unedited parent must not | Yes |
| Coordinate integrity | HXB2 gene ranges asserted against the annotated record | Yes |
| Subtype calling | Reference sequences of known subtype must recover their subtype | Yes |
| Alignment frame | Every retained sequence translates without internal stops | Yes |
| Temporal signal | Root-to-tip R² gate | No — downgrades output |
| Subsampling stability | Rerun tree at 2 seeds; key clade support should be qualitatively stable | No — reported |
| External plausibility | Compare subtype distribution and DRM magnitude against published Kenyan/regional estimates | No — discussed in report |

The last row is a judgement check, not a test: if Virotrace's subtype distribution for Kenya looks
nothing like the published picture, the most likely explanation is a pipeline bug, and that
comparison should happen before publishing, not after a reviewer asks.

## 9. Stated limitations (to be carried verbatim into the report)

1. Convenience sampling — results describe deposited sequences, not the Kenyan epidemic.
2. Pretreatment vs. acquired resistance cannot be separated; treatment status is unknown for most records.
3. Geographic resolution is country-level for most records; subnational inference is limited or absent.
4. HIVdb interpretation is best validated on subtype B; Kenya's epidemic is predominantly non-B.
5. Confirmed recombinants (by the distance screen, per the recombination tree policy) are excluded from the tree. Suspected recombinants and sequences with no subtype call stay in, so the phylogeny is not strictly the non-recombinant subset.
6. Migration reconstruction is confounded by uneven sampling intensity across countries and years.
7. Temporal trends may reflect changing study composition rather than epidemic change.
8. Not a clinical tool; no result should inform individual patient management.
