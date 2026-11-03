# ADR 0002 — Recombinants leave the tree but stay in the DRM layer

**Status:** accepted
**Date:** 2026-11-04

## Context

Kenya's HIV-1 epidemic is predominantly subtype A1 with D, C, and a substantial burden of
circulating and unique recombinant forms. Recombination means different parts of a genome have
different evolutionary histories.

A bifurcating phylogeny assumes a single history for the whole alignment. Genotypic resistance
interpretation makes no such assumption: it reads codons at defined positions and compares them
against an algorithm.

So one analysis is invalidated by recombination and the other is not — and a single
include/exclude decision applied to both would damage one of them.

## Decision

Screen for recombination once, then apply different policies per layer:

- **Phylo layer:** exclude sequences with confirmed breakpoints
  (`recombination.tree_policy: exclude_confirmed`).
- **DRM layer:** retain all sequences (`recombination.drm_policy: keep_all`).

Report the excluded fraction in both outputs, and state in the report that the phylogeny
represents the non-recombinant subset while the resistance estimates do not.

## Consequences

**Positive.** The tree does not silently violate its own model. The resistance estimate is not
biased by dropping a non-random, lineage-correlated chunk of the corpus — recombinants are not
randomly distributed with respect to subtype or region, so excluding them from prevalence would
shift the estimate in an unpredictable direction.

**Negative.** The two layers are computed on different sequence sets, so n differs between the
tree and the prevalence tables. That is a genuine cost: it has to be stated clearly in the report,
and "how many sequences does Virotrace analyse?" has two correct answers.

**Rejected alternatives.**
- *Exclude recombinants everywhere.* Biases the resistance estimate for no methodological gain.
- *Keep recombinants everywhere.* Produces a tree whose topology is partly an artefact.
- *Partition the alignment by breakpoint and build per-partition trees.* Methodologically the
  best answer and out of scope for v1 on a laptop; recorded as a scale-up path.

## Verification

`recombination.tree_policy` and `recombination.drm_policy` are independent config keys, and the
integration test asserts the DRM input set is a superset of the tree input set.
