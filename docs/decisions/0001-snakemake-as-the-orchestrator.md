# ADR 0001 — Snakemake as the orchestrator

**Status:** accepted
**Date:** 2026-11-04

## Context

Virotrace is a multi-stage batch analysis run on one laptop, where reproducibility is the headline
feature and the audience includes reviewers who will spend about ten minutes deciding whether the
work is credible. The orchestration choice has to serve: resumability (a tree run that dies at hour
four must not restart from acquisition), config-driven parameterisation, per-rule environments for
awkward bioinformatics binaries, and legibility to someone reading the repo cold.

## Decision

Snakemake, with rules that contain no analytical logic — each shells out to one `virotrace`
subcommand.

## Alternatives considered

| Option | Why not |
|---|---|
| **Nextflow** | Equally capable and arguably better at cluster execution, but heavier to learn and operate for a solo laptop project, and its channel model buys little when the DAG is a linear-ish file pipeline |
| **Make** | No wildcards worth the name, no config integration, no provenance, no per-rule environments |
| **A Python driver script** | Easiest to start, worst to resume. Every re-run recomputes, or grows ad-hoc "if output exists, skip" logic that is the thing a DAG engine already does correctly |
| **Airflow / Prefect** | Built for scheduled services with a server. Deployment surface with no benefit here |

## Consequences

**Positive.** `--dry-run` is a free, always-current architecture diagram and a cheap CI check that
the DAG wires up. Resumability comes for free. Per-rule conda environments scope a tool upgrade to
one rule. Moving the expensive rules to a cluster later is a profile change, since resources are
already declared per rule.

**Negative.** A second language in the repo, and `.smk` files are harder to unit-test than Python.
The mitigation is the rule that rules hold no logic: anything worth testing lives in
`src/virotrace/` and is reachable from the CLI, so the untestable surface stays thin.

**Enforcement.** A `run:` block longer than three lines becomes a function in the package. Reviewed
on sight.
