"""Virotrace: HIV-1 molecular surveillance for Kenya and Africa.

Two analysis layers over one curated sequence corpus:
  - drm/   genotypic drug-resistance calling and prevalence estimation
  - phylo/ maximum-likelihood and time-resolved phylogenetics with trait reconstruction

NOT a clinical tool. See README.md and docs/01-technical-specification.md.
"""

__version__ = "0.1.0"
SCHEMA_VERSION = 1

__all__ = ["SCHEMA_VERSION", "__version__"]
