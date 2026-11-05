"""Configuration loading and validation.

The single point at which YAML becomes typed objects. Nothing else in the package
reads a config file, and no module holds a hard-coded analytical parameter.

Layering (later overrides earlier, key by key):
    config/config.yaml  ->  config/profiles/<profile>.yaml  ->  CLI overrides
"""

from __future__ import annotations

import copy
import hashlib
import json
import os
from pathlib import Path
from typing import Any, Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from .errors import ConfigError

# --------------------------------------------------------------------------- models


class _Base(BaseModel):
    model_config = ConfigDict(extra="allow", validate_assignment=True)


class GeneRanges(_Base):
    """HXB2 nucleotide coordinates, inclusive, 1-based."""

    pol: tuple[int, int]
    pr: tuple[int, int]
    rt: tuple[int, int]
    in_: tuple[int, int] = Field(alias="in")

    model_config = ConfigDict(extra="allow", populate_by_name=True)

    @model_validator(mode="after")
    def _ranges_are_codon_sized(self) -> GeneRanges:
        for name in ("pol", "pr", "rt", "in_"):
            start, end = getattr(self, name)
            if end <= start:
                raise ConfigError(f"reference.genes.{name}: end must exceed start")
            if (end - start + 1) % 3 != 0:
                raise ConfigError(
                    f"reference.genes.{name}: span {end - start + 1} is not a multiple of 3"
                )
        return self


class ReferenceConfig(_Base):
    accession: str = "K03455"
    record: str = "hxb2.gb"
    genes: GeneRanges


class HypermutConfig(_Base):
    enabled: bool = True
    p_threshold: float = Field(0.05, gt=0, lt=1)
    action: Literal["exclude", "flag_only"] = "exclude"


class QCConfig(_Base):
    min_gene_coverage: float = Field(0.80, ge=0, le=1)
    max_ambiguity_frac: float = Field(0.05, ge=0, le=1)
    exclude_internal_stop: bool = True
    exclude_frameshift: bool = True
    hypermutation: HypermutConfig = HypermutConfig()
    # Minimum nucleotide identity to HXB2 for a sequence to count as pol at all. Real pol
    # sits at roughly 80-95% (even for divergent subtypes); a sequence an aligner forced
    # against the wrong gene sits near 40-60%. Chosen from the observed histogram.
    min_identity_to_reference: float = Field(0.75, gt=0, le=1)
    # Share of a drug class's SDRM positions a sequence must cover to enter that class's
    # denominator. 1.0 = every position: a sequence that does not cover a position says
    # nothing about it, and counting it as "no resistance" would understate prevalence.
    min_class_position_coverage: float = Field(1.0, gt=0, le=1)
    exclude_lab_strains: bool = True
    lab_strain_list: str = "lab_strains.txt"
    # Clone-derived records (a /clone tag, "clone" or "proviral" in the definition) are
    # mostly several PCR clones from one patient, so counting each as a person would let
    # a few heavily cloned patients dominate prevalence. "exclude" removes them from the
    # prevalence and tree layers and keeps them in the DRM table; "keep_all" disables it.
    clone_policy: Literal["exclude", "keep_all"] = "exclude"
    # Several sequences from one person (amplicons, visits) are not independent
    # observations. "one_per_individual" keeps the earliest, least ambiguous sequence of each
    # person the patterns can identify; records they cannot identify count as one person each.
    individual_policy: Literal["keep_all", "one_per_individual"] = "one_per_individual"
    # Regexes applied to the isolate name; group 1 is the individual's id, e.g. RV217 isolates
    # are RV217_<participant>_v<visit>_...  Add a pattern per study whose naming is known.
    individual_id_patterns: list[str] = [r"^(RV217_\d+)_"]
    max_record_error_rate: float = Field(0.25, ge=0, le=1)

    @field_validator("individual_id_patterns")
    @classmethod
    def _patterns_compile_with_one_group(cls, patterns: list[str]) -> list[str]:
        import re

        for pattern in patterns:
            try:
                compiled = re.compile(pattern)
            except re.error as exc:
                raise ValueError(
                    f"individual_id_patterns: {pattern!r} is not a regex: {exc}"
                ) from exc
            if compiled.groups != 1:
                raise ValueError(
                    f"individual_id_patterns: {pattern!r} must have exactly one capture group"
                )
        return patterns


class SubsampleConfig(_Base):
    enabled: bool = True
    strata: list[str] = ["country_code", "year_bin", "subtype_called"]
    year_bin_width: int = Field(2, ge=1)
    max_per_stratum: int = Field(40, ge=1)
    min_per_stratum: int = Field(1, ge=0)
    priority: list[str] = ["date_precision", "gene_coverage", "ambiguity_frac"]

    @model_validator(mode="after")
    def _min_not_above_max(self) -> SubsampleConfig:
        if self.min_per_stratum > self.max_per_stratum:
            raise ConfigError(
                f"subsample.min_per_stratum ({self.min_per_stratum}) exceeds "
                f"max_per_stratum ({self.max_per_stratum})"
            )
        return self


class MixtureConfig(_Base):
    policy: Literal["any", "majority", "strict"] = "any"
    # A codon that could encode more amino acids than this carries no usable information
    # (NNN is the extreme) and is a no-call rather than a mixture.
    max_ambiguous_residues: int = Field(4, ge=1, le=20)


class SDRMConfig(_Base):
    list: str = "sdrm_who.tsv"
    list_version: str = "2009"
    #: Optional second list reported as a sensitivity analysis, never as the headline.
    sensitivity_list: str | None = None
    sensitivity_list_version: str | None = None


class DRMConfig(_Base):
    algorithm: str = "hivdb"
    algorithm_version: str | None = None
    classes: list[str] = ["NRTI", "NNRTI", "PI", "INSTI"]
    require_gene_coverage: bool = True
    mixtures: MixtureConfig = MixtureConfig()
    sdrm: SDRMConfig = SDRMConfig()


class PrevalenceConfig(_Base):
    strata: list[str] = ["year_bin", "drug_class"]
    also_by_subtype: bool = True
    ci_method: Literal["wilson"] = "wilson"
    ci_level: float = Field(0.95, gt=0, lt=1)
    min_stratum_n: int = Field(10, ge=1)
    trend_test: Literal["cochran_armitage", "none"] = "cochran_armitage"
    year_bin_width: int = Field(5, ge=1, le=20)
    # Hand-curated study_key -> treatment context, from reading the papers. Overrides the
    # keyword hint where a study is listed. Empty or absent file = hints only.
    study_annotations: str = "study_annotations.tsv"


class TemporalGateConfig(_Base):
    enabled: bool = True
    min_temporal_r2: float = Field(0.10, ge=0, le=1)
    require_positive_slope: bool = True
    date_randomisation: bool = False


class Admin1Config(_Base):
    enabled: Literal["auto", True, False] = "auto"
    min_resolved_frac: float = Field(0.30, ge=0, le=1)


class TraitsConfig(_Base):
    columns: list[str] = ["country_code", "region"]
    admin1: Admin1Config = Admin1Config()


class PhyloConfig(_Base):
    enabled: bool = True
    temporal_gate: TemporalGateConfig = TemporalGateConfig()
    traits: TraitsConfig = TraitsConfig()


class ReportConfig(_Base):
    min_stratum_n: int = Field(10, ge=1)
    include_limitations: bool = True
    figures: list[str] = []

    @model_validator(mode="after")
    def _limitations_are_mandatory(self) -> ReportConfig:
        # Not configurable by design. See docs/06-configuration.md section 6.
        if not self.include_limitations:
            raise ConfigError(
                "report.include_limitations cannot be false -- a surveillance report "
                "without its limitations section is not publishable"
            )
        return self


class SubtypeConfig(_Base):
    methods: list[str] = ["placement"]
    reference_panel: str = "subtype_references.fasta"
    min_support: float = Field(0.70, ge=0, le=1)  # share of site bootstraps agreeing with the call
    n_bootstrap: int = Field(100, ge=0, le=1000)
    # Panel labels that are not subtyping targets (the LANL set carries chimpanzee, gorilla
    # and HIV-1 groups N/O/P). CRFs are kept unless include_crfs lists the ones to keep.
    exclude_labels: list[str] = ["CPZ", "GOR", "O", "N", "P"]
    include_crfs: list[str] | None = None
    annotate_crfs: bool = True  # second placement with CRFs; annotation only
    min_overlap_sites: int = Field(300, ge=1)  # comparable sites needed for any call
    on_disagreement: Literal["report"] = "report"


class RecombinationConfig(_Base):
    enabled: bool = True
    method: Literal["window_scan"] = "window_scan"
    window: int = Field(300, ge=50)
    step: int = Field(100, ge=10)
    min_window_sites: int = Field(150, ge=20)
    min_margin: float = Field(0.02, ge=0)  # nearest vs next label, in p-distance
    min_segment: int = Field(400, ge=50)  # columns a same-label run must span
    confirm_support: float = Field(0.90, ge=0, le=1)
    confirm_bootstrap: int = Field(100, ge=10, le=1000)
    tree_policy: Literal["exclude_confirmed", "exclude_suspected", "keep_all"] = "exclude_confirmed"
    drm_policy: Literal["keep_all"] = "keep_all"


class Settings(_Base):
    """Fully resolved, validated configuration."""

    version: int
    schema_version: int
    scope: str
    seed: int
    threads: int = Field(4, ge=1)
    snapshot_id: str | None = None

    paths: dict[str, str]
    scopes: dict[str, dict[str, list[str]]]
    reference: ReferenceConfig
    qc: QCConfig
    subsample: SubsampleConfig
    drm: DRMConfig
    prevalence: PrevalenceConfig
    subtype: SubtypeConfig = SubtypeConfig()
    recombination: RecombinationConfig = RecombinationConfig()
    phylo: PhyloConfig
    report: ReportConfig

    # ----------------------------------------------------------------- cross-field
    @model_validator(mode="after")
    def _cross_field_rules(self) -> Settings:
        if self.scope not in self.scopes:
            raise ConfigError(
                f"scope '{self.scope}' is not defined in scopes: {sorted(self.scopes)}"
            )
        if self.drm.require_gene_coverage and self.qc.min_gene_coverage <= 0:
            raise ConfigError(
                "drm.require_gene_coverage is true but qc.min_gene_coverage is 0 -- "
                "per-class denominators would be meaningless"
            )
        return self

    # ----------------------------------------------------------------- derived
    def stage_seed(self, stage: str) -> int:
        """Deterministic per-stage seed.

        Derived from (root seed, stage name) so that adding or reordering a stage
        does not shift another stage's random stream.
        """
        digest = hashlib.sha256(f"{self.seed}:{stage}".encode()).digest()
        return int.from_bytes(digest[:4], "big")

    def config_hash(self) -> str:
        payload = json.dumps(
            self.model_dump(mode="json", by_alias=True), sort_keys=True, separators=(",", ":")
        )
        return hashlib.sha256(payload.encode()).hexdigest()[:12]

    def warnings(self) -> list[str]:
        """Non-fatal configuration smells worth printing loudly."""
        out: list[str] = []
        if self.phylo.traits.admin1.enabled is True:
            out.append(
                "phylo.traits.admin1.enabled is forced true, bypassing the resolved-fraction "
                "gate. Subnational reconstruction on country-only metadata overstates "
                "geographic resolution."
            )
        if self.report.min_stratum_n < 5:
            out.append(
                f"report.min_stratum_n is {self.report.min_stratum_n}; small denominators "
                "will be plotted as if they were estimates."
            )
        if not self.qc.hypermutation.enabled:
            out.append(
                "qc.hypermutation.enabled is false. APOBEC-edited sequences produce "
                "spurious resistance signal; prevalence estimates will be unreliable."
            )
        return out


# --------------------------------------------------------------------------- loading


def deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    """Recursively merge ``override`` into ``base``, returning a new dict."""
    out = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = deep_merge(out[key], value)
        else:
            out[key] = copy.deepcopy(value)
    return out


def _read_yaml(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise ConfigError(f"config file not found: {path}")
    with path.open() as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ConfigError(f"config file is not a mapping: {path}")
    return data


def load_settings(
    config_path: Path | str = "config/config.yaml",
    profile: str | None = None,
    overrides: dict[str, Any] | None = None,
) -> Settings:
    """Load, layer, and validate configuration.

    Layer order, later wins: base file < profile < ``VIROTRACE_SCOPE`` < ``overrides``.
    ``VIROTRACE_PROFILE`` and ``VIROTRACE_SCOPE`` are set by the Snakefile so every
    stage it launches resolves exactly the configuration Snakemake resolved.
    """
    config_path = Path(config_path)
    merged = _read_yaml(config_path)

    profile = profile or os.environ.get("VIROTRACE_PROFILE") or None
    if profile:
        profile_path = config_path.parent / "profiles" / f"{profile}.yaml"
        merged = deep_merge(merged, _read_yaml(profile_path))

    env_scope = os.environ.get("VIROTRACE_SCOPE")
    if env_scope:
        merged = deep_merge(merged, {"scope": env_scope})

    if overrides:
        merged = deep_merge(merged, overrides)

    try:
        return Settings.model_validate(merged)
    except ConfigError:
        raise
    except Exception as exc:  # pydantic ValidationError and friends
        raise ConfigError(f"configuration failed validation:\n{exc}") from exc
