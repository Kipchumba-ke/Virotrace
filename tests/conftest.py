"""Shared fixtures.

All sequence fixtures are synthetic or public reference material. No real
patient-derived sequence data lives in this repository.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from virotrace.drm.sdrm import load_sdrm_list

REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture(scope="session")
def repo_root() -> Path:
    return REPO_ROOT


@pytest.fixture(scope="session")
def config_path() -> Path:
    return REPO_ROOT / "config" / "config.yaml"


@pytest.fixture(scope="session")
def sdrm_list():
    return load_sdrm_list(
        REPO_ROOT / "config" / "reference" / "sdrm_who.PLACEHOLDER.tsv", version="placeholder"
    )


@pytest.fixture
def toy_records() -> list[dict]:
    """Metadata rows exercising the schema edges that break naive code."""
    return [
        {
            "strain": "KE_001",
            "country_code": "KEN",
            "year": 2015,
            "subtype_called": "A1",
            "date_precision": "day",
            "gene_coverage": 0.98,
            "ambiguity_frac": 0.001,
        },
        {
            "strain": "KE_002",
            "country_code": "KEN",
            "year": 2015,
            "subtype_called": "A1",
            "date_precision": "year",
            "gene_coverage": 0.91,
            "ambiguity_frac": 0.02,
        },
        {
            "strain": "KE_003",
            "country_code": "KEN",
            "year": 2016,
            "subtype_called": "A1",
            "date_precision": "month",
            "gene_coverage": 0.95,
            "ambiguity_frac": 0.004,
        },
        {
            "strain": "KE_004",
            "country_code": "KEN",
            "year": 2016,
            "subtype_called": "D",
            "date_precision": "day",
            "gene_coverage": 0.88,
            "ambiguity_frac": 0.01,
        },
        {
            "strain": "UG_001",
            "country_code": "UGA",
            "year": 2016,
            "subtype_called": "A1",
            "date_precision": "day",
            "gene_coverage": 0.99,
            "ambiguity_frac": 0.0,
        },
        {
            "strain": "KE_005",
            "country_code": "KEN",
            "year": 2015,
            "subtype_called": "A1",
            "date_precision": "unknown",
            "gene_coverage": 0.83,
            "ambiguity_frac": 0.04,
        },
    ]
