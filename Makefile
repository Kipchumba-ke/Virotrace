.DEFAULT_GOAL := help
SHELL := /bin/bash
PROFILE ?= dev
CORES ?= 4

.PHONY: help setup lint type test test-all smoke sdrm dag dev kenya africa refresh-golden clean clean-all

help: ## Show this help
	@grep -hE '^[a-z-]+:.*?## ' $(MAKEFILE_LIST) | awk 'BEGIN{FS=":.*?## "}{printf "  \033[36m%-16s\033[0m %s\n", $$1, $$2}'

setup: ## Create the conda environment and install the package
	mamba env create -f environment.yml || conda env create -f environment.yml
	@echo "now run: conda activate virotrace"

lint: ## Lint and format check
	ruff check src tests
	ruff format --check src tests

type: ## Static type check
	mypy

test: ## Fast suite (unit, contract, integration) -- target under 60s
	pytest tests/unit tests/contract tests/integration

test-all: ## Everything including the smoke test
	pytest

smoke: ## Full dev-profile DAG end to end
	pytest -m smoke

sdrm: ## Rebuild the SDRM lists from config/reference/sources/ (checks wild types against HXB2)
	python scripts/sdrm_from_text.py config/reference/sources/sdrm_who_2009.txt --version 2009 \
		--source "Bennett et al. 2009, PLoS ONE 4(3):e4724; transcribed from Stanford HIVdb SDRM page (2026-09-25), original WHO 2009 entries" \
		--out config/reference/sdrm_who.tsv
	python scripts/sdrm_from_text.py config/reference/sources/sdrm_stanford_2026-09-25.txt --version stanford-2026-09-25 \
		--source "Stanford HIVdb SDRM page last updated 2026-09-25 (WHO 2009 + Rhee 2021/2023 + INSTI Tzou 2020/Tao 2023)" \
		--out config/reference/sdrm_stanford_2026-09-25.tsv

dag: ## Validate the DAG for every profile without executing anything
	virotrace validate config --config config/config.yaml --profile dev
	virotrace validate config --config config/config.yaml --profile kenya
	virotrace validate config --config config/config.yaml --profile africa
	pytest tests/smoke/test_dev_build.py -k dag -q

dev: ## Run the toy build (minutes, no network)
	snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 2 \
		--config scope=kenya profile=dev

kenya: ## Run the primary Kenya build
	snakemake -s workflow/Snakefile --configfile config/config.yaml --cores $(CORES) \
		--config scope=kenya profile=kenya

africa: ## Run the Africa context build (long; resumable)
	snakemake -s workflow/Snakefile --configfile config/config.yaml --cores $(CORES) \
		--config scope=africa profile=africa --rerun-incomplete

refresh-golden: ## Regenerate regression baselines -- the diff must be explained in the PR
	pytest tests/regression --snapshot-update || true
	@echo "review the diff before committing; an unexplained refresh blesses a regression"

clean: ## Remove derived results, keep raw snapshots
	rm -rf data/interim data/processed results auspice reports logs .snakemake

clean-all: clean ## Also remove raw snapshots
	rm -rf data/raw
