PACKAGE_NAME = weighted_cardinality_estimation
SRC_DIR = src
BUILD_DIR = build
PIP_WHEELS_DIR := $(HOME)/pip-wheels
VENV := .venv
PYTHON := $(VENV)/bin/python
PIP := $(VENV)/bin/pip
ASV := $(VENV)/bin/asv

.PHONY: asv_publish download-wheels

# Download all deps into a local wheel cache (for offline installs)
download-wheels:
	pip download . scikit-build-core pybind11 pytest -d $(PIP_WHEELS_DIR)

# this command will run regression only for the newest pushed commit and NOT SAVE IT TO DISK
asv_debug:
	$(ASV) run HEAD^! --dry-run --durations --show-stderr

# This command is used to run regression for asv. It runs it since stable commit.
asv_regression:
	$(ASV) run NEW

asv_publish:
	$(ASV) publish
	$(ASV) gh-pages --rewrite

clean: clean_cpp

clean_all: clean_asv clean_cpp clean_mypy

clean_asv:
	rm -rf $(PACKAGE_NAME)
	rm -rf env

clean_cpp:
	rm -rf ./.cache
	rm -rf ./build
	rm -rf $(VENV)

clean_mypy:
	rm -rf ./.mypy_cache
