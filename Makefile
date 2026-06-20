PACKAGE_NAME = weighted_cardinality_estimation
SRC_DIR = src
BUILD_DIR = build
PIP_WHEELS_DIR := $(HOME)/pip-wheels
VENV := .venv
PYTHON := $(VENV)/bin/python
PIP := $(VENV)/bin/pip
ASV := $(VENV)/bin/asv

.PHONY: asv_publish download-wheels venv library-build

# Create venv and install build tools only (no WCE library)
venv:
	uv venv --clear $(VENV)
	uv pip install --python $(VENV) --find-links $(PIP_WHEELS_DIR) --no-index scikit-build-core pybind11 ninja
	uv pip install --python $(VENV) --find-links $(PIP_WHEELS_DIR) --no-index "."

# Build WCE wheel into dist/ and install into its own venv
library-build:
	uv cache clean weighted-cardinality-estimation
	rm -rf dist
	uv build --no-build-isolation \
		-Cbuild-dir=build \
		-Ccmake.args="-DCMAKE_BUILD_TYPE=Release" \
		--wheel .
	uv pip install --python $(VENV) --no-index --reinstall-package weighted-cardinality-estimation dist/*.whl

# Download all deps into a local wheel cache (for offline installs)
download-wheels:
	uv run --no-project --with pip -- pip download . -d $(PIP_WHEELS_DIR)

# this command will run regression only for the newest pushed commit and NOT SAVE IT TO DISK
asv_debug:
	$(ASV) run HEAD^! --dry-run --durations --show-stderr

# This command is used to run regression for asv. It runs it since stable commit.
asv_regression:
	$(ASV) run NEW

asv_publish:
	$(ASV) publish
	$(ASV) gh-pages --rewrite

clean:
	rm -rf ./.cache
	rm -rf ./build
	rm -rf $(VENV)
	rm -rf ./dist
	rm -rf ./*.egg-info
	rm -rf ./.pytest_cache
	rm -rf ./.mypy_cache
	rm -rf ./$(PACKAGE_NAME)
	rm -rf ./env
	find . -name __pycache__ -type d | xargs rm -rf
	rm -rf weighted-cardinality-estimation
	rm -rf asv_results
	rm -rf .benchmarks
