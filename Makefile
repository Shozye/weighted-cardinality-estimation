PACKAGE_NAME = weighted_cardinality_estimation
SRC_DIR = src
BUILD_DIR = build

.PHONY: all help build build_fast test

# Runs with fast compile with -O3 optimisation flag
build_fast:
	pip install . \
		--no-deps \
		--no-build-isolation \
		-Cbuild-dir=$(BUILD_DIR) \
		-Ccmake.args="-DCMAKE_BUILD_TYPE=Release" \
		-vvv

# it will also install dependencies so good for first time
build_with_deps:
	python -m pip install . -Cbuild-dir=$(BUILD_DIR) -vvv

# Runs all unit tests
test:
	python -m pytest tests/

# this command will run regression only for the newest pushed commit and NOT SAVE IT TO DISK
asv_debug:
	asv run HEAD^! --dry-run --durations --show-stderr

# This command is used to run regression for asv. It runs it since stable commit. 
asv_regression:
	asv run NEW

clean: clean_cpp 

clean_all: clean_asv clean_cpp clean_mypy

clean_asv:
	rm -rf $(PACKAGE_NAME)
	rm -rf env

clean_cpp:
	rm -rf ./.cache
	rm -rf ./build

clean_mypy:
	rm -rf ./.mypy_cache
