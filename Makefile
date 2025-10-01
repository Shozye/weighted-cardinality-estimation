.PHONY: build test bench

# i run it in my venv
build:
	python -m pip install -e . -Cbuild-dir=.build -vvv
# -vvv is to throw errors --WExtra --WError during C++ compilation
# --no-deps and --no-build-isolation are to make build faster

build_clean:
	rm -rf .build
	python -m pip install . --no-deps --no-build-isolation -Cbuild-dir=.build -vvv

build_fast:
	python -m pip install . --no-deps --no-build-isolation -Cbuild-dir=.build -vvv

test:
	pytest tests/

# on my PC (ryzen 7 5800X), it runs 83 seconds.
bench:
	pytest benchmarks/bench_sketches.py -q --benchmark-disable-gc --benchmark-warmup=on

# This command is used to run regression for asv. It runs it since stable commit. 
asv_regression:
	asv run NEW