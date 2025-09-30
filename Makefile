.PHONY: build test bench

# i run it in my venv
build:
	python -m pip install -e . -Cbuild-dir=.build -vvv
# -vvv is to throw errors --WExtra --WError during C++ compilation
# --no-deps and --no-build-isolation are to make build faster

build-clean:
	rm -rf _build
	python -m pip install -e . --no-build-isolation -Cbuild-dir=.build -vvv

build_fast:
	python -m pip install -e . --no-build-isolation -Cbuild-dir=.build -vvv

test:
	pytest tests/

# on my PC (ryzen 7 5800X), it runs 83 seconds.
bench:
	pytest benchmarks/bench_sketches.py -q --benchmark-disable-gc --benchmark-warmup=on
