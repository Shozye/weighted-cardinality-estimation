build:
	python -m pip install -e . --no-deps --no-build-isolation -vvv
# -vvv is to throw errors --WExtra --WError during C++ compilation
# --no-deps and --no-build-isolation are to make build faster

test:
	pytest tests/

bench:
	pytest benchmarks/bench_sketches.py -q --benchmark-disable-gc --benchmark-warmup=on
