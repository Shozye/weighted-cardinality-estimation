set dotenv-load

venv_dir := ".venv"

[private]
default:
    just --list

# Create venv and install build tools only (no WCE library)
venv:
    make venv

# Build WCE wheel into dist/ and install into its own venv
library-build:
    make library-build

# Run fast unit tests only (excludes benchmarks)
run-library-test-only-fast:
    {{venv_dir}}/bin/pytest tests/ -m "not benchmark"

# Run benchmarks only
run-library-test-only-benchmarks:
    {{venv_dir}}/bin/pytest tests/ -m "benchmark"

# Remove build artifacts and venv
clean:
    make clean

# Run all ASV benchmarks against HEAD (uses current venv, no rebuild)
asv-run:
    {{venv_dir}}/bin/asv run --python same HEAD^!

# Compare two commits, e.g.: just asv-compare HEAD~1 HEAD
asv-compare from to:
    {{venv_dir}}/bin/asv compare {{from}} {{to}}

# Build HTML report and open it
asv-html:
    {{venv_dir}}/bin/asv publish
    {{venv_dir}}/bin/asv preview
