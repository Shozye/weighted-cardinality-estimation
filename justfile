venv_dir := ".venv"

[private]
default:
    just --list

# Create venv and install build tools only (no WCE library)
venv:
    uv venv {{venv_dir}}
    uv pip install --python {{venv_dir}} --find-links ~/pip-wheels --no-index scikit-build-core pybind11 pytest

# Build WCE wheel into dist/ and install into its own venv
library-build:
    uv venv {{venv_dir}}
    uv pip install --python {{venv_dir}} --find-links ~/pip-wheels --no-index scikit-build-core pybind11 pytest
    uv cache clean weighted-cardinality-estimation
    rm -rf dist
    uv build --no-build-isolation \
        -Cbuild-dir=build \
        -Ccmake.args="-DCMAKE_BUILD_TYPE=Release" \
        --wheel .
    uv pip uninstall --python {{venv_dir}} -y weighted-cardinality-estimation 2>/dev/null || true
    uv pip install --python {{venv_dir}} --no-index dist/*.whl

# Run unit tests using WCE venv
library-test:
    {{venv_dir}}/bin/pytest tests/

# Remove build artifacts and venv
clean:
    rm -rf build {{venv_dir}}
