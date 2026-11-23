#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p data results
shasum -a 256 -c PREREGISTRATION.sha256
export UV_CACHE_DIR="$PWD/.uv-cache"
if [ ! -x .venv/bin/python ]; then
    uv venv --python 3.12.7 .venv
fi
uv pip install --python .venv/bin/python 'numpy==2.5.5' 'scipy==1.18.1' 'gmpy2==2.3.2' 'matplotlib==3.11.2'
gmp_prefix="${GMP_PREFIX:-/opt/homebrew/opt/gmp}"
compiler="${CXX:-clang++}"
for name in sample gates; do
    "$compiler" -O3 -std=c++17 -Wno-deprecated-literal-operator -I"$gmp_prefix/include" "$name.cpp" -L"$gmp_prefix/lib" -lgmpxx -lgmp -o "$name"
done
"$compiler" -O3 -std=c++17 -shared -fPIC score.cpp -o score.dylib
.venv/bin/python grammar.py > catalog_counts.json
.venv/bin/python generate.py > data/generation.log
.venv/bin/python features.py > data/features.log
.venv/bin/python scoring.py > results/scoring.log
.venv/bin/python statistics_run.py > results/statistics.log
.venv/bin/python scaling.py > results/scaling.log
.venv/bin/python controls.py > results/controls.log
.venv/bin/python diagnostics.py > results/diagnostics.log
.venv/bin/python audit.py > results/audit.log
.venv/bin/python distinct.py > results/distinct.log
.venv/bin/python deliver.py > results/export.log
.venv/bin/python population24.py > results/population24.log
.venv/bin/python figure.py > results/figure.log
.venv/bin/python residue_null.py > results/residue_null.log
cat results/scaling.log results/audit.log results/distinct.log
