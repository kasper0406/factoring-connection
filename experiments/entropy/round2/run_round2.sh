#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
p=../.venv/bin/python
export VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1
(cd .. && shasum -a 256 -c round2/prereg.sha256)
shasum -a 256 -c power_addendum.sha256
./build.sh
"$p" make_data.py
"$p" - <<'PY'
from config import *
run(['./sampler',64,0,1024,'cost_panel_v1',DATA/'cost.rows'])
for n in (40,48,56,64,80,96):run(['./sampler',n,0,86,'witness_v1',DATA/f'wit{n}.rows'])
PY
"$p" grammar.py --frozen
shasum -a 256 -c grammar.sha256
"$p" large_data.py
"$p" budget.py
"$p" labels.py
"$p" extract.py DV
"$p" screen.py
shasum -a 256 -c locks.sha256
"$p" confirm_eval.py
"$p" aggregate.py
"$p" supplement.py
"$p" -c 'from aggregate import aggregate; aggregate(64,"S"); aggregate(96,"S")'
"$p" inference.py original
"$p" inference.py augmented
"$p" validation_audit.py
"$p" identity.py
"$p" dyad.py
"$p" decay.py
"$p" conditional.py
"$p" sensitivity.py
"$p" audit_arithmetic.py
"$p" catalog_audit.py
"$p" coverage.py
"$p" certificates.py
"$p" verify_results.py
"$p" mi_null.py
"$p" matched_power.py
"$p" data_manifest.py
"$p" export.py
"$p" figure.py
