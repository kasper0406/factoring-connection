# NFS smoothness experiment

The two preregistrations, exact candidate arrays, fitted models, held-out scores, timing records and report live in this directory. The original cohort uses B=4093. The independent power cohort lives in `power_cohort/` and uses B=262139. Both use y=31.

## Files

| Path relative to each cohort | Contents |
|---|---|
| `PREREGISTRATION.md`, `preregistered.utc` | Frozen protocol and UTC timestamp. SHA-256 records are in the project. |
| `data/manifest.json` | Modulus factors, polynomial coefficients in increasing degree, root m, split, bit size, record counts, file checksum and generation timing. |
| `data/nNN.npz` | Integer rows in generation order. The first `reference_count` are the prior history. Remaining rows are targets. |
| `data/edgesNN.npz` | For each target, the 16 reference indices used. Reference indices are local to that N. |
| `results/dataset_counts.csv`, `results/integer_audit.json` | Population by N, region, stream and label, plus arithmetic and history audit counts. |
| `results/{true,shuffled,planted}/selection_lock.json` | Selection scores, calibration weights, frozen model SHA-256 and chosen variant, written before final scoring. |
| `results/{true,shuffled,planted}/models.joblib` | Fitted size regressions, calibration parameters, all six frozen flexible learners. |
| `results/{true,shuffled,planted}/heldout_scores.npz` | Full precision held-out predictions and target labels, aligned in manifest order and then target row order within N. Columns are joint, rational and algebraic. |
| `results/{true,shuffled,planted}/heldout_summary.json` | Paired family CIs, exact sign-flip probabilities, restricted permutation results and verdict. The permutation NPZs store sampled row indices and every permuted loss. |
| `results/{true,shuffled,planted}/by_modulus_geometry.csv` | Log-loss, AUC, label counts and yield at a fixed confirmation budget, for all models and all targets. |
| `results/family_metrics.csv` | Absolute metrics and paired gains versus logistic, matched learner, and geometry-only where relevant, by family and region. Each eligible modulus has equal weight. |
| `results/pooled_yield.csv` | Counts and ratios aggregated from each run's `yield_ledger.csv`. Bin -1 uses the overall 20% budget. Bins 0–3 independently rank inside each registered u interval. |
| `results/true/unblended_joint.csv` | Descriptive sensitivity using the same trained models and sigmoid calibration, with the convex baseline blend disabled. |
| `results/true/work_seconds.csv`, `work_summary.json` in the same directory | Replayed operating costs and exact confirmed yields by policy. |
| `logs/` | Training and testing output, arithmetic/environment information, and timing records. |

The original cohort also has `CONTROL_SUPPLEMENT.md` and `results/planted_strong/`. The high-bound cohort has `HISTORY_AMENDMENT.md`, with the failed 1024-reference geometry check retained in its logs. Its supplementary global negative control lives in `results/global_shuffle/`, and `results/shuffled/within_cell_placebo.json` records the within-cell label randomization.

The `coords` columns in each data NPZ are `a, b, i, geometry`. Geometry 0 denotes the full box and 1 the special-q lattice. The `valuation` array is indexed by row, rational/algebraic side and ascending prime through 31. The `cofactor_words` array contains unsigned little-endian 64-bit limbs for each absolute stripped cofactor, so `C = int(low) + (int(high) << 64)`. The absolute norm is C times the recorded small prime powers. The stored `sides` columns are rational and algebraic smoothness. For real labels their product is the joint target. Synthetic joint labels are a separate prediction target. Their side labels support the multitask learner.

## Reproduction

Run from this directory in the supplied workspace. The repository Python supplies GMP, PARI and SymPy. The local ML environment is pinned in `requirements-ml.txt`. The replay uses the same two Python environments. Commands below recreate outputs in place.

```sh
ARITH_PY=/Volumes/git/factoring-connection/.venv/bin/python
PROJECT_DIR="$PWD"
ML_PY="$PROJECT_DIR/.venv/bin/python"
export UV_CACHE_DIR="$PROJECT_DIR/.uv-cache"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
[ -x "$ML_PY" ] || uv venv --python "$ARITH_PY" .venv
uv pip install --python "$ML_PY" -r requirements-ml.txt

"$ARITH_PY" generate.py
"$ML_PY" fit_models.py
"$ML_PY" evaluate.py
"$ML_PY" strong_control.py
"$ARITH_PY" audit.py
"$ML_PY" diagnostics.py
"$ML_PY" work_ledger.py
"$ML_PY" verify_results.py

(
  cd power_cohort
  export GUIDE_STUDY=power
  "$ARITH_PY" ../generate.py
  "$ML_PY" ../fit_models.py
  "$ML_PY" ../evaluate.py
  "$ARITH_PY" ../audit.py
  "$ML_PY" ../diagnostics.py
  "$ML_PY" ../work_ledger.py
  "$ML_PY" ../null_diagnostic.py
  "$ML_PY" ../verify_results.py
)
"$ML_PY" summarize.py
```

The model fit uses one BLAS thread and two OpenMP threads. The work replay uses one thread. Its summed phase timings cover both norms and local stripping on every target, score construction, stable ranking, and short-circuit exact PARI confirmation of retained candidates. The comparisons amortize polynomial selection and training across future queries. Times are specific to this Python replay. The source predictor can skip all flexible work when calibration assigns it zero weight.
