# Reproducing the second smoothness experiment

Run from this directory in the shared workspace. The round-one Python environment supplies NumPy, SciPy and scikit-learn. The arithmetic scripts import GMP, PARI and SymPy from the read-only repository environment. The additional full-column comparator uses LightGBM 4.6.0 installed locally in `.extra`. The frozen design is in [PREREGISTRATION.md](PREREGISTRATION.md), including the two cutoffs and conditions for a discovery. The files named `*SUPPLEMENT.md`, `PC2B_DESIGN.md`, `REGION_AUDIT_NOTE.md`, and `ALL_LOCAL_CHECK.md` describe subsequent diagnostics.

## Primary files

| File | Contents |
|---|---|
| `data/manifest.json`, `data/nNN.json` | 48 synthetic semiprimes, prime factors, polynomials, modular root, split, matched parent pairs, ideal roots, reduced bases, match costs and allowed maps/shifts. |
| `data/nNN.npz` | All target and donor probes for one N. Columns of `coords` are a, b, and the removed algebraic special-q (1 in the ordinary box). The first `nt` rows are candidate targets before deduplication. |
| `data/keepNN.npy` | Global target ownership mask. A repeated physical coordinate has a single owner in phase and modulus order, determined before model fitting. |
| `data/gNN.npz` | Coordinate arithmetic, continued fractions, embedding/covariant, raw bit and coset features. |
| `models/` | Cross-fit baselines, family models, plant models, precomputed node offsets, and extra all-column trees. Binary model files are implementation artifacts. Import the defining classes in the loading Python process for objects saved by a script run as main. |
| `results/validation_lock.json`, `lock_provenance.json` | Frozen validation screen, time and preregistration digest. |
| `results/family_summary.csv` | One row per family/cutoff combining physical joint gain, adjusted verdict, and native control sensitivity. |
| `results/A_families.csv` | Primary 108 contrasts. All 8 families plus their combination, both cutoffs and all 6 targets. Includes full-domain gain, uncertain-stratum gain, 99% clustered CI, information per hit, Holm p, narrow-cell tests and final verdict. |
| `results/A_parents.csv`, `match_balance.json` | Direct matched parent comparison, before and after B1, full and uncertain data, all binary/continuous sides. Match-balance columns are log q, log basis-length ratio, and rational and algebraic log-norm envelopes. |
| `results/A_baselines.csv`, `A_covariance.json`, `A_scaling.csv` | Baseline absolute log-loss, centered covariance of the side residuals, and empirical information by thirds of U. |
| `results/B_*.csv`, `Bnull_*.csv`, `crt_*.csv` | Buffered, within-modulus split, Bernoulli(B0) null and raw-bit randomized-CRT null. |
| `results/Afull_*.csv` | All-family robustness analysis with the unscreened prime-power/valuation boosted comparator. |
| `results/all_cell_permutations.csv` | Full cell permutations and exact cell-conditional logistic estimates supplementing the locked disjoint-pair screen. |
| `results/native_sensitivity.csv`, `native_gates_*.csv`, `native_dose_response.csv` | Native twofold plants, validation and microcell gates, and additional 4x/8x plants. Conditional power is across 100 independent test-label draws for fixed trained predictors. Bounds and multiplicity match the described control analysis. |
| `results/pc1_*.csv`, `pc2a_*.csv`, `pc2b*.csv`, `pc2b*match.json` | Gaussian tilt and recovery, forbidden-band leakage fixture, matched arithmetic controls and scaling. |
| `results/audit.json`, `verification.json`, `L1_positive_fixture.json` | Independent arithmetic/geometry checks, independent recomputation of all 108 main gains and Holm values, and exact parallel-shift detector check. |

The main NPZs store exact local valuations as sparse rows (`ptr,vi,ve`), indexed by side then prime through 32,768. The public coordinates determine every projective prime-power class. `sizes,local,w,flags` have shape (nodes, 2 cutoffs, 2 sides). Flags are +1 for forced smooth, −1 for forced rough, 0 for uncertain. `y` contains exact side labels, and `raw` the two normalized log norms. `digits` is available strictly for donor lookups. `leak` is reserved for the deliberate positive control. Signs and integer norms can be reconstructed directly. All permutations and matches use fixed seeds.

For each target the 12 edge slots refer to up to two safe shifts, one Farey neighbour, one admissible projective transform, two componentwise-product factors, two bilinear-sum/difference factors, and four prior distinct-q donors. −1 means unavailable. G4/G5 scores use the domain with an eligible relative, G6 uses special-q children, and G7 has a coset feature everywhere. A fixed projective transform sharing a global polynomial factor has zero resultant and is disallowed, including any exact label-preserving symmetry. Costs of individual probes are visible in the node and edge ledgers.

For per-geometry gains and AUC, see `results/by_regime.csv`. Archived A-split arrays may carry the initial unscaled-angle diagnostic tag; their membership and scores use only the modulus split. The final B-split files use the nominal box skew of 2.2.

The analytic denominator is the numerical Buchstab rough density. B0 approximates a size-local random-integer probability. The registered B1 sees all candidate local columns at screening, then fits depth-four boosting and a 30-knot cubic GAM with out-of-modulus offsets and isotonic calibration. A one-class binary training fold falls back to B0. The separate full-column fit offers every valuation and every projective prime-power level simultaneously. In the locked screen, shared random swaps act on disjoint pairs in common-support N/sign/0.02-nat cells. Each cutoff gets a step-down min-P validation budget of 0.025; the final 108-way Holm budget is 0.01. Full-cell permutations are provided separately.

## Commands

The commands regenerate outputs in place. All Python commands use one BLAS thread. The standalone control and robustness studies reuse the frozen primary design and new synthetic labels or new control moduli.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export UV_CACHE_DIR="$PWD/../.uv-cache"
PY=../.venv/bin/python
clang++ -std=c++17 -O3 -shared -fPIC local.cpp -o local.dylib
uv pip install --python "$PY" --target .extra --no-deps lightgbm==4.6.0
"$PY" rdata.py
"$PY" rfeatures.py
"$PY" raudit.py
for cut in 0 1; do
  "$PY" rbase.py "$cut"
  "$PY" rfit.py "$cut" A
done
"$PY" rstats.py screen
"$PY" rstats.py A
"$PY" rfullcells.py
sh region_run.sh
for cut in 0 1; do
  "$PY" rcontrols.py "$cut"
  "$PY" rplantgate.py "$cut"
  "$PY" rdose.py "$cut"
  "$PY" ralllocal.py "$cut"
  "$PY" rfit.py "$cut" Afull
done
"$PY" rcontrols.py aggregate
"$PY" rdose.py aggregate
"$PY" rstats.py Afull
"$PY" rbandcheck.py
"$PY" rpc2b.py
ROUND2_PC2_POWER=1 "$PY" rpc2b.py
"$PY" rbyregime.py
"$PY" rverify.py
"$PY" rsummary.py
```
