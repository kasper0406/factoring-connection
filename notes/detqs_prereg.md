# detqs pre-registration (frozen before any experiment run)

Agent `detqs`, phase 8, 2026-10-19. Task: polyimpl D1, deterministic fixed-u QS on almost all N.

Frozen code (this commit): `src/detqs/{sieve,linalg,qs,localmodel,harness,analysis,prereg_values}.py`,
tests `tests/test_detqs.py` (15 pass). No sieve run on any experiment instance happened before this commit.
The local model was evaluated beforehand on separately seeded N (`detqs|model|n`) only to size the
predictions. That is model output, not data.

## Algorithm under test (operational, deterministic)
- m = ⌈√N⌉ and Q_t = (m+t)² − N for t = 1..T.
- Factor base: primes ℓ ≤ y with x² ≡ N (mod ℓ) solvable.
- Smoothness is exact: every prime power is divided out.
- Relations are taken in increasing t and added to an incremental exact F₂ elimination.
- Each new kernel vector S is tested with gcd(X_S − Y_S, N). The run stops at the first proper gcd, or
  after DMAX = 30 kernel vectors.
- Parameters: y(n) = 2^⌊n/5⌋ and T(n) = 2^{16,18,18,19,20,21} for n = 40, 48, 56, 64, 72, 80.
- These are small-size parameters. The asymptotic regime y = N^c, c < 1/(6ω), is not reachable at
  40–80 bits. The experiments test *structure*: where the variance of Z_N comes from, and the H2 rank
  event. They do not test asymptotic exponents.

## Instances (seeded, `src/detqs/harness.py`)
- Integers: N uniform in [2^{n−1}, 2^n), seed `detqs|int|n|i`.
- Semiprimes: uniform on S_n by rectangle rejection, with seeds `detqs|sp|n|i` (E1) and `detqs|e2|n|i` (E2).

## E1: variance and lower tail of Z_N = #{t ≤ T : Q_t is y-smooth}
- Samples: n ∈ {40, 48, 56, 64} with 2000 integers and 1000 semiprimes each; n = 72 with 1000 and 500.
- Local model: L(N) = Σ_blocks w_b ∏_{ℓ≤y} λ_{ℓ,u_b}(N) / Σ w_b. Here λ is the exact Hensel distribution of
  v_ℓ(x² − N), weighted by ρ(u − e log ℓ/log y)/ρ(u), against the uniform Z_ℓ law. See
  `localmodel.py`.
- Definitions:
  - Z̄ = sample mean on integers.
  - V = Var(Z/Z̄).
  - V_L = Var(L/L̄).
  - c = ΣZ/ΣL and r_N = Z_N/(c L(N)).
  - V_res = Var(r).
  - V_strat = pooled within-stratum variance of Z/Z̄, with strata (N mod 8, N mod 3, N mod 5).
- Predictions, on integers at every n unless stated otherwise:
  - **P1a** (local excess): V ≥ 10/Z̄.
  - **P1b** (model calibration): 0.4 V_L ≤ V ≤ 2.5 V_L.
  - **P1c**: corr(Z, L) ≥ 0.8.
  - **P1d** (model explains the variance): V_res ≤ 0.25 V.
  - **P1e** (model-free stratification): V_strat ≤ 0.6 V.
  - **P1f** (no lower-tail outliers beyond the local model): min r_N ≥ 0.5, on integers and on semiprimes.
  - **P1g** (H1 on S_n at the frozen T): every sampled semiprime has Z_N ≥ |FB| + 2.
  - **P1h** (the normalisation transfers to S_n): the mean of r_N over semiprimes is within 1 ± 0.15.

## E3: pair correlations by k (from the E1 integer samples)
- Definition: C(k) = [Σ_N Σ_t f_t f_{t+k} / #N] / [Σ_t F_t F_{t+k} / #N²], where F_t = Σ_N f_t and
  k ∈ [T/8, T/2).
- ratio_ℓ = mean C(k : ℓ | k) / mean C(k : ℓ ∤ k), for ℓ = 2, 3, 5.
- Model values, `prereg_values.py`, from the unconditional local pair factor at u(T/4):

| n | ℓ = 2 | ℓ = 3 | ℓ = 5 |
|---|---|---|---|
| 40 | 1.5116 | 1.2609 | 1.2932 |
| 48 | 1.3958 | 1.2056 | 1.2287 |
| 56 | 1.2023 | 1.1082 | 1.1177 |
| 64 | 1.1677 | 1.0901 | 1.0973 |
| 72 | 1.1077 | 1.0583 | 1.0621 |

- Predictions:
  - **P3a**: ratio_ℓ > 1 for every ℓ and n.
  - **P3b**: (ratio_ℓ − 1)/(model_ℓ − 1) ∈ [0.4, 2.0] for every ℓ and n.
- Purpose: the theory (notes/detqs_results.md) says that after summing over t the pair sum is a fat-region,
  BDH-type sum. Its only non-product structure is *local*: congruence classes of k. There is no
  thin-region anomaly. Large k (≈ T/2) is included.

## E2: H2 rank event on S_n
- Samples: n ∈ {40, 48, 56, 64} with 1000 semiprimes each; 400 at n = 72; 200 at n = 80.
- d* is the index of the first fruitful kernel vector in t-order. It is ∞ if none occurs within DMAX.
- Predictions:
  - **P2a**: for Blum N, the Jacobi rule ((X_S Y_S/N) = −1 ⇔ fruitful) holds for every tested vector.
  - **P2b** (fair-coin rank heuristic): for j = 1..5, #{d* > j} is consistent with Binomial(M, 2^{−j}),
    with exact two-sided p ≥ 0.001.
    - Per n: at least 28 of the 30 (n, j) tests pass.
    - Pooled over n: all 5 pass.
  - **P2c**: #{d* > 15} ≤ 2 over all n. The expectation is ≈ 0.14.
  - **P2d**: H1 holds for every run (Z ≥ |FB| + 2), with 0 ACCIDENTAL hits.
  - **P2f** (Lemma P): every tested kernel vector of size ≤ 2 is fruitful.

## Decision rules
- **H1 structure supported:** P1a–P1h and P3a–P3b hold at every n, with at most 2 individual
  (prediction, n) failures, each reported.
- **H2 heuristic supported:** P2a, P2b, P2c and P2f hold.
- Anything else is reported as "mixed", with the failing items.

## Commands
- `PYTHONPATH=. .venv/bin/python -m src.detqs.harness e1`
- `PYTHONPATH=. .venv/bin/python -m src.detqs.harness e2`
- `PYTHONPATH=. .venv/bin/python -m src.detqs.harness summary`
