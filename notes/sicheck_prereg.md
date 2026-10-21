# sicheck pre-registration (phase 20, 2026-10-21)

Agent `sicheck`. Committed **before** any E1/E3/E4 run. Code at the commit of this file: `src/sicheck/{core,diag,model,
run_experiment}.py`, tests `tests/test_sicheck.py` (7 pass).

## 0. Design smoke runs (disclosed; not data)

Before writing this file I ran three smoke scripts (scratchpad, not committed) that only measured the **density of
y-smooth u near kN** and sieve timings, to choose h. They measured nothing about v. The first two (k = 1..10) were
biased upward by inheritance (j·u₁ ∈ I_j when j·v₁ ≤ h); the third sampled 12 random k ∈ [1, K]. Cell values below
come from the third: h ≈ 3/δ̂ (mean count μ ≈ 3 per interval; B10 has μ ≈ 10). μ̂ will be re-measured in E1.

## 1. Theory items (proved in the notes, checked by code; not predictions)

- T1 (exact): κ*³ = 2/9, c³ = 16/9, ν* = 2/(3κ*), trials decreasing at κ* (κ_min³ = 4/9); SS(b) at s = ν* returns
  base κ* and c³ = 16/9; single-interval exponent max(α, 1 − 2α) minimised at 1/3. `src/sicheck/model.py`, tests.
- T2 (**Prop B**, proved in notes): for Blum N (p ≡ q ≡ 3 mod 4) and SI relations u ≡ v (mod N), a dependency S is
  fruitful iff Σ_{i∈S} β_i ≡ 1 (mod 2), β_i = ((Ω_T(u_i) − Ω_T(v_i))/2) mod 2, T = {ℓ : J_N(ℓ) = −1}. Computable
  from N.
- T3 (**Cor. B**, proved): for Blum N, if the factor base consists only of primes with J_N(ℓ) = +1, or if every
  relation has β_i = 0, every dependency is trivial. A constructor that sees only kN and N can therefore make Prop SI
  fail (constructor `beta0`).

## 2. E1: heuristic H_S for deterministic constructors

**Cells** (balanced semiprimes `balanced_semiprime(n, seed)`, seeds 0 and 1; k = 1..K consecutive; I_k = [kN, kN + h]):

| cell | n bits | y | h | K |
|---|---|---|---|---|
| A | 40 | 100 | 6 300 000 | 6000 |
| B | 48 | 200 | 8 000 000 | 6000 |
| B10 | 48 | 200 | 27 000 000 | 3000 |
| C | 56 | 600 | 3 400 000 | 6000 |
| D | 64 | 1500 | 2 800 000 | 4000 |
| E | 80 | 6000 | 6 600 000 | 2000 |

**Constructors** (functions of kN, h, y only; brute-force: the interval is sieved exactly):
- `next`: smallest y-smooth u ∈ I_k; `last`: largest; `rand`: uniform among the smooth u ∈ I_k (rng seeded by
  (seed, kN));
- `cofactor` (greedy product construction): least smooth u₀ ≤ h such that u₀·⌈kN/u₀⌉ is smooth; output that product,
  so v = (−kN) mod u₀.

**Null (size-matched).** p_null(v) = exact density of y-smooth integers in [v − D, v + D], D = max(200, ⌊v/20⌋).
E = Σ p_null(v_k) over the outputs, O = number of outputs with v_k y-smooth. Parity-split nulls p_e(v), p_o(v): the
same density among even (odd) integers of the window.

**Predictions.**
- **P1 (H_S at matched size).** For each constructor, O/E pooled over both seeds lies in [0.85, 1.15] in every cell
  with E ≥ 50, and pooled over all cells in [0.92, 1.08]. (Reason, notes §3: for k running over consecutive integers,
  kN mod M is equidistributed and nearly independent of the constructor's choice of u mod M, so v mod M is uniform for
  M ≪ h.)
- **P2 (local correlation with k, toy-scale).** For `next` and `rand`, let π₀ be the fraction of chosen u that are even
  (measured). Parity model: for odd k, P(v even) = 1 − π₀; for even k, P(v even) = π₀; prediction
  E_par = Σ [P(v even)·p_e(v) + P(v odd)·p_o(v)].
  - P2a: pooled over cells and seeds, O/E (naive null) on odd k is < 1 and on even k is > 1, for both constructors.
  - P2b: O/E_par ∈ [0.85, 1.15] in each k-parity stratum, pooled over cells and seeds, for both constructors.
- **P3 (residue distribution).** Per (cell, seed): δ̂ = total smooth count / (K(h + 1)), μ̂ = δ̂(h + 1).
  - P3a `next`: P(empty) ≥ exp(−μ̂) (over-dispersion from inherited multiples) and P(empty)/exp(−μ̂) ≤ 2.5.
  - P3b `next`: KS distance between v (non-empty intervals) and the exponential law with rate δ̂ truncated to [0, h]
    is ≤ 0.10. `last`: the same for h − v.
  - P3c `rand`: KS distance of v/h from U[0, 1] ≤ 0.05. `cofactor`: KS distance of v/u₀ from U[0, 1] ≤ 0.05.
- **P4 (duplicates, primitivity).** Among outputs with v smooth: the fraction whose primitive pair (u, v)/gcd(u, v)
  was already produced at a smaller k ("duplicate") is larger for `next` than for `rand` (pooled), and lies in
  [0.05, 0.5] for `next`. The fraction with gcd(u, v) > 1 is ≥ 0.2 for every constructor (pooled).

## 3. E3: heuristic R (fruitfulness, rank)

**Instances.** n ∈ {40, 48, 56} with (y, h) = (150, 1 070 000), (300, 1 120 000), (600, 3 400 000). For each n, the
first 6 Blum and the first 6 non-Blum balanced semiprimes in seed order (classification uses p, q: harness only).
k = 1, 2, … shared by all constructors; each constructor collects distinct primitive relations until it has
(factor-base size) + 20, or k reaches 15 000.
Constructors: `next`, `last`, `rand`, `cofactor` (honest); `qsfb` = `next` over the J+ base {ℓ ≤ 2.2y : J_N(ℓ) = +1},
run only when J_N(2) = +1 (so 2 is in it); `beta0` (Blum only) = smallest smooth u ∈ I_k that does not give a relation
with β = 1. For each run: parity rank r, basis of the dependency space, Φ on every basis vector (Φ = 1 iff the square
root of unity Z ∉ {±1}), Φ on 64 uniform random nonzero dependencies, and the prefix sequence (Φ of each new basis
dependency in order of appearance).

**Predictions.**
- **R1.** Honest constructors, runs that reach fb + 20 relations: Φ ≢ 0 in ≥ 95% of runs; the pooled fraction of
  fruitful random dependencies ∈ [0.42, 0.58].
- **R2.** `qsfb` on Blum N: Φ ≡ 0 in 100% of runs with ≥ 1 dependency (T3, proved). On non-Blum N: Φ ≢ 0 in ≥ 50% of
  runs with ≥ 5 dependencies.
- **R3.** `beta0` on Blum N: Φ ≡ 0 in 100% of runs (T3).
- **R4.** Prop B: 0 mismatches between Φ and Σβ on every checked dependency of every Blum run. Non-Blum: mismatches
  occur in ≥ 1 run (the formula is Blum-specific).
- **R5.** Prefix statistic, honest constructors pooled: P(Φ ≡ 0 on the first d dependencies) ∈ [0.5, 2]·2^{−d} for
  d = 1, 2, 3, 4 (each d with ≥ 40 runs).
- **R6.** r ≤ fb − 2 in every run (quadratic characters mod p and q); r = fb − 2 in ≥ 90% of honest runs that reach
  fb + 20 relations.

## 4. E4: end-to-end toy (operational `factor_si`, brute-force constructor)

n ∈ {40, 48, 56, 64, 72, 80} with (y, h) = (150, 1.07e6), (300, 1.12e6), (600, 3.4e6), (1500, 2.8e6), (2500, 5e6),
(6000, 6.6e6); seeds 0, 1, 2; constructor `rand` at all sizes, `next` up to 64 bits; stop at fb + 20 distinct primitive
relations; 64 random dependencies. The algorithm sees only N, y, h, seed.
- **F1.** A nontrivial factor is found in ≥ 95% of runs; every factor is STRUCTURAL (from a dependency); no ACCIDENTAL
  factors (y < p).
- **F2.** Pooled over runs, (v-smooth outputs) / Σ_k p_null(v_k) ∈ [0.85, 1.15].
- **F3.** Descriptive: intervals used, the empty fraction, duplicates, and the Dickman-ρ prediction ρ(log h / log y)
  versus the observed rate (no threshold registered).

## 5. Execution

One heavy job at a time, `memory_pressure` checked first, RSS logged per record (expected < 500 MB).
```
.venv/bin/python -m src.sicheck.run_experiment E1
.venv/bin/python -m src.sicheck.run_experiment E3
.venv/bin/python -m src.sicheck.run_experiment E4
.venv/bin/python -m src.sicheck.score
```
Raw: `experiments/raw/sicheck_E{1,3,4}.jsonl`. Any deviation will be disclosed in the results note.

## Addendum A1 (2026-10-21, after E3, before E3b; committed before running E3b)

E3's R2 was **vacuous** for Blum N: the J+ base (primes ≤ 2.2y with J_N = +1) gave ≤ 14 relations and no dependency
in 15 000 intervals at 48–56 bits. A smoke run (disclosed; it measured relation yield only, on one Blum N per size)
chose smaller instances. **E3b** (`run_experiment E3b`):
- n = 32 (J+ primes ≤ 300, h = 300 000) and n = 36 (J+ primes ≤ 600, h = 1 000 000); the first 8 Blum and the first
  8 non-Blum balanced semiprimes (seed order) with J_N(2) = +1, per size; constructor `next` over the J+ base until
  |base| + 10 distinct relations or k = 20 000; 32 random dependencies. Control: the same Blum N with the full base
  (primes ≤ the same bound), `next`, until |base| + 10 relations.
- **R2b.** Blum, J+ base: Φ ≡ 0 in 100% of runs with ≥ 1 dependency (Cor. B(a), proved), and ≥ 10 such runs in all.
- **R2c.** Blum, full base (control): Φ ≢ 0 in ≥ 90% of runs.
- **R2d.** Non-Blum, J+ base: Φ ≢ 0 in ≥ 50% of the runs with ≥ 5 dependencies (no proof; the obstruction of Cor. B
  needs p ≡ q ≡ 3 mod 4).
