# unbalanced: pre-registration (extend Theorem D′ to unbalanced ratios q/p)

Agent `unbalanced`, 2026-10-18. Written and committed **before** any pre-registered run. Code:
`src/unbalanced/` (operational: `params.py`, `algorithm.py`; harness: `harness.py`,
`experiments.py`, `analyze.py`). Tests: `tests/test_unbalanced.py` (9 pass). Config with frozen
n-only predictions: `experiments/configs/unbalanced_experiments.json`.

## 1. Theoretical claims to be proved in notes/unbalanced_results.md (status now: drafted)

Notation: ratio R = q/p, band j means R ∈ (2^j, 2^{j+1}], S_{n,j} the prime pairs of band j with
2^{n−1} ≤ pq < 2^n, H the weight height of D′.

- **Theorem U (regime I).** D′ with the weight set replaced by ratio-band instances
  (A, B, band j), A ≤ ⌊H2^{−j/2}⌋, B ≤ ⌊H2^{j/2}⌋, gcd = 1, A ≠ B, window computed over the
  band-j p-range, j = 0..J with 2^J ≤ H²/C₀. Every window span stays ≤ (1+√2)H√N + O(1), so the
  primes, γ, the four blocks and the shared baby list are those of D′. Cost ×O(J) = ×O(n);
  exceptional fraction O(n²/H²) **uniformly in every band**. With H = n2^{εn/2}: ratios up to
  N^{ε} at cost N^{1/6+5ε/6}, exceptions O(N^{−ε}); with H = N^{o(1)}: ratios N^{o(1)} at
  N^{1/6+o(1)}.
- **Remark (D′ unchanged).** D′'s own window for (A,B) contains F(p) whenever 1 ≤ N/p² ≤ (B/A)²
  (B ≥ √2A), so D′ itself covers ratio R with multiplicity ≍ H²/√R (exceptional ≍ n²√R/H²).
- **Theorem U_θ (regime II).** For 2^j > H²/C₀: weights (1,B), B ∈ [2^j, 2^{j+1}), ratio window
  (B−Δ, B+Δ], Δ ≍ H², per-band Y_j ≍ Δ^{2/3}N^{1/6}R_j^{−5/6}, K_θ blocks. Cost per band
  ≍ Δ^{2/3}(NR)^{1/6} = Δ^{2/3}q^{1/3}; valid for R < N^{1/5−O(ε)}. Hence ratios up to N^θ
  (θ < 1/5) cost N^{(1+θ)/6+O(ε)} with O(N^{−ε}) exceptions.
- **Proposition O (curvature obstruction).** For any Lehman-linear target (exponent B + AN,
  lock modulus M), one giant of length m structurally covers ≤ 2√(m·x₂³/(|A|N)) + 2 integer
  candidates in [x₁, x₂]. Full coverage of the candidates of ratio R costs ≥ c(NR)^{1/6} = c q^{1/3}.
  So within this family no design reaches N^{1/6} at ratio N^{θ}, θ > 0. Combined with
  Pollard–Strassen (cost p^{1/2}), the worst ratio is R = N^{1/5}, at cost N^{1/5}.

## 2. Designs (all n-only)

- **C1:** D′'s pre-registered finite schedule, `dcltheorem.params.schedule_design(n, 3.0, K=4)`.
- **W1 (regime I):** H, Y and K = 4 from C1; ratio bands j = 0..J1, where J1 is the largest j with
  2^j ≤ H². At n = 48/64/80: H = 9/11/11, Y = 206/1461/7804, J1 = 6.
- **W2 (regime II, per band j):** `schedule_regime2(n, j, 3.0)`. Weights (1,B), B ∈ [2^j, 2^{j+1}),
  ratio window (B−Δ, B+Δ]. Y is chosen on a 2^{1/4} grid by the n-only expected block-0 cost. Δ is
  the least integer with (2Δ+1)Σ_{ℓ∈𝒫₀}1/ℓ ≥ 3. Feasibility requires Y > 2Δ + 1.
- The algorithm is `unbalanced.algorithm.dclw`. With the single band j = 0 it reproduces D′'s target
  list exactly (test `test_band0_reproduces_dprime_targets`).
- Instances are uniform prime pairs from S_{n,j} (`harness.sample_pair_ratio`, rejection from an
  integer box).

## 3. Experiments

| id | what | sizes | seeds |
|---|---|---|---|
| lemma | W1, all 4 blocks recorded, Lemma 1–7 checks on every covering target, certificate | n ∈ {32, 40}, j = 0..J1 | 0–9 |
| oper1 | W1 and C1 on the **same** N, all 4 blocks | n ∈ {48, 64, 80}, j ∈ {0, 2, 4, J1} | 100–123 |
| oper2 | W2, all 4 blocks, lemma checks | (n,j) ∈ {(64,7),(64,9),(80,7),(80,9),(80,11),(96,7),(96,9),(96,11),(96,13)} | 200–215 |
| exc | W1 certificate components only (no stage 3) | n ∈ {64, 80}, j ∈ {0, 3, J1} | 5000–5299 |

Metric E is `dcltheorem.algorithm.metric_E`, unchanged.

**Frozen n-only predictions** (config file):
- W1 block-0 multiplicity λ_j at n = 64: 3.75, 3.20, 3.48, 3.07, 2.97, 2.79, 3.98 for j = 0..6.
- C1 block-0 multiplicity at the band-centre ratio, n = 64: 3.75, 1.10, 0.78, 0.55, 0.37, 0.23,
  0.09 for j = 0..6.
- W2 expected block-0 list:
  - n = 80: 39 722 / 47 318 / 58 370 at j = 7/9/11 (log2 slope ≈ 0.14 per j);
  - n = 96: 259 855 / 300 052 / 373 730 / 466 253 at j = 7/9/11/13 (slope ≈ 0.14);
  - asymptotic theory: 1/6.

## 4. Decision rules (fixed now)

- **R-a (correctness).**
  - 0 incorrect outputs in all runs.
  - 0 violations of L1 (integrality), L2 (window), L3 (mod p), L4 (mod q closed form), L57 (clean
    block: gcd = p ⇔ not FULL) on every covering target in lemma, oper1-W1 and oper2.
  - certified ⇒ success, everywhere.
- **R-b (regime I is uniform in the ratio).**
  - W1 success ≥ 23/24 at every (n, j) of oper1.
  - exc: W1 certified fraction ≥ 0.99 at every (n, j).
  - exc: block-0 uncovered fraction ≤ e^{−λ̄₀} + 0.04, with λ̄₀ the mean per-pair n-only
    multiplicity.
- **R-c (regime I overhead is a constant ~J, not a power).**
  - median E_full(W1) / median E_full(C1) at j = 0 lies in [2, J1 + 2] at each n.
  - That ratio at n = 80, divided by the ratio at n = 48, lies in [0.5, 2].
- **R-d (D′ unchanged degrades with the ratio).**
  - C1's mean block-0 structural coverage (p-side + q-side targets whose true g lies in the
    window) at j = J1 is ≤ 0.3 × its value at j = 0, at each n.
  - C1 success at j = J1 is below W1 success at j = J1 at each n (secondary).
- **R-e (regime II cost law).**
  - The OLS slope of log2(median E_block0) against j lies in [0.08, 0.24] at n = 80 and at
    n = 96 (n-only expected ≈ 0.14; asymptotic 1/6).
  - W2 success ≥ 14/16 at every (n, j).
  - A slope ≤ 0.04 at both n would **refute** the R^{1/6} law at these sizes.
- **R-f (calibration).** The measured median block-0 list (m + giants) is within a factor 1.5 of
  the n-only expected block-0 list at every oper2 (n, j).

**Pollard–Strassen comparison: analytic only (no run).** PS removes every factor ≤ X in
X^{1/2+o(1)}. The note compares the exponent profiles (1+θ)/6 against (1−θ)/4.

## 5. Commands

```
.venv/bin/python -m pytest -q tests/test_unbalanced.py
PYTHONPATH=. .venv/bin/python -m src.unbalanced.experiments config
PYTHONPATH=. .venv/bin/python -m src.unbalanced.experiments lemma
PYTHONPATH=. .venv/bin/python -m src.unbalanced.experiments oper1
PYTHONPATH=. .venv/bin/python -m src.unbalanced.experiments oper2
PYTHONPATH=. .venv/bin/python -m src.unbalanced.experiments exc
PYTHONPATH=. .venv/bin/python -m src.unbalanced.experiments summarize
```
