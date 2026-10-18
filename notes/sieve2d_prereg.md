# sieve2d pre-registration (agent `sieve2d`, 2026-10-18)

Committed before any pre-registered run. Code: `src/sieve2d/{coverage,sampling,experiments}.py`
(harness only: it uses p and q). Tests: `tests/test_sieve2d.py` (13 pass). Config:
`experiments/configs/sieve2d.json` (= `CONFIG` in `experiments.py`).

**Disclosure.** Before writing this file I ran ONE smoke/timing job (not part of the
pre-registered runs): n = 96, block 0, Y = 32768, 10^5 sampled prime pairs, H ∈ {8,…,24}.
Observed uncovered counts were within ≈2% of the independent-residue model for
λ ≤ 8 (e.g. H = 16: 578 obs vs 600 model; H = 20: 27 vs 26.6; H = 24: 1 vs 0.9).
The predictions below were set with that knowledge. The pre-registered runs use 10^6
pairs, so the sample streams differ from the smoke run.

## Question

Theorem D′'s exceptional set E_b consists of pairs (p,q) ∈ S_n that are *uncovered*: no
(v, ℓ) ∈ V × 𝒫_k with ℓ | Bp − Aq. The proved bound is n²/H² (large sieve). The heuristic
is e^{−λ} (Poisson), λ = |V| Σ_{ℓ∈𝒫_k} 1/ℓ. The theory side (notes/sieve2d_results.md)
proves λ^{−2}-type bounds. Here we measure the actual uncovered fraction.

## Designs

- Theorem-shaped band ("thm"): Y = ⌈⌊2^{n/6}⌋/⌊H^{1/3}⌋⌉, blocks k: primes ℓ ≡ 3 (mod 4)
  in (2^kY, 2^{k+1}Y]; V = V(H) of D′.
- Separated band ("sep", exhaustive only): Y = H² + 1.
- Model: p/q mod ℓ uniform on (Z/ℓ)^* independently (prime pairs), or (p,q) uniform on
  (Z/ℓ)² (integer pairs). f_model = ∏(1 − cover prob). λ_eff = −log f_model. Model
  moments of M = #{(ℓ,v): ℓ | Bp − Aq} are exact (tested against residue enumeration).

## Runs

- **E1 (exhaustive).** Every prime pair of S_n, n ∈ {20,22,24,26,28,30}; every
  nondegenerate integer pair (x<y≤2x, 2^{n−1} ≤ xy < 2^n, no Bx = Ay) for
  n ∈ {20,22,24,26}. thm: H = 2..12, blocks 0–3. sep: H = 2..10, blocks 0, 1.
- **E2/E3 (sampled).** 10^6 uniform prime pairs of S_n and 10^6 uniform integer pairs
  of the region, seed 20261018, n ∈ {48, 64, 80, 96}; thm design, H ∈
  {2,3,4,5,6,7,8,10,12,14,16,20,24,28,32} restricted to Y > H² (classes distinct, the
  theorem's regime); blocks 0–3 at n = 48, 64, block 0 at 80, 96.

## Predictions (frozen)

Let E = pairs × f_model (expected uncovered count), obs = observed uncovered count.

- **P1 (model agreement).** E2, prime and integer: every cell with E ≥ 30 has
  obs/E ∈ [0.85, 1.15]; pooled Σobs/ΣE over those cells ∈ [0.95, 1.05]. E1 (small n, where
  P/Y is small): every cell with E ≥ 30 has obs/E ∈ [0.67, 1.5] (descriptive; failures
  reported, not a decision input).
- **P2 (exponential tail, not polynomial).** E2 block 0: pooled over cells with E < 5:
  Σobs ≤ ΣE + 3√ΣE + 3. And every E2 cell with λ ≥ 6 has obs/pairs ≤ 0.1 · 3/λ²
  (at least 10× below the λ^{−2} curve).
- **P3 (slope).** Least squares of log(obs/pairs) on λ_eff over E2 block-0 prime cells
  with obs ≥ 10, weights obs: slope ∈ [−1.15, −0.85]. Same fit with λ (nominal) in place
  of λ_eff: slope ∈ [−1.3, −0.85].
- **P4 (primes vs integers).** Over E2 cells with E ≥ 30 for both kinds (same n, k, H):
  pooled (Σobs_int/ΣE_int)/(Σobs_prime/ΣE_prime) ∈ [0.9, 1.1].
- **P5 (moments; the quantity the theorem controls).** E2 cells with λ ≥ 1:
  Var(M)/model_var ∈ [0.9, 1.1] and m4(M)/model_m4 ∈ [0.8, 1.2] in every cell;
  for integer cells, Qstat/Qmodel ∈ [0.8, 1.2], where Qstat = mean of (1 − M/λ)^4 over
  the nondegenerate pairs (the theorem's weight-free statistic).
- **P6 (uniformity in q).** E1 n = 30, thm block 0, cells with f ∈ [0.01, 0.5]: the
  dispersion index D of per-q uncovered counts (q with ≥ 200 partners) ∈ [0.5, 2].

## Decision rules

- P1 ∧ P2 ∧ P3 ⇒ "uncovered fraction ≈ independent-residue model ≈ e^{−λ}
  (experimentally supported at 48–96 bits)". Any failure is reported as such.
- P4 ⇒ "restricting to primes vs integers changes the rate by < 10%".
- P5 ⇒ "the 4th-moment main term used in the proof is numerically accurate".
- None of these tests the proved bound itself (an upper bound); they test the heuristic
  and the proof's main-term identity.

## Reproduction

```
PYTHONPATH=. .venv/bin/python -m src.sieve2d.experiments config
PYTHONPATH=. .venv/bin/python -m src.sieve2d.experiments exh
PYTHONPATH=. .venv/bin/python -m src.sieve2d.experiments samp
PYTHONPATH=. .venv/bin/python -m src.sieve2d.experiments summarize
```
