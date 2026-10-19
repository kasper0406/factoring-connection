# polyimpl pre-registration (frozen before any run)

Agent `polyimpl`, 2026-10-19. This file and the code it names are committed **before** experiment
E1 is run. Results go to `notes/polyimpl_results.md`. Raw data go to
`experiments/raw/polyimpl_e1.jsonl`, and the summary goes to
`experiments/analysis/polyimpl_e1_summary.json`.

## E1. Residual mass of cyclotomic group-order exposure at the N^{1/6} budget

**Question.** Theorem CB (notes/polyimpl_results.md) shows the following for (Z/N)^* alone.
Annihilation by any B₂-smooth exponents, B₂ ≤ p^{2/3}, never exposes a prime p with
P⁺(p−1) > p^{0.6687}, and Fouvry gives a positive proportion of such p. The natural extension,
Conjecture CB⁺, adds the other cyclotomic tori T_d, d ∈ {2, 3, 4, 6}, whose orders are Φ_d(p).
Proving it needs a joint large-prime-factor statement that is beyond current sieve methods.
E1 measures, at accessible sizes, the fraction of primes, and of balanced pairs, that no
cyclotomic test d ∈ D = {1, 2, 3, 4, 6} exposes under the most generous exposure criterion. It
also tells us whether codex polyimpl_1's Rank-4 route, a deterministic order-exposure portfolio
with a residual-tail bound η ≤ N^{−1/30−δ}, is even heuristically available for cyclotomic
portfolios.

**Definitions (diagnostic code; it uses factorisations of Φ_d(p) and never of N).**
- Bit sizes k ∈ {20, 26, 32}. The primes are p ∈ [2^{k−1}, 2^k). For each k, draw 20 000 distinct
  primes uniformly from the block: random integers from `random.Random(1019·k)`, rejecting
  non-primes (BPSW via gmpy2.is_prime) and repeats.
- N ≈ p², so the budget is X = N^{1/6} ≈ 2^{k/3} (taken as 2^{k/3}, real) and B₂ = X² = 2^{2k/3}.
- Generous exposure: p is **exposed by d** iff P⁺(Φ_d(p)) ≤ B₂. This ignores the costs of
  stage 1, so it over-counts exposure.
- **Residual** r₁(k) = fraction of p exposed by no d ∈ D.
- r₁^{(1)}(k) is the same fraction for D = {1} alone, i.e. P⁺(p−1) > B₂.
- r₁^{(1,2)}(k) is the same fraction for D = {1, 2}.
- **Pair residual** r₂(k) = fraction of pairs (p, q) with both p and q in the residual. The pairs
  are p < q ≤ 2p, formed as consecutive disjoint pairs of a seeded shuffle (seed 7·k) of the
  sample, keeping a pair only if q ≤ 2p. All sampled primes lie in one dyadic block, so every
  pair is balanced.

**Predictions (heuristic: independent Dickman model; frozen).**
- P1. r₁^{(1)} ∈ [0.33, 0.48] at every k. The model value is ln(3/2) = 0.405, up to local and
  finite-size corrections.
- P2. r₁^{(1,2)} ∈ [0.12, 0.22] at every k. The model value is 0.405² = 0.164.
- P3. r₁ ∈ [0.09, 0.20] at every k. The model value is 0.164 · (1 − ρ(3))³ ≈ 0.141.
- P4. r₂ ∈ [0.008, 0.045] at every k. The model value is r₁² ≈ 0.020.
- P5 (stability). r₁(32)/r₁(20) ∈ [0.7, 1.3]: no decay with size.

**Decision rule (frozen).**
- If P3, P4 and P5 hold: record "cyclotomic exposure portfolios leave a size-stable positive
  residual (experimentally supported at 20–32-bit primes)". Codex Rank 4 is then **refuted for
  cyclotomic portfolios (heuristic + data)**, and Conjecture CB⁺ is recorded as experimentally
  supported.
- If r₁(32) < 0.05, or P5 fails by decay: record that the residual shrinks with size, and do
  **not** claim CB⁺ support.
- Either way, report all numbers. Report the binomial standard errors (≈ 0.0025 for r₁ at
  20 000 samples).

**Not tested here.** Fixed-curve ECM lists. Their residual decays with list length under
independence, which is Ward–Maestro Part 5 territory and ecmavg's HIT.

**Command.** `PYTHONPATH=. .venv/bin/python -m src.polyimpl.residual run`
