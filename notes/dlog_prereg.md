# dlog pre-registration: E1/E2 toy runs of the Theorem-DL pipeline

Agent `dlog`, phase 11, 2026-10-20. This file is frozen at the commit that adds it, together with
`src/dlog/{pipeline,harness_e1,params}.py`, `experiments/configs/dlog_e1.json` and `tests/test_dlog.py`.
Calibration used seeds 900–906 (two to three primes each run, 30–44 bits) before freezing. They are not part of E1/E2.

## What is tested

The toy pipeline is a correctness smoke test of the DL post-processing claimed in `notes/dlog_results.md`
(Proposition XK: exact integer dependencies plus archimedean unit elimination, and the hidden-exponent extraction
Lemma HX). It does not test any asymptotic claim. Field selection and special-q descent are imported unchanged from
the Theorem-A toy (`src/thmA/pipeline.py`, design and deviations in `notes/thmA_prereg.md`).

Operational input: (p, g, h, config, seed), where g is a primitive root. The harness (sympy, factor-aware) makes the
instance and checks the answer. The operational code never sees X or the factorisation of p − 1.

## Disclosed toy deviations

1. Relaxed split shape (B-smooth part × ≤ 2 primes in (B, 4B]), Pollard–Brent as the shape test, toy w, H, S. These
   are the thmA toy settings, d = 3, μ = 1/3.
2. Atoms whose norm F_A(P) has a prime factor dividing D = f_d·disc(f) are discarded as descent failures. The
   theorem instead computes local maximal orders at those primes.
3. Miller–Rabin instead of AKS. The flint dense integer nullspace replaces Dixon p-adic lifting.
4. Unit elimination uses arb balls at precision 2·(bits of y) + O(dim) and LLL. A combination counts as torsion
   only if every certified |coordinate| is below half the Blanksby–Montgomery bound.
5. E2 adds NFS pair relations g₁ψ(a − bα) = G(a, b) (rows with T = 1 and no hidden exponents). They exist only to
   make the kernel contain nontrivial algebraic parts γ ≠ 1, which stress-tests unit elimination. They are a
   diagnostic. The theorem's fruitfulness argument does not cover rows without hidden exponents, so P5 uses E1 only.

## Blocks

- E1: 24 primes, bits ∈ {30, 34, 38, 42, 46, 50} × seeds 0–3, pipeline as specified (`extra_rows` = 6).
- E2: 12 primes, bits ∈ {30, 36, 42} × seeds 100–103, plus `pair_rows` = 400 from |a| ≤ 300, 1 ≤ b ≤ 300.

## Predictions (frozen)

- **P1 correctness.** E1: 24/24 and E2: 12/12 return X with g^X = h. Any run that hits a cap counts as a failure.
- **P2 row identity.** 0 violations of T ≡ ±∏p^{e}·g₁^{E}·∏ψ(a − bα)^{e} (mod p), over all target and pair rows.
- **P3 torsion.** For every certified torsion combination x*:
  - W = ∏T_i^{x*_i} ∈ {±1} (μ_K = {±1}, d odd);
  - ψ_q(γ_{x*}) ∈ {±1} (mod q) at 3 auxiliary 40-bit primes q.
  0 violations in E1 ∪ E2.
- **P4 exactness/units.** 0 kernel violations (x^T E ≠ 0) and 0 non-unit kernel vectors (certified Σ Log ∌ 0).
  unit_rank_implied ≤ r₁ + r₂ − 1 in every run.
- **P5 hidden-exponent law (E1 only).** For ℓ ∈ {2, 3, 5, 7} with ℓ | p − 1, take the values B_k mod ℓ over a
  subset of relations that is linearly independent mod ℓ (≤ 20 per run). Pool over E1. The rate of B_k ≡ 0 (mod ℓ) is
  consistent with 1/ℓ: two-sided exact binomial p ≥ 0.001 for each ℓ with ≥ 30 values.
- **P6 unit elimination is needed and works (E2).**
  - Among kernel basis vectors with γ ≠ 1, the fraction with ψ_q(γ) ∈ {±1} at the first auxiliary prime is ≤ 5%,
    pooled.
  - unit_rank_implied = r₁ + r₂ − 1 in ≥ 9/12 E2 runs.
- **P7 (descriptive, no prediction).** dim ker(E mod 2) − dim ker_Z(E), i.e. how often a GF(2) kernel would contain
  vectors with no integer lift (the class-group/unit-type obstruction that Prop XK avoids).

## Commands

```
.venv/bin/python -m src.dlog.harness_e1 experiments/configs/dlog_e1.json experiments/raw/dlog_e1.jsonl
.venv/bin/python -m pytest -q tests/test_dlog.py
```
