# rkrank pre-registration — E1: dual obstructions and rank growth of toy NFS relations (Hypothesis RK)

Agent `rkrank`, phase 8, 2026-10-19. Code: `src/rkrank/{gf2,toynfs,sampler,harness,analyze}.py`, tests
`tests/test_rkrank.py` (frozen with this file). Config `experiments/configs/rkrank_e1.json`.

**Calibration disclosure.** Before writing this file I ran two instances *outside* the registered grid:
(blum, 60 bits, seed 99), using `experiments/configs/rkrank_calib.json`, and (tri, 80 bits, seed 99), using
`rkrank_calib80.json`. I looked at all of their outputs. The thresholds below were chosen after seeing the
60-bit instance. That run also refuted my prior expectation of an archimedean sign obstruction in the monic
toy, and P4 is formulated accordingly. The 80-bit run showed that the box A = 8192, b ≤ 8000 is too small.
The registered box is larger. Neither calibration instance belongs to the registered grid.

## Setting
- Toy field: a monic cubic f = f̂_{N,m} + (x − m)(r₁x + r₀), with |r_i| ≤ 4 and f(m) = N.
- Relations: coprime (a, b) with |a| ≤ A and b ≥ 1, swept until 3.0·cols relations or b = bmax. The
  condition is that a − bm and F(a, b) are both B-smooth, tested exactly. Pairs hitting an index-divisor
  prime are skipped (toy).
- Columns:
  - the rational sign and primes ≤ B;
  - degree-one ideals;
  - exact real-embedding signs [a − bρ_j < 0];
  - 32 quadratic characters.
- The relation list is shuffled (seeded). Prefixes of M = L·cols rows are analysed for L ∈ {1.0, 1.25,
  1.5, 2.0, 3.0}; a level is skipped if not reached.
- Grid: classes {blum, tri} × bits {64 (B = 8192, A = 4096, bmax = 6000), 80 (B = 16384, A = 16384,
  bmax = 12000)} × seeds 0–4, so 20 instances.
- Quantities per level:
  - W = R ∩ (rational block), computed exactly by elimination with the high block first.
  - dim W^⊥ (the rational dual).
  - Gen = the ω functionals w_i(p) = [(p/p_i) = −1]. This is harness-only, since it uses the factors.
  - occ(p) = the number of rows with v_p(a − bm) odd.
  - n_peel = the number of rational primes occurring in no row of the algebraic 2-core. Rows are peeled
    iteratively when they contain an ideal that occurs in exactly one remaining row.
  - c_k = codim of W ∩ F₂^{P_k} in Λ/2 ∩ F₂^{P_k}, where P_k = {p : occ(p) ≥ k} (sign coordinate excluded).
  - Small-square vectors: B-smooth z = x² mod N with |z| ≤ 4N^{2/3}, from the linearised-cell sampler.
    That sampler is Vallée-type, not Vallée's algorithm. Each is tested for z-vector ∈ W.
  - δ_φ = the fraction of rows with φ(row) = 1, for three Liouville-type functionals (all rational prime
    columns, all ideal columns, both) and 10 seeded random dense functionals.
- Definitions:
  - excess_peel = dim W^⊥ − ω − n_peel. It is ≥ 0 by construction when Gen and the peeled e_p are
    independent.
  - sign_in_perp means W ⊆ {u_sign = 0}.

## Predictions
- **P0 (sanity).** W ⊆ Λ/2 (Gen ⊆ W^⊥) in every level-record.
- **P1 (RK modulo the 2-core, i.e. "every obstruction is a peeled prime").**
  - P1a: excess_peel ≤ 4 in every record at levels ≥ 1.5.
  - P1b: excess_peel = 0 in ≥ 80% of records at levels ≥ 2.0.
- **P2 (RK on well-covered support).**
  - P2a: c₃ = 0 in ≥ 90% of records at levels ≥ 1.5, and c₃ ≤ 2 in all of them.
  - P2b: c₁ = 0 in ≥ 80% of records at levels ≥ 2.0.
- **P3 (small squares).**
  - Pooled over records at levels ≥ 1.5, small-square vectors whose support lies in P₃ lie in W at
    rate ≥ 0.95.
  - Negative-z vectors among them, in records without a sign obstruction, lie in W at rate ≥ 0.9.
- **P4 (archimedean obstruction).** In every record, sign_in_perp holds iff some real-sign column equals 1
  on every row. "If" is forced; "only if" is the prediction.
- **P5 (no dense near-obstructions).** At the top level reached, δ ∈ [0.45, 0.55] for the three
  Liouville-type and the 10 random dense functionals, in every instance.

Commands:
```
.venv/bin/python -m src.rkrank.harness experiments/configs/rkrank_e1.json experiments/raw/rkrank_e1.jsonl 10
.venv/bin/python -m src.rkrank.analyze experiments/raw/rkrank_e1.jsonl
```
Any deviation will be recorded as a disclosed deviation in notes/rkrank_results.md.
