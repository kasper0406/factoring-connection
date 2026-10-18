# nfsfruit pre-registration: experiment E1, the fruitfulness functional on toy sieve families

Agent `nfsfruit`, 2026-10-18. This file is committed BEFORE any E1 run. A pilot of at most one
instance per family and size (seed 9001, Blum, 40/44/48 bits) was used only for timing and
box sizes. It is excluded from E1, and E1 uses seeds 0..19.

## Object measured

Let n be composite and fix a relation system: relations i = 1..M, and a GF(2) parity matrix A
(rows are relations; columns are signs, factor-base primes per side, first-degree prime ideals (r,s),
and, for NFS, 40 quadratic characters). Let K0 be the space of true square dependencies. For NFS,
K0 is ker A after every basis vector has been certified a square by an exact square root in Z[α].
For S in K0 the algorithm builds X, Y with X² ≡ Y² (mod n).

**Lemma H (proved in notes/nfsfruit_results.md; the experiment checks it).** S ↦ [X/Y] ∈
μ₂(Z/n)/{±1} ≅ F₂^{ω(n)−1} is a group homomorphism Φ: K0 → F₂^{ω−1}, for every n (Blum or not).
Hence a uniformly random S in K0 is fruitful with probability exactly 1 − 2^{−rank Φ}.

The Lee–Venkatesan conjecture (Conj. 7.1, Blum n) is equivalent to rank Φ ≥ 1 on the full set of
smooth relations. "The factor-detecting row raises the rank" means exactly rank Φ|K0 ≥ 1.

## Families (all exact; src/nfsfruit/families.py)

| family | relation | role |
|---|---|---|
| dixon | x uniform mod n, x² mod n smooth | positive control: hidden root, provably random |
| trivial | x < √n, value x² (no reduction) | negative control: Φ ≡ 0 by construction |
| qs | x = ⌊√n⌋+1+i, x² − n | heuristic, x known |
| lin_base | degree-1 NFS, L1 = x − my, L2 = C1x + C0y, m = ⌊√n⌋ | rational–rational NFS, fixed polynomial |
| lin_rand | degree-1 NFS, m uniform mod n, L1, L2 = Gauss-reduced lattice vectors (+ t·L1, |t| ≤ 2) | randomized polynomial |
| nfs3_base | monic cubic, base-m, m uniform in ((n/2)^{1/3}, n^{1/3}] | standard NFS |
| nfs3_rand | the same plus (x − m)(r1 x + r0), r_i uniform in [−⌊m/4⌋, ⌊m/4⌋] | LV-style randomization, monic-restricted |

Disclosed simplifications relative to LV: degree 3 is fixed (not growing); the polynomial is monic,
so LV's c₀ term is dropped; |r_i| ≤ m/4 (LV need |c_i| ≫ m, which kills the yield at toy sizes);
characters are 40 first-degree primes near 2²⁶ (not exp L(1/3) norms); the relation search stops at
(#columns + 60) instead of exhausting the box.

Parameters (experiments/configs/nfsfruit_e1.json):
- dixon B=300 (40b) / 400 (44b), excess 40;
- trivial B=300, excess 40;
- qs B=400, M=20000;
- lin B=1500, A=600, Bb=200 (40b) / B=2000, A=800, Bb=250 (44b), excess 60;
- nfs3_base B=1500, A=1000, Bb=400 (40b) / B=2500, A=2600, Bb=1300 (44b), nchar 40, excess 60;
- nfs3_rand B=1500, A=1800, Bb=900 (40b) / B=2500, A=3600, Bb=1800 (44b), nchar 40, excess 60.
  (Box enlarged after the seed-9001 pilot showed too few relations; this was fixed before E1.)

Inputs: bits ∈ {40, 44}; classes blum (p ≡ q ≡ 3 mod 4), mixed, both1 (p ≡ q ≡ 1 mod 4), and
tri (3 primes). 20 seeds per (family, class, bits), so 7 × 4 × 2 × 20 = 1120 instances.
Generator: src/nfsfruit/instances.py.

## Per-instance measurements

M, #columns, rank A, dim ker(parity only) (NFS), dim K0, rank Φ, the number of fruitful basis vectors,
4 direct linearity checks Φ(S⊕S′) = Φ(S)+Φ(S′), and 40 random relation subsets (size uniform in
rank A ± 8) recording (dim K0(sub), rank Φ|K0(sub)). Φ on subsets is computed from coordinates in
the full K0 basis; this uses Lemma H, which the linearity checks validate.

## Predictions and decision rules (frozen)

- **P0 (correctness):** 0 failed congruence checks, 0 Lemma-H linearity violations, 0 non-square
  basis vectors after certification. Any violation is a bug or refutes Lemma H. Investigate before
  reporting anything else.
- **P1 (controls):** trivial has rank Φ = 0 in 100% of instances. dixon has rank Φ = ω − 1 in 100%
  of instances with dim K0 ≥ 10.
- **P2 (main, heuristic "random functional" H_rand):** for every non-control family, rank Φ = ω − 1
  in every instance with dim K0 ≥ 10.
  - *Structural degeneracy* is declared if any such instance has rank Φ < ω − 1 (expected count
    under H_rand ≤ 1120·2^{−10} ≈ 1).
  - *Family degeneracy* is declared if ≥ 3 instances in one family do so.
- **P3 (subset statistics):** for semiprimes, pool the subsets with dim K0(sub) = k for k = 1, 2, 3, 4.
  Under H_rand, P(rank Φ|sub = 0) = 2^{−k}. For 3-prime n, P(rank < 2) = 1 − (1 − 2^{−k})(1 − 2^{1−k}),
  the probability that a uniformly random linear map F₂^k → F₂² has rank < 2.
  Pooled within a family over classes and bits, *bias* is declared at |z| > 4 (binomial). Prediction:
  no bias in any family, Blum or not.
- **P4 (Blum vs non-Blum):** no difference in rank-Φ failure rates or subset statistics between blum,
  mixed and both1 (|z| ≤ 4). This tests the claim of Lemma H that the Blum condition is irrelevant to
  the linear-algebraic structure.
- **P5 (characters):** for NFS, record dim ker(parity only) − dim K0, the number of character
  constraints actually binding. No prediction; descriptive only.

Analysis: src/nfsfruit/analyze_e1.py over experiments/raw/nfsfruit_e1.jsonl, with the exact
command recorded in the results note.
