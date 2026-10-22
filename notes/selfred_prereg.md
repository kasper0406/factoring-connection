# selfred pre-registration (phase 22)

Agent `selfred`, 2026-10-22. Committed before any experiment below is run. Both experiments are exact (integer
or rational arithmetic; floats only for display). Code: `src/selfred/`. Raw output: `experiments/raw/selfred_e{1,2}.jsonl`.

## E1. Exact census: what N mod 2^j says about ε (positional floor beyond Siegel–Walfisz)

**Claim being illustrated (Prop FLOOR-m, to be proved in notes/selfred_results.md).** If X, Y are iid on primes and
TV₁ := TV(law of X mod 2^j, uniform on units mod 2^j), then for every function F of N mod 2^j (j ≥ 2), with
N = XY and E = {N ≡ 1 mod 4},

  A*_j := max_F Pr[F(N mod 2^j) = ε | E] − 1/2 ≤ 3·TV₂/(1 − 2·TV₂) ≤ 6·TV₁/(1 − 4·TV₁),   TV₂ ≤ 2·TV₁,

where ε = (−1/X) (= (−1/Y) on E). This is a proved inequality; E1 checks it exactly and shows where A*_j becomes
non-negligible.

**Population.** S_b = all primes in [3·2^{b−2}, 2^b) (top two bits set, as `balanced_semiprime`), b ∈ {14, 16, 18};
n = 2b. Ordered pairs (X, Y) ∈ S_b × S_b, uniform, diagonal included (this is exactly the iid law of the claim).
For each j = 2, …, n: exact A*_j (Bayes-optimal: Σ_c max(n₊(c), n₋(c)) / #E − 1/2) and exact TV₁, TV₂.

**Null (random labels).** For each class c with k_c = n₊(c) + n₋(c) pairs, the expected Bayes advantage if the ε labels
were iid fair coins: E[max(B, k_c − B)] with B ~ Bin(k_c, 1/2), computed exactly. A_null,j = Σ_c E[max]/#E − 1/2.

**Predictions.**
- P1 (proved; must hold): A*_j ≤ 3·TV₂/(1 − 2TV₂) whenever TV₂ < 1/2, for every (b, j). One violation = bug.
- P2 (heuristic, random-model): A*_j / A_null,j ∈ [0.5, 2] for every j with A_null,j ≥ 10^{-3} (i.e. no structure beyond
  random-label fluctuation), and A*_j ≤ 0.02 for all j ≤ n − 2·log₂ b − 2. The GRH bound only guarantees the floor for
  2^j ≲ 2^{b/2}/b²; P2 predicts it holds empirically to j ≈ n − O(log n) (Montgomery-type equidistribution).
- Descriptive: j_{0.05} = least j with A*_j ≥ 0.05, per b.

## E2. Oblivious polynomial maps never move p (Theorem SR1)

**Claim (Theorem SR1, proved in the notes).** For h ∈ Z[x] of degree D (arbitrary height), content c, prime p and
Q = {primes q ∈ (p, 2p)}: G_p := gcd_{q ∈ Q, h(pq) ≠ 0} h(pq) divides p^a · c · S with S an L-smooth integer,
L = 1 + 4Dp/(|Q| ln p) (when |Q| > D(√p + 1)). Refined (Siegel–Walfisz range, proved for p ≥ p₀(h)): G_p / p^{v_p(G_p)}
= U(h) := ∏_{ℓ ≠ p} ℓ^{u_ℓ(h)}, u_ℓ(h) = min_{z: ℓ ∤ z} v_ℓ(h(z)), independent of p.

**Instances.** Seeded (`random.Random("selfred-e2|D|k")`): D ∈ {1, 2, 3, 5, 8}, 8 polynomials per D. Coefficients
uniform in [−2^64, 2^64], multiplied by a random content c = product of 1–3 random primes < 2^20 and, for half of the
polynomials, by a planted fixed-divisor factor (∏_{i=1}^{D}(x − i) · x for D ≥ 2, i.e. divisible by D! on units-ish),
so U(h) is nontrivial. Primes p: the 40 primes following 10^4 and the 40 following 10^5. Q = all primes in (p, 2p).

**Predictions.**
- P1 (proved; must hold): every prime ℓ | G_p with ℓ ≠ p and ℓ ∤ c satisfies ℓ < L. One violation = bug.
- P2 (proved for large p; checked): G_p / p^{v_p(G_p)} = U(h) for every (h, p). Any mismatch is reported with the
  offending prime and p (expected only if some ℓ^{u+1} is too large for (p, 2p) to cover its unit classes).
- Descriptive: number of distinct values of G_p / p^{v_p} per h (prediction: 1).

## Not pre-registered
Unit tests of the reductions (rigidity solver, rad ⇒ squarefree decomposition, valuation split, φ(N) mod ℓ ⇒ s) are
correctness checks, not experiments.
