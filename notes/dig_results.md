# dig: native #E(Z/N) and adaptive interval refinement (2026-11-25)

**Context.** Follow-ups to notes/inp_results.md, run by GPT-6.2-astra via codex:
- Dig A, native point count: session 01a21498-3e5f-7d92-a799-a290d1f75bbb; output `codex/dig_ecount_codex.out`.
- Dig B, interval refinement: session 01a214a3-b118-7541-986c-3c41e1e709a0; output `codex/dig_interval_codex.out`.
- Briefs: `codex/dig_{ecount,interval}.md`. Artifacts: `experiments/dig/`.
- Also an explainer of Theorem EXT and the Schoof gap: `codex/explain_ext_codex.out`.

## A. Native #E(Z/N) of one fixed curve

**Verdict.** No native algorithm, and every known method hits a precise obstruction. Weight conditional on
factoring ∈ P revised from 24% to **about 5%**.

**The minimal target is a scalar, not an operator.** It suffices to compute
u_ℓ(N) = det(1−g_p)·det(1−g_q) mod ℓ for primes ℓ = O(log N), where g_r is Frobenius on E[ℓ]. Equivalently: compute
the N-th Dirichlet coefficient of ζ(s)ζ(s−1)/L(E,s) mod a supplied ℓ ≤ 2 ln N, uniformly, without factoring. CRT over
such ℓ then gives M = #E(Z/N). Constructing the Frobenius *operators* is harder and itself factors N (Prop FROB).

**Method by method** (each uses the local characteristic at an identifiable step):

| Method | Where it needs p |
|---|---|
| Schoof / SEA | Frobenius as x^r mod (ψ_ℓ, r); Elkies eigenvalues; Atkin orbit lengths |
| Pila / Edixhoven–Couveignes | Frobenius at a specified prime |
| Kedlaya / Lauder | an r-adic Frobenius lift; a known base-fibre Frobenius matrix |
| Satoh / AGM | the canonical degree-r Frobenius isogeny |
| Witt vectors | the index-N operation is not componentwise p, q |
| Harvey's average polynomial time | a batch over all primes ≤ B; reaching the hidden prime needs B ~ √N |

**Reformulations that avoid constructing Frobenius.**
- The Weil pairing gives only the public product pq mod ℓ. Projective action can lose the answer (I vs −I on E[3]).
- Lifting to Z/N^k: #E(Z/N^k) = N^{k−1}M, which just repeats the unknown.
- Jacobians and Kummer surfaces raise the degree, and Kummer coordinates can discard the needed sign.
- X₀(N) counts over small fields are factoring-equivalent, but the explicit models have size ~N.
- **Exact identity:** M = N + 1 + (p+q) + a_N − (q+1)a_p − (p+1)a_q, so even a supplied a_N leaves mixed terms.
- CM, supersingular and anomalous identities need promises on the hidden primes.

**Partial information.**
- M is odd iff the cubic has no roots mod p or mod q.
- For an S₃ cubic, (disc/N) = −1 forces M even.
- There is no general evaluator even mod 2 or 3.

**How much of M recovery needs.** For q/p < 4/3, |M − N| ≤ 4.1N^{3/4}, so residues with combined modulus > 2W
(about ¾ of the input bits) recover M. A poly-size list of lifts also suffices.

**Experiment.** Curve y² = x³ − x + 1; 2,400 semiprimes each at 40, 60 and 80 bits. Four N-only readouts: Jacobi
symbols of norms of x^N − x + c in the cubic and 3-division algebras.
- Max within-stratum correlation 0.0155, family-wise p = 0.886.
- Planted control 0.266, against a null cutoff of 0.033.
- No splits.

## B. Adaptive interval refinement

**Verdict.** Coppersmith subdivision is exact out to N^{1/4}·poly(log N). The gap starts at N^{1/4+ε}, and no
reproducible signal was found there. Weight revised from 16% to **about 8%**.

**Beam theory.** With independent per-node errors at rate e, even the optimal width-h^k beam needs
δ = 1/2 − O(1/n), i.e. near-perfect per-step accuracy. Correlated errors can change this. A resamplable
inverse-polynomial edge on each true parent suffices by majority vote. An average edge on fresh moduli does not.

**Scores examined.**
- Exact and weighted divisor counts.
- Truncated and smoothed Voronoi sums: the error is ≫ ½ for polynomial truncation, and smoothing admits
  near-divisors of N + k.
- Primal, dual and shifted lattices: coverage at most KW/X.
- Continued fractions, Lehman/Hart and Fermat tangent tests.
- Reciprocal exponential sums: a uniform Fourier detector needs Ω(√N) frequencies (Parseval).

**Proved: a fixed-exponent coverage barrier for tangent square tests.** A menu of M weight pairs with offsets ≤ J has
verified-hit probability ≪ M(1 + N^{1/4}√(J+1))/X, i.e. poly(n)·N^{−ε} at X = N^{1/4+ε}.

**Experiment.** 384 semiprimes per cell at 48, 72 and 96 bits.

| Width | Correct / 1152 (chance 576) |
|---|---|
| N^{0.26} | 665 (exact hits only, fading with size) |
| N^{0.30} | 616 |
| N^{0.40} | 585 |

A 96-bit N^{0.40} cell scored 58%. Its frozen 3,072-trial replication gave 49.2%, CI [47.5, 51.0], so it did not
replicate.
