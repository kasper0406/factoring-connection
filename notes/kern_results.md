# kern: Frobenius-free counts on the torsion algebra over Z/N; the OR-structure of zero-detection (2026-11-25)

**Process.** GPT-6.2-astra via codex (session 01a21498-3e5f-7d92-a799-a290d1f75bbb). Brief: `codex/kern.md`. Output:
`codex/kern_codex.out`. Artifacts: `experiments/kern/`. The user's sanity check, that the kernel of a ↦ a^N is
trivial, was confirmed.

## Proved

**K1: power maps on A^×.**
- **Orbit data.** Full nonzero vectors of E[ℓ] decide the OR, i.e. whether either side has ℓ-torsion. Roots of ψ_ℓ
  (vectors up to ±1) and projective lines can lose even the zero test: I and −I act identically.
- **Kernel and order.** |A^×[m]| = |A^×/(A^×)^m| = ∏_r ∏_i gcd(m, r^{k_{r,i}} − 1).
  - For m = N this is q^{c_p} p^{c_q}, which is almost always 1: exceptional partners number O(d²).
  - **The exact |A^×| factors N:** φ(N) divides it, and Miller-style square roots of 1 follow.
- **Iterates.** The fixed fraction of a ↦ a^{N^j} is gcd((pq)^j − 1, p^k − 1)/(p^k − 1). For fixed d and j, a return in
  even one component has probability X^{−1+o(1)}.

**K2: the OR.**
- **Square-class form of the target.** For y² = x³ − x + 1 and A = (Z/N)[x]/(x³ − x + 1):
  R(N) = 4 ⟺ M odd ⟺ t_p = t_q = 1 ⟺ |A^×[2]| = |Idem(A)| = 4. Any non-scalar involution or idempotent certifies
  the OR and exposes a factor.
- **Sampling reduction.** Let T = ker(Norm : A^× → (Z/N)^×). Draw a uniform unit a and set b = a³/Norm(a). Then
  Pr[b is a square in A^×] = 2^{2−s}, with s = t_p + t_q. So every draw is a square when R(N) = 4, and at least half
  are non-squares otherwise.
  - That is a one-sided constant-gap test, missing exactly one step: **recognizing non-squares among norm-one units**,
    whose norm-Jacobi symbol is identically +1.
  - It is a quadratic-residuosity problem (cf. Goldwasser–Micali QRP) in the norm-one torus of a cubic algebra.
- **Trace forms lose the distinction.**
  - The Gram matrix of the trace pairing has det −23. Over odd finite fields, dimension and determinant classify
    nondegenerate forms, and the Brauer group vanishes, so Hasse–Witt carries no extra bit. This is consistent with
    Serre's formula.
  - Exhaustive check over 4,186 prime pairs p < q ≤ 499: both local Hasse signs were +1 everywhere. In the stratum
    where both discriminant signs are +1, R = 4 held on 496 pairs and R ≠ 4 on 450. The discriminant is a genuine
    partial filter, but not enough.

## Measured (pre-registered)

1,200 semiprimes each at 40, 60 and 80 bits.
- Whole-element returns: 0.
- Exact |ker μ_N| > 1: 0.
- Proper gcd splits: 6 at 40 bits, 0 at 60 and 80.
- On 1,769 non-split inputs with positive discriminant symbol, the max |correlation| over eight signed readouts was
  0.027, family-wise p = 0.91. Null.

## Sharpest open statement

On inputs with (−23/N) = +1, decide whether every norm-one unit of the cubic algebra (Z/N)[x]/(x³ − x + 1) is a
square: a generalized quadratic-residuosity problem in the torus T. A factor-free certificate of non-squareness on
even an inverse-polynomial fraction of the draws would decide u₂(N). Astra puts this route at about 1%, conditional on
factoring ∈ P.
