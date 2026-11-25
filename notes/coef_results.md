# coef: u_ℓ(N) = N-th coefficient of ζ(s)ζ(s−1)/L(E,s) mod ℓ, without factoring (2026-11-25)

**Process.** GPT-6.2-astra via codex (session 01a21498-3e5f-7d92-a799-a290d1f75bbb). Brief: `codex/coef.md`. Output:
`codex/coef_codex.out`. Artifacts: `experiments/coef/`.

**Motivation.** The user believes rumours that classical poly-time factoring exists. Claude's framing: Jacobi symbols
are computable from N because characters are multiplicative *and* given by congruence conditions. u_ℓ(N) needs
non-abelian Frobenius data at the composite N.

## Proved

**Structure.** u_ℓ(pq) = h(g_p)h(g_q), with h(g) = det(1−g) = 1 + det g − tr g. Expanded:
u_ℓ(pq) = (1+p)(1+q) − (1+q)a_p − (1+p)a_q + a_p a_q mod ℓ.

**Two distinct information losses.**
1. **Combining the primes already loses information, even in an abelian group.** The pairs (1, g²) and (g, g) have
   the same product but give answers 0 and h(g)².
2. **For full image, separate abelian data at each prime is insufficient.** Examples: I vs −I (h = 0 vs 4); in S₃,
   identity vs a 3-cycle (h = 0 vs 1). By Chebotarev this holds even after fixing any finite set of compatible abelian
   data, so no fixed congruence recipe exists.

Arbitrary operations on the integer N remain open.

**Small cases, exhaustively.**
- **ℓ = 2:** u₂ = 1 iff both 2-division cubics are irreducible mod p and mod q (both Frobenius classes are 3-cycles).
  Then u₂ = a_N mod 2. Neither (disc/N) nor the root count of the cubic mod N determines it.
- **ℓ = 3:** enumerating all 48² matrix pairs, all three answers occur in each public determinant class. (I,I) and
  (−I,−I) agree on a_N mod 3 and on all four isogeny-kernel actions, yet give 0 vs 1. A subgroup census shows which
  images need non-abelian local data.
- **Does the decision factor N?** Inside Theorem PT's pure model, an exact evaluator splits a positive fraction of
  inputs (1/18 for u₂, 81/512 for u₃). For unrestricted algorithms no reduction is known: one bit at N alone gives
  nothing.

**A native parity formula (new, verified).** Take E: y² = x³ − x + 1. Its 2-division cubic's splitting field is the
Hilbert class field of Q(√−23), with h = 3. Let R(N) = #{(s,t) ∈ Z² : s² + st + 6t² = N}. Then for distinct good
primes p, q ∤ 46:

  **u₂(pq) = 1 ⟺ R(pq) = 4** (as opposed to 0 or 8).

- Equivalently, u₂(N) is the N-th coefficient of η(z)η(23z) mod 2 on these semiprimes.
- Checked on 1,891 semiprimes.
- The direct algorithm, scanning t, costs O(√N).

**A tempting acceleration fails.** No finite automaton outputs u₂(n) for all n. The argument is density zero versus
≫ X/log X primes; via Christol, the generating series is not algebraic over F₂(z). The semiprime promise is left
open.

**A recovery check.** u(r²) = (1+r)u(r), so an all-index evaluator gives p + q ≡ u_ℓ(N²)/u_ℓ(N) − 1 − N (mod ℓ). A
semiprime-only evaluator, uniform for ℓ ≤ C log N, gives M by CRT, and Theorem EXT then factors N.

## Recent results mined, as a lever on non-abelian data at N

- **Cogentic deterministic L[3/4]:** its tools are abelian cocycles and positivity. Its splitting step needs a
  multiple of r−1, and its short-generator machinery would need an ideal of norm N with chosen primes above p and q.
- **Cogentic L[1/3]:** it factors first.
- **alpo.ge CM derivatives:** the closest positive example, a blind Frobenius-type detection mod N, but only at
  derivative order ≥ p, and inertness is quadratic (abelian).
- **OpenAI FFT:** full transforms, with no random access to a single coefficient.
- **Masullo GRH:** chooses auxiliary primes, with no handle on the divisors of N.

None gives a lever.

## Experiment (pre-registered)

**Method.** 800 balanced semiprimes each at 40, 60 and 80 bits. A sampling estimator built only from N: random
multiples of a point, torsion certification, Jacobi readouts.

**Results.**
- u₂ gain +0.0011, 99% CI [−0.003, +0.005]: no improvement over the baseline.
- u₃ gain −0.0135, slightly worse than the baseline.
- 55/800 accidental splits at 40 bits, 0 at 60 and 80.

Null.

## Sharpest open statement

On good semiprimes N with (−23/N) = +1, decide whether R(N) = #{s² + st + 6t² = N} equals 4, in poly(log N) time.
That is one bit of a fixed weight-1 dihedral form's coefficient, a_f(N) mod 2 for f = η(z)η(23z). Astra puts the
native-coefficient route at about 3%, conditional on factoring ∈ P.
