# Phase 29 synthesis: factoring-equivalent invariants as new views, and their native cost

Orchestrator summary, 2026-10-27, of three notes:
- [`nativeinv`](nativeinv_results.md): the inventory of all 55 FACT rows, plus combinatorial and algebraic deep dives;
- [`nativemod`](nativemod_results.md): fixed-level modular-form coefficients and Hecke traces at composite index;
- [`nativeana`](nativeana_results.md): the analytic precision wall.

**Motivation (the user's point).** A reduction "computing X factors N" shows X is at least as hard as factoring. It
says nothing about how hard X is natively, so each reduction is a new view of the problem. We prefer views whose
natural generalisation is not NP-hard.

## Verdict

None of the 55 factoring-equivalent invariants in the repo has a known native algorithm (one that never factors N)
faster than factoring. The best native costs are a power of N: Õ(N^{1/2}) to Õ(N^{1/3}) for divisor-type invariants,
and N^{1/5} for analytic class-number routes with a √Δ finish. The one subexponential native method is class-group
relation finding, at L[1/2], which is itself a relation method. **No native-time lower bound is proved in general.**
Proved results cover restricted models only:
- Prop LOC: random-oracle local data;
- Theorem LIN: linear character-sum estimators;
- Prop HW: black-box linear extraction;
- Prop LRS: powering in a fixed algebra.

A Theorem-PT-type general lower bound would imply permanent lower bounds, so it is out of reach.

## Why the views collapse

- **Prop TRI (nativeinv, proved).** If f(pq) = F(p)F(q) with quasi-polynomial local data F, then f(N) factors N unless
  F is a single monomial, and otherwise is a fiber bit or an aggregate. This is why 21 of 55 rows are "p + q in disguise":
  φ, σ, ψ, Pillai, #K_{2j−1}(Z/N), class counts of GL₂/SL₂/Heis(Z/N), the lens-space count, necklaces and h(−7N²).
- **Native obstructions met by every construction tried.** Each construction needs one of:
  - p-th powers (Adams operations, necklaces, cat maps);
  - CRT- or divisor-indexed coefficients (class counts, Aff degrees);
  - local element orders (SL₃(Z/N) relations);
  - separate Legendre symbols (lens spaces are exactly QR).
- **The analytic wall (nativeana).**
  - *Local* invariants (σ, r₄, Σ_t H(4N − t²) = 2(N + q)) hold the factor in the local factor at p ≈ √N, which is
    invisible to expansions; even a Coppersmith finish needs relative error N^{−3/4}.
  - *Global* invariants (h, hR, L(1, χ), r₃) converge only like K^{−1/2}, which balances at N^{1/5} with a √Δ finish.
  - The fast L-value algorithms that beat √(conductor) (Hiary for powerful moduli, Vishe, Platt) use the
    factorisation of the conductor.
- **Fixed level does not help (nativemod).**
  - At composite N, the Edixhoven–Couveignes Frobenius step is a FROB instance, or a PRIM program for the literal
    x ↦ x^N.
  - a(N) is a convolution of the two Frobenius classes; there is no "Frobenius at N".
  - Eisenstein coefficients factor N deterministically.
  - Edixhoven–Couveignes Remark 2.2.4 already states the prime-versus-composite dichotomy for σ_{k−1}.

## Genuinely open native questions (none NP-hard)

1. **A single coefficient a_f(N) of a fixed non-CM cusp form.** Equivalently, a fixed cuspidal theta difference,
   (r_{L1} − r_{L2})/4 at level 11, or a 37a1 analogue. Compute it at a balanced N below L[1/3] without separating out
   σ(N). It is in ZPP^FACT, and NP-hardness would give NP = coNP under GRH. No prior work shows it factors N. [nativemod]
2. **Möbius-weighted summatory functions** (Φ, twisted totient) faster than x^{2/3}. [nativeinv S1]
3. **Representations of φ(N) mod 3 by digit automata on semiprimes, and torsion of poly(n)-size Z[N]-matrices
   equal to φ(N).** Partial negatives are proved. [nativeinv S2, S3]
4. **A barrier-side slack:** a finish faster than √Δ given h to precision Δ. A Sutherland-type Δ^{1/3} generic
   bound would leave room down to N^{1/7} in the LIN model. [nativeana]

The subjective probability, per the agents, that any of these yields sub-L[1/3] factoring is below 1%.
