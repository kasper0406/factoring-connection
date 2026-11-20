# multcopp: can the RH-remark's "multiplicity over Z" become a Coppersmith-style factoring algorithm? (2026-11-20)

Question (from grhremark): the remark (alpo.ge/remark.pdf, claimed Lean proof of RH + Dirichlet GRH) says its key new
idea is Stepanov's method "with multiplicity" over Z. Coppersmith/Howgrave-Graham is also a Stepanov-type argument:
auxiliary polynomials vanishing to high p-adic order at a hidden root. Could the remark's mechanism beat the
known-bits barrier (topblock Prop T) and so give factoring below L[1/3]?

Consultants (independent): GPT-6.2-astra xhigh via codex (`codex/multcopp_1.md` → `codex/multcopp_1_codex.out`,
session 01a1f9dc-2f30-7761-9b0c-dccaef23598b) and Gemini 4 Carbon High via agy (`codex/multcopp_1_gemini{,2}.out`,
conversation 789d428f). Both read the Lean sources. Source: the 14-part tarball, sha256
9631b4f020199e580ba86178cb56c6e9b09f690c2369a80dc4334aa93f49c7fc, 175k .lean files. It was unpacked to a scratch
directory and **not built**. Toy tests: `experiments/multcopp/{codex,gemini}/` (both pre-registered).

## Verdict

**Clean kill. No sub-L[1/3] algorithm, and nothing even below N^{1/4}.** Both consultants agree on the mechanism and
on the obstruction. Estimates: P < 1% (codex), < 0.1% (Gemini).

## What the mechanism is (Lean inspection; mathematical correctness of the Lean is UNVERIFIED)

- **Top statements.** `Main.lean` / `Statements.lean` prove:
  - Mathlib `RiemannHypothesis`;
  - Dirichlet GRH for every modulus and every χ, imprimitive included (`DirichletLZerosOnLine`);
  - the primitive version.

  There is no Dedekind or Hecke statement. So, if correct, it discharges exactly the "ERH" of
  detGRH/zroute/dsr, and not nested's X_ERH or dBPW (see grhremark).
- **Target.** The proof bounds the Cesàro von Mangoldt remainder G_z(X) = Σ_{n≤X} Λ(n)(X/n − 1) − X log X +
  (1+γ)X (`Formalism/StatementOC.lean`).
- **p-adic multiplicity.** It comes from derivatives of theta/CM functions on a CM curve at primes p inert in
  K = Q(√−d). There the reduction is supersingular: the Hasse invariant vanishes and D^p ≡ 0 (mod p) for the
  invariant derivation D.
  - Differentiating a times gives valuation Σ_{i odd} ⌊a/p^i⌋ − O(log_p a) ≈ ap/(p²−1).
  - The PDF's "a_p + a_p³ + ⋯" is a text-extraction artifact for a/p + a/p³ + ⋯.
  - Lean: `Formalism/SharpRestriction/Divisibility/CartierLocalBounds.lean` (rho, gam), `Formalism/GaussJetBounds`,
    `Formalism/FormalGroupLog`, `QuadraticExponentIdentity.lean`.
- **The √X chain.** Iterated theta lifts between definite unitary groups,
  `Development/AbelianVarietyTheory/BlockTheory/TwistedDiscreteSeries/IteratedLiftChains.lean`, `ChainParSucc.lean`.
  - Rallis/seesaw make the Gram determinant factor over the m lifts (`Formalism/ThetaLiftIndependence/Statement.lean`:
    DetChain, ChainScalar), so the p-adic gains add over the chain.
  - The archimedean side is paid once, through Bost's slope inequality (`Formalism/LatticeSlopes`) and Chowla–Selberg.
  - m ≍ √X comes from a Chebyshev-node quadrature error O(X/m) (`Formalism/ZetaNonvanishingChebyshevWeights.lean`).
- **The "stop at 2/3 a priori".** An error balance X/m + Xm/k + m + k/n² gives X^{2/3} without twist-family
  "division", and X^{1/2} with an Iwasawa/Behrend family. It is a polynomial exponent in X, unrelated to L-notation
  (Gemini, from `Formalism/SizeWindowBrackets`, `SChain.lean`; codex's reading is similar but labelled speculative).
- **Size of the proof's objects** (Gemini, `CharacterBoxes/…`, `SChain.lean`): m ≈ X^{1/2}, unitary rank ≈ X^{400},
  genus size above exp(X^{402}). This is a proof device only.
- **Stepanov dictionary.** The "points" are all inert primes p ≤ X, weighted by log p. The "auxiliary polynomial" is
  the period/evaluation map between Hermitian lattices of algebraic modular forms. The "degree bound" is the slope
  inequality.

## Why it is not a Coppersmith algorithm

1. **Blind hidden-prime vanishing really happens.** Pick d with (−d/N) = −1. Then one factor p is inert
   (supersingular) and the other split (ordinary), so p | D^a θ(0) for all a ≥ p while the q-part is generically
   nonzero. Hence gcd(E_a mod N, N) = p splits N outright, computed mod N without knowing p.
   - Codex's toy, on y² = x³ − x with 78 semiprimes in [67, 127]: the gcd is 1 at order 31 and equals the inert
     part at order 127, matching the predictions in 78/78.
   - Gemini's toy, on w² = 1 − s⁴ with 30 semiprimes of 14–20 bits: at a ≈ N^{1/4} < p it split 0/30; above √(2N)
     it split 30/30.
   - This is a whole-factor "leak", so the repo's one-bit-leak heuristic (caudit §7 Q4) would misclassify it.
     Its barrier is cost, not information.
2. **The multiplicity is zero below a = p.** Coppersmith gets p^m at degree O(m) because N and x + p₀ are built
   into the ideal (N, f)^m. Here the valuation is ⌊a/p⌋ + ⋯ = 0 for every a < p, in every lift (m·0 = 0). For
   balanced N, informative orders need a ≥ p ≈ √N.
3. **No shortcut is known.**
   - Write D^a x = y^{a mod 2} P_a(x). P_a has leading coefficient (a+1)! (codex, by induction), so the order-a
     coefficient is a unit mod p while a+1 < p.
   - A fast evaluation of that coefficient mod N would binary-search the least prime factor. That is fast modular
     factorials, the Pollard–Strassen / Bostan–Gaudry–Schost route, whose best known cost is N^{1/4}.
   - Fast-forwarding D^a mod a composite N fails: D^ℓ = A_ℓ·D (mod ℓ) uses binom(ℓ,k) ≡ 0 (mod ℓ), which is false
     mod N. This is the Theorem PT / FROB pattern.
4. **Slope ledger at a single prime.** The remark's p-adic gain sums over all ~X/log X inert primes ≤ X, against an
   archimedean height ~ m·b·X log X. One hidden prime contributes ~(log p)/(p log X) of it, which is negligible.
   Injecting N as a level blows the automorphic dimension up to N^{Θ(rank)}. Rewriting Howgrave-Graham through the
   slope inequality is an identity and reproduces the standard bound ξ < 2βθ − θ² ≤ β².
5. **Subexponential X.** Taking X = L_N[β] reaches only L[β]-smooth primes and misses p. This is the size–entropy
   law, Thm SS.

**Single missing lemma** (codex): compute the leading coefficient of P_a mod N (equivalently (a+1)! mod N) for binary
a ≤ √N in L_N[α < 1/3] time. That is a fast-factorial breakthrough and owes nothing to the remark.
