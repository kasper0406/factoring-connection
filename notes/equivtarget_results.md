# equivtarget: "factoring-equivalent ≠ hard": FACT/DL-equivalent primitives as targets (2026-11-20)

**User's point.** "This primitive is factoring-equivalent (or FACT ∧ DL, or solves F_p-DL)" was often used as a kill.
That is only a kill if factoring is hard, which is exactly what we are betting against. An equivalent problem might
have structure, or a known easy regime, that factoring lacks.

Phase 29 (notes/native_results.md) already did this for 55 *invariants*. This note covers the *constructive* and
DL-type targets dropped in phases 18–32.

**Process.**
- Briefs: `codex/equivtarget_common.md`, `codex/equivtarget_1_{codex,gemini}.md`, `codex/equivtarget_2.md`.
- C, GPT-6.2-astra xhigh (codex session 01a1f9ff-4307-7392-8985-9931919161ae): systematic harvest and re-scoring.
  Full report `codex/equivtarget_1_codex_harvest.md`, with a 17-row reduction table and a 537-line concordance over
  77 notes.
- G, Gemini 4 Carbon High (agy conversation 330807ac): transfers from known easy regimes.
- One cross-review round. Toys (pre-registered) in `experiments/equivtarget/`.

## Verdict

**The reframing is right in principle and was worth doing, but it gives no sub-L[1/3] route.** For each target we
asked whether it has an easy regime we can steer N's instances into. Every candidate fails at a *quantitative
conservation law*, not at "≡ FACT". It either collapses to hint generation carrying more than log₂N bits about p, or
lands on one of the two bare walls already isolated in subthird:
- beyond-Johnson 2-D quadratic CRT list decoding (M2);
- multi-cone unit cancellation (LP_exc^rel).

Combined P(α < 1/3) via this direction: < 0.3% (both consultants).

## Why the DL-world easy regimes do not transfer

- **Characteristic obstruction (proved; C, sharpened by G).** Any ring map φ: Z/N → A into an F_ℓ-algebra, unital or
  not, with gcd(N, ℓ) = 1 is zero: Nφ(1) = ℓφ(1) = 0, and Bézout gives φ(1) = 0. Small-characteristic quasi-poly DL
  (BGJT / Kleinjung–Wesolowski), Joux's medium-prime methods, MOV/Frey–Rück, GHS and Weil restriction all preserve
  characteristic. So none of them can receive a factoring instance without a factor.
- **Non-ring transfers exist, but each pays separately (G).**
  - Level-N moduli in characteristic 2 (#X_1(N)(F_2) = φ(N)/2) have model size Θ(N).
  - Group embeddings (Z/N)^× → F_{2^k}^× need λ(N) | 2^k − 1, i.e. k ~ ord ~ N.
  - Positional lifts hit binary carries.
  - Weil restriction stays in characteristic p×q and needs a poly(p) factor base.
  - The Jacobi symbol is a free group map to F_3^× but carries no factor (C).
- **The poly(p) trap (G).** Every sub-L[1/3] or Weil-descent method in characteristic p uses a factor base of size
  ≳ p ~ √N.
- **The ECM–Dickman scissors for BGJT-style descent (G).** Splitting a mod N into L[β]-smooth pieces costs
  max(1−β, β/2) ≥ 1/3.

## Candidates and how they died

| # | Target (equivalence) | Easy regime tried | Outcome |
|---|---|---|---|
| C1 | Smooth ideal in a specified class of Cl(−4kN) (class root ⇔ FACT under GRH) | Steering by public multipliers; wide boxes | **Proved under GRH:** small boxes cover only 1/L[1−a] of classes. Every sublattice gives a cofactor form of the same discriminant, so lattice/slope calls cost ≥ L[max(a, 1−a)] ≥ L[1/2] (G; C agrees). Wide-box search = class-group lattice sieve. Toy 476 vs 70 smooth norms is Knuth–Schroeppel multiplier scoring, a constant-level effect, before selection cost. Survives only as the M2 beyond-Johnson decoder (0.1–0.2%) |
| C2 | One 2-isogeny neighbour mod N from j alone (⇔ FACT: masked Legendre parameter gives a split with probability exactly 4/9, verified 62,752/141,192) | Volcanoes, Kani/SIDH, Frobenius | Root of Φ₂(j, X) mod N = square root + cube root mod N (Cardano): Rabin/RSA-type. Labelled torsion or unique common neighbours would help, but none is constructible from j. < 0.05% |
| C3 | Anomalous elliptic component (trace one at a hidden prime; Cheng's 4p − 1 = Ds²) | Smart/SSSA p-adic log | Real special-family factorer, but D ≤ L[a] covers p^{−1/2+o(1)} of primes (2/1982 near 2²⁴). The genus ≥ 2 variant is worse or unproved. < 0.05% |
| G1 | "CRT-HNP": residues x mod N small mod p and mod q | Multi-sample HNP/ACD lattices | Lattice = the classical approximate-common-divisor simultaneous-approximation lattice (Howgrave-Graham, Cohn–Heninger, GGM 2016). The toy beat σ = 1/4 (20/20 at 64 bits), but the generator must supply (m+1)δ > log₂N total bits about p. Natural N-only biased families such as x_i = (N−1)/d + i are public shifts of one value (C's toy) and give nothing. ≤ 0.1% |
| G2 | Norm-N CM endomorphism ν = π_pπ_q of y² = x³ + 1 (⇔ FACT), via Kani 2^e-gluing over F_{ℓ²} | SIDH break (Castryck–Decru, Robert) | No torsion images are known: free data fix ν mod 3 only (O(1) bits), leaving a torus of ~N candidates (Θ(2^e) scan). C gives an explicit counterexample to G's Kani "iff" at N = 403 (two non-family kernels with F†F = 2^{19}I). Kani adds degree, not information. < 0.05% |
| G3 | LP_exc^rel at unit rank ≥ 2 (Joux-tower / Q(ζ₇)⁺) | Incomparable unit cones S^{e−2} | 12-cone LLL: 0 of 214 full-support vectors below N² (min 3.345). Units preserve norms and valuations, so new parity rank needs inequivalent elements. Subsumed by subthird's LP_exc^rel (0.2%) |

**Also re-scored (C's table, in `codex/equivtarget_1_codex_harvest.md`).** Point counts and orders, Prop SI / promise
CVP, Schnorr–Adleman CVP, hidden-root samplers, NAV / unit-image matching, composite Frobenius, DSR / pinned rows,
period / O1 observables, prime-vector ADD / successor / Zech, regulator multiples and genus membership, factor-moving
products and CRT phase, Chebyshev compositions and S-unit bases, factorial circuits and PPA/PWPP/MCSP. **Zech
addition** is DL-equivalent in its own field, and fixed characteristic makes it genuinely easy (quasi-polynomial
DL). But no map takes a general factoring instance there while keeping the information needed to split N.

## If continuing (both rate it ~0.2%)

- **C.** An output-sensitive quadratic CRT decoder for M2 on blinded, near-uniform target classes. Charge multiplier
  selection, partial assignments visited, and work per independent relation. Cap the first comparison with direct
  sieving at one CPU-hour on 64–96-bit discriminants.
- **G.** Multiplicative S-unit product steering (beyond linear Minkowski LLL) for LP_exc^rel at (e,n,m) = (3,4,3) over
  Q(ζ₇)⁺, about 2 h of Python.
