# cogentic: Google Research "Cogentic" papers vs this repo (2026-11-25)

**Source.** https://github.com/google-research/google-research/tree/master/cogentic_research_papers (posted
2026-11-25). Produced by Cogentic (arXiv:2609.40324), an agentic framework with Gemini 4 as its base model, "with no
human interaction other than the initial prompt". Papers:
- integer_factorization: 54 pp.
- deterministic_factoring_l34: 59 pp.
- discrete_logarithm: 36 pp.
- Plus game theory and learning theory papers, not relevant here.

READMEs are copied to `references/cogentic/`. Comparison by Gemini 4 Carbon (agy 487bb700):
`codex/cogentic_1_gemini{,2}.out`. Brief: `codex/cogentic_1.md`. Not verified by us.

## Claims

1. **Unconditional Las Vegas complete factorization in L_N[1/3, 105]**, with workspace L_N[1/3, 40] and AKS-certified
   output. A Lean 4 formalization is promised "in a companion submission".
   - **Fruitfulness** comes from sealed-root Dixon descent: the same architecture as our Theorem DD / U-MF chain.
   - **Relation supply at 1/3:**
     - two-field rational elimination;
     - Bernoulli acceptance proportional to smooth weight, with backward path telescoping, so only ≤ 0.65 ln ln N ≪ d
       nodes per path need joint control (this sidesteps our nested Lemma OD and Prop BT);
     - a mod-N discrete Laplace coefficient law;
     - a Siegel-free hybrid smooth weight (Pólya–Vinogradov plus Brun);
     - an unconditional LV6′ via large character primes and 3-adic Hensel roots.
   - **Where the constant 105 comes from:** c > 21.15/(1 − 2^{−1/3}) ≈ 102.5. This is the same strong-induction recursion
     blow-up as in our Theorem NR / U-MF.
2. **Unconditional deterministic factoring in O(exp(3(ln N)^{3/4})) = L_N[3/4].** Also: deterministic splitting given
   a multiple of all p−1 in exp(O(√(ln N) ln ln N))·log M̃, and deterministic quadratic nonresidues mod p in L_p[1/2].
   - It implements our Z-route (XK-Z-type elimination over totally real base-m fields) and **bypasses DSR and ERH**
     with three devices:
     - a long-progression perturbation sweep wider than the modulus;
     - Landau–Page repulsion plus Montgomery zero density;
     - subgroup-index saturation.
   - The QNR part uses a multiquadratic Walsh–Hadamard / biquadratic cocycle positivity against a 16th-moment Burgess
     bound. This is the TIED-PARITY / RM(1,K) structure from our hatch round.
3. **Fixed-characteristic F_{p^n} discrete log** in exp((2+o(1))(ln n)²/ln ln ln n), i.e. sub-quasi-polynomial.

## Effect on the repo (if correct)

**Superseded.**
- Theorem A (L[1/2, 0.919]), U-MF (L[2/5]), 2L-final / 2L-MF (L[3/8]) and conditional Theorem NR are all superseded by
  the unconditional L[1/3, 105].
- Worst-case deterministic factoring improves from N^{1/5} to L[3/4].
- Prop PIN / Miller splitting loses ERH.

**Not superseded.**
- Theorem AV: average-case deterministic factoring still beats L[3/4].
- Deterministic L[1/2] (DSR_{1/2}) stays open.
- Theorem DL / dlog25 (prime-field DL at L_p[2/5]) is not superseded by paper 3, which is fixed-characteristic only.
  Paper 1 plus Prop XK might lift prime-field DL to L_p[1/3]; that is unchecked.

**Consistent.**
- The lowerbound model theorem (α ≥ 1/3; floors 1.21/1.53 ≪ 105) and the GRH assessments.
- No proved repo barrier is contradicted:
  - dsr's Thm SB, Thm FN and BOX-∞ concern prescribed targets, and paper 2 sweeps instead.
  - Paper 1 respects nested Prop R2 and escapeR2's Lemma TI / Cor EX.
  - Paper 1 does not prove Conj F for unseeded NFS: every row descends from a blinded random square.

## Sub-L[1/3] relevance (the only goal here)

**None.**
- Paper 1 stays inside all three hatches:
  - its hybrid weight has Dickman-rate density, so H is not violated;
  - its candidates are rank-2 binary forms, so the profile is proved;
  - it uses L[1/3, 18] independent blinded trees with no sharing, so R is intact.
- Paper 2's cocycle device needs a known multiple of λ(N), which is factoring-equivalent, and caps at L[1/2] scale.
- Paper 3 is blocked by the characteristic obstruction: no ring map Z/N → F_p-algebra.

## Credibility (Gemini's read)

- The papers are unusually complete and internally consistent. Unlike the Masullo draft, they have no stubs. Constants
  are transparent, and edge cases are handled (skewness excision, index primes, semiprime active characters).
- **Lean claim: be skeptical until the code is released.** It appears only in the READMEs, not in the papers. Lean
  versions disagree (v4.35.0 vs v4.33.0), and the needed inputs (AKS, effective Chebotarev, Stark, Landau–Page,
  Montgomery zero density, Turing-machine cost) are not in Mathlib.
- **Referee priorities:**
  - P1 Lemma 4.6 + 4.4: measurability in backward path telescoping; whether capped Bernoulli loops preserve
    (2/μ_t)·X_a.
  - P1 Lemma 3.5(4)(a) / 3.6(4): gcd absorption in the hybrid weight; the bad-conductor partition.
  - P2 Lemma 5.1 / 6.3 (semiprime q, Type III 16th-moment fibers) and Prop 9.3 / 10.1 (three-piece split, stopping-time
    peeling).
