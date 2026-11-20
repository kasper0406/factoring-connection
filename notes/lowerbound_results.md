# lowerbound: a unified relation-model theorem "α ≥ 1/3" (2026-11-20)

**Question (user).** Is there literature proving an L[1/3] lower bound for factoring? There is none in any general
model. Factoring ∈ NP∩coNP, so an unconditional superpolynomial classical lower bound would separate P from NP∩coNP.
What exists:
- random-squares thresholds, at L[1/2] for uniform candidates: Pomerance ICM 1994; Croot–Granville–Pemantle–Tetali
  2012; Balister–Bollobás–Morris 2018;
- Pomerance's size principle;
- Zhu–Lv–Liu, DCC 94:51 (2025/26). It shows c ≥ (32/9)^{1/3} ≈ 1.5263 *within the NFS/TNFS/MNFS cost formulae* (an
  elementary closed form, with a pairwise-resultant star constraint for MNFS). That is a parametric bound, not an
  algorithmic one;
- generic-group and generic-ring bounds (Shoup, Maurer–Wolf, Damgård–Koprowski, Aggarwal–Maurer, Brown), which
  cannot express sieving.

The repo's barrier map had no compositional or adaptive theorem (completeness audit). So we tried to prove one.

**Process.**
- C, GPT-6.2-astra xhigh (codex session 01a1fa34-0b6c-7341-960b-40fb5b993e9c): primary prover. Round 1
  `codex/lowerbound_1_codex_proof.md`; round 2 cross-review and merged theorem `codex/lowerbound_2_codex_merged.md`.
- G, Gemini 4 Carbon High (agy conversation 7e0ef22c): independent model M_G, literature and adversaries.
  `codex/lowerbound_1_gemini{,2}.out` and `codex/lowerbound_2_gemini{,2}.out`.
- Read-only Opus critic on the merged theorem. Briefs: `codex/lowerbound_common.md`, `codex/lowerbound_1_*.md`,
  `codex/lowerbound_2.md`.

## Result (MODEL theorem; Opus critic: ACCEPT WITH FIXES)

**Model (priced valuation-probe model).** This is the smooth-number analogue of the generic group model.
- **Computation.** Uniform and charged, replayable from public inputs, coins and paid answers. Only
  relation-certified splits count.
- **Probes.** Smoothness information comes only from paid probes on fresh residual cofactors.
- **H, conditional Dickman.** Every fresh probe succeeds with probability ≤ ρ(s/log B)^{1−o(1)}, given everything
  computed so far, internal state included.
- **Blocks.** These are target, column set, bound B = e^μ and quota h. A block opens before its first credited probe
  and its exceptions (≤ ⌊h/4⌋) are fixed at opening. Each revelation credits one address once. Blocks that share
  columns are pooled, with one exception budget.
- **R, quota.** Extraction needs h distinct fresh successes in some block, with log h ≥ (1−o(1))μ.
- **Credit (after the critic's fix).** Only N-independent precomputation shared across targets may go uncredited.
  Any N-dependent stage (for example a rational sieve with m depending on N) is credited at the address's first touch,
  and success means the joint relation. The critic's patch keeps the supermartingale valid. Reuse budget:
  Σ_i Q_i ≤ e^{o(F)} W.

**Theorem (C, merged).** Suppose every ordinary block profile of credited residual log-sizes satisfies
λ(k) ≥ A√(x log(k+1)) through rank t = e^{cF} (F = x^{1/3}(log x)^{2/3}). Then solving a positive fraction of K
targets costs

  **E[W]/K ≥ L_N[1/3, (4A/3)^{2/3}]**,

with or without an optional single-witness clause, which gives (κ/18)^{1/3}. The theorem allows:
- adaptive block opening and stopping;
- unions and pooling of families;
- batches (targets may be dependent);
- O(log x) advice bits;
- descent under tree incidence (Assumption T).

**Proof.** In each block, exp(θY − (e^θ−1)P) is a nonnegative supermartingale from opening, with the predictable sum
P bounded on every path by the sorted profile sum M. Stopping at d = ⌈3h/4⌉ and applying Markov gives (eM/d)^d. A
union over ≤ t openings follows. The size–entropy calculus is min_μ[μ − log ρ(s/μ)] = (1+o(1))√(2s log s); at
s = A√(x z), z = cF, the threshold is c_A = (4A/3)^{2/3}. The critic verified every step, and the numbers 1.2114,
1.5263, 1.9230 and 0.39894.

## What is proved vs assumed (honest headline, per the critic)

- **Proved profile, all N, no counting heuristic.** Rank-2 binary-form methods whose jointly credited product form has
  N^g | disc (Lemma C2 + Bean/Davenport + Lemma FT; pooled families need log J ≤ η log k), under the explicit
  error-budget condition C_FT·D log(2E) ≤ h/4 (critic MAJOR 2). This gives A = √g:
  - **L[1/3, 1.2114]** with ramified extraction (g = 1);
  - **L[1/3, 1.5263]** with two ring maps (g = 2), which SNFS attains.
- **Not covered by the proved inputs.** Small-quota blocks (e.g. B = x^{1/2}, where D log E ~ x^{2/3} ≫ h) and huge
  boxes (D log E up to B/16).
- **Heuristic profiles only.**
  - GNFS's A = 2, i.e. 1.9230, needs GAE plus typical norm sizes.
  - MNFS 1.9019 and the factory 1.6386 are reproduced as optima of the classical cost programs, not as lower bounds.
    With proper credit they are lower-bounded only at ~1.5263, and that heuristically.
  - Free linear algebra (ω = 0) gives 1.859536 for MNFS and 1.636264 for the factory (both consultants agree).
  - QS/MPQS (exponent 1/2) and towers also rest on typical profiles.
- **KH subclass.** A θ-regular KH subclass gives A = θ/√7 (c ≈ 0.399 at θ = 1/2). It is a checkable sufficient
  condition that no algorithm has yet been shown to satisfy.

**Novelty: real but modest, mostly formalization.** The constant law is Theorem SS(d), i.e. Pomerance's size
principle. What is new:
- adaptive multi-block composition by a stopped exponential supermartingale with a sorted-profile predictable bound;
- explicit credit, exception and reuse accounting;
- all-N rank-2 profile inputs that need no polynomial-counting heuristic.

## What did not survive, and side results

- **G's removal of Assumption T** (analyst cutoff B_* = T²(log N)^C, last-step argument) **fails as stated**.
  - A forced prime can persist in the residual: rows with q and q³ give q² in the last residual.
  - C gave an explicit binary form example: f = X − Y, g = X² + (N−1−q)XY + qY², Res = N, with residual 1 at (0,1)
    on the index-q lattice.
  - Repairs need full-cofactor SM2 plus a non-cusp hypothesis H_cusp on single descent lattices. So general DAG
    descent still needs a separately proved quota or witness.
- **Batch amortization (G's adversary A1).** Unrestricted dense batches really do admit polylog amortized factoring:
  sieving least prime factors up to 2N with K ≈ N/polylog. G's relation-only factory version still needs
  smooth-difference distribution and fruitfulness estimates (C).
  - Inside the pre-registered credit model, no numeric cap is necessary.
  - For target-blind generation, the caps are about K ≤ N^{1−ε} (pure relations) and K ≤ N^{1/2−ε} (with batch
    remainder-tree gcds).
- **G's Axiom NP** (checkable N-publicity) is useful but needs gcd(disc W, N) ∈ {1, N} added. A finite gcd checklist
  cannot detect coefficient-hidden advice such as p+1. The real exclusion of factoring in disguise is
  "relation-certified success only" plus H given the full internal state (critic).
- **MINOR (critic).**
  - R must count rows by distinct fresh revelations, since non-primitive multiples add rows without revelations.
  - Joint-tuple smoothness needs u₁log u₁ + u₂log u₂ ≥ (u₁+u₂)log((u₁+u₂)/2).
  - Blocks must be column-disjoint unless pooled.

## Where sub-L[1/3] must break the model (the useful negative information)

Any sub-1/3 algorithm must violate at least one of these:
1. **H.** It needs a constructor or decoder that selects smooth residuals above the Dickman rate: Prop SI-type
   constructors, or a beyond-Johnson CRT list decoder (M2) on N-public residue tables. Excluding these is a method-class
   restriction, like generic-group models, not an arithmetic theorem. A proof that no such decoder exists would be an
   arithmetic-circuit lower bound.
2. **The profile.** It needs exceptional small norms beyond the error budget at rank ≥ 3: LP_exc^rel, or multi-cone
   unit cancellation in type-D modules. Theorem KH bounds only generator heights.
3. **R or reuse.** It needs thin dependencies (PIN-type pinned coordinates), heavily shared descent DAGs, or one
   witness serving many related inputs.

These are exactly the three walls the subthird/equivtarget threads ended on. The theorem turns the barrier map from a
patchwork into one statement with three named escape hatches.
