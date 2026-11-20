# hatch: attacks on the three escape hatches of the α ≥ 1/3 model theorem (2026-11-20)

**Setup.** notes/lowerbound_results.md says any sub-L[1/3] algorithm must break **H** (Dickman-rate smoothness
probes), **the profile** (LP_exc^rel: exceptional small norms at rank ≥ 3), or **R/reuse** (quota of fresh relations;
shared descent; related inputs).

**Runs (rules `codex/hatch_common.md`).** Each consultant then broke the other's results (`codex/hatch_x.md`).

| Hatch | Consultant | Output and artifacts |
|---|---|---|
| H1 | C, codex | `codex/hatch1_codex.out`, `experiments/hatch/hatch_h1c/` (includes the cross-review `cross_review.md`) |
| H1 | G, Gemini | `codex/hatch1_gemini.out`, `experiments/hatch/hatch_h1g/` |
| H2 | C, codex | `codex/hatch2_codex.out`, `experiments/hatch/hatch_h2c/` |
| H3 | G, Gemini | `codex/hatch3_gemini{,2}.out`, `experiments/hatch/hatch_h3g/` |

Cross-review outputs: `codex/hatch_x_{codex,gemini}.out`. Also included here: the follow-up asking whether
Masullo's orbit counting extends to LP_exc^rel (`codex/masullo_ext_{1,2}.md`, outputs `codex/masullo_ext_*`,
artifacts `experiments/masullo_ext/`).

## Verdict

**No hatch opened. Several new obstructions were proved, and each hatch is now narrower.** Every surviving opening is
an explicit lemma (listed at the end).

## Accepted as proved (after cross-review; internal, unrefereed)

### Hatch 1 (H)

- **Weighted Theorem J (C).** With separate nonnegative weights at every prime-power depth q^k, a profile fixed in
  advance yields at most ρ(log X/log B)^{1−o(1)} smooth witnesses per call on average, even granting every weighted
  Johnson certificate at unit cost. So multiplicity weighting gains nothing at leading order. Adaptive weights are
  not covered.
- **Galois folding fails at split primes (C; G Thm 2b).** A primitive norm vanishes on one conjugate branch only, so
  full-orbit agreement is identically 0.
- **Arithmetic-derivative analogues (G Thms 1a–1d, as corrected by C).** Buium p-derivations, Int(Z), hidden
  derivatives and fixed-order jets do not give a usable integer "multiplicity code".
  - 1a holds: the jet congruences are equivalent to q^m | X.
  - 1b holds only as distinctness: there is no typical-height claim.
  - 1c is proved at fixed order.
  - 1d's Int(Z) index and covolume bounds are correct, but its "universal Johnson wall" is **false**. C gives a
    counterexample: (Y−1)(Y−2) ∈ I* with a small sup norm.
  - Root cause: Z has no derivation and no nontrivial automorphism.
- **Folded O_K CRT codes decode past Johnson when errors respect whole prime orbits (G Thm 2a; affine-subspace
  containment).** A blind replay gave 30/30. C shows that a coordinate-only argument gives the same leading threshold
  for a fixed field, so the Galois action is not needed, and novelty is limited. It is not useful for factoring,
  because of 2b.
- **Measured (C).** Quadratic rational reconstruction recovers witnesses faster *once a divisor is chosen*: 5.96 vs
  10.27 ms per witness at 96 bits. Choosing the branch before enumerating the L[1−β] divisor products is the open step.

### Hatch 2 (profile / LP_exc^rel)

- **Determinant obstruction (C; checked by G).** Take a bounded projective height cutoff u. A candidate whose norm is
  far below the Theorem KH envelope must lie in one of at most D specified proper k-subspaces, or cancel at
  Ω(√x/log x) embeddings. At the exact L[1/4] balance this becomes Θ(D) embeddings.
  - Fixes: allow unequal relative degrees n₁, n₂.
  - It fails at unbounded unit height (h ~ x^{1/4}). Parity rank *inside* the subspaces is not charged.
- **Quantitative Subspace Theorem is too weak (C).** Evertse's cover prefactor has log Θ(r² log D) = x^{1/2+o(1)},
  above budget, and its hypotheses exclude the base-field action. Evertse–Ferretti gives O(m + log log D) per weighted
  frame, but choosing frames costs binom(2n, m)^e.
- **Theorem Rel-BT (G, corrected by C).** Fixed base k, Δ_rel ≠ 0 (relative general position).
  - (i) v_{𝔭}(Δ_rel) ≥ m(m−1)(2n−2)!/(2n−m)! at each 𝔭 | N (relative ramification indices). Hence
    |N Δ_rel| ≥ N^{P(m,2n)}. This answers grdeg's "relative R1" question when Δ_rel ≠ 0, non-power generators included.
  - (ii) Rel-R4: an envelope bound for m independent candidates.
  - (iii) Corrected:
    - The number of unit cones is O_k((1+S)^{r_k}) with S = H/(2n) + e log(C_k m E). G's individual-unit bound is
      false, because generators can be rescaled by units.
    - The per-cone volume is exp(O_k(m log m))·N^{−m(m−1)/(q(q−1))}·Y^{m/q}, with q = 2n. For complex places this
      needs C's spherical-Hölder argument, since the Thunder citation does not prove it.
  - (iv) The boundary error is r·D·(2E+1)^{r−1} per cone.
  - Scope: it closes door (i) at α ≥ 1/3 only inside the box scope H5_r. At the yield-relevant ternary box (E = 1) the
    boundary error exceeds 3^r, so **the L[1/4] regime stays open**.
  - **Masullo's App. A.2 is not needed.** Standard Dirichlet unit theory over fixed k suffices.
- **Experiments.**
  - C: S-unit / multi-cone steering at (e,n,m) = (3,4,3) over Q(ζ₇)⁺ and (3,5,3) at 60 bits gave 0 full-support hits
    in 80,359 tested. An exhaustive audit of 237,608 grid tuples found a best log|F|/log N of 3.47 (44-bit) and 3.01
    (60-bit), never below 2.
  - C, Masullo-extension round: an embedding-threshold cover of exceptional subspaces, tested over Q(√2) at n = 8,
    m = 3. **One 72-bit divided-difference witness escapes it.** Its norms are 108241 and 3175924591, and it raises
    parity rank from 47 to 48 even after all smaller-norm rows in the box are included. One row is not an
    asymptotic yield.

### Hatch 3 (R / reuse)

- **MIN-ROW (G).** Algebraically, one fruitful square row suffices, so "rows > π(B)" is not a linear-algebra
  necessity.
- **PIN-1 (G, repaired).** Miller's chain on a single base splits N given any multiple of its order, with probability
  ≥ 1 − 2^{1−ω} for a uniform base. The deterministic Jacobi-(−1) version needs even ω and odd multiplicities; C's
  counterexample is N = 231 with base −1.
- **PIN-STRIP exponent (corrected).** Use cutoff 13cF_{1/3}, not the original.
- **The Z-route needs no fewer rows.** Pinning one coordinate still needed more than π(B)/2 rows in 20/20 toys.
- **THIN-LEDGER (heuristic).** A genuinely thin dependency family of h = L[η] relations would give L[max(1/5, η)]
  RAM factoring. So R is load-bearing, and the whole hatch reduces to whether thin dependencies exist.
- **DAG-POOL, arithmetic part only.** Sublattice discriminants scale as Q^{(D−1)(D−2)}. There are at most D^{ω(Q)}
  sublattices per index, and the pooled additive profile and cusp budget hold under a degree/size restriction.
  **But it does NOT remove Assumption T.** It presupposes a quota of *fresh* smooth residuals. C's counterexample:
  rows q³b_i with one eliminator qb_0 give independent base rows, yet every stripped parent keeps q².
- **TIED-PARITY (G).** Prime-ideal parity indicators of the degree-p₀ subfields of a (Z/p₀)^r tower span rank p₀^r
  (no sharing), except at p₀ = 2 over F₂. There they form the Reed–Muller code RM(1,r), compressing 2^r columns to
  r+1.
  - The candidate built on it, C-RM, needs ML-RM: small lifts and norms across r coordinates with independent-row
    yield. That is the LP_exc^rel obstacle again, plus class-group 2-torsion and squareclass independence. Speculative.
- **SINGLE-N-REUSE (accounting).** Nine related-input schemes conserve cost for a single N: multipliers in QS/NFS,
  shifted bases and factories, SNFS multipliers, powers N^j, product bundling, dense multiples of N, polynomial maps,
  and exact algebraic M_k.
- **Corrected claims.**
  - Large-prime-cycle thresholds (CGPT/BBM) do **not** follow from marginal prime incidences in this model (C).
  - PLANT's size floor is false in general (4² − 15 = 1).

## Masullo extension (did the orbit counting extend?): no

Masullo's spectral and ray-class machinery counts full-rank objects (ideals in ray classes). The LP_exc^rel family is
a low-rank linear section, a measure-zero subvariety of the torus R_{K/k}G_m, of log-general type since 2n ≫ m.
Bridges tried:
- **B1, Batyrev–Tschinkel torus harmonic analysis.** Needs a fixed toric compactification. The section imposes 2n−m
  additive conditions, which can be inserted exactly via additive characters, but there is no uniform estimate for
  them.
- **B2, one-prime-at-a-time chopping.** The cutoffs are degree-uniform, but it needs closure under π_v-scaling, which
  the module lacks, and it cannot see joint prime incidences of one candidate.
- **B3, ray-class sandwich over KH's degree-T lattice.** Even with Masullo's GRH (square-root errors), the error
  √Y = L[1/2] (or Y^{T/2}) exceeds the 3^r = L[1/4] box. Full-field units also erase the mod-N ray condition.
- **B4, periodic orbits and shells.** These are single prime loops. Base-unit deduplication changes log counts by
  O_k(1).

The only transferable piece turned out to be the standard unit fundamental domain over fixed k (Rel-BT). That is
classical, not Masullo-specific.

## Open lemmas: where the hatches now stand

1. **H, branch selection.** Find, in poly(π(B)) work, a divisor branch whose quadratic reconstruction yields a smooth
   witness, *before* enumerating L[1−β] divisor products.
2. **Profile, growing-rank threshold lemma (C).** For reduced linked modules at the L[1/4] balance, prove
   Σ_i log t_i − log(N𝔠 N𝔠′) ≥ c_k x/T − C_k D m. Here t_i is the least threshold at which ternary vectors with
   |L_i(a)| ≤ t_i span k^m.
   - This would confine every capped candidate to D specified proper subspaces.
   - It still leaves bounding the fresh parity rank those subspaces contribute: uncontrolled recursion gives D^{O(m)}
     families.
   - Gemini's alternative is a rank-3 "two-hyperplane" relative Liouville lemma.
   - The 72-bit witness above is the test case.
3. **R, rank-to-fresh-residual lemma for shared descent (both).** A row introducing an unseen odd large-prime column
   gives a pivot. Bound the remaining rank from rows supported on already-exposed primes or carrying large even
   powers, charging their actual stripped divisors. This is exactly what DAG-POOL needs to remove Assumption T.

## Addendum: Opus critic (read-only) on Rel-BT and the determinant obstruction

**Rel-BT: ACCEPT WITH FIXES.**
- (i), (ii), the corrected (iii) and (iv) were re-derived and are correct. The per-cone volume product is exactly
  N^{−m(m−1)/(q(q−1))}.
- **MAJOR, scope.** "Closes door (i)" holds only under all of:
  - Δ_rel ≠ 0 (relatively frame-less modules stay open);
  - the R3 model ledger (smoothness heuristic, H7 prefactor);
  - the box condition N_units·rD(2E+1)^{r−1} ≤ B/4;
  - a height budget S = poly(x).
- **MAJOR, proof gap.** The complex-place volume lemma (codex's spherical-Hölder sketch) must be written out. C^m as
  R^{2m} is not in Bean–Thunder general position. Only exp(O(m log m)) is needed.
- **MINOR.**
  - (i) is at 𝔭_k = ker(ψ mod p) only.
  - The cones remove O_k^× duplicates, not O_K^×·O_K′^×.
- **Novelty: modest.** (i) and (ii) transfer R1/R4. The new content is the finite-volume unit-cone count over fixed
  k. It is vacuous at E = 1.

**Determinant obstruction: ACCEPT WITH FIXES.**
- All steps verified: the identity, F₀ ≥ 1, the column-replacement/Hadamard contradiction, unit invariance,
  |S₊| ≥ G/(2τ), and the dyadic version with q ≤ D/(2r).
- At the exact balance the count is Θ(D) with constant ≈ 1/(4e²).
- **MAJOR.** It is a structural constraint, not a bound on N_ind or yield. Accurate wording: "few-embedding steering is
  ruled out outside a priori, unit-invariant spaces". The ≤ D subspaces (V_i = V_ī) can have dimension m−1, and their
  parity rank is uncharged.
- **MAJOR.** It holds at bounded projective height only.
- **MINOR.**
  - State KH's hypotheses and log Y ≤ (1−c)H̃.
  - Replace 2n by n₁+n₂.
  - Record the J = N𝔠·N𝔠′ typo fix at subthird_results.md:135.
- **Novelty: incremental.** A Liouville/Schmidt gap principle that quantifies KH's envelope corollary.
