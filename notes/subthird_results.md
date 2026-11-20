# subthird: a focused search for any mechanism below L[1/3] (2026-11-20)

The user's only goal is factoring general N in L[α], α < 1/3 (heuristic is fine). Two consultants, each with a
focus, then a cross-kill round. Common rules: `codex/subthird_common.md`.

- **C**, GPT-6.2-astra xhigh via codex (session 01a1f9cb-c431-7970-9ebc-66f701ae7a77). Focus: a constructive attack
  on the G_r / grdeg type-D crack, with exact experiments. Brief `codex/subthird_1_codex.md`; output
  `codex/subthird_1_codex.out`; proof and data in `experiments/subthird/codex/` (consultant_c.txt, pre-registrations,
  yield.csv, round2_*).
- **G**, Gemini 4 Carbon High via agy (conversation 9681a045). Focus: adaptive constructors, H5, composition, fresh
  paradigms and a 2024–26 literature scan. Output `codex/subthird_1_gemini{,2}.out`; scripts in
  `experiments/subthird/gemini/`.
- Cross-kill brief `codex/subthird_2.md`; outputs `codex/subthird_2_{codex,gemini}.out`.

## Verdict

**No sub-L[1/3] algorithm survives.** The last crack in the barrier map is narrowed further. The c = 0 cheap-generator
mechanism for grdeg door (i), the type-D route to L[1/4], is refuted at the generator-height level (Theorem KH
below). The door itself is not closed: sparse selection inside high-height modules remains open. What
remains is a counting loophole, LP_exc^rel. Estimates for general-N sub-1/3 factoring: sparse type-D selection
0.2–0.3%, the M2 one-witness decoder 0.1–0.2%, M3 0.05–0.1%.

## New results

### Theorem KH: the linear-Kneser height bound (C; reviewed by G; Opus critic in §Critic)

- **Statement.** Take a reduced type-D linked pair: fixed k of degree e, K and K′ of relative degree n, m ≥ 3
  k-independent generator pairs with equal evaluations mod N, and s ≥ m/2 the larger projection dimension. Let
  H = Σ_embeddings log max_j |u_j|. Then
  **H ≥ (x − D log(C_k D)) / (2⌈2n/s⌉ + 1)**, where x = log N.
- **Proof sketch.**
  1. Linear Kneser (Hou–Leung–Xiang; Bachoc–Serra–Zémor) shows that products of degree T ≈ 4n/s of the normalized
     generator span reach K × 0.
  2. Both evaluations agree on monomials and (1,0) maps to 1, so the integral monomial span has index divisible by N.
  3. Hadamard plus a covolume lower bound give exp(T·H) ≳ N.

  No unit-volume input is needed, which is the obstruction that blocked grdeg's "relative R1".
- **Consequence.** At the L[1/4] balance the main term is Ω(mx/n) = x^{3/4+o(1)}, against the x^{1/2+o(1)} the c = 0
  ledger needs, so the bound restores c = 1 (α ≥ 1/3) for generator heights. Non-power generators are included.
  With k = Q it also covers non-subfield frame-less modules (grdeg door (ii)) at the generator-height level.
- **Gaps and rescues.**
  - G found two gaps and gave a rescue for each: (i) unequal reduced degrees n₁ ≠ n₂ (the order of the two
    projections is swapped, giving T ≤ 2⌈2n/s⌉ + 1); (ii) denominator clearing in the subfield reduction (normalize
    by a unit-evaluation element w of the original span and multiply back by w^T).
  - C's caveat: the coincident-side and repeated-subfield reductions need the CC/SUB hypotheses. They are not
    unconditionally WLOG when it comes to extracting a factor.
- **Not covered.** The height-to-yield step is heuristic (typical norm sizes). Exceptional small norms are not
  bounded: C gives an explicit candidate over Q(∛3) with norms 1 and 8, independent of N, inside a high-height
  module.
- **Experiments (C, exact).** 128 boxes on 24 moduli of 60–96 bits, 437,920 coefficient evaluations. Every non-power
  construction (divided differences, full linked-lattice LLL, common sparse polynomials) had **zero incremental
  ideal-parity rank** under the aggregate cap at 84 and 96 bits. Matched power controls had 84. These are finite
  yields, not asymptotic exponents.

### Lemma FT: Farey-trapezoid count (G; C gave a shorter proof)

- **Statement.** #{primitive (a,b) ∈ [−E,E]² : 0 < |F(a,b)| ≤ X} ≤ O(D·A_F·X^{2/D}) + O(D log E), for squarefree F
  of degree D ≥ 3.
- **Consequence.** It removes the H5 box hypothesis at rank 2 up to doubly subexponential boxes, E ≤ 2^{B/(16D)}.
- **Novelty.** No exact prior statement was found. It is close to folklore (Thue-inequality counts: Bombieri–Schmidt,
  Thunder, Stewart) but has no H(F) factor and handles reducible squarefree NFS forms.
- **Check.** Pre-registered, 667/667 thresholds pass.

### Kills and corrections

- G's K2b factor-count bound k ≤ e holds. Its Newton/PTE necessary condition is **false**: C's counterexample is
  {0,4,7,11} vs {1,2,9,10}.
- G's K2c curve pullback needs squarefree pullbacks, which the repo already states (resbarrier Cor C2(d)).
- G's K2d: relative-resultant divisibility holds without a module action, but the pointwise floor fails (C's Q(∛3)
  example).
- G's adaptive-constructor kills (§1A–1D) stand. Single-interval constructors are capped by Prop SI(iii).
  Residue-steered beam search is killed by Prop RP / Lemma RT, measured at 0.60–0.77× random. The Bohr-grid moment
  method blows up to N^{Θ(1)}. DAG descent without Assumption T forces 1+γ ≤ 4α. A linear 2-D Bohr-CRT constructor
  would solve F_p-DL in L_p[β].
- G's survivors M2 and M3 are restatements of open targets.
  - **M2** (class-group positive-definite 2-D CRT capacity constructor) needs a beyond-Johnson one-witness CRT list
    decoder on quadratic norms. That is the pilist/crtdecode open primitive. C's pre-registered M2 slope decoder
    found 25 vs 2 smooth witnesses at 80 bits but lost on work: 554 lattice reductions per hit vs 474 direct norm
    tests. In effect it is a lattice sieve.
  - **M3** (non-exact multi-presentation split pencils) is nonexact O1/O2 = Π_β. C adds a counting kill when
    m = B^{ω(1)}.

## The surviving open statement (LP_exc^rel; G, sharpened by C)

- **Setting.** A reduced linked O_k-module M (relative rank m, degree D = 2en) with candidate set Ω ⊂ Z^r,
  |Ω| = e^C.
- **Counting function.** N_ind(M, Ω; Y) is the number of O_k^×-inequivalent candidates with 0 < |F| ≤ Y whose
  two-sided parity rows are independent modulo (a) candidates in proper k-subspaces and (b) products of
  lower-norm candidates.
- **What sub-L[1/3] needs.** An L[α]-constructible family (α < 1/3, C = O(F_α)) and log Y = o(x^{2/3}) with
  N_ind·ρ(log Y/log B) ≥ π(B), where B = L[α], despite H ≥ Ω(mx/n).
- **Where survivors must live.** By the envelope corollary of Theorem KH, such candidates must be scattered across
  many proper subspaces or many incomparable archimedean profiles. A single volume excess, or one exceptional small
  norm, does not suffice.

## Literature scan (G; spot-checked by both)

- **Scope.** Eleven 2024–26 arXiv items: factorial SLPs, capacity Reed–Solomon list decoding, smooth twins from
  prime lattices, Schnorr-CVP physics heuristics, double index calculus, streaming Regev, small-characteristic DL
  refinements, TNFS. All IDs resolve.
- **Corrections from the spot-check.**
  - Huang v1 vs v4 differ (v4 is a privacy paper).
  - The refined-NFS preprint dates to 2020, not 2024–25.
  - Robinson's orthotope work is a 2022 preprint (the published TNFS lattice-sieving paper is 3D).
  - Lido's quasi-polynomial small-characteristic DL paper dates to 2022.
  - The VQE paper factored 1961.
- **Finding.** None transfers to sub-L[1/3] factoring. Each lands on a repo barrier: slpcirc, crtdecode/pilist,
  siconstruct Prop CV, schnorrcvp, smoothbuild, heurnfs, or periodfind.

## Best next experiments (if continuing)

1. **LP_exc^rel directly (G).** On 64–112-bit split N at (n,m) ∈ {(8,3), (10,4), (12,4)} over Q(√2), run lattice
   enumeration in distinct conjugate-weighting cones with E > 1, targeting multi-profile cancellation
   |F| ≤ 3^D·N^{1/d}. Measure whether the proper-submodule-independent rank stays at 0.
2. **M2 vs a lattice sieve (C).** A held-out comparison forcing the same divisor size, charging total work per fresh
   relation.

## Critic (Opus, read-only): ACCEPT WITH FIXES, no FATAL

1. **Verified.**
   - The Kneser step is correct, with one step of slack.
   - The algebra dichotomy is correct: either k[V] = K×K′, or k[V] is the graph of a k-isomorphism (type A,
     k-coincidence).
   - Γ ⊆ ker is correct.
   - The covolume/Hadamard estimate is correct. Common k-units leave H unchanged, so units give no loophole.
   - At L[1/4] the main term is ≥ mx/(8n) ≈ x^{3/4}(log x)^{5/4}/8 against an error O_k(m log x), a ratio of
     x^{−1/2+o(1)}. So c = 1 at the exponent level.
2. **MAJOR, fixed: unequal degrees.** "V^{t+1} ∋ (z,0)" needs n₁ = n₂ as written. Fix: state the theorem for
   (K₀, K′₀) of unequal degrees with G's Rescue 1. This gives T ≤ 2⌈2n/s⌉ + 1 with n = max(n₁, n₂).
3. **MAJOR, fixed: (1,0) ∈ Λ.** This needs the generators to lie in K₀. G's Rescue 2 is valid, but it bounds the
   height of the original generators, which is inflated by fixed factors such as log N(𝔠) of a common ideal or the
   type-B prefactor. Fix, the **ideal-relative theorem**:
   - let 𝔠, 𝔠′ be the generator ideals, coprime to N (otherwise N splits or the pair is useless);
   - set Λ′ = S ∩ (𝔠^T × 𝔠′^T), which maps onto Z/N via (y,0) with y ≡ 1 mod N;
   - set H̃ = Σ log M_i − log N(𝔠) ≥ 0;
   - the same proof then gives **T·H̃ ≥ x − D log(C_k D)**.

   For type B, take 𝔠 = c𝔟 with a common integer c. This bounds h_v − log N(𝔟), exactly the quantity the ledger
   charges. The repeated-subfield case should be handled by applying this theorem to K₀, K′₀ directly, not by
   literal Lemma SUB substitution.
4. **MAJOR caveat, already flagged by C.** "α ≥ 1/3" is model-level only. It relies on the typical-size heuristic
   (box values ≈ H), and units can defeat that by archimedean cancellation; that gap is LP_exc^rel. Wording adopted
   here: the "c = 0 cheap-generator mechanism is refuted", not "door (i) closed".
5. **MINOR.**
   - The coincident-side case needs no reduction, since the graph case is type A. State that CC's exceptional cases
     either split N or fall into the gcd = N case.
   - The sharper numerator x − en log(2n) should be written out.
   - The envelope corollary needs non-coincident sub-spans.
   - The Bachoc–Serra–Zémor citation was not checked; the Hou–Leung–Xiang statement as used is right.

**Status of Theorem KH:** proved (internal; consultant C + G review + Opus critic), in its ideal-relative,
unequal-degree form above.
