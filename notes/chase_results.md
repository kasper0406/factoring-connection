# chase: the three open lemmas from the hatch round (2026-11-20)

**Process.** Brief rules `codex/chase_common.md`; runs `codex/chase_{L1_gemini,L2_codex,L2_gemini,L3_codex}.{md,out}`;
cross-review `codex/chase_x.md` → `codex/chase_x_{codex,gemini}.out`. Artifacts are in `experiments/chase/`.
- Codex sessions: L2 01a1fa9e-33b0-74a3-b094-028f1b0f0e6d; L3 01a1fa9e-33b1-7051-bc37-1bfeacdb4ea6.
- agy conversations: L2 e4e06b72; L1 7acc38cd.

## Verdict

**No lemma flipped into a construction. No sub-L[1/3] lead.** Each lemma is now partly proved, and what remains is
narrower and stated precisely.

## Lemma 1 (H: branch selection) — partial barrier (Gemini; codex-reviewed)

- **Theorem 1(a–c) (correct, scoped).** For d > 2K², the following reduce to each other in O(π(B) poly(x)):
  - choosing a split-prime CRT branch whose rational reconstruction gives a smooth M2 witness;
  - finding an M2 witness with Q_good > 2K².

  It also gives the exact fiber identity |B_win(Q)| = Σ_{witnesses} τ(Q_good). Random-branch bounds hold pointwise for
  every form. Upper-band branches match box sampling up to subleading order; unrestricted branches lose K^{2β+o(1)}.
  **Neither prices target-dependent selection.**
- **Theorem 2 (correct, averaged).** The two-choice weighted Theorem J for M2: every fixed depth-weight profile,
  including a profile-portfolio fixed independently of the form, gets at most exp[−((1−β)/(2κ) − o(1))F_{1−β}] smooth
  witnesses per call, averaged over balanced forms.
- **Refuted auxiliary claim.** "Tree-blindness" is false. Codex's counterexample: Q = 13u² + 4uv + 45v², N = 581, K = 4,
  B = 5.
- **Measured.** Across 30 unplanted 48–96-bit fixtures, no proxy, beam, meet-in-the-middle, tree-sieve or MCMC selector
  beat windowed random branches at fixed work.
- **Remaining.** Unrestricted computation on the public root table, where a lower bound is P vs NP-level.

## Lemma 2 (profile / LP_exc^rel) — partially proved (codex prover, Gemini adversary; cross-checked)

- **72-bit witness explained.**
  - Codex: its ternary unit representative (√2, −1, 1+√2) is covered by the ternary-box cover. The escape was an
    artifact of the wider archived box [−2,2]⁶.
  - Gemini: its measured pair-prefix (two-hyperplane) cover also captures it.
  - Codex's check: the *analytic* tube bound does NOT contain it (smallest ratio 36.2).
- **Sweep (codex, pre-registered).** 42 boxes on 15 semiprimes. Escaping unit classes fell from 140 (60-bit) to 11
  (72-bit) to 0 at 84, 96 and 112 bits, in the archived, ternary and growing-(n,m) arms alike. At 84+ bits every capped
  point lay in the span of the first two generators.
- **Threshold inequality (codex; proved, verified by Gemini).** For every h,

  Q ≥ (x − D log(C_kD))/g_h − (e/2)Dm log m − O_k(Dm),  with g_h = h(2⌈2n/s_h⌉+1).

  This uses KH applied to V^h plus a Cramer/product-formula loss of (e/2)m log m. The loss is sharp for unrestricted
  forms (Alon–Vũ). Consequences:
  - g_h = o(T), e.g. s_2/m → ∞, gives Q ≥ x/T, so fast-growing product spaces are confined.
  - Any failure needs both slow product growth, s_h = O(hm) for h ≤ cT, and loss ≥ c′m log m at a positive fraction of
    the embeddings.
- **Fresh-rank bound without recursive blow-up (codex; proved under hypotheses).** Call an embedding deep when
  |L_i(a)| < C_k m M_i e^{−m}. Capped candidates with o(m/log x) deep embeddings contribute B^{o(1)} total fresh parity
  rank, across noncoincident spans. Large coincident subfield graphs need a separate bound.
- **Rank-3 two-hyperplane relative Liouville lemma (Gemini; correct with corrected constant).**
  - Statement: for v = L_i × L_j ≠ 0, the tube {|L_i| < s, |L_j| < 2τ/s} spans k-dimension ≤ 2, where
    τ = ‖v‖/(2·3^{(3e+1)/2}(C_kE)^{3e−2}).
  - Scales: O_k(1 + log E) dyadic scales per pair suffice, by Cramer. This is a two-form extension of the determinant
    obstruction, for fixed m only: (m−1)-hyperplane tuples would cost binom(2n, m−1) ≫ 3^r.
- **Adversary constructions (Gemini).** None gives α < 1/3: all four pay c ≥ 1 on full-support candidates.
  - C1: Pisot relative units.
  - C2: multi-cone Lagrange seeds.
  - C3: common-denominator diffuse cancellation.
  - C4: Kummer and CM towers.
  - Corrections from codex: C2's "2^m − 2 small orbits" is wrong in the stated ternary Lagrange box (at most 2m−1
    qualify), though its formal products' parity rank ≤ m stands. C1 and C3/C4 general norm floors are unsupported
    beyond the specific fixtures, and the C4b 74-bit fixture fails the linking hypothesis.
- **Remaining, and both cross-reviewers' best next attack.** Divided-difference modules with total pole degree
  p = O(m). They have exactly the slow product growth (s_h ≤ h(p+1)+1). The question is whether their shared
  denominator permits near-maximal conditioning loss log(M_i/t_i) ≍ m log m at a positive fraction of embeddings.

## Lemma 3 (R: shared descent) — refuted as stated, model repair identified (codex; Gemini-reviewed)

- **Rank formula (proved).** Extracted base rank = rank M − rank M_L (ranks over Q for Z). A row introducing an unseen
  odd large column adds zero base rank at that step.
- **Counterexample (proved, priced).** A coset-table DAG: generators q³b_i mod N; table lookups certify equalities
  between already-factored words. It gains integer rank m (binary ≥ m − ω(N)) with **zero fresh smoothness
  revelations**; every stripped cofactor is 1. Cost N^{1+o(1)}, so it is not an algorithm, but it shows that
  "rank ⇒ fresh probes" is false without accounting for modular-equality certificates. Toy: N = 2773 gives integer rank
  10, binary rank 8, gcd 47.
- **Model repair (natural).** Use r(V) ≤ r(U) + #fresh, where U is generated by other certificates, including modular
  equalities, and charge every certificate's first discovery. A Dickman-type lower bound then additionally needs control
  of r(U) within the budget. That is the precise remaining term.

## Where the hatches stand

- **H.** Equivalent to the original witness search. Fixed or portfolio weighting is proved useless. What remains is
  unrestricted computation (P vs NP-level).
- **Profile.** Proved for fast product growth and for few deep embeddings. Rank 3 is handled. Open: slow-growth,
  poorly conditioned modules; the divided-difference family is the named test.
- **R.** Needs the certificate-accounting amendment. Open: bounding the rank r(U) obtainable from cached factors,
  recognized squares and modular equalities within L[1/3, c] work.
