# inversesieve: can DCL non-coverage be decoded into cheap containers? (r20 idea 2)

Agent `inversesieve`, phase 13, 2026-10-20.
- Pre-registration: `notes/inversesieve_prereg.md`, commit ee39c2a (before any registered run).
- Code: `src/inversesieve/{design,usets,stats,theory,experiments,analyze}.py`; exploratory scripts
  `src/inversesieve/explore_*.py`. Tests: `tests/test_inversesieve.py` (16 pass).
- Raw data: `experiments/raw/inversesieve_e{1,2,3,4}.jsonl`, `inversesieve_x_replicate.json`.
  Summary: `experiments/analysis/inversesieve_summary.json`.
- GPT-6.1-sol: blind round `codex/inversesieve_1.out`; adversarial referee `codex/inversesieve_2.out`.
  Two Opus critics (proofs, literature): §9. All their fixes are applied below.

## 0. Verdict (plain language)

**No constructive compression was found, and the route is refuted for every *forgetful* or *local*
container rule. No deterministic N^{6/31} algorithm comes out of it.**
- What survives is logically a new coverage theorem with explicit exceptions (Prop C). That is the same
  kind of missing lemma that already blocks sub-1/5 coverage. It is not an inverse-sieve statement.
- So for the r20 idea as posed, this is a kill. It is not a proof that no N-dependent container generator
  can exist; that would be a factoring lower bound.

- **The main insight: the information is diffuse.**
  - "The true pair is uncovered" has probability ≈ e^{−λ}, with λ = N^{1/31} at ϑ = 1/5.
  - That improbability is spread evenly over J ≈ N^{5/31} band primes. Each contributes only
    ρ ≈ N^{−4/31} nats.
  - Any bounded set of primes carries o(1) information: it moves the law of those coordinates by total
    variation ≤ Σρ (Prop D).
  - An exhaustive container cover of the locally admissible candidates within the N^{6/31} budget must
    combine the conditions of ≳ J/λ = J^{4/5} = N^{4/31} primes jointly. That is the coverage problem
    again.
  - With boundedly many primes, an inverse theorem exists only at Minkowski/Thue height, and there the
    design costs exactly Harvey's N^{1/5} (Cor. D2).
- **What the inverse theorem really says (Prop C).** "Uncovered ⇒ p ∈ C(N)" is the contrapositive of
  "p ∉ C(N) ⇒ covered".
  - Containers hold ≤ N^{6/31}·N^{1/4} candidates per input.
  - For a forgetful rule (Theorem L), validity *implies* coverage at λ = J^{1/5} for every pair outside an
    explicit set of density n^{O(1)}N^{−7/124} (Cor. L3).
  - The best *proved* exceptional density for the design is N^{−1/31} = N^{−4/124}, from the large sieve.
    It is not efficiently enumerable, and it is a factor N^{3/124} too large (Prop M).
- **Forgetful containers are useless (Theorem L, proved in an aligned-cell model).**
  - Model: a rule that sees N only through N mod M and a power-of-2 size window Δ, and that must hold
    for all pairs uncovered at the primes it sees.
  - Such a rule needs mass ≳ P·min(1, Δ/(PM))/n², averaged over all n-bit N. For Δ ≥ PM that is
    ≳ N^{1/4}/n² Coppersmith containers, which is no gain over brute-force windows in the
    Coppersmith/Theorem-G charging model.
  - Within the 6/31 budget, sub-rules that observe fewer than 69/124 of the bits of N catch only a
    vanishing fraction of the pairs (Cor. L2).
  - *Interpretation, not proved:* a working rule must use N non-forgetfully, which in practice means the
    exact hyperbola pq = N.
- **The inverse-sieve literature gives nothing here (Prop I; literature critic).**
  - Every known inverse-sieve theorem needs, at the primes it uses, one of: sparsity (|S mod ℓ| ≤ αℓ,
    αℓ^k or ℓ^κ); a fixed proportion of omitted classes; or ≥ b√ℓ omitted classes at primes ≍ log T
    (Croot–Mao–Pohoata–Yip 2026, Cor. 5.4).
  - Non-coverage omits only |V| ≈ N^{1/31} ≪ √ℓ classes, on one band of primes ≍ N^{5/31}.
  - All conclusions are size bounds or non-constructive existence statements. For the singleton {p}
    they say only "x − p = 0".
- **Under the random model the exceptions do not exist (Prop E, heuristic).** At ϑ = 1/5 the expected
  number of uncovered pairs, and of uncovered candidates over a dyadic range, tends to 0. The design is
  heuristically a complete cover, and the inverse theorem's whole content would be its proof.
- **Experiments (pre-registered).** 1000 exact U(N) at n = 24–40; 1842 uncovered true pairs; all 172
  dchworst worst pairs. In every tested family, no reproducible useful concentration was found:
  - |U|/model is 0.98–1.01;
  - interval cover numbers at length N^{1/4} are within 0.7% of random;
  - Coppersmith-AP triples are 0.985–1.001 of random;
  - additive energy is 0.99–1.00 of random;
  - no low-degree-image signature mod small primes;
  - Lehman containers: RR 1.02, nominal CI [0.87, 1.20];
  - p's position inside U(N) is generic;
  - worst-pair L* percentile has mean 0.507.

  Four configuration-level frozen predictions failed (P4 ×3, P4b ×1), all in the *dispersed* direction.
  One E4 endpoint was not evaluable. A cf-jump signal in E3 did not replicate on fresh seeds.
- **The cleanest illustration.** In all 88 large-J worst pairs, U(N) = {p} exactly. The uncovered set is
  then a perfect container in the information sense, yet it gives no efficient locator. Prop D and
  Theorem L rule out local and forgetful locators. Nothing better is known.
- **Side finding (experimental; mechanism UNVERIFIED).** At bands with Y ≈ P/16, true pairs are
  uncovered only 0.76–0.87 times as often as the independent-residue model predicts. The cause is
  negative dependence between band halves, i.e. coverage better than random, as in dchworst P3. It does
  not help containers.
- **Exact structure (Theorem S, proved).** For the full signed box, R_ℓ(N) = √(N·F_V(ℓ)), with
  |R_ℓ(N)| = |V| for every N and E×(R_ℓ(N)) = E×(F_V(ℓ)) exactly.
  - F_V is additively structured: a union of H dilated intervals, with 1.4–4.5× random symmetric
    energy.
  - Measured E(R) is 0.99× random symmetric energy. The square-root map hides the structure
    additively (measured, not proved).

## 1. Setting

- Notation is that of dchworst §1 and coverage38 §1. N = pq, p < q ≤ 2p, P = √N.
- The band 𝒫 is J primes ℓ ≡ 3 (mod 4) in (Y, 2Y]. Weights V are primitive (A, B) of height ≤ H.
  λ = |V|Σ1/ℓ, and ρ_ℓ = |R_ℓ(N)|/ℓ.
- F_V(ℓ) = {AB^{−1} mod ℓ}. R_ℓ(N) = {x : Bx² ≡ AN (mod ℓ), some (A,B) ∈ V}.
- U(N) = {x ∈ I : x mod ℓ ∉ R_ℓ(N) ∀ℓ}, with I = [⌈√(N/2)⌉, ⌊√N⌋].
- Since Bp² − AN = p(Bp − Aq), p ∈ U(N) ⇔ the pair is uncovered
  (`test_U_matches_bruteforce_and_true_factor`).
- **Ledger (exact; `theory.ledger`; both reviewers recomputed it).** At ϑ = 1/5:
  - X = N^{6/31} = P^{12/31};
  - Y ≍ J ≍ P^{10/31};
  - λ = P^{2/31} = N^{1/31};
  - ρ ≍ λ/J = N^{−4/31};
  - container mass budget XP^{1/2} = N^{55/124}, a fraction N^{−7/124} of I.
- A Coppersmith container is any set of ≤ N^{1/4} candidates on which Coppersmith runs in poly time: an
  interval, or a residue class mod M ≥ N^{1/4}. Larger known moduli are *easier*.
  - "≤ 3 band primes per modulus" therefore applies only to moduli M ≤ P, where a class still holds
    ≥ 1 candidate per container. What really matters is the mass count (Prop D).

## 2. Theorem S: exact forbidden-residue structure (proved)

Let ℓ ≡ 3 (mod 4) be prime with ℓ ∤ N and ℓ > 2H².

- **(a)** *V = the full signed primitive box* {(A,B): 1 ≤ |A| ≤ H, 1 ≤ B ≤ H, gcd = 1}.
  - Then F_V(ℓ) = ∪_{B≤H} B^{−1}·{A : 1 ≤ |A| ≤ H, gcd(A,B) = 1}, with distinct elements
    (|AB′ − A′B| ≤ 2H² < ℓ). So |F| = |V| and F = −F = F^{−1}.
  - For a general sign-symmetric V only |F| = |V| and F = −F hold. F = F^{−1} failed in 759/900 random
    cases (critic). Counterexample (codex): H = 2, ℓ = 11, V = {(±1,1), (±2,1)}.
- **(b)** *V sign-symmetric* (any).
  - x ∈ R_ℓ(N) ⇔ x ≢ 0 and x²N^{−1} ∈ F. Exactly one of ±Nf is a square, since −1 is a non-residue,
    and it has 2 roots.
  - So **|R_ℓ(N)| = |V| for every N with ℓ ∤ N**, and the independent-model density ∏(1 − |V|/ℓ) is
    N-independent. (`test_theorem_S_exact_size_signed`.)
- **(c)** Coverage at ℓ ⇔ t_ℓ = pq^{−1} mod ℓ ∈ F (ℓ ∤ pq). For the full box this is equivalent to the
  lattice {(A,B) : ℓ | Bp − Aq} having a nonzero vector of height ≤ H.
- **(d)** *V sign-symmetric:* E×(R_ℓ(N)) = E×(F).
  - Let F′ = {f : Nf a square}. Then F = F′ ⊔ (−F′) and R = {±√(Nf) : f ∈ F′}.
  - x₁x₂ = x₃x₄ in R squares to f₁f₂ = f₃f₄ in F′. The case f₁f₂ = −f₃f₄ cannot occur, because N²f₃f₄
    is a square and −1 is not.
  - Each quadruple lifts to 8 of the 16 sign patterns, and the same 8-to-1 count holds for F = ±F′.
  - (`test_theorem_S_mult_energy_identity`; 250/250 E2 cases.)
- **(e) Measured.**
  - E(F)/E(random symmetric set): median 1.52, range 1.08–4.65, rising with ℓ/H⁴ (40S: 4.46).
  - E(R)/E(random symmetric set): median 0.989, range 0.87–1.05.
  - Heuristic: E(F) ≍ H⁴ log H when 4H⁴ < ℓ, since sums of Farey fractions then hold in Q.
- **(f)** Proved: multiplicative structure transports from F to R. Measured only: the additive structure
  of F does not survive into R.

## 3. Logical content and bookkeeping

**Prop C (proved, contrapositive).** For any N-computable C(N), "uncovered ⇒ p ∈ C(N)" ⇔ "p ∉ C(N) ⇒
covered".
- This is a logical equivalence, not an algorithmic reduction.
- A container set of small density inside I does not by itself give a small exceptional fraction of
  inputs. Efficient construction and enumeration remain separate obligations (codex).

**Prop E (heuristic; independent model; design fixed over a dyadic range).**
- E#{uncovered pairs in S_n} ≍ (P²/n²)e^{−Λ}, and E Σ_{N∈[T,2T)}|U(N)| ≍ T^{3/2}e^{−Λ}.
- Here Λ = Σρ_ℓ = N^{1/31+o(1)}, so both tend to 0. Small-n experiments with nonempty U(N) can test only
  conditional structure.

**Prop M (proved arithmetic).**
- Containers hold ≤ N^{55/124+o(1)} candidates per N.
- The large sieve (D′ Lemma 9) gives |U(N)| ≤ (P + 4Y²)/λ, i.e. density N^{−1/31}. The set has an
  explicit membership test but is not efficiently enumerable within the budget.
- A Poisson-quality k-th moment gives density λ^{−k/2} = N^{−k/62}. Reaching N^{−7/124} needs k ≥ 7/2,
  i.e. a 4th moment.
- A theorem covering *all* of U(N) would imply |U(N)| ≤ N^{55/124+o(1)}. That beats the large sieve's
  N^{29/62} by N^{3/124−o(1)} (codex; critic).

## 4. Theorem L: forgetful container rules (proved, in an aligned-cell model)

**Model (hypotheses, not reductions; corrected after both reviews).**
- A rule is a fixed family of sub-rules i, each with an observation (M_i, Δ_i):
  - every prime factor of M_i is < P/2;
  - **Δ_i is a power of 2 with 2M_i ≤ Δ_i ≤ 2^{n−1}**, so cells [kΔ_i, (k+1)Δ_i) are aligned with
    [2^{n−1}, 2^n).
- Sub-rule i outputs C_i(D) ⊂ [P/2, P] as a function of the cell D = (N mod M_i, ⌊N/Δ_i⌋) only.
- Mass m(N) = Σ_i |C_i(D_i(N))| counts outputs with multiplicity.
- Adaptive selection among sub-rules using other data of N is outside the model.
- *Why alignment is a hypothesis.* The draft's "WLOG refine Δ to a power of 2" was **wrong**: refined
  cells need not nest. Unaligned boundary slivers break the bound. Counterexamples: critic, n = 20,
  Δ = 1028; codex, n = 20, Δ = 524280.
  - For general Δ one must add ≤ 2 boundary cells per sub-rule, i.e. a term ≤ Σ_i 2·min(Δ_i, PM_i) to
    the caught-pair counts below.

**Theorem L.** Let 𝒬_i ⊂ S_n be the pairs caught by sub-rule i, i.e. those with p ∈ C_i(D_i(pq)), and
let 𝒬 = ∪𝒬_i. Then

  Σ_{n-bit N} m(N) ≥ Σ_i |𝒬_i|·min(P/8, Δ_i/(4M_i)) ≥ |𝒬|·min_i min(P/8, Δ_i/(4M_i)).

*Proof.*
- **(i)** Restrict C_i(D) to primes, which preserves the caught pairs and lowers the mass. Then
  x ∈ C_i(D) is coprime to M_i. The q with xq ≡ d (mod M_i) and xq ∈ [kΔ_i, (k+1)Δ_i) lie in one class
  mod M_i in an interval of length ≤ 2Δ_i/P. So there are ≤ 2Δ_i/(PM_i) + 1 of them, and
  |𝒬_i| ≤ Σ_D |C_i(D)|(2Δ_i/(PM_i) + 1). (`test_theorem_L_partner_bound`.)
- **(ii)** Every aligned cell contains ≥ ⌊Δ_i/M_i⌋ ≥ Δ_i/(2M_i) n-bit integers N. So
  Σ_N m_i(N) ≥ Σ_D |C_i(D)|Δ_i/(2M_i).
- **(iii)** (Δ/(2M))/(2Δ/(PM) + 1) ≥ min(P/8, Δ/(4M)). Sum over i. ∎

`test_theorem_L_toy_rule_mass` checks the extremal rule exactly on S_20 with aligned cells.

**Corollaries.**
- **L1 (local inverse theorems).** Take 𝒬 = all pairs uncovered at the band primes S dividing ∏M_i,
  with |S| = o(Y/|V|).
  - Theorem N (dchworst §4) applies: 16Y ≤ P^θ with θ < 1/2, invertible weights, |V||S| = o(P). It gives
    |𝒬| = (1 − o(1))|S_n|.
  - If Δ_i ≥ PM_i, the mass averaged over *all n-bit integers* is ≥ cP/n². So there are ≥ cN^{1/4}/n²
    Coppersmith containers on that average, and in the worst case over inputs.
  - In Theorem-G charging, Σ_b √V_b ≥ √(ΣV_b). Each box of volume V with nonzero values < N contributes
    ≤ 2V candidates, since (P/2)³ > N. So the charged worst-case cost is ≥ cN^{1/4}/n.
  - These are statements in the charging model, not runtime lower bounds.
- **L2 (bit budget).**
  - A sub-rule with Δ_i/M_i = P^a can catch at most |𝒬_i| ≤ 8·Σ_N m_i(N)/min(P, Δ_i/M_i) pairs.
  - Within the budget Σ_N m ≤ 2^{n−1}·XP^{1/2}, sub-rules with a < 55/62 − δ together catch only an
    n^{O(1)}P^{−δ} fraction of S_n.
  - Equivalently: almost all pairs must be caught by sub-rules observing ≥ (69/62)log₂P bits of N, i.e.
    ≥ 69/124 of the n bits.
  - Lehman/stationary-point containers observe the top ≈ n/4 bits, up to an O(1) enlargement of the
    window. They fall under L1, consistent with their catch rate ≈ (#containers)·N^{−1/4}.
- **L3 (local triggers imply explicit-exception coverage).** Let validity be required only on globally
  uncovered pairs, with max_N m(N) ≤ XP^{1/2}n^{O(1)} and Δ_i ≥ PM_i/2.
  - By (i) and (ii), the caught set 𝒞 has |𝒞| ≤ Σ_i(4/P + 2M_i/Δ_i)Σ_N m_i(N) ≤ 4P·max m =
    n^{O(1)}XP^{3/2}, a density n^{O(1)}N^{−7/124}.
  - So validity **implies** that every pair outside the explicit forgetful set 𝒞 is covered. The converse
    does not hold within the per-N mass budget (critic).
  - "Covered or Lehman-easy" dichotomies are a special case. coverage38's Theorem U, B1 and Prop L
    concern certifying coverage for the remaining pairs and are unchanged.
- **Scope.**
  - Theorem L is a volume argument of the kind used in below16 T1/T2, applied to containers.
  - It says nothing about rules with Δ/M ≤ P^{55/62} or adaptive rules.
  - That such rules must "use the hyperbola pq = N" is interpretation. The virtual members of U(N) lie on
    the modular hyperbolas xy ≡ N (mod ℓ) for every ℓ, but on no common integer hyperbola. Only the
    exact N distinguishes p from them.

## 5. Diffuseness and the Minkowski endpoint

**Prop D (exhaustive local covering; proved; quantitative form of codex's bounded-k statement).**
- For S ⊂ 𝒫, let A_S = {x ∈ I : x mod ℓ ∉ R_ℓ(N) ∀ℓ ∈ S}. By S(b) and a union bound,
  |A_S| ≥ |I|(1 − Σ_{ℓ∈S}|V|/ℓ) − |S||V|.
- Covering *all* of A_S by intervals of length w, or by classes of modulus ≥ |I|/w (which hold ≤ w + 1
  points), needs ≥ |A_S|/(w + 1) containers.
- Mass ≤ N^{55/124+o(1)} = |I|·N^{−7/124+o(1)} therefore forces Σ_{ℓ∈S}ρ_ℓ ≥ 1 − o(1); the |S||V|
  endpoint term is negligible since |S||V| ≤ J|V| = o(|I|).
- Hence |S| ≥ (1 − o(1))Y/|V| ≍ J/λ = J^{4/5} = N^{4/31+o(1)}. The union bound makes this weak: the
  random model needs Σρ ≳ log N.
- *Scope:* only exhaustive covering of locally admissible candidates, not every method of locating p.
- In the independent model, conditioning on global avoidance moves the law of the coordinates in S by
  total variation ≤ Σ_Sρ_ℓ (codex; proved there).

**Cor. D2 (where inverse theorems exist; proved for this design family).**
- With s = |S| bounded, Prop D needs Σ_Sρ ≥ 1 − o(1), so H² ≳ Y/s on average over S: Thue/Minkowski
  height.
- **Example.** Take signed weights with (1,1) excluded and H = ⌊√ℓ⌋. By Thue's lemma every t ≠ 0, 1 is
  some A/B; checked for all ℓ < 800 (critic) and ℓ ≤ 400 (codex). So "uncovered at ℓ ⇒ t_ℓ = 1 ⇒
  p ≡ ±√N (mod ℓ)".
  - Two such primes give 4 classes mod ℓ₁ℓ₂ > N^{1/4}, which Coppersmith decodes.
  - The dchworst worst pairs show exactly this (E4 P9b).
- The design cost HP/Y² + sY with H ≍ √Y, plus DCL extraction and CRT decoding, is minimised at
  Y ≍ P^{2/5}, giving **N^{1/5}**: Theorem W at c = 0.
- So inverse theorems that use boundedly many primes, in this design family, land exactly on Harvey's
  exponent. This is not a classification of all bounded-prime inverse theorems.

**Prop Q (planted-survivor query model; codex; re-checked by both reviewers).**
- Hide z uniformly in I. Public data are independent of z. The algorithm may ask Q "is x = z?" queries,
  then output containers holding B integers.
- Then Pr[success] ≤ (Q + B)/|I|. With Q = N^{1/5} and B = N^{1/5+1/4}, this is ≤ N^{−1/20+o(1)}.
- This is a model statement. The arithmetic transfer is unproved.

## 6. Literature (Prop I; checked by the literature critic against the primary sources)

| result | hypothesis | met here? |
|---|---|---|
| Green–Harper, GAFA 24 (2014), arXiv 1311.6176, Thm 1.1 | A mod ℓ ⊂ S_ℓ, \|S_ℓ\| = (ℓ+1)/2, ≥ (1/16+δ)ℓ³ additive quadruples, every ℓ ≤ T^{1/2} | no informative application: band exclusions give no such S_ℓ |
| same, Thms 1.2(i), 1.3 | intervals or APs of length ≤ (1−ε)ℓ, for a proportion ε of primes or for all ℓ ≤ T^{1/2}: a fixed fraction omitted | no: only a fraction N^{−4/31} is omitted, on one band |
| Helfgott–Venkatesh, Roth volume (CUP 2009), Thm 1.1 | S ⊂ [0,T]², \|S mod ℓ\| ≤ αℓ for every prime, α fixed | no new conclusion: U × {0} lies on a known line, divisor pairs on the known hyperbola |
| Walsh, Duke 2012 (arXiv 1105.1551), Thms 1.1, 2.4 | fewer than αℓ^k classes (integer k < d); restricted version needs Σ_{ℓ∈𝒫} log ℓ/ℓ ≥ κ log Q | no; one band has prime mass O(1) |
| Walsh, GAFA 24 (2014) (arXiv 1307.0259), Thms 1.2–1.3 | \|S mod ℓ\| ≪ ℓ^κ, real κ < d | no; for {p} the polynomial is x − p |
| Gallagher's larger sieve; Croot–Elsholtz, Acta Math. Hungar. 103 (2004), Thms 2–3 | Σ log ℓ/ν_ℓ > log T; Variant 1 also needs Σ log ℓ < \|A\| | negative denominator, or a trivial gain |
| Croot–Mao–Pohoata–Yip, arXiv 2607.15311 **v2** (17 Sep 2026; v1 is a different Croot–Yip paper), Prop 5.1 | nonempty E ⊂ [T], primes ≡ 3 (mod 4), ∏ℓ ≤ ⌊√T⌋/2 (one band prime fits); ≫ (\|E\|/√T)exp(cΣ√((ℓ−s_ℓ)/(ℓs_ℓ))) points on some m ± x² | applies, but the gain certified by the complement bound is e^{o(1)}, and it is non-constructive |
| same, Cor. 5.4 | ≥ b√ℓ omitted classes at primes ≍ log T | no: \|V\| ≈ N^{1/31} ≪ √ℓ ≈ N^{5/62}, and our primes are ≍ N^{5/31} |

- The nearest literature on survivors of a few-classes-per-prime sieve does not contradict the note:
  Ford–Konyagin–Maynard–Pomerance–Tao, *Long gaps in sieved sets* (JEMS 2021); Balister et al. on the
  Erdős covering problem (Invent. 2022); Jacobsthal-type bounds; the parity phenomenon. These concern
  gaps, density and non-emptiness, not structure. Like the primes, survivors of such sieves carry no
  forced structure.
- Prior art on "inverse theorems for non-covering" in Lehman, Harvey–Hittmeir or McKee–Pinch factoring:
  none found. The searches were inconclusive and did not cover McKee–Pinch specifically.

## 7. Experiments vs frozen predictions (`experiments/analysis/inversesieve_summary.json`)

The CIs and KS p-values below are nominal. In E3 and E4 they treat observations as independent although
some share p (codex).

**E1** (100 N per config; exact U(N)).

| cfg | \|U\| mean | P1 | P2 (z) | P3 | P4 χ² | P5 | P5b | P6 | P8c |
|---|---|---|---|---|---|---|---|---|---|
| 24M | 193 | 0.997 | 1.000 (0.2) | 0.985 | **0.80** | 0.994 | 0.997 | 1.004 | 1.00 |
| 24S | 102 | 0.998 | 1.001 (0.3) | 0.987 | **0.87** | 0.993 | 1.000 | 0.992 | 1.03 |
| 28M | 242 | 1.006 | 1.007 (3.5) | 0.987 | **0.87** | 1.013 | 1.008 | 1.003 | 0.99 |
| 28S | 132 | 1.012 | 0.997 (−0.9) | 0.998 | 0.97 | 1.029 | 1.010 | 1.000 | 0.99 |
| 32M | 171 | 0.983 | 1.003 (0.8) | 0.995 | 0.96 | 0.989 | 0.986 | 0.999 | 0.98 |
| 32S | 210 | 1.002 | 0.999 (−0.3) | 1.001 | 0.99 | 1.009 | 1.002 | 0.997 | 0.98 |
| 36M | 115 | 1.003 | 0.999 (−0.2) | 0.990 | 0.97 | 0.983 | 1.003 | 1.000 | 1.04 |
| 36S | 218 | 0.994 | 1.001 (0.3) | 1.001 | 0.95 | 0.984 | 0.994 | 1.000 | 1.00 |
| 40M | 167 | 1.012 | 1.006 (1.5) | 0.988 | 0.99 | 1.056 | 1.011 | 1.000 | 0.98 |
| 40S | 240 | 0.999 | 0.996 (−1.1) | 0.997 | 0.97 | 1.038 | 0.999 | 1.000 | 1.00 |

Columns: P1 |U|/model; P2 interval cover ratio at L = N^{1/4} (z); P3 Coppersmith-AP triples; P4 χ²
occupancy ratio; P5 pair ratio vs local model; P5b halves; P6 integer energy; P8c Lehman.

- P1, P2, P3, P5 (gap KS p ≥ 0.077), P5b, P6 and P8c pass in all 10 configs. No U(N) was empty.
- No missing classes occurred among the *eligible* (|U| ≥ 4r ln r) occupancy tests. Ineligible ones do
  occur; codex: 24S idx 1 misses class 6 mod 23.
- The pooled survivor trajectory stays within 1.7% of the model at every checkpoint.
- **P4 fails in 24M, 24S and 28M** (χ² ratio 0.80–0.87 < 0.9). U(N) is *more* uniform mod small primes
  than random: the dispersed direction.
  - Candidate explanation, UNVERIFIED: U is I minus a union of APs of difference ℓ, each spread evenly
    mod r.
  - The effect vanishes at n ≥ 32 (0.95–0.99).
- **P4b fails in 40S.** The CRT collision ratio is 0.78: 26 observed vs 33.3 expected, ≈ −1.3σ in a
  Poisson approximation. Dispersed direction; the tolerance was too tight for such counts.
- |U(N)| is over-dispersed across N (var/mean 1.5–3.5). This is a fluctuation across inputs, not
  within-set clustering. Mechanism UNVERIFIED.

**E2.**
- Median E(R)/E_sym = 0.989; median E(F)/E_sym = 1.52 (40S: 4.46). E×(R) = E×(F) in 250/250 cases.
- P7 holds. Disclosed: the frozen P7 Boolean omitted the registered identity clause, which was checked
  separately.

**E3** (1842 uncovered true pairs, 2000 controls). **P8 passes (frozen pooled tolerances).**
- Lehman RR 1.02 [0.87, 1.20]; NN-percentile KS p = 0.13; near-count ratio 0.955; held-out ratio 0.932.
- Unthresholded reports: L* ≤ X_s holds for 100% of both groups (uninformative at these n). Dirichlet
  quality KS p = 0.68; L* KS p = 0.008; cf-jump KS p = 0.004.
- **Follow-ups (exploratory, disclosed).**
  - *Matched controls* (same p, q a uniform prime in the same range; `explore_matched.py`): cf-jump
    p = 0.006, L* p = 0.08, Lehman 13.6% vs 11.0%.
  - *Fresh-seed replication* (`explore_replicate.py`; 3627 uncovered pairs; matched controls):
    - Lehman RR 1.05 [0.92, 1.19];
    - mean Δlog(cf jump) = 0.021 ± 0.030;
    - cf-distribution KS p = 0.31 (codex);
    - dead weights 0.433 vs 0.423.
  - The replication did not save L*, so the L* signal was not independently re-tested.
  - These follow-ups supplied no reproducible useful enrichment. An inverse theorem would need
    enrichment by N^{7/124} with certainty, not ≲ 20%.
- **Held-out side finding (`explore_held_prime.py`; codex reproduced it).**
  - In S-shapes, partner survival on band half 2, given survival on half 1, is 0.90 / 0.86 / 0.81 of the
    independent prediction at n = 32 / 36 / 40. Prime partners give 0.90 / 0.83 / 0.81, and full-band
    uncovered rates are 0.87 / 0.80 / 0.76 of the model.
  - Random classes of the same size give 1.00, and so do the U(N) halves.
  - Same-weight exclusion (|δ_v| < Y² ⇒ at most one band factor) explains only ≈ 0.95. Co-coverage by
    distinct weights is 0.996 of independent (`explore_pair_explore.py`). The rest is UNVERIFIED.
  - Direction: uncovered pairs are *rarer* than Poisson, i.e. coverage better than random.
  - At the real scale Y² = P^{20/31} < P the size-based exclusion is generically absent. That does not
    exclude every occurrence of small δ_v.

**E4** (dchworst worst pairs).
- **P9 (L* endpoint) and P9b pass.** All 172 records are uncovered with p ∈ U(N).
- J ≥ 16 (88 records, 83 distinct pairs): L*-percentile mean 0.507, nominal KS p = 0.31; Lehman 9.1% vs
  control 13.1% / 9.0%.
- **The NN endpoint is not evaluable.** In all 88 J ≥ 16 records, U(N) = {p} exactly: 0 other survivors
  against 0.565 predicted. For J < 16, Σ(|U| − 1) = 40 389 against 40 366. Codex reproduced all 172.
  - The frozen code returned None for this endpoint and the frozen Boolean counted it as a pass. This is
    disclosed: P9 was evaluated on L* alone.
  - The summary's "U_over_model = 156" is selection bias: these N were chosen so that p survives.
- P9b, restricted to the *censored* signed J ≤ 2 records: p ≡ q (mod ℓ₁). This is the (1,1)-exclusion
  artefact, i.e. Cor. D2's Minkowski-endpoint inverse theorem. Some uncensored J = 2 records do not
  have it.

**Decision.**
- Every frozen prediction in the structured (H1) direction holds.
- Four configuration-level predictions failed (P4 ×3, P4b ×1), all dispersed. One endpoint (E4 NN) was
  not evaluable.
- Reading: **no demonstrated useful concentration in any tested container statistic.** This is negative
  evidence for the tested families, not a proof of pseudorandomness.

## 8. Candidate card

| field | content |
|---|---|
| mechanism | DCL design at λ = J^{1/5} + inverse theorem "uncovered ⇒ explicit containers" + Coppersmith |
| payoff if true | deterministic N^{6/31} for every balanced N |
| status | **refuted for forgetful aligned-cell rules (Theorem L) and for exhaustive local covers (Prop D), proved**; reduces to explicit-exception coverage (Prop C, Cor. L3); no reproducible structure in E1–E4 |
| sharpest obstruction | diffuseness: ρ = N^{−4/31} per prime; an exhaustive local cover needs ≳ J^{4/5} primes jointly; bounded-prime inverse theorems sit at Thue/Minkowski height, where the design costs N^{1/5} (Cor. D2) |
| what would revive it | a provable coverage theorem at λ = J^{1/5} outside an efficiently enumerable, non-forgetful exceptional set of density ≤ N^{−7/124}. That is a new coverage theorem, not an inverse-sieve theorem |

## 9. Reviews

- **GPT-6.1-sol referee** (`codex/inversesieve_2.out`). It reproduced the summary, all 16 tests, all 172
  E4 sets, the exponent ledger and the follow-ups. **Verdict: REJECT as written.** Its central objection
  was that the draft's general "kill"/structural-impossibility wording and several "proved" claims were
  stronger than the proofs.
  - Fixes applied:
    - headline rescoped (§0);
    - S(a) and S(c) restricted to the full box;
    - "≤ 3 primes" qualified (§1);
    - Theorem L's dyadic/Δ ≥ 2M "WLOG" turned into hypotheses;
    - L2 restated per sub-rule;
    - L3 "equivalent" → "implies";
    - Prop M ε → o(1);
    - Prop D w + 1 and scope;
    - Cor. D2 scope;
    - P4 mechanism UNVERIFIED;
    - P4b counts corrected (26 vs 33.3);
    - "missing classes" restricted to eligible tests;
    - trajectory 1.7%;
    - P7 Boolean disclosed;
    - nominal statistics;
    - replication wording;
    - E4 NN endpoint "not evaluable";
    - 83 distinct pairs;
    - P9b restricted to censored records;
    - the claim that a locator "needs" sieving withdrawn.
- **Opus proof critic.** It re-ran Theorem S on 900 cases, checked the ledger and Cor. D2 by script, and
  built counterexamples to unaligned Theorem L. **Verdict: ACCEPT WITH FIXES.**
  - Its one substantive error found: the alignment WLOG. Now a model hypothesis, with the boundary term
    stated.
  - Other fixes applied: full box for S(a) and S(c); the S(d) residue step; per-sub-rule form; which
    population the average runs over; L3 one-directional; Thue (not Minkowski) in D2; the test now uses
    aligned Δ and the stated bound.
- **Opus literature critic.** It verified Green–Harper, Helfgott–Venkatesh, both Walsh papers,
  Gallagher/Croot–Elsholtz and Croot–Mao–Pohoata–Yip (v2 only; Prop 5.1 needs primes ≡ 3 mod 4, which
  the band satisfies). It corrected the "few classes" overstatement and added Cor. 5.4, FKMPT, BBMST and
  the parity phenomenon. Prior art: none found. All of this is applied in §0 and §6.

## 10. Commands

```
.venv/bin/python -m pytest -q tests/test_inversesieve.py
.venv/bin/python -m src.inversesieve.experiments e1 --workers 9   # likewise e2, e3, e4
.venv/bin/python -m src.inversesieve.analyze
.venv/bin/python src/inversesieve/explore_matched.py
.venv/bin/python src/inversesieve/explore_replicate.py
.venv/bin/python src/inversesieve/explore_held_prime.py 40 S 30   # also 36 S 40, 32 S 60
.venv/bin/python src/inversesieve/explore_pair_explore.py 40 S
```
