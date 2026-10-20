# ptecons: constructing dense N-dependent split families (the systeq §6 escape)

Agent `ptecons`, phase 16 (route 1, construction side), 2026-10-20. This is the version revised after codex ptecons_2 and two Opus critics (§8).

**Files.**
- Code: `src/ptecons/`.
  - `polys.py`: exact polynomial helpers.
  - `lemmas.py`: diagnostic checkers for every lemma below.
  - `acs.py`: OPERATIONAL aligned-cluster generators. They see only N.
  - `run_experiment.py`: the harness.
- Tests: `tests/test_ptecons.py` (12).
- Pre-registration: `notes/ptecons_prereg.md` (commit 95a669f, made before E1/E2; the pilot is disclosed there).
- Config: `experiments/configs/ptecons.json`.
- Raw results: `experiments/raw/ptecons_E{1,2}.jsonl`.
- Codex (GPT-6.1-sol): blind round `codex/ptecons_1.{md,out}`; adversarial round `codex/ptecons_2.{md,out}`.

**Labels.**
- **PROVED**: proof given here or in the cited file.
- **KNOWN**: published, with a citation.
- **HEURISTIC**: the assumption is named.
- **measured**.
- **refuted for this model**.
- **open**.
- **assessment**.

Novelty is never claimed; "no prior statement found by us" is the most we say.

---

## 0. Verdict

**Plain language.** We found no construction meeting systeq §6 (a)–(c).

We examined seven families: Dickson/Chebyshev, Prouhet–Tarry–Escott (PTE), cyclotomic/Galois orbits, elliptic, group-ring, lattice/LLL, and our own aligned clusters. Every one fails at an identifiable step. The list of failure modes (§4) is a summary of these families, not a classification of all possible constructions. The entries that stay **open** are listed explicitly in §7.

**What the finite-field engine has that Z lacks.**

1. **Linear families.**
   - *In F_q.* The roots b + c·a (c ∈ F_q) move affinely in (a, b), while every middle coefficient of the split product stays zero (Σ_{c∈F_q} c^j = 0). That makes the alignment to the presentation a *linear* condition for a whole family at once.
   - *Over R.* If the roots move affinely in real parameters and the X^{M−2} and X^{M−3} coefficients stay affine, then at most one root moves (**LIN-R**, PROVED). It is a corollary of de Gua's rule, which is KNOWN.
   - *Complex directions.* Vanishing power sums are possible even over Z[i]; codex gives an explicit example. But for "star" families over an imaginary quadratic field, Hajós's lemma plus archimedean non-units force the λ-dependent part of the residual to have degree ≳ M′/12 or height ≳ 2^{M′/4} (**STAR**, PROVED counting + HEURISTIC NC).
2. **Constants are units in F_q but not in Z.**
   - In F_q[X] the "norm" of X − a is a nonzero constant, which costs nothing.
   - In Z the norm ledger ∏ F(a_l) = ±Res(F, P mod F) charges every root with |F(a_l)| ≥ 2 (**LEDGER**, PROVED, elementary).
   - For the Cunningham/F_1 presentation X^n ≡ 1 and even n: ∏_{|a_l|≥2}|a_l| ≤ ‖P mod (X^n − 1)‖₁ (**CUN**, PROVED; the odd-n form was corrected after codex). So a compressed split polynomial has at most log₂‖r‖₁ roots outside {−1, 0, 1}.

**The PTE route specifically.**
- **Ideal classes (one free coefficient) give no common-offset gain.**
  - Kleiman (1975) proved (M − 1)! | C for the constant C = P − Q [KNOWN].
  - So in a cluster sharing one aligned offset, residuals are pairwise spaced by multiples of (M − 1)!. At most 1 + 2V/(M − 1)! members have |v| ≤ V [PROVED].
  - With M ≥ x/(C₀ log B), all but O(1) members have log|v| ≳ (x/(C₀ log B))·log(x/log B).
  - This exceeds x^{2/3} for B = L[α] with α ≤ 1/3; at α = 1/3 the margin is only a log factor (codex).
  - The cost reading (no gain below L[1/3]) uses the random-residual model [HEURISTIC].
- **Near-ideal classes (t ≥ 2 free coefficients).** They obey a proved coefficient-height bound (**CW**, codex ptecons_1; no prior general-t statement found). Their *evaluated* residuals are controlled only under a no-cancellation heuristic NC (§3.5).
- **Generic classes.** Mean-occupancy counting fixes only about √(2M) top moments [HEURISTIC], so the residual is ≈ N^{1−O(M^{−1/2})}.
- **Structured chains** (Prouhet–Thue–Morse, Frolov, Choudhry) fix only O(log M) moments.

**Experiments (pre-registered; two control predictions failed as stated, disclosed in §5).**
- **E1.** Five structured families fixed independently of N show the uniform congruence rate over 200 semiprimes: G = 1.05, 1.07, 1.04, 0.96, 0.99. This is a sanity check close to a tautology. Codex shows that "N-independent ⇒ no excess" is *false* for adversarial U (squares; lcm(range) + 1), so AVG is only a conditional heuristic (§3.1).
- **E2.** Three N-aligned higher-degree split families were tested against COS (degree 2) on 64- and 80-bit semiprimes:
  - No aligned variant had a median residual exponent below COS's by more than 0.01: λ = log|v|/log N of 0.592/0.574 for COS against 0.565–0.606 for the generic quartic.
  - On a post-hoc per-column yield, smooth relations per factor-base column, COS wins at every N: 0.41–0.44 and 0.16–0.17, against ≤ 0.31 and ≤ 0.15.
  - Two-digit sextic alignment is infeasible: 0–4 members against a cap of 20000.
  - These are toy-scale results, confounded by different root bounds and pool truncation (§5).

**Verdict.**
- No candidate meets (a)–(c).
- Refuted for this model (theorem plus the stated heuristics):
  - Dickson;
  - ideal and common-offset PTE;
  - generic and structured PTE chains;
  - cyclotomic/Galois orbits;
  - elliptic preimages;
  - (C₂)^r group rings;
  - real-linear families;
  - imaginary-quadratic stars (under NC);
  - aligned clusters, which are no better than L[1/2].
- **Open** (§7):
  - near-ideal classes whose difference polynomials cancel systematically at an N-aligned point;
  - non-star complex-direction families;
  - an output-sensitive enumerator of split polynomials in the N-lattice {P : P(m) ≡ small mod N}. That enumerator would be a sufficient primitive; it is smoothenum's algebraic or doubly-smooth listing.
  - No candidate is known for any of them.
- Assessment: < 1% that a split-by-identity construction yields sub-L[1/3]. Codex (blind round) estimates about 2%.

---

## 1. Formulation (aligned with systeq §6 and the ptesplit_1 model)

**Inputs.** A relation is a pair of actual integers (u, v) with u ≡ v (mod N), u ≠ v. The data are:
- an evaluation point m (with F(m) ≡ 0 mod N for a presentation F, or simply a residue);
- a split polynomial P = ∏_{l=1}^M (X − a_l) with |a_l| ≤ A;
- u = P(m) = ∏(m − a_l). The pieces are small integers when |m| + A ≤ B, or factor-base symbols otherwise;
- v = (P mod F)(m), or a residue of u mod N.

**Requirements.**
- Smallness of r = P mod F is not smallness of v (ptesplit_1).
- Nontriviality forces M ≥ (x − O(1))/(C₀ log B) when the pieces are ≤ B^{C₀}.
- Cardinality is not rank. Repeated or dependent rows count once, and pieces factored into primes merge columns.

**Cost threshold** [HEURISTIC, random-residual model; codex ptecons_1 and _2, ptesplit_1].
- With log|v| = c·x^γ(log x)^{1−γ} and cheap generation, the optimum is L[γ/2]. A fixed exponent below 1/3 needs γ < 2/3. log|v| = o(x^{2/3}) alone may only improve the constant.
- If v is smooth *by identity*, the size requirement disappears (row 2b). The alignment cost then dominates.

**Cluster + alignment decomposition (our organising device).**
- A *cluster* is a class 𝒞(T) of split polynomials with common top part T: P = T + r_P, deg r_P < t. Equivalently, the members have equal e₁..e_{M−t}: a multigrade PTE class.
- An *alignment* is a condition T(m) ≡ w (mod N) with |w| small.
- One alignment yields |𝒞| relations u_P ≡ w + r_P(m) with a **common offset** w. These are genuinely distinct multiplicative relations: the aligned congruence plus different exact identities u_P − u_{P′} = (r_P − r_{P′})(m). This is the BGJT shape.
- Re-centering each u_P separately mod N forfeits the cluster gain, since each member then needs its own modular condition. (The E2 harness reports the re-centered |v| as a size statistic; λ is unaffected.)

## 2. Candidate table

| # | Family | Density | Residual | N-dependence | Verdict |
|---|---|---|---|---|---|
| 1a | Dickson D_n(X) − D_n(t) over Z | Splits only for n ∈ {1, 2} (all t), n ∈ {3, 6} (t = ±1, ±2), n = 4 (t = 0, ±2) [PROVED; codex; our check n ≤ 30] | Exact identity, e.g. (X² − 4)(X² − 1)² = D₆ − 2 | None: F = P − r is chosen from P (EX) | refuted: (c) |
| 1b | Dilated/shifted D₆: ((X−h)² − 4λ²)((X−h)² − λ²)² | Θ(A²) | One fixed degree-6 F pins h and λ² [PROVED, codex] | Per-member alignment | refuted |
| 1c | D_n − c in a fixed field K | φ(n)/2 ≤ [K:Q]; O((log A)^{r_L}) shifts (τ is a unit) [PROVED, codex; Dirichlet] | — | — | refuted: (a) |
| 2a | Ideal PTE classes (t = 1) | ≤ q = 2A + 1 per class; ideal solutions known for M ≤ 10 and M = 12 | Common offset: spacing (M−1)! (Kleiman) | Alignment ≈ N/V trials | **refuted** (PROVED spacing + HEURISTIC cost) |
| 2b | Ideal class with T(m) ≡ 0: u = ∏(m − a), v = ±∏a, both smooth (codex) | as 2a | Smooth by identity; ≥ C/2 | Needs N \| T(m): a root mod N of a fixed degree-M polynomial | refuted: alignment ≈ N [HEURISTIC] |
| 2c | Generic PTE classes | Mean A^{M−k(k+1)/2}, so k ≲ √(2M) [HEURISTIC]; VMVT for fixed (M, k) [KNOWN, BDG 2016; Wooley 2016] | ≈ N^{1−O(M^{−1/2})} | Prescribed top part | refuted: (b) [HEURISTIC] |
| 2d | Prouhet–Thue–Morse, Frolov lifts, Choudhry chains | Up to A^r parameter tuples, M = 2^r (or 2j^{k−1}) | Only O(log M) moments fixed | none | refuted: (b) |
| 2e | Near-ideal classes, 2 ≤ t ≪ M | ≤ C(q + t − 1, t) per class [PROVED, ptesplit_1/ptecons_1] | Height H ≥ K_s/(t2^{t−1}A^{t−1}) [PROVED, CW]; evaluation only under NC | Prescribed top part | refuted under NC; **open** otherwise |
| 3a | μ_n-orbits ∏_j(X − h − aζ^j) = (X − h)^n − a^n | Θ(A²) | One binomial per orbit (Galois collapse) | — | refuted |
| 3b | Projective μ_n orbits (dX − b)^n − (a − cX)^n | Θ(A⁴) | Mixed terms do not compress; F′ applies [PROVED, codex] | — | refuted |
| 3c | Cunningham N \| Φ_n(m), X^n ≡ 1 | — | CUN: ≈ log₂‖r‖ roots outside {−1, 0, 1} [PROVED] | Special N | refuted: reduces to systeq M13 |
| 3d | Jacobi/Gauss sums (Stickelberger) | O(n²A²) | Ideal relations in Z[ζ_n]; image ≈ N/m | Needs ζ_n ↦ m mod N | refuted: (b) |
| 3e | Complex directions over Z[i] (codex: (X⁴ − 1)⁴(X⁴ + 4), p₁..p₇ = 0; general k: A = O(k³), M ≤ e^{O(k² log k)}) | — | STAR (§3.2): λ-part degree ≳ M′/12 or height ≳ 2^{M′/4} | — | stars refuted under NC; non-stars **open** |
| 4 | Division polynomials, Lattès preimages | O((log A)^{r/2}) points [KNOWN, Néron–Tate]; φ(n) ≤ 2[K:Q] [PROVED, codex; Weil pairing] | Denominators ψ_n² | None | refuted: (a), (c) |
| 5a | Group ring (C₂)^r: ∏_ε(X − Σε_j a_j) | ≈ A^r/(r!)², M = 2^r | ≤ 2 log₂ M leading coefficients can agree [PROVED, codex] | — | refuted: (b) |
| 5b | Abelian group determinants; resultants ∏(a_i − b_j) | — | Cyclotomic norms; translation-invariant | A zero residual reveals a gcd | refuted |
| 6 | Lattice Λ_{N,m} ∩ split variety; Schnorr–Adleman | Split density in the coefficient box A^{−M(M−1)/2} [PROVED, codex] | Expected count ≈ S·2V/N [HEURISTIC] | Genuine | **open** as an algorithm; concrete variants fail (§6) |
| 7a | Aligned clusters ACS (§3.6) | Pool-limited | ≥ N^{1/2} typical [HEURISTIC]; E2 | Base-m digit matching | refuted: no better than L[1/2] |
| 7b | Composition clusters (X² − c)², a² + b² = 2c | r₂(2c) = A^{o(1)} per c | as 7a | none | refuted (E1 G = 0.99) |
| 7c | Real affine-root families | — | LIN-R: one moving root | — | refuted (PROVED) |
| 7d | Block products ∏_j{P_j or Q_j} | 2^r | rank ≤ 2r [PROVED, codex] | — | refuted |

## 3. Lemmas

### 3.1 AVG — conditional equidistribution [HEURISTIC; the universal form is false]

**Exact average.** For N uniform over a set 𝒟 with N > 2V,

  E_N #{u ∈ U, |u| > N : |u mod± N| ≤ V} = |𝒟|^{−1} Σ_{u∈U} Σ_{|r|≤V} #{N ∈ 𝒟 : N | u − r}.

**The uniform rate Σ_u (2V+1)/N needs a hypothesis.** The divisor counts of u − r in 𝒟 must be typical. Two counterexamples (codex ptecons_2, PROVED):
- for semiprime N = pq, t² ≡ 0, ±1 has 5 or 9 solutions mod N, against 3 for a uniform residue;
- u = lcm(𝒟) + 1 succeeds for every N ∈ 𝒟.

**The heuristic we use.** Families uncorrelated with the N-population get the uniform rate. Structure without N-dependence then gives no *systematic* excess. E1 (G = 0.96–1.07) is consistent but nearly tautological.

### 3.2 GAP-R, LIN-R, STAR — rigidity of linear families

**GAP-R** [KNOWN in substance: de Gua's rule (1741) and Popoviciu (1935) §§2, 9, as located by the critics. We re-prove it.]
- *Statement.* Let y₁..y_M be real with e_i(y) = e_{i+1}(y) = 0 for some 1 ≤ i ≤ M − 1. Then at most i − 1 of the y_l are nonzero.
- *Proof.*
  - R(Z) = ∏(Z + y_l) is real-rooted.
  - Its (M − i − 1)-th derivative has Z¹ and Z⁰ coefficients proportional to e_i and e_{i+1}, so it is Z²·(…).
  - A multiple root of f^{(r)}, for real-rooted f, is a root of f of multiplicity ≥ r + μ.
  - So 0 has multiplicity ≥ M − i + 1 in R. ∎
- *Check.* Exhaustive, M ≤ 6 in [−5, 5]: 1364 double-gap cases, 0 violations.

**LIN-R** [PROVED, an immediate corollary of GAP-R; no prior statement found].
- *Statement.* Let M ≥ 3 and P_λ(X) = ∏(X − α_l − L_l(λ)), with α_l ∈ C and L_l real linear forms on R^s. If the X^{M−2} and X^{M−3} coefficients are affine in λ (as identities), then at most one L_l ≢ 0.
- *Proof.* The degree-i homogeneous part of e_i(α + L(λ)) is e_i(L(λ)). Apply GAP-R with i = 2 pointwise, then pick λ outside the kernels of any two nonzero forms. ∎

**Complex directions** (corrected after codex).
- *Wording correction.* e₂ = … = e_k = 0 gives p_j = e₁^j, not p_j = 0, unless e₁ = 0.
- *No growing field is needed.* Vanishing p₁..p_k is possible over a fixed field. Example: (X⁴ − 1)⁴(X⁴ + 4) = X²⁰ − 10X¹² + 20X⁸ − 15X⁴ + 4 splits over Z[i], has house √2, and has p₁ = … = p₇ = 0. Codex (PROVED sketch) builds, for every k, a Gaussian multiset with p₁..p_k = 0, house O(k³) and M ≤ e^{O(k² log k)}.
- So the first draft's claim that complex directions force cyclotomic fields was **false**.
- What remains is a *compression* problem: in these examples M grows much faster than k.

**STAR** (complex directions over an imaginary quadratic K) [PROVED counting + HEURISTIC NC].

*Setup.*
- c₁..c_M ∈ O_K are conjugation-closed, M′ of them nonzero; α, λ ∈ Z with λ ≠ 0.
- P_λ(X) = ∏(X − α − λc_l) ∈ Z[X].
- Write ∏(Y − c_l) = Y^M + Σ_{j≤s} g_j Y^j with g_s ≠ 0. Then

  P_λ(X) = (X − α)^M + Σ_{j≤s} g_j λ^{M−j}(X − α)^j,

  so the λ-dependent part has degree s.

*Counting* [PROVED].
- Let g_low be the lowest nonzero g_j, which is ± the product of the nonzero c_l.
- Non-units of an imaginary quadratic O_K have |c| ≥ √2, so at most 2 log₂|g_low| of the c_l are non-units.
- Each of the w_K ≤ 6 roots of unity has multiplicity ≤ s + 1, by **Hajós's lemma** [KNOWN, Hajós 1953: a nonzero root of a T-term polynomial has multiplicity ≤ T − 1].
- So M′ ≤ 6(s + 1) + 2 log₂|g_low|.
- Checked by exact tests (`star_gaussian_check`); the codex Z[i] example satisfies it (M′ = 20, s = 12).

*Consequence under NC.*
- The λ-dependent residual is ≳ max((m − α)^s, |g_low|).
- Either s ≥ M′/12 − 1 or log₂|g_low| ≥ M′/4. So log|v| ≳ M′·min((log m)/12, (ln 2)/4).
- With M′ ≍ x/log B this is ≫ x^{2/3} for α < 1/3.
- If M′ = o(M), the family is a fixed relation times exact factors and adds no rank.

*Not covered (open).* Non-star complex families (varying α_l, mixed moments), and fields with units of infinite order.

### 3.3 LEDGER and CUN — constants are not units [PROVED, elementary]

**LEDGER.** For monic F of degree d, P = ∏(X − a_l) and r = P mod F:

  ∏_l |F(a_l)| = |Res(F, r)| ≤ ‖r‖₁^d · M(F)^{deg r}.

- This is norm multiplicativity of ∏(θ − a_l) = r(θ) in Z[θ].
- It is not a universal obstruction by itself. For example, F = P − 1 gives F(a_l) = −1 at every root (codex). The difficulty then moves to arranging N | F(m).
- Its content is the contrast with F_q[X], where the norm of X − a is a unit constant.

**CUN.** For n ≥ 3, integer roots a_l and r = P mod (X^n − 1), let z be the number of nonreal n-th roots of unity: z = n − 2 for even n and n − 1 for odd n. Then

  ∏_{|a_l|≥2}|a_l|^{n−2}(1 − 2^{−n})^M ≤ ‖r‖₁^z,

and for even n the sharper ∏_{|a_l|≥2}|a_l| ≤ ‖r‖₁ holds (codex).

- *Proof.* At each nonreal ζ, P(ζ) = r(ζ), so |P(ζ)| ≤ ‖r‖₁. Then:
  - for even n and |a| ≥ 2, ∏_{ζ nonreal}|ζ − a| = |a^n − 1|/|a² − 1| ≥ |a|^{n−2};
  - for odd n it is |a^n − 1|/|a − 1| ≥ |a|^{n−2}(1 − 2^{−n});
  - for a ∈ {−1, 0, 1} the product is ≥ 1. ∎
- *Correction.* The first draft used z = n − 2 for odd n. Codex's counterexample is n = 3 with roots {1, (−2)^84}, failing by a factor of 1.04; it is kept as a regression test. The tests use exact integer comparison.
- *Background.* Res(X^n − 1, P) = ±∏(a_l^n − 1) is the Pierce–Lehmer form [KNOWN: Pierce 1916; Lehmer, Ann. Math. 34 (1933)].
- *Consequence for systeq M13 / F_{1^n}* (even n; odd n needs M = o(n2^n) for the same statement up to 1 + o(1)).
  - A split polynomial with ‖P mod (X^n − 1)‖₁ = h has at most log₂ h roots outside {−1, 0, 1}.
  - Its u-side therefore uses essentially the three symbols m − 1, m and m + 1.

### 3.4 Ideal clusters (Kleiman) and CW (near-ideal height)

**Ideal case** [KNOWN: Kleiman, J. reine angew. Math. 278/279 (1975) 48–51; Rees–Smyth, LNM 1415 (1990) Prop. 2.1; restated as Prop. 2 in Costello–Meyer–Naehrig, EUROCRYPT 2021].
- If P ≠ Q are monic, split over Z, of degree M and with P − Q = C constant, then (M − 1)! | C.
- Short proof (codex): Σ C(a_l, M) − Σ C(b_l, M) = −C/(M − 1)! ∈ Z.
- Every prime of C is ≤ 2A.
- *Spacing consequence* [PROVED]. Under one common offset, residuals of distinct members differ by nonzero multiples of (M − 1)!. So #{members with |v| ≤ V} ≤ 1 + 2V/(M − 1)!.
- *Threshold.* With M ≥ x/(C₀ log B) and log B = x^α(log x)^{1−α}, every member but O(1) has log|v| ≥ M log M − O(M) ≈ x^{1−α}(log x)^α. This is ≫ x^{2/3} for α < 1/3, and still larger by a log factor at α = 1/3 (codex).
- The resulting cost statement, that ideal clusters give no gain below L[1/3], is HEURISTIC (random-residual model). It does not cover 2b, where v is smooth by identity; 2b is killed by its alignment cost instead.
- Our first draft's ½M log M bound (from codex's support argument) is superseded by Kleiman.
- *Measured.* All known ideal solutions we checked (sizes 2–8, 10, 12) satisfy both bounds; for example, M = 12 has |C| = 6.7·10²².

**CW, the near-ideal height bound** [PROVED; codex ptecons_1, `codex/scratch/ptecons_1/ptecons_pte_height_lemma.md`, re-checked step by step by us and by the proof critic; no prior general-t statement found].
- *Statement.* Let P ≠ Q be monic of degree M, split over Z with roots in [−A, A], with deg(P − Q) < t. Put k = M − t, s = ⌊k/2⌋ + 1 and K_s = ⌊s/2⌋!⌈s/2⌉!. Then ‖P − Q‖_∞ ≥ K_s/(t 2^{t−1} A^{t−1}).
- *Proof idea.*
  - Remove the gcd.
  - Each coprime cofactor has ≥ s distinct roots: otherwise W = ∏(X − a)² over one support contradicts the equal moments. This is the classical argument behind "a degree-k solution has size ≥ k + 1".
  - Evaluate (P − Q)^{(h)} at a root of Q's cofactor.

**Barrier side.** The companion results are ptebarrier's Lemma R (`notes/ptebarrier_prereg.md`) and ptesplit's same-remainder theorem.

### 3.5 Near-ideal clusters under NC [HEURISTIC]

**NC (no cancellation).** For a class difference D = P − Q of height H, |D(m)| ≳ max(H, m^{deg D}), where deg D is the actual degree for that pair.

**Under NC, in a common-offset cluster:**
- when (deg D)·log A ≲ (M/4) log M, CW gives log|v| ≳ (M/4) log M;
- otherwise log|v| ≳ deg D · log m.

For pieces ≤ B = L[α], α < 1/3, and m ≥ A^{Ω(1)}, both exceed x^{2/3}.

**Not settled:** NC itself; the case of very small m with t ≈ M; and systematic cancellation of D at an N-aligned m. The last is exactly the open primitive.

### 3.6 Aligned clusters (ACS, our best candidate) and alignment capacity [HEURISTIC + measured]

**Algorithm** (`src/ptecons/acs.py`).
1. Take m = ⌈(kN)^{1/M}⌉ and the offset R₀ = m^M − kN.
2. Greedily choose e₁, e₂, … of a split P with roots in [−A, A] to cancel R₀'s leading base-m digits.
3. Enumerate the class realising them. For the quartic this uses a two-sums meet-in-the-middle, checked against brute force.
4. Set u = ∏(m − a_l) (pieces ≤ m + A) and v = u − kN. The harness reports the centred residue.
5. Smoothness-test v, then do linear algebra.

**Capacity** [HEURISTIC].
- Cancelling the m^{M−2} digit needs |e₂| ≈ m. Since Σa_l² = e₁² − 2e₂, this forces A ≳ √(m/M).
- For a *typical* member under a root-distribution model, the lower coefficients e_i ≈ A^i make |v| ≳ N^{1/2}.
- This is not a deterministic bound. Codex's example (t, −t, 0, …, 0) has large e₂ and e_M = 0. Such members, however, have few distinct pieces, and their residual is ≥ m^{M−J} when only J roots are nonzero.
- Generic realizability caps the matched digits at about √(2M).
- M = 2 is the linear-sieve/COS case, |v| ≈ A√N.

**Cost.** ACS is no better than L[1/2].
- At fixed M with A ≳ m^{1/2}, the factor base would be about N^{1/(2M)}, which is exponential.
- Keeping A subexponential forces M → ∞, and then the √(2M) capacity gives N^{1−O(M^{−1/2})}.
- Adding an algebraic side gives NFS, which RES + G (heurnfs B1) bars from beating 1/3.

### 3.7 Collapse and rank [PROVED, elementary]

- **Galois collapse.** A full μ_n-orbit gives ∏_k(x − ζ^k t) = Res(Z^n − 1, x − tZ) = x^n − t^n: one relation.
- **Disjoint root supports** (Galois/Dickson fibres) give rows ≤ columns over the symbols {m − a}. When pieces are factored into primes, columns merge, so this is a no-excess statement only for symbol factor bases.
- **Block products** (codex) span rank ≤ 2r for 2^r members.
- **Comparison.** In BGJT each linear factor occurs in ≈ q² of the q³ members.

## 4. Failure modes of the families examined (a summary, not a classification)

| Mode | Statement | Status | Applies to |
|---|---|---|---|
| No N-dependence | uniform congruence rate | HEURISTIC (conditional) | 1a, 2d, 3a, 5, 7b |
| Collapse / rank | orbit products are norms; disjoint supports; block rank | PROVED | 1c, 3a, 3d, 7d |
| Fixed-field degree bounds | Dickson φ(n)/2 ≤ D; torsion φ(n) ≤ 2D; Splitting Lemma (systeq §5.4) | PROVED | 1c, 3a, 4 |
| Archimedean non-units | LEDGER / CUN / STAR | PROVED (+ NC for STAR) | 3c, 3e stars |
| Real rigidity | LIN-R (de Gua) | PROVED | 7c |
| Ideal spacing | Kleiman (M−1)! \| C | KNOWN / PROVED | 2a |
| Alignment capacity | √(2M) digits; typical v ≳ N^{1/2} | HEURISTIC, measured | 2c, 7a; 2e, 6 only partially |

Each family examined fails at least one mode. The BGJT orbit in F_q satisfies none of the hypotheses of any mode.

## 5. Experiments (pre-registered, `notes/ptecons_prereg.md`, commit 95a669f)

### E1 — uniform-rate check (`experiments/raw/ptecons_E1.jsonl`)

**Setup.** 200 balanced 24-bit semiprimes, V = 1000, |u| > N only.

| family | size | pairs | observed | expected | G |
|---|---|---|---|---|---|
| generic sextic, [−6, 6], m = 40 | 18564 | 3.71M | 615 | 584.9 | 1.05 |
| symmetric J = 3 clusters | 11575 | 2.32M | 389 | 364.7 | 1.07 |
| PTE-class members (M = 5, e₁..e₃ equal) | 1132 | 226k | 37 | 35.7 | 1.04 |
| Dickson integer-split members | 7826 | 1.48M | 224 | 233.8 | 0.96 |
| composition clusters (a² + b² = 2c) | 17512 | 2.49M | 390 | 392.9 | 0.99 |

**P1** (G ∈ [0.8, 1.25] when expected ≥ 200) **held** for all four eligible families; the PTE row is below the threshold. Per §3.1 this checks the heuristic for these families only.

### E2 — aligned clusters (`experiments/raw/ptecons_E2.jsonl`)

**Setup.** 64-bit semiprimes (5 seeds) and 80-bit (3 seeds); pool cap 20000; B = 2^16. The smoothness control is a random integer of the same bit length.

| variant | median λ (64-bit) | median λ (80-bit) | λ₁₀₀₀ (64 / 80) | members (64-bit) | smooth / control | smooth per column (post hoc; 64 / 80) |
|---|---|---|---|---|---|---|
| cos (M = 2, A = 100) | 0.592 | 0.574 | 0.537 / 0.530 | 20000 | 0.85–0.91 | 0.41–0.44 / 0.16–0.17 |
| quart_g (e₁ exact, e₂ window) | 0.595–0.606 | 0.565–0.583 | 0.556–0.584 / 0.524–0.540 | 3072–13093 | 0.76–0.98 | 0.06–0.27 / 0.11–0.15 |
| quart_s (symmetric quartic) | 0.633 | 0.606 | 0.578–0.581 / 0.562–0.564 | ≈ 19100 | 0.96–1.03 | 0.29–0.30 / 0.01 |
| sext_s (two digits matched) | — | — | — | 1, 1, 1, 1, 4 (80-bit: 0) | — | — |
| sext_s_wide (one digit) | 0.674–0.677 | 0.647–0.649 | 0.695–0.699 / — | 881–1101 | 0.79–1.36 | 0.20–0.31 / 0.01 |
| quart_rand (control) | 0.855 | 0.860 | 0.80 / 0.82 | 20000 | 0.54–0.93 | 0.04 / 0.004 |

The per-column metric is smooth relations divided by π(B) plus the column count. The column count is the number of distinct pieces, or π(max piece) if smaller. It was added after the critics' review and was not pre-registered.

**Evaluation.**
- **P2a** (cos median = 0.5 + log₂100/n ± 0.04; predicted 0.604 and 0.583): measured 0.592 and 0.574. **Held.**
- **P2b** (every aligned variant ≥ cos − 0.02 at every N): **held**. The closest case is quart_g at 80-bit seed 2, 0.565 against 0.574, which is *below* COS by 0.009.
  - So "no variant beat COS" holds only within the registered 0.02 margin.
  - On the post-hoc per-column metric COS wins at every N.
  - *Confounds* (proof critic): root bounds differ by variant; the cap keeps the first 20000 in loop order; quart_s uses 57k–74k distinct pieces for about 19k members.
- **P2c** (λ₁₀₀₀ ≥ 0.45): **held** (minimum 0.524).
- **P2e** (sext_s < 200 members): **held** (0–4).
- **Kill signal:** not triggered.
- **P2d** (smooth ratio ∈ [0.6, 1.6]): **failed** in one cell, quart_rand at 80-bit seed 2: 0.54, from 25 against 46 hits.
  - *Explanation.* With k = 1, v = ∏(m − a) − N is odd with probability 15/16 (and similarly mod 3). This is the local-root effect of systeq E4.
  - COS with two factors shows 0.85–0.91; quart_g, whose k varies, is near 1.
  - The 64-bit quart_rand cells already sit at 0.65–0.75, so 0.54 is this bias plus Poisson noise.
- **P2f** (quart_rand median λ ≥ 0.95): **failed as stated** (0.855 / 0.860).
  - The "unaligned" control is implicitly aligned on its leading digit (m = ⌈N^{1/4}⌉), so |v| ≈ |e₁|m³.
  - With median |e₁| ≈ 0.8A this predicts λ ≈ 0.854 / 0.859, matching the measurement.
  - The prediction was mis-specified.

**Reading.** These are toy-scale measurements. The asymptotic claims rest on §3.6, not on E2. No scaling point beyond 80 bits was run.

## 6. Literature (corrected after the literature critic)

- **The BGJT line.**
  - Joux, SAC 2013 (L[1/4]);
  - Barbulescu–Gaudry–Joux–Thomé, EUROCRYPT 2014 (arXiv 1306.4244);
  - Granger–Kleinjung–Zumbrägel;
  - Kleinjung–Wesolowski, J. AMS 35 (2022).
- **Closest prior remark.** Granger–Kleinjung–Zumbrägel, "Indiscreet logarithms…", AMC 12 (2018) §7: prime fields need "fundamentally new ideas". No published integer-factoring transfer was found.
- **PTE.**
  - Wright, Quart. J. Math. 6 (1935), the conjecture.
  - Borwein–Ingalls, Enseign. Math. 40 (1994).
  - Borwein, *Computational Excursions…* (2002), Ch. 11.
  - Borwein–Lisoněk–Percival, Math. Comp. 72 (2003): symmetric searches; new size-10 solutions.
  - Size 9: Letac (1942). Size 12: Kuosa–Meyrignac–Chen (1999) and Choudhry–Wróblewski (2008).
  - Coppersmith–Mossinghoff–Scheinerman–VanderKam, arXiv 2304.11254: searches to size 16, including Z[i] and Eisenstein integers.
  - Choudhry, arXiv 1908.11192: chains.
- **PTE constants.**
  - Kleiman, J. reine angew. Math. 278/279 (1975);
  - Rees–Smyth, LNM 1415 (1990);
  - Caley, PhD thesis, Waterloo 2012 (bounds on C_n);
  - Caley, Math. Comp. 82 (2013), arXiv 1011.1262, on PTE over Z[i];
  - Filaseta–Markovich (2-adic valuation of C_n).
- **PTE in cryptography.** Costello–Meyer–Naehrig, EUROCRYPT 2021, LNCS 12696 (ePrint 2020/1283): twin-smooth integers from PTE. N-free.
- **VMVT.** Bourgain–Demeter–Guth, Annals 184 (2016); Wooley, Adv. Math. 294 (2016), the cubic case. Constants depend on (M, k).
- **Schnorr–Adleman lattice factoring.**
  - Schnorr, EUROCRYPT '91 (LNCS 547);
  - Adleman 1995 (unpublished manuscript);
  - Vera, arXiv 1003.5461: effective sufficient conditions towards a rigorous reduction, not a failure proof;
  - Ducas 2021 (the SchnorrGate experiments): failure in practice;
  - Aboumrad–Widdows–Kaushik, arXiv 2308.07804: lattice factoring does not scale;
  - Khattar–Yosri, arXiv 2307.09651: a comment on the quantum-assisted sublinear-dimension variant, which fails from about 80 bits.
- **Degree-2 baselines.**
  - The integer-factoring "linear sieve" is Schroeppel's (unpublished; analysed in Pomerance, Math. Centre Tracts 154, 1982).
  - Coppersmith–Odlyzko–Schroeppel, Algorithmica 1 (1986), is the GF(p) discrete-logarithm version. We say "COS" for the shape.
  - Reyneri's cubic sieve (via COS 1986) needs special p (or N).
- **Real roots and fewnomials.**
  - de Gua (1741) and Popoviciu (1935), for GAP-R;
  - Hajós (1953);
  - Lenstra, Number Theory in Progress I (1999), lacunary factors, which is tangential.
- **Elliptic.** Weil pairing; Merel, Invent. Math. 124 (1996); Néron–Tate counts ~c(log A)^{r/2}.
- **Internal.** The Splitting Lemma is systeq §5.4 (internal; not a literature result).

## 7. What remains open (precise)

1. **Near-ideal PTE classes.** Growing M, small t, whose difference polynomials D cancel *systematically* at an N-aligned point m. This is the failure of NC.
2. **Non-star complex-direction families** over a fixed field, and fields with infinite unit groups.
3. **An output-sensitive enumerator** of {P split, |roots| ≤ B} ∩ {P(m) ≡ small mod N} at cost |output|·L[o(1)], with output ≥ L[α] for some α < 1/3. Either log|v| ≤ x^γ with γ < 2/3, or v smooth by identity. This would be a sufficient primitive (smoothenum's algebraic or doubly-smooth listing). Codex notes it need not be necessary: a special arithmetic parametrisation could succeed without solving the general problem.

No candidate is known for any of them.

## 8. Reviews

### 8.1 Codex blind round (ptecons_1)

**Same verdict**, about 2%.

**Adopted:**
- the CW height bound and the support argument;
- the Dickson classification, fixed-field bound and unit density;
- row 2b;
- the (C₂)^r family and its moment rigidity;
- the rank trap;
- split density in coefficient space;
- the γ < 2/3 threshold;
- the elliptic Weil-pairing bound.

**Its experiment** (M = 6, roots in [−12, 12], m = 13) found 12 two-member ideal classes. None had |v| ≤ 1000, because the constants 14400–100800 exceed 2V, as Kleiman predicts.

### 8.2 Codex adversarial round (ptecons_2) — all points applied

- **AVG** was false as a universal statement (squares; lcm + 1). It is now a conditional heuristic.
- **GAP-R/LIN-R** hold; the classical source is Popoviciu 1935.
- **The complex-direction claim** ("forces cyclotomic fields") was false; codex supplied a Z[i] example and a general-k construction. It is replaced by STAR plus an open sub-case.
- **CUN** was false for odd n; corrected, with the even-n bound sharpened to ∏|a| ≤ ‖r‖₁.
- **LEDGER** is not a universal obstruction (F = P − 1).
- **CW's ideal case** is Kleiman 1975, stronger by a factor of 2, with the spacing corollary.
- **HALF** is not a deterministic barrier (zero-root example); A ≳ √(m/M) instead of √m. It is now a typical-member heuristic.
- **"No variant beat COS"** needed a defined criterion (the 0.565 < 0.574 cell). The per-column metric was added.
- **The failure-mode list** is not exhaustive. Reworded, with the open list in §7.
- **The γ < 2/3 threshold** was confirmed. The CW boundary at α = 1/3 is by a log factor.

### 8.3 Opus critics

**Proof/cost critic: ACCEPT WITH FIXES (applied).**
- MAJOR: CUN odd n (fixed).
- MAJOR: PROVED labels on cost consequences (relabelled PROVED bound + HEURISTIC cost).
- MAJOR: §0/§2/§4 inconsistencies about what is open (fixed).
- MAJOR: Kleiman supersedes the ideal CW (applied).
- MAJOR: the E2 metric is confounded (per-column metric and confounds added; claims softened).
- MAJOR: ACS is "no better than" L[1/2] (fixed).
- MINOR, all applied:
  - the P2f formula (now via the median |e₁|);
  - P2d noise;
  - LIN-R needs M ≥ 3;
  - φ is not monotone (removed);
  - NC uses the actual deg D;
  - §3.7 wording;
  - the centred-residue note;
  - the test count and P2a values;
  - the remark on the smoothness-control bucket bias (not re-run).
- The critic verified the CW proof, the SmoothTester, the leakage test and the table numbers against the JSONL.

**Literature critic (applied).**
- Kleiman's year and content.
- The Rees–Smyth reference.
- Caley is Z[i] (the thesis covers the Z bounds).
- Vera mischaracterised; Ducas and Aboumrad–Widdows–Kaushik added.
- Khattar–Yosri scope.
- Adleman is unpublished.
- COS is a DLP paper; Schroeppel's linear sieve; Reyneri special p.
- Size-9 and size-12 attributions.
- Wright 1935.
- BGJT sources and the Granger–Kleinjung–Zumbrägel remark.
- GAP-R = de Gua's rule.
- Pierce–Lehmer for CUN.

No confirmation round was run on this revision.

## 9. Reproduction

```bash
.venv/bin/python -m pytest -q tests/test_ptecons.py
PYTHONPATH=. .venv/bin/python -m src.ptecons.run_experiment E1    # ~2 s
PYTHONPATH=. .venv/bin/python -m src.ptecons.run_experiment E2    # ~3 min
```

The post-hoc per-column metric was computed from the E2 JSONL with π(B) = 6542 and π(max piece).
