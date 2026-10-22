# zroute: deterministic factoring under GRH via the Z route — what is proved, and the one missing lemma

Agent `zroute`, phase 22, 2026-10-22.
- Pre-registration: `notes/zroute_prereg.md` (commit `986f17e`; addendum A1 for E4 in `e66204c`, both before the
  runs they govern).
- Code: `src/zroute/core.py` (operational and deterministic: it sees only N and imports no randomness),
  `src/zroute/harness.py` (harness; the random module is used only by registered controls), `src/zroute/analyze.py`,
  `src/zroute/posthoc_e3.py` (post hoc, labelled). Tests: `tests/test_zroute.py` (10 tests).
- Data: `experiments/raw/zroute_e{1,2,3,4}.jsonl`, `experiments/analysis/zroute_summary.json`; config
  `experiments/configs/zroute.json`.
- GPT-6.1-sol: blind round `codex/zroute_1.{md,out}`; adversarial round `codex/zroute_2.{md,out}` (§9).
- Labels: **proved**; **proved (folklore-level)** for elementary combinations of classical results; **proved from
  A / from U** (conditional on the internal unrefereed chain of `notes/thmA_proof.md` / `notes/lv5unif_results.md`);
  **conditional**; **open**; **observation** (a limitation of a proof template, not a lower bound);
  **experimentally supported**; **heuristic**.

## Verdict (plain language)

**No deterministic subexponential factoring algorithm under GRH was obtained.** The Z route now has everything
except relation supply, and the missing piece can be stated as one explicit algorithmic lemma.

1. **The algebra is done, deterministically [proved].** Lemma XK-Z (§2) turns NFS-type rows (fixed or growing
   degree, composite N, any degree parity) into pure rational relations ∏ q^{e_q} ≡ 1 (mod N). It uses exact
   integer cancellation, archimedean logarithms, LLL and a unit gap. It returns a Q-basis of *everything* that
   exact cancellation can give, not one dependency, with an exact rank formula. It uses no random prime and costs
   L^{O(a)} for B = L^a. Under ERH, Miller's test then factors N deterministically. This needs only the primes below
   2 ln² N to lie in the Q-span of the relations (Prop PIN), not full rank. GPT-6.1-sol proved XK-Z independently in
   its blind round.
2. **Seeds, shape test, recursion and fruitfulness are not needed [proved].** Pure degree-one rows (targets λ with
   known exponent vectors, Gauss split λ ≡ u/v, both halves B-smooth by trial division) feed Miller directly. So do
   seedless two-sided NFS pairs. Number fields are needed only below α = 1/2.
   - In the randomized design, random seeds keep exactly two jobs, supply and rank. Lemma HR: random hidden
     exponent labels make full Q-rank automatic. This is the factoring analogue of dlog's hidden labels.
3. **Rank can be forced deterministically [proved].** Lemma DD1 (§4) takes one row per prime q, with target
   q^H·w and H larger than the total number of prime factors in the row. The relation matrix is then strictly
   diagonally dominant, hence nonsingular, whatever the splits are. Codex's Wronskian-label Lemma WR does the same
   for descent rows after exact elimination. So **rank is not part of the missing lemma.** It remains a hypothesis
   (PIN) only for seedless designs.
4. **The crux, exactly [open].** What is missing is a deterministic subexponential *constructor*:
   - **Hypothesis DSR (deterministic smooth representation).** Given N and a target a = q^H (q ≤ B), find within
     L^{O(1)} time B-smooth u, v ≤ N³ and a small B-smooth multiplier w with u ≡ a·w·v (mod N).
   - **Theorem DZ [proved; ERH]:** DSR ⇒ deterministic L[1/2, O(1)] factoring of every N. The descent version
     DSR-U ⇒ L[2/5, O(1)].
   - **Existence is free under GRH** (Prop GX: every residue class mod N contains an L^a-smooth integer below
     N^{2+ε}). What is missing is constructivity.
   - This mirrors abcsunit ("existence free, constructivity is the failure"). The orchestrator's hope that
     "success is checkable, so a short explicit list suffices" fails for that reason.
5. **Why the known certificates cannot give DSR [observations + one proved template barrier].**
   - **Prop FC (proved, template).** Character-sum certificates of Hölder type certify that an explicit list of m
     candidates hits the good set only in three cases: when m ≳ N^{1/2} (sup on the good set); when m ≳ N/δ
     (Parseval); or when the good set has Fourier algebra norm L^{O(1)}. Subgroup complements have that property,
     which is why Bach-type certificates work for nonresidues; smooth-split sets are not known to.
   - For pair-counting certificates the total lattice enumeration is ≥ 2N/ρ², even with GRH's optimal bias for the
     factor base itself (Prop GX).
   - Averaging over fields, or over N (for "almost all N"), converts the requirement into a level of distribution
     1 − o(1) with only L^{O(1)} elements per residue class. That is beyond Elliott–Halberstam (Observation LV).
   - So "FM/LD* certify only N^{Ω(1)} families" is not a quirk of one proof. It is forced for every
     moment/large-sieve/BV template we know of. It is not a lower bound for all methods.
6. **Literature [§7].** The record is still N^{1/5+o(1)} (Harvey 2021; Harvey–Hittmeir 2022; GRH: Shanks 1971,
   Schoof 1982). Gao–Feng–Hu–Pan (arXiv 2512.19076) improve only log factors for balanced semiprimes;
   Harvey–Hittmeir 2026 (arXiv 2601.11131, ANTS XVII) concerns finding elements of large order.
   - We found no deterministic subexponential algorithm under GRH for factoring, class groups or prime-field
     discrete logarithms, and no hardness result.
   - The nearest GRH-derandomization success is Evdokimov's quasi-polynomial deterministic factoring of polynomials
     over finite fields (1994).
7. **Experiments [§8; toy sizes, correctness only].**
   - E1: the deterministic seedless cubic NFS Z route factored every instance in which PIN was reached, with exact
     relation verification throughout.
   - E2: deterministic orbit targets need 10–19n rows for full rank versus exactly n for random labels.
   - E3: the orbit's yield equals the random yield on average but is over-dispersed (variance/mean 4.3; P3c failed).
     Post hoc, 96% of consecutive orbit hits are the trivially derived split of the previous target, i.e. duplicate
     relations.
   - E4: the dominance design forced full rank and factored every instance it completed.

**Probability assessment (subjective).** Proving DSR, or anything equivalent, with known techniques: ≲ 5% (ours);
codex: ≈ 10% for deterministic L[1/2, O(1)] under GRH by adapting known techniques.

---

## 1. Setting and notation

- N is odd, not a prime power, with no prime factor ≤ max(B, 3 ln² N) (preprocessing: trial division, integer roots,
  AKS). ω(N) ≥ 2.
- B = L^a, where L = L_N[1/2, 1] (or L_N[2/5, 1] where stated).
- FB = (−1, p₁, …, p_{n−1}), the primes ≤ B. Γ = Z^n, ev(e) = ∏ q^{e_q} mod N, Λ = ker ev, U = (Z/N)^×.
- Q₀ := {primes < 2 ln² N} ⊆ FB. Under ERH, Q₀ generates U (Bach 1990; detGRH §1).
- λ(N) is the Carmichael function. "Row" means an integer vector with a certified multiplicative identity mod N.
- detGRH's Propositions Z and EQ are used as stated there (proved, folklore-level).

## 2. Lemma XK-Z: exact elimination of the algebraic side

**Data.**
- f ∈ Z[x] primitive and irreducible of degree d, with leading coefficient t, root α, K = Q(α). An integer m with
  f(m) ≡ 0 (mod N).
- Atoms P = (a_P, b_P), coprime, with F(P) ≠ 0 B-smooth and a_P − b_P m a unit mod N (gcd-checked: otherwise a
  factor).
- Rows i = 1..M, each with a rational part r_i ∈ Z^n, atom exponents y_{iP} and auxiliary exponents c_{ih} (units
  z_h mod N, for example g₁), satisfying
  - ev(r_i)·∏_h z_h^{c_{ih}} ≡ ∏_P ψ(a_P − b_Pα)^{y_{iP}} (mod N),  ψ(a − bα) := a − bm.     (2.1)
  - Seeded rows put the target's exponent vector into r_i. Seedless NFS rows (a − bm and F(a, b) both smooth) have
    r_i = e(a − bm), y_i = e_P, no auxiliary part. Theorem A rows have the g₁ column.
- **Algebraic matrix A** (M × C_alg):
  - the auxiliary columns;
  - the column Σ_P y_{iP} (the exponent of 𝔞 = O_K + αO_K; it is Theorem A's g₁ column);
  - for every prime ideal 𝔯 above r ≤ B, the column Σ_P y_{iP} v_𝔯(𝔍_P), with 𝔍_P := (a_P − b_Pα)𝔞^{−1}.
- Λ₁ := ker_Z(Aᵀ); for x ∈ Λ₁, γ_x := ∏_P (a_P − b_Pα)^{Σ_i x_i y_{iP}}; Λ₂ := {x ∈ Λ₁ : γ_x ∈ μ_K};
  W := span_Q{Σ_i x_i r_i : x ∈ Λ₂}.

**Lemma Ψ (evaluation at composite N) [proved].**
- Because f is primitive, ker(Z[x] → K, x ↦ α) = fQ[x] ∩ Z[x] = fZ[x] (Gauss's lemma). So x ↦ m defines a ring
  homomorphism ψ: Z[α] → Z/N.
- Let S := ψ^{−1}(U). ψ extends uniquely to S^{−1}Z[α] → Z/N (universal property of localisation).
- The group 𝒜 generated by the atoms lies in (S^{−1}Z[α])^×, and ψ|_𝒜 : 𝒜 → U is a homomorphism. So two products
  of atoms that are equal in K have equal images mod N.
- No reduction modulo a prime ideal above N is used, so composite N, repeated roots of f modulo divisors of N and
  non-integral α are all harmless.
- (codex zroute_1 §1.2 states the same; the dlog remark (iii) sketch is the special case.)

**Lemma XK-Z [proved].**
- (i) For x ∈ Λ₂, o(x) := w_d·Σ_i x_i r_i satisfies ev(o(x)) = 1, where w_d := lcm(1, …, 2d²) is a multiple of |μ_K|.
  For odd d, μ_K = {±1}, and o(x) := Σ_i x_i r_i + εe_{−1} with ε ∈ {0, 1} read off from ev(Σ x_i r_i) = ±1.
- (ii) A deterministic algorithm computes integer vectors whose Q-span is W, and their number dim W.
- (iii) dim W = rank_R([V; T]) − rank_R(V), where Z is any integer basis of Λ₁ ⊗ Q, T := RᵀZ and V := (Log γ_z)_{z ∈ Z}
  (codex). In particular dim R^T(Λ₁ ⊗ Q) − r_U ≤ dim W ≤ dim Rᵀ(Λ₁ ⊗ Q), with r_U = r₁ + r₂ − 1 ≤ d − 1.
- (iv) **Cost.** Let M, C_alg ≤ L^{a+o(1)}, entries of A and r_i ≤ (log N)^{O(1)}, and d ≤ log N. Then the deterministic
  cost is L^{(ω+1)a+o(1)} with fast matrix multiplication, or L^{4a+o(1)} classically, plus L^{3a+o(1)} for the unit
  step. Outputs have L^{a+o(1)} bits. This holds at fixed d and at d ≍ (log N/log log N)^{1/5} (L = L_N[2/5, 1]).

*Proof.*
- **(i) Correctness.**
  - For x ∈ Λ₁, γ_x ∈ O_K^×. This is dlog Lemma XK1: the 𝔞-column gives Σ_P y_P = 0, the ideal columns cancel
    every 𝔯 | 𝔍_P, and N(𝔍_P) = |F(P)| is B-smooth, so (γ_x) = O_K. Its proof is degree-free (dlog25 §4).
  - Bad primes r | t·disc f: XK0(iii)/XK0-U computes r-maximal orders and valuations.
    - The only randomness there is factoring f mod r. Replace it by Berlekamp's deterministic algorithm, poly(d)·r
      operations with r ≤ B. There are O(log N) bad primes, so this costs L^{a+o(1)}.
  - Multiply (2.1) over i with exponents x_i. The auxiliary columns are columns of A, so they cancel. Lemma Ψ then
    gives ev(Σ x_i r_i) ≡ ψ(γ_x) (mod N).
  - For x ∈ Λ₂, γ_x = ζ ∈ μ_K ∩ 𝒜, so ψ(γ_x)^{|μ_K|} = 1.
  - If ζ ∈ K has order w then φ(w) ≤ d, and φ(w) ≥ √(w/2) gives w ≤ 2d². So w_d kills μ_K. For odd d, μ_K = {±1}.
- **(ii) Deterministic Q-basis of Λ₁.**
  - Hadamard: every minor of A has |Δ| ≤ Had := (√C_alg·‖A‖_∞)^{C_alg}.
  - Take the first primes p₁ < p₂ < … above 2^{62} until ∏p_k > Had. Then rank_Q A = max_k rank_{p_k} A, because a
    fixed nonzero maximal minor cannot vanish modulo all of them. A p_k attaining the maximum yields a rank profile
    (I, J) with det A_{IJ} ≠ 0.
  - For i ∉ I, put x^{(i)} := det(A_{IJ})·(−A_{i,J}A_{IJ}^{−1}, e_i). This is integral by Cramer, computed by CRT or
    by p_k-adic lifting (deterministic, since p_k ∤ det A_{IJ}), with entries ≤ Had.
  - Cost: L^{a+o(1)} primes × ρ^ω, plus the products: L^{(ω+1)a+o(1)}. Classical Bareiss: L^{4a+o(1)}.
- **(ii) Unit step: Lemma INC [proved; deterministic incremental form of dlog XK4 / dlog25 XK4-U].**
  - Process x^{(i)} in order, keeping a pivot set I′ whose Log vectors are R-independent (|I′| ≤ r_U).
  - For each new j, run XK4's LLL on {Log γ_{x^{(i)}} : i ∈ I′} ∪ {Log γ_{x^{(j)}}} (dimension ≤ 2d), at XK4-U's
    precision κ.
  - **Case 1: a relation exists** (Log γ_{x^{(j)}} ∈ span(I′)). Then a small integer relation exists: the Cramer
    bound R of XK4-U applies, since the vectors lie in the unit lattice and rational dependence equals integer
    dependence there.
    - XK4's argument then shows that LLL's first vector is a relation. Its j-coefficient is nonzero because the
      pivots are independent.
    - It is certified by ‖Σ n Log γ‖_∞ < β_d/2 in ball arithmetic. Here β_d is the Blanksby–Montgomery gap, or
      codex's elementary gap δ_d = 1/(d T_d), T_d = d(6·2^d + 1)^d + 1, with log(1/δ_d) = O(d²).
    - Output z^{(j)} := Σ n_i x^{(i)} + n_j x^{(j)} ∈ Λ₂.
  - **Case 2: no relation exists.** Then every LLL vector has ‖Σ n Log γ‖_∞ ≥ β_d − (rounding) > β_d/2, the
    certificate fails, and j joins I′.
  - The two cases are decided rigorously by the gap.
  - The z^{(j)} (j ∉ I′) are independent, since each contains its own x^{(j)}. There are dim Λ₁ − rank(Log image)
    = dim Λ₂ of them, so they form a Q-basis of Λ₂ ⊗ Q.
  - Cost per call: poly(d)·Õ(M(P + β_x) + κ²), with κ, P = O(dβ_x + d² + d log(kMd log N)) and β_x = L^{a+o(1)} (XK4-U
    ledger). M calls: L^{3a+o(1)}.
  - Root isolation (Schönhage/Pan; dlog25 §4) is deterministic.
- **(iii)** is linear algebra: Λ₂ ⊗ Q = ker of the composite Q^{dim Λ₁} → (Log image).
- **(iv)** sums the above. All poly(d, k, log N) factors are L^{o(1)}. ∎

*Checks.*
- E1 (§8): every one of the output relations verified exactly, with 0 non-unit kernel vectors and 0 undecided
  certifications.
- Test `test_xkz_relations_verify_and_rank_formula`: torsion dimension = dim Λ₁ − pivots.

*Prior art.* This is dlog's Proposition XK (itself partially pre-empted by Buchmann–Pohst LLL unit relations and
Schirokauer's explicit units) applied to all of Λ₂ instead of one vector, with the random prime removed.
GPT-6.1-sol proved the same lemma blind (codex zroute_1 §1). We claim nothing beyond "no prior formulation found by us".

**Proposition PIN (Miller under partial rank) [proved; ERH; folklore-level].**
- *Statement.* Let L_c ⊆ Λ be any relation lattice containing 2e_{−1}, and suppose e_q ∈ L_c ⊗ Q for every q ∈ Q₀.
  Let P be the pivot columns of an HNF basis H of L_c, and D := |det H_P|. Then q^D ≡ 1 (mod N) for all q ∈ Q₀,
  λ(N) | D, and Miller's test with bases a = 2, …, ⌈2 ln² N⌉ splits N.
- *Proof.*
  - e_q = cH over Q, and restricting to P gives c = (e_q)_P H_P^{−1}, whose denominators divide D. So De_q ∈ L_c ⊆ Λ.
  - Q₀ generates U under ERH, so λ(N) | D.
  - The rest is detGRH Prop Z (Miller 1976; Bach 1990). ∎
- *Remarks.*
  - PIN concerns only O(log² N/log log N) coordinates. Full rank on large factor-base primes is unnecessary.
  - Projection rank is not enough (codex): the other coordinates must cancel. Codex also notes that no single
    numerical rank condition is both necessary and sufficient.
  - In E1 and E2, Miller often split N before PIN held (E1: see §8), because the pinned subgroup already contained a
    witness. This is luck, not a guarantee.

## 3. Randomness sources revisited: what the Z route still needs

| detGRH row | randomness | status in the Z route |
|---|---|---|
| 1–3 | field law (g₀, g₁, w, c) | **not needed for α = 1/2**: the degree-one design (§3.2) uses no number field. Needed for α < 1/2 (descent); there it is supply-only (§5) |
| 4a | seed law, supply role | **the crux**: replaced by Hypothesis DSR (§5) |
| 4b | seed law, fruitfulness role | gone (Miller under ERH) |
| — | seed law, **rank role** | random labels give rank automatically (Lemma HR); deterministic labels force it (Lemmas DD1, WR, §4) |
| 6 | shape test | not needed for α = 1/2 (halves are tested for B-smoothness by trial division). For α = 2/5 the descent needs complete factorisations of the halves: recursion, which is deterministic if the outer algorithm is |
| 7, 7′ | kernel sampling, exact LA | deterministic (Lemma XK-Z (ii)) |
| 8, 9 | characters, Hensel root | gone (no square roots) |
| 10–13 | — | deterministic |

### 3.1 Random seeds: exactly two jobs remain

**Lemma HR (hidden labels force rank) [proved].**
- *Setting.* Let s₁, …, s_M be i.i.d. uniform on [0, X)^n with X ≥ N⁶, and T_i := ev(s_i). Let 𝒯 be any random
  variable such that (s_i)_i is conditionally independent of 𝒯 given (T_i)_i. (For example, 𝒯 is all targets,
  descents, field data and the coins of every other step.)
- Let e_i ∈ Z^n be 𝒯-measurable, put r_i := s_i − e_i, and let z^{(1)}, …, z^{(K)} ∈ Z^M (K ≥ n) be 𝒯-measurable and
  linearly independent, with entries of at most N^{1/2} bits.
- *Claim.* Pr[rank_Q{Σ_i z^{(k)}_i r_i}_k < n | 𝒯] ≤ 2/N + 4MnN³/X.
- *Proof.*
  - **Choose ℓ′.** Fix (given 𝒯) a prime ℓ′ ∈ (N, 2N] not dividing a nonzero K × K minor Δ of the z-matrix. One
    exists, since that interval has ≥ N/(2 ln N) primes and Δ has < N^{1/2}·K + K log K bits, so it has fewer prime
    divisors above N.
  - **The labels mod ℓ′ are nearly uniform given 𝒯.** Put m₀ := ℓ′λ(N); gcd(ℓ′, λ(N)) = 1 since λ(N) < N < ℓ′.
    - Let s′ be uniform on [0, m₀⌊X/m₀⌋)^n. Then TV(s, s′) ≤ nm₀/X.
    - s′ mod m₀ is uniform on (Z/ℓ′)^n × (Z/λ(N))^n (CRT), and ev(s′) depends only on s′ mod λ(N) (λ(N) is even).
      So s′ mod ℓ′ is uniform and independent of ev(s′).
    - Pr[ev(s′) = T] = 1/|G_B| ≥ 1/N for T ∈ G_B. Conditioning on ev(s) = T inflates the distance to at most
      2nm₀N/X ≤ 4nN³/X.
    - By conditional independence, given 𝒯 the s_i mod ℓ′ are independent and each is within 4nN³/X of uniform.
  - **Rank.** z mod ℓ′ has full row rank K, so S ↦ ZS mod ℓ′ maps the uniform law on F_{ℓ′}^{M×n} onto the uniform
    law on F_{ℓ′}^{K×n}.
    - A uniform K × n matrix (K ≥ n) has rank < n with probability ≤ Σ_{j≤n} ℓ′^{j−1−K} ≤ 2/ℓ′ < 2/N.
    - Subtracting the fixed matrix ZE does not change this.
    - Rank mod ℓ′ ≤ rank over Q. ∎
- *Remarks.*
  - Codex zroute_1 §3.3 gives an equivalent proof (an affine-hyperplane bound per label coordinate).
  - The labels must stay in the rational block. Treating them as separate auxiliary columns cancels the signal.
  - HR does not give GF(2) rank: conditioning on T can fix parities.

**Proposition ZA (a character-free, fruitfulness-free Theorem A) [proved from A; ERH].**
- *Statement.* In Theorem A, replace:
  - R = r² by T = ev(s), s uniform on [0, N⁶)^n;
  - Steps 8–9 (LV6, characters, Hensel root, Theorem DD) by Lemma XK-Z, Lemma HR and Miller with deterministic bases.
  The result is a Las Vegas algorithm under ERH with expected time L[1/2, c_DL(θ) + o(1)], where
  c_DL(θ) = min_a max(a + 1/√2, a + θ/(2a), 3a) (dlog Theorem DL's constant, = (1 + θ)/√2 for θ ≤ 1/2).
- *Proof sketch (the only new step is supply).*
  - Under ERH, G_B = U. T is then within nN/N⁶ of uniform on U (codex zroute_1 §2.5).
  - dlog Lemma S_p holds at composite N with the same proof:
    - the coprime-pair map (u, v) ↦ uv^{−1} is injective because |uv′ − u′v| < Y² < N;
    - the halves are units because their primes exceed B;
    - the Gauss-reduced vector is the split.
    So Theorem K(ii) holds without the 2^{−ω} factor.
  - Acceptance and the descent depend on s only through T, so Lemma HR applies with M ≥ n + C_col + d + n rows.
    XK-Z costs L^{3a} (Las Vegas) or L^{(ω+1)a} (deterministic).
  - Seed evaluation costs L^{a+o(1)} per trial, which is dominated.
  - The same substitution in Theorem U (L[2/5]) is "proved from U" with dlog25's XK-U.
- *Status.* This is a re-derivation, not an improvement. Its point is diagnostic: in the Z route, random seeds do
  supply and rank, nothing else.

### 3.2 Seedless and field-free designs

- **Degree-one rows [proved].** Take λ ∈ U with known exponent vector x(λ) and the Gauss-reduced (u, v) with
  u ≡ λv. If u and v are B-smooth, then x(λ) − e(u) + e(v) ∈ Λ.
  - The smoothness test is trial division (or a remainder tree) by FB, deterministic, L^{a+o(1)}.
  - There is no shape test, no descent and no recursion. This is the classical Dixon-type row, read over Z.
  - Small explicit targets are useless (detGRH Lemma SR; codex: λ = 2 with |u|, |v| < N/3 forces u = 2v). Targets
    must wrap around modulo N.
- **Seedless two-sided NFS [proved].** A pair (a, b) with a − bm and F(a, b) B-smooth is a row of type (2.1). No
  target, seed, shape test or square root is needed.
  - Rank then has no labels to lean on, and PIN is a genuine hypothesis (E1 shows it failing for lack of
    occurrences; §8).
- **Pooling across fields [proved; codex].** Rational outputs of different fields add (W_total = Σ W_j), but each
  field must be eliminated with its own rows. dim W_j ≤ M_j − rank(A_j) − (unit rank). A union-over-fields
  descendability statement does not by itself give rank: codex gives J fields with one row each, every target
  covered, and all kernels zero.

## 4. Deterministic rank forcing

**Lemma DD1 (diagonal dominance for split rows) [proved].**
- *Construction.* For each prime q ≤ B choose any B-smooth w_q, and suppose the target a_q := q^H·w_q mod N has a
  representation u_q ≡ a_q v_q with u_q, v_q B-smooth. Let K := max_q Ω(w_q) and U* := max_q log₂|u_qv_q|. If
  H > K + U*, then the n − 1 rows ρ_q := He_q + e(w_q) − e(u_q) + e(v_q), together with 2e_{−1}, have full rank n.
  - For Gauss-reduced splits, |u_qv_q| ≤ (u² + v²)/2 ≤ N/√3, so H := K + bitlen(2N) + 2 suffices. For
    representations with |u|, |v| ≤ N³, take H > K + 6 log₂ N.
- *Proof.*
  - On the prime columns, the diagonal entry is ≥ H − v_q(u_q) and the off-diagonal row sum is
    ≤ Ω(w_q) + Ω(u_q) + Ω(v_q) − v_q(u_q).
  - So the block is strictly row diagonally dominant, and nonsingular (Lévy–Desplanques).
  - Adding the row 2e_{−1} multiplies the determinant by ±2. ∎
- *Consequence (ERH).* D := |det| is a multiple of φ(N), and Miller splits N.
- *Credit.* The dominance idea is codex zroute_1 §5, stated there for integer representatives s_q ≡ q^H. The
  extension to split rows with multipliers w is ours.
- *Checks.* E4 (§8): every completed instance was dominant and of full rank, and Miller split all of them.

**Lemma WR (Wronskian labels; codex zroute_1 §3.4) [proved; checked by us].**
- *Construction.* Index the candidate targets by j = 1..P and give target j the label
  x_{j,i} := (n + j)!/(n + j − i)!·t^{n+j−i} (0 ≤ i < n), the i-th derivative of t^{n+j} at an integer t.
- *Claim.* Let t := C + 2, with C := n!(P + 1)^n A₀^n and A₀ := H₀((n + P)^n + PB₀). Then for every accepted subset,
  every kernel basis Z (rank ≥ n, entries ≤ H₀) and every descent part V (entries ≤ B₀), the outputs Zᵀ(X − V) have
  rank n.
- *Proof (checked).*
  - Choose n independent kernel vectors. The output matrix is [p_l^{(r)}(t) − c_{lr}], with
    p_l(t) := Σ_i z^{(l)}_i t^{n+j_i}. The j_i are distinct, so the p_l are independent polynomials.
  - After an invertible constant change of the l-basis their degrees d_l are distinct and > n.
  - The Wronskian has degree D = Σd_l − n(n − 1)/2 and leading coefficient ∏ lc(p_l)·∏_{a<b}(d_a − d_b) ≠ 0.
  - Replacing r derivative entries by constants lowers the degree strictly, because Σ_{replaced} d_l ≥ r(n + 1)
    > r(2n − r − 1)/2.
  - So the determinant is a nonzero integer polynomial with coefficients ≤ C, and t = C + 2 exceeds all its roots
    (Cauchy).
  - The bound is uniform over all data, so V may depend on t through the targets. ∎
- *Notes.*
  - Labels have L^{O(1)} bits when P, n, H₀, B₀ are L^{O(1)}-bounded.
  - It applies after XK-Z, to descent rows, where dominance is lost through the combinations z.
  - Test `test_wronskian_labels_force_rank` checks 300 random instances.

**Consequences for the rank question.**
- In any seeded design the algorithm chooses the labels. Rank is then forced deterministically once the exact torsion
  kernel has dimension ≥ n (enough rows). It is no longer part of the open problem.
- In seedless designs (two-sided NFS, intervals near kN) PIN remains a hypothesis. E1 and E2 show its mechanism at toy
  scale: a prime q is pinned only once it occurs, which takes ≈ q·log(#small primes) rows.
- **Relation to rkrank.**
  - RK (GF(2) rank) ⇒ full Q-rank ⇒ PIN, since rank_Q L ≥ rank_{F₂}(L mod 2).
  - Theorem V_κ's codimension-D^κ deficiency does not imply PIN.
  - Conversely, full Q-rank can hold with zero F₂ information: the doubled lattice (detGRH E1; codex §3.2).
  - The Z route therefore needs strictly less than RK. With labels it needs nothing at all on rank.

## 5. Deterministic relation supply: the crux

### 5.1 The conditional theorem and the missing lemma

**Hypothesis DSR(a, c) (deterministic smooth representation) [open].**
- *Statement.* There is a deterministic algorithm which, given N (preprocessed) and a target a ∈ {q^H mod N : q ≤ B}
  with H := 8⌈log₂ N⌉ and B = L^a, runs in time L_N[1/2, c] and outputs B-smooth integers u, v with |u|, |v| ≤ N³,
  and a B-smooth w with Ω(w) ≤ log₂ N, all with their factorisations, such that u ≡ a·w·v (mod N).
- *Variants.*
  - **DS-T:** among an explicit list of L^{c} targets with Wronskian labels, at least n + 1 split with both halves
    B-smooth.
  - **Smooth compression** (codex zroute_1 §5): v = w = 1 and an integer representative u ≤ N³ of q^H.
  - **DSR-U:** the same at L = L[2/5], with the representation obtained by special-q descent in explicit fields.

**Theorem DZ [proved; ERH; the composition is folklore-level].** ERH + DSR(a, c) ⇒ a deterministic algorithm that
completely factors every N in time L_N[1/2, max(c + a, (ω + 1)a) + o(1)].
- *Proof.*
  - One DSR call per prime q ≤ B gives the rows of Lemma DD1, since Ω(u) + Ω(v) + Ω(w) ≤ 7 log₂ N < H.
  - Full rank, D = |det| (deterministic, L^{(ω+1)a+o(1)}), Prop PIN/Z and Miller split N.
  - Recursion on the factors (AKS for primes, integer roots for powers) takes ≤ log₂ N calls, each at most the bound
    for N.
- Likewise ERH + DSR-U ⇒ deterministic L[2/5, O(1)]. Here the halves of a descent need complete factorisations, by
  recursion on the deterministic algorithm itself (recurse Lemma F_α), and XK-Z handles growing degree. ∎

**Proposition GX (existence is free under GRH) [proved; GRH for Dirichlet L-functions mod N].**
- *Statement.* Let y = L^a with a > 0 fixed. For N large, every a ∈ U has a y-smooth representative u ≤ N^{2+ε}. In
  particular DSR's required output exists, with v = w = 1.
- *Proof (codex zroute_1 §4.5, re-derived).*
  - Let μ be uniform on the primes ≤ y. Under GRH, |Σ_{p≤y} χ(p)| ≪ √y log²(Ny) for χ ≠ χ₀ mod N (from
    ψ(y, χ) ≪ √y log²(Ny); prime powers contribute O(√y)). So |μ̂(χ)| ≤ y^{−1/2+o(1)}.
  - μ^{*k}(a) = φ(N)^{−1}Σ_χ χ̄(a)μ̂(χ)^k ≥ 1/φ(N) − max_{χ≠χ₀}|μ̂(χ)|^k > 0 once k ≥ (2 + ε/2)·ln N/ln y.
  - The product of k primes ≤ y is ≤ y^k ≤ N^{2+ε}. ∎
- *A pair-counting variant (ours).*
  - Under GRH, for y = L^{1/a+δ} ≤ B (a > 1), the number of (n, u, v) with n ≤ y, u, v B-smooth in [1, N/2] and
    u ≡ nv (mod N), weighted by Λ(n), is
    ≥ ψ(y)|S|²/φ(N) − C√y log²(Ny)|S| ≥ (1/2)ψ(y)|S|²/φ(N) = N·L^{O(1)}.
  - *Proof:* orthogonality; the error is ≤ max_{χ≠χ₀}|ψ(y, χ)|·#{u ≡ u′} = max|ψ(y, χ)|·|S|.
  - So the factor base is itself a GRH-optimal ε-biased set: its bias y^{−1/2+o(1)} meets the Parseval floor
    m^{−1/2}.
- **Constructivity is exactly what is missing.** Finding the representative needs a search over N^{1+o(1)} candidates
  (k-fold prime products, or the lattice {u ≡ nv} in [0, N/2]²).

**On how hard DSR is (remark, not a claim).**
- For a prime modulus p, DSR supplies a full-rank relation lattice for F_p^* and representations of prescribed
  residues. These are the relation-collection and individual-descent steps of a deterministic index calculus.
- We found no deterministic subexponential prime-field discrete-logarithm algorithm under GRH (§7). DSR is therefore
  plausibly at least as hard as that open problem.
- We did not prove an implication. Turning the relation lattice into logarithms modulo p − 1 needs its saturation,
  which is a Pohlig–Hellman-type step we did not resolve.
- Prop EQ (detGRH) shows deterministic factoring needs only *arbitrary* full-rank relations. DSR asks for more.

### 5.2 Why known certificates do not give DSR

The orchestrator's question is whether an explicit list of L^{O(1)} fields, seeds or targets can be certified to
contain good ones, and whether detGRH's "N^{Ω(1)} families" is a proof artifact.

**Proposition FC (Fourier-certificate barrier) [proved; a statement about certificate templates].**
- *Setup.* G = U, φ = |G|, W ⊂ G the explicit candidates (|W| = m ≤ φ/2), S ⊂ G the good set (|S| = δφ, δ ≤ 1/4).
  Then #(W ∩ S) = mδ + E with E = φ^{−1}Σ_{χ≠χ₀} Ŵ(χ)\overline{Ŝ(χ)}. A Hölder certificate proves #(W ∩ S) > 0 by
  showing mδ > ‖Ŵ‖′_p‖Ŝ‖′_{p′} (primes: χ₀ omitted, normalised by φ^{−1}).
- (a) **(2, 2) (Parseval only) fails unless m > (3/8)φ/δ.**
  - *Proof:* by Parseval, φ^{−1}Σ_{χ≠χ₀}|Ŵ|² = m − m²/φ ≥ m/2 and φ^{−1}Σ_{χ≠χ₀}|Ŝ|² = δφ(1 − δ) ≥ (3/4)δφ.
- (b) **(1, ∞) (sup on S) fails unless m > (√3/4)(φ/δ)^{1/2}, for every S.**
  - *Proof:* φ^{−1}Σ_{χ≠χ₀}|Ŵ| ≥ φ^{−1}Σ|Ŵ|²/max|Ŵ| ≥ (m/2)/m = 1/2, and max_{χ≠χ₀}|Ŝ| ≥ (δ(1 − δ)φ)^{1/2}
    (Parseval).
- (c) **(∞, 1) (sup on W) succeeds only if the Fourier-algebra norm ‖1_S‖_A = φ^{−1}Σ_χ|Ŝ(χ)| ≤ δ(1 + (2m)^{1/2}).**
  - *Proof:* max_{χ≠χ₀}|Ŵ| ≥ (m/2)^{1/2} (Parseval).
  - A random-like S has ‖1_S‖_A ≍ (δφ)^{1/2}, which would need m ≳ φ/δ.
  - Cosets of subgroups have ‖1_S‖_A = 1. This is exactly why ERH certificates for nonresidues or for subgroup
    escape (Ankeny, Bach) need only polylogarithmically many candidates (codex's correction: "every character-sum
    certificate needs N^{Ω(1)} samples" is false in general).
  - For smooth-split sets no bound ‖1_S‖_A = L^{O(1)} is known. Heuristically (saddle-point estimates for twisted
    smooth sums), |Ŝ(χ)| ≈ |S|e^{−O(u)} for many χ, which suggests a large norm. **[heuristic]**
- (d) **Pair certificates.** Here S = smooth integers in [1, Y] and the hits are pairs (u, v) with u ≡ wv, w ∈ W. The
  (∞, 2) certificate needs max|Ŵ| < m|S|/φ, hence m > φ²/(2Y²ρ_Y²). Enumerating the lattice {u ≡ wv} ∩ [−Y, Y]²
  costs ≥ 4Y²/N per w. So the total work is ≥ 2N/ρ_Y² for every Y ≥ √N. GRH (Prop GX) makes the W side optimal; the
  cost is on the S side. ∎
- *What this does not exclude.* Certificates that use the joint structure of Ŵ and Ŝ, not norms. And non-Fourier
  arguments.

**Observation LV (averaging over fields or over N is a level-of-distribution-1 problem) [observation].**
- **Fields.**
  - For a fixed atom P, FM's field average realises the progression F̂(P) + G(P)·j (modulus t = N^{μ+ν}, length
    X_F/t), and the moments need level ≤ 1/2 in t (LV Lemma 8.5′, Harper, Drappeau).
  - An explicit consecutive family of L^{O(1)} fields realises a progression of length L^{O(1)}: level
    log t/log(tJ) → 1. Individual-progression theorems need q ≤ y^{4√e−δ} (Harper 2012). Averaged ones stop at
    level 1/2 (Fouvry–Tenenbaum; 3/5 for special weights, Drappeau 2015).
  - So a short family needs smooth numbers equidistributed at level 1 − o(1) with L^{O(1)} terms per class. That is
    beyond Elliott–Halberstam, which keeps x/q ≥ x^ε. Friedlander–Granville show EH-type statements for primes fail
    at x/q = (log x)^B; the analogue for smooth numbers at x/q = L^{O(1)} is open.
- **N ("almost all N").**
  - For an explicit integer λ (say 2^e with 2^e > X) and N ∈ [X, 2X], relations are pairs of smooth u, v ≤ Y with
    N | u − λv.
  - The dual large sieve gives #{N ∈ [X, 2X] : no relation} ≤ X(Z + 4X²)/|𝒜|, where 𝒜 is the multiset {u − λv}
    (injective for λ > 2Y) and Z ≈ λY its range.
  - Beating X needs |𝒜| ≫ X², i.e. per-N enumeration Y²/N ≥ X/ρ²: trivial.
  - Per-N cost N^{ε} would need level 1 − ε (EH-type for this bilinear sequence). Subexponential cost needs level
    1 − o(1).
- **Answer to the orchestrator's question.**
  - "Current proofs certify only N^{Ω(1)}-size families" is forced for every moment, large-sieve and BV-type
    template we know of (Prop FC; Observation LV).
  - It is not a lower bound. A proof of DSR would need either a pointwise mechanism or the level-1 statement.
  - At growing degree, codex notes that some families are N^{o(1)} yet still ≫ L[2/5]. The same conclusion holds in
    L-units.

**Observation (checkability) [proved, trivial; detGRH §6].** Global moments are compatible with total failure on any
prescribed list. Any proof must use the arithmetic of the list.

### 5.3 Remarks

- **Elements of large order.** Orbit designs (λ_j = g^j) need ord_N(g) > (number of targets). For general N,
  deterministically finding g of large order is itself the subject of Harvey–Hittmeir 2026 (arXiv 2601.11131).
  - Multi-base boxes avoid the issue: a collision ∏q^{j} ≡ ∏q^{j′} of two targets is itself a relation.
  - E3 found no orbit cycle within 2000 targets for 1000 semiprimes.
- **Redundancy of multiplicatively closed target families** [proved; trivial, Lemma SC-type]. If λ ≡ u/v, then
  qλ ≡ (qu)/v gives the identical row x(λ) + e_q − e(qu) + e(v). E3/post hoc: 96% of consecutive orbit hits were such
  duplicates. Deterministic target families should be multiplicatively spread, e.g. Wronskian or q^H·w labels.
- **Deterministic heuristic L[1/3] (remark, heuristic input).**
  - Base-m NFS with growing degree is deterministic. XK-Z (L^{O(a)} at growing degree) plus labelled targets or PIN
    plus Miller give a deterministic algorithm whose only unproved input is the usual NFS yield heuristic (and PIN
    for seedless rows). No fruitfulness heuristic is used.
  - The dense exact elimination worsens the constant but keeps α = 1/3. We make no claim beyond this; it is
    heuristic.

## 6. Answers to the task items

| item | answer | status |
|---|---|---|
| 1 XK-Z | Lemma XK-Z (deterministic, composite N, all degrees, full Q-basis of the torsion kernel, rank formula, cost L^{O(a)}); Miller under PIN | proved (codex independently) |
| 2a fields | not needed at α = 1/2; at α < 1/2 the short-family statement is a level-1 problem; Fourier/moment templates certify only N^{Ω(1)}-size lists; not a lower bound | observations + Prop FC (proved template) |
| 2b seeds | not needed (degree-one or seedless rows); random seeds do supply + rank only; deterministic targets must wrap around (Lemma SR) and be multiplicatively spread | proved |
| 2c shape test | unnecessary at α = 1/2; at 2/5, recursion is deterministic given a deterministic outer algorithm | proved / conditional |
| 3 rank | PIN suffices; random labels force rank (HR); deterministic labels force rank (DD1, WR); only seedless designs need a rank hypothesis; Z route needs strictly less than RK | proved |
| 4 literature | record N^{1/5}; no deterministic subexponential factoring, class group or F_p DL under GRH found; no hardness known | search result |
| 5 experiments | §8 | experimentally supported (toy) |
| crux | Hypothesis DSR (constructive smooth representation of q^H); existence free under GRH (Prop GX); Theorem DZ: ERH + DSR ⇒ deterministic L[1/2, O(1)] | open / proved implication |

## 7. Literature

Verified by me where stated; codex zroute_1 §4 verified the theorem statements marked (c).
- **Deterministic factoring.**
  - Pollard 1974 / Strassen 1976: N^{1/4+o(1)}.
  - Under GRH: Shanks 1971 and Schoof 1982 §6, N^{1/5+ε} (c; attribution as in detGRH).
  - Unconditional:
    - Hittmeir 2018: N^{2/9+o(1)}.
    - Harvey, Math. Comp. 90 (2021) 2937–2950: O(N^{1/5}log^{16/5}N) (c).
    - Harvey–Hittmeir, Math. Comp. 91 (2022) 1367–1379: a further (log log N)^{3/5} saving (c).
  - Gao–Feng–Hu–Pan, arXiv 2512.19076 (unrefereed): O(N^{1/5}log^{13/5}N/(log log N)^{3/5}) for balanced semiprimes.
    The same exponent; they mention N^{1/6}, N^{1/8} only as possibilities (read by me, abstract and Thm 1.1).
  - Harvey–Hittmeir, "Deterministic methods for finding elements of large multiplicative order", arXiv 2601.11131,
    ANTS XVII (2026): removes the hypothesis D ≥ N^{1/6} in the large-order subroutine; no new factoring exponent
    (read by me, abstract and introduction).
  - **The record N^{1/5+o(1)} is confirmed.** We found no deterministic L[α, O(1)] factoring for any α < 1 under GRH,
    and no hardness or impossibility result.
- **Factoring from exponents (ERH):** Miller, JCSS 13 (1976) 300–317; Bach 1984 §3 and Math. Comp. 55 (1990);
  Shoup §§10.4–10.5. Żrałek 2010 (deterministic L[1/3] given iterated φ; an oracle result).
- **Z route, prior formulations:** McCurley 1988; Stange, arXiv 2211.06821 (heuristic, randomized).
- **Class groups.**
  - Hafner–McCurley, JAMS 2 (1989): Las Vegas, ERH.
  - Biasse–Erukulangara 2023: randomized, GRH.
  - de Boer–Pellet-Mary–Wesolowski, arXiv 2512.01588: rigorous under ERH, randomized (c).
  - Buchmann 1990 and Biasse–Fieker: heuristic.
  - Deterministic: Lenstra 1982 / Schoof 1982 |D|^{1/5+ε} (class number/regulator); Schoof 2008 (Arakelov class
    groups, deterministic exponential) (c).
  - No deterministic subexponential class-group algorithm under GRH found.
- **Discrete logs in F_p^*:** Pomerance 1987 (randomized L[1/2, √2]). No deterministic subexponential algorithm
  under GRH found.
- **Deterministic algorithms under GRH that do exist:**
  - square roots mod p (Ankeny 1952 + Tonelli–Shanks; Bach's constants);
  - Miller's primality test (superseded unconditionally by AKS, Annals 160 (2004));
  - polynomial factoring over F_q: Rónyai 1988/1992, Huang (abelian Galois groups), and **Evdokimov 1994
    (ANTS-I, LNCS 877): deterministic (n^{log n} log q)^{O(1)} under GRH**, later Ivanyos–Karpinski–Saxena. This is
    the closest GRH-derandomization success. It works because the search space (algebraic structures of size
    polynomial in n) is small, unlike smooth representatives mod N.
- **Smooth numbers.**
  - Single progressions: Granville 1993; Soundararajan 2008; Harper, JNT 132 (2012): q ≤ y^{4√e−δ},
    log x/log q → ∞ (c).
  - Short intervals under RH: Soundararajan, arXiv 1009.1591 (intervals ≍ u√x/ρ(u/2) for y ≥ exp(5√(log x log log x)))
    (c); Younis, arXiv 2409.05761 (θ > 1/2) (c).
  - Unconditional short intervals: Matomäki–Radziwiłł (y = x^ε, [x, x + C√x]).
  - All of these concern sets of size ≥ x^{1/2}.
  - Level of distribution: Fouvry–Tenenbaum 1996, Harper 2012 (BV/BDH), Drappeau 2015 (3/5).
  - EH limitations for primes: Friedlander–Granville 1989.
- **Expanders:** Jao–Miller–Venkatesan 2009 (the same prime-character input as Prop GX).
- **Derandomization:** Impagliazzo–Wigderson 1997 ⇒ deterministic L[1/2, O(1)] via Lenstra–Pomerance (polyimpl
  Thm DR; codex re-derived); not known to follow from GRH. Kabanets–Impagliazzo 2004.
- **Diagonal dominance:** Lévy–Desplanques (classical). **Wronskians of independent polynomials:** classical
  (Bôcher).

## 8. Experiments (pre-registered; toy sizes; correctness and mechanism only)

(Filled in below after the runs; see `experiments/analysis/zroute_summary.json`.)

## 9. Reviews

(See below.)

## 10. Candidate cards

```text
ID: zroute-1 deterministic Z route, degree one, dominance labels (Theorem DZ).
Target: every N, deterministic, ERH.
Proved: XK-Z not needed at d = 1; DD1 forces rank; Miller under ERH; complete factorisation by recursion.
Missing: Hypothesis DSR (deterministic smooth representation of q^H); existence free under GRH (Prop GX).
Barrier: Fourier/large-sieve templates need >= N^{1/2} (sup on S) or N/delta (Parseval) candidates, or
  spectral norm L^{O(1)} of the good set (Prop FC); averaging over fields/N needs level 1 - o(1) (Obs LV).
Status: proved reduction; supply OPEN.

ID: zroute-2 seedless deterministic NFS Z route.
Proved: XK-Z (deterministic, all degrees), Prop PIN.
Missing: pointwise NFS yield + PIN (no labels).  E1: factors whenever PIN reached.
Status: reduction proved; supply and rank heuristic.

ID: zroute-3 randomized Z-route Theorem A (Prop ZA).
Proved from A (ERH): character- and fruitfulness-free; hidden labels force rank (Lemma HR).
Status: diagnostic re-derivation; seeds' roles = supply + rank.
```

## 11. Reproduction

```
PYTHONPATH=. .venv/bin/python -m src.zroute.harness e3 experiments/configs/zroute.json experiments/raw/zroute_e3.jsonl
PYTHONPATH=. .venv/bin/python -m src.zroute.harness e2 experiments/configs/zroute.json experiments/raw/zroute_e2.jsonl
PYTHONPATH=. .venv/bin/python -m src.zroute.harness e1 experiments/configs/zroute.json experiments/raw/zroute_e1.jsonl
PYTHONPATH=. .venv/bin/python -m src.zroute.harness e4 experiments/configs/zroute.json experiments/raw/zroute_e4.jsonl
PYTHONPATH=. .venv/bin/python -m src.zroute.analyze > experiments/analysis/zroute_summary.json
PYTHONPATH=. .venv/bin/python -m src.zroute.posthoc_e3          # post hoc, labelled
.venv/bin/python -m pytest -q tests/test_zroute.py
```
