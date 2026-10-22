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

**No deterministic subexponential factoring algorithm under GRH was obtained.** For the direct degree-one design,
everything except relation supply is now proved, and the missing piece is one explicit constructive lemma.

1. **The algebra is done, deterministically [proved].**
   - Lemma XK-Z (§2) turns NFS-type rows into pure rational relations ∏ q^{e_q} ≡ 1 (mod N): fixed or growing degree,
     composite N, any degree parity.
   - It works by exact integer cancellation, archimedean logarithms, LLL and a unit gap, with no random prime.
   - It returns a spanning set of everything that exact ideal cancellation plus removal of the free unit part can
     give, with an exact rank formula. Cost: L^{O(a)} for B = L^a and polylog-height input data.
   - Under ERH, Miller's exponent-based splitting then factors N deterministically once the primes below 2 ln² N lie
     in the Q-span of the relations (Prop PIN). Full rank is not needed.
   - GPT-6.1-sol proved XK-Z independently in its blind round.
2. **For the degree-one design, no seeds, shape test, descent or recursion are needed for correctness [proved].**
   - Rows come from targets λ with known exponent labels: a Gauss split λ ≡ u/v with both halves B-smooth, tested by
     trial division.
   - The existing faster designs (α < 1/2) use number fields and descent, which retain the shape test, met by
     recursion.
   - In the randomized Z route, random seeds keep two jobs, supply and rank. Lemma HR: random hidden exponent labels
     make full Q-rank automatic.
3. **Rank can be forced deterministically in suitably labelled designs [proved].**
   - Lemma DD1 (§4): with one split row per prime q, target q^H·w and H above the row's total number of prime factors,
     the relation matrix is strictly diagonally dominant, hence nonsingular, whatever the splits are.
   - Codex's Wronskian-label Lemma WR forces rank after exact elimination too, given enough within-field torsion
     dependencies.
   - For these designs rank is not part of the missing lemma. Seedless designs, unlabelled seeded designs and thinly
     populated field blocks still need a rank hypothesis (PIN).
4. **The crux, exactly [open].** What is missing is a deterministic subexponential *constructor*:
   - **Hypothesis DSR (deterministic smooth representation):** given N and a target a = q^H (q ≤ B), find within
     L^{O(1)} time nonzero B-smooth u, v with |u|, |v| ≤ N³ and a small B-smooth multiplier w such that
     u ≡ a·w·v (mod N).
   - **Theorem DZ [proved; ERH]:** DSR at B = L[1/2, a] implies deterministic L[1/2, O(1)] complete factorisation of
     every N. DSR at B = L[2/5, a] (same integer outputs) implies L[2/5, O(1)].
   - **The outputs exist under GRH** (Prop GX: every invertible residue class mod N contains an L^a-smooth integer
     below N^{2+ε}). The existence proof supplies no subexponential constructor, and naive searches are exponential.
   - The orchestrator's hope "success is checkable, so a short explicit list suffices" fails for that reason, as in
     abcsunit (existence free, constructivity is the failure).
5. **Why the analysed certificates do not give DSR [one proved template statement + observations].**
   - **Prop FC (proved; norm templates).** A Hölder certificate (p ∈ {1, 2, ∞}, extended to all p by norm monotonicity)
     can show that an explicit list of m candidates meets a good set S of density δ only in these cases:
     - when m > φ(1 − δ) (Parseval only);
     - when m > (√3/4)(φ/δ)^{1/2} (sup on S);
     - when S has Fourier-algebra norm < δ(1 + √(2m)).
   - Cosets and subgroup complements have algebra norm ≤ 2, which is why Bach-type certificates for nonresidues need
     only polylogarithmically many candidates. We conjecture (no proof, no reliable spectral model) that smooth-split
     sets have large norm.
   - Directly applying existing smooth-number distribution theorems to short explicit families, or averaging over
     N, would need level of distribution near 1 with L^{O(1)} terms per class (Observation LV). Existing technology
     stops near 5/8 (Pascadi). This explains why FM/LD* certify only N^{Ω(1)}-size families. It is not a lower bound
     for all methods.
6. **Literature [§7].** The record is still N^{1/5+o(1)}: Harvey 2021 and Harvey–Hittmeir 2022 unconditionally;
   Shanks 1971 and Schoof 1982 under GRH.
   - Gao–Feng–Hu–Pan (IACR ePrint 2025/1004; to appear in Math. Comp.) improve log factors for balanced semiprimes.
   - Harvey–Hittmeir 2026 give deterministic large-order elements, so orbit designs can afford elements of order
     L^{O(1)}.
   - We found no deterministic subexponential algorithm under GRH for factoring, general class groups or prime-field
     discrete logarithms, and no hardness result.
   - The nearest GRH-derandomization success is Evdokimov's quasi-polynomial deterministic polynomial factoring over
     finite fields (1994).
7. **Experiments [§8; toy sizes, correctness and mechanism only; 4 pre-registered predictions failed, all
   diagnosed].**
   - **E1** (seedless deterministic cubic NFS + XK-Z): all 50 780 output relations verified exactly. PIN ⇒ Miller
     split N in 6/6.
   - PIN itself was reached in only 6/20 instances (P1c failed). The cause is the predicted mechanism: rare Q₀ primes
     never occur, and fixed boxes run out. Miller still split all 20 instances at some stage.
   - **E2:** deterministic orbit targets need 10–19n rows for full rank, against exactly n for random labels.
   - **E3:** the orbit yield equals the random yield on average but is over-dispersed (variance/mean 4.3; P3c
     failed). Post hoc, 96% of consecutive orbit hits are trivially derived duplicates.
   - **E4:** the dominance design (DD1) forced full rank and split N in 25/25 instances (40–56 bits).
   - Tries per prime are heavy-tailed (P4d failed: mean 30.7 against a window of 3.65–14.6 at 40 bits). Post hoc, the
     targets a·w inherit the target's own lattice structure, so pointwise cost per prescribed target is the visible
     weak point, exactly where DSR sits.

**Probability assessment (subjective, not a mathematical claim).** Proving DSR, or an equivalent constructive
statement, with known techniques: ≲ 5% (ours); codex ≈ 10% for deterministic L[1/2, O(1)] under GRH.

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
  So w_d·W ⊆ Λ ⊗ Q is spanned by valid relations.
  - For odd d, μ_K = {±1}. The cheaper o(x) := Σ_i x_i r_i + εe_{−1} (ε read off from ev(Σ x_i r_i) = ±1) is also
    valid, but spans W + Qe_{−1} only together with 2e_{−1}. It can drop a direction of W, e.g. a rational image
    −e_{−1} (codex). The code uses this variant and always appends 2e_{−1}.
  - "Everything exact cancellation can give" means W. It excludes modular relations coming from non-torsion units
    whose images happen to equal 1 mod N.
- (ii) A deterministic algorithm computes integer vectors whose Q-span is W. dim W is then obtained by one exact
  rank computation. Alternatively, process a basis of ker Aᵀ adapted to ker Aᵀ ∩ ker Rᵀ first, so that the remaining
  outputs are triangularly independent (codex).
- (iii) dim W = rank_R([V; T]) − rank_R(V), where Z is any integer basis of Λ₁ ⊗ Q, T := RᵀZ and V := (Log γ_z)_{z ∈ Z}
  (codex). In particular dim R^T(Λ₁ ⊗ Q) − r_U ≤ dim W ≤ dim Rᵀ(Λ₁ ⊗ Q), with r_U = r₁ + r₂ − 1 ≤ d − 1.
- (iv) **Cost.** Assume:
  - M, C_alg and the number of atom occurrences are ≤ L^{a+o(1)};
  - entries of A, r_i and y_{iP} are ≤ (log N)^{O(1)};
  - d ≤ log N;
  - f and the atom coordinates have (log N)^{O(1)} bits. Then the deterministic
  cost is L^{(ω+1)a+o(1)} with fast matrix multiplication, or L^{4a+o(1)} classically, plus L^{3a+o(1)} for the unit
  step. Outputs have L^{a+o(1)} bits per coordinate (L^{3a+o(1)} in total). This holds at fixed d and at d ≍ (log N/log log N)^{1/5} (L = L_N[2/5, 1]).

*Proof.*
- **(i) Correctness.**
  - For x ∈ Λ₁, γ_x ∈ O_K^×. This is dlog Lemma XK1: the 𝔞-column gives Σ_P y_P = 0, the ideal columns cancel
    every 𝔯 | 𝔍_P, and N(𝔍_P) = |F(P)| is B-smooth, so (γ_x) = O_K. Its proof is degree-free (dlog25 §4).
  - Bad primes r | t·disc f: XK0(iii)/XK0-U computes r-maximal orders and valuations (Buchmann–Lenstra 1994).
    - The only randomness there is factoring f mod r. Replace it by Berlekamp's deterministic algorithm, poly(d)·r
      operations with r ≤ B. There are O(log N) bad primes, so this costs L^{a+o(1)}.
  - Multiply (2.1) over i with exponents x_i. The auxiliary columns are columns of A, so they cancel. Lemma Ψ then
    gives ev(Σ x_i r_i) ≡ ψ(γ_x) (mod N).
  - For x ∈ Λ₂, γ_x = ζ ∈ μ_K ∩ 𝒜, so ψ(γ_x)^{|μ_K|} = 1.
  - If ζ ∈ K has order w then φ(w) ≤ d, and φ(w) ≥ √(w/2) gives w ≤ 2d². So w_d kills μ_K. For odd d, μ_K = {±1}.
- **(ii) Deterministic Q-basis of Λ₁.**
  - Hadamard: every minor of A has |Δ| ≤ Had := (√C_alg·‖A‖_∞)^{C_alg}.
  - Take the first primes p₁ < p₂ < … above 2^{62} until ∏p_k > 2·Had² (enough both for rank detection, which needs
    ∏ > Had, and for signed CRT after skipping primes dividing det A_{IJ}). Then rank_Q A = max_k rank_{p_k} A, because a
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
    - It is certified by ‖Σ n Log γ‖_∞ < β_d/2 in ball arithmetic, with κ chosen so that
      2^κ β_d > 8(1 + s√d)R′ and R′ = 2^{(s−1)/2}√s(1 + s√d)R.
    - β_d is the Blanksby–Montgomery gap. Codex's elementary gap δ_d = 1/(d T_d), T_d = d(6·2^d + 1)^d + 1, also works,
      but its log(1/δ_d) = O(d²) adds O(d³) to log R. That does not change the L exponent.
    - Output z^{(j)} := Σ n_i x^{(i)} + n_j x^{(j)} ∈ Λ₂.
  - **Case 2: no relation exists.** Then every LLL vector has ‖Σ n Log γ‖_∞ ≥ β_d − (rounding) > β_d/2, the
    certificate fails, and j joins I′. (It is enough that the certificate cannot pass; dependence would have
    guaranteed passage.)
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
  - "Miller" here is the exponent-based splitting procedure with the known multiple D of λ(N), not the strong
    probable-prime test with exponent N − 1.
  - Projection rank is not enough (codex): the other coordinates must cancel. Codex also notes that no single
    numerical rank condition is both necessary and sufficient.
  - In E1 and E2, Miller often split N before PIN held (E1: see §8), because the pinned subgroup already contained a
    witness. This is luck, not a guarantee.

## 3. Randomness sources revisited: what the Z route still needs

| detGRH row | randomness | status in the Z route |
|---|---|---|
| 1–3 | field law (g₀, g₁, w, c) | the degree-one design (§3.2) uses no number field, and DSR (§5) is stated at any α. The existing descent designs for α < 1/2 use fields, where the field law is supply-only |
| 4a | seed law, supply role | **the crux**: replaced by Hypothesis DSR (§5) |
| 4b | seed law, fruitfulness role | gone (Miller under ERH) |
| — | seed law, **rank role** | random labels give rank automatically (Lemma HR); deterministic labels force it (Lemmas DD1, WR, §4) |
| 6 | shape test | not needed in the degree-one design (halves are tested for B-smoothness by trial division). Descent designs need complete factorisations of the halves, by recursion. That needs a quantified recurrence (T(N) ≤ L[α, A] + L[α, p]T(N^θ), c ≥ max(A, p/(1 − θ^α))), not only determinism |
| 7, 7′ | kernel sampling, exact LA | deterministic (Lemma XK-Z (ii)) |
| 8, 9 | characters, Hensel root | gone (no square roots) |
| 10–13 | — | deterministic |

### 3.1 Random seeds: exactly two jobs remain

**Lemma HR (hidden labels force rank) [proved; the random-exponent rank argument is standard in rigorous index
calculus (Pomerance 1987; Hafner–McCurley 1989; Stange 2022). New here, as far as we found: only the conditional form
given the transcript, applied to post-elimination combinations].**
- *Setting.* Let s₁, …, s_M be i.i.d. uniform on [0, X)^n with X ≥ N⁶, and T_i := ev(s_i). Let 𝒯 be a random
  variable that includes (T_i)_i, such that (s_i)_i is conditionally independent of 𝒯 given (T_i)_i. (For example, 𝒯
  is all targets, descents, field data and the coins of every other step.)
- Let e_i ∈ Z^n be 𝒯-measurable, put r_i := s_i − e_i, and let z^{(1)}, …, z^{(K)} ∈ Z^M (K ≥ n) be 𝒯-measurable and
  linearly independent, with K ≤ N^{1/3} and entries of at most N^{1/3} bits. (The application has K, bits ≤ L^{O(1)}.)
- *Claim.* Pr[rank_Q{Σ_i z^{(k)}_i r_i}_k < n | 𝒯] ≤ 2/N + 4MnN³/X.
- *Proof.*
  - **Choose ℓ′.** Fix (given 𝒯) a prime ℓ′ ∈ (N, 2N] not dividing a nonzero K × K minor Δ of the z-matrix. One
    exists, since that interval has ≥ N/(2 ln N) primes, while Δ has ≤ K(N^{1/3} + log K) ≤ 2N^{2/3} bits and so
    fewer than N/(2 ln N) prime divisors above N.
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
  - Codex zroute_1 §3.3 / zroute_2 #8 gives a stronger proof with no bound on K or heights. Normalise n pivot
    columns, condition on everything else, and note that each remaining label coordinate is uniform on an arithmetic
    progression of step ≤ N with ≥ X/N − 1 terms. This gives Pr[failure | 𝒯] ≤ nN/(X − N).
  - The labels must stay in the rational block. Treating them as separate auxiliary columns cancels the signal.
  - HR does not give GF(2) rank: conditioning on T can fix parities.

**Proposition ZA (a character-free, fruitfulness-free Theorem A) [proof sketch; proved from A; ERH].**
- *Statement.* In Theorem A, replace:
  - R = r² by T = ev(s), s uniform on [0, N⁶)^n;
  - Steps 8–9 (LV6, characters, Hensel root, Theorem DD) by Lemma XK-Z, Lemma HR and Miller with deterministic bases.
  The result is a Las Vegas algorithm under ERH with expected time L[1/2, min_a max(a + 1/√2, a + θ/(2a), (ω + 1)a)
  + o(1)].
  - To keep dlog's constant (3a in place of (ω + 1)a), keep only the Q₀ coordinates as hidden outputs, force the
    other rational coordinates to cancel, and compute polynomially many exact dependencies with the few-right-hand-side
    machinery of dlog XK, then use PIN (codex zroute_2 #9). We have not written this out.
- *Proof sketch (the only new step is supply).*
  - Under ERH, G_B = U. T is then within nN/N⁶ of uniform on U (codex zroute_1 §2.5).
  - dlog Lemma S_p holds at composite N with the same proof:
    - the coprime-pair map (u, v) ↦ uv^{−1} is injective because |uv′ − u′v| < Y² < N;
    - the halves are units because the accepted prime set G is gcd-checked against N (a failure is a factor, and
      dropping ≤ ω(N) primes costs nothing);
    - the Gauss-reduced vector is the split.
    So Theorem K(ii) holds without the 2^{−ω} factor.
  - Acceptance and the descent depend on s only through T, so Lemma HR applies. It needs a torsion kernel of
    dimension ≥ n, hence M ≥ C_col + d + n rows (the extra n). XK-Z costs L^{(ω+1)a} (deterministic).
  - Seed evaluation costs L^{a+o(1)} per trial, which is dominated.
  - The same substitution in Theorem U (L[2/5]) is "proved from U" with dlog25's XK-U.
- *Status.* This is a re-derivation, not an improvement. Its point is diagnostic: in the Z route, random seeds do
  supply and rank, nothing else.

### 3.2 Seedless and field-free designs

- **Degree-one rows [proved].** Take λ ∈ U with known exponent vector x(λ) and the Gauss-reduced (u, v) with
  u ≡ λv. If u and v are B-smooth, then x(λ) − e(u) + e(v) ∈ Λ.
  - The smoothness test is trial division (or a remainder tree) by FB, deterministic, L^{a+o(1)}.
  - There is no shape test, no descent and no recursion. This is the classical Dixon-type row, read over Z.
  - With ordinary labels (label = exponent vector of the integer λ), small targets are useless: the relation is an
    integer identity before reduction (detGRH Lemma SR; codex: λ = 2 with |u|, |v| < N/3 forces u = 2v). A small
    residue carrying a large known label, such as q^H mod N, is not covered by this, and DD1 uses exactly that.
- **Seedless two-sided NFS [proved].** A pair (a, b) with a − bm and F(a, b) B-smooth is a row of type (2.1). No
  target, seed, shape test or square root is needed.
  - Rank then has no labels to lean on, and PIN is a genuine hypothesis (E1 shows it failing for lack of
    occurrences; §8).
- **Pooling across fields [proved; codex].** Rational outputs of different fields add (W_total = Σ W_j), but each
  field must be eliminated with its own rows. dim W_j ≤ M_j − rank(A_j) − rank(V_j), with V_j the Log image of field
  j's kernel. rank(V_j) ≤ r_U, so the unit rank bounds the loss from above, not dim W_j. A union-over-fields
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
  extension to split rows with multipliers w is ours ("no prior formulation found by us").
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
- *Hypotheses made explicit (critic).* Every accepted row carries its own label, the indices j are pairwise distinct,
  H₀ bounds the actual kernel vectors produced by INC (≈ Had·R), and for odd d the εe_{−1} correction is part of V.
- *Notes.*
  - Labels have L^{O(1)} bits when P, n, H₀, B₀ are L^{O(1)}-bounded.
  - It applies after XK-Z, to descent rows, where dominance is lost through the combinations z.
  - Test `test_wronskian_labels_force_rank` checks 300 random instances.

**Consequences for the rank question.**
- For the DD1 design (one split per FB prime) and for WR-labelled designs with uniform certificate bounds and a torsion
  kernel of dimension ≥ n, rank is forced deterministically. It is not part of the open problem there.
  - This does **not** cover arbitrary seeded designs, orbit targets or union-over-fields statements. n accepted
    descent rows need not give n torsion dependencies, and a field whose algebraic map is injective contributes
    nothing (codex).
- In seedless designs (two-sided NFS, intervals near kN) PIN remains a hypothesis. E1 and E2 show its mechanism at toy
  scale: a prime q is pinned only once it occurs, which takes ≈ q·log(#small primes) rows.
- **Relation to rkrank.** The notions are incomparable as stated.
  - RK concerns parity cancellation of raw rows (square certificates). The Z route needs exact cancellation and unit
    removal.
  - A row with rational entry 1 and algebraic entry 2 has a rational direction after parity cancellation and no exact
    kernel at all (codex). So RK does not imply full Q-rank of the Z-route outputs. (An earlier draft said it did;
    withdrawn.)
  - Conversely, full Q-rank can hold with zero F₂ information: the doubled lattice (detGRH E1; codex §3.2).
  - For a single integer lattice, rank_Q L ≥ rank_{F₂}(L mod 2) holds, but that is not RK's object.
  - With DD1/WR labels the Z route needs no rank hypothesis.

## 5. Deterministic relation supply: the crux

### 5.1 The conditional theorem and the missing lemma

**Hypothesis DSR_α(a, c) (deterministic smooth representation) [open].**
- *Statement.* There is a deterministic algorithm which, given N and a target a ∈ {q^H mod N : q ≤ B}, runs in time
  L_N[α, c] and outputs nonzero integers u, v, w with their factorisations, such that:
  - all prime factors are ≤ B, where B = L_N[α, a] and H := 8⌈log₂ N⌉;
  - |u|, |v| ≤ N³ and Ω(|w|) ≤ log₂ N;
  - u ≡ a·w·v (mod N).
  It must hold **uniformly for every** preprocessed N (odd, not a prime power, no prime factor ≤ max(B, 3 ln² N)).
  An almost-all version would not give Theorem DZ.
- DSR := DSR_{1/2}.
- *Variants.*
  - **Smooth compression** (codex zroute_1 §5): v = w = 1, i.e. an integer representative u ≤ N³ of q^H. This is a
    special case of DSR's output, so constructing it suffices.
  - **DS-T:** a deterministic list of P = L^{c} targets with Wronskian labels (Lemma WR), with explicit bit
    complexity, such that at least n of them have Gauss splits with both halves B-smooth. Degree one, so XK-Z is not
    needed and WR applies with V = the halves' exponents.
  - **Descent version [not formulated; conditional sketch only].** Replacing u, v by NFS descent certificates would
    need a full specification: the constructor, the row distribution across fields, enough within-field torsion
    dependencies for WR, certificate bounds, and a quantified recursion budget for the shape test. We do not state
    or claim it. (An earlier draft's "DSR-U ⇒ L[2/5]" for a descent version is withdrawn; both reviews.)

**Theorem DZ [proved; ERH; the composition is folklore-level].** For fixed a, c > 0 and α ∈ {1/2, 2/5} (any fixed
α ∈ (0, 1) works the same way):
ERH + DSR_α(a, c) ⇒ a deterministic algorithm that completely factors every integer N in time
L_N[α, max(a + c, (ω + 1)a) + o(1)].
- *Proof.*
  - **Preprocessing,** reapplied to every cofactor in the recursion: remove prime factors ≤ max(B, 3 ln² N) with
    multiplicity; test primality (AKS); detect perfect powers and recurse on their bases. Then FB ⊇ Q₀, and every
    integer below 2 ln² N is a unit.
  - **Rows.** One DSR call per prime q ≤ B, π(B) = L^{a+o(1)} calls in all, costing L^{a+c+o(1)}. Each row is
    He_q + e(w) − e(u) + e(v) (with signs in the −1 column). Since Ω(|u|) + Ω(|v|) + Ω(|w|) ≤ 3 log₂ N + 3 log₂ N
    + log₂ N < H, Lemma DD1 gives full rank n together with 2e_{−1}.
  - **Splitting.** D = |det| (entries O(log N), L^{a+o(1)} bits) is computed deterministically in L^{(ω+1)a+o(1)}.
    Under ERH, Q₀ generates U, so λ(N) | D. Miller's exponent-based splitting with bases a < 2 ln² N then splits N
    (detGRH Prop Z). This also holds when N is not squarefree (e.g. p²q with p > B), because N has ≥ 2 prime-power
    CRT components.
  - **Recursion.** There are O(log N) splitting nodes, each costing at most the bound for N. Factors need not be
    coprime; recursing on both and keeping multiplicities is correct.
- At α = 2/5 the same proof applies verbatim. DSR_{2/5} asks for the same integer outputs with B = L[2/5, a], and
  GX gives their existence at that scale too. No number field, descent or shape test is used. ∎
- *Remark.* The degree-one implication needs neither XK-Z nor a separate rank hypothesis. XK-Z is what any number-field
  formulation of supply would need.

**Proposition GX (the outputs exist under GRH) [proved; folklore: a corollary of Jao–Miller–Venkatesan's GRH
expander; GRH for Dirichlet L-functions mod N].**
- *Statement.* Let y = L_N[α, a] with α ∈ (0, 1) and a > 0 fixed, and ε ∈ (0, 1). For N large, every **invertible**
  residue class mod N contains a y-smooth integer u ≤ N^{2+ε}. In particular DSR_α's output exists, with v = w = 1.
  - The class 0 contains no nonzero B-smooth integer after preprocessing.
- *Proof (codex zroute_1 §4.5, re-derived).*
  - Let μ be uniform on the primes ≤ y not dividing N. Under GRH, |Σ_{p≤y} χ(p)| ≪ √y log²(Ny) for χ ≠ χ₀ mod N (from
    ψ(y, χ) ≪ √y log²(Ny); prime powers contribute O(√y)). So max_{χ≠χ₀}|μ̂(χ)| ≤ y^{−1/2}(log Ny)^{O(1)}
    = y^{−1/2+o(1)}, since log N = y^{o(1)}.
  - μ^{*k}(a) = φ(N)^{−1}Σ_χ χ̄(a)μ̂(χ)^k ≥ 1/φ(N) − max_{χ≠χ₀}|μ̂(χ)|^k > 0 for k := ⌈(2 + ε/2)·ln N/ln y⌉.
  - The product of k primes ≤ y is ≤ y^k ≤ N^{2+ε} for N large. ∎
- *A pair-counting variant (a routine character-sum count; α = 1/2 only).*
  - Take B = L^a with a > 1 and y = L^{1/a+δ} with 0 < δ < a − 1/a (so y ≤ B). Let S be the B-smooth integers in
    [1, N/2], all units, s := |S| = N·L^{−1/(2a)+o(1)}.
  - Weighted by Λ(n), #{(n, u, v) : n ≤ y, u, v ∈ S, u ≡ nv (mod N)} = ψ(y)s²/φ(N) + E, with
    |E| ≤ max_{χ≠χ₀}|ψ(y, χ)|·(s − s²/φ(N)).
  - The ratio error/main is L^{−δ/2+o(1)}, so the count is ≥ (1/2)ψ(y)s²/φ(N) = N·L^{δ+o(1)}.
  - So the factor base is an ε-biased set at the Parseval scale (bias y^{−1/2+o(1)} against the floor m^{−1/2}),
    up to logarithmic factors. This gives many relations but no prescribed-target representations.
- **Constructivity is the gap.** The existence proofs supply no subexponential constructor. Straightforward
  implementations (enumerating k-fold prime products, about N^{2+o(1)} of them, or searching one progression below
  N^{2+o(1)}, about N^{1+o(1)} candidates) are exponential. This is a statement about naive implementations, not a
  lower bound.

**On how hard DSR is (an analogy, not a claim).**
- A strengthened version of DSR has a prime modulus p and arbitrary prescribed targets. It would supply the
  relation-collection and individual-descent steps of a deterministic index calculus in F_p^*.
- We found no deterministic subexponential prime-field discrete-logarithm algorithm under GRH (§7).
- But DSR as stated assumes composite N and only the targets q^H. A full-rank relation sublattice also does not by
  itself determine logarithms: saturation and the large prime-order components of p − 1 remain.
- **We know no reduction in either direction** between DSR and prime-field discrete logarithms.
- Prop EQ (detGRH) shows deterministic factoring needs only *arbitrary* full-rank relations. DSR asks for more
  structured ones.

### 5.2 Why known certificates do not give DSR

The orchestrator's question is whether an explicit list of L^{O(1)} fields, seeds or targets can be certified to
contain good ones, and whether detGRH's "N^{Ω(1)} families" is a proof artifact.

**Proposition FC (norm-certificate templates) [proved; a statement about certificates, not algorithms. Parts (a)
and (c) are standard small-bias facts (Naor–Naor; AGHP; Azar–Motwani–Naor), specialised here. No prior formulation of
the comparison with smooth-split sets found by us].**
- *Setup.*
  - G = U, φ = |G|.
  - W ⊂ G is the explicit candidates, |W| = m ≤ φ/2.
  - S ⊂ G is a good set, |S| = δφ, δ ≤ 1/4. The good set must be specified: for example, the residues whose
    canonical Gauss split has both halves B-smooth. It is not "residues with some smooth representative below N³",
    which is all of U by GX.
  - #(W ∩ S) = mδ + E, with E = φ^{−1}Σ_{χ≠χ₀} Ŵ(χ)\overline{Ŝ(χ)}.
  - A Hölder norm certificate proves #(W ∩ S) > 0 by showing mδ > ‖Ŵ‖′_p‖Ŝ‖′_{p′}, where primes mean χ₀ omitted and
    normalisation by φ^{−1}.
- (a) **(2, 2) (Parseval only) succeeds iff m > φ(1 − δ)** (exact, codex). So it never succeeds for m ≤ φ/2.
  - *Proof:* the condition is mδ > (m(1 − m/φ)·δφ(1 − δ))^{1/2}.
- (b) **(1, ∞) (sup on S) fails unless m > (√3/4)(φ/δ)^{1/2}, for every S.**
  - *Proof:* φ^{−1}Σ_{χ≠χ₀}|Ŵ| ≥ (m/2)/m = 1/2, and max_{χ≠χ₀}|Ŝ| ≥ (3δφ/4)^{1/2} by Parseval.
- (c) **(∞, 1) (sup on W) succeeds only if ‖1_S‖_A := φ^{−1}Σ_χ|Ŝ(χ)| < δ(1 + (2m)^{1/2}).**
  - *Proof:* max_{χ≠χ₀}|Ŵ| ≥ (m/2)^{1/2} by Parseval.
  - Random S has ‖1_S‖_A ≍ (δφ)^{1/2}. The critic measured 0.77–0.89 × (δφ)^{1/2} in (Z/4099)^*.
  - Cosets have ‖1_S‖_A = 1 and subgroup complements have 2(1 − 1/k). This is exactly why ERH certificates for
    nonresidues or subgroup escape (Ankeny, Bach) need only polylogarithmically many candidates. So "every
    character-sum certificate needs N^{Ω(1)} samples" is false (codex).
  - **Conjecture (speculative; no proof and no reliable spectral model):** canonical smooth-split sets do not have
    ‖1_S‖_A = L^{O(1)}.
    - An earlier draft argued this from twisted smooth sums |Ŝ(χ)| ≈ |S|e^{−O(u)} "for many χ". That argument is
      withdrawn. By Parseval at most e^{2Cu}/δ characters can be that large, and twisted smooth-number sums do not
      describe the indicator of a split set.
- (d) **Other exponents and weights.** φ^{−1}Σ_{χ≠χ₀} has total mass < 1, so ‖·‖′_q ≥ ‖·‖′_p for q ≥ p.
  - Every p ≤ 2 therefore inherits (b)'s necessary condition up to constants, and every p ≥ 2 inherits (c)'s.
  - Weighting W does not help, since a probability vector w on m points has ‖ŵ‖₂² ≥ 1/m − 1/φ.
- (e) **Pair certificates (a corrected template calculation).**
  - Let S be the smooth integers in [1, Y], Y < N, s := |S|. The hits are pairs (u, v) with u ≡ wv, w ∈ W.
  - The (∞, 1) certificate applied to |Ŝ|² needs max_{χ≠χ₀}|Ŵ| < ms/(φ − s), hence m > (φ − s)²/(2s²).
  - If the certified family is then searched by exhaustive lattice enumeration of {u ≡ wv} ∩ [−Y, Y]² (≈ 4Y²/N
    points for generic w; area/determinant is **not** a lower bound for every w, cf. codex's N = 1009·1013, w = 1000),
    the work for generic w is Ω(N/ρ_Y²).
  - This is a limitation of this certificate combined with this search. It is not a cost lower bound for pair-based
    algorithms. (Corrected from an earlier, unqualified "≥ 2N/ρ²".) ∎
- *Not excluded.* Certificates that use the joint structure of Ŵ and Ŝ rather than norms, and non-Fourier arguments.

**Observation LV (why direct applications of distribution theorems need level ≈ 1) [observation].**
- **Fields.**
  - For a fixed atom P, FM's field average realises the progression F̂(P) + G(P)·j, with modulus t = N^{μ+ν} and
    length X_F/t. The moments use level ≤ 1/2 in t (LV Lemma 8.5′, Harper, Drappeau).
  - An explicit consecutive family of L^{O(1)} fields realises a progression of length J = L^{O(1)}, with level
    log t/log(tJ) → 1.
  - So a *direct* application of long-progression smooth-number estimates fails. Individual-progression theorems need
    q ≤ y^{4√e−δ} (Harper 2012, sufficient condition). Averaged results with absolute values over moduli and a fixed
    residue reach 3/5 (Drappeau 2015), and 66/107 and 5/8 (Pascadi 2025; arXiv 2505.00653; verified by the literature
    critic), not ≈ 1. Their constants depend on the residue.
  - This does not show that every explicit family needs a level-1 theorem: families can exploit structure that the
    progression model misses.
- **N ("almost all N").** Averaging a fixed explicit design over N ∈ [X, 2X] leads, in the natural formulation, to
  equidistribution of the bilinear sequence {u − λv} in the class 0 modulo N with only L^{O(1)} expected elements per
  modulus.
  - Our attempted dual-large-sieve bound was **vacuous**: injectivity forces the range Z > 2|𝒜|, so the bound exceeds
    2X (codex; critic). It also needs care with repeated reduced fractions at composite N.
  - We record it as a failed template. Nothing here shows that almost-all-N methods need level 1 − o(1).
- **Answer to the orchestrator's question.**
  - "Current proofs certify only N^{Ω(1)}-size families" holds for FM/LD*. It also holds for the norm-certificate
    templates of Prop FC, unless the good set has small Fourier-algebra norm.
  - Direct use of existing distribution theorems fails for short families.
  - We have no argument that it is forced for all methods. In that sense it is a limitation of the known templates, of
    unknown depth.
  - At growing degree, codex notes that some families are N^{o(1)} yet still ≫ L[2/5, O(1)].

**Observation (checkability) [proved, trivial; detGRH §6].** Global moments are compatible with total failure on any
prescribed list. Any proof must use the arithmetic of the list.

### 5.3 Remarks

- **Elements of large order.** Orbit designs (λ_j = g^j) need ord_N(g) > (number of targets).
  - Harvey–Hittmeir 2026 (arXiv 2601.11131v2, Thm 1.1; N ≥ 3, 1 ≤ D < N − 1, which is harmless for D = L^{O(1)})
    deterministically return an element of order > D or a proper factor in O(D^{1/2} log D log N/√(log log D)) bit
    operations. For D = L[α, C] that is
    L[α, C/2 + o(1)], so order acquisition is affordable.
  - Multi-base boxes avoid the issue in any case: a collision ∏q^{j} ≡ ∏q^{j′} of two targets is itself a relation,
    though one collision gives neither PIN nor spanning.
  - E3 found no orbit cycle within 2000 targets for 1000 semiprimes.
- **Redundancy of multiplicatively propagated representations [proved; trivial].** If λ ≡ u/v and the labels satisfy
  x(qλ) = x(λ) + e_q (true for the orbit and for ordinary labels), then the derived representation qλ ≡ (qu)/v gives
  the identical row.
  - Post hoc for E3: 96% of consecutive orbit hits were such duplicates.
  - This favours label schemes where x(qλ) ≠ x(λ) + e_q (WR) or one row per prime (DD1). It does not prove that every
    multiplicatively closed family is rank-deficient.
- **Deterministic heuristic L[1/3] (remark; folklore; heuristic input).** NFS is deterministic apart from the
  linear algebra and the square root.
  - Seedless base-m NFS with growing degree is deterministic. Deterministic enumeration, polynomial exact
    post-processing (XK-Z), PIN and Miller preserve α = 1/3.
  - The unproved inputs are a pointwise simultaneous-smoothness/yield heuristic for the explicit polynomial and box,
    and PIN. No fruitfulness heuristic is used.
  - For WR-labelled targets the usual NFS yield heuristic is not enough: one must also assume prescribed-target
    descents in sufficiently populated field blocks.
  - The dense exact elimination worsens the constant.

## 6. Answers to the task items

| item | answer | status |
|---|---|---|
| 1 XK-Z | Lemma XK-Z: deterministic, composite N, all degrees. Gives a spanning set of W (exact ideal cancellation + removal of the free unit part); rank formula; cost L^{O(a)} under the stated input-height bounds; Miller under PIN | proved (codex independently) |
| 2a fields | the degree-one design needs none, and DSR is stated at any α. FM/LD* and the norm templates certify only N^{Ω(1)}-size lists (unless the good set has small algebra norm); direct use of known distribution theorems needs level ≈ 1. Not a lower bound | Prop FC (proved template statement) + observations |
| 2b seeds | random root seeds are unnecessary (degree-one or seedless rows); random seeds do supply + rank only. With ordinary labels, targets must wrap around (Lemma SR); labels should avoid trivial propagation | proved |
| 2c shape test | unnecessary in the degree-one design. Descent designs need a quantified recursion | proved / conditional |
| 3 rank | PIN suffices; random labels force rank (HR); DD1 and bounded WR designs force it deterministically. Seedless designs, unlabelled seeded designs and thinly populated field blocks still need PIN. RK and Z-route rank are incomparable | proved |
| 4 literature | record N^{1/5}; no deterministic subexponential factoring, class group or F_p DL under GRH found; no hardness known | search result |
| 5 experiments | §8 | experimentally supported (toy) |
| crux | Hypothesis DSR (uniform deterministic smooth representation of q^H). Outputs exist under GRH (Prop GX). Theorem DZ: ERH + DSR_α ⇒ deterministic L[α, O(1)], α ∈ {1/2, 2/5} | open / proved implication |

## 7. Literature

Verified by me where stated. (c) = theorem statement checked by codex zroute_1/zroute_2. (c*) = reported by codex,
not verified by me.
- **Deterministic factoring.**
  - Pollard 1974 / Strassen 1976: N^{1/4+o(1)}.
  - Under GRH: Shanks 1971 and Schoof 1982 §6, N^{1/5+ε} (c).
  - Unconditional:
    - Hittmeir 2018: N^{2/9+o(1)}.
    - Harvey, Math. Comp. 90 (2021) 2937–2950: O(N^{1/5}log^{16/5}N) (c).
    - Harvey–Hittmeir, Math. Comp. 91 (2022) 1367–1379: a further (log log N)^{3/5} saving (c).
  - Gao–Feng–Hu–Pan, IACR ePrint 2025/1004 / arXiv 2512.19076, to appear in Math. Comp. (doi 10.1090/mcom/4188): for
    balanced semiprimes, O(N^{1/5}log^{13/5}N/(log log N)^{3/5}),
    the same exponent. They mention N^{1/6}, N^{1/8} as possibilities. The paper also treats power-divisor and
    special-form problems. (Read by me: abstract and Thm 1.1.)
  - Harvey–Hittmeir, "Deterministic methods for finding elements of large multiplicative order", arXiv 2601.11131,
    ANTS XVII (2026). It removes the hypothesis D ≥ N^{1/6} from the large-order-or-factor subroutine and gives no new
    factoring exponent. (Read by me: abstract.) Its running time is quoted in §5.3 (verified by the literature
    critic).
    - Intermediate steps: Oznovich–Volk, arXiv 2506.07668 (SODA 2026), D ≥ N^{1/6}; Gao et al., D ≥ N^{1/4+o(1)}.
  - Carella, arXiv 1308.2891, claims N^{1/6+ε}. It is unvetted and unpublished, and not counted here.
  - Umans–Wang, arXiv 2511.10851: a combinatorial conjecture (sets S, T of size n^β with entries ≤ exp(n^α) such
    that every i ≤ n divides some s − t) implying deterministic N^{1/6} (verified by the literature critic). It is not
    subexponential and does not follow from GRH.
  - **The record N^{1/5+o(1)} is confirmed.** We found no deterministic L[α, O(1)] factoring for any α < 1 under GRH
    alone, and no non-generic lower bound excluding one. Generic-ring lower bounds (Shoup 1997; Aggarwal–Maurer
    2009) do not apply to index calculus.
  - Classical reductions exist: factoring versus totient or order-multiple information under ERH (Miller), and the
    unrestricted relation-lattice equivalence (detGRH Prop EQ). Further oracle equivalences: Bach–Miller–Shallit 1986
  (σ), Shallit–Shamir 1985, Woll 1987, and Morain–Renault–Smith (AAECC 34, 2023; φ, λ or order oracles). See also
  Adleman–McCurley's open problems (ANTS-I 1994).
- **Factoring from exponents (ERH).** Miller, JCSS 13 (1976) 300–317; Bach 1984 §3 and Math. Comp. 55 (1990);
  Shoup §§10.4–10.5. Żrałek, Math. Comp. 79 (2010): with O(ln n) queries to a φ-oracle, complete factorisation in
  deterministic time exp((1 + o(1))(ln n)^{1/3}(ln ln n)^{2/3}) (an oracle result).
- **Z route, prior formulations.** Stange, arXiv 2211.06821 (Hypothesis 3.1, Thm 3.2). Under ERH and Miller,
  factoring reduces to random relations among residues g^x with random x ∈ [1, n], analysed through coordinate
  covolumes; the runtime is heuristic. This is direct prior art for the hidden-label set-up (Lemma HR) and for PIN's
  per-coordinate idea.
- **Closest classical analogues of a DSR constructor.**
  - CFRAC (Morrison–Brillhart 1975): a deterministic sequence of residues Q_k < 2√N, heuristic L[1/2, √2]. Its rigour
    gap (smoothness along a deterministic sequence) is exactly DSR's gap (Pomerance 1982; Schnorr 1982).
  - Vallée, Math. Comp. 56 (1991): provable generation of small modular squares (randomized).
  - Boneh, STOC 2000 / JCSS 64 (2002): a deterministic CRT-decoding constructor of integers with a large smooth part in
    short intervals.
  - dBPW's sampling of smooth ideals in prescribed classes is the randomized counterpart of DSR.
- **Exact multiplicative relations.** Deterministic polynomial-time computation of multiplicative relations among
  explicit algebraic numbers: Ge 1993; Lenstra–Silverberg, Found. Comput. Math. 18 (2018), Thm 1.11 (also
  arXiv 1509.08843). This is prior art for the extraction in XK-Z, though not with the compressed-product cost ledger.
- **Class groups.**
  - Hafner–McCurley, JAMS 2 (1989): Las Vegas, ERH (imaginary quadratic fields).
  - Biasse–Erukulangara 2023: randomized, GRH (c).
  - de Boer–Pellet-Mary–Wesolowski, arXiv 2512.01588: rigorous under ERH, probabilistic (c).
  - Buchmann 1990 and Biasse–Fieker: heuristic.
  - Deterministic: Lenstra 1982 / Schoof 1982 |D|^{1/5+ε} for class number and regulator (class-group structure is
    |D|^{1/4+ε}); Schoof 2008, Arakelov class groups, deterministic exponential (c).
  - No deterministic subexponential class-group algorithm under GRH found.
- **Discrete logs in F_p^*.** Pomerance 1987: randomized L[1/2, √2]. No deterministic subexponential algorithm under
  GRH found.
- **Deterministic algorithms under GRH that do exist.**
  - Square roots mod p: Ankeny 1952 + Tonelli–Shanks, with Bach's constants.
  - Miller's primality test, superseded unconditionally by AKS (Annals 160 (2004)).
  - Polynomial factoring over F_q: Rónyai 1988/1992; Huang (abelian Galois groups); **Evdokimov 1994 (ANTS-I,
    LNCS 877): deterministic (n^{log n} log q)^{O(1)} under GRH**; later Ivanyos–Karpinski–Saxena. This is the closest
    GRH-derandomization success.
- **Smooth numbers.**
  - Single progressions: Granville 1993; Soundararajan 2008; Harper, JNT 132 (2012), q ≤ y^{4√e−δ},
    log x/log q → ∞ (c). These are sufficient conditions.
  - Short intervals under RH:
    - Soundararajan, arXiv 1009.1591: intervals ≍ u√x/ρ(u/2), y ≥ exp(5√(log x log log x)) (c).
    - Younis, arXiv 2409.05761, Thm 1.3: every fixed θ > 1/2 with y ≥ (log x)^{K(θ)}, which includes y = L^{O(1)},
      with interval length x^θ (c).
  - Unconditional:
    - Matomäki–Radziwiłł: [x, x + C√x] for y = x^ε; also almost-all intervals much shorter than √x for fixed
      smoothness parameter.
    - Jain, arXiv 2502.10530, Quart. J. Math. (2026): almost-all intervals of length
      exp((1 + ε)((11/8)u log u + 4 log log X)) for y ≥ exp((log X)^{2/3+ε}), improving Matomäki. Below √X, but the
      y-range misses L[1/2] and L[2/5] (verified by the literature critic).
  - None gives a pointwise smooth number in every explicit set of size L^{O(1)} near multiples of N.
  - Level of distribution:
    - Fouvry–Tenenbaum 1996; Harper 2012 (BV/BDH);
    - Drappeau 2015: 3/5, with absolute values over moduli and a fixed residue;
    - Pascadi, Compositio 161 (2025) 1923–1974, Thm 1.1: level 66/107; arXiv 2505.00653, Thm 1.5: 5/8. Both concern
      a fixed residue a (implied constants depend on a, so residues varying with N are not covered) and
      (log x)^C ≤ y ≤ x^{1/C} (verified by the literature critic).
  - EH limitations for primes: Friedlander–Granville 1989.
- **Expanders and smooth representatives.**
  - Jao–Miller–Venkatesan 2009: under GRH the Cayley graph on primes ≤ (log N)^{2+δ} is an expander, so every class
    is a product of O(log N) such primes. Prop GX is a corollary (folklore).
  - Unconditional analogues: Harman, Acta Arith. 91 (1999): for cube-free q, every reduced class contains a
    q^{1/(4√e)+ε}-smooth s ≤ q^{9/4+ε}. See also Harman–Shparlinski (IMRN 2016) and Munsch–Shparlinski.
- **Small-bias sets and spectral norm.**
  - Naor–Naor 1993; Alon–Goldreich–Håstad–Peralta 1992; Azar–Motwani–Naor 1998. ε-biased sets fool functions of small
    Fourier-L¹ norm; Prop FC(a)/(c) are instances.
  - Green–Sanders, Annals 168 (2008): the quantitative idempotent theorem (small-norm sets are combinations of
    cosets). Its bounds are useless at norm L^{O(1)}.
  - Shoup 1992: least primitive root under GRH. The generator set has norm ≈ 2^{ω(p−1)}, with polylogarithmic
    certificates.
- **Boundary case.** Kleinjung–Wesolowski (JAMS 2022): rigorous but randomized, fixed characteristic.
- **Derandomization.** Impagliazzo–Wigderson 1997 ⇒ deterministic L[1/2, O(1)] via Lenstra–Pomerance (polyimpl
  Thm DR; codex re-derived it). Not known to follow from GRH. Kabanets–Impagliazzo 2004.
- **Diagonal dominance:** Lévy–Desplanques (classical). **Wronskians of independent polynomials:** classical (Bôcher).

## 8. Experiments (pre-registered; toy sizes; correctness and mechanism only)

Pre-registration: `notes/zroute_prereg.md` (`986f17e`; addendum A1 `e66204c`). Summary:
`experiments/analysis/zroute_summary.json`.

**Disclosed deviations.**
- (1) The E4 run was queued behind E1 by a waiter process that died. It was restarted by hand after E1 finished,
  with identical code and config.
- (2) All E1 fields turned out to have signature (1, 1) (r_U = 1), so P1e tested only that case.
- (3) E4's D computation (below).
- No other deviations.

**E1: seedless deterministic NFS Z route (base-m cubic; 20 instances: 8 × 36-bit and 8 × 40-bit semiprimes, 4 ×
42-bit three-prime products).**

| prediction | result | verdict |
|---|---|---|
| P1a: 0 verification failures, 0 non-unit kernel vectors, 0 undecided certifications | 0 / 0 / 0 over 47 stages (50 780 output relations, each checked ev = 1 exactly) | **pass** |
| P1b: output rank = raw rational rank in ≥ 90% of stages | 46/47 | **pass** |
| P1c: PIN reached in ≥ 16/20 | **6/20** (36-bit: 5/8; 40-bit: 1/8; 42-bit tri: 0/4) | **FAIL** |
| P1d: PIN ⇒ Miller splits | 6/6 (largest base used: 3) | **pass** |
| P1e: pivots = r₁ + r₂ − 1 | 20/20 (all r_U = 1; see deviation 2) | **pass (narrow)** |

- *Why P1c failed (diagnosis from the recorded fields).*
  - In every failing instance the unpinned primes are exactly primes in Q₀ that never occur on the rational side of
    the prefix rows (1–16 primes per instance).
  - In 4 instances the fixed box ran out before the last stage: rows_available < n + C_alg + 7n.
  - This is the occurrence mechanism of Z3/E2, not a defect of XK-Z: output rank equals raw rational rank in 46/47
    stages.
  - My pre-registered window assumed the pilot (seed 999, 40 bits, PIN at multiplier 6.3) was typical. It was not.
- *Descriptive.* Miller split N in **41/41 stages without PIN** as well, and so in 20/20 instances at some stage. The
  pinned subgroup always contained a witness. This is luck at toy size, not a guarantee (Prop PIN needs all of Q₀).
- Cost: 30–340 s and ≤ 1.2 GB RSS per instance.

**E2: rank growth for degree-one rows (10 × 40-bit semiprimes, n = 243).**

| prediction | result | verdict |
|---|---|---|
| P2a: random labels reach full rank at exactly n rows in ≥ 9/10 | 10/10 | **pass** |
| P2b: deterministic orbit needs > 3n rows (or never) in ≥ 8/10 | 10/10; rows_to_full/n = 10.2–19.0 | **pass** |
| P2c: at the cap, n − rank = #non-occurring primes | not evaluable: every orbit run reached full rank before the cap (30n), so the deficit at the cap was 0 = #non-occurring primes, trivially | **vacuous** |
| P2d: PIN ⇒ Miller splits | 20/20 | **pass** |

- Post hoc (E3 below): most orbit rows are duplicates, which explains part of the 10–19n.

**E3: pointwise supply uniformity across N (1000 × 40-bit semiprimes, J = 2000 targets each).**

| prediction | result | verdict |
|---|---|---|
| P3a: mean orbit/mean random ∈ [0.9, 1.1] | 0.99996 (274.04 vs 274.05 hits) | **pass** |
| P3b: min/median ≥ 0.5 | orbit 0.542, random 0.810 | **pass** |
| P3c: variance/mean ≤ 3 | random 1.05; **orbit 4.27** | **FAIL (orbit)** |
| P3d: orbit cycles within J (descriptive) | 0/1000 | — |

- Per-N orbit and random yields are uncorrelated (ρ = 0.017), so no N-level effect is visible at this size.
- *Post hoc (`src/zroute/posthoc_e3.py`; 50 N, labelled).*
  - P(hit_{j+1} | hit_j) = 0.752, against a base rate of 0.132.
  - 9 537 of 9 926 consecutive hit pairs (96%) are the derived split (2u, v) or (u, v/2) of the previous target. They
    give the identical row (§5.3).
  - So the orbit's hits come in runs of duplicates: same mean, inflated variance, fewer independent relations.

**E4: deterministic rank forcing by dominance (Lemma DD1; addendum A1; 10 × 40-, 10 × 48-, 5 × 56-bit semiprimes).**

| prediction | result | verdict |
|---|---|---|
| P4a: complete ⇒ strictly dominant and rank n | 25/25 | **pass** |
| P4b: complete ⇒ Miller splits | 25/25 | **pass** |
| P4c: every prime gets a row (≥ 9/10 at 40, 48; ≥ 4/5 at 56) | 10/10, 10/10, 5/5 | **pass** |
| P4d: mean tries per prime at 40 bits ∈ [0.5/δ, 2/δ] = [3.65, 14.6] (δ = 0.137 from E3) | **30.7** | **FAIL** |

- Tries per prime are heavy-tailed. Medians of the per-instance medians: 17 / 47 / 153 at 40 / 48 / 56 bits. Means:
  30.7 / 80.1 / 273.9. Maxima: 6555 / 1455 / 2949.
- *Post hoc (`src/zroute/posthoc_e4.py`; one 40-bit instance; labelled).*
  - Per-prime median 14 against 1/δ ≈ 7.3; 35 of 242 primes needed > 50 tries.
  - Of 5 569 failed attempts examined, 1 674 (30%) had Gauss-reduced vector equal to the representation inherited
    from the target's own lattice vector, (wu₀, v₀)/g or (u₀, v₀/w).
  - So the multiplier family a·w (small w) is correlated with the target's own lattice. A bad target stays bad for many
    w. This explains part of the tail, not all of it.
  - Mechanistic lesson for DSR: the obvious explicit constructor "a·w for w = 1, 2, …" behaves worse than independent
    random targets, and pointwise cost per target is heavy-tailed.
- *Disclosed deviation (3).* The first E4 run was stopped at 48-bit seed 1 because flint's HNF in `pin_and_exponent`
  took ≈ 20 min on the 327 × 327 dominance matrix (the relation search took 0.6 s). I added a fast path: for a square,
  nonsingular input, D := |det|, a valid multiple of the lattice index. E4 was rerun in full.
  - Rows, tries and completeness are identical on the overlapping instances (`experiments/raw/zroute_e4_hnf_partial.jsonl`).
  - The predictions do not depend on the D computation.

## 9. Reviews

**GPT-6.1-sol blind round (`codex/zroute_1`, xhigh).**
- It reached independently:
  - XK-Z (composite evaluation, deterministic unit kernel, an elementary unit gap, the w_d trick for even degree, the
    rank formula);
  - "PIN suffices";
  - the pooling caveat;
  - hidden-label rank;
  - the GRH smooth-representative theorem (Prop GX).
- It contributed Lemma WR and the dominance idea (smooth compression). These are credited where used.
- Its probability estimate is ≈ 10%.

**GPT-6.1-sol adversarial round (`codex/zroute_2`, xhigh).**
- Verdicts:
  - CORRECT WITH FIXES: Theorem DZ at α = 1/2.
  - OVERCLAIMED: "all designs reduced to supply".
  - WRONG: RK ⇒ Q-rank; FC(d)'s lattice-point bound; the pooling inequality.
  - Vacuous: LV's large-sieve calculation.
- All fixes are applied in this version:
  - DSR made uniform and explicit; DSR-U for descent withdrawn and replaced by the literal DSR_{2/5};
  - rank claims restricted to DD1/WR designs;
  - the RK implication withdrawn;
  - FC(a) given an exact threshold; FC(c)'s heuristic replaced by an explicit speculative conjecture; FC(d)
    rewritten; LV demoted to failed templates plus a direct-application remark;
  - XK-Z given input-height hypotheses, the span-versus-basis issue fixed, and the sign variant qualified;
  - HR given the K bound and conditioning on targets, with codex's stronger hyperplane proof cited;
  - ZA's cost changed to (ω + 1)a and its halves-unit reason fixed;
  - the DSR/DL remark reduced to an analogy with no claimed reduction;
  - GX given the ceiling, invertible classes, the parameter range of the pair variant, and the "naive search"
    wording;
  - literature updated: Pascadi, Jain, Matomäki–Radziwiłł almost-all, Ge/Lenstra–Silverberg, Umans–Wang,
    Harvey–Hittmeir 2026's running time; the Evdokimov explanation removed.

**Opus proof critic (read-only).** ACCEPT WITH FIXES, no FATAL issue.
- Numerical checks: FC Parseval identities on 200 random (W, S); DD1 on 6 semiprimes, every row dominant, max |uv|
  = 0.988·N/√3, Miller split all; WR on 300 instances; HR exact at N = 35.
- One MAJOR issue: the DSR-U ⇒ L[2/5] clause. Fixed by redefining DSR_{2/5} with integer outputs.
- MINOR issues, all applied: ZA labelling and the halves-unit reason; HR's K bound and target measurability; FC
  general p and the (d) work bound; LV composite moduli and vacuity; the pooling inequality; WR's hypotheses
  (distinct labelled indices; H₀ bounding the actual INC output; the εe_{−1} correction folded into V); signed CRT in
  XK-Z(ii); Verdict labels.

**Opus proof critic, confirmation round on v2 (`0e3efca`).** ACCEPT. Every revised item holds. Its two wording nits
(μ restricted to primes ∤ N in GX; "work for generic w" in FC(e)) are applied.

**Opus literature critic (read-only).** OK WITH FIXES, all applied.
- Verified the four items reported only by codex: Harvey–Hittmeir's running time (adding D < N − 1), Pascadi's
  levels and setting, Jain (whose y-range was corrected), and Umans–Wang.
- Corrected:
  - Gao et al. are refereed (to appear in Math. Comp.);
  - Żrałek's oracle result was mischaracterised;
  - the McCurley 1988 citation was ambiguous and has been dropped.
- Added prior art:
  - GX is a corollary of JMV; Harman 1999 is the unconditional analogue;
  - HR's random-exponent rank argument is standard (Pomerance, Hafner–McCurley, Stange);
  - FC(a)/(c) are small-bias facts; Green–Sanders; Shoup 1992;
  - CFRAC, Vallée and Boneh are DSR analogues;
  - oracle equivalences (Bach–Miller–Shallit, Morain–Renault–Smith);
  - generic lower bounds do not apply;
  - Oznovich–Volk; Kleinjung–Wesolowski as the boundary case.

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
