# specialfam: deterministic factoring of explicit special families

Agent `specialfam`, phase 7, 2026-10-19. Pre-registration `notes/specialfam_prereg.md` (commit
`4be7188`, made before any run). Code `src/specialfam/`, tests `tests/test_specialfam.py`, raw data
`experiments/raw/specialfam_e*.jsonl`, summary `experiments/analysis/specialfam_summary.json`.
Starting point: codex `r18_ideas.out` §1–§2.

## 0. Verdict (plain language)

- **Yes, there are rigorous, checked, deterministic algorithms below N^{1/6} for explicit infinite
  families of N, and we prove them.** But below-1/6 special families are *cheap to come by*. A prime
  in a short interval or in a residue class with a large modulus already gives a factoring exponent
  near 0 by known Coppersmith methods (§6). So the value of this note lies in the general mechanism and
  in the careful accounting, not in "beating 1/6".
- **Main statement (Theorem G, proved; no prior formulation found by us).** Pollard–Strassen extends
  from intervals to polynomial values.
  - *Setting.* Φ is any integer polynomial of bounded degree, and the box has a side of length at
    least √(volume) along which Φ has a unit leading coefficient. All values on the (slightly
    extended) box must be nonzero and below N in absolute value.
  - *Conclusion.* If some prime of N divides some value of Φ on the box, a proper divisor of N is
    output deterministically in Õ(√volume).
  - *Engine.* The identity ∏_{s∈S} χ(s) = ±Nm(F_S mod χ) = ±Res(χ, F_S), plus one remainder tree of
    F_S modulo many small monic χ over (Z/N)[Y]. No ring with a norm is needed.
  - *Relation to r18.* r18's (Z/N)[i] construction is the special case Φ = a² + b².
  - *Boxes without a long side.* A norm form over an order (Theorem G2) restores √volume, but only for
    decomposable Φ.
- **Families (proved, via verified analytic inputs).**
  - p = a² + b² with b ≤ p^{0.119+o(1)} (Harman–Lewis): N^{0.15475+o(1)} balanced, N^{0.1238+o(1)} at
    p ≍ N^{2/5}. This confirms r18 §1.
  - Gaussian primes in small polar boxes (Harman–Kumchev–Lewis 0.53): N^{0.1325+o(1)} balanced,
    N^{0.106+o(1)} at N^{2/5}.
  - Coppersmith known-residue families via Linnik constants (Theorem S). The exponent is 3/20 with
    Xylouris's L = 5, which supersedes r18's 2/13 (L = 5.2). It is 1/24 for primorial moduli (Chang)
    and 0.0132 for 2-power moduli (Banks–Shparlinski). These families are "low bits known up to sign"
    in disguise.
  - In the Theorem A/G/G2 families (HL, HKL, polynomial forms), q can be any prime of the right size.
    In the Theorem S families, q is constrained too (q² ≡ 1 mod g).
- **One improvement over a 2026 result (Corollary W, proved; small; the γ < 1/12 part is conditional
  on GFHP Thm 5.3, an unrefereed 2025 preprint).** Let g = gcd(p−1, q−1) = N^γ for
  balanced N = pq. Then N factors deterministically from N alone in N^{E**(γ)+o(1)}, where
  - E** = 1/5 − 3γ/5 on [0, 1/18],
  - E** = 1/6 on [1/18, 1/12],
  - E** = 1/4 − γ on [1/12, 1/6],
  - E** = γ/2 on [1/6, 1/4],
  - and E** equals Ward–Maestro's curve beyond 1/4.

  This is strictly below Ward–Maestro's Prop. 28.2 curve E(γ) on (0, 1/5) \ {1/12}, with equality at
  1/12. It is also strictly below their record min(1/5, E) there. At γ = 1/6 it gives 1/12
  against their 1/6. It is a composition of WM's own recovery of g (Pollard–Strassen on N − 1) with
  published deterministic decoders: Coppersmith subdivision (GFHP 2025 Cor. 3.2, Oznovich–Volk 2025
  Thm 3.2, WM Cor. 14.3) for γ ≥ 1/12, and Gao–Feng–Hu–Pan's Theorem 5.3 for γ ≤ 1/12. The
  composition was not found stated by us.
- **What does not package (in our two models).** Non-decomposable Φ on boxes whose sides are all
  ≪ √volume, e.g. x² + y² + z² on a cube. This is a model statement, not a lower bound.
- **Prior art.**
  - *Theorem G.* No prior formulation for degree ≥ 2 moduli or polynomial boxes was found by us. Every
    piece is known:
    - remainder trees of a product modulo monic f_i (Bernstein, *Fast multiplication and its
      applications*, §18–19);
    - resultant = norm;
    - degree-1 moduli with arbitrary baby/giant sets (Hittmeir 2021 Def. 3.1 "revealing" pairs;
      Harvey 2021 Alg. 4.1);
    - quadratic-extension products in P±1 stage 2 (Montgomery–Kruppa 2008) as a close relative.
  - *Theorem S.* Its decoding step is published. GFHP 2025 Thm 3.1 and Oznovich–Volk 2025 Thm 3.2
    give deterministic known-residue Coppersmith subdivision. O–V's main algorithm has the same
    "learn p mod M, then decode" shape. What remains is only the discovery of the modulus and its
    mixed signs from N² − 1, which we did not find stated.
  - *Corollary W.* Not found, including in Ward–Maestro and GFHP. GFHP's Remark 5.5 compares the two
    decoders for a *known* modulus.
- **Experiments.** 576 instances (40–100 bits) were run as pre-registered, and all seven predictions
  held.
  - Every one of the 516 promised instances was factored; the controls gave 0 wrong factors.
  - Norm-box cost slopes match the predictions: 0.1543 vs 0.15475 (balanced) and 0.1244 vs 0.1238
    (β = 2/5). The wall-time slope is 0.184, against Pollard–Strassen's 0.25.
  - Corollary W's measured slopes at γ = 1/8 and 1/6 (0.161, 0.096) already sit below WM's E(γ)
    (0.1875, 0.1667). Details in §9.
- **Reviews.** GPT-6.1-sol could not run: Codex CLI login expired (401). An independent Opus critic
  (ACCEPT WITH FIXES, all applied) and an Opus prior-art auditor were used instead.

## 1. Setting and cost model

- N ≥ 2 is given. "Deterministic" means no randomness; correctness is unconditional.
- Every output is a verified proper divisor, so false outputs are impossible.
- Bit operations are counted with M(n) = O(n log n) for n-bit integer multiplication (Harvey–van der
  Hoeven 2021). Arithmetic in Z/N costs O(M(log N)).
- A product of two degree-< s polynomials over Z/N costs O(M(s(log N + log s))) by Kronecker
  substitution. For s ≤ N this is O(s log² N).
- Division with remainder of a degree-< 2s polynomial by a *monic* degree-s polynomial over any
  commutative ring A costs O(1) multiplications of size s. Use Newton iteration for the inverse of the
  reversed divisor modulo Y^s. Its constant term is 1, so no element of A is inverted (Bernstein,
  *Fast multiplication and its applications*, §6 and §17; also von zur Gathen–Gerhard, Ch. 9–10).
- Hence a product tree, and a remainder tree, over K monic polynomials of total degree D costs
  O(M(D log N) log K) bit operations. For D ≤ N this is O(D log³ N).
- AGENT_RULES: no factorization or primality routine is used by operational code. Flint is used only
  for multiplication and `%` by monic polynomials (the audited operations, `src/common/README.md`).

## 2. The packaging theorem (Theorem G)

### 2.1 Two ring lemmas (any commutative ring A)

**Lemma 1 (remainder trees; proved).** Let M = M₁M₂ with M₁, M₂ ∈ A[Y] monic. For every F ∈ A[Y],
(F mod M) mod M₁ = F mod M₁.

*Proof.*
- Division by a monic polynomial is unique over any commutative ring. If (Q₁−Q₂)M₁ = R₂ − R₁ with
  deg R_i < deg M₁ and Q₁ ≠ Q₂, the left side has degree deg(Q₁−Q₂) + deg M₁. This holds because its
  leading coefficient is lc(Q₁−Q₂)·1 ≠ 0. That contradicts the degree of the right side.
- Write F = QM + R and R = Q′M₁ + R₁. Then F = (QM₂ + Q′)M₁ + R₁ with deg R₁ < deg M₁. ∎

**Lemma 2 (block products are resultants; proved).** Let χ ∈ A[Y] be monic of degree d ≥ 1, and let
S ⊂ A be a finite list with F_S(Y) = ∏_{s∈S}(Y − s). Then

  ∏_{s∈S} χ(s) = (−1)^{|S|d} · Nm_{A[Y]/(χ)}(F_S mod χ),

where Nm is the determinant of multiplication on the free A-module A[Y]/(χ) with basis
1, Y, …, Y^{d−1}.

*Proof.*
- Both sides are polynomials with integer coefficients in the coefficients of χ and in the elements of
  S. F_S mod χ has integer-polynomial coefficients because χ is monic.
- So it suffices to prove the identity in Z[c, s] ⊂ K := Frac(Z[c, s]) (generic coefficients) and then
  specialise.
- Over an algebraic closure, generic χ has distinct roots θ₁, …, θ_d, and K̄[Y]/(χ) ≅ K̄^d by
  evaluation (CRT). Multiplication by r = F_S mod χ becomes diag(r(θ_i)) = diag(F_S(θ_i)), so
  Nm = ∏_i F_S(θ_i) = ∏_i ∏_s (θ_i − s).
- Also ∏_s χ(s) = ∏_s ∏_i (s − θ_i) = (−1)^{|S|d} ∏_i ∏_s (θ_i − s). ∎

For d = 1 this is plain evaluation. For d = 2 with χ = Y² + c₁Y + c₀ and r = r₀ + r₁Y,
Nm = r₀² − c₁r₀r₁ + c₀r₁² (used in the code). For larger d the code uses a division-free Laplace
determinant (d ≤ 6).

### 2.2 The block-product algorithm

**Proposition P (proved).** Given N, a list S of m integers and monic χ₁, …, χ_K ∈ (Z/N)[Y] of degrees
≤ d, the K values P_k = ∏_{s∈S} χ_k(s) mod N are computed with O((m + dK) log³ N) + O(K·c_d·M(log N))
bit operations (c_d = d! for Laplace, O(d⁴) for Berkowitz). No element of Z/N is inverted.

*Proof.*
- Build F_S (product tree of the monic Y − s), cost O(m log³ N).
- Build the product tree of the χ_k (all monic), cost O(dK log³ N).
- Reduce F_S modulo the root, then walk the remainder tree down. Lemma 1 makes every node equal to
  F_S mod (that node). Each tree level costs O(M((m + dK) log N)).
- Apply Lemma 2 at each leaf.
- In practice the giants are processed in batches of total degree ≈ m. This keeps every tree root at
  degree O(m), with the same total. ∎

### 2.3 The factoring theorem

**Theorem G (polynomial boxes with a long side; proved).**
- *Setup.* Φ(x, w) ∈ Z[x, w₁, …, w_r] has d = deg_x Φ ≥ 1 and x-leading coefficient c(w) ∈ Z[w]. The
  box is Q = I × W with I = [x₀, x₀ + A) ∩ Z and W ⊂ Z^r finite. Write V = A·|W|, and choose
  1 ≤ m ≤ A, I⁺ = [x₀, x₀ + m⌈A/m⌉).
- *Hypothesis (H).* For all x ∈ I⁺ and w ∈ W: 0 < |Φ(x, w)| < N and 0 < |c(w)| < N.
- *Conclusion.* There is a deterministic algorithm which, given N, Φ and Q, outputs a proper divisor
  of N whenever some prime p | N divides Φ(x, w) for some (x, w) ∈ Q. It never outputs anything else.
- *Cost.* With m = ⌈√V⌉ when A ≥ √V, it uses O_d(√V log³ N + K·h) bit operations, where K ≤ 2√V is
  the number of giants. Here h is the cost of writing down the coefficients of Φ(t + Y, w) for one
  (t, w). In general the cost is O_d((m + |W|⌈A/m⌉) log³ N + K·h).

*Proof.*
- **Giants.** Set G = {(v, w) : 0 ≤ v < ⌈A/m⌉, w ∈ W}, and K = |G| ≤ |W|(A/m + 1).
- **Leading coefficients.** For (v, w) ∈ G compute c = c(w) mod N and gcd(c, N).
  - By (H), 0 < |c(w)| < N, so the gcd is not N.
  - If 1 < gcd < N, output it (rule 3).
  - Otherwise c is a unit, and χ_{v,w}(Y) := c^{−1}Φ(x₀ + mv + Y, w) mod N is monic of degree d.
- **Block products.** Proposition P with S = {0, …, m−1} gives
  P_{v,w} = ∏_{j<m} χ_{v,w}(j) = c^{−m} ∏_{j<m} Φ(x₀ + mv + j, w) mod N.
- **Gcds.** Compute g_{v,w} = gcd(P_{v,w}, N) for all (v, w); this is K gcds. Since c is a unit,
  g_{v,w} > 1 iff some prime of N divides some leaf Φ(x₀ + mv + j, w).
  - If some prime p | N divides Φ(x, w) with (x, w) ∈ Q, then x = x₀ + mv + j for a unique (v, j),
    and g_{v,w} > 1.
  - Take the first giant with g > 1. If g < N, output it.
  - If g = N, then N | ∏_j Φ(x₀ + mv + j, w), so every prime of N divides some leaf. Scan the m leaves
    of this block. The first leaf L with gcd(L, N) > 1 gives a proper divisor, because 0 < |L| < N.
    The scan costs O(m) gcds and runs at most once.
- **Cost.** O_d((m + dK) log³ N) for the trees, plus K(h + small determinant + gcd) for the per-giant
  work. When A ≥ √V and m = ⌈√V⌉, K ≤ √V + |W| ≤ 2√V, since |W| = V/A ≤ √V. ∎

**Remarks.**
1. **Degree 1 is Pollard–Strassen.** With Φ = x + 0·w (or the arbitrary-set variant
   χ_k(Y) = Y + t_k), Theorem G is Pollard–Strassen with arbitrary baby and giant sets. The content is
   the degree ≥ 2 moduli: the giant data enter a *small monic modulus*, not a big fixed polynomial,
   so multivariate giant sets cost nothing extra.
2. **The algorithm needs only Φ and the box.** In particular it needs no field, no ring with a norm,
   and no structure of Φ beyond a unit leading coefficient in one long direction.
3. **Every polynomial is unit-leading after a change of variables.**
   - Take a primitive u ∈ Z^k with Φ_top(u) ≠ 0 and complete it to a unimodular basis. In the new
     coordinates the x₁-leading coefficient is the constant Φ_top(u).
   - It must be a unit mod N: if 1 < gcd < N that is a factor, and N | Φ_top(u) is excluded when
     |Φ_top(u)| < N.
   - The real restriction is geometric. The box must be covered by O(√V) lattice lines parallel to
     u, each meeting it in a segment of length ≥ √V.
   - The segments may start at different points. Theorem G extends trivially to a starting point
     x₀(w) that depends on the giant (the code does not implement this).
4. **(H) is the price of the scan.** It is needed on the extended blocks I⁺ × W, not only on Q. It
   holds when Φ > 0 on I⁺ × W with values below N, which is the case in every family below.

### 2.4 Boxes without a long side: norm forms (Theorem G2) and what fails

**Lemma B (box splitting; proved).** Every box Q = ∏_{i≤k}[c_i, c_i + L_i) ∩ Z^k of volume
V = ∏L_i is contained in S₁ + S₂. Here S₁ is a sub-box at the origin and S₂ is a grid, with
|S₁| ≤ 2√V and |S₂| ≤ 2√V.

*Proof.*
- Let κ be the first index with ∏_{i≤κ} L_i ≥ √V, and set m_i = L_i for i < κ.
- Set m_κ = ⌈√V / ∏_{i<κ}L_i⌉ ≤ L_κ.
- Take S₁ = ∏_{i<κ}[0, L_i) × [0, m_κ) and S₂ = c + (m_κ·[0, ⌈L_κ/m_κ⌉)) × ∏_{i>κ}[0, L_i).
- Then |S₁| ≤ √V + ∏_{i<κ}L_i < 2√V.
- Also |S₂| ≤ (L_κ/m_κ + 1)∏_{i>κ}L_i ≤ √V + √V. ∎

**Theorem G2 (norm forms; proved).**
- *Setup.* Let g ∈ Z[T] be monic of degree n, O = Z[T]/(g), Nm: O → Z the determinant of
  multiplication, and L: Z^k → O affine.
- *Hypothesis.* 0 < |Nm L(x)| < N on S₁ + S₂.
- *Conclusion.* For any box Q of volume V, a deterministic algorithm finds a proper factor of N
  whenever some prime of N divides Nm(L(x)) for some x ∈ Q, in O_n(√V log³ N) bit operations.

*Proof.*
- Split Q ⊆ S₁ + S₂ by Lemma B, and write L(s + t) = L₀(s) + L(t) with L₀ the linear part.
- Over the free (Z/N)-algebra O_N = O/NO, form F(X) = ∏_{s∈S₁}(X + L₀(s)) ∈ O_N[X] (a monic product
  tree).
- Evaluate F at the |S₂| points z = L(t) by the remainder tree with monic linear moduli X − z
  (Lemma 1 over A = O_N). This gives F(z) = ∏_s L(s + t).
- Nm is multiplicative and commutes with reduction mod N, so Nm(F(z)) ≡ ∏_s Nm L(s + t) (mod N).
- Finish with gcds and the block scan as in Theorem G.
- Each O_N operation costs O(n²) operations of Z/N, and Nm costs O(n³). ∎

This is r18's mechanism in general form. For n = 2 and O = Z[i] it is literally r18 §1. For Z[i] it
coincides with Theorem G: (Y + mv)² + b² is the characteristic polynomial of −(mv + ib).

**Proposition C (exactly which Φ are norm-packageable with affine L; proved).** For Φ ∈ Z[x₁, …, x_k]
the following are equivalent:
- (a) there are an integer c ≠ 0, a finite product R = ∏ Z[T]/(g_i) of monic orders, and an affine
  L: Z^k → R with c·Φ = Nm_R ∘ L;
- (b) Φ is decomposable: it factors into affine-linear polynomials over Q̄.

*Proof.*
- **(a ⇒ b).** R ⊗ Q̄ is a finite product of local Artinian Q̄-algebras A_j with residue maps π_j.
  On a local A_j, multiplication by y = π_j(y) + n (n nilpotent) has determinant π_j(y)^{dim A_j}.
  So Nm_R(y) = ∏_j π_j(y)^{dim A_j}, and each π_j ∘ L is affine.
- **(b ⇒ a).** By unique factorisation in Q̄[x], the multiset of affine factors of Φ, up to scalars,
  is Galois-stable. In each Galois orbit, normalise one factor ℓ so that one of its coefficients is 1.
  Then ℓ has coefficients in a number field K = Q(θ), θ an algebraic integer, and the orbit's product
  is Nm_{K/Q}(ℓ(x)) times a rational constant.
- Clear denominators: Dℓ has coefficients in Z[θ], and Nm(Dℓ) = D^{[K:Q]}Nm(ℓ).
- Repeated orbits give repeated factors of R.
- The result is ∏ Nm(D_oℓ_o)^{e_o} = ρΦ with ρ = a/b ∈ Q (a, b integers). Add one factor Z to R with
  the constant affine map L ≡ b. Then Nm_R(L(x)) = b·ρΦ = aΦ, so c = a. ∎

Theorem G2 holds verbatim for such products R. Work componentwise in R/NR, take Nm_R = ∏ of the
component norms, and use a common baby set. Each component is a monogenic Z[T]/(g_i); a single
Z[T]/(∏g_i) would not do over Z. If a number field's basis ω_i is not in a monogenic order (e.g.
Maynard's), use Z[θ] with Dω_i ∈ Z[θ]; the norm is multiplied by D^n.

**So which structures package at √V?**

| structure | √V packaging? | via |
|---|---|---|
| interval (Pollard–Strassen) | yes | G (d = 1) |
| one-variable values f(x), x ≤ X | yes (√X) | G, or a_d^{d−1}f(x) = Nm_{Z[θ]}(a_d x − θ) with θ = a_d·(root) |
| a² + b², any box with a long side | yes | G (= r18, G2 with Z[i]) |
| any binary quadratic form ax² + bxy + cy² (any D, also square D) | yes, any box | G2: a·Φ = Nm_{Z[ω]}(ax + ωy) with ω² − bω + ac = 0 |
| binary forms of any degree (a_d ≠ 0, else swap variables) | yes, any box | G2: a_d^{d−1}f(x,y) = Nm(a_d x − (a_dθ)y), monic g(T) = a_d^{d−1}f(T/a_d, 1) |
| incomplete norm forms (Maynard), x³ + 2y³ (Heath-Brown) | yes, any box | G2 (O = Z[θ]) |
| x² + y⁴ (Friedlander–Iwaniec), a² + k ("near a square"), x² + y² + z² **with a long side** | yes | G (not norm forms; G2 does not apply) |
| non-decomposable Φ, no side ≥ √V (e.g. x² + y² + z² on a cube [0,L)³) | **not packaged by G, nor by G2 with affine L**: G costs V/L = V^{2/3}; G2 does not apply (rank-3 forms are not products of linear forms). Non-affine L and non-commutative norms (quaternions) are outside both models | — |

Comment: multivariate multipoint evaluation (e.g. Kedlaya–Umans) does not obviously help. The product
∏_{s∈S}Φ(s, ·) has ≥ |S|^{k₂} coefficients in k₂ ≥ 2 giant variables, which is ≥ V at a balanced
split. This is a model observation only; it is **not a lower bound**.

## 3. Theorem A (norm boxes, r18 §1)

**Theorem A (proved).**
- *Statement.* Fix a rational η ∈ (0, 1/2). There is a deterministic algorithm with public inputs N
  and η, implemented as `normbox_factor`. Suppose N is composite with a prime factor p = a² + b², where
  a ≥ 0 and 1 ≤ b ≤ p^η. (A factor below 11 is returned by trial division.) Then the algorithm returns a proper divisor of N in
  O_η(p^{1/4+η/2} log³ N) bit operations.
- *Corollary.* If p ≍ N^β, the cost is N^{β(1/4+η/2)+o(1)}.

*Proof.*
- **Bands.** Bands k = 4, 5, … use A_k = ⌊2^{k/2}⌋ and B_k = ⌊2^{kη}⌋. Theorem G is applied with
  Φ = a² + b², x = a ∈ [0, A_k], W = [1, B_k] and m = ⌊√((A_k+1)B_k)⌋ + 1. The algorithm stops at the
  first proper divisor. It also stops when (2A_k+1)² + B_k² ≥ N.
- **The promised band.** Let k be the band with 2^{k−1} < p ≤ 2^k. Then a ≤ √p ≤ 2^{k/2} and
  b ≤ p^η ≤ 2^{kη}, so (a, b) lies in band k's box.
- **(H) holds.** If N has a prime factor < 11, trial division returns it. Otherwise N/p > 1 has all
  its prime factors ≥ 11, so 11p ≤ N.
  - Leaves are ≥ 1. The x-range of the blocks is [0, m⌈(A+1)/m⌉) ⊆ [0, 2A+1], because m ≤ A + 1.
  - So leaves are ≤ (2A+1)² + B² ≤ 8p + 4√(2p) + 1 + 2p < 11p for p ≥ 37. Bands 4–5 were checked by
    hand (90 < 143, 146 < 187). The critic checked every band k ≤ 29 against the least eligible p.
  - Every earlier band has smaller leaves. So the loop reaches band k, and Theorem G guarantees success
    there or earlier.
- **Cost.** Band k costs O(√((A_k+1)B_k) + B_k) giants plus babies, which is O(2^{k(1/4+η/2)}) since
  B_k ≤ A_k + 1. Summed over bands ≤ k, this is O_η(2^{k(1/4+η/2)}) = O_η(p^{1/4+η/2}). Multiply by the
  log³ N tree factor. ∎

**Comparison with r18 §1.** All of r18's algorithmic claims check out:
- the cost Õ(m + AB/m);
- monic division without inverses (Lemma 1);
- leaf norms below N. The precise condition is "no prime factor < 11", or a known β.

The gcd = N recovery is the block scan of Theorem G. A product-tree descent is not needed.

## 4. The analytic inputs (verified against the sources)

**Harman–Lewis (HL).** G. Harman and P. Lewis, *Gaussian primes in narrow sectors*, Mathematika 48
(2001) 119–135. The PDF scan was read by us.
- **Theorem 2.** Let X > X₀. For any β, γ with 0 ≤ β ≤ π/2 and X^{−0.381} ≤ γ ≤ π/2, the number of
  Gaussian primes p with β ≤ arg p < β + γ and |p|² ≤ X is > cXγ/log X, with c an absolute positive
  constant.
  - The width is measured against the norm bound X (sector area ≍ X^{0.619}).
  - It is uniform in the position β.
  - It is a lower bound of the right order.
  - The norm range is "≤ X", not dyadic.
- **Theorem 1** as printed reads "infinitely many primes p = m² + n², n < p^θ, θ < 0.119". We do not
  rely on this wording. The proof shows (6) for θ = 0.119; the inequality sign looks like a
  typesetting slip. We derive what we need from Theorem 2 directly.

**Lemma HL (proved from HL Thm 2).**
- *Statement.* There are absolute C, X₀ such that for X ≥ X₀ there are ≥ (c/2)X^{0.619}/log X primes
  p ∈ (X/(C log X), X] with p = a² + b², a ≥ 1 and 1 ≤ b ≤ X^{0.119}.
- *Consequence.* For every η > 0.119, the family 𝓕_η = {p prime : p = a² + b², 1 ≤ b ≤ p^η} meets
  every interval (X/(C log X), X] for X large.

*Proof.*
- Apply Theorem 2 with β = 0 and γ = X^{−0.381}.
- Discard Gaussian primes on the real axis. These are rational primes ≡ 3 (mod 4) with norm p² ≤ X,
  at most √X of them.
- Discard those of norm ≤ X/T. Lattice points with 0 ≤ arg < γ and |z|² ≤ X/T satisfy 0 ≤ b < 2γa,
  so there are ≤ Σ_{a≤R}(2γa + 1) ≤ γX/T + γR + R of them, with R = √(X/T).
- With T = 4 log X/c, at least (c/2)Xγ/log X remain. Here Xγ = X^{0.619} ≫ √X.
- Each remaining π = a + ib has a, b ≥ 1 (1 + i has argument π/4 > γ), so its norm p is a rational
  prime. Moreover b = |π| sin(arg π) ≤ √X·γ = X^{0.119}.
- Finally p > X/T gives b ≤ (Tp)^{0.119} = p^{0.119+o(1)}. ∎

**Harman–Kumchev–Lewis (HKL).** G. Harman, A. Kumchev and P. A. Lewis, *The distribution of prime
ideals of imaginary quadratic fields*, Trans. AMS 356 (2004) 599–620. Read by us.
- **Theorem 2.** Set A = {𝔞 ∈ C : x − y ≤ N𝔞 < x, φ₀ ≤ arg λ(ξ_𝔞) ≤ φ₀ + φ}, with
  λ(ξ) = (ξ/|ξ|)^w. If x ≥ x₀, y ≥ x^{0.765} and φ ≥ x^{−0.235}, then Σ_{𝔭∈A} 1 ≫ φy/log x.
- **Theorem 1.** It yields Q(m − s, n − t) ≪ Q(s, t)^{0.53} + 1.
- **Corollary.** For every z there is a Gaussian prime π with |z − π| ≪ |z|^{0.53} + 1.

**Lemma HKL (proved from HKL Thm 2, Q(i), w = 4).**
- *Statement.* Take x = 2^k and the public centre z₀ = √(x/2)(1 + i). For x large there are
  ≫ x^{0.53}/log x primes p = a² + b² ∈ [x − x^{0.765}, x) with |a + ib − z₀| ≤ 2x^{0.265}.
- *Hence* (a, b) lies in the public square of side 4x^{0.265} + O(1) around (√(x/2), √(x/2)).

*Proof.*
- Apply Theorem 2 with φ₀ = π and φ = x^{−0.235}. The generator in the first quadrant then has argument
  in [π/4, π/4 + φ/4], so its radial deviation from √x is ≤ y/√x = x^{0.265}.
- The arc length is ≤ √x·φ/4 = x^{0.265}/4.
- a, b ≥ 1, so the norm is a rational prime. ∎

**Other prime-producing polynomial families (status of the inputs).**
- Heath-Brown, x³ + 2y³ (Acta Math. 186 (2001)): primes with x, y ≤ X ≍ p^{1/3}. Verified only by
  standard citation.
- Maynard, *Primes represented by incomplete norm forms* (arXiv 1507.05080; Forum Math. Pi 2020):
  N_K(Σ_{i≤n−k} x_iω^{i−1}) has the expected number of prime values when n ≥ 4k. For pure radical
  fields Q(θ^{1/n}) it has infinitely many prime values when n ≥ 22k/7. The abstract was read.
- Friedlander–Iwaniec, x² + y⁴ (Ann. Math. 148 (1998)). Standard citation.
- Harman–Lewis Theorem 4: infinitely many p with {√p} < p^{−0.262}, i.e. 0 < p − ⌊√p⌋² ≤ 2p^{0.238}.

## 5. Theorem S (smooth common part + Coppersmith, r18 §2)

**Coppersmith input (cited).** May, *Using LLL-reduction for solving RSA and factorization problems*
(2010), Theorem 1 (stated in §2, proved in §3, with c ≥ 1).
- *Statement.* Let N have a divisor b ≥ N^β and let f be monic of degree δ. All x₀ with
  f(x₀) ≡ 0 (mod b) and |x₀| ≤ c·N^{β²/δ} are found in time O(c δ⁵ log⁹ N).
- *How it is used.* Only with δ = 1, as repeated calls on shifted windows of width N^{β²}. One call
  per chunk.
- *Code.* The code uses the exact rigorous Howgrave-Graham + LLL condition with lattice dimension
  ≤ 9 (`copper.py`, from trackE).

**Lemma R (root counting; proved).** For H ≥ 1 let ρ(d) = #{r mod d : r² ≡ 1}. Then
Σ_{d|H} ρ(d) ≤ 2τ(H)². If H | N² − 1, this is N^{O(1/log log N)}.

*Proof.*
- ρ is multiplicative.
- For odd ℓ, (r−1)(r+1) ≡ 0 (mod ℓ^e) with gcd(r−1, r+1) | 2 forces r ≡ ±1, so ρ(ℓ^e) = 2.
- ρ(2) = 1, ρ(4) = 2, and ρ(2^e) = 4 for e ≥ 3.
- Hence ρ(d) ≤ 2·2^{ω(d)} ≤ 2τ(d), and Σ_{d|H}τ(d) ≤ τ(H)².
- τ(n) = n^{O(1/log log n)} (Wigert), and H < N². ∎

**Lemma W (divisor window; proved).** If g ≥ D and P⁺(g) ≤ B, then g has a divisor d with
D ≤ d < BD.

*Proof.* Remove prime factors of g one at a time while the quotient stays ≥ D. The final d satisfies
d/ℓ < D for some prime ℓ | d with ℓ ≤ B. ∎

The code tries only d ∈ [D, BD), largest first. The theorem does not need Lemma W.

**Theorem S (proved; the decoding step is published).** The decoder in Step 2 is exactly GFHP 2025
Thm 3.1 / Oznovich–Volk 2025 Thm 3.2 (deterministic known-residue Coppersmith subdivision, with the
β(1−β) dyadic form in O–V). Theorem S adds only Step 1 and the r² ≡ 1 enumeration.
- *Setup.* Let p | N be prime with p ∈ [P, 2P], P = N^β (the
  band is public). Suppose some integer g satisfies
  - g | p² − 1 and g | N² − 1 (for N = pq this is g | gcd(p² − 1, q² − 1)),
  - g ≥ N^γ,
  - P⁺(g) ≤ B = N^η,

  with B and the bound N^γ public.
- *Conclusion.* A deterministic algorithm (`smoothcopp_factor`) outputs a proper divisor of N in
  N^{o(1)}·(N^{η/2} + N^{[β(1−β)−γ]₊}) bit operations.
- *Unknown band.* Scanning the O(log N) dyadic bands gives N^{max(η/2, 1/4−γ)+o(1)}.

*Proof.*
- **Step 1.** Pollard–Strassen block products over Z/M for M = N − 1 and M = N + 1 (Proposition P
  with χ_v = Y + (3 + mv), m ≈ √B). This finds every prime ℓ ≤ B dividing M, with multiplicity, in
  Õ(√B).
  - Scan blocks in order with a shrinking cofactor M_cur.
  - The first x with gcd(x, M_cur) > 1 is the least prime factor of M_cur in range, so no primality
    test is needed.
  - Let H be the B-smooth part of (N−1)(N+1) = N² − 1. Then g | H, because g | N² − 1 and every
    prime factor of g is ≤ B.
- **Step 2.** For every d | H with d ≥ N^γ and every r mod d with r² ≡ 1, write p = r + dx. The x
  range has length ≤ 2P/d + 1.
  - Make the polynomial monic: f(y) = y + (rd^{−1} + x_c) (mod N). Here d is a unit mod N since
    d | N² − 1.
  - Cover the x range by ⌈(2P/d + 1)/N^{β²}⌉ windows, with one Coppersmith call each (δ = 1, b = p ≥ N^β).
- **Correctness.** For d = g (or any divisor of g that is ≥ N^γ) and r = p mod d, which satisfies
  r² ≡ 1 (mod d), some window contains (p − r)/d. The call returns p. Every output is gcd-verified.
- **Cost.** By Lemma R there are ≤ 2τ(H)² = N^{o(1)} pairs. Each costs
  O((N^{β−γ−β²} + 1)·log⁹ N). ∎

**Comparison with r18 §2.**
- r18 stated max{η, …}. The smooth part costs √B with Pollard–Strassen, not B, so the first term is
  η/2.
- The root bound 2τ(H)² is correct.

**Corollary W (arbitrary common factor; proved).**
- *Statement.* Let N = pq be balanced (p < q ≤ Cp) and let g := gcd(p − 1, q − 1) = N^γ. Then
  `common_factor_factor`, given N alone and with no knowledge of g or γ, factors N in
  N^{max(γ/2, 1/4−γ)+o(1)} bit operations.
- *Unbalanced and mixed-sign versions.* The same holds with p ≍ N^β, giving max(γ/2, β(1−β)−γ),
  and with g | gcd(p² − 1, q² − 1) (residues r² ≡ 1).

*Proof.*
- **Dovetailing.** For t = 2, 3, …: take B_t = 2^t and compute the B_t-smooth part H_t of N − 1 in
  Õ(2^{t/2}). For every divisor d of H_t not tried before whose window count is ≤ 2^{⌈t/2⌉}, run the
  Coppersmith windows with r = 1.
- **When it succeeds.** g | N − 1 (because N − 1 = (p−1)q + (q−1)) and every prime of g is ≤ g. So
  once 2^t ≥ g and 2^{t/2} ≥ #windows(g) ≍ √C·N^{1/4}/g, the divisor d = g is tried and p ≡ 1 (mod g)
  is found.
- **Cost.** Stage t costs N^{o(1)}·2^{t/2}, by Lemma R for the number of divisors. The first
  successful stage has 2^{t/2} ≍ max(√g, N^{1/4}/g). ∎

**Against Ward–Maestro 2026, Prop. 28.2 (their N-only curve for the same class).**
- WM's curve is E(γ) =
  - 1/4 − γ on [0, 1/12],
  - 1/8 + γ/2 on [1/12, 1/8],
  - 1/4 − γ/2 on [1/8, 1/6],
  - 1/2 − 2γ on [1/6, 1/5],
  - γ/2 on [1/5, 1/4],
  - 1/4 − γ/2 on [1/4, 3/8],
  - 1/(8√e) + ε on [3/8, 1/2).

  It comes from their Route A (Lemma 26.2 / Cor. 27.1), Route B (McKee–Pinch progression BSGS with
  Harvey–Hittmeir order selection) and Route C (progression trial; Hinek Thm 6.3 for g ≥ N^{1/4}).
- Combining with Corollary W gives E*(γ) = min(E(γ), max(γ/2, 1/4 − γ)), which is
  - 1/4 − γ on [0, 1/6],
  - γ/2 on [1/6, 1/4],
  - WM's branches beyond.
- E* < E exactly on (1/12, 1/5). The gap is largest at γ = 1/6: 1/12 against 1/6.
- With the Coppersmith decoder alone, the threshold for beating 1/5 stays at γ > 1/20.
- WM's report states known-residue Coppersmith subdivision itself: their Corollary 14.3, with cost
  N^{λ−µ−1/4}; λ = 1/2 and µ = γ give 1/4 − γ. Their Part 4 never applies it to the recovered g.
- **Small γ: Gao–Feng–Hu–Pan (GFHP), arXiv 2512.19076, Theorem 5.3 (read by us).**
  - *Statement.* If p ≡ r (mod n) for every prime p | N, then N is factored deterministically in
    O(N^{1/5}log^{16/5}N / n^{3/5} + N^{1/6}log⁴N).
  - *Applies here.* With n = g recovered from N − 1 (cost √g, dovetailed as above) and r = 1, it gives
    max(γ/2, 1/6, 1/5 − 3γ/5).
  - *Comparison.* That is below 1/4 − γ exactly for γ < 1/12 (GFHP Remark 5.5 makes the same
    crossover remark for a known modulus).
- **Final curve.** E**(γ) = min(E(γ), max(γ/2, 1/4−γ), max(γ/2, 1/6, 1/5−3γ/5)), which is
  - 1/5 − 3γ/5 on [0, 1/18],
  - 1/6 on [1/18, 1/12],
  - 1/4 − γ on [1/12, 1/6],
  - γ/2 on [1/6, 1/4],
  - WM beyond 1/4.

  E** < E on (0, 1/5) except at γ = 1/12, where they are equal. E** < 1/5 for every γ ∈ (0, 1/2). WM
  need γ > 1/20 to beat 1/5 and γ > 1/6 to reach below 1/6 (critic correction), while E** ≤ 1/6 for
  γ ≥ 1/18.
- The composition (g recovered from N − 1, then the GFHP or Coppersmith decoder) was not found
  stated by us: "no prior formulation found by us". Its proof is the dovetailing argument above, with
  GFHP Thm 5.3 as a cited, unrefereed (2025 preprint) input for the γ < 1/12 part.

**Families for Theorem S (Linnik-type inputs).**
- *Construction.* Let g be any modulus with P⁺(g) ≤ (log N)^{O(1)}, e.g. g = 2^k or a primorial. Let
  p and q be the least primes ≡ 1 and ≡ −1 (mod g), so q ≠ p. Then g | gcd(p² − 1, q² − 1).
- *Exponent.* If p, q ≪ g^L then N ≪ g^{2L}, so γ ≥ 1/(2L) − o(1). The cost is
  N^{1/4 − 1/(2L) + o(1)} for every β, since β(1−β) ≤ 1/4.
  - **L = 5.2** (Xylouris 2009, arXiv 0906.2749, the value r18 used): 2/13 ≈ 0.1538.
  - **L = 5.18** (Xylouris, Acta Arith. 2011): ≈ 0.1535.
  - **L = 5** (Xylouris, *Linniks Konstante ist kleiner als 5*, Chebyshevskii Sbornik 2018,
    published): **3/20 = 0.15**.
  - **Prime-power moduli ℓ^k, ℓ fixed:** Banks–Shparlinski (J. Anal. Math. 139 (2019) 239–263;
    journal version, which states L < 2.1115) give a lower bound of the right order for π(x; ℓ^k, a)
    when ℓ^k ≤ x^{0.4736}, with constants depending on ℓ. The arXiv v2 does not contain the Linnik
    improvement. Verified via the journal abstract and Baier–Pujahari, arXiv 2107.04348.
    So with g = 2^k there are p, q ≪ g^{2.1115}, γ ≥ 0.2368 and cost **N^{0.0132+o(1)}**.
  - **Under GRH** (L = 2 + ε): N^{ε}.
  - **Smooth moduli:** Chang (J. Anal. Math. 123 (2014)), Cor. 11, verified by the prior-art
    auditor. If log ℓ = o(log q) for every prime ℓ | q, the least prime ≡ a (mod q) is
    < q^{12/5+o(1)}. Primorials qualify, so g = primorial gives γ ≥ 5/24 and cost
    **N^{1/24+o(1)}**.
  - Caveat on L = 5: an erratum note (arXiv 2001.00578) says the L ≤ 5 proof still needs checking.
    5.18 (Acta Arith. 2011) is safe and gives 0.1535.
- *Honest reading.* These families mean "p ≡ ±1 modulo a known large 2-power or smooth modulus", i.e.
  ≈ γ·n low bits of p known up to sign. That is classical Coppersmith territory (Coppersmith 1996;
  Boneh–Durfee–Howgrave-Graham 1999). The only non-classical ingredient is the *discovery* of the
  modulus and of sign patterns from N² − 1.

## 6. Calibration: special families are cheap; where norm boxes sit

The relevant mechanism depends on the shape of the information about p.

| information on p | mechanism | cost |
|---|---|---|
| p in a public interval of length L | Coppersmith subdivision | N^{o(1)}·(L/N^{β²} + 1) |
| p ≡ r (mod d), d ≥ N^γ | Coppersmith | N^{o(1)}·(N^{β−γ−β²} + 1) |
| p divides a polynomial value on a box of volume V | Theorem G / G2 | V^{1/2+o(1)} |

- Baker–Harman–Pintz put primes in every [x − x^{0.525}, x]. So the explicit family
  "p ∈ [2^k − 2^{0.525k}, 2^k]" is factored in N^{0.525β − β² + o(1)}: N^{0.0125} balanced and
  N^{0.05} at β = 2/5. This is far below everything in §7 and is prior art (Coppersmith 1996).
- The norm-box families are interesting because they are *not* localised in size or residue.
  - The HL sector family has p anywhere in (X/(C log X), X].
  - It is contained in {p : p − ⌊√p⌋² ≤ p^{0.238}}, a union of √X short intervals. Coppersmith
    on that superset needs one call per interval, N^{β/2} in total, while Theorem A gives N^{0.3095β}.
- For the HKL family the norm interval [x − x^{0.765}, x] is also known. Coppersmith on it gives
  N^{0.765β−β²}, which ties Theorem G at β = 1/2 (0.1325) and loses at β = 2/5 (0.146 vs 0.106).

**Why T1/T2 (below16) do not apply, and why there is no tension with them.**
- *Different model.* T1/T2 are bounds on the *fraction of S_n* covered by residue-local degree-1
  collision designs in the 𝔐⁺ group model, with labels charged. Theorems A, G and S are not
  collision searches in (Z/N)^*. They are polynomial evaluation over (Z/N)[Y] and lattice reduction.
  So T1/T2 say nothing about them.
- *Inside the envelope anyway.* Even transported formally, they sit strictly inside the T1 envelope.
  - At cost N^{1/8+η/4}, Theorem A covers a fraction ≈ p^{η−1/2} = N^{η/2−1/4} of balanced semiprimes.
    T1 allows n²·N^{3(1/8+η/4)/2−1/4} = N^{3η/8−1/16+o(1)}. Since η/2 − 1/4 < 3η/8 − 1/16 for all
    η < 3/2, we are below it.
  - Theorem S / Corollary W classes have density ≈ Σ_{g ≥ N^γ} g^{−2} = N^{−γ+o(1)}.
    - For γ < 1/4 the cost is N^{1/4−γ}. T1 allows N^{1/8−3γ/2}, and −γ < 1/8 − 3γ/2 iff γ < 1/4.
    - For γ ≥ 1/6 (Corollary W cost N^{γ/2}) T1 allows N^{3γ/4−1/4}, and −γ < 3γ/4 − 1/4 iff
      γ > 1/7.
    - So the classes are inside the envelope throughout (critic correction).
- *Summary.* These are sparse-family results. They do not touch the average-case floor.

## 7. Family table

The cost exponent is log_N(cost). "e_p" is the exponent in p (cost = p^{e_p+o(1)}). The table gives
the balanced value (β = 1/2) and the β = 2/5 value. In the rows using Theorem A/G/G2, q is any prime
of the right size. In the Theorem S / Corollary W rows q is constrained (q² ≡ 1 mod g). All rows are
proved (Theorems A/G/G2/S plus the cited analytic input) unless marked.

| family (explicit, infinite) | analytic input | mechanism | e_p | β = 1/2 | β = 2/5 |
|---|---|---|---|---|---|
| none (worst case) | — | Harvey–Hittmeir | — | 0.2 | 0.2 |
| none | — | Pollard–Strassen | 1/2 | 0.25 | 0.2 |
| p = a² + b², b ≤ p^{0.119+o(1)} | HL Thm 2 | Thm A | 0.3095 | **0.15475** | **0.1238** |
| (a, b) within 2p^{0.265} of public √(x/2)(1+i) | HKL Thm 2 | Thm G | 0.265 | **0.1325** | **0.106** |
| x³ + 2y³ | Heath-Brown 2001 | G or G2 (box X×X) | 1/3 | 0.1667 | 0.1333 |
| Maynard, pure radical, n = 22, k = 7 | Maynard (n ≥ 22k/7); box shape from the abstract, UNVERIFIED | G2 (no long side) | 15/44 | 0.1705 | 0.1364 |
| Maynard, general K, n = 4k | Maynard (n ≥ 4k); box shape UNVERIFIED | G2 | 3/8 | 0.1875 | 0.15 |
| x² + y⁴ | Friedlander–Iwaniec | G (long x side) | 3/8 | 0.1875 | 0.15 |
| p − ⌊√p⌋² ≤ 2p^{0.238} | HL Thm 4 | G (a² + k) | 0.369 | 0.1845 | 0.1476 |
| p ∈ [x − x^{0.525}, x] (comparator) | BHP | Coppersmith | — | 0.0125 | 0.05 |
| p ≡ 1, q ≡ −1 mod g, least primes | Xylouris L = 5 (5.18 safe) | Thm S | — | 0.15 (0.1535) | ≤ same |
| same, g = primorial | Chang Cor. 11 (12/5) | Thm S | — | 1/24 | ≤ 1/24 |
| same, g = 2^k | Banks–Shparlinski 0.4736 | Thm S | — | 0.0132 | ≤ 0.0132 |
| gcd(p−1, q−1) = N^γ, g arbitrary | (promise) | Cor. W (+ GFHP 5.3 for γ < 1/12) | — | E**(γ) | max(γ/2, 6/25−γ) for γ ≥ 1/12 |

Notes on the table:
- The Maynard rows use V = X^{n−k} and p ≍ X^n, so e_p = (n−k)/(2n).
- The HKL and HL rows are the best Gaussian rows we know. We found no "every sector" improvement of
  0.381, and none of 0.53; almost-all results do not give explicit families.

## 8. Prior art

Sources: our WebSearch rounds plus an independent Opus prior-art auditor (§10), who read Hinek's
thesis, Hittmeir 2021, May 2010, HL, HKL, Chang, Banks–Shparlinski, Bernstein's multapps, GFHP,
Oznovich–Volk, BDF and WM §§25–28 in full.

**Pollard–Strassen with arbitrary sets (degree-1 moduli): KNOWN.**
- Strassen 1977 and Pollard 1974.
- Hittmeir, Math. Comp. 90 (2021), Def. 3.1, Alg. 3.2 and Lemma 3.3 ("revealing" pairs B, S; form
  ∏_{b∈B}(X − b) and evaluate at S).
- Harvey, Math. Comp. 90 (2021), Alg. 4.1.
- Bostan–Gaudry–Schost 2007 and Costa–Harvey 2014 need arithmetic-progression shifts.

**Remainder trees of products modulo monic polynomials: KNOWN.** Bernstein, *Fast multiplication and
its applications* (MSRI Publ. 44, 2008): §6 Newton reciprocal over any commutative ring, §17 quotient
and remainder, §18 remainder tree for monic f_i, §19 "small factors of a product" (∏h_j mod f_i). With
h_j = Y − s and f_i = χ_k this is exactly the F_S mod χ_k step of Proposition P.

**Close relative.** Montgomery–Kruppa, *Improved stage 2 to P±1* (ANTS VIII, 2008) uses product and
remainder trees in (Z/N)[X]/(X² − aX + 1) with resultant-style products.

**Theorem G as a factoring statement: no prior formulation found by us.** That is: degree ≥ 2 monic
moduli giving Õ(√V) for primes dividing Φ on a polynomial box; Pollard–Strassen over (Z/N)[i] at
Gaussian points; the norm-form multi-dimensional version; the HL/HKL family exponents. Searches
covered Pollard–Strassen over Gaussian integers, extension rings and resultants; "prime factor of the
form f(x)"; and p = a² + b² with small b.

**Cryptanalysis of structured primes.** Gaussian-integer RSA variants (Elkamchouchi et al.; Teşeleanu
2024) concern the key equation with (p²−1)(q²−1), not structured factors. Cheng's 4p−1 / CM methods
use a different mechanism. Hittmeir 2017 treats a^n ± b^n. We found no attack on factors with a small
Gaussian coordinate, or on factors near a square.

**Known-residue Coppersmith (Theorem S's decoder): KNOWN, deterministic.**
- Coppersmith 1996; Boneh–Durfee–Frankel 1998, Cor. 2.2 (p mod r with r ≥ N^{1/4} gives polynomial
  time); Coppersmith–Howgrave-Graham–Nagaraj, Math. Comp. 2008.
- Gao–Feng–Hu–Pan 2025 (arXiv 2512.19076): Thm 3.1 / Cor. 3.2 give N^{1/4}/m; Thm 5.3 gives
  N^{1/5}/n^{3/5} + N^{1/6} when every prime is ≡ r (mod n).
- Oznovich–Volk 2025 (arXiv 2506.07668): Thm 3.2 gives N^{1/(4r)−α}, with the same
  "learn p mod M, then decode" shape.
- WM Cor. 14.3.
- ROCA (Nemec et al. 2017) enumerates residues modulo a divisor of a primorial and runs
  Howgrave-Graham per guess.
- Boneh–Durfee–Howgrave-Graham 1999 combine exhaustive search with lattices.

**Large gcd(p−1, q−1).**
- Wiener; McKee–Pinch 1996/98 (not retrievable; quoted via Hinek).
- Hinek, PhD thesis (Waterloo 2007): Attack 6.1 gives N^{1/4−γ/2} expected; Thm 6.3 gives polynomial
  time for known g ≥ a + b; Thm 6.4 gives N^{1/4−γ} expected for known g; §6.3.4 seeks g with ECM/NFS.
- Ward–Maestro 2026 Part 4 (deterministic, N only; E(γ)).
- **Not found:** recovering g from the smooth part of N − 1 or N² − 1 and then running a
  known-residue decoder, i.e. Corollary W and the mixed-sign promise of Theorem S.

**Verdict on novelty, at most "no prior formulation found by us", for:**
- Theorem G (and Theorem A as its instance, i.e. r18 §1);
- the HL/HKL family exponents;
- Corollary W's curve E**;
- the N² − 1 mixed-sign promise.

Theorem S itself is a routine corollary of published decoders.

## 9. Experiments (pre-registered; run exactly as frozen)

Commands: `.venv/bin/python -m src.specialfam.run_experiments e1|e2c|e3|e4|e5|summarize`.

All runs used the frozen code at `4be7188`.
- **One disclosed deviation.** The summarizer lacked the P7 statistic (median of max(ps_units, chunks)).
  It was added after the runs; this is analysis-only, and the data are unchanged.
- **Lattice dimension.** The Coppersmith lattice dimension is capped at 9. The rigorous window is then
  N^{0.222} rather than N^{1/4} at β = 1/2, so the chunk slopes in E4/E5 sit above the asymptotic
  exponents by about 1/36 (critic issue 15). The pre-registered predictions use the finite-dimension
  values.

| id | result | prediction | verdict |
|---|---|---|---|
| E1 correctness | 180/180 factored (90 balanced, 90 at β = 2/5); all 180 hits STRUCTURAL-PROMISED (the block of the promised (a, b)); 0 gcd = N blocks | 100%; ≥ 90% promised | P1, P4 hold |
| E2 count slope | β = 1/2: 0.1543 (band units), 0.1550 (total); β = 2/5: 0.1244 / 0.1225 | 0.15475 ± 0.01; 0.1238 ± 0.01 | P2 holds |
| E2 wall slope (β = 1/2) | 0.184 (upper half), 0.197 (all); median 3.5 s at 100 bits (0.31 s at β = 2/5) | [0.13, 0.21]; Pollard–Strassen would be 0.25 | P3 holds |
| E2c control | 60 unstructured balanced semiprimes: 2 factored, both at n = 40, both with p = a² + b², b ≤ p^{0.119} (structural, as expected at 20-bit p); 0 wrong | 0 wrong, ≤ 10% found | P5 holds |
| E3 | HKL-type box 20/20; near-square a² + k 20/20 | 100% | P1 holds |
| E4 (Thm S, mixed CRT signs) | 128/128 factored; chunk slope 0.162 (γ = 0.12), 0.105 (γ = 0.18); ≤ 39 (d, r) pairs | 0.153 ± 0.03; 0.108 ± 0.03 | P6 holds |
| E5 (Cor. W, g = 2·prime) | 168/168 factored; slope of max(PS, chunks): 0.190, 0.161, 0.096, 0.100 for γ = 1/10, 1/8, 1/6, 1/5 | 0.175, 0.153, 0.108, 0.100 (± 0.035); WM E: 0.175, 0.1875, 0.1667, 0.100 | P7 holds |

Reading E5:
- *Where E5 is informative.* The finite-dimension measured slopes at γ = 1/8 and 1/6 (0.161, 0.096)
  are already below WM's asymptotic E(γ) (0.1875, 0.1667). This is consistent with Corollary W.
- *What it is not.* Toy sizes (48–96 bits) are correctness and slope checks only (rule 8). The
  asymptotic exponents are the theorems, not these fits.
- *Not tested.* The GFHP branch of Corollary W (γ < 1/12) was not implemented or tested.

## 10. Reviews

**GPT-6.1-sol rounds could not run.** Three attempts were made: `specialfam_2` prior art twice, and
`specialfam_1` adversarial review once. All three failed at startup with Codex CLI authentication
errors (refresh token reused, then expired; HTTP 401). Logs are in
`codex/scratch/specialfam_{1,2}_attempt*.log`. The prompts are committed as `codex/specialfam_{1,2}.md`
and can be re-run after `codex login`. Two independent Opus subagents substituted for them.

**Opus critic (read-only): ACCEPT WITH FIXES, no fatal issue.**
- *Checks performed.*
  - `block_products` against direct evaluation (300 random cases);
  - `smooth_part_ps` against sympy (300 cases);
  - `unit_roots` against brute force;
  - the Theorem A leaf inequality for bands 4–29;
  - WM's E(γ) recomputed from their Routes A/B/C.
- *Verdicts.*
  - Lemmas 1, 2, B, R, W, Proposition P, Theorems G, G2, A, S and Corollary W: correct.
  - Remark 3 and Proposition C: fixable.
- *Fixes applied.*
  1. (major) "q any prime" was restricted to the A/G/G2 rows.
  2. WM need γ > 1/6, not 1/8, to get below 1/6.
  3. The density of the Theorem S class in the T1 comparison is N^{−γ}, not N^{−2γ}; the case split
     was added.
  4. The Theorem A leaf chain was fixed for p < 37, and the "no factor < 11" hypothesis was folded
     into the algorithm.
  5. Remark 3: unit condition, line-cover geometry, and a giant-dependent x₀.
  6. Proposition C: rational ρ; G2 for products of orders; non-monogenic bases.
  7. Table nits: a_d ≠ 0, monic f, and the scope of the cube example.
  8. Cost K·h.
  9. (H) on I⁺ × W.
  10. §0 wording (rule 5; one divisor output; (H) mentioned).
  11. WM Cor. 14.3 cited.
  12. The HL lattice count.
  13. "Contained in" the near-square set.
  14. Maynard box shape marked UNVERIFIED.
  15. The lattice-dimension cap disclosed in §9.

**Opus prior-art auditor: two findings changed the note.**
- (i) The deterministic known-residue Coppersmith subdivision is published (GFHP 2025 Thm 3.1;
  Oznovich–Volk 2025 Thm 3.2). Theorem S was downgraded to "routine corollary of published decoders".
- (ii) GFHP Thm 5.3 combined with the recovery of g improves the curve on (0, 1/12). This was added to
  Corollary W as E**, after we read Thm 5.3 and Remark 5.5 ourselves.
- It also confirmed HL, HKL, Xylouris L = 5 (with an erratum caveat), Banks–Shparlinski (journal
  version), Chang Cor. 11 and May 2010 Thm 1. It found no every-sector or every-disc improvement of
  0.381 or 0.53.
- The von zur Gathen–Gerhard theorem numbers could not be checked. Bernstein's multapps (§6, 17–19)
  is cited instead.

## 11. Claims and labels

| claim | label |
|---|---|
| Lemmas 1, 2, B, R, W; Proposition P | proved (standard facts, proofs given) |
| Theorem G (any Φ, long side, √V) | proved; no prior formulation found by us |
| Theorem G2 (norm forms, any box) and Proposition C (exact characterisation) | proved |
| Theorem A (norm boxes; r18 §1 confirmed with corrected leaf condition) | proved |
| Lemma HL (b ≤ p^{0.119+o(1)}, from HL Thm 2); Lemma HKL | proved from the cited theorems (read by us) |
| family exponents 0.15475/0.1238 (HL), 0.1325/0.106 (HKL) | proved (conditional only on the cited published theorems) |
| Theorem S (η/2 improves r18's η) | proved; routine corollary of published decoders (GFHP Thm 3.1, O–V Thm 3.2) |
| Corollary W, Coppersmith branch (improves WM Prop. 28.2 on (1/12, 1/5)) | proved; no prior statement found by us |
| Corollary W, GFHP branch (γ < 1/12) | proved conditional on GFHP Thm 5.3 (2025 preprint, read by us, unrefereed) |
| Linnik families: 3/20 (L = 5; 5.18 safe), 1/24 (Chang), 0.0132 (Banks–Shparlinski) | proved from cited results (Xylouris L = 5 has an open erratum note) |
| non-decomposable Φ without a long side is not packaged by G or by G2 with affine L | model observation, not a lower bound |
| Maynard rows (box shape) | UNVERIFIED (abstract only) |
| E1–E5 slopes and correctness | experimentally supported |
