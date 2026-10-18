# nfsfruit: fruitfulness of square congruences in a rigorous randomized NFS

Agent `nfsfruit`, 2026-10-18. Pre-registration: `notes/nfsfruit_prereg.md` (E1 committed in `a3523e9`,
Addendum A in `a3b8852`, both before their runs). Code: `src/nfsfruit/`. Tests: `tests/test_nfsfruit.py`.
Codex: `codex/nfsfruit_1.{md,out}` (blind) and `codex/nfsfruit_2.{md,out}` (adversarial).

## Verdict

**No rigorous L[1/3] factoring theorem was obtained, not even for a meaningful restricted family.** The
fruitfulness gap in Lee–Venkatesan (LV) remains open. The study does narrow it to one precisely stated
statement about a single explicit function β on the set of smooth pairs, and the experiments find no
sign that it fails. §4 gives a short hierarchy of sufficient conditions: the batch-level Conjecture F,
which is the weakest, and the static Lemma M. Feasibility of closing it with current tools: **low**. The needed statement is a
"mod-4 / second-order" equidistribution property that the first-order (mod-2, quadratic-character)
equidistribution tools behind LV do not reach.

What *is* established (label in brackets; the critic status is in §9):

1. **[proved] Lemma H.** For every odd composite n, Blum or not, fruitfulness is a group homomorphism
   Φ: K0 → μ₂(Z/n)/{±1} on the space of true square dependencies. So the argument of LV's conditional
   Theorem 2.3 works for all semiprimes under the single hypothesis "Φ ≢ 0 on the dependency space
   actually sampled". LV's Remark 7.3 says that for non-Blum n
   it is "unclear how (even notionally)" to argue this. Codex reached the same lemma blind, and it is
   probably folklore.
2. **[proved] Bockstein formula.** On K0, Φ equals Σ_{i∈S} β_i for an explicit per-pair function β, which
   is half of a character that vanishes on every relation. For Blum n, in the quadratic sieve and in
   degree-1 NFS, β is computable from n alone via Jacobi symbols. So the "factor-detecting character" of
   r16 §3.3 is a literal, computable extra row there. It also gives a degeneracy theorem for the two-form
   degree-1 sieve: relations whose L1L2 values have T-weighted prime count ≡ 0, 1 (mod 4), where
   T = {ℓ : (ℓ/n) = −1}, can never be fruitful (in particular when all their primes have (ℓ/n) = +1). In QS
   every factor-base prime has (ℓ/n) = +1, and fruitfulness comes instead from J_n(x_i).
3. **[proved] Reduction to a static statement about one function.** Let δ(n,m,f) be the distance from β
   to the nearest function (a,b) ↦ ⟨w, x(a,b)⟩, under the law of kept pairs on the smooth set X_{n,m,f}
   (x = exact square-class vector, including the constant feature when |S| parity is imposed).
   - With all of X used, Φ ≢ 0 ⟺ δ > 0.
   - For a batch of M i.i.d. kept pairs, P(Φ ≡ 0 on the batch kernel) ≤ 2^C(1−δ)^M.
   - The weakest sufficient statement is the batch-level **Conjecture F**: under the law of
     (m, f, batch) produced by LV's search, Pr[Φ ≢ 0] ≥ L[1/3,−o(1)].
   - A static sufficient condition is **Lemma M**: δ ≥ L[1/3,−o(1)] on a set of (m,f) of measure
     L[1/3,−o(1)] under the search-output law. Batches must then be enlarged to ≥ C ln2/δ, at cost L^{o(1)}.
   - Either condition yields LV-time L[1/3,(64/9)^{1/3}+o(1)] factoring of every semiprime.
     [proved implication, modulo LV's non-monic bookkeeping]
   - An x-preserving, β-flipping involution on a non-negligible part of X gives Lemma M; Dixon's
     ε-multiplication is exactly such an involution. Lemma M is pointwise *stronger* than Conj. 7.1. Its
     value is that it removes the batch randomness and asks one static question.
4. **[proved; known in substance] Hidden-root multiplier (degree 1).** Choose the multiplier λ = r² of a
   two-linear-form sieve with r uniform. Then every odd-weight dependency is fruitful with probability
   exactly 1 − 2^{1−ω}. As the critic notes, this is Dixon in disguise: put x_i = r·L2(a_i,b_i); then
   x_i² ≡ L1·L2(a_i,b_i) (mod n), with one shared hidden root. The forms have coefficients of size about √n,
   and the heuristic cost is L[1/2, √2] [assessment]. It illustrates that a hidden root is exactly a
   Dixon-type fibre symmetry. It is not a new algorithm.
   **[proved, counting only]:** a globally uniform hidden root admits a ratio representative λ ≡ c/d with
   |c|,|d| ≤ H with probability ≤ 4H²·2^ω/φ(n). That this excludes the naive transfer to NFS-sized
   coefficients is an [assessment].
   It is *not* a general barrier to hidden-root ideas at L[1/3] (codex round 2 rejected my stronger
   wording).
5. **[experimentally supported]** Toy exact experiments covered Dixon, QS, degree-1 NFS (fixed and
   randomized), degree-3 NFS (base-m and LV-randomized), Blum, mixed, both-1-mod-4 and 3-prime n, at
   40 and 44 bits. Results:
   - 1120 + 1120 instances, 0 structural degeneracies;
   - the factor-detecting functional is nonzero on every instance with dim K0 ≥ 10;
   - prefix spaces of dimension k are fruitless at rate 2^{−k} within |z| ≤ 3.3 in every family;
   - no Blum / non-Blum difference;
   - (the toy "LV-randomized" family is monic with deg R = d − 2, so it is LV-*style*, not LV's distribution);
   - Lemma H: 4248 linearity checks, 0 failures;
   - the Bockstein identity: 102 693 exact checks, 0 failures;
   - the hidden-root algebra: 28 029 checks, 0 failures; odd dependencies fruitful 0.521 (n=238),
     3-prime 0.775 (n=80).

The restricted-family option is met only trivially. Codex points out that Fermat's method factors
n = pq with q − p ≤ n^{1/4}L[1/3, c/2] deterministically in L[1/3, c], and that this family has
≫ X^{3/2}/log²X members. (Coppersmith's method does better on the same family, L[1/3, c/2] — critic.) That is [proved], but it is the close-prime family and says nothing about NFS.
No non-trivial special family with a provable L[1/3] was found.

---

## 1. Lee–Venkatesan, read from the arXiv PDF (1805.08873v1)

### 1.1 The randomized algorithm and its randomness space
The parameters are δ, κ, σ, β, β′ with κ > δ^{-1}, the constraint 2σ > max(β,β′) + δ^{-1}/(3β) + (σδ+κ)/(3β′),
and δ^{-1} < κ + σδ/2. The degree is d = δ(log n/log log n)^{1/3}, and d is **odd** (needed in Lemma 6.10;
it also makes μ(K) = {±1}). The random choices are:
- **m** uniform in [2^{−1/d}L, L], where L = L_n(2/3, δ^{-1}).
- **f** = f̂_{n,m} + Σ_{i=0}^{d−1} c_i (x − m y) x^{d−i−1} y^i. Here f̂ is the base-m expansion of n, and the
  c_i are i.i.d. uniform in a centred interval of length 2L_n(2/3, κ − δ^{-1}). The c_i term dominates
  f̂ because κ > δ^{-1}. There are **no random multiples of n**; the randomness is "Kleinjung-like" but
  with random R.
- **(a,b)** uniform with 0 ≤ a < |b| ∈ [½L(1/3,σ), L(1/3,σ)]. A pair is kept iff a − bm is
  B = L(1/3,β)-smooth and f(a,b) is B′ = L(1/3,β′)-smooth. Smoothness is tested rigorously with the
  hyperelliptic curve method (Cor. 3.19).
- **Stochastic deepening** (Lemma 5.6): a random number 2^i of (m,f) pairs is tried, each searched to depth
  2^{−i}·L(1/3, …). This replaces second-moment bounds.
- **Quadratic characters** come from IdealSampler: a random degree bound k ∈ [d], a random r near x^{1/k} with
  log x = L(1/3, ε), ε = (1/3+o(1))δ, Miller–Rabin, a random simple factor of f mod r. They are large
  primes, chosen to avoid possible Siegel zeros (Lemmas 6.6, 6.10). Their number is
  4d(δκ log n + δ²κ² log 2 · log^{4/3} n (log log n)^{1/3}).
- The kernel is **sampled uniformly** (block Wiedemann). ℓ = (4/3) log log n kernel vectors are produced, and
  two lying over the same element of the small exceptional *set* E_{K,ε} ⊂ H are paired. The §3.3
  pseudocode (checked by the critic) shows that each kernel vector comes from a **fresh batch of about C
  relations for the same (m,f)**. The relation list is emptied after each kernel vector. *Mismatch:* the
  proof of Thm 2.3 instead argues with all of P_S ("running for L(1/3,2σ) finds every generator of P_S").
  So the published conditional proof does not literally match the batch algorithm. Prop R below addresses
  this; batches must be enlarged to ≥ C ln2/δ.

### 1.2 Theorems (as stated)
- **Thm 2.1.** For any n, the Randomised NFS runs in expected time L_n(1/3, (64/9)^{1/3} + o(1)) and
  produces x, y with x² ≡ y² (mod n). The Coppersmith multiple-polynomial variant gives
  L(1/3, ((92+26√13)/27)^{1/3}).
- **Thm 2.3.** For fixed semiprime n with both primes ≡ 3 (mod 4), conditional on Conjecture 7.1, it finds
  x² ≡ y² with x ≢ ±y in expected time L(1/3, (64/9)^{1/3} + o(1)).
- **Thm 2.5** (relations). Under (4.1)/(4.2) it finds, almost surely, an irreducible f of degree d and height
  ≤ L(2/3,κ) with n | f(m), together with L(1/3, max(β,β′)+o(1)) pairs a < |b| ≤ L(1/3,σ) that are smooth on
  both sides. This takes expected time L(1/3,λ) for any λ > max(β,β′) + δ^{-1}/(3β) + (σδ+κ)/(3β′), and
  the failure probability is ≤ L(2/3, κ−δ^{-1})^{−1+o(1)}.
- **Thm 2.6** (squares). For all but an L(2/3, (κ−δ^{-1})/2·(1+o(1)))^{-1} fraction of f, given
  L(1/3, max(β,β′))·Ω(log log n) smooth pairs, a congruence of squares is found in expected time
  L(1/3, 2max(2δ/3, β, β′))^{1+o(1)}.

### 1.3 The conjecture (verbatim content, §7)
Let n = pq with p ≡ q ≡ 3 (mod 4), and χ_n = χ_pχ_q. Define
P = {h ∈ Z[X] : (h(m),n) = 1, h(m) is a square in Z, h(α) = g² with g ∈ Z[α]}, and P_S as the smooth part:
products of linear factors a − bX of *height* L(1/3,σ), smooth on both sides. That height set is not the
same as the sampled box 0 ≤ a < |b| ∈ [½L, L]. Define
χ_P(h) = χ_n(√h(m))·χ_n(φ_{m,α}(√h(α))). It is well defined because χ_n(−1) = 1, and it extends to an
order-2 character on Z[X].
**Conjecture 7.1:** χ_P restricted to P_S cannot be written as a product of the characters χ_𝔭 and
(−1)^{ord_𝔭}. Every such character is trivial on P_S, which consists of squares, so the conjecture is
equivalent to **χ_P ≢ 1 on P_S**, i.e. Φ ≢ 0 on the true dependency space of all smooth linear factors of
that height. This holds under the literal definition P_S ⊂ P. Under the reading LV's proof uses (P_S as the
group of smooth products), it depends on the chosen extension unless the characters cut out K0 exactly.
The proof of Thm 2.3: since χ_P is a nonzero linear functional on the kernel, a uniform kernel element
is fruitful with probability ½. Remark 7.2 supports the conjecture only by counting (dim P_S exceeds the
number of characters). Remark 7.3 notes that if it fails for (n,f), NFS with those parameters never
splits n.

### 1.4 Why the standard argument is not available
In Dixon's algorithm the sampled root x is uniform on its fiber x·μ₂(Z/n), conditional on everything the
algorithm uses (z = x² mod n and its factorization). So each relation carries an independent uniform
"root sign", and any dependency is fruitful with probability 1 − 2^{1−ω} (§2.1 of codex/nfsfruit_1.out).
In QS, linear sieve and NFS there is **no hidden randomness at all**. Both square roots,
u = √∏(a−bm) ∈ Z and v = √∏(a−bα) ∈ Z[α], are canonical up to ±1. So fruitfulness is a deterministic
function of (n, m, f, relation set, S). Randomizing m, f, (a,b) does not help unless some measure-preserving
action on the conditioned data changes the root comparison, and none is known. Our §4 makes this precise:
the missing object is an x-preserving, β-flipping involution.

## 2. Literature, 2018–2026 and the rigorous lineage

Sources: my own searches plus codex/nfsfruit_1.out §1 (citations there, spot-checked).
- **Rigorous NFS.** The only follow-ups found are van Leeuwen's MSc thesis (Bristol 2020, arXiv 2007.02689;
  an exposition of LV, with no new fruitfulness result) and citations in Le Gluher–Spaenlehauer–Thomé 2021
  (heuristic NFS complexity refinement), Dąbrowski–Pomykała–Shparlinski 2023 (oracle factoring) and
  Bernstein–Chou CryptAttackTester 2024 (App. B.3 flags Conj. 7.1 as the remaining condition), plus a
  2025 DCC paper on NFS complexity formulae. **No proof of fruitfulness, or of square-root randomness, for
  NFS / QS / linear sieve was found.**
- **Rigorous L[1/3] for special families.** None found for integer factoring. Adjacent DLP results are
  genuinely rigorous but do not transfer: Granger–Kleinjung–Zumbrägel 2018 and Kleinjung–Wesolowski
  (JAMS 2022, quasi-polynomial DLP in fixed characteristic).
- **Rigorous L[1/2, c].**
  - Dixon 1981;
  - Pomerance 1987: L[1/2, √2];
  - Vallée 1989/91: L[1/2, √(4/3)], rigorous, using quasi-uniform small squares ~n^{2/3};
  - Seysen 1987: L[1/2, √(5/4)] under GRH;
  - A. K. Lenstra 1988: L[1/2, 1] under GRH;
  - **Lenstra–Pomerance 1992: L[1/2, 1+o(1)], unconditional** — still the record;
  - Schnorr–Lenstra 1984 relies on class-number smoothness assumptions.
  **No unconditional c < 1 was found.**
- **Class groups.** de Boer–Pellet-Mary–Wesolowski (arXiv 2512.01588, 2025/26) give ERH-rigorous
  subexponential class groups for arbitrary number fields. This is relevant machinery for S-unit
  computations (our β), but its cost is L_{|Δ_K|}(1/2), not L_n[1/3].
- **Making squares.** Croot–Granville–Pemantle–Tetali (Annals 2012) and Balister–Bollobás–Morris
  (Annals 2018) settle when square products appear, but not root signs.

## 3. Lemma H and the Bockstein formula

Setting: K = Q(α), f monic irreducible of odd degree, f(m) ≡ 0 (mod n), φ: Z[α] → Z/n with α ↦ m.
Relations i are coprime (a_i,b_i) with z_i = a_i − b_i m. Let
K0 = {S : ∏_S z_i ∈ Q^{×2}, ∏_S (a_i − b_iα) ∈ K^{×2}}, an F₂-space (the kernel of the map into
Q^×/Q^{×2} ⊕ K^×/K^{×2}).
For S ∈ K0, u = f′(m)√∏z_i ∈ Z and v = f′(α)√∏(a_i−b_iα) ∈ Z[α], using f′(α)O_K ⊆ Z[α].
Then u² ≡ φ(v)² (mod n), and everything is assumed coprime to n (otherwise a factor is output).

Standing assumptions:
- gcd(n, f′(m)∏z_i) = 1 (otherwise a factor is output);
- "semiprime" means p ≠ q, with p, q > max(B, B′);
- f is monic, for exposition. LV's f is non-monic (leading coefficient C_d + c_0 random). The statements
  carry over with LV's normalisation: f_d(a − bm), f_d(a − bα), the constant |S|-parity column, and a
  BLP-type f′ factor. The constant feature must then be included in x. This is routine but not written out here.

**Lemma H [proved].** Φ(S) := [u/φ(v)] ∈ μ₂(Z/n)/{±1} is well defined and Φ(S△T) = Φ(S) + Φ(T).
*Proof.* u and v are determined up to sign, so the class is well defined. Put I = S∩T. Square roots for
S△T are u(S)u(T)/∏_I z_i and v(S)v(T)/∏_I(a_i−b_iα), up to sign. Since φ(a_i − b_iα) = z_i, the extra
factors cancel under φ. ∎
For odd n with ω prime factors, μ₂(Z/n) ≅ {±1}^ω, so μ₂/{±1} ≅ F₂^{ω−1} (non-canonically), and S is fruitful
iff Φ(S) ≠ 0. A uniform S in K0 (zero included) is therefore fruitful with probability exactly 1 − 2^{−rank Φ}.
A basis vector or a uniform nonzero vector does not have this distribution.
**Consequence:** the argument of Thm 2.3 works for every semiprime under the hypothesis "Φ ≢ 0 on the
dependency space actually sampled"; the Blum condition is not needed. LV's Conj. 7.1 is phrased with
Blum characters, so its non-Blum analogue is this hypothesis, not 7.1 itself. **Prior art:** the same
subgroup mechanism ("dependencies giving trivial factorizations form a proper subgroup", plus uniform
kernel sampling) is explicit in the class-group setting in Seysen, Math. Comp. 48 (1987), p. 777 before
Prop. 8.6 and Thm 8.7 (located by codex round 2). Lemma H is the NFS transcription and is probably folklore. The same
proof works for any two-sided sieve whose relations satisfy r_i ≡ φ(β_i) (degree 1, QS, multiple
polynomials). Codex nfsfruit_1 §2.3 proves the same statement blind and adds the well-definedness
caveat that the root must be evaluable, which holds here via f′(α).

**Proposition G (Bockstein formula) [proved].** Let n = pq. Choose χ_p ∈ \widehat{F_p^×} of order
2^{v₂(p−1)}, so that χ_p(−1) = −1, and likewise χ_q. Let c = max(v₂(p−1), v₂(q−1)) and ζ = e^{2πi/2^c}.
Define ψ(z) = χ_p(z)χ_q(z) on Z coprime to n, and ψ(x) = χ_p(x mod 𝔭)χ_q(x mod 𝔮) on K^× coprime to
𝔭𝔮, where 𝔭 = (p, α−m) and 𝔮 = (q, α−m). Then ψ(−1) = 1, and for x ∈ Z[α], ψ(x) = ψ(φ(x)).
Fix a basis g_1..g_R of U_S/{±1}, where U_S is the S-unit group (S = factor-base primes and infinite places).
Write z_i = ±∏ℓ^{e_iℓ} and a_i − b_iα = ±∏g_j^{c_ij}. Choose integer lifts τ_ℓ, σ_j with
ψ(ℓ) = ζ^{τ_ℓ} and ψ(g_j) = ζ^{σ_j}. Put t_i = Σ τ_ℓ e_iℓ − Σ σ_j c_ij. Then 2^c | t_i for every relation, and
  **Φ(S) = Σ_{i∈S} β_i (mod 2),  where β_i = t_i/2^c mod 2, for all S ∈ K0.**
*Proof.* First, ψ(z_i) = ψ(φ(a_i−b_iα)) = ψ(a_i−b_iα), so ζ^{t_i} = 1. Next, write u/φ(v) ≡ (s_p, s_q) ∈ {±1}²
by CRT. Since χ_p(−1) = −1, ψ(u)/ψ(φ(v)) = χ_p(s_p)χ_q(s_q) = s_ps_q, which is −1 iff fruitful.
The f′ factors cancel because ψ(f′(m)) = ψ(f′(α)). Since S ∈ K0, all E_ℓ = Σ_S e_iℓ and Σ_S c_ij are even,
√∏z = ±∏ℓ^{E_ℓ/2} and √∏(a−bα) = ±∏g_j^{Σ_S c_ij/2} exactly (μ_K = ±1 as d is odd), and ψ(±1) = 1.
Hence ψ(u)/ψ(v) = ζ^{(Σ_S t_i)/2} = (−1)^{Σ_S t_i/2^c}. ∎
Changing the lifts or the basis changes β by an F₂-linear function of the square-class vector
x_i = (e_iℓ mod 2, c_ij mod 2). More generally, any two valid β, including those from a different choice of
χ_p among the characters of order 2^{v₂(p−1)}, agree on K0(X), so by duality they differ by some ⟨w,x⟩ on X.
β is therefore well defined modulo rowspace. β is a fixed function of
(n, m, f, (a,b)). It depends on hidden data (χ_p, χ_q) and on S-unit coordinates.

**Corollary J (Blum, computable rows) [proved; exhaustively checked in E2].** For Blum n we have c = 1 and
ψ = J_n (the Jacobi symbol), so everything below is computable from n:
- (J1) **QS** (z_i = x_i² − n, X_i = x_i). Every odd factor-base prime ℓ | z_i has (n/ℓ) = 1, so J_n(ℓ) = 1.
  Hence Φ(S) = Σ_{i∈S} β_i with β_i = [J_n(x_i) = −1] + [n ≡ 5 (mod 8)]·[x_i odd]. (For n ≡ 5 (mod 8) and
  x odd, v₂(x² − n) = 2; the code uses the equivalent form v₂(z_i)/2.)
  *Proof:* Φ(S) = [J_n(X·Y) = −1] with Y = ∏ℓ^{E_ℓ/2}, and J_n(−1) = 1. The formula is false outside the
  Blum case: for n = 65, x = 9, z = 16, the congruence is fruitful but the expression gives 0 (codex).
- (J2) **Degree-1 NFS** (L1 = c1x+c0y, L2 = d1x+d0y, so d1L1 ≡ c1L2). The algorithm *imposes* |S| even,
  because the congruence U² ≡ λ^{|S|}W² needs λ^{|S|/2}. Here K0 means the even-weight square dependencies;
  this is a restriction, not a consequence of squareness: codex exhibits an odd-weight square dependency for
  n = 77. Here
  Φ(S) = Σ_{i∈S} ⌊t_i/2⌋ mod 2 with t_i = Ω_T(L1L2(a_i,b_i)) and T = {ℓ : J_n(ℓ) = −1}.
  *Proof:* J_n(L1)J_n(L2) = J_n(c1d1) =: κ for every relation, so t_i ≡ [κ = −1] (mod 2).
  Φ(S) = [κ=−1]·|S|/2 + Σ_S t_i/2, and both cases reduce to Σ⌊t_i/2⌋ because |S| is even.
- **Degeneracy:** if ⌊t_i/2⌋ is even for every relation used, in particular if all prime factors have
  J_n(ℓ) = +1, then every congruence is trivial. The same holds in QS for relations with J_n(x_i) = +1
  (when (2/n) = 1). So fruitfulness lives entirely in T-primes appearing with exponent sum ≡ 2 or 3 (mod 4),
  and in the "which root of n mod ℓ" information that the parity matrix does not record.
  In QS terms, J_n(x) = (−z/x) decomposes by reciprocity into the root choices x ≡ ±√n (mod ℓ) at
  ℓ ≡ 3 (mod 4). Those are the ideals above ℓ in Q(√n).

So, for Blum n, the r16 §3.3 "factor-detecting row" is a **literal computable row** in QS and degree-1
NFS. An algorithm may impose c·s = 1 as an inhomogeneous constraint; every solution is then fruitful,
and a solution exists iff Φ ≢ 0. In NFS of degree ≥ 2, β additionally needs S-unit coordinates and ψ on
units (hidden). It is a row only "in principle".

## 4. The reduction and a hierarchy of sufficient conditions

**Proposition R [proved; codex round 2 and critic: correct with fixes].** Fix (m, f). Let X = X_{n,m,f} be the
finite set of kept pairs. In LV's search (a reading of their pseudocode, not a theorem of theirs), pairs are
i.i.d. uniform in the box and kept iff detected smooth. The first M kept relations are therefore i.i.d.
from the *kept law* on X, whatever the depth cut-off. That law is uniform weighted by the HECM detection
probability, and δ is defined with respect to it (or HECM is boosted). If duplicates are discarded the sampling is without
replacement, and the bound below still holds. Let x: X → F₂^C be the exact square-class vector, with
coordinates fixed on the span generated by all of X *before* sampling: rational exponent parities and
sign, U_S/U_S² coordinates, and the constant feature 1 whenever the algorithm imposes |S| even. Then the
sampled true dependency space is {S : Σ_S x_i = 0}.
*Every algorithmic restriction must be included in x.* Codex's n = 77 example has Φ ≠ 0 on the square
kernel but Φ = 0 on its even-weight subspace. For w ∈ F₂^C let δ_w = Pr_{X}[β ≠ ⟨w,x⟩], and let δ = min_w δ_w. Then
  Pr(Φ|K0 ≡ 0) ≤ Σ_w (1−δ_w)^M ≤ 2^C(1−δ)^M.
*Proof.* (ker)^⊥ = rowspace. On K0, Φ = β restricted, so Φ ≡ 0 iff the vector (β_i) lies in the span of
the coordinate vectors (x_i(·))_i, i.e. iff β_i = ⟨w,x_i⟩ for all i, for some w. Take a union bound over w,
using independence of the M samples. ∎
Hence M ≥ (C ln 2 + ln(1/ε))/δ relations make the *whole sampled kernel* degenerate with probability ≤ ε.
One uniform kernel sample is then fruitful with probability ≥ (1−ε)/2, and t samples fail with probability
≤ ε + 2^{−t}.

**Hierarchy of sufficient conditions (after the critic's correction).** Earlier drafts called Lemma M
"the exact missing lemma"; that was wrong. It is one sufficient condition among several:
- **(F) Conjecture F, batch level; the weakest, and essentially equivalent to "LV's algorithm factors".**
  Under the law of (m, f, batch) that LV's stochastic-deepening search actually outputs, with batches of
  ≥ C·L^{o(1)} relations, Pr[Φ ≢ 0 on K0(batch)] ≥ L_n[1/3,−o(1)]. With codex's restart lemma, Lemma H and
  LV's Thms 2.5/2.6, (F) gives L[1/3,(64/9)^{1/3}+o(1)] factoring of every semiprime. The measure matters.
  Deepening favours (m,f) with large |X|, an event that may have uniform probability only
  2^{−i*} ≥ L[1/3, −(2σ−τ)]. A condition that holds on a uniformly-L^{−o(1)} set of (m,f) could miss those
  entirely. (Codex: "favorable fruitfulness parameters could have low relation density".)
- **(E) Full enumeration (LV's proof route).** Enumerate the whole box for (m,f), at cost L(1/3, 2σ), which
  at the optimum equals (64/9)^{1/3}, and use all of X. Then Φ ≢ 0 ⟺ δ > 0 ⟺ Conj. 7.1 at (m,f) (box
  version). So "Conj. 7.1 on a set of (m,f) of measure L^{−o(1)}, *jointly with* enough relation density" suffices.
- **(M) Lemma M, static.** δ(n,m,f) ≥ L_n[1/3,−o(1)] on a set of (m,f) of measure ≥ L_n[1/3,−o(1)] under
  the search-output (|X|-weighted) law. By Prop R, batches of M ≥ (C ln 2 + O(1))/δ relations (an L^{o(1)}
  enlargement of LV's ≈ C) give (F), since 2^C(1−δ)^M = o(1). Lemma M implicitly needs |X| ≫ C/δ; otherwise
  the x-vectors of the batch are independent and β is trivially linear on it. Lemma M is pointwise *stronger*
  than Conj. 7.1. Its point is that it removes the batch randomness and asks a single static
  equidistribution question about one function on the smooth set. That is the form most amenable to
  analytic attack.
- **(I) Involution, structural.** An x-preserving, β-flipping involution on a part of X of mass ≥ δ′ implies
  δ ≥ δ′/2, hence (M). Dixon's ε-multiplication is the case δ′ = 1.

**Involution criterion [proved; a restatement of Dixon's fibre argument].** Suppose some measure-preserving
involution ι acts on a subset of X of mass ≥ δ′, with x∘ι = x and β∘ι = β + 1. Then δ_w ≥ δ′/2 for every w.
In Dixon, ι(x) = εx with the nontrivial ε = (1,−1) acts on *all* of X (ε = −1 does not flip β), which
recovers the classical proof. In QS, what is needed is a smoothness-preserving map that keeps the parity vector
of x² − n and flips the root class x mod ℓ ↔ −x mod ℓ at an odd number of primes ℓ ≡ 3 (mod 4). x ↦ −x flips
them all at once and leaves β unchanged. No such map is known for QS or NFS.

**Why existing tools fall short [assessment].** LV's analytic inputs are first-order:
- smooth numbers in APs (Harper's Bombieri–Vinogradov for smooth numbers; see also Pascadi 2025);
- effective Chebotarev for quadratic extensions;
- Siegel-zero sparsity.

These control x (parities and quadratic-character values). β is second-order: it is the 2^{c+1}-adic
lift of a character that is already trivial mod 2^c on every relation. In QS it is the choice of ideal above ℓ;
in general it is S-unit coordinates mod 4, weighted by hidden characters. Lemma M needs joint equidistribution
of (x, β) over smooth pairs, uniformly over 2^C linear tests. The union bound in Prop. R absorbs the "uniformly"
only if δ is bounded below for *every* w. So a structural (involution-type) argument is needed rather than a
moment bound. Codex's adversarial view of the routes is in §8.

**Remark (QS: the obstruction is already a short character sum) [assessment].** In QS (Blum,
n ≡ 1 mod 8), β(x) = [J_n(x) = −1]. Even the case w = 0 of Lemma M, or of Conjecture F (β ≢ 0 on X), asks
that the Jacobi symbol mod n be non-constant on the smooth-weighted interval x ∈ [√n − M, √n + M] of length
M = L_n[1/2]. That interval is far below the Burgess range n^{1/4+ε}, where no individual-modulus
character-sum bound is known. GRH bounds for sums starting at 0 do not cover shifted intervals. A proof
would therefore have to use averaging over n or the parameters (large sieve for quadratic characters)
together with the smoothness conditioning. Codex round 2 states the same obstacle in terms of *twisted
smooth correlations*, uniform over all prime-sign twists: the most promising route is averaging over
polynomial parameters, and no available estimate controls all twists after simultaneous-smoothness
conditioning.

## 5. Hidden roots: the degree-1 multiplier (Dixon in disguise) and a counting obstruction

**Theorem L1 [proved; checked in E3].** Let r be uniform in (Z/n)^× and λ = r² mod n. Gauss-reduce the lattice
Λ_λ = {(c,d) : c ≡ λd (mod n)} (determinant n) to e1 = (c1,d1), e2. Put (c0,d0) = e2 + t·e1, and
L1 = c1x + c0y, L2 = d1x + d0y (Res = ±n). The coefficients are O(√n) only *with high probability*.
For λ = 1 the reduced basis is (1,1), ((n+1)/2, (1−n)/2) (codex), and in general |e1|·|e2| ≍ n. One may
reject λ with |e2| > T√n; that rule depends on λ only, so it preserves the argument below. It rejects with
probability O(T^{−2}). Then L1(a,b) ≡ λL2(a,b) (mod n) for all (a,b),
because c0 ≡ λd0 and c1 ≡ λd1. For S with ∏_S L1 = U² and ∏_S L2 = W² (no condition on |S|), set X = U and
Y = r^{|S|}W. Then X² ≡ Y². Every computation except Y depends on r only through λ, so conditional on all of it,
r is uniform on r·μ₂(Z/n). For any output S of odd weight, Pr(fruitful | all else) = 1 − 2^{1−ω(n)} exactly. Here "all else" means the
coins independent of the hidden branch; S must not be selected by inspecting Y. For n = p² this is 0, as it
should be.
*Proof.* r ↦ er with e ∈ μ₂ multiplies Y by e^{|S|} = e, which is a bijection of μ₂(Z/n)/{±1}-classes. ∎
The existence of an odd-weight dependency is a pure parity condition (the all-ones functional ≠ 0 on K0); it
held in 318/320 E3 instances. **Scope:** relation generation (smoothness of two random linear forms of size
√n·box) is not claimed rigorous. So this
gives no improvement over Lenstra–Pomerance. **Correction (critic):** both forms are about √n·A and both must be
smooth, so the product is about n·A² and the heuristic cost is L[1/2, √2] (Dixon size), not L[1/2, 1]. In fact
L1 *is* Dixon: x_i := r·L2(a_i,b_i) satisfies x_i² ≡ L1L2(a_i,b_i), with a single shared hidden root. So Theorem L1
is known in substance. It is kept as the cleanest illustration that "hidden root" means "Dixon fibre".

**Proposition B (counting; scope reduced after codex round 2) [proved].** Let Q = φ(n)/2^ω be the number of
squares. If r is uniform and λ = r², then Pr(λ ≡ c/d with 0 < |c|,|d| ≤ H) ≤ 4H²/Q, so
Pr(H(λ) ≤ n^{1/2−ε}) ≤ n^{−2ε+o(1)}. Equivalently, with probability ≥ ½, max(|c1|,|d1|) ≥ (Q/8)^{1/2}. This
guarantees a large coefficient for one form, not both. *What this does not prove:* lower bounds on form values or on smoothness
cost, or anything about samplers that are uniform only within each root fibre (which is all that
fruitfulness needs). It is a barrier to **globally uniform** hidden roots with small ratio
representatives, i.e. to the naive transfer of Theorem L1 to NFS sizes. (My first wording, "hidden
multipliers force L[1/2]", was REJECTED by codex as overclaimed.) Codex nfsfruit_1 §3.1
proves the companion facts:
- per-relation square multipliers leave Φ unchanged;
- appending one uniform Dixon relation smooth over an L[1/3] factor base costs L[2/3, 2/(3β)].
To keep NFS-sized coefficients one would need a quasi-uniform sampler of r over the μ₂-invariant set
{r : r² ≡ c/d, |c|,|d| ≤ L_n[2/3, ·]}. That set has density L[2/3]/n. Vallée-type samplers reach residues of
size about n^{2/3}, and none reaching n^{o(1)} is known. [assessment]

## 6. Experiments (pre-registered; all exact)

Commands:
```
.venv/bin/python -m src.nfsfruit.run_e1 experiments/configs/nfsfruit_e1.json experiments/raw/nfsfruit_e1.jsonl 9
.venv/bin/python -m src.nfsfruit.analyze_e1 experiments/raw/nfsfruit_e1.jsonl
.venv/bin/python -m src.nfsfruit.run_add experiments/raw/nfsfruit_add.jsonl 9
.venv/bin/python -m src.nfsfruit.analyze_add experiments/raw/nfsfruit_add.jsonl
```
NFS squares are certified by an exact square root in Z[α]: inert-prime Hensel lift, then exact squaring.
The precision cap comes from a proved coefficient bound (`cubic.py`). Factors are used only in
`harness*.py` to classify X/Y mod each prime. For semiprimes the operational gcd test is asserted equal to
the diagnostic class.

**E1** (1120 instances: 7 families × {blum, mixed, both1, tri} × {40, 44} bits × 20 seeds).

| family | median M / cols / dim K0 (44b) | eligible (dim K0 ≥ 10) | rank-deficient | basis vectors fruitful (semiprimes) |
|---|---|---|---|---|
| dixon (control) | 119 / 79 / 44 | 160 | 0 | 2557/5203 = 0.491 |
| trivial (neg. control) | 103 / 63 / 103 | 160 | **160** (as designed) | 0 |
| qs | 117 / 79 / 79 | 153 | 0 | 0.501 |
| lin_base | 669 / 609 / 92 | 160 | 0 | 0.495 |
| lin_rand | 669 / 609 / 91 | 160 | 0 | 0.511 |
| nfs3_base | 824 / 764 / 176 | 160 | 0 | 0.495 |
| nfs3_rand | 821 / 770 / 149 | 101 | 0 | 0.499 |

P0: 0 errors, 0 non-square kernel vectors after certification, 4248/4248 Lemma-H linearity checks.
P1, P2 and P4 as predicted. P5 (binding character constraints, dim ker(parity) − dim K0): median 1 (base)
and 0 (rand), max 3, so 40 characters always sufficed. P3 was underpowered for NFS (0 subsets with k ≤ 4).
This was disclosed and replaced by E1b. nfs3_rand at 40 bits often had too few relations (59 ineligible
instances); this is a yield limitation of the LV-style coefficient growth at toy size, disclosed.

**E1b** (prefix nullity, 5 orders per instance). The deficient rate at dimension k matches 2^{−k} (semiprimes)
and 1−(1−2^{−k})(1−2^{1−k}) (3-prime) for k = 1..6 in **every** family. The largest cluster-robust
|z| = 2.6 (lin_rand, k=3); per-class pooled |z| ≤ 3.3 (lin_rand both1); threshold 4. Example (nfs3_base,
semiprimes): 293/600, 147/600, 74/600, 41/600 against 300, 150, 75, 37.5. No Blum / non-Blum effect.
Under these toy conditions the factor-detecting functional behaves like a uniformly random functional.

**E2** (Blum, 240 instances). Identity c(S) = Φ(S) on every basis vector:
- QS: 13 170 checks;
- lin_base: 47 841 checks;
- lin_rand: 41 682 checks;
- **0 disagreements**.

Degeneracy: on the relations with c_i = 0, rank Φ = 0 in all 240 instances (total subspace dimension 13 195).

**E3** (lin_sq, 320 instances). The hidden-root identity held for 28 029 (basis vector, e ∈ μ₂) checks, with
0 failures. A uniform odd-weight dependency was fruitful in 124/238 semiprimes (0.521; proved value 0.5)
and 62/80 3-prime cases (0.775; proved 0.75). Even-weight dependencies, which follow the heuristic: 0.517 and 0.762.
Odd-weight dependencies existed in 318/320 instances.

## 7. Candidate cards

```text
ID: NF-1 averaged transversality (Conjecture F / quantitative Lemma M) for LV randomized NFS.
Target: every semiprime n (no Blum condition, by Lemma H).
Representation: LV's (m, f = fhat + (x-m)R, pairs, characters).
Operation: uniform dependency in K0; gcd(u - phi(v), n).
Why it exposes a factor: Phi != 0 on K0 => half of K0 is fruitful.
Known lemmas: LV Thms 2.1/2.5/2.6; Lemma H; Prop G; Prop R.
New lemma required: Conjecture F (batch-level, under the search-output law), e.g. via the static Lemma M
  (delta >= L[1/3,-o(1)] on an L^{-o(1)} set of (m,f) under the |X|-weighted law; batches >= C ln2/delta).
Correctness: proved conditional on Conjecture F.
Cost: LV's L[1/3,(64/9)^{1/3}+o(1)].
Randomness: LV's coins plus restarts.
Closest prior: LV Conj 7.1 (needs every (m,f), box version). Conjecture F only needs the batch event
  under the search law. Lemma M is pointwise stronger than 7.1 but static.
Smallest experiment: E1/E1b (done; consistent with a uniformly random functional).
Refutation: a family of (n,m,f) with beta exactly linear in x on X (=> Phi == 0); none seen.
Repair: involution criterion / hidden multipliers (NF-2).
Status: OPEN.

ID: NF-2 hidden-root multiplier (Dixon symmetry via lambda = r^2; Dixon in disguise, x_i = r L2(a_i,b_i)).
Target: every n; two-sided sieve.
Operation: odd-weight dependency, Y = r^{|S|} W.
Proved: fruitful w.p. 1 - 2^{1-omega} given an odd-weight dependency.
Cost: coefficients ~ sqrt(n) w.h.p.; heuristic L[1/2,sqrt 2] (Dixon size); no gain over Lenstra-Pomerance.
Barrier: counting only (Prop B): a globally uniform hidden root has small ratio representatives with
  probability <= 4H^2/Q. Fibre-wise-uniform samplers at NFS size are not excluded.
Status: proved (degree 1); refuted for this model at L[1/3] (globally-uniform-r version only).

ID: NF-3 Jacobi row (Blum, QS / degree-1).
Proved identity; turns fruitfulness into rank([A; c]) = rank A + 1 with c computable from n.
Does not prove the rank increase. Status: proved; no complexity gain by itself.
```

## 8. GPT-6.1-sol rounds
- **nfsfruit_1 (blind).** It independently proved Lemma H and the rowspace criterion; gave the Dixon
  root-fibre symmetry; showed per-relation square multipliers leave Φ unchanged and that a single Dixon
  relation over an L[1/3] base costs L[2/3]; gave the Fermat close-prime family as a trivial "special family";
  stated an averaged transversality target equivalent in spirit to Lemma M; and found no 2018–2026 fruitfulness
  result. Its claim "ordinary square multipliers cannot inject the missing randomness" covers per-relation
  multipliers. Theorem L1 escapes it by hiding the root in the *ratio of leading coefficients*, which is
  possible only at size √n.
- **nfsfruit_2 (adversarial)**, verdict "Claims 1–4 sound with fixes; Claim 5 wrong as a general barrier".
  All fixes are applied above:
  - coprimality and distinct-prime hypotheses; a non-canonical identification; the distribution of
    uniform-K0 sampling;
  - the cleaner QS formula and the non-Blum counterexample (n = 65);
  - "|S| even" is an imposed restriction (n = 77 counterexample), so the constant feature must be in x;
  - the union bound is legitimate with coordinates fixed on all of X, and concerns degeneracy of the
    whole kernel; the full-enumeration repair preserves (64/9)^{1/3};
  - Theorem L1's O(√n) holds only w.h.p. (λ = 1 counterexample), with a λ-only rejection rule;
  - Proposition B downgraded to a counting statement;
  - prior art: Seysen 1987 p. 777 for the subgroup mechanism, BLP §8 Lemma 8.2 for square detection, and
    Lenstra–Pomerance §2/§10 for ambiguous forms.
  It found no published statement of the Bockstein formula, of the hidden-root two-form construction, or of
  Conjecture F. Its scratch script `codex/scratch/nfsfruit_2/referee_checks.py` (not committed; git-ignored)
  passed 1008 further kernel-vector tests of Prop. G, including unequal local 2-parts (p ≡ 1 mod 4).

## 9. Review status
- **Codex nfsfruit_2 (adversarial):**
  - Claims 1–4: correct with fixes, applied.
  - Claim 5: wrong as a general barrier; downgraded to counting.
- **Independent read-only Opus critic:**
  - Lemma H: ACCEPT WITH FIXES (monic normalisation; which K0).
  - Prop G: ACCEPT. The critic also ran its own exact non-Blum degree-1 check with c = 4 and 5: 300/300
    agreements per instance, and 2^c | t_i asserted for every relation.
  - J1/J2: ACCEPT. The Verdict's degeneracy wording was REJECTED (false for QS) and is fixed.
  - Prop R inequality and the involution criterion: ACCEPT.
  - "Lemma M is the exact missing lemma": REJECTED (it is only sufficient). Now replaced by the hierarchy
    F/E/M/I with the search-output measure.
  - Theorem L1: ACCEPT WITH FIXES. Its cost is L[1/2, √2], it is Dixon in disguise, and the O(√n) bound holds only w.h.p.
  - Prop B: ACCEPT WITH FIXES (constants; one form).
  - LV reading: the batch structure (fresh batches of ≈ C relations) was corrected. P_S is defined by height.
    The exceptional "set".
  - All fixes are applied in this version.
- **Surviving claims labelled [proved]:** Lemma H, Prop G, Cor J (J1, J2, degeneracy in its stated scope),
  Prop R, the involution criterion, Theorem L1 (known in substance), Prop B (counting), and the implications
  (F) ⇒ L[1/3] factoring and (M) ⇒ (F) (modulo LV's non-monic bookkeeping, which is not written out).
- Folklore / known: Lemma H (cf. Seysen 1987 p. 777), the involution criterion, the Blum–Jacobi fact behind
  Cor J, and L1. Prop G in this packaging: no prior formulation found by us or by either reviewer.
