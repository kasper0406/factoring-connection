# prepread: high-degree algebraic preprocessing, characters, and dense short SLPs (phase 18)

Agent `prepread`, 2026-10-20.
- Pre-registration: `notes/prepread_prereg.md` (commit `fdb1328`; Amendment A1 in the next commit, made before
  any experiment ran).
- Code: `src/prepread/{slp,spectra,dlred,compsearch,run_experiments}.py`, `src/prepread/compsearch.c`.
  Tests: `tests/test_prepread.py`.
- Raw data: `experiments/raw/prepread_e{1,1b,2,3,4,5}.jsonl`.
- GPT-6.1-sol: blind round `codex/prepread_1.{md,out}` (full report in `codex/scratch/prepread_1/report.md`,
  git-ignored); adversarial round `codex/prepread_2.{md,out}`.
- Two read-only Opus critics (proofs; literature). Adjudication in §9.
- Labels: **KNOWN** (with reference), **PROVED** (proof here; "no prior formulation found by us" is the strongest
  novelty claim), **CONDITIONAL**, **HEURISTIC**, **CONJECTURED**, **experimentally supported**,
  **REFUTED for this model**, **UNVERIFIED**.

Notation. N = pq balanced, p, q ∈ (X, 2X]. e(x) = exp(2πix). For Φ ∈ Z[y] (or an SLP), a prime p, a multiplier
c ∈ F_p and a multiplicative character χ of F_p^*,

  A_p(χ, c) = (1/(p−1)) Σ_{y∈F_p^*} χ(y)·e(cΦ(y)/p).

For the composite observable Ψ(y) = e(hΦ(y)/N), CRT gives Ψ̂(χ_pχ_q) = A_p(χ̄_p, h q̄)·A_q(χ̄_q, h p̄), with
q̄ = q^{−1} mod p and p̄ = p^{−1} mod q.

---

## 0. Verdict (plain language)

**No mechanism. The two remainders are not one question: they separate exactly at the order of the character.**

- **O1 (high-degree preprocessing, additive readout) is at least as hard as discrete logarithms in prime fields.**
  - Theorem O1-DL (§2.4, PROVED, routine): an efficiently computable phase that correlates with a character of
    order d at a prime p yields discrete logarithms modulo d at p.
  - Hence a working O1 observable (a 1/poly-heavy peak of order ≥ r/poly, as periodfind's Theorem T requires)
    would give polynomial-time discrete logarithms in F_p^* for a positive proportion of primes p in the range,
    on top of factoring.
  - Prime-field DL is not known to follow from factoring. So O1 is a strictly stronger target than factoring,
    and it is barred by any hardness of prime-field DL (CONDITIONAL barrier, §2.5).
  - Unconditionally, Weil bounds the correlation only below cyclic exponent width √p (§2.1). That is SLP size
    ½·log₂ p, so the unconditional barrier is vacuous for short SLPs, as expected (proving more would give circuit
    lower bounds for DL).
  - A size-and-order barrier with no further hypothesis is false. Pohlig–Hellman is itself a short SLP (§2.3):
    an N-independent polylog program has correlation ≈ 1 with a character of order 2^e at every prime with
    p − 1 = 2^e·(odd B-smooth). That is the p − 1 regime (order annihilation), so the limit bites exactly at
    the largest prime factor of p − 1.
- **The dense-torus question has no per-prime obstruction at all.** At a known prime, root density 1/2 is
  trivial (x^{(p−1)/2} − 1). The difficulty is purely uniformity in p.
  - **Galois theory alone cannot give a barrier** (Prop GAL-0, PROVED). For every interval there are polynomials
    with trivial, cyclic and symmetric Galois groups whose root density is ≈ 1/3 at every prime of the interval.
    Fixed-point proportions of the group are the wrong statistic (§3.1). Effective Chebotarev is silent once the
    degree exceeds ≈ √X/log²X, and dense roots need degree ≥ μX (§3.2). Abelian versus non-abelian is also the
    wrong dividing line: ECM's division-polynomial fields have GL₂ images (§3.4).
  - **Galois theory plus the arithmetic of the roots does give a restricted inverse theorem** (Theorem MIN, §4.4,
    PROVED; no prior formulation found by us). A composition of monic integer polynomials from a fixed finite set
    (e.g. x² + c, |c| ≤ 2), with degree at most 2^Δ·X, has root density < μ at a positive proportion of the
    primes of (X, 2X] for all large X. The proof combines Chebotarev for the finitely many outer prefixes,
    Burnside and Northcott (roots of such compositions have bounded house). So dense compositions need excess
    depth Δ → ∞. The polylogarithmic-depth case stays open (Conjecture QC, §4.6).
  - Together with codex's fixed-prefix and outer-constant obstructions (§4.2–4.3), these close every
    composition class we could state except "intersective words of growing excess depth".
- **Where the limit bites** (§5): for O1, at DL in the subgroup of order P⁺(p − 1); for the dense-torus question,
  nowhere per prime, only through uniformity in p. There, Prop BB (slpcirc) makes a general barrier
  permanent-hard, and Prop GAL-0 shows that no group-theoretic shortcut exists.

**Subjective assessment.** Probability that either remainder yields sub-L[1/3] factoring: < 1%. O1 is now the
less plausible of the two, since it would also break prime-field DL.

---

## 1. The common frame: one two-parameter family of sums

**Prop UF (the two faces). PROVED, trivial.** For Φ ∈ Z[y] and a prime p,

  #{y ∈ F_p^* : Φ(y) ≡ 0}/(p − 1) = E_{c∈F_p} A_p(χ₀, c),   χ₀ = trivial character.

*Proof.* E_{c∈F_p} e(cv/p) = [v ≡ 0]. ∎

So both remainders ask for large values of the same family A_p(χ, c):
- **slpcirc (model U, W_U):** the trivial character, averaged over the readout multiplier c, must lie in
  [μ, 1 − μ] at (almost) every prime. The readout is a gcd, which is why the trivial character is useful.
- **periodfind O1:** a single (χ, c) with |A_p(χ, c)| ≥ τ for χ of order ≥ p/poly, at both p and q (Corollary
  O1-N below). The readout is a Fourier peak, which is why only large-order characters are useful (Theorem T).
- **Small orders in between** reveal only [d | r]-type residuosity bits (periodfind §1.2, §2) and need no new
  analysis.

The order of χ is the dial. Order 1 is easy at each prime and hard only through uniformity (§3–4). Large prime
order is DL-hard at each prime (§2.4).

---

## 2. O1: high-degree preprocessing with an additive readout

### 2.1 The Weil threshold, exactly. KNOWN (Weil; Castro–Moreno; Fu–Wan) + PROVED (width reduction)

*Bound.* Let χ be nontrivial and f ∈ F_p[y, y^{−1}] a Laurent polynomial with exponents in [−a, b], a, b ≥ 0,
a + b ≥ 1, a, b < p. Then |Σ_{y∈F_p^*} χ(y) e(f(y)/p)| ≤ (a + b)√p.
- This is the Weil bound for the rank-one sheaf L_χ ⊗ L_ψ(f) on G_m. Its Swan conductors are a at 0 and b at ∞
  (both prime to p, so no Artin–Schreier degeneracy), and H²_c vanishes.
- Reference: Castro–Moreno, Proc. AMS 128 (2000) 2529–2537; the conductor form for rational f = A/B is Fu–Wan,
  Quart. J. Math. 65 (2014) 1195–1211, Thm 4.6 of the preprint (codex; our §1 reading of the theorem number is
  UNVERIFIED pending the literature critic). Codex quotes the bound [b + max(a − b, 0) + t]√p for f = A/B.

*Effective degree (PROVED).* Only the function y ↦ Φ(y) on F_p^* matters, i.e. Φ mod (y^{p−1} − 1). Lift its
nonzero exponents to integers in a window [−a, b]. The best lift minimises a + b. Call the minimum the **cyclic
width** w_p(Φ): the length of the shortest arc of Z/(p − 1) containing the support and 0 (`slp.weil_width`).
Then

  |A_p(χ, c)| ≤ w_p(Φ)·√p/(p − 1)  for every nontrivial χ and c ≠ 0.

- **Threshold.** The bound is trivial exactly when w_p(Φ) ≥ (p − 1)/√p, i.e. **w ≈ √p**. A τ-correlation needs
  w_p(Φ) ≥ τ(p − 1)/√p.
- For the composite observable this gives |Ψ̂(χ)| ≤ w_p w_q/√N, the D⁴N/r² line of periodfind PF case 2.
- **Short SLPs pass it at once.** w_p ≤ deg Φ ≤ 2^s, so Weil binds only when s < ½·log₂ p − log₂(1/τ). Generic
  SLPs have w_p ≈ p.
- Conductor is the right complexity (codex). y^{p−2} = 1/y has w = 1. A rational representative can have a much
  smaller conductor than the polynomial one. Finding it need not preserve SLP size.

*Compositions do not help Deligne either (HEURISTIC remark).* Pushing L_χ forward through a degree-2 layer gives
Σ_x χ(x)e(G(x² + c)/p) = Σ_u [λ(u − c) + λη(u − c)] e(G(u)/p), where λ² = χ and η is the Legendre symbol. So a
k-layer composition is a sum against a trace function of rank 2^k, and the conductor grows like the degree. A
sheaf bound is again vacuous above √p.

### 2.2 Monomials. PROVED from a KNOWN theorem (codex)

Let d = gcd(m, p − 1) and H = (F_p^*)^m, with L = |H|.
- If χ is nontrivial on ker(y ↦ y^m), then A_p = 0.
- Otherwise A_p(χ, c) = (d/(p − 1))·Σ_{z∈H} θ(z)e(cz/p), with ord χ ≤ L. Expanding the square and using
  Bourgain–Glibichuk–Konyagin (J. London Math. Soc. 73 (2006) 380–398, Thm 6) gives

    |A_p(χ, c)| ≤ (L^{−1} + p^{−η(ε)})^{1/2}  when ord χ ≥ p^ε.

So monomial preprocessing of any degree has power cancellation against characters of polynomially large order.
Concentration at small order L is the order-annihilation regime (periodfind §1.3f).

### 2.3 Prop PH (Pohlig–Hellman is a short SLP). PROVED

**The construction.**
- Fix e ≥ 1 and an odd m (e.g. the odd part of lcm(1..B)). The program has two inputs y, w, ring operations and
  halving (a division by the unit 2, gcd-checked as in slpcirc §1.4).
  - z = y^m, ζ = w^m, ζ⁻¹ := ζ^{2^e − 1}.
  - For i = 0..e−1: s_i = z_i^{2^{e−1−i}} ∈ {±1}; b_i = (1 − s_i)/2;
    z_{i+1} = z_i·((1 + s_i) + (1 − s_i)ζ_i^{−1})/2; ζ_i^{−1} ← (ζ_i^{−1})².
  - Output L = Σ b_i 2^i.
  - Size O(log m + e²). Every exponent and constant is a fixed integer.
- **Claim.** Let p − 1 = 2^e·u with u | m odd, and let w be a quadratic non-residue (probability 1/2). Then ζ
  generates μ_{2^e}, and L(y) = log_ζ(y^m) ∈ [0, 2^e).
- **Correlation.** y ↦ e(L(y)/2^e) is a character χ of order 2^e. With c = round(p/2^e), |A_p(χ, c)| ≥ sinc-type
  ≥ 0.8 for u ≥ 3 (unit test: 0.99 at p = 3·2^{12} + 1, 274 operations).

*Proof.* y ↦ y^m maps F_p^* onto μ_{2^e} because u | m. Peeling bits works because s_i = (−1)^{b_i} when
z_i ∈ ⟨ζ^{2^i}⟩. Then e(cL/p − L/2^e) has phase |L(c/p − 2^{−e})| ≤ 2^e/(2p) < 1/(2u). ∎

**Reading.**
- Codex's prime-specific version: y^{(p−1)/R} with a field constant ζ, plus a Boolean rounding circuit. It gives
  correlation ≥ cos(π/p) with order R = 2^b. Via Xylouris (Linnik constant 5.18), it has R ≥ p^{0.19} infinitely
  often.
- So **"SLP size s and character order ≥ R ⇒ small correlation" is false**, even for N-independent programs.
- The exceptions are exactly the primes where the χ-order part of p − 1 is annihilated by a fixed smooth
  exponent. That is order annihilation (nonsmooth Prop OA, slpcirc row 2).
- The next theorem says this is forced: outside smooth orders, correlation is a DL algorithm.

### 2.4 Theorem O1-DL (correlation ⇒ discrete logarithm). PROVED (routine; the classical analogue of Shor's DL step)

**Statement.** Let p be prime and f: F_p^* → {|z| ≤ 1} computable in time T. Suppose |⟨f, χ⟩| ≥ τ for a
character χ of order d. Given p and a generator g, the following runs in expected time
O((T + log² p)·τ^{−2}·log p·log(log p/δ)) plus one sparse-Fourier search of cost poly(log p, 1/τ)·T:
1. find the frequency j with χ(g) = e(j/(p − 1)) (or one of the ≤ 1/τ² frequencies at least as heavy);
2. for any y, compute log_g y mod d with error probability ≤ δ.

**Proof.**
1. *Frequency.* F(t) = f(g^t) on Z/(p − 1) has F̂(j) = ⟨f, χ⟩ exactly, since g is a generator so restriction
   classes are singletons (Theorem T with i = 1). Find j by sparse Fourier transform over Z_{p−1} with the known
   modulus p − 1 (Mansour 1995; Akavia 2009 Thm 4; cost poly(log p, 1/τ) queries). Its list has ≤ 1/τ² entries;
   try each.
2. *One phase.* For any z, E_s f(z g^s)e(−js/(p − 1)) = F̂(j)·e(j·log_g z/(p − 1)). Substitute t = s + log_g z.
   Sampling m = O(τ^{−2} log(1/δ')) values of s estimates θ(z) = j·log_g z/(p − 1) mod 1 to within 1/16 with
   probability 1 − δ'. F̂(j) itself is estimated once, the same way, at z = 1.
3. *Doubling.* Measure θ_k ≈ frac(2^k θ) at z = y^{2^k} for k = 0..K, K = ⌈log₂(p − 1)⌉ + 3. Going backwards, θ_k
   is the one of θ_{k+1}/2 and (θ_{k+1} + 1)/2 that is nearer the measurement. The two candidates are 1/2 apart,
   so a measurement error < 1/8 picks the right one, and the error in θ halves each step. After K steps θ is
   known to < 1/(2(p − 1)); round to get j·log_g y mod (p − 1). This gives log_g y modulo
   (p − 1)/gcd(j, p − 1) = d.
4. *Union bound* over the K + 1 estimates, δ' = δ/(K + 1). ∎

Notes.
- No generator? A random element is a generator with probability φ(p − 1)/(p − 1) ≫ 1/log log p. Success is
  checkable once the full log is found. The cofactor gcd(j, p − 1) is handled by BSGS in its square root.
- Code: `src/prepread/dlred.py`. Tests: `test_dlog_reduction_planted`, `test_dlog_reduction_noncoprime_frequency`.
- Codex proved the same per-prime statement in its blind round (with explicit access assumptions) and labelled
  the resulting barrier CONJECTURED. Our addition is the N-level corollary below, which removes the access
  assumptions in the factoring setting.

**Corollary O1-N (an O1 factoring observable is a prime-field DL algorithm). PROVED.**
- *Hypothesis.* Fix N-independent programs Φ_X of size s and a computable multiplier map N ↦ h(N), possibly
  randomised. Suppose that for a fraction ≥ μ of balanced N = pq with p, q ∈ (X, 2X], the observable
  Ψ(y) = e(h(N)Φ_X(y)/N) has, for some a of index i = φ(N)/ord(a), a τ-heavy peak of order ≥ r/κ.
- *Conclusion.* For a fraction ≥ μ/2 of primes p ∈ (X, 2X], discrete logarithms in F_p^* are computable in
  expected time poly(s, log X, 1/τ, i, κ, 1/μ).

*Proof.*
1. *Correlation.* Theorem T (T2) gives χ = χ_pχ_q with |Ψ̂(χ)| ≥ √τ/i and ord χ ≥ r/κ. By CRT,
   Ψ̂(χ) = A_p(χ̄_p, h q̄)·A_q(χ̄_q, h p̄), with both factors ≤ 1. So |A_p(χ̄_p, c_p)| ≥ √τ/i, where c_p = h q̄ mod p.
2. *Order.* ord χ = lcm(ord χ_p, ord χ_q) ≤ ord χ_p·(q − 1) and r = (p − 1)(q − 1)/i. So
   ord χ_p ≥ (p − 1)/(κi). The same holds at q.
3. *Averaging.* Call a pair (p, q) good if it satisfies the hypothesis. Then for at least μ/2 of the primes p, at
   least μ/2 of the partners q are good.
4. *Reduction at such a p.* Draw a random prime q′ ∈ (X, 2X] (a primality test is used), set c = h(pq′)·q′^{−1}
   mod p, and run Theorem O1-DL on f = e(cΦ_X/p). It succeeds with probability ≥ μ/2 per draw, and success is
   checkable.
5. *Remaining log.* The log is known modulo a divisor of p − 1 of index ≤ κi. BSGS finishes in poly time. ∎

### 2.5 The barrier, stated three ways

1. **CONDITIONAL (from O1-DL).** Let T_DL(p, ℓ) be the cost of DL in the order-ℓ subgroup of F_p^*. A size-s
   Φ whose phase correlates τ with a character of order divisible by ℓ satisfies

     τ ≤ Õ( (s·log² p / T_DL(p, ℓ))^{1/2} ).

   With ℓ = P⁺(p − 1) (typically ≥ p^{0.6}) and the conjectured L_p[1/3] hardness of prime-field DL,
   τ ≤ L_p[1/3, c]^{−1/2}. That is far below the 1/poly that Theorem T needs.
2. **Sharpness in the order** (Prop PH). Without excluding smooth orders the statement is false. The bound is
   about the largest prime factor of the order, not the order itself.
3. **Unconditional** (Weil, §2.1; BGK for monomials, §2.2). τ ≤ w_p√p/(p − 1). This is vacuous for s ≥ ½·log₂ p.
   - An unconditional bound τ ≤ poly(s)·p^{−c} for all size-s SLPs would give a superpolylogarithmic SLP lower
     bound for weak DL predictors.
   - Codex reports that Coppersmith–Shparlinski (J. Cryptology 13 (2000) 339–360) leave SLP complexity of DL
     explicitly open. Known lower bounds bound degree, sparsity or linear complexity instead: Coppersmith–
     Shparlinski Thms 1, 3; Lange–Winterhof, Acta Arith. 101 (2002); Meidl–Winterhof, IEEE-IT 47 (2001).
   - We expect such a bound to be as hard as the other circuit lower bounds this project keeps meeting (Prop BB).

**What remains open in O1.** Correlation with large-prime-order characters by short N-independent SLPs, without a
DL algorithm hiding inside. O1-DL says such an SLP *is* a DL algorithm (of a restricted, algebraic kind), so the
open question is now "is there a polylog algebraic-circuit DL predictor for F_p^*". We know of no reason to
expect one.

---

## 3. The Galois framing of the dense-torus question

### 3.1 Dictionary, and the right statistic. KNOWN / PROVED (routine)

- For F ∈ Z[x] with lc(F)·disc(rad F) ≢ 0 mod p: Z_p := #{x ∈ F_p : F(x) = 0} = fix(Frob_p) on the distinct roots.
  This holds also when deg F ≥ p (codex).
- For G acting with t orbits: E_G fix = t (Burnside), and Pr_G[fix ≥ M] ≤ t/M (Markov).
- So model U asks that a **specific** element, Frob_p, fix ≥ μp roots for every p in (X, 2X], where
  deg F ≥ μX and t is typically tiny.
- That is a statement about the Frobenius *sequence at small primes*, not about the group.
- The fixed-point proportion FPP(G) = Pr[fix > 0] only says whether there is a root at all. For transitive
  actions the mean is always 1, and FPP ∈ [1/D, 1 − e^{−1} + o(1)] for D-point actions (natural S_D gives
  1 − e^{−1}).
- So the user's "fixed-point proportion in [1/polylog, 1 − 1/polylog]" is plentiful and irrelevant to density.
  The relevant event is "Frob_p ∈ {σ : fix(σ) ≥ μp}", which has Haar measure ≤ t/(μp).

### 3.2 Prop CS (Chebotarev silence, quantified). PROVED (from KNOWN effective Chebotarev)

Assume GRH, and let F have irreducible factors F_1..F_t with fields K_i = Q[x]/F_i of degree n_i. By the
GRH prime ideal theorem (Lagarias–Odlyzko 1977 Thm 1.1; Serre 1981 §2.4),

  Σ_{p∈(X,2X], p∤disc} Z_p(F) = t·(π(2X) − π(X)) + O( Σ_i √X·(log d_{K_i} + n_i log X) ).

- The error is below the main term only if n_i ≲ √X/log² X and log d_{K_i} ≲ √X/log X. Dense roots need
  Σ n_i = deg F ≥ μX.
- So in the whole regime that model U needs, effective Chebotarev is silent by a factor ≥ √X, even on average
  over p.
- Unconditionally it is far worse: Lagarias–Odlyzko Thm 1.3 needs X ≥ exp(10 n_K (log d_K)²) (codex).
- This is the root-counting analogue of the Weil threshold of §2.1. Both faces live above the √-scale where
  every RH-type tool is silent.

### 3.3 Prop GAL-0 (no Galois-only barrier). PROVED (CRT + standard criteria; codex found the S_D case independently)

Fix X and D with X/3 ≤ D < X. Then there are monic F ∈ Z[x] of degree D with Z_p(F) = D, i.e.
δ_p = D/p ∈ (1/6, 1), at every prime p ∈ (X, 2X], whose Galois group is
- (a) trivial: F = ∏_{j≤D}(x − j);
- (b) S_D: F = ∏_{j≤D}(x − j) + M·h(x), with M = ∏_{X<p≤2X} p and h of degree < D chosen by CRT at three auxiliary
  primes ℓ₁, ℓ₂, ℓ₃ ∤ M. Make F Eisenstein at ℓ₁, of type (D − 1, 1) mod ℓ₂, and of type (2, 1, …, 1) with
  distinct factors mod ℓ₃. Transitivity, a (D − 1)-cycle and a transposition give S_D (van der Waerden;
  Dedekind's theorem);
- (c) cyclic of order D: take a prime ℓ ≡ 1 (mod D) such that every p ∈ (X, 2X] is a D-th power residue mod ℓ.
  Such ℓ exist by Chebotarev in Q(ζ_D, p^{1/D} : p ∈ (X, 2X]). Let K ⊂ Q(ζ_ℓ) be the degree-D subfield. Every
  p ∈ (X, 2X] splits completely in K.
  - Choose θ ∈ O_K by CRT so that θ has D distinct residues in the D degree-one primes above each p (possible
    since D < p). θ then generates K, and its minimal polynomial F has D distinct roots mod each p.

To get density ≈ 1/3 rather than ≈ 1, replace D by ⌊X/3⌋ in (a)–(c). Codex's variant prescribes ⌊p/2⌋ linear
factors at each p, giving density ≈ 1/2 with group S_D.

**Consequence.** Every Galois group we tried (trivial, cyclic, symmetric) admits W_U-type root counts on any
interval. So a barrier must use the SLP. And any barrier that covers case (a) implies that the Pochhammer product
has no polylog SLP, which by Bürgisser implies superpolynomial τ(Per) (slpcirc Prop BB(b)). **A barrier of the
form "short SLP + Galois group in class 𝒢 ⇒ no W_U" is permanent-hard whenever 𝒢 contains the trivial group.**

### 3.4 Abelian, non-abelian, and "pinned Frobenius"

- **Prop AB (abelian ⇒ torus torsion). PROVED (Kronecker–Weber).** If F has abelian splitting field of conductor
  m, its roots lie in Q(ζ_m) = Q(G_m[m]), and
  Z_p(F) = Σ_χ n_χ χ(p) for p ∤ m·disc. The sum runs over Dirichlet characters mod m, with n_χ ≥ 0 integers
  (the multiplicity of χ in the permutation character), Σ n_χ = deg F and n_{χ₀} = t.
  - So Z_p depends only on p mod m. Each irreducible factor contributes either its whole degree (p splits
    completely in its field, i.e. p lies in a subgroup H_i ≤ (Z/m)^*) or nothing.
  - p − 1, Williams p + 1 (roots ζ + ζ^{−1}) and Bach–Shallit Φ_k (Gauss periods of order-k subgroups) are all of
    this form.
  - In this precise sense "abelian = G_m-torsion (torus) methods". Prop AB alone gives no size barrier: case
    (c) of GAL-0 is abelian.
- **Non-abelian does not mean "not group order".** Division polynomials ψ_M of a non-CM curve have Galois group
  ⊆ GL₂(Z/M), non-abelian. Their root counts are governed by #E(F_p)[M], i.e. by the characteristic polynomial
  x² − a_p x + p of Frobenius. That is ECM, covered by nonsmooth Prop OA. So the dividing line is not
  abelian/non-abelian.
- **CONJECTURED (Frobenius pinning).** The known dense-root mechanisms in the silent regime are those where Frob_p
  acts through data of size p^{O(1)} that are pinned by p itself:
  - the integer p acting on G_m-torsion;
  - a Frobenius π with |π| = √p acting on torsion of an abelian variety;
  - the identity (all roots rational: Pochhammer, constants).
  Short SLPs of huge degree with pinned Frobenius are torsion polynomials of commutative group schemes (multiples
  [M] have SLP size O(log M) through the group law, Ritt–Eremenko for commuting families), plus the explicit-root
  route. We conjecture that a polylog N-independent SLP with W_U must fall in one of these three classes. We have
  no proof. Prop BB shows that even the third class cannot be excluded without a permanent lower bound.

---

## 4. Quadratic compositions: toward the inverse theorem

Setting: g_i(x) = x² + c_i with c_i ∈ C, a fixed finite set of integers (C = {−2, …, 2} in the experiments), and
F = g_k ∘ ⋯ ∘ g_1. At a prime p, T_k = {0} and T_{j−1} = g_j^{−1}(T_j), so Z_p(F) = |T_0| and

  |T_{j−1}| = Σ_{u∈T_j} (1 + η(u − c_j)) = |T_j|(1 + β_j),   β_j := mean of η(u − c_j) over T_j,

where η is the Legendre symbol (η(0) = 0, one preimage). In Galois terms β_j is the "non-swap" frequency of
Frob_p at the fixed vertices of level j of the preimage tree.

### 4.1 Lemma DEF (deficit). PROVED, elementary

If Z_p(F) ≥ μp with k = log₂(μp) + Δ, then Σ_j (1 − β_j) ≤ 2 ln 2·Δ ≈ 1.39Δ.

*Proof.* Z_p = ∏(1 + β_j), so Σ_j (1 − log₂(1 + β_j)) = k − log₂ Z_p ≤ Δ. By concavity
1 − log₂(1 + β) ≥ (1 − β)/(2 ln 2) on [−1, 1]. ∎

So near-minimal-depth dense compositions need almost every level of the preimage tree to be almost entirely
"square-shifted". The per-prime statement cannot be closed by character-sum methods. The sets T_j are
0-dimensional (roots of a univariate polynomial), and there is no square-root-cancellation theory for character
sums over the roots of a fixed-degree-≥-√p polynomial at a single prime.

### 4.2 Prop OUT (outer constant). PROVED (codex), trivial

A root of F mod p gives a root of g_k. If −c_k is not a square in Z, then g_k has no root at the primes with
(−c_k/p) = −1. That is half the primes asymptotically, uniformly over the finitely many c_k ∈ C. So W_U forces the
outermost map to have a rational root (c_k = −s²).

### 4.3 Prop PFX (fixed-prefix towers). PROVED (codex)

Let F_k = P ∘ R_k with a **fixed** outer prefix P, and suppose FPP(Gal(P_m)) → 0 along a fixed tower of prefixes
P_m (e.g. iterates of x² + c with full wreath image: Stoll 1992, as stated in Jones 2013 Thm 2.3; x² + 1
qualifies). Then for any depths k(X) → ∞, the proportion of primes of (X, 2X] where F has any root tends to 0.

*Proof.* A root at depth k maps to a root of P_m. Chebotarev at fixed depth bounds the proportion by FPP(G_m) +
o(1). Let m → ∞. ∎

FPP for the full binary wreath product satisfies q_{n+1} = q_n − q_n²/2, so q_n ~ 2/n (Odoni 1985). For bounded
degree d it is O_d(1/n) (codex). Our E3 tests this Haar model in the silent regime.

### 4.4 Theorem MIN (no dense compositions of near-minimal degree). PROVED; no prior formulation found by us

**Statement.** Let 𝒢 be a finite set of monic polynomials in Z[x] of degree ≥ 2, and let μ ∈ (0, 1] and Δ ≥ 0.
There are κ = κ(𝒢, μ, Δ) > 0 and X₀ such that, for every X ≥ X₀ and every composition F = g_k ∘ ⋯ ∘ g_1 with
g_i ∈ 𝒢 and deg F ≤ 2^Δ·X, the root density satisfies Z_p(F) < μp for at least a κ-fraction of the primes
p ∈ (X, 2X].

In particular W_U(μ) (density ≥ μ at all but one prime) fails. Since W_U needs deg F ≥ μX, **a dense
composition needs degree ≥ X·2^{Δ(X)} with Δ(X) → ∞**.

**Proof.**
1. *Invariant disc.* Since every g ∈ 𝒢 is monic of degree ≥ 2 with bounded coefficients, there is R = R(𝒢) ≥ 1
   such that |β| ≤ R and g(α) = β imply |α| ≤ R. (For x² + c with |c| ≤ H take R = (1 + √(1 + 4H))/2.
   In general, for |α| > 1, |α|^d ≤ R + H·d·|α|^{d−1} gives |α| ≤ Hd + R/|α|^{d−1}, hence |α| ≤ Hd + √R ≤ R
   once R ≥ (2Hd + 1)², where H bounds the coefficients and d the degrees.) By induction,
   every complex root of every composition over 𝒢 has all conjugates in |z| ≤ R.
2. *Northcott count.* Let N(d) be the number of monic irreducible integer polynomials of degree ≤ d all of whose
   roots lie in |z| ≤ R. N(d) is finite: the coefficients are bounded by binom(d, i)R^i.
3. *Prefix transfer.* Write F = P_m ∘ Q with P_m the outer m layers. If Z_p(F) ≥ μp > μX then, since a layer of
   degree e multiplies preimage counts by at most e,
   #{y ∈ F_p : P_m(y) = 0} ≥ μX/deg Q = μX·deg(P_m)/deg F ≥ μ2^{−Δ}·deg P_m.
4. *Choice of m.* Pick d with (d + 1) > 2^{Δ+1}/μ. Then pick m with 2^m > 2^{Δ+1}N(d)/μ. Every prefix P of m
   layers has deg P ≥ 2^m. Its number t_P of distinct irreducible factors satisfies
   t_P ≤ N(d) + deg P/(d + 1) < μ2^{−Δ} deg P.
5. *Chebotarev on finitely many prefixes.* For each of the finitely many m-layer words P over 𝒢, let G_P act on
   the distinct roots of P. By Burnside, t_P = E_{G_P} fix < μ2^{−Δ} deg P. So the set
   S_P = {σ : fix(σ) ≥ μ2^{−Δ} deg P} has |S_P| < |G_P|.
   By Chebotarev (dyadic form, fixed field), the primes p ∈ (X, 2X] with Frob_p ∉ S_P have proportion
   → 1 − |S_P|/|G_P| ≥ 1/|G_P|. At such p (and p ∤ disc(rad P)), P has < μ2^{−Δ} deg P roots mod p, so by step 3
   Z_p(F) < μp.
6. Take κ < min_P 1/|G_P| over the finitely many prefixes, and X₀ large enough for all of them. Compositions with
   fewer than m layers have degree ≤ (max deg)^m < μX for X large, so they are trivially excluded. ∎

**Remarks.**
- The ingredients are Chebotarev for the finitely many outer prefixes of fixed depth (the uniformity comes from
  finiteness, as in codex's PFX) and Burnside. The step that is *not* group theory is Northcott. Roots of
  compositions over 𝒢 have bounded house, so a prefix has only boundedly many low-degree irreducible factors, and
  most of its Frobenius elements must move most roots. This is exactly what Prop GAL-0 says Galois theory alone
  cannot supply. The Pochhammer product, whose roots 1..D have unbounded height, is the case it excludes.
- *Scope.* The theorem covers depth k ≤ log₂ X + Δ for quadratic 𝒢. It says nothing about k = (log X)^A, where
  the prefix bound in step 3 drops below 1.
- Under GRH, effective Chebotarev on prefixes of depth m ≈ log log X lets Δ grow like c·log log log X (sketch,
  UNVERIFIED constants).
- The theorem is for compositions, not general SLPs. A general SLP has no outer prefix of bounded degree.

### 4.5 Random coefficients. PROVED, trivial

If c_1..c_k are independent and uniform mod p, then E_c Z_p(F) = 1 exactly, since E_c η(u − c) = 0 at each
layer. So Pr_c[Z_p ≥ μp] ≤ 1/(μp). With integer c_i uniform in [0, H), the same holds up to the factor
(1 + p/H)^k. So generic large coefficients are excluded at every single prime, but the interesting case is small
fixed constants.

### 4.6 What remains: Conjecture QC (codex's restricted target, endorsed)

For fixed finite C ⊂ Z and A, B > 0, uniformly over words of length ≤ (log X)^A,
#{p ∈ (X, 2X] : δ_p(F) ≥ (log X)^{−B}} = o(π(2X) − π(X)).

What is proved:
- OUT: the outer constant must be −s²;
- PFX: fixed prefixes with vanishing FPP fail;
- MIN: excess depth must → ∞;
- DEF: near-minimal depth needs almost-total squareness per prime.

The open class is **words of growing excess depth all of whose prefixes are (asymptotically) intersective**,
i.e. have roots at almost every prime. Words whose prefixes carry an integer backbone (a chain 0 ← y_1 ← y_2 ← …
of integers with y_j² + c = y_{j−1}; all such y_j lie in [−2, 2] for C = {−2..2}) are intersective and are what
the E2 maximisers look like.

---

## 5. Where the limit bites (summary)

| Question | Per-prime status | Uniform status | Barrier proved here |
|---|---|---|---|
| O1, χ of order with large prime factor ℓ | **DL-hard**: τ-correlation ⇒ DL mod ℓ (O1-DL) | O1 ⇒ factoring + DL at ≥ μ/2 of primes (O1-N) | conditional on prime-field DL; Weil unconditional only for w < √p |
| O1, χ of smooth order | easy (PH-SLP) | N-independent only when p − 1 is smooth (PH with fixed m) | = order annihilation (OA) |
| O1, monomial Φ | — | — | BGK power saving for ord χ ≥ p^ε (codex) |
| dense-torus, any group | **trivial** (x^{(p−1)/2} − 1) | needs uniformity | none possible from Galois alone (GAL-0); permanent-hard in general (BB) |
| compositions, outer −c non-square | — | fails at ≈ ½ of primes | OUT |
| compositions, fixed prefix tower | — | roots at o(1) of primes | PFX |
| compositions, degree ≤ 2^Δ X | — | fails at a κ-fraction of primes | **MIN** |
| compositions, polylog depth, intersective prefixes | — | open | Conjecture QC; E2/E3 support it |

---

## 6. Experiments

*(to be filled after the runs)*

---

## 7. Candidate cards

*(to be filled)*

---

## 8. Reproduction

```
.venv/bin/python -m pytest -q tests/test_prepread.py
.venv/bin/python -m src.prepread.run_experiments e1   # ~1 h, ≤ 0.4 GB
.venv/bin/python -m src.prepread.run_experiments e2
.venv/bin/python -m src.prepread.run_experiments e3
.venv/bin/python -m src.prepread.run_experiments e4   # needs e1 output
.venv/bin/python -m src.prepread.run_experiments e5
```

## 9. Review record

*(to be filled)*
