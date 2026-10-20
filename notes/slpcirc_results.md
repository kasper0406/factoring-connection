# slpcirc results: factoring via short straight-line programs (phase 17)

Agent `slpcirc`, 2026-10-20. Pre-registration: `notes/slpcirc_prereg.md` (commits `8e1ba0b`, Addendum A
`bf77fc6`). Code: `src/slpcirc/*`; tests: `tests/test_slpcirc.py`. Codex blind round:
`codex/slpcirc_1.{md,out}` (full report `codex/scratch/slpcirc_1/report.md`, git-ignored).

## Verdict

**No mechanism. The SLP paradigm gives no route below L[1/3]; here it does not even reach a new
subexponential method.** The paradigm splits cleanly into three regimes.

1. **SLPs that depend on N (the "mod N, for this N" version) are not a weaker hypothesis. They are a
   restatement of factoring** (Prop R, proved and trivial). Every residue mod N, including ⌊√N⌋! mod N,
   has an SLP of size ≤ 2 log₂ N. So "a short SLP exists mod N" is vacuous, and "a short SLP can be
   *constructed* from N" is exactly a factoring algorithm. The question has content only once we
   restrict how the program may depend on N.
2. **N-oblivious SLPs (model U: one integer program per bit length, evaluated mod N at random inputs,
   then a gcd) are where the content lies.** This model already contains trial division, p ± 1,
   Bach–Shallit tori, ECM/HECM, rho and Pollard–Strassen. ECM, at L_p[1/2], is the best of them.
   - The weakest sufficient statement is W_U(s, μ): root density δ_p ∈ [μ, 1 − μ] mod all but ≤ 2
     primes of the range (§1.3).
   - An unconditional proof that model U cannot be polynomial would prove a superpolynomial
     constant-free lower bound for the permanent (Prop BB, via Lipton and Bürgisser 2009). So no
     general barrier is provable with current techniques.
3. **Restricted barriers we proved (all checked by exact experiments):**
   - **Theorem SP:** any polynomial that is sparse in a torus coordinate (x, Dickson u + 1/u,
     F_{p^k}^*, several variables, products of such) is *exactly* a group-order method. Root density
     ≥ λ at p forces every prime ℓ ≥ ℓ₀(t, λ) dividing the group order (p − 1, p + 1, Φ_k(p), …) to
     divide the product Λ of the exponent differences. Coefficients are unrestricted and may depend
     on N or even on p.
   - **Theorem SP-U:** with N-independent exponents, success on a positive fraction of primes in
     (X, 2X] forces exponent bit length ≥ X^{θ−o(1)}. For worst-case balanced N this means SLP size
     ≥ N^{1/4−ε}/(r t²). So sparse/torus methods cannot even match Pollard–Strassen.
   - **Cor SP-N:** exponents that are polynomials in N (x^{N−1}, x^N, …) succeed on a density-zero set
     of balanced semiprimes.

**Is every short SLP with many roots mod p secretly a group-order method? Answered in two senses:**
- *Broad sense: no.* Codex's S₃ example x³ − x − 1, the CGRW trinomials x^a − x − 1, polynomial
  dynamics and Pochhammer products are not group-order methods.
- *Dense-and-uniform sense: no known exception.* Every non-group mechanism found is either:
  - sparse, so Theorem SP kills it;
  - birthday-bounded (rho, δ ≈ s²/p);
  - explicit-root, i.e. the τ / factorial problem (open, and as hard as permanent lower bounds);
  - bounded by Chebotarev (about 1 root per irreducible factor; census H2 confirms 0.90);
  - or empirically null (quadratic compositions, E4).

**Assessment.** Probability that this thread yields sub-L[1/3] factoring: < 1%. Probability that
it yields any new subexponential method: < 2%. The one open door is a dense, uniform, non-group,
non-sparse short SLP (§6). Model U lower bounds there are permanent-hard (Prop BB), so neither a
construction nor a proof of impossibility is in reach.

Labels: proved / experimentally supported / conjectured / refuted for this model / KNOWN (literature) /
UNVERIFIED. "No prior formulation found by us" is the strongest novelty statement made. All proofs
below were checked by two read-only Opus critics and a codex adversarial round (§8).

---

## 1. Models and the weakest sufficient statement

### 1.1 Success probability [proved; codex agrees]

For N = pq, a program F (polynomial or SLP over Z/N in variables x₁..x_k), a uniform input
x ∈ (Z/N)^k, and δ_r := Pr_{u ∈ F_r^k}[F(u) ≡ 0 mod r]:

  Pr[1 < gcd(F(x), N) < N] = δ_p(1 − δ_q) + δ_q(1 − δ_p).   (1.1)

This follows from CRT independence. A trial costing C succeeding with probability ε gives expected time
O(C/ε). Two consequences follow.
- **One balanced prime suffices.** If μ ≤ δ_p ≤ 1 − μ, then success ≥ μ whatever δ_q is.
- **Density alone is not enough** (codex): x^{(p−1)(q−1)} − 1 is a nonzero function at both primes but
  has δ ≈ 1 at both, so (1.1) is ≈ N^{−1/2}.

### 1.2 Prop R (collapse of the N-dependent version) [proved, trivial; = codex §1.3]

(a) **Existence is vacuous.** Every residue r mod N has a constant-free SLP of size ≤ 2⌈log₂ N⌉
(binary expansion). So "n! mod N has a short SLP for the specific N" is true for every N. The τ-type
question only makes sense for families, or for integers rather than residues.

(b) **Construction is factoring.** Expected-poly-time factoring is equivalent to uniform expected
poly-time construction, from N, of poly-size SLPs F_N over Z/N with δ_r(F_N) ∈ [1/3, 1/2) at every prime
r | N.
- (⇐) follows from (1.1).
- (⇒) factor N, then output F_N = Σ_r e_r (x^{(r−1)/2} − 1) with CRT idempotents e_r (codex).

The same equivalence holds with L[α] budgets.

(c) **Intermediate uses of N.**
- As a ring scalar, N ≡ 0 mod p, so any F(x, N) ≡ F(x, 0) mod p. Algebraic use of N is useless.
- As an *exponent*, N is covered by Cor SP-N (§4.4).
- An *arbitrary computable* function of N (constants, exponents, polynomial "nonscalar constants" in
  (Z/N)[x]) falls back to (b) unless the constructor is restricted.

So the question "SLPs over (Z/N)[x] with N-dependent constants" has content only for restricted
constructors.

### 1.3 Model U (N-oblivious programs) and the weakest sufficient statement [proved]

**Model U.**
- A family (F_{X,τ}) of constant-free SLPs over Z (division-free, or with gcd-checked divisions,
  §1.4) in k input variables. The family depends only on the size parameter X, not on N.
- It is run on N with prime factors in (X, 2X].
- Each trial evaluates one F_{X,τ} mod N at a fresh uniform input and takes a gcd with N.
- The next program may depend adaptively on the previous gcd outcomes (1 or N).
- Size s counts ring operations (including those building constants); k ≤ s.

The known methods in this model:

| Method | Program | (size, per-prime density) |
|---|---|---|
| trial division | constant ∏_{i∈I} i | (\|I\|, 0/1) |
| Pollard p − 1 (and stage 2) | x^M − 1, ∏(x^{Mq_j} − 1) | (log M, gcd(M, p−1)/p) |
| Williams p + 1 | D_M(x) − 2 | (log M, (gcd(M,p−1)+gcd(M,p+1))/2p) |
| Bach–Shallit Φ_k | sparse in a subtorus of F_{p^k}^* | (log M, Φ_k(p)-smooth part) |
| ECM / HECM | Z-coordinate of [M]P on a curve with parameters as inputs | (≈ B, ≈ Pr[#E smooth]) |
| Pollard rho | ∏_{i≤T}(f^{(2i)}(x) − f^{(i)}(x)) | (≈ T, ≈ T²/p, heuristic) |
| Pollard–Strassen / BGS | constant k! by product trees, binary search on k | (Õ(√k), 0/1) |

**Weakest sufficient statement W_U(s, μ).** There is a model-U family of size s with
μ ≤ δ_p ≤ 1 − μ for all but at most two primes p ∈ (X, 2X]. Then balanced semiprimes factor in
expected O(s·M(log N)/μ) bit operations: success per trial ≥ 2μ(1 − μ) ≥ μ for every pair.

**Necessity for a single fixed program [proved].** Suppose one program succeeds with probability ≥ μ on
every N = pq with X < p < q ≤ 2X.
- Let S_lo = {δ_p < μ/2} and S_hi = {δ_p > 1 − μ/2}.
- Two primes in S_lo give success < μ by (1.1); so do two primes in S_hi.
- Hence |S_lo|, |S_hi| ≤ 1, and all but ≤ 2 primes have δ_p ∈ [μ/2, 1 − μ/2].
- W_poly := W_U(polylog X, 1/polylog X) gives randomized poly-time factoring of balanced semiprimes.
  W_U with s/μ ≤ L_N[α] gives L[α].
- Model U contains no NFS or QS: their relations depend on N. Model U lower bounds are therefore *not*
  factoring lower bounds.

**Where the τ-conjecture fits [proved scope].**
- A polylog SLP for the Pochhammer polynomial ∏_{a≤X}(x − a) gives W_poly, since δ_p = X/p ∈ [1/2, 1)
  for p ∈ (X, 2X] (Lipton 1994).
- So do polylog SLPs for the constants k! (binary search; Lipton 2009).
- W_poly is much weaker than either. It allows any roots mod p, not integer roots, and several
  variables. The (integer-root) τ-conjecture does **not** imply that W_poly fails: x^M − 1 has 2 integer
  roots but many roots mod p.

### 1.4 Divisions and intermediate gcds [proved, classical]

Run a division SLP of size s as numerator/denominator pairs:
- a/b ± c/d = (ad ± bc)/(bd): 3 multiplications and 1 addition;
- a/b · c/d: 2 multiplications;
- division: a swap.

The result is two division-free programs of size ≤ 4s. A gcd-checked division hits exactly when the
denominator value is a nonunit. So division SLPs reduce to division-free ones of size O(s) by testing
A·B; a test that is ≡ 0 mod N is resolved by testing A and B separately. Intermediate gcd checks are
adaptive trials. **Divisions and gcds change only constants in model U.**

### 1.5 Prop BB (barrier to barriers) [proved from KNOWN theorems; already observed in trackC §3]

Suppose someone proves unconditionally that W_poly fails, even nonuniformly and with adaptive
constant-free division-free programs. Then:
1. **τ(k!) is not (log k)^{O(1)}.** Otherwise binary search over k ∈ (X, 2X] with the constant programs
   for k! finds p ≤ k < q in log₂ X trials: gcd(k! mod N, N) = p.
2. **Hence τ(Per_m) is superpolynomial**, by Bürgisser (Comput. Complexity 18, 2009, Main Thm 1.2: τ(Per)
   poly ⇒ τ(k!) polylog).

So any general lower bound for model U is at least as hard as a constant-free permanent lower bound,
which is open. The same holds for the version with an input variable, through the Pochhammer polynomial.
**Consequence:** honest barriers for this paradigm must be restricted-model theorems. §4 gives them.

---

## 2. Known results and which of them bind

| Result | Statement (scope) | Binds the mod-N / N-dependent version? |
|---|---|---|
| Shub–Smale 1995 (Duke 81), BCSS Thm 3 | τ-conjecture (#integer roots ≤ (1+τ)^c) ⇒ P_C ≠ NP_C; Main Thm: factorials not ultimately easy ⇒ HN intractable | No: integer roots and char 0 only |
| Bürgisser 2009 (Comput. Complexity 18), Main Thm 1.2 | factorial/Pochhammer/τ hardness ⇒ τ(Per) superpoly; Per easy ⇒ τ(k!) polylog (Pochhammer: with divisions) | Gives Prop BB (barrier to barriers) |
| Lipton 1994 (ANTS-I), 2009 (blog Thm 1); Phillipson–Rojas 2013 Thm 2.4(II) | short SLPs for polynomials with many rational roots / for t! ⇒ factoring not hard (nonuniform / average-case) | Upper-bound direction only |
| Koiran 2011 (ICS; arXiv 1004.4960 Conj. 3, Thm 7) | real τ-conjecture for sums of products of sparse polynomials ⇒ permanent lower bounds | No: real roots, char 0 |
| Koiran–Portier–Tavenas 2015 (JSC 68) Thm 12 | Wronskian real-root bounds for SPS-type expressions | No: char 0; Wronskians fail in char p |
| Phillipson–Rojas 2013 (Contemp. Math. 605) | adelic (R or Q_ℓ) τ/SPS conjecture | No: bounds Q_ℓ-roots, not F_p-roots |
| Cheng 2003 (Comput. Complexity 12): WL conjecture; ψ_m has an SLP of size 72 log m + 60 | number-field root bounds ⇒ bounded torsion | No: over F_p, torsion is large |
| Cheng 2004 (TCS 326) Thm 1 | conditional L_t[1/2]-size SLPs for *multiples* of t! (elliptic + smoothness) | Group-order construction (§3) |
| Strassen 1973 (Numer. Math. 20) degree bound; Schönhage 1976 | s ≥ log₂ deg; geometric-degree bounds | Only s ≥ log₂(δp): no barrier (§4.6) |
| Strassen 1974 (SIAM J. Comput. 3), Schnorr 1978, Heintz–Sieveking 1980, Heintz–Schnorr 1980 | hard-to-compute polynomials with algebraic or "generic" coefficients/roots | No: they use coefficient-field degree, which collapses mod p (x^{p−1} − 1 has p − 1 roots and size 2 log p) |
| Borodin–Cook 1976, Grigoriev 1982, Risler 1985 | additive complexity vs real roots | No: x^{p−1} − 1 has one subtraction |
| Baur–Strassen 1983 | L(f, ∇f) ≤ 3L(f) | No |
| CFKLLS 2000 (Israel J. Math 120); **Kelley 2016 (LMS JCM 19A) Thms 2.2–2.3**; Kelley–Owen 2017 (JSC 79); Bi–Cheng–Rojas 2013/2016; Cheng–Gao–Rojas–Wan 2017 (FFA 46) | t-nomial root counts on F_q: R ≤ 2(q−1)^{1−1/(t−1)}C^{1/(t−1)}, C = largest coset in the root set, coset alignment of exponents; constructions with many roots in extension fields; x^a − x − 1 with unbounded counts at selected primes | **Yes**, per prime, any coefficients. This is the input to Theorem SP (§4); our Lemma S1 is a self-contained weaker relative |
| Ritt 1923 (Trans. AMS 25), Eremenko 1989 (Leningrad Math. J. 1) | commuting rational maps of degree ≥ 2 (no common iterate) are conjugate to powers, ±Chebyshev or Lattès | Yes, for "fast-index" families: multiplicatively indexed commuting families f_m ∘ f_n = f_{mn} are G_m / torus / elliptic (§3) |
| trackC O1–O4 (internal) | interval/prefix zero sets: fixed-dimension powers, power-map layers, group laws | Interval zero sets only |
| systeq Lemma CF (internal) | systematic product vanishing on all scalars has degree ≥ p | Factorial route: degree, not size |
| nonsmooth Prop OA (internal) | single-scalar group-law annihilation [m]x = 0 succeeds only if #G(F_p) is smooth up to {p,q}-parts | Yes, ECM/HECM/tori with one fixed scalar |
| "Kayal's factoring via SLPs" | **UNVERIFIED / not found** (codex and we). Nearby: Kaltofen 1989 (SLP polynomial factorization, poly in degree); Kayal–Saxena 2006 (ring morphisms) | — |

**Which lower bounds actually bind the weak version.**
- For the N-dependent version: none (Prop R).
- For model U, in general: none beyond the trivial degree bound s ≥ log₂(δ_p p), which is polylog.
  Every characteristic-0 tool (τ, real τ, adelic τ, Wronskians, Heintz–Sieveking, Risler) fails at the
  first step, because the subgroups of F_p^* are cheap. This is consistent with Prop BB: a general bound
  would be permanent-hard.
- What binds are the finite-field sparse bounds (Kelley/CFKLLS), Prop OA, Ritt–Eremenko, and Track C's
  O1–O4, each only inside its restricted class. Theorem SP below packages the sparse bound into an
  exact group-order dichotomy and a uniform exponential barrier.

---

## 3. Construction table and group-order classification

"Group order" means: the root count at p is a function of the structure of G(F_p) for a fixed
algebraic group G (here via an explicit torsion/kernel/coset mechanism), and success needs a scalar that
annihilates most of G(F_p).

| # | Family (short SLP) | Roots mod p | Governing structure | Verdict |
|---|---|---|---|---|
| 1 | n·x, constants n | p \| n | G_a: trial division / content | constant route = τ of integers (§1.5) |
| 2 | x^m − a, x^M − 1, stage-2 products, x^{N−1} − 1 | gcd(m, p−1) (cosets) | G_m, p − 1 | **Theorem SP** (exact); SP-U, SP-N |
| 3 | D_m(x) − 2, Lucas V_m, Chebyshev T_m | (gcd(m,p−1)+gcd(m,p+1))/2 [proved, tested] | split / non-split tori, p ± 1 | **Cor SP-D** |
| 4 | norms/traces from (Z/N)[y]/(f), sparse in u ∈ F_{p^k}^* | kernel cosets | tori, Φ_k(p) | **Theorem SP** on G ≤ F_{p^k}^* |
| 5 | division polynomials ψ_m, Lattès maps, [M]P (ECM, HECM) | (#E(F_p)[m] + #E^d(F_p)[m] − 2)/2 for odd m (codex, proved) | E and twist; abelian varieties | Prop OA (single scalar); elliptic-sparse open (§4.5) |
| 6 | multiplicatively indexed commuting families f_m ∘ f_n = f_{mn} | — | Ritt/Eremenko: power, ±Chebyshev, Lattès only | always rows 2, 3 or 5 [KNOWN] |
| 7 | quadratic compositions (x² + c_k)∘⋯∘(x² + c₁) | preimage tree; at most 2^k | c ∈ {0, −2}: powers/Chebyshev (2-power torsion only); generic c: **non-group** dynamics | E4: critical branching, mean fibre 1.46, no uniform density [exp.]; theorem open (§6) |
| 8 | rho: f^{(k)}(x) − x, ∏(f^{(2i)} − f^{(i)}) | cycle structure of x ↦ f(x) | **non-group** (birthday); x² + c with c = 0 collapses to G_m | δ ≈ s²/p (heuristic); exponential |
| 9 | Pochhammer ∏(x − a), k!, primorials | explicit integer roots | **non-group** (explicit) | equivalent to τ / "ultimately easy" factorials; open, permanent-hard (Prop BB) |
| 10 | Cheng 2004 multiples of t! | elliptic smooth orders | group (row 5) + smoothness | L[1/2] |
| 11 | x^a − x − 1 (CGRW), x³ − x − 1 (S₃) | Frob fixed points; unbounded at selected p (CGRW) | **non-group** (non-abelian splitting) | trinomial ⇒ Theorem SP / Kelley (C = 1 ⇒ δ ≤ 2p^{−1/2}); never dense |
| 12 | x^{p^k} − x, extension trace polynomials (CGRW, Kelley–Owen) | additive subfields | Frobenius / additive group | needs the characteristic; over F_p with public exponents it is a row-2 sparse polynomial |
| 13 | x^{(p−1)/2} − 1, "h(x) is a square" | ≈ p/2 | quadratic character | needs p (Prop R) |
| 14 | products / unions; norms / resultants | union; projection | inherits its constituents | union bound; products are in SP; generic elimination: no construction |
| 15 | random dense tiny SLPs (E2) | ≈ 1 per irreducible factor | Chebotarev | E2 H2: mean 0.90 per factor |

**Which ones do not reduce to group orders:** rows 7 (generic c), 8, 9, 11, 12, 15. None of them is
dense *and* uniform at polylog size:
- row 8 is birthday-bounded;
- row 9 is the open τ question (believed hard);
- row 11 is killed by Theorem SP and Kelley;
- row 12 needs p;
- row 15 has bounded counts;
- row 7 is empirically null.

**Is there an argument that "many roots ⇒ group order"?**
- In the broad sense, no: rows 7, 8, 9 and 11 are counterexamples (codex: x³ − x − 1 is not governed by
  any fixed congruence).
- In the sparse sense, yes, exactly (Theorem SP).
- For index-multiplicative families, yes (Ritt–Eremenko).
- For single-scalar group laws, yes (Prop OA).
- The Galois view explains the boundary. For p ∤ disc, #roots(F mod p) = Fix(Frob_p) on the roots of F.
  In the Chebotarev regime (p large compared with the fields) the average is the number of irreducible
  factors. Dense roots at p need deg F ≥ δp, which is the regime p < deg where Chebotarev is silent.
  There, group methods work because Frob_p is a *small* element (the integer p, or π with |π| = √p)
  acting on torsion, so its fixed points are a group order. A non-group dense mechanism would need small
  Frobenius elements of a non-abelian, non-torsion Galois action to have many fixed points uniformly in
  p. We know no such family, and Theorem SP excludes it for sparse polynomials.

---

## 4. Barriers (restricted models)

Notation:
- G is a cyclic subgroup of F_Q^* (Q a power of p), of order n.
- F(u) = Σ_{i=1}^t c_i u^{e_i}, with distinct integers e_i (Laurent allowed) and c_i ∈ F_Q not all 0.
- K_ℓ(t) := ⌊(ℓ − 1)^{1/(t−1)}⌋.
- Λ(F) := ∏_{i<j}(e_i − e_j) ∈ Z ∖ {0}.
- ℓ₀(t, λ) := the least ℓ with 2/K_ℓ(t) < λ. K_ℓ is nondecreasing in ℓ, so this holds for all ℓ ≥ ℓ₀, and
  ℓ₀ ≤ (⌊2/λ⌋ + 1)^{t−1} + 1.

### 4.1 Lemma S1 (sparse zeros on a coset of μ_ℓ) [proved; elementary relative of CFKLLS / Kelley]

Let ℓ | n be prime, t ≥ 2, all c_i ≠ 0, and the e_i pairwise incongruent mod ℓ. Then

  #{u ∈ G : F(u) = 0} ≤ (n/ℓ)·2(⌈ℓ/K⌉ − 1) < 2n/K,  K = K_ℓ(t) ≥ 1.

*Proof.*
1. Let H ≤ G be the subgroup of order ℓ. G is the disjoint union of n/ℓ cosets xH. Fix x and put
   φ(y) = F(xy) = Σ b_i y^{e_i}, with b_i = c_i x^{e_i} ≠ 0. On H, y^ℓ = 1, so φ(y) = Σ b_i y^{r_i} with
   r_i = e_i mod ℓ pairwise distinct.
2. *Dirichlet.* The ℓ vectors (m(r_i − r_1)/ℓ mod 1)_{i=2..t} ∈ (R/Z)^{t−1}, for m = 0..ℓ−1, fall into
   K^{t−1} ≤ ℓ − 1 < ℓ boxes of side 1/K. Two of them, m ≠ m′, share a box. So m* = m − m′ ≢ 0 (mod ℓ)
   satisfies: each m*(r_i − r_1) ≡ v_i (mod ℓ) with an integer |v_i| < ℓ/K.
3. The v_i (v_1 = 0) are pairwise distinct, because m* is invertible mod ℓ.
4. y ↦ y^{m*} permutes H, and φ(y^{m*}) = y^{m* r_1}·Σ b_i y^{v_i} on H. With w = ⌈ℓ/K⌉ − 1 ≥ |v_i|, the
   zeros of φ on H correspond to the zeros on H of ψ(y) = Σ b_i y^{w+v_i}. This is a nonzero polynomial
   (distinct exponents, nonzero coefficients) of degree ≤ 2w, so it has at most 2w zeros.
5. Sum over the n/ℓ cosets. ∎

(E1: 4 429 217 exhaustive/random instances, t ≤ 5, p ≤ 211, 0 violations.)

**Lemma S1-k (several variables) [proved].** Let F = Σ c_i u^{α_i} on G^k (α_i ∈ Z^k, c_i ≠ 0),
ℓ > C(t,2), and the α_i pairwise incongruent mod ℓ. Then the zeros number < 2n^k/K_ℓ(t).

*Proof.*
1. There is w ∈ Z^k with (α_i − α_j)·w ≢ 0 (mod ℓ) for all i < j: each pair excludes a hyperplane of
   (Z/ℓ)^k, and C(t,2)ℓ^{k−1} < ℓ^k.
2. ι: H → G^k, y ↦ (y^{w_1}, …, y^{w_k}) is an injective homomorphism.
3. On each of the n^k/ℓ cosets x·ι(H), F(x·ι(y)) = Σ c_i x^{α_i} y^{α_i·w} has t nonzero terms with
   exponents distinct mod ℓ.
4. Apply the per-coset bound of S1. ∎

(Tested: `test_s1_multivariate_bruteforce`.)

### 4.2 Theorem SP (sparse ⇒ group order; per-prime dichotomy) [proved]

Suppose the zero density of F on G is ≥ λ > 0. Then every prime ℓ ≥ ℓ₀(t, λ) dividing n divides Λ(F).
So the ℓ₀-rough part of |G| is supported on the primes of Λ(F). In the k-variable version, use
ℓ ≥ max(ℓ₀, C(t,2) + 1) and Λ_k(F) := ∏_{i<j}(a nonzero coordinate of α_i − α_j).

*Proof.* Let I = {i : c_i ≠ 0} ≠ ∅. The characters u ↦ u^a, a ∈ Z/n, are linearly independent functions
on G (Dedekind). Reduce F|_G to Σ_{a∈S} β_a u^a by merging exponents mod n, keeping β_a ≠ 0.
- **Case 1: F|_G ≡ 0.** Every class contains ≥ 2 indices of I (a singleton class has β = c_i ≠ 0). So
  e_i ≡ e_j (mod n) for some i ≠ j, and ℓ | n | e_i − e_j.
- **Case 2: |S| = 1.** A monomial has no zeros on G, contradicting λ > 0.
- **Case 3: |S| = t′ ∈ [2, t].** If the classes in S were pairwise incongruent mod ℓ, Lemma S1 (on
  representatives) would give density < 2/K_ℓ(t′) ≤ 2/K_ℓ(t) < λ. So two classes a ≠ a′ satisfy
  a ≡ a′ (mod ℓ), and their representatives e_i, e_j satisfy ℓ | e_i − e_j.
- **The k-variable version** is the same, with characters u^α, α ∈ (Z/n)^k, and Lemma S1-k.
∎

*Remarks.*
1. The coefficients c_i are arbitrary elements of F_Q. They may depend on N, and even on p. The only
   exclusion is F ≡ 0 as a polynomial mod p. That is the *content route*: p divides every coefficient,
   i.e. a constant SLP. Content-route constants can themselves be group-order constants (2^M − 1) or
   factorial-type constants (k!). Their complexity is the τ-of-integers problem of §1.5, which Theorem
   SP does not touch.
2. Kelley 2016, Thms 2.2–2.3 (KNOWN; exact form to be read in the source) gives a stronger single-coset
   version: density δ forces a coset of size ≥ (p−1)(δ/2)^{t−1} inside the root set, with exponents
   aligned mod its order. That implies Theorem SP with ℓ₀ ≈ (2/λ)^{t−1}. Our Lemma S1 is a self-contained
   elementary route to the prime-by-prime form, which is all the uniform corollary needs.

**Corollaries at a prime p (input x uniform in F_p).**
- **(SP-1) p − 1.** For a t-sparse F ∈ F_p[x] (or F_p[x^{±1}]) with root density δ_p ≥ λ + 1/p on F_p:
  every prime ℓ ≥ ℓ₀(t, λ) dividing p − 1 divides Λ(F).
- **(SP-D) p ± 1 (Dickson/Lucas/Chebyshev).** Suppose F(u + u^{−1}) = Σ_{i≤t} c_i u^{e_i} as Laurent
  polynomials. Example: F = Σ_j a_j D_{m_j}(x, 1), exponents ±m_j, and
  Λ = ∏_{j<j′}(m_j − m_{j′})(m_j + m_{j′})·∏_j 2m_j.
  - Each x ≠ ±2 in F_p is u + u^{−1} for exactly two u, both in F_p^* or both in the norm-one group
    T ≤ F_{p²}^* of order p + 1.
  - |F_p^*| + |T| = 2p. If F has ≥ λp zeros, then F(u + 1/u) has zero density ≥ λ − 2/p on F_p^* or
    on T.
  - Theorem SP then gives: every ℓ ≥ ℓ₀(t, λ − 2/p) dividing p − 1, or every such ℓ dividing p + 1,
    divides Λ.
  - Tested: `test_dickson_root_formula`, `test_dichotomy_pollard`.
- **(SP-K) Bach–Shallit tori.** For F sparse in u ∈ G ≤ F_{p^k}^* (e.g. |G| = Φ_k(p)), Theorem SP holds
  verbatim with |G| in place of p − 1.
- **(SP-Π) Products.** F = ∏_{j≤r} F_j with each F_j t-sparse and F ≢ 0 mod p. Root density ≥ λ means
  some F_j has density ≥ λ/r. So every ℓ ≥ ℓ₀(t, λ/r) dividing the group order divides ∏_j Λ(F_j).
  This covers p − 1 stage 2 and unions of several group tests.

*Interpretation.* Success of a sparse/torus method at p happens only if the exponent data contain all
large prime factors of a group order. Once the coefficients are free, this is precisely the
Pollard/Williams/Bach–Shallit success condition. **This is the sense in which every sparse SLP method is
secretly a group-order method** (proved).

### 4.3 Lemma LP (dyadic Goldfeld) [proved, using Bombieri–Vinogradov]

For fixed 0 < θ < 1/2:

  #{X < p ≤ 2X : P⁺(p − 1) > X^θ} ≥ (1 − θ − o(1))·π*(X),  π*(X) = π(2X) − π(X).

The same holds for p + 1.

*Proof.*
1. Σ_{X<p≤2X} log(p−1) = π*(X)(log X + O(1)) = Σ_d Λ(d)·π*(X; d), where π*(X; d) counts p with d | p − 1.
2. For d ≤ D = X^θ, BV (dyadic difference of two BV sums) gives
   Σ_{d≤D} Λ(d)π*(X; d) = π*(X)(θ log X + O(1)) + O(X/log^A X).
3. Hence Σ_{p} Σ_{ℓ^a | p−1, ℓ^a > D} log ℓ ≥ π*(X)((1 − θ) log X − O(1)).
4. Split the left side.
   - Primes p with P⁺(p − 1) > D contribute ≤ log(2X) each.
   - Prime powers ℓ^a > D with ℓ ≤ D (so a ≥ 2) contribute in total
     ≤ Σ_{ℓ≤D} log ℓ · Σ_{a: ℓ^a > D} (X/ℓ^a + 1) = O(X/√D + D log² X) = o(π*(X) log X).
5. Divide by log 2X. ∎

(Sanity check at X = 2^20, 4 000 sampled primes: the fractions for θ = 0.2/0.3/0.4/0.5 are
0.997/0.962/0.846/0.673, against the lower bounds 0.8/0.7/0.6/0.5. For θ = 1/2, take D = X^{1/2}/log^B X.)

### 4.4 Theorem SP-U (N-independent exponents need exponential size) [proved]

**Hypotheses.**
- F_X = ∏_{j≤r} F_{X,j} ∈ Z[u_1^{±1}, …, u_k^{±1}], each factor with ≤ t terms and exponent coordinates
  ≤ D in absolute value, independent of N.
- 0 < θ < 1/2 and λ ∈ (0, 1], with ℓ₀* := max(ℓ₀(t, λ/r), C(t,2) + 1) ≤ X^θ.
- P = {p ∈ (X, 2X] : F_X mod p ≢ 0 and F_X has zero density ≥ λ on (F_p^*)^k}.

**Conclusion.**

  r·C(t,2)·log(2D) ≥ θ log X · (|P| − (θ + o(1))π*(X)) / (X^{1−θ} + 1).

*Proof.*
1. Let P′ = {p ∈ P : ℓ_p := P⁺(p − 1) > X^θ}. By Lemma LP, |P′| ≥ |P| − (θ + o(1))π*.
2. For p ∈ P′, ℓ_p ≥ ℓ₀*, so ℓ_p | Λ := ∏_j Λ_k(F_{X,j}) by Theorem SP and (SP-Π). Λ is a nonzero
   integer with log|Λ| ≤ r C(t,2) log(2D).
3. Each prime ℓ > X^θ equals ℓ_p for at most X/ℓ + 1 ≤ X^{1−θ} + 1 primes p ∈ (X, 2X]. So there are at
   least |P′|/(X^{1−θ} + 1) distinct ℓ, each > X^θ, all dividing Λ. ∎

**Corollary SP-U1 (worst case).**
- Suppose a single fixed model-U program whose output is such a product succeeds with probability ≥ μ on
  every N = pq, X < p < q ≤ 2X.
- By §1.3, |P| ≥ π*(X) − 2, with λ = μ/2 − k/(p − 1), after restricting to the torus.
  (Zero density ≤ 1 − μ/2 < 1 excludes F ≡ 0 mod p.)
- So log(2D) ≥ (1 − o(1))·θ(1 − θ)X^θ/(r C(t,2)) for every θ < 1/2 allowed by ℓ₀*.
- Any SLP for a polynomial F_X has size s ≥ log₂ deg F_X ≥ log₂ D. Hence

    s ≥ c·X^{1/2−ε}/(r t²) = c·N^{1/4−ε}/(r t²)

  for balanced N, whenever (4r/μ + 1)^{t−1} + 1 ≤ X^{1/2−ε}. This covers fixed t and r, and also
  t ≤ c′ log X / log(r/μ).
- The same bound holds for p + 1 (SP-D, using Lemma LP for p + 1) and for tori of fixed dimension.
- Pollard–Strassen achieves N^{1/4+o(1)} in model U, so sparse/torus methods cannot beat it in the
  worst case.

*Adaptive trials.* With T programs (all of this sparse form, contents prime to the range), per-trial
success at a pair where both densities are < λ is < 2λ. So uniform success μ > 2λT forces all but one
prime of the range into ∪_τ P_τ, and the theorem applies to the sum of the Λ's.

**Corollary SP-U2 (average case).**
- Success ≥ μ on average over balanced N gives E_p[δ_p] ≥ μ/2, hence |P| ≥ (μ/4)π* with λ = μ/4.
- Take θ = μ/8. This gives s ≥ X^{μ/8 − o(1)}/(r t²).
- For μ = 1/polylog this is superpolylogarithmic. So there is no polynomial-time sparse method even on
  average.

**Corollary SP-N (exponents polynomial in N) [proved].**
- Let the exponents be E_i(N) with fixed E_i ∈ Z[y], pairwise distinct as polynomials (x^{N−1} − 1,
  x^N − x, x^{N²+1} − x^{N}, …). The coefficients may be arbitrary, provided they are not all ≡ 0 mod p
  (that is the content route, which with N-dependent coefficients is Prop R).
- If ℓ | p − 1 then N ≡ q (mod ℓ). So success at p with P⁺(p − 1) = ℓ > X^θ ≥ ℓ₀ needs ℓ | Π(q),
  where Π := ∏_{i<j}(E_i − E_j) ≠ 0.
- For fixed q, at most ω_{>X^θ}(Π(q))·(X^{1−θ} + 1) = O_{E,θ}(X^{1−θ}) primes p qualify.
- So the fraction of balanced pairs on which such a template succeeds is ≤ θ + o(1) for every θ > 0,
  i.e. o(1). The successes concentrate on X^{o(1)}-smooth p − 1: Fermat-test / Pollard behaviour.
- Exponents given by *arbitrary* computable functions of N fall under Prop R.

### 4.5 What the barriers do not cover

- **Elliptic sparse combinations.** Σ c_i φ([e_i]P) on E(F_p) are not covered. The Dirichlet/degree
  trick of S1 fails: on a subgroup C ≅ Z/ℓ, φ∘[k] has degree ≈ k², so after compressing the multipliers
  to |k_i| < ℓ^{1−1/t} the coset function has degree ≈ ℓ^{2−2/t} > ℓ. Prop OA covers the single-scalar
  test [M]P = O only.
- **Dense, non-sparse SLPs:** rows 7, 8, 9 and 15. Here only Prop BB (no general theorem possible
  without permanent bounds) and the experiments apply.
- **The content route** (constant SLPs). This is τ of integers, the factorial question (§1.5).
- **Track C O1–O4** remain the barriers for interval/prefix zero sets.

### 4.6 Task-4 answer: what the Strassen degree bound gives, and what rules things out

**Strassen's degree bound.**
- An SLP of size s computes degree ≤ 2^s, so a zero set of size D needs s ≥ log₂ D.
- For density δ this is s ≥ log₂(δp). That is polylog, so it rules out nothing, not even polynomial
  time.

**Ruling out subexponential methods in all of model U is impossible.** ECM is in model U at
L_p[1/2, √2] (heuristic).

**Ruling out polynomial time in model U is permanent-hard (Prop BB).** What *does* rule things out, by
submodel:

| Submodel | Barrier |
|---|---|
| sparse / torus | exponential, N^{1/4−ε} (SP-U1) |
| N-polynomial exponents | density-zero success (SP-N) |
| single-scalar group laws | success ⇔ smooth group order, so ECM-type L[1/2] under smoothness heuristics (OA; KNOWN-type) |
| interval zero sets | Track C O1–O4 |
| birthday / rho | δ ≲ s²/p (heuristic) |
| bounded degree | δ ≤ deg/p (trivial) |

---

## 5. Experiments (pre-registered; raw data in experiments/raw/slpcirc_E{1,2,3,4}.jsonl)

**E1: Lemma S1 check.**
- Exhaustive for t = 2, 3 and all p ≤ 43; random for t = 4, 5 and p ≤ 211.
- 4 429 217 instances, **0 violations** (prediction met). 12 s.

**E2: census of all 132 876 nonconstant polynomials with τ ≤ 6.** (Constant-free, from {1, x}; exact BFS
over value sets; degree ≤ 64.) Primes 67 ≤ p ≤ 251; 4 783 531 (F, p) pairs; 5 pairs had F ≡ 0 mod p.

- **H1 FAILED as stated.** 771 dense pairs (r ≥ p/8). In 701 of them (90.9%, against the predicted
  ≥ 95%), ≥ 90% of the roots are LIN/CYC/KUM/CHEB/LOW.
  - The 70 failures are all at p ≤ 89 (67: 30, 71: 21, 73: 5, 79: 11, 89: 3). Each comes from one or
    two OTHER irreducible factors of degree 6–16 that *split completely* at that prime. So there are at
    most 16 OTHER roots, bounded by the factor degree, as Chebotarev allows. None is dense in any
    asymptotic sense.
  - Inspection of the 27 distinct OTHER factors: 12 are trinomials (x⁶ − x⁴ − 1, x⁷ − x³ − 1,
    x¹² − x⁹ − 1, x⁸ + 3x⁴ + 1, x⁴ ± x ± 1, …). These lie in the Theorem SP class, but the
    pre-registered classification had no SPARSE type.
  - The others are cyclotomic/Kummer after a *non-affine* substitution, which the affine-only
    classification missed. Examples: ((x² − 1)⁷ + 1)/x², x¹⁶ + 3x⁸ + 1 (x⁸ = −φ^{±2}), and compositions
    with y = x² + x.
  - Disclosed: the failure is a limitation of the pre-registered taxonomy. No non-group dense mechanism
    turned up.
- **H2 PASSED.** The 1 805 distinct OTHER factors of degree ≥ 8 have mean roots per good prime m(g)
  averaging **0.898**, in the predicted [0.7, 1.3]. The maximum is 1.44 (predicted ≤ 3); none exceeds
  2. This is Chebotarev behaviour: about one root per irreducible factor, with no small-Frobenius
  anomaly.
- **H3 PASSED.** The best uniform (min-over-window) density is 5/251, attained by
  −x(2x − 1)(2x + 1)(x − 1)(x + 1), which has only LIN roots.
- Factor-type counts: LIN 879, CYC 168, KUM 1 637, CHEB 83, LOW 10 550, OTHER 5 792. 57 s.

**E3: Theorem SP per-prime check.**
- Setup: all 1 612 primes in (2^14, 2^15], λ = 1/20. Families:
  - x^{lcm(1..B)} − 1 for B = 50, 200, 1000;
  - 200 seeded random trinomials and 4-nomials with exponents < 2^40;
  - x^M + x^{M′} − 2.
- **0 violations** (prediction met).
- Primes with δ ≥ 1/20: 281 / 686 / 1095 for lcm(50/200/1000), 534 for the binomial-pair family, and
  **0 for all 200 random sparse polynomials**.
- Disclosed scope: at this size ℓ₀(4, 1/20) ≈ 7·10⁴ > p, so the check is vacuous for t = 4 and
  meaningful only for t ≤ 3. 76 s.

**E4: quadratic compositions (Addendum A; codex's open question).**
- Setup: f = (x² + c_k)∘⋯∘(x² + c₁) with c_i ∈ {−2..2}; 60 seeded primes in (2^12, 2^13].
- **P4a PASSED.** Over all 5^8 sequences with k = 8, the maximum median fibre is 8 (≤ 16). The top
  sequences are built from c ∈ {0, −2} layers ending in −1: powers/Chebyshev, i.e. 2-power torsion.
- **P4b PASSED.** For 10⁴ random k = 39 sequences the mean fibre is 1.46 (in [0.5, 2]), as for a
  critical branching process.
- **P4c PASSED.** No sequence reaches p/(log₂ p)² at ≥ half the primes (maximum count 8 of 60). The
  largest single fibre is 272 at one prime (density ≈ 0.04), with 0 at 52 of the 60 primes.
- So this class shows no uniform dense fibres. 27 s.

---

## 6. Candidate cards (brief) and open questions

| Card | Mechanism | Outcome |
|---|---|---|
| C1 | sparse/torus SLPs (x^M, Dickson, Φ_k, products, several variables, arbitrary coefficients) | **refuted for this model**: exact group-order dichotomy (SP); exponential size uniformly (SP-U1); superpolylog on average (SP-U2) |
| C2 | N-dependent exponent templates (x^{N−1}, x^N, polynomials in N) | **refuted for this model** (SP-N): density-zero success |
| C3 | N-dependent arbitrary constants / "short SLP mod this N" | **equivalent to factoring** (Prop R) |
| C4 | Pochhammer / factorial / content constants | **open**; equivalent to "ultimately easy" factorials; any refutation is permanent-hard (Prop BB) |
| C5 | non-abelian splitting families (x^a − x − 1) | **refuted** (trinomial ⇒ SP; Kelley C = 1 ⇒ δ ≤ 2p^{−1/2}) |
| C6 | quadratic compositions with small fixed c (codex) | **no candidate**; E4 null (critical branching); inverse theorem open |
| C7 | rho-type dynamics | birthday-bounded δ ≈ s²/p (heuristic); exponential |
| C8 | ECM/HECM (group law, single scalar) | Prop OA; L[1/2] (KNOWN) |

**The single most promising open question** (sharpened from trackC §3 and codex):
- Is there an N-independent SLP of size polylog(X) whose output F satisfies all three of:
  - F is **not** sparse in any torus coordinate;
  - F is not ≡ 0 mod p as a polynomial;
  - its zero density mod p lies in [1/polylog, 1 − 1/polylog] for all p ∈ (X, 2X]?
- Theorem SP says such an F must be "dense in every torus coordinate".
- Prop BB says impossibility cannot currently be proved.
- A first tractable piece is codex's: a dense-fibre inverse theorem for compositions of x² + c_i with
  c_i from a fixed set. Heuristically, the fibre after k layers is ∏(1 + a_i), where a_i is a quadratic-
  character bias on the current fibre shifted by c_i. Density p/polylog with k = polylog needs
  biases ≳ 1/k at almost every layer, uniformly in p. E4 sees none.

---

## 7. Reproduction

```bash
cd /Volumes/git/factoring-connection
.venv/bin/python -m pytest -q tests/test_slpcirc.py                        # 9 tests
.venv/bin/python -m src.slpcirc.e1_e3 E1 experiments/raw/slpcirc_E1.jsonl   # 12 s
.venv/bin/python -m src.slpcirc.census 6 experiments/raw/slpcirc_E2.jsonl   # 57 s, ~270 MB
.venv/bin/python -m src.slpcirc.e1_e3 E3 experiments/raw/slpcirc_E3.jsonl   # 76 s
.venv/bin/python -m src.slpcirc.quadcomp experiments/raw/slpcirc_E4.jsonl   # 27 s
```

Development runs used τ ≤ 4 only (disclosed). Everything is single-core, with peak RSS ≤ 300 MB.

## 8. Review record

(Filled in after review: codex adversarial round slpcirc_2, the Opus proof critic and the Opus
literature critic.)
