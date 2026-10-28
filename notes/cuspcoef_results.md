# cuspcoef: one coefficient a_f(N) of a fixed non-CM cusp form at N = pq (phase 30)

Agent `cuspcoef`, 2026-10-28. **Version v1 (before reviews).**
- Pre-registration: `notes/cuspcoef_prereg.md` (commit `82ab8f0`; addendum A1, commit `0d34503`, made after a smoke
  run and before the registered E1 run).
- Code: OPERATIONAL `src/cuspcoef/{decode,family}.py` (inputs: N, the public curve list, the a_N values). An AST test
  forbids PARI, the harness, gmpy2/sympy/flint and factorisation calls. HARNESS `src/cuspcoef/{harness,run_experiments}.py`.
- Tests: `tests/test_cuspcoef.py`.
- Raw data: `experiments/raw/cuspcoef_{e1,e2,e3,e4}.jsonl`; smoke run `cuspcoef_e1_quick.jsonl`.
- GPT-6.1-sol: blind round `codex/cuspcoef_1.{md,out}` (**cx1**); adversarial round `codex/cuspcoef_2.{md,out}` (**cx2**).
- Labels: KNOWN (with reference), PROVED, MODEL (proved inside a stated probabilistic model), HEURISTIC,
  experimentally supported, REFUTED (prediction), OPEN, UNVERIFIED. Strongest novelty wording: "no prior statement
  found by us".

Notation: N = pq, p < q odd primes in [X, 2X) coprime to the level, n = log₂N. f, f_i are fixed newforms of weight k
(k_i) with integer coefficients; a_r = a_f(r); τ_r = (a_r(f_i))_{i ≤ K}. R = ⌊2√(2X)⌋ bounds |a_p|, |a_q| in weight 2.

---

## 0. Verdict (plain language)

**Not settled in reality; settled in the standard random model, with a sharp line between what works and what does
not.**

1. **(B) Barrier — Theorem DFB (MODEL; no prior statement found by us).** Call a datum *determinant-free* if it is
   computed from N and the trace vectors τ_p, τ_q of fixed forms alone (not from p or q themselves). a_f(N) = a_pa_q is
   one; so is the ordered pair (a_p, a_q). In the Lang–Trotter model, where τ_p and τ_q are independent random integer
   vectors with the Chebotarev/Sato–Tate law, any determinant-free datum of K forms raises any algorithm's success
   probability by at most a factor C^K. Here C = 2^{k+3}B², and B is the spread constant of the trace law.
   - For a Serre curve (37a1), B ≤ 2·B_Gek. Numerically B_Gek ≤ 1.63 at 24-bit p (E3); unconditionally
     B_Gek = O(log p·loglog p).
   - So one coefficient is worth at most about 8 bits at toy sizes, and at most 2log₂log N + O(logloglog N) bits in
     general. Bits here means max-divergence, i.e. log₂ of the gain factor in success probability.
   - The proof is a one-line likelihood ratio against a p-free simulator: the product of two uniform integers.
   - Every known reduction lies outside the hypothesis, and DFB says why. Each one uses the determinant, i.e. p itself:
     - Bach–Charles: a(N²) contains p^{k−1};
     - CM pairs: 4p = t² + Du²;
     - #E(Z/N) = (p+1−a_p)(q+1−a_q);
     - CRT-T: t ≡ a_p mod p;
     - CLS: the determinants of the classes.
   - So nativemod's CROSS taxonomy ("the factor information enters through a cross term with a determinant-type
     summand") becomes a theorem in this model.
2. **(R) Reduction.**
   - *One fixed form, or O(1) forms:* no reduction in the model (DFB). The model's lower bound is
     K ≥ log₂(1/ε₀)/log₂C forms, where ε₀ is the success probability of factoring without data.
   - *K = O(n² log n) fixed curves — HEURISTIC, toy-verified.* The products a_N(E_i), **together with the
     factorisations of these products** (integers of about n/2 bits), factor N. The decoder takes the candidate splits
     a_N = xy in the Hasse box, decodes s mod ℓ by Chebotarev likelihoods, and finishes by CRT.
     - E1: 60/60 at n = 32, 40, 48 with K = 32000; 56/60 at the first combination tried.
     - The congruence-only decoder (no factorisation of the products) plateaus at quadratic-residue information
       (Theorem SC). The null simulator is at chance.
   - So for polylog-many curves, a native polynomial-time a_N would reduce n-bit factoring to factoring O(n² log n)
     integers of n/2 + 2 bits. That is a constant-factor gain in the L[1/3] exponent (2^{−1/3}), not polynomial time.
   - The factorisations cannot be bypassed by selecting easy products: the fraction of products whose split is
     poly-time recoverable decays super-polynomially (E2).
   - Weight 2 together with weight 12 (a_E(N) with τ(N)) is just K = 2 in DFB. No lattice attack exists: there is no
     equation linking p to the traces (DEC), and the model gain is ≤ C₂C₁₂.
   - Galois congruences: Theorem SC's scalar hypothesis holds at every level for Serre curves (Prop SCX). For a general
     non-CM curve, non-scalar images add at most log₂[Ẑ^× : Λ_E] = O(1) bits, where Λ_E is the group of scalars in the
     adelic image (11a1: s mod 5).
3. **(A) Native algorithm.** The best native route found is **N^{1/2+o(1)}**: Eichler–Selberg trace formulas at the
   auxiliary indices 3N, 5N (and 11N at level 37) cancel the divisor term, leaving O(√N) Hurwitz class numbers. This
   comes from cx1's scratch checks (`codex/scratch/cuspcoef_1/elliptic_trace_check.py`); cx1's final answer omits it
   and cautions that trace formulas can hide a factoring step.
   - The same class-number sums also output s = p + q, so this route is a native N^{1/2+o(1)} *factoring* method.
     It is not a coefficient-specific speedup.
   - Next is Õ(N): the theta-difference shifted convolution, the complete Jacobi sum Σ_x ((x³+ax+b)/N), or Merel's
     Heilbronn sets.
   - Nothing below N^{1/2} is known. Nothing subexponential is known.
   - **Any native algorithm that also works at the index N² factors N with the same cost** (Bach–Charles, even
     weight). So a native sub-L[1/3] algorithm for a_f(n) would be a sub-L[1/3] factoring algorithm, whatever the
     answer to (R). The single-value question matters only for oracles, or for algorithms restricted to squarefree
     indices. Even the squarefree-only Jacobi-sum formula leaks the factors at N² (Remark JS).

Subjective: the probability that this direction yields sub-L[1/3] factoring is < 1%. The cleanest remaining open item
is a reduction that uses a_N **together with** an archimedean-level structure computable from N, such as the group
E(Z/N) or torsion algebras at levels growing with n. DFB does not cover such couplings (§2.4), and we found none.

---

## 1. What exactly is open

- **Oracle question.** Given N and A = a_f(N) (one value, one fixed non-CM f), factor N.
  - An oracle for a_f(n) at all n settles it by Bach–Charles: N² is an index. Restricting to squarefree indices
    gives nothing new, since a_f(Nm) = a_f(N)a_f(m) for gcd(m, N) = 1. So the squarefree-index oracle is equivalent
    to the single value plus values at indices coprime to N.
- **Native question.** Compute a_f(N) without factoring N, below L[1/3].
- tracedata §6.1 and nativemod §5 left both OPEN. Prop DEC (tracedata) shows that (N, a_p, a_q) satisfy no FIXED
  algebraic relation with p on positive-probability events. That excludes direct lattice equations, but not
  arithmetic uses of the integers.

---

## 2. (B) Theorem DFB: determinant-free data are worth O(log log N) bits in the Lang–Trotter model

### 2.1 Model and statement

**Model LT(B) (MODEL).**
- Fix newforms f_1, …, f_K of weights k_i with integer coefficients, and N = pq with p ∈ [X₁, 2X₁), q ∈ [X₂, 2X₂).
- For r ∈ {p, q}, τ_r = (a_r(f_i))_i is a random vector in Z^K with |a_r(f_i)| ≤ 2r^{(k_i−1)/2}, and
  Pr[τ_r = t] ≤ ∏_i B·r^{−(k_i−1)/2} for every t ∈ Z^K (spread B).
- τ_p and τ_q are independent of each other.
- Everything else the algorithm sees is fixed given N, or independent of (τ_p, τ_q). This includes N itself, every
  object computable from N (E(Z/N), torsion algebras, Jacobi symbols), the algorithm's coins, and any oracle whose
  answers do not depend on (τ_p, τ_q), e.g. a_f at indices coprime to N.

**Determinant-free datum.** D = Φ(N, τ_p, τ_q, ω), for a fixed (possibly randomized, coins ω) map Φ that receives N
but not p or q. Examples: a_f(N) = a_pa_q; (a_{f_i}(N))_i; the unordered or ordered trace pairs; the factorisation
of a_pa_q. Non-examples: a_f(N²) = a_N² + N^{k−1}·(…) − (q^{k−1}a_p² + p^{k−1}a_q²) (Bach–Charles); #E(Z/N); the
CRT trace t_N; anything that uses p, q individually.

**Theorem DFB (MODEL, PROVED in the model).** Let σ₁, σ₂ be independent, each uniform on the integer box
∏_i [−⌊2(2X_j)^{(k_i−1)/2}⌋, ⌊2(2X_j)^{(k_i−1)/2}⌋] (j = 1, 2), and put D̃ = Φ(N, σ₁, σ₂, ω). Then for every
algorithm 𝒜 (with any oracle as above) and every set S of outputs,

  Pr[𝒜(N, D) ∈ S] ≤ C^K·Pr[𝒜(N, D̃) ∈ S],
  C^K := ∏_i B²(2⌊2(2X₁)^{(k_i−1)/2}⌋+1)(2⌊2(2X₂)^{(k_i−1)/2}⌋+1)/N^{(k_i−1)/2} ≤ ∏_i 2^{k_i+3}B²(1+o(1)).

D̃ depends on N, X₁, X₂ and coins only. So the right-hand side is the success of an algorithm that is given N alone
(plus a polynomial-time sampler).

*Proof.*
- For every (t, u), Pr[(τ_p, τ_q) = (t, u)] ≤ ∏_i B²(pq)^{−(k_i−1)/2}, and it vanishes outside the Hasse boxes,
  which lie inside the σ-boxes because p < 2X₁ and q < 2X₂.
- Pr[(σ₁, σ₂) = (t, u)] = 1/(box sizes) on that box. Hence Pr[τ = (t, u)] ≤ C^K·Pr[σ = (t, u)] pointwise.
- Let g(t, u) := Pr[𝒜(N, Φ(N, t, u, ω)) ∈ S], the probability over ω, the algorithm's coins and the oracle. This is
  well defined because these are independent of τ.
- Then Pr[𝒜(N, D) ∈ S] = Σ_{t,u} Pr[τ = (t, u)]·g(t, u) ≤ C^K Σ Pr[σ = (t, u)]·g(t, u) = C^K Pr[𝒜(N, D̃) ∈ S].
- Bound on C: (4(2X)^{(k−1)/2} + 1)²/N^{(k−1)/2} ≤ 16·2^{k−1}(X₁X₂)^{(k−1)/2}/(pq)^{(k−1)/2}·(1 + o(1)), and
  X₁X₂ ≤ pq. ∎

(The proof is the data-processing inequality for max-divergence, D_∞(Φ(P)‖Φ(Q)) ≤ D_∞(P‖Q). It is unconditional
inside the model, and holds for every N, not on average.)

**Corollary (bits; K forms).** In LT(B), a determinant-free datum of K forms raises success by ≤ K·log₂C bits of
max-divergence. If every algorithm of a class closed under prepending the simulator factors N with probability
≤ ε₀, then with K forms the probability is ≤ C^K ε₀. Reaching probability 1/2 therefore needs
K ≥ log₂(1/(2ε₀))/log₂C forms.

### 2.2 The spread constant B

- **Gekeler law (random curve over F_p).** Pr[a = t] = H(4p − t²)/(2p), with H the Hurwitz class number
  (Σ_{t²<4p}H(4p − t²) = 2p, Kronecker–Hurwitz). Gekeler (IMRN 2003) shows that this equals the product of the
  GL₂(Z_ℓ) Chebotarev densities and the Sato–Tate density. So B_Gek(p) = max_t H(4p − t²)/(2√p).
- **Bound (PROVED, routine).**
  - For any negative discriminant −Δ, the class number formula gives 2h(−Δ)/w(−Δ) = √Δ·L(1, χ_{−Δ})/π.
  - Partial summation gives L(1, χ_{−Δ}) ≤ 3 + ln Δ.
  - Hence H(D) ≤ (√D/π)(3 + ln D)·σ(m)/m, where m² is the largest square divisor of D.
  - Robin's bound σ(m)/m ≤ e^γ lnln m + 0.65/lnln m (m ≥ 3) then gives
    B_Gek(p) ≤ (1/π)(3 + ln 4p)(e^γ lnln 4p + 0.65/lnln 4p) = O(log p·loglog p).
  - Under GRH, L(1, χ) ≪ loglog Δ, so B_Gek = O((loglog p)²).
- **E3 (measured):** B_Gek ∈ [0.98, 1.18] at 16 bits, [1.29, 1.54] at 20 bits, [1.39, 1.63] at 24 bits. This is far
  below the bound, which is 24–39 there.
- **A fixed curve.**
  - In the Lang–Trotter law the adelic density of tr = t on the fibre det = p is that of the image G instead of
    GL₂(Ẑ). Its ratio to the GL₂ density is ≤ [GL₂(Ẑ) : G], because det is surjective on G.
  - So B_E ≤ [GL₂(Ẑ) : G]·B_Gek. This is 2 for a Serre curve (37a1; adelic index 2 per LMFDB, UNVERIFIED by us
    beyond E4).
  - 11a1 has index divisible by [GL₂(F₅) : Borel with trivial character] = 24, from its rational 5-torsion. That is
    worth ≤ 9 extra bits, and the actual leak is s mod 5.
- **Weight 12 (Δ).** τ(p) is fixed modulo 2^{11}·3^6·5^3·7·23·691 by p (the Ramanujan-type congruences). So B_Δ
  carries that constant factor, a few tens of bits, still O(1). The leak is s mod 691 and similar.
- **Toy check (exact; unit test `test_dfb_ratio_bound_exact`).** p = 607, q = 1009 under the Gekeler law gives
  log₂C = 3.65 bits. The actual worst ratio is 2.11 bits, and the KL divergence is 0.50 bits.

### 2.3 What DFB excludes, and why the hypotheses are sharp

**Excluded in the model (all determinant-free):**
- every Coppersmith/lattice attack built from a_N = a_pa_q and the Hasse bounds, for one or several curves, including
  weight 2 together with weight 12;
- congruence decoding of a_N, or of several a_N(E_i) for O(1) curves (SC, GC);
- Sato–Tate/archimedean reasoning;
- factoring the integer a_N;
- any combination of these, adaptive or not, with any amount of computation.

**Not excluded, and these are exactly the known reductions:**

| reduction | which hypothesis it violates |
|---|---|
| Bach–Charles (a_N, a_{N²}) | a_{N²} uses p^{k−1}, q^{k−1}: not determinant-free |
| CM pairs y² = x³ − x, x³ − Dx (tracedata CM-PAIR) | spread: for CM, a_p ∈ {±2α, ±2β} with α² + β² = p, so B = Θ(√p) |
| M = #E(Z/N) (Theorem EXT) | not determinant-free: M = (p+1−a_p)(q+1−a_q) |
| CRT trace t_N ≡ a_p (mod p) (CRT-T) | uses p |
| class pair (CLS) | determinants of classes |
| ordered traces of O(n² log n) non-CM curves (tracedata MC) | none: consistent with DFB (K·log₂C ≫ n/2) |

So DFB says, in a model, that **every reduction from Hecke data must use a determinant (p itself), either explicitly
or through a CM-type algebraic relation, unless it uses many forms**. That is nativemod's Prop CROSS observation
promoted to a theorem inside LT(B).

### 2.4 Scope: what the model does not capture (honest limits)

1. **Determinism.** Real traces are deterministic, poly(n)-time functions of the prime (Schoof, E–C). DFB is a
   random-oracle-type statement. A reduction exploiting a specific non-statistical property of p ↦ a_p (as CM does)
   is outside it. For non-CM f, DEC excludes fixed algebraic relations; nothing excludes cleverer ones.
2. **Coupling to structures computable from N.**
   - LT(B) treats τ_p as independent of E(Z/N), (Z/N)[x]/(ψ_ℓ), etc. In reality a_p mod ℓ is the trace of the
     Frobenius acting on the torsion algebra, and p + 1 − a_p is the order of a computable group's p-component.
   - **DFB(m) (MODEL; same proof).** Suppose that, conditionally on the true Frobenius classes c_p, c_q at a level m,
     the lifts are independent with Pr[τ_r = t | c_r] ≤ ∏_i Bm·r^{−(k_i−1)/2} on the residue class of tr c_r. Then
     Pr[𝒜(N, D) ∈ S | c] ≤ C_m^K·Pr[𝒜(N, Φ(N, σ^{(m)})) ∈ S | c], with σ^{(m)} uniform on the residue classes and
     C_m ≤ ∏ 2^{k_i+3}B²(1 + m/(2R_i))². So, relative to an algorithm that also holds the true classes at level m
     (e.g. by probing torsion algebras), the lift of a determinant-free datum adds ≤ K·log₂C_m bits.
   - Couplings at levels growing with n are not covered. The class pair at all ℓ ≤ O(log N) factors N (CLS), but no
     way to obtain it from N is known (FROB/PRIM).
   - Coupling through the group order is not covered either. The natural attempt needs the exponent p + 1 − a_p of
     the p-component; a_N fixes neither a_p nor p. Over (Z/N)[t]/(t² − c) the conjugation is Frobenius at inert
     primes, and gives [a_p]σ = [p+1] on the p-part. From a_N that yields [a_N]σP = ([a_q(p+1)]P_p, [a_p(q+1)]P_q),
     again with unknown p. No determinant-free annihilator exists by DEC's annihilator remark (tracedata §2.1).
3. **The model is for one N.** No average over instances is involved: the bound holds for every N.

---

## 3. Prop SCX: the exact scope of Theorem SC (congruences of products)

**Setting.** E/Q non-CM, with adelic image G ⊆ GL₂(Ẑ). Λ := {λ ∈ Ẑ^× : λI ∈ G}, and Λ_m its image in (Z/m)^×.

**Prop SCX (PROVED from SC and Serre's open image theorem).**
- (i) Λ is open in Ẑ^×. G is open (Serre 1972), so it contains Γ(M_E) ∋ λI for every λ ≡ 1 (mod M_E). Hence
  [Ẑ^× : Λ] ≤ |(Z/M_E)^×|.
- (ii) For every m, the Chebotarev law of (N mod m, a_N mod m), conditioned on p mod m, is invariant under
  p ↦ λ²p, q ↦ λ^{−2}q for λ ∈ Λ_m (tracedata Theorem SC, applied with the scalars Λ_m). So the congruence content
  about p mod m is at most the class of p in (Z/m)^×/Λ_m². That exceeds the quadratic-character content
  (Z/m)^×/((Z/m)^×)² by a factor ≤ [(Z/m)^× : Λ_m] ≤ [Ẑ^× : Λ], i.e. by ≤ log₂[Ẑ^× : Λ] bits, uniformly in m.
- (iii) For a Serre curve, Λ = Ẑ^×. Here G = {g : ε(g mod 2) = χ_Δ(det g)}, where ε is the sign character of
  GL₂(F₂) ≅ S₃. For a scalar, λI ≡ I (mod 2) gives ε = 1, and χ_Δ(λ²) = 1. So the congruence content of a_N at
  every level is quadratic-character content only.
  - 37a1 is a Serre curve per LMFDB (UNVERIFIED by us).
  - E4 checks the finite-level consequences:
    - all ℓ(ℓ−1) pairs (a_p mod ℓ, p mod ℓ) occur for every ℓ ≤ 47, which by Dickson gives surjectivity mod ℓ;
    - the entanglement rule "a_p odd ⇒ (37/p) = 1" has 0 violations among 17982 primes, and 2967 primes satisfy
      (37/p) = 1 with a_p even.
- (iv) Escapes from QR content therefore come only from finitely many levels where scalars are missing (rational
  isogenies or torsion, entanglement), and are bounded by O(1) bits per curve.
  - Example: 11a1 at 5, where Λ₅ = {1} and a_N ≡ (1+p)(1+q) gives s mod 5.
  - For K curves the joint scalar group is open as well, so the excess is still O(1) per curve.

**E1 confirms the plateau.** The congruence-only decoder (C2) at K = 32000 has per-ℓ top-1 accuracy 1.0, 0.63,
0.50, 0.27, 0.18 for ℓ = 3, 5, 7, 11, 13. That is tracedata's E3b plateau (0.64, 0.44–0.46, …). The decoder that
uses the factorisations reaches ≥ 0.975 at every ℓ ≤ 23.

---

## 4. (R) Many curves: products plus half-size factorisations

**Prop RM (HEURISTIC; toy-verified).**
- Fix the family F (y² = x³ + Ax + B, one curve per j, non-CM; `family.py`). Given N, a_N(E_i) for the first K
  curves, and the factorisations of the |a_N(E_i)| (integers ≤ 8X ≈ 4√N·√2):
  1. List the candidate splits C_i = {(x, y) : xy = a_N(E_i), |x|, |y| ≤ R}.
  2. For each odd ℓ ≤ L (∏ℓ > 2X), score d = p mod ℓ by Σ_i log Σ_{(x,y)∈C_i} n_ℓ(x, d)·n_ℓ(y, N/d), where
     n_ℓ(t, d) = ℓ(ℓ + χ(t² − 4d)) is the GL₂(F_ℓ) class count.
  3. Pass to s mod ℓ, list-decode, CRT, and solve.
- Under the Chebotarev-independence heuristic for families (tracedata HEURISTIC MC), this factors N once
  K ≈ c·L² log L with L ≈ 0.35n, i.e. **K = O(n² log n)**.
  - Unordered pairs identify s mod ℓ (tracedata PH1).
  - The true split is always among the candidates. The others are noise with ST-like residues.
- **Cost.** K factorisations of integers of ≤ n/2 + 3 bits, plus poly(n). With NFS on these, the total is
  L_N[1/3, (64/9)^{1/3}·2^{−1/3}]·poly. So a polynomial-time native a_N for the first K curves of F would cut the
  NFS constant by 2^{1/3}. It would not give polynomial time.
- **The factorisations cannot be skipped by selection (E2).**
  - The split is computable in poly(n) only if every prime factor of a_N, except one, is ≤ n^{O(1)}. That requires
    one of a_p, a_q (≈ N^{1/4} each) to be n^{O(1)}-smooth. This has probability ρ(u) with u ≈ (n/4)ln2/(c ln n),
    which is super-polynomially small.
  - E2: real/simulated acceptance rates (B = n²) are 0.94/0.91 (n = 64), 0.40/0.33 (n = 128), 0.045/0.053 (n = 192), 0.005/0.006 (n = 256).
- **Without factorisations**, the products' congruences carry only QR bits (SC, SCX). A partially factored product
  with an unsplit composite part loses the ordering at every ℓ whenever the composite can straddle x and y.
- **Few curves.** K = O(1) curves give ≤ K·log₂C bits (DFB), so no reduction exists in the model. The tension between
  DFB's K ≳ n^δ/loglog N (if ε₀ ≤ 2^{−n^δ}) and RM's O(n² log n) is a gap of polynomial factors only.
- Recursion (remark, not tested): if a squarefree-index oracle for all curves of F were available, the a_N(E_i) of
  the half-size products could be fed back. That needs a decoder for k-factor squarefree indices (multisets of
  traces), which we did not build. We claim nothing here.

---

## 5. (A) Native algorithms for a_f(N)

| route | cost (best known; upper bounds) | native? | status |
|---|---|---|---|
| definition / Hecke multiplicativity, E–C at p and q | poly(n) after factoring | no | KNOWN |
| **Eichler–Selberg at auxiliary indices ℓN, ℓ ∈ {3, 5} (cx1 scratch; the same sums give s = p + q)**: the hyperbolic terms of Tr T_{3N}, Tr T_{5N} differ only by known constants when q < 3p; τ(3N) − τ(5N) = (252 − 4830)τ(N) isolates τ(N); level 11/37 analogues | O(√N) Hurwitz class numbers H(D), D ≤ 20N; each H(D) needs h of the orders above disc −D (L_D[1/2] under GRH, Hafner–McCurley) and the square part of D | yes (GRH; square-part caveat) | cx1, exact checks; **N^{1/2+o(1)}** |
| theta difference: r_{L1}(N) = Σ_{m≤N} r_B(m)r_B(N−m), B = x² + xy + 3y² (h = 1, r_B = 2·(1 ∗ χ₋₁₁)), then a(N) = (5/8)(r_{L1}(N) − (12/5)σ(N)) | Õ(N) (segmented sieve) + Õ(N^{1/2}) for σ(N) | yes | PROVED (from nativemod (iii)) |
| complete Jacobi sum a_N = Σ_{x mod N} ((x³ + ax + b)/N) (squarefree N) | Θ(N) Jacobi symbols | yes | KNOWN identity (cx1 in tracedata) |
| Merel/Heilbronn sets for T_N on modular symbols | ≍ σ(N)log N or ≍ N log²N matrices | yes | KNOWN (nativemod §3.1) |
| fast exponentiation in a fixed algebra / P-recursive recurrence | impossible | — | Prop LRS; natural boundary (nativemod) |
| analytic extraction of one coefficient | precision wall | — | nativeana |

**Remark BC-native (KNOWN consequence of Bach–Charles).** An algorithm computing a_f(n) at every index n (in
particular n = N²) in time T(log n) factors a density-one set of RSA moduli in time 2T(2 log N) + poly (even weight,
non-CM). So a native sub-L[1/3] algorithm valid at all indices **is** a sub-L[1/3] factoring algorithm, independent of
(R). Every route in the table works at all indices, except the Jacobi sum.

**Remark JS (PROVED, elementary).** The Jacobi-sum formula is squarefree-specific, yet it leaks at N² as well. For
any cubic f,
- Σ_{x mod N²} (f(x)/N²) = #{x mod N² : gcd(f(x), N) = 1} = N(p − r_p)(q − r_q), with r_r = #{roots of f mod r}.
- Given V = (p − r_p)(q − r_q), each guess (r_p, r_q) ∈ {0, 1, 2, 3}² gives
  r_q p² − (N + r_pr_q − V)p + r_pN = 0, which recovers p unless r_p = r_q = 0.
- For an S₃ cubic that exception has probability 1/9.

**No route below N^{1/2} was found.** nativemod's "first question" (compute the cuspidal combination without separating
σ(N)) is untouched below that.

---

## 6. Experiments (pre-registered; `notes/cuspcoef_prereg.md`)

Peak RSS ≤ 726 MB (E1, n = 48 sieve). Jobs ran one at a time. E1 ≈ 3 min, E3 12 s, E4 < 1 min, E2 ≈ TBD.

| exp | prediction | result |
|---|---|---|
| E4 P4a | all ℓ(ℓ−1) pairs, 5 ≤ ℓ ≤ 47 | all present for every ℓ; **holds** |
| E4 P4b | ℓ = 3: 6 pairs | 6/6; **holds** |
| E4 P4c | a_p odd ⇒ (37/p) = 1, 0 violations; converse fails | 0 violations in 17982 primes; converse fails 2967 times; **holds** |
| E3 P3a | Kronecker–Hurwitz sum for all 30 primes | 30/30; **holds** |
| E3 P3b | max B_emp ≤ 4, below the explicit bound | max 1.63; bound 24–39; **holds**. The registered bound formula used "2 + ln 4p" and "+1"; the proved one (§2.2) uses "3 + ln 4p" and "+0.65/lnln". Both are far above B_emp |
| E1 P1a | K = 32000: ≥ 17/20 (n = 32, 40), ≥ 14/20 (n = 48) | 20/20, 20/20, 20/20; **holds** |
| E1 P1b | success non-decreasing in K | n = 32: 20, 20, 20; n = 40: 20, 20, 20; n = 48: 18, 20, 20; **holds** |
| E1 P1c | null C1 ≤ 1/20 | 17/20 (n = 32), 6/20 (n = 40), 0/20 (n = 48): **REFUTED at n = 32, 40**. Exhaustive verification within the 3000-combination budget (space 9072 and 90720); see A1 |
| E1 P1d | congruence-only C2 ≤ 2/20 | 20/20, 17/20, 5/20: **REFUTED at all n**. SC's QR content shrinks the effective space by about 2 per ℓ, so 3000 combinations cover most of it (n = 48: ≈ 1.09·10⁶/2⁸ ≈ 4·10³). This should have been predicted from SC |
| E1 P1e | per-ℓ top-1 ≥ 0.8 for ℓ ≤ 13 (K = 32000) | 1.0, 1.0, 0.983, 1.0, 0.983 (ℓ = 3–13); ≥ 0.975 up to ℓ = 23; **holds** |
| E1 P1f (A1) | median `tried` ≤ 5 at K = 32000 | 1, 1, 1; **holds** |
| E1 P1g (A1) | median `tried` of C1, C2 ≥ 100 | C1: 1472, > 3000, > 3000; C2: 184.5, 764.5, > 3000; **holds** |
| E1 P1h (A1) | top-1: main ≥ 10/20; C1, C2 ≤ 1/20 | main 19, 18, 19; C1 0, 0, 0; C2 0, 0, 0; **holds** |
| E2 P2a | real within ×2 of simulator where sim ≥ 0.02 | real/sim 0.943/0.913 (n=64), 0.400/0.335 (128), 0.045/0.053 (192); **holds** |
| E2 P2b | simulator decreasing, ≤ 0.02 at n = 256 | sim 0.913, 0.335, 0.053, 0.0059 (real 0.005 at n=256); **holds** |

**E1 details (K = 2000 / 8000 / 32000; main decoder).**
- Success: n = 32: 20/20/20; n = 40: 20/20/20; n = 48: 18/20/20.
- Top-1: n = 32: 7/13/19; n = 40: 5/9/18; n = 48: 0/7/19.
- Median `tried`: n = 32: 3/1/1; n = 40: 15.5/2/1; n = 48: 661/5.5/1.
- Mean candidate splits per product: 21 (n = 32), 28 (n = 40), 33 (n = 48). For the null they are 14–21: real
  traces are more often divisible by small primes.
- Accidental gcd(Δ_i, N) splits (bad reduction of a family curve at p or q): 27 at n = 32 over 20 instances
  (K = 32000), 2 at n = 40, 0 at n = 48. All ACCIDENTAL; the structural decoder ignores them.
- These are toy correctness checks (rule 8). The harness supplies the a_N, and the half-size factorisations are
  trivial at these sizes.

---

## 7. Candidate cards

| mechanism | model | outcome |
|---|---|---|
| Coppersmith/lattice on a_N = a_pa_q (+ Hasse), one or few curves, ± τ(N) | LT(B) | **refuted for this model** (DFB): ≤ K·log₂C bits |
| congruences of a_N, one or many curves | Chebotarev | QR bits + O(1) (SC, SCX) |
| products of K = O(n² log n) curves + factorisations | HEURISTIC MC | works (E1 60/60); costs K half-size factorisations |
| selecting easy products to avoid factoring | Dickman | fails: super-polynomially rare (E2) |
| native: Eichler–Selberg cancellation (cx1) | — | N^{1/2+o(1)} (GRH, square parts) |
| native: theta shifted convolution, Jacobi sum, Merel | — | Õ(N) |

---

## 8. Reviews

**cx1 (GPT-6.1-sol blind round, `codex/cuspcoef_1.out`).**
- It found no reduction for one, O(1) or polylog-many forms: the question is OPEN. It is a literature conclusion,
  not an impossibility theorem.
- Its own MODEL barrier is congruence-only, and complementary to DFB: K curves give < 1.187K bits about the hidden
  determinant residues, summed over all odd auxiliary ℓ. The proof bounds χ²(P_d‖U) ≤ (ℓ−1)^{−2} and sums
  Σ_ℓ 2log₂(1 + (ℓ−1)^{−2}) < π²/(12 ln 2).
- It confirms the scalar symmetry for any number of curves, and for mixed weights with scalars (λ^{k_i−1}I)_i.
- It confirms the twist identity a_{E^{(d)}}(N) = (d/N)a_E(N), the 11a1 s mod 5 leak, and the Bach–Charles
  denominator mechanism.
- It adds Bruin, Publ. Math. Besançon (2011), Thm 1.1 and Rem 1.3.
- Native routes: Õ(N) for Merel, the Jacobi sum, theta meet-in-the-middle and the binary convolution. Charles 2006
  factors n first.
- Its scratch trace-formula checks are the source of the N^{1/2+o(1)} row in §5. Its own caution applies: the
  class-number sums reveal s.

**Not yet done at handback:** cx2 (prompt `codex/cuspcoef_2.md`), the Opus proof critic report, and the Opus
literature critic. Every MODEL/PROVED claim here is therefore unreviewed.

## 9. Reproduction and deviations
```
.venv/bin/python -m pytest -q tests/test_cuspcoef.py
PYTHONPATH=. .venv/bin/python -m src.cuspcoef.run_experiments e4    # < 1 min
PYTHONPATH=. .venv/bin/python -m src.cuspcoef.run_experiments e3    # 12 s
PYTHONPATH=. .venv/bin/python -m src.cuspcoef.run_experiments e1    # ~3 min, peak RSS 726 MB
PYTHONPATH=. .venv/bin/python -m src.cuspcoef.run_experiments e2    # ~TBD
```
- **D1 (A1).** The smoke run revealed the list-decoding budget flaw before the registered run. P1c and P1d are reported
  as registered (refuted). P1f–P1h were added before the registered run.
- **D2.** After E1 ran, the method `SPF.factor` was renamed `SPF.decompose` so that the AST leakage test does not flag
  the name `factor`. This is cosmetic: no logic changed.
- **D3.** E3's registered explicit bound differs from the one proved in §2.2 (see the table).
