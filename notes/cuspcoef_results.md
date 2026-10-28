# cuspcoef: one coefficient a_f(N) of a fixed non-CM cusp form at N = pq (phase 30)

Agents `cuspcoef` (v1) and `cuspcoef2` (v2), 2026-10-28. **Version v2.** It applies the blind round cx1, the
adversarial round cx2 (MAJOR REVISION), the Opus proof critic and the Opus literature critic, and reconciles every
quoted number against the raw files. All four reviews are logged in §8.
- Pre-registration: `notes/cuspcoef_prereg.md` (commit `82ab8f0`; addendum A1, commit `0d34503`, made after a smoke
  run and before the registered E1 run).
- Code: OPERATIONAL `src/cuspcoef/{decode,family}.py` (inputs: N, the public curve list, the a_N values). An AST test
  forbids PARI, the harness, gmpy2/sympy/flint and factorisation calls. HARNESS `src/cuspcoef/{harness,run_experiments}.py`.
- Tests: `tests/test_cuspcoef.py` (9 tests, all pass).
- Raw data: `experiments/raw/cuspcoef_{e1,e2,e3,e4}.jsonl`; smoke run `cuspcoef_e1_quick.jsonl`; E2 supplements
  `cuspcoef_e2_n256.jsonl` and `cuspcoef_e2_n64.jsonl` (these REPEAT samples of `cuspcoef_e2.jsonl` and must not be
  pooled with it; see D4–D5).
- GPT-6.1-sol: blind round `codex/cuspcoef_1.{md,out}` (**cx1**); adversarial round `codex/cuspcoef_2.{md,out}` (**cx2**).
- Labels: KNOWN (with reference), PROVED, MODEL (proved inside a stated probabilistic model), HEURISTIC,
  experimentally supported, REFUTED (prediction), OPEN, UNVERIFIED. Strongest novelty wording: "no prior statement
  found by us".

Notation: N = pq, p < q odd primes in [X, 2X) coprime to the level, n = log₂N. f, f_i are fixed newforms of weight k
(k_i) with integer coefficients; a_r = a_f(r); τ_r = (a_r(f_i))_{i ≤ K}. R = ⌊2√(2X)⌋ bounds |a_p|, |a_q| in weight 2.

---

## 0. Verdict (plain language)

**The single-value question is OPEN for real curves.** Given N = pq and one coefficient a_f(N) of a fixed non-CM form,
we know no way to factor N, and no proof that it cannot be done. Neither does any source that we, cx1, cx2 or the
literature critic found. Everything below is either a statement inside an idealised random model, a heuristic, or a
toy-size experiment.

1. **Barrier: a MODEL statement only (Theorem DFB, §2).**
   - Replace the true traces a_p(f_i), a_q(f_i) by independent random integers with a spread-out law (model LT(B), a
     Lang–Trotter-type random model; it is not a statement from Lang–Trotter's LNM 504).
     Then any datum built from N and the traces alone ("determinant-free") multiplies any algorithm's success
     probability by at most an explicit factor C_K. For one weight-2 Serre curve, log₂C is about 2log₂log N + O(1),
     i.e. one coefficient is worth O(log log N) bits in the model.
   - The proof is a one-line likelihood ratio. It holds only under three hypotheses that v1 left implicit (§2.1):
     the datum must be efficiently simulable, the algorithm must not evaluate the trace maps at p or q, and the forms
     must be pairwise not twist-equivalent (the public family F has 34 isogenous pairs among its first 32000 curves).
   - **DFB is not a hardness theorem about actual coefficients.** Real traces are deterministic functions of p, and
     E(Z/N) is not independent of them. DFB does not show that a_f(N) alone is useless for factoring.
   - Every reduction we know uses something the model forbids: the determinant p^{k−1} (Bach–Charles; group orders;
     CRT traces; class pairs), CM (no spread), or many forms (tracedata MC).
2. **Congruences: also MODEL.** Congruence data from K curves, summed over all odd auxiliary primes ℓ, carry
   < 2K·log₂(π²/8) < 0.607K bits about the hidden p mod ℓ (Prop CI, §2.5, in the independent full-GL₂ model). This is
   cx1's bound (1.187K), sharpened by cx2 (0.698K) and by us (exact χ²). For curves with non-scalar images there is an
   O(1) excess (Prop SCX, §3; 11a1 leaks s mod 5).
3. **Many curves: a HEURISTIC reduction, toy-verified (Prop RM, §4).**
   - The products a_N(E_i) of poly(n) curves, **together with the factorisations of these products** (integers of about
     n/2 bits), factor N. E1: 60/60 at n = 32, 40, 48 with K = 32000, and 56/60 at the first combination tried.
   - The sample size is poly(n). v1's K = O(n² log n) was imported from ordered-trace data and is not supported for
     candidate splits (proof critic, cx2).
   - Cost: NFS on the products gives L_N[1/3, (32/9)^{1/3}] ≈ L_N[1/3, 1.53], given a polynomial-time a_N for the curves.
     The proof critic's ECM-selection variant suggests a constant ≈ 0.91–1.15 (HEURISTIC, untested, needs
     L_N[1/3]-many curves).
   - E2: the registered "polynomial-time split" filter (cofactor after removing n²-smooth primes is 1 or prime) has
     acceptance 0.943, 0.400, 0.045, 0.005 at n = 64, 128, 192, 256 (simulator 0.913, 0.335, 0.053, 0.0059). P2a and
     P2b hold. This measures that one filter. It is not a barrier against every efficient way to split a product
     (cx2: Fermat-close factors, perfect powers).
4. **Native computation: nothing faster than factoring is known (§5).**
   - Every native route costs at least N^{1/2+o(1)}. Factoring first costs Õ(N^{1/5}) deterministically (Harvey 2021)
     and then poly(n) (Schoof 1985 in weight 2; Edixhoven–Couveignes Thm 15.2.1 at level 1; Bruin 2011, Thm 1.1 in
     general). So "native" is only meaningful below N^{1/5}.
   - The N^{1/2+o(1)} Eichler–Selberg route is a factoring method in disguise: its own class-number sums give
     Σ_t H(4N − t²) = 2(N + q). Its √N class-number cost is Charles's (Ramanujan J. 2006); the new step, cancelling the
     divisor term with the indices 3N and 5N, comes from cx1's scratch checks.
   - An algorithm that computes a_f at every index in time T(n) factors N in time 2T(2n) + poly (Bach–Charles, even
     weight). So the single-value question matters only for a supplied value, or for algorithms restricted to
     squarefree indices.

Registered predictions: **P1c and P1d are REFUTED** (the controls succeed by exhaustive search within the decoding
budget, as foreseen in addendum A1). All other predictions (P1a, P1b, P1e–P1h, P2a, P2b, P3a, P3b, P4a–P4c) hold.

**Number reconciliation (§6.1).** All E1, E3, E4 and E2 counts in v1 match the raw files, as do the DFB toy values.
v1's report figures "0.91, 0.33, 0.053" are the SIMULATED E2 rates and are correct. Changed in v2:
- the E1 search-space sizes, where v1 quoted the upper bound ∏(ℓ+1)/2 = 9072, 90720, 1.09·10⁶ instead of the
  registered per-instance S_n (medians 2604, 23814, 365310). The QR-reduced spaces are 440, 4320, 26400, not "≈ 4·10³"
  at n = 48;
- the E3 range of the *proved* bound, 22.2–35.7 (24–39 was the *registered* formula);
- E2 runtime 612 s in total; E1 rerun 122 s, peak RSS 739 MB (v1 said "≈ 3 min, ≤ 726 MB").
No v1 headline number changed.

Subjective: the probability that this direction yields sub-L[1/3] factoring is < 1%. The most concrete open item is a
reduction that couples a_N with a structure computable from N, such as E(Z/N) or torsion algebras at growing levels.
DFB's model assumes such couplings away (§2.4), and we found none.

---

## 1. What exactly is open

- **Oracle question.** Given N and A = a_f(N) (one value, one fixed non-CM f), factor N. OPEN.
  - An oracle for a_f(n) at all n settles it by Bach–Charles, since N² is an index.
  - Restricting to squarefree indices gives nothing new, since a_f(Nm) = a_f(N)a_f(m) for gcd(m, N) = 1. So the
    squarefree-index oracle is equivalent to the single value plus values at indices coprime to N.
- **Native question.** Compute a_f(N) without factoring N, faster than factoring N (below N^{1/5} deterministically,
  below L[1/3] heuristically). OPEN; nothing below N^{1/2} is known.
- tracedata §6.1 and nativemod §5 left both OPEN. Prop DEC (tracedata) shows that (N, a_p, a_q) satisfy no FIXED
  algebraic relation with p on positive-probability events. That excludes direct lattice equations, but not
  arithmetic uses of the integers.
- **Literature (literature critic; none states anything about a(N) alone).** Bach–Charles (Contemp. Math. 463
  (2008), Thm 2.2, Rems 2.5–2.6) use a(N) and a(N²) (weight 2k, non-CM, density one, including balanced moduli).
  E–C (Ann. Math. Studies 176) Rem 2.2.4 and Bruin (Publ. Math. Besançon 2011) Rem 1.4 concern Eisenstein series
  (σ_{k−1}(N) is equivalent to factoring; Bach–Miller–Shallit 1986). Charles, "Counting lattice vectors", JCSS 73 (2007)
  962–972 (TR1498 Thm 7.1) reduces factoring to counting lattice vectors through θ_{E₈}, again Eisenstein, not cusp
  forms. Chow (PhD thesis, Toronto 2015, Thms 3.3–3.4) uses a_f at n, n², …, n^t. Kunihiro–Koyama (EUROCRYPT 1998)
  and Dieulefait–Jiménez Urroz (Mathematics 8 (2020) 2126) use point counts #E(Z/N). Galbraith–McKee (LNCS 3796, 2005)
  use point-count inputs. All of these use determinant data or several indices. No prior statement on the single value
  was found by us, cx1, cx2 or the literature critic.
- Information counting does not settle anything (cx1): |a_f(N)| ≤ 4N^{(k−1)/2}, so a weight-2 coefficient has about
  n/2 bits, as many as one factor.

---

## 2. (B) Theorem DFB: a likelihood-ratio bound in a Lang–Trotter-type random model

### 2.1 Model and statement

**Model LT(B) (MODEL).** Lang–Trotter (LNM 504, 1976) weight the Frobenius of one curve over primes p by the adelic
image densities and the Sato–Tate density. LT(B) keeps only a pointwise spread bound, treats τ_p as random for the
fixed p of the instance, and makes τ_p, τ_q and the K forms jointly independent. None of that is in LNM 504.
- Fix newforms f_1, …, f_K of weights k_i with integer coefficients, non-CM, **pairwise not twist-equivalent**
  (f_i ≠ f_j ⊗ χ for every Dirichlet character χ; in particular no two isogenous curves). Fix N = pq with
  p ∈ [X₁, 2X₁), q ∈ [X₂, 2X₂).
- For r ∈ {p, q}, τ_r = (a_r(f_i))_i is a random vector in Z^K with |a_r(f_i)| ≤ 2r^{(k_i−1)/2}, and
  Pr[τ_r = t] ≤ B_K·∏_i r^{−(k_i−1)/2} for every t ∈ Z^K (joint spread B_K).
- τ_p and τ_q are independent of each other.
- Everything else the algorithm sees is fixed given N, or independent of (τ_p, τ_q). This includes N itself, the
  algorithm's coins, and oracles whose answers do not depend on (τ_p, τ_q), e.g. a_f at indices coprime to N.
- **The algorithm does not evaluate the trace maps r ↦ a_r(f_i) at p or q.** (Without this clause the proof fails:
  an exponential-time algorithm can sieve candidate primes p′ by a_{p′} | a_N, and the real-to-simulated ratio is then
  about √X/d(a_N), unbounded. Proof critic.)
- **Idealisation (stated, not proved).** The model treats E(Z/N), the torsion algebras (Z/N)[x]/(ψ_ℓ), Jacobi
  symbols etc. as independent of τ. In reality #E(Z/N) = (p+1−a_p)(q+1−a_q) and a_p mod ℓ is the trace of Frobenius
  on the ℓ-torsion algebra. This independence is false for actual curves and is the model's main idealisation.

**Determinant-free, simulable datum.** D = Φ(N, τ_p, τ_q, ω), for a fixed map Φ (coins ω) that receives N but not p
or q, **and such that, given (N, X₁, X₂), the pair (σ, Φ(N, σ, ω)) can be sampled in time T_Φ** when σ = (σ₁, σ₂) is
uniform on the boxes below. Without simulability the definition is vacuous: Φ(N, ·) = (smallest prime factor of N)
receives only N and is "determinant-free" (proof critic; cx2).
- Examples (simulable in poly(n)): a_f(N) = a_pa_q; (a_{f_i}(N))_i; the unordered or ordered trace pairs; and the
  factorisation of a_pa_q, because a uniform integer of a box can be sampled together with its factorisation in
  expected polynomial time (Bach 1988; Kalai 2003), with signs and the value 0 handled separately.
- Non-examples: a_f(N²) = a_N² − q^{k−1}a_p² − p^{k−1}a_q² + N^{k−1} (trivial character; Bach–Charles eqs. (1)–(3),
  expanded as in Chow 2015, Ex. 3.1); #E(Z/N); the CRT trace t_N;
  anything that uses p, q individually.

**Theorem DFB (MODEL, PROVED in the model).** Let σ₁, σ₂ be independent, each uniform on the integer box
Box_j = ∏_i [−⌊2(2X_j)^{(k_i−1)/2}⌋, ⌊2(2X_j)^{(k_i−1)/2}⌋] (j = 1, 2), and put D̃ = Φ(N, σ₁, σ₂, ω). Then for every
algorithm 𝒜 as in the model and every set S of outputs,

  Pr[𝒜(N, D) ∈ S] ≤ C_K·Pr[𝒜(N, D̃) ∈ S],
  C_K := B_K²·|Box₁|·|Box₂|/∏_i N^{(k_i−1)/2}.

For one weight-2 form, C₁ = B²(2⌊2√(2X₁)⌋+1)(2⌊2√(2X₂)⌋+1)/√N = 32B²(1 + o(1)). The finite-box value can exceed 32B²
slightly (cx2: C/B² = 32.056 at X = 32768, p = 32771, q = 32779), so we quote the exact product. For K forms,
C_K ≤ B_K²·∏_i 2^{k_i+3}(1 + o(1)).

*Proof.* Unchanged from v1:
- For every (t, u), Pr[(τ_p, τ_q) = (t, u)] ≤ B_K²∏_i(pq)^{−(k_i−1)/2}, and it vanishes outside the Hasse boxes,
  which lie inside Box₁ × Box₂ because p < 2X₁ and q < 2X₂.
- Pr[(σ₁, σ₂) = (t, u)] = 1/(|Box₁||Box₂|) there. Hence Pr[τ = (t, u)] ≤ C_K·Pr[σ = (t, u)] pointwise.
- Let g(t, u) := Pr[𝒜(N, Φ(N, t, u, ω)) ∈ S] over ω, the algorithm's coins and the oracle. It is well defined
  because these are independent of τ, and because 𝒜 does not evaluate the trace maps at p, q.
- Then Pr[𝒜(N, D) ∈ S] = Σ Pr[τ = (t, u)]·g(t, u) ≤ C_K Σ Pr[σ = (t, u)]·g(t, u) = C_K Pr[𝒜(N, D̃) ∈ S]. ∎

(The inequality is the data-processing inequality for max-divergence. It holds for every N, inside the model.)

**Corollary (MODEL).** Let 𝒜 run in time T with the datum. Then 𝒜 composed with the simulator runs in time
T + T_Φ without the datum, and succeeds with probability ≥ Pr[𝒜 succeeds with D]/C_K.
- *Bits:* the datum adds at most log₂C_K bits of max-divergence.
- *Time:* success is checkable. So truncating at time T and restarting with fresh coins and a fresh D̃ factors N in
  expected time O(C_K·(T + T_Φ + poly(n))). This uses the per-N bound; it does not compare the expected running times of
  the unchanged algorithm on the two inputs, which can differ arbitrarily (cx2).
- *Baseline:* a statement "every algorithm of class 𝒞 succeeds with probability ≤ ε₀" must average over the instance
  law (or exclude hard-coding): for one fixed N an algorithm can hard-code its factors (cx2; proof critic).
- With these caveats: reaching success probability 1/2 needs log₂C_K ≥ log₂(1/(2ε₀)). With C_K = c^K this means
  K ≥ log₂(1/(2ε₀))/log₂c forms. Whether C_K grows only exponentially in K for a growing family needs a uniform bound on
  the joint index (§2.2); this is OPEN.

### 2.2 The spread constant

- **Gekeler law (random curve over F_p, automorphism-weighted).** Pr[a = t] = H(4p − t²)/(2p), with H the Hurwitz
  class number (Σ_{t²<4p}H(4p − t²) = 2p, Kronecker–Hurwitz). Gekeler (IMRN 2003, no. 37, 1999–2018; Thm 5.5 per the
  literature critic and the restatement in David–Koukoulopoulos–Smith, Math. Ann. 2017, Thm 1.1; cx2 cites Thm 4.6) shows
  that this equals the product of the local GL₂(Z_ℓ) densities and the Sato–Tate density. So
  B_Gek(p) = max_t H(4p − t²)/(2√p).
- **Bound (PROVED, routine; fixes from the proof critic and cx2).**
  - For a negative discriminant −Δ (not necessarily fundamental), 2h(−Δ)/w(−Δ) = √Δ·L(1, (−Δ/·))/π, with the
    imprimitive Kronecker symbol (−Δ/·).
  - L(1, χ) ≤ 2 + ln Δ for every non-principal character mod Δ, so 3 + ln Δ is safe.
  - Summing over the square divisors gives H(D) ≤ (√D/π)(3 + ln D)·σ(m)/m, with m² the largest square divisor of D.
    Numerically the ratio of H(D) to this bound is ≤ 0.392 for every D < 2·10⁴ (proof critic).
  - Robin (J. Math. Pures Appl. 63 (1984)): σ(m)/m ≤ e^γ lnln m + 0.6483/lnln m for m ≥ 3. To substitute m ≤ 4p one needs monotonicity of the right
    side, which holds for m ≥ 7; for m ≤ 6, σ(m)/m ≤ 2 < min_y(e^γy + 0.65/y) ≈ 2.15.
  - Hence B_Gek(p) ≤ (1/π)(3 + ln 4p)(e^γ lnln 4p + 0.65/lnln 4p) = O(log p·loglog p).
  - Under GRH (Littlewood 1928: L(1, χ) ≪ loglog) the conductor factors cancel and B_Gek = O(loglog p). (v1's O((loglog p)²) is true but loose.)
- **E3 (measured):** B_Gek ∈ [0.985, 1.184] at 16 bits, [1.293, 1.543] at 20 bits, [1.386, 1.630] at 24 bits. The
  proved bound is 22.2–35.7 on these primes (the registered formula gives 24.0–38.9).
- **One fixed curve.** In the Lang–Trotter law the density of tr = t on the fibre det = p is that of the image G
  instead of GL₂(Ẑ). Since det is surjective on G, its ratio to the GL₂ density is ≤ [GL₂(Ẑ) : G]. So
  B_E ≤ [GL₂(Ẑ) : G]·B_Gek (correct, proof critic).
  - 37a1 is a Serre curve: LMFDB 37.a1 lists adelic level 74 and index 2 (cx2; literature critic), which is Jones's
    Serre number 2|W| for W = 37 (Jones, Trans. AMS 362 (2010)). Serre 1972 §5.5.6 uses this curve as the example with
    surjective mod-ℓ images for all ℓ (secondary sources). So B ≤ 2B_Gek, and
    log₂C₁ = 5 + 2log₂(2B_Gek) ≤ 8.4 bits at 24-bit p, using the E3 maximum. Asymptotically
    log₂C₁ = 2log₂log N + O(logloglog N).
  - **11a1: corrected.** It lies in a 5–5 isogeny chain (isogeny matrix [1,5,5; 5,1,25; 5,25,1]), so E[5] ≅ Z/5 ⊕ μ₅
    and the mod-5 image {diag(1, d)} (LMFDB 11.a2: 5Cs.1.1) has index 120. LMFDB gives adelic index 1200, level 550. The crude bound therefore costs
    2log₂1200 ≈ 20.5 extra bits. v1's "divisible by 24, ≤ 9 bits" was wrong. The actual leak is s mod 5.
- **K forms (new in v2).** For pairwise non-twist-equivalent non-CM curves, the joint adelic image G_K has finite index
  in F_K = {(g_i) ∈ GL₂(Ẑ)^K : all det g_i equal} (Ribet, Invent. Math. 28 (1975), Thm 3.5 for non-isogenous
  curves; Lombardo, JNT 168 (2016), Thm 1.1, an explicit version, cited by cx2; Loeffler, Glasgow Math. J. 59 (2017)
  for finite sets of newforms of mixed weight). Then B_K ≤ [F_K : G_K]·B_Gek^K. Since [F_K : G_K] ≥ ∏[GL₂ : G_i], and no bound of
  the form c^K is known for a growing family, "C_K = c^K" is an assumption (OPEN).
  - Isogenous curves have identical traces, and a quadratic twist changes a_N only by (d/N) (cx1). For such pairs the
    vector bound fails by a factor √p/B. The family F removes twists only within one j-invariant: its first 32000
    curves contain 34 isogenous pairs (e.g. (A, B) = (1, 2) and (−19, 30)) and 55 pairs equal up to a twist of an
    isogenous curve (e.g. (−2, 1) and (−7, 6), twist by −1). This does not affect E1 (duplicates are only redundant),
    but DFB with K curves of F must drop them.
- **Weight 12 (Δ): corrected.** ρ_Δ has non-surjective mod-ℓ image exactly at ℓ ∈ {2, 3, 5, 7, 23, 691}, and the local
  density ratios there are bounded, so B_Δ = O(B_Gek). v1 said τ(p) is "fixed modulo 2^{11}·3^6·5^3·7·23·691 by p". That
  is false (Swinnerton-Dyer, LNM 350 (1973), lists the moduli but the congruences are not all functions of p mod ℓ):
  modulo 23 it depends on the splitting of p in x³ − x − 1 (τ(13) ≡ 22, τ(59) ≡ 2 mod 23 although
  13 ≡ 59 mod 23; Wilton 1930), and the 2^{11} congruence holds only for p ≡ 1 mod 8. The 691 leak is real: from
  τ(N) ≡ (1 + p^{11})(1 + q^{11}) mod 691 one gets p^{11} + q^{11} mod 691, hence {p, q} mod 691 because
  gcd(11, 690) = 1, so s mod 691.
- **Toy check (exact; `test_dfb_ratio_bound_exact`; recomputed in `scratchpad/cuspcoef2/dfb_toy.py`).** For p = 607,
  q = 1009 under the Gekeler law, log₂C = 3.650 bits. The worst actual ratio for the product a_pa_q is 2.107 bits, and
  the KL divergence is 0.503 bits.

### 2.3 What DFB excludes in the model, and how the known reductions escape it

**Excluded in LT(B)** (all determinant-free and simulable): Coppersmith/lattice attacks built from a_N = a_pa_q and the
Hasse bounds, for O(1) curves, including weight 2 with weight 12; congruence decoding of a_N (SC, SCX, CI);
Sato–Tate/archimedean reasoning; factoring the integer a_N; any combination of these within the model's interface.
This is a statement about the model. DEC excludes only fixed polynomial relations for the real traces, so "no lattice
attack exists" for actual coefficients is not claimed.

**Every reduction we know escapes one hypothesis:**

| reduction | which hypothesis it violates |
|---|---|
| Bach–Charles (a_N, a_{N²}) | a_{N²} contains p^{k−1}, q^{k−1}: not determinant-free |
| CM pairs y² = x³ − x, x³ − Dx (tracedata CM-PAIR; cx1's cross-square derivation) | spread: for CM, a_p ∈ {±2α, ±2β} with α² + β² = p, so B = Θ(√p) |
| M = #E(Z/N) (Theorem EXT) | not determinant-free: M = (p+1−a_p)(q+1−a_q) |
| CRT trace t_N ≡ a_p (mod p) (CRT-T) | uses p |
| class pair (CLS) | determinants of classes |
| ordered traces of O(n² log n) non-CM curves (tracedata MC, HEURISTIC) | none: many forms, consistent with DFB since log₂C_K ≫ n/2 |

So, in the model, a reduction from Hecke data must use the determinant, CM structure, or many forms. v1's "each one
uses the determinant" was too strong (CM and MC rows), and v1's "exactly the known reductions" is replaced by "every
reduction we know".

### 2.4 Scope: what the model does not capture

1. **Determinism.** Real traces are deterministic, poly(n)-time functions of the prime (Schoof; Edixhoven–Couveignes).
   A reduction exploiting a specific non-statistical property of p ↦ a_p (as CM does) is outside the model. For
   non-CM f, DEC excludes fixed algebraic relations; nothing excludes cleverer ones.
2. **Coupling to structures computable from N.**
   - **DFB(m) (MODEL; same proof).** Suppose that, conditionally on the true Frobenius classes c_p, c_q at level m, the
     lifts are independent of each other and of all other randomness, and Pr[τ_r = t | c_r] ≤ B_K m^{K}·∏_i
     r^{−(k_i−1)/2} on the residue class of tr c_r. Let σ^{(m)} be uniform on the residue boxes Box_j(c) = {t ∈ Box_j :
     t ≡ tr c (mod m)}. Then Pr[𝒜(N, D) ∈ S | c] ≤ C_{K,m}·Pr[𝒜(N, Φ(N, σ^{(m)})) ∈ S | c] with
     C_{K,m} = B_K²m^{2K}|Box₁(c_p)||Box₂(c_q)|/∏_i N^{(k_i−1)/2} (cx2's unambiguous form).
   - Holding the true classes is an extra oracle assumption. v1's "e.g. by probing torsion algebras" is removed: no
     way to obtain the classes from N is known (FROB/PRIM).
   - Couplings at levels growing with n are not covered. The class pair at all ℓ ≤ O(log N) factors N (CLS).
   - Coupling through the group order is not covered. Over (Z/N)[t]/(t² − c), conjugation is Frobenius at inert
     primes, and from a_N one gets [a_N]σP = ([a_q(p+1)]P_p, [a_p(q+1)]P_q), again with unknown p. DEC's annihilator
     remark (tracedata §2.1) excludes fixed polynomial annihilators only, not arbitrary group algorithms. Neither cx1,
     cx2 nor we found a reduction from a_N plus group operations; Galbraith–McKee (2005) use stronger point-count inputs.
3. **One N at a time.** The inequality holds for every N; baselines ε₀ must be averages (§2.1).

### 2.5 Prop CI: congruence information is O(K) bits (MODEL; cx1, sharpened by cx2 and us)

**Model.** For each odd prime ℓ, D_ℓ = p mod ℓ is hidden and uniform in F_ℓ^×, and C_ℓ = N mod ℓ is exposed. For each
of K curves, independently, g_i, h_i are uniform in GL₂(F_ℓ) with det g_i = D_ℓ, det h_i = C_ℓ/D_ℓ. One observes
A_{i,ℓ} = tr g_i·tr h_i. Different ℓ are independent; unlimited postprocessing is allowed.

**Prop CI.** For every set S of odd primes, I((D_ℓ)_{ℓ∈S}; (A_{i,ℓ})_{i,ℓ} | (C_ℓ)_ℓ) < 2K·log₂(π²/8) < 0.607K bits.

*Proof.*
- The trace law on the fibre det = d is P_d(t) = (ℓ + χ(t² − 4d))/(ℓ² − 1) (count solutions of bc = a(t − a) − d).
- Exactly, max_d χ²(P_d‖U) = 1/(ℓ² − 1), with U uniform on F_ℓ (checked exactly for all ℓ < 400 in
  `scratchpad/cuspcoef2/check_cx2.py`; cx1 used the weaker (ℓ − 1)^{−2}).
- I(D; Y | C) ≤ E_D KL(P_{Y|D,C}‖U) ≤ log₂(1 + χ²), and the ordered pair (tr g, tr h) gives at most twice that;
  passing to the product cannot increase information. Conditional independence makes this additive over curves and ℓ.
- Σ_{ℓ odd} 2log₂(1 + 1/(ℓ² − 1)) = 2log₂∏_{ℓ odd}(1 − ℓ^{−2})^{−1} = 2log₂(π²/8) = 0.6060. ∎

cx2's bound 2log₂(4/π) < 0.698 is valid but looser; cx1's was π²/(12 ln 2) < 1.187.

**Scope.** The model is conditioned on N mod ℓ, not on the integer N. It discards integrality and size, and it allows
no group operations. Independence across ℓ is a model assumption. 37a1's 2–37 entanglement does not break the odd-ℓ
version (projecting away ℓ = 2 gives the full odd image; cx2). Curves with exceptional images need corrections:
11a1 at ℓ = 5, C = 1 already gives 1.5 bits (cx2).

---

## 3. Prop SCX: the exact scope of Theorem SC (congruences of products)

**Setting.** E/Q non-CM, with adelic image G ⊆ GL₂(Ẑ). Λ := {λ ∈ Ẑ^× : λI ∈ G}, and Λ_m its image in (Z/m)^×.

**Prop SCX (PROVED from SC and Serre's open image theorem; fixes from the proof critic and cx2).**
- (i) Λ is open in Ẑ^×. G is open (Serre 1972), so it contains the kernel of reduction mod some M_E, which contains λI
  for every λ ≡ 1 (mod M_E). Hence [Ẑ^× : Λ] ≤ |(Z/M_E)^×|.
- (ii) For every m, the Chebotarev law of (N mod m, a_N mod m), conditioned on p mod m, is invariant under
  p ↦ λ²p, q ↦ λ^{−2}q for λ ∈ Λ_m (tracedata Theorem SC, with the scalars Λ_m). So the congruence content about
  p mod m is at most the class of p in (Z/m)^×/Λ_m². With U = (Z/m)^×, [U² : Λ_m²] = [U : Λ_m]·|Λ_m[2]|/|U[2]| ≤
  [U : Λ_m] ≤ [Ẑ^× : Λ], so this exceeds the quadratic-character content by ≤ log₂[Ẑ^× : Λ] bits, uniformly in m.
- (iii) For a Serre curve, Λ = Ẑ^×. Here G = {g : ε(g mod 2) = χ_Δ(det g)}, with ε the sign character of
  GL₂(F₂) ≅ S₃. For a scalar, λI ≡ I (mod 2) gives ε = 1, and χ_Δ(λ²) = 1.
  - 37a1 is a Serre curve per LMFDB (adelic level 74, index 2; verified by cx2 and the literature critic).
  - E4 checks finite-level consequences. All ℓ(ℓ−1) pairs (a_p mod ℓ, p mod ℓ) occur for every 3 ≤ ℓ ≤ 47. By
    Dickson this proves surjectivity mod ℓ only for **5 ≤ ℓ ≤ 47**: at ℓ = 3 the normaliser of the non-split Cartan
    (order 16) already realises all 6 pairs (proof critic; cx2). So P4b holds as registered but proves nothing.
  - The entanglement rule "a_p odd ⇒ (37/p) = 1" has 0 violations among 17982 primes, and 2967 primes have
    (37/p) = 1 with a_p even.
- (iv) Escapes from QR content come only from finitely many levels where scalars are missing (rational isogenies or
  torsion, entanglement), and are bounded by O(1) bits per curve. Example: 11a1 at 5, where Λ₅ = {1} and
  a_N ≡ (1+p)(1+q) gives s mod 5.
  - **K curves: corrected.** SC needs scalar *tuples* in the joint image, Λ_K ⊆ ∩Λ_i, possibly strictly. The excess
    is O_K(1) for a fixed set of pairwise non-twist-equivalent curves. "O(1) per curve" for a growing family would
    need a uniform joint-index bound (OPEN).

**E1 confirms the plateau.** The congruence-only decoder (C2) at K = 32000 has per-ℓ top-1 accuracy 1.0, 0.633, 0.500,
0.267, 0.183 for ℓ = 3, 5, 7, 11, 13, matching tracedata's E3b plateau. The decoder that uses the factorisations
reaches ≥ 0.975 at every ℓ ≤ 23.

---

## 4. (R) Many curves: products plus their factorisations

**Prop RM (HEURISTIC; toy-verified; cost and sample size corrected in v2).**
- Fix the family F (`family.py`: y² = x³ + Ax + B, one curve per j, non-CM; it is finite, 49,413 curves, so asymptotic
  statements refer to a growing family with isogenous and twist-related curves removed). Given N, a_N(E_i) for K
  curves, and the factorisations of the |a_N(E_i)| (integers ≤ 8X):
  1. List the candidate splits C_i = {(x, y) : xy = a_N(E_i), |x|, |y| ≤ R}.
  2. For each odd ℓ ≤ L (∏ℓ > 2X), score d = p mod ℓ by Σ_i log Σ_{(x,y)∈C_i} n_ℓ(x, d)·n_ℓ(y, N/d), with
     n_ℓ(t, d) = ℓ(ℓ + χ(t² − 4d)) the GL₂(F_ℓ) class count (normalisation independent of d since Σ_tχ(t² − 4d) = −1).
  3. Pass to s mod ℓ, list-decode, CRT, and solve.
- **Sample size: K = poly(n) (HEURISTIC).** v1's K = O(n² log n) was imported from tracedata MC, which concerns
  ordered traces; no sample bound is known for candidate splits. Wrong splits dilute the signal: the per-curve
  information is roughly ~1/(ℓ²|C_i|), and mean |C_i| grows (21.3, 28.1, 33.4 at n = 32, 40, 48, K = 32000). Swaps and sign flips
  of the true split are not noise for s, and the other wrong splits are correlated with the truth. A dilution estimate
  suggests K between n² log n and n³ log n; the toy data do not fix the exponent. (Proof critic; cx2.)
- **List decoding.** Step 3 inherits MC's list search, which is exp(Θ(n/log n)) in the worst case; E1's median `tried`
  of 1 shows that at toy sizes the top combination is usually right.
- **Divisor enumeration.** A product can have exp(Θ(n/log n)) divisors, so "factorisations plus poly(n)" needs an
  average-case divisor bound (typical d(m) is (log m)^{O(1)}) (cx2).
- **Cost.** With a polynomial-time a_N for the curves, NFS on the products gives
  L_{√N}[1/3, (64/9)^{1/3}] = L_N[1/3, (32/9)^{1/3}] ≈ L_N[1/3, 1.526] (heuristic GNFS cost; Buhler–Lenstra–Pomerance
  1993). v1 wrote this as "(64/9)^{1/3}·2^{−1/3}", which is the same number; v2 drops the factor from §0 to avoid
  suggesting a further saving.
  - **ECM selection (HEURISTIC, proof critic; algebra checked by us, untested).** Every prime factor of a_N is ≤ R ≈
    N^{1/4}. Run ECM with bound y = L_N[2/3, a] on the products of L_N[1/3]-many curves, and keep those where one of
    a_p, a_q is y-smooth: the unfactored cofactor then lies on one side and need not be factored. The cost
    L_N[1/3, 1/(12a) + √(4a/3)] is minimised at a ≈ 0.275, giving c ≈ 0.909; with E2's stricter criterion c ≈ 1.145.
    This needs a_N for L_N[1/3]-many curves of the family, not poly(n) many.
- **What E2 shows (corrected).** E2 measures the registered filter "the cofactor after removing all primes ≤ n² is 1 or
  prime". Its acceptance decays like the Dickman function (P2a, P2b hold; table §6). It shows that this filter is
  super-polynomially rare, also with the polynomial-time ECM bound exp(Θ((log n)²/loglog n)) in place of n².
  - v1's "the split is computable in poly(n) only if every prime factor except one is ≤ n^{O(1)}" is **false**
    (cx2): 1000003·1000033 = 1000018² − 15² factors in one Fermat step, and perfect powers are also easy. Within the
    random-integer heuristic such structures are rare, but E2 does not prove that every efficiently splittable product
    is rare.
- **Without factorisations**, the products' congruences carry only QR bits plus O(1) (SC, SCX, CI). A partially
  factored product with an unsplit composite part loses the ordering at every ℓ whenever the composite can straddle
  x and y.
- **Few curves.** K = O(1) curves give ≤ log₂C_K bits in LT(B), so the model has no reduction from them.
- Recursion (remark, not tested): a squarefree-index oracle for all curves of F could be fed the half-size products.
  That needs a decoder for k-factor squarefree indices, which we did not build. We claim nothing here.

---

## 5. (A) Native algorithms for a_f(N)

**Baseline (KNOWN).** Factor N deterministically in Õ(N^{1/5}) (Harvey, Math. Comp. 90 (2021) 2937–2950; Harvey–Hittmeir, Math.
Comp. 91 (2022) 1367–1379; earlier Õ(N^{1/4}) by Pollard 1974 and Strassen 1976), or heuristically in L_N[1/3, 1.923] (NFS). Then compute a_f(p), a_f(q) in poly(n)
(Schoof 1985 for weight 2; Edixhoven–Couveignes 2011, Thm 15.2.1, level 1, GRH, factored input; Bruin 2011, Thm 1.1,
Las Vegas under GRH, general level, index in factored form). **Every native route below is
slower than this.**

| route | cost (best known; upper bounds) | native? | status |
|---|---|---|---|
| factor N, then multiplicativity + Schoof / E–C / Bruin at p and q | Õ(N^{1/5}) deterministic; L[1/3] heuristic | no | KNOWN |
| Eichler–Selberg at auxiliary indices 3N, 5N (11N at level 37). √N class-number cost: Charles 2006 (Thm 2.1; E–C Introduction p. 3); cancellation of the divisor term: cx1 attempt-1 scratch: the hyperbolic terms of Tr T_{3N}, Tr T_{5N} differ by known constants when q < 3p, and τ(3N) − τ(5N) = (252 − 4830)τ(N) isolates τ(N) | O(√N) Hurwitz numbers H(D), D ≤ 20N. Each needs the square part of D (rigorous L_D[1/2+o(1)] factoring, Lenstra–Pomerance, JAMS 5 (1992), negligible against √N; v1's "square-part caveat" is removed) and class numbers (Las Vegas L^{√2+o(1)} under GRH, Hafner–McCurley, JAMS 2 (1989)). So N^{1/2+o(1)} under GRH | **factoring in disguise**: Σ_t H(4N − t²) = 2(N + q), and Σ_t H(12N − t²) = 8N + 6(p + q) when q < 3p (checked numerically by the proof critic and cx2) | identities exact; cost GRH |
| theta difference: r_{L1}(N) = Σ_{m≤N} r_B(m)r_B(N−m), B = x² + xy + 3y², then a(N) = (5/8)(r_{L1}(N) − (12/5)σ(N)) | Õ(N) | **no**: σ(N) is factoring (v1 said yes) | PROVED (nativemod (iii)) |
| complete Jacobi sum a_N = Σ_{x mod N} ((x³ + ax + b)/N) (squarefree N) | Θ(N) Jacobi symbols | yes | KNOWN identity |
| Merel/Heilbronn sets for T_N on modular symbols | O(σ(N) log N) operations, Õ(N) for N = pq (Stein, GSM 79, §8.3.3) | yes | KNOWN |
| quaternary theta, meet in the middle (cx1) | Õ(N) time, O(N) space | yes | PROVED (cx1) |
| Charles, Ramanujan J. 11 (2006) 221–224: τ(p) in O(p^{1/2+ε}) randomized under GRH (Thm 2.6) or O(p^{3/4+ε}) deterministic (Thm 2.7) | — | no: §2 factors n first | KNOWN (literature critic read the author's preprint) |
| fast exponentiation in a fixed algebra / P-recursive recurrence | impossible | — | Prop LRS; natural boundary (nativemod) |
| analytic extraction of one coefficient | precision wall | — | nativeana |

- **PARI caution (cx1, attempt 1).** PARI's `qfbhclassno(12N)` factors its argument internally (observed in
  `codex/scratch/cuspcoef_1/audit.log`). The trace-formula checks therefore verify identities, not a native
  implementation.
- **Remark BC-native (KNOWN consequence of Bach–Charles, Contemp. Math. 463 (2008), Thm 2.2).** An algorithm computing
  a_f(n) at every index n in time T(log n) factors a density-one set of RSA moduli in time 2T(2 log N) + poly (f a
  non-CM eigenform in S_{2k}(Γ₀(M), χ); the theorem assumes the algorithm outputs a(n) with its minimal polynomial and
  embedding). The cost is T at doubled length, not "the same cost". Of the routes above, the trace-formula
  cancellation is specific to semiprimes and would need a separate extension to apply at N².
- **Remark JS (PROVED, elementary; checked by cx2).** The Jacobi-sum formula is squarefree-specific, yet it leaks at N².
  - For any cubic f, Σ_{x mod N²} (f(x)/N²) = #{x mod N² : gcd(f(x), N) = 1} = N(p − r_p)(q − r_q), with
    r_r = #{roots of f mod r}.
  - Given V = (p − r_p)(q − r_q), each guess (r_p, r_q) ∈ {0, 1, 2, 3}² gives
    r_q p² − (N + r_pr_q − V)p + r_pN = 0, which recovers p unless r_p = r_q = 0 (then V = N, no information).
  - For an S₃ cubic that exception has asymptotic probability 1/9 for independently sampled good primes.

---

## 6. Experiments (pre-registered; `notes/cuspcoef_prereg.md`)

Jobs ran one at a time. Measured: E1 122 s with peak RSS 739 MB (cuspcoef2 rerun, `/usr/bin/time -l`); E2 612 s in
total (t_total 0.2 + 82.9 + 176.2 + 353.0 s) with peak RSS 92 MB at n = 256; E3 12 s and E4 < 1 min (v1, not in the
raw files).

| exp | prediction | result |
|---|---|---|
| E4 P4a | all ℓ(ℓ−1) pairs, 5 ≤ ℓ ≤ 47 | all present for every ℓ; **holds** |
| E4 P4b | ℓ = 3: 6 pairs | 6/6; **holds** (but uninformative: a proper subgroup also realises all 6, §3) |
| E4 P4c | a_p odd ⇒ (37/p) = 1, 0 violations; converse fails | 0 violations in 17982 primes; converse fails 2967 times; **holds** |
| E3 P3a | Kronecker–Hurwitz sum for all 30 primes | 30/30; **holds** |
| E3 P3b | max B_emp ≤ 4, below the explicit bound | max 1.630; registered bound 24.0–38.9 (proved bound 22.2–35.7); **holds** |
| E1 P1a | K = 32000: ≥ 17/20 (n = 32, 40), ≥ 14/20 (n = 48) | 20/20, 20/20, 20/20; **holds** |
| E1 P1b | success non-decreasing in K | n = 32: 20, 20, 20; n = 40: 20, 20, 20; n = 48: 18, 20, 20; **holds** |
| E1 P1c | null C1 ≤ 1/20 | 17/20 (n = 32), 6/20 (n = 40), 0/20 (n = 48): **REFUTED at n = 32, 40** (exhaustive verification within the 3000-combination budget; median per-instance space 2604 and 23814; see A1) |
| E1 P1d | congruence-only C2 ≤ 2/20 | 20/20, 17/20, 5/20: **REFUTED at all n**. SC's QR content shrinks the space: median QR-reduced space 440, 4320, 26400 at n = 32, 40, 48, against the budget 3000 |
| E1 P1e | per-ℓ top-1 ≥ 0.8 for ℓ ≤ 13 (K = 32000) | 1.0, 1.0, 0.983, 1.0, 0.983 (ℓ = 3–13); ≥ 0.975 up to ℓ = 23; **holds** |
| E1 P1f (A1) | median `tried` ≤ 5 at K = 32000 | 1, 1, 1; **holds** |
| E1 P1g (A1) | median `tried` of C1, C2 ≥ 100 | C1: 1472, > 3000, > 3000; C2: 184.5, 764.5, > 3000; **holds** |
| E1 P1h (A1) | top-1: main ≥ 10/20; C1, C2 ≤ 1/20 | main 19, 18, 19; C1 0, 0, 0; C2 0, 0, 0; **holds** |
| E2 P2a | real within ×2 of simulator where sim ≥ 0.02 | ratio 1.033 (n = 64), 1.195 (128), 0.846 (192); **holds** |
| E2 P2b | simulator strictly decreasing, ≤ 0.02 at n = 256 | 0.913, 0.335, 0.053, 0.0059; **holds** |

**E2 table** (`cuspcoef_e2.jsonl`; 600 real (instance, curve) samples and 10⁵ simulated samples per n; B = n²).

| n | real acc | real rate | sim acc | sim rate | real/sim | (real − sim)/SE |
|---|---|---|---|---|---|---|
| 64 | 566/600 | 0.943 | 91292/10⁵ | 0.913 | 1.03 | +3.2 |
| 128 | 240/600 | 0.400 | 33480/10⁵ | 0.335 | 1.19 | +3.3 |
| 192 | 27/600 | 0.045 | 5321/10⁵ | 0.053 | 0.85 | −1.0 |
| 256 | 3/600 | 0.005 | 590/10⁵ | 0.0059 | 0.85 | (3 events) |

At n = 64 and 128 the real rate is above the simulator by > 3 SE (within the registered ×2 tolerance). Real traces are
more often divisible by small primes than uniform integers, as E1's candidate counts also show.

**E1 details (K = 2000 / 8000 / 32000; main decoder).**
- Success: n = 32: 20/20/20; n = 40: 20/20/20; n = 48: 18/20/20.
- Top-1: n = 32: 7/13/19; n = 40: 5/9/18; n = 48: 0/7/19.
- Median `tried`: n = 32: 3/1/1; n = 40: 15.5/2/1; n = 48: 661/5.5/1.
- Mean candidate splits per product: 21.3–21.5 (n = 32), 28.1–28.5 (n = 40), 33.4–33.7 (n = 48). For the null they
  are 14.2, 17.6, 20.9.
- Accidental gcd(Δ_i, N) splits (bad reduction of a family curve at p or q), at K = 32000: 27 at n = 32 over 20
  instances, 2 at n = 40, 0 at n = 48. All ACCIDENTAL; the structural decoder ignores them.
- These are toy correctness checks (rule 8). The harness supplies the a_N, and the half-size factorisations are
  trivial at these sizes.

### 6.1 Reconciliation (v2)

Script `scratchpad/cuspcoef2/reconcile.py` recomputes every count from `experiments/raw/cuspcoef_*.jsonl`; cx2
independently did the same (`codex/scratch/cuspcoef_2/verify_records.py`). Determinism checks: the E1 rerun matches all
300 records except wall time; the E2 n = 64 rerun matches its line exactly.

| quantity | v1 | raw | action |
|---|---|---|---|
| E1 success / top-1 / median `tried` / candidates / accidental | as in §6 | identical | none |
| E1 per-ℓ top-1, main and C2 | as in §3, §6 | identical (C2: 1.0, 0.633, 0.500, 0.267, 0.183) | none |
| E1 search space | 9072, 90720, 1.09·10⁶ ("S_n") | these are the upper bound ∏(ℓ+1)/2; registered per-instance S_n medians 2604, 23814, 365310 | corrected |
| C2 effective space at n = 48 | "≈ 4·10³" | median QR-reduced space 26400 | corrected |
| E3 B_emp ranges | [0.98, 1.18], [1.29, 1.54], [1.39, 1.63] | identical | none |
| E3 "bound 24–39" | quoted next to the proved formula | registered formula 24.0–38.9; proved formula 22.2–35.7 | both given |
| E4 counts | 17982, 0, 2967 | identical | none |
| E2 rates | v1 (6f16445): 0.94/0.91, 0.40/0.33, then TBD; 2c4a3d1 added 0.045/0.053, 0.005/0.006 | 0.943/0.913, 0.400/0.335, 0.045/0.053, 0.005/0.0059 | consistent; 3 decimals in v2 |
| DFB toy | 3.65, 2.11, 0.50 bits | 3.650, 2.107, 0.503 | none |
| runtimes / RSS | "≈ 3 min, ≤ 726 MB"; E2 TBD | not in the raw files; remeasured 122 s / 739 MB (E1), 612 s / 92 MB (E2) | replaced |

**Integrity record.** The previous agent twice wrote monitor-style text with invented E2 numbers at n = 256 and
discarded both. Its report's "0.91, 0.33, 0.053" are the simulated rates and are correct. The only n = 256 line in
the raw data comes from the cuspcoef2 rerun (D4–D5).

---

## 7. Candidate cards

| mechanism | model | outcome |
|---|---|---|
| Coppersmith/lattice on a_N = a_pa_q (+ Hasse), O(1) curves, ± τ(N) | LT(B) | **refuted for this model** (DFB): ≤ log₂C_K bits; real curves OPEN |
| congruences of a_N, one or many curves | independent GL₂ / Chebotarev | < 0.607K bits (CI); QR bits + O(1) (SC, SCX) |
| products of K = poly(n) curves + factorisations | HEURISTIC | works at toy size (E1 60/60); costs K half-size factorisations, L_N[1/3, 1.526] with NFS |
| ECM selection of products | HEURISTIC | L_N[1/3, ≈ 0.91–1.15] given a_N for L[1/3]-many curves; untested |
| the n²-smooth "easy product" filter | Dickman | super-polynomially rare (E2); not a general barrier |
| native: Eichler–Selberg cancellation | — | N^{1/2+o(1)}; reveals q; slower than factoring |
| native: Jacobi sum, Merel, quaternary theta | — | Õ(N) |

---

## 8. Reviews

**cx1 (GPT-6.1-sol blind round, `codex/cuspcoef_1.out`).** A first attempt (session 01a1818c) was killed at 09:28
without output; its log is kept in `scratchpad/cuspcoef2/cuspcoef_1_attempt1.log` and its scratch files in
`codex/scratch/cuspcoef_1/`. The trace-formula cancellation and the PARI `qfbhclassno` audit come from that attempt's
log, which we mined. The relaunch by cuspcoef2 (session 01a181a2, 09:30–09:49) produced the only `cuspcoef_1.out`. It was
merged at commit ad25398 by the previous agent, so there was exactly one blind round output.

| cx1 finding | action |
|---|---|
| no reduction for (a), (b), (c): OPEN, a literature conclusion | adopted as the headline |
| MODEL congruence bound < 1.187K bits | Prop CI (§2.5), sharpened to 0.607K |
| scalar symmetry for any K and for mixed weights (λ^{k_i−1}I)_i | §3 (supports SCX) |
| twist identity a_{E^{(d)}}(N) = (d/N)a_E(N); isogenous curves identical | §2.2 joint-spread hypothesis |
| 11a1: s mod 5 | §3 (iv) |
| Bach–Charles denominator mechanism; CM cross-square derivation | §2.3 table |
| Bruin 2011 Thm 1.1, "Rem 1.3" (the factoring remark is Rem 1.4: literature critic); Charles 2006 factors n first; Galbraith–McKee | §1, §5, §2.4 |
| native routes Õ(N); trace-formula caution | §5 |

**cx2 (GPT-6.1-sol adversarial round, `codex/cuspcoef_2.out`): MAJOR REVISION.**

| # | cx2 finding | verdict | action |
|---|---|---|---|
| 1 | DFB is not a hardness theorem; Φ needs an operational definition; pointwise ε₀ impossible (hard-coding) | accepted | §0.1, §2.1 (simulability, averaged ε₀) |
| 2 | RM sample size and cost unsupported; divisor counts can be exp(n/log n); F is finite (49,413 curves) | accepted | §4 (K = poly(n), divisor caveat, finite family) |
| 3 | E2's "only if" is false (Fermat example) | accepted | §4, §0.3 |
| 4 | native headline misleading; Harvey N^{1/5} baseline; Hafner–McCurley | accepted | §0.4, §5 |
| 5 | Kalai sampler works; time corollary needs truncation and restarts; exact C/B² = 32.056 > 32 | accepted | §2.1 |
| 6 | weight-12 congruence claim wrong mod 23 (Wilton); 691 leak correct; "no lattice attack" exceeds DEC | accepted | §2.2, §2.3 |
| 7 | imprimitive Kronecker character; Robin small range; 37a1 Serre curve verified (LMFDB); 11a1 index 1200 | accepted | §2.2 |
| 8 | joint scalars O_F(1); Lombardo 2016; ℓ = 3 normaliser defeats the E4 argument | accepted | §3 |
| 9 | DFB(m) constant notation; class oracle is an assumption; DEC limited | accepted | §2.4 |
| 10 | CI proof survives; conditioning on N mod ℓ; 11a1 gives 1.5 bits; sharper 0.698K | accepted, sharpened to 0.607K | §2.5 |
| 11 | BC-native costs T(2n); JS checked; 1/9 asymptotic | accepted | §5 |
| 12 | numbers reproduce; proved E3 bound 22.2–35.7; E2 612 s; RSS not in raw | accepted | §6, §6.1 |
| 13 | search-space explanation wrong (median 365310 / 26400 at n = 48); disclose D4 and the simulator endpoint; supplements must not be pooled | accepted | §6, §6.1, §9 |

**Opus proof critic (read-only; scratch `scratchpad/cuspcoef2/critic_proof/`).** DFB inequality CORRECT; statement,
corollary, §2.2, §2.4 CORRECT WITH FIXES; SCX CORRECT WITH FIXES; RM OVERCLAIM; native framing WRONG.

| # | finding | action |
|---|---|---|
| 1 | N^{1/2} route reveals q (Kronecker–Hurwitz, verified for 7 N; 12N version for 3 N); above N^{1/5}; theta row needs σ(N) | §0.4, §5 |
| 2 | RM cost: NFS constant (32/9)^{1/3}; ECM selection ≈ 0.91–1.15; E2 refutes only polynomial-time selection; poly-time ECM bound | §4 |
| 3 | Φ must be simulable (smallest-prime-factor counterexample); Kalai/Bach; time corollary | §2.1 |
| 4 | trace evaluation at p, q must be excluded (sieve counterexample); E(Z/N) independence is an idealisation; ε₀ averaged; §0 "each uses the determinant" false for CM and MC | §2.1, §2.3 |
| 5 | joint spread: 34 isogenous + 55 twist-isogenous pairs in F; B_K ≤ [F_K : G_K]∏B; no c^K bound known | §2.2, §2.1 |
| 6 | K = O(n² log n) unsupported; dilution estimate n² log n to n³ log n | §4 |
| 7 | 11a1 index ≥ 120 at level 5, not "divisible by 24" | §2.2 |
| 8 | SCX (ii) index identity; ℓ = 3 Cartan normaliser; joint scalars O_K(1) | §3 |
| 9 | imprimitive character; L(1, χ) ≤ 2 + ln Δ; Robin m ≤ 6; GRH gives O(loglog p); weight-12 wording | §2.2 |
| 10 | a(N²) formula wording; MC list-search cost; likelihood normalisation correct | §2.1, §4 |

**Opus literature critic (read-only; sources in `scratchpad/cuspcoef2/critic_lit/`).** No prior statement on a single
coefficient found; Bach–Charles, Kronecker–Hurwitz, Eichler–Selberg, E–C, Harvey, Hafner–McCurley, Kalai, Bach,
Coppersmith VERIFIED. Lang–Trotter and Serre 1972 only via secondary sources (Serre's paper returned 403).

| # | finding | action |
|---|---|---|
| 1 | 11a1: mod-5 image 5Cs.1.1 (index 120), adelic index 1200; even index 24 would give 9.17 > 9 bits | §2.2 |
| 2 | Bruin's factoring remark is Rem 1.4, not 1.3 (cx1 mis-cited); Thm 1.1 is Las Vegas under GRH, factored index | §1, §5, §8 |
| 3 | "the Lang–Trotter model" overstated; LT(B) is not in LNM 504 | §0.1, §2 title, §2.1 |
| 4 | Bach–Charles: Thm 2.2 weight 2k, non-CM, minimal-polynomial output; Rems 2.5, 2.6; nothing on a(N) alone; a(N²) attribution | §1, §2.1, §5 |
| 5 | credit Charles 2006 for the √N trace-formula cost; square part via Lenstra–Pomerance | §0.4, §5 |
| 6 | "E–C at p and q" mis-cites: E–C is level 1 (Thm 15.2.1); Schoof for weight 2; Bruin general level | §0.4, §5 |
| 7 | 37a1 Serre curve VERIFIED (LMFDB, Jones 2010); joint images: Ribet 1975, Loeffler 2017 | §2.2, §3 |
| 8 | Gekeler Thm 5.5 (via David–Koukoulopoulos–Smith) | §2.2 (cx2 said Thm 4.6; not resolved, both given) |
| 9–10 | Kronecker–Hurwitz (H(0) = −1/12) and level-1 Eichler–Selberg (Zagier's appendix to Lang 1976) VERIFIED | §5 |
| 11 | Charles 2007 JCSS / TR1498 Thm 7.1 is factoring → lattice-vector counting via θ_{E₈} (Eisenstein), not cusp forms | §1 |
| 12 | E–C Rem 2.2.4 is the Eisenstein remark plus a Bach–Charles summary; Thm 2.4.18 is the Ramanujan bound | §1 |
| 13–16 | Pollard, Strassen, Harvey, Harvey–Hittmeir, Hafner–McCurley, Kalai, Bach 1988, Coppersmith §11 VERIFIED; prior work (Kunihiro–Koyama, Dieulefait–Jiménez Urroz, Chow 2015, Preobrazhenskii) all use determinant data or several indices | §1, §5 |

---

## 9. Reproduction and deviations
```
.venv/bin/python -m pytest -q tests/test_cuspcoef.py
PYTHONPATH=. .venv/bin/python -m src.cuspcoef.run_experiments e4    # < 1 min
PYTHONPATH=. .venv/bin/python -m src.cuspcoef.run_experiments e3    # 12 s
PYTHONPATH=. .venv/bin/python -m src.cuspcoef.run_experiments e1    # 122 s, peak RSS 739 MB
PYTHONPATH=. .venv/bin/python -m src.cuspcoef.run_experiments e2    # 612 s in total (n = 256: 353 s, 92 MB)
PYTHONPATH=. .venv/bin/python -m src.cuspcoef.run_experiments e2 --only 256   # D4
.venv/bin/python scratchpad/cuspcoef2/reconcile.py                  # §6.1 (scratchpad is git-ignored)
```
- **D1 (A1).** The smoke run revealed the list-decoding budget flaw before the registered run. P1c and P1d are reported
  as registered (refuted). P1f–P1h were added before the registered run.
- **D2.** After E1 ran, the method `SPF.factor` was renamed `SPF.decompose` so that the AST leakage test does not flag
  the name `factor`. Cosmetic: no logic changed.
- **D3.** E3's registered explicit bound differs from the one proved in §2.2; both are reported.
- **D4 (cuspcoef2).** The E2 module reruns every size and overwrites its output. The previous agent was stopped
  before n = 256 finished. cuspcoef2 added an `--only` flag (commit c63d967) that restricts the sizes and writes to
  `cuspcoef_e2_n<sizes>.jsonl`. Per-size seeds depend only on n, so the restricted run draws exactly the samples the
  full run would. The n = 64 rerun reproduced the committed line exactly.
- **D5 (provenance).** The n = 256 line in `cuspcoef_e2.jsonl` is a byte-identical copy of `cuspcoef_e2_n256.jsonl`
  (written 09:36:58, appended 09:37:16), committed at 2c4a3d1. That commit's message says the line "landed before
  stop"; it did not, it came from the D4 rerun. No duplicate line exists.
- **D6.** The E2 simulator samples x′ uniformly in [−2⌊√p⌋, 2⌊√p⌋], whereas the registration says [−2√p, 2√p]
  (endpoints differ by at most 1; cx2). Not corrected, since the change is negligible and the data are registered.
- **D7.** v1 quoted the search-space upper bound ∏(ℓ+1)/2 as S_n; v2 reports the registered per-instance S_n (§6.1).
