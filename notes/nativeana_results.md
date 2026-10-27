# nativeana: the analytic precision wall, mapped (phase 29)

Agent `nativeana`, 2026-10-27.

**Files.**
- Pre-registration: `notes/nativeana_prereg.md`, commit `240730b`, made before any code.
- Code: `src/nativeana/{chi,hankel,harness}.py`. `chi.py` and `hankel.py` are operational and see N or m only.
  `harness.py` holds the ground truth: PARI class numbers, which factor D internally.
- Tests: `tests/test_nativeana.py`.
- Raw data: `experiments/raw/nativeana_E{1,2,3}.jsonl`.
- GPT-6.1-sol blind round: `codex/nativeana_1.{md,out}`. The full report is at
  `codex/scratch/nativeana_1/nativeana_1_report.md` (git-ignored).
- Adversarial round and critics: §9.

**Labels.**
- **KNOWN**: with a reference.
- **PROVED**: the proof is given here.
- **SKETCH**: a proof outline with a named gap.
- **HEURISTIC**, **OBSERVATION**, **experimentally supported**, **UNVERIFIED**.
- "No prior statement found by us" is the strongest priority claim made.

## 0. Verdict (plain language)

**No analytic method we found beats N^{1/5}, and none looks capable of it. For one restricted but natural class of
methods we can now prove a barrier.** This is the outcome of an audit, not a general lower bound.

- **What the wall looks like (OBSERVATION).** Analytic invariants that factor N fall, as far as we found, into three
  regimes.
  - **Local invariants** (σ(N), r₄, the Kronecker–Hurwitz sum A = 2(N + q)). The factor sits in the local factor at the
    prime p ≈ √N. Analytic expansions (singular series, Ramanujan/circle method) are blind below modulus p, and an
    absolute error of N^{1/4} (relative N^{−3/4}) is needed even with a Coppersmith finish.
  - **Global invariants** (h(−4N), hR, L(1, χ_{−4N}), r₃). These can be approximated from small-prime data. In every
    estimator we ran (E1) the error falls like K^{−1/2} in the number K of terms. Then, near the approximate-functional-
    equation length ≈ 2√(4N), it decays rapidly to below rounding precision.
    - The group structure turns precision Δ into a √Δ search, which balances at N^{1/5} (Shanks/Lenstra/Schoof).
    - Relational algorithms (Hafner–McCurley L[1/2]) leave the analytic family.
  - **Coarse conductor readings** (Booker–Hiary–Keating). Cheap, but for a squarefree semiprime they only see the
    exponent pattern of N, which is public.
- **New barrier (Theorem LIN, PROVED; §4.1).** Consider estimators of L(1, χ_{−4m}) that are *linear* in the Kronecker
  symbols χ_{−4m}(n), under these hypotheses:
  - K ≤ Q^{1/4−ε} terms, at positions n ≤ P with PK ≤ Q^{1−ε};
  - weights smooth in m at scale Q, of size ≤ Q^{A₀};
  - the same terms for every m (non-adaptive; a "restricted quadrature" model).

  Each such estimator errs by ≥ c/√K for some m ≡ 1 (mod 4) in every dyadic range [Q, 2Q], i.e. by ≥ c√(Q/K) on
  h(−4m).
  - **Scope.** The hard m need not be a semiprime. The semiprime version is a SKETCH under GRH (§4.3). Euler products
    are log-linear and covered only by a GRH sketch (§4.2), whose rate c/(√K log K) is attained by the truncated
    product.
  - **Consequence (§4.4, model statement).** *If* the finishing step is charged Θ(√Δ) (BSGS/kangaroo), no estimator in
    the model plus that finish beats Q^{1/5}. This is not a factoring lower bound:
    - the finish is imposed, not proved optimal;
    - relational and hybrid methods escape the model (de Haan–Jacobson–Williams compute R in expected D^{1/6+ε} under
      GRH, from a relational candidate);
    - nonlinear and adaptive analytic estimators are not covered.
- **Coefficient extraction (Prop HW, PROVED; §4.5).** It is an elementary Hankel-rank consequence in the AAK/Prony
  circle of ideas. A *linear* rule that uses T point values or jets of a black-box power series cannot separate
  coefficient N from its neighbours within distance T: its moment error is ≥ 1 in ℓ¹ on the window. This formalises
  "smoothing hides the jump" in the black-box model only. Modular forms are not black boxes, and the model says nothing
  about structured algorithms.
- **Amortisation (§3).**
  - Swapping sums in Σ_t H(4N − t²) produces a singular series. For smooth macroscopic weights, its ℓ-factor depends on
    N only through v_ℓ(N) (the classical determinant density). So termwise local-factor evaluation is a divisibility
    test at ℓ = p.
  - The remainder is dominated, for our smooth weight, by the hyperbolic trace-formula term at trace p + q (E3:
    R² ≈ 0.99). Other weights can make cusp terms such as τ(p)τ(q) equally large.
  - We found no faster aggregate evaluation than the weighted hull method (N^{1/3}, codex_1). No computational lower
    bound for amortisation is claimed.
- **Which "fast" algorithms use the conductor's factorisation.**
  - Hiary's power-full method gains only on the square-full part. It permits factoring within its own cost and gives
    no gain for squarefree q.
  - Vishe needs a divisor of size ≈ q^{1/6} (and treats GL(2) twists).
  - Hiary's amortised quadratic L-function algorithm sieves the factorisations of a dense conductor window.
  - Platt's all-character transform uses the prime-power decomposition.
  - Rubinstein's smoothed AFE, Booker 2006, Odlyzko–Schönhage and Lagarias–Odlyzko use no factorisation. None beats
    √(conductor) for a single conductor.
- **Booker–Hiary–Keating (Duke 2015)** is the one analytic algorithm we found that reads an invariant at conductor ≍ N
  in *conjecturally* subexponential time, exp((log N)^{2/3+o(1)}).
  - Even under GRH, θ* < 1 is not proved, and the cost exceeds that of standard factoring algorithms.
  - It certifies squarefreeness. That information is public for semiprimes.

**Nothing here beats N^{1/5}, nor the project's internal (unrefereed) rigorous L[2/5] (Theorem U), nor the published
rigorous L[1/2] (Lenstra–Pomerance).** The single quantitative slack in the analytic class-number route is the
finishing step.
- A generic lower bound of Sutherland's Ω(Δ^{1/3}) type (thesis Thm 2.3, transferred to a window) is only a *candidate
  model* here. It would leave room down to Q^{1/7} in the architecture.
- Whether any finish beats √Δ is open (§7, item 1).

## 1. The wall table

Notation:
- Q = |D| = 4N, with N ≡ 1 (mod 4), so that D = −4N is fundamental. For N ≡ 3 (mod 4) use −N plus the public 2-part
  (N = 15: χ_{−60} is induced from conductor 15).
- Δ: absolute error.
- "Finish": the cheapest known way from that precision to a factor.
- Coppersmith finishing for p + q needs a separation q − p ≍ √N: error E in s = p + q moves q by E·q/(q − p)
  (codex_1).

| Invariant | Factor link | Precision needed + finish | Best native algorithm: unconditional / GRH | Speedup uses factorisation of the conductor? |
|---|---|---|---|---|
| h(−4N) | exact h ⇒ exponent ⇒ Cl[2] ⇒ p (randomised; GRH for sampling; genusleak §2.2, §3.3) | Δ < ½ rounds. Any Δ: window finish √Δ per class | Unconditional: AFE / certified evaluation N^{1/2+o(1)}; Booker 2006 terminates unconditionally. Randomised: Srinivasan's random summation (proved there for *real* quadratic fields, expected d^{1/5+ε}, with an ineffective constant). The imaginary version plus a √Δ finish is our adaptation (E1: retrospective score 0.25–0.29 at 28–44 bits). GRH: Shanks–Lenstra–Schoof N^{1/5}; Booker 2006 \|d\|^{1/4+ε}; Hafner–McCurley L[1/2, √2] (relational, under GRH); 3/√8 for a modified algorithm (Vollmer; Biasse–Erukulangara 2023) | No. Booker's explicit Burgess is for *prime* discriminants. Burgess with r ≥ 4 needs a cube-free modulus, which a semiprime promise supplies |
| L(1, χ_{−4N}) | h = (2√N/π)L | Δ_L = πΔ_h/(2√N): relative √(N)^{−1} to round | as h | AFE needs the *conductor*. For D = −4N with N ≡ 1, 2 (mod 4) it is \|D\| iff N is squarefree; for N ≡ 3 (mod 4) it is N (squarefree N). Squarefreeness is promised here (BHK test it otherwise) |
| R(4N), hR | a certified multiple of R reaches the midpoint (Prop HR, genusleak). It gives a factor only when N(ε) = +1 and the midpoint is useful: an observed 64–97%. N = 731 reflects the public factor 2 | window Δ ⇒ √Δ infrastructure search | AFE N^{1/2}; GRH: N^{1/5} (Lenstra 1982; Schoof 1982); Buchmann; de Boer–Pellet-Mary–Wesolowski, arXiv 2512.01588 (probabilistic subexponential under ERH, any number field; relational). de Haan–Jacobson–Williams, Math. Comp. 76 (2007): computes R with unconditionally certified output in expected D^{1/6+ε} under GRH (a relational candidate plus verification), outside the analytic model | No |
| single H(4N − t²) | none known (discriminant unrelated to p, q individually) | — | as h(4N − t²) | — |
| A = Σ_t H(4N − t²) = 2(N + q) = Σ_{ab=N}max(a, b) over ordered pairs (Kronecker–Gierster–Hurwitz; refined combinatorially by Popa–Zagier, PAMS 145 (2017)) | q = A/2 − N | Δ_A ≤ 2N^{1/4−ε} (Coppersmith on q, uniform) | term-wise via class-number algorithms: ≈ √N·N^{1/5} (a cost of that algorithm, not a lower bound); aggregate: A = 2(F(N) − F(N−1)) with F(x) = Σ_{ab≤x, a≥b} a, by weighted hull counting N^{1/3+o(1)} (codex_1, PROVED there; Sladkey's method with weights) | Hiary 2018 amortised quadratic L: yes (sieved factorisations) |
| Σ_t P_k(t, N)H(4N − t²), k = 4, 6, 8, 10, 14 | = −2(1 + p^{k−1}) (Eichler–Selberg, S_k(1) = 0; E3, nativemod E1a) | p^{k−1} to error ≲ p^{k−2}N^{1/4} | as A | — |
| r₂(N) | ∈ {0, 16}: one QR-type bit; two representations factor | the representations themselves | O(√N) squares; disc hull N^{1/3} | — |
| r₃(N) | 12h(−4N) (N ≡ 1 mod 4); 24h(−N) (N ≡ 3 mod 8) | as h | as h | as h |
| r₄(N) = 8σ(N) | s = p + q = σ − N − 1 | error 8N^{1/4−ε} (with separation) | hyperbola / hull weighted counting N^{1/3}; Li 2026 N^{9/28} for D(N) is circular (factors N); circle method: blind below modulus p (barrier R) | — |
| σ(N) = S(N) − S(N−1) | as r₄ | N^{1/4} on each of S(N), S(N−1) ≈ N² (relative N^{−7/4}) | N^{1/3} (Sladkey; Tao–Croot–Helfgott only *sketch* x^{1/3+o(1)} for Σ τ, §2.1) | — |
| C_N(X) = D_N(X) − D_{N−1}(X) | #{d ≤ X : d \| N}, bisection | exact (< ¼) | N^{1/3} (Hiary, Math. Comp. 85 (2016); hull walk) | — |
| θ(τ) values | carry no factor information; coefficients do | — | polylog by modular transformation (codex_1 §5) | — |
| rep(x)/N phase θ (leakmap GAR) | θN = 1 + jp | (3/4)log₂N bits ⇒ error N^{1/4} | no native estimator | — |
| log conductor (BHK) | squarefree vs not | coarse: O(1) in log | GRH + Conj 3.1: exp((log N)^{2/3+o(1)}); θ* < 1 unproved even under GRH | No; blind for squarefree N |

**The three regimes (OBSERVATION; the proved pieces are §§3–4).**
1. **Local invariants** (σ, r₄, A, φ-type values). Their analytic expressions are singular series or Ramanujan
   expansions, Σ_q c_q(N)a_q. The factor sits at modulus q ∈ {p, …}.
   - Terms below p are identical for N and for any N′ ≡ N modulo lcm(1..p−1) (barrier R, nonsmooth §4.1).
   - Required precision: N^{1/4} absolute, i.e. relative N^{−3/4}.
   - The circle method must at least reach modulus p, since all terms below p are blind. (An earlier draft said
     "N^{3/4} moduli", from the absolute tail bound. That ignores Möbius cancellation and is withdrawn.) As a Fourier
     coefficient at index N, extraction pays Prop HW.
2. **Global invariants** (h, hR, L(1, χ), r₃). These are smooth averages of public symbols with no blindness.
   - In the restricted linear model the error is ≥ c/√K (Theorem LIN, K ≤ Q^{1/4−ε}). Every estimator we ran shows
     K^{−1/2} until near the AFE length ≈ 2√Q.
   - There the truncation tail decays rapidly (E1: the "AFE cliff"). The AFE is exact as an infinite series; its
     truncation tail is what becomes negligible.
   - The group structure turns precision Δ into a √Δ search, so the known wall is N^{1/5}.
   - Relations (Hafner–McCurley) leave the analytic family altogether.
3. **Coarse conductor readings** (BHK). Cheap, but they see only the exponent pattern of N.

## 2. Literature audit: which speedups need a factored conductor

Each statement is about the cited algorithm, not about every possible modification.

- **Hiary, "Computing Dirichlet character sums to a power-full modulus"**, JNT 140 (2014) 122–146 (arXiv 1205.4687),
  Thm 1.1.
  - The cost is about ∏_{ℓ^a‖m} ℓ^{⌈a/3⌉}, plus preprocessing. That is q^{1/3} only for cube-full exponents: q = p²
    gives q^{1/2}.
  - Thm 1.1 assumes the factorisation, but Hiary notes this "is not essential": factoring in q^{1/3+o(1)} (Lehman)
    fits inside the cost. It uses Postnikov's formula χ(1 + ℓ^b x) = e(poly(x)/ℓ^{a−b}) at each prime power.
  - The gain is only on the square-full part. For squarefree pq the cost factor is pq = N: no gain (stated in the
    paper).
  - KNOWN; full text opened by codex_1, and the abstract verified by us.
- **Hiary, "A fast amortized algorithm for computing quadratic Dirichlet L-functions"**, QJM 69(1) (2018) 1–12,
  DOI 10.1093/qmath/hax029.
  - Õ(Q) preprocessing, then Õ(d(q)) per conductor in a window. It uses conductor factorisations, cheaply obtained by
    sieving the whole window.
  - Our family {t² − 4N} has √N members spread over conductors ≍ N. For it the preprocessing alone is ≍ N.
  - KNOWN (codex_1; verified by the literature critic: Thm 1.1 for L(1/2 + it, χ), conductors in [Q, Q + Δ),
    factorisations by a generalised sieve of Eratosthenes).
- **Vishe, "A fast algorithm to compute L(1/2, f×χ_q)"**, JNT 133 (2013) 1502–1524, Thm 1.
  - It takes q = MN with M ≤ N, M = M₁M₂, M₁ | N and (M₂, N) = 1. The cost is O(min(M⁵ + N, q)^{1+ε}), which gains
    only if q has a factor below q^{1/5}. M ≈ q^{1/6} gives q^{5/6+o(1)}.
  - It **needs that useful divisor**, not a full factorisation. A balanced semiprime has none of that size. f is a
    fixed GL(2) cusp form, so this is not a quadratic-character L-function.
  - KNOWN (abstract verified by us; Thm 1 per codex_1).
- **Platt, "Numerical computations concerning the GRH"**, Math. Comp. 85 (2016), Lemma 4.1.
  - It computes all φ(q) characters by FFT, using primitive roots and the prime-power decomposition of q (needs the
    factorisation).
  - The amortisation is over characters, not over one quadratic value.
  - KNOWN (codex_1).
- **Rubinstein**, smoothed AFE (CUP 2005, §3). It needs coefficients, the root number and the conductor.
  - For quadratic χ the root number is public (1 for real primitive characters).
  - The conductor equals |D| iff the input is squarefree up to the 2-part.
  - Cost √(conductor)^{1+o(1)}. No factorisation, given the squarefree promise.
- **Booker, "Artin's conjecture, Turing's method and the RH"**, Exp. Math. 15 (2006). Rigorous FFT-in-height
  evaluation; about √(analytic conductor) coefficients; no factorisation. KNOWN (codex_1).
- **Booker, "Quadratic class numbers and character sums"**, Math. Comp. 75 (2006) 1481–1492.
  - It terminates unconditionally with the correct h, and runs in O(|d|^{1/4+ε}) under GRH.
  - It uses an explicit Burgess bound, stated for prime discriminants.
  - KNOWN (abstract verified by us).
- **Odlyzko–Schönhage**, Trans. AMS 309 (1988), Thm 1.1. It amortises over *heights* at a fixed conductor; there is no
  family of discriminants, and no factorisation is used. KNOWN.
- **Lagarias–Odlyzko**, J. Algorithms 8 (1987). π(x) in x^{1/2+ε}. A smoothed explicit formula plus a local correction
  interval: the smoothing is paid back exactly in a window, the same phenomenon as P2 (conn_ana) and Prop HW. KNOWN.
- **Booker–Hiary–Keating, "Detecting squarefree numbers"**, Duke 164 (2015) 235–275 (arXiv 1304.6937; read via ar5iv).
  - Under GRH, positivity of the explicit formula for L(s, χ_d) (twisted by auxiliary χ_q) gives a computable lower
    bound on log|Δ|, where d = Δℓ².
  - Prop 2.3 gives runtime exp((log N)^{θ*+o(1)}). Their Conj 3.1 (random-matrix heuristics) predicts θ* = 2/3.
  - It certifies squarefreeness (or exposes a square factor) **without factoring**. This is the only analytic algorithm
    we found that reads an invariant at conductor ≍ N in conjecturally subexponential time.
  - The authors' caveats: subexponential time needs θ* < 1, which "seems beyond the current technology to prove, even
    under GRH". The conjectured cost is worse than standard factoring. Their only θ*-free theorem is the RSA
    "not squarefull" certification (Thm 1.1).
  - Whether the twisted prime-sum machinery has *other* useful outputs is not excluded by our observation.
  - What it reads is log-conductor to precision O(log ℓ). For squarefree N that is log N exactly, which is public.
- **Srinivasan**, Math. Comp. 67 (1998) 1285–1308.
  - Random summation: unbiased estimates of L(1, χ) from Õ(ε^{−2}) random symbols. Plus algebraic finishing gives an
    expected d^{1/5+ε}, unconditionally, for real quadratic fields. The abstract notes that the estimate is "not
    effective".
  - Theorem numbers UNVERIFIED (codex_1 could not open the primary paper). Our E1 implementation reproduces the ε^{−2}
    law.
- **Hafner–McCurley**, JAMS 2 (1989) 837–850: L[1/2, √2] under ERH for imaginary class groups, through smooth norms
  (relations; family M). 3/√8 for a modified algorithm: Biasse–Erukulangara, Adv. Math. Commun. 17 (2023); earlier
  Vollmer.
- **Shanks's N^{1/5} under ERH.** Shanks, PSPUM 20 (1971) 415–440, as credited by Harvey (Math. Comp. 90 (2021)).
  Sutherland's thesis (MIT 2007) §11.4 describes it (primes < d^{1/5}); verified by the literature critic.
- **Additional references** (literature critic): Zagier, C. R. Acad. Sci. Paris 281 (1975) 883–886 (the class-number
  Eisenstein series ℋ); Hirzebruch–Zagier, Invent. Math. 36 (1976); Bach, "Improved approximations for Euler
  products", CMS Conf. Proc. 15 (1995) (truncated Euler-product error); Coppersmith, J. Cryptology 10 (1997);
  Schoof, "Quadratic fields and factorization", MC Tract 155 (1982).

## 3. Amortisation (question 2)

### 3.1 Prop AM: the swapped-sum main term (density formula KNOWN; reading ours; E3 P3.2)
Let w be a smooth macroscopic weight in u = t/(2√N), and let M_x(N) be the main term of Σ_t w(t)H(4N − t²) obtained by swapping sums and replacing every
incomplete t-sum by its complete average, keeping the local factors at ℓ ≤ x.

*Claim.* The ℓ-adic factor of the full t-average depends on N only through v_ℓ(N).
- This needs smooth macroscopic averaging. A weight that selects residues of t breaks it. Codex_2's example: t ≡ 0
  (mod 3) gives 6 versus 12 trace-zero matrices mod 3 with det 1 versus det 2, with equal valuations.
- The density of matrices M ∈ M₂(Z_ℓ) with det M = N, summed over traces, is
  #{M mod ℓ^k : det M ≡ N}/ℓ^{3k} = (1 − ℓ^{−2})Σ_{j=0}^{v_ℓ(N)} ℓ^{−j}, for v_ℓ(N) < k.
- This is the classical Siegel local density of ad − bc = N; its product is σ(N)/(Nζ(2)) (KNOWN; e.g. the
  Duke–Rudnick–Sarnak counting circle, cited from memory by the critic).
- *Proof (Hensel).*
  - A matrix that is primitive, i.e. not ≡ 0 mod ℓ, has a non-zero gradient of det mod ℓ. So its lifts are
    equidistributed over det values: each primitive class mod ℓ with det ≡ a mod ℓ has ℓ^{3(k−1)} lifts with
    det ≡ a mod ℓ^k.
  - Primitive matrices mod ℓ number ℓ(ℓ² − 1) with a fixed unit det, and (ℓ² − 1)(ℓ + 1) with det ≡ 0.
  - Non-primitive M = ℓM′ have det = ℓ²det M′, so they recurse with a/ℓ².
  - Summing gives ℓ^{3k}(1 − ℓ^{−2})Σ_{j≤v}ℓ^{−j}.
  - The link from H(4N − t²) to matrix (t, N)-densities, with the ℓ-factor of H as the local orbital density, is
    Gekeler (IMRN 2003 no. 37, 1999–2018) for prime determinant. For composite N the analogue is assumed here, and
    is consistent with the exact identities of E3.
- E3 P3.2 checks the count exactly for ℓ ∈ {3, 5, 7} and k ≤ 3: 0 failures in 118 residues. The critic re-checked
  ℓ = 2 (k ≤ 5), 3 (k ≤ 4), 5 and 11 (k ≤ 2).
- "Only through v_ℓ(N)" refers to this trace-summed density. The t-weighted L(1)-moment factor at ℓ ∤ N also depends
  on the public (N/ℓ).
- Codex_1 independently: Σ_{t mod ℓ}((t² − 4N)/ℓ) = −1 if ℓ ∤ N and ℓ − 1 if ℓ | N. For the raw L(1)-moment the odd-ℓ
  factors are 1 − (1 + (N/ℓ))/(ℓ(ℓ² − 1)) (public) for ℓ ∤ N, and 1 + 1/ℓ for ℓ | N.

*Consequence.*
- M_x(N) equals its value for a prime N as long as x < p.
- The factor at ℓ = p differs from the generic one exactly when p | N. So *termwise local-factor evaluation* performs a
  divisibility test at each ℓ.
- This is not a lower bound: an aggregate algorithm need not evaluate the factors separately.
- Note that 2N·∏_{ℓ|N, ℓ≤x}(1 + 1/ℓ) equals 2(N + q) for p ≤ x < q, exactly the Kronecker–Hurwitz value. This is
  algebra: N(1 + 1/p) = N + q.
- So the amortised analytic main term carries the factor **only through the local factor at modulus p ≈ √N**.

### 3.2 The remainder: hyperbolic and identity terms (identity from Eichler–Selberg; E3)
Expand an even weight w(u), u = t/(2√N), in Chebyshev polynomials U_m. Since U_m(t/2√N) = N^{−m/2}P_{m+2}(t, N),
Eichler–Selberg at every weight gives

  S_w := Σ_t w·H(4N − t²) = 2c₀σ(N) − 2Σ_{m even} c_m N^{−m/2}[1 + p^{m+1} + Tr(T_N | S_{m+2}(1))],

with c_m the U-coefficients of w (rapidly decreasing for a C^∞ bump). Put A_w := (σ(N)/N)(1/π)Σ_t w·√(4N − t²); by
Euler–Maclaurin it equals 2c₀σ(N) up to a super-polynomially small error.
- The p-part sums to −2p·J_w(√(p/q)), where J_w(r) = Σ c_m r^m = (2/π)∫w(u)√(1 − u²)/(1 − 2ur + r²) du. This uses
  Σ U_m(u)r^m = (1 − 2ur + r²)^{−1}.
- The "1"-part is −2J_w(1/√N).
- Equivalently, the p-part is −(2√N/π)∫w(u)√(1 − u²)/(u₀ − u) du, a Stieltjes transform at u₀ = (p + q)/(2√N) > 1:
  the hyperbolic trace p + q, outside the elliptic range.
- The critic verified S_w **exactly** (rational Chebyshev coefficients up to m = 10, including τ(N)) for N = 77, 377,
  10403.
- So the factor appears twice at size √N:
  - in the main term (local factor at p);
  - in a transform of the test function at the hyperbolic class of trace p + q, which lies outside the elliptic
    range. (This is not literal evaluation of the bump outside its support.)
- Cusp traces, e.g. −2c₁₀N^{−5}τ(p)τ(q), are a third contribution. For our bump it is small. For other weights it is
  on the same √N scale: w = U₁₀ makes it the dominant factor-dependent term.
- The term −2J_w(1/√N) is O(1): for w = 1 it is exactly −2, and for the bump it tends to ≈ −0.516.
- **Prior art.** This is a re-reading of known structure, not a new identity:
  - Σ_{dd′=N} min(d, d′)^{k−1} is the textbook hyperbolic term of Eichler–Selberg/Selberg with general test functions
    (Zagier, LNM 627 (1977); Knightly–Li);
  - Chebyshev expansion across weights with the kernel (1 − 2ur + r²)^{−1} is the mechanism behind Serre's μ_q
    (JAMS 10 (1997)) and Conrey–Duke–Farmer (Acta Arith. 78 (1997));
  - swapping and Poisson-summing the elliptic sum, with class numbers handled through the AFE, is Altuğ (Compositio
    151 (2015), and parts II–III).
- **E3 (exploratory, as registered)** with the bump w = exp(−1/(1 − u²)), 60 semiprimes of 16–24 bits:
  - R_w = S_w − A_w has |R_w|/√N = 0.314 (median, all n);
  - the R² against −2pJ_w is 0.985 / 0.994 / 0.995;
  - max |R_w + 2pJ_w|/√N ≤ 0.0044. This is cusp traces, plus −2J_w(1/√N) (≈ −0.52, O(1)), plus discretisation.
    The size is consistent with Σ2|c_m|·dim S_{m+2}·d(N) ≈ 0.009.
- To factor from S_w one needs S_w to absolute error ≈ N^{1/4}. With independent per-term errors that means each of
  ≈ 4√N class numbers to O(1). Correlated errors could do better, and we have no bound for them. We did not check that
  p ↦ S_w is injective.

### 3.3 Other amortisation routes (assessment)
- **Coefficient of θ·ℋ** (Zagier's weight-3/2 mock Eisenstein series). A equals the 4N-th coefficient. Holomorphic
  projection (Mertens, Adv. Math. 301 (2016)) turns it into the σ − λ identity.
  - As an analytic extraction it is subject to Prop HW. As a structural computation it is σ(N).
  - We know no other route, but "every projection algorithm is circular" is not proved (codex_1's caution, adopted).
- **Odlyzko–Schönhage / Hiary-type amortisation over the family t² − 4N.**
  - Each member needs about √N AFE terms, and the family's symbols χ_{t²−4N}(n) are periodic in t mod n.
  - Summing over t first gives back the singular series of 3.1. Evaluated termwise, its ℓ-factors for ℓ ≤ √N decide
    ℓ | N for every such ℓ, which is trial division.
  - Hiary's amortised algorithm needs dense conductor windows and their factorisations.
  - No sub-N^{1/5} aggregate is known; the best verified aggregate is codex_1's weighted hull, N^{1/3}.
- **Verdict (HEURISTIC, supported by 3.1–3.2).** A fast evaluation of the elliptic side is a fast evaluation of the
  divisor/hyperbolic side. Every amortisation we tried re-derives the trace formula and gains nothing. No
  computational lower bound is claimed.

## 4. Model barriers (question 3)

### 4.1 Theorem LIN (PROVED; uses Heath-Brown's quadratic large sieve)
**Setting.**
- Let ε ∈ (0, 1/100) and let η ∈ C_c^∞((0, 1)) be fixed, with η ≥ 0 and ∫η = 1. Put Φ(m) = η((m − Q)/Q). Here Q is
  the scale of m and |D| = 4m (in §1, Q = |D|). Below, φ(·) is Euler's function.
- For m ≡ 1 (mod 4) let χ_m = χ_{−4m} (Kronecker; possibly non-fundamental) and L_m = L(1, χ_m) = Σ_{n odd} χ_m(n)/n.
- An **admissible linear estimator** is L̂(m) = Σ_{n∈S} W_n((m − Q)/Q)·χ_m(n), where:
  - S is a set of odd integers in [1, P] with |S| ≤ K;
  - each W_n ∈ C^∞ satisfies sup_{[0,1]}|W_n^{(j)}| ≤ B·C_j for all j, for a fixed sequence (C_j) and B ≤ Q^{A₀}.

**Theorem LIN.** There are c > 0 and Q₀, depending on ε, A₀, (C_j) and η, with the following property. For Q ≥ Q₀,
if K ≤ Q^{1/4−ε} and PK ≤ Q^{1−ε}, then

  Σ_{m≡1(4)} Φ(m)·|L_m − L̂(m)|² ≥ (c/K)·Σ_{m≡1(4)} Φ(m).

Hence some m ∈ (Q, 2Q), m ≡ 1 (mod 4), has |L_m − L̂(m)| ≥ √(c/K). Since h(−4m) = (2√m/π)L_m, every estimator
ĥ = (2√m/π)L̂ (the factor is smooth in m and absorbed into W_n) has |h(−4m) − ĥ(m)| ≥ c′√(Q/K) at that m.

*Proof.*
1. **Lemma 1 (smooth character sums; Poisson).** (The Poisson terms carry phases e(±hQ/4r), omitted below; only
   their absolute values are used. A Fourier-decay exponent A′ is chosen large enough to give the saving Q^{−A}
   against B ≤ Q^{A₀}.)
   - Let r be odd with r ≤ Q^{1−ε/2}. Let G = η·W with W as above, and ψ_r(m) = (−m/r)·1[m ≡ 1 (4)].
   - Splitting m by residue mod 4r and applying Poisson summation gives
     Σ_m G((m − Q)/Q)ψ_r(m) = (Q/4r)Σ_h Ĝ(hQ/4r)·τ_h, where |τ_h| ≤ r, τ_0 = Σ_{a mod 4r}ψ_r(a), and
     |Ĝ(ξ)| ≪_A B(1 + |ξ|)^{−A}.
   - Since Q/4r ≥ Q^{ε/2}/4, the terms h ≠ 0 total O_A(B·Q^{−A}).
   - If r is not a square, m ↦ (m/r) is a non-principal character mod r, so τ_0 = 0.
   - If r = s², ψ_r(m) = 1[(m, s) = 1]·1[m ≡ 1 (4)], so the sum is (Q/4)(φ(s)/s)∫G + O(BQ^{−A}).
   - For odd n we have χ_m(n) = (−m/n). So χ_m(n)χ_m(k) = ψ_{nk}(m).
2. **Tests.**
   - Let 𝒦 be the set of odd squarefree k ≤ 8K such that nk is not a square for every n ∈ S.
   - The odd squarefree k ≤ 8K number ≈ (4/π²)·8K > 3K.
   - At most K squarefree k are excluded: k is excluded only if it is the squarefree kernel of some n ∈ S.
   - On (4K, 8K] there are ≈ 1.62K odd squarefree k, and the mean of φ(k)/k over them is ≈ 0.86. After the
     exclusions, Σ φ(k)/k over the remaining k in (4K, 8K] is > 0.3K for large K.
   - Each such k has 1/k² ≥ 1/(64K²), so Σ_{k∈𝒦}(φ(k)/k)/k² ≥ c₁/K. Numerically the exact worst case gives
     K·Σ ≈ 0.095 for all K ≤ 3000 (critic).
3. **Near-orthogonality.**
   - Put ψ_k := ψ_k(m) = χ_m(k) and G_{kk′} = Σ_m Φψ_kψ_{k′}.
   - For k ≠ k′ in 𝒦, kk′ ≤ 64K² ≤ Q^{1/2} is not a square, so |G_{kk′}| = O(Q^{−A}).
   - G_kk = (Q/4)(φ(k)/k) + O(Q^{−A}) ≫ Q/log log Q.
   - With u_k = ψ_k/√G_kk, the Gram matrix of the u_k is I + O(KQ^{−A}).
   - For E := L − L̂ the frame bound gives Σ_{k∈𝒦} |⟨E, ψ_k⟩|²/G_kk ≤ (1 + o(1))·Σ_mΦ|E|².
4. **The estimator is orthogonal to the tests.** For n ∈ S and k ∈ 𝒦, nk is not a square and nk ≤ 8PK ≤ 8Q^{1−ε}.
   Lemma 1 with G = ηW_n (8Q^{1−ε} ≤ Q^{1−ε/2} for Q ≥ Q₀) gives |⟨L̂, ψ_k⟩| ≤ K·B·Q^{−A}, with A > A₀ + 2.
5. **The target is not.**
   - Put X = Q^{1−ε}/(8K) and split L_m = Σ_{n≤X} χ_m(n)/n + T_m.
   - In the first part, terms with nk not a square are O(Q^{−A}) by Lemma 1.
   - Terms with nk a square are n = kj². Each is (1/kj²)Σ_mΦ·1[(m, kj) = 1] ≥ 0, and j = 1 contributes G_kk/k (note
     k ≤ 8K ≤ X).
   - So the first part is ≥ G_kk/k − O(Q^{−A}).
   - Second part: |Σ_mΦψ_kT_m| ≤ (Σ_mΦ)^{1/2}·‖η‖_∞^{1/2}(Σ_{m≤2Q}|T_m|²)^{1/2}.
6. **Lemma 2 (tail mean square).** Σ_{m≤2Q, m≡1(4)} |T_m(X)|² ≪_ε Q^{ε/2}(Q/X + Q^{1/2}).
   - For n > Q², Pólya–Vinogradov and partial summation bound the tail by O(√Q log Q/Q²).
   - For a dyadic block n ∈ (Y, 2Y], X ≤ Y ≤ Q², write n = n₁s² and m = m₁t² (n₁, m₁ squarefree). Then
     χ_m(n) = (−m₁/n₁)·1[(m₁t², s) = 1]·1[(n₁, t) = 1].
   - For fixed s and t, quadratic reciprocity (split by residues mod 4) turns the inner sum into
     Σ_{m₁≤2Q/t²}|Σ_{n₁~Y/s²} a_{n₁}(n₁/m₁)|², with |a_{n₁}| ≤ s²/Y and both variables odd squarefree.
   - Heath-Brown (Acta Arith. 72 (1995), Thm 1: Σ*_{m≤M}|Σ*_{n≤N} a_n (n/m)|² ≪_ε (MN)^ε(M + N)Σ|a_n|²) bounds it by
     ≪ Q^{ε/8}(Q/t² + Y/s²)·s²/Y.
   - Summing over t ≤ √(2Q) gives Q^{ε/8}(Qs²/Y + √Q).
   - Since m₁ ≡ 1 (mod 4), reciprocity has no sign: (−m₁/n₁) = (−1/n₁)(n₁/m₁).
   - Bookkeeping (codex_2). For fixed s and t, put 1[(n₁, t) = 1] and (−1/n₁) into the coefficients, keep
     1[(t, s) = 1], and drop the outer restriction (m₁, s) = 1 only after taking absolute squares. This gives
     ‖B_Y‖₂ ≪ Q^{o(1)}Σ_{s≤√(2Y)} s^{−2}(Qs²/Y + √Q)^{1/2} ≪ Q^{o(1)}(√(Q/Y)·log 2Y + Q^{1/4}).
   - Minkowski over s (weights 1/s²) loses a log; Minkowski over the O(log Q) dyadic blocks loses log².
7. **Conclusion.**
   - By Lemma 2, with Q/X = 8KQ^{ε}, the second part is ≪ Q^{1/2}·(Q^{ε/2}(8KQ^{ε} + Q^{1/2}))^{1/2}
     ≪ Q^{1/2+3ε/4}√K + Q^{3/4+ε/4}.
   - Both are ≤ Q^{3/4+ε/4}·O(1) when K ≤ Q^{1/4−ε}.
   - Meanwhile G_kk/k ≥ cQ/(8K log log Q) ≥ cQ^{3/4+ε}/(8 log log Q).
   - So |⟨E, ψ_k⟩| ≥ G_kk/(2k) for Q ≥ Q₀, and step 3 gives
     Σ_mΦ|E|² ≥ (1 − o(1))Σ_{k∈𝒦} G_kk/(4k²) ≥ c·Q·Σ_{k∈𝒦}(φ(k)/k)/k² ≥ c·c₁·Q/K. ∎

*Remarks.*
- (i) The range K ≤ Q^{1/4−ε} is set by the non-squarefree m with large square part (the √Q in Lemma 2). Q^{1/4}
  already covers every K relevant to the 1/5 balance.
  - Restricting to squarefree m does **not** simply extend the range. Writing μ² = Σ_{d²|m}μ(d), the range d > √(Q/r)
    ruins Lemma 1's O(Q^{−A}) (critic). A new orthogonality lemma would be needed, and that is OPEN.
- (ii) Positions up to P = Q^{1−ε}/K are allowed. So any K-term sub-selection of the AFE's n ≤ 4√Q is inside the
  model, and the full AFE (K ≈ 4√Q > Q^{1/4}) is outside it, as it must be.
- (iii) Randomised estimators: the bound holds for each seed, so some m has E_seed|E_m|² ≥ c/K. That is a second
  moment and does not by itself bound the expected cost of a finish. Srinivasan's sampler
  uses positions up to Q^{1+}, which is outside the PK ≤ Q^{1−ε} hypothesis. Its own variance is ≍ (log Q)²/K
  anyway (E1).
- (iv) Smoothness of W_n in m at scale Q is essential. With unrestricted dependence on m a single weight W_1(m) = L_m
  is "linear" (codex_1's objection). Hypothesis (C_j, B ≤ Q^{A₀}) excludes that interpolation, and Lemma 1 is exactly
  where it is used.
- (v) Coverage, algorithm by algorithm:
  - inside: fixed-cutoff partial sums, smooth truncations, short sub-selections of the AFE, and short random samplers
    with cutoff ≲ √(QK) (Pólya–Vinogradov);
  - outside: Euler products at their prime-operation cost (only LIN′), the full AFE, Booker 2006 (algebraic work),
    Srinivasan's sampler (P ≍ Q^{1+}), input-dependent supports and stopping rules, and relational methods.

  This is a restricted quadrature theorem, neither vacuous nor a general lower bound.
- (vi) No prior statement found by us. The ingredients (square-class orthogonality of quadratic characters, the HB
  large sieve, Bessel) are standard in the L(1, χ_d) distribution literature (e.g. Granville–Soundararajan, GAFA 2003).

### 4.2 LIN′: Euler products (SKETCH; GRH)
- Estimators of the form log L̂(m) = Σ_{ℓ∈S} [a_ℓ(m)χ_m(ℓ) + b_ℓ(m)χ_m(ℓ)²], with S a set of primes, |S| = K, and
  smooth weights. This covers truncated Euler products and any re-weighting of them.
- The test argument applies with tests k = primes ℓ ∉ S:
  - ⟨log L, ψ_ℓ⟩ has main term G_ℓℓ/ℓ;
  - composite squarefree k carry no main term for log L;
  - χ_m(ℓ)² is principal-type and orthogonal to every ψ_k with k > 1.
- Result: some m has |log L_m − log L̂(m)| ≥ c·(Σ_{ℓ>p_K} ℓ^{−2})^{1/2} ≍ c/(√K·log K). The odd-power coefficient
  is atanh(1/ℓ), of the same scale. Converting to absolute L-error loses a factor 1/log log Q under GRH. The worst case is S = the
  first K primes. (Corrected after the critic: v1 said c/√(K log K), which is false.)
- The truncated Euler product attains this order, so LIN′ is sharp in its model.
- **Gap:** the tail Σ_{n>X}Λ(n)χ_m(n)/(n log n) needs the HB bound on dyadic blocks up to Q^A, and GRH for the
  pointwise tail beyond. Under GRH, L_m ≥ c/log log Q converts this to an error in L.
- E2's log-linear LS fits match the tail prediction: RMS/(Σ_{odd ℓ>p_K}ℓ^{−2})^{1/2} = 0.94–1.09 for all K = 2…512 in
  both families (post hoc, §6).

### 4.3 Semiprime family (SKETCH; GRH)
The factoring-relevant inputs are m = N = pq with p, q ≍ √Q.
- Lemma 1 over semiprimes: for a non-principal ψ mod r ≤ Q^{1−ε}, Mellin-separating η(pq/Q) gives
  Σ_{pq} η·ψ(p)ψ(q) ≪ Q^{1/2+ε} under GRH.
- The estimator term becomes Σ_{n∈S}|W_n|·Q^{1/2+ε}, which must be ≪ Q/K. That needs BK² ≪ Q^{1/2−ε}.
- Lemma 2 holds, and improves (t = 1 only, so there is no √Q term), since HB applies to any subset of squarefree m.
- Missing: the balanced region, pq ≡ 1 (mod 4), ordering and diagonal removal, imprimitive characters, the family
  mass ≍ Q/log²Q (rebuild the Gram bounds and target correlations with it), the semiprime count and main-term
  asymptotics, the Mellin separation with W_n, and the log² Q loss from
  semiprime density in the constant.
- Hence, under GRH, LIN holds on semiprimes for K ≤ Q^{1/4−ε}/√B.
- Not written in full; E2 family B agrees with family A to within a factor 1.15.

### 4.4 Corollary ARCH: a certified-window trade-off (PROVED given LIN, in the stated model)
**Model assumptions.**
- A K-term estimate costs ≥ K (one symbol per term). This is an assumption, not a theorem.
- S and W_n are fixed for all m (non-adaptive).
- The hypotheses of LIN hold, *including* PK ≤ Q^{1−ε}. Small K alone does not imply this.
- The estimate is linear (LIN). Log-linear Euler products need LIN′ (SKETCH, GRH), which gives the same exponent up to
  logs.
- The algorithm is deterministic and must finish from a window that is certified to contain h(−4m), for every m.
- The finishing cost F(Δ) is **imposed**. F(Δ) = Θ(√Δ) is the BSGS/kangaroo upper bound, not a proved lower bound
  for every finish.

**Statement.**
- The worst-case cost over m ∈ (Q, 2Q] is ≥ min_K [K + F(c′√(Q/K))], since by LIN some m forces Δ ≥ c′√(Q/K).
- If K > Q^{1/4−ε} the cost is already ≥ Q^{1/5}.
- With F(Δ) = √Δ the minimum is ≍ Q^{1/5}, at K = Q^{1/5}.

**Scope.**
- This is PROVED in this model for the class-number problem over all m ≡ 1 (mod 4). The hard m need not be a
  semiprime.
- For factoring inputs it is conditional on §4.3 (GRH, SKETCH). It is **not** a factoring lower bound.
- Randomised algorithms are not covered. A seedwise mean-square bound does not control expected finishing cost,
  because rare large errors can dominate a second moment.

**The finish (candidate model only).**
- Sutherland's thesis Thm 2.3 proves a generic Ω(N^{1/3}) lower bound for order computation, using exactly the
  divisibility count. He also disproved the conjectured Ω(√N) for unbounded order computation: primorial steps give
  o(√N).
- A window analogue would give Ω(Δ^{1/3}). But "a hidden prime component with a prescribed distribution" is a
  *candidate* model, not a theorem about class numbers in analytic windows.
- Conditional on a Δ^{1/3} finish, the balance is min_K max(K, (Q/K)^{1/6}) = Q^{1/7}. No such finish is known.

### 4.5 Prop HW: Hankel-rank uncertainty for coefficient extraction (PROVED)
**Statement.**
- Let m_n = Σ_{i=1}^{T} c_i z_i^n for any complex c_i and z_i. Jets are allowed: a node of multiplicity μ counts μ
  times in T.
- If N ≥ W and T ≤ W, then Σ_{|n−N|≤W}|m_n − δ_{nN}| ≥ 1.
- Consequently, take any linear rule a_N ≈ Σ c_i f(z_i) on (say polynomial) f = Σ a_nq^n, with T ≤ W.
  - It errs by ≥ A on some f whose coefficients are fixed outside [N − W, N + W] and complex with |a_n| ≤ A inside.
  - With *real* coefficients in [−A, A] the worst error is ≥ (2/π)A.
  - The real case cannot give A. Critic counterexample: W = T = 1, rule 0.2·i^{−N}f(i), worst error ≈ 0.894A. Codex_2
    gives another: N = W = T = 1, â₁ = −(8i/29)f(i/2), worst error 5A/√29.
  - Its moment ℓ¹ error, 35/29 > 1, is consistent with the theorem.
- The bound is tight: W + 1 equally spaced nodes on |z| = 1 reproduce δ_N exactly on the window. Both directions are
  checked in `tests/test_nativeana.py`.

*Proof.*
- Consider the (W+1)×(W+1) block H = [m_{j+k}] with 0 ≤ j ≤ W and N − W ≤ k ≤ N.
- H = V₁·diag(c)·V₂ᵀ has rank ≤ T (confluent Vandermonde for jets).
- The same block of δ_N is the anti-identity, with all singular values equal to 1.
- By Weyl, 1 = σ_{T+1}(anti-identity) ≤ σ_{T+1}(H) + ‖H − anti-identity‖ = ‖E‖.
- E = Σ_{|s−N|≤W} e_s·P_s, with each P_s a partial permutation matrix, so ‖E‖ ≤ Σ|e_s|. ∎

**Reading.**
- Smoothing, theta transformations and T-node quadratures of contour integrals, used as *linear* functionals of a
  black-box q-series, cannot isolate index N from its T neighbours. (An exact contour integral does isolate a_N; the
  point is the number of nodes.) This is P2 ("smoothing hides the jump"), now a theorem in the black-box
  linear model.
- It compares with codex_1's exact-recovery version (T > N; approximate error only 2^{−T}). The window version gives
  constant error at window size T.
- Escapes use structure that the model forbids:
  - modularity (relates values, not coefficients);
  - Hecke multiplicativity (a_N = a_p a_q: factoring; Bach–Charles);
  - for L-values, the functional equation (AFE).
- **Scope.** This is not a lower bound for modular coefficients. Theta series and modular forms have no independent
  coefficient variations, and the model ignores precision and modular relations. It is a statement about black-box
  linear information only.
- **Prior art.** This is an elementary finite consequence of the Hankel-rank / AAK mechanism: Kronecker (1881);
  Adamyan–Arov–Krein; Peller, *Hankel operators* (§6 of the MSRI survey); Beylkin–Monzón, ACHA 19 (2005) 17–48;
  Plonka–Pototskaia, arXiv 1609.09603. No prior statement of this exact window corollary was found by us.

### 4.6 What escapes the models (list)
- Relations (Hafner–McCurley, Buchmann). They are non-linear, use smooth norms, and run in L[1/2].
- Group finishing. BSGS uses the group, not the analytic estimate.
- Coarse conductor readings (BHK). These use positivity of zero sums, not linear estimation of L(1).
- Algorithms that are non-smooth in m beyond scale Q (anything that "knows" m arithmetically, e.g. Euclid on m, CF of
  √m) are outside LIN. The infrastructure (SQUFOF, N^{1/4}) is such a method.

## 5. Escapes (question 4)

1. **Conductor or level independent of N.**
   - r_k(N) = [q^N]θ^k, σ(N) = [q^N]Σ dq^d/(1 − q^d) and A = [q^{4N}]θℋ all have fixed level (codex_1 §5).
   - θ values are polylog-cheap, but coefficients are not. Black-box extraction pays Prop HW (window ≈ the index scale
     for exactness). The circle method is blind below modulus p.
   - Hecke structure is factoring (Bach–Charles; Bach–Miller–Shallit for σ).
   - So the "N as index" analytic route is worse than "N as conductor": no functional equation acts on the index.
     This is the hand-off to `nativemod` (algebraic side).
2. **Coarse precision plus cheap pinning.** Exact required precisions:
   - σ / r₄: 8N^{1/4−ε} (separation-qualified);
   - A: 2N^{1/4−ε} (uniform);
   - Σ(t² − N)H = −2(1 + p³): ≈ 3p²N^{1/4} ≈ 3N^{5/4};
   - h: any Δ, with a √Δ finish;
   - hR: likewise;
   - L(1, χ): Δ_h·π/(2√N).

   In the LIN model the cost of Δ_h is ≳ Q/Δ_h², so precision N^{1/2−δ} costs N^{2δ} and finishing costs N^{1/4−δ/2}: the
   balance is δ = 1/10.

   The local invariants need *relative* N^{−3/4}. We know no analytic method that gets there below the AFE length:
   the K^{−1/2} rate is the known behaviour below the AFE length, and it is proved only in the LIN model. They are
   also blind below p.
3. **Combination (BHK-type coarse reading of an asymmetric invariant).**
   - We looked for an L-function built from N alone whose conductor or Γ-factor depends asymmetrically on p and q. For
     squarefree pq, every conductor in the constructions we examined is a product N^a·(public): quadratic and cubic
     fields Q(√±N), Q(∛N) and Q(∛N²) (the same field), twists E ⊗ χ_N.
   - Asymmetry would need a construction that already distinguishes p from q (OBSERVATION; no theorem).
   - So coarse readings yield only exponent patterns: squarefreeness, not factors.

## 6. Experiments (pre-registered; outcomes)

**E1: precision versus cost** (n = 28 / 36 / 44, 20 instances each, N ≡ 1 mod 4, no skipped-seed issues).

| Prediction | Outcome |
|---|---|
| P1.1 EP slope ∈ [−0.65, −0.35] | **met**: −0.567, −0.583, −0.544 |
| P1.2 PS slope | **met**: −0.478, −0.479, −0.559 |
| P1.3 SM slope | **met**: −0.376, −0.419, −0.543 |
| P1.4 AFE exact for T ≥ 4√Q; AFE/PS ratio ∈ [1/3, 3] for T ≤ √Q/16 | **met**: 0 failures; ratios 0.92–1.01 |
| P1.5 RS slope ∈ [−0.6, −0.4]; constant ∈ [0.05, 0.5] | **met**: −0.518, −0.473, −0.477; 0.137, 0.110, 0.107 |
| P1.6 EP log C_min/log N ∈ [0.15, 0.30] | **met**: 0.202, 0.201, 0.188. Others: PS 0.20/0.20/0.20, SM 0.21/0.20/0.21, AFE 0.20/0.20/0.20, RS 0.29/0.26/0.25 (retrospective, uses the true error) |

- **The AFE cliff.** At 44 bits (√Q ≈ 8.4·10⁶) the median error is 1.4·10³ at T = 2^20, 108 at T = 2^22, 2.9 at T = 2^23,
  and 8·10⁻⁶ at T = 2^24.
- Below ≈ √Q/16 the AFE is indistinguishable from the sharp partial sum.

**E2: optimal linear rules** (Q = 2^30; 6000 train / 3000 test; family A: all m ≡ 1 mod 4; family B: balanced
31-bit semiprimes).

| Prediction | Outcome |
|---|---|
| P2.1 slope of test RMS vs K, both families and rules, ∈ [−0.65, −0.35] | **met**: A −0.498 / −0.487; B −0.495 / −0.486 |
| P2.2 rule (ii) (smooth in m) ≥ 0.9 × rule (i) | **met**: ratios 1.00–1.06 (no gain) |
| P2.3 B/A ∈ [0.7, 1.4] | **met**: 1.12–1.18 |
| P2.4 RMS·√K ∈ [0.1, 1.0] | **met**: 0.46–0.50 (A), 0.52–0.60 (B) |
| P2.5 log-linear slope ∈ [−0.65, −0.35] | **NOT MET**: −0.661 (A), −0.670 (B) |

- P2.5 post hoc (corrected after the critic): LIN′'s model predicts RMS ≈ (Σ_{odd ℓ>p_K}ℓ^{−2})^{1/2} ≍ 1/(√K log K).
  The measured ratio RMS/prediction is 0.94–1.09 for every K = 2…512, in both families.
  - The asymptotic local slope of the prediction is about −0.70.
  - An intermediate form, (K log K)^{−1/2}, which we used first, is wrong. RMS·√(K ln K) drifts from 0.375 to 0.318.
- The registered range was set from the linear case without the log factor. That is a specification error, recorded
  as a miss.
- The optimal LS rules sit at RMS ≈ 0.5/√K, the theorem's c/√K shape. Allowing weights to vary with m buys nothing.

**E3: Kronecker–Hurwitz** (n = 16 / 20 / 24, 20 each).

| Prediction | Outcome |
|---|---|
| P3.1 three Eichler–Selberg identities exactly | **met**: 60/60 |
| P3.2 local densities | **met**: 118 residues, 0 failures |
| P3.3 median \|R_w\|/√N ∈ [0.05, 20] | **met**: 0.316, 0.315, 0.314; R² vs hyperbolic 0.985 / 0.994 / 0.995 (reported) |

No factoring is claimed. All "successes" are identity checks.

## 7. Shortlist and NP-hardness notes

**Nothing below beats N^{1/5}.** Ranked by usefulness to the project:
1. **(Barrier, open; the only slack in the analytic route.)** The generic complexity of "compute ord(g), or the group
   exponent, given a window of width Δ containing |G|".
   - LIN + √Δ gives exactly 1/5. A Δ^{1/3} finish would allow 1/7.
   - First question: does a Sutherland-type argument (thesis Thm 2.3) with a window distribution give Ω(√Δ) for
     generic algorithms, or does a primorial-type construction achieve Δ^{1/2−c}? Sutherland's o(√N) primorial steps
     for unbounded order computation show that the answer is not obvious.
   - A matching Ω(√Δ) bound would make N^{1/5} a theorem for "linear analytic estimate + generic finish".
2. **(Barrier, cheap to finish.)** Write LIN′ (Euler products) and the semiprime version in full under GRH.
   Extending LIN beyond K ≤ Q^{1/4} needs a new orthogonality lemma for squarefree m (OPEN; §4.1 remark (i)).
3. **(Probe, low value.)** Nonlinear analytic estimators of L(1, χ_D) outside both models, e.g. polynomial in
   L-values of several twists χ_D·χ_k. In the random model the twisted errors are near-orthogonal, so we expect no
   gain (HEURISTIC).
4. **(Hand-off to nativemod.)** Index-N coefficients of fixed-level forms. Analytically closed by Prop HW and by
   blindness below p; algebraically their question.

**Hardness of natural generalisations.**
- Evaluating σ(n) for arbitrary n is randomised-polytime equivalent to factoring (Bach–Miller–Shallit, SICOMP 15
  (1986)). So computing r₄ is FACT, and so not NP-hard unless NP = coNP.
- Variable-dimension theta coefficients [z^B]∏(1 + z^{a_i}) are #P-complete (codex_1, a parsimonious #3SAT reduction).
  These are not our objects.
- Exact degree-2 floor sums are #P/NP-hard (floorhard). Class numbers themselves are in L[1/2] under GRH, so they are
  not hard in that sense.
- Natural generalisations of the analytic objects are thus either FACT-equivalent or far from our instances.
  Deprioritise none on NP grounds, but none is promising.

## 8. Coordination note for `nativemod` (notes only)
- E3 §3.2 gives the smoothed-weight version of their E1a. The elliptic side with any smooth weight is
  2c₀σ(N) − 2pJ_w(√(p/q)) − 2J_w(1/√N) + (cusp traces at N). The hyperbolic class t = p + q carries the factor
  explicitly.
- Prop HW is the analytic-side statement for their fixed-level coefficients. Linear black-box extraction from values
  needs window-size samples, so any sub-N method must use Hecke/modular structure. That is their side.

## 9. Reviews
- **GPT-6.1-sol blind round (codex/nativeana_1).** Independently:
  - the wall table;
  - the literature audit, including Hiary QJM 2018, de Haan–Jacobson–Williams, de Boer–Pellet-Mary–Wesolowski, and
    the classical Kronecker–Gierster–Hurwitz form A = Σ_{ab=N} max(a, b), refined by Popa–Zagier;
  - the weighted-hull N^{1/3} aggregate for A;
  - the local character sums;
  - an exact-recovery coefficient bound (T > N) and an independent-signs sampling bound (≥ 1/(K+1));
  - the Coppersmith separation caveat;
  - the #P-completeness remark.

  It did *not* find a bound for actual characters, and it flagged smooth-weight models as vacuous unless the smoothness
  is quantified. Theorem LIN answers this with the (C_j, B) hypothesis.
- **Opus proofs critic (read-only): ACCEPT WITH FIXES on all six items; no REJECT.**
  - LIN's proof was checked step by step and is sound. Fixes:
    - nk ≤ 8Q^{1−ε};
    - the test-set counting argument;
    - the ‖η‖_∞ factor;
    - the notation clash φ/η;
    - remark (i) withdrawn (the squarefree extension needs a new orthogonality lemma);
    - scope in §0 (non-squarefree hard m; semiprimes only in the sketch).
  - **WRONG, fixed:** LIN′ rate c/√(K log K) → c/(√K log K). The P2.5 post-hoc explanation is rewritten; the
    tail-sum prediction matches E2 at 0.94–1.09, which we re-ran.
  - **WRONG, fixed:** Prop HW's real-box corollary (counterexample) → complex coefficients, or a (2/π)A constant.
  - **WRONG, fixed:** §5 precision 3pN^{1/4} → 3p²N^{1/4}.
  - **GAP, fixed:** Prop AM now has a Hensel proof. §3.2 was verified exactly by the critic for N = 77, 377, 10403; A_w
    is defined; the −2J_w(1/√N) term is included.
  - Overclaims softened: "no method faster than K^{−1/2}", the N^{3/4} circle-method claim (withdrawn), "not NP-hard
    unless NP = coNP", the conductor hypothesis N ≡ 1, 2 (mod 4), and ARCH framing.
- **GPT-6.1-sol adversarial round (codex/nativeana_2): "LIN's core argument survives; the advertised factoring
  optimality does not follow."** All findings applied:
  - FATAL #8 (ARCH as advertised). ARCH is restated as a deterministic certified-window trade-off on the proved family,
    with the finish cost imposed and the position hypothesis restored. The semiprime and randomised extensions are
    labelled separately. §0 is rewritten.
  - #5: the model is called a restricted quadrature theorem, with algorithm-by-algorithm coverage.
  - #6: the LIN′ log (as the critic).
  - #7: the semiprime gaps are listed.
  - #9: Δ^{1/3} is a candidate model only, citing Sutherland Thm 2.3.
  - #10, #11: the HW counterexample is added; HW is not a modular-coefficient lower bound; AAK/Peller prior art.
  - #12: AM needs smooth macroscopic averaging (the t ≡ 0 mod 3 counterexample), and is reworded as "termwise
    local-factor evaluation performs divisibility tests"; no lower bound.
  - #13: the identity term −2J_w(1/√N) is included, and cusp terms can be √N-sized for other weights (w = U₁₀).
  - #14: de Haan–Jacobson–Williams is a real D^{1/6} regulator computation under GRH (hybrid, relational).
  - #15: Hiary 2014 (q = p² gives q^{1/2}; factoring allowed within its cost), Vishe's hypotheses.
  - #16: the N ≡ 1 (mod 4) hypothesis (N = 15 example) and Prop HR's qualifications (N = 731) are restored.
  - #17: Lagarias–Odlyzko space, Srinivasan real-quadratic.
  - #18: the verdict is now "audit outcome + restricted barrier", and the AFE "cliff" is described as rapid tail decay.
- **Opus literature critic (read-only): OK WITH FIXES.** Every paper exists. Fixed:
  - Srinivasan: real quadratic and ineffective.
  - de Haan–Jacobson–Williams: a computation.
  - Kronecker–Gierster–Hurwitz attribution; Popa–Zagier is a refinement.
  - "Rigorous L[2/5]" is internal and unrefereed.
  - Hiary 2014 ("factorisation not essential"); Hiary 2018 verified; Vishe's hypotheses; BHK caveats.
  - Sutherland §11.4 verified.
  - Gekeler is for prime determinant; Tao–Croot–Helfgott only sketch x^{1/3}.
  - Biasse–Erukulangara and de Boer–Pellet-Mary–Wesolowski (arXiv 2512.01588).
  - PRIOR-ART: §3.2 (Zagier 1977, Serre 1997, Conrey–Duke–Farmer, Altuğ); the Prop AM density (Siegel; KNOWN); Prop HW
    (AAK, Beylkin–Monzón, Plonka–Pototskaia); the Δ^{1/3} argument (Sutherland Thm 2.3).
  - Missing references added.
  - For Theorem LIN the critic searched and found no prior K-term lower bound, so "no prior statement found by us"
    stands.

## 10. Reproduction
```bash
PYTHONPATH=. .venv/bin/python -m pytest -q tests/test_nativeana.py
PYTHONPATH=. .venv/bin/python -m src.nativeana.harness E1   # ~1 min
PYTHONPATH=. .venv/bin/python -m src.nativeana.harness E2   # ~5 s
PYTHONPATH=. .venv/bin/python -m src.nativeana.harness E3   # ~1 min
```
