# analytic: non-relational analytic factoring (phase 17)

Agent `analytic`, 2026-10-20.

**Files.**
- Pre-registration: `notes/analytic_prereg.md`. Commit `7d18952` was made before E1/E2. Amendment A1
  (E3) is in commit `7286ac3`, made before E3 ran.
- Code: `src/analytic/{hiary,hull,ladder,cubic,run_experiments}.py`.
- Tests: `tests/test_analytic.py` (25 pass).
- Raw data: `experiments/raw/analytic_e{1,2,3}.jsonl`.
- Summaries: `experiments/analysis/analytic_summary.json`, `analytic_e3_summary.json`.
- Candidate cards: `notes/analytic_cards.md`.

**Reviews.**
- GPT-6.1-sol blind round: `codex/analytic_1.{md,out}`. The full report is at
  `codex/scratch/analytic_1/analytic_1_report.md` (git-ignored).
- GPT-6.1-sol adversarial round: `codex/analytic_2.{md,out}`. The full report is at
  `codex/scratch/analytic_2/analytic_2_report.md`.
- Two read-only Opus critics: one on correctness and cost, one on literature.
- Every fix from all four is applied in this version; §10 lists them.

**Labels.**
- **KNOWN**: with a reference.
- **PROVED**: the proof is given here. The proofs are elementary. "No prior formulation found by us" is
  the strongest priority claim made.
- **HEURISTIC**, **CONJECTURED**, **experimentally supported**, **REFUTED for this model**, **OPEN**.

## 0. Verdict (plain language)

**No analytic route to factoring was found that is subexponential. None even beats the deterministic
record N^{1/5}.** "Analytic" covers exponential sums, theta functions, L-functions and lattice-point
counting. L[1/3] is far out of reach. This is a verdict about known techniques and the mechanisms we
could build. It is not a lower bound.

What we established:

1. **Hiary's N^{1/3} factoring algorithm contains no exponential sum.**
   - It is a continued-fraction method. For each block [x − H, x + H], one Dirichlet approximation of
     N/x² turns "x + h divides N" into a quadratic equation.
   - Hiary's fast quadratic and cubic exponential-sum algorithms (Annals 2011) belong to his separate work
     on computing ζ. The task's premise conflated the two.
   - The exponent 1/3 comes from a denominator-amplified Taylor error: the block width is H ≍ x·N^{−1/3}.
   - Two independent facts give the same order:
     - the integer hull of {x, y > 0, xy ≥ N} has Θ(N^{1/3} log N) vertices, and every divisor is a
       vertex (KNOWN);
     - Hiary needs "(1/3)·log₂N bits", i.e. 2/3 of the bits of p, the 2-D lattice (Rivest–Shamir) strength. *(Corrected by the orchestrator after topblock: the earlier text said "1/3 of the bits of p".)*
   - These are consistent readings, not one theorem. A claimed universal "curvature law" explaining
     Lehman and Rubinstein as well was rejected on review.
2. **Reduction L_j, the floor-sum ladder (PROVED for each fixed j).**
   - Suppose exact floor sums Σ⌊P(h)⌋ of rational degree-j polynomials were available in polylog time.
     Then all divisors of N up to √N are found with O(N^{1/(j+2)} log² N) oracle calls.
   - j = 1 gives N^{1/3}, operational today with Euclid floor sums. j = 2 gives N^{1/4}. j = 3 gives
     N^{1/5}, the record. j = 4 gives N^{1/6}.
   - This improves the N^{1/(j+1)} of the blind codex construction. The trick is to locate the O(1)
     ambiguous points per block instead of demanding precision 1/x.
   - nonsmooth's N^{1/(2j+2)} has no valid derivation for this ladder.
   - E3 confirms correctness (170/170) and O(1) candidates and oracle calls per block.
   - The missing primitive is exact lattice-point counting under short, nearly parabolic polynomial arcs.
     It is **OPEN**: neither an algorithm nor a hardness result is known. Fourier methods cannot supply it
     (barrier F).
   - Uniform degree would give L[½] only with two further uniformity results that we do not have.
3. **Complete cubic Weyl sums mod N factor N (reduction C3, PROVED; E2 100/100).**
   - The sums are S(a, N) = Σ_{x mod N} e(ax³/N). Nine approximate values factor every N = pq with
     p ≡ q ≡ 1 (mod 3), through Jacobi sums in Z[ω] (randomised Turing reduction).
   - The precision needed is constant absolute error (rounding decoder; measured typical δ* ≈ 3). A
     conic small-root decoder tolerates N^{1/6−ε}; this is codex's argument and is not implemented.
   - The sign test "S(1, N) ≠ 0" decides the fixed-base QR bit "is −3 a square mod N" with absolute
     precision 1/18.
   - Consequence: a total-polylog cubic evaluator for the complete-sum coefficient range would factor
     these N. Hiary's cubic theorem covers only cubic coefficients ≤ K^{−2}.
   - The degree-2/degree-3 line for *pure complete* sums (Jacobi symbol versus factoring-hard) matches the
     known multivariate line of Cai–Chen–Lipton–Lu. A conjectured extension to every degree k is false
     (counterexample k = 4, N = 77).
4. **Dequantisation: no classical route.**
   - The Shor phase h·a^x/N is closed under Weyl differencing (KNOWN: Korobov's method). So differencing
     does not lower its degree, and no useful short dual transform is known.
   - Its bit-cut Schmidt rank equals the (odd part of the) order (KNOWN in substance:
     Dang–Hill–Hollenberg 2019). So tensor contraction pays ≍ r.
   - The QFT half is classically simulable only under input-access hypotheses. The significance of
     modular exponentiation is Yoran–Short's point.
5. **Barrier map (§6).**
   - Hull-walk output size.
   - A PROVED uncertainty principle: a band-limited divisor detector needs ≥ 3d/4 Fourier terms.
   - Ramanujan truncation, in the frozen-weight model.
   - Cubic hardness (C3).
   - A HEURISTIC cost balance for the approximation-access model.
   - None of these is a general factoring lower bound.

**Best idea.** The open question is item 2's primitive: count lattice points exactly under a short
degree-j arc with tiny higher coefficients, in time H^{o(1)}.
- Polylog at j = 2 only ties the old N^{1/4} (Pollard–Strassen, Coppersmith).
- Beating N^{1/5} needs j ≥ 4.
- L[½] additionally needs degree-uniform oracle costs and degree-uniform candidate bounds.
- We know no mechanism for any of these.

## 1. Hiary's algorithms, dissected

### 1.1 The factoring algorithm (KNOWN)

Reference: G. A. Hiary, *A deterministic algorithm for integer factorization*, Math. Comp. 85 (2016),
no. 300, 2065–2069, DOI 10.1090/mcom3037 (arXiv:1408.2608, read in full).

Algorithm 1:
- Trial-divide up to x₀ = ⌈(17N)^{1/3}⌉.
- Then, for blocks [x − H, x + H] with H = ⌊(17N)^{−1/3}x⌋:
  - take the convergent b/q of N/x² with the largest q ≤ 4H (Lemma 3.1);
  - set a = [qN/x];
  - solve (qN − ax) + (bx − a)h + bh² = 0 (Lemma 3.2);
  - test the ≤ 2 integer roots by division.
- Theorem 2.1: O(N^{1/3} log² N) operations on O(log N)-bit numbers.

**What reveals the factor.** Nothing analytic: no sum, no bisection, no precision issue. The exact
identity

  qN/(x+h) = a − bh + E(h),  E(h) = ε₁ − ε₂h/q′ + qNh²/(x²(x+h)),  |ε₁| ≤ ½,  |ε₂| < 1,

makes E(h) an integer at a divisor. The block condition forces |E(h)| < 1, hence E(h) = 0.

**Why 1/3 (PROVED; this is Hiary's estimate (4)).**
- The curvature term needs q·NH²/x³ < ¼, with Dirichlet denominator q generically ≍ H. Hence
  N(H/x)³ ≲ 1/16, i.e. H ≍ xN^{−1/3}.
- Blocks per dyadic range: (17N)^{1/3}(ln 2)/2 = 0.891·N^{1/3}.
- In total: ½(17N)^{1/3} ln(√N/(17N)^{1/3}) = 0.2143·N^{1/3} ln N − 1.21·N^{1/3}.
- **The binding resource is this denominator-amplified linearisation error.** It is not sum length,
  phase degree, precision, or cubic sums.
- Small-denominator ("major arc") approximations allow bigger blocks. Hiary's implementation gains
  about 20% from them.

**Consistent readings of the exponent** (not one theorem):
- **(a) 2-D Diophantine strength.** Hiary notes that the method needs the top (1/3)log₂N bits of p, the
  Rivest–Shamir (EUROCRYPT '85, LNCS 219, 31–34) 2-D strength. Coppersmith needs 1/4.
- **(b) Integer hull (KNOWN).**
  - Take H_N = {x, y > 0 : xy ≥ N}. Its integer hull has Θ(N^{1/3} log N) vertices: Balog–Bárány,
    arXiv 2602.06897, Thm 1.1 (2026 preprint).
  - Earlier, Alcántara–Blanco–Criado–Santos (arXiv 2501.19193) gave the bounds 2N^{1/3} − 2 and
    O(N^{1/3} log N), and an O(log)-per-vertex enumerator.
  - Every divisor is a vertex (KNOWN, Alcántara et al., by strict convexity; equivalently, by AM–GM,
    the functional (N/d)x + dy is minimised only at (d, N/d)).
  - Alcántara et al. already note that hull enumeration factors N deterministically in
    O(N^{1/3} log² N) operations.
  - Any method that visits every vertex in [N^{1/3}, √N] pays Ω(N^{1/3} log N): an output-size bound in
    the unit-cost model. It does not bound early-stopping or non-enumerating methods.
  - E1 measured, exactly, 0.480–0.490·N^{1/3} vertices per full dyadic window, with max/min ≤ 1.03, for
    24–40-bit N. This is numerical evidence only. The padded-window routine is not certified to be
    complete, but all divisors were found as vertices.
  - Hiary blocks / hull vertices was 1.06–1.38 over the measured ranges; per full window it tends to
    0.891/0.485 ≈ 1.84. Hiary's blocks are not hull edges: only the order of magnitude agrees.
- **(c) Affine arc length (HEURISTIC comparison).**
  - The hyperbola's affine arc-length element is κ^{1/3}ds = (2N)^{1/3}dx/x, which gives N^{1/3} per
    logarithmic window.
  - This is a heuristic for local-linearisation methods only. It is not a deterministic local law, since
    thin full-curvature arcs can have O(1) hull vertices (§3).
  - Lehman (Farey multipliers k ≤ N^{1/3}) and Rubinstein (displacement against occupancy, Kloosterman
    second moments) are **not** hull walks. That four methods share the exponent 1/3 is not explained
    by one mechanism (codex analytic_2, C1).
  - For x² + y² = N, the disc hull also has ≍ N^{1/3} vertices (Balog–Bárány, DIMACS Ser. 6 (1991)
    39–44). Most of them lie inside the circle; only the r₂(N) points on it are factor-revealing.

### 1.2 Hiary's exponential-sum algorithms (KNOWN; they compute ζ, not factors)

- **Quadratic.** Hiary, Ann. of Math. 174 (2011) 859–889. Σ_{k≤K} k^r e(ak + bk²) is computed to error ε
  in poly(log(K/ε)) operations, with no precomputation. The method is repeated theta (Poisson / Weil)
  transformations: a continued-fraction recursion in (a, b).
- **ζ(½ + it).** Hiary, Ann. of Math. 174 (2011) 891–946.
  - Quadratic blocks of the Riemann–Siegel sum give t^{1/3}.
  - The main result, t^{4/13}, rests on Theorem 2.2: Σ k^r e(ak + bk² + ck³) with 0 ≤ c ≤ K^{μ−3},
    μ ≤ 1, at polylog cost per query after a K^{4μ} FFT precomputation.
  - Balancing the block queries t^{1/2−β} against K^{4μ} gives β = 5/26. The binding limit is the size
    of the cubic coefficient.
- **Analogy (HEURISTIC).**
  - Riemann–Siegel ↔ trial division.
  - Odlyzko–Schönhage amortisation ↔ Pollard–Strassen / Bostan–Gaudry–Schost multipoint evaluation.
  - Theta-blocking at t^{1/3} ↔ CF-blocking at N^{1/3}.
  - It breaks at exactness. ζ needs an O(1)-error value; factoring needs an exact count, which
    band-limited numerics cannot give (barrier F).

## 2. Candidate table

The cost column is the best known classical cost without factoring first. "Minimal requirement" is what
an evaluator must deliver for the factor to be read off.

| Candidate quantity | Link to a factor | Minimal requirement | Best known evaluation | Below N^{1/3}? N^{1/5}? N^ε? |
|---|---|---|---|---|
| Restricted divisor count C_N(X) = D_N(X) − D_{N−1}(X), with D_N(X) = Σ_{d≤X}⌊N/d⌋ | C_N(X) = #{d ≤ X : d \| N}; bisection, O(log N) queries (folklore) | error < ¼ in each D | N^{1/3}: Hiary; hull walk; ladder j = 1 | no / no / no |
| Full D(N) = Σ_{m≤N} τ(m) | first difference = τ(N): useless | — | N^{1/3}: Tao–Croot–Helfgott, Math. Comp. 81 (2012), Prop. 2.2; Sladkey arXiv 1206.3369. Li, arXiv 2609.29949 (preprint), claims N^{9/28}, using Lenstra–Pomerance factoring of all integers within N^{9/28} of N | not applicable |
| Sawtooth Σ_{d∈I} ψ(N/d), Σ_d e(hN/d) | jump at divisors; the Fourier series loses the endpoint term | all frequencies up to ≍ x | one frequency: Õ((\|h\|N)^{1/4}) by theta blocking (codex); all needed frequencies: ≥ 3d/4 terms (F) | no |
| Voronoi dual sums with smoothing width H | dual length N/H²; de-smoothing re-exposes an H-window (codex) | exact | no restricted-cutoff version below N^{1/3} known | no |
| r₂(N) via θ² | r₂(pq) ∈ {0, 16}: one bit. Two inequivalent representations factor N | the points on x² + y² = N | disc hull walk: N^{1/3} | no |
| r₄(N) = 8σ(N), N odd | s = r₄/8 − N − 1 = p + q (KNOWN) | error < 4 (additive N^{1/4} gives a Coppersmith hint) | extraction needs ≍ N samples; the singular series hits R and also omits cusp parts for higher k | no |
| σ(N) via Ramanujan sums | corrected terms below modulus p vanish (nonsmooth §4.1; R) | moduli ≥ p | √N | no |
| h(−4N) via L(1, χ) | class group + ambiguous forms (Schoof 1982) | \|δL\| < π/(2√Q) | AFE √Q; Shanks' Euler-product + BSGS: Q^{1/5} under ERH (Shanks; Sutherland thesis §11.4) | ties 1/5 under ERH |
| Complete cubic Weyl sum S(a, N) | 9 values → p (C3) | constant absolute error (N^{1/6−ε} with a conic decoder) | O(N); O(p + q) given the factors | no |
| Legendre-symbol interval Σ_{a<N/2}(a/N) | h(−N) (class b) | exact | √N | no |
| Shor sums Σ_x e(h a^x/N) | the order of a ⇒ factor (Miller, randomised) | — | generic ≳ √r; tensor contraction ≍ r (§5) | no |

## 3. Reduction L_j: the floor-sum ladder

**Theorem L_j (PROVED, conditional on an oracle; fixed j).**
- Fix j ≥ 1 and 0 < c ≤ ½. Let B = ⌈N^{1/(j+2)}⌉. The implementation uses ⌊·⌋, which only changes
  constants.
- Suppose an oracle returns F(P; h₁, h₂) = Σ_{h₁≤h<h₂}⌊P(h)⌋ exactly. P ranges over rational
  polynomials of degree ≤ j with O(j log N)-bit coefficients, including rational constant shifts and
  arbitrary subintervals.
- Then the algorithm below finds every divisor of N in [2, √N], or certifies that there is none, using
  trial division below ⌈B/c⌉ plus:
  - O_{j,c}(B log N) blocks;
  - O_{j,c}(B log N) candidate divisions;
  - O_{j,c}(B log² N) oracle calls.

*Proof.*
1. **Blocks.** Use blocks [x, x + H) with H = ⌊cx/B⌋ ≥ 1; there are O(B) per dyadic range.
2. **Taylor polynomial.** Write f(h) = N/(x+h) = P(h) + R(h), with P(h) = Σ_{i≤j}(−1)^i Nh^i/x^{i+1}
   (denominator x^{j+1}) and R(h) = N(−h)^{j+1}/(x^{j+1}(x+h)).
3. **Remainder bound.** For 0 ≤ h < H, R has sign (−1)^{j+1}, and
   |R| ≤ ρ := NH^{j+1}/x^{j+2} ≤ c^{j+1}B/x ≤ c^{j+2} < 1. This uses B^{j+2} ≥ N and x ≥ B/c.
4. **Divisors are candidates.** Let h > 0 with x + h | N.
   - j odd: {P(h)} = 1 − R(h) ∈ [1 − ρ, 1).
   - j even: {P(h)} = |R(h)| ∈ [0, ρ].
   - At h = 0, R = 0, so the block start is tested separately (one division per block).
5. **Counting candidates.** Use ⌊P + ρ⌋ − ⌊P⌋ for odd j. For even j use ⌊P⌋ − ⌊P − ρ − η⌋, with
   0 < η < 1/(DV), where P ∈ D^{−1}Z and ρ = U/V; this gives the closed endpoint. The code takes
   η = 1/(4DV).
6. **Locating candidates.** Bisection costs O((1 + |C_block|) log H) calls per block.
7. **Size of C.** For a near-integer candidate put m = ⌈P(h)⌉ (odd j) or ⌊P(h)⌋ (even j). Then
   |N/(x+h) − m| ≤ ρ, so d = x + h divides N − k with |k| ≤ ρd ≤ c^{j+1}B(1 + c/B) ≤ B/2.
   - Hence the total number of candidates is ≤ Σ_{|k|≤B/2} τ(N − k) ≪_j B log N.
   - This uses Shiu, J. reine angew. Math. 313 (1980) 161–170, Thm 1. Its hypotheses hold for τ, in
     the range y = N^{1/(j+2)} ≥ N^ε, and the constant depends on j.
   - No equidistribution is used. The forced h = 0 tests (one per block) are counted separately. ∎

**Relation to the other two exponents.**
- Codex's blind construction (analytic_1 §1.4) demands precision 1/(8M), so that P reproduces every
  floor. That forces H ≍ xN^{−1/(j+1)} and gives N^{1/(j+1)}. It is a valid but weaker upper bound.
  Precision ρ ≲ 1/H plus candidate location gains one degree.
- At j = 1 the ladder is an N^{1/3} algorithm built only from Euclid floor sums. That agrees with Hiary
  and with the divisor-sum blocking of Tao–Croot–Helfgott (2012) Prop. 2.2; precision 1/x would give
  N^{1/2} there.
- nonsmooth's N^{1/(2j+2)} assumes that remainder < 1 suffices. Then ≍ ρH ≍ H points per block are
  ambiguous, and the total cost is about √N. **That derivation is invalid for this ladder.** This is an
  analytic argument; E3 did not test it. We do not prove that no other reduction could reach it.
- Related prior work: Khomovsky (arXiv 1802.03658) uses degree-1 and degree-2 integer Taylor
  approximations of N/x to search for divisors. Tao–Croot–Helfgott §4 remark that piecewise-polynomial
  approximations "may ultimately be more efficient". Neither states the 1/(j+2) law: no prior formulation
  found by us.

**E3 (amendment A1).** The j = 1 oracle is operational (Euclid). For j = 2, 3 the oracle is simulated by
brute force: only call counts are model costs.

| j | bits | correct | located candidates / block | c^{j+2} | oracle calls / block | slope of ln(blocks/ln N) (forced by construction) |
|---|---|---|---|---|---|---|
| 1 | 24–48 | 70/70 | 0.064 → 0.093 | 0.125 | 2.19–2.68 | 0.336 |
| 2 | 24–40 | 50/50 | 0.038–0.053 | 0.0625 | 2.22–2.46 | 0.252 |
| 3 | 24–40 | 50/50 | 0.013 → 0.034 | 0.031 | 2.10–2.40 | 0.201 |

**What E3 tests.**
- Block counts are fixed by H = ⌊cx/B⌋. So P3.2 and P3.3 ("slope ≈ 1/(j+2)") are consistency checks
  that could only fail through a bug. P3.3's decision rule is tautological.
- The content is:
  - correctness, 170/170, with every divisor caught as a candidate;
  - located candidates per block stay O(1), rising toward the heuristic c^{j+2}·N/B^{j+2};
  - 2–3 oracle calls per block.
- So the total model cost tracks the block count. E3 says nothing about the cost of a j ≥ 2 oracle, or
  about competing constructions.

**Uniform degree and L[½] (conditional; not proved).**
- *If* the whole degree-j reduction cost exp(log N/(j+2) + O(j log j))·polylog N, then
  j ≍ √(log N/log log N) gives L_N[½, O(1)]. This is the same optimisation as codex analytic_1 §4.1.
- It needs **two** results we do not have:
  - (i) a uniform exact floor-sum oracle with T_j = e^{O(j log j)};
  - (ii) a candidate bound uniform in j.
- Shiu's constant depends on ε = 1/(j+2) → 0. The unrestricted majorant Σ_{|k|≤Y} τ(N − k) is genuinely
  false at that scale: for N = M + 1 with M primorial, τ(M) = exp(Θ(log M/log log M)) exceeds every
  L[½, C] (codex analytic_2, C7). A semiprime-promise or algorithm-specific bound might survive; none is
  known.

**The missing primitive (OPEN).**
- At j = 2 near x ≍ √N, the block has H ≍ N^{1/4}, slope ≍ 1, and quadratic coefficient β with
  βH² ≍ √N/x (≍ 1 at x ≍ √N, larger for small x).
- For such full-curvature short arcs, a hull walk *typically* visits ≍ H^{1/3} vertices, giving
  N^{1/4}·N^{1/12} = N^{1/3}: no gain. This is not a universal lower bound.
  - Codex found exact families with O(1) vertices: P = h²/H², and the Taylor family x = H², N = H⁴ + 1,
    where ⌊P₂(h)⌋ = x − h exactly.
- No polylog exact evaluator and no hardness result are known.
- Euclid solves j = 1. One continued-fraction step maps a parabola count to a conic-arc count (square
  roots appear). A fully general polylog conic-sector counter would locate the points of x² + y² = N by
  bisection and so factor N when p ≡ q ≡ 1 (mod 4). The short-regime question is narrower than that.
- Complete quadratic floor sums lie in a different regime (Q below): there,
  Σ_{x<m}⌊ax²/m⌋ = a(m−1)(2m−1)/6 − (m−1)/2 + (a/m)h(−m).
- Barrier F excludes exponential-sum evaluation as a substitute.

**Approximation-access cost balance (HEURISTIC; not a lower bound).**
- The Chebyshev bound E_j(f; [a, b]) ≥ (b − a)^{j+1}min|f^{(j+1)}|/(2^{2j+1}(j+1)!) is correct.
- Suppose further that ambiguities cost ≳ H·E_j per block. Then balancing blocks against ambiguities
  gives ≍_j N^{1/(j+2)} per dyadic range near √N.
- The step from E_j to ambiguous *integer* points is unproved, and the full-curvature families above
  violate it.

## 4. Reduction C3: complete cubic Weyl sums mod N

**Setting.** S(a, N) = Σ_{x mod N} e(ax³/N), with N = pq, p ≠ q primes > 3 and gcd(a, N) = 1.

**Lemma 1 (KNOWN).**
- (i) **CRT.** S(a, pq) = S(aq̄, p)·S(ap̄, q), where q̄q ≡ 1 (mod p) and p̄p ≡ 1 (mod q). This holds
  because q̄q + p̄p ≡ 1 (mod pq).
- (ii) **p ≡ 2 (mod 3).** S(b, p) = 0, since cubing permutes F_p.
- (iii) **p ≡ 1 (mod 3).** S(b, p) = χ̄(b)G(χ) + χ(b)\overline{G(χ)} = 2√p cos(θ_p − 2πk(b)/3), where
  G(χ) = √p e^{iθ_p}, χ(b) = e(k(b)/3), and χ(−1) = 1. The three values are Gauss's cubic periods: the
  roots of s³ − 3ps − pL = 0, with 4p = L² + 27M².
- (iv) **Jacobi sums.** G(χ)³ = p·J(χ, χ), and J(χ, χ) = π_p ∈ Z[ω] has norm p (Ireland–Rosen, Ch. 8 §3
  and Ch. 9 §4). So π_p/√p = e^{3iθ_p}.

**Lemma 2 (PROVED, via Gauss's periods; codex analytic_2).** |S(b, p)| > 1/3 for every unit b. Here L ≠ 0,
because a trace-zero element of Z[ω] has norm 3t² ≠ p. If |s| ≤ 1/3 then |s³ − 3ps| < p ≤ |pL|, a
contradiction.

**Theorem C3 (PROVED; randomised Turing reduction).** Let N = pq, p ≡ q ≡ 1 (mod 3), and let a, b be
random units. Then the nine values S(a^i b^j, N), 0 ≤ i, j ≤ 2, determine p in poly(log N) exact
operations, with probability 48/81 per pair.
- The values may carry an absolute error ≤ c₀ for an absolute constant c₀ > 0. The first-order estimate
  gives coordinate error ≈ 2.7δ, so the ±1 search tolerates δ ≈ ½ (critic).
- Repeating until division verifies the factor gives Las Vegas expected ≈ 81/48·9 ≈ 15 oracle calls.

*Proof.*
1. **The array.** S(a^i b^j, N) = 4√N u_{α(i,j)} v_{β(i,j)}, where u_k = cos(θ_p − 2π(k + s_p)/3) and
   α, β : F₃² → F₃ are the cubic-class functionals at p and q.
2. **Relabeling.** With probability |GL₂(F₃)|/81 = 48/81 the pair (α, β) is invertible. The decoder tries
   all 48 ordered pairs of independent functionals.
3. **A direction.** Take u from the column with the largest |v_l|; a cosine 3-vector has
   max|v_l| ≥ √3/2. Then z_u = Σ_k u_k ω^k = (3/2)e^{iθ_p}ω^{−s_p}·(a real scalar). Hence z_u³ is a real
   multiple of π_p, or of π̄_p after a reflection relabeling.
4. **Two norm-N elements.** γ = z_u³z_v³ and γ′ = z_u³\overline{z_v³} are real multiples of π_pπ_q and
   π_pπ̄_q (up to conjugation and units). Rescaling to norm N and rounding in Z[ω] recovers both.
5. **The factor.** gcd(γ, γ′) = π_p up to a unit, since q splits into non-associate primes. Its norm
   is p.
6. **No false positives.** Every output is verified by division. ∎

**Corollaries.**
- **(i) Fixed-base QR bit (PROVED).** For N ≡ 1 (mod 3), with distinct primes, S(1, N) ≠ 0 iff
  p ≡ q ≡ 1 (mod 3), i.e. iff −3 is a square mod N.
  - Since then |S(1, N)| > 1/9, absolute precision 1/18 decides the bit.
  - At prime squares this fails: S(1, ℓ²) = ℓ for ℓ ≡ 2 (mod 3).
  - This is the nonsmooth §5.5 / Track E minimal-primitive bit. It is not a reduction of factoring to
    that bit.
- **(ii) Barrier for a "cubic Hiary" (PROVED, with scope).**
  - Hiary's cubic Theorem 2.2 needs c ≤ K^{−2} (μ ≤ 1) and allows a K^{4μ} precomputation.
  - The complete sum has K = N and c = a/N mod 1: generic c ≍ 1 (μ = 3) for random a, and μ = 2 only for
    a = O(1).
  - A **total-polylog** evaluator (setup included) to absolute error c₀ factors these N, in either of
    two settings:
    - (a) it works for generic coefficients;
    - (b) it works for pure coefficients 1/M at *all* moduli M, via codex's padding identity
      S(1, a²N) = a·S(a, N) for squarefree a with gcd(a, 3N) = 1. (Alternatively: small a, b with
      independent cubic characters, which exist under GRH.)
  - A huge-precomputation extension is not excluded.
- **(iii) Degree 2 versus 3 for pure complete sums.**
  - Pure quadratic complete sums are ε_N(a/N)√N: Jacobi symbols.
  - Cubic ones factor split semiprimes.
  - The multivariate analogue of this line is KNOWN: Cai–Chen–Lipton–Lu, *On tractable exponential sums*
    (FAW 2010, arXiv 1005.2632). Quadratic sums mod any N are polynomial time; some cubic families are
    #P-hard.
  - **REFUTED:** our conjectured extension "every degree k with gcd(k, p−1), gcd(k, q−1) > 1". For k = 4
    and N = 77, S₄ = S₂ identically.
  - It is also false that every quadratic arithmetic sum is easy: see Q.
- **Prior art.** van Dam–Seroussi (quant-ph/0207131) reduce discrete log to estimating Gauss sums of
  multiplicative characters. That is a different reduction. The ingredients of C3 are classical: Gauss
  periods, Jacobi sums, and the Euler-type gcd of two norm-N representations. For the univariate
  composite-modulus statement: no prior formulation found by us.

**Precision.**
- The rounding decoder needs constant absolute error.
- **This is not a necessity statement.**
  - The approximate γ gives (A₀, B₀) within O(E) of the true point. The true point solves the
    total-degree-2 equation (A₀+u)² − (A₀+u)(B₀+v) + (B₀+v)² = N.
  - Coron's small-root theorem (EUROCRYPT 2004, LNCS 3027, Thm 5) recovers it for E < N^{1/6−ε}
    (codex analytic_2).
  - This is PROVED conditional on that theorem, and not implemented.
  - Gutiérrez–Ibeas–Joux (J. Symb. Comput. 64 (2014)) reach N^{1/4} for the sum-of-two-squares analogue.
- So about (1/6)log₂N of the ≈ ½log₂N significant bits of S can be dropped: the decoder needs only
  ≈ (1/3)log₂N bits. How many the invariant really needs is open.

**E2.**

| bits | p ≡ q ≡ 1: success | structural | accidental | control max \|S\| | control structural | δ* (typical, ≥ 50% at ±δ) |
|---|---|---|---|---|---|---|
| 14 | 20/20 | 15 | 5 | 6·10⁻¹³ | 0 | 10 |
| 16–22 | 80/80 | 80 | 0 | ≤ 2·10⁻¹⁰ | 0 | 3 |

- P2.1 and P2.2 hold.
- **P2.3 missed as specified.** δ* = 10 at 14 bits lies outside [0.03, 3]. The slope −0.170 lies outside
  [−0.15, 0.15]; it comes from one coarse-grid step.
- δ* is a *typical* tolerance of the rounding decoder. The critic measured 60/60 successes up to
  δ = ½ under corner noise.
- The 5/20 accidental gcd hits at 14 bits are a ≈ 1% tail event (expected ≈ 1.2). The structural rate
  per attempt is consistent with 48/81 (16/32 and 0.575 in two re-checks).
- A cosmetic bug: the "all S = 0" early exit compared floats with exact 0. The controls therefore ran all
  6 attempts, which inflated their accidental hits. It does not affect any structural conclusion.

## 5. Dequantisation verdict (classical only)

**The object.** The joint amplitude of observing (k, y) is
A(k, y) = M^{−1}Σ_{x<M, a^x≡y} e(kx/M) = (MN)^{−1}Σ_{h mod N}Σ_{x<M} e(kx/M + h(a^x − y)/N). Conditioning
on y renormalises by √(M·L_y).

- **QFT simulability is conditional (KNOWN).**
  - Griffiths–Niu (PRL 76 (1996) 3228) replace the measured QFT by adaptive single-qubit measurements.
  - The classical simulations need input-access hypotheses: Aharonov–Landau–Makowsky (quant-ph/0611156);
    Browne (NJP 9 (2007) 146); Yoran–Short (PRA 76 (2007) 042321).
  - Yoran–Short (PRA 76 (2007) 060302(R), arXiv 0706.0872) show that simulating the modular-exponentiation
    circuit on the specified inputs, with product measurements, would simulate Shor.
  - The unresolved resource is therefore the coherent correlations that modular exponentiation creates.
    Evaluating a^x mod N is itself easy.
- **Differencing (KNOWN).**
  - Δ_u(αx + h a^x/N) = αu + h(a^u − 1)a^x/N. The exponential part stays in the family. The basis is
    Korobov's method for Σ e(a bⁿ/m); see Vandehey, arXiv 1606.07911.
  - Iterating can only annihilate the coefficient at a period, or expose a gcd. Both need the answer.
  - A stationary-phase dual exists on a real lift. Its dual range is huge, and the lift of modular data is
    the issue. "No useful short transform is known" is the defensible statement, not nonexistence.
- **Bit-cut rank (KNOWN in substance; Dang–Hill–Hollenberg, Quantum 3 (2019) 116).**
  - Write x = 2^s x_hi + x_lo and g(t) = e(h a^t/N) on Z/r. Then the cut matrix is Ψ_hi·diag(ĝ)·Ψ_loᵀ,
    and its rank is ≤ min(2^s, 2^{m−s}, r/gcd(r, 2^s)).
  - Our first statement omitted the gcd. Codex's counterexample: N = 15, a = 2, s = 4 gives rank 1.
  - For odd r, consecutive exponents make both factors Vandermonde matrices with distinct nodes. Since
    2^m ≥ r² ≥ |supp ĝ|², one factor has full column rank on the support. So the rank is exactly
    min(#distinct rows, #distinct columns, |supp ĝ|) (critic: 0 counterexamples in 43,944 cases).
  - That ĝ(ℓ) ≠ 0 for most ℓ (Gauss-type periods at composite N) is HEURISTIC. For prime N, Chebotarëv's
    nonvanishing-minors theorem (Tao, MRL 12 (2005)) gives maximal rank.
  - Exact rank lower-bounds exact MPS bond dimension along this order, ≍ r. It does not lower-bound
    approximate simulation or other contraction orders.
- **Query bounds (KNOWN, black-box only).**
  - Cleve, Inf. Comput. 192 (2004) 162–171: Ω(r^{1/3}/√log).
  - Chakraborty–García-Soriano–Matsliah, ICS 2011, 284–297: Ω̃(√r).
- **Verdict.** All known fast-summation techniques need a low polynomial degree or a low rank: Weyl/van
  der Corput, Poisson/theta, Taylor blocking, Odlyzko–Schönhage, and tensor contraction. The Shor phase
  has neither: it is closed under shifts, and its rank equals the hidden order. **No classical piece of
  the period information is extractable by these techniques.** This is not a theorem against all
  classical algorithms.

## 6. Barrier map

| # | Statement | Scope | Status |
|---|---|---|---|
| H | Visiting every integer-hull vertex of H_N in [N^{1/3}, √N] costs Ω(N^{1/3} log N) | enumerating methods; not early-stopping or selective ones | KNOWN vertex count (Balog–Bárány 2026 preprint) + trivial output bound |
| F | **Uncertainty principle.** Let F(u) = Σ c_k e(ku), F(0) = 1, \|F(r/d)\| ≤ ε for all r ≢ 0 (mod d), for some d ∈ [x, 2x] with (d − 1)ε² < ¼. (This follows from ε < 1/(2W), W ≥ √(2x), the rounding condition for a window sum.) Then F occupies ≥ 3d/4 residue classes of frequencies mod d; for a band-limited F, K ≥ 3x/8 − 1. Proof: DFT of the samples δ₀ + e; at each empty residue ê = −1/d; Parseval gives (d − 2K − 1)/d² ≤ (d−1)ε²/d. Also stable: Σ\|c_h\| ≤ A forces HA ≳ d (codex) | uniform detectors evaluated term by term; not joint algorithms, single-N detectors, or exact-arithmetic methods. For W = Θ(x) it matches trial division in scale | PROVED (elementary; no prior statement with these constants found) |
| R | **Ramanujan truncation (frozen weights).** Σ_{q≤Q} a_q(X)c_q(n) = Σ_{q≤Q} a_q(X)μ(q) whenever P⁻(n) > Q. This applies to σ (nonsmooth §4.1). For r_k with every even k, odd-q singular-series terms are ε_q^k q^{−k/2}c_q(n). The singular series omits cusp contributions for larger k (Rankin 1965) | weights frozen across compared inputs | KNOWN / PROVED |
| C3 | A total-polylog, O(1)-accurate complete cubic Weyl-sum evaluator ⇒ factoring split semiprimes | generic coefficients, or all moduli via padding | PROVED; E2 |
| Q | Σ_{x<m}⌊ax²/m⌋ = a(m−1)(2m−1)/6 − (m−1)/2 + (a/m)h(−m), for prime m ≡ 3 (mod 4), m > 3 | complete periods | KNOWN (Dirichlet; Chamberland–Dilcher arXiv 2510.04387 eq. (5.10)); checked for 154 primes. Note: h(−m) is L[½] under GRH (Hafner–McCurley), so this is not exponential hardness |
| L | Approximation-access balance ≍ N^{1/(j+2)} | that model | HEURISTIC (cost assumption not proved; counterexample families exist) |
| P | Resolving {N/d} = 1/d needs Ω(log N) fractional bits; summing L independently rounded terms needs log₂L more. No universal per-term rule exists | float classifiers | PROVED (codex analytic_1 §5.3) |
| S | σ(n) is not P-recursive (every root of unity is a singularity of Σσ(n)zⁿ) | fixed-order holonomic recurrences | PROVED (codex analytic_1 §4.3; checked) |
| D | Shor phase: closed under differencing; bit-cut rank = order data | Weyl/vdC; exact MPS contraction | KNOWN (Korobov; Dang–Hill–Hollenberg) |
| I | Information: only floor F1 (functions of N mod m carry 0 bits). No information-theoretic sum-length bound exists, since N determines its factors | — | KNOWN (Track E) |

**What is not a barrier.**
- Phase degree by itself: theta transformations make every quadratic sum polylog, and Hiary handles small
  cubic coefficients.
- The obstructions are exactness (F), full-range cubic coefficients (C3), and the unknown exact-counting
  primitive (§3).

## 7. Subexponential analytic ideas: feasibility

| Idea | Calculation | Verdict |
|---|---|---|
| Uniform-degree ladder (§3) | L[½, O(1)] iff T_j = e^{O(j log j)} **and** a j-uniform candidate bound | neither is known; the unrestricted majorant is false at that scale |
| Many-variable theta (Siegel transformations) | enumeration (b + g)^{g/2} (Deconinck–Heil–Bobenko–van Hoeij–Schmies, Math. Comp. 73 (2004) 1417–1442); xy = N needs a quartic penalty | killed |
| Gaussian product isolation M₁/Z (codex analytic_1 §6B) | exact for λ ≈ 4 log N, but the phase is quartic; Fourier linearisation returns an e^{−itN} extraction integral | killed |
| Fermat / square-sieve character sums | detecting one square among T values needs ≳ T characters (second moment); T ≍ √N | exponential |
| Polynomial phases after a factor-encoding substitution | needs the substitution centre ≈ p to N^{1/4}: the Coppersmith hint itself | circular |
| L-values / class numbers | √Q directly; Q^{1/5} under ERH with BSGS | exponential |
| Holonomic / hypergeometric (BGS) | N^{1/4}; σ is not D-finite (S) | exponential |
| Voronoi (Li 9/28) | full D(N) only; correction factors all integers near N with Lenstra–Pomerance | not applicable |

No analytic sum was found whose fast evaluation is plausible and would give L[½] or better.

## 8. Best idea and the sharpest open question

**OPEN.** Decide exactly #{(h, m) : 0 ≤ h < H, m ≤ αh + βh² + γ}. The coefficients are rational, with
denominators H^{O(1)} and βH² ≍ 1. More generally, the arcs have degree j and tiny higher coefficients.
The target is time H^{o(1)}, with arbitrary shifts and subintervals.

Consequences of an answer:
- Polylog at j = 2: deterministic N^{1/4}, a new route to an old exponent.
- j = 3: N^{1/5}, tying the record.
- j ≥ 4: below the record.
- Uniform in j, together with a uniform candidate bound: L[½].

What is known:
- Euclid solves j = 1.
- Hull walking typically costs H^{1/3}, but some families are trivial.
- Complete-period quadratic sums encode class numbers.
- The full-range cubic exponential-sum analogue is factoring-hard (C3).
- Exponential sums cannot substitute (F).

We see no mechanism (HEURISTIC assessment).

## 9. Pre-registration outcomes

**E1.**
- P1.1 140/140. ✓
- P1.2 blocks/prediction ∈ [0.978, 0.992]. ✓
- **P1.3 MISS:** slope 0.3478 against [0.32, 0.345].
- P1.4 ✓: V/(N^{1/3} ln N) = 0.120–0.125, slope 0.331, window max/min ≤ 1.03.
- P1.5 25/25. ✓
- P1.6 ✓: V/blocks = 0.72–0.94.
- **The pre-registered decision rule is NOT MET** (it required P1.3).
- Post-hoc reading: the window ignored the finite-size term −1.21·N^{1/3}, which is 17% at 48 bits. The
  local slope at these sizes is 0.345–0.36 (critic). P1.2, which encodes the exact count, holds.

Hiary means for primes:

| bits | blocks | CF steps/block | blocks/(N^{1/3} ln N) |
|---|---|---|---|
| 24 | 533 | 3.80 | 0.138 |
| 32 | 4 820 | 4.23 | 0.157 |
| 40 | 41 057 | 4.65 | 0.169 |
| 48 | 328 960 | 5.06 | 0.176 |

**E2.**
- P2.1 ✓.
- P2.2 ✓.
- P2.3 miss (§4).

**E3.**
- P3.1 ✓, P3.4 ✓, P3.5 ✓.
- P3.2 and P3.3 hold, but by construction (§3).

**Q check.** All 154 primes 7 ≤ m < 2000 with m ≡ 3 (mod 4), at a = 1: 0 failures. Codex checked
3 coefficients each for 23 primes.

## 10. Reviews and the fixes applied

**GPT-6.1-sol blind round (analytic_1).** Independently found:
- that Hiary uses no sum;
- the zeta balances;
- C_N(X) bisection;
- barriers F (stable form) and P;
- the non-D-finiteness of σ;
- the Li 9/28 caveat;
- two mechanisms, both killed or left open.

Its N^{1/(j+1)} ladder is improved here to N^{1/(j+2)}.

**GPT-6.1-sol adversarial round (analytic_2): verdicts.**
- C1 REJECT as a unified explanation. Fixed: "consistent readings"; Lehman and Rubinstein are not hull
  walks; the disc caveat.
- C2 ACCEPT WITH FIXES. Fixed: ceiling root, h = 0, closed endpoint, the |N − md| ≤ ρd candidate bound,
  and the lower bound downgraded to a heuristic.
- C3 ACCEPT WITH FIXES. Fixed: the "≥ 3d/4 terms" version and the scope.
- C4 ACCEPT WITH FIXES. Fixed: frozen weights, cusp parts.
- C5 ACCEPT WITH FIXES. Fixed: the precision necessity is retracted (N^{1/6−ε} conic decoder); the μ
  range and the padding identity; the QR bit made robust (1/18); the degree-k conjecture REFUTED (k = 4,
  N = 77).
- C6 REJECT as stated. Fixed: the gcd in the rank formula and the counterexample; QFT access hypotheses;
  the B-process wording; amplitude normalisation.
- C7 REJECT as a proved feasibility statement. Fixed: the two uniformity requirements; the primorial
  obstruction; the H^{1/3} lower bound retracted (counterexamples); the primitive relabeled OPEN; the
  βH² scope.

**Opus correctness/cost critic: ACCEPT WITH FIXES.**
- The step-7 constant.
- L[½] needs a uniform short-interval bound.
- "Improves on" codex rather than "corrects".
- The nonsmooth argument is analytic, not E3.
- P3.3 is tautological.
- C3: the μ = 2 hypothesis, the √3/2 bound, δ* is typical.
- "Consistent" rather than "equivalent".
- The Vandermonde step for the rank.
- The E1 decision rule NOT MET.
- The E2 accidental hits are a tail event.
- Q is not exponential hardness.

**Opus literature critic: OK WITH FIXES.**
- Helfgott–Thompson replaced by Tao–Croot–Helfgott 2012.
- Li's de-smoothing described correctly.
- Venue added for Balog–Bárány 1991.
- Hull-walk cost O(N^{1/3} log² N).
- Relabeled KNOWN: divisor = vertex (Alcántara et al.); Q (Chamberland–Dilcher / Dirichlet); Shor
  differencing (Korobov / Vandehey); bit-cut rank (Dang–Hill–Hollenberg).
- Added: Cai–Chen–Lipton–Lu; Yoran–Short 060302(R); Khomovsky; Tao–Croot–Helfgott §4; Gauss's cubic
  periods; Ireland–Rosen chapters; Shanks credit.
- No fabricated references found.

**What survives.**
- Theorem L_j for fixed j (PROVED, oracle-conditional; no prior formulation found by us).
- Theorem C3 and the robust QR-bit corollary (PROVED; univariate composite-modulus form not found in the
  literature).
- Barrier F (PROVED, scoped).
- The Hiary dissection.
- The dequantisation verdict, as a statement about known techniques.

## 11. Reproduction

```bash
.venv/bin/python -m pytest -q tests/test_analytic.py
PYTHONPATH=. .venv/bin/python -m src.analytic.run_experiments e1        # ~6 min
PYTHONPATH=. .venv/bin/python -m src.analytic.run_experiments e2        # ~8 min
PYTHONPATH=. .venv/bin/python -m src.analytic.run_experiments e3        # ~10 min
PYTHONPATH=. .venv/bin/python -m src.analytic.run_experiments summary
PYTHONPATH=. .venv/bin/python -m src.analytic.run_experiments summary_e3
```
