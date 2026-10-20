# analytic: non-relational analytic factoring (phase 17)

Agent `analytic`, 2026-10-20.
- Pre-registration: `notes/analytic_prereg.md`. Commit `7d18952` was made before E1/E2. Amendment A1
  (E3) is in commit `7286ac3`, made before E3 ran.
- Code: `src/analytic/{hiary,hull,ladder,cubic,run_experiments}.py`.
- Tests: `tests/test_analytic.py` (25 pass).
- Raw data: `experiments/raw/analytic_e{1,2,3}.jsonl`. Summaries:
  `experiments/analysis/analytic_summary.json` and `analytic_e3_summary.json`.
- Candidate cards: `notes/analytic_cards.md`.
- GPT-6.1-sol:
  - blind round: `codex/analytic_1.{md,out}`, full report at
    `codex/scratch/analytic_1/analytic_1_report.md` (git-ignored);
  - adversarial round: `codex/analytic_2.{md,out}`.
  - Their findings are applied, see §10.
- Two read-only Opus critics, on correctness/cost and on literature: §10.
- Labels:
  - **KNOWN**, with a reference;
  - **PROVED**, with the proof here. Proofs are elementary; "no prior formulation found by us" at most;
  - **HEURISTIC**;
  - **CONJECTURED**;
  - **experimentally supported**;
  - **REFUTED for this model**;
  - **UNVERIFIED**.

## 0. Verdict (plain language)

**No analytic (exponential-sum, theta, L-function, lattice-point) route to factoring was found that beats
exponential time. None of the analytic ideas examined comes close to the deterministic record N^{1/5},
let alone to L[1/3].** This is a verdict about known techniques and the mechanisms we could build. It is
not a lower bound.

What we established:

1. **Hiary's N^{1/3} factoring algorithm contains no exponential sum.**
   - It is a continued-fraction method. For each block [x − H, x + H], one Dirichlet approximation of
     N/x² turns "x + h divides N" into a quadratic equation in h.
   - Hiary's fast quadratic and cubic exponential-sum algorithms (Annals 2011) belong to his separate work
     on computing ζ. The prompt's premise conflated the two.
   - The exponent 1/3 has three equivalent explanations:
     - a Dirichlet/curvature balance: block width H ≍ x·N^{−1/3};
     - an integer-hull count: the integer hull of {xy ≥ N} has Θ(N^{1/3} log N) vertices, and every
       divisor is one of them;
     - "1/3 of the bits of p", which is the 2-D lattice (Rivest–Shamir) strength.
   - E1 measured Hiary's block count at 0.98–0.99 of the continuum prediction. It also found
     0.48·N^{1/3} hull vertices in every dyadic window.
2. **Reduction L_j, the floor-sum ladder (PROVED; E3 experimentally supported).**
   - Suppose exact floor sums Σ⌊P(h)⌋ of rational degree-j polynomials were available in polylog time.
     Then all divisors are found with O(N^{1/(j+2)} log N) blocks and O(N^{1/(j+2)} log² N) oracle calls.
   - Concretely: j = 1 gives N^{1/3}, operational today with Euclid floor sums. j = 2 gives N^{1/4},
     j = 3 gives N^{1/5}, and j = 4 gives N^{1/6}. Uniform degree j ≍ √(log N/log log N) with
     T_j = e^{O(j log j)} gives L[1/2].
   - The measured slopes were 0.336, 0.252 and 0.201 for j = 1, 2, 3.
   - This corrects two earlier exponents for this construction: N^{1/(j+1)} (GPT-6.1-sol's blind round)
     and N^{1/(2j+2)} (`nonsmooth_results.md` §3).
   - The missing primitive is **exact lattice-point counting under short, nearly parabolic polynomial
     arcs, faster than the hull complexity H^{1/3}**. No such algorithm is known. Fourier/exponential-sum
     evaluation cannot supply it (barrier F below).
3. **Complete cubic Weyl sums mod N factor N (reduction C3, PROVED; E2 100/100).**
   - The sums are S(a, N) = Σ_{x mod N} e(ax³/N). Nine values to absolute precision O(1) factor every
     N = pq with p ≡ q ≡ 1 (mod 3), through Jacobi sums in Z[ω].
   - Whether S(1, N) vanishes decides the fixed-base QR bit "is −3 a square mod N".
   - By contrast, degree-2 complete sums are Jacobi symbols. So the reciprocity / theta-transformation
     wall sits exactly between degree 2 and degree 3.
   - Consequence: Hiary's cubic-sum algorithm (cubic coefficient ≤ K^{−2}) cannot extend to the
     complete-sum range with polylog cost unless a positive proportion of semiprimes are easy to factor.
4. **Dequantisation: no route.**
   - The QFT half of Shor is classically simulable on product inputs. All hardness is in the phase
     x ↦ h·a^x/N.
   - That phase is closed under Weyl differencing, so differencing never lowers its degree. It has no
     stationary-phase (B-process) structure.
   - Its bit-cut matrices have rank ≍ min(2^s, 2^{m−s}, r). The rank is a theorem for odd r; that the
     needed Fourier coefficients do not vanish is heuristic. So tensor-network contraction of Shor's sum
     costs ≳ r.
5. **Barrier map (§6).**
   - Integer-hull output size: Ω(N^{1/3} log N) for hull walks.
   - A band-limited divisor detector needs bandwidth ≥ 3x/8 − 1: the uncertainty principle, PROVED.
   - Ramanujan / singular-series truncations below the least prime factor are N-independent (PROVED).
   - A degree-j approximation-access model costs ≍ N^{1/(j+2)} (HEURISTIC model lower bound).
   - The cubic Weyl-sum hardness above.
   - None of these is a general factoring lower bound.

**Best idea.** Item 2's missing primitive is the one well-specified open summation problem: count lattice points
under a short polynomial arc of degree j with tiny higher coefficients exactly, in time H^{o(1)}. Even a
polylog evaluator at j = 2 only ties Pollard–Strassen and Coppersmith (N^{1/4}). Beating the record needs
j ≥ 4. L[1/2] needs uniformity in j. No mechanism for any of these is known.

## 1. Hiary's algorithms, dissected

### 1.1 The factoring algorithm (KNOWN)

Reference: G. A. Hiary, *A deterministic algorithm for integer factorization*, Math. Comp. 85 (2016),
2065–2069 (arXiv:1408.2608, read in full). Algorithm 1 works as follows.
- Trial-divide up to x₀ = ⌈(17N)^{1/3}⌉.
- Then, for blocks [x − H, x + H] with H = ⌊(17N)^{−1/3}x⌋:
  - take the convergent b/q of N/x² with the largest q ≤ 4H (Lemma 3.1:
    N/x² = b/q + ε₂/(qq′), q ≤ 4H ≤ q′);
  - set a = [qN/x];
  - solve (qN − ax) + (bx − a)h + bh² = 0 (Lemma 3.2);
  - test the ≤ 2 integer roots by division.
- Theorem 2.1: O(N^{1/3} log² N) operations on O(log N)-bit numbers.

**Which quantity reveals the factor:** none in the analytic sense. There is no sum, no bisection and no
precision issue. The exact identity

  qN/(x+h) = a − bh + E(h),  E(h) = ε₁ − ε₂h/q′ + qNh²/(x²(x+h)),  |ε₁| ≤ ½, |ε₂| < 1,

shows that E(h) is an integer at a divisor. The block condition makes |E(h)| < 1, so E(h) = 0.

**Why 1/3 (PROVED; the balance is Hiary's own estimate (4)).** The curvature term must satisfy
q·NH²/x³ < ¼. The Dirichlet denominator q is generically ≍ H, so the condition is N(H/x)³ ≲ 1/16,
i.e. H ≍ x N^{−1/3}.
- Blocks per dyadic range [X, 2X]: ≈ (17N)^{1/3}(ln 2)/2 = 0.891 N^{1/3}, independent of X.
- In total, ½(17N)^{1/3} ln(√N/(17N)^{1/3}) ≈ 0.214 N^{1/3} ln N.
- **The binding resource is the denominator-amplified linearisation error.** It is not sum length, phase
  degree, precision, or the cost of cubic sums.
- When N/x² happens to have a small-denominator approximation, H can be larger ("major arcs"). Hiary's
  implementation exploits this and gains about 20%. Generic blocks cannot.

**Three equivalent readings of the exponent.**
- (a) **2-D Diophantine strength.** Hiary notes that the method needs the top (1/3)log₂N bits of p. That
  is the Rivest–Shamir (EUROCRYPT 1985) 2-D integer-programming strength, against Coppersmith's 1/4.
  Width H ≍ N^{1/6} near √N means √N/N^{1/6} = N^{1/3} blocks.
- (b) **Integer-hull count** (KNOWN theorem + PROVED corollary).
  - The integer hull I(H_N) of {xy ≥ N} has Θ(N^{1/3} log N) vertices: Balog–Bárány, arXiv 2602.06897
    (Feb 2026, preprint), sharpening Alcántara–Blanco–Criado–Santos, arXiv 2501.19193. The latter has
    bounds 2N^{1/3} − 2 and O(N^{1/3} log N), plus a vertex enumerator at O(log) per vertex. Its stated
    motivation is deterministic factoring.
  - Every divisor (d, N/d) lies on the strictly convex boundary of H_N, hence is a vertex (PROVED).
  - So walking the hull factors N in Θ(N^{1/3} log N), and any method that outputs or visits all hull
    vertices costs Ω(N^{1/3} log N) (PROVED, trivially, as output size).
  - E1, exact: the hull vertex count per full dyadic window in [N^{1/3}, √N] is
    (0.480–0.490)·N^{1/3} at every size from 24 to 40 bits. The max/min window ratio is ≤ 1.03, and the
    slope of ln(V/ln N) is 0.331.
  - Hiary's blocks over the same range number 1.06–1.38 × the hull vertices.
- (c) **Curvature scale** (HEURISTIC taxonomy).
  - A convex arc of length L with radius of curvature ρ carries ≍ Lρ^{−1/3} hull vertices; this is the
    Jarník/Andrews/Balog–Bárány scale.
  - For xy = N near √N, L ≍ ρ ≍ √N, giving N^{1/3}. The same holds for x² + y² = N, whose points factor
    N when p ≡ q ≡ 1 (mod 4); the disc hull has ≍ R^{2/3} = N^{1/3} vertices (Balog–Bárány 1991).
  - Lehman's "Farey dissection" and Rubinstein's displacement-versus-occupancy balance (codex
    analytic_1 §1.5) also land at 1/3. That all these 2-D methods agree is suggestive, not a theorem
    (see §10 for the referee's view).

### 1.2 Hiary's exponential-sum algorithms (KNOWN; they are about ζ, not factoring)

- **Quadratic sums.** Hiary, *A nearly-optimal method to compute the truncated theta function, its
  derivatives, and integrals*, Ann. of Math. 174 (2011) 859–889. Σ_{k≤K} k^r e(ak + bk²) is computed to
  error ε in poly(log(K/ε)) operations, with no precomputation, by repeated theta (Poisson / Weil)
  transformations, i.e. a continued-fraction-like recursion in (a, b).
- **ζ(½ + it).** Hiary, *Fast methods to compute the Riemann zeta function*, Ann. of Math. 174 (2011)
  891–946.
  - Blocks of the Riemann–Siegel sum where t log(1 + k/v) is quadratic up to O(1) give t^{1/3}.
  - The main result, t^{4/13}, uses Theorem 2.2. That theorem covers cubic sums Σ k^r e(ak + bk² + ck³)
    with 0 ≤ c ≤ K^{μ−3}, μ ≤ 1 (so |c| ≤ K^{−2}), at polylog cost per query after a K^{4μ}
    FFT precomputation. Balancing t^{1/2−β} block queries against K^{4μ} gives β = 5/26, i.e. t^{4/13}.
  - The limit is the size of the cubic coefficient: the van der Corput iteration keeps the sums cubic
    only for small c.
- **Analogy (HEURISTIC).**
  - Riemann–Siegel (t^{1/2}) ↔ trial division (N^{1/2}).
  - Odlyzko–Schönhage amortised multi-evaluation ↔ Pollard–Strassen / Bostan–Gaudry–Schost multipoint
    evaluation (N^{1/4}).
  - Hiary's t^{1/3} quadratic blocking ↔ Hiary's N^{1/3} linear+CF blocking.
  - The analogy breaks at one point. For ζ, a numerical value with error O(1) is the goal. For
    factoring, an exact count is needed, and band-limited numerics cannot give exactness (barrier F).

## 2. Candidate table

Each "best known evaluation cost" is classical, without factoring first. A "minimal requirement" row is
what an evaluator must deliver for the factor to be read off.

| Candidate quantity | Identity linking it to a factor | Minimal requirement | Best known evaluation | Beats N^{1/3}? N^{1/5}? N^ε? |
|---|---|---|---|---|
| Restricted divisor count C_N(X) = D_N(X) − D_{N−1}(X), with D_N(X) = Σ_{d≤X}⌊N/d⌋ | C_N(X) = #{d ≤ X : d \| N}; bisection, O(log N) queries (PROVED; codex analytic_1 §2.1) | error < ¼ in each D | N^{1/3} (Hiary; hull walk; ladder j = 1) | no / no / no |
| Full D(N) = Σ_{m≤N} τ(m) | only τ(N): useless | — | N^{1/3} (Helfgott–Thompson; Sladkey arXiv 1206.3369). Li, arXiv 2609.29949, *claims* N^{9/28}, but uses Lenstra–Pomerance factoring to de-smooth (UNVERIFIED, preprint) | not applicable |
| Sawtooth Σ_{d∈I} ψ(N/d), or Σ_d e(hN/d) | jump of ½ at a divisor; the Fourier series loses the endpoint term (codex analytic_1 §2.2) | all frequencies \|h\| ≲ x | each single-h sum: Õ((\|h\|N)^{1/4}) by quadratic theta blocking (codex); all frequencies ≳ x sums (barrier F) | no |
| Hyperbola / Voronoi dual sums | dual length N/H² at smoothing width H; de-smoothing re-exposes an H-interval (codex) | exact | no restricted-cutoff version below N^{1/3} known | no |
| r₂(N) via θ² | r₂(pq) ∈ {0, 16}: one QR-type bit only; two inequivalent representations would factor | points on x² + y² = N | hull walk on the disc: N^{1/3} (§1.1(c)) | no |
| r₄(N) = 8σ(N) (N odd) via θ⁴ | s = p + q = r₄/8 − N − 1 (PROVED, KNOWN) | error < 4 (additive N^{1/4} gives a Coppersmith hint) | coefficient extraction needs ~N samples or exact singular series (barrier R) | no |
| σ(N) via Ramanujan sums | every corrected term with modulus < p vanishes (nonsmooth §4.1; general form: barrier R) | moduli ≥ p | √N | no |
| h(−4N) via L(1, χ_{−4N}) | class group + ambiguous forms (Schoof 1982) | \|δL\| < π/(2√Q) | AFE √Q = N^{1/2}; Shanks BSGS hybrid N^{1/5} under GRH (Sutherland thesis §11.4, via codex) | GRH: ties 1/5 |
| Complete cubic Weyl sum S(a, N) | 9 values → p (C3, PROVED) | absolute error O(1) on a quantity of size 4√N | trivial O(N); O(p + q) given the factors | no |
| Interval of Legendre symbols Σ_{a<N/2}(a/N) | h(−N) (class b) | exact | √N | no |
| Shor sums Σ_x e(h a^x/N) | period of a ↦ order ↦ factor (Miller, randomised) | — | ≳ √r generic; ≳ r for tensor contraction (§7) | no |

## 3. Reduction L_j: the floor-sum ladder

**Theorem L_j (PROVED, conditional on an oracle).** Fix j ≥ 1 and c ∈ (0, ½]. Suppose an oracle returns
F(P; h₁, h₂) = Σ_{h₁≤h<h₂} ⌊P(h)⌋ exactly, for polynomials P ∈ Q[h] of degree ≤ j with O(j log N)-bit
coefficients, in time T_j. Then all divisors of N in [B/c, √N], B = ⌊N^{1/(j+2)}⌋, are found with
- O(N^{1/(j+2)} log N) blocks;
- O(N^{1/(j+2)} log² N) oracle calls;
- O(N^{1/(j+2)} log N) trial divisions,

plus trial division below B/c. Total: N^{1/(j+2)}·(log N)^{O(1)}·T_j.

*Proof.*
1. **Blocks.** Use blocks [x, x + H) with H = ⌊cx/B⌋ ≥ 1, x increasing from B/c by H. Per dyadic range
   there are ≤ (B/c)·ln 2 + O(1) blocks, so O(B log N) in total.
2. **Taylor polynomial.** On a block, f(h) = N/(x+h) = P(h) + R(h). Here
   P(h) = Σ_{i≤j}(−1)^i N h^i/x^{i+1} (denominator x^{j+1}) and R(h) = N(−h)^{j+1}/(x^{j+1}(x+h)).
   This follows from the geometric series.
3. **Remainder bound.** For 0 ≤ h < H, R(h) has the constant sign (−1)^{j+1} (R(0) = 0), and
   |R(h)| ≤ ρ := NH^{j+1}/x^{j+2} ≤ c^{j+1}N/(B^{j+1}x) ≤ e·c^{j+1}B/x ≤ e·c^{j+2} < 1. This uses
   N < (B+1)^{j+2} ≤ e·B^{j+2} for B ≥ j + 2, and x ≥ B/c.
4. **Divisors are candidates.** If x + h | N then f(h) ∈ Z, so:
   - j odd: P(h) = f − R gives {P(h)} ∈ [1 − ρ, 1) or h = 0;
   - j even: P = f + |R| gives {P(h)} ∈ [0, ρ].
5. **Counting candidates.** The candidate set C is counted on any sub-range by one difference of two
   oracle calls. For j odd: Σ⌊P + ρ⌋ − Σ⌊P⌋. For j even: Σ⌊P⌋ − Σ⌊P − ρ − η⌋, where η is smaller than the
   spacing of the rationals involved.
6. **Locating candidates.** Bisection lists C with O(|C| log H + 1) calls. Each candidate is tested by
   one division.
7. **Size of C.** Every candidate d = x + h satisfies {N/d} ∈ [1 − ρ, 1) ∪ [0, ρ). This holds because
   f = P ± |R| with |R| ≤ ρ. That is, d | N − k for some |k| ≤ ρd ≤ e·c^{j+1}B(x + H)/x ≤ 2e·c^{j+1}B
   ≤ Y := B (for c ≤ ½). So the total candidate count over all blocks is
   ≤ Σ_{|k|≤Y} τ(N − k) ≪ Y log N. This uses Shiu's Brun–Titchmarsh theorem for τ in short intervals
   (Shiu, J. reine angew. Math. 313 (1980), 161–170: Σ_{x<n≤x+y} τ(n) ≪ y log x for y ≥ x^ε). Hence
   O(B log N) candidates and O(B log² N) calls. ∎

The trivial bound Σ_{|k|≤Y} τ(N − k) ≤ (2Y+1)·max τ = Y·N^{o(1)} already gives N^{1/(j+2)+o(1)}
without Shiu.

**Why this is N^{1/(j+2)} and not N^{1/(j+1)}.**
- Requiring the polynomial to reproduce every floor (precision 1/(2x), codex analytic_1 §1.4) forces
  ρ ≲ 1/x, hence H ≍ xN^{−1/(j+1)}.
- Requiring only ρ ≲ 1/H, and *locating* the expected O(1) ambiguous h per block with the same oracle,
  gains one degree.
- At j = 1 the ladder is an N^{1/3} algorithm built only from Euclid floor sums. This agrees with Hiary,
  with the hull walk, and with every known N^{1/3} divisor-sum method. The precision-1/x construction
  would give only N^{1/2} at j = 1.
- nonsmooth's N^{1/(2j+2)} assumes that remainder < 1 suffices. Then a constant fraction of the h are
  ambiguous and cannot be resolved cheaply. That exponent is **REFUTED for this construction** (E3).

**E3 (pre-registered amendment A1; all predictions held).** The j = 1 oracle is operational (Euclid). The
j = 2, 3 oracles are simulated by brute force; only call counts are model costs.

| j | sizes (bits) | correct | slope of ln(blocks/ln N) vs ln N | 1/(j+2) | codex 1/(j+1) | nonsmooth 1/(2j+2) | located candidates / block | calls / block |
|---|---|---|---|---|---|---|---|---|
| 1 | 24–48 | 70/70 | **0.336** | 0.333 | 0.5 | 0.25 | 0.064–0.093 (c³ = 0.125) | 2.19–2.68 |
| 2 | 24–40 | 50/50 | **0.252** | 0.25 | 0.333 | 0.167 | 0.038–0.053 (c⁴ = 0.0625) | 2.22–2.46 |
| 3 | 24–40 | 50/50 | **0.201** | 0.2 | 0.25 | 0.125 | 0.013–0.034 (c⁵ = 0.031) | 2.10–2.40 |

For primes with n ≥ 32, the ratio blocks/continuum-prediction falls from 1.19 to 1.07; it decays toward 1
as the rounding loss in ⌊cx/B⌋ shrinks.

**What E3 does and does not test.**
- The block count is fixed by construction (H = ⌊cx/B⌋), so the slope column only confirms the counting
  formula.
- The substantive content is threefold:
  - (a) the algorithm stays correct at these block sizes (170/170, with every divisor caught as a
    candidate);
  - (b) located candidates per block stay O(1), as the Shiu-type bound predicts;
  - (c) oracle calls per block stay ≈ 2–3.
- Together these show that the total model cost tracks the N^{1/(j+2)} block count. A construction with
  block sizes from N^{1/(j+1)}, or from nonsmooth's N^{1/(2j+2)}, would show either needlessly many
  blocks or ≍ H candidates per block.
- E3 says nothing about the cost of the j ≥ 2 oracles, which are simulated.

**Uniform degree → L[1/2] (PROVED, conditional; same calculation as codex analytic_1 §4.1, with j+2 in
place of j+1).** If T_j ≤ exp(Cj log(j+2))·poly(log N), the log-cost is
log N/(j+2) + Cj log j + O(log log N). Taking j ≍ √(log N/(C log log N)) gives L_N[½, O(√C)]. If
T_j = exp(O(j^a)), the exponent becomes (log N)^{a/(a+1)}.

**What the missing primitive looks like (PROVED, a calculation).**
- At j = 2, near x ≍ √N, H ≍ N^{1/4}. P has slope ≍ 1 and quadratic coefficient β = N/x³ with
  βH² ≍ 1. So it is a *full-curvature short parabola*.
- The lattice points under such an arc have ≍ H^{1/3} hull vertices. A hull walk costs H^{1/3} per block,
  i.e. N^{1/4}·N^{1/12} = N^{1/3}: exactly Hiary again. The hull law (§1.1) is self-similar across the
  ladder.
- So a gain at j ≥ 2 requires counting lattice points under a curved arc in time below its hull
  complexity. We know no such method. Barrier F excludes the Fourier route.
- Complete quadratic floor sums are in a different, class-number-hard regime. For prime m ≡ 3 (mod 4),
  m > 3, Σ_{x=1}^{m−1}⌊x²/m⌋ = (m−1)(2m−1)/6 − (m−1)/2 + h(−m), via Dirichlet's
  Σ_{r<m} r(r/m) = −m h(−m) (PROVED from KNOWN; checked numerically in §9).

**Model lower bound D_j (HEURISTIC).** Consider methods that see a block only through a degree-j
approximation of N/(x+h), and pay 1 per block plus 1 per ambiguous h.
- By Chebyshev, E_j(f; [x, x+H]) ≥ H^{j+1} min|f^{(j+1)}| / (2^{2j+1}(j+1)!) ≍_j NH^{j+1}/x^{j+2}.
- The ambiguous h number ≈ H·E_j per block under equidistribution of {P(h)}.
- Minimising (x/H) + x·NH^{j+1}/x^{j+2} gives ≍_j N^{1/(j+2)} per dyadic range near √N.
- The equidistribution step is why this is HEURISTIC. It explains why the ladder cannot be improved
  within the model; it says nothing about other models (Coppersmith uses exact mod-p lattices,
  Pollard–Strassen uses multipoint evaluation).

## 4. Reduction C3: complete cubic Weyl sums mod N factor N

**Setting.** S(a, N) = Σ_{x mod N} e(ax³/N), with N = pq, p ≠ q primes > 3 and gcd(a, N) = 1.

**Lemma 1 (KNOWN, standard).**
- (i) S(a, pq) = S(aq̄, p)·S(ap̄, q), with q̄q ≡ 1 (mod p) and p̄p ≡ 1 (mod q). By CRT,
  1/(pq) ≡ q̄/p + p̄/q (mod 1) after multiplying by ax³.
- (ii) If p ≡ 2 (mod 3), x ↦ x³ permutes F_p, so S(b, p) = 0.
- (iii) If p ≡ 1 (mod 3), let χ be a cubic character. Counting cube roots gives
  S(b, p) = χ̄(b)G(χ) + χ(b)G(χ̄). Since χ(−1) = 1, G(χ̄) = \overline{G(χ)}. So
  S(b, p) = 2√p cos(θ_p − 2πk(b)/3), where G(χ) = √p e^{iθ_p} and χ(b) = e(k(b)/3).
- (iv) G(χ)³ = p·J(χ, χ) with J(χ, χ) = π_p ∈ Z[ω] of norm p (Ireland–Rosen, Ch. 9). So
  π_p/√p = e^{3iθ_p}.

**Lemma 2 (PROVED). The entries never vanish.** cos(θ_p − 2πk/3) = 0 would force e^{3iθ_p} = ±i, i.e.
π_p ∈ √−3·Q, and then N(π_p) = 3a² ≠ p.

**Theorem C3 (PROVED).** For N = pq with p ≡ q ≡ 1 (mod 3) and random units a, b, the nine values
S(a^i b^j, N) (0 ≤ i, j ≤ 2) determine p in poly(log N) exact operations, with success probability 48/81 per
pair (a, b). This holds both for exact values and for values with absolute error ≤ c₀, for an absolute
constant c₀ > 0.

*Proof.*
1. **The 3×3 array.** By Lemma 1, S(a^i b^j, N) = 4√N u_{α(i,j)} v_{β(i,j)}, with
   u_k = cos(θ_p − 2π(k + s_p)/3) and α(i, j) = i·k_p(a) + j·k_p(b) linear over F₃ (similarly β for q).
   With probability |GL₂(F₃)|/81 = 48/81, (α, β): F₃² → F₃² is bijective. The decoder tries all 48
   ordered pairs of independent functionals. For the right one, the array is the rank-one matrix
   4√N u vᵀ.
2. **Recovering a direction.** Take u from a column with |v_l| maximal (≥ 1/√2, since Σv_l² = 3/2).
   Then z_u := Σ_k u_k ω^k = (3/2)e^{iθ_p}ω^{−s_p}·(a real scalar). Hence z_u³ is a real multiple of
   π_p. Relabelings in the dihedral group S₃ change π_p to an associate or to π̄_p.
3. **Two norm-N elements.** Likewise z_v³ is a real multiple of π_q. Then γ = z_u³z_v³ is a real
   multiple of π_pπ_q (or of its conjugate), and γ′ = z_u³·\overline{z_v³} is a real multiple of
   π_pπ̄_q (up to conjugation). Rescaling to norm N recovers them exactly: λ² = N/N(γ), then round in
   Z[ω].
4. **The factor.** gcd_{Z[ω]}(π_pπ_q, π_pπ̄_q) = π_p, since q splits into non-associate primes. Its norm
   is p.
5. **Wrong labelings and noise.** Wrong labelings give elements that fail the norm or divisibility test.
   The output is verified by p | N, so there are no false positives.
6. **Precision.** First order: the coordinate error after rescaling is ≤ ~4δ, so δ ≤ c₀ ≈ 0.05
   (worst case) suffices. A ±1 neighbourhood search widens this. ∎

**Corollaries (PROVED).**
- (i) **Fixed-base QR bit.** For N ≡ 1 (mod 3), S(1, N) ≠ 0 iff p ≡ q ≡ 1 (mod 3), i.e. iff −3 is a square
  mod N. This is the Track E / nonsmooth §5.5 minimal-primitive bit.
- (ii) **Hardness of a cubic analogue of Hiary.** A classical algorithm approximating S(a, N) to absolute
  error c₀ in time T factors every such N (about a quarter of balanced semiprimes) in expected
  ≈ 15 oracle calls plus poly(log N). Hiary's cubic theorem has |c| ≤ K^{−2} (μ ≤ 1). The complete sum
  has c = a/N with K = N (μ = 2). **A polylog cubic evaluator for the full coefficient range is therefore
  at least as hard as factoring these N.**
- (iii) **The degree-2 / degree-3 dividing line.** Complete quadratic sums mod N are Jacobi symbols
  times ε_N√N (class c: reciprocity-computable, 0 bits beyond N mod 4). So the line between "closed-form
  via reciprocity / Weil representation" and "factoring-hard" falls exactly between degree 2 and 3.
- **Analogue for degree k (CONJECTURED, sketch).** Degree-k sums at N with ℓ | gcd(p − 1, q − 1) for
  some ℓ | k would work the same way through Jacobi sums in Z[ζ_ℓ]. When gcd(p − 1, q − 1, k) ≤ 2 (e.g.
  safe primes), every such sum reduces to the quadratic case and carries nothing.
- Prior art: van Dam–Seroussi (quant-ph/0207131) reduce discrete log to Gauss-sum estimation for
  multiplicative characters. We found no prior statement of C3 for additive cubic sums at composite N:
  "no prior formulation found by us"; the literature critic's verdict is in §10.

**E2 (pre-registered).** Balanced N, 14–22 bits.

| bits | res1 (p ≡ q ≡ 1) success | structural | accidental | res2 control max \|S\| | res2 structural | δ* (≥ 50% success under uniform noise ±δ) |
|---|---|---|---|---|---|---|
| 14 | 20/20 | 15 | 5 | 6·10⁻¹³ | 0 | 10 |
| 16 | 20/20 | 20 | 0 | 2·10⁻¹² | 0 | 3 |
| 18 | 20/20 | 20 | 0 | 1·10⁻¹¹ | 0 | 3 |
| 20 | 20/20 | 20 | 0 | 4·10⁻¹¹ | 0 | 3 |
| 22 | 20/20 | 20 | 0 | 2·10⁻¹⁰ | 0 | 3 |

- **P2.1 and P2.2 hold.**
- **P2.3 misses as specified.** δ* = 10 at 14 bits lies outside [0.03, 3], and the slope is −0.170,
  just outside [−0.15, 0.15]. Success at δ = 10 falls from 0.6 to 0.17 with size, while δ ≤ 1 stays at
  ≥ 0.97 everywhere.
  - Reading: the decoder needs constant absolute precision, slightly shrinking with N. It needs all but
    ≈ 2 of the ≈ ½log₂N + 2 significant bits.
  - The data rule out the "half the bits" scenario (δ* ∝ N^{1/4}) for this decoder. A lattice decoder
    exploiting the norm-N conic might tolerate N^{c} error; not tried.
- **Accidental gcd hits (disclosed).**
  - At 14 bits there is a small prime pool, and N = 13843 = 109·127 recurs three times. One instance
    hit 5 structural failures in a row; we checked that the structural rate per attempt there is 0.575,
    against the predicted 48/81 = 0.59.
  - A cosmetic bug: the early exit for "all S = 0" compared floats to exact 0. Controls therefore ran all
    6 attempts, which inflated their accidental hits (4/20 at 14 bits). It has no effect on any structural
    conclusion.

## 5. Dequantisation verdict (classical only)

**What would be needed.** Shor's amplitude is A(k) = M^{−1} Σ_{x<M, a^x ≡ y} e(kx/M), which equals
(MN)^{−1} Σ_{h mod N} Σ_{x<M} e(kx/M + h(a^x − y)/N).

- **The QFT is not the obstruction (KNOWN).**
  - Semiclassical QFT: Griffiths–Niu, PRL 76 (1996).
  - Efficient classical simulation of the QFT on product or weakly entangled inputs: Aharonov–Landau–
    Makowsky (quant-ph/0611156); Yoran–Short, PRA 76 (2007) 042321; Browne, NJP 9 (2007) 146.
  - Exact citation details are UNVERIFIED by us beyond the abstracts. The hardness lies in the
    modular-exponentiation phase.
- **No degree reduction (PROVED).** Δ_u(h a^x/N) = h(a^u − 1)a^x/N. The class {αx + h a^x/N} is closed
  under Weyl differencing, which only rescales h (codex analytic_1 §3.1 gives the same identity). The
  A-process of van der Corput therefore never terminates in a quadratic sum. The real-variable phase
  h a^x/N has derivatives growing like (log a)^k a^x/N, so stationary phase (the B-process) has no
  shorter dual. Reduced residues instead introduce carries. This part is HEURISTIC as a statement about
  all transformations.
- **Rank across bit cuts (PROVED for odd r, modulo a non-vanishing HEURISTIC).** Write
  x = 2^s x_hi + x_lo and G = ⟨a⟩ (cyclic, order r).
  - Then e(h a^x/N) = f(g_hi g_lo), with f = e(h·/N)|_G = Σ_{ψ∈Ĝ} f̂(ψ)ψ.
  - So the matrix [f(g_hi g_lo)] equals Ψ_hi·diag(f̂)·Ψ_loᵀ, whose factors are character Vandermonde
    matrices. For odd r the elements a^{2^s t} are distinct, and the rank is
    min(#distinct g_hi, #distinct g_lo, #supp f̂).
  - Generically f̂(ψ) ≠ 0: these are Gauss-type periods, and their non-vanishing is HEURISTIC.
  - So the Schmidt rank of Σ_x e(h a^x/N)|x⟩ at cut s is ≍ min(2^s, 2^{m−s}, r). Any MPS / tensor
    contraction along the bit order needs bond dimension ≍ r, since m ≥ 2 log₂ r in Shor. This is
    consistent with Jozsa–Linden (Proc. R. Soc. A 459 (2003)): entanglement unbounded in block size is
    necessary.
- **Query bounds (KNOWN, black-box only).**
  - Cleve, Inf. Comput. 192 (2004) 162–171: Ω(r^{1/3}/√log) for order finding of a hidden permutation.
  - Chakraborty–García-Soriano–Matsliah (ICS 2011): Ω̃(√r) with polynomially bounded jumps.
  - Explicit modular multiplication is not a black box, but the only non-black-box structures known to
    be exploited are smoothness (p − 1, ECM, index calculus).
- **Structured ranges (KNOWN, codex).** Complete sums over a subgroup of known index d reduce to Gauss
  sums, e.g. d = 2 gives quadratic Gauss sums. This needs the subgroup and its characters, i.e. the
  answer.
- **Verdict.** Classical evaluation fails at a precise point: the phase h·a^x/N is an eigenfunction of
  shifts (closed under differencing) rather than a low-degree polynomial, and its bit-cut rank equals the
  hidden order. Every fast-summation technique we know (Weyl/van der Corput, Poisson/theta, Taylor
  blocking, Odlyzko–Schönhage amortisation, tensor contraction) needs low degree or low rank. **No
  classical piece of the period information is extractable by these techniques.** It is not a theorem
  against all classical algorithms.

## 6. Barrier map

| # | Statement | Scope | Status |
|---|---|---|---|
| H | Hull-walk methods (Hiary, the Alcántara et al. enumerator, hull walks on x² + y² = N) cost Ω(N^{1/3} log N) | methods that visit every integer-hull vertex in [N^{1/3}, √N] | PROVED (output size; vertex count KNOWN, Balog–Bárány 2026 preprint); E1 |
| F | **Uncertainty principle for divisor detection.** If F(u) = Σ_{\|k\|≤K} c_k e(ku), F(0) = 1, \|F(r/d)\| ≤ ε < 1/(2W) for all d ∈ I ⊂ [x, 2x], \|I\| = W ≥ √(2x), r ≢ 0, then K ≥ 3x/8 − 1. Proof: sample at r/d (no aliasing if K < d/2). The DFT of δ₀ + e must vanish at d − 2K − 1 frequencies, each needing \|ê\| = 1/d. Parseval gives (d − 2K − 1)/d² ≤ (d − 1)ε²/d, so 2K + 1 ≥ d − d(d−1)ε² ≥ 3d/4. A stable variant: Σ\|c_h\| ≤ A forces HA ≳ d (codex analytic_1 §5.2) | linear detectors Σ_{d∈I} F(N/d) with a uniform error budget, evaluated frequency by frequency; not exact-arithmetic methods (Hiary, the ladder) and not detectors that exploit error cancellation | PROVED |
| R | **Ramanujan truncation.** If f(N) = Σ_c a_c(\|N\|)c_c(N), the truncation at c < P⁻(N) is Σ_{c<P⁻} a_c(\|N\|)μ(c): the same function of \|N\| as at a prime. This applies to σ (nonsmooth §4.1) and to the odd-q singular-series terms of r_k for even k (G(a, q)^k = ε_q^k q^{k/2} for odd q). For odd k the terms carry (N/q): information through L(·, χ_N) (class b) or N mod q (floor F1) | series whose N-dependence is only through Ramanujan sums | PROVED |
| L | Degree-j approximation-access model: ≍_j N^{1/(j+2)} blocks + ambiguities per dyadic range near √N (§3) | that model | HEURISTIC (equidistribution step) |
| C3 | Polylog O(1)-accurate complete cubic Weyl sums mod N ⇒ factoring of N with p ≡ q ≡ 1 (mod 3) | full coefficient range (μ = 2) | PROVED; E2 |
| Q | Complete quadratic floor sums mod prime m encode h(−m) | complete period | PROVED from KNOWN (Dirichlet); checked |
| P | Precision: resolving {N/d} = 1/d needs Ω(log N) fractional bits; summing L independently rounded terms needs log₂L extra bits. **No** universal Ω(log N)-bits-per-term rule exists (exact arithmetic avoids it) | float classifiers / independent-error summation | PROVED (codex analytic_1 §5.3) |
| S | σ(n) is not P-recursive (every root of unity is a singularity of Σσ(n)zⁿ), so no fixed-order creative-telescoping recurrence exists for σ | holonomic methods | PROVED (codex analytic_1 §4.3; we checked the argument) |
| D | Shor-sum phase closed under differencing; bit-cut rank ≍ min(2^s, 2^{m−s}, r) | Weyl/vdC, MPS contraction | PROVED / HEURISTIC (§5) |
| I | Information: the floor F1 (functions of N mod m carry 0 bits) is the only proved information statement. No sum-length lower bound from information theory exists, because N determines its factors | — | KNOWN (Track E) |

**What is NOT a barrier.**
- Theta transformations make degree-2 sums polylog for every coefficient size, and Hiary's cubic algorithm
  handles small cubic coefficients. The obstruction is *exactness* (F) plus *degree ≥ 3 in the full range*
  (C3). The phase degree alone is not the obstruction.
- Codex correctly notes that complete quadratic *Gauss* sums do not encode class numbers. Our class-number
  statement (Q) concerns *floor* sums.

## 7. Subexponential analytic ideas: feasibility

| Idea | Feasibility calculation | Verdict |
|---|---|---|
| Uniform-degree ladder (§3) | L[½, O(√C)] iff exact degree-j floor sums cost e^{O(j log j)} | open primitive; even j = 2 polylog is unknown; the hull walk gives H^{1/3}; Fourier is barred (F) |
| Many-variable theta sums (Siegel modular transformations) | Enumeration cost (b + g)^{g/2} for genus g (Deconinck et al., Math. Comp. 73 (2004)). Encoding xy = N needs the quartic penalty e^{−λ(xy−N)²} (codex B) | killed: quartic, and coefficient extraction survives the transform |
| Gaussian product isolation M₁/Z (codex analytic_1 §6B) | exact for λ = 4 log N + c, but quartic phase; Fourier linearisation returns an e^{−itN} extraction integral | killed |
| Fermat/square-sieve exponential sums | detecting a square among T values by characters needs ≳ T characters (second moment); pseudosquare wheels give T/2^k with T ≍ √N | exponential |
| Phases polynomial after a factor-encoding substitution | needs the centre of the substitution (≈ p to N^{1/4}): the Coppersmith hint itself (codex §4.1) | circular |
| L-values / class numbers (AFE, Shanks hybrid) | √Q direct; Q^{1/5} with GRH and BSGS; h is class b at primes | exponential |
| Holonomic / hypergeometric (BGS) | N^{1/4} (Pollard–Strassen / BGS); σ is not D-finite (S) | exponential |
| Voronoi / Li 9/28 | full D(N) only; de-smoothing uses smooth-number factoring | not applicable (relational) |

No analytic sum was found whose fast evaluation is *plausible* and would give L[½] or better.

## 8. Best idea and the sharpest open question

**Open question (the minimal analytic primitive).** Count exactly #{(h, m) : 0 ≤ h < H, m ≤ αh + βh² + γ},
with rational α, β, γ, βH² ≍ 1 and denominators ≈ H^{O(1)}, in time H^{o(1)}. More generally, do so for
degree-j arcs with tiny higher coefficients.

What answers to it would give:
- Polylog at j = 2: deterministic N^{1/4} (a new, analytic route to an old exponent).
- Polylog at j ≥ 4: below the N^{1/5} record.
- Uniform in j: L[½].

Known facts:
- Euclid solves j = 1.
- The hull walk costs H^{1/3}.
- Complete-period quadratic floor sums are class-number-hard (Q).
- Exponential-sum evaluation cannot substitute (F).
- The full-range cubic analogue is factoring-hard (C3).

Our assessment (HEURISTIC): a polylog j = 2 counter would be a breakthrough in lattice-point
computation in its own right. We see no mechanism for it.

## 9. Experiments: pre-registration outcomes

**E1** (`run_experiments e1`; 140 Hiary runs, 25 hull runs).
- P1.1 140/140. ✓
- P1.2 blocks/prediction ∈ [0.978, 0.992]. ✓
- **P1.3 MISS.** The slope was 0.3478 against [0.32, 0.345]. The window ignored the finite-size factor
  ln(N^{1/6}/17^{1/3}) = ln N/6 − 0.944. P1.2, which encodes the exact count, holds.
- P1.4: V/(N^{1/3} ln N) = 0.120–0.125, slope 0.331, window max/min ≤ 1.03. ✓
- P1.5 25/25. ✓
- P1.6: V/blocks = 0.72–0.94. ✓
- Decision rule: claim 1 is supported in substance; the P1.3 miss is reported as a window
  mis-specification.

Hiary means for primes:

| bits | blocks | CF steps/block | blocks/(N^{1/3} ln N) |
|---|---|---|---|
| 24 | 533 | 3.80 | 0.138 |
| 32 | 4 820 | 4.23 | 0.157 |
| 40 | 41 057 | 4.65 | 0.169 |
| 48 | 328 960 | 5.06 | 0.176 |

Hull vertices per full dyadic window / N^{1/3}: 0.490, 0.483, 0.480, 0.486, 0.481 (24–40 bits).

**E2.** §4: P2.1 ✓, P2.2 ✓, P2.3 miss (constant absolute precision confirmed qualitatively; the window
was too tight).

**E3.** §3: P3.1–P3.5 all ✓.

**Check of Q** (scratch, exact): for every prime 7 ≤ m < 2000 with m ≡ 3 (mod 4) (154 primes, 0
failures), Σ_{x<m}⌊x²/m⌋ = (m−1)(2m−1)/6 − (m−1)/2 + h(−m), with h from the reduced-form count. The
derivation is Σ(x² mod m) = 2Σ_{QR r} r = m(m−1)/2 − m·h(−m), from Dirichlet's class-number formula.

## 10. Reviews

Filled in after the GPT-6.1-sol adversarial round and the two Opus critics.

## 11. Reproduction

```bash
.venv/bin/python -m pytest -q tests/test_analytic.py
PYTHONPATH=. .venv/bin/python -m src.analytic.run_experiments e1        # ~6 min
PYTHONPATH=. .venv/bin/python -m src.analytic.run_experiments e2        # ~8 min
PYTHONPATH=. .venv/bin/python -m src.analytic.run_experiments e3        # ~10 min
PYTHONPATH=. .venv/bin/python -m src.analytic.run_experiments summary
PYTHONPATH=. .venv/bin/python -m src.analytic.run_experiments summary_e3
```
