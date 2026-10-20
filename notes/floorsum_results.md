# floorsum: fast exact floor sums for the ladder (algorithm side)

Agent `floorsum`, 2026-10-20. DRAFT pending the GPT-6.1-sol adversarial round and two Opus critics
(§10 will list their fixes).

**Files.**
- Pre-registration: `notes/floorsum_prereg.md` (commit `a1a98a4`, before any experiment;
  Amendment A1 for E5 in `6d5b3ab`, before E5, with a disclosed peek).
- Code: `src/floorsum/farey.py` (exact oracle), `src/floorsum/ladder_ops.py` (operational ladder,
  D(N) via the ladder), `src/floorsum/lattice_block.py` (Coppersmith block oracle, KNOWN),
  `src/floorsum/run_experiments.py`.
- Tests: `tests/test_floorsum.py` (18 pass).
- Raw data: `experiments/raw/floorsum_e{1,2,2b,3,4,4b,5}.jsonl`; summary
  `experiments/analysis/floorsum_summary.json`; logs `experiments/logs/floorsum_*.log`.
- Cards: `notes/floorsum_cards.md`.
- GPT-6.1-sol blind round: `codex/floorsum_1.{md,out}`.

**Labels.** KNOWN (with reference), PROVED (proof here), experimentally supported, HEURISTIC,
CONJECTURED, REFUTED for this model, OPEN. "No prior formulation found by us" is the strongest
priority claim made.

## 0. Verdict (plain language)

**No fast exact floor-sum oracle was found. The best exact oracle we have, for every degree j, costs
about H·β^{1/3} per block. Plugged into the ladder this gives N^{1/3} for every j: exactly the trap.
So nothing here beats Hiary, let alone N^{1/4} or the record N^{1/5}.** This is a statement about the
algorithms we could build and the models we could analyse; it is not a lower bound for factoring.

What we established:

1. **The trap is real, and it is proved in two natural models** (§2).
   - Per block the convexity scale is Hβ^{1/3} ≍ cN^{1/3}/B. Over the O(B log N) blocks that is
     Θ(N^{1/3} log N), for every j ≥ 2 (PROVED).
   - For even j, every vertex of the integer hull of {xy ≥ N} is a vertex of some block hull,
     except at candidates (Lemma V, PROVED). So Σ_blocks V_b ≥ V_glob − O(B log N) = Ω(N^{1/3} log N).
   - Two lower bounds follow (PROVED): any algorithm that sees the block region only through
     lattice-membership or lattice-LP queries needs ≥ V_b queries per exact floor sum (Theorem T); and
     any exact method that splits the block into pieces on which ⌊P⌋ is the floor of one rational
     linear function needs ≥ V_b/O(log H) pieces (Prop PL). That covers hull walks, Hiary-type CF
     blocking and our Farey oracle.
   - E1 (pre-registered, all predictions held): for j = 2, 3, 4 the summed block-hull vertex counts
     grow like N^{0.31} (28–44 bits) and are 1.04–1.96 times the global hyperbola hull count.
2. **The best exact oracle: Farey dissection** (§3; operational, any degree, shifts and
   subintervals). Exactness rests on a one-line lemma (Lemma R). Its cost is
   O((1 + Hβ^{1/3}) log H) evaluations (Theorem O, PROVED for blocks whose P″ varies by a bounded
   factor, which includes all ladder blocks). At x ≍ √N this is T_2 ≈ H^{1/3}, T_3 ≈ H^{4/9},
   T_4 ≈ H^{1/2}, i.e. T_j ≈ H^{2(j−1)/(3j)}. The targets were H^{1/3−ε} (j = 2) and H^{1/10} (j = 4).
   E2b measured slopes 0.337, 0.442, 0.506 against 1/3, 4/9, 1/2. The cost is smaller on major arcs
   (slope near a fraction with small denominator).
3. **End to end** (§7, E3). The ladder becomes fully operational for j = 2, 3 (80/80 correct, no
   simulated oracle). Its cost grows like N^{0.29–0.30}/ln²N over 28–44 bits, i.e. N^{1/3} up to logs.
   j = 3 is not cheaper than j = 2 (evaluation ratio 0.85–0.90).
4. **Fourier smoothing cannot help** (§4, Prop F, PROVED in the certified smoothing model).
   - Any degree-K trigonometric approximant of the sawtooth is off by ≥ 1/4 on a set of measure
     ≥ 1/(51K) made of ≤ 2K + 2 intervals (Bernstein + Remez).
   - Every h landing there must be corrected individually. If the block's points {P(h)} have maximal
     gap g, the cost is ≥ 0.099/√g − 3.
   - E4: for ladder blocks the near-integer counts are 0.998–1.003 × 2δH, and gH ≈ 5–8. So the cost is
     ≍ √H per block up to logs (measured min_K(K + n_{1/K}) ≈ 2.9√H).
   - Hence Fourier + exact correction costs ≥ B√H ≈ N^{1/4 + 1/(2(j+2))}, which never beats N^{1/4}
     and is worse than the hull at j = 2. Hiary's theta algorithm supplies each Weyl sum in polylog
     time to inverse-polynomial error (KNOWN); exactness, not the Weyl sums, is the obstruction.
5. **"Euclid for quadratics" does not exist in the obvious sense** (§5).
   - The axis swap maps a quadratic floor sum to a sum of floors of square roots, and swapping back
     returns the original sum (PROVED, identity).
   - Unimodular moves preserve affine arc length λ = ∫|P″|^{1/3}dh (PROVED). Linear Euclid works
     because a line has λ = 0. A ladder block arc has λ ≍ cN^{1/3}/B, which no lattice move reduces.
   - Unrestricted quadratic floor sums are NP-hard (Manders–Adleman, KNOWN, via a two-call reduction
     found independently by codex floorsum_1 and floorhard_1). Complete periods give class numbers
     (KNOWN).
6. **The ladder's real need is factoring with known high bits** (§6, PROVED reformulation).
   - A polylog L_j oracle locates p from its top (1/(j+2))·log₂N bits, i.e. 2/(j+2) of p's bits.
     Read this way, j = 1 is Rivest–Shamir (1/3·log₂N bits) and j = 2 is Coppersmith (1/4·log₂N).
   - j = 2 is therefore already solved for the ladder's purpose by lattice reduction: Coppersmith /
     Howgrave-Graham blocks contain the j = 2 ladder blocks for every x ≤ √N, since
     (β − ½)² ≥ 0. This gives the KNOWN N^{1/4+o(1)}; E5 runs it (12/12).
   - j ≥ 3 would need fewer than (1/4)log₂N known high bits, beyond Coppersmith's bound. That is a
     standing open problem since 1996, not something a floor-sum trick is likely to deliver.
7. **A polylog L_2 oracle would also compute D(N) = Σ_{n≤N} τ(n) in Õ(N^{1/4})** (Prop D, PROVED).
   The best factorisation-free algorithms are N^{1/3} (Sladkey; Tao–Croot–Helfgott); Li's N^{9/28}
   claim uses factorisation of nearby integers. So the oracle is a stronger primitive than the
   factoring application needs. Prop D is implemented and tested for j = 1–4.

**Best open question.** Whether the ladder's actual queries (two fixed shifts, dyadic subintervals,
one N) admit an exact answer below H^{1/3} using N's arithmetic rather than the arc's shape. For
j = 2 the honest answer is "use Coppersmith instead". For j ≥ 3 the question is equivalent in
strength to beating Coppersmith's 1/4-bits barrier on the ladder's intervals. We see no mechanism.

## 1. Setting and two corrections to the premises

- Theorem L_j (notes/analytic_results.md §3): B = ⌈N^{1/(j+2)}⌉, blocks [x, x+H), H = ⌊cx/B⌋,
  P_j(h) = Σ_{i≤j}(−1)^i N h^i/x^{i+1}, and remainder |R| ≤ ρ = NH^{j+1}/x^{j+2}. The ladder needs
  F(P + s; h₁, h₂) only for s ∈ {0, −ρ−η} (even j) or {0, ρ} (odd j) and for bisection subintervals.
- Coefficient sizes (codex floorsum_1, PROVED; our task statement was imprecise here):
  β_i H^i ≍ c^i N/(xB^i), with β_i = N/x^{i+1}. At x ≍ √N:
  - the top term has β_jH^j ≍ c^j N^{(2−j)/(2(j+2))}: O(1) only for j = 2, and it shrinks for j > 2;
  - the curvature excursion is P″H² ≍ N^{(j−2)/(2(j+2))}, e.g. N^{1/6} for j = 4.
  - The curvature is P″(h) = (2N/x³)(1 + O(jh/x)), so on a block it is constant to within a factor
    1 + O(j/B).
- Candidates per block are O(1) only **on average**: Σ over blocks is O(B log N) (Shiu), but a single
  block can hold Ω(√H) of them (codex: N = 4x² + 1, j = 2, where {P(h)} ≤ ρ for h ≲ x^{1/4}). Per-block
  statements below say so where it matters.

## 2. The convexity trap, verified

**Scale (PROVED).** For j ≥ 2 and a ladder block, P″ ≍ N/x³, so
Hβ^{1/3} ≍ (cx/B)(N/x³)^{1/3} = cN^{1/3}/B. The blocks in a dyadic range number ≍ B/c, so the sum
over blocks of the convexity scale is ≍ N^{1/3} per dyadic range and Θ(N^{1/3} log N) in total,
**independently of j**. Raising j shrinks the number of blocks by exactly the factor by which each
block's convexity scale grows.

**Lemma V (PROVED; even j).** Let G = {(d, m) ∈ Z² : d, m ≥ 1, dm ≥ N}. For a block at x let
S_b = {(h, m) ∈ Z² : 0 ≤ h < H, m > P_j(h)} (the lattice points of the open epigraph, shifted to
d = x + h), and let V_b be the number of vertices of conv S_b, i.e. of the lower convex hull of
{(h, ⌊P_j(h)⌋ + 1)}. Every vertex of conv G with d ∈ [x, x+H) that is not at a candidate h is a
vertex of conv S_b. Hence Σ_b V_b ≥ V_glob([x₀, √N]) − #candidates.

*Proof.* For even j, R ≤ 0, so P_j ≥ f = N/(x+h), and S_b ⊆ {m > f} ⊆ G. A vertex v = (d, ⌈N/d⌉) of
conv G lies in S_b unless ⌈f(h)⌉ ≤ P_j(h). That needs an integer in [f(h), P_j(h)], so
{P_j(h)} ≤ P_j − f = |R| ≤ ρ: h is a candidate. (A divisor, f ∈ Z, is the case of equality.) An
extreme point of conv G that lies in S_b ⊆ G is extreme in conv S_b. ∎

E1 measured the closed version, the lower hull of {(h, ⌈P_j(h)⌉)}. The two point sets coincide
unless P_j(h) ∈ Z for some h; among the 18.2 million block points of the E1 instances up to 40 bits
there were none.

For odd j, S_b ⊇ G-block and the extra points sit at candidates. Adding a point can remove hull
vertices, so we have no proof; E1 shows no loss (V_tot/V_glob ≥ 1.11 for j = 3).

With the KNOWN count V_glob = Θ(N^{1/3} log N) (Balog–Bárány, arXiv 2602.06897, Thm 1.1; earlier
bounds by Alcántara–Blanco–Criado–Santos, arXiv 2501.19193) and #candidates = O(B log N)
(Theorem L_j step 7), Σ_b V_b = Ω(N^{1/3} log N) for every even j ≥ 2.

**Theorem T (PROVED; convex-access model).** Let R be a convex region and L = R ∩ Z² finite. Any
algorithm that accesses R only through (i) lattice-membership queries "is z ∈ R?" for z ∈ Z², and
(ii) lattice linear-optimisation queries "argmax_{z ∈ L} ⟨w, z⟩", and outputs |L| correctly for every
convex R, makes ≥ V(conv L) − 1 queries on R.

*Proof.* For each vertex v of conv L, R_v = conv(L ∖ {v}) is convex and R_v ∩ Z² = L ∖ {v} (it lies
in conv L ⊆ R and excludes the vertex v). Its count differs from |L| by one. A membership query
distinguishes R from R_v only if it asks z = v. An optimisation query distinguishes them only if its
answer for R is v. Each query thus rules out at most one alternative R_v. ∎

Applied to the epigraph region of a convex block arc (P convex for j ≥ 2), F(P; 0, H) needs
≥ V_b − 1 such queries (the open epigraph intersected with a box is convex). The ladder evaluates F(P; 0, H) on every block, so a convex-access oracle
costs Ω(N^{1/3} log N) in total for every even j. **Scope:** R_v is not a polynomial-arc region. An
algorithm that evaluates P off the lattice, or uses its coefficients, is outside the model; that is
precisely where an algebraic oracle would have to live.

**Prop PL (PROVED; piecewise-linear-floor model).** Suppose [0, H) is split into intervals I_1, …, I_r
such that on each I_k, ⌊P(h)⌋ = ⌊(p_k h + r_k)/q_k⌋ for integers p_k, q_k ≥ 1, r_k. Then
r ≥ V_b / O(log(H + max q_k)), with V_b as in Lemma V (the lower hull of {(h, ⌊P(h)⌋)}, shifted up
by one).

*Proof.* A hull vertex lying in I_k is a vertex of the lower hull of the points of I_k. Those points
are (u, ⌊ℓ(u)⌋) with ℓ(u) = (p_k u + r_k)/q_k. Since ⌊ℓ(u)⌋ ≥ ℓ(u) − 1 + 1/q_k, they are exactly the
lowest lattice points of the rational polygon {u ∈ I_k, m ≥ ℓ(u) − 1 + 1/q_k}, so their lower hull is
part of that polygon's integer hull. A planar rational polygon with O(1) edges has an integer hull
with O(φ) vertices for input size φ (Cook–Hartmann–Kannan–McDiarmid, Combinatorica 12 (1992)
27–37), here O(log(H + q_k) + bit size of the coefficients). ∎

Every method we know that is exact and "linearises" the arc has this form: hull walks, Hiary-type
continued-fraction blocking of the floor function, the codex cell construction, and our Farey
oracle. So all of them pay Ω(N^{1/3} log N / log N) = Ω(N^{1/3}) over the ladder.

**E1 (pre-registered; every prediction held).** Brute-force block hulls, 5 balanced semiprimes per size.

| j | V_tot / V_glob (min–max) | slope of ln(V_tot/ln N) vs ln N | per-block V / (cN^{1/3}·2^{1/3}/B) at 44 bits |
|---|---|---|---|
| 2 | 1.50–1.96 | 0.309 | 0.72 |
| 3 | 1.11–1.33 | 0.310 | 0.45 |
| 4 | 1.04–1.17 | 0.314 | 0.38 |

At 44 bits: V_tot = 160 770 (j = 2), 121 346 (j = 3), 114 320 (j = 4) against V_glob ≈ 107 000–110 000
and 28 546, 7 570, 3 070 blocks. Raising j from 2 to 4 removes 89% of the blocks and only 29% of the
vertices: the excess is the per-block boundary overhead, and the global hull is the floor. P1.1
(Lemma V check) held on all 50 even-j instances.

## 3. The exact oracle: Farey dissection

**Lemma R (PROVED).** If q ≥ 1, p, r are integers and r ≤ qP(u) − pu < r + 1, then
⌊P(u)⌋ = ⌊(pu + r)/q⌋.

*Proof.* Write qP(u) = n + θ with n = pu + r and 0 ≤ θ < 1, and n = q⌊n/q⌋ + s with 0 ≤ s ≤ q − 1.
Then s + θ < q, so ⌊(n + θ)/q⌋ = ⌊n/q⌋. ∎

**Algorithm O** (`FloorSumConvex`). Input: P = num/den with num ∈ Z[h], P″ of constant sign on
[h_lo, h_hi). Set β̂ = max|second difference|/2 at three points and Q = ⌈β̂^{−1/3}⌉ (capped at H + 1).
Starting at h:
1. take p/q, the last continued-fraction convergent with q ≤ Q of the local slope s = P(h+1) − P(h);
2. take the window [h, h + L) with L = max(1, ⌊1/(qQβ̂)⌋);
3. on the window g(u) = qP(u) − pu is convex (or concave). Split it at its integer extremum, found by
   binary search on the sign of g(u+1) − g(u). On each monotone side G = ⌊g⌋ is monotone; find its
   constant runs by binary search;
4. on each run with G ≡ r, add floor_sum_linear for ⌊(pu + r)/q⌋ (Euclid, O(log)).

Runs are stored with prefix sums, so subinterval queries cost O(log) more, and each (block, shift)
is built once. Shifts are absorbed into the constant coefficient.

**Theorem O (PROVED).** The output is exact for every P with P″ of constant sign. If moreover
max|P″|/2 ≤ κβ̂ on the range, then #windows ≤ 8Hβ̂^{1/3} + 1 and
#runs ≤ 2·#windows + 2(1 + 2κ)Hβ̂^{1/3}, and the cost is O((1 + κ)(1 + Hβ̂^{1/3}) log H) exact evaluations of P on O(j log N + log H)-bit numbers.

*Proof.*
- Exactness: Lemma R on every run; the runs tile the range, because G is monotone on each side of
  the split.
- Windows: L = max(1, ⌊1/(qQβ̂)⌋) ≥ max(1, ⌊1/(Q²β̂)⌋). If β̂ ≤ 1/8 then Q ≤ β̂^{−1/3} + 2 ≤ 2β̂^{−1/3},
  so 1/(Q²β̂) ≥ β̂^{−1/3}/4, and #windows ≤ 8Hβ̂^{1/3} + 1. If β̂ > 1/8 then trivially
  #windows ≤ H ≤ 2Hβ̂^{1/3}. (If Q is capped at H + 1 there is a single window.)
- Runs: on a window, |g′(u)| ≤ |qs − p| + q·max|P″|·L < 1/Q + 2qκβ̂L, since the last convergent with
  q ≤ Q satisfies |qs − p| < 1/Q. The variation of g over the window is therefore
  ≤ L/Q + 2qκβ̂L² ≤ (1 + 2κ)L/Q, using L ≤ 1/(qQβ̂) when L > 1 (a window with L = 1 has one run).
  Each of the two monotone sides has at most (its variation + 1) distinct values of ⌊g⌋, so a window
  has ≤ 2 + 2·variation runs. Summing, #runs ≤ 2·#windows + 2(1 + 2κ)H/Q, and H/Q ≤ Hβ̂^{1/3}.
- Each window and run costs O(log L) evaluations. ∎

For ladder blocks κ = 1 + O(j/B). So T_j(H) = O(Hβ^{1/3} log H) with β = N/x³, which at x ≍ √N is
H^{1/3} (j = 2), H^{4/9} (j = 3), H^{1/2} (j = 4), and in general H^{2(j−1)/(3j)}. Codex floorsum_1
proved the same bound with a fixed-cell variant (cells with ML³ ≤ 1/64, one CF approximation per cell,
and one congruence-class correction).

The bound is an upper bound. On major arcs (slope within 1/(qQ) of p/q with q small) windows are
far longer than Q, and the cost drops. Codex's example P = h + γh²/H² costs O(1).

**E2 (pre-registered).**
- P2.1 (exactness): 0 mismatches in 3 000 random convex polynomials of degree 2–4, and in 5 400
  ladder-block range queries (j = 2, 3, 4; 40, 48, 56 bits; three shifts). ✓
- **P2.2 and P2.3 FAILED as specified, by a design error.** The scaling blocks were placed at
  x = ⌊√N/2⌋ + O(10³). There the slope N/x² is within O(10³/√N) of the integer 4: a major arc for the
  whole block. The cost ratio (windows + runs)/(Hβ^{1/3}) was 0.003–0.39 for j = 2 (slope −0.002 in
  ln H) and 0.08–0.48 for j = 3 (slope 0.316). Only j = 4 behaved as predicted (1.00–1.14, slope 0.493),
  because there the block's slope range ≍ 4/B already exceeds the major-arc width.
- **E2b (post-hoc, disclosed).** The same sizes with x uniform in [√N/4, √N], 8 blocks per size.
  Ratio 0.51–2.88 (median 1.93) for j = 2, 1.82–2.14 for j = 3, 1.49–2.05 for j = 4. Slopes against
  ln H: **0.337, 0.442, 0.506**, against 1/3, 4/9, 1/2. Experimentally supported for generic blocks.

## 4. Fourier and smoothing

**What Hiary computes (KNOWN; codex floorsum_1 checked the theorems).**
- Quadratic: Hiary, Ann. of Math. 174 (2011) 859–889, Thm 1.1. L^{−r}Σ_{n≤L} n^r e(an + bn²) for any
  a, b, to absolute error Aν^κε in ν^{O(1)} operations, ν = (r+1)log(L/ε), no precomputation.
- Cubic: Hiary, Ann. of Math. 174 (2011) 891–946, Thm 2.2. Cubic coefficient 0 ≤ d ≤ L^{μ−3}, μ ≤ 1,
  after L^{4μ}ν^{O(1)} precomputation. For the j = 3 ladder at x ≍ √N the cubic coefficient of
  kP is ≍ kH^{−10/3}, so frequencies k ≲ H^{1/3} have μ ≈ 0.
- So each Weyl sum Σ_h e(kP(h)) costs polylog to inverse-polynomial accuracy. The Weyl sums are not
  the bottleneck.

**Identity.** With ψ(y) = y − ⌊y⌋ − ½ and the midpoint convention ψ₀ (ψ₀(n) = 0),
⌊y⌋ = y − ½ − ψ₀(y) + ½·1_Z(y). The endpoint term must be kept (codex).

**Prop F (PROVED).** Let T be a real trigonometric polynomial of degree K and
G = {y ∈ R/Z : |ψ(y) − T(y)| < 1/4}.
1. G^c is a union of at most 2K + 2 intervals.
2. meas(G^c) ≥ 1/(51K).

*Proof.*
1. On (0, 1), ψ − T ∓ 1/4 = y − ½ − T(y) ∓ ¼. Its derivative 1 − T′ has ≤ 2K zeros, so each function
   has ≤ 2K + 1 zeros. With the jump at 0 this gives ≤ 2K + 2 intervals.
2. Let μ = meas(G^c). On G, |T| < 3/4. If μ ≤ 1/(16πK), the Remez inequality for trigonometric
   polynomials (Erdélyi, J. London Math. Soc. 45 (1992); in the normalisation θ = 2πy:
   ‖T‖_∞ ≤ e^{4K·2πμ} sup_G|T|) gives ‖T‖_∞ ≤ (3/4)e^{1/2} < 1.24.
   - Bernstein's inequality then gives |T′| ≤ 2π·1.24·K < 7.8K.
   - Let δ = 1/(4 + 31.2K). If T(0) ≥ 0, then on (0, δ), T ≥ T(0) − 7.8Kδ ≥ −7.8Kδ while
     ψ < −½ + δ, so |ψ − T| ≥ ½ − δ − 7.8Kδ = ¼. If T(0) < 0, the same holds on (−δ, 0), where
     ψ > ½ − δ and T < 7.8Kδ.
   - So μ ≥ δ ≥ 1/(36K) for K ≥ 1.
   In both cases μ ≥ min(1/(16πK), 1/(36K)) ≥ 1/(51K). ∎

**Corollary F (PROVED; certified smoothing model).** Suppose an algorithm computes
Σ_h ψ(P(h)) = Σ_h T(P(h)) + Σ_{h∈E}(ψ − T)(P(h)) + Err. It uses K frequencies (cost ≥ K), corrects
the explicitly enumerated set E (cost ≥ |E|), and certifies |Err| ≤ Σ_{h∉E}|ψ − T|(P(h)) < ½. Then at
most one h ∉ E lands in G^c. If the points {P(h)} (h in the block) have maximal circular gap g, an
interval of length ℓ contains ≥ ℓ/g − 1 of them, so |E| ≥ 1/(51Kg) − (2K + 2) − 1. The cost
C ≥ max(K, K + |E|) ≥ max(K, 1/(51Kg) − K − 3), and minimising over K (at K = (102g)^{−1/2}) gives
C ≥ 0.099/√g − 3.

*Scope.* It does not cover uncertified schemes that rely on cancellation, aggregated corrections of
resonant clusters (codex), or non-trigonometric smoothing. Finding E is itself a near-integer
location problem of the same kind; the bound grants it for free.

Codex floorsum_1 gave two companion bounds: Ω(√H) for the worst shift (averaging the constant term
makes E_s|{h : ‖P(h) + s‖ < 1/K}| = 2H/K), and Ω(H^{2/3}) for the explicit family
P = h + γh²/H².

**E4 (pre-registered; held).** 200 random j = 2 ladder blocks per size, x uniform in [√N/4, √N].

| bits | mean H | n_δ / (2δH) at δ = 1/⌈√H⌉ | min over K = 2^i of (K + n_{1/K}) / √H |
|---|---|---|---|
| 48 | 1 166 | 0.998 | 2.88 |
| 56 | 4 762 | 1.003 | 2.89 |
| 64 | 19 502 | 0.999 | 2.88 |

**E4b (post-hoc diagnostic).** Median g·H = 5.2, 6.5, 7.9 (≈ ln H), so Corollary F gives
Ω(√(H/log H)) asymptotically. Its explicit constants are vacuous at these sizes: the median bound is
≤ 0.014√H. The measured model cost (≈ 2.9√H) is the realistic figure.

**Consequence.** Per block ≳ √H, so the total is ≥ B√H ≍ N^{1/(j+2)}·N^{(1/2 − 1/(j+2))/2}
= N^{1/4 + 1/(2(j+2))}. That is N^{3/8} at j = 2 (worse than the hull), → N^{1/4} as j → ∞, and never
below N^{1/4}. **Fourier with exact correction is refuted for this model** as a route below N^{1/4}.
Smoothing cannot make K polylog.

**Self-reference (PROVED, elementary).** For the ladder's actual f = N/d, the smoothing correction
set at scale δ is {d : ‖N/d‖ < δ}, i.e. the divisors of N − k with |k| < δd, and it contains the
divisors of N themselves. Li's de-smoothing for D(N) (below) corrects by factoring integers near
N. For restricted divisor sums that would include N itself: circular.

## 5. Euclid for quadratics, reciprocity, and known exact counters

**Axis swap (PROVED).** For Q(h) = (ah² + bh + c)/m increasing and nonnegative on [0, H),
Σ_{h<H}⌊Q(h)⌋ = Σ_{k=1}^{⌊Q(H−1)⌋}(H − t(k)), with t(k) = ⌈(−b + √(b² + 4a(mk − c)))/(2a)⌉ clipped to
[0, H] (codex). Completing the square, the count is #{(y, t) : t² ≤ 4amy + Δ, t ≡ b (mod 2a)}.
Counting this by t gives back the original sum, so the double swap is the identity. Unlike the
linear case, no modulus shrinks: the inverse of a parabola is a square-root curve, not a parabola
with smaller parameters.

**Prop A (PROVED).** The affine arc length λ(Γ) = ∫_Γ κ^{1/3}ds is invariant under equi-affine maps
(det = 1, plus translations), hence under SL₂(Z) ⋉ Z². For a graph y = P(h), λ = ∫|P″(h)|^{1/3}dh.
- One step of linear Euclid (reduce the slope mod 1, swap axes, reflect) is such a map. It preserves
  the lattice and the boundary class.
- A line has λ = 0. A ladder block arc has λ ≍ H(2β)^{1/3} ≍ cN^{1/3}/B.
- So no sequence of lattice moves can bring it to a line. Convex lattice regions with affine
  perimeter λ have ≍ λ hull vertices in the generic case (Balog–Bárány), which is the trap again.
- This explains why the theta transformation exists for e(P(h)) but not for ⌊P(h)⌋. Poisson
  summation is not a lattice automorphism of the point set: it trades the length H for 2βH dual
  terms, and it needs a smooth summand.

**Unrestricted hardness (KNOWN + PROVED reduction).** Manders–Adleman (JCSS 16 (1978) 168–184):
deciding ∃ 0 ≤ h < X with h² ≡ A (mod m) is NP-complete. Since
⌊(h² − A)/m⌋ − ⌊(h² − A − 1)/m⌋ = 1[m | h² − A], two quadratic floor sums decide it. So polylog
quadratic floor sums with arbitrary coefficients and intervals imply P = NP. Codex floorsum_1 and
floorhard_1 found this independently. Floorhard_1 goes further and claims an embedding into the
Taylor family itself, by enlarging N; that is floorhard's result to state, and we only note that it
is consistent with everything here. Either way, the uniform oracle hypothesis of Theorem L_j is far
stronger than what factoring needs (§6).

**Complete periods (KNOWN).** Σ_{x<m}⌊ax²/m⌋ involves h(−m) for primes m ≡ 3 (mod 4) (Dirichlet; see
analytic §6, row Q). Codex gave the elementary form
Σ_{h<p}⌊h²/p⌋ = (p − 1)(p − 2)/3 − (1/p)Σ_{r<p} r(r/p).

**Known exact counters (KNOWN; checked by codex, to be re-checked by the literature critic).**
- Har-Peled, Comput. Geom. 10 (1998) 125–138: output-sensitive discrete convex hulls; a radius-R
  disc in O(R^{2/3} log R), so the exact circle count follows via Pick's theorem.
- Sladkey, arXiv 1206.3369 (2012): exact D(N) in O(N^{1/3}) without factoring, by successive
  approximation of the hyperbola; Tao–Croot–Helfgott, Math. Comp. 81 (2012), Prop. 2.2, the same
  order.
- Li, arXiv 2609.29949 (Sept 2026; claim not audited): exact D(N) in N^{9/28+o(1)}, using
  Lenstra–Pomerance factorisation of the integers near N.
- We found no factorisation-free exact D(N) or lattice-points-under-a-parabola algorithm below the
  hull complexity, and no reciprocity law for incomplete quadratic floor sums.

**Amortisation (HEURISTIC; card F6).** Consecutive Farey windows form Stern–Brocot chains, and a
"universal Euclidean" aggregation could sum a chain in O(log) time. Generic chains have O(1) length
(they are long only near major arcs), so this gains a constant, not an exponent. Sladkey's algorithm
already exploits it.

## 6. Arithmetic structure: what the ladder really needs

**Reformulation (PROVED, elementary; also codex).** At x ≍ √N, an L_j block has length
H ≍ N^{1/2−1/(j+2)}. Knowing p to within H means knowing its top (1/(j+2))log₂N bits, a fraction
2/(j+2) of p's bits. The ladder is "guess the top bits (B guesses per dyadic range), then solve one
block".

| j | block at √N | known high bits / log₂N | ladder cost with a polylog oracle | known block solver |
|---|---|---|---|---|
| 1 | N^{1/6} | 1/3 | Õ(N^{1/3}) | Euclid; Hiary; Rivest–Shamir (EUROCRYPT '85) |
| 2 | N^{1/4} | 1/4 | Õ(N^{1/4}) | Coppersmith (J. Cryptology 10 (1997)); Howgrave-Graham (CaLC 2001) |
| 3 | N^{3/10} | 1/5 | Õ(N^{1/5}) | none known |
| 4 | N^{1/3} | 1/6 | Õ(N^{1/6}) | none known |

*Caveat (codex).* A floor oracle locates a block's candidates in O((1 + C_b) polylog) calls. A
per-interval polylog guarantee from known high bits therefore also needs a per-block candidate bound,
and single blocks can hold Ω(√H) candidates (§1).

**Lattice containment (PROVED from the KNOWN theorem).** Coppersmith/Howgrave-Graham finds the
divisors p ≥ N^α of N with |p − x| ≤ N^{α²} in polynomial time (Cohn–Heninger, arXiv 1008.1284,
Thm 1.1, d = 1 form). The j = 2 block at x = N^α has H ≍ N^{α−1/4} ≤ N^{α²}, because
α² − α + 1/4 = (α − ½)² ≥ 0. So lattice reduction is a polylog block oracle for the j = 2 ladder's
real need, which is locating divisors, not floor sums. That gives the KNOWN deterministic N^{1/4+o(1)}.
For j = 3 the containment holds only for α ≤ (1 − √(1/5))/2 ≈ 0.276; between N^{0.276} and N^{1/2}
the j = 3 blocks are longer than Coppersmith's interval. The lattice method's own cost per dyadic
range is N^{α−α²}, maximal (N^{1/4}) at α = ½.

**E5 (Amendment A1; KNOWN algorithm).** `lattice_factor` with dimension ≤ 16 (m = 7, n = 16 chosen at
√N): 12/12 correct (36–48 bits). At √N its interval X is 0.68–0.78 of the ladder's j = 2 block,
inside the amended [0.2, 2]. The fixed-dimension loss N^{−Θ(1/n)} and the LLL factor 2^{(n−1)/4}
account for the gap; X ≥ H at √N needs n ≍ log N. Python cost: 48 s at 48 bits (52 000 LLL calls of
dimension 16). It is a demonstration, not a new claim.

**Prop D (PROVED; implemented as `divisor_summatory_ladder`, tested j = 1–4).** With an exact L_j
oracle, D(N) = Σ_{n≤N} τ(n) = 2Σ_{d≤√N}⌊N/d⌋ − ⌊√N⌋² is computable with O(B log² N) oracle calls
plus O(B log N) further operations.

*Proof.* Per block, Σ_h⌊N/(x+h)⌋ = F(P_j; 0, H) − Σ_{h∈C}(⌊P_j(h)⌋ − ⌊N/(x+h)⌋).
- For even j, 0 ≤ P − f ≤ ρ < 1, so ⌊P⌋ − ⌊f⌋ ∈ {0, 1}, and it is 1 only if {P} < ρ, i.e. h ∈ C.
- For odd j, f − P ∈ [0, ρ], and ⌊f⌋ ≠ ⌊P⌋ only if {P} ≥ 1 − ρ, so again h ∈ C.
- C is located by bisection (Theorem L_j step 6). The d < 2B are summed directly.
- The total candidate count is O(B log N) (Theorem L_j step 7). ∎

So a polylog L_2 oracle gives D(N) in Õ(N^{1/4}), below every known factorisation-free algorithm
(N^{1/3}) and below Li's factorisation-assisted N^{9/28}. A polylog L_3 oracle would give Õ(N^{1/5}).
This is a second, independent reason to expect the oracle to be hard. It is evidence, not a
reduction from a hard problem.

## 7. End-to-end: the operational ladder (E3, pre-registered; all predictions held)

`ladder_factor_op(N, j)` is analytic's ladder with the simulated oracle replaced by Algorithm O.
Primes of 28–44 bits (3 per size) give the full-range cost; balanced semiprimes (5 per size) check
correctness.

| bits | j = 2 evals | j = 2 blocks | j = 2 windows+runs | j = 3 evals | j = 3 blocks | wall (j = 2) |
|---|---|---|---|---|---|---|
| 28 | 58 953 | 1 164 | 11 085 | 50 115 | 532 | 0.04 s |
| 32 | 150 497 | 2 459 | 27 889 | 129 402 | 996 | 0.11 s |
| 36 | 427 421 | 5 595 | 77 308 | 371 677 | 1 968 | 0.30 s |
| 40 | 1 236 638 | 12 862 | 218 550 | 1 099 327 | 3 922 | 0.85 s |
| 44 | 3 265 620 | 27 211 | 561 917 | 2 952 360 | 7 277 | 2.2 s |

- P3.1: 80/80 correct. ✓
- P3.2: the slope of ln(evals/ln²N) against ln N is 0.290 (j = 2) and 0.296 (j = 3), inside
  [0.28, 0.40]. ✓ (The local exponent sits a little below 1/3 because the log factors are not
  exactly ln²N.)
- P3.3: evals(j = 3)/evals(j = 2) = 0.85–0.90. ✓ Going from j = 2 to j = 3 cuts the blocks by 3.7×
  but leaves the work unchanged.

**End-to-end factoring exponent with the best oracle we have: N^{1/3+o(1)}, for every j.** With
oracles we cannot build: N^{1/(j+2)} (Theorem L_j); and for j = 2 the lattice route already gives
the KNOWN N^{1/4+o(1)}.

## 8. Pre-registration outcomes

- E1: P1.1 ✓ (50/50), P1.2 ✓ (1.04–1.96), P1.3 ✓ (0.309, 0.310, 0.314), P1.4 ✓ (0.72, 0.45, 0.38).
  Decision: C1 supported.
- E2: P2.1 ✓. **P2.2, P2.3 ✗** (major-arc design error, §3). Decision rule NOT MET. Post-hoc E2b
  supports C2 for generic blocks; it is labelled post-hoc.
- E3: P3.1 ✓, P3.2 ✓, P3.3 ✓. Decision: C3 supported.
- E4: P4.1 ✓ (0.998–1.003), P4.2 ✓ (2.88–2.89). Decision: C4's premise supported.
- E5 (amended): P5.1′ ✓ 12/12, P5.2′ ✓ (0.68–0.78).

## 9. Barrier summary (this agent)

| # | Statement | Scope | Status |
|---|---|---|---|
| T | Exact block floor sum needs ≥ V_b − 1 membership/LP queries; Σ_b V_b = Ω(N^{1/3} log N) for even j | convex-access algorithms | PROVED (Lemma V + Balog–Bárány) |
| PL | Piecewise-linear-floor decompositions need ≥ V_b/O(log H) pieces | hull walks, CF blocking, Farey, cells | PROVED |
| A | Lattice moves preserve affine length; Euclid cannot straighten an arc | equi-affine reductions | PROVED |
| F | Certified smoothing: cost ≥ 0.099/√g − 3 per block, ≍ √H for ladder blocks | trigonometric approximants + individual correction | PROVED (+ E4 for g) |
| MA | Unrestricted quadratic floor sums ⇒ P = NP | arbitrary coefficients | KNOWN + 2-call reduction |
| D | Polylog L_2 ⇒ D(N) in Õ(N^{1/4}) | — | PROVED (implication; evidence of strength) |
| HB | Polylog L_j ⇒ factoring from (1/(j+2))log₂N high bits (with a per-block candidate bound) | — | PROVED (reformulation) |

## 10. Reviews and fixes

- GPT-6.1-sol blind round (floorsum_1): reached the same verdict independently (Õ(H^{1/3})
  CF-linearisation oracle, N^{1/3} total for every j). It contributed the coefficient-size correction,
  the per-block candidate counterexample, the Manders–Adleman reduction, the shift-averaging Fourier
  bound and the high-bits table caveat. All are incorporated above.
- GPT-6.1-sol adversarial round: pending.
- Opus correctness/cost critic: pending.
- Opus literature critic: pending.

## 11. Reproduction

```bash
.venv/bin/python -m pytest -q tests/test_floorsum.py
PYTHONPATH=. .venv/bin/python -m src.floorsum.run_experiments e1     # ~4 min
PYTHONPATH=. .venv/bin/python -m src.floorsum.run_experiments e2     # ~1 min
PYTHONPATH=. .venv/bin/python -m src.floorsum.run_experiments e2b    # post-hoc, ~2 min
PYTHONPATH=. .venv/bin/python -m src.floorsum.run_experiments e3     # ~1 min
PYTHONPATH=. .venv/bin/python -m src.floorsum.run_experiments e4     # ~1 min
PYTHONPATH=. .venv/bin/python -m src.floorsum.run_experiments e4b    # post-hoc, ~1 min
PYTHONPATH=. .venv/bin/python -m src.floorsum.run_experiments e5     # ~6 min
PYTHONPATH=. .venv/bin/python -m src.floorsum.run_experiments summary
```
