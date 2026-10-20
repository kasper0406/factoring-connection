# floorsum: fast exact floor sums for the ladder (algorithm side)

Agent `floorsum`, 2026-10-20. Final version. It incorporates the GPT-6.1-sol adversarial round
(floorsum_2: ACCEPT WITH FIXES, major revision), two read-only Opus critics (correctness/cost:
ACCEPT WITH FIXES; literature: OK WITH FIXES), and the orchestrator's top-block request (E6). §10
lists the fixes.

**Files.**
- Pre-registration: `notes/floorsum_prereg.md`. Commit `a1a98a4` was made before the formal E1–E5
  runs; development checks were run before it and are disclosed there. Amendment A1 (E5, with a
  disclosed peek) is `6d5b3ab`; Amendment A2 (E6) is `d40f4a0`, before E6.
- Code:
  - `src/floorsum/farey.py`: the exact oracle;
  - `src/floorsum/ladder_ops.py`: the operational ladder and D(N) via the ladder;
  - `src/floorsum/lattice_block.py`: the Coppersmith divisor block oracle (KNOWN);
  - `src/floorsum/lattice_window.py`: the bivariate lattice for window queries (E6);
  - `src/floorsum/run_experiments.py`.
- Tests: `tests/test_floorsum.py` (19 pass).
- Raw data: `experiments/raw/floorsum_e{1,1int,2,2b,3,4,4b,5,6,6b,6c}.jsonl`. Summary:
  `experiments/analysis/floorsum_summary.json`. Logs: `experiments/logs/floorsum_*.log`.
- Cards: `notes/floorsum_cards.md`.
- GPT-6.1-sol: `codex/floorsum_1.{md,out}` (blind) and `codex/floorsum_2.{md,out}` (adversarial).

**Labels.** KNOWN (with reference), PROVED (proof here), experimentally supported, HEURISTIC,
CONJECTURED, REFUTED for this model, OPEN. "No prior formulation found by us" is the strongest
priority claim made.

## 0. Verdict (plain language)

**No fast exact floor-sum oracle was found.** Our best exact oracle costs about H·β^{1/3} per block,
for every degree j. Plugged into the ladder, that gives N^{1/3} for every j: exactly the trap. So
nothing here beats Hiary's N^{1/3}, let alone N^{1/4} (Pollard–Strassen, Coppersmith) or Harvey's
N^{1/5}. This is a statement about the algorithms we could build and the models we could analyse. It
is not a lower bound for factoring.

**Top blocks with lattices (the orchestrator's question): no.** On top blocks, lattice small-root
methods answer the ladder's positive-width window queries only on sub-blocks of length ≈ xN^{−1/3}
(E6b: measured/predicted = 0.69–1.75 for j = 2, 3). That is Hiary's length scale, whatever j is, so
the total is again N^{1/3}.

What we established:

1. **The trap is real. It is proved in three access models, each with a stated scope** (§2).
   - Per block the convexity scale is Hβ^{1/3} ≍ cN^{1/3}/B. Over the O(B log N) blocks that is
     Θ(N^{1/3} log N) for every j ≥ 2.
   - **Lemma V (PROVED, even j).** Every vertex of the integer hull of {xy ≥ N} is a vertex of some
     block hull, except at candidates.
   - **Theorem T (PROVED).** Deterministic algorithms that see the block's epigraphs only through
     lattice-membership / lattice-LP queries need ≥ V_b queries per block. This holds for each
     individual floor sum, and also for the ladder's difference F(P) − F(P − ρ − η) when that is
     computed through the two epigraphs.
   - **Prop PL (PROVED).** Interval decompositions into linear-floor pieces need ≥ V_b/O(log(H + q))
     pieces.
   - So convex-access and linearising oracles pay Ω(N^{1/3}) in total for every even j.
   - E1 (pre-registered, all held): summed block-hull vertex counts grow like N^{0.31} for j = 2, 3, 4.
2. **The best exact oracle is the Farey dissection** (§3). It is operational for any degree, with shifts
   and subintervals.
   - Exactness rests on a one-line lemma (Lemma R).
   - Cost (Theorem O, PROVED for blocks where P″ varies by a bounded factor, which includes all
     ladder blocks): O((1 + Hβ^{1/3}) log H) evaluations.
   - At x ≍ √N this is T_j ≈ H^{2(j−1)/(3j)}: H^{1/3}, H^{4/9}, H^{1/2} for j = 2, 3, 4. The targets were
     H^{1/3−ε} (j = 2) and H^{1/10} (j = 4).
   - The pre-registered scaling test (E2) failed by a design error: its blocks sat on a major arc. The
     post-hoc, exploratory E2b on generic blocks measured slopes 0.337, 0.442, 0.506.
3. **End to end** (§7, E3). The ladder is fully operational for j = 2, 3, with no simulated oracle and
   80/80 correct.
   - Evaluations grow like N^{0.29–0.30}·ln²N over 28–44 bits, i.e. N^{1/3} up to logs.
   - j = 3 is only 10–15% cheaper than j = 2, with the same exponent.
4. **Fourier smoothing does not help** (§4).
   - **Prop F (PROVED).** A degree-K trigonometric approximant of the sawtooth is off by ≥ 1/4 on a set
     of measure ≥ 1/(51K) made of ≤ 2K + 2 intervals.
   - In the certified individual-correction model (all frequencies up to K charged), that forces a
     cost ≥ 0.099/√g − 3 per block, where g is the block's maximal gap.
   - For ladder blocks g ≈ (ln H)/H (post-hoc E4b, j = 2 only), so the cost is ≍ √H up to logs
     (experimentally supported). Hence the total is ≳ N^{1/4 + 1/(2(j+2))}: never below N^{1/4}, but
     below N^{1/3} for j ≥ 5, so this model alone does not rule out beating Hiary at large j.
5. **"Euclid for quadratics" fails in the obvious form** (§5).
   - The axis swap returns square-root floors, and the double swap is the identity (PROVED).
   - Lattice moves preserve affine arc length (PROVED, classical). This is a heuristic explanation of
     why Euclid-type moves cannot straighten an arc; it is not an impossibility theorem.
   - Unrestricted quadratic floor sums are NP-hard (Manders–Adleman), and #P-hard (floorhard Thm SS).
6. **The ladder's real need is factoring from known high bits** (§6).
   - A polylog L_j oracle, plus a per-block candidate bound, locates p from its top (1/(j+2))·log₂N bits.
   - j = 2 is already served by Coppersmith lattices. Their blocks contain the ladder's blocks because
     (α − ½)² ≥ 0 (E5: 12/12). This is the KNOWN N^{1/4+o(1)}.
   - j ≥ 3 would need fewer than (1/4)log₂N known high bits. No polynomial-time method for that is
     known (Coppersmith EUROCRYPT '96).
   - Window queries at top blocks are not lattice-easy either (item "top blocks" above; §6.4).
7. **A polylog L_2 oracle would compute D(N) = Σ_{n≤N} τ(n) in Õ(N^{1/4})** (Prop D, PROVED; implemented
   and tested for j = 1–4).
   - We found no factorisation-free exact D(N) algorithm below N^{1/3}. Li's N^{9/28} is randomized
     and uses factorisation of the integers near N.
   - So the oracle is a stronger primitive than factoring needs. This agrees with floorhard: the uniform
     oracle is NP-hard and #P-hard.

**Best open question.** Is there a deterministic, factorisation-free top-block oracle for the ladder's
**difference** queries (fixed N, shift −ρ−η, dyadic subintervals) below H^{1/3} (j = 2) or H^{1/10}
(j = 4)? Such an oracle would have to use N's arithmetic rather than the arc's shape, and could not
access the two epigraphs only as convex sets.
- For j = 2 the honest answer is "use Coppersmith".
- For j ≥ 3 a positive answer would beat Coppersmith's high-bits barrier on the ladder's intervals.
- Neither we, nor codex, nor floorhard see a mechanism.
- Note that with factorisation allowed the difference interface is easy (codex; §6.5). The problem is
  circular, because the k = 0 candidates are the divisors of N.

## 1. Setting and corrections to the premises

- **Theorem L_j** (notes/analytic_results.md §3).
  - B = ⌈N^{1/(j+2)}⌉ (the code uses ⌊·⌋; then ρ ≤ c^{j+2}(1 + 1/B)^{j+2} < 1 for c = ½, B ≥ 2).
  - Blocks [x, x+H) with H = ⌊cx/B⌋; P_j(h) = Σ_{i≤j}(−1)^i N h^i/x^{i+1}.
  - |R| ≤ ρ = NH^{j+1}/x^{j+2}.
  - The ladder needs only the differences F(P; I) − F(P − ρ − η; I) (even j) or F(P + ρ; I) − F(P; I)
    (odd j), on [0, H) and on bisection subintervals I. Prop D needs F itself.
- **Coefficient sizes** (codex floorsum_1, PROVED; our task statement was imprecise). With
  β_i = N/x^{i+1}, β_i H^i ≍ c^i N/(xB^i). At x ≍ √N:
  - the top term has β_j H^j ≍ c^j N^{(2−j)/(2(j+2))}: O(1) only for j = 2;
  - the curvature excursion is P″H² ≍ N^{(j−2)/(2(j+2))}, e.g. N^{1/6} for j = 4;
  - P″(h) = (2N/x³)(1 + O(jh/x)), constant on a block to a factor 1 + O(j/B).
- **Candidates per block are O(1) only on average.** The total over blocks is O(B log N) (Shiu). A
  single block can hold Ω(√H) candidates (codex: N = 4x² + 1, j = 2).

## 2. The convexity trap, verified

**Scale.** For j ≥ 2, P″ ≍ N/x³ on a block, so Hβ^{1/3} ≍ (cx/B)(N/x³)^{1/3} = cN^{1/3}/B. A dyadic range
has ≍ B/c blocks, so the convexity scale sums to ≍ N^{1/3} per dyadic range and Θ(N^{1/3} log N) in
total, **independently of j**. Raising j cuts the number of blocks by exactly the factor by which each
block's scale grows. That Hβ^{1/3} is the right per-arc scale is the classical Jarník (Math. Z. 24
(1926)) / Swinnerton-Dyer / Bombieri–Pila count of lattice points on convex arcs (KNOWN); the
computation for the ladder is done here.

**Lemma V (PROVED; even j).**
- Let G = {(d, m) ∈ Z² : d, m ≥ 1, dm ≥ N}.
- For a block at x, let S_b⁼ = {(h, m) : 0 ≤ h < H, m ≥ P_j(h)} (closed epigraph) and
  S_b = {(h, m) : 0 ≤ h < H, m > P_j(h)} (open epigraph), shifted to d = x + h.
- Then every vertex of conv G with d ∈ [x, x+H) that is not at a candidate h is a vertex of both
  conv S_b⁼ and conv S_b.
- Hence Σ_b V_b ≥ V_glob([x₀, √N]) − #candidates.

*Proof.*
1. For even j, R ≤ 0, so P_j ≥ f = N/(x+h). Then S_b ⊆ S_b⁼ ⊆ {m ≥ f} ⊆ G.
2. A vertex v = (d, ⌈N/d⌉) of conv G is missing from S_b⁼ only if ⌈f⌉ < P_j, and from S_b only if
   ⌈f⌉ ≤ P_j. Either way there is an integer in [f, P_j], so {P_j(h)} ≤ P_j − f = |R| ≤ ρ: h is a
   candidate. A divisor (f ∈ Z) is the equality case.
3. An extreme point of conv G that lies in a subset of G is extreme for that subset. ∎

*Restricted range.* Balog–Bárány (arXiv 2602.06897, Thm 1.1) prove Θ(N^{1/3} log N) vertices for the
whole hull. By the symmetry (x, y) ↔ (y, x), about half lie below √N. Removing d < x₀ = O(B) loses
O(B). So V_glob([x₀, √N]) = Θ(N^{1/3} log N); earlier bounds are Alcántara–Blanco–Criado–Santos
(arXiv 2501.19193). With #candidates = O(B log N) (Theorem L_j step 7), Σ_b V_b = Ω(N^{1/3} log N)
for every even j ≥ 2.

*Odd j.* S_b contains the block of G plus extra points at candidates. Adding points can remove hull
vertices, so there is no proof. At the block start the open epigraph also omits (x, N/x) when x | N;
the ladder's forced block-start division handles that. E1 shows no loss: V_tot/V_glob ≥ 1.11 for
j = 3.

*Integral values.* The two hulls coincide unless some P_j(h) is an integer. The post-hoc count
`e1int` found 0 such points among the 73.4 million block points of all E1 instances.

**Theorem T (PROVED; convex-access model).** Fix a finite truncating box Λ.
- Let R be convex and L = R ∩ Λ ∩ Z².
- Consider a deterministic (or zero-error) algorithm that accesses R only through:
  - (i) membership queries "z ∈ R?" for z ∈ Z²;
  - (ii) linear-optimisation queries returning **one** maximiser of ⟨w, z⟩ over L, under a fixed
    deletion-stable tie rule (e.g. lexicographic).
- If it outputs |L| correctly for every convex R, it makes ≥ V(conv L) queries on R.

*Proof.*
1. For each vertex v of conv L, R_v = conv(L ∖ {v}) is convex, lies in R, and has
   R_v ∩ Λ ∩ Z² = L ∖ {v}.
2. A membership query separates R from R_v only at z = v. An optimisation query separates them only if
   its answer for R is v; otherwise the same maximiser stays legal for R_v under the tie rule.
3. Each query therefore rules out at most one of the V alternatives, and all must be ruled out. ∎

*Difference version (critic).* Suppose an algorithm computes Δ = |L′| − |L| for nested convex
R ⊆ R′ (the epigraphs of P and P − ρ − η, truncated) through such queries to R and R′. Replacing R
by R_v keeps R_v ⊆ R′, changes Δ by one, and leaves every R′ query unchanged. So it too needs
≥ V(conv L) queries.

*Application.* With the box [0, H) × [⌊min P⌋, M] (truncated count HM − F(P) + const), each
ladder block needs ≥ V_b queries in this model. By Lemma V the ladder pays Ω(N^{1/3} log N) for every
even j. That holds whether it computes the floor sums, their differences through the two epigraphs,
or (for Prop D) F itself.

*Scope.* R_v is not a polynomial-arc region. An algorithm that evaluates P off the lattice, uses its
coefficients or N's arithmetic, or accesses the non-convex crescent between the arcs directly, is
outside the model. Li (arXiv 2609.29949) already remarks informally that beating 1/3 "requires a
mechanism beyond a straightforward local piecewise-linear resolution". Theorem T and Prop PL formalise
two versions of that remark.

**Prop PL (PROVED; interval linear-floor decompositions).**
- Suppose [0, H) is split into intervals I_1, …, I_r with ⌊P(h)⌋ = ⌊(p_k h + r_k)/q_k⌋ on I_k, where
  p_k, r_k ∈ Z and q_k ≥ 1.
- Then r ≥ V_b / O(log(H + max q_k)), with V_b as in Lemma V.

*Proof.*
1. A vertex of the lower hull of {(h, ⌊P(h)⌋)} that lies in I_k is a vertex of the lower hull of I_k's
   points.
2. After an integer shear and translation (which preserve vertex counts), 0 ≤ p_k, r_k < q_k. Since
   ⌊ℓ(u)⌋ ≥ ℓ(u) − 1 + 1/q_k, those points are the lowest lattice points of the rational polygon
   {u ∈ I_k, m ≥ ℓ(u) − 1 + 1/q_k}.
3. A planar rational polyhedron with O(1) facets has an integer hull with O(φ) vertices, where φ is
   the input size (Cook–Hartmann–Kannan–McDiarmid, Combinatorica 12 (1992) 27–37: at most
   2m^n(12n²φ)^{n−1}). Here φ = O(log(H + q_k)). ∎

*Coverage.* PL covers the Farey oracle (q ≤ H + 1, so the bound is V_b/O(log H)) and any exact
method that outputs such an interval decomposition. Hull walks are covered directly by Theorem T,
since they visit the V_b vertices. Not covered:
- codex's fixed-cell construction, whose per-cell congruence correction is not of this form;
- decompositions into arithmetic progressions, such as Polymath's Prop 2.2.

**E1 (pre-registered; every prediction held).** Brute-force block hulls (closed epigraph), 5 balanced
semiprimes per size.

| j | V_tot / V_glob (min–max) | slope of ln(V_tot/ln N) vs ln N | per-block V / (2^{1/3}cN^{1/3}/B) at 44 bits |
|---|---|---|---|
| 2 | 1.50–1.96 | 0.309 | 0.72 |
| 3 | 1.11–1.33 | 0.310 | 0.45 |
| 4 | 1.04–1.17 | 0.314 | 0.38 |

- At 44 bits V_tot = 160 770, 121 346 and 114 320 for j = 2, 3, 4, against V_glob ≈ 107 000–110 000.
  The block counts are 28 546, 7 570 and 3 070.
- Raising j from 2 to 4 removes 89% of the blocks but only 29% of the vertices. The excess over V_glob
  is per-block boundary overhead; the global hull is the floor.
- P1.1 (the Lemma V check) held on all 50 even-j instances.

## 3. The exact oracle: Farey dissection

**Lemma R (PROVED).** If q ≥ 1 and p, r ∈ Z with r ≤ qP(u) − pu < r + 1, then ⌊P(u)⌋ = ⌊(pu + r)/q⌋.

*Proof.* Write qP(u) = n + θ with n = pu + r, 0 ≤ θ < 1, and n = q⌊n/q⌋ + s with 0 ≤ s ≤ q − 1. Then
s + θ < q, so ⌊(n + θ)/q⌋ = ⌊n/q⌋. ∎

**Algorithm O** (`FloorSumConvex`).
- Input: P = num/den with num ∈ Z[h] and P″ of constant sign on [h_lo, h_hi).
  - The code rejects ranges whose end second differences have opposite signs. That check is
    sufficient for degree ≤ 3; for degree ≥ 4 the caller guarantees the precondition, which holds
    for ladder blocks. Before the check was added, a critic found silent errors on 2 of 3 000 cubics
    with an interior inflection. No ladder block was affected.
- Set β̂ = (max |second difference| at three points)/2, and Q = iroot(⌊1/β̂⌋ + 1, 3) + 1. This is
  ⌈β̂^{−1/3}⌉ or one more, capped at H + 1.
- From the current position h:
  1. p/q := the last continued-fraction convergent with q ≤ Q of the local slope s = P(h+1) − P(h);
  2. window [h, h + L) with L = max(1, ⌊1/(qQβ̂)⌋);
  3. on the window, g(u) = qP(u) − pu is convex (or concave). Split it at its integer extremum,
     found by binary search on the sign of g(u+1) − g(u). On each monotone side G = ⌊g⌋ is monotone,
     and its constant runs are found by binary search;
  4. each run with G ≡ r contributes floor_sum_linear for ⌊(pu + r)/q⌋ (Euclid, O(log)).
- Runs are stored with prefix sums, so a subinterval query costs O(log) more. Each (block, shift)
  pair is built once; the shift is absorbed into the constant term.

**Theorem O (PROVED).**
- The output is exact for every P with P″ of constant sign.
- If max|P″|/2 ≤ κβ̂ on the range, then:
  - #windows ≤ 8Hβ̂^{1/3} + 1;
  - #runs ≤ 4·#windows + 2(1 + 2κ)(1 + Hβ̂^{1/3});
  - the cost is O((1 + κ)(1 + Hβ̂^{1/3}) log H) exact evaluations on O(j log N + log H)-bit numbers.

*Proof.*
- **Exactness.** Lemma R applies on every run. The runs tile the range because G is monotone on each
  side of the split.
- **Windows.** L ≥ max(1, ⌊1/(Q²β̂)⌋).
  - If β̂ ≤ 1/8, then Q ≤ β̂^{−1/3} + 2 ≤ 2β̂^{−1/3}, so 1/(Q²β̂) ≥ β̂^{−1/3}/4 and
    #windows ≤ 8Hβ̂^{1/3} + 1.
  - If β̂ > 1/8, then #windows ≤ H ≤ 2Hβ̂^{1/3}.
  - If Q is capped, there is one window.
- **Runs.**
  - On a window, |g′(u)| ≤ |qs − p| + q·max|P″|·L < 1/Q + 2qκβ̂L. The last convergent with q ≤ Q
    has |qs − p| < 1/Q, and it is 0 if the expansion terminates.
  - So the variation over the window is ≤ L/Q + 2qκβ̂L² ≤ (1 + 2κ)L/Q, using L ≤ 1/(qQβ̂) when
    L > 1. A window with L = 1 has one run.
  - A monotone side with variation V takes at most V + 2 floor values (codex: an integer can be
    crossed with arbitrarily small variation). So a window has ≤ 4 + 2·variation runs.
  - Summing, with H/Q ≤ 1 + Hβ̂^{1/3}, gives the bound.
- **Cost.** Each window and each run costs O(log L) evaluations. ∎

Codex's counterexample to our first constant, P = ((h − 5)² − 4)/10⁶ on [0, 10), produces four runs.
That is within the corrected bound.

**For ladder blocks** κ = 1 + O(j/B), so T_j(H) = O(Hβ^{1/3} log H) with β = N/x³.
- At x ≍ √N this is H^{1/3} (j = 2), H^{4/9} (j = 3), H^{1/2} (j = 4), and in general H^{2(j−1)/(3j)}.
  Then B·T_j ≍ N^{1/3} for every j (critic-checked).
- Codex floorsum_1 proved the same order with a fixed-cell variant.
- The bound is an upper bound. On major arcs (slope within ≈ 1/(qQ) of p/q with small q) the windows
  are much longer and the cost drops. Codex's P = h + γh²/H² costs O(1).

**E2 (pre-registered).**
- **P2.1 ✓.** 0 mismatches in 3 000 random convex polynomials of degree 2–4, and in 5 400 ladder-block
  range queries (j = 2, 3, 4; 40, 48, 56 bits; three shifts).
- **P2.2 and P2.3 FAILED as specified (design error).**
  - The scaling blocks sat at x = ⌊√N/2⌋ + O(10³). There N/x² is within O(10³/√N) of the integer 4:
    the whole block is one major arc.
  - The cost ratio (windows + runs)/(Hβ^{1/3}) was 0.003–0.39 for j = 2 (slope −0.002 in ln H) and
    0.08–0.48 for j = 3 (slope 0.316).
  - j = 4 behaved as predicted (1.00–1.14, slope 0.493), because there the block's slope range ≍ 4/B
    exceeds the major-arc width.
  - The failure invalidates the lower-side scaling prediction, not the upper bound of Theorem O. The
    decision rule is NOT MET.
- **E2b (post-hoc, exploratory; not a confirmatory replacement).** The same sizes, with x uniform in
  [√N/4, √N], 8 blocks per size.
  - Cost ratio 0.51–2.88 (median 1.93) for j = 2, 1.82–2.14 for j = 3, 1.49–2.05 for j = 4.
  - Slopes against ln H: **0.337, 0.442, 0.506**, against 1/3, 4/9, 1/2. These are finite-size fits.

## 4. Fourier and smoothing

**What Hiary computes (KNOWN; checked by codex and the literature critic).**
- **Quadratic.** Hiary, Ann. of Math. 174 (2011) 859–889, Thm 1.1. It computes
  L^{−r}Σ_{n≤L} n^r e(an + bn²) for any a, b, to error Aν^κε, in ν^{O(1)} operations on O(ν²)-bit
  numbers, with ν = (r+1)log(L/ε) and no precomputation.
- **Cubic.** Hiary, Ann. of Math. 174 (2011) 891–946, Thm 2.2. Cubic coefficient c ∈ [0, L^{μ−3}],
  μ ≤ 1, after L^{4μ}ν^{O(1)} precomputation.
- **Coverage for the ladder.**
  - j = 2: every frequency.
  - j = 3: only frequencies k ≲ H^{1/3} have μ ≈ 0. At the optimum K ≍ √H below, μ ≈ 1/6 and the
    precomputation is H^{2/3} unless it is shared across blocks.
  - j ≥ 4: there is no quartic theorem.
- So for j = 2 the Weyl sums are not the bottleneck. For j ≥ 3 they are an additional obstacle.

**Identity.** With ψ(y) = y − ⌊y⌋ − ½ and the midpoint convention ψ₀(n) = 0,
⌊y⌋ = y − ½ − ψ₀(y) + ½·1_Z(y). The endpoint term must be kept (codex).

**Prop F (PROVED).** Let T be a real trigonometric polynomial of degree K and
G = {y ∈ R/Z : |ψ(y) − T(y)| < 1/4}.
1. G^c is a union of at most 2K + 2 intervals.
2. meas(G^c) ≥ 1/(51K).

*Proof.*
1. **Intervals.** On (0, 1), φ = ψ − T = y − ½ − T(y) has derivative 1 − T′, which has ≤ 2K zeros.
   So φ has ≤ 2K + 1 monotone pieces. Each of {φ ≥ ¼} and {φ ≤ −¼} then has ≤ K + 1 components;
   tangencies sit at critical points and are counted the same way. With the jump at 0 identified on
   R/Z, that gives ≤ 2K + 2 intervals.
2. **Measure.** Let μ = meas(G^c). On G, |T| < 3/4.
   - Suppose μ ≤ 1/(16πK). The Remez inequality for trigonometric polynomials (Erdélyi, J. London Math.
     Soc. (2) 45 (1992) 255–264) says ‖T‖_∞ ≤ e^{4K|B|}·sup outside B, for an exceptional set B of
     angular measure |B| < π/2. Here |B| = 2πμ ≤ 1/(8K). So ‖T‖_∞ ≤ (3/4)e^{1/2} < 1.24.
   - Bernstein's inequality (Zygmund, Trigonometric Series, Ch. X; Borwein–Erdélyi, GTM 161) gives
     |T′| ≤ 2π·1.24·K < 7.8K.
   - Let δ = 1/(4 + 31.2K). If T(0) ≥ 0, then on (0, δ), T ≥ −7.8Kδ while ψ < −½ + δ, so
     |ψ − T| ≥ ½ − δ(1 + 7.8K) = ¼. If T(0) < 0, the same holds on (−δ, 0).
   - So μ ≥ δ ≥ 1/(36K).
   - In both cases μ ≥ min(1/(16πK), 1/(36K)) ≥ 1/(51K). ∎

Ganzburg (Constr. Approx. 17 (2001)) sharpens the Remez constant, which would improve 1/(51K).

**Corollary F (PROVED; certified individual-correction model).**
- *Model.* The algorithm writes Σ_h ψ(P(h)) = Σ_h T(P(h)) + Σ_{h∈E}(ψ − T)(P(h)) + Err. It is charged
  K for evaluating all frequencies up to K (an explicit model assumption: degree K need not mean K
  nonzero terms) and |E| for the explicitly corrected set. It certifies
  |Err| ≤ Σ_{h∉E}|ψ − T|(P(h)) < ½.
- Then at most one h ∉ E lands in G^c.
- If the points {P(h)} of the block have maximal circular gap g, an interval of length ℓ contains
  ≥ ℓ/g − 1 of them, so |E| ≥ 1/(51Kg) − (2K + 2) − 1.
- Hence the cost C ≥ max(K, 1/(51Kg) − K − 3) ≥ (102g)^{−1/2} − 3 ≥ 0.099/√g − 3.

*Scope.* The model excludes:
- uncertified schemes that rely on cancellation;
- aggregated corrections of resonant clusters (codex);
- non-trigonometric smoothing.

Finding E is itself a near-integer location problem; the bound grants it for free.

Codex floorsum_1 gave two companion bounds: Ω(√H) for the worst shift (E_s|{h : ‖P(h) + s‖ < 1/K}| =
2H/K), and Ω(H^{2/3}) for the explicit family P = h + γh²/H².

**E4 (pre-registered; held).** 200 random j = 2 ladder blocks per size, x uniform in [√N/4, √N].

| bits | mean H | n_δ / (2δH) at δ = 1/⌈√H⌉ | min over K = 2^i of (K + n_{1/K}) / √H |
|---|---|---|---|
| 48 | 1 166 | 0.998 | 2.88 |
| 56 | 4 762 | 1.003 | 2.89 |
| 64 | 19 502 | 0.999 | 2.88 |

**E4b (post-hoc diagnostic, j = 2).**
- Median g·H = 5.2, 6.5, 7.9 (≈ ln H), max 151 at 48 bits.
- The explicit constants of Corollary F are vacuous at these sizes: the median bound is −0.043√H,
  −0.005√H and +0.014√H.
- The realistic figure is the measured model cost ≈ 2.9√H.
- Gaps are not uniformly small. Codex shows N = x² + 1 (j = 2) has g → 3/4.

**Consequence (experimentally supported for j = 2, HEURISTIC for j ≥ 3).**
- Where g ≲ (log H)/H, Fourier + certified correction costs Ω̃(√H) per block.
- The total is then ≳ B√H ≍ N^{1/4 + 1/(2(j+2))}: N^{3/8} at j = 2 (worse than the hull), tending to
  N^{1/4}, never below it.
- For j ≥ 5 this is below N^{1/3}. There the route is blocked by the missing degree-j Weyl-sum
  algorithms and by the cost of locating E, not by Prop F.
- Smoothing cannot make K polylog.

**Self-reference (PROVED, elementary).** For f = N/d the smoothing correction set at scale δ is
{d : ‖N/d‖ < δ}: divisors of N − k with |k| < δd, including k = 0.
- Li's de-smoothing for D(N) corrects by factoring the integers near N. For restricted sums that
  includes N itself, which is circular.
- Prior art for "the factorisations of nearby integers factor N": Sica, "Factoring with hints"
  (arXiv 0912.1585), which factors N deterministically from those of ≈ N^{1/3+ε} nearby integers.

## 5. Euclid for quadratics, reciprocity, and known exact counters

**Axis swap (PROVED).**
- For Q(h) = (ah² + bh + c)/m increasing and nonnegative on [0, H),
  Σ_{h<H}⌊Q(h)⌋ = Σ_{k=1}^{⌊Q(H−1)⌋}(H − t(k)), with t(k) = ⌈(−b + √(b² + 4a(mk − c)))/(2a)⌉ clipped
  to [0, H] (codex).
- Completing the square, the count is #{(y, t) : t² ≤ 4amy + Δ, t ≡ b (mod 2a)}. Counting by t returns
  the original sum, so the double swap is the identity.
- No modulus shrinks, because the inverse of a parabola is a square-root curve.

**Prop A (affine length; invariance PROVED and classical, interpretation HEURISTIC).**
- For a planar affine map A, λ(AΓ) = |det A|^{1/3}λ(Γ), where λ(Γ) = ∫κ^{1/3}ds = ∫|P″(h)|^{1/3}dh for a
  graph (Blaschke). So λ is invariant under SL₂^±(Z) ⋉ Z². Every linear-Euclid move (slope reduction,
  swap, reflection) is such a map.
- A line has λ = 0; a ladder block arc has λ ≍ cN^{1/3}/B. No invertible affine map turns a
  positive-length arc into a line.
- Generic convex lattice regions have hull-vertex counts governed by affine perimeter: Andrews (Trans.
  AMS 106 (1963)), Balog–Bárány (disc, DIMACS 1991), Bárány–Larman (Math. Ann. 312 (1998)), Bárány
  (J. reine angew. Math. 484 (1997); survey Bull. AMS 45 (2008)).
- This is a geometric analogy for the cost of local linearisation. It is **not** a proof that no
  quadratic Euclidean recurrence exists: an algorithm may split arcs, subtract contributions, or use
  compressed arithmetic representations.
- Likewise, "theta transformations exist for e(P) but not for ⌊P⌋" is a heuristic contrast: Poisson
  summation needs a smooth summand and is not a lattice automorphism.
- No prior statement of the affine-length argument for Euclid was found by us. Sladkey's algorithm
  already uses tangent lines and unimodular changes of coordinates.

**Unrestricted hardness (KNOWN + PROVED reduction; floorhard has the full picture).**
- Manders–Adleman (JCSS 16 (1978) 168–184; Garey–Johnson AN1): deciding whether some x < c has
  x² ≡ a (mod b) is NP-complete.
- Since ⌊(h² − A)/m⌋ − ⌊(h² − A − 1)/m⌋ = 1[m | h² − A], two quadratic floor sums decide it. Codex
  floorsum_1 and floorhard found this independently.
- Floorhard proves much more (notes/floorhard_results.md):
  - #P-hardness (Thm SS);
  - hardness within the Taylor class itself for uniform-in-N oracles (Thm TE), with the hard blocks at
    x ≤ N^{3/8+o(1)};
  - Cor DOM: a uniform oracle with cost H^δ factors m in O(log²m·T(m)).
- Our results are consistent with all of these. None of them touches top blocks.

**Complete periods (KNOWN).** Σ_{x<m}⌊ax²/m⌋ carries h(−m) (Dirichlet; floorhard Prop FP). Codex's
form Σ_{h<p}⌊h²/p⌋ = (p − 1)(p − 2)/3 − (1/p)Σ_{r<p} r(r/p) was checked by the literature critic.

**Known exact counters (KNOWN; verified by the literature critic).**
- **Har-Peled**, Comput. Geom. 10 (1998) 125–138: output-sensitive discrete convex hulls in
  O(|C₀| log δ), so a radius-R disc in O(R^{2/3} log R) (for sufficiently round bodies). The exact count
  then follows from Pick's theorem.
- **Sladkey**, arXiv 1206.3369: exact D(N) in O(N^{1/3}) time and O(log N) space, without factoring.
- **Polymath (Tao–Croot–Helfgott)**, "Deterministic methods to find primes", Math. Comp. 81 (2012)
  1233–1246, §2.1: N^{1/3+o(1)}.
- **Li**, arXiv 2609.29949 (Sept 2026; not audited): N^{9/28+o(1)} expected time. It is randomized only
  through Lenstra–Pomerance factoring of the integers near N.
- **Alcántara–Blanco–Criado–Santos**, arXiv 2501.19193: O(log)-per-vertex hull enumeration, and
  (their Cor 1.4) a deterministic hull-walk factoring algorithm in O(N^{1/3} log² N). This is prior art
  for "a hull walk factors in N^{1/3}", alongside Lehman (Math. Comp. 28 (1974)) and Hiary
  (Math. Comp. 85 (2016), arXiv 1408.2608).
- We found no factorisation-free exact D(N) algorithm or lattice-points-under-a-parabola algorithm
  below the hull complexity. We found no reciprocity law for incomplete quadratic floor sums; the
  known ones cover linear arguments with polynomial weights (universal Euclid; e.g. arXiv 2107.08308).

**Amortisation (HEURISTIC; card F6).** Consecutive Farey windows form Stern–Brocot chains. Chains are
long only near major arcs, so aggregating them gains a constant, not an exponent.

## 6. Arithmetic structure: what the ladder really needs

### 6.1 High bits (reformulation, PROVED; also codex)

At x ≍ √N an L_j block has length H ≍ N^{1/2−1/(j+2)}. Locating p within H means knowing its top
(1/(j+2))log₂N bits, i.e. a fraction 2/(j+2) of p's bits. This holds up to additive constants and
binary alignment; it describes interval uncertainty. The ladder guesses the interval (B per dyadic
range) and then solves one block.

| j | block at √N | known top bits / log₂N | ladder cost with a polylog oracle | known block solver |
|---|---|---|---|---|
| 1 | N^{1/6} | 1/3 | Õ(N^{1/3}) | Euclid; Hiary (2016); Rivest–Shamir (EUROCRYPT '85, LNCS 219, 31–34) |
| 2 | N^{1/4} | 1/4 | Õ(N^{1/4}) | Coppersmith (J. Cryptology 10 (1997) 233–260); Howgrave-Graham (CaLC 2001, LNCS 2146, 51–66) |
| 3 | N^{3/10} | 1/5 | Õ(N^{1/5}) | none known |
| 4 | N^{1/3} | 1/6 | Õ(N^{1/6}) | none known |

- The j = 2 row ("guess the top quarter of p's bits, then run Coppersmith") is cryptographic folklore,
  i.e. KNOWN.
- The classical N^{1/4} is Pollard (1974) / Strassen (1976). The deterministic record is Harvey
  (Math. Comp. 90 (2021)) N^{1/5}, with Harvey–Hittmeir's log-log saving; Hittmeir gave N^{2/9}.
- *Caveat (codex).* A floor oracle locates a block's candidates in O((1 + C_b) polylog) calls, and a
  single block can hold Ω(√H) candidates. So the table is an implication that needs a per-block
  candidate bound, not an equivalence.
- *State of the art (literature critic, codex).* No polynomial-time classical method is known that
  uses only (1/4 − ε)log₂N contiguous high bits of one balanced factor, for fixed ε > 0, without
  further information. Below the Coppersmith bound (EUROCRYPT '96) one brute-forces the missing bits;
  see De Micheli–Heninger's survey, §4.2.2.
  - Special settings do better: p^r q (Boneh–Durfee–Howgrave-Graham 1999); random known bits of p
    and q (Heninger–Shacham 2009); implicit factoring (May–Ritzenhofen 2009); Maurer's oracle
    factoring.
  - None covers the ladder's intervals. There is no impossibility theorem; the barrier is that of the
    lattice method.

### 6.2 Lattice containment for divisors (PROVED from KNOWN)

- **The theorem.** Coppersmith / Howgrave-Graham / May find the divisors p ≥ N^α of N with
  |p − x| ≤ N^{α²} in polynomial time. See Cohn–Heninger, "Ideal forms of Coppersmith's theorem and
  Guruswami–Sudan list decoding" (ICS 2011; arXiv 1008.1284), Thm 1.1 with d = 1.
- **j = 2.** The block at x = N^α has H ≤ ½N^{α−1/4} ≤ ½N^{α²}, since (α − ½)² ≥ 0. So lattice
  reduction is a polylog block oracle for the j = 2 ladder's real need, which is locating divisors,
  not computing floor sums.
- **j = 3.** The containment holds only for α ≤ (1 − √(1/5))/2 ≈ 0.2764.
- **Cost of the lattice method alone.** Its own cost per dyadic range is N^{α−α²}, maximal (N^{1/4}) at
  α = ½.

**E5 (Amendment A1; KNOWN algorithm).**
- `lattice_factor` with dimension ≤ 16, (m, n) = (7, 16) at √N: 12/12 correct, 36–48 bits.
- At √N its interval is 0.68–0.78 of the ladder's j = 2 block. That band was set after a peek and
  carries no evidential weight.
- Fixed (m, n) = (7, 16) gives radius exponent 7/30 rather than 1/4 at √N (codex). The asymptotic
  theorem needs dimension ≍ log N.
- Python cost: 48 s at 48 bits, about 52 000 LLL calls of dimension 16. This is a demonstration, not a
  new claim.

### 6.3 D(N) via the ladder (Prop D, PROVED; `divisor_summatory_ladder`, tested j = 1–4)

**Prop D.** With an exact L_j oracle, D(N) = 2Σ_{d≤√N}⌊N/d⌋ − ⌊√N⌋² takes O(B log² N) oracle calls
plus O(B log N) further operations.

*Proof.*
1. Per block, Σ_h⌊N/(x+h)⌋ = F(P_j; 0, H) − Σ_{h∈C}(⌊P_j(h)⌋ − ⌊N/(x+h)⌋).
2. For even j, ⌊P⌋ − ⌊f⌋ ∈ {0, 1} and is nonzero only on candidates. For odd j it lies in {−1, 0},
   again supported on candidates. At h = 0 the correction is 0.
3. C is located by bisection; the d < 2B are summed directly; the total candidate count is
   O(B log N). ∎

The critic checked it exhaustively for N < 4000, j = 1–6, and on special N up to 10^{11}.

So a polylog L_2 oracle would give D(N) in Õ(N^{1/4}). We found no factorisation-free algorithm below
N^{1/3}, and the oracle would also beat Li's factorisation-assisted N^{9/28}. A polylog L_3 oracle would
give Õ(N^{1/5}). This is evidence that the oracle is a strong primitive, not a reduction from a hard
problem. A sum over fewer than N^{1/3} denominators is cheap by direct evaluation; the issue is the
total cost over all blocks.

### 6.4 Top blocks: window queries by lattice reduction (E6; orchestrator request)

Floorhard (§4.2) shows:
- exact-value queries at top blocks are Coppersmith-easy (T3′);
- for j = 2, divisor detection is known (T2);
- the open part is windows of positive width ρ.

We tested the natural lattice for windows.

**Construction (`lattice_window.py`).**
- Write P_j(h) = num(h)/D with D = x^{j+1}. Then {P_j(h)} ∈ [0, R/D] iff some r ∈ [0, R] satisfies
  f(h, r) := r − g(h) ≡ 0 (mod D), where g = num mod D and R = ⌈ρD⌉.
- Use the basic-strategy lattice h^u f^k D^{m−k}, 0 ≤ k ≤ m, 0 ≤ u ≤ j(m−k). It is triangular, with
  monomials h^a r^b for b ≤ m and a ≤ j(m−b).
- Rows with n‖v‖² < D^{2m} vanish at every root in the box (Howgrave-Graham). Two of them give the
  roots via a resultant, or via a common factor linear in r.
- **Generic determinant bound (PROVED for this lattice, as m → ∞).** The exponent sums are j²m³/6
  (X), jm³/6 (Y) and jm³/3 (D), against D^{m·jm²/2}. This gives X^j·Y < D.

**Pre-registered prediction and its failure.**
- A2 predicted X_max ≈ ρ^{−1/j}, from the generic bound with Y = ρD.
- E6 placed blocks at x = ⌊√N⌋ − O(10⁴). There N/x² = 1 + O(10⁴/√N), again a major arc: the block's
  linear part does not wrap.
- On those blocks the lattice resolves the entire j = 2 block (X_max ≥ H at 64 and 80 bits). In a
  scratch diagnostic on one such 48-bit block, the good rows shared the factor (r − G(h))², where G is
  a non-wrapping integer representative of g. For j = 3 it resolves about H/14.
- P6.1 held but was vacuous (no roots in the boxes). P6.2 and P6.3 FAILED. The decision rule is NOT
  MET.

**E6b (post-hoc, disclosed; generic top blocks, x uniform in [√N/2, √N]; 4 blocks × 3 sizes × j = 2, 3;
m ≤ 4).**
- The lattice does better than the generic bound: log₂(X^jR/D) reaches +3.7 (j = 2) and +12 (j = 3).
- But **X_max/(xN^{−1/3}) ∈ [0.69, 1.75] on all 24 blocks, for both j.**
- In units of the block, ln(H/X_max)/ln H ≈ 0.27–0.33 (j = 2; → 1/3) and 0.38–0.45 (j = 3; → 4/9).
- Raising m from 2 to 4 changes X_max by at most a factor 2.

**E6c (post-hoc).** A root is planted at random h₀ < X_max on 12 generic blocks (j = 2: m = 3;
j = 3: m = 2; 48 and 64 bits). Every good row vanished at it, and h₀ was recovered in 12/12.

**Explanation (HEURISTIC).** At a top block x > N^{1/3}, so for i ≥ 2 the coefficients of g are
c_i = ±N x^{j−i}: small, not reduced mod D. Absorb the nonlinear part s = r − Σ_{i≥2}c_ih^i into a
linear congruence c₀ + c₁h ≡ s (mod D). A 2-D lattice solves that when X·max(R, N x^{j−2}X²) ≲ x^{j+1},
i.e. **X ≲ xN^{−1/3}**. That is Hiary's j = 1 block length (Hiary 2016; Rivest–Shamir). The bivariate
lattice, at the dimensions tried, does no better.

**Answer to the orchestrator.**
- In this construction a lattice window oracle needs ≈ H/(xN^{−1/3}) = cN^{1/3}/B calls per top block.
  That matches the Farey cost: H^{1/3} for j = 2, H^{4/9} for j = 3. The total is N^{1/3}.
- Coppersmith-type methods do not solve positive-width windows below the trap at top blocks. They
  solve exact values (floorhard T3′) and, for j = 2, divisor detection.
- Larger m, other shift sets or other lattices are not excluded. This is a measured statement about
  the standard construction, plus a heuristic explanation.
- Major-arc blocks (slope within ≈ 1/(qH) of p/q with small q) are cheap. By a HEURISTIC count they
  cover only O(Q) of the B top blocks for denominators q ≤ Q.

### 6.5 With factoring allowed the difference interface is easy (codex floorsum_2)

Every candidate gives d = x + h and k = N − dm with |k| ≤ C_jB. So:
1. factor every N − k in that range once (Lenstra–Pomerance, expected L_N[½, 1 + o(1)] each);
2. enumerate the divisors d ≤ √N, assign them to blocks, and check the exact rational candidate
   predicate;
3. answer subinterval counts by rank queries.

Codex's prototype reconstructed the exact candidate sets on 16 instances. The cost is
N^{1/(j+2)+o(1)} expected, below N^{1/3}. But k = 0 means factoring N itself, so this is circular as a
factoring method. It shows only that no unconditional "the difference interface costs the hull"
statement is possible. The meaningful question is the deterministic, factorisation-free one.

## 7. End to end: the operational ladder (E3, pre-registered; all predictions held)

`ladder_factor_op(N, j)` is analytic's ladder with the simulated oracle replaced by Algorithm O. Primes
of 28–44 bits (3 per size) give the full-range cost. Balanced semiprimes (5 per size) check
correctness.

| bits | j = 2 evals | j = 2 blocks | j = 2 windows+runs | j = 3 evals | j = 3 blocks | wall (j = 2) |
|---|---|---|---|---|---|---|
| 28 | 58 953 | 1 164 | 11 085 | 50 115 | 532 | 0.04 s |
| 32 | 150 497 | 2 459 | 27 889 | 129 402 | 996 | 0.11 s |
| 36 | 427 421 | 5 595 | 77 308 | 371 677 | 1 968 | 0.30 s |
| 40 | 1 236 638 | 12 862 | 218 550 | 1 099 327 | 3 922 | 0.85 s |
| 44 | 3 265 620 | 27 211 | 561 917 | 2 952 360 | 7 277 | 2.2 s |

- P3.1 ✓: 80/80 correct.
- P3.2 ✓: the slope of ln(evals/ln²N) against ln N is 0.290 (j = 2) and 0.296 (j = 3), inside
  [0.28, 0.40]. These are finite-size fits.
- P3.3 ✓: evals(j = 3)/evals(j = 2) = 0.85–0.90. Going from j = 2 to j = 3 cuts the blocks by 2.2×
  (28 bits) to 3.7× (44 bits), yet removes only 10–15% of the work.

**End-to-end exponents.**
- With the best oracle we have: N^{1/3+o(1)} for every j (Theorem O; E3).
- With oracles nobody has: N^{1/(j+2)} (Theorem L_j).
- For j = 2, the lattice route already gives the KNOWN N^{1/4+o(1)}.

## 8. Pre-registration outcomes

| Experiment | Outcome | Decision |
|---|---|---|
| E1 | P1.1 ✓ (50/50), P1.2 ✓ (1.04–1.96), P1.3 ✓ (0.309, 0.310, 0.314), P1.4 ✓ (0.72, 0.45, 0.38) | C1 supported |
| E2 | P2.1 ✓. **P2.2, P2.3 ✗** (major-arc design error) | rule NOT MET; post-hoc E2b exploratory only |
| E3 | P3.1 ✓, P3.2 ✓, P3.3 ✓ | C3 supported |
| E4 | P4.1 ✓ (0.998–1.003), P4.2 ✓ (2.88–2.89) | C4's premise supported (j = 2) |
| E5 (amended) | P5.1′ ✓ 12/12; P5.2′ ✓ (0.68–0.78; no evidential weight) | demonstration |
| E6 (A2) | P6.1 ✓ (vacuous); **P6.2, P6.3 ✗** (major arc again; the prediction also ignored the small nonlinear coefficients) | rule NOT MET; post-hoc E6b/E6c |

The two failures share one cause. Blocks placed at x ≈ √N/k have slopes near integers, so they are
unrepresentative. Both are disclosed above.

## 9. Barrier summary (this agent)

| # | Statement | Scope | Status |
|---|---|---|---|
| T | ≥ V_b membership/LP queries per block, also for the two-epigraph difference; Σ_b V_b = Ω(N^{1/3} log N) | deterministic convex-access algorithms; even j | PROVED (Lemma V + Balog–Bárány) |
| PL | Interval linear-floor decompositions need ≥ V_b/O(log(H + q)) pieces | Farey and similar; not AP decompositions | PROVED |
| A | Lattice moves preserve affine length | interpretation as a Euclid obstruction | invariance PROVED (classical); interpretation HEURISTIC |
| F | Certified smoothing costs ≥ 0.099/√g − 3 per block | bandwidth-charged trigonometric approximants + individual correction | PROVED; ≍ √H per block experimentally supported (j = 2) |
| LW | Lattice window oracles (basic strategy, m ≤ 4) reach sub-blocks of length ≈ xN^{−1/3} at top blocks | this construction | experimentally supported (post-hoc E6b/E6c) + HEURISTIC explanation |
| MA | Unrestricted quadratic floor sums ⇒ P = NP | arbitrary coefficients | KNOWN + 2-call reduction (floorhard: #P) |
| D | Polylog L_2 ⇒ D(N) in Õ(N^{1/4}) | — | PROVED implication |
| HB | Polylog L_j + per-block candidate bound ⇒ factoring from (1/(j+2))log₂N top bits | — | PROVED implication |

## 10. Reviews and fixes

**GPT-6.1-sol blind round (floorsum_1).** It reached the same verdict independently: an Õ(H^{1/3})
CF-linearisation oracle and N^{1/3} total for every j. It contributed:
- the coefficient-size correction;
- the per-block candidate counterexample;
- the Manders–Adleman reduction;
- the shift-averaging Fourier bound;
- the high-bits caveat.

**GPT-6.1-sol adversarial round (floorsum_2): ACCEPT WITH FIXES, major revision.** Fixes applied:
- Theorem T: one-optimiser LP interface with a deletion-stable tie rule, deterministic algorithms, a
  truncating box; it was re-scoped to individual sums or two-epigraph differences ("REJECT as scoped").
- Theorem O: the run constant corrected (+2 per side; counterexample reproduced in the text); the
  capped-Q term; the explicit κ hypothesis.
- Prop F: bandwidth charged as a model assumption; Ω̃(√H) only under a gap hypothesis; the N = x²+1
  gap counterexample; the ladder consequence made conditional; Hiary's coverage restricted.
- Prop A: relabelled as a heuristic interpretation.
- High bits: "equivalence" withdrawn; the precise open-problem wording; radius 7/30 at fixed
  dimension.
- Prop D: literature phrased as a search report.
- Pre-registration: wording ("before the formal runs"; E2b exploratory; E5 band).
- The factorisation-assisted bypass, as §6.5.

**Opus correctness/cost critic: ACCEPT WITH FIXES, no FATAL.** It verified every proof and every table
entry, and found 0 mismatches in 19 755 adversarial exactness cases. Fixes applied:
- the difference version of Theorem T (MAJOR);
- labels for the Fourier consequence and the odd-j totals (MAJOR);
- the E3 wording; post-hoc tags in §0; the negative E4b medians; Hiary's scope;
- the convexity precondition check in the code (found silent errors on inflected cubics) and its test;
- Q documentation; the PL normalisation and coverage; the Theorem T count (≥ V);
- the Prop F part 1 proof; the E5 band; the logged integral-point check (`e1int`); the B rounding note.

**Opus literature critic: OK WITH FIXES.** Fixes applied:
- the Polymath pointer (§2.1, N^{1/3+o(1)});
- Li is randomized;
- Alcántara et al.'s Cor 1.4 credited as prior art for hull-walk factoring;
- Hiary 2016 and the factoring baselines (Lehman, Pollard–Strassen, Hittmeir, Harvey, Harvey–Hittmeir)
  cited;
- the Prop A attribution moved to Andrews / Bárány–Larman / Bárány / Jarník;
- Li's informal piecewise-linear remark credited;
- Sica 2009 added;
- the Cohn–Heninger title; the Erdélyi range; a Bernstein reference;
- the special cases of the high-bits problem; the j = 2 row marked as folklore.

**Orchestrator request (after floorhard).** It asked for a top-block window oracle by lattices.
Answered in §6.4 (E6, E6b, E6c): lattices reach Hiary's scale only.

## 11. Reproduction

```bash
.venv/bin/python -m pytest -q tests/test_floorsum.py
for e in e1 e1int e2 e2b e3 e4 e4b e5 e6 e6b e6c summary; do
  PYTHONPATH=. .venv/bin/python -m src.floorsum.run_experiments $e
done   # e1 ~4 min, e5 ~6 min, e6/e6b ~3 min each, the rest ~1 min
```
