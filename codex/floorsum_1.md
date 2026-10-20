# floorsum_1: blind round. FAST EXACT FLOOR SUMS OF DEGREE-j POLYNOMIALS ON SHORT INTERVALS

You are a research mathematician (computational number theory, analytic number theory, lattice
algorithms) with web search. The repository at /Volumes/git/factoring-connection is read-only
context (an AI-assisted factoring research project). Write scratch files only under your current
directory. Read notes/analytic_results.md §3 ("Reduction L_j: the floor-sum ladder") and
src/analytic/ladder.py first; you may cite them instead of re-deriving.

Only CLASSICAL algorithms matter. Be honest: a precise obstruction is as valuable as an algorithm.
Label every claim PROVED (give the proof), KNOWN (give an exact reference you have checked),
HEURISTIC, CONJECTURED, or UNVERIFIED.

## Setting (Theorem L_j, proved in the notes)

Fix j. B = N^{1/(j+2)}, blocks [x, x+H) with H = ⌊cx/B⌋, B/c ≤ x ≤ √N, and
P(h) = Σ_{i=0}^{j} (−1)^i N h^i / x^{i+1} (Taylor polynomial of N/(x+h); denominator x^{j+1}).
If an oracle returns F(P + s; h₁, h₂) = Σ_{h₁≤h<h₂} ⌊P(h) + s⌋ EXACTLY (rational shifts s,
arbitrary subintervals), then all divisors of N up to √N are found with O(B log² N) oracle calls.
At x ≈ √N the block has H ≈ c N^{1/2 − 1/(j+2)}, slope ≈ 1, and the top coefficient has
β_j H^j ≈ c^j (the arc bends by O(1) over the block). Total cost ≈ B · T_j(H).

Targets: j = 2 needs T_2 < H^{1/3} (to beat N^{1/3}); j = 4 needs T_4 < H^{1/10} (to beat the
deterministic record N^{1/5}); polylog T_j for all fixed j gives N^{1/(j+2)+o(1)}.

## Tasks

1. **The convexity trap.** Integer-hull / convex-arc lattice-point counting costs about
   H·(P″)^{1/3} per block. Verify (or refute) that summed over the ladder's blocks this is
   Θ(N^{1/3}) up to logs, independently of j. Is there any hull-type method whose cost depends on
   the arc's degree rather than its curvature?
2. **Fourier / smoothing.** Write ⌊y⌋ = y − 1/2 − ψ(y). Each frequency k gives an incomplete degree-j
   Weyl sum Σ_h e(kP(h)). Check precisely what Hiary (Annals 2011, quadratic; cubic in the ζ paper)
   computes: precision, cost, coefficient ranges. Then analyse smoothing (Vaaler / Beurling–Selberg
   degree-K approximants) plus EXACT correction of the h with ‖P(h)‖ < 1/K. What is the best
   T_j(H) such a scheme can achieve? Give a lower bound in a clearly stated model if you can. Can K
   be polylog?
3. **Euclid for quadratics.** Is there an exact recursive reduction (like the j = 1 Euclid / floor_sum
   algorithm, or the "universal Euclidean algorithm") for Σ_{h<H} ⌊(ah² + bh + c)/m⌋? Analyse the
   axis swap (the inverse of a parabola is a square-root curve), unimodular transformations of
   parabolic regions, and any reciprocity law for quadratic floor sums or for lattice points under
   parabolas (Gauss/Jacobi/theta analogues). Survey known exact algorithms for lattice points under
   parabolas, in circles/ellipses (Gauss circle problem computations), and the divisor summatory
   function D(N) = Σ_{n≤N} ⌊N/n⌋ (is anything below N^{1/3} known WITHOUT factoring?). Give exact
   references.
4. **Arithmetic structure.** P is the Taylor polynomial of N/(x+h); the ladder really only needs the
   O(1) h per block with ‖P(h)‖ ≲ 1/H, i.e. divisors in [x, x+H) of N − k with small |k|. Relate the
   ladder's actual need to known problems (e.g. factoring with known high bits of p, Coppersmith /
   Howgrave-Graham's N^{β²} bound, Rivest–Shamir). For each j, what does a polylog block oracle
   imply about those problems? Can lattice reduction serve as the block oracle for j = 2? For j = 3?
5. **Your best algorithm.** Give the best exact oracle you can for j = 2 in the ladder regime
   (βH² ≈ 1, coefficient denominators ≈ H^{O(1)}) with a cost analysis. Anything below H^{1/3}?
   For j ≥ 3? If you find nothing, state the sharpest obstruction you can prove.

Deliverable: a self-contained report (≤ 4000 words) in your final message, ordered by the tasks,
verdict first.
