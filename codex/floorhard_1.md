# floorhard_1 — blind round: how hard is EXACT polynomial floor-sum evaluation?

You are a research mathematician / complexity theorist. Work rigorously; label every claim
PROVED (give the proof), KNOWN (give a precise reference you are confident exists), HEURISTIC,
CONJECTURED, or OPEN. Do not invent references; if unsure, say "UNVERIFIED". The repository is
read-only for you; write any scratch files (scripts, the full report) under the current working
directory (codex/scratch/floorhard_1/). Python 3 with sympy/gmpy2 may or may not be available;
small exact experiments are welcome. Your final answer (stdout / -o file) should be the full report
(it may be long, ≤ ~6000 words).

## Background (a proved reduction, from this project)

Theorem L_j. Fix j ≥ 1. Let B = ⌈N^{1/(j+2)}⌉. Suppose an oracle returns EXACTLY
  F(P; h₁, h₂) = Σ_{h₁ ≤ h < h₂} ⌊P(h)⌋
for rational polynomials P of degree ≤ j with O(j log N)-bit coefficients (arbitrary rational constant
shifts, arbitrary subintervals), in time T_j. Then N can be completely factored in time
≈ N^{1/(j+2)}·polylog(N)·T_j. The P used are the Taylor polynomials of f(h) = N/(x+h) at blocks
[x, x+H) with H ≍ x·N^{−1/(j+2)}, shifted by rational constants: P(h) = Σ_{i≤j} (−1)^i N h^i / x^{i+1} + c.
The block length is H ≤ N^{1/2 − 1/(j+2)}, and the coefficients have the shape: slope ≍ N/x² , i-th
coefficient ≍ N/x^{i+1}, with |coef_i|·H^i ≍ N H^i/x^{i+1}, so the top coefficient times H^j is
about N^{(j+1)/(j+2)}... (check this yourself) and P is "nearly a hyperbola". Divisors x+h | N are
lattice points (h, m) with N/(x+h) = m, i.e. points within ρ < 1 of the arc. Degree 1 floor sums are
computable in O(log) by the Euclid-like floor_sum recursion, giving an N^{1/3} factoring algorithm.
For j = 2, polylog oracle ⇒ N^{1/4}; j = 3 ⇒ N^{1/5} (the deterministic record, Harvey 2021);
j ≥ 4 would beat the record. Also known in the project: complete quadratic floor sums
Σ_{x<m} ⌊a x²/m⌋ = a(m−1)(2m−1)/6 − (m−1)/2 + (a/m)·h(−m)·(…) for prime m ≡ 3 mod 4
(Dirichlet class-number formula); and complete cubic Weyl sums Σ_{x mod N} e(ax³/N) to O(1)
absolute error factor N = pq with p ≡ q ≡ 1 mod 3.

## The question (hardness side)

What would a FAST exact floor-sum oracle imply? Where exactly does the hardness lie?

1. Reductions FROM fast exact degree-2 floor sums. If Σ_{h₁≤h<h₂} ⌊(a h² + b h + c)/m⌋ (all integers,
   O(log N) bits, arbitrary intervals) were computable in polylog(N) time, what follows?
   Consider: class numbers h(−m) for prime m and composite m (Dirichlet; incomplete character sums
   Σ_{n<m/2}(n/m), Σ (n/m) n); quadratic residuosity mod composite m; factoring m itself (note:
   ⌊an²/m⌋ sums vs residues an² mod m, and counts like #{n ≤ X : an² mod m < Y}); point counting on
   conics; partial Legendre/Jacobi-symbol sums Σ_{n≤X}(n/p) (Burgess-range questions). Decide
   precisely which interval lengths and coefficient sizes each reduction needs. In particular:
   complete-period sums may be easy via closed forms while short or partial sums are hard, or the
   reverse. Is "polylog exact degree-2 floor sums over arbitrary intervals" known-hard (e.g. implies
   computing h(−p) in polylog, or factoring), and, if it implies class numbers, is h(−p) believed
   to need super-polylog time (best known: L[1/2] under GRH, Hafner–McCurley; Q^{1/5} Shanks)?
2. The RESTRICTED setting. Theorem L_j needs only SHORT intervals (H ≤ N^{1/2−1/(j+2)} ≪ m) and the
   special Taylor polynomials of N/(x+h) (nearly hyperbolic arcs, denominators x^{j+1}). Do the
   reductions in (1) survive in this restricted class? If not, what is the sharpest statement of the
   form "an oracle for class C (restricted) in time T would imply X"? Is there any hardness at all for
   exactly the oracle class Theorem L_j needs, beyond the trivial "it factors N in N^{1/(j+2)}T"?
   Conversely, is there a reduction showing the restricted class is as hard as the general one
   (e.g. by Euclid/CF steps or unimodular changes of variables mapping one to the other)?
3. Degree ≥ 3. Relate exact cubic floor sums Σ⌊(a h³ + …)/m⌋ to cubic Weyl/Gauss sums, cubic
   residuosity, class numbers/Hecke characters for Q(√−3) or cubic fields. Does polylog exact
   degree-3 floor-sum computation over full periods mod a composite N factor N? (Note the floor
   sum is not obviously the exponential sum; relate them via Σ⌊f/m⌋ = (Σf − Σ(f mod m))/m and the
   distribution of residues x³ mod m in intervals.)
4. Lower-bound / complexity literature. Exact lattice-point counting under curves: Gauss circle
   exact Σ r₂ algorithms; the Dirichlet divisor problem D(x) = Σ_{n≤x} d(n) exactly (O(x^{1/3})
   known: Vinogradov hull / Tao–Croot–Helfgott / Sladkey; anything better? conditional lower bounds or
   equivalences?); exact prime counting π(x) (Lagarias–Miller–Odlyzko x^{2/3}, Lagarias–Odlyzko
   analytic x^{1/2+ε}, Helfgott–Thompson x^{3/5}, recent?). Are there known equivalences or
   conditional lower bounds like "computing D(x) exactly faster than x^{1/3−ε} would imply …"?
   Counting lattice points in polygons / under convex curves (Barvinok, Lagarias–Odlyzko–type),
   known hardness of counting points under algebraic curves of degree ≥ 2 (e.g., #P-hardness of
   counting solutions to quadratic congruences / Manders–Adleman NP-completeness of
   x² ≡ a mod b with x ≤ c — is that relevant? Note that floor sums count lattice points under
   parabolas: Σ⌊(ah²+bh+c)/m⌋ counts (h, y) with m y ≤ a h² + b h + c, which is related to
   x² ≡ ? mod m with bounded x).
5. Give the sharpest barrier statement you can PROVE, e.g.
   "for j ≥ 2, an exact floor-sum oracle on [class] in time T would imply [X]", with X open or
   believed hard; or the opposite: explain why no such reduction exists for the restricted class
   and give a plausibility argument that a fast restricted oracle might exist. Estimate,
   honestly, the probability that a polylog (or H^{o(1)}) restricted degree-2 oracle exists.

Be concrete; small exact computations to confirm identities are encouraged.
