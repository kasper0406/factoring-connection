# modrec_1: blind attempt at a poly(log) lattice-count residue (builder brief)

You are a research mathematician/programmer. Goal: BUILD an exact evaluator, running in time
polynomial in log T (or at least T^{o(1)}), for one of the residues below, or show precisely why
the mechanism you tried cannot work. Work in your current directory; write code there; python3
with sympy/numpy is fine. Do not use the factorization of any input inside an evaluator (factors
may be used only in independent brute-force checks).

## Targets (any one suffices; P1 is the simplest)
- P1. Octant circle parity: P(T) = #{(x,y) in Z^2 : 0 < x < y, x^2 + y^2 <= T} mod 2, for
  arbitrary integer T. For N = pq = 1 mod 4, P(N) - P(N-1) = [N is a sum of two squares] mod 2
  (r_2(N)/8 is 2 or 0 here, so actually use mod 4 of the octant count: specify which residue
  your method yields and check this claim yourself).
- P2. Hyperbola parity with a line cutoff: H(T,X) = sum_{d<=X} floor(T/d) mod 2. Since
  H(N,X) - H(N-1,X) = #{d <= X : d | N} mod 2, a fast evaluator factors N by bisection on X.
- P3. Fixed-form box count: S(N) = #{0 <= x,y <= isqrt(N) : N | x^2 + x*y + 6*y^2} mod 8
  (one prescribed cutoff; this carries a hidden class/genus bit for discriminant -23).

## Known context (do not rediscover)
- Exact counts of lattice points under conic arcs cost ~T^{1/3}: every exact method walks the
  integer convex hull, which has Theta(T^{1/3}) vertices (Balog-Barany). Tao-Croot-Helfgott
  (Polymath4, Math. Comp. 2012) compute D(x) = sum tau(n) exactly in x^{1/3} and parity of pi(x)
  in x^{1/2-c}. Linear floor sums have exact Euclid reciprocity; complete quadratic Gauss sums
  have exact Landsberg-Schaar reciprocity; incomplete quadratic sums have only approximate
  (Hardy-Littlewood/Hiary) reciprocity with real error terms.
- Mod-2 tricks that DO work: theta(q)^2 = theta(q^2) mod 2; swap symmetry fixes the parity of
  the full circle / full hyperbola count (sum_{d<=T} floor(T/d) = floor(sqrt T) mod 2).
- The hypothesis to test: a count taken only mod 2^k need not respect the hull barrier.

## What to try (pick what you think is most promising; be a builder)
Exact reciprocity/descent with integer bookkeeping kept mod 2^k: SL2(Z) changes of variable,
Euclid/continued-fraction descent of the form or region, piecewise-affine involutions
(Zagier/windmill style: count = #fixed points mod 2), theta inversion done exactly mod 2^k,
modular forms mod 2^k, Barvinok-type short generating functions reduced mod 2, or anything else.
A valid step must map the prescribed instance to a strictly smaller instance plus poly(log)
boundary terms.

## Deliverable (write to your final message, terse)
1. For each mechanism tried: what it computes exactly, verified against brute force on small
   T, and the measured growth of the number of terms/nodes/corrections with T (give numbers at
   several T up to at least 2^36 where feasible). Kill test: corrections growing like T^{c}.
2. Any poly(log) evaluator: code path, brute-force agreement counts, timings at 35-40 bits.
3. The sharpest obstruction or open step you found, with prior-art pointers you actually checked.
Label every claim proved / experimentally supported / conjectured / refuted.
