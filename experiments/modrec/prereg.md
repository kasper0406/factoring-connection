# modrec pre-registration (2026-11-27, before any run)

Question: is there a poly(log N) exact evaluator for a lattice-count residue mod 2^k that the repo
has shown factors N (floor8 short-square query; centred ellipse mod-16 count), via a reduction step
that exploits working mod 2^k rather than exactly?

## E1 Main hypothesis: Euclid/SL2 hull descent with mod-2^k bookkeeping
- Build a factor-blind exact counter for #{(x,y) in box, F(x,y) <= T} (F positive definite binary
  form; box [0,X]^2) by integer-hull (Stern-Brocot) descent with Pick bookkeeping kept mod 2^k.
  Use it to evaluate the floor8 short-square query S(N) mod 8 for x^2+xy+6y^2 and x^2+y^2 as
  differences of region counts at T = jN, jN-1, j <= 8.
- Controls: brute force for N < 10^6; factor-backed lattice_sum from experiments/floor8 on the
  35-bit (33119086993) and 36-bit (38942974393) fixtures and on 20 random balanced semiprimes at
  each of 40, 48, 56 bits.
- Measure: descent nodes (hull vertices + Stern-Brocot steps) vs N, 30..60 bits; log-log slope.
- Mod-2^k pruning test: count nodes whose subtree can be skipped using only mod-2^k data
  (integral affine involution of the node region certified from the real region, or a
  repeated node state mod 2^k with identical subtree residue). Report the fraction pruned.
- PASS (hypothesis alive): pruned node count grows slower than N^{0.1}. KILL: slope >= 0.25
  after pruning (corrections grow like the hull vertex count).

## E2 Uniform-descent barriers (identities, verified numerically)
- E2a: half-period parabola parity sum_{y<p/2} floor(2y^2/p) mod 2 vs h(-p) mod 4 (Dirichlet
  class number formula) and Mordell's sign of ((p-1)/2)! mod p, all primes p = 3 mod 4, p < 2*10^5.
- E2b: whole-disk count sum_{n<=T} r_2(n) mod 16 vs the parity of
  #{n<=T : n = p^e m^2 or 2 p^e m^2, p = 1 mod 4, e = 1 mod 4, p not dividing m}, T <= 10^6.
- Pass = identity holds on every case; any failure is reported and the claim dropped.

## E3 Own idea: modular forms mod 2^k
- The D=-23 hidden bit is a(N) mod 2 for eta(z)eta(23z). Reduce mod 2 by Euler's pentagonal
  theorem and test whether this gives a strictly smaller instance (descent) or an instance of the
  same size (no progress). Check numerically for N < 10^5 against the cubic-splitting formula.
- Serre's density theorem as an obstruction to realising the cutoff 1/(1-q) mod 2^k.
- KILL if the reduced instance has index >= N (no shrinkage).

## Codex
- modrec_1: blind builder attempt on P1/P2/P3 (codex/modrec_1.md). Its claims get checked here.

Every claimed evaluator must be factor-blind, verified at >= 35/36 bits, all operations charged.
