# smoothenum_1 — blind round: output-sensitive enumeration of doubly-smooth NFS pairs

You are a research number theorist / algorithms researcher. You have read-only access to the repo at
/Volumes/git/factoring-connection (scratch work goes in your cwd; you may run python there). Web search is
enabled; use it for literature. Honesty matters more than a positive verdict. Label every statement KNOWN
(with citation), PROVED (give the proof), HEURISTIC (name the exact assumption), or SPECULATIVE.

## Setting
NFS for a general integer N, heuristic L[1/3] model. Notation L[c] = L_N[1/3, c], x = log N.
Base-m polynomial pair: g(X) = X − m, f of degree d = δ (x/log x)^{1/3}, coefficients ≈ m ≈ N^{1/(d+1)}.
Sieve box |a|, |b| ≤ E = L[e]. Rational norm |a − bm| ≈ E·N^{1/d} (size coefficient 1/δ in units of
x^{2/3}(log x)^{1/3}); algebraic norm |b^d f(a/b)| ≈ E^d N^{1/d} (size δe + 1/δ). Factor bases B = L[β]
(rational) and L[β_a] (algebraic). Smoothness of a size-s number at L[β] has probability L[−s/(3β)].
Plain sieving costs E² = L[2e] (the whole box); sparse LA costs (#columns)^2. GNFS = (64/9)^{1/3} = 1.923;
Coppersmith MNFS = 1.9019; Coppersmith's factorization factory reaches 1.6386 per N after a shared
precomputation L[2.0069] that tabulates the rational-smooth pairs (a, b) once for many N sharing m.

## Question
Is there an algorithm that, given (f, g) and the box, outputs most pairs (a, b) with BOTH norms smooth (or
just the rational-smooth pairs) at cost much less than the box area E², ideally ≈ (number of output pairs)
× L^{o(1)}? If so, what single-N L[1/3] constant results? If not, what is the precise barrier?

## Required sections
1. **Target.** Write the MNFS / GNFS optimisation program with a generic enumeration-cost term. Compute:
   (i) the constant if the rational-smooth pairs can be listed at output cost E²·L[−r] (r = rational
   smoothness penalty) — re-derive 1.6386 — and the parameters at that optimum; (ii) the constant if only
   the ALGEBRAIC-smooth pairs can be listed output-sensitively; (iii) the constant if the doubly-smooth pairs
   can be listed output-sensitively; (iv) a parametric family: if listing the rational-smooth pairs costs
   E²·L[−θ r] for 0 ≤ θ ≤ 1, give c(θ) and the minimal θ that beats 1.9019 by a non-negligible margin.
   Numerically minimise explicit programs (show the code).
2. **Mechanisms**, each with a heuristic cost analysis and the effect on the constant:
   (a) bottom-up product enumeration: generate B-smooth integers s of size E·N^{1/d} as products of
   factor-base primes, solve a − bm = s for (a, b) in the box (a linear Diophantine condition: given b,
   a = s + bm, so each s yields ≈ 2E/m... points — compute honestly; note |a−bm| ≤ E·m means for each s
   only O(1) or fewer b's give |a| ≤ E), then test the algebraic side. Honest count vs sieving.
   (b) multi-special-q (special-(q₁, q₂), double special-q on both sides, nested special-q lattices).
   (c) lattice/CRT enumeration of pairs whose norm is divisible by a chosen smooth modulus M of size close
   to the norm (Coppersmith/Howgrave-Graham smooth-multiple small roots; "lattice sieving with a huge
   smooth special-q"; enumerating short vectors of the lattice {(a,b): M | a − bm}). Does it reach
   output-sensitivity? What fraction of the smooth pairs does a single M capture, and how many M are needed?
   (d) factory transfer to a single N: manufacture many related numbers (kN for many k, many m's, many
   rational polynomials) so one rational-side precomputation serves one N. Cost it.
   (e) anything else: sum-product structure, subset-sum / k-list (Wagner) generalized birthday on
   logarithms of primes, Schroeppel–Shamir style, quantum-free tricks, smooth-number generators.
3. **Barriers.** For each failing mechanism a precise statement in the heuristic model. Is there an
   information-theoretic or query-complexity floor (e.g. any algorithm whose only access to smoothness is
   via divisibility/lattice queries needs ≥ X operations)? Is "listing smooth values of a − bm in a box
   below box cost" equivalent to some known hard/open problem? Discuss the trivial lower bound
   (#output) versus the E² sieving cost, and where the true complexity could lie.
4. **Literature.** Is output-sensitive smooth enumeration known? Check: Bernstein–Lange "Batch NFS"
   (2014/921) and the factory; Coppersmith 1993 "Modifications to the number field sieve"; Franke–Kleinjung
   lattice sieving and special-q variants (incl. "double large prime special-q", bucket sieving); Kleinjung
   "special-q with q on both sides"; TNFS / exTNFS; Kirchner 2018; Boudot et al. 2020 (RSA-240); any work on
   enumerating smooth values of polynomials, smooth numbers in short intervals/APs algorithms, Bernstein's
   "how to find smooth parts of integers", "Sieving in function fields", Joux–Lercier, Pomerance's
   comments on whether sieving is optimal. Report precisely what is and is not known.
5. **Best candidate** with full derivation and your own adversarial attack; else rank near misses.

Keep it under ~4500 words. Be numerically careful.
