# Blind attempt: two-variable sieve for prime pairs avoiding Farey classes in a prime band

You are a research mathematician (analytic number theory, sieve methods). Work read-only; you may read files under /home/knielsen/code/factoring-connection (in particular notes/dcltheorem_results.md, Lemma 9 and section 4(iii)) for context, and use web search for literature. Do NOT look at notes/sieve2d_* (blind attempt).

## Problem

Parameters: n large, P = 2^{n/2}, H >= 2 (think H = n^{1+a} or H = N^{eps/4}, N = 2^n), Y = P^{1/3} H^{-1/3} (so P^2 = Y^6 H^2 up to constants).
- V = {(A,B) in [1,H]^2 : gcd(A,B)=1, A != B}, |V| ~ (6/pi^2) H^2.
- Band: 𝒫 = {primes l = 3 mod 4, Y < l <= 2Y}. Since Y > H^2, the |V| residues A B^{-1} mod l are distinct.
- lambda := |V| * sum_{l in 𝒫} 1/l  ~ c H^2 / n.
- S_n = {(p,q) primes: p < q <= 2p, 2^{n-1} <= pq < 2^n}, |S_n| ~ c0 P^2/n^2.
- (p,q) is *uncovered* if no (v,l) in V x 𝒫 has l | Bp - Aq.

Known (the paper's Lemma 9): fixing q and sieving p in an interval by the arithmetic large sieve gives #{uncovered p} <= (P + Q^2)/J with J >= lambda, so the uncovered fraction of S_n is O(n/lambda) = O(n^2/H^2). The Poisson heuristic predicts exp(-lambda) (events l | Bp - Aq nearly independent over l).

## Tasks

1. Derive the strongest bound you can PROVE on #{(p,q) in S_n uncovered}/|S_n| as a function of lambda, H, n (full proof, all error terms explicit). Consider: a two-variable moment method / Bonferroni / Selberg upper-bound sieve on the integer lattice (p,q) in the box (the joint conditions l_1...l_j | B_r p - A_r q define a sublattice of Z^2 of index prod l), lattice point counting with smooth weights (Poisson summation and the dual lattice), handling of degenerate lattices with very short vectors (e.g. all l with the same v: the line Bp = Aq), Selberg presieve to restrict to primes, bilinear / BV / Barban–Davenport–Halberstam inputs, Janson-type inequalities.
2. Identify the level-of-distribution barrier: what is the largest number j of band primes whose joint conditions can be equidistributed in the box with the needed accuracy, and what bound on the uncovered fraction does that imply (e.g. lambda^{-m} for some m)? Is exp(-c lambda), or anything better than a fixed power of 1/lambda, provable by any known method here? If not, say exactly why (accuracy needed vs available) and whether a construction of a "fake" distribution shows the barrier is sharp for level-based methods.
3. Decide whether restricting to primes vs. integers matters (only polylog factors?).
4. Translate into the factoring trade-off: the algorithm costs N^{1/6} H^{5/3} polylog; with the old bound, N^{-eps} exceptions cost N^{5 eps/6}. What does your bound give?

Output: a precise theorem statement, a complete proof (or the precise point where it fails), the barrier analysis, and a list of any literature that already contains this two-variable/lattice moment argument (with exact references). Label every claim proved / heuristic / conjectured.
