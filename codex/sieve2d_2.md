# Adversarial referee: two-variable fourth-moment bound for Theorem D' coverage

You are a referee in analytic number theory. Work read-only on /Volumes/git/factoring-connection.
Read notes/sieve2d_results.md §§1-4 (the claims under review), notes/dcltheorem_results.md
(Theorem D', authoritative; Lemmas 6-10 and §4(iii)), and the code src/sieve2d/bound.py,
tests/test_sieve2d_bound.py. You may write scratch scripts in your scratch directory and run
them with /Volumes/git/factoring-connection/.venv/bin/python (numpy, sympy available).

Check, line by line and adversarially:
1. Theorem 1 and its proof (§2.1-2.5): facts F1-F5, Lemma L (constant 4/sqrt3), Lemma C cases
   (i)/(ii) including the claims r_{v'} >= ceil(j/2) and r <= floor(j/2) under Y >= 16H^2, the
   chord/line count in case (ii), the summations A_j..F_j (are any configurations double-
   counted or missed? is "counted at its own shortest vector" sound?), the falling-factorial
   expansion and the set D. Try to break each inequality numerically on small cases.
2. Corollary 1 and Theorem 2 (scale substitution; hypotheses H <= N^{1/15}/C0; the n^2 loss).
3. Theorem 3: the claim that D''s unchanged algorithm with eps' <= 1/80 fails on O(N^{-2eps'})
   — recheck E_a, E_b, E_c exponents with T = floor(2^{(1/4-eps')n}/n), and all scale
   separations D' needs.
4. Proposition 4 (prime presieve, large sieve) and Proposition 5 (fake distribution via
   Peled-Yadin-Yehudayoff 2011, arXiv 0801.0059; check the exact theorem used, the odd-k
   reduction, the complement/thinning step, and the hypothesis N = n p(1-p) - 1 -> infinity).
5. §4.2-4.4 barrier statements: are any "proved" labels overstated?

Output: for each item CORRECT / GAP (with exact fix) / WRONG (with counterexample), then an overall
verdict (ACCEPT / ACCEPT WITH FIXES / REJECT). Be concise.
