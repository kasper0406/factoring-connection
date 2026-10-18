# Blind + adversarial round: an explicit ECM family as a prefilter for deterministic average-case factoring

You are a research mathematician (analytic number theory, arithmetic of elliptic curves) and an
adversarial referee. Work read-only on /Volumes/git/factoring-connection; you may write scratch
scripts in your scratch directory and run them with
/Volumes/git/factoring-connection/.venv/bin/python (sympy, gmpy2, python-flint, numpy).

Background (read these first):
- codex/r18_ideas.out §5 (the proposal under review; written by a previous GPT run).
- notes/avgtime_results.md (Theorem AV: ladder of D' runs, average N^{1/6} polylog over S_n).
- notes/sieve2d_results.md §3 (Theorem 2: D' coverage exceptional set <= C n^4/H^4, valid for
  2 <= H <= N^{1/15}/C0).
- notes/dchworst_results.md §3 (Theorem R: DCL-wc factors EVERY pair covered by some block; cost
  N^{1/6+o(1)} + list cost, i.e. ~N^{1/6} H^{5/3} for D's design at height H).
S_n = {(p,q): p<q<=2p primes, 2^{n-1} <= pq < 2^n}, uniform. Costs are bit operations, N^x means
up to N^{o(1)}. The algorithm must be deterministic and complete (always factor), and see only N.

The proposal: run the first J curves (with marked points) of a fixed explicit family (e.g.
Suyama sigma = 6,7,8,...) with stage-1 bound B, then a fallback ladder on failures. With
Z_N = #{j <= J : [M_B]P_j vanishes at exactly one of p, q}, the hypothesis (EC2) is
E Z_N >= cJ and Var Z_N <= J N^{o(1)} over S_n; Chebyshev gives failure fraction <~ 1/J, and
codex claimed a complete average cost N^{23/152} with B = N^{1/8}, J = N^{1/38}.

Tasks (be blind: derive things yourself, do not trust §5):
1. Verify or refute the fallback bound (*): E[1_A T_fallback] << N^{1/6} delta^{7/12} + N^{2/15}
   for ANY subset A of S_n of density delta (no independence). Which ladder (AV's D'+detector,
   or Theorem R's DCL-wc) should be used, where must it start, what are the reach sets, and what
   is the best second term once Theorem R replaces D' (no E_a/E_c sets)? Give the sharpest
   statement you can justify.
2. Re-optimise the whole design under (EC2). Is 23/152 optimal? Consider: early stopping
   (run curves sequentially, stop at the first proper divisor), letting B = N^b with b small
   (c = c(b) constant), and the range of J for which (EC2) is assumed. What does a weaker
   hypothesis give (Var <= J^{2-eta}; failure <= J^{-kappa})? What does an ECM stage 2 (BSGS/FFT
   continuation to B2 ~ B^2) change, at the level of exponents?
3. Is (EC2) as stated plausible at all for J = N^{theta}, theta > 0 fixed? Think about
   covariances of the events "ord of P_j mod p is B-smooth" across j, j' over random p,
   including any "common-mode" effects shared by all curves (local data at small primes, the
   location of the Hasse interval, size of p). If (EC2) is false or implausible, state the
   correct hypothesis that the route actually needs (e.g. a single-prime hitting statement)
   and say whether a gcd-staging trick reduces separation to a one-prime statement.
4. Literature (use web search): what is PROVED about smoothness of #E(F_p) or of point orders
   (i) for a fixed curve as p varies, (ii) on average over boxes of curves y^2=x^3+ax+b,
   |a|<=A, |b|<=B, as p varies (Banks-Shparlinski; David-Koukoulopoulos-Smith; Martin-Pollack-
   Smith; David-Smith; others), with the exact box-size ranges, (iii) for smooth numbers in
   Hasse intervals for almost all primes (Matomaki-Radziwill and successors; Pomerance;
   Lenstra's ECM analysis), (iv) conditional results (GRH, EH for number fields:
   Barbulescu-Jouve). Give exact theorem numbers and box/range hypotheses. Which, if any,
   restricted covariance or hitting statement for an EXPLICIT family of size J << p is
   provable (unconditionally, or under GRH/EH)? An almost-all bound with any power saving in J
   would be valuable; so would a proof that such a statement is out of reach of known methods.
5. Adversarial: list every error, gap or overclaim in r18 §5.

Output: numbered sections 1-5; label every claim PROVED / CONDITIONAL (on what) / HEURISTIC /
REFUTED / UNVERIFIED. Exact rational exponents. Be concise (<= 2500 words).
