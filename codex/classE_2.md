# classE_2: ADVERSARIAL review of notes/classE_results.md (v1)

You are a hostile but fair referee (analytic number theory, continued fractions, exponential sums, lattice-point
geometry, statistics). Read /Volumes/git/factoring-connection/notes/classE_results.md (v1) and
notes/classE_prereg.md; code in src/classE/*.py and tests/test_classE.py; raw/analysis data in
experiments/raw/classE_* and experiments/analysis/classE_*.json. Your own blind report is codex/classE_1_report.md
(we credit it as "cx1"). The repo is READ-ONLY; scratch only in the current directory (python:
/Volumes/git/factoring-connection/.venv/bin/python, run with PYTHONPATH=/Volumes/git/factoring-connection).

Attack every claim. In particular:
1. Prop SPLIT (§1.3): the identity, the product bound s_k r_{k+1} >= N/(a_{k+1}+2), and the interpretation "split
   discriminant N^2 => only one character bit (a|p); class E is a CFRAC with relations of size ~N". Is the
   interpretation overstated? Is there any OTHER computable quadratic or higher-character constraint on path data?
2. Theorem F-SQ (§2.2): check the reduction to a step function of y = sqrt(N)/m0 mod 1 (Farey-of-order-2Q
   breakpoints, the boundary layer 1/r, the t-jumps), the double large sieve application (centering, spacing counts,
   (1+XY), the normalisation by pi(I)^2, the conditioning on N = 1 mod 4 and p1 != p2), the Erdos-Turan step and the
   final exponent m0 Q^2 N^{-1/6} log N; and the tb variant (trigonometric approximation of [log p1 - log p2 in J]).
   Is the claimed coverage "path of floor(sqrt N)/N up to convergent denominators N^{7/12-delta}" right? Give the
   correct statement if not. Exact reference/lemma numbers for the double large sieve.
3. §4.2 Heegner/Duke/Michel/Harcos-Michel claims: are the cited results correct for non-fundamental D = -4kN, for
   genus cosets, with level structure? Exact references. Is "Q-leaks are cusp excursions" accurate?
4. §3: FD3-0, the use of your Thm 4, the DSF analogue, the spatial vs frequency split remark, the toy test.
5. §5 battery: registered-rule compliance, deviation D1 (corrected planted expectation 2 rho sqrt(sum S_yy)) --
   is the correction right, and is it legitimate to pass V1 post hoc? Is the tb floor (top 12 bits) adequate? Any
   statistic whose null is uninformative (degenerate-in-disguise)? Recompute a few z from raw data.
6. Verdict wording (§0): overclaims, missing caveats, wrong attributions/novelty.

Output: a markdown referee report: overall recommendation (ACCEPT / ACCEPT WITH FIXES / MAJOR REVISION), then
numbered findings tagged FATAL/MAJOR/MINOR/NIT with exact fixes, then anything you verified correct.
