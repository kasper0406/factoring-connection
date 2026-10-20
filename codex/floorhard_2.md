# floorhard_2 — ADVERSARIAL referee round

You are a hostile but fair referee. Read /Volumes/git/factoring-connection/notes/floorhard_results.md
in full (repository is read-only for you; the code is in src/floorhard/core.py, tests in
tests/test_floorhard.py; background: notes/analytic_results.md §3 Theorem L_j, and your own earlier
blind report codex/floorhard_1.out). Scratch files go under the current directory.

For EACH numbered item give a verdict ACCEPT / ACCEPT WITH FIXES / REJECT, with concrete
counterexamples or corrected statements. Run small exact computations where useful.

1. Theorem SS (parsimonious #SubsetSum → bounded quadratic congruence → floor sums): check every
   inequality (|η| < 1/(2M), separation, wrap-around, window endpoints, u_i ≥ 1). Is the claim
   "the textbook 3SAT → SubsetSum reduction with clause target 4 and slacks {1,2} is parsimonious"
   correct (check: does each satisfying assignment extend to exactly one subset)? Is #SubsetSum
   #P-complete under parsimonious reductions, and is the ETH statement (no 2^{o(√L)}) correctly derived
   (L = O(k²))? Is the bit length right?
2. Prop FP and Cor FP2 (full period ⇒ h(−M) ⇒ factoring): check the formula, the sign structure for
   M = m and M = 3m, the relation h(−4M) = h(−M)(2 − (−M/2)), and whether the decoder's "two distinct
   values average to h(−M)" can fail (e.g. if h(−q) = 0 is impossible but coincidences?).
3. Theorem TE and the prime-N variant: verify the remainder bound, the block-length condition, the
   e-range and the Baker–Harman–Pintz application (exponent of the admissible e-interval vs N).
4. Corollary DOM: is the comparison exponent 2/(j+4) right? Is the claim "L_j is informative only for
   non-uniform oracles" fair?
5. §4.1–4.2 (location of hard blocks; top blocks): (T1) arithmetic; (T2) is it true that every j = 2
   ladder block [x, x + xN^{−1/4}) lies inside the Coppersmith–Howgrave-Graham–May divisor range
   N^{β²} (β = log x/log N), with what losses (constants, N^{o(1)})? (T3) TRY HARD to embed an NP-hard
   (or factoring-hard beyond the ladder) family into TOP blocks x ∈ [√N/κ, √N] of the Taylor class,
   for some fixed j ≥ 2, possibly in a non-vertical (continued-fraction) frame. If you succeed, the
   paper's main scoping claim is wrong — say so with a construction. (T4) is the stated equivalence
   direction correct?
6. §4.3: is it right that the Theorem-SS hard family is solvable in H^{o(1)} time, so that the
   NP/#P/ETH hardness cannot exclude any T_j = N^ε? Is anything stronger (e.g. a linear-size reduction
   giving 2^{Ω(L)}) available from Manders–Adleman's actual construction? If MA's reduction has
   quasi-linear size, does that change §4.3?
7. §5 cubic claims and §6 literature table: errors, missing results, misattributions.
8. The verdict and the sharp barrier statement §7: overclaims? Is the plausibility estimate (2–5% for a
   top-block oracle with T_j = N^{o(1)} at some fixed j ≥ 4) defensible? Give your own.

Output: a structured referee report (≤ 4000 words) with verdicts and required fixes.
