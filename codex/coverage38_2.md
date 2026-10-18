You are an adversarial referee (analytic number theory, geometry of numbers, computational complexity). Read-only access to /Volumes/git/factoring-connection; scratch scripts go in your current working directory.

Review notes/coverage38_results.md (agent coverage38), with code in src/coverage38/ and tests in tests/test_coverage38.py. Background: notes/dchworst_results.md (Theorems W, R; barriers B1-B4) and notes/below16_results.md §1 (cost model M+).

Check, and try hard to BREAK, each claimed result. For each give a verdict CORRECT / CORRECT WITH FIXES / WRONG, the exact fix, and any counterexample (run code where useful):

1. Theorem U (§2): the Cauchy-Schwarz/height-count/Holder proof of X >= (c1 P/(2 sqrt 3))^{2/5} lam^{4/5} J^{-1/5}; the window hypothesis (WIN) and whether it really holds for full windows and Lehman windows (constants, signed weights, the range of S = Bx + AN/x, roots per (v, M)); whether grouping by (v, M) and charging T >= #groups is legitimate in M+; Corollaries U1-U4 (in particular U2's claim that for theta >= 1/4 no design beats 1/5, and U3's 1/6 floor). Is the "single band is extremal" claim correct?
2. The composite-modulus remark: does Theorem R (dchworst §3) really extend to squarefree moduli with prime factors = 3 mod 4 and > H? Is "Harman's sieve is not needed" fair?
3. Proposition L (§3.2) and complements (a)-(c): constants (Hoeffding, union bound, distinctness), the mean-shift bound, the claim that pointwise exponential/character-sum information cannot certify coverage at lam < J^{1/4}, and the Bernstein remark. Is the "Consequence" paragraph overstated?
4. §3.3 literature table (Baker 2012, Bourgain 2005, Fouvry-Shparlinski, Korolev 2016/2019, Bourgain-Garaev 2014 Thm 16, Irving, small products in residue classes / Garaev's remark). Verify citations and ranges with web search; fix anything wrong.
5. Proposition D and the dichotomy (§4), Proposition P (§5), and the conditional table (§6). Is any claim false or overstated? Is there any standard hypothesis we missed that WOULD imply coverage at lam < J^{1/4} (or DCH_pair)?
6. Experiments §7 vs notes/coverage38_prereg.md: were the frozen rules applied correctly (experiments/raw/coverage38_*.jsonl, src/coverage38/analyze.py)?
7. Novelty/prior art: is Theorem U or Prop L already known in some form (e.g. lower bounds for BSGS-type factoring, discrepancy/oracle barriers)?

Output Markdown: overall verdict first (ACCEPT / ACCEPT WITH FIXES / REJECT), then numbered findings with severity (major/minor), then a list of required fixes.
