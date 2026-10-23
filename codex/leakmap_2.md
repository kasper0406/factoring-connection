You are an adversarial referee (research mathematician / theoretical computer scientist). Repo:
/Volumes/git/factoring-connection (read-only for you; scratch files in your current directory). Your job is to find
errors, gaps, overclaims and missing prior art in notes/leakmap_results.md (agent `leakmap`, phase 24). Rigour and
honesty matter more than a positive result. Use web search to check literature claims; give exact citations and mark
anything unverifiable as UNVERIFIED.

Read: notes/leakmap_results.md (the target), notes/leakmap_prereg.md, codex/leakmap_1.out (your own blind round; the
note credits several results to it), src/leakmap/*.py, tests/test_leakmap.py, and as background
notes/factorbits_results.md §3.2 (Theorem PT) and notes/slpcirc_results.md §4.2 (Theorem SP).

Check in particular, with explicit counterexamples where possible:
1. Theorem PT-X: the inequality Pr[correct or split] <= B + 2s + 2 Pr[D and no split]; especially the step
   s0 <= s + Pr[D and no split], the finite distinct-pair version, and whether "any procedure" for c_v is really fine.
2. Lemma NE-C / Theorem NE1 / Cor NE1-bal: the formal reference with base-zero components; the claim that the first
   deviation lies on kappa(x)'s path; the use of the "largest component order"; mu bounds for affine and dyadic
   exponents (gcd(alpha, o) classes; the 2^j term); the factor 4 in the final bound; the small-order mass estimates
   (random base: L tau(x-1)/(x-1); fixed base: L^2 log|y|/log X); whether the classification SMALL/SUBGROUP/UNIVERSAL
   into families M/M/A is honest; the claim that (Q) does not occur in NE1 (Euler criterion a^{(N-1)/2}?).
3. Theorem NE1-R as restated from your blind round; §3.4 (several gates via Kelley/SP): is
   E_j(xq) = alpha_j q + beta_j mod (x-1) right, and the 3 lambda ln X bound?
4. §4: Lemma GAR (rep(x)/N = (x_p - x_q) theta + x_q/N mod 1, theta = j/q + 1/N) and the claim that 3n/4 top bits of
   theta factor N via Howgrave-Graham's partial approximate common divisor algorithm (check the exponent; is it
   3n/4 or something else?).
5. §6 Theorem AC0-eps: re-check your own proof as it is restated (character pair, minor/major arcs, Landau-Page,
   density b^2, LMN). Is anything lost in the restatement? Is the "+O_C(n^{-C})" claim right? Does the major-arc
   partial summation over a rectangle (p <= u, q <= v) with the constraint pq-dependence work (e(beta pq) is not a
   product of a function of p and of q... or is it?). Is there prior art (e.g. results on AC0 and quadratic
   residuosity / factor bits, Shparlinski's book)?
6. Theorem PT-J and Cor PT-J1 (a) unconditional range log|c| <= c0 sqrt(log X): conductor of (c|.) vs Landau-Page
   range; (b) GRH range log|c| <= X^{1/4 - delta}: is the GRH bound for sums over primes of a character with
   conductor f really O(sqrt(X) log^2(fX))? 
7. §5 / E2: is the inference "non-quadratic characters do not leak, quadratic do" sound given the populations
   (NG2 conditions on p = +-1 mod 5; NG1 independence from all quadratic symbols)? Any confound in strata
   (N mod 2^8 m)? Is the post-hoc A1 explanation (square forms that split) convincing?
8. §8 Theorem CLASS and "strongest open model" (PRIM): overclaims? Is there a known result that already classifies
   or breaks PRIM-type models (e.g. generic ring model with a power map; Aggarwal-Maurer; Jager-Schwenk;
   Boneh-Venkatesan; Maurer-Wolf)?
9. Prior art for PT-X/NE1 (generic ring model with N-dependent exponents; "generic RSA" where e is an input).

Output: a numbered list of findings, each tagged WRONG / GAP / OVERCLAIM / MINOR / OK, with the exact fix, then an
overall verdict (ACCEPT / ACCEPT WITH FIXES / MAJOR REVISION). Keep it under ~4000 words.
