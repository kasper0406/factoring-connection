# slpcirc_2: ADVERSARIAL review of notes/slpcirc_results.md

You are a hostile referee (number theory + algebraic complexity), with web search. The repository at
/Volumes/git/factoring-connection is read-only. Read notes/slpcirc_results.md in full (and, as needed,
src/slpcirc/*.py, tests/test_slpcirc.py, experiments/raw/slpcirc_E*.jsonl, notes/slpcirc_prereg.md, and your
own earlier blind report codex/scratch/slpcirc_1/report.md).

Try to BREAK every claim labelled proved. Specifically:
1. Lemma S1 and S1-k (Dirichlet + degree argument on cosets of mu_ell). Find a counterexample or a gap
   (e.g. the per-coset nonvanishing, exponent reduction mod ell vs mod n, the k-variate hyperplane count,
   Laurent exponents, coefficients in F_Q vs F_p, the case where t' < t after merging).
2. Theorem SP (per-prime dichotomy) incl. Cases 1-3, corollaries SP-1, SP-D (Dickson: the u <-> x
   transfer and the density loss), SP-K, SP-Pi (products). Is "every prime ell >= ell0 dividing |G|
   divides Lambda" really what is proved? Is ell0 correct?
3. Lemma LP (dyadic Goldfeld via Bombieri-Vinogradov): check the prime-power error term and whether BV
   applies as used (dyadic, d <= X^theta, theta < 1/2). Is the p+1 version fine?
4. Theorem SP-U and Corollaries SP-U1 (worst case, N^{1/4-eps}), SP-U2 (average case), SP-N
   (exponents polynomial in N). Check the counting (each ell > X^theta arises from <= X^{1-theta}+1
   primes), the reduction from "success on all balanced N" to |P| >= pi* - 2, the torus-vs-F_p^k
   density loss, the adaptive-trials remark, and whether the content route is correctly excluded.
   Is the claimed lower bound on SLP size (s >= log2 deg) correctly applied to Laurent/product forms?
5. Prop R, Prop DIV, Prop BB (does Buergisser 2009 Main Thm 1.2 really give tau(Per) poly => tau(k!)
   polylog, constant-free, division-free? Is the binary-search reduction valid for k! with p <= k < q?).
6. Novelty/literature: is Theorem SP / SP-U already known (Kelley 2016 Thms 2.2-2.3, CFKLLS 2000,
   Bi-Cheng-Rojas, folklore on Pollard p-1)? Give exact statements if so.
7. The verdict and the construction-table classification: anything overclaimed or mislabelled? Any
   construction in the table that is actually dense and uniform? Is the E2 H1 failure honestly reported?

Output: a numbered list of findings, each tagged FATAL / MAJOR / MINOR / OK, with a concrete fix. Then an
overall recommendation (ACCEPT / ACCEPT WITH FIXES / REJECT) per claim (S1, S1-k, SP, SP-D, LP, SP-U,
SP-U1, SP-U2, SP-N, R, DIV, BB).
