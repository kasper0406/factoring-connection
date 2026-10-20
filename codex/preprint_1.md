You are an adversarial expert referee (analytic number theory, computational number theory, algorithms) for a
mathematics preprint. Be skeptical: the paper claims the first rigorous improvement over Lenstra–Pomerance (1992)
for factoring every integer, so assume errors exist and hunt for them.

Files (read-only; you may write scratch files only in your current working directory):
- The preprint: /Volumes/git/factoring-connection/paper2/main.tex (compiled: paper2/main.pdf), bibliography paper2/refs.bib.
- Source notes it was distilled from (consult them only to check for discrepancies or omitted hypotheses; the PAPER is
  what you are refereeing): /Volumes/git/factoring-connection/notes/{thmA_proof,thmA_audit,lemmaM_results,hypD_results,
  pcest_results,lv6_results,ldbound_results}.md.
- Primary sources you should consult on the web where possible: Lee–Venkatesan arXiv:1805.08873 (J. Number Theory 187, 2018),
  Harper arXiv:1208.5992, Lenstra–Pomerance JAMS 5 (1992), Hildebrand–Tenenbaum Trans. AMS 296 (1986), Lagarias–Odlyzko 1977,
  Stark Invent. Math. 23 (1974), Montgomery–Vaughan Mathematika 20 (1973).

Claim: for every admissible fixed (d, k, mu, eta0) a Las Vegas algorithm outputs the complete factorisation of any N in
expected time L_N[1/2, (1+theta)/sqrt2 + o(1)]; hence for every c > 1/sqrt2 expected time O(L_N[1/2, c]).

Please referee the WHOLE paper, section by section, with full rigor. In particular check:
1. Statements: is the main theorem stated honestly (infimum over fixed parameters vs o(1); complete factorisation; prime
   powers; no hypotheses)? Are the numbers in Table 1 correct (recompute exactly)?
2. Section 2 cited facts: are they quoted correctly with the right ranges (CEP, Hildebrand–Tenenbaum Thm 3, Tenenbaum's
   Psi_q, LP main theorem, large sieve, Brun–Titchmarsh)? Is Lemma 2.4 (local behaviour of rho) correctly derived?
3. Section 3: parameter/size lemmas (Lemma 3.1, 3.2, 3.4 sizes with constants), Lemma 3.5 homogeneity, Lemma 3.6
   regularity (constant 200), Lemma 3.7 coprimality, the algorithm (does it match the analysis? any hidden use of r?).
4. Sections 5–7: Theorems S, S' (constants, the 2Y'^2<N shortest-vector argument), Theorem DD (measurability formulation),
   Theorem K and Proposition 7.4 (field selection with factor exits; Wald; harmonic mean; is any conditioning illegitimate?).
5. Section 8: Lemma 8.1 (one-point smoothing) is NEW in this paper and replaces LV Lemmas 5.13/5.18/5.20 — check the kernel
   basis claim, the representation of every |m| <= ||v||_2 with ||x||_inf <= 2V'+1, the total-variation bound, and the
   log-concavity/unimodality/layer-cake argument with great care. Then Lemma 8.3 (good modulus => smooth density), Lemma 8.4
   statement, Corollary 8.5, Lemma 8.6 (fibres, Moebius), Theorem 8.7 (first moment) — check every power of N, L, log N.
6. Section 9: Lemma 9.1 (pairs), Lemma 9.2 (kernel transfer: index, kernel, fibre bound, wedge minor), Lemma 9.3 (Brun–
   Titchmarsh spreading), Lemma 9.4 (divisor switching), Theorem 9.5 (L^2 family bound: classes A/B/C and every exponent
   margin), Theorem 9.6 (second moment: the order of operations, AM-GM step, enlargement, w-average, fibre DS, partner count,
   power count V^4/Y^2), Theorem 9.7 (Paley–Zygmund), Proposition 9.8.
7. Section 10: Lemma IR, rejection events, row identity with signs, Lemmas C, H, discriminant bound, Stark/Landau/LO usage,
   Lemma E (|E|<=1), Lemma X (1/(5d)), SQ algorithm and Lemma SQ, kernel sampling, Proposition 10.15.
8. Sections 11–13: batch smoothness, cost optimisation, uniformity and recursion, dovetailing for small N.
9. Appendix A: Lemma 8.4's proof (the TODO is known; check that Parts 1 and 2 are otherwise right, ranges, eta choice).
10. Attributions and prior art: anything misattributed, overclaimed, or missing (e.g. an existing rigorous c<1 result, or a
    known barrier)? Are bibliography entries correct?

Output format (plain text):
- For each issue: number, location (section/lemma/equation), severity FATAL (breaks the theorem) / MAJOR (needs a new
  argument) / MINOR (local fix) / PRESENTATION, a precise description, and a concrete fix.
- Then a per-section verdict (CORRECT / CORRECT WITH FIXES / WRONG).
- End with: WHOLE-PAPER VERDICT: one of ACCEPT / ACCEPT WITH MINOR FIXES / MAJOR REVISION / REJECT, and WEAKEST LINK: one line.
