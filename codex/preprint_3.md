You are an adversarial expert referee (computational number theory, algorithms, analytic number theory). Read-only: write
scratch files only in your current directory.

Preprint: /Volumes/git/factoring-connection/paper2/main.tex (compiled main.pdf, bibliography refs.bib). Earlier referees have
covered Sections 1–14 and Appendix A. YOUR scope (commit c8e1ee5):
1. Section 15, "The recursive shape test": Lemma (main proof at general a), Lemma OS (oracle substitution, applied to the core
   inside the strong induction), Wald for adaptive sub-calls, Lemma T_rec, Lemma F (solving F(x) <= L^p + L^delta F(sqrt(2x))),
   Theorem RA (c(a) = max(a + theta/(2a), (2+sqrt2)a), c* = sqrt((1+sqrt2) theta)), Remark RA (constants, sanity comparison with
   heuristic NFS 2/sqrt(d+1) and with the random-squares folklore, fixed-depth recursion, prior art), Corollary P (restart schedule;
   parameter sequence d_j = 2j+7, k_j = ceil(d_j^2/4)+j, eta0_j = 1/d_j^2; races against trial division so no N0(j) is needed; "no
   effective rate supplied"). Check every step, including that (A1)–(A3) really hold at a general fixed a, and the restart computation.
2. Section 16, "Application: discrete logarithms in F_p^*": Theorem DL (c_DL = min_a max(a+1/sqrt2, a+theta/(2a), 3a)), Corollary
   (arbitrary base), the literature paragraph (Pomerance 1987; Hafner–McCurley; Biasse et al.; Vollmer 2000 claim retracted in
   ANTS V 2002 corrigendum — check this on the web if possible; Enge–Gaudry; Kleinjung–Wesolowski; Bender–Pomerance; de Boer–
   Pellet-Mary–Wesolowski), Algorithm DLOG, Lemma (transfer with N := p), the row identity, Lemma (ideal bookkeeping), Lemma (units,
   reduction, Lambda_2, torsion), Lemma (unit elimination via Log + LLL + Blanksby–Montgomery; precision and LLL bounds), Lemma
   (hidden exponents), Proposition XK (EXACT_DEPENDENCE; rank profile mod q; Dixon lifting; cost L^{3a}), proof of Theorem DL, Remark
   (interaction with recursion: c_DL = c* with RA; 3 sqrt(theta)/2 with the portfolio), Remark (prior art; what does not transfer).
3. Consistency of the abstract/intro statements about Sections 15–16 with their bodies; honesty of status remarks.
Source notes for comparison: /Volumes/git/factoring-connection/notes/{recurse_results,auditRA_audit,dlog_results}.md.

Output: numbered issues (location; FATAL/MAJOR/MINOR/PRESENTATION; description; concrete fix), then
"VERDICT SECTION 15: ..." and "VERDICT SECTION 16: ..." (CORRECT / CORRECT WITH FIXES / WRONG) and "WEAKEST LINK: ...".
