You are an adversarial expert referee (analytic number theory). Read-only: write scratch files only in your current directory.

Preprint: /Volumes/git/factoring-connection/paper2/main.tex (compiled main.pdf; bibliography refs.bib). A separate referee is
reviewing the whole paper; YOUR scope is the parts revised in commit 22e9a2f:
1. Appendix A ("Proof of Lemma 8.5'"): the general Lemma (regime (R)), the deduction of the main-text Lemma 8.5' (Section 8,
   label lem:85) from it, the cited facts (Drappeau 2015 Lemme 5 = Pascadi 2025 Lemma 3.2; Drappeau Lemme 1(i) from de la
   Bretèche–Tenenbaum; Iwaniec–Kowalski Thm 7.13; Rosser–Schoenfeld; Hildebrand–Tenenbaum Thm 3), the Perron lemma, the
   elementary lemmas (lift averaging, Möbius coprimality, factorisation n = mk, constants, density lower bound, weights over smooth
   moduli), Lemma LS and its corollary, Lemma SM, and the final proof (counting T_Q; Parts 1 and 2). Check every range, constant
   and exponent; in particular that the relative saving survives the small-conductor transfer R_1^{1-beta'} = (log x)^{O(1)}, the
   choice K = 6 + ceil(2 xi_max/a^2), and that the main-text use in Theorem 8.7 (first moment) and Corollary 8.5 still goes
   through with the new epsilon_N.
2. Whether the paper states honestly that Drappeau's Lemme 5 is published but its proof refers to Harper's unpublished preprint
   arXiv:1208.5992 §3.3, and whether any remaining text still claims or relies on Lee–Venkatesan's Lemma 8.5 or Harper's
   Proposition 2. Check the citations on the web where you can (Drappeau arXiv:1307.7554; Pascadi arXiv:2304.11696).
3. Remark rem:noN0 and FACTOR step 3 (alternating trial division with SPLIT; factor 2; algorithm uses no analytic constant).
Source note for comparison: /Volumes/git/factoring-connection/notes/lem85_results.md.

Output: numbered issues (location; FATAL/MAJOR/MINOR/PRESENTATION; description; fix), then
"VERDICT ON REVISED PARTS: CORRECT / CORRECT WITH FIXES / WRONG" and "WEAKEST LINK: ...".
