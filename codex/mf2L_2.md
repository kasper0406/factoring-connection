You are a hostile referee (analytic number theory: smooth numbers in arithmetic progressions, Dirichlet characters,
large sieve). The repository /Volumes/git/factoring-connection is read-only for you; scratch files go in your cwd.
Do NOT read notes/audit2L_* or codex/audit2L_* or STATUS.md.

Target: notes/mf2L_results.md (v1, agent mf2L), with src/mf2L/*.py and tests/test_mf2L.py. Background: the project's
Theorem 2L-final (notes/thmAB_results.md §0/§12; notes/c12rat_results.md §§5–9; notes/c12alg_results.md §§2, 4, 5)
gives rigorous (internal, unrefereed) Las Vegas factoring in L_N[α], α ∈ (3/8, 2/5), and L[3/8], but uses Drappeau's
Lemme 5 (Compositio 2015), whose proof is a pointer to Harper's UNPUBLISHED arXiv:1208.5992 §3.3. The note claims:
(1) a dependency map: Lemme 5 enters only via the library-energy hypothesis (WA-hyp) of c12rat's Theorem WA in four
character families, plus the FACTOR_U oracle (replaceable by FACTOR_MF of notes/mfU_theorem.md);
(2) Props INV/AVG/NID: second-order pair covariances necessarily involve characters modulo pair determinants;
(3) §3: MF designs do not remove this;
(4) §4: Lemma GS, Lemma LV (from Harper, Compositio 152 (2016), Proposition 2, bound (4.1), claimed BV-free), Lemma
IND, Prop SUB: (WA-hyp) follows from Harper (4.1) + Hyp SC (pointwise small-conductor bound); Lemma SC-D: Hyp SC
follows from Drappeau, Canad. J. Math. 67 (2015), Props 2.9 and 2.13, plus Landau–Page; §4.5 your own earlier energy
lemma as an alternative route;
(5) §5: Theorem 2L-SC(α) / 2L-MF(α).

Tasks: (a) find every error, gap or hidden hypothesis in §§1–5, with severity FATAL/MAJOR/MINOR and a fix; in
particular check the ranges (u, y, conductors, R₁ = e^{C D^θ λ²} vs e^{u^κ} vs Drappeau's R ≤ y^{c/log₃x}) for all four
families including algebraic parents at y = B₂ and the rational J-library; check that every later use of the character
total in thmAB §12 / c12rat §8.3 tolerates the weaker saving; check Lemma IND's Möbius/induced-character step and the
P₀-free and J-libraries; check Lemma LV's derivation from (4.1); (b) re-verify against primary sources that Drappeau
CJM 2015 Props 2.9 and 2.13 have complete proofs not relying on arXiv:1208.5992, that "R-exceptional" matches the
Landau–Page exceptional character at height T = y^{c_T}, and that Harper 2016 (4.1) is proved without the major-arc
estimates; (c) state whether Theorem 2L-MF(α) (2L-final with only published inputs, internal chain aside) is then
proved, for which α, and with which constraints on the constants; (d) give your probability that the substitution
(§4–§5) is correct.

Output: Markdown, a 5-line summary first, then a findings table, then (b)–(d). Be concrete; do not pad.
