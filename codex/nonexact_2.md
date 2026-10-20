# nonexact round 2 (ADVERSARIAL): referee the draft on O1 (non-exact / evaluated compressions)

You are an adversarial referee (analytic number theory, geometry of numbers, lattice algorithms, experimental
methodology). Your job is to find errors, gaps, overclaims, mislabelled statuses, and mis-scored experiments.
Be concrete: for each issue give severity (FATAL / MAJOR / MINOR), the exact location, why it is wrong, and a
fix or counterexample. Honesty over politeness. You may read files in this repository (read-only); scratch
computation goes in the current directory.

Read in full:
- notes/nonexact_results.md (the draft under review)
- notes/nonexact_prereg.md (frozen predictions; commit e4e203e, made before the runs)
- src/nonexact/core.py, src/nonexact/run_experiment.py, src/nonexact/posthoc.py, tests/test_nonexact.py
- experiments/raw/nonexact_*.jsonl (re-score every prediction P1-P10 yourself from the raw files)
- background, as needed: notes/ptebarrier_results.md (Lemma R, Theorem W, Prop NX, NX1, Lemma EX usage),
  notes/ptesplit_formulation.md, notes/ptecons_results.md, notes/smoothenum_results.md (T6-T8),
  notes/resbarrier_results.md (sec. 2.6), and your own earlier blind round codex/nonexact_1.out.

Check in particular:
1. Prop EQ / Lemma COL: is "evaluated O1 = Pi" exactly right, and is the scope of the non-exact containment
   (equality only at d = 1) stated correctly? Is Pi_beta (sec. 0 / sec. 6) a faithful residual statement?
2. Prop BRK: are the stated failure points of Lemma R / Theorem W correct and quantitatively right
   (2 log2 N vs M ~ x/log B; step 4 vacuous; J_m contains X - m)?
3. Prop AB: is the Cauchy-Schwarz coincidence count correct and is its reading ("no counting barrier") justified,
   or overclaimed? Are View B coincidences genuinely nontrivial relations?
4. Prop LAT / VA: determinant (Smith form), Henk Thm 1.5 usage with Minkowski II, value-lattice claims.
5. Prop CS: the inequality |R_M(c) - n^M/(p-1)| <= n S_max^{M-2}, the "square-root cancellation needed" reading,
   and the Burgess comparison. Prop BP-NO: is the inapplicability argument sound and correctly labelled?
6. Prop IH: verify (a)-(d) and especially the corrected distinct-relation bound
   ceil((M - 2k + 3)/(2 deg D)). Is the support count M - k + 2 right with multiplicities and common roots? Is the
   claim that common-support points are roots of E right? Is it fair to say IH "closes the gap" for PTE-type
   families with d >= k?
7. Prop NF', Prop COH and Remark PEN (rank-2 pencils have rank 1): proofs and scope.
8. Prop LIN: the model, the inequality D_s/T_s >= B^{s(s-1)/2} e^{-s^2/(2B)}, the optimisation giving
   (N/W) exp(-(1+o(1)) sqrt(2 ln(N/W) ln(2B+s))), and whether "N^{1-o(1)} per relation" is correctly scoped as a
   model statement. Any lattice strategy the model wrongly excludes?
9. Experiments: re-score P1-P10 from raw JSONL against the prereg text. Are the failures (P1/P3 in two E1 cells,
   kill signal, P8 second part, P9(e)) reported honestly, and are the post-hoc diagnoses (multiplicity from
   sampling with replacement; baseline miscalibration; sign identification) correct and properly labelled as
   post hoc? Is the (32,8,10) exhaustive count of 4 vs 0.84 handled correctly?
10. Verdict and probability: overclaims? Is "O1 is not a separate escape" justified? Missing mechanisms?

End with a summary table of issues and an overall recommendation (ACCEPT / ACCEPT WITH FIXES / MAJOR REVISION).
