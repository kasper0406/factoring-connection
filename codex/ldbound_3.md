# ldbound_3: final adversarial re-check of notes/ldbound_results.md (revised version)

You are a hostile referee. The repository /Volumes/git/factoring-connection is read-only for you.
You reviewed an earlier draft as ldbound_2 (your report: /Volumes/git/factoring-connection/codex/ldbound_2.out).
The note has been revised to apply your fixes, two Opus critics' fixes, and your suggestion to let the
base exponent mu (M = N^mu) vary. Read the revised /Volumes/git/factoring-connection/notes/ldbound_results.md
in full (and pcest_results.md as needed; the blind-round report notes/ldbound_codex1_report.md is context).

Check specifically:
1. Every one of your ldbound_2 fixes: applied correctly? Anything still wrong?
2. The mu-generalization: with M = N^mu, mu in [nu + eta0, 1), h >= max(2nu, mu - (d-2)nu, 1 - d mu + nu) + eta0,
   W = N^{mu + nu + eta0/2}: re-verify Lemma GB_w's bound, the sizes paragraph, FM (all steps, Lemma 8.5' ranges,
   fibre counts needing M/V >> N^{o(1)}), Theorem R, Lemma Kern, Lemma 1.3, Lemma W, Theorem L2_fam (every
   margin; does anything secretly need mu = 1/(d+1)?), LD* Steps 0-7, Theorem A'''. Does any pcest result
   use M = N^{1/(d+1)} in a way that breaks for other mu (e.g. Lemma 1.3, irregular-lattice handling,
   Lemma GB's |a_0| bound, Theorem R's covering, the sign filter)?
3. The constants table and the best-mu formula mu* = (1 + (d-1)nu)/(d+1) for k >= (d^2+1)/4; min k values
   (d = 7, 9, 11, 13, 19: 21, 13, 13, 14, 17); A' values. (You may run
   `.venv/bin/python -m src.ldbound.params` from the repo root, read-only.)
4. Props 7.1-7.4 as now stated.
5. Any remaining gap that would prevent the statement "Theorem LD* is proved for (a'''_k; mu), so Theorem A
   holds for this design conditional on LV6 (+ listed caveats)".
Give per-item verdicts CORRECT / CORRECT WITH FIXES / WRONG and a final overall verdict.
