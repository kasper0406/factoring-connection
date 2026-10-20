# snfsrig_2 — adversarial referee

You are an adversarial referee (analytic number theory, NFS). The repository /Volumes/git/factoring-connection is
read-only; write scratch only in your current directory. Read notes/snfsrig_results.md (the note under review) and, as
needed, its base texts: notes/lv5unif_results.md (Theorem U), notes/recurse_results.md §3, notes/nested_results.md §§3–8,
notes/ldbound_results.md, notes/lemmaM_results.md (HECM barrier (c)), and codex/snfsrig_1.out (an earlier blind round).
Code: src/snfsrig/*.py, tests/test_snfsrig.py.

Check, with the aim of finding errors:
1. Lemma RES, HR, PH, RS (proofs, hypotheses, constants). Is Prop DI's table right (T, c1, c2, fields)?
2. Prop UF (universal one-level floor via RES; the (1+p)/(3p+2) formula) — is the model argument correct, and is it
   really "universal" for 2D sieving with one linear side? Does anything (skewness, non-monic, projective roots,
   several rational sides) evade RES's θ_coeff ≥ 1/d?
3. Prop R2-S: the rational-only model computation, the α = 1/3 example (λ = 11/64, margin 5/6), the claim that below
   1/3 targets expand, and the "Why GNFS cannot do this" dilemma paragraph.
4. Prop REC: the recognition floor (HECM 2(1−α)/3 ≤ α ⇒ 2/5; Theorem U recursion; Lemma S″ target-set density
   argument; p−1 density). Is it correct that sub-2/5 for a special family needs a rigorous ECM-speed recogniser? Is
   there a known rigorous ECM statement (Lenstra 1987, Pomerance 1987) that already supplies (O2)? Search the web.
5. Prop RB (relativisation barrier): is the statement precise and true? Is the SNFS numerology (|V| = N^{2ν+o(1)},
   5ν < 1/d) right? Is the contrast with Theorem U consistent?
6. Prop CY / Corollary CY-P: correctness, the bound Φ_a(r) ≤ (r+1)^{φ(a)}, the φ(2e) variant, the o(1) along
   primorials, and the claim that algebraic factors cannot change α.
7. §4: Hypothesis RT, the digit description of a − b r^q (proportion (1 − dν)/(1 + dν)), and the literature comparisons
   (BBDT, Barbulescu–Lachand, Bourgain/Swaenepoel).
8. Anything over-claimed, mislabelled, or missing (e.g. a device or family that does give rigour below 2/5).

Report: for each item CORRECT / CORRECT WITH FIXES / GAP / WRONG with the exact fix; severity FATAL / MAJOR / MINOR.
Verdict first. ≤ 3000 words.
