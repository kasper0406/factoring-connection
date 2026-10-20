# almostall13_2: adversarial referee of notes/almostall13_results.md

You are an adversarial referee. Repository /Volumes/git/factoring-connection is read-only for you (scratch only under
codex/scratch/almostall13_2/). Read notes/almostall13_results.md IN FULL, and as needed notes/nested_results.md (§§3-8),
notes/lv5unif_results.md (§§6, 9, 10), notes/recurse_results.md §3, src/almostall13/*.py, experiments/raw/almostall13_e*.

Check every labelled claim, especially:
1. Lemma DF and Prop NB (i)-(iv): is "descent depends only on the field" true for NEST, including D_A, targets, factor
   exits? Is NB(iv)'s dimension count right (design a''' and plain base-m)? Any design in the project where N enters the
   descent other than through the field?
2. Remark W-AA: the congruence F^(P) = (b g1^{-1})^d N w mod t and its use in LD*-U Step 4.
3. Lemma M1-AA (a)-(c): the Rankin bound, the product over omega(R) smallest primes, the range L <= (log X)^{1.9},
   binning, the B-rough variant, and the claimed consequences for Prop R2.
4. Lemma MP and the consequence (union bound needs relative precision e^{-b1 Q}); is the "moment order within budget"
   remark correct?
5. Prop HR (a)-(c) and Theorem NR-AA: the clocked-recursion strong induction (Markov, r = ceil(E) repetitions, calls on
   m not in K), Lemma S'' on S' = U_A ∩ K, the field-failure probability, (T), the unconditional characters at
   d = delta D' (cost L[1/3, C_chi delta]), the exponent arithmetic of c_NR-AA, the FACTOR_U threshold n_U = exp(t^{4/5}).
6. Definition AA, Lemma TR, and the verdict "best alpha = 2/5 for almost all N".
7. Literature §7 (attributions: Hafner-McCurley 1989 Thms 3.1-3.2; Lenstra-Pomerance 1992 §7 Thms 7.1-7.4 6/7 exponent;
   Lee-Venkatesan; Kleinjung-Wesolowski; Barbulescu-Jouve).
8. E1/E2 numbers against the raw files.

For each issue give severity FATAL / MAJOR / MINOR, location, and a concrete fix. End with a per-claim verdict table
(CORRECT / CORRECT WITH FIXES / INCORRECT) and an overall recommendation. ≤ 4000 words.
