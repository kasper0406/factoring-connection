# pilist round 2 (ADVERSARIAL): attack the draft results note

You are an adversarial referee (analytic number theory, list decoding, lattice methods, additive combinatorics).
Read `/Volumes/git/factoring-connection/notes/pilist_results.md` (DRAFT) and the code it cites in
`/Volumes/git/factoring-connection/src/pilist/` (decoder.py, pencil.py, tails.py, targets.py, posthoc.py,
run_experiment.py), the pre-registration `notes/pilist_prereg.md`, and raw results
`experiments/raw/pilist_E{1,2,3}.jsonl`, `pilist_E2_posthoc.json`, `pilist_E2_coprime_posthoc.json`.
Context notes you may read: notes/smoothenum_results.md, notes/nonexact_results.md §6, notes/ptebarrier_results.md.
The repo is read-only for you; scratch goes in the current directory. Your own blind round is `codex/pilist_1.out`.

Your job is to find errors, gaps, overclaims and mislabels. For each finding give: location, severity
(FATAL / MAJOR / MINOR / NIT), what is wrong, a counterexample or a corrected statement where possible.
Specifically:

1. **Prop JR** (§2.2). Check Lemma S≥V and every regime of the proof: (i) Rankin with s = 2β, the
   rearrangement bound for Σ_{q|P} q^{σ−1}, the "either β²Λ dominates or β ≤ CΛ^{−η/2}" step; (ii) the
   complementary-divisor Rankin with κ = β/2 and its exponent algebra; (iii) the two sub-cases. Are the regime
   boundaries exhaustive? Is V ≤ U legitimately assumed? Is the model H_div (independence; Pr[v_q ≥ j] ≤ ρ_q q^{-j};
   cofactor gain ≤ d^λ) stated so that the bound is meaningful for actual decoders? Is the conclusion
   log Γ = o(F_{1/3}) correct with all terms? Try to construct a P (any product of factor-base prime powers) and a
   call satisfying (J) with Γ = L^{Ω(1)}.
2. **Weighted (GSS) case.** Verify Lemma CS. Try to PROVE the weighted analogue of Prop JR (Γ = L^{o(1)} for every
   weight vector), or find a weight vector that breaks it. Check the GSS decoding condition as stated.
3. **Prop EQ-SEL** (§2.1): is the equivalence correctly stated and proved? Is "LQ is the contrapositive of no
   partial lister" an overclaim?
4. **Prop TAIL** (§3.3): verify the first-moment bound and its adaptive-strategy justification; the H_tail model
   (band (B^{1−ε}, B], ε = C/log u, Θ(u) band primes, p = C'u log u/B^{1−ε}); the claimed consequence that open item
   (f) reduces to H_sel + sparse LA.
5. **O2** (§4): Lemma TL, Lemma PEN (all four parts; the line analysis on the root-pair curve; e ≥ 3 constant 16/3;
   the use of Binyamini–Cluckers–Novikov), Prop O2-count (the quotient linearity; injectivity), and the
   Behrend-type argument that three-term rigidity cannot bound O2.
6. **Cost targets** (§1.3, Table 1.1): check the congruence-instance formula, τ_max, the θ* values, the claim that
   congruence instances need θ → 1 even to beat MNFS.
7. **Experiments** (§5): re-score P1–P12 from the raw files; check the disclosed E1 deviation and the P4/P7
   post-hoc analyses for honesty; check that the 1/(q+1) primitive density is right for a homogeneous linear form.
8. **Literature**: check every citation in §7 for existence and correct attribution (in particular Chinburg–
   Hemenway–Heninger–Scherr 2016, GSS 2000 condition, Herrmann–May Thm 4 formula, Cheng–Wan, Hanrot–Stehlé,
   Behrend/Ruzsa).
9. **Verdict**: is "one wall, three faces" (§6) supported? Overclaimed?

Give an overall recommendation (ACCEPT / ACCEPT WITH FIXES / MAJOR REVISION / REJECT). ≤ 5000 words.
