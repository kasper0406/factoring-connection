# abcsunit round 2 (ADVERSARIAL): attack the draft results note

You are an adversarial referee (Diophantine equations / S-unit equations / abc, smooth numbers and the circle
method, index-calculus factoring, lattices). Read `/Volumes/git/factoring-connection/notes/abcsunit_results.md`
(DRAFT) and the code it cites in `/Volumes/git/factoring-connection/src/abcsunit/` (core.py, model.py,
run_experiment.py, score.py, posthoc.py), the pre-registration `notes/abcsunit_prereg.md`, the raw results
`experiments/raw/abcsunit_{E123,E4,posthoc}.jsonl`, and tests `tests/test_abcsunit.py`. Context you may read:
notes/pilist_results.md §1.1, notes/nonexact_results.md §2.1 and §6, notes/smoothenum_results.md (T7/T8),
notes/systeq_results.md §4–§6. Your own blind round is `codex/abcsunit_1.out`. The repo is read-only for you;
scratch goes in the current directory (exact integer arithmetic for any checks).

Find errors, gaps, overclaims and mislabels. For each finding give: location, severity (FATAL / MAJOR / MINOR /
NIT), what is wrong, and a counterexample or corrected statement where possible. Specifically:

1. **§1 translation.** Is the primitive normalisation (d | k, gcd(k′, u′v′) = 1) correct in all cases (signs,
   u > v vs v > u, N with a prime ≤ y)? Is the NFS "mixed congruence" description accurate (class group, units,
   the degree-one primes above N, the claim that with b smooth it becomes an S-unit equation in K)? Is the
   classification table right?
2. **Prop ABCV and its corollary.** Check every inequality, including rad(k′N) = rad(k′)N/g, the two cases, the
   powerful-part bound k/rad(k) ≤ P(k), #powerful ≤ ζ(3/2)√t, the partial summation, and the handling of g > 1.
   Is "abc is vacuous for relations" an overclaim? Construct relation families (if any) that abc *would* constrain
   and that a factoring method might plausibly use (e.g. relations with highly powerful u, v, or k).
3. **Cor ABC-T, Prop EV, Lemma SC.** Is the reduction to Evertse's fixed-coefficient equation N⁻¹X + N⁻¹Y = 1 in
   Γ_S injective as claimed (signs, gcds)? Is Evertse's bound applied with the right place count? Is Lemma SC's
   injection valid? Is the reading "never bites in any regime" correct — e.g. does a method that needs relations
   of *bounded* height (not just ≥ s relations in total) change anything?
4. **Lemma FR** (§3.2): check items 1–5, especially dim im z = ρ − r, the odd-index statement, and the fruitfulness
   probabilities for ρ − r ∈ {0, 1, 2}. Is it correct that r ≤ s − dim(H/2H) and that dim H/2H = 2 for H = (Z/pq)^*?
   Is the claim that one-sided (QS/Dixon) relations avoid the issue correct? Prior art?
5. **Prop DW** (§3.3): is the description of de Weger's ℓ-adic step accurate? Is Bach's theorem quoted correctly
   (GRH; bound 2 log² m; prime vs integer)? Is "computing det Λ_N factors N" correct as stated? Is the conclusion
   "the step de Weger handles trivially is, for N, at least as hard as factoring" fair, and what remains open?
6. **§4.** Check the LS statements (2011 JTNB, 2012 PLMS), the Drappeau 2013 and Harper 2016 statements, Prop EX
   (exponents: the primitive count N^{1−1/κ+o(1)}; the use of Ψ(x, (log x)^κ) = x^{1−1/κ+o(1)} with y fixed in
   terms of log N rather than log X), the check of your eq. (4.3) (threshold A > κ/(κ−2)), the heuristic
   triply-smooth exponent 2K − 1 − (3K−1)/κ and threshold κ > (3K−1)/(2K−1). Is §4.3's diagnosis ("existence is
   free, constructivity is the failure, a constructive N-version is exactly Π") correct and complete? Is there any
   constructive content in LS / Drappeau / Harper (e.g. via explicit major-arc structure, or via Harper's Bourgain-
   type argument) that could be turned into a search procedure? Is there any known efficient construction of
   polylog-smooth abc triples (with or without N) that we missed?
7. **Experiments.** Check the scoring in `score.py` against the pre-registration (any mis-scoring, any prediction
   scored on the wrong subset). Check the model M1 (`model.py`: σ_ℓ, the shift series, the use of one saddle point,
   the left window near X, the u-subsample) for errors. Is the post-hoc PH1 dispersion model right (Σ m², with
   m = Ψ(X/v′))? Is the PH2 explanation of the two P5 failures (Z-rank < s ⇒ ρ < s) correct? Are the claims in §5
   "Reading" supported? Rerun anything you doubt at small size.
8. **Verdict.** Is the plain-language verdict accurate and not overclaimed? What is the single most important
   correction?

Final answer: a numbered list of findings with severities, then a short overall assessment (ACCEPT / ACCEPT WITH
FIXES / MAJOR REVISION), ≤ 4000 words.
