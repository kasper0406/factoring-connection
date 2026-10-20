# systeq round 2 (ADVERSARIAL): break our barriers for integer analogues of the systematic equation

You are an adversarial referee. Read `notes/systeq_results.md` in this repository (read-only),
and optionally `src/systeq/*.py`, `notes/systeq_prereg.md` and `experiments/raw/systeq_*.jsonl`.
Your blind round is `codex/systeq_1.out`. Scratch work goes in the current directory.

Your job is to find errors and escapes, not to agree. Specifically:

1. **Prop F (§5.1).** Is the proof correct? Is the hypothesis "z ↦ z^s additive on all of A" too strong
   for what BGJT actually uses? BGJT uses additivity only on elements of the form aX + b with a, b
   in the base field k′ and X the generator. Could an exponent s exist with (aX+b)^s = a′X^s + b′ for
   all a, b in a small subset, without global additivity, in some Z/N-algebra, computable from N?
   Construct one or prove that the partial version also forces (p−1)|(s−1) or similar.
2. **Lemma CF (§5.2).** Is the claim that "the only computable systematic equation over Z/N has
   constant field F_p" overstated? Consider S = μ_n ∪ {0} inside (Z/N)[ζ_n] with the product
   T^{n+1} − T, and an abstract automorphism σ instead of a power map. Is there any way to make
   ∏_{c∈S}(z − c) equal to a low-complexity expression in σ(z), z for the elements z one needs?
3. **Lemma DES (§5.3).** Check the case analysis for V_n and the claim "no nondegenerate split
   member for n ≥ 6". Check the Lenstra 1999 citation and statement. Is there a characteristic-0
   family other than fewnomials (e.g. polynomials of the form F(x^q) − G(x), compositions,
   Chebyshev-type, or polynomials over Z[ζ_n] where splitting is free) that keeps E4's "small side"
   while making the big side split into many small pieces? Be concrete; compute a small example.
4. **Rigidity (§5.5) and E5 (§5.6).** Is there an integer setting with presentation freedom we
   missed, e.g. replacing Z/N by a finite Z/N-algebra A = (Z/N)[x]/(F) (where factoring A still
   factors N) and choosing F freely so that "x^n ≡ h(x) mod F" holds with small h? Here the
   presentation IS free (F is ours). Does that give BGJT-like relations whose sides are smooth
   integers, or do the norms N_{A/(Z/N)} reintroduce N-size? Work it out quantitatively: what is
   the size of the integers one must factor? This is our biggest worry, so please be thorough.
5. **Experiments.** Are E1 (positive-control gain G_FF), E2 (Cunningham V-family G ≈ 1.04) and
   E4 fair tests of what they claim? Is the size-bound prediction set up correctly? Any bias?
6. **Literature.** Are any of our "no prior statement found by us" items (Prop F as stated,
   Lemma CF, Lemma DES applied to the BGJT space V, the F_1 PGL₂ remark) known? Cite.

For each item give a verdict: CORRECT / CORRECT WITH FIXES / WRONG, with the fix. Label every
claim KNOWN / PROVED / HEURISTIC / SPECULATIVE. Keep it under ~3000 words.
