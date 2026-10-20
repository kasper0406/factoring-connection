# rank3_2 — adversarial review of notes/rank3_results.md

You are a hostile but fair referee (geometry of numbers, algebraic number theory, NFS complexity). Read-only
access to /Volumes/git/factoring-connection; scratch in your cwd; python available (.venv/bin/python has sympy,
python-flint, gmpy2). Web search enabled.

Read notes/rank3_results.md in full (and, as needed, notes/rank3_prereg.md, src/rank3/*.py, tests/test_rank3.py,
experiments/raw/rank3_E{1,2}.jsonl, notes/resbarrier_results.md §2.3, §2.6, §4.2, notes/ptesplit_results.md §2.4).
Your own blind round is codex/rank3_1.out.

Attack, in particular:
1. **Theorem R1 (ii)** — the inclusion–exclusion proof that ramified linking primes still give
   v_p(Δ_F) ≥ r(r−1)(D−2)!/(D−r)!. Check every step: the set T_K and its size e·f, the valuation 1/e(𝔭_K|p) of
   σ(u_j) − ṽ_j for σ ∈ T_K, the difference-of-rows step, S(t) ≥ (t−1)C(D−2,r−2), and the final inequality
   (e+e′−1)/max(e,e′) ≥ 1. Is anything needed about p | [O_K : O], residue degree > 1, or p ≤ D? Try to break it
   numerically (e.g. Eisenstein / index-divisor examples).
2. **Δ_F ∈ Z** and its well-definedness for our factorisation; Bean–Thunder's exact hypotheses (real vs complex
   coefficients, D > r, Δ ≠ 0, any normalisation) and the exact constant C_r. Verify against the original if you
   can reach it.
3. **Lemma R2 and Corollary R3**: the yield step, H5_r, H7_r, the α ≥ 1/3 proof for growing r, the constants table
   (scoped and unscoped), and whether the claim "rank 2 optimal" is stated with the right scope.
4. **Proposition R4** (cross-polytope argument) and the claim that the resbarrier §4.2 aggregate family is
   impossible "for every N as specified". Is the interpretation of §4.2's premise (family height N^{1/d} meaning
   generator max-heights) fair? Is there a reading of §4.2 that escapes R4?
5. **Split side (§6)**: is R5 really the "best possible two-norm ledger"? Is the O3 statement right?
6. **§7 open list**: anything missing or mis-stated (e.g. is GP really automatic for power bases; can a design be
   in non-general position yet useful; is LP correctly described as weaker than G_r)?
7. **Experiments**: re-score P1–P7 from the raw files; check the disclosed deviations.
Give a verdict (accept / minor / major revision) and a numbered list of issues with severity (MAJOR/MINOR) and the
exact fix. Keep it under ~3000 words.
