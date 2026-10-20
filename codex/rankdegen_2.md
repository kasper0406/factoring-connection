# rankdegen_2 — adversarial review of notes/rankdegen_results.md

You are a hostile but fair referee (geometry of numbers, Diophantine approximation, algebraic number theory, NFS
complexity). Read-only access to /Volumes/git/factoring-connection; scratch in your cwd; python available
(/Volumes/git/factoring-connection/.venv/bin/python has sympy, python-flint, gmpy2, numpy; run with
PYTHONPATH=/Volumes/git/factoring-connection). Web search enabled.

Read notes/rankdegen_results.md IN FULL, and as needed notes/rankdegen_prereg.md, src/rankdegen/*.py,
tests/test_rankdegen.py, experiments/raw/rankdegen_E{1,2,3}.jsonl, notes/rank3_results.md (R1, R2, R3, R5, §7),
notes/resbarrier_results.md §2.2 (Lemma ND). Your own blind round is codex/rankdegen_1.out.

Attack, in particular, and try to break things numerically where possible:
1. **Theorem FH.** (a) The Hölder step with frame weights and the marginal condition. (b) The real-frame identity
   vol{|F_S| ≤ 1} = V₀·Φ_S^{−1/(r+1)} and V₀ = C_r. (c) The complex-to-real step: R_j = Re(e^{iφ_j}L_j), the
   multilinear top coefficient 2^{−r}det, iterated Jensen, and existence of one phase vector for all r+1 minors.
   (d) The claim that uniform weights reproduce Bean–Thunder for totally real forms. Is any of this wrong or known?
2. **Corollary FN (linking).** Galois-orbit integrality of Ψ_O, the valuation ≥ 1 for minors containing the
   congruent cross pair, the use of min over all cross pairs, the symmetric-LP argument, the ramified caveat.
3. **Theorem S**, especially (S2) all-or-nothing coincidence and (S4) "instability only from a whole side". Is the
   classification A/B/C1/C2/D/decoupled exhaustive for Δ_F = 0? Is type D correctly characterised?
4. **Lemma CC** (hypotheses, the finite set C, odd |S|), **Lemma SUB**, and **Lemma PR** (Pólya, |W′(z₀)| ≥ N,
   the use of FH for V_{W̄}).
5. **The Φ program and Proposition PHI**: the reduction to Φ = inf q₁q₂R (AM–GM over the degree split), the
   crossing proof, the claim that pure s suffice, the formal/scoped constant formulas, and whether the program's
   constraint set is really implied by FL for each type in the ledger table (especially C1 both-deficient, C2 with
   ε, ε′, and the balanced-degree Φ_bal values). Is any type mis-costed?
6. **§6 large boxes**: the negative statements (6.1–6.3) and whether LP_exc (6.4) is correctly formulated and
   sufficient. Is there an exceptional family that is neither a sublattice nor a product and could matter?
7. **§7 O3**: Corollary SPL (algebra), the sharpness example, the H_split program and its closed form
   c³ = 8Φ/9, and the verdict.
8. **Experiments**: re-score P1–P9 from the raw files; check the disclosed E3 deviation and anything undisclosed.
Give a verdict (accept / minor / major revision) and a numbered list of issues with severity (MAJOR/MINOR) and the
exact fix. Under ~3500 words.
