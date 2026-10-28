## Own verdict (audit2L, recorded 2026-10-28 ~10:55 CET, BEFORE reading any reviewer output)

Recorded before the codex referee finished (only its process/log size was checked) and before either Opus critic was
launched.

- **(i) α ∈ (3/8, 2/5): ACCEPT WITH FIXES**, relative to the internal one-level chain (paper2 + lv5unif = Theorem U)
  and Drappeau's Lemme 5 (whose proof is in Harper's unrefereed preprint).
- **(ii) α = 3/8: ACCEPT WITH FIXES, at lower confidence** than (i). Everything must be genuinely second order there,
  the acceptability subtraction is closed only by the constant constraint a > 2c₂, and U′ under DW at η₀ = D^{−3/8}
  has had only a leading-order re-check.

Findings of my own targeted review (detail in §1 of the audit note):
1. thmAB round-3 repairs: all four are correct as applied.
   - The marginal lower bound (OPV-A step 6) holds via the domination of ν̂ on interior cells by ν, with no
     circularity.
   - The parent RAT bridge (value-only a_P^{id}; ε_p ≤ D^C/(Bρ_a(q)); perturbation of a relative variance) is correct.
   - The BTS placement (covering interval inside the piece's own dyadic range) satisfies gcd·W ≤ X₀ and log W ≥ ½ log X.
   - The β-drift (D^{ζ−2α}/λ for children, D^{4α−2}/λ for parents) is correct, and its squares fit the budget.
2. c12rat v2 repairs: TYPE, CP, CPW, FIB, DOM, FAM, WA, thin pairs (TF), skewed labels (LS + DOM), PRIM, OPV (Z ≥ 1),
   FO small/mid Rankin, and the Prop DC/RAT bridge (zero component, 1/t corrections, idealised weights) all check.
   - **One accounting flaw, MINOR:** RAT-2's hard-filter step bounds the parent pair sum by "polylog·(δ_e + 𝒱₀)", with
     polylog = e^{O(λ)} = D^{O(1)}. Read literally, that does not give the stated λ^C rate. The spread of the OPV
     constants c_P among parents that pass the filter is in fact 1 + o(1) (factor-64 value range, b-window,
     (64ed)^{γ₂}), and points that never pass contribute 0. So the correct bound is O(𝒱₀ + δ_e²). The fix changes no
     exponent.
3. Exponents (src/audit2L/exponents.py, exact, from scratch): all 46 comparisons hold on the open window.
   - At 3/8 exactly three fail. Two are the first-order budget checks (A/η and √(AB)/η), which second order replaces;
     the second-order checks pass. The third is the acceptability tie 3 − 7α = α, which a > 2c₂ closes.
   - At 2/5 three ties appear (d log V₁ vs 2α, B₂ vs V₁, θ_p vs 3 − 7α), so the window is open at the top.
   - Gaps: children ½ − α, parents 2 − 5α. c(7/18) = 5.387 and c(3/8) = 4.706 reproduce the note's numbers.
4. Drappeau's Lemme 5: I checked it verbatim in the arXiv text (1307.7554v3). It has the 1/φ(q) weight, sums over all
   χ mod q with 1 < cond ≤ x^η, takes Q ≤ x and (log x)^c ≤ y, and its implied constant is effective for A < 1. Its use
   in EN/EN-A (Q ≤ (ξ/δ)^{η_D}) is within these hypotheses.
5. Remaining risks:
   - (a) U′ under DW in the window and at 3/8 is a stated, "routine" transfer of uprime §§3–6 that nobody has re-proved
     line by line. I checked the margins it needs (DS-U′, U1₂, L2 classes, Kern-U, CS′, IR/content rejections); they
     hold at leading order.
   - (b) The internal one-level chain, and Lemme 5 through Harper's preprint.
   - (c) The sheer size of the bookkeeping: about 4000 lines of unrefereed notes.

Weakest link (mine): the glue that transfers U′ to design DW at window and endpoint parameters (c12alg §2.7). It is
unwritten, not known to be wrong. Inherited from below it: Lemme 5 / Harper, and [pp]'s second moment.

Probabilities (mine), correct relative to the internal chain plus published inputs: (i) 70%, (ii) 55%.
