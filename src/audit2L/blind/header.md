# audit2L blind claims: Theorem 2L-final (rigorous Las Vegas factoring in L_N[α], 3/8 ≤ α < 2/5)

This file is a blind excerpt prepared for independent review. It contains the statements and proofs of a chain of
internal (unrefereed) research notes, with all review sections, reviewer verdicts, probability estimates and reviewer
attributions removed. Some sentences still describe how an earlier draft was repaired ("v1 used …"); these are
history, not verdicts. Parenthetical attributions were stripped mechanically, so a few sentences may read abruptly.

**Reviewer instruction.** You may read the underlying source notes (paths below) for statements this excerpt cites
but does not reproduce, but do **not** read the review sections of those notes (sections titled "Reviews", the
"Verdict"/"0. Verdict" sections, and any probability estimates) before you have formed and written down your own
verdict.

## The claim under review

Notation: λ := log log N, D := log N/λ, L[α, c] = L_N[α, c] = e^{(c+o(1))D^αλ}; "polylog" = e^{O(λ)}.

- **Theorem 2L-final(α).** For every fixed α ∈ (3/8, 2/5), a Las Vegas algorithm (FACTOR_2L in design 2L_DW(α))
  outputs the complete prime factorisation of every integer N ≥ 2 in expected time L_N[α, c(α)], where
  c(α) = max(2a, a + s + max(a, σ₁ + σ₂), (a + s)/(1 − 2^{−α})), s = 2κ(1 + C_γ), σ₁ = α/(ad₀),
  σ₂ = (α/a)(4/d₀ + d₀b₂), for any admissible rational design constants a, κ, d₀, b₂ > 0, and C_γ is the polylog
  exponent of the level-2 relation supply γ₂ ≥ (log N)^{−C_γ}.
- **Theorem 2L-final(3/8).** The same at α = 3/8 under the design constraint a > 2c₂, c₂ := 3d₀/(32κb₂), with
  σ₁(3/8) = 3/(8ad₀) + c₂.

Status claimed by the authors: proved relative to the project's internal unrefereed chain (paper2/main.tex,
notes/lv5unif_results.md = Theorem U at L[2/5], and the notes excerpted below) plus published results, of which the
essential one is Drappeau, Compositio Math. 151 (2015), Lemme 5 (stated there without proof as a consequence of
§3.3 of Harper's unrefereed preprint arXiv:1208.5992).

## Logical structure (as claimed)

1. Theorem U's one-level chain (paper2/main.tex, notes/lv5unif_results.md) supplies the split/field-selection
   machinery, LV6 (square roots), fruitfulness, recursion (Lemma F_α) and the second-moment tools.
2. **Hyp U′** (notes/uprime_results.md, excerpt E; with the DW replacements of c12alg §2.7, excerpt F): level-1
   moments μ₁ ≥ (log N)^{−C}, E[N₁²] ≤ (log N)^C(μ₁² + μ₁).
3. **Prop UB** (notes/multidescent_results.md, excerpt D): if the dead parent–child incidence mass is ≤ μ₁/2 then
   Pr[N₂ ≥ 1] ≥ μ₁²/(4E N₁²).
4. **Prop DC / DC_e** (notes/pc3_results.md, excerpt C): dead mass ≤ W(A + √(AB)), A, B the incidence-weighted
   relative variances of the child count and of the label's parent count.
5. **Lemma η, Prop PB, Lemma VD** (notes/c12alg_results.md, excerpt F): W/μ₁ ≤ 3u₂; dead mass ≤ μ₁/2 if
   A ≤ D^{−α}λ^C and B ≤ D^{7α−3}λ^C; A, B split into one-point, rational and algebraic blocks.
6. **Theorem RAT** (notes/c12rat_results.md, excerpt B): the rational blocks at second order.
7. **Theorem AB/AB₂** (notes/thmAB_results.md, excerpt A): the algebraic blocks at second order, and the final
   assembly (§13 of that note, Theorem 2L-final).
8. Prop CLB (excerpt G): one-point child lower bound.

Other source notes that may be consulted: notes/lv5unif_results.md, paper2/main.tex, notes/nested_results.md,
notes/escapeR2_results.md, notes/lem85_results.md.
