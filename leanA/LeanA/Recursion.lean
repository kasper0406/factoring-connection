import LeanA.Defs

/-!
# Lemma F / Lemma F_α (solving the recursion) and the portfolio bookkeeping

[pp] Lemma lem:F (α = 1/2) and lv5unif §10 Lemma F_α (general `0 < α < 1`).

The main statements (MAIN) are copied into `Challenge.lean`; do not change them.
-/

noncomputable section
open Classical

namespace FactoringA

open Real Filter Topology

/-- MAIN. Lemma F_α (lv5unif §10; [pp] Lemma lem:F is `α = 1/2`). Let `F` be nonnegative and
nondecreasing on `[2, ∞)` and suppose, for `x ≥ x₀`,
`F(x) ≤ L_x^{p+ε(x)} + L_x^{δ+ε(x)} F(√(2x))` with `ε(x) → 0`, where `L_x^{b} = LN α b x`.
Then `F(x) ≤ L_x^{c+o(1)}` with `c = max(p, δ/(1 − 2^{−α}))`: for every `η > 0` there is `K` with
`F(x) ≤ K L_x^{c+η}` for all `x ≥ 2`. -/
theorem lemmaF (α : ℝ) (hα0 : 0 < α) (hα1 : α < 1) (F : ℝ → ℝ)
    (hF0 : ∀ x, 2 ≤ x → 0 ≤ F x) (hmono : MonotoneOn F (Set.Ici 2))
    (p δ : ℝ) (hp : 0 ≤ p) (hδ : 0 ≤ δ) (ε : ℝ → ℝ) (hε : Tendsto ε atTop (𝓝 0)) (x₀ : ℝ)
    (hrec : ∀ x, x₀ ≤ x →
      F x ≤ LN α (p + ε x) x + LN α (δ + ε x) x * F (Real.sqrt (2 * x))) :
    ∀ η > 0, ∃ K : ℝ, ∀ x, 2 ≤ x → F x ≤ K * LN α (max p (δ / (1 - 2 ^ (-α))) + η) x := by
  sorry

/-- MAIN. The form of Lemma F used for Theorems RA and U: the recursion inequality of [pp]
Theorem thm:RA, `F(x) ≤ 2 log₂x [P(x) + 2 D(x) F(√(2x))]` with `P(x) ≤ L_x^{p+o(1)}` and
`D(x) ≤ L_x^{δ+o(1)}`, gives `F(x) ≤ L_x^{max(p, δ/(1−2^{−α}))+o(1)}`. -/
theorem lemmaF_logfactor (α : ℝ) (hα0 : 0 < α) (hα1 : α < 1) (F P D : ℝ → ℝ)
    (hF0 : ∀ x, 2 ≤ x → 0 ≤ F x) (hmono : MonotoneOn F (Set.Ici 2))
    (p δ : ℝ) (hp : 0 ≤ p) (hδ : 0 ≤ δ)
    (hP : ∀ η > 0, ∀ᶠ x in atTop, P x ≤ LN α (p + η) x)
    (hD : ∀ η > 0, ∀ᶠ x in atTop, D x ≤ LN α (δ + η) x)
    (hD0 : ∀ᶠ x in atTop, 0 ≤ D x)
    (hrec : ∀ᶠ x in atTop, F x ≤ 2 * Real.logb 2 x * (P x + 2 * D x * F (Real.sqrt (2 * x)))) :
    ∀ η > 0, ∃ K : ℝ, ∀ x, 2 ≤ x → F x ≤ K * LN α (max p (δ / (1 - 2 ^ (-α))) + η) x := by
  sorry

/-- MAIN. [pp] Corollary cor:P, the series of the restart schedule:
`2^{s₀+1} + Σ_{m ≥ 1} 2^{s₀+m} 2^{−m(m+1)/2} ≤ 4·2^{s₀}`. -/
theorem restart_series (s₀ : ℕ) :
    (2 : ℝ) ^ (s₀ + 1) + ∑' m : ℕ, (2 : ℝ) ^ (s₀ + (m + 1)) * (2 : ℝ)⁻¹ ^ ((m + 1) * (m + 2) / 2) ≤
      4 * 2 ^ s₀ := by
  sorry

/-- MAIN. [pp] Corollary cor:P, assembly: if `T_j(N) ≤ L_N[1/2, c_j + o(1)]` for each `j`, with
`c_j → 0`, and for each `j` the portfolio's expected time is eventually at most
`K_j (log N)^{K_j} T_j(N)`, then the portfolio runs in `L_N[1/2, o(1)]`. -/
theorem portfolio_assembly (T : ℕ → ℝ) (Tj : ℕ → ℕ → ℝ) (c : ℕ → ℝ)
    (hc : Tendsto c atTop (𝓝 0)) (hTj : ∀ j, LBound (1 / 2) (c j) (Tj j))
    (hport : ∀ j, ∃ K : ℝ, 0 ≤ K ∧ ∀ᶠ N : ℕ in atTop, T N ≤ K * (Real.log N) ^ K * Tj j N) :
    LBound (1 / 2) 0 T := by
  sorry

end FactoringA
