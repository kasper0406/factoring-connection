import LeanA.Defs

/-!
# (a) Exponent and constant arithmetic

All statements are about the explicit exponent functions of `LeanA/Defs.lean` and carry no
hypotheses beyond their arithmetic side conditions. The main statements (marked MAIN) are copied
into `Challenge.lean`; their statements must not be changed.
-/

noncomputable section
open Classical

namespace FactoringA

open Real

/-- MAIN. [pp] Lemma lem:optmu, lower bound: for `d ≥ 3`, `k ≥ 1` and every `μ ≥ ν`, θ is at least
`θ_min(d, k) + η₀`. -/
theorem theta_ge_thetaMin (P : Params) (hd : 3 ≤ P.d) (hk : 1 ≤ P.k) (hμ : P.ν ≤ P.μ) :
    thetaMin P.d P.k + P.η₀ ≤ P.θ := by
  sorry

/-- MAIN. [pp] Lemma lem:optmu, attainment: at `μ = μ_opt(d, k)`, `θ = θ_min(d, k) + η₀`. -/
theorem theta_at_muOpt (d k : ℕ) (η₀ : ℚ) (hd : 3 ≤ d) (hk : 1 ≤ k) :
    (Params.θ ⟨d, k, muOpt d k, η₀⟩) = thetaMin d k + η₀ := by
  sorry

/-- MAIN. `(1+θ)/√2 < 1 ↔ θ < √2 − 1` ([pp] §1.2). -/
theorem cA_lt_one_iff (P : Params) (hθ : 0 ≤ (P.θ : ℝ)) :
    P.cA < 1 ↔ (P.θ : ℝ) < Real.sqrt 2 - 1 := by
  sorry

/-- MAIN. The constant of every admissible tuple of degree `d` exceeds the infimum
`(1 + 3/(d+1))/√2` (so the infimum is not attained; [pp] Remark rem:claims(a)). -/
theorem cA_gt_inf (P : Params) (hP : P.Admissible) :
    (1 + 3 / ((P.d : ℝ) + 1)) / Real.sqrt 2 < P.cA := by
  sorry

/-- MAIN. At fixed odd `d ≥ 3`, the constants of admissible tuples come arbitrarily close to
`(1 + 3/(d+1))/√2` ([pp] §1.2, display after Theorem thm:main). -/
theorem exists_admissible_fixed_d (d : ℕ) (hd : 3 ≤ d) (hodd : Odd d) (c : ℝ)
    (hc : (1 + 3 / ((d : ℝ) + 1)) / Real.sqrt 2 < c) :
    ∃ P : Params, P.d = d ∧ P.Admissible ∧ P.cA < c := by
  sorry

/-- MAIN. [pp] Corollary cor:main, exponent part: for every `c > 1/√2` some admissible tuple has
`(1+θ)/√2 < c`. -/
theorem exists_admissible_cA_lt (c : ℝ) (hc : 1 / Real.sqrt 2 < c) :
    ∃ P : Params, P.Admissible ∧ P.cA < c := by
  sorry

/-- MAIN. [pp] Table tab:constants, d = 9, k = 13: μ = 17/156, θ = 16/39 + η₀, admissible for
`0 < η₀ ≤ 1/100`, and the constant is below 1 exactly when `η₀ < √2 − 1 − 16/39`. -/
theorem d9_k13 (η₀ : ℚ) (h0 : 0 < η₀) (h1 : η₀ ≤ 1 / 100) :
    (Params.θ ⟨9, 13, 17 / 156, η₀⟩) = 16 / 39 + η₀ ∧ Params.Admissible ⟨9, 13, 17 / 156, η₀⟩ ∧
      (Params.cA ⟨9, 13, 17 / 156, η₀⟩ < 1 ↔ (η₀ : ℝ) < Real.sqrt 2 - 1 - 16 / 39) := by
  sorry

/-- MAIN. [pp] Table tab:constants: at `d = 9` no admissible tuple with `k ≤ 12` has constant
below 1. -/
theorem d9_k_le_12 (P : Params) (hP : P.Admissible) (hd : P.d = 9) (hk : P.k ≤ 12) :
    1 ≤ P.cA := by
  sorry

/-- MAIN. [pp] Prop. prop:cost and the remark after it: for `0 < θ < 1` the cost exponent
`max(a + 1/√2, a + θ/(2a), 2a)` is minimised at `a = θ/√2`, with minimum `(1+θ)/√2`. -/
theorem costExpA_min (θ : ℝ) (h0 : 0 < θ) (h1 : θ < 1) :
    (∀ a > 0, (1 + θ) / Real.sqrt 2 ≤ costExpA θ a) ∧
      costExpA θ (θ / Real.sqrt 2) = (1 + θ) / Real.sqrt 2 := by
  sorry

/-- MAIN. [pp] Theorem thm:RA, optimisation: for `θ > 0`, `c(a) ≥ √((1+√2)θ)` for all `a > 0`,
with equality at `a* = √((√2−1)θ/2)`. -/
theorem costExpRA_min (θ : ℝ) (h0 : 0 < θ) :
    (∀ a > 0, Real.sqrt ((1 + Real.sqrt 2) * θ) ≤ costExpRA θ a) ∧
      costExpRA θ (Real.sqrt ((Real.sqrt 2 - 1) * θ / 2)) = Real.sqrt ((1 + Real.sqrt 2) * θ) := by
  sorry

/-- MAIN. [pp] Remark rem:RA(a): for `0 ≤ θ`, `√((1+√2)θ) < 1 ↔ (1+θ)/√2 < 1`, and for
`0 < θ < √2 − 1`, `√((1+√2)θ) < (1+θ)/√2`. -/
theorem cRA_vs_cA (θ : ℝ) (h0 : 0 ≤ θ) :
    (Real.sqrt ((1 + Real.sqrt 2) * θ) < 1 ↔ (1 + θ) / Real.sqrt 2 < 1) ∧
      (0 < θ → θ < Real.sqrt 2 - 1 →
        Real.sqrt ((1 + Real.sqrt 2) * θ) < (1 + θ) / Real.sqrt 2) := by
  sorry

/-- MAIN. [pp] Corollary cor:P: the portfolio parameter tuples are admissible for every `j ≥ 1`,
satisfy `θ_j < 3/(d_j+1) + 5/d_j²`, and `θ_j → 0`. -/
theorem portfolioParams_spec :
    (∀ j : ℕ, 1 ≤ j → (portfolioParams j).Admissible ∧
      ((portfolioParams j).θ : ℝ) < 3 / ((2 * j + 7 : ℝ) + 1) + 5 / (2 * j + 7 : ℝ) ^ 2) ∧
    Filter.Tendsto (fun j => ((portfolioParams j).θ : ℝ)) Filter.atTop (nhds 0) := by
  sorry

/-- MAIN. For every `c > 0` some admissible tuple has `√((1+√2)θ) < c` (the exponent part of
[pp] Corollary cor:P). -/
theorem exists_admissible_cRA_lt (c : ℝ) (hc : 0 < c) :
    ∃ P : Params, P.Admissible ∧ P.cRAstar < c := by
  sorry

/-- MAIN. [pp] Theorem thm:DL: `c_DL(θ) = min_a max(a + 1/√2, a + θ/(2a), 3a)` equals
`(1+θ)/√2` for `0 < θ ≤ 1/2` (at `a = θ/√2`) and `3√θ/2` for `1/2 < θ < 1` (at `a = √θ/2`). -/
theorem costExpDL_min (θ : ℝ) (h0 : 0 < θ) (h1 : θ < 1) :
    (θ ≤ 1 / 2 → (∀ a > 0, (1 + θ) / Real.sqrt 2 ≤ costExpDL θ a) ∧
        costExpDL θ (θ / Real.sqrt 2) = (1 + θ) / Real.sqrt 2) ∧
      (1 / 2 < θ → (∀ a > 0, 3 * Real.sqrt θ / 2 ≤ costExpDL θ a) ∧
        costExpDL θ (Real.sqrt θ / 2) = 3 * Real.sqrt θ / 2) := by
  sorry

/-- MAIN. The size-model floor of [rec] §3.3 / lv5unif §13: with precision exponent `τ ≥ 0`,
the constraints `k ≳ D^{2(1−2α)}` and `k·d^τ = O(D^α)` are compatible iff
`(2+τ)(1−2α) ≤ α`, i.e. `α ≥ (2+τ)/(5+2τ)`; at `τ = 0` this is `α ≥ 2/5`. -/
theorem size_model_floor (τ α : ℝ) (hτ : 0 ≤ τ) :
    ((2 + τ) * (1 - 2 * α) ≤ α ↔ (2 + τ) / (5 + 2 * τ) ≤ α) ∧
      ((2 : ℝ) * (1 - 2 * α) ≤ α ↔ 2 / 5 ≤ α) := by
  sorry

/-- MAIN. lv5unif §10: `4.12981 ≤ 1/(1 − 2^{−2/5}) ≤ 4.12983`. -/
theorem recursion_factor_two_fifths :
    (4.12981 : ℝ) ≤ 1 / (1 - (2 : ℝ) ^ (-(2 / 5 : ℝ))) ∧
      1 / (1 - (2 : ℝ) ^ (-(2 / 5 : ℝ))) ≤ 4.12983 := by
  sorry

end FactoringA
