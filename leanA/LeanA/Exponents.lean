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

/-! ### Helper lemmas on `θ` -/

lemma hStar_ge (P : Params) : 2 * P.ν ≤ P.hStar ∧ P.μ - ((P.d : ℚ) - 2) * P.ν ≤ P.hStar ∧
    1 - (P.d : ℚ) * P.μ + P.ν ≤ P.hStar := by
  unfold Params.hStar
  exact ⟨le_max_left _ _, le_trans (le_max_left _ _) (le_max_right _ _),
    le_trans (le_max_right _ _) (le_max_right _ _)⟩

lemma theta_eq (P : Params) :
    P.θ = 2 * P.μ + P.hStar + ((P.d : ℚ) - 1) * P.ν + P.η₀ := by
  unfold Params.θ Params.h; ring

/-- The `B`–`C` lower bound (valid for every `μ`). -/
lemma theta_ge_BC (P : Params) (hd : 3 ≤ P.d) :
    3 / ((P.d : ℚ) + 1) + (4 * (P.d : ℚ) - 2) / ((P.d : ℚ) + 1) * P.ν + P.η₀ ≤ P.θ := by
  obtain ⟨_, h2, h3⟩ := hStar_ge P
  rw [theta_eq]
  have hd' : (3 : ℚ) ≤ P.d := by exact_mod_cast hd
  have hpos : (0 : ℚ) < (P.d : ℚ) + 1 := by linarith
  have e : 3 / ((P.d : ℚ) + 1) + (4 * (P.d : ℚ) - 2) / ((P.d : ℚ) + 1) * P.ν
      = (3 + (4 * (P.d : ℚ) - 2) * P.ν) / ((P.d : ℚ) + 1) := by
    field_simp
  rw [e]
  have key : 3 + (4 * (P.d : ℚ) - 2) * P.ν ≤
      (2 * P.μ + P.hStar + ((P.d : ℚ) - 1) * P.ν) * ((P.d : ℚ) + 1) := by
    nlinarith [mul_le_mul_of_nonneg_left h2 (by linarith : (0 : ℚ) ≤ (P.d : ℚ) - 2)]
  have := (div_le_iff₀ hpos).2 key
  linarith

/-- The `A`–`C` lower bound (valid for every `μ`). -/
lemma theta_ge_AC (P : Params) (hd : 3 ≤ P.d) :
    2 / (P.d : ℚ) + ((P.d : ℚ) ^ 2 + P.d - 2) / P.d * P.ν + P.η₀ ≤ P.θ := by
  obtain ⟨h1, _, h3⟩ := hStar_ge P
  rw [theta_eq]
  have hd' : (3 : ℚ) ≤ P.d := by exact_mod_cast hd
  have hpos : (0 : ℚ) < (P.d : ℚ) := by linarith
  have e : 2 / (P.d : ℚ) + ((P.d : ℚ) ^ 2 + P.d - 2) / P.d * P.ν
      = (2 + ((P.d : ℚ) ^ 2 + P.d - 2) * P.ν) / (P.d : ℚ) := by
    field_simp
  rw [e]
  have key : 2 + ((P.d : ℚ) ^ 2 + P.d - 2) * P.ν ≤
      (2 * P.μ + P.hStar + ((P.d : ℚ) - 1) * P.ν) * (P.d : ℚ) := by
    nlinarith [mul_le_mul_of_nonneg_left h1 (by linarith : (0 : ℚ) ≤ (P.d : ℚ) - 2)]
  have := (div_le_iff₀ hpos).2 key
  linarith

lemma theta_case1 (P : Params) (hd : 3 ≤ P.d) (hν : 0 < P.ν)
    (h1 : 1 < ((P.d : ℚ) ^ 2 + 1) * P.ν) (hμ : P.μ = (1 - P.ν) / P.d) :
    P.θ = 2 / (P.d : ℚ) + ((P.d : ℚ) ^ 2 + P.d - 2) / P.d * P.ν + P.η₀ := by
  rw [theta_eq]
  have hd' : (3 : ℚ) ≤ P.d := by exact_mod_cast hd
  have hpos : (0 : ℚ) < (P.d : ℚ) := by linarith
  have hdμ : (P.d : ℚ) * P.μ = 1 - P.ν := by rw [hμ]; field_simp
  have hmax : P.hStar = 2 * P.ν := by
    unfold Params.hStar
    apply max_eq_left
    apply max_le
    · have : (P.d : ℚ) * (P.μ - (P.d : ℚ) * P.ν) < 0 := by nlinarith
      have : P.μ - (P.d : ℚ) * P.ν < 0 := by
        by_contra hc; push_neg at hc; nlinarith
      linarith
    · linarith
  rw [hmax, hμ]
  field_simp
  ring

lemma theta_case2 (P : Params) (hd : 3 ≤ P.d) (hν : 0 < P.ν)
    (h1 : ((P.d : ℚ) ^ 2 + 1) * P.ν ≤ 1) (hμ : P.μ = (1 + ((P.d : ℚ) - 1) * P.ν) / ((P.d : ℚ) + 1)) :
    P.θ = 3 / ((P.d : ℚ) + 1) + (4 * (P.d : ℚ) - 2) / ((P.d : ℚ) + 1) * P.ν + P.η₀ := by
  rw [theta_eq]
  have hd' : (3 : ℚ) ≤ P.d := by exact_mod_cast hd
  have hpos : (0 : ℚ) < (P.d : ℚ) + 1 := by linarith
  have hdμ : ((P.d : ℚ) + 1) * P.μ = 1 + ((P.d : ℚ) - 1) * P.ν := by rw [hμ]; field_simp
  have hmax : P.hStar = P.μ - ((P.d : ℚ) - 2) * P.ν := by
    unfold Params.hStar
    have e : 1 - (P.d : ℚ) * P.μ + P.ν = P.μ - ((P.d : ℚ) - 2) * P.ν := by linarith
    rw [e, max_self]
    apply max_eq_right
    have : ((P.d : ℚ) + 1) * ((P.d : ℚ) * P.ν - P.μ) ≤ 0 := by nlinarith
    have : (P.d : ℚ) * P.ν - P.μ ≤ 0 := by
      by_contra hc; push_neg at hc; nlinarith
    linarith
  rw [hmax, hμ]
  field_simp
  ring

/-- MAIN. [pp] Lemma lem:optmu, lower bound: for `d ≥ 3`, `k ≥ 1` and every `μ ≥ ν`, θ is at least
`θ_min(d, k) + η₀`. -/
theorem theta_ge_thetaMin (P : Params) (hd : 3 ≤ P.d) (hk : 1 ≤ P.k) (hμ : P.ν ≤ P.μ) :
    thetaMin P.d P.k + P.η₀ ≤ P.θ := by
  unfold thetaMin
  split_ifs
  · exact theta_ge_AC P hd
  · exact theta_ge_BC P hd

/-- MAIN. [pp] Lemma lem:optmu, attainment: at `μ = μ_opt(d, k)`, `θ = θ_min(d, k) + η₀`. -/
theorem theta_at_muOpt (d k : ℕ) (η₀ : ℚ) (hd : 3 ≤ d) (hk : 1 ≤ k) :
    (Params.θ ⟨d, k, muOpt d k, η₀⟩) = thetaMin d k + η₀ := by
  have hk4 : (0 : ℚ) < 4 * (k : ℚ) := by
    have : (1 : ℚ) ≤ k := by exact_mod_cast hk
    linarith
  have hν : (0 : ℚ) < 1 / (4 * (k : ℚ)) := by positivity
  unfold thetaMin muOpt
  split_ifs with hc
  · refine theta_case1 ⟨d, k, (1 - 1 / (4 * (k : ℚ))) / d, η₀⟩ hd hν ?_ rfl
    show 1 < ((d : ℚ) ^ 2 + 1) * (1 / (4 * (k : ℚ)))
    rw [mul_one_div, one_lt_div hk4]; exact hc
  · refine theta_case2 ⟨d, k, (1 + ((d : ℚ) - 1) * (1 / (4 * (k : ℚ)))) / ((d : ℚ) + 1), η₀⟩
      hd hν ?_ rfl
    show ((d : ℚ) ^ 2 + 1) * (1 / (4 * (k : ℚ))) ≤ 1
    rw [mul_one_div, div_le_one hk4]; push_neg at hc; exact hc

/-- MAIN. [pp] Lemma lem:optmu, range: for `k > (d+1)/4` the optimal `μ` lies in `(ν, 1)`. -/
theorem muOpt_mem (d k : ℕ) (hd : 3 ≤ d) (hk : ((d : ℚ) + 1) / 4 < k) :
    1 / (4 * (k : ℚ)) < muOpt d k ∧ muOpt d k < 1 := by
  have hd' : (3 : ℚ) ≤ d := by exact_mod_cast hd
  have hk1 : (1 : ℚ) < k := by linarith
  have hk4 : (0 : ℚ) < 4 * (k : ℚ) := by linarith
  set ν : ℚ := 1 / (4 * (k : ℚ)) with hν
  have hν0 : 0 < ν := by positivity
  have hνk : ν * (4 * k) = 1 := by rw [hν]; field_simp
  have hk' : (d : ℚ) + 1 < 4 * k := by linarith
  have hdν : ((d : ℚ) + 1) * ν < 1 := by
    have := mul_lt_mul_of_pos_right hk' hν0
    linarith
  have hν4 : ν < 1 / 4 := by nlinarith
  unfold muOpt
  rw [← hν]
  split_ifs
  · constructor
    · rw [lt_div_iff₀ (by linarith)]; nlinarith
    · rw [div_lt_one (by linarith)]; linarith
  · constructor
    · rw [lt_div_iff₀ (by linarith)]; nlinarith
    · rw [div_lt_one (by linarith)]
      nlinarith [mul_lt_mul_of_pos_left hν4 (by linarith : (0 : ℚ) < (d : ℚ) - 1)]

/-- MAIN. `(1+θ)/√2 < 1 ↔ θ < √2 − 1` ([pp] §1.2). -/
theorem cA_lt_one_iff (P : Params) (hθ : 0 ≤ (P.θ : ℝ)) :
    P.cA < 1 ↔ (P.θ : ℝ) < Real.sqrt 2 - 1 := by
  unfold Params.cA
  rw [div_lt_one (Real.sqrt_pos.2 (by norm_num : (0 : ℝ) < 2))]
  constructor <;> intro h <;> linarith

lemma admissible_facts (P : Params) (hP : P.Admissible) :
    (1 : ℚ) < P.k ∧ 0 < P.ν ∧ 3 / ((P.d : ℚ) + 1) < P.θ := by
  obtain ⟨_, hd, hk, hη, _, _, _⟩ := hP
  have hd' : (3 : ℚ) ≤ P.d := by exact_mod_cast hd
  have hk1 : (1 : ℚ) < P.k := by linarith
  have hν : 0 < P.ν := by unfold Params.ν; exact div_pos one_pos (by linarith)
  refine ⟨hk1, hν, ?_⟩
  have h := theta_ge_BC P hd
  have : 0 < (4 * (P.d : ℚ) - 2) / ((P.d : ℚ) + 1) * P.ν :=
    mul_pos (div_pos (by linarith) (by linarith)) hν
  linarith

/-- MAIN. The constant of every admissible tuple of degree `d` exceeds the infimum
`(1 + 3/(d+1))/√2` (so the infimum is not attained; [pp] Remark rem:claims(a)). -/
theorem cA_gt_inf (P : Params) (hP : P.Admissible) :
    (1 + 3 / ((P.d : ℝ) + 1)) / Real.sqrt 2 < P.cA := by
  obtain ⟨_, _, hθ⟩ := admissible_facts P hP
  have hθr : (3 : ℝ) / ((P.d : ℝ) + 1) < (P.θ : ℝ) := by
    have := (Rat.cast_lt (K := ℝ)).2 hθ
    push_cast at this; exact this
  unfold Params.cA
  rw [div_lt_div_iff_of_pos_right (Real.sqrt_pos.2 (by norm_num))]
  linarith

/-- MAIN. At fixed odd `d ≥ 3`, the constants of admissible tuples come arbitrarily close to
`(1 + 3/(d+1))/√2` ([pp] §1.2, display after Theorem thm:main). -/
theorem exists_admissible_fixed_d (d : ℕ) (hd : 3 ≤ d) (hodd : Odd d) (c : ℝ)
    (hc : (1 + 3 / ((d : ℝ) + 1)) / Real.sqrt 2 < c) :
    ∃ P : Params, P.d = d ∧ P.Admissible ∧ P.cA < c := by
  have hs0 : 0 < Real.sqrt 2 := Real.sqrt_pos.2 (by norm_num)
  rw [div_lt_iff₀ hs0] at hc
  set ε := c * Real.sqrt 2 - 1 - 3 / ((d : ℝ) + 1) with hε
  have hε0 : 0 < ε := by linarith
  obtain ⟨N, hN⟩ := exists_nat_gt (5 / ε)
  obtain ⟨k, hk⟩ : ∃ k : ℕ, k = N + d ^ 2 := ⟨_, rfl⟩
  have hd' : (3 : ℚ) ≤ d := by exact_mod_cast hd
  have hkq : ((d : ℚ)) ^ 2 ≤ k := by rw [hk]; push_cast; linarith [(N.cast_nonneg : (0 : ℚ) ≤ N)]
  have hk9 : (9 : ℚ) ≤ k := by nlinarith
  have hk4 : (0 : ℚ) < 4 * k := by linarith
  set ν : ℚ := 1 / (4 * (k : ℚ)) with hν
  have hν0 : 0 < ν := by positivity
  have hνk : ν * (4 * k) = 1 := by rw [hν]; field_simp
  have hν36 : ν ≤ 1 / 36 := by nlinarith
  have hcase : ¬ (4 * (k : ℚ) < (d : ℚ) ^ 2 + 1) := by push_neg; nlinarith
  have hμ : muOpt d k = (1 + ((d : ℚ) - 1) * ν) / ((d : ℚ) + 1) := by
    unfold muOpt; rw [if_neg hcase]
  have hk1 : 1 ≤ k := by
    have : (1 : ℚ) ≤ k := by linarith
    exact_mod_cast this
  have hθ : Params.θ ⟨d, k, muOpt d k, ν⟩
      = 3 / ((d : ℚ) + 1) + (4 * (d : ℚ) - 2) / ((d : ℚ) + 1) * ν + ν := by
    rw [theta_at_muOpt d k ν hd hk1]; unfold thetaMin; rw [if_neg hcase]
  have hθ' : Params.θ ⟨d, k, muOpt d k, ν⟩ < 3 / ((d : ℚ) + 1) + 5 * ν := by
    rw [hθ]
    have : (4 * (d : ℚ) - 2) / ((d : ℚ) + 1) < 4 := by
      rw [div_lt_iff₀ (by linarith)]; linarith
    nlinarith
  refine ⟨⟨d, k, muOpt d k, ν⟩, rfl, ⟨hodd, hd, ?_, hν0, ?_, ?_, ?_⟩, ?_⟩
  · show ((d : ℚ) + 1) / 4 < k
    nlinarith
  · show ν + ν ≤ muOpt d k
    rw [hμ, le_div_iff₀ (by linarith)]
    have : (d : ℚ) + 3 ≤ 4 * k := by nlinarith
    nlinarith [mul_le_mul_of_nonneg_right this hν0.le]
  · show muOpt d k < 1
    rw [hμ, div_lt_one (by linarith)]
    nlinarith [mul_le_mul_of_nonneg_left hν36 (by linarith : (0 : ℚ) ≤ (d : ℚ) - 1)]
  · show Params.θ _ < 1
    have : 3 / ((d : ℚ) + 1) ≤ 3 / 4 := by
      rw [div_le_div_iff₀ (by linarith) (by norm_num)]; linarith
    linarith
  · unfold Params.cA
    rw [div_lt_iff₀ hs0]
    have h1 : ((Params.θ ⟨d, k, muOpt d k, ν⟩ : ℚ) : ℝ) < 3 / ((d : ℝ) + 1) + 5 * (1 / (4 * (k : ℝ))) := by
      have := (Rat.cast_lt (K := ℝ)).2 hθ'
      rw [hν] at this
      push_cast at this; exact this
    have hkN : (N : ℝ) ≤ k := by
      rw [hk]; push_cast; nlinarith [sq_nonneg (d : ℝ)]
    have hk0 : (0 : ℝ) < k := by
      have : (0 : ℚ) < k := by linarith
      exact_mod_cast this
    have h5 : 5 < ε * k := by
      rw [div_lt_iff₀ hε0] at hN; nlinarith
    have h2 : 5 * (1 / (4 * (k : ℝ))) < ε := by
      rw [mul_one_div, div_lt_iff₀ (by positivity)]; nlinarith
    linarith

/-- MAIN. [pp] Corollary cor:main, exponent part: for every `c > 1/√2` some admissible tuple has
`(1+θ)/√2 < c`. -/
theorem exists_admissible_cA_lt (c : ℝ) (hc : 1 / Real.sqrt 2 < c) :
    ∃ P : Params, P.Admissible ∧ P.cA < c := by
  have hs0 : 0 < Real.sqrt 2 := Real.sqrt_pos.2 (by norm_num)
  have hc' : 1 < c * Real.sqrt 2 := by rwa [div_lt_iff₀ hs0] at hc
  have hε : 0 < c * Real.sqrt 2 - 1 := by linarith
  obtain ⟨n, hn⟩ := exists_nat_gt (3 / (c * Real.sqrt 2 - 1))
  rw [div_lt_iff₀ hε] at hn
  have hlt : (1 + 3 / (((2 * n + 3 : ℕ) : ℝ) + 1)) / Real.sqrt 2 < c := by
    rw [div_lt_iff₀ hs0]
    have : 3 / (((2 * n + 3 : ℕ) : ℝ) + 1) < c * Real.sqrt 2 - 1 := by
      rw [div_lt_iff₀ (by positivity)]; push_cast; nlinarith [(n.cast_nonneg : (0 : ℝ) ≤ n)]
    linarith
  obtain ⟨P, -, hP, hPc⟩ := exists_admissible_fixed_d (2 * n + 3) (by omega) ⟨n + 1, by ring⟩ c hlt
  exact ⟨P, hP, hPc⟩

/-- MAIN. [pp] Table tab:constants, d = 9, k = 13: μ = 17/156, θ = 16/39 + η₀, admissible for
`0 < η₀ ≤ 1/100`, and the constant is below 1 exactly when `η₀ < √2 − 1 − 16/39`. -/
theorem d9_k13 (η₀ : ℚ) (h0 : 0 < η₀) (h1 : η₀ ≤ 1 / 100) :
    (Params.θ ⟨9, 13, 17 / 156, η₀⟩) = 16 / 39 + η₀ ∧ Params.Admissible ⟨9, 13, 17 / 156, η₀⟩ ∧
      (Params.cA ⟨9, 13, 17 / 156, η₀⟩ < 1 ↔ (η₀ : ℝ) < Real.sqrt 2 - 1 - 16 / 39) := by
  have hθ : Params.θ ⟨9, 13, 17 / 156, η₀⟩ = 16 / 39 + η₀ := by
    simp only [Params.θ, Params.h, Params.hStar, Params.ν]
    norm_num [max_def]
    ring
  refine ⟨hθ, ⟨⟨4, by norm_num⟩, by norm_num, by norm_num, h0, ?_, by norm_num, ?_⟩, ?_⟩
  · show 1 / (4 * ((13 : ℕ) : ℚ)) + η₀ ≤ 17 / 156
    norm_num; linarith
  · rw [hθ]; linarith
  · rw [cA_lt_one_iff _ (by rw [hθ]; push_cast; positivity), hθ]
    push_cast
    constructor <;> intro h <;> linarith

/-- MAIN. [pp] Table tab:constants: at `d = 9` no admissible tuple with `k ≤ 12` has constant
below 1. -/
theorem d9_k_le_12 (P : Params) (hP : P.Admissible) (hd : P.d = 9) (hk : P.k ≤ 12) :
    1 ≤ P.cA := by
  have hθ := theta_ge_AC P (by omega)
  obtain ⟨_, _, hk3, hη, _, _, _⟩ := hP
  rw [hd] at hθ hk3
  have hkq : (P.k : ℚ) ≤ 12 := by exact_mod_cast hk
  have hk0 : (0 : ℚ) < P.k := by norm_num at hk3; linarith
  have hν : (1 : ℚ) / 48 ≤ P.ν := by
    unfold Params.ν; rw [div_le_div_iff₀ (by norm_num) (by linarith)]; linarith
  have hθq : (23 : ℚ) / 54 < P.θ := by norm_num at hθ; nlinarith
  have hθr : (23 : ℝ) / 54 < (P.θ : ℝ) := by
    have := (Rat.cast_lt (K := ℝ)).2 hθq
    push_cast at this; exact this
  have hs : Real.sqrt 2 < 77 / 54 := by rw [Real.sqrt_lt' (by norm_num)]; norm_num
  unfold Params.cA
  rw [le_div_iff₀ (Real.sqrt_pos.2 (by norm_num))]
  linarith

/-- MAIN. [pp] Table tab:constants: at `d = 9` every `k ≥ 13` admits an admissible tuple with
constant below 1 (so `k ≥ 13` is necessary and sufficient). -/
theorem d9_k_ge_13 (k : ℕ) (hk : 13 ≤ k) :
    ∃ μ η₀ : ℚ, Params.Admissible ⟨9, k, μ, η₀⟩ ∧ Params.cA ⟨9, k, μ, η₀⟩ < 1 := by
  have hkq : (13 : ℚ) ≤ k := by exact_mod_cast hk
  have hk4 : (0 : ℚ) < 4 * (k : ℚ) := by linarith
  have hk1 : 1 ≤ k := by omega
  have ht : 1 / (4 * (k : ℚ)) ≤ 1 / 52 := by
    rw [div_le_div_iff₀ hk4 (by norm_num)]; linarith
  have ht0 : 0 < 1 / (4 * (k : ℚ)) := by positivity
  generalize htt : 1 / (4 * (k : ℚ)) = t at ht ht0
  have hθm : thetaMin 9 k ≤ 16 / 39 := by
    unfold thetaMin
    rw [htt]
    split_ifs <;> norm_num <;> linarith
  have hμ : t + 1 / 1000 ≤ muOpt 9 k := by
    unfold muOpt
    rw [htt]
    split_ifs
    · rw [le_div_iff₀ (by norm_num)]; norm_num; linarith
    · rw [le_div_iff₀ (by norm_num)]; norm_num; linarith
  have hθ : Params.θ ⟨9, k, muOpt 9 k, 1 / 1000⟩ = thetaMin 9 k + 1 / 1000 :=
    theta_at_muOpt 9 k (1 / 1000) (by norm_num) hk1
  have hθq : Params.θ ⟨9, k, muOpt 9 k, 1 / 1000⟩ ≤ 14113 / 10000 - 1 := by
    rw [hθ]; linarith
  have hθr : ((Params.θ ⟨9, k, muOpt 9 k, 1 / 1000⟩ : ℚ) : ℝ) ≤ 14113 / 10000 - 1 := by
    have := (Rat.cast_le (K := ℝ)).2 hθq
    push_cast at this; exact this
  have hs : (14114 / 10000 : ℝ) < Real.sqrt 2 := by
    rw [Real.lt_sqrt (by norm_num)]; norm_num
  refine ⟨muOpt 9 k, 1 / 1000, ⟨⟨4, by norm_num⟩, by norm_num, ?_, by norm_num, ?_,
    (muOpt_mem 9 k (by norm_num) ?_).2, ?_⟩, ?_⟩
  · show (((9 : ℕ) : ℚ) + 1) / 4 < (k : ℚ)
    norm_num; linarith
  · show 1 / (4 * (k : ℚ)) + 1 / 1000 ≤ muOpt 9 k
    rw [htt]; exact hμ
  · norm_num; linarith
  · show Params.θ ⟨9, k, muOpt 9 k, 1 / 1000⟩ < 1
    linarith
  · unfold Params.cA
    rw [div_lt_one (Real.sqrt_pos.2 (by norm_num))]
    linarith

/-! ### Real-analysis helpers -/

lemma sqrt2_facts : Real.sqrt 2 * Real.sqrt 2 = 2 ∧ 0 < Real.sqrt 2 ∧ 1 < Real.sqrt 2 ∧
    Real.sqrt 2 < 3 / 2 := by
  have h2 : Real.sqrt 2 * Real.sqrt 2 = 2 := Real.mul_self_sqrt (by norm_num)
  have hp : 0 < Real.sqrt 2 := Real.sqrt_pos.2 (by norm_num)
  refine ⟨h2, hp, ?_, ?_⟩ <;> nlinarith

/-- `a + θ/(2a) ≥ (1+θ)/√2` whenever `0 < a ≤ θ/√2` and `θ < 1`. -/
lemma lowerA (θ a : ℝ) (h1 : θ < 1) (ha : 0 < a) (hle : a ≤ θ / Real.sqrt 2) :
    (1 + θ) / Real.sqrt 2 ≤ a + θ / (2 * a) := by
  obtain ⟨hs, hs0, _, _⟩ := sqrt2_facts
  set s := Real.sqrt 2
  have hp : a * s ≤ θ := (le_div_iff₀ hs0).1 hle
  have key : 0 ≤ (θ - a * s) * (1 - a * s) := mul_nonneg (by linarith) (by linarith)
  have e : a + θ / (2 * a) = (2 * a ^ 2 + θ) / (2 * a) := by field_simp
  rw [e, div_le_div_iff₀ hs0 (by positivity)]
  have h3 : s * ((1 + θ) * (2 * a)) ≤ s * ((2 * a ^ 2 + θ) * s) := by
    have e2 : s * ((2 * a ^ 2 + θ) * s) = (s * s) * (2 * a ^ 2 + θ) := by ring
    rw [e2, hs]
    nlinarith
  exact le_of_mul_le_mul_left h3 hs0

lemma costA_lower (θ : ℝ) (h1 : θ < 1) (a : ℝ) (ha : 0 < a) :
    (1 + θ) / Real.sqrt 2 ≤ a + 1 / Real.sqrt 2 ∨ (1 + θ) / Real.sqrt 2 ≤ a + θ / (2 * a) := by
  obtain ⟨hs, hs0, _, _⟩ := sqrt2_facts
  by_cases h : θ / Real.sqrt 2 ≤ a
  · left
    rw [add_div]; linarith
  · right
    exact lowerA θ a h1 ha (le_of_lt (not_le.1 h))

lemma half_div (θ : ℝ) (h0 : 0 < θ) : θ / (2 * (θ / Real.sqrt 2)) = 1 / Real.sqrt 2 := by
  obtain ⟨hs, hs0, _, _⟩ := sqrt2_facts
  rw [div_eq_div_iff (by positivity) hs0.ne']
  field_simp
  nlinarith

/-- MAIN. [pp] Prop. prop:cost and the remark after it: for `0 < θ < 1` the cost exponent
`max(a + 1/√2, a + θ/(2a), 2a)` is minimised at `a = θ/√2`, with minimum `(1+θ)/√2`. -/
theorem costExpA_min (θ : ℝ) (h0 : 0 < θ) (h1 : θ < 1) :
    (∀ a > 0, (1 + θ) / Real.sqrt 2 ≤ costExpA θ a) ∧
      costExpA θ (θ / Real.sqrt 2) = (1 + θ) / Real.sqrt 2 := by
  obtain ⟨hs, hs0, _, _⟩ := sqrt2_facts
  have low : ∀ a > 0, (1 + θ) / Real.sqrt 2 ≤ costExpA θ a := by
    intro a ha
    unfold costExpA
    rcases costA_lower θ h1 a ha with h | h
    · exact le_max_of_le_left h
    · exact le_max_of_le_right (le_max_of_le_left h)
  refine ⟨low, le_antisymm ?_ (low _ (by positivity))⟩
  unfold costExpA
  rw [half_div θ h0]
  refine max_le (le_of_eq (by ring)) (max_le (le_of_eq (by ring)) ?_)
  rw [mul_div_assoc', div_le_div_iff_of_pos_right hs0]
  linarith

/-- MAIN. [pp] Theorem thm:RA, optimisation: for `θ > 0`, `c(a) ≥ √((1+√2)θ)` for all `a > 0`,
with equality at `a* = √((√2−1)θ/2)`. -/
theorem costExpRA_min (θ : ℝ) (h0 : 0 < θ) :
    (∀ a > 0, Real.sqrt ((1 + Real.sqrt 2) * θ) ≤ costExpRA θ a) ∧
      costExpRA θ (Real.sqrt ((Real.sqrt 2 - 1) * θ / 2)) = Real.sqrt ((1 + Real.sqrt 2) * θ) := by
  obtain ⟨hs, hs0, hs1, _⟩ := sqrt2_facts
  have hx : 0 < (Real.sqrt 2 - 1) * θ / 2 := by
    have : 0 < Real.sqrt 2 - 1 := by linarith
    positivity
  have ht0 : 0 < Real.sqrt ((Real.sqrt 2 - 1) * θ / 2) := Real.sqrt_pos.2 hx
  have htt : Real.sqrt ((Real.sqrt 2 - 1) * θ / 2) * Real.sqrt ((Real.sqrt 2 - 1) * θ / 2)
      = (Real.sqrt 2 - 1) * θ / 2 := Real.mul_self_sqrt hx.le
  set s := Real.sqrt 2
  set t := Real.sqrt ((s - 1) * θ / 2)
  have hV : Real.sqrt ((1 + s) * θ) = (2 + s) * t := by
    rw [Real.sqrt_eq_iff_mul_self_eq_of_pos (by positivity)]
    linear_combination (2 + s) ^ 2 * htt + θ / 2 * (s + 3) * hs
  have hθ : θ = 2 * (s + 1) * (t * t) := by
    linear_combination (-(2 : ℝ) * (s + 1)) * htt - θ * hs
  rw [hV]
  unfold costExpRA
  constructor
  · intro a ha
    by_cases h : t ≤ a
    · apply le_max_of_le_right
      exact mul_le_mul_of_nonneg_left h (by linarith)
    · apply le_max_of_le_left
      push_neg at h
      rw [← sub_nonneg]
      have e : a + θ / (2 * a) - (2 + s) * t = (a - t) * (a - (s + 1) * t) / a := by
        rw [hθ]; field_simp; ring
      rw [e]
      apply div_nonneg _ ha.le
      apply mul_nonneg_of_nonpos_of_nonpos <;> nlinarith
  · have e : t + θ / (2 * t) = (2 + s) * t := by
      rw [hθ]; field_simp; ring
    rw [e, max_self]

/-- MAIN. [pp] Remark rem:RA(a): for `0 ≤ θ`, `√((1+√2)θ) < 1 ↔ (1+θ)/√2 < 1`, and for
`0 < θ < 1`, `√((1+√2)θ) < (1+θ)/√2 ↔ θ < √2 − 1`. -/
theorem cRA_vs_cA (θ : ℝ) (h0 : 0 ≤ θ) :
    (Real.sqrt ((1 + Real.sqrt 2) * θ) < 1 ↔ (1 + θ) / Real.sqrt 2 < 1) ∧
      (0 < θ → θ < 1 →
        (Real.sqrt ((1 + Real.sqrt 2) * θ) < (1 + θ) / Real.sqrt 2 ↔ θ < Real.sqrt 2 - 1)) := by
  obtain ⟨hs, hs0, hs1, _⟩ := sqrt2_facts
  set s := Real.sqrt 2
  constructor
  · rw [Real.sqrt_lt' one_pos, div_lt_one hs0, one_pow]
    constructor
    · intro h
      by_contra hc; push_neg at hc
      nlinarith [mul_le_mul_of_nonneg_left hc (by linarith : (0 : ℝ) ≤ 1 + s)]
    · intro h
      nlinarith [mul_lt_mul_of_pos_left (by linarith : θ < s - 1) (by linarith : (0 : ℝ) < 1 + s)]
  · intro hθ0 hθ1
    rw [Real.sqrt_lt' (by positivity), div_pow, sq s, hs, lt_div_iff₀ (by norm_num : (0 : ℝ) < 2)]
    constructor
    · intro h
      by_contra hc; push_neg at hc
      nlinarith [mul_nonneg (by linarith : (0 : ℝ) ≤ θ - (s - 1)) (by linarith : (0 : ℝ) ≤ s + 1 - θ)]
    · intro h
      nlinarith [mul_pos (by linarith : (0 : ℝ) < s - 1 - θ) (by linarith : (0 : ℝ) < s + 1 - θ)]

lemma portfolio_fields (j : ℕ) :
    (portfolioParams j).d = 2 * j + 7 ∧ (portfolioParams j).k = j ^ 2 + 8 * j + 13 ∧
    (portfolioParams j).μ = muOpt (2 * j + 7) (j ^ 2 + 8 * j + 13) ∧
    (portfolioParams j).η₀ = 1 / ((2 * j + 7 : ℕ) : ℚ) ^ 2 := by
  have hD : (2 * j + 7) ^ 2 + 2 * (2 * j + 7) - 11 = 4 * (j ^ 2 + 8 * j + 13) := by
    have : (2 * j + 7) ^ 2 + 2 * (2 * j + 7) = 4 * (j ^ 2 + 8 * j + 13) + 11 := by ring
    rw [this, Nat.add_sub_cancel]
  have hcase : ¬ (4 * ((j ^ 2 + 8 * j + 13 : ℕ) : ℚ) < ((2 * j + 7 : ℕ) : ℚ) ^ 2 + 1) := by
    push_cast; push_neg; nlinarith [(j.cast_nonneg : (0 : ℚ) ≤ j)]
  refine ⟨rfl, ?_, ?_, rfl⟩
  · show ((2 * j + 7) ^ 2 + 2 * (2 * j + 7) - 11) / 4 = _
    rw [hD, Nat.mul_div_cancel_left _ (by norm_num)]
  · show (1 + (((2 * j + 7 : ℕ) : ℚ) - 1) / (((2 * j + 7) ^ 2 + 2 * (2 * j + 7) - 11 : ℕ) : ℚ)) /
      (((2 * j + 7 : ℕ) : ℚ) + 1) = _
    rw [hD, muOpt, if_neg hcase]
    push_cast; ring

lemma portfolio_key (j : ℕ) (hj : 1 ≤ j) : (portfolioParams j).Admissible ∧
    ((portfolioParams j).θ : ℝ) < 3 / ((2 * j + 7 : ℝ) + 1) + 5 / (2 * j + 7 : ℝ) ^ 2 ∧
    0 < ((portfolioParams j).θ : ℝ) := by
  obtain ⟨h1, h2, h3, h4⟩ := portfolio_fields j
  have eta : portfolioParams j = ⟨(portfolioParams j).d, (portfolioParams j).k,
      (portfolioParams j).μ, (portfolioParams j).η₀⟩ := rfl
  rw [eta, h1, h2, h3, h4]
  have hx : (1 : ℚ) ≤ j := by exact_mod_cast hj
  have hcase : ¬ (4 * ((j ^ 2 + 8 * j + 13 : ℕ) : ℚ) < ((2 * j + 7 : ℕ) : ℚ) ^ 2 + 1) := by
    push_cast; push_neg; nlinarith
  have hθ : Params.θ ⟨2 * j + 7, j ^ 2 + 8 * j + 13, muOpt (2 * j + 7) (j ^ 2 + 8 * j + 13),
      1 / ((2 * j + 7 : ℕ) : ℚ) ^ 2⟩ = 3 / ((2 * (j : ℚ) + 7) + 1) +
        (4 * (2 * (j : ℚ) + 7) - 2) / ((2 * (j : ℚ) + 7) + 1) *
          (1 / (4 * ((j : ℚ) ^ 2 + 8 * j + 13))) + 1 / (2 * (j : ℚ) + 7) ^ 2 := by
    rw [theta_at_muOpt _ _ _ (by omega) (by nlinarith)]
    unfold thetaMin; rw [if_neg hcase]; push_cast; ring
  have hμ : muOpt (2 * j + 7) (j ^ 2 + 8 * j + 13) = (1 + ((2 * (j : ℚ) + 7) - 1) *
      (1 / (4 * ((j : ℚ) ^ 2 + 8 * j + 13)))) / ((2 * (j : ℚ) + 7) + 1) := by
    unfold muOpt; rw [if_neg hcase]; push_cast; ring
  set Q : Params := ⟨2 * j + 7, j ^ 2 + 8 * j + 13, muOpt (2 * j + 7) (j ^ 2 + 8 * j + 13),
      1 / ((2 * j + 7 : ℕ) : ℚ) ^ 2⟩ with hQ
  set x : ℚ := (j : ℚ) with hxdef
  have hK : 0 < 4 * (x ^ 2 + 8 * x + 13) := by nlinarith
  have hd0 : 0 < 2 * x + 7 := by linarith
  have hd2 : 0 < (2 * x + 7) ^ 2 := by positivity
  have hsmall : (4 * (2 * x + 7) - 2) / ((2 * x + 7) + 1) * (1 / (4 * (x ^ 2 + 8 * x + 13)))
      < 4 / (2 * x + 7) ^ 2 := by
    rw [div_mul_div_comm, mul_one, div_lt_div_iff₀ (by positivity) hd2]
    nlinarith
  have hpos1 : 0 < (4 * (2 * x + 7) - 2) / ((2 * x + 7) + 1) * (1 / (4 * (x ^ 2 + 8 * x + 13))) :=
    mul_pos (div_pos (by linarith) (by linarith)) (by positivity)
  have hθlt : Q.θ < 3 / ((2 * x + 7) + 1) + 5 / (2 * x + 7) ^ 2 := by
    rw [hθ]
    have : 4 / (2 * x + 7) ^ 2 + 1 / (2 * x + 7) ^ 2 = 5 / (2 * x + 7) ^ 2 := by ring
    linarith
  have hθpos : 0 < Q.θ := by
    rw [hθ]; have : 0 < 3 / ((2 * x + 7) + 1) := by positivity
    have : 0 < 1 / (2 * x + 7) ^ 2 := by positivity
    linarith
  have hνd : 1 / (4 * (x ^ 2 + 8 * x + 13)) ≤ 1 / (2 * x + 7) ^ 2 := by
    rw [div_le_div_iff₀ hK hd2]; nlinarith
  have hν1 : 1 / (4 * (x ^ 2 + 8 * x + 13)) ≤ 1 := by
    rw [div_le_one hK]; nlinarith
  refine ⟨⟨⟨j + 3, by ring⟩, show 3 ≤ 2 * j + 7 by omega, ?_, ?_, ?_, ?_, ?_⟩, ?_, ?_⟩
  · show (((2 * j + 7 : ℕ) : ℚ) + 1) / 4 < ((j ^ 2 + 8 * j + 13 : ℕ) : ℚ)
    push_cast; nlinarith
  · show 0 < 1 / ((2 * j + 7 : ℕ) : ℚ) ^ 2
    push_cast; positivity
  · show 1 / (4 * ((j ^ 2 + 8 * j + 13 : ℕ) : ℚ)) + 1 / ((2 * j + 7 : ℕ) : ℚ) ^ 2 ≤
      muOpt (2 * j + 7) (j ^ 2 + 8 * j + 13)
    rw [hμ]; push_cast
    have hB : 2 / (2 * x + 7) ^ 2 ≤ 1 / ((2 * x + 7) + 1) := by
      rw [div_le_div_iff₀ hd2 (by linarith)]; nlinarith
    have hC : 1 / ((2 * x + 7) + 1) ≤ (1 + ((2 * x + 7) - 1) *
        (1 / (4 * (x ^ 2 + 8 * x + 13)))) / ((2 * x + 7) + 1) := by
      rw [div_le_div_iff_of_pos_right (by linarith)]
      have : 0 ≤ ((2 * x + 7) - 1) * (1 / (4 * (x ^ 2 + 8 * x + 13))) :=
        mul_nonneg (by linarith) (by positivity)
      linarith
    have : 1 / (2 * x + 7) ^ 2 + 1 / (2 * x + 7) ^ 2 = 2 / (2 * x + 7) ^ 2 := by ring
    linarith
  · show muOpt (2 * j + 7) (j ^ 2 + 8 * j + 13) < 1
    rw [hμ, div_lt_one (by linarith)]
    nlinarith [mul_le_mul_of_nonneg_left hν1 (by linarith : (0 : ℚ) ≤ (2 * x + 7) - 1)]
  · show Q.θ < 1
    have : 3 / ((2 * x + 7) + 1) ≤ 3 / 10 := by
      rw [div_le_div_iff₀ (by linarith) (by norm_num)]; linarith
    have : 5 / (2 * x + 7) ^ 2 ≤ 5 / 81 := by
      rw [div_le_div_iff₀ hd2 (by norm_num)]; nlinarith
    linarith
  · have := (Rat.cast_lt (K := ℝ)).2 hθlt
    push_cast at this; exact this
  · exact_mod_cast hθpos

/-- MAIN. [pp] Corollary cor:P: the portfolio parameter tuples are admissible for every `j ≥ 1`,
satisfy `θ_j < 3/(d_j+1) + 5/d_j²`, and `θ_j → 0`. -/
theorem portfolioParams_spec :
    (∀ j : ℕ, 1 ≤ j → (portfolioParams j).Admissible ∧
      ((portfolioParams j).θ : ℝ) < 3 / ((2 * j + 7 : ℝ) + 1) + 5 / (2 * j + 7 : ℝ) ^ 2) ∧
    Filter.Tendsto (fun j => ((portfolioParams j).θ : ℝ)) Filter.atTop (nhds 0) := by
  refine ⟨fun j hj => ⟨(portfolio_key j hj).1, (portfolio_key j hj).2.1⟩, ?_⟩
  have hlim : Filter.Tendsto (fun j : ℕ => 8 * (1 / ((j : ℝ) + 1))) Filter.atTop (nhds 0) := by
    have := (tendsto_one_div_add_atTop_nhds_zero_nat).const_mul (8 : ℝ)
    rwa [mul_zero] at this
  apply squeeze_zero' _ _ hlim
  · filter_upwards [Filter.eventually_ge_atTop 1] with j hj using (portfolio_key j hj).2.2.le
  · filter_upwards [Filter.eventually_ge_atTop 1] with j hj
    have := (portfolio_key j hj).2.1
    have hj0 : (0 : ℝ) ≤ j := Nat.cast_nonneg j
    have e1 : 3 / ((2 * j + 7 : ℝ) + 1) ≤ 3 / ((j : ℝ) + 1) := by
      rw [div_le_div_iff₀ (by positivity) (by positivity)]; nlinarith
    have e2 : 5 / (2 * j + 7 : ℝ) ^ 2 ≤ 5 / ((j : ℝ) + 1) := by
      rw [div_le_div_iff₀ (by positivity) (by positivity)]; nlinarith
    have e3 : 8 * (1 / ((j : ℝ) + 1)) = 3 / ((j : ℝ) + 1) + 5 / ((j : ℝ) + 1) := by ring
    linarith

/-- MAIN. For every `c > 0` some admissible tuple has `√((1+√2)θ) < c` (the exponent part of
[pp] Corollary cor:P). -/
theorem exists_admissible_cRA_lt (c : ℝ) (hc : 0 < c) :
    ∃ P : Params, P.Admissible ∧ P.cRAstar < c := by
  obtain ⟨hspec, hlim⟩ := portfolioParams_spec
  have hev := hlim.eventually (gt_mem_nhds (by positivity : (0 : ℝ) < c ^ 2 / 3))
  obtain ⟨j, hj1, hj2⟩ := ((Filter.eventually_ge_atTop 1).and hev).exists
  refine ⟨portfolioParams j, (hspec j hj1).1, ?_⟩
  unfold Params.cRAstar
  rw [Real.sqrt_lt' hc]
  have hs : Real.sqrt 2 < 2 := by rw [Real.sqrt_lt' (by norm_num)]; norm_num
  have hs0 : 0 ≤ Real.sqrt 2 := Real.sqrt_nonneg 2
  by_cases h : 0 ≤ ((portfolioParams j).θ : ℝ)
  · nlinarith
  · push_neg at h; nlinarith

/-- MAIN. [pp] Theorem thm:DL: `c_DL(θ) = min_a max(a + 1/√2, a + θ/(2a), 3a)` equals
`(1+θ)/√2` for `0 < θ ≤ 1/2` (at `a = θ/√2`) and `3√θ/2` for `1/2 < θ < 1` (at `a = √θ/2`). -/
theorem costExpDL_min (θ : ℝ) (h0 : 0 < θ) (h1 : θ < 1) :
    (θ ≤ 1 / 2 → (∀ a > 0, (1 + θ) / Real.sqrt 2 ≤ costExpDL θ a) ∧
        costExpDL θ (θ / Real.sqrt 2) = (1 + θ) / Real.sqrt 2) ∧
      (1 / 2 < θ → (∀ a > 0, 3 * Real.sqrt θ / 2 ≤ costExpDL θ a) ∧
        costExpDL θ (Real.sqrt θ / 2) = 3 * Real.sqrt θ / 2) := by
  obtain ⟨hs, hs0, _, _⟩ := sqrt2_facts
  constructor
  · intro hθ
    have low : ∀ a > 0, (1 + θ) / Real.sqrt 2 ≤ costExpDL θ a := by
      intro a ha
      unfold costExpDL
      rcases costA_lower θ h1 a ha with h | h
      · exact le_max_of_le_left h
      · exact le_max_of_le_right (le_max_of_le_left h)
    refine ⟨low, le_antisymm ?_ (low _ (by positivity))⟩
    unfold costExpDL
    rw [half_div θ h0]
    refine max_le (le_of_eq (by ring)) (max_le (le_of_eq (by ring)) ?_)
    rw [mul_div_assoc', div_le_div_iff_of_pos_right hs0]
    linarith
  · intro hθ
    have hr0 : 0 < Real.sqrt θ := Real.sqrt_pos.2 h0
    have hrr : Real.sqrt θ * Real.sqrt θ = θ := Real.mul_self_sqrt h0.le
    set r := Real.sqrt θ
    have low : ∀ a > 0, 3 * r / 2 ≤ costExpDL θ a := by
      intro a ha
      unfold costExpDL
      by_cases h : r / 2 ≤ a
      · exact le_max_of_le_right (le_max_of_le_right (by linarith))
      · push_neg at h
        apply le_max_of_le_right (le_max_of_le_left _)
        rw [← sub_nonneg]
        have e : a + θ / (2 * a) - 3 * r / 2 = (2 * a - r) * (a - r) / (2 * a) := by
          rw [← hrr]; field_simp; ring
        rw [e]
        apply div_nonneg _ (by linarith)
        apply mul_nonneg_of_nonpos_of_nonpos <;> linarith
    refine ⟨low, le_antisymm ?_ (low _ (by positivity))⟩
    unfold costExpDL
    have e : r / 2 + θ / (2 * (r / 2)) = 3 * r / 2 := by
      rw [← hrr]; field_simp; ring
    rw [e]
    refine max_le ?_ (max_le le_rfl (le_of_eq (by ring)))
    have : 1 / Real.sqrt 2 ≤ r := by
      rw [div_le_iff₀ hs0]
      by_contra hc; push_neg at hc
      nlinarith [mul_pos hr0 hs0]
    linarith

/-- MAIN. Arithmetic step of the size-model floor (recurse §3.3, lv5unif §13): with precision
exponent `τ ≥ 0`, the constraints `k ≳ D^{2(1−2α)}` and `k·d^τ = O(D^α)` are compatible iff
`(2+τ)(1−2α) ≤ α`. Those constraints themselves are heuristic and are not formalised here; only
the arithmetic equivalence `(2+τ)(1−2α) ≤ α ↔ α ≥ (2+τ)/(5+2τ)` (and its `τ = 0` case
`α ≥ 2/5`) is. -/
theorem size_model_floor (τ α : ℝ) (hτ : 0 ≤ τ) :
    ((2 + τ) * (1 - 2 * α) ≤ α ↔ (2 + τ) / (5 + 2 * τ) ≤ α) ∧
      ((2 : ℝ) * (1 - 2 * α) ≤ α ↔ 2 / 5 ≤ α) := by
  constructor
  · rw [div_le_iff₀ (by linarith)]
    constructor <;> intro h <;> nlinarith
  · constructor <;> intro h <;> linarith

/-- MAIN. lv5unif §10: `4.12981 ≤ 1/(1 − 2^{−2/5}) ≤ 4.12983`. -/
theorem recursion_factor_two_fifths :
    (4.12981 : ℝ) ≤ 1 / (1 - (2 : ℝ) ^ (-(2 / 5 : ℝ))) ∧
      1 / (1 - (2 : ℝ) ^ (-(2 / 5 : ℝ))) ≤ 4.12983 := by
  have hx0 : 0 < (2 : ℝ) ^ (-(2 / 5 : ℝ)) := Real.rpow_pos_of_pos (by norm_num) _
  have hx5 : ((2 : ℝ) ^ (-(2 / 5 : ℝ))) ^ (5 : ℕ) = 1 / 4 := by
    rw [← Real.rpow_natCast, ← Real.rpow_mul (by norm_num)]
    rw [show (-(2 / 5 : ℝ)) * ((5 : ℕ) : ℝ) = ((-2 : ℤ) : ℝ) by push_cast; ring, Real.rpow_intCast]
    norm_num
  set x := (2 : ℝ) ^ (-(2 / 5 : ℝ))
  have hlo : 1 - 1 / (4.12981 : ℝ) ≤ x := by
    by_contra h; push_neg at h
    have := pow_lt_pow_left₀ h hx0.le (by norm_num : (5 : ℕ) ≠ 0)
    rw [hx5] at this
    norm_num at this
  have hhi : x ≤ 1 - 1 / (4.12983 : ℝ) := by
    by_contra h; push_neg at h
    have := pow_lt_pow_left₀ h (by norm_num) (by norm_num : (5 : ℕ) ≠ 0)
    rw [hx5] at this
    norm_num at this
  have hpos : 0 < 1 - x := by norm_num at hhi; linarith
  constructor
  · rw [le_div_iff₀ hpos]
    norm_num at hlo ⊢; linarith
  · rw [div_le_iff₀ hpos]
    norm_num at hhi ⊢; linarith

end FactoringA
