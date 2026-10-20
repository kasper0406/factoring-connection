import LeanA.DefsNFS

/-!
# Elementary asymptotic helpers for the relation-supply argument

* logs of floors and ceilings;
* the `L`-scale `Λ(N) = Lam (1/2) N = √(log N · log log N)` (eventually) and its comparison with
  `log N` and `log log N`;
* a lower bound for `ϱ(X, B)` from CEP in the range used by [pp] Lemma lem:local(iv).
-/

noncomputable section
open Classical

namespace FactoringA

open Real Filter

/-! ### Floors and ceilings -/

lemma one_le_floor_real {x : ℝ} (hx : 1 ≤ x) : (1 : ℝ) ≤ (⌊x⌋₊ : ℝ) := by
  exact_mod_cast Nat.one_le_floor_iff x |>.mpr hx

lemma half_le_floor {x : ℝ} (hx : 1 ≤ x) : x / 2 ≤ (⌊x⌋₊ : ℝ) := by
  have h1 := Nat.lt_floor_add_one x
  have h2 := one_le_floor_real hx
  linarith

lemma log_floor_ge {x : ℝ} (hx : 1 ≤ x) : Real.log x - Real.log 2 ≤ Real.log (⌊x⌋₊ : ℝ) := by
  have hpos : 0 < x / 2 := by linarith
  calc Real.log x - Real.log 2 = Real.log (x / 2) := by
        rw [Real.log_div (by linarith) (by norm_num)]
    _ ≤ _ := Real.log_le_log hpos (half_le_floor hx)

lemma log_floor_le {x : ℝ} (hx : 1 ≤ x) : Real.log (⌊x⌋₊ : ℝ) ≤ Real.log x :=
  Real.log_le_log (by linarith [one_le_floor_real hx]) (Nat.floor_le (by linarith))

lemma log_ceil_ge {x : ℝ} (hx : 1 ≤ x) : Real.log x ≤ Real.log (⌈x⌉₊ : ℝ) :=
  Real.log_le_log (by linarith) (Nat.le_ceil x)

lemma ceil_le_two_mul {x : ℝ} (hx : 1 ≤ x) : (⌈x⌉₊ : ℝ) ≤ 2 * x := by
  have := Nat.ceil_lt_add_one (by linarith : (0:ℝ) ≤ x)
  linarith

lemma log_ceil_le {x : ℝ} (hx : 1 ≤ x) : Real.log (⌈x⌉₊ : ℝ) ≤ Real.log x + Real.log 2 := by
  have h0 : (0:ℝ) < ⌈x⌉₊ := by linarith [Nat.le_ceil x]
  calc Real.log (⌈x⌉₊ : ℝ) ≤ Real.log (2 * x) := Real.log_le_log h0 (ceil_le_two_mul hx)
    _ = _ := by rw [Real.log_mul (by norm_num) (by linarith)]; ring

/-! ### The scale `Λ = Lam (1/2)` -/

lemma Lam_half_nonneg (x : ℝ) : 0 ≤ Lam (1 / 2) x := by
  unfold Lam
  rcases le_or_gt 0 (Real.log x) with h | h
  · have : 0 ≤ max 1 (Real.log (Real.log x)) := le_trans zero_le_one (le_max_left _ _)
    positivity
  · -- negative base: rpow of a negative number is `exp(..) * cos(..)`; avoid by noting `log x`
    -- appears only through `x ≥ 1` in uses. We still need nonnegativity in general.
    have hx : Real.log x ^ ((1:ℝ) / 2) = Real.exp (Real.log (Real.log x) * (1 / 2)) *
        Real.cos ((1 / 2) * Real.pi) := Real.rpow_def_of_neg h _
    rw [hx]
    have : Real.cos ((1 / 2) * Real.pi) = 0 := by
      rw [show (1:ℝ) / 2 * Real.pi = Real.pi / 2 by ring]; exact Real.cos_pi_div_two
    rw [this]; simp

lemma Lam_half_sq {x : ℝ} (hx : 0 ≤ Real.log x) :
    Lam (1 / 2) x ^ 2 = Real.log x * max 1 (Real.log (Real.log x)) := by
  unfold Lam
  have hm : 0 ≤ max 1 (Real.log (Real.log x)) := le_trans zero_le_one (le_max_left _ _)
  rw [mul_pow, ← Real.rpow_natCast, ← Real.rpow_natCast, ← Real.rpow_mul hx, ← Real.rpow_mul hm]
  norm_num

/-- The basic little-o fact: `log t ≤ ε t` eventually. -/
lemma ev_log_le_mul (ε : ℝ) (hε : 0 < ε) : ∀ᶠ t : ℝ in atTop, Real.log t ≤ ε * t := by
  have h := (Real.isLittleO_log_id_atTop).bound hε
  filter_upwards [h, eventually_ge_atTop (1:ℝ)] with t ht ht1
  have h0 : 0 ≤ Real.log t := Real.log_nonneg ht1
  simp only [Real.norm_eq_abs, id] at ht
  rw [abs_of_nonneg h0, abs_of_nonneg (by linarith)] at ht
  exact ht

lemma tendsto_log_nat : Tendsto (fun N : ℕ => Real.log (N : ℝ)) atTop atTop :=
  Real.tendsto_log_atTop.comp tendsto_natCast_atTop_atTop

lemma ev_log_ge (C : ℝ) : ∀ᶠ N : ℕ in atTop, C ≤ Real.log (N : ℝ) :=
  tendsto_log_nat.eventually_ge_atTop C

lemma ev_loglog_ge (C : ℝ) : ∀ᶠ N : ℕ in atTop, C ≤ Real.log (Real.log (N : ℝ)) :=
  (Real.tendsto_log_atTop.comp tendsto_log_nat).eventually_ge_atTop C

lemma ev_nat_ge (C : ℝ) : ∀ᶠ N : ℕ in atTop, C ≤ (N : ℝ) :=
  tendsto_natCast_atTop_atTop.eventually_ge_atTop C

/-- Eventually `Λ² = t log t` with `t = log N`. -/
lemma ev_Lam_sq : ∀ᶠ N : ℕ in atTop,
    Lam (1 / 2) (N : ℝ) ^ 2 = Real.log N * Real.log (Real.log N) := by
  filter_upwards [ev_log_ge 0, ev_loglog_ge 1] with N h0 h1
  rw [Lam_half_sq h0, max_eq_right h1]

/-- `Λ ≥ C` eventually. -/
lemma ev_Lam_ge (C : ℝ) : ∀ᶠ N : ℕ in atTop, C ≤ Lam (1 / 2) (N : ℝ) := by
  filter_upwards [ev_Lam_sq, ev_log_ge (C ^ 2), ev_log_ge 0, ev_loglog_ge 1] with N hsq hC h0 h1
  have hL := Lam_half_nonneg (N : ℝ)
  have : C ^ 2 ≤ Lam (1 / 2) (N : ℝ) ^ 2 := by rw [hsq]; nlinarith
  rcases le_or_gt C 0 with hc | hc
  · linarith
  · nlinarith

/-- `Λ ≤ ε log N` eventually. -/
lemma ev_Lam_le (ε : ℝ) (hε : 0 < ε) : ∀ᶠ N : ℕ in atTop,
    Lam (1 / 2) (N : ℝ) ≤ ε * Real.log N := by
  filter_upwards [ev_Lam_sq, tendsto_log_nat.eventually (ev_log_le_mul (ε ^ 2) (by positivity)),
    ev_log_ge 0] with N hsq hl h0
  have hL := Lam_half_nonneg (N : ℝ)
  have : Lam (1 / 2) (N : ℝ) ^ 2 ≤ (ε * Real.log N) ^ 2 := by
    rw [hsq]; nlinarith
  exact (pow_le_pow_iff_left₀ hL (by positivity) two_ne_zero).mp this

/-- `log log N ≤ ε Λ` eventually. -/
lemma ev_loglog_le (ε : ℝ) (hε : 0 < ε) : ∀ᶠ N : ℕ in atTop,
    Real.log (Real.log N) ≤ ε * Lam (1 / 2) (N : ℝ) := by
  filter_upwards [ev_Lam_sq, tendsto_log_nat.eventually (ev_log_le_mul (ε ^ 2) (by positivity)),
    ev_log_ge 0, ev_loglog_ge 1] with N hsq hl h0 h1
  have hL := Lam_half_nonneg (N : ℝ)
  have : Real.log (Real.log N) ^ 2 ≤ (ε * Lam (1 / 2) (N : ℝ)) ^ 2 := by
    rw [mul_pow, hsq]; nlinarith
  have hpos : 0 ≤ ε * Lam (1 / 2) (N : ℝ) := by positivity
  exact (pow_le_pow_iff_left₀ (by linarith) hpos two_ne_zero).mp this

/-! ### Admissible parameters, as real numbers -/

namespace Params

variable {P : Params}

lemma adm_k (hP : P.Admissible) : 2 ≤ P.k := by
  obtain ⟨_, hd3, hk, _⟩ := hP
  have h3 : (3 : ℚ) ≤ P.d := by exact_mod_cast hd3
  have : (1 : ℚ) < P.k := by linarith
  have : 1 < P.k := by exact_mod_cast this
  omega

lemma nu_real (P : Params) : ((P.ν : ℚ) : ℝ) = 1 / (4 * (P.k : ℝ)) := by
  unfold Params.ν; push_cast; ring

lemma theta_real (P : Params) :
    ((P.θ : ℚ) : ℝ) = 2 * (P.μ : ℝ) + (P.h : ℝ) + ((P.d : ℝ) - 1) * (P.ν : ℝ) := by
  unfold Params.θ; push_cast; ring

lemma h_ge (P : Params) : 2 * (P.ν : ℝ) + P.η₀ ≤ (P.h : ℝ) := by
  have : 2 * P.ν ≤ P.hStar := le_max_left _ _
  have h2 : 2 * P.ν + P.η₀ ≤ P.h := by unfold Params.h; linarith
  exact_mod_cast h2

/-- The real facts about the exponents used below. -/
lemma adm_real (hP : P.Admissible) :
    0 < (P.ν : ℝ) ∧ (P.ν : ℝ) * (4 * P.k) = 1 ∧ 0 < (P.μ : ℝ) - P.ν ∧ 0 < (P.μ : ℝ) ∧
      0 < (P.h : ℝ) ∧ 0 < (P.θ : ℝ) ∧ (P.μ : ℝ) < 1 ∧ 3 ≤ (P.d : ℝ) := by
  have hk := adm_k hP
  obtain ⟨_, hd3, _, hη, hμ, hμ1, _⟩ := hP
  have hk0 : (0 : ℝ) < P.k := by exact_mod_cast (by omega : 0 < P.k)
  have hν : (0 : ℝ) < P.ν := by rw [nu_real]; positivity
  have hη' : (0 : ℝ) < P.η₀ := by exact_mod_cast hη
  have hμ' : (P.ν : ℝ) + P.η₀ ≤ P.μ := by exact_mod_cast hμ
  have hh := h_ge P
  have hd : (3 : ℝ) ≤ P.d := by exact_mod_cast hd3
  refine ⟨hν, ?_, by linarith, by linarith, by linarith, ?_, by exact_mod_cast hμ1, hd⟩
  · rw [nu_real]; field_simp
  · rw [theta_real]; nlinarith

end Params

/-! ### Sizes of the `N`-dependent quantities -/

lemma Lv_rpow (N : ℕ) (c : ℝ) : (Params.Lv N) ^ c = Real.exp (c * Lam (1 / 2) (N : ℝ)) := by
  unfold Params.Lv LN
  rw [one_mul, ← Real.exp_mul, mul_comm]

/-- `B = ⌊L^a⌋`: `log B = aΛ + O(1)` and `B ≥ 2`, eventually. -/
lemma ev_Bsm {a : ℝ} (ha : 0 < a) : ∀ᶠ N : ℕ in atTop,
    2 ≤ Params.Bsm a N ∧ a * Lam (1 / 2) (N : ℝ) - Real.log 2 ≤ Real.log (Params.Bsm a N) ∧
      Real.log (Params.Bsm a N) ≤ a * Lam (1 / 2) (N : ℝ) := by
  filter_upwards [ev_Lam_ge (Real.log 2 / a)] with N hΛ
  have hΛ' : Real.log 2 ≤ a * Lam (1 / 2) (N : ℝ) := by
    rw [div_le_iff₀ ha] at hΛ; linarith
  have h2 : (2 : ℝ) ≤ Real.exp (a * Lam (1 / 2) (N : ℝ)) := by
    calc (2 : ℝ) = Real.exp (Real.log 2) := (Real.exp_log (by norm_num)).symm
      _ ≤ _ := Real.exp_le_exp.mpr hΛ'
  unfold Params.Bsm
  rw [Lv_rpow]
  refine ⟨Nat.le_floor (by exact_mod_cast h2), ?_, ?_⟩
  · have := log_floor_ge (by linarith : (1 : ℝ) ≤ Real.exp (a * Lam (1 / 2) (N : ℝ)))
    rwa [Real.log_exp] at this
  · have := log_floor_le (by linarith : (1 : ℝ) ≤ Real.exp (a * Lam (1 / 2) (N : ℝ)))
    rwa [Real.log_exp] at this

/-- `S = ⌈L^{θ/(4a)}⌉`: `S² ≥ L^{θ/(2a)}`, `S ≥ 1`, `log S ≤ θΛ/(4a) + log 2`. -/
lemma Ssz_bounds (P : Params) {a : ℝ} (ha : 0 < a) (hθ : 0 ≤ (P.θ : ℝ)) (N : ℕ) :
    1 ≤ (P.Ssz a N : ℝ) ∧
      Real.exp ((P.θ : ℝ) / (2 * a) * Lam (1 / 2) (N : ℝ)) ≤ (P.Ssz a N : ℝ) ^ 2 ∧
      (P.θ : ℝ) / (4 * a) * Lam (1 / 2) (N : ℝ) ≤ Real.log (P.Ssz a N) ∧
      Real.log (P.Ssz a N) ≤ (P.θ : ℝ) / (4 * a) * Lam (1 / 2) (N : ℝ) + Real.log 2 := by
  unfold Params.Ssz
  rw [Lv_rpow]
  set z := (P.θ : ℝ) / (4 * a) * Lam (1 / 2) (N : ℝ)
  have hz : 0 ≤ z := by have := Lam_half_nonneg (N : ℝ); positivity
  have h1 : (1 : ℝ) ≤ Real.exp z := Real.one_le_exp hz
  have hc := Nat.le_ceil (Real.exp z)
  refine ⟨le_trans h1 hc, ?_, ?_, ?_⟩
  · have e : (P.θ : ℝ) / (2 * a) * Lam (1 / 2) (N : ℝ) = 2 * z := by
      simp only [z]; field_simp; ring
    rw [e, show (2 : ℝ) * z = z + z by ring, Real.exp_add, sq]
    exact mul_le_mul hc hc (by positivity) (by positivity)
  · have := log_ceil_ge h1; rwa [Real.log_exp] at this
  · have := log_ceil_le h1; rwa [Real.log_exp] at this

lemma log_rpow_nat {N : ℕ} (hN : 1 ≤ N) (c : ℝ) :
    Real.log ((N : ℝ) ^ c) = c * Real.log N :=
  Real.log_rpow (by exact_mod_cast (by omega : 0 < N)) c

/-- `log ⌊N^c⌋ = c log N + O(1)` for `c ≥ 0`, `N ≥ 1`. -/
lemma log_floor_rpow {N : ℕ} (hN : 1 ≤ N) {c : ℝ} (hc : 0 ≤ c) :
    1 ≤ (⌊(N : ℝ) ^ c⌋₊ : ℝ) ∧
    c * Real.log N - Real.log 2 ≤ Real.log (⌊(N : ℝ) ^ c⌋₊ : ℝ) ∧
      Real.log (⌊(N : ℝ) ^ c⌋₊ : ℝ) ≤ c * Real.log N := by
  have h1 : (1 : ℝ) ≤ (N : ℝ) ^ c := Real.one_le_rpow (by exact_mod_cast hN) hc
  refine ⟨one_le_floor_real h1, ?_, ?_⟩
  · have := log_floor_ge h1; rwa [log_rpow_nat hN] at this
  · have := log_floor_le h1; rwa [log_rpow_nat hN] at this

/-- `X = ⌊√(N/2)⌋`: `log X = (1/2) log N + O(1)` and `X ≥ C`, eventually. -/
lemma ev_Xsp (C : ℝ) : ∀ᶠ N : ℕ in atTop,
    C ≤ (Params.Xsp N : ℝ) ∧ 1 ≤ (Params.Xsp N : ℝ) ∧
      Real.log N / 2 - 2 * Real.log 2 ≤ Real.log (Params.Xsp N) ∧
      Real.log (Params.Xsp N) ≤ Real.log N / 2 := by
  filter_upwards [ev_nat_ge (2 * (max C 1 + 1) ^ 2)] with N hN
  have hC1 : 1 ≤ max C 1 := le_max_right _ _
  have hN0 : (0 : ℝ) ≤ N := by positivity
  have hNpos : (0 : ℝ) < N := by
    have : 0 < (max C 1 + 1) ^ 2 := by positivity
    linarith
  have hs : max C 1 + 1 ≤ Real.sqrt ((N : ℝ) / 2) := by
    rw [Real.le_sqrt (by linarith) (by positivity)]; linarith
  have hs1 : (1 : ℝ) ≤ Real.sqrt ((N : ℝ) / 2) := by linarith
  unfold Params.Xsp
  refine ⟨?_, one_le_floor_real hs1, ?_, ?_⟩
  · have := Nat.lt_floor_add_one (Real.sqrt ((N : ℝ) / 2))
    have := le_max_left C 1
    linarith
  · have h := log_floor_ge hs1
    have e : Real.log (Real.sqrt ((N : ℝ) / 2)) = (Real.log N - Real.log 2) / 2 := by
      rw [Real.log_sqrt (by positivity), Real.log_div hNpos.ne' (by norm_num)]
    rw [e] at h
    have : 0 < Real.log 2 := Real.log_pos (by norm_num)
    linarith
  · have h := log_floor_le hs1
    have e : Real.log (Real.sqrt ((N : ℝ) / 2)) = (Real.log N - Real.log 2) / 2 := by
      rw [Real.log_sqrt (by positivity), Real.log_div hNpos.ne' (by norm_num)]
    rw [e] at h
    have : 0 < Real.log 2 := Real.log_pos (by norm_num)
    linarith

/-- If `X ≥ 2^k` then `j ≥ 1`, `(2Y)^k ≤ X < (4Y)^k`. -/
lemma Yj_spec (P : Params) (hk : 1 ≤ P.k) (N : ℕ) (hX : 2 ^ P.k ≤ Params.Xsp N) :
    1 ≤ P.jJ N ∧ (2 * P.Yj N) ^ P.k ≤ Params.Xsp N ∧ Params.Xsp N < (4 * P.Yj N) ^ P.k := by
  have hb : 1 < 2 ^ P.k := Nat.one_lt_two_pow (by omega)
  have hj : 1 ≤ P.jJ N := by
    unfold Params.jJ; exact Nat.log_pos hb hX
  have e2 : 2 * P.Yj N = 2 ^ P.jJ N := by
    unfold Params.Yj
    rw [← pow_succ']; congr 1; omega
  have e4 : 4 * P.Yj N = 2 ^ (P.jJ N + 1) := by
    rw [show 4 * P.Yj N = 2 * (2 * P.Yj N) by ring, e2, pow_succ']
  refine ⟨hj, ?_, ?_⟩
  · rw [e2, ← pow_mul, mul_comm, pow_mul]
    exact Nat.pow_log_le_self _ (by omega)
  · rw [e4, ← pow_mul, mul_comm, pow_mul]
    exact Nat.lt_pow_succ_log_self hb _

/-- `Y = 2^{j−1}`: `j ≥ 1`, `(2Y)^k ≤ X < (4Y)^k`, `log Y = 2ν log N + O(1)`, eventually. -/
lemma ev_Yj {P : Params} (hP : P.Admissible) : ∀ᶠ N : ℕ in atTop,
    1 ≤ P.jJ N ∧ (2 * P.Yj N) ^ P.k ≤ Params.Xsp N ∧ Params.Xsp N < (4 * P.Yj N) ^ P.k ∧
      1 ≤ (P.Yj N : ℝ) ∧
      2 * (P.ν : ℝ) * Real.log N - 3 * Real.log 2 ≤ Real.log (P.Yj N) ∧
      Real.log (P.Yj N) ≤ 2 * (P.ν : ℝ) * Real.log N - Real.log 2 := by
  have hk := Params.adm_k hP
  obtain ⟨hν, hνk, -⟩ := Params.adm_real hP
  filter_upwards [ev_Xsp ((2 : ℝ) ^ P.k)] with N hXs
  obtain ⟨hXC, hX1, hXlo, hXhi⟩ := hXs
  have hX : 2 ^ P.k ≤ Params.Xsp N := by exact_mod_cast hXC
  obtain ⟨hj, h2Y, h4Y⟩ := Yj_spec P (by omega) N hX
  have hY1 : 1 ≤ P.Yj N := Nat.one_le_two_pow
  have hY1' : (1 : ℝ) ≤ P.Yj N := by exact_mod_cast hY1
  have hk0 : (0 : ℝ) < P.k := by exact_mod_cast (by omega : 0 < P.k)
  have hk2 : (2 : ℝ) ≤ P.k := by exact_mod_cast hk
  have hl2 : 0 < Real.log 2 := Real.log_pos (by norm_num)
  have hlY2 : Real.log (2 * (P.Yj N : ℝ)) = Real.log 2 + Real.log (P.Yj N) :=
    Real.log_mul (by norm_num) (by linarith)
  have hlY4 : Real.log (4 * (P.Yj N : ℝ)) = 2 * Real.log 2 + Real.log (P.Yj N) := by
    rw [Real.log_mul (by norm_num) (by linarith),
      show (4 : ℝ) = 2 ^ 2 by norm_num, Real.log_pow]; push_cast; ring
  refine ⟨hj, h2Y, h4Y, hY1', ?_, ?_⟩
  · have h : ((Params.Xsp N : ℕ) : ℝ) < ((4 * P.Yj N : ℕ) : ℝ) ^ P.k := by exact_mod_cast h4Y
    have h' := Real.log_lt_log (by linarith) h
    rw [Real.log_pow, Nat.cast_mul, Nat.cast_ofNat, hlY4] at h'
    -- `k (2 log 2 + log Y) > log N/2 − 2 log 2 = k·2ν log N − 2 log 2`
    have e : Real.log N / 2 = (P.k : ℝ) * (2 * (P.ν : ℝ) * Real.log N) := by
      linear_combination (-(Real.log N) / 2) * hνk
    have h3 : (P.k : ℝ) * (2 * (P.ν : ℝ) * Real.log N - 3 * Real.log 2) ≤
        (P.k : ℝ) * Real.log (P.Yj N) := by nlinarith
    exact le_of_mul_le_mul_left h3 hk0
  · have h : ((2 * P.Yj N : ℕ) : ℝ) ^ P.k ≤ ((Params.Xsp N : ℕ) : ℝ) := by exact_mod_cast h2Y
    have hpos : (0 : ℝ) < ((2 * P.Yj N : ℕ) : ℝ) ^ P.k := by positivity
    have h' := Real.log_le_log hpos h
    rw [Real.log_pow, Nat.cast_mul, Nat.cast_ofNat, hlY2] at h'
    have e : Real.log N / 2 = (P.k : ℝ) * (2 * (P.ν : ℝ) * Real.log N) := by
      linear_combination (-(Real.log N) / 2) * hνk
    have h3 : (P.k : ℝ) * (Real.log 2 + Real.log (P.Yj N)) ≤
        (P.k : ℝ) * (2 * (P.ν : ℝ) * Real.log N) := by linarith
    have := le_of_mul_le_mul_left h3 hk0
    linarith

/-- `V_max = ⌊√(2Y) S⌋`: `log V_max = ν log N + θΛ/(4a) + O(1)`, eventually. -/
lemma ev_Vmax {P : Params} (hP : P.Admissible) {a : ℝ} (ha : 0 < a) : ∀ᶠ N : ℕ in atTop,
    1 ≤ (P.Vmax a N : ℝ) ∧
      (P.ν : ℝ) * Real.log N + (P.θ : ℝ) / (4 * a) * Lam (1 / 2) (N : ℝ) - 3 * Real.log 2 ≤
        Real.log (P.Vmax a N) ∧
      Real.log (P.Vmax a N) ≤
        (P.ν : ℝ) * Real.log N + (P.θ : ℝ) / (4 * a) * Lam (1 / 2) (N : ℝ) + Real.log 2 := by
  obtain ⟨-, -, -, -, -, hθ, -⟩ := Params.adm_real hP
  filter_upwards [ev_Yj hP] with N hYs
  obtain ⟨-, -, -, hY1, hYlo, hYhi⟩ := hYs
  obtain ⟨hS1, -, hSlo, hShi⟩ := Ssz_bounds P ha hθ.le N
  have hl2 : 0 < Real.log 2 := Real.log_pos (by norm_num)
  have h2Y : (1 : ℝ) ≤ 2 * (P.Yj N : ℝ) := by linarith
  have hsq1 : (1 : ℝ) ≤ Real.sqrt (2 * (P.Yj N : ℝ)) := by
    rw [Real.le_sqrt (by norm_num) (by linarith)]; linarith
  have hw1 : (1 : ℝ) ≤ Real.sqrt (2 * (P.Yj N : ℝ)) * (P.Ssz a N : ℝ) := by nlinarith
  have hlw : Real.log (Real.sqrt (2 * (P.Yj N : ℝ)) * (P.Ssz a N : ℝ)) =
      (Real.log 2 + Real.log (P.Yj N)) / 2 + Real.log (P.Ssz a N) := by
    rw [Real.log_mul (by positivity) (by positivity), Real.log_sqrt (by positivity),
      Real.log_mul (by norm_num) (by linarith)]
  unfold Params.Vmax
  refine ⟨one_le_floor_real hw1, ?_, ?_⟩
  · have := log_floor_ge hw1
    rw [hlw] at this; linarith
  · have := log_floor_le hw1
    rw [hlw] at this; linarith

/-- `X_r = 2 M V_max / Y`: `log X_r = (μ−ν) log N + O(Λ)`, eventually. -/
lemma ev_Xr {P : Params} (hP : P.Admissible) {a : ℝ} (ha : 0 < a) : ∀ᶠ N : ℕ in atTop,
    0 < P.Xr a N ∧ ((P.μ : ℝ) - P.ν) / 2 * Real.log N ≤ Real.log (P.Xr a N) ∧
      Real.log (P.Xr a N) ≤
        ((P.μ : ℝ) - P.ν) * Real.log N + ((P.θ : ℝ) / (4 * a) + 1) * Lam (1 / 2) (N : ℝ) := by
  obtain ⟨hν, -, hμν, hμ, -, hθ, -⟩ := Params.adm_real hP
  have hl2 : 0 < Real.log 2 := Real.log_pos (by norm_num)
  filter_upwards [ev_Yj hP, ev_Vmax hP ha, ev_nat_ge 1, ev_Lam_ge (5 * Real.log 2),
    ev_log_ge (4 * Real.log 2 / ((P.μ : ℝ) - P.ν))] with
    N hYs hVs hN1 hΛ hlogN
  obtain ⟨-, -, -, hY1, hYlo, hYhi⟩ := hYs
  obtain ⟨hV1, hVlo, hVhi⟩ := hVs
  have hN1' : 1 ≤ N := by exact_mod_cast hN1
  obtain ⟨hM1, hMlo, hMhi⟩ := log_floor_rpow hN1' hμ.le
  have hΛ0 := Lam_half_nonneg (N : ℝ)
  have hθΛ : 0 ≤ (P.θ : ℝ) / (4 * a) * Lam (1 / 2) (N : ℝ) := by positivity
  have hlog : Real.log (P.Xr a N) =
      Real.log 2 + Real.log (P.Mb N) + Real.log (P.Vmax a N) - Real.log (P.Yj N) := by
    unfold Params.Xr
    unfold Params.Mb at *
    rw [Real.log_div (by positivity) (by positivity), Real.log_mul (by positivity) (by positivity),
      Real.log_mul (by norm_num) (by positivity)]
  have hMlo' : (P.μ : ℝ) * Real.log N - Real.log 2 ≤ Real.log (P.Mb N) := hMlo
  have hMhi' : Real.log (P.Mb N) ≤ (P.μ : ℝ) * Real.log N := hMhi
  have hM1' : (1 : ℝ) ≤ P.Mb N := hM1
  refine ⟨?_, ?_, ?_⟩
  · unfold Params.Xr; positivity
  · rw [hlog]
    have h := (div_le_iff₀ hμν).mp hlogN
    nlinarith
  · rw [hlog]
    have e : ((P.θ : ℝ) / (4 * a) + 1) * Lam (1 / 2) (N : ℝ) =
        (P.θ : ℝ) / (4 * a) * Lam (1 / 2) (N : ℝ) + Lam (1 / 2) (N : ℝ) := by ring
    rw [e]
    nlinarith

/-- `X_F = C′_d M H V_max^d`: `log X_F = (μ+h+dν) log N + O(Λ)`, eventually. -/
lemma ev_XF {P : Params} (hP : P.Admissible) {a : ℝ} (ha : 0 < a) : ∀ᶠ N : ℕ in atTop,
    0 < P.XF a N ∧
      ((P.μ : ℝ) + P.h + P.d * P.ν) / 2 * Real.log N ≤ Real.log (P.XF a N) ∧
      Real.log (P.XF a N) ≤ ((P.μ : ℝ) + P.h + P.d * P.ν) * Real.log N +
        (P.d * ((P.θ : ℝ) / (4 * a)) + 1) * Lam (1 / 2) (N : ℝ) := by
  obtain ⟨hν, -, hμν, hμ, hh, hθ, -, hd⟩ := Params.adm_real hP
  have hl2 : 0 < Real.log 2 := Real.log_pos (by norm_num)
  set Cd : ℝ := (2 : ℝ) ^ (2 * P.d + 5) + 2 * P.d with hCd
  have hCd1 : 1 ≤ Cd := by
    have : (1 : ℝ) ≤ (2 : ℝ) ^ (2 * P.d + 5) := one_le_pow₀ (by norm_num)
    have : (0 : ℝ) ≤ 2 * P.d := by positivity
    linarith
  have hcF : 0 < (P.μ : ℝ) + P.h + P.d * P.ν := by positivity
  filter_upwards [ev_Vmax hP ha, ev_nat_ge 1,
    ev_Lam_ge (Real.log Cd + P.d * Real.log 2),
    ev_log_ge (2 * ((2 + 3 * P.d) * Real.log 2) / ((P.μ : ℝ) + P.h + P.d * P.ν))] with
    N hVs hN1 hΛ hlogN
  obtain ⟨hV1, hVlo, hVhi⟩ := hVs
  have hN1' : 1 ≤ N := by exact_mod_cast hN1
  obtain ⟨hM1, hMlo, hMhi⟩ := log_floor_rpow hN1' hμ.le
  obtain ⟨hH1, hHlo, hHhi⟩ := log_floor_rpow hN1' hh.le
  have hΛ0 := Lam_half_nonneg (N : ℝ)
  have hθΛ : 0 ≤ (P.θ : ℝ) / (4 * a) * Lam (1 / 2) (N : ℝ) := by positivity
  have hM1' : (1 : ℝ) ≤ P.Mb N := hM1
  have hH1' : (1 : ℝ) ≤ P.Hb N := hH1
  have hlog : Real.log (P.XF a N) =
      Real.log Cd + Real.log (P.Mb N) + Real.log (P.Hb N) + P.d * Real.log (P.Vmax a N) := by
    unfold Params.XF
    rw [Real.log_mul (by positivity) (by positivity), Real.log_mul (by positivity) (by positivity),
      Real.log_mul (by positivity) (by positivity), Real.log_pow]
  have hMlo' : (P.μ : ℝ) * Real.log N - Real.log 2 ≤ Real.log (P.Mb N) := hMlo
  have hMhi' : Real.log (P.Mb N) ≤ (P.μ : ℝ) * Real.log N := hMhi
  have hHlo' : (P.h : ℝ) * Real.log N - Real.log 2 ≤ Real.log (P.Hb N) := hHlo
  have hHhi' : Real.log (P.Hb N) ≤ (P.h : ℝ) * Real.log N := hHhi
  have hlC : 0 ≤ Real.log Cd := Real.log_nonneg hCd1
  have hd0 : (0 : ℝ) ≤ P.d := by linarith
  have hVlo' := mul_le_mul_of_nonneg_left hVlo hd0
  have hVhi' := mul_le_mul_of_nonneg_left hVhi hd0
  refine ⟨?_, ?_, ?_⟩
  · unfold Params.XF
    have : (0 : ℝ) < (P.Vmax a N : ℝ) := by linarith
    positivity
  · rw [hlog]
    have h := (div_le_iff₀ hcF).mp hlogN
    nlinarith
  · rw [hlog]
    nlinarith

end FactoringA
