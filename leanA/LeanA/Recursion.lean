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

/-! ## Helper lemmas -/

lemma Lam_nonneg (α x : ℝ) (hx : 1 ≤ x) : 0 ≤ Lam α x := by
  unfold Lam
  apply mul_nonneg
  · exact Real.rpow_nonneg (Real.log_nonneg hx) _
  · exact Real.rpow_nonneg (le_trans zero_le_one (le_max_left _ _)) _

lemma logpow_le_Lam (α x : ℝ) (hα1 : α ≤ 1) (hx : 1 ≤ x) : (Real.log x) ^ α ≤ Lam α x := by
  unfold Lam
  have h1 : 1 ≤ (max 1 (Real.log (Real.log x))) ^ (1 - α) :=
    Real.one_le_rpow (le_max_left _ _) (by linarith)
  have h0 : 0 ≤ (Real.log x) ^ α := Real.rpow_nonneg (Real.log_nonneg hx) _
  nlinarith

lemma LN_mono (α x : ℝ) (hx : 1 ≤ x) {b b' : ℝ} (h : b ≤ b') : LN α b x ≤ LN α b' x :=
  Real.exp_le_exp.2 (mul_le_mul_of_nonneg_right h (Lam_nonneg α x hx))

lemma one_le_LN (α x : ℝ) (hx : 1 ≤ x) {b : ℝ} (hb : 0 ≤ b) : 1 ≤ LN α b x :=
  Real.one_le_exp (mul_nonneg hb (Lam_nonneg α x hx))

lemma LN_mul (α a b x : ℝ) : LN α a x * LN α b x = LN α (a + b) x := by
  unfold LN; rw [← Real.exp_add, add_mul]

/-- Polylogarithms are `L^{o(1)}`. -/
lemma polylog_le_exp (α : ℝ) (hα : 0 < α) (s : ℝ) (hs : 0 < s) (C k : ℝ) :
    ∀ᶠ x in atTop, C * (Real.log x) ^ k ≤ Real.exp (s * (Real.log x) ^ α) := by
  have h1 : (fun w : ℝ => w ^ (k / α)) =o[atTop] fun w => Real.exp (s * w) :=
    isLittleO_rpow_exp_pos_mul_atTop (k / α) hs
  have h2 := h1.def (c := 1 / (|C| + 1)) (by positivity)
  have hlog : Tendsto (fun x => (Real.log x) ^ α) atTop atTop :=
    (tendsto_rpow_atTop hα).comp tendsto_log_atTop
  filter_upwards [hlog.eventually h2, tendsto_log_atTop.eventually (eventually_ge_atTop 0)]
    with x hx hl0
  have e : ((Real.log x) ^ α) ^ (k / α) = (Real.log x) ^ k := by
    rw [← Real.rpow_mul hl0]; congr 1; field_simp
  simp only [e, Real.norm_eq_abs, abs_of_nonneg (Real.rpow_nonneg hl0 _),
    abs_of_pos (Real.exp_pos _)] at hx
  have hp : 0 ≤ (Real.log x) ^ k := Real.rpow_nonneg hl0 _
  calc C * (Real.log x) ^ k ≤ |C| * (Real.log x) ^ k :=
        mul_le_mul_of_nonneg_right (le_abs_self C) hp
    _ ≤ (|C| + 1) * (Real.log x) ^ k := by nlinarith
    _ ≤ (|C| + 1) * (1 / (|C| + 1) * Real.exp (s * (Real.log x) ^ α)) :=
        mul_le_mul_of_nonneg_left hx (by positivity)
    _ = Real.exp (s * (Real.log x) ^ α) := by field_simp

lemma polylog_le_LN (α : ℝ) (hα : 0 < α) (hα1 : α ≤ 1) (s : ℝ) (hs : 0 < s) (C k : ℝ) :
    ∀ᶠ x in atTop, C * (Real.log x) ^ k ≤ LN α s x := by
  filter_upwards [polylog_le_exp α hα s hs C k, eventually_ge_atTop 1] with x hx hx1
  exact hx.trans (Real.exp_le_exp.2 (mul_le_mul_of_nonneg_left
    (logpow_le_Lam α x hα1 hx1) hs.le))

lemma two_le_sqrt_two_mul {x : ℝ} (hx : 2 ≤ x) : 2 ≤ Real.sqrt (2 * x) := by
  rw [Real.le_sqrt (by norm_num) (by linarith)]; linarith

/-- The scaling step: `Lam α (√(2x)) ≤ (2^{-α} + t) Lam α x` eventually. -/
lemma Lam_sqrt_le (α : ℝ) (hα0 : 0 < α) (hα1 : α < 1) (t : ℝ) (ht : 0 < t) :
    ∀ᶠ x in atTop, Lam α (Real.sqrt (2 * x)) ≤ ((2:ℝ) ^ (-α) + t) * Lam α x := by
  have hlog : Tendsto (fun x => (Real.log x) ^ α) atTop atTop :=
    (tendsto_rpow_atTop hα0).comp tendsto_log_atTop
  filter_upwards [hlog.eventually (eventually_ge_atTop ((2:ℝ) ^ (-α) * (Real.log 2) ^ α / t)),
    eventually_ge_atTop 2] with x hx hx2
  set r := (2:ℝ) ^ (-α) with hr
  set y := Real.sqrt (2 * x) with hy
  have hr0 : 0 < r := Real.rpow_pos_of_pos (by norm_num) _
  have hlx : 0 < Real.log x := Real.log_pos (by linarith)
  have hl2 : 0 < Real.log 2 := Real.log_pos (by norm_num)
  have hy2 : 2 ≤ y := two_le_sqrt_two_mul hx2
  have hyx : y ≤ x := by
    rw [hy, Real.sqrt_le_left (by linarith)]; nlinarith
  have hly : Real.log y = (Real.log 2 + Real.log x) / 2 := by
    rw [hy, Real.log_sqrt (by linarith), Real.log_mul (by norm_num) (by linarith)]
  have hlyp : 0 < Real.log y := Real.log_pos (by linarith)
  have hlyx : Real.log y ≤ Real.log x := Real.log_le_log (by linarith) hyx
  have hM : max 1 (Real.log (Real.log y)) ≤ max 1 (Real.log (Real.log x)) :=
    max_le_max le_rfl (Real.log_le_log hlyp hlyx)
  set M := max 1 (Real.log (Real.log x)) with hMdef
  have hM1 : 1 ≤ M := le_max_left _ _
  have hMp : (max 1 (Real.log (Real.log y))) ^ (1 - α) ≤ M ^ (1 - α) :=
    Real.rpow_le_rpow (le_trans zero_le_one (le_max_left _ _)) hM (by linarith)
  have hMp0 : 0 ≤ M ^ (1 - α) := Real.rpow_nonneg (by linarith) _
  have hpow : (Real.log y) ^ α ≤ r * ((Real.log 2) ^ α + (Real.log x) ^ α) := by
    rw [hly, Real.div_rpow (by linarith) (by norm_num), hr, Real.rpow_neg (by norm_num),
      div_eq_inv_mul]
    exact mul_le_mul_of_nonneg_left
      (Real.rpow_add_le_add_rpow hl2.le hlx.le hα0.le hα1.le) (by positivity)
  have hsmall : r * (Real.log 2) ^ α ≤ t * (Real.log x) ^ α := by
    rw [div_le_iff₀ ht] at hx; linarith
  have hlyα : 0 ≤ (Real.log y) ^ α := Real.rpow_nonneg hlyp.le _
  have hlxα : 0 ≤ (Real.log x) ^ α := Real.rpow_nonneg hlx.le _
  have hl2α : 0 ≤ (Real.log 2) ^ α := Real.rpow_nonneg hl2.le _
  unfold Lam
  calc (Real.log y) ^ α * (max 1 (Real.log (Real.log y))) ^ (1 - α)
      ≤ (r * ((Real.log 2) ^ α + (Real.log x) ^ α)) * M ^ (1 - α) :=
        mul_le_mul hpow hMp (Real.rpow_nonneg (le_trans zero_le_one (le_max_left _ _)) _)
          (by positivity)
    _ = r * (Real.log 2) ^ α * M ^ (1 - α) + r * (Real.log x) ^ α * M ^ (1 - α) := by ring
    _ ≤ t * (Real.log x) ^ α * M ^ (1 - α) + r * (Real.log x) ^ α * M ^ (1 - α) := by
        gcongr
    _ = (r + t) * ((Real.log x) ^ α * M ^ (1 - α)) := by ring

/-- Internal variant of Lemma F with the `η′`-slack form of the recursion hypothesis. -/
lemma lemmaF_core (α : ℝ) (hα0 : 0 < α) (hα1 : α < 1) (F : ℝ → ℝ)
    (hF0 : ∀ x, 2 ≤ x → 0 ≤ F x) (hmono : MonotoneOn F (Set.Ici 2))
    (p δ : ℝ) (hp : 0 ≤ p) (_hδ : 0 ≤ δ)
    (hrec : ∀ η > 0, ∀ᶠ x in atTop,
      F x ≤ LN α (p + η) x + LN α (δ + η) x * F (Real.sqrt (2 * x))) :
    ∀ η > 0, ∃ K : ℝ, ∀ x, 2 ≤ x → F x ≤ K * LN α (max p (δ / (1 - 2 ^ (-α))) + η) x := by
  intro η hη
  set r : ℝ := (2:ℝ) ^ (-α) with hr
  have hr0 : 0 < r := Real.rpow_pos_of_pos (by norm_num) _
  have hr1 : r < 1 := Real.rpow_lt_one_of_one_lt_of_neg (by norm_num) (by linarith)
  set c := max p (δ / (1 - r)) with hc
  have hc0 : 0 ≤ c := le_trans hp (le_max_left _ _)
  have hpc : p ≤ c := le_max_left _ _
  have hδc : δ ≤ c - c * r := by
    have h1 : δ / (1 - r) ≤ c := le_max_right _ _
    rw [div_le_iff₀ (by linarith)] at h1
    linarith
  set β := c + η with hβ
  have hβ0 : 0 < β := by linarith
  set g := η * (1 - r) / 4 with hg
  have hg0 : 0 < g := by have : 0 < 1 - r := by linarith
                         positivity
  set t := g / β with ht
  have ht0 : 0 < t := div_pos hg0 hβ0
  have hβt : β * t = g := by rw [ht]; field_simp
  have hev := (hrec g hg0).and ((Lam_sqrt_le α hα0 hα1 t ht0).and
    ((polylog_le_exp α hα0 (2 * g) (by positivity) 2 0).and (eventually_ge_atTop (4:ℝ))))
  obtain ⟨X, hX⟩ := eventually_atTop.1 hev
  set X' := max X 4 with hX'
  have hX'4 : 4 ≤ X' := le_max_right _ _
  set K := max 1 (F X') with hK
  have hK1 : 1 ≤ K := le_max_left _ _
  have claim : ∀ n : ℕ, ∀ x, 2 ≤ x → x ≤ X' + n → F x ≤ K * LN α β x := by
    intro n
    induction n with
    | zero =>
      intro x hx2 hxX
      simp only [CharP.cast_eq_zero, add_zero] at hxX
      have h1 : F x ≤ F X' := hmono (Set.mem_Ici.2 hx2) (Set.mem_Ici.2 (by linarith)) hxX
      have h2 : F X' ≤ K := le_max_right _ _
      have h3 : K ≤ K * LN α β x :=
        le_mul_of_one_le_right (by linarith) (one_le_LN α x (by linarith) hβ0.le)
      linarith
    | succ n ih =>
      intro x hx2 hxn
      by_cases hle : x ≤ X' + n
      · exact ih x hx2 hle
      push Not at hle
      push_cast at hxn
      obtain ⟨hrecx, hscale, hpoly, -⟩ := hX x (by
        have : X ≤ X' := le_max_left _ _
        have : (0:ℝ) ≤ n := Nat.cast_nonneg n
        linarith)
      set y := Real.sqrt (2 * x) with hy
      have hy2 : 2 ≤ y := two_le_sqrt_two_mul hx2
      have hyle : y ≤ X' + n := by
        have : (0:ℝ) ≤ n := Nat.cast_nonneg n
        have hZ : (4:ℝ) ≤ X' + n := by linarith
        have hZ2 : 4 * (X' + n) ≤ (X' + n) ^ 2 := by
          rw [pow_two]; exact mul_le_mul_of_nonneg_right hZ (by linarith)
        rw [hy, Real.sqrt_le_left (by linarith)]
        linarith
      have hFy := ih y hy2 hyle
      have hFy0 := hF0 y hy2
      set A := Lam α x with hA
      have hA0 : 0 ≤ A := Lam_nonneg α x (by linarith)
      have hlogA : (Real.log x) ^ α ≤ A := logpow_le_Lam α x hα1.le (by linarith)
      have h2 : (2:ℝ) ≤ Real.exp (2 * g * A) := by
        rw [Real.rpow_zero, mul_one] at hpoly
        exact hpoly.trans (Real.exp_le_exp.2 (mul_le_mul_of_nonneg_left hlogA (by positivity)))
      -- the second term
      have hLy : LN α β y ≤ Real.exp (β * (r + t) * A) := by
        unfold LN
        apply Real.exp_le_exp.2
        rw [mul_assoc]
        exact mul_le_mul_of_nonneg_left hscale hβ0.le
      have hexp1 : (δ + g) * A + β * (r + t) * A ≤ (β - 2 * g) * A := by
        rw [← add_mul]
        apply mul_le_mul_of_nonneg_right _ hA0
        have : β * (r + t) = c * r + η * r + g := by rw [mul_add, hβt, hβ]; ring
        rw [this]
        have : η - η * r = 4 * g := by rw [hg]; ring
        linarith
      have hT2 : LN α (δ + g) x * F y ≤ K * Real.exp ((β - 2 * g) * A) := by
        calc LN α (δ + g) x * F y ≤ LN α (δ + g) x * (K * Real.exp (β * (r + t) * A)) :=
              mul_le_mul_of_nonneg_left (hFy.trans
                (mul_le_mul_of_nonneg_left hLy (by linarith))) (Real.exp_pos _).le
          _ = K * Real.exp ((δ + g) * A + β * (r + t) * A) := by
              unfold LN; rw [Real.exp_add]; ring
          _ ≤ K * Real.exp ((β - 2 * g) * A) :=
              mul_le_mul_of_nonneg_left (Real.exp_le_exp.2 hexp1) (by linarith)
      have hT1 : LN α (p + g) x ≤ K * Real.exp ((β - 2 * g) * A) := by
        have : LN α (p + g) x ≤ Real.exp ((β - 2 * g) * A) := by
          unfold LN
          apply Real.exp_le_exp.2
          apply mul_le_mul_of_nonneg_right _ hA0
          have : η - η * r = 4 * g := by rw [hg]; ring
          have : 0 ≤ η * r := by positivity
          linarith
        exact this.trans (le_mul_of_one_le_left (Real.exp_pos _).le hK1)
      have hfin : 2 * (K * Real.exp ((β - 2 * g) * A)) ≤ K * LN α β x := by
        have e : LN α β x = Real.exp (2 * g * A) * Real.exp ((β - 2 * g) * A) := by
          unfold LN; rw [← hA, ← Real.exp_add]; ring_nf
        rw [e]
        have : 0 ≤ K * Real.exp ((β - 2 * g) * A) := by
          have := Real.exp_pos ((β - 2 * g) * A)
          positivity
        calc 2 * (K * Real.exp ((β - 2 * g) * A))
            ≤ Real.exp (2 * g * A) * (K * Real.exp ((β - 2 * g) * A)) :=
              mul_le_mul_of_nonneg_right h2 this
          _ = K * (Real.exp (2 * g * A) * Real.exp ((β - 2 * g) * A)) := by ring
      linarith
  refine ⟨K, fun x hx => ?_⟩
  have : x ≤ X' + (⌈x⌉₊ : ℝ) := by
    have := Nat.le_ceil x
    linarith
  exact claim ⌈x⌉₊ x hx this

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
  refine lemmaF_core α hα0 hα1 F hF0 hmono p δ hp hδ ?_
  intro η hη
  filter_upwards [hε.eventually (Iio_mem_nhds hη), eventually_ge_atTop x₀,
    eventually_ge_atTop (2:ℝ)] with x hεx hx₀ hx2
  have hεx' : ε x ≤ η := le_of_lt hεx
  have h1 := hrec x hx₀
  have hFy := hF0 _ (two_le_sqrt_two_mul hx2)
  have h2 : LN α (p + ε x) x ≤ LN α (p + η) x := LN_mono α x (by linarith) (by linarith)
  have h3 : LN α (δ + ε x) x ≤ LN α (δ + η) x := LN_mono α x (by linarith) (by linarith)
  have h4 := mul_le_mul_of_nonneg_right h3 hFy
  linarith

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
  refine lemmaF_core α hα0 hα1 F hF0 hmono p δ hp hδ ?_
  intro η hη
  have hη2 : 0 < η / 2 := by linarith
  filter_upwards [hP _ hη2, hD _ hη2, hD0, hrec,
    polylog_le_LN α hα0 hα1.le (η / 2) hη2 (4 / Real.log 2) 1,
    eventually_ge_atTop (2:ℝ)] with x hPx hDx hD0x hrecx hpoly hx2
  have hl2 : 0 < Real.log 2 := Real.log_pos (by norm_num)
  have hlb : 0 ≤ Real.logb 2 x := Real.logb_nonneg (by norm_num) (by linarith)
  have hpoly' : 4 * Real.logb 2 x ≤ LN α (η / 2) x := by
    rw [Real.rpow_one] at hpoly
    rw [Real.logb]
    calc 4 * (Real.log x / Real.log 2) = 4 / Real.log 2 * Real.log x := by ring
      _ ≤ _ := hpoly
  have hFy := hF0 _ (two_le_sqrt_two_mul hx2)
  have hLP0 : 0 ≤ LN α (p + η / 2) x := (Real.exp_pos _).le
  have hLD0 : 0 ≤ LN α (δ + η / 2) x := (Real.exp_pos _).le
  have hLe0 : 0 ≤ LN α (η / 2) x := (Real.exp_pos _).le
  -- first term
  have hT1 : 2 * Real.logb 2 x * P x ≤ LN α (p + η) x := by
    calc 2 * Real.logb 2 x * P x ≤ 2 * Real.logb 2 x * LN α (p + η / 2) x :=
          mul_le_mul_of_nonneg_left hPx (by positivity)
      _ ≤ 4 * Real.logb 2 x * LN α (p + η / 2) x :=
          mul_le_mul_of_nonneg_right (by linarith) hLP0
      _ ≤ LN α (η / 2) x * LN α (p + η / 2) x := mul_le_mul_of_nonneg_right hpoly' hLP0
      _ = LN α (p + η) x := by rw [LN_mul]; ring_nf
  have hT2 : 2 * Real.logb 2 x * (2 * D x * F (Real.sqrt (2 * x))) ≤
      LN α (δ + η) x * F (Real.sqrt (2 * x)) := by
    have : 4 * Real.logb 2 x * D x ≤ LN α (δ + η) x := by
      calc 4 * Real.logb 2 x * D x ≤ LN α (η / 2) x * D x :=
            mul_le_mul_of_nonneg_right hpoly' hD0x
        _ ≤ LN α (η / 2) x * LN α (δ + η / 2) x := mul_le_mul_of_nonneg_left hDx hLe0
        _ = LN α (δ + η) x := by rw [LN_mul]; ring_nf
    calc 2 * Real.logb 2 x * (2 * D x * F (Real.sqrt (2 * x)))
        = (4 * Real.logb 2 x * D x) * F (Real.sqrt (2 * x)) := by ring
      _ ≤ _ := mul_le_mul_of_nonneg_right this hFy
  calc F x ≤ 2 * Real.logb 2 x * (P x + 2 * D x * F (Real.sqrt (2 * x))) := hrecx
    _ = 2 * Real.logb 2 x * P x + 2 * Real.logb 2 x * (2 * D x * F (Real.sqrt (2 * x))) := by
        ring
    _ ≤ _ := add_le_add hT1 hT2

/-- MAIN. [pp] Corollary cor:P, the series of the restart schedule:
`2^{s₀+1} + Σ_{m ≥ 1} 2^{s₀+m} 2^{−m(m+1)/2} ≤ 4·2^{s₀}`. -/
theorem restart_series (s₀ : ℕ) :
    (2 : ℝ) ^ (s₀ + 1) + ∑' m : ℕ, (2 : ℝ) ^ (s₀ + (m + 1)) * (2 : ℝ)⁻¹ ^ ((m + 1) * (m + 2) / 2) ≤
      4 * 2 ^ s₀ := by
  have hle : ∀ m : ℕ, (2 : ℝ) ^ (s₀ + (m + 1)) * (2 : ℝ)⁻¹ ^ ((m + 1) * (m + 2) / 2) ≤
      2 ^ s₀ * (2 : ℝ)⁻¹ ^ m := by
    intro m
    have hk : 2 * m + 1 ≤ (m + 1) * (m + 2) / 2 := by
      rw [Nat.le_div_iff_mul_le (by norm_num)]; nlinarith [Nat.le_mul_self m]
    have h1 : (2 : ℝ)⁻¹ ^ ((m + 1) * (m + 2) / 2) ≤ (2 : ℝ)⁻¹ ^ (2 * m + 1) :=
      pow_le_pow_of_le_one (by norm_num) (by norm_num) hk
    have h2 : (2 : ℝ) ^ (m + 1) * (2 : ℝ)⁻¹ ^ (m + 1) = 1 := by
      rw [← mul_pow]; norm_num
    calc (2 : ℝ) ^ (s₀ + (m + 1)) * (2 : ℝ)⁻¹ ^ ((m + 1) * (m + 2) / 2)
        ≤ (2 : ℝ) ^ (s₀ + (m + 1)) * (2 : ℝ)⁻¹ ^ (2 * m + 1) :=
          mul_le_mul_of_nonneg_left h1 (by positivity)
      _ = 2 ^ s₀ * (2 : ℝ)⁻¹ ^ m * ((2 : ℝ) ^ (m + 1) * (2 : ℝ)⁻¹ ^ (m + 1)) := by
          rw [show 2 * m + 1 = m + (m + 1) by ring, pow_add, pow_add]; ring
      _ = 2 ^ s₀ * (2 : ℝ)⁻¹ ^ m := by rw [h2, mul_one]
  have hsum : Summable (fun m : ℕ => (2 : ℝ) ^ s₀ * (2 : ℝ)⁻¹ ^ m) :=
    (summable_geometric_of_lt_one (by norm_num) (by norm_num)).mul_left _
  have hsum' : Summable (fun m : ℕ =>
      (2 : ℝ) ^ (s₀ + (m + 1)) * (2 : ℝ)⁻¹ ^ ((m + 1) * (m + 2) / 2)) :=
    Summable.of_nonneg_of_le (fun m => by positivity) hle hsum
  have htsum := hsum'.tsum_le_tsum hle hsum
  rw [tsum_mul_left, tsum_geometric_inv_two] at htsum
  rw [pow_succ]
  linarith

/-- MAIN. [pp] Corollary cor:P, assembly: if `T_j(N) ≤ L_N[1/2, c_j + o(1)]` for each `j`, with
`c_j → 0`, and for each `j` the portfolio's expected time is eventually at most
`K_j (log N)^{K_j} T_j(N)`, then the portfolio runs in `L_N[1/2, o(1)]`. -/
theorem portfolio_assembly (T : ℕ → ℝ) (Tj : ℕ → ℕ → ℝ) (c : ℕ → ℝ)
    (hc : Tendsto c atTop (𝓝 0)) (hTj : ∀ j, LBound (1 / 2) (c j) (Tj j))
    (hport : ∀ j, ∃ K : ℝ, 0 ≤ K ∧ ∀ᶠ N : ℕ in atTop, T N ≤ K * (Real.log N) ^ K * Tj j N) :
    LBound (1 / 2) 0 T := by
  intro ε hε
  have hε3 : 0 < ε / 3 := by linarith
  obtain ⟨j, hj⟩ := (hc.eventually (Iio_mem_nhds hε3)).exists
  obtain ⟨K, hK0, hK⟩ := hport j
  have hpoly : ∀ᶠ N : ℕ in atTop, K * (Real.log N) ^ K ≤ LN (1 / 2) (ε / 3) N :=
    tendsto_natCast_atTop_atTop.eventually
      (polylog_le_LN (1 / 2) (by norm_num) (by norm_num) (ε / 3) hε3 K K)
  filter_upwards [hK, hTj j (ε / 3) hε3, hpoly, eventually_ge_atTop 1] with N hTN hTjN hpN hN1
  have hN1' : (1 : ℝ) ≤ N := by exact_mod_cast hN1
  have hA0 : 0 ≤ K * (Real.log N) ^ K :=
    mul_nonneg hK0 (Real.rpow_nonneg (Real.log_nonneg hN1') _)
  have hL2 : 0 ≤ LN (1 / 2) (c j + ε / 3) N := (Real.exp_pos _).le
  have hprod : K * (Real.log N) ^ K * Tj j N ≤
      LN (1 / 2) (ε / 3) N * LN (1 / 2) (c j + ε / 3) N := by
    rcases le_or_gt 0 (Tj j N) with h | h
    · calc K * (Real.log N) ^ K * Tj j N ≤ LN (1 / 2) (ε / 3) N * Tj j N :=
            mul_le_mul_of_nonneg_right hpN h
        _ ≤ _ := mul_le_mul_of_nonneg_left hTjN (Real.exp_pos _).le
    · have : K * (Real.log N) ^ K * Tj j N ≤ 0 := mul_nonpos_of_nonneg_of_nonpos hA0 h.le
      have : 0 ≤ LN (1 / 2) (ε / 3) N * LN (1 / 2) (c j + ε / 3) N := by
        have := Real.exp_pos ((ε / 3) * Lam (1 / 2) N)
        unfold LN; positivity
      linarith
  rw [LN_mul] at hprod
  have hmono := LN_mono (1 / 2) (N : ℝ) hN1' (b := ε / 3 + (c j + ε / 3)) (b' := 0 + ε)
    (by linarith)
  linarith

end FactoringA
