import LeanA.Recursion
import LeanA.DefsNFS

/-!
# General `L`-notation helpers for the cost assembly (`MainTheorems.lean`)

* monotonicity of `Lam` and `LN` in `x`;
* "constants and polylogarithms are `L^{o(1)}`" along `ℕ`;
* `EvLB Q c f`: `f N ≤ L_N[1/2, c + o(1)]` eventually, for the `N` satisfying `Q`, with its
  calculus (sums, products, constants, division by `p ≥ L^{−o(1)}`);
* conversions to `LBound`.
-/

noncomputable section
open Classical

namespace FactoringA

open Real Filter

/-! ### Monotonicity in `x` -/

lemma mt_loglog_max_mono {x y : ℝ} (hx : 1 ≤ x) (hxy : x ≤ y) :
    max 1 (Real.log (Real.log x)) ≤ max 1 (Real.log (Real.log y)) := by
  have hl0 : 0 ≤ Real.log x := Real.log_nonneg hx
  have hlxy : Real.log x ≤ Real.log y := Real.log_le_log (by linarith) hxy
  rcases hl0.eq_or_lt with h | h
  · rw [← h, Real.log_zero]
    exact max_le (le_max_left _ _) (by linarith [le_max_left 1 (Real.log (Real.log y))])
  · exact max_le_max le_rfl (Real.log_le_log h hlxy)

lemma mt_Lam_mono (α : ℝ) (hα0 : 0 ≤ α) (hα1 : α ≤ 1) {x y : ℝ} (hx : 1 ≤ x) (hxy : x ≤ y) :
    Lam α x ≤ Lam α y := by
  unfold Lam
  have hl0 : 0 ≤ Real.log x := Real.log_nonneg hx
  have hlxy : Real.log x ≤ Real.log y := Real.log_le_log (by linarith) hxy
  apply mul_le_mul (Real.rpow_le_rpow hl0 hlxy hα0)
    (Real.rpow_le_rpow (le_trans zero_le_one (le_max_left _ _)) (mt_loglog_max_mono hx hxy)
      (by linarith))
    (Real.rpow_nonneg (le_trans zero_le_one (le_max_left _ _)) _)
    (Real.rpow_nonneg (by linarith) _)

lemma mt_LN_mono_x (α c : ℝ) (hα0 : 0 ≤ α) (hα1 : α ≤ 1) (hc : 0 ≤ c) {x y : ℝ} (hx : 1 ≤ x)
    (hxy : x ≤ y) : LN α c x ≤ LN α c y :=
  Real.exp_le_exp.2 (mul_le_mul_of_nonneg_left (mt_Lam_mono α hα0 hα1 hx hxy) hc)

lemma mt_LN_pos (α c x : ℝ) : 0 < LN α c x := Real.exp_pos _

/-! ### Constants and polylogarithms along `ℕ` -/

lemma mt_ev_const_le (α : ℝ) (hα : 0 < α) (hα1 : α ≤ 1) (s : ℝ) (hs : 0 < s) (K : ℝ) :
    ∀ᶠ N : ℕ in atTop, K ≤ LN α s N := by
  filter_upwards [tendsto_natCast_atTop_atTop.eventually
    (polylog_le_LN α hα hα1 s hs K 0)] with N hN
  simpa using hN

lemma mt_ev_polylog_le (α : ℝ) (hα : 0 < α) (hα1 : α ≤ 1) (s : ℝ) (hs : 0 < s) (K k : ℝ) :
    ∀ᶠ N : ℕ in atTop, K * (Real.log N) ^ k ≤ LN α s N :=
  tendsto_natCast_atTop_atTop.eventually (polylog_le_LN α hα hα1 s hs K k)

lemma mt_LN_add (α a b x : ℝ) : LN α (a + b) x = LN α a x * LN α b x := (LN_mul α a b x).symm

/-! ### The `Λ`-scale -/

lemma mt_Lam_half_sq {x : ℝ} (hx : 0 ≤ Real.log x) :
    Lam (1 / 2) x ^ 2 = Real.log x * max 1 (Real.log (Real.log x)) := by
  unfold Lam
  have hm : 0 ≤ max 1 (Real.log (Real.log x)) := le_trans zero_le_one (le_max_left _ _)
  rw [mul_pow, ← Real.rpow_natCast, ← Real.rpow_natCast, ← Real.rpow_mul hx, ← Real.rpow_mul hm]
  norm_num

lemma mt_ev_log_le_mul (ε : ℝ) (hε : 0 < ε) (n : ℕ) :
    ∀ᶠ t : ℝ in atTop, Real.log t ^ n ≤ ε * t := by
  have h := (Real.isLittleO_pow_log_id_atTop (n := n)).bound hε
  filter_upwards [h, eventually_ge_atTop (1:ℝ)] with t ht ht1
  have h0 : 0 ≤ Real.log t := Real.log_nonneg ht1
  simp only [Real.norm_eq_abs, id] at ht
  rw [abs_of_nonneg (pow_nonneg h0 _), abs_of_nonneg (by linarith)] at ht
  exact ht

lemma mt_tendsto_log_nat : Tendsto (fun N : ℕ => Real.log (N : ℝ)) atTop atTop :=
  Real.tendsto_log_atTop.comp tendsto_natCast_atTop_atTop

/-- `Λ_N ≤ ε log N` eventually. -/
lemma mt_ev_Lam_le (ε : ℝ) (hε : 0 < ε) : ∀ᶠ N : ℕ in atTop,
    Lam (1 / 2) (N : ℝ) ≤ ε * Real.log N := by
  filter_upwards [mt_tendsto_log_nat.eventually (mt_ev_log_le_mul (ε ^ 2) (by positivity) 1),
    mt_tendsto_log_nat.eventually_ge_atTop 1,
    mt_tendsto_log_nat.eventually_ge_atTop (1 / ε ^ 2)] with N hl h1 hbig
  have h0 : 0 ≤ Real.log N := by linarith
  have hsq := mt_Lam_half_sq h0
  have hL : 0 ≤ Lam (1 / 2) (N : ℝ) := by
    unfold Lam
    exact mul_nonneg (Real.rpow_nonneg h0 _)
      (Real.rpow_nonneg (le_trans zero_le_one (le_max_left _ _)) _)
  have hm : max 1 (Real.log (Real.log N)) ≤ ε ^ 2 * Real.log N := by
    rw [pow_one] at hl
    refine max_le ?_ hl
    rw [div_le_iff₀ (by positivity)] at hbig
    linarith
  have : Lam (1 / 2) (N : ℝ) ^ 2 ≤ (ε * Real.log N) ^ 2 := by
    rw [hsq]; nlinarith
  exact (pow_le_pow_iff_left₀ hL (by positivity) two_ne_zero).1 this

/-! ### `EvLB`: eventual `L`-bounds along a predicate -/

/-- `f N ≤ L_N[1/2, c + o(1)]` eventually, for the `N` satisfying `Q`. -/
def EvLB (Q : ℕ → Prop) (c : ℝ) (f : ℕ → ℝ) : Prop :=
  ∀ ε > 0, ∀ᶠ N : ℕ in atTop, Q N → f N ≤ LN (1 / 2) (c + ε) N

namespace EvLB

variable {Q : ℕ → Prop} {c c' : ℝ} {f g p : ℕ → ℝ}

lemma mono (h : EvLB Q c f) (hc : c ≤ c') : EvLB Q c' f := by
  intro ε hε
  filter_upwards [h ε hε, eventually_ge_atTop 1] with N hN hN1 hQ
  exact (hN hQ).trans (LN_mono _ _ (by exact_mod_cast hN1) (by linarith))

lemma of_le (h : EvLB Q c g) (hfg : ∀ᶠ N : ℕ in atTop, Q N → f N ≤ g N) : EvLB Q c f := by
  intro ε hε
  filter_upwards [h ε hε, hfg] with N h1 h2 hQ
  exact (h2 hQ).trans (h1 hQ)

lemma of_LN (c : ℝ) : EvLB Q c (fun N => LN (1 / 2) c N) := by
  intro ε hε
  filter_upwards [eventually_ge_atTop 1] with N hN1 _
  exact LN_mono _ _ (by exact_mod_cast hN1) (by linarith)

lemma const (K : ℝ) : EvLB Q 0 (fun _ => K) := by
  intro ε hε
  filter_upwards [mt_ev_const_le (1 / 2) (by norm_num) (by norm_num) ε hε K] with N hN _
  simpa using hN

lemma polylog (K k : ℝ) : EvLB Q 0 (fun N => K * (Real.log N) ^ k) := by
  intro ε hε
  filter_upwards [mt_ev_polylog_le (1 / 2) (by norm_num) (by norm_num) ε hε K k] with N hN _
  simpa using hN

lemma add (hf : EvLB Q c f) (hg : EvLB Q c' g) : EvLB Q (max c c') (fun N => f N + g N) := by
  intro ε hε
  filter_upwards [hf (ε / 2) (by linarith), hg (ε / 2) (by linarith),
    mt_ev_const_le (1 / 2) (by norm_num) (by norm_num) (ε / 2) (by linarith) 2,
    eventually_ge_atTop 1] with N h1 h2 h3 hN1 hQ
  have hN1' : (1:ℝ) ≤ N := by exact_mod_cast hN1
  have a1 : LN (1 / 2) (c + ε / 2) N ≤ LN (1 / 2) (max c c' + ε / 2) N :=
    LN_mono _ _ hN1' (by linarith [le_max_left c c'])
  have a2 : LN (1 / 2) (c' + ε / 2) N ≤ LN (1 / 2) (max c c' + ε / 2) N :=
    LN_mono _ _ hN1' (by linarith [le_max_right c c'])
  have e : LN (1 / 2) (max c c' + ε) N =
      LN (1 / 2) (ε / 2) N * LN (1 / 2) (max c c' + ε / 2) N := by
    rw [LN_mul]; congr 1; ring
  rw [e]
  have hpos := mt_LN_pos (1 / 2) (max c c' + ε / 2) N
  have := mul_le_mul_of_nonneg_right h3 hpos.le
  linarith [h1 hQ, h2 hQ]

lemma mul (hf : EvLB Q c f) (hg : EvLB Q c' g) (hf0 : ∀ᶠ N : ℕ in atTop, Q N → 0 ≤ f N) :
    EvLB Q (c + c') (fun N => f N * g N) := by
  intro ε hε
  filter_upwards [hf (ε / 2) (by linarith), hg (ε / 2) (by linarith), hf0] with N h1 h2 h0 hQ
  calc f N * g N ≤ f N * LN (1 / 2) (c' + ε / 2) N := mul_le_mul_of_nonneg_left (h2 hQ) (h0 hQ)
    _ ≤ LN (1 / 2) (c + ε / 2) N * LN (1 / 2) (c' + ε / 2) N :=
        mul_le_mul_of_nonneg_right (h1 hQ) (mt_LN_pos _ _ _).le
    _ = LN (1 / 2) (c + c' + ε) N := by rw [LN_mul]; congr 1; ring

lemma const_mul {K : ℝ} (hK : 0 ≤ K) (hf : EvLB Q c f) : EvLB Q c (fun N => K * f N) :=
  ((const K).mul hf (Eventually.of_forall fun _ _ => hK)).mono (by simp)

lemma div (hf : EvLB Q c f) (hp : ∀ ε > 0, ∀ᶠ N : ℕ in atTop, Q N → LN (1 / 2) (-ε) N ≤ p N) :
    EvLB Q c (fun N => f N / p N) := by
  intro ε hε
  filter_upwards [hf (ε / 2) (by linarith), hp (ε / 2) (by linarith)] with N h1 h2 hQ
  have hp0 : 0 < p N := lt_of_lt_of_le (mt_LN_pos _ _ _) (h2 hQ)
  rw [div_le_iff₀ hp0]
  calc f N ≤ LN (1 / 2) (c + ε / 2) N := h1 hQ
    _ = LN (1 / 2) (c + ε) N * LN (1 / 2) (-(ε / 2)) N := by rw [LN_mul]; congr 1; ring
    _ ≤ LN (1 / 2) (c + ε) N * p N := mul_le_mul_of_nonneg_left (h2 hQ) (mt_LN_pos _ _ _).le

lemma of_lbound (h : LBound (1 / 2) c f) : EvLB Q c f := by
  intro ε hε
  filter_upwards [h ε hε] with N hN _ using hN

lemma toLBound (h : EvLB (fun _ => True) c f) : LBound (1 / 2) c f := by
  intro ε hε
  filter_upwards [h ε hε] with N hN using hN trivial

end EvLB

lemma mt_lbound_mono {α c c' : ℝ} {f : ℕ → ℝ} (h : LBound α c f) (hc : c ≤ c') :
    LBound α c' f := by
  intro ε hε
  filter_upwards [h ε hε, eventually_ge_atTop 1] with N hN hN1
  exact hN.trans (LN_mono _ _ (by exact_mod_cast hN1) (by linarith))

/-- From the conclusion of Lemma F (`F x ≤ K L_x^{c+η}` for `x ≥ 2`) to `LBound`. -/
lemma mt_lbound_of_K (α : ℝ) (hα : 0 < α) (hα1 : α ≤ 1) (F : ℝ → ℝ) (c : ℝ)
    (h : ∀ η > 0, ∃ K : ℝ, ∀ x, 2 ≤ x → F x ≤ K * LN α (c + η) x) :
    LBound α c (fun N => F N) := by
  intro ε hε
  obtain ⟨K, hK⟩ := h (ε / 2) (by linarith)
  filter_upwards [mt_ev_const_le α hα hα1 (ε / 2) (by linarith) K, eventually_ge_atTop 2]
    with N h1 h2
  have h2' : (2:ℝ) ≤ N := by exact_mod_cast h2
  calc F N ≤ K * LN α (c + ε / 2) N := hK N h2'
    _ ≤ LN α (ε / 2) N * LN α (c + ε / 2) N := mul_le_mul_of_nonneg_right h1 (mt_LN_pos _ _ _).le
    _ = LN α (c + ε) N := by rw [LN_mul]; congr 1; ring

end FactoringA
