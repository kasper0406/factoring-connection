import LeanA.Defs

/-!
# (e) Lemma Kern-U and Lemma Kern; field-law algebra (Lemmas GB, homogeneity, forced residue,
pairs)

The main statements (MAIN) are copied into `Challenge.lean`; do not change them.
-/

noncomputable section
open Classical

namespace FactoringA

open Finset Real

/-! ### Lemma GB ([pp] Lemma lem:GB) -/

lemma zmod_natAbs_self (m : ℤ) : ((m : ZMod m.natAbs)) = 0 :=
  (ZMod.intCast_zmod_eq_zero_iff_dvd m m.natAbs).mpr (Int.natAbs_dvd.mpr dvd_rfl)

lemma isUnit_zmod_of_isCoprime {x m : ℤ} (h : IsCoprime x m) : IsUnit (x : ZMod m.natAbs) := by
  obtain ⟨u, v, huv⟩ := h
  have h2 : ((u : ZMod m.natAbs)) * (x : ZMod m.natAbs) = 1 := by
    have := congrArg (fun z : ℤ => (z : ZMod m.natAbs)) huv
    simpa [zmod_natAbs_self] using this
  exact IsUnit.of_mul_eq_one_right _ h2

lemma gb_dvd (g0 g1 : ℤ) (i : ℕ) (T : ℤ) (hcop : IsCoprime g0 g1) :
    g1 ∣ T - gbDigit g0 g1 i T * g0 ^ i := by
  rw [← Int.natAbs_dvd, ← ZMod.intCast_zmod_eq_zero_iff_dvd]
  have hu := isUnit_zmod_of_isCoprime hcop
  simp only [gbDigit, Int.cast_sub, Int.cast_mul, Int.cast_pow, ZMod.coe_valMinAbs]
  rw [mul_assoc, ← mul_pow, ZMod.inv_mul_of_unit _ hu, one_pow, mul_one, sub_self]

lemma gbNext_mul (g0 g1 : ℤ) (i : ℕ) (T : ℤ) (hcop : IsCoprime g0 g1) :
    gbNext g0 g1 i T * g1 = T - gbDigit g0 g1 i T * g0 ^ i :=
  Int.ediv_mul_cancel (gb_dvd g0 g1 i T hcop)

lemma gbT_succ (g0 g1 T : ℤ) (d n : ℕ) :
    gbT g0 g1 T d (n + 1) = gbNext g0 g1 (d - n) (gbT g0 g1 T d n) := rfl

lemma gb_partial (g0 g1 T : ℤ) (d : ℕ) (hcop : IsCoprime g0 g1) :
    ∀ n ≤ d, T = gbT g0 g1 T d n * g1 ^ n +
      ∑ j ∈ range n, gbCoeff g0 g1 T d (d - j) * g0 ^ (d - j) * g1 ^ j := by
  intro n
  induction n with
  | zero => intro _; simp [gbT]
  | succ n ih =>
    intro hn
    have h1 := ih (by omega)
    have hc : gbCoeff g0 g1 T d (d - n) = gbDigit g0 g1 (d - n) (gbT g0 g1 T d n) := by
      rw [gbCoeff, ite_eq_right (by omega), show d - (d - n) = n by omega]
    rw [sum_range_succ, hc, gbT_succ]
    have h2 := gbNext_mul g0 g1 (d - n) (gbT g0 g1 T d n) hcop
    linear_combination h1 - g1 ^ n * h2

lemma gbCoeff_zero (g0 g1 T : ℤ) (d : ℕ) : gbCoeff g0 g1 T d 0 = gbT g0 g1 T d d := by
  simp [gbCoeff]

/-- MAIN. [pp] Lemma lem:GB: `F̂(g₀, g₁) = T`. -/
theorem gb_identity (g0 g1 T : ℤ) (d : ℕ) (hg1 : 2 ≤ g1) (hcop : IsCoprime g0 g1) :
    gbForm g0 g1 T d g0 g1 = T := by
  have h := gb_partial g0 g1 T d hcop d le_rfl
  have hF : gbForm g0 g1 T d g0 g1 =
      ∑ j ∈ range (d + 1), gbCoeff g0 g1 T d (d - j) * g0 ^ (d - j) * g1 ^ j := by
    rw [gbForm, ← sum_range_reflect]
    refine sum_congr rfl fun j hj => ?_
    rw [mem_range] at hj
    rw [show d + 1 - 1 - j = d - j by omega, show d - (d - j) = j by omega]
  rw [hF, sum_range_succ, Nat.sub_self, pow_zero, mul_one, gbCoeff_zero]
  linarith

/-- MAIN. [pp] Lemma lem:GB: `|a_i| ≤ g₁/2` for `1 ≤ i ≤ d`. -/
theorem gb_digit_bound (g0 g1 T : ℤ) (d i : ℕ) (hg1 : 2 ≤ g1) (hi1 : 1 ≤ i) (hid : i ≤ d) :
    2 * |gbCoeff g0 g1 T d i| ≤ g1 := by
  have hne : gbCoeff g0 g1 T d i = gbDigit g0 g1 i (gbT g0 g1 T d (d - i)) := by
    rw [gbCoeff, ite_eq_right (by omega)]
  rw [hne, gbDigit]
  have : NeZero g1.natAbs := ⟨Int.natAbs_ne_zero.mpr (by omega)⟩
  have h := ZMod.valMinAbs_mem_Ioc
    ((gbT g0 g1 T d (d - i) : ZMod g1.natAbs) * ((g0 : ZMod g1.natAbs)⁻¹) ^ i)
  have hn : ((g1.natAbs : ℕ) : ℤ) = g1 := by rw [Int.natCast_natAbs, abs_of_pos (by omega)]
  rw [hn] at h
  obtain ⟨h1, h2⟩ := h
  generalize (((gbT g0 g1 T d (d - i) : ZMod g1.natAbs) *
    ((g0 : ZMod g1.natAbs)⁻¹) ^ i)).valMinAbs = x at h1 h2 ⊢
  rcases abs_cases x with ⟨hx, _⟩ | ⟨hx, _⟩ <;> rw [hx] <;> omega

/-- MAIN. [pp] Lemma lem:GB: if `g₀, g₁ ∈ [M/2, M]` then `|a₀| ≤ 2^d T/M^d + d 2^{d−1} M`. -/
theorem gb_a0_bound (g0 g1 T : ℤ) (d : ℕ) (M : ℝ) (hg1 : 2 ≤ g1) (hcop : IsCoprime g0 g1)
    (hT : 0 ≤ T) (h0 : M / 2 ≤ g0) (h0' : (g0 : ℝ) ≤ M) (h1 : M / 2 ≤ g1) (h1' : (g1 : ℝ) ≤ M) :
    |(gbCoeff g0 g1 T d 0 : ℝ)| ≤ 2 ^ d * T / M ^ d + d * 2 ^ (d - 1) * M := by
  have hg1r : (2:ℝ) ≤ g1 := by exact_mod_cast hg1
  have hg1pos : (0:ℝ) < g1 := by linarith
  have hM : 0 < M := by linarith
  have hg0pos : 0 < (g0:ℝ) := by linarith
  have hM2 : M ≤ 2 * g1 := by linarith
  have h := gb_partial g0 g1 T d hcop d le_rfl
  rw [← gbCoeff_zero] at h
  have hterm : ∀ j ∈ range d,
      |((gbCoeff g0 g1 T d (d - j) * g0 ^ (d - j) * g1 ^ j : ℤ) : ℝ)| ≤
        2 ^ (d - 1) * M * (g1:ℝ) ^ d := by
    intro j hj
    rw [mem_range] at hj
    obtain ⟨k, hk⟩ : ∃ k, d - j = k + 1 := ⟨d - j - 1, by omega⟩
    have hb := gb_digit_bound g0 g1 T d (d - j) hg1 (by omega) (by omega)
    have hbr : 2 * |(gbCoeff g0 g1 T d (d - j) : ℝ)| ≤ g1 := by exact_mod_cast hb
    rw [hk] at hbr ⊢
    generalize (gbCoeff g0 g1 T d (k + 1) : ℤ) = cc at hbr ⊢
    push_cast
    rw [abs_mul, abs_mul, abs_pow, abs_pow, abs_of_pos hg0pos, abs_of_pos hg1pos]
    have e1 : (g0:ℝ) ^ (k + 1) ≤ M * (2 * g1) ^ k := by
      rw [pow_succ']
      exact mul_le_mul h0' (pow_le_pow_left₀ hg0pos.le (by linarith) k) (by positivity) hM.le
    have e2 : (2:ℝ) ^ k ≤ 2 ^ (d - 1) := pow_le_pow_right₀ (by norm_num) (by omega)
    have e3 : (g1:ℝ) ^ d = g1 ^ k * g1 * g1 ^ j := by
      rw [← pow_succ, ← pow_add]; congr 1; omega
    have hcc : |(cc:ℝ)| ≤ g1 / 2 := by linarith
    calc |(cc:ℝ)| * (g0:ℝ) ^ (k + 1) * (g1:ℝ) ^ j
        ≤ (g1 / 2) * (M * (2 * g1) ^ k) * g1 ^ j := by
          gcongr
      _ = M * 2 ^ k * (g1 ^ k * g1 * g1 ^ j) / 2 := by rw [mul_pow]; ring
      _ ≤ M * 2 ^ (d - 1) * (g1 ^ k * g1 * g1 ^ j) / 2 := by gcongr
      _ ≤ 2 ^ (d - 1) * M * (g1:ℝ) ^ d := by
          rw [e3]
          have : 0 ≤ M * 2 ^ (d - 1) * ((g1:ℝ) ^ k * g1 * g1 ^ j) := by positivity
          linarith
  set S := ∑ j ∈ range d, gbCoeff g0 g1 T d (d - j) * g0 ^ (d - j) * g1 ^ j with hS
  have hSb : |(S : ℝ)| ≤ d * (2 ^ (d - 1) * M * (g1:ℝ) ^ d) := by
    rw [hS]; push_cast
    calc _ ≤ ∑ j ∈ range d, |((gbCoeff g0 g1 T d (d - j) : ℝ) * (g0:ℝ) ^ (d - j) * (g1:ℝ) ^ j)| :=
          abs_sum_le_sum_abs _ _
      _ ≤ ∑ j ∈ range d, (2 ^ (d - 1) * M * (g1:ℝ) ^ d) := by
          refine sum_le_sum fun j hj => ?_
          have := hterm j hj
          push_cast at this
          exact this
      _ = d * (2 ^ (d - 1) * M * (g1:ℝ) ^ d) := by rw [sum_const, card_range, nsmul_eq_mul]
  have hTr : (T:ℝ) = (gbCoeff g0 g1 T d 0 : ℝ) * (g1:ℝ) ^ d + S := by exact_mod_cast h
  have hgd : (0:ℝ) < (g1:ℝ) ^ d := by positivity
  have hTr0 : (0:ℝ) ≤ T := by exact_mod_cast hT
  have key : |(gbCoeff g0 g1 T d 0 : ℝ)| * (g1:ℝ) ^ d ≤ T + d * (2 ^ (d - 1) * M * (g1:ℝ) ^ d) := by
    have : (gbCoeff g0 g1 T d 0 : ℝ) * (g1:ℝ) ^ d = T - S := by linarith
    rw [← abs_of_pos hgd, ← abs_mul, abs_of_pos hgd, this]
    calc |(T:ℝ) - S| ≤ |(T:ℝ)| + |(S:ℝ)| := abs_sub _ _
      _ ≤ _ := by rw [abs_of_nonneg hTr0]; linarith
  have hT2 : (T:ℝ) ≤ 2 ^ d * T / M ^ d * (g1:ℝ) ^ d := by
    rw [div_mul_eq_mul_div, le_div_iff₀ (by positivity)]
    have : M ^ d ≤ (2 * g1) ^ d := pow_le_pow_left₀ hM.le hM2 d
    rw [mul_pow] at this
    have := mul_le_mul_of_nonneg_left this hTr0
    linarith
  have : |(gbCoeff g0 g1 T d 0 : ℝ)| * (g1:ℝ) ^ d ≤
      (2 ^ d * T / M ^ d + d * 2 ^ (d - 1) * M) * (g1:ℝ) ^ d := by
    nlinarith
  exact le_of_mul_le_mul_right this hgd

/-- MAIN. [pp] Lemma lem:GB: `a_d ≢ 0 (mod g₁)` when `T` is a unit modulo `g₁`. -/
theorem gb_ad_not_dvd (g0 g1 T : ℤ) (d : ℕ) (hd : 1 ≤ d) (hg1 : 2 ≤ g1)
    (hcop : IsCoprime g0 g1) (hT : IsCoprime T g1) :
    ¬ g1 ∣ gbCoeff g0 g1 T d d := by
  intro hdvd
  have hc : gbCoeff g0 g1 T d d = gbDigit g0 g1 d T := by
    rw [gbCoeff, ite_eq_right (by omega), Nat.sub_self]; rfl
  rw [hc, gbDigit, ← Int.natAbs_dvd, ← ZMod.intCast_zmod_eq_zero_iff_dvd,
    ZMod.coe_valMinAbs] at hdvd
  have hu0 := isUnit_zmod_of_isCoprime hcop
  have huT := isUnit_zmod_of_isCoprime hT
  have hinv : IsUnit ((g0 : ZMod g1.natAbs)⁻¹) :=
    IsUnit.of_mul_eq_one _ (ZMod.inv_mul_of_unit _ hu0)
  have h3 := huT.mul (hinv.pow d)
  rw [hdvd] at h3
  have : Nontrivial (ZMod g1.natAbs) := ZMod.nontrivial_iff.mpr (by omega)
  exact not_isUnit_zero h3

lemma gbForm_homog (g0 g1 T : ℤ) (d : ℕ) (a b : ℤ) (hg1 : 2 ≤ g1) (hcop : IsCoprime g0 g1) :
    g1 ^ d * gbForm g0 g1 T d a b ≡ T * b ^ d [ZMOD Gform g0 g1 a b] := by
  rw [Int.modEq_iff_dvd]
  have hT := gb_identity g0 g1 T d hg1 hcop
  have e : T * b ^ d - g1 ^ d * gbForm g0 g1 T d a b =
      ∑ i ∈ range (d + 1), gbCoeff g0 g1 T d i * (g1 * b) ^ (d - i) *
        ((g0 * b) ^ i - (g1 * a) ^ i) := by
    calc T * b ^ d - g1 ^ d * gbForm g0 g1 T d a b
        = gbForm g0 g1 T d g0 g1 * b ^ d - g1 ^ d * gbForm g0 g1 T d a b := by rw [hT]
      _ = _ := by
        unfold gbForm
        rw [sum_mul, mul_sum, ← sum_sub_distrib]
        refine sum_congr rfl fun i hi => ?_
        rw [mem_range] at hi
        obtain ⟨k, rfl⟩ : ∃ k, d = i + k := ⟨d - i, by omega⟩
        rw [Nat.add_sub_cancel_left]
        ring
  rw [e]
  refine dvd_sum fun i _ => dvd_mul_of_dvd_right ?_ _
  have := sub_dvd_pow_sub_pow (g1 * a) (g0 * b) i
  rw [← neg_sub]
  exact (dvd_neg.mpr this)

/-- MAIN. [pp] Lemma lem:homog: `g₁^d F_A(P) ≡ T b^d (mod G(P))`. -/
theorem homogeneity (g0 g1 T : ℤ) (d : ℕ) (c : Fin d → ℤ) (a b : ℤ) (hg1 : 2 ≤ g1)
    (hcop : IsCoprime g0 g1) :
    g1 ^ d * FA g0 g1 T d c a b ≡ T * b ^ d [ZMOD Gform g0 g1 a b] := by
  have h1 := gbForm_homog g0 g1 T d a b hg1 hcop
  refine Int.ModEq.trans ?_ h1
  rw [Int.modEq_iff_dvd]
  refine ⟨-(g1 ^ d * dotZ c (vP d a b)), ?_⟩
  unfold FA
  ring

/-- MAIN. [pp] Lemma lem:coprime (coprimality and the forced residue). Here `T = N w`,
`t = |G(P)| = ℓ|n|` with `n = G(P)/ℓ` `B`-smooth, `b` a prime with `B < b < ℓ`, `ℓ ∤ N` prime,
every prime factor of `N` exceeds `B`, and `w > t` is prime. Then `gcd(t, b g₁ N w) = 1` and
`F̂(P) ≡ N b^d g₁^{−d} w (mod t)`. -/
theorem forced_residue (N w B ℓ : ℕ) (g0 g1 a b : ℤ) (d : ℕ) (hg1 : 2 ≤ g1)
    (hcop : IsCoprime g0 g1) (hb : (b.natAbs).Prime) (hb0 : 0 < b) (hBb : (B : ℤ) < b)
    (hbℓ : b < ℓ) (hℓ : ℓ.Prime) (hℓN : ¬ ℓ ∣ N) (hdiv : (ℓ : ℤ) ∣ Gform g0 g1 a b)
    (hsm : IsSmooth B (Gform g0 g1 a b / ℓ)) (hNB : ∀ p : ℕ, p.Prime → p ∣ N → B < p)
    (hw : w.Prime) (hwt : |Gform g0 g1 a b| < w) :
    Nat.Coprime (Gform g0 g1 a b).natAbs (b.natAbs * g1.natAbs * N * w) ∧
      ((gbForm g0 g1 ((N : ℤ) * w) d a b : ℤ) : ZMod (Gform g0 g1 a b).natAbs) =
        (N : ZMod (Gform g0 g1 a b).natAbs) * (b : ZMod (Gform g0 g1 a b).natAbs) ^ d *
          ((g1 : ZMod (Gform g0 g1 a b).natAbs)⁻¹) ^ d * (w : ZMod (Gform g0 g1 a b).natAbs) := by
  set G := Gform g0 g1 a b with hGdef
  have hG : G = ℓ * (G / ℓ) := (Int.mul_ediv_cancel' hdiv).symm
  have hG0 : G ≠ 0 := by
    intro h0; apply hsm.1; rw [h0]; simp
  have ht0 : 0 < G.natAbs := Int.natAbs_pos.mpr hG0
  have key : ∀ q : ℕ, q.Prime → q ∣ G.natAbs → q = ℓ ∨ q ≤ B := by
    intro q hq hqt
    have hqG : (q:ℤ) ∣ G := Int.natCast_dvd.mpr hqt
    rw [hG] at hqG
    rcases (Nat.prime_iff_prime_int.mp hq).dvd_or_dvd hqG with h | h
    · left; exact (Nat.prime_dvd_prime_iff_eq hq hℓ).mp (Int.natCast_dvd_natCast.mp h)
    · right; exact hsm.2 q hq h
  have hbabs : ((b.natAbs : ℕ) : ℤ) = b := Int.natAbs_of_nonneg hb0.le
  have cb : Nat.Coprime G.natAbs b.natAbs := by
    apply Nat.coprime_of_dvd
    intro q hq hqt hqb
    have hqe : q = b.natAbs := (Nat.prime_dvd_prime_iff_eq hq hb).mp hqb
    rcases key q hq hqt with h | h
    · have : (q:ℤ) = ℓ := by rw [h]
      rw [hqe, hbabs] at this; omega
    · have : (q:ℤ) ≤ B := by exact_mod_cast h
      rw [hqe, hbabs] at this; omega
  have cg : Nat.Coprime G.natAbs g1.natAbs := by
    apply Nat.coprime_of_dvd
    intro q hq hqt hqg
    have hqpi := Nat.prime_iff_prime_int.mp hq
    have hqg' : (q:ℤ) ∣ g1 := Int.natCast_dvd.mpr hqg
    have hqG : (q:ℤ) ∣ G := Int.natCast_dvd.mpr hqt
    have hq0b : (q:ℤ) ∣ g0 * b := by
      have : g0 * b = g1 * a - G := by rw [hGdef, Gform]; ring
      rw [this]; exact dvd_sub (dvd_mul_of_dvd_left hqg' _) hqG
    rcases hqpi.dvd_or_dvd hq0b with h | h
    · exact hqpi.not_isUnit (hcop.isUnit_of_dvd' h hqg')
    · have hqb : q ∣ b.natAbs := Int.natCast_dvd.mp h
      have h2 := Nat.dvd_gcd hqt hqb
      rw [cb.gcd_eq_one] at h2
      exact hq.one_lt.ne' (Nat.dvd_one.mp h2)
  have cN : Nat.Coprime G.natAbs N := by
    apply Nat.coprime_of_dvd
    intro q hq hqt hqN
    rcases key q hq hqt with h | h
    · exact hℓN (h ▸ hqN)
    · have := hNB q hq hqN; omega
  have cw : Nat.Coprime G.natAbs w := by
    apply Nat.coprime_of_dvd
    intro q hq hqt hqw
    have hqe : q = w := (Nat.prime_dvd_prime_iff_eq hq hw).mp hqw
    have h1 : q ≤ G.natAbs := Nat.le_of_dvd ht0 hqt
    have h2 : ((G.natAbs : ℕ) : ℤ) < w := by rw [Int.natCast_natAbs]; exact hwt
    have h3 : G.natAbs < w := by exact_mod_cast h2
    omega
  refine ⟨Nat.Coprime.mul_right (Nat.Coprime.mul_right (Nat.Coprime.mul_right cb cg) cN) cw, ?_⟩
  have hc := gbForm_homog g0 g1 ((N : ℤ) * w) d a b hg1 hcop
  rw [← hGdef] at hc
  have hz : ((g1 ^ d * gbForm g0 g1 ((N : ℤ) * w) d a b : ℤ) : ZMod G.natAbs) =
      (((N : ℤ) * w * b ^ d : ℤ) : ZMod G.natAbs) := by
    rw [ZMod.intCast_eq_intCast_iff_dvd_sub, Int.natAbs_dvd]
    exact Int.modEq_iff_dvd.mp hc
  have hcopg : IsCoprime g1 G := by
    rw [Int.isCoprime_iff_gcd_eq_one]
    exact Nat.coprime_comm.mp cg
  have hu := isUnit_zmod_of_isCoprime hcopg
  push_cast at hz
  calc ((gbForm g0 g1 ((N : ℤ) * w) d a b : ℤ) : ZMod G.natAbs)
      = ((g1 : ZMod G.natAbs)⁻¹ * (g1 : ZMod G.natAbs)) ^ d *
          ((gbForm g0 g1 ((N : ℤ) * w) d a b : ℤ) : ZMod G.natAbs) := by
        rw [ZMod.inv_mul_of_unit _ hu, one_pow, one_mul]
    _ = ((g1 : ZMod G.natAbs)⁻¹) ^ d * ((g1 : ZMod G.natAbs) ^ d *
          ((gbForm g0 g1 ((N : ℤ) * w) d a b : ℤ) : ZMod G.natAbs)) := by ring
    _ = _ := by rw [hz]; ring

/-! ### Lemma Kern ([pp] Lemma lem:kern) and Lemma Kern-U (lv5unif §3) -/

section KernPoly
open Polynomial

/-- The polynomial `Σ c_i Y^i` attached to a coefficient vector. -/
def kpoly (d : ℕ) (c : Fin d → ℤ) : ℤ[X] := ∑ i : Fin d, C (c i) * X ^ (i : ℕ)

lemma kpoly_coeff (d : ℕ) (c : Fin d → ℤ) (n : ℕ) :
    (kpoly d c).coeff n = if h : n < d then c ⟨n, h⟩ else 0 := by
  rw [kpoly, finsetSum_coeff]
  simp only [coeff_C_mul_X_pow]
  split_ifs with h
  · rw [Finset.sum_eq_single ⟨n, h⟩]
    · simp
    · intro j _ hj; rw [ite_eq_right]; intro hn; exact hj (Fin.ext hn.symm)
    · simp
  · apply Finset.sum_eq_zero; intro j _; rw [ite_eq_right]; intro hn; exact h (hn ▸ j.2)

lemma qmul_coeff (a b a' b' : ℤ) (s : ℤ[X]) (i : ℕ) :
    (((C b - C a * X) * (C b' - C a' * X)) * s).coeff i =
      b * b' * s.coeff i - (a * b' + a' * b) * (if 1 ≤ i then s.coeff (i - 1) else 0) +
        a * a' * (if 2 ≤ i then s.coeff (i - 2) else 0) := by
  have hq : ((C b - C a * X) * (C b' - C a' * X)) * s =
      C (b * b') * s - C (a * b' + a' * b) * (X ^ 1 * s) + C (a * a') * (X ^ 2 * s) := by
    simp only [map_mul, map_add]; ring
  rw [hq, coeff_add, coeff_sub, coeff_C_mul, coeff_C_mul, coeff_C_mul, coeff_X_pow_mul',
    coeff_X_pow_mul']

lemma dot_eval (d : ℕ) (c : Fin d → ℤ) (a b : ℤ) (ha : a ≠ 0) :
    ((dotZ c (vP d a b) : ℤ) : ℚ) =
      (a : ℚ) ^ (d - 1) * ((kpoly d c).map (Int.castRingHom ℚ)).eval ((b : ℚ) / a) := by
  simp only [dotZ, vP, kpoly, Polynomial.map_sum, Polynomial.map_mul, Polynomial.map_pow,
    map_X, eval_finsetSum, eval_mul, eval_pow, eval_X, Int.cast_sum, Int.cast_mul,
    Int.cast_pow, mul_sum, eq_intCast, Polynomial.map_intCast, eval_intCast]
  refine Finset.sum_congr rfl fun i _ => ?_
  have hi : (i:ℕ) ≤ d - 1 := by have := i.2; omega
  have haq : (a:ℚ) ≠ 0 := by exact_mod_cast ha
  have hp : (a:ℚ) ^ (d - 1) = a ^ (d - 1 - (i:ℕ)) * a ^ (i:ℕ) := by
    rw [← pow_add]; congr 1; omega
  have hai : (a:ℚ) ^ (i:ℕ) ≠ 0 := pow_ne_zero _ haq
  rw [hp, div_pow]
  field_simp
  try ring

lemma dot_zero_left (d : ℕ) (hd : 1 ≤ d) (c : Fin d → ℤ) (b : ℤ) :
    dotZ c (vP d 0 b) = c ⟨d - 1, by omega⟩ * b ^ (d - 1) := by
  rw [dotZ, Finset.sum_eq_single ⟨d - 1, by omega⟩]
  · simp [vP]
  · intro j _ hj
    have : d - 1 - (j:ℕ) ≠ 0 := by
      intro h; apply hj; ext; simp only; have := j.2; omega
    simp [vP, zero_pow this]
  · simp

lemma map_lin (a b : ℤ) :
    (C b - C a * X : ℤ[X]).map (Int.castRingHom ℚ) = C (b:ℚ) - C (a:ℚ) * X := by
  simp [Polynomial.map_sub, Polynomial.map_mul]

lemma lin_dvd (d : ℕ) (c : Fin d → ℤ) (a b : ℤ) (hab : IsCoprime a b)
    (h0 : dotZ c (vP d a b) = 0) :
    (C (b:ℚ) - C (a:ℚ) * X) ∣ (kpoly d c).map (Int.castRingHom ℚ) := by
  by_cases ha : a = 0
  · subst ha
    have hbu : IsUnit b := isCoprime_zero_left.mp hab
    have hbq : (b:ℚ) ≠ 0 := by
      have : b ≠ 0 := fun h => by rw [h] at hbu; exact not_isUnit_zero hbu
      exact_mod_cast this
    simp only [Int.cast_zero, map_zero, zero_mul, sub_zero]
    exact (isUnit_C.mpr hbq.isUnit).dvd
  · have haq : (a:ℚ) ≠ 0 := by exact_mod_cast ha
    have hroot : ((kpoly d c).map (Int.castRingHom ℚ)).IsRoot ((b:ℚ) / a) := by
      have := dot_eval d c a b ha
      rw [h0, Int.cast_zero] at this
      exact (mul_eq_zero.mp this.symm).resolve_left (pow_ne_zero _ haq)
    obtain ⟨g, hg⟩ := dvd_iff_isRoot.mpr hroot
    refine ⟨C (-(1 / (a:ℚ))) * g, ?_⟩
    rw [hg, ← mul_assoc]
    congr 1
    have h1 : (a:ℚ) * (-(1 / (a:ℚ))) = -1 := by field_simp
    have h2 : (b:ℚ) * (-(1 / (a:ℚ))) = -((b:ℚ) / a) := by ring
    calc X - C ((b:ℚ) / a) = C ((b:ℚ) * (-(1 / (a:ℚ)))) - C ((a:ℚ) * (-(1 / (a:ℚ)))) * X := by
          rw [h1, h2]; simp; ring
      _ = _ := by simp only [map_mul]; ring

lemma lin_coprime (a b a' b' : ℤ) (hΔ : a * b' - a' * b ≠ 0) :
    IsCoprime (C (b:ℚ) - C (a:ℚ) * X) (C (b':ℚ) - C (a':ℚ) * X) := by
  have hΔq : ((a:ℚ) * b' - a' * b) ≠ 0 := by exact_mod_cast hΔ
  set u : ℚ := (a':ℚ) / (-((a:ℚ) * b' - a' * b))
  set v : ℚ := -(a:ℚ) / (-((a:ℚ) * b' - a' * b))
  have h1 : u * b + v * b' = 1 := by
    simp only [u, v]; field_simp; ring
  have h2 : u * a + v * a' = 0 := by
    simp only [u, v]; field_simp; ring
  refine ⟨C u, C v, ?_⟩
  calc C u * (C (b:ℚ) - C (a:ℚ) * X) + C v * (C (b':ℚ) - C (a':ℚ) * X)
      = C (u * b + v * b') - C (u * a + v * a') * X := by simp only [map_add, map_mul]; ring
    _ = 1 := by rw [h1, h2]; simp

lemma lin_primitive (a b : ℤ) (hab : IsCoprime a b) : (C b - C a * X : ℤ[X]).IsPrimitive := by
  intro r hr
  rw [C_dvd_iff_dvd_coeff] at hr
  have h0 := hr 0
  have h1 := hr 1
  rw [coeff_sub, coeff_C_zero, coeff_C_mul, coeff_X_zero, mul_zero, sub_zero] at h0
  rw [coeff_sub, coeff_C_mul, coeff_X_one, mul_one, coeff_C, ite_eq_right one_ne_zero, zero_sub,
    dvd_neg] at h1
  exact hab.isUnit_of_dvd' h1 h0

lemma kern_poly_dot_zero (d : ℕ) (c : Fin d → ℤ) (a b : ℤ) (ha : a ≠ 0) (g : ℤ[X])
    (h : kpoly d c = (C b - C a * X) * g) : dotZ c (vP d a b) = 0 := by
  have haq : (a:ℚ) ≠ 0 := by exact_mod_cast ha
  have := dot_eval d c a b ha
  rw [h, Polynomial.map_mul, map_lin, eval_mul] at this
  have h0 : eval ((b:ℚ) / a) (C (b:ℚ) - C (a:ℚ) * X) = 0 := by
    simp; field_simp; ring
  rw [h0, zero_mul, mul_zero] at this
  exact_mod_cast this

lemma kern_kernel_aux (d : ℕ) (hd : 3 ≤ d) (a b a' b' : ℤ) (hb : 1 ≤ b) (hb' : 1 ≤ b')
    (hP : IsCoprime a b) (hP' : IsCoprime a' b') (hΔ : a * b' - a' * b ≠ 0) (c : Fin d → ℤ) :
    PhiK d a b a' b' c = (0, 0) ↔
      ∃ e : ℕ → ℤ, (∀ j, d - 2 ≤ j → e j = 0) ∧
        ∀ i : Fin d, c i = b * b' * e i - (a * b' + a' * b) * (if 1 ≤ (i : ℕ) then e (i - 1) else 0) +
          a * a' * (if 2 ≤ (i : ℕ) then e (i - 2) else 0) := by
  have hd1 : 1 ≤ d := by omega
  have hb0 : b ≠ 0 := by omega
  have hb0' : b' ≠ 0 := by omega
  constructor
  · intro hz
    simp only [PhiK, Prod.mk.injEq] at hz
    obtain ⟨hz1, hz2⟩ := hz
    have hdvd : ((C b - C a * X) * (C b' - C a' * X) : ℤ[X]) ∣ kpoly d c := by
      rw [IsPrimitive.Int.dvd_iff_map_cast_dvd_map_cast _ _
        ((lin_primitive a b hP).mul (lin_primitive a' b' hP')), Polynomial.map_mul, map_lin,
        map_lin]
      exact (lin_coprime a b a' b' hΔ).mul_dvd (lin_dvd d c a b hP hz1) (lin_dvd d c a' b' hP' hz2)
    obtain ⟨s, hs⟩ := hdvd
    have hcoef : ∀ n, (kpoly d c).coeff n =
        b * b' * s.coeff n - (a * b' + a' * b) * (if 1 ≤ n then s.coeff (n - 1) else 0) +
          a * a' * (if 2 ≤ n then s.coeff (n - 2) else 0) := fun n => by rw [hs, qmul_coeff]
    -- `c_{d-1} = 0` as soon as `a a' = 0`
    have hcd : a * a' = 0 → c ⟨d - 1, by omega⟩ = 0 := by
      intro haa
      rcases mul_eq_zero.mp haa with h | h
      · rw [h, dot_zero_left d hd1] at hz1
        exact (mul_eq_zero.mp hz1).resolve_right (pow_ne_zero _ hb0)
      · rw [h, dot_zero_left d hd1] at hz2
        exact (mul_eq_zero.mp hz2).resolve_right (pow_ne_zero _ hb0')
    refine ⟨s.coeff, ?_, ?_⟩
    · intro j hj
      by_contra hne
      have hs0 : s ≠ 0 := by rintro rfl; simp at hne
      have hjm : j ≤ s.natDegree := le_natDegree_of_ne_zero hne
      have hlead : s.coeff s.natDegree ≠ 0 := mt leadingCoeff_eq_zero.mp hs0
      have hz : ∀ k, s.natDegree < k → s.coeff k = 0 := fun k hk => coeff_eq_zero_of_natDegree_lt hk
      set m := s.natDegree with hm
      have h2 := hcoef (m + 2)
      rw [kpoly_coeff, dite_eq_right (by omega), hz (m + 2) (by omega), ite_eq_left (by omega),
        ite_eq_left (by omega), show m + 2 - 1 = m + 1 by omega, hz (m + 1) (by omega),
        show m + 2 - 2 = m by omega] at h2
      have haa : a * a' = 0 := by
        have : a * a' * s.coeff m = 0 := by linarith
        exact (mul_eq_zero.mp this).resolve_right hlead
      have hl : (kpoly d c).coeff (m + 1) = 0 := by
        rw [kpoly_coeff]
        split_ifs with hlt
        · have : (⟨m + 1, hlt⟩ : Fin d) = ⟨d - 1, by omega⟩ := by ext; simp only; omega
          rw [this]; exact hcd haa
        · rfl
      have h1 := hcoef (m + 1)
      rw [hl, hz (m + 1) (by omega), ite_eq_left (by omega), Nat.add_sub_cancel, haa] at h1
      have hu : a * b' + a' * b = 0 := by
        have : (a * b' + a' * b) * s.coeff m = 0 := by
          simp only [mul_zero, zero_mul, zero_sub, add_zero] at h1; linarith
        exact (mul_eq_zero.mp this).resolve_right hlead
      apply hΔ
      rcases mul_eq_zero.mp haa with h | h
      · subst h
        have : a' * b = 0 := by linarith
        have ha' : a' = 0 := (mul_eq_zero.mp this).resolve_right hb0
        subst ha'; ring
      · subst h
        have : a * b' = 0 := by linarith
        have ha : a = 0 := (mul_eq_zero.mp this).resolve_right hb0'
        subst ha; ring
    · intro i
      have := hcoef i
      rw [kpoly_coeff, dite_eq_left i.2] at this
      simpa using this
  · rintro ⟨e, he, hc⟩
    let s : ℤ[X] := ∑ j ∈ range (d - 2), C (e j) * X ^ j
    have hs : ∀ j, s.coeff j = e j := by
      intro j
      simp only [s, finsetSum_coeff, coeff_C_mul_X_pow]
      rw [Finset.sum_eq_single j]
      · simp
      · intro k _ hk; rw [ite_eq_right (Ne.symm hk)]
      · intro hj; rw [mem_range, not_lt] at hj; rw [he j hj]; simp
    have hpq : kpoly d c = ((C b - C a * X) * (C b' - C a' * X)) * s := by
      ext n
      rw [qmul_coeff, kpoly_coeff]
      simp only [hs]
      by_cases hn : n < d
      · rw [dite_eq_left hn, hc ⟨n, hn⟩]
      · rw [dite_eq_right hn, he n (by omega)]
        split_ifs <;> simp [he (n - 1) (by omega), he (n - 2) (by omega)]
    have hcd : a * a' = 0 → c ⟨d - 1, by omega⟩ = 0 := by
      intro haa
      rw [hc ⟨d - 1, by omega⟩]
      simp only
      rw [he (d - 1) (by omega), he (d - 1 - 1) (by omega), haa]
      simp
    simp only [PhiK, Prod.mk.injEq]
    constructor
    · by_cases ha : a = 0
      · rw [ha, dot_zero_left d hd1, hcd (by rw [ha]; ring), zero_mul]
      · exact kern_poly_dot_zero d c a b ha _ (by rw [hpq, mul_assoc])
    · by_cases ha : a' = 0
      · rw [ha, dot_zero_left d hd1, hcd (by rw [ha]; ring), zero_mul]
      · exact kern_poly_dot_zero d c a' b' ha ((C b - C a * X) * s)
          (by rw [hpq]; ring)

end KernPoly

lemma dotZ_sub {d : ℕ} (c c' v : Fin d → ℤ) : dotZ (c - c') v = dotZ c v - dotZ c' v := by
  simp [dotZ, sub_mul, Finset.sum_sub_distrib]

/-- MAIN. [pp] Lemma lem:kern(ii): for primitive `P = (a, b)`, `P′ = (a′, b′)` with `b, b′ ≥ 1`
and `Δ ≠ 0`, `ker Φ` consists exactly of the coefficient vectors of `Q·E` with
`Q = (bx − ay)(b′x − a′y)` and `E` a form of degree `d − 3`; in coordinates
`c_i = bb′ e_i − (ab′ + a′b) e_{i−1} + aa′ e_{i−2}` (with `e_j = 0` outside `[0, d−3]`). -/
theorem kern_kernel (d : ℕ) (hd : 3 ≤ d) (a b a' b' : ℤ) (hb : 1 ≤ b) (hb' : 1 ≤ b')
    (hP : IsCoprime a b) (hP' : IsCoprime a' b') (hΔ : a * b' - a' * b ≠ 0) (c : Fin d → ℤ) :
    PhiK d a b a' b' c = (0, 0) ↔
      ∃ e : ℕ → ℤ, (∀ j, d - 2 ≤ j → e j = 0) ∧
        ∀ i : Fin d, c i = b * b' * e i - (a * b' + a' * b) * (if 1 ≤ (i : ℕ) then e (i - 1) else 0) +
          a * a' * (if 2 ≤ (i : ℕ) then e (i - 2) else 0) := by
  exact kern_kernel_aux d hd a b a' b' hb hb' hP hP' hΔ c

/-- MAIN. Lemma Kern-U (lv5unif §3), count: every fibre of `Φ` on `C_H` has at most
`(⌊2H/(bb′)⌋ + 1)^{d−2}` points. -/
theorem kernU_count (d : ℕ) (hd : 3 ≤ d) (a b a' b' : ℤ) (hb : 1 ≤ b) (hb' : 1 ≤ b')
    (hP : IsCoprime a b) (hP' : IsCoprime a' b') (hΔ : a * b' - a' * b ≠ 0) (H : ℕ)
    (z : ℤ × ℤ) :
    (((cube d H).filter (fun c => PhiK d a b a' b' c = z)).card : ℤ) ≤
      ((2 * (H : ℤ)) / (b * b') + 1) ^ (d - 2) := by
  have hβpos : 0 < b * b' := mul_pos (by omega) (by omega)
  have hβne : b * b' ≠ 0 := hβpos.ne'
  set K := 2 * (H : ℤ) / (b * b') with hKdef
  have hK : 0 ≤ K := Int.ediv_nonneg (by positivity) hβpos.le
  let f : (Fin d → ℤ) → (Fin (d - 2) → ℤ) := fun c i => (c ⟨i, by omega⟩ + H) / (b * b')
  have hmaps : Set.MapsTo f ((cube d H).filter (fun c => PhiK d a b a' b' c = z))
      ((Fintype.piFinset (fun _ : Fin (d - 2) => Icc (0:ℤ) K) : Finset _) :
        Set (Fin (d - 2) → ℤ)) := by
    intro c hc
    simp only [coe_filter, Set.mem_ofPred_eq, cube, Fintype.mem_piFinset, mem_Icc] at hc
    rw [Finset.mem_coe, Fintype.mem_piFinset]
    intro i
    simp only [f, mem_Icc]
    have := hc.1 ⟨i, by omega⟩
    exact ⟨Int.ediv_nonneg (by linarith) hβpos.le, Int.ediv_le_ediv hβpos (by linarith)⟩
  have hinj : Set.InjOn f ((cube d H).filter (fun c => PhiK d a b a' b' c = z)) := by
    intro c hc c' hc' hf
    simp only [coe_filter, Set.mem_ofPred_eq] at hc hc'
    have hk : PhiK d a b a' b' (c - c') = (0, 0) := by
      have h1 := hc.2
      have h2 := hc'.2
      simp only [PhiK] at h1 h2 ⊢
      rw [dotZ_sub, dotZ_sub]
      rw [← h2] at h1
      simp only [Prod.mk.injEq] at h1
      rw [h1.1, h1.2]; simp
    obtain ⟨e, he, hce⟩ := (kern_kernel d hd a b a' b' hb hb' hP hP' hΔ (c - c')).mp hk
    have hall : ∀ i, e i = 0 := by
      intro i
      induction i using Nat.strong_induction_on with
      | _ i ih =>
        by_cases hi : d - 2 ≤ i
        · exact he i hi
        · have hlt : i < d := by omega
          have h1 := hce ⟨i, hlt⟩
          simp only [Pi.sub_apply] at h1
          have hp1 : (if 1 ≤ i then e (i - 1) else 0) = 0 := by
            split_ifs
            · exact ih _ (by omega)
            · rfl
          have hp2 : (if 2 ≤ i then e (i - 2) else 0) = 0 := by
            split_ifs
            · exact ih _ (by omega)
            · rfl
          rw [hp1, hp2] at h1
          have hfi := congrFun hf ⟨i, by omega⟩
          simp only [f] at hfi
          have e1 : c ⟨i, hlt⟩ + H = (c' ⟨i, hlt⟩ + H) + b * b' * e i := by linarith
          rw [e1, Int.add_mul_ediv_left _ _ hβne] at hfi
          linarith
    funext i
    have h1 := hce i
    simp only [hall, Pi.sub_apply, mul_zero, sub_zero, add_zero, ite_self] at h1
    linarith
  have hcard := Finset.card_le_card_of_injOn f hmaps hinj
  rw [Fintype.card_piFinset, prod_const, Int.card_Icc, card_univ, Fintype.card_fin] at hcard
  have hcard' : ((((cube d H).filter (fun c => PhiK d a b a' b' c = z)).card : ℕ) : ℤ) ≤
      (((K + 1 - 0).toNat ^ (d - 2) : ℕ) : ℤ) := by exact_mod_cast hcard
  rw [Nat.cast_pow, Int.toNat_of_nonneg (by omega), sub_zero] at hcard'
  exact hcard'

lemma card_cube (d H : ℕ) : ((cube d H).card : ℝ) = (2 * H + 1) ^ d := by
  rw [cube, Fintype.card_piFinset, prod_const, Int.card_Icc, card_univ, Fintype.card_fin]
  have : ((H:ℤ) + 1 - -(H:ℤ)).toNat = 2 * H + 1 := by omega
  rw [this]; push_cast; ring

/-- MAIN. Lemma Kern-U (lv5unif §3), probability: if `bb′ d ≤ 2H + 1` then
`Pr_{c ∈ C_H}[Φ(c) = z] ≤ e/((2H+1)² (bb′)^{d−2})`. -/
theorem kernU_prob (d : ℕ) (hd : 3 ≤ d) (a b a' b' : ℤ) (hb : 1 ≤ b) (hb' : 1 ≤ b')
    (hP : IsCoprime a b) (hP' : IsCoprime a' b') (hΔ : a * b' - a' * b ≠ 0) (H : ℕ)
    (hH : b * b' * d ≤ 2 * (H : ℤ) + 1) (z : ℤ × ℤ) :
    prob (cube d H) (fun c => PhiK d a b a' b' c = z) ≤
      Real.exp 1 / ((2 * (H : ℝ) + 1) ^ 2 * ((b * b' : ℤ) : ℝ) ^ (d - 2)) := by
  have hcount := kernU_count d hd a b a' b' hb hb' hP hP' hΔ H z
  have hβpos : 0 < b * b' := mul_pos (by omega) (by omega)
  have hβne : b * b' ≠ 0 := hβpos.ne'
  have hKle := Int.ediv_mul_le (2 * (H:ℤ)) hβne
  have hK0 : 0 ≤ 2 * (H:ℤ) / (b * b') := Int.ediv_nonneg (by positivity) hβpos.le
  obtain ⟨k, rfl⟩ : ∃ k, d = k + 2 := ⟨d - 2, by omega⟩
  rw [Nat.add_sub_cancel] at hcount ⊢
  rw [prob, card_cube]
  set β : ℝ := ((b * b' : ℤ) : ℝ) with hβ
  set Kr : ℝ := ((2 * (H:ℤ) / (b * b') : ℤ) : ℝ) with hKr
  have hβr : 1 ≤ β := by rw [hβ]; exact_mod_cast hβpos
  have hKr0 : 0 ≤ Kr := by rw [hKr]; exact_mod_cast hK0
  have hKβ : Kr * β ≤ 2 * H := by
    rw [hKr, hβ]; exact_mod_cast hKle
  set x : ℝ := 2 * (H:ℝ) + 1 with hx
  have hx0 : 0 < x := by positivity
  have hHr : β * (k + 2 : ℝ) ≤ x := by
    rw [hβ, hx]; exact_mod_cast hH
  have hc : ((((cube (k + 2) H).filter (fun c => PhiK (k + 2) a b a' b' c = z)).card : ℝ)) ≤
      (Kr + 1) ^ k := by
    rw [hKr]; exact_mod_cast hcount
  have hK1 : Kr + 1 ≤ (x * ((2 * H + β) / x)) / β := by
    rw [mul_div_cancel₀ _ hx0.ne', le_div_iff₀ (by linarith)]; linarith
  have hr : (2 * (H:ℝ) + β) / x ≤ 1 + 1 / ((k:ℝ) + 2) := by
    rw [div_le_iff₀ hx0]
    have : β ≤ x / ((k:ℝ) + 2) := by rw [le_div_iff₀ (by positivity)]; linarith
    have e : (1 + 1 / ((k:ℝ) + 2)) * x = x + x / ((k:ℝ) + 2) := by ring
    rw [e]; rw [hx] at this ⊢; linarith
  have hr0 : 0 ≤ (2 * (H:ℝ) + β) / x := by positivity
  have hexp : ((2 * (H:ℝ) + β) / x) ^ k ≤ Real.exp 1 := by
    calc ((2 * (H:ℝ) + β) / x) ^ k ≤ (Real.exp (1 / ((k:ℝ) + 2))) ^ k := by
          gcongr
          linarith [Real.add_one_le_exp (1 / ((k:ℝ) + 2))]
      _ = Real.exp (k * (1 / ((k:ℝ) + 2))) := (Real.exp_nat_mul _ _).symm
      _ ≤ Real.exp 1 := by
          apply Real.exp_le_exp.mpr
          rw [mul_one_div, div_le_one (by positivity)]; linarith
  calc _ ≤ (Kr + 1) ^ k / x ^ (k + 2) := by
        gcongr
        convert hc
    _ ≤ ((x * ((2 * H + β) / x)) / β) ^ k / x ^ (k + 2) := by
        gcongr
    _ = ((2 * (H:ℝ) + β) / x) ^ k / (x ^ 2 * β ^ k) := by
        rw [div_pow, mul_pow, pow_add]
        field_simp
    _ ≤ Real.exp 1 / (x ^ 2 * β ^ k) := by gcongr

/-- MAIN. [pp] Lemma lem:kern(iii) in the form used for Theorem A (explicit `C_d = 3^{d−2}/2^d`):
if `|a|, |b|, |a′|, |b′| ≤ V` and `H ≥ 4V²` then
`Pr_{c ∈ C_H}[Φ(c) = z] ≤ 3^{d−2}/(2^d H² (bb′)^{d−2})`. -/
theorem kern_prob (d : ℕ) (hd : 3 ≤ d) (a b a' b' : ℤ) (hb : 1 ≤ b) (hb' : 1 ≤ b')
    (hP : IsCoprime a b) (hP' : IsCoprime a' b') (hΔ : a * b' - a' * b ≠ 0) (H V : ℕ)
    (ha : |a| ≤ V) (hbV : b ≤ V) (ha' : |a'| ≤ V) (hbV' : b' ≤ V) (hH : 4 * V ^ 2 ≤ H)
    (z : ℤ × ℤ) :
    prob (cube d H) (fun c => PhiK d a b a' b' c = z) ≤
      3 ^ (d - 2) / (2 ^ d * (H : ℝ) ^ 2 * ((b * b' : ℤ) : ℝ) ^ (d - 2)) := by
  have hcount := kernU_count d hd a b a' b' hb hb' hP hP' hΔ H z
  have hβpos : 0 < b * b' := mul_pos (by omega) (by omega)
  have hβne : b * b' ≠ 0 := hβpos.ne'
  have hKle := Int.ediv_mul_le (2 * (H:ℤ)) hβne
  have hK0 : 0 ≤ 2 * (H:ℤ) / (b * b') := Int.ediv_nonneg (by positivity) hβpos.le
  have hβV : b * b' ≤ (V:ℤ) * V := mul_le_mul hbV hbV' (by omega) (by omega)
  have hHz : (V:ℤ) * V * 4 ≤ H := by
    have : ((4 * V ^ 2 : ℕ) : ℤ) ≤ H := by exact_mod_cast hH
    push_cast at this; nlinarith
  obtain ⟨k, rfl⟩ : ∃ k, d = k + 2 := ⟨d - 2, by omega⟩
  rw [Nat.add_sub_cancel] at hcount ⊢
  rw [prob, card_cube]
  set β : ℝ := ((b * b' : ℤ) : ℝ) with hβ
  set Kr : ℝ := ((2 * (H:ℤ) / (b * b') : ℤ) : ℝ) with hKr
  have hβr : 1 ≤ β := by rw [hβ]; exact_mod_cast hβpos
  have hKr0 : 0 ≤ Kr := by rw [hKr]; exact_mod_cast hK0
  have hKβ : Kr * β ≤ 2 * H := by
    rw [hKr, hβ]; exact_mod_cast hKle
  have hβH : β ≤ H := by
    rw [hβ]
    have : b * b' ≤ (H:ℤ) := by nlinarith
    exact_mod_cast this
  have hH0 : (0:ℝ) < H := by linarith
  have hc : ((((cube (k + 2) H).filter (fun c => PhiK (k + 2) a b a' b' c = z)).card : ℝ)) ≤
      (Kr + 1) ^ k := by
    rw [hKr]; exact_mod_cast hcount
  have hK1 : Kr + 1 ≤ 3 * H / β := by
    rw [le_div_iff₀ (by linarith)]; nlinarith
  calc _ ≤ (Kr + 1) ^ k / (2 * (H:ℝ) + 1) ^ (k + 2) := by
        gcongr
        convert hc
    _ ≤ (3 * H / β) ^ k / (2 * (H:ℝ) + 1) ^ (k + 2) := by
        gcongr
    _ ≤ (3 * H / β) ^ k / (2 * (H:ℝ)) ^ (k + 2) := by
        gcongr
        linarith
    _ = 3 ^ k / (2 ^ (k + 2) * (H : ℝ) ^ 2 * β ^ k) := by
        rw [div_pow, mul_pow, mul_pow, pow_add]
        field_simp
        try ring

/-- MAIN. [pp] Lemma lem:pairs. Let `ℓ` be prime, `g = (g₀, g₁)` coprime, and `P ≠ P′` points
`(a, b)` with `b > 0`, `gcd(a, b) = 1`, `|a|, |b| ≤ V`, both in `Λ_ℓ(g) = {ℓ ∣ g₁x − g₀y}`.
Then `Δ = det(P, P′) ≠ 0`, `ℓ ∣ Δ` and `|Δ| ≤ 2V²` (so `q = |Δ|/ℓ ≤ 2S²` when `V² ≤ ℓS²`). -/
theorem pairs (ℓ : ℕ) (hℓ : ℓ.Prime) (g0 g1 : ℤ) (hg : IsCoprime g0 g1) (V : ℕ)
    (a b a' b' : ℤ) (hb : 0 < b) (hb' : 0 < b') (hP : IsCoprime a b) (hP' : IsCoprime a' b')
    (ha : |a| ≤ V) (hbV : b ≤ V) (ha' : |a'| ≤ V) (hbV' : b' ≤ V)
    (hL : (ℓ : ℤ) ∣ Gform g0 g1 a b) (hL' : (ℓ : ℤ) ∣ Gform g0 g1 a' b') (hne : (a, b) ≠ (a', b')) :
    a * b' - a' * b ≠ 0 ∧ (ℓ : ℤ) ∣ a * b' - a' * b ∧ |a * b' - a' * b| ≤ 2 * (V : ℤ) ^ 2 := by
  refine ⟨?_, ?_, ?_⟩
  · intro h
    apply hne
    have h1 : b ∣ b' := by
      have : b ∣ a * b' := ⟨a', by linarith⟩
      exact hP.symm.dvd_of_dvd_mul_left this
    have h2 : b' ∣ b := by
      have : b' ∣ a' * b := ⟨a, by linarith⟩
      exact hP'.symm.dvd_of_dvd_mul_left this
    have hbb : b = b' := Int.dvd_antisymm hb.le hb'.le h1 h2
    subst hbb
    have : a = a' := mul_right_cancel₀ (ne_of_gt hb) (by linarith : a * b = a' * b)
    rw [this]
  · obtain ⟨u, v, huv⟩ := hg
    have e : a * b' - a' * b = u * (a' * Gform g0 g1 a b - a * Gform g0 g1 a' b') +
        v * (b' * Gform g0 g1 a b - b * Gform g0 g1 a' b') := by
      simp only [Gform]; linear_combination (a' * b - a * b') * huv
    rw [e]
    exact dvd_add (dvd_mul_of_dvd_right (dvd_sub (dvd_mul_of_dvd_right hL _)
      (dvd_mul_of_dvd_right hL' _)) _) (dvd_mul_of_dvd_right (dvd_sub (dvd_mul_of_dvd_right hL _)
      (dvd_mul_of_dvd_right hL' _)) _)
  · calc |a * b' - a' * b| ≤ |a * b'| + |a' * b| := abs_sub _ _
      _ = |a| * b' + |a'| * b := by rw [abs_mul, abs_mul, abs_of_pos hb', abs_of_pos hb]
      _ ≤ V * V + V * V := add_le_add (mul_le_mul ha hbV' hb'.le (by positivity))
          (mul_le_mul ha' hbV hb.le (by positivity))
      _ = 2 * (V : ℤ) ^ 2 := by ring

/-- MAIN. [pp] Lemma lem:pairs, multiplicity: for `P = (a, b)` with `2b > V` and `D ∈ ℤ`, at most
`6` points `P′ ∈ [−V, V]²` have `det(P, P′) = D`. -/
theorem pairs_count (V : ℕ) (a b : ℤ) (hb : (V : ℤ) < 2 * b) (hP : IsCoprime a b) (D : ℤ) :
    ((((Icc (-(V : ℤ)) V) ×ˢ (Icc (-(V : ℤ)) V)).filter
      (fun Q : ℤ × ℤ => a * Q.2 - Q.1 * b = D)).card) ≤ 6 := by
  have hb0 : 0 < b := by omega
  have hbne : b ≠ 0 := hb0.ne'
  have hle : ((((Icc (-(V : ℤ)) V) ×ˢ (Icc (-(V : ℤ)) V)).filter
      (fun Q : ℤ × ℤ => a * Q.2 - Q.1 * b = D)).card) ≤ (Icc (0:ℤ) 3).card := by
    apply Finset.card_le_card_of_injOn (fun Q : ℤ × ℤ => (Q.2 + V) / b)
    · intro Q hQ
      simp only [coe_filter, mem_product, mem_Icc, Set.mem_ofPred_eq] at hQ
      simp only [coe_Icc, Set.mem_Icc]
      refine ⟨Int.ediv_nonneg (by linarith) hb0.le, ?_⟩
      have : (Q.2 + V) / b < 4 := Int.ediv_lt_of_lt_mul hb0 (by linarith)
      omega
    · intro Q hQ Q' hQ' hf
      simp only [coe_filter, mem_product, mem_Icc, Set.mem_ofPred_eq] at hQ hQ'
      simp only at hf
      have hd : b ∣ a * (Q'.2 - Q.2) := ⟨Q'.1 - Q.1, by linarith⟩
      obtain ⟨j, hj⟩ := hP.symm.dvd_of_dvd_mul_left hd
      have e1 : Q'.2 + V = (Q.2 + V) + b * j := by linarith
      rw [e1, Int.add_mul_ediv_left _ _ hbne] at hf
      have hj0 : j = 0 := by linarith
      subst hj0
      have h2 : Q'.2 = Q.2 := by linarith
      have h1 : b * (Q'.1 - Q.1) = 0 := by linear_combination hQ.2 - hQ'.2 + a * h2
      have h1' : Q'.1 = Q.1 := by
        have := (mul_eq_zero.mp h1).resolve_left hbne; linarith
      exact Prod.ext h1'.symm h2.symm
  have h4 : (Icc (0:ℤ) 3).card = 4 := by simp
  omega

end FactoringA
