import LeanA.Defs

/-!
# Orthogonality and lift averaging for smooth numbers in progressions

Source: [pp] (`paper2/main.tex`), Appendix `app:85`, Lemma `lem:E1` (orthogonality) and
Lemma `lem:E2` (lift averaging), the elementary identity behind Part 1 of Lemma 8.5′.

* `S(z, y)` = `smoothSet z y` = the `y`-smooth `n ∈ [1, z]`.
* `Ψ(z, y; q, a)` = `PsiAP`, `Ψ_q(z, y)` = `PsiCop`, `E(z, y; q, a)` = `Eerr`,
  `Ψ(z, y; χ)` = `PsiChi`, `Δ(q)` = `DeltaAll`, `Δ_ℓ(t)` = `DeltaEll`.

The definitions marked TRUSTED and the theorems marked MAIN are copied into `Challenge.lean`.
-/

set_option autoImplicit false

noncomputable section
open Classical

namespace FactoringA

open Finset

/-! ## Definitions -/

/-- TRUSTED. `n` is `y`-smooth: every prime factor `p` of `n` satisfies `p ≤ y`
(used only for `n ≥ 1`). -/
def IsYSmooth (y : ℝ) (n : ℕ) : Prop := ∀ p : ℕ, p.Prime → p ∣ n → (p : ℝ) ≤ y

/-- TRUSTED. `S(z, y)`: the `y`-smooth integers `1 ≤ n ≤ z`. -/
def smoothSet (z y : ℝ) : Finset ℕ := (Finset.Icc 1 ⌊z⌋₊).filter (IsYSmooth y)

/-- TRUSTED. `Ψ(z, y; q, a) = #{n ∈ S(z, y) : n ≡ a (mod q)}`. -/
def PsiAP (z y : ℝ) (q a : ℕ) : ℕ := ((smoothSet z y).filter (fun n => n ≡ a [MOD q])).card

/-- TRUSTED. `Ψ_q(z, y) = #{n ∈ S(z, y) : gcd(n, q) = 1}`. -/
def PsiCop (z y : ℝ) (q : ℕ) : ℕ := ((smoothSet z y).filter (fun n => Nat.Coprime n q)).card

/-- TRUSTED. `E(z, y; q, a) = Ψ(z, y; q, a) − Ψ_q(z, y)/φ(q)`. -/
def Eerr (z y : ℝ) (q a : ℕ) : ℝ := (PsiAP z y q a : ℝ) - (PsiCop z y q : ℝ) / (q.totient : ℝ)

/-- TRUSTED. `Ψ(z, y; χ) = ∑_{n ∈ S(z, y)} χ(n)` for a Dirichlet character `χ` mod `q`. -/
def PsiChi (z y : ℝ) {q : ℕ} (χ : DirichletCharacter ℂ q) : ℂ :=
  ∑ n ∈ smoothSet z y, χ (n : ZMod q)

/-- TRUSTED. `Δ_ℓ(t) = φ(t)⁻¹ ∑_{χ mod t, ℓ ∣ cond χ} |Ψ(z, y; χ)|`. -/
def DeltaEll (z y : ℝ) (ℓ t : ℕ) : ℝ :=
  (1 / (t.totient : ℝ)) *
    ∑ χ ∈ (univ : Finset (DirichletCharacter ℂ t)).filter (fun χ => ℓ ∣ χ.conductor),
      ‖PsiChi z y χ‖

/-- TRUSTED. `Δ(q) = φ(q)⁻¹ ∑_{χ mod q, χ ≠ χ₀} |Ψ(z, y; χ)|`. -/
def DeltaAll (z y : ℝ) (q : ℕ) : ℝ :=
  (1 / (q.totient : ℝ)) *
    ∑ χ ∈ (univ : Finset (DirichletCharacter ℂ q)).filter (fun χ => χ ≠ 1), ‖PsiChi z y χ‖

/-! ## Helper lemmas -/

lemma isUnit_inv_natCast {q a : ℕ} (ha : Nat.Coprime a q) : IsUnit ((a : ZMod q)⁻¹) := by
  have hu : IsUnit (a : ZMod q) := (ZMod.isUnit_iff_coprime a q).2 ha
  exact isUnit_iff_exists_inv.2 ⟨_, ZMod.inv_mul_of_unit _ hu⟩

/-- Full orthogonality: `∑_{χ mod q} χ̄(a) Ψ(z, y; χ) = φ(q) Ψ(z, y; q, a)`. -/
lemma sum_all_chars_PsiChi (z y : ℝ) (q a : ℕ) [NeZero q] (ha : Nat.Coprime a q) :
    ∑ χ : DirichletCharacter ℂ q, χ ((a : ZMod q)⁻¹) * PsiChi z y χ =
      (q.totient : ℂ) * (PsiAP z y q a : ℂ) := by
  have hu : IsUnit (a : ZMod q) := (ZMod.isUnit_iff_coprime a q).2 ha
  simp only [PsiChi, Finset.mul_sum]
  rw [Finset.sum_comm]
  simp_rw [DirichletCharacter.sum_char_inv_mul_char_eq ℂ hu]
  rw [PsiAP, Finset.card_filter, Nat.cast_sum, Finset.mul_sum]
  refine Finset.sum_congr rfl fun n _ => ?_
  have : ((a : ZMod q) = (n : ZMod q)) ↔ n ≡ a [MOD q] := by
    rw [ZMod.natCast_eq_natCast_iff]; exact Nat.ModEq.comm
  by_cases h : n ≡ a [MOD q]
  · simp [this.2 h, h]
  · have h' : ¬ ((a : ZMod q) = (n : ZMod q)) := fun e => h (this.1 e)
    simp [h, h']

/-- `Ψ(z, y; χ₀) = Ψ_q(z, y)`. -/
lemma PsiChi_one (z y : ℝ) (q : ℕ) :
    PsiChi z y (1 : DirichletCharacter ℂ q) = (PsiCop z y q : ℂ) := by
  rw [PsiChi, PsiCop, Finset.card_filter, Nat.cast_sum]
  refine Finset.sum_congr rfl fun n _ => ?_
  by_cases h : IsUnit (n : ZMod q)
  · have hc : Nat.Coprime n q := (ZMod.isUnit_iff_coprime n q).1 h
    simp [MulChar.one_apply h, hc]
  · have hc : ¬ Nat.Coprime n q := fun hc => h ((ZMod.isUnit_iff_coprime n q).2 hc)
    simp [MulChar.map_nonunit _ h, hc]

lemma totient_pos_real {q : ℕ} (hq : 1 ≤ q) : (0 : ℝ) < q.totient := by
  exact_mod_cast Nat.totient_pos.2 (by omega)

/-- Every `y`-smooth `n` is coprime to a prime `ℓ > y`. -/
lemma coprime_of_mem_smoothSet {z y : ℝ} {ℓ n : ℕ} (hℓ : ℓ.Prime) (hyℓ : y < ℓ)
    (hn : n ∈ smoothSet z y) : Nat.Coprime n ℓ := by
  rw [smoothSet, Finset.mem_filter] at hn
  refine Nat.Coprime.symm ((Nat.Prime.coprime_iff_not_dvd hℓ).2 fun h => ?_)
  exact absurd (hn.2 ℓ hℓ h) (not_le.2 hyℓ)

/-- For `t = ℓ s` with `ℓ ∤ s`: `ℓ ∤ cond χ` iff `χ` factors through level `s`. -/
lemma not_dvd_conductor_iff {ℓ s : ℕ} [NeZero (ℓ * s)] (hℓ : ℓ.Prime) (hℓs : ¬ ℓ ∣ s)
    (χ : DirichletCharacter ℂ (ℓ * s)) :
    ¬ ℓ ∣ χ.conductor ↔ χ.FactorsThrough s := by
  have hsd : s ∣ ℓ * s := dvd_mul_left s ℓ
  rw [← DirichletCharacter.mem_conductorSet_iff,
    DirichletCharacter.mem_conductorSet_iff_conductor_dvd χ hsd]
  constructor
  · intro h
    exact Nat.Coprime.dvd_of_dvd_mul_left
      ((Nat.Prime.coprime_iff_not_dvd hℓ).2 h).symm (DirichletCharacter.conductor_dvd_level χ)
  · intro h hd
    exact hℓs (hd.trans h)

/-- The characters mod `ℓ s` with `ℓ ∤ cond` are the lifts of the characters mod `s`. -/
lemma filter_not_dvd_conductor_eq_image {ℓ s : ℕ} [NeZero (ℓ * s)] (hℓ : ℓ.Prime)
    (hℓs : ¬ ℓ ∣ s) :
    (univ : Finset (DirichletCharacter ℂ (ℓ * s))).filter (fun χ => ¬ ℓ ∣ χ.conductor) =
      univ.image (DirichletCharacter.changeLevel (R := ℂ) (dvd_mul_left s ℓ)) := by
  ext χ
  simp only [Finset.mem_filter, Finset.mem_univ, true_and, Finset.mem_image]
  rw [not_dvd_conductor_iff hℓ hℓs χ]
  constructor
  · rintro ⟨_, χ₀, rfl⟩; exact ⟨χ₀, rfl⟩
  · rintro ⟨χ₀, rfl⟩; exact ⟨dvd_mul_left s ℓ, χ₀, rfl⟩

/-- Orthogonality over the characters mod `ℓ s` with `ℓ ∤ cond`, evaluated at a unit `x`. -/
lemma sum_not_dvd_conductor_apply {ℓ s : ℕ} [NeZero s] [NeZero (ℓ * s)] (hℓ : ℓ.Prime)
    (hℓs : ¬ ℓ ∣ s) (x : ZMod (ℓ * s)) (hx : IsUnit x) :
    ∑ χ ∈ (univ : Finset (DirichletCharacter ℂ (ℓ * s))).filter (fun χ => ¬ ℓ ∣ χ.conductor),
        χ x = if (ZMod.cast x : ZMod s) = 1 then (s.totient : ℂ) else 0 := by
  rw [filter_not_dvd_conductor_eq_image hℓ hℓs,
    Finset.sum_image (fun a _ b _ h => DirichletCharacter.changeLevel_injective _ h)]
  have : ∀ ψ : DirichletCharacter ℂ s,
      DirichletCharacter.changeLevel (dvd_mul_left s ℓ) ψ x = ψ (ZMod.cast x) := by
    intro ψ
    have := DirichletCharacter.changeLevel_eq_cast_of_dvd ψ (dvd_mul_left s ℓ) hx.unit
    rwa [IsUnit.unit_spec] at this
  simp_rw [this]
  exact DirichletCharacter.sum_characters_eq ℂ _

/-- The pointwise identity: for `y`-smooth `n` and `(r, ℓ s) = 1`,
`∑_{χ mod ℓs, ℓ ∤ cond χ} χ̄(r) χ(n) = φ(s) · 1[n ≡ r (mod s)]`. -/
lemma sum_not_dvd_conductor_inv_mul {z y : ℝ} {ℓ s r n : ℕ} [NeZero s] [NeZero (ℓ * s)]
    (hℓ : ℓ.Prime) (hyℓ : y < ℓ) (hℓs : ¬ ℓ ∣ s) (hr : Nat.Coprime r (ℓ * s))
    (hn : n ∈ smoothSet z y) :
    ∑ χ ∈ (univ : Finset (DirichletCharacter ℂ (ℓ * s))).filter (fun χ => ¬ ℓ ∣ χ.conductor),
        χ ((r : ZMod (ℓ * s))⁻¹) * χ (n : ZMod (ℓ * s)) =
      if n ≡ r [MOD s] then (s.totient : ℂ) else 0 := by
  have hsd : s ∣ ℓ * s := dvd_mul_left s ℓ
  simp_rw [← map_mul]
  have hri := isUnit_inv_natCast hr
  have hnℓ := coprime_of_mem_smoothSet hℓ hyℓ hn
  by_cases hns : Nat.Coprime n s
  · have hnu : IsUnit (n : ZMod (ℓ * s)) :=
      (ZMod.isUnit_iff_coprime n (ℓ * s)).2 (Nat.Coprime.mul_right hnℓ hns)
    rw [sum_not_dvd_conductor_apply hℓ hℓs _ (hri.mul hnu)]
    congr 1
    apply propext
    rw [ZMod.cast_mul hsd, ZMod.cast_natCast hsd]
    set c : ZMod s := ZMod.cast ((r : ZMod (ℓ * s))⁻¹)
    have hcr : c * (r : ZMod s) = 1 := by
      have hru : IsUnit (r : ZMod (ℓ * s)) := (ZMod.isUnit_iff_coprime r _).2 hr
      rw [← ZMod.cast_natCast hsd (R := ZMod s) r, ← ZMod.cast_mul hsd,
        ZMod.inv_mul_of_unit _ hru, ZMod.cast_one hsd]
    have hc : IsUnit c := isUnit_iff_exists_inv.2 ⟨_, hcr⟩
    rw [← hcr, hc.mul_right_inj, ZMod.natCast_eq_natCast_iff]
  · have h0 : ∀ χ : DirichletCharacter ℂ (ℓ * s),
        χ ((r : ZMod (ℓ * s))⁻¹ * (n : ZMod (ℓ * s))) = 0 := by
      intro χ
      refine MulChar.map_nonunit _ fun hu => hns ?_
      have := (ZMod.isUnit_iff_coprime n (ℓ * s)).1 (isUnit_of_mul_isUnit_right hu)
      exact Nat.Coprime.coprime_mul_left_right this
    have hmod : ¬ n ≡ r [MOD s] := by
      intro hmod
      apply hns
      have hrs : Nat.Coprime r s := Nat.Coprime.coprime_mul_left_right hr
      rw [Nat.Coprime, hmod.gcd_eq]; exact hrs
    simp_rw [h0, Finset.sum_const_zero]
    simp [hmod]

/-- Summed over `S(z, y)`: `∑_{χ mod ℓs, ℓ ∤ cond χ} χ̄(r) Ψ(z, y; χ) = φ(s) Ψ(z, y; s, r)`. -/
lemma sum_not_dvd_conductor_PsiChi (z y : ℝ) {ℓ s r : ℕ} [NeZero s] [NeZero (ℓ * s)]
    (hℓ : ℓ.Prime) (hyℓ : y < ℓ) (hℓs : ¬ ℓ ∣ s) (hr : Nat.Coprime r (ℓ * s)) :
    ∑ χ ∈ (univ : Finset (DirichletCharacter ℂ (ℓ * s))).filter (fun χ => ¬ ℓ ∣ χ.conductor),
        χ ((r : ZMod (ℓ * s))⁻¹) * PsiChi z y χ =
      (s.totient : ℂ) * (PsiAP z y s r : ℂ) := by
  simp only [PsiChi, Finset.mul_sum]
  rw [Finset.sum_comm]
  rw [Finset.sum_congr rfl fun n hn => sum_not_dvd_conductor_inv_mul hℓ hyℓ hℓs hr hn]
  rw [PsiAP, Finset.card_filter, Nat.cast_sum, Finset.mul_sum]
  refine Finset.sum_congr rfl fun n _ => ?_
  split_ifs <;> simp

/-- `Ψ_{ℓ s}(z, y) = Ψ_s(z, y)` for a prime `ℓ > y`. -/
lemma PsiCop_mul_prime (z y : ℝ) {ℓ s : ℕ} (hℓ : ℓ.Prime) (hyℓ : y < ℓ) :
    PsiCop z y (ℓ * s) = PsiCop z y s := by
  rw [PsiCop, PsiCop]
  congr 1
  refine Finset.filter_congr fun n hn => ?_
  have hnℓ := coprime_of_mem_smoothSet hℓ hyℓ hn
  constructor
  · exact Nat.Coprime.coprime_mul_left_right
  · exact Nat.Coprime.mul_right hnℓ

/-! ## Main theorems -/

/-- MAIN. [pp] Lemma lem:E1 (orthogonality): for `q ≥ 1` and `(a, q) = 1`,
`E(z, y; q, a) = φ(q)⁻¹ ∑_{χ ≠ χ₀} χ̄(a) Ψ(z, y; χ)`, where `χ̄(a) = χ(a⁻¹)`. -/
theorem orthogonality_E1 (z y : ℝ) (q a : ℕ) (hq : 1 ≤ q) (ha : Nat.Coprime a q) :
    ((Eerr z y q a : ℝ) : ℂ) =
      (1 / (q.totient : ℂ)) *
        ∑ χ ∈ (univ : Finset (DirichletCharacter ℂ q)).filter (fun χ => χ ≠ 1),
          χ ((a : ZMod q)⁻¹) * PsiChi z y χ := by
  have hfe : (univ : Finset (DirichletCharacter ℂ q)).filter (fun χ => χ ≠ 1) = univ.erase 1 := by
    ext χ; simp
  rw [hfe, Finset.sum_erase_eq_sub (Finset.mem_univ _)]
  have : NeZero q := ⟨by omega⟩
  rw [sum_all_chars_PsiChi z y q a ha, PsiChi_one, MulChar.one_apply (isUnit_inv_natCast ha)]
  have hφ : (q.totient : ℂ) ≠ 0 := by
    exact_mod_cast (Nat.totient_pos.2 (by omega : 0 < q)).ne'
  rw [Eerr]
  push_cast
  field_simp

/-- MAIN. [pp] Lemma lem:E1, second part: `|E(z, y; q, a)| ≤ Δ(q)` for `q ≥ 1`, `(a, q) = 1`. -/
theorem abs_Eerr_le_DeltaAll (z y : ℝ) (q a : ℕ) (hq : 1 ≤ q) (ha : Nat.Coprime a q) :
    |Eerr z y q a| ≤ DeltaAll z y q := by
  have h := congrArg (‖·‖) (orthogonality_E1 z y q a hq ha)
  simp only [Complex.norm_real, Real.norm_eq_abs] at h
  rw [h, DeltaAll, norm_mul]
  have hφ := totient_pos_real hq
  have h1 : ‖(1 / (q.totient : ℂ))‖ = 1 / (q.totient : ℝ) := by
    rw [norm_div, norm_one, Complex.norm_natCast]
  rw [h1]
  refine mul_le_mul_of_nonneg_left ?_ (by positivity)
  refine (norm_sum_le _ _).trans (Finset.sum_le_sum fun χ _ => ?_)
  rw [norm_mul]
  exact mul_le_of_le_one_left (norm_nonneg _) (DirichletCharacter.norm_le_one _ _)

/-- MAIN. [pp] Lemma lem:E2 (lift averaging), exact form: let `ℓ > y` be prime, `s ≥ 1`,
`ℓ ∤ s`, `t = ℓ s`, `(r, t) = 1`. Then
`E(z, y; t, r) − E(z, y; s, r)/(ℓ − 1) = φ(t)⁻¹ ∑_{χ mod t, ℓ ∣ cond χ} χ̄(r) Ψ(z, y; χ)`. -/
theorem lift_averaging_E2_eq (z y : ℝ) (ℓ s t r : ℕ) (hℓ : ℓ.Prime) (hyℓ : y < ℓ)
    (hs : 1 ≤ s) (hℓs : ¬ ℓ ∣ s) (ht : t = ℓ * s) (hr : Nat.Coprime r t) :
    ((Eerr z y t r - Eerr z y s r / ((ℓ : ℝ) - 1) : ℝ) : ℂ) =
      (1 / (t.totient : ℂ)) *
        ∑ χ ∈ (univ : Finset (DirichletCharacter ℂ t)).filter (fun χ => ℓ ∣ χ.conductor),
          χ ((r : ZMod t)⁻¹) * PsiChi z y χ := by
  subst ht
  have : NeZero s := ⟨by omega⟩
  have : NeZero (ℓ * s) := ⟨Nat.mul_ne_zero hℓ.ne_zero (by omega)⟩
  have hall := sum_all_chars_PsiChi z y (ℓ * s) r hr
  rw [← Finset.sum_filter_add_sum_filter_not univ (fun χ : DirichletCharacter ℂ (ℓ * s) =>
    ℓ ∣ χ.conductor), sum_not_dvd_conductor_PsiChi z y hℓ hyℓ hℓs hr] at hall
  have hB := eq_sub_of_add_eq hall
  rw [hB, Eerr, Eerr, PsiCop_mul_prime z y hℓ hyℓ,
    Nat.totient_mul ((Nat.Prime.coprime_iff_not_dvd hℓ).2 hℓs), Nat.totient_prime hℓ]
  have hℓ1 : (1 : ℕ) ≤ ℓ := hℓ.one_lt.le
  have hℓc : ((ℓ : ℂ) - 1) ≠ 0 := by
    have : (1 : ℝ) < ℓ := by exact_mod_cast hℓ.one_lt
    intro h
    have h' : ((ℓ : ℝ) : ℂ) = ((1 : ℝ) : ℂ) := by push_cast; linear_combination h
    exact this.ne' (Complex.ofReal_injective h')
  have hφs : (s.totient : ℂ) ≠ 0 := by
    exact_mod_cast (Nat.totient_pos.2 (by omega : 0 < s)).ne'
  push_cast [Nat.cast_sub hℓ1]
  field_simp
  ring

/-- MAIN. [pp] Lemma lem:E2 (lift averaging): let `ℓ > y` be prime, `s ≥ 1`, `ℓ ∤ s`,
`t = ℓ s`. For every `(r, t) = 1`, `|E(z, y; t, r) − E(z, y; s, r)/(ℓ − 1)| ≤ Δ_ℓ(t)`. -/
theorem lift_averaging_E2 (z y : ℝ) (ℓ s t r : ℕ) (hℓ : ℓ.Prime) (hyℓ : y < ℓ)
    (hs : 1 ≤ s) (hℓs : ¬ ℓ ∣ s) (ht : t = ℓ * s) (hr : Nat.Coprime r t) :
    |Eerr z y t r - Eerr z y s r / ((ℓ : ℝ) - 1)| ≤ DeltaEll z y ℓ t := by
  have h := congrArg (‖·‖) (lift_averaging_E2_eq z y ℓ s t r hℓ hyℓ hs hℓs ht hr)
  simp only [Complex.norm_real, Real.norm_eq_abs] at h
  rw [h, DeltaEll, norm_mul]
  have ht1 : 1 ≤ t := by
    subst ht; exact Nat.one_le_iff_ne_zero.2 (Nat.mul_ne_zero hℓ.ne_zero (by omega))
  have hφ := totient_pos_real ht1
  have h1 : ‖(1 / (t.totient : ℂ))‖ = 1 / (t.totient : ℝ) := by
    rw [norm_div, norm_one, Complex.norm_natCast]
  rw [h1]
  refine mul_le_mul_of_nonneg_left ?_ (by positivity)
  refine (norm_sum_le _ _).trans (Finset.sum_le_sum fun χ _ => ?_)
  rw [norm_mul]
  exact mul_le_of_le_one_left (norm_nonneg _) (DirichletCharacter.norm_le_one _ _)

end FactoringA
