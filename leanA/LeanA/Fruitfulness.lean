import LeanA.Defs

/-!
# (b) Theorem DD (Dixon descent; Rabin/Dixon fibre symmetry), [pp] Theorem thm:DD

Model. The hidden roots `r = (r_1, …, r_n)` are uniform on `((ℤ/N)ˣ)ⁿ` (every row of SPLIT has a
unit root: a non-unit root is a factor exit). Everything the algorithm exposes is a function of the
squares `R_i = r_i²` and of independent coins `ξ`; this is hypothesis (T) of the paper and is built
in here by letting the dependency `S` and the value `Y` be functions of `(R, ξ)` only.
"Conditional on `𝓕 = σ(R, ξ)`" is expressed by counting inside one fibre `{r : r_i² = R_i ∀ i}`
with `ξ` fixed.

The main statements (MAIN) are copied into `Challenge.lean`; do not change them.
-/

noncomputable section
open Classical

namespace FactoringA

open Finset

lemma gcd_val_unit_mul (N : ℕ) [NeZero N] (u : (ZMod N)ˣ) (x : ZMod N) :
    Nat.gcd ((u : ZMod N) * x).val N = Nat.gcd x.val N := by
  rw [ZMod.val_mul, ← Nat.gcd_rec N, Nat.gcd_comm,
    Nat.Coprime.gcd_mul_left_cancel _ (ZMod.val_coe_unit_coprime u)]

lemma properGcd_unit_mul (N : ℕ) [NeZero N] (u : (ZMod N)ˣ) (x : ZMod N) :
    ProperGcd N ((u : ZMod N) * x) ↔ ProperGcd N x := by
  unfold ProperGcd; rw [gcd_val_unit_mul]

lemma sq_eq_one_primePow (p e : ℕ) (hp : p.Prime) (hp2 : p ≠ 2) (x : ZMod (p ^ e))
    (hx : x ^ 2 = 1) : x = 1 ∨ x = -1 := by
  have : NeZero (p ^ e) := ⟨pow_ne_zero _ hp.ne_zero⟩
  set a : ℤ := (x.val : ℤ) with ha
  have hxa : ((a : ℤ) : ZMod (p ^ e)) = x := by simp [a]
  have hdvd : ((p : ℤ) ^ e) ∣ (a - 1) * (a + 1) := by
    have := (ZMod.intCast_zmod_eq_zero_iff_dvd ((a - 1) * (a + 1)) (p ^ e)).mp (by
      push_cast; rw [hxa]; linear_combination hx)
    exact_mod_cast this
  have hpZ : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  have hnot : ¬ ((p : ℤ) ∣ a - 1 ∧ (p : ℤ) ∣ a + 1) := by
    rintro ⟨h1, h2⟩
    have h3 : (p : ℤ) ∣ 2 := by
      have := dvd_sub h2 h1
      have e2 : a + 1 - (a - 1) = 2 := by ring
      rwa [e2] at this
    have h4 : p ∣ 2 := by exact_mod_cast h3
    exact hp2 ((Nat.prime_dvd_prime_iff_eq hp Nat.prime_two).mp h4)
  by_cases h1 : (p : ℤ) ∣ a - 1
  · have h2 : ¬ (p : ℤ) ∣ a + 1 := fun h2 => hnot ⟨h1, h2⟩
    left
    have h := hpZ.pow_dvd_of_dvd_mul_right e h2 hdvd
    have h' := (ZMod.intCast_zmod_eq_zero_iff_dvd (a - 1) (p ^ e)).mpr (by exact_mod_cast h)
    push_cast at h'
    rw [hxa] at h'
    linear_combination h'
  · right
    have h := hpZ.pow_dvd_of_dvd_mul_left e h1 hdvd
    have h' := (ZMod.intCast_zmod_eq_zero_iff_dvd (a + 1) (p ^ e)).mpr (by exact_mod_cast h)
    push_cast at h'
    rw [hxa] at h'
    linear_combination h'

lemma natCard_sqOne_congr {A B : Type*} [CommRing A] [CommRing B] (e : A ≃+* B) :
    Nat.card {u : Aˣ // u ^ 2 = 1} = Nat.card {u : Bˣ // u ^ 2 = 1} :=
  Nat.card_congr ((Units.mapEquiv e.toMulEquiv).toEquiv.subtypeEquiv (fun u => by
    show u ^ 2 = 1 ↔ (Units.mapEquiv e.toMulEquiv u) ^ 2 = 1
    rw [← map_pow, MulEquiv.map_eq_one_iff]))

lemma natCard_sqOne_prod {A B : Type*} [CommRing A] [CommRing B] :
    Nat.card {u : (A × B)ˣ // u ^ 2 = 1} =
      Nat.card {u : Aˣ // u ^ 2 = 1} * Nat.card {u : Bˣ // u ^ 2 = 1} := by
  rw [← Nat.card_prod]
  apply Nat.card_congr
  refine ((MulEquiv.prodUnits).toEquiv.subtypeEquiv
    (q := fun c : Aˣ × Bˣ => c.1 ^ 2 = 1 ∧ c.2 ^ 2 = 1) ?_).trans
    (Equiv.subtypeProdEquivProd (p := fun a : Aˣ => a ^ 2 = 1) (q := fun b : Bˣ => b ^ 2 = 1))
  intro u
  show u ^ 2 = 1 ↔ (MulEquiv.prodUnits u).1 ^ 2 = 1 ∧ (MulEquiv.prodUnits u).2 ^ 2 = 1
  rw [← MulEquiv.map_eq_one_iff (MulEquiv.prodUnits (M := A) (N := B)), map_pow, Prod.ext_iff]
  simp

lemma natCard_sqOne_primePow (p e : ℕ) (hp : p.Prime) (hp2 : p ≠ 2) (he : 0 < e) :
    Nat.card {u : (ZMod (p ^ e))ˣ // u ^ 2 = 1} = 2 := by
  have h3 : 2 < p ^ e := by
    have := hp.two_le
    have : p ≤ p ^ e := Nat.le_self_pow (by omega) p
    omega
  have : Fact (2 < p ^ e) := ⟨h3⟩
  have hne : (1 : (ZMod (p ^ e))ˣ) ≠ -1 := by
    intro h
    have := congrArg Units.val h
    simp only [Units.val_one, Units.val_neg] at this
    exact ZMod.neg_one_ne_one this.symm
  have hS : {u : (ZMod (p ^ e))ˣ | u ^ 2 = 1} = {1, -1} := by
    ext u
    simp only [Set.mem_ofPred_eq, Set.mem_insert_iff, Set.mem_singleton_iff]
    constructor
    · intro hu
      have hu' : (u : ZMod (p ^ e)) ^ 2 = 1 := by
        rw [← Units.val_pow_eq_pow_val, hu, Units.val_one]
      rcases sq_eq_one_primePow p e hp hp2 _ hu' with h | h
      · left; ext; simpa using h
      · right; ext; simpa using h
    · rintro (rfl | rfl) <;> simp
  calc Nat.card {u : (ZMod (p ^ e))ˣ // u ^ 2 = 1}
      = Nat.card ({u : (ZMod (p ^ e))ˣ | u ^ 2 = 1} : Set _) := rfl
    _ = ({1, -1} : Set (ZMod (p ^ e))ˣ).ncard := by rw [hS, Nat.card_coe_set_eq]
    _ = 2 := Set.ncard_pair hne

lemma natCard_sqOne_zmod (N : ℕ) (hodd : Odd N) :
    Nat.card {u : (ZMod N)ˣ // u ^ 2 = 1} = 2 ^ omega N := by
  induction N using Nat.recOnPosPrimePosCoprime with
  | prime_pow p n hp hn =>
    have hp2 : p ≠ 2 := by
      rintro rfl
      exact (Nat.not_odd_iff_even.mpr (Nat.even_pow.mpr ⟨even_two, hn.ne'⟩)) hodd
    rw [natCard_sqOne_primePow p n hp hp2 hn, omega, Nat.primeFactors_prime_pow hn.ne' hp]
    simp
  | zero => exact absurd hodd (by decide)
  | one =>
    rw [omega, Nat.primeFactors_one, Finset.card_empty, pow_zero]
    exact Nat.card_eq_one_iff_unique.mpr ⟨inferInstance, ⟨⟨1, one_pow 2⟩⟩⟩
  | coprime a b ha hb hab iha ihb =>
    have hoa : Odd a := (Nat.odd_mul.mp hodd).1
    have hob : Odd b := (Nat.odd_mul.mp hodd).2
    rw [natCard_sqOne_congr (ZMod.chineseRemainder hab), natCard_sqOne_prod, iha hoa, ihb hob,
      omega, omega, omega, hab.primeFactors_mul,
      Finset.card_union_of_disjoint hab.disjoint_primeFactors, pow_add]

lemma muTwo_card_eq_natCard (N : ℕ) [NeZero N] :
    (muTwo N).card = Nat.card {u : (ZMod N)ˣ // u ^ 2 = 1} := by
  rw [Nat.card_eq_fintype_card, Fintype.card_subtype]
  rfl


lemma properGcd_iff_ne_pm_one_aux (N : ℕ) [NeZero N] (hodd : Odd N) (h1 : 1 < N)
    (ζ : (ZMod N)ˣ) (hζ : ζ ^ 2 = 1) :
    ProperGcd N ((ζ : ZMod N) - 1) ↔ ζ ≠ 1 ∧ ζ ≠ -1 := by
  have hN : 0 < N := by omega
  constructor
  · rintro ⟨hlo, hhi⟩
    refine ⟨?_, ?_⟩
    · rintro rfl
      simp at hhi
    · rintro rfl
      have hu : IsUnit ((((-1 : (ZMod N)ˣ)) : ZMod N) - 1) := by
        have : (((-1 : (ZMod N)ˣ)) : ZMod N) - 1 = -((2 : ℕ) : ZMod N) := by
          push_cast; ring
        rw [this]
        exact ((ZMod.isUnit_iff_coprime 2 N).mpr (Nat.coprime_two_left.mpr hodd)).neg
      obtain ⟨v, hv⟩ := hu
      have := ZMod.val_coe_unit_coprime v
      rw [hv, Nat.Coprime] at this
      omega
  · rintro ⟨hne1, hne2⟩
    have hgpos : 0 < Nat.gcd ((ζ : ZMod N) - 1).val N := Nat.gcd_pos_of_pos_right _ hN
    have hgle : Nat.gcd ((ζ : ZMod N) - 1).val N ≤ N := Nat.le_of_dvd hN (Nat.gcd_dvd_right _ _)
    have hζ' : (ζ : ZMod N) ^ 2 = 1 := by rw [← Units.val_pow_eq_pow_val, hζ, Units.val_one]
    refine ⟨?_, ?_⟩
    · by_contra hcon
      have hg1 : Nat.Coprime ((ζ : ZMod N) - 1).val N := by
        rw [Nat.Coprime]; omega
      have hu : IsUnit ((ζ : ZMod N) - 1) := by
        have := (ZMod.isUnit_iff_coprime _ N).mpr hg1
        rwa [ZMod.natCast_zmod_val] at this
      apply hne2
      have h0 : ((ζ : ZMod N) - 1) * ((ζ : ZMod N) + 1) = 0 := by
        linear_combination hζ'
      have := hu.mul_right_eq_zero.mp h0
      ext
      simp only [Units.val_neg, Units.val_one]
      linear_combination this
    · by_contra hcon
      have hgN : Nat.gcd ((ζ : ZMod N) - 1).val N = N := by omega
      have hdvd : N ∣ ((ζ : ZMod N) - 1).val := by
        have := Nat.gcd_dvd_left ((ζ : ZMod N) - 1).val N
        rwa [hgN] at this
      have hlt := ZMod.val_lt ((ζ : ZMod N) - 1)
      have h0 : ((ζ : ZMod N) - 1).val = 0 := Nat.eq_zero_of_dvd_of_lt hdvd hlt
      rw [ZMod.val_eq_zero] at h0
      apply hne1
      ext
      simp only [Units.val_one]
      linear_combination h0

/-- MAIN. [pp] Theorem thm:DD, exact form. Let `N > 1` be odd. Fix the squares `R` and the coins
`ξ`, and suppose the dependency `S(R, ξ)` is nonempty and `X² = Y²` holds on the fibre, where
`X = ∏_{i ∈ S} r_i`. Then among the roots `r` with `r_i² = R_i`, exactly the proportion
`1 − 2/|μ₂(ℤ/N)|` has `1 < gcd(X − Y, N) < N`. -/
theorem dixon_descent (N : ℕ) [NeZero N] (hodd : Odd N) (h1 : 1 < N) {n : ℕ} {Ξ : Type*}
    (S : (Fin n → (ZMod N)ˣ) → Ξ → Finset (Fin n))
    (Y : (Fin n → (ZMod N)ˣ) → Ξ → (ZMod N)ˣ)
    (R : Fin n → (ZMod N)ˣ) (ξ : Ξ) (hne : (S R ξ).Nonempty)
    (hsq : ∀ r : Fin n → (ZMod N)ˣ, (∀ i, r i ^ 2 = R i) →
      (∏ i ∈ S R ξ, r i) ^ 2 = (Y R ξ) ^ 2) :
    ((univ.filter (fun r : Fin n → (ZMod N)ˣ => (∀ i, r i ^ 2 = R i) ∧
        ProperGcd N (((∏ i ∈ S R ξ, r i : (ZMod N)ˣ) : ZMod N) - ((Y R ξ : (ZMod N)ˣ) : ZMod N)))).card
        : ℝ) =
      (1 - 2 / ((muTwo N).card : ℝ)) *
        ((univ.filter (fun r : Fin n → (ZMod N)ˣ => ∀ i, r i ^ 2 = R i)).card : ℝ) := by
  set T := S R ξ with hT
  set y := Y R ξ with hy
  set F := univ.filter (fun r : Fin n → (ZMod N)ˣ => ∀ i, r i ^ 2 = R i) with hF
  have hmemF : ∀ r, r ∈ F ↔ ∀ i, r i ^ 2 = R i := fun r => by simp [hF]
  have hN3 : 2 < N := by obtain ⟨k, rfl⟩ := hodd; omega
  have : Fact (2 < N) := ⟨hN3⟩
  have hne1 : (1 : (ZMod N)ˣ) ≠ -1 := by
    intro h
    have := congrArg Units.val h
    simp only [Units.val_one, Units.val_neg] at this
    exact ZMod.neg_one_ne_one this.symm
  set ζ : (Fin n → (ZMod N)ˣ) → (ZMod N)ˣ := fun r => (∏ i ∈ T, r i) * y⁻¹ with hζdef
  have hζsq : ∀ r ∈ F, ζ r ^ 2 = 1 := by
    intro r hr
    have := hsq r ((hmemF r).mp hr)
    simp only [ζ, mul_pow, this, inv_pow, mul_inv_cancel]
  have hmaps : Set.MapsTo ζ F (muTwo N) := fun r hr => by
    simp [muTwo, hζsq r hr]
  have hprop : ∀ r, ProperGcd N (((∏ i ∈ T, r i : (ZMod N)ˣ) : ZMod N) - (y : ZMod N)) ↔
      ProperGcd N ((ζ r : ZMod N) - 1) := by
    intro r
    rw [← properGcd_unit_mul N y ((ζ r : ZMod N) - 1)]
    have e : (y : ZMod N) * ((ζ r : ZMod N) - 1) =
        ((∏ i ∈ T, r i : (ZMod N)ˣ) : ZMod N) - (y : ZMod N) := by
      have : y * ζ r = ∏ i ∈ T, r i := by
        simp only [ζ]; rw [mul_comm (∏ i ∈ T, r i), mul_inv_cancel_left]
      rw [mul_sub, mul_one, ← Units.val_mul, this]
    rw [e]
  obtain ⟨i₀, hi₀⟩ := hne
  have hfib : ∀ z ∈ muTwo N,
      (F.filter (fun r => ζ r = z)).card = (F.filter (fun r => ζ r = 1)).card := by
    intro z hz
    have hz2 : z ^ 2 = 1 := (mem_filter.mp hz).2
    set ε : Fin n → (ZMod N)ˣ := Pi.mulSingle i₀ z with hε
    have hεε : ε * ε = 1 := by
      rw [hε, ← Pi.mulSingle_mul, ← sq, hz2, Pi.mulSingle_one]
    have hεsq : ∀ i, (ε i) ^ 2 = 1 := by
      intro i
      by_cases h : i = i₀
      · subst h; simp [ε, hz2]
      · simp [ε, h]
    have hζε : ∀ r, ζ (ε * r) = z * ζ r := by
      intro r
      simp only [ζ, Pi.mul_apply, prod_mul_distrib, ε, prod_pi_mulSingle', hi₀, ↓reduceIte]
      group
    have hεF : ∀ r ∈ F, ε * r ∈ F := by
      intro r hr
      rw [hmemF] at hr ⊢
      intro i
      rw [Pi.mul_apply, mul_pow, hεsq, one_mul, hr]
    apply card_nbij' (fun r => ε * r) (fun r => ε * r)
    · intro r hr
      rw [mem_coe, mem_filter] at hr ⊢
      refine ⟨hεF r hr.1, ?_⟩
      rw [hζε, hr.2, ← sq, hz2]
    · intro r hr
      rw [mem_coe, mem_filter] at hr ⊢
      refine ⟨hεF r hr.1, ?_⟩
      rw [hζε, hr.2, mul_one]
    · intro r _
      simp only [← mul_assoc, hεε, one_mul]
    · intro r _
      simp only [← mul_assoc, hεε, one_mul]
  set c := (F.filter (fun r => ζ r = 1)).card with hc
  set m := (muTwo N).card with hm
  have hFcard : F.card = m * c := by
    rw [card_eq_sum_card_fiberwise hmaps, sum_congr rfl hfib, sum_const, smul_eq_mul]
  have h1mem : (1 : (ZMod N)ˣ) ∈ muTwo N := by simp [muTwo]
  have hm1mem : (-1 : (ZMod N)ˣ) ∈ (muTwo N).erase 1 := by
    rw [mem_erase]; exact ⟨hne1.symm, by simp [muTwo]⟩
  set t := ((muTwo N).erase 1).erase (-1) with ht
  have htcard : t.card = m - 2 := by
    rw [ht, card_erase_of_mem hm1mem, card_erase_of_mem h1mem]; omega
  have hm2 : 2 ≤ m := by
    have := card_le_card (show ({1, -1} : Finset (ZMod N)ˣ) ⊆ muTwo N by
      intro x hx
      simp only [mem_insert, mem_singleton] at hx
      rcases hx with rfl | rfl
      · exact h1mem
      · exact mem_of_mem_erase hm1mem)
    rwa [card_pair hne1] at this
  set G := F.filter (fun r => ζ r ≠ 1 ∧ ζ r ≠ -1) with hG
  have hGeq : univ.filter (fun r : Fin n → (ZMod N)ˣ => (∀ i, r i ^ 2 = R i) ∧
      ProperGcd N (((∏ i ∈ T, r i : (ZMod N)ˣ) : ZMod N) - (y : ZMod N))) = G := by
    ext r
    rw [mem_filter, hG, mem_filter, hmemF]
    simp only [mem_univ, true_and]
    constructor
    · rintro ⟨hr, hp⟩
      refine ⟨hr, ?_⟩
      rwa [hprop, properGcd_iff_ne_pm_one_aux N hodd h1 _ (hζsq r ((hmemF r).mpr hr))] at hp
    · rintro ⟨hr, hp⟩
      refine ⟨hr, ?_⟩
      rwa [hprop, properGcd_iff_ne_pm_one_aux N hodd h1 _ (hζsq r ((hmemF r).mpr hr))]
  have hmapsG : Set.MapsTo ζ G t := by
    intro r hr
    rw [mem_coe, hG, mem_filter] at hr
    rw [mem_coe, ht, mem_erase, mem_erase]
    exact ⟨hr.2.2, hr.2.1, hmaps (mem_coe.mpr hr.1)⟩
  have hGcard : G.card = (m - 2) * c := by
    rw [card_eq_sum_card_fiberwise hmapsG, ← htcard, ← smul_eq_mul, ← sum_const]
    apply sum_congr rfl
    intro z hz
    rw [ht, mem_erase, mem_erase] at hz
    rw [← hfib z hz.2.2, hG, filter_filter]
    congr 1
    apply filter_congr
    intro r _
    constructor
    · exact fun h => h.2
    · intro h; exact ⟨⟨h ▸ hz.2.1, h ▸ hz.1⟩, h⟩
  rw [hGeq, hGcard, hFcard]
  have hmpos : (0 : ℝ) < m := by exact_mod_cast (show 0 < m by omega)
  rw [Nat.cast_mul, Nat.cast_mul, Nat.cast_sub hm2]
  field_simp
  push_cast
  ring

/-- MAIN. `|μ₂(ℤ/N)| = 2^{ω(N)}` for odd `N ≥ 1` ([pp] §sec:notation). -/
theorem card_muTwo (N : ℕ) [NeZero N] (hodd : Odd N) :
    (muTwo N).card = 2 ^ omega N := by
  rw [muTwo_card_eq_natCard, natCard_sqOne_zmod N hodd]

/-- MAIN. Consequently the fruitfulness probability of Theorem DD is `1 − 2^{1−ω(N)} ≥ 1/2`
when `N` is odd with `ω(N) ≥ 2`. -/
theorem dd_prob_ge_half (N : ℕ) [NeZero N] (hodd : Odd N) (hω : 2 ≤ omega N) :
    1 - 2 / ((muTwo N).card : ℝ) = 1 - 2 ^ (1 - (omega N : ℤ)) ∧
      (1 : ℝ) / 2 ≤ 1 - 2 / ((muTwo N).card : ℝ) := by
  rw [card_muTwo N hodd]
  push_cast
  have hpos : (0 : ℝ) < 2 ^ omega N := by positivity
  constructor
  · rw [zpow_sub₀ two_ne_zero, zpow_one, zpow_natCast]
  · have h4 : (4 : ℝ) ≤ 2 ^ omega N := by
      calc (4 : ℝ) = 2 ^ 2 := by norm_num
        _ ≤ 2 ^ omega N := pow_le_pow_right₀ (by norm_num) hω
    have : 2 / (2 ^ omega N : ℝ) ≤ 1 / 2 := by
      rw [div_le_div_iff₀ hpos two_pos]; linarith
    linarith

/-- MAIN. For `ζ ∈ μ₂(ℤ/N)` with `N` odd: `gcd(ζ − 1, N)` is proper iff `ζ ≠ ±1`
(the last step of the proof of Theorem DD). -/
theorem properGcd_iff_ne_pm_one (N : ℕ) [NeZero N] (hodd : Odd N) (h1 : 1 < N)
    (ζ : (ZMod N)ˣ) (hζ : ζ ^ 2 = 1) :
    ProperGcd N ((ζ : ZMod N) - 1) ↔ ζ ≠ 1 ∧ ζ ≠ -1 := by
  exact properGcd_iff_ne_pm_one_aux N hodd h1 ζ hζ

end FactoringA
