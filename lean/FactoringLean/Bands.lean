/-
  Theorem D′, notes/dcltheorem_results.md §3.
  Lemma 5 (one-sided extraction), (e) Lemma 6 (four-band pigeonhole),
  (f) Lemma 7 (clean-band FULL criterion).
-/
import Mathlib

namespace FactoringLean

open Finset

/-! ### Lemma 5: one-sided extraction -/

/-- **Lemma 5 (one-sided extraction).** Let `f(X) = ∏_{i<m} (X − γ^i)`. If the giant `G`
collides with a baby modulo `p` (`p ∣ G − γ^{i*}`, `i* < m`) and with no baby modulo `q`
(`G` not FULL), then `gcd(f(G), p q) = p`. -/
theorem lemma5_one_sided (p q : ℕ) (hp : p.Prime) (hq : q.Prime) (hpq : p ≠ q)
    (G γ : ℤ) (m istar : ℕ) (hi : istar < m)
    (hcol : (p : ℤ) ∣ G - γ ^ istar)
    (hnotfull : ∀ i < m, ¬ (q : ℤ) ∣ G - γ ^ i) :
    Int.gcd (∏ i ∈ range m, (G - γ ^ i)) ((p * q : ℕ) : ℤ) = p := by
  set F := ∏ i ∈ range m, (G - γ ^ i)
  have hpF : (p : ℤ) ∣ F := hcol.trans (dvd_prod_of_mem _ (mem_range.2 hi))
  have hqF : ¬ (q : ℤ) ∣ F := by
    intro h
    obtain ⟨i, hi, hdi⟩ := ((Nat.prime_iff_prime_int.mp hq).dvd_finsetProd_iff _).1 h
    exact hnotfull i (mem_range.1 hi) hdi
  have hcop : Nat.Coprime p q := (Nat.coprime_primes hp hq).2 hpq
  show Nat.gcd F.natAbs (p * q) = p
  rw [Nat.Coprime.gcd_mul _ hcop]
  have h1 : Nat.gcd F.natAbs p = p := Nat.gcd_eq_right (Int.natCast_dvd.mp hpF)
  have h2 : Nat.gcd F.natAbs q = 1 := by
    have : ¬ q ∣ F.natAbs := fun h => hqF (Int.natCast_dvd.mpr h)
    exact Nat.Coprime.symm ((Nat.Prime.coprime_iff_not_dvd hq).2 this)
  rw [h1, h2, mul_one]

/-- The FULL condition in `ZMod q` is the same as the divisibility condition in Lemma 5. -/
theorem not_full_iff (q : ℕ) (G γ : ℤ) (m : ℕ) :
    (∀ i < m, (G : ZMod q) ≠ (γ : ZMod q) ^ i) ↔ (∀ i < m, ¬ (q : ℤ) ∣ G - γ ^ i) := by
  refine forall₂_congr fun i _ => ?_
  rw [← ZMod.intCast_zmod_eq_zero_iff_dvd, Int.cast_sub, Int.cast_pow, sub_eq_zero]

/-! ### (e) Lemma 6: the four-band pigeonhole -/

/-- **Lemma 6 (clean band).** Let `S₀, …, S₃` be pairwise disjoint sets of primes, all
`> Y`. If `0 < n < (Y+1)^4` (in the note `n = q − 1 < q < Y⁴`), then some band `S_k` contains
no prime dividing `n`. -/
theorem lemma6_clean_band (S : Fin 4 → Finset ℕ)
    (hdisj : ∀ i j, i ≠ j → Disjoint (S i) (S j))
    (Y n : ℕ) (hprime : ∀ i, ∀ ℓ ∈ S i, ℓ.Prime) (hbig : ∀ i, ∀ ℓ ∈ S i, Y < ℓ)
    (hn : 0 < n) (hY : n < (Y + 1) ^ 4) :
    ∃ k, ∀ ℓ ∈ S k, ¬ ℓ ∣ n := by
  by_contra h
  push Not at h
  choose f hf hfd using h
  have hinj : Function.Injective f := by
    intro i j hij
    by_contra hne
    exact Finset.disjoint_left.mp (hdisj i j hne) (hf i) (hij ▸ hf j)
  set s := (univ : Finset (Fin 4)).image f
  have hcard : s.card = 4 := by
    rw [card_image_of_injective _ hinj]; simp
  have hdvd : (∏ a ∈ s, a) ∣ n := by
    apply Finset.prod_primes_dvd
    · intro a ha
      obtain ⟨i, -, rfl⟩ := mem_image.1 ha
      exact Nat.prime_iff.mp (hprime i _ (hf i))
    · intro a ha
      obtain ⟨i, -, rfl⟩ := mem_image.1 ha
      exact hfd i
  have hge : (Y + 1) ^ 4 ≤ ∏ a ∈ s, a := by
    rw [← hcard]
    apply Finset.pow_card_le_prod
    intro a ha
    obtain ⟨i, -, rfl⟩ := mem_image.1 ha
    exact hbig i _ (hf i)
  have := Nat.le_of_dvd hn hdvd
  omega

/-- Lemma 6 in the form used: for a prime `q ≤ Y⁴`, some band is clean for `q`, i.e. no
band prime divides `q − 1`. -/
theorem lemma6_clean_band_prime (S : Fin 4 → Finset ℕ)
    (hdisj : ∀ i j, i ≠ j → Disjoint (S i) (S j))
    (Y q : ℕ) (hprime : ∀ i, ∀ ℓ ∈ S i, ℓ.Prime) (hbig : ∀ i, ∀ ℓ ∈ S i, Y < ℓ)
    (hq : q.Prime) (hY : q ≤ Y ^ 4) :
    ∃ k, ∀ ℓ ∈ S k, ¬ ℓ ∣ q - 1 := by
  apply lemma6_clean_band S hdisj Y (q - 1) hprime hbig
  · have := hq.two_le; omega
  · have : Y ^ 4 < (Y + 1) ^ 4 := Nat.pow_lt_pow_left (Nat.lt_succ_self Y) (by norm_num)
    omega

/-! ### (f) Lemma 7: the clean-band FULL criterion -/

/-- **(f) core.** If `gcd(e, d) = 1` then `d ∣ e x ↔ d ∣ x`. -/
theorem coprime_dvd_mul_iff (e d x : ℕ) (h : Nat.Coprime e d) : d ∣ e * x ↔ d ∣ x :=
  ⟨fun hd => (Nat.Coprime.symm h).dvd_of_dvd_mul_left hd, fun hd => dvd_mul_of_dvd_right hd e⟩

/-- **(f) core, integer version.** -/
theorem coprime_dvd_mul_iff_int (e d x : ℤ) (h : IsCoprime e d) : d ∣ e * x ↔ d ∣ x :=
  ⟨fun hd => h.symm.dvd_of_dvd_mul_left hd, fun hd => dvd_mul_of_dvd_right hd e⟩

/-- Cleanliness gives the coprimality used in Lemma 7: if no band prime divides `q − 1`, then
`L = ∏ ℓ²` is coprime to `q − 1`. -/
theorem clean_coprime (S : Finset ℕ) (hprime : ∀ ℓ ∈ S, ℓ.Prime) (n : ℕ)
    (hclean : ∀ ℓ ∈ S, ¬ ℓ ∣ n) : Nat.Coprime (∏ ℓ ∈ S, ℓ ^ 2) n := by
  apply Nat.Coprime.prod_left
  intro ℓ hℓ
  exact Nat.Coprime.pow_left 2 ((Nat.Prime.coprime_iff_not_dvd (hprime ℓ hℓ)).2 (hclean ℓ hℓ))

/-- `gcd(L, q − 1) = 1`, `M ∣ L` and `d ∣ q − 1` give `gcd(M, d) = 1`
(with `M = L/ℓ²` and `d = ord_q(2)`). -/
theorem coprime_of_clean (L M d n : ℕ) (hL : Nat.Coprime L n) (hM : M ∣ L) (hd : d ∣ n) :
    Nat.Coprime M d :=
  (hL.coprime_dvd_left hM).coprime_dvd_right hd

/-- **Lemma 7, group core.** In any group, with `u = x^M`, `γ = x^(M ℓ²)` and
`gcd(M, ord x) = 1`: `u^D = γ^δ ↔ ord x ∣ D − ℓ² δ`. -/
theorem lemma7_core {G : Type*} [Group G] (x : G) (M ℓ D δ : ℤ)
    (hM : IsCoprime M (orderOf x : ℤ)) :
    (x ^ M) ^ D = (x ^ (M * ℓ ^ 2)) ^ δ ↔ (orderOf x : ℤ) ∣ D - ℓ ^ 2 * δ := by
  rw [← zpow_mul, ← zpow_mul, ← div_eq_one, div_eq_mul_inv, ← zpow_sub,
    ← orderOf_dvd_iff_zpow_eq_one]
  have e : M * D - M * ℓ ^ 2 * δ = M * (D - ℓ ^ 2 * δ) := by ring
  rw [e]
  exact coprime_dvd_mul_iff_int _ _ _ hM

/-- **Lemma 7 (FULL in a clean block).** Let the true giant be `γ^{i*} · u^D` mod `q`
(Lemma 4, with `D = (A−B)(p−1)`) and `0 ≤ i* < m`. It is FULL (equals some baby `γ^{i'}`,
`0 ≤ i' < m`) iff `ord ∣ D − ℓ² δ` for some `δ ∈ [−i*, m − 1 − i*]`. -/
theorem lemma7_full_iff {G : Type*} [Group G] (x : G) (M ℓ D istar m : ℤ)
    (hM : IsCoprime M (orderOf x : ℤ)) :
    (∃ i' : ℤ, 0 ≤ i' ∧ i' < m ∧
        (x ^ (M * ℓ ^ 2)) ^ istar * (x ^ M) ^ D = (x ^ (M * ℓ ^ 2)) ^ i') ↔
      (∃ δ : ℤ, -istar ≤ δ ∧ δ ≤ m - 1 - istar ∧ (orderOf x : ℤ) ∣ D - ℓ ^ 2 * δ) := by
  set γ := x ^ (M * ℓ ^ 2)
  have key : ∀ i' : ℤ, γ ^ istar * (x ^ M) ^ D = γ ^ i' ↔ (x ^ M) ^ D = γ ^ (i' - istar) := by
    intro i'
    have hsplit : γ ^ i' = γ ^ istar * γ ^ (i' - istar) := by
      rw [← zpow_add]; congr 1; ring
    rw [hsplit, mul_right_inj]
  constructor
  · rintro ⟨i', h0, hm, h⟩
    refine ⟨i' - istar, by omega, by omega, ?_⟩
    exact (lemma7_core x M ℓ D _ hM).1 ((key i').1 h)
  · rintro ⟨δ, h0, hm, h⟩
    refine ⟨δ + istar, by omega, by omega, ?_⟩
    rw [key]
    have : δ + istar - istar = δ := by ring
    rw [this]
    exact (lemma7_core x M ℓ D δ hM).2 h

/-- Lemma 7 specialised to `(ZMod q)ˣ`: with a clean band (`gcd(L, q − 1) = 1`) and
`M ∣ L`, the coprimality hypothesis of `lemma7_full_iff` holds for every unit `x`. -/
theorem lemma7_coprime_zmod (q : ℕ) [Fact q.Prime] (x : (ZMod q)ˣ) (L M : ℕ)
    (hL : Nat.Coprime L (q - 1)) (hM : M ∣ L) :
    IsCoprime (M : ℤ) (orderOf x : ℤ) := by
  rw [Int.isCoprime_iff_gcd_eq_one, Int.gcd_natCast_natCast]
  exact coprime_of_clean L M (orderOf x) (q - 1) hL hM (ZMod.orderOf_units_dvd_card_sub_one x)

end FactoringLean
