/-
  Theorem D′, notes/dcltheorem_results.md §3.
  (g) Lemma 8: primes q with small ord_q(2) are few.
-/
import Mathlib

namespace FactoringLean

open Finset

/-- A set `s` of primes, each `≥ Q`, all dividing `n ≠ 0`, satisfies `Q^|s| ≤ n`. -/
theorem pow_card_le_of_primes_dvd (s : Finset ℕ) (Q n : ℕ) (hn : n ≠ 0)
    (hprime : ∀ q ∈ s, q.Prime) (hQ : ∀ q ∈ s, Q ≤ q) (hdvd : ∀ q ∈ s, q ∣ n) :
    Q ^ s.card ≤ n := by
  have h1 : (∏ q ∈ s, q) ∣ n :=
    Finset.prod_primes_dvd n (fun q hq => Nat.prime_iff.mp (hprime q hq)) hdvd
  exact (Finset.pow_card_le_prod s id Q hQ).trans (Nat.le_of_dvd (Nat.pos_of_ne_zero hn) h1)

/-- If `q` is an odd prime and `ord_q(2) ≤ T`, then `q ∣ 2^k − 1` for some `1 ≤ k ≤ T`. -/
theorem dvd_two_pow_sub_one_of_orderOf_le (q T : ℕ) (hq : q.Prime) (hq2 : q ≠ 2)
    (hT : orderOf (2 : ZMod q) ≤ T) : ∃ k, 1 ≤ k ∧ k ≤ T ∧ q ∣ 2 ^ k - 1 := by
  have := Fact.mk hq
  have h2 : (2 : ZMod q) ≠ 0 := by
    intro h
    have : q ∣ 2 := by exact_mod_cast (ZMod.natCast_eq_zero_iff 2 q).1 (by exact_mod_cast h)
    exact hq2 ((Nat.prime_dvd_prime_iff_eq hq Nat.prime_two).1 this)
  have hpow : (2 : ZMod q) ^ (q - 1) = 1 := ZMod.pow_card_sub_one_eq_one h2
  have hpos : 0 < orderOf (2 : ZMod q) :=
    orderOf_pos_iff.2 (isOfFinOrder_iff_pow_eq_one.2 ⟨q - 1, by have := hq.two_le; omega, hpow⟩)
  refine ⟨orderOf (2 : ZMod q), hpos, hT, ?_⟩
  rw [← ZMod.natCast_eq_zero_iff, Nat.cast_sub Nat.one_le_two_pow]
  push_cast
  rw [pow_orderOf_eq_one, sub_self]

/-- Primes dividing `2^k − 1` (`k ≥ 1`) and exceeding `2^b`: at most `(k−1)/b` of them, i.e.
`b · #s < k`. -/
theorem card_mul_lt_of_dvd_mersenne (s : Finset ℕ) (b k : ℕ) (hk : 1 ≤ k)
    (hprime : ∀ q ∈ s, q.Prime) (hbig : ∀ q ∈ s, 2 ^ b ≤ q) (hdvd : ∀ q ∈ s, q ∣ 2 ^ k - 1) :
    b * s.card < k := by
  have hne : 2 ^ k - 1 ≠ 0 := by
    have : 2 ≤ 2 ^ k := by
      calc 2 = 2 ^ 1 := by norm_num
        _ ≤ 2 ^ k := Nat.pow_le_pow_right (by norm_num) hk
    omega
  have h := pow_card_le_of_primes_dvd s (2 ^ b) (2 ^ k - 1) hne hprime hbig hdvd
  rw [← pow_mul] at h
  have h' : 2 ^ (b * s.card) < 2 ^ k := lt_of_le_of_lt h (Nat.sub_lt Nat.one_le_two_pow one_pos)
  exact (Nat.pow_lt_pow_iff_right (by norm_num)).1 h'

/-- Every prime `q` in `s` divides some `2^k − 1` with `1 ≤ k ≤ T` ⇒ `s` is covered by the
prime-divisor sets of `2^k − 1`, `k ∈ [1, T]`. -/
theorem card_le_sum_mersenne (s : Finset ℕ) (T : ℕ)
    (h : ∀ q ∈ s, ∃ k, 1 ≤ k ∧ k ≤ T ∧ q ∣ 2 ^ k - 1) :
    s.card ≤ ∑ k ∈ Icc 1 T, (s.filter (fun q => q ∣ 2 ^ k - 1)).card := by
  have hsub : s ⊆ (Icc 1 T).biUnion (fun k => s.filter (fun q => q ∣ 2 ^ k - 1)) := by
    intro q hq
    obtain ⟨k, hk1, hkT, hd⟩ := h q hq
    exact mem_biUnion.2 ⟨k, mem_Icc.2 ⟨hk1, hkT⟩, mem_filter.2 ⟨hq, hd⟩⟩
  exact (card_le_card hsub).trans card_biUnion_le

/-- `(1 + 2 + ⋯ + T) · 2 = T (T + 1)`. -/
theorem gauss_sum_Icc (T : ℕ) : (∑ k ∈ Icc 1 T, k) * 2 = T * (T + 1) := by
  induction T with
  | zero => simp
  | succ n ih =>
    rw [Finset.sum_Icc_succ_top (by omega), add_mul, ih]; ring

/-- **Lemma 8 (sharp form, as in the note).** Let `s` be a set of primes, each `> 2^b`
(`b ≥ 1`, so `b = ⌊log₂ Q₀⌋` for primes `> Q₀`), each with `ord_q(2) ≤ T`. Then
`b · #s ≤ T(T+1)/2`, i.e. `#s ≤ T(T+1)/(2b)`. -/
theorem lemma8_sharp (s : Finset ℕ) (T b : ℕ) (hb : 1 ≤ b)
    (hs : ∀ q ∈ s, q.Prime ∧ 2 ^ b < q ∧ orderOf (2 : ZMod q) ≤ T) :
    b * s.card * 2 ≤ T * (T + 1) := by
  have hmers : ∀ q ∈ s, ∃ k, 1 ≤ k ∧ k ≤ T ∧ q ∣ 2 ^ k - 1 := by
    intro q hq
    obtain ⟨hp, hbig, hord⟩ := hs q hq
    have hq2 : q ≠ 2 := by
      have : 2 ≤ 2 ^ b := by
        calc 2 = 2 ^ 1 := by norm_num
          _ ≤ 2 ^ b := Nat.pow_le_pow_right (by norm_num) hb
      omega
    exact dvd_two_pow_sub_one_of_orderOf_le q T hp hq2 hord
  have h1 := card_le_sum_mersenne s T hmers
  have h2 : ∀ k ∈ Icc 1 T, b * (s.filter (fun q => q ∣ 2 ^ k - 1)).card ≤ k := by
    intro k hk
    have hk1 := (mem_Icc.1 hk).1
    apply le_of_lt
    apply card_mul_lt_of_dvd_mersenne _ b k hk1
    · intro q hq; exact (hs q (mem_filter.1 hq).1).1
    · intro q hq; exact (hs q (mem_filter.1 hq).1).2.1.le
    · intro q hq; exact (mem_filter.1 hq).2
  have h3 : b * s.card ≤ ∑ k ∈ Icc 1 T, k := by
    calc b * s.card ≤ b * ∑ k ∈ Icc 1 T, (s.filter (fun q => q ∣ 2 ^ k - 1)).card :=
          Nat.mul_le_mul_left b h1
      _ = ∑ k ∈ Icc 1 T, b * (s.filter (fun q => q ∣ 2 ^ k - 1)).card := by rw [mul_sum]
      _ ≤ ∑ k ∈ Icc 1 T, k := sum_le_sum h2
  have h4 : (∑ k ∈ Icc 1 T, k) * 2 = T * (T + 1) := gauss_sum_Icc T
  calc b * s.card * 2 ≤ (∑ k ∈ Icc 1 T, k) * 2 := Nat.mul_le_mul_right 2 h3
    _ = T * (T + 1) := h4

/-- **Lemma 8 (form requested): the number of odd primes `q` with `ord_q(2) ≤ T` is at most
`T²`.** -/
theorem lemma8 (s : Finset ℕ) (T : ℕ)
    (hs : ∀ q ∈ s, q.Prime ∧ q ≠ 2 ∧ orderOf (2 : ZMod q) ≤ T) :
    s.card ≤ T ^ 2 := by
  have hmers : ∀ q ∈ s, ∃ k, 1 ≤ k ∧ k ≤ T ∧ q ∣ 2 ^ k - 1 := fun q hq =>
    dvd_two_pow_sub_one_of_orderOf_le q T (hs q hq).1 (hs q hq).2.1 (hs q hq).2.2
  have h1 := card_le_sum_mersenne s T hmers
  have h2 : ∀ k ∈ Icc 1 T, (s.filter (fun q => q ∣ 2 ^ k - 1)).card ≤ T := by
    intro k hk
    have := card_mul_lt_of_dvd_mersenne (s.filter (fun q => q ∣ 2 ^ k - 1)) 1 k
      (mem_Icc.1 hk).1 (fun q hq => (hs q (mem_filter.1 hq).1).1)
      (fun q hq => (hs q (mem_filter.1 hq).1).1.two_le) (fun q hq => (mem_filter.1 hq).2)
    have := (mem_Icc.1 hk).2
    omega
  calc s.card ≤ ∑ k ∈ Icc 1 T, (s.filter (fun q => q ∣ 2 ^ k - 1)).card := h1
    _ ≤ (Icc 1 T).card • T := sum_le_card_nsmul _ _ _ h2
    _ = T ^ 2 := by simp [sq]

end FactoringLean
