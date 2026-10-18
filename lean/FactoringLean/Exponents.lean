/-
  Theorem D′, notes/dcltheorem_results.md §3.
  (a) The exponent identity and its Fermat consequence.
  Lemma 3 (exact collision mod p) and Lemma 4 (discrepancy mod q).
-/
import Mathlib

namespace FactoringLean

/-! ### (a) The exponent identity -/

/-- **Exponent identity.** `B + A N − (B p + A q) = (p − 1)(A q − B)` when `N = p q`. -/
theorem exponent_identity (A B p q : ℤ) :
    B + A * (p * q) - (B * p + A * q) = (p - 1) * (A * q - B) := by
  ring

/-- **Fermat consequence of the exponent identity.** In `ZMod p`, `a^(B + A N) = a^(B p + A q)`
for `N = p q`. This holds for every `a` (no hypothesis `p ∤ a` is needed), since
`x^p = x` in `ZMod p`. -/
theorem pow_exponent_identity_zmod (p : ℕ) [Fact p.Prime] (a : ZMod p) (A B q : ℕ) :
    a ^ (B + A * (p * q)) = a ^ (B * p + A * q) := by
  have hp : ∀ k : ℕ, a ^ (p * k) = a ^ k := by
    intro k; rw [pow_mul, ZMod.pow_card]
  have e1 : A * (p * q) = p * (A * q) := by ring
  rw [pow_add, pow_add, e1, hp, mul_comm B p, hp]

/-- Natural-number congruence form of (a): `a^(B + A N) ≡ a^(B p + A q) (mod p)`. -/
theorem pow_exponent_identity_modEq (p : ℕ) (hp : p.Prime) (a A B q : ℕ) :
    a ^ (B + A * (p * q)) ≡ a ^ (B * p + A * q) [MOD p] := by
  have := Fact.mk hp
  rw [← ZMod.natCast_eq_natCast_iff]
  push_cast
  exact pow_exponent_identity_zmod p (a : ZMod p) A B q

/-- The same statement with the hypothesis `p ∤ a` of the note (Fermat), stated over `ℤ` with
integer exponent difference: if `p ∤ a`, then for integer exponents `e₁ e₂` with
`(p − 1) ∣ e₁ − e₂` the unit `a` satisfies `a^e₁ = a^e₂` in `(ZMod p)ˣ`. -/
theorem units_zpow_eq_of_sub_dvd (p : ℕ) [hp : Fact p.Prime] (u : (ZMod p)ˣ) (e₁ e₂ : ℤ)
    (h : ((p : ℤ) - 1) ∣ e₁ - e₂) : u ^ e₁ = u ^ e₂ := by
  obtain ⟨k, hk⟩ := h
  have h1 : u ^ ((p : ℤ) - 1) = 1 := by
    have hp1 : (1 : ℕ) ≤ p := hp.out.one_lt.le
    have : ((p : ℤ) - 1) = ((p - 1 : ℕ) : ℤ) := by push_cast [hp1]; ring
    rw [this, zpow_natCast, ZMod.units_pow_card_sub_one_eq_one]
  have : e₁ = e₂ + ((p : ℤ) - 1) * k := by linarith
  rw [this, zpow_add, zpow_mul, h1, one_zpow, mul_one]

/-! ### Lemma 3 (exact collision mod p)

Notation of the note: the base `x` (= 2), `M = L/ℓ²`, `u = x^M`, `γ = x^L = u^{ℓ²}`.
The class lock (Lemma 1) gives an integer `g` with `ℓ² g = B p + A q − c`. -/

/-- **Lemma 3, core.** If `ℓ² g = B p + A q − c` then `u^(B + A N − c) = γ^g` modulo `p`,
where `u = x^M` and `γ = x^(M ℓ²)`, for any unit `x` mod `p`. No order hypothesis. -/
theorem lemma3_core (p : ℕ) [Fact p.Prime] (x : (ZMod p)ˣ) (M A B q c ℓ g : ℤ)
    (hg : ℓ ^ 2 * g = B * p + A * q - c) :
    (x ^ M) ^ (B + A * (p * q) - c) = (x ^ (M * ℓ ^ 2)) ^ g := by
  rw [← zpow_mul, ← zpow_mul]
  apply units_zpow_eq_of_sub_dvd
  refine ⟨M * (A * q - B), ?_⟩
  have : M * ℓ ^ 2 * g = M * (B * p + A * q - c) := by rw [← hg]; ring
  rw [this]; ring

/-- **Lemma 3.** With `i* = g − g_lo − j m`, the true giant
`G = u^(B + A N − c) · γ^(−(g_lo + j m))` equals the baby `γ^(i*)` modulo `p`. -/
theorem lemma3 (p : ℕ) [Fact p.Prime] (x : (ZMod p)ˣ) (M A B q c ℓ g glo j m : ℤ)
    (hg : ℓ ^ 2 * g = B * p + A * q - c) :
    (x ^ M) ^ (B + A * (p * q) - c) * (x ^ (M * ℓ ^ 2)) ^ (-(glo + j * m))
      = (x ^ (M * ℓ ^ 2)) ^ (g - glo - j * m) := by
  rw [lemma3_core p x M A B q c ℓ g hg, ← zpow_add]
  congr 1; ring

/-! ### Lemma 4 (discrepancy mod q) -/

/-- **Lemma 4, core.** If `ℓ² g = B p + A q − c`, then modulo the prime `q`,
`u^(B + A N − c) = γ^g · u^((A − B)(p − 1))`. -/
theorem lemma4_core (q : ℕ) [Fact q.Prime] (x : (ZMod q)ˣ) (M A B p c ℓ g : ℤ)
    (hg : ℓ ^ 2 * g = B * p + A * q - c) :
    (x ^ M) ^ (B + A * (p * q) - c) = (x ^ (M * ℓ ^ 2)) ^ g * (x ^ M) ^ ((A - B) * (p - 1)) := by
  rw [← zpow_mul, ← zpow_mul, ← zpow_mul, ← zpow_add]
  apply units_zpow_eq_of_sub_dvd
  refine ⟨M * A * (p - 1), ?_⟩
  have : M * ℓ ^ 2 * g = M * (B * p + A * q - c) := by rw [← hg]; ring
  rw [this]; ring

/-- The exponent difference in Lemma 4 is `(A − B)(p − 1)` modulo `q − 1`:
`(B + A N − c) − (B p + A q − c) ≡ (A − B)(p − 1) (mod q − 1)`. -/
theorem lemma4_exponent_difference (A B p q c : ℤ) :
    ((q - 1) ∣ (B + A * (p * q) - c) - (B * p + A * q - c) - (A - B) * (p - 1)) :=
  ⟨A * (p - 1), by ring⟩

/-- **Lemma 4.** The true giant equals `γ^(i*) · u^((A − B)(p − 1))` modulo `q`. -/
theorem lemma4 (q : ℕ) [Fact q.Prime] (x : (ZMod q)ˣ) (M A B p c ℓ g glo j m : ℤ)
    (hg : ℓ ^ 2 * g = B * p + A * q - c) :
    (x ^ M) ^ (B + A * (p * q) - c) * (x ^ (M * ℓ ^ 2)) ^ (-(glo + j * m))
      = (x ^ (M * ℓ ^ 2)) ^ (g - glo - j * m) * (x ^ M) ^ ((A - B) * (p - 1)) := by
  rw [lemma4_core q x M A B p c ℓ g hg, mul_right_comm, ← zpow_add]
  congr 2; ring

end FactoringLean
