/-
  Theorem D′, notes/dcltheorem_results.md §3.
  (b) The coverage criterion.
  (c) Lemma 1: the modulus-squared class lock (HH 2022 Prop. 3.3).
  Distinctness of the excluded classes used in Lemma 9.
-/
import Mathlib

namespace FactoringLean

/-! ### (b) Coverage criterion -/

/-- **Coverage criterion.** For a prime `ℓ` with `ℓ ∤ p` (in particular when `ℓ ∤ N = p q`),
`ℓ ∣ B p² − A N  ↔  ℓ ∣ B p − A q`. -/
theorem coverage_criterion (ℓ A B p q : ℤ) (hℓ : Prime ℓ) (hp : ¬ ℓ ∣ p) :
    ℓ ∣ B * p ^ 2 - A * (p * q) ↔ ℓ ∣ B * p - A * q := by
  have e : B * p ^ 2 - A * (p * q) = p * (B * p - A * q) := by ring
  rw [e]
  constructor
  · intro h
    rcases hℓ.dvd_or_dvd h with h | h
    · exact absurd h hp
    · exact h
  · intro h; exact dvd_mul_of_dvd_right h p

/-- Coverage criterion with the hypothesis as in the note: `ℓ ∤ N`. -/
theorem coverage_criterion' (ℓ A B p q : ℤ) (hℓ : Prime ℓ) (hN : ¬ ℓ ∣ p * q) :
    ℓ ∣ B * p ^ 2 - A * (p * q) ↔ ℓ ∣ B * p - A * q :=
  coverage_criterion ℓ A B p q hℓ (fun h => hN (dvd_mul_of_dvd_left h q))

/-- The computed roots: if `p ≡ r (mod ℓ)` then `ℓ ∣ B r² − A N ↔ ℓ ∣ B p² − A N`. So
"`p mod ℓ` is a root of `B X² ≡ A N`" is the same as `ℓ ∣ B p² − A N`. -/
theorem root_iff_of_modEq (ℓ A B N p r : ℤ) (h : ℓ ∣ p - r) :
    ℓ ∣ B * r ^ 2 - A * N ↔ ℓ ∣ B * p ^ 2 - A * N := by
  have e : B * p ^ 2 - A * N = (B * r ^ 2 - A * N) + (p - r) * (B * (p + r)) := by ring
  rw [e]
  exact (dvd_add_left (dvd_mul_of_dvd_left h _)).symm

/-! ### (c) Lemma 1: the modulus-squared class lock -/

/-- **Lemma 1 (class lock).** Let `ℓ` be prime, `N = p q`, `r` with `ℓ ∤ r` and
`B r² ≡ A N (mod ℓ)`, `s ≡ N r⁻¹ (mod ℓ²)` (i.e. `ℓ² ∣ r s − N`) and
`c ≡ B r + A s (mod ℓ²)`. If `p ≡ r (mod ℓ)` then `ℓ² ∣ B p + A q − c`, i.e.
`g = (B p + A q − c)/ℓ²` is an integer. -/
theorem lemma1_class_lock (ℓ A B p q r s c : ℤ) (hℓ : Prime ℓ) (hr : ¬ ℓ ∣ r)
    (hroot : ℓ ∣ B * r ^ 2 - A * (p * q))
    (hs : ℓ ^ 2 ∣ r * s - p * q)
    (hc : ℓ ^ 2 ∣ c - (B * r + A * s))
    (hpr : ℓ ∣ p - r) :
    ℓ ^ 2 ∣ B * p + A * q - c := by
  -- ℓ ∤ p since p ≡ r
  have hp : ¬ ℓ ∣ p := by
    intro h; apply hr
    have : r = p - (p - r) := by ring
    rw [this]; exact dvd_sub h hpr
  -- coverage: ℓ ∣ B p − A q
  have hcov : ℓ ∣ B * p - A * q :=
    (coverage_criterion ℓ A B p q hℓ hp).1 ((root_iff_of_modEq ℓ A B (p * q) p r hpr).1 hroot)
  obtain ⟨x, hx⟩ := hpr
  obtain ⟨y, hy⟩ := hcov
  obtain ⟨t, ht⟩ := hs
  obtain ⟨u, hu⟩ := hc
  -- r · (B p + A q − c) is an explicit multiple of ℓ²
  have key : r * (B * p + A * q - c) = ℓ ^ 2 * (x * (y - x * B) - (A * t + r * u)) := by
    have hp' : p = r + ℓ * x := by linarith
    have hc' : c = B * r + A * s + ℓ ^ 2 * u := by linarith
    have hN : A * (r * q) = A * (r * s) - A * ℓ ^ 2 * t - A * (ℓ * x * q) := by
      have : r * q = p * q - ℓ * x * q := by rw [hp']; ring
      rw [this]; linear_combination (-A) * ht
    have hy' : r * B - A * q = ℓ * y - ℓ * x * B := by
      have : r = p - ℓ * x := by linarith
      rw [this]; linear_combination hy
    calc r * (B * p + A * q - c)
        = r * B * p + A * (r * q) - r * c := by ring
      _ = ℓ * x * (r * B - A * q) - ℓ ^ 2 * (A * t + r * u) := by
          rw [hN, hc', hp']; ring
      _ = ℓ ^ 2 * (x * (y - x * B) - (A * t + r * u)) := by rw [hy']; ring
  have hdvd : ℓ ^ 2 ∣ r * (B * p + A * q - c) := ⟨_, key⟩
  have hcop : IsCoprime (ℓ ^ 2) r :=
    ((Irreducible.coprime_iff_not_dvd hℓ.irreducible).2 hr).pow_left
  exact hcop.dvd_of_dvd_mul_left hdvd

/-! ### Distinct excluded classes (used in Lemma 9)

For a prime `ℓ ∤ q`, `ℓ ∣ B p − A q` defines the class `p ≡ A q B⁻¹ (mod ℓ)`. Two weights
`(A, B)`, `(A', B')` with `0 < A, B, A', B' ≤ H`, `H² < ℓ`, and `gcd(A,B) = gcd(A',B') = 1`
give the same class only if they are equal. -/

/-- If both `ℓ ∣ B p − A q` and `ℓ ∣ B' p − A' q` with `ℓ ∤ q` prime, then
`ℓ ∣ A B' − A' B`. -/
theorem same_class_dvd (ℓ A B A' B' p q : ℤ) (hℓ : Prime ℓ) (hq : ¬ ℓ ∣ q)
    (h1 : ℓ ∣ B * p - A * q) (h2 : ℓ ∣ B' * p - A' * q) : ℓ ∣ A * B' - A' * B := by
  have e : q * (A * B' - A' * B) = B * (B' * p - A' * q) - B' * (B * p - A * q) := by ring
  have : ℓ ∣ q * (A * B' - A' * B) := by
    rw [e]; exact dvd_sub (dvd_mul_of_dvd_right h2 _) (dvd_mul_of_dvd_right h1 _)
  rcases hℓ.dvd_or_dvd this with h | h
  · exact absurd h hq
  · exact h

/-- **Distinct classes.** Primitive weights in `[1, H]²` with `H² < ℓ` that cover the same
class modulo a prime `ℓ ∤ q` are equal. Hence the excluded classes number exactly `|V|`. -/
theorem weights_eq_of_same_class (ℓ H : ℕ) (A B A' B' : ℕ) (p q : ℤ)
    (hℓ : ℓ.Prime) (hH : H ^ 2 < ℓ)
    (hA : 0 < A ∧ A ≤ H) (hB : 0 < B ∧ B ≤ H) (hA' : 0 < A' ∧ A' ≤ H) (hB' : 0 < B' ∧ B' ≤ H)
    (hAB : Nat.Coprime A B) (hAB' : Nat.Coprime A' B')
    (hq : ¬ (ℓ : ℤ) ∣ q)
    (h1 : (ℓ : ℤ) ∣ B * p - A * q) (h2 : (ℓ : ℤ) ∣ B' * p - A' * q) :
    A = A' ∧ B = B' := by
  have hd := same_class_dvd ℓ A B A' B' p q (Nat.prime_iff_prime_int.mp hℓ) hq h1 h2
  -- |A B' − A' B| < H² < ℓ, so A B' = A' B
  have hlt1 : A * B' ≤ H ^ 2 := by nlinarith [hA.2, hB'.2]
  have hlt2 : A' * B ≤ H ^ 2 := by nlinarith [hA'.2, hB.2]
  have heq : A * B' = A' * B := by
    have habs : ((A * B' : ℕ) : ℤ) - ((A' * B : ℕ) : ℤ) = 0 := by
      apply Int.eq_zero_of_abs_lt_dvd (m := (ℓ : ℤ))
      · push_cast; exact hd
      · have h1' : ((A * B' : ℕ) : ℤ) < ℓ := by exact_mod_cast lt_of_le_of_lt hlt1 hH
        have h2' : ((A' * B : ℕ) : ℤ) < ℓ := by exact_mod_cast lt_of_le_of_lt hlt2 hH
        rw [abs_lt]; constructor <;> omega
    omega
  -- primitive: A ∣ A' B and gcd(A,B)=1 ⇒ A ∣ A'; symmetric
  have hA_dvd : A ∣ A' := by
    have : A ∣ A' * B := ⟨B', heq.symm ▸ rfl⟩
    exact hAB.dvd_of_dvd_mul_right this
  have hA'_dvd : A' ∣ A := by
    have : A' ∣ A * B' := ⟨B, heq⟩
    exact hAB'.dvd_of_dvd_mul_right this
  have hAA : A = A' := Nat.dvd_antisymm hA_dvd hA'_dvd
  subst hAA
  refine ⟨rfl, ?_⟩
  exact (Nat.eq_of_mul_eq_mul_left hA.1 heq).symm
end FactoringLean
