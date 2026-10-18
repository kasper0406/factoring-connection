/-
  Theorem D′, notes/dcltheorem_results.md §3.
  (h) Lemma 10: the CRT count (FULL is rare), core counting statements.
-/
import Mathlib

namespace FactoringLean

open Finset

/-- **Counting in one residue class.** A set of naturals in `[a, a + len)` whose elements are
pairwise congruent modulo `M` has at most `len / M + 1` elements. -/
theorem card_le_of_pairwise_modEq (S : Finset ℕ) (a len M : ℕ)
    (hS : S ⊆ Ico a (a + len)) (hcong : ∀ x ∈ S, ∀ y ∈ S, x ≡ y [MOD M]) :
    S.card ≤ len / M + 1 := by
  have hmaps : ∀ x ∈ S, (x - a) / M ∈ range (len / M + 1) := by
    intro x hx
    have hx' := mem_Ico.1 (hS hx)
    rw [mem_range, Nat.lt_succ_iff]
    exact Nat.div_le_div_right (by omega)
  have hinj : Set.InjOn (fun x => (x - a) / M) S := by
    intro x hx y hy hxy
    simp only at hxy
    have hx' := mem_Ico.1 (hS hx)
    have hy' := mem_Ico.1 (hS hy)
    have hc : (x - a) ≡ (y - a) [MOD M] := by
      apply Nat.ModEq.add_left_cancel' a
      rw [Nat.add_sub_cancel' hx'.1, Nat.add_sub_cancel' hy'.1]
      exact hcong x hx y hy
    have e1 := Nat.div_add_mod (x - a) M
    have e2 := Nat.div_add_mod (y - a) M
    unfold Nat.ModEq at hc
    rw [hxy, hc] at e1
    omega
  calc S.card ≤ (range (len / M + 1)).card := card_le_card_of_injOn _ hmaps hinj
    _ = len / M + 1 := card_range _

/-- The set of `p` solving the two congruences of Lemma 10 for one triple `(v, ℓ, δ)`:
`ℓ ∣ B p − A q` (coverage) and `d ∣ (A − B)(p − 1) − ℓ² δ` (FULL in a clean band). -/
def lemma10Set (a len ℓ d : ℕ) (A B q δ : ℤ) : Finset ℕ :=
  (Ico a (a + len)).filter
    (fun p : ℕ => (ℓ : ℤ) ∣ B * p - A * q ∧ (d : ℤ) ∣ (A - B) * (p - 1) - ℓ ^ 2 * δ)

/-- The CRT modulus of Lemma 10: `ℓ · d / gcd(d, A − B)`. -/
def lemma10Modulus (ℓ d : ℕ) (A B : ℤ) : ℕ := ℓ * (d / Int.gcd d (A - B))

/-- Any two solutions of the Lemma 10 congruences are congruent modulo `ℓ · d / gcd(d, A − B)`
(CRT: `ℓ ∤ d` since the band is clean, `ℓ ∤ B`). -/
theorem lemma10_pairwise_modEq (a len ℓ d : ℕ) (A B q δ : ℤ) (hℓ : ℓ.Prime)
    (hB : ¬ (ℓ : ℤ) ∣ B) (hd : 0 < d) (hℓd : ¬ ℓ ∣ d) :
    ∀ x ∈ lemma10Set a len ℓ d A B q δ, ∀ y ∈ lemma10Set a len ℓ d A B q δ,
      x ≡ y [MOD lemma10Modulus ℓ d A B] := by
  intro x hx y hy
  obtain ⟨-, hx1, hx2⟩ := mem_filter.1 hx
  obtain ⟨-, hy1, hy2⟩ := mem_filter.1 hy
  set g := Int.gcd d (A - B)
  -- modulo ℓ
  have hℓxy : (ℓ : ℤ) ∣ (x : ℤ) - y := by
    have : (ℓ : ℤ) ∣ B * ((x : ℤ) - y) := by
      have e : B * ((x : ℤ) - y) = (B * x - A * q) - (B * y - A * q) := by ring
      rw [e]; exact dvd_sub hx1 hy1
    rcases (Nat.prime_iff_prime_int.mp hℓ).dvd_or_dvd this with h | h
    · exact absurd h hB
    · exact h
  have hmodℓ : x ≡ y [MOD ℓ] := by
    rw [← Int.natCast_modEq_iff]
    exact (Int.modEq_iff_dvd.2 (by rw [← neg_sub]; exact (dvd_neg).2 hℓxy))
  -- modulo d / g
  have hdxy : (A - B) * x ≡ (A - B) * y [ZMOD d] := by
    rw [Int.modEq_iff_dvd]
    have e : (A - B) * (y : ℤ) - (A - B) * x =
        ((A - B) * (y - 1) - ℓ ^ 2 * δ) - ((A - B) * (x - 1) - ℓ ^ 2 * δ) := by ring
    rw [e]; exact dvd_sub hy2 hx2
  have hdg : (x : ℤ) ≡ y [ZMOD (d : ℤ) / (g : ℤ)] :=
    Int.ModEq.cancel_left_div_gcd (by exact_mod_cast hd) hdxy
  have hmodd : x ≡ y [MOD d / g] := by
    rw [← Int.natCast_modEq_iff]; push_cast; exact hdg
  -- combine (CRT)
  have hcop : Nat.Coprime ℓ (d / g) := by
    have h1 : Nat.Coprime ℓ d := (Nat.Prime.coprime_iff_not_dvd hℓ).2 hℓd
    exact h1.coprime_dvd_right (Nat.div_dvd_of_dvd (Int.gcd_dvd_left (d : ℤ) (A - B) |>
      fun h => by exact_mod_cast h))
  exact (Nat.modEq_and_modEq_iff_modEq_mul hcop).1 ⟨hmodℓ, hmodd⟩

/-- **Lemma 10, single triple.** For a fixed `(A, B, ℓ, δ)` with `ℓ` prime, `ℓ ∤ B`, `ℓ ∤ d`,
`d > 0`, the number of `p ∈ [a, a + len)` with `ℓ ∣ B p − A q` and
`d ∣ (A − B)(p − 1) − ℓ² δ` is at most `len / (ℓ · d / gcd(d, A − B)) + 1`. -/
theorem lemma10_single (a len ℓ d : ℕ) (A B q δ : ℤ) (hℓ : ℓ.Prime)
    (hB : ¬ (ℓ : ℤ) ∣ B) (hd : 0 < d) (hℓd : ¬ ℓ ∣ d) :
    (lemma10Set a len ℓ d A B q δ).card ≤ len / lemma10Modulus ℓ d A B + 1 :=
  card_le_of_pairwise_modEq _ a len _ (filter_subset _ _)
    (lemma10_pairwise_modEq a len ℓ d A B q δ hℓ hB hd hℓd)

/-- The modulus is large: `ℓ · d ≤ (ℓ · d / gcd(d, A − B)) · |A − B|` when `A ≠ B`. With
`|A − B| < H` this is the note's `ℓ d / g′ ≥ ℓ d / H`. -/
theorem lemma10Modulus_ge (ℓ d : ℕ) (A B : ℤ) (hAB : A ≠ B) :
    ℓ * d ≤ lemma10Modulus ℓ d A B * (A - B).natAbs := by
  unfold lemma10Modulus
  set g := Int.gcd d (A - B)
  have hgd : g ∣ d := by exact_mod_cast Int.gcd_dvd_left (d : ℤ) (A - B)
  have hgAB : g ≤ (A - B).natAbs := by
    apply Nat.le_of_dvd (Int.natAbs_pos.2 (sub_ne_zero.2 hAB))
    exact Int.natCast_dvd.mp (Int.gcd_dvd_right _ _)
  calc ℓ * d = ℓ * (d / g) * g := by rw [mul_assoc, Nat.div_mul_cancel hgd]
    _ ≤ ℓ * (d / g) * (A - B).natAbs := Nat.mul_le_mul_left _ hgAB

/-- **Lemma 10 (union bound).** For a finite family of triples `t ∈ I` (weights `A t, B t`,
band prime `ℓ t`, shift `δ t`), each satisfying the hypotheses of `lemma10_single`, the set
`F` of `p ∈ [a, a + len)` that solve the two congruences for *some* triple satisfies
`|F| ≤ Σ_t (len / M_t + 1)`. -/
theorem lemma10_union {ι : Type*} (I : Finset ι) (a len d : ℕ) (q : ℤ)
    (A B δ : ι → ℤ) (ℓ : ι → ℕ)
    (hℓ : ∀ t ∈ I, (ℓ t).Prime) (hB : ∀ t ∈ I, ¬ ((ℓ t : ℕ) : ℤ) ∣ B t) (hd : 0 < d)
    (hℓd : ∀ t ∈ I, ¬ ℓ t ∣ d) :
    ((Ico a (a + len)).filter (fun p : ℕ => ∃ t ∈ I,
        ((ℓ t : ℕ) : ℤ) ∣ B t * p - A t * q ∧
        (d : ℤ) ∣ (A t - B t) * (p - 1) - (ℓ t : ℤ) ^ 2 * δ t)).card
      ≤ ∑ t ∈ I, (len / lemma10Modulus (ℓ t) d (A t) (B t) + 1) := by
  classical
  have hsub : (Ico a (a + len)).filter (fun p : ℕ => ∃ t ∈ I,
        ((ℓ t : ℕ) : ℤ) ∣ B t * p - A t * q ∧
        (d : ℤ) ∣ (A t - B t) * (p - 1) - (ℓ t : ℤ) ^ 2 * δ t)
      ⊆ I.biUnion (fun t => lemma10Set a len (ℓ t) d (A t) (B t) q (δ t)) := by
    intro p hp
    obtain ⟨hpI, t, ht, h1, h2⟩ := mem_filter.1 hp
    exact mem_biUnion.2 ⟨t, ht, mem_filter.2 ⟨hpI, h1, h2⟩⟩
  calc _ ≤ _ := card_le_card hsub
    _ ≤ ∑ t ∈ I, (lemma10Set a len (ℓ t) d (A t) (B t) q (δ t)).card := card_biUnion_le
    _ ≤ _ := sum_le_sum fun t ht =>
        lemma10_single a len (ℓ t) d (A t) (B t) q (δ t) (hℓ t ht) (hB t ht) hd (hℓd t ht)

end FactoringLean
