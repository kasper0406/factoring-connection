/-
  Theorem D′, notes/dcltheorem_results.md §3–4.
  Lemma 9 (coverage via the arithmetic large sieve) and the exceptional-set bound E_b,
  DEDUCED FROM the large-sieve inequality, which is NOT proved here: it enters only as the
  explicit hypothesis `ArithmeticLargeSievePrimeTerms` (a `Prop`, passed as an argument; no `axiom`).
-/
import Mathlib
import FactoringLean.Coverage

namespace FactoringLean

open Finset

/-- **HYPOTHESIS — NOT PROVED IN THIS DEVELOPMENT.**
The arithmetic large sieve, keeping only the `d = 1` and prime `d = ℓ` terms
(Montgomery 1968; Montgomery–Vaughan 1973, constant `N₀ + Q²`):

let `U` be a finite set of integers in `(M₀, M₀ + N₀]`, `Ps` a finite set of primes `≤ Q`, and
for each `ℓ ∈ Ps` let `Ω ℓ ⊆ [0, ℓ)` be a set of `ω(ℓ) < ℓ` residues that `U` avoids modulo `ℓ`.
Then `|U| · (1 + Σ_{ℓ ∈ Ps} ω(ℓ)/(ℓ − ω(ℓ))) ≤ N₀ + Q²`.

This is the **prime-term weakening** of the arithmetic large sieve (Montgomery–Vaughan 1973;
cf. Iwaniec–Kowalski Thm 7.14 / Montgomery–Vaughan *Multiplicative Number Theory I* Thm 19.13,
with Montgomery's lemma for the local factors): of `J = Σ_{d ≤ Q} μ²(d) ∏_{ℓ ∣ d} ω/(ℓ−ω)` only the
`d = 1` and prime `d = ℓ ∈ 𝒫` terms are kept. All dropped terms are `≥ 0`, so this is implied
by the full theorem. For a single dyadic band (Lemma 9: sieving primes `ℓ ∈ (2^kY, 2^{k+1}Y]`,
`Q = 2^{k+1}Y`), any product of two band primes exceeds `(2^kY)² > Q`, so there are no
composite `d ≤ Q` built from band primes and nothing is lost.

Theorems below that use it take `(hLS : ArithmeticLargeSievePrimeTerms)` as an argument. -/
def ArithmeticLargeSievePrimeTerms : Prop :=
  ∀ (M₀ : ℤ) (N₀ Q : ℕ) (U : Finset ℤ) (Ps : Finset ℕ) (Ω : ℕ → Finset ℤ),
    (∀ u ∈ U, M₀ < u ∧ u ≤ M₀ + N₀) →
    (∀ ℓ ∈ Ps, ℓ.Prime ∧ ℓ ≤ Q) →
    (∀ ℓ ∈ Ps, (Ω ℓ).card < ℓ ∧ ∀ w ∈ Ω ℓ, 0 ≤ w ∧ w < (ℓ : ℤ)) →
    (∀ u ∈ U, ∀ ℓ ∈ Ps, ∀ w ∈ Ω ℓ, ¬ (ℓ : ℤ) ∣ u - w) →
    (U.card : ℝ) * (1 + ∑ ℓ ∈ Ps, ((Ω ℓ).card : ℝ) / ((ℓ : ℝ) - (Ω ℓ).card)) ≤ N₀ + (Q : ℝ) ^ 2

/-- The excluded residue of the weight `(A, B)` modulo `ℓ`: the `r ∈ [0, ℓ)` with
`B r ≡ A q (mod ℓ)`. -/
noncomputable def excludedResidue (ℓ : ℕ) (q : ℤ) (v : ℕ × ℕ) : ℤ :=
  ((((v.1 : ℤ) * q : ℤ) : ZMod ℓ) * ((v.2 : ℕ) : ZMod ℓ)⁻¹).val

theorem excludedResidue_spec (ℓ : ℕ) [Fact ℓ.Prime] (q : ℤ) (v : ℕ × ℕ)
    (hB : ¬ ℓ ∣ v.2) :
    0 ≤ excludedResidue ℓ q v ∧ excludedResidue ℓ q v < ℓ ∧
      (ℓ : ℤ) ∣ (v.2 : ℤ) * excludedResidue ℓ q v - v.1 * q := by
  refine ⟨by unfold excludedResidue; positivity, ?_, ?_⟩
  · unfold excludedResidue; exact_mod_cast ZMod.val_lt _
  · rw [← ZMod.intCast_zmod_eq_zero_iff_dvd]
    unfold excludedResidue
    have hB' : ((v.2 : ℕ) : ZMod ℓ) ≠ 0 := by
      rwa [Ne, ZMod.natCast_eq_zero_iff]
    push_cast
    rw [ZMod.natCast_val, ZMod.cast_id', id, mul_comm (((v.1 : ZMod ℓ)) * (q : ZMod ℓ)),
      ← mul_assoc, mul_inv_cancel₀ hB', one_mul, sub_self]

/-- `V ⊆ [1,H]²` has at most `H²` elements. -/
theorem card_weights_le (V : Finset (ℕ × ℕ)) (H : ℕ)
    (hV : ∀ v ∈ V, 0 < v.1 ∧ v.1 ≤ H ∧ 0 < v.2 ∧ v.2 ≤ H) : V.card ≤ H ^ 2 := by
  have : V ⊆ Icc 1 H ×ˢ Icc 1 H := by
    intro v hv
    obtain ⟨h1, h2, h3, h4⟩ := hV v hv
    simp only [mem_product, mem_Icc]; omega
  calc V.card ≤ (Icc 1 H ×ˢ Icc 1 H).card := card_le_card this
    _ = H ^ 2 := by simp [sq]

/-- **Lemma 9 (from the large sieve).** Fix `q`. Let `V` be primitive weights in `[1, H]²`, `Ps`
primes in `(H², Q]` not dividing `q`, and `U` integers in `(M₀, M₀ + N₀]` that no target covers:
`ℓ ∤ B p − A q` for all `(A, B) ∈ V`, `ℓ ∈ Ps`. Assuming the large sieve,
`|U| · |V| · Σ_{ℓ ∈ Ps} 1/ℓ ≤ N₀ + Q²`. -/
theorem lemma9_of_largeSieve (hLS : ArithmeticLargeSievePrimeTerms) (M₀ : ℤ) (N₀ Q H : ℕ) (q : ℤ)
    (V : Finset (ℕ × ℕ)) (Ps : Finset ℕ)
    (hV : ∀ v ∈ V, 0 < v.1 ∧ v.1 ≤ H ∧ 0 < v.2 ∧ v.2 ≤ H ∧ Nat.Coprime v.1 v.2)
    (hP : ∀ ℓ ∈ Ps, ℓ.Prime ∧ ℓ ≤ Q ∧ H ^ 2 < ℓ ∧ ¬ (ℓ : ℤ) ∣ q)
    (U : Finset ℤ) (hU : ∀ u ∈ U, M₀ < u ∧ u ≤ M₀ + N₀)
    (huncov : ∀ u ∈ U, ∀ v ∈ V, ∀ ℓ ∈ Ps, ¬ (ℓ : ℤ) ∣ (v.2 : ℤ) * u - v.1 * q) :
    (U.card : ℝ) * ((V.card : ℝ) * ∑ ℓ ∈ Ps, (1 : ℝ) / ℓ) ≤ N₀ + (Q : ℝ) ^ 2 := by
  classical
  set Ω : ℕ → Finset ℤ := fun ℓ => V.image (excludedResidue ℓ q)
  have hVcard : V.card ≤ H ^ 2 :=
    card_weights_le V H (fun v hv => ⟨(hV v hv).1, (hV v hv).2.1, (hV v hv).2.2.1,
      (hV v hv).2.2.2.1⟩)
  -- ℓ ∤ B for (A,B) ∈ V, ℓ ∈ Ps
  have hBℓ : ∀ ℓ ∈ Ps, ∀ v ∈ V, ¬ ℓ ∣ v.2 := by
    intro ℓ hℓ v hv h
    have := Nat.le_of_dvd (hV v hv).2.2.1 h
    have h2 : v.2 ≤ H ^ 2 := by
      have := (hV v hv).2.2.2.1
      rcases Nat.eq_zero_or_pos H with h0 | h0
      · omega
      · nlinarith
    have := (hP ℓ hℓ).2.2.1
    omega
  -- the excluded residues are distinct, so |Ω ℓ| = |V|
  have hΩcard : ∀ ℓ ∈ Ps, (Ω ℓ).card = V.card := by
    intro ℓ hℓ
    have := Fact.mk (hP ℓ hℓ).1
    apply card_image_of_injOn
    intro v hv v' hv' hvv
    have h1 := (excludedResidue_spec ℓ q v (hBℓ ℓ hℓ v hv)).2.2
    have h2 := (excludedResidue_spec ℓ q v' (hBℓ ℓ hℓ v' hv')).2.2
    rw [← hvv] at h2
    obtain ⟨hA, hB⟩ := weights_eq_of_same_class ℓ H v.1 v.2 v'.1 v'.2
      (excludedResidue ℓ q v) q (hP ℓ hℓ).1 (hP ℓ hℓ).2.2.1
      ⟨(hV v hv).1, (hV v hv).2.1⟩ ⟨(hV v hv).2.2.1, (hV v hv).2.2.2.1⟩
      ⟨(hV v' hv').1, (hV v' hv').2.1⟩ ⟨(hV v' hv').2.2.1, (hV v' hv').2.2.2.1⟩
      (hV v hv).2.2.2.2 (hV v' hv').2.2.2.2 (hP ℓ hℓ).2.2.2
      (by rw [mul_comm]; exact (by simpa [mul_comm] using h1))
      (by rw [mul_comm]; exact (by simpa [mul_comm] using h2))
    exact Prod.ext hA hB
  have hLS' := hLS M₀ N₀ Q U Ps Ω hU (fun ℓ hℓ => ⟨(hP ℓ hℓ).1, (hP ℓ hℓ).2.1⟩)
    (by
      intro ℓ hℓ
      have := Fact.mk (hP ℓ hℓ).1
      refine ⟨by rw [hΩcard ℓ hℓ]; exact lt_of_le_of_lt hVcard (hP ℓ hℓ).2.2.1, ?_⟩
      intro w hw
      obtain ⟨v, hv, rfl⟩ := mem_image.1 hw
      exact ⟨(excludedResidue_spec ℓ q v (hBℓ ℓ hℓ v hv)).1,
        (excludedResidue_spec ℓ q v (hBℓ ℓ hℓ v hv)).2.1⟩)
    (by
      intro u hu ℓ hℓ w hw hdiv
      have := Fact.mk (hP ℓ hℓ).1
      obtain ⟨v, hv, rfl⟩ := mem_image.1 hw
      apply huncov u hu v hv ℓ hℓ
      have h1 := (excludedResidue_spec ℓ q v (hBℓ ℓ hℓ v hv)).2.2
      have e : (v.2 : ℤ) * u - v.1 * q =
          v.2 * (u - excludedResidue ℓ q v) + (v.2 * excludedResidue ℓ q v - v.1 * q) := by ring
      rw [e]; exact dvd_add (dvd_mul_of_dvd_right hdiv _) h1)
  -- ω/(ℓ − ω) ≥ ω/ℓ, so 1 + Σ ω/(ℓ−ω) ≥ |V| Σ 1/ℓ
  have hJ : (V.card : ℝ) * ∑ ℓ ∈ Ps, (1 : ℝ) / ℓ ≤
      1 + ∑ ℓ ∈ Ps, ((Ω ℓ).card : ℝ) / ((ℓ : ℝ) - (Ω ℓ).card) := by
    rw [mul_sum]
    have : ∑ ℓ ∈ Ps, (V.card : ℝ) * (1 / ℓ) ≤
        ∑ ℓ ∈ Ps, ((Ω ℓ).card : ℝ) / ((ℓ : ℝ) - (Ω ℓ).card) := by
      apply sum_le_sum
      intro ℓ hℓ
      rw [hΩcard ℓ hℓ]
      have hlt : (V.card : ℝ) < ℓ := by
        exact_mod_cast lt_of_le_of_lt hVcard (hP ℓ hℓ).2.2.1
      have hV0 : (0 : ℝ) ≤ V.card := by positivity
      rw [mul_one_div]
      apply div_le_div_of_nonneg_left hV0 (by linarith) (by linarith)
    linarith
  calc (U.card : ℝ) * ((V.card : ℝ) * ∑ ℓ ∈ Ps, (1 : ℝ) / ℓ)
      ≤ (U.card : ℝ) * (1 + ∑ ℓ ∈ Ps, ((Ω ℓ).card : ℝ) / ((ℓ : ℝ) - (Ω ℓ).card)) :=
        mul_le_mul_of_nonneg_left hJ (by positivity)
    _ ≤ N₀ + (Q : ℝ) ^ 2 := hLS'

/-- The uncovered set `U_k(q)` of Lemma 9, over the integer interval `(M₀, M₀ + N₀]`. -/
def uncovered (M₀ : ℤ) (N₀ : ℕ) (V : Finset (ℕ × ℕ)) (Ps : Finset ℕ) (q : ℤ) : Finset ℤ :=
  (Ioc M₀ (M₀ + N₀)).filter (fun p => ∀ v ∈ V, ∀ ℓ ∈ Ps, ¬ (ℓ : ℤ) ∣ (v.2 : ℤ) * p - v.1 * q)

/-- **Exceptional set `E_b` (from the large sieve).** Summing Lemma 9 over a finite set `Qs` of
primes `q` (each with its own band `Ps q`, e.g. the least clean band `k(q)`):
`(Σ_{q ∈ Qs} |U(q)|) · J ≤ |Qs| · (N₀ + Q²)`, where `J ≤ |V| Σ_{ℓ ∈ Ps q} 1/ℓ` for all `q`.
The note then uses `J ≥ c₁ H²/n` (PNT mod 4) and `|S_n| ≥ c₀ P²/n²`; those analytic inputs are
NOT formalized. -/
theorem exceptional_Eb_of_largeSieve (hLS : ArithmeticLargeSievePrimeTerms) (M₀ : ℤ) (N₀ Q H : ℕ)
    (V : Finset (ℕ × ℕ)) (Qs : Finset ℤ) (Ps : ℤ → Finset ℕ) (J : ℝ)
    (hV : ∀ v ∈ V, 0 < v.1 ∧ v.1 ≤ H ∧ 0 < v.2 ∧ v.2 ≤ H ∧ Nat.Coprime v.1 v.2)
    (hP : ∀ q ∈ Qs, ∀ ℓ ∈ Ps q, ℓ.Prime ∧ ℓ ≤ Q ∧ H ^ 2 < ℓ ∧ ¬ (ℓ : ℤ) ∣ q)
    (hJ : ∀ q ∈ Qs, J ≤ (V.card : ℝ) * ∑ ℓ ∈ Ps q, (1 : ℝ) / ℓ) :
    (∑ q ∈ Qs, ((uncovered M₀ N₀ V (Ps q) q).card : ℝ)) * J ≤ Qs.card * (N₀ + (Q : ℝ) ^ 2) := by
  rw [sum_mul]
  have : ∀ q ∈ Qs, ((uncovered M₀ N₀ V (Ps q) q).card : ℝ) * J ≤ N₀ + (Q : ℝ) ^ 2 := by
    intro q hq
    have h := lemma9_of_largeSieve hLS M₀ N₀ Q H q V (Ps q) hV (hP q hq)
      (uncovered M₀ N₀ V (Ps q) q)
      (fun u hu => by
        have := (mem_filter.1 hu).1
        rw [mem_Ioc] at this; exact this)
      (fun u hu => (mem_filter.1 hu).2)
    exact (mul_le_mul_of_nonneg_left (hJ q hq) (by positivity)).trans h
  calc ∑ q ∈ Qs, ((uncovered M₀ N₀ V (Ps q) q).card : ℝ) * J
      ≤ ∑ _q ∈ Qs, (N₀ + (Q : ℝ) ^ 2) := sum_le_sum this
    _ = Qs.card * (N₀ + (Q : ℝ) ^ 2) := by rw [sum_const, nsmul_eq_mul]

/-- `E_b` in divided form: if `0 < J`, then `Σ_{q ∈ 𝒬} |U(q)| ≤ |𝒬| (N₀ + Q²) / J`. -/
theorem exceptional_Eb_div_of_largeSieve (hLS : ArithmeticLargeSievePrimeTerms) (M₀ : ℤ)
    (N₀ Q H : ℕ) (V : Finset (ℕ × ℕ)) (Qs : Finset ℤ) (Ps : ℤ → Finset ℕ) (J : ℝ) (hJpos : 0 < J)
    (hV : ∀ v ∈ V, 0 < v.1 ∧ v.1 ≤ H ∧ 0 < v.2 ∧ v.2 ≤ H ∧ Nat.Coprime v.1 v.2)
    (hP : ∀ q ∈ Qs, ∀ ℓ ∈ Ps q, ℓ.Prime ∧ ℓ ≤ Q ∧ H ^ 2 < ℓ ∧ ¬ (ℓ : ℤ) ∣ q)
    (hJ : ∀ q ∈ Qs, J ≤ (V.card : ℝ) * ∑ ℓ ∈ Ps q, (1 : ℝ) / ℓ) :
    ∑ q ∈ Qs, ((uncovered M₀ N₀ V (Ps q) q).card : ℝ) ≤ Qs.card * (N₀ + (Q : ℝ) ^ 2) / J := by
  rw [le_div_iff₀ hJpos]
  exact exceptional_Eb_of_largeSieve hLS M₀ N₀ Q H V Qs Ps J hV hP hJ

/-- Positivity of `J = |V| Σ_{ℓ ∈ 𝒫} 1/ℓ` when `V` and `𝒫` are nonempty and `𝒫` has no `0`. -/
theorem J_pos (V : Finset (ℕ × ℕ)) (Ps : Finset ℕ) (hV : V.Nonempty) (hP : Ps.Nonempty)
    (hP0 : ∀ ℓ ∈ Ps, 0 < ℓ) : 0 < (V.card : ℝ) * ∑ ℓ ∈ Ps, (1 : ℝ) / ℓ := by
  apply mul_pos (by exact_mod_cast hV.card_pos)
  apply Finset.sum_pos (fun ℓ hℓ => by have := hP0 ℓ hℓ; positivity) hP

end FactoringLean
