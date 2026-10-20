import LeanA.Defs
import LeanA.Hypotheses

/-!
# (d), (g) Reverse Markov, Paley–Zygmund, field selection, residue spreading

The main statements (MAIN) are copied into `Challenge.lean`; do not change them.
-/

noncomputable section
open Classical

namespace FactoringA

open Finset Real

lemma card_pos_of_expect_pos {α : Type*} {s : Finset α} {f : α → ℝ} (h : 0 < expect s f) :
    (0 : ℝ) < s.card := by
  rcases s.eq_empty_or_nonempty with rfl | hne
  · simp [expect] at h
  · exact_mod_cast hne.card_pos

lemma prob_eq_div {α : Type*} (s : Finset α) (p : α → Prop) [DecidablePred p] :
    prob s p = ((s.filter p).card : ℝ) / s.card := by
  unfold prob
  congr!

lemma prob_nonneg' {α : Type*} (s : Finset α) (p : α → Prop) : 0 ≤ prob s p := by
  unfold prob; positivity

lemma prob_le_one' {α : Type*} (s : Finset α) (p : α → Prop) : prob s p ≤ 1 := by
  unfold prob
  apply div_le_one_of_le₀ _ (by positivity)
  exact_mod_cast card_filter_le s p

/-- MAIN. [pp] Lemma lem:revmarkov (reverse Markov): if `Z ∈ [0, 1]` and `E Z ≥ β`, then
`Pr[Z ≥ β/2] ≥ β/2`. -/
theorem reverse_markov {α : Type*} (s : Finset α) (Z : α → ℝ) (hZ : ∀ x ∈ s, 0 ≤ Z x ∧ Z x ≤ 1)
    (β : ℝ) (hβ : β ≤ expect s Z) :
    β / 2 ≤ prob s (fun x => β / 2 ≤ Z x) := by
  unfold prob
  unfold expect at hβ
  rcases le_or_gt β 0 with hβ0 | hβ0
  · have : 0 ≤ ((s.filter (fun x => β / 2 ≤ Z x)).card : ℝ) / s.card := by positivity
    linarith
  rcases s.eq_empty_or_nonempty with rfl | hne
  · simp at hβ; linarith
  have hn : (0:ℝ) < s.card := by exact_mod_cast hne.card_pos
  rw [le_div_iff₀ hn] at hβ ⊢
  have key : ∑ x ∈ s, Z x ≤ ∑ x ∈ s, ((if β / 2 ≤ Z x then (1:ℝ) else 0) + β / 2) := by
    apply sum_le_sum; intro x hx; obtain ⟨h0, h1⟩ := hZ x hx
    split_ifs with h <;> linarith
  rw [sum_add_distrib, sum_boole, sum_const, nsmul_eq_mul] at key
  linarith

/-- MAIN. [hypD] Prop PZ / [pp] Theorem thm:Dk, Paley–Zygmund in the form used there:
for a nonnegative integer count `X` with `E X > 0`,
`Pr[X ≥ 1] ≥ (E X)² / (E[X(X−1)] + E X)`. No lower bound `E X ≥ 1` is needed. -/
theorem paley_zygmund {α : Type*} (s : Finset α) (X : α → ℕ)
    (hpos : 0 < expect s (fun x => (X x : ℝ))) :
    (expect s (fun x => (X x : ℝ))) ^ 2 /
        (expect s (fun x => (X x : ℝ) * ((X x : ℝ) - 1)) + expect s (fun x => (X x : ℝ))) ≤
      prob s (fun x => 1 ≤ X x) := by
  have hn := card_pos_of_expect_pos hpos
  unfold expect at *
  rw [prob_eq_div]
  set n : ℝ := (s.card : ℝ) with hn_def
  set S1 := ∑ x ∈ s, (X x : ℝ)
  set S2 := ∑ x ∈ s, (X x : ℝ) * ((X x : ℝ) - 1)
  set K : ℝ := ((s.filter (fun x => 1 ≤ X x)).card : ℝ) with hK_def
  have hS1 : 0 < S1 := by
    by_contra h; push Not at h
    have : S1 / n ≤ 0 := div_nonpos_of_nonpos_of_nonneg h hn.le
    linarith
  have hS2 : 0 ≤ S2 := by
    apply sum_nonneg; intro x _
    rcases Nat.eq_zero_or_pos (X x) with h | h
    · simp [h]
    · have : (1:ℝ) ≤ X x := by exact_mod_cast h
      nlinarith
  have hcs := Finset.sum_mul_sq_le_sq_mul_sq s (fun x => (X x : ℝ))
    (fun x => if 1 ≤ X x then (1:ℝ) else 0)
  have e1 : ∑ x ∈ s, (X x : ℝ) * (if 1 ≤ X x then (1:ℝ) else 0) = S1 := by
    apply sum_congr rfl; intro x _
    split_ifs with h
    · ring
    · have : X x = 0 := by omega
      simp [this]
  have e2 : ∑ x ∈ s, (if 1 ≤ X x then (1:ℝ) else 0) ^ 2 = K := by
    rw [hK_def, ← sum_boole]; apply sum_congr rfl; intro x _; split_ifs <;> simp
  have e3 : ∑ x ∈ s, (X x : ℝ) ^ 2 = S2 + S1 := by
    rw [← sum_add_distrib]; apply sum_congr rfl; intro x _; ring
  rw [e1, e2, e3] at hcs
  have hD : 0 < S2 + S1 := by linarith
  have eq : (S1 / n) ^ 2 / (S2 / n + S1 / n) = S1 ^ 2 / (S2 + S1) / n := by
    field_simp
  rw [eq]
  apply div_le_div_of_nonneg_right _ hn.le
  rw [div_le_iff₀ hD]
  linarith

/-- MAIN. [pp] Theorem thm:Dk, assembly: if `E X ≥ m/P` and `E[X(X−1)] ≤ P m²` with `m, P > 0`,
then `Pr[X ≥ 1] ≥ 1/(P³ + P/m)`. -/
theorem relation_supply_from_moments {α : Type*} (s : Finset α) (X : α → ℕ) (m P : ℝ)
    (hm : 0 < m) (hP : 0 < P)
    (h1 : m / P ≤ expect s (fun x => (X x : ℝ)))
    (h2 : expect s (fun x => (X x : ℝ) * ((X x : ℝ) - 1)) ≤ P * m ^ 2) :
    1 / (P ^ 3 + P / m) ≤ prob s (fun x => 1 ≤ X x) := by
  have hE : 0 < expect s (fun x => (X x : ℝ)) := lt_of_lt_of_le (by positivity) h1
  have hpz := paley_zygmund s X hE
  refine le_trans ?_ hpz
  set E := expect s (fun x => (X x : ℝ))
  set A := expect s (fun x => (X x : ℝ) * ((X x : ℝ) - 1))
  have hA : 0 ≤ A := by
    unfold A expect
    apply div_nonneg _ (by positivity)
    apply sum_nonneg; intro x _
    rcases Nat.eq_zero_or_pos (X x) with h | h
    · simp [h]
    · have : (1:ℝ) ≤ X x := by exact_mod_cast h
      nlinarith
  have hEP : m ≤ E * P := by rwa [div_le_iff₀ hP] at h1
  have k1 : P * m ^ 2 ≤ E ^ 2 * P ^ 3 := by
    have : m ^ 2 ≤ (E * P) ^ 2 := pow_le_pow_left₀ hm.le hEP 2
    nlinarith
  have k2 : E ≤ E ^ 2 * (P / m) := by
    have : 1 ≤ E * P / m := by rw [le_div_iff₀ hm]; linarith
    have : E ^ 2 * (P / m) = E * (E * P / m) := by ring
    rw [this]; nlinarith
  rw [div_le_div_iff₀ (by positivity) (by linarith)]
  nlinarith

/-- MAIN. [pp] Prop. prop:select(ii), the harmonic-mean identity. Let fields `A` be uniform on `s`,
`ς_A ∈ [0, 1]` the per-call row probability, `b = E ς > 0`, and let the kept field have law
`ς_A/(b·|s|)`. Then `E[1/ς_{A*}] = Pr[ς > 0]/b ≤ 1/b` (with `1/0 = 0`). -/
theorem harmonic_mean_identity {α : Type*} (s : Finset α) (ς : α → ℝ)
    (hς : ∀ A ∈ s, 0 ≤ ς A ∧ ς A ≤ 1) (hb : 0 < expect s ς) :
    ∑ A ∈ s, (ς A / (expect s ς * s.card)) * (1 / ς A) = prob s (fun A => 0 < ς A) / expect s ς ∧
      ∑ A ∈ s, (ς A / (expect s ς * s.card)) * (1 / ς A) ≤ 1 / expect s ς := by
  have hn := card_pos_of_expect_pos hb
  set b := expect s ς with hb_def
  have hterm : ∀ A ∈ s, ς A / (b * s.card) * (1 / ς A) =
      if 0 < ς A then 1 / (b * s.card) else 0 := by
    intro A hA
    split_ifs with h
    · field_simp
    · have : ς A = 0 := le_antisymm (not_lt.mp h) (hς A hA).1
      simp [this]
  have hsum : ∑ A ∈ s, ς A / (b * s.card) * (1 / ς A) = prob s (fun A => 0 < ς A) / b := by
    rw [sum_congr rfl hterm, ← sum_filter, sum_const, nsmul_eq_mul]
    unfold prob
    field_simp
  refine ⟨hsum, ?_⟩
  rw [hsum]
  exact div_le_div_of_nonneg_right (prob_le_one' _ _) hb.le

/-- MAIN. [pp] Prop. prop:select(ii), cost bookkeeping. With trial success probability `p ≥ b > 0`
(`b = E ς`), selection cost `C/p`, the kept-field branch reached with probability `b/p` and costing
at most `2(M·C/ς_A + T)` given `A* = A` (two expected DD rounds), the total is
`≤ 3(M·C + T)/p`. -/
theorem split_cost_bookkeeping {α : Type*} (s : Finset α) (ς : α → ℝ)
    (hς : ∀ A ∈ s, 0 ≤ ς A ∧ ς A ≤ 1) (hb : 0 < expect s ς) (p M C T : ℝ)
    (hp : expect s ς ≤ p) (hM : 1 ≤ M) (hC : 0 ≤ C) (hT : 0 ≤ T) :
    C / p + (expect s ς / p) *
        ∑ A ∈ s, (ς A / (expect s ς * s.card)) * (2 * (M * C * (1 / ς A) + T)) ≤
      3 * (M * C + T) / p := by
  have hn := card_pos_of_expect_pos hb
  have H := harmonic_mean_identity s ς hς hb
  set b := expect s ς with hb_def
  set H0 := ∑ A ∈ s, (ς A / (b * s.card)) * (1 / ς A)
  have hH0 : 0 ≤ H0 := by
    apply sum_nonneg; intro A hA
    have := (hς A hA).1
    positivity
  have hbH : b * H0 ≤ 1 := by
    have := H.2
    rw [le_div_iff₀ hb] at this; linarith
  have hb1 : b ≤ 1 := by
    rw [hb_def]; unfold expect
    rw [div_le_one hn]
    calc ∑ x ∈ s, ς x ≤ ∑ x ∈ s, (1:ℝ) := sum_le_sum (fun x hx => (hς x hx).2)
      _ = s.card := by simp
  have hsumpos : 0 < ∑ i ∈ s, ς i := by
    have h' : 0 < expect s ς := hb
    unfold expect at h'
    by_contra hc
    exact absurd h' (not_lt.mpr (div_nonpos_of_nonpos_of_nonneg (not_lt.mp hc) hn.le))
  have hone : ∑ A ∈ s, ς A / (b * s.card) = 1 := by
    rw [← sum_div, hb_def]; unfold expect
    rw [div_mul_cancel₀ _ hn.ne']; exact div_self hsumpos.ne'
  have hexp : ∑ A ∈ s, (ς A / (b * s.card)) * (2 * (M * C * (1 / ς A) + T)) =
      2 * (M * C) * H0 + 2 * T * ∑ A ∈ s, ς A / (b * s.card) := by
    rw [mul_sum, mul_sum, ← sum_add_distrib]
    apply sum_congr rfl; intro A _; ring
  rw [hexp, hone]
  have hp0 : 0 < p := lt_of_lt_of_le hb hp
  have eq : C / p + b / p * (2 * (M * C) * H0 + 2 * T * 1) =
      (C + 2 * (M * C) * (b * H0) + 2 * T * b) / p := by
    field_simp; ring
  rw [eq]
  apply div_le_div_of_nonneg_right _ hp0.le
  have hMC : C ≤ M * C := by nlinarith
  have hMC0 : 0 ≤ M * C := by nlinarith
  nlinarith

/-! ### Published inputs used by Lemma W -/

lemma card_filter_Ioc_add_range (P : ℕ → Prop) [DecidablePred P] {a b : ℕ} (hab : a ≤ b) :
    ((Ioc a b).filter P).card + ((range (a + 1)).filter P).card =
      ((range (b + 1)).filter P).card := by
  rw [← card_union_of_disjoint]
  · congr 1
    ext x
    simp only [mem_union, mem_filter, mem_Ioc, mem_range]
    constructor
    · rintro (⟨⟨h1, h2⟩, h3⟩ | ⟨h1, h3⟩) <;> exact ⟨by omega, h3⟩
    · rintro ⟨h1, h3⟩
      by_cases h : a < x
      · exact Or.inl ⟨⟨h, by omega⟩, h3⟩
      · exact Or.inr ⟨by omega, h3⟩
  · rw [disjoint_left]
    intro x h1 h2
    simp only [mem_filter, mem_Ioc, mem_range] at h1 h2
    omega

lemma primeCounting_eq_card (n : ℕ) :
    Nat.primeCounting n = ((range (n + 1)).filter Nat.Prime).card := by
  simp only [Nat.primeCounting, Nat.primeCounting', Nat.count_eq_card_filter_range]

/-- MAIN. [pp] Lemma lem:W (residue spreading), from Brun–Titchmarsh and Rosser–Schoenfeld: for `w`
uniform among the primes of `(W, 2W]`, every modulus `1 ≤ t < W` and residue `r`,
`Pr[w ≡ r (mod t)] ≤ (24/η₀)/φ(t)`, provided `log W ≤ 2 log N` and `log(W/t) ≥ (η₀/3) log N`. -/
theorem residue_spreading (hBT : BrunTitchmarsh) (hRS : RosserSchoenfeld) (N W t r : ℕ)
    (η₀ : ℝ) (hη : 0 < η₀) (hN : 2 ≤ N) (ht : 1 ≤ t) (htW : t < W) (hW : 21 ≤ W)
    (hlogW : Real.log W ≤ 2 * Real.log N)
    (hgap : η₀ / 3 * Real.log N ≤ Real.log ((W : ℝ) / t)) :
    prob (primesIoc W (2 * W)) (fun w => w ≡ r [MOD t]) ≤ (24 / η₀) / (Nat.totient t : ℝ) := by
  have hW0 : (0:ℝ) < W := by exact_mod_cast (by omega : 0 < W)
  have hlogN : 0 < Real.log N := Real.log_pos (by exact_mod_cast (by omega : 1 < N))
  have hlogW0 : 0 < Real.log W := Real.log_pos (by exact_mod_cast (by omega : 1 < W))
  have ht0 : (0:ℝ) < t := by exact_mod_cast (by omega : 0 < t)
  have hL1 : 0 < Real.log ((W:ℝ) / t) := lt_of_lt_of_le (by positivity) hgap
  have hφ : (0:ℝ) < Nat.totient t := by exact_mod_cast Nat.totient_pos.mpr (by omega)
  have hRHS : 0 ≤ (24 / η₀) / (Nat.totient t : ℝ) := by positivity
  by_cases hcop : Nat.Coprime r t
  swap
  · have hempty : (primesIoc W (2 * W)).filter (fun w => w ≡ r [MOD t]) = ∅ := by
      rw [Finset.filter_eq_empty_iff]
      intro w hw hmod
      simp only [primesIoc, mem_filter, mem_Ioc] at hw
      obtain ⟨⟨hw1, _⟩, hwp⟩ := hw
      have hg : Nat.gcd w t = Nat.gcd r t := hmod.gcd_eq
      rcases Nat.coprime_or_dvd_of_prime hwp t with h | h
      · exact hcop (by rw [Nat.Coprime, ← hg]; exact h)
      · have := Nat.le_of_dvd (by omega) h; omega
    rw [prob_eq_div, hempty, card_empty, Nat.cast_zero, zero_div]
    exact hRHS
  obtain ⟨a, ha⟩ : ∃ a, a = r % t + t := ⟨_, rfl⟩
  have ha1 : 1 ≤ a := by omega
  have har : a ≡ r [MOD t] := by
    unfold Nat.ModEq
    rw [ha, Nat.add_mod_right, Nat.mod_mod]
  have hacop : Nat.Coprime a t := by
    rw [Nat.Coprime, har.gcd_eq]; exact hcop
  have hF : (primesIoc W (2 * W)).filter (fun w => w ≡ r [MOD t]) =
      (Ioc W (2 * W)).filter (fun p => p.Prime ∧ p ≡ a [MOD t]) := by
    unfold primesIoc
    rw [filter_filter]
    apply filter_congr
    intro x _
    constructor
    · rintro ⟨h1, h2⟩; exact ⟨h1, h2.trans har.symm⟩
    · rintro ⟨h1, h2⟩; exact ⟨h1, h2.trans har⟩
  have e1 : ⌊(W:ℝ)⌋₊ = W := Nat.floor_natCast W
  have e2 : ⌊(W:ℝ) + W⌋₊ = 2 * W := by
    rw [show (W:ℝ) + W = ((2 * W : ℕ) : ℝ) by push_cast; ring]; exact Nat.floor_natCast _
  have e3 : ⌊2 * (W:ℝ)⌋₊ = 2 * W := by
    rw [show 2 * (W:ℝ) = ((2 * W : ℕ) : ℝ) by push_cast; ring]; exact Nat.floor_natCast _
  have hcountF : ((Ioc W (2 * W)).filter (fun p => p.Prime ∧ p ≡ a [MOD t])).card +
      primeCountAP W t a = primeCountAP ((W:ℝ) + W) t a := by
    unfold primeCountAP; rw [e1, e2]
    exact card_filter_Ioc_add_range _ (by omega)
  have hcountS : (primesIoc W (2 * W)).card + Nat.primeCounting W =
      Nat.primeCounting (2 * W) := by
    rw [primeCounting_eq_card, primeCounting_eq_card]
    exact card_filter_Ioc_add_range _ (by omega)
  have hBT' := hBT W W t a hW0 hW0 ht ha1 hacop (by exact_mod_cast htW)
  have h205 : (20.5:ℝ) ≤ W := by
    have : (21:ℝ) ≤ W := by exact_mod_cast hW
    have h2 : (20.5:ℝ) = 41 / 2 := by norm_num
    rw [h2]; linarith
  have hRS' := hRS.2 W h205
  rw [e1, e3] at hRS'
  set F := ((Ioc W (2 * W)).filter (fun p => p.Prime ∧ p ≡ a [MOD t])).card with hF_def
  set L1 := Real.log ((W:ℝ) / t)
  set L := Real.log (W:ℝ)
  set φ : ℝ := (Nat.totient t : ℝ)
  have hFreal : (F:ℝ) < 2 * W / (φ * L1) := by
    have : (F:ℝ) = primeCountAP ((W:ℝ) + W) t a - primeCountAP W t a := by
      rw [← hcountF]; push_cast; ring
    linarith
  have hSreal : 3 * W / (5 * L) < ((primesIoc W (2 * W)).card : ℝ) := by
    have : ((primesIoc W (2 * W)).card : ℝ) = Nat.primeCounting (2 * W) - Nat.primeCounting W := by
      rw [← hcountS]; push_cast; ring
    linarith
  have hSpos : 0 < 3 * (W:ℝ) / (5 * L) := by positivity
  rw [prob_eq_div, hF]
  have hWne : (W:ℝ) ≠ 0 := hW0.ne'
  have hφne : φ ≠ 0 := hφ.ne'
  have hL1ne : L1 ≠ 0 := hL1.ne'
  have hLne : L ≠ 0 := hlogW0.ne'
  calc (F:ℝ) / ((primesIoc W (2 * W)).card : ℝ)
      ≤ (2 * W / (φ * L1)) / (3 * W / (5 * L)) :=
        div_le_div₀ (by positivity) hFreal.le hSpos hSreal.le
    _ = 10 * L / (3 * φ * L1) := by field_simp; try ring
    _ ≤ 24 / η₀ / φ := by
        rw [div_div, div_le_div_iff₀ (by positivity) (by positivity)]
        have k1 : L * (η₀ * φ) ≤ 2 * Real.log N * (η₀ * φ) :=
          mul_le_mul_of_nonneg_right hlogW (by positivity)
        have k2 : φ * (η₀ / 3 * Real.log N) ≤ φ * L1 := mul_le_mul_of_nonneg_left hgap hφ.le
        have k3 : 0 < Real.log N * (η₀ * φ) := by positivity
        nlinarith

end FactoringA
