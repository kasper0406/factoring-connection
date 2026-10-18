/-
  Theorem D′, notes/dcltheorem_results.md §3.
  Lemma 2 (window: AM–GM lower end, convexity upper end, window, indexing) and the
  single-target extraction theorem, which chains Lemmas 3, 4, 5 and 7:
  a covering target of a clean band whose true giant is not FULL makes the
  algorithm's gcd equal to p. The wrapper `single_target_from_lock` additionally
  derives its integrality/indexing hypotheses from Lemma 1 and Lemma 2.
-/
import Mathlib
import FactoringLean.Exponents
import FactoringLean.Bands
import FactoringLean.Coverage

namespace FactoringLean

open Finset

/-! ### Lemma 2 (window) -/

/-- **Lemma 2, lower end (AM–GM).** `Lo_v = ⌊2√(A B N)⌋ = ⌊√(4 A B N)⌋ ≤ B p + A q`. -/
theorem lemma2_lower (A B p q : ℕ) : Nat.sqrt (4 * A * B * (p * q)) ≤ B * p + A * q := by
  have : 4 * A * B * (p * q) ≤ (B * p + A * q) ^ 2 := by
    have h : (0 : ℤ) ≤ ((B : ℤ) * p - A * q) ^ 2 := sq_nonneg _
    have : ((4 * A * B * (p * q) : ℕ) : ℤ) ≤ (((B * p + A * q) ^ 2 : ℕ) : ℤ) := by
      push_cast; nlinarith
    exact_mod_cast this
  calc Nat.sqrt (4 * A * B * (p * q)) ≤ Nat.sqrt ((B * p + A * q) ^ 2) := Nat.sqrt_le_sqrt this
    _ = B * p + A * q := Nat.sqrt_eq' _

/-- **Lemma 2, upper end (convexity).** For `0 < a ≤ x ≤ b` and `A, B, N ≥ 0`,
`F(x) = B x + A N / x ≤ max(F(a), F(b))`. -/
theorem lemma2_upper (A B N a b x : ℝ) (hA : 0 ≤ A) (hN : 0 ≤ N) (ha : 0 < a) (hax : a ≤ x)
    (hxb : x ≤ b) :
    B * x + A * N / x ≤ max (B * a + A * N / a) (B * b + A * N / b) := by
  have hx : 0 < x := lt_of_lt_of_le ha hax
  have hb : 0 < b := lt_of_lt_of_le hx hxb
  have hAN : 0 ≤ A * N := mul_nonneg hA hN
  -- F(x) − F(a) = (x − a)(B − AN/(a x)),  F(x) − F(b) = (x − b)(B − AN/(b x))
  have ea : B * x + A * N / x - (B * a + A * N / a) = (x - a) * (B - A * N / (a * x)) := by
    field_simp; ring
  have eb : B * x + A * N / x - (B * b + A * N / b) = (x - b) * (B - A * N / (b * x)) := by
    field_simp; ring
  by_cases h : B ≤ A * N / (a * x)
  · apply le_max_of_le_left
    have : (x - a) * (B - A * N / (a * x)) ≤ 0 :=
      mul_nonpos_of_nonneg_of_nonpos (by linarith) (by linarith)
    linarith
  · apply le_max_of_le_right
    push Not at h
    have hle : A * N / (b * x) ≤ A * N / (a * x) :=
      div_le_div_of_nonneg_left hAN (by positivity) (by nlinarith)
    have : (x - b) * (B - A * N / (b * x)) ≤ 0 :=
      mul_nonpos_of_nonpos_of_nonneg (by linarith) (by linarith)
    linarith

/-- **Lemma 2, window.** With `K = ℓ² > 0`, `g_lo = ⌈(Lo − c)/K⌉ = −⌊(c − Lo)/K⌋` and
`g_hi = ⌊(Hi − c)/K⌋`: if `Lo ≤ K g + c ≤ Hi` (i.e. `Lo ≤ Bp + Aq ≤ Hi`), then
`g_lo ≤ g ≤ g_hi`. -/
theorem lemma2_window (K g c Lo Hi : ℤ) (hK : 0 < K) (hLo : Lo ≤ K * g + c)
    (hHi : K * g + c ≤ Hi) : -((c - Lo) / K) ≤ g ∧ g ≤ (Hi - c) / K := by
  constructor
  · have : -g ≤ (c - Lo) / K := by
      rw [Int.le_ediv_iff_mul_le hK]; linarith
    linarith
  · rw [Int.le_ediv_iff_mul_le hK]; linarith

/-- **Lemma 2, indexing.** With `i_full = g − g_lo`, `j* = ⌊i_full/m⌋`, `i* = i_full − j* m`:
`0 ≤ i* < m` and `g − (g_lo + j* m) = i*`; and if `g_lo ≤ g ≤ g_hi` then
`0 ≤ j* ≤ (g_hi − g_lo)/m`, so the true giant is among those computed. -/
theorem lemma2_index (g glo ghi m : ℤ) (hm : 0 < m) :
    (0 ≤ (g - glo) % m ∧ (g - glo) % m < m ∧ g - (glo + (g - glo) / m * m) = (g - glo) % m) ∧
    (glo ≤ g → g ≤ ghi → 0 ≤ (g - glo) / m ∧ (g - glo) / m ≤ (ghi - glo) / m) := by
  refine ⟨⟨Int.emod_nonneg _ hm.ne', Int.emod_lt_of_pos _ hm, ?_⟩, fun h1 h2 => ⟨?_, ?_⟩⟩
  · have := Int.ediv_mul_add_emod (g - glo) m
    linarith
  · exact Int.ediv_nonneg (by linarith) hm.le
  · exact Int.ediv_le_ediv hm (by linarith)

/-! ### Single-target extraction -/

/-- **Single-target extraction (Lemmas 3, 4, 5, 7 chained).**

Setting: distinct primes `p, q`, `N = p q`; the base reduces to units `xp` mod `p` and `xq`
mod `q` (both are images of `2` in the note); `u = x^M`, `γ = x^(M ℓ²)`, with
`M = L/ℓ²`. The target covers `p`, so (Lemma 1) `ℓ² g = B p + A q − c`; the true giant has
index `j*` and baby index `i*`, `g − g_lo = j* m + i*`, `0 ≤ i* < m`. The integers `G, γ` are
the algorithm's residues (mod `N`) of the true giant and of `γ`.

Hypotheses on `q`: the band is clean (`gcd(M, ord_q x) = 1`, see `lemma7_coprime_zmod`) and
the FULL condition of Lemma 7 fails for every `δ ∈ [−i*, m − 1 − i*]`.

Conclusion: `gcd(f(G), N) = p`, where `f(X) = ∏_{i<m} (X − γ^i)`. -/
theorem single_target_extraction (p q : ℕ) [hp : Fact p.Prime] [hq : Fact q.Prime]
    (hpq : p ≠ q) (xp : (ZMod p)ˣ) (xq : (ZMod q)ˣ)
    (M A B c ℓ g glo jstar istar : ℤ) (m : ℕ)
    (hg : ℓ ^ 2 * g = B * p + A * q - c)
    (hsplit : g - glo = jstar * m + istar) (hi0 : 0 ≤ istar) (him : istar < m)
    (G γ : ℤ)
    (hγp : (γ : ZMod p) = ((xp ^ (M * ℓ ^ 2) : (ZMod p)ˣ) : ZMod p))
    (hγq : (γ : ZMod q) = ((xq ^ (M * ℓ ^ 2) : (ZMod q)ˣ) : ZMod q))
    (hGp : (G : ZMod p) = (((xp ^ M) ^ (B + A * (p * q) - c) *
        (xp ^ (M * ℓ ^ 2)) ^ (-(glo + jstar * m)) : (ZMod p)ˣ) : ZMod p))
    (hGq : (G : ZMod q) = (((xq ^ M) ^ (B + A * (p * q) - c) *
        (xq ^ (M * ℓ ^ 2)) ^ (-(glo + jstar * m)) : (ZMod q)ˣ) : ZMod q))
    (hclean : IsCoprime M (orderOf xq : ℤ))
    (hnotfull : ∀ δ : ℤ, -istar ≤ δ → δ ≤ m - 1 - istar →
        ¬ (orderOf xq : ℤ) ∣ (A - B) * (p - 1) - ℓ ^ 2 * δ) :
    Int.gcd (∏ i ∈ range m, (G - γ ^ i)) ((p * q : ℕ) : ℤ) = p := by
  have hidx : g - glo - jstar * m = istar := by linarith
  -- Lemma 3: G ≡ γ^{i*} (mod p)
  have h3 := lemma3 p xp M A B q c ℓ g glo jstar m hg
  rw [hidx] at h3
  have hcol : (p : ℤ) ∣ G - γ ^ istar.toNat := by
    rw [← ZMod.intCast_zmod_eq_zero_iff_dvd, Int.cast_sub, Int.cast_pow, hGp, h3, hγp,
      ← Units.val_pow_eq_pow_val, ← zpow_natCast, Int.toNat_of_nonneg hi0, sub_self]
  -- Lemma 4 + Lemma 7: G is not FULL mod q
  have h4 := lemma4 q xq M A B p c ℓ g glo jstar m hg
  rw [hidx] at h4
  have hnf : ∀ i < m, ¬ (q : ℤ) ∣ G - γ ^ i := by
    intro i hi hdiv
    rw [← ZMod.intCast_zmod_eq_zero_iff_dvd, Int.cast_sub, Int.cast_pow, sub_eq_zero, hGq, h4,
      hγq, ← Units.val_pow_eq_pow_val, Units.val_inj, ← zpow_natCast] at hdiv
    obtain ⟨δ, hδ0, hδm, hδ⟩ := (lemma7_full_iff xq M ℓ ((A - B) * (p - 1)) istar m hclean).1
      ⟨i, by positivity, by exact_mod_cast hi, hdiv⟩
    exact hnotfull δ hδ0 hδm hδ
  -- Lemma 5
  exact lemma5_one_sided p q hp.out hq.out hpq G γ m istar.toNat (by omega) hcol hnf

/-- **Single-target extraction from the class lock (Lemmas 1, 2, 3, 4, 5, 7 chained).**
Hypotheses of Lemma 1 (the target `(v, ℓ, r)` covers `p`) and of Lemma 2
(`Lo ≤ B p + A q ≤ Hi`) determine `g = (Bp + Aq − c)/ℓ²`, `g_lo = ⌈(Lo − c)/ℓ²⌉`,
`g_hi = ⌊(Hi − c)/ℓ²⌋`, `j* = ⌊(g − g_lo)/m⌋` and `i* = (g − g_lo) mod m`. Then the true giant
is among those computed (`0 ≤ j* ≤ (g_hi − g_lo)/m`), and if the band is clean and the FULL
condition fails, `gcd(f(G), N) = p`. -/
theorem single_target_from_lock (p q : ℕ) [hp : Fact p.Prime] [hq : Fact q.Prime]
    (hpq : p ≠ q) (xp : (ZMod p)ˣ) (xq : (ZMod q)ˣ)
    (M A B r s c ℓ Lo Hi : ℤ) (m : ℕ) (hm : 0 < m)
    (hℓ : Prime ℓ) (hr : ¬ ℓ ∣ r)
    (hroot : ℓ ∣ B * r ^ 2 - A * (p * q)) (hs : ℓ ^ 2 ∣ r * s - p * q)
    (hc : ℓ ^ 2 ∣ c - (B * r + A * s)) (hpr : ℓ ∣ (p : ℤ) - r)
    (hLo : Lo ≤ B * p + A * q) (hHi : B * p + A * q ≤ Hi)
    (g glo ghi jstar istar : ℤ)
    (hgdef : g = (B * p + A * q - c) / ℓ ^ 2) (hglo : glo = -((c - Lo) / ℓ ^ 2))
    (hghi : ghi = (Hi - c) / ℓ ^ 2)
    (hj : jstar = (g - glo) / m) (hi : istar = (g - glo) % m)
    (G γ : ℤ)
    (hγp : (γ : ZMod p) = ((xp ^ (M * ℓ ^ 2) : (ZMod p)ˣ) : ZMod p))
    (hγq : (γ : ZMod q) = ((xq ^ (M * ℓ ^ 2) : (ZMod q)ˣ) : ZMod q))
    (hGp : (G : ZMod p) = (((xp ^ M) ^ (B + A * (p * q) - c) *
        (xp ^ (M * ℓ ^ 2)) ^ (-(glo + jstar * m)) : (ZMod p)ˣ) : ZMod p))
    (hGq : (G : ZMod q) = (((xq ^ M) ^ (B + A * (p * q) - c) *
        (xq ^ (M * ℓ ^ 2)) ^ (-(glo + jstar * m)) : (ZMod q)ˣ) : ZMod q))
    (hclean : IsCoprime M (orderOf xq : ℤ))
    (hnotfull : ∀ δ : ℤ, -istar ≤ δ → δ ≤ m - 1 - istar →
        ¬ (orderOf xq : ℤ) ∣ (A - B) * (p - 1) - ℓ ^ 2 * δ) :
    (0 ≤ jstar ∧ jstar ≤ (ghi - glo) / m) ∧
      Int.gcd (∏ i ∈ range m, (G - γ ^ i)) ((p * q : ℕ) : ℤ) = p := by
  -- Lemma 1
  have hlock := lemma1_class_lock ℓ A B p q r s c hℓ hr hroot hs hc hpr
  have hg : ℓ ^ 2 * g = B * p + A * q - c := by rw [hgdef]; exact Int.mul_ediv_cancel' hlock
  have hK : 0 < ℓ ^ 2 := by
    have : ℓ ≠ 0 := hℓ.ne_zero
    positivity
  -- Lemma 2
  have hw := lemma2_window (ℓ ^ 2) g c Lo Hi hK (by linarith) (by linarith)
  have hmz : (0 : ℤ) < m := by exact_mod_cast hm
  obtain ⟨⟨hi0, him, hidx⟩, hrange⟩ := lemma2_index g glo ghi m hmz
  rw [← hglo, ← hghi] at hw
  refine ⟨by rw [hj]; exact hrange hw.1 hw.2, ?_⟩
  subst hj hi
  exact single_target_extraction p q hpq xp xq M A B c ℓ g glo ((g - glo) / m) ((g - glo) % m) m
    hg (by linarith [Int.ediv_mul_add_emod (g - glo) (m : ℤ)]) hi0 him G γ hγp hγq hGp hGq
    hclean hnotfull

end FactoringLean
