import LeanA.Defs

/-!
# (e) Lemma Kern-U and Lemma Kern; field-law algebra (Lemmas GB, homogeneity, forced residue,
pairs)

The main statements (MAIN) are copied into `Challenge.lean`; do not change them.
-/

noncomputable section
open Classical

namespace FactoringA

open Finset Real

/-! ### Lemma GB ([pp] Lemma lem:GB) -/

/-- MAIN. [pp] Lemma lem:GB: `F̂(g₀, g₁) = T`. -/
theorem gb_identity (g0 g1 T : ℤ) (d : ℕ) (hg1 : 2 ≤ g1) (hcop : IsCoprime g0 g1) :
    gbForm g0 g1 T d g0 g1 = T := by
  sorry

/-- MAIN. [pp] Lemma lem:GB: `|a_i| ≤ g₁/2` for `1 ≤ i ≤ d`. -/
theorem gb_digit_bound (g0 g1 T : ℤ) (d i : ℕ) (hg1 : 2 ≤ g1) (hi1 : 1 ≤ i) (hid : i ≤ d) :
    2 * |gbCoeff g0 g1 T d i| ≤ g1 := by
  sorry

/-- MAIN. [pp] Lemma lem:GB: if `g₀, g₁ ∈ [M/2, M]` then `|a₀| ≤ 2^d T/M^d + d 2^{d−1} M`. -/
theorem gb_a0_bound (g0 g1 T : ℤ) (d : ℕ) (M : ℝ) (hg1 : 2 ≤ g1) (hcop : IsCoprime g0 g1)
    (hT : 0 ≤ T) (h0 : M / 2 ≤ g0) (h0' : (g0 : ℝ) ≤ M) (h1 : M / 2 ≤ g1) (h1' : (g1 : ℝ) ≤ M) :
    |(gbCoeff g0 g1 T d 0 : ℝ)| ≤ 2 ^ d * T / M ^ d + d * 2 ^ (d - 1) * M := by
  sorry

/-- MAIN. [pp] Lemma lem:GB: `a_d ≢ 0 (mod g₁)` when `T` is a unit modulo `g₁`. -/
theorem gb_ad_not_dvd (g0 g1 T : ℤ) (d : ℕ) (hd : 1 ≤ d) (hg1 : 2 ≤ g1)
    (hcop : IsCoprime g0 g1) (hT : IsCoprime T g1) :
    ¬ g1 ∣ gbCoeff g0 g1 T d d := by
  sorry

/-- MAIN. [pp] Lemma lem:homog: `g₁^d F_A(P) ≡ T b^d (mod G(P))`. -/
theorem homogeneity (g0 g1 T : ℤ) (d : ℕ) (c : Fin d → ℤ) (a b : ℤ) (hg1 : 2 ≤ g1)
    (hcop : IsCoprime g0 g1) :
    g1 ^ d * FA g0 g1 T d c a b ≡ T * b ^ d [ZMOD Gform g0 g1 a b] := by
  sorry

/-- MAIN. [pp] Lemma lem:coprime (coprimality and the forced residue). Here `T = N w`,
`t = |G(P)| = ℓ|n|` with `n = G(P)/ℓ` `B`-smooth, `b` a prime with `B < b < ℓ`, `ℓ ∤ N` prime,
every prime factor of `N` exceeds `B`, and `w > t` is prime. Then `gcd(t, b g₁ N w) = 1` and
`F̂(P) ≡ N b^d g₁^{−d} w (mod t)`. -/
theorem forced_residue (N w B ℓ : ℕ) (g0 g1 a b : ℤ) (d : ℕ) (hg1 : 2 ≤ g1)
    (hcop : IsCoprime g0 g1) (hb : (b.natAbs).Prime) (hb0 : 0 < b) (hBb : (B : ℤ) < b)
    (hbℓ : b < ℓ) (hℓ : ℓ.Prime) (hℓN : ¬ ℓ ∣ N) (hdiv : (ℓ : ℤ) ∣ Gform g0 g1 a b)
    (hsm : IsSmooth B (Gform g0 g1 a b / ℓ)) (hNB : ∀ p : ℕ, p.Prime → p ∣ N → B < p)
    (hw : w.Prime) (hwt : |Gform g0 g1 a b| < w) :
    Nat.Coprime (Gform g0 g1 a b).natAbs (b.natAbs * g1.natAbs * N * w) ∧
      ((gbForm g0 g1 ((N : ℤ) * w) d a b : ℤ) : ZMod (Gform g0 g1 a b).natAbs) =
        (N : ZMod (Gform g0 g1 a b).natAbs) * (b : ZMod (Gform g0 g1 a b).natAbs) ^ d *
          ((g1 : ZMod (Gform g0 g1 a b).natAbs)⁻¹) ^ d * (w : ZMod (Gform g0 g1 a b).natAbs) := by
  sorry

/-! ### Lemma Kern ([pp] Lemma lem:kern) and Lemma Kern-U (lv5unif §3) -/

/-- MAIN. [pp] Lemma lem:kern(ii): for primitive `P = (a, b)`, `P′ = (a′, b′)` with `b, b′ ≥ 1`
and `Δ ≠ 0`, `ker Φ` consists exactly of the coefficient vectors of `Q·E` with
`Q = (bx − ay)(b′x − a′y)` and `E` a form of degree `d − 3`; in coordinates
`c_i = bb′ e_i − (ab′ + a′b) e_{i−1} + aa′ e_{i−2}` (with `e_j = 0` outside `[0, d−3]`). -/
theorem kern_kernel (d : ℕ) (hd : 3 ≤ d) (a b a' b' : ℤ) (hb : 1 ≤ b) (hb' : 1 ≤ b')
    (hP : IsCoprime a b) (hP' : IsCoprime a' b') (hΔ : a * b' - a' * b ≠ 0) (c : Fin d → ℤ) :
    PhiK d a b a' b' c = (0, 0) ↔
      ∃ e : ℕ → ℤ, (∀ j, d - 2 ≤ j → e j = 0) ∧
        ∀ i : Fin d, c i = b * b' * e i - (a * b' + a' * b) * (if 1 ≤ (i : ℕ) then e (i - 1) else 0) +
          a * a' * (if 2 ≤ (i : ℕ) then e (i - 2) else 0) := by
  sorry

/-- MAIN. Lemma Kern-U (lv5unif §3), count: every fibre of `Φ` on `C_H` has at most
`(⌊2H/(bb′)⌋ + 1)^{d−2}` points. -/
theorem kernU_count (d : ℕ) (hd : 3 ≤ d) (a b a' b' : ℤ) (hb : 1 ≤ b) (hb' : 1 ≤ b')
    (hP : IsCoprime a b) (hP' : IsCoprime a' b') (hΔ : a * b' - a' * b ≠ 0) (H : ℕ)
    (z : ℤ × ℤ) :
    (((cube d H).filter (fun c => PhiK d a b a' b' c = z)).card : ℤ) ≤
      ((2 * (H : ℤ)) / (b * b') + 1) ^ (d - 2) := by
  sorry

/-- MAIN. Lemma Kern-U (lv5unif §3), probability: if `bb′ d ≤ 2H + 1` then
`Pr_{c ∈ C_H}[Φ(c) = z] ≤ e/((2H+1)² (bb′)^{d−2})`. -/
theorem kernU_prob (d : ℕ) (hd : 3 ≤ d) (a b a' b' : ℤ) (hb : 1 ≤ b) (hb' : 1 ≤ b')
    (hP : IsCoprime a b) (hP' : IsCoprime a' b') (hΔ : a * b' - a' * b ≠ 0) (H : ℕ)
    (hH : b * b' * d ≤ 2 * (H : ℤ) + 1) (z : ℤ × ℤ) :
    prob (cube d H) (fun c => PhiK d a b a' b' c = z) ≤
      Real.exp 1 / ((2 * (H : ℝ) + 1) ^ 2 * ((b * b' : ℤ) : ℝ) ^ (d - 2)) := by
  sorry

/-- MAIN. [pp] Lemma lem:kern(iii) in the form used for Theorem A (explicit `C_d = 3^{d−2}/2^d`):
if `|a|, |b|, |a′|, |b′| ≤ V` and `H ≥ 4V²` then
`Pr_{c ∈ C_H}[Φ(c) = z] ≤ 3^{d−2}/(2^d H² (bb′)^{d−2})`. -/
theorem kern_prob (d : ℕ) (hd : 3 ≤ d) (a b a' b' : ℤ) (hb : 1 ≤ b) (hb' : 1 ≤ b')
    (hP : IsCoprime a b) (hP' : IsCoprime a' b') (hΔ : a * b' - a' * b ≠ 0) (H V : ℕ)
    (ha : |a| ≤ V) (hbV : b ≤ V) (ha' : |a'| ≤ V) (hbV' : b' ≤ V) (hH : 4 * V ^ 2 ≤ H)
    (z : ℤ × ℤ) :
    prob (cube d H) (fun c => PhiK d a b a' b' c = z) ≤
      3 ^ (d - 2) / (2 ^ d * (H : ℝ) ^ 2 * ((b * b' : ℤ) : ℝ) ^ (d - 2)) := by
  sorry

/-- MAIN. [pp] Lemma lem:pairs. Let `ℓ` be prime, `g = (g₀, g₁)` coprime, and `P ≠ P′` points
`(a, b)` with `b > 0`, `gcd(a, b) = 1`, `|a|, |b| ≤ V`, both in `Λ_ℓ(g) = {ℓ ∣ g₁x − g₀y}`.
Then `Δ = det(P, P′) ≠ 0`, `ℓ ∣ Δ` and `|Δ| ≤ 2V²` (so `q = |Δ|/ℓ ≤ 2S²` when `V² ≤ ℓS²`). -/
theorem pairs (ℓ : ℕ) (hℓ : ℓ.Prime) (g0 g1 : ℤ) (hg : IsCoprime g0 g1) (V : ℕ)
    (a b a' b' : ℤ) (hb : 0 < b) (hb' : 0 < b') (hP : IsCoprime a b) (hP' : IsCoprime a' b')
    (ha : |a| ≤ V) (hbV : b ≤ V) (ha' : |a'| ≤ V) (hbV' : b' ≤ V)
    (hL : (ℓ : ℤ) ∣ Gform g0 g1 a b) (hL' : (ℓ : ℤ) ∣ Gform g0 g1 a' b') (hne : (a, b) ≠ (a', b')) :
    a * b' - a' * b ≠ 0 ∧ (ℓ : ℤ) ∣ a * b' - a' * b ∧ |a * b' - a' * b| ≤ 2 * (V : ℤ) ^ 2 := by
  sorry

/-- MAIN. [pp] Lemma lem:pairs, multiplicity: for `P = (a, b)` with `2b > V` and `D ∈ ℤ`, at most
`6` points `P′ ∈ [−V, V]²` have `det(P, P′) = D`. -/
theorem pairs_count (V : ℕ) (a b : ℤ) (hb : (V : ℤ) < 2 * b) (hP : IsCoprime a b) (D : ℤ) :
    ((((Icc (-(V : ℤ)) V) ×ˢ (Icc (-(V : ℤ)) V)).filter
      (fun Q : ℤ × ℤ => a * Q.2 - Q.1 * b = D)).card) ≤ 6 := by
  sorry

end FactoringA
