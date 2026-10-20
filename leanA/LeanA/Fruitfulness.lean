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
  sorry

/-- MAIN. `|μ₂(ℤ/N)| = 2^{ω(N)}` for odd `N ≥ 1` ([pp] §sec:notation). -/
theorem card_muTwo (N : ℕ) [NeZero N] (hodd : Odd N) :
    (muTwo N).card = 2 ^ omega N := by
  sorry

/-- MAIN. Consequently the fruitfulness probability of Theorem DD is `1 − 2^{1−ω(N)} ≥ 1/2`
when `N` is odd with `ω(N) ≥ 2`. -/
theorem dd_prob_ge_half (N : ℕ) [NeZero N] (hodd : Odd N) (hω : 2 ≤ omega N) :
    1 - 2 / ((muTwo N).card : ℝ) = 1 - 2 ^ (1 - (omega N : ℤ)) ∧
      (1 : ℝ) / 2 ≤ 1 - 2 / ((muTwo N).card : ℝ) := by
  sorry

/-- MAIN. For `ζ ∈ μ₂(ℤ/N)` with `N` odd: `gcd(ζ − 1, N)` is proper iff `ζ ≠ ±1`
(the last step of the proof of Theorem DD). -/
theorem properGcd_iff_ne_pm_one (N : ℕ) [NeZero N] (hodd : Odd N) (h1 : 1 < N)
    (ζ : (ZMod N)ˣ) (hζ : ζ ^ 2 = 1) :
    ProperGcd N ((ζ : ZMod N) - 1) ↔ ζ ≠ 1 ∧ ζ ≠ -1 := by
  sorry

end FactoringA
