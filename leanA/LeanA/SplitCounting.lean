import LeanA.Defs

/-!
# (c) Split counting: [pp] Theorems thm:S, thm:Sprime and the combinatorial core of thm:K

`R` is uniform on the unit squares `Q_N` (`unitSquares N`). `SplitEvent N S R` is the event
`u ≡ R v (mod N)` for some `u, v ∈ S`.

The main statements (MAIN) are copied into `Challenge.lean`; do not change them.
-/

noncomputable section
open Classical

namespace FactoringA

open Finset Real

/-- MAIN. `[(ℤ/N)ˣ : Q_N] = |μ₂(ℤ/N)|` (squaring is a homomorphism with kernel `μ₂`), so for odd `N`
the index is `2^{ω(N)}`. -/
theorem card_unitSquares_mul (N : ℕ) [NeZero N] :
    (unitSquares N).card * (muTwo N).card = Fintype.card (ZMod N)ˣ := by
  sorry

/-- MAIN. [pp] Theorem thm:S (split counting; ratio form of the Dixon–Pomerance lemma), probability
bound. -/
theorem split_counting (N : ℕ) [NeZero N] (hodd : Odd N) (Y' : ℝ) (hY : 2 ≤ Y')
    (hYN : 2 * Y' ^ 2 < N) (S : Finset ℕ)
    (hS : ∀ s ∈ S, Y' / 2 ≤ (s : ℝ) ∧ (s : ℝ) ≤ Y' ∧ Nat.Coprime s N) :
    (S.card : ℝ) ^ 4 / (2 * 2 ^ omega N * (Nat.totient N : ℝ) * Y' ^ 2 * (1 + Real.log Y')) ≤
      prob (unitSquares N) (fun R => SplitEvent N S (R : ZMod N)) := by
  sorry

/-- MAIN. [pp] Theorem thm:S, second part: if `2Y'² < N`, `u, v ∈ [1, Y']` and `u ≡ R v`, then
`(u, v)` is an integer multiple of every shortest nonzero vector of `Λ_R`. -/
theorem split_multiple_of_shortest (N : ℕ) (Y' : ℝ) (hYN : 2 * Y' ^ 2 < N) (R : ZMod N)
    (u v : ℕ) (hu : 1 ≤ u) (hv : 1 ≤ v) (huY : (u : ℝ) ≤ Y') (hvY : (v : ℝ) ≤ Y')
    (huv : (u : ZMod N) = R * (v : ZMod N)) (x y : ℤ) (hxy : IsShortestVec N R x y) :
    ∃ j : ℤ, (u : ℤ) = j * x ∧ (v : ℤ) = j * y := by
  sorry

/-- MAIN. [pp] Theorem thm:Sprime (k-products). `D` is a set of primes in a dyadic interval
`(Z, 2Z]`, all coprime to `N`, with `|D| ≥ 4k²`; `S` is a set of products of exactly `k` distinct
primes of `D` in the window `[Y'/2, Y']`, `2Y'² < N`, and `|S| ≥ C(|D|, k)/(2k)`. -/
theorem split_counting_kprod (N : ℕ) [NeZero N] (hodd : Odd N) (k : ℕ) (hk : 1 ≤ k) (Z : ℕ)
    (D : Finset ℕ) (hD : ∀ p ∈ D, p.Prime ∧ Z < p ∧ p ≤ 2 * Z ∧ Nat.Coprime p N)
    (hDcard : 4 * k ^ 2 ≤ D.card) (Y' : ℝ) (hYN : 2 * Y' ^ 2 < N) (S : Finset ℕ)
    (hSsub : S ⊆ kProducts D k) (hSwin : ∀ s ∈ S, Y' / 2 ≤ (s : ℝ) ∧ (s : ℝ) ≤ Y')
    (hScard : (Nat.choose D.card k : ℝ) / (2 * k) ≤ S.card) :
    (S.card : ℝ) ^ 2 / (cK k * 2 ^ omega N * (Nat.totient N : ℝ)) ≤
      prob (unitSquares N) (fun R => SplitEvent N S (R : ZMod N)) := by
  sorry

/-- MAIN. [pp] Theorem thm:K, combinatorial core. Let `N` be odd, `X` with `2X² < N`, and `Y` with
`(2Y)^k ≤ X < (4Y)^k` (as for `J = (Y, 2Y]` in [pp] §sec:params). Let `D ⊆ J` be primes coprime to
`N` with `|D| ≥ 4k²`. For `r` uniform on `[1, N−1]`, with probability at least
`(|D|/k)^{2k} / (4k² c_k 2^{ω(N)} N)` the root `r` is a unit and every shortest nonzero vector
`(u₀, v₀)` of `Λ_{r²}` has `u₀v₀ ≠ 0` and both `|u₀|`, `|v₀|` products of at most `k` distinct
primes of `D` (so the shape test of DESCEND_ONCE passes with targets in `D`). -/
theorem complete_descent_prob (N k X Y : ℕ) [NeZero N] (hodd : Odd N) (hk : 1 ≤ k)
    (hXN : 2 * X ^ 2 < N) (hYX : (2 * Y) ^ k ≤ X) (hXY : X < (4 * Y) ^ k)
    (D : Finset ℕ) (hD : ∀ p ∈ D, p.Prime ∧ Y < p ∧ p ≤ 2 * Y ∧ Nat.Coprime p N)
    (hDcard : 4 * k ^ 2 ≤ D.card) :
    ((D.card : ℝ) / k) ^ (2 * k) / (4 * (k : ℝ) ^ 2 * cK k * 2 ^ omega N * N) ≤
      prob (Icc 1 (N - 1)) (fun r : ℕ => IsUnit (r : ZMod N) ∧
        ∀ x y : ℤ, IsShortestVec N ((r : ZMod N) ^ 2) x y →
          x * y ≠ 0 ∧ IsProdAtMostK D k x.natAbs ∧ IsProdAtMostK D k y.natAbs) := by
  sorry

end FactoringA
