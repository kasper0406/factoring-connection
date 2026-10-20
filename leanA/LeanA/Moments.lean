import LeanA.Defs

/-!
# (d), (g) Reverse Markov, Paley–Zygmund, field selection, residue spreading

The main statements (MAIN) are copied into `Challenge.lean`; do not change them.
-/

noncomputable section
open Classical

namespace FactoringA

open Finset Real

/-- MAIN. [pp] Lemma lem:revmarkov (reverse Markov): if `Z ∈ [0, 1]` and `E Z ≥ β`, then
`Pr[Z ≥ β/2] ≥ β/2`. -/
theorem reverse_markov {α : Type*} (s : Finset α) (Z : α → ℝ) (hZ : ∀ x ∈ s, 0 ≤ Z x ∧ Z x ≤ 1)
    (β : ℝ) (hβ : β ≤ expect s Z) :
    β / 2 ≤ prob s (fun x => β / 2 ≤ Z x) := by
  sorry

/-- MAIN. [hypD] Prop PZ / [pp] Theorem thm:Dk, Paley–Zygmund in the form used there:
for a nonnegative integer count `X` with `E X > 0`,
`Pr[X ≥ 1] ≥ (E X)² / (E[X(X−1)] + E X)`. No lower bound `E X ≥ 1` is needed. -/
theorem paley_zygmund {α : Type*} (s : Finset α) (X : α → ℕ)
    (hpos : 0 < expect s (fun x => (X x : ℝ))) :
    (expect s (fun x => (X x : ℝ))) ^ 2 /
        (expect s (fun x => (X x : ℝ) * ((X x : ℝ) - 1)) + expect s (fun x => (X x : ℝ))) ≤
      prob s (fun x => 1 ≤ X x) := by
  sorry

/-- MAIN. [pp] Theorem thm:Dk, assembly: if `E X ≥ m/P` and `E[X(X−1)] ≤ P m²` with `m, P > 0`,
then `Pr[X ≥ 1] ≥ 1/(P³ + P/m)`. -/
theorem relation_supply_from_moments {α : Type*} (s : Finset α) (X : α → ℕ) (m P : ℝ)
    (hm : 0 < m) (hP : 0 < P)
    (h1 : m / P ≤ expect s (fun x => (X x : ℝ)))
    (h2 : expect s (fun x => (X x : ℝ) * ((X x : ℝ) - 1)) ≤ P * m ^ 2) :
    1 / (P ^ 3 + P / m) ≤ prob s (fun x => 1 ≤ X x) := by
  sorry

/-- MAIN. [pp] Prop. prop:select(ii), the harmonic-mean identity. Let fields `A` be uniform on `s`,
`ς_A ∈ [0, 1]` the per-call row probability, `b = E ς > 0`, and let the kept field have law
`ς_A/(b·|s|)`. Then `E[1/ς_{A*}] = Pr[ς > 0]/b ≤ 1/b` (with `1/0 = 0`). -/
theorem harmonic_mean_identity {α : Type*} (s : Finset α) (ς : α → ℝ)
    (hς : ∀ A ∈ s, 0 ≤ ς A ∧ ς A ≤ 1) (hb : 0 < expect s ς) :
    ∑ A ∈ s, (ς A / (expect s ς * s.card)) * (1 / ς A) = prob s (fun A => 0 < ς A) / expect s ς ∧
      ∑ A ∈ s, (ς A / (expect s ς * s.card)) * (1 / ς A) ≤ 1 / expect s ς := by
  sorry

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
  sorry

/-! ### Published inputs used by Lemma W -/

/-- PUBLISHED INPUT (Brun–Titchmarsh). H. L. Montgomery and R. C. Vaughan, *The large sieve*,
Mathematika 20 (1973), 119–134, Theorem 2: for positive reals `x, y` and coprime positive integers
`q, a` with `y > q`, `π(x+y; q, a) − π(x; q, a) < 2y/(φ(q) log(y/q))`. -/
def BrunTitchmarsh : Prop :=
  ∀ (x y : ℝ) (q a : ℕ), 0 < x → 0 < y → 1 ≤ q → 1 ≤ a → Nat.Coprime a q → (q : ℝ) < y →
    (primeCountAP (x + y) q a : ℝ) - primeCountAP x q a <
      2 * y / ((Nat.totient q : ℝ) * Real.log (y / q))

/-- PUBLISHED INPUT (Rosser–Schoenfeld). J. B. Rosser and L. Schoenfeld, *Approximate formulas
for some functions of prime numbers*, Illinois J. Math. 6 (1962), 64–94: Corollary 1
(`π(z) < 1.25506 z/log z` for `z > 1`) and Corollary 3 (`π(2z) − π(z) > 3z/(5 log z)` for
`z ≥ 20.5`). -/
def RosserSchoenfeld : Prop :=
  (∀ z : ℝ, 1 < z → (Nat.primeCounting ⌊z⌋₊ : ℝ) < 1.25506 * z / Real.log z) ∧
    (∀ z : ℝ, 20.5 ≤ z →
      3 * z / (5 * Real.log z) < (Nat.primeCounting ⌊2 * z⌋₊ : ℝ) - Nat.primeCounting ⌊z⌋₊)

/-- MAIN. [pp] Lemma lem:W (residue spreading), from Brun–Titchmarsh and Rosser–Schoenfeld: for `w`
uniform among the primes of `(W, 2W]`, every modulus `1 ≤ t < W` and residue `r`,
`Pr[w ≡ r (mod t)] ≤ (24/η₀)/φ(t)`, provided `log W ≤ 2 log N` and `log(W/t) ≥ (η₀/3) log N`. -/
theorem residue_spreading (hBT : BrunTitchmarsh) (hRS : RosserSchoenfeld) (N W t r : ℕ)
    (η₀ : ℝ) (hη : 0 < η₀) (hN : 2 ≤ N) (ht : 1 ≤ t) (htW : t < W) (hW : 21 ≤ W)
    (hlogW : Real.log W ≤ 2 * Real.log N)
    (hgap : η₀ / 3 * Real.log N ≤ Real.log ((W : ℝ) / t)) :
    prob (primesIoc W (2 * W)) (fun w => w ≡ r [MOD t]) ≤ (24 / η₀) / (Nat.totient t : ℝ) := by
  sorry

end FactoringA
