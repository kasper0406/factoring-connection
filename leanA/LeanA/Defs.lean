import Mathlib

/-!
# Trusted definitions (part 1)

Every definition in this file is copied verbatim into `Challenge.lean`. Comparator checks that the
two copies agree, so a definition here can only be changed together with the challenge.

Conventions.
* `prob s p` / `expect s f` are probability and expectation under the uniform law on a finite set.
* Exponents are rationals (`ℚ`) where the paper says "rational"; asymptotic statements are about `ℝ`.
* Source: `paper2/main.tex` (cited as [pp]) and the project notes cited in each docstring.
-/

noncomputable section
open Classical

namespace FactoringA

open Real Finset

/-! ## §0. Uniform probability on a finite set -/

/-- `prob s p` = Pr[p(x)] for `x` uniform on the finite set `s` (it is `0` when `s = ∅`). -/
def prob {α : Type*} (s : Finset α) (p : α → Prop) : ℝ :=
  ((s.filter p).card : ℝ) / s.card

/-- `expect s f` = E[f(x)] for `x` uniform on the finite set `s` (it is `0` when `s = ∅`). -/
def expect {α : Type*} (s : Finset α) (f : α → ℝ) : ℝ :=
  (∑ x ∈ s, f x) / s.card

/-! ## §1. L-notation -/

/-- `Lam α x = (log x)^α · (max 1 (log log x))^(1-α)`. For `α = 1/2` and `x ≥ e^e` this is the
usual `√(log x · log log x)`; the `max 1` is the convention of [pp] Lemma F, so that the
quantity is defined and `≥ 0` for all `x ≥ 1`. -/
def Lam (α x : ℝ) : ℝ := (Real.log x) ^ α * (max 1 (Real.log (Real.log x))) ^ (1 - α)

/-- `L_x[α, c] = exp(c · Lam α x)`. -/
def LN (α c x : ℝ) : ℝ := Real.exp (c * Lam α x)

/-- `f(N) ≤ L_N[α, c + o(1)]` as `N → ∞`: for every `ε > 0`, eventually `f N ≤ L_N[α, c + ε]`. -/
def LBound (α c : ℝ) (f : ℕ → ℝ) : Prop :=
  ∀ ε > 0, ∀ᶠ N : ℕ in Filter.atTop, f N ≤ LN α (c + ε) N

/-! ## §2. The fixed parameters of Theorem A ([pp] §1.2, §sec:params) -/

/-- The fixed parameters `(d, k, μ, η₀)`: odd degree `d`, number of large primes `k`,
base exponent `μ`, slack `η₀`. -/
structure Params where
  d : ℕ
  k : ℕ
  μ : ℚ
  η₀ : ℚ

namespace Params

/-- `ν = 1/(4k)`. -/
def ν (P : Params) : ℚ := 1 / (4 * (P.k : ℚ))

/-- `h_*(μ) = max(2ν, μ − (d−2)ν, 1 − dμ + ν)` (kernel, level-½ and window floors). -/
def hStar (P : Params) : ℚ :=
  max (2 * P.ν) (max (P.μ - ((P.d : ℚ) - 2) * P.ν) (1 - (P.d : ℚ) * P.μ + P.ν))

/-- `h = h_*(μ) + η₀`. -/
def h (P : Params) : ℚ := P.hStar + P.η₀

/-- `θ = 2μ + h + (d−1)ν` ([pp] eq:theta). -/
def θ (P : Params) : ℚ := 2 * P.μ + P.h + ((P.d : ℚ) - 1) * P.ν

/-- Admissibility ([pp] §1.2 and §sec:params): `d ≥ 3` odd, `k > (d+1)/4`, `η₀ > 0`, `ν + η₀ ≤ μ < 1`, `θ < 1`. -/
def Admissible (P : Params) : Prop :=
  Odd P.d ∧ 3 ≤ P.d ∧ ((P.d : ℚ) + 1) / 4 < P.k ∧ 0 < P.η₀ ∧ P.ν + P.η₀ ≤ P.μ ∧ P.μ < 1 ∧
    P.θ < 1

/-- The constant of Theorem A, `(1 + θ)/√2`. -/
def cA (P : Params) : ℝ := (1 + (P.θ : ℝ)) / Real.sqrt 2

/-- The constant of Theorem RA at the optimal smoothness exponent, `√((1+√2)θ)`. -/
def cRAstar (P : Params) : ℝ := Real.sqrt ((1 + Real.sqrt 2) * (P.θ : ℝ))

end Params

/-- `θ_min(d, k)` of [pp] Lemma lem:optmu (optimal μ, zero slack). -/
def thetaMin (d k : ℕ) : ℚ :=
  if 4 * (k : ℚ) < (d : ℚ) ^ 2 + 1 then
    2 / (d : ℚ) + ((d : ℚ) ^ 2 + d - 2) / d * (1 / (4 * (k : ℚ)))
  else
    3 / ((d : ℚ) + 1) + (4 * (d : ℚ) - 2) / ((d : ℚ) + 1) * (1 / (4 * (k : ℚ)))

/-- The optimal μ of [pp] Lemma lem:optmu. -/
def muOpt (d k : ℕ) : ℚ :=
  if 4 * (k : ℚ) < (d : ℚ) ^ 2 + 1 then (1 - 1 / (4 * (k : ℚ))) / d
  else (1 + ((d : ℚ) - 1) * (1 / (4 * (k : ℚ)))) / ((d : ℚ) + 1)

/-- The cost exponent of SPLIT at smoothness exponent `a` with Lenstra–Pomerance shape test
([pp] Prop. prop:cost and the remark after it): `max(a + 1/√2, a + θ/(2a), 2a)`. -/
def costExpA (θ a : ℝ) : ℝ := max (a + 1 / Real.sqrt 2) (max (a + θ / (2 * a)) (2 * a))

/-- The cost exponent of the recursive algorithm FACTOR_rec[a] ([pp] Theorem thm:RA):
`c(a) = max(a + θ/(2a), (2+√2)a)`. -/
def costExpRA (θ a : ℝ) : ℝ := max (a + θ / (2 * a)) ((2 + Real.sqrt 2) * a)

/-- The cost exponent of DLOG at smoothness exponent `a` ([pp] Theorem thm:DL):
`max(a + 1/√2, a + θ/(2a), 3a)`. -/
def costExpDL (θ a : ℝ) : ℝ := max (a + 1 / Real.sqrt 2) (max (a + θ / (2 * a)) (3 * a))

/-- The parameter sequence of [pp] Corollary cor:P (portfolio): `d_j = 2j+7`,
`k_j = (d_j² + 2d_j − 11)/4`, `μ_j` the optimum of Lemma 3.2, `η₀,j = 1/d_j²`. -/
def portfolioParams (j : ℕ) : Params :=
  let d : ℕ := 2 * j + 7
  let D : ℕ := d ^ 2 + 2 * d - 11
  { d := d, k := D / 4, μ := (1 + ((d : ℚ) - 1) / D) / ((d : ℚ) + 1), η₀ := 1 / (d : ℚ) ^ 2 }

/-! ## §3. Units modulo N, square roots of unity, proper gcds -/

/-- `ω(N)`, the number of distinct prime factors of `N`. -/
def omega (N : ℕ) : ℕ := N.primeFactors.card

/-- `μ₂(ℤ/N) = {ε ∈ (ℤ/N)ˣ : ε² = 1}`. -/
def muTwo (N : ℕ) [NeZero N] : Finset (ZMod N)ˣ := univ.filter (fun u => u ^ 2 = 1)

/-- `Q_N = ((ℤ/N)ˣ)²`, the unit squares. -/
def unitSquares (N : ℕ) [NeZero N] : Finset (ZMod N)ˣ := univ.image (fun u => u ^ 2)

/-- `x` gives a proper factor of `N`: `1 < gcd(x, N) < N` (with `x` read as its least
nonnegative residue). -/
def ProperGcd (N : ℕ) (x : ZMod N) : Prop := 1 < Nat.gcd x.val N ∧ Nat.gcd x.val N < N

/-! ## §4. Split counting ([pp] §sec:split) -/

/-- `E_S(R)`: `u ≡ R v (mod N)` for some `u, v ∈ S`. -/
def SplitEvent (N : ℕ) (S : Finset ℕ) (R : ZMod N) : Prop :=
  ∃ u ∈ S, ∃ v ∈ S, (u : ZMod N) = R * (v : ZMod N)

/-- `(x, y) ∈ Λ_R = {(x, y) ∈ ℤ² : x ≡ R y (mod N)}`. -/
def InSplitLattice (N : ℕ) (R : ZMod N) (x y : ℤ) : Prop := (x : ZMod N) = R * (y : ZMod N)

/-- `(x, y)` is a shortest nonzero vector of `Λ_R` (Euclidean norm). -/
def IsShortestVec (N : ℕ) (R : ZMod N) (x y : ℤ) : Prop :=
  InSplitLattice N R x y ∧ (x, y) ≠ (0, 0) ∧
    ∀ x' y' : ℤ, InSplitLattice N R x' y' → (x', y') ≠ (0, 0) → x ^ 2 + y ^ 2 ≤ x' ^ 2 + y' ^ 2

/-- Products of exactly `k` distinct elements of `D`. -/
def kProducts (D : Finset ℕ) (k : ℕ) : Finset ℕ := (D.powersetCard k).image (fun T => ∏ p ∈ T, p)

/-- `m` is a product of at most `k` distinct elements of `D` (the shape test's condition). -/
def IsProdAtMostK (D : Finset ℕ) (k : ℕ) (m : ℕ) : Prop :=
  ∃ T ⊆ D, T.card ≤ k ∧ m = ∏ p ∈ T, p

/-- `c_k = 8k² C(2k, k)` of [pp] Theorem thm:Sprime. -/
def cK (k : ℕ) : ℝ := 8 * (k : ℝ) ^ 2 * (Nat.choose (2 * k) k : ℝ)

/-! ## §5. The kernel map of [pp] Lemma lem:kern and Lemma Kern-U -/

/-- `v(P) = (a^{d−1−i} b^i)_{0 ≤ i < d}`. -/
def vP (d : ℕ) (a b : ℤ) : Fin d → ℤ := fun i => a ^ (d - 1 - (i : ℕ)) * b ^ (i : ℕ)

/-- `⟨c, v⟩`. -/
def dotZ {d : ℕ} (c v : Fin d → ℤ) : ℤ := ∑ i, c i * v i

/-- `Φ(c) = (⟨c, v(P)⟩, ⟨c, v(P′)⟩)`. -/
def PhiK (d : ℕ) (a b a' b' : ℤ) (c : Fin d → ℤ) : ℤ × ℤ :=
  (dotZ c (vP d a b), dotZ c (vP d a' b'))

/-- `C_H = [−H, H]^d ∩ ℤ^d`. -/
def cube (d H : ℕ) : Finset (Fin d → ℤ) := Fintype.piFinset (fun _ => Icc (-(H : ℤ)) H)

/-! ## §6. Lemma GB (generalised base expansion, [pp] Lemma lem:GB) -/

/-- The digit `a_i ∈ (−g₁/2, g₁/2]` with `a_i ≡ T_i g₀^{−i} (mod g₁)`. -/
def gbDigit (g0 g1 : ℤ) (i : ℕ) (T : ℤ) : ℤ :=
  (((T : ZMod g1.natAbs) * ((g0 : ZMod g1.natAbs)⁻¹) ^ i)).valMinAbs

/-- `T_{i−1} = (T_i − a_i g₀^i)/g₁`. -/
def gbNext (g0 g1 : ℤ) (i : ℕ) (T : ℤ) : ℤ := (T - gbDigit g0 g1 i T * g0 ^ i) / g1

/-- `gbT g0 g1 T d n = T_{d−n}`, starting from `T_d = T`. -/
def gbT (g0 g1 T : ℤ) (d : ℕ) : ℕ → ℤ
  | 0 => T
  | n + 1 => gbNext g0 g1 (d - n) (gbT g0 g1 T d n)

/-- The coefficient `a_i` (`0 ≤ i ≤ d`): `a_i` is the digit of `T_i` for `i ≥ 1`, and `a_0 = T_0`. -/
def gbCoeff (g0 g1 T : ℤ) (d i : ℕ) : ℤ :=
  if i = 0 then gbT g0 g1 T d d else gbDigit g0 g1 i (gbT g0 g1 T d (d - i))

/-- `F̂(x, y) = Σ_{i=0}^{d} a_i x^i y^{d−i}`. -/
def gbForm (g0 g1 T : ℤ) (d : ℕ) (x y : ℤ) : ℤ :=
  ∑ i ∈ range (d + 1), gbCoeff g0 g1 T d i * x ^ i * y ^ (d - i)

/-- `G(x, y) = g₁ x − g₀ y`. -/
def Gform (g0 g1 x y : ℤ) : ℤ := g1 * x - g0 * y

/-- `F_A(x, y) = F̂(x, y) + G(x, y)·⟨c, v(x, y)⟩` for the field `A = (g, w, c)`, with `T = N w`. -/
def FA (g0 g1 T : ℤ) (d : ℕ) (c : Fin d → ℤ) (x y : ℤ) : ℤ :=
  gbForm g0 g1 T d x y + Gform g0 g1 x y * dotZ c (vP d x y)

/-! ## §7. Smoothness -/

/-- `n` is `B`-smooth and nonzero: every prime factor of `|n|` is `≤ B` ([pp] §sec:notation, `𝒮_B`). -/
def IsSmooth (B : ℕ) (n : ℤ) : Prop := n ≠ 0 ∧ ∀ p : ℕ, p.Prime → (p : ℤ) ∣ n → p ≤ B

/-! ## §8. Prime counting -/

/-- The primes in `(a, b]`. -/
def primesIoc (a b : ℕ) : Finset ℕ := (Ioc a b).filter Nat.Prime

/-- `π(x; q, a)`: the number of primes `p ≤ x` with `p ≡ a (mod q)`. -/
def primeCountAP (x : ℝ) (q a : ℕ) : ℕ :=
  ((Finset.range (⌊x⌋₊ + 1)).filter (fun p => p.Prime ∧ p ≡ a [MOD q])).card

end FactoringA
