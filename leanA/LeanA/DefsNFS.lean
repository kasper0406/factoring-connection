import LeanA.Defs

/-!
# Trusted definitions (part 2): the concrete design of [pp] §sec:params–§sec:candidates

Every definition here is copied verbatim into `Challenge.lean`.

All quantities depend on the fixed parameters `P : Params`, on the smoothness exponent `a > 0`
(`a = θ/√2` for Theorem A; free for Theorem RA, [pp] Lemma lem:generala) and on the integer `N`
handed to SPLIT. Following [pp] §sec:params and Lemma lem:generala:
`B = ⌊L^a⌋`, `S = ⌈L^{θ/(4a)}⌉` (so `S = ⌈L^{1/(2√2)}⌉` at `a = θ/√2`), `X = ⌊√(N/2)⌋`,
`j` maximal with `2^{kj} ≤ X`, `Y = 2^{j−1}`, `J = primes in (Y, 2Y]`, `V_ℓ = ⌊√ℓ S⌋`,
`V_max = ⌊√(2Y) S⌋`, `M = ⌊N^μ⌋`, `H = ⌊N^h⌋`, `W = ⌈N^{μ+ν+η₀/2}⌉`,
`X_F = (2^{2d+5} + 2d) M H V_max^d`, `X_r = 2 M V_max / Y`, `ρ_r = ϱ(X_r)`, `ρ_F = ϱ(X_F)`.
-/

noncomputable section
open Classical

namespace FactoringA

open Real Finset

/-! ## §9. Smooth-number counting -/

/-- `Ψ(x, y) = #{1 ≤ n ≤ x : every prime factor of n is ≤ y}`. -/
def Psi (x y : ℝ) : ℕ :=
  ((Icc 1 ⌊x⌋₊).filter (fun n => ∀ p : ℕ, p.Prime → p ∣ n → (p : ℝ) ≤ y)).card

/-- `ϱ(x, y) = Ψ(x, y)/x`. -/
def rho (x y : ℝ) : ℝ := (Psi x y : ℝ) / x

/-! ## §10. The N-dependent quantities ([pp] §sec:params) -/

namespace Params

variable (P : Params) (a : ℝ) (N : ℕ)

/-- `L = L_N[1/2, 1]`. -/
def Lv : ℝ := LN (1 / 2) 1 N

/-- `B = ⌊L^a⌋`, the smoothness and trial-division bound. -/
def Bsm : ℕ := ⌊(Lv N) ^ a⌋₊

/-- `X = ⌊√(N/2)⌋`. -/
def Xsp : ℕ := ⌊Real.sqrt ((N : ℝ) / 2)⌋₊

/-- `j`, maximal with `2^{kj} ≤ X`. -/
def jJ : ℕ := Nat.log (2 ^ P.k) (Xsp N)

/-- `Y = 2^{j−1}` (for `j ≥ 1`, i.e. all large `N`; at `j = 0` natural subtraction gives `Y = 1`
instead of the paper's `1/2`, which only affects tiny `N` where SPLIT's guards return FAIL). -/
def Yj : ℕ := 2 ^ (P.jJ N - 1)

/-- `J = {ℓ prime : Y < ℓ ≤ 2Y}`. -/
def Jset : Finset ℕ := primesIoc (P.Yj N) (2 * P.Yj N)

/-- `S = ⌈L^{θ/(4a)}⌉`. -/
def Ssz : ℕ := ⌈(Lv N) ^ ((P.θ : ℝ) / (4 * a))⌉₊

/-- `V_ℓ = ⌊√ℓ · S⌋`. -/
def Vl (ℓ : ℕ) : ℕ := ⌊Real.sqrt ℓ * (P.Ssz a N : ℝ)⌋₊

/-- `V_max = ⌊√(2Y) · S⌋`. -/
def Vmax : ℕ := ⌊Real.sqrt (2 * (P.Yj N : ℝ)) * (P.Ssz a N : ℝ)⌋₊

/-- `M = ⌊N^μ⌋`. -/
def Mb : ℕ := ⌊(N : ℝ) ^ (P.μ : ℝ)⌋₊

/-- `H = ⌊N^h⌋`. -/
def Hb : ℕ := ⌊(N : ℝ) ^ (P.h : ℝ)⌋₊

/-- `W = ⌈N^{μ+ν+η₀/2}⌉`. -/
def Wb : ℕ := ⌈(N : ℝ) ^ ((P.μ + P.ν + P.η₀ / 2 : ℚ) : ℝ)⌉₊

/-- `X_F = C′_d M H V_max^d` with `C′_d = 2^{2d+5} + 2d`. -/
def XF : ℝ := ((2 : ℝ) ^ (2 * P.d + 5) + 2 * P.d) * P.Mb N * P.Hb N * (P.Vmax a N : ℝ) ^ P.d

/-- `X_r = 2 M V_max / Y`. -/
def Xr : ℝ := 2 * (P.Mb N : ℝ) * P.Vmax a N / P.Yj N

/-- `ρ_r = ϱ(X_r, B)`. -/
def rhoR : ℝ := rho (P.Xr a N) (Bsm a N)

/-- `ρ_F = ϱ(X_F, B)`. -/
def rhoF : ℝ := rho (P.XF a N) (Bsm a N)

/-! ## §11. Inputs of SPLIT, the field law, candidates and the filtered count -/

/-- The integers handed to SPLIT ([pp] §sec:K, §sec:uniform): `N > 1` odd, with no prime factor
`≤ B`, and not a prime power. -/
def SplitInput : Prop :=
  1 < N ∧ Odd N ∧ (∀ p : ℕ, p.Prime → p ∣ N → Bsm a N < p) ∧ ¬ IsPrimePow N

/-- A field `A = ((g₀, g₁), w, c)`. -/
abbrev Field (P : Params) : Type := (ℤ × ℤ) × ℕ × (Fin P.d → ℤ)

/-- `𝒢 = {(g₀, g₁) ∈ ℤ² ∩ [M/2, M]² : gcd(g₀, g₁) = 1}`. -/
def gSet : Finset (ℤ × ℤ) :=
  ((Icc (((P.Mb N : ℤ) + 1) / 2) (P.Mb N : ℤ)) ×ˢ (Icc (((P.Mb N : ℤ) + 1) / 2) (P.Mb N : ℤ))).filter
    (fun g => Int.gcd g.1 g.2 = 1)

/-- The raw field law `𝒜` ([pp] Definition def:field): `g` uniform on `𝒢`, `w` uniform among the
primes in `(W, 2W]`, `c` uniform on `C_H`, independent; i.e. the uniform law on this product. -/
def fields : Finset P.Field :=
  P.gSet N ×ˢ (primesIoc (P.Wb N) (2 * P.Wb N) ×ˢ cube P.d (P.Hb N))

/-- `𝒫_ℓ = {(a, b) ∈ ℤ² : b prime, V_ℓ/2 < b ≤ V_ℓ, |a| ≤ V_ℓ, gcd(a, b) = 1}`. -/
def cands (ℓ : ℕ) : Finset (ℤ × ℤ) :=
  ((Icc (-(P.Vl a N ℓ : ℤ)) (P.Vl a N ℓ : ℤ)) ×ˢ (Icc 1 (P.Vl a N ℓ : ℤ))).filter
    (fun Q => (Q.2.toNat).Prime ∧ (P.Vl a N ℓ : ℤ) < 2 * Q.2 ∧ Int.gcd Q.1 Q.2 = 1)

/-- The filtered count `N_ℓ(A)` ([pp] §sec:candidates): the number of `Q ∈ 𝒫_ℓ ∩ Λ_ℓ(g)` with
`M V_ℓ/32 ≤ |G(Q)|`, `G(Q)/ℓ` `B`-smooth and `F_A(Q)` `B`-smooth (nonzero); `0` if `ℓ ∣ N`.
Here `F_A = F̂_{g,w} + G·R_c` with `F̂_{g,w}` Lemma GB applied to `N w`. -/
def countN (A : P.Field) (ℓ : ℕ) : ℕ :=
  if ℓ ∣ N then 0 else
    ((P.cands a N ℓ).filter (fun Q =>
      (ℓ : ℤ) ∣ Gform A.1.1 A.1.2 Q.1 Q.2 ∧
      ((P.Mb N : ℝ) * (P.Vl a N ℓ : ℝ)) / 32 ≤ |(Gform A.1.1 A.1.2 Q.1 Q.2 : ℝ)| ∧
      IsSmooth (Bsm a N) (Gform A.1.1 A.1.2 Q.1 Q.2 / ℓ) ∧
      IsSmooth (Bsm a N) (FA A.1.1 A.1.2 ((N : ℤ) * (A.2.1 : ℤ)) P.d A.2.2 Q.1 Q.2))).card

/-- `γ = Pr_{A ∼ 𝒜, ℓ ∈ J}[N_ℓ(A) ≥ 1]` ([pp] §sec:candidates). -/
def gammaD : ℝ := prob (P.fields N ×ˢ P.Jset N) (fun x => 1 ≤ P.countN a N x.1 x.2)

/-- The descendable primes `𝒟_A = {ℓ ∈ J : N_ℓ(A) ≥ 1}`. -/
def descPrimes (A : P.Field) : Finset ℕ := (P.Jset N).filter (fun ℓ => 1 ≤ P.countN a N A ℓ)

/-- The probability, over `r` uniform on `[1, N−1]`, of a complete descent in the field `A`: `r` is a
unit and every shortest vector `(u₀, v₀)` of `Λ_{r²}` has `u₀v₀ ≠ 0` and both halves products of at
most `k` distinct descendable primes ([pp] Theorem thm:K). -/
def descentProb (A : P.Field) : ℝ :=
  prob (Icc 1 (N - 1)) (fun r : ℕ => IsUnit (r : ZMod N) ∧
    ∀ x y : ℤ, IsShortestVec N ((r : ZMod N) ^ 2) x y →
      x * y ≠ 0 ∧ IsProdAtMostK (P.descPrimes a N A) P.k x.natAbs ∧
        IsProdAtMostK (P.descPrimes a N A) P.k y.natAbs)

end Params

end FactoringA
