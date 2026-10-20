import LeanA.DefsNFS
import LeanA.SplitCounting
import LeanA.Moments

/-!
# (d) Relation supply: Theorem D_k ([pp] Theorem thm:Dk) and field selection (i)
([pp] Prop. prop:select(i)), from the moment hypotheses

* PUBLISHED INPUT: `CEP` (Canfield–Erdős–Pomerance 1983).
* INTERNAL HYPOTHESES: `FirstMoment` ([pp] Theorem thm:FM) and `SecondMoment` ([pp] Theorem thm:LD,
  "LD*"). They are statements about the concrete filtered count `Params.countN` of `DefsNFS.lean`.
* PROVED here: `omega_le` ([pp] Lemma lem:omega), `relation_supply` (γ ≥ L^{−o(1)}), and
  `descent_supply` (a non-rejected raw field followed by a uniform root gives a complete descent with
  probability `L^{−o(1)}`).

The main statements (MAIN) and the hypothesis definitions are copied into `Challenge.lean`;
do not change them.
-/

noncomputable section
open Classical

namespace FactoringA

open Real Finset Filter

/-- PUBLISHED INPUT (Canfield–Erdős–Pomerance). E. R. Canfield, P. Erdős and C. Pomerance, *On a
problem of Oppenheim concerning "factorisatio numerorum"*, J. Number Theory 17 (1983), 1–28,
Corollary to Theorem 3.1 ([pp] Fact fact:CEP): for each `ε > 0`, uniformly for `u → ∞` with
`u ≤ (1 − ε) log x / log log x`, `Ψ(x, x^{1/u}) = x u^{−u(1+o(1))}`. -/
def CEP : Prop :=
  ∀ ε > 0, ∀ δ > 0, ∃ U₀ : ℝ, ∀ x u : ℝ, U₀ ≤ u → u ≤ (1 - ε) * Real.log x / Real.log (Real.log x) →
    x * u ^ (-(u * (1 + δ))) ≤ (Psi x (x ^ (1 / u)) : ℝ) ∧
      (Psi x (x ^ (1 / u)) : ℝ) ≤ x * u ^ (-(u * (1 - δ)))

/-- INTERNAL HYPOTHESIS ([pp] Theorem thm:FM, first moment; with [pp] Lemma lem:generala for general
`a`): for split inputs `N ≥ N₀`, `E_{ℓ ∈ J} E_{A ∼ 𝒜} N_ℓ(A) ≥ S² ρ_r ρ_F / (log N)^C`. -/
def FirstMoment (P : Params) (a : ℝ) : Prop :=
  ∃ C N₀ : ℝ, ∀ N : ℕ, N₀ ≤ N → Params.SplitInput a N →
    (P.Ssz a N : ℝ) ^ 2 * P.rhoR a N * P.rhoF a N / (Real.log N) ^ C ≤
      expect (P.fields N ×ˢ P.Jset N) (fun x => (P.countN a N x.1 x.2 : ℝ))

/-- INTERNAL HYPOTHESIS ([pp] Theorem thm:LD, "LD*", second moment; with [pp] Lemma lem:generala):
for split inputs `N ≥ N₀`, `E_{ℓ,A}[N_ℓ(A)(N_ℓ(A) − 1)] ≤ (log N)^C S⁴ ρ_r² ρ_F²`. -/
def SecondMoment (P : Params) (a : ℝ) : Prop :=
  ∃ C N₀ : ℝ, ∀ N : ℕ, N₀ ≤ N → Params.SplitInput a N →
    expect (P.fields N ×ˢ P.Jset N)
        (fun x => (P.countN a N x.1 x.2 : ℝ) * ((P.countN a N x.1 x.2 : ℝ) - 1)) ≤
      (Real.log N) ^ C * (P.Ssz a N : ℝ) ^ 4 * (P.rhoR a N) ^ 2 * (P.rhoF a N) ^ 2

/-- MAIN. [pp] Lemma lem:omega: if `N > 1` has no prime factor `≤ B` (`B ≥ 2`), then
`ω(N) log B ≤ log N`. -/
theorem omega_le (N B : ℕ) (hN : 1 < N) (hB : 2 ≤ B) (hNB : ∀ p : ℕ, p.Prime → p ∣ N → B < p) :
    (omega N : ℝ) * Real.log B ≤ Real.log N := by
  sorry

/-- MAIN. [pp] Theorem thm:Dk (relation supply): from the two moment hypotheses and CEP,
`γ ≥ L^{−o(1)}` uniformly over split inputs. -/
theorem relation_supply (P : Params) (hP : P.Admissible) (a : ℝ) (ha : 0 < a) (hCEP : CEP)
    (hFM : FirstMoment P a) (hLD : SecondMoment P a) :
    ∀ ε > 0, ∃ N₀ : ℕ, ∀ N : ℕ, N₀ ≤ N → Params.SplitInput a N →
      LN (1 / 2) (-ε) N ≤ P.gammaD a N := by
  sorry

/-- MAIN. [pp] Prop. prop:select(i) with Theorem thm:K, concrete form: for any set `𝓡` of fields
with `Pr_𝒜(𝓡) ≤ γ/2`, the probability that a raw field is not in `𝓡` and a uniform `r ∈ [1, N−1]`
gives a complete descent in it is `≥ L^{−o(1)}`, uniformly over split inputs. -/
theorem descent_supply (P : Params) (hP : P.Admissible) (a : ℝ) (ha : 0 < a) (hCEP : CEP)
    (hRS : RosserSchoenfeld) (hFM : FirstMoment P a) (hLD : SecondMoment P a) :
    ∀ ε > 0, ∃ N₀ : ℕ, ∀ N : ℕ, N₀ ≤ N → Params.SplitInput a N →
      ∀ Rej : Finset P.Field, Rej ⊆ P.fields N →
        prob (P.fields N) (fun A => A ∈ Rej) ≤ P.gammaD a N / 2 →
          LN (1 / 2) (-ε) N ≤
            expect (P.fields N) (fun A => if A ∈ Rej then 0 else P.descentProb a N A) := by
  sorry

end FactoringA
