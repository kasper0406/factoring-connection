import LeanA.RelationSupply
import LeanA.Recursion
import LeanA.Exponents

/-!
# (h) Cost assembly: Theorems A, RA, Corollary P, and the assembly steps of Theorems U and DL

## The cost-accounting model

Lean has no running-time semantics for an NFS-type algorithm, and an existence statement
"there is an algorithm with expected time ≤ …" would be vacuous for Lean functions. So the
algorithmic layer is an explicit *cost-accounting model*: a record of real-valued functions of `N`
(expected costs, probabilities, row targets) for the fixed algorithm of [pp] §sec:alg, together
with an *interface* — a list of inequalities between them. Each inequality is one displayed claim
of the paper about that algorithm (cited at the field). The interface never mentions the final
exponent; every bound on a component is either polylogarithmic bookkeeping, a published input
(Lenstra–Pomerance), or the paper's own component statement (Lemma lem:reject, Prop. prop:select(ii),
Prop. prop:BS, Prop. prop:lv6, §sec:uniform). The probability `pTrial` is tied to the concrete
descent probability of `DefsNFS.lean` (interface field `trial_ge`), so the relation supply enters
through the proved theorem `descent_supply`, not through an assumption.

What the Lean theorems certify: given the interface, the published inputs and the two internal
moment hypotheses, the expected cost of FACTOR is `≤ L_N[1/2, c + o(1)]` with the paper's `c`.
What they do NOT certify: that the real algorithm satisfies the interface (that is the content of
[pp] §§sec:alg, sec:lv6, sec:batch and the factor-tree accounting, which are not formalised).

The main statements (MAIN), the model structures and the interface definitions are copied into
`Challenge.lean`; do not change them.
-/

noncomputable section
open Classical

namespace FactoringA

open Real Finset Filter

/-- PUBLISHED INPUT (Lenstra–Pomerance). H. W. Lenstra, Jr. and C. Pomerance, *A rigorous time
bound for factoring integers*, J. Amer. Math. Soc. 5 (1992), 483–516, Theorem 1.1 ([pp] Fact
fact:LP): their Las Vegas algorithm finds the complete prime factorisation of every `n > 1` in
expected time `≤ L_n[1/2, 1 + o(1)]`. Stated for `T n` = the maximum over inputs `≤ n` of that
expected time, the function used by a cost model below. -/
def LenstraPomerance (T : ℕ → ℝ) : Prop := LBound (1 / 2) 1 T

/-- The maximum of `f` over the split inputs `n` with `N₀ ≤ n ≤ x` (`0` if there are none). -/
def supSplit (f : ℕ → ℝ) (a N₀ x : ℝ) : ℝ :=
  sSup (f '' {n : ℕ | N₀ ≤ n ∧ (n : ℝ) ≤ x ∧ Params.SplitInput a n})

/-! ## Theorem A -/

/-- Cost-accounting model of FACTOR with the Lenstra–Pomerance shape test ([pp] §sec:alg). All
entries are functions of the input `N`. -/
structure CostModelA (P : Params) where
  /-- expected running time of FACTOR(N) -/
  factorCost : ℕ → ℝ
  /-- expected running time of SPLIT(N) -/
  splitCost : ℕ → ℝ
  /-- probability that one trial (SAMPLE_FIELD, then DESCEND_ONCE) ends in a row or a factor exit -/
  pTrial : ℕ → ℝ
  /-- the rejection event `𝓡` of [pp] Lemma lem:reject, as a set of raw fields -/
  rej : (N : ℕ) → Finset P.Field
  /-- uniform bound on the expected cost of one trial and of one DESCEND_ONCE call -/
  Ctrial : ℕ → ℝ
  /-- uniform bound on the expected cost of POSTPROCESS -/
  Tpost : ℕ → ℝ
  /-- the row target `M_rows` -/
  Mrows : ℕ → ℝ
  /-- expected cost of the shape test (two Lenstra–Pomerance runs on `|u₀|, |v₀| ≤ √(2N)` and AKS) -/
  Tshape : ℕ → ℝ
  /-- expected cost of one SQ_SEARCH -/
  Tsq : ℕ → ℝ
  /-- maximum over inputs `≤ n` of the Lenstra–Pomerance expected time -/
  TLP : ℕ → ℝ

/-- The interface of the cost model at smoothness exponent `a` (constants `C`, threshold `N₀`). -/
structure CostInterfaceA (P : Params) (a : ℝ) (M : CostModelA P) : Prop where
  /-- One constant and one threshold for all the bookkeeping below. -/
  bounds : ∃ C N₀ : ℝ, 0 ≤ C ∧
    (∀ N : ℕ, N₀ ≤ N → Params.SplitInput a N →
      -- [pp] Prop. prop:select(i) proof / Theorem thm:K proof: a field outside `𝓡` either exits in
      -- SAMPLE_FIELD or is returned, and on the complete-descent event DESCEND_ONCE returns a row
      -- or exits (Prop. prop:BS); so a trial succeeds with at least this probability.
      expect (P.fields N) (fun A => if A ∈ M.rej N then 0 else P.descentProb a N A) ≤ M.pTrial N ∧
      -- [pp] Lemma lem:reject.
      M.rej N ⊆ P.fields N ∧
      prob (P.fields N) (fun A => A ∈ M.rej N) ≤
        Real.exp (-(Real.log (P.Hb N) / (13 * P.d * Real.log (Real.log (P.Hb N))))) +
          2 / (2 * (P.Hb N : ℝ) + 1) + 1 / (P.Wb N : ℝ) ∧
      -- [pp] Prop. prop:select(ii).
      M.splitCost N ≤ 3 * (M.Mrows N * M.Ctrial N + M.Tpost N) / M.pTrial N ∧
      -- [pp] eq:Mrows: `M_rows = C_col + 12 ≤ (d+1)π(B) + O_d(log N)`.
      1 ≤ M.Mrows N ∧ M.Mrows N ≤ (P.d + 1) * (Nat.primeCounting (Params.Bsm a N) : ℝ) + C * Real.log N ∧
      -- [pp] Prop. prop:cost, "One trial": shape test, at most 2k searches, polylog bookkeeping.
      0 ≤ M.Ctrial N ∧ M.Ctrial N ≤ M.Tshape N + 2 * P.k * M.Tsq N + (Real.log N) ^ C ∧
      M.Tshape N ≤ 2 * M.TLP ⌊Real.sqrt (2 * N)⌋₊ + (Real.log N) ^ C ∧
      -- [pp] Prop. prop:BS: `Õ(B + S² log N)`.
      M.Tsq N ≤ (Real.log N) ^ C * ((Params.Bsm a N : ℝ) + (P.Ssz a N : ℝ) ^ 2 * Real.log N) ∧
      0 ≤ M.Tpost N) ∧
    -- [pp] Prop. prop:lv6: `T_post ≤ M_rows² (log N)^{O(1)} + B^{1+o(1)}`.
    (∀ ε > 0, ∃ N₁ : ℝ, ∀ N : ℕ, N₁ ≤ N → Params.SplitInput a N →
      M.Tpost N ≤ M.Mrows N ^ 2 * (Real.log N) ^ C + (Params.Bsm a N : ℝ) ^ (1 + ε)) ∧
    -- [pp] §sec:uniform (proof of Theorem thm:main): at most `log₂N` calls of SPLIT, each raced
    -- against trial division (constant cost below `N₀`), on inputs `≤ N`, plus preprocessing
    -- (primality, perfect powers, trial division to `B`) at every node.
    (∀ N : ℕ, 2 ≤ N → M.factorCost N ≤ C + (Real.logb 2 N + 1) *
      (C * (Params.Bsm a N : ℝ) * (Real.log N) ^ C + 2 * supSplit M.splitCost a N₀ N))

/-- MAIN. [pp] Theorem thm:main (Theorem A), conditional form. For admissible `(d, k, μ, η₀)`, with
`a = θ/√2`: from the published inputs (CEP, Rosser–Schoenfeld, Lenstra–Pomerance), the internal
moment hypotheses FM and LD*, and the cost interface, the expected time of FACTOR is
`≤ L_N[1/2, (1+θ)/√2 + o(1)]`. -/
theorem theoremA (P : Params) (hP : P.Admissible) (hCEP : CEP) (hRS : RosserSchoenfeld)
    (hFM : FirstMoment P ((P.θ : ℝ) / Real.sqrt 2))
    (hLD : SecondMoment P ((P.θ : ℝ) / Real.sqrt 2))
    (M : CostModelA P) (hM : CostInterfaceA P ((P.θ : ℝ) / Real.sqrt 2) M)
    (hLP : LenstraPomerance M.TLP) :
    LBound (1 / 2) P.cA M.factorCost := by
  sorry

/-- MAIN. [pp] Corollary cor:main, conditional form: if the hypotheses of `theoremA` hold for every
admissible tuple, then for every `c > 1/√2` some admissible tuple's FACTOR runs in expected time
`O(L_N[1/2, c])`. -/
theorem corollaryA (hCEP : CEP) (hRS : RosserSchoenfeld)
    (hFM : ∀ P : Params, P.Admissible → FirstMoment P ((P.θ : ℝ) / Real.sqrt 2))
    (hLD : ∀ P : Params, P.Admissible → SecondMoment P ((P.θ : ℝ) / Real.sqrt 2))
    (M : ∀ P : Params, CostModelA P)
    (hM : ∀ P : Params, P.Admissible → CostInterfaceA P ((P.θ : ℝ) / Real.sqrt 2) (M P))
    (hLP : ∀ P : Params, LenstraPomerance (M P).TLP) (c : ℝ) (hc : 1 / Real.sqrt 2 < c) :
    ∃ P : Params, P.Admissible ∧
      ∃ K : ℝ, ∀ N : ℕ, 2 ≤ N → (M P).factorCost N ≤ K * LN (1 / 2) c N := by
  sorry

/-! ## Theorem RA (recursive shape test) -/

/-- Cost-accounting model of FACTOR_rec[a] ([pp] §sec:recursive). -/
structure CostModelRA (P : Params) where
  /-- `F(x) = sup_{2 ≤ n ≤ x} E[cost of FACTOR_rec[a](n)]` -/
  Fsup : ℝ → ℝ
  /-- expected number of DESCEND_ONCE calls in SPLIT^O(N) -/
  tau : ℕ → ℝ
  /-- expected work of SPLIT^O(N) outside the calls of the shape-test oracle `O` -/
  nonO : ℕ → ℝ
  /-- as in `CostModelA` -/
  pTrial : ℕ → ℝ
  rej : (N : ℕ) → Finset P.Field
  Mrows : ℕ → ℝ
  Tsq : ℕ → ℝ
  Tpost : ℕ → ℝ

/-- The interface of the recursive model at smoothness exponent `a`. -/
structure CostInterfaceRA (P : Params) (a : ℝ) (M : CostModelRA P) : Prop where
  bounds : ∃ C N₀ : ℝ, 0 ≤ C ∧
    (∀ N : ℕ, N₀ ≤ N → Params.SplitInput a N →
      -- as in `CostInterfaceA`
      expect (P.fields N) (fun A => if A ∈ M.rej N then 0 else P.descentProb a N A) ≤ M.pTrial N ∧
      M.rej N ⊆ P.fields N ∧
      prob (P.fields N) (fun A => A ∈ M.rej N) ≤
        Real.exp (-(Real.log (P.Hb N) / (13 * P.d * Real.log (Real.log (P.Hb N))))) +
          2 / (2 * (P.Hb N : ℝ) + 1) + 1 / (P.Wb N : ℝ) ∧
      1 ≤ M.Mrows N ∧ M.Mrows N ≤ (P.d + 1) * (Nat.primeCounting (Params.Bsm a N) : ℝ) + C * Real.log N ∧
      M.Tsq N ≤ (Real.log N) ^ C * ((Params.Bsm a N : ℝ) + (P.Ssz a N : ℝ) ^ 2 * Real.log N) ∧
      0 ≤ M.Tsq N ∧ 0 ≤ M.Tpost N ∧ 0 ≤ M.tau N ∧
      -- [pp] Lemma lem:generala (A1): Prop. prop:select(ii) counting DESCEND_ONCE calls.
      M.tau N ≤ 3 * M.Mrows N / M.pTrial N ∧
      -- [pp] Lemma lem:generala (A2): Prop. prop:select(ii) with the oracle's cost removed and at
      -- most 2k searches charged to every DESCEND_ONCE call.
      M.nonO N ≤ 3 * (M.Mrows N * (2 * P.k * M.Tsq N + (Real.log N) ^ C) + M.Tpost N) / M.pTrial N) ∧
    -- [pp] Prop. prop:lv6, as in `CostInterfaceA`
    (∀ ε > 0, ∃ N₁ : ℝ, ∀ N : ℕ, N₁ ≤ N → Params.SplitInput a N →
      M.Tpost N ≤ M.Mrows N ^ 2 * (Real.log N) ^ C + (Params.Bsm a N : ℝ) ^ (1 + ε)) ∧
    -- `F` is finite, nonnegative, nondecreasing ([pp] Theorem thm:RA, "Finiteness").
    (∀ x, 2 ≤ x → 0 ≤ M.Fsup x) ∧ MonotoneOn M.Fsup (Set.Ici 2) ∧
    -- [pp] Theorem thm:RA, "Recursion": factor tree (< 2 log₂x splitting nodes), race against trial
    -- division (factor 2), Lemma OS, Lemma lem:waldrec, and at most two oracle calls per
    -- DESCEND_ONCE call on inputs `≤ √(2x)`.
    (∀ x : ℝ, N₀ ≤ x → M.Fsup x ≤ 2 * Real.logb 2 x *
      (C + C * (Params.Bsm a ⌊x⌋₊ : ℝ) * (Real.log x) ^ C + 2 * supSplit M.nonO a N₀ x +
        2 * (2 * supSplit M.tau a N₀ x) * M.Fsup (Real.sqrt (2 * x))))

/-- MAIN. [pp] Theorem thm:RA (Theorem RA), conditional form. For admissible `(d, k, μ, η₀)` and any
fixed `a > 0`, the expected time of FACTOR_rec[a] is `≤ L_N[1/2, c(a) + o(1)]`,
`c(a) = max(a + θ/(2a), (2+√2)a)`. -/
theorem theoremRA (P : Params) (hP : P.Admissible) (a : ℝ) (ha : 0 < a) (hCEP : CEP)
    (hRS : RosserSchoenfeld) (hFM : FirstMoment P a) (hLD : SecondMoment P a)
    (M : CostModelRA P) (hM : CostInterfaceRA P a M) :
    LBound (1 / 2) (costExpRA (P.θ : ℝ) a) (fun N => M.Fsup N) := by
  sorry

/-- MAIN. [pp] Theorem thm:RA at the optimal `a* = √((√2−1)θ/2)`: expected time
`≤ L_N[1/2, √((1+√2)θ) + o(1)]`. -/
theorem theoremRA_opt (P : Params) (hP : P.Admissible) (hCEP : CEP) (hRS : RosserSchoenfeld)
    (hFM : FirstMoment P (Real.sqrt ((Real.sqrt 2 - 1) * P.θ / 2)))
    (hLD : SecondMoment P (Real.sqrt ((Real.sqrt 2 - 1) * P.θ / 2)))
    (M : CostModelRA P) (hM : CostInterfaceRA P (Real.sqrt ((Real.sqrt 2 - 1) * P.θ / 2)) M) :
    LBound (1 / 2) P.cRAstar (fun N => M.Fsup N) := by
  sorry

/-- MAIN. [pp] Corollary cor:P (portfolio), conditional form. Suppose the hypotheses of
`theoremRA_opt` hold for the portfolio tuples `portfolioParams j` (`j ≥ 1`), and the restart
schedule's expected time `T` satisfies, for each `j ≥ 1`, `T(N) ≤ K_j (log N)^{K_j} F_j(N)`
eventually ([pp] proof of Corollary cor:P: Markov's inequality per phase). Then
`T(N) ≤ L_N[1/2, o(1)]`. -/
theorem corollaryP (hCEP : CEP) (hRS : RosserSchoenfeld)
    (hFM : ∀ j : ℕ, 1 ≤ j → FirstMoment (portfolioParams j)
      (Real.sqrt ((Real.sqrt 2 - 1) * (portfolioParams j).θ / 2)))
    (hLD : ∀ j : ℕ, 1 ≤ j → SecondMoment (portfolioParams j)
      (Real.sqrt ((Real.sqrt 2 - 1) * (portfolioParams j).θ / 2)))
    (M : ∀ j : ℕ, CostModelRA (portfolioParams j))
    (hM : ∀ j : ℕ, 1 ≤ j → CostInterfaceRA (portfolioParams j)
      (Real.sqrt ((Real.sqrt 2 - 1) * (portfolioParams j).θ / 2)) (M j))
    (T : ℕ → ℝ)
    (hport : ∀ j : ℕ, 1 ≤ j → ∃ K : ℝ, 0 ≤ K ∧
      ∀ᶠ N : ℕ in atTop, T N ≤ K * (Real.log N) ^ K * (M j).Fsup N) :
    LBound (1 / 2) 0 T := by
  sorry

/-! ## Theorem U and Theorem DL: assembly steps only -/

/-- MAIN. lv5unif §10 (Theorem U), assembly step: the recursion inequality (10.1),
`F(x) ≤ 2 log₂x [P(x) + D(x) F(√(2x))]` with `P(x) ≤ L_x[2/5, max(2a, a+s+max(a,σ)) + o(1)]` and
`D(x) ≤ L_x[2/5, a+s+o(1)]`, gives `F(x) ≤ L_x[2/5, c + o(1)]` with
`c = max(2a, a+s+max(a,σ), (a+s)/(1−2^{−2/5}))` ((10.3)). The inputs (10.1) are INTERNAL and are not
formalised here. -/
theorem theoremU_assembly (F Pc D : ℝ → ℝ) (a s σ : ℝ) (ha : 0 < a) (hs : 0 ≤ s) (hσ : 0 ≤ σ)
    (hF0 : ∀ x, 2 ≤ x → 0 ≤ F x) (hmono : MonotoneOn F (Set.Ici 2))
    (hP : ∀ η > 0, ∀ᶠ x in atTop, Pc x ≤ LN (2 / 5) (max (2 * a) (a + s + max a σ) + η) x)
    (hD : ∀ η > 0, ∀ᶠ x in atTop, D x ≤ LN (2 / 5) (a + s + η) x)
    (hD0 : ∀ᶠ x in atTop, 0 ≤ D x)
    (hrec : ∀ᶠ x in atTop, F x ≤ 2 * Real.logb 2 x * (Pc x + D x * F (Real.sqrt (2 * x)))) :
    LBound (2 / 5) (max (max (2 * a) (a + s + max a σ)) ((a + s) / (1 - 2 ^ (-(2 / 5 : ℝ)))))
      (fun N => F N) := by
  sorry

/-- MAIN. [pp] Theorem thm:DL, assembly step: with `B = L^a`, if the expected time is
`≤ (M_rows C_trial + T_post)/p` with `M_rows ≤ L^{a+o(1)}`, `C_trial ≤ L^{max(1/√2, θ/(2a))+o(1)}`,
`T_post ≤ L^{3a+o(1)}` and `p ≥ L^{−o(1)}` (all INTERNAL inputs, not formalised here), then the
expected time is `≤ L_p[1/2, max(a + 1/√2, a + θ/(2a), 3a) + o(1)]`; at the optimal `a` this is
`c_DL(θ)` (`costExpDL_min`). -/
theorem theoremDL_assembly (θ a : ℝ) (ha : 0 < a) (T Mrows Ctrial Tpost pT : ℕ → ℝ)
    (hT : ∀ᶠ N : ℕ in atTop, T N ≤ (Mrows N * Ctrial N + Tpost N) / pT N)
    (hM : LBound (1 / 2) a Mrows) (hM0 : ∀ᶠ N : ℕ in atTop, 0 ≤ Mrows N)
    (hC : LBound (1 / 2) (max (1 / Real.sqrt 2) (θ / (2 * a))) Ctrial)
    (hC0 : ∀ᶠ N : ℕ in atTop, 0 ≤ Ctrial N)
    (hTp : LBound (1 / 2) (3 * a) Tpost)
    (hp : ∀ ε > 0, ∀ᶠ N : ℕ in atTop, LN (1 / 2) (-ε) N ≤ pT N) :
    LBound (1 / 2) (costExpDL θ a) T := by
  sorry

end FactoringA
