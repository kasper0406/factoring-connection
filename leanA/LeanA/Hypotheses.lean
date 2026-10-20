import LeanA.DefsNFS

/-!
# Trusted definitions (part 3): the hypotheses and the cost models

Every analytic or algorithmic input of the formalisation, as a `Prop` (passed to theorems as an
explicit hypothesis), and the cost-model records. Docstrings label each one:
* PUBLISHED INPUT: a transcription of a published theorem (citation in the docstring);
* INTERNAL HYPOTHESIS: a claim of the project's own (unrefereed) paper that is not proved in Lean;
* TRUSTED: a plain definition used in statements.

`Challenge.lean` imports this module (together with `Defs` and `DefsNFS`), so this file is part of
the trusted base and must be reviewed as such.
-/

set_option autoImplicit false

noncomputable section
open Classical

namespace FactoringA

section
open Finset Real

/-! ### (from Moments.lean: inputs of Lemma W) -/

/-- PUBLISHED INPUT (Brun–Titchmarsh). H. L. Montgomery and R. C. Vaughan, *The large sieve*,
Mathematika 20 (1973), 119–134, Theorem 2: for positive reals `x, y` and coprime positive integers
`q, a` with `y > q`, `π(x+y; q, a) − π(x; q, a) < 2y/(φ(q) log(y/q))`. (The strict form is also
implied by T. Yamada, arXiv:2312.16090, which has `log(y/q) + 0.8601` in the denominator.) -/
def BrunTitchmarsh : Prop :=
  ∀ (x y : ℝ) (q a : ℕ), 0 < x → 0 < y → 1 ≤ q → 1 ≤ a → Nat.Coprime a q → (q : ℝ) < y →
    (primeCountAP (x + y) q a : ℝ) - primeCountAP x q a <
      2 * y / ((Nat.totient q : ℝ) * Real.log (y / q))

/-- PUBLISHED INPUT (Rosser–Schoenfeld). J. B. Rosser and L. Schoenfeld, *Approximate formulas
for some functions of prime numbers*, Illinois J. Math. 6 (1962), 64–94, p. 69: (3.6) of Corollary 1
after Theorem 2 (`π(z) < 1.25506 z/log z` for `1 < z`) and (3.8) of Corollary 3
(`3z/(5 log z) < π(2z) − π(z)` for `20½ ≤ z`). -/
def RosserSchoenfeld : Prop :=
  (∀ z : ℝ, 1 < z → (Nat.primeCounting ⌊z⌋₊ : ℝ) < 1.25506 * z / Real.log z) ∧
    (∀ z : ℝ, 20.5 ≤ z →
      3 * z / (5 * Real.log z) < (Nat.primeCounting ⌊2 * z⌋₊ : ℝ) - Nat.primeCounting ⌊z⌋₊)

end

section
open Real Finset Filter

/-! ### (from RelationSupply.lean: CEP and the moment hypotheses) -/

/-- PUBLISHED INPUT (Canfield–Erdős–Pomerance). E. R. Canfield, P. Erdős and C. Pomerance, *On a
problem of Oppenheim concerning "factorisatio numerorum"*, J. Number Theory 17 (1983), 1–28,
Corollary on p. 15 (lower half: their Theorem 3.1; upper half: de Bruijn, Indag. Math. 28 (1966),
Part II, Thm 2): if `ε > 0` and `3 ≤ u ≤ (1 − ε) log x / log log x`, then
`Ψ(x, x^{1/u}) = x·exp{−u(log u + log log u − 1 + (log log u − 1)/log u + E(x, u))}` with
`|E(x, u)| ≤ c_ε (log log u)²/(log u)²`. In particular ([pp] Fact fact:CEP)
`Ψ(x, x^{1/u}) = x u^{−u(1+o(1))}` uniformly as `u → ∞` in that range. We state this consequence,
for `x ≥ 16` (so that `log log x > 0`; the source tacitly assumes `x` large — without a lower bound
on `x` the Lean statement would be false, e.g. for negative `x` or `x` just above `e`). -/
def CEP : Prop :=
  ∀ ε > 0, ∀ δ > 0, ∃ U₀ : ℝ, ∀ x u : ℝ, 16 ≤ x → U₀ ≤ u →
    u ≤ (1 - ε) * Real.log x / Real.log (Real.log x) →
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

end

section
open Real Finset Filter

/-! ### (from MainTheorems.lean: cost models) -/

/-- PUBLISHED INPUT (Lenstra–Pomerance). H. W. Lenstra, Jr. and C. Pomerance, *A rigorous time
bound for factoring integers*, J. Amer. Math. Soc. 5 (1992), 483–516, Theorem 1.1 ([pp] Fact
fact:LP; the unnumbered Theorem of §1, p. 483, proved in §10): there is a probabilistic algorithm that
factors any positive integer `n` completely and takes expected time at most `L_n[1/2, 1 + o(1)]`.
Stated for `T n` = the maximum over inputs `≤ n` of that expected time. As a `Prop` about an arbitrary
`T` this is only meaningful together with the MODELLING assumption that a cost model's `TLP` is that
function. -/
def LenstraPomerance (T : ℕ → ℝ) : Prop := LBound (1 / 2) 1 T

/-- The maximum of `f` over the split inputs `n` with `N₀ ≤ n ≤ x` (`0` if there are none). -/
def supSplit (f : ℕ → ℝ) (a N₀ x : ℝ) : ℝ :=
  sSup (f '' {n : ℕ | N₀ ≤ n ∧ (n : ℝ) ≤ x ∧ Params.SplitInput a n})

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

/-- INTERNAL HYPOTHESIS (algorithmic interface for Theorem A) at smoothness exponent `a` (constants `C`,
threshold `N₀`). Each conjunct transcribes one claim of [pp] about its algorithm (cited inline); these
claims are internal and unrefereed and are NOT proved here. In particular the first conjunct carries
[pp] Prop. prop:BS (the search recognises `N_ℓ(A) ≥ 1`) and the SAMPLE_FIELD rejection logic; the
`splitCost` conjunct carries Prop. prop:select(ii) and with it Theorem thm:DD / Remark rem:T (each
round fruitful with probability ≥ 1/2); the `Tpost` conjunct carries Prop. prop:lv6 (with its
Lagarias–Odlyzko / Stark / Landau / Wiedemann inputs); the last conjunct carries the factor-tree and
race accounting of §sec:uniform. Note: these inequalities constrain real numbers; they do not tie the
model to a program, so the theorems below are conditional COST ASSEMBLY statements. -/
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
    (∀ ε : ℝ, 0 < ε → ∃ N₁ : ℝ, ∀ N : ℕ, N₁ ≤ N → Params.SplitInput a N →
      M.Tpost N ≤ M.Mrows N ^ 2 * (Real.log N) ^ C + (Params.Bsm a N : ℝ) ^ (1 + ε)) ∧
    -- [pp] §sec:uniform (proof of Theorem thm:main): at most `log₂N` calls of SPLIT, each raced
    -- against trial division (constant cost below `N₀`), on inputs `≤ N`, plus preprocessing
    -- (primality, perfect powers, trial division to `B`) at every node.
    (∀ N : ℕ, 2 ≤ N → M.factorCost N ≤ C + (Real.logb 2 N + 1) *
      (C * (Params.Bsm a N : ℝ) * (Real.log N) ^ C + 2 * supSplit M.splitCost a N₀ N))

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

/-- INTERNAL HYPOTHESIS (algorithmic interface for Theorem RA) at smoothness exponent `a`: as for
`CostInterfaceA`, plus [pp] Lemma lem:generala (A1), (A2), and the recursion inequality of the proof of
Theorem thm:RA, which carries Lemmas lem:OS, lem:waldrec, lem:Trec, the factor tree and the race against
trial division. Not proved here. -/
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
    (∀ ε : ℝ, 0 < ε → ∃ N₁ : ℝ, ∀ N : ℕ, N₁ ≤ N → Params.SplitInput a N →
      M.Tpost N ≤ M.Mrows N ^ 2 * (Real.log N) ^ C + (Params.Bsm a N : ℝ) ^ (1 + ε)) ∧
    -- `F` is finite, nonnegative, nondecreasing ([pp] Theorem thm:RA, "Finiteness").
    (∀ x, 2 ≤ x → 0 ≤ M.Fsup x) ∧ MonotoneOn M.Fsup (Set.Ici 2) ∧
    -- [pp] Theorem thm:RA, "Recursion": factor tree (< 2 log₂x splitting nodes), race against trial
    -- division (factor 2), Lemma OS, Lemma lem:waldrec, and at most two oracle calls per
    -- DESCEND_ONCE call on inputs `≤ √(2x)`.
    (∀ x : ℝ, N₀ ≤ x → M.Fsup x ≤ 2 * Real.logb 2 x *
      (C + C * (Params.Bsm a ⌊x⌋₊ : ℝ) * (Real.log x) ^ C + 2 * supSplit M.nonO a N₀ x +
        2 * (2 * supSplit M.tau a N₀ x) * M.Fsup (Real.sqrt (2 * x))))

end

section
open Finset

/-! ### (from LiftAveraging.lean: smooth numbers in progressions and characters) -/

/-- TRUSTED. `n` is `y`-smooth: every prime factor `p` of `n` satisfies `p ≤ y`
(used only for `n ≥ 1`). -/
def IsYSmooth (y : ℝ) (n : ℕ) : Prop := ∀ p : ℕ, p.Prime → p ∣ n → (p : ℝ) ≤ y

/-- TRUSTED. `S(z, y)`: the `y`-smooth integers `1 ≤ n ≤ z`. -/
def smoothSet (z y : ℝ) : Finset ℕ := (Finset.Icc 1 ⌊z⌋₊).filter (IsYSmooth y)

/-- TRUSTED. `Ψ(z, y; q, a) = #{n ∈ S(z, y) : n ≡ a (mod q)}`. -/
def PsiAP (z y : ℝ) (q a : ℕ) : ℕ := ((smoothSet z y).filter (fun n => n ≡ a [MOD q])).card

/-- TRUSTED. `Ψ_q(z, y) = #{n ∈ S(z, y) : gcd(n, q) = 1}`. -/
def PsiCop (z y : ℝ) (q : ℕ) : ℕ := ((smoothSet z y).filter (fun n => Nat.Coprime n q)).card

/-- TRUSTED. `E(z, y; q, a) = Ψ(z, y; q, a) − Ψ_q(z, y)/φ(q)`. -/
def Eerr (z y : ℝ) (q a : ℕ) : ℝ := (PsiAP z y q a : ℝ) - (PsiCop z y q : ℝ) / (q.totient : ℝ)

/-- TRUSTED. `Ψ(z, y; χ) = ∑_{n ∈ S(z, y)} χ(n)` for a Dirichlet character `χ` mod `q`. -/
def PsiChi (z y : ℝ) {q : ℕ} (χ : DirichletCharacter ℂ q) : ℂ :=
  ∑ n ∈ smoothSet z y, χ (n : ZMod q)

/-- TRUSTED. `Δ_ℓ(t) = φ(t)⁻¹ ∑_{χ mod t, ℓ ∣ cond χ} |Ψ(z, y; χ)|`. -/
def DeltaEll (z y : ℝ) (ℓ t : ℕ) : ℝ :=
  (1 / (t.totient : ℝ)) *
    ∑ χ ∈ (univ : Finset (DirichletCharacter ℂ t)).filter (fun χ => ℓ ∣ χ.conductor),
      ‖PsiChi z y χ‖

/-- TRUSTED. `Δ(q) = φ(q)⁻¹ ∑_{χ mod q, χ ≠ χ₀} |Ψ(z, y; χ)|`. -/
def DeltaAll (z y : ℝ) (q : ℕ) : ℝ :=
  (1 / (q.totient : ℝ)) *
    ∑ χ ∈ (univ : Finset (DirichletCharacter ℂ q)).filter (fun χ => χ ≠ 1), ‖PsiChi z y χ‖

end

end FactoringA
