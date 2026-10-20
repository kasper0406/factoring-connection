import LeanA.RelationSupply
import LeanA.Recursion
import LeanA.Exponents
import LeanA.AsymptoticsMain

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
fact:LP; the unnumbered Theorem of §1, p. 483, proved in §10): there is a probabilistic algorithm that
factors any positive integer `n` completely and takes expected time at most `L_n[1/2, 1 + o(1)]`.
Stated for `T n` = the maximum over inputs `≤ n` of that expected time. As a `Prop` about an arbitrary
`T` this is only meaningful together with the MODELLING assumption that a cost model's `TLP` is that
function. -/
def LenstraPomerance (T : ℕ → ℝ) : Prop := LBound (1 / 2) 1 T

/-- The maximum of `f` over the split inputs `n` with `N₀ ≤ n ≤ x` (`0` if there are none). -/
def supSplit (f : ℕ → ℝ) (a N₀ x : ℝ) : ℝ :=
  sSup (f '' {n : ℕ | N₀ ≤ n ∧ (n : ℝ) ≤ x ∧ Params.SplitInput a n})

/-! ## Helpers for the cost assembly -/

lemma mt_ev_of_ge {N₀ : ℝ} {R : ℕ → Prop} (h : ∀ N : ℕ, N₀ ≤ N → R N) :
    ∀ᶠ N : ℕ in atTop, R N :=
  (tendsto_natCast_atTop_atTop.eventually_ge_atTop N₀).mono h

lemma mt_rpow_neg_half : (2:ℝ) ^ (-(1 / 2 : ℝ)) = 1 / Real.sqrt 2 := by
  rw [Real.rpow_neg (by norm_num), ← Real.sqrt_eq_rpow, one_div]

lemma mt_adm (P : Params) (hP : P.Admissible) :
    0 < (P.θ : ℝ) ∧ (P.θ : ℝ) < 1 ∧ 0 < (P.h : ℝ) ∧ (P.h : ℝ) ≤ 1 ∧
      0 < (((P.μ + P.ν + P.η₀ / 2 : ℚ)) : ℝ) ∧ (3 : ℝ) ≤ P.d := by
  obtain ⟨_, hν, hθ⟩ := admissible_facts P hP
  obtain ⟨_, hd, _, hη, hμ, _, hθ1⟩ := hP
  have hd' : (3 : ℚ) ≤ P.d := by exact_mod_cast hd
  have hθ0 : (0 : ℚ) < P.θ := lt_trans (by positivity) hθ
  have hh0 : (0 : ℚ) < P.h := by
    have := (hStar_ge P).1
    unfold Params.h; linarith
  have hh1 : P.h ≤ 1 := by
    have e : P.θ = 2 * P.μ + P.h + ((P.d : ℚ) - 1) * P.ν := rfl
    have : 0 ≤ ((P.d : ℚ) - 1) * P.ν := mul_nonneg (by linarith) hν.le
    linarith
  refine ⟨by exact_mod_cast hθ0, by exact_mod_cast hθ1, by exact_mod_cast hh0,
    by exact_mod_cast hh1, ?_, by exact_mod_cast hd⟩
  have : (0 : ℚ) < P.μ + P.ν + P.η₀ / 2 := by linarith
  exact_mod_cast this

lemma mt_Lv_rpow (N : ℕ) (b : ℝ) : (Params.Lv N) ^ b = LN (1 / 2) b N := by
  unfold Params.Lv LN
  rw [← Real.exp_mul]; congr 1; ring

lemma mt_Bsm_le (a : ℝ) (N : ℕ) : (Params.Bsm a N : ℝ) ≤ LN (1 / 2) a N := by
  unfold Params.Bsm
  rw [← mt_Lv_rpow]
  exact Nat.floor_le (Real.rpow_nonneg (Real.exp_pos _).le _)

lemma mt_Ssz_sq_le (P : Params) (a : ℝ) (ha : 0 < a) (hθ : 0 ≤ (P.θ : ℝ)) (N : ℕ)
    (hN : (1 : ℝ) ≤ N) :
    (P.Ssz a N : ℝ) ^ 2 ≤ 4 * LN (1 / 2) ((P.θ : ℝ) / (2 * a)) N := by
  have hq : 0 ≤ (P.θ : ℝ) / (4 * a) := by positivity
  have h1 : 1 ≤ LN (1 / 2) ((P.θ : ℝ) / (4 * a)) N := one_le_LN _ _ hN hq
  have h2 : (P.Ssz a N : ℝ) ≤ 2 * LN (1 / 2) ((P.θ : ℝ) / (4 * a)) N := by
    unfold Params.Ssz
    rw [mt_Lv_rpow]
    have := Nat.ceil_lt_add_one (le_trans zero_le_one h1)
    linarith
  have h0 : (0 : ℝ) ≤ P.Ssz a N := Nat.cast_nonneg _
  have e : LN (1 / 2) ((P.θ : ℝ) / (2 * a)) N =
      LN (1 / 2) ((P.θ : ℝ) / (4 * a)) N * LN (1 / 2) ((P.θ : ℝ) / (4 * a)) N := by
    rw [LN_mul]; congr 1; field_simp; ring
  rw [e]
  nlinarith

lemma mt_primeCounting_le (n : ℕ) : (Nat.primeCounting n : ℝ) ≤ n + 1 := by
  have : Nat.primeCounting n ≤ n + 1 := by
    unfold Nat.primeCounting Nat.primeCounting'
    exact Nat.count_le _
  exact_mod_cast this

/-- The explicit Lemma-reject bound is `≤ L^{−1/2}/2` eventually. -/
lemma mt_reject_small (P : Params) (hP : P.Admissible) :
    ∀ᶠ N : ℕ in atTop,
      Real.exp (-(Real.log (P.Hb N) / (13 * P.d * Real.log (Real.log (P.Hb N))))) +
          2 / (2 * (P.Hb N : ℝ) + 1) + 1 / (P.Wb N : ℝ) ≤ LN (1 / 2) (-(1 / 2)) N / 2 := by
  obtain ⟨-, -, hh0, hh1, hw0, hd⟩ := mt_adm P hP
  have hD0 : (0 : ℝ) < 13 * (P.d : ℝ) := by linarith
  set h : ℝ := (P.h : ℝ) with hh
  set w : ℝ := ((P.μ + P.ν + P.η₀ / 2 : ℚ) : ℝ) with hw
  set D : ℝ := 13 * (P.d : ℝ) with hD
  set κ : ℝ := h ^ 2 / (4 * D ^ 2) with hκ
  have hκ0 : 0 < κ := by positivity
  have hlogt := mt_tendsto_log_nat
  filter_upwards [hlogt.eventually (mt_ev_log_le_mul κ hκ0 3),
    hlogt.eventually (Real.tendsto_log_atTop.eventually_ge_atTop 1),
    hlogt.eventually_ge_atTop (2 * Real.log 2 / h),
    hlogt.eventually_ge_atTop (4 / h),
    mt_ev_Lam_le (h / 2) (by positivity), mt_ev_Lam_le w hw0,
    mt_ev_const_le (1 / 2) (by norm_num) (by norm_num) (1 / 2) (by norm_num) 6,
    eventually_ge_atTop 1] with N E1 E2 E3 E4 E6 E7 E8 hN1
  have hN1' : (1:ℝ) ≤ N := by exact_mod_cast hN1
  have hN0 : (0:ℝ) < N := by linarith
  unfold LN at E8 ⊢
  have hΛ0 : 0 ≤ Lam (1 / 2) (N : ℝ) := Lam_nonneg _ _ hN1'
  have ht0 : 0 ≤ Real.log (N : ℝ) := Real.log_nonneg hN1'
  have hΛsq : Lam (1 / 2) (N : ℝ) ^ 2 = Real.log N * Real.log (Real.log N) := by
    rw [mt_Lam_half_sq ht0, max_eq_right E2]
  -- the sizes of `H` and `W`
  have hNh : (N:ℝ) ^ h = Real.exp (h * Real.log N) := by
    rw [Real.rpow_def_of_pos hN0]; ring_nf
  have hNw : (N:ℝ) ^ w = Real.exp (w * Real.log N) := by
    rw [Real.rpow_def_of_pos hN0]; ring_nf
  have hNh1 : 1 ≤ (N:ℝ) ^ h := Real.one_le_rpow hN1' hh0.le
  have hHdef : (P.Hb N : ℝ) = (⌊(N:ℝ) ^ h⌋₊ : ℝ) := rfl
  have hWdef : (P.Wb N : ℝ) = (⌈(N:ℝ) ^ w⌉₊ : ℝ) := rfl
  have hHlo : (N:ℝ) ^ h / 2 ≤ (P.Hb N : ℝ) := by
    rw [hHdef]
    have h1 := Nat.lt_floor_add_one ((N:ℝ) ^ h)
    have h2 : (1:ℝ) ≤ (⌊(N:ℝ) ^ h⌋₊ : ℝ) := by exact_mod_cast (Nat.one_le_floor_iff _).2 hNh1
    linarith
  have hHhi : (P.Hb N : ℝ) ≤ N := by
    rw [hHdef]
    refine (Nat.floor_le (by linarith)).trans ?_
    calc (N:ℝ) ^ h ≤ (N:ℝ) ^ (1:ℝ) := Real.rpow_le_rpow_of_exponent_le hN1' hh1
      _ = N := Real.rpow_one _
  have hWlo : Real.exp (w * Real.log N) ≤ (P.Wb N : ℝ) := by
    rw [hWdef, ← hNw]; exact Nat.le_ceil _
  generalize (P.Hb N : ℝ) = H at hHlo hHhi ⊢
  generalize (P.Wb N : ℝ) = W at hWlo ⊢
  set t := Real.log (N : ℝ) with ht
  set Λ := Lam (1 / 2) (N : ℝ) with hΛ
  have hH0 : 0 < H := by linarith
  have hlogH : h * t / 2 ≤ Real.log H := by
    have e : Real.log ((N:ℝ) ^ h / 2) = h * t - Real.log 2 := by
      rw [Real.log_div (by positivity) (by norm_num), Real.log_rpow hN0]
    have := Real.log_le_log (by positivity) hHlo
    rw [div_le_iff₀ hh0] at E3
    linarith
  have hlogH1 : 1 < Real.log H := by
    rw [div_le_iff₀ hh0] at E4
    linarith
  have hlogHt : Real.log H ≤ t := Real.log_le_log hH0 hHhi
  have hll0 : 0 < Real.log (Real.log H) := Real.log_pos hlogH1
  have hllt : Real.log (Real.log H) ≤ Real.log t := Real.log_le_log (by linarith) hlogHt
  have hlt0 : 0 ≤ Real.log t := by linarith
  -- term 1
  have hk : (Λ * D * Real.log t) ^ 2 ≤ (h / 2 * t) ^ 2 := by
    have hDt : 0 ≤ D ^ 2 * t := by positivity
    calc (Λ * D * Real.log t) ^ 2 = Λ ^ 2 * D ^ 2 * (Real.log t) ^ 2 := by ring
      _ = D ^ 2 * t * (Real.log t) ^ 3 := by rw [hΛsq]; ring
      _ ≤ D ^ 2 * t * (κ * t) := mul_le_mul_of_nonneg_left E1 hDt
      _ = (h / 2 * t) ^ 2 := by rw [hκ]; field_simp; ring
  have hk' : Λ * D * Real.log t ≤ h / 2 * t :=
    (pow_le_pow_iff_left₀ (mul_nonneg (mul_nonneg hΛ0 hD0.le) hlt0) (by positivity)
      two_ne_zero).1 hk
  have hT1 : Real.exp (-(Real.log H / (D * Real.log (Real.log H)))) ≤ Real.exp (-Λ) := by
    apply Real.exp_le_exp.2
    have : Λ ≤ Real.log H / (D * Real.log (Real.log H)) := by
      rw [le_div_iff₀ (by positivity)]
      have : Λ * (D * Real.log (Real.log H)) ≤ Λ * (D * Real.log t) :=
        mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_left hllt hD0.le) hΛ0
      nlinarith
    linarith
  -- term 2
  have hT2 : 2 / (2 * H + 1) ≤ Real.exp (-Λ) := by
    have h1 : 2 / (2 * H + 1) ≤ 2 / Real.exp (h * t) := by
      rw [← hNh]
      exact div_le_div_of_nonneg_left (by norm_num) (by positivity) (by linarith)
    refine h1.trans ?_
    rw [div_le_iff₀ (Real.exp_pos _), ← Real.exp_add]
    rw [div_le_iff₀ hh0] at E3
    calc (2:ℝ) = Real.exp (Real.log 2) := (Real.exp_log (by norm_num)).symm
      _ ≤ _ := Real.exp_le_exp.2 (by linarith)
  -- term 3
  have hT3 : 1 / W ≤ Real.exp (-Λ) := by
    have hW0 : 0 < W := lt_of_lt_of_le (Real.exp_pos _) hWlo
    rw [div_le_iff₀ hW0]
    calc (1:ℝ) ≤ Real.exp (-Λ + w * t) := Real.one_le_exp (by linarith)
      _ = Real.exp (-Λ) * Real.exp (w * t) := Real.exp_add _ _
      _ ≤ Real.exp (-Λ) * W := mul_le_mul_of_nonneg_left hWlo (Real.exp_pos _).le
  have hfin : 6 * Real.exp (-Λ) ≤ Real.exp (-(1 / 2) * Λ) := by
    have e : Real.exp (-(1 / 2) * Λ) = Real.exp (1 / 2 * Λ) * Real.exp (-Λ) := by
      rw [← Real.exp_add]; ring_nf
    rw [e]
    exact mul_le_mul_of_nonneg_right E8 (Real.exp_pos _).le
  linarith

/-- The success probability of one trial is `≥ L^{−o(1)}` on split inputs: Lemma lem:reject makes
`Pr(𝓡) ≤ γ/2` (with `γ ≥ L^{−1/2}`, `relation_supply`), and `descent_supply` applies. -/
lemma mt_pTrial_lower (P : Params) (hP : P.Admissible) (a : ℝ) (ha : 0 < a) (hCEP : CEP)
    (hRS : RosserSchoenfeld) (hFM : FirstMoment P a) (hLD : SecondMoment P a)
    (pT : ℕ → ℝ) (rej : (N : ℕ) → Finset P.Field)
    (h : ∀ᶠ N : ℕ in atTop, Params.SplitInput a N →
      expect (P.fields N) (fun A => if A ∈ rej N then 0 else P.descentProb a N A) ≤ pT N ∧
      rej N ⊆ P.fields N ∧
      prob (P.fields N) (fun A => A ∈ rej N) ≤
        Real.exp (-(Real.log (P.Hb N) / (13 * P.d * Real.log (Real.log (P.Hb N))))) +
          2 / (2 * (P.Hb N : ℝ) + 1) + 1 / (P.Wb N : ℝ)) :
    ∀ ε > 0, ∀ᶠ N : ℕ in atTop, Params.SplitInput a N → LN (1 / 2) (-ε) N ≤ pT N := by
  intro ε hε
  obtain ⟨N₂, hN₂⟩ := relation_supply P hP a ha hCEP hFM hLD (1 / 2) (by norm_num)
  obtain ⟨N₃, hN₃⟩ := descent_supply P hP a ha hCEP hRS hFM hLD ε hε
  filter_upwards [h, mt_reject_small P hP, eventually_ge_atTop N₂, eventually_ge_atTop N₃]
    with N h1 h2 h3 h4 hQ
  obtain ⟨he, hsub, hpr⟩ := h1 hQ
  have hγ := hN₂ N h3 hQ
  exact (hN₃ N h4 hQ (rej N) hsub (by linarith)).trans he

lemma mt_supSplit_nonneg (f : ℕ → ℝ) (a N₀ x : ℝ)
    (h : ∀ n : ℕ, N₀ ≤ n → Params.SplitInput a n → 0 ≤ f n) : 0 ≤ supSplit f a N₀ x := by
  unfold supSplit
  apply Real.sSup_nonneg
  rintro _ ⟨n, ⟨hn0, -, hs⟩, rfl⟩
  exact h n hn0 hs

/-- `supSplit f ≤ K + L_x^c` once `f n ≤ L_n^c` for all large split inputs. -/
lemma mt_supSplit_le (f : ℕ → ℝ) (a N₀ c : ℝ) (hc : 0 ≤ c)
    (h : ∀ᶠ N : ℕ in atTop, Params.SplitInput a N → f N ≤ LN (1 / 2) c N) :
    ∃ K : ℝ, 0 ≤ K ∧ ∀ x : ℝ, supSplit f a N₀ x ≤ K + LN (1 / 2) c x := by
  obtain ⟨N₁, hN₁⟩ := eventually_atTop.1 h
  have hK0 : 0 ≤ ∑ n ∈ Finset.range N₁, |f n| := Finset.sum_nonneg (fun _ _ => abs_nonneg _)
  refine ⟨∑ n ∈ Finset.range N₁, |f n|, hK0, fun x => ?_⟩
  have hLx := mt_LN_pos (1 / 2) c x
  unfold supSplit
  apply Real.sSup_le _ (by linarith)
  rintro _ ⟨n, ⟨-, hnx, hs⟩, rfl⟩
  have hn1 : (1:ℝ) ≤ n := by exact_mod_cast hs.1.le
  by_cases hn : n < N₁
  · have h1 : f n ≤ |f n| := le_abs_self _
    have h2 : |f n| ≤ ∑ n ∈ Finset.range N₁, |f n| :=
      Finset.single_le_sum (f := fun n => |f n|) (fun _ _ => abs_nonneg _)
        (Finset.mem_range.2 hn)
    linarith
  · push Not at hn
    have h1 := hN₁ n hn hs
    have h2 := mt_LN_mono_x (1 / 2) c (by norm_num) (by norm_num) hc hn1 hnx
    linarith

/-! ### The component bounds on split inputs -/

lemma mt_chain_Mrows (P : Params) (a : ℝ) (ha : 0 ≤ a) (Mr : ℕ → ℝ) (C : ℝ)
    (h : ∀ᶠ N : ℕ in atTop, Params.SplitInput a N →
      Mr N ≤ (P.d + 1) * (Nat.primeCounting (Params.Bsm a N) : ℝ) + C * Real.log N) :
    EvLB (Params.SplitInput a) a Mr := by
  have h1 : EvLB (Params.SplitInput a) (max (max a 0) 0)
      (fun N => ((P.d : ℝ) + 1) * (LN (1 / 2) a N + 1) + C * (Real.log N) ^ (1:ℝ)) :=
    (EvLB.const_mul (by positivity) ((EvLB.of_LN a).add (EvLB.const 1))).add (EvLB.polylog C 1)
  refine (h1.of_le ?_).mono (max_le (max_le le_rfl ha) ha)
  filter_upwards [h] with N hN hQ
  have h2 := hN hQ
  have hπ := mt_primeCounting_le (Params.Bsm a N)
  have hB := mt_Bsm_le a N
  have hd : (0:ℝ) ≤ P.d + 1 := by positivity
  rw [Real.rpow_one]
  have : ((P.d : ℝ) + 1) * (Nat.primeCounting (Params.Bsm a N) : ℝ) ≤
      ((P.d : ℝ) + 1) * (LN (1 / 2) a N + 1) :=
    mul_le_mul_of_nonneg_left (by linarith) hd
  linarith

lemma mt_chain_Tsq (P : Params) (a : ℝ) (ha : 0 < a) (hθ : 0 ≤ (P.θ : ℝ)) (Ts : ℕ → ℝ) (C : ℝ)
    (h : ∀ᶠ N : ℕ in atTop, Params.SplitInput a N →
      Ts N ≤ (Real.log N) ^ C * ((Params.Bsm a N : ℝ) + (P.Ssz a N : ℝ) ^ 2 * Real.log N)) :
    EvLB (Params.SplitInput a) (max a ((P.θ : ℝ) / (2 * a))) Ts := by
  set q := (P.θ : ℝ) / (2 * a) with hq
  have h1 : EvLB (Params.SplitInput a) (0 + max a (q + 0))
      (fun N => 1 * (Real.log N) ^ C *
        (LN (1 / 2) a N + 4 * LN (1 / 2) q N * (1 * (Real.log N) ^ (1:ℝ)))) := by
    refine (EvLB.polylog 1 C).mul ((EvLB.of_LN a).add
      ((EvLB.const_mul (by norm_num) (EvLB.of_LN q)).mul (EvLB.polylog 1 1)
        (Eventually.of_forall fun N _ => mul_nonneg (by norm_num) (mt_LN_pos _ _ _).le))) ?_
    filter_upwards [eventually_ge_atTop 1] with N hN _
    have : (1:ℝ) ≤ N := by exact_mod_cast hN
    have := Real.rpow_nonneg (Real.log_nonneg this) C
    linarith
  refine (h1.of_le ?_).mono (by simp)
  filter_upwards [h, eventually_ge_atTop 1] with N hN hN1 hQ
  have hN1' : (1:ℝ) ≤ N := by exact_mod_cast hN1
  have hl0 : 0 ≤ Real.log N := Real.log_nonneg hN1'
  have hlC : 0 ≤ (Real.log N) ^ C := Real.rpow_nonneg hl0 C
  have hS := mt_Ssz_sq_le P a ha hθ N hN1'
  have hB := mt_Bsm_le a N
  rw [Real.rpow_one, one_mul, one_mul]
  refine (hN hQ).trans (mul_le_mul_of_nonneg_left (add_le_add hB ?_) hlC)
  exact mul_le_mul_of_nonneg_right hS hl0

lemma mt_chain_Tpost (a : ℝ) (ha : 0 ≤ a) (Mr Tp : ℕ → ℝ) (C : ℝ)
    (hMr : EvLB (Params.SplitInput a) a Mr)
    (hMr1 : ∀ᶠ N : ℕ in atTop, Params.SplitInput a N → 1 ≤ Mr N)
    (h : ∀ ε : ℝ, 0 < ε → ∃ N₁ : ℝ, ∀ N : ℕ, N₁ ≤ N → Params.SplitInput a N →
      Tp N ≤ Mr N ^ 2 * (Real.log N) ^ C + (Params.Bsm a N : ℝ) ^ (1 + ε)) :
    EvLB (Params.SplitInput a) (2 * a) Tp := by
  obtain ⟨N₁, hN₁⟩ := h 1 one_pos
  have h1 : EvLB (Params.SplitInput a) (max (a + a + 0) (2 * a))
      (fun N => Mr N * Mr N * (1 * (Real.log N) ^ C) + LN (1 / 2) (2 * a) N) := by
    refine ((hMr.mul hMr ?_).mul (EvLB.polylog 1 C) ?_).add (EvLB.of_LN (2 * a))
    · filter_upwards [hMr1] with N hN hQ
      linarith [hN hQ]
    · filter_upwards [hMr1] with N hN hQ
      have := hN hQ
      nlinarith
  refine (h1.of_le ?_).mono (max_le (by linarith) le_rfl)
  filter_upwards [mt_ev_of_ge hN₁] with N hN hQ
  have h2 := hN hQ
  have hB0 : (0:ℝ) ≤ Params.Bsm a N := Nat.cast_nonneg _
  have hB : (Params.Bsm a N : ℝ) ^ ((1:ℝ) + 1) ≤ LN (1 / 2) (2 * a) N := by
    rw [show (1:ℝ) + 1 = ((2:ℕ):ℝ) by norm_num, Real.rpow_natCast]
    calc (Params.Bsm a N : ℝ) ^ 2 ≤ (LN (1 / 2) a N) ^ 2 :=
          pow_le_pow_left₀ hB0 (mt_Bsm_le a N) _
      _ = LN (1 / 2) a N * LN (1 / 2) a N := by ring
      _ = LN (1 / 2) (2 * a) N := by rw [LN_mul]; congr 1; ring
  rw [one_mul, ← sq]
  exact h2.trans (add_le_add_right hB _)

/-- `TLP(⌊√(2N)⌋) ≤ L_N[1/2, 1/√2 + o(1)]`. -/
lemma mt_TLP_sqrt (Q : ℕ → Prop) (T : ℕ → ℝ) (hT : LBound (1 / 2) 1 T) :
    EvLB Q (1 / Real.sqrt 2) (fun N => T ⌊Real.sqrt (2 * N)⌋₊) := by
  intro ε hε
  set δ := min (ε / 4) 1 with hδ
  have hδ0 : 0 < δ := lt_min (by linarith) one_pos
  have hδ1 : δ ≤ 1 := min_le_right _ _
  have hδε : δ ≤ ε / 4 := min_le_left _ _
  have hm : Tendsto (fun N : ℕ => ⌊Real.sqrt (2 * (N : ℝ))⌋₊) atTop atTop :=
    tendsto_nat_floor_atTop.comp (Real.tendsto_sqrt_atTop.comp
      (tendsto_natCast_atTop_atTop.const_mul_atTop two_pos))
  have hr0 : 0 < 1 / Real.sqrt 2 := by positivity
  have hr1 : 1 / Real.sqrt 2 ≤ 1 := by
    rw [div_le_one (by positivity)]
    exact Real.one_le_sqrt.2 (by norm_num)
  filter_upwards [hm.eventually (hT δ hδ0), hm.eventually (eventually_ge_atTop 1),
    tendsto_natCast_atTop_atTop.eventually
      (Lam_sqrt_le (1 / 2) (by norm_num) (by norm_num) δ hδ0),
    eventually_ge_atTop 1] with N h1 h2 h3 hN1 _
  have hN1' : (1:ℝ) ≤ N := by exact_mod_cast hN1
  have hm1 : (1:ℝ) ≤ (⌊Real.sqrt (2 * (N : ℝ))⌋₊ : ℝ) := by exact_mod_cast h2
  have hmle : (⌊Real.sqrt (2 * (N : ℝ))⌋₊ : ℝ) ≤ Real.sqrt (2 * N) :=
    Nat.floor_le (Real.sqrt_nonneg _)
  have h4 := mt_LN_mono_x (1 / 2) (1 + δ) (by norm_num) (by norm_num) (by linarith) hm1 hmle
  refine h1.trans (h4.trans ?_)
  rw [mt_rpow_neg_half] at h3
  have hΛ0 := Lam_nonneg (1 / 2) (N : ℝ) hN1'
  unfold LN
  apply Real.exp_le_exp.2
  calc (1 + δ) * Lam (1 / 2) (Real.sqrt (2 * N))
      ≤ (1 + δ) * ((1 / Real.sqrt 2 + δ) * Lam (1 / 2) (N : ℝ)) :=
        mul_le_mul_of_nonneg_left h3 (by linarith)
    _ = ((1 + δ) * (1 / Real.sqrt 2 + δ)) * Lam (1 / 2) (N : ℝ) := by ring
    _ ≤ (1 / Real.sqrt 2 + ε) * Lam (1 / 2) (N : ℝ) := by
        apply mul_le_mul_of_nonneg_right _ hΛ0
        nlinarith

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

/-- Unpacking the interface (the only place that depends on its field layout). -/
lemma mt_unpackA {P : Params} {a : ℝ} {M : CostModelA P} (hM : CostInterfaceA P a M) :
    ∃ C N₀ : ℝ, 0 ≤ C ∧
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
        (C * (Params.Bsm a N : ℝ) * (Real.log N) ^ C + 2 * supSplit M.splitCost a N₀ N)) :=
  hM.bounds

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
  obtain ⟨hθ0, hθ1, -, -, -, -⟩ := mt_adm P hP
  obtain ⟨hss, hs0, hs1, -⟩ := sqrt2_facts
  have hq : (P.θ : ℝ) / (2 * ((P.θ : ℝ) / Real.sqrt 2)) = 1 / Real.sqrt 2 := half_div _ hθ0
  set a := (P.θ : ℝ) / Real.sqrt 2 with ha_def
  have ha : 0 < a := by positivity
  have ha1 : a ≤ 1 / Real.sqrt 2 := div_le_div_of_nonneg_right hθ1.le hs0.le
  have hcA : P.cA = a + 1 / Real.sqrt 2 := by unfold Params.cA; rw [ha_def]; ring
  have h2a : 2 * a ≤ P.cA := by rw [hcA]; linarith
  obtain ⟨C, N₀, hC, hb, hTp, hfac⟩ := mt_unpackA hM
  have hb' := mt_ev_of_ge hb
  have hp : ∀ ε > 0, ∀ᶠ N : ℕ in atTop, Params.SplitInput a N → LN (1 / 2) (-ε) N ≤ M.pTrial N :=
    mt_pTrial_lower P hP a ha hCEP hRS hFM hLD M.pTrial M.rej
      (hb'.mono fun N h hQ => ⟨(h hQ).1, (h hQ).2.1, (h hQ).2.2.1⟩)
  have hMr : EvLB (Params.SplitInput a) a M.Mrows :=
    mt_chain_Mrows P a ha.le M.Mrows C (hb'.mono fun N h hQ => (h hQ).2.2.2.2.2.1)
  have hMr1 : ∀ᶠ N : ℕ in atTop, Params.SplitInput a N → 1 ≤ M.Mrows N :=
    hb'.mono fun N h hQ => (h hQ).2.2.2.2.1
  have hTsq : EvLB (Params.SplitInput a) (max a (1 / Real.sqrt 2)) M.Tsq := by
    have := mt_chain_Tsq P a ha hθ0.le M.Tsq C
      (hb'.mono fun N h hQ => (h hQ).2.2.2.2.2.2.2.2.2.1)
    rwa [hq] at this
  have hTpo : EvLB (Params.SplitInput a) (2 * a) M.Tpost :=
    mt_chain_Tpost a ha.le M.Mrows M.Tpost C hMr hMr1 hTp
  have hTsh : EvLB (Params.SplitInput a) (1 / Real.sqrt 2) M.Tshape := by
    have h1 : EvLB (Params.SplitInput a) (max (1 / Real.sqrt 2) 0)
        (fun N => 2 * M.TLP ⌊Real.sqrt (2 * N)⌋₊ + 1 * (Real.log N) ^ C) :=
      (EvLB.const_mul (by norm_num) (mt_TLP_sqrt _ M.TLP hLP)).add (EvLB.polylog 1 C)
    refine (h1.of_le ?_).mono (max_le le_rfl (by positivity))
    filter_upwards [hb'] with N h hQ
    rw [one_mul]; exact (h hQ).2.2.2.2.2.2.2.2.1
  have hCt : EvLB (Params.SplitInput a) (1 / Real.sqrt 2) M.Ctrial := by
    have h1 : EvLB (Params.SplitInput a)
        (max (max (1 / Real.sqrt 2) (max a (1 / Real.sqrt 2))) 0)
        (fun N => M.Tshape N + 2 * P.k * M.Tsq N + 1 * (Real.log N) ^ C) :=
      (hTsh.add (EvLB.const_mul (by positivity) hTsq)).add (EvLB.polylog 1 C)
    refine (h1.of_le ?_).mono (max_le (max_le le_rfl (max_le ha1 le_rfl)) (by positivity))
    filter_upwards [hb'] with N h hQ
    rw [one_mul]; exact (h hQ).2.2.2.2.2.2.2.1
  have hsplit : EvLB (Params.SplitInput a) P.cA M.splitCost := by
    have h1 := (EvLB.const_mul (by norm_num : (0:ℝ) ≤ 3)
      ((hMr.mul hCt (hMr1.mono fun N h hQ => by linarith [h hQ])).add hTpo)).div hp
    refine (h1.of_le ?_).mono (max_le (le_of_eq hcA.symm) h2a)
    filter_upwards [hb'] with N h hQ
    exact (h hQ).2.2.2.1
  -- the factor tree
  intro ε hε
  have hc0 : 0 ≤ P.cA + ε / 2 := by rw [hcA]; positivity
  obtain ⟨K, hK0, hK⟩ := mt_supSplit_le M.splitCost a N₀ (P.cA + ε / 2) hc0
    (hsplit (ε / 2) (by linarith))
  have hlog : EvLB (fun _ => True) 0 (fun N : ℕ => Real.logb 2 N + 1) := by
    refine (((EvLB.polylog (1 / Real.log 2) 1).add (EvLB.const 1)).of_le ?_).mono (by simp)
    filter_upwards with N _
    rw [Real.rpow_one, Real.logb]
    exact le_of_eq (by ring)
  have hlog0 : ∀ᶠ N : ℕ in atTop, True → 0 ≤ Real.logb 2 N + 1 := by
    filter_upwards [eventually_ge_atTop 1] with N hN _
    have := Real.logb_nonneg (b := 2) (by norm_num) (by exact_mod_cast hN : (1:ℝ) ≤ N)
    linarith
  have hB : EvLB (fun _ => True) a (fun N : ℕ => (Params.Bsm a N : ℝ)) :=
    (EvLB.of_LN a).of_le (Eventually.of_forall fun N _ => mt_Bsm_le a N)
  have hg0 := (EvLB.const C).add (hlog.mul (((EvLB.const_mul hC hB).mul (EvLB.polylog 1 C)
    (Eventually.of_forall fun N _ => mul_nonneg hC (Nat.cast_nonneg _))).add
    (EvLB.const_mul (by norm_num : (0:ℝ) ≤ 2) ((EvLB.const K).add (EvLB.of_LN (P.cA + ε / 2)))))
    hlog0)
  have hg : EvLB (fun _ => True) (P.cA + ε / 2) (fun N : ℕ =>
      C + (Real.logb 2 N + 1) * (C * (Params.Bsm a N : ℝ) * (Real.log N) ^ C +
        2 * (K + LN (1 / 2) (P.cA + ε / 2) N))) := by
    refine (hg0.of_le ?_).mono ?_
    · filter_upwards with N _
      exact le_of_eq (by ring)
    · have : a ≤ P.cA + ε / 2 := by rw [hcA]; linarith [hs0, one_div_pos.2 hs0]
      refine max_le hc0 ?_
      rw [zero_add, add_zero]
      exact max_le this (max_le hc0 le_rfl)
  filter_upwards [hg (ε / 2) (by linarith), eventually_ge_atTop 2] with N h1 h2
  have h2' : (2:ℝ) ≤ N := by exact_mod_cast h2
  have hlb : 0 ≤ Real.logb 2 N + 1 := by
    have := Real.logb_nonneg (b := 2) (by norm_num) (by linarith : (1:ℝ) ≤ N)
    linarith
  calc M.factorCost N ≤ C + (Real.logb 2 N + 1) *
        (C * (Params.Bsm a N : ℝ) * (Real.log N) ^ C + 2 * supSplit M.splitCost a N₀ N) :=
        hfac N h2
    _ ≤ C + (Real.logb 2 N + 1) * (C * (Params.Bsm a N : ℝ) * (Real.log N) ^ C +
        2 * (K + LN (1 / 2) (P.cA + ε / 2) N)) :=
        by
          have h3 := hK N
          have h4 := mul_le_mul_of_nonneg_left
            (show C * (Params.Bsm a N : ℝ) * (Real.log N) ^ C + 2 * supSplit M.splitCost a N₀ N ≤
              C * (Params.Bsm a N : ℝ) * (Real.log N) ^ C +
                2 * (K + LN (1 / 2) (P.cA + ε / 2) N) by linarith) hlb
          linarith
    _ ≤ LN (1 / 2) (P.cA + ε / 2 + ε / 2) N := h1 trivial
    _ = LN (1 / 2) (P.cA + ε) N := by congr 1; ring

/-- MAIN. [pp] Corollary cor:main, conditional form: for every `c > 1/√2` there is an admissible
tuple, with constant below `c`, whose FACTOR runs in expected time `O(L_N[1/2, c])` whenever that
tuple's hypotheses hold. -/
theorem corollaryA (hCEP : CEP) (hRS : RosserSchoenfeld) (c : ℝ) (hc : 1 / Real.sqrt 2 < c) :
    ∃ P : Params, P.Admissible ∧ P.cA < c ∧
      ∀ (_hFM : FirstMoment P ((P.θ : ℝ) / Real.sqrt 2))
        (_hLD : SecondMoment P ((P.θ : ℝ) / Real.sqrt 2))
        (M : CostModelA P) (_hM : CostInterfaceA P ((P.θ : ℝ) / Real.sqrt 2) M)
        (_hLP : LenstraPomerance M.TLP),
        ∃ K : ℝ, ∀ N : ℕ, 2 ≤ N → M.factorCost N ≤ K * LN (1 / 2) c N := by
  obtain ⟨P, hP, hlt⟩ := exists_admissible_cA_lt c hc
  refine ⟨P, hP, hlt, fun hFM hLD M hM hLP => ?_⟩
  have hc0 : 0 < c := lt_trans (by positivity) hc
  have hA := theoremA P hP hCEP hRS hFM hLD M hM hLP (c - P.cA) (by linarith)
  obtain ⟨N₁, hN₁⟩ := eventually_atTop.1 hA
  set K₀ := ∑ n ∈ Finset.range N₁, |M.factorCost n| with hK₀
  have hK₀0 : 0 ≤ K₀ := Finset.sum_nonneg (fun _ _ => abs_nonneg _)
  refine ⟨max 1 K₀, fun N hN => ?_⟩
  have hN1 : (1:ℝ) ≤ N := by exact_mod_cast (by omega : 1 ≤ N)
  have hL1 : 1 ≤ LN (1 / 2) c N := one_le_LN _ _ hN1 hc0.le
  by_cases hNN : N < N₁
  · have h1 : M.factorCost N ≤ |M.factorCost N| := le_abs_self _
    have h2 : |M.factorCost N| ≤ K₀ :=
      Finset.single_le_sum (f := fun n => |M.factorCost n|) (fun _ _ => abs_nonneg _)
        (Finset.mem_range.2 hNN)
    have h3 : K₀ ≤ max 1 K₀ := le_max_right _ _
    have h4 : max 1 K₀ ≤ max 1 K₀ * LN (1 / 2) c N :=
      le_mul_of_one_le_right (le_trans zero_le_one (le_max_left _ _)) hL1
    linarith
  · push Not at hNN
    have h1 := hN₁ N hNN
    rw [show P.cA + (c - P.cA) = c by ring] at h1
    have h4 : LN (1 / 2) c N ≤ max 1 K₀ * LN (1 / 2) c N :=
      le_mul_of_one_le_left (by linarith) (le_max_left _ _)
    linarith

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

/-- Unpacking the interface (the only place that depends on its field layout). -/
lemma mt_unpackRA {P : Params} {a : ℝ} {M : CostModelRA P} (hM : CostInterfaceRA P a M) :
    ∃ C N₀ : ℝ, 0 ≤ C ∧
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
          2 * (2 * supSplit M.tau a N₀ x) * M.Fsup (Real.sqrt (2 * x)))) :=
  hM.bounds

/-- MAIN. [pp] Theorem thm:RA (Theorem RA), conditional form. For admissible `(d, k, μ, η₀)` and any
fixed `a > 0`, the expected time of FACTOR_rec[a] is `≤ L_N[1/2, c(a) + o(1)]`,
`c(a) = max(a + θ/(2a), (2+√2)a)`. -/
theorem theoremRA (P : Params) (hP : P.Admissible) (a : ℝ) (ha : 0 < a) (hCEP : CEP)
    (hRS : RosserSchoenfeld) (hFM : FirstMoment P a) (hLD : SecondMoment P a)
    (M : CostModelRA P) (hM : CostInterfaceRA P a M) :
    LBound (1 / 2) (costExpRA (P.θ : ℝ) a) (fun N => M.Fsup N) := by
  obtain ⟨hθ0, -, -, -, -, -⟩ := mt_adm P hP
  obtain ⟨hss, hs0, hs1, -⟩ := sqrt2_facts
  obtain ⟨C, N₀, hC, hb, hTp, hF0, hmono, hrec⟩ := mt_unpackRA hM
  have hb' := mt_ev_of_ge hb
  set q := (P.θ : ℝ) / (2 * a) with hq
  have hq0 : 0 ≤ q := by positivity
  set pE := max (2 * a) (a + q) with hpE
  have hpa : a ≤ pE := le_max_of_le_left (by linarith)
  have hpE0 : 0 ≤ pE := by linarith
  have hp : ∀ ε > 0, ∀ᶠ N : ℕ in atTop, Params.SplitInput a N → LN (1 / 2) (-ε) N ≤ M.pTrial N :=
    mt_pTrial_lower P hP a ha hCEP hRS hFM hLD M.pTrial M.rej
      (hb'.mono fun N h hQ => ⟨(h hQ).1, (h hQ).2.1, (h hQ).2.2.1⟩)
  have hMr : EvLB (Params.SplitInput a) a M.Mrows :=
    mt_chain_Mrows P a ha.le M.Mrows C (hb'.mono fun N h hQ => (h hQ).2.2.2.2.1)
  have hMr1 : ∀ᶠ N : ℕ in atTop, Params.SplitInput a N → 1 ≤ M.Mrows N :=
    hb'.mono fun N h hQ => (h hQ).2.2.2.1
  have hTsq : EvLB (Params.SplitInput a) (max a q) M.Tsq :=
    mt_chain_Tsq P a ha hθ0.le M.Tsq C (hb'.mono fun N h hQ => (h hQ).2.2.2.2.2.1)
  have hTpo : EvLB (Params.SplitInput a) (2 * a) M.Tpost :=
    mt_chain_Tpost a ha.le M.Mrows M.Tpost C hMr hMr1 hTp
  have htau : EvLB (Params.SplitInput a) a M.tau := by
    refine ((EvLB.const_mul (by norm_num : (0:ℝ) ≤ 3) hMr).div hp).of_le ?_
    filter_upwards [hb'] with N h hQ
    exact (h hQ).2.2.2.2.2.2.2.2.2.1
  have hnonO : EvLB (Params.SplitInput a) pE M.nonO := by
    have h1 := (EvLB.const_mul (by norm_num : (0:ℝ) ≤ 3)
      ((hMr.mul ((EvLB.const_mul (by positivity : (0:ℝ) ≤ 2 * P.k) hTsq).add
        (EvLB.polylog 1 C)) (hMr1.mono fun N h hQ => by linarith [h hQ])).add hTpo)).div hp
    refine (h1.of_le ?_).mono ?_
    · filter_upwards [hb'] with N h hQ
      rw [one_mul]
      exact (h hQ).2.2.2.2.2.2.2.2.2.2
    · refine max_le ?_ (le_max_left _ _)
      rw [← max_add_add_left, ← max_add_add_left, add_zero]
      exact max_le (max_le (le_max_of_le_left (by linarith)) (le_max_right _ _))
        (le_max_of_le_left (by linarith))
  -- Lemma F
  have hD0 : ∀ᶠ x : ℝ in atTop, 0 ≤ 2 * supSplit M.tau a N₀ x :=
    Eventually.of_forall fun x => mul_nonneg (by norm_num)
      (mt_supSplit_nonneg _ _ _ _ fun n hn hs => (hb n hn hs).2.2.2.2.2.2.2.2.1)
  have hPc : ∀ η > 0, ∀ᶠ x : ℝ in atTop,
      C + C * (Params.Bsm a ⌊x⌋₊ : ℝ) * (Real.log x) ^ C + 2 * supSplit M.nonO a N₀ x ≤
        LN (1 / 2) (pE + η) x := by
    intro η hη
    obtain ⟨K1, hK10, hK1⟩ := mt_supSplit_le M.nonO a N₀ (pE + η / 2) (by positivity)
      (hnonO (η / 2) (by linarith))
    filter_upwards [polylog_le_LN (1 / 2) (by norm_num) (by norm_num) (η / 2) (by linarith)
      (2 * C + 2 * K1 + 2) C, Real.tendsto_log_atTop.eventually_ge_atTop 1,
      eventually_ge_atTop 1] with x e1 hlx hx1
    have hlC : 1 ≤ (Real.log x) ^ C := Real.one_le_rpow hlx hC
    have hA1 : 1 ≤ LN (1 / 2) (pE + η / 2) x := one_le_LN _ _ hx1 (by positivity)
    have hfl1 : (1:ℝ) ≤ (⌊x⌋₊ : ℝ) := by exact_mod_cast (Nat.one_le_floor_iff _).2 hx1
    have hBa : (Params.Bsm a ⌊x⌋₊ : ℝ) ≤ LN (1 / 2) (pE + η / 2) x :=
      (mt_Bsm_le a _).trans ((mt_LN_mono_x (1 / 2) a (by norm_num) (by norm_num) ha.le hfl1
        (Nat.floor_le (by linarith))).trans (LN_mono _ _ hx1 (by linarith)))
    have hS := hK1 x
    have e : LN (1 / 2) (pE + η) x = LN (1 / 2) (η / 2) x * LN (1 / 2) (pE + η / 2) x := by
      rw [LN_mul]; congr 1; ring
    rw [e]
    set A := LN (1 / 2) (pE + η / 2) x
    set u := (Real.log x) ^ C
    have hB0 : (0:ℝ) ≤ Params.Bsm a ⌊x⌋₊ := Nat.cast_nonneg _
    have hu : 1 ≤ u * A := by nlinarith
    have s1 : C * (Params.Bsm a ⌊x⌋₊ : ℝ) * u ≤ C * u * A := by
      have := mul_le_mul_of_nonneg_left hBa (mul_nonneg hC (by linarith : (0:ℝ) ≤ u))
      nlinarith
    have s2 : C ≤ C * (u * A) := le_mul_of_one_le_right hC hu
    have s3 : K1 ≤ K1 * (u * A) := le_mul_of_one_le_right hK10 hu
    have s4 : A ≤ u * A := le_mul_of_one_le_left (by linarith) hlC
    have s5 : (2 * C + 2 * K1 + 2) * u * A ≤ LN (1 / 2) (η / 2) x * A := by
      exact mul_le_mul_of_nonneg_right e1 (by linarith)
    nlinarith
  have hDb : ∀ η > 0, ∀ᶠ x : ℝ in atTop, 2 * supSplit M.tau a N₀ x ≤ LN (1 / 2) (a + η) x := by
    intro η hη
    obtain ⟨K2, hK20, hK2⟩ := mt_supSplit_le M.tau a N₀ (a + η / 2) (by positivity)
      (htau (η / 2) (by linarith))
    filter_upwards [polylog_le_LN (1 / 2) (by norm_num) (by norm_num) (η / 2) (by linarith)
      (2 * K2 + 2) 0, eventually_ge_atTop 1] with x e1 hx1
    rw [Real.rpow_zero, mul_one] at e1
    have hA1 : 1 ≤ LN (1 / 2) (a + η / 2) x := one_le_LN _ _ hx1 (by positivity)
    have e : LN (1 / 2) (a + η) x = LN (1 / 2) (η / 2) x * LN (1 / 2) (a + η / 2) x := by
      rw [LN_mul]; congr 1; ring
    rw [e]
    have := hK2 x
    have := mul_le_mul_of_nonneg_right e1 (by linarith : (0:ℝ) ≤ LN (1 / 2) (a + η / 2) x)
    nlinarith
  have key := lemmaF_logfactor (1 / 2) (by norm_num) (by norm_num) M.Fsup
    (fun x => C + C * (Params.Bsm a ⌊x⌋₊ : ℝ) * (Real.log x) ^ C + 2 * supSplit M.nonO a N₀ x)
    (fun x => 2 * supSplit M.tau a N₀ x) hF0 hmono pE a hpE0 ha.le hPc hDb hD0 (by
      filter_upwards [eventually_ge_atTop N₀] with x hx using hrec x hx)
  refine mt_lbound_mono (mt_lbound_of_K (1 / 2) (by norm_num) (by norm_num) M.Fsup _ key) ?_
  have hne : (1:ℝ) - 1 / √2 ≠ 0 := by
    have : 1 / √2 < 1 := by rw [div_lt_one hs0]; nlinarith
    linarith
  have e : a / (1 - (2:ℝ) ^ (-(1 / 2 : ℝ))) = (2 + Real.sqrt 2) * a := by
    rw [mt_rpow_neg_half, div_eq_iff hne]
    field_simp
    have hs2 : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
    linear_combination (-1 : ℝ) * hs2
  rw [e]
  unfold costExpRA
  exact max_le (max_le (le_max_of_le_right (by nlinarith)) (le_max_left _ _)) (le_max_right _ _)

/-- MAIN. [pp] Theorem thm:RA at the optimal `a* = √((√2−1)θ/2)`: expected time
`≤ L_N[1/2, √((1+√2)θ) + o(1)]`. -/
theorem theoremRA_opt (P : Params) (hP : P.Admissible) (hCEP : CEP) (hRS : RosserSchoenfeld)
    (hFM : FirstMoment P (Real.sqrt ((Real.sqrt 2 - 1) * P.θ / 2)))
    (hLD : SecondMoment P (Real.sqrt ((Real.sqrt 2 - 1) * P.θ / 2)))
    (M : CostModelRA P) (hM : CostInterfaceRA P (Real.sqrt ((Real.sqrt 2 - 1) * P.θ / 2)) M) :
    LBound (1 / 2) P.cRAstar (fun N => M.Fsup N) := by
  obtain ⟨hθ0, -⟩ := mt_adm P hP
  obtain ⟨-, -, hs1, -⟩ := sqrt2_facts
  have hx : 0 < (Real.sqrt 2 - 1) * (P.θ : ℝ) / 2 := by
    have : 0 < Real.sqrt 2 - 1 := by linarith
    positivity
  have h := theoremRA P hP _ (Real.sqrt_pos.2 hx) hCEP hRS hFM hLD M hM
  rw [(costExpRA_min _ hθ0).2] at h
  unfold Params.cRAstar
  exact h

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
  obtain ⟨hspec, hlim⟩ := portfolioParams_spec
  refine portfolio_assembly T (fun j N => (M (j + 1)).Fsup N)
    (fun j => (portfolioParams (j + 1)).cRAstar) ?_ ?_ ?_
  · have h1 : Tendsto (fun j => ((portfolioParams (j + 1)).θ : ℝ)) atTop (nhds 0) :=
      hlim.comp (tendsto_add_atTop_nat 1)
    have h2 := (h1.const_mul (1 + Real.sqrt 2)).sqrt
    simpa [Params.cRAstar] using h2
  · intro j
    exact theoremRA_opt _ (hspec (j + 1) (by omega)).1 hCEP hRS (hFM _ (by omega))
      (hLD _ (by omega)) (M (j + 1)) (hM _ (by omega))
  · intro j
    exact hport (j + 1) (by omega)

/-! ## Theorem U and Theorem DL: assembly steps only -/

/-- MAIN (ASSEMBLY STEP ONLY; Theorem U itself is NOT formalised). lv5unif §10 (Theorem U), assembly step: the recursion inequality (10.1),
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
  have key := lemmaF_logfactor (2 / 5) (by norm_num) (by norm_num) F Pc (fun x => D x / 2) hF0
    hmono (max (2 * a) (a + s + max a σ)) (a + s) (le_max_of_le_left (by linarith)) (by linarith)
    hP (fun η hη => by
      filter_upwards [hD η hη, hD0] with x h1 h2
      linarith) (by
      filter_upwards [hD0] with x h
      linarith) (by
      filter_upwards [hrec] with x h
      have e : 2 * (D x / 2) * F (Real.sqrt (2 * x)) = D x * F (Real.sqrt (2 * x)) := by ring
      rw [e]; exact h)
  exact mt_lbound_of_K (2 / 5) (by norm_num) (by norm_num) F _ key

/-- MAIN (ASSEMBLY STEP ONLY; Theorem DL itself is NOT formalised; its relation supply at prime modulus is not covered by `FirstMoment`/`SecondMoment`, which are restricted to split inputs). [pp] Theorem thm:DL, assembly step: with `B = L^a`, if the expected time is
`≤ (M_rows C_trial + T_post)/p` with `M_rows ≤ L^{a+o(1)}`, `C_trial ≤ L^{max(1/√2, θ/(2a))+o(1)}`,
`T_post ≤ L^{3a+o(1)}` and `p ≥ L^{−o(1)}` (all INTERNAL inputs, not formalised here), then the
expected time is `≤ L_p[1/2, max(a + 1/√2, a + θ/(2a), 3a) + o(1)]`; at the optimal `a` this is
`c_DL(θ)` (`costExpDL_min`). -/
theorem theoremDL_assembly (θ a : ℝ) (ha : 0 < a) (T Mrows Ctrial Tpost pT : ℕ → ℝ)
    (hT : ∃ K : ℝ, 0 ≤ K ∧ ∀ᶠ N : ℕ in atTop,
      T N ≤ K * (Real.log N) ^ K * (Mrows N * Ctrial N + Tpost N) / pT N)
    (hM : LBound (1 / 2) a Mrows) (hM0 : ∀ᶠ N : ℕ in atTop, 0 ≤ Mrows N)
    (hC : LBound (1 / 2) (max a (max (1 / Real.sqrt 2) (θ / (2 * a)))) Ctrial)
    (hC0 : ∀ᶠ N : ℕ in atTop, 0 ≤ Ctrial N)
    (hTp : LBound (1 / 2) (3 * a) Tpost)
    (hp : ∀ ε > 0, ∀ᶠ N : ℕ in atTop, LN (1 / 2) (-ε) N ≤ pT N) :
    LBound (1 / 2) (costExpDL θ a) T := by
  obtain ⟨K, hK0, hT⟩ := hT
  have hM' : EvLB (fun _ => True) a Mrows := EvLB.of_lbound hM
  have hC' : EvLB (fun _ => True) (max a (max (1 / Real.sqrt 2) (θ / (2 * a)))) Ctrial :=
    EvLB.of_lbound hC
  have hTp' : EvLB (fun _ => True) (3 * a) Tpost := EvLB.of_lbound hTp
  have h1 := (hM'.mul hC' (hM0.mono fun N h _ => h)).add hTp'
  have hpl0 : ∀ᶠ N : ℕ in atTop, True → 0 ≤ K * (Real.log N) ^ K := by
    filter_upwards [eventually_ge_atTop 1] with N hN _
    exact mul_nonneg hK0 (Real.rpow_nonneg (Real.log_nonneg (by exact_mod_cast hN)) _)
  have h2 := ((EvLB.polylog K K).mul h1 hpl0).div (p := pT)
    (fun ε hε => (hp ε hε).mono fun N h _ => h)
  refine EvLB.toLBound ((h2.of_le ?_).mono ?_)
  · filter_upwards [hT] with N h _ using h
  · unfold costExpDL
    rw [zero_add]
    refine max_le ?_ (le_max_of_le_right (le_max_right _ _))
    rw [← max_add_add_left, ← max_add_add_left]
    exact max_le (le_max_of_le_right (le_max_of_le_right (by linarith)))
      (max_le (le_max_left _ _) (le_max_of_le_right (le_max_left _ _)))

end FactoringA
