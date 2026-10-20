import LeanA.DefsNFS
import LeanA.SplitCounting
import LeanA.Moments
import LeanA.Asymptotics
import LeanA.Hypotheses

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

/-! ### CEP in the range of [pp] Lemma lem:local(iv) -/

lemma rho_aux_fin {u lt Λ t a c K ε l2 : ℝ} (ha : 0 < a) (hΛ0 : 0 ≤ Λ)
    (hlt0 : 0 ≤ lt) (hsq : Λ ^ 2 = t * lt) (hK1 : K * lt ≤ ε * a / 2 * Λ)
    (hK2 : (c / a + ε) * l2 ≤ ε * a / 2 * Λ) (hpos : 0 < a * Λ - l2)
    (hu1 : u * (a * Λ - l2) ≤ c * t + K * Λ) : u * lt ≤ (c / a + ε) * Λ := by
  by_contra hcon
  push Not at hcon
  have h1 : (c / a + ε) * Λ * (a * Λ - l2) < u * lt * (a * Λ - l2) :=
    mul_lt_mul_of_pos_right hcon hpos
  have h2 : u * lt * (a * Λ - l2) ≤ (c * t + K * Λ) * lt := by
    have := mul_le_mul_of_nonneg_right hu1 hlt0; linarith
  have e1 : (c / a + ε) * Λ * (a * Λ - l2) =
      c * Λ ^ 2 + ε * a * Λ ^ 2 - (c / a + ε) * l2 * Λ := by field_simp
  have e2 : (c * t + K * Λ) * lt = c * Λ ^ 2 + K * lt * Λ := by rw [hsq]; ring
  rw [e1] at h1; rw [e2] at h2
  have h3 : K * lt * Λ ≤ ε * a / 2 * Λ * Λ := mul_le_mul_of_nonneg_right hK1 hΛ0
  have h4 : (c / a + ε) * l2 * Λ ≤ ε * a / 2 * Λ * Λ := mul_le_mul_of_nonneg_right hK2 hΛ0
  have e3 : ε * a / 2 * Λ * Λ = ε * a * Λ ^ 2 / 2 := by ring
  linarith

lemma rho_aux_end {u lt Λ δ y a c ε : ℝ} (hu0 : 0 < u) (hδ0 : 0 < δ)
    (hlogu : Real.log u ≤ 1 / 2 * lt) (hfin : u * lt ≤ (c / a + ε) * Λ)
    (hδy : (1 + δ) * y = c / (2 * a) + ε) (hy : y = c / (2 * a) + ε / 2) :
    -((c / (2 * a) + ε) * Λ) ≤ Real.log u * -(u * (1 + δ)) := by
  have h1 : Real.log u * (u * (1 + δ)) ≤ (1 / 2 * lt) * (u * (1 + δ)) :=
    mul_le_mul_of_nonneg_right hlogu (by positivity)
  have e3 : (1 / 2 * lt) * (u * (1 + δ)) = (1 + δ) * (1 / 2 * (u * lt)) := by ring
  have e4 : (c / (2 * a) + ε) * Λ = (1 + δ) * (y * Λ) := by rw [← hδy]; ring
  have e5 : 1 / 2 * (u * lt) ≤ y * Λ := by
    have : y * Λ = 1 / 2 * ((c / a + ε) * Λ) := by rw [hy]; field_simp
    rw [this]; linarith
  have h2 : (1 + δ) * (1 / 2 * (u * lt)) ≤ (1 + δ) * (y * Λ) :=
    mul_le_mul_of_nonneg_left e5 (by positivity)
  have e6 : Real.log u * -(u * (1 + δ)) = -(Real.log u * (u * (1 + δ))) := by ring
  rw [e6, e4]; linarith

/-- Lower bound for `ϱ(X, B)` from CEP when `log X = c log N + O(Λ)` and `log B = aΛ + O(1)`,
`Λ = Lam (1/2) N`: `ϱ(X, B) ≥ L^{−c/(2a) − ε}`. -/
lemma rho_lower (hCEP : CEP) {a c K ε : ℝ} (ha : 0 < a) (hc : 0 < c) (hK : 0 ≤ K) (hε : 0 < ε) :
    ∀ᶠ N : ℕ in atTop, ∀ (X : ℝ) (B : ℕ), 0 < X →
      c / 2 * Real.log N ≤ Real.log X →
      Real.log X ≤ c * Real.log N + K * Lam (1 / 2) (N : ℝ) →
      a * Lam (1 / 2) (N : ℝ) - Real.log 2 ≤ Real.log B → Real.log B ≤ a * Lam (1 / 2) (N : ℝ) →
      Real.exp (-((c / (2 * a) + ε) * Lam (1 / 2) (N : ℝ))) ≤ rho X B := by
  set y := c / (2 * a) + ε / 2 with hy
  have hy0 : 0 < y := by positivity
  set δ := (ε / 2) / y with hδ
  have hδ0 : 0 < δ := by positivity
  have hδy : (1 + δ) * y = c / (2 * a) + ε := by
    have : (1 + δ) * y = y + ε / 2 := by rw [hδ]; field_simp
    rw [this, hy]; ring
  obtain ⟨U₀, hU⟩ := hCEP (1 / 2) (by norm_num) δ hδ0
  set U := max U₀ 1 with hUdef
  have hU1 : 1 ≤ U := le_max_right _ _
  have hl2 : 0 < Real.log 2 := Real.log_pos (by norm_num)
  filter_upwards [ev_Lam_sq, ev_log_ge 1, ev_log_ge (3 / c), ev_loglog_ge 0,
    ev_Lam_le 1 one_pos, ev_Lam_le (c / (2 * a * U)) (by positivity),
    ev_Lam_ge (2 * Real.log 2 / a),
    ev_loglog_le (a / 4) (by positivity),
    ev_Lam_ge (4 / a * (|Real.log (c + K)| + Real.log 2)),
    ev_loglog_ge (4 * (c + K) ^ 2 / a ^ 2),
    ev_loglog_le (ε * a / (2 * (K + 1))) (by positivity),
    ev_Lam_ge (2 * (c / a + ε) * Real.log 2 / (ε * a)), ev_log_ge (2 * Real.log 16 / c)] with
    N hsq ht1 ht3 hlt0 hΛt hΛU hΛ2 hlt1 hΛ3 hlt2 hlt3 hΛ4 ht16
  intro X B hX hX1 hX2 hB1 hB2
  set t := Real.log (N : ℝ) with ht
  set Λ := Lam (1 / 2) (N : ℝ) with hΛ
  set lt := Real.log t with hlt
  have hΛ0 : 0 ≤ Λ := Lam_half_nonneg _
  have hΛ2' : 2 * Real.log 2 ≤ a * Λ := by
    have h := hΛ2
    rw [div_le_iff₀ ha] at h; linarith
  have hBΛ : a * Λ / 2 ≤ Real.log B := by linarith
  have hB0 : 0 < Real.log B := by linarith
  have hBpos : (0 : ℝ) < B := by
    rcases Nat.eq_zero_or_pos B with h | h
    · rw [h] at hB0; simp at hB0
    · exact_mod_cast h
  have hlX : 1 < Real.log X := by
    have : c / 2 * (3 / c) ≤ c / 2 * t := by gcongr
    have e : c / 2 * (3 / c) = 3 / 2 := by field_simp
    linarith
  have hllX : 0 < Real.log (Real.log X) := Real.log_pos hlX
  have hXt : Real.log X ≤ (c + K) * t := by
    have : K * Λ ≤ K * t := mul_le_mul_of_nonneg_left (by linarith) hK
    linarith
  have hcK : 0 < c + K := by linarith
  have ht0 : 0 < t := by linarith
  have hllX2 : Real.log (Real.log X) ≤ Real.log (c + K) + lt := by
    rw [← Real.log_mul hcK.ne' ht0.ne']
    exact Real.log_le_log (by linarith) hXt
  have hΛ3' : 4 * (|Real.log (c + K)| + Real.log 2) ≤ a * Λ := by
    have h := hΛ3
    rw [div_mul_eq_mul_div, div_le_iff₀ ha] at h
    linarith
  have habs := le_abs_self (Real.log (c + K))
  have hllB : Real.log (Real.log X) ≤ 1 / 2 * Real.log B := by
    have : lt ≤ a / 4 * Λ := hlt1
    linarith
  set u := Real.log X / Real.log B with hu
  have hu0 : 0 < u := div_pos (by linarith) hB0
  have huB : u * Real.log B = Real.log X := by rw [hu]; field_simp
  have hUu : U₀ ≤ u := by
    have h1 : U * Real.log B ≤ Real.log X := by
      have h2 : U * (a * Λ) ≤ c / 2 * t := by
        have h3 := mul_le_mul_of_nonneg_left hΛU (by positivity : 0 ≤ U * a)
        have e : U * a * (c / (2 * a * U) * t) = c / 2 * t := by field_simp
        rw [e] at h3; linarith
      have h4 := mul_le_mul_of_nonneg_left hB2 (by positivity : 0 ≤ U)
      linarith
    have : U ≤ u := by rw [hu, le_div_iff₀ hB0]; exact h1
    exact le_trans (le_max_left _ _) this
  have hurange : u ≤ (1 - 1 / 2) * Real.log X / Real.log (Real.log X) := by
    rw [hu, div_le_div_iff₀ hB0 hllX]
    have := mul_le_mul_of_nonneg_left hllB (by linarith : 0 ≤ Real.log X)
    linarith
  have hX16 : 16 ≤ X := by
    have h1 : Real.log 16 ≤ Real.log X := by
      have h := ht16
      rw [div_le_iff₀ hc] at h
      linarith
    exact (Real.log_le_log_iff (by norm_num) hX).mp h1
  have hcep := (hU X u hX16 hUu hurange).1
  have hXu : X ^ (1 / u) = (B : ℝ) := by
    rw [Real.rpow_def_of_pos hX]
    have : Real.log X * (1 / u) = Real.log B := by
      rw [hu]; field_simp
    rw [this, Real.exp_log hBpos]
  rw [hXu] at hcep
  have hrho : u ^ (-(u * (1 + δ))) ≤ rho X B := by
    unfold rho
    rw [le_div_iff₀ hX]; linarith
  refine le_trans ?_ hrho
  rw [Real.rpow_def_of_pos hu0]
  apply Real.exp_le_exp.mpr
  -- `u (aΛ − log 2) ≤ c t + KΛ`
  have hu1 : u * (a * Λ - Real.log 2) ≤ c * t + K * Λ := by
    have := mul_le_mul_of_nonneg_left hB1 hu0.le
    linarith
  -- `u² ≤ t`
  have husq : u ^ 2 ≤ t := by
    have h1 : u * (a * Λ) ≤ 2 * ((c + K) * t) := by
      have := mul_le_mul_of_nonneg_left hBΛ hu0.le
      have e : u * (a * Λ / 2) = u * (a * Λ) / 2 := by ring
      linarith
    have h2 : (u * (a * Λ)) ^ 2 ≤ (2 * ((c + K) * t)) ^ 2 :=
      pow_le_pow_left₀ (by positivity) h1 2
    have hsq' : Λ ^ 2 = t * lt := hsq
    have h3 : 4 * (c + K) ^ 2 ≤ a ^ 2 * lt := by
      have h := hlt2
      rw [div_le_iff₀ (by positivity)] at h; linarith
    by_contra hcon
    push Not at hcon
    have h4 : (u * (a * Λ)) ^ 2 = u ^ 2 * (a ^ 2 * lt) * t := by
      rw [mul_pow, mul_pow, hsq']; ring
    have h5 : t * (a ^ 2 * lt) * t ≤ u ^ 2 * (a ^ 2 * lt) * t := by
      have hal : 0 ≤ a ^ 2 * lt := by positivity
      have := mul_le_mul_of_nonneg_right hcon.le hal
      exact mul_le_mul_of_nonneg_right this ht0.le
    have h6 : t * (4 * (c + K) ^ 2) * t ≤ t * (a ^ 2 * lt) * t := by
      have := mul_le_mul_of_nonneg_left h3 ht0.le
      exact mul_le_mul_of_nonneg_right this ht0.le
    have h7 : (2 * ((c + K) * t)) ^ 2 = t * (4 * (c + K) ^ 2) * t := by ring
    -- all the inequalities are equalities; contradiction with strictness
    have h8 : t * (a ^ 2 * lt) * t < u ^ 2 * (a ^ 2 * lt) * t := by
      have hal : 0 < a ^ 2 * lt := by
        have : 0 < 4 * (c + K) ^ 2 := by positivity
        linarith
      have := mul_lt_mul_of_pos_right hcon hal
      exact mul_lt_mul_of_pos_right this ht0
    linarith
  have hlogu : Real.log u ≤ 1 / 2 * lt := by
    have : Real.log (u ^ 2) ≤ lt := Real.log_le_log (by positivity) husq
    rw [Real.log_pow] at this; push_cast at this; linarith
  -- `u log t ≤ (c/a + ε) Λ`
  have hK1 : K * lt ≤ ε * a / 2 * Λ := by
    have h1 : K * lt ≤ K * (ε * a / (2 * (K + 1)) * Λ) := mul_le_mul_of_nonneg_left hlt3 hK
    have e : K * (ε * a / (2 * (K + 1)) * Λ) = (K / (K + 1)) * (ε * a / 2 * Λ) := by
      field_simp
    have h2 : K / (K + 1) ≤ 1 := by rw [div_le_one (by linarith)]; linarith
    have h3 : 0 ≤ ε * a / 2 * Λ := by positivity
    have h4 := mul_le_mul_of_nonneg_right h2 h3
    rw [e] at h1; linarith
  have hK2 : (c / a + ε) * Real.log 2 ≤ ε * a / 2 * Λ := by
    have h := hΛ4
    rw [div_le_iff₀ (by positivity)] at h; linarith
  have hfin := rho_aux_fin ha hΛ0 hlt0 hsq hK1 hK2 (by linarith) hu1
  exact rho_aux_end hu0 hδ0 hlogu hfin hδy hy

/-- MAIN. [pp] Lemma lem:omega: if `N > 1` has no prime factor `≤ B` (`B ≥ 2`), then
`ω(N) log B ≤ log N`. -/
theorem omega_le (N B : ℕ) (hN : 1 < N) (hB : 2 ≤ B) (hNB : ∀ p : ℕ, p.Prime → p ∣ N → B < p) :
    (omega N : ℝ) * Real.log B ≤ Real.log N := by
  have hB0 : (0 : ℝ) < B := by exact_mod_cast (by omega : 0 < B)
  have hpow : B ^ omega N ≤ N := by
    have h1 : B ^ N.primeFactors.card ≤ ∏ p ∈ N.primeFactors, p :=
      Finset.pow_card_le_prod _ _ _ (fun p hp =>
        (hNB p (Nat.prime_of_mem_primeFactors hp) (Nat.dvd_of_mem_primeFactors hp)).le)
    have h2 : ∏ p ∈ N.primeFactors, p ≤ N := Nat.le_of_dvd (by omega) (Nat.prod_primeFactors_dvd N)
    exact le_trans h1 h2
  have h : ((B : ℝ)) ^ omega N ≤ N := by exact_mod_cast hpow
  calc (omega N : ℝ) * Real.log B = Real.log ((B : ℝ) ^ omega N) := by rw [Real.log_pow]
    _ ≤ Real.log N := Real.log_le_log (by positivity) h

/-- MAIN. [pp] Theorem thm:Dk (relation supply): from the two moment hypotheses and CEP,
`γ ≥ L^{−o(1)}` uniformly over split inputs. -/
theorem relation_supply (P : Params) (hP : P.Admissible) (a : ℝ) (ha : 0 < a) (hCEP : CEP)
    (hFM : FirstMoment P a) (hLD : SecondMoment P a) :
    ∀ ε > 0, ∃ N₀ : ℕ, ∀ N : ℕ, N₀ ≤ N → Params.SplitInput a N →
      LN (1 / 2) (-ε) N ≤ P.gammaD a N := by
  intro ε hε
  obtain ⟨hν, hνk, hμν, hμ, hh, hθ, hμ1, hd⟩ := Params.adm_real hP
  obtain ⟨C1, N1, hFM'⟩ := hFM
  obtain ⟨C2, N2, hLD'⟩ := hLD
  set C := max (max C1 C2) 0 with hC
  have hC0 : 0 ≤ C := le_max_right _ _
  have hC1 : C1 ≤ C := le_trans (le_max_left _ _) (le_max_left _ _)
  have hC2 : C2 ≤ C := le_trans (le_max_right _ _) (le_max_left _ _)
  set ε₁ := ε / 8 with hε₁
  have hε₁0 : 0 < ε₁ := by positivity
  set cr : ℝ := (P.μ : ℝ) - P.ν with hcr
  set cF : ℝ := (P.μ : ℝ) + P.h + P.d * P.ν with hcF
  have hcF0 : 0 < cF := by positivity
  have hθc : (P.θ : ℝ) = cr + cF := by rw [Params.theta_real]; ring
  have key : ∀ᶠ N : ℕ in atTop, Params.SplitInput a N → LN (1 / 2) (-ε) N ≤ P.gammaD a N := by
    filter_upwards [rho_lower hCEP ha hμν (K := (P.θ : ℝ) / (4 * a) + 1) (by positivity) hε₁0,
      rho_lower hCEP ha hcF0 (K := P.d * ((P.θ : ℝ) / (4 * a)) + 1) (by positivity) hε₁0,
      ev_Bsm ha, ev_Xr hP ha, ev_XF hP ha, ev_nat_ge N1, ev_nat_ge N2, ev_log_ge 1,
      ev_loglog_le (ε₁ / (C + 1)) (by positivity), ev_Lam_ge (8 * Real.log 2 / ε)] with
      N hr hF hB hXr hXF hN1 hN2 hlog1 hll hΛ
    intro hsplit
    set Λ := Lam (1 / 2) (N : ℝ) with hΛdef
    have hΛ0 : 0 ≤ Λ := Lam_half_nonneg _
    have hρr := hr (P.Xr a N) (Params.Bsm a N) hXr.1 hXr.2.1 hXr.2.2 hB.2.1 hB.2.2
    have hρF := hF (P.XF a N) (Params.Bsm a N) hXF.1 hXF.2.1 hXF.2.2 hB.2.1 hB.2.2
    obtain ⟨-, hS2, -, -⟩ := Ssz_bounds P ha hθ.le N
    set m : ℝ := (P.Ssz a N : ℝ) ^ 2 * P.rhoR a N * P.rhoF a N with hm
    have hm1 : Real.exp (-(2 * ε₁) * Λ) ≤ m := by
      have e : Real.exp (-(2 * ε₁) * Λ) = Real.exp ((P.θ : ℝ) / (2 * a) * Λ) *
          Real.exp (-((cr / (2 * a) + ε₁) * Λ)) * Real.exp (-((cF / (2 * a) + ε₁) * Λ)) := by
        rw [← Real.exp_add, ← Real.exp_add, hθc]; congr 1; field_simp; ring
      rw [e, hm]
      unfold Params.rhoR Params.rhoF
      apply mul_le_mul (mul_le_mul hS2 hρr (by positivity) (by positivity)) hρF
        (by positivity)
      have : 0 ≤ (P.Ssz a N : ℝ) ^ 2 := by positivity
      have : 0 ≤ rho (P.Xr a N) (Params.Bsm a N) := le_trans (by positivity) hρr
      positivity
    have hm0 : 0 < m := lt_of_lt_of_le (Real.exp_pos _) hm1
    set Pw := (Real.log N) ^ C with hPw
    have hlog0 : 0 < Real.log N := by linarith
    have hPw1 : 1 ≤ Pw := Real.one_le_rpow hlog1 hC0
    have hPw0 : 0 < Pw := by linarith
    have hPwe : Pw ≤ Real.exp (ε₁ * Λ) := by
      rw [hPw, Real.rpow_def_of_pos hlog0]
      apply Real.exp_le_exp.mpr
      have hll0 : 0 ≤ Real.log (Real.log N) := Real.log_nonneg hlog1
      have h1 : Real.log (Real.log N) * C ≤ (ε₁ / (C + 1) * Λ) * C :=
        mul_le_mul_of_nonneg_right hll hC0
      have h2 : (ε₁ / (C + 1) * Λ) * C = (C / (C + 1)) * (ε₁ * Λ) := by field_simp
      have h3 : C / (C + 1) ≤ 1 := by rw [div_le_one (by linarith)]; linarith
      have h4 := mul_le_mul_of_nonneg_right h3 (by positivity : 0 ≤ ε₁ * Λ)
      linarith
    have h1 : m / Pw ≤ expect (P.fields N ×ˢ P.Jset N) (fun x => (P.countN a N x.1 x.2 : ℝ)) := by
      refine le_trans ?_ (hFM' N hN1 hsplit)
      exact div_le_div_of_nonneg_left hm0.le (Real.rpow_pos_of_pos hlog0 _)
        (Real.rpow_le_rpow_of_exponent_le hlog1 hC1)
    have h2 : expect (P.fields N ×ˢ P.Jset N)
        (fun x => (P.countN a N x.1 x.2 : ℝ) * ((P.countN a N x.1 x.2 : ℝ) - 1)) ≤ Pw * m ^ 2 := by
      refine le_trans (hLD' N hN2 hsplit) ?_
      have e : (Real.log N) ^ C2 * (P.Ssz a N : ℝ) ^ 4 * (P.rhoR a N) ^ 2 * (P.rhoF a N) ^ 2 =
          (Real.log N) ^ C2 * m ^ 2 := by rw [hm]; ring
      rw [e]
      exact mul_le_mul_of_nonneg_right (Real.rpow_le_rpow_of_exponent_le hlog1 hC2)
        (by positivity)
    have hγ := relation_supply_from_moments (P.fields N ×ˢ P.Jset N)
      (fun x => P.countN a N x.1 x.2) m Pw hm0 hPw0 h1 h2
    unfold Params.gammaD
    refine le_trans ?_ hγ
    -- `Pw³ + Pw/m ≤ 2 exp(3ε₁Λ) ≤ exp(εΛ)`
    have hQ1 : Pw ^ 3 ≤ Real.exp (3 * ε₁ * Λ) := by
      have := pow_le_pow_left₀ hPw0.le hPwe 3
      rw [← Real.exp_nat_mul] at this
      convert this using 2; push_cast; ring
    have hQ2 : Pw / m ≤ Real.exp (3 * ε₁ * Λ) := by
      rw [div_le_iff₀ hm0]
      calc Pw ≤ Real.exp (ε₁ * Λ) := hPwe
        _ = Real.exp (3 * ε₁ * Λ) * Real.exp (-(2 * ε₁) * Λ) := by
          rw [← Real.exp_add]; congr 1; ring
        _ ≤ Real.exp (3 * ε₁ * Λ) * m := mul_le_mul_of_nonneg_left hm1 (by positivity)
    have hl2 : Real.log 2 ≤ ε / 8 * Λ := by
      have h := hΛ
      rw [div_le_iff₀ hε] at h; linarith
    have hQ : Pw ^ 3 + Pw / m ≤ Real.exp (ε * Λ) := by
      have e2 : 2 * Real.exp (3 * ε₁ * Λ) = Real.exp (Real.log 2 + 3 * ε₁ * Λ) := by
        rw [Real.exp_add, Real.exp_log (by norm_num)]
      have h3 : Real.log 2 + 3 * ε₁ * Λ ≤ ε * Λ := by
        have := mul_nonneg hε.le hΛ0
        rw [hε₁]; linarith
      have h4 := Real.exp_le_exp.mpr h3
      linarith
    unfold LN
    rw [show -ε * Lam (1 / 2) (N : ℝ) = -(ε * Λ) by ring, Real.exp_neg, ← one_div]
    exact one_div_le_one_div_of_le (by positivity) hQ
  obtain ⟨N₀, hN₀⟩ := eventually_atTop.mp key
  exact ⟨N₀, fun N hN hs => hN₀ N hN hs⟩

/-! ### Field selection: combinatorics -/

/-- `E_A[|𝒟_A|/|J|] = γ`. -/
lemma expect_descPrimes (P : Params) (a : ℝ) (N : ℕ) :
    expect (P.fields N) (fun A => ((P.descPrimes a N A).card : ℝ) / (P.Jset N).card) =
      P.gammaD a N := by
  unfold expect Params.descPrimes
  rw [Params.gammaD, prob_eq_div]
  have hc : (((P.fields N ×ˢ P.Jset N).filter (fun x => 1 ≤ P.countN a N x.1 x.2)).card : ℝ) =
      ∑ A ∈ P.fields N, (((P.Jset N).filter (fun ℓ => 1 ≤ P.countN a N A ℓ)).card : ℝ) := by
    rw [Finset.natCast_card_filter, Finset.sum_product]
    apply Finset.sum_congr rfl; intro A _
    rw [Finset.natCast_card_filter]
  rw [hc, ← Finset.sum_div, div_div, Finset.card_product, Nat.cast_mul, mul_comm]

/-- Field selection ([pp] Prop. prop:select(i)): if `δ ∈ [0, 1]` has mean `γ > 0`, `Pr(𝓡) ≤ γ/2`,
and `dP ≥ Q` whenever `A ∉ 𝓡` and `δ_A ≥ γ/4`, then `E[dP · 1_{𝓡ᶜ}] ≥ (γ/4) Q`. -/
lemma select_bound {α : Type*} [DecidableEq α] (F Rej : Finset α) (δ dP : α → ℝ) (γ Q : ℝ) (hγ : 0 < γ)
    (hδ : ∀ A ∈ F, 0 ≤ δ A ∧ δ A ≤ 1) (hexp : expect F δ = γ)
    (hRej : prob F (fun A => A ∈ Rej) ≤ γ / 2) (hdP0 : ∀ A ∈ F, 0 ≤ dP A)
    (hQ : ∀ A ∈ F, A ∉ Rej → γ / 4 ≤ δ A → Q ≤ dP A) (hQ0 : 0 ≤ Q) :
    γ / 4 * Q ≤ expect F (fun A => if A ∈ Rej then 0 else dP A) := by
  set Z : α → ℝ := fun A => if A ∈ Rej then 0 else δ A with hZdef
  have hZ : ∀ A ∈ F, 0 ≤ Z A ∧ Z A ≤ 1 := by
    intro A hA; obtain ⟨h0, h1⟩ := hδ A hA
    simp only [hZdef]; split_ifs <;> constructor <;> linarith
  have hF : (0 : ℝ) < F.card := by
    rcases F.eq_empty_or_nonempty with rfl | hne
    · simp [expect] at hexp; linarith
    · exact_mod_cast hne.card_pos
  have hEZ : γ / 2 ≤ expect F Z := by
    have hsum : ∑ A ∈ F, δ A ≤ ∑ A ∈ F, Z A + ∑ A ∈ F, (if A ∈ Rej then (1 : ℝ) else 0) := by
      rw [← Finset.sum_add_distrib]
      apply Finset.sum_le_sum; intro A hA
      obtain ⟨h0, h1⟩ := hδ A hA
      simp only [hZdef]; split_ifs <;> linarith
    rw [Finset.sum_boole] at hsum
    unfold expect at hexp ⊢; rw [prob_eq_div] at hRej
    rw [div_eq_iff hF.ne'] at hexp
    rw [div_le_iff₀ hF] at hRej
    rw [le_div_iff₀ hF]
    linarith
  have hrm := reverse_markov F Z hZ (γ / 2) hEZ
  have hsum : ∑ A ∈ F, (if γ / 2 / 2 ≤ Z A then Q else 0) ≤
      ∑ A ∈ F, (if A ∈ Rej then 0 else dP A) := by
    apply Finset.sum_le_sum; intro A hA
    by_cases hc : γ / 2 / 2 ≤ Z A
    · rw [ite_eq_left hc]
      by_cases hR : A ∈ Rej
      · simp only [hZdef, ite_eq_left hR] at hc; linarith
      · rw [ite_eq_right hR]
        simp only [hZdef, ite_eq_right hR] at hc
        exact hQ A hA hR (by linarith)
    · rw [ite_eq_right hc]; split_ifs
      · exact le_refl _
      · exact hdP0 A hA
  rw [← Finset.sum_filter, Finset.sum_const, nsmul_eq_mul] at hsum
  rw [prob_eq_div] at hrm
  unfold expect
  rw [le_div_iff₀ hF]
  rw [le_div_iff₀ hF] at hrm
  have : γ / 2 / 2 * F.card * Q ≤ ((F.filter (fun A => γ / 2 / 2 ≤ Z A)).card : ℝ) * Q :=
    mul_le_mul_of_nonneg_right hrm hQ0
  have e : γ / 4 * Q * F.card = γ / 2 / 2 * F.card * Q := by ring
  linarith

/-- `|J| ≥ 3Y/(5 log Y)` (Rosser–Schoenfeld). -/
lemma card_primesIoc_ge (hRS : RosserSchoenfeld) (Y : ℕ) (hY : 21 ≤ Y) :
    3 * (Y : ℝ) / (5 * Real.log Y) ≤ (primesIoc Y (2 * Y)).card := by
  have h := hRS.2 (Y : ℝ) (by have : (21 : ℝ) ≤ Y := by exact_mod_cast hY
                              linarith)
  rw [show (2 : ℝ) * (Y : ℝ) = ((2 * Y : ℕ) : ℝ) by push_cast; ring, Nat.floor_natCast,
    Nat.floor_natCast] at h
  have hc := card_filter_Ioc_add_range Nat.Prime (by omega : Y ≤ 2 * Y)
  rw [← primeCounting_eq_card, ← primeCounting_eq_card] at hc
  have : ((primesIoc Y (2 * Y)).card : ℝ) =
      (Nat.primeCounting (2 * Y) : ℝ) - Nat.primeCounting Y := by
    unfold primesIoc
    rw [← hc]; push_cast; ring
  rw [this]; linarith

/-! ### Field selection: arithmetic -/

lemma omega_aux {w L t lt a e l2 lB : ℝ} (hw : 0 ≤ w) (hL : 0 < L) (ha : 0 < a) (he : 0 < e)
    (hl2 : 0 < l2) (hsq : L ^ 2 = t * lt) (ht : 0 ≤ t) (hB : a * L / 2 ≤ lB) (hwB : w * lB ≤ t)
    (hlt : 8 * l2 / (a * e) ≤ lt) : w * l2 ≤ e / 4 * L := by
  have h1 : w * (a * L / 2) ≤ t := le_trans (mul_le_mul_of_nonneg_left hB hw) hwB
  have h2 : 8 * l2 ≤ a * e * lt := by rw [div_le_iff₀ (by positivity)] at hlt; linarith
  -- `w L l2 ≤ (2t/a) l2 ≤ (e/4) t lt = (e/4) L²`
  have h3 : w * L * l2 * a ≤ 2 * t * l2 := by
    have := mul_le_mul_of_nonneg_right h1 hl2.le; nlinarith
  have h4 : 2 * t * l2 * 4 ≤ e * L ^ 2 * a := by
    rw [hsq]; have := mul_le_mul_of_nonneg_left h2 ht; nlinarith
  have h5 : w * l2 * L * a * 4 ≤ e / 4 * L * L * a * 4 := by nlinarith
  have h6 : 0 < L * a * 4 := by positivity
  have : w * l2 * (L * a * 4) ≤ e / 4 * L * (L * a * 4) := by nlinarith
  exact le_of_mul_le_mul_right this h6

lemma descent_arith {lF lγ lJ lY llY ω2 t Λ kR ε ε₁ ε₂ Kc l3 l4 l5 lk lc l2 : ℝ} (hkR : 0 ≤ kR)
    (hlF : lF = lγ - l4 + 2 * kR * (lγ + lJ - l4 - lk) - (lc + ω2 + t))
    (hγ : -(ε₁ * Λ) ≤ lγ) (h1 : (2 * kR + 1) * ε₁ = ε / 4) (hJ : l3 + lY - l5 - llY ≤ lJ)
    (hY : t - 6 * kR * l2 ≤ 2 * kR * lY) (hllY : llY ≤ ε₂ * Λ) (h2 : 2 * kR * ε₂ = ε / 4)
    (hω : ω2 ≤ ε / 4 * Λ)
    (hKc : Kc = l4 + 2 * kR * l4 + 2 * kR * lk + 2 * kR * (l5 - l3) + lc + 6 * kR * l2)
    (hK : Kc ≤ ε / 4 * Λ) : -(ε * Λ) ≤ lF := by
  have a1 := mul_le_mul_of_nonneg_left hγ (by linarith : 0 ≤ 2 * kR + 1)
  have a2 := mul_le_mul_of_nonneg_left hJ (by linarith : 0 ≤ 2 * kR)
  have a3 := mul_le_mul_of_nonneg_left hllY (by linarith : 0 ≤ 2 * kR)
  have e1 : (2 * kR + 1) * -(ε₁ * Λ) = -(ε / 4 * Λ) := by rw [← h1]; ring
  have e2 : 2 * kR * (ε₂ * Λ) = ε / 4 * Λ := by rw [← h2]; ring
  rw [e1] at a1; rw [e2] at a3
  rw [hlF]; rw [hKc] at hK
  nlinarith

/-- `2X² < N` for odd `N` ([pp] Lemma lem:sizes(e)). -/
lemma Xsp_sq_lt {N : ℕ} (hodd : Odd N) : 2 * Params.Xsp N ^ 2 < N := by
  have hs : ((Params.Xsp N : ℕ) : ℝ) ≤ Real.sqrt ((N : ℝ) / 2) := Nat.floor_le (by positivity)
  have h2 : ((Params.Xsp N : ℕ) : ℝ) ^ 2 ≤ (N : ℝ) / 2 := by
    calc ((Params.Xsp N : ℕ) : ℝ) ^ 2 ≤ (Real.sqrt ((N : ℝ) / 2)) ^ 2 :=
          pow_le_pow_left₀ (by positivity) hs 2
      _ = (N : ℝ) / 2 := Real.sq_sqrt (by positivity)
  have h3 : 2 * Params.Xsp N ^ 2 ≤ N := by
    have : (2 * Params.Xsp N ^ 2 : ℕ) ≤ (N : ℝ) := by push_cast; linarith
    exact_mod_cast this
  obtain ⟨m, hm⟩ := hodd
  generalize Params.Xsp N ^ 2 = q at h3 ⊢
  omega

set_option maxHeartbeats 1600000 in
/-- Descendable primes are primes of `J = (Y, 2Y]` coprime to `N` (`N_ℓ = 0` when `ℓ ∣ N`). -/
lemma descPrimes_spec (P : Params) (a : ℝ) (N : ℕ) (A : P.Field) :
    ∀ p ∈ P.descPrimes a N A, p.Prime ∧ P.Yj N < p ∧ p ≤ 2 * P.Yj N ∧ Nat.Coprime p N := by
  intro p hp
  unfold Params.descPrimes at hp
  rw [Finset.mem_filter] at hp
  obtain ⟨hpJ, hc⟩ := hp
  unfold Params.Jset primesIoc at hpJ
  rw [Finset.mem_filter, Finset.mem_Ioc] at hpJ
  obtain ⟨⟨h1, h2⟩, hpr⟩ := hpJ
  have hnd : ¬ p ∣ N := by
    intro hdv
    unfold Params.countN at hc
    rw [ite_eq_left hdv] at hc
    omega
  exact ⟨hpr, h1, h2, (Nat.Prime.coprime_iff_not_dvd hpr).mpr hnd⟩

set_option maxHeartbeats 1600000 in
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
  intro ε hε
  have hk := Params.adm_k hP
  obtain ⟨hν, hνk, hμν, hμ, hh, hθ, hμ1, hd⟩ := Params.adm_real hP
  set kR : ℝ := (P.k : ℝ) with hkR
  have hkR2 : (2 : ℝ) ≤ kR := by rw [hkR]; exact_mod_cast hk
  set ε₁ := ε / (4 * (2 * kR + 1)) with hε₁
  have hε₁0 : 0 < ε₁ := by positivity
  obtain ⟨N₁, hN₁⟩ := relation_supply P hP a ha hCEP hFM hLD ε₁ hε₁0
  set ε₂ := ε / (8 * kR) with hε₂
  have hε₂0 : 0 < ε₂ := by positivity
  set ν : ℝ := (P.ν : ℝ) with hνdef
  have hl2 : 0 < Real.log 2 := Real.log_pos (by norm_num)
  set Kc : ℝ := Real.log 4 + 2 * kR * Real.log 4 + 2 * kR * Real.log kR +
    2 * kR * (Real.log 5 - Real.log 3) + Real.log (4 * kR ^ 2 * cK P.k) + 6 * kR * Real.log 2
    with hKc
  set Ca : ℝ := Real.log 21 + 3 * Real.log 2 with hCa
  set Cd : ℝ := 3 * Real.log 2 + Real.log 4 + Real.log 5 - Real.log 3 + Real.log (4 * kR ^ 2)
    with hCd
  have key : ∀ᶠ N : ℕ in atTop, Params.SplitInput a N →
      ∀ Rej : Finset P.Field, Rej ⊆ P.fields N →
        prob (P.fields N) (fun A => A ∈ Rej) ≤ P.gammaD a N / 2 →
          LN (1 / 2) (-ε) N ≤
            expect (P.fields N) (fun A => if A ∈ Rej then 0 else P.descentProb a N A) := by
    filter_upwards [eventually_ge_atTop N₁, ev_Yj hP, ev_Bsm ha, ev_Lam_sq,
      ev_Lam_ge (2 * Real.log 2 / a), ev_loglog_ge (8 * Real.log 2 / (a * ε)),
      ev_loglog_le ε₂ hε₂0, ev_loglog_le 1 one_pos,
      ev_Lam_le (ν / 2 / (ε₁ + 1)) (by positivity), ev_log_ge (Ca / (2 * ν)),
      ev_log_ge (2 * Cd / ν), ev_Lam_ge (4 * Kc / ε), ev_log_ge 1, ev_Lam_ge 1] with
      N hN1 hY hB hsq hΛB hllω hll2 hll1 hΛν hlogA hlogD hΛK hlog1 hΛ1
    intro hsplit Rej _ hRej
    obtain ⟨hj, h2Y, h4Y, hY1, hYlo, hYhi⟩ := hY
    obtain ⟨hB2, hBlo, hBhi⟩ := hB
    obtain ⟨hN1', hodd, hNB, hpp⟩ := hsplit
    set t := Real.log (N : ℝ) with ht
    set Λ := Lam (1 / 2) (N : ℝ) with hΛdef
    set γ := P.gammaD a N with hγdef
    set Y := P.Yj N with hYdef
    have hΛ0 : 0 < Λ := by linarith
    -- γ ≥ L^{-ε₁}
    have hγ : Real.exp (-(ε₁ * Λ)) ≤ γ := by
      have := hN₁ N hN1 ⟨hN1', hodd, hNB, hpp⟩
      unfold LN at this; rw [show -ε₁ * Lam (1 / 2) (N : ℝ) = -(ε₁ * Λ) by ring] at this
      exact this
    have hγ0 : 0 < γ := lt_of_lt_of_le (Real.exp_pos _) hγ
    have hlγ : -(ε₁ * Λ) ≤ Real.log γ := (Real.le_log_iff_exp_le hγ0).mpr hγ
    -- sizes of `Y` and `J`
    have hνk' : ν * (4 * kR) = 1 := hνk
    have hν8 : ν ≤ 1 / 8 := by
      have : ν * (4 * kR) ≥ ν * 8 := by nlinarith
      linarith
    have hY21 : (21 : ℝ) ≤ Y := by
      have h1 : Real.log 21 ≤ Real.log Y := by
        have h := hlogA; rw [div_le_iff₀ (by positivity)] at h; linarith
      exact (Real.log_le_log_iff (by norm_num) (by linarith)).mp h1
    have hY21' : 21 ≤ Y := by exact_mod_cast hY21
    have hlY0 : 0 < Real.log Y := Real.log_pos (by linarith)
    have hJ := card_primesIoc_ge hRS Y hY21'
    have hJdef : P.Jset N = primesIoc Y (2 * Y) := rfl
    rw [← hJdef] at hJ
    set Jc : ℝ := ((P.Jset N).card : ℝ) with hJc
    have hJ0' : 0 < 3 * (Y : ℝ) / (5 * Real.log Y) := by positivity
    have hJ0 : 0 < Jc := lt_of_lt_of_le hJ0' hJ
    have hlJ : Real.log 3 + Real.log Y - Real.log 5 - Real.log (Real.log Y) ≤ Real.log Jc := by
      have h := Real.log_le_log hJ0' hJ
      rw [Real.log_div (by positivity) (by positivity), Real.log_mul (by norm_num) (by positivity),
        Real.log_mul (by norm_num) (by positivity)] at h
      linarith
    have hllY : Real.log (Real.log Y) ≤ Real.log t := by
      have h1 : 2 * ν * t ≤ 2 * (1 / 8) * t :=
        mul_le_mul_of_nonneg_right (by linarith) (by linarith)
      exact Real.log_le_log hlY0 (by linarith)
    have hYk : t - 6 * kR * Real.log 2 ≤ 2 * kR * Real.log Y := by
      have := mul_le_mul_of_nonneg_left hYlo (by linarith : (0 : ℝ) ≤ 2 * kR)
      have e : 2 * kR * (2 * ν * t - 3 * Real.log 2) = t - 6 * kR * Real.log 2 := by
        linear_combination t * hνk'
      linarith
    -- ω(N) log 2 ≤ (ε/4) Λ
    have hω0 := omega_le N (Params.Bsm a N) hN1' hB2 hNB
    have hΛB' : 2 * Real.log 2 ≤ a * Λ := by
      have h := hΛB; rw [div_le_iff₀ ha] at h; linarith
    have hω : (omega N : ℝ) * Real.log 2 ≤ ε / 4 * Λ :=
      omega_aux (by positivity) hΛ0 ha hε hl2 hsq (by linarith) (by linarith) hω0 hllω
    -- the uniform lower bound `Q`
    set Den : ℝ := 4 * kR ^ 2 * cK P.k * 2 ^ omega N * N with hDen
    have hcK : 0 < cK P.k := by
      unfold cK; have : (0 : ℝ) < Nat.choose (2 * P.k) P.k := by
        exact_mod_cast Nat.choose_pos (by omega)
      positivity
    have hN0 : (0 : ℝ) < N := by exact_mod_cast (by omega : 0 < N)
    have hDen0 : 0 < Den := by positivity
    set G : ℝ := γ * Jc / 4 / kR with hG
    have hG0 : 0 < G := by positivity
    set Q : ℝ := G ^ (2 * P.k) / Den with hQ
    have hQ0 : 0 < Q := by positivity
    -- `γ |J| / 4 ≥ 4k²`
    have hbig : 4 * kR ^ 2 ≤ γ * Jc / 4 := by
      have hpos : 0 < γ * Jc / 4 := by positivity
      rw [← Real.exp_log hpos, ← Real.exp_log (by positivity : (0 : ℝ) < 4 * kR ^ 2)]
      apply Real.exp_le_exp.mpr
      rw [Real.log_div (by positivity) (by norm_num), Real.log_mul hγ0.ne' hJ0.ne']
      have h1 : (ε₁ + 1) * Λ ≤ ν / 2 * t := by
        have h := hΛν
        rw [div_mul_eq_mul_div, le_div_iff₀ (by positivity)] at h; linarith
      have h2 : Cd ≤ ν / 2 * t := by
        have h := hlogD; rw [div_le_iff₀ hν] at h; linarith
      have h3 : 0 ≤ ν * t := by positivity
      rw [hCd] at h2
      have h4 : Real.log t ≤ Λ := by linarith
      linarith
    -- per-field bound from Theorem K
    have hXN := Xsp_sq_lt hodd
    have : NeZero N := ⟨by omega⟩
    have hQA : ∀ A ∈ P.fields N, A ∉ Rej →
        γ / 4 ≤ ((P.descPrimes a N A).card : ℝ) / Jc → Q ≤ P.descentProb a N A := by
      intro A _ _ hA
      have hDc : γ * Jc / 4 ≤ ((P.descPrimes a N A).card : ℝ) := by
        rw [le_div_iff₀ hJ0] at hA; linarith
      have hDcard : 4 * P.k ^ 2 ≤ (P.descPrimes a N A).card := by
        have : ((4 * P.k ^ 2 : ℕ) : ℝ) ≤ ((P.descPrimes a N A).card : ℝ) := by
          push_cast; linarith
        exact_mod_cast this
      have hD := descPrimes_spec P a N A
      have hK := complete_descent_prob N P.k (Params.Xsp N) Y hodd (by omega) hXN h2Y h4Y
        (P.descPrimes a N A) hD hDcard
      refine le_trans ?_ hK
      apply div_le_div_of_nonneg_right _ hDen0.le
      apply pow_le_pow_left₀ hG0.le
      exact div_le_div_of_nonneg_right hDc (by positivity)
    have hsel := select_bound (P.fields N) Rej
      (fun A => ((P.descPrimes a N A).card : ℝ) / Jc) (P.descentProb a N) γ Q hγ0
      (fun A _ => ⟨by positivity, div_le_one_of_le₀ (by
        rw [hJc]; unfold Params.descPrimes; exact_mod_cast Finset.card_filter_le _ _) hJ0.le⟩)
      (expect_descPrimes P a N) hRej (fun A _ => prob_nonneg' _ _) hQA hQ0.le
    refine le_trans ?_ hsel
    -- `L^{-ε} ≤ (γ/4) Q`
    have hpos : 0 < γ / 4 * Q := by positivity
    unfold LN
    rw [show -ε * Lam (1 / 2) (N : ℝ) = -(ε * Λ) by ring, ← Real.exp_log hpos]
    apply Real.exp_le_exp.mpr
    have hlF : Real.log (γ / 4 * Q) = Real.log γ - Real.log 4 +
        2 * kR * (Real.log γ + Real.log Jc - Real.log 4 - Real.log kR) -
        (Real.log (4 * kR ^ 2 * cK P.k) + (omega N : ℝ) * Real.log 2 + t) := by
      rw [Real.log_mul (by positivity) hQ0.ne', Real.log_div hγ0.ne' (by norm_num),
        hQ, Real.log_div (by positivity) hDen0.ne', Real.log_pow, hG,
        Real.log_div (by positivity) (by positivity), Real.log_div (by positivity) (by norm_num),
        Real.log_mul hγ0.ne' hJ0.ne', hDen, Real.log_mul (by positivity) hN0.ne',
        Real.log_mul (by positivity) (by positivity), Real.log_pow]
      push_cast; ring
    have hK : Kc ≤ ε / 4 * Λ := by
      have h := hΛK; rw [div_le_iff₀ hε] at h; linarith
    exact descent_arith (by linarith) hlF hlγ (by rw [hε₁]; field_simp) hlJ hYk
      (le_trans hllY hll2) (by rw [hε₂]; field_simp; ring) hω hKc hK
  obtain ⟨N₀, hN₀⟩ := eventually_atTop.mp key
  exact ⟨N₀, fun N hN hs => hN₀ N hN hs⟩

end FactoringA
