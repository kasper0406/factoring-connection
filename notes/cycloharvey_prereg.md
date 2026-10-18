# cycloharvey: pre-registration (frozen before the search)

Agent `cycloharvey`, 2026-10-18. This file is committed before any run of the
constructive covering search (`src/cycloharvey/search.py`). The calibration gate
(`src/cycloharvey/run_calibration.py`) and two smoke tests of the covering code at
40 bits (supports (−1,1), (1), (−1,1,2), (−2,−1,1,2) with W ∈ {16, 64, 256}) ran before
this file was written. Those smoke tests were used only to debug the code and to set
time budgets.

## 0. Question

Can the Frobenius handles σ_r of R = (Z/N)[X]/(Φ_m), with r ≡ p (mod m) guessed,
support a Lehman/Hittmeir/Harvey-style one-sided collision search with model
exponent < 1/5? The candidate mechanism is homogeneous higher-degree "Taylor
cancellation" in exponents F(p,q).

## 1. Model (frozen)

Identity model. A target τ is a pair (g_τ, window {U_τ, …, U_τ+W_τ−1}).
- **Target function.** g_τ(p) = F_τ(p, N/p), where F_τ ∈ Z[P,Q] has total degree ≤ 3.
  Equivalently, g_τ is a Laurent polynomial Σ_{k∈I} a_k p^k with I ⊆ {−3,…,3}∖{0},
  a_k ∈ Z for k > 0, and a_k ∈ N^{|k|}Z for k < 0. The constant term is irrelevant.
- **Residue guess.** Optionally, p ≡ r (mod M) is guessed, with φ(M) runs. The target
  may then be any element of the integer-valued lattice L_{r,M} = {(F − c)/M^e : F(p,q) ≡ c
  (mod M^e) on the class}.
- **Correctness.** The correctness proof may use only the polynomial identity: every
  real p in the balanced range [√(N/2), √N] (in the class, if guessed) must have g_τ(p)
  in the window of some target. This is how Lehman, Hittmeir, Harvey and HH argue.
- **Cost.** cost = m + Σ_τ ⌈W_τ/m⌉, with a shared baby list of size m and one handle per
  target. This is Harvey's full sharing; HH is the special case of one giant step per
  branch. For the cyclotomic version every ring operation is charged D = φ(m) base
  operations, which is a lower bound since multiplication costs ≥ D. The cost is also
  multiplied by φ(m) guesses. Cost of ord-finding and of collision detection: polylog,
  ignored (as in all calibrated algorithms).
- **One-sidedness.** No cost is charged. A two-sided collision (≡ mod N) is caught by
  sort-and-match (Harvey Alg. 4.2 Step 3 / HH Step 12) and yields the integer u =
  g(p). p is then recovered by solving g(p) = u, a bounded-degree polynomial with
  rational roots. A one-sided collision is caught by the product-tree/Bluestein gcd.
  Both are therefore detected at the same polylog cost. (Codex r7: symmetric F at
  p ≡ q mod m is two-sided. This does not hurt the model.)

## 2. Family / search space (frozen)

- All 63 supports I ⊆ {−3,−2,−1,1,2,3} (every F of total degree ≤ 3). Coefficient
  bound C = N^η, with η ∈ [0, 1/2] optimised. A constant C is the case η = 0.
- Residue guesses: analytic M = N^μ with μ ∈ [0, 1/2]. Constructive M ∈ {1, 5, 12, 60}
  (conductors m ≤ m_max = 60).
- Constructive supports: (1), (−1,1), (1,2), (−1,1,2), (−2,−1,1), (−2,−1,1,2),
  (−3,−1,1,3), (−1,1,2,3), and the full support (−3,−2,−1,1,2,3). The full-support
  run searches all degree-≤3 F simultaneously.

## 3. Objective (frozen)

(A) **Analytic (lattice-volume heuristic).** For each support,
e(I, μ) = S/(K(K+3)) + μ·(1 − 2(K+κ)/(K(K+3))) [+ μ for cyclotomic D-arithmetic],
where K = |I|, S = Σ|i|, and κ = log_M of the index of L_{r,M} over Λ_I. κ is measured
exactly by HNF (`lattice.class_lattice_index`). Objective: min over I and μ.

(B) **Constructive (exact).** `covering.sweep`: a greedy left-to-right cover of the
balanced range. Candidates come from LLL. Every accepted target is verified exactly
(rational arithmetic, rigorous arb root balls for critical points). Uniform window
W = 2^j for integer j ∈ [⌈0.12n⌉, ⌈0.32n⌉]. The cost for each n is the minimum over
W of min_m (m + n_T(W)·⌈W/m⌉), with N = 2^n + 2^{n/2} + 1. Sizes n ∈ {40, 48, 56, 64} use
full sweeps. Sizes n ∈ {72, 80} are sampled: 4 sub-intervals at relative positions 1/8,
3/8, 5/8 and 7/8 of the range, each of relative length 2^{−6}, with n_T extrapolated
linearly. The slope s_I is fitted as log2 cost on n over n ∈ {56, 64, 72, 80}.
Guessed-class runs: n ∈ {48, 56}, one class r (the smallest unit), with cost multiplied
by φ(M).

## 4. Decision rules (frozen)

- **Calibration gate (hard).** The model must reproduce 1/4, 1/3, 2/9 and 1/5 exactly in
  the rational LP (`calibration.exponent_lp`). Its exact-count fits must land within
  ±0.03 of each target, allowing for the polylog drift seen in the fits. If either
  fails, fix the model before the search.
- **Promising.** Some I or μ gives analytic e < 1/5 − 0.01. Alternatively, some
  constructive slope s_I < 0.19, or the full-support slope satisfies
  s_full < s_(−1,1) − 0.02. If anything is promising, run step 4 (exact factoring
  experiments, 32–56 bits, against a Harvey-style baseline).
- **Refuted for this model and space.** No family in §2 reaches analytic exponent
  < 1/5 − 0.01 and no constructive signal is promising. The cyclotomic exponent is then
  the F_p^* exponent plus μ (D-arithmetic), so it is ≥ the F_p^* exponent.
- A constructive slope that differs from the analytic exponent by more than 0.03 is
  reported as a discrepancy, not silently resolved.

## 5. Secondary analysis (orchestrator yardstick, added before the search)

The yardstick is per-branch steps S ~ N^{1/2}/T^δ, giving T = N^{1/(2(1+δ))}. Harvey and
HH have δ = 3/2. For each family we report the *effective* δ_eff = 1/(2e) − 1 after
coefficient inflation (a = S/(2K) > 1/2 for higher supports), φ(m) guesses,
D-dimensional arithmetic and one-sidedness. Without guesses the heuristic gives
δ_eff = K(K+3)/(2S) − 1. The user's hope, δ = 5/2 (second order cancelled with no
coefficient blow-up), corresponds to K = 4 with S = 4, which is impossible because S ≥
1+1+2+2 = 6 for four distinct nonzero exponents. The heuristic δ_eff for K=4 is 4/3.

## 6. Analytic results derived before the search (to be checked, not assumed)

- **D1 (domination).** Every cyclotomic target (handle T_{r,F}, any F, any division by
  M^e) is also a target of the F_p^* model with the same guess r, with the same window.
  The F_p^* handle is a^{F(1,N)}-type with a = b^{M^e}. Hence e_cyc ≥ e_{F_p^*} + log_N D.
- **LB (lower bound).** In the model of §1, for any fixed degree bound d and any guess
  modulus M, cost ≥ c_d N^{1/5} M^{−o(1)}. The proof uses a sublevel-set bound, a
  Wronskian, a lattice count and Hölder. The full proof is in notes/cycloharvey_results.md
  and will be audited before it is claimed as proved.
- **κ lemma.** κ = K(K+1)/2 for primes M > 2·max|i|. Small primes give a bounded excess.
  The HNF table in experiments/analysis/cycloharvey_kappa.json agrees for M ∈ {5, 7, 11,
  13}. With guessing exactly neutral, HH's gain is (φ(M)/M)^{3/5}, the (log log N)^{3/5}
  of their theorem.

## 7. Step 4 (pre-registered regardless of outcome, as calibration only)

- Operational Harvey-style search (randomised base, Alg. 4.1 collision detection) on
  balanced semiprimes of 32–56 bits, with counted list sizes and ring ops. The model's
  count is predicted from N alone.
- Cyclotomic handle check (m ∈ {5, 8, 12}): verify T_{r,F} ≡ α^{F(p,q)} mod p (harness).
  Measure one-sided rates for asymmetric F and symmetric F, with p ≡ q and p ≢ q
  (mod m).
