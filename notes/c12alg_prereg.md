# c12alg pre-registration (E1, E2)

Agent `c12alg`, phase 27, 2026-10-24. Written and committed before any registered run. Code frozen at this commit:
`src/c12alg/dw.py`, `src/c12alg/e1.py`, `src/c12alg/dw_pipeline.py`, `src/c12alg/harness_e2.py`,
config `experiments/configs/c12alg_e2.json`.

**Smoke runs before registration (disclosed).**
- `e1.run_e1b` with 10⁵ samples and one field, and `e1.run_e1a` on 2 pairs at (d, V, H) = (3, 2, 40). This was a
  timing/correctness check. Its outputs (scratchpad only, not committed) were looked at: per-pair J/J_gcd ≈ 1.10 for
  one nearly parallel pair, with archimedean proxy ≈ 0.11; TV = 0 for a pair with b₁b₂ = 1. After that I added the
  condition b₁b₂ ≥ 2 to the E1a pair sampler, so that every registered TV case has a non-trivial representative set.
- `dw_pipeline.factor_dw` on thmA's calibration seeds 900 (semi, 32 bits) and 901 (tri, 40 bits), with
  W_mult = 64. Both were factored correctly. Content rejections were 8–11% of sampled fields. E2 uses fresh seeds
  (≥ 7000).

These are toy correctness checks of exact statements and of the algebra. None says anything about asymptotic
precision.

## E1: two-point law under design DW

**E1a (exact TV of Lemma TP).** Configurations (d, V, H) ∈ {(3,2,40), (3,2,80), (3,2,120), (3,3,80), (3,3,120),
(4,2,40)}, with 12 random pairs each. Each pair has primitive non-proportional points, b_i ≥ 1, max-norm exactly V and
b₁b₂ ≥ 2 (seed `c12alg|E1a|20261024`). The run computes exactly:
- the law of Φ(c) = (R_c(P₁), R_c(P₂)) for c uniform on [−H, H]^d;
- the law of Φ(c) + z, with z uniform on R = Γ ∩ M[0,1)² (columns d−2, d−1). By Lemma TP this second law is the
  explicit measure π ∗ Unif(Γ ∩ Π);
- the exact TV between the two;
- max over z ∈ R of ‖x(z)‖_∞ for the short preimages.

Predictions:
- **P1a.** TV ≤ d(1 + 4V²)/(2H + 1) in every case (72/72; the bound is vacuous where it exceeds 1).
- **P1b.** ‖x(z)‖_∞ ≤ 1 + 4V² for every z in every case.
- **P1c (scaling, descriptive with a band).** In the non-vacuous cases, the median over pairs of TV·(2H + 1)/(dV²) lies
  in [0.05, 4] in each configuration, and the median TV at (3, 2, 120) is smaller than at (3, 2, 40).

**E1b (joint smoothness under DW, Monte Carlo).** d = 3, M = 8, points P = (a, b) with |a| ≤ 3, 1 ≤ b ≤ 3,
gcd(a, b) = 1, filtered by 8|G(P)| ≥ M·|P|_∞ and gcd(G(P), N) = 1. W_m = 2¹⁴, H = 2¹⁹, y = 512 (exact
trial-division smoothness). Six fields: toy N ∈ (2¹², 2¹³), and g primitive in [4, 8]². Each field uses 2·10⁶
samples of (w, c), w uniform in (W_m, 2W_m] and c uniform. For every pair: the marginals M_i, the joint J, and
J_gcd := Σ_{δ|Δ} (φ(Δ/δ)/Δ)·s̄₁(δ)s̄₂(δ). Here s̄_i(δ) is the empirical smooth fraction of F(P_i) in the class
gcd(F(P_i), Δ) = δ, so J_gcd is the lattice principal part with conditional independence given gcd class. The
archimedean proxy is arch := E[Z₁Z₂]/(E Z₁·E Z₂) − 1 with Z_i = ρ_Dickman(log|F_i|/log y), and J_full := J_gcd(1 + arch).

Predictions (statistics over all pairs of all fields; CIs by field bootstrap, 2000 resamples):
- **P1d.** Pooled Σ J/Σ J_gcd ∈ [0.95, 1.08].
- **P1e.** Pooled Σ J/Σ J_full ∈ [0.97, 1.03], CI ⊂ [0.95, 1.05].
- **P1f.** Regressing (J/(M₁M₂) − 1) on (J_full/(M₁M₂) − 1) through the origin gives slope ∈ [0.8, 1.2].
- **P1g.** Residual excess: the mean over pairs of (J − J_full)/(M₁M₂) is within ±0.03.

Reading fixed in advance:
- E1a is a correctness check of an exact lemma.
- E1b tests that, at toy scale, the DW pair law has no correlation beyond (lattice gcd classes) × (archimedean sizes),
  i.e. no unexpected first-order term. It cannot test precision D^{−α}.

## E2: the toy pipeline with DW fields (composite w)

The thmA toy pipeline, with SAMPLE_FIELD replaced by `sample_field_dw` (repairs of notes/c12alg_results.md §4):
- w uniform in (W, 2W] with W = 64·M, composite allowed;
- F = w·F̂_{g,N} + G·R_c;
- reject degree drop and content > 1;
- the IR certificate is unchanged (primes r | w cannot certify).

All other thmA settings are unchanged (`experiments/configs/c12alg_e2.json`). Instances: block E2 has 24 instances
(semi 30/34/38/42 bits × 3, tri 40/44 × 2, unbal 40 × 2, p2q 36 × 2, small_x_semi 40 × 2, prime 36 × 2). Block E2b
(diagnostic, continue_on_accidental) has 8 instances (semi 44 × 4, tri 46 × 4). Seeds start at 7000.

Predictions:
- **P2a.** Correct complete factorisation in 24/24 (E2) and 8/8 (E2b).
- **P2b.** 0 row-identity violations, 0 congruence violations, 0 CU violations (F̂_N(P) a unit mod |G(P)|), 0
  forced-residue violations (F(P) ≡ w·F̂_N(P) mod |G(P)|), and 0 IR_DW claim violations.
- **P2c.** Every tested kernel vector that passes the character test is a square in K (as thmA).
- **P2d.** The content-rejection fraction among sampled fields that pass the gcd checks and the degree check is within
  3 binomial standard errors of P_cont := Pr_{w,c}[gcd(w, c₀, …, c_{d−1}) > 1]. P_cont is computed exactly by Möbius
  from the toy law (w uniform in (W, 2W] for each instance's W, c uniform in [−4, 4]³). Lemma CT predicts that
  content > 1 occurs exactly when some p | w divides every c_i.
- **P2e.** DD fruitfulness in E2b (diagnostic mode): for ω = 2 splits, the fruitful fraction is within |z| ≤ 3 of 1/2;
  for ω = 3, within |z| ≤ 3 of 3/4.

Reading fixed in advance: E2 checks only the algebra of design DW (the field, the row identity, the square root, the
congruence) with composite w. It says nothing about relation supply.
