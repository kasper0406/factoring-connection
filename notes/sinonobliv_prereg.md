# sinonobliv pre-registration (written and committed before any experiment run)

Agent `sinonobliv`, phase 22, 2026-10-22. Code: `src/sinonobliv/` (`core.py`, `msample.py`, `theory.py`,
`run_exp.py`). Commands: `.venv/bin/python -m src.sinonobliv.run_exp E1|E2A|E2B|all`. Raw output:
`experiments/raw/sinonobliv_E{1,2A,2B}.jsonl`. Memory: one job at a time, expected RSS < 1 GB.

Design calculations made before registration (formula evaluations only, no sampling): the expected tuple counts
per call of E2B from `msample.predicted_tuples` (e.g. X = 2^36, h = 2^10, B = 32: m = 1 Johnson 0.464,
m = 2 CH 0.0459, m = 3 CH 0.0133; "bare" d > 2(h+1): m = 2 43.4, m = 3 13.4), used to set the number of calls.

## E1 (route 2: identity families with forced factorisations)

Families (N-oblivious): cp2 = t(t+1); cp3 = t(t+1)(t+2); aur = 4t⁴+1 = (2t²+2t+1)(2t²−2t+1); sq1 = t²+1
(irreducible control). Smoothness bounds y ∈ {2^8, 2^10}.
- cp2: bits b ∈ {40, 48, 56}; t ∈ [t0, t0 + T), t0 = least t with t(t+1) ≥ 2^b, T = min(t0/16, 2^21).
- cp3: b ∈ {48, 56, 60}; t ∈ [t0, 2t0).
- aur: b ∈ {40, 46, 52}; t ∈ [t0, 4t0).
- sq1: b ∈ {40, 48, 56}; t ∈ [t0, t0 + T) as cp2.
Every member's smoothness is decided exactly (all pieces y-smooth). The t-range is split into 16 equal bins.
Independence prediction: pred = Σ_bins n_b·∏_j p̂_j(b), with p̂_j(b) the exact fraction of y-smooth piece-j values in
bin b. R := (#members with all pieces smooth)/pred.
Control: exact y-smooth density of random integers near the member sizes (DensityCurve: windows of 2^21 integers at
half-bit steps, log-linear interpolation); gain := #smooth members / Σ control(member) (descriptive).

- **P1.** For cp2, cp3, aur, in every cell with pred ≥ 30: R ∈ [0.67, 1.5] (forced pieces are smooth
  independently at the level of their marginal densities: heuristic H_piece).
- **P2 (hit law).** cp2 and cp3, b ∈ {40, 48}, 20 random odd N of b bits (seeded), K = 64 multipliers, windows
  I_k = [kN, kN + h], h from siconstruct's toy L-path (`theory.toy_path`). Exact member counts (`count_members_in_window`)
  against the continuous prediction Σ_k (h+1)/F′(t_k), t_k = F^{-1}(kN). Pooled measured/predicted ∈ [0.95, 1.05]
  in every (family, b) cell whose prediction is ≥ 200; cells with smaller predictions are reported only.
- **D1 (descriptive / MODEL).** Per-window expected smooth members of each family (hits per window × measured joint
  rate) against the interval sieve's h·δ; model extrapolation to 512–4096 bits with `theory.consecutive_product_log_hits`.
  No band: toy h is a large power of N (siconstruct §4.2), so toy economics favour the families.

## E2 (route 3: multi-sample joint certificates)

**E2A (Lemma IH, exhaustive).** b ∈ {32, 40}, 10 seeded random odd N each, multipliers k = 1..8 (K = 8), h = 2^9.
For all pairs k_i < k_j and all offsets r_i, r_j ∈ [0, h], g = gcd(k_iN + r_i, k_jN + r_j).
- **P3.** Zero violations: every pair with g > K·h has k_j·r_i = k_i·r_j. (A theorem; the run is a sanity check.)
- D2: report the counts with g > Kh and with h < g ≤ Kh, and the proportional share of the latter.

**E2B (independent windows, Monte Carlo; seeds fixed).** Configs (X, h, B, y): C1 = (2^36, 2^10, 32, 2^8),
C2 = (2^44, 2^12, 41, 2^10); P = lcm(1..B). Window positions independent uniform in [X, 2X). A certified tuple is
(d, z_1..z_m) with d | P, ln d ≥ σ, d | z_i ∈ [x_i, x_i + h]. Cells: m = 1 with σ = Johnson √(tπ); m = 2, 3 with
σ = Cohn–Heninger t^{m/(m+1)}π^{1/(m+1)} (t = ln(h+1), π = ln P); and "bare" m = 2, 3 with σ = ln(2(h+1)).
Calls: C1: 20000 (m = 1, 2, 3), 2000 (bare); C2: 40000 (m = 1, 2, 3), 2000 (bare).
- **P4.** In every cell with predicted tuples ≥ 50: measured/predicted tuple count ∈ [0.85, 1.15], predicted
  = calls·Σ_{d | P, ln d ≥ σ} ((h+1)/d)^m.
- **P5.** In every cell with ≥ 30 tuples with z₁ smooth: (tuples with z₁ smooth)/(model Σ_d ((h+1)/d)^m·δ(1.5X/d))
  ∈ [0.6, 1.6] (cofactor smooth like a random integer of its size).
- **P6.** Per window used, the smooth yield of the CH cells (m = 2, 3; distinct smooth z₁ per call, divided by m)
  is below that of the m = 1 Johnson cell, in both configs.
- D3 (descriptive): the divisor-trick baseline (d = 2^⌊log₂(h+1)⌋, 20000 calls) smooth yield per call, and the bare
  cells' yield per certified tuple.

## Deviations
Any change after the first run will be disclosed in the results note as a deviation, with the reason.

## Addendum A1 (2026-10-22, written after E1, E2A, E2B had run and before E3 was run)

E1, E2A and E2B have been run (results in the results note; P1 failed in 4 cells, P4/P5 in one cell). This addendum
registers one further experiment, motivated by codex's blind round (sinonobliv_1) listing *uncertified decoder outputs*
as the least constrained open item. Code: `src/sinonobliv/uncert.py`; command `run_exp E3`; raw
`experiments/raw/sinonobliv_E3.jsonl`.

**E3 (uncertified roots of a Howgrave-Graham divisor lattice).** f(r) = x + r, modulus P = lcm(1..B); basis
P^{ℓ−j}f^j (j ≤ ℓ), r^i f^ℓ (1 ≤ i ≤ s), scaled by V; python-flint LLL; integer roots in [0, V] of the first three
reduced polynomials. Certified threshold T_cert = (√n·‖b₁‖)^{1/ℓ} (Howgrave-Graham's condition for the first vector).
Configs: U1 = (x uniform in [2^40, 2^41), V = 2^12, B = 40, ℓ = 3, s = 3, y = 2^10); U2 = (x in [2^48, 2^49),
V = 2^10, B = 60, ℓ = 4, s = 4, y = 2^12); 5000 calls each, seeds 6022 / 6023.
- **P7.** Mean number of uncertified roots (gcd(x + r, P) < T_cert) per call ≤ 0.5 in both configs.
- **P8.** If a config has ≥ 30 uncertified roots: their y-smooth fraction divided by the model mean δ(z/gcd(z, P))
  (DensityCurve at the cofactor size; 1 if the cofactor ≤ y) lies in [0.33, 3] (no smoothness bias beyond the
  known divisor).
- **P9.** Zero violations: every r ∈ [0, V] with gcd(x + r, P) ≥ T_cert is a root of the first reduced polynomial.
- D4 (descriptive): certified roots per call and their smooth fraction against the same model.
