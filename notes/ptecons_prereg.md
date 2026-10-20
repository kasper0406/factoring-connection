# ptecons pre-registration (2026-10-20)

Agent `ptecons`, phase 16. This file is committed before E1 and E2 are run.
- Code: `src/ptecons/{polys,lemmas,acs,run_experiment}.py`.
- Config: `experiments/configs/ptecons.json`.
- Raw output: `experiments/raw/ptecons_E{1,2}.jsonl`.

**Disclosed pilot.** One pilot ran before this file was written: `run_experiment pilot`, at 48 and 64 bits, with seeds 900 and 901. It was used only to choose window and cap parameters, and to find and fix two bugs: the sextic window had the wrong sign range, and the quartic inner loop was slow. The fast quartic loop is checked against brute force in `tests/test_ptecons.py`. The pilot's seeds are excluded below.

Pilot readings, which inform the predictions:
- 48-bit COS median λ = 0.62;
- quart_g 0.63;
- quart_s 0.68;
- the two-digit sextic gave 24 members (48 bits) and 1 member (64 bits) against a cap of 20000.

## Hypotheses

**H1 (AVG).** Fix a structured family of integers u, built from split polynomials at fixed evaluation points and independent of N. Over seeded N, |u mod± N| ≤ V occurs at the uniform rate (2V+1)/N. Structure gives no excess congruence yield unless N enters the construction.

**H2 (alignment capacity).** Families aligned to N on their top coefficients — degree-4 generic, degree-4 symmetric and degree-6 symmetric clusters — give residuals that are typically no smaller than COS's (degree 2) at an equal pool budget. Matching two digits with three square-roots, as the sextic does, is generically infeasible.

## E1 (AVG check)

**Setup.**
- 200 balanced 24-bit semiprimes (`balanced_semiprime(24, s)` for s = 0..199), with V = 1000.
- Families (config E1):
  - generic sextic multisets in [−6, 6] at m = 40;
  - symmetric J = 3 clusters;
  - members of nontrivial M = 5 PTE classes (equal e₁..e₃);
  - Dickson integer-split members D_n(X) − D_n(t), with n ∈ {3, 4, 6}, at m ∈ [40, 2000);
  - double-quadratic composition clusters (a² + b² = 2c).
- Only pairs with |u| > N are counted.
- G = observed / expected, where expected = Σ (2V+1)/N.

**Prediction P1.** G ∈ [0.8, 1.25] for every family whose expected count is ≥ 200.

## E2 (aligned-cluster sieve variants)

**Setup.**
- Balanced semiprimes: 64 bits with seeds 1–5, and 80 bits with seeds 1–3.
- Pool cap: 20000 members per (variant, N).
- B = 2^16 for the smoothness test of v. Each v is paired with a control: a uniform random integer of the same bit length.
- λ = log₂|v| / log₂ N.

**Variants (config E2):**
- cos (A = 100);
- quart_g;
- quart_s;
- sext_s (two digits);
- sext_s_wide (s₂ window effectively open);
- quart_rand (unaligned control).

**Predictions:**
- **P2a.** cos: median λ = 0.5 + log₂(100)/n ± 0.04.
- **P2b.** At every N, each of quart_g, quart_s and sext_s_wide has median λ ≥ (cos median λ) − 0.02.
- **P2c.** λ at rank 1000 (λ₁₀₀₀) ≥ 0.45 for every variant with ≥ 1000 members.
- **P2d.** smooth_ratio ∈ [0.6, 1.6] for every (variant, N) with ≥ 30 control-smooth hits.
- **P2e.** sext_s has < 200 members at every N.
- **P2f.** quart_rand has median λ ≥ 0.95, i.e. unaligned values are uniform mod N.

**Kill signal for H2.** At two or more N, some aligned variant has median λ below cos − 0.05. That would be a structured higher-degree family beating COS, and it would trigger a scaling study.

Both experiments are toy-scale checks of heuristics. They are not factoring results.
