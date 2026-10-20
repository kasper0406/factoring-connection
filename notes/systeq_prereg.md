# systeq pre-registration (agent `systeq`, phase 16, 2026-10-20)

Committed before any run of the experiments below.

## Disclosed pre-registration activity

The following was run before this file was committed:
- A code smoke test of the positive control for q = 4, 5, 7. It gave G_FF = 9.1, 38.2 and 1265 and verified every relation mod I.
- A bit-length and timing check of the E2 parameters: no smoothness counts were taken.
- A tiny Descartes check (n = 7, H = 2): 0 nondegenerate split members.

The q ≤ 7 positive-control rows are therefore *not* predictions. Only q ∈ {8, 9, 11, 13, 16} are predicted.

## Question

The small-characteristic engine has a structural feature: the systematic equation makes one side of a relation smooth BY CONSTRUCTION. Does the closest integer analogue give any such smoothness gain? The gain is measured against the size-bound prediction, i.e. random elements of the same size.

**Gain** G = (observed relation yield) / (yield predicted if each side were a random element of its size).

## E1: positive control, F_{q²}[X]/(I) with I | h₁X^q − h₀ (BGJT/KW model)

- **Candidates:** all of P(V)(F_{q²}), V = span(X^{q+1}, X^q, X, 1), enumerated exactly.
- **Relation:** f splits over F_{q²}, and r_f = A X h₀ + B h₀ + C X h₁ + D h₁ (degree ≤ 3, h₁f ≡ r_f mod I) also splits.
- **Size prediction:** P_rand(q+1)·P_rand(3). Here P_rand(D) = C(Q+D, D)/((Q^{D+1}−1)/(Q−1)), the exact fraction of split binary forms of degree D.
- **Model M_KW:** frac_split(V) ≈ 8/q³ and P(r_f splits | f split) ≈ 1/6. This gives the predicted G_FF* = (8/q³)/6 / (P_rand(q+1)P_rand(3)):

  | q | G_FF* |
  |---|---|
  | 8 | 2.7e3 |
  | 9 | 2.0e4 |
  | 11 | 1.5e6 |
  | 13 | 1.7e8 |
  | 16 | 3.8e11 |

- **Predictions:**
  - (P1a) |log₁₀ G_FF − log₁₀ G_FF*| ≤ 1 for each q ∈ {8, 9, 11, 13, 16}.
  - (P1b) G_FF is increasing in q.
  - (P1c) #relations ≥ Q + 1 (the size of the linear factor base) for every q.
  - (P1d) The congruence h₁f ≡ r_f mod I holds for every relation (asserted in code).

## E2: integer analogue, the "F_{1^n}" V-family on Cunningham-type N = Φ_n(m)

- **Setup:**
  - m^n ≡ 1 (mod N), so μ_n = ⟨m⟩ is an explicit "constant field" in Z/N, and X^n ≡ 1 plays the role of X^q ≡ h₀/h₁.
  - For f = A x^{n+1} + B x^n + C x + D with |A|, …, |D| ≤ H = 32 and A ≠ 0, f(m) ≡ r_f := (A+C)m + (B+D) (mod N). This is checked on every candidate.
- **Relation:** f(m) and r_f are both B-smooth, with B = 2^16.
- **Instances:** (n, m) ∈ {(5, 520), (7, 130), (11, 19), (13, 11)}, giving N of 37–43 bits. There are 200 000 seeded candidates each.
- **Size prediction:** Σ base(bitlen f(m))·base(bitlen r_f). Here base(ℓ) is the exact empirical B-smooth rate of 200 000 seeded uniform integers in [2^{ℓ−1}, 2^ℓ).
- **Predictions:**
  - (P2a) G ∈ [0.5, 2] for every instance.
  - (P2b) **Kill signal for the barrier:** G ≥ 4 on any instance would mean that structured integer smoothness exists, and would trigger investigation.
  - (P2c) The f(m) side has log-size ≥ log N. This is the Dixon/L[1/2] regime, which the barrier predicts.

## E3: Descartes count

- **Setup:** exhaustive over A x^{n+1} + B x^n + C x + D ∈ Z[x], coefficients in [−4, 4], A ≠ 0, n = 2..10. Count the members that split completely over Q, excluding the degenerate x^n(Ax + B) (C = D = 0).
- **Prediction (P3, PROVED by Descartes' rule, so this is an exact check):** 0 nondegenerate split members for every n ≥ 7.
- For n ≤ 6 the counts are reported without prediction. n = 2 is all cubics, so the count is positive.

## E4: cyclotomic calibration (no N)

- **Setup:** m^n − 1 = ∏_{d|n} Φ_d(m), with n = 6, m ∈ [2^20, 2^20 + 30000), and n = 12, m ∈ [4096, 16384); B = 2^16.
- **Per-piece prediction:** Σ_m ∏_d base(bitlen Φ_d(m)).
- **Prediction (P4):** R_piece = observed/predicted ∈ [0.1, 1.5] for both instances. Cyclotomic splitting gives at most per-piece-size smoothness, and the residue-class restrictions (Zsigmondy) push it down. It never gives smoothness by construction.

## Command

```
.venv/bin/python -m src.systeq.run_experiment all
```

Config: `experiments/configs/systeq.json` (seed 0). Raw: `experiments/raw/systeq_E{1,2,3,4}.jsonl`. Memory is under 300 MB; it runs as one job.
