# hypD pre-registration (frozen before any prediction-testing run)

Agent `hypD`, phase 7, 2026-10-19. Code: `src/hypD/{split,specialq,rank,harness,analyze}.py`, tests
`tests/test_hypD.py`, configs `experiments/configs/hypD_e{1,2,3a,3b}.json`. Raw output goes to
`experiments/raw/hypD_e*.jsonl`. All arithmetic is exact. These are toy correctness and plausibility
checks of the model; **no complexity claim** is attached to any of them.

Calibration disclosure: before this file was written I ran yield and timing checks for E2 and E3 on
instance seeds 1000–1002, which are disjoint from the registered seeds 0–4. I also ran the unit tests,
including one exhaustive check of the Theorem S bound at N = 1009·1013. No registered seed was run.

Commands:
```
.venv/bin/python -m src.hypD.harness e1 experiments/configs/hypD_e1.json experiments/raw/hypD_e1.jsonl 8
.venv/bin/python -m src.hypD.harness e2 experiments/configs/hypD_e2.json experiments/raw/hypD_e2.jsonl 8
.venv/bin/python -m src.hypD.harness e3 experiments/configs/hypD_e3a.json experiments/raw/hypD_e3a.jsonl 8
.venv/bin/python -m src.hypD.harness e3 experiments/configs/hypD_e3b.json experiments/raw/hypD_e3b.jsonl 8
.venv/bin/python -m src.hypD.analyze e1 experiments/raw/hypD_e1.jsonl
.venv/bin/python -m src.hypD.analyze e2 experiments/raw/hypD_e2.jsonl
.venv/bin/python -m src.hypD.analyze e3 experiments/raw/hypD_e3a.jsonl experiments/raw/hypD_e3b.jsonl
```

## E1 — split counting (Theorem S)
- Setting:
  - Y = ⌊√(N/2)⌋ and J = (2^j, 2^{j+1}] with j = ⌊(log₂Y)/k⌋.
  - G is a random subset of the primes in J, each kept independently with probability g.
  - S is the set of products of exactly k distinct primes of G lying in [Y/2, Y].
  - R = r² mod N with r uniform. E_S is the event that R ≡ u/v with u, v ∈ S; it is decided by
    `split_in_S`.
  - Grid: classes {blum, tri} × bits {40, 48} × seeds 0–4 × k ∈ {2, 3} × g ∈ {1.0, 0.3}, 4000 r-trials
    each.
- **P1 (theorem check).** No record has an empirical rate below the Theorem S bound
  |S|⁴/(2·2^ω·N·Y²·(1+log Y)) minus 3 standard errors.
- **P2 (model).** Under the coprime-pair model, Pr[E_S] ≈ |S|²/φ(N): characters cancel and the 2^ω
  cancels.
  - The pooled ratio of observed to model hits lies in [0.8, 1.25] for each (class, k).
  - At least 95% of records have |z| ≤ 3.

## E2 — special-q counts in toy LV-style fields
- Setting:
  - Field A = (m, r0, r1). The cubic is the monic base-m polynomial plus (x − m)(r1x + r0), with r_i
    uniform in [−4, 4] and m uniform among values with m³ ≤ N < 2m³.
  - 40 primes ℓ are drawn from J = (2^14, 2^15] once per N. 60 fields per N.
  - B = 8192. Each box has |i| ≤ 32 and 0 ≤ j ≤ 32 in Gauss-reduced special-q lattice coordinates.
  - N_ℓ(A) counts primitive box pairs with (a − bm)/ℓ and F(a, b) both B-smooth.
  - Grid: classes {blum, tri} × bits {64, 72} × seeds 0–3.
- Measured per N over all (A, ℓ) cells:
  - μ = E N_ℓ and γ = Pr[N_ℓ > 0];
  - PC = E[N(N−1)]/μ²;
  - PC_alg = ΣN(N−1)/(Σ n_rat(n_rat−1)·p̂²), where n_rat is the number of rational-smooth box points and
    p̂ = ΣN/Σn_rat. PC_alg is the pair-correlation ratio of the algebraic side given the rational side.
- **P3.** PC_alg ∈ [0.67, 2.0] for every N.
- **P4.** |γ − mean(1 − exp(−n_rat·p̂_A))| ≤ 0.07 for every N, where p̂_A is the field's own pooled
  F-smooth rate. This is a Poisson model with field-dependent quality.
- **P5 (cross-target / complementary-types check).**
  - The between-field overdispersion of g_A = (fraction of the 40 ℓ with N_ℓ > 0), relative to the
    binomial variance, lies in [0.5, 4].
  - Pr_A[g_A ≥ γ/2] ≥ 0.8.
  - The same-field pair ratio Pr[N_ℓ > 0 and N_ℓ′ > 0]/γ² lies in [0.9, 1.3].

## E3 — rank hypothesis RK (Design V)
- Setting:
  - A toy cubic NFS with the same field shape as E2. Its rows have columns for the sign, the rational
    primes ≤ B, the degree-one algebraic primes, and 16 quadratic characters.
  - The box is swept until rows ≥ 1.5·cols or the box ends. An instance is *valid* if rows ≥ 1.2·cols.
  - Dixon vectors are the exponent vectors of B-smooth x² mod N for uniform x, 20 per instance.
  - Grid: 48 bits with B = 1200, |a| ≤ 6000, b ≤ 600; and 56 bits with B = 3000, |a| ≤ 8000, b ≤ 1000.
    Classes {blum, tri} × seeds 0–4.
- **P6 (RK).** Over valid instances, at least 80% of Dixon vectors (rational support only, zero algebraic
  and character part) lie in the GF(2) row space of the NFS rows.
- **P7 (square-class structure).** Write Q_N for the squares mod N; the label p ∈ Q_N is harness-only.
  - P7a. Single-prime vectors e_p with p ∉ Q_N are in the span at rate ≤ 0.05.
  - P7b. Those with p ∈ Q_N that occur oddly in ≥ 3 rows are in the span at rate ≥ 0.7.
  - Reason: span ∩ (rational-only vectors) ⊆ {v : ∏p^{v_p} ∈ Q_N} whenever character-trivial
    even-exponent combinations are true squares in K. RK says equality holds.

Any deviation from these commands or parameters will be recorded as a disclosed deviation in
notes/hypD_results.md.

## Addendum A (2026-10-19, after the E1 run; frozen before E1x is run)
- **E1 as registered is uninformative for P2. This is a design error, disclosed.** The registered
  sampling grid had pooled model-expected hits of about 0.3: Pr[E_S] ≈ 10⁻⁶ against 4000 trials per record.
- The frozen E1 code also crashed on an ACCIDENTAL gcd(r, N) > 1. This happened in a 40-bit `tri`
  instance with a 13-bit factor. The fix resamples r and counts such events, and E1 was rerun
  unchanged otherwise. The crashed partial output is kept as `experiments/raw/hypD_e1_crashed_partial.jsonl`.
- **E1x (exact).** `src/hypD/e1exact.py`, config `experiments/configs/hypD_e1x.json`.
  - For each record it enumerates S × S and computes |{u/v} ∩ Q_N| / |Q_N| exactly.
  - Grid: classes {blum, tri} × bits {36, 40} × seeds 0–4 × k ∈ {2, 3}, G-modes random g = 1.0,
    random g = 0.5, and "qr" (only ℓ ∈ Q_N).
  - Command: `.venv/bin/python -m src.hypD.e1exact experiments/configs/hypD_e1x.json experiments/raw/hypD_e1x.jsonl 8`
- Predictions:
  - **P1′.** exact ≥ Theorem S bound in every record.
  - **P2′.** For the random modes with |S| ≥ 200: exact/model ∈ [0.7, 1.3].
  - **P2″.** For the qr mode with |S| ≥ 200: exact/(2^ω·model) ∈ [0.7, 1.3]. This is the positivity
    excess: when every element of S lies in Q_N, every u/v is a square.
