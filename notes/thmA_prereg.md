# thmA pre-registration: E1 / E1b, an end-to-end toy run of the Theorem A pipeline

Agent `thmA`, phase 10, 2026-10-19. This file is committed before any E1/E1b run.

**Purpose.** A correctness smoke test of the whole final pipeline (notes/thmA_proof.md §3) at toy scale:
- the random design (a‴_k; μ) with a random prime multiplier w, so F(g₀, g₁) = N·w;
- Dixon descent of uniform random squares: split R = u₀/v₀, shape test, special-q relation for every large prime,
  batch smoothness per box;
- LV6 columns: sign, rational primes, ONE parity column, first-degree algebraic primes r ∤ f_d·disc f, and 40
  quadratic characters;
- GF(2) kernel, Hensel square root in Z[β], congruence, gcd.

This is **not an asymptotic test**. Toy sizes say nothing about L[1/2] constants.

**Code (frozen at this commit).**
- `src/thmA/{arith,numfield,smooth,pipeline,params}.py`: operational. The pipeline sees only N, the public
  config and a seed.
- `src/thmA/harness_e1.py`: harness. It knows the factorisations and uses them only to check outputs and to
  annotate ω.
- Config `experiments/configs/thmA_e1.json`. Tests `tests/test_thmA.py`.
- Command:
  `.venv/bin/python -m src.thmA.harness_e1 experiments/configs/thmA_e1.json experiments/raw/thmA_e1.jsonl`,
  then `.venv/bin/python -m src.thmA.analyze experiments/raw/thmA_e1.jsonl`. The analysis script is written
  after the run. It implements exactly the statistics defined below.

**Calibration (disclosed).**
- Seeds 900–906 (balanced semiprimes at 32/40/48/50 bits, one 3-prime N, one unbalanced N, one prime power)
  were used to choose S = 12, min_large = 1 and the B tiers.
- They are not in E1/E1b.
- Calibration observations: every run factored correctly. The special-q success rate was about 5–10%.
  Accidental gcd(r, N) hits ended the 32-bit and the 3-prime 48-bit runs early.

**Toy deviations from Theorem A (disclosed).**
- d = 3, μ = 1/3.
- w is uniform among primes in (2M, 4M], not (W, 2W] with W = N^{μ+ν+η₀/2}. So Lemma W's analytic role is not
  tested. The algebra with F(g₀, g₁) = Nw is tested.
- H = 4. The h/Kern/level/window floors are not imposed.
- Admissible split: each half is (B-smooth part) × (at most k_large = 2 distinct primes in (B, 4B]). Each row
  must contain at least one large prime (min_large = 1), so every row uses a special-q atom. Theorem A's shape
  is "≤ k distinct primes of the dyadic J ≈ N^{1/(2k)}, no smooth part". That shape is infeasibly rare at
  30–50 bits: polylog^{−k} with large constants.
- Shape test by Pollard–Brent on the split halves, not Lenstra–Pomerance. It is never applied to N.
- Primality by deterministic-base Miller–Rabin, not AKS.
- Polynomial factorisation and root finding mod p by python-flint.
- Toy character primes of 62 bits; t_χ = 40, not O(d log N).
- Kernel by Gaussian elimination, not Wiedemann. Uniform kernel vectors are random XORs of a basis.
- Row target: (#odd columns that occur) + 40 + 20.
- Dependencies are tested by lv6 §3.7's refinement: up to 8 pre-sampled uniform kernel vectors per round, each
  skipped if it lies in the span of those already tested.
  - The algorithm's answer is the first fruitful one.
  - For statistics, all 8 are tested (stats_all_deps).
  - A round with no fruitful dependency is followed by fresh rows in the same field (≤ 6 rounds).
- Composites below N₀ = 2^20 are factored by trial division.
- E1b only (diagnostic): `continue_on_accidental` records an accidental factor, treats the object as invalid
  and continues on the same N. The purpose is to measure DD statistics on 3-prime N. E1 runs the algorithm as
  specified: any proper gcd ends the split.

**Instances.**
- E1, 48 instances:
  - balanced semiprimes at 30/34/38/42/46/50 bits × seeds 0–2, plus 50 bits × seeds 3–4 (20);
  - 3-prime N at 42/45/48/50 bits × seeds 0–2 (12);
  - unbalanced N with an 8- or 10-bit factor (4);
  - small prime × balanced semiprime (4);
  - p²q (4);
  - p², p³, (pq)², a prime (4).
- E1b, 12 instances (diagnostic mode):
  - 3-prime N at 46/50 bits × seeds 10–13;
  - semiprimes at 40 bits × seeds 10–13.

**Predictions (registered).**
- **P1 (correctness).** In E1, all 48 instances output the correct complete factorisation: the product equals N
  and every factor is prime. Any cap or exception counts as a failure. Prediction: 48/48. Kill: ≥ 1 wrong output
  (a failure from a cap is reported separately).
- **P2 (row identity).** 0 ROW_IDENTITY_VIOLATION over all rows: r² ≡ q·t^{−E}·∏ψ(z_P)^{e}, with the sign of
  u₀v₀ in q.
- **P3 (congruence).** 0 CONGRUENCE_VIOLATION: X_S² ≡ Y² (mod N) for every square dependency.
- **P4 (characters).** Pooled over E1 + E1b, at least 90% of tested kernel vectors are squares in K.
- **P5 (Theorem DD).** Statistics over tested square dependencies, pooled, E1 + E1b, by ω of the split input:
  - the fruitful fraction lies within |z| ≤ 3 of 1 − 2^{1−ω}, i.e. 1/2 for ω = 2 and 3/4 for ω = 3;
  - binomial z;
  - by DD and lv6 §3.7 the outcomes within a round are i.i.d.
  - The registered test needs ≥ 30 dependencies per class. Fewer means "not testable".
- **P6 (accidental factors).** Factor discoveries are classified STRUCTURAL_DD (the final gcd) or ACCIDENTAL_*
  (gcd with r, a split half, ℓ, G, F, field data).
  - In E1, the observed number of ACCIDENTAL split exits is compared with the null
    λ = Σ_splits (3·r_trials + 2·sq_candidates + 4·fields_sampled)·Σ_{p | n} 1/p.
  - Prediction: observed ≤ λ + 3√λ + 1.
- **P7 (preprocessing).** The small factors of the unbalanced and small × semi instances are found by trial
  division. p², p³ and (pq)² are handled by the perfect-power step. The prime is recognised. All 12 succeed.
- **Descriptive (not predicted).** Per split: r-trials, shape acceptance, descent failures, special-q success
  rate (γ̂), fields sampled, rows, columns, kernel dimension, time.
