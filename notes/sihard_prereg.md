# sihard pre-registration (phase 20)

Agent `sihard`, 2026-10-21. Written and committed before any experiment below was run. The experiments are small,
exact sanity checks of the hardness analysis. None of them can confirm an asymptotic statement.

Code: `src/sihard/` (to be written after this file is committed). Raw output: `experiments/raw/sihard_E{1,2,3,4}.jsonl`.
Memory budget ≤ 2 GB RSS, one job at a time, `memory_pressure` checked first.

## E1. Smooth-number counts in intervals near kN versus random intervals

*Question.* Prop SI's heuristic H_S treats the position kN as random. Do y-smooth counts in I_k = [kN, kN + h] behave
like counts in intervals at random positions of the same size?

*Design.*
- N = `balanced_semiprime(48, seed)`, seeds 0..9.
- y ∈ {2^10, 2^12}. Target mean μ ∈ {2, 8, 32}.
- h = round(μ / ρ(u_mid)), with u_mid = log(100.5·N)/log y and ρ the Dickman function (numerical; design parameter only).
- Near-kN sample: k = 1..200, intervals [kN, kN + h].
- Control sample: 200 intervals [x, x + h], x uniform in [N, 201N], from `random.Random(f"sihard-E1-{seed}-{y}-{mu}")`.
- Smoothness is decided exactly: numpy uint64 cofactor arrays, divided by every prime ≤ y until not divisible
  (all values < 2^63).

*Predictions (pass/fail).*
- **P1.** In each (y, μ) cell, pooled over the 10 N, the ratio of the near-kN mean count to the control mean count lies
  in [0.9, 1.1].
- **P2.** In each cell the index of dispersion (variance/mean) of the near-kN counts lies in [0.7, 1.6], and differs
  from the control's by at most 0.3.
- **P3.** In the μ = 2 and μ = 8 cells, the fraction of empty near-kN intervals differs from the control fraction by at
  most 0.05.
- *Exploratory (no pass/fail).* The ratio of the observed control mean to the Dickman prediction μ.

## E2. Johnson-radius certificate decoders: exact averaged counts

*Question.* Theorem J (to be proved in the results note) bounds the number of y-smooth *certified* outputs of a
Coppersmith/Boneh/Howgrave-Graham call, averaged over the interval position x. Check the rigorous inequality exactly,
and see which moduli P do best at toy scale.

*Design.*
- X = 2^40, y = 2^8, V ∈ {2^8, 2^10, 2^12}. Average over x ∈ [X, X + 2^16).
- Moduli P (products of primes ≤ y, exponents allowed): the family registered in code before running:
  (i) single smooth d ≈ V (the divisor trick, π := log P = log V);
  (ii) primorial-type P = ∏_{p ≤ z} p with π ≈ c·log V, c ∈ {1.5, 2, 4, 8};
  (iii) highly-composite-type P = ∏_{p ≤ z} p^{⌊log B/log p⌋} with π ≈ c·log V, c ∈ {2, 4, 8, 16};
  (iv) 5 random subsets of the factor base with π ≈ 4 log V.
- S from the idealised Johnson condition: log S = √(log V · log P) (rounded up to an integer S).
- Measured: A = average over x of #{z ∈ [x, x+V] : gcd(z, P) ≥ S, z y-smooth}, computed exactly.
- Rigorous bound: B = ((V+1)/M)·#{z ∈ [X, X+M+V] : gcd(z, P) ≥ S, z y-smooth}, M = 2^16, together with the
  divisor-sum form ((V+1)/M)·Σ_{d|P, d≥S} #{m : dm ∈ [X, X+M+V], m y-smooth}.
- Baseline: dens = density of y-smooth integers in [X, X + M + V].

*Predictions.*
- **P4 (exact).** A ≤ B ≤ (divisor-sum form) in every cell.
- *Exploratory.* The best A/dens over the registered family, versus the divisor-trick value (i).

## E3. Oblivious two-list merging

*Design.*
- X = 2^36, h ∈ {2^6, 2^10}.
- Lists A and B: 2000 distinct 2^8-smooth integers each, drawn from [2^18, 2^18·√2) by `random.Random("sihard-E3-…")`, so
  every product lies in [2^36, 2^37).
- Measured: the average over x ∈ [X, 2X) of #{(a, b) : ab ∈ [x, x + h]}, computed exactly by counting products.
- Bound: |A||B|(h+1)/X.
- Bohr version: N uniform in [2^30, 2^31), target {ab : ab mod N ∈ [0, h], ab ≥ N}, average over 4096 seeded N.
  Bound 2(h+1)(1 + ln K)|A||B|/X_N, K = ⌈max(ab)/2^30⌉.

*Predictions.*
- **P5 (exact).** The measured average never exceeds the bound.
- **P6.** In the interval version, measured/bound ∈ [0.9, 1.0].

## E4. Explicit-factor-base reduction (X3C → S-smooth integer in an interval)

*Design.* 300 random X3C instances: q ∈ {2, 3, 4} (universe 3q), m ∈ {q+1, …, q+4} random 3-subsets, seeded. Build
(S, x, h) by the reduction of the results note (random primes found with a seeded RNG and `sympy.isprime`). Count the
S-smooth integers in [x, x + h] by exhaustive enumeration of all exponent vectors with product ≤ x + h, independently of
the reduction logic. Count exact covers by brute force.

*Prediction.* **P7 (exact).** The two counts are equal in all 300 instances.

## Deviations
Any deviation from this file will be recorded in the results note as a disclosed deviation.
