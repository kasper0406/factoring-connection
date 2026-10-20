# analytic: pre-registration (phase 17)

Agent `analytic`, 2026-10-20. Written and committed BEFORE E1/E2 were run. Code at the commit of
this file: `src/analytic/{hiary,hull,cubic,run_experiments}.py`; tests `tests/test_analytic.py`
(18 pass). The tests ran small instances only (N ≤ 1.7·10^7 for the cubic decoder, ≤ 40-bit for
Hiary); no E1/E2 data were produced.

## Why these experiments

The analysis (to be written in `notes/analytic_results.md`) makes two claims that these
experiments check.

1. **The exponent 1/3 of Hiary's algorithm is a convex-hull count.** Hiary's Algorithm 1
   (arXiv 1408.2608) is not an exponential-sum method. It tests a block [x − H, x + H] for divisors
   with one continued-fraction step and one quadratic equation, and H = ⌊(17N)^{−1/3}x⌋ is the
   scale over which the hyperbola y = N/x stays within a 2-D Diophantine (Dirichlet) window. Each
   dyadic range [X, 2X] therefore needs ≈ (17N)^{1/3} (ln 2)/2 blocks, whatever X is. The integer
   hull of {xy ≥ N} has Θ(N^{1/3} log N) vertices (Balog–Bárány, arXiv 2602.06897;
   Alcántara–Blanco–Criado–Santos, arXiv 2501.19193). Every divisor is a hull vertex. So
   "walk the hull" algorithms cost Ω(N^{1/3} log N), and Hiary's block count has the same order.
   E1 checks the block count against the continuum prediction, the per-dyadic-window constancy
   of both counts, and the divisor ⇒ vertex fact.
2. **Complete cubic Weyl sums mod N factor N (reduction R3).** For N = pq with p ≡ q ≡ 1 (mod 3),
   the nine values S(a^i b^j, N) = Σ_{x mod N} e(a^i b^j x³/N) determine p by exact Z[ω]
   arithmetic (proof in `src/analytic/cubic.py` docstring; full proof in the results note).
   S ≡ 0 when a prime factor is ≡ 2 (mod 3). The reduction needs S to absolute precision O(1),
   i.e. all ≈ ½log₂N + O(1) significant bits. So the cubic analogue of Hiary's polylog theta-sum
   algorithm cannot exist in the full coefficient range unless factoring is easy. E2 checks
   correctness and measures the noise threshold.

## E1: Hiary blocks and hull vertices (`run_experiments e1`)

- Hiary arm: N ∈ {prime_instance(n, s), balanced_semiprime(n, s)}, n ∈ {24, 28, …, 48},
  s = 0..9. The algorithm is `hiary_factor(N)` (operational, sees N only).
- Hull arm: balanced_semiprime(n, s), n ∈ {24, 28, 32, 36, 40}, s = 0..4. Vertices of the
  integer hull of {xy ≥ N} with x ∈ [⌈N^{1/3}⌉, ⌊√N⌋] are computed by brute force
  (`hull_vertices`, diagnostic), together with Hiary's block count over the same range
  (trial division off).
- **Predictions.**
  - P1.1: Hiary is correct on 140/140 (the factor is found for semiprimes; None for primes).
  - P1.2: for primes, blocks / predicted_blocks(N) ∈ [0.95, 1.05] for every instance, where
    predicted_blocks = ½(17N)^{1/3} ln(√N/(17N)^{1/3}).
  - P1.3: the least-squares slope of ln(mean blocks / ln N) against ln N (primes) is in
    [0.32, 0.345].
  - P1.4: the hull vertex count V satisfies V/(N^{1/3} ln N) ∈ [0.05, 1] at every size. The
    slope of ln(V/ln N) against ln N is in [0.30, 0.37]. For full dyadic windows with X ≥ 2N^{1/3},
    max/min of the per-window counts is ≤ 1.6 at each size, averaged over instances.
  - P1.5: every divisor of N in [⌈N^{1/3}⌉, ⌊√N⌋] is a hull vertex (25/25 instances).
  - P1.6 (report only): V / (Hiary blocks over the same range) ∈ [0.3, 3].
- **Decision rule.** E1 supports claim 1 iff P1.1, P1.2, P1.3 and P1.5 hold, and P1.4 holds up to at
  most one band miss. A miss is reported as a refutation of the stated constant or exponent.

## E2: cubic Weyl-sum reduction (`run_experiments e2`)

- Instances: balanced N = pq of n bits, n ∈ {14, 16, 18, 20, 22}, 20 per size with
  p ≡ q ≡ 1 (mod 3) (arm res1), and 20 per size with p ≡ 2 (mod 3) (arm res2, control). The
  harness generator is `semiprime_mod3`, seeded.
- Oracle: `weyl_oracle` (brute force, O(N), float64). The decoder `cubic_decode` sees N and
  oracle values only, and uses at most 6 attempts (one random pair (a, b) per attempt, 9 oracle
  calls each).
- Classification (rule 3): STRUCTURAL = factor obtained from the Z[ω] gcd; ACCIDENTAL = gcd(a, N)
  or gcd(b, N) nontrivial. The accidental null is ≈ 2·(1/p + 1/q) per attempt.
- **Predictions.**
  - P2.1: res1 success ≥ 98/100 overall, all STRUCTURAL except accidental gcd hits. An attempt
    fails only if the class map (i, j) ↦ (χ_p, χ_q)(a^i b^j) is not onto μ₃², with probability
    33/81. Six attempts fail with probability 0.0046.
  - P2.2: res2 controls: max |S| < 10^{−6} for every oracle call, and 0 STRUCTURAL successes.
  - P2.3: noise. 10 res1 instances per size; each oracle value gets uniform noise in [−δ, δ],
    δ ∈ {0.01, 0.03, 0.1, 0.3, 1, 3, 10}, with 3 noise seeds. δ*(n) = the largest δ with
    success rate ≥ 50%. Prediction: δ* ∈ [0.03, 3] at every size, and the slope of ln δ* against
    ln N is in [−0.15, 0.15]. That is, constant absolute precision: the decoder needs all but
    O(1) of the ≈ ½log₂N significant bits of S.
- **Decision rule.** E2 supports the reduction iff P2.1 and P2.2 hold. P2.3 calibrates the
  precision claim; a slope outside the window is reported as a correction to it.

## Disclosed exception to the exact-arithmetic rule

The oracle S(a, N) is a real analytic quantity. It is computed in float64, with rounding error
≤ N·2^{−45} < 10^{−6} for N < 2^{22}, then converted exactly to a Fraction. Everything
downstream (the decoder, Z[ω] Euclid, rounding to norm N) is exact. Floats elsewhere are used for
fits only.

## Budget

E1 ≈ 10 min, E2 ≈ 10 min, single process, < 200 MB RSS. Run one at a time after checking
`memory_pressure`.
