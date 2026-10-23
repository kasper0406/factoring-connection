# genusleak pre-registration (phase 24)

Agent `genusleak`, 2026-10-23. Written and committed before any experiment code exists or runs.

## Question
Can an ambiguous class of Cl(−4kN) (or an ambiguous/square form of discriminant 4kN) be reached in subexponential
time without h, its smooth part, or an element's order? The theory (notes/genusleak_results.md) carries the main
weight. The experiments below are exact toy computations. They check proved lemmas against real computation and
measure the constants of the exponential routes. Toy sizes cannot show asymptotics; any factoring success is a
correctness check only (rule 8).

## Common setup
- Instances: `src/common/instances.py::balanced_semiprime(n, seed)`, seeds 0, 1, 2, … in order. Where a congruence
  filter is stated (e.g. N ≡ 1 mod 4), instances failing it are skipped and the next seed is used. The skipped seeds
  are recorded.
- Operational code (`src/genusleak/forms.py`, `cf.py`, the algorithms in `cayley.py`, `sqrtred.py`) sees N only.
  Ground truth (p, q, class numbers from PARI via cypari2) lives in `src/genusleak/harness.py` and is never imported
  by operational code. PARI's class-group routines factor D internally, so they are harness-only.
- Exact integer/rational arithmetic, except E5, which uses mpmath (60 digits) for L-value approximations and fits.
- Raw output: `experiments/raw/genusleak_E{1..5}.jsonl`. Memory ≤ 2 GB; one job at a time.

## E1: Cayley-graph walks (return and collision statistics), D = −4N, N ≡ 1 (mod 4)
- n ∈ {20, 24, 28}, 20 instances each. S = prime forms f_ℓ = (ℓ, b, c) (b the least nonnegative root) for the first
  6 odd primes ℓ with (D/ℓ) = 1. The walk is the simple random walk on the subgroup H = ⟨S⟩ ⊂ Cl(D), with steps
  f_ℓ^{±1}, each with probability 1/12. H is enumerated by BFS on reduced forms.
- W_t = number of closed walks of length t at the identity (exact integer DP on H), for t ≤ 40.
- W_t^free = the same count in Z^6 (exact).
- Hidden genus signal (harness): G_t = (1/|H|)·[(Σ_s χ_p(s))^t + (Σ_s χ_q(s))^t]/12^t, with χ_p(f_ℓ^{±1}) = (ℓ/p).
  P_t = W_t/12^t.
- Derived sample sizes (signal-to-noise 1, even t only):
  - M_ret(t) = P_t/G_t²;
  - M_coll(t) = √(2·P_{2t})/|G_{2t}|.
- **P1 (girth; proved lemma, implementation check).** W_t = W_t^free for every t < g₀ = ⌈log N / log ℓ_max⌉, in
  every instance.
- **P2 (return estimator; proved lower bound + predicted range).** min over even t ≤ 40 of M_ret(t)/|H| is ≥ 1/4 in
  every instance (proved) and ≤ 100 in the median instance per n (predicted).
- **P3 (collision estimator).** min over t ≤ 20 of M_coll(t)/√|H| is ≥ 1/√2 in every instance (proved) and ≤ 50 in
  the median instance per n (predicted).
- Reported, not predicted: the first t with W_t > W_t^free, compared with g₀; and whether H = Cl(D) (harness).

## E2: SQUFOF density versus infrastructure jumps (principal cycle of discriminant 4N)
- n ∈ {28, 32, 36}, 30 instances each, k = 1, no congruence filter. Continued fraction of √N over one full period
  (exact recurrence).
- A *proper square* is an even index i with Q_i = r² for which the reverse cycle from the square-root form reaches an
  ambiguous form giving a nontrivial factor of N. The classification uses the actual reverse-cycle walk.
- Statistics:
  - i₁ = the index of the first proper square (the SQUFOF forward cost);
  - ρ = (# proper squares in one period)/(period length).
  - Uniform random jumps into the cycle (an idealisation of giant steps to random distances, indexed by position)
    need 1/ρ forms in expectation.
- **P4.** The median over instances of i₁·ρ lies in [0.5, 2] for every n. Sequential SQUFOF and random jumps cost the
  same up to a constant.
- **P5.** The least-squares slope of log₂(median i₁) against n (three points) lies in [0.18, 0.32]. Heuristic
  exponent 1/4.

## E3: the square-root reduction (oracle → factor), D = −4N, N ≡ 1 (mod 4)
- n ∈ {24, 32}, 10 instances each, 40 queries each.
- S = the first 10 split odd primes. Sample x = ∏ f_ℓ^{e_ℓ}, with e_ℓ uniform in [0, N²).
- Oracle (harness, brute force on H = ⟨S⟩): given y = x², it returns the root z of y that is least in the
  lexicographic order of reduced (a, b, c). This oracle is deterministic and does not see x. The reduction outputs
  the reduced form of z·x⁻¹ and reads a factor from it if it is a non-public ambiguous form.
- **P6.** Pooled over instances, the per-query success rate is within 0.5 ± 0.08 at each n. The predicted value
  1 − |A_pub|/|Cl[2]| = 1/2 holds when H = Cl(D); instances with H ≠ Cl(D) are reported separately.
- **P7.** Every claimed factor is correct (100%).

## E4: the midpoint of the principal cycle (distance R/2)
- n = 28: 100 instances; n = 32: 50 instances. Multipliers k ∈ {1, 2, 3, 5, 6, 7, 10, 11}. CF of √(kN) over one
  period of length L. If L is even, g = gcd(Q_{L/2}, N).
- **P8.** The fraction of N with a nontrivial g at k = 1 lies in [0.2, 0.8] at each n.
- **P9.** The fraction of N with a nontrivial g for some k in the list is ≥ 0.9 at each n.

## E5: what an analytic approximation of h buys, D = −4N
- n ∈ {32, 40, 48}, 20 instances each.
- ĥ(x) = (√|D|/π)·∏_{ℓ≤x}(1 − (D/ℓ)/ℓ)^{-1}, for x = 2^8, 2^9, …, 2^20. Exact h from PARI `qfbclassno` (harness).
- Δ(x) = |h − ĥ(x)|. Cost model C(x) = π(x) + 2√Δ(x): Euler factors, plus BSGS on a window of width 2Δ.
- **P6′.** The median over instances of the fitted slope d log Δ / d log x (least squares over the grid) lies in
  [−0.75, −0.25] at each n.
- **P7′.** The median over instances of log(min_x C(x))/log|D| lies in [0.15, 0.30] at each n.

## Decision rules
- Each P is reported as met / not met; nothing is re-run to rescue it.
- A failed "proved" prediction (P1, or the lower-bound halves of P2/P3, or P7) means an implementation bug or a wrong
  lemma. It is investigated and disclosed.
- None of E1–E5 can support a subexponential claim. They can refute a stated lemma or calibrate constants.
