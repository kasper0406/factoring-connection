# floorhard: pre-registration

Agent `floorhard`, 2026-10-20. Written and committed BEFORE any of E1–E3 was run. No code for
E1–E3 existed when this file was written, so no data have been seen.

## Claims under test (analysis in notes/floorhard_results.md)

Notation: an exact floor-sum oracle returns F(f, D; h₁, h₂) = Σ_{h₁≤h<h₂} ⌊f(h)/D⌋ for an integer
polynomial f and integer D ≥ 1. Window identity (W): for integers 0 ≤ u < v ≤ D,
#{h ∈ [h₁, h₂) : f(h) mod D ∈ [u, v)} = F(f − u, D) − F(f − v, D) (two calls).

- **(A) Full-period closed form.** For odd squarefree m and gcd(a, m) = 1,
  Σ_{h=0}^{m−1} ⌊ah²/m⌋ = a(m−1)(2m−1)/6 − (m−1)/2 + Σ_{d | m, d ≡ 3 (mod 4)} (a/d)·h*(−d),
  with h*(−3) = 1/3 and h*(−d) = h(−d) for d > 3 (d squarefree, so −d is fundamental).
  Consequence: a full-period oracle yields h(−m) (m ≡ 3 mod 4) or h(−3m) (m ≡ 1 mod 4) from a few
  multipliers a, and then a factor of m via ambiguous forms (randomized).
- **(B) #SubsetSum ⇒ quadratic congruences ⇒ degree-2 floor sums.** For weights w₁..w_k and target
  T, an explicit (a, b = ∏ p_i, window) makes #{x in window : x² ≡ a (mod b)} equal the number of
  subsets with sum T (parsimonious).
- **(C) Taylor embedding.** For N = A·b²·t³ (+ e, 0 ≤ e < E) and x = bt, the window counts of the
  degree-j Taylor polynomial P_j of N/(x+h) on [0, c') equal #{h < c' : A h² ≡ r (mod b)} for every
  residue r, once t ≥ t₀(A, b, c', j). Hence O(log m) calls of the Taylor-class oracle factor any m
  (bisection on h² ≡ 1 mod m).

## E1: full-period formula and class-number factoring (`run_experiments e1`)

- E1a: every odd squarefree m with 3 ≤ m < 3000; a ∈ {1, 2, 3, 5, 7, 11} with gcd(a, m) = 1.
  Left side by brute force; class numbers by `src.nonsmooth.forms.class_number` (reduced forms).
- E1b: 30 balanced semiprimes m = pq (seeded, p ≠ q odd primes > 3), 15 with m ≡ 3 (mod 4) and 15
  with p ≡ q ≡ 3 (mod 4), at 16–22 bits. Oracle = brute-force full-period floor sum (numpy, exact
  int64). Decoder: random multipliers a with (a/M) = 1 and (a/3) = 1 when M = 3m; recover h(−M) as
  the mean of the two observed values of the class-number part; check it against
  `forms.class_number`; factor with ambiguous forms (`src.nonsmooth.equiv.reduce_class_number` on
  disc −4M with h(−4M) = h(−M)(2 − (−M/2))).
- **Predictions.**
  - P1.1: 0 failures of (A) in E1a.
  - P1.2: decoded h(−M) equals the reduced-form count in 30/30.
  - P1.3: a nontrivial factor of m is found in ≥ 29/30 (randomized; ambiguous-form step).
  - P1.4: median number of full-period oracle calls per instance ≤ 4, max ≤ 12.

## E2: #SubsetSum reduction (`run_experiments e2`)

- E2a (brute-force floor sums): 40 random instances, k = 3, weights in [1, 6], target uniform in
  [0, S]; the two floor sums of (W) computed by brute force over the window.
- E2b (structural core): 200 random instances, k ∈ {4, …, 10}, weights in [1, 2^k]; count via
  enumeration of the 2^k CRT square roots (diagnostic) and compare with brute-force subset-sum count.
- **Predictions.** P2.1: equality in 40/40 (E2a). P2.2: equality in 200/200 (E2b).

## E3: Taylor embedding and bisection factoring (`run_experiments e3`)

- E3a: for j ∈ {2, 3, 4}, 60 random (A, b, c') with b ≤ 200 (squarefree odd and non-squarefree
  allowed), A ≤ 5, c' ≤ 3b; t = t₀ from the proof; N = A b² t³; x = bt. Check that the block length
  H = ⌊x/(2B)⌋ (B = ⌈N^{1/(j+2)}⌉) is ≥ c', and that for every residue r mod b the window count of
  P_j equals #{h < c' : A h² ≡ r (mod b)} (brute force).
- E3b: prime variant: N = A b² t³ + e prime (random e < E from the proof), same check, 20 instances
  (j = 2).
- E3c: 40 odd m ≤ 20 000 with at least two distinct prime factors (seeded), j = 2:
  factor m by bisection on h² ≡ 1 (mod m) using only window counts of the embedded Taylor class
  (oracle evaluated by brute force: model cost only). Count oracle calls.
- **Predictions.** P3.1: 60/60 (E3a). P3.2: 20/20 (E3b). P3.3: 40/40 factored, oracle calls
  ≤ 2⌈log₂ m⌉ + 4 for each.

## Disclosed limitations

- All oracles are brute force: these experiments test identities and reductions, never costs.
- Sizes are toy. Toy factoring is a correctness check only.
