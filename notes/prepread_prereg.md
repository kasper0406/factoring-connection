# prepread pre-registration (agent `prepread`, phase 18, 2026-10-20)

Frozen before any experiment below is run. Code: `src/prepread/*`. Raw output:
`experiments/raw/prepread_e{1,2,3,4,5}.jsonl`. Memory budget ≤ 2 GB RSS, one heavy job at a time.

## Questions

1. **periodfind O1.** Ψ(y) = e(h·Φ(y)/N) with Φ a short SLP of huge degree. By CRT and Theorem T, a heavy
   large-order peak needs, at p (and at q), a large value of
   A_p(χ, c) = (1/(p−1)) Σ_{y∈F_p^*} χ(y) e(cΦ(y)/p) for a character χ of large order.
2. **slpcirc dense-torus question.** Root densities mod p of short N-independent SLPs, in particular
   compositions (x² + c_k)∘⋯∘(x² + c_1).

## Analytic claims to be proved or refuted in notes/prepread_results.md (not by these experiments)

- Weil threshold: |A_p(χ, c)| ≤ w·√p/(p−1), w = cyclic exponent width of Φ mod (y^{p−1} − 1); vacuous at
  w ≳ √p.
- Reduction: a computable f on F_p^* with |⟨f, χ⟩| ≥ τ for χ of order d gives discrete logs modulo d in
  time poly(cost(f), log p, 1/τ). Hence O1 success ⇒ DL in F_p^* for a positive fraction of primes.
- Galois framing: root count = fixed points of Frob_p; the dense regime lies outside effective Chebotarev;
  Galois groups alone impose no constraint (constructions with trivial, cyclic and S_D groups).

## E1: large-prime-order character correlations of high-degree short SLPs (float64 FFT, disclosed)

- Primes: 4 seeded primes in (2^13, 2^14) (all c ∈ F_p^*) and 2 seeded primes in (2^15, 2^16) (256 seeded c).
- For p: ℓ = P⁺(p−1), v = v_ℓ(p−1), J_p = {j : ℓ^v ∤ j} (characters whose order is divisible by ℓ).
  g = least primitive root (HARNESS uses the factorisation of p−1; this is the prime p's own data).
- Statistic: M_p(Φ) = max_{c} max_{j∈J_p} |A_p(j, c)| with A_p(j, c) = (1/(p−1)) Σ_t e(cΦ(g^t)/p) e(−jt/(p−1)).
  Baseline B_p = median of M_p over 10 seeded uniformly random maps F_p^* → F_p.
- Families:
  - F1: 20 seeded random SLPs each of size 10, 14, 18 (ops +, −, × on earlier nodes; leaves x, 1, 2, 3);
  - F2: compositions (x² + c)^{∘k}, c ∈ {1, −1, 3}, k ∈ {8, 12, 16, 20};
  - F3: x^{2^k} + x (k = 10, 16, 20), D_{2^k}(x) Dickson (k = 10, 16), x^{M} with M = lcm(1..30);
  - F4 (controls, HARNESS): planted-DL phase y ↦ ⌊p·log_g(y)/(p−1)⌋; Pohlig–Hellman SLP
    PH_e(y^{m}) (§E4) at primes where it applies.
- A member is "nondegenerate at p" if |Φ(F_p^*)| ≥ p/4. Predictions:
  - **H1.** Every nondegenerate F1–F3 member at every prime has M_p(Φ)/B_p ≤ 2.
  - **H2 (Weil sanity).** For members with w < √p (only possible for small degree), |A_p| ≤ w√p/(p−1) for
    all j, c: 0 violations.
  - **H3.** F4 planted-DL control: M_p/B_p ≥ 20 at every prime.
- **Kill criterion (for the barrier picture):** some non-control F1–F3 member with M_p/B_p ≥ 10 at ≥ half
  of the primes. That would be a candidate O1 observable.

## E2: exhaustive search over quadratic compositions (C helper, exact integer arithmetic)

- C = {−2, −1, 0, 1, 2}; f_c = (x² + c_k)∘⋯∘(x² + c_1); fibre Z_p(c) = #{x ∈ F_p : f_c(x) = 0}, by backward
  preimage propagation from {0}, all 5^k sequences.
- Per-prime: for b ∈ {7, 8, 9, 10}, 4 seeded primes in (2^b, 2^{b+1}), k = b + 2. Record M*(p) = max_c Z_p(c)
  and the argmax sequences.
- Uniform: for b ∈ {9, 10}, 8 seeded primes in (2^b, 2^{b+1}) jointly, k ∈ {b + 2, b + 4}. Record
  U* = max_c min_p Z_p(c).
- Predictions:
  - **H4a.** M*(p) ≤ max(2k², 4·2^{max(v₂(p−1), v₂(p+1))}) at every prime.
  - **H4b.** median_p M* at b = 10 is ≤ 3 × median_p M* at b = 7 (critical-branching extreme-value
    heuristic ≈ 1.8×; a density-bounded-below mechanism would give ≈ 8×).
  - **H5.** U* ≤ k in every uniform run.

## E3: Chebotarev regime versus silent regime for x² + 1 (exact)

- All primes in (2^14, 2^16). For k = 1..24, Z_p(k) = #{x : f^{∘k}(x) = 0}, f = x² + 1.
- Haar-wreath prediction: fraction of primes with Z_p > 0 equals p_k, p_0 = 1, p_k = (1 − (1 − p_{k−1})²)/2;
  E Z = 1.
- **H6.** For every k ≤ 24: |frac_k − p_k| ≤ 3√(p_k(1 − p_k)/#primes) + 0.01, and mean Z ∈ [0.5, 2]. (For
  k ≥ 15 the degree 2^k exceeds p: this is the Chebotarev-silent regime; a pass means Frobenius at these
  small primes still looks Haar-distributed on the tree.)

## E4: operational DL-from-correlation reduction (positive controls)

- Reduction code sees only p, g, y and query access to f; it never computes logs itself.
- (a) Planted f(y) = e(⌊p·log_g(y)/(p−1)⌋/p) (HARNESS-built), 6 primes from E1, 20 random targets each.
- (b) The same reduction applied, at each E1 prime, to the F1–F3 member with the largest M_p/B_p (its best
  (c, j) pair is taken from E1), 20 random targets.
- **H7.** (a) 20/20 recoveries at every prime. (b) 0/20 recoveries at every prime (chance level is
  ≈ 20/(p − 1)).
- The Pohlig–Hellman SLP construction of the results note (correlation 1 with a character of 2-power order
  at primes whose odd part of p − 1 divides the SLP's exponent) is checked by a unit test, not here.

## E5: random-SLP root-density census (exact)

- 2000 seeded random SLPs of size 20 (as F1); 12 seeded primes in (2^13, 2^14); δ_p = #roots/p.
- **H8.** No SLP has min_p δ_p ≥ 1/64 over the 12 primes, except SLPs whose output vanishes identically
  at some of these primes or is reported as a structured hit (classified by hand in the results note).

## Disclosure rules

Any deviation (prime counts, sizes, time-outs) will be recorded as a disclosed deviation. Post-hoc analyses
will be labelled as such.
