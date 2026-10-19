# specialfam pre-registration (phase 7)

Agent `specialfam`, 2026-10-19. Written and committed BEFORE any run below. Code and tests are frozen
at this commit: `src/specialfam/{polyprod,normbox,copper,smoothcopp,harness,run_experiments}.py`,
`tests/test_specialfam.py` (28 tests, all pass), config `experiments/configs/specialfam.json`.

**Disclosed pre-commit activity.** (i) Smoke tests at 39–80 bits while writing the code (a handful of
instances; not data). (ii) One timing probe at 100 bits (one normbox instance, 3.1 s; one Algorithm-D
instance, 0.3 s), used only to choose sizes. (iii) The finite-lattice predicted slopes below were
computed from `best_params` on synthetic N = 2^{n−1}+2^{n/2}+1 (no instances, no factors).

## Algorithms (operational: N + public family parameters only)

- **A (norm box, Theorem A).** Public: N, η, and the family's size exponent β (used only to start the
  dyadic band scan at k = ⌊β·bitlen(N)⌋ − 3). Bands k: box a ∈ [0, ⌊2^{k/2}⌋], b ∈ [1, ⌊2^{kη}⌋];
  block products Π_{j<m}((mv+j)² + b²) by one monic remainder tree (resultant packaging, `polyprod`);
  gcd = N resolved by scanning the block. No inversion, no primality test.
- **B (Theorem S).** Public: N, smoothness bound B, lower bound Dmin on g, the balanced band. PS finds
  the B-smooth part H of N²−1; divisors d ∈ [Dmin, B·Dmin) of H, all r² ≡ 1 (mod d); rigorous
  Coppersmith chunks (lattice dimension ≤ 9).
- **D (Corollary W).** Public: N, balanced band. g = gcd(p−1, q−1) arbitrary. Dovetail B = 2^t with a
  chunk budget 2^{⌈t/2⌉} per divisor, residue r = 1.

## Experiments

| id | what | cells | instances |
|---|---|---|---|
| E1 | A on p = a²+b², b ≤ p^{0.119} (b uniform), q random prime; β = 1/2 (p<q<2p) and β = 2/5 | n ∈ {40,48,…,96,100} | 10 per (β, n) = 180 |
| E2 | cost slopes from E1 | — | — |
| E2c | A on unstructured balanced semiprimes (control), bands up to n/2+1 | n ∈ {40,56,72} | 20 per n = 60 |
| E3 | Theorem G on two more families: HKL-type box (public centre, side 2^{⌈0.265 kp⌉}); near-square p = a²+k, k ≤ 2^{⌈0.238 kp⌉} | n ∈ {48,64} | 10 per (fam, n) = 40 |
| E4 | B, g = 2^a·3·5 ≈ N^γ, g ∣ gcd(p²−1, q²−1) with random independent CRT signs per prime power | γ ∈ {0.12, 0.18}, n ∈ {48,…,96} step 8 (+100) | 8 per cell |
| E5 | D, g = 2·(prime) ≈ N^γ (not smooth), p ≡ q ≡ 1 mod g | γ ∈ {1/10, 1/8, 1/6, 1/5}, n ∈ {48,…,96} step 8 | 6 per cell |

Command: `.venv/bin/python -m src.specialfam.run_experiments e1|e2c|e3|e4|e5|summarize`. Raw:
`experiments/raw/specialfam_e*.jsonl`; summary `experiments/analysis/specialfam_summary.json`. Slopes are
least-squares of log2(median per n) against n (floating point used only for the fit).

## Predictions and decision rules (frozen)

- **P1 (correctness).** E1, E3, E4, E5: every promised instance is factored (100%). Any failure is a
  bug or a theorem error; it will be reported and investigated, not removed.
- **P2 (E1 cost, counts).** Slope of log2(units of the successful band) per bit of N:
  β = 1/2: predicted 1/8 + η/4 = 0.15475, accept [0.145, 0.165];
  β = 2/5: predicted β(1/4 + η/2) = 0.1238, accept [0.114, 0.134].
  (Counts follow a deterministic formula; this is a sanity check of the implementation.)
- **P3 (E1 wall time).** β = 1/2 wall-time slope over the upper half of sizes in [0.13, 0.21].
  Pollard–Strassen on the same p would have slope 0.25.
- **P4 (classification).** ≥ 90% of E1 hits are STRUCTURAL-PROMISED (the block holding the promised
  (a, b)); the rest are STRUCTURAL-OTHER. Null expectation for a hit at another leaf ≈ V/p per prime.
- **P5 (control).** E2c: 0 wrong factors. Found rate ≤ 10%, and every find should have p (or q) with a
  representation b ≤ p^{0.119} (structural) or a leaf divisible by p (accidental, rate ≈ V/p).
- **P6 (E4 chunks).** Slope of log2(median chunks) within ±0.03 of the finite-lattice prediction
  (dimension ≤ 9): γ = 0.12 → 0.153, γ = 0.18 → 0.108. Asymptotic (dimension → ∞): 1/4 − γ.
- **P7 (E5).** Slope of log2(median max(ps_units, chunks)) within ±0.035 of max(γ/2, s_fin(γ)), with
  s_fin = 0.175, 0.153, 0.108, 0.081 for γ = 1/10, 1/8, 1/6, 1/5 (so predictions 0.175, 0.153, 0.108,
  0.100). Ward–Maestro Prop. 28.2 E(γ): 0.175, 0.1875, 0.1667, 0.100. The comparison of interest is
  γ = 1/8 and 1/6, where the prediction is below E(γ) even at dimension 9.

Toy-size factoring here is a correctness and slope check only (AGENT_RULES 8). It is not evidence of
asymptotic cost beyond the stated slopes.
