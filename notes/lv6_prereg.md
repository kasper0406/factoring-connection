# lv6 pre-registration: E1, the LV6 post-processing pipeline on toy non-monic fields

Agent `lv6`, phase 9, 2026-10-19. Written and committed BEFORE any run on the registered seeds.
Calibration runs used seeds 100 and 101 only (not registered; disclosed).

## What is tested
The toy checks the *algebraic* claims of LV6 (notes/lv6_results.md), not any asymptotic or cost claim:
- (C1) the column set suffices: sign, rational primes ≤ B, ONE parity column (number of atoms, which squares
  both g₁ and f_d), first-degree algebraic primes (r, s) with r ∤ D = disc(f)·f_d only (bad primes get NO
  column), plus t quadratic characters at random first-degree primes of norm ≈ 2^128;
- (C2) Algorithm SQ (Hensel lifting at one unramified prime p > B, enumeration of the 2^{t_p − 1} sign
  classes, exact check η² = Γ) decides squareness correctly and is complete under the coefficient bound A;
- (C3) the explicit ψ-construction of Y (with f_d^{−E/2}, ψ(γ⁻)^{−1}, ψ(g′(β))^{−1}) gives X² ≡ Y² (mod N);
- (C4) the residual group detected by characters on the parity-only kernel is nontrivial (so characters are
  needed), and after t = 40 characters the full kernel consists of squares.

## Setup (frozen code: src/lv6/{numfield,sqrt,chars,pipeline,e1}.py; config experiments/configs/lv6_e1.json)
- Cells (d, bits of N): (3, 40), (3, 48), (5, 36), (5, 40); seeds 0, 1, 2, 3 → 16 instances. N balanced
  semiprimes from a seeded generator in the harness.
- Field: design (a″_k) toy: (g₀, g₁) random coprime in [M/2, M], M = ⌊N^{1/(d+1)}⌋, R with coefficients in
  [−2^h, 2^h] (h = 6 for d = 3, h = 2 for d = 5), f = F̂ + (g₁x − g₀)R (src/pcest/geometry.field).
  Fields with f_d = 0, disc f = 0, gcd(f_d g₁ disc f, N) ≠ 1, or f reducible are resampled (events logged).
- Relations: coprime (a, b), 1 ≤ b ≤ A, |a| ≤ A (A = 1200 for d = 3, 500 for d = 5), with G = g₁a − g₀b and
  F(a, b) both 8192-smooth (exact), up to 6000 relations. Rows: disjoint groups of 1..gmax relations with
  random signs e = ±1 (gmax = 3 for d = 3, 1 for d = 5), mimicking u/v Dixon-descent rows.
- Characters: t = 40, first-degree primes sampled uniformly with norm in (2^127, 2^128].
  **Disclosed deviation from the theorem:** the theorem takes log x ≍ C_d (log N)² and t = O_d(log N) from
  the explicit bound on |H|; the toy uses fixed 128-bit primes and t = 40.
- Tests per instance: the first K = 6 basis vectors of the full left kernel (parity + characters) and of the
  parity-only left kernel; then the LV6 combination search over the 2^4 − 1 nonzero combinations of the first
  m′ = 4 full-kernel vectors (first square in canonical order).
- Non-square certificates: residue non-square at p (Gamma not a square in some O_K/𝔓, valid since p is
  unramified and prime to the index), or a character witness among ≤ 200 random first-degree primes of
  40 bits.

## Predictions (frozen)
- **P1 (C3, exact).** Every certified square (η² = Γ exactly) yields X² ≡ Y² (mod N): 0 violations, and no
  "ψ(g′(β)) non-unit" events.
- **P2 (C2, exact).** Every "not a square" output of Algorithm SQ is certified: 0 UNCERTIFIED.
- **P3 (C1/C4).** Eligible instances: ker_full_dim ≥ 4. Holds if ≥ 15/16 instances are eligible and in all but
  at most one eligible instance, all 6 tested full-kernel basis vectors are squares.
- **P4 (LV6 combination search).** In every eligible instance the search finds a square.
- **P5 (C4).** The character rank on the parity-only kernel is ≥ 1 in ≥ 15/16 instances, and the fraction of the
  6 tested parity-only kernel vectors that are squares is ≤ 1/2 in ≥ 14/16 instances.
- **P6 (reported, no prediction).** Fraction of nontrivial gcd(X − Y, N) among squares (plain NFS-type
  congruences, no hidden roots: Theorem DD is not tested here); log₂‖η‖/log₂A; the number t_p of prime factors
  of g mod p and sign tries; bad-prime incidences; log₂|H| bound.

Commands:
```
.venv/bin/python -m src.lv6.e1 experiments/configs/lv6_e1.json experiments/raw/lv6_e1.jsonl
.venv/bin/python -m src.lv6.e1 --analyze experiments/raw/lv6_e1.jsonl > experiments/raw/lv6_e1_summary.json
```
Factor discoveries: gcd(F·G, N) ≠ 1 during relation collection is logged as ACCIDENTAL; gcd(X − Y, N) outcomes
are reported in P6 (classified after the fact by the harness).
