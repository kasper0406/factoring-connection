# periodfind pre-registration (agent `periodfind`, 2026-10-20)

Frozen before any experiment below is run. Code: `src/periodfind/*`. Raw output:
`experiments/raw/periodfind_e{1,2,3,4}.jsonl`.

## Main question (orchestrator focus)

For which cheaply computable g: (Z/N)^* → C does F_g(t) = g(a^t mod N) put a large fraction of its
L² mass on frequencies at (or near) multiples of 1/r, with r = ord_N(a)? If such a g exists, sparse
Fourier learning (KM / Akavia SFT / GMS) returns a sample j/r from polylog evaluations.

## Analytic predictions (to be proved in notes/periodfind_results.md, not by these experiments)

- **T (transfer identity).** Let H = ⟨a⟩, r = |H|, i = φ(N)/r, and ĝ(χ) = E_{y∈(Z/N)^*} g(y)χ̄(y).
  On Z/r, F̂_g(j) = E_t F_g(t)e(−jt/r) = Σ_{χ: χ(a) = e(j/r)} ĝ(χ). Hence |F̂_g(j)|² ≤ i·Σ_{χ(a)=e(j/r)}|ĝ(χ)|².
  A heavy peak of order d forces g to correlate with a multiplicative character of order ≥ d.
- **W (additive/Weil-class flatness).** If g(y) = Σ_h b_h ψ_h(y) e(h y/N), with h ∈ (Z/N)^*, ψ_h ∈ {1, (·/N)}
  and Σ|b_h| ≤ B, then every j whose η_j lies outside {1, (·/N)|_H} has |F̂_g(j)| ≤ B√N/r.
- **GR (generic ring).** If g is computed from y by ring operations mod N with zero tests (any size),
  the mass on frequencies j with (r/gcd(r_p, r_q)) ∤ j is ≤ 8·ε_split. Here ε_split = Pr_t[some zero test
  on the run splits N].
- **Comb.** For g = [y = c], every coefficient has relative mass exactly 1/r.

## E1: spectrum census (exact orbit enumeration, float FFT with a stated error bound)

- Instances: `src/common/instances.py` balanced semiprimes, bits n ∈ {22, 26}, seeds 0..9.
  For each, 3 random bases a (seeded), keeping a with r ≥ 2^{n−6}. (HARNESS: r is computed from the
  factors.)
- Observables (panel P1–P12):
  - P1 comb [y = c];
  - P2 additive e(hy/N);
  - P3 interval [y < N/2] (±1-centred);
  - P4 LSB (−1)^y;
  - P5 Jacobi (y/N);
  - P6 shifted Jacobi ((y+1)/N);
  - P7 Jacobi × interval;
  - P8 LSB(y² mod N);
  - P9 LSB(y^{-1} mod N);
  - P10 B-smooth indicator of the integer y (B = 50), centred;
  - P11 OA test [y^M = 1] with M = lcm(1..B₁), B₁ = 30;
  - P12 pseudorandom ±1 hash of y (baseline).
- Metric: M_peak = max over j with η_j ∉ {1, Jacobi} of |F̂(j)|² / ‖F‖², and H_mass = the total relative
  mass on j of order ord(η_j) ≥ √r.
- Predictions:
  - P1: M_peak = 1/r exactly.
  - Every observable except P11 (and except Jacobi at its exceptional frequency): M_peak ≤ 30·ln(r)/r.
    This is a random-function baseline with slack, and fails nowhere. For the Weil-type panel members
    (P2, P6), M_peak also respects the bound of prediction W.
  - P11: heavy (M_peak ≥ 1/100) exactly when r/gcd(M, r) ≤ 100. Otherwise it obeys the 30·ln(r)/r bound.
- Kill criterion for the barrier picture: any non-character observable with M_peak ≥ 1/100 not explained by
  OA (P11) refutes "only characters concentrate" for that class.

## E1b: transfer identity check (diagnostic)

N ≤ 2000, all characters enumerated from the factorisation. Check T to within 1e-9 for every j and every
panel observable. Prediction: 0 violations.

## E2: SFT runs (operational part sees only N, a and query access)

- Algorithm: KM-style low-bit bucketing on Z_Q, Q = 2^m ≥ r², with an unbiased bucket-weight estimator of
  2 queries per sample. Threshold τ ∈ {1/4, 1/16} (relative to the estimated ‖F‖²). Hoeffding sample sizes;
  query counts recorded.
- Functions: positive control g = e(log_a(y)/r) (HARNESS-built pure tone); Jacobi; additive; interval; comb.
- Predictions:
  - The positive control recovers r by continued fractions in ≥ 9/10 runs.
  - Jacobi outputs only the frequency ≈ Q/2 (and 0).
  - Additive, interval and comb output no frequency other than ≈ 0 (comb: nothing, since the estimated
    ‖F‖² is ≈ 0 from fewer than r samples).
  - Queries ≤ 10^7 per run.

## E3: generic-ring theorem regression (diagnostic)

Random decision trees of depth ≤ 3 with zero tests y^m − c. Here m is drawn from {random, multiples of
r_p, multiples of r_q, multiples of r} (factor-aware HARNESS choices, to make tests non-constant). Check
"visible mass ≤ 8 ε_split" on every tree. Prediction: 0 violations (it is a theorem; this checks the code).

## E4: single-peak estimation from random samples ("only part of the pattern")

For F = additive phase e(h a^t/N), estimate F̂(k) at the true strongest peak from m uniform samples,
m/r ∈ {1/16, 1/4, 1, 4}. Detect the peak if the |estimate|² exceeds the 99% null quantile at a random
non-peak k. Prediction: detection rate ≤ 0.2 at m ≤ r/4, and rising only once m ≳ r.
