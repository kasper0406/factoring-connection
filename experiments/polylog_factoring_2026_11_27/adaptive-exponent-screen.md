# Adaptive exponent screen for the scaled quarter bit

The fixed-power cancellation proofs leave exponents that depend on the scale open. This screen tests seven such rules on balanced semiprimes. No rule produced a clear signal, and the small positive skew for the exponent equal to the scale did not replicate on a larger independent sample. The experiment tests the specified rules and sampled populations. It does not give an asymptotic bound. [Fixed-power result](quarter-all-powers.md)

## Design

The target is the sign of the centered quarter window for $cu\bmod N$, with $u$ a mixed root. Each candidate applies the same public quarter rule to $c^{e(c)}\bmod N$. The exponent choices are the representative $c$, its modular inverse, $\lfloor N/c\rfloor$, $N\bmod c$, $2^{\text{Euclidean step count}}+1$, $2^{\text{popcount}(c)}+1$, and the largest convergent denominator for $c/N$ at most $\lfloor\sqrt N\rfloor$. Exponent 65,537 is a fixed control. A factor-planted control reports the correct sign with probability 55%.

Use eight fresh seeded balanced semiprimes at each of 24, 36, 48 and 64 bits. Draw 8,192 uniform unit scales independently for each modulus, for 262,144 sampled scales. Use the factors only to construct the scoring label. The geometric baseline has accuracy at most $1/2+o(1)$ on these balanced moduli. A union bound with Hoeffding's inequality gives a simultaneous 95% half-width of 0.745 percentage points over the 36 tier-by-rule estimates, conditional on the sampled moduli.

## Results

Values are signed percentage points above 50% accuracy.

| Bits | Scale | Inverse | Quotient | Remainder | Euclid | Popcount | Convergent | Fixed | Planted |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 24 | 0.319 | 0.333 | 0.072 | 0.101 | 0.586 | 0.111 | -0.459 | 0.160 | 5.377 |
| 36 | 0.227 | -0.377 | 0.363 | -0.136 | -0.156 | 0.430 | -0.505 | -0.255 | 4.875 |
| 48 | 0.134 | -0.089 | -0.026 | 0.024 | 0.040 | 0.031 | -0.089 | -0.076 | 4.913 |
| 64 | 0.305 | -0.269 | -0.204 | -0.102 | -0.153 | -0.197 | -0.031 | -0.079 | 4.765 |

The scale exponent was positive at every initial tier. The independent confirmation tested this rule alone on 32 fresh moduli, half at 48 bits and half at 64, with 16,384 samples each. It scored 262,137 correct out of 524,288, or 49.9987%. The pooled 95% Hoeffding half-width is 0.188 percentage points. Smaller, more concentrated, or input-family-specific effects remain compatible with the experiment.

Run `python adaptive_exponent_screen.py` and `python adaptive_exponent_screen.py --confirm-scale`. The raw counts and seeds are saved with this note. [Screen records](adaptive-exponent-results.json) · [Confirmation records](adaptive-scale-confirmation.json)
