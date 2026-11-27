# Held-out carry and Euclidean features for the scaled quarter bit

This experiment tests predictors that use integer representatives and carry information. These operations fall outside the fixed sparse-power cancellation theorem. The target is the hidden-pair sign on a uniform unit scale. [Weak prediction reduction](quarter-fourier.md) · [Earlier adaptive exponent screen](adaptive-exponent-screen.md)

## Design

For each of 36 and 64 bits, generate 24 balanced semiprimes with a fixed seed. Use 16 entire moduli for training and eight different moduli for testing, with 4,096 independent unit scales per modulus. The factors enter only the scoring labels. The native feature code takes $N$ and the public scale.

The 106 features use the scale, its modular inverse, its square and cube, and three modular powers. Each representative contributes normalized size, distance to zero, bit count, square-root scale, quarter sign, two shifted Jacobi symbols and six low bits. Extra features use the integer reduction carries, Euclidean quotients, convergent denominator below $\sqrt N$, and public bits of the modulus. Fit a standardized logistic regression and a fixed histogram gradient-boosted model separately in each tier. The planted control appends one feature with a 55% correct label, independently generated for each row. The constant positive predictor measures sampling imbalance. There is no held-out tuning.

Compute simultaneous 95% half-widths across the ten tier-by-model estimates using Hoeffding's inequality conditional on the sampled moduli and fitted models. The intervals quantify uncertainty from sampled scales on those moduli. They do not bound generalization to every balanced modulus.

## Results

Held-out accuracy in percent, using 32,768 scales per tier.

| Bits | Constant | Logistic | Boosted | Logistic with planted feature | Boosted with planted feature |
|---:|---:|---:|---:|---:|---:|
| 36 | 49.481 | 49.799 | 49.649 | 54.868 | 54.233 |
| 64 | 50.415 | 49.991 | 49.475 | 54.800 | 54.404 |

The simultaneous 95% half-width is 0.956 percentage points. The native feature models showed no clear advantage at this resolution. Both models detected the planted signal. The result leaves smaller effects, different feature combinations, and modulus families outside the sampled distribution open.

Run `python carry_feature_screen.py` with NumPy and scikit-learn. This run used NumPy 2.3.5 and scikit-learn 1.8.0. The raw counts give the test result on every held-out modulus. [Raw results](carry-feature-results.json)
