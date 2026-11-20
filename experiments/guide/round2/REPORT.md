# Round 2 — do smooth parents seed smoother children?

**The matched descent experiment found no reliable inheritance beyond child norm size and local divisibility.** It compared 320 smooth-parent/rough-parent lattice pairs across 20 fresh held-out semiprimes. Each algebraic special-q exceeded B=262,139 and was removed once. A smooth parent had both remaining sides B-smooth. Within each modulus, the arms shared the same polynomial and small-prime roots. Pairs were matched on q, reduced-basis geometry and the norm envelope before scoring any children. The median absolute difference in log q was 0.031.

| Local cutoff y | Smooth-parent hits / uncertain children | Rough-parent hits / uncertain children | Raw difference, percentage points | B1-adjusted difference, percentage points [99% CI] |
|---:|---:|---:|---:|---:|
| 2,000 | 236 / 93,603 | 250 / 93,609 | −0.0243 | −0.0169 [−0.1006, +0.0621] |
| 32,768 | 15 / 47,177 | 11 / 46,718 | +0.0286 | +0.0294 [−0.0294, +0.1096] |

The raw differences weight matched pairs and moduli equally. Before restricting to uncertain children, smooth parents produced 734 hits and controls 659, each over 112,620 children. The apparent advantage on that broader population also had a 99% interval spanning zero. The twofold parent plant was detectable on algebraic-side outcomes at both cutoffs. On joint outcomes, a supplementary fourfold parent plant passed the held-out gain test at both cutoffs. [Parent estimates for both sides and W](results/A_parents.csv) · [Match balance](results/match_balance.json)

## Eight feature families

The run used 1,077,330 distinct target coordinates, cubic and quartic polynomials on 48 fresh 60–110-bit semiprimes, and 447,884 final test targets. Round one tested coarse geometry with q=31. Round two used genuine parent-generated algebraic q lattices, all eight specified feature families, frozen baseline offsets and train-only spline residualization. The primary stratum required both stripped cofactors to lie above the deterministic smooth/rough boundaries. It contained 758 joint hits at y=2,000 and 40 at y=32,768. [Preregistration](PREREGISTRATION.md) · [Frozen validation screen](results/validation_lock.json)

Gains below are **joint log-loss reductions in units of 10⁻⁵ nats per candidate**, with 99% two-stage modulus/block bootstrap intervals. Positive favors the feature family. The last two columns give the percentage of 100 synthetic test-label draws detecting a twofold odds contrast on that family's own feature (Monte Carlo SE at most five percentage points). Each pair lists y=2,000 then 32,768. This measures held-out gain sensitivity with the trained model fixed, before the validation and microcell gates.

| Family | Gain at 2,000 [99% CI] | Gain at 32,768 [99% CI] | Joint plant % | Algebraic plant % |
|---|---:|---:|---:|---:|
| G1 coordinate S-units | −2.1 [−7.6, +2.6] | +2.0 [−1.4, +5.6] | 0 / 4 | 98 / 45 |
| G2 continued fractions and roots | −7.9 [−18.7, −0.3] | +0.2 [−9.2, +11.5] | 0 / 0 | 100 / 97 |
| G3 embeddings, chambers and covariants | −4.1 [−19.6, +6.0] | +3.8 [−2.3, +12.3] | 0 / 0 | 0 / 0 |
| G4 shifts and Farey relatives | −7.0 [−18.5, +1.6] | +1.5 [−11.0, +18.2] | 0 / 0 | 0 / 0 |
| G5 multiplicative and bilinear relatives | −5.9 [−19.1, +4.7] | +1.7 [−1.4, +9.4] | 0 / 0 | 8 / 0 |
| **G6 parent seed** | **−1.6 [−4.6, +0.8]** | **+1.7 [−2.4, +6.5]** | **25 / 1** | **100 / 95** |
| G7 history, cosets and donor digits | −2.2 [−6.3, +1.3] | −1.0 [−14.6, +11.9] | 0 / 1 | 100 / 0 |
| G8 raw bits, carries and MLP | −21.6 [−44.9, −6.5] | +0.8 [−21.0, +19.4] | 0 / 0 | 100 / 49 |

**Every family missed the adjusted discovery rule.** The 108 comparisons covered eight families and their combination, both cutoffs, joint and single-side smoothness, both continuous valuation fractions W, and residual cross-side dependence. Binary discovery required a lower 99% bound above 0.02 nats per hit, Holm p<0.01 and the conditional size-cell check. No hypothesis survived the Westfall–Young validation screen. For joint smoothness, G1 at y=32,768 had Holm p=0.416 and the other 15 family/cutoff comparisons had p=1.000. Within-block side residual covariance was +0.000039 [−0.000248, +0.000325] and +0.000033 [−0.000201, +0.000257], respectively. [All gains, CIs and adjusted p-values](results/A_families.csv) · [Box and lattice results](results/by_regime.csv) · [Covariance](results/A_covariance.json)

## What the controls establish

B0 used the finite-prime saddle-point approximation divided by the Buchstab rough density, with polynomial root-density correction. B1 used cross-fitted 30-knot cubic splines, isotonic calibration and depth-four boosting. At y=2,000, joint uncertain-stratum log-loss was 0.015450 for B0 and 0.015322 for B1. At y=32,768 it was 0.003205 and 0.003538. One single-class training fold required the analytic fallback. The registered B1 screened local columns inside training folds. An additional comparator offered all 955 or 10,630 inputs simultaneously, including every projective prime-power class and valuation. Refitting the families gave the same adjusted verdict. [Baselines](results/A_baselines.csv) · [Full-column check](results/Afull_families.csv)

The independent audit checked 2,275,404 probe nodes, all 5,645,883 retained relative links and the parent determinants. L1 rejected 211 of 432 proposed polynomial/shift combinations. Deliberately exposing the forbidden-band bit improved joint log-loss by 0.000920 and 0.000426 at the two cutoffs. The buffered within-modulus split quarantined donor labels by angle and q band. Its Bernoulli(B0) and randomized-CRT controls yielded no adjusted discovery. [Arithmetic and leakage audit](results/audit.json) · [L1 fixture](results/L1_positive_fixture.json) · [Regional null](results/Bnull_families.csv)

No native twofold **joint** plant passed the held-out gain rule in its original draw. For the algebraic side, G1, G2, G6, G7 and G8 passed at y=2,000, and G2, G6 and G8 passed at y=32,768. None of the native twofold plants cleared every validation, narrow-cell and test gate. Across the entire final sample, the full 0.02-nat conditional analysis had only four joint-outcome-discordant cells at y=2,000 and zero at y=32,768. The Gaussian plant also missed its overall 80–115% recovery requirement. At θ=0.02, recovery was 71% and 49%. These measurements leave economically meaningful joint effects unresolved, particularly for G3–G5. [Native sensitivity](results/native_sensitivity.csv) · [Stronger plants](results/native_dose_response.csv) · [Full-cell tests](results/all_cell_permutations.csv)

The registered arithmetic norm-splitting control had only 12 and zero matched hits. A supplementary control used 24 more fresh moduli and smaller coefficients, preserving exact polynomial root counts and matching the exact algebraic small-prime part, sign and stripped log-size cell. At y=2,000 forced splitting raised algebraic smoothness from 1.36% to 5.14%. G8 detected it (Holm p=0.000676). G2 and G3 missed it. No real-data family supplied evidence of a gain growing with u or an exponent change. [Arithmetic controls](results/pc2b_power.csv) · [Scaling](results/A_scaling.csv) · [Reproduction and data dictionary](README.md)
