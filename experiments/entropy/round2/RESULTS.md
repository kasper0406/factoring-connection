# Entropy round 2

**Zero tests met the preregistered polynomial-lead rule.** All 20 locked test–target pairs passed the split-clean association test at 64 bits. All five taken forward passed Holm correction at both 80 and 96 bits. The five contain three distinct comparisons. Their fitted exponential rates are close to one factor of two per four additional product bits.

## Confirmation

The table gives signed centered advantage in **parts per million** after the registered split exclusions. Centering removes the agreement expected from predictor and target frequencies within public residue and size strata. The last column is twice the signed-data log likelihood of the polynomial fit minus twice that of the free exponential fit.

| Frozen comparison using the first 64 CF steps of √(7N) | 64 bits | 80 bits | 96 bits | Fitted β in 2⁻ᵝⁿ | 2Δlog L |
|---|---:|---:|---:|---:|---:|
| Ranks 1–2, parity of the even-index certificate count vs. the factor's quadratic character of 2 | 450.81 | 29.26 | 1.885 | 0.2463 | −6,416 |
| Ranks 3–4, low-byte parity of packed counts vs. the same character | 448.66 | 31.55 | 1.709 | 0.2506 | −3,875 |
| Rank 5, Jacobi readout of the packed counts vs. the smaller factor being 1 mod 3 | 136.03 | 8.67 | 0.493 | 0.2499 | −6,086 |

The first four locks apply when N ≡ 1 mod 8. On that fiber both factors have the same character of 2, making the smaller-factor and random-factor labels identical. The primary uniform-sample analysis gave likelihood ratios of −1,626, −1,000 and −1,995 in the same order. The augmented analysis incorporates an outcome-blind, fixed-size top-up on the eligible population. It overlaps the primary cohort. Prime-cluster errors, dropping the 40-bit anchor, and allowing polynomial powers through 32 all preserved the preference for exponential decay.

At 96 bits, rank 1 had 959 positive readouts in 224,971,804 clean rows. Its exact fixed-margin association p was 2.55 × 10⁻⁴⁴. The 95% centered-effect interval under a common conditional odds ratio was [1.634, 2.125] ppm. The distribution-free centered interval was [−9.164, 9.164] ppm. Its signed **raw** agreement advantage was −54.22 ppm with exact 95% Clopper–Pearson interval [−119.56, 11.12] ppm. These intervals answer different questions.

For this readout the positive rate fell from 1,117 ppm at 64 bits to 4.26 ppm at 96 bits. A qualifying CF denominator has squarefree core dividing 210. The convergent identity A² − 7NB² = (−1)ⁱQᵢ forces its signed core to be a square modulo both factors on clean rows. All 223 rank-1 positive rows in the primary 96-bit sample passed independent certificate checks. Including flagged rows raised the augmented rank-1 centered effect to 2.086 ppm. Of 224,998,647 eligible rows, 26,843 had at least one flag, including 562 with a proper gate gcd. The per-reason counts overlap.

## Power and controls

At 96 bits the uniform primary cohort recovered the planted 1/n² agreement leak at +78.86 ppm, exact 95% CI [12.78, 144.93], one-sided p = 0.00967. On the augmented N ≡ 1 mod 8 cohort the plant measured +157.19 ppm [91.85, 222.53], against a planted 108.51 ppm. Power at the single-test 1% cutoff was 82.34%, or 64.68% at the first five-test Holm cutoff. The 1/n plant measured +0.010399 [0.010333, 0.010464] with power effectively 100%. The weak plant's signed estimates alone did not resolve its decay law. The independent low-bit sentinel gave p = 0.430. The within-pair shuffles had zero discoveries and the public-floor sanity checks had zero centered signal.

The descriptive MI audit found substantial sparse-cell bias. For rank 1 at 96 bits, stratified plug-in MI was 4.249 × 10⁻⁶ bits per eligible row. Its exact fixed-margin null expectation was 4.254 × 10⁻⁶. Miller–Madow alone left 4.242 × 10⁻⁶. The excess after subtracting the null expectation was −4.86 × 10⁻⁹ bits. The signed, pooled association test can detect a shared direction across strata that individually contain only one rare positive.

Restricting the synthetic control to rank 1's exact clean rows shows the effect of the top-up directly. At 96 bits its single-test power rose from 23.67% to 82.34%, and its observed p moved from 0.135 to 1.20 × 10⁻⁶. At 64 bits the corresponding augmented power was 78.42% and the realized p was 0.0312, above the 1% threshold. matched_plants.tsv records each lock's actual rows, outcomes and powers.

## Scope and verification

The run covered all seven grammar families. Of 53 primitive variants, 51 met the measured 50 µs cutoff at 64 bits. Their readouts and compositions yielded 20,000 truth tables distinct up to complementation on 720,000 discovery products. Screening 66 smaller- and random-factor targets across four screens gave 5,280,000 comparisons. Bonferroni sign tests and prime-cluster sensitivity promoted 30 pairs to independent validation. All 20 locks came from multiplier CFs. No PNT or GL target locked. The reduced 40-bit universe and split filters left only 766 clean discovery rows there.

The data comprise 599,830,000 sampled products. Discovery and validation used unique products and mutually disjoint prime pools, also disjoint from confirmation. Confirmation used independent draws with replacement. Among 168 million 64-bit, 40 million 80-bit, and 390 million 96-bit draws, the product-repeat counts were 88, zero and zero. Before any confirmation estimates were inspected, a second registration added 128 million uniform 64-bit draws and 170 million 96-bit draws on N ≡ 1 mod 8. This delivered at least 40.64 million clean rows per 64-bit lock and at least 219.97 million per 96-bit lock.

An independent numerical audit caught a round-to-nearest conversion in the initial PNT target encoding. We corrected it to the specified floor, replayed every discovery and validation comparison, and verified the original promotion and lock file hashes. Independent checks passed for 48 PNT quadratures, 10,240 regenerated high-size inputs and predictions, and all 456 agreement intervals. The runbook documents the implemented proxies and omissions, split-budget visibility, model bounds and replay commands.

The run finished within one hour of the three-hour wall budget. The complete numerical records are in the [confirmation estimates](results/confirmation.tsv), [fixed-margin tests and intervals](results/conditional_exact.tsv), [signed likelihood fits](results/scaling.tsv), [split counts](results/split_flags.tsv), [matched plant measurements](results/matched_plants.tsv) and [MI audit](results/mi_null.tsv). The [figure](results/decay_and_power.png) shows the measured decay and recovery of the weak plant. See the [registration](PREREGISTRATION.txt), [power addendum](POWER_ADDENDUM.txt) and [replay guide](RUNBOOK.txt) for the design and implementation.
