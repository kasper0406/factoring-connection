# Can geometry or search history predict NFS smoothness?

**Neither experiment found a replicated improvement over stripped cofactor size and small-prime information.** At the higher smoothness bound, the same pipeline recovered a planted twofold conditional odds ratio in both polynomial families. At the lower bound even a fourfold plant missed the decision threshold, so that experiment has limited sensitivity.

## Design and arithmetic

The study used 104 distinct balanced semiprimes of 60–110 bits and irreducible base-m cubics and quartics. The 3,610,064 labeled targets comprised 1,814,188 ordinary box points and 1,795,876 special-q lattice points. A further 534,528 prior points supplied history. Each target was primitive. The absolute norms were the homogeneous polynomial value and the rational norm |a−mb|. The rational special-q was 31, inside the local cutoff y=31. A relation required both stripped cofactors to be B-smooth.

The first cohort used 40 moduli and B=4,093. After its sensitivity controls failed, a separately preregistered cohort used 64 fresh moduli and B=262,139. Fit, calibration, model selection and final test used separate N and disjoint eight-row tiles in each coordinate system. The 48 final moduli contributed 1,677,744 targets. Modulus sizes spanned the same range at fit and test. [Original protocol](PREREGISTRATION.md) · [Independent cohort protocol](power_cohort/PREREGISTRATION.md)

PARI factorization and repeated prime-product gcd agreed on all 8,289,184 side labels, including references. A separate audit reconstructed both norms and local factors for every saved row, repeated full PARI factorization on 29,952 sampled sides, and verified all 57,761,024 target–reference links. Each link had nonzero determinant coprime to every prime in (31,B], excluding equal projective residue classes. The higher bound initially left 381 targets on the first fit modulus with fewer than 16 eligible references. Before fitting, the reference pool grew from 1,024 to 8,192 per modulus. Every target and exclusion was retained. [Arithmetic audit](results/integer_audit.json) · [Higher-bound audit](power_cohort/results/integer_audit.json) · [Reference-pool amendment](power_cohort/HISTORY_AMENDMENT.md)

## Held-out prediction

The logistic baseline used piecewise functions of stripped log sizes and exact prime valuations through 31, followed by sigmoid calibration. Boosting and a two-hidden-layer MLP each tested base inputs alone, added coarse geometry, and added geometry plus prior 4- and 16-reference outcome rates. The comparison was paired within region and averaged equally across regions and moduli. A signal required greater than three modulus-cluster SE against both logistic and the same learner using base inputs, replicated in both families, with the registered sign-flip and permutation checks.

| B | Family | Held-out relations / targets | Selected joint log-loss | Gain over logistic, 95% t interval | Gain over matched model | Feature permutation p |
|---:|---|---:|---:|---:|---:|---:|
| 4,093 | Cubic | 82 / 279,624 | 0.00213300 | +0.00000106 [−0.00000925, +0.00001137] | +0.00000108 | 0.74 |
| 4,093 | Quartic | 121 / 279,624 | 0.00284757 | −0.00000782 [−0.00001420, −0.00000145] | +0.00000946 | 0.16 |
| 262,139 | Cubic | 12,691 / 559,248 | 0.09303785 | 0 [0, 0] | 0 | 1.00 |
| 262,139 | Quartic | 15,236 / 559,248 | 0.10778315 | 0 [0, 0] | 0 | 1.00 |

Gains use natural logs per candidate. Both cohorts selected geometry boosting. At B=262,139 calibration assigned zero joint weight to all six flexible fits, making their predictions identical to logistic. The zero intervals describe that identity. Scoring the frozen models without blending also gave worse joint log-loss than logistic for all six fits in both families. The lower-bound selected AUC gains versus logistic were +0.00086 and −0.00235. Only four cubic and two quartic moduli had both outcomes in both regions. Higher-bound joint AUCs were 0.67080 and 0.68159.

The largest apparent single-side gain of the selected model was on higher-bound quartic rational norms, +0.000228 against logistic. Against base-only boosting the gain was −0.0000114. No tested geometry or history variant improved either side by three SE against both baselines in both families. [All lower-bound metrics and intervals](results/family_metrics.csv) · [All higher-bound metrics and intervals](power_cohort/results/family_metrics.csv) · [Higher-bound unblended fits](power_cohort/results/true/unblended_joint.csv)

## Controls

At B=4,093 the original odds-ratio 2.01 plant and the pre-announced 4.06 plant both missed. At B=262,139 the 2.01 plant improved log-loss by 0.00113726 on cubics (3.89 SE) and 0.00139322 on quartics (5.42 SE) against both baselines. The exact one-sided sign-flip probabilities were 0.0000305 and 0.0000153. Restricted-feature permutation p was 0.01 in each family. At the 20% budget, 12,368 planted positives were recovered versus 11,557 under logistic, a 7.02% increase. [Original plant](results/planted/heldout_summary.json) · [Stronger plant](results/planted_strong/heldout_summary.json) · [Higher-bound plant](power_cohort/results/planted/heldout_summary.json)

Both registered shuffled-label runs missed the replication rule. At the higher bound the shuffle preserved the event rate of each N/region. It passed the quartic comparison, and feature permutation gave p=0.01 in both families. Re-randomizing target labels within the same cells showed an expected quartic gain of 0.000806 against logistic, versus 0.000770 observed. The cubic within-cell excess gave p=0.020, quartic p=0.772. A supplementary global shuffle erased those rate differences and gave exactly zero gain in both families. The restricted-feature permutation alone is therefore insufficient evidence for candidate-specific prediction. [Stratified shuffle](power_cohort/results/shuffled/heldout_summary.json) · [Within-cell randomization](power_cohort/results/shuffled/within_cell_placebo.json) · [Global shuffle](power_cohort/results/global_shuffle/heldout_summary.json)

## Yield and work

At B=4,093, 111,840 confirmations yielded 136 real relations for the selected model versus 137 for each baseline. The one-thread replay measured 62.0 versus 85.1 relations/s for selected versus logistic, including all-candidate norm stripping, scoring, ranking and exact confirmation. At B=262,139 all three policies returned 11,536 relations in 223,680 checks. The replay measured 2,368 versus 2,392 relations/s. An implementation can skip the zero-weight model.

When reranking separately within the four registered u bins, selected/logistic yields at B=4,093 were 19/19, 73/73, 18/19 and 2/2. At the higher bound they were 6,985/6,985, 3,485/3,485 and 67/67. The highest bin was empty. These observations give no evidence of gain growing with u or of an exponent change. The tested features and finite ranges leave the general conditional-smoothness hypothesis open. [Lower-bound ledger](results/pooled_yield.csv) · [Higher-bound ledger](power_cohort/results/pooled_yield.csv) · [Reproduction and data dictionary](README.md)

The run took 15.9 minutes from first preregistration through final numerical checks. Instrumented stages used 17.5 CPU-minutes. [Timing ledger](results/run_accounting.json)
