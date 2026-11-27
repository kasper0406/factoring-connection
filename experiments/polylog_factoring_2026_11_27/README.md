# Classical factoring research checkpoint

The objective is a classical factoring algorithm with polynomial bit complexity, or progress toward a general algorithm below the number field sieve scale. The work here establishes stronger oracle reductions, a faster evaluator for one fixed scalar, and barriers for specific representations. The native polynomial-bit-time evaluation step remains open. [Original November 26 roadmap](https://github.com/kasper0406/factoring-connection/blob/9514cbfa013cd9ae3f030044db1dc3fb6bb4a9ed/notes/research_roadmap_nov26.md)

## Priorities

1. **General input: evaluate the one-sided quadratic gap.** One randomized comparison between two public counts on a bounded ten-coordinate polytope with one rank-two quadratic equation supplies each sampled digit. The reduction includes repeated factors and complete recursive factorization. The independent audit charges $O(n^4)$ comparisons per constant-success round under the conservative sampler, where $n=\lceil\log_2(N+1)\rceil$. The missing step is a polynomial-bit evaluator resolving an additive gap of order $N$ against an ambient count of order $N^3$. [Count reduction](single-quadratic-count.md) · [Independent audit](single-quadratic-audit.md)
2. **General exact control: keep the unscaled fixed window.** Deterministic bisection uses the parity of mixed roots of $u^2=1\pmod N$ in translated half-windows. Signed disk counts modulo 4 realize the same queries at one cutoff and discriminant. Exact evaluation for arbitrary affine images is #P-hard even when the modulus factorization is supplied, so affine isolation is a reduction tool rather than a plausible generic evaluator. The hardness proof leaves the unscaled and doubling-only query families open. [Fixed window](fixed-window-split.md) · [Common discriminant](fixed-discriminant-split.md) · [Affine hardness](affine-root-hardness.md)
3. **Balanced input: retain the fixed absolute moment as a control.** Tangent batching computes the exact scalar in expected $N^{2/5+o(1)}$ bit cost. Direct floor reciprocity still takes $\Theta(N^{1/2}\Delta^{1/4})$ terms in the generic branch. Theta, Hecke and finite Poisson transforms preserve the hidden slopes or the high coefficient. These are verified boundary results for future global transforms. [Tangent hybrid](tangent-progressions.md) · [Floor reciprocity](fixed-cap-carries.md) · [Theta and Poisson audit](theta-hecke-transform.md)

A direct lattice-block evaluator now computes central shifted-window parity in $\sqrt N\operatorname{poly}(\log N)$ bit operations. It improves the direct reference exponent for this scalar and still covers a power-of-$N$ number of blocks. [Lattice block evaluator](lattice-window.md)

For any number of distinct odd prime divisors, nonadaptive affine isolation takes $O(\log N)$ weighted queries, or $O((\log N)^3)$ unweighted modulo-4 queries, with constant success probability. These are alternatives when a candidate evaluator can process a batch more cheaply than adaptive queries. [Affine isolation](affine-pair-bound.md) · [Unweighted reconstruction](affine-parity.md)

The exact arbitrary-affine window problem should not be used as the native evaluator target. A parsimonious reduction from counting Subset Sum makes it #P-hard, and its parity is $\oplus$P-hard, even for a square-free modulus whose complete factorization is supplied. The construction scales roots of unity before applying one cyclic half-window. It does not cover the deterministic unscaled family, powers-of-two scaling alone, or the one-sided gap comparison. [Proof and exhaustive checker](affine-root-hardness.md)

On inputs with two distinct prime divisors, an inverse-polynomial prediction advantage on the scaled quarter-window bit also suffices. Fourier localization replaces exact bit reconstruction in that setting. [Weak prediction reduction](quarter-fourier.md)

## Approximate counts and analytic targets

The general odd-input reduction also works with a one-bit gap test on ordinary quadratic counts in a fixed-dimensional domain. The inverse phase modulo $N^2$ gives exact periods for all nonstationary blocks. A lower-input interval isolates a mixed root, and a second affine interval keeps the nonstationary variation small. On a good trial, the gap bit has a constant-size Fourier coefficient at the hidden frequency. Repeated prime factors are included. [Input isolation](square-input-isolation.md) · [Affine interval filter](affine-interval-variance.md) · [One-bit comparison](count-comparison.md) · [Single quadratic count](single-quadratic-count.md)

For inputs with exactly two prime divisors, a constant-precision normalized Salié or square-modulus Kloosterman evaluation would also suffice. One signed phase-interval sum per scale preserves the hidden frequency. Single-term sampling has high variance both before and after grouping by phase, so this route needs aggregation across distinct phase values. [Salié sum](salie-quarter.md) · [Square-modulus sum](square-kloosterman.md) · [Signed interval](salie-phase-bins.md) · [Grouped variance](salie-grouped-variance.md)

The Fourier decoder reuses one sample pool across its interval search. It recovered planted noisy signals on inputs up to 36 bits. The public predictor tests found no clear signal for the specified adaptive exponent rules or for the carry and Euclidean feature models on held-out moduli. [Decoder and calibration](quarter-fourier.md) · [Exponent screen](adaptive-exponent-screen.md) · [Held-out feature screen](carry-feature-screen.md)

## Acceptance test for an evaluator

Use the unfactored input and public window parameters. Compare the exact residue with a direct row scan on both sides of a transition. For the first check, use the two semiprime controls 11,009 and 38,942,974,393, then add repeated-factor inputs such as 315 and 1,575. Include queries with no mixed roots, queries containing multiple mixed roots, and shifts for which the public-root or axis correction changes the raw count. Charge preprocessing, empty blocks, endpoint corrections and any auxiliary modulus construction. After the finite check, give a bound in the bit length for the entire split, including every query the decoder needs.

Full-period totals, local CRT counts taken independently, and representations supplied after finding the hidden root are correctness controls. A successful evaluator must account for which two local residues occur in the same clipped integer interval. On the class formulation, require counts only at the powers actually queried by the reduction. A full-orbit evaluator solves a stronger task.

## Restricted approaches already checked

| Approach | Established restriction |
|---|---|
| Explicit separation of the cutoff into functions of the two prime residues | The unfiltered row can have rank $\Omega(\sqrt N)$. The mixed filtered table can still have rank $\Theta(N^{1/8})$, even when one direction suffices for a closed-form sum. [Full boundary](fixed-cutoff-crt-rank.md) · [Filtered table](filtered-crt-rank.md) |
| Stationary linear state updates | The shift-indexed sequence requires $\Omega(\sqrt N)$ dimension on balanced semiprimes. The full odd-class orbit needs order greater than $L/4$. Read-once binary matrix products can also have exponential width on general inputs. [Shift](window-threshold.md) · [Class power](class-exponent-period.md) · [Binary width](binary-prefix-rank.md) |
| Uniform-input modular samplers | Low-degree rational maps have small success under the square decoder. On independent CRT outputs a pairwise gcd captures rare inverse-polynomial success. Equality branches can be instrumented with gcd probes. These arguments leave integer-representative operations open. [Rational maps](uniform-rational-sampling.md) · [Power maps](power-sampling.md) · [Independent outputs](independent-sampling.md) · [Equality branches](arithmetic-branching.md) |
| Universal shift certificates for mixed roots | Removing the public roots leaves the low-degree ideal of vanishing certificates unchanged. The argument is for every integer lift and leaves interval-specific certificates open. [Proof](mixed-root-certificates.md) |
| Fixed finite families of unweighted cyclic binomial folds | Infinitely many balanced semiprimes make all individual folds units although the exact middle coefficient splits. The proof uses bounded gaps in a fixed progression and does not address growing width families or combined outputs. [Proof](fixed-cyclic-folds.md) |
| Exact arbitrary affine root windows | Counting the selected roots is #P-hard and parity is $\oplus$P-hard even with the modulus factorization supplied. The reduction does not cover unscaled or doubling-only windows. [Proof](affine-root-hardness.md) |
| Genus characters and fixed character quotients on the cyclic class family | Genus characters are constant on every odd power. Finite quotient collisions occur before the queried cutoff on three controls. The collision experiment is not an asymptotic lower bound. [Audit](cyclic-class-evaluator-audit.md) |
| Literal floor-sum reciprocity for the fixed-cap moment | Swapping the summations retains a square boundary and the longer axis. The direct exact sum has $\Theta(N^{1/2}\Delta^{1/4})$ terms on the generic branch. [Identity and checker](fixed-cap-carries.md) |

## Reproduce the principal results

Run these commands in this directory. The lattice and tangent programs use SymPy. The direct disk enumerators are controls with power-of-$N$ running time.

```sh
python fixed_discriminant_split.py
python cyclic_class_queries.py
python affine_parity.py
python bounded_quadratic_roots.py
python tangent_progressions.py 38942974393
python check_affine_root_hardness.py
python check_single_quadratic_audit.py
python check_cyclic_character_quotients.py
python fixed_cap_carries.py
python theta_hecke_transform.py
```
