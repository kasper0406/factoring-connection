# Classical factoring research checkpoint

The objective is a classical factoring algorithm with polynomial bit complexity, or progress toward a general algorithm below the number field sieve scale. The work here establishes stronger oracle reductions, a faster evaluator for one fixed scalar, and barriers for specific representations. The native polynomial-bit-time evaluation step remains open. [Original November 26 roadmap](https://github.com/kasper0406/factoring-connection/blob/9514cbfa013cd9ae3f030044db1dc3fb6bb4a9ed/notes/research_roadmap_nov26.md)

## Priorities

1. **General input: evaluate one clipped root count.** Use the parity of mixed roots of $u^2=1\pmod N$ in a translated window of length $(N-1)/2$. It suffices for deterministic bisection on any odd modulus with at least two distinct prime divisors. The same task is realized by signed disk counts modulo 4 at one cutoff and one discriminant. The forms are powers of a public class of known order, and query exponents take polynomial bit time to construct. The simple one-row formulation is the smaller test of a new evaluator. The class-group formulation supplies an alternative if the evaluator can exploit the order structure. [Fixed window](fixed-window-split.md) · [Common discriminant](fixed-discriminant-split.md) · [Cyclic class encoding](cyclic-class-queries.md)
2. **Balanced input: retain the fixed absolute moment as a control.** The implemented hybrid uses roots modulo the public discriminant and tangent batching. It computes the exact scalar with expected $N^{2/5+o(1)}$ bit cost. The 36-bit control has trace 69,529,082 and divisor 254,119. The bound concerns the scalar and remains larger than the complexity of established general factoring algorithms. Further geometric tuning is a reference-evaluator track. [Progression evaluator](discriminant-progressions.md) · [Tangent hybrid](tangent-progressions.md)
3. **Test partial information before requiring exact output.** The quarter-window reduction splits inputs with exactly two distinct prime divisors without adaptive queries. A lattice finishes from a bounded interval or progression, so a fixed-precision evaluator can target just enough bits to narrow a root to $N^{1/2-\varepsilon}$ candidates. Uniformly chosen windows of this size hit mixed roots with exponentially small probability. A candidate must explain how it obtains the partial coordinate from the public input. [Quarter queries](nonadaptive-quarter.md) · [Small-root decoder](bounded-quadratic-roots.md)

A direct lattice-block evaluator now computes central shifted-window parity in $\sqrt N\operatorname{poly}(\log N)$ bit operations. It improves the direct reference exponent for this scalar and still covers a power-of-$N$ number of blocks. [Lattice block evaluator](lattice-window.md)

For any number of distinct odd prime divisors, nonadaptive affine isolation takes $O(\log N)$ weighted queries, or $O((\log N)^3)$ unweighted modulo-4 queries, with constant success probability. These are alternatives when a candidate evaluator can process a batch more cheaply than adaptive queries. [Affine isolation](affine-pair-bound.md) · [Unweighted reconstruction](affine-parity.md)

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

## Reproduce the principal results

Run these commands in this directory. The lattice and tangent programs use SymPy. The direct disk enumerators are controls with power-of-$N$ running time.

```sh
python fixed_discriminant_split.py
python cyclic_class_queries.py
python affine_parity.py
python bounded_quadratic_roots.py
python tangent_progressions.py 38942974393
```
