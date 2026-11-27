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

## Degree-eight benchmark and one reciprocal quotient range

The smooth upper-digit moment uses $P_8(t)=B_8(t)/B_8(0)=1-20h^2-40h^3-30h^4$, where $h=t(1-t)$. Average $P_8(z/N)$ over the lifts of each unit lower digit $a$ in the isolation window and subtract the two public stationary blocks. For $d=\gcd(a^2-c_0^2,N)$ the upper digits run through $s+d\mathbb Z$ with multiplicity $d$. Bernoulli multiplication gives the exact block average $(d/N)^8P_8(s/d)$. The remaining stationary block on a singleton trial contributes $P_8(((z_0+2uj)\bmod N)/N)$. [Benchmark](bernoulli_moment.py) · [Isolation](affine-pair-bound.md)

The local sum of eighth-power weights at $p^e$ is at most $2/(1-p^{-7})$ by the valuation counts in the input-isolation proof. Write $Q=2^{\omega(N)}$. Subtract the $Q$ stationary weights of one and average over the translated window. The expected supremum of the remaining signal is at most $(wQ/N)/1457$, using $\sum_{p\ge3}p^{-7}\le1/1458$. At the isolation width $wQ/N\le 1/8192$, Markov's inequality leaves a constant probability of a singleton with remainder at most $1/64$. A singleton in a width-one window has zero remainder. The Fourier series of $P_8$ gives $\|P_8-\cos(2\pi\cdot)\|_\infty<1/96$. With a $1/64$ comparison gap about zero, the bit can disagree with the cosine sign only on $|\cos|<1/24$. This union of two arcs occupies at most $N/24+2$ digits. The clean sign has coefficient $2\cos(\pi/(2N))/(N\sin(\pi/N))$ at frequency $2u$. The corrupted coefficient remains above $1/4$ for $N\ge15$. The query uses a fixed number of quadratic counts because a $z^k$ weight can be counted by adding $k$ variables in $[0,z)$. [Local valuation counts](square-input-isolation.md) · [Fourier recovery](quarter-fourier.md)

For the reciprocal benchmark let $M=N^2$, $A=(c_0+Nj)^2\bmod M$, and $I=[L,H)$ be one interval of centered denominators with the same nearest quotient $q$. Inverse reciprocity rewrites
$$e_M(A/x)=e_{|x|}(-\operatorname{sgn}(x)A\overline M)\,e^{2\pi i A/(Mx)}.$$
Omitting the smooth correction over the centered range costs $O(\log N/N)$ in the normalized sum. Within the bin, $r=M-q|x|$ runs in a progression and a second inverse reciprocity step replaces the reciprocal modulus by $|r|$ when $q$ is a unit. Recursing on each remaining modulus individually retains the cost of visiting its denominator. [Reciprocity checks and evaluator](reciprocal_range.py)

The prototype evaluates the whole bin at once in the original coordinate. For each unit lower digit $a\in W$, put $b(a)=1-c_0^2/a^2\bmod N$ and $Y_a=\{y:L\le a+Ny<H\}$. The normalized $j$-signal has coefficient
$$ e_M(a+c_0^2/a)\,\frac1N\sum_{y\in Y_a} e_N(b(a)y) $$
at the distinct frequency $2c_0/a\bmod N$. This gives the exact scan of $w$ lower digits and removes the upper-lift enumeration. For $b\ne0$, the finite geometric sum divided by $N$ has magnitude at most $1/(2|b|_N)$. Each derivative value has at most $Q$ preimages since it specifies one unit square. By orthogonality in $j$, truncating to $|b|_N\le K$ thus has mean squared error at most $wQ/(2NK)$ over uniform translations of $W$ and uniform $j$, at the actual bin endpoints. The stored spectrum serves every query digit. A fixed macroscopic interval retains the stationary coefficient with magnitude $|Y_a|/N=|I|/M+O(1/N)$.

The truncated evaluator generates the retained digits by solving the $2K+1$ public congruences $(1-b)a^2=c_0^2\pmod N$, skipping nonunit coefficients. The same certified windows of size $\Theta(\sqrt N)$ used for mixed-root parity solve these general monic quadratics. Its preprocessing costs $O((K+1)(1+w/\sqrt N)\operatorname{poly}(\log N))$ bit operations, including empty blocks. Each subsequent digit uses only the retained coefficients. Fixed accuracy permits a fixed $K$ at the isolation width. Central windows still require exponentially many blocks in the bit length. [Certified block solver](lattice-window.md)

The degree-eight control recovered a verified divisor on six synthetic singleton windows, including repeated factors. On the $N=1575$ fixture it used 22 progression weights in place of 34,650 direct lift weights per query. The bin diagnostic matched 90 literal sums against both reciprocity steps, checked the spectrum by Parseval, and compared the truncated solver with a full lower-digit scan. At $N=471311$ and width 3682, 688 lattice blocks retained one coefficient, with RMS truncation error $0.004342$. These figures measure work units and approximation error on controls. [Moment benchmark](bernoulli_moment.py) · [Bin benchmark](reciprocal_range.py)

On the 36-bit control $N=38942974393$, the positive $q=2$ bin at start 12391793379 and width 131072 used 96 lattice blocks and retained one coefficient at $K=8$. The lower-digit scan took about 1.09 seconds and the truncated evaluation about 0.10 seconds in a run with FLINT. The normalized RMS truncation error was $5.61\times 10^{-6}$. The fixture uses known factors to place a singleton window. The evaluators take only the public window parameters. SymPy's Python LLL backend failed an internal assertion on one basis in this run. To use the FLINT backend, run the following from this directory.

```sh
uv run --no-project --with sympy --with python-flint==0.9.0 python reciprocal_range.py 38942974393 --start 12391793379 --width 131072
```

## Recovering a factor from the full stationary magnitude

Let $M=N^2=\prod_i q_i$, where $q_i$ are the complete prime-power components, and let $e_i$ be the CRT idempotents modulo $M$. For unit $t$, stationary phase gives the full normalized Kloosterman sum
$$F(t)=K_M(1,t^2)/N=\sum_{u^2=1\bmod M}e_M(2ut)=\prod_i2\cos(4\pi e_i t/M).$$
The product extends the signal to all $t$. A proper $\gcd(t,N)$ already splits the input. At the $1/N$ fraction where $t\equiv0\pmod N$, using zero in the transformed signal changes each Fourier coefficient by at most $O(\log N/N)$. [Stationary-phase identity](square-kloosterman.md)

Use the smooth transform $H(t)=\log(1+|F(t)|)$. On uniform $t\bmod M$, put $C_i=\cos(8\pi e_it/M)$ and $X_i=\frac12\log(2+2C_i)$. The $X_i$ are independent, and $\mathbb E X_i=\log 2/q_i$ by taking absolute values in the product of $1+e_{q_i}(r)$ over $r\bmod q_i$. Their variances are bounded below and their centered fourth moments bounded above by absolute constants. In fact each of $C_i\ge1/2$ and $C_i\le-1/2$ contains a constant fraction of the odd grid. Also $|2\cos(2\pi r/q_i)|\ge2/q_i$, and counting grid points near the two cosine zeros gives $\Pr(X_i<-s)\le(3/2)e^{-s}$ for $s\ge0$. These estimates are uniform in the prime powers.

Fix $i$ and write $S=\sum_{k\ne i}X_k$. If $Y=S-\mathbb ES$ and $\sigma^2=\mathbb EY^2$, independence and the fourth-moment bounds give $\mathbb EY^4\le C\sigma^4$ for an absolute $C$. Hölder gives $\mathbb E|Y|\ge\sigma/\sqrt C$. Since $\mathbb E Y=0$, Cauchy–Schwarz then gives $\Pr(Y>0)\ge1/(4C)$. Hence $\Pr(S\ge0)$ is bounded below. For fixed $S$, $H=\log(1+e^{S+X_i})$ is increasing in $C_i$. When $S\ge0$, its values on the two bands differ by at least $\log((1+\sqrt3)/2)$. The independent-copy formula for covariance now gives $\mathbb E[H C_i]\ge c$ for an absolute $c>0$. This is a real Fourier coefficient of $H$ at $4e_i$, and $\gcd(4e_i,N)=N/\sqrt{q_i}$ is proper whenever the input has at least two distinct prime divisors.

The transform is 1-Lipschitz in the trace magnitude and is $O(\log N)$ in size. A sufficiently small constant RMS error in the normalized magnitude therefore leaves an inverse-polynomial heavy coefficient after rescaling and clipping. The sampled Fourier decoder works directly on uniform $t$ and can verify factors by gcd. Even one-bit magnitude comparisons suffice. Choose a public uniform threshold on a constant-mesh grid for $H$ and test the corresponding threshold $e^h-1$ on $|F|$. The mean answer approximates rescaled $H$, and a constant raw additive gap perturbs it only over an interval of thresholds of constant length. No lower-digit isolation window is needed. Computing the complete normalized magnitude without factors remains open. [Fourier recovery](quarter-fourier.md)

The diagnostic verified the product identity on seven synthetic grids up to 50,625 points and checked 18 literal square-modulus Kloosterman sums. Quantizing the raw trace to a grid of step $1/64$ preserved the predicted nonunit coefficients. Dense Fourier recovery using the quantized values found a proper divisor on all six grids where it was run. Known prime powers generated the controls. [Diagnostic](stationary_log.py)

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
python bernoulli_moment.py
python reciprocal_range.py
python stationary_log.py
python tangent_progressions.py 38942974393
```
