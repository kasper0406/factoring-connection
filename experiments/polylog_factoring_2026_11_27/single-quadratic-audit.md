# Independent audit of the single-quadratic comparison

The randomized reduction is algebraically correct on its stated oracle model. Given a polynomial-bit-time evaluator for the public count threshold, it gives expected polynomial-bit-time complete factorization. It covers repeated odd prime powers and any number of distinct odd prime divisors. It does not give a polynomial-time evaluator for the count.

## Encoding audit

Put $M=N^2$ and make the linear change $v=t-x$. The only quadratic constraint becomes

$$xv=c_s^2+r+Mq.$$

When the complement flag is zero, $r=0$ and this is a clipped modular hyperbola $xv\equiv c_s^2\pmod M$. When the flag is one, $1\le r<M$ counts the complement. For fixed values of the other variables, Euclidean division supplies a unique $(r,q)$. The affine half-interval selector supplies a unique $(\eta,k)$, and the four xor inequalities supply a unique $\lambda$. Thus auxiliary variables introduce no multiplicity.

For each lower representative, lift, phase representative and upper digit, exactly one of the two values of $s$ is complemented. Its ambient population is

$$B=wN^2(N+1)/2.$$

The remaining case has the desired sign in the four-term contrast. This proves $K_\sigma=B+\sigma D_{\rm raw}$. The quadratic form on the independent linear coordinates $(x,t)$ has matrix determinant $-1/4$, hence rank two.

All ten coordinates have $O(\log N)$ bits. In particular, $|k|\le N+1$. The ranges $|x|<N^2$, $0\le t<(N^2+N)/2$ and $c_s<N^2$ give $|q|<8N^2$. The count and its threshold are below a fixed power of $N$, so their binary encodings also have $O(\log N)$ bits.

The independent checker enumerates the free coordinates, reconstructs every auxiliary coordinate, verifies every inequality and compares the result with four separately counted congruence slices. It passed 30 fresh encodings on $N=9,15,21,25,45$, including prime powers and products with repeated factors. It checked 604 accepted modular-hyperbola points. The largest observed ratios were $|k|/(N+1)=0.625$ and $|q|/N^2=1.809$. [Independent check](check_single_quadratic_audit.py)

## Reduction and total work

Let $n=\lceil\log_2 N\rceil$ and retain the fixed constant $C=8192$. Trial division reaches $2^{18}C(n+1)=O(n)$. This has an enormous fixed constant and polynomial bit cost.

There are $O(n)$ public dyadic isolation widths. At the useful width, a trial has an isolated mixed root and small nonstationary remainder with probability at least $1/(16C)$. The random polarity preserves the signal with probability one half. The resulting comparison bit differs from a half-interval bit on at most $1/32$ of its upper digits, leaving a constant Fourier coefficient at twice a nontrivial square root of unity.

The reusable Fourier estimator needs $O(n^3)$ predictor pairs at constant threshold and constant failure probability. This follows from $O(n)$ distance shells, squared bounded-weight concentration and a union bound over the full binary interval tree. Each predictor value is one public quadratic-count threshold. Trying a fixed multiple of $C$ trials at all widths therefore gives a crude $O(Cn^4)$ count-query bound per constant-success round. All frequency candidates are checked with gcd and exact division. Repeating unsuccessful rounds gives expected polynomial time under the evaluator assumption.

Even factors are removed directly. A leaf with one distinct odd prime is recognized by polynomial-time primality and exact perfect-power tests. Every other odd node has a mixed square root of unity and is split recursively. Repeated exponents do not change the local root count, so the proof covers arbitrary positive integers.

## Why fixed dimension does not evaluate the count

Barvinok's fixed-dimensional algorithm applies to lattice points in rational polyhedra and their linear projections. The hyperbola equation above is not a linear projection constraint. [Barvinok and Woods](https://arxiv.org/abs/math/0211146)

After the change $v=t-x$, an exact evaluator must count divisor pairs of $c_s^2+N^2q$ in coupled clipped intervals while $q$ ranges over a power-sized interval. Listing $q$, its divisors, lower input blocks or upper digits loses polynomial bit complexity. The factors of $N$ do not by themselves separate the ordinary interval representatives, so a polynomial-time converse from factoring to this count has not been established.

General binary quadratic feasibility is already NP-complete for equations as restricted as $\alpha x^2+\beta y=\gamma$. This prevents an appeal to fixed dimension and degree alone. It is not an NP-hardness proof for the structured modular-hyperbola family used here. [Manders and Adleman](https://doi.org/10.1016/0022-0000(78)90044-2)

The exact count is a #P function because a coordinate tuple has polynomial length and its constraints are checked in polynomial time. The threshold family is randomized-factoring-hard by the reduction above. No #P-hardness result or factoring-equivalence in both directions is proved here.

The required resolution is also far beyond a generic relative approximation. The useful gap has width $N/2$, while the ambient count is $B=wN^2(N+1)/2$. Its relative scale is $\Theta(1/(wN^2))$, between $N^{-2}$ and $N^{-3}$ on the queried widths. An FPRAS whose cost is polynomial in the reciprocal relative error would therefore take a power of $N$. A candidate needs to preserve the signed cancellation around the public baseline rather than approximate the two large counts independently.

Standard modular-hyperbola discrepancy estimates reach an additive error of order $M^{1/2+o(1)}$ at modulus $M$. Here $M=N^2$, so this is $N^{1+o(1)}$, while the required gap is below $N/2$. Applying a Weil or incomplete-Kloosterman bound to the positive counts separately therefore misses the decision scale. This is an analytic comparison rather than a lower bound on algorithms. [Shparlinski, Theorem 13](https://arxiv.org/abs/1103.2879)

## Verdict

Keep this as the broadest current reduction to one fixed-dimensional public count comparison. The strongest next test is an exact or gap-preserving algorithm for clipped modular-hyperbola counts that couples lower blocks before enumeration. Polyhedral counting, CRT separation and enumeration over the quotient $q$ do not supply that algorithm.

Run the independent check from the repository root.

```sh
python experiments/polylog_factoring_2026_11_27/check_single_quadratic_audit.py
```
