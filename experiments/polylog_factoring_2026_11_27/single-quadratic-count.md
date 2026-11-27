# One quadratic count with a public threshold

The gap oracle in the count comparison reduction can be stated as one ordinary lattice-point count at a public threshold. Its input is a bounded rational polytope in ten coordinates and one integer quadratic equation whose quadratic part has rank two. The reduction needs one count comparison for each sampled upper digit. A polynomial-time evaluator for this decision problem remains open. [Gap comparison and isolation](count-comparison.md)

## Pack the signed contrast by complementation

Let $M=N^2$, $H=(N-1)/2$, and $c_1=c_0+Nj$. Let $b_0$ be the isolation offset, $w$ the width, and $\alpha,\beta$ the affine sign parameters. The input and phase representatives are
$$
 x=s_0-b_0+Ny,\qquad t=t_0+Nz,
 \qquad 0\le s_0<w,\quad 0\le y,t_0<N,\quad 0\le z\le H .
$$
Choosing $y$ over a full period gives each allowed $x\bmod M$ once, including when $s_0-b_0<0$. A selector $\eta\in\{0,1\}$ records the sign class by a linear constraint with a unique integer quotient $k$,
$$
 (H+1)\eta\le \alpha(s_0-b_0)+\beta-Nk
 \le H+(N-1-H)\eta .
$$
Sign class zero has positive sign. An independent selector $s\in\{0,1\}$ chooses the reference or query digit. Let $C_{s,\eta}$ be the count satisfying $xt-x^2\equiv c_s^2\pmod M$ in the corresponding region. The signed contrast before the public stationary correction is
$$
 D_{\rm raw}=C_{1,0}-C_{1,1}-C_{0,0}+C_{0,1}.
$$

Choose a polarity $\sigma\in\{1,-1\}$. In one polytope, encode $\lambda=s\mathbin{\mathrm{xor}}\eta$ with four linear inequalities,
$$
 \lambda\ge s-\eta,\quad \lambda\ge\eta-s,\quad
 \lambda\le s+\eta,\quad \lambda\le2-s-\eta.
$$
For $\sigma=1$ set the complement flag $\nu=1-\lambda$, and for $\sigma=-1$ set $\nu=\lambda$. Here $\nu$ is an affine expression. Add an integer residual $r$ with $\nu\le r\le(M-1)\nu$, and impose
$$
 xt-x^2-c_0^2-s(c_1^2-c_0^2)-r=Mq .
$$
For each choice of the other coordinates, there is exactly one possible residue $r$ and hence one quotient $q$ before applying the residual bounds. When $\nu=0$, the constraint counts the original congruence. When $\nu=1$, it counts its complement, including pairs with nonunit $x$. All those pairs belong to the known ambient population. The quotient bounds $|k|\le N+1$ and $|q|\le8M$ suffice with parameters chosen in their residue ranges.

For either polarity there is one complemented selector case for each sign class. Across those two cases the entire ambient population has size
$$
 B=wN^2(H+1).
$$
Thus the single count is $K_\sigma=B+\sigma D_{\rm raw}$. The ten integer coordinates are $s_0,y,t_0,z,k,s,\eta,\lambda,r,q$. The quadratic part $xt-x^2$ has rank two because $x,t$ are independent affine coordinates. Every coefficient and bound has $O(\log N)$ bits.

## Use a one-sided public threshold

Let $P_j-P_0$ be the public signed stationary correction, computed from the blocks of $x\equiv\pm c_0\pmod N$. Then $D(j)=D_{\rm raw}-(P_j-P_0)$. The count comparison can answer zero for
$$
 K_\sigma\le B+\sigma(P_j-P_0)+N/4
$$
and one for
$$
 K_\sigma\ge B+\sigma(P_j-P_0)+3N/4,
$$
with either answer in the gap. Rounding the lower threshold down and the upper threshold up gives integer thresholds.

Choose $\sigma$ uniformly once per isolation trial. In a singleton trial the signed stationary change, whenever nonzero, has one fixed sign because the reference membership bit and the affine sign of the isolated block are fixed. The polarity agrees with that sign with probability one half. The small-remainder analysis then makes the one-sided bit agree with the hidden membership change except on at most $1/32$ of the upper digits on a good aligned trial. Its heavy Fourier coefficient and gcd decoder are unchanged. All conclusions hold with repeated prime factors. [Fourier recovery](quarter-fourier.md)

## Checks

The diagnostic verified 30 packed counts and threshold comparisons against independently counted congruences for both polarities on three small moduli, including a repeated-prime input. It counts nonunit pairs in the complement. Run `python check_single_quadratic_count.py`.

## Comparison with the finite-subspace identity in grh-lean

At commit 72fead1 of lsernau/grh-lean, New.q_binomial_submodule
sums over subspaces U of an m-dimensional space over a field of
order q. Each dimension-i subspace has weight (-X)^i q^(i(i-1)/2).
The total equals the product of (1-q^j X) for 0<=j<m. Grouping by
dimension gives the Gaussian binomial coefficients. In the Hecke
transfer application an integrality constraint restricts U to the
subspaces of a single residue kernel. The same product then uses
the kernel dimension. At X=1 this sum is the indicator of zero
nullity, including the empty-product case.

The analogous step for the present comparison would be a public
encoding that collapses its weighted fibers while preserving the
threshold gap. Even over a known odd prime field and before applying
the interval filters, x(t-x)=c^2 for c nonzero is not an affine
subspace: it contains (c,2c) and (-c,-2c) but omits their midpoint.
Thus rank-nullity on the existing coordinates does not yield that
encoding. This leaves auxiliary-variable representations open.
The product identity computes a full weighted sum; specializing
the Gaussian coefficient polynomial at q=1 recovers ordinary choose
numbers and supplies no new fast individual-coefficient procedure.

[Subspace identity](https://github.com/lsernau/grh-lean/blob/72fead102f663d5c7a407964a13dc0cf9d977207/grh/src/Lib/General/Combinatorics.lean#L89) · [Kernel restriction](https://github.com/lsernau/grh-lean/blob/72fead102f663d5c7a407964a13dc0cf9d977207/grh/src/Lib/Weil/Unramified/HeckeTransfer.lean#L166)
