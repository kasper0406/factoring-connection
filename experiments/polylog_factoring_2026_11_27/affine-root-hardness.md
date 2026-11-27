# Exact affine root windows contain counting problems

The general affine-window oracle is substantially broader than the fixed unscaled window used by the deterministic bisection. Exact counts in the affine model are #P-hard, and their parities are ⊕P-hard, even when the complete factorization of the modulus is public. This redirects evaluator work toward the unscaled half-window, the prescribed doubling scales, or an approximate statistic. It does not rule out any of those restricted targets.

## Parsimonious encoding

Take a Subset Sum instance with positive weights $w_1,\ldots,w_k$, total $S$, and target $T$. The repository's deterministic quadratic-congruence construction sets $M=4S+4$, chooses pairwise coprime odd prime powers $q_i>8kM$, and puts $B=\prod_iq_i$. It constructs units $v_i\bmod q_i$ and

$$
\theta_i=(B/q_i)v_i,\qquad c=\sum_i\theta_i\pmod B.
$$

The element $c$ is a unit because modulo $q_i$ only $\theta_i$ survives. Every root of $X^2=c^2\pmod B$ is

$$
x_\epsilon=\sum_i\epsilon_i\theta_i\pmod B,\qquad \epsilon_i\in\{-1,1\}.
$$

Multiplication by $c^{-1}$ gives a bijection to the roots $u_\epsilon^2=1\pmod B$. Consequently a count of roots of unity whose public affine image $cu\bmod B$ lies in an interval is exactly the quadratic-root count already used in the parsimonious reduction. A polynomial-time exact evaluator for arbitrary public scales and intervals would give polynomial-time #SubsetSum counting. Taking the result modulo two gives the corresponding reduction from parity Subset Sum.

## Hardness already holds at a half-period window

The construction guarantees, for the signed sum $s_\epsilon=\sum_i\epsilon_iw_i$,

$$
\left\|\frac{x_\epsilon}{B}-\frac{s_\epsilon}{M}\right\|_{\mathbb R/\mathbb Z}<\frac{3}{16M}.
$$

All signed sums lie strictly between $-M/4$ and $M/4$, and adjacent possible sums differ by at least two. Put $s_T=2T-S$ and place one boundary at $(s_T+1)/M$. A cyclic interval of length $(B-1)/2$ ending at

$$
E_T=\left\lfloor B\frac{s_T+1}{M}\right\rfloor
$$

and take the residues $E_T-(B-1)/2+1,\ldots,E_T$ cyclically modulo $B$. This interval contains $x_\epsilon$ exactly when the subset represented by $\epsilon$ has sum at most $T$. A qualifying signed sum is at least $13/(16M)$ below the upper boundary, and the next possible signed sum is at least that far above it. The opposite boundary lies outside the occupied quarter-circle arc. Since $B>8M$, rounding the endpoints cannot change membership. Thus one affine half-window count computes the cumulative Subset Sum count. Two adjacent targets recover the exact-target count. The modulus and every query parameter have polynomial bit length, and its full factorization is part of the construction.

Translate the cyclic interval to $[0,(B-1)/2)$ to put the query exactly in the form $z=au+b\bmod N$. Thus the result applies directly to a proposed evaluator for arbitrary affine maps, such as the strongest isolation oracle. It also applies to arbitrary public-square root counts through $z=cu$. It does **not** apply to the deterministic fixed-window query with scale $a=1$, or to the doubling-only family $a=2^i$. Those restrictions prevent the construction from choosing the local CRT weights independently.

## Check

The check constructs the prime-power moduli without a large-prime search, enumerates every CRT sign vector, verifies $u^2=1\pmod B$, and compares the half-window count with the cumulative subset count. It covers 10,546 target instances, including ten weights and moduli up to 176 bits.

```sh
python experiments/polylog_factoring_2026_11_27/check_affine_root_hardness.py
```

[Underlying quadratic-count proof](../../notes/floorhard_results.md#21-arbitrary-intervals) · [Affine isolation target](affine-isolation.md) · [Fixed unscaled window](fixed-window-split.md)
