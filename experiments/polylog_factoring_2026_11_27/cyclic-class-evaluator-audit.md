# What the cyclic class encoding does and does not compress

The common-discriminant queries are powers of one public class of known
power-of-two order. This makes every queried **form** computable in polynomial
bit time. It does not give a composition rule for its clipped disk count. Two
natural uses of the group structure can be excluded precisely. A calculation
with class characters cannot distinguish the queries at all, and retaining
only a low-order quotient of the cyclic group loses the required bit on the
published controls. A nonstationary fixed-exponent evaluator remains open.
[Class encoding](cyclic-class-queries.md) · [Full-orbit recurrence bound](class-exponent-period.md)

## Genus characters are constant on every query

Write the cyclic subgroup as \(G=\langle g\rangle\), of order \(L=2^k\).
Every query in the fixed-window reduction is \(g^e\) for an odd exponent
\(e\). A genus character is a homomorphism from the proper class group to
\(\{\pm1\}\). On this subgroup,

$$
\chi(g^e)=\chi(g)^e=\chi(g).
$$

Thus all genus characters, including those coming from the many prime factors
of the nonfundamental discriminant, have exactly the same value on every query.
Computing a genus vector or any other exponent-two quotient cannot supply the
window parity. This conclusion uses only that the exponents are odd. It does
not assume that the full class group is cyclic.

Higher-order characters retain more exponent bits. They still need enough
resolution to separate queries with different counts. On the three controls
below, two actual queried powers agree in the indicated quotient and have
opposite corrected disk parities.

| \(N\) | input bits | common quotient | query starts | mixed parities |
|---:|---:|---:|---:|---:|
| 11,009 | 14 | exponent modulo \(2^{13}\) | 436, 2,726 | 1, 0 |
| 28,829 | 15 | exponent modulo \(2^{15}\) | 1,213, 6,924 | 1, 0 |
| 471,311 | 19 | exponent modulo \(2^{19}\) | 2, 109,055 | 1, 0 |

These are finite counterexamples to the proposed small-quotient evaluator, not
an asymptotic lower bound. The diagnostic uses the supplied factors only to
label the true mixed-root parity. Form construction and exponent recovery use
the public algorithms.

## Why class powering does not power the count

Ideal multiplication gives a short basis for \(g^{e+f}\) from bases for
\(g^e\) and \(g^f\). Multiplying lattice elements multiplies their norms. It
therefore maps two inequalities of the form \(Q\le NB_0\) to a product bound
of order \((NB_0)^2\), rather than to the original two-dimensional disk.
Multiplying the corresponding theta series likewise produces a
four-dimensional convolution. Neither operation gives the clipped
two-dimensional count at \(e+f\). Repeated squaring is consequently a
polynomial-time form constructor, as already implemented, without being a
count evaluator.

A baby-step/giant-step table of **counts** has no verified combine operation.
If it instead realizes a fixed stationary linear state over the entire odd
orbit, the existing characteristic-two recurrence theorem requires more than
\(L/4\) states. Since \(L=\Theta(N^2)\), that route is exponential in the bit
length. This theorem does not cover an algorithm specialized to the
\(O(\log N)\) adaptively selected exponents.

The Hecke reformulation has the same missing input. Acting at index \(N\) on a
binary theta series introduces the ideals of norm \(N\). The public roots give
two of their local pairings, while the other pairings are the hidden mixed
classes. Standard factor-by-factor Hecke evaluation would therefore assume the
decomposition being sought. This is a diagnostic about that implementation,
not a proof that every Hecke algorithm must factor first.

## Decision

Do not prioritize genus vectors, fixed low-order character quotients, or a
baby-step/giant-step table whose entries are only class powers. Reopen the
cyclic route for a formula that computes the clipped count at one supplied odd
power while explicitly charging its boundary work. Such a formula would also
be a direct evaluator for the equivalent root-window problem.

Reproduce the finite quotient collisions from this directory with:

```sh
python check_cyclic_character_quotients.py
```
