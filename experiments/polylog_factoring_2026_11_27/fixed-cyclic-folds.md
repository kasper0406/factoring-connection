# Fixed cyclic folds miss the central binomial split infinitely often

For every finite collection of positive widths, there are infinitely many balanced semiprimes for which all the unweighted cyclic folds of the middle binomial coefficient are units. The exact middle coefficient still gives a proper divisor on every one of these inputs. This extends the roadmap's two finite counterexamples to arbitrary fixed width collections. [Roadmap arithmetic control](https://github.com/kasper0406/factoring-connection/blob/9514cbfa013cd9ae3f030044db1dc3fb6bb4a9ed/notes/direct_split_poly.md)

## Prime pairs

Let $L$ be the least common multiple of the widths. Bounded gaps in the fixed progression $1\bmod 2L$ give infinitely many pairs $p<q$ from that progression with one fixed positive difference $2h$. In particular $L\mid h$. For this standard consequence of the Maynard sieve, take a sufficiently large fixed admissible tuple of forms $2L(t+i k!)+1$ for $1\le i\le k$. Primes dividing $2Lk!$ never obstruct the tuple, and larger primes cannot exclude every residue. The fixed-form case of Maynard's theorem gives infinitely many $t$ for which at least two forms are prime. There are finitely many differences in the tuple, so one difference recurs infinitely often. [Maynard, Dense clusters of primes in subsets, Theorem 3.1](https://arxiv.org/abs/1405.2593)

For large $p$, put $m=p+h$ and $N=pq=m^2-h^2$. We have $B=\lfloor\sqrt N\rfloor=m-1$, $q<2p$ and $B\equiv0\pmod d$ for every tested width $d$. Let $F_d$ be the constant coefficient of $(1+x)^{2B}\bmod(x^d-1)$. This is exactly the unweighted middle-position fold.

## Residue at the smaller prime

Modulo $p$ in the quotient ring, $(1+x)^p=1+x^p=1+x$ since $p\equiv1\pmod d$. Thus
$$(1+x)^{2B}=(1+x)^{2p+2h-2}=(1+x)^{2h}.$$
The constant coefficient of the last polynomial in the cyclic quotient is a fixed positive integer, depending only on $d,h$. Its reduction is nonzero for sufficiently large $p$.

## Residue at the larger prime

The same Frobenius identity at $q$ gives $(1+x)^{2B}=(1+x)^{q-1-2h}$ because $2B=q+(q-2h-2)$. For large $q$ the exponent on the right is positive. Since $q\equiv1\pmod d$, the $d$-th roots of unity are distinct in $\mathbb F_q$ and the constant coefficient is the average of the evaluations there. The evaluation at $-1$, when present, is zero. Every other evaluation is $(1+\zeta)^{-2h}$ by Fermat's little theorem.

Consider the fixed characteristic-zero constant
$$C_{d,h}=\frac1d\sum_{\substack{\zeta^d=1\\\zeta\ne-1}}(1+\zeta)^{-2h}.$$
The Galois group permutes the summands, so this is rational. For $\zeta=e^{2\pi i j/d}$ the summand is
$$(2\cos(\pi j/d))^{-2h}e^{-2\pi i j h/d}=(2\cos(\pi j/d))^{-2h}>0$$
because $d\mid h$. Hence $C_{d,h}>0$. For all sufficiently large $q$, reduction of the cyclotomic integers with $d$ and all $1+\zeta$ for $\zeta\ne-1$ inverted gives the same average in $\mathbb F_q$. The nonzero rational then reduces to a nonzero field element. Therefore $F_d$ is a unit for every width in the collection once the primes are large enough.

Meanwhile $p\le B<q$ and $2p\le2B<3p$, so $p$ does not divide the exact central coefficient. We also have $q\le2B<2q$, so $q$ does divide it. Its gcd with $N$ therefore equals $q$.

## Scope and checks

The theorem rules out a universal split by taking individual gcds from a fixed finite collection of unweighted middle-position folds. It leaves input-dependent widths, combinations of outputs and other readouts open. The counterexample family has fixed prime gap and admits immediate Fermat factoring for sufficiently large primes.

The check script computes the rational average with exact fractions in the cyclic quotient, and compares its two predicted local residues against modular cyclic powering. It tested 284 pairs across the width collections $(1,2,3)$, $(2,4,8)$ and $(3,5,6)$. All tests agreed and every folded output was a unit. Run `python check_cyclic_folds.py`.
