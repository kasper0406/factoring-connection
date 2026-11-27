# Several public congruences for shifted count prediction

The short-summary bound leaves room for a long CRT transcript of $N$. On the same balanced distribution, independent random auxiliary moduli provide little additional information even when their combined modulus is close to $N$ on an exponent scale. The proof applies to the distribution of the *products*, so it does not require either individual prime to be equidistributed over the combined CRT modulus. [Dense bit and recovery](dense-count-prediction.md)

## Adding one random modulus

Take $X$ dyadic and draw two distinct primes uniformly from $[X,5X/4]$. Their unordered product is uniform on a set $\mathcal S$ of size $T\asymp X^2/(\log X)^2$. Fix any already observed modulus $B$. Let $\lambda(N)$ be a label taking at most $V$ values. For example, retain a vector of $t$ dense count signs at shifts chosen independently of the next modulus, giving $V=2^t$. Keep the actual joint law of the cells $(a,\lambda)$ with $a=N\bmod B$.

Choose an independent prime $D$ uniformly from a pool $\mathcal D$ in $[R,2R]$, coprime to $B$. In each cell write $N=a+Bj$. Let $v_j$ be its indicator in the integer $j$ interval of length $H=O(1+X^2/B)$. The histogram of $N\bmod D$ in this cell is a permutation of the histogram of $j\bmod D$. Let $T_c$ be the number of elements in cell $c$. Orthogonality gives
$$D\sum_{r\bmod D}\left(\sum_{j\equiv r\ (D)}v_j-T_c/D\right)^2=\sum_{h=1}^{D-1}\left|\sum_j v_j e_D(hj)\right|^2.$$
Write $E_D$ for the sum on the right over all cells. By Cauchy–Schwarz, the total variation distance between the actual law and the law with an independent uniform residue modulo $D$ obeys $4T^2\mathrm{TV}_D^2\le BV E_D$. All nonzero fractions $h/D$ in the pool are reduced and spaced by at least $1/(4R^2)$. The additive large sieve gives $\sum_D E_D\le (H+4R^2)T$. Thus, when $|\mathcal D|\asymp R/\log R$,
$$\mathbb E_D\,\mathrm{TV}_D\ll\sqrt{\frac{BV(H+R^2)}{T|\mathcal D|}}\ll\sqrt V\,\log X\sqrt{\log R\left(\frac1R+\frac{BR}{X^2}\right)}.$$
The last inequality is used with $B\le O(X^2)$ and $R\ge2$. It bounds the increase in success of any prediction rule, even if computation on these residues is unrestricted. A uniform extra residue adds no information in the comparison law. [Additive large sieve](https://personal.science.psu.edu/rcv4/personal/Publications/large_sieve.pdf)

## Consequences at the lattice scale

Start with $B=\ell$ for an auxiliary prime $\ell\asymp\sqrt X$, large enough for native lattice recovery. For one fixed shift, $V=2$. With $R=X^{3/2-\eta}$, $0<\eta<3/2$, the combined residue modulus is of order $X^{2-\eta}$. The increase in mean success is at most
$$O\!\left((\log X)^{3/2}\left(X^{-3/4+\eta/2}+X^{-\eta/2}\right)\right).$$
The bound is uniform over fixed shifts, so it also holds for uniform shifts sampled independently of $D$ and $N$. Together with the first-modulus bias bound, advantage remains exponentially small in the input bit length. For $t$ such shifts and an arbitrary decision on their joint label, the increase in success over using $N\bmod\ell$ alone is bounded by the same expression times $2^{t/2}$. The comparison preserves correlation between the different signs.

The effect also composes for several independently chosen primes, sampled without replacement. Condition on the previous prime choices and use $B$ equal to their product times $\ell$. Group the *actual* products by the previously visible residues. A hybrid replacing the newest residue with an independent uniform one incurs the bound above, and replacing all the residues incurs the sum of the bounds. The next prime is independent of the product and of the fixed label conditional on the previous prime choices. For example, choose the extra primes at a fixed scale $R=(\log X)^C$ and stop after $\lfloor(3/2-\eta)\log X/\log(2R)\rfloor$ draws. The combined modulus has size $X^{2-\eta-o(1)}$ and stays below $O(X^{2-\eta})$. For $t=O(\log\log X)$ and sufficiently large constant $C$, the extra residues give less than any prescribed inverse-polynomial advantage in $\log X$. The prime pool has enough unused primes at each step. Each CRT residue takes polynomial bit time to compute, and the number of draws is polynomial as well.

If the shifts may be chosen *after* seeing the extra residue, keep the full slope $(p+q)\bmod\ell$ as the label. This uses $V=\ell$ and preserves all possible signs, including choices that depend on the public transcript. With $R=X^{1-\eta}$ and $0<\eta<1/2$, the comparison still has power-saving error and the combined modulus is $X^{3/2-\eta}$. The full-slope comparison has a separate limit. In the independent field model from the preceding note each pair $(N\bmod\ell,(p+q)\bmod\ell)$ has probability at most $2/(\ell-1)^2$. After adjoining an independent uniform residue modulo $D$, the probability assigned to the actual support of at most $T$ triples is at most $2T/((\ell-1)^2D)$. Including the preceding note's mean error for the field model shows why this stronger comparison cannot extend to $\ell^2D\gg T$.

These are distributional results for randomly chosen exogenous moduli. The bounds leave open deterministic arithmetic chosen using the full input, such as a modulus tailored to a relation involving the two factors.

## Finite checks

The diagnostic builds exact residue histograms on two balanced prime ranges, using one and three shifts. For five choices of the previously retained modulus, it checks the total variation bound for optimal joint prediction, the squared-energy formula, and the additive large-sieve budget over every prime in the auxiliary pool. These finite tables check the inequalities. Their prediction gains are measured at the small fixture sizes. [Residue table checks](congruence_summary_check.py)

```sh
uv run --no-project --with sympy --with python-flint==0.9.0 --with gmpy2 python experiments/polylog_factoring_2026_11_27/congruence_summary_check.py
```
