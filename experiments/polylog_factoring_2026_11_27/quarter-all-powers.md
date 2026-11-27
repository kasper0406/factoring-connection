# Fixed power quarter predictors reduce to direct gcd tests

On balanced semiprimes, this rules out a distinct advantage from any fixed exponent in the predictor $g_0(c^e\bmod N)$. A fixed exponent here may depend on the public modulus, but is held fixed across the unit scale $c$. If this predictor has inverse-polynomial correlation with the hidden quarter bit, direct gcds involving $c^{e-1}$ have inverse-polynomial success. The conclusion uses a known sparse exponential-sum theorem and holds asymptotically. [Previous large-kernel argument](quarter-power-predictor.md)

Let $N=pq$ with $p<q\le2p$, $u$ a mixed root, and $g_0$ the public quarter sign. Truncate the Fourier expansion of $g_0(cu)$ to $0<|k|\le K$. The discarded part costs $O(K^{-1/2}+1/N)$ on uniform units. Expand $g_0(c^e)$ in its entire Fourier series, of absolute coefficient sum $O(\log N)$. The average of each resulting pair factors over the two primes, with local terms of the shape
$$\frac1{r-1}\sum_{x\in\mathbb F_r^\times} e_r(A x^e+B x),\qquad r=p,q.$$
For the small nonzero target harmonics, $B$ is always nonzero unless a gcd of a tested $k$ already splits $N$. The coefficient $A$ may be zero.

## Two cancellation cases

1. If at either prime $r$ we have $\gcd(e,r-1)\ge r^{3/4}$, fix the value of $x^e$. Its fiber is a multiplicative coset of that size. The linear phase $e_r(Bx)$ has average over the coset of magnitude at most $\sqrt r/\gcd(e,r-1)\le r^{-1/4}$ by the Gauss-sum bound. The $x^e$ phase is constant on this fiber, so this bounds the whole local sum.
2. If at either prime $r$ both $\gcd(e,r-1)$ and $\gcd(e-1,r-1)$ are below $r^{3/4}$, Bourgain's sparse exponential-sum estimate gives a local bound $O(r^{-\gamma})$ for an absolute $\gamma>0$. The two exponents are $1$ and the reduction of $e$ modulo $r-1$. The relevant gcds are those of the two exponents and their difference with $r-1$. If $A=0$ the linear sum has magnitude $1/(r-1)$ directly. [Bourgain, Mordell type exponential sum estimates](https://www.numdam.org/articles/10.1016/j.crma.2004.06.013/) · [Statement of the gcd hypotheses](https://onlinelibrary.wiley.com/doi/full/10.1002/mana.202300426)

In either case the other local factor is bounded by one. Summation over the harmonics costs only $O(\log(2K)\log N)$. For $K$ polynomial in the inverse correlation and $\log N$, a power saving at a balanced prime is exponentially small in the bit length. Thus these cases cannot sustain inverse-polynomial correlation.

## Remaining exponents

If neither cancellation case holds, both local kernels of $x\mapsto x^{e-1}$ have size at least $r^{3/4}$. The previous resonance bound then has $\eta\le p^{-1/4}$. Choose the truncation bounds $J,K$ polynomial in the inverse promised correlation and the bit length. The small resonance event
$$h c^{e-1}\equiv ku\pmod N,\qquad 0<|h|\le J,\quad 0<|k|\le K,$$
must occur with inverse-polynomial probability. Test each bounded pair by $\gcd(hc^{e-1}-k,N)$, after checking the small integers themselves by gcd. A resonance with unit $h,k$ splits because $u$ has opposite signs at the two primes. The public procedure need not identify the local kernel case and does not call the quarter predictor.

This covers the single fixed-power quarter rule. It leaves exponent choices depending on the scale, compositions of several powers, and more general modular maps open. It establishes no factoring method that selects an exponent giving a useful correlation.

## Checks

The existing exhaustive planted controls exercise the direct gcd outcome. Identity and arbitrary exponents serve as finite controls only, since the sparse-sum theorem supplies the asymptotic cancellation. Run `python check_quarter_power_predictor.py`.
