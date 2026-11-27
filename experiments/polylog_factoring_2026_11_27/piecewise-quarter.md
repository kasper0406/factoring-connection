# Piecewise sparse predictors and direct gcd tests

On balanced semiprimes, the sparse-polynomial predictor result extends to polynomially many public arithmetic progression cells. A predictor may choose a different fixed-size sum of modular powers on each cell. Inverse-polynomial correlation still implies inverse-polynomial success for direct gcd tests. In particular this covers choosing an exponent from $\lfloor N/c\rfloor$. [Sparse predictor proof](sparse-quarter.md) · [Scale-dependent exponent screen](adaptive-exponent-screen.md)

## Add the cell frequencies

Suppose the unit scales are partitioned by at most $C$ public cells of the form $[a,b]\cap(r\bmod d)$, with polynomially bounded $C,d$ and a fixed bound on terms in each polynomial. Check any nonunit step by gcd. On each cell the prediction is the quarter rule applied to the cell's sparse polynomial $P(c)$.

Use bounded Fejér approximations of the target and of each cell indicator. Keep the nonzero target frequencies $ku$, $0<|k|\le K$, omitting its $O(1/N)$ mean. This costs $O(K^{-1/2}+1/N)$ in correlation. In the progression coordinate a cell is an interval modulo $N$. With harmonic bound $L$, its indicator approximation has $L^2$ error $O(L^{-1/2})$ and Fourier absolute sum $O(\log(2L))$. Its scalar frequencies are $\ell=j d^{-1}\bmod N$ for $|j|\le L$. Replacing all indicators costs $O(C/\sqrt L)$ in correlation. These estimates apply over uniform units by $N/\varphi(N)=O(1)$.

On each cell use the local subgroup grouping of the sparse proof. After expanding the predictor's quarter rule, the phase is $hP(c)+\ell c-kuc$. The subgroup lemma gives a power saving unless the phase is constant on both local cosets. For the local exponent-1 subset $S_r$ this requires
$$hV_{S_r}(c)=ku-\ell\pmod r,$$
where $V_S(c)=\sum_{i\in S}a_i c^{e_i-1}$.

## Bound the tail and test for factors

The two targets cannot both vanish. Otherwise $j\equiv dku\pmod N$. For sufficiently large balanced inputs with $dK+L<p$, the small $j$ would equal both $dk$ and $-dk$ as integers. Therefore any resonant $h$ belongs to one prescribed residue class modulo at least one prime. Its quarter coefficients at $|h|>J$ sum in absolute value to $O(1/J+\log N/p)$. Bound the first term on each side by $1/J$, then sum reciprocals on a progression of step at least $p$.

The total contribution of nonconstant phases and large resonances is at most
$$O\left(C\log(2L)\log(2K)\left(p^{-\gamma}\log N+1/J+\log N/p\right)\right)$$
for fixed $\gamma>0$ from the sparse bound. Take $K,L,J$ polynomial in $C$, inverse promised correlation and the bit length. The discarded terms can be made smaller than half an inverse-polynomial signal on sufficiently large inputs. As the absolute coefficient sum for each bounded harmonic triple is polynomial, a small resonance among the cells and harmonics must occur with inverse-polynomial probability on a uniform scale.

At a small resonance, orient $u$ so the local targets are $k-\ell$ and $-k-\ell$. They differ by the unit $2k$, after testing small integers by gcd. Try the subset gcds with these shifted targets. If both local exponent-1 subsets fail by matching their target at both primes, the cross-polynomial point bound makes nonconstant coincidences exponentially rare, even over the polynomially many choices. If both cross polynomials are constant, the intersection subset realizes the targets in reverse order and gives a gcd. The public procedure tests every cell's polynomial on each sampled scale. It does not need to identify the hidden cosets.

## Quotient corollary

For promised correlation $\delta$, discard $c<\delta N/16$. Their unit probability is at most $\delta/8$ up to an endpoint correction. On the remaining scales, $\lfloor N/c\rfloor\le16/\delta$, so its values partition the domain into polynomially many intervals. Any polynomial-time rule choosing a fixed-size sparse polynomial from that quotient falls within the result. The full representative, modular inverse, and highly fragmented bit-count rules remain outside it.

The diagnostic exercises subset-intersection recovery at a nonzero progression harmonic $\ell=2/3\bmod N$ as well as at zero. Run `python check_sparse_quarter.py`.
