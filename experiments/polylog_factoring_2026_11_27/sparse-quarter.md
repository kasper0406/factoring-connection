# Sparse polynomial quarter predictors reduce to direct gcd tests

Fix a number of terms $s$. On balanced semiprimes, a predictor $g_0(P(c)\bmod N)$ with $P(c)=\sum_{i=1}^s a_i c^{e_i}$ and inverse-polynomial correlation with the hidden quarter bit yields a direct gcd procedure with inverse-polynomial success. Exponents and coefficients may depend on $N$. Each sampled scale costs $O(2^s)$ subset sums times polynomially many harmonic tests, with modular powering for the terms. The proof does not require bounded polynomial degree. [Single-power case](quarter-all-powers.md)

## Separate exponent groups locally

For a prime $r\mid N$, start with $d=r-1$, and include the public exponents 0 and 1 alongside the $s$ exponents. Whenever two exponents differ modulo $d$ but $\gcd(e-f,d)>d/r^\varepsilon$, replace $d$ by that gcd. Here $\varepsilon=1/(4(s+1))$. Each replacement merges at least two residue classes, keeps every previously merged pair merged, and loses a factor less than $r^\varepsilon$. Thus there are at most $s+1$ replacements, and the terminal $d=d_r$ is at least $(r-1)/r^{1/4}$. In particular 0 and 1 stay distinct for large enough $r$.

Condition the local scale on a coset of the subgroup of order $d_r$. Its uniform parameter can be written as $y^{(r-1)/d_r}$ for uniform $y\in\mathbb F_r^\times$. Group polynomial terms by exponent modulo $d_r$. For the resulting exponents in $y$, the gcd of any difference between distinct groups with $r-1$ is at most $(r-1)r^{-\varepsilon}$. Including the group of exponent 0 gives the same bound for each nonconstant exponent individually. Bourgain's estimate therefore gives a normalized additive sum of magnitude $O(r^{-\gamma_s})$ whenever any nonconstant coefficient survives. In that case the value distribution also puts probability at most $1/r+O(r^{-\gamma_s})$ on any one field element, by additive Fourier inversion. The same estimates apply to subsets of the exponents translated by minus 1, since the original grouping included 1. [Bourgain's sparse exponential-sum theorem](https://www.numdam.org/articles/10.1016/j.crma.2004.06.013/) · [Gcd hypotheses](https://onlinelibrary.wiley.com/doi/full/10.1002/mana.202300426)

## Extract a small resonance

Truncate the hidden quarter sign at harmonics $0<|k|\le K$ and expand the predictor in its full additive Fourier series. The truncation costs $O(K^{-1/2}+1/N)$. The conditioned phase $hP(c)-kuc$ has cancellation unless it is constant on both local cosets. With $S_r=\{i:e_i\equiv1\pmod{d_r}\}$, this constancy requires
$$hV_{S_r}(c)\equiv ku\pmod r,\qquad
V_S(c)=\sum_{i\in S}a_i c^{e_i-1}.$$
The values $V_{S_r}$ are constant on their respective cosets. Check the bounded nonzero $k$ for proper gcds first, so they can be assumed units. For each $k$ and coset pair there is at most one $h$ modulo $N$ satisfying both conditions, because then both local $V_{S_r}$ must be nonzero. Resonances with $|h|_N>J$ cost $O(\log(2K)/J)$ by the quarter coefficient bound. The contribution of nonconstant phases is $O(p^{-\gamma_s}\log N\log(2K))$. For $J,K$ polynomial in the inverse promised correlation and the bit length, the argument forces inverse-polynomial probability of a resonance with $0<|h|\le J$, $0<|k|\le K$.

## Convert resonances even when the local groups differ

Orient the hidden root as $u\equiv1\pmod p$ and $-1\pmod q$. The public test enumerates all subsets $S$ and bounded signed $h,k$, checking $\gcd(hV_S(c)-k,N)$. If neither $S_p$ nor the test with reversed $k$ at $S_q$ splits at a resonance, the cross values must satisfy
$$hV_{S_p}\equiv k\pmod q,\qquad hV_{S_q}\equiv-k\pmod p.$$
These cross equalities are rare if the cross polynomial is nonconstant on its corresponding coset. The point bound above and a union bound over bounded harmonics charge at most $O(JK(p^{-\gamma_s}+q^{-\gamma_s}))$ probability. This remains exponentially small in the bit length on balanced inputs for fixed $s$.

If both cross polynomials are constant on the cosets, their constant terms come exactly from the intersection $S_p\cap S_q$, because the exponents $e_i-1$ are constant at a prime precisely for its own subset. Hence $hV_{S_p\cap S_q}$ equals $-k$ at $p$ and $k$ at $q$. The intersection test gives a proper gcd. The enumerator needs neither the prime factors nor the subgroup orders. Negative exponents in $V_S$ are evaluated on units.

This leaves growing term counts, exponent choices that depend on the scale, and more general arithmetic circuits outside the result. The diagnostic verifies direct gcd recovery on three planted two-term predictors and a separate case where the local subsets differ and the intersection is needed for the stated harmonic pair. Run `python check_sparse_quarter.py`.
