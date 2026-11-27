# Sparse predictors with Fourier approximations to their branch masks

The direct-gcd result for piecewise sparse quarter predictors extends to any polynomial list of public branch masks with polynomially bounded Fourier absolute sums. Frequencies of those masks may be large and need not be supplied. This includes fixed-depth trees comparing modular affine values with interval endpoints. [Arithmetic progression case](piecewise-quarter.md) · [Sampled frequency search](quarter-fourier.md)

## Construct the public approximations

A branch mask $I(c)$ indicates which formula is selected. Assume it is evaluable in polynomial bit time and that $\sum_f|\widehat I(f)|\le B$ for a public polynomial bound $B$. If all coefficients of magnitude below $\tau$ are discarded, Parseval gives squared error at most $\tau B$. The shared-sample Fourier search returns every coefficient of magnitude at least $\tau$, with a polynomial-size list and bounded failure probability. Estimate the returned coefficients by their bounded additive averages. Estimating each to precision proportional to the desired $L^2$ error divided by the square root of the list size suffices. This constructs a Fourier approximation with polynomially many known frequencies and polynomial absolute coefficient sum. The same argument applies when such an approximation is available by another method.

For a polynomial list of masks, make the sum of their approximation errors smaller than a fixed fraction of the promised correlation. The error over unit scales differs by at most a constant factor on balanced semiprimes. In the correlation calculation it is multiplied by two bounded signs. The approximations themselves need not take values in $[0,1]$.

## Handle any mask frequency

As in the progression proof, a public mask frequency $\ell$ contributes the term $\ell c$ to the local phase. The target uses a small nonzero frequency $ku$, where $u$ is the mixed root. If $\ell=ku\bmod N$, then $\gcd(\ell-k,N)$ gives a factor. Test all public mask frequencies against the bounded target harmonics first. A nonunit harmonic gives a direct gcd as well.

In the remaining cases, at least one prime has a nonzero local target $ku-\ell$. Thus the possible resonant predictor harmonics form one residue class modulo at least one prime. Their large-frequency weight is bounded as before by $O(1/J+\log N/p)$. The sparse subgroup cancellation and the intersection test for the two local exponent subsets are unchanged. All harmonic absolute sums of the branch approximation are polynomial. Consequently inverse-polynomial correlation implies enough small resonances for the direct gcd procedure to split with inverse-polynomial probability, unless the public frequency precheck already splits.

## Application and checks

An interval test on $a c+b\bmod N$ has Fourier absolute sum $O(\log N)$ when $a$ is a unit, by the geometric-sum bound and permutation of frequencies. A nonunit $a$ is either zero modulo $N$ or gives a factor by gcd. A fixed-depth tree of such tests has masks which are products of a fixed number of interval indicators or their complements. Fourier convolution bounds their absolute sums by the products of the individual bounds. The resulting theorem covers scale-dependent formula choices using comparisons of arbitrary public modular affine representatives.

The diagnostic built products of two affine interval masks. On three semiprimes the sampled search found a planted hidden frequency and the public frequency gcd recovered a divisor. It also verified the tail-energy and product bounds with exact spectra used as controls. The planted multiplier uses fixture factors. These checks exercise the reduction, not a factor-blind source of a hidden frequency. Run `python check_spectral_branch_masks.py`.
