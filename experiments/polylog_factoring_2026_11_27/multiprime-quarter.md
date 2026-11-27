# Weak parity prediction for a fixed number of prime divisors

For an odd input with a fixed number $k\ge2$ of distinct prime divisors, weak prediction of scaled quarter parity also gives factoring, provided the correlation exceeds an explicit public Fourier contribution. The procedure works with repeated prime factors. Its bounds are polynomial in the bit length and inverse promised correlation for each fixed $k$. [Two-prime case](quarter-fourier.md)

## Approximate the product of quarter signs

In the proof write $m_1,\ldots,m_k$ for the complete prime-power components of $N$. Choose one root from each opposite-sign pair by requiring its first CRT sign to be positive. The public pair corresponds to all positive signs. The remaining $R=2^{k-1}-1$ roots have sign vectors $\sigma\in\{\pm1\}^k$ whose first coordinate is positive and whose coordinates are not all positive. Write $u_\sigma$ for the corresponding root. After removing the public pair, the hidden parity sign is
$$G(c)=\prod_\sigma g_0(cu_\sigma),$$
where $g_0$ is the public quarter-window sign.

Use the degree-$K$ Fejér approximation to $g_0$, with coefficients
$$A_h=(1-|h|/(K+1))\widehat g_0(h),\qquad |h|\le K<N.$$
The approximation has absolute value at most one because its kernel is nonnegative and normalized over the grid. Its $L^2$ error is $O(K^{-1/2})$ by the $1/|h|_N$ coefficient bound, including the Fejér attenuation of retained coefficients. Substituting it in all $R$ factors gives a polynomial $T(c)$ with $\|T\|_\infty\le1$ and $\|G-T\|_2=O(R/\sqrt K)$, by telescoping the product and using that multiplication by each $u_\sigma$ permutes the residues.

## Public and factor-revealing frequencies

Expanding the product, every Fourier term has CRT frequency vector
$$v_i=\sum_\sigma h_\sigma\sigma_i,\qquad |v_i|\le U=RK.$$
The coefficients at these integer vectors can be computed without factoring by convolution of the abstract sign vectors. There are at most $(2U+1)^k$ states. The actual scalar frequency is the CRT lift of $v$. If any prime factor is at most $2U$, bounded gcd tests already split the input. Otherwise the vector-to-frequency map is injective on this box.

Vectors of the form $(t,\ldots,t)$ have the public scalar frequency $t$. Let $B(c)$ be their contribution to $T(c)$. It is computable to polynomial precision from $N$, $k$ and $K$ without the primes. Every other vector has unequal coordinates. For its scalar frequency $f$, a test of $\gcd(f-t,N)$ splits whenever $t$ is one of its coordinates. The search only needs to try the integers in $[-U,U]$. Public low frequencies can be discarded before these tests.

## Decode centered correlation

For a predictor with expected sign $F(c)$, extend $F$ by zero on nonunits, using proper gcds opportunistically. Suppose
$$\left|\mathbb E_{c\bmod N}F(c)(G(c)-B(c))\right|\ge\delta.$$
Choose $K=O(R^2/\delta^2)$ with a sufficient constant so the product approximation error is at most $\delta/2$. Then $F$ has correlation at least $\delta/2$ in absolute value with $T-B$. Its Fourier support has size at most $D=(2U+1)^k$ and its $L^2$ norm is at most one by orthogonal projection of $T$. Thus some factor-revealing frequency in its support has $|\widehat F(f)|\ge\delta/(2\sqrt D)$. Sampled Fourier localization and the bounded gcd tests recover a divisor. Computing the public contribution and the search both take polynomial time for fixed $k$. Trying all $k$ up to a fixed promised upper bound lets the decoder omit advance knowledge of $k$. False candidates are checked by gcd.

The public subtraction matters because products of hidden harmonics can have all their local signed sums equal. A raw prediction accuracy by itself need not certify correlation outside these public frequencies. This result gives a testable centered target and leaves growing $k$ and native predictors open.

The centered target is nonvacuous for fixed $k$ with large prime components. In the continuous torus limit, the product of quarter signs depends on $k$ independent coordinates, and is not a function of their sum alone. Choose a generic point on a sign-change hyperplane of one mixed sign vector, away from the other sign-change hyperplanes. Moving across it in a direction preserving the coordinate sum flips exactly one factor. This happens on a positive-measure family of nearby slices. Hence the variance after projection onto diagonal frequencies is positive. Riemann sums and the Fejér approximation transfer this to sufficiently large grids. This observation concerns exact or planted predictors, not a way to compute a useful predictor from $N$.

## Checks

The diagnostic built the abstract vector convolutions for three- and four-prime inputs, including one repeated-prime input. It checked that every nondiagonal vector gave a proper gcd and that the expanded spectra reproduced the direct products at sampled scales. Run `python check_multiprime_quarter.py`.
