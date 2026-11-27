# Pointwise Salié coefficients and quadratic residuosity

Grouping the public Salié summands by their additive phase gives an exact coefficient formula. On odd squarefree semiprimes, pointwise evaluation of these coefficients is equivalent to the ordinary quadratic residuosity decision problem, allowing gcd shortcuts. This identifies the arithmetic needed for pointwise evaluation. It leaves open how to compute the oscillatory aggregate from that information. [Normalized Salié target](salie-quarter.md)

For a unit $c\bmod N$, put $b=c^2/4$ and write
$$W_c(t)=\sum_{m\in U_N\,:,m+b/m=t}(m/N).$$
Then $S_N(1,b)=\sum_t W_c(t)e_N(t)$. For $N=pq$ squarefree and a mixed root $u$, the exact coefficient is
$$
W_c(t)=\chi_N(t+c)+\chi_N(t-c)+\chi_N(t+cu)+\chi_N(t-cu).
$$

To prove this, at one prime $r$ write $a=c/2$. The equation is $m^2-tm+a^2=0$, so its roots are units. At a root, $(m+a)^2/m=t+2a$. If $t\ne\pm2a$, the two signs $\chi_r(t+2a)$ and $\chi_r(t-2a)$ are equal precisely when the discriminant is a square. When equal, there are two roots and each contributes the common sign. When opposite, there is no root and the two signs sum to zero. At the endpoints, the double root contributes the one nonzero sign among $\chi_r(t\pm2a)$ and is counted once. Thus the local coefficient is their sum. Multiply the two local coefficients by CRT and expand. The two mixed translations use $u$ because their choices of sign disagree at the primes.

## Exact decision equivalence

When $d=t^2-c^2$ is a unit, the coefficient is either $4\chi_N(t+c)$ if $d$ is a square modulo $N$, or zero otherwise. A Jacobi sign of $-1$ settles the zero case publicly. The remaining gate is precisely the usual quadratic residuosity decision for Jacobi-positive units. If $d$ is a proper nonunit, its gcd factors the semiprime and the coefficient can be evaluated locally. If $d=0\bmod N$, check $t-c$ and $t+c$ by gcd. Either a gcd factors the input or $t=\pm c\bmod N$, in which case the unique weighted root is $m=t/2$.

Conversely, given a Jacobi-positive unit $d$, take $t=(d+1)/2$ and $c=(d-1)/2$ modulo $N$. This gives $t^2-c^2=d$. If $c$ is zero then $d=1$. If $c$ is a proper nonunit its gcd factors the semiprime. Otherwise the coefficient is nonzero exactly for a square. These transformations take polynomial bit time.

The equivalence is for a single grouped phase coefficient. It does not supply the normalized complex sum to constant precision. In the root expression, that sum compresses all coefficient values to one hidden cosine. A native way to perform that compression is still needed.

## Checks

Independent enumeration verified 245,980 phase coefficients on three semiprimes, including the Jacobi-positive square and nonsquare cases. Run `python check_salie_phase_weights.py`.
