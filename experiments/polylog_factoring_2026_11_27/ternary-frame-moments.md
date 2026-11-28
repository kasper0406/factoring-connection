# Quartic moments of the 72 ternary Gram frames

For the mixed Gram $A=x^2+1,\ B=-Nx,\ C=N^2$, one can read all three quartic moments that identify a factor from a single scan of $1\le b<N$. Each step uses at most two integer square tests. For $\gcd(x,N)=1$ exactly three positive-coordinate rows pass, and each passing $b$ already gives a proper factor by $\gcd(b,N)$. This costs $O(N\,\mathrm{polylog}\,N)$ bit operations for $0\le x<N$. At $x=1$, one fixed weight-three CM coefficient gives the same factor. [Mixed Gram construction and count](mixed-ternary-embeddings.md) · [Exact coordinate scan and independent checks](check_ternary_frame_moments.py)

Assume $N=dr$ with distinct odd primes and $\bigl(\frac{-A}{N}\bigr)=-1$, naming $r$ the prime where $-A$ is a square. There are 24 coordinate frames and 48 other frames. The latter form one orbit under all signed coordinate permutations: this group has order 48, and any nonidentity element fixing both marked vectors would be a reflection through their plane. The possible primitive normals to reflection planes of the signed permutation group have lengths $1,\sqrt2$. The primitive normal to a noncoordinate frame has length $r$.

## The three traces

For each frame put
$$
 P_0=\sum_i v_i^4,\qquad P_1=\sum_i u_i v_i^3,\qquad P_2=\sum_i u_i^2v_i^2 .
$$
Subtract $24(N^4,-xN^3,x^2N^2)$ from the vector of full traces and divide by 48 **in the integers**. The three resulting integers $M_0,M_1,M_2$ are the values on a single noncoordinate frame by orbit invariance. They satisfy $\gcd(N,M_0,M_1,M_2)=d$.

Indeed $v=d w$ with $w$ primitive. Modulo $r$, $w\cdot w=0$, $u\cdot w=0$ and $u\cdot u=A\ne0$. Thus $u,w$ span $w^\perp$. If $M_0=M_1=0\bmod r$, then $(w_1^3,w_2^3,w_3^3)$ belongs to $(w^\perp)^\perp=\langle w\rangle$. All its nonzero coordinate positions therefore have $w_i^2=\lambda$ for the same $\lambda\ne0$. Isotropy forces $r$ to divide the number of such positions. Consequently $r=3$ and every position occurs. In that case $M_2=d^2\lambda A\ne0\bmod r$. The first two traces thus suffice once the factor 3 is excluded. One random linear combination $M_0+yM_1+zM_2$, for independent uniform coefficients modulo $N$, succeeds with probability at least $1-1/r$.

These are weighted Fourier coefficients of the genus-two theta series for the ternary sum of squares. Write a frame as the matrix with columns $u,v$, and grade its term by its Gram matrix. The harmonic quartic $H(y)=\sum_i y_i^4-\frac35(y\cdot y)^2$ on $\mathbb R^3$ gives all three weights by taking coefficients of degrees $0,1,2$ in $s$ in $H(v+s u)$. The radial part is explicit from the Gram. Each coefficient is pluriharmonic because every pair of column Laplacians on $H(v+s u)$ is a scalar multiple of $\Delta H$. In degrees below four, averaging a polynomial over signed coordinate permutations gives only constants and linear combinations of $u\cdot u,u\cdot v,v\cdot v$. These finite traces then depend only on the prescribed Gram and the known number 72. Pall gives the ordinary representation formula. Pluriharmonic theta weights are classical; the calculation here identifies which weights carry a proper gcd on this particular Gram. [Pall (1949), §3 and Theorem 4](https://doi.org/10.4153/CJM-1949-032-9) · [Kashiwara and Vergne (1978), §5](https://www.cmls.polytechnique.fr/perso/vergne/publications2/Remplacer/78onthesegalshaleweil.pdf) · [Ibukiyama on vector-valued theta constructions](https://doi.org/10.1080/10586458.2017.1353454)

## A row scan for the exact coefficients

The three columns $v,\ Nu+xv,\ u\times v$ form $N$ times an orthogonal matrix of determinant $-1$ after division by $N$. Each row $(b,t,c)$ satisfies
$$
b^2+t^2+c^2=N^2,\qquad t\equiv xb\bmod N .
$$
If $b=0\bmod N$, the row is an axis of length $N$. On any other row $c^2=-A b^2\bmod N$, so $d$ divides $b,t,c$ and $\gcd(b,N)=d$. This already recovers a factor from one row without constructing a full frame.

To compute the moments, scan positive $b<N$, try $t=xb\bmod N$ and its translate by $-N$, and retain $c=\sqrt{N^2-b^2-t^2}$ when it is a positive integer. Each retained row contributes exactly
$$
 (b^4,\ ab^3,\ a^2b^2),\qquad a=(t-xb)/N,
$$
to the normalized traces. The remaining argument gives the multiplicity, including rows with an axis in one column.

After division by $d$, a nonaxis row has norm squared $r^2$ and reduces modulo $r$ to one of the two isotropic lines $\ell=\langle(1,x,\pm s)\rangle$, $s^2=-A\bmod r$. Lift a generator to $l\in\mathbb Z^3$ with $l\cdot l=0\bmod r^2$, using a nonzero coordinate to correct it by a multiple of $r$. Put
$$
 L_\ell=\{z\in\mathbb Z^3:z\cdot l=0\bmod r\}+\mathbb Z(l/r).
$$
This is an integral lattice of covolume one. Its first summand has index $r$, and adjoining $l/r$ has index $r$ in the reverse direction. A positive definite integral unimodular ternary lattice is isometric to $\mathbb Z^3$, so it has exactly six unit vectors. For any primitive $y$ of squared norm $r^2$ reducing into $\ell$, write $y=hl+rz$, $r\nmid h$; expanding its norm gives $l\cdot z=0\bmod r$. Thus $y/r$ is one of its six unit vectors. Conversely every nonintegral unit vector yields such a row. The integral unit vectors of $L_\ell$ are exactly the signed coordinate directions perpendicular modulo $r$ to $l$. There are two if $r\mid x$ and zero otherwise; the first and third coordinates of the line never vanish.

An orthonormal row basis of $L_\ell$ with determinant $-1$, scaled by $N$, reconstructs a marked frame: read $v$ from column one and $u$ from the difference between columns two and $xv$, divided by $N$. Its third column is $u\times v$. Each nonaxis first row admits four such completions, corresponding to the oriented signed bases in its orthogonal coordinate plane. The two lines share no nonaxis row. Permuting row positions supplies a factor of three, and each row with $b,c>0$ represents four choices of the joint sign of $(b,t)$ and the sign of $c$. Hence $4\cdot3\cdot4=48$ occurrences per retained monomial, exactly the normalization above. Any integral row in a noncoordinate frame lies along its second coordinate in this row lattice and has $b=0$, contributing zero. The count with $b,c>0$ is three when $r\nmid x$ and two otherwise. This derivation is a direct ternary neighbor calculation, requiring no prime factor to perform the scan.

## A fixed classical coefficient when $x=1$

Suppose $\bigl(\frac{-2}{N}\bigr)=-1$. Fix $u=(1,1,0)$. The frame conditions give $v_1+v_2=-N$ and, with $Y=v_1-v_2$ and $Z=v_3$, reduce to $Y^2+2Z^2=N^2$. In the Euclidean ring $\mathbb Z[\sqrt{-2}]$ write $r=\pi\bar\pi$; the other prime $d$ is inert. The noncoordinate solutions satisfy $Y+Z\sqrt{-2}=\pm d\pi^2$ or its conjugate. In particular
$$
 2M_2=N^2+d^2(\Re\pi^2)^2 .
$$
The Hecke character $\psi((\alpha))=\alpha^2$ gives the normalized weight-three form $f=\eta(z)^2\eta(2z)\eta(4z)\eta(8z)^2$. The ideals of norm $N^2$ give
$$
 a_f(N^2)=d^2(\pi^4+r^2+\bar\pi^4)=8M_2-5N^2,\qquad
 \gcd(a_f(N^2),N)=d .
$$
Modulo the prime ideal $(\pi)$ of $\mathbb Z[\sqrt{-2}]$, only $\bar\pi^4$ survives inside the parentheses, proving the gcd. One can evaluate the coefficient without factors as half the sum of $Y^2-2Z^2$ over the signed solutions of $Y^2+2Z^2=N^2$, in $O(N)$ square tests. The eta product identifies an exact coefficient at a composite index; a polylogarithmic evaluation here would factor every modulus in this subfamily. The general composite-index computation of eigenforms already has factoring reductions in the literature. [The weight-three newform](https://www.lmfdb.org/ModularForm/GL2/Q/holomorphic/8/3/d/a/) · [Ahlgren, Ono and Penniston, weight-three eta products](https://uva.theopenscholar.com/files/ken-ono/files/065_8.pdf) · [Bach and Charles on composite indices](https://arxiv.org/abs/0708.1192)

## Checks

Run `python experiments/polylog_factoring_2026_11_27/check_ternary_frame_moments.py` from the repository root. Across 175 mixed Grams for eight semiprimes, the script found 510 positive rows. The first quartic missed 88 instances and the first two missed 32, all at $r=3$. The three traces returned the proper factor every time. On 24 Grams their exact values agreed with a separate enumeration of the full marked representations. Six independent binary theta coefficients agreed with the $x=1$ formula.

At $N=77,x=1$, the positive rows, written $(a,b,c)$, are $(-1,14,42),(0,42,49),(-1,63,42)$. The normalized traces are $(18903073,-252791,4165)$, and $a_f(77^2)=3675$; the gcd is 7. [Reproduce the checks](check_ternary_frame_moments.py)
