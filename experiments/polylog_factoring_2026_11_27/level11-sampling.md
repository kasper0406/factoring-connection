# Uniform sampling from the level-11 quaternary norm form

For every positive integer $n$, there is an exact uniform sampler for the signed ordered integer solutions of

$$
Q(x,y,z,w)=x^2+xy+3y^2+z^2+zw+3w^2=n .
$$

It uses $n$ alone and takes expected time polynomial in its bit length. A second sampler draws uniformly from the primitive solutions, whose coordinates have gcd one. This shell is nonempty precisely when $121\nmid n$. The output is one representation. The number of representations and the prime factors of $n$ remain separate computational questions. [Implementation](level11_sampler.py)

The earlier level-11 experiment computed the theta identities for this form and its two companion lattices. The earlier four-square sampler used the Hurwitz order, where every right ideal is principal. The discriminant-11 order has two ideal classes. Here the theta identities control the cost of rejecting nonprincipal ideals and the weights needed to sample all possible coordinate gcds. The published SQIsign2D–West algorithm already samples uniform primitive ideals of prescribed norm. Its §3.3 proves the projective-line correspondence used below. [Level-11 theta identities](../../notes/nativemod_results.md#L466-L502) · [Four-square sampler](four-square-incidence.md) · [SQIsign2D–West, §3.3](https://www.math.u-bordeaux.fr/~damienrobert/pro/publications/articles/sqisign2d.pdf#page=14)

## Sampling a primitive ideal and recognizing its generators

Put $\omega=(1+\sqrt{-11})/2$. The order $\mathcal O=\mathbb Z[\omega]+\mathbb Z[\omega]j$, with $j^2=-1$ and $j\omega=(1-\omega)j$, has reduced norm $Q$ in the basis $1,\omega,j,\omega j$. Its four units are $\pm1,\pm j$. First suppose $11\nmid n$, so the order splits modulo $n$. The construction of the splitting below uses no hidden primes.

Draw $u,v\bmod n$ uniformly until $\gcd(u,v,n)=1$. Their span is a uniform line in $\mathbb P^1(\mathbb Z/n)$, since each line has the same number of generators. Its matrices with both columns in that line form a right ideal modulo $n$. Lift it to an integer lattice $I$ in $\mathcal O$ by adjoining $n\mathcal O$. Its index is $n^2$, and every element has norm divisible by $n$ because its determinant vanishes modulo $n$.

For $\alpha\in I$ of norm $n$, the lattice $\alpha\mathcal O$ already has index $n^2$, hence equals $I$. Its generators of that norm are exactly the four unit associates. They are primitive because the line is a direct summand at each prime dividing $n$. Conversely, a primitive norm-$n$ element maps to a matrix with at least one unit entry locally at each $r^e\parallel n$. Since its determinant vanishes modulo $r^e$, its columns span a unique direct-summand line there. CRT gives the unique global line. Thus returning a uniform generator whenever the list is nonempty samples every primitive point with equal probability.

The algorithm computes a basis of $I$ by Hermite normal form and applies exact four-dimensional LLL. For its reduced row basis $B$, enumerate the integer coefficients $c$ of $x=cB$ satisfying

$$
|c_j|\le\sqrt{2n}\,\|(B^{-1})_{\cdot j}\|_2 .
$$

This includes every point of norm $n$ since $Q(x)\ge\|x\|_2^2/2$. The box has bounded size independently of $n$. Indeed $Q(x)\le4\|x\|_2^2$, so each nonzero ideal vector has length at least $\sqrt n/2$. The LLL product bound and $|\det B|=n^2$ then give every basis row length $O(\sqrt n)$. Cofactors give $\|(B^{-1})_{\cdot j}\|_2=O(1/\sqrt n)$. HNF, LLL and enumeration therefore take polynomial bit time.

## Why at least a third of the ideals succeed

The three previously verified theta identities imply $3r_1+2r_2=12\sigma$ and $3r_2+2r_3=12\sigma$ away from 11. Here $r_1$ counts $Q$, and $r_2,r_3$ count the two companion forms. Remove common square divisors from all three counts by Möbius inversion. The resulting $u_i(n)$ count primitive points. The same operation converts $\sigma(n)$ into $\psi(n)=n\prod_{r\mid n}(1+1/r)$, the number of projective lines. Elimination gives

$$
9u_1(n)=12\psi(n)+4u_3(n).
$$

Consequently $u_1(n)/(4\psi(n))\ge 1/3$ for every $n$ coprime to 11, including powers of two and nonsquarefree inputs. Möbius inversion is only used in the proof. The sampler needs at most three ideal draws on average. [Earlier theta identities](../../notes/nativemod_results.md#L483-L502)

## Splitting the order using the public modulus

On the odd part $s>1$ of $n$, find $a^2+b^2=-11\bmod s$. Set $j=\left(\begin{smallmatrix}0&1\\-1&0\end{smallmatrix}\right)$ and $\sqrt{-11}=\left(\begin{smallmatrix}a&b\\b&-a\end{smallmatrix}\right)$. These give an integral coordinate isomorphism modulo $s$, since 2 and 11 are units there. To find $a,b$, draw $x,y\bmod s$ and set $t=-(x^2+y^2)/11\bmod s$ in $[0,s)$. If $t$ is a unit and a verified sum $z^2+w^2=t$ is found, the quotient $(x+iy)/(z+iw)$ has norm $-11\bmod s$. The case $t=1$ works immediately. For primes $t\equiv1\bmod4$ the implementation finds a square root of $-1$ with probability $1/2$ and uses Cornacchia's algorithm. Every accepted answer is checked by exact arithmetic.

Each unit $t$ has $s\prod_{r\mid s}(1-(\frac{-1}{r})/r)\ge\varphi(s)$ preimages. The prime number theorem modulo 4 gives $\gg s/\log s$ eligible primes for sufficiently large $s$, after omitting divisors of $s$. Thus the expected number of trials is polynomial in $\log s$. This is the Gaussian norm-ratio method of Pollack and Treviño with the right-hand side changed to $-11$. [Pollack and Treviño, §4](https://math.colgate.edu/~integers/sjs15/sjs15.pdf)

For $2^e\parallel n$, use $W=\left(\begin{smallmatrix}0&-3\\1&1\end{smallmatrix}\right)$ for multiplication by $\omega$ on $(1,\omega)$ and $C=\left(\begin{smallmatrix}1&1\\0&-1\end{smallmatrix}\right)$ for conjugation. Lift $b=0\bmod2$ to a root of $b^2+b+4=0\bmod 2^e$ one bit at a time. Then $J=(bI+W)C$ squares to $-I$ and satisfies $JW=(I-W)J$. Modulo 2 the algebra generated by $W$ is $\mathbb F_4$ and $C$ acts by its nonidentity automorphism. Thus $I,W,J,WJ$ form a basis of the full matrix ring modulo 2, hence modulo $2^e$. Combine the two coordinate splittings by CRT at the public coprime moduli $s,2^e$.

## Full shells and the prime 11

For $11\nmid n$, write a point as $g x$ with $x$ primitive of norm $m=n/g^2$. Propose each $g\le\sqrt n$ with mass proportional to $1/g^2$, retaining it when $g^2\mid n$. This can be done by drawing a dyadic block $[2^k,2^{k+1})$ with mass proportional to $2^{-k}$, drawing uniformly within it and accepting with probability $2^{2k}/g^2$. The implementation rejects proposals above $\sqrt n$.

Write $\eta(m)=\psi(m)/m$. Compute $h=n/\gcd(n,m^{\operatorname{bitlength}(n)})$ by modular powering. Draw a jointly primitive pair modulo $h$ and require its first coordinate to be a unit. Locally the probability is $1/(1+1/r)$ at each prime present in $n$ and absent from $m$. This inserts the factor $\eta(m)/\eta(n)$. Make **one** uniform ideal draw at $m$, restarting the whole proposal on an empty ideal. Each individual primitive generator then has unconditioned probability $1/(4\psi(m))$ in that draw. Multiplying the three weights gives a per-point mass proportional to $1/(4n\eta(n))$, independent of $g$. The choice $g=1$ alone gives total acceptance at least $1/6$. [Earlier content correction](four-square-incidence.md#sampling-when-kne2)

For general $n=11^e m$, remove the public power of 11 and multiply the sampled element by $(\sqrt{-11})^e$. Multiplication by $\sqrt{-11}$ sends each binary pair $(x,y)$ to $T(x,y)=(-x-6y,2x+y)$, with $T^2=-11I$. This maps the entire norm-$m$ shell bijectively to the norm-$11m$ shell in one step. Indeed $4Q=(2x+y)^2+(2z+w)^2+11(y^2+w^2)$. If $11\mid Q$, the fact that $-1$ is a nonsquare modulo 11 forces $2x+y$ and $2z+w$ to vanish modulo 11. Both applications of $T^{-1}=-T/11$ are integral. For $e=1$ primitivity is preserved, while for $e\ge2$ every output has all coordinates divisible by 11.

## Finite checks

The exhaustive diagnostic checked all 1,826 projective lines at 34 unramified norms, including powers of two and repeated odd primes. It compared the full outputs to 4,868 points obtained by independent direct enumeration. All 58 content strata had exactly equal per-point sampler mass at their respective norms. At $n=81$, the three contents $1,3,9$ have 264, 24 and 4 points. Twelve ramified norms supplied another 172 points, all matched by the multiplication bijection. Twenty-four seeded draws through 254 bits checked both public samplers, and a further 256 full draws at 81 exercised all three contents. These checks verify exact finite identities and output validity. [Runnable diagnostic](check_level11_sampler.py)

## Can the ideal test distinguish the hidden primes?

The sampler supplies a public yes/no test as well as points. Fix its matrix coordinates and let $C_n(L)$ be 1 if the ideal of a supplied projective line $L$ has a norm-$n$ generator. For $N=pq$ with distinct primes away from 11, a uniform line has independent uniform components $L_p,L_q$. The computation sees $N$ and $L$ alone. We use the factors to analyze its output.

Put $\pi=3/5$ and $\lambda_r=a_r/(r+1)$, where $a_r$ is the prime coefficient of $\eta(z)^2\eta(11z)^2$, the weight-two newform of level 11. In the coordinate splittings induced from $N$,
$$
\mathbb E[C_N\mid L_p]=\pi+\lambda_q(C_p(L_p)-\pi).
$$
Indeed the ideal at the specified $p$-line has $p$-primary quotient in the order. At $q$ it is the full local order. Its $q+1$ subideals are exactly its intersections with the ideals of the $q$-lines. Passing to its ideal-class representative preserves the classes of these subideals. The normalized two-by-two Brandt matrix at $q$ has stationary masses $(3/5,2/5)$, from the inverse unit-group sizes. Its eigenvalues are 1 and $\lambda_q$ by the level-11 identities. Its principal-column entries are therefore $\pi+(1-\pi)\lambda_q$ and $\pi-\pi\lambda_q$, giving the formula. Interchanging the primes gives the other marginal. [Brandt matrices and the level-11 eigenform](../../notes/nativemod_results.md#4-q4--brandt-matrices-and-quaternary-theta-series)

For any statistic $f(L_p)$ of standard deviation one, $|\operatorname{Cov}(C_N,f)|\le |\lambda_q|/2\le\sqrt q/(q+1)$ by the Hasse bound for the level-11 elliptic curve. If $a_q=0$, conditioning on acceptance leaves the $p$-line exactly uniform. For $551=19\cdot29$, both coefficients vanish. The computation accepts 360 of the 600 lines. Every row accepts 18 of 30 and every column accepts 12 of 20. The resulting $20\times30$ zero-one matrix has rank 20 over the rationals, as checked by exact elimination. A single polynomial equality on the line coordinates modulo $N$ would give a rectangle of accepting CRT pairs and a matrix of rank at most one. The integer lattice test computes a predicate outside that class. Either prime component by itself leaves its probability unchanged.

**A factor-producing test.** Choose a reference line $T$ using $N$ and the chosen splitting before the sample. Draw a uniform principal line $L$, for example by rejection with the public predicate. Compute the gcd of $N$ and the determinant of the two primitive pairs representing $T$ and $L$. This returns a proper divisor exactly when one prime component agrees and the other differs. With $\mu=\mathbb E C_N=\pi+(1-\pi)\lambda_p\lambda_q$, its exact success probability is
$$
\frac{1}{\mu}\left[
\frac{\pi+\lambda_q(C_p(T_p)-\pi)}{p+1}
+\frac{\pi+\lambda_p(C_q(T_q)-\pi)}{q+1}
-\frac{2 C_N(T)}{(p+1)(q+1)}
\right].
$$
The subtraction removes the double hit, which returns $N$. Uniformly for $p,q\in[R,2R]$, the probability is $(1/(p+1)+1/(q+1))(1+O(R^{-1/2}))$. Even selecting a reference line using its publicly computable class gives order $\sqrt N$ expected trials on this family. At 551 a nonprincipal reference gives probability $1/12$, and a principal reference gives $7/90$. The test verified the marginal and collision formulas on all 4,032 cells of 21 semiprimes. [Public predicate](level11_sampler.py) · [Independent enumeration and gcd checks](check_level11_marginals.py)
