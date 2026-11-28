# A prescribed three-square class with the full factorization of its squared radius

Let $N$ be the product of two distinct odd primes and $x$ a short public parameter with Jacobi symbol $(-(x^2+1)/N)=-1$. One can reduce factoring $N$ to integer-point decisions on a single sphere **and give the evaluator the complete prime factorization of its squared radius**. Writing $A=x^2+1$, that radius squared is $M=A(p_1p_2)^2$. The two auxiliary primes are distinct, comparable in size and coprime to $2AN$. For fixed or polylogarithmic $x$, the randomized preprocessing has expected polynomial bit complexity with a polynomial-time primality test. It uses only the prime number theorem modulo four. [Construction](factored_gaussian_packing.py) · [Executable checks](check_factored_gaussian_packing.py) · [Short-parameter domain](short-jacobi-grams.md)

The earlier packing produced a near-square $M$ of unknown factorization. The Gaussian identity below puts $M$ in a specified squareclass. The extra step here constructs the two prime factors of its conductor $T=p_1p_2$. Thus even the exact total numbers of ordinary and primitive representations of $M$ can accompany every individual decision. The construction supplies the factorization of $M$. It gives no general factorization guarantee for the congruence modulus. [Original packing](packed-ternary-decision.md) · [Fixed-field packing and global samplers](fixed-field-packing.md)

## Finding the two primes

Set $J=8(x+1)$, $H=(1+\lfloor\sqrt{NJ}\rfloor)\lceil\log_2(N+1)\rceil$. Seek Gaussian primes $\pi=a+ib$, $\rho=c+ie$ with

$$
0<b<a,\quad 0<e<c,\quad H^2/2<|\pi|^2,|\rho|^2\le H^2 .
$$

Require that their coordinate pairs span the same projective line modulo $N$, have the same parity in the real coordinate, and have $b/a,e/c$ in the same one of $J$ equal intervals in $[0,1)$. For two distinct such points, define

$$
B=ac+be,\qquad C=|ae-bc|=Nt,\qquad T=B^2+C^2=|\pi|^2|\rho|^2 .
$$

The primitive vectors have different slopes. The congruence and parity imply $t$ is a positive even integer. Also $C/B<1/J$. An odd prime $1\bmod4$ has exactly one prime vector in the chosen octant, so the two norms are distinct primes. Exclude primes dividing $2AN$.

**Why a random line works.** The prime number theorem in the fixed progression $1\bmod4$ and the two-square theorem give $\gg H^2/\log H$ eligible prime vectors across all the lines, before accounting for the parity and slope labels. Excluding prime factors of $2AN$ removes at most $O(\log(2AN))$ vectors. Each primitive vector belongs to exactly one line. There are $N\prod_{s\mid N}(1+1/s)<2N$ lines for the semiprimes here. A line containing no usable pair accounts for at most $2J$ vectors. For $x<N$, $\log H=O(\log N)$, and the full prime count eventually dominates $4NJ$.

Every line contains at most $O(1+H^2/N)$ primitive vectors in the disk, uniformly even for a very narrow lattice. To see this, Gauss-reduce its inverse-image lattice of determinant $N$ to a basis $u,v$. Write $s=\|u\|^2$. For $i u+j v$ in the disk,

$$
(s i+j\,u\cdot v)^2+(jN)^2\le sH^2 .
$$

When $j=0$, only $i=\pm1$ can be primitive. For nonzero $j$ there are $O(H\sqrt s/N)$ choices of $j$ and $O(1+H/\sqrt s)$ choices of $i$ for each. Since $s=O(N)$ for a reduced basis, the claim follows. These same bounds give an explicit output enumeration with polynomial bit work per pair of coefficients.

It follows by dividing the remaining prime count by the uniform bound that a fraction $\gg1/\log H$ of all lines contain a pair. This only uses the *total* prime count. We do not require information about primes in any specified progression with modulus $N$. To sample a uniform line, draw a uniform pair modulo $N$ until it is unimodular. Every line has exactly $\varphi(N)$ unimodular generators. From a generator $(a,b)$ let $g=\gcd(a,N)$, $m=N/g$. A basis is $(g,b(a/g)^{-1}\bmod m),(0,m)$, taking the residue to be zero if $m=1$. This also works for nonunit first coordinates.

At this height, $p_i=O(N(x+1)\log^2 N)$ and $T=O(N^2(x+1)^2\log^4 N)$. Positivity of the packing coefficient below requires $B>(x+\sqrt A)C$, giving $T>(x+\sqrt A)^2N^2t^2$. The construction comes within powers of $\log N$ of this size bound. For every sufficiently large $N$ the initial height suffices. The implementation widens its search after $16\operatorname{bitlength}(H)^2$ unsuccessful trials to cover small inputs. The probability of a widening is $\exp(-\Omega(\log H))$ at any sufficiently large height. Repeating the bound at successive doubled heights gives expected polynomial work. With a strict output-size budget one can retry at the initial height in the asymptotic regime.

## Using the known factorization

For $C=Nt$ set

$$
P=B^2-C^2-2xBC,\qquad Q=2ABt.
$$

The identity $(x+i)(B-iC)^2=(Px+QN)+iP$ gives $(Px+QN)^2+P^2=AT^2$. The angular condition ensures $P>0$. Squarefreeness of $T$ implies $\gcd(B,C)=1$. Thus $P$ is a unit at all primes dividing $2Bt$. If a prime dividing $A$ also divided $P$, then $P=(B-xC)^2-AC^2$ would force $B=xC$ there and hence $T=0$ there. This contradicts the exclusions. Therefore $\gcd(P,Q)=1$, and also $\gcd(T,Q)=1$.

For a possible first vector $\|u\|^2=A$, put $W=Pu-Qv$. The points of squared norm $M$ with $W=Pu\bmod Q$ are precisely the partners with $\|v\|^2=N^2$, $u\cdot v=-Nx$. Indeed

$$
Q(\|v\|^2-N^2)=2P(u\cdot v+xN),\qquad
ABt\mid u\cdot v+xN,\qquad |u\cdot v+xN|\le B/t .
$$

The inequality follows from $\|v\|\le(P+T)\sqrt A/Q=(B-xC)/(t\sqrt A)$. Since $At^2>1$, the right side of the norm identity vanishes. Every resulting $W$ is primitive. A prime divisor of its content divides $Q$ or the minors of $(u,v)$, by $u\times W=-Q(u\times v)$. In the first case it would divide the primitive $u$, and in the second it would divide $N$ and also $M$, which are coprime. Refining the class modulo $Q\ell$ therefore gives the original polynomial menu of existence decisions. The factorization of $A$ can be supplied by trial division in time polynomial in $x+1$.

For prime factors $s=p_1,p_2$ write $\chi_s=(\frac{-A}{s})$. The full sphere has

$$
r_3^*(M)=r_3^*(A)\prod_s(s-\chi_s),\qquad
r_3(M)=r_3(A)\prod_s(1+s-\chi_s).
$$

The base counts are enumerable in polynomial time for short $x$. An exact uniform sampler is available even on the unrestricted sphere. Relative to the entire primitive sphere the chance of one of the 48 factor-revealing points is exactly $48/r_3^*(M)$. Known counts and fast global sampling leave the individual residue question. [Counting identities and samplers](fixed-field-packing.md)

## An occupied and an empty class with the same supplied arithmetic

At $N=77,x=1$ a seeded run chose Gaussian primes $(237,10)$ and $(205,8)$. Their norms are $56\,269$ and $42\,089$. It gave $B=48\,665,t=2,Q=389\,320$ and

$$
M=2(56\,269\cdot42\,089)^2=11\,217\,746\,060\,351\,790\,962.
$$

The primitive count is $28\,419\,501\,120$ and the ordinary count is $28\,420\,681\,428$. At the same modulus $17Q=6\,618\,440$ the following two queries reduce to the same first vector $(1,1,0)$. They also lie in the same proper spinor genus. Their cosets coincide locally at the primes of $Q$. At 17 their common norm is a unit, and a proper local rotation of square spinor norm carries one residue to the other. [Local rotation proof](ternary-coset-spinor-genus.md)

| Residue triple | Number of integral points | Result of direct enumeration |
|---|---:|---|
| $(2\,555\,529,1\,776\,889,3\,503\,880)$ | 1 | The point gives $v=(-14,-63,42)$ and $\gcd(v)=7$. |
| $(4\,112\,809,219\,609,5\,061\,160)$ | 0 | The three-square congruence still holds. |

The checker constructs complete squared-radius factorizations on twelve mixed challenges, including a 388-bit $M$ from an 81-bit $N$. It independently checks line enumeration on 2,016 directions. For $N=15$ a full reduction recovered 3 with 132 direct integer-point decisions, including 60 empty classes. The examples through 64-bit auxiliary primes use deterministic Miller–Rabin in its proved range. Above that range the default Python backend certifies primes by Lucas's criterion after factoring the auxiliary candidate minus one. That certificate generator has no polynomial worst-case guarantee. Substituting a polynomial-time primality test gives the stated theoretical bound. [Checks](check_factored_gaussian_packing.py) · [Agrawal, Kayal and Saxena](https://doi.org/10.4007/annals.2004.160.781)

Sardari studies the least *prime* represented by binary forms of fundamental discriminant. Zaman gives bounds for primes and almost primes represented by positive definite binary forms. Here two already identified primes suffice in a projective line chosen at random, and their product is represented by $B^2+N^2t^2$ at the required angle. Earlier repository searches also used Gauss-reduced two-dimensional lattices for norm relations. The application here is to supply all the auxiliary norm factors to the evaluator. [Sardari](https://ems.press/journals/jems/articles/17390) · [Zaman](https://doi.org/10.1093/qmath/hay028) · [Earlier norm-lattice calculation](check_near_square_prime_norms.py)
