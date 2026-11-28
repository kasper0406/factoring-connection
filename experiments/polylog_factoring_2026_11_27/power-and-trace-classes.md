# Which prime residue classes does a public power certificate constrain?

A residue $a$ with $a^{2^h}\equiv D\bmod N$ certifies that $D$ is a $2^h$-th power at every prime dividing $N$. For any squarefree integer $D$ divisible by an odd prime, increasing $h$ adds no universal exclusions of prime residue classes beyond the quadratic character of $D$. It can change their densities. A different recurrence, $x\mapsto x^2-2$, can certify that each prime lies in just *two* dyadic classes. The second certificate has a strong restriction on the inputs where it can exist. The comparison below concerns congruence information implied by the existence of a certificate. Its numerical value can also be used in further computations. [Trace verifier and class search](cyclic_trace_classes.py) · [Root counts and controls](check_cyclic_trace_classes.py)

## Ordinary power roots and congruence classes

Fix $h\ge1$, $d=2^h$, a squarefree $D$ with an odd prime divisor, and a modulus $m$. Exclude primes dividing $2mD$. Write $F=\mathbf Q(\alpha)$ for $\alpha^d=D$. Put $E=\mathbf Q(\zeta_d)$ and $L=E(\alpha)$. The polynomial is Eisenstein at a prime of $E$ over any odd prime dividing $D$, since that prime is unramified in $E$ and $D$ has valuation one. Thus $[L:E]=d$ and the translations $\alpha\mapsto\zeta_d^t\alpha$ form a cyclic group of order $d$. Every unit automorphism of $E$ extends fixing $\alpha$, using the basis $1,\ldots,\alpha^{d-1}$. The full group is the affine group of translations and odd multipliers on $\mathbf Z/d$.

The commutator subgroup of this group consists of even translations: multipliers act on translations by their unit multiplier, and $-1$ generates all even translation differences. The same statement at $d=2$ gives the trivial subgroup. Hence the maximal abelian subfield of $L$ is $E(\sqrt D)$. Inside $F$, its intersection is precisely $\mathbf Q(\sqrt D)$ by taking the invariants of the unit automorphisms that fix $\alpha$. In particular $F\cap\mathbf Q(\zeta_m)$ is either rational or this quadratic field.

Let $r$ be a unit residue modulo $m$. There are positive-density primes $s\equiv r\bmod m$ such that $X^d-D$ has a root modulo $s$ exactly when the cyclotomic automorphism for $r$ fixes $F\cap\mathbf Q(\zeta_m)$. Indeed the image on $\mathbf Q(\zeta_m)$ of the subgroup fixing $\alpha$ in the compositum $L\mathbf Q(\zeta_m)$ fixes precisely that intersection. Chebotarev realizes its elements as Frobenius classes. Conjugates of such elements have the same restriction to the abelian cyclotomic field, and the Frobenius fixed-point criterion detects a root. Thus if $m$ contains the quadratic conductor, exactly the classes with quadratic sign $+1$ are possible. If the quadratic field is absent, every unit class is possible. The assertion fixes $D,d,m$ before considering primes. [Stevenhagen and Lenstra, *Chebotarëv and his density theorem*](https://pub.math.leidenuniv.nl/~lenstrahw/PUBLICATIONS/1996d/cheb.pdf)

For a dyadic modulus $m=2^k$ with $k\ge h$ one can also give the conditional density *within* each class. Eisenstein at the odd prime works over $\mathbf Q(\zeta_{2^k})$ too. With the Frobenius multiplier fixed at $r$, the translation $t\bmod d$ is uniform in the Chebotarev count. A fixed root requires $(r-1)u=-t\bmod d$ for some $u$. The density is $1/\gcd(r-1,d)$, positive in every class. Information from this density depends only on the truncated valuation of $r-1$. The hypotheses on $D$ matter. For instance a $2^h$-th root of $-1$ forces $s\equiv1\bmod 2^{h+1}$ by the order of the root.

## A compact cyclic certificate

Let $C_1(x)=x$ and $C_{2t}(x)=C_t(x)^2-2$ on dyadic indices. Fix $j\ge1$, $M=2^{j+2}$. For any odd prime $s$, the number of solutions of $C_{2^j}(x)=0\bmod s$ is exactly
$$
 \begin{cases} 2^j,&s\equiv1\ \text{or}\ -1\bmod M,\\
 0,&\text{otherwise}.\end{cases}
$$
To prove this, write $x=z+z^{-1}$ in a quadratic extension of $\mathbf F_s$. By induction $C_t(x)=z^t+z^{-t}$. Vanishing forces $z$ to have exact order $M$, and Frobenius must send $z$ to $z$ or $z^{-1}$ for its trace to lie in the base field. Conversely every primitive $M$-th root has base-field trace under either congruence. Pairing inverses gives $\varphi(M)/2=2^j$ distinct traces. This is the elementary cyclotomic trace behind Lucas sequences and the $p\pm1$ group-order methods. [Earlier Lucas and Chebyshev work](../../notes/trackA_analysis.md)

Given a candidate $x$, $j$ modular squarings and a gcd can therefore give a factor or certify that **both** factors belong to $\{1,-1\}\bmod M$. If the full certificate holds, the list algorithm completes the search at its stated cost. Its validity holds for every $N$ in the balance promise. For a randomly chosen $x\bmod N$, conditional on both distinct primes satisfying the congruence, CRT gives exactly $4^j$ successful residues. The certificate can only exist on inputs whose two primes meet the congruence.

For a predetermined list of small integer traces there is also a coverage bound *without* assuming randomness of the trace values. Take a prime uniformly from $[U,2U]$. For a given integer $x$, the nonzero integer $C_{2^j}(x)$ has $O(M\log(2+|x|)/\log U)$ possible distinct prime divisors there, by its size. It is nonzero over the integers since $C_{2^j}(x)=0$ would require an integer square equal to two at the final step. Independently, the two progressions contain $O(1+U/M)$ integers in the interval. Dividing the smaller count by $\asymp U/\log U$ gives a bound
$$
 O\left(\min\left\{\frac{M\log(2+|x|)}{U},\
                       \frac{\log U}{U}+\frac{\log U}{M}\right\}\right).
$$
Consequently any list of polynomially many predetermined $(x,j)$ of polynomial bit length in $\log U$ misses both primes on a $1-U^{-1/2+o(1)}$ fraction of pairs. The count also covers choosing dynamically among the same list. It does not address choosing new values from the specific integer $N$.

The checker verified all 2,416 prime/level root counts for odd primes below 2,000 and $1\le j\le8$. Four planted products exercised the full certificate and two-class completion. At $N=5662847\cdot6803327$, $x=18$ at $j=5$ certified the classes $\pm1\bmod128$. The fixture primes were selected from the trace value. On 64 independent products at each of 16, 24 and 32 bits per prime, scanning public $x=3,\ldots,34$ through those respective levels gave zero splits or full certificates. [Checks](check_cyclic_trace_classes.py) · [Native quadratic character rows](cfrac-character-rows.md) · [An input-specific cyclic cubic](shanks-split-seeds.md)
