# What a correlated cubic relation actually gives

The conductor solver records which combination of Shanks fields has small conductor. Keeping the modular roots during that calculation yields more. A nontrivial combination of conductor $C$ gives a modular cube root for **one** rational Shanks parameter of norm at most $9C$. A globally trivial combination can give a factor directly. On the three earlier examples the second calculation recovered a proper gcd from 90 of the 135 seed cases previously recorded as trivial. [Conductor solver](shanks-conductor-relations.md) · [Witness fusion and tests](shanks_witness_fusion.py)

## Compress the witness along with the field

Let $E=\mathbf Q(\omega)$, $\omega^2+\omega+1=0$, and $\gamma(\beta)=\beta/\bar\beta$. A public root $x_i$ for a good rational parameter $a_i/b_i\bmod N$ supplies
$$
 u_i=(x_i-\omega)/(x_i-\bar\omega),\qquad
 u_i^3=\gamma(a_i-3b_i\omega)\pmod N.
$$
Given a relation $e_i\in\{0,1,2\}$, put $U=\prod_i u_i^{e_i}$ modulo $N$ and $G=\prod_i\gamma(a_i-3b_i\omega)^{e_i}$ exactly in $E$. Suppose the conductor solver reports tame conductor $g$ together with its factors. At each $r\mid g$ pick a generator $\pi_r$ of one of the two primes above $r$. Select $\pi_r$ if $v_{\pi_r}(G)=1\bmod3$ and $\bar\pi_r$ otherwise. The product $\beta_0$ of the selected elements has norm $g$. The quotient $G/\gamma(\beta_0)$ has every valuation divisible by three, including the unique prime above 3, since conjugation negates the valuations of $G$. In the Euclidean ring $\mathbf Z[\omega]$ it is therefore a global cube times a unit. Units are $\mu_6$ and their cubes are $\{1,-1\}$. Rotating $\beta_0$ by $1,\omega,\omega^2$ exhausts the three unit classes in the ratio $\gamma(\beta_0)$. For one rotation $\beta$ we obtain
$$
 G=\gamma(\beta)h^3,\quad h\in E,\quad h\bar h=1.
$$
Thus $U/h\bmod N$ witnesses the short representative. All divisions are at units over every prime of $N$ by the valid input certificates and the coprimality of the resulting conductor with $N$. The rotation search also certifies its answer by exact cubing. To extract $h$ from the norm-one number $z=h^3$, seek a rational root $T$ of $T^3-3T=z+\bar z$. Then $h=(z+T)/(T^2-1)$, with $z=\pm1$ handled directly. Fixed-degree rational polynomial factorization is polynomial in bit length. Small split primes can be oriented by trial modulo $r$ and by enumerating Eisenstein vectors of norm $r$. The cost is polynomial in the original sizes and the *small prime cutoff*, with no large integer factorization. [Lenstra, Lenstra and Lovász, rational polynomial factorization](https://doi.org/10.1007/BF01457454)

Write $\beta=A+B\omega$. If $B\ne0$ the reduced parameter $(a,b)=(3A,-B)/\gcd(3A,B)$, up to overall sign, has
$$
 Q(a,b)=a^2+3ab+9b^2\le 9g\le9C.
$$
The primary certificate is the norm-one element $z'=U/h$ with $(z')^3=\gamma(\beta)\bmod N$. Where $z'-1$ is a unit, $x=(z'\bar\omega-\omega)/(z'-1)$ is a rational residue and solves the homogeneous Shanks equation. A global $\omega$ rotation of $z'$ changes the choice of root. A denominator nonunit at just one input prime exposes a factor. For balanced inputs with $C$ polylogarithmic and $B\ne0$, $b$ is automatically a unit for sufficiently large $N$ by its size. The resulting parameter is then a good ordinary Shanks parameter. For a trivial field one can take $\beta=1$ and use the gcd test below.

On $N=1247$, the original certificates at seeds 63 and 66 fuse to the short parameter $-2$ and root 1152. At $N=2813$, seeds 78 and 111 give the same parameter and root 2200. Its norm and conductor both equal 7. These are roots constructed from the two inputs to the relation, without supplying either factor of $N$. A multirow search can still be an effective way to discover this root. The equivalence identifies what it must produce. [Implementation](shanks_witness_fusion.py)

## A trivial relation can split the input

For $N=pq$, $p,q>3$ distinct, each local norm-one group in the Eisenstein algebra has exactly three cube roots of 1. In the split case it is a copy of $\mathbf F_p^\times$, and in the inert case it is the norm-one subgroup of $\mathbf F_{p^2}^\times$. In both cases these roots are $1,\omega,\omega^2$. Thus a modular cube root of 1 has two labels. If the labels differ, comparing its two coefficients to a global candidate $\omega^j$ by gcd recovers a prime. If $G=h^3$ exactly, $U/h$ supplies such a root.

For independent uniform good seeds, condition on their resulting public parameters and choose any **nonzero** relation using only these parameters. For each fixed accepted parameter there are nine good seed preimages, three independently at each prime. Multiplication by a nonzero $e_i\bmod3$ permutes the local labels. Hence if the relation is globally trivial, $U/h$ has uniform labels in $\mu_3^2$ and the gcd splits $N$ with probability exactly $2/3$. This argument also applies to one uniformly sampled seed when a lift and relation are chosen from its parameter alone. Correlated generation can change the label law completely. The underlying extraction by disagreement between roots is the standard congruence-of-powers mechanism. The earlier hidden-prime cube-extractor note used an extraction oracle on an asymmetric promise. Here the actual supplied split witnesses and the exact cube of a global product supply the comparison. [Earlier cube-root extraction](bounded-quadratic-roots.md) · [Root-fiber count](shanks-split-seeds.md)

At $N=713$, the public seed $10$ gives the shortest rational lift $1/6$. This cubic is already split over $\mathbf Q$. Two of its rational roots are $-3/2$ and $-1/3$. The first agrees with 10 modulo 23 and the second modulo 31. Directly computing $\gcd(2\cdot10+3,713)=23$ and $\gcd(3\cdot10+1,713)=31$ exposes both. Over all previously recorded trivial seed cases at $N=713,1247,2813$, 90 gave a proper gcd and 45 had diagonal labels.

**How frequent is the one-row event?** For balanced $p,q$ tending to infinity, the probability that the shortest rational lift is globally split is $\Theta(N^{-2/3})$ for a uniform seed, and the same order holds for factor recovery from that event. To prove this, homogenize the rational Shanks map at a primitive direction $X/Y$:
$$
 F=X^3-3XY^2-Y^3,\qquad D=XY(X+Y),\qquad
 F-3D\omega=(X-Y\omega)^3,\qquad Q(F,D)=(X^2+XY+Y^2)^3.
$$
The pair $(F,D)$ is primitive, since reducing $F$ modulo $X,Y,X+Y$ gives a unit at every prime of $D$. A globally split ordinary parameter has at least one rational root, so it comes from this map. The shortest-lift bound $Q\le3N$ requires $|X|,|Y|=O(N^{1/6})$. This gives $O(N^{1/3})$ candidate parameters and nine modular seeds per parameter.

For the lower bound take primitive $1\le X,Y\le cN^{1/6}$ with a sufficiently small constant $c$. There are $\Theta(N^{1/3})$ such pairs. All the seed unit checks pass for large balanced inputs by size. Their $(F,D)$ vectors have squared length strictly smaller than the Eisenstein lift lattice covolume $3\sqrt3N/2$. Any independent lattice vector of at most that length would span an area smaller than the covolume. Since $(F,D)$ is primitive, it is the unique shortest lift up to sign. The parameters are distinct modulo $N$: if their shortest vectors shared a lattice they would agree up to sign, and the three global roots of a positive $X/Y$ have exactly one positive member. Thus each parameter has a nine-element modular seed fiber, six of which give proper factors. In particular the deliberately chosen *short rational seed itself* is diagonal and yields no split through this comparison. These controls account for both the observed strict gain and its cost on large uniform inputs.

## Which correlated root transformations only reproduce the original field?

Over $E$, use $u=(x-\omega)/(x-\bar\omega)$ and parameter coordinate $v=u^3$. Suppose a rational transformation of the known root produces $H(u)$, and its new cubed value depends rationally on the public parameter alone, so $H(u)^3\in E(v)$. The quotient $H(\omega u)/H(u)$ has cube 1 in a field and is identically $\omega^j$ for a fixed $j$. Therefore $H(u)/u^j$ is invariant under $u\mapsto\omega u$. Since the invariant field is $E(v)$, **all** such transformations are
$$
 H(u)=u^j R(v),\qquad H(u)^3=v^jR(v)^3,\qquad j\in\{0,1,2\}.
$$
This covers rational maps of any degree, including Möbius maps. The only resulting global field is the original one or the trivial field. A pointwise-valid identity of such maps whose product cubed class vanishes also gives a diagonal cube-root label: the resulting ratio is a single constant in $\mu_3$ in the function field. It cannot force different labels after good specialization at the two primes. For several formally independent roots $u_i^3=v_i$, the same argument applied to each separate deck transformation gives $H=R(v_1,\ldots,v_m)\prod u_i^{j_i}$. These statements concern transformations descending to the parameter field. Input-dependent integer lifts, reductions by $N$, and restricted specializations need not descend. [Kummer theory](https://www.jmilne.org/math/CourseNotes/CFT.pdf) · [Hoshi on simplest cubic field isomorphisms](https://arxiv.org/abs/0810.3374)

After these tests, a productive correlation would have to give either a short root certificate for a nontrivial small-conductor character or a globally trivial relation with a detectable mixed root. The open question is how to force one of those outcomes from $N$ more often, while tracking the actual independence of the two local labels.
