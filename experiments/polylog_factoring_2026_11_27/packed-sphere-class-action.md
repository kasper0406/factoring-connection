# Which ideal classes return to the packed sphere residues?

For the mixed ternary Gram construction, **exactly six proper ideal classes** send the public basepoint to a packed point, up to proper integral rotation. Two classes give coordinate frames. The four useful classes have reduced leading coefficient equal to the hidden split prime or twice that prime. Supplied class powers can be evaluated in polynomial bit time. This gives an exact target set for a walk and a way to recognize the hidden factor on a successful return. [Class action](packed_sphere_class_action.py) · [Executable checks](check_packed_sphere_class_action.py)

Let $N=dr$ be a product of distinct odd primes and take a nonnegative $x$ with $A=x^2+1$ and Jacobi symbol $(-A/N)=-1$. Label the prime at which $-A$ is a square by $r$. Let $Q$ be the least power of 7 at least $8(x+1)(N+1)$, and put
$$
H=QN+x,\qquad M=H^2+1,\qquad W_0=(H,1,0).
$$
A return means $W\equiv u\pmod Q$ for **some** integral $u$ of squared length $A$. This union is invariant under the 24 proper signed permutations. A query prescribing one particular $u$, or a further residue modulo another prime, can retain fewer points. Packing identifies the union bijectively with the 72 frames having Gram $(A,-Nx,N^2)$ by $W=u-Qv$. There are 24 coordinate frames and 48 informative frames. Every resulting point is primitive. If a prime divides $W$, it divides every coordinate of $W\times v=u\times v$, whose squared length is $N^2$. It must divide both $N$ and $M$, contradicting $\gcd(N,M)=\gcd(N,A)=1$. [Mixed frame count](mixed-ternary-embeddings.md) · [Packing bijection](packed-ternary-decision.md)

## Evaluating a supplied class power

View $W_0=Hi+j$ in the Hurwitz order $\mathcal H$. The intersection of $\mathbb Q(W_0)$ with $\mathcal H$ is exactly $O=\mathbb Z[W_0]$ of discriminant $-4M$. Indeed the fourth coordinate of $a+bW_0$ is zero, so its four Hurwitz coordinates must be integers. The scalar and $j$ coordinates then force $a,b\in\mathbb Z$. Thus the action uses the proper invertible ideals of this order, including when $M$ has square factors.

For a primitive positive definite form of discriminant $-4M$, use binary composition and Gauss reduction to compute its $e$-th power $[a,b,c]$. Use the ideal convention $I=(a,b/2+W_0)$. Generate the left Hurwitz ideal $\mathcal H I=\mathcal H g$ by conjugating its two generators, taking a Hurwitz left gcd, and conjugating back. Return the point
$$
W=gW_0\bar g/a.
$$
Here $\operatorname{nrd}(g)=a$. Locally $I$ is generated over $O$ by an element of the quadratic field, so extension to $\mathcal H$ squares the index $[O:I]=a$. Also $\mathcal H I$ is stable under right multiplication by $W_0$, making $W$ integral. Locally, $g$ differs on the left from a quadratic generator by a Hurwitz unit, so the rotated embedding is still optimal and $W$ primitive. Hurwitz generators differing by a left unit give properly rotated points. Form reduction bounds $a\le 2\sqrt{M/3}$. Binary powering and Euclidean quaternion gcd take polynomial time in the bit lengths of the form, $e$ and $M$. The class number and a factorization of $M$ are unnecessary. [Existing form arithmetic](../../src/nonsmooth/forms.py) · [Existing Hurwitz gcd](four_square_sampler.py) · [Mosunov, §§4.3–4.4](https://arxiv.org/html/1311.3379v2#S4.SS3)

The classical action on primitive sphere points modulo proper integral rotations is transitive. Its stabilizer here is $\{1,[\mathfrak t_2]\}$, where $\mathfrak t_2$ is the unique prime above 2. Indeed $M$ is $1\bmod 4$ or $2\bmod 8$, so $2$ ramifies in the field and the conductor of $O$ is odd. At odd primes an optimal embedding has the projective local order units as its integral torus stabilizer. At two the projective quaternion division algebra is compact and preserves the lattice. The local quadratic uniformizer, in addition to the units and rational scalars, therefore stabilizes the proper rotation orbit. Its class squares to $[(2)]=1$ and is nontrivial since no element $s+tW_0$ has norm $s^2+Mt^2=2$. Khayutin's Remark 8.15 describes precisely this extra quotient. Equivalently, the primitive three-square formula here gives $12h(-4M)$ points. Every proper rotation orbit has size 24 because a nonidentity cube rotation fixes only an axis of primitive squared length 1, 2 or 3. [Khayutin, Proposition 8.8 and Remark 8.15](https://ems.press/content/serial-article-files/32576?nt=1#page=53) · [Voight, Chapter 30, optimal embeddings and the primitive three-square formula](https://link.springer.com/chapter/10.1007/978-3-030-56694-4_30)

## Identifying all six classes

Since $-M\equiv-A\pmod N$, the two conjugate invertible primes $\mathfrak r_+,\mathfrak r_-$ above $r$ have norm $r$. For each one let $g$ generate its extended left Hurwitz ideal. Conjugation by $g$ is a rotation of reduced denominator exactly $r$. Its norm is $r$ and it cannot have a rational prime divisor as a Hurwitz element. The denominator formula for primitive integral Euler quaternions gives precisely $r$ after clearing any halves. [Euler denominator proof](ternary-frame-rotations.md)

Take $u_0=(x,1,0)$ and $v_0=(-N,0,0)$. The rotation sends $v_0$ to an integer vector since $r\mid N$, and it sends $W_0$ to an integer vector by the ideal construction. Its image of $u_0=W_0+Qv_0$ is therefore integral as well. The resulting frame is informative. If the image of $v_0$ were axial, the first column of the rotation would be integral. Subtracting $x$ times that column from the image of $u_0$ gives the integral second column, and their cross product gives the integral third column. This contradicts denominator $r$.

The two conjugate ideals yield distinct proper rotation orbits. Equality would force $[\mathfrak r_+]^2$ to be $1$ or $[\mathfrak t_2]$. The former would give a principal ideal $\mathfrak r_+^2$ of norm $r^2$. Every nonreal element of $O$ has norm at least $M>2r^2$, and the real elements $\pm r$ generate $\mathfrak r_+\mathfrak r_-$. The latter would make $\mathfrak t_2\mathfrak r_+^2$ principal of norm $2r^2$, below $M$ and not an integer square. Each orbit contains 24 points, exhausting the 48 informative frames. The base orbit accounts for the 24 coordinate frames. In the Picard group the six return classes are consequently
$$
1,\quad [\mathfrak t_2],\quad
[\mathfrak r_+],\quad[\mathfrak r_-],\quad
[\mathfrak t_2\mathfrak r_+],\quad[\mathfrak t_2\mathfrak r_-].
$$

For the four informative ideals the indices $a=r,2r$ are squarefree, so each has a primitive form $[a,2B,(B^2+M)/a]$. Translating $B$ makes $-a<2B\le a$. As $M>(2r)^2$, the last coefficient exceeds $a$, and the form is already reduced. A gcd of its leading coefficient with $N$ gives $r$. For a uniform Picard class the exact informative probability is $4/h(-4M)$, irrespective of how many coordinates or orientations are tested. This statement concerns the uniform distribution. A power of a chosen class lies in its own cyclic subgroup.

There is a deterministic restriction on words in small primes. If an integral invertible ideal $J$ of norm coprime to $N$ gives an informative orbit, multiplication by the conjugate of the appropriate norm-$r$ ideal and, if needed, by $\mathfrak t_2$ gives a principal integral ideal. Its norm is $r\operatorname{Norm}J$ or $2r\operatorname{Norm}J$. Its $r$-adic valuation is odd, so its generator is nonreal and the norm is at least $M$. Thus $\operatorname{Norm}J\ge M/(2r)$. Binary powering can pass this bound. The preceding point-transporter analysis already obtained the same scale for the reduced denominator of an individual rational rotation and computed $[r,2B,(B^2+M)/r]$ **from a supplied target**. The class action computes an orbit in the forward direction and identifies the entire returning subset of the Picard group. Mosunov’s factorization connection uses ambiguous classes and factors of the represented norm $M$. Here $M$ is coprime to $N$, and the two informative prime classes have distinct inverse orbits even after quotienting by $[\mathfrak t_2]$. [Earlier transporter and prime-neighbor calculation](packed-point-rotations.md) · [Earlier Gaussian conductor and perpendicular classes](ternary-shapes.md) · [Mosunov, §5](https://arxiv.org/html/1311.3379v2#S5) · [Ellenberg, Michel and Venkatesh, sphere class action](https://arxiv.org/abs/1001.0897)

## A large exponent does not force a return

Take $N=77$, $x=1$, $Q=2401$ and $M=34179874885$. Exponentiation by
$$
e=|\mathrm{GL}_2(\mathbb Z/7^4)|=7^{12}(7^2-1)(7^2-7)=27904034997216
$$
kills every invertible matrix modulo $Q$. A generator of $\mathcal H I^e$ need not be the ordinary $e$-th power of a generator of $\mathcal H I$. With the public split prime 19, coprime to $NQ$, the executable class action gives
$$
[19,14,(49+M)/19]^e=[5981,5498,5716006],\qquad
W=(-53506,174585,-28932).
$$
The sorted least absolute residues of $W\bmod Q$ are $(120,684,688)$. Every allowed $u$ has squared length 2 and gives $(0,1,1)$. Even starting at the informative norm-11 class $[11,6,(9+M)/11]$, the same exponent gives $[166934,36086,206701]$ and residues $(384,917,925)$. This disproves the proposed universal-return exponent. Selecting a useful power from $N,x$ alone remains open.

Run from the repository root:

```bash
python experiments/polylog_factoring_2026_11_27/check_packed_sphere_class_action.py
```

The checks compare the six classes to all 72 independently enumerated frame images on seven inputs. The fixture $N=143,x=0$ includes $25\mid M$. A full enumeration at $N=15,x=0$ gives 2,760 primitive classes and 1,380 proper rotation orbits with exactly two classes apiece. Exactly six classes meet the packed union.
