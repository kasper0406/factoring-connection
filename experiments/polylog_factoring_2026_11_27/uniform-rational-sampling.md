# Low degree rational sampling for the square residue decoder

A low-degree rational map applied to independent uniform residues can increase the square-residue decoder's success probability by at most a polynomial factor in the degree, except when coefficient gcds or a deterministic constant already give the split. The statement allows coefficients to depend on the unfactored modulus.

## Fiber bound

Let $N=pq$ for odd primes $p<q\le2p$ and $B=\lfloor\sqrt N\rfloor$. The decoder first tries $\gcd(z,N)$. For a unit $z$ it checks whether $z^2\bmod N=t^2$ as ordinary integers and tries $\gcd(z-t,N)$. The successful unit set has exactly $2(B-1)$ elements. For every unit $1\le t\le B$ there are two mixed CRT signs, and the only excluded $t$ is $p$. The nonzero nonunits add $p+q-2$ successes. Write $S$ for their union.

Take $z=F(U)/G(U)$ for $U$ independent uniform coordinates in $\mathbb Z/N\mathbb Z$, with total degrees at most $D$. When the denominator is a nonunit, allow a gcd attempt on the denominator. Assume the reductions of $G$ are nonzero polynomials and $F/G$ is formally nonconstant over each prime field.

For a fixed desired residue $a\bmod r$, where $r=p$ or $q$, the equation $F(U)-aG(U)=0$ is a nonzero polynomial equation of total degree at most $D$. The Schwartz–Zippel bound gives probability at most $D/r$. The two reductions of uniform $U$ are independent by CRT. Hence each $a\bmod N$ can be output with probability at most $D^2/N$, counting only unit denominators. The probability that the denominator is a nonunit is at most $D/p+D/q$, again by the polynomial zero bound. In total the probability of returning a divisor in one trial is at most
$$\min\left(1,\frac{D^2(2(B-1)+p+q-2)}{N}+\frac{D}{p}+\frac{D}{q}\right).$$

The bound is $O((D^2+D)/\sqrt N)$ on balanced inputs. It also applies conditionally to a map chosen from previous trials if the new inputs are again uniform and independent after conditioning on that history.

## Constant reductions and zero denominators

For explicitly supplied coefficients, first take gcds of the denominator coefficients with $N$. A proper nonunit coefficient splits $N$. Otherwise either the denominator is identically zero modulo $N$, which makes the rational expression undefined, or some coefficient $g_j$ is a unit. Write $f_j$ for the numerator's coefficient at the same monomial, including zero. Form the coefficient list of $H=g_jF-f_jG$. Over either prime field, $F/G$ is formally constant exactly when every coefficient of $H$ vanishes. If this holds at just one prime, a nonzero coefficient of $H$ modulo the other prime supplies a proper gcd. If it holds over both primes, the map on unit denominators equals the single public constant $f_j/g_j$. Apply the decoder once to this constant. Denominator gcds remain subject to the $D/p+D/q$ bound when both denominator reductions are nonzero.

For polynomially many explicitly described maps of polynomial degree, the union bound keeps the random success probability exponentially small in the bit length after the coefficient and constant cases are resolved. This leaves room for nonuniform samples, degree exponential in the bit length with a succinct evaluation circuit, or a direct construction of a successful constant. [Roadmap decoder and uniform baseline](https://github.com/kasper0406/factoring-connection/blob/9514cbfa013cd9ae3f030044db1dc3fb6bb4a9ed/notes/research_roadmap_nov26.md)
