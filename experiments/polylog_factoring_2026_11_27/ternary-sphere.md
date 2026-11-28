# Uniform primitive points on a sphere of known odd radius

Given any odd $N\ge1$, one can sample exactly uniformly from all signed,
ordered primitive triples $x^2+y^2+z^2=N^2$ in expected polynomial time in
$\log(2N)$. The algorithm uses $N$ alone and includes points with a zero
coordinate. The square radius is what lets a four-square sampler supply
the inputs. For a general squarefree squared norm $n$, the three-square
count involves a quadratic class number. Pollack and Schorn give a
randomized representation algorithm for admissible $n$ with a polynomial
runtime bound under ERH. [Earlier class-number investigation](../../notes/geometry_results.md) ·
[Pollack and Schorn](https://doi.org/10.1090/mcom/3349) ·
[Primitive sampler](ternary_sphere.py)

## Four Gaussian-unit lifts per point

Draw an exactly uniform primitive four-square representation of $N$ using
the existing public norm sampler. Interpret its four coordinates as two
Gaussian integers $u,v$. The classical parametrization of Pythagorean
quadruples gives the Hopf point

$$
x=|u|^2-|v|^2,\qquad y+iz=2u\bar v,\qquad N=|u|^2+|v|^2.
$$

Reject when $\gcd(x,y,z)>1$. Choose one of the three cyclic permutations
of an accepted triple uniformly. [Classical parametrization and counts](https://www.maths.tcd.ie/EMIS/journals/JIS/VOL18/Huerlimann/huerli6.pdf) ·
[Four-square sampler and proof](four-square-incidence.md)

For odd $N$, exactly one of $|u|^2,|v|^2$ is odd, so $x$ is odd and $y,z$
are even. A common Gaussian divisor of $u,v$ contributes its norm to
every coordinate. Conversely, if an odd rational prime $r$ divides all
three coordinates, it also divides $N$ and both Gaussian norms. For an
inert $r$, this forces $r\mid u,v$. For $r=\pi\bar\pi$ split, the two
zero norms and $u\bar v=0\bmod r$ force a common Gaussian prime divisor.
Indeed each nonzero residue of norm zero in
$\mathbf Z[i]/r\simeq\mathbf F_r\times\mathbf F_r$ has just one nonzero
component, and $u\bar v=0$ requires the same component for both. This
argument includes a zero residue by using any prime dividing the other
Gaussian norm. Thus primitivity is equivalent to $u,v$ being coprime in
$\mathbf Z[i]$.

Every primitive sphere point has exactly one odd coordinate since
$N^2\equiv1\bmod4$. First put it in the $x$ position. For $x\ne N$,
write $(y+iz)/(N-x)=u/v$ in lowest Gaussian terms. Stereographic
projection gives the same normalized point $ (x,y,z)/N $ as the Hopf
formula. If both Gaussian norms were odd, both would be 1 modulo 4.
The raw Hopf coordinates would then have gcd exactly 2 by the odd-prime
argument and the fact that $u\bar v$ has odd norm. After division by 2
their first coordinate would be even. The chosen point therefore comes
from opposite norm parities, with raw gcd one and raw radius exactly $N$.
The pole $x=N$ is primitive only for $N=1$ and comes from $(u,v)=(1,0)$.
The reduced fraction determines its pair up to a common unit
$1,-1,i,-i$, also at the poles. This proves exactly four lifts per
first-coordinate-odd point. The three cyclic shifts cover disjoint possible
positions for the odd coordinate, establishing the uniform law.

## Counting accepted inputs

Put $\chi(r)=(-1/r)$ at odd primes. The number of coprime Gaussian pairs
of squared norm sum $N$ is

$$
A(N)=8N\prod_{r\mid N}(1-\chi(r)/r).
$$

To see this, sort all pairs by their Gaussian gcd ideal, choosing any
generator for each ideal. Dividing by that generator gives a unique
coprime pair up to the fixed choice. Thus the ordinary count $8\sigma(N)$
is the Dirichlet convolution of $A$ with Gaussian ideal norm counts.
At $r$, the power-series factor for $\sigma$ is
$1/((1-t)(1-rt))$ and the ideal count has factor
$1/((1-t)(1-\chi(r)t))$, by splitting or inertia. Their quotient is
$(1-\chi(r)t)/(1-rt)$, yielding the formula at every odd prime power.
There are consequently $3A(N)/4$ primitive triples, in agreement with
the classical primitive three-square formula. [Cooper and Hirschhorn, equation (1.19)](https://doi.org/10.1007/s11139-006-0240-6)

The input sampler is uniform on the $8N\prod_{r\mid N}(1+1/r)$ rationally
primitive quaternion representations. Each accepted Gaussian pair is in
that population. Hence the acceptance probability is exactly
$$
\prod_{\substack{r\mid N\\r\equiv1\ (4)}}\frac{r-1}{r+1}
\ \ge\ (\varphi(N)/N)^2 .
$$
The expected number of quaternion draws is $O((1+\log\log(3N))^2)$.
Computing the required split of the quaternion algebra once per radius
allows it to be reused across rejection attempts.

## Perpendicular forms and factor information

A primitive point $P=(x,y,z)$ gives the full perpendicular lattice
$P^\perp\cap\mathbf Z^3$ directly. For $g=\gcd(x,y)>0$, choose
$sx+ty=g$. The vectors $e=(-y/g,x/g,0)$ and $f=(-sz,-tz,g)$ satisfy
$e\times f=P$, so span the full saturated lattice. The axial case
$g=0$ gives the ordinary coordinate plane. Its binary norm form
$[e\cdot e,2e\cdot f,f\cdot f]$ has discriminant $-4N^2$.
It is primitive. At an odd prime a nonsingular ternary quadratic space
has no totally isotropic two-dimensional subspace. Modulo 2, a plane on
which every vector has even norm would equal $(1,1,1)^\perp$, which would
force three odd coordinates in $P$.

Reduction of the form is public arithmetic. It need not be the principal
class. For example, $(1,2,2)$ at radius 3 gives the class of $[2,2,5]$,
which represents 2 and does not represent 1. A general split from the
corresponding class would need another operation, such as finding a
factor-revealing ambiguous class. On a squarefree input, inspecting the
content of a rejected Hopf point itself can reveal a prime $r\equiv1\bmod4$
only on a fraction $2/(r+1)$ of raw draws per such prime. Indeed either
Gaussian prime above $r$ must divide both inputs, each event has exactly
$8\sigma(N/r)$ lifts, and the two events are disjoint. A prime
$r\equiv3\bmod4$ never divides that content for squarefree $N$.
[Prior ambiguous-form analysis](../../notes/codex_r4_sos_classgroup.md)

At the 50 odd radii through 99, exhaustive checks compared 24,824
quaternions to independent three-coordinate enumeration. They verified
all 5,498 four-element unit fibers, all 16,494 primitive ternary points
and their binary form discriminants, plus the identity between Hopf content
and Gaussian gcd norm for every input. There were also 64 public sampler
calls, including repeated prime powers and radius $1009\cdot1013$.
[Reproduce the checks](check_ternary_sphere.py)
