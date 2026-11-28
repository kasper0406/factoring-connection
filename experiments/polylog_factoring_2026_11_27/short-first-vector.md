# Fixing the first vector in the ternary Gram search

For a product $N=pq$ of distinct odd primes, the four public choices
$x=0,1,2,6$ cover a limiting $15/16$ of pairs sampled from a common dyadic
prime interval. On every covered input a noncoordinate ternary Gram
witness exists **with first vector prescribed exactly** as $(x,1,0)$.
Two other small first vectors raise this density to $63/64$.
The denominator of a rational parametrization then identifies the
splitting prime. Finding that denominator from $N$ remains the search
problem. [Ternary Gram and its 48 factor-revealing witnesses](mixed-ternary-embeddings.md)

## Four choices and their exact domain

Take $A=x^2+1$ and test the Jacobi symbol $(-A/N)$. A value $-1$ guarantees
that exactly one of $p,q$, denoted $r$, splits in the sense that $-A$ is a
nonzero square modulo $r$. If a test gives zero, compute $\gcd(A,N)$.
For $x=0,1,2,6$, the norms $A$ are $1,2,5,37$. In each case **every**
integer vector of squared norm $A$ is a signed permutation of $(x,1,0)$,
as can be checked by listing the integer squares up to 37. An ambient
signed permutation can therefore move any of the 48 witnesses to one
with the specified first vector, preserving the gcd of its partner.

On inputs coprime to 1480, the four Jacobi tests give a negative value
exactly when $p$ and $q$ have different signatures under the characters
of $-1,-2,-5,-37$. These four outcomes are independent on the reduced
residue classes modulo $1480=8\cdot5\cdot37$. In fact they determine the
two independent signs of $-1,2$ modulo 8 and the quadratic characters
modulo 5 and 37. The Chinese remainder theorem gives 36 classes per
signature. Thus coverage is certain on every pair with different
signatures. The prime number theorem in these fixed progressions gives
limiting probability $15/16$ for distinct uniform primes in $[R,2R)$.
This fraction measures when the test supplies a Gram with a short
prescribed first vector. It does not measure factor recovery.

## The remaining binary norm equation

The same geometry permits any coprime integers $a,b$ and $c,d$ with
$ad-bc=1$. Use the coordinate witness
$$
u=(a,b,0),\quad v_0=N(c,d,0),\quad
A=a^2+b^2,\quad B=N(ac+bd),\quad C=N^2(c^2+d^2).
$$
The determinant is $AC-B^2=N^2$. The usual family takes $(a,b)=(x,1)$
and $(c,d)=(-1,0)$. For every other integral partner of this fixed $u$,
write $v=(Nc+bt,Nd-at,z)$. Its remaining norm constraint is
$$
z^2=t(2N-At).
$$
The case $t=0$ is $v_0$. For $t>0$, put $z/t=k/\ell$ in lowest terms
with $\ell>0$. Then
$$
D=k^2+A\ell^2,\qquad t=\frac{2N\ell^2}{D},\qquad z=\frac{2Nk\ell}{D}.
$$
Because $\gcd(D,\ell)=1$, this partner is integral exactly when $D\mid 2N$.
For a negative Jacobi test, the nonsplitting prime cannot divide $D$,
since $\ell$ is a unit at any prime dividing $D$. Hence
$D\in\{1,2,r,2r\}$. The first two values, when attained, give
$\gcd(v)=N$. Either of the latter gives $\gcd(D,N)=r$ and $\gcd(v)=N/r$.
Indeed $N/r$ divides all coordinates while $r\nmid t$, and $\gcd(a,b)=1$.
The partner has no other common prime outside $N$ since its minors
have squared sum $N^2$. Any successful slope has
$\max(|k|,\ell)\ge\sqrt{r/(A+1)}$. Even with constant $A$, the height is
at least of order $N^{1/4}$ on balanced pairs.

There are 72 frames in the mixed instance for this general coordinate
choice as well. The gcd $g$ of the minors divides $N$ and equals
$\gcd(v)$ by $Av_j-Bu_j=\sum_i u_i(u_i v_j-u_j v_i)$ and $\gcd(A,N)=1$.
Dividing $v$ by $g$ gives the primitive Gram with $m=N/g$, determinant
$m^2$ and the same unit $A$. Its frame count is
$24\#\{s\bmod m:s^2\equiv -A\bmod m\}$. Only $m=1,r$ contribute,
with 24 and 48 frames respectively. For a single first-vector orbit
of size $h$, exactly $48/h$ factor-revealing partners remain at the
prescribed $u$ by signed-permutation symmetry. [Primitive frame count](mixed-ternary-embeddings.md)

At $N=11009$ the tests for $x=0,1,2$ all give $+1$, and $x=6$ gives
$-1$. The slope $k/\ell=8/1$ has $D=8^2+37=101$, yielding the split
$11009=101\cdot109$. This verifies the parametrization on a known input.
At $x=0$ the norm equation reduces to the familiar two-square search.
The other norms replace it by a different fixed positive binary form.
[Earlier circle analysis](../../notes/geometry_results.md)

## Two more first vectors and completeness

Use $(a,b)=(2,3)$ of norm 13 and $(3,7)$ of norm 58, with complements
$(c,d)=(-1,-1)$ and $(-1,-2)$ respectively. Each norm again has a single
signed-permutation orbit. The two new characters $-13,-58$ supply the
independent quadratic signs of 13 and 29. The six signatures occupy
all 64 possibilities equally among the 193,536 reduced residue classes
modulo $557960=8\cdot5\cdot13\cdot29\cdot37$. Thus this menu covers all
different-signature pairs, of limiting density $63/64$.

Up to signs and permutation, the complete list of norms of primitive
coordinate vectors with a single three-square orbit is $1,2,5,10,13,37,58$.
For a vector $(a,b,0)$ with $0<a<b$ the orbit has 24 elements. If an
odd square divides $A=a^2+b^2$, the residual integer is still $1$ or $2$
modulo 4 and has a three-square representation by Legendre's theorem.
Scaling it gives a separate imprimitive orbit. Otherwise $A$ is
squarefree, and Gauss's formula gives $r_3(A)=12h(-4A)$ for $A>2$.
The even fundamental discriminants of class number two are
$-20,-24,-40,-52,-88,-148,-232$. The candidates 6 and 22 contain an
inert prime to odd exponent and cannot themselves be sums of two squares.
The exceptional smaller orbits come from $A=1,2$. Restricting to
$A=x^2+1$ gives $x=0,1,2,3,6$, and the character at $x=3$ introduces
no independent sign. [Kenku's classification of the even discriminants](https://doi.org/10.1112/plms/s3-22.4.734) ·
[Gauss's three-square formula](../../notes/geometry_results.md)

The checker exhausts the 16 and 64 signatures, verifies all seven orbits
directly, and compares direct partner enumeration with the complete
binary parametrization on 223 mixed Grams, containing 1,089 partners at
the prescribed first vectors. [Reproduce the calculations](check_short_first_vector.py)
