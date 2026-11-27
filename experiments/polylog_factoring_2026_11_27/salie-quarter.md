# An approximate Salié sum suffices for factoring

For an odd input with exactly two distinct prime divisors, a constant-precision approximation to one normalized Salié sum on uniform unit scales suffices for randomized polynomial-time factoring. Each summand is computable from public input by modular inversion and the Jacobi symbol. The sum has substantial cancellation. Single-term importance sampling needs order $N$ independent terms for constant mean-square precision on balanced semiprimes. [Fourier recovery](quarter-fourier.md)

## Exact identity

Write $e_N(t)=\exp(2\pi i t/N)$ and $\chi(t)=(t/N)$. For a unit $a$ set
$$S_N(a,b)=\sum_{m\in(\mathbb Z/N)^\times}\chi(m)e_N(am+bm^{-1}).$$
Let $G_N=\sum_{x\bmod N}e_N(x^2)$ be the quadratic Gauss sum, of value $\sqrt N$ for $N\equiv1\pmod4$ and $i\sqrt N$ otherwise. This definition also handles repeated factors. Then
$$\frac{S_N(a,b)}{G_N}=\chi(a)\sum_{v^2\equiv4ab\bmod N}e_N(v).$$

To derive the identity, take the additive Fourier transform in $b$ of the root sum. Its unnormalized coefficient at $z$ is $\sum_v e_N(v-zv^2/(4a))$. If $z$ is a nonunit, translation by $N/r$ for any prime $r\mid(z,N)$ preserves the quadratic term and multiplies the linear phase by a nontrivial $r$-th root of unity. The coefficient thus vanishes. For unit $z$, complete the square and evaluate the quadratic Gauss sum. This gives $G_N\chi(-z/(4a))e_N(a/z)$. Fourier inversion followed by $m=z^{-1}$ gives the identity, using $G_N^2=\chi(-1)N$.

Take $a=1$, $b=c^2/4\bmod N$, and let $u$ be a mixed root. For a unit $c$, the four roots are $\pm c$ and $\pm cu$. Hence
$$\frac{S_N(1,c^2/4)}{G_N}
=2\cos(2\pi c/N)+2\cos(2\pi cu/N).$$

## Recover a factor from approximate values

Suppose a public algorithm approximates the real normalized sum with root-mean-square error at most $1/2$ on uniform units, including its internal randomness. Subtract the known cosine, divide by two, and clip to $[-1,1]$. The resulting value has mean-square distance at most $1/16$ from $g(c)=\cos(2\pi cu/N)$. Convert it to a randomized sign of that expected value if needed by the sampler, and let $F(c)$ be the expected sign. By Jensen the same error bound holds for $F$.

For an odd modulus with two distinct prime divisors, the unit average of $g^2$ is at least $1/2$, by the Ramanujan sum at the unit frequency $2u$. Cauchy–Schwarz gives $\mathbb E_{c\in U} F(c)g(c)\ge1/4$. Extend $F$ by zero on nonunits. A uniform gcd probe can supply a factor directly. Otherwise the correlation on all residues is at least $1/8$, since $\varphi(N)/N\ge1/2$. This is the real part of the Fourier coefficient at $u$. The sampled heavy-frequency search recovers $u$ or $-u$ among its list. Test each output frequency minus 1 by gcd. This avoids a harmonic search over an unknown multiple of $u$.

## Variance of sampling individual terms

The normalized complex summands each have absolute value $N^{-1/2}$. For any positive sampling probabilities $q_m$, the unbiased one-term importance estimator has second moment $\sum_m |a_m|^2/q_m\ge(\sum_m |a_m|)^2=\varphi(N)^2/N$ by Cauchy–Schwarz. On the two-prime-divisor family the target magnitude is at most 4, so the variance is $\Omega(N)$ for large inputs.

Even estimating just the real part has average variance $\Omega(N)$ on balanced semiprimes. Its absolute term magnitudes are $N^{-1/2}|\cos\theta_{m,c}|$, where the phase includes division by $G_N$. Since $|\cos\theta|\ge\cos^2\theta$, it suffices to average the latter. For each fixed $m$, the nonconstant part in $\cos^2\theta_{m,c}$ is a unit quadratic phase in $c$. The complete Gauss bound and deletion of nonunits give an average $O(N^{-1/2})$ over unit $c$ on balanced semiprimes. Thus the average sum of absolute real magnitudes is $\Omega(\sqrt N)$. Jensen and the same importance-estimator bound give the stated average variance. Averaging $L$ independent one-term samples needs $L=\Omega(N)$ for constant mean-square error within this model. A method that aggregates oscillatory terms algebraically is outside this variance calculation.

## Checks

The diagnostic compared the public normalization with the factor-supplied cosine identity on 138 unit scales across five moduli, including two repeated-factor inputs. It also computed the minimum possible variance of a real one-term estimator for each tested scale. Run `python check_salie_quarter.py`.

## Fixed sparse linear combinations on square arguments

Suppose N=pq has comparable distinct odd prime factors. Even after
the public cosine is supplied exactly, fitting arbitrary fixed
complex weights to polynomially many individual summands leaves
mean-square error close to 2 on uniform unit shifts. The weights and
chosen summand indices may depend on N but are fixed across shifts.
This is a statement about this linear approximation class and the
uniform unit distribution, not the particular short decoder query list.

Indeed the shift dependence of each summand is a chirp
f_a(h)=e_N(a h²) with unit coefficient a=(4m)^(-1). Absorb all
constant phases and normalization in the fitted weights. The target
after public subtraction is g(h)=e_N(uh)+e_N(-uh). Use the inner
product given by averaging on unit h. Put B_r=(sqrt(r)+1)/(r-1)
for r=p,q. Each target-chirp inner product has magnitude at most
2 B_p B_q: CRT factors each quadratic-linear sum into two prime
sums, and deleting zero costs at most one from the complete Gauss
bound sqrt(r). For two distinct chirps the same argument applies
wherever their coefficients differ modulo a prime. The local
factor is one where they agree. For sufficiently large p,q every
off-diagonal Gram entry therefore has magnitude <= mu=max(B_p,B_q).

For k distinct chirps with (k-1)mu<1 the least Gram eigenvalue
is at least 1-(k-1)mu by Gershgorin. The squared norm of the
orthogonal projection of g onto their span is consequently at most

    4 k B_p² B_q² / (1-(k-1)mu).

The target squared norm is 2+2/phi(N) by the unit Ramanujan sum
at 2u. Thus for k=poly(log N) the best mean-square error tends
to 2. Taking real parts of fitted sums only doubles the dictionary
size by adding conjugate chirps; duplicate coefficients can be
merged. This rules out fixed sparse reweighting as a uniform
unit-scale approximation at absolute error 0.2, even with optimal
weights. Shift-dependent selection and nonlinear aggregation remain
outside the argument.

## Polynomial combinations of a fixed set of summands

Take k fixed quadratic phases f_a(h)=e_N(a h²), allow their complex
conjugates, and fit a polynomial of total degree d with arbitrary
coefficients independent of h. Each monomial is again one quadratic
phase, since multiplication adds its coefficients modulo N. Padding
products with ones shows that at most K=(2k+1)^d distinct phases
occur, including the constant. This covers polynomials in real and
imaginary parts through their expressions using conjugates.

The sparse projection estimate above applies to all these phases.
One extension needs checking: a combined coefficient a can be a
nonunit or zero. Locally where a vanishes, its correlation with
either linear target is -1/(r-1) since the target frequency is a
unit. This satisfies the same B_r bound. The off-diagonal estimate
already allows a difference divisible by one prime. Thus if the
actual number of distinct phases is M and (M-1)mu<1, the projection
bound holds with k replaced by M.

In particular, if k is poly(log N) and d=o(log N/log log N), then
M<=K=N^o(1). For balanced growing primes the root-mean-square error
of the best such fit to the hidden cosine pair tends to sqrt(2).
The public cosine pair may be added separately without changing
this residual bound. The shift distribution here is uniform on
units. This degree argument does not bound general arithmetic
circuit size: repeated squaring produces large formal degree with
few operations, and dense phase support can have a short expression.

## What a square argument determines

Write A(h)=2 cos(2 pi h/N), g(h)=A(uh), T=A+g. Each unit
square fiber consists of ±h and ±uh, so its conditional mean of
g is T/2. For uniform units the minimum mean-square error among
*all* functions of h² estimating g is exactly N/phi(N). Indeed
the conditional residual is (g-A)/2. Both A and g have second
moment 2+2/phi(N). Their cross moment is
-2/(p-1)-2/(q-1): the four additive frequencies ±1±u are
divisible by exactly one prime and each normalized unit sum is
-1/(r-1) for the other prime r. Substitution gives N/phi(N).
Enumeration over square fibers confirmed this for all 45 pairs of
distinct odd primes through 31.

Consequently learning the hidden summand solely as a function of
the square has irreducible error even with unlimited degree.
The full normalized Salié sum T is constant on every fiber and
has no such ambiguity. To apply the sparse bound to T directly,
sum the four target-chirp correlation bounds (at frequencies
±1 and ±u). The projection numerator becomes 16 M B_p² B_q²
with the same Gram denominator. The full target second moment is
4+4/phi(N)-4/(p-1)-4/(q-1). Thus a fit with M=N^o(1) phases has
optimal RMS error tending to 2 for balanced growing primes. The
decoder uses h separately when subtracting A from an estimate of
T. These calculations distinguish that procedure from fitting g
on square arguments alone.

## A public gcd test for exceptional Gram entries

For an explicit list of M different quadratic coefficients modulo
N=pq, check gcd(a-b,N) on each pair. Any nonunit difference gives
a proper factor immediately. Otherwise every off-diagonal Gram
entry has modulus at most B_p B_q rather than max(B_p,B_q), since
both prime Gauss bounds apply at once. With beta=B_p B_q the full
target projection bound therefore improves to

    16 M beta² / (1-(M-1)beta), whenever (M-1)beta<1.

In particular for balanced primes tending to infinity,
beta=(1+o(1))/sqrt(N). For any fixed epsilon>0 and
M<=(1-epsilon)sqrt(N), a list passing the gcd test captures only
O_epsilon(1/sqrt(N)) of the target's squared norm. Its optimal
root-mean-square error still tends to 2 on uniform unit shifts.
An explicit list in that size range giving error <=0.2 on all
units must therefore contain a pair revealing a factor for all
sufficiently large balanced inputs. Checking all pairs is
polynomial in the bit length only for polynomial-size lists
(or under an independently established faster bound); this result
does not construct such lists. Repeated coefficients modulo N
must first be merged so that gcd=N cannot mask a proper collision.

## Sign-only prediction on uniform units

For N=pq with distinct odd primes write h=qi+pj mod N and choose
u=1 mod p, -1 mod q. Put C(h)=cos(2 pi uh/N). The full normalized
sum T(h) factors as 4 cos(2 pi i/p) cos(2 pi j/q). On uniform
units i,j are independently uniform nonzero residues. The sign
never vanishes. Multiplication by u exchanges the public and
hidden cosines while preserving T. Consequently

    E[sign(T(h)) C(h)] = E|T(h)|/4 = m_p m_q,

where m_r=(csc(pi/(2r))-1)/(r-1). To see the finite sum, reduce
2i modulo r to the centered representative in [-(r-1)/2,(r-1)/2].
Then |cos(2 pi i/r)| is cos(pi times that representative/r).
Summing this geometric progression gives csc(pi/(2r)); remove
i=0. Hence m_3=1/2, and m_r>=1/2 for r>=5 by csc(x)>=1/x.
As both primes grow the product tends to 4/pi².

Any sign predictor with error probability rho on uniform units
(including internal randomness) therefore has hidden cosine
correlation at least m_p m_q - 2 rho. Extend the predictor by
zero on nonunits. The real part of its normalized Fourier
coefficient at u is at least phi(N)/N times this bound. At 90%
accuracy this is uniformly positive for every such semiprime.
For sufficiently large balanced inputs even 80% accuracy suffices,
since 4/pi²>0.4. The sampled heavy-frequency algorithm described
above can list the hidden frequency; a proper gcd verifies the
split. This trades numerical approximation for a single binary
classification target. No efficient predictor is supplied.

## Beating the public sign baseline

For every fixed epsilon>0, a predictor of the full sign with unit
accuracy >=3/4+epsilon on sufficiently large balanced semiprimes
suffices for randomized polynomial-time splitting on those inputs.
The public baseline sign(cos(2 pi h/N)) has accuracy tending to
3/4. Exhaustive unit-scale enumeration gave 0.745098 for 101*103
and 0.749506 for 1009*1013.

Here is a Fourier proof. Put s(x)=sign(cos(2 pi x)) on the circle.
The target in CRT coordinates x=i/p,y=j/q is s(x)s(y), and
h/N=x+y mod one. In the continuous limit the mean target at fixed
x+y=z is the convolution L(z)=1-4||z||, by intersecting two half
circles. Its sign is s(z), and its absolute integral is 1/2.
The CRT grid averages converge by Riemann integrability, including
for the public baseline's product with the target. Discarding
nonunits costs o(1). This proves the baseline accuracy.

Extend the predictor by zero on nonunits, using its expected output
if randomized. Its correlation on the full grid is >=1/2+2epsilon-o(1).
Take the symmetric partial Fourier sum P_K of s through orders ±K.
Its L2 error tends to zero by Parseval, so P_K(x)P_K(y)
approximates the target in L1 on the torus by Cauchy-Schwarz.
For fixed K this error passes to the growing CRT grids by Riemann
integrability. Choose K depending only on epsilon to make the error
less than epsilon/4 for sufficiently large inputs.

In this finite expansion each mode (a,b) is the additive frequency
k mod N with k=a mod p and k=b mod q. The diagonal a=b contributes
a function of z alone with coefficients s_hat(a)^2. This diagonal
polynomial converges uniformly to L since the squares of the real
Fourier coefficients are absolutely summable. For large enough K
its absolute grid mean is <=1/2+epsilon/4. Its correlation with
any bounded predictor obeys the same bound. Hence the off-diagonal
modes carry correlation >=epsilon for large inputs. By the triangle
inequality some such k has predictor Fourier magnitude at least
epsilon divided by the finite sum of absolute product coefficients
(an epsilon-dependent constant). Sampled heavy-frequency recovery
lists k in polynomial time at fixed epsilon. Test gcd(k-a,N) for
every listed k and |a|<=K. When p,q>2K and a differs from b,
the matching gcd equals p. This reduction requires no knowledge
of either prime to perform the search. It assumes access to the
predictor on random unit shifts and supplies no predictor itself.

## Quantitative advantage and unbalanced semiprimes

The previous reduction has bit cost polynomial in log N, 1/epsilon
and the cost of querying the predictor. No balance promise on
distinct odd p,q is necessary. Here are bounds making the dependence
on epsilon uniform. Assume 0<epsilon<=1/4 and a known lower bound
3/4+epsilon on accuracy at uniformly random units.

Set K=ceil(C/epsilon²) for a sufficiently large absolute constant C.
Trial-divide through ceil(C' K²/epsilon), with a sufficiently large
absolute C'. If this fails then both primes exceed that bound. For
r=p,q let alpha_r(a) be the normalized DFT coefficient of the public
quarter sign on the r-grid. The geometric sum of its interval gives
|alpha_r(a)|<=1/|a| for nonzero centered a. Truncating at K costs
O(1/sqrt(K)) in normalized L2 by Parseval. On the product grid
the product of the two truncated series approximates the target
in L1 with error O(1/sqrt(K)) by telescoping the product and
Cauchy-Schwarz. Each truncated factor has L2 norm at most one.

For |a|<=K compare alpha_r(a) with the continuous coefficient
alpha(a). On grid cells avoiding the two sign jumps, the exponential
varies by at most 2 pi |a|/r; the jump cells cost at most 4/r.
Hence the coefficient error is O((1+|a|)/r). As all coefficients
have magnitude at most one, summing errors in diagonal products
costs O(K²/min(p,q)) uniformly in z. The remaining continuous
diagonal tail costs O(1/K) by the squared coefficient decay.
The absolute grid mean of the triangle L differs from 1/2 by
O(1/N) by Lipschitz continuity. Extending the predictor by zero
on nonunits costs O(1/min(p,q)) in correlation. The parameter
choices make all errors together smaller than epsilon. Therefore
the off-diagonal correlation is at least epsilon.

The sum of the absolute truncated product coefficients is
O(log²(2K)), by the harmonic sum bound. Some off-diagonal
frequency thus has magnitude Omega(epsilon/log²(2K)). The sampled
interval search in quarter-fourier.md lists it in time polynomial
in log N and 1/epsilon. Trying the 2K+1 possible first coordinates
by gcd then splits N. In particular any known inverse-polynomial
advantage in bit length above 75% suffices. This is an oracle
reduction for each promised input, not an estimator or a practical
sample-budget claim.

## Public intervals with Jacobi-dependent decisions

On balanced growing semiprimes, partition the public residue h
into B contiguous intervals, and in each choose an arbitrary sign
rule on chi(h). Even the best such table for that N, with B
polynomial in log N, cannot attain 3/4 plus an inverse-polynomial
advantage in the bit length on uniform units.

Indeed write its expected output on units as A(h)+chi(h)D(h),
where A,D are constant on each interval and bounded by one.
Extend it by zero on nonunits. The unmasked A on all residues
has normalized Fourier coefficient O(B/|k|_N) at k nonzero
by the geometric sum. Removing nonunits changes this by at most
1/p+1/q. The chi D term is already zero on nonunits. At every
frequency its coefficient is O(B log(N)/sqrt(N)): the normalized
Fourier transform of the primitive Jacobi character has sup norm
1/sqrt(N), and the normalized Fourier l1 norm of an interval is
O(log N). Multiplication by the interval convolves these spectra.

If accuracy reached 3/4+epsilon, the preceding quantitative proof
would yield a coefficient of magnitude Omega(epsilon/log²(2K))
with residues a,b of size <=K=O(epsilon^-2), a different from b.
Its centered representative must have magnitude >=min(p,q)-K:
otherwise the two congruences force it to equal both a and b.
For inverse-polylog epsilon and polylog B the preceding upper
bound there is polylog(N)/sqrt(N), contradicting the required
coefficient for sufficiently large balanced N. The result allows
arbitrary interval endpoints but does not cover general functions
of the bits of h. Exhaustive majority fitting for N=1009*1013
with 32 equal-width bins and both Jacobi signs gave 0.74950593
accuracy, the public baseline to the displayed precision.

## Interval tests on a shift and its modular inverse

Take a predictor made of R rectangular cells I x J in the public
pair (h, inverse(h) mod N), with bounded weights of magnitude <=1
and zero output off units. For a partition this includes an arbitrary
sign table. For balanced N=pq its normalized Fourier coefficient
at a nonzero k has magnitude at most

    O(R (1/|k|_N + log²(N)/sqrt(min(p,q)))).

For one cell expand both interval indicators in additive Fourier
series. At k the complete unit sum in modes a,b, divided by N,
has phase (a-k)h+b/h. When (a,b)=(k,0) this contributes at most
the magnitude of the k-th interval coefficient, O(1/|k|_N).
Otherwise at least one prime has a nonzero coefficient pair.
At that prime the normalized complete Kloosterman sum is
O(1/sqrt(r)) by Weil if both coefficients are nonzero, and of
magnitude 1/r if precisely one is nonzero. The other CRT factor
has magnitude at most one. Each interval has normalized Fourier
l1 norm O(log N), giving the stated error after summing all modes.

For R polynomial in log N the bound at the off-diagonal frequencies
from the quantitative reduction is negligible in bit length:
|k|_N >= min(p,q)-K there. Consequently this model cannot reach
75% plus an inverse-polynomial advantage on uniform unit shifts for
large balanced semiprimes. This is specific to a short rectangle
description, not general arithmetic with modular inverses.

For p=1009, q=1013, exhaustive majority fitting with 32 equal
bins on each axis got 764584 of 1020096 unit shifts correct
(74.95216%). These are best in-sample labels using known-factor
ground truth; they do not constitute a public training method.

## Exact finite-modulus score of the public cosine rule

With known distinct odd primes p,q put P=(p-1)/2 and Q=(q-1)/2.
The public rule is correct on exactly

    2 P Q + 2 floor_sum(P, 2p, 2q, q)

of the 4PQ unit shifts. The repository floor_sum(n,m,a,b) sums
floor((ai+b)/m) for 0<=i<n by Euclidean recursion. Thus one can
calibrate planted benchmarks with polynomial bit cost.

To derive the count, subtract the nearest half-integer from each
CRT angle i/p,j/q (in turns). The residuals are signed v/(2p),
w/(2q), with each 1<=v<=P and 1<=w<=Q taking both signs once:
doubling and centering permutes the nonzero residues. The cosine
sign picked up from the two half-integers is exactly the product
of the local cosine signs. Agreement therefore asks whether the
sum of residuals has magnitude below 1/4. Opposite signs always
agree. Equal signs agree iff v/p+w/q<1/2. At fixed v there are
floor(q(p-2v)/(2p)) such w, with no endpoint tie. Reverse the
v order to obtain the displayed sum. Reflection of the rectangle
by (v,w)->(P+1-v,Q+1-w) sends every point below the line to one
strictly above it. Hence the exact public score never exceeds 75%,
even at small or unbalanced moduli.

Verified against full unit enumeration for all 45 pairs of distinct
odd primes through 31 and for 1009*1013. On the latter benchmark
the public rule gets 764568 correct, exactly 16 fewer than the
best in-sample shift/inverse bin table.

## A square-only sign predictor needs any nonnegligible advantage

Suppose a predictor gets N and y=h² mod N (equivalently y/4),
where h is a uniform unit, and predicts the sign of the full
normalized sum with probability at least 1/2+epsilon. This sign
is well defined on y. Then N=pq for distinct odd primes can be
split in randomized time polynomial in log N, 1/epsilon and the
predictor cost. In particular the 75% benchmark for an algorithm
given h does not apply to this square-only model.

Define f(h) as the expected answer on h² and set it to zero on
nonunits. Since u²=1, f(uh)=f(h) exactly. Its Fourier coefficients
at k and uk therefore agree, including for randomized prediction.
The unit correlation with the target is at least 2epsilon.
Trial-divide as in the quantitative proof above, so both primes
are large compared with K=C/epsilon² and the loss on the full grid
and the truncated-product L1 error together cost less than epsilon/2.
The (0,0) term is negligible: each local sign mean on the full
odd prime grid is ±1/r, giving coefficient magnitude 1/N.
Consequently some nonconstant truncated frequency with residues
(a,b), |a|,|b|<=K, has predictor Fourier magnitude at least
Omega(epsilon/log²(2K)) by the same coefficient l1 bound.

If a differs from b the gcd search already succeeds. If a=b
is nonzero, multiplying the frequency by u gives residues (a,-a)
with exactly the same predictor magnitude. These differ because
both primes exceed 2K. The bounded-error heavy-frequency list
therefore contains a splitting frequency in either case; no
knowledge of u is needed to find it. Run the list search with
fresh predictor calls on squared sample inputs and test all small
offsets by gcd. A square-only algorithm with inverse-polylog
advantage is sufficient. A predictor given a particular square
root can break the required symmetry, as the public baseline does.
