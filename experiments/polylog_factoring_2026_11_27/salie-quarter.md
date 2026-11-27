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

## Grouping summands with the same phase

On the distinct odd semiprime promise, put b=c²/4 for a unit c and
let h(m)=m+b/m modulo N. The public partner b/m has exactly the
same complex summand as m. More generally h(t)=h(m) is equivalent
to (t−m)(tm−b)=0 modulo N. At each prime the choices are m and
b/m, and they have identical local quadratic character because b is
a square. Thus every summand in a phase fiber is identical. If a
third distinct element t is found outside the two public choices,
it agrees with m at exactly one prime, so gcd(t−m,N) factors N.

Negation conjugates the normalized summand because χ(−1)G_N equals
the conjugate of G_N. Grouping the negative phase as well therefore
preserves the sum of absolute real parts. For real orbit weights w,
the optimal single-orbit unbiased importance estimator has second
moment (∑|w|)², exactly the same as the optimum for individual real
terms. Thus the independent-sampling variance bound below persists
even if these complete groups are supplied. Exhaustive checks on
N=35,77,143 and all square unit targets verified constant Jacobi sign
on 1,692 fibers and proper gcds on all 5,856 ordered pairs outside
the two public choices.

The follow-up question of computing the signed weight at a phase is
already settled as an oracle equivalence on odd squarefree semiprimes.
[Pointwise weights and quadratic residuosity](salie-phase-weights.md)
gives both reductions. [Grouped sampling variance](salie-grouped-variance.md)
covers independent sampling even with pointwise access. A sufficient
open target is the [signed phase interval](salie-phase-bins.md).

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

## Fourier tests on the square argument itself

For uniform unit h let Y=h² mod N and S be the full sign. Every
additive Fourier mode on Y satisfies

    |E[S e_N(aY)]| = O(log(2p) log(2q) / sqrt(N))

uniformly in a, including zero and nonunits. Indeed write
h=qi+pj. The mixed quadratic term vanishes modulo N, so the
expectation factors into averages of s(i/p)e_p(aq i²) and
s(j/q)e_q(ap j²) on nonzero local residues. At a prime r with
nonzero quadratic coefficient, expand the square wave in its
additive Fourier series. Its coefficient l1 norm is O(log(2r)),
and each completed quadratic sum has magnitude sqrt(r); removing
zero changes the unnormalized sum by one. With zero quadratic
coefficient the unnormalized sign sum on nonzero inputs has
absolute value at most two. Both cases meet the stated bound.

For a bounded predictor B(Y), Fourier expansion therefore bounds
its absolute sign correlation by its normalized Fourier l1 norm
times the displayed bound. A predictor constant on L intervals has
l1 norm O(L log(2N)), even with arbitrary fitted interval boundaries.
This also applies after relabeling Y by a public invertible affine
map, which preserves that norm. Polylogarithmically many intervals
give negligible advantage in bit length for growing semiprimes.
The same test applies to any proposed square-only predictor whose
Fourier l1 norm can be bounded, regardless of its implementation.

Exhaustive majority fits using 32 equal intervals gave 5628/10200
for 101*103 and 514124/1020096 for 1009*1013. These retrospective
scores use known factors to choose the labels independently for each
modulus. Summing over h weights every unit square equally.

## Shifted Jacobi symbols on the square input

For Y=h² and k fixed public shifts d_j, allow an arbitrary decision
function on the k Jacobi symbols (Y+d_j | N). On balanced semiprimes
its absolute sign correlation on unit h is bounded by

    O(k/min(p,q) + 2^(k/2) (k+1)^2 log(2p) log(2q) / sqrt(N)).

To see this, expand the Boolean decision function on nonzero answers
in Walsh monomials. Its multilinear extension stays in [-1,1] and
its coefficient l1 norm is <=2^(k/2). Each zero event has probability
<=2/(p-1)+2/(q-1) by the quadratic root count. On a monomial the
sign correlation splits by CRT into two local averages of the quarter
sign times the quadratic character of a product of <=k polynomials
of the form (other_prime*i)^2+d_j. Remove even root multiplicities
at a cost of O(k) omitted terms per unnormalized sum. If the remaining
polynomial is nonconstant, complete the quarter sign in additive
Fourier modes. The Weil bounds for the squarefree polynomial character
sum, including a linear additive twist, give O((k+1)sqrt(r)) per
completed sum at prime r. The sign has Fourier l1 norm O(log(2r)).
If the polynomial was a constant times a square, the character is
constant away from its roots and the signed sum is O(k+1). This
also covers the empty monomial. Multiplying local bounds and summing
Walsh coefficients proves the estimate without any condition on
collisions among the shifts.

For fixed eta>0 and k<=(1-eta)log2 N the advantage is negligible
in the bit length on balanced inputs. Using d=1,2,3,4 and allowing
separate labels for zero answers, exact best in-sample accuracy
was 52.2353% at 101*103 and 50.3811% at 1009*1013. To enumerate,
tabulate signed counts by the four local Legendre symbols at each
prime, then multiply counts while multiplying the symbols
coordinatewise. The computed labels use the factors.

## Public labels for square-only learning and held-out scoring

Sample a uniform unit h, give the predictor only (N,h² mod N),
and label the example with the readily computed sign of cos(2 pi h/N)
by an integer quarter comparison. For any fixed square-only predictor
its accuracy on these public labels is EXACTLY its accuracy on the
hidden quarter labels sign(cos(2 pi uh/N)): substitute h -> uh,
which permutes units and leaves the square fixed. This concerns
population accuracy; individual labels may differ. Hence even an
inverse-polynomial held-out population advantage on this public task
suffices by the quarter Fourier reduction. Neither training examples
nor validation examples require factoring. Fix the trained model
before testing and use independent draws to avoid selection bias.
The label task is distinct from predicting the full Salié sign.

There is also a graded label: the conditional mean of cos(2 pi h/N)
given the square is T(h)/4 by averaging the four roots. Its
correlation with any bounded square-only predictor equals correlation
with the hidden cosine, giving a directly estimable Fourier statistic.

check_square_only_score.py tests a simple 32-bin majority model with
100000 training draws and 100000 separate scoring draws per modulus.
The fixed-seed diagnostics returned 52.862% for 10403 and 49.920%
for 1022117. For independent random validation draws the two-sided
99% Hoeffding radius at this sample size is 0.515 percentage points
per model. Fixed seeds are for reproducibility. These results do
not establish sustained advantage as bit length grows.

## Maximum public-label accuracy from the square

Within each unit square fiber the two root pairs contribute two
copies each of signs s(h/N) and s(uh/N). Thus the label distribution
is either unanimous or exactly tied. The full Salié sign attains
the best possible public-label accuracy: on a unanimous fiber the
sum of the two cosine pairs has their common sign, and on a tied
fiber either prediction is optimal. Its public-label score is
precisely the public cosine rule's full-sign score (the same two
signs are compared). Hence the exact floor sum above also computes
the Bayes ceiling for any square-only public-label predictor,
including one with unlimited computation. This ceiling is <=75%.

More precisely, writing U for the probability of a unanimous
fiber and a for a model's conditional accuracy on those fibers,
its public-label score is 1/2+U(a-1/2). Its decisions on tied
fibers are invisible to this score. At N=1009*1013 exactly
127260 of 255024 unit squares are unanimous, giving ceiling
74.95059%. Exhaustive fiber enumeration confirmed the floor-sum
ceiling for the 45 products of distinct odd primes through 31 and
also for 101*103 and 1009*1013.

## Recognizing unanimous fibers is the original quarter-parity task

For a unit square y=h² modulo a distinct odd semiprime let Q(y)
count its roots in the first quarter. Each of the two sign pairs
contributes one iff its folded root has positive cosine. Consequently
Q is 0 or 2 on unanimous fibers and 1 on tied fibers. With +1
denoting unanimity, the recognition target is exactly

    (-1)^Q(y) = s(h/N) s(uh/N).

Thus a proposed recognizer on y (or one given h as well) can be
turned into a hidden-quarter predictor by multiplying its output
by the computable s(h/N). This preserves correctness on every
individual query and hence preserves average advantage exactly.
The reverse operation also works when h is given. Exact recognition
is the full original parity evaluator on unit scales in sign form;
the previous heavy-frequency reduction already covers a noisy version.

A method may abstain on difficult scales. If it returns a recognition
answer on a fraction alpha of uniform unit scales and has conditional
accuracy 1/2+gamma there, filling abstentions with independent fair
signs gives advantage alpha*gamma. No independence between coverage
and the true class is needed. A public efficiently recognizable subset
with both coverage and conditional advantage inverse-polynomial in
bit length would therefore suffice. The public single-root labels
from the training experiment do not certify recognition accuracy:
on tied fibers they cannot reveal which of the two classes applies
from the observed label alone.

## Coverage of searching integer square lifts

For a uniform unit square y represented in [0,N), test whether
y+kN is an integer square at L distinct nonnegative levels chosen
from N alone. For consecutive levels starting at zero, the covered
squares are exactly the distinct unit residues x² mod N for
1<=x<=isqrt(LN-1). For any fixed choice of levels the count of
covered residues is at most sqrt(LN)+L: each interval
[kN,(k+1)N) has at most sqrt(N)(sqrt(k+1)-sqrt(k))+1 integer
square roots. The widths decrease with k, so their sum is at most
the sum for the first L levels. Repeated residues can only reduce
coverage. Divide by phi(N)/4 to get the hit probability. Since
phi(N)/N>=8/15 for distinct odd primes, polylog L gives negligible
coverage in the bit length, even if the chosen levels are large.

On a hit, return a root x selected solely from y and independent
coins. The sampled h conditioned on y is uniform among the four
roots, so gcd(h-x,N) splits N with probability exactly one half:
the two roots outside ±x each agree with x at precisely one prime.
Thus root witnesses would bypass label prediction if they covered
an inverse-polynomial fraction. This level bound is specific to
levels chosen independently of y.

Exact covered-residue counts for L=1,16,64 respectively were
100,388,744 out of 2550 at N=101*103, and 1009,4024,8016 out
of 255024 at N=1009*1013, by enumerating and deduplicating the
unit square residues up to the stated root cutoff.

### Tight order for quotient lists

For distinct odd primes and 1≤L≤N, the best target-independent list has hit probability of order sqrt(L/N), uniformly in the modulus. The preceding upper bound gives this order since L≤sqrt(LN). For the matching lower bound use levels 0 through L−1, and put R=isqrt(LN−1)<N. Of the integers from 1 to R at least R(1−1/p−1/q)≥7R/15 are units. Each unit square can occur at most four times in that interval since it has four roots modulo N. Dividing by the number of unit squares gives the claimed lower bound. Thus choosing scattered or very large levels cannot improve the asymptotic order.

If the target chooses from K menus of at most L levels, their union has at most KL levels. Even perfect menu selection therefore gives O(sqrt(KL/N)) average success for KL≤N. To reach inverse-polylogarithmic success this family needs KL at least N/polylog(N). This counts distinct quotient choices and applies even when selection uses the full target. It leaves open how to compute a useful quotient from a large implicit family.

## Integer-affine choices of the lift level

Let the input representative range over 0<=y<N. An integer-affine
choice k(y)=ay+b gives x²=(1+aN)y+bN. Put A=|1+aN| (nonzero)
and B=bN. The nonnegative roots x lie in one interval of length
at most sqrt(AN), since the nonnegative radicands span at most AN.
Also x²=B mod A. We count roots of this last congruence as a
union of progressions even when B has large square factors.

For a prime power r^e exactly dividing A, if r^e divides B then
x is zero modulo r^ceil(e/2). Otherwise existence forces the
valuation of B to be even, say 2t<e. Divide x by r^t: the
result is a unit square root modulo r^(e-2t). There are at most
two unit roots for odd r and at most four for r=2. This gives
at most that many classes for x modulo r^(e-t). Every local
period is at least r^(e/2). CRT thus gives at most 2^(omega(A)+1)
classes with a common period >=sqrt(A), also valid for A=1.
Each class supplies at most sqrt(N)+1 values of x in the interval.
Each x gives at most one input y. This bound is independent of b.

For |a|<=N^C with fixed C, omega(A)=O(log N/log log N) since
omega(A)!<=A. Even a menu of poly(log N) such rules, with any
input-dependent selection among them, covers at most N^(1/2+o(1))
residues. Dividing by phi(N)/4 gives hit probability at most
N^(-1/2+o(1)) for uniform unit squares on distinct odd semiprimes.
This includes piecewise-affine predictors with polynomially many
pieces. It does not bound arbitrary polynomial-bit slopes of
superpolynomial magnitude, or a nonlinear rule with no short
affine menu. Prime-power decompositions are used only in the proof.

### Extension of the same coverage test to polynomial selectors

The fixed-level and affine tests above address the same witness
mechanism. For k=P(y) with integer coefficients, degree and coefficient
bit lengths polynomial in log N, the existing additive large-sieve
coverage inequality (see ../../notes/codex_r15_lean_faithfulness.md)
also gives O(sqrt(N) log N) hits per polynomial. Here is the extra
hypothesis check needed to apply it. F(X)=X+NP(X) reduces to X
at every prime dividing N, so cannot be a constant times a square
over Q: write any such representation as an integer constant times
the square of a primitive integer polynomial and reduce it modulo
that prime. Thus F=c H² G with integer factors and a nonconstant
squarefree primitive G. Outside primes dividing
2c lc(H) lc(G) disc(G), G stays nonconstant and squarefree.
The exceptional integer has polynomial bit length by the polynomial
factor height bound and the determinant formula for the discriminant.
At good primes r between sqrt(N)/2 and sqrt(N), the quadratic
character bound for G gives at most r/2+O(deg(F) sqrt(r)+deg(F))
residues on which F is square, counting zeros of H separately.
This is <=3r/4 for sufficiently large N. Only poly(log N) primes are
exceptional, leaving order sqrt(N)/log N good primes. For a hit set
of size M each such prime contributes at least M²/3 nonzero Fourier
energy by Parseval and Cauchy-Schwarz on occupied residue classes.
The additive large sieve bounds their sum by 2NM, proving the claim.
A polylogarithmic menu still covers a negligible fraction of unit
squares. This reuses the earlier sieve and is a coverage bound only.

## Distinct-prime dependence of two sign reductions

For odd N with s>=2 distinct primes, let d_i be their full prime
powers. On unit h the normalized Salié value is the sum of e_N(x)
over the roots of h². Write h=sum_i (N/d_i) j_i by CRT. Independent
local sign choices give a product of 2 cos(2 pi j_i/d_i). Thus
the exact sign, extended by zero off units, has a Fourier coefficient
at every mixed root of one equal to the product of local quantities
(d_i)^(-1) sum_{j unit mod d_i} |cos(2 pi j/d_i)|. Imaginary parts
cancel by reflection. Each factor is >=1/6: the unit density is
>=2/3 and the unit average of cosine squared is >=1/4 by the
Ramanujan sum at frequency two. This gives a threshold 6^(-s)
for the sampled Fourier search. Independent of any factorization,
one can try decreasing thresholds and verify gcds. The resulting
oracle reduction is polynomial in log N and 6^s, hence polynomial
in bit length for s=O(log log N). A predictor with error rho on
unit shifts changes each normalized coefficient by at most 2 rho.
This bound for the full Salié sign degrades with s.

There is a different conclusion for the sampled-root public label.
Given only h², any predictor for sign(cos(2 pi h/N)) with advantage
epsilon yields a split in time polynomial in log N and 1/epsilon,
for any s>=2. Indeed extend the expected predictor by zero off
units. Its correlation with the public square wave on all residues
is >=2 epsilon phi(N)/N >=2 epsilon/(s+1), since the i-th distinct
prime is at least i+1. The one-dimensional Fourier truncation in
quarter-fourier.md yields a heavy frequency with nonzero small
centered index a (trial division handles small N). If gcd(a,N)>1
it already splits N. Otherwise invariance under h -> uh guarantees
the same coefficient at au for every mixed root u. The sampled
list and gcd(k-a,N) recover a split without knowing u or s, using
s<=log2 N to set thresholds. This promise concerns predicting the
sampled root's label, not the full Salié sign on all its conjugates.
It supplies no predictor with that advantage.

## Information limit for the sampled-root label at many prime factors

Put s=omega(N). Even an unrestricted square-only predictor has
population advantage at most O((1+s)^2 2^(-s/4)) for predicting
the sampled root's public quarter sign. This qualifies the preceding
sufficient condition: along families with s/log log N tending to
infinity, the required inverse-polynomial advantage is impossible.

For the proof extend the bounded predictor f(h²) by zero on nonunits.
Its normalized Fourier transform is constant on each orbit under
multiplication by all square roots of unity. For nonzero centered
frequency a with |a|<=K=floor(2^(s/2)), at most log2 K components
can have their full prime power dividing a. Every other component
has two distinct sign choices on a, since N is odd. The orbit thus
contains at least 2^(s/2) frequencies. Parseval bounds the absolute
value of each corresponding coefficient by 2^(-s/4). The public
quarter wave has constant coefficient of magnitude 1/N; its
nonzero coefficients have magnitude <=1/|a| by the interval sum.
The low-frequency contribution to correlation is therefore
O((1+s)2^(-s/4)), and the tail costs O(K^(-1/2)) by Cauchy-Schwarz
and Parseval. Convert full-grid correlation to unit correlation
using N/phi(N)<=s+1, and divide by two for advantage.

Exact majority scores obtained by enumerating every square fiber
at N=15,105,1155,15015,255255 were respectively 4/8, 32/48,
280/480, 3224/5760, 49684/92160. These examples are diagnostics,
not a monotonicity claim. At N=15 the label is exactly balanced
on every unit square. This obstruction concerns predicting a random
preimage label; the aggregate window parity is deterministic on the
square, and the bound does not transfer to that target.

## Factored evaluation for all odd moduli and scales

For the normalized complete sum T(c)=S_N(1,c²/4)/G_N, suppose the
prime powers r=p^e exactly dividing N are supplied. Put v_r=(N/r)⁻¹
modulo r. The root-sum identity and CRT give T(c) as the product
over these r of the following real numbers.

* If p does not divide c, the factor is 2 cos(2π(c v_r mod r)/r).
* If p divides c and e=1, the factor is 1.
* If p divides c and e≥2, the factor is 0.

Indeed the additive phase factors over CRT; locally the roots are
±c in the first case and only 0 in the second. In the last case all
roots are divisible by p. Translation by r/p preserves the root
set, and multiplies its sum by a nontrivial p-th root of unity. The
sum is therefore zero. This covers nonunit scales and repeated
prime factors without enumerating roots.

With k distinct prime divisors, absolute error 2^(−t) costs polynomial
bit time in log N and t. For example enclose each local factor to
error 2^(−t)/(k·3^k), so the telescoping product error is at most
2^(−t) even allowing approximate factors of magnitude up to 3,
using exact dyadic multiplication afterwards. A direct FFT of the
original summands agreed on every scale of every odd N from 3 through
251 (15,875 values, maximum absolute discrepancy below 5·10^(−15)).
The supplied factors determine the CRT phases; obtaining them from the
unfactored index remains the computational issue for this formula.

## Factoring from the logarithm of the magnitude

A public evaluator for log|T(c)| with pointwise absolute error 1/32
at unit scales suffices for randomized polynomial-time factoring of
arbitrary odd N. No sign or phase of T is needed. With the factors
known, the preceding cosine product supplies this precision in
polynomial time. The log-error promise differs from additive error
in the original sum when its magnitude is small.

For the proof extend the unit-scale log magnitude to all c by the
function L(c)=∑_{r} log(2|cos(2π c v_r/r)|), with r and v_r as
above. The minimum angular distance from a cosine zero on an odd
grid of size r is at least π/(2r). Thus each factor lies between
1/r and 2 and |L|≤log N. In particular even an additive approximation
to |T| within 1/(64N) suffices on units.

For each exact prime power r, consider k=2v_r(N/r) mod N, with
gcd(k,N)=N/r. Every other local summand has zero Fourier coefficient
there by CRT. To bound the real part of the surviving coefficient,
take x uniform mod r, z=cos(4πx/r), and f(z)=log(2+2z)/2.
Then E z=0 and Var z=1/2. Since f'(z)≥1/4 on the sampled range,
Cov(f(z),z)≥Var(z)/4=1/8 by the independent-copy covariance identity.
This is the desired Fourier coefficient's real part. Consequently
each component gives a frequency of constant unnormalized signal
strength (normalized Fourier transform of the unscaled function L),
regardless of the number of distinct primes.

Simulate bounded-function Fourier queries by checking gcd(c,N)
first. A proper gcd ends the search successfully. At c=0 return
zero, changing any coefficient by at most log(N)/N. On units clip
the evaluator to [−log N,log N]. Equivalently analyze a hypothetical
function filled with exact L on proper nonunits: the simulation
agrees until such a query immediately yields a divisor. The surviving
coefficient after scaling by log N has size Ω(1/log N) for sufficiently
large N. Use the [sampled heavy-frequency algorithm](quarter-fourier.md)
and check gcd(k,N) for all returned frequencies. Small N and prime
powers can be handled directly; recurse on proper splits and amplify
success. FFTs of the proof extension for every odd N≤2000 with two
or more distinct prime factors checked all predicted frequencies on
676 moduli. The smallest observed real coefficient exceeded 0.2310.
These checks use the factors and test the reduction.
