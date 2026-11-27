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
