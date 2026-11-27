# One quadratic count with a public threshold

The gap oracle in the count comparison reduction can be stated as one ordinary lattice-point count at a public threshold. Its input is a bounded rational polytope in ten coordinates and one integer quadratic equation whose quadratic part has rank two. The reduction needs one count comparison for each sampled upper digit. A polynomial-time evaluator for this decision problem remains open. [Gap comparison and isolation](count-comparison.md)

## Pack the signed contrast by complementation

Let $M=N^2$, $H=(N-1)/2$, and $c_1=c_0+Nj$. Let $b_0$ be the isolation offset, $w$ the width, and $\alpha,\beta$ the affine sign parameters. The input and phase representatives are
$$
 x=s_0-b_0+Ny,\qquad t=t_0+Nz,
 \qquad 0\le s_0<w,\quad 0\le y,t_0<N,\quad 0\le z\le H .
$$
Choosing $y$ over a full period gives each allowed $x\bmod M$ once, including when $s_0-b_0<0$. A selector $\eta\in\{0,1\}$ records the sign class by a linear constraint with a unique integer quotient $k$,
$$
 (H+1)\eta\le \alpha(s_0-b_0)+\beta-Nk
 \le H+(N-1-H)\eta .
$$
Sign class zero has positive sign. An independent selector $s\in\{0,1\}$ chooses the reference or query digit. Let $C_{s,\eta}$ be the count satisfying $xt-x^2\equiv c_s^2\pmod M$ in the corresponding region. The signed contrast before the public stationary correction is
$$
 D_{\rm raw}=C_{1,0}-C_{1,1}-C_{0,0}+C_{0,1}.
$$

Choose a polarity $\sigma\in\{1,-1\}$. In one polytope, encode $\lambda=s\mathbin{\mathrm{xor}}\eta$ with four linear inequalities,
$$
 \lambda\ge s-\eta,\quad \lambda\ge\eta-s,\quad
 \lambda\le s+\eta,\quad \lambda\le2-s-\eta.
$$
For $\sigma=1$ set the complement flag $\nu=1-\lambda$, and for $\sigma=-1$ set $\nu=\lambda$. Here $\nu$ is an affine expression. Add an integer residual $r$ with $\nu\le r\le(M-1)\nu$, and impose
$$
 xt-x^2-c_0^2-s(c_1^2-c_0^2)-r=Mq .
$$
For each choice of the other coordinates, there is exactly one possible residue $r$ and hence one quotient $q$ before applying the residual bounds. When $\nu=0$, the constraint counts the original congruence. When $\nu=1$, it counts its complement, including pairs with nonunit $x$. All those pairs belong to the known ambient population. The quotient bounds $|k|\le N+1$ and $|q|\le8M$ suffice with parameters chosen in their residue ranges.

For either polarity there is one complemented selector case for each sign class. Across those two cases the entire ambient population has size
$$
 B=wN^2(H+1).
$$
Thus the single count is $K_\sigma=B+\sigma D_{\rm raw}$. The ten integer coordinates are $s_0,y,t_0,z,k,s,\eta,\lambda,r,q$. The quadratic part $xt-x^2$ has rank two because $x,t$ are independent affine coordinates. Every coefficient and bound has $O(\log N)$ bits.

## Use a one-sided public threshold

Let $P_j-P_0$ be the public signed stationary correction, computed from the blocks of $x\equiv\pm c_0\pmod N$. Then $D(j)=D_{\rm raw}-(P_j-P_0)$. The count comparison can answer zero for
$$
 K_\sigma\le B+\sigma(P_j-P_0)+N/4
$$
and one for
$$
 K_\sigma\ge B+\sigma(P_j-P_0)+3N/4,
$$
with either answer in the gap. Rounding the lower threshold down and the upper threshold up gives integer thresholds.

Choose $\sigma$ uniformly once per isolation trial. In a singleton trial the signed stationary change, whenever nonzero, has one fixed sign because the reference membership bit and the affine sign of the isolated block are fixed. The polarity agrees with that sign with probability one half. The small-remainder analysis then makes the one-sided bit agree with the hidden membership change except on at most $1/32$ of the upper digits on a good aligned trial. Its heavy Fourier coefficient and gcd decoder are unchanged. All conclusions hold with repeated prime factors. [Fourier recovery](quarter-fourier.md)

## Checks

The diagnostic verified 30 packed counts and threshold comparisons against independently counted congruences for both polarities on three small moduli, including a repeated-prime input. It counts nonunit pairs in the complement. Run `python check_single_quadratic_count.py`.

## Comparison with the finite-subspace identity in grh-lean

At commit 72fead1 of lsernau/grh-lean, New.q_binomial_submodule
sums over subspaces U of an m-dimensional space over a field of
order q. Each dimension-i subspace has weight (-X)^i q^(i(i-1)/2).
The total equals the product of (1-q^j X) for 0<=j<m. Grouping by
dimension gives the Gaussian binomial coefficients. In the Hecke
transfer application an integrality constraint restricts U to the
subspaces of a single residue kernel. The same product then uses
the kernel dimension. At X=1 this sum is the indicator of zero
nullity, including the empty-product case.

The analogous step for the present comparison would be a public
encoding that collapses its weighted fibers while preserving the
threshold gap. Even over a known odd prime field and before applying
the interval filters, x(t-x)=c^2 for c nonzero is not an affine
subspace: it contains (c,2c) and (-c,-2c) but omits their midpoint.
Thus rank-nullity on the existing coordinates does not yield that
encoding. This leaves auxiliary-variable representations open.
The product identity computes a full weighted sum; specializing
the Gaussian coefficient polynomial at q=1 recovers ordinary choose
numbers and supplies no new fast individual-coefficient procedure.

[Subspace identity](https://github.com/lsernau/grh-lean/blob/72fead102f663d5c7a407964a13dc0cf9d977207/grh/src/Lib/General/Combinatorics.lean#L89) · [Kernel restriction](https://github.com/lsernau/grh-lean/blob/72fead102f663d5c7a407964a13dc0cf9d977207/grh/src/Lib/Weil/Unramified/HeckeTransfer.lean#L166)

## A two-dimensional kernel and the missing aggregate

For each candidate pair put A=[[x,c],[c,t-x]] over Z/(N^2), with
c a unit. Its determinant is exactly the tested residual. Using c as
a unit pivot reduces A by invertible operations to diag(1,det A)
up to a unit on the second diagonal. Hence the number of kernel
vectors equals gcd(det A,N^2), and the candidate passes exactly when
this number equals N^2. This converts the pointwise predicate into a
kernel-size test over the ring. It still leaves the sum over candidates.

Even the complete collection of residue-field ranks of this matrix
loses the lift digit: for c=x=1 compare t=2 and t=2+N. Both give
rank one at every prime dividing N. Their kernel sizes modulo N^2
are respectively N^2 and N. Both phase digits are allowed for odd
N>=3. The same information loss occurs for any matrix polynomial in
the undivided coordinates reduced modulo N. An encoding using the
upper digits separately need not suffer that loss.

A sufficient alternative target is a polynomial-size matrix built
from public input whose ranks at two prime divisors of N differ.
Scan the remaining matrix modulo N for a nonzero entry. A nonunit
gives a proper factor by gcd. A unit permits row and column elimination
and decreases both field ranks by one. Reaching a zero remainder
would imply equal initial ranks. Thus constructing a separating matrix
already suffices for a split without computing any subspace sum.
No such construction from the aggregated comparison is supplied here.

## Projecting the linear lift onto the upper phase digit

Fix a unit lower digit a and a unit scale c. Let b be the least
residue of a+c²/a modulo N, and h=(a(b-a)-c²)/N. Substituting
x=a+Ny, t=b+Nz and scale c+Nj into the quadratic congruence
gives exactly (b-2a)y+az=2cj-h modulo N. The full row has rank
one at every prime divisor since a is a unit. The projection onto
z is more informative: with g=gcd(b-2a,N)=gcd(a²-c²,N), its
image is one residue class z=z_j modulo g, where
z_j=a^(-1)(2cj-h) modulo g. Each element of the image has g
preimages. For H=(N-1)/2 and 0<=z_j<g, the exact clipped count
is g floor((H+g-z_j)/g). Thus the N upper lifts of one lower
digit cost only polynomially many bit operations in log N.

If g=1 the clipped count is independent of j and cancels against
the reference. If 1<g<N the computation already supplies a divisor.
If g=N then a²=c² modulo N. Outside the two publicly corrected
digits a=±c modulo N, gcd(a-c,N) is proper: odd prime powers
dividing N must divide one of a-c and a+c exclusively.
Consequently *any* explicitly located lower digit with nonzero
corrected contrast yields a verified divisor. This explains why
collapsing the inner lifts still leaves the outer search. It gives
no bound on the work of locating such a digit in a wide window.

## Window product as a factor-or-zero certificate

Remove the two public digits ±c from the cyclic interval I, and
compute P=product_{a remaining in I}(a-c)(a+c) modulo N. If
gcd(P,N)=1, every remaining unit block has g=1 and cancels, and
nonunits have no solutions. Thus the corrected contrast is exactly
zero for every upper shift and affine sign on this window.

If the gcd is proper it gives a divisor directly. If it equals N,
bisect the product, descending into a child with gcd greater than
one until a proper gcd appears or a single digit remains. At a
leaf with gcd N, neither linear factor vanishes modulo N because
the two public digits were removed. Their product vanishes, hence
gcd(a-c,N) is proper. This branch can give a factor even on a
window whose signed contrast happens to cancel.

Splitting at the cyclic wrap and removed points gives a constant
number of ordinary intervals. Each of their two shifted factorial
products can be computed with the classical square-root block
product and monic multipoint evaluation. Bisection on one path
keeps the total at softly O(sqrt(|I|)) modular operations by the
geometric series for the square roots of the subinterval lengths.
The certificate bypasses upper-digit sampling on a successful
window but still requires a fast product or a different detector
for a wide interval. It supplies a factor-or-zero interface, not
the value of a nonzero comparison.

## An exact factorial specialization of the window detector

For odd composite N take 1<=w<N/2. If gcd(w,N)>1 we already have
a proper divisor. Otherwise set 2c=w modulo N and choose the cyclic
window a=c+i for 1<=i<=w. Neither public digit is present, because
a-c=i and a+c=w+i modulo N and both representatives lie in [1,2w].
Thus the window product is exactly (2w)! modulo N. Its nonunit bit
tests whether the least prime divisor p is at most 2w. In particular
the first positive width is (p+1)/2. Binary search up to (N-1)/2
therefore recovers p in O(log N) queries, or splits N early via
gcd(w,N). Recurse and primality-test to factor completely.

Conversely, given the distinct prime divisors r of N, the general
window-product nonunit bit is easy. Split the cyclic window at wrap
and the excluded public digits. On each remaining interval [L,U],
test whether it contains an integer congruent to c or -c modulo r
by comparing floor quotients at U and L-1. This gives a deterministic
polynomial-time oracle equivalence for the Boolean detector on arbitrary
public windows. It does not assert the same encoding inside a promised
distribution of isolated signed comparisons: the factor-or-zero
detector can succeed where the actual signed contrast vanishes.

## Reflection cancels the half interval except at a boundary digit

For a unit lower digit a let C_j(a) be its unsigned clipped lift
count. Keep the notation b,g,z_j,H from the projection calculation.
Define e=H when b is nonzero and e=0 otherwise. Then

    C_j(a)+C_j(-a) = N + g [z_j = e modulo g].

Indeed negating x negates the phase t=x+(c+Nj)^2/x modulo N^2.
If the low phase digit b is nonzero, its upper digit z becomes
N-1-z. The two half-interval indicators sum to one except at H,
where both are one. If b=0 the reflected digit is -z modulo N
and the extra unit is at zero instead. The projected image supplies
exactly g preimages of the boundary digit when it occurs.

After reference subtraction an equally weighted reflected pair
therefore contributes just the difference of two boundary counts.
At a stationary digit g=N, the shift slope 2c/a is a unit, so the
boundary is reached at exactly one j modulo N. After dividing by N,
the reflected contrast is a point indicator minus a constant;
every nonzero normalized Fourier coefficient has magnitude 1/N.
Thus enforcing equal reflection weights on an isolated stationary
pair loses the constant-size Fourier signal used by the existing
decoder. This calculation does not cover opposite weights or an
unpaired digit in the actual isolation window.

## Opposite reflection weights preserve the primitive harmonic

For a stationary unit digit a (a²=c² modulo odd N), put s=2c/a
modulo N. Its lower phase digit b=2a modulo N is nonzero. Its
upper phase z=z_0+sj runs through all residues as j varies. The
oppositely weighted pair (C_j(a)-C_j(-a))/N is +1 at z<H, zero
at z=H and -1 at z>H. Centering at H and taking the normalized
Fourier transform in j gives magnitude at frequency s equal to

    (2/N) sum_{k=1}^H sin(2 pi k/N) = cot(pi/(2N))/N.

Reference subtraction affects only frequency zero. The magnitude
therefore tends to 2/pi. In contrast, every nonstationary raw
block is periodic in j with period g=gcd(a²-c²,N)<N by the lift
formula. Its Fourier transform is supported on multiples of N/g,
so it contributes exactly zero at all frequencies coprime to N,
even with a j-independent window or sign. Since s is a unit, an
isolated stationary pair has no nonstationary contamination at s
in this linear statistic. Other stationary pairs can interfere.

This statement precedes the nonlinear threshold comparison.
Orthogonality in the full Fourier sum does not imply small
pointwise error or cheap access to an individual comparison bit.

## Sparse rational shifts that cancel a prime-period lift block

Suppose g=gcd(a²-c²,N)=p is a proper odd prime divisor. For a
fixed lower digit, the clipped count as a function of j is a
constant plus p times the indicator of an interval in z_j modulo p.
The interval length is (H+1) mod p=(p+1)/2, and z_j runs affinely
through all residues with unit slope. Its first nonconstant harmonic
is nonzero: summing a nontrivial p-th root over a proper nonempty
interval gives a nonzero geometric series. Subtracting the reference
only changes the constant.

Consider an explicit linear combination of t<p distinct shifts
modulo N with nonzero rational weights. If the resulting block
vanishes for every j, its shift polynomial F must vanish at a
primitive p-th root. Reduce exponents modulo p and collect weights
in each class. The resulting polynomial of degree at most p-1 is
a rational multiple of 1+X+...+X^(p-1) by cyclotomic irreducibility.
Since t<p, some residue class is empty, so all collected weights
are zero. In particular any chosen shift has a distinct companion
congruent to it modulo p. Taking gcds of its differences with all
other shifts exposes a proper divisor in at most t-1 tests: such a
difference is divisible by p and is nonzero modulo N.

This argument needs exact cancellation of even one actual block of
prime period, for all shifts of the input digit j. It allows signed
weights and does not assume cancellation of arbitrary periodic
functions. Approximation, cancellation only after summing blocks,
and nonlinear comparisons remain outside the statement.

## One-step differences suppress the aggregate semiprime remainder

Let N=pq for distinct odd primes. Take the difference at j+1 and j
before thresholding. At any nonstationary block with g=r in {p,q},
the upper class moves by 2c/a=±2 modulo r because a²=c² modulo r.
The clipped count is a constant plus r times a half-interval
indicator on these classes. Its difference has magnitude r on at
most four shifts per period r and vanishes elsewhere. The expected
absolute normalized difference is therefore at most 4/N for a
uniform shift j. Blocks with g=1 vanish identically.

There are at most 2N/r lower digits having g=r, even before
windowing. Consequently for *any* fixed lower window and weights
of absolute value at most one independent of j, the nonstationary
sum R, normalized by N, satisfies

    E_j |R(j+1)-R(j)| <= 8/p + 8/q.

This uses only the triangle inequality, without assuming independent
phases or signs. Markov's inequality bounds the fraction exceeding
1/4 in magnitude by 32/p+32/q. A stationary block has a unit slope
s=2c/a modulo N. Its fundamental Fourier coefficient is multiplied
in magnitude by 2|sin(pi s/N)| under the same difference. That
factor can be small, so the remainder estimate alone does not give
a uniform decoder. Nor does sparsity per block locate the active
digits in polynomial bit time.

## Dyadic lags and small stationary roots

Suppose N=pq with distinct odd primes p<q<2p. Take an isolated
nonpublic stationary digit with slope s=2c/a modulo N and let
m=min(s mod N, -s mod N). Then m²=4 modulo N and m is not 2.
For fixed integer B>=1, if m<=B sqrt(N), testing 4+kN for squares
for 1<=k<=B² recovers a nontrivial root and a factor via gcd(m-2,N)
at the matching k. All these tests use public input only.

Otherwise consider the first power of two h with hm>=N/4.
It satisfies N/4<=hm<=N/2, including when h=1. The linear
difference at lag h multiplies the magnitude of the stationary
fundamental by 2|sin(pi hs/N)|>=sqrt(2). Telescoping the one-step
bound gives mean absolute normalized nonstationary error at most
8h(1/p+1/q)<=4(p+q)/m<12/B, since p+q<3 sqrt(N).
Thus for any fixed desired error level one may choose a constant B,
run the small-square tests and try O(log N) public dyadic lags.
At least one lag has both guarantees if the square tests did not
already split N. For B=1024 the mean error is below 0.012.

These are guarantees for the linear difference and an isolated
stationary signal. The absolute value or threshold of a difference
has a different spectrum. Computing the signed difference from
unfactored input in polynomial bit time is still open.

## One-sided comparison at the selected lag

Under the comparable-semiprime and singleton conditions above,
write the normalized corrected lag contrast as U(j)+E(j), with
U(j)=±(I(z+hs)-I(z)), z=z_0+sj and I the upper half-window
indicator used in the count. At the good dyadic lag let d be
the least absolute residue of hs. Then N/4<=d<=floor(N/2).
The positive part {U=1} is a single cyclic interval of length d
in the z coordinate, for either stationary sign. Its normalized
Fourier magnitude at s is sin(pi d/N)/(N sin(pi/N)), hence at
least 1/(sqrt(2) pi).

Use a one-sided comparator returning 0 below or at 1/4 and 1
above or at 3/4, arbitrary in between. Whenever |E|<=1/4 it
returns exactly the indicator of U=1. By the preceding mean bound
with B=1024, disagreement occupies at most 48/1024 of all shifts.
Even adversarial answers there change a normalized coefficient by
at most that fraction. Thus its magnitude at s exceeds 0.178.
This gives a constant heavy frequency for the existing sampled
decoder after thresholding. The comparison is a single polarity
of the packed count: take the two scales c+Nj and c+N(j+h)
(reducing shifts modulo N) and subtract the public stationary
correction for their difference. The lower digits and signed
window are unchanged. Uniform polynomial bit cost for this
one-sided threshold remains the required evaluator assumption.

## Window-sensitive lag bound for arbitrary semiprime ratios

For N=pq, distinct odd primes p<q, each prime-period class g=r
lies in the two progressions a=±c modulo r. Any cyclic interval of
width w contains at most 2(w/r+1) such digits. Applying the
one-step estimate only to these digits, then telescoping, gives

    E_j |R(j+h)-R(j)| <= (8h/N)(w/p+w/q+2).

For a nonpublic stationary slope, its least absolute representative
m is congruent to ±2 modulo q and m!=2. Hence m>=q-2.
The good dyadic lag above satisfies h<=N/(2m), giving a mean
error bound of 4[(w/N)(p+q)+2]/(q-2), with no factor-ratio
hypothesis. For instance if q>=2048 and w/N<=1/8192 this is
below 0.005. Smaller q can be handled by trial division through
2048. For singleton isolation under these width conditions the
one-sided comparator therefore disagrees on at most 0.02 of the
shifts; its normalized coefficient at s exceeds 0.20. The bound is
uniform in the window location and the fixed weights of magnitude
at most one (with stationary weight ±1). This removes the need for
the small-square preprocessing on this semiprime argument. The
actual comparison still needs a polynomial-time evaluator.
