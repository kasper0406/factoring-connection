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

## One-sided lag comparisons on general good isolation trials

The [general isolation argument](count-comparison.md) gives, after
its stated preprocessing, a constant-probability singleton trial
with normalized mean square E_j |R(j)|² <= 1/512. For every lag h
on this same trial, translation invariance and (x-y)²<=2x²+2y² give

    E_j |R(j+h)-R(j)|² <= 1/128.

Hence at most 1/8 of shifts have difference error larger than 1/4.
For the nonpublic stationary unit frequency s, choose the first
dyadic h with h min(s mod N,-s mod N)>=N/4. As before its clean
positive pulse occupies one interval in phase and its normalized
coefficient at s is at least 1/(sqrt(2) pi). The one-sided gap bit
therefore has coefficient magnitude at least 1/(sqrt(2) pi)-1/8
> 0.10. No semiprime assumption enters this deduction. The
mean-square guarantee holds on the same trial for every candidate
lag; no independence between lag selection and the window is needed.

This is a corollary of the existing general reduction for the
variant comparing two public shifts. It retains that reduction's
assumptions and preprocessing, including separate primality and
perfect-power handling. The remaining computational problem is
evaluating the packed one-sided comparison on unfactored inputs.

## All upper shifts act on one fixed congruence by a phase shear

Fix unit c and M=N². For shift j set k=j/c modulo N and
u=1+Nk modulo M. Multiplication by u on both x and t maps
x(t-x)=c² modulo M bijectively to x'(t'-x')=(c+Nj)² modulo M:
u c=c+Nj modulo M and u is a unit. This bijection fixes both
lower digits modulo N. In particular the public input window and
its affine weight are unchanged. For t=b+Nz with 0<=b,z<N the
new phase digits are exactly (b,z+kb modulo N).

Thus on the single base graph at c, each lag contrast can be
computed using the signed mask I(z+k' b)-I(z+k b), where k'=(j+h)/c
modulo N. The mask consists of two modular linear interval tests.
The two public stationary blocks can equivalently be removed at
the lower-input filter. Even the base graph can be normalized to
uv=1 modulo M by writing u=x/c, v=(t-x)/c; the input window is
then tested on cu modulo N and t=c(u+v) modulo M.

This gives a common modular inverse graph for all scales and
shifts at fixed N. Batching the masks against a succinct
representation of its windowed phase distribution is the open
evaluation step. Materializing the graph costs phi(N²) entries.

## Swap symmetry gives an exact interval parity target

There is a variant on modulus N itself. For an odd N, a public
unit c and interval I, count ordered pairs x,y in I with
xy=c² modulo N, taking residues once. Swap is an involution
preserving the window on both coordinates. Modulo two the count
therefore equals the number of a in I with a²=c² modulo N.
XOR away the contributions of the two public roots ±c. The
remaining bit is exactly the parity of the nonpublic roots in I.
Requiring *both* coordinates in I makes the cancellation exact.

On a singleton isolation trial from square-input-isolation.md the
corrected bit is one. More generally any odd interval suffices.
Halve it and query one child: corrected root parity is additive
on disjoint children, so the other child's parity follows from the
parent. Descend to a singleton in O(log N) queries. Its residue a
gives a proper divisor gcd(a-c,N) for odd N, since the two signs
partition the odd prime powers dividing N and a is nonpublic.
This handles repeated prime factors when there are at least two
distinct primes. The earlier isolation routine and perfect-power
preprocessing supply the remaining cases. This is a reduction to
an exact parity oracle for a different count, without a phase
digit or Fourier decoding. The swap does not itself evaluate the
parity in polynomial bit time. Direct enumeration still scales
with the window width.

## Deterministic starting interval with scale one

Fix c=1. For odd N having k>=2 distinct prime factors, CRT and
the two roots ±1 on each odd prime power give exactly 2^k roots
of a²=1 modulo N. Negation pairs them without fixed points. Thus
[1,(N-1)/2] contains exactly 2^(k-1) roots. Removing 1 leaves
the public interval [2,(N-1)/2] with odd parity, independently
of the factorization. Every root there yields a proper gcd(a-1,N).

Consequently exact inverse-pair parity on square boxes [L,U]^2,
2<=L<=U<=floor(N/2), suffices for deterministic splitting with
one query per bisection level. Indeed swaps cancel the cross terms
in parity when an interval is bisected. Even prefix squares
[2,B]^2 alone suffice: the endpoint prefix parities start at
zero (empty prefix) and one; keep the half whose endpoints differ.
The prefix parity can oscillate, so the invariant is differing
endpoints rather than monotonicity. Primality testing, removing
factors of two and perfect-power detection handle the other cases;
recursing on proper splits gives complete factorization with a
polynomial number of queries in log N. Direct parity scans still
have width-dependent bit cost.

## Exact prefix evaluation near either known endpoint

For an odd N with at least two distinct prime factors let F(B)
be the inverse-pair parity on [2,B], 1<=B<=H=(N-1)/2. From the
left endpoint enumerate integers k>=1 with 1+kN<=B² and count
the perfect squares among them modulo two. Each square has exactly
one positive root in [2,B]. From the right endpoint start at F(H)=1.
A root a in (B,H] corresponds bijectively to an odd positive
t=N-2a<N-2B with t²=4 modulo N; invertibility of 2 gives the
converse. Enumerate squares 4+kN for 0<=k and 4+kN<(N-2B)²,
counting only odd square roots, and XOR their parity with one.
There are no negative k for N>=15. Both procedures use integer
square roots on O(log N)-bit operands and need no factors of N.

Choosing the nearer endpoint in these coordinates costs
O(1+min(B²,(N-2B)²)/N) square tests. Thus both endpoint regions
of width O(sqrt(N)*poly(log N)) admit exact polylogarithmic bit
time. At B near N/4, including the first midpoint in the earlier
bisection, the cheaper bound still scales linearly in N.

## A fixed quarter box and nonadaptive dyadic scales

Suppose odd N has exactly two distinct prime divisors, with arbitrary
positive exponents. Its square roots of one are ±1 and ±u, where
1<u<N/2 picks a unique representative. Fix I=[1,floor(N/4)].
At scale c=2^j mod N count ordered pairs in I² with xy=c² mod N.
XOR out the known diagonal contributions ±c. By swap cancellation,
the resulting bit is one precisely when the fractional part of
2^j u/N lies in the first or last quarter of [0,1). No endpoint
ties occur because N is odd and u is a unit.

This bit says exactly that binary digits j+1 and j+2 of u/N agree.
The first digit is zero. With L=bit_length(N), scales j=0,...,L-2
recover the first L digits by keeping or flipping the preceding bit.
Their dyadic interval has length 2^(-L)<1/N, so ceiling its lower
endpoint times N recovers u. Its gcd with N after subtracting one
splits N. All boxes and scales are fixed from N before the answers
arrive. If N has more distinct primes the bit is instead a parity
over multiple hidden pairs, so this single-root decoder does not
apply. Uniform fast evaluation of even these fixed-box parities is open.

## Symmetric binary hashing for nonadaptive decoding on general inputs

For odd N with k>=2 distinct prime divisors there are M=2^(k-1)-1
nonpublic root pairs. Represent each by v=min(u,N-u), using L bits
where L=bit_length(N). Sample an affine binary hash H(v)=Av+b with
t output bits and independent fair entries. On both coordinates x
of the quarter box at scale c impose the same extra test: normalize
r=x/c mod N, require v=min(r,N-r)>1 and H(v)=0. The pair count
remains swap invariant. Its parity counts surviving diagonal roots
and the hash retains or discards both signs together at every scale.

Try t=1,...,L+1. At t=ceil(log2(2M)) the mean number of surviving
pairs is lambda=M/2^t in (1/4,1/2]. Pairwise independence gives
Pr(exactly one)>=E[X]-E[X(X-1)]>=lambda-lambda²>=3/16.
On such a trial all dyadic-scale bits decode that one root by the
previous argument. Generate the hash once per trial and keep it
across scales. Independently sampled trials can be batched before
any answers arrive; accept only a candidate yielding a verified
proper divisor. This uses a fixed outer box with an additional
polynomial-bit predicate inside it. The unfiltered arbitrary-prefix
reduction remains simpler for adaptive queries. No efficient oracle
for these hash-filtered counts is implied by efficient hash evaluation.

## Continued fractions can stop the quarter decoder early

For a recovered prefix of l binary digits of u/N, let x be the
midpoint of its dyadic cell. For odd N and a nonpublic square root
u of one, gcd(u-1,N) and gcd(u+1,N) are complementary coprime
factors. The reduced denominators of (u-1)/N and (u+1)/N are
therefore these same two factors in opposite order. Call the
smaller denominator D. Compute the continued fractions of both
public rationals x ± 1/N and inspect every convergent denominator
for a proper gcd with N. At 2^l>D² the midpoint error is strictly
less than 1/(2D²), so the corresponding reduced fraction is a
convergent by Legendre's criterion. A divisor is found at or before
that precision, without knowing D. The l digits use l-1 quarter
queries because the first binary digit is given. Checking after
each prefix allows early stopping while keeping a predetermined
sequence of query scales. This argument also applies to any one
root pair isolated by the hash variant. It changes query count,
not the computational cost of a quarter-parity evaluation.

## Short translated intervals via univariate small roots

For an interval centered at a public integer A, use the monic
quadratic f(T)=(A+T)^2-1 modulo N. Univariate Coppersmith finds
all offsets of size below a fixed constant times sqrt(N) in
polynomial bit time for fixed degree. A constant number of
overlapping translates handles an interval of width sqrt(N),
filtering endpoints and deduplicating recovered integers.
Consequently parity on any such interval, or a difference of
prefix parities across it, is computable from N without its
factorization. The argument also works with a public unit scale
by replacing 1 with c².

After l=ceil(bit_length(N)/2) binary digits of a hidden root u/N,
the dyadic cell confines u to an interval of width N/2^l<=sqrt(N).
Run the translated quadratic solver and verify candidates by gcd.
Thus at most l-1 quarter queries suffice on a singleton trial,
even for comparable factors, for which small-denominator recovery
alone can need nearly the full bit length. This is a theorem-based
application; the wide prefix from 2 to floor(N/4) still spans
Theta(sqrt(N)) such intervals if covered directly.

[Don Coppersmith, Finding Small Solutions to Small Degree Polynomials, sections 1–3](https://cr.yp.to/bib/2001/coppersmith.pdf)

### Fixed-dimension recovery diagnostic

The reproducible check_quarter_lattice.py uses 24 seeded semiprime
draws with primes near 2^10 to 2^11. Known factors generate the
target and an interval containing it. Recovery receives only N,
the center and radius. Its six lattice rows come from T^i f^j N^(2-j)
for i=0,1 and j=0,1,2, scaling column k by radius^k. It inspects
integer roots of the LLL-reduced row polynomials and verifies the
congruence and gcd. Planted-root recoveries were 24/24, 24/24 and
14/24 at radii floor(sqrt(N)/32), floor(sqrt(N)/8) and
floor(sqrt(N)/2) respectively. Each radius uses fresh offsets.
This fixed dimension is an empirical control; a miss cannot certify
that the interval is empty. The diagnostic supplies no prefix oracle.

### Certifying empty and nonempty intervals

For a reduced lattice row encoding g(T), the sum of absolute scaled
coefficients bounds |g(T)| throughout |T|<=radius. If this sum
is strictly below N², every modular root in the interval is an
integer zero of g: each generating row takes a multiple of N² at
such a root, and LLL takes integer linear combinations. The row
is nonzero since the original lattice has full rank. Listing its
integer roots and checking the original congruence thus certifies
completeness, including empty intervals, for the requested
nonpublic roots. If no row meets the inequality the result is
inconclusive.

On 48 seeded fixtures per radius (half planted, half randomly
centered), all 48 lists were certified at each of sqrt(N)/32 and
sqrt(N)/8, with 24 empty and 24 nonempty in each set. None were
certified at sqrt(N)/2. All certified lists matched an exhaustive
interval scan. These are fixed-dimension, small-input measurements.

### Bounded nonzero levels for a completeness certificate

Let S be the sum of absolute scaled coefficients of a nonconstant
reduced row g. Throughout the interval every modular root must obey
g(T)=j N² for an integer |j|<=floor(S/N²). Solving every integer
level in that range and checking the modular congruence certifies
the complete list. A constant row must be skipped since its matching
level would vanish identically. Setting a fixed maximum level K
bounds the number of degree-at-most-five equations by 2K+1 per
accepted row.

On exactly the same 48 widest intervals in the diagnostic, K=8
certified 14 empty and 20 nonempty lists compared with zero
certificates for K=0. The other 14 remained inconclusive. All
certified lists matched exhaustive enumeration. These results
concern small generated inputs and a fixed six-row basis.

### Level enumeration on a wide interval

Let g be any nonconstant integer combination of the six diagnostic
basis polynomials before column scaling. In degree d=2j+i, the
leading coefficient is a nonzero multiple of N^(2-j): the six
basis degrees are distinct and f is monic. Write W for the span
of the integer interval and h=floor(W/d). The d-th forward
difference on d+1 equally spaced points in the interval has
absolute value d! h^d |lc(g)|. The positive and negative weights
each sum to 2^(d-1), so the difference is bounded above by
2^(d-1) times the range (maximum minus minimum) of g there.
Consequently the range divided by N² is at least

    d! h^d / (2^(d-1) N^floor(d/2)).

For W=floor(N/4)-O(1) and 1<=d<=5 this grows at least linearly
in N. Even exact integer-domain extrema and an optimally chosen
row therefore leave Omega(N) candidate levels for any method that
enumerates every multiple of N² between the extrema. This argument
does not require every intervening multiple to be attained. Using
congruences to skip levels or computing their combined parity
directly would be a different method.

### Survivors of a periodic square sieve at the quarter cutoff

At B=floor(N/4) and scale one, the left-endpoint test has
V=floor((B²-1)/N) candidate levels k=1,...,V. For distinct odd
auxiliary primes r coprime to N, filter k by requiring 1+kN to
be a square modulo every r, including zero. Let Q be their product.
In each complete period of length Q the count of survivors is
exactly product_r ((r+1)/2): the affine map permutes residues at
each prime and CRT combines the tests. If Q<=V then the complete
periods alone retain at least V/(2*2^t) candidates, with t the
number of primes. Since t!<=Q<=V, t=O(log V/log log V).
Thus scanning surviving levels still takes V^(1-o(1)) visits
in this periodic sieve model. This includes primes selected in
response to N as long as the product bound holds. It makes no
survivor claim for Q>V and gives no lower bound for algorithms
that determine parity without visiting the survivors.

### A conditional exact sieve with polynomial description length

Assume Dirichlet GRH. For each positive nonsquare D=1+kN <= B²,
the Kronecker symbol (D/u) defines a nonprincipal quadratic
character on the units modulo 8D. Bach's GRH bound gives generators
of that group among primes coprime to 8D of size O(log²(8D)).
At least one generator r has character value -1. Thus D fails the
square-residue test modulo r, and r cannot divide N since D=1 mod r
for any r dividing N. A square passes every test, allowing residue
zero. For B<=N the same list of odd primes up to C log²(8N²), for
a sufficiently large absolute C, therefore singles out exactly the
square levels throughout the interval. This list is computable in
polynomial time in log N. Its combined period can exceed V.

This reduces the evaluator to the parity of an interval intersection
of explicitly given modular residue sets. Testing membership at one
level already takes polynomial time by integer square root. The open
step is to aggregate the intersection parity in polynomial time in
log N. The generator bound supplies no such aggregation algorithm.

Reference: Eric Bach, “Explicit bounds for primality testing and
related problems”, Mathematics of Computation 55 (1990), 355–380.

### Binary prefixes isolate the interval condition

Use n=ceil(log2 N) bits for the root representative a. Decompose
[0,B+1) into disjoint aligned binary blocks, one per set bit of
B+1, largest first. Each block fixes high bits and leaves the low
bits free, hence is affine over F_2. XOR the root-count parities
of the blocks and remove the known contribution at a=1 for B>=1.
This handles the cutoff with O(n) affine restrictions. For the
quarter cutoff the same decomposition works without knowing factors.

Even on one such block, linear equations over F_2 with auxiliary
bits cannot in general represent exactly the root witnesses under
projection onto their ordinary bits. For N=6545, B=1636, the block
[1024,1535] occurs in the decomposition of [0,1637). Exhaustive
testing of a² mod 6545 on those 512 inputs gives precisely
1121, 1189, 1429 as roots. Every nonempty projection of an affine
solution space is affine and has power-of-two cardinality. The
same obstruction holds if each valid input must have an odd number
of lifts and every invalid input an even number, because nonempty
fibers of a linear projection all have the same size. This example
does not exclude a representation of aggregate parity alone or
a nonlinear encoding.

### Degree of the prefix indicator

For an aligned block [L,L+2^d), write a=L+sum_i 2^i b_i.
Represent the root predicate by its unique multilinear polynomial
over F_2 in the d free bits. The block parity is exactly the
coefficient of the product of all d bits: a monomial on s bits
contributes 2^(d-s) ones to the truth-table sum. Thus a degree
bound strictly below d can certify even parity only.

The full degree occurs even on blocks containing just one nonpublic
root, with arbitrarily many free bits. Take H=2^d, d>=2, x=2H-2,
N=x²-1 and the block [H,2H). Every a here has 0<a²-1<2N:
the upper bound follows from (x+1)²-1<2(x²-1) for x>=6.
Hence divisibility forces a=x. Its low-bit pattern consists of
ones except at bit zero, giving indicator (1+b_0) times the
product of b_1 through b_(d-1). This block lies below the quarter
cutoff for d>=3. Its degree is d, of order log N, despite its
short factored expression. Exhaustive checks agree for d=2..12.
The example concerns arbitrary aligned subblocks and does not
assert balanced semiprimality or a circuit-size lower bound.

For general [L,U] with L>=2 the exact work of direct level scanning
is at most floor((U²-1)/N)-ceil((L²-1)/N)+1 integer square tests.
This provides a polynomial-time base case whenever that count is
polynomial in log N. It leaves long blocks at macroscopic L open.

### Quarter parity for moduli divisible by four

Suppose 4 divides N, gcd(c,N)=1, and count all solutions to
a²=c² mod N in [1,N/4]. For N>4 the maps a -> -a and
a -> a+N/2 preserve the roots and partition them into quartets.
Each quartet meets each open quarter once: the endpoints are
nonunits, and the four representatives have the form
r, N/2-r, N/2+r, N-r with 0<r<N/4. Hence the count is the total
number of roots divided by four, independently of the unit scale.
If the evaluator omits the public pair ±c, subtract its explicitly
testable contribution to the interval.

Write N=2^s m with m odd. CRT and multiplication by c give
2^omega(m) times 2 roots for s=2, or times 4 for s>=3,
using the standard square roots of unity modulo prime powers.
For s=2 and m>1 the quarter count is odd exactly when m is
a prime power. For s>=3 it is odd exactly when m=1. For N=4
the inclusive endpoint contributes one. Integer root extraction
and primality testing decide these cases in polynomial bit time.
Direct enumeration verified the parity for all multiples of four
through 512 and all unit scales 1<=c<min(N,18).

For odd N, any translation by t preserving the roots must send
both c and -c to roots. Subtracting the resulting congruences
gives 4ct=0 mod N, forcing t=0. This rules out the translation
step of this quartet argument on odd semiprimes.

### Doubling an odd modulus preserves the quarter-parity problem

Let Q(n,c) count roots of a²=c² mod n in [1,floor(n/4)].
For odd m>1 and unit d mod m, let c be the unique odd lift of
2d mod m to mod 2m. There is an exact count identity

    Q(2m,c) + Q(m,d) = 2^(omega(m)-1).

For each sign pair of roots at scale d modulo m, take its
representative r in [1,m/2). Doubling induces a bijection on sign
pairs onto those at scale c mod m. The new lower-half representative
is 2r for r<m/4 and m-2r for r>m/4. Only the latter is odd.
Oddness is precisely the extra condition for a root mod m to be
a root mod 2m at this odd scale. The first quarter of 2m spans
the lower half of m. Each original sign pair therefore contributes
exactly once across the two counts. CRT gives the right-hand side.

Modulo two the correction is one exactly for prime powers m,
which is decidable in polynomial bit time. The scale conversion
is invertible by halving modulo m. Consequently a polynomial-time
unit-scale quarter-parity evaluator on either family gives one on
the other with a single oracle call. When excluding the public
pair on both sides the right-hand side decreases by one, by the
same sign-pair correspondence. Verified the exact full-count
identity for all 13,230 unit scales at odd moduli 3 through 255.

### Square lifts of a unit-scale instance

For odd m>1, gcd(d,m)=1 and any positive h, set N=h²m and c=hd.
These choices need no coprimality between h and m. A solution has
h² dividing a², hence h dividing a. Substituting a=hb identifies
the lifted quarter count with the count of roots b²=d² mod m
in 1<=b<=floor(hm/4). Write R=2^omega(m) and Q=Q(m,d).
Each complete period contributes R. With j=h mod 4, the remaining
partial period contributes respectively 0, Q, R/2, R-Q for
j=0,1,2,3. The last two cases follow by reflection b -> m-b,
and no roots occur at the fractional endpoints since m is odd.

Thus for odd h the full lifted parity equals Q mod 2. For h
divisible by four it vanishes; for h=2 mod 4 it indicates that
m is a prime power. In particular h=2 places the problem on a
modulus divisible by four but erases the unknown quarter bit.
Any exclusion of the two public residues is handled on the lifted
instance by checking those residues against its actual cutoff.
Direct enumeration verified the exact counts in all 10,464
instances with odd 3<=m<=65, unit residues d, and 1<=h<=12.

### Arbitrary scales and the square part of a gcd

For any N>=2 and integer c, put D=gcd(N,c²), m=N/D. Write
D=s²t with t squarefree, and put h=st. Then multiplication by h
bijects the following finite sets, including when m=1:

    {1<=b<=floor(sm/4) : b²=(c/h)² mod m}
    {1<=a<=floor(N/4) : a²=c² mod N}.

Indeed D divides c² and must divide a², precisely forcing h|c and
h|a respectively. After substitution and division by D the
congruence has common coefficient t. For p|m the valuation of D
at p equals that of c² and is even, so p does not divide t and
c/h is a unit modulo m. Cancellation is therefore exact. The
cutoff follows from N/h=sm. If m=1 the count is simply floor(s/4).
Full counts can be adjusted for public roots by checking distinct
public residues in the original interval.

For promised squarefree N this reduction uses gcd arithmetic alone:
D=gcd(N,c), s=1, h=D, and the cutoff is exactly the original
quarter cutoff at the smaller modulus m. For unrestricted N the
formula requires the largest square divisor of D, for which the
identity supplies no polynomial-time algorithm. The unit case D=1
is unchanged. Verified the witness bijection for every residue c
at each N from 2 through 179, including zero scales.

### Jacobi symbol at uniform unit scales

For distinct odd primes N=pq, fix a nonpublic square root u of 1.
The nonpublic quarter bit H(c) is 1 iff the least absolute residue
of cu has magnitude <=N/4. Take c uniform among units. When N=3
mod 4, H is exactly balanced and independent of the Jacobi symbol
chi(c). Indeed H(c)=H(-c) while chi(-c)=-chi(c). Multiplication
by u permutes the units. The number of units in [1,N/4] is
floor(N/4)-floor(p/4)-floor(q/4)=phi(N)/4 in this case, by the
opposite residues of p and q modulo 4. Reflection doubles this
number for the folded interval. Hence each joint cell has phi(N)/4
members. Enumeration checked this identity among all applicable
pairs of distinct odd primes below 48.

For N=1 mod 4 the reflection doubles the Jacobi sum on [1,N/4].
In general E[(2H-1)chi(c)] equals
2 chi(u)(1+chi(-1)) sum_{1<=a<=N/4} chi(a) / phi(N).
The character is primitive modulo the odd squarefree N, so the
Pólya–Vinogradov inequality bounds this by O(log(N)/sqrt(N))
in absolute value. The constant predictor also has only O(1/N)
advantage by the same floor count. Every deterministic predictor
on a single sign is constant or plus/minus that sign, so this
bounds the best advantage from this symbol alone. This statement
is about uniform unit scales; it leaves adaptive scales and joint
features of c open.

### Several Jacobi probes at related scales

For unit c, every multiplicative probe chi(kc)=chi(k)chi(c).
Thus any number of fixed multipliers, including the dyadic sequence,
provides only the original sign. This also holds for multipliers
selected using only N, independent randomness and preceding signs.

Additive shifts require a different argument. For balanced distinct
odd primes p,q of comparable size, fix k>=1 public shifts t_i
pairwise distinct modulo each prime. For uniform unit c, a predictor
using only the vector of chi(c+t_i) has advantage at most

    O(k/sqrt(N) + 2^(k/2) (k+1)^2 log(N)/sqrt(N))

in predicting H(c). The bound is uniform in the shifts and in the
decision function, for a fixed balance ratio. In particular it tends
to zero exponentially in the input length for k <= (1-epsilon)log2 N
with any fixed positive epsilon (absorbing polynomial factors).

Proof: extend the folded quarter sign F=2H-1 to all residues by
using the same interval test on x=uc. Conditioning on units costs
O(1/sqrt(N)). A zero probe has probability O(k/sqrt(N)) by a union
bound; set arbitrary Boolean signs there. Expand the Boolean
decision function in Walsh characters on k signs. Its coefficient
l1 norm is at most 2^(k/2) by Parseval. Off the zero set each Walsh
monomial is a product of Jacobi symbols. To avoid multiplying the
zero-event cost by the Walsh norm, evaluate the multilinear Walsh
extension directly at the raw signs including zero: it is an average
over Boolean completions and stays in [-1,1].

For a nonempty subset S, change variables to x=uc. The sum involves
the Jacobi character of the product of x+u t_i over S, up to a
constant sign. At each prime this polynomial has |S| distinct
simple roots. The standard Weil bounds for complete multiplicative
and mixed additive-multiplicative character sums bound each local
Fourier sum by O((|S|+1)sqrt(p)) or its q analogue, including
frequency zero. CRT multiplies the local bounds. Completing the
folded quarter interval costs O(log N) in normalized Fourier l1
norm by the geometric-series bound; the same holds for F. Dividing
by N gives the stated bound per nonconstant Walsh monomial. The
constant monomial contributes O(1/N) by interval length. This
argument leaves shifts chosen from the input scale, growing adaptive
decision trees, and predictors using other features of c open.

### Adaptive shifted Jacobi probes with bounded depth

Keep balanced N=pq and a uniform unit scale. Allow the next shift
to depend on N, independent randomness and all prior Jacobi
answers. Assume along each nonzero-answer path the distinct shifts
are distinct modulo both primes (exact repeats may be cached).
A difference with a proper nonunit gcd would itself reveal a
factor. Suppose the nonzero-answer binary tree has L leaves and
depth k. Then the advantage for predicting the nonpublic quarter
bit using only the transcript is

    O(L (k+1)^2 log(N) / sqrt(N)).

Fix the randomness. Replace a zero answer by a fresh fair sign and
keep following the binary tree. This changes the prediction only
on inputs encountering a zero. A union bound over the at most L-1
binary query nodes gives probability O(L/sqrt(N)) for comparable
primes, even before conditioning on reaching any given node.
For a fixed leaf with prescribed signs e_i and shifts t_i its
probability conditional on c in the modified process is exactly

    product_i (1 + e_i chi(c+t_i))/2.

Its expansion has coefficient absolute sum one. The preceding
mixed character-sum estimate bounds its correlation with the
folded quarter sign by O((k+1)^2 log(N)/sqrt(N)), initially over
all residues. Sum this bound over leaves with their output signs.
Conditioning the resulting bounded predictor on unit inputs costs
O(1/sqrt(N)). Averaging over the independent randomness preserves
the estimate.

Since L<=2^k, depth at most (1/2-epsilon)log2(N) gives exponentially
small advantage in bit length. The leaf bound also applies to
polynomial-size binary trees of polynomial depth. It does not give
a comparable result for arbitrary polynomial-time branching
programs, whose unfolded trees can have exponentially many leaves.

### A signed return in the dyadic quarter sequence

For odd N>1 write b(x)=1 if the residue x mod N exceeds N/2 and
zero otherwise. For nonzero x let H(x) indicate that its least
absolute residue has magnitude below N/4. Then over F_2,

    H(x) = 1 + b(x) + b(2x).

This is the adjacent-digit identity above. Summing d consecutive
bits cancels the interior terms. Given a signed return
2^d = e mod N with e equal to +1 or -1, the XOR is therefore
(d mod 2) XOR [e=-1], for *every* nonzero starting residue.
Computing that output and verifying the supplied return need only
modular exponentiation and parity of d. Finding such a d efficiently
is a separate question.

In particular the two sign pairs of a unit-scale semiprime give
the same cycle parity, so the full quarter sequence (public and
nonpublic together) has XOR zero over a return. Dropping the public
pair gives the formula above. For any odd modulus the corresponding
aggregate simply repeats the per-pair parity once per sign pair.
This explains why the cycle checksum alone cannot distinguish the
hidden semiprime pair from the public one. A different return,
where 2^d squares to one but is neither sign of one, already splits
the modulus by gcd(2^d-1,N). Verified the formula at first signed
returns for all 6,566 units at odd moduli from 3 through 179.

### Distribution of short dyadic transcripts

Fix N=pq with distinct odd primes, a nonpublic root u of unity,
and a uniform unit scale c. Let W be k consecutive hidden quarter
bits starting at c. For every prescribed word w,

    | #{c : W=w} - phi(N)/2^k | <= 8.

Indeed x=cu mod N permutes the units. The k adjacent-digit equality
tests of x/N specify its first k+1 digits once its first digit is
chosen. The two choices give complementary binary prefixes. Each
prefix corresponds to a half-open interval of length N/2^(k+1)
in [0,N). Inclusion-exclusion over 1,p,q,pq counts its units;
each multiple count differs from interval length divided by its
divisor by at most one. Summing the two interval errors gives eight.
Thus the total variation distance from k fair independent bits is
at most 4*2^k/phi(N).

A fresh uniform unit and its public quarter sequence sample exactly
the same distribution without u. This reproduces the marginal law
but gives no coupling to the correct transcript at a specified c.
Histograms alone cannot verify an evaluator. Exhaustively checked
the count bound for all 45 semiprimes formed from odd primes
through 31, for k=1 through the bit length of each modulus.

### Joint law of public and hidden dyadic words at the same scale

For N=pq, distinct odd primes, write c=qi+pj mod N. Uniform
units correspond to independent nonzero residues i mod p, j mod q.
Choose the nonpublic u=1 mod p, -1 mod q. Then the public and
hidden fractional parts are i/p+j/q and i/p-j/q modulo one.

For k consecutive quarter bits at each coordinate, the joint law
has total variation distance from 2k fair bits at most
(1+2^(k+3))*(1/p+1/q). To see this, take independent continuous
uniform A,B on the circle. A+B,A-B are independent uniform:
the integer matrix is a surjective endomorphism of the two-torus
and hence preserves normalized Haar measure. Round A,B down to
the p- and q-grids. The two output coordinates each move by at
most delta=1/p+1/q in circular distance. The quarter word can
change only within delta of a multiple of 1/2^(k+1), by the
adjacent-digit description. A union bound over both coordinates
costs at most 2^(k+3)*delta. Conditioning the grid inputs to be
nonzero costs at most delta. The continuous words are independent
and uniform, since each word has two complementary prefixes.

For balanced primes and k below (1/2-epsilon)log2(N), observing
the public word gives exponentially small advantage for predicting
a specified hidden bit within the word, averaged over unit scales.
This is a restriction on the observed information; c itself is
available to the general evaluator. For p=1009,q=1013,k=4 exact
enumeration gave joint total variation distance 0.006657 (rounded).

### Additive shifts recover folded root multiples from disagreement counts

For N=pq with distinct odd primes, write |x|_N for the least
absolute residue magnitude, choose a nonpublic root u of unity,
and extend the hidden unit-scale bit to all c by
F(c)=[4|cu|_N<N]. On nonunits this extension equals the public
folded-quarter bit, since cu equals c or -c modulo N whenever
p or q divides c. In particular F(0)=1 by convention.

Let D(h) count all c mod N with F(c) different from F(c+h).
Then D(h)=2|hu|_N. Multiplication by u turns this into the size
of the symmetric difference of a circular interval with its shift.
The interval has 2 floor(N/4)+1 points, so the smaller of its
length and its complement length is (N-1)/2. For a shift of folded
length d <= (N-1)/2, exactly d points leave and d enter.

In particular D(1)/2 recovers the canonical hidden root and hence
a factor by gcd. More weakly D(2c) modulo 4 equals zero exactly
when F(c)=1: doubling the folded representative gives an even
representative in the first quarter, while folding back from the
second quarter subtracts an even number from odd N. Thus a fast
mod-four disagreement counter alone suffices for every unit-scale
bit in the dyadic decoder. D(h) modulo two always vanishes.
This reduction does not supply a fast mod-four counter.
Direct enumeration checked all shifts for the 21 semiprimes formed
from distinct odd primes through 19.

### Reflection shows the mod-four counter is point evaluation

For every Boolean F on residues modulo odd N with F(-x)=F(x),
let D_F(h) be the number of disagreements between F and its shift
by h. Then

    D_F(2c) mod 4 = 2 (F(0) XOR F(c)).

Write D_F(2c)/2 as the support size of F minus the overlap
sum_x F(x)F(x+2c). The support size has parity F(0) by negation.
Reflection x -> -2c-x preserves the overlap summand and fixes
exactly -c. Thus overlap parity is F(-c)F(c)=F(c). This also
follows by squaring the support polynomial in the binary cyclic
group algebra: reflection identifies autocorrelation with convolution
and Frobenius doubles every exponent.

For the extended hidden quarter bit F(0)=1, so the mod-four
disagreement oracle and the point oracle are interchangeable with
one call in either direction, since two is invertible modulo N.
The exact integer disagreement count still carries the magnitude
information in the previous identity. Reflection alone supplies no
new way to compute the surviving fixed-point value. Checked all
15,364 choices of odd modulus 3 through 17, symmetric Boolean
predicate and shift parameter in that range.

### Constant additive accuracy suffices for disagreement rates

For an odd semiprime suppose queries at arbitrary h approximate
D(h)/N to additive error 1/100. Halving gives estimates R(h) of
||hu/N||, distance to the nearest integer, within delta=1/200.
Clip the estimates to [0,1/2]. The following deterministic reduction
takes O(log N) queries and polynomial bit time.

Find a dyadic t with R(t) in [0.1,0.4]. Some true dyadic distance
enters [1/8,3/8] within O(log N) steps: above 3/8 the next distance
is below 1/4, and below 1/8 it doubles from at least 1/N.
Choose the global sign of u so that a=tu/N modulo one is positive
in the centered interval. Then |a-R(t)|<=delta and a is at least
0.095 from both endpoints 0,1/2.

For each h=2^j, j=0 through ceil(log2 N), choose the sign of R(h)
minimizing the discrepancy between ||R(t) +/- R(h)|| and R(t+h).
The signed estimate is within 7 delta on the circle of hu/N.
Indeed the correct sign has residual at most 3 delta. If the wrong
sign wins, the two true predictions differ by at most 6 delta
(each estimated prediction differs from its true value by at most
2 delta, and each residual is at most 3 delta; more directly the
wrong true prediction differs from the true observation by at most
2 delta + 3 delta + delta). Writing r=||hu/N||, the separation is
2 min(a,1/2-a,r,1/2-r), by folding a+r and a-r. Thus r lies within
3 delta of 0 or 1/2, and even the wrong signed estimate is within
7 delta of the correct phase. Use 7 delta as the uniform bound.

Starting from the last signed estimate, repeatedly take the half
(on the circle, two choices) nearest to the preceding signed
estimate. The two candidates are half a turn apart, and the true
branch is within (7 delta)/2 + 7 delta < 1/4 of the preceding
estimate. Each step halves the error. At j=0 the error is below
1/(2N); rounding after multiplication by N gives the root up to
global sign. A gcd gives a factor. All choices can use rational
arithmetic with O(log N) bits for fixed-precision oracle answers.
This reduces the precision requirement to a constant but still
requires input-specific estimates, not just a marginal average.

The reference implementation in check_disagreement_decoder.py takes
only N and a rate callback. It makes 2(bit_length(N)+1) calls.
A separate known-factor oracle tested 46 semiprimes, including
(2^61-1)(2^127-1), under four error profiles: rounding to the
nearest 1/50, fixed errors of either sign at 1/100, and seeded
mixed signs at that limit. All 184 runs recovered a verified
proper factor. The rounded mode exposes only 51 possible answers
per call. The implementation uses exact rational phase arithmetic.

### Feeding public Salié sums to the coarse decoder

The existing salie-quarter.md identity gives a second input path
for check_disagreement_decoder.py. At a unit shift h subtract
2 cos(2 pi h/N) from the real normalized Salié value and divide
by two. Clipping produces a cosine estimate z. The rate is
arccos(z)/pi. At nonunit h the extended rate is simply 2|h|_N/N.

An absolute error of 0.0001 in the normalized sum yields at most
0.00005 cosine error and 0.005 rate error. Indeed two angles in
[0,pi] separated by v differ in cosine by at least
2 sin²(v/2) >= 2(v/pi)². Rounding the rate to the nearest 0.01
still meets the 0.01 decoder tolerance (with additional precision
or slack when implementing the elementary functions numerically).

Run python check_salie_disagreement.py. Direct summation with
public inputs recovered a proper factor on each of seven moduli
from 15 to 10403, using floating point for this small diagnostic.
This is an end-to-end functional baseline with linear-scale
enumeration per query, not a polylogarithmic sum evaluator.

### Direct phase decoding tolerates 0.2 normalized Salié error

Write C(h)=cos(2 pi hu/N). An additive 0.2 approximation of the
real normalized sum on every requested unit shift yields an estimate
of C(h) within 0.1 by subtracting the public cosine and dividing
as above. Nonunit shifts have the public cosine. Clip estimates
to [-1,1]. The following reduction uses O(log N) such queries.

Choose a dyadic t with estimated cosine of absolute value <=0.81.
The earlier dyadic distance argument gives a shift with distance in
[1/8,3/8], where the true absolute cosine is at most sqrt(1/2).
Thus a qualifying shift exists. At the chosen shift the true cosine
has magnitude <=0.91, so the reference sine has magnitude >0.41.
Orient the global root sign to make that sine positive.

At a dyadic shift h use arccos of the estimated C(h) for the
unsigned angle. For its sign use the difference of estimates
C(h-t)-C(h+t), whose true value is twice the product of the two
sines. If the sign is correct, angular error is at most
arccos(0.9)<0.46 radians (the largest inverse-cosine change over
an interval of width 0.1 occurs at an endpoint). A wrong sign
requires |sin(2 pi hu/N)| <= 0.2/(2*0.41)<0.244. The true angle
then lies within arcsin(0.244) of 0 or pi modulo sign. The unsigned
estimate is within arccos(sqrt(1-0.244²)-0.1) of that same endpoint.
Their sum is below 0.77 radians. Either way the circular error is
strictly below pi/3, enough for nearest-half backtracking: the
correct half differs from the preceding estimate by less than
1.5*(pi/3)=pi/2, while the two halves are pi apart. At bit_length(N)
doublings the final error in turns is below 1/(2N). Rounding and
gcd finish the reduction. The strict margins permit fixed-precision
angle estimates followed by exact rational backtracking. This is a
pointwise guarantee distinct from the prior randomized RMS reduction.

check_salie_disagreement.py now also rounds the publicly summed real
value to the nearest multiple of 0.4 before querying the direct
phase decoder. Its interval [-4,4] permits 21 answers. The decoder
rounds angles to 1/4096 of a turn before exact rational backtracking.
Both public-sum paths recovered verified factors on 15, 77, 143,
323, 1147, 2021 and 10403. The earlier 184 rate-oracle simulations
also passed after sharing the backtracking routine.

### Kernel dimension for a clipped inverse involution

For unit c, index a binary vector space by the integers x in the
interval I. Define M e_x=e_y if x is a unit and y=c²/x mod N
lies in I, and zero otherwise (use representatives in [0,N)).
Set D=1+M. The exact number of roots in I is

    2 nullity(D) - nullity(D²).

A retained inverse pair gives a two-dimensional swap block and
contributes nullities (1,2); a fixed root contributes (1,1);
an unmatched or nonunit basis element contributes (0,0). Thus the
formula counts exactly fixed roots. In particular the parity is
nullity(D²) mod 2. But D²=1+M² is diagonal and its nullity counts
all retained elements, which is the earlier ordered inverse-pair
count. Simply squaring this operator returns us to the known
cancellation. Testing binary ranks for all 1,964 unit scales at
moduli 3 through 80 agreed with direct quarter counts.

The representation respects the interval and has one basis vector
per interval element. Constructing its columns on demand is cheap,
but the quarter matrix has floor(N/4) columns. A genuine advance
here would compress the aggregate into polynomial dimension or
give a polynomial-bit nullity algorithm for this specific operator.

### Rank of the translated-window answers across square inputs

Can all translations of a single window share a short binary feature
vector for each unit square y? Suppose the answer at translation t
is a binary dot product of that vector with a vector depending on
t and N alone. For odd N>1 and 1<=w<N coprime to N, the required
dimension is at least the number of distinct unit squares mod N.

Indeed the root indicator columns v_y have disjoint nonempty
supports and vanish at zero. Applying the cyclic window-sum map W
to each gives the answer matrix. If Wv=0, adding consecutive window
answers gives v(t)=v(t+w) for every t. Coprimality makes v constant.
Thus W is injective on the column span, since that span vanishes at
zero. The answer matrix has full column rank over F_2. A shared
dot-product representation of dimension d factors it through F_2^d,
so d bounds its rank. In particular w=(N-1)/2 qualifies, and for
distinct odd primes N=pq the rank is (p-1)(q-1)/4. Binary elimination
verified ranks 2,3,5,6,15,30,72 at N=15,21,33,35,77,143,323.

This applies to full root counts varying both y and translation.
A single fixed cutoff, the general splitting reduction at fixed
scale, nonlinear readout, and a kernel built jointly from the query
parameters are different models. In particular it gives no lower
bound on their evaluation time.

### A two-dimensional quarter determinant over a large number field

For distinct odd primes N=pq and a unit h, put A=2cos(2 pi h/N),
B=2cos(2 pi hu/N) for a mixed root u of one. Each positive
eigenvalue among A,B contributes one first-quarter root. Neither
vanishes, so sign(AB)=(-1)^Q for the full quarter count Q.
A real two-dimensional companion matrix with characteristic
polynomial X²-T(h)X+D has the desired determinant sign. Here
T(h)=A+B is the normalized Salié sum and
D=AB=(T(h)²-T(2h)-4)/2 by double-angle. This is a determinant-sign
encoding, not a nullity formula. The entries have public sum
expressions; the known direct evaluation still scans the units.

In particular D itself has rational algebraic degree phi(N)/4.
Indeed AB is the sum of the two cosine pairs at frequencies
h(1+u) and h(1-u); these have respective conductors p,q in some
order. Write their values as x in Q(zeta_p) and z in Q(zeta_q).
The two cyclotomic fields have intersection Q since their
compositum has degree (p-1)(q-1). An automorphism fixing x+z
has its change in x equal to minus its change in z, hence rational.
A finite-order automorphism cannot translate x by a nonzero rational
(apply it repeatedly). Thus both terms are fixed separately.
Within each prime field the stabilizer of a primitive cosine pair
consists exactly of the signs, since the two roots of the associated
quadratic are the primitive root of unity and its inverse. The
stabilizer in the product Galois group has size four, proving the
degree claim. Consequently explicitly expanding a dense rational
minimal polynomial for D would take linear-scale output. This
does not constrain a succinct expression or direct sign test.

### Extending the quarter determinant to any odd modulus

For odd N>1, unit h and m=2^(omega(N)-1) opposite pairs of roots
of unity, take one representative r per pair and let lambda_r be
cos(2 pi h r/N). The quarter count Q is the number of positive
lambda_r. Therefore its parity is read from the sign of their
product and the parity of m. Write t_j=sum_r cos(2 pi j h r/N),
with t_0=m. For 1<=j<=m and gcd(j,N)=1, twice t_j is the public
normalized Salie sum at scale jh. A proper gcd(j,N) itself splits
the modulus in a factoring reduction. Here m<N, so no j in this
range is a multiple of N. This qualification matters: a nonunit
scale cannot simply be substituted into the unit-scale formula.

Given m and t_1,...,t_m, the i-th power sum of the lambdas is
2^(-i) sum_{j=0}^i binom(i,j)t_|i-2j| by expanding cosine.
Newton identities then give the product e_m in O(m²) arithmetic
operations. Absolute errors delta in the t_j propagate to at most
delta in each power sum since the displayed weights sum to one.
To bound the final error, use the generating series for elementary
symmetric functions: perturbing the power sums multiplies it,
through degree m, by exp(sum_{i=1}^m (-1)^(i-1) Delta_i z^i/i).
The sum of absolute coefficients of the perturbing factor minus
one is bounded by exp(delta H_m)-1, with H_m the harmonic sum.
The original polynomial has coefficient l1 norm <=2^m. Each
|lambda_r|>=sin(pi/(2N))>=1/N because N is odd. Consequently
delta=(4N)^(-m)/(4m) suffices to preserve the product sign.
Dyadic input rounding to this tolerance followed by exact rational
arithmetic uses polynomially many bits in m and log N (denominators
from binomial averaging and Newton divisions also have polynomial
bit length). Numerical checks covered 9,45,105,315,385,1155 and
all unit h below min(N,30), using root lists as reference data.

This is a conditional postprocessing bound. It needs the pair
count m or a promise fixing it, plus the indicated traces. For
omega(N)=O(log log N) the processing and requested precision are
polynomial in log N. At unrestricted omega(N) even this matrix
dimension can exceed any polynomial in the input bit length.

### Recovering the missing zeroth trace from a pair-count bound

Given an integer upper bound M>=m, use the preceding t_j for
1<=j<=2M. If 2M>=N, direct enumeration already costs polynomial
time in M. Otherwise a nonunit j yields a proper gcd in the
splitting reduction. With the actual cosine traces available,
build a real symmetric Toeplitz matrix B of size 2M+1 with zero
diagonal and off-diagonal entry t_|i-j|. Then

    m = -lambda_min(B).

Indeed B+mI is one half of V V*, where V has one column for each
of the 2m roots x at scale h and entry exp(2 pi i jx/N) in row j.
It is positive semidefinite and has a kernel since there are more
rows than columns. Removing the common diagonal shifts its lowest
eigenvalue to -m. This proof needs no spacing estimate for the roots.

Absolute trace errors <=1/(16M) perturb any row sum in absolute
entry norm by <=1/8 and hence perturb the lowest eigenvalue by
<=1/8. Computing that eigenvalue to another 1/8 and rounding its
negative recovers m exactly. Thus a known polynomial upper bound
replaces the need for the exact number of distinct primes in this
postprocessing step. The determinant sign still requires finer
traces as above. Tests with seeded perturbations recovered m in
93 unit-scale instances across 9,45,105,315,385,1155,5005, using
root-enumerated traces. No public fast trace evaluator is assumed.

### Randomized linear summaries for translated windows

The exact shared-rank obstruction above changes if each fixed query
may fail with probability delta. Let a square have at most R roots,
and use n=bit_length(N) bits (including for endpoint N). At each bit
position j choose a random binary matrix mapping the more significant
prefix to b bits, with 2^b>=2Rn/delta. Store one parity per (j,bucket):
the number mod 2 of roots x with x_j=0 whose higher prefix hashes to
that bucket. This vector depends on the square and the public coins.

For cutoff z, XOR the bucket of its higher prefix whenever z_j=1.
With exact prefix comparisons each root x<z contributes at its first
differing bit; all others contribute zero. Each false prefix match
has probability 2^-b, since the difference is nonzero. For [a,b')
XOR the answers for both endpoints. A union bound over at most 2Rn
comparisons gives error <=delta per fixed query. Wrapping windows
also need the total parity, one extra coordinate. Thus the binary
dot product uses O(R n²/delta) coordinates, polynomial for a unit
square modulo a semiprime at constant error. Fresh coins handle
adaptive queries; correctness for all translations simultaneously
is a stronger requirement. This differs from symmetric hashing for
root isolation: the target here is the parity of the whole interval.

Each bucket is an aggregate over roots satisfying affine bit tests.
Forming it from an explicit root list is easy; an efficient way to
obtain it directly from the square would be needed for factoring.
A seeded root-list diagnostic on N=101*103 with 2048 buckets per
position had 27 errors in 2000 random interval queries (union bound
112/2048 per query). This tests decoding only.

### A single public binary prefix suffices on distinct odd semiprimes

This specializes the existing interval Fourier reduction to remove
the arbitrary endpoint from the evaluator. Let n=bit_length(N),
B=2^(n-3), N=pq for distinct odd primes. Ask for the parity of
roots of h² with ordinary n-bit representative starting in 000,
i.e. 0<=x<B. For unit h, remove the publicly testable contribution
of ±h. The remaining bit is [|uh|_N<B] for a mixed root u.
For nonunit h this last bit is publicly computable because uh=±h.
Thus exact single-prefix evaluations supply the full sign function
2[|uh|_N<B]-1 on the cyclic group.

At frequency u the magnitude of its normalized Fourier coefficient
is 2 sin(pi (2B-1)/N)/(N sin(pi/N)) by a geometric sum. This is
uniformly bounded below: N>=15 and 1/4-1/N <= (2B-1)/N < 1/2,
so the magnitude exceeds 0.3. The existing sampled heavy-frequency
search lists u and gcd(u-1,N) splits the modulus. This also tolerates
a predictor with sufficiently small constant error at uniform unit
shifts: an error rate rho changes the coefficient by at most 2 rho.
An arbitrary tiny advantage over chance is insufficient for this
argument because the shorter interval can give biased labels.

Hence the outstanding evaluator may focus on a deterministic
codimension-three binary subspace intersected with the modular
root condition; there are no prefix hashes or variable thresholds
in this semiprime sufficiency statement. This remains a reduction,
not a fast method to compute that intersection parity. Exhaustive
unit-shift checks on 45 products of distinct odd primes through
31 verified removal of the public pair; the smallest fundamental
coefficient magnitude there was 0.376945.

### Splitting the free bits of the 000 prefix

Let R=F_2[Z]/(Z^N-1). On the prefix 000 let the remaining d=n-3
bits give x=i+2^a j, with 0<=i<2^a and 0<=j<2^(d-a). Summing
Z^(x²) in R and reading coefficient y gives the desired parity.
Could the individual monomials first be represented as a sum of
a few products of separate low- and high-bit functions over R?

Any exact pointwise decomposition requires at least
min(2^a,2^(d-a)) products. Evaluate Z at a primitive N-th root
zeta in an algebraic closure of F_2 (N is odd). Each separate
product specializes to a rank-one matrix at most. But after
removing nonzero row and column factors the matrix has entries
(zeta^(2^(a+1)j))^i. Its nodes are distinct since 2 is coprime to
N and the number of columns is at most N. Its leading square minor
of the claimed size is Vandermonde and invertible. A balanced
split thus requires order sqrt(N) products in this representation.

For a *specified* unit square y on a distinct odd semiprime, the
binary coefficient matrix [x²=y mod N] on the same split has at
most two ones: x<B<N/2 and the four roots form two opposite pairs.
Thus its rank is at most two. This low rank by itself gives no
efficient way to locate or count its entries. The argument rules
out compressing every pointwise phase across this split before
coefficient extraction; it gives no lower bound on computing the
single summed coefficient, even via a different phase-sum identity.

### Existence alone on prefix 0000 suffices for a semiprime split

Let N=pq for distinct odd primes, n=bit_length(N), B=2^(n-4).
Suppose an evaluator answers only whether a unit square y has any
root below B. No count or parity is supplied. Put A(h)=[|h|_N<B]
and let u be a mixed root of one. The answer on y=h² is
F(h)=A(h) OR A(uh). Extend F by A(h) at nonunits, where uh=±h.
This extension is computable and invariant under h -> uh.

For uniform h on all residues let alpha=(2B-1)/N. Since A<=F
and E F<=2 alpha, the covariance with the public interval satisfies
Cov(F,A)=alpha(1-E F)>=alpha(1-2 alpha). Here 1/16<=alpha<=1/4
for N>=15, so this is at least 1/32. Centering removes the constant
Fourier mode. The normalized interval coefficients at nonzero
centered frequencies a have magnitude <=1/|a| by the geometric
sum. By Parseval and Cauchy-Schwarz the correlation tail past K is
<=sqrt(2/K). Fix e.g. K=32768, making the tail <1/64. For N<=2K+1
use trial division. Otherwise some 0<|a|<=K has Fourier magnitude
for F at least 1/(128(1+log K)), by the harmonic sum bound.
Invariance gives the same magnitude at au. If gcd(a,N)>1, trial
division through K splits N; otherwise gcd(au-a,N) is proper.
The existing sampled Fourier list search followed by gcd checks
for all small offsets therefore yields randomized polynomial-time
splitting given a polynomial-time existence test.

The fourth zero bit keeps the union probability bounded away from
one without assumptions on the overlap of the two root pairs.
A public polynomial-size kernel with positive nullity precisely
on these nonempty fibers would suffice. The unresolved part is
building such a kernel. Enumeration at all residues of the 45
semiprimes formed from distinct odd primes through 31 verified
the decision-function identity, including its nonunit extension.

### Removing two-cycle kernel vectors for the existence test

Use the clipped inverse matrix above with a sign on its edges.
For unit y set z=y/x mod N on unit x in the prefix. Put J e_x=e_z
when z is retained and x<=z, J e_x=-e_z when retained and x>z,
and zero on every other column. Over F_3 or the reals let L=I-J.
An unmatched index contributes [1], a fixed root contributes [0],
and a retained two-cycle contributes a 2-by-2 block of determinant
2. Thus the kernel has precisely the coordinate vectors of the
retained roots as a basis. The use of odd characteristic matters.

Over the reals L^T L is diagonal with respective entries 1,0,2.
This yields an exact unit gap between zero and the positive singular
values. It also exhibits the limitation: the root projector
(I-L^T L)(2I-L^T L)/2 is just the diagonal of pointwise root
indicators. The sign removes spurious vectors but leaves one
coordinate per candidate; knowing how to compute a column or the
spectral gap does not supply a small global kernel. A direct diagonal
encoding of the root predicate would have the same dimension. Checked
the Gram identity and zero coordinates at 310 unit inputs across
seven odd moduli, including nonsquare unit inputs.

### Choosing the zero-prefix length with more prime factors

Suppose the existence evaluator accepts a prefix length as input.
For odd N with s>=2 distinct prime factors, including repeated
powers, set R=2^s and B=2^(n-s-2) for n=bit_length(N). This
requests s+2 leading zeros, giving alpha=B/N between 1/(4R) and
1/(2R). The exponent is nonnegative since N>=3*5*3^(s-2).
For a proof on all residues put A(h)=[0<=h mod N<B] and
F(h)=max_u A(uh), with u over square roots of unity. At units
F is exactly the requested existence bit. At zero F=1; if a
sampled nonzero h is a nonunit, its gcd already finishes the split.
Thus the Fourier search can simulate this full-grid signal or stop
successfully, without having to solve a nonunit existence query.

The same covariance argument now gives Cov(F,A)>=alpha(1-R alpha)
>=1/(8R). The interval Fourier tail beyond K costs <=sqrt(2/K).
Choose K a sufficiently large multiple of R². For N<=2K+1 use
trial division; otherwise a nonzero |a|<=K has predictor Fourier
magnitude Omega(1/(R log(2K))). If a is a nonunit, trying small
divisors succeeds. For unit a, invariance gives the same magnitude
at au for any mixed root u; gcd(au-a,N) then splits N. Searching
the heavy list and trying offsets takes time polynomial in R and
log N with a unit-square prefix-existence oracle.

No supplied s is needed. Try guesses t=2,3,... with prefix length
t+2 and the corresponding search threshold and explicit work cap.
For an expected-time guarantee, dovetail repeated independent
bounded runs for guesses: at stage j run each guess t<=j for
j-t+1 fresh trials (skip invalid prefix lengths). Raise confidence
per trial at the true t to a sufficiently large constant, or note
that failure through stage s+k is exp(-Omega(k²)); even exponentially
growing per-stage budgets then have finite expected overhead.
Verify every returned gcd. On the promised inputs this costs
poly(log N,2^s) in expectation assuming polynomial-time existence
queries. Thus the prefix length is logarithmic in bit length when
s=O(log log N). The reduction gives no efficient existence evaluator.
Checks of the orbit mask, unit existence bit and covariance bound
covered five moduli from 45 through 15015, including repeated powers.

### Recovering a root by changing the query modulus

The exponential dependence on the number of prime factors in the
same-modulus existence reduction disappears if a uniform oracle
can be queried on auxiliary odd moduli with variable prefix length.
Remove factors 2 and 3 first. For a unit square y modulo the remaining
N>1, set B=2^bit_length(N). Seek a root x<B with 3 not dividing x.
One exists since the opposite representatives r,N-r are both below
B and cannot both be multiples of 3.

Maintain a unit a modulo d=3^j such that an eligible x is congruent
to ±a modulo d. Start at d=3,a=1. The three classes up to sign
above it are represented by b=a+td, 0<=t<3. For each candidate
form M=3dN and, by CRT, the unit square Y congruent to y modulo N
and b² modulo 3d. Query for a root of Y below the *same absolute*
B by requiring bit_length(M)-bit_length(N) leading zeros. Since
an odd prime power has only ±b as roots of b² for unit b, a yes
answer is exactly feasibility of this digit. Take the first yes
and update a=b,d=3d. Once d>2B, exactly one representative from
±a can lie in [0,B); min(a,d-a) is that root. Verify its square.
This takes at most three exact decisions per ternary digit, hence
O(log N) calls on moduli below 12N², plus polynomial bit work.

This deterministic root selection uses only N,y. When y=h² for a
uniform unit h, comparison with the selected root x splits N by
gcd(h-x,N) with probability 1-2/2^omega(N), at least one half
when N has two or more distinct prime factors. For each prime
power the sign either matches or differs, independently. Primality
and perfect-power tests and recursion handle all integer inputs.
This is an all-input randomized polynomial-time reduction conditional
on the existence evaluator. It does not justify replacing that oracle
by the fixed four-bit oracle, or by a splitter which simply returns
the known factor 3 of an auxiliary modulus. The executable
check_prefix_root_recovery.py verified all 181 unit squares across
six moduli using a direct, exhaustive prefix-existence oracle.

### A probability bound sufficient for binary root recovery

The binary decoder only needs one occupied terminal progression. Let L be its depth, and let a polynomial-time predictor assign complementary rational probabilities to the two children of each visited node using only public inputs. Sample a complete path according to those probabilities and run the existing short-root solver at its leaf. Check any returned root by squaring modulo N. If the total probability of occupied leaves is at least 1/P(log N) for every promised input, independent trials give expected polynomial time when P is polynomial. The predictor may be wrong at particular nodes. A single occupied path with cumulative negative log probability O(log log N) suffices, as does the same total mass spread over many paths.

There is also a deterministic version under the single-path guarantee. For a fixed threshold 1/K, expand only prefixes whose probability product is at least 1/K. There are at most K such prefixes per depth because their probabilities sum to at most one. Every prefix of a leaf of weight at least 1/K survives. Thus O(K(L+1)) node visits suffice, including the terminal root checks. Doubling K removes the need to know the threshold in advance. Exact products of L rational predictions of polynomial bit length still have polynomial bit length.

A constant per-digit probability of following a particular occupied path can yield exponentially small path mass. The required estimate is on the complete path or on total occupied mass. A public predictor meeting either bound on general inputs remains open.

### How much target information must reach the leaf selector?

For each fixed odd N, it suffices to bound the mean success probability on a uniform unit square. For a given unit square y, sample a fresh uniform unit u and apply the predictor to yu² modulo N. Squaring maps the unit group onto its square subgroup with equal fiber sizes, so this input has precisely the required distribution. A verified root z gives the root z/u modulo N. Independent restarts attain expected inverse-mean trial count. This standard square blinding only randomizes the target at the given modulus.

Suppose the selector sees the target only through a summary with at most K values. Its distribution over terminal leaves can depend arbitrarily on that summary and N. Write X for the decoder's interval limit. A terminal leaf contains at most 2(X+1) integer candidates, counting both signs. Each candidate squares to one target. Thus for any fixed summary and leaf, at most 2(X+1) targets with that summary can succeed. Average over the leaf distribution for that summary, then sum over summaries and divide by the number Q of unit squares. Mean success is at most min(1, 2K(X+1)/Q). This bound concerns leaf selection; the final lattice solver still gets the entire target.

For N=pq with distinct odd primes, Q=(p−1)(q−1)/4 is of order N. The implemented X is of order N^(2/5), giving mean success O(K/N^(3/5)). In particular, inverse-polylogarithmic success requires at least (3/5)log₂ N − O(log log N) bits of target summary capacity. This counting condition does not bound the cost of computing a sufficiently rich summary or establish its sufficiency.

### Interval cost of removing a known unit scale

For a public unit h and target h² modulo odd N, multiplication by h sends the roots of 1 bijectively to the target roots. A cyclic interval I of width 0<w<N becomes the selector S={u:hu mod N belongs to I}. Its exact number of cyclic runs in canonical u order is min(h,N−h,w,N−w), taking h in [1,N). A run begins when hu is in I and hu−h is outside. Multiplication permutes residues, so the number of starts is |I minus (I+h)|. For an interval, this is min(h,N−h,w,N−w): replace I by its complement or reverse the shift without changing the count, then both the width and shift are at most N/2 and the overlap has size max(0,w−h). Enumeration verified the formula for every unit and nontrivial width on odd moduli below 80; the interval origin does not affect it.

For small positive h one can generate the pullback of a nonwrapping interval [L,U) directly: for j=0 through h−1 take [ceil((jN+L)/h),ceil((jN+U)/h)), clipping to [0,N). Negative small representatives use reflection. Ordinary intervals decompose into O(log N) binary blocks; blocks above the lower half can first be reflected using root symmetry. Thus polylogarithmic signed scale permits polynomially many fixed-target prefix counts. If the signed scale and both interval widths w,N−w are of order N, a representation of the selector itself by intervals requires order N runs, even if XOR of intervals is allowed (each contributes at most two boundaries). This bound concerns pointwise selector conversion; a formula preserving just the root parity may be shorter.

### Wide intervals with a short exceptional set after reflection

For odd N>1, a unit square target y, and a cyclic interval I in canonical residue order, put H=[1,(N−1)/2] and S={x in H: exactly one of x,N−x lies in I}. Summing each root pair's contribution gives the parity of the roots in S. Alternatively use H minus S and toggle its root parity exactly when N is a prime power: H contains 2^(ω(N)−1) roots, and prime-power recognition is polynomial in bit length. This requires only the promise that y is a unit square. The sets and their lengths can be formed from the endpoints of I and its reflection, splitting wrapping intervals at zero and sorting the resulting constant number of boundaries.

Let L be the smaller of |S| and |H minus S|. Translate the quadratic (a+t)²−y on each component of the smaller set and run the existing square-root-scale solver in blocks. This evaluates the original interval parity in polynomial bit time when L is at most √N times a fixed polynomial in log N, regardless of the width of I or the size of a supplied scale. Both parity identities matched exhaustive counts on 88,406 (N,y,I) triples for odd N from 3 through 37. For a block wholly inside H the criterion is simply min(|I|,|H|−|I|), so it still leaves blocks occupying about half of H outside this guarantee.

### Frequency of short reflection exceptions

Write H=(N−1)/2 for the number of opposite pairs of nonzero residues, I=[a,a+w) cyclically, and t=min(w,N−w). Its reflection starts at 1−a−w. Thus with d=(1−2a−w) mod N the disagreement count on the lower half is exactly s=min(t,d,N−d). Indeed the full symmetric difference has twice the size of I minus its translate by d, and zero is never in the difference. Uniform start a makes d uniform since N is odd.

For an integer threshold K with H>2K, the number of starts satisfying min(s,H−s)≤K is N if t≤K, is 2K+1 if K<t<H−K, and is 4K+3 if t>K and t≥H−K. In the middle case only d at circular distance at most K contributes. In the last case distances H−K through H add 2(K+1) starts. Exhaustive membership tests checked every start and width for odd N from 25 through 101 with K=floor(√N).

Consequently the implemented criterion covers a fraction (2 floor(√N)+1)/N of uniform cyclic starts at quarter width for sufficiently large N. For intervals wholly within the lower half and away from zero, s=w, so translating such an interval within that half cannot improve its classification. This is coverage of the reflection shortcut, not an estimate of the distribution of queries reached by the adaptive decoder.

### Combining reflection-evaluable selectors at several scales

Fix the H=(N−1)/2 opposite pairs of nonzero residues. A folded selector passing the reflection test with threshold K is constant outside at most K pairs. Pulling it back by any unit multiplier preserves that property because multiplication permutes these pairs. For t such selectors, every point outside the union of their exceptional sets sees the same t-bit vector. Therefore every pointwise Boolean combination of them is constant outside at most tK pairs. A target selector accepting w pairs requires t≥ceil(min(w,H−w)/K). This covers XOR formulas using the known parity of the whole half as a free constant. For a lower-half prefix of the pair representatives [1,w], the bound is attained: partition the shorter of that prefix and its suffix complement into contiguous blocks of at most K and XOR their counts, using the total parity for the complement. These are cyclic-interval queries, with no dyadic alignment restriction in the upper bound.

For w near H/2 and K=floor(√N), even changing scales and complementing requires order √N windows to reproduce the selector pointwise. A count identity using special properties of the actual square roots can evade this argument: it only needs the parity on the root set to agree.

### Two evaluations classify polynomial maps on the roots of one

Let N>1 be odd. To test a public polynomial P over ℤ/Nℤ as a map from the roots of 1 to themselves, evaluate v=P(1) and w=P(−1). The test v²=w²=1 mod N is necessary and sufficient. At any root r the idempotents (1±r)/2 are orthogonal with sum one, giving P(r)=(1+r)v/2+(1−r)w/2. Squaring this expression proves sufficiency. If v or w differs from both ±1, gcd(v−1,N) or gcd(w−1,N), respectively, is a proper factor. Otherwise the four possible pairs give only the two constant maps and r↦±r. In particular any different root permutation provided by a polynomial circuit immediately splits N by two public evaluations. This also applies to rational circuits whose division gates are units at both test inputs: their values at a root are the same idempotent combinations, and inverses glue componentwise. The argument includes repeated odd prime powers.

This test concerns maps preserving the entire root set. Bit-dependent selection, a map only promised on some hidden roots, or an aggregate parity identity can require different information.

### Quadratic products lose prefix parity at bounded precision

For distinct odd primes with N=pq, omit the known root 1 and form
B(t)=∏_{x=2}^t (x²−1). Every block of q consecutive integers contains
an x equal to 1 modulo q, and likewise for p. Thus B(t) is zero modulo
N^k once t−1 ≥ k max(p,q). Its residue then stays fixed even if a
mixed root enters the prefix.

For an explicit infinite family of collisions, take p and q of comparable
size tending to infinity, and k polynomial in log N. Let r be the larger
representative in the opposite pair of mixed roots. Then r>N/2, so both
B(r−1) and B(r) vanish modulo N^k for sufficiently large N. The two
prefix parities differ. For N=77 and k=1, the same collision already
occurs at t=33 and 34, since 34²−1=15·77. This rules out reading parity
from the product residue alone. It leaves open methods that also use the
endpoint, subdivide the interval, or recover more than the truncated product.

For the marked product F_t(Z)=∏_{x=2}^t (Z+x²−1), the same obstruction
covers coefficients through degree d. If t−1 ≥ (k+d) max(p,q), then
at least k+d constant factors are divisible by p, and likewise by q.
A term of degree j ≤ d omits only j constants. Every such term is
divisible by N^k. Taking t=r−1 and t=r as above gives opposite parities
with identical low coefficients for polynomial-in-log-N choices of k,d.

The first surviving coefficient can itself encode a factor. On the full
residue system let G(Z)=∏_{x=0}^{N−1}(Z+x²−1). Modulo p its order
of vanishing at zero is exactly 2q: there are q copies of each sign
class and every remaining constant is a unit modulo p. Modulo q the
order is 2p. Thus for p<q the first nonzero coefficient modulo N is
at degree 2p, and its gcd with N is p. This gives an exact reduction
from semiprime factoring to finding that degree for this implicitly
specified product. Scanning to it costs order √N on balanced inputs.

### A single value of the full product splits a semiprime

For N=pq with distinct odd primes, choose a unit u with Jacobi symbol
(u/N)=−1. Exactly half the units qualify and the Jacobi test uses only
N. Set Z=1−u in G. A prime factor of N divides G(1−u) exactly when
x²=u has a solution modulo that prime. The two Legendre symbols are
opposite, so gcd(G(1−u),N) returns the prime where u is a square.

Conversely, given p,q one can evaluate every point in polynomial time.
In F_p the one-period product ∏(x²−u) is
−u (u^((p−1)/2)−1)², since every nonzero square occurs twice.
Raise this to q for the full product modulo p, interchange p and q
for the other residue, and apply CRT. This proves randomized
polynomial-time equivalence of full-product point evaluation and
factoring on distinct odd semiprimes. It removes the coefficient scan
from the reduction but supplies no public fast product algorithm.

More generally, for any odd N with h≥2 distinct prime factors, a random
unit gives independent uniform local Legendre signs by CRT. A mixed
sign pattern occurs with probability 1−2^(1−h) and guarantees a proper
gcd, even with repeated prime factors. Sampling a nonunit already
yields a factor when sampling from 1 through N−1. Prime powers can be
handled by exact integer root and primality tests.
