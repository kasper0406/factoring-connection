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

### Empty corners in the inverse-pair box

For the second-quarter interval, test rectangles near its endpoints
N/b and N/d for b,d in {2,4}. If xy≡1 mod odd N>16, then
(bx−N)(dy−N)≡bd mod N. Both factors on the left are odd, while bd
is even. Consequently their product cannot equal bd as integers,
and its absolute value is at least N−bd. A whole rectangle can be
discarded whenever the product of its two maximum absolute endpoint
offsets in these scaled coordinates is smaller than N−bd. This works
directly on public interval bounds. Across the 55 pairs of distinct
primes from 7 through 43, all 7,156 combinations of endpoint choices
and inverse pairs inside the discrete second quarter passed the bound.

This batch test still leaves the interior square [N/3,3N/8]² for
large N: all four endpoint offsets there are at least N/4 in
absolute value. In particular a row scan using just these corner
exclusions still has Θ(N) interior rows to handle. The test prunes
rectangles but supplies no asymptotic improvement for the wide query.

### Coverage of rational-center root exclusions

Write H=⌊√(N−1)⌋. For a reduced fraction a/b in [0,1], with
b≤H coprime to N, consider |bx−aN|≤H. At a root of unity x
the integer (bx−aN)²−b² is divisible by N and has magnitude <N.
It vanishes, so x≡±1 modulo N. This excludes mixed roots throughout
the strip. The public roots are checked separately.

Even adaptively chosen distinct centers have limited explicit coverage.
For M≥1 centers let B=⌈√M⌉. For each b≤B there are at most b+1
choices of a, each covering at most 2H/b+1 integers. Summing gives
O(HB+B²). Each selected center with b>B adds at most 2H/B+1
positions. Hence total coverage is O(√(NM)+M). Covering a fixed
positive fraction of the second quarter by these individual strips
requires Ω(N) centers. This argument does not address implicit
summation over many centers or a stronger exclusion predicate.

There is also an exact coverage statement when all denominators are
allowed. Put $h=\lfloor\sqrt N\rfloor$. For every candidate $x$,
place the $h+1$ residues $jx$ for $0\le j\le h$ on the circle modulo
$N$. Two consecutive residues have cyclic gap at most $N/(h+1)<\sqrt N$
(duplicates give gap zero). Subtract their indices, reversing sign if
needed. This supplies $1\le b\le h$ and an integer $a$ with
$|bx-aN|<\sqrt N$. For a mixed root, every such denominator shares a
proper factor with $N$: writing $v=bx-aN$ gives $N\mid v^2-b^2$.
Its magnitude is less than $N$, except possibly when $b^2=N$ and
$v=0$. That exception cannot occur for a root since $x$ is a unit
and $0<b<N$. Thus $v=\pm b$, and an invertible $b$ would imply
$x=\pm1$. The unconstrained strips cover every candidate. Finding
which intervals escape the coprime-denominator strips still requires
handling the denominators that share a hidden factor. The converse
that every escaping candidate is a mixed root is not asserted.
Exhaustive checks found a covering denominator at all 22,800 residues
of odd moduli 3 through 301 and verified the nonunit conclusion at
every mixed root (even allowing the closed integer radius there).

The converse fails even on unit residues that pass the elementary gcd
checks. Take primes $p>3$, $q>9p$, $N=pq$, and $x=jq+2$ for
$1\le j<p$. For any coprime denominator $1\le b\le\sqrt N$,
the least distance of $bjq$ from a multiple of $N$ is at least $q$,
because $p$ does not divide $bj$. Hence the least distance of $bx$
is at least $q-2b>\sqrt N$. In contrast $b=p$ covers this candidate
since $2p<\sqrt N$. Every such candidate is a nonroot already modulo
$q$. Excluding the at most three indices where $x=0,1,-1$ modulo
$p$ leaves at least $p-4$ candidates with $\gcd(x(x^2-1),N)=1$.

Choose $9p<q<18p$ (Bertrand's postulate). These survivors occupy
$\Omega(\sqrt N)$ distinct blocks in the canonical partition of width
$2\lfloor\sqrt N/16\rfloor+1$, since they are spaced at least $q$
apart. At most four blocks contain any genuine root. Consequently even
an ideal test of complete coprime-strip coverage leaves power-many
empty blocks uncertified on this infinite family. For $p=101,q=1009$,
exhaustive denominator checks verified all 97 gcd-clean examples.
Locating even one uncovered candidate still gives a possible factoring
route: any short covering denominator for that candidate is nonunit.
The obstruction concerns certifying empty blocks by complete coverage.

**Decoding and density.** On a semiprime N=pq with p<q, set
H=floor(√N) and K=floor(H/p). A short nonunit denominator must be
pj for 1≤j≤K. Its signed remainder is pt with |t|≤K, so jx=t modulo
q. Conversely each such congruence gives a short nonunit denominator.
Since j<q, each pair covers one residue modulo q. Hence at most
K(2K+1)/q ≤ 3/p of all candidates admit even one such denominator.
Uniform independent trials take at least p/3 draws in expectation.
For 9p<q<18p the examples above give an order-matching lower bound
on the success probability.

Given x, run Euclid starting at N,x until the first remainder r≤H,
and take its signed coefficient b of x. Initially b=1 for r=x.
At later steps the consecutive-remainder identity
r_prev |b| + r |b_prev| = N gives |b|<√N since r_prev>H.
This supplies a short denominator in polynomial bit time. Every
candidate escaping all coprime strips therefore yields a proper gcd.
For 101·1009, 1,515 of 101,909 residues admitted some short nonunit
denominator and this Euclidean choice split on 1,072 residues.

For a promised bounded ratio q/p<C and p>√C, a simpler decoder tries
gcd(jx−t,N) for 1≤j≤floor(√C), |t|≤floor(√C). At an escaping
candidate the true pair above has gcd divisible by q. It cannot have
gcd N, since j<p is a unit and jx=t modulo N with |t|≤K<√N
would itself cover the candidate. Thus on this family the hint already
supplies one of a bounded list of linear gcd witnesses. Generating the
hint from N remains open.

**Exact random-input baseline.** The residue classes modulo q above
are precisely the distinct reduced fractions t/j in the signed box
|t|≤K, 1≤j≤K. Two fractions equal modulo q have cross difference
of magnitude at most 2K²≤2q/p<q, so they agree as rationals.
The positive coprime pairs in the K by K square number
2·sum(phi(i), i=1..K)−1 by reflection across its diagonal. Including
negative fractions and zero gives exactly
S=4·sum(phi(i), i=1..K)−1 classes, with phi(1)=1. Uniform sampling
therefore hits the set with probability S/q. This counts the existence
of some short nonunit denominator; a specific decoder can miss it.
The formula matched exhaustive class counts on all 561 pairs of distinct
odd primes below 150. A single gcd(x,N) splits on (p+q−2)/(pq) of
uniform residues, already at least 1/p, so the earlier 3/p bound limits
the per-draw gain of this denominator event to a factor of three.

More generally let x=F(z) modulo N for a uniform residue z and a public
polynomial whose reduction modulo q is nonconstant of degree d. Each
of the S classes has at most d preimages modulo q, so the hit probability
is at most min(1,dS/q). Polynomial degree in the bit length cannot make
this a polynomial-success sampler on balanced semiprimes. This argument
leaves open high-degree succinct maps, distributions on nonuniform
seeds, and reductions that become constant at a hidden prime.

**Succinct powers.** For a uniform unit z and public positive exponent
e, the image of z^e modulo q is uniform on a subgroup of order
m=(q−1)/gcd(q−1,e). If the denominator event has probability at least
δ, the preceding S favorable residues give m≤S/δ. For q/p≤C with
known constant C, use S≤3C, and set B=ceil(3C/δ). Then q−1 divides
E=e·lcm(1,...,B), since m divides the lcm. Given an inverse-polynomial
lower bound δ and a polynomial-bit exponent e, E has polynomial bit
length and is constructible in polynomial bit time. The same argument
covers a polynomial of degree D applied after powering if it remains
nonconstant of degree D modulo q: replace S by DS in the image bound.

This E suffices to split N without the window. Write E=u·2^s with u
odd. For a fresh uniform unit w test gcd(w^(u·2^i)−1,N) for i=0..s.
The q component reaches one by the last step, and its first hitting
index is the base-two logarithm of the order of a uniform element in
the Sylow 2-subgroup modulo q. Each index has probability at most 1/2
since q is odd. The independent p component either never reaches one
or reaches it at a fixed index conditional on its local input. With
probability at least 1/2 the hitting indices differ and some gcd splits.
This is the standard group-order splitting mechanism. Exhausting the
unit bases on 330 choices of small semiprime and exponent multiple of
q−1 confirmed the bound. The argument supplies a decoder for a promised
concentrating exponent, not a way to construct that exponent from N.

**Chebyshev specialization.** Let D_e be the normalized Chebyshev
polynomial D_e(z+1/z)=z^e+z^(−e). For X uniform modulo q, every
input trace X arises twice in the disjoint union of the split cyclic
group of order q−1 and the norm-one group of order q+1 in the quadratic
extension. Indeed the two roots of z²−Xz+1 are inverse, and at X=±2
the single root belongs to both groups. Each output trace has at most
two lifts in each group. Consequently for any S target residues the hit
probability is at most S·(gcd(e,q−1)+gcd(e,q+1))/q. Enumeration
verified the single-target bound for all exponents 1 through 100 and odd
primes below 160 (3,600 prime/exponent pairs).

If the hit probability is at least δ, one of the two groups has power
image order at most 4S/δ. On q/p≤C, the public multiple
E=e·lcm(1,...,ceil(12C/δ)) is therefore divisible by q−1 or q+1.
The previous scalar decoder handles q−1. For q+1 choose a uniform
unit D modulo N and work in the quadratic algebra t²=D. With probability
1/2, D is nonsquare modulo q. Sample a uniform algebra unit z by
rejecting nonunit norms (a proper gcd already splits N). Then
w=z/conjugate(z) is uniform in the q-side norm-one group by cyclicity
of the extension-field units, independently of its p component. Apply
the same successive squarings as in the scalar decoder, testing for
local identity by taking the gcd of N and both coefficients of w^k−1.
Conditional on the nonsquare choice the mismatch probability is at least
1/2 by the same hitting-index argument. Trying both decoders thus costs
expected polynomial bit time under the stated concentration promise and
a polynomial-bit exponent. This connects this specific proposed sampler
to the classical multiplicative and norm-one order methods. It provides
no exponent with the concentration promise.

**Polynomial rules for the degree.** Let p<q be primes in [B,2B],
and allow a shared menu of polylog(B) nonzero integer polynomials F of
uniformly bounded degree and coefficient magnitudes at most a fixed
power of B. The menu can depend on B; selection within it may depend on
N. Use positive exponents |F(N)|, discarding zero. For either sign
ε=±1 we have F(pq)=F(εp) modulo q−ε. If the power image in that group
has size m≤L, then q−ε divides mF(εp). Fix p,F,ε and an integer
L=polylog(B). Unless F(εp)=0, each m≤L allows at most B^o(1)
choices of q, since the nonzero integer |mF(εp)| is polynomially
bounded in B and its divisor count is B^o(1). The exceptional p with
F(εp)=0 number at most the degree per polynomial and sign. Summing
gives at most B^(1+o(1)) covered pairs out of B^(2−o(1)) prime pairs.
Thus even input-dependent selection from this menu misses a fraction
tending to one for the necessary small-image condition. This bounds
the specific concentration mechanism, not the power of arbitrary
algorithms using these degrees.

As a finite check use the menu N−1, N+1, N²−1 and L=8. Enumerating
all p<q in [B,2B) gave 117/210 pairs at B=100, 670/9,045 at B=1,000
and 6,625/533,028 at B=10,000 for which at least one q-side group
image had size at most eight. These are group-size checks with known
factors, not measured factoring success rates.

**Larger degree menus.** The same obstruction extends to a shared
polylogarithmic-size menu of dense nonzero integer polynomials whose
degrees and coefficient bit lengths are polynomial in log B. Choose
a fixed A≥1 such that every nonzero M=|mF(±p)| in the preceding
argument satisfies log(2M)=O((log B)^A), uniformly at p in [B,2B]
and for the allowed polylogarithmic m. Such an A exists by direct
evaluation. Although the full divisor count can now be large, only
divisors at most 3B matter.

For any fixed 0<s<1, their number is at most
(3B)^s times the product of (1−r^(−s))^(−1) over distinct primes
r dividing M. This follows by summing (3B/d)^s over divisors.
There are O((log B)^A) such primes. Bounding the i-th smallest
by i+1, the logarithm of the product is O_s((log B)^(A(1−s))).
Take s=1−1/(2A). The exponential correction is B^o(1).
Summing over signs, menu entries, and m costs only B^o(1).
Exceptional primes p with some F(±p)=0 number at most the sum of
twice the polynomial degrees. Relative to all distinct prime pairs,
the fraction covered is therefore at most B^(−1/(2A)+o(1)).
This uses the same Rankin weighting principle as the earlier
input-independent exponent analysis in quadratic-chirp.md. Here
one freezes p and counts the possible partners after polynomial
specialization. In particular multiplying a menu entry by a shared
polynomial-bit smoothness exponent still falls under this bound.
Constructing an arbitrary polynomial separately from the digits of
each N need not yield a shared menu of this size.

**Same-input trace collisions.** For two distinct public degrees
0≤a<b put u=b−a and v=b+a. Modulo an odd prime r, let h_r be the
fraction of field inputs x with D_a(x)=D_b(x). The exact formula is

    2r h_r = sum over t=r−1,r+1 of [gcd(t,u)+gcd(t,v)−gcd(t,u,v)].

The trace equality holds just when z^u=1 or z^v=1. Inclusion-exclusion
counts this union in each cyclic group, and each trace input has two
lifts across the disjoint union. Direct recurrence evaluation verified
all 3,264 prime/degree-pair cases for odd primes below 100 and distinct
degrees from 0 through 16.

For a uniform input modulo N=pq the two equality events are independent,
so the exact proper-gcd probability is h_p(1−h_q)+h_q(1−h_p).
If it is at least δ, one h_r is at least δ/2 and one of the four gcds
in its positive terms is at least rδ/4. Thus one of the cyclic groups
of order t=r±1 has index t/gcd(t,w)≤8/δ for w either u or v.
The public multiple uv·lcm(1,...,ceil(8/δ)) then covers that group
order. The earlier scalar or norm-one decoder applies without a
balance promise. For b=N, D_N equals D_q as a function modulo p
and D_p modulo q by Frobenius and composition. That crossed degree
identity alone gives no collision probability guarantee. The formula
concerns fixed degrees independent of the sampled trace.

**Local series at trace two.** For odd distinct primes p<q, set N=pq.
As a formal polynomial modulo p, D_N(2+X)=D_q(2+X^p), by
composition and D_p(Y)=Y^p in characteristic p. The linear coefficient
of D_q(2+X) is q² (differentiate the trace recurrence at 2).
Consequently D_N(2+X)−2 has order exactly p modulo p, and order q
modulo q. Over the composite ring its first nonzero coefficient is
therefore at degree p, with gcd q against N. For an explicit example
N=77 first differs from 2 in degree seven: the coefficient is 2
modulo 7 and zero modulo 11.

For any N and 1≤k≤N the exact coefficient c_k satisfies
k c_k=N binom(N+k−1,2k−1), by the Chebyshev recurrence
(equivalently differentiate and expand the second-kind polynomial).
Hence any individual nonconstant coefficient nonzero modulo N
already has an index sharing a factor with N. At 1≤k<N this gives
a proper gcd directly from k. The explicit binomial formula agreed
with the predicted first term on 4,870 coefficient checks across
all distinct odd prime pairs below 80. A succinct test of whether the
series equals 2 modulo (N,X^(T+1)) would decide p>T on semiprimes
and recover p by bisection. Dense truncation reaches order p; the
binary degree of D_N makes ordinary evaluation cheap, but does not
supply this high-order zero test.

**Exact first degree for arbitrary odd N.** For each exact prime
power p^e dividing N set a=floor((e+1)/2), b=ceil((e+1)/2).
The first positive degree surviving modulo p^e is (p^a+p^b)/2.
Take the minimum over the prime-power components for the first degree
modulo N. For example, at N=45 the local degrees are 6 modulo 9
and 5 modulo 5, giving first degree 5 and gcd(5,45)=5.

To prove the local formula substitute X=T²/(1+T). The trace identity
turns D_N(2+X)−2 into ((1+T)^N−1)²/(1+T)^N.
The denominator is a unit and the substitution doubles the first
nonzero degree, also over the ring modulo p^e. From
i·binom(N,i)=N·binom(N−1,i−1), the valuation of binom(N,i) is at
least e−v_p(i). At i=p^j for 0≤j≤e equality holds because the
bottom j base-p digits of N−1 are all p−1 (Lucas' theorem).
A surviving product in the square therefore requires indices i,l
with v_p(i)+v_p(l)≥e+1. The smallest possible sum i+l is p^a+p^b:
at fixed exponent sum, balancing the exponents minimizes the sum of
their powers, with equality only at the balanced pair up to swapping.
At this first sum the minimizing product has valuation e−1. There
are one or two orderings, which cannot cancel for odd p. This proves
the formula. Direct coefficient recurrence verified all 1,828 local
thresholds for every odd N from 3 through 1,999.

For composite odd N the first degree k is below N (also for a single
prime power of exponent at least two), and the preceding coefficient
identity gives 1<gcd(k,N)<N. Binary search using a succinct prefix-zero
test thus splits every odd composite; recurse and handle even factors
directly. Conversely the displayed local formula evaluates the test
in polynomial bit time given the factors, regardless of their number.
This gives polynomial-time Turing equivalence for this particular
scalar on general inputs. The simpler binomial series (1+T)^N−1
already has first nonzero degree modulo N equal to the least prime
factor, by the same valuation identity at i=p. Chebyshev changes the
threshold through squaring and substitution but supplies no known fast
method to compute either threshold from unfactored N.

As a check, using every reduced center with b≤bitlength(N), the
numbers of covered positions in [⌊N/4⌋,⌊N/2⌋) were 7,869 / 366,015
for N=1,464,059, 103,475 / 36,282,880 for N=145,131,521, and
1,258,233 / 3,625,283,755 for N=14,501,135,021. Direct squaring
verified the exclusion on every covered integer.

### A compact rational filter for the second quarter

A uniformly bounded trigonometric polynomial making a constant-size
transition across W residue positions needs degree Ω(N/W) by Bernstein's
derivative inequality on the circle. Even granting a √N-wide guard
band leaves degree Ω(√N) in this bounded, consecutive-frequency model.
The degree requirement changes when division is available.

For odd N put L=bitlength(N), k=⌈log_9(8N)⌉, and define
A(t)=∏_{j=0}^L(t+2^(−j))² and F(t)=A(t)^k/(A(t)^k+A(−t)^k).
The denominator is positive for every real t. For 1/N≤t≤1 choose
a dyadic s with t/2≤s≤t. Its factor in A(−t)/A(t) is at most 1/9;
the rest are at most one. Hence F(t)≥1−9^(−k), and by symmetry
F(−t)≤9^(−k). The rational degree is O(log² N). For each nonzero
residue x let θ=2πx/N. Both |sin θ| and |cos θ| are at least 1/N,
by the grid spacing for odd N and the sine chord bound. Thus
F(sin θ)F(−cos θ) approximates the strict second-quarter indicator
N/4<x<N/2 with pointwise error at most 2·9^(−k). Summing over all
roots of one incurs error at most 1/4, so rounding gives the exact
count. For intervals specified using floored endpoints, individually
check the at most two differing residues by modular squaring.

The expression is small even in bit complexity for known x: on
either sign one product A(±t) has all its unsquared factors of
magnitude at least 2^(−L), giving a denominator lower bound
2^(−2kL(L+1)). Polynomially many precision bits suffice.
Independent direct evaluation and rounding agreed with enumeration
for all 127 odd N from 3 through 255. These checks supplied the
roots by enumeration. The missing computation is summing the
rational weight over the implicit roots. This filter removes the
need for a high-degree polynomial mask; it gives no fast evaluation
of that sum from N alone.

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

### Adaptive modular equality tests on mixed roots

Suppose a public decision procedure on an unknown mixed root uses
addition, subtraction, multiplication modulo odd N, branches on equality
of two registers modulo N, and Boolean leaves. Public constants and the
next operation may depend on earlier branch answers. If the procedure
gives different answers on two mixed roots, it yields a deterministic
factor finder with constant many modular operations and gcds per step.

Simulate one path, maintaining two residues per register by substituting
+1 and −1 for the input along that same path. These are formal
substitutions even when neither public input would follow the simulated
branches. At an equality gate subtract its two registers and take the gcd
of each substituted difference with N. A proper gcd splits N. Otherwise
each difference is either zero or a unit. Every mixed root uses the plus
substitution on at least one prime-power component and the minus
substitution on at least one other. Thus equality holds at every mixed
root exactly when both differences are zero, and fails at every mixed
root otherwise. Continue along this common branch. If simulation reaches
a Boolean leaf, all mixed roots share its answer, contradicting the
distinguishing promise. The reasoning allows repeated odd prime powers.

In particular the lower-half predicate separates a mixed root from its
negative. A polynomial-bit-time modular equality decision tree for that
predicate, promised correct only on mixed roots, would already factor
the modulus without receiving a root. The restriction to modular
equality matters: comparison of canonical integer representatives reads
the lower-half bit immediately on a supplied root. This argument does
not classify algorithms that evaluate an aggregate directly.

### Parity of a modular equality test summed over root tuples

Consider the XOR of a Boolean leaf over all ordered tuples of m roots
of one, using the same arithmetic and branching model as above. Write
k for the number of distinct prime factors of odd N. Each prime-power
component chooses a sign pattern in {+1,−1}^m for the whole tuple.
Group the tuples by the set S of patterns used across components.

Simulate a group by running one formal substitution per member of S
along a common path. At each equality branch compute the gcd of every
substituted difference with N. A proper gcd gives a factor. Otherwise
each difference is either zero or a unit, so the common answer is true
exactly when all substituted differences vanish. The leaf is then
constant across the group. Its size is the number of surjections from
k labeled components onto S: permuting the pattern labels acts freely,
so this size is divisible by |S| factorial. Every group with at least
two patterns therefore contributes even parity. The remaining XOR is
over the public tuples of signs, one tuple per singleton group. The same
cancellation holds on the restricted domain of tuples having at least one
mixed coordinate, where all singleton groups are absent.

For two distinct prime factors, only singleton and pair groups occur.
Simulating every such group costs O(4^m) times the per-path cost and
requires no factorization. This gives an effective dichotomy for
m=O(log log N): either the simulation exposes a factor or the full
parity equals the XOR on public sign tuples. For arbitrary k the
cancellation proof persists but enumerating all possible pattern sets
can be expensive. This result is specific to arithmetic with modular
equality branches and Boolean leaves. Ordered interval comparisons,
weights computed outside the ring, and other aggregate representations
require their own analysis.

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

### Full-product evaluation given an arbitrary factorization

Write P_N(u)=∏_{x=0}^{N−1}(x²−u). For each exact odd prime power
p^e dividing N, compute the least residue b of −u(u^((p−1)/2)−1)²
modulo p. Then P_N(u) ≡ b^(N/p) modulo p^e.

To prove the lift, let A_j be the product on one period of length p^j.
For j≥1 group the next period by x=a+t p^j with 0≤t<p.
Modulo p^(j+1) the quadratic term in t drops out of each x²−u.
In the product over t, terms with at least two linear corrections
also vanish. The one-correction terms sum to a multiple of
p^j ∑t, divisible by p^(j+1) for odd p. Thus A_(j+1) ≡ A_j^p
modulo p^(j+1), including when constants are nonunits. Raising
congruent integers modulo p^j to the pth power preserves congruence
modulo p^(j+1) by the binomial theorem. Induction from the one-period
formula A_1 ≡ b modulo p gives A_e ≡ b^(p^(e−1)) modulo p^e.
There are N/p^e full periods. For 2^e exactly dividing N, half the
factors are even and N/2 ≥ 2^(e−1) ≥ e, so the local residue is zero.

CRT and modular exponentiation now evaluate the full product at any
point in polynomial bit complexity given the prime factorization.
Together with recursive random-unit splitting above, this establishes
expected polynomial-time oracle equivalence for every N≥2.

### Normalizing the product gives the local quadratic characters

For odd N and a unit u, put C=(−4u)^N modulo N. At every exact
prime power p^e dividing N, P_N(u) is zero if u is a square modulo p
and equals C otherwise. Indeed N−N/p is divisible by φ(p^e), and
replacing the base residue b by −4u in the nonsquare case is valid
after raising to N/p by the lifting congruence above. Consequently
E=P_N(u)/C is the CRT idempotent with value one on the nonsquare
components. The root of unity S=1−2E has local values (u/p) lifted
as the integers ±1 modulo p^e.

Can a single exponent k give S=u^k for every unit? Modulo p^e the
unit group is cyclic of order m=φ(p^e). Testing a generator forces
k ≡ m/2 modulo m, and this condition suffices. These congruences
are compatible exactly when all v₂(p−1) agree: necessity follows
from v₂(k)=v₂(m)−1. For sufficiency take k=lcm(m)/2, since the lcm
divided by each m is then odd. This exponent uses the factorization.
For N=15 the requirements already conflict modulo 2 and 4. Thus
even an exponent chosen separately for each N cannot implement the
normalized sign on every modulus and every unit.

### The exponents of a sparse sign polynomial suffice to factor

Suppose N is odd and F(X)=∑ c_i X^{d_i} realizes the normalized sign
S on every unit modulo N. For each prime p dividing N, every nonzero
residue modulo p lifts to a unit modulo N. Thus F agrees on F_p^×
with X^((p−1)/2). Reduce the exponents modulo p−1 and collect
coefficients. Both resulting polynomials have degree at most p−2
and agree on p−1 points, so their coefficients agree. In particular
some positive d_i equals (p−1)/2 modulo p−1.

Hence M=N·lcm{2d_i : d_i>0} is a multiple of every φ(p^e) for
p^e exactly dividing N. Its bit length is bounded by the input bit
length plus the sum of the listed exponent bit lengths (up to one bit
per exponent). Even a superset of the support works. The coefficients
are unnecessary.

To split an odd divisor D of N with at least two distinct prime factors,
write M=2^s t with t odd. For a uniform random unit a modulo D, take
gcd(D,a^(2^j t)−1) for j=0 through s, using repeated squaring.
The local orders of a^t are independent powers of two. On each odd
prime-power component every individual order has probability at most
1/2 by cyclicity of the unit group. With probability at least 1/2 two
chosen components have different orders, and some gcd separates them.
A sampled nonunit can already give a factor. Recurse with the same M,
handling primes and prime powers by primality and exact root tests.
Thus a polynomial-bit explicit exponent list covering a correct sign
polynomial is itself sufficient for expected polynomial-time factoring.
This argument requires an explicit list and gives no bound for circuits
whose expanded support is large.

The exponent test also covers a sparse rational expression A/B for S
when B is invertible on all units. More generally, assume A(u)=S(u)B(u)
on all units and that B is not identically zero as a function on
F_p^× for any prime p dividing N. Modulo X^(p−1)−1 over F_p, the
reduced numerator equals the reduced denominator shifted cyclically by
(p−1)/2. Choose any nonzero reduced denominator coefficient. Its shifted
numerator coefficient is nonzero as well, so some original numerator
exponent a and denominator exponent b satisfy a−b ≡ (p−1)/2 mod p−1.
Thus the list of nonzero absolute differences |a−b| passes the preceding
exponent test. Its length is at most the product of the support sizes,
and it uses neither the coefficients nor an expansion of 1/B.

The same support extraction tolerates errors. Fix A and B, assume B is
a unit on every unit input, and let T be their total number of terms.
If A/B returns S modulo N with error strictly below 1/T on a uniform
unit, the preceding exponent list still factors N. For each p dividing N
the local discrepancy D=A−X^((p−1)/2)B on F_p^× has error probability
at most the global one, because projection of a uniform unit is uniform.
Collect exponents modulo p−1. If D has t>0 surviving terms, then in the
cyclic sequence D(g^j) for a generator g no t consecutive values vanish:
such values would give a homogeneous Vandermonde system on t distinct
nodes g^d with nonzero coefficients multiplied by units. Each cyclic
window of length t therefore contains an error. Counting windows gives
error fraction at least 1/t ≥ 1/T, a contradiction. Thus the local
function identity holds at every prime, enough for exponent extraction
even if some values still disagree at prime-power precision.

For N=pq with distinct odd primes, the same proof works when uniform
sampling is restricted to Jacobi symbol −1. Every nonzero residue at
p then has exactly (q−1)/2 allowed partners at q, and conversely, so
both projections remain uniform. The error bound concerns the actual
distribution for one fixed formula. A finite run of successful tests
alone does not certify it.

### Weak predictions on Jacobi-negative inputs

For a fixed distinct odd semiprime N, a random unit u of Jacobi symbol
−1 always gives a mixed root S(u). If a trial outputs either S(u) or
−S(u) with probability ε, gcd(output−1,N) finds a proper factor with
probability at least ε. Independent trials therefore suffice in expected
1/ε attempts. Every gcd is verifiable even when the other outputs are
arbitrary. No sparsity hypothesis enters this test.

There is also a one-bit version. The two possible mixed roots are a
single opposite pair modulo N. Exactly one representative in [0,N)
is even; call it r. This r depends only on N. Encode it in n=bitlength(N)
bits. Suppose a predictor given N, a fresh uniform Jacobi-negative unit
u, and an independent uniform n-bit mask m guesses the binary inner
product of m with r modulo two with probability at least 1/2+δ.
Regard the predictor as a randomized oracle on m, resampling u and
its internal coins each call. Its expected signed output then has
Walsh coefficient at r at least 2δ. The Goldreich–Levin list decoder
recovers a list containing r with constant probability in time polynomial
in n and 1/δ (use independent oracle samples for its correlation
estimates). Test list entries by gcd with N after subtracting one.
This suffices at inverse-polynomial advantage for each N, even though
predictions of individual fixed bits need not have that advantage.
The assumption concerns random *masks* of the even representative;
an ordinary bit guess with a skewed marginal does not imply it.

### The low bit selects one hidden prime character

For N=pq with distinct odd primes, let s be the representative of S(u)
in [0,N). Encode odd as +1 and even as −1. On all units, this signed
low bit equals (u/ℓ) for one fixed ℓ in {p,q}. Indeed, equal local
signs give 1 or N−1, with the required parities. Of the two mixed
roots exactly one is odd. Take ℓ to be the prime where that odd mixed
root equals +1. The two mixed sign patterns then give the same rule.
Multiplication by the public Jacobi symbol gives the other local
character. In particular the Jacobi symbol alone predicts the low bit
with exactly half accuracy on uniform units, even conditional on either
Jacobi sign, by independence of the two local signs.

For N=77 the odd mixed root is 43, equal to +1 modulo 7, so this bit
is the Legendre character modulo 7. An exact low-bit oracle therefore
gives a quadratic-residuosity test on units (both recovered signs must
be positive). This observation supplies no classical polynomial-time
procedure for extracting the hidden prime from character values. The
random-mask decoding argument above requires queries on many masks
of a fixed root, which this single-bit oracle does not provide.

For a public unit scale a, the signed low bit of [a S(u)]_N (least
nonnegative residue), divided by the signed low bit of [a]_N, again
equals one of the two local Legendre characters. The signs at the
two equal-sign inputs are +1 and −1 after normalization; the two mixed
roots receive opposite signs because their scaled representatives add
to odd N. Thus the selected prime can depend on the scale. Changing u
in the unscaled character never changes which prime it selects.

With a scaled-bit evaluator, fix any Jacobi-negative u and set r=[S(u)]_N.
For j≥1, the ordinary parity bit of [2^j r]_N equals
floor(2^j r/N) modulo two, since 2^j r is even and N is odd.
These are consecutive binary digits of r/N. After
ceil(bitlength(N)/2) nonadaptive queries the compatible interval has
width at most √N, so the earlier shifted quadratic solver recovers the
mixed root among its candidates. Equivalently, tracking whether the
selected character changes relative to scale one suffices: at a
Jacobi-negative u a switch flips the normalized sign. Guessing the
single initial orientation costs only two candidate prefixes. This
reuses the existing binary-prefix decoder; the additional assumption
is access to the scale-dependent bit or character selection.

### Exact shifted low-bit correlation returns the complementary factor

Let f be the signed low bit on units, extended by zero elsewhere. For
N=pq write p for the prime selected by f. For every unit shift h,
the integer sum of f(x)f(x+h) over x modulo N is 2−q. CRT gives
q−2 choices avoiding 0 and −h at q. At p the sum of the Legendre
symbol of x(x+h) is −1: away from x=−h it equals the character of
x/(x+h), which runs once over all residues except 1. The omitted
summand is zero. Thus one exact sum at h=1 gives q by subtraction.

Conditional on x and x+h both being units, a uniform sample of the
sign product has mean −1/(p−2). On balanced semiprimes, distinguishing
this Bernoulli sign distribution from a fair sign from independent
products alone needs Ω(p²) samples for constant advantage: its
single-sample relative entropy from fair is O(1/p²), and Pinsker's
inequality bounds the product-distribution advantage. Every unit shift
has this same law. Shifts with a proper nonunit gcd already split N.
This only bounds using the sampled sign products as independent bits;
the positions, correlated queries, and structured summation remain
available. The exact identity gives no fast way to sum the primitive.

### Counting consecutive units already contains the factors

For distinct odd primes N=pq the number C(2) of x modulo N with x
and x+1 both units equals (p−2)(q−2). Thus p+q=(N+4−C(2))/2.
In fact C(2)=2φ(N)−N+2, so the unsigned support count here
is affinely equivalent to the totient. Membership is public by gcd.

More generally C(k) counts runs of k units, the Schemmel totient of
order k. If k does not exceed the smallest prime divisor, CRT gives
C(k)=N∏_{p|N}(1−k/p), since k distinct residue classes are forbidden
at each p. This identity yields a deterministic polynomial-time
reduction from arbitrary integer factoring to exact run counts of
logarithmic length. Trial-divide out all primes ≤ L=bitlength(N);
let M be the remaining cofactor. When M>1 every remaining prime
exceeds L, and there are at most L distinct primes. The counts on M
at k=0,...,L (with C(0)=M) interpolate the integer polynomial
(M/rad(M))∏_{p|M}(p−k). Its integer roots are precisely the distinct
remaining primes. Polynomial factorization over the rationals and
exact division recover them and their multiplicities in polynomial
bit time. The interpolant has degree at most L and polynomial-bit
coefficients, for example by bounding its coefficient norm by
(M/rad(M))∏(p+1) ≤ M·2^L. This identifies the cost as counting
points exactly even for a public predicate; it gives no shortcut
for computing a Schemmel totient without the factors.

For odd N>1 let p be its least prime divisor. The parity at even
lengths is precisely C(2t) mod 2 = [2t < p]. Each valid start x
pairs with 1−2t−x (negate the run and reverse its order). The unique
fixed start satisfies 2x=1−2t modulo N. Doubling its entries gives
the odd numbers from 1−2t through 2t−1, all units exactly when
2t−1<p. Since p is odd this is equivalent to 2t<p. For odd lengths
the fixed run contains zero and the parity always vanishes.
Consequently the first even length with even count is p+1. Binary
search on even lengths recovers p with O(log N) parity queries,
including for repeated prime factors. Below the first prime all even
lengths give the same bit, so restricting lengths to polynomial in
the bit length loses the threshold for balanced semiprimes.
The fixed-point predicate is exactly coprimality of N with the odd
double factorial (2t−1)!!, or equivalently with (2t)! for odd N.
This brings the symmetry approach back to the factorial gcd test.

### Wide intervals via an extracted prime power

For the target y=1 on odd N, trial-divide by the primes up to bitlength(N)
and select the largest exact power s of one of these primes dividing N.
Write m=N/s. Each root is in one of the two disjoint progressions ±1 mod s.
In a half-open interval let f be the first point of such a progression.
For x=f+sk the remaining condition modulo m is the monic quadratic
(k+f/s)²−1/s²=0, with division by s performed modulo m. When m=1
only the progression length is needed. This calculation handles an
unfactored composite m and counts all matches, including when the two
prime factors of m select opposite signs.

Apply this to the same constant number of intervals selected by
reflection. If their total length is at most s⌊√m⌋, each progression
has at most ⌊√m⌋ entries per interval. A constant number of translated
short-root calls suffices, with degrees and lattice ranks polynomial in
the bit length under the existing solver bound. Trial division through
the bit length is likewise polynomial. The old cutoff was ⌊√N⌋.
When s and m are comparable the new threshold is Θ(N^(3/4)).
The improvement relies on a publicly extractable divisor and does not
apply to balanced semiprimes with two large factors.

At N=3^8·101·103=68,254,083, intervals of width 524,288 starting at
0, 524,288, and 6,291,456 gave parities 1, 0, and 1 respectively.
All three exceeded the old reflection cutoff. At width 32,768 the new
path resolved all 2,082 previously unresolved grid starts j·32,768<N,
including seven odd answers. The same width resolved 1,050 such starts
at N=5^5·101·109, also with seven odd answers. Reference values used
all roots computed independently by SymPy, while the evaluator used
only N, the interval, the target one and the existing half-parity bit.

### A common Jacobi sieve for all candidate factor sums

Let b be the bit length of odd N. Sample m uniformly among the odd
integers in [1,8N]. Reject a candidate sum s when the Jacobi symbol
((s²−4N)/m) is −1. Every nonnegative square discriminant passes,
including when the symbol is zero. For a fixed nonsquare integer
1≤D≤N, quadratic reciprocity gives a nonprincipal real character on
the units modulo Q=4D, taking m to (D/m). Indeed, after removing
square factors its primitive character has the fundamental discriminant
of the nontrivial squarefree part of D as conductor, dividing 4D.
Exactly half the units of a period give −1. Complete periods cover
at least half of [1,8N], and all these witnesses are odd. Hence one
test rejects D with probability at least φ(Q)/(2Q). If Q has h
distinct prime divisors, its i-th smallest prime is at least i+1;
multiplying the totient factors gives φ(Q)/Q≥1/(h+1)≥1/(b+3).

With T=4(b+3)b independent choices the probability any particular
nonsquare survives is at most exp(−2b). A union bound over *every*
D≤N gives failure probability ≤N exp(−2b)<1/N. Thus the same short
list of unfactored public moduli works simultaneously throughout
the balanced factor-sum range (where 0≤s²−4N≤N/2). This is an
unconditional randomized statement. Computing each Jacobi symbol
costs polynomial bit time. In the 2,989,601,243 benchmark, scanning
the 3,316 even sums with m=2·randrange(4N)+1 from Python Random
seed 20261127 left 1,636 after one test, 197 after four, 19 after
eight and just 112,716 after fifteen.

This gives a local-property description with a short public list.
A proposed sum already admits a deterministic integer square test.
The outstanding search cost is locating a sum satisfying the tests
without scanning the interval. For large composite m the acceptable
classes are implicit; the small-prime wheel cannot list them in
polynomial space just from the short description.

### Contiguous blocks rejected by a single small modulus

For odd N, index the even candidate sums by s=2t. For every odd
modulus m>1, the class t=(N+1)/2 modulo m always passes the
square-discriminant test, since (N+1)²−4N=(N−1)². It also passes
the Jacobi filter, which accepts zero. Thus every m consecutive
indices include a survivor of that particular test. A certificate
that covers W consecutive indices by contiguous blocks, each rejected
throughout by one odd modulus at most L, needs at least
ceil(W/(L−1)) blocks, even with overlaps and adaptive modulus choices.
For a balanced search of width proportional to √N, taking L polynomial
in log N leaves exponentially many blocks in the input bit length.
This bound concerns that certificate format. A joint congruence
argument, interleaved residue classes, or larger moduli require
separate analysis.

For the benchmark, 112716²−4·2989601243=27322², giving the verified
split 42697·70019=2989601243.

### Certifying an empty interval from rejection moments

Use the first fifteen Jacobi moduli generated above and s=2t with
54678≤t≤56335 (the lower half). Let r(t) be the number of tests
rejecting t. Direct enumeration gives the following frequencies for
r from 0 through 15 in order:

    0 2 0 17 77 141 278 322 318 250 159 68 18 7 1 0

The ordinary even-order Bonferroni upper bound on the survivor count
is the true count plus the sum of binom(r−1,d) over candidates with
r>d, by the alternating partial binomial identity. Orders 8, 12 and
14 give 12463, 20 and 0 respectively.

A different certificate uses only moments through order eight on this
instance. Any real polynomial P with P(0)=1 gives the survivor upper
bound Σ_t P(r(t))². Put G_ij=Σ_t r(t)^(i+j) for 0≤i,j≤4, with
zeroth powers equal to one. G is positive definite because the observed
support has more than four distinct values. Constrained quadratic
minimization gives the optimal coefficients G⁻¹e₀/(G⁻¹)₀₀ and the
exact upper bound 1/(G⁻¹)₀₀ = 279106006052799305/333564824161743793, below one.
Thus these nine integer moments imply emptiness without requiring a
table of individual passes. At degree three the optimum is above 2.51.

The moments here came from scanning all 1658 candidates. An arithmetic
evaluator would need to obtain adequate moments directly from N and
the endpoints. Each moment of fixed order k expands into intersection
counts involving at most k tests. This example supplies neither a fast
intersection counter nor a uniform bound on the degree needed.

### Computing fixed-order moments from small prime filters

The [CRT moment calculation](check_factor_sum_moments.py) uses the ten
odd primes through 31 on the same lower half. The kth binomial moment
is the sum, over k-test subsets, of their simultaneous rejection counts.
Split a subset into coprime products A and B and form the two lists of
rejecting residues. For a left residue a write t=a+Aj. The allowed
values of j+a/A modulo B are precisely the right residues multiplied by
the inverse of A modulo B. Two binary searches in that sorted list,
with whole periods accounted for, count hits between the endpoints.
Here a/A denotes modular division. The empty side has modulus one.

With T primes at most L, moments through fixed order k cost at most
O((T+1)^k L^ceil(k/2) poly(log N, k, log L)) bit operations for
endpoints of O(log N) bits, independent of interval width except for
endpoint arithmetic. Thus fixed-order moments are polynomial-time
computable when T and L are polynomial in the input bit length.

The ten-prime benchmark gave a squared-polynomial survivor bound of
1.54719 at degree four and 0.55330 at degree five (rounding to five
decimal places), certifying emptiness at moment order ten. Direct
enumeration verified every binomial moment. The largest half-list had
31,104 residues for an interval of only 1,658 candidates. Here order
ten uses every subset of the ten tests, so the certificate gives no
degree saving over exact inclusion–exclusion. A general algorithm
would need a choice of small primes that separates the target and a
uniformly sufficient, efficiently computable moment order.

### Exact moment threshold on the ten-prime benchmark

For the lower half in the preceding CRT experiment, the aggregate
integer moments through order nine suffice to prove emptiness, and
order eight leaves an integral ambiguity. For r in {0,…,10} use the
polynomial F(r)=∏_{j=2}^{10}(1−r/j). It equals one at zero, 1/10 at
one, and zero elsewhere on this domain. Its sum therefore bounds the
survivor count from above. Equivalently compute the sum as
Σ_{j=0}^9 (−1)^j (10−j) M_j/10 where M_j=Σ_t binom(r(t),j).
The coefficients are the forward differences of F at zero.
The measured moments give exactly 7/10.

For the converse, the actual rejection histogram and an alternative
nonnegative integer histogram indexed from zero through ten are

    0 7 40 128 293 370 397 267 124 30 2
    1 0 58 116 251 496 229 399 61 47 0

Their generating polynomials differ by (1−z)^9(1+2z).
The first eight derivatives and the values agree at z=1, so the two
histograms give identical binomial moments M_0 through M_8. Any
inference using only these aggregate integers, including parity
and congruences, must allow both zero survivors and one. The second
histogram is a witness to moment ambiguity. Realizability under the
same arithmetic filters is an extra constraint outside this argument.

### A degree lower bound for small-prime interval certificates

Let X_i(t) be the rejection bit of a filter modulo a distinct prime
p_i≤L. Discard constant-pass filters and assume each remaining filter
both passes and rejects some residue. This holds for the nonconstant
factor-sum tests since s=N+1 always passes locally. Write T for the
number of filters. Consider any real polynomial F of degree at most k
in their bits, nonnegative on the Boolean cube and satisfying F(0)≥1.
Its sum over W consecutive t upper-bounds the survivor count. With
D equal to the sum of binom(T,j) for 0≤j≤min(k,T), this upper bound
is itself at least

$$
\max\{0,(W-D L^{2k})/(D L^k)\}.
$$

To see this, take expectation μ over a full CRT period. The bits are
independent Bernoulli variables with means α_i in [1/L,1−1/L].
Their centered, variance-one versions have absolute value at most √L.
Products φ_S of these centered variables form an orthonormal basis.
Multilinearizing F leaves only |S|≤k. Each coefficient
c_S=E[F φ_S] satisfies |c_S|≤μ(F)L^(k/2) by nonnegativity.
Evaluating at zero gives 1≤F(0)≤D μ(F)L^k.
Each nonconstant basis function has mean zero over its own period
of length at most L^k. Its sum on the interval therefore has magnitude
at most L^(3k/2), by discarding whole periods and bounding the tail.
Consequently the interval sum of F is at least (W−D L^(2k))μ(F),
giving the displayed inequality whenever its numerator is positive.

In particular, a certificate strictly below one needs
W < D(L^(2k)+L^k). For T,L polynomial in b=log N and
W=N^δ with fixed δ>0, D≤(T+1)^k for k≥1, hence k=Ω(b/log b).
This covers squared-polynomial moment bounds and even arbitrary
nonnegative low-degree polynomials using the labeled bits separately.
It allows the filters and coefficients to depend on the interval.
The result concerns a single such majorant summed over the interval.
It leaves open other deductions using arithmetic facts about which
bit patterns occur, as well as efficient evaluation at growing order.

### A least-nonresidue problem inside uniform small-prime separation

Take a prime r≥17 congruent to 1 modulo 8. Set
N=(4r−1)(4r−9), t=4r−3, and s=2t. The two displayed odd factors
have ratio less than two. The true sum of these factors is 8r−10.
The candidate s=8r−6 is four larger, and lies in the balanced
interval since t²−N=16r and 128r≤N. It is false for every factor
pair because its discriminant is 64r, a nonsquare. Also gcd(t,N)=1
since gcd(t,16r)=1.

An odd prime p rejects this s in the square-discriminant sieve
exactly when the Legendre symbol (r/p) is −1. At p=r it is zero.
Quadratic reciprocity gives (r/p)=(p/r) for distinct odd primes,
and 2 is a quadratic residue modulo r. Thus the first rejecting
prime is precisely the least positive quadratic nonresidue modulo r
(the least nonresidue must be prime by multiplicativity).
As log N is of order log r, a uniform polylogarithmic prime cutoff
separating every false sum on every balanced odd composite would
give a polylogarithmic least-nonresidue bound for every such r.

This is an implication for complete pointwise separation by the
odd-prime discriminant filters. The constructed factors need not
be prime. The implication does not require a factoring algorithm
to reject this candidate, and a direct integer square test rejects
it immediately. For r=1009, N=16248945 and s=8066, the first
rejecting prime is 11.

### The size of an explicit survivor list on general inputs

For arbitrarily large even k, take X=2^k and partition [X,2X) into k
equal subintervals. By the prime number theorem and averaging, one
subinterval contains at least k primes: the mean is of order X/k².
Choose k distinct primes there and take their product N, an odd
squarefree integer of bit length b=Θ(k²). Its smallest prime factor
is at least 2^k, eventually exceeding any fixed polynomial in b.

The ratio of the largest to the smallest chosen prime is at most
1+1/k. For any split into two sets of k/2 primes, the ratio of the
larger product to the smaller is at most (1+1/k)^(k/2)<2. All these
pairs are balanced. Unique factorization makes their unordered pairs
distinct. A sum and the common product determine the unordered pair,
so there are at least binom(k,k/2)/2 distinct genuine balanced sums.
This grows as exp(Ω(√b)). Each discriminant is a perfect square and
every sound square-discriminant filter passes every one.

For k=20, the first twenty primes above 2²⁰ run from 1048583 to
1048837. Their product has 401 bits. The product
of the ten largest is less than twice the product of the ten smallest,
verifying at least 92,378 balanced sums directly.

Thus a uniform algorithm explicitly listing all surviving balanced
sums can exceed polynomial output length before false positives enter
the analysis. For general inputs an interval evaluator may return one
verified witness or guide a split without listing all witnesses. A
polynomial list bound remains meaningful under a semiprime promise.

### Two odd inputs suffice to balance an arbitrary split

Suppose an evaluator decides whether an interval contains the sum of
two positive factors a≤b≤2a of a given odd integer. Remove factors
2 and 3 and test primality first. For a remaining composite N, let
3^k be the largest power of three at most N/3. Query the full balanced
sum interval on N·3^k and N·3^(k−1). These are odd inputs below N².
At least one query is occupied, and every balanced pair on either
input exposes a proper divisor of N.

For existence, take any nontrivial split N=uv with 5≤u≤v.
Then v/u≤N/25<3^k. Some integer j in [0,k] has 3^j within a
multiplicative factor √3 of v/u, by rounding the logarithm.
Choose h from {k−1,k} with the parity of j. Then h≥j and the factors

    u · 3^((h+j)/2)   and   v · 3^((h−j)/2)

have product N·3^h and ratio (in either order) at most √3.
For soundness, if either factor of a balanced pair were coprime to N,
it would divide the padding 3^h≤N/3, and its companion would be at
least N. Their ratio would be at least three. If one factor contained
all of N, its companion would divide the padding, giving the same
contradiction. Thus taking gcd with either factor splits N.

On an occupied input M, bisect the integer interval from ceil(2√M)
through floor(3√(M/2)), then recover the factors from the sum using
its square discriminant. All endpoints can be computed with integer
square roots. This costs O(log N) interval queries per split and
recurses on proper divisors. The reduction requires a correct empty
answer on unoccupied queries and supports repeated prime factors in
the padded inputs. Exhaustive enumeration checked existence and both
gcd conclusions for all 571 composites coprime to six through 3000.

### Balancing with one squarefree odd multiplier

For odd composite N≥3465 let K be the largest initial product of odd
primes 3·5·7··· at most N/3. Let p be the next prime. Trial divide N
through 6p, returning a proper divisor if found. Smaller inputs can be
handled by trial division. This preprocessing is polynomial in log N:
the index of p is O(log N), and the nth prime is O(n log(2n)).
Now gcd(N,K)=1 and K>N/(3p).

There is a chain of divisors of K covering [3,K/3] with consecutive
ratios at most 5/3. For K=1155 the chain is
3,5,7,11,15,21,33,35,55,77,105,165,231,385.
On appending a prime q, merge the chain and its q-multiple. Their
ranges overlap if q≤K/9, and the merged chain retains the ratio
bound. This inequality holds initially for q=13. It persists since
the following prime is less than 2q by Bertrand's postulate while
the new primorial is Kq.

For any nontrivial split N=uv with u≤v, the ideal divisor for padding
u is x=√(Kv/u). We have x≥√K>3 and x≤K/3 because
u>6p>3√(3p)>3√(N/K). Choose the chain divisor d nearest to x on a
logarithmic scale. Then ud and vK/d differ by a factor at most 5/3:
their ratio is (d/x)² and d is within a factor √(5/3) of x.
The unknown split and d are used only to prove that the window is
occupied. The algorithm constructs K from N alone.

Query balanced interval occupancy on the single odd input NK and
bisect to recover a sum. Every returned balanced pair splits N by
gcd as in the preceding reduction because K≤N/3 and gcd(N,K)=1.
The input NK is squarefree whenever N is squarefree. For
N=1009·1000003 the multiplier is 111546435 and a verified balanced
pair is (326232907,345001035), with gcds 1009 and 1000003.
The padding introduces many known prime factors. An oracle promised
only semiprime inputs is outside this reduction.

### Fermat search after primorial balancing

Consider distinct primes p,q in [B,2B], N=pq. Fix any menu of
poly(log B) odd multipliers K≤B^C with C fixed, allowing the choice
within the menu to depend arbitrarily on N. For 1≤H≤B, the fraction
of these semiprimes admitting a split via the first H ordinary Fermat
midpoints on some KN is at most B^o(1) √(H/B). In particular it
tends to zero for H≤B^(1−ε) with fixed ε>0. All primorials used by
the preceding reduction on this input range form such a menu.

Indeed, for each K discard inputs sharing a prime with K, counting
them all as successes. K has O(1) prime factors in [B,2B], so their
combined fraction over the menu is B^(−1+o(1)). For coprime inputs,
a pair exposing a factor by gcd must be a=pd, c=qK/d with d dividing
K, up to exchange. Set M=KN and u=log(a/c)/2. Its Fermat midpoint is
(a+c)/2=√M cosh(u). Reaching it in H tests starting at ceil(√M)
requires √M (cosh(u)−1)≤H. Since cosh(u)−1≥u²/2 and √M≥B,
this forces |log(q/p)−log(d²/K)|≤√(8H/B).

For fixed p,K,d the allowed q in [B,2B] therefore occupy an interval
of length O(√(BH)), by the mean value theorem for the exponential
on log(q/p) in [−log 2,log 2]. Counting all integers there gives
O(1+√(BH)) possibilities. The classical divisor bound gives
τ(K)=B^o(1) uniformly for K≤B^C. Sum over divisors, multipliers and
the O(B/log B) primes p, and divide by the Θ((B/log B)²) prime pairs.
This proves the bound, including rounding and the discarded inputs.

This result is specific to consecutive Fermat midpoints and to a
small shared menu. It permits choosing a different menu at each scale.
A method that selects among many more multipliers from the arithmetic
of N, or skips large sets of midpoints, needs separate analysis.
The balancing reduction itself uses interval occupancy and still holds.

### Certifying a short interval of Fermat midpoints

For odd M and integer midpoints A≤t≤B with A≥ceil(√M), set
g(t)=t+√(t²−M). The possible larger factors lie in the exact integer
interval [ceil(g(A)), floor(g(B))]. Conversely every divisor in that
interval gives an integer midpoint in [A,B], since M is odd and
g is increasing. Intersect with the balanced factor range if required.

This permits the standard degree-one unknown-divisor lattice already
used in [the affine decoder](../../src/trackE/coppersmith.py).
For exact occupancy, retain every integer root of a certified short
polynomial and test whether the corresponding candidate divides M and
lies inside the requested interval. Returning a gcd from an outlying
candidate could give an incorrect interval answer.

For M=2989601243 we tested [56352,56360] (larger factor 70019)
and [56340,56348] (empty), using one five-dimensional lattice per
interval. With c the floor of the center of the integer factor
interval and Y its covering radius, scale coefficient i by Y^i in
the basis M², M(c+x), (c+x)², x(c+x)², x²(c+x)². The first LLL
row in each case had sum of absolute entries below the square of the
lower factor endpoint. Its unscaled polynomial must vanish at every
valid offset: the value is divisible by the square of its factor
and has strictly smaller absolute value. Exact integer-root extraction
then gives both answers.

A fixed positive exponent slack ε in the standard small-root theorem
permits factor intervals of width O(M^(1/4−ε)) in polynomial bit time
for 0<ε<1/4, since the larger factor is at least √M. For a fixed
factor ratio strictly above one, g has bounded derivative, so the
same scale of midpoint widths suffices there. Near √M the derivative
diverges; the exact transformed width controls that case. Covering
the entire balanced factor interval [√M,√(2M)] at this scale still
takes order M^(1/4+ε) blocks. This provides a rigorous local stopping
condition for an interval search.

### Avoiding precision loss when decoding low sum bits

Suppose odd N>30 has a balanced pair a≤b≤2a. Handle a common
factor with 15 directly. Set M=N when N≡3 mod 4, and M=15N
otherwise. In the second case 5a and 3b give a balanced pair:
their ratio 5a/(3b) lies in [5/6,5/3]. Every balanced pair on M
splits N by gcd since a factor coprime to N would divide 15,
leaving a companion at least N>30. The multiplier also preserves
squarefreeness when applied after the coprimality check.

Now M≡3 mod 4, so the two factors of any pair have opposite residues
modulo four. Their sum is zero modulo four. Given its residue s modulo
2^k for k≥3, use the even half h=s/2 of the canonical representative.
The larger factor satisfies (x−h)²≡h²−M modulo 2^k, and the right
side is odd. For a genuine hint it has exactly four roots, each a
single residue at full stride 2^k (lift by testing both next bits).
Thus this preprocessing eliminates the singular cases in the existing
low-bit decoder, even when the hint may describe any balanced pair.

The larger-factor interval has width O(√M). In each residue class its
index width is O(√M/2^k). For any fixed ε in (0,1/4), taking
2^k≥M^(1/4+ε) puts that width within the degree-one unknown-divisor
lattice regime (use slack ε/2 to absorb constants). Its affine
polynomial is monic after division by 2^k modulo the odd input.
Check roots for exact divisibility and balance before taking gcds.
This needs a supplied residue of a genuine sum on the transformed
input. Multiplying by 15 does not compute it.

Exhaustive direct congruence checks at k=8 on all 109 eligible
balanced inputs from 31 through 3999 coprime to 15 verified the
existence of a balanced pair and four residue classes containing the
larger factor for every such pair.

### Combining an archimedean bin with the low product bits

Let B≤p,q≤2B be primes with p≡1 and q≡3 modulo four, N=pq.
Choose uniformly among these pairs. Set m=2^k with 8≤m≤B. A guide
knows N modulo m and which of H intervals contains N. The interval
partition may even depend on the observed residue. If the guide
returns at most L guesses for (p+q) modulo m, its success rate is

    O(L (log B)² (1/m + H/B)).

Fix a product residue n, one interval of length Δ, and one guess s.
The low-bit decoding argument gives at most four possible classes for
p modulo m. Each has at most B/m+1 integer representatives. For each
such p, q must lie in the single class n/p modulo m and an interval of
length at most Δ/B. Thus there are at most 4(B/m+1)(1+Δ/(Bm))
pairs for this guess, even counting composite candidates. Within a
product residue the interval lengths sum to at most 3B² after clipping
to [B²,4B²]. Summing over the m residues and L guesses and dividing
by Θ(B²/(log B)²) prime pairs gives the bound. The denominator uses
the prime number theorem only in the two fixed reduced classes modulo
four. Randomized lists satisfy the same bound by conditioning.

For example m of order B^(1/2+δ), 0<δ<1/2, is large enough for the
preceding lattice threshold with a suitable fixed slack. Even H as
large as B^(1−η), fixed η>0, gives vanishing success for any
polylogarithmic list. This conclusion assumes the guide observes only
the bin and residue. A statistic computed from the full N can behave
differently within a bin.

As a finite check take B=8192 and m=128, and divide [B²,4B²) into
H equal intervals. Exhaustive counts on all 190080 eligible pairs
give optimal single-guess hit counts 15352, 16405, and 19821 for
H=1,8,64 respectively. In each (bin,residue) cell the best guess
is the most frequent sum residue. These are exact finite optima.

### The first ambiguous sum bit as quadratic residuosity

For distinct odd primes N=pq≡3 mod 4, label p≡1 and q≡3 modulo
four. Let a=−2 when N≡3 mod 8 and a=2 when N≡7 mod 8. Then

    p+q ≡ N+1     (mod 16)  if a is a square modulo N,
    p+q ≡ N+1+8   (mod 16)  otherwise.

Indeed N+1−(p+q)=(p−1)(q−1). The second factor has exactly one
factor of two. The product vanishes modulo 16 precisely for p≡1
modulo eight. The supplementary laws give Jacobi symbol (a/N)=+1
in both branches. Thus the two local Legendre symbols agree, and
(a/p)=(2/p) because −1 is a square modulo p. This common sign
is positive exactly when p≡1 modulo eight. CRT proves the formula.

Another expression for the same bit is the divisor sum of the
quadratic character χ(d)=(a/d) over positive d dividing N. The sum
is (1+χ(p))(1+χ(q)), hence four or zero. In the negative-sign choice
a=−2, the classical norm form x²+2y² represents N precisely in the
four case: each split prime is a norm in the Euclidean ring Z[√−2],
whereas an inert prime divides any norm to even order. For example
187=17·11 and 91=13·7 are both 11 modulo 16. The former equals
13²+2·3² and has sum 28; the latter has no such representation
and has sum 20. Exhaustively checking squares modulo N verified the
bit formula on all 143 eligible prime pairs below 100.

This identifies one missing bit on the semiprime promise with a
specific residuosity decision on the public input. A method to evaluate
the bit quickly, or to obtain the further bits needed for the lattice,
does not follow from the Jacobi symbol or the divisor-sum identity.
The statement does not automatically transfer to 15N with more than
two distinct prime factors.

### A single quartic test for the first ambiguous bit

The [earlier rational quartic audit](../../notes/newfamily_results.md)
already identifies the hidden-prime character in a quartic test on
squares. In the preceding semiprime case the exact specialization is

    p+q ≡ N+1 (mod 16)  iff  x⁴ ≡ 4 (mod N) has a solution.

For any unit t, locally x⁴=t² requires x²=t or x²=−t. At q≡3
modulo four exactly one sign is square. At p≡1 modulo four both signs
are squares precisely when (t/p)=+1, otherwise neither is. Taking t=2
and applying the preceding sum-bit formula proves the specialization.
Exhaustive fourth-power tables verified the general t identity on
23936 units across all eligible prime pairs below 50.

For the fixed target there are exactly eight roots on each yes input.
Four choose distinct signs at p and q, so gcd(x²−2,N) splits N for
exactly half of the roots. This extraction needs a root with suitable
signs. A decision bit alone does not provide it. For general t, an
ordinary quadratic-residuosity oracle on N can answer the quartic
question: choose u from ±t to have Jacobi symbol +1 (possible since
(−1/N)=−1). Its two local signs agree and (u/p)=(t/p). Conversely
the quartic test on t² decides whether t is a square modulo N whenever
(t/N)=+1; a negative Jacobi symbol already excludes a square. These
oracle reductions concern variable t. They imply no such equivalence
for the single target 4.

### Where the norm conic first detects the missing bit

Keep the preceding semiprime promise and sign a. Consider primitive
solutions of x²−Ny²=a z² modulo m, meaning gcd(x,y,z,m)=1.
On the yes branch there is a primitive solution for every m. On the
no branch the exact obstruction is p² dividing m or q² dividing m.

For an odd prime ℓ away from N, set z=1. The two sets of residues
x² and a+Ny² each have (ℓ+1)/2 elements, so they intersect. At a
solution x and y cannot both vanish, and a nonzero partial derivative
lifts the solution through all powers of ℓ. At 2 take y=z=1:
N+a≡1 modulo eight, so it is a square at all binary precisions.

At a prime r dividing N, the two local Legendre symbols agree by the
sign choice for a. On the yes branch a square root of a lifts, giving
a solution with y=0 and z=1. On the no branch reduction modulo r
forces x=z=0 modulo r, so (0,1,0) is still a primitive solution at
that precision. Modulo r² the equation now forces N y²=0 modulo r².
Since r divides N exactly once, this forces r to divide y as well.
CRT combines the local criteria. There is also a real point in both
cases. The Hasse–Minkowski theorem therefore identifies rational
solubility with the yes branch.

For N=91 and a=−2, every squarefree modulus passes. The first failing
modulus is 49. Direct enumeration checked the criterion on 91,187,39,55
at moduli 2 through 35 and at 49 (140 checks). One can form m=N²
without the factors, and primitive solubility there decides the bit
exactly. Efficiently testing that composite modulus remains open here.
This calculation concerns existence; counting or constructing points
requires separate analysis.

### Counting the affine conic with opposite local signs

For N the product of two distinct odd primes, choose a public unit a
whose Jacobi symbol modulo N is −1. Sampling units finds one in two
trials on average. For N≡3 mod 4 the fixed choice a=−1 works.
Write r for the prime factor where a is a square, and u for the other
prime. In the chart y=1 of the norm conic, reduction modulo N gives
x²=a z². Count all ordered pairs of residues (x,z), including zero.

Modulo r the equation consists of two distinct lines through the
origin, giving 2r−1 points. Modulo u it has only the origin. Thus the
CRT count C modulo the composite (i.e. the integer number of points
over Z/NZ) is exactly 2r−1. Moreover, every nonzero pair has z=0
modulo u and z nonzero modulo r, so gcd(z,N)=u. The zero pair is
always available. A routine returning any other point already factors
the input. An exact count gives r=(C+1)/2. Knowing C modulo 2^(k+1)
gives r modulo 2^k. On balanced inputs k>(1/4+ε)log₂N suffices for
the affine divisor lattice at fixed ε>0, using smaller slack to absorb
the constant in r≥√(N/2).

For N=187 and a=−1 there are 33 pairs, with all 32 nonzero pairs
revealing 11 by gcd. Exhaustive enumeration checked the formula and
every nonzero witness on 58 combinations of semiprime and test unit,
using primes p<20, p<q<25 and a in {−1,2,−2,3} with Jacobi −1.
Drawing z uniformly hits a nonzero eligible z with probability
(r−1)/N, which is exponentially small in the bit length on balanced
inputs. An independent evaluation of the count or a way to find such
a z remains necessary.

### The conic Fourier transform away from its origin

For the preceding split-sign conic let F(h,k) be the unnormalized
sum of exp(2πi(hx+kz)/N) over its points. Write r for the split
prime and u for the inert prime. The inert component is the origin.
Over r each line contributes r when its character is trivial and zero
otherwise, by geometric summation. CRT rescales both local frequency
coordinates by the same unit. Subtracting the doubly counted origin
therefore gives the exact integer values

* 2r−1 if h=k=0 modulo r,
* r−1 if k²=a h² modulo r and (h,k) is nonzero modulo r,
* −1 otherwise.

Every nonzero frequency modulo N can be evaluated in polynomial bit
time from the public inputs. First take gcd(h,k,N). If proper, use
the resulting factorization and the formula. Otherwise this gcd is
one. Take gcd(k²−a h²,N); a proper gcd again suffices, and a gcd of
one gives −1. The gcd cannot be N in this branch since the equation
at the inert prime would force both coordinates to vanish there.
Thus every nonzero frequency with value different from −1 supplies
a factor through these gcds. Their fraction among all frequencies,
including the origin as exceptional, is (2r−1)/r².

At (0,0) the value is the desired count. Its nearest nonzero grid
neighbor (1,0) always has value −1. Continuity on the underlying
real frequency axis does allow recovery from a much finer point:
using representatives 0≤x<N, the real part of F(1/N,0) differs
from F(0,0) by at most (2r−1)·2π²/N² by 1−cos(t)≤t²/2.
But that phase is not a character modulo N, so the off-origin
integer-frequency evaluator does not apply there. For N=187, a=−1,
a direct numerical 2D FFT reproduced 30976 entries at −1, 3872 at
16 and 121 at 33, with maximum deviation below 10⁻¹³.

### Exact fractional frequencies on a coordinate axis

Take canonical coordinates in [0,N). At the inert prime u both
coordinates must vanish. At the split prime r the two distinct lines
give exactly two choices of z for each nonzero x, and one at x=0.
Consequently the projection to x has multiplicity two on the nonzero
multiples of u and multiplicity one at zero. For all real or complex t,

    F(t,0) = 2 Σ_{j=0}^{r−1} exp(2πitj/r) − 1.

For t outside r times the integers this is
2(1−exp(2πit))/(1−exp(2πit/r))−1. At a removable singularity the
value is 2r−1. In particular F(1/2,0)=1+2i cot(π/(2r)).
The real part at this half frequency is constant; the imaginary part
determines r. Its derivative with respect to real r is
π csc²(π/(2r))/r² ≥ 4/π, by sin(x)≤x. Thus even absolute error
1/2 in the imaginary part allows unique integer recovery, e.g. by
monotone bisection with sufficient numerical guard precision.
On N=187 direct summation gave 1+21.5834373i at t=1/2 and real
part 32.9941563 at t=1/N, consistent with r=17.

For this axis the apparent analytic continuation problem is exactly
the transform of an arithmetic progression of hidden spacing N/r.
All integer probes before the first nonzero multiple of r return −1.
The geometric series evaluates the transform given r. A useful advance
here would need an independent computation from N and a; changing
frequency or sharpening the continuity estimate does not supply one.

### Uniform sampling on the norm-one conic and the accuracy needed

Keep (a/N)=−1 with split prime r and inert prime u, and replace zero
by one in x²−az²=1. This is the [previously audited norm-one
group](../../notes/codex_r1_mechanisms.md).
Its local orders are r−1 and u+1: in the split algebra a pair of units
must have product one; in the quadratic field the norm has kernel of
size u+1. Hence the cardinality C is N−1+r−u.

Exact uniform sampling uses random v,w modulo N, keeping pairs for
which d=v²−aw² is invertible. Output ((v²+aw²)/d,2vw/d) modulo N.
This is the quadratic-algebra homomorphism α↦α/conjugate(α). Its
kernel has (r−1)(u−1) elements, the scalar units. Its image has
(r−1)(u+1) elements by division of the unit-group size, so every
norm-one point occurs equally often. Acceptance probability is
(1−1/r)²(1−1/u²), bounded below by a positive constant. On N=187
with a=−1 exhaustive enumeration gives 192 outputs, each with
160 preimages.

An independently computed confidence interval for C of width W gives
an interval of width at most W for r, even when the primes are very
close. In fact D=C−(N−1)=r−u and
r=(D+√(D²+4N))/2 has derivative strictly between zero and one
everywhere. The same Lipschitz bound holds for u as a function of D.
On balanced inputs W=O(N^(1/4−ε)) for fixed 0<ε<1/4 therefore
suffices for the existing unknown-divisor lattice (use smaller slack
for fixed constants). Exact uniform sampling alone supplies no such
estimate. In particular the probability of any repeated output in m
independent draws is at most m(m−1)/(2C), negligible at polylogarithmic
m for balanced inputs. This addresses collision-based size estimation
only. The earlier norm-one audit already treats public exponentiation.

### Interval bias of sampled norm-one coordinates

Take a uniform point (x,z) on the norm-one conic over N=ru above.
For any public h,k with gcd(h,k,N)=1 and any interval I of canonical
residues, the probability that hx+kz modulo N lies in I differs
from |I|/N by O(log N / √N), uniformly in the interval and projection.

Here is the Fourier estimate. At an odd prime ℓ dividing N, write
e(v)=exp(2πiv/ℓ). For a nonzero frequency pair (H,K), insert additive
orthogonality for x²−az²=1. The zero multiplier contributes zero.
Completing the two quadratic Gauss sums leaves, up to a factor of
absolute value one, the sum over nonzero v of

    e(−v − (H²−K²/a)/(4v)).

The classical Weil bound for a prime-field Kloosterman sum bounds its
magnitude by 2√ℓ. If the reciprocal coefficient vanishes the sum is
−1 directly. At zero frequency the local count is ℓ±1. CRT factors
the global Fourier coefficients with invertible local rescalings.

In the one-dimensional interval Fourier expansion the zero frequency
gives C|I|/N, where C is the conic size. For frequencies coprime to N
the conic coefficient has magnitude at most 4√N, since the projection
is primitive at both primes. The sum of interval coefficient magnitudes
is O(N log N), by the geometric series bound. After the inversion
factor 1/N this gives O(√N log N) count error. For frequencies rj,
1≤j<u, the conic coefficients are O(r√u), and the interval coefficient
is O(u/min(j,u−j)), independently of the number of full periods in I.
These terms contribute O(√u log u) after inversion. Interchanging r,u
handles the remaining frequencies. Divide by C=(r−1)(u+1)=Θ(N).

Consequently an observation giving one of H interval bins of one such
projection has total variation distance O(H log N / √N) from the bin
of a uniform residue. For m independent draws with the next projection
and bins chosen from the preceding bin labels, coupling bounds the
transcript distance by O(mH log N / √N). Thus polylogarithmically
many such observations can be simulated with independent uniform
residues up to negligible error. This does not cover multiple
correlated projections of the same point or tests on exact coordinates.

For N=187,a=−1, direct enumeration over the 192 points gave maximum
absolute interval deviations 1315/17952, 967/8976, and 39/1496 for
projections x, x+z, and 3x+5z respectively, taking the range of the
centered cumulative counts.

### The square-class bit of a sampled norm-one point

Keep the sampler's unit α=v+w√a with norm d=v²−aw², and write
g=α/conjugate(α)=x+z√a. Since x+1=2v²/d and x−1=2aw²/d,
for unit v,w the Jacobi symbols of 2(x+1) and 2(x−1) are respectively
(d/N) and (a/N)(d/N). They are opposite here. The label (d/N)
also works when a coordinate vanishes, and changing the sampler
preimage multiplies d by a scalar square.

Over each local ring at r or u, g is a square in the norm-one group
exactly when d is a scalar square. In fact this gives an explicit
bijection over the composite too: for each scalar t with t²=d,
the group element α/t has norm one and square g. Conversely if h
has norm one and h²=g, then α/h equals its conjugate and hence is
a scalar (2 is invertible); its square is d. Thus extracting halves
of a sampled point with its preimage in hand amounts to extracting
ordinary modular square roots of its norm.

On a uniform point the two local square-class signs are independent
fair bits: CRT makes the local components independent and each local
cyclic group has even order. Their product is the public Jacobi label.
Restricting to either label by rejection sampling takes two trials
on average and leaves the two compatible local sign pairs equiprobable.
Conditioning on the label alone can increase the probability of any
event by at most two. The complete point can contain other information.

On N=187 and a=−1, direct group squaring and enumeration of sampler
preimages gave 48 distinct points in each of the four local sign
pairs and agreed with both local norm characters for every preimage.

### Prescribing a norm by choosing the conic

For any unit d modulo an odd semiprime N=pq with distinct primes,
choose v uniformly modulo N and set a=v²−d. Keep a when its Jacobi
symbol is −1. Then α=v+√a is a unit of norm exactly d, and
g=α/conjugate(α) is a publicly constructible norm-one point. The
preceding bijection sends square roots of d to halves of g and back.

Write e_ℓ for the number (zero or two) of square roots of d modulo
ℓ. The exact number of successful v is ((p−e_p)(q−e_q)−1)/2.
Indeed the number giving nonzero Jacobi symbol is (p−e_p)(q−e_q).
The sum of the Jacobi symbols over v is +1 by CRT: each local
character sum is −1 since (v−y)(v+y)=d has exactly ℓ−1 solutions.
Thus the probability is at least ((1−2/p)(1−2/q)−1/N)/2 ≥ 1/15,
using p≥3, q≥5 after ordering. Construction takes at most fifteen
trials in expectation. Enumeration for all 160 units d modulo 187
gave 67, 76, 82 or 93 successful v out of 187, matching the formula.

Conversely a point g=x+z√a presented without its sampler preimage
admits α=1+g whenever gcd(2+2x,N)=1, since its norm is 2+2x
and α/conjugate(α)=g. A proper gcd splits N. If x=−1 modulo N
then z=0 since N is squarefree, and α=√a works. Thus, with the
conic parameter allowed to vary, square decision for norm-one points
and ordinary modular quadratic residuosity reduce to each other in
randomized polynomial time (after a split, local decisions are easy).
An algorithm returning a half on every square input would also factor
N: plant d=t² for uniform unit t, recover a root t′ from the returned
half, and compute gcd(t−t′,N). Conditional on d, the constructed
instance is independent of which of its four roots was planted, so
the gcd is proper with probability one half. The decision equivalence
alone gives no such factoring reduction. Fixed-a algorithms have a
different promise and the forward construction need not query them.

### Halving uniform squares on the fixed circle

For search, even the single circle x²+z²=1 (mod N) suffices for all
odd composites with at least two distinct prime factors. Use complex
multiplication modulo N as the group law. Plant a uniform point H and
give a solver only G=H². This is uniform in the subgroup of squares.

Uniform sampling still works with prime powers. Choose random α=v+iw
with invertible norm and output α/conjugate(α). To see surjectivity
over each odd prime power, use α=1+H if its norm is a unit, otherwise
use α=i(1−H). The two norms are 2(1+x) and 2(1−x) for H=x+iz,
so at least one is a unit locally. Both choices have quotient H.
CRT yields surjectivity globally, and equal fibers give uniformity.
For each prime p dividing N, the excluded norm-zero pairs modulo p
occupy either one point or two lines through the origin. Acceptance
is therefore at least the product of (1−1/p)², hence at least
1/(ω(N)+1)² by ordering the primes and telescoping. This bound
suffices for expected polynomial bit time.

For any verified half H′ of G, the ratio H′/H squares to one. Since
its inverse equals its conjugate, it is a scalar t, with t²=1 mod N.
There are 2^ω(N) such scalars, from independent signs on the prime
powers. Conditional on G and the solver output, the planted H is
uniform in the fiber, so gcd(t−1,N) is proper with probability
1−2^(1−ω(N)). If the solver succeeds with probability η on uniform
squares for this N, one trial splits N with probability exactly
η times that factor (ignoring any earlier accidental gcd). Detect
primes and perfect powers separately, strip powers of two, and recurse
on proper splits. A per-input inverse-polynomial η thus suffices for
randomized expected polynomial-time factoring if the solver itself is
polynomial time.

Enumeration with a deterministic canonical returned half for each G
gave 8/16, 24/48, 96/128 and 120/240 successful planted points for
N=15,45,105,225 respectively. This is the classical planted-root
reduction applied to the circle. Producing some point of a prescribed
norm and halving a specified group point are different tasks; the
former alone gives no half of G.

### Forcing a split from a unique odd-order half

For two distinct odd primes N=pq, sample a uniformly among units with
Jacobi symbol (−a/N)=+1. Exactly half these parameters make the square
subgroup of the norm-one group odd in order. Indeed at ℓ=p or q the
cyclic local order ℓ−(a/ℓ) is 2 modulo 4 exactly when (−a/ℓ)=−1.
The public filter leaves the two patterns (+,+) and (−,−) equally
likely. Recognizing the good pattern itself is quadratic nonresiduosity
of −a under the positive Jacobi promise; recognizing it is unnecessary
for the following randomized reduction.

On a good parameter let S be the subgroup of squares. Each G in S
has a unique half H₀ inside S, since squaring is a bijection there.
Plant a uniform H conditioned on the negative public square-class
label from the sampler above, then submit G=H². At each local prime
the scalar −1 is a nonsquare in its local group (order 2 modulo 4).
Consequently exactly two of the four halves of each G have negative
label. They differ from H₀ by the two mixed CRT signs. This proves
both that the submitted G is uniform in S and that comparing H₀ with
the planted H gives a proper gcd every time.

Thus a polynomial-time solver returning the unique odd-order half
with probability η averaged over good parameters and uniform targets
at each semiprime gives success probability at least η/2 per trial.
Allow arbitrary outputs on bad parameters and verify the gcd before
accepting. This is a promise on the solver's correctness, with a
polynomial runtime cap on all calls. It restricts the search problem
to uniquely specified outputs while retaining a factoring reduction.
It supplies no algorithm to invert squaring inside S.

For N=187 exactly 40 of 80 eligible parameters are good. At a=3 the
group has 180 points and S has 45. Exhaustive squaring found precisely
two negative-label planted points per target and all 90 comparisons
with the unique half split N.

### Coverage of affine exponents for unique halving

On a favorable conic the local square subgroup at p is cyclic of odd
order A=(p+σ)/2, where σ=1 for p≡1 mod 4 and −1 otherwise.
For uniform G in the square subgroup, G^e is its correct half at p
with probability gcd(A,2e−1)/A. Thus an exponent that works for
every G must invert 2 modulo both local square-subgroup orders.

Take N=pq uniformly over distinct primes in [B,2B]. At each scale
allow a shared menu of poly(log B) formulas e=(cN+d)/s with integer
coefficients of magnitude at most B^K for fixed K and s>0. Use a
formula only on inputs where it gives an integer. The probability
that ANY formula gives the correct half at even ONE prime component
of a uniform square on a favorable conic is at most B^(−1+o(1)).
Selection from the menu may use N and the observed target.

Proof for one formula and one labeled prime: replace 2e−1 by its
integer multiple s(2e−1) to upper-bound the gcd. Modulo A this becomes
Cq+D, with C=−2cσ and D=2d−s. If C=D=0 the formula is never integral.
Otherwise put g=gcd(|C|,|D|)>0. Expand gcd(A,Cq+D) as a sum of φ(h)
over common divisors. For each h≤B+1 and fixed σ, at most O(B/h+1)
integers p in range satisfy h|A, since p≡−σ mod 2h. The congruence
Cq≡−D mod h has no solutions unless gcd(C,h)|D, and otherwise
at most O(B gcd(g,h)/h+1) integers q in range solve it. Counting all
integers here gives a valid upper bound on prime pairs. The weight
φ(h)/A is O(h/B). Summing the resulting products gives
O(B log B τ(g)+B), since Σ_{h≤B+1} gcd(g,h)/h ≤ O(τ(g) log B).
Now τ(g)=B^o(1) by the classical divisor bound. Divide by the
Θ(B²/log² B) prime pairs and sum over signs, components and formulas.

As a numerical check, for e=(N+1)/2 the mean probabilities of a
correct half at either component were 0.021805, 0.005398 and 0.001358
at B=128,512,2048, from all 253,2775,32385 unordered pairs respectively
and the exact cyclic-group probabilities. The argument concerns these
affine exponent menus; it leaves other ways to compute the unique
half open, including exponents obtained from further arithmetic of N.

### An interval logarithm for the signed prime difference

Let r,u be the split and inert primes of a mixed norm-one conic,
with max(r,u)≤2 min(r,u). Write D=r−u. Its group order N−1+D
gives the public target relation G^D=G^(1−N) for every point G.
The same holds for a tuple of independent points. We have
|D|≤W=floor(√(N/2)), by the balance promise. A returned candidate
d can be checked by testing d²+4N for a square and verifying the
resulting factors. Checking both signs is inexpensive.

Two independent uniform points, viewed as one tuple in a product
group, have order equal to λ=lcm(r−1,u+1) with probability at least
1/ζ(2)=6/π². For each prime ℓ dividing λ choose a local cyclic
component attaining the maximum ℓ-adic valuation in λ. At least one
of the two local samples has full ℓ-part there with probability
1−1/ℓ². The events across ℓ are independent by primary decomposition.
Their product bounds the success probability from below. This proof
does not require computing λ or recognizing successful tuples.

On success the tuple order exceeds W, since λ≥max(r−1,u+1)≥√N−1
and √N−1>√(N/2) for N≥15. Thus each of the two signed windows
[−W,0], [0,W] contains at most one solution, and one contains D.
An interval-log solver may report any solution in each window or
report emptiness. Verify both answers, and resample on failure.
For N=187,a=−1 the local orders are 16 and 12, the target D is 6,
and 24576 of 36864 ordered sample pairs have full exponent 48.

Classical baby-step giant-step on the tuple gives a reference search
of O(√W) group operations and stored entries per window, using
canonical residue coordinates for equality. This is the N^(1/4)
scale and supplies no polylogarithmic discrete logarithm. The tuple
argument isolates the computational cost from possible ambiguity of
the signed difference.

### A known order on both components shortens the interval logarithm

Suppose a conic point has known exact order m with known prime
factorization, on the mixed conic for N=ru from the preceding
interval-log reduction. For each prime ℓ dividing m, raise the point
to m/ℓ. This cannot be the global identity. If it is the identity in
one component, gcd(N,x−1,z) of the resulting point gives a split.
Otherwise both local orders contain the full ℓ-part of m. After all
tests succeed without a split, m divides both r−1 and u+1. Hence

    m² divides (r−1)(u+1) = N−1+D, where D=r−u.

Choose d₀≡1−N modulo m² and write D=d₀+m²j. For the public tuple G,
solve (G^(m²))^j=G^(1−N−d₀). Each signed half-window now has
O(1+√N/m²) candidate indices. On the full-exponent tuple event
the original uniqueness in that window also gives uniqueness after
substitution. Baby-step giant-step uses O(1+N^(1/4)/m) group
operations after logarithmic-cost exponent setup.

For an example requiring no search for torsion, take a=−3 on odd
semiprimes N≡2 mod 3 coprime to 3. The point (1/2,1/2) has order six
in both components: it satisfies T²−T+1=0 and its cube is −1, with
neither component of order two. The split prime is 1 modulo three by
quadratic reciprocity. Thus D≡1−N modulo 36. At N=221 only −4 lies
in the balanced difference window with that residue, recovering 13
and 17. All 223 prime pairs 5≤p<150, p<q<2p with pq≡2 mod 3 verified
the congruence.

Large common orders have a distributional limitation independent of
how the conic is selected. For primes in [B,2B], any common order
belongs to one of the four gcds gcd(p±1,q±1). For each 2≤T≤2B+1
the fraction of pairs admitting any common order at least T is
O(log²(B)/T): union-bound over signs and divisors T≤m≤2B+1, with
O((B/m+1)²) integer pairs per divisor, then divide by the
Θ(B²/log² B) distinct prime pairs. This rules out a power-size
common order on most balanced semiprimes. It leaves open how often
one can construct a point of known order whose local orders differ,
which is the immediate-split case above.

### Reducing a characteristic-zero torsion point

A primitive complex mth root of unity satisfies Φ_m. At every
prime ℓ not dividing m, its image in any residue field where it is
reduced still has exact order m. Indeed X^m−1 is squarefree in that
characteristic. Its cyclotomic factors Φ_d for d dividing m therefore
have disjoint roots, and a root of Φ_m cannot be a root of X^d−1
for any proper divisor d. The same proof works after adjoining
the root formally and projecting to any field component. It does not
depend on which embedding or conjugate one chooses. Local splitting
degrees may still differ.

On a rational quadratic norm-one conic a torsion unit lies in a
quadratic field or in the split algebra. In the field case φ(m)≤2,
hence m is 1,2,3,4 or 6. In the split case rational roots of unity
are ±1. Thus reducing one such point, with coordinates defined at the
primes and nondegenerate conic there, places it in the equal-order
case whenever gcd(m,N)=1. This explains why the explicit order-six
point gives a congruence but cannot trigger the order-mismatch gcd
on the coprime inputs. It places no restriction on a modular point
constructed directly using N and only known to satisfy X^m=1.

For N=221, the modular point (137,85) on x²+3z²=1 illustrates the
distinction. It equals the primitive sixth root (1/2,1/2) modulo 13
and the identity modulo 17; gcd(137−1,85,221)=17. This example was
constructed using the factors, so it supplies a check only. Polynomial
gcd checks verified disjointness from every proper-order factor for
288 coprime (m,ℓ) pairs with 2≤m≤23 and odd primes ℓ<50.

### Projecting the odd-prime sum sieve onto low bits

For m=2^k the locally possible endings are the image of f(x)=x+N/x
on odd x modulo m. For an input x let v be the 2-adic valuation
of x−N/x, truncated at k. The fiber equation is
(y−x)(y−N/x)=0 modulo m. When 2v<k it gives two disjoint classes
modulo 2^(k−v), hence 2^(v+1) inputs in the fiber. This follows
since the two differences have separation of valuation v: if their
product vanishes modulo m, one has valuation v and the other at
least k−v. When 2v≥k there is just one class modulo 2^ceil(k/2),
of size 2^floor(k/2). The number of x with v≥j equals 2^(k−j)
times the number of square roots of N modulo 2^j, since x is odd.
Summing the reciprocal fiber sizes gives the image size in O(k)
integer operations. At k=26 the counts for N≡1,3,5,7 modulo 8
are respectively 1,398,104, 8,388,608, 4,194,304, and 8,388,608.

For the 97-bit benchmark N=118464083331556819902753234409, *all*
1,398,104 endings survive the earlier odd-prime tests through 31
even after imposing the balanced sum interval. Here is a certificate.
Write s=a+2^26 t for a proposed ending a. Modulo the odd sieve
period P=100280245065, testing s is equivalent to testing
2^26(t+a/2^26), using inversion modulo P. Thus the same set of
acceptable positions is translated cyclically for each a. Each
ending has at least 622,224 consecutive choices of t in the sum
interval. Partition P into 322,329 adjacent blocks whose lengths
are at most 311,112. The checker finds a passing position in each
block. Every cyclic interval of 622,224 positions contains an entire
block, including at wraparound. Consequently each ending has a
passing lift. This is a benchmark-specific certificate and does
not apply to extra predicates or narrower sum intervals.

### Which ranks can see a split quadratic algebra?

For odd N=pq and (a/N)=−1, the algebra A=R[T]/(T²−a) splits
at exactly one prime. Suppose a linear map constructed from the
unmarked algebra is invariant up to change of bases under algebra
isomorphism, and its construction commutes with extension of the
base field. At each prime ℓ its rank equals that for the split
algebra F_ℓ×F_ℓ: after extending to the algebraic closure the
algebras are isomorphic and matrix rank is unchanged. A rank gap
between the two characteristics in such a construction already
occurs for the public split-algebra control. Multiplication and
trace give a concrete example: the trace pairing has Gram matrix
diag(2,2a), of rank two on both sides. Its determinant square class
does distinguish the local algebras, with its orientation hidden.

There is also a bound permitting a marked generator and explicit
dependence on a. Fix N before sampling a. Let M(U) be a public
matrix with polynomial entries in U of degree at most d. Suppose
the generic ranks over F_p(U) and F_q(U) agree and equal r. At
each prime choose one generically nonzero r-minor (r=0 is trivial).
Its degree is at most rd, so the specialized rank can drop on at
most rd field elements. A uniform unit a conditioned on Jacobi −1
has uniform projection to each local unit group, since either
local sign has equally many partners of the required opposite sign.
Thus the probability of unequal ranks is bounded by
min(1, rd/(p−1)+rd/(q−1)). Fixed lists of matrices obey the union
bound. For unequal generic ranks, unconditioned random unit
specialization preserves both generic ranks except with probability
bounded by the corresponding sum of minor-degree bounds; that route
already supplies a rank gap without selecting split versus inert.

For multiplication by T−1 the determinant is 1−a. The rank drops
precisely at a=1 locally, and the two drops are mutually exclusive
under the Jacobi condition. Its exact gap probability is therefore
1/(p−1)+1/(q−1). Direct enumeration verified this on all 21 distinct
semiprimes with primes from 3 through 19. These arguments do not
bound the degree of a compact exponentiation circuit or cover a
matrix whose entries encode a ground-field point search.

### Scalar squares in the quadratic algebra and the hidden character

Fix a Jacobi-negative unit a for a distinct odd semiprime N, and
label the split prime p and the inert prime q. For a scalar unit b,
the class of b is a square in A=(Z/N)[T]/(T²−a) exactly when
(b/p)=+1. Over p the two coordinates of the scalar are (b,b).
Over q every scalar unit is a square in F_(q²), since q+1 is even.
This gives the oriented local character as a ground-field property.

There is an exact public reduction to ordinary quadratic residuosity
on Jacobi-positive inputs. Form t=b if (b/N)=+1, otherwise t=ab.
Then (t/N)=+1 and (t/p)=(b/p), since (a/p)=+1. Thus t is a
square modulo N exactly when b is a square in A. Conversely a
Jacobi-positive input is a square modulo N exactly when it is a
scalar square in A. This is decision equivalence; it does not
supply a polynomial-time algorithm for quadratic residuosity.
Direct enumeration of pairs x,y with 2xy=0 modulo N checked
3,192 (a,b) cases across the six semiprimes formed from 3,5,7,11.

Extraction is stronger in a concrete sense. The element 1 has eight
square roots in A, four choices at the split prime and two at the
inert prime. Each root x+yT other than the two public scalars ±1
gives a factor. If y=0, use gcd(x−1,N) on the mixed scalar root.
If y≠0, the inert-field root is ±1, hence q divides y and gcd(y,N)
is proper. Uniform sampling from these eight roots would split with
probability 3/4. The decision test by itself supplies no such sampler.

### How often an empty-interval certificate must be constructible

For N=pq with distinct odd primes, the even sums 2√N≤s<N+1 contain exactly one square discriminant: a square gives the integer factor pair (s±√(s²−4N))/2, and the upper bound excludes (1,N). A complete deterministic emptiness test is stronger than needed. Bisect this list of even sums and run independent polynomial-time certificate attempts on both children. Accept only certificates whose sound verifier certifies emptiness, then keep the other child. It suffices that on each empty sibling on this particular bisection path a fresh attempt succeeds with probability at least ε=1/poly(log N). At most ceil(log₂ N) splits and at most 2/ε expected attempts per split give a Las Vegas reduction. The last sum gives both factors by integer square root. Conversely the two primes are themselves a short, deterministically checkable emptiness certificate for any empty interval under this promise. Certificate length alone therefore cannot demonstrate progress. The outstanding requirement is a public constructor with a success guarantee on those siblings. Average success over unrelated empty intervals need not give this guarantee.

### Adaptive rational slopes in a contiguous divisor scan

Consider the explicit scheme covering [X,2X], with X comparable to √N, by intervals on each of which the full real range of F(t)=a t+qN/t has length at most K≥1, for independently chosen integers a and q≥1. This permits testing O(K+1) integer levels by solving quadratics, as in the [Hiary implementation](../../src/analytic/hiary.py). Even adaptive choices require Ω(X^(2/3)/K^(2/3)) blocks for 1≤K≤√X, with constants depending on the balance bounds.

Indeed F'' is comparable to q/X. The midpoint second difference on a block of width w forces qw²/X=O(K). The mean value theorem supplies a point with |F'|≤K/w, and variation of F' across the block is O(qw/X)=O(K/w). Thus throughout any block of width ≥H the slope N/t² is within O(K/(qH)) of a/q and q≤Q=O(KX/H²). This slope ranges over a bounded interval, so there are O(q) eligible numerators per denominator (take H≥K). The union of these rational neighborhoods has length O(KQ/H)=O(K²X/H³). Since |(N/t²)'| is comparable to 1/X, taking H a sufficiently large constant times (K²X)^(1/3) leaves at least half the original interval outside *every* admissible block of width ≥H. Covering the remainder takes Ω(X/H) short blocks. The same argument applies on any fixed proportional balanced subinterval. This makes the generic-denominator heuristic in [the analytic notes](../../notes/analytic_results.md) unnecessary for this specific cover model. Arithmetic filters that collectively skip integer levels or distant blocks are outside the argument.

### Several weighted sums from the same modular factor pair

For an odd modulus m, a candidate sum s passes the exact square-discriminant sieve precisely when u+v=s and uv=N have a solution modulo m (take u,v=(s±√(s²−4N))/2). The same pair simultaneously supplies every weighted sum t=au+bv and a square root of its discriminant: t²−4abN=(au−bv)². Consequently adding any number of such congruences with freely chosen witness values t leaves exactly the original candidates, even when all weights must share one ordered pair. CRT gives the same conclusion across coprime moduli. On the balanced benchmark N=2989601243, exhaustive enumeration through prime 31 gives survivors 112716, 114144 and 115446. The latter two have nonsquare integer discriminants and pass every weighted test of this form at the same moduli. Independently measured residues or interval bounds on t impose additional conditions and are outside this equivalence. See the [prior weighted sieve comparison](../../notes/signsplit_results.md).

### Coverage of the small Lehman multiplier benchmark

Take primes X≤p≤2X and 1.35p≤q≤1.65p. Search every multiplier k≤K at the first H integer centers starting at ceil(2√(kpq)), where K,H are positive integers polynomial in log X. A successful proper gcd from a square at such a center implies k=ij and center=ip+jq for positive integers i,j: the two roots of the associated monic quadratic multiply to kpq, and a proper split puts p and q in different roots (for large X, k<p; odd factors make extracting the gcd with twice the root harmless). The excess above 2√(kpq) is p(√i−√(j q/p))² < H. Hence i is comparable to j, j=O(√K), and q/p lies within O(√(H/(Xj))) of i/j. For each j there are O(j) numerators. For fixed p, counting all integers q in these intervals gives O(√(HX) K^(3/4)+K) possibilities, including endpoint rounding. The prime number theorem gives Θ(X/log X) prime candidates in the original q interval uniformly for these p. The hit fraction is therefore O(log X [√(H/X) K^(3/4)+K/X]), tending to zero. This also covers input-dependent choices restricted to this same menu. The square test can succeed for planted ratios, as documented in the [benchmark](check_lehman_window_baseline.py). The counting bound addresses small multipliers and short consecutive windows only.

**Sparse menus of larger multipliers.** On the same family, fix a menu shared by all inputs at scale X. Allow M multipliers, each at most X^C for a fixed C, and H consecutive centers per multiplier as above. For polylogarithmic M,H the fraction split is at most X^(−1/2+o(1)). Indeed each k has at most τ(k)=X^o(1) ordered decompositions ij=k by the classical divisor bound, uniformly for k≤X^C. For fixed p, each decomposition supplies O(1+√(HX/j)) candidate integers q from the preceding ratio estimate. Summing over the menu and dividing by the number of candidate primes gives O(M max τ(k) (1+√(HX)) log X / X). Count all cases where a menu entry shares a prime with N as successes separately: the menu contains at most CM prime divisors ≥X, and each occurs in O(X/log X) input pairs out of Θ(X²/log² X). Random menus independent of the particular N obey the same average bound by conditioning. More generally the counting estimate shows that constant coverage with polylogarithmic H requires at least X^(1/2−o(1)) distinct entries across the family, even if an input-dependent selector tests only one of them. This constrains the range size of such a selector, not the computational cost of selecting an entry.
