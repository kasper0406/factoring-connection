# QRP, hidden conductors and interval-count targets (2026-11-25–26)

Consolidated from the conversation and five recovered early archives. The full early
archive inventory is in [experiment provenance](../experiments/nov26_consolidation/README.md).
Results from later scratch work are explicitly marked REPORTED when the underlying
proof or experimental record was not recovered. None of these gives a general fast
factoring algorithm or a reduction from factoring to ordinary QRP.

## 1. Exact character access, and the missing conductor decoder

For N=rq with distinct odd primes and a promised Jacobi-positive unit x,
QRP asks whether both local Legendre symbols are +1. Pick a public unit t
with Jacobi(t,N)=-1 and label r by (t/r)=+1. A QRP oracle gives (a/r)
for all units: query a if its Jacobi sign is +1 and ta otherwise.
Conversely either designated local character answers QRP. This is an equivalence
of decision capabilities, not a factor-extraction algorithm. A useful inverse-
polynomial advantage on uniformly random promised inputs at a fixed N can be
amplified via multiplication by random squares (random self-reduction).

Folded permutation proof (mixed moduli p=1 mod 4, q=3 mod 4): folded multiplication
by x on unit representatives in [1,(N-1)/2] has sign (x/p). The same sign arises
on the group of squares under multiplication by x^2 (the local square groups
have even and odd orders respectively). For an odd lift A of x, write W(k) for
the number of folded sums and differences from unordered pairs of unit
representatives. If C(k) counts z mod N with z,z+k units, pairing z with -z-k
gives 2W(k)=C(k)-1_{gcd(k,N)=1}. The CRT values for C(k) imply
W(k) mod 2 = 1_{p does not divide k}. Hence the permutation sign exponent is
sum_{1<=k< N/2, p∤k} floor(Ak/N) mod 2. Full floor-sum reciprocity leaves the
unknown denominator q. For the prefix with A=N+4 and K<N/4, subtract from
K(K+1)/2 to obtain T(floor(K/p)) mod 2, where T(n)=n(n+1)/2.
On the first doubling bracket crossing p this is a threshold, so chosen-prefix
parity factors N. Handle q=3,7 separately. No reduction from the complete sign
to these prefixes is known.

Conductor wheel: reciprocity turns queries (ell/r) into constraints
(r/ell) for known small primes ell, after calibration at -1. A product W
of sufficiently many small primes identifies r under Dirichlet GRH
(pseudosquare uniqueness); a polynomial-size random large-prime wheel
also has a finite-height uniqueness guarantee (uniform effective Chebotarev).
Neither is a decoding algorithm. Explicit consecutive-prime wheel enumeration
has B^{1-o(1)} lifts for an interval of length B. In the early toy experiment,
30/30 14–30-bit factor transcripts selected one integer before divisibility
verification; median explicit lifts grew 53.5, 393.5, 2957, 22420, 171523.
The public complement character is determined by N and one hidden character.
Random large-prime transcripts can be simulated by choosing symbols and testing
public Jacobi consistency; the prior simulation argument is for that random
wheel, not the consecutive small-prime case. See also
[wheelgap](wheelgap_results.md) and [factorbits](factorbits_results.md).
The planted high-agreement CRT/Coppersmith cases (including an 81-bit example)
are valid conditional decoders, not a generic algorithm: they rely on atypical
agreement of residues of p and q. A uniform consecutive-prime decoder remains open.

Measured QRP screen (original preregistration and results recovered): 27,648
held-out samples, 72 distinct held-out moduli across 40/60/80 bits and three
mod-4 strata, 230 public features. Logistic accuracy 0.50420, extra trees
0.49504, global stump 0.49877, within-modulus stump 0.50145; all adjusted tests
null. Planted 55% feature control achieved 0.54427. This rules out the tested
large signals, not arbitrary inverse-polynomial asymptotic advantages.
See [archived run summary](../experiments/nov26_consolidation/qrp_screen.json).

## 2. Genus, theta and higher-residuosity views

Let N=rq, r=1 mod 8, q=3 mod 8. Supplementary quadratic characters give the
mod-8 class, not the next bit. At r, (-1)^((r-1)/8)=( (2+sqrt(2))/r )
(the sign is independent of the root). A fixed order-four character is needed
to distinguish r=1,9 mod 16 among rational residue classes; this does not
bound N-dependent quadratic queries. See [hidclass](hidclass_results.md).

Early archived results: for squarefree m=3 mod 4, m>3,
S(m)=sum_{0<a<m/2} Jacobi(a,m)=(2-Jacobi(2,m))h(-m).
Genus theory at three/four prime discriminants recovers QRP by a constant
number of padded half-sum residues mod 8 or 16. The four-prime comparison:
choose b=1 mod 4 prime with Jacobi(b,N)=-1; choose distinct eligible primes
c_i=-N mod 4 with (c_i/b)=-1 and a common Jacobi(c_i,N) sign. Then
S(Nbc_0)=S(Nbc_1) mod 16 iff c_1/c_0 is square mod N.
The archived proof uses the Schur complement of the two-auxiliary-prime
block of the Redei matrix. Sampling random square cosets and eligible primes
gives the stated one-sided preprocessing; audit quantitative prime density
assumptions separately from the unconditional algebra.

The fifth dissection is
    4 sum_{0<a<m/5} Jacobi(a,m)=(5-Jacobi(5,m))h(-m)-h(-5m)
for 5∤m. With (5/N)=-1 and appropriate c_i, two fifth-sums mod 4 answer
the comparison; dropping the twist without its divisibility hypothesis is
wrong (m=5863). More generally set k=[(5/N)=-1], g=[(c_i/5)=-1] and hold
the public Jacobi sign j=[(c_i/N)=-1] fixed. The archived difference rule is
(Q(Nc_1)-Q(Nc_0))/2 = kg [(c_1c_0/p)=-1] mod 2.
This shortens the interval by a constant, not its asymptotic complexity.

Other reported reductions (not freshly audited): elliptic-curve twist point
counts mod powers of two and ternary theta counts recover low factor bits;
in particular U(N)=#{0<x<y<z: x^2+y^2+z^2=N^2} satisfies 8U=N+2r-21 on
this family (also in hidclass). A uniformly supplied cyclotomic remainder
modulo a small auxiliary prime's cyclotomic polynomial can recover p+q via
CRT; a single fixed finite catalog gives only bounded residue data.
No native fast evaluator. A shifted hidden-character oracle has a flat
Fourier spectrum of size r^{-1/2}; the large-coefficient positional variant
would give a factor, but no mechanism for deriving it was found.

## 3. Unspecified-root Jacobi bit

Take a unit c with d=c^2+2 a square unit at both primes, and assume the
four roots of s^2-cs-1/2=0 mod N are distinct units. Both primes split -2.
The product of the two local s-roots is -1/2, a square at both primes,
so their local signs agree. J(c)=Jacobi(s,N) is independent of the root.
For any global root t^2=-2 mod N,
    s(c+t)=(s+t/2)^2,
hence J(c)=Jacobi(c+t,N) if c+t is a unit (its norm is d).
This is a decision target with no root supplied; evaluating it by finding
a new independent CRT root would itself split N with constant probability.
Indeed roots of the same quadratic differ in exactly one CRT component
with probability 1/2 under uniform independent sampling.

Reported public randomization: multiply c+omega by U^2 in
(Z/N)[omega], omega^2=-2, write the result A+B omega and normalize c'=A/B
when B is a unit and the original admissibility conditions hold. Then
J(c')=J(c) Jacobi(B,N), by taking the sign at a root of -2.
The exact uniformity proof depends on the admissible gate and should be
reconstructed before claiming a self-reduction for this problem.

Reported class-group lifts, pending recovery/audit of full proofs:
- Lift c+omega to a generator with prime norm ell, appropriate 2-adic
  congruences, and use an unramified dihedral extension of Q(sqrt(-2ell)).
  The Artin sign at a norm-N ideal is J(c).
- An explicitly supplied ambiguous ideal/form (ell,0,2Nh) of discriminant
  -8ellNh has fourth-power membership equivalent to a prescribed J sign
  in a certifying auxiliary-prime branch. Constructing a half yields
  a square root modulo N, hence can factor; decision is the intended target.
- A modified lift ell=3 mod 8, auxiliary h of Jacobi +1 and prescribed
  sign e, gives the reported ordinary-class-number test
  h(-8ellNh)=32 mod 64 iff h is nonsquare at both hidden primes and J(c)=e;
  otherwise 0 mod 64. The claimed sampler has success 1/4 per trial,
  using Dirichlet GRH plus a justified fixed-field prime-sampling estimate.
  This is a reported conditional reduction, not an unconditional
  complexity claim. Recover the exact 2-adic/local sign conventions and
  Redei rank calculation before using it.
- Conversely the reported one-query reduction from ordinary QRP sets
  c=1/2 (known root s=1) and chooses h in the challenge's residue class,
  then enforces the lift sign. The construction, not just the sign
  identity, needs independent audit.

Fixed x=1 (c=-2), r=1 and q=19 mod 24: reported fixed theta form
3u^2+64v^2=N has 8 or 0 signed representations, giving J=- or +.
A CM curve E:y^2=u(u^2+4u+2) gives a_N(E)=-J mod 3 via its 3-division
polynomial. The prior modular-form congruence was
theta(3z)theta(64z) coefficient = sigma_1(n)+a_n(E) mod 3 on n=19 mod 24,
with a finite Sturm check. The exact theta normalization and finite
certificate have not been recovered. Keep this as a target rather than
a certified implementation. Fixed x=1 alone does not give arbitrary
factor bits.

## 4. Exact interval identities and moving-cutoff reduction

For N=rq, r=1 mod 8, q=3 mod 8, define
    K_N(y)=#{1<=a<=y : Jacobi(a,N)=-1}.
The number of units up to N/4 is phi(N)/4 by inclusion-exclusion.
Since chi(2)=chi(-1)=-1, the signed first-quarter sum vanishes:
the even terms in the first half contribute -S(N/4), and the odd terms,
by reflection, contribute S(N/2)-S(N/4). Thus
    K_N(N/4)=phi(N)/8.
The b low bits of this count are efficiently interconvertible with
r mod 2^(b+2) for b>=1: for g(a)=(N+1-a-N/a)/8 on a=1 mod 8,
v2(g(a)-g(a'))=v2(a-a')-2, since N=3 mod 8. Bit lifting is unique.
A known coprime multiplier c=1 mod 8 only scales the quarter count by
phi(c) on the analogous multiplicative family; it is not an extra precision
oracle. This must not be mistaken for computing phi(N)/8 cheaply.

For ANY N=pq=3 mod 4, p<q, set C_N(y)=K_N(y) mod 2. For
0<=y<=floor(sqrt(N)),
    y+C_N(y)+C_N(N-1-y) = floor(y/p) mod 2.
Proof: reflection a->N-a exchanges positive and negative signs; the full
negative count phi(N)/2 is even. The two reflected counts therefore give
unit parity in the prefix. Below sqrt(N), the only nonunits are multiples
of p. Double y (cap at sqrt(N)) until the first nonzero bit; the bracket
has y<2p, so binary search gives p. The proof works for any coloring of
opposite unit pairs with exactly one marked member, not just Jacobi signs.
It does not provide an efficient aggregate.

For odd squarefree M=1 mod 8 the finite class-number dissection is
    h(-8M)=4 S_M(M/8)-2 S_M(M/4)
where S_M(y)=sum_{1<=a<=y} Jacobi(a,M), subject to the usual small-
discriminant unit conventions. For squarefree N=19 mod 24, putting M=3N
and separating the residues mod 3 gives
    h(-24N)=4 sum_{5N/24<a<7N/24} Jacobi(a,N).
This is a useful exact fixed-window target. On the fixed c=-2 family the
earlier class-number argument reported J=+ iff h(-24N)=16 mod 32 on
(q/r)=-1. With K0=K_N(N/4), K=negative signs in the central interval,
this is equivalently K0-K=2 mod 4; the products of the nonzero signs
of both intervals are always +1. The reported auxiliary h=1 mod 24 prime
test uses n=Nh, Jacobi(h,N)=+1, and counts modulo 4 in the half interval
and modulo 8 in the central interval. Its class-number sign and sampling
proof (small h=O(log^2 N) under Dirichlet GRH) must be recovered before
deployment. A factor-free evaluator for variable-endpoint C_N(y) would
already factor mixed semiprimes. Fixed endpoints alone do not imply that
capability. Existing class-group algorithms give subexponential benchmarks;
a constant-factor window shortening is not a breakthrough.

## 5. Research disposition

Keep these as precise tests for genuinely new reciprocity/aggregation
mechanisms, not a request for more equivalent expressions. Reopen conductor
decoding only for a structured algorithm using r|N, and interval counts
only for a bounded-state or shrinking-argument evaluator. The current
primary direction is the [canonical index / ellipse count](shared_collision_nov26.md),
which has a different geometry but a related native count obstruction.
