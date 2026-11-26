# Multiplication-network and shared collision search (2026-11-26)

Status: AI-assisted research consolidation, based on the conversation following commit
2db67a9190340630afaef7a023a9bff981bdcce9. **No new end-to-end factoring exponent is established.**
Elementary identities below have been rederived; full asymptotic claims from the conversation
are recorded as pending independent proof audit when the original scratch proof is unavailable.
Original late-stage scratch scripts/data could not be recovered in the consolidation session.
See [validation and provenance](../experiments/nov26_consolidation/README.md) and the
[roadmap](research_roadmap_nov26.md). Do not interpret reported toy checks as complexity proofs.

## 1. Multiplication counting modulo 3

Rudich–Rousseau, [arXiv:2610.00600](https://arxiv.org/abs/2610.00600), give a recursive
multiplication-network count with deterministic subexponential (N^{o(1)}) runtime.
This is not by itself a polynomial-time algorithm or a comparison with L[1/3].

For distinct odd prime factors, remove the known assignments 1 and N from an ordered
divisor-prefix count. There are at most two other assignments. Thus the exact normalized
count modulo 3 is sufficient; a count with a known factor 2^E can be normalized over F_3.
Correction: the published implementation can duplicate branches with multiplicity divisible
by 3. It cannot simply be run mod 3 unchanged. Enumerate each genuine branch exactly once,
or otherwise prove its multiplicity is invertible and known. The expense of intermediate
tables survives. A fresh audit of the modified contraction and its complexity is required.

### Cut-rank barrier (proved in the specified flattening)

Let m=2^k, k>=3. Write operands d=u+ms, e=v+mt, with each half in [0,m).
For d e = N (mod m^2), N odd, the low halves are odd and v=N/u (mod m).
The allowed high halves form the affine line
t=b_u-(N/u^2)s (mod m), where b_u=(N-uv)/(mu) interpreted modulo m (divide by m first).
There are m/8 distinct slopes: the squares among odd units modulo 2^k are exactly
the residues 1 mod 8. Over an extension of F_3 with a primitive m-th root, the Fourier
character perpendicular to one line has nonzero sum on that line and zero sum on all
lines of other slopes. Therefore the low/high flattening has rank at least m/8 over F_3.
This is a barrier for low-rank representations across this cut, not arbitrary symbolic
representations, elimination orders, or tensor networks.

### Masked product counts (local construction)

For a fixed integer pair of positive weights a,b, set H=ap+bq and Z=ap-bq. Then
Z^2=H^2-4ab pq; p=(H+Z)/(2a), q=(H-Z)/(2b).
For a specified H and interval of possible pq, the constraint is an interval (or two
intervals) for Z, plus linear integrality and operand bounds. Requiring pq=N mod 2^k
gives Z^2=H^2-4abN mod (4ab 2^k). In the dyadic coprime cases studied, solving the
2-adic square congruence yields a bounded number of progressions (including compressed
families for nonunits); intersect with the interval and linear constraints.
For fixed operands or fixed levels the count is polynomial in the bit length.
Account for odd factors of a,b separately; the 2-adic argument alone does not solve
an arbitrary composite-modulus root problem.

Reported local algorithms: near sqrt(N), a block of O(N^{1/4}) operand positions with
O(N^{1/2}) product width; arbitrary specified comparable center, O(N^{1/6}) positions
with O(N^{1/3}) product width, using a small rational approximation to the slope.
11,000 two-operand toy cases reportedly agreed with enumeration.
Reported separate-level cover: Theta(N^{1/3}) for full coverage of [sqrt(N)/2,sqrt(N)]
by positive weighted levels, under the particular listing model. Reported dyadic
filter bound: every 128 consecutive levels contains one compatible with low product
bits alone in the specified coprime-weight regime. Original proofs missing: do not
cite the last two estimates as general lower bounds until reconstructed.
Next: batch across levels or change the cut; fixed-level compression is not enough.

## 2. Structured square-weight collision evaluation

Existing baseline: [dcltheorem](dcltheorem_results.md), Theorem D' and Corollary D'',
internally reviewed, give N^{1/6} polylog(N) deterministic time on all but an inverse-
polylogarithmic fraction of bounded-ratio prime pairs (may return FAIL).
It remains exponential in log N and is not an improvement over heuristic GNFS.
See also [multimod](multimod_results.md) on the collision-stage bottleneck.

For public units x,y,z modulo N, consider g(u,v)=x^(u^2)y^(v^2)z^(uv).
A target T is hit locally if T=g(u,v) modulo a hidden prime. For a list of m unit
targets let Q(t)=prod_j(1-t/T_j). Handle nonunits by gcd. In a block starting at v0,
g(u,v0+i)=A B^i y^(i^2), A=g(u,v0), B=y^(2v0)z^u.
Build a shared polynomial
    F(A,B)=prod_{0<=i<L} Q(A B^i y^(i^2)).
Evaluate it at the block anchors. The original construction invokes fast
multivariate multipoint evaluation over Z/NZ ([Bhargava et al., arXiv:2205.00342](https://arxiv.org/abs/2205.00342));
the claimed grid bound is G^{3/4+o(1)} polylog(N) for m=1 and
G^{3/4+o(1)} sqrt(m) polylog(N) for m<=sqrt(G). Treat these bounds as pending
audit of the actual degree/representation, ring hypotheses, and evaluation cost.
The block identity itself is elementary. Reported 700 direct checks agreed.
Retain row-block products separately: the full grid product may vanish at both
factors, whereas a block gcd can be proper. Example from the conversation:
N=3053, a block value 2580 gives gcd 43; the target was factor-backed for this check.

A short-list evaluator alone does not give a faster factorer. Need a list of
<=G^{1/2-epsilon} public targets that retains a one-sided hit often enough, and charge
list construction, all failures, and extraction. The reported uniform-shortlist
bound is O(m/Y) at the representative band Y~N^{1/6} under the stated sampling
hypotheses; full proof not recovered. Do not count factor-backed target selection
as an algorithm.

## 3. Canonical one-index target

Suppose N=pq, p<q<=2p, p,q~P. Let ell be a prime not dividing 2N for which
4N has a square root modulo ell^2. Choose a root branch and its largest positive
representative c<2sqrt(N), and let Delta=(4N-c^2)/ell^2>0. Assume gcd(N,ell Delta)=1,
and work with the quadratic order of discriminant -Delta. Put
    alpha=(c+ell sqrt(-Delta))/2.
This is integral of norm N. The invertible prime ideals (p,alpha),(q,alpha)
are coprime, and their product is (alpha); their classes C,C^{-1} have equal
minimum A_C. All class conventions below are up to inversion.

If gcd(u,v)=1 and 2pv=c u (mod ell), the integer
    j=(p v^2+q u^2-cuv)/ell^2
is represented primitively by the class C (or its inverse). Indeed write
v=k u+ell w, k=c/(2p) mod ell. The restricted form divided by ell^2 is
[(p k^2+q-ck)/ell^2, (2pk-c)/ell, p], of discriminant -Delta.
Primitivity in (u,v) implies primitivity in (u,w). This form represents p
and corresponds to one of the two prime ideals over p; the choice of sign does
not affect the minimum.

**Elementary lemma.** If 4j^2<Delta, then j=A_C. For two independent integral
vectors in the same positive definite binary form, the product of their values
is at least Delta/4 (the determinant of their Gram matrix is Delta/4 times the
square of an integer determinant). A vector of no greater value must therefore
be collinear with the primitive vector, hence a multiple of it.

**Extraction, proved.** Given a candidate A_C, for each grid weight set
H=ell^2 A_C+cuv. At a true hit H=p v^2+q u^2 and
    H^2-4Nu^2v^2=(p v^2-q u^2)^2.
Test whether the discriminant is a nonnegative integer square and test the
divisibility of (H +/- sqrt(discriminant))/(2v^2). No candidate is trusted
unless it is an exact proper divisor of N.

**Reported band-coverage theorem, pending full audit:** K~log^2(P);
primitive u in [K,2K), v in [2K,3K), prime ell in [Y,2Y].
At Y~P^{1/4} log^4(P), all but O(log^{-2} P) of bounded-ratio prime pairs
have a canonical hit. At Y<=P^{1/4}/log^4(P), at most O(1/log P) have such a
hit (same grid). The upper proof was described as the large-sieve coverage
lemma in dcltheorem plus a Gauss-sum estimate ruling out exceptionally small
Delta. The lower proof uses the necessary slope proximity
|v-u sqrt(q/p)| << Y/P^{1/4}. Reconstruct and audit both, including the
number of distinct square-weight residues, c choice, and distribution of
exceptional small discriminants; the original proof file is unavailable.

**Conditional payoff only:** if the upper coverage theorem is valid and
A_C is evaluable from public (N,ell,c) in polylog time, scanning the band and
trying all polylog weights costs N^{1/8} polylog(N) on almost all balanced
prime pairs. No evaluator for A_C and no end-to-end speedup have been found.
Reported experiment: 192 semiprimes, 51 with grid hits, canonical minimum
recovered a factor on 43 of those 51, accounting for 52/64 locked weights.
These minima were factor-backed.

A primitive representation j cannot be divisible by an odd prime l coprime
to Delta for which (-Delta/l)=-1: the reduced binary form is anisotropic mod l,
forcing both coordinates to be divisible by l. Thus inert auxiliary primes
can exclude candidate indices, but no proof of adequate shortlist sparsity.

## 4. Public ellipse count and squared-class decoding

Let Q0 be the known principal form of discriminant -Delta<-4, and a_E(j) count
integral invertible ideals of norm j in class E. To avoid conductor issues,
use primitive invertible ideals and the appropriate form/ideal correspondence
(or initially restrict to fundamental discriminants); explicitly audit the
nonfundamental case. If gcd(j,N)=1, factorization at the unramified primes
above p,q gives
    a_1(Nj)=2 a_1(j)+2 a_{C^2}(j).
The four ideal choices have classes 1,1,C^2,C^{-2}, and conjugation gives
equal counts for the last two. With only the units +/-1, the norm representation
count up to sign is the corresponding principal ideal count.

Define U_N(B)= (1/2) #{(x,y) in Z^2: 0<Q0(x,y)<=NB, N | Q0(x,y)} mod 8.
A query returns one count modulo 16 before division by 2. For B<min(p,q) and
in the range where the principal form has only axis vectors,
    U_N(B) = 2 floor(sqrt(B)) + 2 V(B) (mod 8),
where V(B) counts nonzero vectors up to sign in the C^2 form through B.
More generally the principal term is known and can be evaluated by direct
form-specific means; charge it. For the promised large-ell balanced branches,
Delta is small enough that the decoding cutoffs are below p,q; verify the
precise endpoint and principal-term hypotheses before invoking the decoder.

For a reduced form [a,b,d] with |b|<=a<=d, the first doubling endpoint
reaching a lies below 2a, where there are at most three nonzero vector pairs.
The count mod 4 is therefore exact there. Binary search gives a.
Subtract the floor(sqrt(B/a)) multiples of (1,0); at B<=floor((Delta+a^2)/(4a))
there are at most two further vector pairs. Indeed d<=that bound, y=+/-1
is necessary for any such nonaxis vector (complete the square), and the
allowed interval of x has length <=1 after centering, with endpoint ties
handled by the count. The first additional value gives d, and |b|=sqrt(4ad-Delta).
If the principal term is cheap and the query conditions above hold,
O(log N) counts recover the reduced class C^2 up to inversion. The direct
enumeration test reported 30,719/30,719 reconstructed reduced forms.
The count is not currently known to be cheap.

**Class-square ambiguity:** the minimum of C does not in general follow from C^2.
Conversation counterexample: Delta=20, ell=11, (N,c)=(11009,204),(28829,336);
both squared classes principal, minima 1 and 2. Public N,c differ, so this
does not exclude algorithms using them. A known odd multiple e of the class
order would give C=(C^2)^((e+1)/2); alternatively enumerate square roots of the
class and verify candidate minima. Efficient lifting without factoring the
discriminant, and with manageable total number of lifts, is open. Known
class-group square-root algorithms require extra data; see
[genusleak](genusleak_results.md).

## 5. Next steps

1. Independently reproduce the band coverage proof and the multipoint exponent.
2. Build a factor-free evaluator for the modular ellipse count; the known
   principal norm-N representation supplies only the two same-orientation
   ideals. The cross-oriented CRT root would directly split N by gcd.
3. Analyze and implement class-square lifting for the actual discriminants
   in the band (including nonfundamental and small-discriminant exceptions).
4. Charge end-to-end branch generation, class lifting and failures; compare
   with D'' and with GNFS, not just separate-level enumeration.

No factorization in P or sub-L[1/3] factorization follows from these results.
