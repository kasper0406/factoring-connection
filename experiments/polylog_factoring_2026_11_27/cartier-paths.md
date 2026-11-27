# Coefficient matrices as finite-state paths

The [earlier scalar Hasse coefficient](../../notes/nativemod_results.md) mixes
a prime-index coefficient with its partner's coefficient. In higher genus the
order of the two factors matters. For odd m define B_m(i,j) as the coefficient
of x^(mi−j) in f(x)^((m−1)/2), 1≤i,j≤g, with deg(f)=2g+1. These are
Cartier–Manin/Hasse–Witt coefficient matrices at good prime indices, up to
transpose convention. The composite index denotes only the indicated
coefficient calculation.

For distinct odd primes p,q, expansion in characteristic p gives
B_pq = B_q B_p mod p; expansion in characteristic q gives B_p B_q mod q.
Indeed the first factor f^((p−1)/2) supplies an exponent pk−j and
the pth power of the second supplies p(qi−k). The first exponent bounds
force 1≤k≤g (for p>2g+1), giving precisely matrix multiplication.

For f=1+x^d, d=2g+1, the support of B_r for prime r>d coprime to d
is the partial permutation i ↦ ri mod d on states 1 through g. Exponents
in the required congruence class automatically lie between 0 and
d(r−1)/2. Each nonzero weight is a binomial with top (r−1)/2.
When d<p<q<2p all those tops are smaller than either prime. Each path
has at most one intermediate state, and every surviving weight is a
unit at both primes. Therefore an entry of B_N splits N by gcd exactly
when just one of the two ordered paths survives. For d=5 the residues
2 and 3 give opposite single arrows. At N=91 the two diagonal entries
reduce to 70 and 39 and their gcds are 7 and 13.

[The check](check_cartier_paths.py) evaluated the actual coefficients
and verified the predicted local zero status, entry by entry, for all
110 pairs 13≤p<q<100 with q<2p. Degrees 5, 7, 11 split 17, 26, 52
inputs respectively, and their union split 73. This is a small exact
support check. Computing the individual binomials by their definition
uses indices of order N. The state count alone gives no shortcut for
computing the weights modulo the composite.

Research pointer: seek a coefficient or transfer system with the same
order-sensitive factorization modulo the hidden primes and an explicit
binary-index recurrence of polynomial bit cost over Z/N. Verify both
the recurrence cost (including coefficient growth or nonunit divisions)
and a one-sided matrix entry on unplanted inputs. If the proposed
shortcut requires p,q as transition labels it assumes away the main
difficulty. The curve and Cartier–Manin operator are established tools;
no literature-priority claim is made for this interpretation.

## Scope of the finite-state analogy

The arrows above read hidden prime labels, not digits of the public N.
They describe coefficient support and do not constitute a factoring
automaton. The [earlier coefficient analysis](../../notes/coef_results.md)
rules out a fixed digit automaton for its u₂ coefficient on all integers
and explicitly leaves the semiprime promise open. [Proposition
AUT](../../notes/nativeinv_results.md) also excludes integer-weighted
digit automata for specified multiplicative functions on all integers.
The curve example supplies no exception to either result.

There is a useful arithmetic control for the order identity. In the
cubic quotient T³=T+1, define C_m to have columns the coordinates of
1, T^m and T^(2m). Compute C_N modulo N by binary powering. Modulo p
the map with matrix C_p is Frobenius and is linear on base-field
coefficients, so C_N=C_p C_q. Modulo q the order reverses. This
uses modular arithmetic registers whose number of possible contents
grows with N. On the 200 larger unplanted controls, the [check](check_power_substitution.py) verified both identities and found zero
proper factors from gcds of the entries or determinant of C_N.
Thus fast computation and the order identity alone do not give a
one-sided zero. No speedup follows from the finite-state analogy.

## Minimum degree for a local rank drop in power substitution

For any monic squarefree f of degree d over F_p, let C_N have columns
1,T^N,...,T^((d−1)N) in the power basis modulo f. Over a splitting
field, evaluation at the roots α of f converts this into a Vandermonde
matrix with nodes α^N. Its rank is precisely the number of distinct
nodes. Raising to p is injective in the splitting field, so for N=pq a
collision requires distinct nonzero α,β with α/β of order q. Since
the order of the ratio divides the lcm of the two element orders,
at least one element has order divisible by q. That element lies in
an extension of degree e≤d. Hence q divides p^e−1 and ord_q(p)≤d.

The least degree allowing any such f is exactly max(2,ord_q(p)).
For order t>1, take the minimal polynomial of a primitive q-th root
of unity over F_p: its t distinct conjugates all map to 1. For t=1
take two distinct q-th roots in F_p as the roots of a quadratic.
This sharpness argument uses the hidden primes; it is not a public
construction. Interchange p,q for the other local field.

Repeated factors do not automatically help. For d≤min(p,q), in the
local algebra over a splitting field T=α+ε on each repeated-root
block with ε^d=0. Thus T^N=α^N on that block. Below both order
thresholds the rank simply counts distinct roots, namely
d−deg gcd(f,f′) in these characteristics. A difference of those
degrees is already detectable in polynomial Euclidean division over
Z/N: a nonunit leading coefficient yields a factor, and otherwise
both field computations have the same degree sequence.

The smallest of the two local orders over all 100 unplanted 32-bit
prime controls was 5,987,065, and over all 100 unplanted 48-bit prime
controls it was 133,632,043,168. These are factor-assisted diagnostics
using multiplicative orders, not timings or all-moduli lower bounds.
They exclude every smaller-degree squarefree polynomial for this
full-matrix rank test on the corresponding cohort. Entries, selected
minors, and other matrix constructions have different conditions.

## Do cheap coefficient valuations predict the path type?

At degree five, restrict to balanced p<q<2p with N≡1 mod 5.
The possibilities modulo five are (1,1), (4,4), and the mixed pair
{2,3}. Only the mixed pair has opposite single arrows. Equivalently
this is a fixed quadratic-residuosity question: by quadratic
reciprocity, 5 is a nonsquare at both hidden primes precisely in
the mixed class. It has Jacobi symbol +1 in all three cases.

The two diagonal entries are binomial((N−1)/2, i(N−1)/5) for i=1,2.
Their valuations at any *public* small prime cost O(log N) floor
divisions via Legendre's formula. A [held-out check](check_cartier_valuations.py)
used the twelve valuations from the primes 2,3,5,7,11,13. Fit one
centered least-squares linear score on 2,048 generated semiprimes,
then froze it for a separately seeded test set with 2,048 mixed and
2,048 other inputs. Both sets used disjoint prime-size ranges for p
and q to keep the factor sizes comparable and the larger prime below
twice the smaller. The score classified 994 mixed and 1,061 other
inputs correctly: 2,055/4,096 (50.17%; balanced guessing gives 50%).
This measures one predictor on one distribution. It provides no
evidence of useful predictive gain and gives no impossibility
result for nonlinear statistics or other arithmetic observables.

## Guaranteed witness menu from adjacent binomial carries

Changing the coefficient family removes the need for a favorable degree-five
residue pair. Suppose N=pq with distinct odd primes p<q<2p. For auxiliary
prime d≥5 and j from 1 to g=(d−1)/2, set m_j=⌊jN/d⌋ and take the
gcd of N with binomial(m_j,m_(j−1)). When d<p, m_j<N/2<p².
The low base-p digit of m_j is ⌊p(jq mod d)/d⌋, a strictly increasing
function of jq mod d. Its high digit never decreases with j. By Lucas's
theorem the binomial vanishes modulo p exactly when jq mod d descends
from the preceding residue. At q the corresponding multiplier is p.

For u in [1,d−1] there are exactly ⌊gu/d⌋ descents in the first g steps
of repeatedly adding u modulo d, since each wrap subtracts d. If the
least positive residues of p and q differ by at least three, these two
counts differ (3g>d). Some j therefore gives a proper gcd.

The auxiliary primes can be selected from N alone. Handle q−p=2 by
checking whether N+1 is a square. Otherwise the nonzero integer
∏_{h=−2}^2(q−p+h) has absolute value at most (⌊√N⌋+2)^5.
List consecutive primes starting at 5 until their product exceeds this
bound. One of them divides none of the five terms; for it the two
residues differ by at least three. If the list reaches p, trial division
already splits N. By Chebyshev's primorial estimate the largest d is
O(log N), so there are O(log² N) public binomial indices to test.
This proves witness coverage for the modified family, not for the
original fixed curve matrices. It does not improve the evaluation
cost: some tops have order N, and the argument gives no polylogarithmic
algorithm for their residues modulo N. More elementary binomial
witnesses already exist, for example the central binomial with lower
index ⌊√N/2⌋+1 for separated balanced primes.

The [check](check_adjacent_carries.py) compared the carry predictions
with direct integer binomial gcds on the small cohort of 110 pairs,
up to its first split or trial-division cutoff. On the 200 larger
unplanted controls, factor-assisted path labels predicted a split on
every input by d≤17. Large coefficient residues were not evaluated.

## The coefficient beside the Gaussian CM middle term

For N=pq with p≡1 and q≡3 modulo four, put k=(N−1)/2 and
j=(N−3)/4. In (X³−X)^k = X^k (X²−1)^k, the coefficient at
X^(N−1) is identically zero by parity. The coefficients at X^(N−2)
and X^N both have magnitude C=binom(k,j). This polynomial comes
from the cubic for y²=x³−x; at a composite index the calculation is
ordinary coefficient extraction.

Modulo p, the low digits of k and j are (p−1)/2 and 3(p−1)/4,
respectively. Lucas's theorem gives C≡0 mod p. Modulo q the exact
first-digit decomposition gives

    C ≡ binom((q−1)/2,(q−3)/4) · binom((p−1)/2,(p−1)/4) mod q.

The first factor is a unit. The second is a unit whenever (p−1)/2<q,
which includes all balanced inputs of this residue type. Hence
gcd(C,N)=p on that family. Outside this range, the same identity
characterizes the extra zero by a carry when adding two copies of
(p−1)/4 in base q. At N=91 the unsigned coefficient is 39 modulo 91
and yields 13. Direct integer binomials checked all 65 balanced pairs
of this type with both primes below 100.

This is a d=4 instance of the adjacent-binomial mechanism: the lower
index is floor(N/4) and the upper one is floor(2N/4). The special
feature is its position beside the missing middle term of this CM
cubic. The earlier growing menu covers more inputs. The Lucas proof
uses the hidden primes to explain the coefficient; a fast computation
from N alone is still needed.

## Formal multiplication on the Gaussian CM curve

On E: y²=x³−x at odd N=pq as above, use a local parameter T at
infinity and let F=[N](T) be its formal multiplication series. The
standard height theorem for elliptic formal groups gives order p at
the ordinary prime p≡1 mod 4 and q² at the supersingular prime
q≡3 mod 4. Multiplication by the other prime is invertible on the
formal group locally and preserves the order. Equivalently the
inseparable parts of the respective isogenies factor through the
first and second relative Frobenius. Consequently the entire reduced
series belongs to F_p[[T^p]] and F_q[[T^(q²)]], respectively.

For balanced factors put B=floor(√(2N)). Then p≤B<q², so the gcd
of N and *all* coefficients through degree B is exactly q. In fact
the first nonzero coefficient modulo N has index p, itself a factor.
In contrast the ordinary formal derivative vanishes identically modulo
N. For the jth Hasse derivative, Lucas's theorem on monomials shows
that it too vanishes on the p component unless p divides j, and on
the q component unless q² divides j. In particular no positive
derivative order coprime to N exposes either component this way.

Ordinary truncated series through degree less than p vanish on both
components in the balanced case. More generally, in a unital algebra
of dimension d over either residue field, any nilpotent has index
at most d, by its multiplication operator. Evaluation on a nilpotent
in dimension at most min(p,q²) therefore gives zero on both sides.
The argument concerns explicit nilpotent evaluations. It gives no
lower bound for extracting order or content directly from a succinct
description of the multiplication map.

As an independent check at N=91 use the formal coordinate with
x=T^(−2) and y=−T^(−3)√(1−T⁴); its change from −x/y has unit
linear term and only powers of two in its denominators. The normalized
logarithm has coefficient binom(2j,j)/(4^j(4j+1)) at degree 4j+1,
by integrating the invariant differential. Solving L(F)=91L(T)
coefficient by coefficient over the rationals through degree 49 and
then reducing gave first terms 3T^13 modulo 13 and T^49 modulo 7.
The coefficient at degree 13 is 42 modulo 91, yielding gcd 7.
This cutoff concerns the elliptic multiplication map and is distinct
from the normalized correction series in the Lean comparison.

## Scalar evaluations at multiples of the modulus

In the same coordinate define A(Z)=L(NZ)/N. Its degree-d coefficient
for d=4j+1 is N^(d−1) binom(2j,j)/(4^j d). For every prime ℓ|N
the valuation is at least d−1−v_ℓ(d), since N is odd and squarefree.
Thus every nonlinear coefficient is divisible by N in the localized
coefficient ring, and all terms of degree d>2K+1 vanish modulo N^K.
The reduced rational coefficients have denominators coprime to N
and can be computed with integer arithmetic without its factors.

The relation L([N](T))=NL(T) gives exactly

    [N](Nz)/N = A⁻¹(N A(z)).

Modulo N^K, solve A(w)=y starting with w=y and iterating
w←y−(A(w)−w). The nonlinear part has coefficients divisible by N,
so each iteration gains at least one digit modulo N. Truncation and
K iterations cost polynomial bit time in K and log N. At each hidden
prime A preserves valuation since it is the identity plus nonlinear
terms divisible by that prime. Consequently a unit z gives
v_ℓ([N](Nz))=2 at both primes; more precisely [N](Nz)/N²≡z mod N.
For the latter congruence use A(z)≡z mod N and A(w)≡w mod N³ when
w is a multiple of N.

In the unscaled multiplication series all coefficients are integral
at the odd primes (the formal coordinate change is integral there).
After T=Nz, every term of degree at least p is therefore invisible
modulo N^K for K≤p. In particular a scalar evaluation at that
precision cannot read off the degree-p term from the preceding
section. This observation is specific to this substitution and
precision, and does not preclude arithmetic on the compact map.

For N=91 and K=16, the scaled-log calculation for z=1,2,11,23
verified A(w)=91A(z) mod 91^16 and w/91≡z modulo 91 in all four
cases. Here w=[91](91z)/91.

## Multiplication on actual points in degrees one and two

For distinct odd primes N=pq with max(p,q)<2 min(p,q), multiplication
by N on E: y²=x³−x is a bijection on E(F_ℓ) and E(F_(ℓ²)) for
each ℓ dividing N. This applies in particular to the mixed ordinary
and supersingular case above. It prevents a nonidentity component
from becoming the identity under the final multiplication in either
degree. Intermediate computations and higher extensions are separate.

Write t for the Frobenius trace at ℓ. The degree-two group size is
(ℓ+1−t)(ℓ+1+t). Both factors are multiples of four: the curve and
its quadratic twist each have full rational two-torsion. Both are at
most ℓ+1+2√ℓ by Hasse. For ℓ≥7 this is less than 2ℓ, hence less
than four times either hidden prime by balance. For ℓ=3,5 the Hasse
bound is already less than 12. A positive multiple of four in this
range cannot be divisible by p or q. Thus the degree-two order is
coprime to N; its degree-one subgroup inherits the conclusion.
This argument also applies to any good elliptic curve with full
rational two-torsion at both primes.

Direct prime-field point counts and the trace formula checked all 781
balanced odd prime pairs below 300. Degree three already behaves
differently: at N=299=13·23 the trace at 23 is zero, giving
#E(F_(23³))=12168, divisible by 13. At 13 the trace has absolute
value six and the degree-three group size is coprime to 299
(by the formula ℓ³+1−t³+3ℓt). This supplies an actual one-sided
kernel in degree three, without a general method for finding its
points efficiently. The torsion at 23 in this example has order 13;
it is distinct from the characteristic-23 formal height phenomenon.

## How often small extensions can contain N-torsion

Let p and q be distinct uniformly chosen primes in [B,2B]. Even allowing
the elliptic curve to depend on N=pq, suppose it has good reduction and
a rational point of order two at each hidden prime. The probability
that *any* such curve at either component has nonidentity N-torsion
over *any* extension of degree d≤D is O(D² log B / √B). In particular
this probability is negligible in the input bit length for polylogarithmic
D. This is an existence bound; it grants the curve selector the factors.

At base prime p let t be the Frobenius trace. It is even, and Hasse
gives |t|≤2√p. Define S₀=2, S₁=t, S_d=t S_(d−1)−p S_(d−2).
The extension group size H=p^d+1−S_d is positive. By Cauchy's
theorem a nontrivial N-kernel requires and suffices for p|H or q|H.

* For p|H one needs p|t^d−1, since S_d≡t^d modulo p for d≥1.
  For each even integer t in the Hasse range across the scale and each
  d this is a nonzero integer of magnitude at most (2√(2B))^d+1.
  It has O(d) prime divisors at least B. Summing over O(√B) traces
  and d≤D gives O(√B D²) eligible base primes.
* For q|H fix p and allow every trace in its Hasse range. For each
  trace that occurs for a curve and each d, the bound
  H≤(p^(d/2)+1)² again permits O(d) prime divisors at least B.
  This gives O(√B D²) eligible partner primes per p, even for
  adaptive curve choice. For the fixed Gaussian CM curve there is
  just one trace per p, improving this part alone to O(D²) partners.

Divide by the Θ(B/log B) available primes and allow both component
orientations. Excluding equal primes changes only the constant.
The even-trace hypothesis matters in the first case: trace one would
make t^d−1 zero. The result concerns actual points in finite residue
fields under multiplication by N; formal coefficients, nonreduced
algebras and other multiplication scalars are outside the argument.

For the fixed curve y²=x³−x, exact prime-field counts and the trace
recurrence gave the following numbers of pairs admitting any nonzero
kernel (including a kernel on both sides):

| B | Pairs | Degrees 1–3 | Degrees 1–8 |
|---|---:|---:|---:|
| 128 | 253 | 0 | 31 |
| 512 | 2775 | 5 | 31 |
| 2048 | 32385 | 29 | 379 |

## The trace-one exception already has an operational CM splitter

The [Cheng implementation](../../src/nonsmooth/cheng.py) and
[earlier nonsmooth audit](../../notes/nonsmooth_results.md) already
cover this exception. When 4p−1=Dv², a suitable CM twist of
discriminant −D can have exactly p points at the hidden prime p.
Multiplication by N then kills that component. The earlier audit
reports 180/180 special inputs split and 0/180 random controls.
As a replay, cheng_factor with seed 20261127 and d_list=(3,) returned
19 from 437 in one try, 37 from 1517 in eight, and 61 from 4331 in
three. Only N, the seed and public configuration entered the solver.

One can bound coverage by the *number* of discriminants even when
their sizes are unrestricted. At scale p in [B,2B], take a shared
menu of H distinct positive integers D. For each D there are at most
floor(√(8B/D)) positive v, each determining at most one p. Ordering
the menu as D₁<⋯<D_H gives D_i≥i, so the total number of covered
primes is O(√(BH)). Dividing by Θ(B/log B) primes, and allowing
either hidden factor to be covered, gives fraction
O(log B · √(H/B)) on balanced pairs at this scale. This bound even
allows perfect selection from the menu after seeing N. It applies to
the trace-one discriminant condition, independently of the cost of
constructing the curves. A discriminant computed from N outside a
small shared menu or a different trace mechanism is outside its scope.
