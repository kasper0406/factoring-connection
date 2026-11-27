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
