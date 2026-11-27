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
