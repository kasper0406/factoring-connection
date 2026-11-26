# Computing the quartic parity with modular symbols at level 59

For a promised $N=pq=t^2+3t+17$ on the canonical discriminant $-59$ branch, the parity of the roots of $X^4-2tX^2-8X+t-5$ equals a sum of ordinary continued-fraction paths. The level and the linear algebra are fixed. Each path costs $O(\log N)$ Euclidean steps; our exact evaluator visits all $N$ endpoints. A shortcut using seven paths, including two whose offsets grow with $t$, takes polylogarithmic time but gives exactly the same output on examples with zero, one and sixteen quartic roots. [Quartic and input promise](cubic59.md)

The repository already computes fixed-level modular symbols by continued fractions, lists the cost of enumerating the full Hecke correspondence, and tested a resultant on a related norm-one quartic. We reused its projective-line conventions and independently compared its Heilbronn matrices with the endpoint implementation. The branch identity allows us to remove the unknown-divisor and coprimality filters from the sum. [Fixed-level endpoint method](../../src/newfamily/modsym.py) · [Composite-index Hecke work](../../notes/nativemod_results.md) · [Earlier quartic resultant](../../notes/buildbit_results.md)

## A four-dimensional quotient that reads the right bit

Write $e=\{\infty,0\}$ in the relative weight-two modular-symbol space for $\Gamma_0(59)$ over $\mathbb F_2$. A path between rational endpoints is represented by its unimodular continued-fraction edges. Its 60 possible edge labels are the bottom rows in $\mathbb P^1(\mathbb F_{59})$. The two-term and three-term Manin relations give dimension 11. At this level there are no elliptic points of orders two or three in the projective group. The integral presentation therefore introduces no elliptic torsion on passage to characteristic two. [Manin relations and Hecke operators, chapter 3](https://www.williamstein.org/books/modform/modform/weight_two.html)

For each good odd prime $r$, let $\epsilon_r=1$ if $Y^3+2Y+1$ is irreducible modulo $r$, and 0 otherwise. Form the quotient $V$ by the images of $T_r-\epsilon_r$ for $r=3,5,7,11,13$. Finite row reduction gives $\dim V=4$, $e\ne0$, and $T_2^2=1$ on $V$. The relations are Hecke-stable by commutativity. The checks verify the matrix identities, including preservation of the relation space by $T_2$. The following argument explains why the five small primes suffice for the good odd primes of any size.

1. **The quotient comes from cuspidal homology.** There are two cusps. On their degree-zero boundary $T_5$ acts by $1+5$, which vanishes modulo two. Thus $T_5$ maps all relative symbols into the cuspidal subspace. Since $\epsilon_5=1$, this subspace surjects onto $V$. The boundary group is torsion-free. The integral cuspidal Hecke algebra acts on $V$ by compatibility of the integration pairing on absolute homology.
2. **Two explicit weight-two forms certify the odd eigenvalues.** For $K=\mathbb Q(\sqrt{-59})$ let $\Theta_0,\Theta_1$ be the full lattice theta series of its principal class and either nonprincipal class. The two nonprincipal classes have equal theta series. The normalized dihedral eigenform $g=(\Theta_0-\Theta_1)/2$ has weight one and quadratic character $\chi_{-59}$. This follows by grouping the ideal-character series by the three classes: the two values on the nonprincipal classes sum to $-1$, and each ideal corresponds to two vectors. Its good prime coefficient is the trace of the two-dimensional dihedral representation, reducing to $\epsilon_r$. The prime 2 is inert in $K$ and $a_2(g)=0$. Both $g\Theta_0$ and $g^2$ are integral classical weight-two cusp forms with trivial character at level 59. Their reductions are $g(q)$ and $g(q^2)$ because $\Theta_0\equiv1\bmod2$. They are independent by their first exponents. Every good odd $T_r$ acts on both as $\epsilon_r$ by the Fourier coefficient recurrence. In weight two modulo two $T_2$ selects the even coefficients, so it swaps them. Here the weight-one recurrence gives $U_2g=g(q^2)$ modulo two. [Ideal-character theta series](https://www.math.mcgill.ca/darmon/pub/Articles/Research/79.DHRV/michigan.pdf)
3. **Transfer the calculation to $V$.** The integral cuspidal Hecke algebra in weight two and level 59 is generated additively by $T_n$ through $n=\lfloor 2(59+1)/12\rfloor=10$. The product and prime-power recurrences then give algebra generators $1,T_2,T_3,T_5,T_7$. Express any good odd $T_r$ in these generators, reduce modulo two, and substitute the calibrated scalars for the three odd generators. The resulting polynomial in $T_2$ acts as $\epsilon_r$ on the two forms above. The swap matrix has minimal polynomial $X^2-1$ in characteristic two. The polynomial is therefore $\epsilon_r$ modulo $X^2-1$, so it gives $\epsilon_r$ on $V$ as well. [Agashe–Stein, integral Hecke generation bound](https://www.williamstein.org/papers/generating_hecke/appdx.pdf)

## An exact sum with no factorization

Put $P(b,d)=\{\infty,b/d\}$ and $A_d=\sum_{0\le b<d}P(b,d)$ in $V$. The determinant-$n$ Hecke representatives send $e$ to these paths, giving $T_ne=\sum_{d\mid n}A_d$ for $n$ coprime to 59. In particular $A_1=e$. On the branch the primes have the same cubic type, so $\epsilon_p=\epsilon_q=B$, the quartic root parity. Consequently $A_p=A_q$ and $T_{pq}e=T_pT_qe=Be$. This gives the public identity

$$
\sum_{b=0}^{N-1}\{\infty,b/N\}=(1+B)e .
$$

All numerators enter including nonunits. Accumulate by XOR in $V$ and return 1 for a zero sum, 0 for $e$. Setup constructs fixed matrices of dimension 60 once. Each integer division uses operands of $O(\log N)$ bits. The complete scan has $N$ terms and at most $O(N\log N)$ Euclidean edges. The comparison below counts all edges. This benchmark trades the quartic for uniform small-state summands; it costs more than the earlier direct norm-count verifier for fixed discriminant.

| $N$ | Quartic roots from independent local checks | Returned bit from $N$ alone | Euclidean edges |
|---:|---:|---:|---:|
| 31,877 | 1 | 1 | 323,567 |
| 252,521 | 0 | 0 | 3,011,625 |
| 29,014,397 | 16 | 0 | 462,654,503 |

For the first example we also summed literal uncompressed edge words and reduced only at the end. For six small primes the 360 Hecke columns agree with the repository's independent Heilbronn construction in the **full** relative space. Every good odd prime from 17 through 97 passed a separate scalar-action check.

## Why seven fast endpoints miss the bit

The branch gives the public slope $t/N$. Consider the paths to $(t+s)/N$ for fixed $s$. Set $m_s=s^2-3s+17$. The first division satisfies

$$N=(t+s)(t+3-s)+m_s.$$

When $t+s>m_s$ and $0<t+s<N$, the next quotient is $\lfloor(t+s)/m_s\rfloor$ and the rest of the Euclidean algorithm starts with the remainder $(t+s)\bmod m_s$ (or terminates). Thus the **entire list of projective-line edge labels**, including its length, depends only on $t\bmod 59m_s$. The same proof works for a different fixed prime level by replacing 59. The theorem concerns the reduced edge labels; the full integer quotient still contains information about $t$.

For $s=-2,-1,0,1,2$ the five paths have joint period dividing 947,835. We also evaluated the growing offsets $s=\pm\lfloor\sqrt t\rfloor$, where the fixed-period argument does not apply. In the progression $t=5385+947835k$ we sought balanced inputs whose seven values in $V$ match. The three rows below have exactly the same four-bit value at each corresponding endpoint. The five fixed paths additionally match edge for edge before projection.

| $t$ | $N$ | Quartic roots | Target bit |
|---:|---:|---:|---:|
| 5,385 | 29,014,397 | 16 | 0 |
| 80,573,563,065 | 6,492,099,065,231,252,883,437 | 1 | 1 |
| 80,639,911,515 | 6,502,795,329,388,949,329,787 | 0 | 0 |

The witnesses were selected using control factorizations. The fast endpoint calculations use only $N$. The collision excludes every postprocessing function of those seven $V$ values alone. It leaves open other endpoints, growing-level symbols, unreduced words at the two growing offsets and algorithms using $N$ in other ways.

**Next blocker.** Sum the interior endpoints collectively while charging the construction of every group. A concrete starting point is a recurrence for the parity of visits to a given projective edge across all continued fractions $b/N$, with $N=t^2+3t+17$. It would have to avoid enumerating either the $N$ numerators or a power-sized set of intermediate states. Use the matched triple before drawing conclusions from another small selection of endpoints. Independent limitations on how these bits transfer across discriminants are recorded in [the growing-discriminant experiment](growing59.md).

From the repository root, with python-flint, NumPy and SymPy installed:

    python experiments/ellipse_disk_mod8/manin59.py 31877 --complete
    python experiments/ellipse_disk_mod8/manin59.py 31877
    python experiments/ellipse_disk_mod8/manin59_check.py --output rerun.json

[Public evaluator](manin59.py) · [Independent checks](manin59_check.py) · [Recorded results](manin59_results.json)
