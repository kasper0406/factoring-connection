# What survives when all quadratic norm directions are combined?

For $N=pq$ with distinct odd primes $p<q<2p$, suppose a short public $x$ gives $A=x^2+1$ and Jacobi symbol $(-A/N)=-1$. Exactly one hidden prime splits in $K=\mathbf Q(\sqrt{-A})$. Call it $r$. The other prime, $d$, is inert. For coprime $k,j$, if $d\mid k^2+Aj^2$ then $j$ must be a unit modulo $d$, giving a square root of $-A$ there, a contradiction. This is the public test that treats the factors differently. The earlier work shows how to use a supplied slope and how to cover nonprincipal binary forms. Here we ask whether multiplying contributions from *all* directions gives an easier computation. [Marked slopes](marked-conic-slopes.md) · [Coverage by reduced forms](pooled-conic-batch.md)

## The full product gives a factor

Put $X=\lfloor\sqrt{2N}\rfloor$. Then $q\le X<2p$, so $p$ and $q$ themselves are the only integers up to $X$ sharing a nontrivial gcd with $N$. Write $\Delta$ for the fundamental discriminant of $K$, and $\chi_\Delta$ for its primitive quadratic Dirichlet character. Given short $x$, trial division supplies the squarefree part of $-A$ and hence $\Delta$. The elementary product
$$
U(X)=\prod_{\substack{1\le n\le X\\\chi_\Delta(n)=1}}n\pmod N
\qquad\text{satisfies}\qquad \gcd(N,U(X))=r .
$$
At $N=77,x=1$, the field has discriminant $-8$. Its positive character classes are $1,3\bmod8$. Below $X=12$ the product runs over $1,3,9,11$ and its gcd with $77$ is $11$.

One can instead multiply the norms of **all integral ideals of the maximal order** with norm up to $X$. The number $a(n)$ of these ideals of norm $n$ is $\sum_{e\mid n}\chi_\Delta(e)$, by the Euler factors for split, inert and ramified primes. In particular $a(r)=2$ and $a(d)=0$. Thus this larger product has exactly the same gcd with $N$, including in fields with nontrivial class group. Its integer value can be rearranged as

$$
Z(X)=\prod_{n\le X}n^{a(n)}
=\prod_{e\le X}\left(e^{\lfloor X/e\rfloor}\lfloor X/e\rfloor!\right)^{\chi_\Delta(e)}.
$$

The last equality is a rational identity with an integer value. Negative exponents can involve nonunits modulo $N$, so direct modular inversion is unjustified. Their gcds agree even though the products can have different residues and multiplicities. Representations by one nonprincipal or principal binary form need not cover $r$; using all ideals covers both degree-one primes over $r$.

## Cost and the limit of a residue-only filter

The mask in $U$ is periodic modulo $|\Delta|$. Its product splits into ordinary arithmetic progression products. On each progression of length $L$, build the polynomial for a block of $\lceil\sqrt L\rceil$ consecutive factors and evaluate it at all block starts with monic multipoint evaluation modulo $N$. A short tail is direct. Summed over the at most $|\Delta|$ classes this costs $\widetilde O(|\Delta|+\sqrt{X|\Delta|})$ modular operations. For fixed $A$ this is the established quarter-power scale. General deterministic factoring attains $N^{1/5+o(1)}$ bit operations. The full ideal product supplies an interpretation of the known factor-separating mask. It supplies no demonstrated improvement on either benchmark. [Pollard–Strassen and fast progression products](https://cs.uwaterloo.ca/~eschost/publications/pollard.pdf) · [Harvey–Hittmeir](https://doi.org/10.1090/mcom/3708) · [Earlier factorial algorithms](../../src/trackC/baselines.py)

More characters alone cannot make a universally valid mask arbitrarily sparse on a complete period. Given any $m$ coprime to $N$, put $c=N\bmod m$. The possible unit residue pairs obey $b=c/a\bmod m$. The involution $a\mapsto c/a$ partitions the units into singletons and pairs. A subset guaranteed to contain an endpoint of *every algebraically compatible pair* has minimum cardinality $(\varphi(m)+f)/2$, where $f$ counts the square roots of $c$ modulo $m$. Each singleton and exactly one endpoint of each two-element orbit attain the bound. In the mixed character example with $m=|\Delta|$, there are no singletons and the positive-character mask attains equality. This is a statement about an entire residue period. Extra information about the actual factors or a bounded interval can change the question. The earlier information-floor result already identified these fibers. [Earlier residue information analysis](../../notes/trackE_analysis.md)

The checker compared the two gcds on 475 mixed balanced fixtures across seven field discriminants. At $N=19\cdot23$ and $A=17$, the principal form misses 23, while two ideals of norm 23 exist and both products recover it.

The relevant next result would be an algorithm to evaluate an input-dependent product or other aggregate with a lower **total** bit cost, or a public way to select a productive subset cheaply. A public character that distinguishes the two primes already exists. [Exact checks](check_collective_norm_product.py)
