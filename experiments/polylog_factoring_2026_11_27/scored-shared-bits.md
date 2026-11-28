# How accurate would an input-only shared-bit test need to be?

For distinct odd primes $p<q<2p$, the fixed singleton binomials yield $A=(p\mathbin{\&}q)\bmod 2^k$. Suppose a cheaper procedure reports a $k$-bit word $M$ which can contain errors. We can rank candidate residues without knowing their locations or their number. Let $w=k-1$, excluding the known units bit, and let $e$ be the number of discrepancies on these positions. On a density-one set of balanced prime pairs the smaller prime's residue is reached after exploring at most
$$
 2^{\gamma k+o(k)}(w+1)\binom{w}{e}
$$
nodes up to subexponential factors, **simultaneously for every report**. Here $\gamma=1-H_2(1/4)=0.1887218755\ldots$. A batch containing this residue yields a verified factor. This is a cost for *using* side information. [Exact fiber theorem](shared-and-fiber-entropy.md) · [Ranker and factoring algorithm](scored_shared_bits.py) · [Reproducible tests](check_scored_shared_bits.py)

## Ranking by exact integer arithmetic

Use blocks of length $\lfloor k^{1/3}\rfloor$ (at least one) above bit zero. At a prefix for the first factor $x\bmod 2^j$, enumerate its possible extensions $x'$ to the next block. For each, compute $y'=N/x'$ modulo the new power of two. Group children by their actual AND block $D$. Let $f(D)$ be the number of children in that group. This enumeration takes $2^{O(k^{1/3})}\operatorname{poly}(\log N)$ per expanded node.

Assign each possible error *word* an initial probability by drawing its weight uniformly from $0,\ldots,w$ and then its support uniformly at that weight. The probability of a particular $l$-bit prefix with $h$ errors is $1/((l+1)\binom lh)$. For example this follows by integrating $t^h(1-t)^{l-h}$ over $0\le t\le1$. Conditional on a node, give a child group the conditional probability of its implied error block $M_{\rm block}\mathbin{\rm XOR}D$. Split that mass equally among the $f(D)$ children. Groups with no children receive no mass. The outgoing mass sums to at most the mass of the node.

The inverse mass of a prefix is the integer
$$
 (l+1)\binom{l}{h}\prod_{\rm blocks} f(D).
$$
Its magnitude never decreases down a path. A priority queue expands least keys first. At each depth at most $K$ nodes have key at most $K$. On the true path the arity product is **exactly the same product for every reported word**, since each $D$ there is the true AND block and the node is the true residue prefix. The local block estimate in the exact-fiber proof bounds it by $2^{\gamma k+o(k)}$ on typical inputs. Shifting the block partition by the fixed units bit still leaves $O(k^{1/3})$ initial and final bits to charge directly. This proves the bound for every error pattern on the same set of true inputs. It also bounds the expense of generating children and comparing keys, since each key has polynomially many bits. The queue enumerates all odd residues if exhausted.

When the exact AND is supplied, the simpler unpruned one-bit tree has the same exponent. The exact checker verifies each compatible population by direct enumeration, and the density-one theorem controls the maximum over all prefix layers. The implementation offers this as the complete-list option. [Direct tree](hh_residue_lists.py) · [Exact whole-fiber count](shared-and-fiber-entropy.md)

## From rank to factoring work

Send successive batches of $1,2,4,\ldots$ leaves to the deterministic residue-list completion. If $R$ bounds the true leaf rank, the sizes tested sum to $O(R)$, and the largest is $O(R)$. The bit cost up to subexponential factors is bounded by
$$
 R+\min\big\{N^{1/5}(R/2^k)^{3/5},\;N^{1/4}R/2^k\big\}
$$
using the corresponding inverse-mass bound for the queue work. Equivalently, one may use the common upper bound $2^{\gamma k+o(k)}(w+1)\binom we$ on rank and queue work throughout. For an error ceiling $\rho w$, $0\le\rho\le1/4$, put $\alpha=\gamma+H_2(\rho)$. Choose $k$ near $\log_2N/4$ when $\alpha\le1/2$, or near $\log_2N/(3+2\alpha)$ otherwise. These choices give the following exponents on almost all distinct prime pairs drawn uniformly from $[U,2U]$. The error ceiling holds on each input of interest, independently of how the mistaken positions were selected. [Completion and balanced sample model](hh-residue-lists.md)

| Maximum errors | Upper exponent of deterministic factoring work after receiving the report |
| ---: | ---: |
| 0% | 0.047181 |
| 5% | 0.118780 |
| 10% | 0.152410 |
| 20% | 0.188881 |
| 25% | 0.200000 |

Entries report $\min(\alpha/4,\alpha/(3+2\alpha))$ to six decimals, with an asymptotically vanishing additional exponent. The unlabelled true AND is zero on about 75% of positions in this sample. Reporting zero everywhere reaches that accuracy already. A public estimator must also pay for generating its report. A cheap reduction in the error ceiling below 25% would improve this deterministic completion bound. A polynomial-time factoring algorithm still requires a stronger source of information.

## A matching limitation for lists using just the low product and report

For $0<\rho<1/4$ choose $\rho'<\rho$ with $\rho'\to\rho$ and $(\rho-\rho')\sqrt w\to\infty$. Draw each reported nonunit bit from a Bernoulli distribution of parameter $(1/4-\rho')/(1-2\rho')$ and independently each error bit with parameter $\rho'$. Take their XOR as the true AND. It is Bernoulli with parameter $1/4$. For a one, take ordered factor bits $11$. For a zero, choose uniformly among $00,01,10$. The resulting pair of odd words is uniform. On the vanishing set where the error count exceeds $\rho w$, use the exact report instead.

Conditional on the initially drawn report, a particular input pair with $h$ errors and $z$ zero AND bits has probability $(\rho')^h(1-\rho')^{w-h}3^{-z}$. With probability $1-o(1)$, $h=\rho'w+o(w)$ and $z=3w/4+o(w)$. For any list rule receiving only the low product and report and keeping at most $L$ first residues per observation, at most $2^w L$ ordered pairs are admitted at each fixed report: there are at most $2^w$ possible odd products, and fixing a product and one unit fixes its partner. Hence lists of size $2^{(\gamma+H_2(\rho)-\epsilon)w}$ for fixed $\epsilon>0$ have vanishing success probability under this reporting scheme. The switch to an exact report on the exceptional set does not change that conclusion. For error ceilings at or above $1/4$ one can take $\rho'\to1/4$ from below to obtain exponent one. For $\rho=0$ the exact-fiber theorem gives the lower bound. This sharpness statement is for the low-residue observation model. The full integer $N$ gives additional arithmetic to a factoring algorithm.

The ranker was checked on 682 small rankings (18,724 exact integer keys), 208 factorizations with altered or unrelated reports, and a separately labelled 69-bit fixture. The older capped-tree method also handles erasures and one-sided guarantees. [Checks](check_scored_shared_bits.py) · [Capped-tree alternative](noisy-singleton-lists.md)
