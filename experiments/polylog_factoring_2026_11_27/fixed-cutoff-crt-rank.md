# CRT rank of the weighted boundary at the fixed ellipse cutoff

At a single prescribed ellipse cutoff, a formula that first represents the entire weighted row boundary as a sum of products of functions of the two prime residues requires exponentially many terms in the input bit length on balanced semiprimes. The bound permits coefficients and local functions that depend on the factors and on the cutoff. Its scope is the full residue table. At a fixed row, the congruence restricts the evaluator to at most four entries, so the relevant open question is whether the boundary and the congruence can be evaluated jointly across rows. [Roadmap](https://github.com/kasper0406/factoring-connection/blob/9514cbfa013cd9ae3f030044db1dc3fb6bb4a9ed/notes/research_roadmap_nov26.md)

## Interval rank

Let $N=pq$ for odd primes $p<q$. Represent $s\bmod N$ by an integer in $[0,N)$, and define a $p\times q$ matrix $M_K$ over $\mathbb F_p$ with entry $[s<K]$ at the pair $(s\bmod p,s\bmod q)$. If $q\le K\le (p-1)q$, then
$$\operatorname{rank}_{\mathbb F_p} M_K\ge\min(p,\lceil q/2\rceil).$$

Write $K=kq+r$ with $0\le r<q$. Representatives in column $j$ are $j+qt$. Among them the interval keeps exactly the first $k+1$ values of $t$ when $j<r$, and the first $k$ otherwise. After reindexing the rows by multiplication by $q^{-1}$ modulo $p$, the column is a circular shift, by $q^{-1}j$, of the indicator of an interval of length $k+1$ or $k$ starting at zero.

One of the two length classes contains at least $\lceil q/2\rceil$ consecutive columns. It therefore contains at least $\min(p,\lceil q/2\rceil)$ distinct shifts. For each occurring length $d$, we have $1\le d<p$. The $p$ shifts of the corresponding indicator form a basis. As polynomials in $\mathbb F_p[X]/(X^p-1)$, they are the images of the standard monomial basis under multiplication by $1+X+\dots+X^{d-1}$. This map is invertible because $X^p-1=(X-1)^p$ and the value of the multiplier at 1 is $d\ne0$ in $\mathbb F_p$.

For the weighted matrix with entries $(a+2s)[s<K]$, multiplying by the weight modulo $p$ scales each row. Exactly one row is zeroed, so the rank drops by at most one. A sum of $R$ terms $f_i(s\bmod p)g_i(s\bmod q)$ has rank at most $R$ over $\mathbb F_p$. Consequently this weighted boundary requires at least $\min(p,\lceil q/2\rceil)-1$ such terms.

## Application at the public cutoff

Use the fixed-cap parameters $A=\lfloor\sqrt{\Delta/3}\rfloor$ and $w^2+\Delta y^2\le4NA$. For $y=1$, write the positive allowed parity coordinates as $w=a+2s$ for $a=1$ when $\Delta$ is odd and $a=2$ otherwise. The upper bound gives exactly
$$K=1+\left\lfloor\frac{\lfloor\sqrt{4NA-\Delta}\rfloor-a}{2}\right\rfloor$$
allowed values if the upper bound is at least $a$. Apply the rank bound whenever $q\le K\le(p-1)q$. On the near-square branch, $A=O(N^{1/4})$ and $K=O(N^{5/8})$, so the upper restriction holds for sufficiently large balanced inputs. If $A\ge3$, then $K\ge q$ for sufficiently large balanced inputs since $q\le\sqrt{2N}$ and $\Delta=O(\sqrt N)$. Even one row then forces $\Omega(\sqrt N)$ terms if its full weighted boundary is materialized in this representation. [Fixed-cap decoder](https://github.com/kasper0406/factoring-connection/blob/9514cbfa013cd9ae3f030044db1dc3fb6bb4a9ed/experiments/ellipse_aggregate/weighted_recurrence_single.md)

## Checks

The check script compares the interval and weighted rank bounds with modular elimination for 412 intervals. For the two historical collision controls and the larger public branch, it uses the fixture primes for the rank certificates and ordinary integer shell enumeration for the trace. The two collision inputs use their original auxiliary branches.

| $N$ | Cap | Eligible first-row coordinates $K$ | Certified weighted rank lower bound | Actual first-row matches | Trace modulo $N$ | Decoded divisor |
|---:|---:|---:|---:|---:|---:|---:|
| 11,009 | 2 | 148 | 69 | 0 | 288 | 109 |
| 28,829 | 2 | 240 | 126 | 0 | 630 | 227 |
| 38,942,974,393 | 358 | 3,733,843 | 153,246 | 1 | 69,529,082 | 254,119 |

The table reports lower bounds for the full weighted boundary, obtained by counting the distinct shifts in the larger length class. The actual matches show the difference in scope between this matrix and the trace computation.

Run with `python check.py` using the Python standard library. [Historical fixtures](https://github.com/kasper0406/factoring-connection/blob/9514cbfa013cd9ae3f030044db1dc3fb6bb4a9ed/experiments/ellipse_row_batches/results.json)

## Roadmap consequence

Drop explicit low-rank separation of the full cutoff mask from the fixed-cap track. A replacement should exploit the congruence throughout the computation. For example, require a proposed batching rule to compute the contribution of a block of rows directly from public data on the modular matches in that block. Measure block construction and empty blocks, and check that the total work over the full cap is polynomial in the bit length. The rank proof gives no complexity bound for that joint contraction or for a factor-blind randomized evaluator of the fixed scalar.
