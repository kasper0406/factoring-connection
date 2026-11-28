# Predicting a factor's residues from a balanced semiprime

Let $p,q$ be two distinct primes sampled uniformly without replacement from $[R,2R)$, let $N=pq$, and put $X=R^2$. Fix $A>1$ and $0<t<1/2$. At each prime $\ell\le B=(\log R)^A$, a public rule sees $N\bmod\ell$ and predicts a single residue $a_\ell$ for $p\bmod\ell$. Count an agreement with weight $\log\ell$. Even the best such rules reach total weight $t\log N$ with probability at most $X^{-t/A+o(1)}$. The exponent is sharp for the agreement event. Independent randomized public rules attain it for every fixed balanced input, over their own randomness. The comparison concerns the information supplied to the weighted CRT decoder. It grants the scorer the factor to check which coordinates agreed. [Order-hint decoder and determinant](factoring-with-an-order-hint.md) · [Runnable checks](check_local_order_scores.py)

## Counting local agreements

For each $\ell$, allow any function $a_\ell:(\mathbb Z/\ell)^\times\to\mathbb Z/\ell$. The functions can vary with $R$. A portfolio of $X^{o(1)}$ complete rule tables can also select a table after inspecting the whole input. Discarding coordinates after inspection only decreases the score. Conditioning on independent random choices gives the same upper bound. Either factor may be graded, at a cost of a factor two.

For a squarefree product $d$ of tested primes, fix the product residue $u=pq\bmod d$. Agreement throughout $d$ prescribes one residue for $p\bmod d$ by CRT. If this is a unit, it fixes $q\bmod d$ as well. Thus at most $\varphi(d)$ of the ordered residue pairs modulo $d$ qualify. Every residue class has at most $\lceil R/d\rceil$ representatives in the interval. For $d\le R$ the number of qualifying ordered integer pairs is at most $4R^2/d$. The prime number theorem supplies $\Theta(R^2/\log^2 R)$ ordered distinct prime pairs. This argument uses no distribution of primes in progressions of growing modulus.

Let $D$ be the product of the agreeing primes. If $D\ge N^t$, multiply its primes until the product $d$ first reaches $X^t$. Then $X^t\le d<BX^t<R$ for all sufficiently large $R$. A union bound and Rankin's inequality give, for $0<\sigma<1/A$,
$$
 \Pr(D\ge N^t)
 \ll \log^2 R \sum_{\substack{d\ge X^t\\d\mid\prod_{\ell\le B}\ell}}\frac1d
 \le \log^2 R\,X^{-t\sigma}\prod_{\ell\le B}(1+\ell^{\sigma-1}).
$$
The logarithm of the product is $O_\sigma(B^\sigma)=o(\log R)$, even on summing over all integers through $B$. Letting $\sigma$ approach $1/A$ gives the asserted exponent. The same proof applies to a threshold $N^{t-o(1)}$. At $t=1/4$, it bounds the probability of attaining the critical auxiliary score on full balanced intervals. For a prime coordinate the contribution $(\log\gcd(\ell,p-a_\ell))^2/\log\ell$ is exactly its agreement weight.

## A public construction attaining the exponent

Use every prime in $[B/2,B]$. Draw independent uniform units $c_\ell\bmod\ell$ and predict $a_\ell(u)=c_\ell u^2\bmod\ell$. For every fixed $N=pq$ in the cohort and large enough $R$, the indicators of agreement with $p$ are independent and have respective probabilities $1/(\ell-1)$. Let $m$ be the number of band primes and set $k=\lceil t\log N/\log(B/2)\rceil$. Every $k$-subset has product at least $N^t$.

The probability of exactly $k$ agreements is the coefficient of $z^k$ in $\prod_\ell((\ell-2)+z)/(\ell-1)$. The events specifying the exact subset are disjoint. Their combined probability is at least
$$
 \binom{m}{k} B^{-k}\prod_\ell(1-1/(\ell-1)).
$$
The last product tends to one. The prime number theorem gives $m\asymp B/\log B$, whereas $k=(t+o(1))\log N/\log B$, so $m/k\to\infty$. Using $\binom{m}{k}\ge(m/k)^k$ gives log probability at least $-k(\log k+\log\log B+O(1))=-(t/A+o(1))\log N$. For the matching upper bound on each fixed input, $\mathbb E[D^\sigma]=\prod_\ell(1+(\ell^\sigma-1)/(\ell-1))=\exp(o(\log N))$ when $\sigma<1/A$, and Markov's inequality applies. Identifying the agreeing coordinates from $N$ remains a separate algorithmic step.

## Which local data enter the count

Two variations follow from the same residue-fiber calculation.

- For $t(1+1/A)<1/2$, including $t=1/4$, take any predeclared pool of at most $(\log R)^A$ distinct primes up to $X^{1-t/A}$. A pool prime equals a hidden prime with probability at most $X^{-1/2+o(1)}$ over the whole pool. Away from these inputs, any pool prime $\ell\ge X^{t/A}$ has match probability $O(\log^2R(1/\ell+1/R+\ell/R^2))=X^{-t/A+o(1)}$. Take the union over these primes. Among the others the first-crossing product is below $X^{t+t/A}<R$. The Rankin product has log at most $O_\sigma((\log R)^{A\sigma})$, by bounding the $i$-th prime below by $i+1$. This proves the same exponent for pools containing large, widely separated primes.
- For tested primes in $[B/2,B]$, allow each anchor to read $N$ modulo its prime and one further preassigned prime at most $B$. At $t=1/4$ let $d'$ be the least common multiple of all moduli consulted by a first-crossing $d$. Then $d'\le d B^{\omega(d)}\le d^2 X^{o(1)}\le R^{1+o(1)}$, using $\omega(d)\le\log d/\log(B/2)$. Fixing the product modulo $d'$ fixes the predicted first factor modulo $d$, with at most $\varphi(d')/\varphi(d)$ unit lifts. Hence at most $\varphi(d')^2/\varphi(d)$ ordered unit pairs modulo $d'$ qualify. Counting interval representatives now gives match probability $X^{o(1)}/d$, since $d/\varphi(d)=X^{o(1)}$. Summing yields the original exponent even with this extra input to every table.

The earlier CRT audit studied agreement of integers in smooth-number windows, and the earlier exponent-list bound tested bit-length-only exponents against $p\pm1$. Here the observation is the actual product at each permitted modulus, and the target is the residue of one of its two prime factors. [CRT audit, §1](../../notes/crtdecode_results.md) · [Johnson certificate analysis, §5.2](../../notes/sihard_results.md) · [Previous exponent-list bound](quadratic-chirp.md#bit-length-only-exponent-lists-miss-almost-all-balanced-prime-pairs)

## Whole-input anchors near the square root

A second construction uses $s=\lfloor\sqrt N\rfloor$. Fix any menu $J=J(R)$ of at most $(\log R)^A$ integer shifts, and form $a_j=s+j$ and $e_j=|N-a_j^2|$. An input-dependent choice may retain any pairwise coprime collection of these full moduli with $e_j\ge2$. For a target factor $r\in\{p,q\}$ put $g_j=\gcd(e_j,r-a_j)$. The probability of finding a proper factor by $\gcd(N,a_j)$, or of attaining
$$
 \sum_j\frac{(\log g_j)^2}{\log e_j}
 \ \ge\ (1/4-o(1))\log N
$$
with unit anchors, is at most $X^{-1/4+o(1)}$. The statement concerns this specific family and score. In particular its moduli are the full Fermat remainders.

For the preliminary gcd, fix $r$ and $j$. The value $\lfloor\sqrt{rq}\rfloor$ lies in an interval of length $O(R)$, so $r\mid s+j$ selects $O(1)$ possible values. The derivative with respect to $q$ is bounded above and below by positive constants on the balanced range. Each selected value therefore has $O(1)$ integer preimages. The total probability at either prime, over the whole menu, is $X^{-1/2+o(1)}$. Also, integers within $X^{1/4}$ of any square in the input range number $O(RX^{1/4})$. After normalization they have probability $X^{-1/4+o(1)}$. On the remaining inputs all the $e_j$ are at least $X^{1/4}$.

For a unit anchor $a$, reduction modulo $r-a$ gives
$$
 \gcd(|N-a^2|,r-a)=\gcd(q'-r,r-a),\qquad q'=N/r.
$$
Indeed $N-a^2\equiv a(q'-a)$, and $\gcd(a,r-a)=1$. Thus all $g_j$ divide $V=|q'-r|<R$. They are pairwise coprime on a retained collection. Set $\gamma_j=\log g_j/\log X$. Their sum is at most $1/2$. Since $\log e_j\ge(\log X)/4$, a qualifying score forces $\sum_j\gamma_j^2\ge1/16-o(1)$. For fixed $\eta>0$, coordinates smaller than $\eta$ contribute at most $\eta/2$ to this last sum. The remaining set $S$ has at most $1/(2\eta)$ members. Consequently
$$
 d=\prod_{j\in S}g_j\ \ge\ X^{\sqrt{1/16-\eta/2-o(1)}},\qquad d\mid V ,
$$
for small enough $\eta$. The number of possible sets $S$ is $X^{o(1)}$ for fixed $\eta$.

Fix such a set of cardinality $C$ and a possible $d$. Write $v=q'-r$, and restrict $|v|$ to a dyadic band $[T,2T)$. For every prime power dividing $d$ exactly, $r-s$ must equal one of the $C$ selected shifts modulo that power. There are therefore at most $C^{\omega(d)}$ acceptable residue classes modulo $d$ by CRT. The function $f(r)=r-\sqrt{r(r+v)}$ satisfies
$$
 |f'(r)|=\frac{v^2}{2\sqrt{r(r+v)}\,(2r+v+2\sqrt{r(r+v)})}
          \asymp T^2/R^2 .
$$
Its range on the allowed interval has length $O(T^2/R)$ and $r-s=\lceil f(r)\rceil$. For one class modulo $d$, at most $O(1+T^2/(Rd))$ possible integers in that class meet the rounded range. Each has $O(R^2/T^2)$ preimages because $f$ is monotone with the displayed derivative. This bounds the number of possible $r$ by $O(R/d+R^2/T^2)$ for each class.

Only $O(T/d)$ nonzero $v$ in this band are multiples of $d$, since $d\le |v|<2T$. Multiplying the two counts and dividing by $R^2$ gives $O((T/R+d/T)/d^2)=O(1/d^2)$ per acceptable class. Restoring the prime-pair normalization and summing the $O(\log R)$ bands costs $X^{o(1)}$. For fixed $C$, $C^{\omega(d)}=d^{o(1)}$, so the sum over $d\ge X^u$ costs at most $X^{-u+o(1)}$. Apply the lower bound on $d$ above, then let $\eta$ decrease to zero. This proves the claimed probability bound, including input-dependent selection of a coprime subcollection.

## Finite checks

The check exhausted 685 unit product fibers across several small squarefree moduli. For primes in $[128,256)$, 138 ordered distinct pairs met a test threshold, compared with 176 first-crossing witness incidences and an elementary integer-count upper bound of 2,540. On a fixed input, all 2,880 possible randomized tables at primes $5,7,11,13$ gave probability $137/1440$ for exactly two hits and $101/960$ for agreement product at least 35. For the square remainders, it also checked 3,414 gcd identities and the finite derivative-based bound on 2,816 rounded residue classes.

For each prime size $b$, a separate seeded check sampled 256 distinct balanced inputs, used all primes in $[(2b)^2/2,(2b)^2]$, and graded agreement product at least $N^{1/4}$ with the smaller prime. The cubic table uses $a_\ell(u)=u^3$. The randomized table uses independent multipliers in $a_\ell(u)=c_\ell u^2$. The constructor sees only the public product and the table randomness.

| Bits per prime | Band primes | Cubic table hits | Randomized table hits |
| ---: | ---: | ---: | ---: |
| 16 | 75 | 27/256 | 13/256 |
| 24 | 152 | 0/256 | 1/256 |
| 32 | 255 | 2/256 | 0/256 |
| 40 | 382 | 2/256 | 1/256 |

These counts grade raw agreement. The agreeing prime coordinates were read using the known synthetic factor. [Executable checks](check_local_order_scores.py)
