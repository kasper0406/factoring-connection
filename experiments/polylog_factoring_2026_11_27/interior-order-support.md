# Prime factors of all binary substrings

Choose distinct primes $p,q$ uniformly from $[R,2R)$, with $R$ an integer tending to infinity. Write $X=R^2$ and $N=pq$. Let $S(N)$ contain every prime dividing any positive proper contiguous substring of the binary expansion of $N$, read as an integer. For fixed $0<t\le1/4$, the chance that these primes account for at least $N^t$ of either $r^2-1$, $r\in\{p,q\}$, is at most
$$
 \exp\!\big[-(\sqrt{t/6}-o(1))\sqrt{\log X\,\log\log X}\big].
$$
This is smaller than every inverse polynomial in the input bit length. All powers of a supported prime in $r^2-1$ count, even if its first power appears only once in one word. A separate elementary count shows that a direct gcd of any proper word with $N$ splits at most $X^{-1/2+o(1)}$ of the inputs. Both estimates allow inspecting every word before choosing any of them. [Endpoint words and their sharper bound](digit-order-support.md) · [Sparse-order decoder](sparse-order-support.md) · [Checks](check_interior_order_support.py)

## Direct gcd with any proper word

Fix the position and width of a word. With powers of two $B,C$ and $H=BC$, write $z=\lfloor N/B\rfloor\bmod C$ and $J=N-Bz$. For a positive proper word, $J>0$. A discarded low bit includes the final one of odd $N$, and a discarded top bit includes its leading one. In particular, $J$ has the form $Hk+v$ with $0\le v<B$ and $0\le k<4R^2/H$. There are $O(R^2/C+B)=O(R^2/C)$ possible $J$ since $H\le8R^2$ for a window in the binary input.

If $p\mid z$, then also $p\mid J$. Each $J<4R^2$ has $O(1)$ prime divisors in the factor interval. For each such $p$, the prime partner has the form $q=J/p+Bj$ with $1\le j<C/p$. There are at most $C/R$ choices. Multiplication gives $O(R)$ possible ordered pairs per word position. Summing the $O(\log^2 R)$ positions and dividing by the $\Theta(R^2/\log^2 R)$ distinct prime pairs proves the bound.

## A large prime in an order and in a word

For $2\le Y\le2R+1$, the probability that *some* $\ell>Y$ dividing $p^2-1$ also divides a positive word at one fixed position is
$$
 \ll (\log(2R))^{O(1)}Y^{-1/3}.
$$
There are $O(\log R)$ candidate primes for a given $p$. For each word value range $0<z<C$, the set $E_p$ of values divisible by at least one of those primes therefore has size $O(C\log R/Y)$.

First take $C\le RY^{1/3}$. For fixed odd $p\in[R,2R)$, let $n_p(z)$ count all integer partners in that interval which give $z$. We use the collision bound
$$
 \sum_{p\ {\rm odd}}\sum_z n_p(z)^2
 \ll (R^2+R^3/C)(\log(2R))^{O(1)} ,
 \qquad R\le p<2R .
$$
Equal outputs with $q_1\ne q_2$ imply $pv$ is within distance less than $B$ of a multiple of $H$, where $v=|q_1-q_2|<R$. Each $(p,v)$ accounts for at most $2R$ ordered partners. We bound such $(p,v)$ by counting divisors of integers $n\le2R^2$ within these strips. The diagonal accounts for $O(R^2)$.

When $B\ge R^{1/4}$, Shiu's mean-value bound for the divisor function on intervals of length comparable to $B$ gives $O(B(\log R)^{O(1)})$ per strip. The interval near zero also satisfies this by the ordinary summatory bound. There are $O(1+R^2/H)$ strips, giving $O((B+R^2/C)(\log R)^{O(1)})$. The hypotheses use fixed exponents since the strip endpoints are $O(R^2)$ and $B\ge R^{1/4}$.

When $B<R^{1/4}$, we have $C\le RY^{1/3}\le 2R^{4/3}$, so $H\le 2R^{19/12}$. For each nonzero offset $b$ with $|b|<B$, write $g=\gcd(b,H)$. Divide the congruence $n=b\bmod H$ by $g$ and apply Shiu to $\tau(n/g)$ in the resulting coprime class modulo $H/g$ up to $2R^2/g$. The length divided by the modulus is $2R^2/H\gg R^{5/12}$. Since $g$ is a power of two, $\tau(n)\le\tau(g)\tau(n/g)$ costs $O(\log R)$, and $\varphi(H/g)\ge(H/g)/2$. The resulting bound is $O((R^2/H)(\log R)^{O(1)})$ for each $b$. At $b=0$, odd $p$ implies $H\mid v$, so there are $O(R^2/H)$ pairs directly. Summing the offsets proves the collision bound. [Shiu's theorem](https://doi.org/10.1515/crll.1980.313.161) · [Koukoulopoulos, *The Distribution of Prime Numbers*, §14](https://dms.umontreal.ca/~koukoulo/documents/publications/primes.pdf)

By Cauchy–Schwarz, the number of pairs with $z\in E_p$ is at most $(\sum_p |E_p|)^{1/2}(\sum_{p,z}n_p(z)^2)^{1/2}$. Dividing by the number of prime pairs bounds their probability by a polylogarithmic factor times $\sqrt{(1+C/R)/Y}$. This gives the claimed $Y^{-1/3}$ bound in this range.

For $C>RY^{1/3}$, fix $p$ and one candidate $\ell>Y$. Group the $q$ by $k=\lfloor pq/H\rfloor$. Each group is an interval. In it, $\ell\mid z$ requires $pq\bmod B\ell$ to lie among $B$ residues, since $\lfloor pq/B\rfloor\equiv Ck\bmod\ell$. Multiplication by $p$ permutes residues modulo $B\ell$. Counting over all groups gives $O(R/\ell+B(1+R^2/H))=O(R/\ell+R^2/C)$ partners. Summing over candidate primes and normalizing proves the same claim. These arguments include word positions whose exact availability depends on $N$ by bounding a larger set with positions fixed for the size band.

## Pooling all word positions

The polynomial number of positions in the bit length costs a factor $O(\log^2R)$ in the large-prime estimate. Put $L=\log X$, $W=\sqrt{L\log L}$ and $Y=\exp(aW)$ for fixed $a>0$. If there is no supported prime $>Y$ and the covered part reaches $X^t$, multiply its prime factors with repetitions until first reaching $d\ge X^t$. Then $d\le X^t Y<R$ and $p^2=1\bmod d$. The latter has at most $2\cdot2^{\omega(d)}$ roots. Counting all integers $p$ per root class and normalizing gives a bound of a polylogarithmic factor times
$$
 \sum_{\substack{d\ge X^t\\P^+(d)\le Y}}\frac{2^{\omega(d)}}d
 \le X^{-t\sigma}\prod_{\ell\le Y}\left(1+
              \frac{2\ell^{-1+\sigma}}{1-\ell^{-1+\sigma}}\right),
 \qquad 0<\sigma<1.
$$
Choose $\sigma=(\tfrac12\log L-2\log\log L)/\log Y$. The log of the Euler product is $O(Y^\sigma\sum_{\ell\le Y}1/\ell)=o(W)$, using the prime harmonic sum. The smooth part therefore costs $\exp[-(t/(2a)-o(1))W]$, while a supported large prime costs $\exp[-(a/3-o(1))W]$. Taking $a=\sqrt{3t/2}$ proves the assertion. Symmetry allows either factor.

If at most $s=O(\log\log X)$ supporting primes must suffice, there is also a stronger bound $(\log X)^{O(1)}X^{-t/(4s)}$. Set $U=X^{t/s}$. Some one supported prime power must reach $U$. A supported prime larger than $U^{3/4}$ costs at most a polylogarithmic factor times $U^{-1/4}$ by the large-prime estimate. For each smaller prime, choose its first power reaching $U$, which is at most $U^{7/4}<R$. The condition that this power divide $p^2-1$ gives at most four classes, each containing $O(R/U)$ integers. There are at most $U^{3/4}$ bases to try. This yields the other contribution of the same order.

### Shifting or multiplying the input before reading the bits

The bounds also hold after pooling words from every $aN+b$ with $1\le a\le(\log R)^K$ and $|b|\le(\log R)^K$, for any fixed $K$. Here $H\ll aR^2$ for large $R$. In the collision count the shift $b$ cancels. The product tested for proximity to a multiple of $H$ is $apv\le2aR^2$. Its number of representations for fixed $a$ is still at most its divisor count. The longer interval and larger number of strips cost at most a power of $\log R$ in both applications of Shiu. For the zero offset, $H\mid apv$ gives $H/\gcd(H,a)\mid v$. For large $C$, the map $q\mapsto apq+b\bmod B\ell$ has fibers of size at most $a$ since $p$ is a unit. The number of partners is therefore $O(aR/\ell + B(1+aR^2/H))$ in the interval-group count, allowing a further polylogarithmic factor in the endpoint term if needed. The union over coefficients and word positions is still polynomial in $\log R$.

For the direct gcd, $J=aN+b-Bz$ again has $O(aR^2/C)$ possible values. A factor $p$ dividing $z$ must divide $J-b$. If $J=b$, then $N\mid Bz$, hence $N\mid z$, which yields the whole modulus in the gcd. Otherwise each nonzero $J-b$ allows $O(1)$ choices of $p$ in the prime interval. Writing $z=ph$ gives $a q=Bh+(J-b)/p$ and at most $C/R$ positive choices of $h$. The number of pairs across all such transformations is $R(\log R)^{O(1)}$.

## Finite checks

Each row uses the same 128 unplanted prime pairs as the endpoint-word test. A covered part reaches the threshold when it is at least $N^{1/4}$ for either factor. The first three columns use the known test factors for grading. Only the last column is a direct factoring operation on $N$.

| Bits per prime | Best single proper word | Best two bases from all words | All supported bases from all words | Direct gcd |
| ---: | ---: | ---: | ---: | ---: |
| 16 | 127 | 121 | 127 | 0 |
| 32 | 60 | 38 | 113 | 0 |
| 48 | 4 | 2 | 80 | 0 |
| 64 | 0 | 0 | 61 | 0 |

At 64 bits per prime, pooling raised the graded count from 18 for end cuts to 61 for all contiguous words. A factorization of every word and all matching hidden valuations were available to the grader. The experiment also checked 238 size-and-position combinations at $R=16,32,64$, including the integer collision comparison and all direct-gcd witnesses. [Reproducible checks](check_interior_order_support.py)
