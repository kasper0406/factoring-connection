# How much can a short list of primes tell us about a factor?

The October experiment already combined prime powers from $p-1$ and $p+1$, tried the possible residues, and used Coppersmith to recover $p$. It split all 128 planted mixed-sign cases. Its large divisor had to divide both $p^2-1$ and $q^2-1$, making it visible in $N^2-1$. Here we remove the condition on $q$ when a short list of candidate primes is available. We also bound how often a polynomial in $N$ can supply such primes on random balanced inputs. [October Theorem S, §§5 and 9](../../notes/specialfam_results.md) · [Original implementation](../../src/specialfam/smoothcopp.py)

## Recovery with one structured prime

Let $N=pq$ with primes $p<q\le2p$. For a set $S$ of primes write $[m]_S=\prod_{\ell\in S}\ell^{v_\ell(m)}$. Suppose $S$ is given, $s=|S|=O(\log\log N)$, and for either factor $r$,
$$
 [r^2-1]_S>N^{1/4+\epsilon}
$$
for fixed $\epsilon>0$. We can factor $N$ deterministically in polynomial bit time. The inputs are $N,S,\epsilon$. The relevant valuations and signs are unknown. Only bases at most $2N$ are needed. Reading a list with larger entries incurs its representation cost.

Check the bases for a common factor with $N$. Use the integer bit quantum $b=\lfloor\epsilon\lfloor\log_2N\rfloor/(2s)\rfloor$. For the true prime powers $\ell_i^{v_i}$ in $r^2-1$, take $w_i=\lfloor\log_2(\ell_i^{v_i})/b\rfloor$ when $b\ge1$. Their sum is at most $K=\lfloor\log_2(2N)/b\rfloor=O_\epsilon(s)$. Enumerate the nonnegative vectors with sum at most $K$. There are $\binom{K+s}{s}=\exp(O_\epsilon(s))$. Replace each $w_i$ by the smallest nonnegative exponent $e_i$ with $\ell_i^{e_i}\ge2^{bw_i}$. For the true vector $e_i\le v_i$, and the total loss in the logarithm of the product is less than $bs$ bits. Thus one enumerated $d=\prod_i\ell_i^{e_i}$ divides $r^2-1$ and exceeds $N^{1/4}$. When $b=0$, enumerate the exact exponents with sum at most $\lfloor\log_2(2N)\rfloor=O_\epsilon(s)$. All comparisons can use integers.

For each $d\ge N^{1/4}$, try the roots $a^2=1\bmod d$. There are at most $2^{s+1}$, including the four possible roots at a power of two. For each $a$, the possible $r=a+dy$ put $y$ in an interval of length $O(\sqrt N/d+1)=O(N^{1/4})$. A bounded number of certified windows for the linear polynomial $y+ad^{-1}\bmod N$ covers the interval. Verify every output by division. This step reuses the October root enumerator and the established affine lattice algorithm. Rounding the exponent vectors allows us to use them even when $N^2-1$ supplies no large common divisor. [Published lattice bounds](https://cr.yp.to/ntheory/smallheight-20080502.pdf) · [Reader](smooth_order_cover.py) · [End-to-end checks](check_smooth_order_cover.py)

For example, $N=761987=701\cdot1087$. For $S=\{3,5\}$, the supported parts of $701-1$ and $701+1$ are 25 and 27. Both lie below $N^{1/4}$, while their product 675 exceeds $N^{3/8}$. The full common divisor $\gcd(701^2-1,1087^2-1)$ is just 24. The reader recovered 701 from $N$ and the two bases at its default $\epsilon=1/8$. It also recovered factors of planted 47-bit and 135-bit products with the same bases and common divisor 24, using the affine lattice in both cases. In the exhaustive small check it recovered a factor in all 476 qualifying combinations of an input and a base list. Only the larger prime met the promise in 136 of these.

## Can we obtain the bases from the input?

Select distinct primes $p,q$ uniformly in $[R,2R)$, put $X=R^2$, and fix $0<t\le1/4$. Then $X\le N<4X$. These estimates count inputs satisfying the divisibility promise. A decoder can have other reasons for splitting an input.

**A preselected pool.** Allow a pool of at most $(\log R)^A$ arbitrary primes for fixed $A$, declared before drawing the inputs. Even selecting the best $s=O(\log\log R)$ of them after seeing both factors attains $[r^2-1]_S\ge N^t$ with probability at most $X^{-t+o(1)}$. Indeed there are only $\exp(O((\log\log R)^2))=X^{o(1)}$ choices of bases and exact exponents, since $r^2-1<4R^2$. Every resulting $d\ge X^t$ has at most $2^{s+1}$ roots of unity, and each class accounts for at most $1+R/d$ integers. Divide by the number of primes in the interval and allow a factor two for the choice of $r$. For a fixed finite $S$, Bugeaud, Evertse and Győry already give the order of magnitude of the integer count for values of $x^2-1$ with large $S$-part. [Published fixed-set density theorem, Theorem 2.3](https://doi.org/10.4064/aa170828-7-3)

**Primes dividing a polynomial in the public product.** Fix a primitive integer polynomial $F$ of positive degree $h$. Grant every prime dividing $F(N)$ for free. Let $D$ contain every power from $r^2-1$ at those primes, even if a prime divides $F(N)$ only once. For fixed $h$,
$$
 \Pr\{D\ge N^t\text{ for either }r\}\le X^{-t+o(1)}.
$$
This covers $F(N)=N^2-1$ and products of any fixed number of fixed-degree polynomial values. The same bound holds for degree polynomial in $\log R$ when $O(\log\log R)$ supporting primes can account for the threshold. Taking the product of a polynomial-size menu of preselected primitive polynomials gives the same statement for its combined pool. Coefficients may vary with the size band and must be preselected for that band. Primitivity ensures a nonzero reduction at each prime. The at most $h$ roots with $F(N)=0$ have probability $X^{-1+o(1)}$ and can be counted separately.

Set $T=X^t$ and consider $r=p$. If a supporting prime $\ell>T$ occurs, first fix $p$. There are $O(\log R)$ prime factors of $p^2-1$, and at each such $\ell$ at most $h$ residues of $q\bmod\ell$ satisfy $F(pq)=0\bmod\ell$, since $p$ is a unit. Counting all integers in those classes and dividing by the number of prime partners bounds the probability by $X^{-t+o(1)}$.

Otherwise multiply supported prime factors of $p^2-1$, allowing repetitions, up to the first product $d\ge T$. Then $T\le d\le T^2\le R$, $p^2=1\bmod d$, and $\operatorname{rad}(d)\mid F(pq)$, where the radical is the product of the distinct primes in $d$. The first congruence gives at most $2\cdot2^{\omega(d)}$ classes for $p\bmod d$. For each such $p$, the second gives at most $h^{\omega(d)}$ classes for $q\bmod\operatorname{rad}(d)$. Both moduli fit within the interval length $R$. The number of qualifying ordered integer pairs is therefore
$$
 O\!\left(R^2\,\frac{(2h)^{\omega(d)}}{d\,\operatorname{rad}(d)}\right).
$$
For fixed $h$, the numerator weight is $X^{o(1)}$ uniformly for $d\le R$. The same holds for $h=(\log R)^{O(1)}$ if $\omega(d)=O(\log\log R)$. For fixed $0<\sigma<1$,
$$
 \sum_{d\ge T}\frac{1}{d\,\operatorname{rad}(d)}
 \le T^{-\sigma}\prod_{\ell\ {\rm prime}}\left(1+
 \frac{\ell^{-2+\sigma}}{1-\ell^{-1+\sigma}}\right).
$$
The product converges. Divide the integer counts by the $\Theta(R^2/\log^2 R)$ ordered distinct prime pairs, let $\sigma$ approach one, and apply symmetry. At the quarter threshold even the entire prime support of $N^2-1$ qualifies on at most $X^{-1/4+o(1)}$ of these inputs. October's density discussion explicitly labeled the estimate for the narrower common-divisor family heuristic. [October density discussion, §6](../../notes/specialfam_results.md) · [Finite witness and residue-class checks](check_sparse_order_coverage.py)

**The exponent at the quarter threshold is attained.** Already $F(N)=N-1$ gives $X^{-1/4+o(1)}$ qualifying inputs, with an existential choice of the supporting primes. Take $K=\sqrt R/(\log R)^B$ for a sufficiently large fixed $B$. Bombieri–Vinogradov in the interval $[R,2R)$ gives $\gg R/\log R$ incidences $(p,k)$ with $K\le k\le2K$ and $p=1\bmod k$, since $\sum_{K\le k\le2K}1/\varphi(k)\gg1$. Each $p$ has at most $\tau(p-1)=R^{o(1)}$ such $k$. If sparse support is required, discard $p$ with $\omega(p-1)>C\log\log R$. Their incidences total at most $(\log R)^{-C\log2}\sum_{m\le2R}\tau(m)2^{\omega(m)}\ll R(\log R)^{3-C\log2}$, using $\tau(m)2^{\omega(m)}\le\tau_4(m)$. Choose $C$ sufficiently large. [Bombieri's large sieve theorem](https://doi.org/10.1112/S0025579300005313)

To each of the remaining $R^{1-o(1)}$ primes assign one $d=(p-1)/k$. These $d$ have order $\sqrt R(\log R)^B$. There are $R^{1/2+o(1)}$ possible values. Cauchy–Schwarz therefore gives $R^{3/2-o(1)}$ ordered pairs assigned the same $d$, even after deleting equal pairs. Each has $d\mid p-1,q-1,N-1$ and $d>N^{1/4}$. This counts the availability of a support. The lower bound at this endpoint gives a logarithmic margin and does not supply the fixed exponent margin of the recovery theorem.

## What happens on unplanted inputs?

Each row uses 512 seeded pairs. The threshold is $N^{5/16}$. For the first four columns $S=\{2,3,5,7\}$. “Shared” uses its part of $\gcd(p^2-1,q^2-1)$. “One sign” takes the best result among $p-1,p+1,q-1,q+1$. “Both signs” takes the better of $p^2-1,q^2-1$. “Additional” counts inputs qualifying only for the third test among those three. The last two columns grant the best four of the first 64 primes, or every supported power from $F(N)=(N^2-1)(N^2+1)(N^2+N+1)$, respectively. All columns are promise counts scored with known synthetic factors.

| Bits per prime | Shared | One sign | Both signs | Additional | Best four of 64 | Polynomial support |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 16 | 20 | 109 | 303 | 189 | 500 | 339 |
| 32 | 0 | 1 | 3 | 2 | 252 | 16 |
| 48 | 0 | 0 | 0 | 0 | 5 | 0 |
| 64 | 0 | 0 | 0 | 0 | 0 | 0 |

The polynomial-support check exhausts ordered distinct prime pairs for $R=16,64,256$, identifies a large prime or a first-crossing divisor for every qualifying pair, and verifies the finite bound for each divisor. [Coverage and witness checks](check_sparse_order_coverage.py)

**Research decision.** The mixed-sign decoder is available as a finishing step. A general factoring proposal along this route needs a way to choose bases from $N$ that reaches the threshold frequently. Replacing one preselected polynomial by another, granting the prime support of its value, or adjusting the signs falls within the bounds above.
