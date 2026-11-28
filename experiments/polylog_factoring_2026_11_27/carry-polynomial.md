# Factoring from few changes in binary carries

An $m$-bit composite $N$ can be split in deterministic polynomial time for fixed $h$ if some factorization $N=ab$, $a,b>1$, has at most $h$ changes in its sequence of outgoing binary carries. The same conclusion holds if the *slopes* of the carries change at most $h$ times. The algorithm takes only $N,h$ and the choice of difference order. It can return a proper factor on other inputs as well.

Write $L_a(X)$ for the polynomial whose coefficients are the binary digits of $a$, least significant first. Let $u_i=[X^i]L_aL_b$, $n_i=[X^i]L_N$, and $c_i$ be the outgoing carry of column $i$. Extend $c_i$ by zero on both ends. With $C=\sum_i c_iX^i$, the column equation gives

$$
u_i+c_{i-1}=n_i+2c_i,\qquad L_aL_b=L_N+(2-X)C .
$$

For order $k=1$ or $2$, enumerate the positions and nonzero coefficients of $J=(1-X)^kC$. At most $h$ positions occur in a range of length $m+1$. Every $u_i\le m$, so induction gives $0\le c_i\le m$. The jump magnitudes are therefore at most $m$ for first differences and $2m$ for second differences. Enumerating positions and heights gives $m^{O(h)}$ candidates. Enforce $J(1)=0$, and also $J'(1)=0$ for the second difference. These identities determine the final one or two jump heights after their positions are chosen. Integration recovers $C$.

Factor each admissible $F=L_N+(2-X)C$ over $\mathbb Z$. For every irreducible $G\mid F$, test $\gcd(N,G(2))$. At the true profile, the monic digit polynomial $L_a$ has some irreducible factor whose value at 2 has absolute value greater than one and divides $a$. Its gcd is therefore a proper divisor of $N$, including when $a,b$ are composite. Dense integer polynomial factorization has a deterministic polynomial-time algorithm. The executable uses FLINT when available and SymPy otherwise. [Lenstra, Lenstra and Lovász](https://doi.org/10.1007/BF01457454) · [Implementation](carry_polynomial.py)

### What a short description of the carries permits

**One plateau, unbounded total carry.** Let $L\ge3$, choose $H>L+1$, and take $0\in I\subset[0,H-1)$ so successive elements of $I\cup\{H-1\}$ differ by more than $L$. Form
$$
a=2^H+2^{H-1}+\sum_{i\in I}2^i,\qquad b=2^H+2^L-1.
$$
The low $L$ bits of $b$ translate each set bit of $a$ into a length-$L$ interval. The intervals starting below $H-1$ are disjoint. The two high intervals overlap from $H$ through $H+L-2$. The leading bit of $b$ adds one more term at $H$, and its other translates miss all intervals. Thus $c_i=1$ exactly for $H\le i<H+L$. There are two changes, while both the propagated sum $\sum c_i$ and the raw overlap $\sum_i\max(u_i-1,0)$ equal $L$. Both operands have $H+1$ bits, and their Hamming weights can both grow by taking $|I|$ and $L$ to infinity. This constructs integer products for every such choice. Infinitude of prime pairs in the family is open here.

**Growing carry height with four slope changes.** For $a=2^s-1,\ b=2^t-1,\ 2\le s\le t$, the first differences of the carry are $+1$ from 1 through $s-1$, zero through $t$, and $-1$ from $t+1$ through $s+t-1$, with zeros outside. Hence $J=(1-X)^2C$ has four coefficients $+1,-1,-1,+1$ at $1,s,t+1,s+t$. The carry reaches height $s-1$ and has total mass $t(s-1)$. Some such inputs already split through the unmodified digit polynomial.

**Frequency of few first differences.** At a downward step $c_{i-1}\le 2c_i+1$. With $h\ge2$ changes including the first rise, a carry must descend to zero in at most $h-1$ drops, so its peak is at most $2^{h-1}-1$. Evaluating the column identity at $X=1$ gives
$$
\operatorname{wt}(a)\operatorname{wt}(b)=\operatorname{wt}(N)+\sum_i c_i\le m\,2^{h-1}.
$$
For primes $p,q\in[R,2R]$, at least one thus has $o(m)$ set bits whenever $2^h=o(m)$. The number of such binary strings is $\exp(o(m))=R^{o(1)}$ by the binomial bound. The prime number theorem gives $\asymp R/\log R$ primes in the interval, so at most $R^{-1+o(1)}$ of the pairs satisfy this first-difference promise. The zero-carry case obeys the same bound directly. A separate low-column analysis covers sparse changes of either order. [Carry prediction bounds](carry-prediction-bounds.md)

### Precedents and checks

Landau and Immerman illustrate undoing carries on 1729 in §5 and explain why unrestricted coefficient splitting produces exponentially many candidates (pp. 14–15 of the July 1993 manuscript). Mackall's Algorithm 1 factors the unmodified digit polynomial at successive bases. A 2024 patent also treats carries using polynomial correction. The 2026 binary convolution height paper defines raw overlap on the unpropagated convolution. The two-change family above separates that statistic and the total outgoing carry from the number of changes. These comparisons establish a difference in promises without establishing priority for the sparse-difference enumeration. [Landau and Immerman](https://web.cs.umass.edu/publication/docs/1992/UM-CS-1992-066.pdf#page=14) · [Mackall](https://eoinmackall.github.io/files/papers/factorization.pdf) · [US 12,034,849](https://patents.google.com/patent/US12034849B2/en) · [Binary convolution height](https://doi.org/10.1080/27684830.2026.2727395)

| Input | Successful reconstruction | Result |
|---|---|---|
| $29{,}014{,}226{,}011{,}817{,}451{,}384{,}471{,}551$, 85 bits | Two changes at 42 and 60, total carry 18 and raw overlap 18. Tried 2,602 profiles and factored 1,454 polynomials. | $4{,}398{,}046{,}773{,}247 \times 6{,}597{,}070{,}815{,}233$. Fermat reaches this pair after 111,073,297,014 increments. |
| $1{,}099{,}503{,}239{,}183$, 40 bits | Four unit slope changes at 3, 20, 23 and 40, peak 17, total carry 340 and raw overlap 322. Tried 1,027 profiles and factored 474 polynomials. | $1{,}048{,}571\times1{,}048{,}573$. Fermat also finds this twin pair on its first test. |

The unmodified digit polynomial yielded no factor on either fixture. The 85-bit fixture's smaller prime also lies $2^{18}-1$ above $2^{42}$, within the usual known-high-bits range for balanced factors. The examples separate carry parameters and show no advantage over all applicable factoring methods. The check also verifies 182 schoolbook reconstructions and 155 promised input-only splits on small prime products. [Runnable checks](check_carry_polynomial.py) · [Carry recurrences](carry-recurrence.md) · [Carry parity](carry-parity.md)
