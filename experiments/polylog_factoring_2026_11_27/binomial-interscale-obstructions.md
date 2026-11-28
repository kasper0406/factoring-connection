# Can one dyadic binomial level be updated from the preceding levels?

Write $B_t=\binom{2t}{t}$. The fixed-multiplier decoder asks for the divisibility of $B_{\lfloor hN/D\rfloor}$, where $D$ varies over powers of two and $h=1,3$. Two obstructions apply to proposed updates using the numerical coefficients themselves. The first works modulo $N$ and allows the complete preceding residues as inputs. The second concerns exact rational duplication formulas with a fixed number of auxiliary quantities. Previous work bounded the spacing of neighboring indices and identified intervals on which individual shift ratios are units. [Two-stream decoder](two-fixed-and-streams.md) · [Neighboring-index calculation](binomial-transition-windows.md)

## The complete preceding residues can all vanish

Fix a power of two $D\ge16$ and any additional fixed modulus $L\ge1$. There are arbitrarily large balanced semiprimes $N=pq$, $p<q<2p$, of each of the following types. All of them satisfy $N\equiv1\bmod\operatorname{lcm}(D,L)$.

| Residues of both primes modulo $D$ | Every earlier eligible binomial modulo $N$ | Current $B_{\lfloor N/D\rfloor}$ and $B_{\lfloor3N/D\rfloor}$ |
|---|---|---|
| $-1$ | zero | both zero modulo $N$ |
| $D/2-1$ | zero | both units modulo $N$ |

Here an earlier eligible binomial may have **any** power-of-two denominator $d\mid D/2$ and any positive odd numerator $g$ with $4g\le d$. Even receiving all of these coefficients as residues, along with $N\bmod L$, leaves the next Boolean answers ambiguous on these two families.

To construct them, write $M=\operatorname{lcm}(D,L)$ and choose a unit $a\bmod M$ reducing to the indicated class modulo $D$. Such a lift exists by choosing a unit at every odd prime of $M$ and an odd lift to its 2-part. Both indicated classes square to one modulo $D$. Hence $a^{-1}\bmod M$ reduces to the same class as $a$ modulo $D$. The prime number theorem for fixed arithmetic progressions provides $p\equiv a\bmod M$ in $[X,6X/5]$ and $q\equiv a^{-1}\bmod M$ in $[7X/5,8X/5]$ for all sufficiently large $X$.

For $s\in\{p,q\}$ the established valuation criterion tests whether $g(N/s)\bmod d$ is in the upper half. The size hypotheses hold once $p>D/2$. Each lower denominator sees $-g\bmod d$ in both families. At $D$ the two possible residues are $D-h$ and $D/2-h$ for every fixed positive odd $h$ with $4h\le D$. One is in the upper half and the other is in the lower half. [Central-binomial valuation criterion](two-fixed-and-streams.md#exact-search-and-the-binomial-observations)

In particular, a fixed polynomial expression over the rationals in $N$, the preceding coefficients, and the index $x=(N-1)/D$ cannot give the next coefficient's zero-versus-unit behavior by reduction modulo $N$. On both families that expression reduces to its value at $N=0$, all previous coefficients zero, and $x=-1/D$. A fixed nonzero rational value is a unit for sufficiently large $p,q$. A zero rational value vanishes on both families. The argument also applies to rational expressions whose denominators are units on these inputs. It leaves open computations using the full integer $N$ in other ways, or measuring more than the preceding residues.

## No finite exact rational duplication state for the central binomial

For any distinct positive integers $m_1,\ldots,m_\ell$, the sequences $B_{m_1 n},\ldots,B_{m_\ell n}$ are algebraically independent over $\mathbb C(n)$ in the following precise sense. No nonzero polynomial in these sequences with rational-function coefficients in $n$ vanishes for every sufficiently large integer $n$. The result also holds after adjoining $4^n$ to the coefficient field and after replacing any $B_{mn}$ by $B_{mn+r}$ with fixed integer offset $r$.

For the proof, set $U_n=4^{-n}B_n$. Its ratio after a unit increment is $(n+1/2)/(n+1)$. Therefore

$$
 s_m(n)=\frac{U_{m(n+1)}}{U_{mn}}
       =\prod_{j=0}^{m-1}\frac{n+(2j+1)/(2m)}{n+(j+1)/m}.
$$

For a nonzero rational function of $n$, sum its signed zero and pole multiplicities in the orbit $-1/(2m)+\mathbb Z$. The result for $s_m$ is one. Every $s_{m'}$ with $m'<m$ has orbit sum zero because all its zero and pole denominators divide $2m'$ or $m'$. The orbit sum of $R(n+1)/R(n)$ is zero for every nonzero rational $R$, by translation. It follows that if any exponent $e_m$ is nonzero, then $4^a\prod_m s_m(n)^{e_m}$ fails the rational shift-ratio condition at the largest $m$ with $e_m\ne0$. If all $e_m=0$ and $a$ is a nonzero integer, the constant $4^a$ is also impossible because $R(n+1)/R(n)\to1$.

Consider a shortest rational-function linear relation among distinct monomials in $4^n$ and the $U_{mn}$, valid for all sufficiently large $n$. Divide by the coefficient of one monomial, shift by one, and subtract a multiple of the original relation to eliminate that monomial. Minimality forces every remaining coefficient to vanish. Thus for two distinct monomials, the ratio of their shift ratios must equal $R(n+1)/R(n)$ for a nonzero rational $R$. The preceding paragraph excludes that possibility. A one-term relation is excluded by nonvanishing. This proves the assertion about the $U_{mn}$ and $4^n$. Since $B_{mn}=(4^n)^mU_{mn}$, the assertion holds for $B_{mn}$ as well. A fixed index shift changes $B_{mn}$ by a nonzero rational function of $n$ by the ordinary binomial recurrence.

Two consequences specify the scope of the obstruction.

- For a **fixed finite** collection of dyadic levels, fix a residue class $N=Dn+e$ for their common finest denominator. Then each center $\lfloor hN/d\rfloor$ is $mn+r$ with $m=hD/d$, $r=\lfloor he/d\rfloor$. The distinct slopes for $h=1,3$ give algebraically independent exact coefficients as $n$ varies. Restricting a putative universal identity to this class suffices to refute it.
- Suppose an exact algorithm uses a fixed state of $\ell$ scalar quantities at $n$, reads $B_n$ from them by a rational function, and updates to the state at $2n$ by rational functions of $n$ and that state, with defined divisions. Composing its reads and updates makes $B_n,B_{2n},\ldots,B_{2^\ell n}$ rational functions of $\ell$ formal state variables over $\mathbb C(n)$. Any $\ell+1$ such functions are algebraically dependent. Their resulting identity on the actual states contradicts the theorem. This excludes fixed-width exact rational duplication even with arbitrarily chosen auxiliary scalar functions. It gives no complexity lower bound for growing states, modular computations, or a procedure that outputs only divisibility bits.

The companion check verifies the rational shift quotient, its distinguishing orbit at every slope tested, and both modular families against literal central binomials. [Finite checks](check_binomial_interscale_obstructions.py)
