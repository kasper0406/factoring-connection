# Fixed cap moment by discriminant progressions

The whole fixed-cap moment can be evaluated using public modular roots in the discriminant. Combining this method with the repository's rational-tangent shell enumeration gives a worst-case expected $N^{4/9+o(1)}$ bit bound for this balanced-semiprime evaluator. This improves the bound for the fixed scalar relative to the two evaluators described in the weighted disk README. The exponent still gives exponential time in the input bit length and is larger than the exponent of standard general factoring methods. [Weighted disk README](https://github.com/kasper0406/factoring-connection/blob/9514cbfa013cd9ae3f030044db1dc3fb6bb4a9ed/experiments/ellipse_aggregate/README.md)

## Progression evaluator

Let $A=\lfloor\sqrt{\Delta/3}\rfloor$ be the public cap. For even $\Delta$ put $x=w/2$, $D=\Delta/4$ and $H_j=Nj$. For odd $\Delta$ put $x=w$, $D=\Delta$ and $H_j=4Nj$. Each shell solves $x^2+D y^2=H_j$ in positive integers. In the even case the original equation forces even $w$. In the odd case the equation already enforces its required parity.

Factor $D$ by trial division once. At each level $j$ solve $x^2=H_j\pmod D$ and enumerate the resulting disjoint progressions up to $\lfloor\sqrt{H_j-D}\rfloor$. For each $x$ test whether $(H_j-x^2)/D$ is a square and add the corresponding original weight. This uses the targeted-root idea from the upper-band evaluator. Here it evaluates all levels of the exact absolute moment, and represents nonunit roots as progressions rather than listing every lift modulo $D$. [Upper-band evaluator](https://github.com/kasper0406/factoring-connection/blob/9514cbfa013cd9ae3f030044db1dc3fb6bb4a9ed/experiments/ellipse_disk_mod8/conductor_band.md)

For a prime power $\ell^e\parallel D$, let $v=\min(e,v_\ell(H_j))$. When $v=e$, the roots form one residue class modulo $\ell^{\lceil e/2\rceil}$. When $v<e$ is odd, there are no roots. Otherwise divide out $\ell^v$, find at most four unit roots modulo $\ell^{e-v}$, and scale by $\ell^{v/2}$. This gives at most four classes modulo $\ell^{e-v/2}$. Unit roots are found by modular powering or randomized Tonelli–Shanks and digit lifting, with the binary case lifted directly. CRT combines the classes with at most $4^{\omega(D)}=D^{o(1)}$ work per level up to polynomial bit factors.

The combined progression step at level $j$ is $D/s_j$, where $s_j^2\mid\gcd(D,H_j)$. Since $D$ is coprime to $N$, each possible $s$ divides $\prod_{\ell^e\parallel D}\ell^{\lfloor e/2\rfloor}$ and satisfies $s^2\mid4j$. Consequently
$$\sum_{j\le A}s_j\le\sum_{s^2\mid D} s\left\lfloor\frac{A\gcd(4,s^2)}{s^2}\right\rfloor\le 4A\sum_{s^2\mid D}\frac1s \le A D^{o(1)}.$$

Counting even empty progressions, the total cost is therefore
$$N^{o(1)}\left(\sqrt D+A+\frac{A\sqrt{NA}}D\right)=N^{o(1)}\left(\sqrt\Delta+\sqrt N\,\Delta^{-1/4}\right).$$
This includes discriminant factorization, all root queries, CRT, and candidate square tests. Randomness affects the time to find prime-field nonresidues. Each returned root is exact.

## Crossover with geometric shells

On the other branch, the repository's rational-tangent evaluator covers the whole shell in $O(1+(Nj/\sqrt\Delta)^{1/3})$ blocks. This follows from dividing each of its two graph lengths by its rational-tangent block length and charging the terminal blocks. A constant number of lines is tested per block. Deduplicate the intersection of the two graph halves and accumulate positive $w$ for the prescribed scalar. Summing over $j\le A$ gives $N^{o(1)}(\sqrt\Delta+N^{1/3}\sqrt\Delta)$. [Rational-tangent implementation](https://github.com/kasper0406/factoring-connection/blob/9514cbfa013cd9ae3f030044db1dc3fb6bb4a9ed/experiments/ellipse_row_batches/arc_batches.py)

Use geometric shells if $\Delta^9<N^2$ and progressions otherwise. At the crossover both main terms are $N^{4/9}$, and $\sqrt\Delta=O(N^{1/4})$ throughout the public near-square branch. Decode the resulting scalar by the original fixed-cap formula. Neither evaluator requires the factors of $N$. The crossover bound follows from the two per-method bounds. The standalone script exercises the progression method.

## Checks

The progression code matched brute-force congruences for 7,260 residue queries with moduli up to 120. The end-to-end checks covered 560 small balanced semiprimes, including any direct branch gcds, and matched independent integer shell sums on the eligible disks. The two historical auxiliary branches gave traces 288 and 630. At $N=38{,}942{,}974{,}393$ the standalone program used 68 trial divisor probes, 358 levels, 734 prime-power queries, 480 progressions, and 5,967 square tests. It computed the published residue 69,529,082 and returned 254,119.

Run `python progression_moment.py 38942974393` with SymPy available.
