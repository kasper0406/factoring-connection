# A rectangle evaluator for the modular ellipse

The [47-bit cutoff pair](class_selection.md#next-blocker-and-acceptance-test) has signed residues 0 at 581 and 4 at 582. A public projection through the discriminant reproduces both. Grouping integer multiples of a candidate before testing its square saves work specifically in parity. An additional local test eliminates candidates on whole norm levels. On 582 the new evaluator used 18,042 candidate square tests for positive parity and 19,089 for the full signed residue. The reference used 7,971,513. The new method still scans at least $B$ residues on a full disk; the improvement in this instance does not give a better factoring exponent.

## Move the sum from rows to congruence classes

Use the repository's public principal form and eligible branch. Write its positive-quadrant equation as $x^2+D y^2=L j$. For $4\mid\Delta$ take $(L,D)=(N,\Delta/4)$ and $x=w/2$. For odd $\Delta$ take $(L,D)=(4N,\Delta)$ and $x=w$, where $w=|2x_{\rm original}+(\Delta\bmod 2)y_{\rm original}|$. The equality enforces the original parity. In both cases $\gcd(L,D)=1$. In the decoder range $0\le B\le\lfloor\sqrt{\Delta/3}\rfloor$ the punctured disk has no axis points. Each positive pair contributes four signed vectors. Consequently positive counts modulo 4 give the complete signed residue modulo 16, and positive parity is a separate one-bit output. If counting modulo simultaneous negation, divide the signed multiplicity by two.

For a rectangle of norm levels and positive coordinates, first apply the public projection
$$j\equiv L^{-1}x^2\pmod D.$$
The input supplies neither a factor of $N$ nor a mixed square root. Factor only the public $D$. Choose the smallest product $m$ of full prime-power factors of $D$ at least as large as the width of the level interval. The decoder cap guarantees $D$ is large enough. For $x=r\bmod m$, there is at most one admissible level. Scan residues $r$ once, retaining only those whose implied level lies in the interval. A rejected residue discards every row in that class for the entire interval. For the retained levels, scan each remaining prime-power modulus once and combine compatible roots by CRT **modulo $D$**. These operations never split the quadratic congruence modulo $N$.

Each surviving level $t$ and residue $r\bmod D$ leaves a short progression $z=r+kD$ within $1\le z\le\lfloor\sqrt{Lt-D}\rfloor$. The remaining square test is on $v=(Lt-z^2)/D$. For a plain query one would simply test these values and clip the positive row to the rectangle. The following grouping further reduces the number of tests.

## Cancel multiples before the square test

Let $E$ be the product of primes whose squares divide $D$. For any candidate $(x,j)$ take the largest $h$ coprime to $E$ with $h\mid x$ and $h^2\mid j$. Then $x=h z$, $j=h^2t$, and no prime outside $E$ divides $z$ with its square dividing $t$. For an actual solution $h^2\mid D y^2$. Each prime of $h$ occurs at most once in $D$, so $h\mid y$. This proves that every solution has a unique base solution $(z,y_0,t)$, and that every multiplier $h$ coprime to $E$ restores a valid solution $(hz,hy_0,h^2t)$. The extraction is defined using $x,j$ alone. Square factors of $D$ are the reason for restricting the multipliers.

For $j_0\le j\le B$, write $h_-=1+\lfloor\sqrt{(j_0-1)/t}\rfloor$, $h_+=\lfloor\sqrt{B/t}\rfloor$. The base contributes with coefficient
$$M=\sum_{e\mid E}\mu(e)\left(\left\lfloor\frac{h_+}{e}\right\rfloor-
                         \left\lfloor\frac{h_--1}{e}\right\rfloor\right).$$
On an unmasked disk this depends only on $t$. A zero coefficient in the output ring rejects all its progressions before the square test or the CRT lift. A base can come from $t<j_0$, so the ray version begins the projection at 1 even for a thin requested interval.

For a row window $a\le y\le b$, intersect the multiplier interval with $a^2\le h^2 v\le b^2$ using integer square roots, then apply the same inclusion-exclusion. This is valid even before knowing whether $v$ is a square. For a surviving base test $v$ once. A nonsquare contributes zero for every positive multiplier. A square contributes $M$ and gives the exact sharp boundaries with no enumeration of its multiples.

Two further filters apply to the base. First, $z$ must be coprime to the product of primes $p\nmid E$ with $p^2\mid t$. Second, suppose an odd prime $p\mid t$, $p\nmid LD$, has $(-D/p)=-1$. Reduction of $z^2+D y_0^2=Lt$ modulo $p$ forces $p$ to divide both coordinates. This forces $p^2\mid t$ and contradicts the defining base condition. Thus a single such prime rejects the whole level. Both tests use factorizations of public candidate levels. The implementation caches them and records that work.

## What the full query costs

The direct scan of the first modulus costs $m$ residue evaluations. Preparing the other local tables costs $\sum q_i$ residue evaluations for the remaining prime powers. Target lookups, their lists, inclusion-exclusion, factoring $D$ and the encountered levels, and CRT combinations are extra. For a surviving level $t$ let $R_t$ be the number of compatible complete residues modulo $D$. Then the number of actual coordinate positions visited is bounded by
$$\sum_t R_t\left(1+\frac{\sqrt{Lt}}{D}\right).$$
More precisely each residue $r$ contributes the number of its positive representatives up to $\lfloor\sqrt{Lt-D}\rfloor$. The code charges those positions even if the primitive check or the row clip later removes them.

For an explicit upper bound up to bit-length factors, put $s=\omega(D)$, $J=2^{\omega(E)}$, $V$ equal to the number of distinct levels at the seed, $K=\sum R_t$, and $H$ equal to the number of visited positions. In addition to factoring $D$ and the separately counted public base levels, one query costs
$$\widetilde O\left(2^s+m+\sum q_i+(s+J)V+(1+s)K+(1+J)H\right).$$
Here $JH$ allows for row clipping, and the exact square-root cutoffs account for the ellipse boundary. Only the modulus decomposition and base-level arithmetic are reused across queries in the decoder; the reported decoder totals charge rebuilding the residue tables. All timings include construction. The raw counters include polynomial evaluations during public factorizations separately from the final integer square tests.

There are two visible limitations. For a full disk $m\ge B$ by construction, and for prime $D$ and $B>1$ the implementation scans $m=D$. For every level $B/4<t\le B$ the only multiplier on a full disk is 1. Parity cannot cancel any primitive candidate in that band by scaling. Both limitations concern this algorithm.

## Controls and measurements

| $B$ on the 47-bit input | Positive disk | Signed residue mod 16 | Plain projection tests | Exact ray tests, no inert filter | Parity tests, no inert filter | Parity tests, filtered | Mod-4 tests, filtered |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 65 | 8 | 0 | 3,922 | 1,569 | 1,157 | 456 | 586 |
| 145 | 12 | 0 | 11,464 | 6,856 | 6,314 | 1,417 | 1,417 |
| 257 | 16 | 0 | 26,064 | 16,692 | 15,513 | 3,163 | 3,619 |
| 401 | 20 | 0 | 48,374 | 30,559 | 26,688 | 6,989 | 7,680 |
| 581 | 24 | 0 | 80,675 | 52,152 | 46,409 | 16,471 | 17,518 |
| 582 | 25 | 4 | 82,246 | 53,723 | 47,980 | 18,042 | 19,089 |

At 582 the mod-4 run scanned 2,809 seed residues plus 46 other local residues, factored 91 base levels, tested 93 Legendre symbols and visited 30,457 coordinate positions. The reference's 7,971,513 square tests and the new 19,089 are counts of the final integer-square checks in two different algorithms. The per-stage counts above and [raw results](rectangle_rays_results.json) allow setup costs to be compared independently. The modulus shape matters. At $N=19{,}498{,}647{,}323$, $B=71$ the parity run visited 63,667 setup residues for 463 final square tests.

The independent direct area enumeration on the 194 repository branches and ten additional ties agreed on all 22,788 exact or modular rectangle comparisons. These include narrow row and level ranges, empty cutoffs, both discriminant parities, and 49 branches with a square factor in $D$. All 204 recovered a verified prime divisor, including the collision pair and 12 tied first minima. Fifteen N-only scaling cases also agreed with reduced-form enumeration at their first mixed cutoff. On eight further fixed inputs from the earlier held-out cohort, 58 parity and mod-4 comparisons agreed with an independent count from the principal and squared prime classes. That last reference uses fixture factors; neither the rectangle evaluator nor its decoder takes them.

**Factoring from the rectangle interface.** The known principal contribution is $\lfloor\sqrt B\rfloor$ positive points up to the cap. Use the [existing doubling-and-bisection routine](../ellipse_row_batches/row_batches.py) on the remaining count modulo 4. Before twice the first mixed norm there are at most three mixed pairs, so positivity is detectable there, including a tie. At the first norm $a$, bisect in the positive $y$ coordinate using the rectangle $[a,a]\times[1,y_0]$. Subtract the known principal row when $a$ is square. The remaining monotone count is 0, 1 or 2. For the first qualifying $y$, verify $w^2=4Na-\Delta y^2$ and take $\gcd((cy-\ell w)/2,N)$. On $140{,}000{,}556{,}000{,}551$ this input-only run located $a=582$ and returned $14{,}000{,}029$ in 42 queries and 251,868 candidate square tests. This uses row-clipped mod-4 queries in addition to the radial signed count.

## Next test

At 582, **16,546 of the 18,042 filtered parity tests** were at surviving nonsquare levels $B/4<t\le B$. The very same upper-band count was 14,975 tests at 581. Only one genuine mixed point enters at 582. Local compatibility is easy to exhibit even below the first hit: $39389^2\equiv 7N\pmod{1455062}$, but the positive disk has no mixed point at level 7.

The next target is the parity of this primitive nonsquare band as one aggregate. The residual quadratic in a terminal progression is
$$v(k)=(Lt-r^2)/D-2rk-Dk^2.$$
Test a rule that combines these progressions across different $(t,r)$ before generating all of them. Check the zero and one outcomes at the adjacent cutoffs, then the other inputs. To improve the exponent it must also avoid the $m\ge B$ seed scan. Further optimizations limited to multipliers of one primitive direction leave the observed upper band untouched.

Run from the repository root with Python, python-flint, SymPy and NumPy. The benchmark also imports the existing disk tests and their dd dependency.

    python experiments/ellipse_disk_mod8/rectangle_rays.py 140000556000551 --cutoff 582
    python experiments/ellipse_disk_mod8/rectangle_rays.py 140000556000551
    python experiments/ellipse_disk_mod8/rectangle_rays_check.py --output rectangle-rerun.json
