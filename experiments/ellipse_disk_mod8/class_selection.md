# Can the public middle coefficient select the prime class?

On the near-square $\ell=1$ branch, short rational approximations to the ratio of the factors give a public shortlist before any classes are enumerated. The rule tested here recovered 14 of 96 fresh balanced semiprimes. Spending the same number of class-minimum probes on the smallest possible indices recovered 13. The full list recovered all 96. A counting argument explains the loss of coverage and applies to other rules confined to small prime-class minima. We would stop tuning this shortlist and put the next experiment on bulk evaluation of a single nonsquare cutoff.

## Selection rule

Let $N=pq$, $p<q\le2p$, $c=\lfloor2\sqrt N\rfloor$ with the [usual downward adjustment](../ellipse_aggregate/factor.py) if $4N-c^2\le4$. Write $\delta=2\sqrt N-c\in(0,2)$ and $\Delta=4N-c^2$. The prime class up to inversion has representative $[p,-c,q]$. Write $m$ for its reduced minimum. This is the minimum of $C$ itself. The index $s$ in the [whole-prefix decision](class_prefix.md) is the reduced minimum of $C^2$.

For positive coprime weights $u,v$, put $H=pv^2+qu^2$. Then

$$H-2uv\sqrt N=p(v-u\sqrt{q/p})^2,\qquad
A=H-cuv=[p,-c,q](v,u).$$

A good approximation $v/u$ puts the integer $H$ just above $2uv\sqrt N$. The search uses $K=\lfloor N^{1/16}\rfloor$, $1\le u\le K$, $u\le v\le\lceil\sqrt{2}u\rceil$, and the first four integers $H>2uv\sqrt N$. It keeps each distinct $A=H-cuv$ in $[1,\lfloor\sqrt{\Delta/3}\rfloor]$. All endpoints use integer square roots.

1. Each proposed $H$ also gives a direct test. If $H^2-4Nu^2v^2$ is a square, divide $H\pm\sqrt{H^2-4Nu^2v^2}$ by $2v^2$, requiring exact division, a proper divisor $p$ of $N$, and $p\le N/p\le2p$. This certifies the weight represents $A$ in the true prime class.
2. For each proposed $A$, enumerate roots $b^2=-\Delta\pmod{4A}$ and construct the primitive reduced $[A,|b|,d]$. Run the existing public [batched verifier](class_prefix.py) on the resulting forms. This can find a factor even when a proposal hit the true minimum through an unrelated weight. It can miss when the true minimum is outside the shortlist. The command returns an empty answer in that case.
3. The controls use the same root enumeration and verifier on (i) the equally many lowest indices and (ii) the full cap. Since this experiment searches $C$ directly, these controls can enumerate reduced forms at the given leading indices without lifting their squares.

The proposal costs $O(K^2 g)$ integer-arithmetic operations for gap length $g\ge1$, before roots and verification. A cheap proposal is useful only if it covers the true class on enough inputs. For any weight searched, $0<A<\delta uv+g\le 2K(\sqrt2K+1)+g$. Thus every successful shortlist lookup has $m=O(K^2+g)$. A direct certificate gives $m\le A$, so the same bound applies even when the certified vector is longer than a shortest vector.

## How often can the prime class have a small minimum?

Choose the two distinct primes in $[P,2P]$. For $1\le T\le\sqrt P$ the fraction of pairs whose near-square prime class has $m\le T$ is

$$\ll (\log P)^2\left(\frac{T}{\sqrt P}\right)^{2/3}.$$

Here is a counting proof. Split at a parameter $0<\eta\le1$. When $\delta<\eta$, the product $N$ lies between $c^2/4$ and $(c+\eta)^2/4$. The $O(P)$ possible integers $c$ each give $O(P\eta+1)$ integers $N$. A semiprime determines its unordered pair, so this accounts for $O(P^2\eta+P)$ pairs.

For $\delta\ge\eta$, take a vector attaining $m\le T$. For sufficiently large $P$, $T<P$, so neither coordinate can vanish. The minimum can be taken with $u,v>0$ since $c>0$. With $r=\sqrt{q/p}$,

$$p(v-ru)^2+\delta uv=m\le T.$$

It follows that $|v-ru|\le\sqrt{T/P}$, $v\asymp u$, and $u\ll\sqrt{T/\eta}$. For a fixed $p,u,v$, the admissible values of $r$ form an interval of length at most $2\sqrt{T/P}/u$. Since $q=pr^2$ and $1\le r\le\sqrt2$, this gives $O(1+\sqrt{PT}/u)$ integer possibilities for $q$. There are $O(u)$ choices of $v$ per $u$. Summing over $u$ and all $O(P)$ integers $p$ gives $O(PT/\eta + P^{3/2}T/\sqrt\eta)$ further pairs.

Divide by the $\asymp P^2/(\log P)^2$ distinct prime pairs in the interval and take $\eta=(T/\sqrt P)^{2/3}$. This proves the bound, including the exceptional inputs where $2\sqrt N$ is very close to an integer. The only prime distribution input is the [prime number theorem](https://mathworld.wolfram.com/PrimeNumberTheorem.html).

In particular a search confined to $m\le N^{1/4-\epsilon}$ covers a vanishing fraction for any fixed $0<\epsilon<1/4$. For this four-position rule $m=O(N^{1/8})$ on success, giving a fraction at most $N^{-1/12+o(1)}$ on these dyadic balanced inputs. This concerns the public $\ell=1$ family and small minima of $C$. Sparse proposals near $N^{1/4}$, aggregates over unlisted classes, and other public branches are open.

## Blind benchmark

The four scales use 24 distinct uniform unordered pairs of primes per scale, sampled without replacement within each tier by rejection from odd integers. The seed, parameter rule and gap are fixed in the [benchmark](class_selection_check.py). The selector and both control searches receive only $N$. The test uses the factors afterwards to label $m$ by Gauss reduction and check every output. The benchmark runs both paths even on a direct hit to measure class retention separately. Each listed search includes modular-root generation and the product verifier on misses.

| Size of $N$ | Direct certificate | Retained true class | Combined | Same number of lowest indices | All indices |
|---|---:|---:|---:|---:|---:|
| 21–22 bits | 7/24 | 8/24 | 8/24 | 8/24 | 24/24 |
| 29–30 bits | 2/24 | 4/24 | 4/24 | 4/24 | 24/24 |
| 37–38 bits | 0/24 | 2/24 | 2/24 | 1/24 | 24/24 |
| 45–46 bits | 0/24 | 0/24 | 0/24 | 0/24 | 24/24 |

In the largest tier the proposal visited 512 weight pairs and tried 1,200 integer positions in total. It generated 218 classes from 715 distinct-minimum probes summed across inputs, followed by 22,512 polynomial evaluations. Exhaustive search visited 46,805 minima, generated 10,243 classes and used 5,274,635 evaluations. The two totals pay for substantially different success rates. The [per-input results](class_selection_results.json) include setup, build-plus-search time and work counts for all three methods.

Reusing the 194 repository semiprimes with newly constructed $\ell=1$ branches gave 124/194 successes, the same total as the lowest-index control. The original historical branches use different public data. On the two collision integers the shortlist factors 11,009 and misses 28,829. Exhaustive search factors both. The class builder also agreed with independent signed-form enumeration on the full cap and a sparse subset for each of 248 discriminants through 500.

## Next blocker and acceptance test

Compute the signed punctured-disk residue $U_N(B)\bmod16$ for a growing nonsquare cutoff by combining contributions over an interval of norm levels. One concrete starting point is to evaluate the parity of the number of positive solutions $x^2+d y^2=Nj$ over a rectangle of $(j,y)$, with $j\le B$ and $d=\Delta/4$ for the even-discriminant fixture below. A general evaluator needs a rule for combining these blocks that keeps the congruence and the sharp boundary exact. Count the cost of constructing and combining every block. For this fixture the axes are empty, so four times the positive count gives the signed count; its full residue needs the positive count modulo 4. Parity already detects the change in this test.

A specific floor-sum starting point is $L(t)=\lfloor\sqrt{\max(t,0)}\rfloor$. With $Y=\lfloor\sqrt{(NB-1)/d}\rfloor$ the positive count is exactly
$$\sum_{j=1}^{B}\sum_{y=1}^{Y}\bigl(L(Nj-dy^2)-L(Nj-dy^2-1)\bigr).$$
The difference is one precisely at a positive square. The first open step is a reciprocity or subdivision rule for this *difference on a whole rectangle*, in which interior contributions cancel before listing the shells or tracing their boundaries one at a time. The displayed identity on its own has the original enumeration cost.

| Public $N$ | $c$ | $\Delta$ | $B$ | Positive disk | $U_N(B)\bmod16$ |
|---:|---:|---:|---:|---:|---:|
| 140,000,556,000,551 | 23,664,366 | 5,820,248 | 581 | 24 | 0 |
| 140,000,556,000,551 | 23,664,366 | 5,820,248 | 582 | 25 | 4 |

The reference uses only public inputs and enumerates the shells with a small square-residue sieve. It performs 7,971,513 integer square tests through 582. The new positive point has $(x,y)=(175811278,186427)$ and gives $\gcd(cy/2-x,N)=14000029$. These coordinates are a check on the reference. The proposed evaluator gets only the public branch and requested cutoff. Its first milestone is to reproduce both residues and the existing small exact cases while showing how its work grows across increasing nonsquare cutoffs. Matching these two numbers alone gives no complexity guarantee.

For reproduction with Python, SymPy, NumPy and python-flint installed, run

    python experiments/ellipse_disk_mod8/class_selection.py 11009
    python experiments/ellipse_disk_mod8/class_selection_check.py --output results-rerun.json
