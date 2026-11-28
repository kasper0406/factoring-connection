# How often can any subset of short fractions give a relation?

Combining fractions sometimes gives a relation that every singleton misses. At $N=247$ and $t=23$, three short fractions have outside-base squarefree parts $5,17,85$ and multiply to $-1/768$. The question for factoring is how often **any** such subset exists. Averaging over the input modulus lets us count all subset sizes at once. [Earlier example and polynomial-time evaluator](short-rational-relations.md#exact-construction)

## Experiment being bounded

Choose an unordered pair of distinct primes $p,q$ uniformly from $[X,2X]$ and set $N=pq$, $M=4X^2$, $L=\sqrt{\log M\log\log M}$. All logarithms are natural. For each $N$, fix a base $\mathcal B_N$ of at most $(\log M)^A$ primes, for a fixed $A$, before drawing targets. One row for a unit target $t$ is a reduced fraction $u/v$ with $v>0$, $\gcd(v,N)=1$, $u\equiv tv\bmod N$ and $0<|uv|\le N$. The collector gets **every** such fraction. A dependency is a nonempty subset of distinct target-indexed rows in which every prime outside the base has even total exponent in the product of the absolute heights. This necessary condition allows either sign and either parity of the number of rows. Repeating a row twice within one target contributes the zero vector and is removed.

| Target supply | Upper bound on the probability of any dependency with $(\log M)^{O(1)}$ draws | Larger budgets for which success still tends to zero |
|---|---|---|
| Independent uniform unit squares | $\exp[-(1/\sqrt{2}-o(1))L]$ | $\exp(aL)$ targets for every fixed $a<1/2$ |
| Independent pairs of consecutive unit squares from the conic sampler | $\exp[-(1/2-o(1))L]$ | $\exp(aL)$ pairs for every fixed $a<1/(2\sqrt{2})$ |

Both probabilities average over the choice of $N$ and the draws at that $N$. In either model, the average bound precludes inverse-polynomial success in the bit length on **every** balanced semiprime. Finding rows and running faster linear algebra cannot raise the probability that a qualifying subset exists. Other target distributions, heights and sources of modular relations require separate analysis.

## Counting over the modulus

Choose constants $c,d>0$ with $cd<1/2$, and set

$$Y=\exp(cL),\qquad \sigma=1-d\sqrt{\log\log M/\log M}.$$

Enlarge the base to $S=\mathcal B_N\cup\{\ell\text{ prime}:\ell\le Y\}$. Call an integer eligible if every prime outside $S$ has even exponent in its absolute value.

1. **Eligible fractions.** For $W\ge1$, the number of eligible positive integers up to $W$ is at most $W^\sigma Z$, where $Z=\zeta(2\sigma)\prod_{\ell\in S}(1+\ell^{-\sigma})$. This follows by writing the integer as a square times a squarefree product of primes in $S$. Uniformly in the original base, $\log Z=o(L)$. Indeed, $\sum_{\ell\le Y}\ell^{-\sigma}\le Y^{1-\sigma}\sum_{\ell\le Y}1/\ell=O((\log M)^{cd}\log\log Y)$ by the prime harmonic sum bound, and base primes above $Y$ add at most $(\log M)^A Y^{-\sigma}$. In a primitive fraction, eligibility of $|uv|$ implies eligibility of both coordinates. Rankin's bound applied to the two coordinates together therefore counts at most $2M^\sigma Z^2$ signed eligible fractions of height at most $M$. Each fixes at most one unit target modulo $N$. Any target law with point probabilities $O(1/M)$ hits one with probability at most $\exp[-(d-o(1))L]$.
2. **Repeated primes within one target.** Over all signed fractions of height at most $M$, at most $O((M/\ell)(1+\log M))$ have $\ell\mid uv$. To see this, divide the divisible coordinate by $\ell$ and sum the number of pairs of positive integers of product at most $M/\ell$. Thus the number of ordered pairs of fractions sharing any prime $\ell>Y$ is $O(M^2(1+\log M)^2/Y)$. If two distinct normalized fractions $u/v,w/h$ represent the same target modulo $N$, then $N\mid uh-wv$. This nonzero integer has absolute value at most $2M^2$. It has only $O(1)$ prime divisors at least $X$ and therefore allows only $O(1)$ candidate moduli $pq$. For each modulus the first fraction specifies the target. There are $\Theta(M/(\log M)^2)$ possible moduli by the prime number theorem, and each conditional target point has mass $O(1/M)$ for uniform unit squares. The average probability of any such repeated prime is $O((1+\log M)^4/Y)$.
3. **Repeated primes across independent targets.** For each fixed $N$ a prime $\ell$ occurs in some fraction of a single target with probability $O((1+\log M)/\ell)$ by the row count and the point-probability bound. Multiply the bounds for two independent draws and sum over $\ell>Y$. The probability for a pair of draws is $O((1+\log M)^2/Y)$. This includes repeated target values. [Earlier fixed-modulus cross-target bound](short-rational-relations.md#which-cancellations-can-move-between-independent-targets)

If a selected row is ineligible, it contains a prime outside $S$ to odd exponent. A dependency requires another selected row containing that prime. Consequently an eligible row, a within-target repeat, or a cross-target repeat must occur. There is no union bound over subsets. For $m$ independent targets the result is at most

$$m\exp[-(d-o(1))L]+O\big((m+m^2)(1+\log M)^4/Y\big).$$

For polynomially many draws, take $c=d$ arbitrarily close to $1/\sqrt{2}$ from below. For $m\le\exp(aL)$, choose $d>a$ and $c>2a$ with $cd<1/2$. This is possible for every fixed $a<1/2$. The same argument applies to independent uniform units. It also transfers to a joint target law within total variation distance $\epsilon$ of such draws with an additive $\epsilon$ loss. [Labelled dictionary sampler and its distance bound](short-rational-relations.md#generating-near-uniform-targets-together-with-labels)

## What survives for consecutive square targets?

For each fresh accepted unit $s$, form $x=(s-s^{-1})/2$, $y=(s+s^{-1})/2$ and keep unit $x,y$. Collect all short rows for $t=x^2$ and $t+1=y^2$. Each marginal has point probabilities $O(1/M)$ by the [conic fiber count](short-rational-relations.md#a-coverage-bound-for-any-one-fraction-from-each-target). Distinct sampled pairs are independent conditional on $N$.

Group rows from one draw by their rational value after subtracting the target offset, taking $u/v$ for a row of $t$ and $w/h-1$ for a row of $t+1$. Each group has at most two members. For a mixed pair of rows, $N$ divides $E=(u+v)h-wv$. When $E\ne0$, its absolute value is at most $3M^2$, so the repeated-prime count over possible moduli in step 2 still applies. The same-target determinant handles rows with equal offsets. With high probability, every prime above $Y$ that occurs in multiple rows appears only within one group of one draw. Cross-draw repeats cost at most four times the bound in step 3 per pair of draws.

When $E=0$, uniqueness of the reduced fraction with positive denominator gives $h=v$, $w=u+v$. Away from repeated primes between groups, each selected group must cancel its own odd exponents outside $S$. A singleton must be eligible. In a doubleton the product of the two absolute heights is $|u(u+v)|v^2$. Since $\gcd(u,u+v)=1$, both $|u|$ and $|u+v|$ must be eligible. The two starting rows give $v\max(|u|,|u+v|)\le M$.

If both magnitudes are at most $\sqrt M$, there are $O(M^\sigma Z^2)$ choices for the two signed integers. Otherwise put the larger magnitude in a dyadic band $[W,2W)$ starting at $W\ge\sqrt M$. Choose the larger coordinate, then $1\le v\le M/W$. This gives $O(MW^{\sigma-1}Z)$ choices per band. Each choice fixes $t$ whenever the denominator is a unit. The point-probability bound and $O(\log M)$ bands give probability at most $\exp[-(d/2-o(1))L]$ per draw for an eligible doubleton. Thus the pooled estimate holds with $d/2$ in the single-draw exponent. Take $d=2c$ with $c<1/2$ for the polynomial-budget bound. Taking $c>2a$, $d>2a$ and $cd<1/2$ gives the stated larger-budget result. This covers subsets spanning many rows of both members of many independent pairs.

**Nonunit checks within the height bound.** For a unit target, if a short unreduced congruent row exposes $p\mid v$, then $p\mid u$ as well. Writing $u=pu'$, $v=pv'$ gives $0<|u'v'|\le q/p\le 2$. There are only constantly many possible fractions, each fixing the target modulo $q$. The probability is $O(1/X)$ per draw under either marginal bound, including the symmetric case with $q\mid v$. Conic or unit sampling itself also encounters a proper factor with probability $O(1/X)$ per attempt on these inputs. These bounds cover gcds of row coordinates and sampler inputs.

**Increasing the height.** The dependency exponents in the table also hold with height cap $N\exp(o(L))$, including $N$ times any fixed power of the bit length. In the row counts use $H=M\exp(o(L))$ as the bounding height, keeping the modulus count and target probabilities unchanged. The fraction counts gain at most $\exp(o(L))$, the determinant is $O(H^2)$ and still has only constantly many prime divisors above $X$, and the dyadic split for an exact pair moves to $\sqrt H$. All these changes fit in the stated error terms.

## Finite check

For all six products of distinct primes from $\{11,13,17,19\}$, there are 289 unit-square targets. With base $\{2,3\}$, 254 admit an odd dependency and 234 have an eligible singleton at the original base. Thus 20 require multiple rows. Enumerating by target and independently enumerating by integer determinant give identical sets of repeated-prime fraction pairs, with 961, 171 and 65 pairs at thresholds 3, 7 and 13. These counts verify the incidence argument on small inputs. The asymptotic conclusion comes from the bounds above. [Exhaustive check](check_same_target_determinants.py)

There are also 44 possible consecutive-square targets on these moduli. At $Y=3$, the check enumerated 106,906 nonempty subsets avoiding all three events (eligible singleton, eligible exact pair, repeat between distinct groups). Every subset retained at least one outside-base prime to odd exponent. The exact pair is necessary in the decomposition. For $N=209$ and $t=168$, the rows $4/5$ and $9/5$ both retain the prime 5, and their product is the square $36/25$.
