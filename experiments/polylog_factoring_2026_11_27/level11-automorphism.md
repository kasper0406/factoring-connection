# Denominators of units in a sampled ideal's left order

Given \(n\) coprime to 11, sample a uniform projective line and its right ideal \(I\) of norm \(n\) in the level-11 order \(\mathcal O\). Compute a noncentral unit \(u\) of its left order \(\mathcal O_L(I)\). The entries of \(nu\) in the original basis are integers. Return \(\gcd(n,(nu)_1,\ldots,(nu)_4)\) when it is a proper divisor. The target order determines which units are available. On \(n=pq\), a uniform trial succeeds with probability at most \(2/(p+1)+2/(q+1)\). On some infinite families its success probability is zero. The computation and verification both use the input integer alone. [Public probe](level11_automorphism.py) · [Earlier level-11 sampler](level11-sampling.md)

To compute the units, form the column HNF of the products \(a\bar b\) over basis elements \(a,b\) of \(I\). This gives \(I\bar I=n\mathcal O_L(I)\), as seen locally by writing \(I=\alpha\mathcal O\). Enumerate all its points of norm \(n^2\) with the fixed-dimensional lattice routine. The column covolume is \(n^4\), since all maximal orders of discriminant 11 have the same covolume. The squared shortest length is \(\gg n^2\) since the reduced norm is integral on the left order. The LLL enumeration box therefore has bounded size as in the norm sampler. There are four or six points. In the four-unit class every noncentral unit is \(\pm i\). In the six-unit class every noncentral unit is a sign or an integral scalar translate of a primitive cube root of unity. Thus all eligible choices give the same integrality divisor on that line. [Brandt matrices and their ideal classes](https://jvoight.github.io/quat-book.pdf) · [Level-11 class counts](../../notes/nativemod_results.md#L466-L502)

## Exact split probability

Let \(p,q\) be distinct primes away from 11 and \(n=pq\). Label the four-unit source class 0 and the six-unit class 1. Write \(B_s(0,j)\) for the number of lines at \(s\) whose ideal reaches class \(j\). These sum to \(s+1\). The earlier theta calculation gives \(B_s(0,0)=(3(s+1)+2a_s)/5\) where \(a_s\) is the coefficient of the level-11 weight-two newform. Define

\[
t_0(r)=1+\left(\frac{-4}{r}\right),\qquad
t_1(r)=1+\left(\frac{-3}{r}\right)
\]

using Kronecker symbols. These are the numbers of projective lines fixed by a noncentral unit of each respective order. They count eigendirections in the local two-dimensional space. At \(r=2\) a Gaussian unit is nonscalar because \(i-1\) has norm 2, so it has one fixed line. A primitive cube root has none there. At \(r=3\) the latter has one. This is the local horizontal-isogeny criterion for the two quadratic orientations. [Colò and Kohel, §2](https://doi.org/10.1515/jmc-2019-0034)

Write \(E_p\) for the event that the chosen target unit is integral in the original order at \(p\). Delete the \(p\)-step of the sampled ideal, obtaining its intermediate ideal \(I_q\). Its left order equals the original order at \(p\) and the final left order at all primes away from \(p\). Hence \(E_p\) occurs precisely when the intermediate order has a noncentral unit preserving the final \(p\)-line. Indeed a preserving unit also belongs to the target order, and all noncentral units there agree on integrality. Conversely, a target unit integral at \(p\) already belongs to the intermediate order globally. These two classes match when preservation occurs, since the two unit groups have orders four and six and neither contains the other.

Condition on the \(q\)-line first. At an intermediate order of class \(j\), exactly \(t_j(p)\) out of its \(p+1\) choices work. Thus the number of lines with \(E_p\) and target class \(j\) is \(B_q(0,j)t_j(p)\), yielding \(\Pr(E_p)\le2/(p+1)\). If both \(E_p\) and \(E_q\) occur, the target unit belongs to the original source order globally. The original Gaussian unit must stabilize both lines, giving exactly \(t_0(p)t_0(q)\) possibilities. A proper gcd occurs on the symmetric difference. Its exact count out of \((p+1)(q+1)\) is

\[
\sum_{j=0}^{1}\bigl[B_q(0,j)t_j(p)+B_p(0,j)t_j(q)\bigr]-2t_0(p)t_0(q).
\]

For \(p,q\equiv11\bmod12\), with neither prime equal to 11, all four local line counts vanish. This conclusion holds for every choice of the sampled line.

## Other endomorphisms chosen after seeing the line

Allow any noncentral element \(\beta\) of the target order of norm at most \(K\). If it is integral in the source at \(p\), it belongs to the intermediate order globally and preserves the \(p\)-line. The element \(2\beta-\operatorname{tr}(\beta)\) has trace zero and norm at most \(4K\). The intermediate order's rank-three trace-zero lattice has minimum nonzero norm at least one. Packing bounds the number of these elements by \(O((1+\sqrt K)^3)\), uniformly in the ideal. If \(\beta\) reduces to a scalar at \(p\), its nonzero discriminant \(\operatorname{tr}(\beta)^2-4\operatorname{nrd}(\beta)\) is divisible by \(p^2\). This is impossible for \(p^2>4K\). Each remaining direction fixes at most two lines. For odd \(p\), the trace-zero part determines the same eigendirections as \(\beta\).

Consequently, even inspecting every such target endomorphism adaptively has probability \(O((1+\sqrt K)^3/(p+1))\) of local integrality at \(p\) when \(p^2>4K\). In particular polynomial norm bounds in the bit length give exponentially small coverage at a balanced prime. For an intermediate menu of \(M\) elements already specified before the last line is drawn, the same count is at most \(2M/(p+1)\) provided all reductions are nonscalar. These bounds concern the indicated denominator test and degree constraints.

## Checks

The diagnostic computed the actual units and gcds for all 2,234 lines at 16 semiprimes. All 224 proper gcds matched the counts by prime and target class predicted independently from the eta product. It also verified the action of representative units on an integral basis of each ideal class. At \(n=65\), 24 of 84 lines gave a factor. At \(n=119\), 24 of 144 did. At \(n=23\cdot47=1081\), all 1,152 lines failed, as the inert-prime criterion requires. The fixture primes only enumerate the complete test population.

In 96 separately seeded trials given only the input integer, the routine returned 16 verified proper divisors. For a concrete six-unit example, seed 2 at \(n=91\) chose the public line \((30,29)\). One target unit had numerator \((42,7,-14,-42)\) over 91, of norm \(91^2\). Its coordinate gcd with 91 gave 7. [Reproducible diagnostic](check_level11_automorphism.py)
