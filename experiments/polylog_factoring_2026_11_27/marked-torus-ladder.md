# What a public cubic character buys in an order test

At $N=143=11\cdot13$ the two local norm-one groups have order 12. The public cubic Euler power sends every point to 1, so its three gcds never split this input. The full sequence of successive cubes splits **all 96** points with nonzero cubic Jacobi mark, out of 144 points total. Below is the exact law for both tests. It identifies the group-order promise behind the improvement and bounds how often the simplest version applies. [Implementation](marked_torus_ladder.py) · [Exhaustive checks](check_marked_torus_ladder.py) · [Earlier Lucas exponent analysis](../../notes/primality_analysis.md)

Let $N=pq$ for distinct primes greater than three. Write $\epsilon_r=1$ for $r=1\bmod 3$ and $-1$ otherwise. The norm-one subgroup $T_r$ in $\mathbf F_r[\omega]$ is cyclic of order $s_r=r-\epsilon_r$. Let $z$ be uniform in $T_N=T_p\times T_q$. Its hidden mark $u_r\in\mathbf F_3$ is given by $z_r^{s_r/3}=\omega^{u_r}$. Its public mark $c=-\epsilon_p u_p-\epsilon_q u_q$ is computed in polynomial time by cubic reciprocity. In particular, when $N=2\bmod3$, $c\ne0$ certifies that the hidden marks differ. To get a factor we need a public element with different *actual* local endpoints. At any query $w$, we compare its two Eisenstein coefficients with those of each global root $1,\omega,\omega^2$, taking a gcd of both differences with $N$. [Mark identity, sampling and explicit reciprocity algorithm](eisenstein-torus-labels.md)

## One Euler power

Put $m_r=s_r/3$, $g=\gcd(m_p,m_q)$, $e=(N-\epsilon_p\epsilon_q)/3$ and $h_r=g/m_r$. The sign product is public from $N\bmod3$. The probability of a proper gcd after the single power $z^e$, **conditional on any of the three public marks**, is

$$h_p+h_q-(1+\sigma)h_ph_q,\qquad
\sigma=\begin{cases}1&3\mid e/g,\\1/3&3\nmid e/g.\end{cases}$$

Indeed $e=\epsilon_qm_p+\epsilon_pm_q+3m_pm_q$ and $\gcd(e,m_r)=g$. In a local generator whose $m_r$-th power is $\omega$, an input reaches a cube root of unity precisely when its exponent is $(m_r/g)v_r$, an event of probability $h_r$. On a hit $v_r\bmod3$ is uniform, the hidden mark is $(m_r/g)v_r$ and the actual endpoint has exponent $(e/g)v_r$. The other prime's uniform mark preserves the individual hit probability on each public class. On two hits the public mark is still uniform, since $m_p/g,m_q/g$ are coprime. If $3\nmid e/g$, the diagonal $v_p=v_q$ maps bijectively onto the three public classes because $\epsilon_p m_p/g+\epsilon_q m_q/g\equiv\epsilon_p\epsilon_q e/g\bmod3$. This proves the formula. Conditioning a uniform draw on the certified disagreement cannot improve this particular test.

## The full sequence of cubes

Write $s_r=3^{a_r}b_r$ with $3\nmid b_r$, $G=\gcd(b_p,b_q)$ and $C_r=b_r/G$. For an integer $B\ge1$ let $M$ be the least common multiple of all positive integers at most $B$ coprime to 3. Let $D$ be $M$ times the 3-free part of $N-\epsilon_p\epsilon_q$. Both quantities are public. Test $z^D,z^{3D},z^{9D},\ldots$ against the three global roots. After $O(\log N)$ terms the local 3-parts have vanished. Computing and using $D$ costs polynomial time in $B$ and $\log N$.

The probability that a component ever reaches the three roots is

$$\kappa_r=\frac{\gcd(C_r,M)}{C_r}.$$

Indeed $\gcd(b_r,N-\epsilon_p\epsilon_q)=G$, by reduction modulo $b_r$. Uniformity on the prime-to-three component gives the formula. These two events are independent of each other and of both hidden marks. If only one occurs, a gcd splits the input. When both occur, a split fails precisely when the two 3-primary components have the same order $3^j$ and the same last nonidentity cubic root. Give the identity the tag $(0,0)$.

For nonzero $u_r$ the tag is deterministically $(a_r,D b_r^{-1}u_r\bmod3)$. For zero $u_r$ its law assigns probability $3^{1-a_r}$ to $(0,0)$, and $3^{j-a_r}$ to each $(j,1),(j,2)$ for $1\le j<a_r$. This follows by counting elements of each order in the index-three subgroup. Given $c$, average the collision probability of these independent tag laws over the three equiprobable pairs $(u_p,u_q)$ with $-\epsilon_pu_p-\epsilon_qu_q=c$. Calling the result $\sigma_c$, the success probability is exactly

$$\kappa_p+\kappa_q-(1+\sigma_c)\kappa_p\kappa_q.$$

Without mark selection the tag collision probability is $(9^{\min(a_p,a_q)}+3)/(4\cdot3^{a_p+a_q})\le1/3$, by summing the two tags of every positive height. In particular if $C_p\mid M$ or $C_q\mid M$, a uniform point gives a factor with probability at least $2/3$. This includes $\min(C_p,C_q)\le B$ and permits an arbitrarily large common part $G$.

There is a sharper public filter when $a_p=a_q=a$ and $3^{a+1}\mid N-\epsilon_p\epsilon_q$. Then $\epsilon_p b_p+\epsilon_q b_q=0\bmod3$. A zero and a nonzero hidden mark have different heights. When both marks are nonzero, equal terminal roots require $u_p/b_p=u_q/b_q\bmod3$, which forces $c=0$. Thus $\sigma_c=0$ for either nonzero class. Under the same cofactor promise **every point with nonzero mark splits**. The promised orders are an analytic condition on the input. The gcd algorithm itself uses only $N,z,B$. The example $N=143$ has $a_p=a_q=1$ and $C_p=C_q=1$. Conversely $z=\omega$ can have different hidden marks while every common power is a global root. Selecting this particular point is distinct from drawing uniformly in a mark class.

## Which correlations can change the cofactor cost?

Take independent uniform torus points $z_1,\ldots,z_k$ and choose a word $\prod_j z_j^{v_j}$ from their public marks. At either prime the unresolved prime-to-three group after the power $D$ has size $d_r=C_r/\gcd(C_r,M)$. Conditional on the entire mark transcript its $k$ input coordinates are still independent uniform elements of a cyclic group of order $d_r$. The probability that the word kills this component is exactly

$$\gcd(d_r,v_1,\ldots,v_k)/d_r.$$

Thus a menu of $T$ primitive exponent vectors (entries of gcd one) chosen from the marks has split probability at most $T(\kappa_p+\kappa_q)$. Asking also for the marks of multiplicative combinations adds only linear combinations of the known marks. Nonprimitive exponents can improve the hit rate by the displayed gcd, recovering ordinary group-order annihilation. This statement leaves room for operations using other information about the integer representatives.

For $B=1$ and distinct primes sampled uniformly in $[R,2R]$,

$$\mathbb E(1/C_p+1/C_q)\ll (\log R)^4/R.$$

To see this, fix $a_p,a_q$ and the two signs. Set $X=R/3^{a_p}$, $Y=R/3^{a_q}$. The $b_p$ lie between constant multiples of $X$ and the $b_q$ between constant multiples of $Y$. Including all integers in these ranges only increases the sum. Expanding $\gcd(i,j)=\sum_{d\mid i,j}\varphi(d)$ bounds the sum of gcds by $O(XY\log R)$. Dividing by $b_p\gg X$ gives $O(Y\log R)$. Sum geometrically over $a_q$ and over the $O(\log R)$ possible $a_p$, then normalize by $\gg R^2/(\log R)^2$ prime pairs. Symmetry gives the other term. In particular the basic ladder with polylogarithmically many marked primitive words has vanishing average coverage. The simpler promise $\min C_r\le B$ itself holds on at most $O(B(\log R)^4/R)$ of the pairs. The latter count does not bound the larger set $C_r\mid M$ for growing $B$.

The unconditional small-menu count concerns this torus, its common-power tests and the specified information available for choosing words. Prime-transcript arguments in a generic ring differ from a cubic symbol computed on explicit integer coordinates. The symbol supplies a comparison of local powers. Computing the CRT element representing those powers on ordinary inputs would meet the earlier mixed-root extraction criterion immediately. [Pure-ring model](../../notes/factorbits_results.md) · [Earlier group-order sufficient conditions](../../notes/suffmap_results.md) · [Extraction with retained witnesses](shanks-witness-fusion.md)

## How a cubic endpoint gives the original square-root witness

There is a direct conversion to the object in the repo's window reductions. Suppose we hold an actual $h\in T_N$ with $h^3=1$ and different local labels. For each $j=0,1,2$ compute

$$I_j=\big(1+\omega^{-j}h+(\omega^{-j}h)^2\big)/3 \pmod N.$$

At either prime $I_j$ is 1 precisely when the label equals $j$, and otherwise it is zero. Hence $I_j$ has zero $\omega$ coefficient modulo $N$ and lies in $\mathbf Z/N$. At each of the two occupied labels $t_j=2I_j-1$ is a mixed square root of 1. Conversely, from any such square root $t$ put $I=(1+t)/2$ and $h=1+(\omega-1)I$. Its local values are $1$ and $\omega$. These conversions use a constant number of ring operations. For $N=143$ and $z=1+66\omega$, the public cubic mark is 1 and $z^{16}=12+77\omega$. Here $I_0=78$ and $t_0=12$, with $12^2=1\bmod143$ and $\gcd(12-1,143)=11$. The mark alone is a single element of $\mathbf F_3$. The projector takes the actual ring element as its argument. [Fixed-window reduction](fixed-window-split.md) · [Weak quarter-window predictor](quarter-fourier.md)
