# Factoring from a list of quadratic-form classes

**Research decision.** The norm-$N$ bisector is a useful fast path. For growing nonsquare cutoffs, test whether public data can select a short list of quadratic-form classes cheaply. There is now an executable verifier for such a list and an end-to-end baseline that constructs it exhaustively. The baseline still pays for each possible minimum in the requested interval.

## A prefix decision sufficient for factoring

Use an eligible public branch $(N,\ell,c,\Delta)$ with $N=pq$, $p<q\le2p$, $c^2+\Delta\ell^2=4N$, $\gcd(N,\ell\Delta)=1$ and $\Delta>4$. In the [repository's four-line reduction](../../notes/shared_collision_nov26.md#4-public-ellipse-count-and-squared-class-decoding), let $C$ be the proper ideal class of a prime over $p$, up to inversion. The public principal element of norm $N$ pairs the primes over $p$ with those over $q$. Thus the unordered prime classes at either rational prime are $\{C,C^{-1}\}$. The first mixed norm divided by $N$ is the minimum $s$ of the reduced square $C^2$; its inverse has the same minimum. This also works for nonmaximal orders under the stated good-prime hypothesis.

For a public interval $I$ of possible minima inside the usual cap, define
$$\mathcal L(I)=\{F\text{ a primitive reduced class of discriminant }-\Delta:
                       \min(\mathrm{red}(F^2))\in I\}/(F\sim F^{-1}).$$

A class on this list represents a proper prime divisor of $N$ **if and only if** $s\in I$. In one direction, $C$ represents $p$ and is on the list when $s\in I$. In the other direction, a representation of $p$ or $q$ corresponds to one of the two invertible ideals of that prime norm, so its class is $C$ or $C^{-1}$. This gives a factor-or-empty answer on the whole prefix $I=[1,B]$ and on intervals. It uses the existence of a mixed point in that prefix. The numeric signed count modulo 16 is a different output. [Forms and ideal classes](https://math.mit.edu/~drew/sutherland-phd.pdf)

The simple rational-anchor extension of reflection would seek a point of norm $j$ in the field and multiply it by the public point. Clearing a denominator $h$ requires $u^2+\Delta v^2=4j h^2$. For nonsquare $j$, $v\ne0$, so $h\ge\sqrt{\Delta/j}/2$. In particular an integral multiplier is unavailable below $\Delta/4$. This does not exclude cancellation of denominators in another algorithm. The class verifier avoids needing such an anchor.

## Verifying a public list in one batch

Let $F=[A,b,d]$ be a reduced representative. Set $Q=\lfloor\sqrt{2N}\rfloor$. If $F(x,y)=p$ or $q$, reduction gives
$$|x|\le X_F=\lfloor\sqrt{2Q/A}\rfloor,\qquad
|y|\le Y_F=\lfloor\sqrt{2Q/d}\rfloor.$$
Indeed $F(x,y)\ge(Ax^2+d y^2)/2$ by $|b|\le A\le d$. Check axes directly; a prime on an axis requires $A$ or $d$ to equal that prime. On the other points $4AF=(2Ax+by)^2+\Delta y^2$. The unit public root $r=c/\ell\bmod N$ therefore supplies the linear tests
$$i-\rho y,\qquad i=|x|,\ y>0,\qquad
 \rho\in\{\pm(b+r)/(2A),\ \pm(b-r)/(2A)\}\pmod N.$$
The public cap guarantees $A<p$. Any prime representation yields a nonunit in at least one test.

Pool all these rectangles. Let $W=\sum X_FY_F$ and $T=\sum Y_F$, counting the distinct sign routes used for each form. Take long-side block widths of order $\min(X_F,\sqrt W)$ rounded down to powers of two. Build a single polynomial $P_k(z)=\prod_{1\le i\le k}(z-i)$ for each width and share it across all forms. The initial number of evaluations is $O(\sqrt W+T)$. Split nonunit blocks; at a given level each route and row has at most two genuine nonunit positions including padding for all sufficiently large promised inputs, since $2X_F=O(N^{1/4})<p$. The remaining levels each have $O(T)$ evaluations. At a leaf, try both signs and require the ordinary form value to divide $N$ strictly and exactly. This last test excludes incidental modular hits and makes negative prefix results meaningful.

With fast monic product and remainder trees the verifier costs $\widetilde O(K+\sqrt W+T)$ in bit operations for $K$ already supplied classes. For a single class this is $\widetilde O(1+(N/\Delta)^{1/4})$, since $Ad$ is between $\Delta/4$ and $\Delta/3$. No hidden mixed root or supplied prime factor enters the verifier. The bound excludes *building the list*.

For example, at $\Delta\asymp\sqrt N$ and for similarly shaped forms with $A,d\asymp\sqrt\Delta$, the row term is of order $K N^{1/8}$. Beating $N^{1/5}$ by a fixed power in this model requires $K\le N^{3/40-\eta}$ for some $\eta>0$, with list generation also below the target. Skew forms have different row costs, so report actual $W,T$, and setup rather than using a universal limit on $K$. Even a single supplied typical class leaves a power of $N$ in this verifier's bound.

## What the current list builder costs

For each $a\in I$ the [baseline](class_prefix.py) enumerates roots $b^2=-\Delta\bmod 4a$ using SymPy, forms $[a,b,(b^2+\Delta)/(4a)]$, and enforces reduction and primitivity. Factoring the public discriminant reuses the deterministic block routine from the [bisector experiment](bisector.md). The baseline rules out classes missing a locally represented unit square at a prime of $\Delta$. At odd primes the reduced form has rank one modulo that prime and the square character of a nonzero diagonal entry suffices. At 2 it tests the residue 1 modulo 8. It then runs the [existing conic and Gauss-duplication lift](../ellipse_global/README.md#lifting-with-a-supplied-factorization-of-δ) on the retained candidates. PARI reports an obstruction where the rational conic is insoluble. For a square class the earlier unit repair yields a square representation prime to the discriminant, and composing the seed with all ambiguous classes gives the full fiber. Each returned root is checked by duplication.

All of the discriminant work, modular-root enumeration (including any factorization of $4a$ inside SymPy), conics, lift production and deduplication is included in the measured build times. No claim of polylogarithmic list generation follows. In particular, this builder performs at least $|I|$ iterations. The [near-square coverage bound](bisector.md#how-frequently-can-small-cutoffs-work) shows why scanning all small $a$ cannot support an almost-all $N^{1/8}$ factorer on the $\ell=1$ family.

## Which shortcuts deserve a test?

For every class on the list the two known index-$N$ sublattices give back the original form after rescaling. Explicitly, let $t=(c-\ell b)/2$ and
$$M=\begin{pmatrix}t&-\ell d\\ \ell A&t+\ell b\end{pmatrix}.$$
Then $\det M=N$ and $F(Mv)=N F(v)$. Conjugating the public element gives the other sign. Each image is the full corresponding public congruence line, as seen from the index and the public root. The shortest-vector and full theta-series data on those unmixed lines therefore reproduce the candidate's own data. Ranking by the transformation itself supplies no new observation about the hidden prime class. Similarly, order-two class characters are trivial on every class square.

A measurable target is to use the public middle coefficient $c$ *before* listing the classes. On the $\ell=1$ branch the true class is represented up to inversion by $[p,c,q]$. Equivalently a true candidate $F=[A,b,d]$ has a unimodular basis $(x,y),(z,t)$ with
$$2Axz+b(xt+yz)+2dyt=c,\qquad xt-yz=\pm1,$$
and its diagonal norms lie in the balanced factor interval. The determinant identity then forces their product to be $N$.

There is an immediate filter on proposed selection methods. **The frame condition modulo any modulus is soluble for every primitive candidate.** At each prime power choose a vector $v$ of unit norm among $(1,0),(0,1),(1,1)$ and a complementary vector $w$ of determinant 1. Replacing $w$ by $w+k v$ shifts the middle coefficient by $2k F(v)$. At odd primes this attains every residue. At 2 the middle coefficient already has the required parity $b\equiv c\pmod 2$, so it attains every admissible residue there as well. Combine by CRT. Independently, the real balanced-frame condition is soluble for every candidate: the positive Gram matrix with diagonal entries $\sqrt N$ and off-diagonal $c/2$ has the same determinant as $F$, so an area-preserving real change of basis gives that matrix. Consequently, separate modular solubility and continuous feasibility cannot prune this list. The finite tests construct local frames at 8, 32, 9, 25 and 49 explicitly.

A range test would have to detect or exclude such a basis for **many candidates at once**. A concrete success condition is retaining the true class on unseen inputs with total build and verification work below the chosen benchmark. Charge empty intervals and failed attempts. This is a proposed selection mechanism; no cheap range test has been established.

## Controls

* On 194 repository branches and ten deliberately selected tie cases, the empty prefix through $s-1$, the prefix through $s$, the singleton interval $[s,s]$, and the full cap all gave the expected answers. Seventy first minima were nonsquare. The full cap needs no supplied index.
* On all 248 order discriminants with $4<\Delta\le500$, both a full and an upper-range class list agreed with independent enumeration and forward Gauss duplication, covering 1,548 proper classes.
* The original 15 N-only scaling inputs were rerun with the full cap and with an explicitly supplied first index. On N-only inputs the factors were used in the controls to label $s$ by squaring the reduced form $[p,c,q]$. The blind runs received neither $s$ nor the factors. Both public-line similarity identities were also checked on the resulting lists.
* Two additional controls use the first primes after $10^7,1.4\cdot10^7$ and after $10^9,1.4\cdot10^9$. The 47-bit case was run both ways; the 61-bit case is a supplied-index test only.

| $N$ | $s$ | Public cap | Classes given $s$ / blind cap | Search polynomial values given $s$ / blind cap |
|---|---:|---:|---:|---:|
| 6,494,011 | 28 | 47 | 4 / 16 | 315 / 948 |
| 38,942,974,393 | 205 | 358 | 8 / 90 | 2,010 / 19,899 |
| 140,000,556,000,551 | 582 | 1,392 | 16 / 436 | 15,924 / 383,854 |
| 1,400,000,032,800,000,161 | 6,253 | 38,935 | 8 / unmeasured | 13,212 / unmeasured |

These are evaluations during the final search. [Raw results](class_prefix_results.json) separately record setup, root enumeration, conic and class counts, search shape, polynomial construction and wall times. [Reproduction and independent controls](class_prefix_check.py).

Use Python with python-flint, NumPy and SymPy, and PARI/GP on PATH or supplied via the GP environment variable. From the repository root run

    python experiments/ellipse_disk_mod8/class_prefix.py 6494011
    python experiments/ellipse_disk_mod8/class_prefix.py 6494011 --only-minimum 28
    python experiments/ellipse_disk_mod8/class_prefix_check.py --output results-rerun.json
