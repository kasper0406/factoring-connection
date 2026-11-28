# How many low factor pairs survive their product and bitwise AND?

Let $P,Q$ be independent uniform odd residues modulo $2^k$. Observe the low product $n=PQ\bmod 2^k$ and the bitwise AND $A=P\mathbin{\&}Q$. An ordered prefix pair is uniquely specified by its first coordinate, so the exact number of compatible pairs is

$$
 F_k(n,A)=\#\{x\ {\rm odd}\bmod 2^k:
            x\mathbin{\&}(n/x\bmod 2^k)=A\}.
$$

Write $H_2(t)=-t\log_2t-(1-t)\log_2(1-t)$ and

$$
 \gamma=1-H_2(1/4)=\tfrac34\log_2 3-1=0.18872187554\ldots.
$$

**The typical fiber has $2^{(\gamma+o(1))k}$ ordered pairs.** More precisely, there are absolute constants $c,C>0$ such that for all sufficiently large $k$,

$$
 \Pr\!\left(\left|\log_2 F_k(n,A)-\gamma k\right|>Ck^{8/9}\right)
 \le \exp(-c k^{1/9}).
$$

The probability weights observations by how often an actual random pair produces them. The earlier entropy count proved the lower exponent. The upper bound below uses short blocks in which the inverse-product map is affine. It also bounds the work of enumerating the entire prefix tree from the supplied bits. [Earlier lower bound and binomial interpretation](adaptive-binomial-lift.md#how-much-does-the-fixed-singleton-menu-leave-ambiguous)

## An exact average for a block

Fix a block length $b$, write $M=2^b$, and take an odd residue $h\bmod M$. For unsigned $b$-bit words $r,s$, let $f_h(r,s)$ count the pairs $r',s'$ with

$$
 s'-h r'=s-h r\pmod M,\qquad r'\mathbin{\&}s'=r\mathbin{\&}s .
$$

For $0\le d<b$, let $z_d$ be the number of zero bits of $r\mathbin{\&}s$ at positions strictly above $d$. With $r,s$ fixed and $h$ uniform over odd residues, the exact average is

$$
 \mathbb E_h f_h(r,s)=1+
 \sum_{\substack{0\le d<b\\r_d\ne s_d}}
      3^{z_d}2^{-(b-d-1)}. \tag{1}
$$

The identity pair contributes one. For an alternative, a possible odd slope requires that $r'-r$ and $s'-s$ have the same 2-adic valuation $d<b$. Thus all bits below $d$ are fixed and both bits at $d$ flip. This preserves their AND exactly when $r_d\ne s_d$. Above $d$ there are three possible ordered bit pairs for each zero of the given AND and one for each one. For each resulting alternative pair, division of the slope congruence by $2^d$ prescribes precisely one odd residue of $h\bmod 2^{b-d}$. Its probability is $2^{-(b-d-1)}$, including probability one at $d=b-1$.

Now choose $r,s$ independently and uniformly. In a suffix of length $\ell=b-d-1$, the number $z_d$ of zero AND bits is binomial with parameter $3/4$. Hoeffding's inequality and a union bound over all suffixes give

$$
 z_d\le\tfrac34(b-d-1)+b^{2/3}\quad\hbox{for every }d
$$

except with probability at most $b\exp(-2b^{1/3})$. The empty suffix always satisfies the bound. On this event, (1) gives

$$
 \mathbb E_h f_h(r,s)\le 1+b\,2^{\gamma b+(\log_2 3)b^{2/3}}.
$$

Markov's inequality now gives, for independent uniform $r,s$ and uniform odd $h$, and sufficiently large $b$,

$$
 \Pr\!\left(\log_2 f_h(r,s)>\gamma b+3b^{2/3}\right)
 \le \exp(-c_0b^{1/3}) \tag{2}
$$

with an absolute $c_0>0$. Keeping the observed AND fixed in (1) retains the typical three choices per zero. Averaging the number of choices over the AND first gives the weaker exponent $\log_2(5/4)$.

## From block counts to the whole fiber

Use $b=\lfloor k^{1/3}\rfloor$ for sufficiently large $k$. Read the first $b$ bits together and split the remaining positions into full $b$-bit blocks and one possible shorter block. On a compatible prefix of length $j\ge b$ write $x,y$ for its least nonnegative coordinates. Appending digits $u,v\bmod 2^b$ gives the exact constraint

$$
 (x+2^j u)(y+2^j v)\equiv n\pmod{2^{j+b}}
 \quad\Longleftrightarrow\quad
 v\equiv x^{-1}\frac{n-xy}{2^j}-(y/x)u\pmod{2^b}.
$$

On the path of the actual $P,Q$, the slope $h=-Q/P\bmod 2^b$ depends only on the first $b$ bits and is uniform over the odd residues. In each later full block, the actual digits $(r,s)$ are independent uniform $b$-bit words, independently of $h$. Since the actual child lies on the line, its intercept is $s-hr$. Thus the number of children respecting that block's product and AND equals $f_h(r,s)$. Every full block on the true path satisfies the upper bound in (2) except on a combined event of probability at most $(k/b)\exp(-c_0 b^{1/3})$. The union bound permits the same slope to occur in all blocks.

For each fixed observation $(n,A)$, give the tree root mass one and divide each node's mass equally among its viable children at that block. A node with no children contributes zero mass to the leaves. The first and last blocks have at most $2^b$ children apiece, since $x$ determines $y$. For each actual path outside the exceptional event, minus the logarithm of its leaf mass is therefore at most

$$
 2b+(k/b)(\gamma b+3b^{2/3})\le B_k,\qquad
 B_k=\gamma k+C_1 k^{8/9}.
$$

The total mass of the leaves is at most one. Consequently each observation has at most $2^{B_k}$ leaves of mass at least $2^{-B_k}$. To bound **all** leaves, use the input law once more. Given its observation, the actual pair is uniform on the fiber. If a fiber exceeds $2^{B_k+1}$ leaves, the conditional probability of a leaf lighter than $2^{-B_k}$ is greater than $1/2$. The probability of such large fibers is at most twice the exceptional path probability, hence at most $\exp(-c k^{1/9})$ for a smaller absolute $c$.

For the matching lower bound, the $w=k-1$ nonunit bits of the AND are independent Bernoulli variables of parameter $1/4$. Except with probability $\exp(-\Omega(w^{1/3}))$, their total weight differs from $w/4$ by at most $w^{2/3}$. The number of such AND words is at most $2^{w H_2(1/4)+O(w^{2/3})}$ by the binomial entropy bound. There are $2^w$ product residues and $2^{2w}$ equiprobable input pairs. Summing the mass of fibers smaller than $2^{\gamma w-C_2 w^{2/3}}$ over these observations gives probability $\exp(-\Omega(w^{2/3}))$ for sufficiently large $C_2$. This proves the claimed lower tail together with the exceptional AND probability.

## Enumeration and the balanced-prime model

The usual one-bit prefix enumerator can use the upper bound without implementing the block decomposition. From an odd $x$ at depth $j$, try $x$ and $x+2^j$ modulo $2^{j+1}$, compute the odd partner by inversion and check the new AND bit. The live list contains exactly $F_j$ elements. For $j<\sqrt k$ its size is at most $2^{\sqrt k}$. At every larger $j$ the marginal observation has the law of a fresh $j$-bit instance. A union bound puts every such list below $2^{\gamma k+O(k^{8/9})}$ except with probability $\exp(-\Omega(k^{1/18}))$. Thus enumeration takes $2^{(\gamma+o(1))k}$ bit operations with probability tending to one. It is exact for every input, including inputs that take longer.

The bounds also hold when $P,Q$ are residues of distinct primes chosen uniformly in $[R,2R]$, with $k\to\infty$ and $2^k\le R^{1-\delta}$ for fixed $\delta>0$. Brun–Titchmarsh bounds the number of primes in each odd residue class by $O_\delta(R/(2^k\log R))$. By the prime number theorem, the probability of each residue pair is at most a constant depending on $\delta$ times its uniform-model probability, also after conditioning on distinctness. Since the full fiber and the sizes of each enumerated prefix list are functions of the residue pair, their exceptional probabilities transfer. One may order the primes afterward by symmetry. If $k=\lfloor\theta\log_2(PQ)\rfloor$ for fixed $\theta<1/2$, only a bounded number of depths occur within this dyadic sample and the same union bound applies. [Montgomery and Vaughan, *The large sieve*](https://doi.org/10.1112/S0025579300004708)

For an odd balanced semiprime $N=pq$ with its actual $k$ low AND bits supplied, an explicit residue list of size $S$ can be tested using $\widetilde O(S(1+N^{1/4}/2^k))$ bit operations by affine divisor recovery. On the density-one set above, $k=\lfloor(2/9)\log_2N\rfloor$ gives time $N^{1/4-(2/9)H_2(1/4)+o(1)}=N^{0.06971597\ldots+o(1)}$. At $k=\lfloor(1/4)\log_2N\rfloor$ the exponent is $\gamma/4=0.04718047\ldots$. These statements give runtime on a set of relative density tending to one. The estimate assumes the true singleton answers are supplied, and leaves the cost of evaluating them from $N$ open. [Affine divisor recovery](coupled-rank-windows.md) · [Exact block checks](check_shared_and_fiber_entropy.py)
