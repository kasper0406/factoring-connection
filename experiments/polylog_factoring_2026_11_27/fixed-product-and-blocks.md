# What a fixed dyadic product reveals about upper factor bits

Fix **any** odd residue $n\bmod 2^k$. Choose $X$ uniformly over the odd residues and put $Y=n/X\bmod 2^k$. The upper $b$ bits of the two endpoints jointly have total variation distance
$$
O(k^2 2^{b-k/2})
$$
from two independent uniform $b$-bit words. Thus for $b\le k/2-\omega(\log k)$ their AND is close to iid Bernoulli $1/4$, uniformly in the product. If $2b>k-1$, support alone gives a distance of at least $1-2^{k-1-2b}$ for the joint endpoint words. As one consequence, at position $j\ge1$ the probability of a positive singleton reading is $1/4+O(j^2 2^{-j/2})$ for **every** specified product residue at precision at least $j+1$. Reduce modulo $2^{j+1}$ to apply the result. The same statement holds for $(hX)\mathbin{\&}(hY)$ at any odd multiplier by replacing $n$ with $h^2n$. [Earlier exact forced-answer classification and Bayes decision](cofactor-bit-completions.md) · [Modular hyperbola survey](https://arxiv.org/abs/1103.2879)

Two earlier list-size conclusions can also be made uniform in $n$. With one supplied AND the typical list has $2^{(1-H_2(1/4)+o(1))k}$ endpoints. With both supplied AND streams at multipliers $1,3$ or $1,5$, the terminal lists are tight and all prefix lists are polynomial-size with probability tending to one. These statements previously averaged over the product as well as the endpoints. They hold for **each** odd product residue when the endpoint is sampled uniformly within that residue's fiber. This still permits exceptional endpoints at any one product and does not give a worst-case guarantee for the actual prime factors of a fixed integer. [Single-stream theorem](shared-and-fiber-entropy.md) · [Two-stream theorem](two-fixed-and-streams.md)

## A uniform block comparison, including side information

Write $M=2^K$, $D=2^b$, $L=M/D$. Also permit a prescribed lower word $X=a\bmod 2^l$ with odd $a$ and $\max(1,l)+b\le K$; $l=0$ specifies only the original oddness. For any product $n$ and any compatible $a$, each upper word is individually uniform on $\mathbf Z/D$. Their joint distance to two independent upper words satisfies
$$
\mathrm{TV}\le C K^2\frac{2^lD}{\sqrt M}.\tag{1}
$$
The comparison measure also gives an efficient approximate sampler for any statistic depending solely on these two upper words.

For $1\le r<D$, the additive character of an upper word is $g_r(x)=e_D(r\lfloor x/L\rfloor)$, where $e_m(t)=\exp(2\pi i t/m)$. In its expansion on $\mathbf Z/M$, only frequencies $u=r\bmod D$ occur. Their coefficients are $L^{-1}\sum_{j=0}^{L-1}e_M(-uj)$ and have magnitude at most $\min(1,D/(2|u|_M))$, where $|u|_M$ is the least unsigned distance to zero. Summing the harmonic series in this one residue class gives an $O(K)$ bound on the Fourier $\ell^1$ norm.

The normalized dyadic Kloosterman sum satisfies
$$
\left|\mathbb E_{x\text{ odd}} e_M(ux+vn/x)\right|\le 8\sqrt{2\gcd(u,v,M)/M}.
$$
To verify it, cancel the common dyadic factor from $u,v,M$ and set the resulting modulus to $2^e$. For $e\ge2$ write $x=z+2^s j$ with $s=\lceil e/2\rceil$. Its inverse modulo $2^e$ is $z^{-1}-2^s jz^{-2}$. Summing over $j$ kills the term unless $uz^2=vn\bmod2^{e-s}$ with the reduced coefficients. If just one coefficient is even this is impossible. If both are odd there are at most four square roots modulo $2^{e-s}$, giving at most $4\cdot2^{2s-e}$ surviving $z\bmod2^s$. Division by the number of odd residues yields the bound. The two smallest moduli obey it trivially. This is the standard complete-sum estimate at a power of two.

Conditioning on $X=a\bmod 2^l$ inserts a normalized additive character filter with $\ell^1$ norm at most $2^l$. Its frequency shifts are multiples of $M/2^l$, hence multiples of $D$. Every frequency in the expansion of $g_r(X)$, after any such shift, therefore retains the exact valuation $v_2(r)$. For $r,s\ne0$ the joint upper-word character mean satisfies
$$
\left|\mathbb E[g_r(X)g_s(Y)\mid X=a\bmod2^l]\right|
\le C K^2 2^l\sqrt{\gcd(r,s,D)/M}.
$$
If exactly one index is zero, the mean vanishes by the exact marginal. By Cauchy–Schwarz and Parseval on the $D$-by-$D$ grid, its total variation from uniform is at most half the square root of the sum of squares of the nonconstant character means. The sum of $\gcd(r,s,D)$ over $1\le r,s<D$ is at most $\sum_{i=0}^{b-1}2^i(D/2^i)^2<2D^2$, proving (1). The Parseval step controls all bins jointly. Summing a separate discrepancy bound over each rectangle would lose a further factor of $D$. [Prior prime-modulus use of rectangle discrepancy](growing-level-intersection.md)

## Fixed-product two-stream collision bounds

The following applies to any of the supplied AND schedules whose carry matrix gives a compatible-quadruple bound $C_h\lambda^b$ with $\lambda<8$ uniformly over a $b$-bit interval. In particular it includes the full schedules for $h=3,5$ and the corresponding periodic schedules of densities $2/3,3/5$. [Finite carry certificates](sparse-fixed-multipliers.md)

Take an alternative first differing at bit $d\ge4$, so both coordinates first differ there. Use the next $b=\min(k-d,\lfloor d/4\rfloor)$ bits. Their high-word differences satisfy $\Delta v=c\Delta u\bmod2^b$, with odd differences and $c=-n/a^2\bmod2^b$, where $a$ is the actual low $b$-bit endpoint. For fixed $n$, the slope need not be uniform. Each slope has at most four preimages among the $2^{b-1}$ equally likely odd $a$.

Conditional on $a$, the true high words have distance at most $O((d+b)^2 2^{(3b-d)/2})$ from independent words by (1) with $K=d+b,l=b$. Consider the event that these high words admit **at least one** alternative satisfying the slope, the selected observations, and some common incoming pair of multiplier carries. Taking this union over all carries removes their dependence on the unused lower digits. Replacing the true high words by independent ones costs at most the stated total variation on this event indicator. One cannot transfer an unbounded count of alternatives at the same cost.

In the independent comparison model, bound the indicator by the number of matching alternatives. Each possible quadruple of high words with odd differences prescribes just one slope. Averaging over $a$ costs at most $4/2^{b-1}$ for this slope and $4^{-b}$ for the actual independent word pair. The carry bound, summed over at most $h^2$ incoming pairs, gives probability $O_h((\lambda/8)^b)$. The additional variation error is $O((b+1)^2 2^{-b/2})$ uniformly for $d\ge4b$. Thus the probability of a complete alternative first differing at $d$ decays exponentially in $b$.

For $t\le d\le k-t$, one has $b\ge\lfloor t/4\rfloor$, and any given $b$ occurs for at most five $d$. Summing the geometric tail bounds the probability of a middle alternative by $O_h(e^{-c_h t})$ uniformly in $n,k$. As before, each observation has at most $4^t$ endpoints without a middle alternative. Given $n$ and the observed ANDs the true endpoint is uniform on its fiber, so $\Pr(F_k>2\cdot4^t\mid n)\le O_h(e^{-c_h t})$. Choose $t$ constant for tightness, or a sufficiently large multiple of $\log k$ and union over all prefix depths for polynomial maximum width with high probability. Fixed omissions of the lowest few observations are absorbed into the constant lower cutoff on $d$.

## Fixed-product one-stream entropy

Take $b=\lfloor k^{1/3}\rfloor$ and use the existing single-stream block argument after a free initial $4b$ bits. At each later full $b$-bit block the true high words conditioned on $a=X\bmod2^b$ are within $O((b+1)^2 2^{-b/2})$ of independent words by (1), uniformly over its starting depth $d\ge4b$. The slope distribution over $a$ is at most four times uniform on odd residues. Hence the earlier high-probability bound on the **true-path** compatible child count for one block transfers with at most four times its exceptional probability plus the variation error. This applies to the indicator that the count is excessive. A union bound controls every full block, after which the earlier equal-child mass proof bounds the entire list for fixed $n$. Its initial and final free blocks cost at most $5b$ bits of mass.

For the matching lower bound, in the independent comparison model each full block's AND weight differs from $b/4$ by at most $b^{2/3}$ except with probability $\exp(-\Omega(b^{1/3}))$. The same variation and union bounds transfer this event to every fixed $n$. The overall AND weight then differs from $(k-1)/4$ by $O(k^{8/9})$. There are at most $2^{H_2(1/4)(k-1)+O(k^{8/9})}$ such AND words, and $2^{k-1}$ equally likely endpoints for this fixed product. Counting small fibers gives the lower tail exactly as in the original entropy proof. Uniformly in $n$, the probability that $\log_2 F_k$ differs from $(1-H_2(1/4))k$ by more than $Ck^{8/9}$ is at most $\exp(-c k^{1/9})$ for all sufficiently large $k$. No sampling assumption about primes conditioned on a specific full integer is used. [Original one-block probability and leaf mass argument](shared-and-fiber-entropy.md#an-exact-average-for-a-block)

## Checks

The checker verified 218 literal upper-word Fourier expansions and the dyadic complete-sum inequality at small moduli. At $k=20,b=4$ the largest unconditioned upper-word distance among the three products $1,-1,82791$ was $0.011902$. With a prescribed lower four-bit word $1$ or $3$, the largest distance was $0.048828$. At $k=16$ and the fixed product $1$, exhaustively grouping all 32,768 ordered endpoints by each sparse two-stream observation gave a pair-weighted median fiber of eight in both cases. The largest fibers had 128 and 130 endpoints. This difference between typical and exceptional endpoints is present even when the product is fixed. [Reproduce](check_fixed_product_and_blocks.py)
