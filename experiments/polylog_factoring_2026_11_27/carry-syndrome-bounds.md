# Accidental factors from a binary carry template

A wrong carry template can still factor a semiprime: its corrected binary polynomial need only be divisible by the polynomial encoding the smaller factor. For distinct primes sampled uniformly from the same $b$-bit interval, the probability that **any mask in a predeclared menu of $2^{o(b)}$ masks** meets this weaker condition is at most $N^{-0.08592775+o(1)}$. The zero template with any fixed Hamming error radius is included. On 1,024 fresh products at each size, radius two split 114 inputs at 12 bits per prime, 15 at 16 bits, two at 20 bits, and none at 24, 32, 48 or 64 bits. Most finite successes used incorrect carry masks. [Public decoder](carry-parity.md) · [Earlier bound for the true carry mask](carry-prediction-bounds.md) · [Executable checks](check_carry_syndromes.py)

Write $U_z(T)=\sum_i z_iT^i\in\mathbb F_2[T]$ for the binary digits of a nonnegative integer $z$. Its inverse encoding reads the polynomial coefficients as bits. For $p<q<2p$, the smaller-factor decoder with mask $E$ factors $U_N+U_E$ and tests its binary-polynomial divisors as candidate integers in the balanced interval below $\sqrt N$. Without a divisor cap, it returns $p$ exactly when
$$
 s_p(N):=U_N\bmod U_p=U_E\bmod U_p.
$$
This polynomial remainder is the syndrome. The implementation's divisor cap can only remove successes. Since binary encoding is unique, a divisor read as the integer $p$ must be $U_p$.

## A bound uniform over the target syndromes

Let $R=2^{b-1}$, $0\le k\le b-2$, and divide $[R^2,4R^2)$ into pieces of the aligned blocks of width $2^{b+k}$. In each block allow **up to $M$ target syndromes**, which may depend arbitrarily on $p$ and the block. When $p,q$ are distinct random primes in $[R,2R)$, the probability that $s_p(pq)$ belongs to the allowed set is
$$
 \ll b^2\left(M2^{-k}+2^k(3/4)^{b-k}\right).
$$
The constant is absolute. The same estimate up to a constant holds when $p$ denotes the smaller prime after sorting the pair. In particular, masks may depend on the product's high prefix $\lfloor N/2^{b+k}\rfloor$. A rule allowed to see all of $N$ and choose from a predeclared menu is covered by taking all its possible syndromes. An unrestricted rule that may generate a fresh target from the full input could use the factorization to supply the true mask.

**Short polynomial multiples.** Fix an odd $b$-bit $p$, initially allowing every integer $q\in[R,2R)$. For $q<q'$ with their products in the same block and the same syndrome, put $d=q'-q$ and $W=(pq)\mathbin{\mathrm{XOR}}(pq')$. The high product bits cancel. Therefore
$$
 U_W=U_p V,\qquad 0\ne V\in\mathbb F_2[T],\quad \deg V\le k,\qquad
 0<d<2^{k+1}.
$$
The degree follows from $\deg U_p=b-1$ and $\deg U_W<b+k$. The bound on $d$ follows from $pd<2^{b+k}$.

**Counting ordinary integer gaps.** Fix $p,V$, which fixes $W$. The exact integer identity $x\mathbin{\mathrm{XOR}}W=x+W-2(x\mathbin{\mathrm{AND}}W)$, applied at $x=pq$, gives
$$
 2((pq)\mathbin{\mathrm{AND}}W)=W-pd.
$$
Let $a$ be the number of ones in the lowest $k$ bits of $W$. Modulo $2^{k+1}$, the doubled AND on the left has at most $2^a$ possible values. Since $p$ is odd, at most $2^a$ of the allowed gaps $d$ occur. For each fixed $d$, the same equation fixes every bit of $pq$ at which $W$ has a one. If $w$ counts the ones of $W$ in its lowest $b-1$ bits, at most $R2^{-w}$ values of $q$ are possible. Indeed $pq\bmod R$ runs through all residues exactly once on the full partner interval. For this $V$ the number of collision pairs in all blocks combined is at most $R2^{a-w}$.

**Average over the first factor.** Now take $p$ uniformly among odd $b$-bit integers. For fixed $V\ne0$, write $V=T^vV_0$ with $v\le k$ and $V_0(0)=1$. Multiplication by $V_0$ permutes the odd polynomials modulo $T^{b-1-v}$. The input $U_p$ is uniform there, as the fixed top bit has already been discarded. In particular the coefficients of $U_W$ at indices $k+1,\ldots,b-2$ are independent fair bits. All of them are included in $w-a$. Hence $\mathbb E_p[2^{a-w}]\le (3/4)^{b-k-2}$. Summing over fewer than $2^{k+1}$ nonzero $V$ bounds the expected number of within-block equal-syndrome pairs by $O(R2^k(3/4)^{b-k})$.

**Count members of the selected sets.** There are $O(R2^{-k})$ possible blocks. Within any block, each selected syndrome accounts for one occupied place plus at most the number of unordered collisions in its class, since $h\le1+\binom h2$. Apply this simultaneously to at most $M$ classes. Their collision terms sum to at most the total for that block, irrespective of the choice of targets. Summing the blocks and averaging over $p$ gives $O(MR2^{-k}+R2^k(3/4)^{b-k})$ successful partners per first factor. Divide by $R$, then condition on $p,q$ both prime and distinct, which costs $O(b^2)$ by the prime number theorem. This proves the bound without requiring equidistribution in a growing arithmetic progression.

With $A=\log_2(4/3)$ and $M=2^{o(b)}$, choose $k=\lfloor Ab/(2+A)\rfloor$. Both terms are $2^{-Ab/(2+A)+o(b)}$. Since $\log_2 N=2b+O(1)$, the exponent in base $N$ is $A/(4+2A)=0.0859277546\ldots$. More generally $M=2^{\rho b+o(b)}$ gives exponential decay for $\rho<A/(1+A)=0.2933049474\ldots$ by taking $k$ near $(A+\rho)b/(2+A)$. For instance, even an exhaustive Hamming ball of radius $\theta\,\operatorname{bitlength}(N)$ about a fixed word has exponentially small success when $2H_2(\theta)<A/(1+A)$ and $0\le\theta<1/2$, which includes every $\theta<0.02093389$. These are upper bounds. They do not identify the sharp exponent. A separate calculation treats a missing word at an interior position. [Interior-word extension](carry-syndrome-windows.md)

## What the finite experiment measured

The offline screen receives the smaller prime $p$. It computes $s_p(N)$ and the remainders of $T^j\bmod U_p$ at the permitted mask positions $1\le j<\operatorname{bitlength}(N)$. Radius zero succeeds exactly for zero syndrome. Radius one adds the individual monomial remainders. Radius two also admits their pairwise XORs at distinct positions. The implementation keeps a set of earlier remainders, so the pair test needs only a linear number of set lookups. The separate true-mask count compares $N$ with the coefficient word of $U_pU_q$.

| Bits per prime | Inputs | Uncapped split counts at radii 0, 1, 2 | True-mask counts at radii 0, 1, 2 |
| ---: | ---: | --- | --- |
| 12 | 1,024 | 2, 17, 114 | 0, 2, 4 |
| 16 | 1,024 | 0, 2, 15 | 0, 0, 0 |
| 20 | 1,024 | 0, 0, 2 | 0, 0, 0 |
| 24 | 1,024 | 0, 0, 0 | 0, 0, 0 |
| 32 | 1,024 | 0, 0, 0 | 0, 0, 0 |
| 48 | 1,024 | 0, 0, 0 | 0, 0, 0 |
| 64 | 1,024 | 0, 0, 0 | 0, 0, 0 |

Every one of the 131 positive inputs was replayed at all three radii with the public decoder, which received $N$ alone and used its default divisor cap. All results matched the screen. For example, $1{,}870{,}182{,}893=35{,}449\cdot52{,}757$ split at radius one even though its true mask had weight 17. A disjoint unfiltered public check used four inputs per size through 24 bits and one per size at 32, 48 and 64 bits. None of these 19 split. The decoder skipped no masks in any of the 450 public calls. The 7,168-input table counts an exact factor-supplied screen with public verification of all its positives.

The test also independently enumerates masks on the 78 smaller prime pairs in the original experiment. Across 368 fixed-first-factor and block-size cases it verifies the short-multiple and AND identities, bounds on the possible gaps, the exact average of the pinned-bit weights and the target-set bound. It encounters 9,425 collision pairs in those checks. The script's default runs the finite proofs and the seven-cohort screen. Passing the option $\texttt{--replay}$ also runs the public decoder. [Reproduce the experiment](check_carry_syndromes.py)

## Related results

Devevey, Sakzad, Stehlé and Steinfeld encode integer arithmetic carries by polynomials in their analysis of integer polynomial learning with errors. Lemire and Kaser use binary polynomial products and remainders in universal hashing, where collision probabilities are averaged over a separately chosen key. Here the same $p$ determines the integer progression $pq$ and its binary polynomial modulus $U_p$. [Carry-polynomial arithmetic](https://research.monash.edu/en/publications/on-the-integer-polynomial-learning-with-errors-problem/) · [Carryless multiplication and hashing, §4](https://arxiv.org/abs/1503.03465)

Lipmaa, Wallén and Dumas analyze differences of XOR under modular integer addition. Schmid studies the joint distributions of binary digit sums of integer multiples. The estimate above combines an exact XOR difference at gap $pd$ with a bounded-degree polynomial multiple of $U_p$, and counts separate syndrome classes within an aligned block. [Additive differential probability of XOR](https://mathexp.eu/dumas/publications/Lipmaa_Wallen_Dumas_2004_ADP.pdf) · [Joint binary digit sums of multiples](https://matwbn.icm.edu.pl/ksiazki/aa/aa43/aa4347.pdf)
