# Local error tolerance for the two fixed binomial streams

For $h=3$ or $5$, take both fixed readings on each eligible bit level of a balanced semiprime. A report is a pair of bits. A mismatch in either component charges one error at that position. Let $k$ be the decoding depth and $t=\lceil C\log k\rceil$ for a sufficiently large constant $C$. An evaluator may produce wrong reports at up to a fraction $\delta$ of the positions in **each contiguous interval of length at least $t$** wholly inside the reported levels. The guarantee allows adversarial errors subject to this local budget, simultaneously over all such error patterns for each good input. It covers $\delta=1/90$ for $h=3$ and $\delta=1/60$ for $h=5$. The capped decoder takes polynomial time and recovers the factor on a density-one set in the dyadic prime model. Independent report errors with fixed probability strictly below the stated rate satisfy the local budget with high probability by the Chernoff bound and a union bound over at most $k^2$ windows, on increasing $C$ if needed. A total error budget with unlimited long bursts would give a different promise. [Exact case and prime-sample transfer](fixed-binomial-streams.md) · [Decoder](noisy_fixed_streams.py)

## Counting nearby blocks

Let $T$ be the four-carry transition matrix in the [exact proof](two-fixed-and-streams.md), and let $U$ admit all 16 digit quadruples. Weight an agreeing position by one and a disagreeing position by $1/32$. Write $M=T+(U-T)/32$ and $v=(32M)^{11}\mathbf1$. Integer arithmetic gives the following componentwise bounds $Mv\le\lambda v$. The last inequality is checked using rational numbers before taking roots.

| $h$ | Rational upper bound $\lambda$ | $\delta$ | Strict inequality |
|---:|---:|---:|---|
| 3 | $37/5$ | $1/90$ | $(37/5)^{90}32^2<8^{90}$ |
| 5 | $71/10$ | $1/60$ | $(71/10)^{60}32^2<8^{60}$ |

It follows that $M^b\mathbf1$ is componentwise at most $(\max v/\min v)\lambda^b$. All word quadruples disagreeing on at most $2\delta b$ positions have weight at least $32^{-2\delta b}$. Their count from any carry state is consequently $O_h((\lambda32^{2\delta})^b)$. [Exact integer certificates](check_noisy_fixed_certificate.py)

Fix a prefix of length $j\le k$ and let $d$ be the index of the first difference between the true first coordinate and a possible alternative. On the next $b=\min(d,j-d)$ bits, the two odd differences specify one odd slope modulo $2^b$. Its uniform probability is $2^{1-b}$, independent of the two actual high blocks, which have joint probability $4^{-b}$. A union over possible incoming multiplication carries accounts for their dependence on the slope. Hence the probability of **any** such alternative block within $2\delta b$ disagreements of the true block is $O_h(\rho^b)$, with $\rho=\lambda32^{2\delta}/8<1$. This event is defined without reference to a chosen noisy report.

For each $j$, summing over $t\le d\le j-t$ costs $O_h(\rho^t)$ since each $b$ occurs at most twice. Sum over $j\le k$ and transfer to primes at full depth with the earlier $O(k^2)$ density factor. Except on a set of probability $O_h(k^3\rho^t)$ there is no middle first difference of this kind, for any prefix. Choosing $C$ large enough gives the stated conclusion. Fixed bottom omissions all precede $t$ for sufficiently large $k$.

## Why a capped search tolerates adversarial errors

Try each low $t$-bit odd residue as a separate seed. Extend it by the publicly known product congruence, retaining children only when every completed interval respects the discrepancy budget against the supplied reports. Drop a seed entirely if its width exceeds $2^t$. The true seed survives. Indeed every admissible rival in that seed starts with the same $t$ bits. If it first differed in the middle, its length-$b$ block and the true block would both be within $\delta b$ of the same report and hence within $2\delta b$ of each other, excluded on the good-input event above. Thus only the top $t$ bits of the current prefix can vary in that seed. When $j<2t$ there are at most $2^t$ possible extensions directly. The cap is safe for the correct seed for every admissible report on the same good input. False seeds may die or be discarded at any depth. Final candidates are verified by integer division.

Both the number of seeds and the per-seed cap are $O(2^t)$, polynomial in $k$ for fixed $C$. The reference checker uses small explicit $t$ and covers cases with deliberate corruptions and false-seed abandonment. Its finite tests and the asymptotic guarantee use separate parameters. Constructing reports of the requisite quality from the unfactored input remains open. [Checks](check_noisy_fixed_streams.py) · [Earlier error models for one stream](scored-shared-bits.md)
