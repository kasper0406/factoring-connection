# Can coarse square-root phases predict a deeper block of binomial readings?

For distinct primes sampled uniformly in $[R,2R)$ with $R=2^k$, a square-root phase of $N=pq$ gives a weak forecast at an upper dyadic level. Arbitrarily combining polylogarithmically many coarse measurements of these phases, even with the entire low word $N\bmod R$, gives no constant improvement over the all-zero forecast on a long middle block. In particular this source cannot generate the two accurate report streams needed by the decoder. The restriction concerns which features of $N$ enter the forecast. The full integer can be used for final verification. [Upper-level calculation](upper-bit-square-root-phase.md) · [Uniform fixed-product block theorem](fixed-product-and-blocks.md) · [Required reporting accuracy](noisy-fixed-streams.md)

## How many joint observations can the phases give?

Consider the phase $2\sqrt{uvN}/D+\phi$ modulo one, with positive integers $u,v\le A$. This includes applying the original square-root phase to publicly scaled inputs $uvN$. Use dyadic periods $D_{\min}\le D\le 4R$. At each triple let the observations collectively change at at most $B$ boundary points per unit period. For example a sign test has at most two, and a quantizer with $s$ uniformly spaced bins has $s$. Multiple fixed offsets and harmonics are allowed by charging their combined boundaries. For each given low product word the entire menu, including its offsets, can be different.

The full menu has at most

$$ C=O\big(1+B A^2(1+\log R)+B A^3 R/D_{\min}\big) $$

joint cells on $R\le z=\sqrt N<2R$. Indeed a single phase passes through at most $O(B(1+\sqrt{uv}\,R/D))$ boundaries. Sum over $u,v,D$, using $\sum_{u,v\le A}\sqrt{uv}\le A^3$ and $\sum_D 1/D\le 2/D_{\min}$. Allowing boundary points to give special output only doubles the partition size. Adaptive queries drawn from this menu and arbitrary postprocessing still give at most $C$ final predictions for a fixed low word. Thresholds chosen adaptively from outside the enumerated menu would require their own count. At high enough precision even the unscaled square root determines $N$.

## Prime-pair bound on middle-block forecasts

Fix an odd multiplier $h$, a starting bit $j$ and a length $b$ with $j+b\le k$. The target word consists of bits $j,\ldots,j+b-1$ of $(hp)\mathbin{\&}(hq)$. For the binomial interpretation use $h=1$ or $3$ and the eligible dyadic levels of the decoder. The forecasts may also use $N\bmod 2^{k+r}$, $r\ge0$. For $0\le\delta<1/4$ define $T_b(\delta)=\Pr[\operatorname{Bin}(b,1/4)\le\lfloor\delta b\rfloor]$. Then

$$
\Pr_{p,q}\big[\text{block Hamming error}\le\delta b\big]
\ll k^2\left(2^r C T_b(\delta)+(j+b)^2 2^{(b-j)/2}\right). \tag{1}
$$

Start with independent uniform odd integers $X,Y$ in $[R,2R)$, and condition on any odd $n=XY\bmod R$. The residue of $X$ is uniform and $Y=n/X\bmod R$. Reduced to precision $2^{j+b}$, the product for $hX,hY$ is $h^2n$ and the first residue is again uniform. The fixed-product block theorem gives total variation $O((j+b)^2 2^{(b-j)/2})$ between the two upper $b$-bit words at that precision and two independent uniform words. Taking their bitwise AND therefore gives independent Bernoulli($1/4$) bits up to the same error. [Block theorem, equation (1)](fixed-product-and-blocks.md#a-uniform-block-comparison-including-side-information)

For each fixed $n$, there are at most $2^r$ refinements of the low word and $C$ possible full-menu cells per refinement. Consequently the predicted word belongs to a set of at most $2^r C$ fixed words, regardless of any statistical dependence between the outputs and the factors. For iid Bernoulli($1/4$) bits, the radius-$\lfloor\delta b\rfloor$ Hamming ball of greatest mass is centered at zero. Every predicted one replaces one mismatch probability $1/4$ by $3/4$. The mass of the union of all balls is thus at most $2^r C T_b(\delta)$. Apply the variation bound once to the union, average over $n$, and condition on distinct primes, an event of probability $\asymp k^{-2}$ by the prime number theorem. This proves (1), also for randomized forecasts by fixing their independent random coins first.

For instance take $j=\lfloor k/2\rfloor$, $b=\lfloor k/4\rfloor$. All these bits are eligible for $h=1,3$ for large $k$. If $A,B,R/D_{\min}$ grow at most polynomially in $k$ and $r=O(\log k)$, then $2^r C$ is polynomial in $k$. For every fixed $\delta<1/4$, the binomial tail and the right side of (1) decay exponentially in $k$. Consequently even the *mean* fraction of correct middle-block bits cannot exceed $3/4+o(1)$ using these inputs. A constant-zero prediction achieves $3/4+o(1)$ in mean on the prime sample, by applying both binomial tails in the same comparison and conditioning argument. Arbitrary further modular arithmetic on $N$, higher-precision phase digits and genuinely new input transformations lie outside this bound.

## How many additional observations would accurate reports require?

The same counting argument yields a stronger requirement when the number $L=2^r C$ of possible predicted middle words can be exponential. Fix $0\le\delta<1/4$. To predict the $b$-bit block with at most $\delta b$ errors with probability bounded away from zero, the permitted additional observations must allow at least
$$
\log_2 L\ \ge\ b\,[H_2(1/4)-H_2(\delta)-o(1)]
$$
in the example $j=\lfloor k/2\rfloor$, $b=\lfloor k/4\rfloor$. Here $H_2$ denotes binary entropy. At 1% error this is at least $0.7304b$ extra bits of possible output. This counts distinguishable words. It does not assert that an arbitrary collection of that size contains the right prediction.

Indeed for an iid Bernoulli($1/4$) word, restrict to weights within $b^{2/3}$ of $b/4$. The omitted mass is $\exp(-\Omega(b^{1/3}))$, and every retained word has probability at most $2^{-bH_2(1/4)+O(b^{2/3})}$. A ball of radius $\lfloor\delta b\rfloor$ has at most $2^{bH_2(\delta)}$ elements, by the elementary binomial bound. Apply the total variation error once to the union and then charge the $O(k^2)$ cost of conditioning on true prime pairs. If the claimed inequality is missed by a fixed positive fraction of $b$ the success probability tends to zero. The weaker full-ball probability in (1) has the advantage of giving exponential tails near 75% accuracy at no typical-weight cutoff. The same argument works for any public features with a bounded number of possible outputs at each fixed low product.

## Exact finite check

The checker constructs all cut points as rationals in the variable $N$ for $u,v\le3$, $D\in\{R/2,R,2R\}$ and four phase bins. It verifies the cuts against direct integer square roots and finds 154 real cells at every tested dyadic scale. It then fits the Hamming-optimal word separately in every cell for two fixed odd products at $k=12,15,18$, enumerating **all odd endpoints** at each product. At $k=18$, for bits 9 through 14, the oracle-fitted mean errors were $0.230258$ at product residue 1 and $0.251500$ at residue 82791. The corresponding true-word variation distances from iid Bernoulli($1/4$) were $0.079590$ and $0.048462$. The bound for a union of predicted Hamming balls was also checked with exact probabilities. [Reproduce](check_quantized_tangent_transport.py)
