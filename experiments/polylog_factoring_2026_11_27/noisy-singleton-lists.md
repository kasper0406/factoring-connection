# What partial accuracy on the fixed singleton bits buys

For a balanced odd semiprime $N=pq$, the fixed dyadic singleton tests read the bits of $p\mathbin{\&}q$. Suppose a proposed evaluator supplies the first $w$ **nonunit** bits, allows erasures, and can make mistakes. The decoder below builds one list that accounts for both the unknown orientations of the two factors and the errors in the supplied answers. It needs no error locations. On typical balanced prime pairs with 90% accuracy on these bits and no erasures, the deterministic work **after receiving the readings** is at most $N^{0.182334+o(1)}$. With 95% accuracy it is at most $N^{0.164355+o(1)}$. The exact-answer exponent was $1/8$. These are conditional guarantees about computing a factor from readings. [Exact-answer list completion](hh-residue-lists.md) · [Factorial ratios giving the bits](cofactor-bit-binomials.md) · [Implementation](noisy_singleton_lists.py)

## One tree for the two budgets

Let $m=2^{w+1}$ and let $s$ positions be marked as erasures. For a candidate ordered pair of odd endpoints $a,b$ with $ab=N\bmod m$, define $B(a)$ as the count of **observed** positions where $a,b$ differ. Let $E(a)$ be the Hamming error count between $a\mathbin{\&}b$ and the nonerased readings. For specified upper bounds $B,E$, the exact list of endpoints with $B(a)\le B$ and $E(a)\le E$ has size at most
$$
F(w,A,E)=2^A V(w-A,E),\qquad
A=\min(w,s+B),\qquad V(l,E)=\sum_{i=0}^{\min(l,E)}\binom li.
$$
The list can be built in $F\,\mathrm{poly}(\log N)$ bit operations for $w+1\le\log_2 N$. If every reported one is certified, enforce these $t$ positions exactly and replace $w$ by $w-t$ in both $F$ and $A$. The remaining $E$ possible mistakes are missing true ones.

To prove the bound, start from $a=b=1\bmod2$. At the next position the product and the current endpoints already determine whether the new two bits must agree. Changing the lift of $a$ flips both new bits. At an erased position there are two free choices. At an observed position with differing bits there are two orientations, both giving AND zero. A reported one there incurs a compulsory error, or kills the state when ones are certified. At an observed position with equal bits the two choices are $00,11$. One agrees with the reading, one consumes an error. A certified one at such a position forces $11$.

On a surviving path there are at most $A$ free choices. The decision about which kind of node comes next may depend on the entire prefix. With $l$ positions left and $0\le a\le l$ free choices remaining, $F(l,a,e)$ still bounds the number of leaves. At a free node the two subtrees together fit $2F(l-1,a-1,e)=F(l,a,e)$ when $a>0$. A free node with zero remaining capacity discards its path. At an error-scored node with $a<l$, the two child bounds sum to $F(l,a,e)$ by Pascal's rule, with $V(l,-1)=0$. With $a=l$ use the trivial bound $2^l$. Forced errors and pruned paths only reduce these counts. For the actual implementation the cap counts differing *observed* positions and every erased position can use one free choice. The same reasoning applies to each prefix, whose length and erasures do not exceed their final values, so its width also never exceeds $F(w,A,E)$. Certified-one positions either discard the path or give one compulsory child and can be omitted in this count.

Enumerating all subsets of up to $E$ wrong readings first would charge $V(w-s,E)$ even though opposite bits already force a zero on that candidate path. The joint tree pays $V(w-A,E)$ instead. Its bound is never greater than $2^w$. Given a candidate list of size $S$, the deterministic list completion uses
$$
\left(S+\min\{N^{1/5}(S/m)^{3/5},\ N^{1/4}S/m\}\right)\mathrm{poly}(\log N)
$$
bit work. If neither cap is known, try all $O(w^2)$ pairs in increasing order of $F$. The true caps suffice. Every earlier trial has at most the corresponding bound on candidates, and every returned divisor can be checked directly. [Shared list completion and its order setup](hh-residue-lists.md)

## Accuracy, erasures and the factoring exponent

Let the two distinct primes be sampled in $[R,2R]$. For a predetermined erasure pattern with $s=\sigma w+o(w)$ and $w-s$ growing proportionally to $w$, the true observed orientation count is $(w-s)/2+O(w^{2/3})$ on a density-one set. On independent uniform odd endpoints modulo $m$, this count is binomial with parameters $w-s,1/2$. The tail is exponentially small in $w^{1/3}$. For $m\le R$, comparing residue counts in the integer square to the number of prime pairs costs $O(\log^2 R)$. The argument applies for $w$ proportional to $\log R$, including the boundedly many possible precisions at $N\in[R^2,4R^2]$. Erasures selected independently of the primes give the same estimate by conditioning on the pattern. For input-dependent erasures, the deterministic tree bound still uses the actual observed orientation count.

Suppose at most a fraction $\rho$ of nonerased readings are wrong. For fixed $0\le\rho<1/4$ the bound is $S\le 2^{(\alpha+o(1))w}$ with
$$
\alpha=\frac{1+\sigma}{2}+\frac{1-\sigma}{2}H_2(2\rho).
$$
At $\rho\ge1/4$ use $\alpha=1$ for this guarantee. Take $w=\lfloor\log_2 N/(3+2\alpha)\rfloor$. Balancing candidate generation with the shared collision term gives exponent $\alpha/(3+2\alpha)$. This is strictly below $1/5$ whenever $\sigma<1$ and $\rho<1/4$. It counts individual inputs in the stated density-one set, without asserting a mean runtime bound. The error promise may be worst-case over the positions. An expected success rate by itself supplies no such per-input guarantee.

For an example with one-sided readings, retain each true nonunit one independently with probability $r$, report every other bit as zero, and use no erasures. With high probability, $t=rw/4+o(w)$ certified ones are reported, $E=(1-r)w/4+o(w)$ ones are missed, and the orientation count is $w/2+o(w)$. For fixed $r\in[0,1]$, the resulting list exponent is
$$
\alpha=\tfrac12+(\tfrac12-r/4)H_2\!\left(\frac{1-r}{2-r}\right).
$$
For $r=1/2$ this gives $\alpha\approx0.844361$, improving the $0.875$ obtained by simply fixing the reported ones. The sampler here assumes access to true bits and illustrates the guarantee.

## Why an all-zero predictor earns no exponential list saving

In the uniform odd-residue model each nonunit AND bit is one with probability $1/4$, so a predictor always reporting zero has 75% expected accuracy. A counting statement explains what happens even when the true orientation and error counts are supplied for the particular instance. Fix any odd product residue and any readings determined by that residue. Its $2^w$ ordered endpoints split into at most $(w+1)^2$ classes by their exact pair $(B(a),E(a))$. For a uniform endpoint the probability that its **exact-count** class has size below $\delta 2^w/(w+1)^2$ is at most $\delta$. The list from upper bounds contains that class. Taking $\delta=2^{-\sqrt w}$ shows that even the per-instance true caps leave $2^{w-o(w)}$ candidates with high probability. The prime comparison above transfers this conclusion at the logarithmic precisions used here. It applies to the low-product candidate list. Algorithms using other bits or arithmetic properties of the full integer $N$ require separate analysis. The verified-one variant here requires the reports of one to be valid, as they automatically are for the all-zero baseline.

The tests exhaustively compared 41,790 capped trees with direct enumeration of the entire product fiber. All 314 end-to-end runs recovered a factor from deliberately degraded native factorial-ratio readings, with 314 erasures and 248 flipped answers across the runs. On a 20-position planted transcript with two erasures and two errors, the decoder retained 2,932 residues. Its joint bound was 57,344, versus 176,128 for choosing the error subset independently. The all-zero baseline on 196,608 endpoints at 48 fixed product residues had 74.69% bit accuracy. Even with each endpoint's true caps, the median candidate list held 972 of the 4,096 endpoints. [Reproduce](check_noisy_singleton_lists.py)
