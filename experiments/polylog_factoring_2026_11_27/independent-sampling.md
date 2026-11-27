# Pairwise gcd for independent CRT sampling

Let $N=pq$ for distinct odd primes. Draw a unit $Z$ from any distribution whose residues at $p$ and $q$ are independent. Write $s$ for its probability of splitting $N$ by the square-residue decoder. For two independent draws $Z,Z'$ from the sampler, put
$$C=\Pr(1<\gcd(Z-Z',N)<N).$$
Then either $C\ge s^2/4$, or the pair of most probable local residues is itself a successful decoder input and occurs with probability above $1/4$. This includes nonlinear nonuniform samplers whose local outputs are independent. The gcd needs two samples and no factor-supplied anchor.

## Weighted forest bound

Write $\alpha_i,\beta_j$ for the two marginal probabilities, $m_p=\max_i\alpha_i$ and $m_q=\max_j\beta_j$. The graph of successful residue pairs is a forest. Root each component and orient its edges from parent to child. For a left child $i$, its edge contributes at most $\alpha_i m_q$ to $s$. Each left vertex is a child at most once, so the sum over these edges is at most $m_q$. Right children contribute at most $m_p$. Hence $s\le m_p+m_q$. [Forest proof](power-sampling.md)

## Pairwise gcd

Put $A=\sum_i\alpha_i^2$ and $B=\sum_j\beta_j^2$. The two independent samples agree at the first prime with probability $A$ and at the second with probability $B$. The agreement events are independent, so $C=A+B-2AB$.

* If $A,B\le1/2$, then $C\ge(A+B)/2\ge(m_p^2+m_q^2)/2\ge s^2/4$.
* If $A\ge1/2\ge B$ or vice versa, then $C\ge1/4$.
* If $A,B>1/2$, then $m_p\ge A>1/2$ and $m_q\ge B>1/2$. If the joint mode is a successful edge, its probability exceeds $1/4$. Otherwise $s\le1-m_pm_q\le (1-m_p)+(1-m_q)\le(1-A)+(1-B)$. Applying the first case to $1-A,1-B$ gives $C\ge s/2$.

All quantities involving the primes are used only to prove the bound. A public trial can simply run the decoder on both samples and take the gcd of their difference.

## Approximate independence

If the sampler is within total variation distance $\delta$ of a product distribution on the unit pairs, its one-sample success rate differs from that of the product distribution by at most $\delta$. The two-sample laws differ by at most $2\delta$ by coupling each draw. Thus either the actual square-decoder probability exceeds $1/4-\delta$ by the mode case, or the actual pairwise gcd probability is at least
$$\frac{(\max(0,s-\delta))^2}{4}-2\delta.$$
A rare but inverse-polynomial success rate with a negligible independence error therefore already gives inverse-polynomial gcd success. The alternative is a frequent decoder hit on a single output. The proof gives no corresponding guarantee for samplers with substantial dependence between the local residues.

## Checks

The check script constructs the successful edges by CRT and direct decoding, then evaluates the two probabilities as exact fractions on 2,000 seeded marginal pairs over five semiprimes. Run `python check_independent_sampling.py`.
