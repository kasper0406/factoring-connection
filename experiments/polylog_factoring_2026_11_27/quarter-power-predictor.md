# Direct gcd for quarter predictors based on large power kernels

Consider a predictor that applies the public quarter rule to $c^e\bmod N$ instead of the hidden $cu\bmod N$. High exponents are cheap to evaluate and fall outside the low-degree feature results. In a regime where the local kernels of $c\mapsto c^{e-1}$ are large, any useful correlation is already detected by simple gcds on $c^{e-1}$. The exponent can be chosen in any public way. [Weak prediction reduction](quarter-fourier.md) · [Previous power sampling comparison](power-sampling.md)

Let $N=pq$ with odd $p<q$, let $u$ be a mixed root, and use the quarter sign $g_0$ from the Fourier reduction. For a uniform unit $c$, put $v=c^{e-1}\bmod N$. Conditional on $v$, the CRT components of $c$ are independent uniform multiplicative cosets of subgroups of sizes
$$d_p=\gcd(e-1,p-1),\qquad d_q=\gcd(e-1,q-1).$$
For a prime $r$ and a nonzero additive frequency, the average of the additive character on such a coset has magnitude at most $\sqrt r/d_r$. Expand the coset indicator into the multiplicative characters trivial on the subgroup and apply the Gauss-sum bound to each. Set $\eta=\max(\sqrt p/d_p,\sqrt q/d_q)$.

## Fourier resonance

On a fixed coset, the predictor is $g_0(vc)$ and the target is $g_0(uc)$. Truncate the Fourier expansion of the target to absolute frequencies $0<|k|\le K$. Its omitted part contributes $O(K^{-1/2}+1/N)$ to correlation over uniform units, by Parseval, the quarter-rule coefficient bound, and $N/\varphi(N)=O(1)$. Expand the predictor with its full Fourier series, whose absolute sum is $O(\log N)$. Each resulting term uses the additive frequency $hv-ku$ for some predictor harmonic $h$. If this is nonzero modulo $N$, at least one local coset bounds its conditional expectation by $\eta$. The total nonresonant contribution is therefore $O(\eta\log(2K)\log N)$.

For each target frequency $k$ and value $v$, the resonant predictor frequency is unique modulo $N$, namely $h=ku/v$. Its quarter coefficient is at most $1/|h|_N$. Resonances with $|h|_N>J$ contribute at most $O(\log(2K)/J)$. Write $\rho$ for the probability over $c$ that there is a resonance with $0<|h|\le J$ and $0<|k|\le K$. The remaining resonances contribute $O(\rho\log(2K))$. In total,
$$
\left|\mathbb E\,g_0(c^e)g_0(cu)\right|
=O\left(K^{-1/2}+N^{-1}+\eta\log(2K)\log N+
\frac{\log(2K)}{J}+\rho\log(2K)\right).
$$

## Recover a divisor directly

Assume $J,K<N/2$. Precheck the bounded nonzero integers by gcd. If any is a nonunit, this already gives a proper factor. Otherwise a small resonance $hv\equiv ku\pmod N$ gives $\gcd(hv-k,N)$ equal to the prime where $u\equiv1$. Try all bounded pairs on each sampled $v$. The procedure uses neither the predictor's output nor the hidden $u$.

For an inverse-polynomial correlation, choose $J,K$ polynomial in its inverse and the bit length so that the truncation terms are smaller than the signal. In the regime where the $\eta$ term is also smaller, the direct gcd has inverse-polynomial success. The subgroup-size hypothesis is part of the theorem. Small local kernels and more general high-degree maps remain outside this argument.

The diagnostic exhausts unit scales for identity, arbitrary and factor-planted exponents on three moduli. The planted cases use $e-1=\mathrm{lcm}((p-1)/2,(q-1)/2)$, with both primes $3\bmod4$. Half the scales then give $v=\pm u$ and the direct gcd succeeds on each. Run `python check_quarter_power_predictor.py`.
