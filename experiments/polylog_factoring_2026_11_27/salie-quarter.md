# An approximate Salié sum suffices for factoring

For an odd input with exactly two distinct prime divisors, a constant-precision approximation to one normalized Salié sum on uniform unit scales suffices for randomized polynomial-time factoring. Each summand is computable from public input by modular inversion and the Jacobi symbol. The sum has substantial cancellation. Single-term importance sampling needs order $N$ independent terms for constant mean-square precision on balanced semiprimes. [Fourier recovery](quarter-fourier.md)

## Exact identity

Write $e_N(t)=\exp(2\pi i t/N)$ and $\chi(t)=(t/N)$. For a unit $a$ set
$$S_N(a,b)=\sum_{m\in(\mathbb Z/N)^\times}\chi(m)e_N(am+bm^{-1}).$$
Let $G_N=\sum_{x\bmod N}e_N(x^2)$ be the quadratic Gauss sum, of value $\sqrt N$ for $N\equiv1\pmod4$ and $i\sqrt N$ otherwise. This definition also handles repeated factors. Then
$$\frac{S_N(a,b)}{G_N}=\chi(a)\sum_{v^2\equiv4ab\bmod N}e_N(v).$$

To derive the identity, take the additive Fourier transform in $b$ of the root sum. Its unnormalized coefficient at $z$ is $\sum_v e_N(v-zv^2/(4a))$. If $z$ is a nonunit, translation by $N/r$ for any prime $r\mid(z,N)$ preserves the quadratic term and multiplies the linear phase by a nontrivial $r$-th root of unity. The coefficient thus vanishes. For unit $z$, complete the square and evaluate the quadratic Gauss sum. This gives $G_N\chi(-z/(4a))e_N(a/z)$. Fourier inversion followed by $m=z^{-1}$ gives the identity, using $G_N^2=\chi(-1)N$.

Take $a=1$, $b=c^2/4\bmod N$, and let $u$ be a mixed root. For a unit $c$, the four roots are $\pm c$ and $\pm cu$. Hence
$$\frac{S_N(1,c^2/4)}{G_N}
=2\cos(2\pi c/N)+2\cos(2\pi cu/N).$$

## Recover a factor from approximate values

Suppose a public algorithm approximates the real normalized sum with root-mean-square error at most $1/2$ on uniform units, including its internal randomness. Subtract the known cosine, divide by two, and clip to $[-1,1]$. The resulting value has mean-square distance at most $1/16$ from $g(c)=\cos(2\pi cu/N)$. Convert it to a randomized sign of that expected value if needed by the sampler, and let $F(c)$ be the expected sign. By Jensen the same error bound holds for $F$.

For an odd modulus with two distinct prime divisors, the unit average of $g^2$ is at least $1/2$, by the Ramanujan sum at the unit frequency $2u$. Cauchy–Schwarz gives $\mathbb E_{c\in U} F(c)g(c)\ge1/4$. Extend $F$ by zero on nonunits. A uniform gcd probe can supply a factor directly. Otherwise the correlation on all residues is at least $1/8$, since $\varphi(N)/N\ge1/2$. This is the real part of the Fourier coefficient at $u$. The sampled heavy-frequency search recovers $u$ or $-u$ among its list. Test each output frequency minus 1 by gcd. This avoids a harmonic search over an unknown multiple of $u$.

## Variance of sampling individual terms

The normalized complex summands each have absolute value $N^{-1/2}$. For any positive sampling probabilities $q_m$, the unbiased one-term importance estimator has second moment $\sum_m |a_m|^2/q_m\ge(\sum_m |a_m|)^2=\varphi(N)^2/N$ by Cauchy–Schwarz. On the two-prime-divisor family the target magnitude is at most 4, so the variance is $\Omega(N)$ for large inputs.

Even estimating just the real part has average variance $\Omega(N)$ on balanced semiprimes. Its absolute term magnitudes are $N^{-1/2}|\cos\theta_{m,c}|$, where the phase includes division by $G_N$. Since $|\cos\theta|\ge\cos^2\theta$, it suffices to average the latter. For each fixed $m$, the nonconstant part in $\cos^2\theta_{m,c}$ is a unit quadratic phase in $c$. The complete Gauss bound and deletion of nonunits give an average $O(N^{-1/2})$ over unit $c$ on balanced semiprimes. Thus the average sum of absolute real magnitudes is $\Omega(\sqrt N)$. Jensen and the same importance-estimator bound give the stated average variance. Averaging $L$ independent one-term samples needs $L=\Omega(N)$ for constant mean-square error within this model. A method that aggregates oscillatory terms algebraically is outside this variance calculation.

## Checks

The diagnostic compared the public normalization with the factor-supplied cosine identity on 138 unit scales across five moduli, including two repeated-factor inputs. It also computed the minimum possible variance of a real one-term estimator for each tested scale. Run `python check_salie_quarter.py`.
