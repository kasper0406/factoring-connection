# Sampling variance after grouping Salié phases

Pointwise access to the signed phase coefficients still leaves a variance obstruction to estimating the signed half-interval sum by independent importance sampling. This holds even with optimal sampling probabilities, on balanced squarefree semiprimes. [Pointwise coefficient oracle](salie-phase-weights.md) · [Signed interval target](salie-phase-bins.md)

At a prime $r$ the local signed coefficient is $L_r(t)=\chi_r(t+c)+\chi_r(t-c)$ for unit $c$. Its squared sum over $t$ is $2r-4$. The two square terms contribute $r-1$ each, and the cross term uses $\sum_t\chi_r(t^2-c^2)=-1$. Every nonzero additive Fourier coefficient of $L_r^2$, under the unnormalized convention, is $O(\sqrt r)$. The square terms differ from constants at single points, and the cross term satisfies the standard mixed quadratic character-sum bound.

By CRT, $W_c^2$ has mean $(2-4/p)(2-4/q)$. Its unnormalized Fourier coefficients off zero are $O(N^{3/4})$ on balanced inputs. If both local frequencies are nonzero use $O(\sqrt p\sqrt q)$. If only one is nonzero use its square-root bound and the other prime's linear-size bound. The absolute Fourier coefficient sum of an interval indicator is $O(\log N)$ in the normalized convention. Hence on the fixed half interval $I$,
$$\sum_{t\in I} W_c(t)^2=(2-4/p)(2-4/q)|I|+O(N^{3/4}\log N)=\Omega(N),$$
uniformly for unit $c$. Since $|W_c(t)|\le4$, the sum of absolute magnitudes over the interval is also $\Omega(N)$. The signed sum is $O(\sqrt N\log N)$ by the earlier interval Fourier bound.

For the normalized target $\sum_{t\in I} W_c(t)/\sqrt N$, the minimum second moment of an unbiased one-term importance estimator is the square of the sum of absolute normalized magnitudes. It is therefore $\Omega(N)$. Subtracting the squared mean only removes $O(\log^2 N)$. Independent repetitions need $\Omega(N)$ samples for constant mean-square error in this model. This applies after grouping by phase, so separate exact cancellation within a phase bucket does not improve the sampling exponent. Structured aggregation across distinct phases lies outside the proof.
