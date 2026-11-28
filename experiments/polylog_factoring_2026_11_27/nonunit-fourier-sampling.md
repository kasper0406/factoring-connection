# Fourier samples can reveal a factor while all interval masses look uniform

Let $N=pq$ for distinct odd primes $p<q<2p$. Take $f(x)=1$ for nonunits modulo $N$ (including zero) and zero elsewhere, and put $s=p+q-1=\sum_x f(x)$. An ideal Fourier sample of the normalized state $f/\sqrt s$ gives a factor with probability tending to one. Yet the exact largest discrepancy from uniform over **all contiguous frequency intervals** is only $(p+q-2)/N$. The example distinguishes samples, additive interval estimates, and conditional branch estimates. Shor uses the state obtained by coherent modular exponentiation. This calibration uses a different state whose preparation already has factoring power. [Shor](https://doi.org/10.1137/S0097539795293172) · [Earlier Fourier-side order-finding audit](../../notes/periodfind_results.md) · [Existing quarter-window predictor reduction](quarter-fourier.md)

The earlier ledger already gave the spectrum and its splitting probability. The interval error and prefix precision below quantify what a classical estimator would have to supply. [Earlier nonunit-mask calculation](https://github.com/kasper0406/factoring-connection/blob/6457469a8c353e2e24d45955b515457f6c52c4ee/experiments/polylog_factoring_2026_11_27/README.md#L140-L164)

## Exact distribution and interval discrepancy

Write $\widehat f(k)=\sum_x f(x)e^{-2\pi ikx/N}$ and $\mu(k)=|\widehat f(k)|^2/(Ns)$. Add the indicators of the multiples of $p$ and $q$, then subtract the extra copy of zero. The two subgroup sums give $\widehat f(k)=p\,1_{p\mid k}+q\,1_{q\mid k}-1$, so
$$
|\widehat f(k)|^2=1+p(p-2)1_{p\mid k}+q(q-2)1_{q\mid k}+2N1_{k=0}.
$$
An interval $I$ of length $L$ contains either $\lfloor L/d\rfloor$ or $\lceil L/d\rceil$ multiples of $d$ for $d=p,q,N$. Their counts therefore differ from $L/d$ by at most $1-1/d$. Subtracting $L/N$ from $\mu(I)$ gives a sum of precisely these three count discrepancies, weighted by their positive coefficients in the display divided by $Ns$. All three attain their positive maxima at $I=\{0\}$. Thus
$$
\max_I\left|\mu(I)-\frac{|I|}{N}\right|=\mu(0)-\frac1N=\frac{s-1}{N}.
$$
This also holds for cyclic intervals. Returning $|I|/N$ is an $N$-only oracle that meets any fixed inverse-polynomial tolerance in the **bit length** on all sufficiently large balanced inputs.

The outcomes coprime to $N$ have weight one, and zero has weight $s^2$. All other outcomes reveal a proper gcd. Hence
$$
\Pr_{k\sim\mu}(1<\gcd(k,N)<N)=\frac{\varphi(N)}{N}\left(1-\frac1s\right).
$$
For balanced $p\ge5$ all nonunits have mass at least $1/N$, while the units have mass $1/(Ns)$. Total variation from uniform is then also exactly the displayed split probability. Both quantities tend to one as $p\to\infty$. Ordinary basis measurement of the *input* state itself gives a factor with probability $1-1/s$. Filtering a uniform state by the gcd predicate accepts with probability only $s/N$. Thus access to a normalized nonunit state is already a substantive resource.

## A smaller absolute error suffices by bisection

Assume $N\ge1024$ and $1\le t\le\lfloor\sqrt N\rfloor$. This range includes neither $q$ nor $2p$. The first positive nonunit is $p$, giving
$$
\mu([1,t])=\frac{t+p(p-2)1_{t\ge p}}{Ns}.
$$
For $t<p$ this is at most $1/N$. For $t\ge p$ it is at least $1/(6\sqrt{2N})$, since $p-2\ge p/2$, $s<3p$, and $p>\sqrt{N/2}$. Consequently any prefix estimate with absolute error at most $1/(40\sqrt N)$ decides whether $p\le t$ by comparison with $1/(16\sqrt N)$. The gap holds even at $N=1024$. Binary search recovers $p$ in $O(\log N)$ queries using only $N$ and those estimates. Trial division handles smaller inputs. The tolerances have descriptions of $O(\log N)$ bits. An evaluator whose running time is polynomial in the *inverse* tolerance could still take power-size time. [Exact-rational decoder and checks](check_nonunit_fourier_sampling.py)

The easy uniform estimate and the accuracy sufficient to factor are both on the $N^{-1/2}$ scale, with separated constants. Conditional estimates give the probability of choosing the left half given that the sample lies inside a particular interval of positive mass. If each answer has absolute error at most $\eta$, clip to $[0,1]$ and draw one branch accordingly. Couple that draw to the true one using the same uniform variate. As long as the paths agree, they separate at the next step with probability at most $\eta$. A balanced binary tree of depth $\lceil\log_2N\rceil$ then gives a sample within $\eta\lceil\log_2N\rceil$ in total variation. This requires a guarantee at each visited positive-mass interval. Coarse absolute estimates on the original probability scale need not meet it after division by a small parent mass.

## Checks

| $N$ | Largest interval discrepancy | Total variation | Proper factor from a Fourier sample |
|---:|---:|---:|---:|
| $35$ | $0.285714$ | $0.623377$ | $0.623377$ |
| $143$ | $0.153846$ | $0.802676$ | $0.802676$ |
| $2{,}989{,}601{,}243$ | $3.7702\cdot10^{-5}$ | $0.99995343$ | $0.99995343$ |
| $494{,}727{,}028{,}433{,}071{,}069$ | $2.8592\cdot10^{-9}$ | $0.99999999643$ | $0.99999999643$ |

The checker computed literal complex Fourier amplitudes and independent integer prefix extrema on eight small semiprimes. For the four tabulated inputs it also checked 500 randomly chosen interval masses and subdivisions against the exact formulas. The 59-bit fixture uses the seeded primes 633,064,633 and 781,479,493. On the two larger inputs, estimates perturbed toward the wrong decision at every call still yielded the smaller factor after 16 and 30 prefix questions. The factor-supplied estimator is for calibration. No efficient classical fine-precision evaluator is established. [Reproduce](check_nonunit_fourier_sampling.py)
