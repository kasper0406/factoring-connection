# Affine interval filters instead of bit masks

The general square-modulus reduction can replace its binary inner product filter by one additional modular interval. The count now uses two lower-input interval conditions and one upper-phase interval, with no Jacobi weights and no bit mask. It still works for repeated prime factors. [Input digit isolation](square-input-isolation.md)

## Exact covariance of the new signs

Let $s(t)=1$ for residues $0\le t\le(N-1)/2$ and $-1$ otherwise. Choose independent uniform $a,b\bmod N$ and use the sign $s(ax_0+b)$ for the lower input digit $x_0$. Hold these choices fixed over the entire upper-digit query trial. A nonunit $a$ encountered during construction can be tested by gcd.

For two lower digits $x,y$, let $g=\gcd(x-y,N)$ and $m=N/g$. Averaging over $a,b$ makes the two affine values uniform over pairs in a common residue class modulo $g$. Conditional on that class, the two values are independent and uniform in it. Since $m$ is odd, each class contains either $(m-1)/2$ or $(m+1)/2$ points of the positive half. Its mean sign is therefore $\pm1/m$. This gives the exact covariance identity, including coincident inputs,
$$\mathbb E_{a,b}s(ax+b)s(ay+b)=\left(\frac{\gcd(x-y,N)}N\right)^2.$$
Here covariance means the uncentered product expectation used in the variance expansion.

## Bound the variance on all odd inputs

Write $w_x=\gcd(x^2-c_0^2,N)/N$. The normalized reference subtraction changes a block by at most $w_x$. As the independent isolation window translates, the probability that it retains both $x,y$ is at most $w/N$, where $w$ is the window width. The sign identity therefore bounds the expected squared norm by $w/N$ times
$$\sum_{x,y\in U_N}w_x w_y\left(\frac{\gcd(x-y,N)}N\right)^2.$$
This sum factors over prime powers. For the local factor at $p^e$, use
$$\gcd(x-y,p^e)^2=\sum_{p^k\mid(x-y,p^e)}J_2(p^k),$$
where $J_2(1)=1$ and $J_2(p^k)=p^{2k}-p^{2k-2}$ otherwise. At depth $k$, group the local weights by residues modulo $p^k$ and let $S_k$ be the squared sum of their class masses. The two classes of $\pm c_0$ at depth $k\ge1$ have mass $1+(e-k)(1-1/p)$ each. All other classes have the same weights as they would at exponent $k$, with total squared mass at most $2/p$ by the earlier squared-weight formula. At depth zero the single mass is at most $2e+1$.

The contribution of $k=e$ is at most $2+2/p$. The rest is bounded by
$$4\sum_{j\ge1}(j+1)^2/p^{2j}\le24/p^2.$$
Here $j=e-k$ and $p\ge3$. Thus the local factor is at most $2(1+5/p)$, uniformly in $e$. For $Q=2^{\omega(N)}$, the global factor is at most $Q(1+\omega(N))^5$, since $1+5/p\le(1+1/p)^5$ and the latter product telescopes after ordering the primes.

At the isolation width, $wQ/N\le1$. The expected squared norm is therefore polynomial in the input bit length. The same Markov and clipping argument, with a larger polynomial threshold, preserves inverse-polynomial probability of a singleton trial with a usable hidden frequency. The normalized RMS error allowance remains $1/8$. The Fourier decoder and verified gcd tests then split the input in randomized polynomial time given the count oracle.

## Count formulation

Each sign can be queried as a separate unweighted count. Apart from the interval filters, the phase condition on $x,t$ is the single quadratic congruence
$$xt\equiv x^2+c^2\pmod{N^2}.$$
There is no need to test whether $x$ is a unit inside the count: a nonunit cannot satisfy this congruence because $c$ is a unit. For a valid $x$ exactly one $t$ is counted. Lower input and upper phase digits can be expressed with bounded integer variables, and the affine sign filter adds one linear congruence and interval. This gives a fixed-dimensional modular quadratic counting target. A polynomial-time evaluator for it remains open.

The diagnostic verified the exact covariance by enumeration for 30 modulus-divisor choices and compared the local class expansion with direct pair enumeration on seven prime powers, up to $3^6$. Run `python check_affine_interval_variance.py`.
