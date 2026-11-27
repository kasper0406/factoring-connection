# Weak prediction for translated half windows

For any odd input with exactly two distinct prime divisors, predict the parity contributed by the hidden root pair in a half-length interval. The scale is uniform over units and the interval start is uniform over all residues. A predictor with accuracy $3/4+\delta$ on this joint distribution gives a randomized split in time polynomial in $\log N$, $1/\delta$, and the predictor's running time. The threshold $3/4$ reflects a public shift-only baseline. Accuracy above $1/2$ alone is insufficient for this reduction. [Fixed window exact reduction](fixed-window-split.md) · [Quarter Fourier decoder](quarter-fourier.md)

## Exact public baseline

Let $M=(N-1)/2$, $I_a=a+[1,M]\pmod N$, and let $u$ be a mixed root. Write $G_a(t)=(-1)^{[t\in I_a]+[-t\in I_a]}$. At scale $c$ the hidden sign is $G_a(cu)$ after removing the predictable pair $\pm c$. Over uniform $t\bmod N$ the mean is
$$\mu(a)=1-4(M-J(a))/N,\qquad
J(a)=\max(0,M-|((-2a)\bmod N)-(M+1)|).$$
Indeed $J(a)=|I_a\cap(-I_a)|$ counts pairs of integers in $[1,M]$ with sum congruent to $-2a$. Their sums lie between 2 and $N-1$. The symmetric difference of the intervals has size $2(M-J(a))$.

Across all shifts, the overlap 0 occurs twice, each overlap from 1 through $M-1$ occurs twice, and $M$ occurs once. Summing the resulting arithmetic progressions of absolute values yields
$$\mathbb E_a |\mu(a)|=\frac12+\frac1{2N^2}.$$
Consequently the best shift-only predictor under uniform scales over all residues has accuracy exactly $3/4+1/(4N^2)$. Its label is the sign of $\mu(a)$. For unit scales the shift-only accuracy changes by at most the nonunit fraction.

## Decode excess accuracy by a one dimensional Fourier search

Let $F_a(c)$ be the expected predicted sign, extended by zero from units to all residues. Put $\nu=1-\varphi(N)/N$. If $\nu\ge\delta/4$, uniform gcd probes already split with probability $\Omega(\delta)$, since at most half of the nonunit mass is zero. Run these probes as a parallel option. If $N$ is polynomial in $1/\delta$, trial division suffices. Otherwise the claimed accuracy gives
$$\mathbb E_{a,c}F_a(c)(G_a(cu)-\mu(a))
\ge 2\delta-\nu-1/(2N^2)\ge\delta$$
in the remaining case. This follows by bounding the shift-only contribution by $\mathbb E_a|\mu(a)|$. Since each shift's centered correlation is at most 2, a uniform shift has centered correlation at least $\delta/2$ with probability $\Omega(\delta)$.

For each fixed shift, the template $G_a(t)$ has at most four jumps on the cyclic line. Taking its first difference in the additive Fourier transform bounds every nonzero coefficient by $2/|k|_N$. Centering removes the zero frequency exactly. A cutoff $K=O(\delta^{-2})$ makes the Fourier tail small enough by Parseval. As in the quarter proof, a good shift forces a heavy coefficient of $F_a$ at some frequency $ku\bmod N$, with $0<|k|\le K$. Run the sampled frequency search with the shift held fixed, then test the bounded integers by gcd against the recovered frequencies. A correct frequency gives a verified factor. Sampling polynomially many shifts suffices. The search can also use a direct promise of inverse-polynomial centered correlation in place of the accuracy premise.

The two-prime-divisor hypothesis ensures that there is just one hidden pair. The calculation of a shift-only baseline alone gives no predictor with excess accuracy on unfactored inputs.

## Checks

The diagnostic computed each shift mean by direct enumeration on 78 odd moduli and verified the exact baseline. On three semiprimes it ran sampled Fourier recovery at a fixed shift with planted predictions independently flipped 12% of the time. It recovered a factoring frequency and a divisor in every case, using 25,525 to 32,618 predictor calls. These tests verify the decoder on a planted signal. Run `python check_shift_prediction.py`.
