# Discrepancy for interval and periodic quarter predictors

On balanced semiprimes, any quarter-window predictor whose expected output is periodic with small period on each of a small number of integer intervals has exponentially small advantage. The bound holds even when its interval endpoints and residue labels depend arbitrarily on the factors. It concerns the hidden-pair bit after subtracting the public pair. [Weak prediction reduction](quarter-fourier.md)

Let $N=pq$ with $p<q\le2p$. For a mixed root $u$ write $u=aq+s$ with $s=\pm1$. Then $a$ is a nonzero residue modulo $p$, since $u$ has the opposite sign at $p$. The quarter bit in sign form is $g(c)=1$ if the distance of $cu/N$ to the nearest integer is at most $1/4$, and $-1$ otherwise.

In any progression $c=c_0+dj$ with $p\nmid d$, take a consecutive block of $p$ terms. Starting at the block's first phase, their phases modulo one have the form
$$\theta + \frac{da j}{p}+\frac{dsj}{N},\qquad 0\le j<p.$$
The first two terms describe a rotated uniform $p$-point grid. Its signed sum for the quarter rule has absolute value at most 2. The final term moves any grid point by less than $d/q$. Its sign can change only if the grid point is within $d/q$ of a quarter-window endpoint. The two endpoint neighborhoods contain at most $4dp/q+4$ grid points, so the actual block sum has absolute value at most $10+8dp/q$.

Now partition the integer scales in $[0,N)$ into at most $S$ intervals. On each interval let a predictor's expected sign depend only on $c\bmod d$, where the period can vary by interval and $1\le d\le D<p$. In an interval of length $T$, summing over the residue classes uses at most $T/p$ complete blocks and at most $dp$ positions in incomplete blocks. Multiplying the class sums by their respective predicted signs in $[-1,1]$ yields a correlation numerator of magnitude at most
$$10T/p+8dT/q+dp.$$
Sum over the intervals. Extending the predictor to nonunits for this bound and then removing them changes the numerator by at most $p+q-1$. The correlation on unit scales therefore has magnitude
$$O\left(\frac{q+SDp}{\varphi(N)}\right)=O\left(\frac{SD}{\sqrt N}\right).$$
Prediction advantage is half that correlation. Polynomially many intervals and polynomially bounded periods in the input bit length cannot supply the inverse-polynomial advantage needed by the reduction. The conclusion leaves representative-dependent operations with high-period or nonperiodic structure open.

## Adding the Jacobi symbol

The same conclusion holds if each label also sees the Jacobi symbol $\chi(c)=(c/N)$. On unit scales, decompose a label into its average and half-difference at the two Jacobi signs. The half-difference is multiplied by $\chi(c)$. For each arithmetic progression cell, let $I$ be its indicator. Its normalized additive Fourier coefficients have total absolute value $O(\log N)$ by the geometric-sum bound, since its step is a unit and permutes the frequencies. The same bound holds for $g$, since multiplying its argument by $u$ permutes the quarter rule's frequencies. The Jacobi character is primitive for a product of two distinct primes. The Gauss-sum bound gives $|\widehat\chi(k)|\le N^{-1/2}$ at every additive frequency. Thus the normalized correlation on a cell with the Jacobi factor is at most $O(N^{-1/2}\log^2 N)$ by Fourier convolution. The nonunit terms vanish because $\chi=0$ there. Summing over at most $SD$ cells gives advantage at most $O(SD\log^2 N/\sqrt N)$.

The diagnostic checked the complete-block discrepancy by direct evaluation on three semiprimes. Run `python check_quarter_periodic_features.py`.
