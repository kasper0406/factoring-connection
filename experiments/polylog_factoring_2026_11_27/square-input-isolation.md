# Unweighted input digit counts for general factoring

Moving the isolation window and pairwise independent signs to the lower input digit extends the square-modulus counting reduction to repeated prime factors. It now covers every odd input with at least two distinct primes. The oracle still supplies two ordinary counts for each upper digit, using simple public bit and interval filters. Estimating the corrected difference to additive RMS error $N/8$ suffices for a randomized polynomial-time split. [Previous phase-digit construction](square-digit-isolation.md)

## Filter the input blocks

Write $x=x_0+Ny$ modulo $N^2$. Keep a block when $(x_0+b)\bmod N<w$, for a public random offset $b$, and give it a pairwise independent sign $\epsilon(x_0)$. Independently choose a uniform unit $c_0$ and set $c=c_0+Nj$. In each block count phases of $x+c^2/x$ whose upper base-$N$ digit is at most $(N-1)/2$. The oracle can count the positive and negative blocks separately. Take their difference, remove the public stationary blocks $x\equiv\pm c_0\pmod N$, and subtract the same corrected value at $j=0$. Divide by $N$ to obtain $F(j)$.

For a fixed block, set $d=\gcd(x_0^2-c_0^2,N)$. Its upper digit runs over one class modulo $d$, with multiplicity $d$, and this class depends only on $j\bmod d$. The normalized difference between two shifts is at most $d/N$. If $d<N$, the contribution has no Fourier coefficient at a unit frequency. This remains true when $d$ contains a partial prime power.

Stationary blocks are $x_0=c_0u\bmod N$ for roots $u^2=1$. The lower window now isolates these input representatives directly. On a singleton trial the remaining stationary block contributes a signed half-interval indicator of a phase with upper-digit slope $2u$. Its coefficient at $2u\bmod N$ has magnitude $1/\pi+O(1/N)$. The affine isolation probability is the same as before because $c_0$ is uniform over units. [Isolation probability](affine-pair-bound.md)

## The squared weights tolerate repeated factors

Signs are independent by input block, so different preimages of the same phase contribute no cross terms to the expected squared norm. Even with the reference subtraction, we have
$$\mathbb E_{b,\epsilon}\|F\|_2^2\le \frac wN\sum_{x_0\in U_N}\left(\frac{\gcd(x_0^2-c_0^2,N)}N\right)^2.$$

The sum factors over complete prime powers. At $p^e$, there are two stationary inputs of weight one. For $1\le v<e$ there are $2(p-1)p^{e-v-1}$ inputs with $v_p(x_0^2-c_0^2)=v$, because the two linear factors have at most one divisible by $p$. There are $(p-3)p^{e-1}$ remaining unit inputs with valuation zero. The local squared sum is therefore
$$2+\frac2p-\frac{p+3}{p^{e+1}}\le2(1+1/p).$$
Consequently the global variance bound is at most $(w/N)Q(1+\omega(N))$, where $Q=2^{\omega(N)}$. This is uniform in the exponents.

## Public trial and decoding

Use the isolation width from the prior proof, with constant $C=8192$: the largest dyadic width at most $N/(CQ)$ if that ratio is at least one, and width one otherwise. The singleton probability is at least $1/(8C)$. In either case $wQ/N\le1$, so the expected squared norm is at most $n+1$ for $n=N.\mathrm{bit\_length}()$. Markov's inequality gives simultaneous singleton isolation and squared norm at most $16C(n+1)$ with probability at least $1/(16C)$.

Clip the normalized oracle answer at $64C(n+1)$. On a good trial the Fourier coefficient changes by at most $1/16$, by the squared norm bound. The oracle's RMS error contributes at most $1/8$ more. A constant-size hidden coefficient survives before dividing by the clipping threshold. The bounded signal then has an inverse-polynomial heavy coefficient. Sampled Fourier recovery and verified gcds of frequencies minus 2 give a split. Trying all the dyadic widths avoids needing $Q$ in advance. Removing even factors and testing for primality and perfect powers handles the other branches of complete factorization.

The reference evaluator can visit just the $w$ lower input blocks, counting the admissible upper digits by floors. It still has a power-of-$N$ cost for wide windows. A polynomial-bit-time evaluator for the noisy count target is unresolved.

## Checks

The diagnostic verified the local squared-weight formula and the unit-frequency identity on five inputs, four with repeated factors. It recovered proper divisors on all five, including products with three and four distinct primes. Direct enumeration modulo $N^2$ cross-checked four upper digits on each small fixture. Run `python check_square_input_isolation.py`.
