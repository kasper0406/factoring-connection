# Nonadaptive unweighted disk reduction for general factoring

An odd input with at least two distinct prime divisors can be split with constant probability by $O(n^3)$ nonadaptive signed disk queries modulo 4, with $n=N.\mathrm{bit\_length}()$. Every query for that input uses one common numerical cutoff and a primitive positive form. The reduction uses affine isolation and reads an isolated coordinate by CRT. It replaces the weighted scalar at the cost of more queries. [Affine isolation probability](affine-pair-bound.md) · [Common cutoff construction](affine-isolation.md)

## Read the isolated coordinate modulo small primes

An affine trial maps roots $u^2=1\pmod N$ to $z=au+b\bmod N$ for a public unit $a$. It retains the positions $0\le z<w$. The isolation theorem gives a batch of $O(n)$ such triples $(a,b,w)$ with constant probability that at least one trial retains exactly one mixed root.

Choose distinct auxiliary primes $\ell$ coprime to $N$ whose product exceeds $N$. They can all be chosen below $O(n)$. By the prime number theorem, the sum of logarithms of primes up to a sufficiently large constant times $n$ exceeds $2\log N$, and skipping primes dividing $N$ loses at most $\log N$. Trial enumeration and primality testing of these auxiliary integers take polynomial bit time.

For each trial and each $\ell$, query the parity of the retained mixed roots in every residue class $r\bmod\ell$. When the trial has one root, exactly one parity is odd for each prime. These residues recover its coordinate by CRT with no ambiguity in $[0,w)$. Then $u=a^{-1}(z-b)\bmod N$, and $\gcd(u-1,N)$ verifies a split. A trial with multiple retained roots is allowed to produce a candidate. The same checks guarantee that every returned divisor is proper.

## Encode a residue class by a shorter affine window

Write $z=r+\ell v$ in the residue class. If $r<w$, the range is
$$0\le v<w_r,\qquad w_r=1+\lfloor(w-1-r)/\ell\rfloor.$$
Set $a'=a\ell^{-1}\bmod N$ and $b'=(b-r)\ell^{-1}\bmod N$. The query is just a length-$w_r$ affine window for $v=a'u+b'\bmod N$. Because $w_r\le w$, the same padded form construction and numerical cutoff apply, including when $w_r$ is non-dyadic. Empty residue classes need no query.

Use one reference form with the same leading coefficient and cutoff and only an axis in its disk. Subtract its signed count modulo 4 from the query count. The $y=1$ and $y=-1$ rows then contribute twice the affine root count. Correct the predictable endpoints and roots $u=\pm1$ as in the affine weighted construction, using counts in place of weights. The axis cancels for repeated factors too.

At most $\sum_\ell \ell=O(n^2)$ disk queries are needed per trial, plus the one reference reused across all trials. All query choices depend only on the public modulus, samples and auxiliary primes. The implementation can stop once it verifies a factor. Repeating batches improves the success probability, and recursion with primality and exact perfect-power tests gives complete factorization given a polynomial-bit-time disk evaluator. That evaluator remains open.

## Checks

Direct signed disk enumeration exercised 1,920 seeded trials on ten inputs, including repeated factors and products of four distinct primes. The decoder recovered and verified a factor from every one of the 335 singleton trials. The run made 14,763 signed count calls with early rejection and checked the shared cutoff for every planned query. Run `python affine_parity.py`.
