# A fixed cutoff suffices for the general disk reduction

Complete integer factorization reduces to signed disk counts modulo 4 at one public numerical cutoff per odd input. Each query changes the form while retaining that cutoff. The underlying parity oracle uses a window of fixed length and changes only its starting point. The same construction handles moduli with repeated factors. Computing the disk residues in polynomial bit time remains open. [Earlier fixed-form and varying-cutoff reduction](https://github.com/kasper0406/factoring-connection/blob/9514cbfa013cd9ae3f030044db1dc3fb6bb4a9ed/experiments/ellipse_disk_mod8/native_dyadic.md)

## Bisection with a fixed window length

For odd $N\ge5$, set $M=(N-1)/2$. Let $h(x)$ indicate a root of $x^2=1\pmod N$ other than $\pm1$. Let $P(t)=\sum_{x=1}^t h(x)\bmod2$ and ask for
$$W(a)=\sum_{x=a+1}^{a+M}h(x)\bmod2.$$
Reflection $h(x)=h(N-x)$ gives $W(a)=P(a)+P(M-a)\pmod2$ for $0\le a\le M$. If $N$ has at least two distinct prime divisors, then $W(0)=1$ because the roots pair off under negation, leaving $2^{\omega(N)-1}-1$ mixed roots in the lower half.

At $a=\lfloor M/2\rfloor$, the two prefixes cancel if $M$ is even. If $M$ is odd, their difference is just $h(a+1)$, which can be tested directly. If that test succeeds, $\gcd(a,N)$ splits $N$. Otherwise the two endpoints have different parity. Bisect without a monotonicity assumption, retaining endpoints of different parity. For adjacent endpoints $a,a+1$, the parity changes exactly when one of $a+1$ and $a+M+1$ is a mixed root. Test both and take the gcd of the root minus 1 with $N$. This takes $O(\log N)$ queries. Prime powers have no mixed root and can be handled by exact integer root tests and primality testing. Remove even factors directly and recurse on verified divisors.

## Signed disks at one numerical cutoff

Put $K=\lfloor M/2\rfloor$, $A=(N+1)^4$, $R=N(N+1)$ and $C_0=A(K+1)^2$. For window start $a$, set $s=a+K+1$. Choose the cross coefficient $b$ in $[0,R)$ by CRT so that $b\equiv2s\pmod N$ and $b\equiv1\pmod{N+1}$. Choose $c\in[C_0,C_0+N)$ so that $c\equiv s^2-1\pmod N$. Use
$$Q(x,y)=Ax^2+bxy+cy^2,\qquad Q'(x,y)=Q(x,y)+2NB_0y^2,$$
with the common cutoff
$$B_0=\left\lceil\frac{C_0+AK^2+RK+N}{N}\right\rceil.$$
Both forms are primitive because $\gcd(A,b)=1$. Their coefficients and cutoff have $O(\log N)$ bits.

Writing $H=NB_0$, we have $H<2C_0$, since $A(2K+1)>RK+2N$. At $y=1$, every $|x|\le K$ has $Q(x,1)\le H$. If $|x|\ge K+1$, then $Q(x,1)\ge C_0+A(K+1)^2-R(K+1)>H$, since $(A-R)(2K+1)>2N$. For $|y|\ge2$, completing the square gives $Q(x,y)\ge4C_0-R^2/A>H$. The second disk has only the axis, which cancels exactly against the same axis in the first count. This avoids any assumption on the nonzero axis solutions for nonsquarefree $N$. The cutoff is below the usual $\sqrt{|\operatorname{disc}|/3}$ bound for both forms, since $B_0<2C_0/N<A(K+1)$ and $4Ac-b^2>3A^2(K+1)^2$.

Modulo $N$, both forms reduce to $(x+sy)^2-y^2$. The remaining rows $y=\pm1$ therefore contribute twice the number of roots in $[a+1,a+2K+1]$, by the symmetry $(x,y)\mapsto(-x,-y)$. The signed count difference modulo 4 yields that parity. If $M$ is even, check and remove the extra endpoint at $a+M+1$. Remove the known roots $\pm1$ when present to obtain $W(a)$. Every split uses just the same $B_0$ for its current modulus.

## Checks

The reference enumerator computed signed disks directly, including the axis. The check covered 6,478 windows over all 158 odd integers from 5 through 319. It verified the common cutoff, primitive coefficients and discriminant bound, split all 84 inputs with two or more distinct prime factors, and returned no split for the 74 primes or prime powers. Run `python fixed_window_split.py`. Reference enumeration serves as a correctness control and has a power-of-$N$ cost.
