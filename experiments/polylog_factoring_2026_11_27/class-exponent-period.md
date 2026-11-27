# Recurrence size over the full odd class orbit

The fixed-discriminant family has a cyclic subgroup of order $L$. Index its odd powers by $2b+1$ and let $g(b)$ be the signed disk count with the axis and the two public roots removed, divided by two and reduced modulo two. If the odd input has at least two distinct prime factors, every stationary recurrence over $\mathbb F_2$ for this full periodic sequence has order greater than $L/4$. The reduction itself asks only for a subset of this orbit, so the bound is for a stronger evaluator. [Class encoding](cyclic-class-queries.md) · [Disk parameters](fixed-discriminant-split.md)

## Half an orbit exposes a boundary point

Take the centered odd shear $j\bmod L$ for each query form $(Lx+jy)^2+Dy^2$. Let $H=NB_0$, $M=(N-1)/2$ and $K=\lfloor M/2\rfloor$. The construction gives
$$LK\le\sqrt{H-D}<LK+L/8.$$
Indeed $H-D$ lies between $(LK+2N)^2$ and $(LK+2N)^2+3N$, and $L>16N^2$. Only $y=\pm1$ can contribute away from the axis.

Adding $L/2$ to the class exponent adds $L/2$ to the ideal parameter $t$ modulo $L$, because that is the unique element of order two and the parameter addition law reduces to ordinary addition for this element. Inverting an odd parameter to get the shear also adds $L/2$ to the shear modulo $L$. Choose $j\in[L/8,L/4)$ odd. The centered paired shear is $j'=j-L/2$. At $y=1$ the allowed integers run from $-K$ through $K-1$ for $j$, and from $-K+1$ through $K$ for $j'$, by the radius bound.

Reading $x+j/L$ modulo $N$ means $x+jL^{-1}\bmod N$. The two intervals of root residues have length $2K$ and start at $a=jL^{-1}-K$ and $a+(M+1)$ modulo $N$. The possible $j$ in the chosen interval realize every $a$ modulo $N$ because $L/8>2N$. The parity difference of the intervals equals the mixed-root parity on their complement, since the total number of mixed roots is even.

If $M$ is even, the complement is just $a+M$. A shift placing a mixed root there makes the parity difference one. If $M$ is odd, the complement consists of offsets $M-1,M,2M$ from $a$. Multiplying by 2 modulo $N$ turns those offsets into the three consecutive residues $N-3,N-2,N-1$. If the parity vanished for every $a$, the induced binary root indicator in the doubled coordinate would satisfy $f(x)+f(x+1)+f(x+2)=0$ everywhere. Subtract consecutive equations to get period 3. When $3\nmid N$, this would make it constant. When $3\mid N$, it vanishes at $0,2,-2$, because zero and the two public roots were removed from the undoubled indicator. These three residues cover all classes modulo 3. In either case $f$ would vanish, contradicting the existence of mixed roots. Thus some pair has different parities.

## Stationary recurrence consequence

The sequence in $b$ has period $P=L/2$, and the differing pair shows that $P/2$ is not a period. Its minimal shift annihilator divides $X^P-1=(X+1)^P$ over $\mathbb F_2$. If its degree were at most $P/2$, it would divide $(X+1)^{P/2}=X^{P/2}-1$, forcing period $P/2$. Hence its degree, and therefore the dimension of a stationary linear representation, exceeds $P/2$. Since $L=\Theta(N^2)$ for the stated construction, the bound is exponential in the bit length.

The check script independently enumerated signed disks and verified exponent differences to find witnesses on all 84 odd inputs with mixed roots through 319. The other 74 moduli had no mixed roots. Run `python check_class_exponent_period.py`.
