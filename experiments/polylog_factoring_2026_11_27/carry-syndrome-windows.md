# Predicting the carry syndrome with an interior word missing

Knowing the bits on **both sides** of a missing word of the product still
leaves little chance of finding a factor with the binary-polynomial
decoder on a random balanced prime pair. This estimate counts every mask
that passes the actual divisor test, including masks far from the true
carry word. It allows the predictor to know which prime is being targeted
when it constructs its table. [Decoder](carry-parity.md) ·
[Bound when only the high prefix is given](carry-syndrome-bounds.md)

## Statement

Write $R=2^{b-1}$, take distinct uniform primes $p,q\in[R,2R)$, and pad
$N=pq$ to $2b$ bits. For $0\le k\le b-2$, hide the $b+k$ bits at positions
$s,\ldots,s+b+k-1$, with $0\le s\le b-k$. For each $p$ and each assignment
of the outside bits, allow an arbitrary set of at most $M$ target residues
modulo $U_p$, the binary digit polynomial of $p$ over $\mathbf F_2$.
The success event is $U_N\bmod U_p$ lying in that set. Put

$$
A=\log_2(4/3),\qquad D=2^k(3/4)^{b-k}.
$$

Uniformly in the window position and the choice of tables, the probability
of success is at most
$$
R^{o(1)}\bigl(M2^{-k}+D^{2/3}\bigr).
$$
The nontrivial regime has $D=o(1)$. For $M=2^{o(b)}$, choose
$k=\lfloor 2Ab/(5+2A)\rfloor$. The bound becomes
$N^{-A/(5+2A)+o(1)}=N^{-0.0711890\ldots+o(1)}$. Here $\log_2N=2b+O(1)$.
A union over the $O(b)$ possible window positions, or over the choice
of either prime, preserves the exponent. A menu of corrections produces
at most as many target residues as corrected masks. Extra information
taking at most $J$ values can be included by multiplying $M$ by $J$.
For the lowest window, the high-prefix argument gives a stronger exponent.
[High-prefix theorem](carry-syndrome-bounds.md)

## Collisions inside a window

First take $p$ uniform over odd $b$-bit integers and allow all integer
partners $R\le q<2R$. Fix $p$. For every outside word and every selected
syndrome that occurs, set aside its largest successful partner. This sets
aside at most $M2^{b-k}$ partners. Every remaining successful $q$ has a
larger partner $q'$ in the same class. Put $B=2^s$ and let

$$
K=\bigl((pq)\mathbin{\mathrm{XOR}}(pq')\bigr)/B,\qquad h=(q'-q)/B.
$$

Both divisions are exact. The shared low bits give $B\mid p(q'-q)$, and
$p$ is odd. The shared high bits give $pBh<B2^{b+k}$, so
$1\le h<2^{k+1}$. Equal syndromes imply $U_p\mid T^s U_K$. Since
$U_p(0)=1$, we have $U_K=U_pV$ for some nonzero polynomial $V$ of degree
at most $k$. Writing $x=\lfloor pq/B\rfloor$ yields the exact integer
identity

$$
2(x\mathbin{\mathrm{AND}}K)=K-ph .
$$

For fixed $p,V$, let $a$ count the ones in the lowest $k$ bits of $K$,
and let $w$ count those in its lowest $b-1$ bits. The left side modulo
$2^{k+1}$ takes at most $2^a$ values, so there are at most $2^a$ possible
$h$, since $p$ is odd. Each feasible $h$ fixes $x$ at every set bit of $K$.
Consequently the possible $z=x\bmod R$ lie in a set $E_{p,V}$ of at most
$R2^{a-w}$ residues. Pool the sets over all $V$ in $E_p$ *before* counting
how often partners land in them.

For fixed $V=T^vV_0$, $v\le k$ and $V_0(0)=1$. Multiplication by $V_0$
permutes the odd polynomials modulo $T^{b-1-v}$. Thus as the middle digits
of $p$ vary, the bits of $K$ in positions $k+1,\ldots,b-2$ are independent
fair bits. Each contributes $3/4$ to the expectation of $2^{a-w}$. Summing
over fewer than $2^{k+1}$ nonzero $V$ gives the position-independent bound
$$
\sum_{p\ \mathrm{odd}} |E_p| \ll R^2 D .
$$

## How often a shifted product word occurs

Define $n_p(z)=\#\{q\in[R,2R):\lfloor pq/B\rfloor\equiv z\pmod R\}$.
For any integer $T\ge5$, the number of odd $p$ having some $n_p(z)\ge T$
is $O(R^{1+o(1)}/T^2)$. The window constraint gives $1\le B\le2R$.

To prove the bound, put $H=BR$. These partners are lattice points
$(q,t)$ with $t\equiv pq\bmod H$ in the rectangle
$[R,2R)\times[Bz,B(z+1))$. The lattice has determinant $H$ and the
rectangle has area $H$. A noncollinear set of $t_0$ points can be
triangulated into at least $t_0-2$ nondegenerate triangles using all its
points. Each lattice triangle has area at least $H/2$. Hence five or more
points in this rectangle must be collinear.

Order at least $T$ such points by their first coordinates. Collinearity
makes the shortest successive increment $(v,u)$ satisfy
$$
1\le v\le R/(T-1),\qquad 0<|u|\le B/(T-1),\qquad pv\equiv u\pmod H .
$$
Indeed $u=0$ would give $H\mid v$, since $p$ is odd and $H$ is a power
of two. The nonzero increment yields $B\ge T-1$. Write $m=pv$.
Then $R\le m\le 2R^2/(T-1)$ and $m$ lies within $B/(T-1)$ of a
*positive* multiple of $BR$ because $B/(T-1)<R$. The existence of that
multiple implies $B\ll R/T$. There are therefore at most
$$
O\bigl((1+R/(BT))(1+B/T)\bigr)=O(R/T^2)
$$
possible integers $m$. Each has $R^{o(1)}$ possible divisors $p$ by the
standard divisor bound, as $m\le2R^2$. This proves the occupancy estimate.

Each exceptional $p$ has at most $R$ partners. For every other $p$ there
are at most $T|E_p|$ partners whose shifted word belongs to $E_p$. Restoring
the set-aside partners and dividing by the number of ordered prime pairs
$\gg R^2/b^2$, we obtain
$$
O(b^2)\left(M2^{-k}+TD+R^{o(1)}/T^2\right).
$$
Taking $T=\max(5,\lceil D^{-1/3}\rceil)$ proves the statement for $D=o(1)$.
The trivial probability bound covers larger $D$. No assumption about
equidistribution of primes in short intervals is needed.

## Finite checks

The script exhausts all odd $p$ and all integer $q$ in the indicated
binary ranges, at every permissible offset. It checks the gaps, the AND
identity, the cardinality of each pooled set's component, the multiplicative
bit average and the Cauchy–Schwarz incidence inequality
$(\sum_p\sum_{z\in E_p}n_p(z))^2\le(\sum_p |E_p|)(\sum_{p,z}n_p(z)^2)$.
It also computes the best possible one-target table at each position.
The stronger collinearity estimate above is proved analytically.
[Run the checker](check_carry_syndrome_windows.py)

| $b$ | $k$ | Positions | Successive collision pairs checked | Min–max hits of best table | All $(p,q)$ per position |
|---:|---:|---:|---:|---:|---:|
| 5 | 1 | 5 | 107 | 83–128 | 128 |
| 6 | 1 | 6 | 404 | 288–512 | 512 |
| 7 | 2 | 6 | 2,292 | 820–1,182 | 2,048 |
| 8 | 2 | 7 | 8,323 | 2,852–4,518 | 8,192 |
