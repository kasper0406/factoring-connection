# Marked ideals at discriminant 11

For an odd input $N$ coprime to 11, a projective line modulo $N$ specifies a right ideal $I$ of norm $N$ in the fixed quaternion order $\mathcal O$. The classifier computes its ideal class and the matrix modulo $2^k$ of an equivalence from a fixed representative. Its expected cost is polynomial in $\log N+k$. Both the splitting and the ideal calculation use the public integer. [Order, splitting, and line construction](level11-sampling.md) · [Implementation](level11_marking.py)

## Computing the transport

Use the basis $(1,\omega,j,\omega j)$ of the four-unit order $\mathcal O$. Take $\mathcal O$ as the principal representative $P$. For the six-unit class $A$, take the norm-three right ideal $J$ defined by $x+y+z=z+w-x=0\bmod 3$. In the splitting $\omega=\operatorname{diag}(1,0)$, $j=\left(\begin{smallmatrix}0&1\\-1&0\end{smallmatrix}\right)$ at 3, this is the line $(1,1)$. All norm-three elements of $\mathcal O$ have image on a coordinate axis.

Let $n$ be the absolute norm of the input ideal $I$. A principal ideal has four generators of norm $n$. Otherwise enumerate the norm-$3n$ elements $\beta\in I\bar J$. Local principality gives $\bar J J=3\mathcal O$. Thus $\beta J\subseteq 3I$, and equality of norms gives $(\beta/3)J=I$. There are six such transports, differing by right multiplication by units of the left order of $J$. This also works when $3\mid n$.

For either target class, write its representative norm as $d\in\{1,3\}$. The transport lattice has covolume $(dn)^2$ in order coordinates, and each nonzero vector has norm divisible by $dn$. The norm form is bounded above and below by fixed multiples of squared Euclidean length. These bounds and four-dimensional LLL put the shell of norm $dn$ into a coefficient box of constant size. Computing its integer basis and all transports takes polynomial bit time. [Lattice bound and splitting cost](level11-sampling.md)

The default base is $P$. To start at $J$, use $\gamma=(1,-1,0,1)\in J$, of norm 6. The index of $\gamma\mathcal O$ in $J$ is four, so $\gamma$ generates $J$ at every odd prime. The ideal $\gamma I+NJ$ gives the relative norm-$N$ subideal of $J$ specified by the input line. Its absolute norm is $3N$, and its lattice at 2 is $J_2=\mathcal O_2$.

## The five states modulo 2

The three nonzero vectors over $\mathbb F_2$ are permuted by the transport. The projective units of $P$ fix $(1,1)$, leaving three possible images. The projective units of the left order of $J$ form the alternating subgroup, leaving two possible permutation signs. Number the principal images by $0,1,2$ in the order $(1,0),(0,1),(1,1)$, and the two other states by $3,4$ for even and odd. Each marked state has stationary mass $1/5$.

For a good prime $r$, let $B_r(i,j)$ count the $r+1$ line choices taking state $i$ to state $j$. Write $a_{11}(r)$ for the coefficient of $\eta(z)^2\eta(11z)^2$, $a_{44}(r)$ for the trace of the newform of level 44, and $s=(0,0,0,1,-1)$. Then

$$
B_r(i,j)=\frac{r+1}{5}+a_{11}(r)\left(\mathbf 1_{i=j}-\frac15\right)
+\frac{a_{44}(r)-a_{11}(r)}2 s_i s_j.
$$

The three individually labelled 2-neighbors of a terminal class have principal-bit pattern $111$ at $A$, or one of the three unit vectors at $P$. Indeed the unmarked Brandt rows at 2 are $(1,2)$ and $(3,0)$. If $C$ is the central principal bit and $C_i$ the bit at neighbor $i$, the four pattern indicators are $1-C$ and $C+C_i-1$. Coprime ideal modifications commute for each individual line. Each labelled neighbor therefore obeys the same single-prime marginal as $C$ from its transported base. This makes the three-dimensional space of mean-zero functions of the four patterns an eigenspace with eigenvalue $a_{11}(r)$. [Unmarked marginal](level11-sampling.md)

The remaining contrast is $s$. At 2, conjugation by $\operatorname{diag}(1,2)$ takes the principal congruence subgroup $K(2)$ onto $K_0(4)$. The cuspidal Brandt module therefore corresponds to weight two at level $44=4\cdot11$, new at 11. This space has dimension four, with three translates of the level-11 newform and one level-44 new eigenline. Thus $s$ has eigenvalue $a_{44}(r)$. Projecting onto constants, the old space, and $s$ gives the displayed matrix. [Pizer's Brandt correspondence](https://wstein.org/msri06/refs/Pizer-An_Algorithm_for_Computing_Modular_Forms_on_Gamma0N.pdf) · [Level-44 newform](https://www.lmfdb.org/ModularForm/GL2/Q/holomorphic/44/2/a/a/)

At $N=19\cdot29$, both level-11 coefficients vanish. The level-44 elliptic curve gives $a_{44}(19)=8$ and $a_{44}(29)=0$. Given a 29-line of state 3, its twenty lifts have state counts $(4,4,4,8,0)$. State 4 interchanges the last two counts. Every 19-line has six lifts of each state. Starting instead from state 3, the expected signed value on a uniform $pq$-line is $a_{44}(p)a_{44}(q)/((p+1)(q+1))$ by successive applications of the two prime transitions. [Curve for the point counts](https://www.lmfdb.org/EllipticCurve/Q/44/a/1)

## Integral frames at higher powers of two

For $M=2^k$, $k\ge2$, reduce each transport $\beta/d$ in the fixed splitting $\mathcal O/M\simeq M_2(\mathbb Z/M)$. Return the target class and the lexicographically least row-major matrix. Enumerating the four or six witnesses gives the full right coset of the reference unit group. The coordinate splitting uses $O(k)$ bit lifts, and the four-by-four matrix work is polynomial in $k$.

Every output matrix has determinant $n/d\bmod M$. Within a fixed absolute determinant sector $n\bmod M$, the principal class contributes $|\mathrm{SL}_2(\mathbb Z/M)|/4$ states and the other class contributes $|\mathrm{SL}_2(\mathbb Z/M)|/6$. The units are already distinct modulo 4. All marked stabilizers are trivial, so these states have equal mass. The sector size is $5M^3/16$, using $|\mathrm{SL}_2(\mathbb Z/2^k)|=3M^3/4$.

From the repository root run `python experiments/polylog_factoring_2026_11_27/check_level11_marking.py`. The checks compare every five-state output with independent integral witnesses on 1,040 CRT cells, test the signed sums from the nonprincipal base, compare full frames with an independent 2-adic splitting, and check determinants at $k=128$. The public entry points are `marked_line_classifier(n)` and `framed_line_classifier(n,k)`. Each returns a callable on a primitive coordinate pair $(u,v)$ and also accepts `base="nonprincipal"`. [Exact checks](check_level11_marking.py)
