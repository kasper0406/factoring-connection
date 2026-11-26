# A centered recurrence on quadratic forms

An unmasked signed count modulo 16 suffices to factor a balanced semiprime in $O(\log N)$ queries **if the evaluator accepts the explicitly transformed forms below**. Every queried region is centered at zero. The congruence and the ellipse change together. The construction gives a working evaluator and an adaptive reduction. The native evaluator still uses power-of-$N$ work.

The [earlier decoder](README.md) held the form fixed and added predicates on bits of $y$. Here the query is always the same scalar

$$R_N(F,B)=\sum_{\substack{z\in\mathbb Z^2\\0<F(z)\le NB\\N\mid F(z)}}1\pmod {16}.$$

The forms along a path can be nonprincipal and their discriminants can grow. This condition matters if an implementation only covers the original fixed principal form. The [previous obstruction](translation.md) concerned translates of one periodic weight. The recurrence below changes the lattice, at the same center.

## Exact count over the whole disk

Write $F=[A,b,C]$ for $Ax^2+bxy+Cy^2$. Every primitive integer direction has one of three slopes modulo 2. If $y$ is odd, use $x\equiv r y\bmod 2$, with $r=0,1$. Otherwise $x$ is odd, and use $y\equiv 0x\bmod 2$. At stride $m=2^k$ each occupied slope lifts to $r$ or $r+m$ modulo $2m$, in the same chart. The two children partition the *primitive directions*. Their full lattices can intersect at nonprimitive vectors.

In the first chart take the basis $(m,0),(r,1)$ of the centered lattice $L$. The form on that lattice is

$$F_L=[Am^2,\ m(2Ar+b),\ Ar^2+br+C].$$

The other chart swaps the two physical coordinates. Before primitive normalization the discriminant is $m^2\operatorname{disc}(F)$. At each node Gauss reduction gives the first two squared minima $a_L,c_L$ and a shortest vector $v$ in physical coordinates. Put $H=NB$.

* If $a_L>H$, the node is empty.
* If $a_L\le H<c_L$, every nonzero point of this disk in $L$ lies on the line of $v$. When $v$ is primitive in the ambient lattice, the entire contribution of this direction is
  $$2\left\lfloor\frac{\lfloor\sqrt{H/a_L}\rfloor}{N/\gcd(N,a_L)}\right\rfloor.$$
  Indeed $N\mid k^2 a_L$ iff $N/\gcd(N,a_L)\mid k$ for squarefree $N$. When $v$ is imprimitive, that node owns no ambient primitive direction meeting the disk.
* Otherwise visit both slope lifts, or all three charts at the root. Choose strictly $c_L>H$ for the ray test so that equal minima and boundary points survive.

Every nonzero integer vector is assigned to the unique leaf containing its ambient primitive direction. That leaf counts **all** its signed multiples at once. This also accounts for even vectors belonging to several centered lattices. The [implementation](centered.py) uses only the public form, $N$ and $B$ and computes the residue without constructing a mixed root or looping over norm indices.

**Cost of the actual evaluator.** Independent vectors of norm at most $H$ have absolute determinant at most $2H/\sqrt{|\operatorname{disc}F|}$, by the Gram determinant identity. Their determinant in a stride-$m$ node is a multiple of $m$. Consequently no node branches beyond this bound, and there are at most $3m/2$ charts at any positive depth. The total node count is $O(1+H/\sqrt{|\operatorname{disc}F|})$. Arithmetic per node costs polynomially many bit operations in the input size. The depth-first stack uses logarithmically many nodes.

For the public principal form of discriminant $-\Delta$ this node bound has the right order in the large N-only branch regime. Inside the disk take the rectangle with $1\le x\le\lfloor\sqrt H/4\rfloor$, $1\le y\le\lfloor\sqrt{H/\Delta}/2\rfloor$. A positive constant fraction of its points are primitive: the share with a common divisor is at most $\sum_{d\ge2}1/d^2<1$. They are distinct directions. Each requires its own ray leaf. Thus the native geometric traversal takes $\Theta(H/\sqrt{\Delta})$ nodes when both sides grow. Since $\Delta<4\sqrt N+4$ on the N-only branch, this is already $\Omega(N^{3/4})$ at $B=1$. This bound is for the implemented traversal. It includes unsuccessful rays and does not rely on a timing fit.

## Adaptive reduction using the unmasked scalar

Start with the repository's [public branch](../ellipse_aggregate/README.md) constructed from $N$ alone, or a supplied eligible branch. Its principal form is $Q=[1,b,(\Delta+b^2)/4]$, $b=\Delta\bmod2$. Let $s$ be the first mixed norm divided by $N$. On the original form, double and bisect to find $s$ as in the [unmasked radial reduction](../ellipse_global/README.md). The public contribution is $4\lfloor\sqrt B\rfloor$ in the search range. Before $2s$ at most three mixed sign/conjugation orbits are present. A nonempty prefix there has signed mixed residue 4, 8 or 12.

At $B=s$ there are at most eight signed mixed vectors, including ties. Each is primitive. A common prime divisor of its coordinates either divides $N$, impossible at $s<\min(p,q)$, or can be divided out to give a smaller positive mixed norm. Thus each of these vectors belongs to exactly one primitive child of the current node.

For each candidate child request $R_N(Q_L,s)$. Subtract all public signed points in that lattice through $s$. For each public direction $(x,y)$ the positive multiples in the first chart have step $m/\gcd(m,x-r y)$, and similarly with $x,y$ exchanged in the second chart. This gives an exact public subtraction by two integer divisions and the two antipodal signs. The residual is between 0 and 8 as an integer, so its residue modulo 16 detects occupancy. Choose one nonempty child. Its siblings need no further descent. Once only one geometrical direction is possible, Gauss reduction gives the primitive mixed vector $(x,y)$ directly. The final gcd is

$$\gcd\left((cy-\ell(2x+by))/2,\ N\right).$$

The determinant bound gives $O(\log N)$ levels with at most three queries on the first level and two thereafter. The search for $s$ also costs $O(\log N)$ queries. All remaining arithmetic has polynomial bit cost. This proof uses the public eligibility and cap from the existing reduction.

**What a general evaluator must accept.** If $g$ is the content of $Q_L$, it is a power of 2 because the substitution is invertible at every odd prime and $Q$ is primitive. Divide by $g$, query at $\lfloor s/g\rfloor$, and optionally Gauss-reduce the result. This preserves the scalar since $N$ is odd. A public square root of the new discriminant modulo $N$ is $m c/(g\ell)$ up to sign. The queried forms have $O(\log N)$-bit coefficients and discriminants on this branch. They need not have a known vector of norm $N$, nor lie in the principal class. An algorithm that requires those extra promises cannot be substituted without a further reduction. The CLI implements the complete factorer using the native evaluator for each request.

## Can equivalence or local representation remove most children?

There is a concrete obstruction already at the single cutoff $B=1$. Consider an N-only branch with $\Delta\equiv3\bmod8$ and $\Delta>4$. Put $\omega=(1+\sqrt{-\Delta})/2$, $\mathcal O=\mathbb Z[\omega]$, $r_0=(c-1)/2\bmod8$. Take a power of two $m\ge8$ of order $\sqrt N/2$ (the benchmark rounds down). Use the first chart $r\equiv r_0\bmod8$. For all sufficiently large such inputs the following $m/8$ children survive three independent tests.

1. **Geometry.** The two independent vectors $(r,1),(r-m,1)$ are in the disk for $0\le r<m$ since $m\le\sqrt N/2$ and $\Delta=O(\sqrt N)$.
2. **Every prime power.** At an odd prime power the basis substitution is invertible, and the parent already has a public vector of norm $N$. At 2 the odd ratio $N/Q(r,1)$ is 1 modulo 8, since $Q((c-1)/2,1)=N$. Such a unit is a square at every 2-power: from precision $2^j$, $j\ge3$, changing an odd root by $2^{j-1}$ toggles its square's next bit. Multiplying $(r,1)$ by the resulting local root solves the norm equation.
3. **Gauss equivalence, including reflection.** The lattice is $I_r=\mathbb Z(r+\omega)+m\mathcal O$. Since the norm of $r+\omega$ is odd and $m$ is a 2-power, $\mathcal O I_r=\mathcal O$ by Bézout. Any isometry between two such norm lattices extends over their common rational plane to multiplication by a norm-one element, possibly followed by conjugation. This follows by dividing the isometry by its value at 1: the resulting isometry fixes 1 and acts by sign on its perpendicular line. Taking $\mathcal O$-spans forces the multiplier to be a unit of $\mathcal O$. Its only units are $\pm1$ since $\Delta>4$. Thus equivalence forces the slopes to be $r$ or $-1-r$ modulo $m$. The second is outside the selected residue modulo 8. All $m/8$ forms are inequivalent.

On the 35-bit input $N=33,119,086,993$, $\Delta=3243$ and $m=65,536$. All 8,192 selected forms pass. Exactly two represent $N$ over the integers. One has the public representation. The other contains the mixed one, at $r=12,122$. This fixture was deliberately selected for a mixed hit at $B=1$. On the pre-existing 36-bit input $N=38,942,974,393$, the corresponding 8,192 forms contain only the public hit. The global labels come from a separate exact square test on the single audit shell and are never supplied to a factorer. Finite checks include direct reduction of every form and an independent mod-2048 test. The all-prime-power claim follows from the proof. This result concerns these pruning rules; an aggregate algorithm could compute many occupancy bits together.

## Verification

The [controls](centered_check.py) compare every distinct requested residue against a separate integer row enumeration. They check the cutoffs $s-1,s$ explicitly, and the end-to-end factor output. [Recorded results](centered_results.json).

| Cohort | Branches | Distinct residue checks within branches | Largest number of scalar queries for factoring | Native nodes across factor runs |
|---|---:|---:|---:|---:|
| Two historical collision inputs plus 192 replacement branches | 194 | 4,271 | 32 | 16,366,332 |
| Deliberately selected tied minima | 10 | 278 | 33 | 3,623,634 |

All 204 returned the verified factor pair. Of 4,549 distinct requests counted within their branches, 3,584 had positive cutoff and a nonprincipal reduced form. At the two consecutive cutoffs the mixed signed residue was $0\to4$ in 192 singleton cases and $0\to8$ in 12 ties. Forty extra direct comparisons cover nonprimitive forms, axes, singular multiples and a sheared form. The larger local-solubility fixtures measure one level of the tree; they are not end-to-end timings.

From the repository root run

    uv run --no-project --with python-flint==0.9.0 python experiments/ellipse_disk_mod8/centered.py 11009
    python experiments/ellipse_disk_mod8/centered_check.py --census 33119086993
    python experiments/ellipse_disk_mod8/centered_check.py

**Next decision.** A method to combine a whole block of these inequivalent forms at one center, or to detect which half contains a global representation using less work than examining its slopes, would address the missing evaluator. The two opposite-outcome fixtures give a direct test at $B=1$. Any method using only prime-power solubility gives the same answer on every selected child. The geometric recursion is useful as a reference implementation for a proposed aggregate.

The [equal-norm bisector search](bisector.md) distinguishes these two groups with preprocessing included. It also bounds how frequently a shallow prefix on this near-square branch can contain a mixed point.
