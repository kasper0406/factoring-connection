# What the packed congruence does to a quaternion rotation

For the fixed-field packing, a quaternion rotation satisfying a fine mark splits into two binary norm equations. One has short norm $4A$. The other has norm $4AN^2$, and its four integer coordinates have exactly the same gcd with the original input $N$ as the hidden frame partner. The two auxiliary primes disappear from these norms. This explains why their factorization facilitates counting by Hecke steps while leaving a binary norm problem involving the original primes. [Fixed-field construction](fixed-field-packing.md) · [Individual Hecke counts](factored-ternary-hecke.md) · [Earlier fixed-first-vector binary norm](short-first-vector.md)

Write $A=x^2+1$, $C=tN$, $T=B^2+C^2$, $P=B^2-C^2-2xBC$ and $Q=2ABt$. Identify vectors with pure Hamilton quaternions. Let $u$ be the first vector at the mark, let $z$ be any base vector of norm $A$, and let $h$ have norm $T$. The target is $W=h z\bar h=Pu-Qv$, where the packed Gram gives $u\cdot v=-Nx$ and $|v|^2=N^2$. In particular $u\cdot W=A(B^2-C^2)$.

The identity $Wh=T h z$ and the congruences $P\equiv-T\bmod B$, $P\equiv T\bmod t$ give integral quaternions
$$
 U=(uh+hz)/B,\qquad V=(uh-hz)/t .
$$
The packing conditions make $T$ a unit at both denominators. Multiplying out gives $uU=Uz$ and $uV=-Vz$. These two rational spaces each have dimension two, since the centralizer of a nonzero pure quaternion $z$ is $\mathbb Q(z)$. The Euclidean inner product of $uh$ and $hz$ is $u\cdot W$. Their individual squared lengths are $AT$, so
$$
 |U|^2=4A,\qquad |V|^2=4AN^2,\qquad
 h=-u(BU+tV)/(2A).
$$
For a converse enumeration, take integral points on these two planes of the stated norms, reconstruct $h$, and verify integrality and the original mark. The planes depend only on the two short vectors. The informative gcd follows directly from
$$
 V\bar h=(T u-W)/t=2N(tN+xB)u+2ABv .
$$
Right multiplication by $\bar h$ is invertible modulo $N$ and $2AB$ is a unit there. Thus $\gcd(N,V_0,V_1,V_2,V_3)=\gcd(N,v_1,v_2,v_3)$. At $x=0$, after aligning the source and target axes, the two planes are Gaussian lines. The long component amounts to a Gaussian integer of norm $N^2$.

For comparison, the generic modular linearization is $c h\equiv T h z\bmod m$. When a unit intertwiner exists, its full solution module is an invertible translate by multiplication of the centralizer of $z$. For primitive $z$, the commutator with $a+\mathbf b$ is $2\mathbf b\times z$. The cross-product map has two nonzero Smith invariants, both one. Hence the modular intertwining lattice has exact index $(m/\gcd(2,m))^2$. At an odd fine test prime $\ell\nmid A$, its fixed-determinant fiber has $\ell-(-A/\ell)$ elements at level one and $\ell^{k-1}(\ell-(-A/\ell))$ modulo $\ell^k$. This is the unit norm map in the split quadratic algebra or its unramified quadratic extension. At an odd prime dividing $A$ the relative determinant must be a square modulo $\ell$, and when it is there are $2\ell^k$ lifts: the centralizer units reduce to $a+bz$ with $a\ne0$, and the norm image is exactly the units square modulo $\ell$. Thus fine targets of the same unit norm have the same local multiplicity at the varying test prime. The remaining distinction comes from the integral norm equation.

## Recovering the possible inputs from the scalar query

Suppose the native reader receives only $(M,m,c)$ and the supplied factorization of $T$ into one or two distinct Gaussian-split primes. Recover $A=M/T^2$ and test $A=x^2+1$. Cornacchia at each known prime and Gaussian multiplication give every positive ordered decomposition $T=B^2+C^2$, at most four for two primes. For each $2AB\mid m$ put
$$
 t=\gcd(C,m/(2AB)),\qquad N=C/t,\qquad \ell=m/(2ABt).
$$
These identities follow from $C=tN$, $m=2ABt\ell$ and $\ell\nmid N$. Check the packing parameters, recover the centered first-vector mark modulo $Q$, and check the fine Gram conditions at $\ell$. One pair of scalar parameters has two genuine interpretations: $M=7501^2$, $T=13\cdot577$, $m=7548$ admits $(B,C,t,N,\ell)=(74,45,3,15,17)$ and $(51,70,2,35,37)$. The second has negative $P$.

An occupied fine class has a unique interpretation. Any two candidate first vectors satisfy $P_1u_1\equiv P_2u_2\bmod 2A$, so $2A$ divides each coordinate of $u_1\times u_2$. The cross product has length at most $A$, forcing $u_2=\pm u_1$. An actual witness gives $u_i\cdot W=A(B_i^2-C_i^2)$, so the legs $(B_i,C_i)$ agree or swap. A swap in $B_1t_1\ell_1=B_2t_2\ell_2$ would give $N_2\ell_1=N_1\ell_2$. But $N_1$ is coprime to both $N_2$ and $\ell_1$, since each Gaussian decomposition is primitive and the test prime excludes $N_1$. This rules out the swap.

The reference evaluator retains every compatible interpretation. It factors each candidate $N$ as a baseline, constructs the coordinate and noncoordinate frame orbits for each eligible semiprime, and unions points after independent checks of squared norm, residue and primitivity. The originating interpretation gives completeness on promised packed queries. Independent verification gives soundness for every branch. For $x=0,t=1$ the fine Gram check excludes the extra imprimitive coarse solution. This evaluator needs up to a constant number of original-scale factorizations.

The code uses a randomized square root of $-1$ and an injectable factor routine whose default is SymPy. For a theoretical deterministic comparison, general factoring costs $N^{1/5+o(1)}$. Deterministically finding a nonresidue at a supplied prime $s$ takes $s^{1/(4\sqrt{\mathrm e})+o(1)}$ operations by the Burgess–Vinogradov bound, sufficient to run Cornacchia at a Gaussian-split prime or Tonelli–Shanks in the frame reconstruction. In the balanced auxiliary construction $s$ and every candidate $N$ are at most the original $N^{1+o(1)}$ for short $x$. Thus the resulting deterministic scalar baseline is also $N^{1/5+o(1)}$ in this regime. For arbitrary oversized conductors the scale is set by the supplied primes and the reconstructed candidates. [Harvey and Hittmeir](https://doi.org/10.1090/mcom/3708) · [Least quadratic nonresidue bound](https://arxiv.org/abs/1311.7556) · [Frame reconstruction](prescribed-ternary-embedding.md)

At original $N=15$, the cases $(x,B,T)=(0,52,29\cdot101)$ and $(1,62,13\cdot313)$ each have 72 marked primitive points. The component $V$ gives gcd $N$ on the 24 coordinate frames. On all 48 other frames it gives 3 at $x=0$ and 5 at $x=1$. [Component and scalar reader](packed_quaternion_collapse.py) · [Reproducible checks](check_packed_quaternion_collapse.py)

Sardari's strong approximation result treats a prescribed four-square problem under a factoring assumption and a conjecture about sums of two squares. The hardness construction in §2 encodes separate choices in a product of Gaussian prime norms. Applying it to a single-prime ternary conjugation orbit would require an additional reduction. [Sardari, *Complexity of Strong Approximation on the Sphere*](https://doi.org/10.1093/imrn/rnz233)
