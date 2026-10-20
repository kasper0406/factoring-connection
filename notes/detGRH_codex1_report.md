# Deterministic factoring under GRH: certificate lattices, square roots, and the missing relation theorem

Independent report, 20 October 2026.

**Conclusion — OPEN.** The proposed argument does not establish deterministic subexponential complete factorization, either for every integer or for a density-one set of integers. Its character-theoretic conclusion is valid under a strong, precisely formulated hypothesis. That hypothesis is impossible for many moduli, and the supplied NFS/descent construction does not establish it for the remaining moduli.

**PROVED.** There is a useful alternative reduction: a full-rank lattice of verified multiplicative relations for generators of \((\mathbb Z/N)^\times\) gives an exponent multiple, and under GRH that exponent multiple yields deterministic polynomial-time factoring. The relation lattice need not be the full kernel and need not be saturated. Finding such a full-rank lattice deterministically in subexponential time is the unresolved step.

**Scope and labels.** PROVED means proved below, or an explicitly identified published theorem under its stated assumptions. FALSE means an asserted implication has a counterexample. OPEN includes missing algorithmic or analytic statements; it does not mean their negations are proved. HEURISTIC means an explicitly unproved model. Literature absence is reported as a search finding, not as a theorem of impossibility.

I used the permitted background notes and the indicated FM/LD* proof sections. I did not use other researchers' detGRH notes, repository status reports, or material in `codex/`. The randomized Theorem A and its internal analytic inputs remain unrefereed premises here: this report does not certify their complete proof chain. The elementary obstructions and reductions below are independent of their validity.

## 1. Notation and the two different kinds of relations

Write

\[
N=\prod_{i=1}^{\omega}p_i^{e_i},\qquad
G=(\mathbb Z/N)^\times,
\quad T=G[2],\quad D=\{1,-1\}\subseteq T.
\]

We may first remove powers of 2, recognize primes, and recognize perfect powers deterministically. At a splitting step, assume that \(N\) is odd and \(\omega\ge2\). Repeated odd prime factors are allowed. The Chinese remainder theorem gives

\[
T\simeq\{\pm1\}^{\omega},\qquad T/D\simeq\mathbb F_2^{\omega-1}.
\]

For \(z\in T\), a nonzero class in \(T/D\) is equivalent to

\[
1<\gcd(z-1,N)<N.
\tag{1.1}
\]

Indeed, the gcd is the product of precisely those prime powers at which the CRT component of \(z\) is \(+1\).

Throughout, \(L=L_N[1/2,1]=\exp(\sqrt{\log N\log\log N})\). Thus \(L^{O(1)}\) is subexponential in the bit length of \(N\), whereas \(N^c\), for fixed \(c>0\), is exponential.

### 1.1 Certificate group versus modular presentation

**PROVED: these are different kernels.** Let \(C\) be a multiplicative certificate group with an evaluable homomorphism

\[
\operatorname{ev}_C:C\longrightarrow G.
\]

For Dixon one can use nonzero rational numbers whose numerator and denominator are coprime to \(N\). For NFS one uses a suitable subgroup of \(\mathbb Q^\times\times K^\times\), with the reduction at the chosen root of the field polynomial. Only elements at which that reduction is defined and invertible are admitted. Leading coefficients, denominator corrections, and derivative factors must be included in the actual identity; they cannot be discarded as harmless constants.

Choose certificate atoms and let

\[
\Gamma=\mathbb Z^m\xrightarrow{\ c\ }C,
\qquad \operatorname{ev}=\operatorname{ev}_C\circ c,
\qquad \Lambda=\ker(\operatorname{ev}).
\tag{1.2}
\]

This free formal presentation avoids silently quotienting out signs or other torsion. Let \(\Lambda_0\subseteq\Lambda\) be the lattice of verified modular relations actually found. Then

\[
\Gamma/\Lambda\simeq H:=\operatorname{ev}(\Gamma),
\qquad
\Gamma/\Lambda_0\twoheadrightarrow H.
\tag{1.3}
\]

By contrast,

\[
\ker c\subseteq\Lambda,
\qquad
\Sigma:=\{a\in\Gamma:c(a)\in C^2\}
\tag{1.4}
\]

are respectively exact certificate identities and certificates that are squares. Generally \(\Sigma\ne2\Gamma\): a certificate can have a square root in \(C\) that is absent from the chosen atom list.

An NFS row \(r_i^2=\operatorname{ev}_C(c_i)\) is a modular relation only after adding a formal generator for \(r_i\): it gives

\[
2a_i-b_i\in\Lambda,
\quad\operatorname{ev}(a_i)=r_i,
\quad c(b_i)=c_i.
\tag{1.5}
\]

It is not normally an exact identity in \(C\). Also, adding one new seed generator and one row does not by itself bring the modular relation lattice closer to full rank: both the number of generators and the number of relations increase.

### 1.2 The half-relation map

**PROVED.** There is a canonical isomorphism

\[
\Theta:\frac{\Lambda\cap2\Gamma}{2\Lambda}\longrightarrow H[2],
\qquad [2a]\longmapsto\operatorname{ev}(a).
\tag{1.6}
\]

It is well defined because \(2a\in\Lambda\) implies \(\operatorname{ev}(a)^2=1\). Its kernel consists precisely of \(a\in\Lambda\), and every element of \(H[2]\) has a lift \(a\). The quotient in (1.6) is therefore an \(\mathbb F_2\)-space. Composing with \(T\to T/D\) gives the factor-detecting half-relation map.

For a target \(r\), choose a lift \(a_r\), and suppose a descent gives \(b_r\) with

\[
\operatorname{ev}(b_r)=r^2.
\]

If \(b_r=2h_r\) in this presentation, the certificate square root gives

\[
y_r=\operatorname{ev}(h_r),\qquad
\chi_r=\frac{y_r}{r}
=\operatorname{ev}(h_r-a_r)
=\Theta([b_r-2a_r]).
\tag{1.7}
\]

If the root lies outside the original presentation, adjoining that computed root gives the same formula in an enlarged \(\Gamma\). Formula (1.7), rather than a Dirichlet character assertion, is what is always justified.

If a different descent has \(b'_r=b_r+2z\), then \(2z\in\Lambda\), and its root differs by \(\operatorname{ev}(z)\in T\). A nondiagonal difference already factors \(N\). Thus independence of the root from the descent certificate is a substantive condition, not a consequence of Gaussian elimination.

### 1.3 Where multiplicativity actually holds

**PROVED.** Multiplication holds for evaluated certificates and for formal products of chosen roots. If square roots of a free list of square certificates are fixed, multiplying their roots defines a homomorphism on the free exponent group. This gives a homomorphism on modular residue classes only if every relation between two representations is respected by the chosen roots.

Concretely, to define a root choice on \(\Gamma/\Lambda\), every \(\lambda\in\Lambda\), not merely every found \(\lambda\in\Lambda_0\), must have trivial root defect. Consistency on \(\Lambda_0\) defines at most an object on the larger presented group \(\Gamma/\Lambda_0\). Identifying that object with a function of the residue modulo \(N\) is unjustified.

In NFS, usually only products of rows are squares in \(C\). There is then no separately defined \(y(r_i^2)\) for each row. The algorithm produces

\[
Y_S^2=\prod_{i\in S}r_i^2,
\]

after choosing a square dependency \(S\). Treating this as evaluation of an already constructed root homomorphism on every individual target reverses the construction.

**PROVED: choices are deterministic if ordered, but not intrinsic.** A fixed row order, pivot rule, and setting of free variables give a reproducible dependency. The dependency space is canonical for the specified matrix; a selected nonzero vector, a complement, and a solution to an inhomogeneous descent system are not canonical without these extra choices. Binary linear algebra enforces parity, not compatibility with every modular relation or with all carries from multiplying representations.

For a fixed true NFS square dependency, the rational and number-field square roots are each unique up to the relevant field's square roots of unity. In the odd-degree setting those are \(\pm1\), so their quotient class in \(T/D\) is well defined. Different dependencies can give different classes and hence different factors.

## 2. What a multiplicative modular square-root selector would imply

The correct domain of such a selector is the square subgroup

\[
Q=G^2.
\]

Asking for a square root of every element of \(G\) is impossible because \(G\) contains nonsquares. The intended hypothesis is a homomorphism

\[
y:Q\longrightarrow G,
\qquad y(q)^2=q.
\tag{2.1}
\]

### 2.1 Exact existence criterion

**PROVED: (2.1) exists if and only if every \(p_i\equiv3\pmod4\).**

Suppose some \(p_i\equiv1\pmod4\). Its local unit group contains an element of order 4. By CRT there is \(a\in G\) of order 4. Then \(q=a^2\in Q\) has order 2. A homomorphism has \(y(q)^2=y(q^2)=y(1)=1\), contradicting \(y(q)^2=q\ne1\).

Conversely, if all \(p_i\equiv3\pmod4\), then

\[
G\simeq T\times Q,
\qquad |Q|\text{ odd}.
\]

Squaring is an automorphism of \(Q\). Its inverse, included into \(G\), is a selector satisfying (2.1). It is the unique homomorphic selector: two selectors differ by a homomorphism \(Q\to T\), and such a homomorphism is trivial because \(|Q|\) is odd.

An equivalent proof exposes the character obstruction. If

\[
\chi(r)=y(r^2)r^{-1},
\tag{2.2}
\]

is a homomorphism, then \(\chi:G\to T\) is a retraction: for \(\varepsilon\in T\), \(\chi(\varepsilon)=\varepsilon\). Every homomorphism to \(T\) kills \(G^2\), so \(T\cap G^2\) must be trivial. A local prime \(1\pmod4\) produces a nontrivial element of this intersection.

For example, modulo \(65\), \(8^2=-1\) and \(-1\) has order 2. Multiplicativity would simultaneously require \(y(-1)^2=1\) and \(y(-1)^2=-1\).

**PROVED: passing to the diagonal quotient does not remove this obstruction when \(\omega\ge2\).** Suppose only that the root ratio \(\bar\chi(r)=[y(r^2)/r]\in T/D\) is a homomorphism. Its value at the identity forces \(y(1)\in D\), and for \(\varepsilon\in T\) it must satisfy \(\bar\chi(\varepsilon)=[\varepsilon]\). If some \(p_i\equiv1\pmod4\), take \(\varepsilon\) to be \(-1\) at that prime-power component and \(+1\) elsewhere. This is a square in \(G\), yet is nondiagonal. Every homomorphism to an exponent-2 group kills squares, contradicting \(\bar\chi(\varepsilon)=[\varepsilon]\ne0\). Thus a global character argument cannot be restored merely by allowing sign ambiguity in the selected roots.

**PROVED: in the possible case the components are prescribed, not arbitrary.** Put

\[
m_i=\frac{p_i^{e_i-1}(p_i-1)}2,
\]

which is odd. The inverse of squaring on the local square group sends \(r^2\) to \(r^{m_i+1}\). Hence

\[
\chi(r)\bmod p_i^{e_i}=r^{m_i}
=\left(\frac r{p_i}\right)\in\{\pm1\}.
\tag{2.3}
\]

The last equality follows because the left side is a root of 1 and has the indicated Legendre symbol value after reduction modulo \(p_i\). Thus \(\chi\) is the vector of the independent local Legendre characters. In particular its image is all of \(T\), so for \(\omega\ge2\) it cannot land in \(D\).

Every real character of \(G\) is a product of the local Legendre characters: each odd prime-power unit group is cyclic and has exactly one nontrivial character to \(\{\pm1\}\). This classification applies only after a homomorphism on \(G\) has been established.

Including the seed \(-1\) imposes \(y(1)=1\) and \(\chi(-1)=-1\), the diagonal element. It does not itself yield a factor. In particular, treating \(-1\) as an extra freely adjustable sign is invalid.

### 2.2 The valid GRH implication

**PROVED, conditional on GRH.** Suppose (2.1) holds and is evaluable on squares of the bases being tested. If \(\bar\chi:G\to T/D\) is nonzero, some unit integer

\[
2\le r<2(\log N)^2
\tag{2.4}
\]

has \(\bar\chi(r)\ne0\), unless the preliminary gcd tests have already found a factor.

Choose a linear functional on \(T/D\) nonzero on the image of \(\bar\chi\). Its composition with \(\bar\chi\) is a nonprincipal real Dirichlet character modulo \(N\). The Ankeny–Bach least-character-witness theorem under GRH gives (2.4). Testing all these bases requires no knowledge of the character or of the factorization. A computed root detects the desired event through (1.1). See [Bach, 1990](https://doi.org/10.1090/S0025-5718-1990-1023756-8).

**FALSE as an unconditional implication from multiplicativity alone.** A homomorphism whose image is diagonal need not have a fruitful generator. The trivial homomorphism is an immediate counterexample if no root-selector/fibre condition is imposed. The root-selector condition changes the situation: it makes a global multiplicative selector impossible at primes \(1\pmod4\), and forces a nontrivial quotient character when all primes are \(3\pmod4\) and \(\omega\ge2\).

**PROVED, conditional on the full Dirichlet GRH.** The units among the positive integers below the Bach bound generate \(G\); otherwise a nonprincipal character of the quotient by their generated subgroup contradicts the corresponding small-character-witness theorem. Consequently, saying that (2.1) holds on the square subgroup of the group generated by *all* these bases does not avoid the existence obstruction. Under GRH that group is \(G\).

For a proper subgroup \(H\), a homomorphic inverse of squaring \(H\to H^2\) exists exactly when its Sylow 2-subgroup has exponent at most 2. This follows by applying the order-4 argument and splitting a finite abelian group into its odd and 2-primary parts. Multiplicativity merely on a free monoid of short integer words says much less and has no GRH witness consequence.

### 2.3 The proposed small seeds fail the stated Gauss/shape construction

**PROVED.** For every \(r\ge2\) satisfying

\[
r^4+1<N,
\tag{2.5}
\]

the shortest nonzero vector in

\[
\mathcal L_{r^2}=\{(u,v)\in\mathbb Z^2:u-r^2v\in N\mathbb Z\}
\]

is \(\pm(r^2,1)\).

For a vector with \(u-r^2v=kN\ne0\), Cauchy–Schwarz gives

\[
\|(u,v)\|\ge\frac N{\sqrt{1+r^4}}>\sqrt{1+r^4}.
\]

For \(k=0\), the vector is an integer multiple of \((r^2,1)\). This proves the claim.

Condition (2.5) holds eventually for every \(r\le2\log^2N\), and indeed for an initial segment bounded by \(L^{o(1)}\). The Theorem A shape test requires \(|u|\) and \(|v|\) to be squarefree products of primes in \(J\), where the primes of \(J\) have size \(N^{2\nu+o(1)}\). But \(u=r^2\) is not squarefree and its prime factors are much smaller than \(J\). Thus these seeds are rejected. Replacing the factoring subroutine by recursion leaves this rejection unchanged.

If instead one accepts their obvious factor-base certificates, \(r^2<N\) is an ordinary integer square, its rational square root is \(r\), and the resulting congruence is tautological. In (1.7) this is \(b_r=2a_r\), so \(\chi_r=1\).

**FALSE: formal multiplicativity supplies the needed modular selector.** The positive rational square-root function is multiplicative on positive rational squares and returns \(r\) on these targets. It does not descend to the finite modular square subgroup. For illustration, modulo \(77\) the unique homomorphic selector on \(Q\) sends \(4\) to \(9\), whereas the rational square root sends the integer \(4\) to \(2\). The roots differ fruitfully: \(\gcd(9-2,77)=7\).

**PROVED: actual algorithms can satisfy the hypothesis in the possible case if the exponent information is available.** Given any multiple \(M=2^st\) of \(\lambda(N)\), with \(t\) odd, for an all-prime-\(3\pmod4\) modulus the odd number \(t\) annihilates \(Q\). Thus \(y(q)=q^{(t+1)/2}\) is a computable homomorphic selector. Factoring \(N\) first also supplies this information. This proves satisfiability by an actual algorithm, without supplying a fast way to obtain the information.

**OPEN.** Producing a homomorphic selector in subexponential time for all-prime-\(3\pmod4\) moduli would suffice under GRH. Its abstract existence is proved above; constructing it without already obtaining the required arithmetic information is not supplied by NFS descent. For other moduli, exact multiplicativity of this kind is algebraically impossible.

## 3. Structural fruitlessness and the functional \(\Phi\)

### 3.1 A general proof of the homomorphism statement

**PROVED.** Suppose rows \((x_i,c_i)\) satisfy

\[
x_i^2=\operatorname{ev}_C(c_i),\qquad x_i\in G.
\]

Assume \(\operatorname{ev}_C(C[2])\subseteq D\); this holds for the rational/odd-degree-field certificate group. Let

\[
K_0=\left\{s\in\mathbb F_2^t:\prod_i c_i^{s_i}\in C^2\right\}.
\]

For \(s\in K_0\), choose a computed certificate root \(g_s\), and put

\[
\Phi(s)=\left[\frac{\prod_i x_i^{s_i}}
{\operatorname{ev}_C(g_s)}\right]\in T/D.
\tag{3.1}
\]

Different roots differ by \(C[2]\), so this is well defined. For two subsets \(s,t\) with intersection \(I\), a root for their symmetric difference is

\[
g_sg_t/\prod_{i\in I}c_i.
\]

The corresponding quotient of the products of \(x_i\)'s has the correction \(\prod_{i\in I}x_i^2\). These corrections cancel under evaluation. Hence

\[
\Phi(s+t)=\Phi(s)+\Phi(t).
\tag{3.2}
\]

This independently verifies the substantive content of Lemma H with its necessary evaluability and torsion assumptions. If \(\rho=\operatorname{rank}\Phi\), the proportion of all vectors of \(K_0\), including zero, giving a factor is exactly

\[
1-2^{-\rho}.
\tag{3.3}
\]

### 3.2 Nonzero \(\Phi\): no seed and no GRH are needed

**PROVED.** Given an explicit basis \(s_1,\ldots,s_k\) of the actual square-dependency space and evaluable square roots, if \(\Phi\ne0\), some basis vector is fruitful. Otherwise linearity makes \(\Phi\) zero everywhere. Test the basis by exact square roots and gcds. This costs at most \(k\) root computations, not \(2^k\) trials.

This derandomizes *kernel sampling*, with a possible increase in the leading \(L\)-constant. It does not prove \(\Phi\ne0\), and a kernel obtained from incomplete parity/character information need not yet be \(K_0\).

**FALSE: \(\Phi\ne0\) alone implies a small fruitful individual seed.** Its domain is a space of dependencies, not \(G\). A target may not even have an individual root certificate. Even where targets do, a selected branch can be fruitless on every small target while another dependency is fruitful. For example, modulo \(77\), the certificate \(4\) can be paired with stored roots \(2\) or \(9\); both are valid rows, and the second gives a factor. A separate target rule that returns the rational root on every small integer square still gives \(\chi_r=1\) throughout that small range. No assertion about small arguments of a Dirichlet character follows from the nonzero dependency functional.

### 3.3 An exact structural criterion

**PROVED.** Fix the entire allowed relation family and its actual dependency space \(K_0\). Structural fruitlessness means exactly

\[
\Phi|_{K_0}=0.
\tag{3.4}
\]

Sampling fruitlessness means that the vectors tested so far lie in \(\ker\Phi\), although (3.4) is false. Structural fruitlessness is relative to a specified relation family; enlarging the family can change it.

For an exact square-class matrix \(E:\mathbb F_2^t\to V\), \(K_0=\ker E\). Extend each coordinate of \(\Phi\) to a linear row \(\beta_j\) on \(\mathbb F_2^t\). Then

\[
\Phi=0\quad\Longleftrightarrow\quad
\beta_j\in\operatorname{rowspace}(E)\quad\text{for every }j.
\tag{3.5}
\]

This is the identity \((\ker E)^\perp=\operatorname{rowspace}(E)\). It is an exact characterization, not an equidistribution theorem. Dimension or the number of dependencies alone does not rule it out: choose any \(\beta_j\) in that rowspace.

There is no general prohibition against (3.4). The family \(c_i=r_i^2\) in the rational certificate group, with stored roots \(x_i=r_i\), has \(\Phi=0\) identically for every \(N\). It is a legitimate, completely fruitless family of certificates.

On the other hand, a single-valued selector \(y:Q\to G\) cannot have \(y(r^2)/r\in D\) for *every* \(r\in G\), even without multiplicativity. For fixed \(R\), the quotients as \(r\) ranges over its root fibre exhaust \(T\). This does not contradict fruitlessness of a deterministic relation family that sees only one chosen root class from each encountered fibre.

### 3.4 The second-order formula does not manufacture a character of the seed

**PROVED, in the semiprime NFS setting.** Let \(N=pq\). Choose characters \(\eta_p,\eta_q\) of orders \(2^{v_2(p-1)}\), \(2^{v_2(q-1)}\), respectively. Then \(\eta_p(-1)=\eta_q(-1)=-1\). Put \(c=\max(v_2(p-1),v_2(q-1))\), \(\zeta=e^{2\pi i/2^c}\), and \(\eta=\eta_p\eta_q\), interpreted also on field elements by their two reductions. In particular \(\eta(-1)=1\).

Use exact rational and field multiplicative coordinates for a finite list of ordinary NFS rows \(q_i\equiv\alpha_i\pmod N\). Write

\[
q_i=\pm\prod_a u_a^{e_{ia}},\qquad
\alpha_i=\pm\prod_b v_b^{f_{ib}},
\]

and choose lifts \(\eta(u_a)=\zeta^{A_a}\), \(\eta(v_b)=\zeta^{B_b}\). Then

\[
t_i=\sum_a A_ae_{ia}-\sum_bB_bf_{ib}\equiv0\pmod{2^c}.
\]

On a square dependency all summed coordinates are even, so for its two roots \(U,V\),

\[
\eta(U)/\eta(V)
=\zeta^{\sum_{i\in S}t_i/2}
=(-1)^{\sum_{i\in S}t_i/2^c}.
\tag{3.6}
\]

For a CRT sign pair \((\varepsilon_p,\varepsilon_q)\), its \(\eta\)-value is \(\varepsilon_p\varepsilon_q\), which detects the nonzero element of \(T/D\). Thus the per-row values \(\beta_i=t_i/2^c\bmod2\) give \(\Phi\) on dependencies. This verifies the Bockstein mechanism in its stated setting. It requires higher 2-power character data and exact field coordinates; it does not assert that the resulting per-row function is a quadratic Dirichlet character of a small seed.

### 3.5 The B1 obstruction and what GRH does remove

**PROVED, with the restricted scope.** Consider a Blum semiprime and a degree-1 two-form construction with

\[
L_1(P)\equiv\lambda L_2(P)\pmod N,
\]

and an even-size dependency \(|S|=2h\). If \(U^2=\prod_SL_1(P)\) and \(V^2=\prod_SL_2(P)\), the resulting root of 1 is

\[
z=U/(\lambda^hV).
\]

Suppose every prime in the factor base has Jacobi symbol \(+1\). Signs also have Jacobi symbol \(+1\) for a Blum semiprime. Therefore \(J_N(U)=J_N(V)=1\). Each admitted row forces \(J_N(\lambda)=1\), so \(J_N(z)=1\). On \(T\) for a Blum semiprime, the two nondiagonal sign pairs have Jacobi symbol \(-1\). Hence every such dependency is fruitless.

Under GRH, a nonprincipal Jacobi character has a nonresidue below \(2\log^2N\). Thus a sufficiently large factor base cannot satisfy this particular absence-of-nonresidues condition. **The converse is FALSE:** putting a nonresidue in the factor base does not force any admitted relation to use it with the exponent pattern needed for (3.6), or force \(\beta\) outside the rowspace. An unused nonresidue added to a tautological relation family changes no dependency's fruitfulness.

**OPEN.** GRH alone has not here been shown to exclude (3.4) for deterministic QS or NFS relation families. B1 is a necessary obstruction for the particular two-form mechanism, not a general explanation of QS fruitfulness or a sufficient condition after it disappears. For example, \(7^2-33=16\) gives \(\gcd(7-4,33)=3\), although \((2/33)=1\).

### 3.6 Why the hidden-root proof cannot be reused on small known roots

**PROVED.** If \(r\) is uniform in \(G\), conditional on \(R=r^2\), it is uniform on a coset of \(T\). If the full transcript, including acceptance and dependency choice, is invariant under independently multiplying each stored descent root by \(T\), a dependency containing at least one such root has uniformly distributed quotient in \(T\). Exactly two of its \(2^\omega\) elements are diagonal, giving success probability \(1-2^{1-\omega}\).

This proves the stated DD conclusion under its transcript-invariance assumptions. A deterministic small integer is not uniform on its root fibre. Fixing a reproducible rule for constructing the other root does not restore that conditional distribution. Earlier factor tests or root-dependent selection can also destroy it, so the argument is per fresh invariant trial, not an automatic statement after arbitrary adaptive failures.

**PROVED: fixed-test averaging is insufficient in principle.** For a concrete model, index \(2^t\) inputs by bit strings \(v\in\mathbb F_2^t\), set the relation labels to the basis vectors \(e_j\), and set \(\beta_v(j)=v_j\). The predictor \(w=v\) fits every label exactly. For each fixed \(w\), however, the disagreement fraction for uniform \(v\) has mean \(1/2\) and variance \(1/(4t)\); Chebyshev gives an exceptional fraction at most \(1/(4t\varepsilon^2)\) for disagreement outside \(1/2\pm\varepsilon\). Thus every fixed test is balanced for almost all inputs, while an input-dependent test annihilates every input. This is the logical obstruction illustrated by the internal Theorem B. Balance against fixed tests cannot exclude the input-dependent rowspace in (3.5).

## 4. Integer relations give a route without NFS square-dependency fruitfulness

### 4.1 Full rank is enough; full kernel is unnecessary

**PROVED.** Let \(g_1,\ldots,g_m\) generate \(G\). Suppose an integer matrix \(A\in\mathbb Z^{m\times m}\) has nonzero determinant and every row is verified to satisfy

\[
\prod_{j=1}^m g_j^{A_{ij}}=1\pmod N.
\tag{4.1}
\]

Then

\[
M=|\det A|
\]

is a multiple of \(|G|=\varphi(N)\), and therefore of the Carmichael exponent \(\lambda(N)\).

**Proof.** The row lattice \(\Lambda_0\) is full rank, so \(|\mathbb Z^m/\Lambda_0|=|\det A|\). Evaluation induces a surjection of this finite group onto \(G\). Its order is therefore divisible by \(|G|\). Equivalently, the adjugate shows that \(Mg_j=0\) in additive group notation for every generator, which already proves the exponent assertion. No equality \(\Lambda_0=\Lambda\), saturation, or Smith normal form is required. ∎

If \(g_j\) generate only \(H\), the same proof gives only an annihilator of \(H\). Under Dirichlet GRH one can ensure \(H=G\) by including all unit bases up to the Bach bound. A gcd with a nonunit base gives a factor instead.

If \(|A_{ij}|\le2^b\), Hadamard's inequality gives

\[
\log_2|\det A|\le mb+\frac m2\log_2m.
\tag{4.2}
\]

Determinants, relation verification, and subsequent modular exponentiation take time polynomial in the encoded matrix size and \(\log N\). Thus \(L^{O(1)}\)-time construction of an \(L^{O(1)}\)-size matrix suffices for subexponential factoring.

### 4.2 Factoring from an exponent multiple: explicit proof

**PROVED, conditional only on quadratic Dirichlet GRH.** Given \(N\) as above and a positive multiple \(M\) of \(\lambda(N)\), one can split \(N\) in time polynomial in \(\log N+\log M\).

Write \(M=2^st\), \(t\) odd. For each integer \(2\le a<2\log^2N\), first compute \(\gcd(a,N)\), and then inspect the chain

\[
a^t,a^{2t},\ldots,a^{2^st}=1\pmod N.
\tag{4.3}
\]

Take gcds with 1 subtracted from the chain entries. Any proper gcd splits \(N\).

To prove that some chain splits, put \(s_i=v_2(p_i-1)\), and \(S=\max_i s_i\). We have \(s\ge S\), and \(t\) kills every local odd-order component.

If not all \(s_i\) are equal, choose a prime with \(s_i=S\) and another with smaller \(s_j\). The nonprincipal character \((\cdot/p_i)\) has a small nonresidue \(a\) under GRH, or a small nonunit already factors \(N\). In the chain, at exponent \(t2^{S-1}\), the \(p_i^{e_i}\)-component is \(-1\), whereas every component with \(s_j<S\) is \(+1\). Every other component is also \(\pm1\). This is a nondiagonal root of 1, so its gcd is proper.

If all \(s_i=S\), choose two distinct prime divisors \(p_i,p_j\). The character

\[
\left(\frac\cdot{p_i}\right)
\left(\frac\cdot{p_j}\right)
\]

is nonprincipal. A small base where it is \(-1\) has different signs at these two components of \(a^{t2^{S-1}}\). At every component this entry is \(\pm1\), so (1.1) applies. This proves the algorithm. ∎

The prime-power case is handled by the preliminary perfect-power test; recursive splitting handles the full factorization. This is the group-exponent mechanism behind [Miller, 1976](https://www.cs.cmu.edu/~glmiller/Publications/Papers/Mi76.pdf), with the small-base bound supplied by [Bach, 1990](https://doi.org/10.1090/S0025-5718-1990-1023756-8).

This route does not use an \(\mathbb F_2\) NFS square dependency, its functional \(\Phi\), or an exact square root in the number field. The final Miller step can still expose a nontrivial modular square root of 1; the distinction is the source of the exponent, not an assertion that every operation avoids the number 2.

### 4.3 What exact NFS dependencies provide, and what they do not

**PROVED as a conditional algebraic statement.** Exact cancellation of rational exponents and all algebraic ideal valuations leaves units. Exact cancellation of the units modulo roots of unity then gives genuine multiplicative relations modulo \(N\), up to explicitly removable torsion. In odd-degree fields that torsion is \(\pm1\), which can be removed by doubling the relation.

Archimedean logarithms can certify this second cancellation: for fixed field degree, a non-torsion unit has a positive degree-dependent lower bound on the size of its logarithm vector; sufficiently accurate certified approximations separate zero from nonzero unit logs. Integer linear algebra and LLL can therefore find and certify suitable exact relations when enough independent unit relations are supplied. This explains the role of the proposed XK mechanism. Its modular-prime rank search can be replaced by deterministic rational/integer linear algebra at a polynomial cost in the explicit input size.

Here is the relevant fixed-degree derivation. Suppose \(s>r_U\) independent integer dependencies leave units \(u_1,\ldots,u_s\), where \(s\le d\) and \(r_U\le d-1\). Their vectors \(v_j=\operatorname{Log}(u_j)\) lie in the unit logarithm lattice. A degree-dependent Mahler-measure lower bound, together with Kronecker's theorem, gives a computable \(\delta_d>0\) such that every nonzero vector of that lattice has norm at least \(\delta_d\). For fixed \(d\), this also follows from finiteness of algebraic integers of bounded degree and bounded Mahler measure, with an effective enumeration near 1.

Take an LLL-reduced basis of the unit lattice for the proof, without requiring it as algorithmic input. Its Gram–Schmidt lengths are bounded below by a fixed multiple of \(\delta_d\). Back substitution bounds the integer coordinates of every \(v_j\) by \(O_d(1+V/\delta_d)\), where \(V=\max_j\|v_j\|\). A nonzero vector of minors of a maximal-rank coordinate submatrix therefore supplies an integer relation \(n\ne0\) with

\[
\sum_jn_jv_j=0,\qquad
\|n\|\le R:=O_d(1+V/\delta_d)^d.
\]

If all vectors are zero, any coordinate vector suffices. To find a relation, apply LLL to the integer vectors
\((e_j,\lfloor2^\kappa v_j\rceil)\). The true relation gives a lattice vector of norm \(O_d(R)\), so LLL returns one of norm at most \(R'=O_d(R)\). If its coefficient vector were not a relation, its last coordinates would have norm at least \(2^\kappa\delta_d-O_d(R')\). Choose \(2^\kappa\delta_d\) larger than a suitable fixed multiple of \(R'\); this contradicts the LLL length bound. Certified approximations at \(\kappa+O_d(\log R)\) bits suffice. The output unit has zero logarithm vector and hence is a root of unity. Its modular image is \(\pm1\) in odd degree; doubling gives an exact modular relation. When the starting integer coefficients have polynomially bounded bit length in the explicit input, \(\log V\), \(\log R\), and the required precision do too. This justifies deterministic polynomial postprocessing at fixed degree, while saying nothing about full rank of its projected relations.

**FALSE: one primitive exact dependency supplies a full-rank modular presentation.** An exact dependency is one vector in the kernel. Even many more rows than parity columns need not yield independent relations on a fixed generating set of \(G\). For example, relations only on \(g_1\) produce a rank-one sublattice in \(\mathbb Z^m\), regardless of how many copies are collected. Saturation of the exact-dependency kernel does not imply that its image in seed exponent coordinates has full rank.

**OPEN.** To use §4.1, deterministic NFS/descent would have to produce enough verified relations whose vectors on a fixed small generating set span \(\mathbb Q^m\). The supplied relation-count theorems and the existence of primitive exact dependencies do not imply this spanning statement. Random hidden exponent labels can supply useful distributions in discrete-log algorithms; deterministic labels require a new argument.

## 5. Can field averages be turned into a short explicit list?

### 5.1 What Paley–Zygmund actually says

Let \(\mathcal A_N\) be the proposed field family, \(F=|\mathcal A_N|=N^{\Theta(1)}\), and let \(Z(A)\ge0\) be an observable count. First and second moments give

\[
\Pr[Z>0]\ge\frac{(\mathbb EZ)^2}{\mathbb EZ^2}.
\tag{5.1}
\]

For the actual two-stage formulation, FM/LD* first average over both \(A\) and \(\ell\). If \(g_A\in[0,1]\) is the fraction of descendable primes and \(\mathbb E g_A\ge\gamma\), then

\[
\Pr[g_A\ge\gamma/2]\ge\gamma/2.
\tag{5.2}
\]

Indeed \(\mathbb Eg_A\le\gamma/2+\Pr[g_A\ge\gamma/2]\). A bound \(\gamma\ge L^{-o(1)}\) thus gives many globally useful fields, conditional on the asserted moment inputs.

**PROVED: global moments do not imply a prescribed short list hits them.** Let \(T_N\subset\mathcal A_N\) have size \(L^{O(1)}\). Choose a set of good fields of density \(\gamma\), disjoint from \(T_N\), and put \(Z=a\) on that set and zero elsewhere. Then

\[
\mathbb EZ=a\gamma,\qquad
\mathbb EZ^2=a^2\gamma,
\]

and (5.1) is exact. This example can even have fixed positive \(\gamma\). Every global moment is compatible with complete failure on the prescribed list.

With no information about the position of a good set of size \(\gamma F\), a list that must intersect every such set needs more than \(F-\lceil\gamma F\rceil\) distinct elements. Random sampling instead needs \(O(\gamma^{-1}\log(1/\varepsilon))\) samples for failure probability \(\varepsilon\). These are different assertions.

Counts observed during a successful computation certify success. A lower bound on \(g_A\), however, is not automatically cheaply certifiable: its definition ranges over \(N^{2\nu+o(1)}\) target primes. Positive verification does not itself give an efficient pessimistic estimator or a deterministic hitting set.

### 5.2 Short coefficient boxes: the existing proof's scale

**PROVED about the hypotheses of the stated method; OPEN beyond them.** In the supplied fixed-degree design,

\[
M=N^\mu,\quad H=N^h,\quad
V=N^{\nu+o(1)},\quad
t=|G(P)|\asymp MV=N^{\mu+\nu+o(1)},
\]

and \(w\) ranges over primes at scale

\[
W=N^{\mu+\nu+\eta_0/2}.
\]

The displayed proof uses:

* rational fibres of length \(M/V\), with a power-of-\(N\) margin over the coprimality errors;
* a coefficient window satisfying the first-moment kernel/window conditions;
* for the second moment, \(H\gtrsim V^2=N^{2\nu+o(1)}\), to count translates of a kernel lattice generated by products of two linear forms;
* \(\log(W/t)\gg\log N\) in the multiplier residue-spreading bound.

The last bound follows directly by dividing Brun–Titchmarsh's estimate for primes in a class modulo \(t\) by the number of primes in \((W,2W]\):

\[
\Pr[w\equiv a\pmod t]
\ll\frac{\log W}{\varphi(t)\log(W/t)}.
\tag{5.3}
\]

For the coefficient pair map, the kernel contains forms divisible by
\((bx-ay)(b'x-a'y)\); its natural basis has coefficients of size \(V^2\). The box-volume comparison used in LD* therefore has an error of the same scale as the side length if that side is shortened far below \(V^2\). This is a failure of this particular proof estimate, not a universal lower bound on all possible short-box arguments.

The family size in this design is

\[
|\mathcal A_N|\asymp M^2H^d\frac W{\log W}
=N^{3\mu+dh+\nu+\eta_0/2+o(1)}.
\tag{5.4}
\]

Thus the supplied positive margins give power-size families, not \(L^{O(1)}\)-size boxes. Even the \(c\)-box alone has at least \(N^{2d\nu+o(1)}\) points at the second-moment scale. Translating or shortening the boxes, freezing \(g\) or \(w\), or taking a diagonal list requires new estimates; the average over the omitted parameters cannot simply be retained.

### 5.3 Initial segments and second moments

**PROVED.** For a particular explicit list \(A_1,\ldots,A_Q\), valid first/second moments over that list imply existence by exactly the same inequality (5.1). There is no probabilistic-versus-deterministic objection to this use of the second-moment method. The missing theorem is the list-specific first moment and its error control.

There is no numerical “minimal \(Q\)” determined by global FM/LD*. If a proposed local estimate is

\[
\sum_{j\le Q}Z(A_j)\ge Q\mu-\mathcal E(Q),
\]

its existence threshold is \(Q\mu>\mathcal E(Q)\). A global error bound provides no suitable \(\mathcal E(Q)\) for an arbitrary short initial segment. Calling the segment explicit does not prove this inequality.

### 5.4 Individual progressions and character sums

**PROVED as limitations of the indicated estimates.** Smooth numbers in a single prescribed progression are a different problem from averages over moduli and residue classes. The identity

\[
\Psi(X,B;t,a)=\frac1{\varphi(t)}
\sum_{\chi\bmod t}\overline{\chi(a)}
\sum_{\substack{n\le X\\P^+(n)\le B}}\chi(n)
\tag{5.5}
\]

shows the missing uniform smooth-character-sum input. The ordinary GRH prime-counting theorem in progressions cannot be substituted for the inner smooth sum.

For comparison, even when a square-root-size error \(E\asymp X^{1/2}\operatorname{polylog}(Xt)\) is available for a counting problem with density \(\rho\), a main term \(X\rho/\varphi(t)\) dominates only at

\[
X\gg (\varphi(t)/\rho)^2\operatorname{polylog}(Xt).
\tag{5.6}
\]

This is a conditional error-balancing calculation, **not** a claimed general GRH theorem for smooth numbers. For primes \(w\) in an individual progression, the standard GRH error similarly gives a sufficient scale \(W\gg\varphi(t)^2\operatorname{polylog}(Wt)\), still a power of \(N\) at the moduli in this construction.

Ordinary completion of a nontrivial character sum over a prime modulus \(q\) has an error of order \(q^{1/2}\) up to degree and logarithmic factors. To exploit it on an interval with expected good density \(\gamma\), the main term must exceed the aggregate error, including the Fourier mass of the good-field predicate. That predicate has not been shown to be a single character, a bounded-complexity collection of character tests, or to have bounded Fourier mass. Consequently no \(L^{O(1)}\) list bound follows from that approach as stated.

The Ankeny exception is specific and powerful: a nontrivial character has a small witness. A smoothness/large-relation-count event need not be the complement of a proper character subgroup. Replacing one by the other is precisely the unproved step.

### 5.5 GRH expanders

**PROVED: mixing a random walk does not directly enumerate a short hitting list.** Suppose, even granting a suitable transitive graph on the \(F\) fields, its normalized second eigenvalue is \(\lambda<1\), and its good set has density \(\gamma\). From a fixed start, the standard spectral estimate gives

\[
|P^t1_{\rm good}(x)-\gamma|
\le\lambda^t\sqrt{F\gamma}.
\tag{5.7}
\]

A sufficient mixing time for positive hitting probability satisfies

\[
t>\frac{\log(F/\gamma)}{2\log(1/\lambda)}.
\]

Enumerating all words of length \(t\) in a degree-\(D_0\) graph costs \(D_0^t\); for a fixed-degree fixed-gap expander this is a power of \(F\). Independently of spectral estimates, the union of endpoints of a list that must hit *every* density-\(\gamma\) subset must contain more than \((1-\gamma)F\) vertices. A deterministic chosen path can remain in a bad region.

GRH expansion of Cayley graphs in residue or class groups is useful for randomized sampling. Here there is additionally no established identification of the field parameter family with a transitive group action preserving the needed count statistics. A special low-complexity hitting theorem could change the conclusion; the generic expander method does not supply one.

### 5.6 What “existence of a short list” means

**PROVED.** For each individual \(N\), some good field exists if the asserted positive moment is true. A nonuniform list of length one could name it. The algorithm needs a *uniform, efficiently constructible* list. Even a length-dependent advice list working for all inputs of a given bit length is a nonuniform derandomization assertion, not a uniform GRH algorithm.

Under the Impagliazzo–Wigderson assumption that some language in \(E=\mathrm{DTIME}(2^{O(n)})\) requires circuits of size at least \(2^{\varepsilon n}\) for every sufficiently large \(n\), for some fixed \(\varepsilon>0\), a polynomial-size hitting set for the bounded successful computation follows by a pseudorandom generator. To see the reduction, truncate a Las Vegas factoring algorithm at twice a uniform expected-time bound \(T\), verify outputs, view success as a circuit of size \(T^{O(1)}\), and enumerate all generator seeds of length \(O(\log T)\). The total time is \(T^{O(1)}\). This validates the logic of the internal DR route, conditional on that additional complexity assumption; GRH is not that assumption. See [Impagliazzo–Wigderson, STOC 1997, 220–229](https://doi.org/10.1145/258533.258590).

**OPEN.** None of the supplied FM/LD*, individual-progression, generic character-sum, or expander arguments proves the desired uniformly explicit \(L^{O(1)}\) field list. The scale requirements above are requirements of specified methods or abstract information bounds, not impossibility theorems for arithmetic derandomization.

## 6. Inventory of randomness and deterministic costs

The table concerns the fixed-degree Theorem A architecture. “Removable” preserves a subexponential bound where stated; it need not preserve its leading constant. All gcd and output verification steps are deterministic.

| Source or apparent source | Status and deterministic replacement | Cost or unresolved point |
|---|---|---|
| \((g_0,g_1)\), coefficient vector \(c\), and selection of a field | **OPEN** at \(L^{O(1)}\) cost under GRH. Exhaust the stated finite family for an exponential fallback. | Existing family has size \(N^{\Theta(1)}\); global moments do not produce a short explicit hitting list. |
| Prime multiplier \(w\) | **OPEN** to fix it or restrict it to a short list while retaining LD*. Finding an individual prime in the prescribed dyadic interval is **PROVED** deterministically under GRH by scanning a sufficiently long GRH-guaranteed prime-containing interval. | Scanning can cost \(W^{1/2}\operatorname{polylog}W\), exponential here. Constructing a suitable multiplier with the needed residue-spreading statistics is stronger than merely finding a prime. Exhausting all \(w\)'s costs \(N^{\Theta(1)}\). |
| Rejection sampling for coprimality, primality, content, and irreducibility certificates | **PROVED** removable as tests: use gcds, deterministic primality, and deterministic polynomial arithmetic. | Enumerating rejected choices does not derandomize finding a surviving good field; a global small rejection probability says nothing about a particular short list. |
| Dixon seeds \(r\) | **OPEN** at deterministic subexponential cost with guaranteed useful rows and fruitfulness. | Enumeration of all residues costs \(N^{1+o(1)}\); the suggested small seeds fail §2.3. There are separate supply and root-branch problems. |
| Gauss reduction and tie breaking | **PROVED** already deterministic after a fixed convention. | Polynomial in \(\log N\). It supplies no replacement for hidden-root randomness. |
| Special-\(q\) candidate search | **PROVED** already deterministic in the stated algorithm. Enumerate its reduced lattice box, primality-test the permitted second coordinates, and batch-test smoothness. | Per target \(\widetilde O(B+S^2)\) with the given batch method. What remains unknown is that the deterministic field/targets have successful candidates. |
| Lenstra–Pomerance factoring inside the shape test | **OPEN** to derandomize with a known subexponential worst-case bound under GRH. **PROVED** deterministic exponential fallback. | Inputs are at most \(\sqrt{2N}\); the general \(n^{1/5+o(1)}\) algorithm gives \(N^{1/10+o(1)}\) per such input. LP's randomized bound scales as \(L_N[1/2,1/\sqrt2+o(1)]\). |
| Recursive replacement of the shape test | **PROVED** removes the call to a separate randomized factoring black box, conditional on the recursive algorithm's analysis. | It repeats the same unresolved field/seed searches on smaller inputs. It does not turn a randomized algorithm into a deterministic one. |
| Wiedemann/Kaltofen–Saunders preconditioners and random kernel vectors | **PROVED** removable unconditionally on an explicit matrix. Compute a basis by deterministic elimination and test it as in §3.2. | Dense binary elimination on \(M_r=L^{a+o(1)}\) rows costs \(M_r^{3+o(1)}\) elementary bit operations up to standard representation factors; root testing costs a polynomial in the explicit sizes. Nonzero \(\Phi\) remains necessary for that relation family. |
| Random quadratic-character prime ideals | **PROVED** removable under the needed Hecke GRH at \(L^{O(1)}\) cost for fixed degree and the bounded-height relations here. Enumerate enough small unramified degree-one primes to separate every relevant nonsquare class, and retain independent character constraints. | A common discriminant/height bound gives a prime cutoff polynomial in the total explicit bit size, hence \(L^{O(1)}\). Excluding atom supports can enlarge the cutoff beyond polylogarithmic size. The original leading constant is not claimed. See the proof below. |
| Hensel prime for the number-field square root | **PROVED** already deterministic as specified: take the least admissible prime above \(B\), with a fixed tie rule. | Sieve/search in a subexponential-size interval. Since \(\log|f_d\operatorname{disc}f|=O_d(\log N)\) while \(B\) is subexponential, not every prime in \((B,2B]\) can divide that fixed integer for large \(N\). |
| Roots and polynomial factors modulo small factor-base/Hensel primes | **PROVED** removable unconditionally at \(L^{O(1)}\) for fixed degree by exhaustive field arithmetic or root enumeration; under GRH use deterministic finite-field factoring. | Brute force at primes of size \(O(B)\) can cost \(B^{O_d(1)}\). For the original very large character primes, brute force is not acceptable; use GRH finite-field factoring or the preceding small-prime replacement. |
| Hensel lifting, selection of local square-root sign vectors, exact root verification | **PROVED** deterministic. Enumerate the finitely many local sign choices and verify the exact equation. | Fixed degree gives at most \(2^d\) sign patterns and polynomial work in the expanded coefficient bit sizes. This computes a root of a genuine field square; it does not address modular fruitfulness. |
| Random primes/rank profiles in exact \(\mathbb Z\)-linear algebra, such as XK | **PROVED** removable unconditionally by deterministic rational linear algebra, HNF, or fraction-free methods; certified real approximation and LLL are deterministic. | Polynomial in explicit matrix dimensions and bit sizes, hence \(L^{O(1)}\) here. This does not establish full rank of the resulting seed-relation lattice. |
| Restarts, retry clocks, caching, and stochastic deepening | **PROVED** deterministic scheduling is possible; **OPEN** to obtain the same successful short search. | Dovetail fields, targets, and work budgets. The schedule only helps if a successful bounded deterministic branch is known to lie in the enumerated set. Random expected-time estimates alone do not prove that. |

**PROVED: justification for deterministic character separation under GRH.** For any nonsquare certificate class \(u\in K^\times/K^{\times2}\), the quadratic extension \(K(\sqrt u)/K\) has a nontrivial quadratic Hecke character. A GRH effective Chebotarev/least-prime theorem finds an unramified degree-one prime where that character is \(-1\), with norm polynomial in \(\log|\Delta_{K(\sqrt u)}|\) and the logarithmic size of excluded supports, for fixed degree. The exclusion avoids undefined row evaluations.

Every product of at most \(M_r\) bounded-height atoms has height with logarithm bounded by a polynomial in the total expanded input size. The extension discriminants and excluded supports therefore have one uniform logarithmic bound of size \(L^{O(1)}\). Enumerating all degree-one primes up to the corresponding uniform cutoff separates *all* nonsquare classes represented by possible dependencies; there is no union bound over their exponentially many names. Character evaluation is linear on the parity dependency space. Deterministic elimination can keep a basis of its restrictions. The resulting kernel contains exactly the square classes needed for exact root extraction. Enumeration, factorization over the auxiliary residue fields, and linear algebra all cost \(L^{O(1)}\). The effective GRH prime mechanism is the same kind used in [Schoof, Theorem 6.1](https://www.mat.uniroma2.it/~schoof/mcpaper.pdf); deterministic finite-field factoring is available with fixed-degree polynomial bit complexity under [Evdokimov, 1994](https://doi.org/10.1007/3-540-58691-1_58).

**PROVED: exact algebraic root extraction itself is not a fundamental random step.** Once the radicand is given explicitly and is a square in a fixed-degree field, a deterministic algebraic algorithm can recover and verify its root in time polynomial in the explicit bit size; the described bounded-prime Hensel/sign enumeration is one such approach in this setting. The difficult square-class kernel and the subsequent modular branch are separate issues.

**OPEN, not “irreducible” as a lower bound.** No listed source of randomness is proved necessary for all algorithms. Within this proof architecture the unresolved arithmetic tasks are deterministic field/multiplier supply, useful target supply, and the fruitfulness or full-rank-relation property. Recursive shape testing couples back to those tasks. Kernel sampling, local root routines, and routine random linear algebra are dispensable implementation choices at an \(L^{O(1)}\) cost.

## 7. Literature check and precise distinctions

### 7.1 Deterministic factoring bounds

**PROVED — historical conditional theorem.** René J. Schoof, “Quadratic fields and factorization,” in H. W. Lenstra Jr. and R. Tijdeman, eds., *Computational Methods in Number Theory, Part II*, Mathematical Centre Tracts **155** (1982), 235–286, gives GRH-conditional deterministic \(O_\varepsilon(N^{1/5+\varepsilon})\) factoring algorithms. They use quadratic class groups and real-quadratic infrastructure; see Algorithms 6.4 and 6.5. These are exponential-time algorithms in the input length, not deterministic subexponential index calculus. [Author's scan](https://www.mat.uniroma2.it/~schoof/mcpaper.pdf).

**PROVED — modern unconditional improvement of the baseline.** David Harvey, “An exponent one-fifth algorithm for deterministic integer factorisation,” *Mathematics of Computation* **90** (2021), 2937–2950, proves \(N^{1/5+o(1)}\) deterministic factoring without GRH. David Harvey and Markus Hittmeir, “A log-log speedup for exponent one-fifth deterministic integer factorisation,” *Mathematics of Computation* **91** (2022), 1367–1379, give

\[
O\!\left(\frac{N^{1/5}\log^{16/5}N}{(\log\log N)^{3/5}}\right)
\]

bit operations. Thus “the current one-fifth bound requires GRH” is FALSE. [Harvey](https://arxiv.org/abs/2010.05450), [Harvey–Hittmeir](https://arxiv.org/abs/2105.11105).

**OPEN — literature search finding.** I found no established uniform deterministic subexponential complete-factoring algorithm under GRH for every \(N\), or for a density-one set of integers without an additional oracle or hypothesis. The one-fifth exponent remains the general deterministic benchmark in the sources checked. This is not a proof that such an algorithm is impossible. “Almost all” must specify the ensemble and whether the output is one divisor or complete factorization; semiprimes have density zero among all integers, so these formulations cannot be interchanged.

**PROVED: one-divisor factoring for almost all unrestricted integers is already easy.** Trial-divide by all primes up to \(\log N\). This takes polynomial time in \(\log N\) and finds a proper divisor for a density-one set of integers. To prove density one without a uniform sieve estimate, fix \(B\). The density of integers divisible by none of the primes at most \(B\) is \(\prod_{p\le B}(1-1/p)\), by periodicity. Eventually \(\log N\ge B\), so the upper density of inputs on which this trial division finds no divisor is at most that product. Letting \(B\to\infty\) makes it zero. Finitely many tiny inputs do not matter. This observation supplies neither complete factorization nor a result for balanced semiprimes or the sifted residuals in Theorem A. The nontrivial almost-all objective must use complete factorization or a specified hard input ensemble.

The analytic and group-theoretic references used in the reductions are N. C. Ankeny, “The least quadratic non residue,” *Annals of Mathematics* (2) **55** (1952), 65–72; E. Bach, “Explicit bounds for primality testing and related problems,” *Mathematics of Computation* **55**, no. 191 (1990), 355–380, [DOI](https://doi.org/10.1090/S0025-5718-1990-1023756-8); and G. L. Miller, “Riemann's hypothesis and tests for primality,” *Journal of Computer and System Sciences* **13**, no. 3 (1976), 300–317, [author's PDF](https://www.cs.cmu.edu/~glmiller/Publications/Papers/Mi76.pdf).

Classical precursors to Schoof's treatment are D. Shanks, “Class number, a theory of factorization, and genera,” *Proceedings of Symposia in Pure Mathematics* **20** (1971), 415–440, [AMS volume](https://doi.org/10.1090/pspum/020/0318071), and H. W. Lenstra Jr., “On the calculation of regulators and class numbers of quadratic fields,” in *Journées Arithmétiques 1980*, LMS Lecture Note Series **56** (1982), 123–150, [author's bibliography](https://pub.math.leidenuniv.nl/~lenstrahw/PUBLICATIONS/pub.html). Small-generator bounds and baby-step/giant-step or infrastructure calculations in this lineage do not constitute a deterministic subexponential relation-collection theorem.

### 7.2 What the randomized and class-group papers establish

The following are separate results; there is no single deterministic “Pomerance–Seysen” theorem implied by their combination.

| Exact reference | Established statement relevant here |
|---|---|
| J. D. Dixon, “Asymptotically fast factorization of integers,” *Math. Comp.* **36** (1981), 255–260. [DOI](https://doi.org/10.1090/S0025-5718-1981-0595059-1). | **PROVED:** rigorous randomized random-squares factoring. The root-fibre argument is probabilistic. |
| C.-P. Schnorr, “Refined analysis and improvements on some factoring algorithms,” *Journal of Algorithms* **3** (1982), 101–127. [Publisher](https://www.sciencedirect.com/science/article/pii/0196677482900128). | **PROVED:** randomized and heuristic variants, with a refined rigorous analysis of Dixon's Monte Carlo method. It does not establish deterministic subexponential factoring. |
| C.-P. Schnorr and H. W. Lenstra Jr., “A Monte Carlo factoring algorithm with linear storage,” *Math. Comp.* **43** (1984), 289–311. [Author's PDF](https://pub.math.leidenuniv.nl/~lenstrahw/PUBLICATIONS/1984c/artc.pdf). | **PROVED:** a Monte Carlo class-group/order method. Its proposed fast behavior is not a deterministic GRH theorem. |
| M. Seysen, “A probabilistic factorization algorithm with quadratic forms of negative discriminant,” *Math. Comp.* **48** (1987), 757–780. [DOI](https://doi.org/10.1090/S0025-5718-1987-0878705-X). | **PROVED under ERH:** randomized class-group factoring, with \(L_N[1/2,\sqrt{5/4}]\)-type running time. In particular the word “probabilistic” is substantive. |
| A. K. Lenstra, “Fast and rigorous factorization under the generalized Riemann hypothesis,” *Indagationes Mathematicae (Proceedings)* **91** (1988), 443–454; also cited as *Indag. Math.* **50**. [Publisher](https://www.sciencedirect.com/science/article/pii/S1385725888800222). | **PROVED under GRH:** randomized \(L_N[1/2,1+o(1)]\) factoring, combining Seysen's method with a smoothness test. Not deterministic. |
| H. W. Lenstra Jr. and C. Pomerance, “A rigorous time bound for factoring integers,” *J. Amer. Math. Soc.* **5** (1992), 483–516. [PDF](https://community.ams.org/journals/jams/1992-05-03/S0894-0347-1992-1137100-0/S0894-0347-1992-1137100-0.pdf). | **PROVED:** randomized expected \(L_N[1/2,1+o(1)]\) factoring without GRH. Exhausting its random tapes is not a subexponential derandomization. |
| J. L. Hafner and K. S. McCurley, “A rigorous subexponential algorithm for computation of class groups,” *J. Amer. Math. Soc.* **2** (1989), 837–850. [Author institution](https://research.ibm.com/publications/a-rigorous-subexponential-algorithm-for-computation-of-class-groups). | **PROVED under ERH:** a Las Vegas algorithm for imaginary-quadratic class-group structure in expected \(L_{|\Delta|}[1/2,\sqrt2+o(1)]\) time. It is not deterministic. |
| J. Buchmann, “A subexponential algorithm for the determination of class groups and regulators of algebraic number fields,” in *Séminaire de Théorie des Nombres, Paris 1988–1989*, Progress in Mathematics **91** (1990), 27–41. | **HEURISTIC running-time analysis:** a probabilistic generalization to number fields, not a deterministic theorem from GRH alone. The historical distinction is reviewed in [de Boer–Pellet-Mary–Wesolowski, §12.2](https://koendeboer.com/publication/rigorousmethods/RigorousMethods.pdf). |
| K. de Boer, A. Pellet-Mary and B. Wesolowski, “Rigorous methods for computational number theory,” arXiv:**2512.01588**, 2025, Theorem 1.3. [Preprint](https://arxiv.org/abs/2512.01588). | **PROVED under ERH, as stated in this recent preprint:** a probabilistic rigorous class/unit-group algorithm, including fixed-degree \(L_{|\Delta|}[1/2,O(1)]\) bounds. It removes old heuristic gaps in the randomized setting, not randomness. |

**PROVED as a reduction.** A computed quadratic class-group presentation with explicit ideal/form representatives can expose its 2-torsion. Ambiguous forms encode factorizations of the discriminant; this is the mechanism used by Schoof. Small generators under GRH do not by themselves compute their orders or the relation lattice. Therefore small-generation theorems do not turn Hafner–McCurley or Buchmann-style relation collection into deterministic subexponential algorithms.

**OPEN — literature search finding.** I found no general deterministic subexponential algorithm under GRH alone for the full class-group structure problem as discriminant grows, or for discrete logarithms in arbitrary prime fields \(\mathbb F_p\). The classical rigorous prime-field index-calculus algorithm is randomized: C. Pomerance, “Fast, rigorous factorization and discrete logarithm algorithms,” in *Discrete Algorithms and Complexity*, Academic Press (1987), 119–143. [Author's scan](https://math.dartmouth.edu/~carlp/disclog.pdf). Results for special groups, smooth group orders, bounded characteristic, or finite-field *polynomial factorization* do not answer that prime-field discrete-log question.

### 7.3 NFS/QS fruitfulness under GRH

**PROVED statement of what the cited paper claims.** Jonathan D. Lee and Ramarathnam Venkatesan, “Rigorous Analysis of a Randomised Number Field Sieve,” *Journal of Number Theory* **187** (2018), 92–159, DOI 10.1016/j.jnt.2017.10.019, analyzes a randomized polynomial family. Its square-congruence theorem is distinct from its factoring theorem; the latter, for Blum semiprimes, assumes an additional character nondegeneracy conjecture. That conjecture is not stated as a consequence of GRH. [Paper](https://arxiv.org/abs/1805.08873), [full text, §7](https://arxiv.org/html/1805.08873v1).

**OPEN — literature search finding.** I found no published theorem proving the required deterministic NFS or general QS dependency fruitfulness merely from Dirichlet/Hecke GRH. The class-group algorithms above have their own proved group-theoretic mechanisms; those cannot be imported into an arbitrary NFS dependency matrix without proving the corresponding relation-spanning or nondegeneracy assertion.

### 7.4 “Factoring reduces to modular square roots”: which problem?

**PROVED — unconditional randomized reduction.** M. O. Rabin, “Digitalized signatures and public-key functions as intractable as factorization,” MIT/LCS/TR-212 (1979), uses a root algorithm independent of the chosen hidden root. Choose uniform \(r\), query a root of \(r^2\), and compare. The fibre calculation in §3.6 proves success probability \(1-2^{1-\omega}\). No distributional assumption on the root algorithm's output is needed. [MIT record](https://dspace.mit.edu/entities/publication/15016449-ca26-4d45-82ea-7c200d82d777).

**PROVED — deterministic GRH reduction to a total root-or-factor problem.** E. Jeřábek, “Integer factoring and modular square roots,” *Journal of Computer and System Sciences* **82** (2016), 380–394, Theorem 3.8, gives deterministic reductions under quadratic GRH involving

\[
\mathrm{FacRoot}:\quad (a/N)=1
\ \Longrightarrow\ \text{produce a root of }a\text{ or a proper factor of }N.
\]

Jacobi \(+1\) allows nonsquares. The reduction finds such a small nonsquare, forcing a factor. The total problem Root may return a verifiable nonsquare witness, also giving a proper divisor on these inputs. [Theorem and definitions](https://arxiv.org/html/1207.5220v3).

**OPEN here:** no verified source in this search establishes the unrestricted promise-only deterministic equivalence. A promise to work only on actual squares is weaker than FacRoot. Neither that total-problem theorem nor the small-seed selector theorem supplies an arbitrary-selector deterministic reduction.

### 7.5 Barriers versus impossibility

**PROVED methodological barriers:** the absence of a global homomorphic selector at primes \(1\pmod4\); the small-seed shape rejection; the rowspace criterion for structural fruitlessness; the failure of fixed-test quantifiers; and the inability of global moments or generic expansion to locate a prescribed short list. These are exact statements with specified scopes.

**OPEN — no general impossibility claim.** I found no impossibility theorem ruling out deterministic subexponential Dixon-type factoring in the ordinary integer bit model under GRH. Generic-group bounds or lower bounds for particular proof techniques would not give such a theorem. Conversely, a sufficiently strong uniform circuit-hardness assumption derandomizes the known Las Vegas factoring algorithm as in §5.6. Thus none of the obstructions above licenses the claim that randomness is inherently necessary.

**Finite verification.** The companion [verification script](verify_report.py), using only Python's standard library, checks 150,498 selector multiplication identities for 58 odd moduli up to 250, the order-4 obstruction for the other 66 moduli in that range, 1,508 shortest-vector instances up to 1,000, and 2,028 exponent-multiple splitting runs up to 2,000, including repeated prime factors and extra even/odd multipliers. All passed. These checks illustrate the elementary statements; the arguments above are their proofs.

## 8. Final claims, verdict, and one remaining hypothesis

### (i) Claims considered PROVED

1. Certificate squarehood, exact certificate relations, and the finite modular kernel are different objects. The half-relation map (1.6) precisely identifies the root defect when the relevant roots are represented.
2. Formal multiplicativity descends to residue classes only after consistency with all relevant modular relations. A deterministic pivot convention does not prove that consistency.
3. A homomorphic square-root selector \(G^2\to G\) exists exactly when every odd prime divisor of \(N\) is \(3\pmod4\). In that case it is unique and its root-ratio components are the local Legendre symbols.
4. Under GRH, a nonzero quotient character arising from such a selector has a fruitful unit base below \(2\log^2N\). For \(\omega\ge2\) the genuine global homomorphic selector necessarily has a nonzero quotient character.
5. For \(r^4+1<N\), the stated shortest-vector split is \(\pm(r^2,1)\). The proposed small initial segment therefore fails Theorem A's shape condition; accepting its rational-square certificates instead gives tautologies.
6. With the stated evaluability/torsion conditions, \(\Phi\) is a homomorphism on the true dependency space. If it is nonzero, testing a basis finds a factor deterministically, without GRH or a new seed.
7. Structural fruitlessness is exactly \(\Phi=0\), equivalently the appropriate factor-detecting rows belong to the square-class rowspace. The B1 nonresidue condition is only a restricted sufficient obstruction to fruitfulness, not its converse.
8. The DD probability follows from root-fibre uniformity and transcript invariance. Deterministic small roots do not have that distribution.
9. A full-rank lattice of verified modular relations on generators yields a multiple of \(\varphi(N)\), hence an exponent multiple, even when it is a proper unsaturated sublattice of the full kernel.
10. Given an exponent multiple, the explicit small-base chain in §4.2 factors deterministically under quadratic GRH in polynomial time in the input and exponent bit lengths.
11. A primitive exact NFS dependency alone does not imply full rank on a fixed set of modular generators.
12. Global FM/LD*/Paley–Zygmund bounds do not imply that a prescribed subexponential list contains a good field. The supplied box and multiplier estimates operate at power-of-\(N\) scales; restricted-list estimates require new input.
13. Kernel sampling, ordinary random linear algebra, bounded-prime local root routines, and—under Hecke GRH—quadratic-character separation can be made deterministic at an \(L^{O(1)}\) cost in the stated fixed-degree setting.
14. The modern deterministic \(N^{1/5+o(1)}\) factoring bound is unconditional. The cited rigorous subexponential factoring/class-group algorithms are randomized; the cited deterministic GRH root-or-factor reduction has a stronger input/output specification than a promise-square selector.
15. Finding one proper divisor for almost all unrestricted integers is already deterministic polynomial time by small-prime trial division. The substantive density-one question concerns complete factorization or an explicitly restricted ensemble.

### (ii) Verdict

**OPEN: deterministic subexponential complete factoring under GRH does not follow from the proposed argument or the supplied average relation-supply statements.** The multiplicative-selector idea is valid only after imposing compatibility that is impossible for many \(N\), and unconstructed in the remaining case. For the actual descent, the suggested small seeds fail even before fruitfulness. A short good-field list would solve one search problem, but would not by itself prove useful deterministic targets, nonzero \(\Phi\), or full rank of exact modular relations. The same gaps remain for a complete-factoring almost-all claim until an input ensemble and the requisite exceptional-set estimates are proved. If the objective instead means just one proper divisor on a density-one subset of all integers, §7.1 already gives that unconditionally in polynomial time.

### (iii) The single most promising remaining hypothesis

**OPEN — deterministic full-rank modular relation supply.** There is an absolute constant \(C\) and a uniform deterministic algorithm with the following property. On input an odd composite \(N\) that is not a prime power, let \(b_1,\ldots,b_m\) consist of \(-1\) and all positive unit bases at most \(\lceil2\log^2N\rceil\); a nonunit gcd may instead be returned as a factor. In time at most

\[
\exp\!\bigl(C\sqrt{\log N\log\log N}\bigr),
\]

the algorithm either finds a proper factor or constructs, from exact NFS/descent certificates or another explicitly specified relation procedure, an integer matrix \(A\in\mathbb Z^{m\times m}\) of total encoded length at most the same bound such that

\[
\det A\ne0,
\qquad
\prod_{j=1}^m b_j^{A_{ij}}\equiv1\pmod N
\quad(1\le i\le m).
\]

Every output condition is directly verifiable. This hypothesis plus Dirichlet GRH proves deterministic \(L_N[1/2,O(1)]\) complete factoring by §§4.1–4.2; no separate NFS fruitfulness hypothesis or saturation theorem is needed. It is a precise research target, not a claimed consequence of the current moment theorems, and it is not advertised as weaker than factoring: under GRH the reduction makes its computational strength essentially equivalent to the desired result. Its value is to replace an ambiguous modular square-root choice by one explicit, checkable missing relation-spanning theorem.
