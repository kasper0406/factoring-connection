**Nested descendability, ERH, and the \(L[1/3]\) factoring question**

Independent research report, 20 October 2026.

**Result: OPEN.** I do not prove an ERH-conditional \(L_N[1/3,c]\) factoring algorithm, or an algorithm for every integer with any fixed exponent \(\alpha<1/2\). The strongest general factoring bound I can support here is the unconditional expected bound \(L_N[1/2,1+o(1)]\) of Lenstra–Pomerance. The principal obstruction is a quantitative supply theorem for short, two-sided relations whose factors have completed descent certificates in the same field. dBPW supplies a different, full-dimensional ideal sampler; its discriminant-dependent norm scale prevents a direct substitution.

There are positive results below: an adaptive one-sided consequence of dBPW, an elementary restricted-ideal counting lemma, a uniform split-counting argument, the adaptive-tree version of transcript invariance, and conditional cost recurrences. None establishes the missing relation supply.

**Scope and conventions.** I read the permitted theorem roadmap and all of recurse_results §§3.1–3.5. I did not read the prohibited project reports or other agents' work. I used the original dBPW PDF, version 2, dated February 2026, through the browser's PDF fetch and extraction. A local curl download failed because the shell could not resolve arxiv.org; there is no local PDF copy. The source is [dBPW, arXiv:2512.01588v2](https://arxiv.org/pdf/2512.01588v2).

PROVED means either a proof is supplied below, or an explicitly identified external theorem is invoked under its stated hypotheses. It does not mean that I independently audited every line of a 154-page preprint. HEURISTIC marks a distributional assumption or an algorithmic prediction. OPEN marks a missing proof. A proof of an implication from additional hypotheses is not a proof that those hypotheses hold. In particular, the internal FM, LD*, and growing-degree LV6 claims are not assumed verified.

Put
\[
 t=\log N,\qquad D=t/\log t,\qquad
 E=t^{1/3}(\log t)^{2/3},\qquad
 Q=t^{2/3}(\log t)^{1/3},\qquad
 D'=Q/E=D^{1/3}.
\]
Thus \(E=D'\log t\), \(Q=(D')^2\log t\), and \(D'Q=t\). Write
\[
 L_N[\alpha,c]=
 \exp\!\big((c+o(1))t^\alpha(\log t)^{1-\alpha}\big).
\]
The letter \(D\) here is not a field discriminant.

**1. What dBPW actually supplies**

The following formula register records the relevant hypotheses and conclusions of Definitions 2.11–2.13, Theorems 1.1, 7.1, 9.5, 17.7, and Remark 9.15. Throughout, \(n=[K:\mathbb Q]\), \(\Delta=|\Delta_K|\), and \(\rho_K=\operatorname{res}_{s=1}\zeta_K(s)\). The degree is unrestricted; the algorithmic lattice block size is \(2\le b\le n\). The representation uses a monic irreducible defining polynomial of bit length polynomial in \(\log\Delta\), an LLL-reduced integral basis, and rational coordinates for elements and ideal bases. These are external ERH-conditional results, except for the Haar-average assertion in Theorem 7.1. [dBPW, pp. 3, 10–13, 28, 47–50, 57, 87–88](https://arxiv.org/pdf/2512.01588v2).

**Density and geometry.**
\[
 \delta_{\mathcal S}[x]
 =\min_{x/e^n\le v\le x}
 \frac{\#\{\mathfrak c\in\mathcal S:N\mathfrak c\le v\}}
      {\rho_K v}.
 \tag{1}
\]
\[
 K_{\mathbb R}\simeq\mathbb R^{r_1}\times\mathbb C^{r_2},\quad
 n=r_1+2r_2,\quad
 r\mathcal B_\infty=\{(z_\sigma):|z_\sigma|\le r\},\quad
 \Gamma_K=\sup_{x\mathfrak a}
       \frac{\lambda_n(x\mathfrak a)}{\lambda_1(x\mathfrak a)}
 \le\Delta^{1/n}.
 \tag{2}
\]
Here the supremum is over ideal lattices; complex coordinates contribute twice to the norm.

**Theorem 1.1, as printed in the introduction.**
\[
\begin{gathered}
 \mathfrak a\subseteq\mathcal O_K,\quad
 \text{LLL-reduced basis of }\mathcal O_K\text{ supplied},\\
 r\ge48\,b^{2n/b}n^{7/2}\Delta^{3/(2n)}
                    (N\mathfrak a)^{1/n},\\
 \Pr\big[\beta\in\mathfrak a,\
     (\beta)\mathfrak a^{-1}\in\mathcal S\mathcal S_{B_0}\big]
 \ge\delta_{\mathcal S}[r^n]/3-\varepsilon,\\
 B_0=(\log\Delta+\log\log(1/\varepsilon))^{O(1)},\\
 \mathbb E T=\operatorname{poly}
   (\text{input length},\log\Delta,\log N\mathfrak a,
                  \log(1/\varepsilon),b^b).
\end{gathered}
\tag{3}
\]
The useful accuracy range has \(0<\varepsilon<1\).

**Theorem 7.1.** Let \(P=\operatorname{Pic}^0_{K\mathfrak m}\), \(G\le P\) have finite index \(I\), and let an Arakelov divisor \(\mathfrak d\) have class Haar-uniform on \(G+[\mathfrak b]\). Write \(\mathfrak m=\mathfrak m_0\mathfrak m_\infty\), \(N\mathfrak m=2^{|\mathfrak m_{\mathbb R}|}N\mathfrak m_0\). Then
\[
\begin{gathered}
 r\ge8n^2\Gamma_K\Delta^{1/(2n)}(N\mathfrak m)^{1/n},\\
 (\mathcal S,\mathfrak m_0)=1,\qquad
 [d^0(\mathcal S)]\subseteq
       G+[\langle\tau\rangle]-[\mathfrak b],\\
 \mathbb E_{\mathfrak d}\Pr_{\alpha\leftarrow
                 \operatorname{Exp}(\mathfrak d)\cap r\mathcal B_\infty}
 \left[
 (\alpha)\operatorname{Exp}(-\mathfrak d)\in\mathcal S
 \ \middle|\
 \alpha\operatorname{Exp}_\infty(-\mathfrak d_\infty)
                         \in\tau K_{\mathfrak m,1}
 \right]
 \ge \frac{N\mathfrak m_0}{\phi(\mathfrak m_0)}
             \frac I3\,\delta_{\mathcal S}[r^n].
\end{gathered}
\tag{4}
\]
The conditioning prescribes finite residues and real signs.

**Theorem 9.5.** Inputs additionally include a subgroup-membership oracle, rationally represented archimedean distortion \(y\), \(\omega\ge1\), and
\[
\begin{gathered}
 (\mathfrak a,\mathfrak m_0)=1,\quad
 \tau\in K_{\mathfrak m},\quad[\langle\tau\rangle]\in G,\quad
 0<\varepsilon<\min(1,20/n),\\
 (\mathcal S,\mathfrak m_0)=1,\qquad
 [d^0(\mathcal S)]\subseteq G-[d^0(\mathfrak a)],\\
 r=48\omega b^{2n/b}n^{7/2}\Delta^{3/(2n)}
                  (N\mathfrak m_0)^{1/n},\\
 B_0=\widetilde O\!\left(
 I^2n^2\{n^2(\log\log(1/\varepsilon))^2+
                     (\log(\Delta N\mathfrak m))^2\}\right),\\
 T_{\rm rw}=\left\lceil
 7n+\log N\mathfrak m+\log|P|-\log I+
                 2\log(1/\varepsilon)+1\right\rceil.
\end{gathered}
\tag{5}
\]
Its joint success guarantee is
\[
\begin{gathered}
 \beta\in\mathfrak a,\quad
 (\beta)\mathfrak a^{-1}\in
                    \mathcal S\mathcal S_{B_0}^{G},\quad
 \beta\in\tau K_{\mathfrak m,1},\\
 |N_{K/\mathbb Q}\beta|
       \le N\mathfrak a\,B_0^{T_{\rm rw}}r^n,\\
 p_{\rm success}\ge
 \frac{N\mathfrak m_0}{\phi(\mathfrak m_0)}
             \frac I3\,\delta_{\mathcal S}[r^n]-\varepsilon.
\end{gathered}
\tag{6}
\]
Here \(\mathcal S_{B_0}^{G}\) uses primes of norm at most \(B_0\), coprime to \(\mathfrak m_0\), with normalized Arakelov classes in \(G\). Complexity is polynomial in \(b^b,\log\Delta,\operatorname{size}(\mathfrak a),\log(1/\varepsilon),\log N\mathfrak m,I,\operatorname{size}(y),\operatorname{size}(\tau),\operatorname{size}(\omega)\); oracle calls are bounded by \(I\operatorname{poly}(\log\Delta,\log N\mathfrak m)\).

**Theorem 17.7.** With the algorithm's small-prime exclusion modulus \(\mathfrak m_0\), let
\[
\begin{split}
 \rho_K^{\rm cut}
 &=\min\left(\rho_K,
 \max\left\{
 e^{(\log\Delta)^{2/3}(\log\log\Delta)^{2/3}},
 e^{n^{2/3}(\log n)^{4/3}}
 \right\}\right),\\
 B_{\max}&=\operatorname{poly}
       (L_\Delta[1/2,O(1)],L_{n^n}[2/3,O(1)],\rho_K^{\rm cut}).
\end{split}
\tag{7}
\]
If \(S\) contains every prime ideal of norm at most \(B_{\max}\) coprime to \(\mathfrak m_0\), and \((\mathfrak a,\mathfrak m_0)=1\), Algorithm 3 returns
\[
 \alpha\in\mathfrak a,\qquad
 \mathfrak a=(\alpha)\prod_{\mathfrak p\in S}\mathfrak p^{-v_{\mathfrak p}},
 \quad v_{\mathfrak p}\ge0,
\tag{8}
\]
in expected time polynomial in the quantities in (7), input sizes, and \(\log\Delta\).

**Interpretation and consequences — PROVED from the register.**

First, (1) counts *all integral ideals up to varying norm bounds*. It is not a relative density among rational primes, degree-one primes, ideals in one dyadic window, or values of a binary form. In particular, a polynomial relative density of acceptable primes is not the hypothesis \(\delta_{\mathcal S}[r^n]\ge t^{-C}\) for their multiplicative monoid.

Second, the sampling is in real dimension \(n\). The rank-two lattice
\[
 \{a-b\alpha:a,b\in\mathbb Z\}\subset K
\tag{9}
\]
has smaller rank for \(n>2\); it is not a full ideal lattice. A nonzero \(\mathcal O_K\)-ideal has rank \(n\), because it contains \(c\mathcal O_K\) for a nonzero integer \(c\). Thus (9) cannot be used as the input ideal in (3) or (5). The exception \(n=2\) does not give the growing-degree NFS construction.

Conditioning a full-dimensional sample on lying in (9) does not preserve its success bound. In a coefficient cube of radius \(H\), the slice has \(O(H^2)\) points versus \(\Theta(H^n)\) points in total. This example explains the loss under conditioning; it is not an asserted exact probability formula for the distorted dBPW boxes.

Third, (3) and (5) have different normalizations. The introduction's radius includes \((N\mathfrak a)^{1/n}\). In the stronger detailed theorem, \(r\) is normalized; the actual lattice box is scaled by \((N\overline{\mathfrak a})^{1/n}\), with \(\overline{\mathfrak a}=\mathfrak a\prod\mathfrak p_i\). One should use (5)–(6) for a norm calculation, rather than silently identify the two radii.

Taking norms in (6) gives the explicit bound
\[
 N((\beta)\mathfrak a^{-1})
 \le B_0^{T_{\rm rw}}(48\omega)^n
 b^{2n^2/b}n^{7n/2}\Delta^{3/2}N\mathfrak m_0.
\tag{10}
\]
The cancellation of \(N\mathfrak a\) is significant. A small target ideal does not make the certified cofactor bound proportionally small.

More precisely, the proof constructs
\[
 (\beta)\mathfrak a^{-1}
       =\mathfrak c\prod_{i=1}^{T_{\rm rw}}\mathfrak p_i,\qquad
 N\mathfrak c\le r^n,
\tag{11}
\]
on the success event, with the walk primes already known and \(B_0\)-smooth. The full bound in (10) can be larger than a fixed power of \(N\), since \(T_{\rm rw}\log B_0\) can be \(O(t\log t)\). The unknown bulk in (11) has the cleaner \(r^n\) bound. Discarding the known walk factors does not discard the \(\Delta^{3/2}\) term.

Equation (4) is a geometric averaging statement, conditional on an already uniform Arakelov class. It does not require ERH. The algorithm obtains approximate uniformity using the ERH random walk, then pays for lattice sampling. Neither this sequence of operations nor (8) is a tower of contracting special-\(q\) descents.

The printed input assumptions should also be respected. However, computing the maximal order is not an unavoidable separate obstruction to *Part 1*: §9.5 gives an order variant away from its conductor. Its radius uses
\[
 |\Delta(R)|^{3/(2n)},\qquad
 |\Delta(R)|=[\mathcal O_K:R]^2\Delta,
\tag{12}
\]
and requires the modulus to include the conductor and the appropriate coprimality conditions. This can avoid maximal-order computation, at a worse discriminant scale. It does not make a rank-two slice full-dimensional. [dBPW, Lemma 9.14 and Remark 9.15](https://arxiv.org/pdf/2512.01588v2#page=57).

**The NFS norm calculation — PROVED upper bounds; typical sizes HEURISTIC.**

Consider first a monic LV polynomial of degree \(n=\delta D'\) and height \(H\le e^{\kappa Q}\). A discriminant bound follows directly from Hadamard's inequality for the Sylvester determinant of \(f,f'\):
\[
 \log|\operatorname{disc}f|
       =O(n\log H+n\log n)=O(t).
\tag{13}
\]
Since \(\Delta_K\) divides the discriminant of the power-basis order, \(\log\Delta_K=O(t)\). Nonmonic representations require a specified integral order; any extra leading-coefficient or index contribution must be included. The monic case already exhibits the problem.

A known root \(m\bmod N\) gives an evaluation map on the defining order. It does not give \(\log\Delta_K=O(Q)\). The equality \(\log\Delta_K\asymp t\) for the particular random fields is a typical-size expectation, not a lower bound established here.

For any field with \(\log\Delta_K=\zeta t\), \(\zeta>0\), and trivial ray modulus, the central bound is
\[
 \log(r^n)
 =\frac32\zeta t+
       \frac{2n^2}{b}\log b+\frac72n\log n+n\log(48\omega).
\tag{14}
\]
Taking \(b=2\) makes the extra terms \(o(t)\) at \(n\asymp D'\). Thus the certified unknown cofactor scale is \(N^{3\zeta/2+o(1)}\), even for a target prime of size \(e^{O(Q)}\). This is an upper-bound scale, not a claim that each sample has that norm.

Even improving the discriminant exponent from \(3/2\) to \(1/2\), as the geometric radius might suggest in a particularly well-shaped field, would leave a term of order \(t\). To certify a usual \(L[1/3]\) smoothness cost by this route, one would need the unknown logarithmic norm scale to be \(O(Q)\), together with uniform ideal-density bounds. A uniform construction of the required known-root fields at that effective discriminant scale is **OPEN**.

**2. The strongest nested statement obtained here**

**Adaptive one-sided ideal lemma — PROVED under ERH, using Theorem 9.5.**

Fix \(K\), its computational representation, and the parameters of (5). Given a history \(\mathcal H\), allow both the target ideal \(\mathfrak a_{\mathcal H}\) and a family \(\mathcal S_{\mathcal H}\) to be arbitrary functions of that history. Suppose all hypotheses of (5) hold and
\[
 \delta_{\mathcal S_{\mathcal H}}[r^n]\ge\eta
\tag{15}
\]
for every history being considered. Use fresh sampling coins and \(\varepsilon\le\eta/6\). For the full subgroup and trivial modulus, the conditional success probability is at least \(\eta/6\), uniformly in the history. If membership is decidable at cost \(C_{\rm test}\), repetition costs at most \(6(C_{\rm sample}+C_{\rm test})/\eta\) in expectation.

*Proof.* Condition on the entire history. The field, ideal and family are now fixed admissible inputs. Equation (6) gives \(\eta/3-\varepsilon\ge\eta/6\). Repeating with fresh coins gives the geometric bound. The same argument works at every stopping history. No uniform distribution of the endogenous input ideal is needed. \(\square\)

If \(\mathcal S\) is the monoid generated by a desired set of primes containing the walk primes, then \(\mathcal S\mathcal S_{B_0}=\mathcal S\). In that case this is a genuine one-sided result for adaptively arising children.

The missing hypothesis is (15), at the norm scale (14), for the recursively defined prime set. The lemma does not prove it from a relative prime density.

**An elementary restricted-ideal density bound — PROVED.**

Let \(\mathcal P\) consist of \(m\) distinct prime ideals of norm at most \(y\), where
\[
 y\ge2,\quad m\ge\delta y/\log y,\quad0<\delta\le1.
\]
Let \(\mathcal S=\langle\mathcal P\rangle\), \(x\ge e^n y\), and \(u=\log x/\log y\). Then
\[
 \delta_{\mathcal S}[x]\ge
 \frac1{\rho_Ky}
       \left(\frac{\delta}{u\log y}\right)^u.
\tag{16}
\]
No ERH is required once the prime count is supplied.

*Proof.* For \(v\in[x/e^n,x]\), set \(h=\lfloor\log v/\log y\rfloor\ge1\). The products of exactly \(h\) primes from \(\mathcal P\), with repetition, give
\[
 \#\mathcal S(v)\ge {m+h-1\choose h}\ge(m/h)^h.
\]
Distinct exponent vectors give distinct ideals. Since \(v<y^{h+1}\),
\[
 \frac{\#\mathcal S(v)}{\rho_Kv}
 \ge\frac1{\rho_Ky}
       \left(\frac{\delta}{h\log y}\right)^h
 \ge\frac1{\rho_Ky}
       \left(\frac{\delta}{u\log y}\right)^u.
\]
For the last inequality, \(z\log(z\log y/\delta)\) is increasing for \(z\ge1\), \(y\ge2\), \(0<\delta\le1\). Take the minimum over \(v\). \(\square\)

This bound handles arbitrary sets of prime ideals and exposes a quantitative issue: its \(y^{-1}\) loss alone costs \(e^{\Theta(Q)}\) when \(y=e^{\Theta(Q)}\). It is insufficient for the desired \(e^{O(E)}\) bound. Stronger multiscale combinatorics can improve such estimates; even an ideal analogue with an ideal Dickman bound would leave the thin-slice and discriminant problems.

**A precise sufficient two-sided lemma — OPEN.**

There are two kinds of nodes. At a rational node \(q\), require \(q\mid G_A(a,b)\), and the children come from \(G_A(a,b)/q\) and the ideal of \(a-b\alpha\). At an algebraic node \(\mathfrak q=(q,\alpha-s)\), require \(a-bs\equiv0\bmod q\); the children come from \(G_A(a,b)\) and \((a-b\alpha)\mathfrak q^{-1}\). These are different special-\(q\) lattices. One cannot identify the two conditions just because both norms are \(q\).

Work away from leading-coefficient, discriminant and conductor exceptions, with an explicit procedure for those exceptions. For primitive \((a,b)\), the ordinary nonexceptional algebraic prime factors are degree one: if \(p\nmid b\) divides \(F_A(a,b)\), the relevant root is \(a/b\bmod p\).

Fix a finite search budget at each node and a contraction factor \(0<\lambda<1\). Define \(\mathcal H_0(A)\) to contain the handled base primes and ideals. Define \(\mathcal H_j(A)\) to consist of targets admitting a *completed certificate*: a searched point for which every nonbase child has norm at most \(q^\lambda\) and belongs to \(\mathcal H_{j-1}(A)\). Use both node types.

A sufficient analytic conclusion is
\[
 \mathbb E_A
 \frac{|\mathcal H_h^{\rm rat}(A)\cap J|}{|J|}
                  \ge t^{-C},
 \qquad h=O(\log t),
\tag{17}
\]
with \(C\) independent of \(N\) and the level, and with all certificate construction, testing, field exceptions and square-root work bounded as in §4 below. Equivalently, one could establish appropriately normalized first and second moments of *completed-tree counts* yielding (17).

It is essential that these are completed-tree good sets. An average density estimate for the next level's primes, followed by an unproved assertion that a parent lands in that set, is circular. Also, a bound \(t^{-C_j}\) with \(C_j\) growing with \(j\) is not the uniform statement (17).

**Why a bare density principle cannot supply the missing lemma — PROVED limitations.**

There is a simple class obstruction to an unrestricted assertion. Let \(H\) be a proper subgroup of \(\operatorname{Cl}(K)\), and allow only primes whose classes lie in \(H\). If \([\mathfrak q]\notin H\), no identity
\[
 (\beta)=\mathfrak q\prod_{\mathfrak p\in\mathcal P}
                          \mathfrak p^{e_{\mathfrak p}}
\tag{18}
\]
can hold: its class would assert \([\mathfrak q]^{-1}\in H\). When \(H\) has bounded index, this is compatible with positive relative prime density; prime-ideal equidistribution under Hecke ERH makes that explicit. Allowing all sufficiently small base primes can remove this particular obstruction. It does not prove a two-sided smoothness theorem.

An even more general distinction concerns arbitrary *ideal families*. For a finite searched slice, let \(\mathcal D\) be its set of output cofactor ideals, with \(|\mathcal D|\le M\), and put \(\mathcal S\) equal to all integral ideals except \(\mathcal D\). Then the slice never hits \(\mathcal S\), while
\[
 \delta_{\mathcal S}[x]\ge
       \delta_{\rm all}[x]-\frac{e^nM}{\rho_Kx}.
\tag{19}
\]
This follows by deleting at most \(M\) ideals from every cumulative count in (1). Whenever the last term is small, density in the ambient ideal space says little about that slice. This disproves a general transference from arbitrary ideal families to a thin search. It does not disprove (17), whose families have additional multiplicative structure.

Likewise, marginal probabilities cannot establish a joint estimate: indicators \(X=\mathbf1_E\) and \(Y=\mathbf1_{E^c}\) have positive means but \(XY=0\). Averaging over the field does not repair an unproved correlation between acceptable children and selected parent values.

**3. Assessment of the proposed designs**

**(a) Algebraic ideal-class descent with independently randomized rational \(g\): OPEN as a factoring design.**

The ideal-class part can use the adaptive lemma if (15) is available. The proposed independence on the rational side is not present automatically. Given \(K\), \(\alpha\), and an evaluation map \(\psi(\alpha)=m\bmod N\),
\[
 \psi(a-b\alpha)=a-bm
     =g_1^{-1}(g_1a-g_0b)\pmod N
\tag{20}
\]
requires \(g_0\equiv mg_1\bmod N\). Thus \(g\) is constrained to a congruence lattice once the algebraic map is fixed. In the LV construction the field itself changes with \(g\). Initial independence of some coefficient coins is not conditional independence of rational values after the field and selected points have been exposed.

For a general \(\beta=\sum_i c_i\alpha^i\), the residue is \(\sum_i c_im^i\bmod N\). There is no identity supplying a small linear rational value with the same residue. Ideal-class reduction alone does not solve this lifting problem.

**(b) Good sets defined inductively: the definition is valid; the density induction is OPEN.**

The completed-certificate definition above is implementable by capped searches and memoization. It eliminates logical ambiguity about a child's meaning. To prove its density, however, one needs the parent count with factors in that exact field-dependent set. Neither ordinary smoothness moments nor ERH prime counts give this implication.

A quantitative issue persists even if a separate estimate is obtained at each level. If a proof loses a factor \(t^{-c}\) at each of \(h\asymp\log t\) levels, the resulting lower bound is \(e^{-O((\log t)^2)}\). Raised to \(2k\), \(k\asymp D'\), this costs \(e^{O(E\log t)}\), not \(L_N[1/3,c]\) with finite \(c\). A uniform exponent, a summable loss, or a better global selection argument is required.

**(c) dBPW at each level: PROVED one-sided sampling, OPEN useful contracting factoring descent.**

The direct substitution pays (10) at every node. The field discriminant does not shrink with the target prime. In the \(N^{\Theta(1)}\) discriminant regime, this is not a certified geometric contraction.

Hidden-square rows *can* absorb a general algebraic element if its residue is prescribed. Suppose an evaluation map extends to the order used for sampling and has kernel \(\mathfrak n\), with quotient \(\mathbb Z/N\mathbb Z\). Take a unit \(r\), put \(R=r^2\bmod N\), and use the ray condition
\[
 \beta\equiv R\pmod{\mathfrak n}.
\tag{21}
\]
Then \(r^2=\psi(\beta)\) is a valid one-sided row. If enough such \(\beta\)'s factor over a manageable algebraic base and their product can be certified square, the argument in §5 proves fruitfulness. A small rational counterpart is unnecessary for this redesigned row.

This is an algebraic construction, not a complexity theorem. Since \(N\mathfrak n=N\), the normalized radius bound acquires
\[
 r_{\rm box}^{\,n}\ \supseteq\ \Delta^{3/2}N.
\tag{22}
\]
The notation means that the specified radius formula contains this factor, not that every sample's norm is bounded below by it. The known smooth walk contribution must also be included.

If instead one samples an unconstrained \(\beta\) and tries to choose \(r\) afterwards, one has not produced a known random square root of \(\psi(\beta)\). Adding the representative of \(\psi(\beta)\) as a new rational column usually adds one new column per row. Factoring that representative can cost a call at size \(N\), with no fixed size contraction. Gauss-splitting it creates another smoothness/descent problem.

The known-\(\tau\) ray theorem is the appropriate statement for (21). One should not import the known-factorization assumption on the auxiliary modulus in Theorem 13.2 as a proof that every use of Theorem 9.5 is circular: the prescribed residue is already available here. Order/conductor compatibility, the norm scale and relation supply remain to be proved.

There is a useful cost diagnosis even under an optimistic smoothness law. Suppose the unknown bulk has the certified scale \(x=e^{ct+o(t)}\), while the final base is \(y=e^{b\Lambda_\alpha}\). A Dickman-type law would give
\[
 u=\frac{ct}{b\Lambda_\alpha}(1+o(1)),\qquad
 \log(1/p_{\rm smooth})\sim u\log u
 =\left(\frac{c(1-\alpha)}b+o(1)\right)
             t^{1-\alpha}(\log t)^\alpha.
\]
The smoothness cost is then at the \(L_N[1-\alpha]\) scale; balancing it with a base and linear algebra at \(L_N[\alpha]\) requires \(\alpha\ge1/2\). At a final \(L[1/3]\) base, this calculation gives \(L[2/3]\) smoothness work. Early levels with a much larger cutoff can be cheaper, but repeating the same discriminant-sized sampler does not lower the terminal scale.

This is **HEURISTIC cost analysis of the certified envelope**, not a complexity lower bound for dBPW's actual output distribution or for factoring. The positive theorem (7)–(8) instead gives an ERH one-step ideal reduction at its stated \(1/2\)-type discriminant scale. No smaller exponent follows from substituting NFS parameters into that theorem.

**(d) Several fields or a fresh field at each level: OPEN.**

Independent fields can help rational primes: a rational target can be searched in any field sharing the rational side. To amplify success, one needs conditional bounds given that shared side and the history selecting the target. A bound averaged jointly over \(g\) and fields is insufficient if failure is concentrated on some \(g\)'s.

An algebraic child cannot simply move to a fresh field. A prime ideal \(\mathfrak p\subset\mathcal O_{K_1}\) is not an ideal of \(K_2\). Replacing it by its rational norm \(p\) loses information, since
\[
 p\mathcal O_{K_1}=\prod_{\mathfrak P\mid p}
                              \mathfrak P^{e(\mathfrak P/p)}
\tag{23}
\]
contains the other primes above \(p\). Recovering the chosen prime from \(p\) needs relations for those other factors or a separate transfer construction.

Shared rational residues can glue explicit *elements* in different fields. They do not give canonical transfers of arbitrary ideal classes. A fresh-field strategy therefore needs a fully specified relation identity, a bound on all extra columns and fields, and a proof of its conditional success rate. None is supplied by the present hypotheses.

**(e) Jensen, Paley–Zygmund, whole trees and field randomness.**

**PROVED:** if \(g(A)\) is the density of targets with completed certificates, then
\[
 \mathbb E g(A)^{2k}\ge(\mathbb E g(A))^{2k}.
\tag{24}
\]
This is useful after (17) has been proved. It does not prove (17).

**PROVED:** for a nonnegative completed-tree count \(Z\),
\[
 \Pr(Z>0)\ge(\mathbb EZ)^2/\mathbb EZ^2.
\tag{25}
\]
Using it requires moments of completed trees; moments of isolated smooth values are different quantities.

**OPEN:** a suitable tree-moment estimate for the LV field law. A random degree-\(d\) coefficient vector has only \(d\) coordinate directions. Evaluation constraints at sufficiently many independent points determine that vector; one cannot multiply fresh-coefficient probabilities at every node of a much larger tree. This identifies a failure in a proposed independence proof, not an information-theoretic impossibility theorem. One random object can support many correlated successful searches.

**4. Costs, tree size, split counting and recursion**

**Uniform split counting — PROVED.**

Let \(N\) be odd, \(h=2^{\omega(N)}\), and \(\mathcal P\) a set of \(m\) primes coprime to \(N\). Suppose every product of \(k\) distinct primes of \(\mathcal P\) is at most \(X\), with \(2X^2<N\), and \(m\ge2k^2h\). Let \(\mathcal U\) be all these products. For \(R\) uniform on the squares in \((\mathbb Z/N\mathbb Z)^\times\), the Gauss split has both positive coordinates in \(\mathcal U\), up to simultaneous sign, with probability at least
\[
 \frac{|\mathcal U|^2}{2\phi(N)}.
\tag{26}
\]

*Proof.* Partition \(\mathcal U\) into the \(h\) square classes. At least \(|\mathcal U|^2/h\) ordered pairs have the same class. For each \(u\), at most \(k{m-1\choose k-1}=(k^2/m)|\mathcal U|\) partners share a prime with it. At least \(|\mathcal U|^2/(2h)\) pairs are therefore coprime and in the same class.

Their ratios are distinct modulo \(N\): equality implies \(uv'-u'v\) is a multiple of \(N\) of absolute value less than \(N\), hence zero; reduced positive fractions are unique. For such a pair, a shortest lattice vector \(s\) has norm at most \(\|(u,v)\|\). Its determinant with \((u,v)\) is a multiple of \(N\) of absolute value at most \(2X^2<N\), hence zero. Primitivity forces \(s=\pm(u,v)\). Divide the pair count by \(\phi(N)/h\). \(\square\)

This proof is uniform in \(k\) and avoids an exponential loss in \(\omega(N)\). It establishes the needed split-counting implication without accepting the internal proof label.

For the \(1/3\) size model, take \(k=\kappa D'+O(1)\), \(X=\sqrt N/3\), and \(Y=X^{1/k}/2\). For a good set of relative density \(g\ge t^{-C}\) in \((Y,2Y]\), prime counting and
\({m\choose k}\ge(m/k)^k\) give
\[
 \log|\mathcal U|
 \ge \tfrac12t-\kappa(C+1)E+o(E).
\]
Consequently,
\[
 \Pr(\text{acceptable split})\ge
       \exp[-2\kappa(C+1)E+o(E)].
\tag{27}
\]
Here \(\log\log Y+\log k=\log t+O(1)\), which includes the growing-\(k\) combinatorial loss.

After trial division to \(e^{bE}\), \(\omega(N)\le t/(bE)\). Meanwhile \(\log m\asymp Q\), which is larger than \(\log(2^{\omega(N)})=O(t/E)\) by a factor of order \(\log t\). Thus the hypothesis \(m\ge2k^2h\) holds for these good fields. This observation matters for inputs with many distinct prime factors.

**Multi-level norm arithmetic — HEURISTIC algorithmic conclusion, PROVED conditional bookkeeping.**

Take \(d=\delta D'\), \(B=e^{bE}\), and search inflation \(S=e^{sE}\). For a target with \(z=\log q\), the proposed degree-one geometry gives upper envelopes of the form
\[
 \log X_F\le\kappa_FQ+\frac d2z+dsE+o(t),\qquad
 \log X_r\le m_0Q+\frac12z+sE+o(Q+z).
\tag{28}
\]
For an algebraic target, division by its norm can improve the first envelope; not using that improvement is conservative.

If \(z\ge bE\) and all children are required to be at most \(e^{\lambda z}\), then
\[
 u_F+u_r\le
 \left(\frac{\kappa_F+m_0}{\lambda b}
       +\frac{\delta}{2\lambda}
       +\frac{\delta s}{\lambda b}+o(1)\right)D'.
\tag{29}
\]
Under the additional restricted-smoothness hypothesis, a polynomial prime-density loss contributes \(O(D'\log t)=O(E)\) to the logarithmic search cost. Also \(D'\log D'=E/3+o(E)\). Choosing \(b\) large makes the coefficient of \(s\) in such a search inequality smaller than the available \(2s\), leaving finite constants.

These statements explain why the size model is plausible. They do not supply the restricted-smoothness hypothesis.

With contraction \(z\mapsto\lambda z\), a top bound \(z=O(Q)\) reaches \(bE\) in
\[
 O(\log(Q/E))=O(\log t)
\]
levels. In a completed certificate, every nonbase child has norm exceeding \(B\), and the total logarithmic norm at a node is \(O(t)\). Hence there are at most \(O(t/\log B)\) such children, counting multiplicity, under (28). The entire completed tree has at most
\[
 \exp(O((\log t)^2))=\exp(o(E))
\tag{30}
\]
nodes. Base-factor multiplicities are stored in exponent vectors. This is a bound for a completed tree; it does not justify discarding failed search trees without charging their work. Capped searches or a uniform conditional expected-work bound are needed.

At the largest target, \(\log q\le b_1Q\), (28) permits algebraic testing inputs of size at most
\[
 N^{\tau_F+o(1)},\qquad \tau_F=\delta b_1/2.
\tag{31}
\]
Choose \(\delta b_1/2<1\). Other terms in (28) are \(o(t)\), including \(dsE=O(Q)\). Thus recursive complete factorization tests are strictly smaller than the original input. Factoring split halves has contraction \(1/2+o(1)\).

**Recursion lemma — PROVED, conditional on outer work bounds.**

Suppose an outer procedure on \(n\) makes at most
\(L_n[\alpha,a+o(1)]\) expected auxiliary calls on integers at most \(n^{\tau+o(1)}\), for fixed \(\tau<1\), and its remaining expected work is \(L_n[\alpha,p+o(1)]\). Suppose auxiliary outputs are canonical complete factorizations, its logical decisions do not depend on their running times, and auxiliary calls use fresh coins. Then its self-recursive implementation has expected cost \(L_n[\alpha,c+o(1)]\) for every
\[
 c>\max\{p,\ a/(1-\tau^\alpha)\},
\tag{32}
\]
provided the outer success and work bounds hold uniformly and finite inputs are handled.

*Proof.* The unique canonical factorization makes the outer logical transcript independent of the chosen correct factorization oracle. If \(F(x)\) is the worst expected factorization cost up to \(x\), conditional expectation at every adaptive call gives
\[
 F(n)\le e^{(p+o(1))\Lambda_\alpha(n)}
      +e^{(a+o(1))\Lambda_\alpha(n)}
                      F(n^{\tau+o(1)}).
\]
Polynomial factors also absorb the at most \(O(\log n)\) splits needed for complete factorization after a successful proper split. Since
\(\Lambda_\alpha(n^\tau)=(\tau^\alpha+o(1))\Lambda_\alpha(n)\),
strong induction closes with any strict margin in (32). For finite exceptional inputs, dovetail a certified trial-division algorithm. \(\square\)

For the proposed tree, take \(\tau>\max(1/2,\delta b_1/2)\), still below 1. If (17), the search budgets and growing-degree postprocessing were proved, (27), (30), and the row target \(e^{bE+o(E)}\) would give finite \(a,p\). Equation (32) would then give a finite \(L_N[1/3,c]\) constant. This is a conditional reduction, not a proof of (17).

**One-level intermediate targets — PROVED algebra of a conditional model; algorithms OPEN.**

For a growing-degree one-level design, the stated size constraints imply
\[
 \theta\gtrsim1/d,\qquad\theta\gtrsim d/k,\qquad
 \theta\lesssim D^{2\alpha-1}.
\]
Thus \(d\gtrsim D^{1-2\alpha}\), \(k\gtrsim D^{2(1-2\alpha)}\). If the good-prime loss is
\(\log(1/\gamma)=O(d^\xi\log t)\), split acceptance requires \(kd^\xi=O(D^\alpha)\). Therefore
\[
 (2+\xi)(1-2\alpha)\le\alpha,\qquad
 \alpha\ge\frac{2+\xi}{5+2\xi}.
\tag{33}
\]

| Assumed uniform density loss | Conditional one-level threshold |
|---|---:|
| \(t^{-O(1)}\) | \(2/5\) |
| \(t^{-O(d)}\), or comparable \(d^{O(d)}\) loss | \(3/7\) |
| \(t^{-O(d^2)}\) | \(4/9\) |

None of these rows supplies a rigorous algorithm: the required uniform special-\(q\) moments, field exceptions and postprocessing are unproved here. In the multi-level exponent arithmetic, \(\alpha>1/3\) allows
\(\lceil(1-2\alpha)/(3\alpha-1)\rceil\) levels; this gives 1, 2, and 6 at \(\alpha=0.4,0.38,0.35\). At \(1/3\), constants rather than logarithmic exponents must contract, as in (28)–(30).

Finally, a portfolio of fixed-degree bounds with constants tending to zero at exponent \(1/2\) would give \(\exp(o(\sqrt{t\log t}))\), assuming those individual algorithms were proved. It gives no fixed exponent below \(1/2\). For example, a saving factor \(1/\log\log t\) in the exponent still dominates every \(t^{1/2-\epsilon}\) scale.

**5. Transcript invariance through an adaptive multi-level descent**

**PROVED, with explicit conditions.**

Preprocess \(N\) to be odd and not a prime power. Let
\[
 U=(\mathbb Z/N\mathbb Z)^\times,\quad
 H=\ker(u\mapsto u^2),\quad |H|=2^{\omega(N)}.
\]
Sample independent uniform \(r_i\in U\) and reveal \(R_i=r_i^2\). Conditional on the \(R_i\), write
\[
 r_i=s_i(R_i)\epsilon_i,
\tag{34}
\]
where the sections \(s_i\) are used only in the proof and the \(\epsilon_i\) are independent uniform elements of \(H\).

Suppose every exposed quantity is a function of the \(R_i\) and coins independent of the \(r_i\): field choices, splits, canonical factorizations, special-\(q\) points, all failed and successful child searches, their order, cache contents, primality certificates, characters, linear algebra, square-root choices, and stopping rules. This condition allows dependence on previous exposed quantities. By induction over computational steps, the entire transcript \(\mathcal T\) is measurable with respect to \(\sigma((R_i),C)\). It therefore leaves the conditional joint law of the \(\epsilon_i\) unchanged.

Suppose the algorithm then selects a nonempty row subset \(I=I(\mathcal T)\) and a unit \(Y=Y(\mathcal T)\) such that
\[
 Y^2=\prod_{i\in I}R_i\pmod N.
\]
Put \(X=\prod_{i\in I}r_i\). Choose one \(i_0\in I\), and condition also on all \(\epsilon_i\) except \(\epsilon_{i_0}\). Multiplication by that remaining uniform kernel element makes \(X/Y\) uniform on \(H\). Thus
\[
 \Pr(1<\gcd(X-Y,N)<N\mid\mathcal T)
      =1-\frac2{|H|}
      =1-2^{1-\omega(N)}\ge\frac12.
\tag{35}
\]
For odd prime powers, the only unit square roots of 1 are \(\pm1\), so the same proof includes nonsquarefree \(N\). A nonempty coefficient vector modulo 2 is what matters if rows are combined with general integer exponents.

The hidden roots may be multiplied to form \(X\) only after the dependency and \(Y\) have been selected. After a trivial final gcd, discard the rows and use fresh roots; conditioning repeatedly on failure while reusing the same hidden signs would require a separate argument.

Recursive tests do not compromise this proof. Their inputs depend on the outer \(R_i\), and their fresh random tapes are independent of the outer hidden signs. Inner executions may have their own hidden-root arguments. A cache keyed by field, ideal and exposed history remains invariant.

This establishes (T) as a structural property of a correctly designed multi-level algorithm. It does not establish that the required certificates or square roots can be produced within the claimed time.

**6. ERH: uses, possible weakening, and field averages**

In dBPW, the random-walk input is ERH for Hecke \(L\)-functions over the fixed field with the relevant modulus, including the characters controlling the archimedean directions. The mechanism is a bound on prime-ideal character sums and therefore on the nontrivial eigenvalues of the Arakelov walk. The geometric Haar-average result (4) and exact lattice arithmetic are separate. Dedekind-zeta RH also enters the small-prime product estimates used to control the residue contribution and smooth-ideal counts. See the authors' explicit discussion in [dBPW §21.3, pp. 112–113](https://arxiv.org/pdf/2512.01588v2#page=112).

**PROVED logical distinction:** Dirichlet GRH alone is not the full hypothesis of the general Hecke-walk theorem. Nor does that theorem assert pointwise restricted smoothness of a binary form in a special-\(q\) box.

For the present task, effective Chebotarev or Hecke prime-ideal counting can supply enough degree-one primes and separating characters. It does not count simultaneous products of prescribed field-dependent prime sets in two linked values. Likewise, a square-root-error estimate for ordinary prime counts does not itself give a positive correlation estimate for those two restricted smoothness events.

**Conditional replacement — PROVED implication, OPEN hypothesis.** It would suffice to replace ERH by quantitative mixing and prime-count estimates of the required strength for the actual fields and moduli used. For an averaged argument, one needs bounds that survive the algorithm's selection and use of fields. For example, if \(\mathcal B\) is the bad-field set and \(Z\) is a relation or completed-tree count, the relevant comparison can be
\[
 \mathbb E[Z\mathbf1_{\mathcal B}]
                  =o(\mathbb EZ),
\tag{36}
\]
not merely \(\Pr(\mathcal B)=o(1)\).

The LV known-root distribution is an \(N\)-dependent constrained family of growing-degree fields. No average spectral theorem or restricted two-sided smoothness theorem sufficient for (17) has been established here for that family. Averages over unrelated field families do not automatically apply.

In a size-biased field-selection argument, even unconditional marginal field averages must be weighted correctly. If \(q_A\) is a field's row success probability, retaining the first successful field gives density proportional to \(q_A\), and
\[
 \mathbb E_{\rm retained}(1/q_A)
       =\Pr(q_A>0)/\mathbb E_Aq_A.
\tag{37}
\]
This follows by integrating \(q_A(1/q_A)\). It is useful bookkeeping, but still needs a lower bound on \(\mathbb E_Aq_A\). With factor exits, include those as absorbing successes rather than condition away their contribution.

**7. Literature and what its theorems do not imply**

**General integer factoring.** Lenstra and Pomerance, “A rigorous time bound for factoring integers,” *J. Amer. Math. Soc.* **5** (1992), 483–516, prove expected complete factorization time \(L_N[1/2,1+o(1)]\). This is unconditional. [Original paper, author-hosted PDF](https://math.dartmouth.edu/~carlp/PDF/paper85.pdf).

Their Theorem 6.1, p. 499, concerns cumulative counts \(\Psi(x,y;\mathcal P)\) for arbitrary prime sets. It requires both reciprocal prime mass near \(y\) and a lower prime count at a smaller scale. In its notation,
\[
\begin{gathered}
 0\le\eta\le1,\quad
 2\le y\le\exp((\log x)^{1/2}(\log\log x)^\eta),\\
 u=\log x/\log y,\quad
 v=y^{1-1/\log u},\quad w=v^{(\log u)^{-\eta}},\\
 \sum_{\substack{p\in\mathcal P\\v<p\le y}}\frac1p
       \ge\frac1{A\log u},\qquad
 \pi(w;\mathcal P)\ge\frac{w}{B_1\log w},\\
 \Psi(x,y;\mathcal P)\ge x\exp\!\left[
 -u\{\log u+12(\log u)^\eta+\log\log u+
                2(\log u)^{\eta-1}\log B_1+\log A\}\right].
\end{gathered}
\tag{38}
\]
There is also an effective lower threshold on \(x\), and \(A,B_1\ge1\). This is a cumulative integer count, not a count in each short interval, progression or polynomial-value slice. The displayed range also prevents direct use with \(x=e^{\Theta(t)}\), \(y=e^{\Theta(Q)}\). These distinctions are independent obstacles to substituting it into the proposed descent. [Theorem 6.1, original PDF](https://math.dartmouth.edu/~carlp/PDF/paper85.pdf#page=17).

**Randomized NFS.** Jonathan D. Lee and Ramarathnam Venkatesan, “Rigorous Analysis of a Randomised Number Field Sieve,” *J. Number Theory* **187** (2018), 92–159, DOI 10.1016/j.jnt.2017.10.019. Their Theorem 2.1 produces congruences of squares in expected \(L_N[1/3,(64/9)^{1/3}+o(1)]\). Theorem 2.3 adds a fruitful-factor conclusion for semiprimes with both factors \(3\bmod4\), conditional on their Conjecture 7.1. That conjecture is an additional character assertion, not identified with ERH. Consequently this paper is not an ERH-only factoring theorem for every \(N\). [Paper, Theorems 2.1–2.3 and §7](https://arxiv.org/pdf/1805.08873).

**Finite-field discrete logarithms.** There is an important positive exception to a blanket “no rigorous sub-\(1/2\)” claim. Thorsten Kleinjung and Benjamin Wesolowski, “Discrete logarithms in quasi-polynomial time in finite fields of fixed characteristic,” *J. Amer. Math. Soc.* **35** (2022), 581–624, DOI 10.1090/jams/985, prove expected time
\[
 (pn)^{\,2\log_2 n+O(1)}
\tag{39}
\]
in \(\mathbb F_{p^n}\). For every fixed characteristic this is quasipolynomial, hence asymptotically better than \(L[1/3,c]\). It is unconditional. It is not a uniform \(L[1/3]\) theorem for all characteristics growing with field size; for \(n=1\), (39) is merely polynomial in \(p\). [Theorem 1.1 and published reference](https://doi.org/10.1090/jams/985).

Large-characteristic NFS-DL \(L[1/3]\) analyses and the usual function-field-sieve smoothness analyses must be distinguished from this theorem. I found no ERH-only \(L[1/3]\) result covering arbitrary prime fields. The fixed-characteristic result exploits special finite-field elimination geometry; it is not a number-field arbitrary-good-set lemma.

**Number-field class groups and units.** Jean-François Biasse's “An \(L(1/3)\) algorithm for ideal class group and regulator computation in certain number fields,” arXiv:0912.1927, and “An \(L(1/3)\) algorithm for discrete logarithm computation and principality testing in certain number fields,” arXiv:1204.1292, concern restricted defining-polynomial families and use heuristics. In particular, the latter explicitly assumes both relation generation and smoothness heuristics. [First paper](https://arxiv.org/abs/0912.1927), [second paper, Heuristics 1–2](https://arxiv.org/pdf/1204.1292).

Jean-François Biasse and Claus Fieker, “Subexponential class group and unit group computation in large degree number fields,” *LMS J. Comput. Math.* **17**, Special Issue A (2014), 385–403, DOI 10.1112/S1461157014000345, explicitly assumes GRH **and additional smoothness heuristics**. Its general subexponential result and its restricted \(L[1/3]\) regimes are not ERH-only theorems. [Original article](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/4387ACB036E3358143A563F196E386CB/S1461157014000345a.pdf).

dBPW's 2025 preprint, revised February 2026, supplies a genuine ERH-only general class/unit-group result with input integral basis:
\[
 \operatorname{poly}\left(
 \text{input length},L_\Delta[1/2,O(1)],
 L_{n^n}[2/3,O(1)],
 \min(\rho_K,L_\Delta[2/3+o(1),O(1)])
 \right).
\tag{40}
\]
This is Theorem 1.3, not a general \(L_\Delta[1/3]\) theorem. Its ideal-class logarithm and principality applications use the decomposition theorem (8). [dBPW, Theorem 1.3](https://arxiv.org/pdf/2512.01588v2#page=4).

For the NFS degree \(n\asymp D'\), the degree term in (40) is \(e^{o(E)}\); the standard residue upper bound \(\rho_K\le(\log\Delta)^{O(n)}\), with \(\log\Delta=O(t)\), contributes at most \(e^{O(E)}\). This residue upper bound itself does not require ERH. The \(L_\Delta[1/2]\) term remains at exponent \(1/2\) when \(\log\Delta\asymp t\). This is a direct parameter substitution, not a lower bound on every implementation.

**Jacobians.** Andreas Enge, Pierrick Gaudry and Emmanuel Thomé, “An \(L(1/3)\) Discrete Logarithm Algorithm for Low Degree Curves,” arXiv:0905.2177 (2009), explicitly uses Heuristic 5 for smoothness of restricted principal divisors, including the special-place descent, and Heuristic 6 for relation generation. The \(L[1/3]\) claim is therefore HEURISTIC. The earlier Enge–Gaudry paper is in *EUROCRYPT 2007*, LNCS **4515**, 379–393, with the \(L[1/3+\epsilon]\) formulation. [EGT, §3 and Heuristics 5–6](https://arxiv.org/pdf/0905.2177), [earlier paper](https://www.iacr.org/archive/eurocrypt2007/45150379/45150379.pdf).

Claus Diem's “An index calculus algorithm for non-singular plane curves of high genus,” ECC 2006, also labels the relevant relation and descent statements heuristic. It should not be cited as a rigorous \(L[1/3]\) theorem. [Author's slides, Heuristic Results 1–2](https://www.math.uni-leipzig.de/~diem/preprints/ecc2006.pdf).

Diem also has a separate rigorous result: “On the discrete logarithm problem in class groups of curves,” *Math. Comp.* **80** (2011), 443–475, DOI 10.1090/S0025-5718-2010-02281-1, gives expected \(\widetilde O(q^{2-2/g})\) for every fixed \(g\ge2\). This is a fixed-genus power bound, not a uniform growing-genus \(L[1/3]\) theorem. Constants depending on fixed genus cannot be discarded by allowing genus to grow. [Original journal abstract](https://www.ams.org/journals/mcom/2011-80-273/S0025-5718-2010-02281-1/home.html).

There are rigorous subexponential Jacobian algorithms at the \(1/2\) scale or \(1/2+\epsilon\) formulations under their curve and presentation hypotheses: Andreas Enge, “Computing discrete logarithms in high-genus hyperelliptic Jacobians in provably subexponential time,” *Math. Comp.* **71** (2002), 729–742, DOI 10.1090/S0025-5718-01-01363-1; Jean-Marc Couveignes, “Algebraic groups and discrete logarithm,” in *Public-Key Cryptography and Computational Number Theory* (2001), 17–28, DOI 10.1515/9783110881035.17; Florian Hess, “Computing relations in divisor class groups of algebraic curves over finite fields,” preprint, cited as 2004 in EGT.

Couveignes' Theorem 1 assumes a supplied rational origin and a polynomial-size plane presentation; it is a subexponential theorem, not an \(L[1/3]\) theorem. [Original chapter preview, pp. 19–20](https://api.pageplace.de/preview/DT0400.9783110881035_A19986600/preview-9783110881035_A19986600.pdf#page=31). Enge's author presentation distinguishes the proven \(L[1/2]\) work from the heuristic \(L[1/3]\) construction. [Enge, ECC 2006](https://cacr.uwaterloo.ca/conferences/2006/ecc2006/enge.pdf).

The Hess preprint's historical URLs were not retrievable in this session; its citation is recorded in EGT's bibliography. I do not certify exact constants or all hypotheses from that unavailable text. No inference about factoring relies on it.

**2024–2026 check.** The material development relevant here is dBPW's 2025/2026 rigorous ideal sampler and class/unit-group analysis. Andreas-Stephan Elsenhans and John Voight, “Computing class groups and unit groups in Magma,” arXiv:2510.05501 (2025), carefully distinguishes rigorous correctness, GRH-conditional correctness and heuristic regimes; it is not a new general \(L[1/3]\) runtime theorem. [Original preprint](https://arxiv.org/abs/2510.05501).

Searches for rigorous or ERH-only general \(L[1/3]\) factoring did not identify such a theorem as of this report. This is a literature-search finding, not a proof that no unpublished or overlooked result exists. The affirmative finite-field result (39) and the general ideal result (40) should both be acknowledged without transferring their hypotheses or conclusions to integer factoring.

**8. Final claim ledger**

The following are the claims I consider **PROVED**, with qualifications preserved:

1. The stated dBPW theorem implications hold under their external ERH and representation hypotheses. Their density is (1), their lattice is full-dimensional, and their explicit cofactor bound is (10); these do not assert rank-two two-sided sampling.
2. The adaptive one-sided lemma follows from the pointwise ideal-sampling theorem when the local density hypothesis (15) holds.
3. Equation (16) is an unconditional combinatorial lower bound for the monoid generated by an arbitrary sufficiently numerous set of prime ideals.
4. Class compatibility and (19) show why a bare ambient density argument cannot establish a universal slice-hitting theorem.
5. The uniform split probability (26) holds, including the growing-\(k\) estimate (27) when the stated good-prime density is available.
6. Geometric contraction gives the completed-tree size (30), and the recursion implication (32) closes under its uniform work and success hypotheses.
7. Transcript invariance (T), and the fruitful-gcd probability (35), hold for adaptive multi-level trees satisfying the specified exposure rule.
8. Equation (33) gives the one-level thresholds \(2/5,3/7,4/9\) within the stated conditional size and density-loss model. It does not establish algorithms at those exponents.
9. Lenstra–Pomerance provides unconditional complete factoring in expected \(L_N[1/2,1+o(1)]\). Kleinjung–Wesolowski provides the separate unconditional fixed-characteristic finite-field result (39).

**HEURISTIC:** \(L[1/3]\) multi-level factoring in the proposed size model, the needed restricted smoothness law for the endogenous good sets, and typical discriminant-size assertions not proved for the specified field distribution.

**OPEN:** (17), a usable substitute for it, a fresh-field transfer that handles algebraic children at the required cost, and the necessary uniform growing-degree special-\(q\) moment and postprocessing estimates. No fixed \(\alpha<1/2\) factoring bound is established here, even under ERH.

**Best supported general factoring exponent: \(\boxed{\alpha=1/2}\), with unconditional constant \(1\) from Lenstra–Pomerance.**

**Single sharpest barrier:** prove a uniform quantitative supply of short two-sided special-\(q\) relations whose factors lie in the recursively certified good sets of the same field. ERH-controlled full-dimensional Arakelov sampling does not currently provide this supply at contracting NFS norm scales.
