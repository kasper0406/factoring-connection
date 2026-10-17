> Source: GPT-6.1-sol (codex exec, web search, high effort), round 1, 2026-10-17. Prompt: codex/r1_tal_bgs.md. Raw: codex/r1_tal_bgs.out. Curated copy; claims not yet independently re-verified by orchestrator unless noted.

**Tal’s exponent saving does not directly transfer to an unfactored RSA modulus.** The saving comes from reconstructing prime-specific factorial ratios; its general interval-product primitive still has square-root cost. BGS already works over \(\mathbb Z/N\mathbb Z\), including factor discovery through failed inversions. Harvey and Harvey–Hittmeir go further by batching collisions modulo hidden factors, giving a deterministic \(N^{1/5+o(1)}\) algorithm for your predicate through factoring.

I read the full papers, including the technical appendices. The versions used are:

| Reference | Version consulted |
|---|---|
| Tal, *Computing Modular Factorials Below the Square-Root Barrier* | September 2026, ECCC TR26-211, 69 pages; [PDF](https://eccc.weizmann.ac.il/report/2026/211/download) |
| Bostan–Gaudry–Schost, *Linear Recurrences with Polynomial Coefficients and Application to Integer Factorization and Cartier–Manin Operator* | SIAM J. Comput. 36(6), 2007, pp. 1777–1806; [PDF](https://mathexp.eu/bostan/publications/BoGaSc07.pdf) |
| Harvey, *An Exponent One-Fifth Algorithm for Deterministic Integer Factorisation* | arXiv:2010.05450v1, 12 October 2020; [PDF](https://arxiv.org/pdf/2010.05450v1) |
| Harvey–Hittmeir, *A Log-Log Speedup for Exponent One-Fifth Deterministic Integer Factorisation* | arXiv:2105.11105v1, 24 May 2021; [PDF](https://arxiv.org/pdf/2105.11105v1) |

Below, write \(N=P Q\), with unknown primes \(P<Q\). This avoids confusing the RSA factors with Tal’s **known prime \(p\)** and **auxiliary conductor \(q\)**. Put \(n=K-1\), \(\ell=\lceil\log_2N\rceil\), and let \(\mathsf M(b)\) denote integer-multiplication cost.

---

**1. What is already available for \(U(K)\)**

For a semiprime,

\[
U(K)\iff P\ge K.
\]

Consequently, a deterministic bound immediately available from these papers is

\[
T_U(N,K)\le
\min\left\{
T_{\mathrm{BGS}}(N,K-1),\
O\!\left(
\frac{N^{1/5}\ell^{16/5}}{(\log\ell)^{3/5}}
\right)
\right\},
\]

where

\[
T_{\mathrm{BGS}}(N,b)=
O\!\left(
\mathsf M(\sqrt b\,\ell)
+\log b\,\mathsf M(\ell)\log\ell
\right).
\]

The first branch tests for a prime divisor at most \(K-1\); it need not compute the entire factorial residue. The second factors \(N\) and compares its smaller factor with \(K\). The cited results are BGS Lemma 13 and HH Theorem 1.1. [BGS, §5](https://mathexp.eu/bostan/publications/BoGaSc07.pdf), [HH](https://arxiv.org/pdf/2105.11105v1).

Ignoring logarithms, HH wins over the square-root-in-\(K\) branch when \(K>N^{2/5}\). Near the interesting threshold for balanced RSA, \(K\asymp\sqrt N\), the comparison is \(N^{1/5}\) versus \(N^{1/4}\).

This does **not** establish that unit testing is as difficult as factorial evaluation. It establishes that a uniformly fast exact unit predicate is a factoring primitive: binary search needs \(O(\log N)\) calls.

One interval distinction matters throughout:

\[
\gcd\!\left(\prod_{A\le j<B}j,N\right)>1
\]

means that the interval contains a **multiple of a prime factor**, not necessarily a divisor of \(N\). For example, \([14,15)\) detects a nonunit modulo \(77\), although \(14\nmid77\).

---

**2. Tal: inputs, identities, reconstruction, and costs**

Tal’s bounds are

\[
\widetilde O\!\left(q^c+\sqrt{p/q}\right)
\]

for \(q\mid p-1\), modulo \(p\) or \(p^2\), and

\[
\widetilde O_r\!\left(
k^{O_r(1)}
\left[
q^{c_r}+
\left(\frac{\Phi_r(p)}q\right)^{1/(2\varphi(r))}
\right]
\right)
\]

for fixed odd prime-power \(r\), \(q\mid\Phi_r(p)\), modulo \(p^k\). These are Theorems 1.1–1.3. **The paper does not give numerical values for \(c,c_r\), or the precision-overhead exponent**; making them explicit is an open problem in §8. Thus a more explicit per-step polynomial exponent cannot honestly be extracted from this paper. [Tal](https://eccc.weizmann.ac.il/report/2026/211/download).

**2.1. What the setup actually requires — §§1.3, 2.1–2.2**

For each conductor \(m\), define

\[
L_m=\mathbb Q(\zeta_m),\qquad
\mathcal O_m=\mathbb Z[\zeta_m],\qquad
f_m=\operatorname{ord}_m(p).
\]

Choose \(\eta_m\in\mathbb F_{p^{f_m}}^\times\) of exact order \(m\). This specifies

\[
\rho_{\mathfrak p_m}:\mathcal O_m\longrightarrow
\mathbb F_{p^{f_m}},
\qquad
\zeta_m\longmapsto\eta_m,
\]

and the prime ideal \(\mathfrak p_m=\ker\rho_{\mathfrak p_m}\).

At divisors \(m\mid q\), choices must be compatible:

\[
\zeta_m=\zeta_q^{q/m},\qquad
\eta_m=\eta_q^{q/m},\qquad
\chi_m=\chi_q^{q/m}.
\]

The character convention is

\[
\rho_{\mathfrak p_m}(\chi_m(x))
=x^{-(p^f-1)/m},
\qquad x\ne0,
\]

with every character power extended by zero at zero.

| Construction | Required arithmetic data | What is **not** required |
|---|---|---|
| Split, \(q\mid p-1\) | Known \(p,q\); an exact-order-\(q\) element of \(\mathbb F_p^\times\); compatible character/reduction map | A generator of all \(\mathbb F_p^\times\); the full factorization of \(p-1\) when \(q\) is supplied |
| Higher cyclotomic | Known \(p,q,r\); a representation of \(\mathbb F_{p^r}\); exact-order-\(q\) element; compatible maps | Full factorization of \(p^r-1\); a generator of its multiplicative group |
| \(p+1\) construction | Quadratic extension \(\mathbb F_{p^2}\); character of order \(2q\) | A full multiplicative-group generator |

Factoring **\(q\)** suffices to certify the chosen element’s order:

\[
\eta^{q/s}\ne1
\quad\text{for every prime }s\mid q.
\]

With \(q\) supplied, its factorization costs \(q^{O(1)}\), already absorbed into reconstruction. Random sampling followed by powering has expected \(q/\varphi(q)\) trials.

Factoring **\(p-1\)** or **\(\Phi_r(p)\)** is used to *select a useful \(q\)*, rather than as an additional prerequisite once \(q\) is supplied. Proposition 1.6 gives approximate divisor selection from a known factorization. For fixed \(r\), randomized subexponential factoring of \(\Phi_r(p)\) costs \(p^{o(1)}\). The split deterministic setup uses a \(p^{1/4+o(1)}\) search for suitable power nonresidues, absorbed by its stated bound.

The character can be represented through this order-\(q\) subgroup; no discrete logarithm in the entire field multiplicative group is required.

**2.2. Gauss–Jacobi identities — Proposition 2.1, Appendix B**

For \(Q=p^f\),

\[
G_Q(\theta)=\sum_{x\in\mathbb F_Q}\theta(x)\psi_Q(x),
\]

\[
J_Q(\theta_1,\ldots,\theta_s)
=\sum_{x_1+\cdots+x_s=1}\prod_i\theta_i(x_i).
\]

When the characters and their product are nontrivial,

\[
J_Q(\theta_1,\ldots,\theta_s)
=\frac{\prod_iG_Q(\theta_i)}
{G_Q(\prod_i\theta_i)},
\qquad
G_Q(\theta^p)=G_Q(\theta).
\]

For a nondegenerate binary Jacobi sum,

\[
J_Q(\chi^u,\chi^v)\,
\overline{J_Q(\chi^u,\chi^v)}
=Q.
\]

The degenerate values are explicit:

\[
J_Q(\chi^u,\chi^v)=
\begin{cases}
Q-2,&m\mid u,\ m\mid v,\\
-1,&\text{exactly one of }u,v\text{ is }0\bmod m,\\
-\chi(-1)^u,&u+v=0\bmod m,\ m\nmid u.
\end{cases}
\]

In the split case, with \(J_a=J_p(\chi,\ldots,\chi)\),

\[
J_a=\frac{G_p(\chi)^a}{G_p(\chi^a)},\qquad
J_a\overline{J_a}=p^{a-1},
\]

\[
J_{u+v}=J_uJ_vJ_p(\chi^u,\chi^v).
\]

**Algorithmic use:** an addition chain constructs \(J_a\) from \(O(\log a)\) binary sums. These formulas are not implemented by summing \(p\) or \(p^f\) terms.

**Requirements:** the finite-field character identities, Frobenius automorphism, and nontriviality hypotheses. Their proofs use finite-field changes of variables and character orthogonality.

**2.3. Reconstructing a Jacobi sum without enumerating the field — §§2.4–2.8, §3**

For representatives \(s\in(\mathbb Z/m\mathbb Z)^\times\), put

\[
c_s(u,v)=
\left\lfloor\frac{(u+v)s}{m}\right\rfloor
-\left\lfloor\frac{us}{m}\right\rfloor
-\left\lfloor\frac{vs}{m}\right\rfloor.
\]

For \(H_m=\langle p\rangle\), the orbit exponents are

\[
e_t(u,v)=
\sum_{j=0}^{f_m-1}c_{[tp^j]_m}(u,v).
\]

Stickelberger gives the **prime-ideal factorization**

\[
\bigl(J_{p^{f_m},m}(u,v)\bigr)
=
\prod_{t\in(\mathbb Z/m\mathbb Z)^\times/H_m}
\mathfrak p_{m,t}^{\,e_t(u,v)}.
\tag{Tal (15)}
\]

In the split case \(f_m=1\). For the multifold sum, telescoping yields

\[
(J_a)=
\prod_{t\in(\mathbb Z/q\mathbb Z)^\times}
\mathfrak p_{q,t}^{\lfloor at/q\rfloor}.
\tag{Proposition 2.4; Appendix C}
\]

The reconstruction pipeline is:

| Step | Operation and cost | Exact information required |
|---|---|---|
| Build cyclotomic arithmetic | Construct \(\Phi_m\), basis, multiplication and conjugation; \(m^{O(1)}\) bit operations | Integer \(m\) |
| Construct ideals | Form the above ideals in integer-basis/Hermite form; polynomial in \(m,\log p\) | Known \(p\); the selected factors of \(\Phi_m\bmod p\), determined by \(\eta_m\); Frobenius orbits |
| Recover polarized generator | Given \(I\overline I=w\mathcal O_m\), recover \(\beta\) with \((\beta)=I,\ \beta\overline\beta=w\); polynomial in input size | The actual ideal \(I\), its polarization \(w=p^f\), and CM-order structure |
| Orient the root of unity | Test at most \(2m\) candidates using proper-conductor sums; polynomial in \(m,f,\log p\) | Correct sums at every proper divisor of \(m\) |
| Assemble multifold sum | Addition chain; \(O(\log a)\) binary factors, stored with binary exponents | Previously reconstructed binary sums |
| Evaluate locally | Polynomial in \(m,k,\log p\) and representation size | Chosen prime above known \(p\), Hensel lifts, and known local valuations |

The generator-recovery step is Tal’s Theorem 2.7, invoking Lenstra–Silverberg. The ambiguity is only a root of unity because two polarized generators have quotient of absolute value one at every embedding.

The orientation test is explicit, not an unspecified sign choice. In Lemma 2.8, polynomial CRT combines each candidate with proper-conductor components into

\[
A(X)=\sum_{j=0}^{m-1}a_jX^j.
\]

The correct candidate has integral coefficients and satisfies

\[
A'(1)\equiv
\begin{cases}
0\pmod m,&m\text{ odd},\\[2mm]
\dfrac{Q-1}{2}\bigl(\gcd(m,u)+\gcd(m,v)\bigr)\pmod m,
&m\text{ even}.
\end{cases}
\]

Processing conductors in increasing order supplies all proper-conductor inputs.

Summing polynomial costs over \(m\mid q\), and including addition-chain factors, gives \(\widetilde O(q^c)\) in the split case and \(\widetilde O_r(q^{c_r})\) in the higher case (§3.3).

**Local divisions deserve attention.** Lemma 2.5 handles a factor \(A/b\) by finding \(v_p(b)\), increasing precision, dividing out the matching power of \(p\), and then inverting the remaining unit. A normalized Gauss sum of valuation \(e\) requires reduction to precision \(p^{e+k}\) before removing \(p^e\). This uses a **known uniformizer \(p\)** and a **specified local prime**, not merely a ring inversion.

**2.4. Split factorial decoder — Proposition 4.1, §4.2, Appendix A**

Let

\[
B=\frac{p-1}{q}.
\]

The central identity is

\[
(aB)!
\equiv
(-1)^{a-1}\rho_{\mathfrak p}(J_a)(B!)^a
\pmod p,
\qquad 1\le a<q.
\tag{Tal (2)}
\]

Algorithm:

1. If necessary, replace \(n\) by \(m=p-1-n\) using
   \[
   n!\equiv(-1)^{m+1}(m!)^{-1}\pmod p.
   \]
2. Choose the nearest multiple \(aB\) to the lower-half target.
3. Compute \(B!\) using the interval-product routine.
4. Reconstruct \(J_a\), then use (2).
5. Multiply or divide by the intervening interval, of length at most \(B/2\).

Costs:

\[
B!:\widetilde O(\sqrt B),\qquad
J_a:\widetilde O(q^c),\qquad
\text{correction}:\widetilde O(\sqrt B),
\]

plus logarithmic-cost powering and inversion.

The proof in Appendix A uses

\[
1_{y=0}=1-y^{p-1}
\]

and

\[
\sum_{x\in\mathbb F_p}x^d=
\begin{cases}
-1,&d>0,\ p-1\mid d,\\
0,&p-1\nmid d,
\end{cases}
\]

so that exactly one multinomial term survives. This is one precise place where the **field characteristic and group order** enter.

**2.5. Higher cyclotomic decoder — §§3.2, 5; Appendices B–D**

Let

\[
r=s^e,\qquad h=s^{e-1},\qquad
d=\varphi(r)=(s-1)h,\qquad
D=\frac{\Phi_r(p)}q.
\]

Lemmas 2.2–2.3 remove the at-most-one factor \(s\) from \(q\), obtaining, for every \(m>1\) dividing the remaining conductor,

\[
\gcd(m,p^h-1)=1,\qquad \operatorname{ord}_m(p)=r.
\]

Define

\[
\Theta_{r,m}=J_{p^r}(\chi_m,\chi_m^{p^h-1}).
\]

Gauss–Jacobi and Frobenius invariance give

\[
\sigma_u(\Theta_{r,m})
=G_{p^r}(\chi_m^{(p^h-1)u}),
\qquad
\Theta_{r,m}\overline{\Theta_{r,m}}=p^r.
\tag{Tal (12)}
\]

Thus one reconstructed family \((\Theta_{r,m})_{m\mid q}\) supplies all recursive Gauss sums. For a requested \(a\),

\[
g=\gcd(a,q),\quad m=q/g,\quad a'=a/g,\quad
u\equiv a'(p^h-1)^{-1}\pmod m,
\]

\[
G_{p^r}(\chi_q^a)=\sigma_u(\Theta_{r,m}).
\tag{Tal (17)}
\]

The lower-conductor exponent is unchanged:

\[
a\frac{p^r-1}{q}=a'\frac{p^r-1}{m}.
\]

The factorial identity follows from Gross–Koblitz. With Teichmüller character \(\omega\), \(\pi^{p-1}=-p\), and base-\(p\) digit sum \(s_p(a)\),

\[
G_{p^f}(\omega^{-a})
=-\pi^{s_p(a)}
\prod_{i=0}^{f-1}
\Gamma_p\!\left(\left\langle\frac{ap^i}{p^f-1}\right\rangle\right),
\tag{Tal (13)}
\]

whose first-order consequence is

\[
G_{p^f}(\omega^{-a})
\equiv
-\frac{\pi^{s_p(a)}}{\prod_i a_i!}
\pmod{\pi^{s_p(a)+p-1}}.
\tag{Tal (14)}
\]

Choose digits satisfying

\[
aD=\sum_{i=0}^{d-1}x_ip^i,\qquad
1\le x_u\le x_{u+h}\le\cdots\le x_{u+(s-2)h}<p.
\]

Set

\[
y_0=x_0-1,\quad y_u=x_u\ (0<u<h),\quad
y_i=x_i-x_{i-h}\ (i\ge h),
\]

\[
Y_0=x_{d-h}-1,\qquad
Y_u=x_{d-h+u}\ (u>0),\qquad
Y=\sum_uY_u=\sum_i y_i.
\]

Appendix D proves that the digits of \(a(p^r-1)/q\) are

\[
p-1-y_0,\ldots,p-1-y_{d-1},Y_0,\ldots,Y_{h-1}.
\]

Their sum is \(d(p-1)\). Consequently the Gauss sum has known valuation \(d\). Define its normalized unit by

\[
G_{p^r}(\chi_q^a)\equiv-p^dU_a
\pmod{\mathfrak p_q^{d+1}}.
\]

Proposition 5.1 gives

\[
\prod_{u=0}^{h-1}
\binom{Y_u}{y_u,y_{u+h},\ldots,y_{u+(s-2)h}}
=(-1)^Y U_a^{-1}
\quad\text{in }\mathbb F_p.
\]

Wilson’s identity

\[
(p-1-y)!\,y!\equiv(-1)^{y+1}\pmod p
\]

converts the complementary digit factorials into this expression.

**Lattice step, Lemma 5.2.** The lattice

\[
\Lambda_D=
\left\{x\in\mathbb Z^d:
\sum_i x_ip^i=0\pmod D
\right\}
\]

has determinant \(D\). Multiplication by \(\zeta_r\) preserves it. A short vector \(v\), followed by

\[
v,Tv,\ldots,T^{d-1}v,
\]

gives independent directions of length \(O_r(D^{1/d})\). Fixed-dimensional lattice reduction and rounding cost \(\operatorname{poly}_r(\log p)\).

Put \(E=\lceil C_rD^{1/d}\rceil\). The target makes one total \(Y_0=n+O_r(E)\), its \(s-1\) parts \(n/(s-1)+O_r(E)\), and every other chain \(O_r(E)\).

**Recursion, Theorem 5.3.** Compute one part factorial recursively; obtain the other large part factorials by short ratios; compute all small factorials directly; use Proposition 5.1; correct from \(Y_0!\) to \(n!\). There is one child of size below \(3n/4\), and all other products have length \(O_r(E)\).

Therefore:

\[
\text{per level}:
\widetilde O_r(E^{1/2})
+\text{polynomial local evaluation},
\]

\[
\text{total}:
\widetilde O_r\!\left(q^{c_r}+D^{1/(2d)}\right).
\]

For \(r=3\), this becomes particularly concrete:

\[
aD=pA+B,\qquad
\binom{A-1}{B-1}
\]

is reconstructed, with \(A-1\approx n\) and \(B-1\approx n/2\). Only \((B-1)!\) is recursive. The correction lengths are \(O(\sqrt D)\), hence their costs are \(\widetilde O(D^{1/4})\).

**Requirements:** known \(p\) controls the digits, the Gauss exponent, the valuation, Wilson complementation, and all factorial nonvanishing conditions. The lattice rounding itself is an integer algorithm and requires no finite-field division.

**2.6. Divisors of \(p+1\) — Proposition 4.2, Appendix G**

Assume \(q\ge4\), \(q\mid p+1\), and \(B=(p+1)/q\) odd. Choose an order-\(2q\) character on \(\mathbb F_{p^2}\); its restriction to \(\mathbb F_p^\times\) is quadratic.

For odd \(a<q\),

\[
E_a=\sum_{\operatorname{Tr}(x)=1}\chi(x)^a,\qquad
k_a=(aB-1)/2,\qquad t=(p-1)/2.
\]

The reconstruction uses

\[
G_{p^2}(\chi_m)=E_mG_p(\lambda),\qquad
E_m\overline E_m=p
\]

when the restriction is quadratic. Its ideal is explicitly the product of primes indexed by upper-half representatives \(u>m/2\), as in (37). At proper conductors \(m\mid p+1\),

\[
E_m=(-1)^{(p+1)/m+1}.
\]

Appendix G’s root selection uses integral polynomial CRT and

\[
A'(1)\equiv0\pmod m.
\]

The decoder is

\[
\rho_{\mathfrak p}(E_a)=\binom{t}{k_a},\qquad
\binom{2k_a}{k_a}
=(-4)^{k_a}\rho_{\mathfrak p}(E_a).
\]

Hence, with \(z=aB-1=2k_a\),

\[
z!=(-4)^{z/2}\rho_{\mathfrak p}(E_a)\bigl((z/2)!\bigr)^2.
\]

One-child recursion and \(O(B)\)-length corrections give

\[
\widetilde O(q^C+\sqrt{p/q}).
\]

Again, the field, trace, quadratic character, and known \(p\) are essential.

**2.7. Modulo \(p^2\) — §6, Appendix F**

Proposition 6.1 lifts the split identity to

\[
(aB)!\equiv
(-1)^{a-1}J_a^{(2)}(B!)^a
\left(1+\frac{ap}{q}(H_{aB}-H_B)\right)
\pmod{p^2},
\]

where \(J_a^{(2)}=\rho_{\mathfrak p^2}(J_a)\) and

\[
H_t=\sum_{j=1}^t j^{-1}\pmod p.
\]

Appendix F first specializes Young’s congruence to

\[
\frac{
\binom{a(p+1)B}{(p+1)B,\ldots,(p+1)B}
}{
\binom{aB}{B,\ldots,B}
}
\equiv(-1)^{a-1}J_a^{(2)}\pmod{p^2}.
\]

The block expansion used is

\[
\prod_{\substack{1\le j\le(p+1)t\\p\nmid j}}j
\equiv ((p-1)!)^t t!(1+tpH_t)\pmod{p^2}.
\]

For the harmonic correction, Hensel-lift \(\eta\) to \(\widetilde\eta\bmod p^2\), and calculate

\[
Q_t=
\frac{(1-\widetilde\eta^t)^{p-1}-1}{p}\pmod p.
\]

Lemma 6.2 gives

\[
H_{aB}-H_B
=\sum_{t=1}^{q-1}(\eta^{at}-\eta^t)Q_t\pmod p.
\]

This uses a roots-of-unity filter, the binomial congruence

\[
\binom pm\equiv p(-1)^{m-1}/m\pmod{p^2},
\]

and \(H_{p-1}=0\bmod p\).

**Cost:** \(q-1\) exponentiations modulo \(p^2\), or \(\widetilde O(q)\) bit operations, absorbed by reconstruction. Everything else retains the split bound.

**Requirements:** known \(p\); roots lifted relative to \(p\); inverses in \(\mathbb F_p\); Fermat divisibility making each integer quotient \(Q_t\) well-defined.

For larger indices,

\[
(p+t)!\equiv-pt!\pmod{p^2}\quad(0\le t<p),
\]

and \(n!\equiv0\bmod p^2\) for \(n\ge2p\).

**2.8. Modulo \(p^k\) — §7, Appendix E**

The returned recursive object is

\[
P_M(Z)=\prod_{j=1}^{M}(j+pZ)\pmod{p^k},\qquad M<p.
\]

Its \(Z^i\)-coefficient is divisible by \(p^i\), so \(\deg P_M<k\) modulo \(p^k\).

Morita gamma satisfies

\[
\frac{\Gamma_p(x+1)}{\Gamma_p(x)}
=
\begin{cases}
-x,&x\in\mathbb Z_p^\times,\\
-1,&x\in p\mathbb Z_p,
\end{cases}
\qquad
\Gamma_p(x)\Gamma_p(1-x)=(-1)^{a_0(x)}.
\]

Cyclic Gross–Koblitz arguments satisfy

\[
p\gamma_{i+1}=\gamma_i+p-1-y_i\quad(i<d),
\]

\[
p\gamma_{d+u+1}=\gamma_{d+u}+Y_u.
\]

Set \(z_i=\gamma_{i+1}-1\) and \(w_u=-\gamma_{d+u+1}\). Then

\[
\sum_i z_i=\sum_u w_u.
\]

For

\[
\mathcal G_p(z,w)=
\frac{\prod_i\Gamma_p(pz_i)}
{\prod_u\Gamma_p(pw_u)},
\]

Proposition 7.1 is

\[
U_a=(-1)^Y\mathcal G_p(z,w)
\frac{\prod_iP_{y_i}(z_i)}
{\prod_uP_{Y_u}(w_u)}.
\]

The linear gamma term cancels because the sums of shifts agree. The remaining expansion is

\[
\log_p\mathcal G_p(z,w)
=
-\sum_{m\ge1}
\frac{\lambda_mp^{2m+1}}{2m(2m+1)}
\left(\sum_i z_i^{2m+1}-\sum_u w_u^{2m+1}\right).
\]

For \(p\ge5\), \(\mathcal G_p=1\bmod p^3\). Only \(O(k)\) coefficients matter modulo \(p^k\).

The actual algorithm:

| Step | Cost and conditions |
|---|---|
| Determine gamma coefficients, Lemma 7.2 | A grid of \((k+1)^d\) small lattice points supplies a full-rank linear system; only \(O(k)\) reconstructed Gauss values are evaluated. Polynomial in \(k\), with short products of length \(O_r(kE)\). |
| Find interpolation points, Lemma 7.3 | Translate a rounded lattice point along a short vector; polynomial in \(k,\log p\). Produces \(k\) arguments distinct modulo \(p\). |
| Recursive call | One call returning \(P_{\lfloor n/(s-1)\rfloor}(Z)\). |
| Recover parent values | \(O_r(k)\) shifted interval products of length \(O_r(kE)\), plus polynomial work in \(k\). |
| Interpolate | Degree below \(k\); Vandermonde differences are units modulo \(p^k\). Polynomial in \(k,\log p\). |

The branch requiring gamma work ensures \(p>2k+1\), so retained logarithm/exponential denominators are units. If size conditions fail, direct interval products have the required polynomial-in-\(k\) overhead.

For \(n\ge p\), §7.5 first uses

\[
v_p(n!)=\sum_{j\ge1}\left\lfloor n/p^j\right\rfloor.
\]

If this is below \(k\), then \(n<kp\). For \(p>k\), write \(n=ap+m\), obtaining

\[
n!=p^aa!\left(\prod_{b=0}^{a-1}P_{p-1}(b)\right)P_m(a).
\]

This requires two shifted-factorial polynomial computations and \(O(k)\) evaluations.

**Requirements:** known \(p\), its valuation, \(p\)-adic precision, local normalized Gauss sums, and unit interpolation differences. These are stronger requirements than access to modular inversion.

**2.9. Uniform logarithmic improvement — Proposition 4.3, Appendix H**

This separate algorithm has no useful-divisor assumption:

\[
O\!\left(\sqrt p\,(\log p)^{3/2}\sqrt{\log\log p}\right)
\]

bit operations modulo \(p\) or \(p^2\).

Its identity is modular halving:

\[
\begin{aligned}
\prod_{a<j\le b}(j+pz)
\equiv{}&2^{b-a}
\prod_{\lfloor a/2\rfloor<j\le\lfloor b/2\rfloor}(j+pz/2)\\
&\times
\prod_{\lfloor(a+p)/2\rfloor<j\le\lfloor(b+p)/2\rfloor}
(j+p(z-1)/2)
\pmod{p^2}.
\end{aligned}
\]

The routing permutation is \(x\mapsto x/2\) for even \(x\), and \(x\mapsto(x+p)/2\) for odd \(x\).

Two batching results do the arithmetic:

- **Lemma H.1:** at most \(H^{3/8}\) requested prefix endpoints up to \(H\), including \(H\), cost \(O(\sqrt H(\log p)^2)\). It applies BGS Theorem 15 to a two-dimensional recurrence for the coefficients of \(P_m(Z)\bmod p^2\).
- **Lemma H.2:** \(s\) disjoint intervals of total length \(V\), with separate shifts, cost
  \[
  O((s+1+\sqrt V)(\log p)^6).
  \]

Take \(t\asymp\log p\) halving rounds, \(B\asymp\log p/\log\log p\), and \(H\asymp p/B\). Binary-expansion avoidance bounds make the remaining total interval length small. Prefix evaluation dominates the stated bound.

This is a logarithmic improvement measured against the **modulus \(p\)**. It is not a general \(K^{1/2-\delta}\) interval-product algorithm.

---

**3. Exactly what survives modulo an unfactored composite**

The following conclusions are my transfer analysis of the preceding arguments.

| Tal component | Status over \(R_N=\mathbb Z/N\mathbb Z\) |
|---|---|
| Product/remainder trees with monic divisors, Lemma 2.6 | **Survives.** Gives \(\widetilde O(\sqrt M)\) algebra operations for an \(M\)-term interval. |
| Fixed-dimensional cyclotomic lattice rounding, Lemma 5.2 | **Survives as an integer algorithm**, replacing its integer parameter \(p\) by \(N\). The resulting lattice point has no established factorial decoder modulo \(N\). |
| Applying a supplied cyclotomic evaluation map | **Survives when the map exists** and required denominators are units. |
| Hensel lifting a supplied simple root/factor | **Survives under unit-derivative conditions.** It does not construct the missing finite-field identities. |
| Polarized-ideal recovery over a known CM-order | **Remains an available number-field algorithm.** What fails is the justification that the constructed ideal/polarization represents the desired character sum for the composite-modulus factorial problem. |
| Frobenius, character norm, and finite-field orthogonality identities | **Do not transfer with \(p=N\).** |
| Wilson complementation and factorial-ratio decoders | **Do not transfer.** |
| Removing known \(p\)-adic valuations, §§2.6, 7 | **Requires the prime/local component.** Failed modular inversion does not replace valuation removal. |
| Modular-halving routing, Appendix H | Much of the routing **does survive for odd \(N\)**, but its scale becomes \(N\), not the unknown least factor or the short prefix length \(K\). |

There are several distinct breaks.

**First: the necessary order data are wrong.**

\[
q\mid N-1
\]

does not imply \(q\mid P-1\) or \(q\mid Q-1\). Similarly,

\[
q\mid\Phi_r(N)
\]

does not give the corresponding divisibility for either hidden prime. A global exact-order-\(q\) element may also have smaller order in one CRT component.

You can test stronger componentwise root conditions through

\[
\gcd(\eta^{q/s}-1,N),
\]

but passing these tests still does not establish the finite-field factorial identities with \(N\) replacing \(p\).

**Second: Appendix A fails without an inversion failure.**

For \(N=15,y=3\),

\[
1-y^{N-1}=1-3^{14}\equiv7\pmod{15},
\]

although \(y\ne0\). The zero-indicator identity is false. Likewise, \(x\mapsto x^N\) need not be a ring endomorphism, and Wilson gives no replacement: \(14!\equiv0\bmod15\).

These are false identities even in division-free expressions. The “failed inversions expose factors” strategy cannot repair them.

**Third: the higher construction needs characteristic-\(P\) digits and valuations.**

The digit calculation in Appendix D, normalization by \(P^d\), and Wilson conversion all concern a particular hidden prime. Digits in base \(N\) do not encode those calculations.

**Fourth: even unusually favorable roots do not fix the decoder.**

One might successfully construct composite-modulus cyclotomic maps or polarized generators under extra hypotheses. That still leaves the missing assertion that their residues encode the desired factorial ratios. Successful algebraic reconstruction is not a proof of the finite-field decoding congruence.

**Fifth: failed inversion is not automatically factor discovery.**

If \(\gcd(a,N)\) is proper, it gives a factor. If it equals \(N\), a product-tree descent may locate a useful constituent, but an individual zero denominator can remain uninformative.

**Answer to the requested sub-step question:** Tal provides no proved generic operation computing an arbitrary long consecutive product, or its unit predicate, over unfactored \(R_N\) in \(K^{1/2-\delta}\operatorname{polylog}N\) time. Its relevant generic product subroutines are Lemma 2.6 and Lemma H.2, both with square-root dependence on the applicable total length. The sub-square-root work reconstructs **special prime-dependent ratios**, followed by shorter BGS products.

Expanding \(\prod_{j\in S}(x-j)\) densely also has output length \(\Omega(|S|)\) coefficients. A scalar value or compressed representation avoids that output bound, but no such general faster primitive is proved here.

---

**4. BGS: exactly what is batched**

**The elementary construction.** For block length \(b\), construct

\[
F_b(X)=\prod_{i=1}^{b}(X+i)
\]

and evaluate it at \(0,b,2b,\ldots\). These values are consecutive block products. Product trees and ordinary multipoint evaluation cost

\[
O(\mathcal M_R(b)\log b)
\]

ring operations when the number of points is \(O(b)\). Choosing \(b\asymp\sqrt n\) gives the usual factorial algorithm. Monic remainder trees work over arbitrary commutative rings.

**BGS’s improvement is value shifting, not a smaller polynomial degree.** Theorem 5 converts

\[
F(r),F(r+\beta),\ldots,F(r+d\beta)
\]

into the values on the shifted progression \(r+\alpha+i\beta\), in

\[
2\mathcal M_R(d)+O(d)
\]

operations. Lagrange interpolation becomes a convolution involving \(1/(a+i-j)\).

Its precise condition is that

\[
\beta,\ 2,\ldots,d,\ \alpha-d\beta,\ldots,\alpha+d\beta
\]

are units, with the inverse of their product supplied.

**Polynomial matrices.** Theorem 8 batches

\[
M_b(X)=M(X+\alpha b)\cdots M(X+\alpha)
\]

at \(X=0,\beta,\ldots,b\beta\). The recursion uses

\[
M_{2b}(X)=M_b(X+\alpha b)M_b(X),
\]

with an extra linear matrix in the odd case. For matrix dimension \(a\), its cost is

\[
O\!\left(\mathsf{MM}(a)b+a^2\mathcal M_R(b)\right).
\]

Theorem 14 computes a recurrence term at index \(n\) in

\[
O\!\left(
\mathsf{MM}(a)\sqrt n+a^2\mathcal M_R(\sqrt n)
\right).
\]

Theorem 15 obtains the same bound for selected indices \(n_1<\cdots<n_t=n\), provided

\[
t<n^{1/2-\epsilon}
\]

for fixed \(\epsilon>0\). That qualification matters for Tal’s Appendix H. [BGS, §§3–7](https://mathexp.eu/bostan/publications/BoGaSc07.pdf).

**The actual factoring application is slightly different from ordinary consecutive blocks.** In Lemma 13, take \(c=\lfloor\sqrt b\rfloor\) and

\[
F(X)=\prod_{j=0}^{c-1}(X+cj+1).
\]

Evaluate at \(X=0,\ldots,c-1\). These are products over residue classes modulo \(c\), collectively covering \(1,\ldots,c^2\). A short tail completes the range to \(b\).

Lemma 12 batches invertibility checks using a scalar product tree and locates a nonunit by descent. Auxiliary inversion conditions concern explicit integers below \(N\); a failure therefore exposes a nontrivial gcd, followed by factor recovery.

The result is Lemma 13’s bound quoted at the start: **find a prime divisor at most \(b\), or certify none exists**, without knowing the factorization of \(N\). Theorem 11 gives complete factorization in

\[
O(\mathsf M(N^{1/4}\log N))
\]

bits.

For your task, the final multiplication of all block values into a factorial is unnecessary. The expensive batched evaluations remain.

---

**5. Harvey and HH: hidden-factor collisions, rather than factorials**

**5.1. The common polynomial primitive**

Given giant-step residues \(v_1,\ldots,v_t\in R_N\), construct

\[
f(X)=\prod_{j=1}^{t}(X-v_j).
\]

Evaluate at the geometric progression

\[
1,\alpha,\ldots,\alpha^{m-1},
\qquad \alpha\in R_N^\times.
\]

Harvey Lemmas 2.3–2.4 give

\[
\text{construction}:O(t\ell^3),\qquad
\text{geometric evaluation}:O((t+m)\ell^2).
\]

Bluestein’s convolution uses

\[
ij=\binom i2+\binom{-j}2-\binom{i-j}2.
\]

This evaluation needs **only an invertible \(\alpha\)**, not a field, primitive root, known group order, or factorization of \(N\).

A gcd of \(f(\alpha^i)\) with \(N\) detects a collision modulo \(P\) or \(Q\). Exact matches modulo \(N\) are handled separately; after removing them, a gcd equal to \(N\) can be resolved by inspecting constituent differences.

This is Harvey Algorithm 4.1/Proposition 4.1, and HH Algorithm 4.1/Proposition 4.2. [Harvey](https://arxiv.org/pdf/2010.05450v1), [HH](https://arxiv.org/pdf/2105.11105v1).

**5.2. Harvey’s search**

The key congruence is

\[
\alpha^{aQ+bP}\equiv\alpha^{aN+b}\pmod P.
\]

It follows because the two exponents agree modulo \(P-1\). The algorithm computes the right-hand residue modulo \(N\); it does not need to know \(P\).

Lemma 3.3 says that if

\[
(N/r)^{1/2}\le P<\sqrt N,
\]

there exist positive \(a,b\), \(ab\le r\), such that

\[
0\le aQ+bP-\sqrt{4abN}
<
\frac{\sqrt N}{4r\sqrt{ab}}.
\]

Algorithm 4.2 forms

\[
t_{a,b}=\alpha^{aN+b-\lceil\sqrt{4abN}\rceil},
\]

\[
v_{a,b,j}=\alpha^{-jm}t_{a,b},
\qquad
0\le j<\frac{\sqrt N}{4rm\sqrt{ab}},
\]

and searches against all baby steps \(\alpha^i,\ 0\le i<m\).

The total giant-step count is

\[
O\!\left(
\left[\frac{\sqrt N}{\sqrt r\,m}+r\right]\log N
\right).
\]

Proposition 4.2 gives

\[
O\!\left(
\left[\frac{\sqrt N}{\sqrt r\,m}+r\right]\ell^4
+m\ell^2
\right).
\]

A guessed value \(u=aQ+bP\) is checked using

\[
Y^2-uY+abN,
\]

whose roots must be \(aQ,bP\).

The full algorithm chooses

\[
r\asymp N^{1/5}\ell^{-4/5},
\qquad
m\asymp N^{1/5}\ell^{6/5}.
\]

Small-factor elimination establishes the lower-bound hypothesis. Proposition 2.7 constructs a sufficiently large-order element or finds a factor; this setup costs \(O(N^{1/5}\ell^2)\) with \(D\asymp N^{2/5}\). No factorization of \(\varphi(N)\) is supplied.

**There is no special input assumption \(N\equiv c\bmod m\)** in Harvey’s theorem. The useful congruence concerns hidden factors and exponents; all enumerated residues are computed from \(N\).

**5.3. HH’s lattice and sieve refinement**

HH uses a small-prime product \(m\), checks \(\gcd(m,N)\), and enumerates possible reduced residues

\[
P\equiv\sigma\pmod m.
\]

It also enumerates multiplicative cells

\[
\sigma_0\le P<(1+1/m_0)\sigma_0.
\]

Proposition 3.3 computes \(a,b\), using only \(N,m,m_0,\sigma,\sigma_0\), such that the matching hidden factor satisfies

\[
\left|aQ+bP-\left(aN/\sigma_0+b\sigma_0\right)\right|
\le
\frac{4\sqrt{Nm}}{m_0^{3/2}},
\]

and

\[
aQ+bP\equiv aN/\sigma+b\sigma\pmod{m^2}.
\]

The lattice condition is

\[
b\equiv(N/\sigma^2)a\pmod m.
\]

After scaling coordinates, two-dimensional lattice reduction finds the required vector in \(O(\ell^3)\) bits. The factor \(P\) is used only in the correctness proof.

The congruence modulo \(m^2\) reduces the baby-step range to

\[
\lambda=
\left\lceil\frac{4\sqrt N}{(mm_0)^{3/2}}\right\rceil,
\]

using baby steps \(\beta^{m^2i}\), \(0\le i\le2\lambda\). There are

\[
O(\varphi(m)m_0\log N)
\]

giant steps. Proposition 4.4 gives

\[
O\!\left(
m\ell(\log\ell)^2+
\varphi(m)m_0\ell^4+
\frac{\sqrt N\,\ell^2}{(mm_0)^{3/2}}
\right).
\]

Taking a primorial with \(\varphi(m)/m\asymp1/\log\log N\) and

\[
mm_0\asymp
N^{1/5}(\log\log N)^{2/5}/(\log N)^{4/5}
\]

gives Theorem 1.1’s log-log improvement. All hypotheses—coprimality and sufficient order—are established algorithmically. [HH, §§3–4](https://arxiv.org/pdf/2105.11105v1).

These methods do compute polynomials of the requested shape \(\prod(X-v_j)\), but the \(v_j\) are engineered exponentiation residues. Their speedup comes from covering a much larger candidate space with a small collision table, not from multiplying \(K\) consecutive integers in sub-\(\sqrt K\) time.

---

**6. Three concrete directions for the unit predicate**

**Idea 1 — Restrict HH to a certified factor interval. Status: proved adaptation; no novelty or optimality claim.**

Suppose a previous predicate result certifies \(P\ge A\), and we want to decide whether

\[
P\in[A,B),\qquad B\le\sqrt N+1.
\]

Set HH’s sieve modulus \(m=1\), and let \(s=m_0\). Restrict multiplicative cells to this interval. Their number is

\[
J=O(1+s\theta),\qquad \theta=\log(B/A).
\]

For each cell start \(x\), Proposition 3.3 computes \(a,b\) such that

\[
\left|aQ+bP-(aN/x+bx)\right|
\le4\sqrt N/s^{3/2}
\]

whenever \(P\) belongs to that cell.

Set

\[
\lambda=\left\lceil4\sqrt N/s^{3/2}\right\rceil,\qquad
\kappa=2\lambda+1.
\]

Use the same collision polynomial, now with only \(J\) giant steps and \(\kappa\) baby steps.

There is a way to avoid paying the full deterministic high-order setup for every narrow interval. Choose random units \(\beta\), and test

\[
\gcd(\beta^i-1,N)=1,\qquad 1\le i<\kappa.
\]

A proper gcd factors \(N\). A gcd equal to \(N\) causes rejection. Passing proves both component orders are at least \(\kappa\).

If \(\kappa\le A-1\), primitive elements in both CRT components pass. Their probability is \(\Omega((\log\log N)^{-2})\), so this is a Las Vegas setup with only logarithmic expected repetition overhead.

**Cost:**

\[
\widetilde O\!\left(
(1+s\theta)\ell^3+
\left(1+\frac{\sqrt N}{s^{3/2}}\right)\ell^2
\right)
\]

expected bit operations, suppressing the additional logarithmic factors.

Ignoring logarithms, balancing gives

\[
s\asymp N^{1/5}\theta^{-2/5},
\qquad
T=\widetilde O\!\left((1+N^{1/5}\theta^{3/5})\ell^3\right),
\]

subject to the order-test and parameter conditions above.

For \(A\asymp\sqrt N\), interval width \(W=B-A\), and \(\theta\asymp W/\sqrt N\),

\[
T=\widetilde O\!\left(
\left[1+\frac{W^{3/5}}{N^{1/10}}\right]\ell^3
\right).
\]

**Why it is exact:** if \(P\) is inside the interval, one enumerated cell guarantees a hidden-factor collision. Exhausting all cells without a factor therefore excludes \(P\) from the interval. If any collision finds a factor outside the proposed interval, factoring lets us compare \(P\) directly.

**Where it fails to give a general breakthrough:** the lower bound \(P\ge A\) is essential to the cheap randomized order certification. Starting with an interval of constant relative width near \(\sqrt N\) still costs \(N^{1/5+o(1)}\). Subsequent narrowing does not erase that first cost. This is also an interval test for the **prime factor itself**, not a test for arbitrary nonunits in the interval.

I checked the lattice/collision construction on three small semiprimes; the promised collision occurred in each case. The justification is the argument above, not those examples.

**Idea 2 — Compress the zero information of a character product. Status: exact reduction proved; faster evaluation conjectural; simple sampling shortcut refuted.**

For odd squarefree \(N\), define

\[
C_N(K)=\prod_{j=1}^{K-1}\left(\frac jN\right),
\]

using the Jacobi symbol, including its zero value. Then

\[
U(K)\iff C_N(K)\ne0.
\]

This discards virtually all factorial information: the output has only three possibilities. Every individual character value is computable from \(N,j\) alone, without factors, generators, or group-order information.

An equivalent count is

\[
Z_N(n)=
\sum_{j=1}^{n}
\left(1-\left(\frac jN\right)^2\right).
\]

It counts nonunits exactly. For \(n<N\),

\[
Z_N(n)=\left\lfloor\frac nP\right\rfloor+
\left\lfloor\frac nQ\right\rfloor,
\]

and \(U(K)\iff Z_N(K-1)=0\).

**Concrete proposed operation:** evaluate this character-prefix zero count, or only whether it vanishes, without evaluating all its terms.

**Established costs:**

- Individual-symbol enumeration: \(\widetilde O(K\ell)\) bits.
- BGS product followed by a Jacobi symbol or gcd: square-root-in-\(K\) cost.
- A hypothetical \(\widetilde O(K^{1/2-\delta}\ell^c)\) prefix evaluator would give the desired improvement; this is an unproved primitive.

Explicit complete Gauss-sum evaluation does not supply it. For the primitive quadratic character modulo squarefree \(N\),

\[
\tau_N(a)=\left(\frac aN\right)\tau_N(1),\qquad
\tau_N(1)^2=\left(\frac{-1}N\right)N.
\]

Thus complete Gauss sums merely re-express the already easy *individual* symbol. They do not evaluate its long prefix product or the incomplete principal-character count above.

**Sampling obstruction:** when \(n=P<Q\), exactly one index in \(1,\ldots,n\) is nonunit. Uniform independent index sampling needs \(\Theta(P)\) samples for constant detection probability. This refutes that sampling shortcut, not arbitrary algorithms exploiting the explicit arithmetic of \(N\).

This direction isolates a genuinely smaller output, but the missing operation remains substantial: fast incomplete zero-counting for an explicitly evaluable character.

**Idea 3 — Replace factorial nonunits by a smooth-order factor certificate. Status: proved restricted-family method; this is Pollard’s classical mechanism, not a new algorithm.**

Choose a bound \(B\) and construct

\[
E_B=\operatorname{lcm}(1,\ldots,B).
\]

For a random unit \(a\bmod N\), compute

\[
g=\gcd(a^{E_B}-1,N).
\]

If

\[
P-1\mid E_B,\qquad Q-1\nmid E_B,
\]

then \(P\mid g\), and the probability that \(Q\nmid g\) is

\[
1-\frac{\gcd(E_B,Q-1)}{Q-1}\ge\frac12.
\]

A successful trial factors \(N\), answering every \(U(K)\).

**Cost:** \(\log E_B=O(B)\). Sieve construction and modular exponentiation cost

\[
\widetilde O(B\ell)
\]

bits, or \(O(B)\) modular multiplications plus preprocessing. It is below square-root-in-\(K\) cost when \(B\ll\sqrt K\), especially for polylogarithmic \(B\).

**Failure:** it detects smoothness of a factor’s group order, not its location below \(K\). A failed trial gives no negative certificate for \(U(K)\); both factors may also yield \(g=N\). This makes a useful inexpensive factor-discovery branch, followed by BGS or HH for an exact answer.

The most concrete research continuation here is Idea 1’s interval specialization. For a new general prefix primitive, Idea 2 states the unresolved operation precisely: compute whether an incomplete principal-character count is zero, with less work than the corresponding modular product.