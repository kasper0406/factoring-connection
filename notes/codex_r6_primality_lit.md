> Source: GPT-6.1-sol round 6 (codex/r6_primality_lit.md → .out), 2026-10-18. Curated copy; not yet independently verified.

Primality tests **do sometimes expose factors and partial information**. What they generally lack is a guarantee of doing so efficiently on every composite input. The distinction is precise: a failed congruence can fail modulo **both** hidden primes, whereas a useful gcd requires a value that vanishes modulo exactly one.

**Status convention:** **PROVED** covers cited theorems and elementary derivations given below; it does not imply novelty. **CONJECTURED** identifies heuristic assumptions. **UNVERIFIED** identifies a proposed extraction capability or probability bound not established here.

**1. Common setup and the exact CRT criterion — PROVED.** Let \(N=pq\), with distinct odd primes, and sample bases uniformly from \((\mathbb Z/N\mathbb Z)^\times\). For any computed scalar \(z\),

\[
1<\gcd(z,N)<N
\iff
[p\mid z]\mathbin{\mathrm{XOR}}[q\mid z].
\]

Thus \(z\ne0\pmod N\) proves failure of an asserted identity, but need not expose a factor: \(z\) might be nonzero in both CRT components.

For any positive exponent \(e\), define

\[
A_p(e)=\frac{\gcd(e,p-1)}{p-1},\qquad
A_q(e)=\frac{\gcd(e,q-1)}{q-1}.
\]

Cyclicity of the two multiplicative groups and independence under CRT give the exact formula

\[
\boxed{\Pr\bigl(1<\gcd(a^e-1,N)<N\bigr)
=A_p(e)+A_q(e)-2A_p(e)A_q(e).}
\tag{1}
\]

For polynomial or extension-ring computations, apply the same criterion to each recorded coefficient, pivot, denominator, or resultant. A ring residual that vanishes entirely modulo one prime but not the other guarantees a factor from some coefficient gcd. Merely being nonzero modulo \(N\) does not.

All probabilities below condition on unit bases. If bases are uniform in \(1,\ldots,N-1\), an initial gcd finds a factor with probability

\[
\frac{p+q-2}{N-1};
\]

combine this with the conditional probabilities using the unit fraction \(\varphi(N)/(N-1)\). Sampling \(2,\ldots,N-2\) instead removes two always-passing bases and changes the normalization.

**2. Fermat, Miller–Rabin, and Solovay–Strassen — PROVED exact probabilities.** Write

\[
N-1=2^t d,\quad d\text{ odd},\qquad
u=v_2(p-1),\quad v=v_2(q-1),\quad k=\min(u,v).
\]

Put

\[
g_p=\gcd(N-1,p-1),\qquad g_q=\gcd(N-1,q-1),
\]

\[
A=\frac{g_p}{p-1},\quad B=\frac{g_q}{q-1},\quad
F=AB,\quad
C_k=\frac{4^k+2}{3\cdot4^k}.
\]

For semiprimes there is a useful simplification:

\[
\boxed{g_p=g_q=\gcd(p-1,q-1)=g.}
\tag{2}
\]

Also \(t=k\) if \(u\ne v\), while \(t>k\) if \(u=v\).

| Computation or event | Exact probability |
|---|---:|
| Fermat passes: \(a^{N-1}=1\pmod N\) | \(F\) |
| Fermat residual gcd exposes a factor | \(A+B-2F\) |
| Miller–Rabin passes | \(FC_k\) |
| MR finds a nontrivial square root of \(1\) | \(F(1-C_k)\) |
| MR with the additional gcd checks specified below exposes a factor | \(A+B-F(1+C_k)\) |
| Solovay–Strassen passes | \(F/2\) |
| SS residual gcd exposes a factor | \((A+B)/2-F\) |
| SS checks gcds of both \(b-1\) and \(b+1\) | \(A+B-2F\) if \(u=v\); \(A+B-\tfrac32F\) otherwise |

These are factor-extraction probabilities, distinct from the usual compositeness-detection bounds. The original probabilistic tests are [Rabin, *J. Number Theory* 12 (1980), 128–138](https://doi.org/10.1016/0022-314X(80)90084-0), and [Solovay–Strassen, *SIAM J. Computing* 6 (1977), 84–85](https://doi.org/10.1137/0206006), with its [1978 correction](https://doi.org/10.1137/0207009). The displayed extraction formulas are derived here.

- **Fermat** computes \(b=a^{N-1}\pmod N\). The gcd \(\gcd(b-1,N)\) factors precisely when the local order divides \(N-1\) at exactly one prime. Equation (1) gives the table entry.
- **Miller–Rabin** computes \(x_i=a^{2^i d}\pmod N\), accepting when \(x_0=1\), or \(x_i=-1\) for some \(0\le i<t\). A nontrivial-root factor occurs when \(x_i\not\equiv\pm1\pmod N\) but \(x_{i+1}=1\): locally \(x_i=(1,-1)\) or \((-1,1)\).
- **Augmented MR**, as counted in the table, checks \(\gcd(x_i-1,N)\) and \(\gcd(x_i+1,N)\) for \(0\le i<t\), and \(\gcd(x_t-1,N)\). This also catches local \(+1\) or \(-1\) when the other component is not either sign. Arbitrary additional gcds inside modular exponentiation constitute a different experiment.
- **Solovay–Strassen** computes \(b=a^{(N-1)/2}\pmod N\) and \(J=(a/N)\in\{\pm1\}\), accepting iff \(b=J\pmod N\). Its residual gcd factors precisely when \(b=J\) at one prime but not the other.

For SS, each local matching event has probability \(A/2\) or \(B/2\): the other prime supplies an independent Legendre-symbol sign. Their intersection has probability \(F/2\), giving the residual-gcd formula.

A concrete counterexample to “rejection yields a factor” is \(N=35,a=2\). Fermat, MR, and SS reject, yet the Fermat residual and all the MR/SS sign residuals above have gcd \(1\) with \(35\).

**Balance alone does not determine the probabilities.** If \(p,q=\Theta(\sqrt N)\) and \(g\) is bounded by a polynomial in \(\log N\), augmented MR has probability \(O(g/\sqrt N)\), while nontrivial-root extraction has probability \(O(g^2/N)\). These are exponentially small in the input length. Large common divisors of \(p-1,q-1\) can change this substantially, even for balanced factors.

**Carmichael numbers are different — PROVED.** If \(N\) is Carmichael, every unit satisfies \(a^{N-1}=1\). Consequently, every MR rejection produces a nontrivial square root of \(1\), and therefore a factor if the preceding value is retained.

More precisely, let \(N=\prod_{j=1}^r p_j\), \(r\ge3\), and \(u_j=v_2(p_j-1)\). Define

\[
w_u(0)=2^{-u},\qquad
w_u(j)=2^{j-1-u}\quad(1\le j\le u).
\]

Then

\[
\Pr(\text{no factor from MR})
=\sum_{j=0}^{\min u_i}\prod_{i=1}^r w_{u_i}(j)
\le 2^{1-r}\le\frac14.
\]

Thus one unit base factors a Carmichael number with probability at least \(3/4\); independent repetitions amplify this. This is the nontrivial-root mechanism underlying [Miller’s 1976 paper, §2](https://www.cs.cmu.edu/~glmiller/Publications/Papers/RiemannsHypothesisandTestsforPrimality.pdf), combined with the strong-test counting argument associated with Rabin and [Monier, *Theoretical Computer Science* 12 (1980), 97–108](https://doi.org/10.1016/0304-3975(80)90007-9).

**3. Lucas and strong Lucas tests — PROVED, with the sampling model specified.** For parameters \(P,Q\), let

\[
D=P^2-4Q,\qquad \gcd(2QD,N)=1,
\]

and let \(\rho,\sigma\) be the roots of \(X^2-PX+Q\). Then

\[
U_m=\frac{\rho^m-\sigma^m}{\rho-\sigma},
\qquad V_m=\rho^m+\sigma^m.
\]

Set \(z=\rho/\sigma\). This is a norm-one element in the quadratic algebra. At a prime \(r\mid N\), put \(\chi_r=(D/r)\). The local norm-one torus is cyclic of order

\[
h_r=r-\chi_r.
\]

Exactly,

\[
U_m=0\pmod r\iff z^m=1,
\qquad
V_m=0\pmod r\iff z^m=-1.
\tag{3}
\]

The ordinary Lucas test computes \(U_L\), where \(L=N-(D/N)\). The strong test writes \(L=2^t d\) and checks \(U_d=0\), or \(V_{2^i d}=0\) for some \(i<t\). These definitions and local identities appear in [Grantham, *Frobenius Pseudoprimes*, §2, Theorems 2.1–2.2](https://arxiv.org/pdf/1903.06820).

Consequently:

- \(\gcd(U_m,N)\) exposes a factor iff the local order of \(z\) divides \(m\) at exactly one prime.
- \(\gcd(V_m,N)\) exposes a factor iff the local power is \(-1\) at exactly one prime.

For **uniform independent elements of the two full tori**, all Fermat/MR formulas above apply with \(p-1,q-1,N-1\) replaced by \(h_p,h_q,L\). In particular,

\[
\Pr(\gcd(U_L,N)\text{ proper})
=
\frac{\gcd(L,h_p)}{h_p}
+\frac{\gcd(L,h_q)}{h_q}
-
2\frac{\gcd(L,h_p)\gcd(L,h_q)}{h_ph_q}.
\tag{4}
\]

This is **not automatically the distribution obtained by random Lucas parameters**.

An implementable alternative fixes \(D\) and samples \(P\pmod N\), conditioned on \(P^2-D\) being a unit, setting \(Q=(P^2-D)/4\). The map

\[
P\longmapsto z=\frac{P+\sqrt D}{P-\sqrt D}
\]

is locally a bijection onto the torus with its identity removed. Hence replace the local probability of \(z^m=1\) by

\[
\frac{\gcd(m,h_r)-1}{h_r-1}.
\]

For a strong chain, the local masses for \(\operatorname{ord}(z^d)=2^j\) are

\[
f_r(0)=\frac{\gcd(d,h_r)-1}{h_r-1},\qquad
f_r(j)=\frac{\gcd(d,h_r)2^{j-1}}{h_r-1}
\quad(1\le j\le v_2(h_r)).
\]

Let \(S_r=\sum_{j=0}^{t}f_r(j)\), taking out-of-range masses as zero. Monitoring both \(U_{2^i d}\) and \(V_{2^i d}\) for \(i<t\), plus \(U_L\), gives exact factor probability

\[
\boxed{S_p+S_q-S_pS_q-\sum_{j=0}^{t}f_p(j)f_q(j).}
\tag{5}
\]

These sampling formulas are elementary cyclic-group counts, derived here.

**Baillie–PSW.** The common strong variant combines base-2 MR with a strong Lucas test using Selfridge parameters: search \(D=5,-7,9,-11,\ldots\) for \((D/N)=-1\), then \(P=1,Q=(1-D)/4\), with square/gcd precautions. Factors arise through those gcds or the MR/Lucas conditions just described. For a fixed \(N\), the standard procedure is deterministic, so there is no random-base factor probability to multiply between its stages. The stages are not independent. See [Baillie–Fiori–Wagstaff, *Strengthening the Baillie–PSW primality test* (2021), §§2–3](https://homes.cerias.purdue.edu/~ssw/bfw.pdf).

**UNVERIFIED:** a general correctness theorem or universal factor-extraction probability for standard BPSW. The cited paper discusses heuristic arguments predicting counterexamples; excellent observed detection is not a factoring theorem.

**4. Grantham’s Frobenius tests.** The general test, for monic \(f\) of degree \(d\), computes

\[
f_0=f,\quad
F_i=\operatorname{gcmd}(X^{N^i}-X,f_{i-1}),\quad
f_i=f_{i-1}/F_i.
\]

It requires \(f_d=1\), checks \(F_i(X^N)\equiv0\pmod{F_i}\), and checks

\[
(-1)^{\sum_{2\mid i}\deg(F_i)/i}
=\left(\frac{\operatorname{disc}f}{N}\right).
\]

**PROVED factor condition:** a nonunit leading coefficient during polynomial Euclid exposes a factor. Different local gcd degrees guarantee this obstruction. More generally, coefficient gcds expose factors according to the common CRT criterion. A failed Frobenius identity may instead have unit residual coefficients. [Grantham, §3, especially Proposition 3.5 and the test definition](https://arxiv.org/pdf/1903.06820).

For the distinct **random quadratic Frobenius test**, \(f=X^2-bX-c\), admissible parameters satisfy

\[
\left(\frac{b^2+4c}{N}\right)=-1,\qquad
\left(\frac{-c}{N}\right)=1.
\]

After trial division and square detection, it checks that \(X^{(N+1)/2}\) is scalar, that \(X^{N+1}=-c\), and a strong chain for \(N^2-1=2^r s\), with \(-1\) checked through index \(r-2\). Grantham proves an error bound below \(1/7710\) for the prescribed randomized procedure. **That bounds passing, not failure to reveal a factor.** [Grantham, *A Probable Prime Test With High Confidence*, *JNT* 72 (1998), 32–47, §§2–3](https://arxiv.org/pdf/1903.06823).

**UNVERIFIED:** transferring that \(1/7710\) bound into a factor-extraction bound. A numerical extraction probability requires a specified polynomial distribution and specified gcd instrumentation.

**5. AKS, APR-CL, ECPP, and Pocklington–Pratt.**

| Method | What is computed | Exact circumstances exposing factors; probability qualification |
|---|---|---|
| **AKS** | Perfect-power detection; auxiliary \(r\) with sufficiently large \(\operatorname{ord}_r(N)\); small gcds; \((X+a)^N\stackrel?=X^N+a\pmod{N,X^r-1}\). | Perfect-power/small-gcd branches can expose factors. A failed polynomial check exposes one through coefficient gcds exactly when some coefficient vanishes at one hidden prime only. Fixed input and prescribed parameters give a deterministic outcome, not a random-base probability. |
| **APR-CL** | Jacobi sums in cyclotomic rings, prescribed exponentiations and saturation checks, followed by a restricted divisor search. | Coefficient/denominator gcds obey CRT. If the algorithm reaches the final stage with its divisor restrictions established, a composite is factored by that final search. Early rejection need not supply those restrictions. No canonical random-base probability. |
| **ECPP** | CM constructions provide a curve and candidate order \(m=c\ell\); scalar multiplications check \([m]P=O\) and that \([c]P\) is nonzero modulo every prime divisor; recursively certify \(\ell\). | A slope denominator or projective coordinate vanishing at exactly one prime exposes a factor. Failure of an order check can occur at both primes and reveal none. Exact fixed-curve probabilities are given below. |
| **Pocklington** | Factor \(F\mid N-1\), \(F>\sqrt N\); certify its prime factors \(\ell\); check \(a_\ell^{N-1}=1\) and \(\gcd(a_\ell^{(N-1)/\ell}-1,N)=1\). | A proper gcd directly factors; gcd \(N\) merely invalidates that witness. For any specified exponent, equation (1) gives the exact probability for a uniform unit base. |
| **Pratt** | Complete factorization of \(N-1\), recursive certificates, and a witness of order \(N-1\), checked by exponentiation and inequalities for every prime \(\ell\mid N-1\). | Optional gcds have equation (1)’s probability. A failed inequality or exponentiation need not factor. Short certificates do not imply polynomial-time certificate construction. |

Primary sources are [AKS, *PRIMES is in P*, §4](https://www.cse.iitk.ac.in/users/manindra/algebra/primality_v6.pdf); [Cohen–Lenstra, *Math. Comp.* 42 (1984), 297–330, §§2–3 and 7–12](https://pub.math.leidenuniv.nl/~lenstrahw/PUBLICATIONS/1984a/art.pdf); [Morain, ECPP presentation, §3](https://arxiv.org/pdf/math/0502097); [Grégoire–Théry–Werner, formal Pocklington certificates](https://www-sop.inria.fr/members/Benjamin.Gregoire/Publi/pock.pdf); and [Pratt, *SIAM J. Computing* 4 (1975), 214–220](https://doi.org/10.1137/0204018).

Some technical details matter:

- **AKS — PROVED.** A failed check supplies a short tuple identifying a nonzero coefficient in the *folded* residual. It need not identify an original binomial coefficient. For comparison, before folding, if \(p\) is the least prime divisor, then
  \[
  \binom Np\equiv N/p\pmod N,
  \]
  which does expose a factor. Locating that degree efficiently is a different task from finding an AKS mismatch.
- **APR-CL — PROVED.** Its underlying prime-modulus Gauss-sum identity is
  \[
  \tau(\chi)^N\equiv \chi(N)^{-N}\sigma_N(\tau(\chi))\pmod N,
  \]
  where \(\sigma_N\) acts on character roots of unity. Jacobi sums implement the checks in smaller rings. Ultimately it can establish that every divisor \(r\mid N\) satisfies \(r\equiv N^i\pmod S\) for some \(i<t\), with \(S>\sqrt N\); enumerating these residues then finds every possible small divisor. Cohen–Lenstra explicitly distinguish this successful final-stage restriction from the early elimination of ordinary composites.
- **ECPP — PROVED.** A standard sufficient certificate has prime \(\ell>(N^{1/4}+1)^2\). For a fixed smooth curve and uniform independent CRT points, write
  \[
  E(\mathbb F_r)\simeq C_{d_{r,1}}\times C_{d_{r,2}}.
  \]
  Then
  \[
  H_r(e):=\Pr([e]P=O\bmod r)
  =\frac{\gcd(e,d_{r,1})\gcd(e,d_{r,2})}{\#E(\mathbb F_r)}.
  \]
  A suitable projective \(Z\)-coordinate of \([e]P\) exposes a factor with exact probability
  \[
  H_p(e)+H_q(e)-2H_p(e)H_q(e).
  \]
  This is a fixed-curve, uniform-point statement; it is not an overall ECPP factor probability. Polynomial generation-time estimates for practical fast ECPP are heuristic, while certificate soundness is proved.

**6. The discrete-log view — PROVED, with an ambiguity caveat.** For a unit \(a\), let \(s=p+q\). Since

\[
\varphi(N)=N-s+1,
\]

Euler’s theorem gives

\[
\boxed{a^{N-1}=a^{s-2},\qquad
a^N=a^{s-1},\qquad
a^{N+1}=a^s\pmod N.}
\tag{6}
\]

Thus the primality-test residue supplies a discrete-log equation for a short hidden exponent. It determines that exponent **modulo \(\operatorname{ord}_N(a)\)**, not necessarily as an integer. For example, \(a=1\) supplies no discrimination.

If \(q/p\le C\), then

\[
2\sqrt N\le s\le
(\sqrt C+1/\sqrt C)\sqrt N.
\]

The interval has width \(O(\sqrt N)\), so baby-step/giant-step searches it in \(O(N^{1/4})\) group operations and comparable storage. Uniqueness follows if the order exceeds the interval width. Otherwise collisions require candidate validation, further bases, or order-handling arguments. Given the correct \(s\), factor using

\[
p,q=\frac{s\mp\sqrt{s^2-4N}}2.
\]

For factors whose ratio merely stays bounded, \(s=\Theta(\sqrt N)\); the sharper \(s\sim2\sqrt N\) requires \(q/p\to1\).

**Quadratic tori — all four cases, PROVED.** Continue using the norm-one element \(z=\rho/\sigma\), and set

\[
\chi_p=(D/p),\quad\chi_q=(D/q),\quad J=\chi_p\chi_q.
\]

Locally \(z^p=z^{\chi_p}\), \(z^q=z^{\chi_q}\). Hence

\[
H=(p-\chi_p)(q-\chi_q)
=N-\chi_qp-\chi_pq+J
\]

annihilates \(z\), yielding

\[
z^{N+1}=z^{\chi_qp+\chi_pq+1-J},
\qquad
z^{N-J}=z^{\chi_qp+\chi_pq-2J}.
\tag{7}
\]

| \((\chi_p,\chi_q)\) | \(z^{N+1}\) | \(z^{N-J}\), the Lucas-test exponent |
|---|---|---|
| \((+1,+1)\) | \(z^{p+q}\) | \(z^{p+q-2}\) |
| \((-1,-1)\) | \(z^{-(p+q)}\) | \(z^{-(p+q)-2}\) |
| \((+1,-1)\) | \(z^{q-p+2}\) | \(z^{q-p+2}\) |
| \((-1,+1)\) | \(z^{p-q+2}\) | \(z^{p-q+2}\) |

These identities follow from the local Frobenius relations underlying Lucas tests; the complete table is derived here.

**Normalization matters.** An ordinary Lucas root \(\rho\) has norm \(Q\), and instead satisfies

\[
\rho^r=Q^{(1-\chi_r)/2}\rho^{\chi_r}.
\]

Its \(N+1\) identities therefore have extra factors \(Q^K\), with

\[
K=0,\ s,\ p-1,\ q-1
\]

in the four table rows, respectively. The displayed simple exponent identities apply directly to norm-one elements, including roots with \(Q=1\). From computed Lucas values one can recover the normalized power efficiently:

\[
z^m=\frac{V_m+\sqrt D\,U_m}{V_m-\sqrt D\,U_m}.
\]

The same interval-DLP interpretation therefore applies, including its order ambiguity.

**7. Prior art using precisely this connection.**

| Source | What is established or proposed |
|---|---|
| **Shanks (1971), “Class number, a theory of factorization, and genera,” *Proc. Sympos. Pure Math.* 20, 415–440.** | BSGS/order-search machinery and class-group factoring. This is the search primitive; credit for the explicit semiprime identity should be distinguished from credit for BSGS. Pollard’s paper below explicitly describes Shanks’s meet-in-the-middle method. |
| **Pollard (1978), “Monte Carlo methods for index computation (mod \(p\)),” *Math. Comp.* 32, 918–924, §3, p. 923.** | Explicitly proposes solving \(a^{2x}=a^{N+1}\pmod N\), with \(x=(p+q)/2\), assuming \(M<p<q<N/M\), \(M>N^{1/3}\). Gives an approximate \(O(\sqrt{N/M})\) kangaroo analysis and acknowledges possible failure. This is the requested short-exponent connection, equivalently a totient search. [Primary-paper transcription](https://studylib.net/doc/10510783/monte--carlo--methods--for--index-...-by-j.--m.-pollard). |
| **McKee–Pinch (1996), “Old and new deterministic factoring algorithms,” ANTS-II, LNCS 1122, 217–224.** | Removes the Monte Carlo element from Pollard’s discrete-log factoring method. Harvey identifies its bound as \(N^{1/3+o(1)}\). [Publisher abstract](https://doi.org/10.1007/3-540-61581-4_57), [Harvey §1](https://arxiv.org/pdf/2010.05450). |
| **Hittmeir (2018), “A babystep-giantstep method for faster deterministic integer factorization.”** | Searches \(a^X=a^{N+1}\pmod N\) while sieving possible \(p+q\) residues using \(pq=N\). Proves a deterministic bound of the form \(N^{1/4}\exp[-\Omega(\log N/\log\log N)](\log N)^{O(1)}\). [Primary text, Theorems 1.2 and 7.1](https://arxiv.org/pdf/1608.08766). |
| **Hittmeir (2020 preprint; 2021 publication), “A time-space tradeoff for Lehman’s deterministic integer factorization method.”** | Replaces the unweighted sum search by Lehman-type weighted sums and proves deterministic \(N^{2/9+o(1)}\) factoring. [Theorem 1.1](https://arxiv.org/abs/2006.16729). |
| **Harvey (2020 preprint; 2021 publication), “An exponent one-fifth algorithm for deterministic integer factorisation.”** | Batches weighted searches and collisions modulo hidden factors, proving \(O(N^{1/5}\log^{16/5}N)\) bit complexity. [Theorem 1.1, §3, Algorithms 4.1–4.3](https://arxiv.org/pdf/2010.05450). |

The weighted identity is elementary and explains the mechanism. For weights \(A,B\) and unit base \(a\),

\[
a^{Aq+Bp}\equiv a^{AN+B}\pmod p.
\]

Thus a candidate \(w=Aq+Bp\) can create a collision modulo \(p\) without creating one modulo \(q\), enabling a gcd. Once \(w\) is correct,

\[
w^2-4ABN=(Aq-Bp)^2
\]

validates it and recovers the factors. The improvement comes from organizing and batching these searches, rather than extracting factors from an ordinary primality-test rejection.

**CONJECTURED/heuristic:** the usual square-root running-time analysis of an implementation of Pollard’s kangaroo rests on random-walk assumptions. Deterministic BSGS bounds and the cited deterministic factoring bounds do not require that heuristic.

**8. The logical implication — PROVED, but weaker than the proposed slogan.** Suppose a polynomial-time algorithm, on every composite \(N\), produces \(h\) of polynomial length, and a polynomial-time decoder \(D(N,h)\) always returns a proper divisor. Then factoring is polynomial-time: use deterministic primality testing at each node and recursively split composites. The factor tree has \(O(\log N)\) nodes.

An inverse-polynomial success probability with verifiable divisors similarly gives randomized expected polynomial-time factoring by repetition.

What does **not** follow is that compositeness certificates cannot encode factors cheaply. A factor is already an \(O(\log N)\)-bit certificate, efficiently checked. The missing requirement is **efficient uniform production together with efficient decoding**. Consequently:

- **PROVED:** a universally efficient factor-producing primality test would give efficient factoring.
- **CONJECTURED:** classical factoring is not randomized polynomial-time.
- **UNVERIFIED:** any unconditional prohibition against recovering factors from primality-test transcripts.

AKS illustrates the distinction: its failed folded-polynomial identity is efficiently checkable evidence of compositeness; its correctness proof does not require extracting a factor from that failure.

For comparison, factoring computes \(\varphi(N)\) and Carmichael’s exponent \(\lambda(N)\) efficiently. Conversely, an explicitly given positive, polynomial-bit-length multiple of \(\lambda(N)\) permits randomized splitting with probability at least \(1/2\) per unit trial, using repeated squaring and a nontrivial root of \(1\). Therefore computing either function is randomized polynomial-time equivalent to factoring. For a semiprime, knowing \(\varphi(N)\) gives \(p+q=N+1-\varphi(N)\) directly.

Historical qualification: [Miller (1976), §3](https://www.cs.cmu.edu/~glmiller/Publications/Papers/RiemannsHypothesisandTestsforPrimality.pdf) gives deterministic equivalences conditional on ERH. [Long, “Random Equivalence of Factorization and Computation of Orders,” Princeton TR 284 (1981)](https://www.lix.polytechnique.fr/Labo/Francois.Morain/Introuvables/long-orders.pdf) establishes the randomized order-computation equivalence. Computing a power \(a^N\) is much weaker than supplying its order or a group-wide annihilating exponent.

**9. Exactly what MR/SS leak about the 2-adic structure — PROVED.** There is a complete sampling law for experiments. Write

\[
p-1=2^u m_p,\quad q-1=2^v m_q,\qquad
c_r=\frac{\gcd(d,m_r)}{m_r}.
\]

Independently at each prime:

1. The odd-order component of \(a^d\) vanishes with probability \(c_r\).
2. Independently, its 2-primary order is \(2^{K_r}\), where
   \[
   \Pr(K_r=0)=2^{-u_r},\qquad
   \Pr(K_r=j)=2^{j-1-u_r}\quad(1\le j\le u_r).
   \]
3. If the odd component survives, no squaring-chain value is \(\pm1\). Otherwise \(x_i=1\) locally iff \(K_r\le i\), and \(x_i=-1\) iff \(K_r=i+1\).
4. The local Legendre symbol of the original base is negative iff \(K_r=u_r\).

This generates the joint distribution of every MR sign pattern and the base’s Jacobi symbol. **The valuations \(u,v\) alone do not specify the unconditional distribution:** the odd-component probabilities also matter.

For the observable global sign flags, the formulas simplify to

\[
\Pr(x_i=1)=F\,4^{\min(i,k)-k},\qquad 0\le i\le t,
\]

\[
\Pr(x_i=-1)=
\begin{cases}
F\,4^{i-k},&0\le i<k,\\
0,&k\le i<t.
\end{cases}
\tag{8}
\]

First-hit probabilities follow by subtracting successive \(+1\) probabilities. Conditional on Fermat passing, \(K_p,K_q\) are independent with distribution \(w_k\); this also proves the MR probabilities above.

For the project’s bit, impose the necessary promise \(N=pq\equiv1\pmod4\). Then

\[
\epsilon=\left(\frac{-1}{p}\right)
=\left(\frac{-1}{q}\right),
\]

and

\[
\boxed{
s\equiv\epsilon(N+1)\pmod8,\qquad
\epsilon=(-1)^{(N+1-s)/4}.}
\tag{9}
\]

Thus, given \(N\), \(s\bmod8\) determines \(\epsilon\). Under this promise, \(\epsilon=+1\) is equivalent to both primes being \(1\bmod4\), to \(-1\) being a square modulo \(N\), and to \(N\) being a sum of two integer squares.

MR gives a detectable positive certificate:

\[
x_i=-1\pmod N,\ i\ge1
\quad\Longrightarrow\quad
x_{i-1}^2=-1\pmod N
\quad\Longrightarrow\quad \epsilon=+1.
\]

Its exact per-base probability is

\[
\boxed{
\Pr(\text{MR supplies this square root of }-1)
=
F\,\frac{4^k-4}{3\cdot4^k}.}
\tag{10}
\]

It is zero when \(\epsilon=-1\), since then \(u=v=k=1\). When \(\epsilon=+1\), it lies between \(F/4\) and \(F/3\). The certificate is cheap to recognize, but its occurrence can be exponentially rare.

SS exposes another exact structural event. With \(b=a^{(N-1)/2}\),

\[
\boxed{
\Pr(b=-1\pmod N)=
\begin{cases}
0,&u=v,\\
F/4,&u\ne v.
\end{cases}}
\tag{11}
\]

Therefore observing \(b=-1\) certifies \(u\ne v\). Under \(N\equiv1\pmod4\), this also certifies \(\epsilon=+1\). Absence of this event does not certify \(\epsilon=-1\); in particular it never occurs when \(u=v\ge2\).

Jacobi symbols alone do not resolve the hidden signs:

\[
(-1/N)=(-1/p)(-1/q),\qquad
(2/N)=(2/p)(2/q).
\]

Their values are already determined by \(N\bmod4\) and \(N\bmod8\). For a random unit base, its Jacobi symbol is an unbiased sign regardless of the hidden primes’ residue classes.

**UNVERIFIED:** a polynomial-time decoder of \(\epsilon\), or a growing prefix of \(p+q\), from arbitrary rejected MR/SS residues without factoring. The rare certificates above establish partial leakage, not such a decoder.

**10. Experiment specification.** Record separately: rejection, initial gcd, Fermat residual gcd, nontrivial-root extraction, augmented sign gcds, and positive partial-information certificates. Stratify semiprimes by \(u,v\) **and** \(g\), and state every parameter distribution. For interval-DLP experiments, retain and validate collisions instead of assuming the first matching exponent is \(p+q\).

As consistency checks, I exhaustively verified the MR/SS probability formulas over every unit base for all **276 pairs of odd primes below 100**, and checked both torus exponent identities in all four Legendre-sign cases on **1,781 valid small parameter choices**. These checks supplement the derivations; they do not establish claims about uninstrumented implementations.