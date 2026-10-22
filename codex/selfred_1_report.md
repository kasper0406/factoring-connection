# Factoring self-reductions and partial information: phase 22 blind report

Agent `selfred`; 22 October 2026. Repository read only; this report and its checks are in the scratch directory. Literature was checked using web search and accessible papers. This is a research report, not a claim of a new factoring algorithm.

**Labels.** **KNOWN** means supported by the linked reference; **PROVED** means a proof is supplied below; **HEURISTIC**, **CONJECTURED**, **OPEN**, and **FALSE** have their usual meanings. **UNVERIFIED** is an additional audit flag, not a mathematical status. An **OPEN** literature status below means that this search found neither the indicated algorithm nor the indicated reduction; it is not a theorem that none exists. Complexity is classical bit complexity, with ordinary binary integer inputs and outputs. Randomized and deterministic reductions are distinguished.

## 1. Verdict

1. **OPEN — search finding.** I found no efficient transformation that moves the unknown large factors to varied, publicly specified nonmonomial functions of those factors and supplies the analogue of RSA message blinding. This is not an impossibility theorem and is not equivalent to saying that factoring has been proved non-random-self-reducible.
2. **PROVED.** A height-free, finite-sample strengthening of PSR holds. For nonzero \(h\in\mathbb Z[X]\) of degree \(d\), and \(s>d\) distinct \(q\in(p,2p)\), every common divisor of the \(h(pq)\), after removing coefficient content and powers of \(p\), divides a Vandermonde determinant built from \(d+1\) of the \(q\)'s in an interval of length at most \(dp/(s-d)\). Thus its primes are \(O(d\log p/\delta)\)-smooth at relative prime density \(\delta\), and its bit length is \(O(d^2\log(d\log p/\delta))\). See Theorem H.
3. **PROVED, with corrected hypotheses.** A public, nonconstant, bounded-degree rational expression \(R_N(p)=N'\) gives a polynomial-time factoring algorithm by clearing denominators and finding rational roots. Polynomial separated products have only the monomial constant case. Rational separated products can cancel identically without being monomials. Small *enumerable* unknown cofactors or errors are allowed; unrestricted approximation is not covered. See Theorem R.
4. **PROVED.** An unconditional exclusion of all such transformations would imply \(\mathrm{FACTOR}\notin\mathrm{FP}\). An unconditional exclusion of polynomial-time Turing reductions from factoring to squarefreeness or \(\epsilon\) would do so too: a factoring algorithm could ignore its oracle.
5. **OPEN — search finding.** For arbitrary integers, the reverse reductions \(\mathrm{FACTOR}\le_T\omega,\mu,\lambda,\mathrm{SF},\operatorname{rad}\), and for balanced semiprimes \(\mathrm{FACTOR}\le_T\epsilon\), were not located. These targets are easy given factors. On the distinct-semiprime promise, the signature-only quantities are constants; unrestricted oracle queries remain a different issue.
6. **PROVED / KNOWN.** Radical and full squarefree decomposition are polynomial-time Turing equivalent; at most \(\lfloor\log_2 N\rfloor\) radical queries suffice. The parity squarefree kernel is also Turing equivalent to them. Shallit–Shamir give reductions *to* the divisor-counting function, not a proof that these functions factor integers. Buchmann–Lenstra equate finding a maximal order with computing the radical, not with factoring.
7. **KNOWN, correcting the background.** Allender–Saks–Shparlinski cover **squarefreeness**, primality, and gcd: Theorem 1 embeds \(\mathrm{TC}^0\) in the corresponding oracle circuit classes; Corollary 2 excludes each exact problem from \(\mathrm{AC}^0[\ell]\) for every prime \(\ell\). These are not semiprime prediction lower bounds. [ASS]
8. **PROVED conditional on GRH; PROVED unconditionally in a smaller range.** With independent dyadic prime factors, a function of \(N\bmod m\) alone has vanishing advantage for \(\epsilon\), under GRH, throughout \(m\le X^{1-\eta}\), where \(p,q\asymp X\) and \(\eta>0\) is fixed. This is wider than \(m\le N^{1/4-o(1)}\). For powers of two an unconditional range \(2^j\le N^c\), some absolute \(c>0\), follows from Bourgain's prescribed-digit theorem. Distributional details and proofs are in §4.
9. **PROVED.** There is an honest, exponentially tiny positive advantage: one Fermat square test factors a fraction \(\Omega(N^{-1/4})\) of independently sampled dyadic semiprimes. For an exactly fair mixture of the two same-mod-4 factor populations, predicting \(\epsilon\) on those inputs and guessing otherwise has that order of advantage over the optimal constant. This is a factoring event, not a nonsplitting leak or a near-prime generator artefact. No new efficient nonsplitting predictor is claimed.

## 2. Self-reduction candidates and algebraic barriers

### 2.1 What counts as moving a factor

**PROVED — distinction.** A statement \(g_r(p)\mid N'\) is weaker than saying that \(g_r(p)\) is prime or that it is the complete factorization of \(N'\). For example, \(p+1\) is composite for odd \(p>2\). The barriers below already apply to moving *divisors*, so they also apply to the stronger prime-factor requirement when its hypotheses hold. The polynomial descriptions below are public symbolic descriptions; the numerical value \(g_r(p)\) is not assumed known.

**PROVED — distinction.** Factor movement, inversion-preserving self-reduction, random self-reducibility over a distribution of moduli, and information-theoretic instance hiding are four separate requirements. Multiplication can provide the first kind of trivial inversion-preserving relation without producing an independently random RSA modulus or hiding \(N\) from an unbounded observer.

### 2.2 Candidate table

| Candidate | Classification | Exact reason / limitation |
|---|---|---|
| \(N'=kN\), or \(kN^j\) | **PROVED: works only as a multiplier / monomial move** | The old primes remain \(p,q\); extra factors come from \(k\). Factoring (kN) factors \(N\), but the secret has not become \(p+r\), for example. If \(\gcd(k,N)\) is proper, the computation already splits \(N\). |
| Fixed \(h(N)\), polynomial \(h\ne0\) | **PROVED: PSR barrier** | If fixed \(g(p)\mid h(pq)\) for all sufficiently large primes \(p\) and all prime \(q\in(p,2p)\), then \(g=cX^j\). Corollary H1 gives a proof. This is a common-divisibility assertion, not a classification of all factors of \(h(N)\). |
| Families \(h_n\), degree polynomial in \(n\), arbitrary heights; random coefficients independent of the instance | **PROVED: Theorem H applies to each fixed family member / random seed** | Positive-density success in \(q\) allows only content, \(p\)-powers, and a small smooth remainder. Arbitrary coefficient heights are allowed in the theorem; polynomial-time output additionally requires polynomial bit size. A probabilistic success assertion does not justify fixing a seed that works everywhere. |
| A completely general \(N\)-dependent polynomial map | **OPEN outside the stated hypotheses** | Theorem H cannot compare different polynomials \(h_{pq}\) as \(q\) changes. Theorem R applies when an explicit separated algebraic identity is supplied. |
| \(\operatorname{Norm}_{K/\mathbb Q}(\alpha(N))\), fixed field and polynomial coordinates | **PROVED: polynomial barrier** | The norm is the determinant of multiplication by \(\alpha\), hence a polynomial in \(N\). Theorem H applies if a moving divisor is asserted uniformly. |
| Norm with an \(N\)-dependent field or algebraic description | **OPEN in general; Theorem R when a public separated identity exists** | A norm is not intrinsically factoring-equivalent. Polynomial field degree and polynomial-size structure constants permit determinant evaluation; arbitrary succinct, exponentially large degrees do not come with that guarantee. |
| Resultants / discriminants with polynomial dependence on \(N\) | **PROVED: polynomial barrier in that setting** | Determinant formulas give polynomials in \(N\). A shared root of *two* polynomials modulo \(N\) makes their resultant divisible by \(N\); a repeated root of \(f\) makes \(\operatorname{Res}(f,f')\) divisible by \(N\). One known simple root of \(f\) alone says neither. Off-\(N\) factors are not automatically functions of \(p\). |
| \(\Phi_m(N)\), \(m>1\) | **PROVED: not a direct \(p\)-divisibility move; H applies** | \(\Phi_m(0)=1\), so \(p\nmid\Phi_m(N)\). Other proposed uniform \(g(p)\) divisors face H. Its ordinary output length is \(\Theta(\varphi(m)\log N)\). Exponential \(m\) has exponential degree; polynomial degree permits polynomial-size evaluation in the ordinary index representation. |
| \(\Phi_m(x)\), with \(x\) related to \(N\) | **PROVED / OPEN depending on the relation** | If \(x\equiv0\pmod N\), the preceding constant-term argument applies. Known roots of unity modulo \(N\) yield order tests and gcds. If \(x=x(N)\) is polynomial, compose and use H. No general barrier is asserted for arbitrary efficiently computed integer representatives \(x(N)\). |
| \(\#E(\mathbb Z/N\mathbb Z)\), good reduction | **PROVED: chosen-twist cardinality oracle factors \(N\)** | CRT gives \(\#E(\mathbb F_p)\#E(\mathbb F_q)\), not a publicly evaluable low-degree polynomial in \(p,q\). Four twist counts give a randomized factoring reduction, proved in §2.6. Computing one specified curve's count is not established here to be factoring-equivalent. |
| Torus order \((p+1)(q+1)\) | **PROVED: factoring-equivalent on semiprimes** | Subtract \(N+1\) to obtain \(p+q\); solve a quadratic. |
| \((p^2+p+1)(q^2+q+1)\) | **PROVED: factoring-equivalent on semiprimes** | With \(s=p+q\), the value is \(s^2+(N+1)s+N^2-N+1\). Its positive \(s\) is recoverable by a quadratic equation. A supplied value \(p^2+p+1\) alone gives \(p\) directly. |
| Lucas / recurrence computations for orders involving \(p+1\) or \(p^2+p+1\) | **PROVED: distinguish residue computations from order computation** | Fast recurrence evaluation modulo \(N\) does not compute the integer group order. Smooth-order annihilation and component separation can split \(N\). There is no identity here giving the full moved order cheaply from \(N\). |
| \(a^N-1\), fixed integer \(a\ge2\) | **PROVED: exponential output length** | \(a^p-1\mid a^{pq}-1\), but the output has \(\Theta(N\log a)\) bits. Even the moved divisor has \(\Theta(p\log a)\) bits. Modular exponentiation computes a residue, which need not retain the integer divisibility. |

**PROVED — cyclotomic size detail.** From \(\Phi_m(z)=\prod_{d\mid m}(z^d-1)^{\mu(m/d)}\), take logarithms for real \(z\ge2\). The main term is \(\varphi(m)\log z\); the remaining absolute sum is at most \(O(2^{\omega(m)}/z)\). Also \(2^{\omega(m)}\le2\varphi(m)\), so for growing \(z=N\) the bit length is \(\Theta(\varphi(m)\log N)\), uniformly in \(m\). The elementary product formula gives

\[
\frac{\varphi(m)^2}{m}=\prod_{\ell^a\parallel m}\ell^{a-2}(\ell-1)^2\ge\frac12.
\]

Every factor is at least one except possibly \(\ell=2,a=1\), which contributes \(1/2\). Thus \(\varphi(m)\ge\sqrt{m/2}\). Conversely \(\varphi(m)\le\operatorname{poly}(n)\) forces \(m\le2\varphi(m)^2=\operatorname{poly}(n)\), permitting trial division of the index and direct polynomial-size evaluation. Allowing an exponential-length symbolic value changes the computational problem.

### 2.3 Theorem H: a height-free, finite-sample PSR bound

**PROVED.** Let \(p\) be a positive prime, \(h\in\mathbb Z[X]\setminus\{0\}\), \(d=\deg h\), and \(c=\operatorname{cont}(h)>0\). Let \(Q\subset(p,2p)\cap\mathbb Z\) contain \(s>d\) distinct integers; they need not be prime. Let \(G\ge1\) satisfy

\[
G\mid h(pq)\qquad(q\in Q).
\]

Write \(z_{p'}\) for the largest divisor of \(z\) coprime to \(p\), and put

\[
R=\left(\frac{G}{\gcd(G,c)}\right)_{p'}.
\]

If \(d=0\), then \(R=1\). If \(d\ge1\), there are \(q_0<\cdots<q_d\) in \(Q\), all in an interval of integer length

\[
H\le \left\lfloor\frac{dp}{s-d}\right\rfloor,
\]

such that

\[
R\mid\Delta:=\prod_{0\le i<j\le d}(q_j-q_i).
\tag{H}
\]

Consequently \(R\) is \(H\)-smooth and

\[
\log_2 R\le\frac{d(d+1)}2\log_2 H.
\tag{H-size}
\]

*Proof.* Order all samples as \(u_1<\cdots<u_s\). Each gap \(u_{j+1}-u_j\) occurs at most \(d\) times in

\[
\sum_{i=1}^{s-d}(u_{i+d}-u_i)<dp.
\]

One block of \(d+1\) consecutive samples therefore has diameter at most the stated \(H\). Use that block as the \(q_i\).

Let \(f(Y)=h(pY)/c\in\mathbb Z[Y]\), with coefficient vector \(a\). For every prime \(\ell\ne p\), some coordinate of \(a\) is a unit modulo \(\ell\): \(h/c\) is primitive, and multiplying its coefficient of degree \(i\) by \(p^i\) does not change this property. If \(\ell^t\mid R\), then all \(f(q_i)\) are divisible by \(\ell^t\). Let \(V=(q_i^j)_{0\le i,j\le d}\). The identity

\[
\operatorname{adj}(V)Va=(\det V)a
\]

shows that \(\ell^t\) divides every coordinate of \(\Delta a\). A unit coordinate of \(a\) forces \(\ell^t\mid\Delta\). This holds for every prime power in \(R\), proving \(H\). Each positive difference in \(\Delta\) is at most \(H\), giving smoothness and (H-size). ∎

**PROVED — density corollary.** Define relative density as

\[
s\ge\delta\bigl(\pi(2p)-\pi(p)\bigr).
\]

The prime number theorem gives \(s\gg\delta p/\log p\). Provided \(s\ge2d\),

\[
H\ll d\log p/\delta,
\qquad
\log_2R=O\!\left(d^2\log(2+d\log p/\delta)\right).
\]

For \(d,\delta^{-1}\le\operatorname{poly}(n)\), with \(p\) an \(O(n)\)-bit prime, both the smoothness bound and the remainder's bit length are polynomial in \(n\). The only analytic input is a lower bound of the right order for the total number of primes in \((p,2p)\); fixed-modulus PNT is sufficient, for example [MV] (§13.1) with the unconditional PNT in its background.

**PROVED — necessary exceptions.** The hypotheses \(h\ne0\) and \(s>d\) matter. A polynomial vanishing at all sampled \(pq\)'s gives no restriction on a common divisor. The content cannot be bounded without a height hypothesis. Powers of \(p\) cannot be excluded because \(h(X)=X^j\) supplies them. If \(\delta\) is so small that \(s\le d\), the proposed polynomial smoothness conclusion need not hold.

**PROVED — information interpretation, and its limit.** For fixed public \(h,p,Q\), \(R\) can be described by at most the number of bits in (H-size), and every prime of \(R\) can be found by trial division up to \(H\). This is a size/smoothness bound, **not** an \(o(1)\)-information theorem about \(p\). Small smooth values can carry residue information. For example, with \(h(X)=X-1\) and \(Q\) consisting of primes \(q\equiv1\pmod4\), divisibility by \(4\) distinguishes \(p\equiv1\pmod4\) from \(p\equiv3\pmod4\). The theorem also does not supply an algorithm to extract a hidden \(G\) from an arbitrary output integer. If a displayed component contains a lone \(p\)-power, its gcd with \(N\) already factors \(N\).

**PROVED — Corollary H1 (fixed PSR).** Suppose fixed nonzero \(g,h\in\mathbb Z[X]\) satisfy \(g(p)\mid h(pq)\) for all sufficiently large primes \(p\) and all primes \(q\in(p,2p)\). Write \(g(X)=X^j g_0(X)\), \(g_0(0)\ne0\). For \(p>|g_0(0)|\), \(p\nmid g_0(p)\). Theorem H with fixed \(d\) gives

\[
|g_0(p)|\le c\,(O(\log p))^{d(d+1)/2}.
\]

A nonconstant fixed polynomial grows as a positive power of \(p\), contradicting that bound. Thus \(g_0\) is constant, so \(g=c_0X^j\). The same proof works for any fixed positive relative density. ∎

**PROVED — scope for randomized families.** Fixing a public seed and a size class fixes \(h\), so H is height-free even for randomized coefficients. If success holds only on a subset \(Q_r\), use its actual size. Averaging success probabilities permits statements about the distribution of \(|Q_r|\), not a seed good for every \(p,q\). An \(N\)-dependent choice \(h_N\) destroys the shared-polynomial premise unless separately controlled. The theorem itself makes no independence assumption; it requires the *same* polynomial at the displayed samples.

### 2.4 Theorem R: publicly described algebraic movement reveals the factors

**PROVED, using KNOWN polynomial factorization [LLL].** Input consists of \(N=pq\), a rational number \(Y\), and a public rational function

\[
R_N(X)=U(X)/V(X)\in\mathbb Q(X),
\]

represented by explicit polynomials of degree and coefficient bit length polynomial in \(n=\lceil\log_2N\rceil\). Suppose \(V(p)\ne0\), \(R_N(p)=Y\), and \(R_N\) is nonconstant. Then \(N\) can be factored in deterministic polynomial time from this tuple.

*Proof.* Clear coefficient denominators and form the nonzero polynomial \(F(X)=U(X)-YV(X)\). It has \(p\) as an integer root. Deterministic polynomial factorization over \(\mathbb Q\) [LLL] finds all its linear factors in time polynomial in its explicit encoding. Test their positive integer roots by division into \(N\); \(p\) is among them. ∎

**PROVED — separated product formulation.** Given explicit nonzero rational functions \(A,B\), a *known nonzero scalar* \(C\), and an exact identity

\[
N'=C\,A(p)B(q),
\]

take \(R_N(X)=C A(X)B(N/X)\). Degrees remain polynomial; clearing negative powers introduces only polynomially many factors of \(N\) into the coefficients, so coefficient bit lengths remain polynomial. All denominators must be nonzero at the actual factors. If a rational expression \(C(X,N/X)\) is used instead of a scalar, apply the theorem to the **whole** expression, including that cofactor; cancellation must be checked there.

**PROVED — exact exception for polynomials.** For nonzero polynomials \(A,B\), \(A(X)B(N/X)\) is constant iff \(A=aX^j\), \(B=bX^j\), with the same exponent. Indeed, any nonzero complex root of \(A\) would be a root of the constant Laurent polynomial, which is impossible; likewise for \(B\). Thus both have roots only at zero and are monomials; their exponents must agree. The product is \(abN^j\), the multiplier/monomial case.

**FALSE — rational exception claimed to be only monomial.** For rational functions, the constant case is more general. Take

\[
A(X)=X+1,\qquad B(Y)=\frac{Y}{N+Y}.
\]

Then \(A(X)B(N/X)=1\) identically, even though neither function is monomial. The correct exception is **rational-function identity**, tested after cancellation. This example concerns rational factor expressions; it does not claim two positive integer moved divisors.

**PROVED — polynomially bounded unknown cofactor.** If \(N'=K R_N(p)\), \(1\le K\le n^c\), enumerate \(K\), form \(R_N(X)-N'/K\), and apply the proof. The same works for \(R_N(p)=K N'\). The bound concerns the **value** of \(K\), not its bit length: an \(O(n)\)-bit cofactor has exponentially many possibilities. If \(R_N\) is nonconstant, every enumerated equation is nonzero.

**PROVED — approximation that is actually covered.** If \(R_N(p)\) is integral and 

\[
|R_N(p)-N'|\le n^c,
\]

enumerate the integer errors. More generally the set of candidate exact values must be publicly enumerable in polynomial time. A tiny interval for rational values is not enough without denominator control.

**OPEN — unqualified approximation version.** Relative error, even exponentially small in \(n\), does not provide such a theorem. Already \(N'=N\) approximates \((p+1)(q+1)\) with relative error \(O(N^{-1/2})\) on balanced inputs. A theorem extracting factors from that approximation with no further input would itself prove factoring polynomial-time. This is a counterexample to the proposed *algebraic argument*, not an unconditional proof that the conclusion is false (factoring in P has not been excluded).

**PROVED — other qualifications.** The data is \((N,N',A,B,C)\), not just \((N,N')\), unless the descriptions can be recomputed from those inputs. For randomized maps, the necessary coins or descriptions must be auxiliary output. Explicit degree and coefficient bounds are needed: an arithmetic circuit can denote a polynomial of exponential degree or height. Zero factors, poles, and identically zero equations must be excluded. A prime factor merely hidden somewhere in a large cofactor gives no exact equation of this form.

**PROVED — escapes from R.** A modular move \(u\equiv ap\pmod M\) introduces an unknown carry in \(ap=u+kM\). If it has exponentially many possible values, R does not permit enumerating it. Boundedly many carries restore the argument. Floor functions, bit operations, arbitrary computable \(g\), implicit polynomials of exponential degree, and large unknown cofactors are not covered. No impossibility result for those escapes is established here.

### 2.5 Why RSA blinding escapes this argument

**PROVED.** RSA changes \(c=x^e\bmod N\) to \(c'=c r^e\bmod N\). The corresponding secret is the residue \(x'=rx\bmod N\), so \(x'=rx-kN\) with an unknown carry. The public relation

\[
(rx-kN)^e\equiv c r^e\pmod N
\]

is a consequence of \(x^e\equiv c\pmod N\), not an independent nonzero integer polynomial in \(x\). Even fixed-degree root finding **modulo an unknown composite factorization** is different from rational root finding over \(\mathbb Q\). Lifting the congruence to an equality introduces another unknown quotient. These variables range over exponentially many possibilities in general. When \(\gcd(e,\varphi(N))=1\), multiplication by uniform \(r\in(\mathbb Z/N\mathbb Z)^*\) makes \(xr\) uniform in that group, while preserving an efficiently reversible relation between the messages.

**PROVED.** In factoring, substituting \(q=N/X\) into \(pq=N\) gives the identity \(X(N/X)=N\). That equation alone reveals nothing new. But a *new exact nonconstant* separated expression evaluated at \(p\) gives the independent equation used in R. Multipliers leave the identity constant. This explains the algebraic difference; it does not prove that a group action on a different encoding of the factoring secret cannot exist.

**KNOWN.** Alexi–Chor–Goldreich–Schnorr's RSA/Rabin bit reductions concern message inversion; Håstad–Näslund concern individual RSA/Rabin message bits. They are not statements about a fixed bit of a factor. Håstad–Näslund **Theorems 8.1 and 9.1** were checked in the author-hosted manuscript (RSA and their injective Rabin variant, respectively). ACGS's exact theorem numbering was not independently checked here. [ACGS], [HN]

### 2.6 Explicit order-value reductions

**PROVED — torus / cyclotomic order.** If \(T=(p+1)(q+1)\), then \(s=T-N-1=p+q\), and \(p,q\) are the roots of \(X^2-sX+N\). For \(T=(p^2+p+1)(q^2+q+1)\), expand to obtain

\[
T=s^2+(N+1)s+N^2-N+1.
\]

It is strictly increasing for \(s>0\), so the correct integer \(s\) is found by solving that quadratic and checking. Thus these exact integer outputs, if computable, would factor \(N\).

**PROVED — chosen elliptic twist counts.** Let \(N=pq\) with odd distinct primes. Choose a nonsingular short Weierstrass curve \(E\) modulo \(N\); a proper gcd of its discriminant with \(N\) already factors \(N\). For each unit \(d\), its quadratic twist has local counts

\[
\#E^{(d)}(\mathbb F_r)=r+1-\chi_r(d)t_r\quad(r=p,q).
\]

This identity follows by counting solutions in \(dy^2=f(x)\): the contribution at each \(x\) is \(1+\chi_r(d)\chi_r(f(x))\). CRT multiplies the local counts.

Choose independent uniform units \(d_1,d_2\); their Legendre-sign vectors in \(\{\pm1\}^2\) are independent uniform vectors. With probability \(6/16=3/8\), they generate all four vectors. In that event, the counts for \(E,E^{(d_1)},E^{(d_2)},E^{(d_1d_2)}\) sum to

\[
4(p+1)(q+1).
\]

Use the resulting candidate \(s\) to solve \(X^2-sX+N\), and verify any factors. Repeat on failure. This is a constant expected number of cardinality queries; one count can be reused across attempts. No assumption on nonzero traces is needed. Uniform units are sampled by rejection and gcd. This proves a reduction for an oracle allowing selected twists, not for the count of one specified curve. ∎

**KNOWN / PROVED — converse.** Given \(p,q>3\), Schoof's deterministic polynomial-time finite-field point-counting algorithm computes both local counts, and CRT gives the composite count. Handle the finitely many smaller primes directly. Thus the chosen-curve/twist cardinality problem and semiprime factoring are equivalent under randomized polynomial-time Turing reductions. René Schoof, *Elliptic curves over finite fields and the computation of square roots mod p*, Mathematics of Computation **44**(170), 483–494 (1985), introduction. [Schoof]

### 2.7 Meta-barrier and random self-reducibility literature

**PROVED — meta-barrier.** If factoring were in FP, compute \(p,q\) and output \((p+r)(q+r)\) with public random positive \(r\) of polynomial bit size. This is a nonmonomial moved-divisor map of ordinary polynomial output size. Consequently, an unconditional theorem excluding *all* such polynomial-time maps implies factoring not in FP. If “known prime factors” literally requires prime outputs, one can instead sample and certify primes in intervals depending on \(p\), after factoring; the simpler shifted-divisor example is enough for the intended formulation. If transformations are required to be independently random moduli or information-theoretically hide \(N\), those additional definitions must be analyzed separately.

**PROVED — no-reduction version.** If factoring is in FP, it has a polynomial-time oracle reduction to every total oracle problem, using zero queries. Therefore “no polynomial-time Turing reduction from factoring to \(\epsilon\)/squarefreeness” implies factoring not in FP. The corresponding exclusion of randomized reductions implies factoring not in randomized polynomial time. A lower bound only in a restricted query model does not have this implication.

**KNOWN.** Martín Abadi, Joan Feigenbaum, Joe Kilian, *On hiding information from an oracle*, JCSS **39**(1), 21–50 (1989), §6 **Question 1**, explicitly asks whether FACTORING is encryptable leaking at most input length, and whether it is random-self-reducible. Its general negative results concern NP-hard functions under its hiding definitions, not a factoring impossibility theorem. [AFK]

**KNOWN.** Joan Feigenbaum, Lance Fortnow, *On the random-self-reducibility of complete sets*, SICOMP **22**(5), 994–1005 (1993): **Theorem 3.1** puts an NP language with a nonadaptive polynomial-query random self-reduction in their \(\mathrm{AM}_{poly}\); **Corollary 3.2** improves the advice to logarithmic; **Corollary 3.3** gives a third-level PH collapse for NP-complete languages; **Theorem 3.5** treats adaptive \(O(\log n)\)-query reductions for NP-complete languages. The conclusion does not forbid a reduction for factoring. [FF]

**PROVED, using KNOWN primality certificates [Pratt].** For example, the language “\(N\) has a prime divisor at most \(b\)” lies in NP∩coNP: a complete prime factorization certifies a no-answer, and a small prime divisor certifies a yes-answer. Polynomial-size primality certificates can be supplied. It is polynomial-time Turing equivalent to factoring by threshold search. Thus the NP-*completeness* hypothesis in the PH-collapse application is not available for factoring. Membership in NP∩coNP alone is not an obstruction to random self-reducibility.

**KNOWN.** Donald Beaver, Joan Feigenbaum, *Hiding instances in multioracle queries*, STACS 1990, LNCS **415**, 37–48, provides multioracle instance hiding / locally random reductions, with the noncommunication requirement on the oracles. This is not a map to independent factoring instances; the oracle function can change. The general construction is also described in [FF] (introduction, refs. 7–8). [BF]

**KNOWN.** Andrej Bogdanov, Luca Trevisan, *On worst-case to average-case reductions for NP problems*, FOCS 2003, 308–317; journal version SICOMP **36**(4), 1119–1159 (2006), extends the NP-complete barrier to nonadaptive self-correctors for samplable distributions. It does not establish non-random-self-reducibility of factoring or distinguish factoring from any of the natural bits in §3. [BT]

**PROVED — fixed-modulus Rabin distinction.** For a fixed \(N\), if \(y\) is a quadratic residue, choose uniform unit \(r\), ask for a square root of \(yr^2\), and divide the answer by \(r\). This randomizes the residue while leaving the modulus fixed. To factor \(N=pq\), choose uniform unit \(x\), submit \(x^2\), and compare the returned root \(z\) with \(x\). Conditional on the square, the unknown \(x\) is uniform among four roots; with probability \(1/2\), \(\gcd(x-z,N)\) is proper. An oracle succeeding on a fraction \(\alpha\) of uniform unit squares gives success \(\alpha/2\). This is the mechanism underlying the factoring equivalence of modular square roots; see Michael O. Rabin, *Digitalized signatures and public-key functions as intractable as factorization*, MIT/LCS/TR-212 (1979). [Rabin]

**OPEN — search findings.** I found no worst-case-to-average-case random self-reduction that changes a general factoring modulus to an independent balanced semiprime modulus, including the Blum-integer restriction. I also found no theorem providing this through the hard nonsmooth part while only randomizing smooth parts. A polynomially bounded smooth part is directly computable by trial division; multiplying by random known smooth cofactors preserves the hard primes and is not such a result. The accurate qualified statement is: **factoring is not known, in this literature search, to be random self-reducible over the modulus**. Per-modulus residue randomization should always be named separately.

## 3. Partial-information functions: status and reduction directions

### 3.1 Definitions and status table

**PROVED — conventions.** For \(M=\prod r^{e_r}\), define

\[
\omega(M)=\#\{r:e_r>0\},\quad \Omega(M)=\sum e_r,\quad
\lambda(M)=(-1)^{\Omega(M)},\quad \operatorname{rad}(M)=\prod_{e_r>0}r.
\]

Here \(\lambda\) is **Liouville's** function, not Carmichael's exponent. Write \(\mathrm{SF}(M)\) for squarefreeness. There are three different squarefree outputs:

* radical \(\operatorname{rad}(M)\), the largest squarefree divisor;
* parity kernel \(s(M)=\prod_{e_r\text{ odd}}r\), with \(M=u^2s(M)\);
* full decomposition \(M=\prod_{i\ge1}b_i^i\), where the \(b_i\) are squarefree and pairwise coprime.

**PROVED.** Complete factoring computes every row below in polynomial time. The reverse directions concern an oracle on its **full stated domain**, not merely the value at the original semiprime.

| Target oracle / data | Value on distinct \(N=pq\) | Best supported classification in this search |
|---|---|---|
| \(\omega(M)\) | \(2\) | **OPEN:** no general polynomial-time algorithm or reduction from factoring found. **PROVED / KNOWN:** \(\omega\le_T d\) by signature recovery, §3.3; \(\omega\le_T\mathrm{FACTOR}\). C9/O9 in the readable FOCS list [Adl]. |
| Möbius \(\mu(M)\) | \(+1\) | **OPEN** for computation and factoring reduction on arbitrary inputs. **PROVED:** \(\mathrm{SF}(M)=[\mu(M)\ne0]\); \(\mu(M)=\lambda(M)\mathrm{SF}(M)\). **KNOWN:** one \(d\)-oracle call suffices [SS]; proof below. |
| Liouville \(\lambda(M)\) | \(+1\) | **OPEN:** no factoring reduction or polynomial-time algorithm found. **PROVED:** determined by the prime signature and thus reducible to \(d\); not the factoring-equivalent Carmichael function. |
| Squarefreeness \(\mathrm{SF}(M)\) | yes | **OPEN:** no polynomial factoring reduction found. **KNOWN:** deterministic reduction to Euler's full totient [Woll]. **KNOWN:** exact squarefreeness is outside \(\mathrm{AC}^0[\ell]\) for prime \(\ell\) [ASS]. BHK's analytic certificates and their costs are discussed below. |
| \(\operatorname{rad}(M)\) | \(N\) | **OPEN:** no polynomial factoring reduction found. **PROVED:** equivalent under deterministic polynomial-time Turing reductions to full squarefree decomposition and the parity kernel (§3.2). **KNOWN:** radical computation is polynomially equivalent to maximal-order construction [BL] (Theorem 1.3). |
| Parity kernel \(s(M)\), or largest square divisor \(u^2=M/s(M)\) | \(N\), respectively \(1\) | **PROVED:** Turing equivalent to radical / decomposition (§3.2). This does **not** say the single numerical output \(s(N)\) reveals a semiprime's factors. C7 is this parity-kernel problem [Adl]. |
| Full squarefree decomposition | \(b_1=N\) | **PROVED:** same Turing equivalence; on a squarefree input it returns the input as a single block. |
| Smaller factor \(p\bmod4\) | nonconstant | **OPEN:** factoring reduction and general prediction amplification not found. On \(N\equiv1\pmod4\), it is exactly the \(\epsilon\) bit; on \(N\equiv3\pmod4\), it is an orientation bit (the two factors have opposite classes). |
| \(\epsilon=(-1/p)\), \(N\equiv1\pmod4\) | \(+1\) for \(p,q\equiv1\pmod4\), otherwise \(-1\) | **OPEN:** no factoring reduction or polynomial-time algorithm found. **PROVED:** equals fixed-base QR of \(-1\), and is determined by \(\varphi(N)\bmod8\). Seurin formulates the corresponding Blum/pseudo-Blum distributional assumption [Seurin]. |
| \(\varphi(N)\bmod\ell\), fixed odd prime \(\ell\) | depends on \(p+q\bmod\ell\) | **PROVED:** equivalent to \(p+q\bmod\ell\), because \(\varphi(N)=N+1-p-q\). **OPEN:** no general factoring reduction or polynomial-time algorithm found for a fixed small nontrivial modulus. The \(\Phi\)-hiding assumption is a distributional assumption about a divisibility predicate, not an unconditional equivalence [CMS]. |
| \(\varphi(N)\bmod2\), \(\bmod4\), \(\bmod8\) | \(0,0\), and \(0\) or \(4\) | **PROVED:** first two are trivial for odd distinct semiprimes; the third determines \(\epsilon\) on \(N\equiv1\pmod4\). Thus calling every fixed totient residue a nontrivial open bit is inaccurate. |
| Fixed-base quadratic residuosity \(\mathrm{QR}_a(N)\) | depends on local Legendre symbols | **PROVED:** factoring suffices; a known integer-square base is trivial, a Jacobi \(-1\) base is a nonresidue, and \(a=-1\) gives \(\epsilon\) under the promise. **OPEN:** no general reverse reduction found for nontrivial fixed bases. Random-base QR assumptions must not be silently substituted for a theorem about each fixed base. |

**KNOWN.** The classical quadratic-residuosity cryptographic assumption is associated with Shafi Goldwasser and Silvio Micali, *Probabilistic encryption*, JCSS **28**(2), 270–299 (1984). This is an assumption, not a proof of a reverse factoring reduction. [GM]

**KNOWN.** Christian Cachin, Silvio Micali, Markus Stadler, *Computationally private information retrieval with polylogarithmic communication*, EUROCRYPT 1999, LNCS **1592**, 402–414, introduces \(\Phi\)-hiding, involving a small odd prime dividing the totient of a composite with unknown factors. Yannick Seurin, *On the lossiness of the Rabin trapdoor function*, PKC 2014, LNCS **8383**, 380–398, introduces the \(2\)-\(\Phi/4\)-hiding extension distinguishing Blum from pseudo-Blum moduli in its specified distributions. Neither reference proves that an arbitrary \(\epsilon\) oracle factors. [CMS], [Seurin]

### 3.2 Radical, parity kernel, and squarefree decomposition

**PROVED — at most \(\lfloor\log_2 M\rfloor\) radical calls.** Put \(X_0=M\). Until \(X_i=1\), query \(R_i=\operatorname{rad}(X_i)\) and set \(X_{i+1}=X_i/R_i\). If \(M=\prod r^{e_r}\), then

\[
X_i=\prod r^{\max(e_r-i,0)},\qquad R_i=\prod_{e_r\ge i+1}r.
\]

For \(t=\max e_r\), there are \(t\le\lfloor\log_2M\rfloor\) calls. Put \(R_t=1\) and \(b_j=R_{j-1}/R_j\), \(1\le j\le t\). These are precisely the pairwise coprime products of primes of exact exponent \(j\), so \(M=\prod b_j^j\). Conversely the decomposition gives the radical by multiplying the nonunit blocks. Query lengths never exceed the input length. ∎

**PROVED — parity kernel suffices too.** Given an oracle \(s(X)\), recursively set \(X_{i+1}=\sqrt{X_i/s(X_i)}\). Its prime exponents are \(\lfloor e_r/2^{i+1}\rfloor\); \(s(X_i)\) is the product of primes whose exponent has binary digit \(i\) equal to one. Therefore

\[
\operatorname{rad}(M)=\operatorname{lcm}_{i\ge0}s(X_i).
\]

There are at most \(1+\lfloor\log_2\max e_r\rfloor\) nonterminal calls. Radical, in the other direction, yields the decomposition by the preceding algorithm and hence the parity kernel. A largest-square-divisor oracle yields \(s(X)\) by division, and vice versa. These are function-oracle reductions; squarefree **decision** does not automatically give these numeric outputs. ∎

**PROVED — \(p^2q\) example.** On the promise \(M=p^2q\), distinct primes, squarefreeness is constantly false, \(\mu=0\), \(\lambda=-1\), \(\omega=2\), and \(d=6\). In contrast, \(\operatorname{rad}(M)=pq\) gives \(p=M/\operatorname{rad}(M)\), while \(s(M)=q\) gives \(p=\sqrt{M/q}\). Thus these *function values* factor that promised family immediately. This does not produce an algorithm for computing those values, or a reduction from factoring balanced \(pq\) to factoring \(p^2q\).

**OPEN — special-case search finding.** I found no polynomial-time algorithm for general \(p^2q\) with two comparable unknown primes. Small factors are removed by trial division, perfect powers are detected by integer root tests, and sufficiently strong extra congruence/size information can invoke small-root methods. None justifies calling the general \(p^2q\) problem easy. The separate literature on \(p^rq\) for large \(r\) is not being asserted to settle \(r=2\).

### 3.3 What a divisor-counting oracle supplies

**KNOWN.** Jeffrey Shallit, Adi Shamir, *Number-theoretic functions which are equivalent to number of divisors*, IPL **20**(3), 151–153 (1985), proves that \(\mu\) needs one \(d\)-oracle call, and that functions depending only on the exponent multiset can be obtained with at most \(\log_2M\) such calls. The author-affiliated publication record states these claims explicitly. [SS]

**PROVED — one-call construction for \(\mu\).** Let \(b=\lfloor\log_2M\rfloor\), and choose a prime \(r>b+1\) of size \(O(b)\) (searching a short polynomial-sized range suffices). Query

\[
D=d(M^{r-1})=\prod_i((r-1)e_i+1).
\]

Since \(1\le e_i\le b<r\), a factor \((r-1)e_i+1\) is divisible by \(r\) exactly when \(e_i=1\), in which case it equals \(r\). For \(e_i>1\) it is greater than one and coprime to \(r\). Thus \(D\) is a pure \(r\)-power iff \(M\) is squarefree. In that case \(D=r^{\omega(M)}\), giving \(\mu(M)=(-1)^{\omega(M)}\); otherwise \(\mu(M)=0\). The query \(M^{r-1}\) has only \(O(b^2)\) bits. Handle \(M=1\) directly. ∎

**PROVED — full signature with at most \(b\) calls.** The polynomial

\[
P(T)=\prod_i(1+e_iT)
\]

has degree at most \(b\), \(P(0)=1\), and \(P(j)=d(M^j)\) for \(j=1,\ldots,b\). Interpolate over the rationals. For each integer \(1\le e\le b\), repeatedly divide by \(1+eT\) to find its multiplicity. This recovers the exponent multiset and hence \(\omega,\lambda,d,\mu\) and squarefreeness. All queries and arithmetic have polynomial bit size. ∎

**PROVED — limitation.** This implication is **not** \(d\le_T\mathrm{FACTOR}\) in reverse. These polynomial roots are the small exponents, not the large prime bases. Recovering the signature \({1,1}\) of a semiprime supplies no factor. No converse factoring reduction from \(d\) was found.

### 3.4 What signature-only queries can do

**PROVED — extraction lemma.** For a positive, polynomial-bit-length query integer \(X\), write the \(N\)-supported part formally as \(p^a q^b\). If \(a\ne b\), \(N\) is factorable in polynomial time from \(N,X\), without the oracle: while \(\gcd(X,N)=N\), divide \(X\) by \(N\). At the first unequal exponent, its gcd with \(N\) is \(p\) or \(q\). There are at most \(\log_2X\) iterations. If both exponents are zero it is already coprime. Equivalently, \(\gcd(X,N^k)\), with \(k\ge\log_2X\), represents the supported part, but a proper gcd with \(N\) can emerge only after removing the common copies.

**PROVED — restricted-query barrier.** Consider queries made only by multiplication and powering from \(N\) and completely factored public integers, with polynomial output length. Either a public integer already shares a proper gcd with \(N\), or \(p,q\) occur with the same exponent. In the latter case the whole signature is known from the construction (merge public factors with gcd/division); a signature oracle supplies no additional information about the unknown primes.

**OPEN — no unrestricted extrapolation.** An arbitrary query \(X=f(N)\), even one coprime to \(N\) such as a suitable shifted integer, can have an unknown signature dependent on all bits of \(N\). There is no proof here that its oracle answer is uninformative. The extraction lemma only covers unequal exponents of the **old** primes. A reduction could exploit the other factors of transformed queries, and might use adaptive queries. Constancy of the target on semiprimes, multiplicativity on simple products, or the gcd lemma does not rule out this possibility.

### 3.5 What the requested historical references actually establish

**KNOWN / UNVERIFIED distinction.** Leonard M. Adleman, *Algorithmic number theory—the complexity contribution*, FOCS 1994, 88–113, contains the following readable problem list on printed pp. **98–99**. I checked a text reproduction of the paper; its mathematical reduction symbol has imperfect OCR. The definitions show that the random reduction is the authors' randomized oracle-reduction notion. [Adl]

| Number in the readable FOCS paper | Content, in modern notation |
|---|---|
| **C7** | From \(M\), return \(u,s\in\mathbb N\) with \(M=u^2s\), \(s\) squarefree. |
| **O7a** | Decide whether C7 has a deterministic polynomial-time algorithm. |
| **O7b** | Ask whether integer factoring C5 has a randomized polynomial-time reduction to C7. The exact printed subscript glyph is **UNVERIFIED** from OCR. |
| **C8 / O8** | Decide squarefreeness / ask whether that decision is in P. |
| **C9 / O9** | Compute the number of distinct prime factors / ask whether that function is in P. |

**KNOWN.** The FOCS text also records the exponent-multiset/divisor-count equivalence of Shallit–Shamir and reductions to \(d\), not a factoring equivalence for C9. It records no progress in its 1994 remarks on C8/C9. Those are historical statements, not evidence that no later result exists. [Adl]

**KNOWN metadata; exact ANTS statements UNVERIFIED.** Leonard M. Adleman, Kevin S. McCurley, *Open problems in number theoretic complexity, II*, ANTS-I, LNCS **877**, 291–322 (1994), DOI 10.1007/3-540-58691-1_70, is verified from the publisher and McCurley's publication list. Its accessible PDF has broken font encoding, and the screenshot interface did not expose readable page images in this session. C7/O7a/O7b are independently cited by the primary ESIGN report [ESIGN] (§3.1); I could not directly verify C8/O8/C9/O9 and their exact typography in the **ANTS** version. Do not treat verification of the FOCS list as direct verification of the ANTS list. [AM]

**KNOWN.** Heather Woll, *Reductions among number theoretic problems*, Information and Computation **72**(3), 167–179 (1987), DOI 10.1016/0890-5401(87)90030-7: the accessible abstract identifies a new **deterministic reduction from squarefreeness to Euler's totient**, along with an unrelated order/discrete-log reduction. This is the easy-to-hard direction. Exact internal theorem numbering was not accessible. [Woll]

**KNOWN, with partial verification.** Susan Landau, *Some remarks on computing the square parts of integers*, Information and Computation **78**(3), 246–253 (1988), DOI 10.1016/0890-5401(88)90028-4: the abstract confirms a deterministic polynomial-time Turing reduction of its largest-squarefree-factor function to Euler's totient and another reduction to Carmichael's function. The displayed definitions of its \(\pi,\theta\) were omitted by the accessible rendering, so the exact \(\pi/\theta\) notation and theorem numbers are **UNVERIFIED** here. Du–Volkovich, FSTTCS 2021, §1.1 independently confirms the Landau/Woll direction toward totient oracles. [Landau], [DV]

**KNOWN.** Johannes A. Buchmann, Hendrik W. Lenstra Jr., *Approximatting rings of integers in number fields* (the published title uses “Approximatting”), JTNB **6**(2), 221–260 (1994): **Theorem 1.3** gives polynomial transformations between constructing the ring of integers of a number field and finding the largest squarefree divisor of an integer. **Theorem 6.12** equates maximality recognition for an order with squarefreeness decision. **§6.11** explains that their construction improves Chistov's Turing reduction to a single-field transformation. Their remarks about obtaining maximal orders from discriminant data are carefully distinguished from factoring. [BL]

**KNOWN by the checked bibliography and exposition in [BL]; original theorem UNVERIFIED.** A. L. Chistov, *The complexity of constructing the ring of integers of a global field*, Doklady Akademii Nauk SSSR **306** (1989), 1063–1067; English translation Soviet Mathematics Doklady **39** (1989), 597–600. I verified the citation and the description of the reduction through [BL], not the original paper or its theorem number.

**KNOWN metadata; exact book locations UNVERIFIED.** Eric Bach, Jeffrey Shallit, *Algorithmic Number Theory, Volume I: Efficient Algorithms*, MIT Press (1996), has Chapter **14**, “Reductions Among Number-Theoretic Problems,” as verified on Shallit's book pages. I did not obtain the chapter text, so I do not assign theorem numbers or claim it proves an equivalence absent from the papers above. [BS]

**OPEN — separation search finding.** I found no ordinary relativized/oracle separation proving factoring strictly harder than the listed natural arithmetic functions. Generic-ring results, the restricted multiplicative-query argument, and the circuit lower bounds concern explicitly restricted models. They should not be reported as an unrestricted oracle separation between the actual integer problems.

### 3.6 Booker–Hiary–Keating: result and cost

**KNOWN.** Andrew R. Booker, Ghaith A. Hiary, Jon P. Keating, *Detecting squarefree numbers*, Duke Mathematical Journal **164**(2), 235–275 (2015), uses explicit formulas for quadratic Dirichlet \(L\)-functions. **Theorem 1.1**, under GRH, proves RSA-210, RSA-220, RSA-230, RSA-232 **not squarefull**: each has a prime factor of exponent one. It does **not** certify those four numbers squarefree. The paper reports months of desktop computation for RSA-210. [BHK]

**KNOWN / HEURISTIC.** **Proposition 2.3** gives, under quadratic GRH, an algorithm returning a nontrivial square factor or a conditional certificate of squarefreeness; for squarefree inputs its cost is \(\exp((\log M)^{\theta^*+o(1)})\), with \(\theta^*\) defined by extreme low-zero gaps among twists. The estimate \(\theta^*=2/3\) is **HEURISTIC**, based on additional random-matrix/extreme-gap reasoning, not a consequence of GRH alone. The introduction gives a rigorous elementary squarefreeness bound \(M^{1/6+o(1)}\) via Pollard–Strassen. None is polynomial in \(\log M\). [BHK]

### 3.7 Totient bits, class-number congruences, and \(\epsilon\)

**PROVED — exact tiny totient output.** For odd primes with \(N\equiv1\pmod4\), either both are \(1\pmod4\), giving \(16\mid(p-1)(q-1)\), or both are \(3\pmod4\), giving \(\varphi(N)\equiv4\pmod8\). Hence

\[
\epsilon=+1\iff\varphi(N)\equiv0\pmod8,
\qquad
\epsilon=-1\iff\varphi(N)\equiv4\pmod8.
\]

This is an exact reduction from \(\epsilon\) to three low totient bits, not a way to compute those bits from \(N\).

**KNOWN.** Genus theory gives \(r_2(\mathrm{Cl}(D))=\omega(D)-1\) for negative fundamental discriminants, and Rédei theory computes the \(4\)-rank from residue-symbol data at the prime-discriminant factors. One checked source is Étienne Fouvry, Jürgen Klüners, *On the 4-rank of class groups of quadratic number fields*, Inventiones Mathematicae **167**, 455–513 (2007), introduction and its algebraic discussion. The \(4\)-rank formula is also explicitly recorded in [Chan] (§1.4). [FK]

**PROVED — why \(h(-4N)\bmod8\) is not simply \(\epsilon\).** For \(N\equiv1\pmod4\), \(D=-4N\) is fundamental and has three ramified rational primes. Thus the \(2\)-rank is two, \(4\mid h(D)\), and \(h(D)\equiv4\pmod8\) iff its \(4\)-rank is zero. The Rédei matrix depends on more than \(\chi_4(p)\), including residue symbols between the factors. In fact the following balanced pairs show that neither possible residue determines \(\epsilon\), even after supplying \(N\bmod8\):

| \(p,q\) | \(N\bmod8\) | \(\epsilon\) | \(h(-4N)\) |
|---|---:|---:|---:|
| \(13,17\) | \(5\) | \(+1\) | \(16\) |
| \(7,11\) | \(5\) | \(-1\) | \(8\) |
| \(17,29\) | \(5\) | \(+1\) | \(12\) |
| \(19,23\) | \(5\) | \(-1\) | \(20\) |

**PROVED / reproducible finite verification, using KNOWN form reduction.** The class numbers in the table are obtained by enumerating primitive reduced positive binary quadratic forms: \(b^2-4ac=D\), \(-a<b\le a\le c\), with the usual positive-boundary convention and \(\gcd(a,b,c)=1\), and \(a\le\sqrt{|D|/3}\). The complete enumeration is implemented in `check_claims.py` beside this report. Each class has one reduced representative with this convention; the class-number computation by reduced-form enumeration is also documented in the [Magma handbook](https://magma.maths.usyd.edu.au/magma/handbook/text/374). This verifies the finite counterexamples, not a statistical independence assertion.

**OPEN — search finding.** These congruences can encode factor-dependent information in restricted families, but I found no polynomial-time procedure to compute the relevant class-number residue from a balanced semiprime that thereby predicts \(\epsilon\) beyond the residue floor, and no theorem that \(h(-4N)\bmod8\) generally determines \(\epsilon\). Computing the Rédei matrix from already known prime factors is not an algorithm for obtaining those factors or their bit.

## 4. PT's open classes: what follows, and what does not

### 4.1 Specify the probability space first

**PROVED — necessary formulation.** “Prediction advantage” is undefined without an instance distribution and a baseline. An adversarial distribution supported on \(p,q\equiv1\pmod4\) makes \(\epsilon\) constant. In this section, take independent uniform primes \(P,Q\in[X,2X]\), discard \(P=Q\), order them as \(p<q\), and, for \(\epsilon\), condition on \(pq\equiv1\pmod4\). A public modulus \(m=m(X)\) is fixed for the size class or chosen independently of the instance. The observation is only \(C=N\bmod m\), not all of \(N\). For a general residue target, advantage is relative to its best constant prediction, not necessarily \(1/2\).

**PROVED — exact floor versus probabilistic independence.** For every distribution and every \(F\),

\[
I((p\bmod m,q\bmod m);F(N\bmod m)\mid N\bmod m)=0.
\]

This tautological conditional-information statement is exact, but by itself it says nothing about whether the observed product residue predicts a factor residue. The latter needs equidistribution. If \(U,V\) are independent uniform in \((\mathbb Z/m\mathbb Z)^*\), then \(UV\) is independent of \(U\); the conditional factor pairs lie on \(\{(u,c/u)\}\). The same argument with \(\operatorname{lcm}(4,m)\) handles the \(N\equiv1\pmod4\) promise and \(\epsilon\).

### 4.2 A GRH theorem, with a stronger range than requested

**KNOWN — analytic input.** Hugh L. Montgomery, Robert C. Vaughan, *Multiplicative Number Theory I: Classical Theory*, Cambridge University Press (2007), **Theorem 13.7**, equations **(13.19)–(13.21)**, assumes GRH for characters modulo \(m\) and gives

\[
\psi(x,\chi)=E_0(\chi)x+O(\sqrt{x}\log x\log(mx)),
\]
\[
\pi(x,\chi)=E_0(\chi)\operatorname{li}(x)+O(\sqrt{x}\log(mx)).
\tag{GRH-primes}
\]

Here \(E_0\) is one for the principal character and zero otherwise. **Corollary 13.8** gives the corresponding pointwise progression estimate; the mean-square orthogonality calculation appears immediately afterward. The chapter text was checked directly. [MV]

**PROVED conditional on GRH — Theorem L.** Let \(m<X\), \(L=\operatorname{lcm}(4,m)<4X\), and suppose no prime in the sampling interval divides \(L\); this holds, in particular, if \(m<X\) and \(X>2\). Let \(\mu_L\) be the residue law of a uniform prime in \([X,2X]\), and \(u_L\) uniform on its unit group. Under GRH modulo \(L\),

\[
\|\mu_L-u_L\|_{TV}
\ll\sqrt{\frac{\varphi(L)}{X}}\,\log X\log(XL).
\tag{L}
\]

Every predictor \(F(N\bmod m)\) of \(\epsilon\), including an unbounded or randomized predictor, has advantage at most a constant times the right side plus \(O(\log X/X)\). In particular it is \(o(1)\) uniformly for \(m\le X^{1-\eta}\), every fixed \(\eta>0\), and even whenever \(\sqrt{m/X}\log^2X=o(1)\).

*Proof.* Let \(T=\pi(2X)-\pi(X)\asymp X/\log X\). For a nonprincipal character \(\chi\) modulo \(L\), (GRH-primes) gives

\[
\left|\widehat\mu_L(\chi)\right|
=\left|T^{-1}\sum_{X\le r\le2X\atop r\text{ prime}}\chi(r)\right|
\ll\frac{\log X\log(XL)}{\sqrt X}.
\]

Character orthogonality and Cauchy–Schwarz give

\[
\sum_{a\in(\mathbb Z/L)^*}|\mu_L(a)-1/\varphi(L)|^2
=\frac1{\varphi(L)}\sum_{\chi\ne\chi_0}|\widehat\mu_L(\chi)|^2,
\]

and hence (L). The product pair law differs from \(u_L\otimes u_L\) by at most twice that total variation. In the uniform model, for every \(c\), the first component \(U\) conditional on \(UV=c\) is uniform; \(\chi_4(U)\) is a fair sign independent of \(UV\). Conditioning on \(UV\equiv1\pmod4\), an event of probability \(1/2\), preserves this. Conditioning perturbs total variation by at most a constant factor when the actual event probability is bounded below. Discarding the diagonal changes it by \(O(1/T)\). On this promise \(\chi_4(P)=\chi_4(Q)\), so ordering the factors does not alter the target. A measurable predictor cannot increase total variation, proving the bound. ∎

**PROVED — the requested GRH range is valid in this model.** Since \(N\asymp X^2\), every \(m\le N^{1/4-o(1)}\) is well inside Theorem L's range. A pointwise-error proof would reach roughly \(m\ll\sqrt X\) with logarithmic losses; the mean-square / total-variation proof is stronger. There is no information-theoretic phase transition at \(N^{1/4}\) for observing the **product residue alone**.

**PROVED — general functions of the smaller factor residue.** For a fixed \(\eta>0\), the same \(o(1)\) independence conclusion holds for every \(f(p\bmod m)\) in the range \(m\le X^{1-\eta}\), relative to its best constant predictor. To accommodate the ordering, partition \([X,2X]\) into \(K\) consecutive equal intervals. For fixed \(K\), PNT and the same character bound show approximate uniformity of residues in each interval; the error is \(O(K\sqrt{m/X}\log^2X)\). On pairs in different intervals, the smaller factor is identified by the interval labels, independently of the approximately uniform residue pair. The same-interval probability is \(O(1/K)\). Let \(K\to\infty\) sufficiently slowly. Thus the joint law of \(p\bmod m\) and \(N\bmod m\) is \(o(1)\) from the independent uniform model. The proof also works after the \(N\equiv1\pmod4\) conditioning by adjoining modulus \(4\). ∎

**PROVED — limits.** This theorem concerns the natural independent prime law, not every distribution of balanced semiprimes, not individual fixed \(N\), not a modulus selected using arbitrary bits of the instance, and not multiple residues whose **combined** lcm exceeds the range. It says nothing about an algorithm allowed the high bits of \(N\) in addition to \(N\bmod m\).

### 4.3 Unconditional low bits: a verified route using Bourgain

**KNOWN.** Jean Bourgain, *Prescribing the binary digits of primes, II*, Israel Journal of Mathematics **206** (2015), 165–182, **unnumbered Theorem in §1**, equations **(1.1)–(1.2)**: for \(Z=2^b\), any set \(A\subset\{1,\ldots,b-1\}\) of \(r<\gamma b\) digit positions, some absolute \(\gamma>0\), and any assignments, the prime count is asymptotic to \(2^{-r}Z/\log Z\), with the odd least bit fixed. The paper and the allowed range were checked directly. [Bourgain]

**PROVED — unconditional consequence.** Take \(X=2^{b-1}\), \(m=2^j\), and prescribe the top bit and the \(j-1\) low nonconstant bits of a prime. There are \(r=j\) prescribed positions. For \(j\le\gamma' b\), any fixed \(0<\gamma'<\gamma\), the theorem gives uniformly

\[
\#\{r\in[X,2X]:r\text{ prime},\ r\equiv a\pmod{2^j}\}
=(1+o(1))\frac{T}{\varphi(2^j)}
\]

for every odd \(a\). Summing the relative errors yields total-variation distance \(o(1)\) from uniform. Apply the pair/fiber argument in Theorem L, adjoining the two mod-4 bits when necessary. Therefore any function of the low \(j\) bits of \(N\) has advantage \(o(1)\) for \(\epsilon\) for

\[
2^j\le N^c
\]

with some absolute \(c>0\) (choose any fixed \(2c<\gamma\) with slack for endpoints). This is an existence-of-constant range; I do not assign a numerical value to Bourgain's unspecified constant. Fixed leading-bit refinements also handle the interval bins used for the smaller-factor version. ∎

**KNOWN, but not used to claim an unchecked power-of-two theorem.** M. B. Barban, Yu. V. Linnik, N. G. Tshudakov, *On prime numbers in an arithmetic progression with a prime-power difference*, Acta Arithmetica **9** (1964), 375–390; P. X. Gallagher, *Primes in progressions to prime-power modulus*, Inventiones Mathematicae **16** (1972), 191–201. The checked paper of Stephan Baier and Sudhir Pujahari, *A Bombieri–Vinogradov-type theorem with prime power moduli*, Acta Arithmetica **204**(2) (2022), introduction, **(1.4)**, records the fixed **odd-prime** modulus asymptotic

\[
\pi(x;\rho^k,a)=\frac{\operatorname{li}(x)}{\varphi(\rho^k)}
\{1+O_{\rho,A,\eta}((\log x)^{-A})\},
\]

in the ranges \(\rho^k\le x^{3/8-\eta}\) for BLT and \(\rho^k\le x^{2/5-\eta}\) for Gallagher, with the cited later \(5/12\) improvement. Original theorem numbers are **UNVERIFIED** here. [BLT], [Gallagher], [BP]

**KNOWN / qualification.** Henryk Iwaniec, *On zeros of Dirichlet's L series*, Inventiones Mathematicae **23** (1974), 97–104, Theorems **1–2**, extends the zero-free-region work to general moduli controlled by their radical; powers of two are permitted. I did not reconstruct a precise prime-progression asymptotic with a specified exponent for \(2^j\) from those theorems in this report. Using Bourgain avoids silently extending the odd-prime formulations above to \(2\). [Iwaniec]

### 4.4 Why this does not contradict known-factor-bit attacks

**KNOWN.** Dan Boneh, Glenn Durfee, Yair Frankel, *An attack on RSA given a small fraction of the private key bits*, ASIACRYPT 1998, LNCS **1514**, 25–34, **Corollary 1, printed p. 28**, states that an \(n\)-bit RSA modulus can be factored in polynomial time given \(r\ge2^{n/4}\) and \(p\bmod r\). Its balance condition is **(1), p. 27**, \(\sqrt N/2<q<p<2\sqrt N\) (their factors are ordered oppositely to ours). If \(\gcd(r,N)>1\), factor by that gcd; otherwise their proof uses \(q\bmod r=N/p\bmod r\). It explicitly handles the boundary with a constant-bit exhaustive search. The claim was checked in a scan of the proceedings. Thus the conference citation is **Corollary 1**, not “Corollary 2.2”; whether the latter labels a different full version is **UNVERIFIED**. [BDF]

**KNOWN, with original numbering UNVERIFIED.** The underlying small-root method is Don Coppersmith, *Small solutions to polynomial equations, and low exponent RSA vulnerabilities*, Journal of Cryptology **10**(4), 233–260 (1997). Its original theorem numbers were not checked; the exact residue threshold needed here is independently stated in the checked BDF corollary. In particular \(\log_2r\ge n/4+O(\log n)\) is more than enough; polynomially many enumerated missing bits also preserve polynomial cost. [Coppersmith], [BDF]

**PROVED — no matching information threshold.** Knowing \(p\bmod M\) *together with the full integer \(N\)* and observing \(N\bmod M\) alone are different experiments. For the first, write \(p=a+Mt\); sufficiently large \(M\) makes \(t\) a small root modulo an unknown divisor of \(N\). For the second, the pair lies on a product fiber and \(a\) is not given. Theorems L and Bourgain's corollary concern the second experiment. The apparent meeting at \(N^{1/4}\) is therefore not a tight secrecy/recovery threshold.

### 4.5 What the circuit and automata papers say

| Reference | Checked result | Consequence for the present problem |
|---|---|---|
| Allender–Saks–Shparlinski, JCSS 2001 | **KNOWN:** Theorem **1**, \(\mathrm{TC}^0\subseteq\mathrm{AC}^0[\mathrm{Primes}],\mathrm{AC}^0[\mathrm{GCD}],\mathrm{AC}^0[\mathrm{Square\!-!Free}]\). Corollary **2** excludes exact computation in \(\mathrm{AC}^0[\ell]\) for every prime \(\ell\). [ASS] | **PROVED distinction:** no lower bound for \(\epsilon\) on the promised semiprime slice follows. |
| Anna Bernasconi, Igor Shparlinski, *Circuit complexity of testing square-free numbers*, STACS 1999, LNCS **1563**, 47–56 | **KNOWN:** superpolynomial size is needed for bounded-depth, unbounded-fanin squarefreeness circuits; checked author-affiliated abstract. [Bernasconi] | **PROVED distinction:** exact squarefreeness is a different target, constant on our semiprimes. |
| Ben Green, *On (not) computing the Möbius function using bounded depth circuits*, CPC **21**(6), 942–951 (2012) | **KNOWN:** **Theorem 1** gives \(2^{-b}\sum_{x<2^b}\mu(x)F(x)=o(1)\) for polynomial-size fixed-depth \(\mathrm{AC}^0\) families; the paper gives a quantitative bound. Its **Theorem 2** restates LMN. [Green] | **PROVED distinction:** \(\mu(pq)=1\), so this is not orthogonality to \(\chi_4(p)\) on semiprimes. Restricting a uniform average to a sparse slice is not licensed. |
| Konieczny–Lemańczyk–Müllner, *Multiplicative automatic sequences*, Math. Z. **300**, 1297–1318 (2022) | **KNOWN:** published **Theorem 1** (arXiv v2 **1.1**) writes \(a(n)=f_1(v_\rho(n))f_2(n/\rho^{v_\rho(n)})\), with fixed prime \(\rho\), eventually periodic \(f_1,f_2\), and specified multiplicativity/vanishing properties. [KLM] | **PROVED:** on sufficiently large semiprimes with \(p,q\ne\rho\), this is \(f_2(N)\), hence depends on \(N\) modulo a fixed period. The residue floor applies. It does not classify arbitrary automatic predictors. |
| Clemens Müllner, *Automatic sequences fulfill the Sarnak conjecture*, Duke Math. J. **166**(17), 3219–3290 (2017) | **KNOWN:** **Theorem 1.2** establishes Möbius disjointness for the automatic symbolic systems; **Theorem 1.3** gives prime-letter frequencies under its strong-connectivity and leading-zero hypotheses. [Müllner] | **PROVED distinction:** the stated Möbius theorem alone is not a semiprime factor-bit theorem. Its type-II machinery does give a more relevant restricted implication below. |
| Christian Mauduit, Joël Rivat, *Sur un problème de Gelfond: la somme des chiffres des nombres premiers*, Ann. Math. **171**(3), 1591–1646 (2010) | **KNOWN:** **Theorem 1** gives a power saving for \(\sum_{n\le x}\Lambda(n)e(\alpha s_b(n))\) when \((b-1)\alpha\notin\mathbb Z\); **Theorem 3** gives the digit-sum residue distribution, with the necessary congruence obstruction. [MR] | **PROVED distinction:** bits/digit sums of a prime are not the weighted digits of the product of two primes. Type-II estimates are the relevant ingredient, not an automatic transfer of the final theorem. |
| Bourgain, Israel J. Math. 2015 | **KNOWN:** the prescribed-bit count in §4.3. [Bourgain] | **PROVED:** enough for the low-bit product floor; not by itself a theorem about arbitrary \(\mathrm{AC}^0\) functions of all product bits. |

**PROVED — a nontrivial restricted extension using Müllner's type-II bound.** The checked arXiv v3 of [Müllner], §5.4.2, **Proposition 5.5**, equation **(29)**, bounds

\[
\sum_{u\asymp U}\sum_{v\asymp V}a_u b_v f(uv)e(\theta uv)
\]

uniformly for \(|a_u|,|b_v|\le1\), subject to its carry property (Definition 4.1), Fourier property (Definition 4.2), and balanced log-length condition (28). In particular, if a scalar-valued sequence \(f\) satisfies those hypotheses with linear Fourier decay \(\gamma(t)\ge c t\), its bound is \(O(X^{2-\kappa}(\log X)^C)\) when \(U,V\asymp X\), for some \(\kappa>0\).

Take \(a_u=1_{u\text{ prime}}\chi_4(u)\), \(b_v=1_{v\text{ prime}}\). This yields a negligible normalized correlation against \(\epsilon f(pq)\). For the promise \(pq\equiv1\pmod4\), use

\[
\chi_4(p)1_{\chi_4(p)=\chi_4(q)}=
\tfrac12(\chi_4(p)+\chi_4(q)).
\]

Each resulting sum has bounded separate coefficients of the permitted kind. Normalization by \(\asymp X^2/\log^2X\) preserves \(o(1)\). The excluded diagonal is negligible. Thus these carry/Fourier-decay classes admit a rigorously justified floor. This is a deduction from the stated type-II estimate, not a claim that every automaton output meets its hypotheses. In the absence of linear decay, merely \(\gamma(t)\to\infty\) is not enough for this normalized conclusion without checking the logarithmic losses. ∎

**OPEN / UNVERIFIED blanket extension.** I did not prove or verify a theorem covering *every* finite-automaton predictor of \(\epsilon\) on balanced semiprimes. Establishing the required decomposition with the correct weights and quantitative bounds may extend the preceding argument, but the Sarnak headline theorem alone is insufficient. I likewise found no full \(\mathrm{AC}^0\) floor for this target. This is a search/proof limitation, not a theorem that such extensions are impossible or absent from all literature.

### 4.6 What LMN would actually require

**KNOWN.** Nathan Linial, Yishay Mansour, Noam Nisan, *Constant depth circuits, Fourier transform, and learnability*, JACM **40**(3), 607–620 (1993), supplies Fourier concentration for \(\mathrm{AC}^0\). One checked form is [Green] (Theorem 2): a depth-\(d\), size-\(S\) Boolean circuit satisfies

\[
\sum_{|A|>t}|\widehat F(A)|^2\le2S\,2^{-t^{1/d}/20}.
\tag{LMN}
\]

**PROVED — sufficient correlation criterion with normalization.** Let \(\mathcal S\) be the semiprime slice in a \(b\)-bit cube, with uniform density \(\rho\), and let \(y\) be \(\epsilon\), centered by its mean on the slice. Suppose

\[
\beta_t=\max_{|A|\le t}
\left|\mathbb E_{N\in\mathcal S} y(N)(-1)^{\sum_{i\in A}N_i}\right|.
\]

If the Fourier tail of \(F\) above degree \(t\) is at most \(\eta\), then

\[
\left|\mathbb E_{\mathcal S} yF\right|
\le \sqrt{\sum_{j\le t}\binom bj}\,\beta_t
+2\sqrt{\eta/\rho}.
\tag{AC-criterion}
\]

*Proof.* Expand \(F=F_{\le t}+F_{>t}\). For the low part, Cauchy–Schwarz bounds \(\sum_{|A|\le t}|\widehat F(A)|\) by the square root of the number of terms, since \(\sum_A|\widehat F(A)|^2=1\). For the tail,

\[
\mathbb E_{\mathcal S}|F_{>t}|
\le \sqrt{\mathbb E_{\mathcal S}|F_{>t}|^2}
\le \sqrt{\eta/\rho},
\]

and \(|y|\le2\). This proves the criterion. The usual same-length dyadic semiprime slice has \(\rho\asymp1/\log^2(2^b)\), by counting distinct unordered prime pairs in the chosen dyadic interval. ∎

**PROVED — why merely \(o(1)\) for each parity is inadequate.** LMN needs \(t=\operatorname{polylog}(b)\), and the number of parities is \(\exp(O(t\log b))\). The criterion requires a bound beating that loss, as well as the slice-density loss. “All low-degree correlations tend to zero” without a uniform rate is insufficient.

**OPEN — weighted bilinear problem.** The needed parity sums can be expanded into additive Fourier sums with typical terms

\[
\sum_{p,q\asymp X}\chi_4(p)e(bpq/2^k),
\]

but controlling the expansion costs, major arcs, high digit positions, and uniform rates is essential. The low-bit theorems handle \(k\) in their modulus range. A parity of a few arbitrary bit positions may involve a denominator as large as the full product range; “few positions” does not mean a small low-bit modulus. I have not proved the uniform estimates needed for (AC-criterion). The type-II bounds in the cited digital papers suggest a route for some classes, but do not, as quoted, establish the full \(\mathrm{AC}^0\) claim.

### 4.7 Algorithms that use the integer structure

**PROVED — scope of PT.** In the supplied pure-tree model, a test integer is fixed by the public size parameters, independent coins and prior zero-test answers. If it vanishes modulo exactly one of \(p,q\), its gcd splits \(N\); otherwise the composite run follows their common prime transcript. That is the partition mechanism. It does not cover unrestricted uses of the binary expansion, \(N\)-dependent exponents, chosen integer representatives, continued fractions of \(\sqrt N\), forms of discriminant \(-4N\), or lattice constructions depending on \(N\).

**OPEN — unrestricted leakage search finding.** I found no general theorem excluding factor-bit prediction by those \(N\)-dependent methods, and no generic amplification from their weak fixed-bit prediction to factoring. Some such algorithms actually factor on special subfamilies; §5 gives an explicit example. A proof about a generic-ring model should not be extended to a program allowed to read all integer bits.

## 5. A positive leak, its size, and an artefact check

### 5.1 One Fermat test gives an exponentially tiny real advantage

**PROVED.** On an odd semiprime \(N=pq\), compute \(a=\lceil\sqrt N\rceil\). If \(a^2-N\) is a square \(b^2\), test \(a-b\) and \(a+b\) as factors. The computation has polynomial bit complexity. If

\[
q-p\le\sqrt X,\qquad p,q\in[X,2X],
\]

then this one test succeeds. Indeed,

\[
\frac{p+q}2-\sqrt{pq}
=\frac{(q-p)^2}{2(\sqrt p+\sqrt q)^2}
\le\frac18<1.
\]

The left endpoint \((p+q)/2\) is an integer strictly above \(\sqrt N\), so it is \(a\), and \(b=(q-p)/2\).

**PROVED — unconditional density lower bound for independent prime sampling.** Partition \([X,2X]\) into \(B=O(\sqrt X)\) bins of length at most \(\sqrt X\). For each residue class \(r=1,3\pmod4\), let \(t_{i,r}\) be the number of primes in bin \(i\), and \(T_r=\sum_i t_{i,r}\asymp X/(2\log X)\), by PNT in the fixed progressions. The number of same-class, distinct unordered pairs within bins is

\[
\sum_i\binom{t_{i,r}}2
\ge\frac12\left(T_r^2/B-T_r\right)
\gg T_r^2/\sqrt X.
\]

All such pairs satisfy the one-test condition. Dividing by the total number \(\asymp X^2/\log^2X\) of same-class pairs shows a factoring success probability \(\Omega(X^{-1/2})=\Omega(N^{-1/4})\), under the independently sampled distribution conditioned on \(N\equiv1\pmod4\). No prime-pair conjecture or short-interval lower bound is used: Cauchy–Schwarz over the bins suffices. ∎

**PROVED — predictor and the baseline qualification.** For an exactly unbiased target, use the standard balanced mixture: choose a fair hidden sign; choose two independent uniform primes from \([X,2X]\) in class \(1\bmod4\) for sign \(+1\), or in class \(3\bmod4\) for sign \(-1\); discard equality. The density proof applies separately to each population, so its average factoring success is still \(\Omega(X^{-1/2})\). Run the test; on success compute \(\epsilon\) from the recovered factor; otherwise output an independent fair sign. If success probability is \(s\), accuracy is exactly \(1/2+s/2\), exceeding the optimal constant by \(\Omega(N^{-1/4})=2^{-\Theta(n)}\).

**PROVED — finite-prime bias is not being hidden.** For uniformly sampled primes conditioned on \(N\equiv1\bmod4\), without the fair-mixture adjustment, the same algorithm still has accuracy \(1/2+s/2\). Its improvement over that distribution's *best constant* need not follow from this calculation, because the two finite prime-class counts can differ. I do not claim such an improvement for that law. The fair-mixture formulation removes this baseline issue explicitly. The advantage remains exponentially small, with no amplification claim.

**PROVED — artefact check.** Within each stated residue class the two primes are sampled independently and uniformly; no closeness is imposed. The class-choice label is hidden, and every non-random answer is supported by recovered and verified factors. The advantage does not arise from a supplied label, a near-prime generation rule, variable prime lengths, training/test overlap, or a fitted classifier. It is leakage **through factoring**, not a new nonsplitting statistic. Allowing polynomially many Fermat steps enlarges the near-square success set but does not change that distinction.

### 5.2 Other positive information must be charged to its input

**PROVED.** A radical value for \(p^2q\), a totient residue modulo \(8\), an exact torus order, or selected elliptic twist counts can all reveal information as shown above. They are extra oracle outputs. Treating their arithmetic definitions as polynomial-time procedures would be a leakage artefact. Similarly, observing a known Legendre symbol of \(p\) is extra factor information; computing a Jacobi symbol of \(N\) only gives a product of local symbols.

**OPEN — search finding.** This investigation produced no new inverse-polynomial-advantage, nonsplitting, polynomial-time predictor of \(\epsilon\) on independently generated balanced semiprimes.

## 6. Verification limits and unresolved work

**UNVERIFIED — exact source details.** The following are deliberately not promoted to verified literature facts:

1. Exact C8/O8/C9/O9 numbering and literal wording in the **ANTS** Adleman–McCurley paper. Its metadata and C7 citation are confirmed, and the FOCS C7–C9 list is readable, but the ANTS PDF extraction is corrupted.
2. The printed randomized-reduction subscript in FOCS O7b. The meaning is supported by its definitions; the OCR glyph is ambiguous. No literal transcription is supplied.
3. Internal theorem numbers and full technical hypotheses in Woll and Landau beyond the reductions explicitly stated in the accessible abstracts. In particular the rendered definitions of Landau's \(\pi,\theta\) were missing.
4. Chistov's original theorem number/text, and exact theorem/page references within Bach–Shallit Chapter 14. Buchmann–Lenstra Theorems 1.3 and 6.12 were directly checked instead.
5. Original BLT/Gallagher theorem numbers, and a specific numerical power-of-two progression exponent deduced directly from their work or Iwaniec's extension. The odd-prime statements are checked through Baier–Pujahari's explicit historical statement; the unconditional powers-of-two result proved here uses Bourgain.
6. Coppersmith's original theorem numbers and whether a separate BDF version labels its residue result “Corollary 2.2.” The conference **Corollary 1** and its exact hypotheses are verified.
7. ACGS theorem numbers. The RSA/Rabin-message scope is checked; Håstad–Näslund 8.1/9.1 are verified.
8. A full finite-automaton or \(\mathrm{AC}^0\) lower bound for predicting \(\epsilon\) on the specified semiprime distribution. The restricted carry/Fourier-decay implication and the LMN sufficient criterion are given, but their missing hypotheses are not assumed away.
9. An unrestricted oracle separation for any natural bit in §3, or a modulus-changing random self-reduction for factoring/Blum factoring. The report records **no result found**, not nonexistence.
10. Any universal theorem excluding modular/nonalgebraic factor movement or large unknown cofactors, or permitting arbitrary approximate algebraic identities. The proofs here do not cover them.

**PROVED — audit boundary.** Finite computational checks accompanying this report are sanity checks of identities and examples, not evidence for asymptotic lower bounds. The principal proofs are H, R, the oracle decompositions, the residue total-variation argument, the LMN criterion, and the Fermat density argument. No experiment is used to settle an open complexity question.

**PROVED — reproducibility record.** Running **python3 check_claims.py** passed 750 finite instances of Theorem H's divisibility and size assertions; radical/parity-kernel and one-call Möbius identities for integers 1–1000; divisor-oracle signature interpolation for 1–200; all four displayed class-number enumerations; the exact \(3/8\) twist-sign probability; and the stated one-test Fermat sufficient condition on dyadic prime pairs at \(X=256,512,1024\). Full output is in **check_claims_output.txt**. These checks are finite verification only.

## References and checked source links

The identifiers below link to the paper, publisher, or author-hosted copy used. Bibliographic entries with partial verification are qualified in §6. Search-result dates were not used as publication dates.

* **[ACGS]** Werner Alexi, Benny Chor, Oded Goldreich, Claus P. Schnorr, *RSA and Rabin functions: Certain parts are as hard as the whole*, SIAM Journal on Computing **17**(2), 194–209 (1988). DOI [10.1137/0217013](https://epubs.siam.org/doi/abs/10.1137/0217013).
* **[HN]** Johan Håstad, Mats Näslund, *The security of all RSA and discrete log bits*, Journal of the ACM **51**(2), 187–230 (2004). [Author manuscript](https://people.kth.se/~johanh/hnrsaacm.pdf), §§8–9, Theorems 8.1 and 9.1.
* **[LLL]** Arjen K. Lenstra, Hendrik W. Lenstra Jr., László Lovász, *Factoring polynomials with rational coefficients*, Mathematische Annalen **261**, 515–534 (1982). DOI 10.1007/BF01457454. [Paper copy](https://www.math.ucdavis.edu/~deloera/MISC/LA-BIBLIO/trunk/Lovasz/LovaszLenstraLenstrafactor.pdf).
* **[AFK]** Martín Abadi, Joan Feigenbaum, Joe Kilian, *On hiding information from an oracle*, Journal of Computer and System Sciences **39**(1), 21–50 (1989). DOI 10.1016/0022-0000(89)90018-4. [Author-hosted paper](https://www.cs.yale.edu/homes/jf/AFK-JCSS89.pdf), §6 Question 1.
* **[FF]** Joan Feigenbaum, Lance Fortnow, *On the random-self-reducibility of complete sets*, SIAM Journal on Computing **22**(5), 994–1005 (1993). DOI 10.1137/0222061. [Author-hosted paper](https://lance.fortnow.com/papers/files/rsr.pdf).
* **[BF]** Donald Beaver, Joan Feigenbaum, *Hiding instances in multioracle queries*, STACS 1990, LNCS **415**, 37–48 (1990). DOI [10.1007/3-540-52282-4_30](https://doi.org/10.1007/3-540-52282-4_30). The general-construction description is additionally checked in [FF]'s introduction.
* **[BT]** Andrej Bogdanov, Luca Trevisan, *On worst-case to average-case reductions for NP problems*, FOCS 2003, 308–317; SIAM Journal on Computing **36**(4), 1119–1159 (2006). [Journal abstract](https://epubs.siam.org/doi/10.1137/S0097539705446974).
* **[Rabin]** Michael O. Rabin, *Digitalized signatures and public-key functions as intractable as factorization*, MIT/LCS/TR-212 (1979). [MIT repository](https://dspace.mit.edu/entities/publication/15016449-ca26-4d45-82ea-7c200d82d777).
* **[AM]** Leonard M. Adleman, Kevin S. McCurley, *Open problems in number theoretic complexity, II*, ANTS-I, LNCS **877**, 291–322 (1994). [Publisher](https://link.springer.com/chapter/10.1007/3-540-58691-1_70); [PDF with unreadable extraction](https://www.cse.iitk.ac.in/users/manindra/CS681/open-problems-in-number.pdf).
* **[Adl]** Leonard M. Adleman, *Algorithmic number theory—the complexity contribution*, 35th FOCS, 88–113 (1994). DOI 10.1109/SFCS.1994.365702. [Readable paper reproduction](https://paperzz.com/doc/9179771/algorithmic-number-theory-the-complexity-contribution), printed pp. 98–99.
* **[Woll]** Heather Woll, *Reductions among number theoretic problems*, Information and Computation **72**(3), 167–179 (1987). [Publisher](https://www.sciencedirect.com/science/article/pii/0890540187900307), DOI 10.1016/0890-5401(87)90030-7; abstract accessible through search, full-text fetch blocked.
* **[Landau]** Susan Landau, *Some remarks on computing the square parts of integers*, Information and Computation **78**(3), 246–253 (1988). [Publisher abstract](https://www.sciencedirect.com/science/article/pii/0890540188900284), DOI 10.1016/0890-5401(88)90028-4.
* **[SS]** Jeffrey Shallit, Adi Shamir, *Number-theoretic functions which are equivalent to number of divisors*, Information Processing Letters **20**(3), 151–153 (1985). DOI 10.1016/0020-0190(85)90084-5. [Author-affiliated record and abstract](https://weizmann.elsevierpure.com/en/publications/number-theoretic-functions-which-are-equivalent-to-number-of-divi/).
* **[BL]** Johannes A. Buchmann, Hendrik W. Lenstra Jr., *Approximatting rings of integers in number fields*, Journal de Théorie des Nombres de Bordeaux **6**(2), 221–260 (1994). DOI 10.5802/jtnb.113. [Full paper](https://www.numdam.org/item/JTNB_1994__6_2_221_0.pdf), Theorems 1.3, 6.12, §6.11 and bibliography [6] for Chistov.
* **[BS]** Eric Bach, Jeffrey Shallit, *Algorithmic Number Theory, Volume I: Efficient Algorithms*, MIT Press (1996). [Author book page](https://cs.uwaterloo.ca/~shallit/ant.html); [contents](https://cs.uwaterloo.ca/~shallit/antdesc.html), Chapter 14.
* **[DV]** Yang Du, Ilya Volkovich, *Approximating the number of prime factors given an oracle to Euler's totient function*, FSTTCS 2021, LIPIcs **213**, article 17, 17:1–17:10 (2021). [Paper](https://drops.dagstuhl.de/storage/00lipics/lipics-vol213-fsttcs2021/LIPIcs.FSTTCS.2021.17/LIPIcs.FSTTCS.2021.17.pdf), §1.1.
* **[ESIGN]** *Security of \(N=P^2Q\) in ESIGN*, CRYPTREC external report, dated 28 January 2002; no individual author identified on the checked title page. [Report](https://www.cryptrec.go.jp/exreport/cryptrec-ex-1010-2001.pdf), §3.1. Used only to corroborate the ANTS C7/O7a/O7b citation.
* **[BHK]** Andrew R. Booker, Ghaith A. Hiary, Jon P. Keating, *Detecting squarefree numbers*, Duke Mathematical Journal **164**(2), 235–275 (2015). DOI 10.1215/00127094-2856619. [arXiv v2](https://arxiv.org/pdf/1304.6937), Theorem 1.1, Proposition 2.3, §3.
* **[GM]** Shafi Goldwasser, Silvio Micali, *Probabilistic encryption*, Journal of Computer and System Sciences **28**(2), 270–299 (1984). [Publisher](https://www.sciencedirect.com/science/article/pii/0022000084900709), DOI 10.1016/0022-0000(84)90070-9.
* **[CMS]** Christian Cachin, Silvio Micali, Markus Stadler, *Computationally private information retrieval with polylogarithmic communication*, EUROCRYPT 1999, LNCS **1592**, 402–414 (1999). DOI 10.1007/3-540-48910-X_28. [Author-hosted paper](https://people.csail.mit.edu/silvio/Selected%20Scientific%20Papers/Private%20Information%20Retrieval/Computationally%20Private%20Information%20Retrieval%20with%20Polylogarithmic%20Communication.pdf).
* **[Seurin]** Yannick Seurin, *On the lossiness of the Rabin trapdoor function*, PKC 2014, LNCS **8383**, 380–398 (2014). DOI 10.1007/978-3-642-54631-0_22. [Full version, ePrint 2013/256](https://eprint.iacr.org/2013/256.pdf), dated 15 April 2014.
* **[FK]** Étienne Fouvry, Jürgen Klüners, *On the 4-rank of class groups of quadratic number fields*, Inventiones Mathematicae **167**, 455–513 (2007; online 2006). DOI 10.1007/s00222-006-0021-2. [Author-hosted manuscript](https://math.uni-paderborn.de/fileadmin/mathematik/AG-Computeralgebra/Publications-klueners/ranks.pdf).
* **[Chan]** Yik Tung Chan, *On the 2-part of class groups and Diophantine equations*, PhD thesis, University College London (10 August 2020), §1.4. [Thesis](https://discovery.ucl.ac.uk/id/eprint/10107372/1/thesis_chan_ucl.pdf). Used for its explicitly displayed Rédei rank formula; original Rédei theorem numbering not verified.
* **[MV]** Hugh L. Montgomery, Robert C. Vaughan, *Multiplicative Number Theory I: Classical Theory*, Cambridge University Press (2007). [Chapter 13](https://personal.science.psu.edu/rcv4/personal/Publications/MNTI/17.0_pp_419_451_Conditional_estimates.pdf), Theorem 13.7, (13.19)–(13.21), Corollary 13.8.
* **[BLT]** M. B. Barban, Yu. V. Linnik, N. G. Tshudakov, *On prime numbers in an arithmetic progression with a prime-power difference*, Acta Arithmetica **9** (1964), 375–390. [Publisher](https://www.impan.pl/en/publishing-house/journals-and-series/acta-arithmetica/all/en/publishing-house/journals-and-series/acta-arithmetica/all/9/4/95481/on-prime-numbers-in-an-arithmetic-progression-with-a-prime-power-difference), DOI 10.4064/aa-9-4-375-390.
* **[Gallagher]** Patrick X. Gallagher, *Primes in progressions to prime-power modulus*, Inventiones Mathematicae **16** (1972), 191–201. DOI [10.1007/BF01425492](https://doi.org/10.1007/BF01425492). The quoted progression range is checked through [BP].
* **[BP]** Stephan Baier, Sudhir Pujahari, *A Bombieri–Vinogradov-type theorem with prime power moduli*, Acta Arithmetica **204**(2) (2022). [Checked manuscript v5](https://arxiv.org/pdf/2107.04348), introduction, (1.4), and historical discussion.
* **[Iwaniec]** Henryk Iwaniec, *On zeros of Dirichlet's L series*, Inventiones Mathematicae **23** (1974), 97–104, Theorems 1–2. [Digitized journal volume](https://gdz.sub.uni-goettingen.de/download/pdf/PPN356556735_0023/PPN356556735_0023.pdf).
* **[Bourgain]** Jean Bourgain, *Prescribing the binary digits of primes, II*, Israel Journal of Mathematics **206** (2015), 165–182. [Paper](https://arxiv.org/pdf/1307.0398), unnumbered Theorem, §1, (1.1)–(1.2).
* **[Coppersmith]** Don Coppersmith, *Small solutions to polynomial equations, and low exponent RSA vulnerabilities*, Journal of Cryptology **10**(4), 233–260 (1997). DOI 10.1007/s001459900030. [Author-affiliated record](https://research.ibm.com/publications/small-solutions-to-polynomial-equations-and-low-exponent-rsa-vulnerabilities).
* **[BDF]** Dan Boneh, Glenn Durfee, Yair Frankel, *An attack on RSA given a small fraction of the private key bits*, ASIACRYPT 1998, LNCS **1514**, 25–34 (1998), Corollary 1, p. 28. [Author page](https://cs.stanford.edu/people/dabo/pubs/abstracts/bits_of_d.html); checked against the [proceedings scan](https://www.nzdr.ru/data/media/biblio/kolxoz/Cs/CsLn/Advances%20in%20Cryptology%20-%20ASIACRYPT%20%2798%28LNCS1514%2C%20Springer%2C%201998%29%28ISBN%203540651098%29%28434s%29.pdf).
* **[ASS]** Eric Allender, Michael Saks, Igor Shparlinski, *A lower bound for primality*, Journal of Computer and System Sciences **62**(2), 356–366 (2001). DOI 10.1006/jcss.2000.1725. [Revised ECCC paper](https://eccc.weizmann.ac.il/report/1999/010/comment/1/download/), Theorem 1 and Corollary 2.
* **[Bernasconi]** Anna Bernasconi, Igor Shparlinski, *Circuit complexity of testing square-free numbers*, STACS 1999, LNCS **1563**, 47–56 (1999). DOI 10.1007/3-540-49116-3_4. [Author-affiliated abstract](https://researchers.mq.edu.au/en/publications/circuit-complexity-of-testing-square-free-numbers/).
* **[Green]** Ben Green, *On (not) computing the Möbius function using bounded depth circuits*, Combinatorics, Probability and Computing **21**(6), 942–951 (2012). DOI 10.1017/S0963548312000284. [Paper](https://arxiv.org/pdf/1103.4991), Theorems 1–2.
* **[KLM]** Jakub Konieczny, Mariusz Lemańczyk, Clemens Müllner, *Multiplicative automatic sequences*, Mathematische Zeitschrift **300**, 1297–1318 (2022). [Published Theorem 1](https://link.springer.com/article/10.1007/s00209-021-02834-3); arXiv v2 labels it Theorem 1.1.
* **[Müllner]** Clemens Müllner, *Automatic sequences fulfill the Sarnak conjecture*, Duke Mathematical Journal **166**(17), 3219–3290 (2017). DOI 10.1215/00127094-2017-0024. [arXiv v3](https://arxiv.org/pdf/1602.03042), Theorems 1.2–1.3 and Proposition 5.5, (29).
* **[MR]** Christian Mauduit, Joël Rivat, *Sur un problème de Gelfond: la somme des chiffres des nombres premiers*, Annals of Mathematics **171**(3), 1591–1646 (2010). [Journal PDF](https://annals.math.princeton.edu/wp-content/uploads/annals-v171-n3-p04-p.pdf), Theorems 1 and 3.
* **[LMN]** Nathan Linial, Yishay Mansour, Noam Nisan, *Constant depth circuits, Fourier transform, and learnability*, Journal of the ACM **40**(3), 607–620 (1993). DOI [10.1145/174130.174138](https://doi.org/10.1145/174130.174138); the quantitative form quoted here is checked in [Green] (Theorem 2).
* **[Pratt]** Vaughan R. Pratt, *Every prime has a succinct certificate*, SIAM Journal on Computing **4**(3), 214–220 (1975), Theorems 1–2. [Publisher text](https://epubs.siam.org/doi/10.1137/0204018).
* **[Schoof]** René Schoof, *Elliptic curves over finite fields and the computation of square roots mod p*, Mathematics of Computation **44**(170), 483–494 (1985). [Paper copy](https://pages.cs.wisc.edu/~cs812-1/Schoof85.pdf), abstract and introduction.

[ACGS]: https://epubs.siam.org/doi/abs/10.1137/0217013
[HN]: https://people.kth.se/~johanh/hnrsaacm.pdf
[LLL]: https://www.math.ucdavis.edu/~deloera/MISC/LA-BIBLIO/trunk/Lovasz/LovaszLenstraLenstrafactor.pdf
[AFK]: https://www.cs.yale.edu/homes/jf/AFK-JCSS89.pdf
[FF]: https://lance.fortnow.com/papers/files/rsr.pdf
[BF]: https://doi.org/10.1007/3-540-52282-4_30
[BT]: https://epubs.siam.org/doi/10.1137/S0097539705446974
[Rabin]: https://dspace.mit.edu/entities/publication/15016449-ca26-4d45-82ea-7c200d82d777
[AM]: https://link.springer.com/chapter/10.1007/3-540-58691-1_70
[Adl]: https://paperzz.com/doc/9179771/algorithmic-number-theory-the-complexity-contribution
[Woll]: https://www.sciencedirect.com/science/article/pii/0890540187900307
[Landau]: https://www.sciencedirect.com/science/article/pii/0890540188900284
[SS]: https://weizmann.elsevierpure.com/en/publications/number-theoretic-functions-which-are-equivalent-to-number-of-divi/
[BL]: https://www.numdam.org/item/JTNB_1994__6_2_221_0.pdf
[BS]: https://cs.uwaterloo.ca/~shallit/antdesc.html
[DV]: https://drops.dagstuhl.de/storage/00lipics/lipics-vol213-fsttcs2021/LIPIcs.FSTTCS.2021.17/LIPIcs.FSTTCS.2021.17.pdf
[ESIGN]: https://www.cryptrec.go.jp/exreport/cryptrec-ex-1010-2001.pdf
[BHK]: https://arxiv.org/pdf/1304.6937
[GM]: https://www.sciencedirect.com/science/article/pii/0022000084900709
[CMS]: https://people.csail.mit.edu/silvio/Selected%20Scientific%20Papers/Private%20Information%20Retrieval/Computationally%20Private%20Information%20Retrieval%20with%20Polylogarithmic%20Communication.pdf
[Seurin]: https://eprint.iacr.org/2013/256.pdf
[FK]: https://math.uni-paderborn.de/fileadmin/mathematik/AG-Computeralgebra/Publications-klueners/ranks.pdf
[Chan]: https://discovery.ucl.ac.uk/id/eprint/10107372/1/thesis_chan_ucl.pdf
[MV]: https://personal.science.psu.edu/rcv4/personal/Publications/MNTI/17.0_pp_419_451_Conditional_estimates.pdf
[BLT]: https://doi.org/10.4064/aa-9-4-375-390
[Gallagher]: https://doi.org/10.1007/BF01425492
[BP]: https://arxiv.org/pdf/2107.04348
[Iwaniec]: https://gdz.sub.uni-goettingen.de/download/pdf/PPN356556735_0023/PPN356556735_0023.pdf
[Bourgain]: https://arxiv.org/pdf/1307.0398
[Coppersmith]: https://research.ibm.com/publications/small-solutions-to-polynomial-equations-and-low-exponent-rsa-vulnerabilities
[BDF]: https://cs.stanford.edu/people/dabo/pubs/abstracts/bits_of_d.html
[ASS]: https://eccc.weizmann.ac.il/report/1999/010/comment/1/download/
[Bernasconi]: https://researchers.mq.edu.au/en/publications/circuit-complexity-of-testing-square-free-numbers/
[Green]: https://arxiv.org/pdf/1103.4991
[KLM]: https://link.springer.com/article/10.1007/s00209-021-02834-3
[Müllner]: https://arxiv.org/pdf/1602.03042
[MR]: https://annals.math.princeton.edu/wp-content/uploads/annals-v171-n3-p04-p.pdf
[LMN]: https://doi.org/10.1145/174130.174138
[Pratt]: https://epubs.siam.org/doi/10.1137/0204018
[Schoof]: https://pages.cs.wisc.edu/~cs812-1/Schoof85.pdf
