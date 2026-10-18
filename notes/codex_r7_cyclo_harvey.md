> Source: GPT-6.1-sol round 7 (codex/r7_cyclo_harvey.md → .out), 2026-10-18. Evaluates user's cyclotomic-Frobenius/Harvey idea and corrections. Curated copy.

The Frobenius workaround is valid. **A speedup from it is UNVERIFIED.** It supplies powers with hidden exponents, but neither shorter search intervals nor sufficiently many one-sided collisions follow automatically. I read both project files; no files were changed.

**1. What is computable — proved.** Assume \(N=pq\), with distinct primes, and first check \(\gcd(m,N)=1\). In
\[
R=(\mathbb Z/N\mathbb Z)[X]/(\Phi_m(X)),\qquad D=\varphi(m),
\]
every \(r\in(\mathbb Z/m\mathbb Z)^*\) defines the computable automorphism
\(\sigma_r:X\mapsto X^r\). For the correct guess \(r\equiv p\pmod m\), its reduction modulo \(p\) is Frobenius on the entire algebra, which may be a product of fields.

Consequently, for a unit \(\alpha\),
\[
\alpha^{p^k}\equiv\sigma_r^k(\alpha),\qquad
\alpha^{q^k}\equiv\sigma_r^{-k}(\alpha^{N^k})\pmod p.
\]
More generally, for \(F(P,Q)=\sum c_{ij}P^iQ^j\), the computable element
\[
T_{r,F}=\prod_{i,j}\sigma_r^{\,i-j}(\alpha)^{c_{ij}N^j}
\]
satisfies \(T_{r,F}\equiv\alpha^{F(p,q)}\pmod p\). Negative coefficients require inverses.

Thus the proposed handles exist. For fixed degree and coefficients of polynomial bit length, powering takes polynomially many ring operations. The costs are the \(D\) residue guesses, arithmetic on \(D\) coefficients, and constructing elements with useful orders. Using \(\alpha=\zeta_m\) alone merely reproduces the guessed residues.

There are two important qualifications:

- **Proved:** modulo \(\Phi_d(p)\), the powers satisfy the polynomial relation \(\Phi_d(p)=0\); they are not \(d\) independent handles. In the \(p^2+1\) subgroup, \(p^2\equiv-1\), so \(p^3\equiv-p\) and \(q^2\equiv-N^2\).
- **Proved:** the subgroup can sometimes be constructed without knowing \(p\). If \(\sigma\) has order four, set \(\alpha=\sigma^2(\beta)/\beta\). For a correct Frobenius guess, \(\alpha^{p^2+1}=1\pmod p\). But obtaining degree-four residue fields for every possible factor is an additional obligation; guessing residues does not guarantee that splitting type.

This is established cyclotomic machinery. **It is not exactly the same algorithm as APR-CL:** APR-CL uses special Gauss/Jacobi-sum identities to restrict possible divisor residues, typically to powers of \(N\) modulo an auxiliary modulus. See Cohen–Lenstra, *Primality Testing and Jacobi Sums*, **Math. Comp. 42 (1984), 297–330**, and Lenstra, *Galois theory and primality testing*, **LNM 1142 (1985), 169–189**. ([Original APR-CL paper](https://pub.math.leidenuniv.nl/~lenstrahw/PUBLICATIONS/1984a/art.pdf), [Lenstra chapter](https://doi.org/10.1007/BFb0074800))

The residue-class decoding step is separate: Lenstra’s 1984 algorithm handles modulus \(M>N^{1/3}\); Coppersmith–Howgrave-Graham–Nagaraj construct all divisors in a specified class for \(M\ge N^{1/4+\epsilon}\), fixed \(\epsilon>0\). These results **consume residue information; they do not obtain the correct residue for free**. ([Lenstra, **Math. Comp. 42 (1984), 331–340**](https://doi.org/10.1090/S0025-5718-1984-0726007-1), [CHGN, **77 (2008), 531–545**](https://doi.org/10.1090/S0025-5718-07-02007-8))

**2. Which identities can produce factors? — proved, with explicit limitations.**

The identity
\[
\sigma(x^N)=\sigma(x)^N
\]
holds globally and supplies no separating defect.

By contrast, searching integer candidates \(u\) for
\[
z=\alpha^u-T_{r,F}
\]
can give \(z=0\pmod p\) at \(u=F(p,q)\). A gcd of its coefficient vector with \(N\) returns \(p\) precisely when \(z\ne0\pmod q\).

That second condition is substantive. For example, if \(p\equiv q\pmod m\), the guessed \(\sigma\) is Frobenius at both primes, and
\[
T_{r,F}\equiv\alpha^{F(q,p)}\pmod q.
\]
Writing \(L_q=\operatorname{ord}(\alpha\bmod q)\), the collision is one-sided exactly when
\[
L_q\nmid F(p,q)-F(q,p).
\]
For \(F(P,Q)=AQ+BP\), this becomes
\(L_q\nmid(A-B)(q-p)\). For symmetric \(F\), it automatically fails: the collision holds at both primes. Such a collision may still identify a useful integer candidate, but it cannot be counted as a separating gcd event.

Two concrete cost checks expose the missing speedup:

- **Proved operation count, conditional on suitable element orders:** search \(\alpha^p=\sigma_r(\alpha)\), writing \(p=r+mt\). For balanced factors, \(t\) ranges over \(L=O(\sqrt N/m)\) integers. Sharing \(B\) baby steps across all \(D\) guesses gives \(B+DL/B\) list elements. Even granting fast detection of collisions modulo an unknown factor, dense algebra arithmetic costs approximately
  \[
  \widetilde O\!\left(D\sqrt{DL}+D^2\right).
  \]
  With \(m=N^\mu\), \(D=N^{\mu+o(1)}\), this is
  \[
  N^{\max(1/4+\mu,\,2\mu)+o(1)}.
  \]
  This particular scheme reaches exponent \(1/4\) at constant conductor, before proving a deterministic order-selection guarantee. Ordinary hashing detects equality modulo \(N\), so it cannot replace the unknown-factor collision machinery.

- **Proved failure of a concrete higher-order cancellation:** let \(x=Aq,\ y=Bp,\ K=ABN,\ s=x+y\). The homogeneous cubic
  \[
  P=x^3+y^3-9K(x+y)=s^3-12Ks
  \]
  has fourth-order flatness near \(x=y\). Nevertheless, with \(\delta=s-2\sqrt K\),
  \[
  P+16K^{3/2}=\delta^2(s+4\sqrt K).
  \]
  Its centered deviation divided by the original deviation \(\delta\) equals
  \[
  (x-y)^2\frac{s+4\sqrt K}{s+2\sqrt K}\ge1
  \]
  whenever the positive integers \(x,y\) differ. It cancels more Taylor terms while **widening** the integer exponent interval. For \(A=B=1\), the respective scales are \((q-p)^4/\sqrt N\) and \((q-p)^2/\sqrt N\).

**Conjectured:** some other polynomial family and batching method could improve this balance. **UNVERIFIED:** any resulting exponent below \(1/5\), including a proof of coverage, coefficient bounds, element orders, and collision extraction. The calculations above are scoped failures, not a general lower bound against extension methods.

**3. Prior literature located.**

| Primary source | Relevance |
|---|---|
| Bach–Shallit, [*Factoring with cyclotomic polynomials*](https://doi.org/10.1090/S0025-5718-1989-0947467-1), **Math. Comp. 52 (1989), 201–219**; preliminary FOCS 1985, 443–450 | Direct prior art for automorphisms, group-ring exponents, and \(\Phi_d(p)\) factoring methods. Primarily exploits known multiples/smoothness, rather than higher-order Lehman cancellation. |
| Hittmeir, [*A babystep-giantstep method for faster deterministic integer factorization*](https://arxiv.org/abs/1608.08766), **Math. Comp. 87 (2018), 2915–2935** | Sieves possible \(p+q\) residues; exponent \(1/4\) with a superpolynomial saving. |
| Hittmeir, [*A time-space tradeoff for Lehman’s deterministic integer factorization method*](https://arxiv.org/abs/2006.16729), **Math. Comp. 90 (2021), 1999–2010** | Proved exponent \(2/9\). |
| Harvey, [*An exponent one-fifth algorithm for deterministic integer factorisation*](https://arxiv.org/abs/2010.05450), **Math. Comp. 90 (2021), 2937–2950** | Proved exponent \(1/5\), through more extensive collision batching. |
| Harvey–Hittmeir, [*A log-log speedup for exponent one-fifth deterministic integer factorisation*](https://arxiv.org/abs/2105.11105), **Math. Comp. 91 (2022), 1367–1379** | The 2021 preprint and 2022 publication are the same work. Proves \(O(N^{1/5}\log^{16/5}N/(\log\log N)^{3/5})\). |
| Harvey–Hittmeir, [*A deterministic algorithm for finding \(r\)-power divisors*](https://arxiv.org/abs/2202.12401), **Res. Number Theory 8 (2022), article 94** | Lattice method for \(p^r\mid N\); not an extension-field Lehman method. |
| Hittmeir, [*Integer factorization as subset-sum problem*](https://arxiv.org/abs/2205.10074), 2022 preprint; **JNT 249 (2023), 93–118** | Further sieving/Fermat refinements. |
| Hales, [*Divisors in Residue Classes Revisited*](https://arxiv.org/abs/2410.05030), 2024 preprint | Extends Lenstra’s residue-class algorithm to certain number rings; does not establish the proposed collision speedup. |

**UNVERIFIED:** an identifiable result called “Lehman method in number fields,” or a published cyclotomic/higher-degree Harvey–Hittmeir construction with the proposed Taylor-cancellation advantage. I found no matching primary source; this is not a novelty proof.

**4(a). Close-prime BSGS — correction proved.** Put \(\Delta=q-p\), \(s=p+q\). Then
\[
s-2\sqrt N=\frac{\Delta^2}{s+2\sqrt N}.
\]
Thus \(a^{N+1-2\lfloor\sqrt N\rfloor}=a^{s-2\lfloor\sqrt N\rfloor}\pmod N\) has a discrete logarithm in an interval of width \(O(1+\Delta^2/\sqrt N)\). With a suitable base and candidate verification, BSGS costs
\[
O(1+\Delta/N^{1/4})
\]
group operations, plus powering/setup. An unknown gap bound can be handled by doubling.

It is **not instant throughout \(\Delta\ll N^{3/8}\)**: that condition permits costs approaching \(N^{1/8}\). Constant search work occurs at the \(N^{1/4}\) gap scale, where Fermat already takes \(O(1)\) square tests. Coppersmith also handles sufficiently accurate factor approximations. ([Coppersmith, **J. Cryptology 10 (1997), 233–260**](https://doi.org/10.1007/s001459900030))

Prior art includes Pollard’s interval-log factoring proposal, **Math. Comp. 32 (1978), 918–924, §3**; McKee–Pinch, *Old and new deterministic factoring algorithms*, **ANTS II, LNCS 1122 (1996), 217–224**; and Hittmeir 2018 above. McKee’s separate *Speeding Fermat’s factoring method*, **Math. Comp. 68 (1999), 1729–1737**, should not be conflated with a new BSGS result. ([Pollard](https://doi.org/10.1090/S0025-5718-1978-0491431-9), [McKee–Pinch](https://pure.royalholloway.ac.uk/en/publications/old-and-new-deterministic-factoring-algorithms), [McKee](https://doi.org/10.1090/S0025-5718-99-01133-3))

Your modulo correction is generally right, **with an exception**. If \(p=r+Mt\), then
\[
p+q\equiv r+Nr^{-1}+Mt(1-Nr^{-2})\pmod{M^2}.
\]
Usually only the residue modulo \(M\) is fixed. If \(r^2\equiv N\pmod M\), the last term vanishes and the sum *is* fixed modulo \(M^2\). More generally, choosing \(Br^2\equiv AN\pmod M\) fixes \(Aq+Bp\) modulo \(M^2\): this is precisely the weighted phenomenon in Harvey–Hittmeir, Proposition 3.3. ([Primary paper](https://arxiv.org/html/2105.11105))

**4(b). ECM curve plus twist — correction proved.** For a unit RHS ratio \(r\),
\[
(r/N)=(r/p)(r/q)=-1
\]
means the two local symbols differ. The points cover opposite twist classes at **exactly one** hidden prime, which need not be \(p\).

Counterexample: \(N=7\cdot11\), \(E:y^2=x^3+x+1\), and \(x=0,9\). The RHS values are \(1,739\); their ratio has Jacobi symbol \(-1\), yet both are squares modulo \(7\).

The orders \(p+1-t\) and \(p+1+t\), and x-only arithmetic that works on either twist, are established machinery; see Montgomery, *Speeding the Pollard and elliptic curve methods of factorization*, **Math. Comp. 48 (1987), 243–264, §10.3.1**. ([Original paper](https://wstein.org/edu/124/misc/montgomery.pdf))

Trying both gives two correlated opportunities for smooth point orders, at additional arithmetic cost. The union bound permits at most a factor-two probability gain when the marginals agree; it supplies no smoothness guarantee or change of asymptotic exponent. **UNVERIFIED:** priority for a particular paired-twist ECM variant, or a special speedup beyond those extra trials.