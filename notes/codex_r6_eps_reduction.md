> Source: GPT-6.1-sol round 6 (codex/r6_eps_reduction.md → .out), 2026-10-18. Curated copy; not yet independently verified.

**UNVERIFIED / open:** I found no classical randomized polynomial-time reduction establishing the requested implication, even for an exact \(\varepsilon\)-oracle. The closest literature is **Seurin’s \(2\)-\(\Phi/4\)-Hiding problem**, which concerns precisely this predicate. Adleman–McCurley O11 concerns general quadratic residuosity.

I read the requested project sections; no files were changed. Below, **PROVED** includes elementary derivations given here; **UNVERIFIED** means no supporting theorem was located. Assume \(\delta\ge1/\mathrm{poly}(n)\); otherwise a runtime polynomial in \(1/\delta\) need not be polynomial in \(n\).

**1. Literature map**

| Source | Exact statement and scope |
|---|---|
| **Adleman–McCurley (1994), C11/O11** | C11 takes \((a,m)\) and decides whether \(\gcd(a,m)=1\) and \(x^2\equiv a\pmod m\) is solvable. **O11a:** is C11 in \(P\)? **O11b:** does complete integer factoring C5 reduce to C11 by randomized polynomial-time reduction? They record C11 \(\le_P\) C5 and no progress in 1994. These are unrestricted problems, with variable \(a\), rather than fixed-\(-1\) average-case SOS. |
| **Its neighbours** | C9 counts distinct prime factors; they discuss estimating this through QR frequencies. C10 extracts \(e\)-th roots subject to \(\gcd(e,\varphi(m))=1\); factoring reducibility is O10, while removing that restriction permits Rabin’s reduction. C12 finds a prime-modulus nonresidue. C13 finds the least prime with a prescribed Legendre-symbol signature. C14 extracts roots modulo a prime; Rem14\(_{94}\) records Schoof’s deterministic polynomial algorithm for fixed \(a\). [Author’s original text](https://www.mccurley.org/papers/open.ps.gz) |
| **Seurin (PKC 2014), §3.1** | **CONJECTURED:** indistinguishability of uniformly sampled balanced Blum integers \(p,q\equiv3\pmod4\) and pseudo-Blum integers \(p,q\equiv1\pmod4\), called **\(2\)-\(\Phi/4\)-Hiding**. This is deciding \(\varepsilon\). For an equal mixture, prediction advantage \(\delta\) corresponds to distinguishing advantage \(2\delta\). His prime windows and mixture weights must be matched before transferring a statement to the project’s \(D_n\). [Paper](https://eprint.iacr.org/2013/256), [author’s slides, slide 13](https://www.iacr.org/workshops/pkc2014/slides/7-3-lossiness-Rabin.pdf) |
| **Jager–Schwenk, 2012 full version, Theorem 3** | **PROVED:** for each \(m=pq\), distinct odd primes, a generic-ring algorithm deciding whether a hidden uniform \(x\in J_m=\{x:(x/m)=1\}\) belongs to \(QR_m\), with advantage \(\alpha\) using \(t\) ring operations, yields a factor with probability at least \(\alpha/[2(t^2+4t+3)]\). The reduction runs it once, adds at most \(2t\) ring operations and \((t+2)^2\) gcds, and samples from \(J_m\) and \(\mathbb Z_m^*\). **This is a restricted-model result with a variable hidden element.** [Full version, §5.2](https://eprint.iacr.org/2009/621.pdf) |
| **Schoof (1985), §4** | **PROVED:** for fixed integer \(a\), a square root modulo a known prime \(p\), when one exists, is computable deterministically in \(O((\log p)^9)\), with an \(a\)-dependent constant. This does not cover unfactored composite moduli. [Original paper](https://pages.cs.wisc.edu/~cs812-1/Schoof85.pdf) |
| **Rybalov (2020)** | **PROVED, conditional:** assuming worst-case hardness of SOS decision and \(P=BPP\), there is no polynomial strongly generic algorithm—one succeeding outside an exponentially vanishing exceptional fraction. This concerns generic-case complexity over integers, not weak prediction on balanced semiprimes. [Paper and abstract](https://www.mathnet.ru/eng/pdm707) |

**PROVED:** SOS decision and the three fixed-element QR decisions lie below factoring. For odd prime \(\ell\), their local conditions are:

| Element | Local square condition |
|---|---|
| \(-1\) | \(\ell\equiv1\pmod4\) |
| \(2\) | \(\ell\equiv1,7\pmod8\) |
| \(-2\) | \(\ell\equiv1,3\pmod8\) |

CRT requires the condition at **every** prime divisor. Jacobi symbol \(-1\) settles a negative instance; Jacobi \(+1\) generally leaves the hidden local signs unresolved. **UNVERIFIED:** no standard-model factoring equivalence or worst-case-to-\(D_n\)-average-case theorem for these fixed elements was located.

General QR self-reductions randomize the **element while retaining the modulus**; they do not randomize a hard modulus into \(D_n\). This distinction is explicit in [Feigenbaum–Kannan–Nisan, §3](https://www.cs.yale.edu/homes/jf/FKN-Structures90.pdf).

**Rabin versus the fixed root.** **PROVED:** for fixed odd semiprime \(N\), if an algorithm extracts a valid root on an \(\alpha\) fraction of uniform quadratic residues, choose uniform \(x\in\mathbb Z_N^*\), submit \(x^2\), receive \(y\), and compute \(\gcd(x-y,N)\). Conditional on valid extraction, \(x\) is uniform among four roots, giving factoring probability \(\alpha/2\). Factoring supplies root extraction in the converse direction. [Rabin’s 1979 report](https://publications.csail.mit.edu/lcs/pubs/pdf/MIT-LCS-TR-212.pdf)

That proof does **not** establish factoring from one root of the fixed element \(-1\):

\[
\boxed{\text{one SOS representation of squarefree }N
\;\longleftrightarrow\;
\text{one root of }-1\bmod N}
\]

Both conversions are deterministic polynomial time:

- \(N=a^2+b^2\) gives \(r=ab^{-1}\pmod N\); squarefreeness ensures \(b\) is invertible.
- A root \(r\) gives a Gaussian gcd \(\gcd_{\mathbb Z[i]}(N,r+i)\) of norm \(N\).

Two roots \(r,t\) with \(t\not\equiv\pm r\) split \(N\). **Independent uniform** roots achieve this with probability \(1/2\); independent calls to an arbitrary root algorithm need not. It may always return the same root.

**2. Reduction map**

Write \(S(m)=[m\text{ is SOS}]\), \(Q(m)=[-1\text{ has a root mod }m]\), and \(\Omega(m)\) for prime factors counted with multiplicity.

**PROVED promise obstruction:**

\[
\Omega(kN)=2+\Omega(k)>2\qquad(k>1).
\]

Thus \(kN,Nr,Nr^2,N^2\) are outside **every** balanced-semiprime support. Their query distributions have total-variation distance **1** from every \(D_j\), regardless of congruence or bit length. Allowing arbitrary queries gives no accuracy guarantee on them.

| Attempt | Exact failure |
|---|---|
| **(a) Small-prime or known-factor multipliers \(kN\)** | All nontrivial queries violate the promise. Even granting a perfect global SOS oracle, after checking \(\gcd(k,N)\), coprimality gives \(S(kN)=S(k)S(N)\): the response merely repeats the target bit or masks it. Also, \(k\equiv1\pmod4\) does **not** imply \(S(k)=1\): \(k=21\) is a counterexample. Multipliers whose primes are all \(1\bmod4\) simply repeat \(S(N)\). |
| **(b) \(Nr^2\)** | Outside the promise, with \(\Omega(Nr^2)=4\) for prime \(r\). **Exactly:** \(S(Nr^2)=S(N)\), since square multiplication preserves every valuation parity. A perfect \(Q\)-oracle adds only the known constraint at \(r\). |
| **(c) Multiplicativity** | The computable Jacobi symbol satisfies \((-1/N)=\varepsilon^2=1\), hence \((-1/kN)=(-1/k)\): it loses the desired sign. SOS’s Boolean indicator combines by **AND** on coprime inputs, not signed multiplication. For example, \(21\) and \(209\) are coprime non-SOS semiprimes, and their product remains non-SOS, although multiplying their \(-1\) labels gives \(+1\). |
| **(d) Jager–Schwenk-style rerandomization** | Their proof exploits dependence on a hidden variable \(x\) and CRT variation among its inputs. Here \(A\) receives only \(N\). Replacing \(-1\) by \(-u^2\bmod N\) would require an oracle accepting a variable element; \(A(N)\) receives no such element. With \(x=-1\), the polynomial \(X+1\) vanishes globally and exposes no CRT disagreement. Public-input ring handles can be tracked by ordinary integer representatives. |
| **(e) Search-to-decision through varied moduli** | Standard witness descent requires a predicate such as “does a root exist with this prefix?” Neither SOS decision nor \(Q(N)\) supplies it. Encoding those prefix constraints into appropriately distributed semiprimes is the missing reduction. A tempting query \(S(N-a^2)\) tests a remainder’s being **two** squares, not one: \(65-2^2=61=5^2+6^2\), but \(61\) is not a square. |
| **(f) Quadratic forms / genus theory** | Queries corresponding to discriminants \(-4kN\) usually reduce to querying \(kN\), outside the promise. Genus characters describe a quotient \(G/G^2\); they do not produce square roots of classes, class orders, or useful ambiguous forms. Computing the usual character matrices from the hidden discriminant factorization is circular. |

Additional related-modulus checks:

- **PROVED:** \(S(N^2)=1\) for every \(N\), whereas \(Q(N^2)=Q(N)\) for odd \(N\), by prime-power lifting. Thus these two global decision predicates must not be conflated.
- **PROVED:** shifts do not preserve the label even within the project’s prime windows:
  \[
  3\,120\,349=1601\cdot1949,\quad\varepsilon=+1,
  \]
  \[
  3\,120\,353=1723\cdot1811,\quad\varepsilon=-1.
  \]
  Both have 22 bits and both factors’ top two bits set. This refutes label preservation under \(N\mapsto N+4\); it does not exclude a subtler reduction.
- **UNVERIFIED:** no suitable reduction using SOS decisions on nearby moduli was located. Classical multiplier-discriminant factoring methods require richer information: [Schoof, §6](https://www.mat.uniroma2.it/~schoof/mcpaper.pdf) computes useful ambiguous forms through class-group/regulator machinery.

For the “Ireland–Rosen-style” route, **PROVED:** \(\Delta=-4N\) is fundamental with three prime-discriminant factors, so
\[
\dim_{\mathbb F_2}\mathrm{Cl}(\Delta)[2]=2
\]
for either value of \(\varepsilon\). Moreover,
\[
(1,0,N),\qquad (2,2,(N+1)/2)
\]
are publicly available ambiguous forms that reveal no proper factor of \(N\). Producing **an** ambiguous form therefore does not suffice; producing a useful one is the extraction problem. I verified no Ireland–Rosen polynomial-time reduction overcoming this distinction.

**Average-case accuracy creates a second obstruction.** **PROVED:** support membership alone does not preserve advantage. For a query law \(Q\),

\[
\Pr_Q[A(M)=\varepsilon(M)]
\ge \tfrac12+\delta-\operatorname{TV}(Q,D_n).
\]

An allowed deterministic oracle may be wrong throughout any subset \(E\) with \(D_n(E)\le1/2-\delta\), yet satisfy the original guarantee. Repeating the same query cannot repair those errors. Adaptive or correlated queries require a distributional argument beyond their being valid semiprimes.

**3. Separations and barriers**

**UNVERIFIED:** I found no published relativized or black-box separation for the stated concrete arithmetic problem on mixed \(D_n\).

There is a closely relevant **conditional** barrier: Seurin assumes that factoring remains hard with a \(2\)-\(\Phi/4\)-Hiding decision oracle when proving limitations on tight black-box reductions for deterministic Rabin–Williams signatures. This is an additional hardness assumption, rather than a proof that the oracle cannot help factoring. [Paper, introduction and gap-one-way definition](https://eprint.iacr.org/2013/256.pdf), [author’s slides 23–24](https://www.iacr.org/workshops/pkc2014/slides/7-3-lossiness-Rabin.pdf)

**PROVED observation:** an unconditional oracle world where ordinary binary-integer factoring is hard would already establish unrelativized factoring hardness: an algorithm can ignore the oracle. Abstract hidden-ring separations consequently do not automatically settle this arithmetic question.

**4. A realistic 1–3-day target**

A rigorous restricted obstruction is plausible; resolving the unrestricted reduction is not a realistic expectation.

**ONE lemma to formalize and audit — elimination of multiplicative queries.**

> Let \(R\) be a polynomial-time oracle algorithm whose queries on target \(N\) all have the form \(kN^e\), with positive integers \(k,e\). If \(R\) factors with probability \(\eta(n)\) for **every** oracle satisfying the stated promise, then an oracle-free algorithm factors \(D_n\) with probability at least \(\eta(n)/2\).

**PROVED sketch:** choose the allowed oracle
\[
A_*(M)=
\begin{cases}
\varepsilon(M),&M\text{ is a valid balanced semiprime},\\
+1,&\text{otherwise}.
\end{cases}
\]
Among the permitted queries, only \(k=e=1\) can be valid. Guess \(\varepsilon(N)\) once, consistently simulate answers at \(N\), and answer \(+1\) elsewhere. The correct guess reproduces the entire adaptive transcript. This is restricted oracle elimination, not a computational separation.

**Experiment:** exhaustively enumerate small project-window semiprimes and adaptive query trees using \(k\le100,\ e\le3\). Compare \(A_*\)’s transcript with both guessed-bit simulations; require exact agreement on the correct branch and zero nontrivial multiplier queries inside the promise. I already checked 10,413 small cases of the underlying SOS and fixed-\(-1\) product/power identities without a failure.