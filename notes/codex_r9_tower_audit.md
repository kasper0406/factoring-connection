> Source: GPT-6.1-sol r9_tower_audit (codex/r9_tower_audit.md → .out), 2026-10-18. Curated copy. Scripts in /tmp/claude-r9 (ephemeral).

**T1, T2, T3, T4, and both 2-part isomorphisms survive the audit. P5’s conclusion is repairable, but its written proof contains a false Frobenius identification and omits essential conditioning. T4′ proves a public XOR identity, not its claimed absence of information. Substantial prior art exists.**

I kept the repository read-only. Independent scripts and results are under `/tmp/claude-r9/`.

| Claim | Verdict | Finding |
|---|---|---|
| T1: quartic character suffices when \(N\equiv1\pmod8\) | **CORRECT** | The expansion, product identity, and reduction from \(a\bmod8\) to \(a\bmod4\) are valid. |
| T1(iv): ambiguity criterion and explicit \(h\) formulas | **CORRECT** | The change under \(a\mapsto a+2\) is computed correctly. |
| T1(v): quartic character is insufficient when \(N\equiv5\pmod8\) | **CORRECT** | Explicit collision: \((p,q)=(7,11)\) and \((23,59)\). |
| T2, with its stated \(N\equiv1\pmod8,\ k\ge6\) hypotheses | **CORRECT** | The difference factors into terms divisible by \(2^{k-3}\) and \(8\). |
| R1: \(h=[2,2,p^*]\) on undetermined QR₂-true pairs | **CORRECT** | Admissibility and the conductor-16 evaluation are valid. |
| T3, including the Rédei linearity step | **CORRECT** | All required symbols are defined; minimal ramification gives the required class-group interpretation. |
| T4: fourth-power membership decision rule | **CORRECT** | Both quartic characters must vanish; their changes are exactly \((h,h)\). |
| T4′: \(z_4\oplus z_8=[N\equiv9\pmod{16}]\) | **CORRECT** | Its linearity argument is valid. |
| T4′: “therefore these depths carry nothing about \(h\)” | **CORRECT WITH GAP** | A public XOR does not establish independence of either remaining bit from \(h\). |
| P5: asymptotic independence for fixed \(a,b\) | **CORRECT WITH GAP** | The statement holds with a repaired Chebotarev argument below. The supplied proof does not establish it. |
| P5’s identification of \([a,b,p^*]\) with ordinary \(\operatorname{Frob}_p\) | **WRONG** | Counterexample: \([5,5,-11]=1\), while ordinary Frobenius at \(11\) in \(\mathbf Q(\zeta_5)\) is trivial. |
| \(\mathrm{Cl}(-12N)[2^\infty]\cong\mathrm{Cl}(-3N)[2^\infty]\) | **CORRECT** | The conductor-2 extension map has kernel of order \(3\). |
| \(\mathrm{Cl}(-28N)[2^\infty]\cong\mathrm{Cl}(-7N)[2^\infty]\) | **CORRECT** | The extension map has trivial kernel; here even the full class groups are isomorphic. |

**1. T1: audit of the modular derivation**

References below are to [tower_results.md](/home/knielsen/code/factoring-connection/notes/tower_results.md:43).

- **Lines 43–46:** Correct. The decomposition
  \[
  (\mathbf Z/32\mathbf Z)^\times=\{\pm1\}\times\langle5\rangle
  \]
  gives the stated even characters. Given \(\varepsilon\), \(\psi(p)\) determines \(p\bmod16\).

- **Lines 62–63:** Correct. After dividing the binomial expansion by \(4\), all omitted terms vanish modulo \(8\):
  \[
  u\equiv\varepsilon\bigl(a+4\binom a2\bigr)
   =\varepsilon a(2a-1)\pmod8.
  \]

- **Lines 64–65:** Correct. The potentially questionable simplification is valid:
  \[
  4x^2y^2-4xy=4xy(xy-1)\equiv0\pmod8.
  \]
  Substituting \(x+y\equiv c\pmod8\) therefore gives
  \[
  uv\equiv(5-2c)a(c-a)\pmod8.
  \]

- **Line 66:** Correct. Replacing \(a\) by \(a+4\) changes \(a(c-a)\) by
  \[
  4c-8a-16.
  \]
  This vanishes modulo \(8\) precisely when \(c\) is even. Thus the quartic character suffices when \(N\equiv1\pmod8\).

- **Lines 67–71:** Correct. Under \(a\mapsto a+2\), the change is
  \[
  2c-4a-4\pmod8.
  \]
  It is nonzero exactly when \(2a\equiv c\pmod4\), equivalently
  \(\psi(p)=\psi(q)\). The four printed formulas for \(h\) check out.

- **Line 72:** Correct. For odd \(c\), the change under \(a\mapsto a+4\) is \(4\pmod8\), so \(s\bmod128\) changes by \(64\).

A concrete demonstration of the last point is:

| \((p,q)\) | \(N\bmod128\) | \(\varepsilon\) | \(\psi(p)\) | \(s\bmod128\) |
|---|---:|---:|---:|---:|
| \((7,11)\) | 77 | −1 | −1 | 18 |
| \((23,59)\) | 77 | −1 | −1 | 82 |

Thus the same quartic-character key gives different sums. The order-8 character distinguishes them.

**2. T2: the general lemma**

[Lines 74–78](/home/knielsen/code/factoring-connection/notes/tower_results.md:74) are correct under their stated hypotheses. In \(\mathbf Z/2^k\mathbf Z\),
\[
f_N(x')-f_N(x)
=(x'-x)\frac{xx'-N}{xx'}.
\]
The denominator is a unit. Since \(k\ge6\), the assumed congruence implies \(x'\equiv x\pmod8\), whence
\[
xx'-N\equiv x^2-1\equiv0\pmod8.
\]
The product consequently vanishes modulo \(2^k\).

The \(N\equiv1\pmod8\) hypothesis is essential. The preceding \((7,11)\), \((23,59)\) example also refutes the corresponding unqualified statement for \(k=7\): the distinguished primes agree modulo \(16\), but the sums differ modulo \(128\).

**3. T3: admissibility, linearity, and the class-group dictionary**

Stevenhagen’s theorem requires **all pairwise Hilbert symbols to be \(+1\), including at infinity**, and
\[
\gcd(\Delta(a),\Delta(b),\Delta(c))=1.
\]
This is a *three-way* gcd condition, not pairwise coprimality. Definitions 7.6 and 7.8 specify minimal ramification and the Artin interpretation; Theorem 7.9 supplies symmetry and linearity. [Stevenhagen, §§7–8](https://arxiv.org/pdf/1806.06250v3#page=18)

My check of [the T3 proof](/home/knielsen/code/factoring-connection/notes/tower_results.md:229) is as follows.

1. **Choice of \(r\): correct.** Quadratic reciprocity makes exactly one of \((q/p),(p/q)\) positive.

2. **Both \(R_4\) computations: correct.** With rows ordered by the dyadic discriminant, \(-r,-r'\), both matrices are
   \[
   R_4=
   \begin{pmatrix}
   0&1&1\\
   0&0&0\\
   0&1&1
   \end{pmatrix}.
   \]
   Hence \(r_4=1\).

3. **Identification of the surviving ambiguous class: correct.** For \(-4N\), the relation is \(P_rP_{r'}=1\) in the class group. For \(-8N\), it is \(P_rP_{r'}=T_8\). In both cases,
   \[
   C[2]\cap2C=\langle T\rangle,
   \]
   with \(T\ne1\).

4. **Second-kind decompositions: correct.** The three required splitting tests hold for both \((-r,4r')\) and \((-r,8r')\).

5. **Conversion to \(R_8\): correct.** Writing
   \[
   A=[-r,r',2],\qquad A'=[-r,2r',2],
   \]
   gives
   \[
   r_8(-4N)=1-A,\qquad r_8(-8N)=1-A'.
   \]
   The complements cancel in their XOR.

6. **Admissibility of the linearity calculation: correct.** For the three triples
   \[
   (-r,r',2),\quad(-r,2,2),\quad(-r,2r',2),
   \]
   the local checks printed in the note are sufficient. At \(2\), use
   \[
   -r\equiv1,\qquad r'\equiv7\pmod8.
   \]
   At infinity, no pair consists of two negative arguments. All unlisted odd places are harmless. The three-way discriminant gcd is \(1\), even for the repeated argument \(2\).

7. **Minimal ramification: correct.** For the pairs defining \(A,A'\), the quadratic discriminants are coprime and \(\Delta(-r)\equiv1\pmod8\). Thus the minimally ramified extensions are unramified over the relevant imaginary quadratic fields. Their symbols at the norm-2 ideal agree with the class-group symbols being used.

Consequently,
\[
\begin{aligned}
A'+A
&=[-r,2,2]\\
&=[2,2,-r]\\
&=\psi_2(-r)\\
&=[r\equiv7\pmod{16}].
\end{aligned}
\]
**There is no forbidden use of linearity here.**

The negative argument in the final symbol is also legitimate: \(-r\equiv1\pmod8\), and the allowed quadratic twists evaluate trivially there. Therefore T3 and C3 follow.

For example, my independent computation for \(p=23,q=71\) gives
\[
\mathrm{Cl}(-6532)[2^\infty]\cong C_8\times C_2,\qquad
\mathrm{Cl}(-13064)[2^\infty]\cong C_4\times C_2.
\]
Here \(r=23\), the rank XOR is \(1\), and
\[
-(1633+1)+64\equiv94=23+71\pmod{128}.
\]

**4. T4 and T4′**

[Lines 271–281](/home/knielsen/code/factoring-connection/notes/tower_results.md:271) are correct. “Every entry is \(+1\)” refers to multiplicative residue symbols; the associated \(\mathbf F_2\)-matrix is zero.

The two second-kind genus characters are independent. Their quartic lifts therefore test \(T\in4C\) by two simultaneous vanishing conditions. The required triples are admissible, and linearity gives
\[
(A',B')=(A,B)+(h,h).
\]
The printed decision rule follows immediately.

The \((0,0)\) outcome really can leave ambiguity **even after conditioning on \(K_0\)**:

| \((p,q)\) | \(N\bmod128\) | \((z_4,z_8)\) | \(h\) |
|---|---:|---|---:|
| \((17,257)\) | 17 | \((0,0)\) | 0 |
| \((73,137)\) | 17 | \((0,0)\) | 1 |

Both satisfy the T4 hypotheses; I independently computed both memberships.

T4 concerns the **particular class \(T\)**. When \(r_4=2\), its fourth-power membership cannot generally be replaced by the scalar 8-rank.

For T4′, the identity at lines 287–289 is correct. The information-theoretic conclusion is not established: a public XOR only gives
\[
I(h;z_4,z_8\mid K_0)=I(h;z_4\mid K_0),
\]
not zero. Proving zero requires a separate distribution argument. The identity also says nothing by itself about deeper memberships beyond \(4C\).

**5. P5: explicit proof error and a repair**

The failure at [line 317](/home/knielsen/code/factoring-connection/notes/tower_results.md:313) is concrete.

Take \(a=b=5,\ p=11\). The symbol is defined, and
\[
[5,5,-11]=\left(\frac{-11}{5}\right)_4=1\in\mathbf F_2.
\]
But \(11\equiv1\pmod5\), so its ordinary Frobenius in
\(F_{5,5}=\mathbf Q(\zeta_5)\) is trivial. The difference is the infinite-prime factor required when \(p^*=-p\). [Stevenhagen, Definition 7.8](https://arxiv.org/pdf/1806.06250v3#page=18)

There are two further gaps:

- Disjointness from \(\mathbf Q(\zeta_{16})^+\) alone does not justify independence **after conditioning on a class modulo \(8\)**.
- The finite-family paragraph needs to control the compositum and the admissibility conditions explicitly.

Here is a repair, which also proves the finite-family version.

Let \(F_i=F_{a_i,b_i}\) for finitely many fixed symbols, and put
\[
E=\prod_i\mathbf Q(\sqrt{a_i},\sqrt{b_i}).
\]
Choose a finite multiquadratic field \(M\) containing \(\mathbf Q(\zeta_8)\) and the quadratic fields needed to encode all remaining admissibility conditions. Put \(L=M\prod_iF_i\).

Every admissible prime splits in \(E\). If it is inert in \(\mathbf Q(\sqrt2)\), then necessarily
\[
\sqrt2\notin E.
\]
Moreover, \(\operatorname{Gal}(L/E)\) has exponent \(2\): its restrictions to every \(F_i/E_i\) and to \(M/\mathbf Q\) have exponent \(2\).

Let \(C=\mathbf Q(\zeta_{16})^+\). If \(C\subset L\), the square of its order-4 character would vanish on \(\operatorname{Gal}(L/E)\). Thus the quadratic subfield \(\mathbf Q(\sqrt2)\) of \(C\) would lie in \(E\), a contradiction. Since \(M\) already contains \(\sqrt2\),
\[
C\cap L=\mathbf Q(\sqrt2).
\]

The nontrivial automorphism of \(LC/L\) consequently exchanges the two possible values \(\psi(p)=i,-i\), while preserving:

- the class modulo \(8\);
- admissibility;
- every symbol value.

For negative \(p^*\), the infinite-prime correction is fixed on the chosen class and does not affect this pairing. Chebotarev gives equal limiting densities for the paired conjugacy classes.

Thus **P5 is true for fixed symbols, conditional on a nonempty admissible prime population**, but its proof should be replaced. This argument does not establish C5 for class-group data whose defining fields vary with the second prime.

**6. The conductor-2 isomorphisms**

[Lines 123–126](/home/knielsen/code/factoring-connection/notes/tower_results.md:123) are correct in the stated setting \(p,q>7\).

For \(d=-3N\) or \(-7N\), these are fundamental discriminants, and the order of discriminant \(4d\) has conductor \(2\). Extension of ideals gives a surjection
\[
\mathrm{Cl}(4d)\longrightarrow\mathrm{Cl}(d)
\]
with kernel order \(2-(d/2)\); the unit index is \(1\). The underlying exact sequence is standard. [Conrad, Theorems 5.2–5.3](https://kconrad.math.uconn.edu/blurbs/gradnumthy/conductor.pdf#page=14)

Here \(d\equiv5\pmod8\) gives kernel order \(3\), while \(d\equiv1\pmod8\) gives kernel order \(1\). Restriction to the 2-primary subgroups is therefore an isomorphism.

The qualification matters: the **full** groups for \(-12N\) and \(-3N\) are not isomorphic. For \(N=11\cdot19\), my class numbers are \(12\) and \(4\), respectively; both 2-parts are \(C_2^2\).

**7. Independent numerical reproduction and code audit**

My [exact checker](/tmp/claude-r9/audit_exact.py) imports **nothing from the repository**. It enumerates primitive reduced forms, composes through coprime equivalent forms and CRT, and enumerates entire power images \(C^2,C^4,C^8,\ldots\). It uses Python integers and no random subgroup generation.

| Independent check | Population | Failures |
|---|---:|---:|
| T1 identities, \(p<q<2000,\ p\equiv q\pmod4\) | 22,666 pairs | 0 |
| T1(iv), \(N\equiv1\pmod8\) | 11,288 pairs | 0 |
| T2, \(k=6,\ldots,11\) | All applicable pairs above | 0 ambiguous keys |
| T3, \(11\le p<q<300\) | 105 pairs | 0 |
| T4 under all its hypotheses, \(q<300\) | 14 pairs; 4 decided | 0 |
| T4′, \(q<300\) | 38 pairs | 0 |
| Conductor isomorphisms | 792 comparisons | 0 |
| Single-prime statements, \(p<20{,}000\) | 556 primes \(1\bmod8\), 565 primes \(7\bmod8\) | 0 |

I also reproduced the **128 quartic-character keys**, **64 \(K_0\) keys with 32 ambiguous**, and **5,572 undetermined pairs**. For \(N\equiv5\pmod8\), all 128 quartic-character keys were ambiguous; none of the 256 order-8-character keys were.

Independent prime enumeration reproduced the bound-4000 populations: **37,030 total**, **18,390 undetermined**, **9,453 T3**, **2,064 T4**, and **4,217 T4′**. I did **not** independently recompute every class group at that bound.

The small computations matched **2,376 stored class-group records** for six discriminants. Results: [main checks](/tmp/claude-r9/exact_results.json), [supplementary checks](/tmp/claude-r9/supplement_results.json), [membership examples](/tmp/claude-r9/extra_results.json).

For `src/tower/`:

- The reduced-form enumeration, Sylow stopping criterion, torsion-to-rank conversion, and membership iteration are mathematically sound.
- The existing tests share the imported composition implementation. Their agreement alone was not independent validation; the new checker avoids that dependency.
- The permutation calculations are finite-sample diagnostics. Null agreement cannot prove asymptotic independence or C5.
- `symbols.py` checks consequences through class groups and elementary residue symbols; it does not independently implement general Rédei symbols.

**8. Prior art**

**Yes: 8-ranks of these discriminants expressed through residue symbols and Artin conditions are established literature.** Priority for the exact paired XOR deserves a narrower statement.

| Source | Relevance and limitation |
|---|---|
| **Rédei–Reichardt (1934); Rédei (1939)** | The quartic-extension framework and subsequent trilinear-symbol theory are classical. Stevenhagen explains their historical relationship to \(R_4\) and \(R_8\). This establishes prior art for the mechanism, rather than a located verbatim T3 formula. [Stevenhagen, §§1,4](https://arxiv.org/pdf/1806.06250v3#page=7) |
| **Stevenhagen, *Ray class groups and governing fields* (1989)** | Establishes governing fields for quadratic 8-rank problems with one prime varying. Fixing \(q\) puts these discriminants in that framework. [Paper](https://www.numdam.org/item/PMB_1989___1_A1_0/) |
| **Milović, arXiv:1512.08034; IJNT 14 (2018), 2165–2193** | Directly studies \(-4pq,-8pq,8pq\), with \(p,q\equiv1\pmod4\). Proposition 3 gives an explicit residue-symbol criterion for an unramified cyclic octic extension through a specified genus extension, including the dyadic condition. It directly overlaps T4’s setting; its hypotheses exclude T3’s \(7\bmod8\) primes. [Full paper, §2](https://discovery.ucl.ac.uk/id/eprint/10073324/1/dpq8rank.pdf) |
| **Wu–Yue (2007), *8-ranks of Class Groups of Some Imaginary Quadratic Number Fields*** | Directly treats \(D=-4pq\) under T4’s \(p,q\equiv1\pmod8,\ (p/q)=1\) hypotheses, giving an \(r_8=2\) criterion involving quartic symbols and quadratic forms. This is clear prior art for T4-type statements. I checked the publisher’s abstract, not the paywalled full proof. [Article](https://doi.org/10.1007/s10114-007-0965-1) |
| **Chan–Koymans–Milović–Pagano, arXiv:1908.01752** | Sections 2.1–2.2 develop Rédei-symbol linearity and Artin-pairing cancellation identities. These are close methodological precedents, although their principal distribution problem concerns positive discriminants. [Paper](https://arxiv.org/pdf/1908.01752) |
| **Koymans–Pagano, *Higher Rédei reciprocity and integral points on conics*** | Develops substantially broader reciprocity and Artin-pairing machinery. It is relevant background, not a located source for the exact T3 XOR. [Paper](https://arxiv.org/pdf/2005.14157) |
| **Yue (2011), *8-ranks of class groups of quadratic number fields and their densities*** | Gives explicit congruence, quartic-symbol, and form criteria for mixed \(3\bmod4/1\bmod4\) prime pairs. Relevant neighboring literature, but its stated prime hypotheses do not directly cover T3 or T4. [Article](https://doi.org/10.1007/s10114-011-8273-1) |

I did **not locate the exact T3 XOR stated verbatim** in the inspected sources. It is a valid short consequence of established Rédei reciprocity. The note’s blanket “New: T1–T4” should therefore become a precise claim about this particular corollary and its application to \(p+q\bmod128\), with priority left unestablished.