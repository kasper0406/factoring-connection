# Adversarial referee report on `detGRH_results.md`, §§1–6

Date: 20 October 2026. The draft and the requested context notes were read. Additional checks of the randomness inventory used `lv6_results.md` and `recurse_results.md`. Repository files were not modified; `codex/detGRH_1.out` was not opened.

The central algebraic reduction survives: full-rank **kernel** relations on a generating factor base give a multiple of the group exponent, and Miller's method then factors under ERH. The principal problems are the interpretation of certificate-selected roots, an unsupported transfer from unrestricted relation production to DS-NFS, an incomplete deterministic cost audit, and assumptions smuggled into the analytic observations.

Exact computations are reproducible with:

```sh
/Volumes/git/factoring-connection/.venv/bin/python referee_checks.py
```

The checks covered MR on 250 odd moduli through 501, FL on 180 collected lattices, and Miller's subgroup argument on 440 modulus/exponent pairs and 47,080 bases outside the proposed subgroup. These finite checks support the proofs; they do not verify the analytic hypotheses. Results are in `referee_checks.json`.

## 1. Lemma MR(a), MR(b), uniqueness, and Corollary MR

**Verdict: CORRECT WITH FIXES.** Both parts of MR, including uniqueness, are correct. The corollary needs the precise small-prime qualification described below.

Write

\[
U=(\mathbb Z/N)^\times=\prod_i C_{n_i},\qquad
n_i=p_i^{e_i-1}(p_i-1).
\]

For (a), a section of squaring over the square subgroup exists exactly when each cyclic square subgroup has odd order, namely when every \(n_i/2\) is odd. This is equivalent to every \(p_i\equiv3\pmod4\). In that case squaring is an automorphism of \(Q\), so its inverse is the unique section. Repeated prime factors cause no problem.

For (b), define the map being split explicitly:

\[
\sigma:U/\{\pm1\}\longrightarrow Q,\qquad [u]\longmapsto u^2.
\]

The assertion must mean \(\sigma(y(s))=s\), not merely the weaker equality \([y(s)]^2=[s]\) in the quotient. With this definition, the proof using

\[
\chi(r)=y(r^2)[r]^{-1}\in U[2]/\{\pm1\}
\]

is valid. Its restriction to \(U[2]\) is the quotient map. A local sign is a nonzero class when \(\omega\ge2\), and cannot be the image under an exponent-two homomorphism of an even power. Thus all \(n_i/2\) must be odd. Uniqueness follows because \(Q\) then has odd order. For \(\omega=1\), \(\sigma\) is an isomorphism for every odd prime power, as claimed.

For the principal root,

\[
\frac{y(r^2)}r=\bigl((r/p_i)\bigr)_i
\]

up to the diagonal sign. Its quotient kernel is proper. A non-diagonal sign gives a proper gcd for either representative of the root class.

**Bach qualification.** The bound \(2(\log N)^2\) permits a prime dividing \(N\); the general quoted bound for a prime **coprime** to \(N\) outside a proper subgroup is \(3(\log N)^2\). Thus either:

- assume that trial division has removed every prime below \(2(\log N)^2\), as in Proposition EQ; or
- state “a small prime gives a factor by gcd, or is a unit with nontrivial \(\chi\)”; or
- use the unit bound with constant 3 and change the honesty cutoff accordingly.

The distinction is explicitly recorded in [Lamzouri–Li–Soundararajan, introduction §1.1](https://arxiv.org/pdf/1309.3595) and [Shoup, §10.5, printed p. 325](https://shoup.net/ntb/ntb-v2.pdf). Here and below, logarithms in these constants are natural logarithms. Dirichlet GRH suffices for this subgroup step.

Finally, “computing principal roots already factors \(N\)” is a reduction statement: it is unconditional with randomized inputs, or deterministic here with the repaired GRH small-base argument. It is not a claim that the multiplicative map cannot exist. It does exist for every product of odd prime powers whose primes are all 3 modulo 4.

## 2. Proposition FL(i)–(iii)

**Verdict: CORRECT**, in the stated free presentation \(\Gamma=\mathbb Z^n\), with \(2e_{-1}\in\Lambda_c\).

For (i), the isomorphism follows directly by halving in the torsion-free group \(\Gamma\):

\[
(\Lambda\cap2\Gamma)/2\Lambda\simeq G_B[2].
\]

Injectivity and surjectivity in the draft are correct. The sign generator contributes a nonzero class mapping to \(-1\). In particular, \(W'\) has dimension \(\dim G_B[2]\), not always \(\omega\).

For (ii), a class represented both by \(v\in\Lambda_c\) and by \(w\in2\Gamma\) differs by an element of \(2\Lambda\), hence \(v\) itself is even. This proves the image equality. Because the intersection already contains the sign line, it contains a fruitful class exactly when its dimension is at least 2.

For (iii), the dimension calculation is valid:

\[
\dim(V_c\cap W')\ge\dim W'-c.
\]

Full rank makes the stated finite quotient interpretation legitimate. Thus \(c\le\dim W'-2\) is a sufficient condition. For a semiprime with \(G_B[2]=U[2]\), this sufficient condition is exactly that \([\Lambda:\Lambda_c]\) be odd.

**Do not turn this sufficient condition into an equivalence.** Even index can be fruitful. For \(N=77\), \(FB=(-1,2,3)\), the group generated is all of \(U\); taking \(\Lambda_c=\Lambda\cap2\Gamma\) gives index 2 in \(\Lambda\), yet its \(\Phi\)-image is all four square roots of 1.

The generality sentence about arbitrary finitely generated certificate groups needs qualification. Halving is not unique in the presence of torsion. For example, \(\Gamma=C_2\) mapping isomorphically to \(\{\pm1\}\) has \(\Lambda=0\), whereas \(G_B[2]=C_2\). The displayed unquotiented isomorphism then fails. Use a finite free presentation, or formulate the appropriate quotient by the image of \(\Gamma[2]\) and account for the other torsion explicitly. The main proposition already uses a free presentation and is unaffected.

## 3. Proposition CD and the descent-to-residues criterion

**Verdict: CORRECT WITH FIXES.** The ambiguity calculation is correct. The “iff” is correct for the ambiguity relation, or for the stated canonical completion, but is false for an arbitrary selected root function.

For a fixed vector \(c\), two decompositions give

\[
2(g-g')=\lambda'-\lambda\in\Lambda_c\cap2\Gamma,
\qquad
\frac{\operatorname{ev}(g)}{\operatorname{ev}(g')}
=\Phi(\lambda'-\lambda).
\]

Every value in this image is attainable by varying the decomposition. Therefore the root class is intrinsically well defined for a fixed certificate when \(\Lambda_c\) is fruitless, and in that case it is a homomorphism on \(C\) modulo signs.

There are two distinct notions of “certificate” in the text:

1. a vector \(c\), whose decomposition has not yet been selected;
2. a vector together with a selected decomposition, or a fixed completion convention.

Keep them distinct. A convenient complete fix is to specify the canonical convention before CD. Choose relation vectors \(b_j\in\Lambda_c\) whose reductions form a basis of the image of \(\Lambda_c\) in \(\Gamma/2\Gamma\). For \(c\in C\), write its parity uniquely as \(\sum a_j\bar b_j\), with \(a_j\in\{0,1\}\), and set

\[
Y(c)=\operatorname{ev}\left(\frac{c-\sum a_jb_j}{2}\right).
\]

Then \(Y\) is actually a homomorphism on the integer group \(C\): addition of parity coefficients introduces only carries \(2b_j\), whose halved evaluations are 1. In particular,

\[
Y(2g)=\operatorname{ev}(g),\qquad
Y(c+2h)=Y(c)\operatorname{ev}(h).
\]

With this convention, CD(b)'s ratio formula and its full range \(G_B[2]\) are valid, and descent independent of certificate vectors holds iff \(G_B[2]\subseteq\{\pm1\}\). Alternatively, retain all decompositions and say that **all admissible root values** of a residue agree modulo signs iff this condition holds.

**Counterexample to the unrestricted selected-function reading.** Let \(N=77\), \(FB=(-1,2,3)\), and \(\Lambda_c=\Lambda\). Then \(G_B=U\), and define \(Y(c)\) to be the principal root of \(\operatorname{ev}(c)\) for every \(c\in C\). This root is admissible: choose a lift \(g\) of that root and put \(\lambda=c-2g\in\Lambda\). The selected function descends to residues even though \(|G_B[2]|=4\).

Concretely, the same vector \(c=2e_2\) permits

\[
g=e_2,\ \lambda=0,\ Y=2;
\qquad
g=2e_3,\ \lambda=2e_2-4e_3,\ Y=9.
\]

Both roots square to 4 modulo 77; their ratio is 43, a nontrivial square root of 1. Selecting the principal value consistently is different from allowing all admissible decompositions or using the normalized canonical convention. The draft must state which object it means.

## 4. Theorem M(S), M(G), vacuity, and Lemma SC

**Verdict: GAP as an application to Theorem A; the abstract dichotomy is CORRECT WITH FIXES.**

**Case (S).** The first three bullets are correct: the full-lattice image is diagonal, every collected lattice is fruitless, and the certificate root class is intrinsic and multiplicative on \(G_B^2\). For every \(r\in G_B\), its ratio to the certified root is diagonal.

The Miller conclusion is also correct for unit bases in \(G_B\) and a multiple of \(\exp G_B\), but its given explanation is false. **The whole Miller sequence need not lie in \(G_B[2]\).** For example,

\[
N=65,\quad FB=(-1,2),\quad G_B=\langle2\rangle\simeq C_{12},
\quad G_B[2]=\{1,64\}.
\]

With base 2 and exponent 12, the sequence is \(8,64,1\). Its first element has order 4. The correct explanation is that a nontrivial square root of 1 extracted at the last nonidentity step would belong to \(G_B[2]\); that group contains only the diagonal signs. Equivalently, the cyclic 2-primary part of \(G_B\) projects injectively into every local 2-primary factor, so the sequence cannot acquire different local terminal signs.

The B1 instance is correct for a **Blum semiprime \(N=pq\)** with distinct \(p,q\equiv3\pmod4\), in the scope of the cited even-weight, two-form degree-one sieve. On \(U[2]\),

\[
J_N(\epsilon_p,\epsilon_q)=\epsilon_p\epsilon_q.
\]

Since \(-1\in G_B\), inclusion in \(\ker J_N\) gives equality \(G_B[2]=\{\pm1\}\). This does not extend just from “all prime divisors are 3 modulo 4” to arbitrary multiplicities: \(J_{p^eq^f}(\epsilon)=\epsilon_p^e\epsilon_q^f\). Nor does B1 apply to arbitrary QS or hidden-root relations. Additional fixed coefficient factors in a two-form certificate must be included or separately shown to lie in the Jacobi kernel.

The final ERH bullet is correct after small-factor preprocessing and the Bach qualification in §5 below.

**Case (G), as a mathematical implication.** Fix an actual root selection \(F(R)=Y(c(R))\) and require

\[
F(r^2)=r\pmod{\{\pm1\}}
\]

for every small prime unit in the chosen cutoff. With the canonical convention of §3, \(c(r^2)=2e_r\) does imply this equality. Without that convention, certificate-vector honesty alone does not imply root-value honesty: the example \(c=2e_2\), \(Y=9\) modulo 77 shows the distinction.

With this correction and a valid small-unit bound, the contradiction with MR is sound. Nonmultiplicativity produces \(R,R'\) with

\[
\delta=F(RR')/(F(R)F(R'))\in U[2]\setminus\{\pm1\},
\]

and the displayed gcd is a proper factor. For non-Blum moduli, MR already rules out multiplicativity without honesty.

**This corrected implication is not vacuous.** For example, take \(N=1019\cdot1031=1{,}050{,}589\). Its least prime factor exceeds \(2\log^2N\approx384.47\), and \(-1,2\) already generate all of \(U\): the local orders of 2 are 1018 and 515. Define a root function to equal \(r\) on each small-prime square \(r^2\), and to equal the principal root elsewhere. These prescriptions are consistent, and this function is not multiplicative. Computation gives

\[
F(4)=2,\quad F(9)=3,\quad F(36)=1{,}049{,}564,
\quad\gcd(F(36)-6,N)=1031.
\]

Every such root choice has a certificate in \(2\Gamma\). One may take \(\Lambda_c=2\Lambda+\mathbb Z(2e_{-1})\), a fruitless collected lattice containing the sign control. Thus this is a genuine example of the abstract mechanism.

**The application to Theorem A is not established.** Its actual `DESCEND_ONCE` is a partial procedure. For a small prime \(r<B\), the split \((r^2,1)\), even when it is the shortest split, fails the requirement that the halves be products of distinct primes from \(J\), whose primes exceed \(B\). Theorem A does not return the alleged honest certificate there. It also does not certify every square. A direct smooth-square shortcut can be added, but it is an additional algorithmic rule. A successful extension to all squares, with the relevant cost, is not supplied by Theorem A. The inequality \(r^2<\sqrt N\) also needs an asymptotic threshold and separate treatment of small inputs.

**Lemma SC itself is correct:** ordinary integer factorizations satisfy

\[
e(RR')-e(R)-e(R')=0
\]

when the ordinary product is below \(N\). For the homomorphic canonical certificate convention, this also makes the root defect trivial. For arbitrary root selections, the vector identity alone does not imply that conclusion.

The following deductions are too strong:

- The word-length lower bound requires specified generators of size at most \(2\log^2N\), specified lifts, and a consistent multiplicative certification before reduction. A larger factor base or arbitrary sections is not covered by SC. The numerical size bound must use the actual generators.
- Cocycles of an arbitrary section need not generate the entire kernel. For example, for \(\mathbb Z\to\mathbb Z/3\), choose section values \(0,4,2\). All its additive defects generate \(6\mathbb Z\), whereas the kernel is \(3\mathbb Z\). Generation needs additional hypotheses on the lifts, or the usual Schreier generators for a fixed generating set.
- A Schreier enumeration is an available construction, not a proof that every search needs an \(N\)-size set. No algorithmic equivalence between finding a fruitful pair and the original derandomization problem is proved.

Finally, replace both assertions “\(\Phi\not\equiv0\) implies fruitful” and “structural fruitlessness is \(\Phi_{\rm full}\equiv0\)” by assertions about the **quotient** \(\bar\Phi\) in \(U[2]/\{\pm1\}\). The unquotiented map is never trivial here: it sends \(2e_{-1}\) to \(-1\).

## 5. Proposition Z: Bach, Miller's subgroup, and the determinant

**Verdict: CORRECT WITH FIXES.** The algebraic and Miller arguments are correct under the standing preprocessing. The blanket formulation of Bach's bound must be repaired.

The precise usable statements under Dirichlet GRH are:

- A proper subgroup \(H<U_m\) omits a prime \(\ell<2(\log m)^2\), **allowing \(\ell\mid m\)**.
- If \(\ell\nmid m\) is required, Bach's quoted general constant is 3.
- Consequently, after removal of prime factors below \(2(\log m)^2\), the unit primes below that bound generate \(U_m\).

These are the distinctions in [Lamzouri–Li–Soundararajan, §1.1](https://arxiv.org/pdf/1309.3595). The draft's generator conclusion is therefore valid in its rough-input domain. Without preprocessing, test gcds first, or use the unit bound with constant 3. Specify whether endpoints are included consistently; taking ceilings in the algorithm is harmless. This is not a smallness condition on the index of \(H\).

Let \(D=[\Gamma:\Lambda_c]>0\). Then

\[
\Gamma/\Lambda_c\twoheadrightarrow\Gamma/\Lambda\simeq G_B,
\qquad
\exp G_B\mid |G_B|\mid D.
\]

Once \(G_B=U\), this gives \(\lambda(N)\mid D\). For a nonzero square minor use its **absolute value**; it is the index of a full-rank sublattice contained in \(\Lambda_c\), hence a multiple of \(D\). No knowledge of the factorization is used to compute such a minor or its determinant.

For Miller, write \(D=2^st\), \(t\) odd. The maximal \(k\) exists because \((-1)^t=-1\), and no \(k\ge s\) is possible. The preimage

\[
H=\{a\in U:a^{2^kt}\in\{\pm1\}\}
\]

is a subgroup. If \(b^{2^kt}=-1\), CRT applied to \(b\) on one prime-power component and 1 on the others proves that \(H\) is proper. This works equally well when \(N\) has repeated prime factors.

For \(a\notin H\), the sequence from \(a^{2^kt}\) ends at 1. Its last nonidentity element \(x\) satisfies \(x^2=1\). It cannot be \(-1\): an occurrence after level \(k\) contradicts maximality, and an occurrence at level \(k\) puts \(a\) in \(H\). Thus every such base produces a proper factor. The proof does not require computing \(H\), \(b\), or \(k\); the algorithm tests all squaring levels for each small base.

The stated modular-exponentiation cost is a conservative valid bound. Deterministic polynomial-time determinant/HNF computation in the input bit size suffices for Z. This claim alone does **not** justify an \(L^{3a}\) bit bound for all the larger exact-dependency computations in §5 of the draft.

This exponent-multiple factoring step is classical: [Miller 1976, Lemma 5 and Theorems 3–4](https://www.cs.cmu.edu/~glmiller/Publications/Papers/Mi76.pdf); see also [Shoup, §§10.4–10.5](https://shoup.net/ntb/ntb-v2.pdf).

## 6. Proposition EQ and preprocessing

**Verdict: CORRECT WITH FIXES.** The stated equivalence for unrestricted kernel-relation production is correct. There is no circularity in its proof.

For relation production \(\Rightarrow\) splitting, extract independent rows and apply Z. Binary-encoded exponents are essential: the running time is polynomial in their **bit lengths**, not their magnitudes. Verification by modular exponentiation also has polynomial cost in that representation.

For splitting \(\Rightarrow\) relation production, recursive splitting gives the prime factorization and hence

\[
\lambda(N)=\operatorname{lcm}_{p^e\parallel N}p^{e-1}(p-1).
\]

The vectors \(\lambda(N)e_q\), with \(2e_{-1}\) if desired, are full-rank kernel relations and have polynomial output size for the small factor base. This computes a multiple of the needed orders without factoring any \(p-1\).

The domain restrictions are handled as follows:

- remove powers of 2;
- handle 1 and primes by deterministic primality testing;
- detect perfect powers by integer-root computations; recurse on their bases, which need not themselves be prime;
- trial-divide at the current cofactor's cutoff before every oracle call;
- call the relation/splitting procedure only on the remaining admissible odd inputs.

Only \(O(\log N)\) calls are needed. The displayed time comparison should use a nondecreasing upper bound \(T^*(N)=\max_{m\le N}T(m)\), or assume monotonicity. In particular, all-input splitting and complete factoring are equivalent up to this overhead. An almost-all-input splitting guarantee would not automatically survive recursion on the cofactors.

**The important limitation:** EQ proves equivalence with producing **arbitrary** full-rank kernel relations. It does not prove equivalence with **DS-NFS**, where the relations must arise from the prescribed bounded fields, seeds, and descents. Factoring \(N\) lets one print \(\lambda(N)e_q\); the proof does not realize these vectors as those descents or their exact dependencies. Remove “by Proposition EQ, DS-NFS is equivalent” unless a separate reverse construction is given. The proved implication is

\[
\mathrm{ERH}+\mathrm{DS\text{-}NFS}\Longrightarrow\text{deterministic factoring},
\]

with a separately established exact-postprocessing cost. The same issue affects the opening verdict's assertion that this particular NFS hypothesis is equivalent to deterministic \(L[1/2]\) factoring.

## 7. Consequences Z1–Z5

**Verdict: GAP in the stronger conclusions; Z1–Z2 and the underlying observation in Z3 are correct.**

**Z1: CORRECT, as a sufficient replacement architecture.** Given the stipulated kernel relations, Z needs neither hidden roots nor an \(\mathbb F_2\) fruitfulness argument. This does not show that a relation source already proved for Theorem A supplies the required relations.

**Z2: CORRECT WITH FIXES.** The power-map subgroup argument is Miller's argument. Higher 2-adic information survives in integer relations, and doubling a full-rank collected lattice preserves its usefulness for Z while destroying its non-diagonal \(\mathbb F_2\) root image. If the standing sign-control assumption is retained, double the other relations and reinsert \(2e_{-1}\); otherwise the literal doubled lattice omits that stipulated control. “Bockstein tower” is an interpretation, not an additional theorem established here.

**Z3: the full-lattice fact is CORRECT; the source-level conclusion is a GAP.** Since \(\Gamma/\Lambda\) is finite, \(\Lambda\otimes\mathbb Q=\mathbb Q^n\), so no nonzero rational linear functional annihilates the **whole** kernel. This excludes a rank obstruction for the unrestricted lattice.

It does not exclude a structural rank obstruction for the **allowed source** of smooth NFS relations or their projected exact dependencies. For example, over \(N=77\), \(FB=(-1,2,3)\), the restricted source

\[
\langle (2,0,0),(0,30,30)\rangle\subseteq\Lambda
\]

uses every column, yet has rational rank 2 rather than 3. This example is not a claim about the actual NFS source; it disproves the inference from finite ambient index alone. Proving that the permitted source spans enough directions is precisely part of the missing theorem. “Only sampling” is defensible only if sampling includes the unresolved arithmetic of which relations are obtainable.

A prime absent from all rows detects one particular rank defect. Enlarging the factor base adds columns as well as possible relations and does not automatically fix the defect. Say “enlarge/change the relation search and prove that the enlarged source supplies the required directions.”

**Z4: CORRECT WITH FIXES.** A usual QS row has

\[
2e_{x_i}-e_{Q(x_i)}
\]

only after adjoining \(x_i\) as a generator. If every \(x_i\) is a new independent column, its coefficient 2 forces its row coefficient to vanish in an exact rational elimination of the \(x_i\)-columns. The ordinary QS rows therefore do not automatically give the pure factor-base kernel relations needed by Z.

But “QS relations are not kernel relations” is too categorical. Sometimes \(x_i\) itself factors over the base. For example,

\[
N=77,\quad x=9=3^2,\quad Q(x)=4=2^2
\]

gives the genuine kernel relation \(3^4 2^{-2}\equiv1\pmod{77}\) on \((-1,2,3)\). Reused roots or known multiplicative dependencies among roots also defeat the “one new independent column per relation” premise. The correct conclusion is that Z does not automatically apply to **the usual QS relation supply**, not that it can never apply to a QS-derived relation set.

The contrast with NFS is conditional too. NFS rows initially live in a rational/algebraic certificate group. Exact elimination can produce pure rational kernel relations, but their projected rank and extraction cost must be proved. Prop XK supplies a method for one torsion dependency in its stated setting; it does not prove all the full-rank supply assertions being made here.

**Z5: correctly labelled conditional/open and unverified.** It proposes replacing the fruitfulness conjecture by another substantial hypothesis, namely full rank of the projected exact-relation lattice, plus adequate extraction cost. It does not discharge the LV conjecture. Retain that qualification in the opening verdict too.

## 8. §5 randomness inventory and claimed derandomizations

**Verdict: GAP.** The inventory is a useful outline, but it omits an explicit random choice in the proposed replacement and confuses several deterministic arithmetic costs.

| Rows | Referee finding / required correction |
|---|---|
| 1–3 | These are the field-law randomness in \((g_0,g_1,w,c)\). The cited proof averages over them. No small deterministic replacement is proved. Replace **irreducible** by **unresolved in the present proof**. Neither GRH nor the moment argument proves a lower bound against derandomization. |
| 4a | Uniform hidden roots make the square input uniform and are used for relation supply. Replacing them by factor-base products changes the input distribution. Theorem S′/K's success guarantee does not transfer to that list. CS does not prove a seed-list lower bound. Again, **unresolved**, not irreducible. |
| 4b | With known exponent vectors, a successfully descended product seed can contribute to an eventual kernel relation after exact algebraic elimination. This is a conditional reduction requiring full projected rational rank. It is not a drop-in derandomization of the original random-square supply. |
| 5 | Special-q search is already deterministic given the field and target prime. Correct. |
| 6 | Recursion removes the *separate* randomized shape-test subroutine if a deterministic outer algorithm with the requisite uniform cost/supply guarantees has already been constructed. Strong induction then applies to strictly smaller auxiliary inputs, with small inputs handled separately. This is a conditional closure property, not an independent deterministic factoring algorithm. |
| 7 | Deterministic elimination over \(\mathbb F_2\) is standard. Its cubic bit cost is \(L^{3a+o(1)}\); fast field elimination costs \(L^{\omega_{\rm MM}a+o(1)}\). Once true square dependencies are available, testing a basis finds a fruitful one whenever the reduced \(\Phi\) is nonzero. Exact integer/rational elimination has a different bit-cost analysis. |
| 8 | GRH can provide deterministic separating quadratic characters in fixed degree, subject to the exclusions and construction details below. It needs GRH for the relevant number-field extensions, a stronger formulation than just Dirichlet GRH used in Z. The Z route avoids these columns conditional on its own relation-supply hypothesis. |
| 9 | The Hensel prime in the cited algorithm is already the least suitable prime, selected by deterministic sieving. Random finite-field factorization and Tonelli–Shanks choices at that prime still need deterministic replacements. |
| 10 | Exhaustive evaluation finds the linear roots of \(f\bmod r\), at total \(L^{2a+o(1)}\) cost for primes \(r\le B\) and fixed \(d\). Correct. |
| 11 | Z removes trivial-gcd retries once full-rank kernel relations are supplied. It does not remove retries needed to obtain those relations. |
| 12–13 | The listed arithmetic tests, parameter computation, and AKS are deterministic. Correct. |

**Missing random prime in XK.** `dlog_results.md` §4, EXACT_DEPENDENCE step 2, explicitly picks a random prime

\[
q\in[Q,2Q],\qquad Q=M(\log p)^3.
\]

Its expected cubic cost uses a probability bound for preservation of a rational rank profile. This randomness is absent from the inventory even though rows 4b and 7 cite XK as deterministic. LLL and subsequent exact verification do not remove that choice. A deterministic rational rank computation and a deterministically chosen suitable lifting prime will fix determinism in polynomial time. Merely trying all candidate primes does not retain the asserted cubic bound: there can be \(O(M\log\log p)\) bad candidates by the supplied minor bound.

The distinction is real in the literature: the often-cited cubic certified solver of [Mulders–Storjohann (2004)](https://cs.uwaterloo.ca/~astorjoh/jscdense.pdf) is explicitly randomized. Deterministic normal-form methods do exist, but “polynomial time” does not mean the claimed exponent 3. For instance, [Storjohann, *Algorithms for Matrix Canonical Forms*, integer complexity bound in the introduction and §6.2](https://cs.uwaterloo.ca/~astorjoh/diss2up.pdf), gives a safe bound of the form \(\widetilde O(M^{\omega_{\rm MM}+1}\log\|E\|)\) for the relevant dense exact tasks. With small original entries, that is \(L^{(\omega_{\rm MM}+1)a+o(1)}\), not automatically \(L^{3a+o(1)}\). The draft must specify an algorithm and account for the number of independent projected relations required.

**The slack arithmetic itself is correct, with a missing range restriction.** In RA,

\[
\frac{a^*}{c^*}=\frac1{2+\sqrt2}=1-\frac1{\sqrt2},
\qquad
\frac{3a^*}{c^*}=0.8786796564\ldots<1.
\]

Thus a genuinely deterministic cubic postprocessing term fits. A deterministic exact bound with exponent \(C<2+\sqrt2\) also fits; a bound \(C=\omega_{\rm MM}+1\), using a matrix multiplication exponent below \(1+\sqrt2\), is a possible repair once the complete required extraction is implemented and analysed.

For nonrecursive Theorem A, \(a=\theta/\sqrt2\), so

\[
3a<(1+\theta)/\sqrt2\iff\theta<1/2.
\]

But Theorem A is stated for **all \(\theta<1\)**. “Which holds” is false in that full range. It does hold in the regime beating \(L[1/2,1]\), where \(\theta<\sqrt2-1\). At \(\theta=3/4\), the proposed cubic exponent is about 1.591, exceeding the original exponent 1.237. At \(\theta=1/2\), equality still preserves the leading \(L\)-constant, though the term is not slack. With exact cost \(C=\omega_{\rm MM}+1\), the analogous strict condition is \(\theta<1/\omega_{\rm MM}\).

**Character-prime details.** A uniform bound on \(\log|\Delta_{K(\sqrt h)}|\), as supplied in lv6, allows effective Chebotarev to distinguish each nonzero residual square class. It does not justify undefined residue-symbol entries: the usual evaluation \((a-b\alpha/\mathfrak p)\) requires \(\mathfrak p\) to avoid that row element. All prime ideals of polylogarithmic norm can meet the supports of the smooth atoms. Excluding only ramified primes is insufficient.

A repair must either construct appropriate local unit-part characters, or enumerate sufficiently many primes outside the factor base and bad supports, with the exclusion set in the effective bound. The latter can change “polylogarithmically many columns” into \(B\,L^{o(1)}\) columns. At fixed degree, GRH prime-counting estimates can support such an enumeration with suitable slack, but the row count and deterministic elimination costs must then be redone. Also distinguish “some prime ideal” from the required **first-degree** prime ideals; in fixed degree one can discard higher-degree primes using the corresponding counting estimate. A bare invocation of least-prime Chebotarev does not do all this bookkeeping.

**Other finite-field coins.** The literal source uses randomized selection of a simple root for the character primes, randomized polynomial factorization at the Hensel prime, and square-root/nonresidue routines in extension fields. Linear-root exhaustive evaluation alone does not describe all these replacements. For fixed degree they are removable without creating a new supply problem: Berlekamp-style deterministic searches and exhaustive shifts can be budgeted at \(p\,\mathrm{poly}(d,\log p)\) at the Hensel prime \(p\asymp B\). For example, for a field generator \(\zeta\), a quadratic nonresidue can be found among \(\zeta+a\), \(a\in\mathbb F_p\), by the character-sum bound for its irreducible norm polynomial, when \(p\) exceeds a degree-dependent constant. Supply that argument or cite a deterministic fixed-degree algorithm. It does not justify exhaustive search up to the enormous original character-prime range.

Finally, DS-NFS must charge the cost of **executing** the descents, shape tests, exact dependencies, and processing all fields, not just the cost of listing field/seed pairs. Bounds on total rows, coefficient bit lengths, and the number of fieldwise eliminations are needed before the claimed \(\max(c_0,C)\) running-time constant follows. The final determinant cost must use the bit lengths of the extracted relations, which can be much larger than the original entries of \(E\). Prop XK explicitly ends with a qualification that its fuller factoring extraction cost was not established; the draft cannot cite that sketch as a finished deterministic full-rank extraction theorem.

## 9. Observations FS and CS

**Verdict: WRONG as stated.** FS is a valid warning about the current full-progression proof, not a proved family-size floor. CS's Hölder calculation is correct **conditional on an additional spectral bound**, but that bound does not follow from GRH for the target set in the statement.

**FS.** Put \(t=N^{q_0+o(1)}\), \(X_F=N^{\xi+o(1)}\), with the cited proof requiring \(q_0<\xi/2\). Using its existing full-progression input requires

\[
X_F/t\ge t^{1+\varepsilon}
\]

for some fixed margin in the chosen parameters. Reproducing that argument by enumerating an entire consecutive progression is polynomial-size in \(N\). This does explain why directly shrinking the current random coefficient box is not justified by the existing argument.

The stronger inferences do not follow:

- Bombieri–Vinogradov is not “meaningless” for every shorter progression. The quoted theorem has a particular level of distribution; failure to fall within it is not a lower bound on all possible moment proofs.
- A sparse explicit subfamily need not be a short consecutive progression. It can spread samples throughout the original range. A second moment over an arbitrary initial segment of an arbitrary enumeration is not automatically a short coefficient box.
- The implication \(J\ll t\Rightarrow tJ\ll X_F^{1/2+o(1)}\) is arithmetically false. For instance, take \(t=X_F^{0.4}\), \(J=t^{0.9}\). Then \(J=o(t)\), but \(tJ=X_F^{0.76}\). The correct condition for the claimed short-interval regime is \(J\ll X_F^{1/2}/t\). A subexponential \(J=N^{o(1)}\) does satisfy that condition with the fixed strict level margin.

Relabel FS as **“limitation of the present full-progression argument; no lower bound for general explicit subfamilies.”** Assertions that no GRH theorem can help, or that a pointwise theorem is the unique possible remedy, exceed this observation.

**CS: recomputation.** Work in the unit group, of size \(g=\varphi(N)\). Assume first that \(E\) is a set of **distinct residues**, \(m=|E|\), and \(S\) is a set of units, \(s=|S|\). Let

\[
U_S=\max_{\chi\ne1}|\widehat S(\chi)|.
\]

Character orthogonality gives, up to the harmless conjugation convention,

\[
Z=\frac1g\sum_\chi\widehat E(\chi)|\widehat S(\chi)|^2,
\qquad \mathrm{Main}=\frac{ms^2}{g}.
\]

For conjugate exponents \(1\le p\le2\), \(2\le q\le\infty\), Parseval and interpolation give

\[
\sum_\chi|\widehat E(\chi)|^p\le g m^{p/2},
\qquad
\sum_{\chi\ne1}|\widehat S(\chi)|^{2q}
\le U_S^{2q-2}gs.
\]

Consequently,

\[
\boxed{\frac{|\mathrm{Error}|}{\mathrm{Main}}
\le\frac{g}{\sqrt m}\,s^{-2+1/q}U_S^{2-2/q}.}
\]

Under the **extra assumption** \(U_S\le sN^{-1/4+o(1)}\), and with \(s=N^{1/2+o(1)}\), \(g=N^{1+o(1)}\), this is

\[
m^{-1/2}N^{1/2+o(1)}.
\]

The cancellation of \(q\) is correct. With only \(U_S\le s\), the best choice in this range is \(q=2\), giving

\[
m^{-1/2}N^{3/4+o(1)}.
\]

For \(p>2\), using \(\sum|\widehat E|^p\le gm^{p-1}\) gives the stated conditional bound \(m^{-1/p}N^{1/2+o(1)}\), which is worse. These exponent calculations are not the error in CS.

There are three substantive defects:

1. **The GRH input is unsupported.** GRH bounds a character sum over the full prime interval \(J\). It does not bound the sum over an arbitrary field-dependent subset of good primes, or over the descendable halves formed from that subset. A character-biased filter can remove all cancellation. For example,
   \[
   S=\{1\le a\le\sqrt N:(a,N)=1,\ \chi(a)=1\}
   \]
   has \(\widehat S(\chi)=s\), and in balanced examples \(s\asymp\sqrt N\). For the modulus in §4, the script finds 513 such units in \([1,1025]\), with character sum 513. An asymptotic family follows by taking \(N=pq\), \(p<q<2p\), and \(\chi=(\cdot/p)\); the interval contains a full period modulo \(p\). A new theorem about the actual descendable set is required.
2. **A family of seeds can collide.** If \(E\) is an indexed family of products, Parseval gives \(\sum|\widehat E|^2=g\mathcal C_E\), where \(\mathcal C_E=\sum_x m_x^2\), not an upper bound \(gm\). The draft writes a lower bound and subsequently uses it as an upper bound. The general displayed ratio has \(\sqrt{\mathcal C_E}/m\) in place of \(m^{-1/2}\). If all seeds coincide, there is no saving from the number of indices.
3. **A sufficient absolute-error certificate is not a necessary seed count.** Even with the extra spectral assumption, these inequalities merely fail to certify positivity at small \(m\). They do not prove that a smaller explicit list fails. Indeed \(E=\{1\}\) meets \(S/S\) whenever \(S\ne\varnothing\), by \(u=v\); such a hit need not produce any useful rank. Thus a hitting statement also falls short of DS-NFS's rank requirement.

Retain the Hölder calculation as a conditional observation with distinctness/energy hypotheses. Delete the asserted GRH consequence for the actual target set and replace “needs” by “this particular absolute-value bound would only certify.” The analogous assertion that every deterministic use of expander mixing must enumerate all walks is also not a proved lower bound.

## 10. Overall verdict, prior art, and deterministic subexponential factoring

**Verdict: GAP in the draft's stronger assessments.** Its statement that this work has not constructed a deterministic subexponential factoring algorithm is justified. Several claims advertised as proved are stronger than their arguments.

The revisions with the greatest effect are:

1. State the certificate convention and root-value honesty explicitly; do not present Theorem M(G) as a property already satisfied by Theorem A's partial descent.
2. Retain EQ for unrestricted kernel relations. Remove its claimed converse for the restricted DS-NFS hypothesis.
3. Separate deterministic finite-field elimination from deterministic exact bit complexity; include XK's random prime and the full projected-rank extraction cost.
4. Scope “no structural obstruction” to the unrestricted lattice, and replace “irreducible” by “unresolved.”
5. Replace FS and CS's lower-bound language by statements of the limitations of the specified proof inputs. CS's spectral assumption is currently missing.

The defensible positive conclusion is: **under ERH, full-rank kernel relations on a suitably generating factor base suffice for deterministic factoring; finding such a relation set within the proposed NFS restrictions and time budget remains unproved.** This identifies a sufficient supply problem. It does not prove that every proposed relation source has no further structural obstruction, or that field/seed randomness is necessary.

**Prior art.** The determinant-to-exponent observation is an immediate finite-group presentation argument. The subsequent factoring reduction and its ERH derandomization are classical. The following sources support the relevant antecedents; they should not all be cited as if they state the draft's exact DS-NFS claim.

| Source | What it supports |
|---|---|
| [G. L. Miller, “Riemann's hypothesis and tests for primality,” JCSS 13 (1976), 300–317](https://www.cs.cmu.edu/~glmiller/Publications/Papers/Mi76.pdf) | Lemma 5 and Theorem 3 give ERH factoring reductions from suitable exponent functions; Theorem 4 treats the period problem. This is the primary origin of the power-map/gcd part of Z. |
| [E. Bach, *Discrete Logarithms and Factoring*, UCB/CSD-84-186 (1984)](https://www2.eecs.berkeley.edu/Pubs/TechRpts/1984/CSD-84-186.pdf) | §3 reduces factoring to obtaining exponents from discrete logarithms, and explicitly notes that the reduction becomes deterministic under ERH by finding a small element outside the failure subgroup. It does not supply an efficient deterministic relation collector. |
| [E. Bach, “Explicit bounds for primality testing and related problems,” Math. Comp. 55 (1990), 355–380](https://doi.org/10.1090/S0025-5718-1990-1023756-8) | The explicit small-element constants used to make the Miller enumeration concrete. The unit/nonunit distinction is as in §5 of this report. |
| [K. S. McCurley, “A key distribution system equivalent to factoring,” J. Cryptology 1 (1988), 95–105](https://www.mccurley.org/papers/equiv.pdf) | §3.1 discusses extracting a nonzero order multiple from exponent collisions and using Miller-type ideas to factor, with hypotheses on the chosen base. Relevant precedent, not a theorem asserting full-rank deterministic relation supply. |
| [K. S. McCurley, “The discrete logarithm problem,” Proc. Symp. Appl. Math. 42 (1990), 49–74](https://www.mccurley.org/papers/dlog.pdf) | §5 surveys factor bases, multiplicative relations, and index calculus. |
| [J. L. Hafner and K. S. McCurley, “A rigorous subexponential algorithm for computation of class groups,” JAMS 2 (1989), 837–850](https://research.ibm.com/publications/a-rigorous-subexponential-algorithm-for-computation-of-class-groups) | Classical relation lattices and integer normal forms for finite abelian groups. This is a Las Vegas algorithm, not a deterministic subexponential factoring result. |
| [V. Shoup, *A Computational Introduction to Number Theory and Algebra*, §§10.4–10.5](https://shoup.net/ntb/ntb-v2.pdf) | Explicitly treats factoring from a nonzero multiple of \(\lambda(N)\), and credits the algorithm and its ERH deterministic version to Miller. It also records the two Bach constants. |
| [E. Bach and J. Shallit, *Algorithmic Number Theory*, Vol. 1, MIT Press (1996)](https://mitpress.mit.edu/9780262526296/algorithmic-number-theory-volume-1/) | Standard background on exponent/order methods, ERH, and factoring. I verified the publication, but did not verify an exact numbered statement there of the draft's EQ formulation; do not cite an invented theorem number. |
| [R. Crandall and C. Pomerance, *Prime Numbers: A Computational Perspective*, 2nd ed. (2005)](https://link.springer.com/book/10.1007/0-387-28979-8), and [Pomerance, “A tale of two sieves” (1996)](https://math.dartmouth.edu/~carlp/PDF/paper109.pdf) | Standard sources on congruence-of-squares factoring and sieve relation collection. They provide context, not a verified citation for the precise deterministic full-rank equivalence. |
| [K. E. Stange, “Factoring using multiplicative relations modulo \(n\): a subexponential algorithm inspired by the index calculus,” 2022/2023, especially §2 and Theorem 3.2](https://arxiv.org/pdf/2211.06821) | A particularly close prior formulation: multiplicative relations modulo a composite, rational linear algebra, and extraction of an order or order multiple. Its relation-distribution and subexponential conclusions use a heuristic hypothesis, and its reduction is randomized. It should be discussed alongside Miller and Bach, even though it is not the identical deterministic EQ statement. |

I did not locate the exact small-factor-base search-problem equivalence labelled EQ in those words. Its proof is elementary from the cited classical ingredients, so “folklore-level reformulation” is a fair attribution; an assertion of substantive novelty would not be supported. Stange is sufficiently close that the literature discussion should not omit it.

**What is known about deterministic subexponential factoring?** The search did not locate a proved classical deterministic subexponential algorithm for **complete prime factorization of every integer**, under GRH alone or unconditionally. The established worst-case bound found is exponential in input length, \(N^{1/5+o(1)}\): [Harvey (2020 preprint/2021 publication)](https://arxiv.org/abs/2010.05450), with the explicit logarithmic improvement of [Harvey–Hittmeir (2021 preprint/2022 publication)](https://arxiv.org/abs/2105.11105). A [November 2025 paper by Chris Umans and Siki Wang](https://arxiv.org/abs/2511.10851) still describes \(1/5\) as the current deterministic exponent and proposes a further-conjecture route to \(1/6\). These results are not subexponential in \(\log N\).

The historical GRH \(N^{1/5+\varepsilon}\) comparison is consistent with [Schoof, *Quadratic fields and factorization* (1982)](https://ir.cwi.nl/pub/13014/13014D.pdf) and the later literature. It is not evidence that no stronger GRH theorem could exist. The randomized GRH algorithm of [A. K. Lenstra (1988)](https://www.sciencedirect.com/science/article/pii/S1385725888800222) must not be misclassified as deterministic.

**“Almost all \(N\)” needs a different answer depending on the task.** For finding **one proper factor**, there is already an elementary unconditional deterministic polynomial-time algorithm on a density-one set of composite integers: trial-divide up to \(\log N\). The integers with no prime divisor below any growing cutoff have density tending to zero. For a proof without a uniform sieve estimate, fix a constant \(K\); their upper density is at most \(\prod_{p\le K}(1-1/p)\), and then let \(K\to\infty\). Primes themselves have density zero and can be handled by AKS. Thus an unqualified claim that no deterministic subexponential *splitting* algorithm is known for almost all integers is wrong.

For **complete factorization**, this observation does not suffice: the remaining large cofactor may still be difficult. I found no theorem giving deterministic subexponential complete factorization on a density-one set under GRH alone in the checked sources. State this as a literature-search finding, specify the averaging family, and do not infer it from a worst-case record. Results on [oracle factoring of almost all integers](https://www.sciencedirect.com/science/article/abs/pii/S0885064X23000109) assume additional information and do not settle this question. Factoring polynomials over finite fields is also a different problem.

The draft's negative accomplishment claim is therefore defensible when it refers to complete factoring. Its sweeping barrier, equivalence-to-DS-NFS, and complete-derandomization claims require the revisions above.

## Final table of verdicts

| Claim | Verdict | Decisive point |
|---|---|---|
| 1. MR and Corollary MR | **CORRECT WITH FIXES** | MR and uniqueness hold; define quotient squaring precisely and repair the small-unit cutoff. |
| 2. FL(i)–(iii) | **CORRECT** | Valid for the free presentation and included sign control; the codimension condition is sufficient. |
| 3. CD | **CORRECT WITH FIXES** | Fix canonical/transport-compatible completion, or speak of all admissible roots; the unrestricted selected-function “iff” is false. |
| 4. M, G, SC | **GAP** | Correct abstract nonvacuous implication after root-honesty fixes; actual Theorem A descent lacks the asserted totality and small-square honesty. Miller explanation and Schreier claims need correction. |
| 5. Z | **CORRECT WITH FIXES** | Determinant divisibility and Miller subgroup proof hold; Bach's constant 2 permits nonunits, repaired by preprocessing. |
| 6. EQ | **CORRECT WITH FIXES** | Unrestricted relation production is equivalent; preprocessing is sound. The claimed converse for DS-NFS is not proved. |
| 7. Z1–Z5 | **GAP** | Ambient full rank does not prove source full rank; usual QS limitation needs qualification; NFS extraction remains conditional. |
| 8. Randomness table | **GAP** | XK random prime omitted; exact deterministic cubic cost and character-prime details are not supplied; “irreducible” is unsupported. |
| 9. FS and CS | **WRONG** | FS overstates a method limitation; CS assumes a GRH spectral consequence that is false for arbitrary target sets. Conditional Hölder exponents are correct. |
| 10. Overall/literature | **GAP** | Classical reduction survives; stronger structural and DS-NFS claims do not. Almost-all splitting is already easy and differs from complete factoring. |
