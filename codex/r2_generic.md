Context: a research project on classical polynomial-time integer factorization (brief summarized below). Round 1 (by you, in a separate session) produced the attached candidate cards. Our working synthesis is below. Your job: check it hard, then go beyond it.

SYNTHESIS (to be checked):
(S1) Every operational mechanism tried so far in Tracks A/B/D (binary sections tested by random-monic-quotient PIT, quotient sketches F_N mod (X^r-c), rank/pivot extraction in poly-dimensional quotient algebras, lifted/divided derivatives), and the product-tree methods in Track C, are instances of ONE model: a uniform poly-size straight-line program over Z/NZ (possibly in a poly-dimensional R_N-algebra given by structure constants), with random inputs, followed by gcd tests of computed elements with N. Call this the "generic ring + gcd" (GR+gcd) model. In that model a factor appears iff some computed element vanishes in some but not all CRT components.
(S2) Hence, in GR+gcd, the question reduces to: construct uniformly from N a poly-size SLP C (in poly-dim algebra) such that the distribution of "C vanishes mod p" differs noticeably from "C vanishes mod q". All known positive instances draw their skew from (i) algebraic-group orders (p-1, p+1, Phi_k(p), #E(F_p), Jacobians, tori) or (ii) polynomials with many integer roots (Lipton 1994/Shub–Smale tau/Cheng), (iii) integer-root products as in BGS (cost sqrt).
(S3) F_N itself: F_N mod p = F_q(X)^p. For x in F_p, with t=(1+x)/x, roots satisfy t^q - (t-1)^q = 1, where q acts as q mod (p-1). So the root density of F_N (and its sections/sketches) in each component is governed by Fermat-curve-type/Jacobi-sum counts at characters of order gcd(q, p-1)-ish, i.e. it's another order-coincidence source, not a new one.

TASKS
1. Verify or refute S1–S3 precisely. For S3 derive the exact root count of F_q over F_p (and over F_{p^k}) in terms of character sums, and the probability that a random scalar/random low-degree quotient exposes a factor. Correct anything wrong.
2. Prior art (use web search, primary sources): what is known about factoring in the generic ring model? Check Aggarwal–Maurer (Eurocrypt 2009, "Breaking RSA generically is equivalent to factoring"), Jager–Schwenk (Asiacrypt 2009 / JoC, Jacobi symbol in the generic ring model), Damgård–Koprowski (Eurocrypt 2002, generic lower bounds for root extraction), Boneh–Venkatesan (1998), Maurer–Wolf, Lipton 1994, and any result giving (a) an unconditional lower bound for factoring by generic ring algorithms, or (b) an equivalence "generic poly-time factoring iff short SLPs with skewed root densities exist", or (c) a barrier saying such a lower bound would imply tau-conjecture/VP≠VNP-type consequences. State exact theorems, with numbers.
3. Prove the cleanest version of the reduction "randomized poly-time GR+gcd factoring of balanced semiprimes ⇔ existence of a uniform poly-size SLP family with root-density skew ≥ 1/poly" (or explain why the equivalence fails, e.g. adaptivity, branching on gcd outcomes, equality tests). This is scaffolding, not novelty.
4. The main creative task: propose at least FOUR concrete mechanisms that step OUTSIDE GR+gcd, i.e. that use operations which are not ring operations on residues: integer size/order (comparisons, floor, rounding, carries, bit extraction of canonical representatives), lattice reduction, continued fractions, Coppersmith-type small-root decoding, collision search on integer intervals (Harvey/HH-style), or mixing representatives of residues across different moduli (e.g. reduce a residue mod N as an integer mod an auxiliary prime). For each: a full candidate card (template in the brief §4), the precise way it uses non-generic structure, the smallest missing lemma, and the smallest exact experiment. Include at least one Track E card that specifies an exact partial-information format, cites the exact decoder theorem that consumes it (Coppersmith Theorem number, parameters), and proposes a non-oracle source for it. Be explicit about where each one would collapse back to sqrt/N^{1/4}/N^{1/5} cost.
5. Adversarial pass on your own task-4 cards, then a ranking.

Output Markdown. Label every claim: proved / experimentally supported / conjectured / refuted for this model / UNVERIFIED.

=== ROUND 1 CARDS (yours, for reference) ===
## 1. Mechanisms as candidate cards

The strongest surviving conditional reduction is **A1: binary support search with ordinary arithmetic circuits**. Randomized formal identity testing over \(\mathbb Z/N\mathbb Z\) is already available; the missing operation is repeated section construction without circuit growth.

The most concrete experimental mechanism is **B1: compare multiplication ranks in weighted quotient rings**. It survives a useful repair cycle, but its degree-two version has exact counterexamples, and the experiments show no speedup.

I independently checked the seed identities on all 831 composites through \(1000\), and checked every scalar evaluation of \(F_{561}\). One clarification to the supplied Lucas statement: for \(F_{pq}\bmod p\), the digit-dominated support excludes the two removed endpoints \(j=0,q\). For example,
\[
\operatorname{supp}(F_{35}\bmod5)=\{5,10,25,30\}.
\]
No candidate below uses hidden-characteristic digits operationally.

Write \(\widehat F_N=(1+X)^N-1-X^N\in\mathbb Z[X]\) for the integer lift. Let \(\mathsf M(n)\) denote integer multiplication cost, and \(\mathsf P(r)\) polynomial multiplication cost in ring operations. All candidate computations reduce coefficients modulo \(N\).

### A1 — Random-quotient identity tests guide binary support extraction

**ID and one-sentence mechanism:**  
**A1 — conditional reduction.** Construct ordinary circuits for successive even/odd sections, and use randomized formal identity tests to select a nonzero child.

**Target problem / supported input family:**  
Every composite \(N\), including prime powers and integers with several prime factors. Even inputs can be split immediately.

**Representation and how it is built from \(N\):**  
Start with the repeated-squaring circuit for \(F_N\), of size \(O(n)\). At depth \(d\), the current polynomial is
\[
S_{M,r}(Y)=\sum_{j\ge0}[X^{r+Mj}]F_N\,Y^j,\qquad M=2^d.
\]
The representation must be an ordinary arithmetic circuit, with shared gates. Merely attaching a “section” operator does not provide an efficient evaluator.

**Exact operation that extracts information:**  
Apply the even/odd gate rules, construct both child circuits, and test them modulo independently sampled monic polynomials \(h(Y)\) of degree \(e=O(n)\). A nonzero remainder certifies a nonzero formal section.

Repeat the quotient test sufficiently often before declaring a section zero.

**Why that information exposes a factor:**  
After \(O(n)\) selections, the accumulated residue bits identify an interior support exponent \(k\). Then
\[
[X^k]F_N\ne0\quad\Longrightarrow\quad 1<\gcd(k,N)<N.
\]

A tiny diagnostic example is
\[
F_{15}=5X^3+3X^5+10X^6+10X^9+3X^{10}+5X^{12}\pmod{15}.
\]
Selecting odd twice gives
\[
S_{4,3}(Y)=5.
\]
The recovered residue \(3\) already gives \(\gcd(3,15)=3\).

**Known lemmas being used:**  
The coefficient-annihilator identity; the section gate rules; randomized formal PIT over finite rings. The latter does **not** require constructing a field extension or knowing a factor of \(N\).

**New lemma required:**  
A uniform algorithm must construct the ordinary circuits encountered along the selected path in total polynomial time, with maximum size \(L(n)=n^{O(1)}\).

This is the sole remaining primitive for this randomized reduction. It is substantial: repeated constant-factor overhead is insufficient.

**Correctness argument or conjecture:**  
Assuming that construction lemma, use fresh randomness at every adaptive test. Amplify each possible false-zero decision to error at most \(\delta/O(n)\). The complete path is correct with probability at least \(1-\delta\). Verify the final GCD; restarting after failure gives an expected polynomial-time splitter under the construction lemma.

**Cost in \(n\), including construction and extraction:**  
One circuit evaluation modulo a degree-\(e\) monic polynomial costs \(O(L(n)\mathsf P(e))\) ring operations. The elementary PIT proof below needs \(O(e\log(n/\delta))\) repetitions per decision. Thus a conservative search bound is
\[
O\!\left(nL(n)e\mathsf P(e)\log(n/\delta)\right)
\]
ring operations, plus circuit construction. Here \(e=O(n)\), coefficients have \(n\) bits, and straightforward storage is \(O(L(n)en)\) bits.

The existing dense quotient implementation instead costs \(O(n\mathsf P(2^d))\) at depth \(d\).

**Randomness and success probability:**  
Only the identity tests require randomness. Polynomial amplification suffices. The conditional algorithm has verifiable outputs.

**Closest prior method and the actual difference:**  
AKS uses quotient-ring identities for primality; the proposed output here is a support exponent. More decisively, Arvind–Mukhopadhyay–Srinivasan, Section 6, Algorithm 2 and Theorem 6.6, already give randomized circuit PIT over finite rings using random monic quotients, crediting Agrawal–Biswas. That primitive is known, rather than a new conjecture. [AKS](https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/primality_journal.pdf), [finite-ring PIT](https://arxiv.org/pdf/0801.0514).

Generic coefficient extraction is a different problem: \([X^T]\prod_i(1+X^{w_i})\) counts subset-sum solutions. Koiran–Perifel study the zero-coefficient problem, including characteristic-zero hardness under strong nondeterministic Turing reductions. Those results do not establish hardness for this \(F_N\), or for randomized formal PIT. [Koiran–Perifel](https://doi.org/10.1016/j.tcs.2007.08.008).

**Smallest informative experiment:**  
Hash-cons the section transformations on repeated-squaring circuits and follow one diagnostic support path and one PIT-selected path. Measure unique reachable gates after every depth, including gates needed to evaluate the chosen child.

**Counterexample or scaling behavior that would refute this version:**  
The obvious constant-coefficient linear-state compression is already refuted: a single scalar section observable can require exactly \(2^d\) states. Proof below. Exponential DAG growth would likewise refute the tested circuit construction, though not every possible construction.

**Possible repair if that happens:**  
Seek nonlinear circuit identities specific to binomial powers, or an efficiently evaluable representation with section operations stronger than ordinary circuits. The repair must specify its evaluator; a nested-section description alone does not suffice.

---

### B1 — Detect different hidden-field ranks in weighted quotient algebras

**ID and one-sentence mechanism:**  
**B1 — executable, experimentally weak.** Compute multiplication by a weighted quotient sketch of \(F_N\), and expose a nonunit when its two hidden-field ranks differ.

**Target problem / supported input family:**  
The rank-gap guarantee is stated for distinct-prime semiprimes. Every reported nontrivial GCD remains a valid certificate on arbitrary inputs.

**Representation and how it is built from \(N\):**  
Choose public \(r,c\), check parameter GCDs, and compute
\[
h=X^r-c,\qquad f=F_N\bmod h.
\]
Store \(r\) coefficients. Build the \(r\times r\) matrix \(B_f\) for multiplication by \(f\).

First inspect multiplication by
\[
W=X(X+1)(X^2+X+1).
\]
For \(\gcd(N,6)=1\), its roots are public roots of \(\widehat F_N\). Classify factors found from this small polynomial as baseline discoveries. Investigate \(F_N\) only when \(W\) is a unit in the quotient.

**Exact operation that extracts information:**  
Perform elimination using only unit pivots. Before inversion, compute the GCD with \(N\). Full pivot scanning either finds a useful nonunit or proceeds by valid ring operations.

Over a field,
\[
\operatorname{rank}B_f=r-\deg\gcd(h,f).
\]
Different CRT ranks therefore force a useful nonunit somewhere in such elimination.

**Why that information exposes a factor:**  
A pivot or inspected entry divisible by exactly one hidden prime gives a factor by GCD. A rank difference guarantees this phenomenon; equal ranks do not exclude it.

A repaired example is
\[
N=47959,\quad r=2,\quad c=46575,
\]
for which exact exponentiation gives
\[
f=9209+16167X.
\]
Both coefficients and the public-mask resultant are units modulo \(N\), but
\[
\det B_f
=9209^2-46575\cdot16167^2
\equiv37810\pmod{47959},
\]
and
\[
\gcd(37810,47959)=199.
\]
The diagnostic ranks are \(1\) modulo \(199\) and \(2\) modulo \(241\).

**Known lemmas being used:**  
Multiplication rank in a polynomial quotient; unit-pivot elimination; the perfect coefficient pairing proved below.

For analysis only,
\[
F_{pq}\bmod p=F_q(X)^p.
\]
The candidate computes neither \(p\) nor \(F_q\).

**New lemma required:**  
For an explicitly constructible distribution of polynomially bounded \(r,c\), prove an inverse-polynomial probability of unequal CRT ranks on the intended independent-prime distribution.

For uniform accepted weights, this has a precise finite-field formulation. If \(\pi_{\ell,j}\) is the conditional probability of rank \(j\) in the \(\ell\)-component, then CRT independence gives
\[
\Pr(\text{rank gap})=1-\sum_{j=0}^{r}\pi_{p,j}\pi_{q,j}.
\]
The missing result is a lower bound for this quantity, not merely a nonzero-sketch theorem.

**Correctness argument or conjecture:**  
Rank-gap extraction is proved. The required frequency of rank gaps is unproved. A uniform positive guarantee for the repaired degree-two rank model is **refuted**.

**Cost in \(n\), including construction and extraction:**  
Per probe,
\[
O(n\mathsf P(r)+r^3)
\]
ring operations and \(O(r^3)\) GCD checks suffice with straightforward elimination; storage is \(O(r^2n)\) bits. Thus polynomially bounded \(r\) gives polynomial cost per attempt. Repetition count is unresolved.

**Randomness and success probability:**  
Weights can be sampled uniformly from \(R_N\), rejecting public-mask failures. Exact degree-two success probabilities below range from \(0\) to positive constants. No general lower bound is established.

**Closest prior method and the actual difference:**  
AKS establishes compositeness through low-degree quotients; this mechanism asks for an imbalance in quotient-algebra ranks that forces factor extraction. It does not use sparse interpolation: the expanded support is exponentially large on balanced semiprimes, so finite-field sparse-interpolation guarantees do not supply the decoder. [AKS](https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/primality_journal.pdf), [van der Hoeven–Lecerf](https://www.texmacs.org/joris/ffsparse/ffsparse.html).

**Smallest informative experiment:**  
For each small semiprime and degree \(r\), enumerate weights separately in the two diagnostic fields, obtain the complete rank distributions, and calculate the exact CRT rank-gap probability. This is more discriminating than counting occasional successful GCDs.

**Counterexample or scaling behavior that would refute this version:**  
For \(N=101\cdot113\), every admissible quadratic weight gives rank \(2\) in both components. The same holds for \(107\cdot193\) and \(1009\cdot1481\). Hence varying \(c\) alone cannot rescue degree two.

Moreover, at \(N=35\), avoiding \(W\) eliminates rank gaps for **every monic quotient**, regardless of degree. This is a rank obstruction, not an assertion that all matrix entries are units.

**Possible repair if that happens:**  
Vary quotient degree or quotient shape. Analyze equal-rank nonunits separately, with a matched random-GCD baseline. Dividing out one copy of \(W\) is insufficient: hidden-characteristic multiplicities persist.

---

### C1 — Test collisions of consecutive roots through a compressed discriminant

**ID and one-sentence mechanism:**  
**C1 — exact reduction, extraction primitive unresolved.** Replace factorial evaluation by a squarefreeness/unit-discriminant test for a consecutively rooted polynomial given by a shift grammar.

**Target problem / supported input family:**  
The least-prime-factor prefix predicate for every composite \(N\).

**Representation and how it is built from \(N\):**  
For binary \(K\), describe
\[
P_K(T)=\prod_{j=0}^{K-1}(T-j)
\]
through
\[
P_{a+b}(T)=P_a(T)P_b(T-a).
\]
Binary recursion gives an \(O(\log K)\)-size **parameterized shift grammar**. It is not automatically an ordinary arithmetic SLP: evaluating the same subexpression at different arguments requires additional work.

**Exact operation that extracts information:**  
Determine whether
\[
D_K=\operatorname{disc}(P_K)
\]
is a unit modulo \(N\), using the grammar without materializing \(K\) coefficients or \(K\) roots.

A concrete recursive target is the cross-resultant in
\[
D_{2m}=D_m^2\operatorname{Res}\!\bigl(P_m(T),P_m(T-m)\bigr)^2.
\]

**Why that information exposes a factor:**  
Directly,
\[
D_K=\prod_{d=1}^{K-1}d^{\,2(K-d)}
    =\prod_{j=1}^{K-1}(j!)^2.
\]
Thus
\[
D_K\in R_N^\times\quad\Longleftrightarrow\quad K\le p,
\]
where \(p\) is the smallest prime divisor. Binary search recovers \(p\).

For \(N=35\),
\[
D_5\equiv29,\qquad D_6\equiv15,\qquad D_8\equiv0\pmod{35}.
\]
Their GCDs with \(35\) are \(1,5,35\), respectively.

**Known lemmas being used:**  
The discriminant product formula; resultants of products; the factorial-unit equivalence.

The exact recurrence
\[
D_{K+1}D_{K-1}=K^2D_K^2
\]
also holds. It does not permit cancellation after nonunits appear.

**New lemma required:**  
A faster unit-resultant algorithm for these translation grammars. A meaningful intermediate target is
\[
K^{1/2-\eta}\operatorname{poly}(n,\log K)
\]
for some fixed \(\eta>0\); polynomial time in \(n,\log K\) would solve the main objective.

**Correctness argument or conjecture:**  
The reduction and recurrences are proved. Faster extraction is conjectural. At present this is an exact specification of a different operation, without an algorithmic advantage.

**Cost in \(n\), including construction and extraction:**  
Grammar construction costs \(O(\log K)\). The unit predicate already has a roughly \(\sqrt K\) product/block baseline. Dense discriminant computation is unnecessary and can be more expensive.

If the proposed predicate cost were \(T(n,K)\), binary-search factor extraction costs \(O(n)\) queries. An improvement over the square-root product primitive would not automatically beat the best deterministic factorization bound.

**Randomness and success probability:**  
The stated grammar/resultant route targets deterministic computation. A randomized unit test would need amplification across the adaptive search.

**Closest prior method and the actual difference:**  
Pollard/Strassen batch consecutive products; BGS Lemma 13 already finds a prime divisor bounded by \(b\), or certifies none, in
\[
O\!\left(\mathsf M(\sqrt b\log N)+
\log b\,\mathsf M(\log N)\log\log N\right)
\]
bit operations. BGS Theorem 11 gives complete factorization in \(O(\mathsf M(N^{1/4}\log N))\). Here the proposed operation is a unit test on a discriminant grammar, rather than evaluation of block products. No faster implementation is established. [BGS](https://mathexp.eu/bostan/publications/BoGaSc07.pdf).

Harvey obtains \(O(N^{1/5}\log^{16/5}N)\); Harvey–Hittmeir divide this by \((\log\log N)^{3/5}\). These remain the relevant deterministic comparison, rather than trial division. [Harvey](https://arxiv.org/html/2010.05450v1), [Harvey–Hittmeir](https://arxiv.org/html/2105.11105v1).

Lipton’s 1994 preview concerns short programs for polynomials with many rational roots; it should not be quoted as an unrestricted uniform theorem about factorial descriptions. The \(\tau\)-conjecture concerns integer roots versus ordinary constant-free circuit size; it does not apply directly to a stronger shift grammar. [Lipton preview](https://page-one.springer.com/pdf/preview/10.1007/3-540-58691-1_45), [\(\tau\)-conjecture formulation](https://arxiv.org/pdf/math/0304100).

Tal’s Theorem 1.1 takes a **known prime** \(P\), a supplied divisor \(d\mid P-1\), and \(m<P\), obtaining \(\widetilde O(d^c+\sqrt{P/d})\) modular-factorial computation. This does not provide the unknown-factor-modulus unit operation here. [Tal, Theorems 1.1–1.3](https://eccc.weizmann.ac.il/report/2026/211/download).

**Smallest informative experiment:**  
Implement the cross-resultant recursion while retaining shift grammars, instrumenting every conversion into dense coefficients, evaluation points, or subresultant states. Compare only the Boolean unit answer against a transparent block-product implementation.

**Counterexample or scaling behavior that would refute this version:**  
A purported short recurrence that divides by \(D_K\) after \(K>p\) is invalid. A recursion that expands \(\Theta(K)\) shifted evaluations has lost its compression.

**Possible repair if that happens:**  
Use unit certificates or fraction-free identities that remain meaningful at collision points. Any repair must reduce actual work, rather than retain only a short symbolic description.

---

### D1 — Divide the lifted derivative before reducing modulo \(N\)

**ID and one-sentence mechanism:**  
**D1 — exact restricted-family mechanism; generic scalar version fails.** Retain information erased by ordinary differentiation through
\[
H_N(X)=\widehat F_N'(X)/N
      =(1+X)^{N-1}-X^{N-1},
\]
then seek nonunit scalar values.

**Target problem / supported input family:**  
Distinct odd-prime semiprimes. The success formula depends explicitly on \(\gcd(p-1,q-1)\).

**Representation and how it is built from \(N\):**  
A repeated-squaring circuit of size \(O(n)\). No large integer coefficients or actual lift to modulus \(N^2\) are needed to evaluate the resulting expression.

**Exact operation that extracts information:**  
Choose \(a\in R_N\), evaluate \(H_N(a)\), and compute its GCD with \(N\).

A repaired version removes the public factor:
\[
H_N=(2X+1)Q_N.
\]
The quotient \(Q_N\) also has a divide-free \(O(n)\)-gate circuit.

**Why that information exposes a factor:**  
Let
\[
d=\gcd(p-1,q-1).
\]
For \(a\ne0,-1\) modulo \(p\), substitute \(t=(a+1)/a\). Then
\[
H_N(a)=0
\iff t^{q-1}=1.
\]
The map excludes \(t=1\), so there are exactly \(d-1\) roots modulo either factor. Consequently uniform scalar sampling splits with probability
\[
\boxed{\frac{(d-1)\bigl(p+q-2(d-1)\bigr)}{pq}}.
\]

For \(N=91\), \(d=6\), giving probability \(50/91\).

**Known lemmas being used:**  
Fermat reduction of exponents; cyclicity of finite-field multiplicative groups; CRT independence.

**New lemma required:**  
For a broader algorithm, an efficiently constructible exponent or sampling operation must expose roots beyond public torsion with inverse-polynomial probability. Uniform sampling and affine changes of variables do not supply this.

The initially tempting lemma—“the divided derivative gives a useful scalar signal on all semiprimes”—is false.

**Correctness argument or conjecture:**  
The exact probability is proved. The public root \(a=-1/2\) is simple in both components. Removing it leaves exactly \(d-2\) roots in each.

If \(d=2\), \(Q_N(a)\) is a unit for every \(a\in R_N\). Moreover,
\[
\gcd(H_N(a),N)=\gcd(2a+1,N).
\]
On that family the original scalar procedure is precisely a public linear-polynomial GCD test.

**Cost in \(n\), including construction and extraction:**  
\(O(n)\) modular multiplications and one GCD per attempt; polynomial storage. Expected time is this cost divided by the displayed success probability.

**Randomness and success probability:**  
Uniform \(a\) gives the exact formula. For balanced factors and bounded \(d\), success is \(O(N^{-1/2})\). Favorable order-correlated families can have constant success.

**Closest prior method and the actual difference:**  
The surviving root condition is a multiplicative-order condition, close to Pollard \(p-1\) arithmetic. Cheng explicitly discusses the distinction between modular exponentiation and modular factorial computation. The proposed difference is the compact exponentiation observable obtained by dividing the integer derivative; its analysis shows that it does not escape the order bottleneck. [Cheng](https://qcheng2023.github.io/paper/factorial.pdf).

**Smallest informative experiment:**  
Exhaustively compare GCDs of \(H_N(a)\), \(Q_N(a)\), and \(2a+1\) on semiprimes classified diagnostically by \(d\).

**Counterexample or scaling behavior that would refute this version:**  
Every \(d=2\) input refutes the stripped scalar version: it never returns a factor, regardless of the scalar sampler.

**Possible repair if that happens:**  
Change the algebra or observable, rather than repeat scalar evaluation. Higher-dimensional evaluation returns to the rank/sketch questions of B1; a bounded lift alone does not change divisibility of the original value.

---

### E1 — Count short-interval intersections of CRT discriminant masks

**ID and one-sentence mechanism:**  
**E1 — concrete filtering representation, two unresolved extraction requirements.** Encode possible factor sums through small-prime quadratic-residue masks, then search their intersection with a short interval without enumerating CRT classes.

**Target problem / supported input family:**  
Balanced semiprimes; for a direct Coppersmith transfer, initially restrict to \(q-p\ge\delta p\) for fixed \(\delta>0\).

**Representation and how it is built from \(N\):**  
For public odd auxiliary primes \(\ell\nmid N\), store
\[
S_\ell=\{s\bmod\ell:s^2-4N\text{ is a square, including }0\}.
\]
Each mask costs polynomial time in \(\ell,n\). The combined predicate is a short conjunction of local masks, without an expanded list modulo \(L=\prod\ell\).

**Exact operation that extracts information:**  
Compute the exact number of admissible integers in an interval,
\[
C([A,B])=
\#\{s\in[A,B]:s\bmod\ell\in S_\ell\text{ for every }\ell\},
\]
then recursively descend into occupied intervals.

The concrete computational object is a CRT interval-counting problem, equivalently a structured multiple-choice sum problem.

**Why that information exposes a factor:**  
The true \(s=p+q\) survives every mask because
\[
s^2-4N=(q-p)^2.
\]
An exact surviving sum is verified by an integer square test. Alternatively, a polynomial-size list of sufficiently narrow factor intervals permits small-root postprocessing.

The local mask size is exactly
\[
\boxed{|S_\ell|=\frac{\ell+\left(\frac N\ell\right)}2}.
\]

**Known lemmas being used:**  
The discriminant identity; quadratic character sums; CRT.

Across a complete period, the exact density is
\[
\rho=\prod_{\ell\mid L}\frac{|S_\ell|}{\ell}.
\]
Averaged over translations, a length-\(W\) interval has \(W\rho\) survivors.

**New lemma required:**  
Two independent obligations remain:

1. Output-sensitive or polynomial-time exact interval counting for these masks.
2. A useful bound on false survivors in the particular interval containing the planted factor sum.

There is presently no honest reduction of both obligations to a smaller established primitive.

**Correctness argument or conjecture:**  
Mask soundness and density are proved. Efficient decoding and adequate false-candidate elimination are conjectural.

For converting sum intervals into factor intervals,
\[
p(s)=\frac{s-\sqrt{s^2-4N}}2,\qquad
\left|\frac{dp}{ds}\right|=\frac p{q-p}.
\]
Thus a short sum interval does not uniformly give an equally short factor interval near \(p=q\).

**Cost in \(n\), including construction and extraction:**  
For polynomially many polynomial-size auxiliary primes, mask construction and storage are polynomial in \(n\). Explicit CRT expansion costs
\[
\prod_\ell |S_\ell|,
\]
which is prohibitive. Direct interval scanning costs \(\Theta(\sqrt N)\) for a general balanced interval. The proposed compressed counting cost is unknown.

**Randomness and success probability:**  
The basic masks are deterministic. Randomizing auxiliary primes does not, by itself, prove either required bound.

**Closest prior method and the actual difference:**  
Hittmeir’s hyperbolic sieve and reduction to multiple-choice subset sum are the immediate prior construction. The only proposed difference is an output-sensitive exact interval-count operation on the compressed masks; the local sieve itself is not new. [Hittmeir, Section 5](https://arxiv.org/pdf/2205.10074).

Coppersmith’s Theorem 4 recovers factors from the high-order \(\frac14\log_2N\) bits of a factor. For balanced semiprimes, intervals of width \(N^{1/4-\varepsilon}\) provide a safe interface after boundary guesses. A compressed conjunction of masks is not yet that information. [Coppersmith, Theorem 4](https://nymity.ch/anomalous-tor-keys/bibliography/pdf/Coppersmith1997a.pdf).

**Smallest informative experiment:**  
Measure both survivor counts and the size of an exact interval-counting representation as masks are added. Record them separately.

For \(N=11413\), the tested interval contains \(108\) integers. Eight masks leave one survivor, but the expanded CRT list contains \(332\,640\) classes. Fourteen masks still leave one survivor while that list contains \(13\,881\,466\,368\,000\) classes.

**Counterexample or scaling behavior that would refute this version:**  
Near-half local densities do not imply independence on one short interval. Few survivors do not imply cheap counting. Either exponentially many false leaves or exponentially growing counting states defeats this version.

**Possible repair if that happens:**  
Aim first for a restricted-family interval theorem, or a structured approximation that outputs actual factor intervals with certified coverage. Do not treat a small local-mask description as already decoded information.

---

### F1 — Correlate quadratic norm-one group orders using \(N-\operatorname{Jacobi}(D,N)\)

**ID and one-sentence mechanism:**  
**F1 — independent connection, proved restricted-family splitter.** Exponentiate random elements of a quadratic norm-one group by an exponent constructed directly from \(N\), seeking identity in exactly one CRT component.

**Target problem / supported input family:**  
Distinct odd-prime semiprimes whose shifted prime orders have a large common divisor and unequal sizes. This route does not use \(F_N\), coefficient sections, or factorials.

**Representation and how it is built from \(N\):**  
Choose \(D\) with \(\gcd(D,N)=1\), and use the rank-two algebra
\[
A_N=R_N[T]/(T^2-D).
\]
Sample \(\alpha=u+vT\), checking that its norm \(u^2-Dv^2\) is a unit. Set
\[
z=\alpha/\bar\alpha,\qquad
J=\left(\frac DN\right),\qquad m=N-J.
\]

**Exact operation that extracts information:**  
Compute
\[
z^m=a+bT
\]
by binary exponentiation and return a nontrivial value of
\[
\gcd(N,a-1,b).
\]
Failed norm inversions are separately valid GCD opportunities.

**Why that information exposes a factor:**  
Let
\[
\epsilon_p=\left(\frac Dp\right),\quad
\epsilon_q=\left(\frac Dq\right),\quad
A=p-\epsilon_p,\quad B=q-\epsilon_q,\quad d=\gcd(A,B).
\]
The local norm-one groups are cyclic of orders \(A,B\), and
\[
\gcd(m,A)=\gcd(m,B)=d.
\]
Conditional on \(\alpha\) being a unit, \(z\) is uniform in their CRT product. Therefore
\[
\boxed{\Pr(\text{identity in exactly one component})
=\frac dA+\frac dB-\frac{2d^2}{AB}}.
\]

**Known lemmas being used:**  
Quadratic finite-field/split-algebra norm-one groups; the elementary Hilbert-90 map \(\alpha\mapsto\alpha/\bar\alpha\); cyclic-group root counts.

**New lemma required:**  
A new, efficiently constructible order-selective multiplier or algebra family must create substantially more favorable local orders. Choosing \(D\) alone accesses only four order pairs.

**Correctness argument or conjecture:**  
The stated probability is proved. If \(A\mid B\) and \(B\ge2A\), it is at least \(1/2\). If \(A=B\), it is zero.

Examples:
\[
N=35,\ D=2:\quad A=B=6,\quad P=0;
\]
\[
N=77,\ D=2:\quad A=6,\ B=12,\quad P=\frac12.
\]

**Cost in \(n\), including construction and extraction:**  
Polynomial time per sample: a Jacobi symbol, norm GCD/inversion, and \(O(n)\) operations in a rank-two algebra. Expected splitting time is this cost divided by the displayed probability.

**Randomness and success probability:**  
Uniform \((u,v)\), conditioned on unit norm, gives the stated uniform distribution. Rejection has bounded expected cost; norm failures may themselves split \(N\).

**Closest prior method and the actual difference:**  
This uses the same quadratic/Lucas arithmetic underlying Williams’ \(p+1\) method. The proposed variant uses the explicit exponent \(N-J\), rather than a smoothness-bound exponent, and its success condition is the common-divisor formula above. It is not a new general group-factorization method. [Williams’ method as implemented in GMP-ECM](https://members.loria.fr/PZimmermann/ecmnet/cov/ecm/pp1.c.gcov.html).

**Smallest informative experiment:**  
Compute the four diagnostic order-pair GCDs and compare the exact success formula with exhaustive small norm-one groups.

**Counterexample or scaling behavior that would refute this version:**  
Twin primes with \(\epsilon_p=-1,\epsilon_q=1\) give \(A=B\), hence zero identity-separation probability. For \(p=10007,q=10037\), the four possible GCDs are only \(2,2,4,6\); every twist has small success probability.

**Possible repair if that happens:**  
Smooth multipliers recover familiar \(p\pm1\) methods. An actual advance requires a new way to construct selective multipliers or local orders.

## 2. Adversarial self-review

A separate mathematical critic and complexity auditor checked the draft identities and costs. The following is the subsequent adversarial synthesis.

| Card | Assumption most likely to fail | Mathematical reason |
|---|---|---|
| **A1** | Repeated sections retain small, efficiently constructible ordinary circuits. | For \(M=2^d,c=2\), one scalar section sequence has minimal constant-coefficient recurrence order exactly \(M\). Linear-state compression cannot supply the missing lemma. |
| **B1** | Small quotient degrees yield rank gaps frequently. | After excluding public roots, complete quadratic weight distributions can be concentrated at rank \(2\) in both hidden fields. Nonzero sketches and occasional nonunit entries do not imply rank separation. |
| **C1** | Unit-only output makes compressed discriminants substantially easier. | Its predicate is exactly the factorial-unit predicate. The short discriminant recurrence requires invalid cancellation precisely after the informative threshold. |
| **D1** | Dividing the lifted derivative reveals new scalar information generically. | When \(d=2\), \(H_N=(2X+1)Q_N\) with every \(Q_N(a)\) a unit. Its GCD is exactly that of the public linear factor. |
| **E1** | Near-half local masks give a cheaply searchable sparse global set. | Complete-period density supplies neither short-interval sparsity nor an efficient interval counter. These are separate requirements. |
| **F1** | Varying the quadratic parameter creates many useful local orders. | It changes only the signs in \(p\pm1,q\pm1\). Equal orders give zero separation; bounded common divisors give small probabilities. |

The most consequential repair was to **B1**.

Initially,
\[
F_{35}\bmod(X^2-2)=33+4X
\]
looked useful: both coefficients are units, whereas its multiplication determinant has GCD \(7\) with \(35\).

But
\[
\operatorname{Res}(X^2-c,X^2+X+1)=c^2+c+1.
\]
Thus \(c=2\) already supplies \(7\) through a tiny public polynomial. The analogous initial successes at \(N=77,c=2\) and \(N=143,c=3\) are explained by the same baseline.

The repaired experiment excludes these public-root singularities. It then finds a genuine additional rank phenomenon at \(N=47959\), but also exact failures:

| \(N\) | Local quadratic rank counts \((\#0,\#1,\#2)\), after filtering | Exact rank-gap probability |
|---:|---|---:|
| \(101\cdot113\) | \((0,0,99)\), \((0,0,111)\) | \(0\) |
| \(107\cdot193\) | \((0,0,105)\), \((0,0,189)\) | \(0\) |
| \(199\cdot241\) | \((0,6,189)\), \((0,0,237)\) | \(2/65\) |
| \(251\cdot401\) | \((1,8,240)\), \((0,6,393)\) | \(1661/33117\) |
| \(1009\cdot1481\) | \((0,0,1005)\), \((0,0,1479)\) | \(0\) |

Therefore **the repaired quadratic rank-only model has no positive uniform success guarantee**, even if every weight is tried. Higher degrees and equal-rank nonunits remain distinct questions.

There is also a useful distinction between rank and elimination behavior. An invertible matrix can contain nonunits—for example,
\[
\begin{pmatrix}1&p\\0&1\end{pmatrix}\pmod{pq}.
\]
Consequently, an equal-rank factor discovery is not evidence for the rank-gap mechanism.

## 3. Ranking and experiment specifications

The scores below are research priorities: clarity \(5\) means a fully specified operation; missing-lemma score \(5\) means a small remaining gap. A zero denotes a proposed general extension already contradicted in its stated model.

| Rank | Card | Clarity | Missing-lemma score | Product | Assessment |
|---:|---|---:|---:|---:|---|
| 1 | **A1** | 5 | 2 | 10 | One precise remaining construction lemma; formal randomized zero testing is known. |
| 2 | **B1** | 5 | 2 | 10 | Cheap exact tests; varying degree is essential after quadratic counterexamples. |
| 3 | **C1** | 4 | 2 | 8 | Clear unit-resultant interface, but no demonstrated computational advantage. |
| 4 | **E1** | 4 | 1 | 4 | Two substantial obligations remain. |
| 5 | **F1** | 5 | 0 | 0 | Valid restricted-family splitter; twist-only generalization fails. |
| 6 | **D1** | 5 | 0 | 0 | Exact calibration mechanism; stripped scalar model can have no signal. |

The first three merit the next exact experiments. None currently has a small technical gap sufficient for a polynomial-time factorizer.

**A1 — the most discriminating experiment: pathwise circuit growth.**

- **Inputs:** \(143,323,899,11413,47959\), plus \(49,125,561\) and prime \(101\).
- **Operation:** Freeze a hash-consed ordinary-circuit implementation of the section gate rules, including explicit constant simplification. Follow a diagnostic path to a known support exponent and an operational path selected by random-quotient PIT.
- **Measure:** Reachable unique gates, construction operations, evaluation work, and peak storage at every depth. Count transformed dependencies needed by the chosen child.
- **Refutation:** Growth proportional to \(2^d\), or another exponential in depth, refutes this implementation as the required construction lemma—even if its printed description stays short. It does not refute all possible nonlinear compression.

**B1 — the most discriminating experiment: exact rank distributions before scaling.**

- **Inputs:** Exhaustive balanced semiprimes with distinct primes in \([101,251]\).
- **Public configuration:** \(r\in\{2,3,4,5,7,8\}\); uniform weights, with public parameter and \(W\)-resultant checks.
- **Operation:** The diagnostic harness enumerates weights locally to calculate every \(\pi_{p,j},\pi_{q,j}\). The candidate itself receives only \(N,r,c\).
- **Measure:** Exact rank-gap probability, equal-rank nonunit frequency, and improvement over matched random-GCD work.
- **Refutation:** An input whose two rank distributions agree as the same point mass for every degree in the menu refutes that menu’s rank-only splitting guarantee. Exact probabilities below a predeclared bound refute that quantitative bound.

**C1 — the most discriminating experiment: does unit-resultant recursion avoid state expansion?**

- **Inputs:** \(N=101\cdot113\) and \(107\cdot193\); binary cutoffs \(K\le512\), including cutoffs around the diagnostic least factor.
- **Operation:** Implement the cross-resultant recursion on shift grammars, returning only unit/nonunit or a certified factor.
- **Measure:** Distinct shifted states, dense coefficients ever materialized, resultant/subresultant dimensions, ring operations, and GCDs. Compare against a block-product unit test.
- **Refutation:** Requiring \(\Omega(\sqrt K)\) state or arithmetic work defeats the proposed sub-square-root implementation. A recurrence that cancels nonunits fails algebraically, independent of timing.

I ran smaller exact diagnostics for all six cards. The quotient scaling test used 12 independently sampled semiprimes at each of 8, 12, and 16 **factor bits**, with 48 public probes per instance:

| Factor bits | Total probes | Probes passing public filters | Factors from subsequent matrix inspection | Actual CRT rank gaps |
|---:|---:|---:|---:|---:|
| 8 | 576 | 512 | 123 | 2 |
| 12 | 576 | 574 | 10 | 0 |
| 16 | 576 | 576 | 2 | 0 |

Matched independent random-GCD trials, using the observed GCD budgets, have expected successful-trial counts approximately \(302.9,54.6,3.50\). Thus these measurements provide **no evidence of an advantage**. Most matrix discoveries were equal-rank events.

Additional checks included:

- All seed identities through \(1000\), including prime powers and multifactor composites.
- All 561 scalar evaluations of \(F_{561}\).
- 1,791 operational probes on inputs \(2\) through \(200\), checking every reported factor.
- Exact derivative root counts, discriminants, CRT masks, and small norm-one groups.
- Random-quotient PIT: \(F_{561}\) was detected in all 128 trials despite being zero at every scalar. The formally zero \(F_{17}\) produced zero in all trials. Nilpotent-content examples modulo \(9\) and \(25\) were also detected.

These are correctness and mechanism diagnostics, not asymptotic evidence. Counters instrument explicit quotient arithmetic, GCD calls, and inversions; they do not expand the internal bit work of Python’s built-in arithmetic.

Reproduce with:

```bash
python3 -B phase1_experiments.py > phase1_summary.json
```

Saved artifacts: [experiment code](/home/knielsen/code/factoring-connection/codex/phase1_experiments.py), [raw results](/home/knielsen/code/factoring-connection/codex/phase1_results.json), [mathematical notes](/home/knielsen/code/factoring-connection/codex/PHASE1_NOTES.md), [status](/home/knielsen/code/factoring-connection/codex/STATUS.md).

## 4. Proved scaffolding and precisely delimited obstructions

### 4.1 Random monic quotients supply formal PIT without factoring

**Status: proved; known specialization, not a novelty claim.**

Let \(P\in R_N[Y]\) be represented by an arithmetic circuit, with degree at most \(D\ge1\). Choose \(e\ge2\) with \(2^e\ge4D\), and sample a uniform monic degree-\(e\) polynomial \(h\) over \(R_N\).

If \(P\ne0\), then
\[
\Pr(P\bmod h\ne0)\ge\frac1{4e}.
\]

**Proof.** Choose a prime power \(p^a\mid N\) for which \(P\) is nonzero. Locally, write
\[
P=p^vQ,\qquad v<a,\qquad \overline Q\ne0\in\mathbb F_p[Y].
\]
This choice is used only in the proof.

Reduction modulo a monic polynomial is linear and commutes with coefficient reduction. Hence a nonzero coefficient of \(Q\bmod h\) modulo \(p\) survives multiplication by \(p^v\) modulo \(p^a\).

The reduction \(\bar h\) is uniform among the \(p^e\) monic degree-\(e\) polynomials over \(\mathbb F_p\). At least \(p^e/(2e)\) are irreducible. At most \(D/e\) distinct irreducibles of degree \(e\) divide \(\overline Q\). Thus
\[
\Pr(P\bmod h\ne0)
\ge\frac1{2e}-\frac{D}{ep^e}
\ge\frac1{4e}.
\]
No irreducibility test is required. Circuit evaluation in \(R_N[Y]/h\) uses only monic reduction. Repetition gives polynomial-time amplification. \(\square\)

The useful consequence for A1 is precise: **efficient ordinary-circuit section construction would suffice; a separate unknown randomized formal zero-test lemma is unnecessary.**

### 4.2 A binary section observable requires exponentially many linear states

**Status: proved; obstruction only for constant-coefficient linear realizations.**

For \(0\le r<M\), define
\[
s_t=[X^r]\bigl((1+X)^t\bmod(X^M-c)\bigr).
\]
Then
\[
\boxed{\sum_{t\ge0}s_tz^t
=\frac{z^r(1-z)^{M-1-r}}{(1-z)^M-cz^M}}.
\]

**Proof.** In the quotient algebra,
\[
\sum_{t\ge0}(1+X)^tz^t
=\frac1{1-z-zX}.
\]
Multiply numerator and denominator by
\[
\sum_{j=0}^{M-1}(zX)^j(1-z)^{M-1-j}.
\]
Using \(X^M=c\), the denominator becomes \((1-z)^M-cz^M\). Extracting coefficient \(X^r\) gives the formula. \(\square\)

For \(M=2^d,c=2\), over every odd field the denominator has degree \(M\). Its values at \(z=0,1\) are nonzero, so it is coprime to the numerator. The reduced rational generating function therefore has denominator degree \(M\).

Any constant-coefficient linear realization producing the entire sequence has dimension at least \(M=2^d\). Nested extension coordinates cannot conceal this cost: a rank-\(M\) state over a degree-\(e\) algebra has \(Me\) base coordinates.

This does not rule out nonlinear circuits or a special operation computing only the single time \(t=N\).

### 4.3 Rank extraction has a perfect coefficient-pairing formulation

**Status: proved; standard algebra used as an extraction interface.**

For monic \(h\) of degree \(r\), let
\[
H_{ij}=[X^{r-1}]\bigl(X^{i+j}f\bmod h\bigr),
\qquad 0\le i,j<r.
\]
Then \(H\) has the same field rank as multiplication by \(f\).

**Proof.** Define
\[
G_{ij}=[X^{r-1}](X^{i+j}\bmod h).
\]
Entries below the anti-diagonal in degree order are zero, and anti-diagonal entries are one. Thus \(\det G=\pm1\). Moreover,
\[
H=GB_f.
\]
Therefore their ranks agree. For \(h=X^r-c\), \(G\) is exactly the anti-diagonal permutation matrix.

Over a field, if \(g=\gcd(f,h)\), multiplication by \(f\) has a kernel of dimension \(\deg g\), giving
\[
\operatorname{rank}B_f=r-\deg g.
\]
Finally, if elimination over \(R_{pq}\) encounters only units and zeros, every pivot contributes equally to both CRT ranks. Unequal final ranks therefore force a nontrivial nonunit. \(\square\)

### 4.4 Public cyclotomic roots explain the first sketch successes

**Status: proved; elementary scaffolding, novelty not asserted.**

For \(N\ge3\), the only cyclotomic factors of the integer polynomial \(\widehat F_N\) are:

- \(\Phi_2=X+1\), exactly for odd \(N\), with multiplicity one.
- \(\Phi_3=X^2+X+1\), for \(N\equiv1,5\pmod6\), with multiplicities two and one, respectively.

**Proof.** The assertions at \(-1\) follow by evaluation and differentiation. At a primitive cube root \(\omega\), use \(1+\omega=-\omega^2\). Direct substitution gives precisely the stated residue classes. In the double-root case,
\[
\widehat F_N''(\omega)=N(N-1)\ne0,
\]
so the multiplicity is exactly two.

For \(m\ge4\), take \(\zeta=e^{2\pi i/m}\). A zero would imply
\[
(2\cos(\pi/m))^N=|1+\zeta^N|\le2.
\]
But its left side is at least \((\sqrt2)^3>2\). Finally, \(\widehat F_N(1)=2^N-2\ne0\). \(\square\)

Thus \(W=X(X+1)\Phi_3\) captures the public zero sets relevant to the repair, without their multiplicities.

For \(N=35\), the obstruction is stronger:
\[
F_{35}\bmod5=2^5[X(X+1)\Phi_3^2]^5,
\qquad
F_{35}\bmod7=5^7W^7.
\]
If \(W\) is a unit modulo any monic quotient \(h\), then \(F_{35}\) is a unit in both component algebras. Their multiplication ranks are both \(\deg h\).

This rules out rank mismatch after filtering at \(N=35\); it does not rule out nonunit entries of the resulting invertible matrices.

### 4.5 Removing the derivative’s public root can remove all scalar signal

**Status: proved; refutes the stripped scalar model when \(d=2\).**

For distinct odd primes \(p,q\), let \(N=pq\) and \(d=\gcd(p-1,q-1)\). Then
\[
H_N=(2X+1)Q_N,
\]
where \(Q_N\) has an \(O(\log N)\)-size divide-free circuit and exactly \(d-2\) roots in each hidden field.

**Proof.** Write \(N-1=2m\), \(A=(1+X)^2\), \(B=X^2\). Then
\[
Q_N=\sum_{j=0}^{m-1}A^{m-1-j}B^j.
\]
Binary geometric-sum composition computes this with \(O(\log m)\) gates; for example,
\[
G_{2u}=(A^u+B^u)G_u.
\]
This is an identity over \(\mathbb Z\), requiring no division.

The substitution in D1 gives \(d-1\) roots for \(H_N\). The root \(-1/2\) is simple because, with \(e=N-1\),
\[
H_N'(-1/2)=2e(1/2)^{e-1}\ne0
\]
modulo either factor: \(e\equiv-1\) in both characteristics. Dividing by \(2X+1\) therefore removes exactly one root in each field. \(\square\)

When \(d=2\), every \(Q_N(a)\) is a unit modulo \(N\), proving
\[
\gcd(H_N(a),N)=\gcd(2a+1,N).
\]

### 4.6 The CRT mask density is exact, but applies only to an average interval

**Status: proved; useful calibration for E1.**

For odd \(\ell\nmid N\),
\[
|S_\ell|=\frac{\ell+\chi_\ell(N)}2.
\]

**Proof.** The equation
\[
y^2=s^2-4N
\]
is equivalent to
\[
(s-y)(s+y)=4N.
\]
There are \(\ell-1\) ordered nonzero factor pairs, and division by \(2\) is valid. Hence
\[
\sum_s\chi_\ell(s^2-4N)=-1.
\]
There are \(1+\chi_\ell(N)\) zero discriminants. Including zero among the squares gives the formula.

CRT multiplies the local densities. Averaging a length-\(W\) interval over all translations modulo \(L\), each admissible residue is counted exactly \(W\) times, giving mean survivor count \(W\rho\). \(\square\)

The proof supplies no bound for the planted interval.

### 4.7 The quadratic-torus probability reduces exactly to shifted-order correlation

**Status: proved; restricted-family guarantee.**

Using the notation of F1,
\[
m=N-\epsilon_p\epsilon_q
\equiv\epsilon_p(q-\epsilon_q)=\epsilon_pB\pmod A.
\]
Therefore \(\gcd(m,A)=\gcd(A,B)=d\), and symmetrically for \(B\).

The map \(\alpha\mapsto\alpha/\bar\alpha\) is a surjective homomorphism onto each local norm-one group, with equal-sized fibers. Uniform unit sampling therefore produces uniform group elements. Exactly \(d\) elements in either cyclic group are killed by \(m\). CRT independence now gives
\[
\Pr(\text{exactly one identity})
=\frac dA\left(1-\frac dB\right)
+\frac dB\left(1-\frac dA\right).
\]

This identifies the precise missing information for any repair: an efficiently constructed operation must alter the unfavorable local order correlation, rather than merely describe the same rank-two algebra more compactly.

=== BRIEF ===
# Research agent brief: find a simple classical connection for integer factorization

**Prepared:** 17 October 2026  
**Mode:** Active mathematical research, with exact computational experiments  
**Primary target:** A classical polynomial-time factorization algorithm, or a concrete new mechanism that materially advances toward one

## 1. Mission

Investigate whether a simple connection between explicit integer arithmetic and compressed algebraic representations can produce an efficient classical factorization algorithm.

The working hypothesis is that the eventual mechanism may be conceptually simple and may not resemble sieving or classical simulation of a quantum algorithm. Treat this as a productive hypothesis, not an established fact. Do not dismiss a direction because factoring has been studied for a long time. Equally, do not mistake a short identity, a short program describing an enormous object, or successful toy examples for an efficient algorithm.

**Do actual research.** Read relevant primary sources, formulate candidate mechanisms, derive consequences, implement the most testable constructions, find counterexamples, and revise. A literature survey or a statement that the problem is unresolved is not the requested deliverable.

Begin with the polynomial lead below, but do not become attached to it. You may replace it with a better representation if you explain what new operation becomes possible.

The central question is:

> Which factor-revealing feature has a compact exact representation, and what operation would let us extract that feature without expanding the representation?

### Problem and computational model

Input is an integer `N > 1`, in binary. Write `n = N.bit_length()` and `R_N = Z/NZ`.

The initial hard family is `N = p*q`, with independently generated, distinct primes of similar bit lengths. The factors are not given to the candidate algorithm. Also test prime powers, unbalanced factors, and composites with more than two prime factors.

The first algorithmic goal is to return any `d` with `1 < d < N` and `N % d == 0`. Recursive splitting and primality testing turn a sufficiently uniform splitting algorithm into full factorization; account for those costs when making a final claim.

Use ordinary classical bit complexity. Deterministic polynomial time or randomized expected/bounded-error polynomial time would satisfy the main objective; state which one a proposal targets. Count preprocessing, repetitions, memory, bit precision, coefficient lengths, and representation conversion. A polynomial bound in `N`, `p`, a polynomial degree near `sqrt(N)`, or a group order is not a polynomial bound in `n`.

A restricted-family improvement, a new constructive lemma, a useful reduction with a smaller missing primitive, or a carefully delimited obstruction is worthwhile intermediate progress. Label it accurately and use it to choose the next experiment.

## 2. Mathematical starting point

The elementary claims in this section are supplied with proofs. They are research scaffolding, **not claims of novelty**. Check them independently and search for prior formulations before attributing originality.

### 2.1 The polynomial and its first surviving coefficient

Define

\[
F_N(X)=(1+X)^N-1-X^N\in R_N[X].
\]

The prime/composite polynomial identity is the starting point of the AKS primality test; see [AKS], Lemma 2.1 and Section 2. AKS uses suitable low-degree quotient rings to obtain a primality test. It does not supply the factor-extraction operation sought here.

Let `N` be composite and let `p` be its smallest prime divisor. For `1 <= k < p`,

\[
\binom Nk=\frac{N(N-1)\cdots(N-k+1)}{k!}\equiv0\pmod N,
\]

because `k!` is invertible modulo `N`.

At `k = p`, use the exact integer identity

\[
\binom Np=\frac Np\binom{N-1}{p-1}.
\]

All denominators in the product for the second binomial coefficient are units modulo `N`, so

\[
\binom{N-1}{p-1}
=\prod_{j=1}^{p-1}\frac{N-j}{j}
\equiv(-1)^{p-1}\pmod N.
\]

For odd `p`, the sign is positive. For `p = 2`, `-N/2` and `N/2` are the same residue modulo `N`. Therefore the least nonnegative residue is exactly

\[
\boxed{[X^p]F_N=N/p\pmod N.}
\]

Consequently,

\[
\boxed{\min\{k:[X^k]F_N\ne0\}=p.}
\]

This holds for every composite `N`, including prime powers, not just squarefree semiprimes. Since `p < N`, the removed endpoint `X^N` causes no issue. For prime `N`, `F_N` is the zero polynomial instead.

Examples:

| N | First surviving degree p | Coefficient residue N/p |
|---|---:|---:|
| 9 | 3 | 3 |
| 15 | 3 | 5 |
| 35 | 5 | 7 |
| 45 | 3 | 15 |
| 49 | 7 | 7 |
| 143 | 11 | 13 |

The expression for `F_N` has an arithmetic circuit with `O(n)` gates using repeated squaring. This is a statement about circuit description size, not expanded polynomial size or the cost of coefficient extraction.

### 2.2 A weaker extraction task already suffices: any surviving index

Write `a_k = binom(N,k) mod N` for an interior coefficient. The integer identity

\[
k\binom Nk=N\binom{N-1}{k-1}
\]

gives

\[
k a_k=0\quad\text{in }R_N.
\]

Thus, for `0 < k < N`,

\[
a_k\ne0\Longrightarrow1<\gcd(k,N)<N.
\]

**We do not need the first coefficient, or even its value. Finding any index of a nonzero coefficient gives a factor immediately.**

Target the weaker problem whenever it permits a simpler construction:

> Given the compact circuit for `F_N`, find one exponent in its nonzero support without enumerating the support.

### 2.3 Why individual coefficient queries do not provide the missing information

Suppose an algorithm asks for `a_k` at an explicitly chosen interior index `k`.

If `gcd(k,N) = 1`, its answer is already known to be zero. Otherwise, computing that GCD has already found a nontrivial factor.

Therefore an individual-coefficient oracle can be simulated up to factor discovery: compute the GCD first; either return zero or terminate with a factor. This does not prove factoring hard. It shows that **merely making coefficient evaluation fast is insufficient for this polynomial**. The useful operation must act on many indices at once, reveal a support index, or exploit different arithmetic information.

This is a particularly important filter for proposals based on random coefficient sampling, a binomial-coefficient library, or a fast high-order derivative at one guessed index.

### 2.4 Prefix testing and its factorial connection

For `K >= 1`, define

\[
Z(N,K):=[F_N(X)\equiv0\pmod{N,X^K}].
\]

For composite `N`,

\[
Z(N,K)\Longleftrightarrow K\le p
\Longleftrightarrow\gcd\bigl(N,(K-1)!\bigr)=1.
\]

For finding `p`, it suffices to use `2 <= K <= floor(sqrt(N))+1`. A uniform polynomial-time prefix predicate would recover `p` by `O(n)` binary-search queries. A randomized predicate must have sufficiently amplified error for the adaptive search.

This is an exact connection to bounded divisibility and modular factorials. It is not, by itself, a speedup. Investigate whether the polynomial representation makes the unit/nonunit predicate cheaper than evaluating the factorial, and whether factorial or recurrence techniques suggest a better compression. Read [BGS] and the related classical constructions [Harvey], [HH].

Do not require full factorial evaluation unnecessarily: the predicate asks only whether it is a unit modulo `N`. Conversely, do not treat an efficient unit predicate as already available.

### 2.5 Ordinary differentiation loses the signal

\[
F_N'(X)=N(1+X)^{N-1}-NX^{N-1}=0\quad\text{in }R_N[X].
\]

Every positive-order ordinary derivative is therefore zero. A characteristic-zero method based on the logarithmic derivative or derivative ratios cannot simply be transferred to this ring.

Hasse derivatives retain coefficient information:

\[
D^{[k]}F_N(0)=[X^k]F_N.
\]

But evaluating them at individually selected orders has the limitation in Section 2.3. A potentially useful derivative construction must aggregate orders, preserve additional lifted information, or access a different invariant.

### 2.6 Scalar evaluations and sparse recovery have specific traps

A nonzero polynomial over `R_N` can define the zero function on all scalar residues. For example, `F_561(a) = 0 mod 561` for every scalar `a`: modulo each of `3, 11, 17`, exponentiation to `561` acts as the identity. Yet `F_561` is not the zero polynomial, by Section 2.1.

Therefore scalar random evaluation in `R_N` is not automatically a sound polynomial identity test. Distinguish formal polynomial identity from equality as functions on a finite ring.

Nor is the expanded support polynomially small in general. For `N = p*q`, `p < q`, reduction modulo `q` gives

\[
F_N(X)\equiv(1+X^q)^p-1-X^{pq}\pmod q.
\]

For every `1 <= j < p`, `binom(p,j)` is nonzero modulo `q`, because its prime divisors are at most `p`. Thus `F_N` has at least `p-1` nonzero coefficients. On balanced semiprimes this is exponentially many in `n`.

The relevant opportunity is structured or partial support recovery, not an unsupported assumption of ordinary sparsity.

## 3. Priority research tracks

These are starting directions, not a compulsory architecture. First investigate A and B, connect them to C and D, and keep one genuinely independent alternative alive. The proposed modifications below are hypotheses, not results established by the cited papers.

### Track A — Binary coefficient sections and compressed witness extraction

For any polynomial `P`, define

\[
E(P)(Y)=\sum_j[X^{2j}]P\,Y^j,
\qquad
O(P)(Y)=\sum_j[X^{2j+1}]P\,Y^j.
\]

Then

\[
P(X)=E(P)(X^2)+X\,O(P)(X^2).
\]

If `P` is nonzero, at least one section is nonzero. Repeatedly selecting a nonzero section reconstructs the bits of a support index in `O(log degree(P))` levels.

There are elementary gate rules. Writing `f=f_e(X^2)+Xf_o(X^2)` and similarly for `g`,

\[
E(fg)=f_e g_e+Yf_o g_o,
\qquad
O(fg)=f_e g_o+f_o g_e.
\]

Addition is componentwise. A single section transformation has constant-factor circuit overhead with appropriate sharing. Repeating it `O(n)` times need not: a constant-factor increase at every level can still become exponential.

**Precise conditional target:** Construct a representation for all sections reached along one search path that remains polynomial in `n`, together with a polynomial-time formal zero test for those representations over `R_N`. This would turn the support-index reduction into a factorizer.

Tasks:

1. Derive and implement section transformations on arithmetic circuits, preserving shared subexpressions. Measure growth in unique nodes, not printed expression length alone.
2. Exploit the special repeated-power structure of `(1+X)^N`. Investigate matrix recurrences, rational generating functions, and reusable formulas for sections; do not default to a generic circuit expansion.
3. Determine what information a zero test actually needs. Test whether a weaker certificate or a small family of quotient evaluations can select a nonzero child without reconstructing either child.
4. Study coefficientwise products with interval or bit-pattern masks. The mask may have a short description, but prove or measure the cost of applying it.
5. Look for a class closed under both sectioning and zero testing. It need only contain the required descendants, not arbitrary arithmetic circuits.

Record separately the cost of constructing a section, representing it, evaluating it, and testing its formal identity. Making only one of these cheap is partial progress, not the full algorithm.

**First experiment:** Follow both diagnostic paths chosen with known support information and operational paths chosen without factors. Diagnostic paths measure representation growth only; they are not factoring results.

### Track B — Quotient sketches that retain useful exponent information

Compute compact projections such as

\[
F_N(X)\bmod(X^r-1),\qquad
F_N(X)\bmod(X^r-c),
\]

with polynomially bounded `r` and cheaply chosen `c`. Repeated squaring evaluates each projection without expanding degree `N`.

Modulo `X^r-1`, a coefficient at position `j` is the sum of original coefficients at `j, j+r, j+2r, ...`. A unit-valued sum can mix contributions divisible by different hidden prime factors. A nonzero projected coefficient is therefore not automatically a factor-revealing coefficient of the original polynomial.

Modulo `X^r-c`, folded terms acquire weights `c^t` according to the quotient of their exponent by `r`. Investigate whether several weights or several quotient rings expose a useful exponent, a nonunit, or enough partial information for another method.

Concrete questions:

- Can paired sketches isolate one hidden characteristic's contribution without knowing that characteristic?
- Can small-degree resultants, subresultants, annihilators, ranks, or failed inversions reveal a nontrivial GCD?
- Can a small family of weighted projections identify one support exponent modulo enough small moduli without enumerating exponentially many matching choices?
- Can the relation `F_N' = 0` constrain the information discarded by reduction? Differentiation before and after reduction need not commute: retain the quotient terms in the derivation.
- Can a sketch answer a prefix or binary-section question rather than reconstruct the full polynomial?

Read [AKS] for the distinction between proving compositeness and extracting a witness, and [Sparse] for interpolation mechanisms and their assumptions. Characteristic, number of terms, roots of unity, discrete-log steps, and field-versus-ring behavior must all be audited before transfer.

Do not use generic sketch collisions to claim an impossibility theorem for this family: the decoder receives `N` and can exploit its explicit arithmetic. A negative result must identify the restricted sketch/decoder model it actually addresses.

**First experiment:** Exhaustively inspect several small families of weighted projections, their coefficient GCDs, and basic invariants on tiny semiprimes. Derive the reason for any successful pattern before scaling it.

### Track C — Factorial unit tests and compressed recurrence evaluation

The prefix predicate is precisely a factorial-unit predicate. This links the polynomial lead to algorithms for products of consecutive terms and linear recurrences with polynomial coefficients [BGS].

Investigate whether unit/nonunit information admits a simpler recurrence than the factorial value itself. Candidate tools include product trees, block products, multipoint evaluation, cancellation-aware recurrences, and transformations that isolate an interval containing a divisor.

Every proposed recurrence should state whether its dimension, polynomial degree, or number of blocks grows with `K` rather than `log K`. If a recurrence uses division modulo `N`, compute the GCD of each denominator before inversion: a nonunit can be useful output, not merely an implementation failure.

A particularly relevant new lead is Yann Tal's September 2026 preprint [Tal]. Its abstract describes modular-factorial improvements through Jacobi sums, reconstruction, and recursive reductions. **Its setting supplies a known prime modulus, or specified prime-power settings, and appropriate divisor conditions.** It is not a general unknown-factor-modulus algorithm. Inspect the full construction and ask which identities survive without knowing the prime or its group structure.

A useful adaptation must avoid assuming a factorization of `N`, its totient, auxiliary quantities of comparable unknown difficulty, or a known prime underlying `R_N`. Character sums are useful here only if the proposed way to construct or evaluate them is explicit and affordable.

**First experiment:** Implement a transparent block-product baseline, then one recurrence or invariant that attempts to answer only the unit predicate. Compare actual arithmetic counts and largest representation sizes.

### Track D — Hidden characteristic, lifted information, and annihilators

For `N = p*q`, the polynomial has different Frobenius structure in its two hidden components:

\[
F_N\bmod p=(1+X^p)^q-1-X^{pq},
\]

\[
F_N\bmod q=(1+X^q)^p-1-X^{pq}.
\]

The goal is not to perform these reductions using known factors. It is to find an operation over `R_N` that separates the two structures without knowing them.

Investigate Hasse derivatives, divided-power constructions, formal-series sections, finite differences, recurrence singularities, and annihilating ideals. Consider modest lifts such as modulus `N^t` where `t` is polynomially bounded, but derive what survives the lift and how it is extracted. Ordinary differentiation becomes zero after reduction modulo `N`; differentiation before reduction or in a lift retains other data, whose usefulness must be demonstrated rather than assumed.

A promising failure mode of a field algorithm is encountering a nonunit pivot. Instrument such algorithms to return `gcd(pivot,N)`. However, selecting a pivot guaranteed to become a useful nonunit is part of the research problem.

Avoid circular use of Lucas/Kummer digits in base `p`, a `p`-Frobenius map, or Cartier operators parameterized by an unknown factor. They can guide identities in the proof, but an operational algorithm must construct its own parameters from `N`.

**First experiment:** Build small lifted or sectioned representations and classify the denominators, pivots, and annihilators that arise. Separate structural factor discovery from accidental divisibility by a small tested integer.

### Track E — Transfer partial information into lattices or hidden relations

A sketch does not have to reveal a complete factor. Determine whether it yields a rigorously specified interval, congruence, or correlated set of bits that existing small-root methods can use. Coppersmith's work [Coppersmith] is the relevant reference, but extract the precise theorem and parameter requirements before applying it.

For hidden-relation approaches, use [Regev], [Pilatte], and [EG] to identify the information that classical postprocessing consumes. A short-relation existence theorem is not a procedure for finding those relations; a decoder is not a sampler.

One exact interface is multiplicative characters. With generators `a_i`, define

\[
\mathcal L=\{z\in\mathbb Z^d:\prod_i a_i^{z_i}=1\pmod N\}.
\]

If a computable character satisfies `chi(a_i)=exp(2*pi*i*w_i)`, then

\[
\langle w,z\rangle\in\mathbb Z\quad\text{for every }z\in\mathcal L.
\]

The question is whether explicit arithmetic supplies useful character phases or partial phase constraints efficiently. For squared generators `a_i=b_i^2`, ambient characters taking values only in `{+1,-1}` become trivial. A proposed source of phases must avoid that degeneracy and avoid constructing the unknown group decomposition first.

Also inspect [Gao] for the neighboring classical idea of deliberately avoiding uninformative shortest vectors. Transfer a concrete objective or construction, not merely the instruction to use lattice reduction.

Treat [Cleve] as a model-specific warning about generic oracle order finding, not as a lower bound for explicit integer arithmetic.

**First experiment:** Specify an exact partial-information format and demonstrate its downstream usefulness with a diagnostic oracle. Then separately attempt to produce it from `N` alone. Never combine those two stages into an unsupported end-to-end speed claim.

### Track F — One independent connection

Reserve some effort for a route not suggested by the polynomial lead. Examples include a different low-description-length generating function, a nonunit-producing algebraic identity, or an arithmetic invariant of a low-dimensional object.

The proposal must specify an executable operation and why it sees the hidden factors. Merely restating `pq=N` as optimization, SAT, a spectral problem, or a dynamical system is not enough. Show where the representation gives a cheaper operation or a smaller search problem.

## 4. Research procedure

### Phase 1: establish the mathematical and literature baseline

Reproduce the elementary lemmas and tests in this brief. Read the core sources needed for the two most promising tracks; do not read every listed paper before doing mathematics.

Create a compact map of the nearest prior art. For each relevant result record its exact problem, input representation, theorem assumptions, complexity variables, and the primitive it provides. Use primary sources. Verify title, authors, date, and theorem statement; earlier conversational summaries are leads, not evidence.

Keep discovery work separate from verification. In the first generation pass, propose at least six mechanisms, including one outside the supplied tracks. Then select two or three with the clearest operation and smallest missing lemma.

### Phase 2: write a candidate card before substantial coding

Use this template for each serious candidate:

```text
ID and one-sentence mechanism:
Target problem / supported input family:
Representation and how it is built from N:
Exact operation that extracts information:
Why that information exposes a factor:
Known lemmas being used:
New lemma required:
Correctness argument or conjecture:
Cost in n, including construction and extraction:
Randomness and success probability:
Closest prior method and the actual difference:
Smallest informative experiment:
Counterexample or scaling behavior that would refute this version:
Possible repair if that happens:
```

Reject or repair proposals that name the missing operation without implementing or bounding it. Do not abandon a promising representation merely because the first naive implementation expands it.

### Phase 3: alternate construction and criticism

For each selected mechanism:

1. Derive the strongest statement you can justify.
2. Implement the smallest exact experiment that distinguishes it from the obvious baseline.
3. Search for boundary cases and counterexamples.
4. Audit representation growth and total work.
5. Produce a revised mechanism or a narrow obstruction identifying what a revision must change.

Use a separate mathematical critic and complexity auditor if real subagents are available. Give the critic the candidate and code, not a persuasive narrative about expected success. Without subagents, perform separately documented review passes. Do not claim independent review unless it actually occurred.

The critic must state which assumption fails and give an example, derivation, or missing bound. “Factoring is believed hard” is not a technical critique.

Do not stop after a survey or after one failed proposal. Complete at least one repair cycle on the best candidate, unless a complete proof and independently checked implementation are already obtained.

### Phase 4: keep the result reproducible

Preserve proofs, candidate versions, failed constructions, seeds, raw measurements, and exact commands. A negative result about one representation is useful when its scope is explicit and it points to a different operation to try.

Allocate most effort to mathematical construction and discriminating experiments, not ever-expanding bibliographies or GPU optimization. Use exact integer/ring arithmetic initially. Floating-point exploration is allowed only with a later precision analysis.

## 5. Experimental design and leakage controls

### Input families

Begin with exhaustive small composites for algebraic correctness. Then use reproducibly generated balanced semiprimes at increasing bit lengths, for example `16, 24, 32, 48, 64`, extending only as resources permit. Successful factoring of small inputs is a correctness check, not evidence of polynomial scaling.

Include independent-prime balanced semiprimes, prime powers, composites with several distinct factors, highly unbalanced semiprimes, Carmichael numbers, and prime inputs. Near-equal primes, special algebraic forms, and deliberately smooth `p-1` or `q-1` belong in separately labeled diagnostic families.

### Separate algorithm inputs from diagnostic truth

The candidate receives only `N`, public configuration, and a random seed. A separate test harness may retain factors to compute labels and understand failures.

Never leak factors through filenames, seeds deterministically encoding factors, candidate-specific parameter choices, precomputed group orders, smoothness-selected instances, library factorization calls, or a simulator initialized with the factors.

Use actual random generation of the intended family and record the sampling procedure. Freeze candidate choices before the holdout evaluation. Special test families are useful but must not silently replace the intended distribution.

### Record meaningful costs

For each run record at least:

- Correctness, result, seed, wall time, timeout, and peak memory.
- Ring operations and GCD calls, with operand bit lengths.
- Largest polynomial degree actually materialized and number of stored coefficients.
- Circuit nodes, recurrence dimension, matrix size, section depth, and rank where relevant.
- Number of sketches, candidate combinations, restarts, and rejected trials.
- Numerical precision and error bounds for any approximate component.

Separate preprocessing, information acquisition, decoding, and verification. Count all work, including unsuccessful attempts. Distinguish parallel depth from total classical work.

Compare against transparent trial division and Pollard-rho implementations on small tests, a relevant product/block baseline for Track C, and the nearest prior construction for the actual proposed operation. Use controlled same-language arithmetic comparisons where possible. Do not claim superiority to GNFS from a small-input benchmark or an improvement over trial division.

Empirical fits are diagnostics, not complexity proofs. Show raw counters versus both `n` and the actual structural parameters. In particular, a method polynomial in a hidden factor can look innocuous on tiny bit lengths.

## 6. Mandatory complexity and correctness audit

For every proposed breakthrough, answer the following before calling it an algorithmic improvement.

**Representation.** Is the object genuinely stored in polynomial space, or merely described by a short program that later expands? Are all required operations efficient on that representation?

**Coefficient access.** Does the method require sequential coefficients up to `p`, a dense truncation at `K`, or exponentially many terms after sectioning? Section 2.3 rules out individual coefficient evaluation as the missing advantage here.

**Field assumptions.** Does a theorem require a known prime, an algebraically closed field, a characteristic bound, or invertible denominators? Prove the needed ring version or report factor-revealing nonunits.

**Hidden inputs.** Is any step using `p`, `q`, `phi(N)`, `lambda(N)`, an unknown order, a hidden lattice basis, or informative samples supplied by an oracle? State exactly how it is obtained.

**Numerical cost.** An error tolerance `2^(-poly(n))` uses polynomially many bits. An error tolerance `2^(-p)` needs order-`p` bits. A tiny error written with a short exponent is not automatically cheap to achieve. Conversely, do not call every exponentially small error an exponential-precision requirement.

**State-space cost.** A short diffusion rule, tensor description, or matrix-power expression is not enough. Bound the represented state, ranks, contractions, and the cost of the final observable.

**Degree and root bounds.** Check polynomial identity-test assumptions when formal degree is exponential. Over a finite ring, nonzero polynomials can vanish at every scalar point. Auxiliary prime moduli also change coefficient divisibility: reducing modulo an unrelated prime does not preserve vanishing modulo `N` for free.

**Randomness.** Bound the success probability and cost per attempt, including constructing the distribution and rejecting bad outputs. Distinguish worst-case, average-case, and restricted-family claims.

**Certificates.** A candidate factor is cheaply verifiable. That does not bound search time or certify an average-case/polynomial-time theorem. A decoder's correctness with ideal samples does not validate the source of those samples.

**Novelty.** Compare with the closest specific method. A new notation, a different implementation of a known identity, or a known reduction rediscovered here is not an asymptotic breakthrough. It can still be useful scaffolding for the next modification.

Use these checks to improve proposals, not to replace proposal generation.

## 7. Workspace and required deliverables

Suggested layout:

```text
factoring-connection/
  README.md
  STATUS.md
  references.md
  notes/
    seed_lemmas.md
    prior_art_map.md
    candidate_001.md
    candidate_002.md
    review_and_repairs.md
  src/
    seed_checks.py
    ring_polynomials.py
    circuit_sections.py
    quotient_sketches.py
    candidates/
  tests/
  experiments/
    configs/
    raw/
    analysis/
  REPORT.md
  NEXT_STEPS.md
```

Use Python with exact arithmetic for the first prototypes; use existing computer-algebra packages where appropriate, but inspect any operation that may internally factor a modulus. Optimize only after identifying an informative mechanism. Do not buy cloud resources, publish results, or contact authors without authorization.

A candidate interface should accept `N` and public configuration only, return a factor or explicit failure, and emit an auditable metrics record. Reference factorization functions belong in the harness, never inside the purported new factorizer.

The final report must contain:

1. **Best surviving mechanism:** its precise new operation and the strongest supported claim.
2. **Mathematics:** proved lemmas, conjectures, counterexamples, and any remaining proof obligation, clearly separated.
3. **Experiments:** reproducible commands, raw results, input distributions, failures, and scaling diagnostics.
4. **Comparison:** closest prior art and exactly what changed.
5. **Continuation:** the next most decisive mathematical question or experiment, with a concrete proposed attack.

Use status labels such as `proved`, `experimentally supported`, `conjectured`, `refuted for this model`, and `not yet tested`. A useful report may contain a novel mechanism plus a sharply identified unresolved lemma. Do not inflate that into a completed factorizer, and do not stop at “no breakthrough found” when the work supports a more specific next step.

## 8. Primary-source reading list

These are selected starting points, not a claim of an exhaustive literature search. Bibliographic/abstract pages were checked when preparing this brief; the AKS seed identity was also checked in the paper. Except for the elementary derivations above, do not assume that every cited theorem or implementation has been independently verified here. Open the relevant full papers, record the version consulted, and verify exact conditions before relying on them.

**[AKS]** Manindra Agrawal, Neeraj Kayal, Nitin Saxena. *PRIMES is in P*. Annals of Mathematics 160 (2004), 781–793. [Paper](https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/primality_journal.pdf). Focus on Lemma 2.1 and the quotient-ring reduction in Section 2; distinguish primality detection from factor extraction.

**[BGS]** Alin Bostan, Pierrick Gaudry, Éric Schost. *Linear Recurrences with Polynomial Coefficients and Application to Integer Factorization and Cartier–Manin Operator*. SIAM Journal on Computing 36(6) (2007), 1777–1806. [Author-hosted paper](https://mathexp.eu/bostan/publications/BoGaSc07.pdf). This directly connects distant recurrence terms, product computations, and classical factorization. Audit the dependence on the recurrence index.

**[Harvey]** David Harvey. *An exponent one-fifth algorithm for deterministic integer factorisation*. arXiv:2010.05450 (2020). [Record](https://arxiv.org/abs/2010.05450). A rigorous `N^(1/5+o(1))` classical factoring construction; a useful structured-search baseline, not polynomial in input bit length.

**[HH]** David Harvey, Markus Hittmeir. *A log-log speedup for exponent one-fifth deterministic integer factorisation*. arXiv:2105.11105 (2021). [Record](https://arxiv.org/abs/2105.11105). Study the arithmetic structure and the actual batching operation, not just its asymptotic headline.

**[Sparse]** Joris van der Hoeven, Grégoire Lecerf. *Sparse polynomial interpolation: Faster strategies over finite fields*. Author-hosted final version dated 2 May 2024. [Text](https://www.texmacs.org/joris/ffsparse/ffsparse.html). Relevant to projections, support reconstruction, and circuit access. Its field, sparsity, characteristic, and root-of-unity assumptions require explicit checking for our setting.

**[Tal]** Yann Tal. *Computing Modular Factorials Below the Square-Root Barrier*. ECCC TR26-211, September 2026 preprint. [Record and paper access](https://eccc.weizmann.ac.il/report/2026/211/). A recent lead involving known prime/prime-power moduli, Jacobi sums, and recursive factorial computations. Treat the claims as preprint results to inspect, and preserve the known-modulus/divisor conditions. No general unknown-factor-modulus speedup is asserted by this brief.

**[Coppersmith]** Don Coppersmith. *Small Solutions to Polynomial Equations, and Low Exponent RSA Vulnerabilities*. Journal of Cryptology 10 (1997), 233–260. [Publisher record](https://link.springer.com/article/10.1007/s001459900030). Use to specify exactly what partial factor information would suffice, rather than casually assuming that arbitrary correlated bits are usable.

**[Regev]** Oded Regev. *An Efficient Quantum Factoring Algorithm*. arXiv:2308.06572 (2023; consult current version). [Record](https://arxiv.org/abs/2308.06572). A source for the hidden-relation formulation and the separation between quantum information acquisition and classical postprocessing.

**[Pilatte]** Cédric Pilatte. *Unconditional correctness of recent quantum algorithms for factoring and computing discrete logarithms*. arXiv:2404.16450; Forum of Mathematics, Pi 14 (2026), e5. [Record](https://arxiv.org/abs/2404.16450). Structural existence/correctness results for suitable modified constructions; do not equate existence of short relations with constructive classical access to them.

**[EG]** Martin Ekerå, Joel Gärtner. *Extending Regev's factoring algorithm to compute discrete logarithms*. arXiv:2311.05545 (2023; consult current version). [Record](https://arxiv.org/abs/2311.05545). Relevant to what approximate lattice information suffices for classical recovery. Extract the exact noise, distribution, and precision requirements.

**[Gao]** Yiming Gao, Yansong Feng, Honggang Hu, Yanbin Pan. *On Factoring and Power Divisor Problems via Rank-3 Lattices and the Second Vector*. arXiv:2512.19076 (December 2025 preprint). [Record](https://arxiv.org/abs/2512.19076). A concrete neighboring classical lead emphasizing useful rather than merely shortest vectors. Reproduce and inspect the construction before generalizing it.

**[Cleve]** Richard Cleve. *The query complexity of order-finding*. arXiv:quant-ph/9911124 (1999). [Record](https://arxiv.org/abs/quant-ph/9911124). Read the oracle model carefully. Its lower-bound conclusions do not automatically apply to algorithms exploiting the explicit representation of modular integers.

For additional searches, combine terms such as “arithmetic circuit coefficient sections,” “order of vanishing straight-line programs,” “polynomial identity testing composite characteristic,” “Hadamard product compressed polynomials,” “factorial unit test composite modulus,” “holonomic recurrence distant terms,” and “constructive AKS factor extraction.” Treat these as search directions, not as claims that a particular theorem already exists.

## 9. Small exact smoke test

Save this as `src/seed_checks.py`. It is reference/testing code using trial division and explicit coefficients. It is deliberately not a fast factoring algorithm and must not appear as the hidden implementation of a candidate.

```python
"""Exact tests for the elementary seed lemmas; small instances only."""
from math import comb, gcd, isqrt


def least_divisor_reference(value: int) -> int:
    if value < 2:
        raise ValueError("value must be at least 2")
    for divisor in range(2, isqrt(value) + 1):
        if value % divisor == 0:
            return divisor
    return value


def check_seed_lemmas(limit: int = 200) -> int:
    if not 4 <= limit <= 2000:
        raise ValueError("use a small test limit between 4 and 2000")
    composites = 0
    for modulus in range(2, limit + 1):
        first = least_divisor_reference(modulus)
        # Only interior coefficients: the endpoints of F_N are removed.
        coeffs = [comb(modulus, k) % modulus
                  for k in range(1, modulus)]
        if first == modulus:
            assert not any(coeffs), (modulus, "prime identity")
            continue

        composites += 1
        assert not any(coeffs[:first - 1])
        assert coeffs[first - 1] == modulus // first
        for k, coefficient in enumerate(coeffs, start=1):
            assert k * coefficient % modulus == 0
            if coefficient:
                assert 1 < gcd(k, modulus) < modulus

        # product = (K-1)! mod modulus at each iteration.
        product = 1
        for cutoff in range(1, isqrt(modulus) + 2):
            prefix_zero = not any(coeffs[:cutoff - 1])
            assert prefix_zero == (cutoff <= first)
            assert prefix_zero == (gcd(product, modulus) == 1)
            product = product * cutoff % modulus

    for point in range(561):
        assert (pow(1 + point, 561, 561)
                - 1 - pow(point, 561, 561)) % 561 == 0
    return composites


if __name__ == "__main__":
    count = check_seed_lemmas()
    print(f"Seed checks passed for {count} small composites.")
```

During preparation, the first-index/value and coefficient-annihilator identities were also checked directly for all 831 composites up to 1,000, and the scalar-zero property of `F_561` was checked at all 561 residues. These checks supplement the proofs; they are not evidence of a fast extraction procedure.

## 10. Start here

Read Sections 1–2, run the smoke tests, and create `STATUS.md`. Then investigate binary sectioning and weighted quotient sketches in parallel or as separate passes. Produce explicit candidate cards, implement the most discriminating small experiment for each, and perform a repair cycle after criticism.

The aim is to find the connection, not merely describe why existing approaches are expensive. Keep asking what information a representation exposes cheaply that was not cheaply accessible before, and make that operation concrete.
