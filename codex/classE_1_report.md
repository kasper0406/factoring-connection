| Question | Result | Status / limitation |
|---|---|---|
| Q1: factor-bearing convergents | A proper denominator divisible by \(p\) must equal \(p\), with Euclidean remainder \(p\). A congruence and an additional cylinder inequality give an exact iff. | **Proved.** The proposed small-ratio congruence alone is insufficient. |
| Q1: square-root argument | For a unit \(t=\lfloor\sqrt N\rfloor\), every proper convergent denominator and positive remainder is coprime to \(N\). | **Proved.** The nonunit case is exactly \(q=p+2\), where \(t=p\) already factors \(N\). |
| Q1: two ends | The reversed digit word represents \(\min(\bar a,N-\bar a)/N\). Two prefixes of denominator at most \(\sqrt N\) can leave a single junction digit. | **Proved.** Its cost need not be bounded; smaller thresholds can leave a growing middle word. |
| Q2: uniform argument | The entire one-ended prefix through denominator \(T\) has a universal reference law within \(O(T^2/N)\) in total variation for balanced semiprimes. | **Proved.** At \(T=N^{1/2-\delta}\), error \(O(N^{-2\delta})\), including arbitrary postprocessing of that prefix. |
| Q2: joint/full path | Endpoint rectangle discrepancy is \(O(N^{-1/2}\log^2N)\); an elementary joint-prefix bound is \(O(T^4N^{-1/2}\log^2N)\). | **Proved bounds.** Fine endpoint paths cannot generally be modeled as independent continuous endpoints. An arbitrary full class-E floor is **open**. |
| Q2: structured arguments | Target-conditional endpoint equidistribution gives a floor when its resolution controls the event sets of the readout. | **Conditional theorem.** Qualitative equidistribution alone is insufficient. |
| Q3: local sums and congruences | Explicit public congruences for \(12M F_2\) and \(24M F_3\) give a DSF decoder for the unweighted local sum under a nondegeneracy condition. | **Proved here.** Root-order local parts have different weights; some local parts are identically zero. |
| Q3: three-factor computation | An exact three-variable weighted lattice count evaluates the sum in polynomial time in the bit length. | **Proved reduction** to the established Barvinok algorithm, with pseudocode. |
| Q4: reduced coefficients | Indefinite reduction is a quadratic-irrational CF path with unbounded state. Arbitrary coefficient costs/composition paths do not thereby become finite-state class E plus genus. | Classical correspondence **proved**; general classification and a new nonquadratic leak **open**. |
| Q5: battery | Structured endpoint/junction features, shifted cone sums, and coefficient readouts, with geometry-aware floors and one-sided certificates. | **Proposed**, not run. No scalable new leak established. |

# Blind report: continued fractions at structured arguments

**Scope.** Background context came from the specified DS/DSF, class-E, GAR, class-H, and open-class sections of newfamily_results.md, leakmap_results.md, and genusleak_results.md. No notes/classE_* file was read. Computation was restricted to small exact checks in this scratch workspace. “Open” below means unresolved in this report, not a lower bound.

Two qualifications matter. The integer \(N\) determines its factors mathematically, so a prediction theorem must restrict the readout and identify a reference baseline. Also, if “real size” means the exact real number \(N\), the floor is tautological. Operationally it must mean specified size information or a specified smooth size baseline, together with fixed congruences. The distributional theorems below are stated independently of this ambiguity; converting them to target floors requires a separate bound for the baseline.

## Q1. Exact structure and factor-bearing conditions

### 1.1 Notation — proved identities

Use the canonical expansion

\[
\frac aN=[0;b_1,\ldots,b_r],\qquad b_r\ge2,\quad 0<a<N,\quad(a,N)=1.
\]

Write convergents \(u_i/v_i\), with \((u_{-1},u_0)=(1,0)\), \((v_{-1},v_0)=(0,1)\). Set \(R_{-1}=N,R_0=a\) and \(R_i=R_{i-2}-b_iR_{i-1}\). Thus \(R_{r-1}=1,R_r=0,v_r=N\). Induction gives

\[
av_i-Nu_i=(-1)^iR_i,\qquad
v_iR_{i-1}+v_{i-1}R_i=N. \tag{1}
\]

For every proper positive-index convergent,

\[
v_iR_i<N, \tag{2}
\]

because \(R_{i-1}>R_i\). This is also the usual estimate \(|a/N-u_i/v_i|<1/v_i^2\).

### 1.2 \(a=\lfloor\sqrt N\rfloor\) — proved

Put \(t=\lfloor\sqrt N\rfloor\), \(D=N-t^2=jt+d\), \(0\le d<t\). Here \(1\le D\le2t\), \(j\in\{0,1,2\}\), and the unit condition forces \(d>0\). The exact path is

\[
\boxed{\frac tN=[0;t+j,\ \text{the digit word of }t/d].} \tag{3}
\]

Indeed, \(N=(t+j)t+d\), followed by Euclid on \(t,d\).

For \(z=\sqrt N=t+\eta\), \(0<\eta<1\),

\[
\eta=\frac D{z+t},\qquad
\sqrt N=[t;\lfloor(z+t)/D\rfloor,\ldots]. \tag{4}
\]

This quadratic recursion differs from the rational recursion using \(D\bmod t\). Two exact approximation identities explain apparent similarities:

\[
\frac tN=\frac1{\sqrt N}-\frac\eta N,\qquad
t+\frac Nt=2\sqrt N+\frac{\eta^2}{t}. \tag{5}
\]

The second expression is twice the Newton approximation to \(\sqrt N\), within \(1/t\) of \(\sqrt{4N}\). Neither estimate guarantees a growing common CF prefix: distance to cylinder boundaries also matters. For example,

\[
N=53\cdot59=3127,\quad t=55,\quad
55/3127=[0;56,1,5,1,7],
\]

while \(\sqrt N\) and \(\sqrt{4N}\) begin respectively

\[
[55;1,11,2,3,\ldots],\qquad [111;1,5,4,1,\ldots].
\]

Thus identifying these paths by a universal shift or truncation is false. Equations (3)–(5) are exact relations without such a claim.

In factor-dependent notation, let \(h=(p+q)/2,e=(q-p)/2\). Then

\[
N=h^2-e^2,\qquad D=2t(h-t)+(h-t)^2-e^2. \tag{6}
\]

This introduces the unknown Fermat representation; it does not extract it. Also, every nonzero convergent of \(t/N\) is at most \(1/\lfloor N/t\rfloor=O(N^{-1/2})\), whereas \(p/q>1/2\). None approximates \(p/q\) accurately.

**Stronger gcd obstruction.** For odd distinct primes in the stated range,

\[
(t,N)>1\ \Longleftrightarrow\ q=p+2,
\]

and then \(t=p\). Proof: \(p\le t<q<2p\), so a nontrivial gcd requires \(t=p\); \(pq<(p+1)^2\), together with the even prime gap, is equivalent to \(q=p+2\).

In particular, this nonunit case has \(N\equiv3\pmod4\); it never occurs on the specified \(\epsilon\)-target population \(N\equiv1\pmod4\).

In every remaining case, **all proper convergent denominators and all positive Euclidean remainders for \(t/N\) are coprime to \(N\)**. The theorem in §1.7 proves this: a factor-bearing denominator would require \(t\equiv\pm1\pmod q\), whereas \(p<t\le q-2\). This excludes direct gcd exposure, not statistical leakage from all digit costs.

There is nevertheless a shared quadratic datum: \(D\) is the first complete-quotient denominator in the usual \(\sqrt N\) recursion. Since \(t^2\equiv-D\pmod N\), the rational inverse obeys \(\bar t\equiv-tD^{-1}\pmod N\). For its convergents, (1) also gives the exact norm identity

\[
R_i^2+Dv_i^2
=N\bigl((v_i-tu_i)^2+Du_i^2\bigr).
\]

At unit \(v_i\), this restates the already public constraint \((-D/p)=(-D/q)=1\). It does not turn each remainder into the principal-cycle norm \(Q_i\) of \(\sqrt N\), or give the individual remainder-symbol constraints of that path.

### 1.3 \(a=\lfloor N^{2/3}\rfloor\) — proved

Let \(w=N^{1/3}\), \(A=\lfloor w^2\rfloor\), \(B=\lfloor w\rfloor\), \(C=\lfloor N/A\rfloor\), \(E=N-CA\). For \(N\ge8\), \(C\in\{B,B+1\}\): writing \(A=w^2-\xi\), \(0\le\xi<1\), gives

\[
0\le N/A-w=\frac{w\xi}{w^2-\xi}<1.
\]

For a unit argument, \(E>0\), and

\[
\frac AN=[0;C,\ \text{the digit word of }A/E]. \tag{7}
\]

This is a public Euclidean pair, not a factor-related quadratic norm identity. No analogue of the square-root gcd obstruction is proved for this argument.

### 1.4 \(a=\lfloor N/k\rfloor\pm c\) — proved for fixed \(k,c\)

Assume \(k\ge2\), \(0<a<N\). Write \(N=kt+s\), \(0\le s<k\), \(a=t+e\), \(e=\pm c\). Then

\[
N=ka+d,\qquad d=s-ke=O_{k,c}(1).
\]

* For \(0<d<a\): first digit \(k\), then Euclid on \(a,d\).
* For \(d=-h<0\), \(a>2h\): first digits \(k-1,1\), then the digit word of \((a-h)/h\).
* For \(d=0\): the unit requirement fails except for trivial small cases.

After a single large quotient, the remaining tail has bounded modulus and bounded length. Its data depend on congruences modulo constants depending on \(k,c\); the large quotient is an affine floor expression in \(N\). This is a structural simplification, not automatically a target floor for arbitrary \(F\): unrestricted arithmetic on the large digit may retain essentially the public input.

### 1.5 \(a=\lfloor\omega N\rfloor\) — proved structure; growing-prefix claims need hypotheses

For a fixed computable \(0<\omega<1\),

\[
\frac{\lfloor\omega N\rfloor}N=\omega-\frac{\{\omega N\}}N. \tag{8}
\]

A word remains a prefix exactly when this rational stays in its CF cylinder. Every fixed prefix of an irrational \(\omega\) eventually occurs. No uniform growing-prefix assertion holds for arbitrary \(\omega\).

If all partial quotients of \(\omega\) are bounded by \(A\), its position within each cylinder stays an \(A\)-dependent positive proportion away from both endpoints. Cylinder lengths are comparable to \(v^{-2}\). Thus, for some \(c_A>0\), every prefix of \(\omega\) with denominator at most \(c_A\sqrt N\) remains a prefix in (8): its boundary distance exceeds \(1/N\). In particular the golden argument has a long substantially deterministic forward word, not a random initial word.

For rational \(\omega=h/k\) in lowest terms, put \(r=hN\bmod k\). Then

\[
ka=hN-r,\qquad a/N=h/k-r/(kN). \tag{9}
\]

For large \(N\), \(r>0\), the rational \(h/k\) is a convergent by the \(1/(2k^2)\) criterion. There is an appropriate side expansion of \(h/k\), a large digit, and a bounded-modulus tail. At \(r=0\) the unit requirement usually fails. A fixed real must also come with an algorithm computing its floor at the required precision; polynomial-time computability is not automatic.

### 1.6 Exponent residues — proved CRT formulas

For \(a=2^{(N-1)/2}\bmod N\),

\[
a_p=\left(\frac2p\right)2^{(q-1)/2}\pmod p,\qquad
a_q=\left(\frac2q\right)2^{(p-1)/2}\pmod q. \tag{10}
\]

Indeed,

\[
(pq-1)/2=(p-1)(q-1)/2+(p-1)/2+(q-1)/2,
\]

and Euler's criterion applies. For \(b=2^{N-1}\bmod N\),

\[
b_p=2^{q-1}\pmod p,\quad b_q=2^{p-1}\pmod q,\quad b=a^2\pmod N. \tag{11}
\]

These generally are **not** CRT combinations of two signs. The half-exponent components are signs precisely when \(b=1\pmod N\). Opposite signs then factor by \(\gcd(a\pm1,N)\); matching signs give \(a=1\) or \(N-1\). More generally, \(\gcd(b-1,N)\) factors if exactly one component in (11) equals \(1\). For instance \(b_p=1\) means \(\operatorname{ord}_p(2)\mid q-1\).

After CRT and taking the integer representative, the CF is exactly the ordinary Euclidean path in (1). Formulas (10)–(11) describe its factor-dependent inputs, but do not give a short quotient-word formula or an equidistribution theorem over the specified semiprimes.

### 1.7 Exact factor-bearing denominator criterion — proved

Exclude the final denominator \(v_r=N\), which trivially contains both factors.

**Theorem 1.** For \(p<q<2p\), \((a,pq)=1\), a proper convergent denominator divisible by \(p\) exists iff there are \(\sigma\in\{\pm1\}\), \(1\le u\le p-1\), such that

\[
\boxed{a=qu+\sigma,\quad
w=(-\sigma u^{-1}\bmod p)\in\{1,\ldots,p-1\},\quad w<q-p.} \tag{12}
\]

Its denominator and remainder both equal \(p\). No proper denominator is divisible by \(q\).

**Necessity.** If \(v_i=ps\), (1) forces \(R_i=pt\), \(s,t>0\). By (2),

\[
p^2st<N,\qquad st<q/p<2,
\]

so \(s=t=1\). Equation (1) gives \(a=qu+\sigma\), \(\sigma=(-1)^i\). The preceding denominator \(w\) satisfies \(uw\equiv-\sigma\pmod p\). Equation (1)'s second identity gives \(q=R_{i-1}+w>p+w\). The same product bound excludes a proper denominator divisible by \(q\).

**Sufficiency.** A reduced \(u/v\in(0,1)\), on specified side \(\sigma\), has a unique side predecessor denominator \(w\in\{1,\ldots,v-1\}\) with \(uw\equiv-\sigma\pmod v\). Its continuation cylinder extends from \(u/v\) to the adjacent mediant on that side, with length \(1/[v(v+w)]\). In (12),

\[
a/N-u/p=\sigma/N.
\]

It lies strictly in that cylinder iff \(1/N<1/[p(p+w)]\), exactly \(w<q-p\). Equality cannot occur here: the determinant of the candidate and \(a/N\) is \(p>1\), whereas a boundary mediant would have determinant \(1\). ∎

**The proposed congruence is not sufficient.** Take \(p=29,q=31,N=899,a=63=2q+1\). Both primes have their top two bits set. Although \(a\equiv1\pmod q\), with product \(rs=1\), its side predecessor is \(w=-2^{-1}\bmod29=14\not<2\). The proper denominators are \(14,43,57,157,371\), none divisible by \(29\).

Outside the balanced range, an exact general version is available. For a candidate reduced \(u/(ps)\in(0,1)\), \(ps<N\), write

\[
as-qu=\sigma t,\quad t>0,\quad
w=-\sigma u^{-1}\bmod ps,\quad 1\le w<ps.
\]

It is a convergent iff

\[
\boxed{t(ps+w)<q.} \tag{13}
\]

Thus \(st<q/p\) and \(a\equiv\sigma t/s\pmod q\) are necessary. The denominator \(s<q\) is invertible modulo \(q\). The product bound alone omits the side-cylinder condition. Equality in (13) is impossible for these candidates because the remaining determinant is \(pt>1\).

**Exact uniform-unit probability.** For each sign, inversion permutes \(u=1,\ldots,p-1\), giving \(q-p-1\) acceptable values. The two families are disjoint, hence

\[
\Pr_a[\text{proper factor-bearing denominator}]
=\frac{2(q-p-1)}{(p-1)(q-1)}=O(N^{-1/2}). \tag{14}
\]

Also \(\gcd(v_i,N)=\gcd(R_i,N)\), from (1) and \((a,N)=1\). Positive-remainder gcd exposure is precisely the same event.

### 1.8 Approximate factors and Coppersmith — proved implications

A denominator merely of order \(p\) gives no useful precision. A certified estimate \(P_0\) with

\[
|P_0-p|<N^{1/4-\eta}
\]

for fixed \(\eta>0\) is different: \(P_0+X\) has a small root modulo the unknown divisor \(p=N^{1/2+o(1)}\), recoverable by the unknown-divisor Coppersmith method. A polynomial-size list containing one such estimate also suffices. This is the standard \(N^{\beta^2}\) bound for a divisor of size \(N^\beta\); see [Coppersmith's account, integer factorization section](https://cr.yp.to/bib/2001/coppersmith.pdf). Exponent slack avoids endpoint constants. It is not a lower bound.

Consequently any proposed “near-\(p\)” mechanism should report absolute error on the \(N^{1/4}\) scale. Exact divisibility already factors by gcd.

As worded, “any computable \(a(N)\)” has a trivial positive answer: factor first and take \(a=q-1\). It is a unit and \(1/p\) is a convergent, since \(w=1<q-p\). The substantive question is a construction that avoids prior factoring. A polynomial-time routine always producing a unit argument with a factor-bearing proper denominator would be a polynomial-time factoring algorithm: scan its \(O(\log N)\) convergents and take gcds. None is supplied here.

For a convergent to approximate \(p/q\) accurately, the input \(a/N\) must itself be close to that ratio, except for coarse early convergents. Being a rational CF supplies no such factor-ratio approximation.

### 1.9 The two-ended decomposition — proved exactly

Let \(A(b)=\begin{pmatrix}b&1\\1&0\end{pmatrix}\). Then

\[
A(b_1)\cdots A(b_r)=
\begin{pmatrix}N&v_{r-1}\\a&u_{r-1}\end{pmatrix}.
\]

Reversal transposes this product, so

\[
\boxed{[0;b_r,\ldots,b_1]=v_{r-1}/N,\quad
v_{r-1}\equiv(-1)^{r-1}\bar a\pmod N.} \tag{15}
\]

Since \(b_r\ge2\), \(v_{r-1}<N/2\), and therefore \(v_{r-1}=\min(\bar a,N-\bar a)\). If \(b_1=1\), merge the reversed word's final \(1\) with its preceding digit to canonicalize it. The affected proper denominator exceeds \(N/2\), so small-denominator reverse prefixes are unaffected.

Write \(w_j=K(b_r,\ldots,b_{r-j+1})\) for reversed continuants, \(w_0=1,w_{-1}=0\). At a split after digit \(k\),

\[
\boxed{N=v_k w_{r-k}+v_{k-1}w_{r-k-1}.} \tag{16}
\]

Choose \(k\) maximal with \(v_k\le\sqrt N\). Applying (16) at \(k+1\) shows \(w_{r-k-1}<\sqrt N\). Thus the word consists of:

1. a forward prefix of denominator at most \(\sqrt N\);
2. the one digit \(b_{k+1}\);
3. the reversal of a prefix of (15), also of denominator below \(\sqrt N\).

This is deterministic, not independence. A single digit can have a very large cost. At threshold \(N^{1/2-\delta}\), a word of ones shows that the uncovered middle can have length proportional to \(\delta\log N\). For state-dependent costs, reversing the suffix also requires transport of matrix residues using the terminal matrix and suffix matrix modulo \(m_0\), together with parity. Bare reversal of the costs without transporting states is incorrect.

### 1.10 GAR and visibility of the hidden phase

GAR gives, with \(d=x_p-x_q\),

\[
\operatorname{rep}(x)/N\equiv d\theta+x_q/N\pmod1,\qquad
\theta=(q^{-1}\bmod p)/p. \tag{17}
\]

Every CF readout is a function of that expression. If \(d=0\), this particular representative has no phase dependence. If \(d\ne0\), CF cylinder crossings can respond to \(\theta\), but \(d,x_q\) are unknown and correlated with it.

Holding those quantities fixed, a cylinder of denominator \(V\) has typical scale \(V^{-2}\), corresponding locally to phase resolution \(1/(|d|V^2)\), away from its boundaries. Multiplication by \(d\) on the circle also introduces branches. This is a resolution observation, **not a phase-recovery algorithm**. If \(x_q\) were known as an integer and \(d\not\equiv0\pmod p\), then \(\gcd(\operatorname{rep}(x)-x_q,N)=q\) already factors. Treating local CRT components as free known coefficients is circular.

Write \(N\theta=1+jp\). A phase approximation of error \(\Delta\) gives an approximate multiple of \(p\) with error \(O(N\Delta+1)\). The same small-root theorem factors when \(\Delta<N^{-3/4-\eta}\), for any fixed \(\eta>0\), explaining the \(3n/4\) precision scale with endpoint qualifications. No theorem here says that every useful class-E bit predictor must recover that precision, or any explicit approximation to \(\theta\).

## Q2. Floor theorems with precise quantifiers

### 2.1 A one-ended total variation theorem — proved

Let \(\mathcal P_T(x)\) be the CF word stopped at the last convergent denominator at most \(T\), **excluding the next, overshooting digit**. For rational \(a/N\), assume \(T<N\). The empty prefix is allowed.

**Theorem 2.** For \(a\) uniform among units modulo \(N\), and \(U\) uniform on \((0,1)\),

\[
d_{\rm TV}\bigl(\mathcal P_T(a/N),\mathcal P_T(U)\bigr)
\ll\frac{\tau(N)T^2}{\varphi(N)}. \tag{18}
\]

For balanced semiprimes this is \(O(T^2/N)\), and hence \(O(N^{-2\delta})\) at \(T=N^{1/2-\delta}\). Every function of that prefix, including its matrices modulo fixed \(m_0\), inherits the bound by data processing. This includes probabilities of arbitrary residue-class or discontinuous output events.

**Proof.** Each stopping word ends at a reduced rational \(u/v\), \(v\le T\). There are at most two representations of each such rational, so \(O(T^2)\) possible cells. For a prefix ending with denominators \(w,v\), the next denominator exceeds \(T\) iff

\[
b_{\rm next}\ge L=\left\lfloor\frac{T-w}{v}\right\rfloor+1.
\]

The union of these continuation cylinders is an interval, from \(u/v\) to the side mediant with continuation digit \(L\). There are \(O(T^2)\) boundary points. For any interval \(I\), Möbius inversion gives

\[
\#\{1\le a<N:(a,N)=1,\ a/N\in I\}
=\varphi(N)|I|+O(\tau(N)). \tag{19}
\]

Sum the normalized discrepancies over the partition and handle boundary atoms with the same estimate. ∎

The theorem also permits fixed numerator-residue labels. If \((m,N)=1\), CRT and the same Möbius argument give

\[
\#\{a/N\in I:(a,N)=1,\ a\equiv c\pmod m\}
=\frac{\varphi(N)}m|I|+O(\tau(N)).
\]

Thus the joint prefix and \(a\bmod m\) law is within \(O(m\tau(N)T^2/\varphi(N))\) of the continuous-prefix law with an independent uniform label. Given \(N\bmod m\), (1) then determines every remainder residue from that label and the convergent matrices. This covers the intended fixed-modulus remainder counts too. For fixed \(m\), sufficiently large balanced primes are automatically coprime to \(m\); a small shared factor is detected directly.

Inversion permutes units, so the same result holds for a prefix of \(\bar a/N\). For the **literal reversed path**, use the folded variable \(\min(U,1-U)\), with density \(2\) on \((0,1/2)\), as reference. Folding changes constants only. The overshooting digit is deliberately excluded: it can reveal much finer information than the stopped prefix.

This is a theorem about the entire prefix law, stronger than a mean-value statement. It does not include the junction, the full path cost, or the Dedekind endpoint correction involving the inverse.

### 2.2 Joint endpoint equidistribution and the gcd caveat — proved from Weil

Set

\[
\mu_N=\frac1{\varphi(N)}\sum_{a\in(\mathbb Z/N)^\times}
\delta_{(a/N,\bar a/N)}.
\]

Its Fourier coefficient at \((h,k)\) is \(S(h,k;N)/\varphi(N)\). CRT gives

\[
S(h,k;pq)=S(hq^{-1},kq^{-1};p)\,
S(hp^{-1},kp^{-1};q), \tag{20}
\]

with the inverses taken modulo the respective primes. At a prime, the factor is \(p-1\) when both arguments vanish, \(-1\) when exactly one vanishes, and bounded by \(2\sqrt p\) otherwise. The last assertion is [Weil's Kloosterman bound](https://www.renyi.hu/~gharcos/weil.pdf). Consequently,

\[
|S(h,k;N)|\le4\sqrt{N\gcd(h,k,N)}. \tag{21}
\]

One must not drop the gcd factor at arbitrary frequencies: for \(p\mid h,k\), the normalized coefficient can have scale \(N^{-1/4}\), instead of \(N^{-1/2}\).

For rectangle discrepancy, choose the Erdős–Turán–Koksma cutoff \(H=c\sqrt N<p\), uniformly possible because \(p>\sqrt{N/2}\). Nonzero frequencies in this box avoid simultaneous divisibility by either prime. Summing the Fourier weights gives

\[
\boxed{\sup_{I,J}|\mu_N(I\times J)-|I||J||
\ll N^{-1/2}\log^2N.} \tag{22}
\]

A disjoint union of \(R\) rectangles has error at most \(O(RN^{-1/2}\log^2N)\); folding the inverse coordinate changes constants only.

Both prefix partitions have \(O(T^2)\) cells, so for independent uniform \(U,V\),

\[
d_{\rm TV}\left(
(\mathcal P_T(a/N),\mathcal P_T(\bar a/N)),
(\mathcal P_T(U),\mathcal P_T(V))\right)
\ll T^4N^{-1/2}\log^2N. \tag{23}
\]

This explicit sufficient bound works at \(T\le N^{1/8-\eta}\), with error \(O(N^{-4\eta}\log^2N)\). At \(T=N^{1/2-\delta}\), this proof gives exponent \(c=4\delta-3/2>0\), requiring \(\delta>3/8\). No optimality is claimed.

### 2.3 A precision obstruction for two fine endpoint paths — proved

If \(T\gg N^{1/4}\) by a power, the joint prefix law is asymptotically at total variation distance **one** from the independent continuous reference in (23). This is not merely a weak discrepancy estimate.

A stopping cell with last denominator \(v\) has length

\[
\frac1{v(Lv+w)}\le\frac1{vT}.
\]

The combined measure of cells with \(v<\epsilon T\) is \(O(\epsilon+T^{-1})\), summing over at most \(2\varphi(v)\) words per denominator. All remaining cells have length at most \(1/(\epsilon T^2)\). The discrete pair law occupies at most \(\varphi(N)\le N\) cell pairs. The independent continuous measure of their union is bounded by

\[
O(\epsilon+T^{-1})+\frac{N}{\epsilon^2T^4}.
\]

When \(T^4/N\) grows by a power, choose \(\epsilon\to0\) with the last expression tending to zero. The union has discrete probability one, proving the assertion. The folded version is identical up to constants.

Thus for \(0<\delta<1/4\), two \(N^{1/2-\delta}\)-prefixes cannot jointly be replaced by independent real CF prefixes in total variation. This does not prove factor-bit leakage: their coupling includes a public determinant constraint. It does show why ordinary endpoint equidistribution cannot dismiss every fine-path readout.

### 2.4 What the fixed-denominator literature actually proves — verified

**Mean length.** With \(\kappa=12\log2/\pi^2\), Porter's theorem is

\[
\frac1{\varphi(N)}\sum_{(a,N)=1}r(a/N)
=\kappa\log N+C_P-1+O_\varepsilon(N^{-1/6+\varepsilon}), \tag{24}
\]

where \(C_P\) is Porter's constant. Bykovskii–Frolenkov, Theorem 2, improve the unnormalized error to \(O_\varepsilon(N^{5/6-7/174+\varepsilon})\). For balanced semiprimes, division by \(\varphi(N)\asymp N\) gives

\[
\kappa\log N+C_P-1+O_\varepsilon(N^{-6/29+\varepsilon}). \tag{25}
\]

Their Theorem 1 gives normalized error \(O_{\varepsilon,x}(N^{-4/21+\varepsilon})\) for the mean Gauss–Kuzmin tail-count statistic at fixed \(x\in(0,1]\). These exponents and Porter's universal constant are verified in [Bykovskii–Frolenkov, pp. 645–646](https://www.mathnet.ru/php/getFT.phtml?jrnid=sm&paperid=8718&what=fullteng). These are mean estimates, not arbitrary-law comparisons.

The same pages record Ustinov's predecessor bound for that tail count: unnormalized error \(O_{\varepsilon,x}(N^{5/6}(\log N)^{7/6+\varepsilon})\), or normalized error \(O_{\varepsilon,x}(N^{-1/6}(\log N)^{7/6+\varepsilon})\) for balanced semiprimes.

**Concentration and tails.** For \(S_1=\sum b_i\), Aistleitner–Borda–Hauke prove, for fixed \(C>0\), all \(N\ge3\), \(0<t\le(\log N)^C\),

\[
\Pr\left[\left|S_1-\frac{12}{\pi^2}\log N\log\log N\right|
\ge t\log N\right]\ll_Ct^{-1}, \tag{26}
\]

and

\[
\Pr[\max b_i\ge t\log N]\le\frac{12}{\pi^2t}
+O_C\left(\frac{(\log\log N)^3}{t\log N}\right). \tag{27}
\]

See [Theorems 1 and 3](https://arxiv.org/pdf/2210.14095). They also treat Dedekind sums and other digit statistics. These do not assert polynomial total variation comparison of full laws; (26) itself has only logarithmic-strength consequences for polylogarithmic \(t\).

**Polynomial large deviations.** David–Kim–Mor–Shapira, Theorem 1.2, prove: for each fixed \(\epsilon>0\) and fixed finite word \(\mathbf w\), there exist positive \(\alpha_\epsilon,\alpha_{\epsilon,\mathbf w}\) such that for all sufficiently large fixed \(N\),

\[
\Pr[|r/\log N-\kappa|>\epsilon]<N^{-\alpha_\epsilon},
\quad
\Pr[|D_{\mathbf w}(a/N)-D_{\mathbf w}|>\epsilon]
<N^{-\alpha_{\epsilon,\mathbf w}}. \tag{28}
\]

Here \(D_{\mathbf w}(a/N)\) is its overlapping word frequency. The paper explicitly says its methods **do not estimate these exponents**. Numerical exponent values would be unverified. See [Theorem 1.2 and its following paragraph](https://arxiv.org/html/2401.15586v2). Fixed-tolerance concentration does not control parity, unnormalized cost residues, or every bounded postprocessing.

These results do not prove the requested floor for every full class-E functional and every bounded readout. Such a theorem needs a metric, normalization, reference law, and restriction on dependence on \(N\). For the full class this report leaves it **open**. Equations (18), (22), and (23) give unconditional replacements with precise scopes and exponents.

### 2.5 Turning a reference law into a target floor — proved implication

Let \(Y(N)\) be a target, centered if desired, \(C\) the designated public floor, and \(g(C,S)\in[-1,1]\). Suppose a reference law \(\nu_C\) satisfies

\[
\left|\mathbb E_a g(C(N),S_N(a))
-\int g(C(N),s)\,d\nu_{C(N)}(s)\right|\le\rho_N.
\]

Define \(B(C)=\int g(C,s)d\nu_C(s)\). Conditioning on \(N\) gives

\[
\left|\mathbb E[Yg(C,S_N(a))]-\mathbb E[YB(C)]\right|
\le\mathbb E[|Y|\rho_N]. \tag{29}
\]

The baseline target bound must still be supplied independently. Equation (18) gives \(\rho_N=O(N^{-2\delta})\) for every allowed one-ended-prefix rule whose reference baseline belongs to the designated floor. If unrestricted \(N\) is available to the predictor, the baseline may itself factor, making the proposed floor vacuous. Also, a bound for the mean of a raw feature is not a bound for every nonlinear predictor of that feature.

### 2.6 Structured arguments: a quantitative conditional theorem

Take an explicit semiprime ensemble within public stratum \(C\), with target \(B\in\{0,1\}\), and \(Z_N=(a(N)/N,\bar a(N)/N)\).

**Hypothesis H(\(\Delta\)).** For each target subclass \(B=b\), its conditional endpoint law differs on every axis-parallel rectangle from the **same** public reference measure \(\nu_C\) by at most \(\Delta\). Include needed public size parameters or formulate this for the associated public family of rectangles. If a readout uses terminal matrix residues or parity, require the same assertion jointly on those finitely many labels.

Unconditional equidistribution over both targets is insufficient: two biased subclasses can average to an equidistributed population.

**Theorem 3 — conditional.** Suppose each prediction event admits inner and outer approximations, each a disjoint union of at most \(R\) rectangles in endpoint coordinates and finite labels, with reference mass of their difference at most \(\eta\). Under H(\(\Delta\)), its probability in each target subclass differs from reference by at most \(R\Delta+\eta\). The difference between the two subclass probabilities is at most \(2(R\Delta+\eta)\).

**Proof.** Sum the rectangle discrepancies for the inner and outer approximations and sandwich the event. Then subtract the two target-class bounds. ∎

For two \(T\)-prefixes, take \(R=O(T^4)\), \(\eta=0\), with endpoint atoms controlled separately. Thus \(\Delta=O(X^{-\kappa})\), \(T=X^\beta\), \(4\beta<\kappa\) imply a polynomial floor for every function of those prefixes. Fixed coarse partitions or polynomially many rectangles need much less resolution than full paths.

The weakest hypothesis tailored to a battery is equidistribution against just its actual bounded test functions. It suffices but becomes nearly the desired conclusion if that class includes every full class-E readout. The useful noncircular version is quantitative rectangle/cylinder equidistribution **at the resolution and complexity required by the readout**.

For exponent residues, no such target-conditional theorem is proved here. For \(\lfloor\omega N\rfloor\), the first endpoint tends to \(\omega\), so a uniform two-dimensional reference is already false. A Dirac reference at \(\omega\) describes weak convergence only: uniform rectangle discrepancy from it stays one when the endpoint differs from \(\omega\). Instead use appropriately restricted coarse tests, or rescale to \(\omega N-\lfloor\omega N\rfloor\) and check the event-set complexity in those coordinates. The square-root and two-thirds endpoints tend to zero and likewise need rescaling, for example to their fractional-part displacements. A hypothesis with an incorrect uniform reference does not become plausible by conditioning on factor bits.

## Q3. Lattice-cone sums with two and three denominator factors

### 3.1 Normalization and local decomposition — proved

Write

\[
F_d(n;\mathbf a;M)=\frac1M\sum_{k=1}^{M-1}
\frac{e^{2\pi i kn/M}}{\prod_{j=1}^d(1-e^{2\pi i k a_j/M})},
\quad (a_j,M)=1. \tag{30}
\]

This is exactly the question's sign convention; compare the equivalent denominator-sign conventions in [Beck–Diaz–Robins](https://arxiv.org/pdf/math/0204035). If a parameter is a nonunit, first compute its gcd with \(N\); otherwise the displayed root sum can have poles.

Partitioning nontrivial \(N\)-th roots by exact orders \(p,q,pq\) gives

\[
\boxed{F_d(n;\mathbf a;N)
=q^{-1}F_d(n;\mathbf a;p)+p^{-1}F_d(n;\mathbf a;q)
+U_d(n;\mathbf a;N),} \tag{31}
\]

where \(U_d\) sums over primitive \(N\)-th roots with normalization \(1/N\). The weights \(1/q,1/p\) are essential.

By contrast, the **sawtooth-index** identity for Dedekind sums is

\[
s(a,N)=s(a,p)+s(a,q)+U_{\rm saw}(a,N).
\]

Its unweighted local part is not the root-order local part of the cotangent representation \(s(a,N)=(4N)^{-1}\sum\cot(\pi k/N)\cot(\pi ak/N)\). The two unit terms differ. More generally, an unnormalized periodic-kernel sum splits unweighted over multiples of the primes; normalization by \(1/N\) introduces the weights in (31).

There is a general local/unit decomposition, but no normalization-independent theorem that every object called its local part factors \(N\).

### 3.2 Explicit public congruences — proved

All divisions on the right below mean modular inversion of the unit parameters.

**Theorem 4.** For every \(M>1\), units \(a,b,c\), and integer \(n\),

\[
12M F_2(n;a,b;M)\in\mathbb Z,
\]

\[
\boxed{12M F_2(n;a,b;M)\equiv
C_2(n;a,b)=\frac{-a^2-b^2-3ab+6n(a+b)-6n^2}{ab}\pmod M.} \tag{32}
\]

Put \(A=a+b+c,\ B=ab+ac+bc\). Then

\[
24M F_3(n;a,b,c;M)\in\mathbb Z,
\]

\[
\boxed{24M F_3(n;a,b,c;M)\equiv
C_3(n;a,b,c)=
\frac{4n^3-6An^2+2(A^2+B)n-AB}{abc}\pmod M.} \tag{33}
\]

These residues are public aggregates using ring arithmetic and inverses. Bare \(M^kF_d\bmod M\) is not an integer residue until denominators are cleared. For these sums, after multiplying by \(12\) or \(24\), extra powers \(M^k\), \(k\ge2\), give the trivial zero residue.

**Proof.** The one-factor kernel is

\[
f_M(t)=F_1(t;1;M)=
\begin{cases}
(M-1)/(2M),&t\equiv0,\\
(2\operatorname{rep}_M(t)-M-1)/(2M),&t\not\equiv0.
\end{cases} \tag{34}
\]

One can derive it by taking a first difference and using the geometric series. In particular \(2Mf_M(t)\) is integral and

\[
2Mf_M(t)\equiv2t-1\pmod M. \tag{35}
\]

The identity \((1-e^{2\pi it})^{-1}=\tfrac12+\tfrac i2\cot(\pi t)\), with opposite roots paired, gives

\[
F_2(0;a,b;M)=\frac{M-1}{4M}-s(ab^{-1},M), \tag{36}
\]

\[
F_3(0;a,b,c;M)=\frac{M-1}{8M}
-\frac12\{s(ab^{-1},M)+s(ac^{-1},M)+s(bc^{-1},M)\}. \tag{37}
\]

Integrality of \(12Ms(h,M)\) and its classical congruence \(12Ms(h,M)\equiv h+h^{-1}\pmod M\) prove the two assertions at \(n=0\).

For general \(n\), cancel one denominator:

\[
F_2(n;a,b;M)-F_2(n+b;a,b;M)=f_M(a^{-1}n), \tag{38}
\]

\[
F_3(n;a,b,c;M)-F_3(n+c;a,b,c;M)=F_2(n;a,b;M). \tag{39}
\]

Take \(t=b^{-1}n\bmod M\) in \(\{0,\ldots,M-1\}\) and sum (38) along \(0,b,\ldots,(t-1)b\). Multiplying by \(12M\), (35) adds the residue

\[
-6ba^{-1}t(t-1)+6t.
\]

Add (36) and substitute \(t\equiv b^{-1}n\pmod M\) to obtain (32).

Next take \(t=c^{-1}n\bmod M\) and sum (39). Multiplication by \(24M\) adds \(-2\sum_{j=0}^{t-1}C_2(jc;a,b)\), or

\[
2t(3+a/b+b/a)-6c(1/a+1/b)t(t-1)
+2c^2(ab)^{-1}t(t-1)(2t-1).
\]

Substitute \(t\equiv c^{-1}n\), add the residue from (37), and simplify to (33). The integer sum formulas remove all divisions by \(2,3\) before reducing modulo \(M\). Equations (38)–(39) propagate integrality as well. ∎

### 3.3 A DSF decoder, with the precise local oracle and exceptions — proved

For \(d=2,3\), put \(D_2=12,D_3=24\), and

\[
L_d=F_d(n;\mathbf a;p)+F_d(n;\mathbf a;q),
\qquad X=D_dNL_d\in\mathbb Z.
\]

The congruences imply

\[
X\equiv C_d(n;\mathbf a)(p+q)\pmod N. \tag{40}
\]

Indeed \(X=q(D_dpF_d(p))+p(D_dqF_d(q))\); reduce at each prime. If \(\gcd(C_d,N)=1\),

\[
\boxed{p+q=XC_d^{-1}\bmod N,} \tag{41}
\]

with the least nonnegative representative giving the exact sum because \(p+q<N\). Solve the resulting quadratic. If \(\gcd(C_d,N)\) is proper, the public coefficient already factors \(N\). If \(C_d\equiv0\pmod N\), this decoder says nothing.

This is a DSF analogue for the **unweighted** sum of local values. The root-order local part in (31) is weighted and is a different oracle; (41) does not automatically decode it.

Some actual weighted local parts do factor. For example,

\[
F_2(0;1,1;M)=-\frac{(M-1)(M-5)}{12M}.
\]

If \(L_{\rm root}=F_2(p)/q+F_2(q)/p\), \(H=p+q\), then

\[
X_{\rm root}=12NL_{\rm root}=-H^2+2N+6H-10,
\quad H=3+\sqrt{2N-1-X_{\rm root}}.
\]

But a universal assertion without nondegeneracy is false. Changing \(k\mapsto M-k\) in (30) yields

\[
F_d(n;\mathbf a;M)=(-1)^dF_d(a_1+\cdots+a_d-n;\mathbf a;M).
\]

Thus for every odd \(M\),

\[
\boxed{F_3(2;1,1,2;M)=0.} \tag{42}
\]

Both weighted and unweighted local parts are identically zero and carry no nontrivial oracle value.

For generalized Bernoulli/Todd sums, fixed-weight denominator-clearing and Laurent-polynomial corrections are established by [Chae–Jun–Lee](https://arxiv.org/pdf/1406.3563). Such results require the correct weight, cone, and normalization, rather than a blanket substitution into every shifted sum. Equations (32)–(33) avoid an unverified translation of those conventions.

### 3.4 Is the non-aggregate content class E?

For two factors, convolution of (34) gives a one-dimensional weighted floor sum. Euclidean reciprocity computes it; the unshifted case (36) is literally Dedekind plus a public term. See [Beck–Robins' polynomial-time formulas for Dedekind–Rademacher sums](https://arxiv.org/abs/math/0111329). Fixed rational shifts can be transported in finite residue state after clearing their fixed denominators. This proves finiteness of the shift state; a literal class-E representation additionally requires the reciprocity costs to telescope into the prescribed digit/state form, with any endpoint terms retained.

For arbitrary \(n=n(N)\), however, the moving shift can have denominator growing with \(N\). The recursion retains an \(O(N)\)-size shift, not just a state modulo fixed \(m_0\). “Computable by Euclid” alone therefore does not establish membership in the literal narrow class E. Likewise higher-weight two-dimensional cone formulas may depend on unbounded convergent coefficients. A finite-state additive representation needs a separate proof for the selected functional.

For three factors, the exact formula below is a two-variable floor sum, or a three-dimensional lifted cone count. Fixed-dimensional unimodular cone decomposition computes it. That is a higher-dimensional geometric analogue of CF reciprocity, **not a proved reduction to one ordinary CF additive functional**. Polynomial-time computation is neither a floor theorem nor a new factor leak.

Unshifted products of an odd number of cotangents vanish by root pairing. Thus the unshifted three-cotangent Zagier sum is zero; the next nonzero higher unshifted product has four factors. See [Zagier's original definition and parity observation](https://people.mpim-bonn.mpg.de/zagier/files/doi/10.1007/BF01351173/HigherDedekindSums.pdf). Shifts or derivatives can remove that cancellation but change the precise computational class.

### 3.5 Exact polynomial-time three-factor algorithm — proved reduction

For odd \(N\), set \(h=(N-1)/2\). Equation (34) becomes

\[
f_N(t)=\frac{\operatorname{rep}_N(t-1)-h}{N}.
\]

The inverse Fourier transform of a product is cyclic convolution. Substitute first variables \(a(x+1),b(y+1)\), and put

\[
\alpha=\operatorname{rep}_N(-ac^{-1}),\quad
\beta=\operatorname{rep}_N(-bc^{-1}),\quad
\gamma=\operatorname{rep}_N(c^{-1}(n-a-b)-1).
\]

Then

\[
F_3(n;a,b,c;N)=\frac1{N^3}
\sum_{0\le x,y<N}(x-h)(y-h)
\bigl(\operatorname{rep}_N(\gamma+\alpha x+\beta y)-h\bigr). \tag{43}
\]

Replace the representative by \(L-N\lfloor L/N\rfloor\). Since \(\sum_{x=0}^{N-1}(x-h)=0\), the affine part \(L-h\) sums to zero. Hence

\[
\boxed{F_3(n;a,b,c;N)=-\frac1{N^2}
\sum_{0\le x,y<N}(x-h)(y-h)
\left\lfloor\frac{\gamma+\alpha x+\beta y}{N}\right\rfloor.} \tag{44}
\]

This is a weighted lattice count in the bounded polytope

\[
P=\{(x,y,z):0\le x,y\le N-1,\ z\ge1,\
Nz\le\gamma+\alpha x+\beta y\}. \tag{45}
\]

Each integer pair contributes precisely its nonnegative floor number of integer \(z\)'s. The weight is the fixed degree-two polynomial \((x-h)(y-h)\).

~~~text
FourierDedekind3(N, a, b, c, n):
    require N odd and gcd(a*b*c, N) = 1
    ci    := inverse_mod(c, N)
    alpha := (-a*ci) mod N
    beta  := (-b*ci) mod N
    gamma := (ci*(n-a-b)-1) mod N
    h     := (N-1)/2
    P     := {0 <= x,y <= N-1, z >= 1,
              N*z <= gamma + alpha*x + beta*y}
    Z1, Zx, Zy, Zxy := WeightedLatticeSums(P; 1, x, y, x*y)
    H := Zxy - h*(Zx+Zy) + h*h*Z1
    return exact_rational(-H, N*N)

WeightedLatticeSums(P; requested monomials):
    G(X,Y,Z) := BarvinokShortGeneratingFunction(P)
    for each x^i*y^j:
        apply (X*d/dX)^i * (Y*d/dY)^j to G
        evaluate the removable specialization X=Y=Z=1
        by a generic one-parameter Laurent expansion
    return the exact integer sums
~~~

Dimension, inequality count, and weight degree are fixed, and coefficients have \(O(\log N)\) bits after modular reduction. Barvinok produces a polynomial-size rational generating function. A fixed number of Euler derivatives preserves polynomial complexity. A generic exponential substitution avoids zero denominator directions; constant Laurent coefficients evaluate the specialization after pole cancellation. Output bit length is polynomial. This proves exact polynomial-time evaluation without enumerating \(N^2\) pairs or factoring \(N\).

The general theorem is [Barvinok, 1994](https://pubsonline.informs.org/doi/10.1287/moor.19.4.769); polynomial-weight generation by differential operators appears explicitly in [De Loera–Hemmecke–Köppe–Weismantel, Lemma 3](https://www.math.ucdavis.edu/~deloera/RECENT_WORK/mixed-approximation.pdf). A direct modular cone implementation is discussed by [Xin–Xu](https://arxiv.org/pdf/2303.01185); the available draft describes itself as an announcement. The algorithm above relies on the established general theorem, not an assumption that their implementation is complete.

The corresponding two-factor formula is

\[
F_2(n;a,b;N)=\frac1{N^2}\sum_{0\le x<N}(x-h)
\bigl(\operatorname{rep}_N(\gamma+\alpha x)-h\bigr), \tag{46}
\]

with \(\alpha=-ab^{-1}\bmod N,\ \gamma=b^{-1}(n-a)-1\bmod N\). Euclidean weighted floor-sum descent evaluates it. Formula (44) supplies concrete additional geometric readouts, not a reason to identify them with narrow finite-state class E.

**Small verification only.** Exact rational convolution checks on 252 small odd-modulus parameter tuples verified (32), (33), (42), and (44), including composite moduli and moduli divisible by \(3\). Another 125 tuples checked (32)–(33) including even moduli. Exhaustive checks on 43,948 unit pairs for small balanced semiprimes verified (12), (14), and the square-root gcd obstruction. These checks validate signs and formulas; the proofs are the justification. No large battery was run.

## Q4. Reduced-form coefficients beyond genus

### 4.1 What reduces to quadratic continued fractions — classical/proved

For positive nonsquare discriminant \(D\), an indefinite form \((a,b,c)\) determines the quadratic endpoint

\[
\xi=\frac{-b+\sqrt D}{2a}.
\]

Reduction acts by integral Möbius transformations, and a quadratic-irrational CF gives its reduction word, with the standard convention adjustments. In the \(\sqrt{kN}\) notation, complete quotients satisfy

\[
\xi_i=\frac{\sqrt{kN}+P_i}{Q_i},\quad t_i=\lfloor\xi_i\rfloor,\quad
P_{i+1}=t_iQ_i-P_i,\quad
Q_{i+1}=\frac{kN-P_{i+1}^2}{Q_i}. \tag{47}
\]

The form coefficients and complete-quotient states recover one another under the usual conventions. A coefficient statistic along a reduction segment is therefore a statistic of a class-H path **with its unbounded state retained**.

For the principal path, convergents give

\[
A_{i-1}^2-kNB_{i-1}^2=(-1)^iQ_i. \tag{48}
\]

If \(p\nmid Q_i\), then \(p\nmid A_{i-1}\), and

\[
\left(\frac{(-1)^iQ_i}{p}\right)=1,
\]

and similarly at \(q\). This proves the principal-cycle genus restriction. On a nonprincipal class, local admissibility may carry its class character; it is not automatically the same \(+1\) condition.

Square coincidences, congruences of squares assembled from norms, and ambiguous endpoints can expose factors. These are established mechanisms. The form/CF correspondence and square-form algorithm are treated by [Gower–Wagstaff](https://citeseerx.ist.psu.edu/document?doi=1becf939066e33f2e3def14347b81ccd8d22934d&repid=rep1&type=pdf). Their distributional runtime analysis uses assumptions, rather than a theorem that genus restrictions are the only statistical constraints.

### 4.2 Why a general finite-state/genus classification does not follow

Three distinctions block the proposed formal implication:

1. **Unbounded state.** \(P_i,Q_i\), leading coefficients, and infrastructure distances are not determined by a digit and a matrix modulo fixed \(m_0\). Reconstructing them from the complete word and discriminant permits global postprocessing; it does not express every coefficient cost as \(\sum F(b_i,\gamma_i\bmod m_0)\).
2. **Composition history.** Composition makes new coefficients using products, gcds, and congruence solving, then reduces. Individual reduction segments have quadratic CF descriptions, but an arbitrary composition tree retains intermediate integer coefficients and global relations. No fixed-state additive representation follows.
3. **Negative discriminant.** At \(-4kN\), the roots are nonreal. Positive-definite Gauss reduction is modular reduction in the upper half-plane, using translations and swaps depending on the full coefficient state. It is not the real quadratic CF (47), nor automatically a rational \(a/N\) path. Rational Euclidean subroutines inside an algorithm do not prove that all its coefficient statistics belong to class E.

These are failures of a claimed implication, not separation theorems between computational classes. Enlarging “class H” to mean every function of all quadratic-form states and arbitrary composition histories makes inclusion largely definitional and supplies no floor.

Moreover, the coefficients together with their discriminant retain \(N\). Unrestricted integer postprocessing can factor and compute any nonquadratic target. A meaningful barrier needs a statistic grammar, path-length budget, specification of access to full coefficients, and restriction on postprocessing.

### 4.3 Status of the open class

No new nonquadratic leakage mechanism for squarefree balanced \(pq\) is proved. Equally, there is no theorem here exhausting nonmultiplicative coefficient statistics by genus characters, norm smoothness/squareness, collisions, period/order information, or approximate reconstruction.

A scoped **conjecture** could assert no residual signal for a particular preregistered polynomial-length family after removing certified splits and known genus witnesses. It would concern that family, not all coefficient arithmetic. A leading coefficient with a proper gcd, a useful factor approximation, or a recovered period/midpoint is a legitimate route to non-genus target information, but already invokes an established reconstruction mechanism. Additional claims should identify the numerical information and explain why those mechanisms do not account for it.

## Q5. Proposed battery and one-sided positive controls

### 5.1 Population and baseline

Use **200,000 semiprimes total**, equally divided among product bit lengths \(n=48,64,80,104,128\). For \(b=n/2\), sample distinct primes independently from

\[
L=3\cdot2^{b-2}\le p,q<U=2^b,
\]

sort them, and reject squares and duplicate products. Then products have exactly \(n\) bits and \(q<2p\) automatically. Record generation choices. For \(\epsilon\), retain \(N\equiv1\pmod4\); this approximately halves its population. If the intended budget is 200,000 **per size**, the same design scales, with different power.

Use floor covariates \(n\), \(N\bmod256\), and normalized product size with a frozen smooth baseline and an explicit correction for factor top-bit support. A useful continuous prime-density approximation for the archimedean floor is

\[
p_{\min}(N)=\max(L,N/U),\quad p_{\max}(N)=\sqrt N,
\]

\[
b_{\rm arch}(N)=
\frac{\displaystyle\int_{\max(p_{\min},\sqrt{7N/8})}^{\sqrt N}
\frac{du}{u\log u\log(N/u)}}
{\displaystyle\int_{p_{\min}}^{\sqrt N}
\frac{du}{u\log u\log(N/u)}}. \tag{49}
\]

The numerator is zero if its lower limit exceeds its upper limit. The factor \(1/u\) is the Jacobian for conditioning on a product. Formula (49) is a **heuristic smooth prime-density floor**, not an exact finite-sample or residue-conditioned theorem. Calibrate it by training-only residual fits in normalized size and \(N\bmod256\). Fit the same public covariates for \(\epsilon\), rather than assuming exact \(1/2\) balance.

Partition prime pools before creating training and confirmation products, so the final confirmation set shares no prime with training. Within a split, account for reuse in either prime role with dyadic/multiway cluster uncertainty. Freeze arguments, features, tests, and multiplicity correction before confirmation. Residualize both labels and features on public covariates with cross-fitting. Size-bin permutations alone can leave smooth size confounding; show robustness to finer bins and to the geometry-aware floor.

With 40,000 products per size and half for confirmation, an illustrative \(z\)-threshold \(4.5\) corresponds to residual correlation about \(0.032\) on 20,000 independent observations, or \(0.045\) on 10,000. Multiplicity, imbalance, rarity, and prime dependence worsen these limits. Nulls cannot rule out smaller advantages.

### 5.2 Argument menu

Compute a preliminary gcd for each argument and record splits before a unit-only analysis.

| Family | Frozen arguments | Purpose |
|---|---|---|
| Square-root neighborhood | \(\lfloor\sqrt N\rfloor+j,\ j=-2,-1,0,1,2\); next-square gap \(g=(\lfloor\sqrt N\rfloor+1)^2-N\) | Compare the exact no-gcd theorem at \(j=0\) with nearby arguments and junction data. |
| Cubic-size neighborhood | \(\lfloor N^{2/3}\rfloor+j,\ j=-1,0,1\) | Test rounding sensitivity and first remainders. |
| Fixed irrational slopes | \(\lfloor\omega N\rfloor\), \(\omega=(\sqrt5-1)/2,\sqrt2-1,\sqrt[3]2-1\) | Separate a mostly fixed forward word from variable reverse/junction data; certify floors numerically. |
| Near-rational controls | \(\lfloor N/k\rfloor\pm1,\ k=3,5,7,17\) | Exercise bounded-remainder simplification; primarily public-structure controls. |
| Exponent representatives | \(b^{(N-1)/2}\bmod N,\ b^{N-1}\bmod N,\ b=2,3,5\) | Target CRT/order interactions; record gcd/sign events first. |
| Randomized references | Two fresh uniform units drawn by rejection; one frozen salted-hash argument separately | Exact uniform-argument calibration; the deterministic hash is only a pseudorandom heuristic reference. |

Use the inverse through exact reversal instead of counting the same path as an independent observation. Do not silently replace nonunits by neighboring units: that changes the distribution and can censor informative events.

### 5.3 CF readouts and factor diagnostics

Take \(T=\lfloor N^u\rfloor\), \(u=1/8,1/4,3/8,1/2\), and the one-digit split in §1.9. Freeze:

* Each end: length; counts in digit bins \(1,2,[3,7],[8,31],[32,\infty)\); \(\sum\min(b_i,16)\), \(\sum\min(b_i,64)\), \(\sum\lfloor\log_2b_i\rfloor\).
* Finite-level costs: signed digit sums and a small fixed list of occupancies \(\sum1_{b_i\in I}1_{(v_{i-1},v_i)\bmod m\in C}\), \(m=3,5,8\), with frozen bins \(I\) and classes \(C\). Test selected total-cost residues modulo \(3,5,7\).
* Junction: central digit, bit length, residue modulo \(3,5,8\), adjacent digit bins; a few frozen products of forward/reverse costs; gap length between the \(N^{3/8}\)-prefixes.
* Exact Dedekind/Rademacher data: alternating digit sum, its residues, and signed/absolute normalized values. Record the public endpoint correction separately.

Truncation, maxima, position-selected junctions, and products of costs are **derived path statistics**; not all are literally one fixed-\(F\) additive class E. Label them separately.

Record \(\gcd(v_i,N)\), \(\gcd(R_i,N)\) as public exposure diagnostics, outside the fixed-modulus additive class. On random references, compare their rate with the known-factor oracle probability (14). For the unit square-root argument it must be zero. Hidden factors may measure \(\min_i|v_i-p|\) on the \(N^{1/4}\) scale during analysis, but this is an **oracle diagnostic**, never an input feature.

Report overall effects and effects after excluding certified splits and known genus witnesses. Show exclusion prevalence and covariate shifts; exclusion changes the population. If an effect disappears under exclusion, an existing mechanism accounts for it at the tested resolution.

### 5.4 Two-/three-factor cone features

For a small frozen selection of arguments, compute

\[
F_2(n;1,a;N),\quad F_3(n;1,a,a+1;N),
\quad n\in\{0,\lfloor N/2\rfloor\},
\]

subject to unit checks. Reserve \(F_3(n;1,a,a^2+1;N)\) for a separately budgeted secondary family. Validate the exact evaluator on a small label-blind pilot and measure cost before committing to all instances: polynomial complexity does not guarantee a fast implementation.

For each tuple remove the explicit aggregate:

\[
J_2=\frac{12NF_2-\operatorname{rep}_N(C_2)}N,\qquad
J_3=\frac{24NF_3-\operatorname{rep}_N(C_3)}N.
\]

These are exact integers. Test sign, bit length, residues modulo \(3,5,7,8\), and a frozen robust normalization, with the same floor residualization. Include the zero identity (42) as an arithmetic control. Do not put \(F_d(p)+F_d(q)\) into the public feature table: by (41) it is a factor-oracle control.

### 5.5 Reduced-form readouts: a separate exploratory family

For \(k=1,3,5\), inspect at most 128 and 512 steps of the principal \(\sqrt{kN}\) path. Freeze counts/clipped sums of \(Q_i/\sqrt{kN}\), residues of \(P_i,Q_i\) modulo \(3,5,8\), parity-separated digit bins, and a few successive-coefficient ratio bins. Keep this family separate from rational class E. Record square norms, any actually collected smooth relations, and gcd exposures as explanatory controls. A finite null does not prove a genus-only classification.

### 5.6 One-sided controls with public certificates

1. **\(\epsilon=+1\).** A unit \(z\) with \(z^2\equiv-1\pmod N\) proves both primes are \(1\bmod4\). An odd-index principal norm \(Q_i=r^2\), \((r,N)=1\), gives \(z=A_{i-1}r^{-1}\bmod N\) by (48). If bounded prefixes yield too few witnesses, use hidden factors to generate a valid witness on a preregistered random fraction of positive-label cases, solely for pipeline calibration. Verify it publicly. Failure supplies no negative label.
2. **Archimedean positive.** A Fermat witness \(h^2-N=d^2\), \(d>0\), with \(h\ge15d\), certifies \(p/q=(h-d)/(h+d)\ge7/8\). Equivalently \(N\ge224d^2\). Verify the nontrivial factors \(h\pm d\). A short public Fermat scan gives a natural rare-event control; oracle-generated positive witnesses calibrate power. A failed scan proves nothing about the target.
3. **Local-part oracle.** Compute \(L_2\) or \(L_3\) with hidden factors for nondegenerate tuples and run (41). This validates exact rational arithmetic and the decoder and recovers both labels. Gate it to a random subset of positive cases if a sparse one-sided detection control is wanted. It is unavailable auxiliary information.
4. **Floor/type-I controls.** Generate synthetic labels from the fitted public floor, or synthetic features depending only on those covariates. Residual tests should find no beyond-floor effect. Blind label permutations separately calibrate multiplicity.

Store oracle controls separately from public features. Their purpose is to verify detection of a real certificate at stated prevalence, not to count known-factor constructions as new leaks.

## Ranked mechanisms worth further work

1. **Partial local-value separation in two-/three-factor sums — algebraically strongest, algorithmically unresolved.** Equations (32)–(41) specify a decisive residue exposing \(p+q\). Search for coarse public readouts correlated with it after removing the aggregate. Exact separation already is factoring; the decoder is an obstruction and a precise target, not a bypass.
2. **Junction/reverse data for structured exponent residues — the most direct class-E frontier.** Their CRT components are not mere signs, and uniform-numerator theorems do not apply. A replicated target-conditional deviation may be real. First exclude gcd/order events and floor errors. Ordinary endpoint equidistribution cannot dismiss fine junction costs.
3. **Three-factor weighted cone costs — a computable extension beyond narrow fixed-state class E.** Formula (44) supplies concrete readouts with public modular aggregates and no factor oracle. A small validated experiment is justified. Neither a factor-separation mechanism nor reduction to ordinary fixed-state CF costs is proved.
4. **Full reduced-form coefficient statistics after known certificates are removed — broad and weakly constrained.** A finite-state classification does not follow from the reduction-word correspondence. Narrowly specified families may reveal an unmodeled bias, but no new nonquadratic mechanism or exclusion theorem is established.
5. **Direct factor-bearing rational convergents — low priority as a new leak.** Their uniform probability is \(O(N^{-1/2})\), and they are absent for the unit square-root argument. Forcing them efficiently is already factoring. Approximate variants must achieve the \(N^{1/4}\) error scale to have an established decoder.
