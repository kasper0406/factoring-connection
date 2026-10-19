# ldbound_1: a multiplier design, and the distinction between two residue averages

2026-10-19. All asymptotics are as N tends to infinity. The repository outside this scratch directory was read only.

**PROVED in this report:** a restricted all-residue second-moment upper bound, including non-reduced classes; bounded multiplicity of the admissible fractions in the original filtered design; and LD* for the multiplier design defined below. **CONJECTURED / unresolved here:** LD and LD* for the original fixed-resultant design. The factoring conclusion for the multiplier design still has the repository's LV6 assumption (and the stated smoothness-testing assumption, if that version of the cost is used).

The new design is

\[
 F_{g,z,c}=z\widehat F_g+G_gR_c,\qquad
 1\le z\le Z_0=N^h,\quad |c_i|\le C_d Z_0,
 \tag{M}
\]

with z and c independent and uniform, and g, candidate points, and the filter as in (a''_k). The fixed constant C_d is chosen sufficiently large in the first-moment proof. Thus

\[
 F_{g,z,c}(g_0,g_1)=Nz.
\]

As in the notes, a nontrivial gcd(g_1,N) terminates with a factor. For statements about actual number fields one may equivalently restrict to (g_1,N)=1 and normalize the g law on that subset. Its size is still asymp M^2. The fibre lower bound with this restriction is checked in the FM proof below. Analytic upper bounds are preserved when this restriction is imposed.

Take

\[
 s=1/(d+1),\quad \nu=1/(4k),\quad h>s+\nu.
 \tag{M-range}
\]

The resulting budget is **PROVED arithmetic**:

\[
 \theta_M=2s+h+(d-1)\nu,
 \qquad \inf_{h>s+\nu}\theta_M=3s+d\nu.
 \tag{M-cost}
\]

In particular the infimum over fixed k and positive height margins is unchanged: (1+3/(d+1))/sqrt(2). This is an infimum, not an assertion that a fixed positive height margin is an o(1).

## 0. Conventions and analytic inputs

Put Lambda=sqrt(log N log log N), L=exp(Lambda), and B=exp(a Lambda), with a>0 fixed. All parameters d,k,a and all polynomial margins are fixed. A bound exp(O(sqrt(log N))) or exp(O(sqrt(log N/log log N))) is L^{o(1)}. A fixed power of L is not. Polynomially small additive terms below are compared explicitly with the desired main scale before they are discarded.

Write Psi(x)=Psi(x,B), rho(x)=Psi(x)/x, and S_B for the positive B-smooth integers. Signed smooth integers exclude zero. Let

\[
 A_m(a;x)=\#\{1\le v\le x:v\in S_B,\ v\equiv a\pmod m\}.
\]

The following are the **KNOWN** inputs used in the new analytic proof. The precise primary source is [Harper, arXiv:1208.5992](https://arxiv.org/pdf/1208.5992), with the locations specified here. No GRH or Siegel lower bound is used.

* Smooth Numbers Result 1, p. 6 (Hildebrand--Tenenbaum): for x>=y>=2 the saddle-point expression has relative error O(1/log(u+1)+1/log y). In our range its consequence is
  \[
  x^{\alpha}\zeta(\alpha,B)\ll (\log N)^C\Psi(x),\quad
  \alpha=1-\frac{\log(u\log u)+O(1)}{\log B}\longrightarrow1.
  \tag{SP}
  \]
  We use this only for x=N^{c+o(1)}, with c in a fixed compact subset of (0,infinity). The displayed alpha expansion is valid for log x<B<=x^{1/3}, as recorded on p. 7, and our B satisfies this.
* Proposition 1, p. 9: log^{1.1}x<=y<=x; nonprincipal chi of conductor r<=x^{d_H}, modulus m<=x; zero-free region Re(s)>1-epsilon, |Im(s)|<=H; C_H/log y<epsilon<=alpha(x,y)/2; y^{0.9 epsilon}log^2x<=H<=x^{d_H}; and either y>=(Hr)^{C_H} or epsilon>=40 log log(myH)/log y. Its conclusion is
  \[
  |\Psi(x,y;\chi)|\ll\Psi(x,y)\sqrt{\log x\log y}
       (x^{-0.3\epsilon}\log H+H^{-0.02}).
  \tag{H1}
  \]
* Zeros Result 1(iv), p. 8, and its application on p. 15: the primitive characters with R<r<=2R having a zero in Re(s)>299/300, |Im(s)|<=r^{100} satisfy
  \[
  \sum_{\chi^*\ {m exceptional}}\frac1{\varphi(r)}\ll R^{-1/10}.
  \tag{ZD}
  \]
  This is effective, for R>=3. Consequently their number is O(R^{9/10}).
* Multiplicative Large Sieve 1, p. 18: for any complex sequence on an interval of length X,
  \[
  \sum_{r\le Q}\frac r{\varphi(r)}\sum_{\chi^*\bmod r}
       |\sum a_v\chi^*(v)|^2\le (X+3Q^2)\sum|a_v|^2.
  \tag{LS}
  \]

These statements include imprimitive characters where indicated in H1; LS itself is used only for primitive characters. The constants d_H,C_H are Harper's constants, not the polynomial degree.

The first-moment proof also uses **KNOWN / checked repository input** Lemma 8.5' of pcest_results, §2.2, with x=N^{xi+o(1)}, Q=N^{gamma+o(1)}, 1/(2k)+eta_1<gamma<xi/2. Its proof uses [Lee--Venkatesan, §8, Lemma 8.5 and (8.2)](https://arxiv.org/pdf/1805.08873): (log x)^K<y<x^{1/log log x}, x^epsilon<=Q<=sqrt(Psi(x,y)), omega->infinity and omega=y^{O(1)}. In the fixed-k regime the proof bounds the character majorant W, not merely a residue maximum. This input is used for FM only; it is not squared to prove the theorem below.

Other elementary or repository inputs are proved when they are needed, or identified by their exact lemma. Prime counting on (Y,2Y] uses the classical prime number theorem, pi(2Y)-pi(Y)~Y/log Y, Y->infinity (see [NIST DLMF, §27.12](https://dlmf.nist.gov/27.12)). Nothing requires a prime in a subpolynomial interval. If ell divides N, a factor has already been found; throughout the moment analysis ell does not divide N, as in the notes.

The elementary estimates tau(v)<=C_epsilon v^epsilon and A^{omega(v)}<=C_{A,epsilon}v^epsilon, for each fixed epsilon>0, can be proved by separating the finitely many primes below 2^{1/epsilon}, respectively A^{1/epsilon}. For the other primes use e+1<=2^e<=p^{epsilon e}, respectively A<=p^epsilon; for each of the finitely many remaining primes absorb (e+1)p^{-epsilon e} in its finite supremum. We also use v/phi(v)<<log(2v). To prove this weaker sufficient bound, expand v/phi(v)=sum_{d|v}mu(d)^2/phi(d). The identity
sum_{v<=T}1/phi(v)=sum_{d<=T}mu(d)^2/(d phi(d)) sum_{j<=T/d}1/j<<log(2T)
follows from the convergent Euler product for sum mu(d)^2/(d phi(d)); the stated bound follows in turn from the divisor expansion. These suffice everywhere a logarithmic loss is allowed.

## 1. Q1: the Farey premise needs an admissibility correction

**PROVED.** In the original filtered design, if t=|G(P)|=ell n, n is B-smooth, then

\[
 (b,g_1)=(g_1,t)=1.
 \tag{1.1}
\]

Indeed b is prime, b>B, b<ell, and gcd(a,b)=1. If b divides g_1 then b divides t, which is impossible. Since (g_0,g_1)=1,

\[
 (g_1,t)=(g_1,g_0b)=(g_1,b)=1.
\]

**PROVED (bounded primitive representations).** Fix a filtered t and u mod t. The number of admissible (P,g) with b/g_1=u mod t and G(P)=+t or -t is bounded by an absolute constant.

To see this, first consider (b,g_1). Distinct such primitive pairs give distinct real fractions b/g_1. Any two satisfying the same modular equation have nonzero determinant divisible by t, and hence their real fractions differ by at least t/M^2. All those fractions lie in [V/(2M),2V/M]. Their number is at most

\[
 1+\frac{3MV}{2t}\le49,
\]

because t>=MV/32. For each pair (b,g_1), the equation g_1a-g_0b=+t, or -t, has solutions (a,g_0) on a translate of Z(b,g_1). The g_0 interval has length M/2<=g_1, so there are at most two possibilities for each sign. This proves the assertion with, for example, 196.

Thus **FALSE for the actual admissible fractions** is the assertion that arbitrarily large K-fold Farey multiplicities remain in this filtered, prime-b design. The argument also explains why removing the filter or the prime-b/reduced-fraction condition changes the issue.

The passage u -> N u^d may still identify several u's. Its multiplicity is at most

\[
 J_d(t):=\#\{u\in(\mathbb Z/t\mathbb Z)^*:u^d=1\}\le 2d^{\omega(t)}.
 \tag{1.2}
\]

For odd prime powers the unit group is cyclic and the kernel has size at most d; for powers of 2 its cyclic factors give at most 2d. CRT proves (1.2). A pointwise assertion J_d(t)=L^{o(1)} for all polynomial-size smooth t is not made.

**PROVED averaged replacement.** For fixed A>=1 and R=N^{r+o(1)},

\[
 \sum_{n\le2R,\ n\in S_B} A^{\omega(n)}
 \le \Psi(R)\exp(O_A(u_R+\log\log B))(\log N)^C
 =\Psi(R)L^{o(1)}.
 \tag{1.3}
\]

Proof: use Rankin with alpha_R in (SP). The weighted Euler product divided by zeta(alpha_R,B) is

\[
 \prod_{p\le B}(1+(A-1)p^{-\alpha_R}).
\]

Chebyshev's prime bound and partial summation give

\[
 \sum_{p\le B}p^{-\alpha_R}
 \ll \log\log B+
 \frac{B^{1-\alpha_R}}{(1-\alpha_R)\log B}
 \ll \log\log B+u_R.
\]

For the last inequality use B^{1-alpha_R} asymp u_R log u_R and (1-alpha_R)log B asymp log u_R from (SP). The finitely many small u_R cases are irrelevant here. This proves (1.3).

The weighted version of Theorem AR below, proved at the end of §3, makes this replacement valid even when lift energy is attached to n. Therefore the **unweighted marginal sampling** of the original structured residues causes no L-scale loss. This does not include the partner's rational-smoothness indicator; that distinction is essential.

**PROVED model barrier without primitivity.** Take t prime, MV=c t with a fixed 0<c<1, and let K->infinity with K log V=o(V) and V<=M. Allow all numerators and denominators in [V/2,V] times [M/2,M]. Choose primitive (b_0,g_0) in [0.55V/K,0.60V/K] times [0.55M/K,0.60M/K]. Mobius inversion gives Theta(MV/K^2) choices: its error is O((M/K)log(V/K)), smaller than the area under the condition on K. Each has Theta(K) scalar multiples in the original box. The determinant of two of the smaller vectors has magnitude <t, so distinct primitive vectors give distinct modular slopes. Moreover the determinant of a smaller vector and any vector in the original box has magnitude O(MV/K)<t. Thus another vector in that box with its modular slope must be a scalar multiple; the representation count is Theta(K), not merely bounded below by K. The counting measure of high-multiplicity residues is Theta(K^{-2}), while their representation-weighted measure is Theta(K^{-1}). An energy of size K^2 on them has bounded all-residue mean and representation-weighted mean Theta(K). Taking K=L^epsilon with epsilon>0 small makes this an impermissible loss and still satisfies the geometric size conditions. This is a countermodel to an inference from box area and all-residue L^2 alone. It is not a counterexample to smooth-number dispersion and it does not satisfy (1.1).

## 2. Q2: what the original averaging does and does not remove

**PROVED identity.** Fix P=(a,b), P'=(a',b'), Delta=epsilon ell q with epsilon=+1 or -1, and G(P)=ell n. Integer points of the g fibre have the form

\[
 (g_0,g_1)=(g_0^*,g_1^*)+j(a,b).
\]

Consequently

\[
 n'(j)=n'_0-\epsilon qj,
 \qquad
 r_0(j)=N\left(\frac b{g_1^*+bj}\right)^d\pmod{\ell n},
 \tag{2.1}
\]

and, equivalently,

\[
 bn'=b'n-\epsilon qg_1.
 \tag{2.2}
\]

The j interval has length O(M/V). Thus the partner indicator 1_B(n'(j)) and the coarse residue use the same j. The coefficient randomness c changes lift values but does not change r_0(j). Averaging c therefore does not create an independent coarse-residue variable.

**PROVED location of the loss.** After a marginal all-residue bound one can discard 1_B(n') by positivity. This replaces its expected density rho_r by 1 and loses

\[
 \rho_r^{-1}=L^{(s-\nu)/(2a)+o(1)}.
 \tag{2.3}
\]

Alternatively, keeping the partner density by a fibre divisor-switching bound makes one take the supremum of the attached lift energy over that fibre. Controlling the unweighted average of that energy does not bound this supremum or its correlation with 1_B(n'). In the max-coarse-class version this can cost up to the coarse modulus t; it cannot be called L^{o(1)}.

**PROVED information-theoretic countermodel.** On each fibre let a binary partner weight have density p=rho_r. Give its support lift energy 1/p and give the complement energy 1. The unweighted mean energy is below 2; the weighted mean is 1, whereas the desired independent scale is p. The loss is 1/p. This can be repeated at every outer index ell,q,n. Even exact uniformity of coarse-class totals is compatible with the model: inside the selected coarse classes put their fixed total mass on only a p fraction of the fine lifts. Since rho_F<p and |Delta| is polynomial, the model respects the per-lift capacity after rounding. It is a model of what marginal estimates cannot imply, not an arithmetic disproof of LD*.

**CONJECTURED for (a''_k).** The actual mixed correlation in (2.1), with the smooth-number lift energies, satisfies LD*. Neither this model nor the Farey calculation settles that assertion. In particular, outer averaging is not by itself a proof of independence.

## 3. Q3: a restricted all-residue upper bound

### Theorem AR [PROVED]

Let

\[
 x=N^{\xi+o(1)},\quad R=N^{r+o(1)},\quad Y=N^{\lambda+o(1)},
 \quad r,\lambda>0,\quad \xi>r+3\lambda,
 \tag{3.1}
\]

with fixed positive margins. Let 1<=Z<=L^C, C fixed, and set m=ell^2 n q, where ell is prime in (Y,2Y], R<n<=2R is B-smooth, and 1<=q<=Z. This is a multiset indexed by (ell,n,q). Then

\[
 \sum_\ell\sum_{n\sim R,\ n\in S_B}\sum_{q\le Z}
       \sum_{a\bmod m} A_m(a;x)^2
 \le L^{o(1)}\frac{\Psi(x)^2\rho(R)}{Y\log Y}.
 \tag{AR}
\]

The omitted factor log(2Z) on the right is L^{o(1)}. It can be included explicitly throughout the proof. The same assertion holds with any window contained in [-x,x], excluding zero. The centered variance also satisfies (AR) when, for e=(a,m), the subtracted term is

\[
 M_m(a)=
 \begin{cases}
 \Psi_{m/e}(x/e,B)/\varphi(m/e),&e\in S_B,\\
 0,&e\notin S_B.
 \end{cases}
 \tag{3.2}
\]

This is a main-square upper bound within L^{o(1)}. It is not a vanishing relative-variance theorem. It retains the smooth-n sparsity factor rho(R), as well as the polynomial sparsity from ell^2.

### 3.1 A large-conductor lemma [PROVED from LS]

For 1<=U<=Q<=x,

\[
 \mathcal H(Q,U):=
 \sum_{v\le Q}\frac1{\varphi(v)}
 \sum_{\substack{\chi\bmod v\\\operatorname{cond}\chi\ge U}}
 |\Psi(x,B;\chi)|^2
 \ll\Psi(x)(Q+x/U)(\log(2x))^C.
 \tag{3.3}
\]

This includes imprimitive characters. For completeness, if chi is induced from chi* mod r,

\[
 \Psi(x,B;\chi)=
 \sum_{\substack{d\mid v/r\\d\in S_B}}
 \mu(d)\chi^*(d)\Psi(x/d,B;\chi^*).
\]

Cauchy gives the factor tau(v/r). For r in a dyadic interval [D,2D], summing over v with rd|v uses

\[
 \sum_{v\le Q,\ rd\mid v}
 \frac{\tau(v/r)}{\varphi(v)}
 \ll\frac{\tau(d)}{\varphi(r)\varphi(d)}
      (\log(2Q/(rd)))^2.
\]

Indeed tau(ds)<=tau(d)tau(s), phi(rs)>=phi(r)phi(s), and sum_{s<=T}tau(s)/phi(s)<<log^2(2T), by expanding tau and using sum_{s<=T}1/phi(s)<<log(2T). The latter follows from n/phi(n)=sum_{h|n}mu(h)^2/phi(h), with the resulting convergent reciprocal sum.

Apply LS to the primitive characters r~D and the sequence 1_{S_B} on [1,x/d]. The resulting bound is

\[
 \ll\Psi(x)\sum_{d\le Q/D}\frac{\tau(d)}{\varphi(d)}
       (x/(dD)+D)(\log(2x))^2
 \ll\Psi(x)(x/D+D)(\log(2x))^4.
\]

Here sum tau(d)/(d phi(d)) converges and sum_{d<=T}tau(d)/phi(d)<<log^2(2T). Summing the dyadic D proves (3.3). No unreduced additive frequency has been put into a primitive large sieve.

### 3.2 Reduced residues: principal ell-character

By character orthogonality, the uncentered reduced sum is

\[
 \sum_{(a,m)=1} A_m(a;x)^2
   =\frac1{\varphi(m)}\sum_{\chi\bmod m}|\Psi(x,B;\chi)|^2.
 \tag{3.4}
\]

For large N, ell>B and ell>Z, hence (ell,nq)=1. A character mod ell^2 nq factors uniquely. If its ell^2-component is principal, its smooth-number sum is exactly the sum of the corresponding character mod nq. The prefactor is 1/phi(ell^2).

Choose a fixed eta>0 sufficiently small that

\[
 100\eta<d_H,\quad \eta<1/4000,
 \quad \eta\xi<\min(r,\lambda)/10,
\]

and set m_0=exp(eta sqrt(log x)), U=x^eta. The constants can be reduced once more without affecting the argument. All ensuing choices are uniform on compact ranges with the margins in (3.1).

**Small conductors, including 1.** There are O(m_0^2) primitive characters of conductor <=m_0. Bound their smooth sums by Psi(x). Since

\[
 \sum_{n\sim R,\ n\in S_B}\sum_{q\le Z}\frac1{\varphi(nq)}
 \ll \rho(R)(\log N)^C\log(2Z),
\]

their contribution after summing ell is

\[
 \ll\frac{\Psi(x)^2\rho(R)}{Y\log Y}
       m_0^2(\log N)^C\log(2Z).
 \tag{3.5}
\]

This is admissible because log(m_0^2)=O(sqrt(log N))=o(Lambda). This deliberately does not attempt to prove an o(1) variance. A possible Siegel character is included without an ineffective constant.

**Generic medium conductors, m_0<r<=U.** Call a primitive character generic if it has no zero in the region of (ZD). Apply H1 with

\[
 \epsilon=\min(1/300,10\log r/\log B),\quad H=r^{100}.
\]

The following checks are explicit. B>=log^{1.1}x and B<=x. Moduli nq<=2RZ<x. The conductor and H are <=x^{d_H}; epsilon<=alpha/2 eventually and epsilon>C_H/log B. Also H>=B^{0.9 epsilon}log^2x. The second alternative epsilon>=40 log log(nq B H)/log B holds because log r>=eta sqrt(log x)>>log log x. If epsilon is uncapped, x^{-0.3 epsilon}=r^{-3u_x}; if capped, x^{-0.3 epsilon}=x^{-0.001}. Our eta makes its logarithmic multiplier smaller than r^{-2} for r<=U. Thus

\[
 |\Psi(x,B;\chi)|\ll\Psi(x)(\log N)^C r^{-2}.
 \tag{3.6}
\]

The same bound holds for every induced character to nq: Euler factors added to L have zeros on Re(s)=0, and H1 is stated for imprimitive characters. For each modulus, sum phi(r)r^{-4} over r>m_0; it is bounded by an absolute constant. Summing its 1/phi(nq) prefactor proves an acceptable bound of the scale in (3.5), even without the m_0^2 factor.

**Exceptional medium conductors.** Put alpha_R=alpha(R,B). For a fixed primitive conductor r, define

\[
 W_r=\sum_{q\le Z}\sum_{\substack{n\sim R,\ n\in S_B\\r\mid nq}}
      \frac1{\varphi(nq)}.
\]

With h=r/(r,q), the n count is zero unless h is B-smooth, and otherwise is at most Psi(2R/h). Rankin and (SP) give, uniformly in h,

\[
 \Psi(2R/h)\ll(\log N)^C\Psi(R)h^{-\alpha_R}.
\]

Therefore

\[
 W_r\ll\rho(R)(\log N)^C r^{-\alpha_R}
       \sum_{q\le Z}\frac{(r,q)^{\alpha_R}}q
 \ll\rho(R)(\log N)^C\log(2Z)\tau(r)r^{-\alpha_R}.
 \tag{3.7}
\]

The last step expands over divisors h|(r,q), and uses sum_{q<=Z,h|q}1/q<=log(2Z)/h and alpha_R<=1. By (ZD), the number of exceptional characters with r~D is O(D^{0.9}). For large N, alpha_R>0.99 and tau(r)<=r^{0.01} for r>=m_0. Hence their total, using the trivial bound Psi(x), is

\[
 \ll\Psi(x)^2\rho(R)(\log N)^C\log(2Z)
       \sum_{D\ge m_0\ {m dyadic}}D^{-0.08}.
 \tag{3.8}
\]

After the ell prefactor this is smaller than the required scale. This is exactly where the smooth-n restriction is retained. Treating its density as 1 would fail.

**Large conductors, r>=U.** For fixed v=nq, its number of representations is <=tau(v). For every fixed epsilon_0>0 it is <=N^{epsilon_0} for large N. We use this only in a term with a fixed polynomial saving, not as an L^{o(1)} factor. By (3.3) the contribution is

\[
 \ll \frac{N^{\epsilon_0}}{Y\log Y}\Psi(x)
        (RZ+x/U)(\log N)^C.
 \tag{3.9}
\]

Relative to the target, its two terms are N^{r-xi+epsilon_0+o(1)} and N^{-eta xi+epsilon_0+o(1)}. Choose epsilon_0 smaller than half of both fixed gaps. They are polynomially small.

### 3.3 Nonprincipal ell-character: the essential ell^2 split

There are two cases. A character mod ell^2 has conductor ell or ell^2. On B-smooth integers the former is the same character modulo ell. Its prefactor satisfies

\[
 1/\varphi(\ell^2 nq)=\ell^{-1}/\varphi(\ell nq).
\]

Every reduced modulus ell nq has exactly one prime factor in (Y,2Y] and the pair (n,q) has at most tau(nq)<=N^{epsilon_0} representations. Apply (3.3) with Q=O(YRZ), U=Y, and the extra factor 1/Y. This gives

\[
 \ll N^{\epsilon_0}\Psi(x)(RZ+x/Y^2)(\log N)^C.
 \tag{3.10}
\]

For conductor divisible by ell^2, apply (3.3) directly with Q=O(Y^2RZ), U=Y^2:

\[
 \ll N^{\epsilon_0}\Psi(x)(Y^2RZ+x/Y^2)(\log N)^C.
 \tag{3.11}
\]

Divide (3.10)--(3.11) by Psi(x)^2 rho(R)/(Y log Y). Their polynomial exponents are respectively at most

\[
 r+\lambda-\xi+\epsilon_0,
 \quad-\lambda+\epsilon_0,
 \quad r+3\lambda-\xi+\epsilon_0.
\]

All are negative by (3.1), with epsilon_0 chosen small enough. The factors Psi(x)/x, rho(R), Z, and logarithms are N^{o(1)} here and are absorbed against those fixed polynomial gaps. They have not been inserted as multiplicative L^{o(1)} losses into a main term.

This completes the reduced-residue proof. In particular, the loss from forgetting the sparse family in the large sieve is paid with polynomial room, after distinguishing conductor ell from conductor ell^2. Treating all nonprincipal ell-components as merely conductor >=ell would erase that room.

### 3.4 Non-reduced classes

For e=(a,m), the count is zero if e is not B-smooth; otherwise it is the reduced-class count at (x/e,m/e). Since ell>B, such e has no ell factor. Choose a small fixed delta>0, less than r/4 and xi/4.

For e>N^delta, every original class has at most x/m+1 integers. Therefore its entire gcd layer has square sum at most (x/m+1)Psi(x/e). By Rankin,

\[
 \sum_{\substack{e\mid m\\e>N^\delta}}\Psi(x/e)
 \ll\Psi(x)(\log N)^C
       N^{-\delta\alpha_x/2}\sum_{e\mid m}e^{-\alpha_x/2}
 \le\Psi(x)N^{-\delta/3+o(1)}.
\]

Because m<x with a fixed polynomial margin and x/Psi(x)=N^{o(1)}, summing these layers is polynomially small relative to the right side of (AR).

For e<=N^delta, write u=(e,n), v=e/u. Then u|n, v|q, and

\[
 n=un_1,\quad q=vq_1,\quad
 m/e=\ell^2n_1q_1,\quad x/e=x/(uv).
\]

Drop the additional conditions of this canonical decomposition and sum over all smooth u,v. Apply the reduced theorem with R'=R/u, x'=x/(uv), Z'=Z/v. Its hypotheses hold uniformly: R' has exponent >=3r/4, and

\[
 \xi'-r'-3\lambda=(\xi-r-3\lambda)-\log v/\log N
 =\xi-r-3\lambda+o(1)>0
\]

since v<=Z. Using Rankin twice, the resulting right side relative to the original right side is bounded by a polylogarithm times

\[
 u^{1-\alpha_R-2\alpha_x}v^{-2\alpha_x}.
\]

Both exponents are below -1.9 for large N. Summing u and v costs an absolute constant. This proves (AR) for all residues. The centering (3.2) is the mean within each gcd layer, so its centered square sum is no larger than its uncentered square sum. For signed windows use (A(a)+A(-a))^2<=2A(a)^2+2A(-a)^2. Zero must remain excluded. This proves all assertions of AR.

### 3.5 Ranges, weights, and the precise strength

**PROVED application to the original parameters.** There r=s-nu, lambda=2nu, xi=s+h+dnu, so

\[
 \xi-r-3\lambda=h+(d-5)\nu>0
\]

under h>=2nu+eta_0, for every d>=3. The moduli have size N^{s+3nu}L^{O(1)}, exactly the requested exponent. No half-level assertion for these enlarged moduli is necessary.

Define the all-coarse-class averaged lift energy

\[
 \mathcal E_{\ell,n,q}(x)=\frac{\ell^2nq}{\Psi(x)^2}
              \sum_{a\bmod\ell^2nq}A_{\ell^2nq}(a;x)^2.
\]

Since n~R and ell~Y, multiplying (AR) by O(Y^2R)/Psi(x)^2 gives

\[
 \frac1{\pi_J}\sum_\ell\sum_{n\sim R,\ n\in S_B}
       \sum_{q\le Z}\frac1q\mathcal E_{\ell,n,q}(x)
 \le L^{o(1)}\Psi(R).
 \tag{3.12}
\]

This includes the natural harmonic q-weight. It is the all-residue assertion used in §5. A bound uniform in each q is not claimed: (3.7) used the q average. The ell average is used in (3.10)--(3.11); at the full stated range the proof does not supply a uniform fixed-ell assertion. Averaging smooth n is used in (3.7), and cannot be discarded there.

The large-sieve diagonal itself has relative size at most Y^2RZ/Psi(x)=N^{r+2lambda-xi+o(1)}. It is harmless. The stronger r+3lambda<xi is the cost of enlarging to all moduli in the ell^2-conductor contribution, not the actual diagonal of the sparse family.

**PROVED weighted extension.** (AR) and (3.12) remain valid with a factor A^{omega(n)}, fixed A, inside their sums. Here are the extra checks, to avoid multiplying unrelated averages. Replace the Rankin product in (SP) for the n count by the product in (1.3), at cost E_A=exp(O_A(u_R+log log B))(log N)^C=L^{o(1)}. When h|n, submultiplicativity gives A^{omega(n)}<=A^{omega(h)}A^{omega(n/h)}, so (3.7) gains E_A A^{omega(h)}. The latter is <=C_epsilon h^epsilon for every fixed epsilon>0. With epsilon=0.01 and tau(r)<=r^{0.01}, the exceptional dyadic exponent is still < -0.05. In the large-conductor terms the maximum weight is N^{epsilon_0} with epsilon_0 chosen smaller, and is absorbed by the same polynomial gaps. In §3.4 the small-layer u sum gains A^{omega(u)} and still converges, since its exponent is <-1.9; the large-layer weight can again be absorbed in a polynomial saving. Small and generic conductor contributions gain E_A. This proves the extension with a uniform L^{o(1)} bound.

Together with §1, this proves the unweighted marginal structured-residue assertion at the relevant averaged normalization: if K_{t,D}(r)=Dt^2 Psi(x)^{-2}sum_j A_{tD}(r+tj;x)^2, then

\[
 \frac1{\pi_J}\sum_\ell\sum_{n\sim R,\ n\in S_B}\sum_{q\le Z}\frac1q
 \frac1t\sum_{(P,g):\,|G(P)|=t}K_{t,\ell q}(N(b/g_1)^d)
 \le L^{o(1)}\Psi(R).
 \tag{3.13}
\]

Indeed each coarse residue has at most 196 J_d(t) preimages, and (1/t)sum_r K_{t,D}(r)=mathcal E_{ell,n,q}(x). The weighted AR extension, rather than an independent use of (1.3), pays for J_d(t). Partner rational weights are absent from (3.13).

## 4. Q4: what changes can remove the missing average?

**PROVED obstruction to a coefficient-only repair at fixed resultant.** If two integral homogeneous forms of degree d have the same value N at the primitive point (g_0,g_1), their difference is divisible in Z[x,y] by G=g_1x-g_0y. Over Q this follows by substituting x=g_0y/g_1; integrality of the quotient follows from Gauss's lemma because G is primitive. Thus the existing R family already spans every such coefficient perturbation. Its coarse congruence remains N(b/g_1)^d, regardless of the size or law of R.

This does not prove that every alternative rational map fails. It identifies exactly why more coefficients in the same fixed-resultant affine lattice cannot supply a new independent coarse residue.

**PROVED valid multiplier design.** Use (M). The identity of Lemma 1.2 becomes

\[
 g_1^d F(P)\equiv Nz b^d\pmod{G(P)}.
 \tag{4.1}
\]

The common-root map mod N remains well-defined: F(g_0,g_1)=Nz=0 mod N. Its relation image remains G(P)/g_1. There is no z factor in that image and no new z parity column. The ordinary leading-coefficient bookkeeping and the g_1 column remain as in the repository's LV6 problem.

The rational law, candidate lattice, and filter have not changed, so Lemma Q, Lemma DS, Lemma 1.3, and Theorem R apply exactly. Conditional on z the random perturbation is still G R_c, with the same v(P), image lattice Gamma, and kernel. Lemma Kern and Prop. 4.2 therefore apply with coefficient height H_R=C_d Z_0 and offset z Fhat. The requirement H_R>=V^2 follows from h>s+nu>2nu. Norm bounds are

\[
 |G/\ell|\le O(MV/\ell),\qquad |F(P)|\le O_d(MZ_0V^d),
 \tag{4.2}
\]

so the exponents are r=s-nu and xi=s+h+dnu. No factor involving z is omitted from xi.

### First moment for (M) [PROVED]

An elementary scalar-kernel lower bound is useful here and avoids treating a large new offset as a small error.

Let P be primitive, p=max(|a|,b) asymp V, and let c be uniform on [-H_R,H_R]^d, H_R/V->infinity polynomially. There are constants c_d,c'_d>0 such that

\[
 \Pr[\langle c,v(P)\rangle=m]\ge
       \frac{c'_d}{H_R V^{d-1}}
 \quad (|m|\le c_d H_R V^{d-1}).
 \tag{4.3}
\]

Proof: exchange coordinates and signs so the larger coordinate is p and the other is q, 0<=q<=p, gcd(p,q)=1. Successively solve m modulo p using the last monomial q^{d-1}, divide the remainder by p, and repeat. This gives coefficients c_i, i>0, with |c_i|<=p/2 and

\[
 |c_0|\le |m|/p^{d-1}+(d-1)p/2.
\]

For sufficiently small c_d and H_R>>_d V, this is a representative inside [-H_R/4,H_R/4]^d. The scalar kernel consists of multiples of bx-ay, with the d-1 linearly independent monomial multiples as an integral basis, each of height <=V. Take the coefficients in that basis in [-W_0,W_0], W_0=floor(H_R/(8dV)). All these distinct kernel vectors have sup norm <=H_R/4. There are at least c_d H_R^{d-1}/V^{d-1} of them. Add them to the representative and divide by (2H_R+1)^d. This proves (4.3).

For eventual number-field use one may also require the leading coefficient of F to be a unit modulo N. This retains (4.3), with a logarithmic loss, whenever the prescribed F-value is B-smooth and nonzero. In the original coordinates only the first kernel-basis parameter changes the leading coefficient, with step g_1 b. For each prime p|N, p!=b, it excludes one class of that parameter modulo p. If b|N, smoothness of F(P) and (a,b)=1 already force its leading coefficient to be nonzero modulo b, since F(P) is that leading coefficient times a^d modulo b. Mobius inversion on the first basis parameter leaves at least W_0 phi(N)/N-2^{omega(N)} choices, up to an absolute constant for interval length. This is >>W_0/log N, since W_0 grows polynomially and 2^{omega(N)}=L^{o(1)} after preprocessing. The other parameters remain free. Thus the FM argument below can count only fields with unit leading coefficient. Their unconditional probability is 1-o(1): for fixed g,z the leading coefficient is z a_hat_d+g_1 c_0, and the union bound for its prime divisors in N is O(sum_{p|N}1/p+omega(N)/H_R)=o(1).

On the filter t=|G|>=MV/32,

\[
 |z\widehat F(P)|/t\le C'_d Z_0 V^{d-1}.
\]

Choose the fixed ratio H_R/Z_0=C_d sufficiently large. Then this offset is at most one eighth of the central m-window in (4.3). Regardless of the sign of G, the positive F-interval [tW/4,tW/2], W=c_d H_R V^{d-1}, is covered by that central m-window. All these F-heights are within fixed factors of X_F=O_d(MH_R V_max^d).

Call t good if its reduced-class smooth counts on a common finite multiplicative grid covering those fixed factors have relative error at most a sufficiently small fixed constant. A grid with ratio, for example, 1.01 contains a subinterval inside [tW/4,tW/2]. Hildebrand's estimate, in the precise form quoted as LV Fact 3.13, applies here: log B>=(log log x)^{5/3+epsilon} and the interval parameter 1<=z_int<=B^{5/12}; take a fixed z_int=100. Thus the smooth count of that subinterval is comparable to its length times rho_F. For the reduced-class principal term we use the coprimality comparison quoted as LV Corollary 3.17 / pcest §2.1(T): its condition omega(t)<=B^c/log(1+u) holds uniformly since omega(t)<=log t=O(log N), and its relative error is o(1). Choose the goodness constant below the fixed interval comparison constant. The grid is common to ell,P,g,z. These facts and their applications are given in [LV §§3 and 8](https://arxiv.org/pdf/1805.08873).

For good t and (z,t)=1, (4.1), the unit property of Fhat(P), and (4.3) give

\[
 \Pr_c[F(P)\in S_B]\ge \rho_F/(\log N)^C.
\]

Among 1<=z<=Z_0, the fraction coprime to t is at least

\[
 \varphi(t)/t-t/Z_0\gg1/\log N,
\]

because Z_0/t=N^{h-s-nu+o(1)} has a fixed positive exponent. Lemma 8.5' makes bad t rare under the exact fibres exactly as in pcest FM Steps 1 and 4; these fibres and their coprimality counts did not change. Its checks are

\[
 \gamma=s+\nu>2\nu+\eta_1,
 \quad \gamma<\xi/2
 \ \Longleftrightarrow\ h>s-(d-2)\nu,
\]

both with fixed margins. Its y,omega and large-conductor ranges are consequently those already checked in pcest §2.2. The optional restriction (g_1,N)=1 does not destroy the exact fibre lower bound. On g=g*+j(a,b), primitivity removes one j-class modulo each prime dividing t. Requiring (g_1,N)=1 removes one j-class modulo each prime p|N with p!=b; if b|N, the fixed value of g_1 modulo b is already nonzero because b does not divide t. Since (t,N)=1 these restrictions combine by CRT. The resulting lower bound on a central fibre of length asymp M/V is its length times (phi(t)/t)(phi(N)/N), up to a harmless factor for b|N, with error at most 2^{omega(t)+omega(N)}. The latter is N^{o(1)} and is absorbed additively by M/V=N^{s-nu+o(1)}. The density factors cost only logarithms. Also the number of excluded g in the parameter square is o(M^2): the union bound is O(M^2 sum_{p|N}1/p+M omega(N)), and N has no prime factor <=B. Removing bad t costs o(1) of the rational smooth mass, while good z costs only a logarithmic factor. Summing the fibre counts yields

\[
 \mathbb E N_\ell\ge S^2\rho_r\rho_F/(\log N)^C.
 \tag{FM-M}
\]

This is a first-moment proof for the new offset, not an assertion that the old small-offset step applies unchanged.

### Irreducibility and degree [PROVED elementary estimate]

These are not needed for the analytic LD* inequality, but are useful for checking the field design. Uniformly in fixed g,z, the probability that F has degree less than d or is reducible over Q is at most H_R^{-1/4+o(1)}.

Here is a finite-field sieve proof. Take primes p<=Q=H_R^{1/4}, excluding p dividing g_1Nz and a fixed finite set depending on d. For each remaining p, varying c mod p gives the uniform affine hyperplane of degree-at-most-d polynomials with prescribed nonzero value Nz/g_1^d at g_0/g_1. The fraction irreducible of degree d is >=c_d>0. To verify this fraction, fix a nonzero leading coefficient, translate the prescribed point to 0, and count monic irreducibles with a prescribed nonzero constant term. The norm map from F_{p^d}^* to F_p^* has (p^d-1)/(p-1) elements in each fibre. At most O_d(p^{d/2}) elements lie in proper subfields. Divide the remaining full-degree elements by d. This gives asymptotically p^{d-1}/d polynomials for each constant term, uniformly in that term. Summing over the p-1 leading coefficients proves the fraction assertion for d>=3.

There are K asymp Q/log Q usable primes: only O(log N) were omitted. The c cube is uniform modulo p, or modulo pp', up to total variation O_d(p/H_R), respectively O_d(pp'/H_R). Thus the sum of the finite-field irreducibility indicators has mean >=c_d K and variance O_d(K+K^2 Q^2/H_R). Chebyshev gives probability of no irreducible reduction O_d(1/K+Q^2/H_R)=H_R^{-1/4+o(1)}. A reducible rational polynomial, or one of smaller degree, has no degree-d irreducible reduction. This proves the estimate.

Rejecting these fields removes at most O(S^2)H_R^{-1/4+o(1)} from the mean, which is negligible against FM-M; the second-moment upper bound is preserved by positivity. The unit-leading-coefficient restriction was handled directly in the FM fibre lower bound above, so its possibly larger exceptional probability is not incorrectly subtracted as a tiny additive error. Conditioning on these acceptances changes upper bounds by only 1+o(1). If a polynomial has integer content U, then U divides Nz. For unit-leading-coefficient fields, (U,N)=1 automatically, and dividing by U preserves the common root. Every counted smooth value has smooth U and remains smooth after division. The usual leading-coefficient norm bookkeeping can therefore use the primitive part. A proof of LV6 for this law is still a separate requirement, not supplied here.

### Constant [PROVED arithmetic; factoring implication conditional on LV6]

For every fixed odd d>=7, choose k large enough and then h=s+nu+eta with eta>0 sufficiently small. Then

\[
 \theta_M=\frac3{d+1}+\frac d{4k}+\eta<\sqrt2-1,
 \qquad c_M=\frac{1+\theta_M}{\sqrt2}<1.
 \tag{4.4}
\]

At the zero-margin boundary the smallest integer k and the infimum constants are:

| d | smallest k allowing a positive sufficiently small margin | boundary c_M at that k | infimum over fixed k, eta |
|---:|---:|---:|---:|
| 7 | 45 | 0.999770 | 0.972272 |
| 9 | 20 | 0.998788 | 0.919239 |
| 11 | 17 | 0.998268 | 0.883883 |
| 19 | 18 | 0.999770 | 0.813173 |

These figures use the repository's shape-test/smoothness-test cost model. With the Pollard--Strassen replacement its optimization is unchanged after substituting theta_M, and its limiting d=9 constant is 0.967000, still conditional on LV6. Across odd d and fixed k the formal infimum is 1/sqrt(2); no growing-d or growing-k uniform theorem is asserted.

## 5. Q5: full LD* proof for (M)

All expectations in this section include the independent z. Let D=|Delta|=ell q, t=|G(P)|, t'=|G(P')|, and T=dH_R V^{d-1}. Define the lift counts and kappa exactly as in pcest §4.2, with Fhat replaced by z Fhat. Both points are filtered. Put

\[
 K_P(z)=\frac{D\sum_{j\bmod D}\Psi_{I_z}(tD,z\widehat F(P)+Gj)^2}
                   {4T^2\rho(Tt)^2},
 \tag{5.1}
\]

and define K_{P'} analogously, with the lift permutation from Lemma Kern. That permutation is a unit permutation: P and P' are primitive, and b,b' are invertible modulo Delta because b,b'>q and b,b'<ell. Cauchy--Schwarz and the arithmetic-geometric mean give, pointwise in z,

\[
 \kappa_{\rm alg}(z)\le\sqrt{K_P(z)K_{P'}(z)}
                      \le(K_P(z)+K_{P'}(z))/2.
 \tag{5.2}
\]

All windows I_z lie inside one common [-x,x], x=C_d M H_R V_max^d, and x asymp Tt asymp Tt' by the filter. Hence rho(Tt),rho(Tt') and rho(x) are comparable within polylogarithms, uniformly in ell and the pair. Dominate the windows by this envelope.

Since (Fhat(P),t)=1, multiplication by Fhat(P) permutes the residues mod t. Each residue is hit by at most Z_0/t+1 values of z. Because Z_0>=t for all filtered t for large N,

\[
 \begin{split}
 \mathbb E_z K_P(z)
 &\le (\log N)^C\frac{D}{t T^2\rho(x)^2}
                 \sum_{a\bmod tD}A^{\pm}_{tD}(a;x)^2\\
 &\le (\log N)^C\frac{tD}{\Psi(x)^2}
                 \sum_{a\bmod tD}A^{\pm}_{tD}(a;x)^2
 =:(\log N)^C E_{\ell,|n|,q}.
 \end{split}
 \tag{5.3}
\]

Here A^pm counts nonzero signed smooth integers in the envelope. Constants such as 4 and the comparability x/(Tt) are absorbed in the fixed constant. Crucially E depends only on ell,|n|,q, not on b,g_1, the position on the g fibre, or the partner.

This averaging did not require z to be independent in the two marginals. The one shared z suffices because (5.2) is applied before its average. It did not condition z on being smooth or coprime; non-reduced classes are the reason §3.4 was necessary.

It remains to retain the partner density. Fix ell,P,P',n, with det(P,P')=epsilon ell q. Along G(P)=ell n, (2.1) makes n' a progression of step q. Enlarge the g-parameter interval to length O(M/V). Its n' interval has length O(qM/V), bounded below at the same scale by using a fixed enlargement, and at most O(X_r), since q<=2S^2 and X_r asymp MV/ell. The filter restricts |n'| to [X_r/64,X_r] up to fixed constants. Reflect negative parts and cover them by a bounded number of intervals of length comparable to qM/V, within a fixed enlargement of this filtered range. If necessary shorten the constant in the lower endpoint and enlarge the length by a constant; both operations have only constant cost.

Lemma DS applies: q<=L^{O(1)}=B^{O(1)}, and qM/V>=X_r B^{-C} for a sufficiently large fixed C because M/V asymp X_r/S^2. It includes non-reduced residues. We obtain, uniformly in every fixed n,

\[
 \#\{g\text{ on this fibre}:n'\in S_B,\ \text{both filters}\}
       \le (M/V)\rho_r(\log N)^C.
 \tag{5.4}
\]

Primitivity and square-boundary restrictions on g only decrease this count. This is the step that could not be combined with a moving energy in the old design. In (5.3) that energy is now constant along the fibre.

For fixed P,q and sign epsilon there are O(1) possible P': solutions of det(P,P')=epsilon ell q lie on a translate of ZP, whose spacing is at least V/2 inside a box of diameter O(V). There are O(V^2) choices for P and |G| asymp M^2 choices for the field parameter g, where |G| here means the cardinality of the primitive parameter set, not the polynomial value. To avoid that overloading, write this cardinality as |mathcal G| asymp M^2 in the following estimate.

Use (5.2), symmetry of ordered pairs, (5.3), and (5.4). Up to polylogarithms, the entire LD* left side is at most

\[
 S^2\frac{V}{M}\rho_r\,
 \frac1{\pi_J}\sum_\ell\sum_{n\in S_B\cap\text{filtered range}}
                \sum_{q\le2S^2}\frac1q E_{\ell,|n|,q}.
 \tag{5.5}
\]

There are only O(1) dyadic n ranges, with R asymp MV/ell, and their rho(R) are comparable to rho_r within polylogarithms. Apply (3.12), including its signed version, to each. Its hypotheses are

\[
 r=s-\nu>0,\quad\lambda=2\nu>0,\quad
 \xi-r-3\lambda=h+(d-5)\nu>0.
\]

All hold with fixed margins under (M-range), also for d=3 since s>nu. Consequently (5.5) is at most

\[
 L^{o(1)}S^2\frac VM\rho_r\,(MV/\ell)\rho_r
       =L^{o(1)}S^4\rho_r^2.
 \tag{LD*-M}
\]

This is precisely LD* for (M), with its actual rational weights and actual lift normalization. No product of separately averaged rational and algebraic estimates was used. The displayed calculation uses the independent construction law before eligibility rejections. If those rejections are imposed, their joint-law density is bounded by 1+o(1) relative to that law, and the nonnegative LD* sum has the same upper bound. The coefficient-kernel upper bound is likewise dominated by the unconditioned bound; the FM lower bound for eligible fields was checked separately in §4.

At the search size S^2=(log N)^{C_3}/(rho_r rho_F), the same proof and the rational first-moment lower bound imply the relative **filtered** LD formulation. Indeed if R_f is the rational filtered count, E R_f>>S^2 rho_r/polylog tends to infinity, so E[R_f(R_f-1)]>=(E R_f)^2-E R_f>>S^4 rho_r^2/polylog. Divide LD*-M by this lower bound. The original unfiltered LD, with arbitrarily small |G| and its varying algebraic density, is not claimed.

### Consequence [PROVED implication, with stated assumptions]

Choose one common S(N) through S^2=(log N)^{C_3}/(rho_r rho_F), with common dyadic envelope arguments, and choose C_3 above FM-M's logarithmic loss. Then FM-M gives E N>=1. Lemma Kern, Prop. 4.2, the filter comparison, and LD*-M give

\[
 \mathbb E[N(N-1)]\le L^{o(1)}S^4\rho_r^2\rho_F^2,
 \qquad PC\le L^{o(1)}.
\]

Paley--Zygmund gives the required D_k success parameter L^{-o(1)}. The subsequent factoring implication uses the repository's Theorem R, split/field-selection reductions, and LV6. With the stated smoothness-testing cost it has c_M=(1+theta_M)/sqrt(2). This report does not prove LV6 or turn that implication into an unconditional factoring theorem.

## 6. Referee attack and precise limits

1. **Possible objection: “you squared an L1 estimate.”** No. AR uses (3.4), the primitive large sieve with induced-character multiplicities, and H1/ZD. Its small-conductor loss is explicitly exp(O(sqrt(log N))), acceptable for an upper bound but too large for an o(1) variance. The old B5 barrier remains valid.
2. **Possible objection: “you forgot the ell factor.”** Equations (3.10)--(3.11) distinguish conductor ell and ell^2. The worst exponent is r+3lambda-xi, explicitly negative. For the original parameters this is -(h+(d-5)nu). Without that split the naive large-sieve estimate need not work.
3. **Possible objection: “non-reduced classes change density polynomially.”** They are divided by their exact gcd. Large gcd layers are polynomially negligible by a capacity bound. The small layers have the absolutely summable exponents in §3.4. Replacing all those main terms by Psi(x)/m would be incorrect.
4. **Possible objection: “your extra variable is the same for both points.”** Equation (5.2) is pointwise and uses one shared z. Only a one-coordinate residue average is required afterwards. Assuming independent multipliers for the two points would be an invalid field design; it is not done here.
5. **Possible objection: “you lost the partner density after averaging residues.”** Equation (5.3) removes the moving residue before (5.4) is applied. This is precisely why the multiplier was introduced. In the old design the analogous energy still depends on the fibre parameter in (2.1).
6. **Possible objection: “z is restricted to units / smooth numbers.”** It is uniform on the entire interval. Units are used only as a subset for the FM lower bound. The upper bound includes every z. Conditioning z on being smooth would invalidate (5.3).
7. **Possible objection: “a bigger shift broke FM.”** The new scalar lower bound (4.3) and fixed ratio H_R/Z_0 repair this. The original tiny-offset argument is not used.
8. **Possible objection: “you hid a divisor bound in L^{o(1)}.”** Polynomial-size tau and A^{omega} bounds appear only in terms with fixed polynomial savings, or with the decaying exceptional-conductor power. At the main scale they are handled by a smooth weighted Euler product, with explicit exp(O(u_R)) loss.
9. **Possible objection: “AR implies a maximum over coarse classes.”** It does not. Equations (3.12) and (3.13) are averages. LV* with a maximum over the coarse class remains unproved. The multiplier proof avoids that maximum.
10. **Remaining limitation:** the independent proof in this report has not been externally refereed. It proves analytic moments for the specified law. Its factoring consequence retains LV6, including primitive-part/leading-coefficient and transcript bookkeeping, and the repository's chosen smoothness-test cost assumption or its Pollard--Strassen replacement.

No appeal to the newer BV or BDH theorems is needed for AR. In particular [Harper 2025, Corollary 1](https://arxiv.org/pdf/2412.19644) is an all-moduli asymptotic with log^C x<=y<=x and x^{0.51}<=Q<=x, rather than the sparse-family assertion (AR); its Q Psi diagonal does not by positivity give the missing partner-weight correlation. The de la Breteche--Drappeau, Drappeau, and Pascadi results are therefore not invoked at an unchecked range.

## 7. Short verdict

* **Original (a''_k):** LD, LD*, and the partner-weighted version of piece (ii) remain **CONJECTURED** here. The stated Farey multiplicity premise is **FALSE for its admissible filtered fractions**. Their marginal structured-residue bound (3.13) is **PROVED**; the missing mixed correlation is explicitly (2.1) with its partner smoothness indicator.
* **Piece (i): PROVED** as (AR)/(3.12), an L^{o(1)} upper main-square bound, averaged over ell, smooth n, and q. It includes imprimitive characters and non-reduced classes, retains both sparsity factors, and does not assert an o(1) variance or a coarse-class maximum.
* **Multiplier design (M): LD* and filtered LD are PROVED.** FM, Theorem R and Lemma Kern remain valid by the checks above. The new finite-k budget is theta_M=3/(d+1)+d/(4k)+eta, eta>0. The limiting constant remains (1+3/(d+1))/sqrt(2), including 0.919239 at d=9.
* **Smallest remaining input for the multiplier design's factoring conclusion:** the stated LV6 adaptation for this field law and heights (and the stated smoothness-testing guarantee if its faster cost is used). No lift-dispersion hypothesis remains in the analytic second-moment argument for (M).
