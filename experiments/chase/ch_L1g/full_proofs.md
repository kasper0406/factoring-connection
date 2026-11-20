# Chase Lemma 1 (H / Branch Selection): Complete Proofs and Method Analysis

Workspace: `/Volumes/git/factoring-connection/codex/scratch/ch_L1g/`

## 1. Setup and Notation

Let $N$ be a composite integer with $x = \log N$, $F_\alpha = x^\alpha (\log x)^{1-\alpha}$, and $L_N[\alpha, c] = \exp((c + o(1)) F_\alpha)$.
Fix $0 < \beta < 1/3$ (e.g. $\beta = 1/4$), smoothness bound $B = L_N[\beta, \kappa]$ ($\kappa > 0$), and positive primitive search box
$$\Omega_{\text{prim}} = \{(u, v) \in [1, K]^2 : \gcd(u, v) = 1\}, \qquad K = L_N[1-\beta, \nu], \quad t_0 := \log(K^2) = 2\nu F_{1-\beta}.$$
In M2 (`notes/subthird_results.md:76-79`, `codex/subthird_1_gemini.out:22-30`, `experiments/hatch/hatch_h1c/report.md:34-39,54-59`), $Q(u, v) = a u^2 + b u v + c v^2$ is a positive-definite binary quadratic form of discriminant $b^2 - 4ac = -4N$ with $a, c \asymp \sqrt{N}$ and $|b| \le a$. Thus on $\Omega_{\text{prim}}$, $\log Q(u, v) = \frac{1}{2}x + O(F_{1-\beta}) = (\frac{1}{2} + o(1))x$, $u_0 := \frac{x}{2\log B} \sim \frac{1}{2\kappa}(x/\log x)^{1-\beta}$, and Dickman density is $\rho(u_0) = L_N[1-\beta, -\frac{1-\beta}{2\kappa}]$. For $2\nu > \frac{1-\beta}{2\kappa}$, the heuristic expected witness supply is $|\Omega_{\text{prim}}| \rho(u_0) = L_N[1-\beta, 2\nu - \frac{1-\beta}{2\kappa}] \gg 1$.

Let $\mathcal{P}_{\text{good}}(Q) = \{q \le B \text{ odd prime} : q \nmid aN,\ (\frac{-N}{q}) = +1\}$, $Q_{\text{good}}(u, v) = \prod_{q \in \mathcal{P}_{\text{good}}(Q)} q^{v_q(Q(u,v))}$, and $\mathcal{W}(Q) = \{(u, v) \in \Omega_{\text{prim}} : Q(u, v) \text{ is } B\text{-smooth}\}$.
For any modulus range $(D_0, D_1]$ with $D_0 \ge 2K^2$, define the set of CRT root branches:
$$\mathcal{B}_{(D_0, D_1]}(Q) = \big\{ (d, r) : D_0 < d \le D_1,\ \text{supp}(d) \subseteq \mathcal{P}_{\text{good}}(Q),\ 0 \le r < d,\ Q(r, 1) \equiv 0 \pmod d \big\}.$$

---

## 2. Theorem 1: Exact Bipartite Fiber Identity, Pointwise Random-Branch Bounds, and Reduction [PROVED]

**Theorem 1.** Let $Q(u, v) = a u^2 + b u v + c v^2$ be any positive-definite integral form of discriminant $-4N < 0$, and let $D_1 \ge D_0 \ge 2K^2$.
1. **(Exact Fiber Identity & Uniqueness).** Partial Euclidean rational reconstruction $\text{rec}(r, d)$ (stopping at the first remainder $\le K$) defines a partial map $\text{rec}: \mathcal{B}_{(D_0, D_1]}(Q) \rightharpoonup \Omega_{\text{prim}}$ in $O(\log^2 d)$ bit operations. Distinct points in $\Omega_{\text{prim}}$ have disjoint preimage fibers, and for every $(u, v) \in \Omega_{\text{prim}}$,
   $$|\text{rec}^{-1}(u, v)| = \tau_{(D_0, D_1]}\big(Q_{\text{good}}(u, v)\big) := \#\big\{d \mid Q_{\text{good}}(u, v) : D_0 < d \le D_1\big\}.$$
   In particular, the exact number of winning branches in $\mathcal{B}_{(D_0, D_1]}(Q)$ is
   $$|\mathcal{B}_{\text{win}}(Q)| = \sum_{(u, v) \in \mathcal{W}(Q)} \tau_{(D_0, D_1]}\big(Q_{\text{good}}(u, v)\big),$$
   and if $D_1 \ge B D_0$, every $(u, v) \in \mathcal{W}(Q)$ with $Q_{\text{good}}(u, v) > D_0$ satisfies $|\text{rec}^{-1}(u, v)| \ge 1$.
2. **(Polynomial-Time Search Equivalence).** Given a winning branch $(d, r) \in \mathcal{B}_{\text{win}}(Q)$ with $d > 2K^2$, $\text{rec}(r, d)$ and trial division by $q \le B$ output a certified $B$-smooth witness $(u, v) \in \mathcal{W}(Q)$ in $O(\log^2 d + \pi(B)\text{poly}(x))$ bit operations. Conversely, given any witness $(u, v) \in \mathcal{W}(Q)$ with $Q_{\text{good}}(u, v) > 2K^2$, greedy multiplication of prime-power factors of $Q_{\text{good}}(u, v)$ yields $d \in (2K^2, 2BK^2]$ and $r := u v^{-1} \bmod d$ with $(d, r) \in \mathcal{B}_{\text{win}}(Q)$ in $O(\pi(B)\text{poly}(x))$ bit operations.
3. **(Pointwise Random-Branch Upper Bounds — Individual Form $Q$, Zero Averaging).** Suppose $a, |b|, c \le C_0 \sqrt{N}$, $0 < \beta < 1/2$, $B = L_N[\beta, \kappa]$, $K = L_N[1-\beta, \nu]$.
   - **(i) Upper-band $m$-prime branches:** Let $\eta = 1/\log x$, $\mathcal{P}_{\text{up}} = \{q \in \mathcal{P}_{\text{good}}(Q) : B^{1-\eta} < q \le B\}$, and $m = \lfloor \frac{\log(2K^2)}{(1-\eta)\log B} \rfloor + 1$. Assuming $|\mathcal{P}_{\text{up}}| = B^{1-o(1)}$, uniform sampling from $\mathcal{B}_m^{\text{up}}(Q) = \{(q_1\cdots q_m, r) : q_i \in \mathcal{P}_{\text{up}} \text{ distinct}\}$ satisfies
     $$\Pr_{(d, r) \leftarrow \mathcal{B}_m^{\text{up}}(Q)}\big[\text{rec}(r, d) \in \mathcal{W}(Q)\big] \le \frac{|\mathcal{W}(Q)|}{K^2} \exp\big(O(F_{1-2\beta}) + o(F_{1-\beta})\big).$$
   - **(ii) Unrestricted prime-power branches:** Assuming $|\mathcal{P}_{\text{good}}(Q) \cap (B/4, B/2]| = B^{1-o(1)}$, uniform sampling from the full prime-power branch set $\mathcal{B}_{(2K^2, 2BK^2]}(Q)$ satisfies, for every individual form $Q$,
     $$\Pr_{(d, r) \leftarrow \mathcal{B}_{(2K^2, 2BK^2]}(Q)}\big[\text{rec}(r, d) \in \mathcal{W}(Q)\big] \le |\mathcal{W}(Q)| \cdot K^{-2(1-\beta) + o(1)}.$$

**Proof.**
(1) If $(u_1, v_1), (u_2, v_2) \in \Omega_{\text{prim}}$ both satisfy $u_i \equiv r v_i \pmod d$ with $d > 2K^2$, then $d \mid (u_1 v_2 - u_2 v_1)$ while $|u_1 v_2 - u_2 v_1| < K^2 < d$, forcing $u_1 v_2 = u_2 v_1$ and hence $(u_1, v_1) = (u_2, v_2)$ by primitivity and $v_i > 0$. Moreover, if $u - r v = -k d$ with $(u, v) \in \Omega_{\text{prim}}$ and $0 \le r < d$, then $0 \le k < v$ and $|\frac{r}{d} - \frac{k}{v}| = \frac{u}{d v} \le \frac{K}{d v} < \frac{1}{2 K v} \le \frac{1}{2 v^2}$ since $d > 2K^2$. By Legendre's continued-fraction theorem, $k/v$ is a convergent of $r/d$; in the Euclidean remainder sequence $(R_j, V_j)$ starting at $(d, 0), (r, 1)$ with $R_j \equiv r V_j \pmod d$, the unique convergent with $V_j = v \le K$ is the first remainder $R_j = u \le K$ (as adjacent convergents satisfy $R_{j-1} V_j + R_j V_{j-1} = d > 2K^2$, so $R_{j-1} > K$). Hence $\text{rec}(r, d) = (u, v)$. Conversely, for $(u, v) \in \Omega_{\text{prim}}$ and $d \mid Q_{\text{good}}(u, v)$ with $D_0 < d \le D_1$, no $q \mid d$ can divide $v$ (else $q \mid Q(u,v) - (bu + cv)v = a u^2$ with $q \nmid a$, forcing $q \mid u$). At each $q^e \parallel d$, Hensel lifting of the two simple roots of $Q(X, 1) \bmod q$ gives $Q(u, v) \equiv a(u - r_{q^e,0} v)(u - r_{q^e,1} v) \pmod{q^e}$ with $r_{q^e,0} - r_{q^e,1} \in (\mathbb{Z}/q^e\mathbb{Z})^\times$, so $q^e \mid Q(u, v)$ iff $u v^{-1} \equiv r_{q^e,\sigma} \pmod{q^e}$ for a unique $\sigma \in \{0, 1\}$. CRT assembles a unique $r = u v^{-1} \bmod d$ with $(d, r) \in \mathcal{B}_{(D_0, D_1]}(Q)$. Summing fiber sizes over $(u, v) \in \mathcal{W}(Q)$ gives the exact identity. If $Q_{\text{good}}(u, v) > D_0$ and $D_1 \ge B D_0$, multiplying prime factors of $Q_{\text{good}}(u,v)$ ($\le B$ each) until the product first exceeds $D_0$ lands in $(D_0, B D_0] \subseteq (D_0, D_1]$.

(2) Immediate from (1).

(3)(i) Every $(d, r) \in \mathcal{B}_m^{\text{up}}(Q)$ satisfies $d > B^{(1-\eta)m} \ge 2K^2$. Since $m \sim \frac{2\nu}{\kappa}(x/\log x)^{1-2\beta} = B^{o(1)}$ while $|\mathcal{P}_{\text{up}}| = B^{1-o(1)}$, $\log|\mathcal{B}_m^{\text{up}}(Q)| = m\log(2|\mathcal{P}_{\text{up}}|/m) + O(m) = \frac{2\nu}{1-\eta} F_{1-\beta} - o(F_{1-\beta}) = 2\nu F_{1-\beta} + o(F_{1-\beta})$ for $\eta = 1/\log x$. For any $(u, v) \in \Omega_{\text{prim}}$, $\log Q(u, v) \le \frac{1}{2}x + 2\log K + O(1)$, so $\omega_{\text{up}}(Q(u, v)) \le \frac{\log Q(u,v)}{(1-\eta)\log B} = O(u_0)$. Because $u_0 / m \asymp (x/\log x)^\beta$, each witness lies on at most $\binom{\omega_{\text{up}}}{m} \le (e \omega_{\text{up}} / m)^m = \exp(O(m\log x)) = \exp(O(F_{1-2\beta}))$ branches in $\mathcal{B}_m^{\text{up}}(Q)$.

(3)(ii) For any integer $M \le 3 C_0 \sqrt{N} K^2$, $\omega(M) \le (1 + o(1))\frac{\log M}{\log\log M} = O(x/\log x)$. Letting $p_1 < p_2 < \dots$ be all primes in increasing order, $p_{\omega(M)} = O(x)$ and Rankin's bound at exponent $\beta \in (0, 1/2)$ gives
$$\tau_{\le 2BK^2}(M) \le (2BK^2)^\beta \prod_{p \mid M} (1 - p^{-\beta})^{-1} \le (2BK^2)^\beta \prod_{i=1}^{\omega(M)} (1 - p_i^{-\beta})^{-1} \le K^{2\beta} \exp\Big(\beta \log(2B) + O\big(x^{1-\beta}/\log x\big)\Big) = K^{2\beta + o(1)},$$
since $x^{1-\beta}/\log x = F_{1-\beta} (\log x)^{-1-\beta} = o(F_{1-\beta})$. For the denominator, partition $(B/4, B/2]$ into $J = \lceil 4 \log(2K^2)/\log B \rceil + 4 = O(x^{1-2\beta})$ log-equal intervals $I_j = (e^{\mu_j - \delta}, e^{\mu_j}]$ of log-width $\delta = \frac{\log 2}{J}$. By pigeonhole, some $I_{j^*}$ contains $S \ge |\mathcal{P}_{\text{good}}(Q) \cap (B/4, B/2]| / J = B^{1-o(1)}$ primes. Choosing $m^* = \lfloor \frac{\log(2K^2)}{\mu_{j^*} - \delta} \rfloor + 1$ distinct primes in $I_{j^*}$ gives $\log d > m^*(\mu_{j^*} - \delta) \ge \log(2K^2)$ and $\log d \le m^* \mu_{j^*} \le \log(2K^2) + \mu_{j^*} + m^*\delta < \log(2K^2) + \log(B/2) + \log 2 = \log(2BK^2)$ (for $x$ large enough that $m^*\delta = m^* \frac{\log 2}{J} < \log 2$). Thus $|\mathcal{B}_{(2K^2, 2BK^2]}(Q)| \ge 2^{m^*} \binom{S}{m^*} = K^{2 - o(1)}$, yielding the claim. $\square$

---

## 3. Theorem 2: Two-Choice Projective Weighted Theorem J for Quadratic Norms [PROVED]

Theorem 1(3)(ii) leaves a $K^{2\beta}$ worst-case gap if a hypothetical witness were a primorial of tiny primes, and applies to uniform branch sampling. We now prove the full M2 extension of Weighted Theorem J (`experiments/hatch/hatch_h1c/report.md:9-27`, `notes/sihard_results.md:428`), obtaining the **exact Dickman exponent $-\frac{1-\beta}{2\kappa}$** uniformly across all fixed nonnegative depth-weight profiles $(h_{qk})_{q \le B, k \ge 1}$ — including simultaneously certifying **all $2^{\omega(d)}$ root branches** of every modulus $d \ge K^2$ in any $\exp(o(F_{1-\beta}))$-sized repertoire!

**Theorem 2.** Fix $0 < \beta < 1$, $\kappa > 0$, $B = L_{N_0}[\beta, \kappa]$, $x = \log N_0$, $F_a = x^a(\log x)^{1-a}$, and $2 \le K \le \exp(O(F_{1-\beta}))$ with $t_0 = \log(K^2)$. Let $a_0 = \lceil \frac{1}{2}\sqrt{N_0} \rceil$ and consider the balanced positive-definite form family
$$\mathcal{F}(N_0) = \big\{ Q(u, v) = a u^2 + 2 b_0 u v + c v^2 : a \in [a_0, 2a_0),\ b_0 \in (-a/2, a/2],\ c \in [4a_0, 8a_0) \big\},$$
where every $Q \in \mathcal{F}(N_0)$ has negative discriminant $-4N = 4(b_0^2 - ac)$ with $N \in [\frac{15}{4} a_0^2, 16 a_0^2] \subset [\frac{15}{16} N_0, 5 N_0]$ and $a \asymp c \asymp \sqrt{N}$.
For any nonnegative weights $(h_{qk})_{q \le B, k \ge 1}$ with $E = \sum_{q, k} h_{qk}^2 \log q > 0$, define $S(z) = \sum_{q \le B} \log q \sum_{k \le v_q(z)} h_{qk}$, and grant at unit cost every $(u, v) \in \Omega_{\text{prim}}$ with $S(Q(u, v)) \ge \sqrt{E t_0}$ (strictly below the two-choice projective Johnson benchmark $J_2 = \sqrt{2 E t_0}$). Uniformly over all weight profiles:
$$\mathbb{E}_{Q \leftarrow \mathcal{F}(N_0)}\big[ \#\big\{ (u, v) \in \Omega_{\text{prim}} : Q(u, v) \text{ is } B\text{-smooth and } S(Q(u, v)) \ge \sqrt{E t_0} \big\} \big] \le \exp\left( -\left(\frac{1-\beta}{2\kappa} - o(1)\right) F_{1-\beta} \right).$$

**Corollary 2.1 (Simultaneous $2^{\omega(d)}$-Branch and Portfolio Bound).** For any $B$-smooth modulus $d = \prod_{q \le B} q^{e_q} \ge K^2$, setting $h_{qk} = \mathbf{1}_{k \le e_q}$ gives $E = \log d \ge t_0$, so every $(u, v) \in \Omega_{\text{prim}}$ lying on *any* of the $2^{\omega(d)}$ CRT root branches of $Q \bmod d$ satisfies $S(Q(u, v)) \ge \log d \ge \sqrt{E t_0}$. Consequently, any fixed or adaptive algorithm selecting from a portfolio of $|\mathcal{R}| \le \exp(o(F_{1-\beta}))$ moduli $d \ge K^2$ or weight profiles yields at most $\exp[-((1-\beta)/(2\kappa) - o(1)) F_{1-\beta}]$ smooth witnesses per call in expectation over $\mathcal{F}(N_0)$, even if all $2^{\omega(d)}$ branches of each modulus are reconstructed at unit cost.

**Proof of Theorem 2.**
Step 1 (Point-value multiplicity bound across $\mathcal{F}(N_0)$). We have $|\mathcal{F}(N_0)| = \sum_{a=a_0}^{2a_0-1} a \cdot 4a_0 \ge 4 a_0^3$. For $Q \in \mathcal{F}(N_0)$ and $(u, v) \in \Omega_{\text{prim}}$ (with $K \le a_0$ for large $N_0$), $|2b_0 u v| \le a u v \le \frac{1}{2}a(u^2 + v^2)$ and $c \ge 4a_0 \ge 2a$, so
$$\tfrac{1}{2}a_0(u^2 + 3v^2) \le Q(u, v) \le 2a_0 K^2 + 2a_0 K^2 + 8a_0 K^2 = 12 a_0 K^2 =: Z_{\max}.$$
Fix an integer $z \in [2a_0, Z_{\max}]$, $a \in [a_0, 2a_0)$, and $v \in [1, K]$. Let $\mu_{a,v}(z)$ be the number of triples $(u, b_0, c)$ with $u \in [1, K]$, $\gcd(u, v) = 1$, $b_0 \in (-a/2, a/2]$, $c \in [4a_0, 8a_0)$, and $a u^2 + 2b_0 u v + c v^2 = z$.
First, for fixed $u$ with $\gcd(u, v) = 1$, $v \mid (z - a u^2)$ is necessary, and any two solutions $(b_0, c), (b_0', c')$ of $2u b_0 + v c = (z - a u^2)/v$ satisfy $(v / \gcd(2, v)) \mid (b_0' - b_0)$; since $b_0$ ranges over an interval of length $a$ and $v \le K \le a$, there are at most $\lfloor \frac{a}{v/2} \rfloor + 1 \le \frac{3a}{v}$ pairs $(b_0, c)$ per $(u, v)$.
Second, $z - a u^2 = 2b_0 u v + c v^2 \in [4a_0 v^2 - a K v,\ 8a_0 v^2 + a K v]$, so $u^2$ lies in an interval $[U_0, U_0 + W_v]$ of width $W_v = \frac{4a_0 v^2 + 2a K v}{a} \le 4v^2 + 2Kv \le 6Kv$.
- If $z \ge 20 a_0 K v$, then since $v \le K$, $8a_0 v^2 + a K v \le 10 a_0 K v \le z/2$, so $U_0 \ge \frac{z}{2a} \ge \frac{z}{4a_0}$ and the interval of positive $u$ with $u^2 \in [U_0, U_0 + W_v]$ has length $\Delta u = \sqrt{U_0 + W_v} - \sqrt{U_0} \le \frac{W_v}{2\sqrt{U_0}} \le \frac{6Kv}{\sqrt{z/a_0}}$.
- If $z < 20 a_0 K v$, then $\Delta u \le \sqrt{W_v} \le \sqrt{6Kv} < \sqrt{6 K \cdot \frac{20 a_0 K^2}{z}}$? More directly, $v > \frac{z}{20 a_0 K}$ implies $\frac{\sqrt{6Kv}}{v} = \sqrt{\frac{6K}{v}} < \sqrt{\frac{120 a_0 K^2}{z}} = \frac{\sqrt{120}\, K}{\sqrt{z/a_0}}$.
In both cases, $\frac{\Delta u}{v} \le \frac{12 K}{\sqrt{z/a_0}}$. Now let $g = \gcd(a, v)$ and $m' = v/g$. Necessity of $a u^2 \equiv z \pmod v$ with $\gcd(u, v) = 1$ requires $g \mid z$ and $u^2 \equiv (a/g)^{-1}(z/g) \pmod{m'}$, which restricts $u \bmod m'$ to at most $4 \cdot 2^{\omega(m')} \le 4 \cdot 2^{\omega(v)}$ residue classes. An interval of length $\Delta u$ contains at most $(\frac{\Delta u}{m'} + 1) \cdot 4 \cdot 2^{\omega(v)} \le 4 \cdot 2^{\omega(v)} (1 + \frac{12 K g}{\sqrt{z/a_0}})$ valid $u$. Multiplying by $\frac{3a}{v} \le \frac{6a_0}{v}$ and summing over $a \in [a_0, 2a_0)$ (using $\sum_{a=a_0}^{2a_0-1} \gcd(a, v) \le \sum_{d \mid v} d (\frac{a_0}{d} + 1) \le 2 a_0 \tau(v)$ since $v \le a_0$) and $v \in [1, K]$ yields total multiplicity $\mu(z) := \sum_{a, v} \mu_{a,v}(z)$:
$$\mu(z) \le 24 a_0^2 \sum_{v=1}^K \frac{2^{\omega(v)}}{v} + \frac{576 a_0^{5/2} K}{\sqrt{z}} \sum_{v=1}^K \frac{\tau(v) 2^{\omega(v)}}{v} \ll a_0^2 (\log K)^4 \left(1 + \frac{K}{\sqrt{z/a_0}}\right).$$

Step 2 (Weighted Rankin and exponential Markov). Rescale $(h_{qk})$ so $E = t_0$, making the certificate condition $S(z) \ge t_0$. Put $s = \log a_0 = \frac{1}{2}x + O(1)$, $u_0 = s / \log B$, and $\theta = 1 - \frac{\log u_0 - 2\log\log x}{\log B} \to 1$ (so $\theta > 3/4$ for large $N_0$). For every $z \in [2a_0, Z_{\max}]$ ($Z_{\max} = 12 e^{s + t_0}$):
- First term of $\mu(z)/|\mathcal{F}(N_0)|$: $\frac{a_0^2 (\log K)^4}{4 a_0^3} \le (\log K)^4 e^{-s} \cdot \left(\frac{12 e^{s+t_0}}{z}\right)^\theta \ll (\log K)^4 e^{-(1-\theta)s + \theta t_0} z^{-\theta}$.
- Second term of $\mu(z)/|\mathcal{F}(N_0)|$: $\frac{a_0^{5/2} K (\log K)^4}{4 a_0^3 z^{1/2}} = \frac{1}{4}(\log K)^4 e^{-s/2 + t_0/2} z^{-1/2} \cdot \left(\frac{12 e^{s+t_0}}{z}\right)^{\theta - 1/2} \ll (\log K)^4 e^{-(1-\theta)s + \theta t_0} z^{-\theta}$, where $+t_0/2$ cancels $-t_0/2$ identically!

Hence, multiplying by $e^{2\theta(S(z) - t_0)} \ge 1$ on certified $z$ and summing over all $B$-smooth $z$,
$$\log \mathbb{E}_{Q \leftarrow \mathcal{F}(N_0)}[\#\text{hits}] \le -(1-\theta)s + \theta t_0 - 2\theta t_0 + \sum_{q \le B} \log M_q + O(\log x), \qquad M_q = \sum_{j \ge 0} q^{\theta(2\sum_{k\le j} h_{qk} - j)}.$$
By `experiments/hatch/hatch_h1c/report.md:24-26` (valid verbatim with $F = F_{1-\beta}$ and $\sum_q D_q \log q = E = t_0 = O(F_{1-\beta})$), $\sum_{q \le B} \log M_q \le \theta t_0 + o(F_{1-\beta})$. The $t_0$ terms cancel identically ($\theta t_0 - 2\theta t_0 + \theta t_0 = 0$), leaving $-(1-\theta)s + o(F_{1-\beta}) = -\left(\frac{1-\beta}{2\kappa} - o(1)\right) F_{1-\beta}$. $\square$

---

## 4. Structural Obstructions for the Four Adaptive Mechanisms

### 4.1 Meet-in-the-Middle over Root Branches (`MITM_COLLIDE`) — Proposition MITM [PROVED in random-branch / fixed-split model]

Split $\mathcal{P}_{\text{good}}(Q) = \mathcal{P}_L \sqcup \mathcal{P}_R$ into two halves and sample $M$ branches in $\mathcal{B}_m^{\text{up}}(\mathcal{P}_L)$ and $M$ branches in $\mathcal{B}_m^{\text{up}}(\mathcal{P}_R)$ ($d_L, d_R > 2K^2$). Each colliding candidate $(u, v) = \text{rec}(r_L, d_L) = \text{rec}(r_R, d_R)$ indeed has forced divisor $d_L d_R > 4K^4$. However, for every witness $(u, v) \in \mathcal{W}(Q)$, Theorem 1(3)(i) bounds its left-degree in each half by $\exp(O(F_{1-2\beta}))$, so its probability of appearing in $\mathcal{L}_L$ (and independently in $\mathcal{L}_R$) is at most $M \cdot K^{-2} \exp(O(F_{1-2\beta}) + o(F_{1-\beta}))$. Summing over all $|\mathcal{W}(Q)|$ witnesses,
$$\mathbb{E}[\#\text{witness collisions}] \le |\mathcal{W}(Q)| \cdot M^2 K^{-4} \exp\big(o(F_{1-\beta})\big) = \left(M \frac{|\mathcal{W}(Q)|}{K^2}\right) \cdot \frac{M}{K^{2 - o(1)}}.$$
For every budget $M \ll K^2 = L_N[1-\beta, 2\nu]$ (in particular $M = L_N[\beta, O(1)]$), the factor $M / K^{2-o(1)} \ll 1$ makes two-sided branch MITM **strictly inferior in exponent** to testing the $M$ left-branches directly (mirroring `notes/crtdecode_results.md:443` PH-$\theta$). Conversely, splitting the full norm $Q(u,v) = d_L d_R$ at $d_L, d_R \approx N^{1/4} K$ puts $d_L \gg K^2$, where in-box branch density is $O(K^2 / d_L) = N^{-1/4+o(1)}$, requiring exponential list sizes $M = N^{1/4-o(1)}$ (`notes/sihard_results.md:111` Prop MM).

### 4.2 Lattice Enumeration with Partial-Norm Pruning (`BEAM_SV`, `BEAM_QNORM`) — Lemma Non-Visibility [PROVED]

**Lemma (Witness Non-Visibility Before $\frac{3}{4} K_0^2$).** Let $(u^*, v^*) \in \Omega_{\text{prim}}$ satisfy $\|(u^*, v^*)\|_\infty \ge K_0$, and let $(d_k, r_k)$ be any ancestor branch of $(u^*, v^*)$ with $d_k < \frac{3}{4} K_0^2$. Then the shortest Euclidean nonzero vector $\mathbf{b}_1$ of $\Lambda_{d_k, r_k}$ satisfies $\|\mathbf{b}_1\|_\infty \le \|\mathbf{b}_1\|_2 \le (4/3)^{1/4}\sqrt{d_k} < K_0 \le \|(u^*, v^*)\|_\infty$, so $\mathbf{b}_1 \ne \pm(u^*, v^*)$ and $\mathbf{b}_1, (u^*, v^*)$ are $\mathbb{Z}$-linearly independent.

**Consequence.**
1. Pruning intermediate branches by smallest $\|\mathbf{b}_1\|_\infty$ (`round2_crt.py`) only retains branches where an *already tiny* vector $\mathbf{w}$ ($\|\mathbf{w}\|_\infty \le K_{\text{eff}} \ll K$) happens to divide successive batch primes; it is therefore disguised presieving of the tiny sub-box $[1, K_{\text{eff}}]^2$. Asymptotically at $2\nu = \frac{1-\beta}{2\kappa} + 2\varepsilon$, any sub-box $[1, K^{1-\delta}]^2$ has expected smooth count $K^{2(1-\delta)}\rho(u_0) = L_N[1-\beta, 2\varepsilon - 2\delta\nu] \to 0$ once $\delta > \varepsilon/\nu$, and even on nonempty boxes $2\log K = o(x)$ changes $\log Q \sim \frac{1}{2}x$ only at subleading order $L_N[1-2\beta]$.
2. Pruning by the exact reduced partial cofactor $A_k = \min_{\Lambda_{d_k, r_k}\setminus 0} Q(u, v)/d_k$ (`BEAM_QNORM`) computes the minimum of the reduced quadratic form in the ideal class $[\mathfrak{a}\mathfrak{d}_k] \in \mathrm{Cl}(-4N)$. Across all reduced forms of discriminant $-4N$, at most $\sum_{A_1 \le Y} 2^{\omega(\gcd(A_1, N))+1} = Y N^{o(1)}$ classes out of $N^{1/2\pm o(1)}$ have $A_1 \le Y$; hence any $L_N[\beta]$ beam encounters only $A_k = N^{1/2-o(1)}$, and before $d_k \approx K^2$ the minimizing vector is not $(u^*, v^*)$.

### 4.3 Sieving-in-the-Branch-Tree (`TREE_SIEVE`) — Lemma Tree-Blindness [PROVED]

**Lemma (Zero Lookahead Before $K^2/B$).** Let $(d_k, r_k)$ be any branch with $d_k \le K^2 / B$. For *every* good split prime $q \le B$ with $q \nmid d_k$ and *both* roots $\sigma \in \{0, 1\}$ of $Q(X, 1) \bmod q$, the extended lattice $\Lambda_{d_k q, r_{k,\sigma}'}$ has determinant $d_k q \le K^2$ and hence (by Minkowski's convex-body theorem on $[-K, K]^2$, area $4K^2 \ge 4 d_k q$) contains a nonzero primitive point in $[-K, K]^2$. Moreover, whenever $\lambda_1(\Lambda_{d_k q, r'}) \ge L$, standard lattice-point counting in $[1, K]^2$ gives $\#(\Lambda_{d_k q, r'} \cap \Omega_{\text{prim}}) = \frac{6}{\pi^2}\frac{K^2}{d_k q} + O(\frac{K \log K}{L} + \log K)$.

**Consequence.** At every depth $d_k \ll K^2/B$, *no* good prime $q \le B$ is eliminated by box-emptiness and every $q \nmid d_k$ hits $\frac{2}{q}(1 + o(1))$ of the primitive box points of $\Lambda_{d_k, r_k}$ independently of $r_k$. Non-trivial lookahead begins only at penultimate depth $d_{m-1} \in (2K^2/B, 2K^2]$, where enumerating and sieving $\Lambda_{d_{m-1}, r_{m-1}} \cap \Omega_{\text{prim}}$ ($\Theta(B)$ points) is an ordinary special-$d_{m-1}$ lattice sieve (`notes/smoothenum_results.md:198` Prop LQ), costing $\Theta(B \log\log B)$ work for $\Theta(B \rho(u_0 + 1))$ expected hits, i.e. $\rho(u_0)^{1+o(1)}$ hits per unit work.

### 4.4 Proxy Sampling & Local Branch-Walk (`PROXY_TILT`, `MCMC_PROXY`) — Lemma CRT-Scramble [PROVED]

**Lemma (CRT-Scramble on Penultimate Neighbors).** Let $(d', r')$ be any sub-branch with $d' > 2K^2 / B$, and let $(u_1, v_1) \ne (u_2, v_2)$ be two distinct points in $\Lambda_{d', r'} \cap \Omega_{\text{prim}}$.
1. $u_1 v_2 - u_2 v_1 = k d'$ for a nonzero integer $k$ with $1 \le |k| < \frac{K^2}{2K^2/B} = \frac{B}{2}$.
2. If a good split prime $p \nmid d'$ vanishes on the *same* root $r_{p,\sigma}$ at both $(u_1, v_1)$ and $(u_2, v_2)$, then $p \mid ((u_1 - r_{p,\sigma} v_1)v_2 - (u_2 - r_{p,\sigma} v_2)v_1) = u_1 v_2 - u_2 v_1 = k d'$, so $p \mid k$.
3. In particular, **zero** good split primes $p \in [B/2, B] \setminus \text{supp}(d')$ can vanish on the same root branch at $(u_1, v_1)$ and $(u_2, v_2)$, and at most $\omega(|k|) \le (1+o(1))\frac{\log B}{\log\log B} \ll u_0$ primes $< B/2$ outside $d'$ can share a branch.

**Consequence.** Single-prime Markov swaps $(q_i, \sigma_i) \to (q_i', \sigma_i')$ on $\mathcal{B}_m$ either stay on the identical $(u, v)$ (yielding zero fresh candidates) or jump to $(u_2, v_2)$ with $|k| < B/2$, scrambling all prime-branch incidences in $[B/2, B]$ and destroying any smooth local gradient toward $\mathcal{W}(Q)$. Meanwhile, factorized prime-weight proxies $w(d) = \prod w(q_i)$ at $d \asymp 2K^2$ leave the untested cofactor scale $Q(u,v)/d \asymp \sqrt{N}$ unchanged and remain bounded by Theorem 2 / Dickman density $\rho(u_0)^{1-o(1)}$.

---

## 5. Pre-Registered Experimental Record and Disclosed Post-Hoc Diagnosis

Files: `prereg.md` (written prior to running), `experiment.py`, `results.jsonl`, `summary.json`, `posthoc_diag.py`, `posthoc.json`.

### 5.1 Pre-Registered Scorecard (`summary.json`)
- **P1 (40-bit exact bipartite identity):** **PASS (`4/4` exact, 0 discrepancies)**: winning branch counts `0, 6, 0, 14` equal summed good-divisor fiber counts `0, 6, 0, 14` across `11,698` traversed branches.
- **P2 (72/84/96-bit fixed-work main-shell comparison at `W_budget = 4000` per form, `72,000` pooled work):**
  - Against naive unwindowed `RND_BRANCH` (`3` hits / `3` main-shell hits, `34.08 ms`): `PROXY_TILT` (`3/3`), `BEAM_SV` (`4/4`, `1.33x`), `BEAM_QNORM` (`5/4`, `1.33x`), `MCMC_PROXY` (`2/2`, `0.67x`), and `MITM_COLLIDE` (`1/0`, `0` collisions across all 30 instances) all stayed well below `2.0x`.
  - **Disclosed P2 failure as scored:** `TREE_SIEVE` scored `13/12` (`4.0x` naive `RND_BRANCH`), caused by (a) an omitted per-prime sublattice progression setup charge in `experiment.py:444` (`TREE_SIEVE` spent `217.49 ms` vs `RND_BRANCH`'s `34.08 ms`) and (b) modulus overshoot up to `2 B K^2` in unwindowed `RND_BRANCH` (`431` in-box points).
- **P3 (`BEAM_SV` vs `TINY_BOX_CTRL` at 72–96 bits):** **Disclosed failure as scored (`4` vs `2` hits = `2.0x > 1.25x`).** Because `2K^2 > B^2` prevented `BEAM_SV` from collapsing into a tiny sub-box (`4/4` hits had `max(u,v) > K/4`), unpresieved `TINY_BOX_CTRL` lacked the `d > 2K^2` divisor boost (`0.0028%` vs `0.696%` per test), while `BEAM_SV` (`4` hits, `135.31 ms`) merely matched presieved `RND_BRANCH` (`3` hits, `34.08 ms`).
- **P4 (Step-1 partial-cofactor vs terminal smoothness correlation, 72–96 bits):** **PASS (`r = -0.00257`, `|r| < 0.10`, `n = 1,559`).**
- **P5 (CRT-Scramble check on `p > B/2`):** **PASS (`0` violations in `14,788` distinct penultimate-neighbor pairs).**

### 5.2 Disclosed Post-Hoc Check (`posthoc.json`)
- Charging 1 progression-setup op per prime per sublattice (`TREE_SIEVE_HONEST`) and windowing only the final prime of random branch choice into `(2K^2/d_{m-1}, 4K^2/d_{m-1}]` (`RND_WINDOW`) at identical `72,000` pooled work (72–96 bits):
  - `RND_WINDOW`: **`10` hits (`9` main-shell)**, `2,216` tested, `57.84 ms`
  - `BEAM_QNORM`: `5` hits (`4` main-shell), `102.51 ms`
  - `BEAM_SV`: `4` hits (`4` main-shell), `131.60 ms`
  - `TREE_SIEVE_HONEST`: `1` hit (`1` main-shell), `26.46 ms`
- At penultimate depth `m - 1`, `r = -0.05031` (`n = 1,674` overall) and `r = -0.09806` (`n = 725` with `d_{m-1} < 3(K/4)^2/4`, where `0 / 725` shortest vectors equaled the terminal point).
