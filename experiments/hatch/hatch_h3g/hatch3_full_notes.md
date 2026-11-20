# Hatch 3: R / Thin Dependencies, Shared Descent DAGs, and Related-Input Reuse — Full Proofs & Ledgers

Date: 2026-11-20
Workspace: `/Volumes/git/factoring-connection/codex/scratch/hatch_h3g`
Context references:
- `notes/lowerbound_results.md` (Priced valuation-probe model: H, R, Credit, Reuse; Theorem C merged; Hatch 3 at lines 36–41, 87–98, 117)
- `codex/lowerbound_1_codex_proof.md`, `codex/lowerbound_2_codex_merged.md`, `codex/lowerbound_2_gemini.out`
- `notes/zroute_results.md` (Lemma XK-Z, Prop PIN, Prop ZA, Lemmas HR, DD1, WR)
- `notes/heurnfs_results.md` (§3.1 Assumption T bottom bound, §3.4–3.5 multipliers and rational sides)
- `notes/detGRH_results.md` (Prop FL, Prop Z, Prop EQ)
- `notes/fruit_results.md` (Lemma HY, Theorem HYB, Prop IV-LV)
- `notes/smoothbuild_results.md` (Lemmas KI, LP, AT; CGPT/BBM)
- `notes/smoothenum_results.md` (Table 1 T1–T8, Prop LQ, Prop FT)
- `notes/escapeR2_results.md` (Lemma TI, Cor EX, Example TB tied fields)
- `notes/selfred_results.md` (Theorems SR1, SR2, Lemma SIG)

Status labels used throughout: **[proved]** (unconditional mathematical proof), **[model-proved]** (exact mathematical deduction inside the named heuristic/probe model), **[heuristic]**, **[speculative]**, **[measured]** (pre-registered toy Python run).

---

## 1. Precise Statement Attacked

In `notes/lowerbound_results.md:36-41, 117`, the unified $L_N[1/3, (4A/3)^{2/3}]$ relation-model lower bound rests on three structural clauses in Escape Hatch 3:
1. **Quota (Axiom R, `notes/lowerbound_results.md:36`, `codex/lowerbound_1_codex_proof.md:15`):** Certified factor extraction requires $h$ distinct fresh probe successes in an opened block at smoothness bound $B = e^\mu$ with $\ln h \ge (1 - o(1))\mu - O(1)$.
2. **Tree Incidence (Assumption T, `notes/lowerbound_results.md:53, 87-92`, `notes/heurnfs_results.md:184-216`):** Multi-level descent / special-$q$ elimination assumes each large prime $> B^{O(1)}$ has degree $\le 2$ (one parent, one eliminator); Round 2's single-witness fallback (`codex/lowerbound_2_codex_merged.md:28`) dropped the constant to $(\kappa/18)^{1/3}$ and needed an extra non-cusp hypothesis `H_cusp`.
3. **Reuse Budget (`notes/lowerbound_results.md:40`, `codex/lowerbound_1_codex_proof.md:11, 54, 110`):** Target-credited probes satisfy $\sum_{i=1}^K Q_i \le e^{o(F_{1/3})} W$, ruling out a single witness serving many related inputs.

We attack all three parts — **(a) thin dependencies**, **(b) shared descent DAGs**, and **(c) related-input reuse for a single given $N$** — asking in each case whether an algorithm can break the hypothesis to reach $L_N[\alpha]$ with $\alpha < 1/3$, or whether a sharp new theorem closes that part of Hatch 3.

---

## 2. Attack on (a): Thin Dependencies and True Minimal Row Count

### 2.1 True Algebraic Minimal Row Count (`Lemma MIN-ROW` and `Lemma PIN-1`) [proved]

**Is $h > \pi(B)$ ("rows > active columns") a pure linear-algebra necessity?**
**No [proved].** Algebraically, the minimal row count for a certified relation split is **$h_{\min} = 1$** both in the $\mathbb{F}_2$ congruence-of-squares route and in the randomized $\mathbb{Z}$-route!

- **$\mathbb{F}_2$ route ($h_{\min} = 1$):** One relation whose rational and algebraic components are already squares ($U^2 \equiv \psi(V)^2 \pmod N$ with $U \not\equiv \pm \psi(V)$, e.g. Fermat $x^2 - N = y^2$ or a single square-norm NFS pair) is a rank-0 external-column dependency of weight $1$. More generally, if $h$ rows touch a set $\mathcal{E}$ of external columns that must vanish mod $2$, rank-nullity gives $h \ge \operatorname{rank}_{\mathbb{F}_2}(R_{\mathcal{E}}) + 1$; $\operatorname{rank}_{\mathbb{F}_2}(R_{\mathcal{E}})$ can be $0$ even when $B$ is large.

- **Randomized 1-Coordinate $\mathbb{Z}$-Route (`Lemma PIN-1`, [proved]; sharpens `notes/zroute_results.md:211` Prop PIN from $|Q_0| \sim 2\ln^2 N/\ln\ln N$ coordinates to $1$ coordinate):**
  Let $N = \prod_{i=1}^\omega p_i^{a_i}$ be odd with $\omega \ge 2$ distinct prime factors.
  1. *(Blum / equal 2-adic totient valuations, deterministic 1-coordinate PIN):* Suppose $v_2(\varphi(p_1^{a_1})) = \dots = v_2(\varphi(p_\omega^{a_\omega})) = e \ge 1$ (in particular every Blum semiprime $p \equiv q \equiv 3 \pmod 4$, where $e = 1$) and $a_1, \dots, a_\omega$ are odd. Let $q_0$ be **any single integer** with Jacobi symbol $(q_0 / N) = -1$. Then for **every** positive integer $D$ satisfying $q_0^D \equiv 1 \pmod N$ (or $q_0^D \equiv \pm 1 \pmod N$, after replacing $D \mapsto 2D$), Miller's power chain on base $q_0$ alone splits $N$ with probability $1$.
  2. *(General odd composite $N$, randomized 1-coordinate PIN):* Let $b \in (\mathbb{Z}/N\mathbb{Z})^\times$ be drawn uniformly at random. Then with probability $\ge 1 - 2^{1 - \omega} \ge 1/2$ over $b$, **every** positive integer $D$ satisfying $b^D \equiv 1 \pmod N$ (even if chosen adversarially knowing the factorization of $N$) splits $N$ via Miller's power chain on base $b$ alone.

  *Proof of `Lemma PIN-1`.* Let $r_i = \operatorname{ord}_{p_i^{a_i}}(b)$ and $k_i = v_2(r_i)$. Write any $D > 0$ with $b^D \equiv 1 \pmod N$ as $D = 2^s t$ with $t$ odd (note doubling $D \mapsto 2D$ increments $s$ and leaves $t$ unchanged). Since $r_i \mid 2^s t$ and $r_i / 2^{k_i}$ is odd, $r_i / 2^{k_i} \mid t$, hence $\operatorname{ord}_{p_i^{a_i}}(b^t) = 2^{k_i}$ — **independent of which multiple $D$ of $\operatorname{ord}_N(b)$ is supplied**. Whenever $k_1, \dots, k_\omega$ are not all equal, let $k_{\max} = \max_i k_i \ge 1$; then $y := b^{2^{k_{\max}-1} t}$ satisfies $y \equiv -1 \pmod{p_i^{a_i}}$ for $k_i = k_{\max}$ (since $(\mathbb{Z}/p_i^{a_i})^\times$ is cyclic for odd $p_i$) and $y \equiv +1 \pmod{p_j^{a_j}}$ for $k_j < k_{\max}$, so $1 < \gcd(y - 1, N) < N$.
  For (1), $(q_0/N) = -1$ forces $(q_0/p_i) = -1$ (so $k_i = e$) for some $i$ and $(q_0/p_j) = +1$ (so $k_j \le e - 1$) for some $j$, giving $k_i \ne k_j$ deterministically.
  For (2), $k_1, \dots, k_\omega$ are independent under CRT. Let $e_i = v_2(\varphi(p_i^{a_i})) \ge 1$ and $e_{\min} = \min_i e_i$. For each $1 \le c \le e_{\min}$, $\Pr[\forall i, k_i = c] = \prod_i 2^{-(e_i - c + 1)} \le 2^{-\omega(e_{\min} - c + 1)}$, and $\Pr[\forall i, k_i = 0] = \prod_i 2^{-e_i} \le 2^{-\omega e_{\min}}$. Summing gives $\Pr[k_1 = \dots = k_\omega] \le \frac{2^{-\omega} - 2^{-\omega(e_{\min}+1)}}{1 - 2^{-\omega}} + 2^{-\omega e_{\min}} = \frac{1 + (2^\omega - 2)2^{-\omega e_{\min}}}{2^\omega - 1} \le \frac{1 + 2^\omega - 2}{(2^\omega - 1) 2^{\omega - 1}} \dots$ — check directly: for $e_{\min} \ge 1$, $2^{-\omega e_{\min}} \le 2^{-\omega}$, so numerator $\le 1 + (2^\omega - 2)2^{-\omega} = 2 - 2^{1-\omega} = 2^{1-\omega}(2^\omega - 1)$, giving $\Pr[k_1 = \dots = k_\omega] \le 2^{1-\omega}$. $\blacksquare$

### 2.2 Cost Ledger of a Hypothetical Thin-Dependency Algorithm Under H (`Prop THIN-LEDGER`) [model-proved]

Why is R **load-bearing** in `notes/lowerbound_results.md:36`? Suppose an algorithm could genuinely extract a split from $h = L_N[\eta, c_h]$ ($0 \le \eta < 1/3$) $B$-smooth relations on a standard rank-2 GNFS polynomial pair ($\lambda(r) \sim 2\sqrt{2 x \ln r}$), where $B = e^\mu$ is a free parameter!

1. **Thinness by restricting prime support (`Fixed-Column Thinness`, [model-proved]):**
   Suppose $h \ll \pi(B)$ is achieved by requiring every kept relation to use only a fixed set $S$ of $|S| < h$ primes/atoms (for instance, $S = Q_0$ of size $\sim 2\ln^2 N/\ln\ln N$ for pure $Q_0$-PIN). By `notes/smoothbuild_results.md:219` **Lemma AT**, any $|S|$ multiplicatively independent atoms generate at most $\Psi(e^s, p_{|S|})$ integers $\le e^s$. Under the natural subset-Dickman refinement of H, each probe succeeds with probability $\le \rho(s / \ln p_h)^{1-o(1)}$, which is identical to setting $\mu_{\text{eff}} = (1 + o(1))\ln h$. Thus **R holds automatically at $\mu_{\text{eff}}$**, forcing $W \ge L_N[1/3, (4A/3)^{2/3}]$ — and if $h = (\ln N)^{O(1)}$ (e.g. pure $Q_0$-PIN with $B_{\text{eff}} = 2\ln^2 N$, $\mu_{\text{eff}} \sim 2\ln x$), balancing $\ln W \ge \frac{2\sqrt{x \ln W}}{2\ln x}\ln x = \sqrt{x \ln W}$ forces **$\ln W = \Omega(x)$, i.e. exponential time $N^{\Omega(1)}$ ($\alpha = 1$)**!

2. **Thinness at large $B$ (`Unrestricted-Support Thinness`, [model-proved]):**
   Suppose instead a relation family allowed **all** primes $\le B = e^\mu$ ($\mu$ freely chosen large!) and still produced a valid dependency from $h = L_N[\eta, c_h]$ ($0 \le \eta < 1/3$) $B$-smooth relations:
   - *In the pure valuation-probe model* (where 1 smoothness probe costs $1$ regardless of $B$, `notes/lowerbound_results.md:30`): setting $\mu = s$ makes $\rho(s/\mu) = 1$, collapsing probe cost to $O(h) = L_N[\eta]$ (even $O(1)$ if $h = O(1)$!). This proves Axiom R is **indispensable** in the pure probe model.
   - *In real RAM time* (charging ECM cofactoring $T_{\text{ECM}}(\mu) = \exp((1+o(1))\sqrt{2\mu \ln \mu})$ per probe + linear algebra $h^\omega$):
     Optimize probe count $k = \exp(c_1 F_\alpha)$ and smoothness log-bound $\mu$ subject to $k \cdot \rho(A\sqrt{x\ln k}/\mu) \ge h$:
     * **For $0 \le \eta \le 1/5$** (including $h = 1$ or $h = O(\ln^2 N)$): set $\alpha = 1/5$, degree $d \asymp (x/\ln x)^{2/5}$, norm log-size $s = A\sqrt{c_1} x^{3/5}(\ln x)^{2/5}$, and $\mu = c_2 x^{2/5}(\ln x)^{3/5}$ ($B = L_N[2/5, c_2]$). Then:
       - ECM per probe costs $\exp(\sqrt{2 c_2 x^{2/5}(\ln x)^{3/5} \cdot \frac{2}{5}\ln x}) = \exp(\sqrt{\frac{4}{5}c_2} F_{1/5})$;
       - Smoothness penalty is $-\ln \rho(s/\mu) = (1+o(1))\frac{A\sqrt{c_1}}{5 c_2} F_{1/5}$;
       - Yield $k \rho(s/\mu) \ge h$ holds whenever $c_1 > \frac{A\sqrt{c_1}}{5 c_2} \iff c_2 > \frac{A}{5\sqrt{c_1}}$;
       - Minimizing $c_1 + \sqrt{\frac{4}{5} c_2}$ over $c_1 > 0$ gives $c_1 = (A/100)^{2/5}$ and **total RAM time $L_N\!\left[1/5,\ (\sqrt{5} A / 4)^{2/5}\right]$** ($= L_N[1/5, (\sqrt{5}/2)^{2/5} \approx 1.0456]$ for GNFS $A = 2$, and $L_N[1/5, (5/8)^{1/5} \approx 0.9103]$ for $A = \sqrt{2}$)!
     * **For $1/5 < \eta < 1/3$**: set $k = h \cdot L_N[\eta, o(1)]$ ($s \asymp x^{(1+\eta)/2}(\ln x)^{(1-\eta)/2}$) and $\mu \asymp x^{(1-\eta)/2 + \varepsilon}$ ($B = L_N[(1-\eta)/2 + \varepsilon]$). Then $-\ln \rho(s/\mu) \asymp x^{\eta - \varepsilon}\ln x = o(F_\eta)$, while ECM costs $\exp(x^{(1-\eta)/4 + \varepsilon/2}) = o(L_N[\eta])$ (since $(1-\eta)/4 < \eta \iff \eta > 1/5$). Total RAM time is **$L_N[\eta, O(c_h)]$**!
     * **Combined Thin-Dependency RAM Ledger:** **$W_{\text{RAM}}(\eta) = L_N[\max(1/5, \eta)]$** — strictly sub-$L_N[1/3]$ for every $\eta < 1/3$! (Indeed, today's GNFS + ECM ALREADY generates $L_N[1/4]$ completely factored $L_N[1/2]$-smooth relations in $L_N[1/4]$ RAM time; it fails to factor $N$ **solely** because those $L_N[1/4]$ rows have rank $L_N[1/4]$ inside $\pi(L_N[1/2])$ columns!)

### 2.3 Why None of the Candidate Sub-Mechanisms Breaks R (Four Obstructions)

1. **Z-route PIN (`Lemma PIN-STRIP`, [model-proved]):**
   In `notes/zroute_results.md:211, 223`, PIN requires vectors supported *purely* on $\{-1\} \cup Q_0$ ($Q_0 = \{p < 2\ln^2 N\}$): all active columns in $\mathcal{E} = (2\ln^2 N, B]$ must still cancel over $\mathbb{Q}$ (`rank-nullity`: $\dim W_{Q_0} = \operatorname{rank}(R) - \operatorname{rank}(R_{\mathcal{E}}) \le h - \operatorname{rank}(R_{\mathcal{E}})$). Could $Q_0$ absorb a non-negligible share of candidate norm sizes?
   *Claim (`Lemma PIN-STRIP`):* Under H's partial-divisibility bound (and unconditionally for uniform/rational cofactors by Rankin's method at saddle $\alpha = 11/12$: $\sum_{m \ge e^\sigma, P^+(m) \le 2x^2} m^{-1} \le e^{-\frac{11}{12}\sigma} \prod_{p \le 2x^2}(1 - p^{-11/12})^{-1} = \exp(-\frac{11}{12}\sigma + O(x^{1/6}/\ln x))$), across any $t = \lfloor e^{c F_{1/3}} \rfloor$ probes, with probability $1 - o(1)$ **every** probed cofactor $C_r$ has $Q_0$-smooth part $\le \exp(2 c F_{1/3}) = \exp(o(x^{2/3}(\ln x)^{1/3}))$.
   Thus stripping $Q_0$ for free shrinks $\lambda(r)$ by $o(\lambda(r))$ and leaves every row needing exact $\mathbb{Z}$-cancellation on $(2\ln^2 N, B]$ (which is strictly harder than $\mathbb{F}_2$ parity cancellation whenever even prime powers $p^2 \parallel C_r$ occur!).

2. **Large-prime cycles ($k$-LP hypergraphs) [known + model-proved]:**
   For unforced cofactors whose primes $> B_0$ obey conditional marginal incidence $\le q^{-1+o(1)}$, Croot–Granville–Pemantle–Tetali (2012) and Balister–Bollobás–Morris (2018) (`notes/smoothbuild_results.md:236`) prove the 2-core is empty below $(e^{-\gamma} - o(1)) \min_\mu \pi(e^\mu)/\rho(s/\mu)$, recovering $\ln t \ge (1-o(1))\sqrt{2 s \ln s}$ even when **all** primes ($B = \infty$, unbounded $k$-LP cycles) are allowed. Polynomial root-count second moments (`notes/smoothbuild_results.md:251` Prop GR) inflate pair collisions by at most $r_G - 1 \le d - 1 = L_N^{o(1)}$.

3. **Quadratic characters / Schirokauer maps (`Lemma CHAR-TRAP`, [proved]):**
   In $\mathcal{O}_K$, quadratic characters $(\gamma / \mathfrak{q})$ and $\ell$-adic Schirokauer maps only detect obstructions inside the finite-rank quotient $H_{K,\ell} = \{x \in K^\times : (x) = \mathfrak{a}^\ell\}/K^{\times \ell}$ (`notes/rkrank_results.md:125`, `notes/zroute_results.md:207`). Without cancelling prime-ideal valuations $v_{\mathfrak{p}}(\gamma) \equiv 0 \pmod \ell$ across all $\mathfrak{p}$ first, local power residues mod $\mathfrak{q}_1, \dots, \mathfrak{q}_m$ do not lift to a global power in $K^\times$ (`notes/fruit_results.md:350`). Thus they *add* $O_d(\log N)$ columns, never *replace* prime-ideal columns.

4. **Planted dependencies (`Lemma PLANT`, [proved + model-proved]):**
   Suppose $h$ algebraic elements $\gamma_1(X), \dots, \gamma_h(X) \in \mathbb{Z}[X]$ have a planted square product $\prod_{i=1}^h \gamma_i(\alpha) = \Gamma(\alpha)^2$ in $K = \mathbb{Q}[X]/(f(X))$ ($f(m) \equiv 0 \pmod N$). Writing $\prod_{i=1}^h \gamma_i(X) - \Gamma(X)^2 = Q(X) f(X)$ in $\mathbb{Z}[X]$:
   - If $Q(X) = 0$, then $\prod \gamma_i(m) = \Gamma(m)^2$ holds identically over $\mathbb{Z}$, yielding only tautological/fruitless congruences $U = \pm \Gamma(m)$ (`notes/smoothbuild_results.md:180`).
   - If $Q(X) \ne 0$ (e.g. $h = 1$, $d = 2k$, $\Gamma(X) = X^k + \dots$, $R_\Gamma(X) = \Gamma(X)^2 - f(X)$), evaluating at $X = m$ gives $R_\Gamma(m) = \Gamma(m)^2 - N$, collapsing to MPQS ($\ln|R_\Gamma(m)| \ge \frac{1}{2}x$, $\alpha = 1/2$).
   - Even granting an oracle that plants **free** algebraic squares for every rational linear form $(a - bm)$ (`notes/smoothenum_results.md:126` Table 1 row T6), cancelling the rational primes of $a - bm \approx E N^{1/d}$ still requires $\pi(B_{\text{rat}})$ rows and costs $\ge L_N[1/3, (4/9)^{1/3} = 0.7631]$ — still $\alpha = 1/3$!

---

## 3. Attack on (b): Shared Descent DAGs (`Theorem DAG-POOL`)

### 3.1 Why Assumption T Was Thought Load-Bearing
In `notes/heurnfs_results.md:188-201`, the descent bottom bound `#bottom >= #net_rows` used Assumption T (`deg(q) <= 2`). When one eliminator $e_q$ serves $m \gg 2$ parents, two cases occur in the peeled core $(\mathcal{R}', \mathcal{L}')$ above base bound $B_0$:
- **Star of 1-LP rows:** If $e_q$ and its $m$ parents $p_1, \dots, p_m$ each contain *only* $q$ above $B_0$, then every parent $p_i$ is **itself** a 1-LP bottom relation (`R_1 = m + 1 > m = net_rows`), so `#bottom >= #net_rows` already holds without Assumption T!
- **Multi-prime internal sharing / special-$q$ pools:** In general, eliminating $|\mathcal{L}'|$ independent columns above $B_0$ to leave $h_0$ independent rows over $B_0$ still requires $|\mathcal{R}'| \ge h_0 + |\mathcal{L}'|$ distinct fresh relations, and under SM2 above $B_* = W^2(\ln N)^C$ at most one hit per prime $> B_*$ is unforced. However, in `codex/lowerbound_1_codex_proof.md:73`, pooling $s$ distinct sublattices lost $-\frac{D}{2}\ln s$, which vanishes once $s = L_N[1/3]$!

### 3.2 `Theorem DAG-POOL` (Sublattice Discriminant Cancellation Closes Heavy Sharing) [proved in model]

**Statement (`Theorem DAG-POOL`).** Let $W(X, Y) \in \mathbb{Z}[X, Y]$ be a primitive squarefree binary form of degree $D \ge 4$ with $|\operatorname{disc}(W)| \ge N^k$ ($k \in \{1, 2\}$). Let $B_0 = e^{\mu_0}$ with $\mu_0 \ge \beta_0 F_{1/3}$ ($\beta_0 > 0$). Let $\Lambda_1, \dots, \Lambda_s \subset \mathbb{Z}^2$ ($1 \le s \le W_{\text{work}} \le e^{c F_{1/3}}$) be **any** $s$ distinct primitive-containing divisor sublattices of squarefree indices $Q_1, \dots, Q_s$ satisfying $P^-(Q_j) > B_0$ and $Q_j \mid W(\Lambda_j)$, with stripped integer forms $G_j(u, v) := W(M_j(u, v)) / Q_j$ ($\det M_j = Q_j$). Then:
1. **[proved, exact algebra]** $|\operatorname{disc}(G_j)| = Q_j^{(D - 1)(D - 2)} |\operatorname{disc}(W)| \ge N^k Q_j^{(D - 1)(D - 2)}$.
2. **[proved, sublevel sum]** Across **all $s$ sublattices pooled together** (plus $\Lambda_0 = \mathbb{Z}^2$), after discarding at most $(s + 1) C_{\text{FT}} D \ln(2E)$ Farey-cusp points, the $r$-th smallest stripped residual log-size satisfies:
   $$\lambda_{\text{pool}}(r) \ge \frac{k x}{2(D - 1)} + \frac{D}{2}\ln r - O\!\left(D \ln D + \frac{D \ln D}{\mu_0}\ln s + \ln s\right) = (1 - o(1))\left[\frac{k x}{2(D - 1)} + \frac{D}{2}\ln r\right]$$
   for all relevant ranks $r$, because $D \ln D, \ln s = O(F_{1/3}) = o(x^{2/3}(\ln x)^{1/3})$!
3. **[model-proved consequence]**
   - **Heavily shared DAGs ($s \le \frac{h_0}{8 C_{\text{FT}} D \ln(2E)}$):** When $s$ shared eliminators/sublattices serve $h_0 \ge B_0^{1-o(1)}$ net $B_0$-smooth cofactor probes, total cusp exceptions are $\le h_0 / 8 < \lfloor h_0 / 4 \rfloor$, so **neither Assumption T nor `H_cusp` is needed**: C's supermartingale proves the full $L_N[1/3, (16k/9)^{1/3}]$ floor ($1.2114$ at $k = 1$; $1.5263$ at $k = 2$).
   - **Arbitrary $s \le W_{\text{work}}$:** Under `H_cusp` (bounding Farey-cusp hits across the $s$ sublattices by $\le h_0/4$), the exact same $1.2114 / 1.5263$ floor holds — upgrading Round 2's $(k/18)^{1/3}$ ($0.3816 / 0.4807$) whenever a $B_0$-quota $h_0$ is reached! And if **no** $B_0$-quota is used ($B_0 = O(x^2)$), either all primes lie $\le B_* = W^2 x^C$ (falling under R at scale $\le B_*$) or some prime $> B_*$ occurs (where Round 2's `H_cusp` single-witness bound $(k/18)^{1/3}$ already forces $\alpha \ge 1/3$).

*Proof of `Theorem DAG-POOL`.*
(1) $\operatorname{disc}(W \circ M_j) = (\det M_j)^{D(D-1)} \operatorname{disc}(W) = Q_j^{D(D-1)} \operatorname{disc}(W)$. Since $\operatorname{disc}$ is homogeneous of degree $2D - 2$ in coefficients, dividing by $Q_j$ multiplies by $Q_j^{-(2D-2)}$, leaving exponent $D(D-1) - (2D-2) = (D-1)(D-2)$.
(2) By Bean's inequality + `Lemma FT` (`codex/lowerbound_1_codex_proof.md:69-82`), sublattice $j$ contributes $\le C_0 |\operatorname{disc}(G_j)|^{-\frac{1}{D(D-1)}} e^{2\sigma/D} + C_{\text{FT}} D \ln(2E)$ primitive box points with $0 < |G_j| \le e^\sigma$, where $|\operatorname{disc}(G_j)|^{-\frac{1}{D(D-1)}} \le N^{-\frac{k}{D(D-1)}} Q_j^{-\frac{D-2}{D}}$.
Since $W$ is primitive of degree $D$, mod each prime $p \mid Q$ ($Q$ squarefree) $W \bmod p \ne 0$ vanishes on $\le D$ lines of $\mathbb{P}^1(\mathbb{F}_p)$; by CRT at most $D^{\omega(Q)}$ primitive-containing index-$Q$ sublattices satisfy $Q \mid W(\Lambda)$. Since $P^-(Q) > B_0$, $\omega(Q) \le \frac{\ln Q}{\mu_0}$, so $D^{\omega(Q)} \le Q^{\delta_x}$ with $\delta_x := \frac{\ln D}{\mu_0}$.
Split the $s$ sublattices at $Q_j \le s^2$ vs $Q_j > s^2$:
- For $Q_j > s^2$ ($D \ge 4 \implies \frac{D-2}{D} \ge \frac{1}{2}$): $\sum_{Q_j > s^2} Q_j^{-\frac{D-2}{D}} \le s \cdot (s^2)^{-1/2} = 1$.
- For $Q_j \le s^2$: each integer $n \le s^2$ occurs at most $n^{\delta_x} \le s^{2\delta_x}$ times, giving $\sum_{Q_j \le s^2} Q_j^{-(1 - 2/D)} \le s^{2\delta_x} \sum_{n=1}^{s^2} n^{-1 + 2/D} \le \frac{D}{2} s^{2\delta_x + 4/D}$.
Inverting $r \le C_0 N^{-\frac{k}{D(D-1)}} e^{2\sigma/D} (2 + \frac{D}{2} s^{2\delta_x + 4/D})$ gives $\sigma \ge \frac{kx}{2(D-1)} + \frac{D}{2}\ln r - O(D\ln D + D\delta_x \ln s + \ln s)$. At $D \asymp (x/\ln x)^{1/3}$ and $\mu_0 = \beta_0 F_{1/3}$, $D\delta_x \ln s = \frac{D\ln D}{\beta_0 F_{1/3}}\ln s \le \frac{c}{\beta_0} D\ln D = O(F_{1/3}) = o(x^{2/3}(\ln x)^{1/3})$. $\blacksquare$

---

## 4. Attack on (c): Reuse Across Related Inputs for a Single Given $N$

### 4.1 Complete Ledger of Auxiliary-Instance & Witness-Reuse Schemes for One $N$

**Crucial distinction (`Theorem SINGLE-N-REUSE`, [model-proved]):** In `notes/lowerbound_results.md:40`, the reuse budget $\sum_{i=1}^K Q_i \le e^{o(F)} W$ is needed **only** when amortizing across $K \gg 1$ distinct input targets (where dense-batch Adversary A1 gets $\operatorname{poly}(\ln N)$ per target at $K \approx N/\operatorname{poly}(\ln N)$, `codex/lowerbound_1_gemini.out:59`). When given **one** input $N$ ($K_{\text{input}} = 1$), manufacturing $K$ auxiliary instances $M_1, \dots, M_K$ does *not* divide total work $W$ by $K$: every probe is charged once to $N$ ($Q_N \le W$, reuse $S = 1$), and all auxiliary presentations simply form a pooled or adaptive menu of blocks for $N$!

| # | Related-Input / Reuse Mechanism for Single $N$ | Shared Witness / Base | Single-$N$ Cost Ledger | Status & Obstruction |
|---|---|---|---|---|
| C1 | Multipliers $k N$ in QS / MPQS ($x^2 - kN$) | Shared rational FB $\le B$ | $L_N[1/2, 1]$ | **[known/model-proved]** Residue size $\ge 2E\sqrt{kN} \ge 2E\sqrt{N}$ (`notes/heurnfs_results.md:75`). |
| C2 | Multipliers $k N$ in NFS, shared base $m$ ($f_k(m) = kN$) | 100% shared rational sieve $a-bm$ & rational FB; disjoint algebraic FBs | $L_N[1/3, 1.9019]$ ($\omega=2$); $1.8595$ ($\omega=0$) | **[model-proved]** Identical to Coppersmith MNFS (`notes/smoothenum_results.md:329`); `escapeR2` Lemma TI blocks cross-$k$ ideal cancellation when $K_{k_1}\cap K_{k_2}=\mathbb{Q}$. |
| C3 | Shifted bases $m_k = m_0 + k$ or universal factory $m_0 = 2^t$ | $\le 75\%$ rational overlap (`m_0+k`) or 1 chargeable table (`2^t`) | $\ge L_N[1/3, 1.9019]$ ($2.0069$ at factory pt) | **[model-proved]** Table precomputation $2e$ is charged to the single $N$ (`notes/smoothenum_results.md:322-345`). |
| C4 | Multiplier search $1 \le k \le K$ for SNFS-small height $H_k$ (`Lemma SNFS-MULT`) | One low-height pair $(f_{k^*}, g_{k^*})$ for $k^* N$ | $\ge L_N[1/3, 1.9230]$ under G ($1.5263$ Bean/FT floor) | **[proved under G / Bean]** Pigeonhole $\ln K + d\ln H \ge x$ inflates $\ln m_k \approx \frac{x+\ln K}{d}$, conserving $\frac{\ln(kN)}{d}+\ln H_k \ge \frac{2x}{d}$! |
| C5 | Powers $N^j$ ($j \ge 2$) / $N$-adic lifts $\bmod N^j$ | One relation $\bmod N^j$ | $L_N[1/3, j^{1/3} c]$ ($> c$) | **[model-proved]** $N^j \mid \operatorname{Res}$ scales norm $\frac{2jx}{D}$; kernel $(\mathbb{Z}/N^j)^\times \to (\mathbb{Z}/N)^\times$ has rough order $N^{j-1}$. |
| C6 | Product bundling $M = M_1 \cdots M_K$ (1 witness splits $K$ coprime moduli) | $100\%$ probe reuse ($S = K$) | $\frac{1}{K} L_N[1/3, K^{1/3} c] \to \infty$ | **[model-proved]** $K$ enters *inside* $\exp(c K^{1/3} F_{1/3})$ via $\ln M = K x$, blowing up faster than $1/K$. |
| C7 | Dense-batch difference factory A1 on $\{k N : 1 \le k \le K\}$ | Smooth pairs $u - v = k N$ ($u \equiv v \bmod N$) | $N^{1/2+o(1)}$ ($\alpha = 1$) | **[proved]** Multiples of $N$ in $[1, KN]$ have scale-invariant density $\frac{K}{KN} = \frac{1}{N}$ (birthday floor $\sqrt{N}$). |
| C8 | Polynomial/ring maps $M_k = h_k(N)$ ($\Phi_m(N)$, $N \pm 1$) | Any batch method on $\{h_k(N)\}$ | Carries $0$ split of $N$ | **[proved]** `notes/selfred_results.md:127` Thm SR1: divisors of $h_k(pq)$ depend only on $pq = N$, moving no large factor of $p$. |
| C9 | Exact algebraic moved-divisor $M_k = A_k(p)B_k(q)$ | Any solver on $\{M_k\}$ | Poly-time *if constructible* | **[proved]** `notes/selfred_results.md:226` Thm SR2: constructing one non-monomial $M_k$ already factors $N$ via Sturm/LLL. |
| C10 | **Best Candidate (`Candidate C-RM`):** $F$-relative Reed–Muller tied tower $L = F(\sqrt{\delta_1}, \dots, \sqrt{\delta_r})$ over base-$m$ field $F = \mathbb{Q}(\theta)$ ($f_0(m)=N$) | Shares rational $m$ **AND** compresses $V = 2^r - 1$ extensions into $(r+1)\pi_F(B_a)$ ideal-parity columns! | Unconditional $\delta_i$: $\alpha = 1/2$; under **ML-RM** ($\theta_{\text{norm}} \to 0$): constant improvement at $r=O(1)$, **$\alpha = 1/4$ only if multi-coordinate norm-sharing + unit control hold** | **[new theorem + open ML-RM]** See §4.2–4.3 below! |

### 4.2 New Theorem on Tied-Field Column Sharing (`Theorem TIED-PARITY`) [proved]

`notes/escapeR2_results.md:211-234` and `notes/completeness_audit_results.md:125` left open whether tied subfields in an abelian Galois tower can share prime-ideal parity columns across $V$ fields to evade MNFS's $V \cdot B_a$ column growth. We settle the exact column rank for all elementary abelian $p_0$-towers:

**Theorem TIED-PARITY [proved].** Let $p_0$ be prime, $r \ge 2$, and $L/F$ a Galois extension with $\operatorname{Gal}(L/F) \cong (\mathbb{Z}/p_0\mathbb{Z})^r$, containing $V = \frac{p_0^r - 1}{p_0 - 1}$ relative degree-$p_0$ subfields $K_w = L^{\ker w}$ ($[w] \in \mathbb{P}((\mathbb{F}_{p_0}^r)^*)$). Let $\mathfrak{p}$ be any prime of $F$ unramified and splitting completely in $L$, and let $\mathcal{S}_{p_0, r}(k) \subseteq k^{\mathbb{F}_{p_0}^r}$ be the $k$-span of the pullback valuation vectors of all $V \cdot p_0$ prime ideals of all $K_w$ above $\mathfrak{p}$.
1. **Over $k = \mathbb{Q}$ (for every prime $p_0 \ge 2$) and over $k = \mathbb{F}_2$ (for every odd prime $p_0 \ge 3$):**
   $\dim_k \mathcal{S}_{p_0, r}(k) = p_0^r = 1 + V(p_0 - 1)$. **Zero column compression ($\tau = 0$) occurs across the $V$ subfields!**
2. **Over $k = \mathbb{F}_2$ at $p_0 = 2$ (multiquadratic $L = F(\sqrt{\delta_1}, \dots, \sqrt{\delta_r})$, $V = 2^r - 1$):**
   $\mathcal{S}_{2, r}(\mathbb{F}_2) = \mathrm{RM}(1, r)$ (the first-order Reed–Muller code), so **$\dim_{\mathbb{F}_2} \mathcal{S}_{2, r}(\mathbb{F}_2) = r + 1$** — compressing $2^r$ columns into $r + 1$ ($\tau = 2^r - r - 1$ cross-field parity transfers per split prime)!

*Proof.* Identify the primes of $L$ above $\mathfrak{p}$ with $\mathbb{F}_{p_0}^r$; the primes of $K_w$ pull back to the affine hyperplane indicators $1_{\{\langle w, x \rangle = c\}}$. First, $\mathbf{1} = \sum_{c \in \mathbb{F}_{p_0}} 1_{\{\langle w, x \rangle = c\}} \in \mathcal{S}_{p_0, r}(k)$. For any $x_0 \in \mathbb{F}_{p_0}^r$, exactly $V = \frac{p_0^r - 1}{p_0 - 1}$ hyperplanes pass through $x_0$, and for $y \ne x_0$ exactly $V' = \frac{p_0^{r-1} - 1}{p_0 - 1}$ pass through both $x_0$ and $y$. Hence $\sum_{H \ni x_0} 1_H = V \delta_{x_0} + V'(\mathbf{1} - \delta_{x_0}) = p_0^{r-1} \delta_{x_0} + V' \mathbf{1}$. Since $p_0^{r-1}$ is invertible in $\mathbb{Q}$ (all $p_0$) and $p_0^{r-1} \equiv 1 \pmod 2$ (all odd $p_0 \ge 3$), every point indicator $\delta_{x_0}$ lies in $\mathcal{S}_{p_0, r}(k)$, proving (1). For $p_0 = 2$ over $\mathbb{F}_2$, $1_{\{\langle w, x \rangle = c\}} = 1 + c + \sum_{i=1}^r w_i x_i$, which spans $\mathrm{RM}(1, r)$ of dimension $r + 1$, proving (2). $\blacksquare$

### 4.3 Best Candidate (`Candidate C-RM`), Exponent Ledger, and Missing Lemma (`ML-RM`)

By `Theorem TIED-PARITY`, growing odd-prime degrees ($p_0 \to \infty$) have $\tau = 0$, while $p_0 = 2$ over $F = \mathbb{Q}$ fixes degree $2$ ($\alpha = 1/2$). The **unique** algebraic geometry combining growing degree $d$ with exponential $\mathrm{RM}(1, r)$ column compression is **`Candidate C-RM`**:
- Choose base-$m$ polynomial $f_0(X)$ of degree $d$ ($f_0(m) = N$, $F = \mathbb{Q}(\theta)$, $\psi_0(\theta) = m \bmod N$).
- Pick $r$ elements $\delta_1, \dots, \delta_r \in \mathbb{Z}[\theta]$ independent in $F^\times / F^{\times 2}$ with **known modular square roots** $\psi_0(\delta_i) \equiv u_i^2 \pmod N$.
- Form the $V = 2^r - 1$ quadratic extensions $K_w = F(\sqrt{\Delta_w})$ ($\Delta_w = \prod_{i: w_i=1} \delta_i$), compatible with $\Psi_L(\sqrt{\Delta_w}) = u_w := \prod_{i: w_i=1} u_i \bmod N$. By `Theorem TIED-PARITY`, all $V = 2^r - 1$ fields share a **single** $(r + 1)\pi_F(B_a)$-column ideal-parity base (plus $\le 2^r d$ unit/character columns of $H_L$).

**Why our unconditional instantiation hits $\alpha = 1/2$:**
Setting $u_i = \lfloor\sqrt{m}\rfloor + i$ and $\delta_i = (u_i^2 - m) + \theta$ gives $\psi_0(\delta_i) = u_i^2$ in $O(1)$ time with $|u_i| \approx N^{1/(2d)}$, **but** $|a_i| \approx N^{1/(2d)}$ makes $|N_{F/\mathbb{Q}}(\delta_i)| = |f_0(-a_i)| \approx m |a_i|^{d-1} \approx N^{\frac{1}{2} + \frac{1}{2d}}$, so $|N_{F/\mathbb{Q}}(\Delta_w)|$ grows as $N^{|w|(\frac{1}{2} + \frac{1}{2d})}$! Even optimizing linear $b_i \theta + a_i$ with $b_i m + a_i = u_i^2$ ($|a_i|, |b_i| \approx m^{1/3}$) leaves $|N_{F/\mathbb{Q}}(\delta_i)| \gtrsim m (m^{1/3})^d = N^{1/3 + 1/d}$ per factor.

**Exact Missing Lemma (`ML-RM`, [open; conjectured false for generic $N$]):**
Given $N$, construct in time $L_N[\alpha]$ ($\alpha < 1/3$) a degree-$d$ field $F = \mathbb{Q}(\theta)$ ($\psi_0 : \mathbb{Z}[\theta] \to \mathbb{Z}/N\mathbb{Z}$) and $r$ $F^\times/F^{\times 2}$-independent $\delta_i \in \mathcal{O}_F$ with $\psi_0(\delta_i) \equiv u_i^2 \pmod N$ such that for a family of $M = L_N[\alpha]$ elements $\gamma_w \in K_w$:
1. Both $|\psi_{K_w}(\gamma_w) \bmod N|$ and $|N_{K_w/\mathbb{Q}}(\gamma_w)|^{1/2}$ have log-size $o(x^{2/3}(\ln x)^{1/3})$ (specifically $\lesssim x^{2\alpha}$ at $d \asymp x^{1-2\alpha}$);
2. The unit/Selmer obstruction $\dim_{\mathbb{F}_2} \operatorname{span}\{\gamma_w\} \bmod \text{squares}$ inside $H_L$ is $\le L_N[\alpha]$ (requiring $2^r \le L_N[\alpha]$, i.e. $r = O(F_\alpha)$ — note at $2^r = L_N[\alpha]$, $\mathrm{RM}(1, r)$ compresses $L_N[\alpha]$ subfields into $O(F_\alpha \cdot \pi(B_a))$ columns!).
Wait: notice that condition (1) requires simultaneously small rational lift and small algebraic norm across an $r$-dimensional space — reducing directly back to Hatch 2's **`LP_exc^rel`** (`notes/subthird_results.md:83`)!

---

## 5. Pre-Registered Toy Experiments (`prereg.md` SHA-256 `ddaa0200...`) [measured]

Code & raw data in `/Volumes/git/factoring-connection/codex/scratch/hatch_h3g/`:
- Pre-registration: `prereg.md` (written & hashed before any experimental script existed)
- Operational + harness runner: `run_hatch3_experiments.py` → `hatch3_results.json` (`63.24 s`, exact integer / `Fraction` / `GF(2)` arithmetic, zero factor leakage in operational routines)
- Post-hoc diagnostic: `posthoc_e2.py` → `posthoc_e2.json`

### Summary of Registered Predictions (7 / 8 PASS, 1 / 8 FAIL — diagnosed)

| ID | Tested Statement | Registered Criterion | Measured Result | Verdict |
|---|---|---|---|---|
| **P1a** | `Lemma PIN-1`: 1-coordinate $\mathbb{Z}$-route Miller split ($B=600, \pi(B)=109$, 40-bit $N$) | Blum Jacobi-$-1$ $q_0$: `10/10`; general random unit $b$: `[45, 85]/100` | Blum **`10/10`**; general **`55/100`** | **PASS** |
| **P1b** | Thin-row obstruction: rows $h_{Z1}$ to pin 1 coordinate $q_0$ vs $h_{F2}$ | $h_{Z1} > 54$ in `20/20`; $h_{Z1}/h_{F2} \ge 0.80$ in $\ge 18/20$ | $\min h_{Z1} = \mathbf{101}$ (`20/20`); median $106$ vs $100$; $h_{Z1}/h_{F2} \ge \mathbf{1.019}$ in **`20/20`** | **PASS** |
| **P2a** | `Theorem DAG-POOL` (1): exact sublattice discriminant scaling ($D=4$, 60-bit $N$) | $\operatorname{disc}(G_j) = q_j^6 \operatorname{disc}(W)$ in `40/40` | **`40/40`** exact integer matches (`0` violations) | **PASS** |
| **P2b** | `Theorem DAG-POOL` (2): raw pooled $\lambda(200; s=64) - \lambda(200; s=1)$ | $\ge -2.50$ nats in `10/10` (vs naive $-8.318$) | **`-15.09` to `-11.17`** (mean `-12.28`) | **FAIL** *(diagnosed below)* |
| **P2c** | Shared-DAG rank conservation without Assumption T (`m ∈ {2,5,20}`, stars & 2-level diamonds) | `net_B0 <= fresh` and (under SM2) `B0_cofactor >= net_B0` in `30/30` | **`30/30`** and **`30/30`** ($\min \text{ratio} = 1.000$) | **PASS** |
| **P3a** | `Theorem TIED-PARITY`: exact $\mathbb{F}_2$ & $\mathbb{Q}$ ideal ranks in $(\mathbb{Z}/p_0\mathbb{Z})^r$ | $p_0=2 \implies \operatorname{rk}_{\mathbb{F}_2} = r+1$ (`4/4`); odd $p_0 \in \{3,5\}$ & all $\mathbb{Q} \implies \operatorname{rk} = p_0^r$ (`10/10`) | $p_0=2, r=2..5$: $\operatorname{rk}_{\mathbb{F}_2} = \mathbf{3,4,5,6}$ ($\tau = 1,4,11,26$); odd $p_0$ & $\mathbb{Q}$: **`100%` full rank ($\tau = 0$)** | **PASS** |
| **P3b** | `Candidate C-RM` ($d=3$, $r=4$, $V=15$ extensions, 48-bit $N$): square congruence & norm slope | `150/150` exact $\psi_0(\Delta_w) \equiv u_w^2 \pmod N$; slope $\in [0.45, 0.75]$ in `10/10` | **`150/150`** exact; slopes **`0.571`–`0.639`** (mean `0.610` vs asymp. `2/3 = 0.667`) | **PASS** |
| **P3c** | `Lemma SNFS-MULT` ($d=3$, $1 \le k \le 1024$, 48-bit $N$): $(\ln m_k + \ln H_k)/(\frac{2}{3}\ln N)$ | $\ge 0.85$ in `10/10` | $\min$ across $10,240$ pairs = **`0.9007`** (`0.901`–`0.993` across `10/10`) | **PASS** |

### Post-Hoc Diagnosis of `P2b` (`posthoc_e2.py`, labelled post-hoc; registered failure stands)
Why did raw $\lambda_{\text{pool}}(200; s=64) - \lambda(200; s=1)$ equal $-12.28$ nats? Two distinct effects compounded in our mis-calibrated baseline:
1. **`s=1` at $q_1 > 500$ vs $Q=1$ base box $\Lambda_0$:** At $q_1 \approx 503$, the single sublattice $G_1$ enjoys a $+\frac{D-2}{2}\ln q_1 = +\ln 503 \approx +6.22$ nat size *increase* above the unforced $Q=1$ form $W$ (`raw_r200_s1` $\approx 40.5$ vs `Q1_base_r200` $\approx 34.1$). Pooling $s = 64 \ll B_0 = 500$ primes in $[503, 900]$ trades $+\ln 500$ against $-2\ln 64$, landing near the $Q=1$ floor rather than staying $+6.2$ nats above it!
2. **Small-quota cusp overflow (`r/s = 200/64 = 3.1` points/sublattice):** Taking the top $200$ points across $64$ sublattices selects only the $3$ most skewed cusp vectors per sublattice (directly illustrating `Theorem DAG-POOL` (3)'s cusp condition $r \gg s D \ln(2E)$ and Round 2's small-$h$ cusp overflow `codex/lowerbound_2_gemini.out:5`), dropping another $\approx 3.9$ nats (`raw_r200_s64` $\approx 28.0$ vs 8-cusp-stripped `ft8_r200_s64` $\approx 31.9$).
3. **Direct check of `Theorem DAG-POOL`'s proved inequality:** Against `Q1_base_r200` with 8 cusp points stripped per lattice, `ft8_r200_s64 - Q1_base_r200` is **`-2.16` to `-2.73` nats** (well inside $-\ln 64 = -4.16$, beating naive $-8.32$ by `~5.8` nats), and against the proved Bean floor $\frac{x}{3} + 2\ln(1000)$ at $r = 1000$, `excess_over_bean_r1000_s64` is **`+4.69` to `+6.03` nats strictly positive** across every instance.

---

## 6. Overall Verdict on Hatch 3

- **Did we break $\alpha = 1/3$?** **No.** Every candidate mechanism across (a), (b), and (c) either stays at $\alpha \ge 1/3$ or reduces back to Hatch 1 (`beyond-Johnson CRT decoding M2`) or Hatch 2 (`LP_exc^rel`).
- **Did we make a dent on the lower bound's Hatch 3?** **Yes — four concrete mathematical dents (closing large parts of Hatch 3 and sharpening the remainder):**
  1. **`Lemma PIN-1` + `Lemma PIN-STRIP` (closes the PIN sub-hatch of R):** Although randomized $\mathbb{Z}$-route splitting needs only **1** pinned coordinate (`Lemma PIN-1`, proved & verified: `10/10` Blum, `55/100` general), PIN still requires exact $\mathbb{Q}$-cancellation of all primes $> 2\ln^2 N$. By Rankin/H (`Lemma PIN-STRIP`), stripping $Q_0$ changes norm log-sizes by $O(F_{1/3}) = o(\lambda(r))$, leaving $h_{Z1} > h_{F2}$ (`min ratio 1.019` in E1). Meanwhile, `Prop THIN-LEDGER` proves R is genuinely load-bearing: a hypothetical family needing $h = L_N[\eta]$ ($0 \le \eta < 1/3$) $B$-smooth GNFS rows would run in **$L_N[\max(1/5, \eta)]$ RAM time** via ECM at $B = L_N[\max(2/5, (1-\eta)/2+\varepsilon)]$.
  2. **`Theorem DAG-POOL` (closes heavily shared descent DAGs without Assumption T or `H_cusp`):** Sublattice discriminant growth $Q_j^{(D-1)(D-2)}$ collapses the $s$-sublattice pooling penalty from $-\frac{D}{2}\ln s$ (`codex/lowerbound_1_codex_proof.md:73`) to $-O(F_{1/3}) = o(\lambda(r))$ across *all* $s \le W$ squarefree-index sublattices above $B_0$. Consequently, heavy sharing ($s \ll h_0 / (D\ln E)$) **eliminates both Assumption T and `H_cusp`**, restoring the full $1.2114 / 1.5263$ floor; only *unshared* cusp-hunting ($s \gtrsim h_0 / (D\ln E)$) still needs `H_cusp`.
  3. **`Theorem SINGLE-N-REUSE` + `Lemma SNFS-MULT` (closes single-$N$ auxiliary amortization):** For $K_{\text{input}} = 1$, probe reuse $S = 1$ holds automatically regardless of auxiliary fan-out, and multiplier search $k N$ conserves $\frac{\ln(kN)}{d} + \ln H_k \ge \frac{2x}{d}$ (`min ratio 0.901` over $10,240$ pairs in E3b).
  4. **`Theorem TIED-PARITY` + `Candidate C-RM` (narrows `escapeR2`'s open tied-field problem):** Odd-prime abelian towers ($p_0 \ge 3$) have **zero** $\mathbb{F}_2$ or $\mathbb{Q}$ ideal-column sharing ($\tau = 0$, proved & verified `7/7`), while $2$-elementary towers compress $2^r - 1$ extensions to $r + 1 = \dim \mathrm{RM}(1, r)$ columns mod $2$ (`Candidate C-RM`), whose norm blow-up $N^{|w|(1/2 + 1/(2d))}$ isolates `ML-RM` $\implies$ `LP_exc^rel`.
- **Single sharpest surviving crack in Hatch 3 & Next Step:**
  1. *Proving Axiom R (2-core emptiness below $\pi(B)$) unconditionally or from polynomial-value correlations* — specifically ruling out structured polynomial families with sub-$B$ $\mathbb{F}_2$-cycles (`Prop THIN-LEDGER` shows any $h = L_N[\eta < 1/3]$ break immediately gives $L_N[\max(1/5, \eta)]$).
  2. *Removing `H_cusp` in the unshared regime $s \gtrsim h_0 / (D\ln E)$* via uniform Thue–Mahler / subspace bounds on distinct continued-fraction cusps across $s$ prime sublattices of a single degree-$D$ binary form.
