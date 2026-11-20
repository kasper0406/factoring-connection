# Independent Experimental & Adversarial Specification (`guide_spec.md`)
Role: DESIGNER + ADVERSARY (No Heavy Compute)
Context: `/Volumes/git/factoring-connection/codex/guide_common.md`, `notes/lowerbound_results.md`, `notes/hatch_results.md`, `notes/chase_results.md`, `notes/cadobench_results.md`.

## 0. Core Question
Under the priced valuation-probe model (`notes/lowerbound_results.md`), any adaptive factoring strategy obeys the $L_N[1/3]$ lower bound provided **Assumption H** holds: conditional on all prior probes and internal state, a fresh residual cofactor of log-size $s$ is $B$-smooth with probability $\le \rho(s / \log B)^{1 - o(1)}$. Sieving/divisor tricks up to $y$ trade 1-for-1 (Theorem J, weighted Theorem J in `notes/hatch_results.md`), and norm size is governed by lattice geometry (Theorem KH). This specification tests the lone remaining empirical loophole in H: **does any global, constructive, or history-dependent feature of $(a,b)$ predict $B$-smoothness of the $y$-rough cofactor beyond exact cofactor log-size and exact local congruence information up to $y$?**

---

## 1. Exact Target, Saturated Baseline, and Candidate Global Features

### 1.1 Candidate Domain and Exact Targets
- **Primitive points:** Canonical pairs $(a, b) \in \mathbb{Z}^2$ with $\gcd(a, b) = 1$ and $b > 0$ (or $(1,0)$ excluded), globally deduplicated by $(a,b)$ before any split.
- **Two sampling regimes:**
  - **Regime U (Unconditioned Skewed Box):** Primitive $(a,b) \in [-A, A] \times [1, H]$ with $A/H \approx \text{skew}$; norms $N_0(a,b) = |G(a,b)|$, $N_1(a,b) = |F(a,b)|$.
  - **Regime Q (Special-$q$ Sublattice, strictly $q > B_1$):** Primitive $(a,b) = i(a_0,b_0) + j(a_1,b_1)$ on $\Lambda_{q,\rho}$ with **prime $q > B_1$**; $N_0(a,b) = |G(a,b)|$, $N_1(a,b) = |F(a,b)|/q$. *(Adversarial rule: CADO's default `qmin` satisfies `qmin < B_1` in 6/6 `cadobench` testbeds `c60`–`c110`; all $q \le B_1$ must be excluded or shifted above $B_1$.)*
- **Cutoffs and exact arithmetic (`gmpy2`):** Fix smoothness bounds $B_s = 2^{\text{lpb}_s}$ ($s \in \{0,1\}$) and test **two pre-registered small-prime cutoffs** $y \in \{2000,\; 32768\}$ ($y < \min(B_0, B_1)$). Compute via Bernstein product/remainder trees (`src/cadobench/smoothcount.py`):
  $$P_y = \prod_{p \le y} p, \qquad S_{s,y}(a,b) = \prod_{p \le y} p^{v_p(N_s(a,b))}, \qquad C_{s,y}(a,b) = \frac{N_s(a,b)}{S_{s,y}(a,b)}.$$
- **Deterministic strata (must be separated):** Since $C_{s,y}$ has no prime factors $\le y$, let $p_{\text{next}}(y)$ be the least prime $> y$.
  - If $C_{s,y} \le B_s$: $C_{s,y}$ is **deterministically $B_s$-smooth** ($Y_s = 1$).
  - If $B_s < C_{s,y} < p_{\text{next}}(y)^2$: $C_{s,y}$ cannot factor into $\ge 2$ primes $> y$, so it is **deterministically non-smooth** ($Y_s = 0$).
  Primary evaluation is restricted to the **uncertain stratum** $\mathcal{U} = \{C_{0,y} \ge \max(B_0 + 1, p_{\text{next}}(y)^2)\} \cap \{C_{1,y} \ge \max(B_1 + 1, p_{\text{next}}(y)^2)\}$ (and $\mathcal{U}_s$ per side); full-box metrics are reported alongside with deterministic flags inside the baseline.
- **Primary targets (binary):**
  1. $Y_{\text{joint}}(a,b) = \mathbf{1}[C_{0,y}(a,b) \text{ is } B_0\text{-smooth and } C_{1,y}(a,b) \text{ is } B_1\text{-smooth}]$;
  2. $Y_0(a,b) = \mathbf{1}[C_{0,y}(a,b) \text{ is } B_0\text{-smooth}]$ and $Y_1(a,b) = \mathbf{1}[C_{1,y}(a,b) \text{ is } B_1\text{-smooth}]$ ($10^2$–$10^3\times$ higher hit counts, giving sharp single-side power).
- **Pre-registered high-power secondary targets:**
  3. Continuous $(y, B_s]$ valuation fraction $W_s(a,b) = \frac{\log \text{smooth\_part}_{B_s}(C_{s,y})}{\log C_{s,y}} \in [0, 1]$ ($O(1)$ variance on every point);
  4. Cross-side residual covariance $\mathrm{Cov}(Y_0 - \hat{p}_{0,\text{base}},\, Y_1 - \hat{p}_{1,\text{base}})$ testing cross-side coupling beyond size and primes $\le y$.

### 1.2 Baseline Hierarchy
- **Baseline B0 (Analytic Rough-Smooth Saddle-Point):** Exact conditional probability $P_{\text{rough}}(s_s; y, B_s)$ that a random $y$-rough integer of log-size $s_s = \log C_{s,y}$ is $B_s$-smooth (Buchstab-Dickman Laplace saddle point over primes in $(y, B_s]$, NOT plain $\rho(s_s/\log B_s)$), shifted by the polynomial's constant $(y, B_s]$ root-density offset $\alpha_{s,(y,B_s]}(F_s)$ computed solely from $F_s$'s coefficients.
- **Baseline B1 (Saturated Size + Local-$y$ Comparator — Mandatory Null):** Out-of-fold cross-fitted ensemble (30-knot monotone/cubic B-spline GAM + isotonic calibration + depth-4 GBDT) using **only**:
  1. Exact log-sizes $s_0 = \log C_{0,y}$, $s_1 = \log C_{1,y}$, raw $\log N_0, \log N_1$, $\log S_{0,y}, \log S_{1,y}$, interaction $s_0 \otimes s_1$, and exact boundary indicators;
  2. Exact local-$y$ sieve state: $(v_p(N_0), v_p(N_1))_{p \le y}$, projective residues $(a:b) \bmod p^{e_p}$ ($p^{e_p} \le y$), and signs $(\mathrm{sgn}(G), \mathrm{sgn}(F))$ (which determine $|F| \bmod 4$ given $F \bmod 4$);
  3. Polynomial-level constant offset (computed from public $(F,G)$ root counts or absorbed by per-$N$ centering) so polynomial-quality variance across $N$ cannot masquerade as within-$N$ candidate guidance.

### 1.3 Creative Catalogue of Candidate Global Features (8 Families)
Every feature below is pre-sanitized by the rules in Section 2 and orthogonalized against $(s_0, s_1, \log S_{0,y}, \log S_{1,y})$:

- **G1. $y$-Truncated Coordinate Arithmetic & $S$-Unit Structure of $(a,b)$:**
  $\log S_y(|a|)/\log|a|$, $\log S_y(|b|)/\log|b|$, indicators $\mathbf{1}[|a| \text{ is } y\text{-smooth}]$, $\mathbf{1}[|b| \text{ is } y\text{-smooth}]$ ($S$-unit coordinates), $\log S_y(|a \pm b|)/\log|a \pm b|$, squarefree core ratio $\log \mathrm{rad}_y(ab)/\log S_y(|ab|)$. *(Tests $abc$-type / $S$-unit smoothness repulsion or attraction without inspecting primes $> y$ dividing $a, b$.)*
- **G2. Continued-Fraction & Diophantine Approximation Structure of $a/b$:**
  CF depth $L(a/b)/\log b$, max partial quotient $\log(\max_k a_k)/\log b$, Khinchin mean $\frac{1}{L}\sum \log a_k$, Gauss–Kuzmin histogram KL, shared CF prefix length and size-residualized Liouville exponent $-\log|a/b - \theta_j|/\log b$ relative to each real root $\theta_j$ of $F(x,1)=0$ and $-g_0/g_1$, plus distance to low-denominator Farey fractions.
- **G3. Minkowski Conjugate-Embedding Balance, Real Chambers & Binary-Form Covariants:**
  Zero-sum embedding vector $\ell_j(a,b) = \log|a - \theta_j b| - \frac{1}{d}\sum_{i=1}^d \log|a - \theta_i b|$ ($\sum_j \ell_j \equiv 0$, orthogonal to $\log|F(a,b)|$): $\|\ell\|_2$, $\min_j \ell_j$, $\max_j \ell_j$, entropy of softmax$(-\ell)$, real-root chamber signature $\sigma(a,b) = (\mathrm{sgn}(a - \theta_j b))_{j=1}^{r_1} \in \{\pm 1\}^{r_1}$ *(tests narrow class group / Artin reciprocity parity constraints)*, scale-free Jacobian $|b F_a(a,b)|/|F(a,b)|$, and Hessian invariant $\log|H_F(a,b)| - \frac{2(d-2)}{d}\log|F(a,b)|$.
- **G4. Sanitized Shifts, Unimodular/Farey Neighbours & Galois/Automorphism Relatives:**
  Residualized outcomes $(Y_{s'} - \hat{p}_{\text{base}}, W_{s'} - \hat{w}_{\text{base}})$ at:
  (i) **Safe shifts** $(a+u, b+v)$ explicitly verified to satisfy $\gcd(G(u,v)F(u,v), P_{(y,B]}) = 1$;
  (ii) **Bezout / Farey neighbours** $(a', b')$ with $a b' - a' b = \pm 1$ (e.g. penultimate CF convergent of $a/b$, which lies far away in Euclidean norm yet satisfies $\gcd(a b' - a' b, P_{(y,B]}) = 1$ identically!), or in Regime Q ($q > B_1$) unimodular lattice steps $i j' - i' j = \pm 1$ checked for $\gcd(\mathrm{Res}_x(F, F \circ M), P_{(y,B]}) = 1$;
  (iii) **Galois / involution / $SL_2(\mathbb{Z})$ images** when $F$ has algebraic symmetries or subfield structure, resultant-checked.
- **G5. Multiplicative Relations Among Candidates:**
  On composite grid points $(a,b) = (a_1 a_2, b_1 b_2)$ and bilinear/Brahmagupta combinations $(a_1 b_2 \pm a_2 b_1, b_1 b_2)$: residualized smoothness of factor candidates $(a_1, b_1), (a_2, b_2)$ filtered to enforce $\gcd(a b_r - a_r b, P_{(y,B]}) = 1$.
- **G6. Descent / Parent-Seed Smoothness Inheritance ("Planting Seeds"):**
  In Regime Q ($q > B_1$), compare **Seeded-$q$** lattices where $(q, \rho)$ arose as a large prime factor of an otherwise $B$-smooth parent relation $(a_{\text{par}}, b_{\text{par}})$ against **Control-$q$** lattices matched on $(\log q$, Gauss-reduced basis skewness/norm envelope, and local $\le y$ roots) from rough parents. Evaluate on $\Lambda_{q,\rho}$ after deleting $(a_{\text{par}}, b_{\text{par}})$ and any point with $\gcd(a b_{\text{par}} - a_{\text{par}} b, P_{(y,B]}) > 1$, reporting the effect **before** (positive via basis-shortness selection) and **after** Baseline B1 residualization.
- **G7. Search-History Trajectories, Ultra-$B^+$ Cosets & Other-Candidate Cofactor Digits:**
  (i) Exponentially weighted moving averages of $(Y - \hat{p}_{\text{base}})$ along adaptive/sequential walks across previously visited disjoint $q' > B_1$ lattices;
  (ii) Coset hashes / Legendre symbols modulo primes $q^* > \max(N_0, N_1)$ *(strictly above maximum norm so $q^* \mid N_s$ is impossible)*;
  (iii) Bit popcount, run-length entropy, and mantissa statistics of **other** sanitized candidates' cofactors $C_{s,y}(a', b')$ ($\gcd(a b' - a' b, P_{(y,B]}) = 1$; never this candidate's $C_{s,y}$).
- **G8. Raw $(a,b)$ Bit-Planes, Carry Cascades & Digital Complexity:**
  Binary/ternary Hamming weights, NAF weight, Lempel–Ziv complexity, carry count in $a + b$ and $g_1 a + g_0 b$, plus learned MLP/1D-CNN/GBDT embeddings of raw $(a,b)$ bits — **evaluated strictly on held-out $N$ (Split A) without polynomial mod-$p$ inputs**.

---

## 2. Adversarial Leakage Checklist, Residualization, Splitting, and Multiple Testing

### 2.1 Seven Non-Obvious Leakage Pathways (Audited on `cadobench_testbed_c60..c110`)
1. **L1 — Shift & Parallel-Root-Line Leakage (`a ± 1` is LIVE-LEAKING in 5/6 `cadobench` testbeds!):**
   Excluding shifts $(u,v) \equiv (0,0) \pmod p$ ($p \in (y,B]$) is **not sufficient**. Because $G(a+u, b+v) = G(a,b) + G(u,v)$, **every prime $p \in (y, B_0]$ dividing $G(u,v)$ leaks 100% identically**; and on side 1, **every prime $p \in (y, B_1]$ dividing $F(u,v)$ leaks along the entire root sublattice $a \equiv (u v^{-1}) b \pmod p$**!
   - *Measured audit (`verify_spec_numbers.py`):* For $(u,v) = (\pm 1, 0)$ ($a \pm 1$) at $y = 2000$, $G(1,0) = Y_1$ in Kleinjung `polyselect` contains prime factors **inside $(2000, B_0]$ in 5 of 6 testbeds**: `c70` ($2843$), `c80` ($13267, 15083$), `c90` ($10459, 12073$), `c100` ($11383, 12457$), `c110` ($13219, 14867$).
   - *Rule:* Any shift $(u,v)$ MUST satisfy $\gcd(G(u,v), P_{(y,B_0]}) = \gcd(F(u,v), P_{(y,B_1]}) = 1$ (and general maps $\Phi$ must satisfy $\gcd(\mathrm{Res}(F, F\circ\Phi), P_{(y,B_1]}) = 1$).
2. **L2 — Projective-Line Collision in "Far Away" Candidates (`p | (a b' - a' b)`):**
   Euclidean distance $\|(a-a', b-b')\| \gg B$ does **not** prevent congruence sharing! Two points share the exact projective class in $\mathbb{P}^1(\mathbb{F}_p)$ (hence identical divisibility by $p$ on both sides) iff $p \mid (a b' - a' b)$. In $\Lambda_{q,\rho}$, $a b' - a' b = \pm q(i j' - i' j)$. Moreover, $B \gg I, J$ in `c60`–`c110`, so spatial box averages do not average out $p$-lattices.
   - *Rule:* Every reference point $(a', b')$ used in $(a,b)$'s features MUST satisfy $\gcd(a b' - a' b, P_{(y,\max B_s]}) = 1$.
3. **L3 — Nonlinear / Buchstab-Cusp Size Leakage (`log C_{s,y}` curve shape):**
   A $y$-rough integer does **not** follow plain $\rho(\log C/\log B)$: it has exact step discontinuities ($1$ for $C \le B$, $0$ for $B < C < p_{\text{next}}(y)^2$) and Buchstab derivative cusps at $B y^k, B^2, \dots$. Any nonlinear coordinate feature ($|a|, |b|, a/b$, CF length) will beat a parametric $\rho(u)$ baseline by fitting size curvature.
   - *Rule:* Restrict primary evaluation to $\mathcal{U}$, compare strictly against nonparametric Baseline B1, and verify any candidate hit under **narrow-bin size permutation** ($\Delta s_0, \Delta s_1 \le 0.02$ nats).
4. **L4 — Coordinate-Divisor & Root-at-Infinity Leakage (`p | a` or `p | b`):**
   If $p \in (y, B]$ divides $c_d, c_0, g_1,$ or $g_0$, then $p \mid b$ or $p \mid a$ forces $p \mid F$ or $p \mid G$. Even when $\gcd(c_0 c_d g_0 g_1, P_{(y,B]}) = 1$, $p \mid a b$ forces $(a:b) \bmod p \in \{(0:1),(1:0)\}$ (non-roots), shifting smoothness by $-O(\frac{\log p}{p})$.
   - *Rule:* Coordinate arithmetic features of $a, b, a\pm b$ may inspect primes $\le y$ ONLY (`S_y(|a|)`), never primes in $(y, B]$.
5. **L5 — Special-$q \le B_1$ & $B^+$-Coset Roughness Leakage (LIVE in 6/6 `cadobench` testbeds!):**
   - *Measured audit:* `qmin < B_1` in **all 6** testbeds (`c60`: $61961 < 2^{19}$ through `c110`: $240000 < 2^{27}$), and `cadobench` observed 13–25% cross-$q$ duplicate relations. Membership in another CADO $q_2$-coset with $q_2 \le B_1$ directly reveals $q_2 \mid C_{1,y}$! Conversely, exposing whether $(a,b)$ hits a root coset mod $q^* \in (B_s, C_{s,y}]$ reveals $q^* \mid C_{s,y} \implies Y_s = 0$!
   - *Rule:* Enforce $q > B_1$ in Regime Q; restrict auxiliary coset moduli to $q^* > \max N_s$ or mask polynomial root lines.
6. **L6 — Raw-Bit Cunningham Linear-Form Memorization & Self-Cofactor Digit Sums:**
   Base-$2^k$ chunk sums and alternating sums of $C_{s,y}$ compute $C_{s,y} \bmod (2^k \mp 1)$, leaking Mersenne/Fermat/Cunningham primes in $(y, B]$. Even on raw $(a,b)$ bits, **11 primes** dividing $2^m \pm 1$ ($m \le 24$) sit inside $(2000, 2^{21}]$ ($2113, 2731, 5419, 8191, 43691, 61681, 65537, 131071, \dots$), allowing an MLP/GBDT on within-$N$ splits to learn $a - r b \equiv 0 \pmod p$ directly from raw bits!
   - *Rule:* Forbid ALL features of this candidate's $C_{s,y}$ beyond $\log C_{s,y}$; evaluate raw-$(a,b)$-bit models strictly on **held-out $N$** (Split A) without polynomial mod-$p$ inputs, plus a randomized-CRT-label null on Split B.
7. **L7 — Cross-$N$ Polynomial-Alpha Confounding & Cross-$q$ Duplicates:**
   Different $N$'s have different $(y, B]$ root densities $\alpha_{(y,B]}(F_N)$; any $N$-level geometry proxy (skew, box ratio) can predict polynomial quality across $N$ with zero within-$N$ guidance value.
   - *Rule:* Deduplicate $(a,b)$ globally; include exact polynomial-level $\alpha_{(y,B]}(F_N)$ in Baseline B1 and score **mean within-$N$ log-loss gain** $\overline{\Delta\mathrm{LL}}_N$ and within-$N$ size-stratified $\Delta\mathrm{AUC}$.

### 2.2 Residualization Protocol
1. Cross-fit Baseline B1 out-of-fold on training data to obtain calibrated baseline logit $\eta_{\text{base}}(x_i) = \mathrm{logit}\,\hat{p}_{\text{base}}(x_i)$.
2. Pre-orthogonalize every continuous global feature $z_{ij}$ against $(s_0, s_1, \log S_{0,y}, \log S_{1,y})$ via a 20-knot spline regression fit strictly on train: $\tilde{z}_{ij} = z_{ij} - \hat{g}_j(s_0, s_1, \log S_{0,y}, \log S_{1,y})$.
3. Fit the test model with **frozen baseline offset**: $\mathrm{logit}\,P(Y_i = 1) = \eta_{\text{base}}(x_i) + f_\theta(\tilde{z}_i)$ (penalized GLM + shallow monotonic-constrained/regularized GBDT).
4. **Mandatory sanity check for any apparent positive:** Permute $\tilde{z}_i$ within joint $(s_0, s_1)$ micro-cells of width $0.02 \times 0.02$ nats (preserving $N$ and $\mathrm{sgn}(F,G)$); a genuine signal must vanish under cell-permutation and survive cell-stratified conditional logistic regression.

### 2.3 Train / Validation / Test Splitting & Multiple-Testing Control
- **Split Axis A — Held-out $N$ (Primary Decision Split):** Generate $K \ge 16$ independent semiprimes across at least 3 size tiers (e.g. 45-, 60-, 70/80-digit balanced semiprimes, independent seeds): `Train-N` (8), `Val-N` (4), `Test-N` (4, locked until final evaluation).
- **Split Axis B — Held-out Region within $N$ (Secondary Polynomial-Geometry Split):** Partition each $N$ by **disjoint special-$q$ bands** ($q > B_1$, interleaved even/odd $\log q$ strata after global $(a,b)$ deduplication) AND by **buffered angular sectors** in $\arctan(a / (\mathrm{skew}\cdot b))$ separated by a $\ge 5\%$ guard band. Pair Split B with a **synthetic local-free Bernoulli control** $Y^{\text{null}}_i \sim \mathrm{Bernoulli}(P_{\text{B0}}(s_i))$ to confirm $0$ false-positive lift from pipeline artifacts.
- **Multiple-Testing Control:** Pre-register 8 family-level omnibus models (G1–G8), 1 all-feature omnibus model, and $J \le 32$ individual single-feature tests. Screen on `Val-N` using **Westfall–Young step-down min-$P$ permutation** (1,000 within-$(N, \text{size-cell})$ permutations, FWER $\le 0.05$). Lock surviving hypotheses and evaluate once on untouched `Test-N` under **Holm–Bonferroni FWER $\alpha = 0.01$**, with two-sided 99% block-bootstrap CIs clustered at the $(N, q)$-block level.

---

## 3. Planted-Signal Positive Controls & $u \log u$ Effect-Size Scaling

### 3.1 Connecting Log-Loss / KL Gain to $u \log u$, Assumption H, and $L_N[\alpha, c]$
Let $u_s = \log C_{s,y}/\log B_s$ and $U := \sum_s u_s \log u_s$ (so baseline smoothness scales as $p_0 \approx \exp(-(1+o(1)) U)$). Let $r(G) = P(G \mid Y=1)/P(G)$ be the likelihood ratio of global features $G$ ($\mathbb{E}[r]=1$). For rare smooth events ($p_0 \ll 1$), expected per-candidate binary log-loss improvement satisfies:
$$\Delta\mathrm{LL} = \mathrm{LL}(p_0) - \mathrm{LL}(p_1) = p_0 \, D_{\mathrm{KL}}\big(P(G \mid Y=1) \,\|\, P(G)\big) + O(p_0^2).$$
Define the **base-rate-normalized information per smooth hit**:
$$\mathcal{I}_{\text{hit}} := \frac{\Delta\mathrm{LL}}{p_0} \approx \mathbb{E}\!\left[\log \frac{p_1(G)}{p_0} \,\middle|\, Y = 1\right] \quad \text{(in nats / relation)}.$$
*(Note: without constructive sublinear navigation into high-$r(G)$ sets, post-hoc rejection filtering of $\Omega$ pays $1/\alpha_{\text{sel}}$ candidates to gain $\le 1/\alpha_{\text{sel}}$ smoothness by Bayes' rule, yielding zero speedup. Even granting **free $O(1)$-cost constructive navigation** into a guided set achieving log-enrichment $\Lambda(u) \asymp \mathcal{I}_{\text{hit}}(u)$, the size-entropy calculus `notes/lowerbound_results.md:57` yields four regimes:)*

1. **Bounded gain $\Lambda(u) = O(1)$ as $U \to \infty$ (e.g. $\mathcal{I}_{\text{hit}} = 0.2$ nats $\iff 1.22\times$ enrichment):** Saves an **$O(1)$ constant factor** only; both $\alpha = 1/3$ and $c$ in $L_N[1/3, c]$ are unchanged.
2. **Sublinear growth $\Lambda(u) = o(U)$ (e.g. $\Theta(\log u)$ or $\Theta(u)$, as in fixed-$k$ norm splitting $\rho(u/2)^2 / \rho(u) \approx e^{u \log 2}$):** Absorbed into $o(1)$ of $L_N[1/3, c + o(1)]$; neither $\alpha = 1/3$ nor $c$ changes asymptotically.
3. **Proportional $u \log u$ gain $\Lambda(u) = \theta \, U$ for fixed $\theta \in (0, 1)$:** Multiplies the Dickman exponent by $1 - \theta$ ($P \approx \rho(u)^{1 - \theta}$), **refuting Assumption H** (`ρ(u)^{1-o(1)}`) and reducing the **exponential constant** in $L_N[1/3, c]$ from $c_0$ to $c_0 (1 - \theta)^{2/3}$ (re-optimizing polynomial degree $d$ and $\log B$).
4. **Full sub-$1/3$ complexity exponent drop ($\alpha < 1/3$ at standard norm profiles):** Requires $\theta(u) \to 1$, i.e. $\Lambda(u) = (1 - o(1))\, U$ leaving residual log-miss-rate $o(U)$.

### 3.2 Two Mandatory Planted-Signal Positive Controls & Exact Power Table
- **PC1 — Calibrated Gaussian-Tilted Oracle Feature (Tests pipeline sensitivity across $\theta \in \{0.01, 0.02, 0.05, 0.10, 0.25\}$):**
  On the real candidates and real labels $Y$, draw synthetic feature $G_\theta \mid (Y=0) \sim \mathcal{N}(0,1)$ and $G_\theta \mid (Y=1) \sim \mathcal{N}(\sqrt{2\theta U},\, 1)$, passed through nonlinear nuisancewarps + noise dimensions. Exact Bayes KL is $D_{\mathrm{KL}} = \theta U$ nats/hit.
  - *Exact Quadrature Power (`verify_spec_numbers.py`, 99% power at Holm FWER $0.01$ over 20 tests):*
    - At `c60–c70` joint scale ($u_0=3.2, u_1=4.4 \implies U = 10.24$, $p_0 = 5\times 10^{-4}$):
      - $\theta = 0.10$ ($\mathcal{I}_{\text{hit}} = 1.02$ nats/hit): needs **$100$ smooth hits** ($2.0 \times 10^5$ test candidates);
      - $\theta = 0.05$ ($\mathcal{I}_{\text{hit}} = 0.51$ nats/hit): needs **$167$ smooth hits** ($3.3 \times 10^5$ test candidates);
      - $\theta = 0.02$ ($\mathcal{I}_{\text{hit}} = 0.205$ nats/hit): needs **$366$ smooth hits** ($7.3 \times 10^5$ test candidates);
      - $\theta = 0.01$ ($\mathcal{I}_{\text{hit}} = 0.102$ nats/hit): needs **$698$ smooth hits** ($1.4 \times 10^6$ test candidates).
    - On single-side targets ($p_0 \approx 10^{-2}$–$5\times 10^{-2}$), $1.4 \times 10^6$ test candidates yield $>14{,}000$ hits, giving $>99\%$ power down to **$\theta \le 0.002$** ($\mathcal{I}_{\text{hit}} \approx 0.01$ nats/hit).
  - **Pre-registered decision threshold:** Pipeline must recover $\hat{\theta}/\theta \in [0.80, 1.15]$ on PC1 at $\theta \ge 0.02$. A real feature is declared **non-null** iff Holm-$p < 0.01$ and lower 99% CI of $\mathcal{I}_{\text{hit}} > 0.02$ nats/hit on Split A; it is an **exponent/H-breaking lead** only if $\hat{\mathcal{I}}_{\text{hit}}(U) / U$ stays $\ge \theta_{\min} = 0.02$ and non-decreasing across $\ge 3$ increasing $U$-bins!
- **PC2 — Arithmetic Planted Signals (No label peeking):**
  - **PC2a ($O(1)$ leakage detector):** Hold primes in $(y, 2y]$ out of Baseline B1 and feed a shift/bit feature coupled to $(y, 2y]$; verify pipeline flags an $O(1)$-nats/hit ($\theta(U) \propto 1/U \to 0$) constant-factor leak.
  - **PC2b (Algebraic Norm-Splitting Plant, $\Theta(u)$ growth):** Mix irreducible degree-4 forms with Sophie-Germain / reducible forms $F_{\text{split}}(a,b) = Q_1(a,b)Q_2(a,b)$ (or points on a norm-splitting curve) matched in $\log C_{1,y}$ and $\le y$ root counts; verify G2/G3/G8 detect the $\Delta\log P \approx +u_1\log 2$ jump on held-out $N$.

---

## 4. Theoretical Prior: Literature & Feature-by-Feature Predictions

### 4.1 Governing Analytic Number Theory Results
1. **Friable Equidistribution in Progressions & Characters (Hildebrand–Tenenbaum 1986; Fouvry–Tenenbaum 1991/1996; Granville 1993; Soundararajan 2008; Harper 2012; Drappeau 2015):** $\Psi(X, B; q, a)$ satisfies Bombieri–Vinogradov / Elliott–Halberstam-type equidistribution up to $q \le X^{1/2-\varepsilon}$ (and $X^{3/5-\varepsilon}$ for well-factorable weights), with no non-local modulus bias once small primes $\le y$ and size are fixed.
2. **Pretentious Multiplicative Independence & Short Shifts (Granville–Soundararajan 2001/2007; Hildebrand 1985/1989; Matomäki–Radziwiłł 2016):** Bounded multiplicative indicators correlate only through archimedean $n^{it}$ (pure size, absorbed by B1) or low-conductor Dirichlet characters $\chi$ (absorbed by local-$y$ unless a character conductor divides $\Delta_F$ above $y$). Consecutive/shifted values $n, n+k$ or independent linear forms obey the **Hardy–Littlewood local singular-series $\times$ Dickman product law**.
3. **Sparse Algebraic Smooth Tuples & $S$-Unit Equations (Balog–Wooley 1998; La Bretèche–Tenenbaum 2005; Lagarias–Soundararajan 2011):** Balog–Wooley construct infinitely many $k$-tuples of consecutive smooth integers by exploiting cyclotomic/algebraic polynomial factorizations $n^m - 1$ on **density-$O(X^{-1+1/m})$ thin sets**, confirming that super-Dickman clustering arises *only* where norms algebraically factor. For $S$-unit sums $a + b = c$ ($a, b$ $y$-smooth), La Bretèche–Tenenbaum prove counts match the local singular-series saddle-point heuristic (no $\Omega(u\log u)$ anomaly; $abc$ bounds extreme smooth collisions).
4. **Polynomial Values (Murphy 1999; Martin 2002; Dartyge–Martin–Tenenbaum):** For irreducible $F \in \mathbb{Z}[a,b]$, smoothness depends asymptotically only on $|F(a,b)|$ and local root densities $\alpha_F$, independent of Minkowski embedding skewness or CF digits.

### 4.2 Summary Prior Table (Expected Clean-Pipeline Outcomes)

| Feature Family | Predicted $\Omega(u\log u)$ Signal | Predicted Finite-$y$/$u$ Residual (Clean Pipeline) | Mechanism / Citation |
|---|---|---|---|
| **G1** ($y$-smooth $a,b$, $a\pm b$) | **NULL** ($\theta = 0$) | Tiny **non-null** $O(\frac{\log y}{y}) \lesssim 4\times 10^{-3}$ nats/hit | $p \nmid ab$ ($p>y$) avoids non-root $(1:0),(0:1) \in \mathbb{P}^1(\mathbb{F}_p)$; La Bretèche–Tenenbaum (2005). |
| **G2** (CF of $a/b$, root Liouville) | **NULL** ($\theta = 0$) | **NULL** ($0$ within noise once $\log C_{s,y}$ splined) | Norm reduction already in $\log C_{1,y}$; Martin (2002), Murphy (1999). |
| **G3a** (Conjugate balance $\ell_j$, $H_F$) | **NULL** ($\theta = 0$) | **NULL** on irreducible $F$ (**NON-NULL** $\Theta(u)$ on PC2b reducible $F$) | Prime ideal norms see only $\prod |a-\theta_j b|$; splits only if $F$ factors over $\mathbb{Q}$. |
| **G3b** (Real-root sign chambers $\sigma$) | **NULL** ($\theta = 0$) | Tiny **non-null** at $u \in (1,2)$, decaying as $O(e^{-cu})$ | Narrow-class / Artin reciprocity parity at ramified $q \mid \Delta_F, q > y$; Chebyshev bias $O(1)$. |
| **G4 / Cross-side** (Safe shifts, Bezout, $Y_0 \leftrightarrow Y_1$) | **NULL** ($\theta = 0$) | Tiny **negative non-null** $-O(\frac{u^2}{y\log y}) \lesssim -10^{-3}$ nats/hit | Coprime $(y,B]$ singular-series tail ($p \mid G(a',b') \Rightarrow p \nmid G(a,b)$); Hildebrand, Granville–Soundararajan. |
| **G5** (Multiplicative $(a,b)$ products) | **NULL** ($\theta = 0$) | **NULL** (unless $F$ is multiplicative / $d=2$ norm composition) | No homomorphism from $(a_1 a_2, b_1 b_2)$ to $F$ for $d \ge 3$. |
| **G6** (Smooth-parent descent $q > B_1$) | **NULL** ($\theta = 0$) | **Positive BEFORE B1**, **NULL AFTER B1** | $\Lambda_{q,\rho} = \mathfrak{q}\cap \mathbb{Z}[\theta]$ is memoryless of parent $(a_{\text{par}},b_{\text{par}})$ once reduced-basis norm size is conditioned out! |
| **G7 / G8** (Walk history, $q^* > N_{\max}$, raw bits on Split A) | **NULL** ($\theta = 0$) | **NULL** on Split A (**LARGE FALSE POSITIVE** on Split B if L1/L5/L6 violated!) | Fouvry–Tenenbaum, Drappeau (2015); Mauduit–Rivat digital independence. |
