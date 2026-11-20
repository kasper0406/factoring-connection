# Hatch 1 Technical Report: Why Capacity-Achieving List Decoders Exist for Polynomial Codes (`F_q[T]`) But Not for Integer CRT Smooth-Search (`Z`), and an Assessment of Buium `p`-Derivations

Date: 2026-11-20
Workspace: `/Volumes/git/factoring-connection/codex/scratch/hatch_h1g/`
Context references:
- `notes/lowerbound_results.md:25-60, 107-120` (Hatch 1: breaking conditional Dickman `H` via beyond-Johnson CRT list decoders on `N`-public residue tables)
- `notes/crtdecode_results.md:109-274` (CRT code `C_y`, valuation metric, `Prop F1-F4`, capacity row 262, 2026 plain-RS capacity rows 263-264)
- `notes/sinonobliv_results.md:141-196` (`Lemma IC`, `Prop RD`, multi-sample ACD route 3)
- `notes/pilist_results.md:251-365` (`Prop JR`, `H_sel`, beyond-Johnson crux)
- `notes/sihard_results.md:402-548` (`Theorem J`: averaged `rho(u_0)^{1-o(1)}` ceiling for single-integer Johnson certificates)
- `notes/subthird_results.md:75-80` (`M2`: class-group positive-definite 2D-CRT constructor at `beta = 1/4`)
- `notes/multcopp_results.md:19-81` (supersingular CM multiplicity over `Z` requires derivation order `a >= p ~ sqrt(N)`)
- `notes/conn_alg.md:140-146` and `notes/nonsmooth_results.md:313-350` (`Prop FROB` and Buium's arithmetic derivative mod `N`)

Labels used throughout: `KNOWN` (published/cited) · `PROVED` (exact proof given here) · `PROVED in MODEL` · `HEURISTIC` (named assumption) · `SPECULATIVE` · `MEASURED` (pre-registered toy run).

---

## 1. Unified Algebraic Dissection: How Polynomial Codes Cross the Johnson Wall (`sqrt(R)` -> `R`)

### 1.1 Why Nonlinear Univariate/Bivariate Interpolation Stops at Johnson (`KNOWN` + unified view)

Let a code encode messages from a set `M` of entropy `b = ln |M|` nats into a codeword of total log-modulus `W` nats (`W = n ln q` for Reed-Solomon over `F_q`; `W = ln P = sum_{q <= B} e_q ln q` for integer CRT). Suppose a received word agrees with a message `m in M` on log-weight `A` nats (`alpha_agr = A / W`, rate `R = b / W`).

In **Guruswami-Sudan** (for RS) and **Boneh / Howgrave-Graham / Guruswami-Sahai-Sudan** (for CRT, `notes/crtdecode_results.md:256-259`), the decoder builds a degree-`L` polynomial in a single evaluation channel `Y` vanishing to multiplicity `r <= L` on the ideal `(P, Y - w)^r`:
1. **Constraint / index cost:** Vanishing to order `r` in dimension `1 + 1` (or `P^{max(0, r-j)} (Y - w)^j` in `Z[Y]`) consumes `(r(r+1)/2) W` nats of determinant.
2. **Degree / archimedean cost:** Because the basis uses **nonlinear powers** `1, Y, Y^2, ..., Y^L`, the `j`-th term has size `j * b` (`deg(f^j) = j(k-1)` or `ln|m^j| = j b`), accumulating `(L(L+1)/2) b` nats in the scaled determinant.
3. **Root-forcing balance:** Dividing the `(L+1)`-dimensional average norm by `r` gives
   `A  >  [ (r+1) / (2(L+1)) ] W  +  [ L / (2 r) ] b`.
   Optimizing `theta = r / L in (0, 1]` by AM-GM hits the **Johnson floor**:
   `A  >=  sqrt(b W)`   (`<==>` `alpha_agr >= sqrt(R)`).

### 1.2 The Three Algebraic Pillars That Let Polynomial Decoders Reach Capacity (`1 - R - eps`)

Every capacity-achieving polynomial list decoder replaces the **height-growing powers** `(1, f, f^2, ..., f^s)` (`size = j * b`) with **`s` linearly-related views** `(L_0(f), L_1(f), ..., L_{s-1}(f))` whose degrees/heights **all stay bounded by `1 * b`**! When `Q_0(T) + sum_{j=1}^s Q_j(T) L_{j-1}(f)(T) = 0` is interpolated with multiplicity `1` across `s + 1` linear channels (packing `m >= s` evaluations per local ideal so `m - s + 1` equations hold per agreeing block):
- Constraint cost is `((m - s + 1)/m) W` shared across `s + 1` polynomials `Q_0, ..., Q_s`.
- Degree of `E(T) = Q_0(T) + sum_{j=1}^s Q_j(T) L_{j-1}(f)(T)` is bounded by `D + k - 1`, paying **`1 * b`** instead of `(s/2) b`!
- Dividing by `(m - s + 1)/m` yields the **Linear-Algebraic Capacity Law**:
  `A  >  W / (s + 1)  +  [ m / (m - s + 1) ] b`   (`--> b + eps W` as `1 << s << m`).

For this law to hold and yield a small list, `R = F_q[T]` provides **Three Algebraic Pillars**:
- **Pillar P1 (Global Prime-Independent Operator + Local Jet Compatibility):** Each `L_j: M -> R` is defined globally on `R` independent of the prime ideal `p_i = (T - alpha_i)`, yet `f mod p_i^m` determines `L_j(f) mod p_i^{m - s + 1}` simultaneously at every `p_i`.
  *Root cause in `F_q[T]`:* Equal characteristic (`F_q subset F_q[T]`). Every prime ideal `(T - alpha_i)` is a translate of **one global coordinate `T`** by a **constant `alpha_i` killed by `d/dT`** (`d(T - alpha_i) = dT` is independent of `i`!).
- **Pillar P2 (`F_q`-Linearity + Degree/Height Non-Expansion):** `L_j` is `F_q`-linear and `deg(L_j(f)) <= deg(f) < k` (`size(L_j(f)) <= b`, never `j * b`).
- **Pillar P3 (Bounded-Dimension Kernel / Wronskian Non-Degeneracy):** Once `E(T) = 0` globally, the operator equation `Q_0 + sum_{j=1}^s Q_j L_{j-1}(f) = 0` over `F_q(T)` confines `f` to an affine `F_q`-subspace of dimension `<= s - 1` (size `q^{s-1}`, independent of `k`, pruned to `O(1/eps)` via subspace-evasive sets or subspace designs).

### 1.3 Exact Algebraic Ingredient Table Across All Seven Literature Families (`KNOWN` + synthesis)

| Decoder Family | Code Changed? | Exact Algebraic Ingredient (P1 + P2) | List-Extraction / Pruning Ingredient (P3) | Why No `Z`-CRT Smooth-Search Transfer |
|---|---|---|---|---|
| **Parvaresh-Vardy** (FOCS 2005) | Yes (`s`-tuple) | Power `h_j(T) = f(T)^{a^j} mod E(T)` reduced modulo an irreducible `E(T)` of deg `k` so `deg(h_j) < k` (`size = b`, not `a^j b`) | Quotient field `F_q[T]/(E(T)) =~ F_{q^k}` where `h_j = f^{a^j}`, turning `Q(T, Y, Y^a, ...) = 0 mod E` into degree-`< a^s` univariate root-finding over `F_{q^k}` | `Z` has no degree-reducing quotient `mod E` that commutes with `mod q_i`: `(m^{a^j} mod M_0) mod q_i` cannot be computed from `m mod q_i` because of integer wraps `floor(m^{a^j}/M_0)`! |
| **Folded RS** (Guruswami-Rudra 2006/2008; Guruswami-Wang 2013; KRSW 2018) | Yes (`m`-folded orbits) | Global ring automorphism `T |-> gamma T` (`gamma in F_q^*` of order `>= n`), preserving `deg f < k` (`size = b`) and sliding across `(alpha, gamma alpha, ..., gamma^{m-1} alpha)` | Diagonal monomial action `(T^i) |-> gamma^{i j} T^i` + periodic subspace designs (Guruswami-Kopparty 2016) / KRSW subspace-evasive pinning (`dim <= s - 1`) | `Aut(Z) = {id}`! On `O_K`, `Gal(K/Q)` fixes `Z` and permutes split prime ideals `q_1 != sigma(q_1)`; primitive smooth norms never vanish mod `q O_K` (Theorem 2 below). |
| **Multiplicity Codes** (Kopparty-Saraf-Yekhanin 2011; Kopparty 2015) | Yes (`m` Hasse derivs) | Global `F_q`-linear Hasse derivations `partial^{(j)}` (`f(T+Z) = sum f^{(j)}(T) Z^j`) with `deg f^{(j)} = deg f - j < k` + Leibniz rule `E^{(u)}(alpha) = 0` for `0 <= u <= m - s` | Order-`(s-1)` algebraic linear ODE solved by Taylor/Hensel lifting from `s-1` initial jets + subspace-evasive sets (DKSS 2009) / subspace designs | `Z` has no non-zero global additive derivation (`d(1) = 0 ==> d(n) = 0` by `1 + ... + 1`); `q`-adic digits depend on `q` and carries break additivity; also `s <= k` forces `alpha_agr >= 1/(k+1) >> 1/pi(B)` (Prop 1d). |
| **BHKS Ideal Framework** (Bhandari-Harsha-Kumar-Sudan, RANDOM 2021 / IEEE IT 2024, arXiv:2103.07930) | Yes (affine FRS ideals) | Coprime ideals `I_i subset F_q[T]` generated by orbits/jets of a **commuting global degree-preserving linear operator algebra** (affine maps `a T + b` + Hasse derivations) | Generalized operator Wronskian matrix `det(L_u(f_v)) != 0` bounding kernel dimension `<= s - 1` via auxiliary minimum-distance codes | Even in `F_q[T]`, BHKS does *not* decode arbitrary coprime ideals to capacity — only ideals invariant under global degree-preserving operators on `F_q[T]`, which `Spec Z` lacks. |
| **Plain-RS Capacity** (Brakensiek-Chen-Putterman-Zhang-Zheng, arXiv:2609.08005, Sept 2026) | **No** (`m = 1` plain RS, prime `q`) | **Hidden Hasse derivatives** `(Y_0, Y_1, ..., Y_d) = (f, f^{(1)}, ..., f^{(d)})` as formal jet variables coupled by total derivative `D_T` over a **prime field `F_q`** (`q >= n`) | Prime characteristic `q > k` prevents Frobenius kernel `F_q[T^p]`; Kolchin-Ostrowski Wronskian independence bounds hidden-jet curve/subspace dimension; time `q^{O(eps^{-12/theta})}` | Needs global derivation `partial^{(j)}` over one prime field `F_q`, AND `d = ceil(eps^{-3/theta}) < k`. At `eps ~ 1/pi(B)`, `d = pi(B)^{3/theta} >> k = poly(x)` fails (`f^{(d)} == 0`) and `q^{O(d^4)}` exceeds sieving (`crtdecode_results.md:263`). |
| **Plain-RS Capacity & Proximity Gaps** (Jeronimo, arXiv:2609.05870, Sept 2026) | **No** (`m = 1` plain RS, prime `q`) | Reparameterized hidden-derivative interpolation + exact-support mutual correlated agreement (MCA) over prime `F_q` (`q >= n`) at fixed additive slack `gamma = alpha_agr - R > 0` | Field-size-independent list bound `n^{O_gamma(1)}` via curve proximity gaps; runtime `q^{C(gamma, L)}` | Same missing global derivation over `Z`, AND additive slack `gamma = alpha_agr - R ~ 1/pi(B) -> 0` blows up `C(gamma, L) -> infinity` (`crtdecode_results.md:264`). |

---

## 2. New Obstruction Results for Arithmetic Analogues (`PROVED`)

### 2.1 Theorem 1 (`Theorem Int-J` & `Theorem Buium-Lin`): Why Buium `p`-Derivations and `Int(Z)` Cannot Beat Johnson

Recall Joyal/Buium's `p`-derivation on `Z` (Buium, *Invent. Math.* 1995; *Arithmetic Differential Equations*, AMS 2005): for a prime `p`,
`delta_p(x) := (x - x^p) / p : Z -> Z`,   `delta_p^{(0)}(x) = x`,   `delta_p^{(u+1)}(x) = delta_p(delta_p^{(u)}(x))`.
On `Z`, `delta_p` satisfies the nonlinear carry rule `delta_p(x + y) = delta_p(x) + delta_p(y) + (x^p + y^p - (x+y)^p)/p` and jet compatibility: if `x == w (mod p^m)` (`m >= 1`), then `x^p == w^p (mod p^{m+1})`, so `delta_p^{(u)}(x) == delta_p^{(u)}(w) (mod p^{m - u})` for `0 <= u < m`.

We prove four exact obstructions closing every way of assembling `delta_p` into a CRT multiplicity decoder over `Z`:

- **Proposition 1a (Local `w = 0` Redundancy, `PROVED`).**
  *For any prime `q`, integer `z != 0` with `v_q(z) = m >= 1`, and `0 <= u <= m`, we have `v_q(delta_q^{(u)}(z)) = m - u` and `q^u delta_q^{(u)}(z) == z (mod q^{m + 1})`. Consequently, the `s` jet congruences `delta_q^{(u)}(z) == 0 (mod q^{m - u})` (`0 <= u < s <= m`) are identically equivalent to the single divisibility `q^m | z`.*
  **Proof.** Write `z = q^m t` with `q \nmid t` and `m >= 1`. Since `m q >= m + 1` (for all `q >= 2, m >= 1`), `z^q = q^{m q} t^q == 0 (mod q^{m + 2})` when `m >= 2` (and `== 0 mod q^2` when `m = 1, q = 2`). Specifically, `q * delta_q(z) = z - z^q = q^m t (1 - q^{m(q-1)} t^{q-1})`, where `1 - q^{m(q-1)} t^{q-1} == 1 (mod q^{m(q-1)})` and `m(q-1) >= 1`. Thus `v_q(delta_q(z)) = m - 1` and `q * delta_q(z) == z (mod q^{m + 1})` whenever `m(q-1) >= 1` (`m + m(q-1) = m q >= m + 1`). Induction on `u` gives `v_q(delta_q^{(u)}(z)) = m - u` and `q^u delta_q^{(u)}(z) == z (mod q^{m + 1})` (or `q | z` at the last step `m - u = 1`). Converse: `q^{m-u} | delta_q^{(u)}(z)` at `u = 0` already gives `q^m | z`. QED.

- **Proposition 1b (Local-to-Global Prime Incoherence, `PROVED`).**
  *Unlike `f^{(1)}(T) in F_q[T]`, `delta_q(j) = (j - j^q)/q` depends explicitly on `q`. For any `j in Z \ {-1, 0, 1}` (`j = -1` gives `delta_q(-1) = 0` for odd `q` and `-1` for `q = 2`), `delta_{q_1}(j) != delta_{q_2}(j)` for all distinct primes `2 <= q_1 < q_2` with `|j| >= 2`.*
  **Proof.** For `|j| >= 2`, the function `g(t) = (|j|^t - |j|) / t` is strictly increasing on `t >= 2` because `(d/dt)(e^{a t} - e^a)/t > 0` for `a = ln|j| >= ln 2` and `t >= 2` (`(a t - 1) e^{a t} + e^a > (2 ln 2 - 1) 4 + 2 > 0`). For `j >= 2`, `delta_q(j) = -g(q)` is strictly decreasing in `q`; for `j <= -2` and odd `q`, `delta_q(j) = +g(q)` is strictly increasing in `q` (while `delta_2(j) = (j - j^2)/2 < 0 < g(q_2)`). Thus no global integer equals `delta_q(j)` across two primes, and CRT-gluing `(delta_q(j) mod q^{e_q - 1})_{q <= B}` produces an integer of generic size `P / rad(P) = exp((1 - o(1)) W) >> H`. QED.

- **Theorem 1c (`Theorem Buium-Lin`: Fixed-Prime Linear Buium Jet Interpolation Degrades Exponentially, `PROVED`).**
  *Fix any prime `p_0 >= 2` and order `s >= 2`, and consider the `(s+1)`-dimensional linear Buium jet lattice*
  `Lambda_{p_0, s} = { (c_0, c_1, ..., c_s) in Z^{s+1} : c_0 + sum_{u=0}^{s-1} c_{u+1} delta_{p_0}^{(u)}(w) == 0 (mod P') }`
  *(`P' | P`, `gcd(p_0, P') = 1`, `ln P' = W'`). On `[0, H]` (`H = e^b >= 4`), `delta_{p_0}^{(u)}(Y)` is a polynomial of exact degree `p_0^u` with leading coefficient `+- p_0^{-(p_0^u - 1)/(p_0 - 1)}`, so `M_u := max_{0 <= j <= H} |delta_{p_0}^{(u)}(j)| = H^{p_0^u} p_0^{-(p_0^u - 1)/(p_0 - 1)} (1 - O(H^{1 - p_0}))`. Hence the scaled lattice covolume satisfies*
  `covol^{1/(s+1)} = exp( [ W' + ((p_0^s - 1)/(p_0 - 1)) b - O(p_0^s ln p_0) ] / (s + 1) )`.
  *For every `p_0 >= 2` and `s >= 3`, `(p_0^s - 1)/((p_0 - 1)(s + 1)) >= (2^s - 1)/(s + 1) > s / 2`, which is **strictly worse** than ordinary `r = 1` Coppersmith `(1, Y, ..., Y^s)` and diverges as `p_0^s / s -> infinity` rather than approaching `1`.*
  **Proof.** Induction on `u` shows `deg(delta_{p_0}^{(u)}) = p_0^u` and `|[Y^{p_0^u}] delta_{p_0}^{(u)}| = p_0^{-(1 + p_0 + ... + p_0^{u-1})} = p_0^{-(p_0^u - 1)/(p_0 - 1)}`. Since `[Z^{s+1} : Lambda_{p_0, s}] = P'` (`c_0` coefficient is `1`), multiplying by `prod_{u=0}^{s-1} M_u` sums the exponents `sum_{u=0}^{s-1} p_0^u = (p_0^s - 1)/(p_0 - 1)`. For `s >= 3`, `2^s - 1 > s(s + 1)/2` (`s = 3: 7 > 6`; `s = 4: 15 > 10`; induction by `2^s > s + 1`). QED.

- **Theorem 1d (`Theorem Int-J`: Exact Capacity Barrier on the Entire Pólya Ring `Int(Z)`, `PROVED`).**
  *Every polynomial expression in Buium's `p`-derivations `delta_p^{(u)}(Y)` (over ANY finite collection of primes `p` and orders `u`, with integer coefficients and total polynomial degree `<= L`) lies in the Pólya ring `Int(Z)_{<= L} = { F in Q[Y]_{<= L} : F(Z) subset Z } = span_Z { binom(Y, k) : 0 <= k <= L }`. Fix any `P >= 2` (`W = ln P`), `X in Z`, window `[0, H]` (`b = ln H`), multiplicity `1 <= r <= L`, and `H >= 2 L`. Define*
  `J_{P, r, L} := { F in Int(Z)_{<= L} : forall u in Z, P^r | F(-X + P u) }`,
  `I^*_{P, r, L} := { F in Int(Z)_{<= L} : forall d | P, forall u in Z, d^r | F(-X + d u) } subset J_{P, r, L}`.
  *Then:*
  1. `[ Int(Z)_{<= L} : I^*_{P, r, L} ] >= [ Int(Z)_{<= L} : J_{P, r, L} ] >= P^{r(r+1)/2}`; *and for `r = 1` and squarefree `P`, `[ Int(Z)_{<= L} : I^*_{P, 1, L} ] = prod_{q | P} q^{floor(L/q) + 1}`.*
  2. *Let `S = { j_u = floor(u H / L) : 0 <= u <= L } subset {0, ..., H}` and `Psi_S(F) = (F(j_0), ..., F(j_L)) in Z^{L+1}`. Then `Psi_S` is injective and*
     `covol(Psi_S(J_{P, r, L}))  >  P^{r(r+1)/2} * ((H - L) / L)^{L(L+1)/2}  >=  P^{r(r+1)/2} * (H / (2 L))^{L(L+1)/2}`.
  3. *Consequently, the Minkowski root-forcing threshold `T = covol(Psi_S(J_{P, r, L}))^{1/(L+1)} <= d^r` (`ln d >= (1/r) ln T`) requires:*
     `ln d  >  [ (r + 1) / (2(L + 1)) ] W  +  [ L / (2 r) ] (b - ln L - ln 2)  >=  sqrt( W * (b - ln L - ln 2) )`   *(whenever `H > 2 L`).*
  *In the smooth-search regime (`b = Theta(F_{1 - beta})`, `beta < 1/2`), any polynomial-in-`pi(B)` lattice dimension `L = L_N[beta, O(1)]` has `ln L = O(F_beta) = o(b)`, so `ln d >= (1 - o(1)) sqrt(b W)` — **the exact Johnson barrier**! Breaking `sqrt(b W)` via `Int(Z)` denominators `1/k!` requires `ln L >= b - O(1)`, i.e., `L >= Theta(H)`, which is slower than sieving `[0, H]`.*

  **Proof.**
  - *Inclusion:* If `f in Int(Z)`, then `(f(n) - f(n)^p)/p in Z` for all `n in Z` by Fermat's Little Theorem, so `Int(Z)` is closed under `delta_p` for every prime `p`. By Pólya (1915), `{binom(Y, k)}_{k=0}^L` is a `Z`-basis of `Int(Z)_{<= L}`, as is `{binom(Y + X, k)}_{k=0}^L`.
  - *(1) General `r`:* Define `pi: Int(Z)_{<= L} -> (Z / P^r Z)^r` by `pi(F) = (F(-X + P u) mod P^r)_{u=0}^{r-1}`. Clearly `J_{P, r, L} subset ker(pi)`. Restrict `pi` to the rank-`r` direct summand `V_r = span_Z { binom(Y + X, k) : 0 <= k < r }`. The `r x r` integer evaluation matrix `E_{u, k} = binom(P u, k)` (`0 <= u, k < r`) has degree `k` in `u` with leading coefficient `P^k / k!`. By the Vandermonde identity `det(u^k)_{0 <= u, k < r} = prod_{k=0}^{r-1} k!`, we obtain `det(E) = (prod_{k=0}^{r-1} P^k / k!) (prod_{k=0}^{r-1} k!) = P^{r(r-1)/2}`. Letting `s_0, ..., s_{r-1} >= 1` be the Smith invariant factors of `E` (`prod s_i = P^{r(r-1)/2}`), we have
    `[Int(Z)_{<= L} : J_{P, r, L}] >= |pi(V_r)| = prod_{i=0}^{r-1} (P^r / gcd(P^r, s_i)) >= prod_{i=0}^{r-1} (P^r / s_i) = P^{r^2 - r(r-1)/2} = P^{r(r+1)/2}`.
  - *(1) Exact `r = 1` formula:* By CRT across `q | P`, shift `Z = Y - w_q`. By Lucas's theorem, `binom(q u, k) == binom(u, k / q) (mod q)` if `q | k` and `0 (mod q)` if `q \nmid k`. Since `{binom(u, m)}_{m=0}^{floor(L/q)}` is unitriangular on `u = 0, ..., floor(L/q)`, `q | sum_{k=0}^L a_k binom(q u, k)` for all `u in Z` iff `q | a_{q m}` for `0 <= m <= floor(L/q)`, giving index `q^{floor(L/q) + 1}`.
  - *(2) & (3):* In basis `{binom(Y, k)}_{k=0}^L`, `[Y^k] binom(Y, k) = 1/k!`, so `det(Psi_S) = [prod_{0 <= u < v <= L} (j_v - j_u)] / [prod_{k=0}^L k!]`. For `j_u = floor(u H / L)` and `u < v`, `j_v - j_u > (v - u)(H / L) - 1 >= (v - u)((H - L)/L) >= (v - u)(H / (2 L))` since `H >= 2 L`. Taking the product over `0 <= u < v <= L` yields `((H/(2L))^{L(L+1)/2} prod_{k=0}^L k!) / prod_{k=0}^L k! = (H/(2L))^{L(L+1)/2}`. Multiplying by `P^{r(r+1)/2}`, taking the `r(L+1)`-th root, and applying AM-GM along with `L(r+1)/(r(L+1)) >= 1` (since `L - r >= 0`) gives `sqrt(W (b - ln L - ln 2))`. QED.

### 2.2 Proposition 1e (Universal Derivative-Order Barrier at Vanishing Agreement, `PROVED`)
Even over `F_q[T]` (Kopparty multiplicity codes and BCPZZ/Jeronimo hidden-derivative plain-RS codes), every `Hasse-derivative` operator `L_j(f) = f^{(j)}` vanishes identically on `F_q[T]_{<k}` once `j >= k`. Hence any non-degenerate derivative interpolation must choose order `s <= k`, which in `A > W/(s+1) + ...` forces agreement fraction `alpha_agr = A / W > 1 / (k + 1)`. In smooth search (`k = b / ln B = poly(x)`, `alpha_agr ~ 1 / pi(B) = L_N[beta, -kappa]`), `1 / (k + 1) >> 1 / pi(B)` by an `L_N[beta]` factor — **so derivative/multiplicity order collapses before reaching smooth-search agreement even in equal characteristic!**

---

## 3. Galois Automorphisms on `O_K` (`Theorem 2`, `PROVED`)

- **Theorem 2a (Positive Coding-Theory Result: Beyond-Johnson Linear-Algebraic Decoding of Galois-Folded `O_K` CRT Codes, `PROVED`).**
  Let `K/Q` be cyclic Galois of degree `m` with generator `sigma` of `Gal(K/Q)`, rational primes `q_1, ..., q_N` (`P = prod q_i`, `W = m ln P`), and messages `alpha in O_K` with `max_{tau: K -> C} |tau(alpha)| <= B_arch` (`b approx m ln(2 B_arch)`). Given `w mod P O_K` agreeing with `alpha` modulo `D O_K` (`D | P`, full-orbit agreement `A_orbit = m ln D`), LLL on the rank-`m(s+1)` `Z`-lattice `Lambda = { (c_0, ..., c_s) in O_K^{s+1} : c_0 + sum_{j=1}^s c_j sigma^{j-1}(w) == 0 (mod P O_K) }` (`1 <= s <= m`) runs in `poly(m, ln P, ln B_arch, ln|disc K|)` time and forces `E(alpha) = c_0 + sum_{j=1}^s c_j sigma^{j-1}(alpha) = 0` in `K` whenever
  `A_orbit  >  W / (s + 1) + [ s / (s + 1) ] m ln B_arch + (1/2) ln|disc K| + (m^2 (s + 1) / 2) ln 2`.
  By the **Artin-Moore Galois Wronskian Lemma** (`det(sigma^j(x_v))_{0 <= j < s, 1 <= v <= s} != 0` iff `x_1, ..., x_s in K` are `Q`-linearly independent, proved by induction via `K^sigma = Q`), `L(x) = sum_{j=1}^s c_j sigma^{j-1}(x)` has `dim_Q ker(L) <= s - 1`, confining `alpha` to an explicit `(s-1)`-dimensional `Q`-affine subspace containing at most `(2 floor(B_arch) + 1)^{s - 1}` points of `O_K`.
  *Regime note:* For fixed `m >= 4`, `s = m - 1`, `R in (1/s^2, 1)`, and large primes `ln P >> m^2 + ln|disc K|`, `A_orbit / W < sqrt(R)` genuinely beats the Johnson bound (verified experimentally in E3a!), though the worst-case kernel count is `V^{(s-1)/m}` unless either `B_arch` is small (`R -> 0`) or `>= 2` independent short lattice vectors simultaneously vanish (which collapses the intersection to `{alpha}`, as occurs in 100% of E3a hits!).

- **Theorem 2b (Negative Smooth-Search Result: Primitive Smooth Norms Have Zero Full-Orbit Agreement, `PROVED`).**
  Let `alpha in O_K \ {0}` be primitive (`alpha not in q O_K` for every rational prime `q`). Then **no rational prime `q` satisfies `q O_K | alpha`**, so **`A_orbit = 0`** identically! At every split prime `(q) = q_1 ... q_m` dividing `N_{K/Q}(alpha)`, `alpha == 0 (mod q_1)` holds on a non-Galois-invariant ideal `sigma(q_1) = q_2 != q_1`, leaving `sigma(alpha) mod q_1 = sigma(alpha mod q_m)` completely unknown. Restricting the congruence to single prime ideals `q_{i, 1}` divides the modulus norm `N(D_ideal)` by `m` relative to `P O_K` (or reduces `dim` by `m`), restoring the single-ideal Biasse-Quintin / Cohn-Heninger Johnson threshold `sqrt(b W)`.

---

## 4. Exponent Ledger If a Beyond-Johnson Smooth-Search CRT Decoder Worked (`PROVED in MODEL`)

Let `x = ln N`, `F_alpha = x^alpha (ln x)^{1 - alpha}`, `B = L_N[beta, kappa]`, `n = pi(B) = L_N[beta, kappa]` (up to `F_beta^{-1}`), and suppose a one-witness CRT decoder finds one `B`-smooth element in a window/box of size `V = e^b` (`b = ((1 - beta) c_X / kappa) F_{1 - beta} + O(ln x)`, where target candidate size is `ln X = c_X x`, rate `R ~ ((1 - beta)/kappa) F_{1 - beta} / (x pi(B))`, agreement fraction `alpha_agr ~ 1 / pi(B)`) in time `T_dec = pi(B)^a = L_N[beta, a kappa]` (`a >= 1` constant):

1. **1D Interval `Prop SI` (`u in [k N, k N + V]`, `c_X = 1`, `notes/sinonobliv_results.md:88-90, 191-193`):**
   - Residue `v = u - k N <= V` must also be `B`-smooth (tested by ECM in `L^{o(1)}`), costing `1 / rho(b / ln B) = exp(Theta(F_{1 - 2 beta}))` trials.
   - Balancing `max(beta, 1 - 2 beta)` forces **`beta = 1/3`**, with total factoring constant
     `c(a) = min_{kappa > 0} max( (1 + a) kappa + 2 / (9 kappa^2), 2 kappa ) = (3/2) (4/9)^{1/3} (1 + a)^{2/3}` (for `a >= 1/3`):
     - `a < 1.14156`: beats MNFS `1.9019`;
     - `a = 1` (linear in residue table): **`L_N[1/3, 6^{1/3}] = L_N[1/3, 1.8171]`**;
     - `a -> 0`: **`L_N[1/3, (16/9)^{1/3}] = L_N[1/3, 1.2114]`** (and `F_p` DLP in `L_p[1/3, 1.232]`, `sihard_results.md:180`).

2. **2D Class-Group Norm `M2` (`Q_m(k, v)` of disc `-4N`, `c_X = 1/2`, `subthird_results.md:76-79`) or Doubly-Split `Pi_beta` (`pilist_results.md:203-206`):**
   - **No second side to test!** Each `B`-smooth primitive witness `Q_m(k, v)` (sampled on privately squared random class walks to preserve `2`-torsion/genus entropy, `codex/subthird_2_codex.out:19`) immediately yields a certified `Cl(-4N)` relation!
   - For **any `beta in (0, 1/3)`**:
     - `pi(B) = L_N[beta, kappa]` relations, each decoded in `L_N[beta, a kappa]` -> relation collection `L_N[beta, (1 + a) kappa]`.
     - Sparse GF(2) Wiedemann/Lanczos on `pi(B)` columns -> `L_N[beta, 2 kappa]`.
     - **Total factoring complexity: `L_N[beta, max((1 + a) kappa, 2 kappa)]`, achieving `alpha = beta < 1/3` for ANY fixed polynomial exponent `a >= 1`!** (E.g. at `beta = 1/4, a = 1, kappa = 1/2`: `L_N[1/4, 1]`.)

---

## 5. The Precise Missing Lemma

To cross the Johnson bound on `Z` (or `Z^2` for `M2`) at smooth-search parameters (`b = Theta(F_{1 - beta})`, `W = c_X x pi(B)`, `A = c_X x`), a linear-algebraic / arithmetic-jet decoder requires:

**Missing Lemma (Bounded-Height Prime-Coherent Jet Lifting at Order `s > pi(B)`):**
*Construct, in `poly(pi(B), x)` time from the residue table `(r_{q^e})_{q <= B}`, `s >= pi(B) + 1` maps `T_0, T_1, ..., T_{s-1}: M -> Z` on the message box `M` (`ln |M| = b`) and compressed moduli `d_q | q^{e_q}` (`sum_{q <= B} (ln d_q / (e_q ln q)) >= 1 - o(1)`) such that:*
1. *(Local determinacy)* Whenever `Q(m) == 0 (mod q^{e_q})`, the values `(T_0(m), ..., T_{s-1}(m)) mod d_q` are polynomial-time computable from `r_{q^e}` alone;
2. *(Height non-expansion)* `max_{m in M} ln |T_u(m)| <= C_0 * b + o(x)` uniformly for all `0 <= u < s` (preventing the `j * b` Coppersmith sum and `p_0^u * b` Buium sum);
3. *(Kernel compression)* Nonzero solutions `(c_0, ..., c_s)` to the resulting congruence lattice induce `c_0 + sum_{u=1}^s c_u T_{u-1}(m) = 0` whose solution set in `M` can be intersected/pruned to `<= poly(pi(B))` candidates without `s <= k` degeneracy (`k = b / ln B << pi(B)`).

Our Theorems 1a-1e and 2b prove that **none** of `q`-adic digits, Buium `p`-derivations `delta_p`, integer-valued polynomials `Int(Z)`, Hasse/hidden derivatives, or `O_K` Galois automorphisms satisfies (1)+(2)+(3).


---

## 6. Toy Experiments (`MEASURED` + Disclosed `POST-HOC`)

Files in `/Volumes/git/factoring-connection/codex/scratch/hatch_h1g/`:
- Pre-registration (written before any script): `preregistration.md`
- Main experiment script (exact Cohen Alg. 2.6.7 integral LLL, zero float rounding in LLL): `toy_experiments.py`
- Pre-registered raw JSON output: `results.json`
- Disclosed post-hoc diagnostic script: `posthoc_diag.py`

### 6.1 Pre-Registered Scorecard (`results.json`, 17.2s runtime)

| ID | Pre-registered prediction | Result | Status |
|---|---|---|---|
| **P1a** | `[Int(Z)_{<= L} : I^*_{P, r, L}] >= [Int(Z)_{<= L} : J_{P, r, L}] >= P^{r(r+1)/2}` in 100% of `(r, L)`, and `r = 1` matches Lucas's `prod_{q \| P} q^{floor(L/q)+1}` | 5/5 cells exact (`(r,L) = (1,2),(1,4),(2,4),(2,6),(3,6)`); e.g. at `(2,6)`, `ln idx(I^*) = 77.7516 >= 67.7712 = 3 ln P` | **PASS** |
| **P1b** | `Int(Z)_{<= 4}` (`r = 2`) minimum decoded `ln d` on `P = prod_{q<=23} q, H = 2^{12}` is `>= 0.85 * Z[Y]` and `> 0.70 * sqrt(b W)` | `sqrt(b W) = 12.64`; `Z[Y]` min/mean `11.20 / 13.35` vs `Int(Z)` min/mean `11.55 / 13.26` (`1.03x` min, `0.993x` mean) | **PASS** |
| **P1c** | Linear Buium `delta_2` (`s = 3, 4`) at `H = 2^{10}, P_odd = prod_{3<=q<=29} q` requires larger `ln d` than `r = 1` Coppersmith | `s = 3`: Buium `12.84` vs Copp `13.85`; `s = 4`: Buium `12.58` vs Copp `13.66` | **FAIL (as scored; diagnosed below)** |
| **P2** | `v_q(delta_q^{(u)}(z)) == m - u` in 100% of checks; symmetric CRT lift `Y_1` mod `M = 9,699,690` has `min \|Y_1\| > 256` and median `\|Y_1\|/M in [0.15, 0.35]` | `3360 / 3360` exact valuations; `min \|Y_1\| = 30,907 >> 256`; median `\|Y_1\|/M = 0.2516` | **PASS** |
| **P3a** | Over `K = Q(zeta_5)` (`m = 4`), `11` primes `103..157`, `R = 0.2191` (`sqrt(R) = 0.4681`), `5/11` primes agreeing (`A/W <= 0.4630 < sqrt(R)`): `s = 3` Galois-folded LLL hits in `>= 80%` of 30 trials vs `0%` for `s = 1` | `s = 3` hit **30/30 (100%)** (and intersecting the first 4 LLL basis vectors uniquely pinpointed `alpha` in **30/30**!), vs **0/30** for `s = 1` | **PASS** |
| **P3b** | 100% of primitive `163`-smooth elements in `[-6, 6]^4 \ {0}` of `Z[zeta_5]` have full-orbit agreement `A_orbit == 0` | `11,456` primitive `163`-smooth elements found; **`0 / 11,456` (`0%`)** have `A_orbit > 0` | **PASS** |
| **P4** | On 40 planted 28-bit `N = p q` smooth-norm windows (`k^2 + N` is `47`-smooth, `H = 512`), full-modulus `Z[k]`, `Int(Z)`, and linear `delta_2` recover `0/40` | `ln val = 19.18 << sqrt(ln P * ln H^2) = 58.49`: `0/40` (`Z[k]`), `0/40` (`Int(Z)`), `0/40` (linear `delta_2`) | **PASS** |

### 6.2 Disclosed Post-Hoc Diagnosis of P1c (`POST-HOC`, `posthoc_diag.py`)

Why did `P1c` give `12.84 < 13.85` at `(W' = 21.90, b = 6.93)`? Printing the shortest LLL vectors in `posthoc_diag.py` revealed two instructive facts that actually **reinforce** Theorem 1c:
1. **Degenerate `c_3 = c_4 = 0` hyperplane collapse:** At `W' = 21.90, b = 6.93`, the full-rank dimension-`s+1` Minkowski threshold `(W' + sum ln M_u)/(s+1)` increases once `ln M_u > 14.23` (`s >= 3`), so LLL sets **`c_3 = c_4 = 0` identically** in *both* lattices (e.g. `j_0 = 868`: `c_Buium = [-43834, 484, 1, 0, 0]`, `c_Copp = [-87668, 969, -1, 0, 0]`). Neither `delta_2^{(2)}` (`deg 4`) nor `delta_2^{(3)}` (`deg 8`) was used at all! Within `s = 2`, `(1, Y, (Y - Y^2)/2)` saved a single constant `ln 2 = 0.69` nat solely because `q = 2` was omitted from `P_odd` (when `2 | P` in `E1b`, that `ln 2` is immediately reclaimed by the `q = 2` Lucas index `2^{floor(L/2)+1}`, giving `11.55 >= 11.20`), amplified by uniform `j_0 in [3, H)` where `(j_0 / H)^2 < (j_0 / H)`.
2. **Active `s = 3, 4` test (`W' = 69.55, b = 5.55`, `j_0 in [3H/4, H)`):** When `W' / b` is large enough that `c_3, c_4 != 0` actually activate (`P_odd = prod_{3<=q<=79} q`, `H = 2^8`), Theorem 1c's ordering holds cleanly:
   - `s = 2`: `mean Buium = 27.685` vs `mean Copp_r1 = 27.668`;
   - `s = 3` (`c_3 != 0` active): **`mean Buium = 25.901 > mean Copp_r1 = 24.973`**;
   - `s = 4` (`c_4 != 0` active for Copp, while Buium still zeroes `c_4 = 0` because `deg 8` is too costly): **`mean Buium = 25.963 > mean Copp_r1 = 24.765`**.
