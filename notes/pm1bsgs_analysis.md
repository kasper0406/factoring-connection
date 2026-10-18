# pm1bsgs: Pollard p−1 stage 1 followed by a BSGS search for b^N ≡ b^x (mod p)

Agent `pm1bsgs`, 2026-10-18. Idea source: `notes/user_ideas/2026-10-18_pm1_frobenius_bsgs.md`.
Code: `src/pm1bsgs/` (operational code in `operational.py`, cost formulas in `costmodel.py`, factor-aware harness in `harness.py`). Tests: `tests/test_pm1bsgs.py` (17 pass). Raw data: `experiments/raw/pm1bsgs_*`. Derived data: `experiments/analysis/pm1bsgs_*.json`.

---

## 0. Verdict for the user (plain language)

**What is right.** Every step of the idea is mathematically correct, once two small points are made precise.
- ord_p(b) divides R, provided S is defined as gcd(p−1, E) (§2, L1).
- N ≡ q (mod R) (§2, L2).
- A baby-step/giant-step search over a window of length X finds p whenever ord_p(b) ≤ X, at a cost of about √X operations up to log factors. This holds whether R is prime or composite (§2, L3–L5).

The idea also beats the classical prime-only stage 2. That stage covers primes up to B2 at a cost of about B2 operations, while a BSGS window covers them at about √B2. This is a real, polynomial-size improvement in how the cost grows with the required success rate (§6.2). Your instinct that stage 2 should be a square-root search, not a prime-by-prime scan, is correct and important.

**What is already known.** This continuation is, in substance, the "FFT continuation" of p−1:
- Pollard 1974 sketched it at the end of his p−1 paper.
- Montgomery–Silverman 1990 made it practical.
- Montgomery–Kruppa 2008 and Kruppa's 2010 thesis refined it.
- For ECM, Brent 1986 gave the birthday-paradox version.

Those algorithms already test every integer in the window that is coprime to a fixed modulus P, composites included. Kruppa's thesis says so explicitly (§7). Two of your points are not new advantages:
- **Using N.** The congruence b^N ≡ b^x (mod p) just says that ord_p(b) divides N − x. So N only fixes where the search window sits. Searching for a multiple of ord_p(b) in [1, X] (plain order-finding) works equally well. So does a window around a random number T, which we checked by experiment (§5.3).
- **Composite cofactors.** Covering composite R never helps at the best parameter choice. Raising B1 to √X costs no more than the stage-2 search itself. Once B1 ≥ √X, any cofactor R ≤ X with no prime factor ≤ B1 must be 1 or a prime (§3, Prop. D). In our budget optimisation, every best parameter choice had B1² ≥ X: 73 of 73 levels in the measured model.

**True complexity.**
- **One group (p−1 only).** New random bases do not help: p−1 does not change, and in our data 100% of instances either succeed for every base or fail for every base. The expected cost for a balanced N is **N^{1/4}·Θ(1/log² N) times the stage-2 polylog factor**. That is **N^{1/4−o(1)}**: a gain of Θ(log² N) over plain BSGS, and not polynomial or subexponential.
- **For a constant fraction of N**, the cost is a smaller power of N. For example, half of all balanced N cost N^{0.151+o(1)} (§6.2). This is the known p−1 success profile.
- **Many groups (elliptic curves, i.e. ECM with this stage 2).** The cost stays L_p[1/2, √2]. The stage 2 only improves the o(1) term, by a factor polynomial in log p, as Brent 1986 and Brent–Kruppa–Zimmermann 2017 already state (§6.3).

**Measured.** On balanced semiprimes of 32–80 bits (500 holdout instances per size), the hybrid's expected cost to factor at 80 bits is about 6.6·10^6 ring operations in our tree-based implementation. Plain BSGS needs 3.1·10^8 and Pollard rho 2.1·10^6. The hybrid's *median* (4.6·10^5) is 4.7× below rho's (§5.4).

---

## 1. Setup and notation

N = pq with p ≠ q odd primes. a is a random unit mod N, B is the stage-1 bound, and E = E_B = ∏_{ℓ≤B} ℓ^{e_ℓ} with e_ℓ = ⌊log_ℓ B⌋. Stage 1 computes b = a^E.

S_p := gcd(p−1, E) and R_p := (p−1)/S_p. o_p := ord_p(a) and d_p := ord_p(b). The same notation is used for q.

"BSGS over window W" means: find t ∈ W with d_p | t by comparing b^{u_i} and b^{v_j} through gcd(∏(b^{u_i} − b^{v_j}), N). We never hash: the values mod p are hidden.

---

## 2. Correctness (each item: status, then proof)

**L1 (proved).** d_p = o_p / gcd(o_p, E), and d_p | R_p. Moreover, b is uniform on the unique subgroup of order R_p of F_p^*, so P(d_p = d) = φ(d)/R_p for each d | R_p.

*Proof.* (a^E)^k = 1 ⇔ o_p | Ek ⇔ o_p/gcd(o_p,E) | k. For each prime ℓ, v_ℓ(d_p) = max(0, v_ℓ(o_p) − e_ℓ) ≤ max(0, v_ℓ(p−1) − e_ℓ) = v_ℓ(R_p). The map x ↦ x^E on a cyclic group of order p−1 has image of order (p−1)/gcd(p−1,E) = R_p, and every fibre has the same size. ∎

*Precision note.* L1 is false if "S" means the full B-smooth part of p−1 and some ℓ^{v_ℓ(p−1)} > B. Example: p = 17, B = 2, E = 2. The full 2-smooth part is 16, so R = 1, but b = a² can have order 8. The user's note does say "S | E", which is the right condition. Test: `test_lemma1_order_divides_R` (exhaustive for p < 400).

**L2 (proved).** N ≡ q (mod p−1), hence modulo every R | p−1. So a^N ≡ a^q (mod p) for every a: since p ≡ 1 (mod p−1), we get N = pq ≡ q. Test: `test_lemma2_N_congruent_q`.

**L3 (proved: the equivalence that removes N).** b^N ≡ b^x (mod p) ⇔ d_p | N − x.

So the search "find x ∈ [0,X] with b^N ≡ b^x" is the search "find a multiple of d_p in [N−X, N]". Note that x = N is always a (useless) solution: the task is to find a *second* exponent with the same residue mod d_p, i.e. a nonzero multiple of d_p. That is order finding. Any interval of length ≥ d_p contains a multiple of d_p, so the windows [N−X, N], [1, X] and [T−X, T] (T arbitrary) all carry the same guarantee "d_p ≤ X+1 ⇒ hit". They also cost the same (one extra exponentiation b^T). Test: `test_lemma3_equivalence_with_order`.

*What N does add (experimentally supported, §5.3).* When X < d_p, the window [N−X, N] still contains a multiple of d_p with probability about X/d_p, because N mod d_p is "random". [1, X] never does. This makes the hybrid's mode "N" roughly 1.6–2× cheaper in expected cost than the unshifted order-finding window. A random shift T gives the same effect: a_shift ≈ a_N at every n. And a random T is re-drawn on every restart, whereas N mod d_p is fixed per instance. So N contributes nothing that a random number does not.

**L4 (proved: gcd detection).** For integers u ≠ v, p | b^u − b^v ⇔ d_p | u − v, and q | b^u − b^v ⇔ d_q | u − v. So a single test gcd(b^N − b^x, N) returns N exactly when lcm(d_p, d_q) | N − x.

For the user's specific test x* = N mod d_p, the failure probability is

  P_fail = Σ_{d_p | R_p, d_q | R_q} φ(d_p)φ(d_q)/(R_p R_q) · [d_q | N − (N mod d_p)]

over a uniform base. Conditioned on d_p, d_q > 1 (stage 1 did not finish either prime), the denominators become (R_p − 1)(R_q − 1). The formula is implemented as `harness.user_first_x_failure_prob`.

Measured over 200 instances per n at B1 = 2^10 (`experiments/analysis/pm1bsgs_failure.json`): P_fail is 0 for almost every instance. There was one exception at n = 56. There R_p = 884171 is prime and R_q = 8, and ⌊N/R_p⌋ happens to be ≡ 0 (mod 8), so the single test fails for *every* base. The failure condition is a congruence, not only "d_p = d_q".

**L5 (proved: full-backtracking criterion).** The operational algorithm takes gcds per block, then per point, then per pair. It returns a proper factor iff the window W contains some t divisible by exactly one of d_p, d_q (finding q also factors N).

Checking the first two multiples of d_p in W and the first two of d_q decides this. *Proof:*
- If d_p = d_q, no pure hit exists.
- If d_q ∤ d_p, two consecutive multiples of d_p cannot both be divisible by d_q, since d_q would then divide their difference d_p.
- If d_q | d_p and d_q ≠ d_p, then d_p ∤ d_q, and the same argument applies with p and q swapped.
- If W holds fewer than two multiples, the ones present are checked directly. ∎

Consequently, once |W| ≥ 2·max(d_p, d_q), the algorithm fails iff d_p = d_q. That happens with probability Σ_{d | gcd(R_p,R_q)} φ(d)²/(R_p R_q) (`harness.prob_equal_orders`, test `test_equal_order_probability`). Conditioned on both orders > 1, it is 0 for all 200 sampled instances at every n (gcd(R_p, R_q) = 1 for random balanced p, q is the typical case).

The stage-1 analogue is "both orders are completed by the same prime-power step", which gives gcd = N. Its exact mean probability is 9.7·10^-3 at n = 32, 4.5·10^-5 at n = 40 and < 2·10^-7 for n ≥ 48.

Tests: `test_interval_pure_hit_bruteforce`, and `test_exact_probs_exhaustive`, which enumerates every base for N = 103·107, 211·223 and 157·181 and matches the exact outcome distribution (Fractions) to the operational runs.

**L6 (proved: cost).** For a window of size X = mJ:
- The hash-free BSGS with the plain product of all differences costs mJ = X multiplications (`naive=True`), which is linear.
- With subproduct and remainder trees it costs O(M(m) log m · J/m), i.e. Õ(√X) for m ≈ √X (Pollard 1974 / Strassen 1976 style).
- Both the baby set and the giant set are geometric progressions. So the root-polynomial build and the evaluation can each be done in O(M(m)) (chirp-z / Bostan–Schost). BKZ 2017 states O(M(d)) for p−1 (§7). We did not implement this; it is charged as the assumed-constant variant `a_N_gp`.

The user's cost claim "O(N^{1/4}/√Σ) when S ≥ Σ" is correct up to this polylog factor and the factor 2 for "p or q".

---

## 3. Why composite cofactors do not matter at the optimum

**Prop. D (proved).** Fix X and let B ≤ √X. Let H(B, X) be the hybrid: stage 1 to B, then a BSGS window of length X. Let H(B′, X) be the same with B′ = ⌈√X⌉.
1. Every (p, a) on which H(B, X) finds p is also handled by H(B′, X). Either stage 1 finds a factor, or d′_p | d_p ≤ X and the window still hits. Separation failures (d_p = d_q) are excluded here, since L5 shows they are negligible.
2. The cost of H(B′, X) is at most cost(H(B, X)) + c₁√X, where c₁ ≈ 2.2 multiplications per unit of B.
3. Under H(B′, X), the B′-rough part of R_p is 1 or a single prime: a product of two primes > √X already exceeds X. R_p can additionally contain excess powers ℓ^{v−e_ℓ} of primes ℓ ≤ B′.

So the "regardless of whether R is composite" feature buys at most a factor 2 in cost. It is never needed when stage 1 is pushed to √X, and pushing it there costs no more than the stage-2 search itself.

The one genuine extra over a prime-only stage 2 is R = ℓ^k·r, where ℓ ≤ B has excess valuation. Its probability per prime is about Σ_ℓ 1/(φ(ℓ)ℓ^{e_ℓ}) = O(1/B^{1/2}), which is negligible.

**Experimentally supported.** For each budget 2^c, we took the parameter choice that is best on the training split. At every level with train success strictly between 2% and 99%, B1² ≥ X held:
- `a_N`: 73 of 73 levels;
- `a_shift`: 76 of 76;
- `a_N_gp`: 60 of 62. The two exceptions sit at the saturation edge, with train success 0.96–0.97.

Composite coverage only enters at *saturation*, where one insists on success for every instance. There the cheapest choice is a small B1 with X ≈ p/S_min. That is the user's regime, and its cost is ≈ √(p/2), i.e. plain BSGS up to a constant.

---

## 4. Implementation (operational: sees only N)

`src/pm1bsgs/operational.py`:
- **Base and stage 1.** `choose_base(N, seed)` with gcd(a, N) checked; an ACCIDENTAL factor is recorded if the gcd is > 1. `stage1` computes b = a^E prime power by prime power, takes gcd(b−1, N) every 64 steps, and replays step by step on gcd = N.
- **(a) `pm1_bsgs(..., mode="N")`, the user's hybrid.** Baby values u_i = b^{N+i} for i < m. Giant values v_j = b^{jm} for j < J. The window is [N−(J−1)m, N+m−1]. F(T) = ∏(T−u_i) comes from a subproduct tree. ∏_j F(v_j) is computed by a remainder tree using only flint `*` and `%` by monic polynomials, which is the audited subset. gcd is taken per block; on gcd = N we backtrack per point, then per pair. `naive=True` gives the hash-free product of all differences (linear cost).
  - `mode="order"`: u_i = b^i, v_j = b^{jm} (j = 1..J), window [1, Jm].
  - `target=T`: the random-shift control.
- **(b) `pm1_prime`.** Classical stage 2: ∏_{B1<r≤B2} (b^r − 1) with a prime-gap table, 2 multiplications per prime, gcd per 256 primes, and per-prime replay.
- **(c) `pm1_birthday`.** Brent 1986's walk x_1 = b, x_{j+1} = x_j² b^{c_j}, with public random bits c_j (value-independent). All K(K−1)/2 pairs are tested at once via ∏_i F′(x_i) = ±∏_{i<j}(x_i−x_j)², which is Brent's §7 identity. The derivative is built by hand.
- **(d) Plain BSGS**: `pm1_bsgs` with B1 = 1.

**Cost unit `eq_ops`** = scalar mults mod N + gcds + polynomial work. Polynomial work uses the model M(d) = d⌈log₂ d⌉ per product and 2M(·) per remainder; a power x^e is charged the binary-method count. `costmodel.py` reproduces the instrumented counts exactly (tests) for runs without early exit or backtracking.

Wall-time calibration (`experiments/raw/pm1bsgs_walltime.jsonl`): the tree BSGS runs at 49–50 ns per eq_op for m ≥ 4096, and the Python scalar modmul loop at 70 ns per op. So in this implementation eq_ops track wall time to within a factor of 1.5 across methods.

**Validation.** `experiments.py validate` ran 2800 operational runs at n = 32..80 with random parameters and bases:

| Method | Runs | Factors found | Of which runs with gcd = N backtracking |
|---|---|---|---|
| bsgs | 1680 | 926 | 298 |
| prime | 560 | 261 | 25 |
| birthday | 560 | 330 | 110 |

**All 2800 runs agree with the exact harness predictor (0 disagreements).** All 1517 factors were STRUCTURAL and 0 were ACCIDENTAL, against a null expectation for gcd(a, N) > 1 of about 2^{-n/2+1} per run. That validation is what allows the large-sample statistics below to use the predictor.

---

## 5. Measurements

The instance family is `balanced_semiprime(n, seed)` for n ∈ {32, 40, 48, 56, 64, 72, 80}, with 1000 seeds per n.
- Seeds 0–499 are the training split, used to pick the best parameter at each budget. Seeds 500–999 are the holdout split, on which everything below is reported.
- Exact outcomes come from the validated predictor, using the factorisations of p−1 and q−1 (harness only).
- Rho uses `src/common/baselines.pollard_rho_brent` on the holdout seeds and is actually run.

### 5.1 The smooth part S_B(p−1) of random balanced p (`pm1bsgs_smooth.json`, `pm1bsgs_asymptotics.json`)

There are 8000 primes per n, using both p and q.

- **E[S_B^{-1/2}]** (the quantity that controls the expected stage-2 window), at n = 80 (p ≈ 2^40):

  | B | Measured | Fixed-B limit Π_B (see §6.1) | Measured × ln B |
  |---|---|---|---|
  | 2^8 | 0.158 | 0.151 | 0.88 |
  | 2^12 | 0.112 | 0.103 | 0.93 |
  | 2^16 | 0.085 | 0.077 | 0.94 |

  Π_B·ln B → 0.857. **E[S^{-1/2}] ≈ 0.86–0.94 / ln B** (experimentally supported).
- **(p−1)/2 prime** (the hardest case: S_B = 2 for every B < (p−1)/2). Measured fraction against the Hardy–Littlewood prediction C₂/ln p (C₂ = 0.6602):

  | n | 32 | 40 | 48 | 56 | 64 | 72 | 80 |
  |---|---|---|---|---|---|---|---|
  | Measured | 0.056 | 0.047 | 0.041 | 0.033 | 0.035 | 0.030 | 0.025 |
  | C₂/ln p | 0.060 | 0.048 | 0.040 | 0.034 | 0.030 | 0.027 | 0.024 |

- **Semiprime version**, E[max(S_p, S_q)^{-1/2}] at n = 80: 0.165 (B = 2^4), 0.055 (2^8), 0.024 (2^12), 0.013 (2^16), 0.010 (2^20). This decays like 1/ln² B, as §6.1 predicts.
- The median of log P₁(p−1)/log p is 0.55–0.59.

### 5.2 Success probability vs. budget, at equal eq_ops (holdout; one base per instance)

Each column is the parameter that is best on the training split for that method at budget ≤ 2^c. `a_N` is the user's hybrid (stage 1 + BSGS on b^N ≡ b^x). `a_order` is plain order-finding BSGS after stage 1. `a_shift` uses a window at a random T. `a_N_gp` is `a_N` charged under the geometric-progression cost model, with assumed constants (2M build, 3M per evaluation block). `b_prime` is classical p−1 with prime-only stage 2. `c_bday` is the Brent birthday stage 2. `d_plain` is BSGS without stage 1. Rho is the fraction of actual runs finishing within the budget.

**n = 64**

| Budget | a_N | a_order | a_shift | a_N_gp | b_prime | c_bday | d_plain | rho |
|---|---|---|---|---|---|---|---|---|
| 2^12 | 0.20 | 0.11 | 0.18 | 0.22 | 0.26 | 0.11 | 0.00 | 0.00 |
| 2^14 | 0.46 | 0.36 | 0.47 | 0.55 | 0.47 | 0.32 | 0.00 | 0.00 |
| 2^16 | 0.66 | 0.57 | 0.66 | 0.80 | 0.66 | 0.56 | 0.00 | 0.14 |
| 2^18 | 0.84 | 0.79 | 0.84 | 0.93 | 0.79 | 0.73 | 0.00 | 0.89 |
| 2^20 | 0.95 | 0.92 | 0.95 | 0.98 | 0.87 | 0.89 | 0.06 | 1.00 |
| 2^22 | 0.98 | 0.97 | 0.98 | 1.00 | 0.94 | 0.98 | 0.35 | 1.00 |
| 2^24 | 1.00 | 1.00 | 1.00 | 1.00 | 0.97 | 0.99 | 1.00 | 1.00 |

**n = 80**

| Budget | a_N | a_order | a_shift | a_N_gp | b_prime | c_bday | d_plain | rho |
|---|---|---|---|---|---|---|---|---|
| 2^14 | 0.13 | 0.10 | 0.14 | 0.19 | 0.15 | 0.11 | 0.00 | 0.00 |
| 2^16 | 0.32 | 0.27 | 0.32 | 0.44 | 0.33 | 0.25 | 0.00 | 0.00 |
| 2^18 | 0.56 | 0.47 | 0.55 | 0.71 | 0.47 | 0.42 | 0.00 | 0.01 |
| 2^20 | 0.76 | 0.70 | 0.76 | 0.88 | 0.65 | 0.66 | 0.00 | 0.11 |
| 2^22 | 0.88 | 0.83 | 0.87 | 0.95 | 0.77 | 0.83 | 0.01 | 0.88 |
| 2^24 | 0.97 | 0.95 | 0.97 | 0.99 | 0.86 | 0.90 | 0.02 | 1.00 |
| 2^26 | 0.99 | 0.99 | 0.99 | 1.00 | 0.92 | 0.97 | 0.23 | 1.00 |
| 2^28 | 1.00 | 1.00 | 1.00 | 1.00 | 0.95 | 0.99 | 0.87 | 1.00 |

n = 48 and the full curves are in `pm1bsgs_budget.json`.

Reading (experimentally supported):
- **Low budgets.** Below about √p, (a) ≈ (b): the cheapest wins come from p−1 being smooth, and both methods get them.
- **High budgets.** (a) pulls clearly ahead of (b), because BSGS covers a prime R ≤ X at √X instead of X/ln X.
- **Birthday.** (c) tracks (a) at a modestly higher cost, consistent with its random rather than structured coverage.
- **Plain BSGS.** (d) is useless until its window reaches p.
- **Rho.** Rho has a sharp threshold near √p ≈ 2^20 (n = 80).

### 5.3 What N contributes (a_N vs a_order vs a_shift)

Mean expected cost with restarts (§5.4) for `a_N` / `a_shift` / `a_order`:

| n | 48 | 64 | 80 |
|---|---|---|---|
| a_N | 3.6e4 | 7.2e5 | 6.6e6 |
| a_shift | 3.1e4 | 7.3e5 | 7.1e6 |
| a_order | 6.5e4 | 1.4e6 | 1.2e7 |

The random shift reproduces the N-window. The unshifted window is about 2× worse, which is the "lottery" effect of L3. Earlier, a training-saturation artifact made `a_N` fail on one instance, a Sophie-Germain-type instance with p−1 = 4r and q−1 = 2r′. Because N mod d_p is fixed, every restart lands on the same miss. We fixed this by preferring larger coverage at saturated levels. The shifted window does not have this failure mode.

### 5.4 Expected total cost to factor, with restarts (holdout, 500 instances per n)

The schedule is: at budget levels 2^8, 2^9, …, 2^36, use the train-best parameter for that level with a **new base** each time, and sum all attempts until success. A cell reads mean / median eq_ops; "cens" is the fraction not factored by the 2^36 level, counted at the cost spent so far (a lower bound).

| n | a_N (user) | a_order | a_N_gp | b_prime | c_bday | d_plain | rho (measured) |
|---|---|---|---|---|---|---|---|
| 32 | 1.75e3 / 716 | 3.25e3 / 1.49e3 | 1.15e3 / 710 | 1.09e3 / 215 | 3.6e3 / 1.55e3 | 9.09e3 / 7.19e3 | 659 / 773 |
| 40 | 7.84e3 / 3.37e3 | 1.55e4 / 3.49e3 | 4.36e3 / 1.57e3 | 7.02e3 / 1.47e3 | 1.70e4 / 6.37e3 | 6.92e4 / 5.27e4 | 2.18e3 / 2.31e3 |
| 48 | 3.64e4 / 7.29e3 | 6.54e4 / 1.35e4 | 1.57e4 / 6.95e3 | 5.10e4 / 7.41e3 | 7.04e4 / 1.42e4 | 4.04e5 / 3.92e5 | 8.09e3 / 8.47e3 |
| 56 | 1.24e5 / 2.85e4 | 2.12e5 / 2.84e4 | 4.45e4 / 1.41e4 | 3.32e5 / 1.49e4 | 2.18e5 / 3.04e4 | 2.39e6 / 1.72e6 | 3.36e4 / 3.39e4 |
| 64 | 7.20e5 / 5.99e4 | 1.43e6 / 1.19e5 | 2.13e5 / 2.91e4 | 6.79e6 / 6.17e4 | 1.68e6 / 1.13e5 | 1.22e7 / 1.25e7 | 1.31e5 / 1.36e5 |
| 72 | 2.51e6 / 2.22e5 | 3.89e6 / 2.21e5 | 5.62e5 / 1.14e5 | ≥4.05e8 / 2.41e5 (2.8% cens) | 4.75e6 / 2.43e5 | 6.15e7 / 5.01e7 | 5.19e5 / 5.41e5 |
| 80 | 6.55e6 / 4.57e5 | 1.20e7 / 9.24e5 | 1.80e6 / 2.36e5 | ≥7.86e8 / 9.03e5 (5.4% cens) | 1.19e7 / 9.38e5 | 3.12e8 / 2.40e8 | 2.10e6 / 2.16e6 |

**New bases do not help (experimentally supported, `pm1bsgs_basedep.json`).** At the train-best parameters for budgets 2^18 and 2^22, for both a_N and b_prime and n ∈ {48, 64, 80}, **100%** of holdout instances either succeed for all 10 bases or fail for all 10. The instance is either favourable or not. All of the restart schedule's benefit comes from the growing budget.

**Scaling fits (diagnostic only).** Slope of log₂(cost) vs n, in "N^slope" units; theory says every method tends to 1/4 up to log factors, except b_prime.

| Method | Mean, all n | Mean, n ≥ 56 | Median |
|---|---|---|---|
| a_N | 0.253 | 0.237 | 0.192 |
| a_shift | 0.252 | 0.245 | 0.202 |
| a_order | 0.250 | 0.236 | 0.192 |
| a_N_gp | 0.222 | 0.218 | 0.177 |
| c_bday | 0.250 | 0.235 | 0.184 |
| d_plain | 0.311 | 0.293 | 0.312 |
| rho | 0.244 | 0.248 | 0.241 |
| b_prime | 0.43 | 0.49 | 0.241 |

b_prime's means are censored lower bounds; its mean slope near 1/2 is the p−1 worst case P₁(p−1) ~ p. **Ratio a_N/d_plain** (the gain over plain BSGS) is 5.2, 8.8, 11, 19, 17, 25 and 48 for n = 32…80, growing roughly like log² p, as §6.1 predicts.

The slopes below 1/4 at these sizes are log-factor effects (Rule 8: they are not evidence of a smaller exponent). The medians are genuinely smaller powers (§6.2).

---

## 6. Asymptotics

### 6.1 Single group, expected cost

**(proved; fixed B, p → ∞).** S_B(p−1) = gcd(p−1, E_B) depends only on p mod E_B. By the prime number theorem in arithmetic progressions, p mod E_B is equidistributed over units. So

  lim_{p→∞} E[S_B^{-s}] = Π_B(s) := ∏_{ℓ≤B} Σ_{k=0}^{e_ℓ} P_ℓ(k) ℓ^{-ks}, with P(ℓ^k | p−1) = 1/φ(ℓ^k).

Mertens' theorem gives Π_B(s) ≍ 1/ln B for every fixed s > 0, since the factor at ℓ is 1 − 1/ℓ + O(ℓ^{-1-s}). Numerically, Π_B(1/2)·ln B → 0.857; the measured values are in §5.1.

**(conjectured; B = p^θ, needs uniformity beyond PNT-in-AP).** E[S_B^{-1/2}] ≍ 1/ln B continues to hold for B up to p^{1/2}. A Brun–Titchmarsh/Selberg sieve upper bound plausibly gives the "≪" direction; we did not check this (UNVERIFIED). The lower bound "≫" follows from P((p−1)/2 prime) ~ C₂/ln p, which is the Hardy–Littlewood conjecture (measured match in §5.1).

**Consequence (conjectured, under the above).**
- *One prime.* The expected cost of the best (B, doubling-X) strategy is Θ̃(√p / ln p) times the stage-2 polylog.
- *Balanced N = pq*, where p and q are independent so either may be found: E[min(√R_p, √R_q)] ≍ √p · E[max(S_p,S_q)^{-1/2}] ≍ √p/ln² B. The bound uses max(S_p,S_q)^{-1/2} ≤ S_p^{-1/4}S_q^{-1/4} for "≤" and P(S_p = S_q = 2) for "≥". Hence

  **E[cost] = Θ(N^{1/4} / log² N) × (stage-2 polylog)** = N^{1/4−o(1)} with o(1) = Θ(log log N / log N).

  The stage-2 polylog is log N for geometric-progression/FFT evaluation and log² N for generic trees. Plain BSGS has E[cost] = Θ(N^{1/4}) × the same polylog. **The gain is exactly a Θ(log² N) factor.** Against rho (Θ(N^{1/4}) multiplications, heuristic), the hybrid gains about log N in the GP/FFT model and is roughly level in the tree model. Measured at n = 80: a_N_gp/rho = 0.86 and a_N/rho = 3.1.
- **Not polynomial, not subexponential.** Under the Hardy–Littlewood conjecture, a fraction ≍ 1/log² N of balanced N has both (p−1)/2 and (q−1)/2 prime. In the generic-group model, such N need Ω(√p) operations for this method family (order finding for a prime-order element; Shoup/Nechaev-type bound, cited from memory, UNVERIFIED).

### 6.2 Single group, cost at a fixed success rate (polynomial gain on a constant fraction)

**Model.** The B-smooth part of an integer of size B^w has log S/log B with density ρ(u)ω(w−u) on [0, w−1], plus an atom ρ(w) at w (ρ is Dickman, ω is Buchstab). The total mass checks numerically to 1 − 10^{-7}. As w → ∞ this tends to e^{-γ}ρ(u); Pinsky 2019, Thm 1 gives the e^{-γ}ρ law for the closely related p_N-smooth measure (VERIFIED by the literature subagent). Applying the model to p−1 for prime p is a heuristic (H–T 1993, p. 468: only positive-proportion results are proven).

**Optimisation.** Write cost = max(B, √X) and drop polylogs. A grid scan over (β, ξ), with B = p^β and X = p^ξ, puts the optimum at ξ = 2β. That is Prop. D again: X = B², so R is prime. The cost exponent c(δ) (cost = p^{c}, i.e. N^{c/2}) at success rate δ:

| δ | 0.05 | 0.10 | 0.25 | 0.50 | 0.75 | 0.90 | 0.99 |
|---|---|---|---|---|---|---|---|
| BSGS/FFT stage 2, one prime | 0.230 | 0.251 | 0.292 | 0.344 | 0.406 | 0.456 | 0.495 |
| BSGS/FFT stage 2, semiprime | 0.214 | 0.231 | 0.262 | **0.301** | 0.344 | 0.387 | 0.456 |
| prime-only stage 2, one prime | 0.335 | 0.379 | 0.471 | 0.607 | 0.779 | 0.905 | 0.99 |
| prime-only stage 2, semiprime | 0.303 | 0.336 | 0.402 | 0.493 | 0.607 | 0.729 | 0.905 |
| plain BSGS / rho | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 |

Status: conjectured, i.e. derived from the ρ·ω model. So **half of balanced N cost about N^{0.151}** with a square-root stage 2, against about N^{0.247} with a prime-only stage 2. This polynomial advantage of the √B2 continuation over the prime-only one is the substantive, correct content of the user's idea. It is the classical Montgomery–Silverman gain, now quantified as an exponent.

### 6.3 Re-rolled groups (ECM with the same BSGS stage 2)

**Status: conjectured**, i.e. derived under the standard ECM heuristic that #E(F_p) behaves like a random integer in the Hasse interval for smoothness. That is the same status as ECM's own L[1/2, √2] bound.

Per curve, take stage-1 bound B and window X = B^k, with cost ≈ B + Õ(B^{k/2}). Let u = ln p/ln B. Success means the B-rough cofactor of #E is ≤ B^k, which has probability between ρ(u) and ρ(u−k)·O(1). Since ρ(u−k)/ρ(u) = u^{k(1+o(1))}, the window multiplies the success probability by only (log p)^{O(k)}.

Optimising B = L_p[1/2, 1/√2] gives **L_p[1/2, √2 + o(1)] for every fixed k.** If X instead grows polynomially in p, the stage-2 cost p^{ξ/2} is exponential in log p and dominates.

So the BSGS stage 2 changes ECM only inside the o(1). This agrees with Brent 1986 ("speedup is of order log(p) … not very impressive") and BKZ 2017 ("can be absorbed into the o(1) term"), both VERIFIED quotes (§7).

---

## 7. Prior art

The search used WebSearch/WebFetch through a literature subagent. AMS and HAL PDFs were blocked (403).

| Source | Status | What it says |
|---|---|---|
| J. M. Pollard, "Theorems on factorization and primality testing", Math. Proc. Camb. Phil. Soc. 76(3) (1974) 521–528, doi:10.1017/S0305004100049252 | bib VERIFIED; content PARTIAL | p−1 is in §4. Its stage 2 tests each prime q ∈ (B1, B2] at 2 mults per prime (via Kruppa thesis Ch. 2). The paper ends by suggesting an FFT continuation (Kruppa). Its deterministic O(N^{1/4+ε}) factoring evaluates products along geometric progressions with chirp-z (Bernstein–De Feo–Leroux–Smith 2020, §3). **To check:** Pollard's own wording and theorem numbers. |
| V. Strassen, "Einige Resultate über Berechnungskomplexität", Jber. DMV 78 (1976/77) 1–8 | bib via secondary citations only (Costa–Harvey 2014; Hittmeir 2016) | Streamlined version of Pollard's N^{1/4+ε}. **To check:** zbMATH/primary record. |
| R. P. Brent, "Some integer factorization algorithms using elliptic curves", Australian Computer Science Communications 8 (1986) 149–163; arXiv:1004.3366 | VERIFIED (full text) | §6 birthday-paradox phase 2: walk Q_{j+1} = Q_j² or Q_j²·Q (multiplicative notation), collisions via ∏(x_i − x_j) and gcd, success ≈ 1 − exp(−r²/n₁). §7: D² = ∏P′(x_j), no FFT. §9.2: "can be used as a second phase for Pollard's 'p−1' algorithm". Speedup "of order log(p)". Analysis assumes the cofactor n₁ is prime. **Our method (c) implements exactly this walk and identity.** |
| P. L. Montgomery, "Speeding the Pollard and elliptic curve methods of factorization", Math. Comp. 48 (1987) 243–264, doi:10.1090/S0025-5718-1987-0866113-7 | bib VERIFIED; text not accessed | BSGS pairing stage 2, still O(π(B2) − π(B1)) operations (Kruppa's description). |
| P. L. Montgomery, R. D. Silverman, "An FFT extension to the P−1 factoring algorithm", Math. Comp. 54 (1990) 839–854, doi:10.1090/S0025-5718-1990-1011444-3 | bib + abstract VERIFIED; full text UNVERIFIED | Evaluates a polynomial along a geometric progression by circular convolution in O(√B2 log B2) (Kruppa §2.2). **To check:** the exact exponent set covered, and any remark on composite cofactors. |
| P. L. Montgomery, A. Kruppa, "Improved stage 2 to P±1 factoring algorithms", ANTS VIII, LNCS 5011 (2008) 180–195; A. Kruppa, PhD thesis, Nancy 2010 | VERIFIED (thesis PDF) | Stage 2 "checks s1(ℓmax−s1) candidates, not necessarily prime but including all primes in ]B1,B2]" (§2.4). Covers "all values coprime to P, in particular all primes" (§2.6). The correctness argument assumes only b1^q ≡ 1 with gcd(q, 2P) = 1, so **composite cofactors coprime to 2P in the window are detected**. The probability model (§5.3.3) counts only prime cofactors. |
| P. Zimmermann, B. Dodson, "20 years of ECM", ANTS VII, LNCS 4076 (2006) 525–542 | VERIFIED (author PDF) | §3.1: the birthday continuation hits "most primes in [B1,B2], and usually other larger primes"; the standard continuation guarantees all primes in [B1,B2], for O(√B2) curve operations. §3.2: Eq. (4) "computes many x_σ − x_τ that do not correspond to prime values of σ−τ". |
| R. P. Brent, A. Kruppa, P. Zimmermann, "FFT extension for algebraic-group factorization algorithms", Ch. 8 in *Topics in Computational Number Theory Inspired by Peter L. Montgomery* (CUP 2017); maths-people.anu.edu.au/brent/pd/rpb264.pdf | VERIFIED | Improvements "can be absorbed into the o(1) term" of L(p)^{√2+o(1)}. §8.2.4: "We also assume o is prime; if it is composite, then it is very likely too large to be found". Cost O(M(d)) for p−1 with geometric-progression points. |
| J. M. Pollard, "Monte Carlo methods for index computation (mod p)", Math. Comp. 32 (1978) 918–924 | bib VERIFIED | Kangaroo/rho for discrete logs; listed as the "Pollard 1978" lead. It is not a factoring stage 2. |
| Smooth-part distribution: Knuth–Trabb Pardo, TCS 3 (1976) 321–348; Hildebrand–Tenenbaum, JTNB 5 (1993) 411–484; R. G. Pinsky, arXiv:1606.02965 (2019 Springer Proc.); Bach–Peralta, Math. Comp. 65 (1996) 1701–1715 | bib VERIFIED; Pinsky Thm 1 VERIFIED | **To check:** a published statement of the e^{-γ}ρ limit for the y-smooth part of a uniform integer (Tenenbaum's book III.5 was not checked); Bach–Peralta content not opened. Our fixed-B statement (§6.1) is proved directly from PNT in APs. |

**Placement of the user's idea:**

| Component | Status |
|---|---|
| Stage 2 as a √X-cost collision search in ⟨b⟩, detected by gcd | Known: Pollard 1974 (sketch), Brent 1986 (birthday), Montgomery–Silverman 1990 (FFT/BSGS) |
| "Also handles composite R" | Implicit in the known algorithms (Kruppa §2.4–2.6). Not emphasised in the literature, and provably irrelevant at the optimum (Prop. D) |
| "Use b^N, because N ≡ q (mod R)" | Correct, but equivalent to order finding with an arbitrary window offset (L3). No prior formulation that uses N this way was found by us; it adds nothing a random offset does not |
| Expected-cost and quantile analysis | Our derivation. The ingredients are standard; the expected-cost Θ(log² N) statement for balanced semiprimes was not found stated anywhere by us |

---

## 8. Candidate card

```text
ID and one-sentence mechanism: PM1-BSGS — stage-1 p−1, then hash-free BSGS (multipoint gcd) for a multiple of ord_p(b) in a window placed at N.
Target problem / supported input family: any N; favourable when p−1 or q−1 has a large B-smooth part.
Representation and how it is built from N: b = a^E mod N; polynomials over Z/N (flint, mul and monic % only).
Exact operation that extracts information: gcd(∏_j F(v_j), N), F = ∏(T − b^{N+i}).
Why that information exposes a factor: p | b^u − b^v ⇔ ord_p(b) | u − v (L4); separation fails only if ord_p(b) = ord_q(b) (L5).
Known lemmas being used: Fermat/Lagrange, CRT, subproduct/remainder trees (Pollard–Strassen), Brent's F′ identity.
New lemma required: none. L1–L5 and Prop. D are elementary and proved here.
Correctness argument or conjecture: proved (§2); predictor ≡ operational on 2800 runs.
Cost in n, including construction and extraction: B + Õ(√X); expected Θ(N^{1/4}/log² N)·polylog (conjectured, §6.1); quantile N^{c(δ)/2} (§6.2).
Randomness and success probability: the base barely matters (100% all-or-nothing); success is a property of p−1 and q−1.
Closest prior method and the actual difference: p−1 with FFT continuation (Montgomery–Silverman 1990; Montgomery–Kruppa 2008). The difference is the window offset N, which is equivalent to an arbitrary offset.
Smallest informative experiment: a_N vs a_order vs a_shift at equal budget (§5.3).
Counterexample or scaling behaviour that would refute this version: (p−1)/2 and (q−1)/2 both prime ⇒ cost Θ̃(√p) (seen at n = 72, seed 982).
Possible repair if that happens: re-roll the group (ECM). Then L_p[1/2, √2], and the stage 2 only touches the o(1).
```

Mechanism verdict: **refuted for this model as a route to a better-than-N^{1/4−o(1)} expected cost.** It is correct and useful as a stage 2, but that stage 2 is known.

---

## 9. Reproduction

```bash
cd /home/knielsen/code/factoring-connection
.venv/bin/python -m pytest -q tests/test_pm1bsgs.py                     # 17 passed
.venv/bin/python -m src.pm1bsgs.experiments smooth   --count 4000       # raw/pm1bsgs_smooth.jsonl.gz
.venv/bin/python -m src.pm1bsgs.experiments phase1   --count 1000       # raw/pm1bsgs_phase1.jsonl.gz
.venv/bin/python -m src.pm1bsgs.experiments rho      --count 500 --start 500   # raw/pm1bsgs_rho.jsonl
.venv/bin/python -m src.pm1bsgs.experiments validate --count 40         # raw/pm1bsgs_validate.jsonl (0 disagreements)
.venv/bin/python -m src.pm1bsgs.experiments walltime                    # raw/pm1bsgs_walltime.jsonl
.venv/bin/python -m src.pm1bsgs.analysis smooth|asym|validate|budget|basedep|fail
```

- Configs are in `experiments/configs/pm1bsgs_*.json`.
- The prime sieve (≤ 2^30) cache is in `experiments/analysis/pm1bsgs_sieve_cache.json`.
- Floating point is used only for the ρ/ω model evaluation and the fits; everything else is exact integers or Fractions.

## 10. Open items / what remains to check

1. **Montgomery–Silverman 1990 full text.** Check the exponent set covered and any remark on composite cofactors.
2. **Pollard 1974.** Check his exact stage-2 and FFT remarks.
3. **Strassen 1976.** Get the primary bibliographic record.
4. **Smooth-part law.** Find a published statement of the e^{-γ}ρ law for the y-smooth part of a uniform integer (Tenenbaum III.5?).
5. **Uniform-in-B estimate.** A uniform-in-B (B ≤ p^{1/2}) version of E[S_B(p−1)^{-1/2}] ≍ 1/log B, i.e. a sieve upper bound plus a Hardy–Littlewood-free lower bound. This would turn the §6.1 Θ(log² N) statement from conjectured into proved (upper bound) and conditional (lower bound).
6. **Generic lower bound.** A precise generic-model lower bound for order finding of a prime-order element, to cite in §6.1.
