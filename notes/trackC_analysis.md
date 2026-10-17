# Track C analysis: the factorial unit predicate

> **r3 audit applied 2026-10-17 (agent `repair`).** Corrections are marked "Corrected 2026-10-17 after r3 audit". The adjudication table is in `notes/repair_adjudication.md`. Reproductions: `src/repair/audit_trackC.py`.

Owner: factC. Date 2026-10-17. Companion notes: `notes/trackC_tal_audit.md` and `notes/trackC_candidates.md`.

**Predicates.**
- `U(K) := [gcd((K−1)!, N) = 1]`, which for composite N is equivalent to `K ≤ p`, the least prime factor.
- `I(A,B) := [gcd(∏_{A≤j<B} j, N) > 1]`, i.e. `[A, B)` contains a multiple of some prime factor.

Write `L = B − A`. `Z_L := {A mod p : [A, A+L) contains a multiple of p} = {0, −1, …, −(L−1)}`.

---

## 1. Baselines (instrumented)

**Code.**
- `src/trackC/polyback.py`: dense polynomials over Z/NZ, with a flint backend and a pure-Python Karatsuba backend that counts R_N multiplications exactly. Monic division is by Newton inversion, so no scalar inverse is ever taken inside polynomial code.
- `src/trackC/baselines.py`:

| Method | What it does | Cost |
|---|---|---|
| (a) `naive_interval` | running product; one gcd per block of 512; Lemma-12 descent to locate the nonunit | L mults |
| (b) `ps_interval` | Pollard–Strassen: subproduct tree for `F(x) = ∏_{i<b}(x+A+i)`, `b = ⌊√L⌋`, then a monic remainder tree at `x = 0, b, 2b, …` | `O(M(√L) log L)` |
| (c) `bgs_interval` | BGS Thm-5 shifting of evaluations: keep `g_d(x) = ∏_{j<d}(A + v x + j)` at `x = 0..d`; doubling `d → 2d` uses 3 Lagrange shifts (one middle product each); increment `d → d+1` direct | `O(M(√L))` |

For (c), every needed inverse is gcd-checked by batch inversion. A hit is reported as `aux_hit` with a witness integer, and `interval_predicate` then falls back to PS. The aux integers are `1..2v+3` and `d + v·t` with `|t| ≤ 2d+1`, all at most about `3L`.

**Flint audit** (complements `src/common/README.md`):
- `nmod_poly` divmod by a polynomial whose leading coefficient is a nonunit mod composite N **aborts the process** with a C-level "Impossible inverse" error that Python cannot catch.
- `fmpz_mod_poly` raises `ValueError` without reporting the gcd.
- So operational code never calls flint division. Multiplication and `mul_low` agree with the Karatsuba backend (tests).

**Tests.** `tests/test_trackC_baselines.py` covers:
- exhaustive `N < 400` with random intervals;
- 10 composite families (prime powers, Carmichael, 3- and 4-prime products, primes), on both backends;
- the prefix predicate, and binary-search recovery of `p`.

All pass: 49 tests in the repo, of which 12 are trackC baselines and 5 trackC mechanisms.

**Scaling.** Interval `[1, K)` with `K < p` (full work), balanced 64-bit and 128-bit `N` (seed 0). Raw data: `experiments/raw/trackC_baselines.jsonl`; fits: `experiments/analysis/trackC_baselines_fit.json`. Selected 64-bit rows:

| log2 K | naive rn_mul / s | PS pmul_model / s | BGS pmul_model / s | BGS gcds |
|---:|---|---|---|---:|
| 16 | 65 663 / 0.003 | 126 230 / 0.007 | 18 288 / 0.002 | 24 |
| 20 | 1 050 623 / 0.053 | 748 596 / 0.033 | 91 953 / 0.008 | 30 |
| 24 | 16 809 983 / 0.83 | 4 165 722 / 0.15 | 442 086 / 0.033 | 36 |
| 28 | — | 22 134 920 / 0.74 | 2 064 015 / 0.14 | 42 |
| 30 | — | 50 331 810 / 1.60 | 4 423 263 / 0.29 | 45 |

`pmul_model = Σ n·⌈log2 n⌉` over polynomial products with `n = len a + len b`. That is an FFT cost model. Wall times are from flint.

Log-log slopes for `K ∈ [2^16, K_max]`:

| Method | wall | exact rn_mul (Karatsuba backend) | pmul_model |
|---|---:|---:|---:|
| naive | 1.00 | 1.000 | — |
| PS | 0.56 | 0.815 | 0.617 |
| BGS | 0.52 | 0.763 | 0.563 |

These match the expected exponents: `K^{1/2}` times logs for an FFT multiply, and `K^{0.79}` for Karatsuba. The PS/BGS ratio of model cost grows like `log K`: 6.9× at `2^16`, 11.4× at `2^30`. That is the BGS log-factor saving. 128-bit `N` gives identical slopes.

The largest materialized polynomial has `√K + 1` coefficients for PS and `≈ 3√K` for BGS (the middle-product convolution).

**End-to-end** (E7): `least_factor_by_search` with BGS `U(K)` found the correct least factor on all 13 instances:
- balanced 32/48/64-bit;
- unbalanced 64/20;
- `p²`, `p³`;
- 3- and 4-prime products;
- Carmichael 6601 = 7·23·41;
- a 48-bit prime (correct "no factor").

Each needs `≈ log2 √N` predicate calls. The 64-bit balanced case takes 33 calls, `pmul_model ≈ 2.7·10^8`, 29 s. This is a correctness check only: the cost is `Θ̃(N^{1/4})`, as expected.

**Status.** `experimentally supported`: the three baselines agree, and their counters fit exponents 1, 1/2+log and 1/2.

---

## 2. Unit-only mechanisms: definitions and outcomes

Full cards are in `notes/trackC_candidates.md`. Summary of operations and experiments:

**M1. Random-shift probes.** Unit status of `∏_{i<L}(s+i)` at random `s` (BGS per probe).
- Proved: by CRT, `p | ∏(s+i)` iff `s mod p ∈ Z_L`. That event is independent of the target `A`'s residue, so a probe carries no information about `I(A, A+L)` beyond the null probability `L/p + L/q`.
- E1 (12 balanced 32-bit N, `L = 2048`, 40 probes per target): hit rate 32/480 = 0.067 when `I(A)` is true and 37/480 = 0.077 when it is false. The null is 0.062–0.077.
- All 69 hits are explained by the probe's own interval. They are STRUCTURAL for the probe interval and ACCIDENTAL for the target.
- `refuted for this model`.

**M2. Power sums.** `S_e(A,B) = Σ_{A≤j<B} j^e mod N` by Faulhaber, with Bernoulli numbers mod N via gcd-checked `1/(n+1)`. Cost `O(e²)` ring operations plus `O(e)` powers. Results:
- Proved (degree): for fixed `L`, `S_e(A, A+L) mod p` is a polynomial in `A` of degree `≤ e`. It has at most `e` zeros over all `p` residues, while a test needs zeros exactly on `Z_L` (`L` points). So `e ≥ L` is necessary (for `L < p`). ~~Also `e < p−1` is needed for the formula to be a degree-`e` polynomial at all.~~ **Corrected 2026-10-17 after r3 audit:** direct binomial expansion `Σ_i (A+i)^e = Σ_k C(e,k) A^k Σ_i i^{e−k}` is always a polynomial of degree `≤ e`, for every `e` (checked at `p = 7, e = 10`). The denominator obstruction below concerns the implemented Faulhaber/Bernoulli evaluator, not every evaluator. For `e ≡ 0 mod (p−1)`, `S = L − #multiples`, which is never `0 mod p` for `L < p−1`.
- Proved (cost circularity, *for the Faulhaber evaluator*): Faulhaber of degree `e` needs `1/k` for all `k ≤ e+1`. With `e ≥ L`, these inversions are exactly a unit test of `(e+1)!`. A failure exposes a prime `≤ e+1` (trial-division-type), and the work is at least `e ≥ L`.
- E2 is an exact count over all `A mod p` for 8 primes `p ∈ [500, 4000]`, `L ∈ {5, 50, p/4}`, and `e ∈ {1,2,3,5,8,13, p−2, p−1, p}`. The number of zeros never exceeds `e`. The zero rate is `≤ e/p` both inside and outside `Z_L`.
- One explained bias: odd `e` show a higher rate inside `Z_L`. That is because `Σ_{|t|≤h}(c+t)^e` with odd `e` is divisible by the center `c` (tested). The "signal" is literally `gcd(center, N)`: a single-index query, the BRIEF §2.3 trap. I classify these as STRUCTURAL-TRIVIAL.
- `refuted for this model`. (M1 scope, r3 audit: the independence argument is for independent uniform shifts. It does not cover correlated shifts or arbitrary combinations.)

**M2′. Fermat-type sums `Σ j^{N−1}`.** Mod `p`, `N−1 ≡ q−1 (mod p−1)`, so this is the case `e = q−1 mod (p−1)` of M2: a generic exponent. Direct evaluation costs `L` powerings. No signal (E2 rows `e = p−2, p` behave like generic `e`). `refuted for this model`.

**M3. Truncated `F_N`.** `Z(N,K) = [(1+X)^N ≡ 1 mod (N, X^K)]` (`truncated_FN_test`).
- E3 checked agreement with `U(K)` exhaustively for composite `N < 600` plus four Carmichael numbers: 0 mismatches.
- Cost: `⌈log2 N⌉ + popcount` products of length-`K` truncated series. At 64 bits: 93 products, `fn_model = 2.07·10^8` at `K = 2^16`, versus BGS `1.8·10^4`, 11 000× more. The slope is 1 in `K`, against 1/2 for BGS.
- Dutta–Jindal–Pandey–Sinhababu (CCC'21) Thm 30, read by factC in the LIPIcs PDF: polylog constant-free circuits for `trunc(1/(1+x)^{d+1}, m)`, `m ∈ {d−1, d}`, imply `τ(n!) = O(log^{c+1} n)`. Their Thm 32: small truncations of all `(1+k²x)^{i/k}` give nonuniform polynomial-time factoring.
- ~~So a short-circuit *truncation* is at least as hard as a short factorial program.~~ **Corrected 2026-10-17 after r3 audit:** the Dutta et al. theorems concern specified characteristic-zero truncation families. They do not establish a hardness reduction for this modular `F_N` truncation, which would need another reduction (`not yet tested`). What is `proved` here is only the equivalence of the predicate `Z(N, K)` with `U(K)`. `refuted as a speedup` (the evaluator is linear in `K`).

**M4. Jacobi-symbol prefix product.** `C_N(K) = ∏_{j<K}(j/N)` (codex idea 2). It is zero iff `U(K)` fails (tested, odd `N < 400`).
- It discards all information except three values. But the only known evaluator is `K` Jacobi symbols, or BGS followed by one Jacobi symbol.
- `Σ_{j<K}(j/N)²` is the coprime count `K − 1 − ⌊(K−1)/p⌋ − ⌊(K−1)/q⌋ + …`.
- The characters mod N computable without factoring are the Jacobi symbol and its powers. The relevant one, `(·/N)²`, is the principal character itself, so character theory gives nothing beyond the original count.
- `not yet tested` for a fast evaluator; no candidate operation is known.

**M5. Compressed discriminant** (codex card C1). `disc(∏_{j<K}(T−j)) = ±∏_{k<K} k!²` is a unit iff `U(K)`. The recursion `D_{2m} = D_m² Res(P_m(T), P_m(T−m))²` needs the cross-resultant `∏_{j<m} P_m(j+m)`. That is `m` evaluations of a degree-`m` polynomial, `Õ(m)`. The total is `Õ(K)`, worse than BGS's `Õ(√K)` for the same predicate. `refuted as a speedup` *for this evaluator*. `Õ(K)` describes the proposed evaluator; it is not a lower bound against all evaluators (r3 audit). The cross-resultant may differ by a sign, which is irrelevant after squaring. The missing primitive is a sub-linear resultant of shifted Pochhammers, which is the same problem again.

**M6. Short-circuit group-law identities.** `x^K − 1`, Dickson/Lucas, division polynomials. See §3, O4: for `x^K − 1`, containing a long prefix forces `(p−1) | K` (`proved`, lemma below; known prior art: Pollard `p−1`). **Corrected 2026-10-17 after r3 audit:** the "(resp. the group order)" analogues for Dickson/Lucas and division polynomials are **not proved**. An integer prefix has no specified group-compatible embedding into a torus or a curve. Moreover, annihilating a point needs its **order**, not `#E(F_p)`, to divide `K` (e.g. on `y^2 = x^3+2x+3` over `F_101`, `#E = 96`, and the point `(3,6)` has order 12, so `12·P = O` while `96 ∤ 12`). Status of those analogues: `conjectured` / `not yet tested`.

---

## 3. Structural question: can a unit-only invariant have log-size dimension or degree?

**Model conventions.** An invariant is *faithful at p* if, for the prime `p | N` in question, its readout is `0 mod p` exactly when the predicate's `p`-component is nonunit. This is automatic for `N = p^k`. It is also automatic for any invariant built componentwise by CRT-uniform ring operations whose correctness argument is per prime, as all known methods are.

Not covered: "mixed" invariants for `N = pq` that are correct only because the `q`-component happens to vanish where the `p`-component does not.

Everything below holds even if all coefficients, matrices and exponents may depend on `p` (non-uniform).

**O1 (C-finite / finite-algebra powers).** `proved`.

*Statement.* Let `u_K = ℓ(M^K v)` with `M` an `R×R` matrix over `F_p`. Equivalently, `u_K = λ(α^K)` for `α` in a commutative `F_p`-algebra of rank `R`, or `u` satisfies a constant-coefficient linear recurrence of order `R`. Then:
- (a) if `u_K ≠ 0` for `1 ≤ K ≤ p` and `u_K = 0` for `p < K ≤ p+R`, then `R ≥ p+1`;
- (b) if `u_A = 0 ⟺ A mod p ∈ Z_L` for `1 ≤ A < 2p`, with `1 ≤ L ≤ p−1`, then `R ≥ L+1`. This is tight: `u_A = ∏_{i<L}(A+i)` has order `L+1`.

*Proof.* Write the characteristic polynomial of the recurrence as `x^t g(x)` with `g(0) ≠ 0` and `t ≤ R`. For `K ≥ t`, `(u_K)` satisfies the invertible recurrence `g` of order `R − t`. An invertible recurrence that has `R − t` consecutive zeros vanishes identically, forwards and backwards.
- In (a), there are `R ≥ R − t` consecutive zeros at `p+1..p+R`, so `u_K = 0` for all `K ≥ t`. In particular `u_R = 0`. If `R ≤ p` this contradicts `u_R ≠ 0`.
- In (b), the run `A = p−L+1..p` has `L` zeros. If `R ≤ L`, the same argument gives `u_A = 0` for all `A ≥ t`, including `A = p+1`. But `[p+1, p+L]` contains no multiple of `p`. Contradiction. ∎

*Corollaries.*
- Any invariant of the form `K ↦ (1+X)^K mod (N, f(X))` with `deg f = r`, or `K ↦ M^K` with a linear readout, needs dimension `r ≥ p+1` to carry `U(K)`. All such invariants have `O(R^ω log K)` jump-ahead, which is exactly why they cannot carry it.
- **Corrected 2026-10-17 after r3 audit (polynomial readouts).** For a polynomial readout of degree `D` of a dimension-`R` state, O1 bounds the **lifted dimension** `C(R+D, D)`, not `R`. Its hypothesis also needs zeros on `p < K ≤ p + C(R+D, D)`. The earlier claim "needs `R ≥ p+1`" is **wrong**. Counterexample, reproduced for `p ∈ {7, 11, 13, 31, 101, 257}`: `α` primitive in `F_{p^2}`, state `z_K = α^K` (`R = 2`), degree-4 readout `Norm(z − α^{p+1})·Norm(z − α^{p+2})`. It is nonzero for `1 ≤ K ≤ p` and zero at `K = p+1, p+2`. Interval version: state `(A, 1)` (`R = 2`) with readout `∏_{i<L}(A+i)` carries `I(A, A+L)` exactly (`p = 101, L = 7`). `audit_trackC.o1_poly_readout`.
- AKS-style quotient sketches indexed by the exponent are in this class.
- E4: an exhaustive check over all order-`R` sequences for `p = 5, 7` (`R ≤ 3`) and `p = 11` (`R ≤ 2`) finds 0 sequences with either pattern. The tightness example `p = 11, L = 3, R = 4` passes.

**O2 (one power-map layer).** `proved`.

*Statement.* Let `P(x) = Q((x+u)^e)` over `F_p`, with `p > 3`, arbitrary `e` (even huge or `p`-dependent), and arbitrary `u, Q`. If the zero set of `P` is exactly `Z_L` with `2 ≤ L ≤ p−2`, then `Q` has at least `⌈(L−1)/2⌉` distinct roots, so `deg Q ≥ (L−1)/2`.

*Proof.* Let `g = gcd(e, p−1)`. The nonzero fibres of `y ↦ y^e` are cosets of `μ_g`. The set `T = Z_L + u`, which is `L` consecutive residues, must be a union of fibres, so `T \ {0}` is `μ_g`-invariant.
- If `g ≥ 3`, then for `k ∈ {1, 2}` the sum of `y^k` over each coset is 0. So `Σ_T y = Σ_T y² = 0`.
- For `T = {a, …, a+L−1}`, the first equation gives `a = −(L−1)/2`. Substituting, the second becomes `L(L−1)(L+1)/12 ≡ 0`. This is impossible for `p > 3` and `2 ≤ L ≤ p−2`.
- Hence `g ≤ 2`. Then fibres have size `≤ 2`, so `Q` needs at least `⌈(L−1)/2⌉` roots. ∎
- **Strengthened 2026-10-17 after r3 audit (proved):** every fibre has at most two points, including the singleton fibre `{0}`. So the `L` points of `T` lie in at least `⌈L/2⌉` fibres, and **`deg Q ≥ ⌈L/2⌉`**. Exhaustive check, all primes `5 ≤ p < 100`, all `e, u, L`: 1 693 434 exact one-layer configurations, 0 violations, and the bound `⌈L/2⌉` is attained for all 986 `(p, L)` with an exact fold (`audit_trackC.o2_strengthening`).

So power maps buy at most the classical `∏(x² − i²)` factor of 2. E5 (exhaustive, all primes `p < 140`, all `4 ≤ L ≤ p−2`, all `g | p−1`, all `u`; 1962 pairs) finds 0 violations. The minimum ratio is 1.007, so the bound is essentially attained by `g = 2` with a centered interval.

**O2′ (k layers, fibre bound).** `proved`, elementary. If `φ` is any composition of translations and power maps `y ↦ y^{e_i}`, its fibres have size `≤ ∏ g_i` with `g_i = gcd(e_i, p−1)`. So `#roots(Q) ≥ L / ∏ g_i`.

~~A large compression therefore needs exponents sharing a large factor with `p−1`.~~ **Corrected 2026-10-17 after r3 audit:** O2′ does **not** exclude large compression from many layers with bounded gcds. With `g_i ≤ G` the bound is `#roots(Q) ≥ L/G^k`, which becomes trivial at `k = ⌈log_G L⌉ = O(log L)` layers (e.g. `G = 2, L = 1000, k = 10`: bound `< 1`). With exponents independent of `p` such as `lcm(1..B)`, `g_i = gcd(lcm(1..B), p−1) = ∏_{ℓ ≤ B} ℓ^{min(v_ℓ(p−1), ⌊log_ℓ B⌋)}`, the truncated prime-power-smooth part. This is **not** always the `B`-smooth part: for `p = 257, B = 10` they are 8 and 256, and they differ for 536 of the 2 261 primes below 20 000. It is still Pollard `p−1` structure.

E8 (exhaustive two-layer search, `p ≤ 31`) shows the one-layer factor 2 is not a theorem for two layers:
- With `gcd(e_i, p−1) ≤ 2`, two layers reach about `0.6–0.8 × ⌈(L−1)/2⌉` roots, for example `p = 31, L = 29: 11` versus 14.
- Allowing `e ≡ 0` or `(p−1)/2 mod (p−1)` (Fermat indicator, Legendre symbol) gives `O(1)` roots by small-`p` coincidences of character runs. For large `p`, runs of quadratic residues are `≪ p^{1/4+ε}` (Burgess, UNVERIFIED citation), so these cannot persist to `L ≈ √p`.
- ~~Status of a general k-layer bound beyond O2′: `conjectured` (`#roots(Q) ≥ c^{−k} L` for exponents with bounded gcd).~~ **Corrected 2026-10-17 after r3 audit:** in that form it is already `proved` by O2′ with `c = G` (`g_i ≤ G` ⇒ `#roots ≥ G^{−k} L`). The open question is a better constant, e.g. whether two gcd-2 layers can beat `3L/8`. The C-2R cycle (`notes/repair_C2R_results.md`) answers it for the tested family. Only the character exponent `E ≡ (p+1)/2` achieves it, and that exponent encodes `p−1`.

**O3 (bounded formal degree).** `proved`, trivial. If the test value is `P(A)` with `deg P < p`, then `deg P ≥ L`, because `P` needs `L` roots mod `p`. ~~Any escape needs formal degree `≥ p`.~~ **Corrected 2026-10-17 after r3 audit: WRONG.** O3 only forces degree `≥ L`, and degree is not circuit size. A hypothetical short circuit for the Pochhammer polynomial `∏_{i<L}(x+i)` has degree exactly `L < p` and the exact zero set `Z_L` (checked at `p = 101, L = 10`). Large formal degree is therefore *not* a necessary feature of an escape. The cheap sources of degree `≥ p` listed (power maps, algebra powers, group laws) remain the known ones.

**O4 (group-law identities contain long prefixes only via the group exponent).** `proved`.

*Lemma.* If `H < F_p^*` is a proper subgroup, its least positive non-member is `< √p + 1`.

*Proof.* Let `n ∉ H` be minimal and `m = ⌈p/n⌉`. Then `0 < mn − p < n`, so `mn ≡ mn − p ∈ H`, hence `m ∉ H`, hence `m ≥ n`. So `p/n + 1 > n`. ∎

*Consequences.*
- If `{1..L−1} ⊂ μ_K(F_p)` with `L − 1 ≥ √p + 1`, then `μ_K = F_p^*`, i.e. `(p−1) | K`. A short-circuit polynomial `x^K − 1` divisible by `∏_{i<L}(x−i)` mod `p` is exactly the Pollard `p−1` condition.
- E6 checks all primes `p < 20000`: the maximal `L` with `[1, L]` inside a proper subgroup is `≤ 0.69(√p + 1)`.
- The analogous statements for Dickson/Lucas (`p+1`) and elliptic division polynomials are `conjectured` by the same pattern, not proved here (r3 audit: treating them as proved is wrong). For curves, the correct quantity is the **order of the point**, not `#E` (example in M6). There is also no specified group-compatible embedding of an integer prefix into a curve.
- **Novelty calibration (r3 audit, verified 2026-10-17 against the primary source).** The `√p + 1` lemma is much weaker than known results. Martin–Parvardi, *Subproducts of small residue classes*, Canad. Math. Bull. 65(1) (2022) 1–8, Corollary 3: `y(p) ≪_ε p^{1/(4√e)+ε}`, where `y(p)` is the least `y` such that subset products of `{1,…,y}` cover `(Z/pZ)^*`. They also cite Burthe: the integers up to `O_ε(p^{1/(4√e)+ε})` generate the group. Either result shows that no proper subgroup contains `[1, L−1]` once `L ≫ p^{0.152}`. So O4's consequence holds at a much shorter prefix, and O4 is **not a new asymptotic obstruction**. Numerically, the least `y` with `⟨1..y⟩ = F_p^*` is at most 29 for all `p < 20 000`, against `√p ≈ 141` (`audit_trackC.generation_bound`).

**Connection to the literature.**
- Shub–Smale 1995 (statement via Bürgisser 2009, ECCC TR06-113, read by factC): the τ-conjecture implies `P_C ≠ NP_C`. Their proof needs only that no nonzero multiples `m_n n!` have polylog `τ`.
- Their Theorem 2 (per `notes/prior_art_map.md`, not re-read by factC) says a polylog BSS decision procedure for `z ∈ {1, …, t}` makes `t!` ultimately easy. That is precisely the interval-zero-set problem that O1–O4 address in restricted models.
- Lipton 1994 (abstract and preview only; full theorem UNVERIFIED): polylog SLPs for polynomials with many rational roots ⇒ factoring is not hard on average.
- Lipton 2009: polylog SLPs for `t!` ⇒ polynomial-size factoring circuits by gcd binary search, which is exactly our prefix reduction.
- Dutta et al. Thm 30/32 (read): the truncation analog.

Our obstructions are unconditional, but only in restricted models. The general SLP model is equivalent to these open problems. Proving `τ`-type lower bounds there would imply new algebraic-complexity separations (Bürgisser Thm 1.1: hardness of `n!` ⇒ `τ(Per_n)` superpolynomial). So an unconditional obstruction for general SLPs is out of reach.

**What an escape must have.** All of the following, simultaneously. (Item 1 withdrawn 2026-10-17 after r3 audit: O3 gives only degree `≥ L`. Item 4 is proved only for `x^K − 1`.)
1. ~~formal degree `≥ p` (O3)~~ formal degree `≥ L` if the test is a single polynomial of degree `< p` (O3);
2. a state that is not a fixed-dimension power `M^K` (O1);
3. not a one-layer power map, and, if a composition of `k` power maps, not one with `∏ g_i ≪ L` (O2/O2′; `k = O(log L)` layers with bounded gcd are *not* excluded by O2′);
4. a zero set mod `p` that is additively ordered (an interval), not a union of cosets of a subgroup of some algebraic group (O4);
5. uniformity over all hidden `p ≤ √N`, since `p`-dependent exponents or characters are unavailable.

The known operations with polylog jump-ahead (matrix powers, power maps, group laws: Lucas, elliptic) all produce subgroup-structured zero sets (`proved` only for `x^K − 1` and linear readouts; heuristic for the others). Interval-structured zero sets are known only from explicit products (`√`-cost via BGS) or from `p`-specific Gauss/Jacobi sums (Tal).

The single most precise open question: **is there a family of polylog-size division-free SLPs over `Z[x]`, uniform in `L`, whose value at `x` is `0 mod p` iff `x mod p ∈ Z_L`, for every prime `p` in a dyadic range above `L`?** This is the mod-`p` form of Shub–Smale's Theorem 2.

---

## 4. Reproduction

```bash
cd /home/knielsen/code/factoring-connection
.venv/bin/python -m pytest -q tests/test_trackC_baselines.py tests/test_trackC_mechanisms.py
.venv/bin/python -m src.trackC.exp_baselines --config experiments/configs/trackC_baselines.json   # ~10 min, 1 core
for e in E1 E2 E3 E4 E5 E6 E7 E8; do .venv/bin/python -m src.trackC.harness_mechanisms --exp $e; done
E8_VARIANT=no_zero_exp .venv/bin/python -m src.trackC.harness_mechanisms --exp E8
E8_VARIANT=gcd_le2     .venv/bin/python -m src.trackC.harness_mechanisms --exp E8
```

**Outputs.**
- `experiments/raw/trackC_baselines.jsonl` and `experiments/analysis/trackC_baselines_fit.json`.
- `experiments/raw/trackC_mechanisms.jsonl`, with all E1–E8 records. The first E8 run, without a `variant` field, is variant `all`.

**Instances.** `src/common/instances.py` generators. Seeds are given in the scripts: `balanced_semiprime(bits, seed)` and so on.

**CPU.** Everything is single-core. The largest jobs are about 10 min for the baselines and about 5 min for E8.
