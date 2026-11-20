# Pre-Registration: Consultant G Sub-L[1/3] Falsification & Theorem Checks

Timestamp: 2026-11-20 (written strictly BEFORE coding or running any experiment).
Directory: `/Volumes/git/factoring-connection/codex/scratch/subthird_g`

## Experiment E1: Adaptive Constructor Steering vs GRH Flatness & Smooth-Kernel Lattice (`s * a^2`)
- **Target:** Test Candidate 1A (Task 1):
  1. Adaptive residue-guided prime-product beam search (`beam_width = 64`, primes `p <= B`) steering modular distance `min(u mod N, N - (u mod N))` or Gauss-split log-size `ln(|u'|*|v'|)` toward `B`-smooth residues vs oblivious random `B`-smooth products with identical smoothness-test budget (`M = 1024` evaluated leaves per seed).
  2. Smooth-squarefree-kernel constructor (`smoothbuild` W12: `s * a^2 = v (mod N)` via 2-D Euclidean reduction of `Lambda_{a^2}`) vs standard random-target `t` Gauss splits (`M = 1024` splits per seed).
- **Sizes:** Balanced semiprimes `N = p*q` of `n in {32, 40, 48}` bits, `3` deterministic seeds per size; `B = round(2^(n/4))`.
- **Pre-registered Hypotheses & Pass/Fail Criteria:**
  - **H1a (Kill prediction — beam search):** Pooled smooth-residue yield ratio `beam / random` lies in `[0.65, 1.45]` at `48` bits (no super-Dickman steering). *Sub-1/3 lead criterion:* ratio `> 2.5` at `48` bits AND strictly increasing across `32 -> 40 -> 48` bits.
  - **H1b (Kill prediction — `s * a^2` lattice):** `100%` of nontrivial Euclid splits of `a^2 mod N` satisfy `|s| * |v| >= N / 2` (`ln|s| + ln|v| >= ln(N/2)`), and smooth-pair yield ratio `sq_kernel / random_gauss` lies in `[0.65, 1.45]`. *Sub-1/3 lead criterion:* any nontrivial split with `|s|*|v| < N^(0.8)` or yield ratio `> 2.5`.

## Experiment E2: Exact Verification of Farey–Trapezoid Lemma FT (Closing Rank-2 H5)
- **Target:** Verify our new **Lemma FT**: for any squarefree binary form `F(x, y) in Z[x, y]` of degree `D >= 3` (`disc(F) != 0`), `E >= 2`, `X > 0`:
  `#{primitive (a, b) in [-E, E]^2 : 0 < |F(a, b)| <= X} <= 384 * D * 9^((D-1)/D) * A_F * X^(2/D) + 4 * D * (floor(log2(E)) + 2)`.
  (And test tighter practical constant vs Davenport's `2*(D+1)*(2*E + 1)` in the cusp regime where `E >> B`.)
- **Instances:** `12` binary forms `F = F_1 * F_2` (`D in {4, 5}`, GNFS base-`m`, SNFS, and `SL_2(Z)`-skewed `x -> x + 17y` pairs for `N` of `28, 36` bits), evaluated exhaustively on primitive `(a, b)` in `[-E, E]^2` for `E = 256` across dyadic `X` thresholds.
- **Pre-registered Pass/Fail Criteria:**
  - **H2a (Theorem check):** `0` violations of Lemma FT across all `(instance, X)` pairs (`A_F` bounded via Bean's `16 * |disc(F)|^(-1/(D*(D-1)))` AND exact numerical `A_F`).
  - **H2b (Cusp-regime sharpness vs Davenport):** At the smallest threshold `X_min = min_{primitive} |F(a,b)|`, the actual primitive count is `<= 4 * D * (floor(log2(E)) + 2)` in `12/12` instances (whereas Davenport's error term `2*(D+1)*(2*E+1)` exceeds `4 * D * (floor(log2(E)) + 2)` by `> 10x` at `E = 256`).

## Experiment E3: Non-Subfield Block-Degenerate Frame-Less Modules (`grdeg` Door (ii) Closure Test) & Door (i) Relative-Rank-3 Toy Check
- **Target:**
  1. Construct codex-type `(1, s, t, t^2)` vs `(1, s', t, t^2)` frame-less non-`O_k`-module families (`K_0 = Q(t)` cubic `t^3 = 3`, `K = K_0(s)` with `s^2 = c_1(t)`, `K' = K_0(s')` with `s'^2 = c_2(t)` sharing a root mod `N` via LLL-minimized relative polynomials over `O_{K_0}`, `r = 4, D = 12`) for split semiprimes `N = p*q` across `3` size tiers (`~20, ~30, ~40` bits).
  2. Also construct genuine **Door (i)** relative-rank `m = 3` modules over fixed `k = Q(i)` (`e = 2`, `r = 6`, `D = 12`, `K = k(s), K' = k(s')` of relative degree `3` with non-power / LLL-reduced linked generators) across the same `N` tiers to test whether `min_{genuine a in {0,+-1}^6} log|F(a)| / log N` stays bounded near the `c = 0` (`e/d = 1/6`) prediction or grows like `(r-1)/(D-1) = 5/11` (relative R1).
- **Pre-registered Hypotheses & Pass/Fail Criteria:**
  - **H3a (Door (ii) closure):** All constructed `(1, s, t, t^2)` instances have `0` frames (`Delta_F = 0`), `N | Norm_{K_0/Q}(Res_{K_0}(f, g))` holds `100%` of the time, `a_2 = 0` gives coincident norms `W(a) == W'(a)`, and for `a_2 != 0` (`a in {0, +-1}^4`), `min log|F(a)| / log N >= 0.45` (matching our `d_0 / d = 1/2` RR2 floor!).
  - **H3b (Door (i) relative-R1 falsification test):** On the `k = Q(i)`, `m = 3` (`r = 6, D = 12`) LLL-linked modules, measure `slope = min_{genuine a in {0,+-1}^6} log|F(a)| / log N` at `~40`-bit `N`. *Pass criterion for `c < 1` (sub-1/3 signal):* `slope < 0.30` (below `(m-1)/(2n - 1) = 2/5 = 0.40` and approaching `e/d = 1/6 = 0.167`). *Fail (relative-R1 holds empirically):* `slope >= 0.35`.
