# cadobench: CADO-NFS on an Apple M4 as a measuring instrument for the L[1/3] cost model

Agent `cadobench`, phase 15, 2026-10-20. Pre-registration: `notes/cadobench_prereg.md` (commit 8657101,
before every validation run). Code: `src/cadobench/`. Tests: `tests/test_cadobench.py`. Raw data:
`experiments/raw/cadobench_*`. GPT-6.1-sol: `codex/cadobench_1.{md,out}`.

## Verdict

**The heuristic smoothness model is accurate at real sizes. The sieve cost model is not.**

- **Smoothness model (experimentally supported, c60–c110).** Dickman ρ, evaluated at the *actual* norm of
  each point and shifted by Murphy's α, predicts:
  - **Single-side smoothness: within ±11%** (measured/predicted 0.89–1.10) at all six sizes.
  - **Joint smoothness: within the counting noise.** Exhaustive windows give joint/(product of marginals)
    = 1.03, 0.98, 0.98 (c60, c70, c80). Random samples, pooled over c60–c110: 111 doubly smooth points
    observed against 120 predicted.
  - So heuristic H's "joint events multiply" holds, and its ρ-at-actual-size estimate is good to ~10% once
    α is included. Without α the algebraic side is off by a factor 2.0–2.7.
- **Sieve-cost model (refuted for this model at real sizes).** las CPU does **not** scale with box area:
  - **CPU ∝ area^0.52–0.71** when I goes ±1 at fixed special-q.
  - **26–47% of las CPU is a fixed per-special-q cost.** Bucket fill is mostly a per factor-base-prime
    lattice set-up, independent of the area.
  - **Relations ∝ area^0.52–0.58**, because larger regions have larger norms.
  - The "E²" box term is thus only ~50–70% of las time. Cofactorization (P−1/P+1/ECM chains plus trial
    division) is 19–43%.
- **Assumption T (heurnfs §3.1) is about half right.**
  - After singleton removal, 51–59% of large ideals occur in exactly two relations.
  - The rest occur 3–36 times; the mean is 2.6–3.3.
  - 13–25% of all las relations are duplicates found by two special-q.
- **Linear algebra stays small.** It is 9–19% of wall time, the matrix has 5–12% as many rows as there are
  unique relations, and the LA share grows with size. The empirical L[1/3] constant fitted over c70–c110 is
  **c ≈ 1.13**: o(1) terms dominate at these sizes, so the asymptotic constant 1.923 says nothing quantitative
  here.
- **For smoothenum:**
  - las spends **0.18–0.69 ms of single-core CPU per relation** (c70–c110), visiting 16k–47k box points
    per relation at 10–17 ns per point.
  - las misses almost nothing inside its own definition: **recall 99.5–100%**.
  - That definition (S_CADO, which caps the cofactor at 2^mfb) keeps only **57%, 64% and 94%** of the
    doubly-smooth set S_B at c60, c70 and c80.
  - An output-sensitive enumerator must beat ~0.3–0.7 ms per relation per core. Its yield ceiling over
    las at the same B and box is S_B/S_CADO.

Status labels follow AGENT_RULES §5. Everything below is **experimentally supported** at these sizes on
this machine unless marked otherwise. None of it is an asymptotic statement.

## 0. Scope change: task 4 withdrawn, unrun

On the orchestrator's instruction, all constant-factor work was withdrawn before anything was
pre-registered or run:
- cofactorization chains;
- lpb/mfb/lambda tuning;
- NEON/SME/Metal offload;
- bwc tuning.

Nothing in this note claims a speed-up. The work was refocused on (1) baselines, (2) exporting testbeds,
(3) validating the heurnfs cost model, and (4) later, comparisons against smoothenum.

## 1. Build notes (CADO-NFS master 692ecb7e, 2026-10-02, in ~/src/cado-nfs; not committed)

- **Dependencies (Homebrew):** gmp 6, hwloc, libomp, mpfr, wget. cmake 4.3.1 and Apple clang 21 were already
  present. Python deps (flask, requests) come from CADO's `scripts/setup-venv.sh` (venv `cado-nfs.venv`, Python 3.14.7).
- **No source patches were needed.** `local.sh`:
  ```
  eval `/opt/homebrew/bin/brew shellenv`
  export LDFLAGS="-L/opt/homebrew/opt/libomp/lib"; export CPPFLAGS="-I/opt/homebrew/opt/libomp/include"
  CFLAGS="-O2 -DNDEBUG -I/opt/homebrew/opt/libomp/include"; CXXFLAGS="-O2 -DNDEBUG -I/opt/homebrew/opt/libomp/include"
  GMP=/opt/homebrew/opt/gmp; HWLOC=/opt/homebrew/opt/hwloc
  export PATH="$HOME/src/cado-nfs/cado-nfs.venv/bin:$PATH"   # must come AFTER brew shellenv
  ```
- **Pitfalls hit:**
  1. `brew shellenv` puts Homebrew's python3 ahead of the venv, so cmake fails with "Importing the flask
     Python module failed". The fix is the PATH line after it.
  2. Homebrew's fmt is too old for C++20, so CADO falls back to its embedded copy automatically.
  3. hwloc reports "Topology is not symmetric" (4 P-cores + 6 E-cores), so CADO applies no CPU binding.
     This is harmless.
- **Validation.**
  - All six factorizations below returned exactly the generator's factor pair.
  - The harness checks this in `numbers.factors_ok`; harness only, the factors are never passed to CADO.
  - CADO's own `make check` suite was not run. The task allowed "test suite *or* smoke factorization".
- **Numbers.** RSA-style `balanced_semiprime(bits, seed=0)` from `src/common/instances.py`, with bits =
  198/231/264/297/331/364 for c60–c110 (`src/cadobench/numbers.py`).
- **Memory discipline.**
  - Every heavy job ran under `src/cadobench/memwatch.py`, with at most one heavy job at a time. Its
    kill rules:
    - tree RSS > 8 GB (las: 6 GB);
    - available memory < 2 GB;
    - pressure rising to warn.
  - CADO's clients run as `--daemon` and leave the process tree. memwatch therefore also tracks every
    process whose command line matches the CADO paths.
  - Peak RSS ≤ 1.86 GB at every size, so c110 was allowed (gate: c100 peak < 6 GB; measured 1.40 GB).
  - **Watchdog deviation.** From ~16:30 the system sat at pressure level 2 (warn) with ~5 GB available,
    even with none of my jobs running (other agents, Docker VM). The watchdog then killed a V6 analysis and
    two las runs.
    - I made the pressure rule baseline-aware: kill on critical, on a rise above the level at job start,
      on < 2 GB available, or on the RSS cap.
    - The killed runs were redone under the new rule. No job came close to the RSS or available-memory
      limits.

## 2. Phase timings (wall seconds, `-t 10`; from CADO's timestamped log, `run_cado.parse_log`)

| digits | total wall (s) | polyselect | sieve | filter | linalg | sqrt | other | CPU (s) | peak RSS (MB) |
|---|---|---|---|---|---|---|---|---|---|
| c60 | 16.2 | 3.3 (20%) | 4.8 (30%) | 2.0 (13%) | 3.8 (23%) | 1.1 (7%) | 1.0 | 27 | (n/a)¹ |
| c70 | 19.6 | 3.3 (17%) | 9.9 (50%) | 2.6 (13%) | 2.4 (12%) | 0.9 (5%) | 0.4 | 70 | 101 |
| c80 | 69.3 | 6.7 (10%) | 48.0 (69%) | 6.5 (9%) | 6.5 (9%) | 1.3 (2%) | 0.4 | 344 | 139 |
| c90 | 220.1 | 10.1 (5%) | 167.3 (76%) | 17.8 (8%) | 21.3 (10%) | 2.3 (1%) | 1.1 | 1236 | 1089 |
| c100 | 477.5 | 25.6 (5%) | 267.4 (56%) | 114.1 (24%) | 60.5 (13%) | 2.8 (1%) | 7.1 | 2905 | 1398 |
| c110 | 1453.0 | 101.9 (7%) | 960.4 (66%) | 96.4 (7%) | 271.1 (19%) | 5.6 (0%) | 17.1 | 8266 | 1858 |

Notes:
- "linalg" = bwc + characters. "filter" = dup1 + dup2 + purge + merge, summed over every filtering
  round.
- At **c100**, CADO ran singleton removal **10 times**: its default relation target was too low, so each
  purge requested more sieving. Hence the 24% filtering share. CADO's own dup2/purge real-time totals
  (63.3 s / 42.4 s) match our parser.
- Repeat runs of c70–c90 (`*_run2`, `*_run3`) are in §2.1. The c60 row is the smoke run.
- ¹ The smoke run predates the memwatch fix that tracks CADO's detached clients, so its RSS (107 MB) misses las.
- c110 overlapped for ~10 min with niced single-core Python (V1–V3), so its wall time may be inflated by
  a few per cent (disclosed).

### 2.1 Run-to-run spread

See the addendum at the end of this note (filled in after the repeats finished).

## 3. Sieve profile (las, special-q on side 1)

Parameters are CADO's defaults (`parameters/factor/params.cNN`, as snapshotted). Yields are from the full
runs. The funnel is per special-q, from a single-threaded las run on the first 50 special-q at qmin.

| size | I | lim0/lim1 | lpb0/lpb1 | mfb0/mfb1 | special-q sieved | raw rels/q | points/q | after sieve | enter cofact. | relations/q | cofact. survival |
|---|---|---|---|---|---|---|---|---|---|---|---|
| c70 | 11 | 343k/244k | 20/21 | 19/42 | 1488 | 120 | 2.1M | 8739 (4.2e-3) | 172 | 128 | 0.74 |
| c80 | 11 | 293k/340k | 21/21 | 41/42 | 6045 | 52 | 2.1M | 17102 (8.2e-3) | 67 | 46 | 0.68 |
| c90 | 11 | 404k/811k | 23/23 | 46/46 | 19558 | 55 | 2.1M | 18174 (8.7e-3) | 107 | 51 | 0.48 |
| c100 | 11 | 650k/800k | 25/26 | 48/51 | 54380 | 79 | 2.0M | 3551 (1.8e-3) | 147 | 90 | 0.61 |
| c110 | 12 | 1.3M/2.0M | 26/27 | 50/52 | 62063 | 129 | 7.7M | 7571 (9.8e-4) | 243 | 166 | 0.68 |

Reading the table:
- "After sieve" are the points whose sieved log-norms pass the threshold on both sides.
- Trial division then rejects 96–99.6% of them by the size of the leftover cofactor.
- **48–74% of the candidates entering cofactorization become relations.**
- Relations per special-q at qmin are 52–166, falling with q over the run.

**Where las's CPU goes.** Single thread, CADO's I, first 50 special-q. These are `-T` self-times grouped by
`validate.tree_shares`:

| size | CPU/q (ms) | bucket fill | per-point (norms, small sieve, apply buckets, resieve, survivor scan) | trial division | cofactorization (P−1/P+1/ECM) | other |
|---|---|---|---|---|---|---|
| c70 | 23 | 33% | 28% | 11% | 22% | 6% |
| c80 | 29 | 29% | 24% | 12% | 31% | 4% |
| c90 | 35 | 34% | 21% | 11% | 31% | 4% |
| c100 | 27 | 49% | 24% | 7% | 14% | 6% |
| c110 (I=12) | 78 | 45% | 34% | 9% | 10% | 3% |

Caveats:
- **las files its fuzzy `cofactoring_sync` timer under category −1 ("uncategorized").** So las's own
  category totals show only 8–15% cofactoring. The shares above count that timer as cofactorization, which
  is what it is (`las-process-bucket-region.cpp:599`).
- **Multi-threaded CPU seconds are not a work measure on this chip.**
  - The same 50 special-q at c90 cost 1.71 CPU-s at `-t 1` and 2.74 CPU-s at `-t 10` (0.48 s elapsed).
  - E-cores inflate CPU seconds by ~60%.
  - Every work breakdown here therefore uses `-t 1`.
- No sampling profiler (Instruments/xctrace) was used. The scope change removed the offload work it was
  meant to inform.

## 4. Cost-model validation against the pre-registration

| # | Prediction (prereg) | Measured | Verdict |
|---|---|---|---|
| V1a | single-side measured/ρ_α in [0.8, 3.0], both sides, c60–c100 | side 0: 0.887, 0.896, 0.894, 0.918, 0.910; side 1: 1.100, 1.089, 1.068, 0.991, 1.048 (extra: c110 0.935/1.011) | **confirmed** |
| V1b | uncorrected ρ worse than ρ_α on the algebraic side at ≥ 4/5 sizes | 5/5 (uncorrected ratios 2.75, 2.31, 2.21, 2.01, 2.29) | **confirmed** |
| V2 | joint/(marginal product) in [0.67, 1.5], c60–c80 (≥ 30 hits) | S_B: 1.028 (246 hits), 0.978 (616), 0.979 (143); S_CADO: 0.957, 0.998, 1.044 | **confirmed** |
| V3 | las recall of S_CADO ≥ 85% | 139/139, 394/396, 135/135 (100%, 99.5%, 100%); las output ⊂ S_B always | **confirmed** |
| V4a | las CPU per q vs area exponent in [0.85, 1.25] (c90; c100) | `-t 10`: 0.61 (c90), 0.54 (c100); `-t 1`: 0.59, 0.52; other sizes 0.61–0.71 | **refuted** |
| V4b | per-point share (norms + sieving + bucket fill) ≥ 60% and cofactorization ≤ 30% (c90, c100) | c90: 55% and 42% (trial division + cofact.); c100: 73% and 21%. Literal las categories: c90 51%/13%, c100 69%/13% | **refuted at c90** (both parts); **confirmed at c100** |
| V4c | relations per q vs area exponent < 0.9 | 0.58 (c90), 0.52 (c100); 0.52–0.58 at every size | **confirmed** |
| V5a | LA share < 25% at c80–c100 | 9.4%, 9.7%, 12.7% (c110: 18.7%) | **confirmed** |
| V5b | matrix N / unique relations in [0.05, 0.5] | 0.117, 0.078, 0.059 (c70 0.083, c110 0.052) | **confirmed** |
| V5c | LA grows faster than sieve from c80 to c100 | LA ×9.3 vs sieve ×5.6 (ratio 1.67); c100→c110 ×4.5 vs ×3.6 | **confirmed** |
| V6 | after purge: fraction of large ideals with exactly 2 occurrences < 50% **and** mean ≥ 3 (c80–c100) | exactly 2: 51.2%, 50.5%, 55.4% (c110 55.5%); mean 2.85, 3.01, 3.26 (c110 2.97) | **refuted** (first part fails at every size) |
| V7 | fitted L[1/3] constant c in [1.0, 1.9], c70–c100 (c110 included) | total wall 1.13; sieve 1.17; LA 1.24 (local slopes 0.85–1.33) | **confirmed** (descriptive only) |

Details.

**V1 (`validate.v1`; 2·10^5 seeded points per size over 20 special-q spread across CADO's q-range).**
- Murphy α: rational side 0.570 (= Σ log p/(p²−1)); algebraic side −4.81, −4.39, −4.08, −4.43, −5.46, −7.09.
- Our α matches CADO's printed value (−4.81 for c60).
- **The rational side is consistently 9–11% *below* ρ_α**, and the algebraic side 0–10% above. These
  residuals are expected at the second order: the Ψ(x, y)/x vs ρ(u) correction, and the spread of u over
  the box (u0 ≈ 3.0–3.5, u1 ≈ 4.3–4.8).
- The pooled random-sample joint S_B counts are 31, 35, 9, 12, 18, 6, against 31.9, 34.6, 8.0, 10.1, 23.7,
  11.8 predicted by mean ρ_α(u0)·ρ_α(u1), i.e. 111 vs 120. This was not pre-registered; it adds
  confirmation at c90–c110, where V2 windows would be too sparse.

**V2/V3 (`validate.v23`).** Exhaustive windows for 5 special-q per size:
- c60: the full region 1024 × 512.
- c70 and c80: i ∈ [−512, 512), j ∈ [0, 1024), i.e. 2^20 lattice points (≈ 637k coprime) per q.

Exact smoothness is by Bernstein's batch test with remainder trees over all primes ≤ 2^lpb. las was run on
exactly that (q, ρ) with the production parameters. **S_CADO / S_B = 0.57, 0.64, 0.94.** At c60 and c70,
CADO's mfb0 < lpb0 limits the rational side to one large prime and discards ~40% of S_B by design.

**V4 (`validate.v4`).** 200 special-q from qmin at `-t 10` (3 repeats; CPU sd ≤ 5.4%), plus a disclosed
deviation: 50 special-q at `-t 1`, 3 repeats, sd ≤ 4.8%.
- **Post-hoc fit (not pre-registered), CPU = c0 + c1·area.** The fixed part c0 is **26%, 35%, 40%, 47%,
  35%** of CPU at CADO's I (c70–c110), i.e. 6–27 ms per special-q.
  - Its origin is the per-prime p-lattice set-up inside `fill_in_buckets`: π(lim) ≈ 50k–250k primes per
    special-q.
- **The pre-registered regression CPU = c_area·area + c_out·survivors is misspecified.** It has no
  intercept, and survivors are collinear with area. It is kept in the JSON (`fit_cpu_area_out`) but not
  used.

**V6 (`relstats.incidence_stats`).** Our own singleton removal is iterated to a fixed point over large
ideals (p > lim) of all unique relations. CADO additionally removes cliques, so our post-purge sets
(100k / 161k / 479k / 910k / 1.81M relations) are upper bounds on CADO's.

| size | raw → unique (dup %) | large ideals | singletons before purge | relations surviving purge | after purge: exactly 2 / ≥ 3 / mean / max |
|---|---|---|---|---|---|
| c70 | 179k → 156k (13%) | 104k | 44% | 64% | 59% / 41% / 2.63 / 9 |
| c80 | 314k → 258k (18%) | 209k | 32% | 63% | 51% / 49% / 2.85 / 12 |
| c90 | 1.08M → 891k (17%) | 820k | 34% | 54% | 51% / 50% / 3.01 / 18 |
| c100 | 4.28M → 3.22M (25%) | 3.55M | 47% | 28% | 55% / 45% / 3.26 / 36 |
| c110 | 8.03M → 6.05M (25%) | 6.53M | 47% | 30% | 56% / 44% / 2.97 / 20 |

- Relations carry 0–5 large primes. The mode is 3 at c90–c110.
- **At c100 and c110, 70% of unique relations die in singleton removal.**

**V7 (`analyze.v7`).**
- The fit is log T = log A + c·(ln N)^{1/3}(ln ln N)^{2/3} over c70–c110, on means of runs.
- Five points with an implicit o(1) make this purely descriptive. The local slopes jump around
  (c90→c100 = 0.85) because CADO's per-size parameter files and its filtering rounds change discretely.

## 5. What this means for the asymptotic programme (and for smoothenum)

1. **heurnfs's smoothness inputs can be trusted at practical sizes, but only with α.**
   - Heuristic H (ρ at actual norm size, independent sides) predicts yields to ~10%.
   - Any cost estimate that uses box-maximum norms without α is off by the observed factor ~2 on the
     algebraic side, and more from the box maximum.
2. **The "E²" term is not the whole sieve cost at real sizes.**
   - A fixed cost per special-q (∝ π(lim)) is 26–47% of las CPU.
   - Yield per special-q grows only as area^0.55.
   - So the region per special-q that is cheapest in practice is set by this fixed cost, not by E² alone.
   - An output-sensitive enumerator that avoids the box must *also* avoid a per-lattice factor-base pass,
     or pay for it per special-q exactly as las does.
3. **Cheaper smoothness tests are capped.** Cofactorization plus trial division is 19–43% of las CPU, and
   the cap falls with size: c100–c110 spend 21% and 19%.
   - This matches heurnfs's model-PROVED statement that the test cost has slack. Even a free test would
     save at most that constant fraction.
4. **The bar an enumerator must clear (single M4 core, CADO defaults, start of the q-range):**

   | | c70 | c80 | c90 | c100 | c110 |
   |---|---|---|---|---|---|
   | CPU per relation (ms) | 0.18 | 0.63 | 0.69 | 0.30 | 0.47 |
   | box points per relation | 16k | 46k | 41k | 22k | 47k |
   | ns per box point | 11 | 14 | 17 | 13 | 10 |

   - Its yield ceiling over las on the same box and bound is |S_B|/|S_CADO|: ×1.8 (c60), ×1.6 (c70),
     ×1.06 (c80).
   - Any further gain must come from *not visiting* the box.
5. **Assumption T.**
   - Half the large ideals in a real collection behave as T assumes: exactly two incidences.
   - The other half are shared by 3–36 relations, and that tail lengthens with size.
   - On top of this, 13–25% of relations are rediscovered by a second special-q.
   - heurnfs says sharing in special-q collections is "not covered by T". The data show it is the typical
     case for half of the large primes, so any descent-style lower bound that leans on T should state that
     it excludes the observed regime.
6. **LA is not yet in balance with sieving** (sieve : LA wall ≈ 3.5 : 1 at c110). Its share grows from
   9% (c80) to 19% (c110). Under the model, sieve = LA only asymptotically.

## 6. GPT-6.1-sol round (`codex/cadobench_1.out`, high effort, web search)

Summary of its findings; citations are in the .out file, and we did not re-verify them:
- **Already implemented in CADO:**
  - Kleinjung polyselect with size/root optimisation;
  - bucket sieving with multilevel buckets;
  - mixed-representation ECM chains (Bouvier–Imbert, PKC 2020);
  - batch cofactorization (`batch`, `batchlpb`, `batchmfb`);
  - parallel merge (Bouillaguet–Zimmermann 2021);
  - MPI block Wiedemann.
- The default parameter files were tuned with OPAL/NOMAD.
- **The only recent sieve-side experiment it found** is Bouillaguet–Fleury–Fouque–Kirchner, ePrint
  2023/801: batch extraction instead of small-prime sieving, ~1.10× local throughput, ~5% overall.
- It also cites a 2026 RSA-260 GPU report: ~17% from workunit tuning. That is UNVERIFIED by us.
- Its ranked levers for this laptop are all constant-factor tuning (5–20% estimates):
  - thread layout;
  - sieve geometry;
  - lpb/mfb/λ;
  - polyselect budget;
  - one NEON kernel.
- It could not substantiate five reliable ≥ 5% wins.
- This agrees with the scope change: nothing known offers an asymptotic lever, and the constant-factor
  levers are mostly already implemented.

## 7. Testbed data format (for smoothenum)

Directory `experiments/raw/cadobench_testbed_c{60,70,80,90,100,110}/` (from `python -m cadobench.export`):

- **`meta.json`**
  - N, generator call, CADO commit (692ecb7e…);
  - `f0_rational_low_to_high` and `f1_algebraic_low_to_high` (coefficient strings, low degree first);
  - skew, the poly file, CADO parameters (lim, lpb, mfb, lambda, I, ncurves), qmin, α;
  - `conventions`:
    - F(a, b) = Σ c_k a^k b^{d−k};
    - special-q on side 1;
    - **a = i·a0 + j·a1, b = i·b0 + j·b1, then (a, b) → (−a, −b) if b < 0**;
    - i ∈ [−I/2, I/2), j ∈ [0, J);
    - n0 = |G(a, b)|, n1 = |F(a, b)|/q;
    - the S_B and S_CADO definitions.
- **`qlattices.jsonl`** (c70–c110): one line per special-q for the 200 consecutive special-q from qmin
  (the V4 set):
  - `q, rho, a0, b0, a1, b1, I, J` exactly as las printed them;
  - `las_relations`: how many relations las output for that (q, ρ).
- **`las_relations.jsonl.gz`** (c70–c110): every relation las output on that set with production
  parameters at CADO's I (`-t 10`). Fields:
  - `q, rho, a, b, i, j`;
  - `line`, the CADO relation `a,b:hex primes of G:hex primes of F` with multiplicity, q included.
  - A relation divisible by two special-q of the set appears once per (q, ρ).
- **`las_summary.json`**: las CPU and elapsed seconds (3 repeats), funnel, `-T` categories and the exact
  command line for that set. Use it as the las yardstick on identical inputs.
- **`truth_windows.json`** + **`truth_hits.jsonl.gz`** (c60, c70, c80; the V2/V3 windows):
  - the window per special-q and exact counts (marginals, |S_B|, |S_CADO|, las recall);
  - every S_B point of the window: `q, rho, i, j, a, b, SB, SCADO, las`.
  - Anything not listed is exactly *not* S_B. The test is exact, Bernstein batch over all primes ≤ 2^lpb.

To compare a prototype, run it on the same (f, g), the same q-lattice (or the same (a, b) box) and the same
2^lpb. Report:
- the set of pairs found, scored against `truth_hits` (recall of S_B and of S_CADO);
- single-core CPU per pair;
- the same numbers for las from `las_summary.json`, or `python -m cadobench.lasbench --testbed c90_run1
  --q0 <q> --q1 <q+1> -t 1`.

## 8. Harness and reproduction

- `src/cadobench/run_cado.py`: full factorization under memwatch; per-phase wall times; key parameters.
  Command: `python -m cadobench.run_cado --digits 90 --out ...`.
- `src/cadobench/lasbench.py`: las on a fixed q-range with `--set key=value` overrides and repeats. It
  parses relations, funnel, `-T` categories and tree, CPU, elapsed, peak memory, and reports relations per
  elapsed second and per CPU second.
- `src/cadobench/smoothcount.py`: exact batch smoothness (product and remainder trees), Murphy α (exact
  root counts mod p^k, p ≤ 2000), Dickman ρ (float, model values only), q-lattice maps.
- `src/cadobench/relstats.py`: accepted relation files, per-q lattices and yields, incidence and
  singleton removal (numpy), filtering and matrix numbers from the log.
- `src/cadobench/validate.py`: V1, V2/V3, V4. `src/cadobench/analyze.py`: phase table, V5, V6, V7.
  `src/cadobench/export.py`: testbeds.
- **Tests:** `.venv/bin/python -m pytest -q tests/test_cadobench.py`. 11 tests:
  - exact smoothness against brute force;
  - ρ values;
  - α against CADO's printed value;
  - norms against a real CADO relation;
  - lattice round trip;
  - log and las parsers;
  - the memwatch kill;
  - deterministic las yield, which needs the local build.

Commands (from the repo root, `PYTHONPATH=src`):
```
python -m cadobench.run_cado --digits D --threads 10 --tag run1 --out experiments/raw/cadobench_cD_run1.json   # D = 70..110
python -m cadobench.validate v1  --digits D --out experiments/raw/cadobench_v1_cD.json                          # D = 60..110
python -m cadobench.validate v23 --digits D --out experiments/raw/cadobench_v23_cD.json                         # D = 60, 70, 80
python -m cadobench.validate v4  --digits D --out experiments/raw/cadobench_v4_cD.json                          # -t 10, 200 q
python -m cadobench.validate v4  --digits D --nq 50 --threads 1 --out experiments/raw/cadobench_v4t1_cD.json    # -t 1, 50 q
python -m cadobench.analyze --v6                                                                                # V5, V6, V7, tables
python -m cadobench.export --digits 60 70 80 90 100 110
```
Each heavy command ran on its own; the c60 baseline is the smoke run (`~/cadobench-work/c60_smoke`).

## 9. Disclosed deviations from the pre-registration

1. **V4 rerun at `-t 1`.** The pre-registered runs were at `-t 10`. On P+E cores, multithreaded CPU
   seconds inflate by ~60%, so the work breakdowns use the `-t 1` reruns (50 special-q). Both data sets are
   kept and agree on the exponents (0.52–0.71).
2. **V4b category definition.** las mis-files the cofactorization timer (§3). Under the literal categories
   the cofactorization share is 13–15% (≤ 30%, "confirmed"); under correct accounting it is 21–42%. The
   verdict above uses the correct accounting: refuted at c90, confirmed at c100.
3. **V4 post-hoc intercept fit.** It was added because the pre-registered no-intercept fit is
   misspecified.
4. **Extra sizes beyond the pre-registration:**
   - V1, V4, V6 at c110;
   - V4 at c70 and c80;
   - V5/V7 including c110;
   - the pooled V1 joint count.
5. **memwatch pressure rule** made baseline-aware mid-session (§1). The killed runs were redone.
6. **c110 timing** overlapped with niced single-core Python work (§2).
