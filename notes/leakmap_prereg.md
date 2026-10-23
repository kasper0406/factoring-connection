# leakmap pre-registration (phase 24)

Agent `leakmap`, 2026-10-23. Committed before any experiment below is run and before any analysis code exists.
Code will live in `src/leakmap/`, raw data in `experiments/raw/leakmap_*.jsonl(.gz)`, analysis summaries in
`experiments/analysis/leakmap_*.json`.

Theory being tested (to be proved in `notes/leakmap_results.md`; stated here so the predictions are fixed):

- **PT-X (decoupling).** For any tree with gcd/zero tests mod N = pq: the run splits iff the two diagonal transcripts
  τ_p, τ_q differ, and otherwise ends in τ_p. For any map κ from primes to transcripts and deviation event
  D = {τ_p ≠ κ(p)} ∪ {τ_q ≠ κ(q)}, with p, q iid: Pr[correct or split] ≤ max_b Pr[F = b] + 2s + 2 Pr[D ∧ no split].
- **NE1 (one-layer N-exponent trees, element tests).** With κ = formal transcript (gate outputs y^{E(N)} replaced by a
  formal variable), Pr_q[node v deviates at p | p] ≤ Σ_components D_v · μ_p(o), where o is the order of the gate base
  in that component of the unit group at p and μ_p(o) = max_t Pr_q[E(N) ≡ t (mod o) | p]. A deviation at p is one of:
  SMALL (o ≤ L), UNIVERSAL (exp(G_p) | E(N) − t: every base is fooled; an archimedean pin), SUBGROUP (otherwise;
  a cross gcd coincidence gcd(E(N) − t, |G_p|) ≥ o).
- **H-GEN (conjecture, class H).** Everything the continued-fraction prefix of √N reveals about the factors, beyond the
  fiber floor and beyond splitting, is genus information: the quadratic symbols (ℓ|p) of small ℓ (incl. ℓ = −1, 2)
  constrained through ((−1)^i Q_i | p) = 1 on small Q_i, plus square relations (CFRAC) and square forms (SQUFOF).
  In particular it reveals nothing about NON-quadratic characters of p.
- **PT-J (pure trees + Jacobi symbols of pure values).** In the population N ≡ 1 (mod 4), ε(p) = ε(q), a single Jacobi
  test (c|N) = (c|p)(c|q) of an N-independent integer c carries information on ε iff β(c)β(−c) ≠ 0, where
  β(c) = E_X[(c|X)] over the prime law.

## E1: instrumented N-exponent trees (checks of proved statements + quantification)

Populations (all with p < q, both b-bit primes with the top two bits set, b = n/2):
- `iid`: all unordered pairs at n = 24 (b = 12) and n = 28 (b = 14); 200 000 uniformly sampled pairs at n = 32 (seed
  `leakmapE1|32`).
- `fermat`: pairs with 0 < q − p ≤ H = 2^{b/2}, all of them, b = 12, 14, 16.

Trees (deterministic, fixed coins; definitions in `src/leakmap/netrees.py` once written):
- `FERMAT`: element tests 2^{N−1} = 1, then (adaptive, on the zero branch only) 3^{N−1} = 1.
- `MR2`: Miller–Rabin ladder for base 2, nodes b_j = 2^{(N−1)/2^{t−j}} with tests b_j = 1 and b_j = −1, j = 0..t,
  t = v₂(N−1) (the tree branches on the public t, PT-m style).
- `LUCAS`: in Z/N[x]/(x² − 3x + 1) (α = x, norm 1): tests α^{N+1} = 1 and α^{(N+1)/2} = ±1.
- `QFROB`: in Z/N[x]/(x² − x − 3): test α^N = σ(α) = 1 − α.
- `ARCH`: nonadaptive tests 2^N = 2^k for k = 1..32 (archimedean probe).
- `REPAIR`: test (2^{N−1} − 1)·(3^{M} − 1) with M = lcm(1..64) (a pure p−1 factor times a Fermat factor).
- `PM1`: pure control 2^{M} − 1, M = lcm(1..64) (no gate; PT case).

Targets: ε = (−1/p) on N ≡ 1 (mod 4) pairs; [p ≡ 1 mod 3] on N ≡ 1 (mod 3) pairs; [q − p < X/8] (archimedean).

Registered predictions:
- **P1 (diagonal lemma):** for every pair and tree, the instrumented run mod N splits iff τ_p ≠ τ_q, and otherwise
  ends with transcript τ_p. 0 mismatches.
- **P2 (PT-X):** with κ = formal transcript and Bayes-optimal leaf labels computed exactly from the population, the
  inequality holds in every (n, population, tree, target) cell.
- **P3 (NE1 counting bound):** for every p in the population and every tree, the exact deviation frequency
  Pr_q[τ_p ≠ κ(p) | p] is ≤ Σ_{v on κ(p)'s path} Σ_components D_v·μ_p(o_v), with μ computed exactly by enumeration
  over the population's partners of p. 0 violations.
- **P4 (exact two-sided law, random bases):** for `FERMAT` and `MR2` variants with a uniformly random base per pair
  (iid at n = 28, 20 000 pairs × 20 bases), the empirical frequency of two-sided deviations (liars) matches the mean of
  gcd(E(N) − t, p−1)gcd(E(N) − t, q−1)/((p−1)(q−1))-type predictions (order-pattern lemma) within 4 SE.
- **P5 (classification, descriptive with one directional prediction):** in `iid` populations, UNIVERSAL deviations
  occur only for `ARCH` (pins q − p = k − 1) and their count equals the number of pairs with q − p + 1 ∈ [1, 32]
  hitting the pin (exactly, by the proof); in `fermat` populations the UNIVERSAL share of all deviations is larger
  than in `iid` at every b.
- Every factor found by an instrumented run is classified STRUCTURAL (a deviation or reference split explains it) or
  ACCIDENTAL; ACCIDENTAL count expected 0 for element tests (CRT criterion).

## E2: does the CF prefix leak only genus information? (test of H-GEN)

Instances: `balanced_semiprime(n, seed)` from `src/common/instances.py`, seeds 9.5·10^12 + n·10^9 + i
(disjoint from every earlier range), i < 2·10^5 at n = 32 and n = 40. Fresh-seed replication range
9.6·10^12 + n·10^9 + i, i < 2·10^5, run only for statistics that hit.

Statistics (16): the CF-prefix menu of `src/eps/cf.py::cf_features` at T ∈ {64, 1024}: n_a1, n_a2, sumbits, maxbits,
n_Qeven, n_sq_even, n_Q1mod4_odd, n_2sq_odd. Computed from N only (imported read-only from src/eps).

Targets (diagnostic labels from p, q; p = smaller factor), each on its own sub-population:
- Genus (quadratic) targets — H-GEN allows leakage:
  - G1: ε = (−1|p), on N ≡ 1 (mod 4).
  - G2: [p ≡ 1 (mod 3)] = (−3|p) indicator, on N ≡ 1 (mod 3).
  - G3: [p ≡ ±1 (mod 5)] = (5|p) indicator, on N ≡ 1 (mod 5).
  - G4: [p ∈ {1,2,4} (mod 7)] = (−7|p) indicator, on N ≡ 1 (mod 7).
- Non-genus targets — H-GEN predicts NO leakage:
  - NG1 (cubic character mod 7): [p ≡ ±1 (mod 7)], on N ≡ 1 (mod 7). Independent of every quadratic symbol of p
    under uniform residues.
  - NG2 (quartic character mod 5): [p ≡ 1 (mod 5)], on the sub-population N ≡ 1 (mod 5) and p ≡ ±1 (mod 5) (so the
    quadratic symbol (5|p) = +1 is fixed).
- Archimedean target (exploratory, no directional prediction): A1 = [q/p < 1 + 1/32].

Test: stratified permutation z for the Pearson correlation of statistic and target, strata = N mod (2^8·m)
(m = target modulus: 1, 3, 5, 7, 7, 5 for G1–NG2) and, for A1, (N mod 2^8, top 6 bits of N). Primary z at n = 32 is
the two-way cluster-robust z (clusters = p and q; prime reuse), at n = 40 the plain permutation z (also report
cluster-robust). Hit iff |z| ≥ 5.0 (112 tests per size; Bonferroni-type family-wise error ≈ 6·10^{−5} per size).
Replication on fresh seeds: same sign and |z| ≥ 3.5.

Registered predictions:
- **Q1 (replication):** at n = 32, G1 has ≥ 1 confirmed hit (eps found 9).
- **Q2 (genus reach):** at n = 32, G2 ∪ G3 ∪ G4 have ≥ 1 confirmed hit (moderate confidence; a null here is a power
  failure, not a refutation).
- **Q3 (H-GEN, the decisive one):** NG1 and NG2 have 0 confirmed hits at both n = 32 and n = 40, and their maximum |z|
  over the 16 statistics is < 4 at each size. A confirmed NG hit refutes H-GEN as stated.
- **Q4:** at n = 40, the number of genus hits is ≤ the number at n = 32 (decay).
- Positive controls (one-sided, proved): PCℓ_even: some even i ≤ 1024 with Q_i = ℓk² (gcd(k,N)=1) ⇒ (ℓ|p) = +1;
  PCℓ_odd: some odd i ≤ 1024 with Q_i = ℓk² ⇒ (−ℓ|p) = +1; ℓ ∈ {3, 5, 7}. Prediction: 0 violations.
- MDE (80% power at |z| = 5): reported for every target; NG nulls are claimed only down to that MDE.

## E3: Legendre-symbol bias of compressed pure values (search for a PT-J leak)

Primes: all primes X in [2^23, 2^24) (≈ 5·10^5). For each pure value c in the menu below compute
z₁ = √π·β̂(c), z₂ = √π·β̂(−c), where β̂(c) = mean over X of (c|X) (0 when X | c; such X are counted as splits and
reported), π = number of primes used.

Menu (N-independent SLPs; exact definitions in `src/leakmap/jacobi.py`):
- F1: 2^{2^k} + 1, k = 5..60. F2: 2^{k!} − 1, k = 5..30. F3: 3^{lcm(1..k)} + 2, k = 5..60 step 5.
- F4: V_{k!}(3) − 2 and V_{k!}(3) + 2 (Lucas), k = 5..30. F5: ∏_{x≤B}(x²+1), B ∈ {10, 30, 100, 300, 1000}.
- F6: B!, B ∈ {10, 100, 1000, 10000}; B! + 1. F7: 300 random SLPs of length 40 (ops: square, multiply by an
  earlier value, add an earlier value, add 1), seed `leakmapE3`.
- Excluded as structural: c = ±k², ±2k² (checked symbolically for the explicit families; for random SLPs by testing
  whether c mod X is a square for 64 random X — never expected).

Hit iff min(|z₁|, |z₂|) ≥ 4.5. **Prediction R1: 0 hits** (H0: every compressed pure value in the menu has
β(c)β(−c) = 0 up to noise; i.e. pure+Jacobi trees built from these values leak nothing about ε). Also report the
largest |z₁ z₂|.

## Resource limits
≤ 2 GB RSS per job, one heavy job at a time, `memory_pressure` checked before each heavy run.

## Deviations
Any deviation from this file will be listed in `notes/leakmap_results.md` under "Deviations".
