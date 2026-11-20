# Pre-Registration: Chase Lemma 1 (H / Branch Selection)

Timestamped before generating any fixtures or running any experiment.
Workspace: `/Volumes/git/factoring-connection/codex/scratch/ch_L1g/`
Interpreter: `python3` (exact integer arithmetic only; zero external dependencies).

## 1. Hypotheses & Pre-Registered Predictions

- **P1 (Exact Bipartite Identity, 40-bit exhaustive check):**
  On 4 unplanted 40-bit semiprimes (`B = 128, K = 32`, deterministic seeds `910000 + rep`), exhaustive traversal of all good-split-prime-power CRT branches `(d, r mod d)` with `2K^2 < d <= 2 B K^2` satisfies the exact identity
  `#winning_branches(Q) == sum_{(u,v) in W_good(Q)} #{d | Q_good(u,v) : 2K^2 < d <= 2 B K^2}`
  with 0 discrepancies across all 4 instances (`4/4` exact).

- **P2 (Fixed-Work Main-Shell Comparison at 72, 84, 96 bits):**
  Across 5 bit lengths `(bits, B, K) in [(48, 256, 256), (60, 384, 384), (72, 512, 512), (84, 1024, 768), (96, 2048, 1536)]` (6 unplanted semiprime fixtures per size, deterministic seed `920000 + 1000*bits + rep`; notice `2K^2 > B^2` at every size so every branch has `>= 3` prime factors), compare all branch selectors at an identical charged work budget `W_budget = 4000` elementary operations per form:
  - `RND_BRANCH`: uniform random good-prime branch until `d > 2K^2`, Euclid reconstruction, smooth test on fresh primitive `(u,v) in [1,K]^2`.
  - `RND_UPPER`: uniform random branch restricted to upper-band primes `q in (B/4, B]`.
  - `PROXY_TILT`: importance-sampled branch choosing primes with weight `q^{-0.85}` (Dickman small-prime tilt).
  - `BEAM_SV`: width-4 shuffled-batch beam pruning by Euclidean Gauss-reduced coordinate sup-norm `max(|u1|, |v1|)` (`round2_crt` mechanism).
  - `BEAM_QNORM`: width-4 shuffled-batch beam pruning by exact reduced partial quadratic norm `min_{Lambda_{d_k,r_k}\setminus 0} Q(u,v)/d_k`.
  - `MITM_COLLIDE`: split-base meet-in-the-middle intersecting reconstructed `(u,v)` across left/right branch tables (`d_L, d_R > 2K^2`), testing collisions first and spending residual budget on non-colliding reconstructed candidates.
  - `MCMC_PROXY`: Metropolis-Hastings single-prime swap chain over `d > 2K^2` branches scored by small-prime cofactor log-gcd `log2 gcd(Q(u,v)//d, prod_{p<=31} p^8)` at temperature `T = 2.0`.
  - `TREE_SIEVE`: penultimate-branch (`2K^2/B < d' <= 2K^2`) sublattice enumeration + prime line sieving, charging every setup op, lattice point, and prime incidence.
  - Controls: `TINY_BOX_CTRL` (uniform primitive sampling in `[1, K_eff]^2` matched to `BEAM_SV`'s 90th-percentile output coordinate `K_eff` on that form) and `BOX_SIEVE` (random-row sieve charging setup + cells + incidences).
  **Prediction P2:** Pooled over the three largest sizes (`72, 84, 96` bits) in the main shell `max(u, v) > K // 4` (`93.75%` of the box, where `1 - o(1)` of asymptotic witnesses reside), **no** proxy/MITM/beam/MCMC selector (`PROXY_TILT`, `BEAM_SV`, `BEAM_QNORM`, `MITM_COLLIDE`, `MCMC_PROXY`, `TREE_SIEVE`) achieves `>= 2.0x` the fresh smooth witness count of `RND_BRANCH` at fixed work `W_budget = 4000`, and `MITM_COLLIDE` achieves `<= 1.0x` `RND_BRANCH`.

- **P3 (Tiny-Box Control Dissection of `BEAM_SV`):**
  Unconditional `BEAM_SV` (including `max(u, v) <= K // 4`) achieves `<= 1.25x` the pooled fresh smooth witness count of `TINY_BOX_CTRL` across `72, 84, 96` bits at equal work `W_budget = 4000`, confirming any raw per-norm-test enrichment comes from coordinate shrinkage into `[1, K_eff]^2` rather than branch-level cofactor prediction.

- **P4 (Partial-Product Non-Predictiveness / Witness Non-Visibility):**
  Across pooled in-box branches (`max(u, v) > K // 4`, up to 200 sampled valid branches per form at `72, 84, 96` bits), the Pearson correlation `r` between first-step exact reduced partial log-cofactor `-log2(min_{Lambda_{q1,r1}\setminus 0} Q / q1)` and terminal residual log-smoothness fraction `log2(B-smooth part of Q(u,v)//d) / log2(Q(u,v)//d)` satisfies `|r| < 0.10`.

- **P5 (Exact CRT-Scramble Verification):**
  Across all distinct primitive box-point pairs `(u1, v1) != (u2, v2)` lying on the same penultimate sub-branch `(d', r')` (`d' > 2K^2 / B`) encountered in `TREE_SIEVE` across all 30 instances, **0** pairs vanish on the same root branch modulo any good split prime `p in (B // 2, B]` with `p nmid d'` (`0` violations).
