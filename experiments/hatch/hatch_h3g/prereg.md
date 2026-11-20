# Pre-Registration: Hatch 3 (R / Thin Dependencies, Shared Descent DAGs, Cross-Input Reuse)

Frozen before any experimental code is written or executed.
Date: 2026-11-20
Workspace: `/Volumes/git/factoring-connection/codex/scratch/hatch_h3g`
Global RNG seed: `20261120`
Rules: Exact integer / rational / GF(2) arithmetic for all algebraic, rank, discriminant, and splitting checks (no factor leakage in operational routines).

## Experiment E1: Thin Dependencies, 1-Coordinate Z-Route Splitting (Lemma PIN-1), and Support Growth (Lemma PIN-STRIP)
- **Instances:** 20 deterministic seeded 40-bit balanced semiprimes `N = p * q` (`p, q` primes in `[2^19 + 1, 2^20)`): 10 Blum (`p ≡ q ≡ 3 mod 4`) and 10 general odd semiprimes.
- **Operational relation family:** Factor base `FB` of primes `≤ B = 600` (`π(600) = 109`). Degree-1 two-linear-form relation generator on `λ = q_0^H mod N` (with `H = 84 > 2 log2 N` so diagonal dominance / non-tautological wrap holds, or random `B_0`-smooth multipliers `λ = ∏_{q ≤ 19} q^{e_q} mod N`), Gauss-reduced to `(c1, d1), (c0, d0)` and evaluated on primitive `(a, b)` to collect certified `B`-smooth multiplicative relations `∏_{q ≤ B} q^{v_q} ≡ ±1 (mod N)`.
- **Registered Predictions:**
  - **P1a (Lemma PIN-1 splitting):**
    (i) On all 10 Blum semiprimes, given ANY single integer relation pinning the smallest prime `q_0` with Jacobi `(q_0 / N) = -1` (`q_0^D ≡ 1 mod N`, `D > 0`), Miller's power chain splits `N` in `10/10` instances.
    (ii) Across `10 × 10 = 100` trials of uniform random units `b ∈ (Z/N)^×` on the 10 general semiprimes, pinning `b^D ≡ 1 mod N` splits `N` via Miller's chain in `[45, 85]` of `100` trials.
  - **P1b (No thin-row escape from column cancellation):** Let `h_F2` be the number of chronological `B`-smooth relations until the first nontrivial `GF(2)` parity dependency on `FB`, and `h_Z1` the number of relations until the first nonzero `Q`-combination cancelling all primes `q > q_0` (pinning `q_0` up to sign). Prediction: `h_Z1 > 54` (`> π(B)/2`) in `20/20` instances, and `h_Z1 / h_F2 ≥ 0.80` in `≥ 18/20` instances.

## Experiment E2: Shared Descent DAGs and Sublattice Discriminant Pooling (Theorem DAG-POOL)
- **Instances:** 10 deterministic seeded 60-bit balanced semiprimes `N`. Monic base-`m` cubic `f(X)` (`m = floor(N^(1/3))`, `f(m) = N`), product form `W(a, b) = (a - b m) F(a, b)` of total degree `D = 4`.
- **Sublattices:** For each `N`, find the first `s = 64` primes `q_j > B_0 = 500` admitting a root `f(r_j) ≡ 0 (mod q_j)` with `q_j ∤ disc(W)`, form the index-`q_j` sublattice `Λ_j = {(a, b) : a ≡ r_j b (mod q_j)}`, Gauss-reduce to basis `M_j` (`det M_j = q_j`), and form the stripped integer form `G_j(u, v) = W(M_j(u, v)) / q_j`.
- **Registered Predictions:**
  - **P2a (Exact discriminant scaling):** On `40/40` tested pairs (`10` instances × `4` sublattices each), exact polynomial resultant/discriminant arithmetic gives `disc(G_j) == q_j^6 * disc(W)` (`q_j^((D-1)(D-2))` at `D = 4`), `0` violations.
  - **P2b (Theorem DAG-POOL vs naive union loss):** Enumerate primitive `(u, v) ∈ [-32, 32]^2` (`u >= 0`) on `s = 1` vs pooled `s = 64` sublattices, sorting positive `log(|G_j(u, v)|)`. Whereas the naive union bound (`lowerbound_1_codex_proof.md:73`) subtracts `(D/2) ln(64) = 2 ln(64) ≈ 8.318` nats at rank `r = 200`, `Theorem DAG-POOL` predicts `λ_pool(200; s=64) - λ_pool(200; s=1) >= -2.50` nats in `10/10` instances.
  - **P2c (DAG rank conservation without Assumption T):** Across `30` constructed shared-eliminator hypergraphs/DAGs (star degrees `m ∈ {2, 5, 20}` and 2-level shared diamonds over `B_0`), verify `0` violations of `net_B0_rows <= distinct_fresh_relations` and `B0_cofactor_relations >= net_B0_rows` when the `> B_0` submatrix has full column rank `|L|`.

## Experiment E3: Cross-Input Reuse, Abelian Tied Subfields (Theorem TIED-PARITY), Candidate C-RM, and Multiplier Height Conservation (Lemma SNFS-MULT)
- **Part E3a (Theorem TIED-PARITY exact incidence ranks):** Construct the exact `{0, 1}` matrix whose rows are all `V * p_0` hyperplane indicators `{x ∈ F_{p_0}^r : <w, x> = c}` across the `V = (p_0^r - 1)/(p_0 - 1)` degree-`p_0` subfields of `G = (Z/p_0 Z)^r` for `(p_0, r) ∈ {(2,2), (2,3), (2,4), (2,5), (3,2), (3,3), (5,2)}`.
  - **P3a:** Over `GF(2)` at `p_0 = 2`, `rank_F2 == r + 1` in `4/4` cases (`3, 4, 5, 6`). Over `GF(2)` at odd `p_0 ∈ {3, 5}` (`3/3` cases) AND over `Q` in all `7/7` cases, `rank == p_0^r` (`τ = 0`, zero column compression).
- **Part E3b (Candidate C-RM and Lemma SNFS-MULT on 10 seeded 48-bit semiprimes):**
  - **P3b (Candidate C-RM identity and norm slope):** For monic base-`m` cubic `f_0(m) = N` (`d = 3`) and `r = 4` planted linear elements `δ_i = a_i + θ` (`u_i = isqrt(m) + i`, `a_i = u_i^2 - m`), verify `ψ_0(∏_{i ∈ w} δ_i) ≡ (∏_{i ∈ w} u_i)^2 (mod N)` in `150/150` non-empty `w ⊆ {1..4}` (`10` instances × `15` subsets), and verify that the OLS slope of `ln|N_{F/Q}(Δ_w)| / ln N` against Hamming weight `|w| ∈ {1, 2, 3, 4}` lies in `[0.45, 0.75]` in `10/10` instances (matching theoretical `1/2 + 1/(2d) = 2/3` up to finite-size leading-coefficient constants).
  - **P3c (Lemma SNFS-MULT multiplier conservation):** On the same `10` 48-bit semiprimes (`d = 3`), scan multipliers `1 ≤ k ≤ 1024`, setting `m_k = floor((k N)^(1/3))` and `H_k = max_{0 ≤ j ≤ 2} c_{j,k}` in the base-`m_k` expansion `k N = m_k^3 + c_{2,k} m_k^2 + c_{1,k} m_k + c_{0,k}`. Prediction: `min_{1 ≤ k ≤ 1024} (ln m_k + ln H_k) / ((2/3) ln N) >= 0.85` in `10/10` instances.
