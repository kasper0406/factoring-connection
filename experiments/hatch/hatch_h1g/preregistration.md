# Pre-Registration: Hatch 1 — Buium p-Derivations, Int(Z) Lattices, and Galois-Folded CRT Decoding

Timestamp: 2026-11-20 (written strictly prior to writing or executing any experiment script).
Workspace: `/Volumes/git/factoring-connection/codex/scratch/hatch_h1g/`

## Hypotheses Under Test

1. **H1 (`Int(Z)` and Buium `p`-Derivation Capacity Wall — Theorem Int-J & Theorem Buium-Lin):**
   - Replacing standard Coppersmith/Boneh polynomials `Z[Y]_{<= L}` by the full Pólya ring `Int(Z)_{<= L} = span_Z { binom(Y, k) : 0 <= k <= L }` (which contains all polynomials in Buium's iterated `p`-derivations `delta_p^{(u)}(Y)` of total polynomial degree `<= L` for all primes `p`) satisfies the exact index lower bound `[Int(Z)_{<= L} : I^*_{P, r, L}] >= [Int(Z)_{<= L} : J_{P, r, L}] >= P^{r(r+1)/2}` (`1 <= r <= L + 1`), and at `r = 1` with squarefree `P` satisfies `[Int(Z)_{<= L} : I^*_{P, 1, L}] = prod_{q | P} q^{floor(L/q) + 1}`.
   - Empirically, LLL decoding on `Int(Z)_{<= L}` does not beat Boneh/Howgrave-Graham's Johnson-bound agreement radius by more than finite-degree `ln L / ln H` terms (`<= 10%` at `H = 2^12, L = 4`), while linear Buium jet interpolation `c_0 + sum_{u=0}^{s-1} c_{u+1} delta_2^{(u)}(j) == 0 (mod P')` performs strictly *worse* than `r = 1` Coppersmith `(1, j, ..., j^s)` for `s >= 3` due to degree growth `2^u`.

2. **H2 (Local `q`-Derivation Redundancy at `w = 0` & CRT-Gluing Height Blowup):**
   - At `z` with `v_q(z) = m >= 1`, `v_q(delta_q^{(u)}(z)) == m - u` for all `0 <= u <= m`.
   - Across coprime primes `q <= 19` (`M = prod_{q <= 19} q = 9,699,690`), the symmetric CRT lift `Y_1 in (-M/2, M/2]` of `(delta_q(j) mod q)_{q <= 19}` for `j in [2, 256]` has no small-height representative (`min_{2 <= j <= 256} |Y_1| > 256`, median `|Y_1| / M` in `[0.15, 0.35]`).

3. **H3 (Galois-Folded `O_K` Decoding Beyond Johnson vs Smooth-Norm Split-Prime Failure):**
   - **(a) Full-orbit beyond-Johnson decoding:** Over `K = Q(zeta_5)` (`m = 4`, cyclic Galois `sigma: zeta_5 |-> zeta_5^2`, `disc(K) = 125`), with `N = 11` rational primes `q_i` (`103 <= q_i <= 163`, `W = 4 sum ln q_i approx 212.5` nats) and messages `alpha = sum_{i=0}^3 a_i zeta_5^i` (`a_i in [-60, 60]`, `b = 4 ln 121 = 19.18` or `a_i in [-350, 350]`, `b = 4 ln 701 = 26.21`, rate `R approx 0.123` to `0.24`), specifically at `a_i in [-1800, 1800]` (`b = 4 ln 3601 = 32.76` or `a_i in [-25000, 25000]`, `b = 4 ln 50001 = 43.28`, `R = 43.28 / 212.5 = 0.204`, Johnson fraction `sqrt(R) = 0.451`), wait — set `a_i in [-60000, 60000]` so `b = 4 ln 120001 = 46.78`, **`R = 46.78 / 212.5 = 0.220`, Johnson threshold `sqrt(R) = 0.469` (`> 5/11 = 0.4545`)**, with `5` of `11` primes agreeing (fraction `5/11 = 0.4545 < sqrt(R)` and `< 0.50` unique decoding): `s = 3` Galois-folded LLL (`rank 16`) satisfies `c_0 + c_1 alpha + c_2 sigma(alpha) + c_3 sigma^2(alpha) == 0` in `K` with `dim_Q ker <= 2` in `>= 80%` of 30 random trials, whereas `s = 1` recovers `0%`!
   - **(b) Primitive smooth-norm failure:** For all non-zero primitive vectors `a in [-6, 6]^4` (`gcd(a_0, a_1, a_2, a_3) = 1`) in `Z[zeta_5]` with `163`-smooth norm `|N_{K/Q}(alpha)|`, the full-orbit agreement `A_orbit = 4 sum_{q <= 163 : q O_K | alpha} ln q` equals `0` in 100% of instances.

4. **H4 (Toy Smooth-Search Wall on Quadratic Norms `Q(k, 1)`):**
   - On 40 planted smooth-norm instances of `Q(k, 1) = k^2 + N_0` / reduced discriminant `-4 N` forms (`N` 28-bit semiprimes, `B = 47`, `H = 512`), standard Boneh `Z[k]` (`L = 4, r = 2`), `Int(Z)` (`L = 4, r = 2`), and linear Buium `(1, Q(k,1), delta_2(Q(k,1)))` all certify `0/40` witnesses without trial-dividing cofactors when run on the full factor-base modulus `P = prod_{q <= 47} q^{e_q}` (`ln X << sqrt(ln H * ln P)`).

## Pre-Registered Pass/Fail Criteria

- **P1a**: 100% of `(r, L)` cells satisfy `[Int(Z)_{<= L} : I^*_{P, r, L}] >= [Int(Z)_{<= L} : J_{P, r, L}] >= P^{r(r+1)/2}`, and `r = 1` matches `prod_{q | P} q^{floor(L/q)+1}` exactly.
- **P1b**: Empirical minimum decoded `ln d` over 30 planted trials (`P = prod_{q <= 23} q`, `H = 2^12`, `L = 4, r = 2`) for `Int(Z)` (`J` and `I^*`) is `>= 0.85 *` standard `Z[Y]` Boneh and `> 0.70 * sqrt(ln H * ln P)`.
- **P1c**: Linear Buium `delta_2` (`s = 3` and `s = 4`, `H = 2^10`, `P = prod_{q <= 29} q`) requires strictly larger minimum `ln d` than `r = 1` Coppersmith `(1, j, ..., j^s)`.
- **P2**: `v_q(delta_q^{(u)}(z)) == m - u` holds in 100% of tested `(q, m, u)` tuples; `min_{2 <= j <= 256} |Y_1| > 256` and `median(|Y_1|/M) in [0.15, 0.35]`.
- **P3a**: At `5/11` agreeing rational prime ideals (`A/W <= 0.46 < sqrt(R)`), `s = 3` Galois-folded LLL over `Z[zeta_5]` yields an exact annihilating affine operator `E(alpha) == 0` (`dim_Q ker <= 2`) in `>= 80%` of 30 trials, vs `0%` for `s = 1`.
- **P3b**: `A_orbit == 0` in 100% of primitive `163`-smooth elements in `[-6, 6]^4 \ {0}` in `Z[zeta_5]`.
- **P4**: `0/40` full-modulus smooth-search recoveries for `Z[k]`, `Int(Z)`, and linear `delta_2` at 28-bit `N`, `B = 47`, `H = 512`.
