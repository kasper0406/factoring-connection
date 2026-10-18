# Final report: searching for a simple classical connection for integer factorization

**Project:** `factoring-connection` (brief: `BRIEF.md`, 17 October 2026). **Report date:** 2026-10-17 (revised the same day
after the GPT-6.1-sol report review, `codex/r5_report_review.out`; see §6).
**Bottom line:** no polynomial-time factorizer, and no asymptotic speedup over known methods, was found. The project
produced:
* a set of exact obstructions, each scoped to a restricted model;
* one conditional reduction whose missing primitive is the brief's own target;
* one classification conjecture with partial proofs;
* one partial-information interface (bits of `φ(N) mod 2^k` into Coppersmith) with an unimplemented, growing family of
  predicates. Its first two non-trivial predicates are identified with classical residuosity decisions (one also with a
  class-group 4-rank). The next level is shown not to be determined by the corresponding 2-power class-group ranks.

Status labels follow BRIEF §7: `proved`, `experimentally supported`, `conjectured`, `refuted for this model`,
`not yet tested`. Text produced by GPT-6.1-sol (via `codex exec`) is labelled as such. Where an agent did not reproduce
such a claim, it is marked "reported by GPT-6.1-sol, not independently reproduced". Where a GPT-6.1-sol note conflicts
with a later agent-verified note, the agent-verified note is used. `notes/repair_adjudication.md` is authoritative for
every claim it corrects, and the corrected form is used throughout.

**How costs are stated.** Four kinds of statement are kept apart:
1. exact operation counts of a specific implementation;
2. empirical averages over the tested instance distributions;
3. heuristic extrapolations;
4. proved bounds within a restricted model.

Only the last kind is a lower bound, and only inside its model. Measured scaling exponents are fits, not asymptotic theorems.

Notation: `F_N = (1+X)^N − 1 − X^N` over `R_N = Z/NZ`; `n = log2 N`; `N = pq` with `p < q` in the balanced family;
`S_{r,M}` is the coefficient section of class `r mod M`; `s = p + q`; `ε = (−1/p)`. The **information floor** is relative to
`N mod m`: a quantity "beyond the floor" is one not determined by `N mod m` (`trackE_analysis.md` §5). It is not an
information-theoretic statement conditional on the full integer `N`, which determines its factors.

---

## 1. Best surviving mechanism

No candidate survives as an algorithm. Every mechanism the project implemented falls into one of two groups.
* Its exact operation count grows at least linearly in `p`, in `q − p`, or in `√K` for a prefix cutoff `K ≈ √N`
  (measured on the implementations).
* Its per-attempt success probability, averaged over the tested balanced instances, is of order `1/p`. This is a measured
  average, not a uniform per-instance rate: the weighted-rank family has success probability exactly 0 on 65–86% of
  tested balanced instances.

The restricted-model bounds of §2 explain several of these measurements. None is a lower bound against all algorithms.
The three items below are what is left, ranked by how much they constrain or advance the search for a polynomial-time method.

| Rank | Item | Kind | Strongest supported claim | Why this rank |
|---|---|---|---|---|
| 1 | **φ(N) mod 2^k bit tower feeding Coppersmith** (Track E + GPT-6.1-sol r4/r5) | partial-information interface with an unimplemented, growing family of predicates | **Proved:** given `s mod 2^k`, the roots of `x² − sx + N mod 2^k` are two cosets mod `2^{k−v}` (`v = v_2(q−p)`). Whenever the exact Howgrave-Graham test `rigorous_ok` holds at modulus `2^{k−v}`, two lattice decodes factor `N`. **Experimentally supported:** the fitted threshold `k ≈ n/4 + n/(4d) − 2 + v`. **Proved:** `s mod 8` ⇔ `ε`. **Proved:** given `ε`, `s mod 16` adds nothing. **Proved:** for `N ≡ 1 (mod 8)`, `s mod 32` adds exactly `QR_2 = [2 is a square mod N]`, and `s mod 64` adds nothing beyond `QR_2`. GPT-6.1-sol's 4-rank form of `QR_2` was verified exhaustively for `p < q < 1000`. **Refuted for these invariants:** `s mod 128` is not determined by `N mod 128`, `ε`, `QR_2` and the `r4`, `r8` of `Cl(−4N)` and `Cl(−8N)` | The decoder half is complete and measured. The interface is not refuted in any model and is not order-revealing. But the decoder needs about `n/4` bits beyond the floor. Only two predicates (`ε`, `QR_2`) are identified, both are fixed-element residuosity decisions with no known polynomial algorithm, and no reduction obtains the needed prefix from them. |
| 2 | **Conditional section-SLP reduction** (escape E-1 + Agrawal–Biswas PIT) | conditional reduction | If a uniform constructor, restricted to ring operations on residue handles, answers in `poly(n)` time the child queries of the binary descent through `F_N` (`poly(n)`-size division-free SLPs over `R_N[Y]` for both children of the current section at every level), then randomized polynomial-time factoring follows (**proved**, using known PIT) | The reduction predates the project (the brief's Track A target, codex r1 card A1, AMS/Agrawal–Biswas PIT). Every tested instantiation of the primitive is refuted (gate-local growth ×3.0–3.7 per level, measured; isolating algebra observables need rank `≥ min(M, p+1)`, proved). With an unrestricted constructor the hypothesis is close to tautological. |
| 3 | **C-2R classification conjecture** (repair cycle) | partially proved obstruction | In the pre-registered run, every pass of the threshold `d ≤ 3L/8` used the character exponent `E ≡ (p+1)/2`. That exponent satisfies `(p−1) \| 2(E−1)`, so it already factors `N`. On the 16 largest instances, no non-character exponent gives any compression (`d = m + 1`). Squaring exponents provably give no constant-factor gain. Non-character exact compression below the threshold does occur at small `p` (`p = 257, E = 70`: `d = 7 < 9`; verified) | The most novel mathematics in the project (no prior formulation found by us), but it points away from factoring. The character fold's empirical ratio `d/m ≈ 3/4` (for `m ≤ 942`) is only a constant factor on a `√p`-degree object. |

### 1.1 Item 1 in detail: the φ(N) mod 2^k bit tower

**Operation.** Produce, from `N` alone, information about `s = p + q` (equivalently `φ(N) = N + 1 − s`) modulo `2^k` beyond
what `N mod 2^k` fixes.

**Decoder half (`notes/trackE_analysis.md` §1, §4.1).**
* Prop. E2 (**proved**): with `v = v_2(q − p)` and `k > 2v`, the roots of `x² − s x + N (mod 2^k)` are exactly the two cosets
  `p, q (mod 2^{k−v})`. Each coset is decoded by the Howgrave-Graham lattice. Success is **guaranteed** whenever the exact
  integer inequality `rigorous_ok` (`trackE_analysis.md` §1.2) holds at modulus `2^{k−v}` ([May10] Thms 2–3).
  Brute-force checks: 200 random pairs in the test suite.
* Fitted threshold (**experimentally supported**; 570 cells, n = 64…512, 5 seeds): decoding succeeded from about
  `k_min = n/4 + n/(4d) − 2` known low bits (instance variation about ±0.7 bit). The rigorous condition is 1–2 bits higher.
  In the `sum2k` experiment (40 instances) success occurred exactly when `k − v ≥ k_min`. That is a finite observation, not
  a theorem. At n = 512, LLL at the threshold took 0.003 s (d = 9) to 52 s (d = 65).

**Identified predicates.** Assume `p ≡ q (mod 4)` (i.e. `N ≡ 1 mod 4`) and write `p = ε + 4u`, `q = ε + 4v`.
* Exact identity (GPT-6.1-sol r4; elementary, re-derived by the writer): **`s = ε(N+1) − 16 ε u v`**.
  *Derivation:* `N = 1 + 4ε(u+v) + 16uv`, so `ε(N+1) = 2ε + 4(u+v) + 16εuv = s + 16εuv`.
* **Modulo 8.** `s ≡ 2 (mod 4)` always. `s mod 8` (its weight-4 bit) is equivalent to `ε`, i.e. to "−1 is a square mod N", i.e. to
  "N is a sum of two squares" (Prop. E2′, **proved**).
* **Modulo 16.** `s ≡ ε(N+1) (mod 16)`, so the weight-8 bit adds nothing once `ε` is known (**proved**; table verified).
* **Modulo 32.** `s ≡ ε(N+1) − 16ε·[(2/p) = −1]·[(2/q) = −1] (mod 32)` (**proved**, verified). The new information is the
  *correction* relative to the known baseline `ε(N+1)`, not the raw weight-16 bit. For `N ≡ 1 (mod 8)`, `(2/p) = (2/q)`, so the
  correction is exactly `QR_2`. For `N ≡ 5 (mod 8)` the correction is 0.
* **Modulo 64.**
  * `N ≡ 5 (mod 8)`: `s ≡ ε(N+1) − 32ε·(((r − ε)/8) mod 2) (mod 64)`, with `r ∈ {p, q}` the prime with `(2/r) = +1`
    (**proved**, verified). Again the indicator is the correction, not the raw weight-32 bit. Example: `p = 37, q = 41`
    has indicator 1, yet `s = 78` has weight-32 bit 0.
  * `N ≡ 1 (mod 8)`, with `b = [2 is a non-square mod N]`: `uv ≡ b(1 + 2[N ≡ 1 mod 16]) (mod 4)`, so
    `s ≡ ε(N+1) − 16ε(uv mod 4) (mod 64)` is determined by `(N mod 64, ε, QR_2)` (**proved**; GPT-6.1-sol r5 item 3,
    re-derived by the writer, verified).
  * *Proof of the `N ≡ 1 (mod 8)` case.* `b = 0` forces `p, q ≡ ±1 (mod 8)`, so `u, v` are even and `uv ≡ 0`. `b = 1` forces
    `u, v` odd. Then `N ≡ 1 + 4ε(u+v) (mod 16)` fixes `u + v mod 4`. `u + v ≡ 0` gives `uv ≡ −u² ≡ 3`, and `u + v ≡ 2` gives
    `uv ≡ u² ≡ 1 (mod 4)`. `N ≡ 1 (mod 16)` is the case `u + v ≡ 0`.
* **Modulo 128** (`N ≡ 1 mod 8`; needs `uv mod 8`). **Not determined** by `(N mod 128, ε, QR_2)`: 32 of 64 keys are ambiguous.
  **Not determined** after adding `r4` and `r8` of both `Cl(−4N)` and `Cl(−8N)`: 73 of 248 keys are ambiguous, covering 1 199 of
  3 418 pairs (`src/repair/verify_r5.py`; exhaustive, `p < q < 1000`). So the class-group description of the tower does not
  continue through these 2-power ranks (**refuted for these invariants**). ~~The missing quantity is the individual value at `p` of
  an order-8 Dirichlet character mod 32~~ *Corrected 2026-10-18 (`tower128`, `notes/tower_results.md`):* for `N ≡ 1 (mod 8)` the
  missing quantity is `ψ(p)` for the even quartic character `ψ = ω²` mod 16 (i.e. `p mod 16`); the order-8 value is needed only
  when `N ≡ 5 (mod 8)`. The bit is undetermined given `(N mod 128, ε, QR_2)` iff `p ≡ q (mod 16)`. Extended test (exhaustive,
  `11 ≤ p < q < 4000`): `(r4, r8)` of `Cl(−4N)` and `Cl(−8N)` together **determine** the bit in the cell `p ≡ q ≡ 7 (mod 8)`
  (Thm T3, via Rédei linearity), and 4th-power membership of the public class decides 3/8 of the cell `p ≡ q ≡ 1 (mod 8)`,
  `(p/q) = 1` (T4). In the `QR_2`-false cells every tested invariant of 6 distinct discriminants is at its null (conjecture C5).
  None of this is an algorithm: `r8` is not known to be computable without the class group, whose computation already factors `N`.
* **Class-group form of `QR_2`** (`N ≡ 1 mod 8`): `QR_2 ⇔ r4 Cl(−8N) > 0` if `ε = +1`, and `QR_2 ⇔ r4 Cl(−4N) > 0` if `ε = −1`.
  The proof (GPT-6.1-sol r4 §3) evaluates genus maps on ambiguous classes via Rédei matrices `M_{−4N}`, `M_{−8N}`. **That proof
  was not re-derived by an agent**; its consequences were verified (below).
* **Public class** `T = [(2, 2, (N+1)/2)]` of discriminant `−4N`: `T ∈ Cl² ⇔ (2/p) = (2/q) = +1` (GPT-6.1-sol; verified).

**Verification of the r4/r5 claims (`src/repair/verify_r4.py`, `src/repair/verify_r5.py`; tests `tests/test_repair_r4.py`,
`tests/test_repair_r5.py`).** Exhaustive over all pairs of distinct odd primes. Class groups are computed by brute force:
enumerate reduced primitive forms, compose by Dirichlet/Shanks composition, and count elements of order dividing 2, 4 and 8.
The form code is validated against the analytic class-number formula `h(D) = −|D|^{−1} Σ (D/a) a`, the group axioms, and
groups with known structure (`Cl(−164) ≅ Z/8`).

| Claim | Scope | Pairs `p<q<300` | Pairs `p<q<1000` | Failures |
|---|---|---:|---:|---:|
| `s = ε(N+1) − 16εuv` | `p ≡ q (mod 4)` | 902 | 6 901 | 0 |
| mod-16 table | `p ≡ q (mod 4)` | 902 | 6 901 | 0 |
| mod-32 formula | `p ≡ q (mod 4)` | 902 | 6 901 | 0 |
| `N ≡ 1 (8)`: `(2/p) = (2/q)`; the mod-32 correction equals `[2 is a non-square mod N]` | `N ≡ 1 (mod 8)` | 442 | 3 418 | 0 |
| **4-rank criterion** (`Cl(−8N)` if `ε=+1`, `Cl(−4N)` if `ε=−1`) | `N ≡ 1 (mod 8)` | 442 | 3 418 | 0 |
| `uv mod 4 = b(1 + 2[N ≡ 1 mod 16])`; `s mod 64` determined by `(N mod 64, ε, QR_2)` | `N ≡ 1 (mod 8)` | 442 | 3 418 | 0 |
| mod-64 formula | `N ≡ 5 (mod 8)` | 460 | 3 483 | 0 |
| `T ∈ Cl(−4N)² ⇔ (2/p)=(2/q)=1` | `p ≡ q (mod 4)` | 902 | 6 901 | 0 |
| `r4 = 2 − rank_{F_2} M` for both printed genus matrices | `N ≡ 1 (mod 8)` | 442 | 3 418 | 0 |
| Milović Prop. 1 as quoted (`ε=+1`: `r4 Cl(−4N) = 2 ⇔ p≡q≡1 (8), (p/q)=1`) | `ε = +1`, `N ≡ 1 (8)` | 202 | 1 569 | 0 |
| `r2 = 2` for `Cl(−4N)`, `Cl(−8N)` | as above | 902 / 442 | 6 901 / 3 418 | 0 |
| control: 4-rank criterion with the **other** discriminant | `N ≡ 1 (mod 8)` | agrees on 186/442 | agrees on 1 569/3 418 | — |
| `s mod 128` determined by `(N mod 128, ε, QR_2, r4/r8 of Cl(−4N), Cl(−8N))`? | `N ≡ 1 (mod 8)` | no: 26/165 keys ambiguous | no: 73/248 keys ambiguous | — |

The control row shows that the `ε`-dependent choice of discriminant is essential. Pairs with `p ≢ q (mod 4)` (928 of 1 830
below 300) are outside the claims' scope and were not tested. Largest class number met: `h(−8N) = 576` (`< 300`) and
`2 496` (`< 1000`). Runtimes: 0.7 s and 45 s (r4), 39 s (r5, bound 1000). Raw output: `experiments/raw/verify_r4.json` and
`experiments/raw/verify_r5.json`.

**What item 1 does not give.** It computes nothing. Deciding `ε` for `N ≡ 1 (mod 4)` is a fixed-element quadratic
residuosity problem that reduces to factoring. No converse reduction and no polynomial algorithm is known (GPT-6.1-sol r4 §1,
citing Adleman–McCurley 1994; reported, not independently reproduced). Jager–Schwenk's generic-ring equivalence concerns a
random input element. It does not cover this modulus-only problem (§5). The class-group form of `QR_2` gives no algorithm
either: computing `r4` from an unfactored discriminant is not known to be easier than computing the class group. Bits needed
beyond `s mod 64` are not identified with any predicate.

### 1.2 Item 2: the conditional section-SLP reduction

* **Statement (proved, conditional).** Run the binary descent on `F_N = (1+X)^N − 1 − X^N` (endpoints removed). Descending
  through `(1+X)^N` itself could select the index `0` or `N`, whose gcd with `N` is not proper.
  * *Assumption.* A uniform constructor, using only ring operations on residue handles, produces in `poly(n)` time
    division-free SLPs over `R_N[Y]` of size `poly(n)` for both children of the current section at every level actually
    reached. Short circuits that merely exist along some unidentified good path are not enough.
  * *Conclusion.* Randomized formal PIT over `R_N` (Agrawal–Biswas 2003, Thm 4.8; Arvind–Mukhopadhyay–Srinivasan,
    arXiv:0801.0514 §6) selects a nonzero child at each level. After `≤ n` levels the descent reaches an interior support
    index `i`, and `1 < gcd(i, N) < N` (seed lemma L2).
  * Sources: `trackA_analysis.md` §3.4 E-1, §4, §9; `prior_art_map.md` §28, §39.
* **Where zero testing is easy (odd balanced semiprimes, `M = 2^k`).** At levels with `M ≤ p − 1`, both children are nonzero
  (C1), so no test is needed. At levels with `M > p`, a nonzero section has at most 3 terms (C2) and is never a scalar false
  zero (C4).
* **Caveat on the constructor (GPT-6.1-sol r2 §1).** If the constructor may run arbitrary polynomial-time bit computations on
  `N`, any polynomial-time factoring algorithm can be compiled into it. The reduction is informative only for restricted
  constructors.
* **Evidence against the primitive.**
  * The gate-local E/O transform grows ×3.6–3.7 (naive) and ×3.0 (Karatsuba) per level on diagnostic and operational paths
    alike. This is a measurement of the implementation; at level 9 a ≤ 3-term section of a 16-bit `N` has 1.46·10^6 nodes.
  * Any isolating algebra observable needs rank `D ≥ min(M, i+1)` (proved), hence `D ≥ p + 1` at every level `2^k > p` in
    the balanced family (Corollary 3′, proved). Under a dense-storage assumption this is `Ω(p)` work per level.
  * No binomial-specific identity was found (Lucas/Chebyshev doubling, roots-of-unity filter, holonomic recurrence,
    Frobenius; `trackA_analysis.md` §6).

### 1.3 Item 3: C-2R

Precise operation: fold the interval zero set `Z_L` (`L = 2m + 1 ≈ √p`) through `φ(x) = ((x + m)^E + v)^2`, then read
`Q(φ(x))` with `deg Q = d`. The one-layer optimum is `d = m + 1 = ⌈L/2⌉`. See §2.3 (Conjecture C2R-X) and
`notes/repair_C2R_results.md`.

Strongest supported claims:
* **Pre-registered rule met.** On the 16 largest instances (`p ≥ 4099`), 0/16 non-character exponents reach `d ≤ 3L/8`. On all
  16, the non-character optimum equals `m + 1` (experimentally supported, exhaustive).
* **All 14/32 observed passes used `E ≡ (p+1)/2`** with `p ≡ 3 (mod 4)`. That exponent satisfies `(p−1) | 2E − 2`, so it
  already factors `N` (50/50 demo).
* **Non-character exact compression below the threshold exists at small `p`.** `p = 257, m = 8, E = 70, v = 186` is exact with
  `d = 7 < 9`. Here `gcd(E, p−1) = 2` and there is no `k ≤ 12` with `256 | k(E−1)` (GPT-6.1-sol r5 item 10; verified with
  two independent implementations, `src/repair/verify_r5.py` and `src/repair/c2r_search.brute_fold`). This matches the arm-3
  column of the results table (`d = 7` at `p = 257`).
* **The character fold's ratio.** `d/m ≈ 3/4` is experimentally supported for `m ≤ 942` only.

---

## 2. Mathematics

### 2.1 The unified gap theorem (formulation: GPT-6.1-sol r3 §2; components proved in `trackA_analysis.md` §3 and `trackC_analysis.md` §3 O1)

**Model.** `R` a commutative unital ring; `V` an `R`-module generated by `D` elements; `T ∈ End_R(V)`; `v ∈ V`;
`λ : V → W` an `R`-linear map into an `R`-module. Put `u_t = λ(T^t v)` and `L(f) = λ(f(T)v)`. An *observable* is exactly
this linear `L`. The state may be computed by any means, including repeated squaring, towers, CRT products and nonlinear
intermediate arithmetic. But `V, T, v, λ` stay fixed as the index `t` (or the coefficients of `f`) vary.

**Theorem (proved; determinant trick plus forward propagation plus Fitting decomposition).**
1. `u_t` satisfies a monic constant-coefficient recurrence of order `≤ D`.
2. `D` consecutive zeros force all later values to vanish.
3. The same holds along arithmetic progressions, with `T` replaced by `T^b` (nonnegative exponents, or `T` invertible).
4. Over a field, if the nilpotent part has dimension `e`, such a block forces `u_t = 0` for all `t ≥ e`. The proof gap
   (zero block starting before `e`) was found by the r3 critic and closed (`trackA_analysis.md` §3.2).
5. For scalar outputs over a field, the sharp parameter is the minimal realization dimension (Hankel rank).

A polynomial readout of degree `d` linearizes onto `C(D+d, d)` coordinates. The theorem then bounds that lifted dimension,
**not** `D` (corrected O1 corollary; counterexample in §2.4).

**Instances in this project.**
* Track A Theorem A1, corrected: an observable that isolates class `r mod M` on `[0, N]` and is nonzero on `X^i` needs
  `D ≥ min(M, i+1)`. Balanced Corollary 3′: `D ≥ min(M, p+1)`.
* Track C O1: fixed-dimension powers carrying the prefix predicate need order `≥ p + 1`; for the interval predicate,
  `≥ L + 1` (tight).
* Also: the Hadamard-mask corollary E1, fixed-quotient projected coefficient sequences, and fixed-matrix Krylov/Hankel sequences.

**Exact boundary against known algorithms** (r3 §2 table; checked against the cited algorithm descriptions by GPT-6.1-sol):

| Algorithm | Relation to the model |
|---|---|
| Pollard `p−1` | Its scalar power is *inside* (`a^t − 1` has a dimension-2 realization). But its zero set is periodic (an order condition), not an exact prefix or class mask, so the gap bound does not apply. Boundary check: `gcd(2^6 − 1, 77) = 7` with one scalar power. A "unified" theorem demanding `D ≥ 7` here would have smuggled in a false faithful-mask assumption. |
| ECM | Curve arithmetic is nonlinear in the coordinates. Point annihilation plus gcd extraction does not supply a faithful linear zero pattern. Outside. |
| BGS / Pollard–Strassen | The factorial recurrence `u_{K+1} = (K+1)u_K` has an index-dependent transition. Product trees, shifts and multipoint evaluation are outside the fixed-transition hypothesis. |
| Harvey / Harvey–Hittmeir | Collision searches over growing collections of powers plus integer arithmetic. Individual power sequences can be inside; the decoder is outside. |
| QS / NFS | Smoothness testing, integer or number-field relations and parity dependencies are not one fixed linear-power observable with the required zero pattern. Outside. |
| Coppersmith | Integer lattices and size bounds recover small roots. Modular evaluation subroutines do not capture that extraction. Outside. |

The theorem says nothing about gcds, pivots, adaptive decoders or nonlinear predicates. It is an obstruction for specific
access models, not a factoring lower bound.

For general straight-line programs, the literature offers several *differently quantified implications*, not an
equivalence with the project's uniform modular predicate:
* Shub–Smale Theorem 2: a polylogarithmic BSS procedure deciding `z ∈ {1, …, t}` would make `t!` "ultimately easy" (some
  nonzero multiples `m_t·t!` have polylog constant-free complexity).
* Lipton 2009: polylog SLPs for `t!` give polynomial-size (nonuniform) factoring circuits by gcd binary search.
* Lipton 1994: short SLPs for polynomials with many rational roots would undermine factoring hardness. The full theorem is
  unverified.
* Bürgisser (Thm 1.1 ECCC / 1.2 journal): factorial hardness implies constant-free permanent lower bounds.
* Dutta et al. Thms 30/32: implications for characteristic-zero truncation families.

An unconditional obstruction for general SLPs would require substantially stronger lower-bound results than those
established here (`trackC_analysis.md` §3).

### 2.2 Proved lemmas

"GPT" = proved by GPT-6.1-sol; the status column says whether an agent reproduced it.

| ID | Statement (short) | Status | Source |
|---|---|---|---|
| L1–L4 | First surviving coefficient `[X^p]F_N = N/p`; `k a_k = 0`; prefix test ⇔ `K ≤ p` ⇔ `(K−1)!` unit; `F_N' = 0`, Hasse derivatives keep coefficients | proved; checked all composites `N ≤ 1000` | `seed_lemmas.md` L1–L4 |
| L5 | `F_N ≡ 0` as a function on `R_N` ⇔ `N` prime or Carmichael | proved | `seed_lemmas.md` L5 |
| L6–L7 | Lucas/Kummer description of `supp F_N`; balanced `|supp| = 2q − p − 1` | proved; checked `N ≤ 1000` | `seed_lemmas.md` L6–L7 |
| G | Gap lemma (LFSR) and Fitting form (gap closed after r3) | proved; exhaustive `F_3`, `D = 3` | `trackA_analysis.md` §3.2 |
| A1 | Dimension conservation, corrected: `D ≥ min(M, i+1)`; support-aware (ii); multilinear (iii) with separate isolation hypothesis; towers/families with hypotheses attached | proved; orchestrator checked; exhaustive `F_2`, `F_3` | `trackA_analysis.md` §3.3 |
| Cor 3′ | Balanced `N`, `M = 2^k`, isolating `L` with `L(F_N) ≠ 0` ⇒ `D ≥ min(M, p+1)`; tight up to least support index | proved | `trackA_analysis.md` §3.3; `repair_adjudication.md` Part 1 |
| C1, C2 | Balanced: all classes hit for `M ≤ p−1`; ≤ 3 nonzero terms for `M > p`, `gcd(M,N) = 1` | proved; 217 772 coprime cases | `trackA_analysis.md` §4 |
| C3 | Pure classes need `M > min(p−1, q−p)`; gcd is proper when the value is nonzero mod the supporting prime | proved (necessity); threshold experimental | `trackA_analysis.md` §0, §5 |
| C4 | Odd balanced, `M = 2^k ≥ 2`: no scalar false zeros (repaired span proof) | proved; 0 in 44 110 and 0 in 5 246 | `trackA_analysis.md` §4 |
| E1 (mask) | Isolating a mask with a zero-run `g` before a nontrivial index needs rank `≥ g + 1`; field prefix/interval bounds need an observed trailing zero block | proved; 38/38 tight | `trackA_analysis.md` §7 |
| O1 | (a) prefix predicate needs order `R ≥ p + 1`; (b) interval needs `R ≥ L + 1`, tight | proved; orchestrator checked | `trackC_analysis.md` §3 |
| O2 | One power-map layer: `deg Q ≥ ⌈L/2⌉` (strengthened after r3) | proved; 1 693 434 configurations, bound attained for all 986 `(p, L)` | `trackC_analysis.md` §3 |
| O2′ | `k` layers: `#roots(Q) ≥ L/∏ gcd(e_i, p−1)`; trivial at `k = O(log L)` with bounded gcds | proved | `trackC_analysis.md` §3 |
| O3 | Single polynomial of degree `< p` needs degree `≥ L` (not `≥ p`) | proved (corrected) | `trackC_analysis.md` §3 |
| O4 | Least non-member of a proper subgroup of `F_p^*` is `< √p + 1`; for `x^K − 1` a long prefix forces `(p−1) \| K`. Much weaker than Martin–Parvardi 2022 / Burthe; not a new obstruction | proved (only for `x^K − 1`) | `trackC_analysis.md` §3 |
| M1–M5 | Random shifts carry no target information (independent uniform shifts only); power sums need `e ≥ L` and the Faulhaber evaluator is circular; `Z(N,K) ⇔ U(K)`; Jacobi prefix product reformulation; discriminant reformulation | proved (scoped) | `trackC_analysis.md` §2 |
| Tal | No step of Tal's construction gives a sub-√ unit test over `Z/NZ` with unknown factors; modular halving gives zero saving for `L ≪ N` | proved / refuted for this model | `trackC_tal_audit.md` §2, §5 |
| B1–B2 | `F_N ≡ F_m^{p^k} (mod p, f)`; readouts reduce to the vanishing pattern `(1+α)^m = 1 + α^m` | proved; 0 violations on all composites `< 5000` | `trackBD_analysis.md` §2 |
| B3 | Squarefree `N`: column-gcd elimination returns a factor ⇔ pivot profiles differ (conditional on finding unit pivots; failure bound `(1−φ(N)/N)^{64}`); profile split = rank split only for Krylov-type matrices | proved; 99 000 matrices | `trackBD_analysis.md` §2 |
| B4 | Vanishing depends on `q mod L_α`; Frobenius coincidence sufficient, not necessary | proved | `trackBD_analysis.md` §2 |
| B5 | Sporadic roots outside `U = {0, −1, ω, ω²}` form anharmonic 6-orbits (classical: Cauchy–Mirimanoff, Boyd–Martin–Thom Lemma 5.2) | proved (`p > 3`) | `trackBD_analysis.md` §2 |
| D1–D2 | Divided lifts are unit multiples of `G_k` or determined by `S`; `G_k(α) = 0` is a Fermat-type order condition; scalar `G_1` split probability `(g−1)(p+q−2(g−1))/(pq)` | proved; census agrees | `trackBD_analysis.md` §4 |
| D3 | LC of `c ↦ F_N(c)` is `e_p = q mod (p−1)` mod `p`; **new (repair):** BM over `Z/N` hits a nonunit by index `p + e_p − 1`, and any nonunit at index `≥ 2e_p` reveals `p`; Stirling window `Δ^K F_N(0) = (K+1)! S(N, K+1)` | proved; 60/60 runs, max ratio 0.47 | `trackBD_analysis.md` §4 |
| E1 (decoder) | Exact rigorous Howgrave-Graham condition (integer test `rigorous_ok`) | proved ([May10] Thms 2–3) | `trackE_analysis.md` §1.2 |
| E2, E2′ | `s mod 2^k` root cosets; `s mod 8` ⇔ sum of two squares | proved | `trackE_analysis.md` §4.1 |
| E3, E4, E5 | `s mod M` has `2^{ω'}` ambiguity, removed exactly by the unordered pair; Legendre label lemma; `p mod ℓ` CRT; `N mod (p−1)` is complete (Pollard exponent) | proved | `trackE_analysis.md` §4 |
| F1, F2 | Information floor relative to `N mod m`: a function of `(p mod m, q mod m)` is informative beyond `N mod m` iff non-constant on the fibre `{(u, c/u)}`; Jacobi symbols and CRT discriminant masks carry 0 bits beyond `N mod m` | proved (bit counts under the uniform-residue model) | `trackE_analysis.md` §5 |
| C2R-1 | Reflection criterion: `d = |A| − r*(c)`; exact folds satisfy `4r*(c) = |A| + S_A(c) + [c∈A] − 2[c/2∈A]` | proved | `repair_C2R_results.md` §3.1 |
| C2R-2 | `E ≡ (p+1)/2`, `p ≡ 3 (4)`: first layer is `y ↦ χ(y) y`, image `H ∩ T`; exactness ⇔ explicit boundary sign conditions | proved (mechanism); ratio `d/m ≈ 3/4` experimental | `repair_C2R_results.md` §3.2 |
| C2R-3 | Squaring exponents (`E ≡ 2`; `E ≡ 2 + (p−1)/2`, `p ≡ 1 (8)`): `r*(c) ≤ τ(c_0)` resp. `2τ(|c_0|)`, so no constant-factor compression, exact or not | proved; 779 + 186 primes | `repair_C2R_results.md` §3.3 |
| r4-a | `s = ε(N+1) − 16εuv`; mod-16 table; mod-32 formula; mod-64 formula for `N ≡ 5 (8)` | proved (GPT; elementary, re-derived in §1.1); verified exhaustively `p < q < 1000` | `codex_r4_sos_classgroup.md` §2; `src/repair/verify_r4.py` |
| r4-b | `N ≡ 1 (8)`: `QR_2` ⇔ `r4 > 0` of `Cl(−8N)` / `Cl(−4N)` by `ε`; `T ∈ C² ⇔ (2/p) = (2/q) = 1`; printed genus matrices give `r4 = 2 − rank M` | proved by GPT, proof not re-derived by an agent; verified exhaustively `p < q < 1000` | `codex_r4_sos_classgroup.md` §3; `src/repair/verify_r4.py` |
| r5-a | `N ≡ 1 (8)`: `uv ≡ b(1 + 2[N ≡ 1 mod 16]) (mod 4)`, so `s mod 64` is a function of `(N mod 64, ε, QR_2)` | proved (GPT-6.1-sol r5; re-derived in §1.1); verified for 3 418 pairs `p < q < 1000` | `codex/r5_report_review.out` item 3; `src/repair/verify_r5.py` |
| r4-c | One root of `−1 mod N` ↔ one primitive representation `N = a² + b²` (deterministic poly time); two roots `t ≢ ±r` factor | proved by GPT; reported by GPT-6.1-sol, not independently reproduced | `codex_r4_sos_classgroup.md` §1 |
| r2-a | Exact all-character root count `R_{p,k}(q)` (Jacobi sums); exact exposure formulas (2)–(5) for scalar, quotient, resultant and rank-gap events | proved by GPT; 65 exact checks in GPT-6.1-sol's own script `codex/round2_experiments.py`; not independently reproduced by an agent | `codex_r2_generic_ring.md` §2 |
| r2-b | Nonadaptive ring+gcd splitting ⇔ a sampled observable with `Pr(Z_p ⊕ Z_q) ≥ 1/poly`; constant-only SLP splitter needs `L ≥ s + O(1)`; restricted permanent barrier | proved by GPT; reported by GPT-6.1-sol, not independently reproduced | `codex_r2_generic_ring.md` §1, §3 |
| r1-F1 | Norm-one torus with exponent `N − (D/N)`: split probability `d/A + d/B − 2d²/(AB)` | proved by GPT; reported by GPT-6.1-sol, not independently reproduced | `codex_r1_mechanisms.md` F1 |

### 2.3 Conjectures and open questions

| ID | Statement | Status | Evidence | Source |
|---|---|---|---|---|
| C2R-X | For every `ε > 0` there are `p_0, K` such that for `p > p_0`, an exact two-layer fold (both gcds 2) with `≥ εL` internal pairs has `(p−1) \| k(E−1)` for some `k ≤ K` | conjectured | D2 = 0 at every instance `p ≥ 4099`; no twisted identity gave an exact fold with even one pair; proved for squaring exponents. Small-`p` non-character compression exists (`p = 257, E = 70`), consistent with the existential `p_0` | `repair_C2R_results.md` §3.5 |
| C2R ratio | Best exact character fold has `d/m → 3/4` | experimentally supported for `m ≤ 942`; asymptotic form conjectured | `m = 70 … 942`: 0.720–0.787 | `repair_C2R_results.md` §3.2 |
| Twisted linearity | Twisted identities (`k ≥ 3`) give `≈ m/k²` non-exact pairs | experimentally supported; proof open (Weil error too large at length `√p`) | table §3.4 | `repair_C2R_results.md` §3.4 |
| BD-(iii) | Apart from public roots and Frobenius coincidences, public algebraic relations among sketch values hold with probability `O(#conditions/|F_{p^d}|)` | conjectured | every B/D experiment; measured per-sketch split rate ≈ 2.3–2.6/p on average | `trackBD_analysis.md` §5 |
| O4-analogues | Lucas (`p+1`) and elliptic analogues of O4 | conjectured; the elliptic version must use point orders, not `#E` | — | `trackC_analysis.md` §3 |
| E-1 | Binomial-specific poly-size SLPs for `S_{r,2^k}(F_N)` | open (not yet tested beyond the gate-local transform and the identities in §6) | — | `trackA_analysis.md` §3.4, §9 |
| Shub–Smale mod p | Polylog division-free SLPs over `Z[x]`, uniform in `L`, vanishing mod `p` exactly on `Z_L` for all `p` in a dyadic range | open | — | `trackC_analysis.md` §3 |
| Track E X | No mechanism of Tracks A–D yields a non-floor coset bit at sub-`N^{1/4}` cost | conjectured / model-conditional | rests on known-algorithm costs | `trackE_analysis.md` §0, §6 |
| r4 hardness | Fixed-character decisions (`ε`, `QR_2`, `ψ(p)`) are hard on suitable semiprime distributions | conjectured by GPT-6.1-sol | — | `codex_r4_sos_classgroup.md` §1–2 |
| Tower via 2-power ranks | `s mod 128` is determined by `(N mod 128, ε, QR_2)` and `r4`, `r8` of `Cl(−4N)`, `Cl(−8N)` | refuted for these invariants (73/248 keys ambiguous, `p < q < 1000`) | `experiments/raw/verify_r5.json` | `src/repair/verify_r5.py` |
| Tower via other public class groups | `s mod 128` determined by 2-power ranks of `Cl(−4kN)`, small `k` | partly: determined in cell p≡q≡7 (8) by (r4,r8) of Cl(−4N),Cl(−8N) (T3, proved); 3/8 of cell p≡q≡1 (8),(p/q)=1 (T4); QR_2-false cells at null (C5, conjectured) | `notes/tower_results.md` | `src/tower/` (2026-10-18) |
| Coarse predicates | Track C obstructions for predicates correct only when `\|K − p\| > N^{1/4−ε}` | not yet tested | — | `trackE_analysis.md` §6 |
| C-HHloc | Harvey–Hittmeir restricted to a certified interval near `√N`, cost `Õ(W^{3/5}N^{−1/10})` (GPT-6.1-sol "proved adaptation") | not yet tested; reported by GPT-6.1-sol, not independently reproduced; not polynomial in any case | — | `trackC_candidates.md` C-HHloc |

### 2.4 Counterexamples (every one reproduced with exact arithmetic by the `repair` agent unless noted)

| Claim refuted | Counterexample | Source |
|---|---|---|
| A1 without the index hypothesis; "selector needs rank `2^{k+1}`" | `N = 35, M = 16, r = 5`: rank-6 observable "coefficient 5 of `R_35[X]/(X^6)`" isolates the class and returns `a_5 = 7` | adjudication A2 |
| A-3 "rank-`D` selection is a coin flip" | `f ↦ f(2)` over `F_101` separates `1` and `X` but isolates neither class; `BudgetSelector` was a literal coin | adjudication A4 |
| C2 for non-coprime `M` | `N = 35, M = 7, r = 0`: terms 7, 14, 21, 28 | adjudication A6 |
| C3 "pure class always gives a factor" | `N = 35, M = 8, r = 1`: section `21Y^3`, gcd 35 at `c = 0` | adjudication A7 |
| C4 as printed; even `N` | `N = 6`: `S_{0,2} = 3Y + 3Y²` is the zero function on `Z/6` | adjudication A8 |
| "Scalar traps only for Carmichael at `M = 1`" | `S_{0,4}(F_51) = 34Y³ + 34Y⁶ + 17Y⁹ + 17Y¹²` vanishes on `Z/51`; 744 traps among odd squarefree `N ≤ 1500` | adjudication A9 |
| Field prefix bound without trailing block | `F_5`, `u_t = 1 − 2^t`: values `0,4,2,3,0`, rank 2 `< 4` | adjudication A10 |
| Circuit growth "proved" | `X^N` chain has one nonzero coordinate; `(1+X)^{2^j}` mod `(35, X^8 − 3)` has 0–5 zero coordinates | adjudication A11 |
| Rank ⇒ bit cost | `f = X^N`: rank-4096 isolating observable `c^{⌊N/M⌋}[N ≡ r]` costs 40 multiplications | adjudication A12 |
| O1 polynomial-readout corollary | `α` primitive in `F_{p²}`, `z_K = α^K`, degree-4 norm readout zero exactly at `K = p+1, p+2` (`p ∈ {7, …, 257}`); interval: state `(A, 1)` with readout `∏(A+i)` | adjudication C1 |
| O3 "escape needs degree `≥ p`" | Pochhammer `∏_{i<10}(x+i)` at `p = 101`: degree 10, zero set exactly `Z_10` | adjudication C6 |
| O2′ consequences | `gcd(lcm(1..10), 256) = 8` vs 10-smooth part 256 (`p = 257`); differ for 536 of 2 261 primes `< 20 000` | adjudication C4–C5 |
| O4 elliptic analogue with `#E` | `y² = x³ + 2x + 3` over `F_101`: `#E = 96`, `(3, 6)` has order 12 | adjudication C7 |
| M2 "`e < p−1` needed" | `p = 7, e = 10`: direct expansion is a degree-`≤ e` polynomial | adjudication C9 |
| B3 universal `2^{−64}` | `N = 105`: `(1 − φ/N)^{64} ≥ 1.0·10^{−17}` | adjudication B1 |
| Rank split ⇔ profile split in general | `[3 5]` over `R_15`: rank 1 both sides, pivots `{1}` vs `{0}`, returns 3 | adjudication B2 |
| "Every non-Frobenius count ≡ 0 mod 6" (literal definition) | `p = 7, q = 11`: roots `2, 4 = ω, ω²`, count 2; fails in 1 039/5 260 pairs | adjudication B3 |
| Degree independence of `P_gap` | `N = 209`, `c = −1`: `P_gap = 0` at `r = 2`, `1/9` at `r = 6` (extension roots) | adjudication B4 |
| M3b necessity (point coincidence) | `N = 247`, `f = X³ + 5`: `Q_N = diag(1, −1, 1)` mod 13, rank split | adjudication B5 |
| `Δ^K` zero predicate = prefix predicate | `N = 391`: gcds 391, 17, 391 at `K = 12, 13, 14` | adjudication B7 |
| Model M has a lower bound | `N = 77`, constant matrix `[2^6 − 1] = [63]` returns 7 | adjudication B9 |
| "Non-character exponents have `o(L)` pairs" without exactness | twisted identities `(p−1) \| k(E−1)`, `k = 4, 6, 8`: `≈ m/k²` pairs | `repair_C2R_results.md` §3.4 |
| "Only the character exponent compresses" (first version of this report; true only on the 16 large instances / at the `3L/8` threshold) | `p = 257, m = 8, E = 70, v = 186`: exact, `d = 7 < 9`, not order-revealing for `k ≤ 12` (GPT-6.1-sol r5; verified by `writer`) | `src/repair/verify_r5.py` |
| "Class-group ranks continue the bit tower" (first version of NEXT_STEPS) | `s mod 64` needs no rank (r5-a); `s mod 128` is not determined by `r4`, `r8` of `Cl(−4N)`, `Cl(−8N)` (73/248 ambiguous keys) | `src/repair/verify_r5.py` |
| "The indicator is the raw new bit" | `p = 37, q = 41` (`N ≡ 5 mod 8`): mod-64 indicator 1, weight-32 bit of `s = 78` is 0; the indicator is the correction relative to `ε(N+1)` (GPT-6.1-sol r5; checked by hand) | `codex/r5_report_review.out` item 4 |
| Uniform degree-2 rank-gap guarantee | `N = 101·113, 107·193, 1009·1481`: `P_gap = 0` for all `r ≤ 16` (GPT-6.1-sol r1; reproduced by `sketchBD`) | `trackBD_analysis.md` §6.5 |
| Character-order restriction in r2 S3 | `R_{17,1}(19) = 2`, `R_{17,1}(41) = 5` with equal gcds (reported by GPT-6.1-sol, not independently reproduced) | `codex_r2_generic_ring.md` §2 |
| SOS ⇔ `−1` QR for arbitrary integers | `45 = 3² + 6²`, `−1` not a square mod 45 (GPT-6.1-sol; elementary) | `codex_r4_sos_classgroup.md` §1 |

### 2.5 Remaining proof obligations

1. Conjecture C2R-X. The missing tool is a bound on `|A_E ∩ (c − A_E)|` or on `S_{A_E}(c)` for sets of size `√p`.
2. The `3/4` ratio in Prop. C2R-2 needs `Σ_{a∈I} χ(a(c−a)) = o(|I|)` on intervals of length `√p`. This is beyond
   Pólya–Vinogradov; whether Burgess covers it was not checked.
3. The r4 4-rank criterion: an agent-checked proof (currently GPT-6.1-sol's genus-matrix argument plus exhaustive verification).
4. B/D (iii): a proof of the `1/|F_{p^d}|` sporadic rate, or a counterexample family.
5. The `p+1` and elliptic analogues of O4 (with point orders).
6. A lower-bound model for Track B/D that excludes Pollard `p−1` (restricted readout degree and constants); none is proved.
7. The standard-model reduction question for `ε` on an explicit balanced distribution (Section 5).

---

## 3. Experiments

### 3.1 Environment and conventions
Python `.venv` (CPython 3.13.5; sympy 1.14, gmpy2 2.3.2, python-flint 0.9.0, numpy, fpylll 0.6.4). Exact integer or
`Z/NZ` arithmetic throughout; floating point only in fits. Candidates receive only `N`, public configuration and a seed;
diagnostic and harness code is separate (`AGENT_RULES.md` §2). Flint audit: `nmod_poly` division by a nonunit leading
coefficient aborts the process, so operational code never calls flint division (`trackC_analysis.md` §1).

### 3.2 Input distributions
* `balanced_semiprime(bits, seed)` (`src/common/instances.py`): `p, q` independent uniform primes of `⌈n/2⌉, ⌊n/2⌋` bits
  **with the top two bits set**, resampled until `p ≠ q` and `N` has exactly `n` bits.
  * For even `n` (every size used in the experiments) both primes lie in `[3·2^{b−2}, 2^b)`, so **`q/p < 4/3`** (measured max
    1.31 over 200 seeds at 16, 24, 32 bits). For odd `n` the halves differ by a bit and `q/p` reaches about 2.6.
  * This is an RSA-style restriction. Consequences: `e_p = q mod (p−1) = q − p + 1`, `e_p/p ≤ 0.33`, and `q mod (p−1)` is
    not equidistributed (`trackBD_analysis.md` §3).
  * At 16 bits only about 66 distinct semiprimes of this shape exist, so 16-bit rows contain repeats;
    `balanced_semiprime(16, 0)` and `(16, 1)` coincide.
* Other families: unbalanced semiprimes, prime powers, 3- and 4-prime products, Carmichael numbers, primes,
  `diag_near_equal` (diagnostic), `diag_bigg` (`q = 1 + t(p−1)/s`, diagnostic, Track B/D).
* Exhaustive ranges:
  * all composites `N ≤ 1000` (seed lemmas);
  * all balanced `N ≤ 5000` (Track A census);
  * all composites `4 ≤ N < 5000` (Track B/D tiny run);
  * all primes `5 ≤ p < 100`, `< 140` and `< 20 000` (Track C O2, E5 and E6);
  * 32 fixed primes 257–24 631 (C-2R);
  * all odd prime pairs `< 300` and `< 1000` (r4/r5).

### 3.3 Reproduction commands (from the notes; run from the project root)

```bash
.venv/bin/python -m pytest -q                                   # full suite: 168 tests, ~8 s
.venv/bin/python src/seed_checks.py                              # BRIEF §9 smoke test + extended checks to N <= 1000

# Track A (notes/trackA_analysis.md §8; < 10 CPU-min)
.venv/bin/python -m src.trackA.run_circuits --bits 16 24 32 --seeds 2 --levels 9 --out experiments/raw/trackA_circuits.jsonl
.venv/bin/python -m src.trackA.run_descent --bits 16 20 24 28 32 --seeds 5 --max-rank 1048576 --out experiments/raw/trackA_descent.jsonl
.venv/bin/python -m src.trackA.zero_tests --maxN 5000 --out experiments/raw/trackA_zero_tests.jsonl
.venv/bin/python -m src.trackA.theory_checks --out experiments/raw/trackA_theory.jsonl
.venv/bin/python -m src.trackA.run_masks --out experiments/raw/trackA_masks.jsonl
.venv/bin/python -m src.trackA.analyze                           # -> experiments/analysis/trackA_summary.md

# Track C (notes/trackC_analysis.md §4; baselines ~10 min, E8 ~5 min)
.venv/bin/python -m src.trackC.exp_baselines --config experiments/configs/trackC_baselines.json
for e in E1 E2 E3 E4 E5 E6 E7 E8; do .venv/bin/python -m src.trackC.harness_mechanisms --exp $e; done
E8_VARIANT=no_zero_exp .venv/bin/python -m src.trackC.harness_mechanisms --exp E8
E8_VARIANT=gcd_le2     .venv/bin/python -m src.trackC.harness_mechanisms --exp E8

# Track B/D (notes/trackBD_analysis.md §6; ~2-3 core-hours)
.venv/bin/python -m src.trackBD.exp_tiny --limit 5000 --R 200 --M 100 --procs 12 --out experiments/raw/trackBD_tiny.jsonl
.venv/bin/python -m src.trackBD.exp_scalar --bits 16 20 24 28 32 --inst 200 --K 4000 --procs 12 --out experiments/raw/trackBD_scalar.jsonl
.venv/bin/python -m src.trackBD.exp_scalar --bits 36 40 --inst 60 --K 4000 --procs 12 --out experiments/raw/trackBD_scalar.jsonl
.venv/bin/python -m src.trackBD.exp_scalar --family diag_bigg --bits 24 32 --s 1 2 4 8 16 32 --inst 30 --K 4000 --out experiments/raw/trackBD_scalar_diag.jsonl
.venv/bin/python -m src.trackBD.exp_scalar --family carmichael --inst 33 --K 2000 --out experiments/raw/trackBD_scalar_carm.jsonl
.venv/bin/python -m src.trackBD.exp_expect --bits 8 10 11 12 --count 12
.venv/bin/python -m src.trackBD.exp_rankdist --procs 8
.venv/bin/python -m src.trackBD.exp_rank --bits 16 20 24 28 32 40 --inst 150 --procs 8
.venv/bin/python -m src.trackBD.exp_bm --bits 16 20 24 --inst 30 --procs 8
.venv/bin/python -m src.trackBD.exp_sketch_bal --bits 16 20 24 28 32 --inst 100 --K 100 --procs 10
.venv/bin/python -m src.trackBD.exp_sketch_bal --bits 36 40 --inst 30 --K 100 --procs 10
.venv/bin/python -m src.trackBD.analyze > experiments/analysis/trackBD_tables.txt

# Track E (notes/trackE_analysis.md §8; threshold ~25 min on 16 workers)
.venv/bin/python -m src.trackE.exp_threshold experiments/configs/trackE_threshold.json
TRACKE_WORKERS=8 .venv/bin/python -m src.trackE.exp_formats pm1 sum2k legendre highnoise sumM crtnoise hsbits flips cost
.venv/bin/python -m src.trackE.info_floor
.venv/bin/python -m src.trackE.analyze                           # -> experiments/analysis/trackE_tables.md

# Repair: r3 audit reproductions and the pre-registered C-2R cycle (notes/repair_C2R_results.md §5; ~1.45 CPU-h total)
.venv/bin/python -m src.repair.audit_trackA  --out experiments/raw/repair_audit_trackA.jsonl
.venv/bin/python -m src.repair.audit_trackC  --out experiments/raw/repair_audit_trackC.jsonl
.venv/bin/python -m src.repair.audit_trackBD --out experiments/raw/repair_audit_trackBD.jsonl
.venv/bin/python -m src.repair.c2r_search --procs 32 --out experiments/raw/repair_C2R.jsonl
.venv/bin/python -m src.repair.c2r_mixing --procs 32 --out experiments/raw/repair_C2R_mixing.jsonl
.venv/bin/python -m src.repair.c2r_theory --procs 16 --out experiments/raw/repair_C2R_theory.jsonl

# r4 / r5 verification (0.7 s at bound 300 and 45 s at bound 1000 for r4; 39 s at bound 1000 for r5)
.venv/bin/python -m src.repair.verify_r4 --bound 300 --out experiments/raw/verify_r4.json
.venv/bin/python -m src.repair.verify_r5 --bound 1000 --out experiments/raw/verify_r5.json
```

Raw results: `experiments/raw/*.jsonl` and `*.json` (two files gzipped). Configs: `experiments/configs/`. Aggregated tables:
`experiments/analysis/`. GPT-6.1-sol prompts, raw outputs and its own scripts: `codex/`.

### 3.4 Key measurements and scaling diagnostics

Every factoring success below is a correctness check, not evidence of polynomial scaling. The last column labels the kind of
statement: **count** (exact operation count of the implementation), **avg** (empirical average over the tested
distribution), **fit** (exponent fitted to finite data) or **heur** (heuristic extrapolation).

| Track | Measurement | Key numbers | Diagnosis |
|---|---|---|---|
| A | Hash-consed E/O section circuits | Unique nodes at level 0/3/6/9 (32-bit, naive): 97 / 2 563 / 119 819 / 6 411 324. Ratio ×3.72 per level (naive), ×3.04 (Karatsuba); all five paths agree within 3% | count + fit: growth ≈ `3^k`–`4^k` per level for this transform |
| A | Exact quotient selector, 25 runs, 16–32 bits | 25/25 paths alive at every level; the first scalar always certifies. Cost ratio exactly 2.00 per level (32-bit level 20: 4.7·10^7 rank-weighted products, 4.8 s) | count: rank `2^{k+1}` at level `k`, so `Θ(N log N)` to full depth for this selector |
| A | Random descent | paths die 1–3 levels after `bitlen(p)` (25/25), matching support density `≈ 2/p` | avg |
| A | Pure-class gcds | 338 STRUCTURAL vs 3 ACCIDENTAL (null ≈ 2.05); first hit at `2M ∈ (q−p, 4(q−p)]` in 25/25 | count + avg: Fermat-type, `≈ (q−p) log N` |
| A | Scalar zero-test census | 0 false zeros in 44 110 nonzero sections (all balanced `N ≤ 5000`) | exhaustive |
| A | Prefix mask | nonzero ⇔ `K > p` in 36/36 (24–48 bits); 48-bit rank `2^20`: 7.9·10^7 products, 9.6 s | count: linear in rank |
| C | Interval-product baselines, 64-bit `N`, `K ≤ 2^30` | Log-log slopes (wall / exact Karatsuba `rn_mul` / FFT model): naive 1.00 / 1.000 / –; Pollard–Strassen 0.56 / 0.815 / 0.617; BGS 0.52 / 0.763 / 0.563. Identical at 128 bits | fit: consistent with the known `Õ(√K)` |
| C | End-to-end least factor by BGS `U(K)` + binary search | 13/13 instances (balanced 32/48/64, unbalanced, `p²`, `p³`, 3- and 4-prime, Carmichael 6601, a prime). 64-bit: 33 calls, 29 s | count: `Θ̃(N^{1/4})` (known BGS bound) |
| C | Truncated `F_N` prefix test (M3) | 0 mismatches vs `U(K)` for composites `< 600` plus 4 Carmichael; at 64-bit, `K = 2^16`: 11 000× BGS cost | count: linear in `K` |
| C | Random-shift probes (M1) | hit rate 0.067 (`I` true) vs 0.077 (false); null 0.062–0.077 | avg: no signal |
| C | O2/O4 exhaustive | O2: 0 violations, all primes `< 140`; strengthened bound attained for all 986 `(p, L)`, `p < 100`. O4: max ratio 0.69 for `p < 20 000`; least generating prefix ≤ 29 | exhaustive |
| B/D | Theorem B1 | 0 violations on every irreducible component of every sketch of all composites `< 5000` | exhaustive |
| B/D | Scalar sketches, `K = 4000` weights, 16–40 bits | The exact criterion predicts hit counts within Poisson error at every size. About 70% of hits are public (`c ≡ 0, −1, ω`). Non-public rate ≈ null (32 bits: 33 SPOR vs 27 NULL). No FROB hit at ≥ 28 bits | avg: order `1/p` per sketch on average; per-instance rate is exact and often 0 |
| B/D | 99 000 rank-split matrices | "factor ⇔ profile split": 0 exceptions; M1 ranks = B1 prediction in 19 800/19 800. Split counts fall ≈ ×4 per +4 bits; 0 at 40 bits | avg |
| B/D | Exact rank distributions | mean `P_gap·p ≈ 2.3–2.6` at 24 bits; `P_gap = 0` exactly for 65–86% of balanced instances | avg (exact per instance) |
| B/D | BM over `Z/N` (`c ↦ F_N(c)`) | `LC_p = e_p` in 90/90; first structural nonunit at `2e_p` (58) or `2e_p − 1` (1); 31/90 met an accidental nonunit first | count: `≈ 2(q − p)` terms (proved bound `p + e_p − 1`); a weak Fermat-gap method |
| B/D | Degree-`r` sketches, 16–40 bits | SPLIT counts track the prediction (e.g. 24-bit `r = 2`: 13 vs 10). Per-sketch split rate ≈ 10^{−5} at 32 bits; 0 at 36–40 bits | avg |
| E | Coppersmith threshold, n = 64…512, d = 3…65 | `k_min = n/4 + n/(4d) − 2 + δ`, `δ = −0.05 ± 0.72` (low bits). Decode 0.003 s (d = 9) to 52 s (d = 65) at n = 512 | fit; polynomial guarantee via `rigorous_ok` |
| E | Formats | `φ(N) mod 2^k`: success iff `k − v ≥ k_min` (40 instances). `N mod (p−1)`: 80/80. `s mod M`: decodes ≈ `2^{ω'}/2`. Legendre bits alone: list `2^{105}` at n = 512. Random known bits, density `δ ≥ 0.5`: 40/40 | avg (finite) |
| C-2R | Pre-registered two-layer folding, 32 primes | Arm 1 / 2 / 3: 0/32, 14/32, 0/32 (= audit). 16 largest: non-character optimum = one-layer degree. X1: 0/16. Best exact character fold `d/m ∈ [0.72, 0.79]` for `m = 70 … 942`; non-character `d = 7 < 9` at `p = 257` | exhaustive per instance; ratio fit |
| r4/r5 | Exhaustive class groups and 2-adic identities | §1.1 table: all r4 claims hold; wrong-discriminant control agrees on only 42% / 46%; `s mod 64` determined without ranks; `s mod 128` not determined by `r4`/`r8` | exhaustive |

### 3.5 Failures and implementation defects found

* `BudgetSelector` returned `rng.getrandbits(1)` once `2M > D`, so its "coin-flip" measurement was by construction (r3 audit A4).
* The first Euclid readout lacked the gcd check on the final constant remainder and missed every `r = 2` split; fixed and
  regression-tested (`trackBD_analysis.md` §2 B3).
* Flint `nmod_poly` division by a nonunit leading coefficient aborts the process (not catchable); worked around.
* Track E's static leakage test is a substring audit only; operational modules were also read by hand.
* GPT-6.1-sol r2's CFRAC pilot rule found 0/4 factors at 24 and 32 factor bits (reported by GPT-6.1-sol, not independently reproduced).
* The first version of NEXT_STEPS proposed a mod-64 class-group experiment that was algebraically vacuous (r5 review item 3).
* Every mechanism card in `notes/*_candidates.md` ends in "refuted for this model", "correct but exponential", or "not yet tested".

---

## 4. Comparison with prior art

Novelty standard: we write "no prior formulation found by us" at most. Most results are reformulations or bookkeeping of
known facts. The prior-art map (`notes/prior_art_map.md`) is GPT-6.1-sol output. Sources marked "agent-read" were checked
against the primary PDF by an agent.

| Result | Closest prior method or fact | What changed | Novelty verdict |
|---|---|---|---|
| Seed lemmas L1–L7 | AKS Lemma 2.1 proof; OEIS A020639 (Sykora 2017); Lucas, Kummer | Exact residue `N/p`, exact balanced support count, Carmichael function-trap characterization | Known / bookkeeping |
| Unified gap theorem, A1, O1, Cor 3′, mask corollary | LFSR theory, determinant trick, Massey 1969 | Applied to coefficient access, quotient sketches and prefix predicates; the balanced bound `D ≥ min(M, p+1)` | Bookkeeping in new language (r3 §3) |
| C1–C4 section structure | Lucas; linear independence of characters | Exact sparsity and no-false-zero facts for binary sections of balanced `N` | Known tools; useful instance facts |
| Conditional section reduction | BRIEF Track A; Agrawal–Biswas 2003 (Thm 4.8); AMS 0801.0514 §6; codex r1 A1 | Measured the gate-local transform and bounded algebra-observable selectors | Prior formulation exists |
| Track C baselines and M3 | BGS 2007 (Thm 5, 8, Lemma 12–13; agent-read); Pollard–Strassen; Lipton 2009 | Instrumented re-implementation | Known |
| M3 truncation | Dutta–Jindal–Pandey–Sinhababu CCC 2021 Thm 30/32 (agent-read) | Modular `F_N` truncation equivalence; no hardness transfer | Known equivalence; transfer not proved |
| M4 Jacobi prefix product, M5 discriminant | codex r1 C1 / Idea 2; resultant formulas | Reformulations of the factorial-unit predicate | Reformulation |
| O2 one-layer bound `deg Q ≥ ⌈L/2⌉` | Power-map fibres; `∏(x² − i²)` halving; sparse root bounds (Kelley; Canetti et al.) | Interval vs multiplicative-fibre moment argument | **No prior formulation found by us** (GPT-6.1-sol r3: "did not find its exact formulation") |
| O2′, O3 | Fibre composition; root counting | — | Elementary |
| O4 | Pollard `p−1`; Burthe; Martin–Parvardi 2022 Cor. 3 (agent-read) | Weaker than prior art | Not new |
| Tal audit | Tal ECCC TR26-211 (agent-read) | Step-by-step transfer table to unknown modulus: no step survives | Audit, no new mechanism |
| B1–B2 exact sketch criterion | Frobenius; AKS quotient rings | Readouts reduce exactly to `(1+α)^q = 1 + α^q` at public points | Bookkeeping |
| B3 elimination ⇔ profile split | Unit-pivot elimination folklore; codex r1 B1 | Exact iff for squarefree `N`, with corrected failure bound | Elementary |
| B4 Frobenius vs sporadic | Pollard / Bach–Shallit order coincidences | `q mod L_α` dependence; Frobenius sufficient not necessary | Bookkeeping |
| B5 six-orbits | Cauchy–Mirimanoff; Helou; Boyd–Martin–Thom 2015 Lemma 5.2 (agent-read) | Cross-prime polynomial `F_e mod p`; census data | Symmetry classical; census new data only |
| D1–D2 lifted derivatives | Fermat-witness gcds; Miller–Rabin on Carmichael; Grantham Frobenius test | Lifts carry nothing beyond `G_k` and `S` | Known phenomenon |
| D3 BM bound `p + e_p − 1` | Massey's lemma; Fermat's method | Proved nonunit index for BM over `Z/N` on `F_N(c)` | Elementary application |
| E1 decoders, E5 | Coppersmith 1997 Thm 4/5; Howgrave-Graham; May 2010 survey (agent-read); Pollard `p−1` | Calibrated price list: fitted `k_min ≈ n/4 + n/(4d) − 2` | Known |
| E2 `φ(N) mod 2^k` cosets | Coppersmith low bits; 2-adic Hensel | Exact coset structure with `v_2(q−p)` loss | Elementary |
| E2′, r4, r5-a bit tower | Fermat two-squares; genus theory; Rédei matrices (Stevenhagen arXiv:1806.06250); Milović Prop. 1; Adleman–McCurley 1994 | Identifies the new information in `s mod 8`, `16`, `32`, `64` with `ε`, nothing, `QR_2` (⇔ `r4 > 0`), nothing; shows `r4`/`r8` do not determine `s mod 128` | Each ingredient known; **the packaging as a tower of predicates feeding Coppersmith: no prior formulation found by us** (largely GPT-6.1-sol's derivation; search not exhaustive) |
| E3, F1–F2 information floor | CRT; Hittmeir's hyperbolic sieve; Stănică–Mulder–Hittmeir 2026 | Exact statement that functions of `N mod m` carry 0 bits beyond `N mod m`; masks of codex E1/E2 carry none | Elementary |
| C-2R (Lemma C2R-1, Prop C2R-2, Thm C2R-3, Conj C2R-X) | Quadratic-character folding; Pollard `p−1`; sums of two squares | Exact reflection criterion; classification of threshold-passing compression as order-revealing | **No prior formulation found by us** |
| r2 exposure formulas, nonadaptive equivalence | Aggarwal–Maurer 2009; Jager–Schwenk 2013; Brown 2010 (GPT-6.1-sol reports) | Exact CRT-disagreement formulation | GPT-6.1-sol output; not independently checked |

**Items with no prior formulation found by us:**
1. the O2 one-layer interval bound and its `⌈L/2⌉` strengthening;
2. the C-2R package (reflection criterion, character-fold mechanism, squaring-exponent theorem, Conjecture C2R-X);
3. the packaging of the low bits of `φ(N) mod 2^k` as a tower of predicates (`ε`, `QR_2` ⇔ 4-rank positivity, then
   undetermined by 2-power ranks at `mod 128`), as an interface to Coppersmith.

All three are small. None is a factoring advance. Our literature search was limited to the sources listed in `references.md`.

---

## 5. Continuation: the single most decisive next question

**Correction to the first version.** The first version posed this question in a "ring-plus-gcd (generic ring)" model. That
was not well-posed (GPT-6.1-sol r5 items 1, 2, 8).
* With `N` public and no hidden input element, every handle can be tracked by its canonical representative. Ring
  operations, equality and gcd tests can then be simulated in polynomial time, so a "generic" `ε`-decider is just an
  ordinary one. Jager–Schwenk's result re-randomizes a hidden input element and does not apply to a modulus-only problem.
* The model also does not contain every mechanism of the project. It contains:
  * quotient sketches and section evaluations (Tracks A, B);
  * elimination, Euclid and Berlekamp–Massey over `Z/N`;
  * power-map folds, `G_k` evaluations and BGS interval products.
* It does not contain:
  * computation mod `N²` followed by division by `N` (Track D). D1 shows that the specific quotients constructed equal
    unit multiples of `G_k`, which is computable mod `N`, but the lift itself is not an `R_N`-algebra operation;
  * Jacobi symbols, integer representatives, continued fractions, lattice reduction and integer root finding (Track E).
* Finally, "a proof would close the program" did not follow. A factoring-equivalence transfers hardness. It does not rule out
  a polynomial-time factoring algorithm, and it does not show that every partial-information format must compute `ε`.

**The question, restated in the standard model.** Let `D_n` be the output distribution of `balanced_semiprime(n, ·)`
(`n` even, uniform seeds), conditioned on `N ≡ 1 (mod 4)`. On this distribution `ε = (−1/p) = (−1/q)` takes both values.

> **Is there a randomized polynomial-time algorithm `A` with `Pr_{N←D_n}[A(N) = ε(N)] ≥ 1/2 + 1/poly(n)` above the best
> predictor that uses only `N mod 2^k`?**

**Why it is decisive for this project: a gate.** Among the coset-shaped formats Track E priced, `φ(N) mod 2^k` is the cheapest
(about `n/4 + v + O(1)` bits, ≤ 2 decodes), and `ε` is the first predicate it needs beyond the floor. Any algorithm that
produces `s mod 8` on `D_n` predicts `ε` (Prop. E2′, exact). Hence:
* **Negative answer.** If `ε` is not predictable on `D_n` (for instance, if factoring reduces to `ε`-prediction on `D_n`),
  then the tower route is no easier than factoring itself. That conclusion concerns this route only; other routes are untouched.
* **Positive answer.** It would be one predicate. The decoder needs about `n/4` bits. The next predicates are `QR_2` (at
  `mod 32`) and an unidentified one at `mod 128`, which 2-power ranks of `Cl(−4N)` and `Cl(−8N)` do not determine (§1.1).
  A positive answer is necessary for the route, not sufficient.

**Status.** Open. It is the average-case form of the sum-of-two-squares decision for semiprimes. That decision reduces to
factoring; the converse is a long-open question (Adleman–McCurley 1994, as reported by GPT-6.1-sol r4, not independently
reproduced). The project should not expect to settle it by proof. Hence the attack is mostly empirical, with each part able
to change the plan.

**Concrete attack.**
1. **Correlation battery with a positive control and a power analysis** (NEXT_STEPS item 1).
   * *Menu, frozen in advance.* Track A section values, Track B sketch coefficients, ranks and residues, BM lengths, Jacobi
     symbols of sketch values, continued-fraction **prefixes** of `√N`, and counts of small reduced forms of `Cl(−4N)`.
   * *Run.* `10^5` instances of `D_n` at `n = 32, 48, 64`.
   * *Positive control.* "The period of `√N` closes within `T` steps and is odd." That implies negative Pell is solvable,
     hence `ε = +1`. It is a known one-sided signal, not leakage. The run is valid only if the battery detects it at its
     measured frequency.
   * *Power analysis.* At `S = 100` statistics the corrected threshold is `|z| ≥ 6.0`; the minimum detectable correlation at
     80% power is `|ρ| ≈ 0.022`.
   * *Interpretation.* A null result excludes effects of that size for that menu only. It does not show that `ε` is
     unpredictable.
2. **The `s mod 128` level** (NEXT_STEPS item 2). Test whether 2-power ranks of `Cl(−4kN)` for small public `k` determine
   `s mod 128`. This decides whether any public class-group data continues the tower past `QR_2`. It costs minutes.
3. **Reduction map** (NEXT_STEPS item 3). Record exactly where known self-reductions fail for the fixed element `−1` on `D_n`.
   Distinguish two things: refuting a particular extraction lemma (by an explicit oracle on which its construction fails), and
   a computational separation, which such examples can neither establish nor refute.

---

## 6. Process and review

**Who did what.**

| Work | Done by |
|---|---|
| Orchestration, ownership, proof checks of A1 and O1/O2/O4, curation of GPT-6.1-sol notes, STATUS log, placement of this report | main Claude session (orchestrator) |
| Shared infrastructure, seed lemmas, Track A (sections, A1, C1–C4, selectors, masks) | Opus 5.5 agent `infraA` |
| Track C (baselines, M1–M8, O1–O4, Tal audit) | Opus 5.5 agent `factC` |
| Track B/D (B1–B7, D1–D3, sketch experiments) | Opus 5.5 agent `sketchBD` |
| Track E (decoders, formats, information floor) | Opus 5.5 agent `trackE` |
| r3 adjudication (36 items), audit reproductions, pre-registered C-2R repair cycle | Opus 5.5 agent `repair` |
| Final deliverables, `src/repair/verify_r4.py`, `src/repair/verify_r5.py` and their tests, this revision | Opus 5.5 agent `writer` |
| Round 1: six mechanism cards (A1, B1, C1, D1, E1, F1), prior-art map, Tal/BGS transfer analysis, phase-1 diagnostics | GPT-6.1-sol (`codex exec`, web search) |
| Round 2: generic-ring synthesis corrections, exposure formulas, cards E2/E3/H2/F2 | GPT-6.1-sol |
| Round 3: adversarial critic audit of Tracks A, B/D, C; unified gap theorem; C-2R design | GPT-6.1-sol |
| Round 4: sum-of-two-squares bit, higher bits of `s`, class-group cards F1–F3 | GPT-6.1-sol |
| Round 5: review of the first version of this report (14 items) | GPT-6.1-sol |

**Independent reviews that actually occurred.**
* Orchestrator proof checks: Theorem A1 (infraA) and O1, O2, O4 (factC) (STATUS log).
* Independent read-only Opus critic on Track B/D: 15 findings, all incorporated (`trackBD_analysis.md` header).
* Independent read-only Opus critic on Track E: rigorous condition, determinant, E2, E2′, Legendre lemma, F2 counts, drift
  arithmetic and leakage separation judged correct; seven fixes applied (`trackE_analysis.md` §7a).
* GPT-6.1-sol r3 adversarial audit of Tracks A, B/D and C, adjudicated item by item by the `repair` agent with exact
  reproductions (`src/repair/audit_track*.py`, 22 tests). **36 WRONG/GAP items; the critic was right on 36/36**: claim
  changed in 27, scope or wording clarified in 9, critic wrong on 0. Every number the critic gave reproduced exactly.
* `factC` checked GPT-6.1-sol's r1 Tal/BGS claims against the Tal and BGS PDFs (all spot checks confirmed; Harvey/HH details
  not read by the agent).
* `sketchBD` reproduced GPT-6.1-sol r1 claims (a) cyclotomic factors, (c) `G_1` split probability, and the `P_gap` values for
  `101·113`, `107·193`, `1009·1481`, `47959`.
* `writer` reproduced the GPT-6.1-sol r4 numeric claims exhaustively (`src/repair/verify_r4.py`; §1.1). GPT-6.1-sol had
  reported checks on 238 prime pairs below 200; the agent check covers 6 901 pairs below 1000.
* GPT-6.1-sol r5 review of the first version of this report: 14 items. The orchestrator agreed with all 14. `writer`
  verified the two computational items independently: item 3 (`uv mod 4` formula, 3 418 pairs, 0 failures) and item 10
  (exact non-character fold at `p = 257`, two implementations). All 14 items are addressed in this revision; the index is in
  `notes/review_and_repairs.md`.
* GPT-6.1-sol r1 and r2 describe internal critic/auditor passes within the same GPT-6.1-sol run. Those are not counted here as
  independent reviews.
* Pre-registration: the C-2R decision rule, instances and frozen extension X1 were committed (`a8920b7`) before any instance
  run. It was a replication pre-registration (the critic had already run arms 1–3); X1 was new.

**Not reviewed independently:** GPT-6.1-sol r2's theorems (exposure formulas, nonadaptive equivalence, constant-only SLP bound,
permanent barrier), r1's F1 torus formula, r4's genus-matrix proof (only its consequences were checked numerically), the
GPT-6.1-sol prior-art map entries not re-read by an agent, and the HH certified-interval adaptation.

The chronological review and repair index is `notes/review_and_repairs.md`.
