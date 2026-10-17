# Track C candidate cards (BRIEF §4 template)

Owner: factC. Date 2026-10-17. Evidence, proofs and reproduction commands are in `notes/trackC_analysis.md`. Throughout, `L = B − A`, `p` is the least prime factor, and `Z_L = {0, −1, …, −(L−1)} mod p`.

---

## C-B0 — BGS evaluation-shifting interval product (reference baseline)

| Field | Content |
|---|---|
| **ID and one-sentence mechanism** | C-B0: compute `∏_{A≤j<B} j mod N` by maintaining values of `g_d(x) = ∏_{j<d}(A + v x + j)` on `x = 0..d`. Lagrange shifting (BGS Thm 5) doubles `d`. |
| **Target problem / supported input family** | `I(A,B)` and `U(K)`; any N. |
| **Representation and how it is built from N** | `d+1` values in Z/NZ; all integers are public. |
| **Exact operation that extracts information** | `gcd(product, N)`; then Lemma-12 descent to a witness `j`. |
| **Why that information exposes a factor** | Direct: the product is a nonunit iff some `j` shares a prime with N. |
| **Known lemmas being used** | BGS Thm 5, Thm 8, Lemma 12, Lemma 13. |
| **New lemma required** | None. |
| **Correctness argument or conjecture** | `proved`; tests pass. |
| **Cost in n, including construction and extraction** | `O(M(√L))` ring ops plus `O(log L)` gcds. Measured slope 0.52 (wall), 0.563 (FFT model), 0.763 (exact Karatsuba). `Θ̃(N^{1/4})` per prefix query near `√N`. **Not polynomial.** |
| **Randomness and success probability** | Deterministic. |
| **Closest prior method and the actual difference** | It is BGS. |
| **Smallest informative experiment** | Done (`exp_baselines`). |
| **Counterexample or scaling behavior that would refute this version** | n/a |
| **Possible repair if that happens** | n/a |

---

## C-M1 — Random-shift transfer

| Field | Content |
|---|---|
| **ID and one-sentence mechanism** | C-M1: hope that nonunit status of `∏_{i<L}(s+i)` at random shifts `s` (or random linear combinations of such values) informs `I(A, A+L)`. |
| **Target problem / supported input family** | `I(A,B)`; semiprimes. |
| **Representation and how it is built from N** | Values `P_L(s)` mod N computed by BGS. |
| **Exact operation that extracts information** | `gcd(P_L(s), N)`, or `gcd(Σ c_s P_L(s), N)`. |
| **Why that information exposes a factor** | Claimed: correlation across shifts. |
| **Known lemmas being used** | CRT. |
| **New lemma required** | — |
| **Correctness argument or conjecture** | **Refuted** (`proved`): `p \| P_L(s)` iff `s mod p ∈ Z_L`, independent of `A`. Linear combinations are `0 mod p` with probability `≈ 1/p` unless all terms are. |
| **Cost in n, including construction and extraction** | `√L` per probe. |
| **Randomness and success probability** | Per probe `≈ L/p + L/q`, the accidental-gcd null. |
| **Closest prior method and the actual difference** | Random-interval probing, i.e. birthday-free trial division. |
| **Smallest informative experiment** | E1: hit rate 0.067 when `I(A)` is true and 0.077 when false; null 0.062–0.077. All 69 hits are explained by the probe's own interval. |
| **Counterexample or scaling behavior that would refute this version** | Observed: no conditional dependence. |
| **Possible repair if that happens** | Shifts correlated with `A` modulo the unknown `p` would need `p`. Shifts `A + tN` are trivial (same residues). None found. |

---

## C-M2 — Power-sum / Faulhaber unit test

| Field | Content |
|---|---|
| **ID and one-sentence mechanism** | C-M2: test `gcd(S_e(A,B), N)` with `S_e = Σ_{A≤j<B} j^e`, computed in `poly(e, log N)` by Faulhaber with Bernoulli numbers mod N. |
| **Target problem / supported input family** | `I(A,B)`. |
| **Representation and how it is built from N** | `B_0..B_e mod N` by the recurrence `Σ_k C(n+1,k)B_k = 0`, with gcd-checked `1/(n+1)`. |
| **Exact operation that extracts information** | `gcd(S_e, N)`; also failed denominators. |
| **Why that information exposes a factor** | Claimed: the interval-dependent value mod `p`. |
| **Known lemmas being used** | Faulhaber; von Staudt–Clausen (denominators); Fermat (`e ≡ 0 mod p−1` gives a count). |
| **New lemma required** | — |
| **Correctness argument or conjecture** | **Refuted** (`proved`). For fixed `L`, `S_e(A, A+L) mod p` is a polynomial in `A` of degree `≤ e`, so it has at most `e` zeros; detection needs zeros exactly on `Z_L`, hence `e ≥ L`. But degree `e ≥ L` Faulhaber needs `1/k` for `k ≤ e+1`, i.e. a unit test of `(e+1)!`: circular, with cost `≥ e ≥ L`. For `e ≡ 0 mod (p−1)` (unknown), `S ≡ L − #multiples ≠ 0`. |
| **Cost in n, including construction and extraction** | `O(e²)` ring ops; useful only for `e ≥ L`. |
| **Randomness and success probability** | Zero rate `≤ e/p` regardless of `I` (E2: exact count over all residues, 8 primes; zeros never exceed `e`). |
| **Closest prior method and the actual difference** | — |
| **Smallest informative experiment** | E2. The odd-`e` bias inside `Z_L` is explained: centered odd power sums are divisible by the center, i.e. `gcd(center, N)`. STRUCTURAL-TRIVIAL (BRIEF §2.3). |
| **Counterexample or scaling behavior that would refute this version** | Shown. |
| **Possible repair if that happens** | Fermat-quotient sums `Σ q_p(j)` (Lerch: `≡ W_p`) need `p`. Mod-N analogues need `φ(N)`. No repair found. |

---

## C-M3 — Truncated `F_N` prefix test

| Field | Content |
|---|---|
| **ID and one-sentence mechanism** | C-M3: `Z(N,K) = [(1+X)^N ≡ 1 mod (N, X^K)]`, which equals `U(K)`, evaluated by repeated squaring of truncated series. |
| **Target problem / supported input family** | `U(K)`; all composite N (E3: 0 mismatches for `N < 600` and 4 Carmichael numbers). |
| **Representation and how it is built from N** | Length-`K` series. |
| **Exact operation that extracts information** | First nonzero coefficient index `k`; then `gcd(k, N)`. |
| **Why that information exposes a factor** | BRIEF §2.1–2.4. |
| **Known lemmas being used** | `[X^p]F_N = N/p`. |
| **New lemma required** | A sub-linear truncation representation. |
| **Correctness argument or conjecture** | `proved` equivalent to `U(K)`. |
| **Cost in n, including construction and extraction** | `O(M(K) log N)`: slope 1 in `K`, about 11 000× BGS at `K = 2^16` (64-bit N). |
| **Randomness and success probability** | Deterministic. |
| **Closest prior method and the actual difference** | Dutta et al. CCC'21 Thm 30/32 (read): short constant-free circuits for related truncations imply easy `n!` and factoring. |
| **Smallest informative experiment** | E3. |
| **Counterexample or scaling behavior that would refute this version** | Linear scaling (observed). |
| **Possible repair if that happens** | Would need a short circuit for the truncation, which is conditionally as hard as factoring (Dutta Thm 30). Cross-track: Track A/B. |

---

## C-M4 — Jacobi-symbol prefix product / coprime count (from codex r1 idea 2)

| Field | Content |
|---|---|
| **ID and one-sentence mechanism** | C-M4: `C_N(K) = ∏_{j<K}(j/N)` is 0 iff `U(K)` fails. Equivalently, the coprime count `Σ_{j<K}(j/N)² < K−1`. |
| **Target problem / supported input family** | `U(K)` and `I(A,B)` for odd N. |
| **Representation and how it is built from N** | Jacobi symbols, each computable in `poly(n)`. |
| **Exact operation that extracts information** | A fast incomplete principal-character count. |
| **Why that information exposes a factor** | Direct. |
| **Known lemmas being used** | Jacobi reciprocity; Gauss sums of the quadratic character mod N (they only re-express single symbols). |
| **New lemma required** | A sub-`√K` evaluator of `Σ_{j<K} χ_0(j)` for the principal character mod N. |
| **Correctness argument or conjecture** | Reformulation `proved` and tested. |
| **Cost in n, including construction and extraction** | Known: `K` symbols, or BGS. |
| **Randomness and success probability** | — |
| **Closest prior method and the actual difference** | Class-number formulas give `Σ_{j<p/2}(j/p)` at special endpoints. Class-group methods (Schnorr–Lenstra, SQUFOF) are known factoring routes (UNVERIFIED as a link to this exact sum). |
| **Smallest informative experiment** | None proposed. The characters computable from N alone are powers of `(·/N)`, and `(·/N)²` is `χ_0` itself, so character theory adds nothing. |
| **Counterexample or scaling behavior that would refute this version** | — |
| **Possible repair if that happens** | Status: `not yet tested`. The open primitive is stated precisely. |

---

## C-M5 — Compressed discriminant (codex card C1)

| Field | Content |
|---|---|
| **ID and one-sentence mechanism** | C-M5: `disc(∏_{j<K}(T−j)) = ±∏_{k<K} k!²` is a unit iff `U(K)`. Use the doubling `D_{2m} = D_m² Res(P_m(T), P_m(T−m))²`. |
| **Target problem / supported input family** | `U(K)`. |
| **Representation and how it is built from N** | Shift grammar `P_{a+b}(T) = P_a(T)·P_b(T−a)`. |
| **Exact operation that extracts information** | Unit test of the cross-resultants. |
| **Why that information exposes a factor** | Direct. |
| **Known lemmas being used** | Resultant product formula. |
| **New lemma required** | A sub-`√m` resultant of shifted Pochhammers. |
| **Correctness argument or conjecture** | **Refuted as speedup.** `Res(P_m(T), P_m(T−m)) = ∏_{j<m} P_m(j+m)` costs `Õ(m)` by multipoint evaluation, giving `Õ(K)` overall versus BGS `Õ(√K)` for the same predicate. The required primitive is the original problem. |
| **Cost in n, including construction and extraction** | `Õ(K)`. |
| **Randomness and success probability** | — |
| **Closest prior method and the actual difference** | — |
| **Smallest informative experiment** | Analytic. |
| **Counterexample or scaling behavior that would refute this version** | — |
| **Possible repair if that happens** | None. |

---

## C-M6 — Fixed-dimension recurrences / quotient-ring powers (C-finite invariants)

| Field | Content |
|---|---|
| **ID and one-sentence mechanism** | C-M6: carry unit-only information in `K ↦ M^K v` (or `(1+X)^K mod (N, f)`, or any rank-`R` algebra power) with `O(R^ω log K)` jump-ahead. |
| **Target problem / supported input family** | `U(K)`, `I(A,B)`. |
| **Representation and how it is built from N** | `R×R` matrices over Z/NZ, which may even depend on `p` non-uniformly. |
| **Exact operation that extracts information** | Nonunit test of a linear or polynomial readout. |
| **Why that information exposes a factor** | — |
| **Known lemmas being used** | Cayley–Hamilton; zero runs of LRS. |
| **New lemma required** | — |
| **Correctness argument or conjecture** | **Refuted** (`proved`, theorem O1). A faithful readout needs order `R ≥ p+1` for the prefix predicate and `R ≥ L+1` for the interval predicate (tight). E4 exhaustive small cases: 0 sequences. Not covered: "mixed" `N = pq` readouts whose correctness relies on `q`-zeros filling gaps. |
| **Cost in n, including construction and extraction** | `poly(p)` at minimum, so no gain. |
| **Randomness and success probability** | — |
| **Closest prior method and the actual difference** | AKS-type sketches indexed by the exponent; Lucas/`p±1`. |
| **Smallest informative experiment** | E4. |
| **Counterexample or scaling behavior that would refute this version** | — |
| **Possible repair if that happens** | Leave the C-finite class: P-recursive (BGS, `√`) or nonlinear dynamics without a commutative group-law jump. |

---

## C-M7 — Power-map compositions `Q(φ(x))`, `φ` = powers and translations

| Field | Content |
|---|---|
| **ID and one-sentence mechanism** | C-M7: replace the degree-`L` Pochhammer by `Q((x+u)^e)` with small `deg Q`. Huge `e` costs only `log e`. |
| **Target problem / supported input family** | `I(A, A+L)` as a function of `A`. |
| **Representation and how it is built from N** | Exponents and constants are public. |
| **Exact operation that extracts information** | `gcd(Q((A+u)^e), N)`. |
| **Why that information exposes a factor** | It would if the zero set mod `p` equalled `Z_L`. |
| **Known lemmas being used** | Fibres of `y ↦ y^e` are `μ_g`-cosets, `g = gcd(e, p−1)`. |
| **New lemma required** | — |
| **Correctness argument or conjecture** | **One layer refuted** (`proved`, O2): `deg Q ≥ (L−1)/2`; only the `±` symmetry survives. E5: 0 violations over all primes `< 140`; the bound is attained. **k layers:** `#roots(Q) ≥ L/∏ gcd(e_i, p−1)` (`proved`, O2′). ~~So compression needs `p−1` to share large factors with uniform exponents.~~ **Corrected 2026-10-17 after r3 audit:** O2′ does not exclude `k = O(log L)` layers with bounded gcds. The bounded-gcd form `≥ G^{−k} L` is already proved. The C-2R repair cycle (`notes/repair_C2R_results.md`) shows that for two gcd-2 layers the only exponent beating `3L/8` in the tested range is the character exponent `(p+1)/2`, which encodes `p−1`. E8 (`p ≤ 31`): two layers with `gcd ≤ 2` reach `≈ 0.6–0.8` of the one-layer bound. `e ≡ 0` or `(p−1)/2` gives `O(1)` roots only through small-`p` character coincidences. A general k-layer bound beyond O2′ with a better constant is open (the bounded-gcd form with constant `G` is proved). |
| **Cost in n, including construction and extraction** | `deg Q + O(Σ log e_i)`. |
| **Randomness and success probability** | — |
| **Closest prior method and the actual difference** | `∏(x² − i²)` halving; Pollard `p−1`. |
| **Smallest informative experiment** | E5, E8. |
| **Counterexample or scaling behavior that would refute this version** | Done for one layer. |
| **Possible repair if that happens** | The next experiment is E8 at larger `p`, with sampled (not exhaustive) exponents and with additive mixing of two branches (fewnomials). That leaves the power-map class. |

---

## C-M8 — Group-law identities with many roots (`x^K − 1`, Dickson, division polynomials)

| Field | Content |
|---|---|
| **ID and one-sentence mechanism** | C-M8: use short-circuit polynomials with many roots so that `∏_{i<L}(x−i)` divides them mod `p`. Then the unit status of the evaluation reveals `p`. |
| **Target problem / supported input family** | `U(K)`. |
| **Representation and how it is built from N** | Repeated squaring, Lucas chains, EDS doubling. |
| **Exact operation that extracts information** | `gcd(f(a), N)`. |
| **Why that information exposes a factor** | — |
| **Known lemmas being used** | — |
| **New lemma required** | — |
| **Correctness argument or conjecture** | `proved` for `G_m`: if `[1, L−1] ⊂ μ_K(F_p)` and `L−1 ≥ √p+1`, then `(p−1) \| K` (least non-member of a proper subgroup is `< √p+1`; E6 confirms with max ratio 0.69 for `p < 20000`). So the identity works only when the group order divides `K`: it is Pollard `p−1`, and similarly `p+1` and ECM (those analogues are `conjectured` here). (r3 audit, 2026-10-17: for curves the relevant quantity is the point's order, not `#E`. The `√p` lemma is much weaker than Burthe / Martin–Parvardi `p^{1/(4√e)+ε}`, so it is not a new asymptotic obstruction.) |
| **Cost in n, including construction and extraction** | That of the smoothness bound, not `L`. |
| **Randomness and success probability** | — |
| **Closest prior method and the actual difference** | `p−1`, `p+1`, ECM; Cheng 2004 (`m·n!` via division polynomials; extra factors break exactness). |
| **Smallest informative experiment** | E6. |
| **Counterexample or scaling behavior that would refute this version** | — |
| **Possible repair if that happens** | None within group laws. Zero sets are cosets or torsion, not intervals. |

---

## C-HHloc — Harvey–Hittmeir restricted to a certified interval (from codex r1 idea 1)

| Field | Content |
|---|---|
| **ID and one-sentence mechanism** | C-HHloc: given `P ≥ A`, decide `P ∈ [A, B)` near `√N` by HH collision search over cells inside the interval. |
| **Target problem / supported input family** | Location of the **prime** `P` (not arbitrary multiples). Balanced semiprimes. |
| **Representation and how it is built from N** | As in HH. |
| **Exact operation that extracts information** | Collision gcds. |
| **Why that information exposes a factor** | — |
| **Known lemmas being used** | HH Prop. 3.3 (UNVERIFIED by factC). |
| **New lemma required** | — |
| **Correctness argument or conjecture** | Codex claims a "proved adaptation". `not yet tested` by factC. |
| **Cost in n, including construction and extraction** | Claimed `Õ((1 + W^{3/5} N^{−1/10}) n^3)` for width `W`, which beats BGS's `√W` on intervals near `√N`. The first wide step still costs `N^{1/5}`. **Not polynomial.** |
| **Randomness and success probability** | Las Vegas order certification. |
| **Closest prior method and the actual difference** | HH 2021. |
| **Smallest informative experiment** | Implement HH cell search on 40–56-bit semiprimes and compare with BGS on `I(A,B)` for `A ≈ √N`. |
| **Counterexample or scaling behavior that would refute this version** | — |
| **Possible repair if that happens** | Out of track C's unit-only scope. Flagged for the orchestrator. |
