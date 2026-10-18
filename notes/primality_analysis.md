# Primality tests and the factors of N: exact leak analysis (agent `primality`, 2026-10-18)

**Question (user).** "How come there is no way to get insights into which primes divide N from a primality test, if N is
composite?" Follow-ups: can the failure residue be iterated or recursed on, and are there more favourable exponents or groups
than `a^(N−1)`, and can such checks be combined?

**Short answer.**
1. A test computes a value `z` in a ring `R_N ≅ R_p × R_q` (CRT). A gcd reveals a factor **iff `z` vanishes in exactly one
   component** (Thm 0).
2. For every power-type test (Fermat, Miller–Rabin, Solovay–Strassen, Lucas, Frobenius), "vanishes at `p`" is the
   group-order coincidence `ord_p(x) | E` for the public exponent `E`. For `E` built from `N` this has probability
   `≈ g/p` with `g = gcd(p−1, q−1)`: about `2^{−n/2}·polylog` on balanced semiprimes (Thms 1–3; measured, E1).
3. **The output does contain the factors**, as a short discrete logarithm. The Fermat residue is
   `a^{N−1} = a^{p+q−2}`, and Lucas residues are `z^{±(p+q)}` or `z^{±(q−p)+2}`. Reading the factors off costs an interval
   discrete log, `~N^{1/4}` group operations (Thm 5; measured slope 0.24, E2). This is Pollard's 1978 method; it is the
   baseline, not a new algorithm.
4. Partial information (2-adic pattern, Euler values, the sum-of-two-squares bit `ε`) leaks only through the same rare
   liar events. A transcript certifies `ε = +1` (by exhibiting `√−1 mod N`) with probability `≈ g²/(3N)` per base; the
   continued-fraction control certifies it far more often at fixed `T` (Thm 8, Prop 9; E4).
5. Logically this is forced: Fermat/MR witnesses can be sampled from `N` alone. So any efficient map "witness → factor"
   **is** an efficient factoring algorithm (Thm 10).

Status labels follow AGENT_RULES rule 5. GPT-6.1-sol's literature note (`notes/codex_r6_primality_lit.md`) was checked
claim by claim (§9). Every formula it states that we could test agrees exhaustively with our independent derivation.

---

## 0. Setting and the CRT leak criterion

`N` is odd and squarefree, `N = ∏ p_i`; mostly `N = pq` with `p < q`. A base is uniform in the relevant group
`G(Z/N) = ∏ G(F_{p_i})` (CRT, independent components).

**Theorem 0 (CRT criterion; proved, elementary).** Let `R_N = Z/N[X]/(f)` (`f` monic, possibly `f = X`) and `z ∈ R_N` with
coefficient vector `(z_0,…,z_{d−1})`. Then `gcd(z_0,…,z_{d−1}, N)` is a proper divisor `⇔` `z ≡ 0` in some but not all of
the components `R_N/p_i R_N`. A *single* coefficient gcd `gcd(z_j, N)` can also be proper while `z` vanishes in no
component.
* The first kind is called **structural**.
* The second kind is **accidental**; it has probability about `1/p` per coordinate when the coordinates behave randomly.
  §3 (AKS) shows a regime where they do not.

*Proof.* `p_i | gcd(z_0,…,N) ⇔ p_i | z_j ∀j ⇔ z ≡ 0 mod p_i`. ∎ Every test in `src/primality/ops.py` reports both kinds
separately.

So **"the test rejects" means `z ≢ 0 mod N`**, i.e. `z ≢ 0` modulo *at least one* prime. A factor needs the asymmetric
event: `≡ 0` modulo one prime, `≢ 0` modulo another. For a random base, rejection is almost always **symmetric**: `z ≢ 0`
at both primes. Exactly (Fermat, `N = pq`):

`P(fails at both | witness) = (1−A)(1−B)/(1−AB)`, with `A = g/(p−1)`, `B = g/(q−1)`.

## 1. Order-pattern lemma and exact leak formulas

**Theorem 1 (order-pattern lemma; proved).** Let `x` be uniform in `G = ∏_i C_{n_i}` (cyclic components). A check `(E, ζ)`,
`ζ = ±1`, has gcd observable `gcd(x^E − ζ, N) = ∏{p_i : x_i^E = ζ}`. The indicator `[x_i^E = ζ]` depends only on
`d_i = ord(x_i)`:
* `[d_i | E]` for `ζ = 1`;
* `[d_i | 2E, d_i ∤ E]` for `ζ = −1`, since `−1` is the unique element of order 2.

The `d_i` are independent with `P(d_i = d) = φ(d)/n_i` (or `φ(d)/(n_i − 1)` for `d > 1` when the identity is excluded). A
family of checks reveals a proper divisor iff the primes' indicator vectors are not all equal:

`P(leak) = 1 − Σ_v ∏_i P_i(v)`.

*Proof.* A cyclic group has exactly `φ(d)` elements of order `d`, and `x^E = 1 ⇔ d | E`; CRT gives independence. ∎
Implemented as `formulas.indicator_dist / leak_prob`.

**Theorem 2 (MR ladder law; proved).** Let `E = 2^t m` (m odd) and let `τ_i` be the least `j` with `x_i^{m2^j} = 1`. If the
odd part of `ord(x_i)` does not divide `m`, set `τ_i = ∞`. Then `τ_i` is independent across `i`, with

`P(τ_i = j) = (h_i/m_i)·w_{s_i}(j)`,

where `n_i = 2^{s_i} m_i`, `h_i = gcd(m, m_i)`, `w_s(0) = 2^{−s}` and `w_s(j) = 2^{j−1−s}` (`1 ≤ j ≤ s`). The ladder values
`b_j = x^{m2^j}` satisfy `b_j ≡ 1 (mod p_i) ⇔ τ_i ≤ j` and `b_j ≡ −1 (mod p_i) ⇔ τ_i = j+1`. Hence:
* **sqrt** (classical "nontrivial √1"): `max τ ≤ t` and the `τ` are not all equal;
* **single** (Fermat-type split): some but not all `τ_i ≤ t`;
* **prefix** (some `gcd(b_j − 1, N)` proper, `j ≤ t`) = sqrt ⊔ single;
* **pm** (also `gcd(b_j + 1, N)`, including `b_t + 1`): prefix plus "one `τ_i = t+1`, the others `> t+1`";
* **strong liar**: all `τ_i` equal and `≤ t`.

`formulas.tau_dist`, `ladder_events`.

**Theorem 3 (closed forms, two cyclic components; proved).** For orders `n_p, n_q` and exponent `e = 2^t m`, put
* `D_i = gcd(e, n_i)`, `κ_i = D_i/n_i`;
* `u_i = min(v₂ n_i, t)`, `u = min(u_p, u_q)`.

Then:

| event | exact probability |
|---|---|
| single | `κ_p(1−κ_q) + κ_q(1−κ_p)` |
| sqrt | `κ_pκ_q [1 − (4^u + 2)/(3·2^{u_p+u_q})]` |
| prefix | single + sqrt |
| strong liar | `odd(D_p)·odd(D_q)·(4^u + 2) / (3 n_p n_q)` |

*Proof.* Use `κ_i = (h_i/m_i)·2^{u_i−s_i}` and `Σ_{j≤u} w_{s_p}(j)w_{s_q}(j) = 2^{−s_p−s_q}(4^u+2)/3`. ∎

**Corollary 3a (Fermat / MR on `N = pq`; proved).**
* `N − 1 ≡ q − 1 (mod p−1)`, so `D_p = D_q = g := gcd(p−1, q−1)`.
* `t = v₂(N−1)` equals `k := min(v₂(p−1), v₂(q−1))` iff `v₂(p−1) ≠ v₂(q−1)`, and `t > k` otherwise.

With `A = g/(p−1)`, `B = g/(q−1)`, `F = AB = g²/φ(N)`, `h = odd(g)`:

| event | exact probability per unit base |
|---|---|
| Fermat `gcd(a^{N−1}−1, N)` splits | `A + B − 2F` |
| MR nontrivial `√1` | `(2/3)(1 − 4^{−k})·F` |
| MR prefix (= GPT's "augmented MR") | `A + B − F(1 + C_k)`, `C_k = (4^k+2)/(3·4^k)` |
| strong liar (Monier) | `h²(4^k+2)/(3φ(N)) = F·C_k` |
| MR exhibits `√−1` (`b_i = −1` for some `1 ≤ i < t`) | `h²(4^k − 4)/(3φ(N)) = F(4^k−4)/(3·4^k)` |

**Corollary 3b (Carmichael; proved).** If `λ(N) | N−1`, every `A_i = 1` and the Fermat leak is 0. The MR prefix factor
probability is `1 − Σ_{j≤min s_i} ∏_i w_{s_i}(j) ≥ 1 − 2^{1−r}` (`r` = number of primes).
* `561 = 3·11·17`: exactly `31/32`. Base 2 gives the ladder `263, 166, 67, 1, 1`, and `gcd(67 − 1, 561) = 33`.
* Chernick numbers at 40–64 bits: the minimum exact rate is `7/8`.

**Corollary 3c (Solovay–Strassen; proved).** `b = a^{(N−1)/2}` is `b_{t−1}` on the MR ladder. With
`d' = gcd((N−1)/2, p_i−1)`:

`P(gcd(b ∓ 1, N) proper) = 1 − [d'_p d'_q + (g−d'_p)(g−d'_q)]/φ − (1−A)(1−B)`.

This equals GPT's `A + B − 2F` (`u = v`) and `A + B − (3/2)F` (`u ≠ v`). The "residual" gcd `gcd(b − (a/N), N)` alone is
`(A+B)/2 − F`.

**Corollary 3d (Lucas / strong Lucas; proved).** Take `z = α/β` (norm one) in `T_N = T_p × T_q`, with `|T_r| = r − (D/r)`.
* `U_k ≡ 0 (mod r) ⇔ z^k = 1` and `V_k ≡ 0 ⇔ z^k = −1` (since `α − β = √D` is a unit). So Thm 3 applies with
  `n_r = r − (D/r)` and `e = N + 1`.
* For `N = pq`: `gcd(N+1, p − δ_p) = gcd(q + δ_p, p − δ_p)`.
* Uniform `P` (with `Q = (P²−D)/4` a unit) corresponds to `z` uniform on `(T_p∖{1}) × (T_q∖{1})`. The identity-excluded
  law of Thm 1 applies. GPT's formula (5) is the same statement.

**Quadratic Frobenius (proved).** Let `α` be uniform in `(Z/N[√Δ])^*` and `z = α^N − σ(α)`. The local kernel sizes are
* inert `p`: `gcd(q−1, p²−1)` out of `p²−1`;
* split `p`: `gcd(q²−1, p−1)` out of `(p−1)²`.

So the structural leak is `κ_p(1−κ_q) + κ_q(1−κ_p)` with these ratios.

**Lemma 4 (proved).** For odd `N > 1` and any unit `a`, `a^{N−1} ≢ −1 (mod N)`.
*Proof.* It would force `v₂(ord_{p}(a)) = t+1` for every `p | N`, hence `2^{t+1} | p−1` for all `p | N`, hence
`N ≡ 1 (mod 2^{t+1})`. That contradicts `t = v₂(N−1)` (prime powers: `ord | p^{k−1}(p−1)`, same 2-part). ∎

**Verification (all exhaustive over every base; 0 failures).** `src/primality/verify.py`,
`experiments/raw/primality_verify.json`.
* MR (prefix, sqrt, pm, liar, liar levels), Fermat (final and full binary chain), SS and Lemma 4:
  * all 3 819 odd semiprimes `< 20 000`;
  * 1 796 three-prime `N < 20 000`;
  * the 22 Carmichael numbers `≤ 200 000`.
* Strong Lucas (Selfridge `D`, every `P`; prefix, pm, liar and the plain `U_{N+1}` test): 1 623 semiprimes `< 12 000`.
* Frobenius (every `α ∈ R^*`): 76 semiprimes `< 700`.
* Closed forms (Cor 3a, 3c, Thm 3) agree with the engine on every instance.

## 2. Measured leak rates (E1, `experiments/analysis/primality_tables.md`)

Balanced semiprimes (`balanced_semiprime`), random unit bases, exact per-instance predictions summed:

| bits | instances × bases | Fermat exp / obs | MR prefix exp / obs | MR sqrt exp / obs | median log₂ rate |
|---:|---:|---|---|---|---:|
| 32 | 20 000 × 5 | 54.8 / 57 | 55.1 / 58 | 0.28 / 1 | −12.9 |
| 48 | 5 000 × 5 | 0.068 / 0 | 0.068 / 0 | 4e−6 / 0 | −20.9 |
| 64 | 5 000 × 5 | 2.5e−4 / 0 | 2.5e−4 / 0 | 7e−11 / 0 | −28.9 |
| 96 | 2 000 × 3 | 7e−10 / 0 | 7e−10 / 0 | — | −44.9 |
| 128 | 2 000 × 3 | 1e−14 / 0 | 1e−14 / 0 | — | −61.6 |

* The median rate is `≈ 2^{−n/2 + 3}`, i.e. `~2g/√N` (mean `g ≈ 15–20` at every size). **Experimentally supported.**
* At 32 bits:
  * MR with all `±1` gcds: 73.7 expected, 78 observed.
  * Strong Lucas (P-sampled, U and V): 85.3 expected, 87 observed.
  * BPSW leaked on 21 of 20 000 N. A heuristic expectation that treats base 2 and `P = 1` as random gives 31.8, so the
    random-base model is only approximate for the fixed BPSW parameters.
* **Diagnostic family `diag_shared_gcd`** (`p−1`, `q−1` share a prime `G ≈ p/2^{10}`), where leaks are measurable at large
  sizes:

  | bits | MR prefix exp / obs | MR sqrt exp / obs |
  |---:|---|---|
  | 64 | 584 / 621 | 3.2 / 3 |
  | 96 | 487 / 508 | 17.1 / 23 |
  | 128 | 370 / 365 | 3.4 / 4 |

  Leaks are classified STRUCTURAL: the mechanism is `ord_p | N−1` with large `g`. This family is a deliberate
  `p−1`/`q−1` coincidence and is not the intended distribution.
* **Chernick Carmichael numbers** (40/52/64 bits): 3 618–3 624 expected vs 3 633–3 654 observed, out of 4 000 bases.
* Null comparison (rule 3): nonunit bases gave 4 hits in 100 000 draws at 32 bits (ACCIDENTAL; null
  `≈ 2/√N·10^5 ≈ 3`). No other accidental hits arise in scalar tests.

**Tie-in with Track B/D.** Every leak above is a group-order coincidence. `ord_p(x)` must divide a public exponent `E`;
since `ord_p(x) | |G(F_p)|`, it must divide `gcd(E, |G(F_p)|)`. This is exactly Track B/D Prop D2: `G_1(α) = 0 ⇔ β^{N−1} = 1`
has `g − 1` roots. It is also B4: vanishing depends on `q` only modulo the local order. AKS is the polynomial instance (§4).

## 3. The discrete-log view: the residue contains `s = p + q`

**Theorem 5 (proved).** For every unit `a`: `a^{N−1} = a^{s−2}`, `a^N = a^{s−1}`, `a^{N+1} = a^s (mod N)`.
*Proof.* `N − 1 = φ(N) + (s − 2)`. ∎

For the norm-one torus `T_N` of `Z/N[√D]`, with `δ_r = (D/r)`, `|T_r| = r − δ_r`:

`(p−δ_p)(q−δ_q) = N − δ_q p − δ_p q + δ_pδ_q` annihilates `T_N`, so **`z^N = z^{δ_q p + δ_p q − δ_pδ_q}`**, i.e.

| `(δ_p, δ_q)` | `z^{N+1}` | Lucas-test exponent `z^{N−(D/N)}` |
|---|---|---|
| `(+1,+1)` | `z^{s}` | `z^{N−1} = z^{s−2}` |
| `(−1,−1)` | `z^{−s}` | `z^{N−1} = z^{−s−2}` |
| `(+1,−1)` | `z^{q−p+2}` | `z^{N+1} = z^{q−p+2}` |
| `(−1,+1)` | `z^{p−q+2}` | `z^{N+1} = z^{p−q+2}` |

* The standard Lucas test needs `(D/N) = −1`, so for `N = pq` it is **always a mixed row**: it encodes `q − p`, not
  `p + q`.
* For a general (non-norm-one) `α`:
  * mod an inert `p`: `α^N ≡ σ(α)^q`;
  * mod a split `p`: `α^N ≡ α^q`;
  * if both primes are inert: `α^{N+1} = σ(α)^{s}` (σ is the public conjugation).

  So the Frobenius-test residue is a pure power of the *public* element `σ(α)` in that case.
* GPT's normalisation: for `ρ` of norm `Q`, `ρ^{N+1} = Q^K ρ^{(table)}` with `K = 0, s, p−1, q−1` in the four rows.
* **Proved; verified** on 4 696 (prime pair, `D`) cases covering all four patterns equally, 20 random elements each, plus
  1 248 cases for the `Q^K` form (0 failures).

**Demonstration (E2, operational `src/primality/dlog.py`).** The code
1. takes only `N`;
2. computes the Fermat residue `y = a^{N−1}`;
3. solves `a^x = y` for `x ∈ [2√N − 2, 3√N − 2]` (kangaroo or BSGS; uniqueness holds when `ord(a)` exceeds the width);
4. factors from `s = x + 2` via `s² − 4N = (q−p)²`.

The Lucas variant solves `z^c = z^{N+1}` on `|c − 2| ≤ √(5N)` and factors from `q − p`.

| method | sizes | success | fit of log₂(mean mults) vs n (n ≥ 24) |
|---|---|---|---|
| Fermat residue, kangaroo | 16–68 bits, 30 each | 420/420 | slope **0.241** |
| Fermat residue, BSGS | 16–60 bits | 360/360 | slope **0.236** (≈ `N^{1/4}`; 30 744 mults at 60 bits vs `N^{1/4}` = 30 488) |
| Lucas residue, kangaroo | 16–56 bits | 180/180 | slope **0.229** |

Example: `N = 10403 = 101·103`, `a = 2`. The test outputs `y = 2^{10402} mod 10403 = 9296`, and `gcd(9295, N) = 1`. The
kangaroo finds `x = 202` (65 multiplications), so `s = 204`, `s² − 4N = 4`, giving `101·103`.

**Conclusion (proved + measured).** The rejected test's output determines the factors (for `ord(a) > √N`). It stores them
as a discrete logarithm in an interval of width `Θ(√N)`. Generic interval-DLP algorithms (BSGS deterministic; kangaroo
heuristic) need `Θ(√W) = Θ(N^{1/4})` operations. Exponential in `n`; **not polynomial** (rule 8).

**Prior art (from `notes/codex_r6_primality_lit.md`, citations checked against its quoted statements, not re-read).**
* Pollard 1978 (*Math. Comp.* 32, §3) proposed exactly this: kangaroo on `a^{2x} = a^{N+1}`, `x = (p+q)/2`.
* Shanks's BSGS is the deterministic search.
* McKee–Pinch 1996 derandomise it.
* Hittmeir 2018 sieves `p+q` residues and reaches `N^{1/4}·exp(−Ω(log N/log log N))`, deterministic.
* Lehman-type weighted sums lead to Hittmeir 2021 (`N^{2/9}`) and Harvey 2021 (`N^{1/5}`, deterministic).

Nothing in this section is new.

**AKS (Thm 7; proved, verified).** AKS checks `Δ_a = (X+a)^N − X^N − a ≡ 0` in `R_N[X]/(X^r − 1)`. Modulo `p`, with
`τ_p : X ↦ X^{p mod r}` (= Frobenius on `F_p[X]/(X^r−1)`):

`Δ_a(N) ≡ τ_p( (X+a)^q − X^q − a )  (mod p)`.

So the mod-`p` residue is the AKS defect **of exponent `q` in characteristic `p`**, twisted by an automorphism that
depends on the unknown `p mod r`. It is not a power of one public element: modulo `p` the residue is
`τ_p((X+a)^q) − …` and modulo `q` it is `τ_q((X+a)^p) − …`, two different twists.
* Its mod-`p` content depends on `q` only modulo `L = lcm(r, ord(X+a))`, where `L | p^{ord_r p} − 1`. This is Track B/D's
  `L_α`.
* `Δ_a ≡ 0 (mod p)` iff `(X+a)^q = X^q + a` in every component. A Frobenius coincidence `q ≡ p^j (mod L)` suffices
  (B4).
* Taking the norm `Z/N[X]/(Φ_r) → Z/N` (a product of `r−1` public automorphic images) collapses the twist. It returns the
  Fermat residue `c_a^N = c_a^{s−1}` of the integer `c_a = ±((−a)^r − 1)`. So AKS's output contains at least the Fermat
  discrete log.

Verified: twist identity on 300 random `(p, q, r, a)` cases, and period `L` on all cases where `ord_r p ≤ 4` (0 failures).
* **AKS coefficient leak, a Kummer regime (new to this note; mechanism = Track A L6–L7 support).** If `q < 2p`, then mod
  `p`, `(X+a)^q = (X^p + a)(X+a)^{q−p}`. Its folded support is the union of two arcs of length `q−p+1` in `Z/r`. If the
  arcs do not cover `Z/r`, then **`r − cov` coefficients of `Δ_a` vanish mod `p` but not mod `q`**, and their gcds return
  `p`.

  | bits | `r` | N in Kummer regime | hits | predicted `64(r − cov)` | other N: hits / null `64r(1/p+1/q)` |
  |---:|---:|---:|---:|---:|---|
  | 20 | 397 | 40/40 | 689 027 | 686 720 | — |
  | 24 | 563 | 27/40 | 386 822 | 385 536 | 273 / 270 |
  | 28 | 773 | 18/40 | 309 689 | 309 120 | 174 / 152 |
  | 32 | 1013 | 4/40 | 144 193 | 144 064 | 96 / 83 |
  | 40 | 1579 | 1/40 | 95 489 | 95 424 | 12 / 8.7 |

  The prediction comes from the first factor only; the small surplus is accidental zeros. **Experimentally supported;
  mechanism proved.**

  The regime requires `q − p ≲ r = Õ(log⁵ N)`. That is exactly where Fermat's difference-of-squares method factors
  instantly; it vanishes on balanced semiprimes as `n` grows (1/40 at 40 bits). Whole-polynomial (structural) leaks: 0
  everywhere. AKS's trial-division step finds nothing for balanced `N ≥ 2^{20}` (`p > r`).

## 3b. Iterating / recursing on the failure residue (orchestrator follow-up (3))

**Proposition 6 (proved).** Group operations and powering by public integers, applied to the public `a` and to
`y = a^{N−1}`, only produce elements `a^{E}` with `E` a *public* integer.
* Since `y = a^{N−1}`, any word in `a, y` is `a^{c₀ + c₁(N−1) + c₂(N−1)² + …}`.
* Reduced modulo `ord(a)`, the exponent is `f(s−2)` for the known polynomial `f(X) = Σ c_i X^i`. For example
  `y^{N−1} = a^{(s−2)²}`, and `y^{N+1}·a^{−2} = a^{(s−2)s − 2}`.
* Combining bases (`a^{s−2}·b^{s−2} = (ab)^{s−2}`) gives the same unknown `s` in another base.
* Mixing Fermat with Lucas gives `s` and `q−p`. These are not independent (`s² − (q−p)² = 4N`).

So every derived residue is an interval-DLP-type equation `g^{f(s)} = Y` with one hidden integer `s` ranging over
`Θ(√N)` values. Two consequences:
1. **Degree ≥ 2 makes it worse.**
   * BSGS and kangaroo need `f` linear: a step `s → s+1` must multiply by a fixed element.
   * For `deg f = 2` the exponent interval has width `Θ(N) ≳ ord(a)`. The DLP solution stops being unique, and any
     collision `a^x = a^{(N−1)²}` gives a multiple of `ord(a)`. That is **order finding**, `~√ord(a) ≈ N^{1/2}`
     generically.
   * Measured (E3, 16–36 bits, 20 each):
     * squared-residue kangaroo: median mults grow with slope **0.52**, success 8–15/20 via the order-multiple ladder
       split, and **0** runs recovered `s` itself;
     * plain residue: slope 0.25;
     * a direct scan over `s` with finite differences: linear in `s − 2√N`, i.e. Fermat's method.
2. **Combining residues** gives no reduction of the search space. The pair `(a^{s−2}, z^{q−p+2})` is determined by `s`.
   * Measured BSGS costs: the Fermat residue alone (width `√N`) is cheaper than the Lucas residue (width `2√(5N)`) by
     ≈ 3× at every size.
   * No joint algorithm below `√W` was found or is known to us. Since the second equation is a deterministic function of
     the first unknown, we **conjecture** that generic joint algorithms need `Ω(N^{1/4})`. Shoup-type generic lower
     bounds are for prime-order groups and single DLPs, so this is not covered by a theorem we know.

The known way to make "combining" pay is to use collisions that are **one-sided** (Lehman/Hittmeir/Harvey weighted sums,
§5b (d)).

## 4. Partial information: 2-adic patterns, Euler values, and `ε`

**Theorem 8 (transcript law; proved).**
* For `N = pq`, the ±1-pattern of the MR ladder (`S0`, `L1`, …, `L_t`, `NONE`) and the class of the Euler value
  (`+1`, `−1`, other) have exact laws given by Thm 2.
* The pattern `NONE` (no `±1` anywhere) has probability `≥ 1 − (A + B)`, with `A + B = g/(p−1) + g/(q−1)`.
* All other patterns, every certificate, and every leak live on this complementary event.

Hence, for any two semiprimes, or any two hidden-parameter classes (e.g. `ε = ±1` given `N mod 8`), the total-variation
distance between the per-base pattern laws is at most `E[A + B | class 1] + E[A + B | class 2]`. With `k` bases it is at
most `k` times that. On balanced inputs this is `≈ 2^{−n/2+4}`.

* GPT's (8), (10) and (11) give the marginal laws: `P(x_i = ±1)`, `P(√−1 exhibited)`, `P(b = −1)`. We verified them
  exhaustively (§9).
* **The pattern does carry the hidden data**: liar levels `L_j` occur iff `j ≤ min(v₂(p−1), v₂(q−1))`. In particular
  `L_j` with `j ≥ 2` occurs only if `ε = +1` (for `N ≡ 1 mod 4`). But all of it is on the rare liar event.

**Exact conditional laws (E4b; all prime pairs in `[2^{b−1}, 2^b)`, `N ≡ 1 mod 4`, exact rationals).**
* `TV(MR pattern | ε=+1 vs −1, N mod 8)`: `6·10^{−3}` (b = 7) down to `3·10^{−4}` (b = 11).
* SS: similar.
* `P(L_{≥2} | ε = −1) = 0` exactly; `P(L_{≥2} | ε = +1) ≈ TV`.
* The leak probability itself depends on `ε`, e.g. 0.029 vs 0.007 at b = 11. The reason is that `ε = +1` forces
  `4 | g`. It is still `O(g/p)`.

**One-sided `ε` certificates vs `n` (E4a, 4 000 instances per size, `N ≡ 1 mod 4`).** MR (`b_i = −1`, `i ≥ 1`) and SS
(`b = −1`) both exhibit `x` with `x² ≡ −1 (mod N)`. That proves `ε = +1`, and Hermite–Serret / Brillhart gives
`N = u² + v²` (`dlog.sos_from_sqrt_minus1`, tested).

| bits | MR rate/base (mean, ε=+1) | SS rate/base | CF PC1 hits (T=1024) | CF PC2 hits (T=1024) | per ~1 000 ε=+1 N |
|---:|---:|---:|---:|---:|---|
| 16 | 2.0e−3 | 1.5e−3 | 389 | 627 | 748 |
| 32 | 6.1e−7 | 2.4e−7 | 29 | 977 | 990 |
| 48 | 8.7e−11 | 6.0e−11 | 1 | 252 | 1 013 |
| 64 | 5.7e−15 | 3.7e−15 | 0 | 19 | 1 009 |
| 128 | 3.3e−35 | 2.3e−35 | 0 | 0 | 960 |

* PC1: odd CF period of `√N` within `T`. PC2: odd `i ≤ T` with `Q_i` a square. Both use our own recurrence,
  `experiments.cf_controls`.
* The MR/SS rate is `≈ F/3 ≈ g²/(3N)`: about `2^{−n}`, so useless.
* The continued-fraction control fires with rate `~T·N^{−1/4}` (PC2 is SQUFOF-like). So it is the far better
  `ε`-certificate per unit work; it is still `N^{1/4}`-type, the same cost as factoring by SQUFOF.
* No certificate fired on an `ε = −1` instance (consistent with proofs).

**Pre-registered battery (E4c; config `experiments/configs/primality_battery.json` committed in `be7723d` before the run).**
* **Design.** 14 transcript features from MR/SS/Lucas: ±1 patterns, certificate events, low bits of residues, Jacobi
  symbols of shifted residues. 40 000 seeds per size at 32/48/64 bits (≈ 20 000 usable, `N ≡ 1 mod 4`), 4 bases each.
  Within-stratum (`N mod 8`) Welch z. Bonferroni threshold `|z| ≥ 3.32`; MDE `|ρ| ≈ 0.029` at 80% power.
* **Result.** **0 hits**; max `|z| = 2.16`.
* The ±1-pattern and certificate features were **DEGENERATE** at all sizes (never fired), as predicted from Thm 8.
* Controls detected:

  | control | z (32 / 48 / 64 bits) |
  |---|---|
  | natural CF-PC2 (T = 256) | +154 / +27.5 / +7.2 |
  | planted ρ = 0.06 | +8.0 / +8.0 / +8.2 |

* **Verdict:** null for this menu at these sizes. It excludes correlations ≥ 0.029; it does not show `ε` is unpredictable.

**Proposition 9 (2-adic extraction needs roots of unity; proved).** Let `N = pq ≡ 1 (mod 4)`. Reading `s mod 2^k` (k ≥ 3)
out of `y = a^{s−2}` by Pohlig–Hellman needs the 2-Sylow projection of `y` or a known element of 2-power order.
(i) **Projection is factoring-hard.** Projecting `y` and `a` to the 2-Sylow subgroup needs an odd multiple `M` of
    `odd(λ(N))`. Then `M·2^n` is a multiple of `λ(N)`, which factors `N` in randomized polynomial time (Miller/Long;
    Thm 3 with `κ = 1` gives success `≥ 1/2` per base).
(ii) **Known 2-power elements are `±1`, `√−1`, or factors.** Let `z` have order `2^j`, `j ≥ 2`. Squaring `j−2` times gives
    `w` of order 4. Either `w² = −1`, so `w` is a square root of −1: it exists iff `ε = +1`, and it yields
    `N = u² + v²` in polynomial time. Or `w² ≠ −1` is a nontrivial `√1`, and `gcd(w² − 1, N)` is a proper divisor.
(iii) **Corollary.** Without factoring, the only extractor of 2-adic information is a square root of −1. Its existence *is*
    `ε = +1`, and finding one is the sum-of-two-squares search problem (REPORT §5). The MR/SS certificates above are this
    route, at rate `≈ g²/(3N)`.

*Proof.* (i) is Long 1981 / Miller 1976 as cited by GPT. (ii): an element of order 2 is `−1` or a nontrivial
`e ≡ (1, −1)`; for the latter `gcd(e−1, N) = p` or `q`. ∎

## 5. The meta-argument

**Theorem 10 (proved, elementary; the formulation is ours, the ingredients are standard).**
(a) **Recursion.** Suppose a randomized algorithm `T` runs in time `poly(n)` and, on every composite input, outputs with
    probability `≥ 1/poly(n)` a string `h` such that a `poly(n)`-time `D(N, h)` is a proper divisor. Then complete
    factorisation runs in expected polynomial time. Recurse on `d` and `N/d`, stopping at primes certified by AKS (or MR
    plus verification). The tree has at most `2 log₂ N` nodes, and divisors are checked by one division.
(b) **Samplable witnesses are useless for factoring unless factoring is easy.** MR witnesses have density `≥ 3/4` in
    `(Z/N)^*` for every odd composite `N` (Rabin–Monier). For any PPT extractor `X`:

    `Pr_{w ← MR-witnesses(N)}[X(N, w) ∈ {p, q}] ≤ (4/3)·Pr[X(N, S(N)) ∈ {p, q}]`,

    where `S` samples a uniform unit and keeps it if it is a witness. The right side is a PPT factoring success
    probability. The same holds for Fermat witnesses on non-Carmichael `N` (density `≥ 1/2`), for SS, and for Lucas or
    Frobenius parameters, with their densities. It holds for every transcript value computed from `(N, w)`.
(c) **What does *not* follow.** "Certificates of compositeness cannot encode factors" is false as a slogan: `p` itself is a
    short certificate. The correct statement concerns **efficiently producible** certificates. They are produced from `N`
    by public randomness, so they carry no more polynomial-time-accessible information about the factors than `N` does
    (GPT §8, same point).
(d) **CRT reading.** A Fermat witness certifies `∃ r | N : a^{N−1} ≢ 1 (mod r)`, an existential statement checked through
    the CRT product without naming `r`. For a random witness it holds for **both** primes with probability
    `(1−A)(1−B)/(1−AB) = 1 − O(g/p)`, so nothing in the witness singles out one prime. The test needs only the
    "OR"; a factor needs the "XOR".
(e) **Partial hints.**
    * If such a `T` produced `≥ n/4 + O(1)` low bits of `s` (Track E, `s mod 2^k`), Coppersmith/Howgrave-Graham decoding
      factors.
    * So any primality-type computation leaking a growing prefix of `s` would be a factoring algorithm.
    * The project's first such bit is `ε` (E2′): the question of REPORT §5, untouched by these tests (§4).

## 5b. Exponent and group menu (orchestrator follow-up: "more favourable checks? combined?")

**Lemma M1 (proved; checked by enumeration).** For `G(F_r) ≅ C_{d_1} × … × C_{d_k}` and uniform `x`,

`A_r(E) := P(x^E = 1) = ∏_i gcd(E, d_i) / |G(F_r)|`.

The gcd leak probability for `N = pq` is `A_p + A_q − 2A_pA_q`. This handles non-cyclic groups (split quadratic units
`C_{p−1}²`, elliptic curves `C_{d_1} × C_{d_2}`). The elliptic formula was checked by enumerating all points on random curves
mod every prime `< 200`: 1 870 (curve, `E`) pairs, 0 failures (`menu.check_ec_formula`).

**Lemma M2 (N-derived exponents are q-driven; proved).** Write `|G(F_p)| = Ψ(p)` with `p` satisfying `p^k ≡ 1` (or
`p ≡ ±1`) modulo `Ψ(p)`. Then `E(N) = E(pq)` reduces to a polynomial in `q`:
* `gcd(N^k − 1, p^k − 1) = gcd(q^k − 1, p^k − 1)`;
* `gcd(N^k − 1, Φ_k(p)) = gcd(q^k − 1, Φ_k(p))`;
* `gcd(N+1, p+1) = gcd(q−1, p+1)`;
* `gcd(N−1, p+1) = gcd(q+1, p+1)`;
* `gcd(N, p−1) = gcd(q, p−1)`;
* `gcd(N²−1, p−1) = gcd(q²−1, p−1)`.

Checked for all `p < 120`, `pq < 2000`, `k ∈ {1,2,3,4,6}` (0 failures).

So `A_p = gcd(f(q), Ψ(p))/Ψ(p)`, a "coincidence between two independent primes". Using `gcd(x, n) = Σ_{d | gcd} φ(d)` and
equidistribution of `q` mod `d`, `E_q[gcd(f(q), n)] = Σ_{d|n} ρ*_f(d)`, where `ρ*_f(d)` is the number of unit roots of `f`
mod `d`.
* For `f = X − 1` this is `τ(n)`: the mean `g` is 15.4 vs mean `τ(p−1)` = 21 at 32 bits. **Heuristic** at fixed size,
  because `q` is not equidistributed modulo divisors `d ≈ q`.
* So **every fixed-degree `N`-derived exponent has `A_p = p^{−1+o(1)}` on average**.

By contrast `E = L_B = lcm(1..B)` gives `A_p = (B-powersmooth part of Ψ(p))/Ψ(p)`. This equals 1 with probability ≈ the
Dickman-type smoothness probability of `Ψ(p)`, a function of `p` alone. That is Pollard `p−1`, Williams `p+1`,
Bach–Shallit `Φ_k`, and ECM. N-derived exponents encode the *primality hypothesis*: they test `ord | N − 1`, which holds
for primes by Fermat. Smooth exponents encode a *smoothness hypothesis about the hidden group*.

**Measured (E6; exact `A_p, A_q`; 400/300/300 balanced instances at 32/48/64 bits; full table in `primality_tables.md` §E6).**

*N-derived exponents: median log₂ of the leak probability (32 / 48 / 64 bits).* `P(leak ≥ ½) = 0` in every cell.

| group | `N−1` | `N+1` | `N` | `N²−1` | `N⁶−1` |
|---|---|---|---|---|---|
| `(Z/N)^*` | −13.6 / −20.8 / −28.9 | −12.7 / −20.6 / −28.8 | −14.8 / −22.8 / −30.8 | −11.4 / −19.2 / −27.5 | −10.3 / −18.2 / −26.4 |
| torus `T_D` (Selfridge D) | −12.5 / −20.7 / −28.6 | −12.9 / −20.9 / −29.7 | −14.8 / −22.8 / −30.8 | −11.1 / −19.3 / −27.4 | −10.0 / −18.3 / −26.1 |
| `(Z/N[√D])^*` | −26.7 / −42.5 / −58.9 | −27.1 / −43.0 / −59.0 | — | — | — |
| `F_{r^6}^*` (ideal) | −89.4 / −137.9 / −186.0 | −90.2 / −137.8 / −186.1 | — | −85.9 / −133.9 / −182.3 | −81.5 / −129.3 / −178.1 |

*Smooth exponents: `P(leak ≥ ½)` and, in brackets, `P(A_p = 1)` (32 / 48 / 64 bits).*

| group | `L_1000` | `L_100000` | `(N−1)·L_1000` | `(N²−1)(N⁶−1)·L_10000` |
|---|---|---|---|---|
| `(Z/N)^*` | 0.48 / 0.34 / 0.07 [0.61 / 0.22 / 0.04] | 0.00 / 0.44 / 0.48 [1.0 / 0.72 / 0.37] | 0.48 / 0.34 / 0.07 | 0.17 / 0.56 / 0.28 |
| torus `T_D` | 0.42 / 0.35 / 0.07 | 0.00 / 0.39 / 0.49 | 0.41 / 0.35 / 0.07 | 0.17 / 0.47 / 0.28 |
| `T_6` (ideal, order `Φ_6(r)`) | 0.06 / 0.00 / 0.00 | 0.45 / 0.09 / 0.00 | 0.06 / 0.00 / 0.00 | 0.24 / 0.01 / 0.00 |
| `F_{r^6}^*` (ideal) | 0.00 / 0.00 / 0.00 | 0.16 / 0.00 / 0.00 | 0.00 | 0.02 / 0.00 / 0.00 |

Reading the tables:
* **N-derived exponents.** `log₂ A ≈ −(n/2)·deg Ψ + O(log n)` in every group. Bigger groups (`F_{r^k}^*`, `|G| ≈ p^k`)
  make it **worse**: the coincidence `gcd(f(q), Ψ(p))` stays small while `|G|` grows. Raising the degree of `E`
  (`N^6 − 1` vs `N − 1`) gains about 3 bits.
* **Smooth exponents** split `N` exactly when one hidden order is smooth and the other is not. At 32 bits, `L_100000`
  makes *both* `A = 1`, which gives gcd `= N` and leak 0. Pollard's method fixes this by taking gcds stage by stage.
* **Products** `(N^a − 1)·L_B` are no better than `L_B` alone (to within 0.01).
* **Random elliptic curves** at 32 bits (1 200 curves, cyclic approximation `gcd(E,#E)/#E`): `N ± 1` give a median leak of
  `2^{−13}`, while `L_1000` and `(N+1)·L_1000` give `P(leak ≥ ½) = 0.50`. Each new curve is a fresh smoothness trial.

**Ranking of ways to "combine checks" (all prior art; nothing here is new).**
(a) **Multiply exponents** (Pollard `p−1` stage 1/2, Williams `p+1`, Bach–Shallit cyclotomic): success iff the hidden
    group order is smooth. Cost `B` per group; subexponential only with (b).
(b) **Re-roll the group** (Lenstra ECM): each curve gives an independent smoothness trial. `L_p[1/2, √2]`, the best method
    when `p` is small.
(c) **Combine hidden-`s` DLP residues** (Fermat + torus): same unknown, no gain (§3b). Cost `N^{1/4}` (Pollard 1978;
    Hittmeir 2018 shaves a subexponential factor).
(d) **Weighted one-sided collisions** (Lehman; Hittmeir 2021; Harvey 2021): `a^{Aq+Bp} ≡ a^{AN+B} (mod p)` holds mod `p`
    for all `A, B`, but mod `q` only if `ord_q(a) | (A−B)(p−1)`.
    * Verified on 263 random cases (`verify_codex.check_weighted`, 0 failures). The collision was one-sided in 202/263;
      the two-sided cases are small-`q` order coincidences.
    * One-sided collisions are exactly what the CRT criterion needs. Batching them gives `N^{1/5+o(1)}` deterministic
      (Harvey), the best known *deterministic* bound.
(e) **Relation combining** (congruence of squares: QS `L_N[1/2, 1]`, NFS `L_N[1/3, 1.923]`): many partial facts (smooth
    relations) combined by linear algebra mod 2 into one asymmetric event (`x² ≡ y²`, `x ≢ ±y`).
    * This is the only family that beats exponential time for balanced `N`.
    * It does not use Fermat-type residues at all.

## 6. Candidate cards (rule 9)

**PT-1: factor from the primality residue by interval DLP.**
* Target: balanced semiprimes. Representation: `y = a^{N−1}`.
* Operation: kangaroo/BSGS on `[2√N, 3√N]`, which exposes `s`. Lemmas: Thm 5.
* Cost: `Θ(N^{1/4})` (measured slope 0.24). Prior: Pollard 1978.
* Refuting scaling: slope 0.25 (observed).
* Status: **works; exponential**; baseline only.

**PT-2: iterate or combine residues.**
* Mechanism: `a^{f(s)}`, `f` known.
* Result: degree ≥ 2 turns it into order finding (slope 0.52, measured); combining Fermat + Lucas has the same single
  unknown.
* Status: **refuted as an improvement over PT-1** in the generic model; no sub-`N^{1/4}` combination is known to us.

**PT-3: 2-adic / `ε` leakage from transcripts.**
* Mechanism: MR/SS liar levels.
* Exact rate `≈ g²/(3N)` (proved); battery null (pre-registered).
* Status: **refuted for this model and menu**. The only extractor is `√−1` (Prop 9).

**PT-4: AKS coefficient leak.**
* Mechanism: Kummer sparsity of `(X^p+a)(X+a)^{q−p}` when `q − p ≲ r`.
* Status: **real but confined** to the Fermat-method regime; vanishes for balanced `N` (1/40 instances at 40 bits).

**PT-5: favourable exponent/group.**
* N-derived exponents are q-driven (Lemma M2), giving `A ≈ p^{−deg Ψ + o(1)}`.
* Smooth exponents reduce to `p−1` / `p+1` / ECM.
* Status: **no improvement over known methods**.

## 7. Open question (single most promising)

Is there a **one-sided** (weighted, Lehman/Harvey-type) collision structure inside a primality-test computation? It would
need intermediate values that collide modulo `p` but not `q` with probability `≫ g/p`. Every test analysed here produces
only symmetric residues (functions of `s`, or of `±(q−p)`). Harvey's batching is the known way to exploit one-sided
collisions (`N^{1/5}`); whether a primality-style identity gives a better collision family is **not yet tested**.

## 8. Reproduction (from the project root; 32 cores; total < 5 min)

```
.venv/bin/python -m pytest -q tests/test_primality.py                       # 20 tests, ~2 s
.venv/bin/python -m src.primality.verify --bound 20000 --three-bound 20000 --carm-max 200000 \
    --lucas-bound 12000 --frob-bound 700 --torus-bound 3000                  # ~31 s -> experiments/raw/primality_verify.json
.venv/bin/python -m src.primality.verify_codex --bound 12000 --lucas-bound 8000   # -> primality_verify_codex.json
for c in epsexact aks kangaroo iterate epscert leakrates battery menu; do
  .venv/bin/python -m src.primality.experiments $c --workers 30; done      # raw -> experiments/raw/primality_*.jsonl
gzip -f experiments/raw/primality_{battery,leakrates,menu,epscert}.jsonl
.venv/bin/python -m src.primality.analyze                                   # -> experiments/analysis/primality_tables.md
```

Operational vs harness split:
* `ops.py` and `dlog.py` see only `N`, parameters and seeds.
* `formulas.py`, `verify*.py`, `experiments.py`, `menu.py` and `analyze.py` use factors and are harness/diagnostic code.
* `dlog` uses `gmpy2.is_square` only on candidate `s² − 4N`.

## 9. Verification of GPT-6.1-sol r6 (`notes/codex_r6_primality_lit.md`)

All items were checked exhaustively over every unit base for all 2 319 odd prime pairs with `pq < 12 000`
(`src/primality/verify_codex.py`, `experiments/raw/primality_verify_codex.json`):

| GPT item | result |
|---|---|
| (1) gcd formula for `e ∈ {N−1, (N−1)/2, N+1, m}` | agrees, 0 failures |
| (2) `g_p = g_q = g`; `t = k` iff `u ≠ v` | agrees |
| §2 table: Fermat pass/gcd; MR pass, `√1`, augmented; SS pass, residual gcd, both-gcd (`u = v` / `u ≠ v`) | all 10 entries agree |
| (8) `P(x_i = 1)`, `P(x_i = −1)` | agrees (all `i`) |
| (10) MR `√−1` probability `F(4^k−4)/(3·4^k)` | agrees |
| (11) `P(b = −1) = 0` (`u = v`), `F/4` (`u ≠ v`) | agrees |
| Carmichael `P(no factor) ≤ 2^{1−r}` | holds for all 33 listed + 12 Chernick (via the exhaustively verified engine); equals the strong-liar probability |
| (5) P-sampled strong Lucas | agrees on 1 069 pairs `< 8 000` |
| (7) torus table + `Q^K` normalisation | agrees on 1 248 cases, all 4 patterns |
| `z^m = (V_m + √D U_m)/(V_m − √D U_m)` | agrees (200 random) |
| weighted identity mod `p` | agrees (263 cases) |

**Notes and minor disagreements.**
1. GPT's "augmented MR" omits `gcd(x_t + 1, N)`. Including it (our `pm`) adds the event "`a^{N−1} ≡ −1` at exactly one
   prime" (`101·103`: 299/5 100 vs 199/5 100). By Lemma 4 this can never happen at both primes.
2. GPT's Carmichael statement is for `r ≥ 3`. For 2-prime `N` with `λ | N−1`, which cannot occur for squarefree
   Carmichael, the formula is the same.
3. GPT's citations (Pollard 1978 p. 923, Hittmeir arXiv 1608.08766, Harvey arXiv 2010.05450, Long TR 284) were **not**
   re-read by this agent; they are reported as given.

## 10. Status summary

| claim | status |
|---|---|
| Thm 0–3, Cor 3a–3d, Frobenius kernels, Lemma 4 | proved; verified exhaustively |
| Thm 5 torus table, `Q^K` normalisation | proved; verified |
| Thm 7 AKS twist / period / norm reduction | proved; twist + period verified |
| AKS Kummer-regime leak count `64(r − cov)` | mechanism proved; counts experimentally supported |
| Leak rate `≈ 2g/√N` on balanced semiprimes | proved (formula) + experimentally supported (32–128 bits) |
| `N^{1/4}` cost of residue → factors | experimentally supported (slope 0.24); generic bound is standard |
| Iterated residues are no easier | Prop 6 proved; measured slope 0.52 |
| Joint Fermat + Lucas DLP needs `Ω(N^{1/4})` generically | conjectured |
| Transcript patterns carry `≤ O(g/p)` TV about `ε` | proved (Thm 8) |
| Battery null for 14 features at 32/48/64 bits | experimentally supported (MDE 0.029) |
| Prop 9 (2-adic extraction ⇒ `√−1` or factor) | proved |
| Thm 10 meta-argument | proved |
| Lemma M1, M2; `A_p ≈ τ(Ψ(p))/Ψ(p)` average | M1/M2 proved; average heuristic, experimentally supported |
| No one-sided collision family inside primality tests beyond weighted sums | not yet tested |
