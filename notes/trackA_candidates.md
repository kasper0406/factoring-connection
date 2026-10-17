# Track A candidate cards (agent infraA)

Template: BRIEF §4. Proofs, data and commands are in `notes/trackA_analysis.md` (cited as §n).
No novelty claims. Status labels per AGENT_RULES.

---
### A-1. Gate-local E/O sectioning of the repeated-squaring circuit + formal PIT
- **Mechanism:** Section the `O(n)`-gate circuit of `F_N` with the E/O gate rules (hash-consed), select a nonzero child by random-quotient PIT, and repeat `n` levels to reach a support index.
- **Target / family:** any composite `N`; tested on balanced semiprimes 16/24/32 bits.
- **Representation:** hash-consed DAG over `R_N[Y]` (`src/circuit_sections.py`).
- **Extraction:** a nonzero child certified by evaluation, then `gcd(index, N)`.
- **Why a factor:** L2 (`k a_k = 0`).
- **Known lemmas:** L1, L2, L7; PIT over `R_N` (Agrawal–Biswas; AMS 0801.0514).
- **New lemma required:** a poly-size bound on the section circuits along a path. **Refuted for this transform:** nodes grow ×3.6–3.7 per level (naive gate, → 4) and ×3.0 (Karatsuba gate) on all paths (§2).
- **Correctness:** proved (tests vs explicit coefficients).
- **Cost in n:** `Θ(3^k)`–`Θ(4^k)` per squaring gate at level `k`, so exponential.
- **Randomness:** PIT only.
- **Closest prior:** the BRIEF conditional reduction and codex r1 A1. Difference: we measured it.
- **Smallest experiment:** done (`run_circuits`).
- **Refuting behaviour:** a constant growth factor > 1 per level (observed, 3.0–3.7).
- **Repair:** a non-gate-local, binomial-specific identity (escape E-1, §3.4, §9). None found (§6). **Status: refuted for this model.**

### A-2. Exact quotient-evaluation selector (rank 2M)
- **Mechanism:** At level `k`, compute `(1+X)^N mod (X^{2M} - c)` for random `c`. Coefficients `r` and `r+M` are the two children at `c`. A nonzero value certifies that child.
- **Target:** balanced semiprimes 16–32 bits (`run_descent`).
- **Representation:** dense element of `R_N[X]/(X^{2M}-c)` (flint).
- **Extraction / why a factor:** the path reaches a support index, then gcd (L2).
- **Known lemmas:** C1, C2, C4 (§4).
- **New lemma required:** none for correctness. ~~For cost it would need a rank `< 2M`, which is impossible in the algebra-observable model (Theorem A1, proved).~~ **Corrected 2026-10-17 after r3 audit:** A1 (corrected form) only gives rank `≥ min(2M, i+1)` for an isolating certificate with support index `i`, i.e. `≥ p+1` in the balanced family (Cor. 3′). A cheaper isolating certificate of rank about `p+1` at high levels is not excluded, but it is still `Ω(p)` under dense storage.
- **Correctness:** experimentally 25/25 paths alive at every level. The first scalar always certified (C4 explains why).
- **Cost:** rank `2^{k+1}` at level `k` *for this selector*. The cost per level doubles exactly (ratio 2.00), total `Θ(N log N)`, so exponential and worse than trial division. (The doubling is a property of this selector, not a theorem.)
- **Closest prior:** the direct coefficient sketch / AKS quotients.
- **Refuting behaviour:** cost doubling, observed. **Status: correct but exponential (experimentally supported + proved lower bound in model).**

### A-3. Fixed-budget (non-isolating) selector, rank D independent of level
- **Mechanism:** use rank-`D` evaluations for every level.
- **Result:** ~~for `2M > D` the two child classes coincide modulo any rank-`D` observable's resolution (A1). The measured success equals a coin flip.~~ **Corrected 2026-10-17 after r3 audit:** A1 shows only that a rank-`D` observable cannot *isolate* the child classes (nonzero index `≥ 2M`, or `≥ p+1` in the balanced family). Failure to isolate is not indistinguishability. Over `F_101`, `f ↦ f(2)` (rank 1) separates `1` from `X` without isolating either class. A1 gives no coin-flip bound. The implemented `budget_1024` selector flips a literal coin once `2M > D` (`BudgetSelector` calls `rng.getrandbits(1)`), so its measured coin-flip rate (dies 1–4 levels after `max(log2 D, log2 p)`) is by construction. **Status: refuted for the implemented selector (experiment). The general claim "every rank-`D` selection rule is a coin flip" is not proved.**

### A-4. Parity / descent through `X^{2M} - s^2 = (X^M - s)(X^M + s)`
- **Mechanism:** a child vanishes iff the parent is odd/even in `Y`. Test this from the parent at `±s` (two rank-`M` algebras) or in the tower `R_N[s]/(s^2-c)[X]/(X^M - s)`.
- **Result:** proved identities (`descent_split`, NonUnit-checked). The rank is `M + M = 2M`. A1 Cor. 2 (corrected 2026-10-17) gives only `≥ min(2M, i+1)` (`≥ p+1` in the balanced family), so this construction is not proved optimal. Obtaining `P(-c)` from `P(c)` needs a primitive `2M`-th root of unity (rank `M` extension). **Status: refuted for this model (no saving).**

### A-5. Section-value gcd ("pure classes")
- **Mechanism:** `gcd(S_{r,M}(c), N)` is a proper factor when class `r` is hit by only one prime's Lucas support.
- **Why:** the value is `≡ 0` mod the other prime (C3).
- **Lemma:** pure classes need `M > min(p-1, q-p)` (proved). The first hit is at `2M ∈ (q-p, 4(q-p)]` (experimentally, 25/25). Counts: 338 STRUCTURAL vs 3 ACCIDENTAL (null ≈ 2.05).
- **Cost:** `Θ((q-p) log N)` products of size `≈ q-p`. For random balanced primes `q-p = Θ(p)`, so exponential. For near-equal primes it is dominated by Fermat's method.
- **Closest prior:** Track B coefficient gcd sketches; Fermat. **Status: works structurally, no gain (exponential).** Handed to sketchBD as cross-reference.

### A-6. Masks via Hadamard product (prefix / interval / bit-pattern / class)
- **Mechanism:** extract `P ⊙ m` for a short-description mask `m` by a quotient evaluation.
- **Result:** Corollary E1 (proved). Rank `≥` (zero-run length + 1): class `≥ M`, bit `j` `≥ 2^j+1`, interval `≥ a+1`. Over fields, prefix/interval need rank `≥` last used index + 1, provided a trailing zero block of length `≥` rank is observed (corrected 2026-10-17 after r3 audit; F_5 counterexample `u_t = 1 − 2^t`). The bounds are tight (exhaustive F_2/F_3) and attained by `X^M - c` and truncation. The prefix mask is informative only at `K ≈ p` (36/36), so it costs `≈ p`. **Status: refuted for this model.**

### A-7. Early stop + gcd scan of the progression `r + tM`
- **Mechanism:** descend `k` levels, then scan `gcd(r + tM, N)` for `t ≤ N/M`.
- **Critique:** for `M < q` every progression contains `≈ q/M ≥ 1` multiples of `p` regardless of `r`, so a hit is not evidence for the descent. The balanced cost `min_k (2^k + N/2^k) ≈ √N` matches trial division. **Status: refuted (trivial).**

---
**Surviving open direction:** E-1 only. A poly-size, constructible, binomial-specific program for `S_{r,2^k}((1+X)^N)`, which must not be an algebra observable in the sense of A1. Proposed first experiment: §9 of the analysis.
