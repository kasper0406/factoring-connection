# avgtime: complete deterministic factoring with average cost N^{1/6}·polylog on S_n

Agent `avgtime`, 2026-10-18. Inputs: notes/dcltheorem_results.md (Theorem D′, authoritative),
codex/r16_strategy.out §3 item 2 and codex/r16_strategy_lemmas.md. Code: `src/avgtime/`
(`exponents.py` exact-fraction audit and LPs; `detector.py`, `ladder.py` operational).
Tests: `tests/test_avgtime.py` (21 tests). No experiments were run, so there is no
pre-registration; every number below is exact rational arithmetic or a proof.

## 0. Verdict (plain language)

**GPT-6.1-sol's 2/11 is correct for the design it analyses, but that design is not optimal.
Re-running D′ with a growing weight height (a "ladder"), plus a cheap detector for primes q
with small ord_q(2), gives a complete deterministic algorithm whose average cost over S_n
is N^{1/6}·(log N)^{O(1)}. That is the cost of a single D′ run, so the average-case exponent
is 1/6, the floor of the mechanism.** The worst case is still Harvey–Hittmeir's N^{1/5}.

- **Theorem AV (proved; it uses D′ with H as a free parameter, plus three short lemmas R, S, C).**
  There is a deterministic algorithm LAD with three properties.
  - It outputs the complete prime factorisation of every N ≥ 2.
  - Its worst-case cost is HH(N) + O(N^{7/36}·polylog), where
    HH(N) = O(N^{1/5}log^{16/5}N/(log log N)^{3/5}) [Harvey–Hittmeir 2022].
  - Its average cost over uniform pairs of S_n = {pq : p < q ≤ 2p} is
    O(N^{1/6}(log N)^{c+1}·log log N). Here c is the polylog exponent of D′'s cost
    bound. Using D′'s own n²/H² coverage bound instead of the n/H² refinement, the
    factor is (log N)^{c+2}·log log N.
- **Why the ladder works.** One D′ run with weight height H costs ≈ N^{1/6}H^{5/3}. It
  fails on at most ≈ n/H² of S_n. The failure fraction decays faster in H (exponent 2)
  than the cost grows (exponent 5/3). So rerunning on FAIL with H ← H^{6/5} adds only
  O(N^{1/6}·polylog) per rung, over O(log log N) rungs. One HH run at the top costs
  N^{1/5}·N^{−1/30} = N^{1/6}.
- **Why the detector is needed (within the proved bounds).** Without it, the bound on the
  two-sided set E_c caps the failure exponent of one D′ run.
  - Even with the sharper dyadic recount of E_c (Lemma D, due to the Opus critic), the cap
    is e ≤ 1/4 − 3h − μ, which is 1/40 for the balanced design.
  - The no-detector ladder then has infimum **N^{17/100+o(1)}**. With D′'s one-threshold
    E_a/E_c accounting it is N^{19/110}. Both are proved within the stated bounds, by an
    elementary LP and a closed form. They are limits of these bounds, not lower bounds on
    the algorithm.
  - The detector is a baby-step giant-step search for the smallest k with
    gcd(2^k − 1, N) > 1, run up to T = N^{3/10}. It costs N^{3/20}.
  - It factors every pair with min ord ≤ T except ord_p 2 = ord_q 2. Those coincidences
    number O(T³/P²) = O(N^{−1/10})·|S_n|.
  - The remaining pairs have ord_q(2) > N^{3/10} > N^{1/4}, which pushes the E_c count down
    to H³m/T.
- **Task item 1: the audit of ε.**
  - D′'s written proof (T = 2^{(1/4−ε)n}) **certifies every ε ∈ (0, 1/52)**. The binding
    estimate is E_c1: −1/12 + 13ε/3 < 0. The other constraints are Y > H² (ε < 1/7),
    Y⁴ > √2P (ε < 1/4) and E_c2 (ε < 1/16). E_a never binds.
  - These are ranges certified by the estimates. They do not show that D′ fails outside
    them.
  - So ε = 1/55 is admissible, and D′ + HH gives E[T] ≤ N^{2/11}polylog. **GPT-6.1-sol's
    claim is confirmed.**
  - **Re-choosing T** (an analysis-only parameter) certifies ε < 1/46.
  - **Lemma D**, which sums the E_c count dyadically over ord_q(2), certifies the
    *unchanged* D′ for **ε < 1/40**. It also removes E_a from D′ entirely.
  - **With the detector**, ε < 1/22. That covers the whole useful range, since D′ beats HH
    only for ε < 1/25.
- **Task item 2: re-optimisation.**
  - For "one D′ run, then HH", **2/11 is the exact optimum** under the proved bounds, even
    with H, Y, T and the baby length m all free and the detector included. The cost is
    ≥ 1/6 + 5h/3 for any Y and m; the coverage term is N^{−2h}. All of this is within the
    proved upper-bound model, not a lower bound on run time.
  - Laddering is the change that matters:

    | design | exponent |
    |---|---|
    | single HH fallback | 2/11 = 0.1818 |
    | ladder, no detector | 17/100 (19/110 with the one-threshold count) |
    | ladder + detector | **1/6** (polylog) |

- **Task item 3: conditional tails.**
  - For the single-fallback design, a coverage tail λ^{−j} gives 1/6 + 1/(30 + 36j):
    3/17 for j = 2 and 4/23 for j = 3. A tail e^{−cλ} gives 1/6.
  - For the ladder + detector the exponent is already 1/6 with the **proved** j = 1. Any
    tail λ^{−j} with j > 5/6 suffices, with a modified schedule (§6). A stronger tail
    (sieve2d) improves only the polylog factor of the average.
  - With j ≥ 2 the detector also becomes unnecessary (§6).
  - Its value lies elsewhere: the fixed-cost almost-all statement (D′'s trade-off N^{1/6+δ}
    vs N^{−f(δ)} exceptions) and the worst case.
- **Task item 4: prior art.** A deterministic complete algorithm with proved average cost
  below N^{1/5} on balanced semiprimes was **not found by us** (own WebSearch plus
  GPT-6.1-sol avgtime_2, 20 items, §8). The ingredients are standard:
  - restarting with growing parameters;
  - BSGS order search;
  - D′ itself.
- **Not claimed.**
  - No worst-case improvement.
  - No result for ratios q/p > 2.
  - No lower bound: we do not show that 1/6 is optimal for average-case factoring.
  - No novelty claim.

## 1. Setting and the free-H form of D′

Notation as in notes/dcltheorem_results.md: n bits, P = 2^{n/2}, S_n = {(p,q) : p < q ≤ 2p
primes, 2^{n−1} ≤ pq < 2^n}, |S_n| ≥ c₀P²/n². Exponents are in units of n = log₂N, so "x"
means N^x up to (log N)^{O(1)}.

DCL_{n,H} denotes D′'s algorithm (§2 of dcltheorem_results) with weight height H as a free
n-only integer and Y = ⌈2^{n/6}H^{−1/3}⌉, four bands, ℓ ≡ 3 (mod 4), base 2, m = m̄_k (or the
endpoint-competitive rule). For a threshold T (analysis only) define, as in D′ §4(iii):

- E_a(T) = {(p,q) ∈ S_n : ord_q(2) ≤ T};
- E_b(H) = {ord_q(2) > T, p ∈ U_{k(q)}(q)} (coverage failure; independent of T as a set
  bound);
- E_c(H,T) = {ord_q(2) > T, p ∈ F_{k(q)}(q)} (two-sided true giants).

**Proposition U (proved: D′'s proof with every "n ≥ n₀(ε)" step re-checked uniformly in H).**
There are absolute constants n₀, H₀, C, c such that for n ≥ n₀, H₀ ≤ H ≤ 2^{n/50} and
1 ≤ T ≤ P/2:

- (i) DCL_{n,H} outputs a proper divisor of N or FAIL on every input.
- (ii) It costs ≤ C n^c H^{5/3} 2^{n/6} bit operations on every n-bit input.
- (iii) Every FAIL pair (p,q) ∈ S_n with ord_q(2) > T lies in E_b(H) ∪ E_c(H,T), and
  - |E_b(H)| ≤ C (n²/H²)|S_n| (D′, Lemma 9), and ≤ C (n/H²)|S_n| by Lemma R below;
  - |E_c(H,T)| ≤ C n (H^{14/3}P^{1/3}/T + H^{10/3}/P^{1/3})|S_n| (D′, Lemma 10).

*Proof (uniformity check).* D′'s proof uses largeness of n only in the following places.
Each holds uniformly for H₀ ≤ H ≤ 2^{n/50} and n ≥ n₀:

| use in D′ | requirement | holds because |
|---|---|---|
| \|V\| ≥ H²/2 | H ≥ H₀ | coprime pairs (6/π² + o(1))H² |
| Lemma 6 (a clean band exists) | Y⁴ > √2P | Y⁴ ≥ 2^{2n/3}H^{−4/3} ≥ 2^{(2/3−4/150)n} > 2^{n/2+1} |
| Lemma 9 (classes distinct) | Y > H² | H^{7/3} ≤ 2^{7n/150} < 2^{n/6}/2 |
| Lemma 9 (PNT mod 4) | Σ_{ℓ∈𝒫_k}1/ℓ ≥ c′/n | Y ≥ 2^{n/6−n/150} → ∞ |
| Lemma 9 (sieve level) | Q² = 4^{k+1}Y² ≤ 256Y² = O(P) | y ≤ 1/6 |
| Lemma 1 (ℓ ∤ p) | 16Y < P/2 | y < 1/2 |
| Lemma 10 / m̄ bound | m̄_k ≤ 3HP/(4^kY²) | HP/Y² ≥ H^{5/3}2^{n/6} ≫ 2 |
| cost table | H²n, Y, nT_k ≤ n^{O(1)}H^{5/3}2^{n/6} | H^{1/3} ≤ 2^{n/6} |
| counts | #q ≤ 3P/n, \|S_n\| ≥ c₀P²/n² | PNT |

Lemmas 1–5, 7, 8 and 10 have no size hypotheses beyond these. In D′ the constant ε enters
only through H = n2^{εn/2} and T = 2^{(1/4−ε)n}. Both are now free, and the bounds of D′
§4(iii) are stated before substituting them. ∎

**Lemma R (coverage O(n/H²); codex r16 §3.1; verified here).** Under the hypotheses of
Proposition U, |E_b(H)| ≤ C(n/H²)|S_n|.

*Proof.*
- Fix q and the clean band k = k(q). Let U be the set of **primes** p ∈ (P/2, P) covered
  by no target of band k. Put Q = ⌊P^{1/2}⌋ and Z = ⌊Q/(32Y)⌋. Then Z ≥ 2^{n/12}/64 → ∞,
  and Z < Y because Y² > 2^{n/4}/32.
- U avoids the s = |V| distinct classes AqB^{−1} mod each ℓ ∈ 𝒫_k (D′ Lemma 9; s < ℓ). It
  also avoids the class 0 mod every prime r ≤ Z, since p > P/2 > Z.
- Apply the arithmetic large sieve, in Montgomery's form for squarefree composite moduli
  (Montgomery 1968; the constant N₀ + Q² is Montgomery–Vaughan 1973), with N₀ = P/2 + 1 and
  moduli ≤ Q:

    |U| ≤ (N₀ + Q²)/J, J = Σ_{d≤Q} μ²(d)∏_{r|d} ω(r)/(r − ω(r)), so |U| ≤ (1.5P + 1)/J.

  D′'s Lemma 9 used only prime moduli. Lemma R needs the composite d = ℓa.
- Keep only the terms d = ℓa of J, with ℓ ∈ 𝒫_k and a ≤ Z squarefree. These are squarefree
  because every prime of a is ≤ Z < Y < ℓ, and d ≤ 16Y·Q/(32Y) ≤ Q. Hence
  J ≥ (Σ_ℓ s/(ℓ−s))·Σ_{a≤Z}μ²(a)/φ(a) ≥ λ·ln Z, with λ = sΣ_ℓ 1/ℓ ≥ (H²/2)(c′/n).
  Here Σ_{a≤Z}μ²(a)/φ(a) = Σ_{rad(m)≤Z}1/m ≥ Σ_{m≤Z}1/m ≥ ln Z. The identity holds because
  1/φ(a) = ∏_{r|a}(1/r + 1/r² + …).
- Since ln Z ≥ c″n, J ≥ c‴H² and |U| ≤ C P/H².
- Sum over the ≤ 3P/n primes q and divide by |S_n| ≥ c₀P²/n². This gives
  |E_b|/|S_n| ≤ C n/H². ∎

The refinement only changes the polylog in Theorem AV, never an exponent.

**Proposition U′ (free baby length and band scale; proved, Opus critic item 9).**
- **Statement.** Fix a baby length m ≤ m̄_k (n-only) instead of m̄_k or the
  endpoint-competitive rule, and any n-only Y with H² < Y, Y⁴ > √2P and 256Y² ≤ P.
  Proposition U then holds with:
  - cost ≤ C n^c (m + T_k(1 + m̄_k/m));
  - the same E_b bound;
  - the E_c bound with m̄_k replaced by m. That is
    |E_c|/|S_n| ≤ C n (H³m/T + H²Ym/P), in exponents 3h + μ − t and 2h + y + μ − 1/2.
- *Proof.*
  - The cost table of D′ §4(ii) charges m babies and Σ_τ⌈W_τ/m⌉ ≤ T_k(1 + m̄_k/m) giants.
  - Lemma 7's δ-range is contained in (−m, m), so Lemma 10 holds with m̄_k replaced by m.
  - The separations are exactly those used in Proposition U. ∎
- This is the model that the LPs of §5 optimise (`exponents._stage`, which also enforces
  m ≤ m̄, i.e. μ ≤ h + 1/2 − 2y).

**Lemma D (dyadic recount of E_c; no E_a; proved, Opus critic item 6).**
- **Setting.** Drop the threshold T. Put F := {(p,q) ∈ S_n : p ∈ F_{k(q)}(q)}, where F_k(q)
  is defined as in D′ Lemma 10 with d = ord_q(2) arbitrary.
- **Statement.** Every FAIL pair of DCL_{n,H} lies in E_b ∪ F. With A := H³m,

    |F|/|S_n| ≤ C n (A/√P + A²/P + H²Ym/P).

*Proof.*
- **Success needs no order hypothesis.** D′'s success implication (§4(iii)) uses Lemmas 5
  and 7 on the clean band k(q). Neither needs ord_q(2) > T: Lemma 7 uses only
  gcd(L, q − 1) = 1. So FAIL ⊆ E_b ∪ F.
- **Per-q count.** Lemma 10's count for a fixed d holds for every d:
  |F_k(q)| ≤ min(P, 2|V||𝒫_k|m(PH/(2^{k+1}Yd) + 1)) ≤ min(P, AP/d) + 2^{k+1}H²Ym.
- **Sum over q, split by d_q.** By Lemma 8, #{q : d_q ≤ D} ≤ 2D²/(n−1), and the total
  number of q is ≤ 3P/n.
  - q with d_q ≤ A contribute ≤ P·2A²/(n−1).
  - The dyadic shells D < d_q ≤ 2D with A ≤ D ≤ √P contribute
    ≤ Σ 8D²/(n−1)·AP/D ≤ 16A√P·P/(n−1).
  - The shells with D > √P contribute ≤ Σ (3P/n)·AP/D ≤ 6AP·√P/n.
  - The "+1" terms contribute ≤ (3P/n)·16H²Ym.
- **Normalise.** Divide by |S_n| ≥ c₀P²/n². ∎

**Consequences.**
- For D′ itself (m = m̄, μ = 1/6 + 5ε/6, h = ε/2) the three terms have exponents
  −1/12 + 7ε/3, −1/6 + 14ε/3 and −1/6 + 5ε/3. Each is < −ε iff ε < 1/40, 1/34 and 1/16
  respectively.
- So the **unchanged** D′ algorithm meets its Theorem D′ guarantee O(N^{−ε}) for every
  ε < 1/40, and its proof needs no E_a term.
- In the ladder with the detector, Lemma D is not needed. There the pairs with
  ord_q(2) ≤ T = N^{3/10} are factored by the detector, which is a stronger cut than the
  dyadic one (T > √P).

**Lemma S (small-order detector; `src/avgtime/detector.py`; elementary, standard BSGS).**
Let N = pq be odd, T ≥ 1, m = ⌈√T⌉, J = ⌈T/m⌉, and d_r = ord_r(2).
- **The algorithm.** Form f(X) = ∏_{0≤i<m}(X − 2^i) over Z/N. For 1 ≤ j ≤ J compute
  G_j = gcd(f(2^{jm}), N) by a product tree and a remainder tree. At the first j with
  G_j > 1, compute gcd(2^k − 1, N) for every k ∈ ((j−1)m, jm].
- **What it finds.** If min(d_p, d_q) ≤ T and d_p ≠ d_q, it outputs the prime with the
  smaller order.
- **Cost.** O(√T) ring operations plus trees, which is Õ(√T·n) bit operations. The output
  is validated.

*Proof.*
- G_j > 1 iff some r ∈ {p, q} has d_r | jm − i for some 0 ≤ i < m. Since 2 is a unit,
  2^{jm} − 2^i = 2^i(2^{jm−i} − 1). So G_j > 1 iff the block ((j−1)m, jm] contains a
  positive multiple of d_p or d_q.
- The least such multiple is d_min = min(d_p, d_q). So the first hit is the block
  j₀ = ⌈d_min/m⌉ ≤ J.
- At k = d_min, gcd(2^k − 1, N) is divisible by the prime of order d_min. It is not
  divisible by the other prime, whose order d_max > d_min does not divide d_min. ∎

The test `test_detector_exhaustive_small_semiprimes` checks the exact success
characterisation, using the scanned range J·m ≥ T, on all pq with 3 ≤ p < q < 260 and
T ∈ {1, 7, 30, 100}. There were 0 mismatches.

**Lemma C (coincidences).** For T ≥ 1 let C(T) = {(p,q) ∈ S_n : ord_p 2 = ord_q 2 ≤ T}.
Then |C(T)| ≤ C (T³/P²)|S_n|.

*Proof.* For each k ≤ T, 2^k − 1 has fewer than k/log₂(P/2) prime factors > P/2 (D′
Lemma 8). So the number of pairs with both orders equal to k is ≤ (2k/(n−2))². Summing over
k ≤ T gives ≤ 4T³/(n−2)² (for T ≥ 1). Dividing by |S_n| ≥ c₀P²/n² gives the claim. ∎

## 2. The algorithm LAD_n (operational reference: `src/avgtime/ladder.py`)

All parameters are n-only.

- T = 2^{⌊3n/10⌋}.
- H_1 = ⌈√n⌉, H_{i+1} = min(⌈H_i^{6/5}⌉, H_max), with H_max = n²·2^{⌈n/60⌉}. The ladder
  stops at H_k = H_max. Then k = O(log n), and H_{i+1}^{5/3} ≤ 2^{5/3}H_i² for every i
  (`exponents.ladder_schedule`, tested for n = 256 … 4096).
- n₀ is absolute and is chosen so that Proposition U applies to every rung. That needs
  H_1 ≥ H₀ and H_max ≤ 2^{n/50}; the second holds only for n ≳ 8000. So the schedule tests
  at n ≤ 4096 check the arithmetic of the schedule, not this hypothesis.

Input N:
- If n < n₀, run HH and stop. This is a finite set of sizes, so the cost is O(1).
- If N < 4, the output is trivial.
- If N is even, output 2 together with HH's factorisation of N/2. The worst case is
  unchanged, and S_n contains no even N.
0. Run the detector (Lemma S) with T.
1. For i = 1, …, k, run DCL_{n,H_i}.
2. If neither step 0 nor step 1 produced a divisor, run Harvey–Hittmeir 2022, which outputs
   the complete factorisation of N.

**Finish.** As soon as a step outputs a proper divisor d, test d and N/d for primality with
AKS (deterministic, polylog). If both are prime, output {d, N/d}. Otherwise run HH on N.
- On S_n the finish never calls HH, because d ∈ {p, q}.
- In the reference implementation the fallback is Lehman's deterministic N^{1/3} method and
  primality is gmpy2's BPSW-type test. These are stand-ins for correctness tests at
  32–40 bits only.
- The implementation's toy ladders (H ∈ {2, 3, 4}) are below n₀. They test the plumbing,
  not the theorem.

## 3. Theorem AV and proof

**Theorem AV.** Let c be as in Proposition U(ii). There are absolute n₀ and C such that:

- (i) **Completeness.** For every N ≥ 2, LAD_n (n = bit length of N) outputs the complete
  prime factorisation of N. It is deterministic and sees only N.
- (ii) **Worst case.** On every n-bit input the cost is
  ≤ HH(N) + C n^{c+4}2^{7n/36} + C n^{O(1)}2^{3n/20}. That is
  O(N^{1/5}log^{16/5}N/(log log N)^{3/5}), Harvey–Hittmeir's bound.
- (iii) **Average over S_n.** For n ≥ n₀,

    (1/|S_n|)·Σ_{(p,q)∈S_n} cost(pq) ≤ C·n^{c+1}·log n·2^{n/6} = N^{1/6}(log N)^{c+1}(log log N)·O(1).

  With D′'s original n²/H² coverage bound instead of Lemma R, the factor is n^{c+2}log n.
- (iv) **Other distributions.** The same bound, times K, holds for any distribution on S_n
  with point masses ≤ K/|S_n|.
  - This includes the RSA-type distribution for even n: two independent uniform n/2-bit
    primes, conditioned to be distinct. It induces uniform unordered pairs on its support,
    which lies in S_n ∪ S_{n−1}, with bounded mass ratio by the PNT.
  - Apply (iii) separately to the two shells; LAD uses each input's own bit length.
  - (iv) says nothing about uniform n-bit integers, or about key-generation rules that
    restrict p and q further.

*Proof.*
- (i) Every step outputs only validated proper divisors. On S_n a proper divisor is p or q.
  Otherwise HH factors N completely.
- (ii) Each stage costs ≤ C n^c H_i^{5/3}2^{n/6} ≤ C n^c (2n²2^{n/60})^{5/3}2^{n/6}, with
  H_i ≤ H_max ≤ 2n²·2^{n/60}. There are k = O(log n) stages, so their total is
  ≤ C n^{c+10/3}log n·2^{7n/36} ≤ C n^{c+4}2^{7n/36}, and 7/36 < 1/5. The detector costs Õ(2^{3n/20}), AKS is
  polylog, and HH runs at most once.
- (iii) **Reach sets.** Write R_i for the pairs of S_n that reach stage i, and R_HH for
  those reaching HH.
  - A pair passes the detector only if min(d_p, d_q) > T or d_p = d_q ≤ T (Lemma S). So it
    has d_q > T, or it lies in C(T).
  - A pair with d_q > T that fails stage i − 1 lies in E_b(H_{i−1}) ∪ E_c(H_{i−1}, T)
    (Proposition U(iii)).
  - Hence, for i ≥ 2, R_i ⊆ C(T) ∪ E_c(H_{i−1},T) ∪ E_b(H_{i−1}), and
    R_HH ⊆ C(T) ∪ E_c(H_k,T) ∪ E_b(H_k).
  - Let W ≤ C n^{16/5}2^{n/5} be the worst-case total cost from (ii). Bound the average by
    the five terms below.
  - (a) **Detector and stage 1** (paid by everyone): Õ(2^{3n/20}) + C n^c H_1^{5/3}2^{n/6}
    = O(n^{c+5/6}2^{n/6}).
  - (b) **Stage i ≥ 2, on E_b(H_{i−1})**:
    C n^c H_i^{5/3}2^{n/6}·min(1, C n/H_{i−1}²) ≤ C n^c·2^{5/3}H_{i−1}²·min(1, Cn/H_{i−1}²)·2^{n/6}
    ≤ C′ n^{c+1}2^{n/6}. Summed over k = O(log n) stages: O(n^{c+1}log n·2^{n/6}).
  - (c) **HH on E_b(H_k)**: C n^{16/5}2^{n/5}·Cn/H_max² ≤ C n^{16/5+1−4}2^{n/5−n/30}
    = O(n^{1/5}2^{n/6}).
  - (d) **Pairs in C(T) ∪ E_c(H_i,T)** for some i: these pay at most W.
    - By Lemma C, |C(T)|/|S_n| ≤ C 2^{(9/10−1)n} = C 2^{−n/10}.
    - By Proposition U(iii), for H ≤ H_max: H^{14/3}P^{1/3}/T ≤ n^{28/3}2^{(14/180 + 1/6 − 3/10)n+O(1)}
      = n^{O(1)}2^{−n/18}, and H^{10/3}/P^{1/3} ≤ n^{O(1)}2^{−n/9}.
    - A union over the k = O(log n) stages gives a fraction ≤ n^{O(1)}2^{−n/18}.
    - Contribution: ≤ W·n^{O(1)}2^{−n/18} = n^{O(1)}2^{(1/5−1/18)n} = n^{O(1)}2^{13n/90}
      = o(2^{n/6}).
  - (e) **Finish** (primality tests): polylog.
  - The total is O(n^{c+1}log n·2^{n/6}), provided c ≥ 0. With n²/H², step (b) gives
    n^{c+2}. ∎
- (iv) is immediate from (iii), since the cost is a nonnegative function of N.

**Remarks.**
- The ladder ratio 6/5 is not special. Any H_{i+1} ≤ (H_i²/n^{o(1)})^{3/5}·O(1) works. What
  matters is that the proved tail exponent 2 (in H) exceeds the cost exponent 5/3.
- For a coverage tail (n/H²)^j with 5/6 < j < 1, the ladder needs a modified schedule
  (codex avgtime_1 item 9):
  - H_{i+1} = ⌈H_i^ρ⌉ with 1 < ρ < 6j/5;
  - H_max = n^A·2^{n/(60j)} with a large fixed A;
  - the same T.

  Then 5ρ/3 − 2j < 0, 1/(60j) < 1/50, and the E_c estimates stay below N^{−1/30}.
- For j < 5/6 a geometric ladder cannot reach 1/6, because the rung contributions grow.
  The boundary j = 5/6 still works with a ratio-2 ladder over O(n) rungs (Opus critic).
- The log n factor (number of rungs) can be removed with a slightly faster schedule. That
  is not pursued here.

## 4. Task item 1: audit of every ε-constraint in D′ (exact; `exponents.dprime_requirements`)

D′ as written uses h = ε/2, y = 1/6 − ε/6, t = 1/4 − ε and m = m̄ (μ = 1/6 + 5ε/6). The
target is an exceptional fraction O(N^{−ε}). E_b equals N^{−ε} exactly by the choice of H.
Each remaining requirement is strict, since polylog factors must be absorbed.

| constraint (D′ location) | exponent condition | admissible ε |
|---|---|---|
| Y > H² (Lemma 9, distinct classes) | 2h − y < 0 ⇔ −1/6 + 7ε/6 < 0 | ε < 1/7 |
| Y⁴ > √2P (Lemma 6, clean band) | 1/2 − 4y < 0 ⇔ −1/6 + 2ε/3 < 0 | ε < 1/4 |
| sieve level 256Y² ≤ P, band < P/2 (Lemmas 9, 1) | y < 1/4 | all ε |
| E_a: n·T(T+1)/P (Lemma 8) | 2t − 1/2 + ε = −ε < 0 | all ε > 0 |
| E_c1: n·H^{14/3}P^{1/3}/T (Lemma 10, main term) | −1/12 + 10ε/3 + ε < 0 | **ε < 1/52 (binding)** |
| E_c2: n·H^{10/3}/P^{1/3} (Lemma 10, "+1" term) | −1/6 + 5ε/3 + ε < 0 | ε < 1/16 |
| (usefulness, not validity) cost < HH | 1/6 + 5ε/6 < 1/5 | ε < 1/25 |
| *Lemma D (dyadic recount, replaces E_a + E_c1):* n·H^{14/3}P^{1/3}/√P | −1/12 + 7ε/3 + ε < 0 | ε < 1/40 |
| *Lemma D:* n·A²/P | −1/6 + 14ε/3 + ε < 0 | ε < 1/34 |

- **Conclusion.** D′'s written estimates certify Theorem D′ for every fixed ε ∈ (0, 1/52),
  with n₀(ε) → ∞ as ε → 1/52.
  - The stated ε ≤ 1/60 is conservative; D′'s own "where the constants come from"
    paragraph already says so.
  - With Lemma D, the unchanged algorithm is certified for every ε ∈ (0, 1/40).
  - None of these thresholds is shown to be sharp for the algorithm. At ε = 1/52 the E_c1
    estimate keeps a factor n^{17/3}, which is a failure of the estimate, not of D′.
- **At ε = 1/55.** The exponents are cost 2/11, E_a = −2/55, E_b = −1/55, E_c1 = −1/44 and
  E_c2 = −3/22. All of them are < −1/55 except E_b, which equals −1/55 by construction.
- **So for D′ + HH** the average satisfies
  E[T] ≤ C n^c 2^{(2/11)n} + HH(N)·C2^{−n/55} = N^{2/11}(log N)^{O(1)}.
  **GPT-6.1-sol's claim is confirmed.** It needs ε < 1/52; the margin at 1/55 is
  −1/44 vs −1/55.
- **With T re-chosen.** E_a + ε < 0 needs t < 1/4 − ε/2. E_c1 + ε < 0 needs
  t > 1/6 + 10ε/3. Both hold iff **ε < 1/46**.
- **With the detector** (E_a replaced by C(T), whose exponent is 3t − 1). Now
  1/6 + 10ε/3 < t < (1 − ε)/3 iff **ε < 1/22**. E_c2 (1/16) and Y > H² (1/7) do not bind,
  and the detector costs t/2 < 1/6. So the detector-augmented D′ holds on all of
  (0, 1/25), the useful range.

## 5. Task item 2: re-optimising all free parameters

### 5.1 One D′ run, then HH: 2/11 is optimal (proved, within the proved bounds)

The free parameters are h, y, μ (baby length m = N^μ ≤ m̄) and t.
- **Cost.** The proved bound m + T_k(1 + m̄/m) is ≥ max(T_k, 2√(T_k m̄)). It is tight up to
  constants, because W_τ ≍ m̄ for most v. In exponents,
  cost ≥ max(2h + y, (3h + 1/2 − y)/2), which is minimised at y = 1/6 − h/3, giving
  1/6 + 5h/3, for every y and μ.
- **Coverage.** The coverage term N^{−2h} does not depend on y, μ or t: J ≍ H²/n for any
  band scale.
- **Hence** E[T] exponent ≥ min_h max(1/6 + 5h/3, 1/5 − 2h) = 2/11, attained at h = 1/110.
  This is ε = 1/55, where E_a and E_c are slack.
- **LP check.** The exact LP (`exponents.single_fallback`, simplex over Fractions) returns
  exactly 2/11, with or without the detector. Choosing T, Y or m differently cannot help,
  and neither can a cheaper E_a fallback, because E_a never binds here.

### 5.2 Ladders (proved; exact LP values)

A k-rung ladder needs:
- cost_i ≤ z + e_{i−1}, because rung i is reached by ≤ N^{−e_{i−1}} of S_n;
- HH: 1/5 − e_k ≤ z.

Rung i has failure exponent e_i ≤ min(2h_i, E_c terms, E_a or C(T) term). The model is
Proposition U′. There are three accountings of the failure set:
- "one-threshold" is D′'s E_a/E_c split;
- "dyadic" is Lemma D;
- "detector" uses C(T) and E_c at the shared t.

| k rungs | no detector, one-threshold | no detector, dyadic (Lemma D) | with detector |
|---|---|---|---|
| 1 | 2/11 = 0.181818 | 2/11 | 2/11 |
| 2 | 35/198 = 0.176768 | 16/91 = 0.175824 | 16/91 = 0.175824 |
| 3 | 223/1278 = 0.174491 | 213/1228 = 0.173453 | 116/671 = 0.172876 |
| 4 | 0.173535 | 0.171962 | 0.171146 |
| 6 | 0.172850 | 0.170527 | 0.169247 |
| 8 | 0.172744 | 0.170120 | 0.168269 |
| 12 | 0.172728 | 0.170003 | 0.167354 |
| ∞ | 19/110 = 0.172727 | **17/100** | **1/6** |

- **No-detector limit, one-threshold count (closed form; branch 3 added after
  avgtime_1 item 1).** For a rung with failure exponent e, put h = e/2 and optimise t, μ
  and y. The minimum cost C(e) is:
  - 1/6 + 5e/6 for e ≤ 1/46 (the balanced design);
  - 1/8 + 11e/4 for 1/46 ≤ e ≤ 1/14. Here E_a and E_c1 force μ ≤ 1/4 − 3e, with
    t = (1/2 + 3h + μ)/3;
  - 9e/2 for 1/14 ≤ e ≤ 1/12, where y ≤ 1/4 binds.

  `exponents.stage_cost_closed` equals the LP at 8 test points.
- **Reachability, for z > 1/6, by induction.** Since C is increasing, e_{i−1} ≤ e*
  implies e_i ≤ e*, where e* solves C(e*) − e* = z. Conversely, every e < e* is reachable
  in finitely many rungs.
  - The final condition 1/5 − z ≤ e* = 4(z − 1/8)/7 gives **z ≥ 19/110**.
  - At the limit, e* = 3/110, h = 3/220, μ = 37/220, y = 19/110 and t = 13/55. The last
    rungs cost N^{1/5}, as much as HH. The infimum is not attained.
- **No-detector limit, dyadic count (Lemma D; Opus critic item 6).** Lemma D removes the
  E_a constraint and replaces E_c1 by 1/4 − 3h − μ (and 1/2 − 6h − 2μ, which is implied).
  - With h = e/2: C(e) = 1/6 + 5e/6 for e ≤ 1/40, and C(e) = 1/8 + 5e/2 for
    1/40 ≤ e ≤ 1/20 (`stage_cost_closed_dyadic`, which equals the LP).
  - Then C(e) − e = 1/8 + 3e/2 = 1/5 − e gives e* = 3/100 and **z = 17/100**.
  - So **19/110 is an artefact of the one-threshold accounting**. The best no-detector
    value within the proved lemmas is 17/100 (an infimum).
- **With the detector,** the per-rung cap is e ≤ 2h up to h = 1/50 (with t = 3/10: E_c1 needs 1/6 + 20h/3 ≤ 3/10, and C(T) gives −1/10). Hence C(e) − e =
  1/6 − e/6 < z for every z > 1/6, and the LP values tend to 1/6 geometrically (ratio
  ≈ 5/6 per rung).
- **Polylog form.** The finite-k LP is an exponent statement. Theorem AV takes
  k = O(log n) rungs and obtains N^{1/6}·polylog directly.

### 5.3 Small-order primes (E_a): cheaper fallbacks

| fallback | cost | covers | use |
|---|---|---|---|
| HH on E_a (D′ as written) | N^{1/5}·n T²/P | — | never binding for one run |
| gcd(2^k − 1, N), k ≤ T sequentially | T | min ord ≤ T, ord_p ≠ ord_q | T = N^{1/6} only |
| **BSGS (Lemma S)** | √T | same | T up to N^{1/3−} for cost ≤ N^{1/6} |

The coincidence set C(T) of the detector has size T³/P². So T = N^{3/10} gives
N^{−1/10}, which is below the N^{−1/30} needed at the top rung. The E_c1 term then has
exponent 1/6 + 14h/3 − 3/10 ≤ −1/18 for h ≤ 1/60 + o(1).

The detector is what lets the E_c count use T = N^{3/10} > √P. Without it, even the
dyadic count is limited by the q with ord_q(2) ≈ √P, which caps a balanced rung at
e = 1/40 < 1/30.

## 6. Task item 3: conditional coverage tails

Let λ = H²/n. "Tail" means |E_b| ≤ polylog·λ^{−j}|S_n|. The proved value is j = 1 (D′:
j = 1 with n² in place of n).

| tail | single HH fallback | ladder, no detector | ladder + detector |
|---|---|---|---|
| λ^{−1} (proved) | 2/11 | 17/100 (19/110 one-threshold) | **1/6, polylog n^{c+1}log n** |
| λ^{−2} | 3/17 = 0.17647 | 1/6 (proved; dyadic LP 0.166684 at k = 8) | 1/6, polylog |
| λ^{−3} | 4/23 = 0.17391 | 1/6 | 1/6, polylog |
| λ^{−j} | 1/6 + 1/(30 + 36j) | 1/6 for j ≥ 2 | 1/6 for any j > 5/6 (modified schedule, §3) |
| e^{−cλ} | 1/6, polylog (H ≍ n) | 1/6 | 1/6, polylog |

- **Single fallback, λ^{−j}.** Balance 1/6 + 5h/3 = 1/5 − 2jh. The E_a/E_c cap
  1/18 − 28h/9 ≥ 2jh holds at h = 1/(50 + 60j) for every j ≥ 1, so it never binds.
  Formula proved; LP agreement tested for j = 1, 2, 3.
- **Ladder without detector, j ≥ 2 (proved).** With h = e/(2j) and the dyadic count, the
  balanced design's cap is e ≤ 1/4 − 3h − μ.
  - For j = 2 that is e ≤ 1/26, which exceeds the 1/30 needed at the top. So
    C(e) − e = 1/6 + 5e/12 − e < 1/6 + δ along the whole ladder, and the exponent is 1/6.
  - Larger j only lowers h.
  - (With the one-threshold count the balanced cap at j = 2 is 1/32 < 1/30. The
    unbalanced branch still gives C(e) − e = 1/8 + 3e/4 < 1/6 at e = 1/30; Opus critic
    item 11.)
- **Ladder + detector, any j > 5/6.** Use the modified schedule of §3: ratio ρ < 6j/5 and
  H_max = n^A·2^{n/(60j)}. The rung contribution H_{i+1}^{5/3}λ_i^{−j} then stays polylog.
- **Upshot for sieve2d.** For the average-cost theorem, a better tail buys only polylog
  factors. What it buys substantively:
  - D′'s fixed-cost trade-off (cost N^{1/6+δ} with exceptions N^{−2jδ·3/5});
  - the single-run statement;
  - the worst-case program.

## 7. What is and is not proved

- **Proved.**
  - Proposition U and U′.
  - Lemmas R, D, S and C.
  - Theorem AV (i)–(iv).
  - The ε audit of §4: 1/52 as written, 1/46 with T re-chosen, 1/40 with Lemma D, 1/22
    with the detector.
  - The single-fallback optimum 2/11, within the proved upper-bound model.
  - The no-detector ladder infima: 19/110 (one-threshold) and 17/100 (dyadic). These are
    within the proved bounds and are not lower bounds on the algorithm.
  - The conditional table of §6.
- **Exact computation.** All LP values (`src/avgtime/exponents.py`, a two-phase simplex in
  Fractions; tests in `tests/test_avgtime.py`).
- **Not proved.** Any worst-case improvement; ratios q/p > 2; any lower bound on average
  deterministic factoring; tails better than λ^{−1}.
- **Dependencies.** Theorem AV inherits D′'s correctness (three referee passes, Lean
  formalisation of Lemmas 1–8 and 10 by `leanform`). It adds only the free-H uniformity
  check, the elementary Lemmas S and C, and the r16 refinement R (polylog only).

## 8. Task item 4: prior art (own WebSearch + GPT-6.1-sol `codex/avgtime_2.out`)

- **The question.** A deterministic, always-complete algorithm with proved **average** cost
  N^c, c < 1/5, over balanced semiprimes, or an almost-all statement below 1/5.
- **What was searched.** Own WebSearch: three queries on average-case and almost-all
  deterministic factoring and the HH statement. GPT-6.1-sol checked 20 items.
- **Not found.** No such theorem was found in:
  - McKee–Pinch 1996: §5 is heuristic N^{1/4+ε};
  - McKee 1999: heuristic;
  - Lehman 1974;
  - Hittmeir 2018, 2021, 2023;
  - Harvey 2021; Harvey–Hittmeir 2022; Harvey–Hittmeir 2026 (arXiv 2601.11131);
  - Gao–Feng–Hu–Pan (arXiv 2512.19076): balanced worst case N^{1/5}, log factor only;
  - Stănică–Mulder–Hittmeir (arXiv 2606.13018);
  - Gower–Wagstaff 2008 (SQUFOF): average analysis heuristic, N^{1/4};
  - Pollard rho: heuristic;
  - class-group methods;
  - fixed-curve ECM.
- **Different scope.**
  - Knuth–Trabb Pardo 1976 and Hafner–McCurley 1989 give genuine input-distribution
    analyses, but for trial-division-type algorithms on all integers.
  - Źrałek 2010 gives deterministic p − 1 under smoothness or totient hypotheses.
  - Crişan–Erban 2021: almost all *unrestricted* semiprimes have a small factor.
- **Positive proportion only.** GPT-6.1-sol cites a posted report, Ward–Maestro (Sept 2026,
  "Floors for deterministic integer factoring", Part 5, Prop. 35.3). It proves base-2 p−1
  success on a positive proportion of pairs at exponents near 1/(4√e) ≈ 0.1516, and
  explicitly does not give an average for an always-completing algorithm. *We have not read
  it ourselves (UNVERIFIED).* A positive proportion does not give an average below 1/5,
  since the failure fraction is constant.
- **Carella 1308.2891.** Claims worst-case 1/6. The proof defect (polynomial content
  divisible by B) was found independently by r13, r16 and avgtime_2.
- **p − 1 obstruction (our item 6).** Confirmed for smoothness-based analyses: at B = N^{1/6}
  the success probability is a constant, so any fallback costs a constant fraction of
  N^{1/5}. GPT-6.1-sol notes that no universal impossibility theorem exists for every
  stage-2 design. Compare pm1bsgs: E = N^{1/4−o(1)}.
- **Verdict.** Not found by us. The ingredients are standard:
  - D′;
  - restarting with growing parameters ("doubling");
  - BSGS order search.

  The combination's content is the observation that D′'s proved tail (H^{−2}) beats its
  cost growth (H^{5/3}). Novelty judgments are left to the orchestrator.

## 9. Review

Both reviewers were read-only. Every requested fix has been applied.

- **GPT-6.1-sol adversarial referee** (`codex/avgtime_1.out`, xhigh). **Verdict: ACCEPT
  WITH FIXES.**
  - **Confirmed CORRECT:**
    - the ε audit and the 2/11 corollary;
    - Proposition U's uniformity, including Lemma 10 for general T;
    - Lemma R, with the μ²/φ identity now written out;
    - Lemmas S and C;
    - Theorem AV(iii): the reach-set inclusions, with no independence needed, every
      exponent, and the HH term;
    - (iv);
    - prior art. It also checked Ward–Maestro Prop. 35.3 itself: positive proportion only.
  - **Fixes applied:**
    - (1) ERROR: the C(e) envelope was missing the branch 9e/2 for e ≥ 1/14, where y ≤ 1/4
      binds. Added and tested. The value 19/110 is unaffected, since e* = 3/110 lies in
      the valid branch. Reachability is now stated for z > 1/6, by induction.
    - (2) "Valid exactly" became "certified by these estimates".
    - (3) n₀ now also covers H_1 ≥ H₀ and H_max ≤ 2^{n/50}.
    - (7) (iv) is now applied per shell, and its scope is stated.
    - (8) The m ≤ m̄ constraint was added to the LP; all values are unchanged. The
      optimum is called an optimum of the upper-bound model.
    - (9) GAP: j > 5/6 needs a modified schedule. Codex's schedule was added to §3 and §6.
    - (10) GAP: an explicit "n < n₀ ⇒ HH" step was added.
- **Independent Opus critic** (read-only subagent). Its own LP encoding (sympy), a
  brute-force detector check on 9000 pairs with 15 values of T (0 misses, 0 spurious
  outputs), and hand re-derivations. **Verdict: ACCEPT WITH FIXES.**
  - **Confirmed CORRECT:** Proposition U, Lemmas R, S and C, Theorem AV(iii)/(iv), the
    single-fallback 2/11, and the LP encoding (for the stated model).
  - **Fixes applied:**
    - (6) ERROR (overclaim): "the no-detector ladder stalls at 19/110" was wrong. D′'s
      success proof never uses ord_q(2) > T. Summing Lemma 10 dyadically over ord_q(2)
      (now **Lemma D**) gives the no-detector infimum 17/100, verified by the LP and in
      closed form. The wording now says "within the one-threshold accounting".
    - (7) The unchanged D′ is certified for ε < 1/40 (Lemma D). Added to §0 and §4.
    - (9) Proposition U′ (free m, y) was stated as a separate lemma.
    - (10) Same as codex (1).
    - (11) The j ≥ 2 no-detector limit 1/6 is now proved, no longer "experimentally
      supported".
    - (12) The j = 5/6 boundary was corrected.
    - (13) The §11 parenthetical was replaced.
    - (14) The even-N step now runs HH on N/2.
    - (15) "Three lemmas"; the log log N factor was added.
    - (2) Montgomery's composite-moduli large sieve is now cited explicitly in Lemma R.

## 10. Reproduction

```
.venv/bin/python -m pytest -q tests/test_avgtime.py          # 21 tests, ~5 s
.venv/bin/python -c "from src.avgtime.exponents import *; print(dprime_admissible(), dprime_dyadic_sup()); \
  print([ladder(k,1,False)[0] for k in (1,2,3)]); print(single_fallback(2))"
```

## 11. Candidate card

- **Mechanism.** A ladder of D′ runs, DCL_{n,H_i} with H_{i+1} = ⌈H_i^{6/5}⌉, behind a
  BSGS small-order detector (T = N^{3/10}), with an HH fallback at the top.
- **Outcome.**
  - Proved: average cost N^{1/6}polylog on S_n; worst case = HH; complete; D′ certified
    for ε < 1/40 (Lemma D).
  - Refuted for this model: none.
  - Superseded: the 2/11 single-fallback corollary (correct but not optimal).
- **Scoped obstruction (within the proved bounds).**
  - Without small-order handling, the E_c count is limited by the primes q with
    ord_q(2) ≈ √P. That caps a balanced rung at failure exponent 1/40, and the ladder's
    infimum is 17/100.
  - A better coverage tail cannot lower the average exponent below 1/6 within D′, since
    every run costs ≥ N^{1/6}.
- **Open question.** Prove a matching lower bound in some model: is N^{1/6} the right
  average-case exponent for deterministic factoring of balanced semiprimes? Or does some
  input-adaptive method get below 1/6 on average? (Cf. below16; multimod Thm A covers only
  full-candidate covering.)
