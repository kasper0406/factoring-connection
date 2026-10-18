# unbalanced: Theorem D′ beyond q ≤ 2p — ratio bands, the ratio gap, and a curvature barrier

Agent `unbalanced`, 2026-10-18 (restarted after the first machine was lost). Pre-registration:
`notes/unbalanced_prereg.md` (commit `b3beb0d`, before any pre-registered run). Code:
`src/unbalanced/` (operational `params.py`, `algorithm.py`; harness `harness.py`,
`experiments.py`, `analyze.py`; diagnostic `obstruction.py`). Tests:
`tests/test_unbalanced.py`, `tests/test_unbalanced_obstruction.py`. Raw data:
`experiments/raw/unbalanced_*.jsonl`. Summary: `experiments/analysis/unbalanced_summary.json`.
Inputs: notes/dcltheorem_results.md (Theorem D′, authoritative), notes/multimod_results.md
§3.1–3.2 and §3.9, GPT-6.1-sol `codex/unbalanced_1.out` (blind) and `codex/r16_strategy.out` §2.

## 0. Verdict (plain language)

**D′ extends cleanly to every ratio q/p ≤ N^{o(1)}, but not to all semiprimes. There is a
real ratio gap. For ratios N^{θ} with 0 < θ < 1/3, our best cost is above N^{1/6}. The
worst case is θ = 1/5 (p ≈ N^{2/5}), where we only reach N^{1/5}, the known Harvey /
Harvey–Hittmeir exponent. We prove that no scheme in a precise and fairly broad model can do
better at full coverage.**

- **Theorem U (proved; regime I).** D′ with ratio-band weight rectangles
  A ≤ H/√R, B ≤ H√R (R = 2^j ≤ H²/16) has the same exponent. It costs N^{1/6+5ε/6}·polylog
  and fails on O(N^{−ε}) of the pairs with q/p ≤ N^{ε}. Every window span stays
  ≤ (1+√2)H√N, so the primes, γ, the four blocks and the baby list are D′'s. Only the target
  count grows, by the number of bands, O(n). The large-sieve, clean-band and FULL-count
  lemmas carry over with H replaced by H² in one gcd bound.
- **Theorem U_θ (proved; regime II).** For H²/16 < R ≤ N^{θ}, θ ≤ 1/5 − 3ε, use weights
  (1, B) with 2 ≤ B ≤ H² (a construction suggested by GPT-6.1-sol, simpler than our
  pre-registered one) and a per-band prime scale Y_j ≍ (NR_j)^{1/6}H^{−2/3}. A band costs
  H^{4/3}(NR)^{1/6}, so the exponent is (1+θ)/6 + 2ε/3. Our pre-registered fine-band
  design W2 (weights (1, B) with B ≈ R) has the same exponent (Variant U_θ′, proved).
- **Profile (proved, with Pollard–Strassen as prior art).** At ratio N^{θ} the best exponent
  we have is e(θ) = min((1+θ)/6, (1−θ)/4) (+O(ε) on the DCL side). So:
  - it exceeds 1/6 for every 0 < θ < 1/3;
  - it is below 1/5 for every θ ≠ 1/5;
  - its maximum is 1/5, at θ = 1/5.
- **Proposition O′ (proved in its model; the strongest negative result here).** Take any
  one-dimensional collision search: any number of baby lists with *arbitrary* exponent sets,
  giants with non-constant Lehman-linear targets αp + βq + c₀ ((α, β) ≠ (0, 0)), any lock
  moduli (prime, composite, CRT),
  either side extracted, multipliers, rational weights. Covering a fraction f of the integer
  candidates of ratio band R costs ≥ c·min(√(fp), f^{2/3}q^{1/3}).
  - At R = N^{1/5} this is ≥ c f^{2/3}N^{1/5}. **In this model, at full (constant-fraction)
    coverage, the 1/5 peak cannot be removed.**
  - The bound is attained inside the model by explicit searches, so it is sharp there. The
    witnesses of §5.2 use relaxed (non-lock) coverage. Real
    algorithms reach the same profile in the almost-all sense: Theorem U_θ on the q^{1/3}
    branch and Pollard–Strassen on the √p branch (§5.2).
  - Its scope is FULL coverage of the candidates, not average-case success. An almost-all
    scheme could beat it only if its covered sets favoured the true factor beyond their
    density. Non-linear targets, lattice methods and multi-dimensional searches are also
    outside it (§5.3).
- **Populations.**
  - For the bounded-ratio population S_n, and per ratio band up to N^{o(1)}: N^{1/6+o(1)}.
  - Uniform over semiprimes with both primes ≥ N^{1/3}: the exponent θ has a continuous law
    on [0, 1/3]. At cost N^{1/6+η} our methods (profile e(θ)) handle only ≈ 30η of them,
    and they give **no** almost-all exponent below 1/5. Exponent 1/5 − η′ still handles all
    but ≈ 30η′.
  - Counted by N ≤ X over all semiprimes, almost all have p < X^{1/3}. Pollard–Strassen
    already gives N^{1/6} for 1 − o(1) of them, so that measure is trivial. But it gives no
    power-saving exceptional set, since a ≍ 1/log log N fraction has θ near 1/5.
- **Beating R^{1/6}: no escape found.** We tried composite locks, q-side extraction,
  multipliers kN, rational or anisotropic weights, finer sub-band lists, composite moduli for
  multiplicity, Lucas/cyclotomic handles, Coppersmith windows, PS on lock classes, and the
  McKee–Pinch pigeonhole. Each is either inside Proposition O′'s model, or costs more (§6).
- **Experiments (pre-registered, run exactly as frozen; §7).** All six frozen decision rules
  hold.
  - 0 wrong outputs and 0 lemma violations, on 6625 covering targets.
  - Regime I (W1) succeeds 24/24 in every cell, uniformly over q/p up to 2^7, at a constant
    ≈ 7× D′ overhead.
  - Unchanged D′ loses ≈ 95% of its structural coverage at q/p ≈ 2^6.
  - Regime II's (W2) cost grows with slope 0.152–0.158 in log₂R at n = 80/96, against
    n-only 0.14 and asymptotic 1/6. The R^{1/6} law is seen, not refuted.
  - Lists match their frozen n-only predictions within 3%.
- **Extraction remark (GPT-6.1-sol, labelled conditional).** Harvey–Hittmeir 2026 gives a
  base of order > m mod N deterministically. With seven blocks and exact-equality recovery,
  that would remove D′'s order and FULL exceptional sets entirely, leaving only coverage
  (§3.4).

## 1. Populations and notation

- n is the bit length, N ∈ [2^{n−1}, 2^n).
- **Ratio band j:** S_{n,j} = {(p, q) primes : 2^j p < q ≤ 2^{j+1}p, 2^{n−1} ≤ pq < 2^n}.
  Write R = R_j = 2^j. Then p ≍ P_j := √(N/R) and q ≍ Q_j := √(NR) (within factors √2).
- S_n^{θ} = ∪_{j ≤ θn} S_{n,j} is the set of pairs with q/p ≤ N^{θ} roughly.
  - j = 0 is D′'s population S_n, up to the boundary q = 2p.
- **Counts (PNT, uniform for j ≤ n/3):** |S_{n,j}| ≥ c₀ 2^n/n². Choose p ∈ [x, 1.1x],
  q ∈ [y, 1.1y] with xy = 0.6·2^n and y/x = 1.5·2^j. Both x, y ≥ 2^{n/3}, so each
  progression has ≫ x/n, y/n primes.
- Exceptions are **counted by pairs**, per band. Because |E_j| ≤ C2^{−εn}|S_{n,j}| for
  every band, the union over bands is also an O(2^{−εn}) fraction of S_n^{θ}.

## 2. Theorem U (regime I: R ≤ H²/16)

### 2.1 Algorithm DCL-wide-I_{n,ε} (= `algorithm.dclw` with `params.theorem_design_U`)

All parameters are D′'s (dcltheorem §2.1):
- H = ⌈n2^{εn/2}⌉;
- Y = ⌈2^{n/6}H^{−1/3}⌉;
- four blocks 𝒫_k = {ℓ ≡ 3 (4) prime, 2^kY < ℓ ≤ 2^{k+1}Y};
- base 2, L_k = ∏ℓ².

**Change.** Let J be the largest j with 2^j ≤ H²/16. For each j ≤ J let

  V_j = {(A, B) : 1 ≤ A ≤ a_j := ⌊H2^{−j/2}⌋, 1 ≤ B ≤ b_j := ⌊H2^{j/2}⌋, gcd(A,B) = 1, A ≠ B}.

- The instance (A, B, j) carries the ratio interval (2^j, 2^{j+1}].
- Its window is computed over the integer p-range x ∈ [x_lo, x_hi] = {x : 2^j < N/x² ≤ 2^{j+1}}
  (`ratio_range`, exact integers).
- Lo = ⌊2√(ABN)⌋ and Hi = max_{x∈{x_lo,x_hi}} ⌈Bx + AN/x⌉. The rest of each block is D′
  verbatim, with one shared baby list per block for all instances of all bands:
  - targets for every ℓ and every instance;
  - roots by a^{(ℓ+1)/4};
  - c, g_lo, W;
  - giants u_ℓ^{B+AN−c}γ^{−(g_lo+jm)};
  - m = m̄_k, or the endpoint-competitive grid;
  - product/remainder trees, gcd, validated outputs.
- With J = 0 the target list is exactly D′'s (test `test_band0_reproduces_dprime_targets`).

### 2.2 Lemmas (changes to D′'s proof only)

Fix (p, q) ∈ S_{n,j}, j ≤ J. Write R = 2^j.

- **U1 (span).** For x in band j's p-range,
  F(x) = Bx + AN/x ≤ b_j√(N/R) + a_j√(2RN) ≤ (1+√2)H√N, and Lo ≥ 0.
  So Hi − Lo ≤ (1+√2)H√N + 2 ≤ span̄. This is D′'s span̄ (5/2)H(2^{n/2}+1) + 2H + 3.
  Hence m̄_k is unchanged. (GPT-6.1-sol: with the exact interval minimum the real span is
  ≤ (√2−1)H√N; we use the cruder bound, as the code does.)
- **U2 (weights).**
  - a_j b_j ≤ H²;
  - |AB′ − A′B| ≤ a_j b_j ≤ H² < Y < ℓ;
  - |V_j| ≥ c₂H². For a_j ≥ 4: a_j ≥ H/(2√R) and b_j ≥ H√R/2. Coprime pairs in [1,a]×[1,b]
    number ≥ b·Σ_{A≤a}φ(A)/A − O(a log a) ≥ c·ab, and A = B removes ≤ a.
  - |A − B| ≤ b_j ≤ H√R ≤ H².
- **Lemma 1 (lock), Lemmas 2–5, Lemma 7 (FULL ⇔ d | (A−B)(p−1) − ℓ²δ) hold verbatim.**
  - They only use ℓ ∤ ABN, ℓ < p and the window containing F(p).
  - ℓ ≤ 16Y ≪ √(N/(2R)) ≤ p, and ℓ > Y > H² ≥ A, B.
  - The window contains F(p) because N/p² = q/p ∈ (2^j, 2^{j+1}].
- **Lemma 6 (clean block)** needs Y⁴ > q. Here q ≤ √(2^{n}·2^{J+1}) ≤ H·2^{n/2}, and
  Y⁴ ≥ 2^{2n/3}H^{−4/3} > H2^{n/2} for ε < 1/14 and n ≥ n₀.
- **Lemma 9 (coverage, large sieve) holds uniformly in j.**
  - Fix a prime q and let I_q = {p ∈ Z : q/2^{j+1} ≤ p < q/2^j}, of length ≤ q/(2R) + 1.
  - By U2, the |V_j| classes AqB^{−1} mod ℓ are distinct. So the uncovered
    U_k(q) ⊂ I_q avoids |V_j| classes mod every ℓ ∈ 𝒫_k.
  - Montgomery–Vaughan (prime terms only, exactly as in D′) gives
    |U_k(q)| ≤ (|I_q| + 1 + 4^{k+1}Y²)·n/(c₁H²).
- **Lemma 10 (FULL count)** holds with g′ = gcd(A−B, d) ≤ H² (instead of ≤ H):
  |F_k(q)| ≤ 2|V_j||𝒫_k| m̄_k (|I_q|H²/(2^kY·T) + 1).

### 2.3 Theorem U

**Theorem U (proved).** Fix ε ∈ (0, 1/60]. For n ≥ n₀(ε), DCL-wide-I_{n,ε} has these
properties.
- It outputs a proper divisor or FAIL on every input.
- It costs O(n^{c}·2^{(1/6+5ε/6)n}) bit operations on every n-bit input.
- It fails on at most C·2^{−εn}|S_{n,j}| pairs of each band j ≤ J, where 2^J ≥ H²/32 ≥ 2^{εn}.

In particular it fails on O(N^{−ε}) of the pairs with q/p ≤ N^{ε}.

*Proof.*
- **(i) Correctness.** Identical to D′.
- **(ii) Cost.** Per block, the targets number ≤ 2Σ_j|V_j||𝒫_k| ≤ 2(J+1)H²·2^kY with
  J + 1 ≤ 2log₂H + 1 = O(n). The list length is m̄_k + T′ with T′ ≤ T + Σ W/m. Every other
  row of D′'s cost table scales by at most J + 1. So the cost is O(n^{O(1)}(HP/Y² + nH²Y)),
  i.e. O(n^{O(1)}H^{5/3}2^{n/6}).
- **(iii) Exceptions, per band j**, with T = ⌊2^{(1/4−ε)n}⌋. The failing pairs lie in
  E_a ∪ E_b ∪ E_c, as in D′. For each, sum over the primes q of the band, at most
  C√(NR)/n of them, and divide by |S_{n,j}| ≥ c₀2^n/n².
  - **E_a (ord_q(2) ≤ T).** By Lemma 8, ≤ T(T+1)/(n−1) bad q, each with ≤ |I_q| + 1 ≤
    C√(N/R) partners. Fraction ≤ C n T²/√(NR) ≤ C n 2^{−2εn}.
  - **E_b (uncovered in the clean block).** Σ_q |U| ≤ (C√(NR)/n)(C√(N/R) + 256Y²)·n/(c₁H²)
    ≤ C 2^n/H², because Y² ≤ √(N/R) when R ≤ H². Fraction ≤ C n²/H² ≤ C 2^{−εn}.
  - **E_c (FULL).** Using |𝒫_k| ≤ 2^kY, m̄_k ≤ 3H√N/(4^kY²) and |I_q| ≤ C√(N/R):
    Σ_q|F_k| ≤ (C√(NR)/n)·2H²·2^kY·m̄_k·(C√(N/R)H²/(2^kYT) + 1). Dividing by |S_{n,j}|:

      ≤ C n (H⁵√N/(Y²T) + H⁴/Y) ≤ n^{O(1)}(2^{(−1/12 + 23ε/6)n} + 2^{(−1/6+3ε)n}),

    using √R ≤ H, Y² ≥ 2^{n/3}H^{−2/3}, T ≥ 2^{(1/4−ε)n−1} and H ≤ 2n2^{εn/2}. Compared
    with D′, the extra factor is one H, coming from g′ ≤ H². For ε ≤ 1/60,
    −1/12 + 23ε/6 ≤ −ε − 1/360, so this is O(2^{−εn}) for n ≥ n₀(ε) (see the constant
    note).
- **Union.** The total is O(2^{−εn}) per band. ∎

*Constant note.* D′'s own proof needs ε < 1/52 at this step, and ours needs ε < 1/58. We
state Theorem U for ε ∈ (0, 1/60], where 29ε/6 ≤ 29/360 < 1/12. Then the first term is
n^{O(1)}2^{−(1/12−29ε/6)n}·2^{−εn}, and 1/12 − 29/360 = 1/360 > 0 absorbs the
polynomial factor for n ≥ n₀(ε). So no change of range is needed.

**Corollary U″ (polylog).** With H = ⌈n^{1+a}⌉, as in D″, the algorithm covers ratios up to
n^{2+2a}/32. It costs N^{1/6}(log N)^{O(1+a)} and fails on O((log N)^{−2a}) per band. The
proof is the same.

## 3. Theorem U_θ (regime II: H²/16 < R ≤ N^{θ}, θ ≤ 1/5 − 3ε)

### 3.1 The coarse construction (theorem design; GPT-6.1-sol `unbalanced_1.out` §2)

Put h = H² (D′'s H). For each regime-II band j (J < j ≤ J_θ := ⌊θn⌋):
- **Weights.** V′ = {(1, B) : 2 ≤ B ≤ h}, with the band's ratio interval (2^j, 2^{j+1}].
  - So |V′| = h − 1 ≥ h/2.
  - The classes qB^{−1} mod ℓ are distinct for ℓ > h.
  - |A − B| = B − 1 < h.
- **Prime scale.** Y_j = ⌈(2^{(n+j)/2}/h)^{1/3}⌉. This is D′'s balance Y ≍ (span/|V|)^{1/3}.
  Blocks 𝒫_{j,k} = {ℓ ≡ 3 (4) : 2^kY_j < ℓ ≤ 2^{k+1}Y_j}, k = 0..3, each with its own
  L, γ and baby list.
- **Everything else is D′:** roots, c, windows over the band's p-range, giants, trees, gcds,
  validation.
- Not run as a frozen experiment. It is `params.WideDesign` with
  `insts = [Inst(1, B, 2^j, 1, 2^{j+1}, 1, j) for B in 2..h]`.

**U3 (span; proved; corrected after the Opus critic).**
- Lo = ⌊2√(BN)⌋ ≤ F by AM–GM, and Hi bounds F on the band by convexity.
- Using B ≤ h and R > h/16:
  Hi − Lo ≤ F(x₁) + O(1) ≤ B√(N/R) + √(2NR) + O(1) ≤ (h/R + √2)√(NR) + O(1) ≤ 18Q_j.
- For R > h the minimiser √(N/B) lies right of the band, so F is decreasing there and the
  real span is ≤ (√2−1)√(NR) (GPT-6.1-sol).
- For h/16 < R ≤ h (about 4 bands) F need not be monotone. The critic measured widths up to
  5.95·Q_j at n = 48, h = 64, R = h/16.
- We charge span̄_j := 18Q_j + O(1) with Q_j := 2^{(n+j+1)/2}, and m̄_{j,k} = ⌊span̄_j/(4^kY_j²)⌋ + 2.
  These are public (bit-length-based) bounds; Lemma 10's δ-range uses them.
  GPT-6.1-sol `unbalanced_2` independently gives 20√(2^nR) + 2 and also offers an
  alternative: keep only 2 ≤ B ≤ min(h, R), which preserves monotonicity.

**Lemmas 1–7 hold verbatim** (same checks as §2.2):
- ℓ > Y_j > h ≥ B. Indeed Y_j ≥ (2^{n/2}/h)^{1/3} > h for ε < 1/8.
- ℓ ≤ 16Y_j < x₁. Indeed Y_j ≤ (2^{(n+θn)/2})^{1/3} ≪ 2^{(n−θn)/2}/4 for θ < 1/2.
- Lemma 6: Y_j⁴ > √2 q. Here Y_j⁴ ≥ (Q_j/(2h))^{4/3} > 2Q_j because Q_j^{1/3} > 4h^{4/3}.

**Theorem U_θ (proved).** Fix ε ∈ (0, 1/60] and θ ∈ (0, 1/5 − 3ε]. For n ≥ n₀(ε, θ), the
algorithm DCL-wide_{n,ε,θ} has these properties. It runs regime I (§2) and then the coarse
construction for bands J < j ≤ ⌊θn⌋, on every input.
- (i) It outputs only proper divisors or FAIL.
- (ii) It costs n^{O(1)}·(H^{5/3}2^{n/6} + H^{4/3}2^{(1+θ)n/6}), i.e.
  N^{1/6 + max(5ε/6, θ/6 + 2ε/3)}·(log N)^{O(1)} bit operations, on every n-bit input.
- (iii) It fails on O(2^{−εn})·|S_{n,j}| pairs of each band j ≤ θn, hence on O(N^{−ε}) of
  S_n^{θ}.

*Proof.* Regime-I bands are Theorem U. Fix a regime-II band j; write R = 2^j, Y = Y_j,
P_j = √(N/R), Q_j ≍ √(NR).

- **(ii) Cost per block.**
  - Targets ≤ 2h·2^{k}Y.
  - The list is m̄_{j,k} + T′ ≤ 3·span̄_j/(4^kY²) + 2 + 2h2^{k+1}Y. This is
    ≤ C(Q_j/Y² + hY) = C h^{2/3}Q_j^{1/3}.
  - The exponent tree costs O(Y log Y) mults mod N. The remaining rows of D′'s table are
    Õ(targets + list).
  - Summing the geometric series over j ≤ θn gives C h^{2/3}(2^{n(1+θ)/2})^{1/3}n^{O(1)}.
- **(iii) E_a.** T := ⌊2^{(n+j)/4}2^{−εn}⌋, i.e. ≍ Q_j^{1/2}N^{−ε}.
  - Lemma 8 gives ≤ T(T+1)/(n−1) primes q > 2^{n/2−1} with ord_q(2) ≤ T, each with
    ≤ |I_q| + 1 ≤ 2P_j partners.
  - So |E_a| ≤ 2P_jT(T+1)/(n−1). Dividing by c₀2^n/n² gives ≤ C n T²/Q_j ≤ C n 2^{−2εn}.
- **(iii) E_b (large sieve).** For a prime q of the band let I_q = [q/2^{j+1}, q/2^j) ∩ Z,
  with |I_q| ≤ P_j + 1.
  - In the clean block k(q), the uncovered p ∈ I_q avoid h − 1 distinct classes mod every
    ℓ ∈ 𝒫_{j,k}.
  - So |U| ≤ (|I_q| + 1 + 4^{k+1}Y²)/((h−1)Σ1/ℓ) ≤ C(P_j + Y²)n/h.
  - Now Y² ≤ P_j ⇔ (NR)^{1/3}h^{−2/3} ≤ (N/R)^{1/2} ⇔ R^{5/6} ≤ N^{1/6}h^{2/3}, which holds
    for θ ≤ 1/5.
  - Summing over ≤ CQ_j/n primes q and dividing by c₀2^n/n²: fraction ≤ C n²/h = C2^{−εn}.
- **(iii) E_c (FULL; D′ Lemma 10 with g′ = gcd(B−1, d) ≤ h).**
  - Σ_q |F_k(q)| ≤ (CQ_j/n)·2h·2^kY·m̄_{j,k}·(CP_j h/(2^kY T) + 1).
  - Divide by c₀2^n/n², using Q_jP_j ≍ N, m̄_{j,k} ≤ C Q_j/(4^kY²) and Y³ ≍ Q_j/h:

      fraction ≤ C n (h² m̄/T + hR/Y + hQ_jY/N)
               ≤ C n (h^{8/3}Q_j^{−1/6}2^{εn} + h^{4/3}R^{5/6}N^{−1/6} + h^{2/3}N^{(2θ−1)/3}).

    (The third term is ≤ N^{−1/5+2ε/3}; the critic restored a dropped h^{2/3}.)

  - **First term.** With Q_j ≥ 2^{(n−1)/2} and h ≤ n²2^{εn}, it is
    ≤ n^{O(1)}2^{(−1/12 + 11ε/3)n}. That is O(2^{−εn}) for ε < 1/56.
  - **Second term.** It is ≤ n^{O(1)}2^{4εn/3}R^{5/6}2^{−n/6}, which is O(2^{−εn}) when
    R ≤ 2^{(1/5 − 14ε/5)n − O(log n)}. That holds for θ ≤ 1/5 − 3ε and n ≥ n₀(ε, θ).
    **This is where θ < 1/5 enters the success proof.** For larger R the two-sided "+1" CRT
    term, multiplied by the number of (B, ℓ, δ), is too large.
- **Union** over bands as in §1. ∎

**What fails at θ ≥ 1/5** (corrected after the critic).
- E_c's second term exceeds 1, and the cost (1+θ)/6 overtakes PS's (1−θ)/4. Both turn at
  θ = 1/5.
- The large sieve (E_b) holds until R ≈ N^{1/5}h^{4/5}, i.e. Y_j² ≤ P_j. So E_b is **not**
  the binding constraint.
- With the Harvey–Hittmeir extraction of §3.4, E_c disappears. Then only cost and E_b limit
  θ, and both turn at 1/5 + O(ε).

### 3.2 Variant U_θ′: the pre-registered fine-band design W2 (proved; what oper2 runs)

The pre-registered design is `params.regime2_insts`.
- Weights (1, B) with B ∈ [2^j, 2^{j+1}), each with ratio window (B − Δ, B + Δ].
- **Span (U4).** For N/x² ∈ (B−Δ, B+Δ] and B ≥ 2Δ:
  F(x) − 2√(BN) = (Bx − N/x)²/(Bx + N/x + 2√(BN)) ≤ Δ²x²/(4√(BN)) ≤ Δ²√N/(2B^{3/2}).
  This is `regime2_span_bound`; tested in `test_regime2_span_bound_holds`.
- **Cost.** One baby list per band and block, shared by all R values of B. Targets ≍ RY,
  list ≍ Δ²√N R^{−3/2}Y^{−2}. With Y ≍ Δ^{2/3}N^{1/6}R^{−5/6} the per-band cost is
  Δ^{2/3}(NR)^{1/6}: the same exponent as §3.1 with Δ in place of h.
- **Feasibility.** Y > 2Δ needs R ≤ N^{1/5}(8Δ)^{−2/5}.
- **Coverage (Lemma 9′).**
  - Partition I_q into I_{q,b} = {p : b < q/p ≤ b+1}, b ∈ [2^j, 2^{j+1}).
  - On I_{q,b} every B ∈ [b+1−Δ, b+Δ] ∩ [2^j, 2^{j+1}) is admissible. There are ≥ Δ − 1 of
    them, and at most one is ≡ 0 (mod ℓ) because 2Δ < ℓ. So ω(ℓ) ≥ Δ − 2 distinct classes
    qB^{−1}.
  - Hence |U ∩ I_{q,b}| ≤ (|I_{q,b}| + 1 + Q²)·Cn/Δ. Summing over b costs an extra R·Q²,
    which is ≪ |I_q| because R^{3/2}Y² ≪ √N.
  - Exceptional fraction: C n²/Δ.
- **FULL (Lemma 10′, roles swapped).** Here |A − B| = B − 1 ≍ R, so D′'s per-q count would
  lose a factor R. Instead fix (p, B, ℓ, δ) and count q.
  - In a clean block, FULL means d = ord_q(2) divides X := (B−1)(p−1) + ℓ²δ.
  - For R ≥ 4Δ and n ≥ n₀, 0 < X < 3q, because ℓ²m̄ ≤ 4·span̄ + 8·4^KY² < (R/2)p.
  - So q ≡ 1 (mod e) for some divisor e > T of X, with ℓ ∤ e (clean), q ≡ Bp (mod ℓ), and
    q ∈ (p(B−Δ), p(B+Δ)]. That gives ≤ τ(X)(2Δp/(Tℓ) + 1) values of q.
  - Summing over δ, ℓ, B and p, and dividing by |S_{n,j}|, the fraction is
    ≤ Cn²τ*(m̄Δ/T + m̄RY/√(NR)), with τ* = max_{X<N} τ(X) = N^{O(1/log log N)}.
  - At Y ≍ Δ^{2/3}N^{1/6}R^{−5/6} and Δ ≍ H², this is
    ≤ N^{−1/12+11ε/3+o(1)} + N^{−1/6+4ε/3+o(1)}.
  - No averaging of gcd(B−1, d) is needed. (r16 suggested that averaging; swapping the roles
    avoids it.)
- **Clean block.** Y_j can be as small as ≈ 2Δ = N^{ε+o(1)}. Lemma 6 then needs
  K = ⌈log q/log Y_j⌉ + 1 = O(1/ε) dyadic blocks, which costs a constant factor 2^{O(1/ε)}.
  - The frozen W2 uses K = 4. At the oper2 sizes Y_j⁴ < q for some (n, j); for example
    n = 96, j = 13 has Y = 373, so Y⁴ ≈ 2^{34} against q ≈ 2^{55}.
  - A clean block is then not guaranteed. `certified` checks cleanliness directly.
- **Variant U_θ′ (proved).** Take T_j = ⌊2^{(n+j)/4}2^{−εn}⌋ as in U_θ, and K blocks as above.
  Then for θ ≤ 1/5 − 3ε the statement of Theorem U_θ holds with W2 in
  place of §3.1, with Δ = ⌊h/(64·2^{K})⌋, K = ⌈3/(5ε)⌉ + 2 and a constant C_ε = 2^{O(1/ε)}.
  The coarse design is simpler and has no such constant, so it is the main theorem.
- GPT-6.1-sol `unbalanced_2`: four blocks suffice for W2 asymptotically only when
  θ < (1+16ε)/23. An order-1/ε block count is needed near the top of the range.

### 3.3 Exponent profile with Pollard–Strassen

**Pollard–Strassen (prior art).**
- Deterministic. It finds every prime factor ≤ z of N in O(z^{1/2}(log N)^{O(1)}) bit
  operations (Pollard 1974; Strassen 1976; Bostan–Gaudry–Schost 2007; Harvey 2021,
  Prop. 2.5).
- A pair of ratio ≥ N^{θ} has p ≤ N^{(1−θ)/2}. So PS costs N^{(1−θ)/4} in the **worst
  case**, with no exceptions.

**Corollary P (proved).** For each θ ∈ [0, 1] there is a deterministic algorithm that
factors all but O(N^{−ε}) of the pairs with ratio in (N^{θ}/2, 2N^{θ}] at cost
N^{e(θ)+O(ε)}, where

  e(θ) = min((1+θ)/6, (1−θ)/4).

The DCL branch is valid only for θ ≤ 1/5 − 3ε; PS is valid for any θ.
- e(θ) > 1/6 on (0, 1/3).
- max e = 1/5, at θ = 1/5.
- e(θ) ≤ 1/5 − |θ − 1/5|/6 elsewhere.

**A whole range of ratios** (all q/p ≤ N^{θ₀}) costs the maximum over θ ≤ θ₀:
- for θ₀ ≤ 1/5 − 3ε, (1+θ₀)/6 + O(ε) (Theorem U_θ);
- for θ₀ ≥ 1/5, 1/5 + O(ε). That is no better than Harvey 2021 / Harvey–Hittmeir 2022,
  which give N^{1/5} in the worst case for all N.

**The ratio gap is real for every scheme we have.** For 0 < θ < 1/3 the cost exceeds
N^{1/6} by N^{min(θ/6, (1/3−θ)/4)}.

### 3.4 Remark: deterministic extraction via Harvey–Hittmeir 2026 (proved given their Thm 1.1; not used in the theorems)

**The cited theorem.** Harvey–Hittmeir, arXiv 2601.11131 v2 (5 Jun 2026), Thm 1.1, checked
against the arXiv text: "input integers N ⩾ 3 and D ⩾ 1 with D < N−1 … outputs either
some α ∈ ℤ_N^* with ord_N(α) > D or a nontrivial divisor of N", in
O(D^{1/2}log D (log log D)^{−1/2} log N).

**The construction (GPT-6.1-sol `unbalanced_1.out` §2).**
- Take D = m and 7 disjoint prime blocks with Y⁷ > N.
- (p−1)(q−1) < N has ≤ 6 prime factors > Y. So some block's L is coprime to
  ord_N(α) | λ(N), and then ord_N(γ) = ord_N(α^L) > m.
- Test gcd(γ^i − 1, N) for 1 ≤ i < m. A proper gcd factors N. Suppose ord_N(γ) > m and no gcd is
  proper. Then ord_p(γ), ord_q(γ) ≥ m: otherwise gcd(γ^{ord_p(γ)} − 1, N) is p or N, and
  N would force ord_N(γ) ≤ ord_p(γ) < m.
- Find exact giant–baby equalities mod N by sorting, and test each directly: solve
  Bz² − Tz + AN = 0 for the proposed T = c + ℓ²(g_lo + jm + i).
  - For the true giant the index is unique mod p (ord_p(γ) ≥ m), so the test succeeds.
  - A failed test proves the giant is not the true one.
- Otherwise, f(G) ≡ 0 (mod N) with no exact equality means there is an i with p | G − γ^i
  and q ∤ G − γ^i. Scanning that one giant (O(m)) then gives a proper factor.

**Consequence (proved, given Thm 1.1).**
- On the clean block, every covering true giant yields a factor.
- E_a and E_c disappear. The exceptional set of Theorems D′, U and U_θ reduces to E_b
  (coverage) alone.
- The extra cost is O(√m·polylog) for α, three more blocks, and O((m+G)log) for the sort.
- What is derived is ord_p(γ), ord_q(γ) ≥ m, and that is all that is used. Babies are then
  distinct mod p, so each giant needs at most one exact-equality test.
- The good block depends on (p−1)(q−1), not on q alone. So E_b must be union-bounded over
  all 7 blocks, at a constant factor (critic's finding 10).
- The θ restriction then comes only from E_b (θ ≤ 1/5) and from cost.

**Operational completion (GPT-6.1-sol `unbalanced_2`).**
- Use D ≥ max_k m_k, so one base α serves every block.
- If some baby gcd equals N, skip the block. This keeps repeated babies on bad blocks from
  creating quadratically many exact matches.
- Test and remove every exact match.
- On a product gcd N, scan that giant's m differences.

We did not implement this, and it is not part of the frozen experiments.

## 4. Which populations make the result meaningful

Profile: e(θ) = min((1+θ)/6, (1−θ)/4), with θ = log(q/p)/log N. The functions
`obstruction.e_profile`, `frac_theta_below` and `frac_covered_at` compute it.

1. **Bounded ratio (S_n, D′) and any fixed band up to N^{o(1)}.** Exponent 1/6 + o(1) with
   N^{−ε}-type exceptions per band (Theorems D′ and U). This is the meaningful
   almost-all-1/6 statement.
2. **Uniform over semiprimes with both primes ≥ N^{1/3}**, counted by pairs (equivalently by
   N ≤ X). By Mertens, the density of p ≈ N^{a} is ∝ da/(a(1−a)). With a = (1−θ)/2, θ has
   density ∝ 1/(1−θ²) on [0, 1/3], so Pr[θ ≤ t] = log₂((1+t)/(1−t)) (proved asymptotically;
   `frac_theta_below`).
   - At cost N^{1/6+η}, only θ ≤ 6η (DCL) and θ ≥ 1/3 − 4η (PS) are handled. That is a
     fraction ≈ 21η/ln 2 ≈ 30η: 0.151 at η = 0.005 and 0.301 at η = 0.01. **So our methods give no
     almost-all N^{1/6+o(1)} statement for this population.** That no method can is not
     claimed.
   - Exponent 1/5 − η′ fails only on θ ∈ (1/5 − 6η′, 1/5 + 4η′), a fraction ≈ 30η′. For our profile,
     every exponent below 1/5 therefore holds for a positive proportion < 1, and none holds
     for 1 − o(1).
   - The median θ is (√2−1)/(√2+1) = 0.172, where e = 0.195.
3. **Log-uniform ratio on [1, N^{1/3}]** (θ uniform). Same conclusion, with constants
   10η·3 = 30η.
4. **All semiprimes counted by N ≤ X.**
   - Almost all have p < X^{δ} for every fixed δ: the fraction with p > X^{δ} is
     ≈ ln((1−δ)/δ)/ln ln X (critic's correction). So PS gives exponent δ/2 with o(1) exceptions, and the
     question is trivial.
   - With **power-saving** exceptions (N^{−ε}) nothing below 1/5 follows. The pairs with
     θ ∈ [1/5 ± η] are a ≍ η/log log X ≫ X^{−ε} fraction.
   - Counted this way, "all but N^{−ε} of all semiprimes at N^{1/6+O(ε)}" is **not achieved
     by our methods**. In the full-cover sense of §5 (covering all but an N^{−ε} fraction of
     the candidates for each N), it is **impossible for
     every scheme in Proposition O′'s model** (full-cover sense), and not achieved by ours.

## 5. The curvature obstruction (Proposition O′)

### 5.1 Model 𝔏 (one-dimensional Lehman-linear collision searches)

Fix N and a ratio band R₀. The **candidates** are the integers p ∈ [x₁, x₂] with
R₀ ≤ N/p² ≤ 2R₀, so x₂ = ⌊√(N/R₀)⌋ and x₂ − x₁ ≥ x₂/4. A search consists of:
- baby lists Λ_1, …, Λ_s, where Λ_i = {γ_i^e : e ∈ S_i} for an **arbitrary** set S_i ⊂ Z
  with |S_i| = m_i;
- giants G, each compared against one list i(G), and each carrying a target
  τ = (α, β, c₀, M, r, g₀) with integers α, β, c₀, g₀, (α, β) ≠ (0, 0), and modulus M ≥ 1.
  A target with α = β = 0 does not depend on p. Its collision would be mod N, two-sided, and
  could never be extracted, so it is excluded (the critic's finding 5).

The target's real-valued index on its class is

  g_τ(p) = (αp + βN/p + c₀)/M²,   p ≡ r (mod M).

G **covers** a candidate p in its class iff ⌊g_τ(p)⌋ ∈ g₀ + S_{i(G)}.
- For the true factor g_τ is an integer (a lock), so this is exactly "G collides at p if p
  is the factor". For other candidates it is a relaxation, which only enlarges coverage.
- For β = 0, refine the class to the sub-class on which g_τ is integral. A true factor
  outside that sub-class cannot collide, so nothing is lost.

The **cost** is X = Σ_i m_i + #giants.

Nothing is assumed about how giants are computed, the group or base, the modulus (prime,
composite, CRT products), the residue class, or which factor is extracted.
- A q-side target αq + βp + c₀ is, in the p-parametrisation, βp + αN/p + c₀: the same form.
- Multipliers kN give integer weights.
- **Rational normalisations** (critic's finding 6). A target may be (αp + βN/p + c₀)/D with
  D | βM², up to units. Integrality on the class forces this, since the x² coefficient is
  βN·M²·r^{−3}. Then g″ = 2βNM²/(Dp³). We assume D ≤ |β|M², so κ ≥ 2N/x₂³, and (a)–(c)
  hold verbatim with |β| replaced by |β|M²/D ≥ 1.

### 5.2 Statement and proof

**Proposition O′ (proved).**
- (a) If β ≠ 0, a giant with |S| = m covers at most 2m + 8√(m x₂³/(|β|N)) ≤ 2m + 8√(m x₂/R₀)
  candidates.
- (b) If β = 0, it covers at most m.
- (c) Hence if a search in 𝔏 covers at least a fraction f of the candidates, then

    X ≥ min( √(f x₂)/4, (f/64)^{2/3}(x₂R₀)^{1/3} ) ≍ min(√(fp), f^{2/3}q^{1/3}),

  where x₂ ≍ p and x₂R₀ ≍ q, both up to a factor √2.

*Proof of (a).*
- Put p = r + Mx. Then g(x) := g_τ(r+Mx) satisfies g''(x) = 2βN/p³, a chain-rule
  cancellation of M². WLOG β > 0 (otherwise negate g; the half-open floor intervals flip,
  and nothing changes). So g'' ≥ κ := 2βN/x₂³ > 0 on the band.
- Let x₀ minimise g on the band (possibly an endpoint) and m₀ = g(x₀). Consider the right
  part x ≥ x₀; the left part is symmetric.
  - There g′(x) ≥ κ(x − x₀) ≥ 0.
  - Since g′ is increasing, g(x) − m₀ ≤ (x − x₀)g′(x) ≤ g′(x)²/κ, so
    **g′(x) ≥ √(κ(g(x) − m₀))**.
- For an integer v put u = v − m₀ and P_v = {x ≥ x₀ : v ≤ g(x) < v + 1}, an interval.
  - It is empty if u ≤ −1.
  - Its length is ≤ √(2(u+1)/κ), from g − m₀ ≥ κ(x−x₀)²/2.
  - For u > 0 its length is also ≤ 1/√(κu), because g′ ≥ √(κu) on P_v.
  - It contains ≤ length + 1 integers.
- The m values v ∈ g₀ + S give distinct u with integer spacing. At most one lies in (−1, 0]
  (length ≤ √(2/κ)), at most one in (0, 1) (length < 2/√κ), and the i-th value ≥ 1 is ≥ i.
- So the right part covers ≤ m + (√2 + 2 + Σ_{i≤m} i^{−1/2})/√κ ≤ m + 5.5√(m/κ) integers.
- Doubling for the two parts, with √(m/κ) = √(m x₂³/(2βN)):
  2m + 11/√2·√(m x₂³/(βN)) ≤ 2m + 8√(m x₂³/(βN)). Finally x₂³/N ≤ x₂/R₀. ∎(a)

*Proof of (b).*
- On the refined class, g(x) = g(0) + kx with k = αM′/M² a nonzero integer, where M′ is the
  refined modulus. So each value has at most one preimage. ∎(b)

*Proof of (c).*
- Let G_i be the number of giants on list i. The coverage is
  ≤ Σ_i G_i(2m_i + 8√(m_i x₂/R₀)) ≤ 2X² + 8√(x₂/R₀)X^{3/2}, since m_i, Σ G_i ≤ X.
- It must be ≥ f x₂/4. If both X² < f x₂/16 and X^{3/2} < f√(x₂R₀)/64, the total is
  < f x₂/4, a contradiction. ∎

**Corollary O″ (proved, in 𝔏).**
- At ratio R₀ = N^{θ}, full or constant-fraction coverage costs
  ≥ c·N^{min((1−θ)/4, (1+θ)/6)}.
- That is ≥ c N^{1/5} at θ = 1/5, with any number of moduli. Compare:
  - cycloharvey's LB is N^{1/5} with a single modulus;
  - multimod Thm A is N^{1/6} for balanced pairs with many moduli. Prop O′ is its
    ratio-dependent form, and it now allows **arbitrary baby sets**.
- **Sharpness inside 𝔏** (witnesses supplied by the critic; proved by direct count):
  - (i) **√p branch.** β = 0, one list of m = ⌈√x₂⌉ consecutive exponents, giants with
    α = 1, M = 1, each covering m consecutive candidates. Cost ≍ √p.
  - (ii) **q^{1/3} branch.** For each unit ratio sub-band b ∈ [R₀, 2R₀) take α = b, β = 1, a
    single modulus M ≍ N^{1/6}R₀^{−5/6}, and every residue class r mod M. The cost is
    m + R₀M with m ≍ √N R₀^{−3/2}/M², giving full coverage at ≍ (NR₀)^{1/6} ≍ q^{1/3}.
  - Witness (ii) uses non-lock targets, i.e. classes whose giants cannot actually be
    computed from N. So optimality inside 𝔏 is about **relaxed** coverage. The lower bound
    is unaffected (it only becomes stronger).
- **Real algorithms reach the profile only in the almost-all sense.**
  - Theorem U_θ makes the true factor covered for all but N^{−ε} of the pairs at
    q^{1/3}N^{O(ε)}.
  - PS gives √p in the worst case, but it is not an 𝔏-search.
  - That U_θ covers a constant fraction of the *candidates* for a fixed N is **not proved**.
    It would need control of ω(ℓ) = 2#{B ≤ h : χ_ℓ(B) = χ_ℓ(N)}, i.e. short character sums
    with h < ℓ^{1/4}, below Burgess. It is heuristic: Poisson, as observed in D′ §8.3.
- **Conclusion.** No 𝔏-search beats e(θ) at full coverage, and our real algorithms achieve
  e(θ) + O(ε) for almost all pairs of each band.

Tests: `test_prop_o_bound_random` (interval windows, best g₀), `test_prop_o_prime_arbitrary_sets`
(best m-set for one giant = the m most popular floor values). The bounds hold on 700
random targets, with the worst ratios 0.50 (interval) and > 0.1 (sets). For an interval at
the band edge, a single side is in range, so the factor-2 slack is expected.

### 5.3 Scope (what Proposition O′ does NOT say)

- **Full cover, not average case.** It bounds how many candidates a search can cover for
  one N. An almost-all algorithm needs only the *true* factor covered for most N.
  - It could beat the bound only if its covered sets C(N) contained the true factor with
    probability well above |C(N)|/#candidates.
  - All lock-class designs have the density property: their classes are N-determined, and
    Lemma 9 is a density statement. But no general theorem excludes correlated designs.
  - Average-case lower bounds: **not proved** (GPT-6.1-sol r16 and unbalanced_1 make the
    same point).
- **Accidental collisions** (mod a reduced group order) are not counted. Heuristically they
  occur at rate ≈ X²/ord_p(γ) ≤ X²/p^{1−o(1)}, which matters only when X ≳ √p, the other
  branch.
- **Outside the model:**
  - targets non-linear in (p, q): exponent polynomials of degree ≥ 2 in N, e.g. N² ≡ q²
    (mod p−1), whose curvature in p is larger (multimod Thm A′ handles the balanced
    degree ≥ 2 case);
  - multi-dimensional collision searches;
  - lattice / Coppersmith decoding;
  - polynomial-evaluation methods (PS);
  - class groups;
  - subexponential methods (heuristic or rigorous-randomised).
- **Cost model.** X counts babies plus giants, i.e. window comparisons. A computation that
  processed a structured family of windows implicitly is not excluded.

## 6. Attempts to beat the R^{1/6} penalty (candidate cards, condensed)

All are evaluated at intermediate R = N^{θ}, 0 < θ < 1/3.

| # | Modification | Outcome | Label |
|---|---|---|---|
| 1 | Anisotropic weights (A, B) with B/A ≈ R and A = a > 1 (rational weights, multipliers kN) | span ÷ a but instances × a. Cost Δ^{2/3}q^{1/3}a^{1/3}, minimised at a = 1. In Prop O′, β = A ≥ 1 | refuted for this model |
| 2 | Composite / CRT locks ℓ₁ℓ₂ \| Bp − Aq | M cancels in g″ (Prop O′(a) is M-free). Each giant covers the same number of candidates of a smaller class | refuted for this model |
| 3 | q-side extraction (exponent 1 + RN ≡ q + Rp mod q−1) | Same target value Rp + q, same window. In q-units the per-giant coverage is ×R and so is the candidate count | refuted for this model |
| 4 | Finer ratio sub-bands, one list each | Per sub-band Δ^{4/3}(p/R)^{1/3}, × R/Δ sub-bands ⇒ Δ^{1/3}(Rq)^{1/3}, worse | refuted (computed) |
| 5 | Coarse small weights (1, B ≤ H²), full band window (GPT-6.1-sol) | Same exponent H^{4/3}q^{1/3}, simpler. Now the main construction | proved (§3.1), no gain |
| 6 | Composite moduli for multiplicity (Σ_M 1/M ≍ 1 per dyadic range vs 1/n for primes) | Coverage per target is still ≈ 1/M. Re-optimising gives λ^{2/3}q^{1/3} independent of the modulus density | refuted (computed; consistent with Prop O′) |
| 7 | Lucas (p+1), cyclotomic Φ_k(p) handles | p+1: N ≡ −q, target −(Bp + Aq), same family. Cyclotomic: cycloharvey D1 (no new targets) | refuted for this model |
| 8 | Exponents of degree ≥ 2 in N (targets with q², …) | g″ ≍ q²/p² ≫ R/p, larger curvature. No construction beats linear | not proved in general (outside 𝔏) |
| 9 | Coppersmith / lattice windows | Need p mod M with M ≳ p^{1/2}. Locks give p mod ℓ only, and for R ≥ ℓ every class is hit (no information). Guessing costs ≥ √p | no gain found (heuristic accounting) |
| 10 | Pollard–Strassen on lock classes | Must cover a union of classes of density ≍ 1 ⇒ √p | no gain (computed) |
| 11 | McKee–Pinch pigeonhole (λp ≡ ±μq mod k, λμ ≤ k) | Guaranteed lock at k ≈ N^{1/3} ⇒ deterministic N^{1/3}. Worse | no gain |
| 12 | HH 2022 local lattices | Accounting S + √N/S^{3/2} with A = O(√(S/R)) recovers N^{1/5} at best (GPT-6.1-sol) | no gain |
| 13 | Structured (non-interval) baby sets exploiting the common curvature g″ = 2N/p³ of all A = 1 targets | Prop O′ allows arbitrary sets and gives the same bound | refuted for this model (new in O′) |

**Further qualifications** (GPT-6.1-sol `unbalanced_2`).
- With finite H the DCL/PS cost crossover is R_* = N^{1/5}H^{−16/5}.
- Factor-dependent structure beats the profile on **special** inputs. For example, if p − 1
  is z-smooth, Pollard p − 1 costs Õ(z), which is below the peak when z < N^{1/5}. But there
  is no density theorem making that almost-all.
- A fixed power saving N^{−ε} requires a fixed ε. Letting ε → 0 to get exponent 1/6 + o(1)
  does not keep a fixed power saving.

**Conclusion.** Every attempt either lies in 𝔏, where Prop O′ makes q^{1/3} optimal, or
costs more. Beating the penalty needs a mechanism outside 𝔏: non-linear targets with
cancelling curvature, multi-dimensional collisions, lattices, or an average-case selection
that favours the true factor. We found none.

## 7. Experiments (pre-registered; notes/unbalanced_prereg.md, commit b3beb0d)

All four frozen experiments were run exactly as specified: designs, sizes, seeds and metric
E. The deviations are execution-only (§9). Raw: `experiments/raw/unbalanced_{lemma,exc,oper1_W1,oper1_C1,oper2}.jsonl`.
Summaries: `experiments/analysis/unbalanced_summary.json`, `unbalanced_decisions.json`.
The algorithms saw only N; classes and certificates are computed afterwards by the
harness. These are toy-size correctness and cost-law checks, not evidence of practical speed.

### 7.1 lemma (W1, n = 32, 40, j = 0..6, 10 pairs each; all 4 blocks recorded)

- 140/140 correct and factored. 140/140 certified, and certified ⇒ success.
- **0 violations** of L1 (integrality), L2 (window), L3 (mod p), L4 (mod q closed form) and
  L57 (clean block: gcd = p ⇔ not FULL) on **1905 covering targets**, 5 of them FULL.

### 7.2 oper1 (W1 vs C1 = D′ on the same N; 24 pairs per cell; all 4 blocks run)

| n | j | W1 success (block 0) | W1 med E_full | C1 success (block 0) | C1 med E_full | W1/C1 | C1 block-0 struct. coverage (p+q side) | W1 block-0 p-coverage |
|---|---|---|---|---|---|---|---|---|
| 48 | 0 | 24 (24) | 2.27e7 | 24 (24) | 2.91e6 | 7.79 | 5.79 | 9.92 |
| 48 | 2 | 24 (24) | 2.28e7 | 24 (20) | 2.94e6 | 7.76 | 1.96 | 10.33 |
| 48 | 4 | 24 (24) | 2.23e7 | 24 (15) | 2.88e6 | 7.76 | 0.46 | 9.42 |
| 48 | 6 | 24 (24) | 2.27e7 | 22 (17) | 2.92e6 | 7.79 | 0.25 | 4.83 |
| 64 | 0 | 24 (24) | 2.56e8 | 24 (24) | 3.73e7 | 6.86 | 6.17 | 8.79 |
| 64 | 2 | 24 (24) | 2.57e8 | 24 (20) | 3.74e7 | 6.86 | 1.67 | 9.25 |
| 64 | 4 | 24 (24) | 2.57e8 | 24 (13) | 3.74e7 | 6.88 | 1.33 | 10.58 |
| 64 | 6 | 24 (24) | 2.56e8 | 23 (19) | 3.72e7 | 6.88 | 0.33 | 6.17 |
| 80 | 0 | 24 (24) | 1.56e9 | 24 (24) | 2.36e8 | 6.61 | 4.75 | 7.79 |
| 80 | 2 | 24 (24) | 1.56e9 | 24 (18) | 2.35e8 | 6.62 | 1.25 | 9.08 |
| 80 | 4 | 24 (24) | 1.55e9 | 24 (10) | 2.34e8 | 6.63 | 0.54 | 7.58 |
| 80 | 6 | 24 (24) | 1.56e9 | 21 (7) | 2.35e8 | 6.63 | 0.25 | 4.50 |

- W1: 288/288 correct, 288/288 certified, 0 lemma violations on 3314 covering targets.
- First-hit classes (W1, all cells): STRUCTURAL_P 198, STRUCTURAL_Q 36, ACCIDENTAL 54.
  ACCIDENTAL hits are concentrated at 48 bits, where the orders mod p are small. No null
  rate was computed for this run (rule 3 of AGENT_RULES asks for one; **not done**). As in
  D′, any gcd in (1, N) is a valid factor.
- C1 (unchanged D′) still factors 21–24/24 at j = 6. Its block-0 coverage collapses
  (0.25–0.33 vs 4.75–6.17 at j = 0), so it succeeds through later blocks and ACCIDENTAL
  hits: 7–17 of the first hits at j = 6. Those are not covered by any theorem.

### 7.3 exc (W1 certificate components only; 300 pairs per cell)

| n | j | certified | block-0 uncovered | e^{−λ̄₀} + 0.04 (bound) | λ̄₀ | block-0 clean | uncovered in all 4 blocks | FULL / covering |
|---|---|---|---|---|---|---|---|---|
| 64 | 0 | 1.000 | 0.017 | 0.063 | 3.75 | 0.967 | 0 | 0 / 3970 |
| 64 | 3 | 1.000 | 0.063 | 0.087 | 3.07 | 0.953 | 0 | 0 / 3238 |
| 64 | 6 | 1.000 | 0.010 | 0.059 | 3.98 | 0.957 | 0 | 0 / 4114 |
| 80 | 0 | 0.997 | 0.050 | 0.087 | 3.07 | 0.947 | 1 | 0 / 3312 |
| 80 | 3 | 1.000 | 0.093 | 0.122 | 2.51 | 0.953 | 0 | 0 / 2657 |
| 80 | 6 | 1.000 | 0.020 | 0.079 | 3.25 | 0.960 | 0 | 0 / 3530 |

Uncovered fractions sit near the Poisson value e^{−λ̄₀}, uniformly in the ratio band, as
Theorem U's uniformity predicts. The proved rate needs H = n^{1+a}; these schedules have
H ≍ √n.

### 7.4 oper2 (W2, fine-band regime II; 16 pairs per cell; all 4 blocks run)

| n | j | Y | Δ | success | certified | med list block 0 | n-only expected list | ratio | med E_block0 | first hits P/ACC |
|---|---|---|---|---|---|---|---|---|---|---|
| 64 | 7 | 353 | 26 | 16 | 16 | 6 138 | 6 178 | 0.99 | 1.65e6 | 15/1 |
| 64 | 9 | 93 | 23 | 16 | 16 | 7 504 | 7 484 | 1.00 | 2.05e6 | 13/3 |
| 80 | 7 | 2666 | 35 | 16 | 16 | 40 747 | 39 722 | 1.03 | 1.43e7 | 14/2 |
| 80 | 9 | 593 | 29 | 16 | 16 | 47 565 | 47 318 | 1.01 | 1.62e7 | 16/0 |
| 80 | 11 | 156 | 22 | 16 | 16 | 58 934 | 58 370 | 1.01 | 2.22e7 | 15/1 |
| 96 | 7 | 23951 | 45 | 16 | 16 | 258 355 | 259 855 | 0.99 | 1.27e8 | 15/1 |
| 96 | 9 | 5334 | 38 | 16 | 16 | 298 094 | 300 052 | 0.99 | 1.41e8 | 15/1 |
| 96 | 11 | 1412 | 31 | 16 | 16 | 373 631 | 373 730 | 1.00 | 1.80e8 | 13/3 |
| 96 | 13 | 373 | 26 | 16 | 16 | 453 458 | 466 253 | 0.97 | 2.35e8 | 13/3 |

- 144/144 correct and factored, 144/144 certified, 0 lemma violations on 1406 covering
  targets.
- **Slopes of log₂(median) against j:**

  | n | E_block0 | list (block 0) |
  |---|---|---|
  | 64 | 0.157 | 0.145 |
  | 80 | 0.158 | 0.133 |
  | 96 | 0.152 | 0.138 |

  The n-only prediction is ≈ 0.14 and the asymptotic value is 1/6 ≈ 0.167.

### 7.5 Frozen decision rules

| rule | criterion | observed | verdict |
|---|---|---|---|
| R-a | 0 incorrect, 0 lemma violations, certified ⇒ success | 0 incorrect in 572 W runs + 288 C1 runs; 0 violations on 6625 covering targets; 0 certified failures | **holds** |
| R-b | W1 success ≥ 23/24 per oper1 cell; exc certified ≥ 0.99; uncovered ≤ e^{−λ̄₀} + 0.04 | 24/24 in all 12 cells; certified ≥ 0.997; every uncovered fraction below its bound | **holds** (regime I is uniform in the ratio) |
| R-c | W1/C1 E_full at j = 0 in [2, J1+2] = [2, 8]; drift 80/48 in [0.5, 2] | 7.79 / 6.86 / 6.61; drift 0.85 | **holds** (overhead ≈ J1+1 = 7, constant). The 48-bit value is near the upper end |
| R-d | C1 block-0 structural coverage at J1 ≤ 0.3 × at j = 0; secondary C1 success < W1 at J1 | ratios 0.043 / 0.054 / 0.053; C1 22/23/21 vs W1 24/24/24 | **holds** (both) |
| R-e | slope of log₂ med E_block0 vs j in [0.08, 0.24] at n = 80, 96; success ≥ 14/16; refutation if ≤ 0.04 | 0.158, 0.152; 16/16 everywhere | **holds**; the R^{1/6} law is **not refuted** |
| R-f | median block-0 list within ×1.5 of the n-only expectation | 0.97–1.03 | **holds** |

**Frozen numeric predictions.**
- **W1 λ_j at n = 64.** The frozen values 3.75 / 3.20 / 3.48 / 3.07 / 2.97 / 2.79 / 3.98
  are design quantities. The exc means at j = 0, 3, 6 (3.75, 3.07, 3.98) reproduce them
  exactly. That checks the bookkeeping only.
- **C1 multiplicity at band centre.** The predicted collapse 3.75 → 0.09 is mirrored by the
  measured structural coverage, which falls to ≈ 5% (R-d).
- **W2 lists.** All nine within 3% of the frozen n-only values.

**Reading.**
- Regime I is uniform in the ratio at a constant overhead ≈ J + 1, as Theorem U says.
- Regime II's cost grows with the ratio at slope 0.15–0.16 in log₂R. That is the R^{1/6}
  law at these sizes, the penalty Prop O′ says is unavoidable in 𝔏.
- Everything is correct end to end.
- There is no practical speed claim. At these sizes W1 costs ≈ 7× D′ on the same N, and D′
  itself is 2–5× slower than HH (dcltheorem §8).
- The coarse U_θ construction was **not** among the frozen runs. Only the critic's
  exploratory n = 48 checks exist for it (0 violations).

## 8. Reviews

- **GPT-6.1-sol blind, `codex/unbalanced_1.out`** (relaunched by the orchestrator).
  - It confirms the rectangle regime (with the sharper span (√2−1)H√N) and the fine-band
    cost H^{4/3}N^{1/6}R^{1/6}.
  - It supplied the coarse (1, B ≤ H²) construction, now Theorem U_θ.
  - It corrects the finite-H crossover to R_* = N^{1/5}H^{−16/5}.
  - It confirms the curvature bound for full coverage by consecutive windows, with the
    qualification that this is not an average-case bound.
  - It proposed the Harvey–Hittmeir 2026 extraction (§3.4).
  - It checked the literature. Lehman, McKee, McKee–Pinch, Hittmeir 2020, Harvey 2021, HH
    2022, Gao et al. 2512.19076 and HH 2026 contain no sub-1/5 intermediate-ratio theorem.
    It repeats the Carella 1308.2891 content defect.
- **GPT-6.1-sol r16 strategy** (§2 "unbalanced"). Same envelope, same full-cover caveat. It
  warned about gcd(B−1, d) and about the N ≤ X population trap; both are handled (§3.2, §4).
- **Independent Opus critic** (read-only; its own numerical checks).
  - Coarse design at n = 48 (h = 24 and 64, several bands): all pairs factored, 0 lemma
    violations.
  - Prop O′(a): 400 random targets, 0 violations, worst ratio 0.28.
  - §4 numbers reproduced.
  - **Verdicts:** Theorem U ACCEPT; U_θ ACCEPT WITH FIXES; U_θ′ ACCEPT; Prop O′ ACCEPT WITH
    FIXES. All fixes were applied:
    - span constant 18Q_j for h/16 < R ≤ h;
    - the dropped h^{2/3} restored;
    - the "what fails at 1/5" wording;
    - T_j stated in U_θ′;
    - (α, β) ≠ (0, 0);
    - rational normalisation D | βM²;
    - the **mislabelled** tightness witness (U_θ does not provably cover a constant
      fraction of candidates for fixed N), replaced by explicit 𝔏 witnesses;
    - the N ≤ X constant ln((1−δ)/δ);
    - §0/§4 overclaims reworded;
    - §3.4 precisions.
- **GPT-6.1-sol adversarial referee, `codex/unbalanced_2.out`** (xhigh).
  - **Verdicts:** Theorem U ACCEPT; U_θ ACCEPT WITH FIXES (the same span repair, found
    independently); U_θ′ ACCEPT (it confirms the swapped FULL count and X > 0); Prop O′
    **REJECT as written**.
  - The rejection is for the same constant-target counterexample (α = β = 0: cost 2, full
    coverage) and the same optimality-quantifier gap. Its repair, (α, β) ≠ (0, 0), is the
    one applied, and with it the proof is CORRECT per that referee.
  - Further fixes applied:
    - p ≥ √(N/(2R)) on line "Lemma 1";
    - public m̄;
    - §3.4 operational completion;
    - R_*;
    - fixed-ε vs ε → 0;
    - special smooth inputs;
    - four blocks suffice for W2 only when θ < (1+16ε)/23.
- **Post-fix status.**
  - Theorem U, U_θ, U_θ′: proved, both reviewers concurring after fixes.
  - Prop O′ in its repaired form (non-constant targets, floor relaxation, per-N
    constant-fraction coverage): proved; both reviewers' repairs applied. The "optimal in 𝔏"
    claim is now limited to relaxed coverage, with explicit witnesses.

## 9. Reproduction and disclosed deviations

```
.venv/bin/python -m pytest -q tests/test_unbalanced.py tests/test_unbalanced_obstruction.py
PYTHONPATH=. UNB_PROCS=8 UNB_MAXTASKS=0 UNB_RESUME=1 .venv/bin/python -m src.unbalanced.experiments lemma
PYTHONPATH=. UNB_PROCS=8 UNB_MAXTASKS=0 UNB_RESUME=1 .venv/bin/python -m src.unbalanced.experiments exc
PYTHONPATH=. UNB_PROCS=8 UNB_MAXTASKS=0 UNB_RESUME=1 .venv/bin/python -m src.unbalanced.experiments oper1
PYTHONPATH=. UNB_PROCS=8 UNB_MAXTASKS=0 UNB_RESUME=1 .venv/bin/python -m src.unbalanced.experiments oper2
PYTHONPATH=. .venv/bin/python -m src.unbalanced.experiments summarize   # also writes decisions
```

The frozen `config` step was **not** re-run: it would overwrite the frozen predictions file.

**Disclosed post-freeze changes.** None changes any design, instance, algorithm or metric.
1. `experiments.run_wide` crashed on `(q/p).bit_length()` (a float). It is now
   `(q//p).bit_length() − 1`, a display field only. This was found before any oper run;
   the lemma runs were started after the fix.
2. Execution knobs:
   - `UNB_PROCS`: 8 instead of 30 on this 10-core shared machine.
   - `UNB_MAXTASKS`: worker recycling disabled. With the frozen value 4, every recycled
     worker recomputed D′'s n-only schedule (~200 s; profile: `params._expected_cost`).
   - `UNB_RESUME`: on restart, run only the missing (design, n, j, seed) jobs.
   - Runs were interrupted and resumed twice: exc at 1151/1800 records (18:53), and oper1
     at 80/288 W1 records (20:41). The second time, the killed pool's orphaned workers kept
     computing and were killed by hand; their results were never written.
   - Every record is a deterministic function of its job key, and the final files have
     1800 / 288 / 288 / 144 distinct keys.
3. `analyze.decisions()` evaluates the frozen rules R-a…R-f. Its only interpretation:
   R-d's "p-side + q-side structural coverage" is `cover_p + cover_q` of block 0, as
   recorded by `_struct_cover`.
4. `obstruction.py` and `tests/test_unbalanced_obstruction.py` are new diagnostics.

## 10. Prior art (curated; claims limited to "no prior formulation found by us")

- Lehman 1974: rational ratio approximations, deterministic N^{1/3}.
- McKee–Pinch 1996 §5: weighted locks λp ≡ ±μq (mod k), O(max(λp, μq)/k²) per lock; a
  pigeonhole N^{1/3}; a heuristic balanced N^{1/4}. Their per-lock cost already shows the
  ratio dependence max(λp, μq).
- Harvey 2021 and Harvey–Hittmeir 2022: deterministic N^{1/5} for all N (all ratios at
  once).
- Pollard 1974 / Strassen 1976 / Bostan–Gaudry–Schost 2007: deterministic z^{1/2} for
  factors ≤ z.
- Harvey–Hittmeir 2026 (arXiv 2601.11131): deterministic large-order elements, used only in
  the §3.4 remark.
- No intermediate-ratio result below N^{1/5} was found by GPT-6.1-sol (two sweeps) or by us.
  Theorem U_θ's (1+θ)/6 for θ < 1/5 is, as far as we found, new as an almost-all statement.
  We claim only "no prior formulation found by us".
- The curvature obstruction is the ratio-dependent form of multimod Thm A (1)/Lemma C, with
  arbitrary baby sets added.

## 11. Candidate cards (summary)

- **DCL-wide-I (Theorem U).**
  - Representation: (Z/N)^*, base 2.
  - Factor feature: g = (Bp + Aq − c)/ℓ² in the ratio-band window.
  - Proved: the exponent of D′ for ratios ≤ N^{ε}.
  - Experimentally supported: §7.
- **DCL-wide-II (Theorem U_θ / U_θ′).**
  - Proved: almost-all exponent (1+θ)/6 + O(ε) for θ ≤ 1/5 − 3ε.
  - Scoped obstruction: Prop O′. In model 𝔏, q^{1/3} is necessary for full coverage.
- **Refuted for this model:** composite locks, q-side, anisotropic or rational weights,
  multipliers, sub-band lists, Lucas/cyclotomic handles, structured baby sets (§6).
- **Most promising open question.** An average-case (non-density) mechanism, or a
  non-linear / two-dimensional target family, that beats the curvature barrier near
  p ≈ N^{2/5}. Equivalently: is there any deterministic almost-all exponent below 1/5 for
  pairs with p ≈ N^{2/5}?
