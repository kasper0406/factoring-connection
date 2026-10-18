# dchworst: the worst case of the D′ scheme

Agent `dchworst`, 2026-10-18 (restarted on the new machine).
- Pre-registration: `notes/dchworst_prereg.md` (commit bde6c64; addendum A in commit 08176bb,
  before the e6 runs).
- Code: `src/dchworst/{coverage,theory,experiments,wcalg,experiments_wc}.py`.
- Tests: `tests/test_dchworst.py`, `tests/test_dchworst_wc.py`.
- Raw data: `experiments/raw/dchworst_e{1..6}*.jsonl`.
- Blind GPT-6.1-sol: `codex/dchworst_1.out`. Review round: `codex/dchworst_2.out`. Opus critic: §8.

## 0. Verdict (plain language)

**The worst case of D′ now rests on coverage alone. Coverage itself is still out of reach.**

- **Obstacles (b) and (c) of D′ §6 are removed (Theorem R, proved).**
  - (b) was small ord_q(2); (c) was all-FULL pairs.
  - DCL-wc is a modified D′. It is deterministic, costs N^{1/6+o(1)} plus the design's list
    cost, and factors **every** balanced N = pq whose pair is covered by at least one block.
  - It needs no clean band, no order hypothesis modulo the hidden q, and no FULL exclusion.
  - Three standard ingredients do the work:
    - a base-selection dichotomy, an elementary one-sided adaptation of Harvey–Hittmeir 2026
      (Lemmas 2.1, 2.3, 2.4, §3.1);
    - a γ-order pre-step on the baby list;
    - the exploration of every two-sided collision down the baby product tree. An exact
      collision there reveals Bp + Aq, hence p. Exact-match handling is itself prior art:
      Harvey 2021, Alg. 4.2, steps 1 and 3–4.
  - Both reviews accept Theorem R with fixes (§8); the fixes are applied. The theorem is
    about the fast variant. `src/dchworst/wcalg.py` is a correctness prototype whose cost
    does not meet (iii); see §3.
  - The orchestrator asked whether small ord_q(2) ≤ T is exploitable. It is: BSGS costs
    Õ(T^{1/2}) per base and factors N unless ord_p(b) = ord_q(b). That equal-order case is
    handled by other bases and the smooth-subgroup argument (§3, Lemma O).
- **Consequences.**
  - **Theorem W, unconditional cost criterion.** If coverage of every pair can be proved
    with H ≤ √Y·J^{−c}, the deterministic worst-case exponent is
    min(1/5, (1−c)/(7/2−c)). This is below 1/5 iff c > 3/8.
  - Minkowski gives c = 0, so this framework reproduces Harvey's 1/5 worst case (not new).
  - The pair form of DCH (a polylogarithmic H suffices for every pair) now implies a
    deterministic worst-case N^{1/6+o(1)} for balanced semiprimes.
  - Theorem D′'s exceptional set shrinks to the pairs uncovered by every block. That set has
    fraction O(n²/H²); E_a and E_c disappear.
- **The worst-case coverage gain is not proved (direction B: no).** I found no design that
  provably covers every pair with c > 0, nor even a logarithmic gain. §5 gives precise
  barriers:
  - positivity arguments whose error budget is ≳ √J reach only c = 1/4, and c > 3/8 is needed
    (a heuristic barrier on methods);
  - single-modulus pigeonhole gives only c = 0;
  - for any proof that also works for coprime non-prime pairs, coverage implies a covering
    statement U(x, H) for a structured set of ≍ H² = Y^{1/4} integers. We know no method
    for sets this sparse.
- **Negative direction (A): proved only in the range H ≤ c√log N.**
  - Brun–Titchmarsh plus a union bound shows that a positive proportion of pairs is
    uncovered (Theorem N).
  - Forcing one fixed ratio at many band primes is impossible by size (proved). Whether CRT
    forcing can construct uncovered pairs beyond the union-bound regime is UNVERIFIED. An
    earlier "at most 6 primes" density argument was refuted by GPT-6.1-sol (§4).
  - For √log N ≪ H ≪ log N, uncovered pairs are only predicted by the Poisson model, which
    the data support. A proof needs a lower-bound sieve at density λ → ∞ over one dyadic
    band, which bounded-level Bonferroni cannot give.
- **Experiments (pre-registered e1–e5, plus e6 for Theorem R).**
  - The exhaustive worst case over all of S_32 and S_36 (6.2M and 78M pairs) follows the
    random (Poisson) extreme: the fitted ĉ matches in all 8 cases (P4).
    - Frozen P1 failed through a mis-stated prediction (the excluded (1,1) weight).
    - Frozen P2 failed in 4 of 125 cells. That is the number the model itself predicts
      (4.47); one cell is a lone outlier with no structure found.
    - Frozen P3 failed almost only on the side of *better-than-random* coverage (120 of 121
      out-of-range cells).
  - The sampled tails at n = 60–84 (e2), the structured families (e3, ≤ 6% deviations), the
    exhaustive D′-design coverage (e4) and the all-blocks search at n = 64 (e5) all pass.
  - e6, Theorem R: **0** failures on covered pairs, exhaustively over S_26 and S_28 and on 36
    adversarial cells. D′ failed on 42 + 105 covered pairs of S_26, S_28 and on 6 adversarial
    covered pairs; DCL-wc factored every one. Details are in §7.
- **DCH repair (orchestrator item 1, verified).**
  - The paper's algorithmic DCH (all N coprime to 6, every integer candidate, A ≠ B) is false
    on squares: for N = t², p = t, Bp² − AN = (B − A)t² has no band divisor.
  - The repaired form is in §6.

## 1. Setting and notation

- N = pq with p < q ≤ 2p primes, n bits. S_n is the population of such pairs, as in D′, and
  P = 2^{n/2}.
- A **weight** is a primitive (A, B) with 1 ≤ |A|, B ≤ H.
  - *Positive* weights form D′'s V (A ≥ 1, A ≠ B).
  - *Signed* weights also allow A < 0, and we allow the diagonal (1,1).
  - The height is max(|A|, B).
- A **block** is a set 𝒫 of primes ℓ ≡ 3 (mod 4) with H < ℓ < √(N/2). Theorem R needs
  only H < ℓ. The large-sieve statements need ℓ > H² for positive weights and ℓ > 2H² for
  signed weights, so that the classes are distinct.
- The pair (p,q) is **covered** by (𝒫, V) iff ℓ | Bp − Aq for some (A,B) ∈ V and ℓ ∈ 𝒫.
- For ℓ ∤ pq put t_ℓ = p q^{−1} mod ℓ. Then ℓ | Bp − Aq ⇔ A B^{−1} ≡ t_ℓ. So coverage
  depends only on (t_ℓ)_ℓ.
- μ_ℓ(t) is the least height of a weight with A B^{−1} ≡ t (mod ℓ), the first minimum of
  the lattice Λ_ℓ = Z(p,q) + ℓZ² over the weight set.
- h_L(p,q) = min_{ℓ∈L} μ_ℓ(t_ℓ).

## 2. Theorem W (cost criterion) — proved

Let J ≤ (1/2 + o(1))Y/ln Y prime moduli ℓ ≡ 3 (mod 4) in (Y, 2Y], so that
j := log_Y J ≤ 1 − o(1) and 0 ≤ c ≤ 1/2, cover every pair with weights of height H = √Y·J^{−c}. The DCL
list cost is X = HP/Y² + JH², up to N^{o(1)}. Theorem R below supplies the extraction
unconditionally. With J = Y^j, j ∈ [0, 1], balancing the two terms gives the exponent of N

  E(j, c) = (1 + j(1−2c)) / (5 + 2j(1−c)).

∂E/∂j has the sign of (1−2c)·5 − 2(1−c) = 3 − 8c. So the optimum is

  min(1/5, (1−c)/(7/2−c)),

which is below 1/5 iff c > 3/8. Values: c = 1/4 and c = 3/8 give 1/5; c = 2/5 gives 6/31;
c = 1/2 gives 1/6. (`theory.worst_case_exponent` and the brute-force check in
`test_cost_exponent_formula`; GPT-6.1-sol dchworst_1 recomputed this independently.)

*Remarks.*
- X is the list cost at m = W. With m optimised the cost is ≈ √(T′W) + T′. Its minimum over
  y = log_P Y matches E(j, c) whenever E < 1/5 (critic), so the formula is the achievable
  exponent of this design family.
- The floor 1/5 is the j = 0 case, i.e. Minkowski (Corollary R1).

For J ≍ Y/log Y, the target c > 3/8 means: for every pair, some primitive (A, B) with
|A|, B ≤ Y^{1/8−δ} has Bp − Aq divisible by a prime in (Y, 2Y], where P ≍ Y^{25/8}.

## 3. Theorem R (worst-case extraction from coverage) — proved

**Algorithm DCL-wc** (`src/dchworst/wcalg.py`; operational: it sees only N and the public
parameters H, V, blocks 𝒫_1..𝒫_K, T and B0).

0. If N < 4, output FAIL. If N is even, output 2. If N = s² is a perfect square, output s.
1. **Base selection (Lemma O).** For b = 2, 3, …, B0:
   - check gcd(b, N);
   - compute r_b, the least r ∈ [1, T_scan] with gcd(b^r − 1, N) > 1, where
     s = ⌊√T⌋ + 1 and T_scan = s⌈T/s⌉ = T + O(√T). Use baby-step giant-step by polynomial
     evaluation: f(X) = ∏_{i=1}^{s}(X − b^{−i}), evaluated at X_j = b^{js}. Then scan the first block j whose value has
     gcd > 1.
   - If that gcd is proper, output it. If no r ≤ T_scan exists, **select b** and go to step 2.
     Otherwise (gcd = N) record r_b.

   If every b ≤ B0 was recorded, put Λ = lcm_b r_b and test gcd(1 + jΛ, N) for
   1 ≤ j ≤ (⌊√N⌋ − 1)/Λ.
2. For each block k, with L = ∏_{ℓ∈𝒫_k} ℓ², u_ℓ = b^{L/ℓ²} and γ = b^L:
   - **Targets.** Build the D′ targets τ = (A, B, ℓ, r, c, g_lo, W), with roots of
     Br² ≡ AN (mod ℓ) and the class lock mod ℓ². This is D′ §2.2 step 2 verbatim. For A < 0
     the window bounds are Lo = B a₀ + ⌊AN/a₀⌋ and Hi = B b₀ + ⌈AN/b₀⌉, since
     x ↦ Bx + AN/x is increasing.
   - Choose the baby length m as in D′.
   - **γ pre-step (Lemma G).** Over the babies γ^i (0 ≤ i < m), find the least r ∈ [1, m)
     with gcd(γ^r − 1, N) > 1, using a product tree of the γ^r − 1 and a binary descent.
     - If that gcd is proper, output it.
     - If it equals N, put o = r and run the **equal-order split**. E = L·o is a multiple of
       ord_N(b) whose factorisation is known: the band primes, and o by trial division. For
       each prime power r^e ∥ E, compute y_r = b^{E/r^e} (all at once by the all-but-one
       exponent tree), then gcd(y_r^{r^j} − 1, N) for j = 0..e, and output any proper gcd.
       Otherwise D = ∏ r^{v_r}, with v_r the least j where y_r^{r^j} ≡ 1 (mod N); test
       gcd(1 + jD, N) for 1 ≤ j ≤ (⌊√N⌋ − 1)/D.
   - **Giants.** G_{τ,j} = u_ℓ^{B+AN−c}·γ^{−(g_lo + jm)} for 0 ≤ j < ⌈W_τ/m⌉. Compute
     f(G) = ∏_i(G − γ^i) at all giants by a remainder tree.
     - Output any proper gcd(f(G), N).
     - Collect the **N-giants**, those with gcd = N.
   - **Exploration (Lemma X).** For each N-giant, descend the subproduct tree of f. At a
     node whose value has gcd = N, evaluate both children: output any proper gcd, and
     descend into the children whose gcd is N.
   - **Exact leaf.** When leaf i (the factor X − γ^i) is reached with G ≡ γ^i (mod N), put
     S = c + ℓ²(g_lo + jm + i) and test whether S² − 4ABN = D₀² is a perfect square. If it
     is, output a proper gcd((S ± D₀)/2, N).
3. Output FAIL.

All outputs are validated (1 < d < N, d | N).

**Theorem R.** Let N = pq ∈ S_n, and let the parameters satisfy H < ℓ < √(N/2) for every
band prime. Then:
- (i) DCL-wc outputs a proper divisor of N or FAIL;
- (ii) if (p,q) is covered by (𝒫_k, V) for at least one k, DCL-wc outputs a proper
  divisor;
- (iii) its cost is

  Õ( Σ_k (m_k + T′_k + |V||𝒫_k| + |𝒫_k| log Y) + B0·T^{1/2} + √N/T + √N/Ψ(√(N/2) − 1, B0) )

  ring operations mod N plus Õ(|V||𝒫_k|) word operations for target construction. Here
  T′_k = Σ_τ⌈W_τ/m_k⌉ is the giant count.
  - The |V||𝒫_k| term charges the rejected (weight, prime) combinations. It is needed for an
    arbitrary V (GPT-6.1-sol: an adversarial V with 9 592 tested weights and 2 targets).
  - For D′'s V and for sign-symmetric signed families it is ≤ H²|𝒫_k|, already part of the
    list cost. With ℓ ≡ 3 (mod 4), exactly one of ±AN B^{−1} is a square, so accepted
    targets pay for the tests.

  With T = ⌈N^{1/3}⌉ and B0 = ⌈n⁴⌉ the last three terms are N^{1/6}·n^{O(1)}, so the cost
  is N^{1/6+o(1)} + Õ(Σ_k(m_k + T′_k + Y)).

*Proof.*

**Lemma O (one-sided order dichotomy).**
- Statement. For b with gcd(b, N) = 1, write o_p = ord_p(b) and o_q = ord_q(b). If
  min(o_p, o_q) ≤ T_scan, the scan returns r_b = min(o_p, o_q), possibly r_b ∈ (T, T_scan].
  The gcd at r_b is proper iff o_p ≠ o_q. If no r ≤ T_scan is found, then
  o_p, o_q > T_scan ≥ T.
- Proof of the dichotomy.
  - b^r ≡ 1 (mod p) iff o_p | r, so the least r with gcd(b^r − 1, N) > 1 is
    min(o_p, o_q).
  - Say r = o_q ≤ o_p. If o_p ≠ o_q then o_p > r, so p ∤ b^r − 1 and the gcd is q.
  - If o_p = o_q, the gcd is N.
  - The BSGS scan visits r = js + i in increasing order inside the first block j with a
    nontrivial value, because b^{js+i} − 1 = b^{i}(b^{js} − b^{−i}).
- Proof of the all-recorded case. Suppose every b ≤ B0 is recorded, so
  ord_p(b) = ord_q(b) = r_b for all b ≤ B0.
  - G_p = ⟨2, …, B0⟩ ≤ F_p^× is cyclic of order lcm_b ord_p(b) = Λ. Likewise mod q.
  - G_p contains every B0-smooth integer in [1, p−1], and these are distinct mod p. So
    Λ ≥ Ψ(p−1, B0), and Λ | p − 1.
  - Hence p = 1 + jΛ with j ≤ (⌊√N⌋ − 1)/Λ (p < √N), and the enumeration finds p.
  - The elementary bound Ψ(x, (log x)^A) ≥ x^{1−1/A}/(log x)^{A−1} follows from counting
    products of k = ⌊log x/(A log log x)⌋ primes ≤ y = (log x)^A: there are
    ≥ (π(y)/k)^k ≥ (log x)^{(A−1)k} of them. With A = 4 and B0 = ⌈n⁴⌉ ≥ (ln p)⁴, this gives
    √N/Λ ≤ P^{1/4}n³.
- The method is the one Harvey–Hittmeir use for large order mod N (arXiv 2601.11131 v2,
  Lemmas 2.1, 2.3 and §3.1, smooth-number subgroup). The one-sided reading "first r
  ⇒ factor unless the two orders are equal" is the only change. ∎

**Lemma G (γ pre-step).** By the same argument applied to γ over r < m:
- either a proper divisor is output;
- or ord_p(γ) = ord_q(γ) = o < m;
- or **both ord_p(γ) ≥ m and ord_q(γ) ≥ m**.

In the equal case, γ^o ≡ 1 (mod N) gives o_p(b), o_q(b) | E = L·o.
- For each prime r | E, y_r = b^{E/r^e} has order r^{v_r(o_p)} mod p and r^{v_r(o_q)}
  mod q.
- If v_r(o_p) ≠ v_r(o_q), then at j = min of the two, gcd(y_r^{r^j} − 1, N) is p or q.
- Otherwise o_p = o_q = D, and D is computed. Since b was selected in step 1,
  D = o_p(b) > T. Also D | p − 1, so the enumeration finds p, at cost ≤ √N/T.
- (In the test mode `force_base`, where step 1 is skipped, the enumeration is still complete
  but its cost is not bounded.) ∎

**The true target.** Let (A, B) ∈ V, ℓ ∈ 𝒫_k with ℓ | Bp − Aq.
- Then a = A N B^{−1} ≡ p² (mod ℓ) is a nonzero square, and p mod ℓ is one of the two
  computed roots r*. This is D′'s coverage criterion; the sign of A plays no role.
- Lemma 1 (class lock): g = (Bp + Aq − c)/ℓ² ∈ Z. Lemma 2 (window): g_lo ≤ g ≤ g_hi, also
  for A < 0 by monotonicity.
- Let j* = ⌊(g − g_lo)/m⌋ and i* = g − g_lo − j*m.
- Mod p: B + AN − (Bp + Aq) = (p−1)(Aq − B), so u_ℓ^{B+AN−c} ≡ u_ℓ^{ℓ²g} = γ^g. Hence
  G* := G_{τ*,j*} ≡ γ^{i*} (mod p) (D′ Lemma 3; no order hypothesis needed).

**Lemma X (exploration).** Assume both ord_p(γ) ≥ m and ord_q(γ) ≥ m (Lemma G).
- (a) Every giant G has at most one *p-root* i ∈ [0, m) with G ≡ γ^i (mod p), because the
  γ^i are distinct mod p.
  - A node of the subproduct tree with gcd(f_c(G), N) = N contains a p-root.
  - Hence at most one such node per level.
  - The exploration therefore evaluates at most 2⌈log₂ m⌉ node polynomials per N-giant.
- (b) If G* is an N-giant, its exploration outputs a proper divisor.
  - Every node c on the root-to-leaf-i* path satisfies p | f_c(G*), so gcd(f_c(G*), N) ∈
    {p, N}.
  - By induction along the path, either some node yields p (output), or leaf i* is reached
    with G* ≡ γ^{i*} (mod N).
  - In the latter case S = c + ℓ²(g_lo + j*m + i*) = c + ℓ²g = Bp + Aq.
  - Then S² − 4ABN = (Bp − Aq)² is a perfect square, and {(S ± D₀)/2} = {Bp, Aq}.
  - gcd(Bp, N) = p because 1 ≤ B < q, and gcd(|Aq|, N) = q because 1 ≤ |A| < p. ∎

**Proof of (ii).** Suppose block k covers (p,q).
- Steps 0–1 either output a divisor or select b.
- In block k, Lemma G either outputs a divisor or gives orders ≥ m mod p and q.
- Then G* is computed.
  - If gcd(f(G*), N) = p, it is output.
  - Otherwise gcd = N and Lemma X(b) applies.
- Earlier blocks can only output divisors or continue. Their FAIL branches, the
  enumerations, are complete by Lemmas O and G, so they cannot abort. ∎

**Proof of (iii).**
- Step 1 costs Õ(B0 T^{1/2}) plus the Λ-enumeration.
- Per block:
  - the exponent tree costs O(|𝒫_k| log Y) mults;
  - targets, giants and the remainder tree cost Õ(m + T′) (D′ §4);
  - Lemma G costs Õ(m): a product tree plus a binary descent;
  - the equal-order split costs O(|𝒫_k| log Y) plus √N/T.
- Exploration, by Lemma X(a), batched level by level, costs Σ_levels Õ(m + 2#N-giants) =
  Õ(m + T′). Multipoint evaluation of each node polynomial is done at its ≤ 2#N-giants
  assigned points.
- ∎ for the fast variant.

**The implementation is a correctness prototype, not a cost witness** (both reviews). Three
places fall short of the theorem's bound:
- `run_block` computes each u_ℓ = b^{L/ℓ²} separately, Θ(J² log Y) mults instead of the
  all-but-one tree, which D′'s `all_but_one` implements;
- `equal_order_split` computes each b^{E/r^e} separately;
- Lemma G scans r linearly, and the exploration evaluates each node at each giant
  separately, O(m) per N-giant.

None of these changes an output. The two "FAIL" aborts in `select_base` and
`equal_order_split` are unreachable for N = pq by Lemmas O and G, as noted in the code.

**Corollaries (proved).**

- **R1 (Minkowski, unconditional worst case N^{1/5+o(1)}).**
  - Design: signed weights with (1,1), H = ⌊√(2Y)⌋, and a block of 1 prime ℓ ∈ (Y, 2Y].
  - If t_ℓ = 1, the weight (1,1) covers. Otherwise pigeonhole on (a, b) ∈ [0, ⌊√ℓ⌋]² gives a
    primitive (A, B) with |A|, B ≤ √ℓ, B ≥ 1 and A ≡ t_ℓB. So every pair is covered.
  - The cost is HP/Y² + O(Y) + N^{1/6+o(1)}, minimised at Y = P^{2/5}: N^{1/5+o(1)}.
  - This recovers **Harvey's proved exponent 1/5** (Harvey 2021, Math. Comp.) inside the
    DCL framework at c = 0; it is not new.
  - McKee–Pinch 1996 §5 is the weighted-congruence precursor, but its balanced method is
    heuristic N^{1/4+ε}, and its deterministic construction is N^{1/3+ε} (GPT-6.1-sol, from
    the scans).
- **R2 (DCH_pair ⇒ worst case 1/6).** Suppose the following, call it DCH_pair(κ). For all
  large n and every (p,q) ∈ S_n, some block of D′'s design with H = ⌈n^κ⌉ covers (p,q),
  with weights from V or with signed weights. Then DCL-wc factors every N ∈ S_n
  deterministically in N^{1/6}·n^{O(1)}·H^{5/3} bit operations, i.e. N^{1/6+o(1)}.
  - DCH_pair is a hypothesis. Theorem N shows it needs κ ≥ 1/2.
- **R3 (D′ sharpened).** With D′'s own design, DCL-wc fails only on pairs uncovered by every
  block. That set is contained in {(p,q): p ∈ U_0(q)}, uncovered in block 0, without any
  condition on ord_q(2), so it is not D′'s E_b, which carries d > T and the clean block.
  - D′ Lemma 9 applied to the fixed block 0 gives an exceptional fraction O(n²/H²).
  - The r16/sieve2d prime-presieve denominator (sieve2d Prop. 4; checked by GPT-6.1-sol)
    gives O(n/H²). That refinement is not mine.
  - E_a (small ord_q(2)) and E_c (FULL) are gone, and so is D′'s restriction ε ≤ 1/60. Only
    Y > H² (the sieve's distinct classes, ε < 1/7) and H < ℓ remain.
- **R4 = Theorem W made unconditional.** §2 applies with extraction supplied by Theorem R.

Answer to the orchestrator's cost question:
- BSGS for ord_q(2) ≤ T costs Õ(T^{1/2}) ring ops.
- D′'s E_a threshold is T = N^{1/4−ε}, which costs N^{1/8}.
- It finds a factor unless ord_p(2) = ord_q(2). That case needs Lemma O's other bases and
  the smooth-subgroup argument (≤ n⁴ bases, N^{1/6}n^{O(1)} with T = N^{1/3}). It is not
  automatic.

## 4. Negative direction (A)

**Theorem N (proved; union bound and Brun–Titchmarsh).** Let L be any set of primes in
(Y, 16Y] with 16Y ≤ P^{θ}, θ < 1/2, and let V be any weight set. Put
λ = |V|·Σ_{ℓ∈L} 1/(ℓ−1). Then

  #{covered pairs of S_n} ≤ (2/(1−θ) + o(1))·λ·|S_n| + o(|S_n|)    (when |V|·|L| = o(P)).

*Proof.*
- For fixed (A, B, ℓ) and a prime p, the partners q with ℓ | Bp − Aq form one class mod ℓ
  in the interval I_p of q-values.
- Montgomery–Vaughan's Brun–Titchmarsh gives ≤ 2|I_p|/(φ(ℓ) log(|I_p|/ℓ)).
- The p with |I_p| < P/n² lie within O(P/n²) of the endpoints of the p-range. For each
  condition (A, B, ℓ) they contribute O(P/n² · (P/(n²ℓ) + 1)) by the trivial count. Summed
  over all conditions, this is O(λ·|S_n|/n² + |V||L|P/n²) = o(λ|S_n|) + o(|S_n|).
  (Fix by the critic and GPT-6.1-sol: the earlier P^{0.99} cutoff does not give the
  constant.)
- For the others, log(|I_p|/ℓ) ≥ (1 − θ − o(1)) log P.
- Summing over p and comparing with |S_n| = (1 + o(1))Σ_p |I_p|/log P (PNT on each I_p)
  gives the bound. ∎

**Consequence.** With θ ≈ 1/3 the constant is 3. For D′'s four blocks and positive weights,
λ ≈ (6/π²)H²·(ln 16/2)/ln Y ≈ 7.3H²/n. Signed weights give 14.6H²/n. So
**a positive proportion of S_n is uncovered by every block for H ≤ 0.2√n (positive
weights) or H ≤ 0.14√n (signed weights)**, and DCH_pair(κ) fails for κ < 1/2.

Note that "uncovered" means Theorem R gives no guarantee. DCL-wc can still succeed through
accidental collisions; this is frequent at toy sizes (e6: uncovered_ok) and heuristically
negligible at real sizes.

**T2 (CRT forcing).**
- *Refuted as first drafted.* The draft said "prescribing t_ℓ at F needs ∏(ℓ−1) ≤ |S_n|,
  hence |F| ≤ 6". Density below 1/|S_n| does not exclude a pair (GPT-6.1-sol dchworst_2).
  - Counterexample: (2617, 5233) satisfies its own ratios at F = {19, 23, 31, 43, 47, 59, 67}.
  - Yet ∏(ℓ−1) ≈ 8.8·10¹⁰ > P².
- *Heuristic remark.* Each forced prime changes λ by only |V|/(ℓ−1) = O(H²/Y). Turning
  forcing into a construction needs equidistribution at level ∏ℓ, so this route is
  UNVERIFIED.
- *Proved obstruction (GPT-6.1-sol).* Forcing one fixed ratio r = A/B at every prime of a set
  F with ∏_F ℓ > 2HP is impossible: ∏_F ℓ | Bp − Aq ≠ 0 and |Bp − Aq| ≤ 2HP.

**The range √log N ≪ H ≪ log N (conjectured: uncovered pairs exist).**
- The exact local density of uncovered classes is Δ = ∏_ℓ(1 − |V|/(ℓ−1)). With signed
  weights and D′'s blocks, log Δ = −(36 ln 16/π²)H²/ln N + o(1) (GPT-6.1-sol, exact).
- The predicted uncovered count is therefore |S_n|e^{−λ} with λ = (144/π²)H²/n. It tends to
  ∞ iff H < (π√ln 2/12 − o(1))·n ≈ 0.218n (critic's correction of "H = o(n)"). The
  experiments e1/e2/e5 agree with this model (§7).
- A proof needs a lower bound of size e^{−λ} with λ → ∞. Truncated inclusion–exclusion to
  depth r needs equidistribution modulo products of r band primes, of size Y^r.
- Even granting level P² = Y^{6+o(1)} for pairs, only r ≤ 6 is available, and odd
  Bonferroni truncations of depth ≤ 5 stay positive only for λ < 2.18.
- Fundamental-lemma and β-sieve bounds use sieve dimension ≍ |V| ≍ H² on the band, not the
  intensity λ ≍ H²/n. Their sifting limit exceeds the available s = log D/log Y ≤ 6 once H
  is a large constant.
- So standard sieve technology proves nothing in this range. This is a method barrier, not
  an impossibility theorem.

## 5. Positive direction (B): no provable gain; precise barriers

**What would suffice.**
- c > 3/8 (Theorem W). At J = Y: for every pair, a primitive weight of height ≤ Y^{1/8−δ}
  with a band divisor of Bp − Aq.
- Even c ∈ (0, 3/8] gives no exponent gain. A log-factor gain in H gives only a log-factor
  speedup, of the same kind as Harvey–Hittmeir 2022 or Gao et al. 2025.
- None of these is proved.

**B1. Positivity barrier (proved conditional calculation; heuristic as a barrier).**
- The covering count M = Σ_{v,ℓ} 1[ℓ | f(v)] has mean λ ≍ JH²/Y.
- *Proved algebra.* If an argument "M = λ + O(E)" has error budget E ≳ J^{1/2}, then
  positivity forces H ≫ Y^{1/2}J^{−1/4}, i.e. c ≤ 1/4. That is exponent 1/5 (Theorem W).
- At c > 3/8 the main term is λ ≤ J^{1/4}. Such an argument needs E = o(J^{1/4}).
- *Model, not theorem.* J^{1/2} is what square-root cancellation over the J moduli gives in
  each Fourier mode e(h(Bp − Aq)/ℓ). A per-modulus bound gives more.
- We do not claim every counting method must incur this error. The true fluctuation of M is
  only ≍ √λ. (GPT-6.1-sol r16 states the same diagnostic: α(1/4) = 3/13 without the 1/5
  floor.)

**B2. Single-modulus pigeonhole (proved, for that strategy).**
- Minkowski/pigeonhole certifies a short vector in one sublattice
  Λ_M = Z(p,q) + MZ², i.e. divisibility by one fixed M. Coverage is a union over ℓ, which is
  not a lattice statement.
- The single-modulus guarantee is optimal over classes. Weights of height ≤ h hit at most
  W(h) = (12/π² + o(1))h² classes, so for h < 0.9√M some class has μ_M(t) > h
  (0.9 = √(π²/12)).
- For composite M even √M fails for primitive vectors (GPT-6.1-sol). M = 4k with
  t = 2k − 1 forces height ≥ M/4. The pair (p, q) = (7, 11) realises this at M = 12.
- Hence this strategy gives c = 0. We do not claim that every geometry-of-numbers argument
  is so limited.

**B3. Interval content (proved conditional implication).**
- *Premise.* The proof is assumed to work for all coprime pairs (p, q), not only prime pairs;
  natural sieve and lattice methods use only ℓ ∤ pq. Prime-pair coverage alone cannot be
  specialised to (x, x+1).
- Under that premise, take (p, q) = (x, x+1). The values become
  Bx − A(x+1) = bx − a with (A, B) = (a, a+b).
- So coverage at height H implies **U(x, H)**: for every x ≍ P, some integer of
  {bx − a : 1 ≤ |b| ≤ 2H, |a| ≤ H, (a, a+b) a weight} has a prime factor in the band.
- U(x, H) is a covering statement for a structured set of ≍ H² integers, a union of ≍ H
  short intervals of length ≤ 2H around multiples of x.
- At c > 3/8 the set has H² = Y^{1/4} integers near x ≍ Y^{25/8}, against band primes ≈ Y.
- *Context, not a lower bound.* The single-interval statement "[x − H, x] contains a
  multiple of a band prime" implies U, and the converse is not shown.
  - The best unconditional single-interval results we found need length x^{1/3−η} at
    x ≍ Y³, i.e. Y^{1−3η}, for prime factors (Jutila 1974 via Friedlander–Lagarias 1987,
    Lemma 2.2).
  - For all divisors in (Y, 2Y] they need x^{27/97+ε} ≈ Y^{0.84} (Graham–Kolesnik
    exponent pairs; Weingartner arXiv 2105.13568).
  - Melfi's MathOverflow question 437849 (2023) asks the divisor version at length log²x,
    Y = x^{1/3}. We found it unanswered.
- We know no method that proves covering statements for sets as sparse as U(x, Y^{1/8}).
  We did not find that U implies a *named* conjecture, and neither did GPT-6.1-sol.

**B4. Other routes considered (no gain found; assessments, not theorems).**
- Composite or smooth moduli: B2 (primitivity) and the root/class cost cancel the gain
  (r16, Thm C neutrality).
- Divisors in intervals (Ford 2008) are distributional statements.
- Continued fractions of q/p turn the values into a rank-2 GAP {aD_i + bD_{i−1}} over a
  parallelogram. The problem is SL₂(Z)-invariant, so this changes coordinates, not
  difficulty.
- 3D Minkowski: the equation Bp − Aq = Mt is bilinear in (M, t), not a lattice condition.
- N-dependent moduli: t_ℓ is unknown, and choosing ℓ from N alone cannot see the lattice
  Λ_ℓ.
- Umans–Wang's prime-cover theorem (arXiv 2511.10851, Thm 5.5) is conditional on SPDC.
  He–Sahai (arXiv 2608.06681) refute the AP form of SPDC but not DCH.

**Sharpest barrier statement.**
- *Sufficient (proved; Theorems R and W).* A proof that every balanced pair has a band-prime
  divisor of some Bp − Aq with height Y^{1/2}J^{−c}, c > 3/8, gives a deterministic
  worst-case exponent below 1/5.
- *Necessary* only for proofs that go through Theorem R's coverage guarantee with this
  design family. Accidental collisions and other designs are not excluded.
- For such a proof:
  - (i) a counting argument needs an error budget o(J^{1/4}), beyond square-root
    cancellation over the moduli (B1, heuristic);
  - (ii) single-modulus pigeonhole is capped at c = 0 (B2), so a union-type existence
    argument of a new kind is needed;
  - (iii) if the proof also covers coprime non-prime pairs, it proves U(x, Y^{1/8}) (B3).

No general factoring lower bound prevents < 1/5. The barrier is specific to this method.

## 6. DCH: which version, and the repair (for paperwriter; orchestrator item 1)

- **Verified counterexample.** paper/main.tex Conjecture DCH (algorithmic form, ~line 1119)
  quantifies over every integer candidate p ∈ [√(N/2), √N], all N coprime to 6, and
  requires A ≠ B.
  - Take N = t² with t prime and p = t. Then Bp² − AN = (B − A)t² with 0 < |B − A| < H < ℓ
    < t, so no band prime divides it.
  - The conjecture is false as stated, for every κ. Signed weights do not help when
    ℓ > 2H, which D′'s asymptotic design satisfies: (B + |A|)t² then has no band divisor.
    With only H < ℓ they can: for H = 4, ℓ = 7, the weight (−3, 4) gives 7t².
- **multimod §3.6 DCH(κ)** allows A = B: coprime A, B ≤ (log N)^κ, so (1,1) is included.
  For N = t², (1,1) gives 0, which every ℓ divides. The notes' version is not refuted by
  squares.
- **Recommended statement.** Theorem R needs only the true factor. Use the *pair* form:
  - **DCH_pair(κ).** For all large n and every (p, q) ∈ S_n, some block 𝒫_k
    (k = 0..3) of D′'s design with H = ⌈n^κ⌉ contains ℓ with ℓ | Bp − Aq for some
    primitive (A, B), 1 ≤ |A|, B ≤ H. The diagonal (1,1) is allowed.
  - With Theorem R, DCH_pair(κ) ⇒ deterministic worst-case N^{1/6+o(1)} on S_n. A ≠ B is no
    longer needed: exact identification handles the two-sided weight (1,1).
- The candidate-wise form (every integer p) is the 𝔐⁺ full-covering statement of
  multimod. If kept, it must either exclude N = square or allow A = B (two alternative repairs).


## 7. Experiments (pre-registered; notes/dchworst_prereg.md)

All runs used the code frozen at commit bde6c64, for e1–e5, or at 08176bb with addendum A,
for e6. The analysis is in `src/dchworst/analyze.py` and
`experiments/analysis/dchworst_summary.json`.

### 7.1 e1: exhaustive worst case of h_J over all of S_32 (6 185 272 pairs) and S_36 (78 295 606 pairs)

Band ℓ ≡ 3 (mod 4) in (Y, 2Y], with nested prefixes. The table shows h*_J, the worst case
over all pairs, against h_pred(J), the least H with |S_n|·Pr_Poisson[h > H] < 1.

| n | Y | weights | J=16 | J=64 | J=128 | J=\|L\| (\|L\|) |
|---|---|---|---|---|---|---|
| 32 | 2^13 | signed | 67 / 65 | 38 / 40 | 29 / 29 | 19 / 19 (440) |
| 32 | 2^13 | positive | 92 / 93 | 56 / 56 | 40 / 42 | 28 / 26 (440) |
| 36 | 2^13 | signed | 71 / 68 | 41 / 43 | 33 / 32 | 19 / 19 (440) |
| 36 | 2^13 | positive | 99 / 98 | 61 / 60 | 46 / 45 | 28 / 28 (440) |
| 32 | 2^11 | signed | 34 / 33 | 23 / 21 | 17 / 17 (128) | |
| 36 | 2^11 | signed | 35 / 35 | 23 / 23 | 17 / 18 (128) | |

The frozen design-only h_pred values in the pre-registration are reproduced exactly: for
n = 32 signed, J = 16 → 65, …, 440 → 19.

**Verdicts on the frozen predictions.**

- **P1: FAILS as stated.**
  - In every (n, Y), h*_1 exceeds ⌊√ℓ₁⌋.
  - The pairs above the bound are **exactly** the pairs with p ≡ q (mod ℓ₁): 11 314,
    2 560 and 428 at n = 32; 147 929, 36 271 and 8 153 at n = 36. The weight (1,1) was
    excluded from the signed set, and for t = 1 every other lattice vector has height
    ≥ ℓ/2.
  - The correct statement, with (1,1) allowed as in Theorem R, holds for all pairs. That is
    Corollary R1.
  - This is an error in the prediction, which labelled itself "proved; sanity" without the
    exclusion. It is not a structured obstruction.
- **P2: FAILS in 4 of 125 cells**, all with positive weights and all with h* > h_pred + 3:
  - (32, 2^13, J = 12): 111 vs 105, pair (45007, 72977);
  - (32, 2^13, J = 96): 55 vs 47, pair (37957, 70181);
  - (36, 2^11, J = 12): 65 vs 58, pair (166613, 222149), with 4p − 3q = 5;
  - (36, 2^13, J = 12): 121 vs 115, pair (154669, 245849).

  The decision rule requires analysing these pairs.
  - Under the Poisson model, the expected number of P2 cells with at least one pair above
    h_pred + 3 is Σ_cells min(1, E[#{h > h_pred+3}]) = **4.47**. The ±3 tolerance was too
    tight for the model's own extreme-value fluctuation, and 4 exceedances is what the
    model predicts.
  - Three of the four worst pairs have E[#{h > h* − 1}] = 0.06–0.27, ordinary extremes.
  - The J = 96 pair is a lone outlier, with E[#{h > 54}] = 0.004. We found no structure in
    it:
    - no small Bp ± Aq (A, B ≤ 12) has a band divisor;
    - its signed/positive split (79 of 128 primes where the positive minimum exceeds the
      signed one) and its mean μ are typical (random pairs: 58–75);
    - q/p has no unusually good approximation (min |Bp − Aq| = 3559 for A, B ≤ 12).

    By J = 128 it is covered at h ≤ 39.
  - **Verdict: no structured obstruction identified.** The prediction failed through a
    too-narrow tolerance (3 cells) plus one unexplained single-pair outlier.
- **P3: FAILS in 4 of 12 (n, Y, sign) groups.**
  - Over all groups, 121 of 8 241 cells fall outside [0.7, 1.4]. All but one have **fewer**
    uncovered pairs than Poisson predicts; the exception is one signed cell with 39 observed
    against 26.5 expected.
  - In the failing groups every outside cell is on the low side. The structure therefore
    helps coverage:
  - With signed weights, Minkowski's cutoff makes μ_ℓ ≤ √ℓ deterministic except at t = 1.
  - Two simultaneous t = 1 events would need ℓ₁ℓ₂ | q − p with ℓ₁ℓ₂ > q − p, which is
    impossible.
  - The Poisson model, which treats the primes as independent, puts mass there.
  - The remaining groups pass (≥ 98.3% of cells within [0.7, 1.4]).
- **P4: PASSES in all 8 (n, Y ∈ {2^11, 2^13}, sign) cases.** ĉ_obs − ĉ_pred lies in
  [−0.014, +0.025]; ĉ ≈ 0.32–0.40.

**Decision.** e1's frozen rule asked: P2–P4 all hold ⇒ "the random-model extreme; no
structured obstruction". As frozen, this did not hold. Our reading after the analysis above:
- the exhaustive worst case is the random-model extreme up to one single-pair outlier;
- the only systematic deviations make coverage *better* than random.

Read the label as *experimentally supported, with the disclosed P1/P2/P3 failures*.

**Local exponent.** h*_J ∝ J^{−ĉ} with ĉ ≈ 0.35–0.40 over J ∈ [16, 440]. The model predicts
the same. It is a finite-range slope, not the asymptotic c: for fixed Y and J → Y/log Y the
model gives h ≍ √(log|S|·Y/J). It is consistent with c → 1/2 (DCH) and says nothing about
the provable c.

### 7.2 e2: sampled tail at the D′ scaling (n = 60, 72, 84; 200k / 100k / 50k pairs; Hcap = 64)

**P2 passes in all 30 cells**: |h* − h_pred| ≤ 2 everywhere. **P3 passes in all 30 cells**:
100% of cells within [0.5, 2], ≥ 95.2% within [0.7, 1.4]. Worst case over all four blocks:

| n | signed h* / pred | positive h* / pred |
|---|---|---|
| 60 | 8 / 8 | 11 / 11 |
| 72 | 8 / 8 | 12 / 11 |
| 84 | 9 / 9 | 12 / 12 |

### 7.3 e3: structured families (n = 72, block 0, 231 primes; 20 000 pairs per family; baseline 100 000)

Uncovered fraction at H_ref (positive H = 8, signed H = 6), as a ratio to the random
baseline:

| family | positive | signed |
|---|---|---|
| twin_gap | 0.997 | 0.955 |
| near_2p | 1.005 | 1.031 |
| near_3_2 | 1.021 | 1.054 (z = 3.4) |
| near_5_4 | 1.048 (z = 3.3) | 1.046 |
| cong_p_M30 | 0.987 | 0.980 |
| cong_mp_M30 | 1.032 | 1.061 (z = 3.9) |
| crt_worst2 | 1.030 | 1.006 |

**The prediction passes**: every ratio lies in [0.5, 2], and crt_worst2 lies in [0.8, 1.25].
No family is a structured adversary (ratio > 2 with a 3σ excess).
- The small but significant +5–6% excesses come from one dead weight each.
- For near-rational q/p, the residual of the approximating weight is tiny and has no band
  divisor.
- For q ≡ −p (mod ∏_{r≤29} r), the value p + q = M30·(small) has no band divisor, so
  (−1,1) is dead. Expected effect ≈ λ/|V| ≈ 4%.

### 7.4 e4: exhaustive D′-design coverage (positive weights, 4 blocks, Y = 2^{n/6}/⌊H^{1/3}⌋)

| n | H | uncovered (any block) | Poisson | uncovered in every *clean* block |
|---|---|---|---|---|
| 32 | 4 | 297 659 | 290 102 | 416 331 |
| 32 | 5 | 19 954 | 18 646 | 50 612 |
| 32 | 6 | 4 455 | 4 121 | 18 994 |
| 32 | 7 | 33 | 29.4 | 1 551 |
| 32 | 8 | **0** | 0.13 | 120 |
| 28 | 7 | 1 | 1.9 | 50 |
| 24 | 6 | 1 | 2.0 | 32 |

**Passes.** All 14 cells with expected ≥ 20 have ratio in [0.886, 1.133]. The minimal
full-coverage H* is 7, 8, 8 for n = 24, 28, 32, equal to the least H with expected < 0.5.

The last column shows what Theorem R buys at these sizes. D′'s proof needs coverage in a
block clean for q; Theorem R needs any block. At n = 32, H = 8 that is 120 pairs against 0,
and it stays 4 at H = 10.

### 7.5 e5: random search for pairs uncovered by all four D′ blocks (n = 64, 10⁶ pairs, positive weights)

| H | Y | uncovered in all 4 blocks: observed / Poisson | per-block ratios obs/exp |
|---|---|---|---|
| 7 | 1625 | 4 721 / 4 701.6 (1.004) | 1.003, 1.003, 1.001, 1.003 |
| 9 | 813 | 78 / 75.0 (1.04) | 1.000, 1.003, 1.001, 1.000 |

**Passes**: the all-4 ratio is within [0.5, 2], and every block is within [0.8, 1.25]. The
blocks behave as independent: the all-4 count is the product of the block rates.

Explicit uncovered pairs at n = 64, the "constructions" of direction A, found by search:
- H = 7: (3569033173, 4719025931);
- H = 9: (2329183463, 4009445411).

These are worst-case witnesses for D′'s design at these H. By Theorem R, DCL-wc has no
guarantee on them.

### 7.6 e6: Theorem R (addendum A)

**e6a, exhaustive over all of S_26 (146 536 pairs) and S_28 (504 015 pairs)**, H = 3, 4 blocks:

| n | design | mode | covered | covered_fail | uncovered (succeeded anyway) | D′ fails on covered → DCL-wc ok |
|---|---|---|---|---|---|---|
| 26 | pos | full | 137 189 | **0** | 9 347 (4 966) | 42 → 42 |
| 26 | pos | force2 | 137 189 | **0** | 9 347 (4 821) | |
| 26 | sgn | full / force2 | 146 536 | **0 / 0** | 0 | |
| 28 | pos | full | 471 769 | **0** | 32 246 (16 624) | 105 → 105 |
| 28 | pos | force2 | 471 769 | **0** | 32 246 (16 395) | |
| 28 | sgn | full / force2 | 504 015 | **0 / 0** | 0 | |

Routes in (28, sgn, full):
- 412 119 one-sided collisions;
- 24 760 exact leaves, mostly the (1,1) weight;
- 63 889 by Lemma O's first_r, small orders;
- 2 809 by Lemma G;
- 438 by exploration.

All improper-output counts are 0. max(frontier/N-giants) = 1.0 in every cell.

**e6b, adversarial families, n = 40, 48, 56** (both designs, both modes; 36 cells): covered_fail = **0** in every cell and improper = 0.
- small2_q (q | 2^k − 1) and small2_both (p, q | 2^k − 1, 2^{k′} − 1): Lemma O (full) or
  Lemma G (force2) factors almost all of them before any DCL step.
- eqord, n = 48: D′ failed on **5 of 5** covered pairs. DCL-wc factored all 5, by exact leaf
  or collision, and under force2 two of them via the equal-order split and D-enumeration.
- small2_both, n = 56: D′ failed on 1 covered pair; DCL-wc factored it.

**Verdicts on the frozen predictions.**
- F1 (covered_fail = 0 in every cell): **passes** (8 + 36 cells).
- F2 (improper = 0): **passes**.
- F3 (frontier ≤ 1 per N-giant): **passes**.
- F4 (equal-order runs > 0 for eqord under force2 at every n): holds at n = 40 (27 runs,
  pos) and n = 48 (2 pos, 1 sgn). **Untested at n = 56**: the generator found no eqord pair
  in 3·10⁶ tries, so the family is empty and the cell absent.

Both reviewers ran further adversarial checks (§8). Their positives:
- the critic's q = 2p − 1 FULL family: 101 D′ failures, all factored by DCL-wc;
- GPT-6.1-sol's equal orders for bases 2 and 3: all factored.

**Post-run code change (disclosed, output-neutral).** After e6, `wcalg.run_block` and
`equal_order_split` were switched to D′'s all-but-one exponent tree, as the reviews asked.
An exhaustive comparison over S_24 (42 622 pairs × 2 designs) gives 0 differences in
(factor, status) against the e6 version (commit 08176bb).

## 8. Review

Theorem R and §§2, 4–6 went to two independent read-only reviewers. Both returned **ACCEPT
WITH FIXES** for Theorem R, and neither found a covered input that fails. All fixes are
applied above.

**Opus critic** (read-only subagent; scratch in the session scratchpad).
- *What it verified.* Every lemma, the smooth-number bound with natural logs, and the
  batched exploration bound.
- *Exhaustive runs with an instrumented copy.* n = 22 and 24 over several designs, both
  modes and forced base 3: 0 covered failures, 0 Lemma X(a) violations, 0 invalid outputs.
  All of the D_enum, order_split, lambda_enum, exact_leaf and explore paths were hit.
- *A FULL-forcing family, q = 2p − 1 at n = 44–52.* D′ failed on 101 covered pairs; DCL-wc
  factored every covered pair, 100 of the 101 D′ failures via exact_leaf.
- *Fixes applied:*
  - R3 scope: "uncovered by every block", which is not D′'s E_b;
  - R2 cost written as n^{O(1)}H^{5/3};
  - R1 relabelled "recovers", with Harvey 2021 attribution;
  - Õ(m) for Lemma G;
  - the code-vs-theorem disclosure;
  - Theorem W: prime moduli ≡ 3 (mod 4), j ≤ 1 − o(1), and the achievable-exponent remark;
  - Theorem N: cutoff P/n² and separate positive/signed constants;
  - T2 relabelled heuristic;
  - the H < 0.218n threshold;
  - "equivalent" replaced by sufficient / conditionally necessary;
  - B3's arithmetic and scale comparison (the "25×" claim was wrong; the Jutila comparison
    holds only at x ≍ Y³) and rephrasing;
  - B1 labelled a heuristic barrier;
  - §6 "equivalently" dropped.

**GPT-6.1-sol, codex/dchworst_2.out** (xhigh, read-only, own adversarial scripts in
`codex/scratch/dchworst_2/`).
- *Runs.* 31 592 single-weight covering runs, 1 409 full and 1 399 forced-base runs, and
  40 runs at default parameters (T = ⌈N^{1/3}⌉, B0 = n⁴) at 39–80 bits, all succeeded. They
  included:
  - equal orders of bases 2 and 3;
  - different orders of b collapsing to equal γ-orders;
  - band primes shared with p − 1 and q − 1;
  - negative-weight exact leaves.
- *Fixes applied:*
  - **(major)** the code does not realise the cost bound: u_ℓ and the b^{E/r^e} are not
    computed by exponent trees, and the exploration is unbatched. It is now labelled a
    correctness prototype;
  - **(major)** rejected (weight, prime) tests must be charged: the |V||𝒫_k| term;
  - T_scan defined (first_r may return r ∈ (T, T_scan]; correct but undocumented);
  - Theorem N's cutoff;
  - **(major)** T2's density inference refuted by the explicit pair (2617, 5233), and
    withdrawn;
  - B1–B3 relabelled, and "equivalent" dropped;
  - prior art: exact-match handling is Harvey 2021 Alg. 4.2; the class lock is HH 2022
    Prop. 3.3; Lemma O adapts HH 2026 Lemmas 2.1, 2.3, 2.4 and §3.1;
  - McKee–Pinch's balanced method is heuristic N^{1/4+ε}, not the source of 1/5;
  - the signed-square qualification ℓ > 2H, and ℓ > 2H² for signed large-sieve classes.
- *Novelty assessment (theirs, which we adopt).* Lemma O is "an elementary adaptation of
  existing machinery". Theorem R as a whole is a combination of known tools. We claim only
  "no prior formulation found by us" for the statement "coverage by any block ⇒ worst-case
  extraction".

**GPT-6.1-sol blind attempt, codex/dchworst_1.out.** It found no worst-case gain either
(UNVERIFIED). Inputs we used:
- the composite-modulus primitive counterexample (B2);
- the interval literature (B3);
- the exact signed density constant 36 ln 16/π² and the CRT fixed-ratio obstruction (§4);
- an exact check: all 53 917 prime partners of q = 1 500 007 against 118 primes in
  (100, 1600] are covered for H ≥ 6.

## 9. Reproduction

```
PYTHONPATH=. .venv/bin/python -m pytest -q tests/test_dchworst.py tests/test_dchworst_wc.py
PYTHONPATH=. .venv/bin/python -m src.dchworst.experiments e1   # ~25 min, ~4.5 GB RAM (n = 36)
PYTHONPATH=. .venv/bin/python -m src.dchworst.experiments e2
PYTHONPATH=. .venv/bin/python -m src.dchworst.experiments e3
PYTHONPATH=. .venv/bin/python -m src.dchworst.experiments e4
PYTHONPATH=. .venv/bin/python -m src.dchworst.experiments e5
PYTHONPATH=. .venv/bin/python -m src.dchworst.experiments_wc e6a   # Theorem R, exhaustive n = 26, 28
PYTHONPATH=. .venv/bin/python -m src.dchworst.experiments_wc e6b   # Theorem R, adversarial families
PYTHONPATH=. .venv/bin/python -m src.dchworst.analyze             # verdicts -> experiments/analysis/dchworst_summary.json
```

## 10. Candidate card

- **Name.** DCL-wc: worst-case extraction for divisor-cover Lehman, plus the worst-case
  coverage problem.
- **Mechanism.** D′'s targets and class lock. Base selection by a one-sided order dichotomy
  (Lemma O, after Harvey–Hittmeir 2026). A γ-order pre-step (Lemma G). Exploration of
  two-sided collisions down the baby subproduct tree, with exact identification of
  Bp + Aq (Lemma X).
- **Proved.**
  - Theorem R: coverage by any block ⇒ success, in N^{1/6+o(1)} plus the list cost, for
    every balanced N.
  - Theorem W: the cost criterion min(1/5, (1−c)/(7/2−c)).
  - Corollary R1: Minkowski ⇒ worst case N^{1/5+o(1)}, which recovers Harvey's exponent.
  - Corollary R2: DCH_pair ⇒ worst case N^{1/6+o(1)}.
  - Corollary R3: D′'s failures are only uncovered pairs.
  - Theorem N: a positive proportion of pairs is uncovered for H ≤ 0.2√n.
  - The fixed-ratio CRT obstruction. The general T2 density claim was refuted and withdrawn.
- **Experimentally supported.**
  - The worst case over all of S_32 and S_36 is the Poisson extreme, with one unexplained
    single-pair outlier.
  - Structured families deviate by ≤ 6%.
  - DCL-wc: 0 failures on covered pairs in e6.
- **Refuted for this model.** Single-modulus pigeonhole (c = 0) cannot beat 1/5, and
  neither can positivity with error budget ≳ √J (c ≤ 1/4; heuristic as a barrier).
- **Open.**
  - Worst-case coverage with c > 3/8 (would give < 1/5).
  - DCH_pair(κ) for any κ (would give 1/6).
  - Existence of uncovered pairs for √log N ≪ H ≪ log N.
- **Most promising open question.** A union-type existence theorem: for every balanced
  pair, a band-prime divisor of some Bp − Aq with height Y^{1/8−δ}. The first test case is
  the union-of-short-intervals statement U(x, H) for x = Y³ and H = Y^{1/8}.
