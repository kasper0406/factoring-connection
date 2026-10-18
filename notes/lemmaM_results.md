# lemmaM: fruitfulness, averaging over N, involutions, and a fruitfulness-free route (Dixon descent)

Agent `lemmaM`, phase 6, 2026-10-18.
- Pre-registration: `notes/lemmaM_prereg.md`, commit `e529d50`, made before any E1/E2 run.
- Code: `src/lemmaM/`. Tests: `tests/test_lemmaM.py`. Raw data: `experiments/raw/lemmaM_e{1,2}.jsonl`.
- GPT-6.1-sol: `codex/lemmaM_1.{md,out}` (blind round) and `codex/lemmaM_2.{md,out}` (adversarial round).
- Independent Opus critic: ACCEPT WITH MAJOR FIXES. All fixes are applied here; see §9.

## Verdict

**We obtained no rigorous L[1/3] factoring, no rigorous L[1/2, c < 1], and no proof of Lemma M or Conjecture F.** The
best outcome is a reduction of type 4/5.

- **Theorem DD [proved; a folklore-level argument].** Take a uniformly random square R = r² mod N and descend it to the
  factor base, as NFS discrete-log descent does: split R ≡ u/v by rational reconstruction, then find one special-q relation
  for each large prime of uv. Any dependency that contains such a relation is then fruitful with probability exactly
  1 − 2^{1−ω}.
  - This holds for every N that is not a prime power.
  - There is no Blum condition, no LV Conjecture 7.1, no Conjecture F, and no condition on the other relations.
  - The required hypothesis is that the exposed transcript is invariant under r ↦ εr.
  - The argument is the Rabin/Dixon fibre symmetry. nfsfruit's Theorem L1 is its single-shared-root special case.
  - What is new here is only the combination: hidden roots supplied by descent, so the relations have NFS size rather than
    Dixon size. In nfsfruit's terms, this is a *new relation set* that carries the β-flipping ε-action by construction
    (criterion (I) with δ′ = 1, restricted to descent relations). It says nothing about LV's own relation set.
- **What remains is relation supply: Hypothesis D (§3), rigorous one-level descent of uniform squares in an LV-randomized
  NFS.** It involves no Bockstein or 2-adic data. The only character condition left is that u/v be a square, which is the
  classical Dixon/Pomerance kind. Hypothesis D is nonetheless open and **not** close to known tools. Three ingredients are
  missing:
  1. A first moment. At fixed degree the progression modulus t = a − bm is |F|^{1−o(1)}, beyond every known level of
     distribution for smooth numbers (Harper, Pascadi 5/8). A design that respects level 1/2 changes the constant; see
     regime (a′).
  2. A *complete-target-list* average E_A q(A): the probability that **all** k special-q targets of one split succeed in
     one field. This is a k-point joint event, and its natural proof needs pair correlations with very small loss. LV avoid
     exactly these.
  3. A counting lemma for splits whose large primes fall in a field-dependent good set.
- **Conditional consequences [proved implications from Hypothesis D, with costs in the standard smoothness model].**
  - **(a)** Fixed odd degree d.
    - The formal constant is L[1/2, 2/√(d+1)] (d = 5: 0.816).
    - With LV-type coefficient inflation H ≥ M^{2+ε}, which keeps the moduli within level 1/2, it is
      **L[1/2, √(6/(d+1))]**. This is < 1 for odd d ≥ 7 (d = 7: 0.866; d = 9: 0.775).
  - **(b)** One-level descent with degree ≍ (log N)^{1/5} gives L[2/5]. The one-stage L[2/5] phenomenon is old: Gordon 1993
    for NFS discrete logs.
  - **(c)** This architecture **cannot reach rigorous L[1/3]**. Rigorous HECM splitting at B1 = L[2/3] alone costs
    L_{B1}[2/3] = L[4/9]. In general, 2β₁/3 ≤ α together with 1 − β₁ ≤ α forces α ≥ 2/5, so an ECM-speed rigorous
    smoothness test would be needed.
  - **(d)** The QS analogue, with special-q relations x² − N = ℓ·c and a QR filter on split primes, is a heuristic
    L[1/2, 1] randomized QS whose dependencies are provably fruitful. It removes hypothesis H2 of below16 O-4 at the price of
    randomness. It does not touch the deterministic question.
- **Averaging over N (angle a) [refuted for this model as a route to Lemma M].**
  - The decoupling parametrizations exist (§4.1).
  - The classical large sieve is loss-free at these lengths: Q ≥ H² (§4.2).
  - Codex's Theorem A gives balance of Jacobi signs on a fixed smooth set for almost all Blum semiprimes.
  - Codex's Theorem B shows that fixed-test averaging is compatible with δ = 0 for every N. So no fixed-test averaging
    argument yields Lemma M, and a structural control of the adaptive predictor is needed.
- **Involutions (angle b).**
  - The canonical x-preserving, β-flipping involution for QS/CFRAC is multiplication by the fundamental unit of Z[√N].
    - It is nontrivial for every Blum N [proved, classical Pell argument; E2: 80/80].
    - It displaces relations by a whole period, P ≈ N^{1/2+o(1)}/h typically (toy fit: slope 0.511).
  - An explicitly computable β-flipping pair factors N at once (codex square-twist identity). So a usable involution must
    be either a factoring algorithm in its own right or a non-constructive symmetry of the sampling law. Dixon and DD use
    the latter.
  - Odd-degree fields have no order-2 automorphisms; degree-1 forms have no GL₂ symmetry except ±I (codex).
  - Möbius symmetries permute coordinates and fix the root ratio (§5).
  - **(B1) [proved, scope-restricted].** In even-|S| two-form degree-1 sieves (nfsfruit Cor. J2, no known root of λ), every
    congruence is trivial whenever the least prime ℓ with (ℓ/N) = −1 exceeds B. So Conjecture F for that mechanism, for all
    Blum N, implies a least-nonresidue bound L[1/2], far beyond Burgess.
- **Class groups (angle c).** We agree with r18 and codex: under uniform class sampling, b + 1/(4b) ≥ 1, and norms are
  ≥ √|Δ|. r18's rejection of "hidden square multipliers appended to NFS" applies to one-step smoothness and does not apply
  to descent (§6).
- **Literature (angle d).** For 2018–2026 we found no rigorous fruitfulness result, no rigorous L[1/3], and no rigorous
  L[1/2, c < 1]. The closest prior art for DD is Stange (arXiv 2211.06821; random exponents g^x in index calculus mod n,
  heuristic NFS version), together with NFS-DL descent (Joux–Lercier, Commeine–Semaev, Gordon, Schirokauer, Al Aswad–Pierrot,
  all heuristic). We found no statement of DD itself, or of the conditional constants √(6/(d+1)) and 2/√(d+1) for factoring.
  This is a search result, not a novelty certification.
- **Experiments (pre-registered, exact):** E1, 960 toy QS-descent records at 48–72 bits; E2, 240 continued-fraction full
  periods at 28–40 bits. All pre-registered predictions held (§7). They are correctness checks only.

---

## 1. Theorem DD (Dixon descent)

**Setting.**
- N is odd with ω ≥ 2 distinct prime factors. Repeated factors are allowed: μ₂(Z/N) ≅ {±1}^ω since each (Z/p^e)^× is
  cyclic. Prime powers are removed first.
- A relation is a pair (X_i, c_i) with X_i ∈ (Z/N)^× and c_i in a certificate group 𝒢, with a homomorphism
  ev: 𝒢 → (Z/N)^× such that X_i² = ev(c_i).
  - Dixon/QS: 𝒢 = Q^× and ev is reduction mod N.
  - NFS: 𝒢 = Q^× × K^× and ev(q, γ) = q·φ(γ)^{−1} on elements coprime to N, with the LV/BLP f′ and f_d factors.
- If ∏_S c_i = g² and g has been computed, put X_S = ∏_S X_i and Y = ev(g). Then X_S² = Y², and Y is a unit.

**Descent relations.**
1. Draw r uniform in (Z/N)^×. A non-unit r is a factor of N: output it if 1 < gcd < N, and resample if gcd = N. The code
   draws r from [2, N−2]. This breaks exactness only on the fibre of R = 1, which has probability 2^ω/N.
2. Put R = r² mod N.
3. A procedure 𝒟(R; ξ) runs with fresh coins ξ and run-wide coins Ξ (for example the field (m, f)). Both are independent
   of r. It either fails or returns (W, c) with W²R = ev(c).
4. Set X = rW.

**Theorem DD [proved].** Let 𝒯 be the exposed transcript. It contains N; Ξ; all ξ_i, R_i, W_i, c_i; all other relations
(X_j, c_j); the choice of S; and g. Assume:
- (T) 𝒯 is invariant under every replacement r_i ↦ ε_i r_i (ε_i ∈ μ₂), where only the stored roots X_i = r_i W_i transform;
- S contains at least one descent relation.

Then Pr[1 < gcd(X_S − Y, N) < N | 𝒯] = 1 − 2^{1−ω} ≥ 1/2.

*Proof.*
1. Conditional on R_i, the root r_i is uniform on r_i·μ₂, and the draws are independent.
2. By (T), conditioning on 𝒯 leaves the fibre elements i.i.d. uniform on μ₂.
3. Fix any descent relation i₀ ∈ S. Then X_S/Y = ε_{i₀}·Z with Z 𝒯-measurable given the other ε's, so X_S/Y is uniform
   on μ₂.
4. Exactly two elements of μ₂ are ±1. For any other ζ, gcd(ζ − 1, N) is the proper divisor ∏_{ζ≡1 (p^e)} p^e. ∎

**Remarks.**
- (R1) **(T) is necessary.** For example, a relation with root r^{−1} and certificate R^{−1} makes the product
  identically 1. The selection of S, g or acceptance must not use X_i, Jacobi symbols of X_i, or earlier gcd outcomes.
  The statement is per dependency: after an observed trivial gcd, take fresh hidden roots or non-adaptive kernel samples.
- (R2) **No leak.** The split, the HECM or trial-division tests, the special-q searches, and LV's character and
  square-root steps act on R, u, v, ℓ, (a, b) and the field only. A factorisation of an integer computed from R is a
  function of R.
- (R3) **Consistency with the earlier barriers.** The rejected strategies needed the uniform square itself to be
  factor-base-smooth (nfsfruit_1 §3.1; lemmaM_1 (22): cost L[2/3]), or a uniform λ = r² with small ratio representatives
  (nfsfruit Prop B; lemmaM_1 (24)). Descent needs neither. "Re-randomising a known root is factoring" is respected: descent
  never sees a root.
- (R4) **Prior art.**
  - Rabin 1979 (square roots and factoring); Dixon 1981 and Pomerance's random-squares analyses (fibre symmetry);
    nfsfruit Theorem L1 (the same exact statement for one shared hidden root).
  - Stange 2022: random exponents in index calculus mod n, with order finding and rational linear algebra under her
    Hypothesis 3.1. Her NFS version is heuristic.
  - NFS-DL individual-log descent.
  - DD only packages these for an NFS factor base.

## 2. Cost of Dixon descent (heuristic smoothness model; implications from Hypothesis D)

**One-level descent.**
- *Split.* (u, v) is the Gauss-reduced shortest vector of {u ≡ Rv (mod N)}, so |u|, |v| ≤ √(2N). Accept iff u, v are
  B1-smooth.
- *Special-q.* For each odd-valuation prime ℓ ∈ (B, B1] of uv, search the special-q box of Λ_ℓ = {a ≡ bm (mod ℓ)}. This
  needs a regular reduced basis ≍ √ℓ; for random m the irregular case, m/ℓ close to a low-height rational, is rare. The
  search looks for primitive (a, b) with (a − bm)/ℓ and F(a, b) both B-smooth.
- *Substitution.* Then ℓ ≡ φ(a − bα)/c with c = (a − bm)/ℓ.
- *Rows and cost.* Use C·L^{o(1)} pure-descent rows; LV §6 needs extra character columns and combinations. Every nonzero
  dependency then contains a descent relation. Total cost ≈ C·L^{o(1)}·(k·T_SQ + T_split) + T_LA.

**(a) Fixed odd d, L[1/2] scale.**
- Use a non-monic f of degree d with m = N^{s+o(1)}, ‖f‖ = N^{s+o(1)}, s = 1/(d+1). Take B = L[1/2, a] and
  B1 = exp((log N)^{3/4}).
- Rigorous HECM (Lenstra–Pila–Pomerance) tests B1-smoothness in L_{B1}[2/3] = exp(O((log N)^{1/2}(log log N)^{1/3})) =
  L[1/2, o(1)].
- The split succeeds with probability ρ(u₁)² = L^{−o(1)}, u₁ = (log N)^{1/4}/2. Each split has k ≍ (log N)^{1/4} + O(log log N)
  large primes.
- The ℓ^{(d−1)/2} and box factors are N^{o(1)}. So both sides have size N^{θ+o(1)}, θ = 2s, and T_SQ = L[1/2, θ/(2a)].
- Optimum a = 1/√(d+1), giving **L[1/2, 2/√(d+1)]**: d = 3: 1; d = 5: 0.816; d = 7: 0.707. This is the formal constant of
  plain fixed-degree NFS.
- **(a′) Level-of-distribution-compatible design.**
  - With LV randomization, F(a, b) ≡ N b^d (mod t) with t = a − bm. At the sizes above, t = |F|^{1−o(1)} (codex lemmaM_2),
    which no BV-type theorem for smooth numbers reaches.
  - Inflate the coefficients of R to H ≥ M^{2+ε}. Then t ≤ |F|^{1/2−ε}, the range of LV's own §8 analysis.
  - This costs F ≍ N^{2s} and the rational cofactor ≍ N^{s}, so θ = 3/(d+1) and c = √(6/(d+1)): < 1 for odd d ≥ 7.
  - The moduli are still multiples of ℓ with a moving residue N b^d. Their first moment is part of Hypothesis D and is open.

**(b) One level, growing degree.**
- B = L[α], log B1 ≍ (log N)^{β₁}(log log N)^{1−β₁}, d ≍ (log N/log log N)^{(1−β₁)/2}.
- The split forces β₁ ≥ 1 − α. The special-q values have log-size ≍ (log N)^{(1+β₁)/2}(log log N)^{(1−β₁)/2}, which forces
  β₁ ≤ 4α − 1. Hence α ≥ 2/5, with equality at β₁ = 3/5, d ≍ (log N/log log N)^{1/5}, and log-sizes
  ≍ (log N)^{4/5}(log log N)^{1/5}.
- HECM on B1 = L[3/5] costs L_{B1}[2/3] = L[2/5].
- So L[2/5, c] for an unoptimised c. The one-stage 2/5 is in Gordon's NFS-DL analysis (1993). Only the factoring and
  fruitfulness packaging is ours.
- With k ≍ (log N/log log N)^{2/5}, a per-target success probability (log N)^{−O(1)} still gives L[2/5, c′].

**(c) No rigorous L[1/3] from this architecture.**
- Rigorous HECM at B1 = L[2/3] costs L[4/9]. More generally, HECM splitting needs 2β₁/3 ≤ α, and with 1 − β₁ ≤ α this
  again forces α ≥ 2/5 (codex lemmaM_2).
- One-level descent of L[2/3] primes costs L[1/2] at the re-optimised degree d ≍ (log N)^{1/6}, or L[2/3] at the L[1/3]-NFS
  degree.
- Multi-level descent would also make the special-q primes field-dependent. Calling this "strictly harder" is an
  assessment.

**(d) QS analogue (the toy of E1).**
- Special-q relation x² − N = ℓ·c with x ≡ ±√N (mod ℓ), scanned near √N. The cofactor ≍ 2√N·J is independent of ℓ.
- Splits are rejected if a large prime is inert in Q(√N). Example: N = 77, R = 16 ≡ 3/5, and both 3 and 5 are inert. The
  rejection costs 2^{−k} = L^{−o(1)}.
- Heuristic L[1/2, 1], with fruitfulness from DD.

## 3. Hypothesis D and what it needs

**Hypothesis D (open).** Use LV-randomized fields A = (m, f): m uniform; f = f̂_{N,m} + (x − m)Σ c_i x^{d−1−i}, c_i
i.i.d. uniform. Let q(A) be the probability that a capped one-level descent of a uniform square succeeds in field A. Then
Pr_A[q(A) ≥ L^{−o(1)}] ≥ L^{−o(1)}, with time cap T_D.

**Corollary DD [proved implication, modulo an adaptation of LV §6 that is not written out].** Hypothesis D implies that
every N is factored in expected time (C·T_D + T_LA)·L^{o(1)}.
- The algorithm: good-field search with restarts, capped descents, pure-descent rows, LV's character columns and square
  root, then Theorem DD.
- The LV §6 adaptation covers characters, square roots and the height bounds. LV's Thm 2.6 assumes pairs of height
  ≤ L[1/3]; ours have height S√ℓ.
- Regime (a′): L[1/2, √(6/(d+1)) + o(1)]. Regime (b): L[2/5, c].

**The three missing ingredients [assessment; codex lemmaM_2 and the critic concur].**
1. **First moment at the right moduli.** E_A[q(A)] needs smooth values of F(a, b) in progressions N b^d mod t, with t = ℓc
   a multiple of a fixed special-q prime, averaged over special-q boxes.
   - At fixed d with minimal height, t ≍ |F|^{1−o(1)}, so no known estimate applies.
   - Under (a′), the level 1/2 is respected, but the family of moduli (multiples of ℓ, moving residues) is not LV's
     averaging family.
   - Pascadi's 5/8 needs a fixed residue.
2. **The complete-target-list average.** The needed statement is Pr_{A,Z}[N_ℓ(A) > 0 for every target ℓ of the split Z]
   ≥ L^{−o(1)}, under the actual accepted-split law.
   - Reverse Markov then gives good fields, so stochastic deepening *does* apply once this average is proved. An earlier
     draft said otherwise; that was wrong (critic, codex).
   - The marginal "Lemma SQ" of an earlier draft does not imply it. Two field types could each descend complementary
     prime types (codex).
   - The natural proof is a second moment of pair-correlation type, with loss exp(o((log N)^{1/4})) in (a) and polylog in
     (b).
3. **Split counting.** The split law has saddle weights ℓ^{−σ}, not 1/ℓ. It carries the square-class condition
   u/v ∈ (Z/N)^{×2}, presumably handled by the Dixon/Pomerance class-pairing trick, but coprimality complicates it.
   Large primes must land in a field-dependent good set.

**Research programme, not a reduction.** Two randomized algebraic sides f_j = f̂_j + (x − m)R_j remove the deterministic
rational side. Descend only primes that split completely in K₁.
- The joint density of (⟨c, v⟩, ⟨c, v′⟩) is ≤ 1/(2H²‖v ∧ v′‖) for continuous cubes (Ball's slicing bound).
- But the coefficients are integers, so lattice, index and congruence effects must be handled.
- Arithmetic restrictions survive.
- The target set becomes field-dependent.
- The splitting density needs an S_d Galois group and uniform Chebotarev.

**Why this is still a different problem from Conjecture F [assessment].** Conjecture F needs equidistribution of the
Bockstein half of a hidden 2-power character, which is already sub-Burgess for QS and needs a least-nonresidue bound for
degree 1 (§5). Hypothesis D needs smooth-value statements of LV's kind, but per-target and jointly. It needs no hidden
2-adic data. It is neither shown easier nor shown harder; it moves the difficulty from character theory to multiplicative
number theory.

## 4. Averaging over N (angle a)

**4.1 Decoupling parametrizations [proved; elementary].**
- **QS.** N = x² − s and (x/N) = (−s/x) for odd x and Blum N (codex (8)). The weight becomes 1_{P⁺(s) ≤ B}. The second
  moment over N needs 1_sm(s)·1_sm(s + x′² − x²), a binary problem (codex (10)).
- **LV.** For fixed m, (N, R) ↦ f is injective (codex (11)). The semiprime condition f(m) ∈ 𝓕 remains.
- **Degree 1.** The determinant couples the coefficients (codex (12)).

**4.2 Lengths [proved, codex lemmaM_1 (5)].**
- Σ_m |Σ_t a_t(t/m)|² ≪ min(Q² + H, Q + H²)Σ|a_t|² (squarefree restrictions; absolute constant). So for Q ≥ H² there is no
  ε-loss. Heath-Brown's (QH)^ε matters only near Q ≍ H.
- Liu (arXiv 2505.09637; cited by codex, unverified by us) makes the ε-dependence explicit, exp₄(C/ε).
- The obstruction is the N-dependence of the weights and the adaptive test w, not the ε.

**4.3 Results and limits.**
- **Theorem A (codex lemmaM_1 §5) [proved; checked by us and by the critic].** For fixed weights on odd t ≤ H < √N and a
  balanced Blum family of size ≍ N/(log N)², the exceptional proportion is ≪ (log N)²E/(η²A²). In particular, Jacobi signs
  on the odd L[1/3]-smooth t ≤ L[2/3] are balanced for almost all n.
- **Theorem B (codex) [proved].** With t a uniform prime ≤ B, x(t) = (e_t, 1) and β_n(t) = [(t/n) = −1], every fixed test
  cancels for almost all n, yet min_w δ_w = 0 for every n.
- **Counting observation [observation].** The LV parameter space has size N^{δκ−1+o(1)}, which is ≪ 2^C. Bounds on bad
  sets for each fixed test cannot be summed over tests.
- **Fibre mixing (codex (28)–(30)) [proved].** δ ≥ η* := Σ_x min(Pr[x, β = 0], Pr[x, β = 1]). This is stronger than Lemma M.
  On descent relations DD supplies η* = 1/2 conditionally.
- **QS fixed-test first moment [sketch].** Via N = x² − s, reciprocity and (5) with moduli ≍ √X and length ≍ √X·H, one gets
  Σ_N Z_w(N) = o(Σ_N |X_N|) when Hρ → ∞. This gives only an average δ_w ≥ 1/2 − o(1) for each w, which is useless for Lemma M
  by Theorem B.

**Verdict on (a):** refuted for this model, i.e. fixed-test averaging. A route via averaging needs an inverse theorem or
another way to control the adaptive predictor. DD moves the factoring question elsewhere; it does not resolve this one.

## 5. Involutions (angle b) and barriers

- **Unit twist [proved, classical; E2].**
  - In the continued-fraction relations of √N, k ↦ k + P is multiplication by ε = a + b√N, the fundamental unit of
    Z[√N] (the Pell order).
  - It preserves the Q_k and multiplies the root by a mod N, with a² ≡ (−1)^P. Signed norms change when P is odd.
  - For Blum N, a is always a nontrivial √1. The argument (codex, critic): there is no negative Pell solution since some
    p ≡ 3 (mod 4); if a ≡ ±1 (mod N), factoring (a − 1)(a + 1) = Nb² gives a smaller positive solution or a negative one.
    E2: 80/80. The critic also checked 780/780 Blum N with p, q < 400.
  - Other classes: 58/80 nontrivial for p ≡ 1, q ≡ 3 (mod 4); for p ≡ q ≡ 1 (mod 4), 31 nontrivial and 49 of norm −1.
  - Displacement is one period, P ≈ N^{1/2+o(1)}/h typically (fit at 28–40 bits: slope 0.511; not an asymptotic proof).
- **Square-twist identity (codex (15)) [proved].** If q′ = s²q and α′ = t²α with s, t explicit, the root ratio of the pair
  is s/φ(t), so a β-flipping pair gives gcd(s − φ(t), N). Consequence: an efficiently constructible β-flipping involution is
  a factoring algorithm. Otherwise the symmetry must act on the sampling law, as in Dixon and DD.
- **Galois / QS conjugation / degree-1 GL₂ [proved, codex].**
  - Odd-degree K has no order-2 automorphism.
  - b ↦ −b preserves the modular root.
  - Primitive odd-determinant form pairs have only M = ±I.
- **Möbius symmetries [proved, ours; corrected by the critic].**
  - Suppose γ = (p q; r s) ∈ GL₂(Q) satisfies F∘γ = λF and fixes m ∈ P¹(Q), so the rational form satisfies G∘γ = λ′G.
  - Then (a′ − b′α) = (p − rα)·σ(a − bα), with σ(α) = γ^{−1}·α ∈ Aut(K), φ∘σ = φ, and λ′ = φ(p − rα).
  - For even |S| the root ratio is preserved exactly. For odd |S| it changes by an explicit constant in μ₂, which by (15)
    is itself a factorisation.
  - σ permutes the ideal coordinates of x, so this is not an instance of criterion (I).
- **MNFS [assessment, codex §3.5].** No identification of exact square classes across fields; no opposite-label proof.
- **(B1) Least-nonresidue barrier [proved; scope: even-|S| two-form degree-1 sieves with no known root of λ, i.e.
  nfsfruit Cor. J2].**
  - For Blum N, Φ(S) = Σ_S ⌊Ω_T(L1L2)/2⌋, where T is the set of primes with (ℓ/N) = −1.
  - If no prime ℓ ≤ B has (ℓ/N) = −1, then T ∩ FB = ∅ and every congruence is trivial.
  - Hence Conjecture F for this mechanism for all Blum N implies n_J(N) ≤ B = L[1/2]. Unconditionally only
    N^{1/(4√e)+ε} is known (Burgess–Vinogradov, primitive (·/N)); under GRH, 2 log² N (Bach).
  - The barrier does **not** apply to hidden-root two-form sieves (nfsfruit L1), nor to QS in general (N = 33: 7² − 33 = 2⁴,
    (2/33) = +1, yet gcd(7 − 4, 33) = 3).
  - For almost all N in natural families the least nonresidue is small. The averaging family must be stated (Pollack's
    results); we do not claim a specific family here.

## 6. Class groups and the r18 rejections (angle c)
- **[proved in the model, codex (16)].** Uniform-class-sampling relation generation with norms N^{1/2+o(1)} costs
  max(b + 1/(4b), 2b) ≥ 1. Lenstra–Pomerance remark the same limitation of their Algorithm 10.1.
- **Geometric obstruction (codex (17)).** Nonrational integers of norm < D/4 do not exist in an imaginary quadratic field
  of discriminant −D.
- **Counting obstruction (codex (19)–(20)).** A nearly uniform class sampler finds representatives of norm ≤ D^α
  (α < 1/2) with probability ≤ D^{α−1/2+o(1)}.
- So LV-style randomization does not lower c for class groups. We agree with r18.
- **On r18's rejection of uniform hidden square multipliers.** It is correct for one-step smoothness of the square, or for
  small ratio representatives. Descent instead splits the square into two √N-halves that are B1-smooth with B1 ≫ B, and
  pays special-q costs, as NFS-DL individual-log descent does. This is a genuine change in relation supply,
  heuristically. Its rigorous status is Hypothesis D (§3), which is open.
- **Modern inputs.** Pascadi's 5/8, Harper's BV, Matomäki–Radziwiłł and Jain's short intervals bear on the first moment of
  Hypothesis D. None of them reaches the needed moduli directly.

## 7. Experiments (pre-registered in `e529d50`; exact arithmetic)

Commands:
```
.venv/bin/python -m src.lemmaM.harness_e1 experiments/configs/lemmaM_e1.json experiments/raw/lemmaM_e1.jsonl 8
.venv/bin/python -m src.lemmaM.analyze e1 experiments/raw/lemmaM_e1.jsonl
.venv/bin/python -m src.lemmaM.run_e2 experiments/configs/lemmaM_e2.json experiments/raw/lemmaM_e2.jsonl 4
.venv/bin/python -m src.lemmaM.analyze e2 experiments/raw/lemmaM_e2.jsonl
```

**E1** (toy QS-descent; 4 classes × {48, 56, 64, 72} bits × 20 seeds × {pure, mixed, plain}; 960 records). The numbers were
re-verified by the critic.

| prediction | result |
|---|---|
| P1: X² ≡ Y² always | 0 non-squares |
| P2 pure, semiprimes, 1/2 | 7165/14400 = 0.4976, cluster z = −0.63 (blum 0.493, mixed 0.504, both1 0.496) |
| P2 pure, 3-prime, 3/4 | 3498/4740 = 0.738, z = −1.73 |
| P3 mixed, descent-containing | semiprimes 0.504 (z = 0.83); 3-prime 0.756 (z = 0.87) |
| mixed, QS-only (reported) | 0.501 (z = 0.34); 0.745 (z = −0.97) |
| P4 flip check | 32 092/32 092 |
| P5 r-trials per relation, plain/descent ≥ 5 | 9.5, 13.2, 17.0, 22.7 (48 to 72 bits) |

- **P4 is a determinism and bookkeeping check, not evidence for DD.** Since (εr)² = r², it confirms that
  `split_and_descend` is a stateless function of R and that X_S·ε flips the class. The probabilistic content of DD is tested
  only by P2 and P3.
- **P5 counts r-trials only.** Counting all smoothness tests (r-trials + special-q tests), the plain/descent ratio is
  1.13 / 1.68 / 2.39 / 2.34 at 48 / 56 / 64 / 72 bits. Special-q relations are cached across descents, which is harmless for
  DD.
- Special-q searches succeeded 100% (1.3–1.7 distinct ℓ per descent).
- Two 3-prime 48-bit instances hit an ACCIDENTAL factor: a split prime equal to a 16-bit factor of N, at the null rate
  ≈ #splits·2/p.
- The toy's special-q values are QS-sized. **No complexity claim** is made.

**E2** (unit twist; {blum, mixed, both1} × {28, 32, 36, 40} bits × 20 seeds).
- P6: identities in 240/240.
- P7: slope 0.511 ∈ [0.40, 0.55].
- P8: blum 80/80 nontrivial (predicted ≥ 0.8; it is in fact a theorem); mixed 58 nontrivial / 22 trivial; both1
  31 nontrivial / 49 norm −1.

## 8. Candidate cards

```text
ID: DD-1 Dixon descent (hidden root injected through the coins).
Target: every N not a prime power.
Operation: r uniform, R = r^2, split R = u/v, special-q descent of large primes, relation (rW)^2 = ev(c).
Why it exposes a factor: Theorem DD (proved, folklore-level): Pr = 1 - 2^{1-omega} per dependency, under transcript invariance.
New lemma required: Hypothesis D = (first moment at special-q moduli) + (complete-target-list average) + (split counting).
Cost (heuristic / conditional on Hyp. D): L[1/2, 2/sqrt(d+1)] formal; L[1/2, sqrt(6/(d+1))] level-1/2 design (c<1 for odd d>=7);
  L[2/5] one-level (Gordon-type); no rigorous L[1/3] with HECM splitting (L[4/9] barrier).
Closest prior: Rabin; Dixon/Pomerance; nfsfruit L1; Stange 2022; NFS-DL descent (Gordon 1993, Joux-Lercier, Commeine-Semaev).
Status: fruitfulness PROVED; relation supply OPEN (not near known tools).

ID: DD-2 QS + descent (QR filter on split primes). Heuristic L[1/2,1], provably fruitful; removes H2 of below16 O-4 (randomized).
Status: proved (fruitfulness), heuristic (cost); E1 toy check passed.

ID: IV-1 unit-twist involution (QS/CFRAC). Exists; beta-flipping for every Blum N (proved); displacement one period (~N^{1/2}/h).
Status: refuted for this model at subexponential scale.

ID: AV-1 averaging over N. Loss-free large sieve; fixed-test cancellation for almost all N (codex Thm A).
Status: refuted for this model as a route to Lemma M (codex Thm B).

ID: B1 least-nonresidue barrier for even-|S| two-form degree-1 sieves. Status: proved (scope-restricted).
```

## 9. Reviews
- **GPT-6.1-sol lemmaM_1 (blind).** It supplied Theorems A and B, the loss-free bound (5), the odd-degree, conjugation, GL₂
  and square-twist results, the class-group obstructions, the literature survey (Pascadi, de Boer–Pellet-Mary–Wesolowski,
  Al Aswad–Pierrot, Barbulescu–Jouve, Liu), and the fibre-mixing condition. It did **not** find the descent route. It left
  "a generator that chooses its own smooth target with a fibre-uniform root" as the open possibility, and DD answers that
  possibility heuristically.
- **GPT-6.1-sol lemmaM_2 (adversarial).**
  - DD: CORRECT WITH FIXES (transcript invariance, per-dependency).
  - (a) arithmetic: CORRECT WITH FIXES (non-monic height, lattice regularity, C·L^{o(1)} rows).
  - (b) L[2/5]: CORRECT WITH FIXES (log factors; Gordon prior art; the HECM 2β₁/3 constraint).
  - QS analogue: WRONG as stated. It is fixed by the QR filter, which the code already had.
  - Lemma SQ ⇒ Hyp. D: WRONG, replaced by the complete-target-list statement.
  - "First moment of LV type": WRONG (modulus range), fixed via (a′) and §3.
  - Two-sided repair: WRONG as a completed reduction, now labelled a research programme.
  - B1: CORRECT WITH FIXES (scope).
  - B2: CORRECT WITH FIXES, with the elementary Blum proof added.
- **Independent Opus critic: ACCEPT WITH MAJOR FIXES.** All applied:
  - DD hypothesis, novelty wording and prior art (Rabin, L1, Gordon);
  - (c) corrected from "L[2/3]" to L[1/2] at the re-optimised degree, plus the L[4/9] HECM barrier;
  - Corollary DD row count and height adaptation;
  - Verdict wording ("first moment of LV type", "first-order", "stochastic deepening");
  - split square-class condition;
  - B1 scope (conflict with L1 and QS);
  - "moot" deleted;
  - unit twist stated as a theorem;
  - Möbius formula and the odd-|S| caveat;
  - P4 described honestly and the P5 total-test ratio added;
  - labels ("observation", "sketch").
  Its scratch checks (780 Blum Pell cases; LPP HECM statement) are in the session scratchpad, not committed.
