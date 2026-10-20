# inversesieve: can DCL non-coverage be decoded into cheap containers? (r20 idea 2)

Agent `inversesieve`, phase 13, 2026-10-20.
- Pre-registration: `notes/inversesieve_prereg.md`, commit ee39c2a (before any registered run).
- Code: `src/inversesieve/{design,usets,stats,theory,experiments,analyze}.py`; exploratory scripts
  `src/inversesieve/explore_*.py`. Tests: `tests/test_inversesieve.py`.
- Raw data: `experiments/raw/inversesieve_e{1,2,3,4}.jsonl`, `inversesieve_x_replicate.json`.
  Summary: `experiments/analysis/inversesieve_summary.json`.
- GPT-6.1-sol: blind round `codex/inversesieve_1.out`; adversarial referee `codex/inversesieve_2.out`.
  Opus critics: §9.

## 0. Verdict (plain language)

**Kill. Non-coverage does not carry usable structure, and no deterministic N^{6/31} algorithm comes out of
this route.** The idea fails for a structural reason, not for lack of a clever lemma.

- **The main insight: the information is too diffuse.**
  - The event "the true pair is uncovered" is extremely unlikely: probability ≈ e^{−λ}, λ = N^{1/31}.
  - That improbability is spread over J ≈ N^{5/31} band primes. Each prime contributes only
    ρ ≈ N^{−4/31} nats.
  - A Coppersmith container can absorb a modulus of size ≲ N^{1/2}, i.e. about 3 band primes. Those carry
    o(1) information.
  - To concentrate enough information to cover the survivors, one must combine the conditions of
    ≳ J^{4/5} = N^{4/31} primes jointly (Prop D). That is the coverage problem again, in disguise.
  - With boundedly many primes, an inverse theorem exists only at Minkowski height. There it reproduces
    Harvey's 1/5 exactly (Cor. D2).
- **What the inverse theorem really says (Prop C, proved).** "Uncovered ⇒ p ∈ C(N)" is the same
  statement as "every pair outside C(N) is covered".
  - The containers have total mass ≤ N^{6/31}·N^{1/4}. So the statement is a coverage theorem at
    λ = J^{1/5} for all pairs except an *explicitly enumerable* set of density ≤ N^{−7/124+o(1)}.
  - No known argument produces explicit exceptional sets. The best proved exceptional density is not
    explicit, and it is a further factor N^{3/124} too large (Prop M).
- **Containers that ignore most of N are provably useless (Theorem L, proved).**
  - Take a container rule that sees N only through a residue N mod M and a size window (cells with
    Δ/M ≥ P^{55/62+δ}), and must hold for pairs uncovered at the primes it sees.
  - Such a rule needs on average ≳ P·min(1, Δ/(PM))/n² candidate mass, i.e. ≥ N^{1/4}/n² Coppersmith
    containers when Δ ≥ PM.
  - To fit the N^{6/31} budget, a rule must use more than 69/124 ≈ 56% of the bits of N in a non-forgetful
    way. That means the exact hyperbola pq = N, which is what codex's kill test excludes as a mere
    rediscovery.
- **The literature does not apply (Prop I).** Every inverse-sieve theorem (Green–Harper,
  Helfgott–Venkatesh, Walsh, Gallagher/Croot–Elsholtz, Croot–Mao–Pohoata–Yip 2026) concerns sets
  occupying *few* residue classes. U(N) occupies all but a fraction N^{−4/31} of the classes. For the
  singleton {p} the theorems apply only vacuously: the "curve" is x − p = 0.
- **Under the random model the exceptions do not even exist (Prop E, heuristic).** At ϑ = 1/5 the expected
  number of uncovered pairs, and of uncovered integer candidates over a dyadic range, is
  N^{O(1)}e^{−N^{1/31}} → 0. There is no population of exceptions whose structure could be exploited.
- **Experiments (pre-registered; 1000 exact U(N) at n = 24–40; 1842 uncovered true pairs; all 172
  dchworst worst pairs).** U(N) and the uncovered true factors are pseudorandom in every tested container
  family:
  - |U|/model is 0.98–1.01;
  - interval covering numbers at length N^{1/4} are within 0.7% of random;
  - Coppersmith-AP triples are 0.985–1.001 of random;
  - class occupancy mod small primes shows no low-degree-image signature;
  - additive energy is 0.99–1.00 of random;
  - Lehman containers: RR 1.02, CI [0.87, 1.20];
  - p's position inside U(N) is generic (KS p = 0.13);
  - worst-pair L* percentile has mean 0.507 (KS p = 0.31);
  - in all 88 large-J worst pairs, U(N) = {p} exactly. The uncovered set is then a perfect container in
    the information sense, but finding its single point is the sieve itself.

  All frozen predictions in the structured (H1) direction pass. The two failures are in the *dispersed*
  direction (§7).
- **A side finding (experimentally supported; mechanism UNVERIFIED).** At bands with Y ≈ P/16, true pairs
  are uncovered *less* often than the independent-residue model predicts: 0.76–0.87 of it. The cause is
  negative dependence between band halves. Coverage is better than random there, which matches dchworst
  P3. It does not help containers.
- **One weak positive (proved, Theorem S).** The exact forbidden sets are R_ℓ(N) = √(N·F_V(ℓ)).
  - For sign-symmetric weights, |R_ℓ(N)| = |V| for every N.
  - E×(R_ℓ(N)) = E×(F_V(ℓ)) exactly.
  - F_V(ℓ) is additively structured: a union of H dilated intervals, with 1.4–4.5× random symmetric
    energy. The square-root map erases that: E(R) is 0.99× random symmetric energy.

## 1. Setting

- Notation is that of dchworst §1 and coverage38 §1. N = pq, p < q ≤ 2p, P = √N.
- The band 𝒫 is J primes ℓ ≡ 3 (mod 4) in (Y, 2Y]. Weights V are primitive (A, B) of height ≤ H.
  λ = |V|Σ1/ℓ, and ρ_ℓ = |R_ℓ(N)|/ℓ.
- F_V(ℓ) = {AB^{−1} mod ℓ}. R_ℓ(N) = {x : Bx² ≡ AN (mod ℓ) for some (A,B) ∈ V}.
- U(N) = {x ∈ I : x mod ℓ ∉ R_ℓ(N) ∀ℓ}, with I = [⌈√(N/2)⌉, ⌊√N⌋].
- The pair is covered iff p mod ℓ ∈ R_ℓ(N) for some ℓ (Bp² − AN = p(Bp − Aq)). So p ∈ U(N) ⇔ uncovered
  (`test_U_matches_bruteforce_and_true_factor`).
- **Ledger (exact; `theory.ledger`).** The attaining single band at λ = J^ϑ (coverage38 U2) has
  X = N^{(1+ϑ)/(6+ϑ)}. At ϑ = 1/5:
  - X = N^{6/31} = P^{12/31};
  - Y ≍ J ≍ P^{10/31};
  - λ = P^{2/31} = N^{1/31};
  - H² ≍ λ (up to logs);
  - ρ ≍ λ/J = N^{−4/31};
  - Coppersmith containers hold ≤ P^{1/2} integers, so the container mass fraction is
    X·P^{1/2}/P = P^{−7/62} = N^{−7/124};
  - log P/log Y = 31/10, so a modulus M ≤ P holds ≤ 3 band primes.

## 2. Theorem S: exact forbidden-residue structure (proved)

Let ℓ ≡ 3 (mod 4) be prime with ℓ ∤ N and ℓ > 2H². Let V be sign-symmetric:
(A,B) ∈ V ⇔ (−A,B) ∈ V, including (±1,1).

- **(a)** F_V(ℓ) = ∪_{B≤H} B^{−1}·{A : |A| ≤ H, gcd(A,B) = 1}. This is a union of H dilates of pieces of
  [−H, H]. Its elements are distinct, because |AB′ − A′B| ≤ 2H² < ℓ. So |F| = |V|, and F = −F = F^{−1}.
- **(b)** x ∈ R_ℓ(N) ⇔ x ≢ 0 and x²N^{−1} ∈ F. Since −1 is a non-residue, exactly one of ±Nf is a square
  for each f, which gives 2 roots. So **|R_ℓ(N)| = |V| for every N.** The local density
  ∏(1 − |V|/ℓ) is therefore N-independent. (`test_theorem_S_exact_size_signed`.)
- **(c)** Coverage of the true pair at ℓ ⇔ t_ℓ = pq^{−1} mod ℓ ∈ F. Equivalently, the lattice
  Λ_ℓ = {(A,B): ℓ | Bp − Aq} has a nonzero vector in the height-H box.
- **(d)** E×(R_ℓ(N)) = E×(F).
  - *Proof.* Let F′ = {f ∈ F : Nf is a square}. Then F = F′ ⊔ (−F′) and R = {±√(Nf) : f ∈ F′}.
  - A solution of x₁x₂ = x₃x₄ in R squares to f₁f₂ = f₃f₄ in F′.
  - Conversely, each such quadruple lifts to 8 of the 16 sign patterns. The same 8-to-1 count holds for
    F = ±F′.
  - (`test_theorem_S_mult_energy_identity`; all 250 E2 cases.)
- **(e) Measured.** E(F)/E(random symmetric set) has median 1.52, range 1.08–4.65. It grows with ℓ/H⁴,
  reaching 4.46 at 40S: the Farey set's rational structure. E(R)/E(random symmetric) has median 0.989,
  range 0.87–1.05.
  - Heuristic order: E(F) ≍ H⁴ log H when 4H⁴ < ℓ, since sums of Farey fractions then hold in Q.
- **(f) Consequence.** The structure of the forbidden slopes is multiplicative and Farey-like. Squaring
  transports it to the candidate side only multiplicatively. Additively, the classes p must avoid look
  random.

## 3. Logical content and moment bookkeeping

**Prop C (proved, trivial; codex agrees).** For any N-computable C(N), "uncovered ⇒ p ∈ C(N)" ⇔ "p ∉ C(N)
⇒ covered". The inverse theorem is a coverage theorem off the containers, plus Coppersmith inside them.

**Prop E (heuristic; Poisson model).**
- At ϑ = 1/5, E#{uncovered pairs in S_n} ≍ (P²/n²)e^{−Λ} with Λ = N^{1/31+o(1)}.
- Likewise E Σ_{N∈[T,2T)}|U(N)| ≍ T^{3/2}e^{−Λ}. Both tend to 0.
- The design is heuristically a *complete* covering. The inverse theorem would then be true for empty
  reasons. Its whole content is the proof.
- Small-n experiments with nonempty U(N) can test only *conditional* structure.

**Prop M (proved arithmetic).**
- Containers hold ≤ XP^{1/2} = N^{55/124} candidates per N, a fraction N^{−7/124} of I.
- The large sieve (D′ Lemma 9) gives |U(N)| ≤ (P + 4Y²)/λ, i.e. exceptional density N^{−1/31} =
  N^{−4/124}. That density is not explicit.
- A Poisson-quality k-th moment would give density λ^{−k/2} = N^{−k/62}. Reaching N^{−7/124} needs
  k ≥ 7/2, i.e. a 4th moment at Poisson quality (`theory.moment_order_needed` = 7/2). The resulting set
  {M = 0} is still not explicit.
- If the theorem covered *all* of U(N) (not only p), it would imply |U(N)| ≤ N^{55/124−ε}. That beats the
  large sieve's N^{29/62} by N^{3/124+ε} (codex).

## 4. Theorem L: forgetful container rules (proved)

**Model.**
- A rule is a fixed finite family of sub-rules i. Sub-rule i has an observation (M_i, Δ_i), where every
  prime factor of M_i is < P/2. It outputs a set C_i(D) ⊂ [P/2, P] depending on N only through the cell
  D = (N mod M_i, ⌊N/Δ_i⌋).
- WLOG Δ_i is a power of 2 with 2M_i ≤ Δ_i ≤ 2^{n−1}. Refining cells only gives the rule more
  information, at a loss ≤ 2 in Δ.
- Mass: m(N) = Σ_i |C_i(D_i(N))|.
- Adaptive selection among sub-rules using other data of N is *outside* the model (codex's caveat).

**Theorem L.** Suppose that for every (p,q) ∈ 𝒬 ⊂ S_n, p ∈ ∪_i C_i(D_i(pq)). Then

  2^{−(n−1)} Σ_{N n-bit} m(N) ≥ (2|𝒬|/P²)·min_i min(P/8, Δ_i/(4M_i)).

*Proof.*
- (i) Fix i, a cell D and x ∈ C_i(D). Then x is coprime to M_i (p prime > factors of M_i; drop the other
  x). The q with xq ≡ d (mod M_i) and xq ∈ [kΔ_i, (k+1)Δ_i) lie in one class mod M_i inside an interval of
  length Δ_i/x ≤ 2Δ_i/P. So there are ≤ 2Δ_i/(PM_i) + 1 of them (`test_theorem_L_partner_bound`).
  - Summing: |𝒬| ≤ Σ_i Σ_D |C_i(D)|·(2Δ_i/(PM_i) + 1).
- (ii) Every cell (aligned, interior) contains ≥ ⌊Δ_i/M_i⌋ ≥ Δ_i/(2M_i) n-bit integers N. So
  Σ_N m(N) ≥ Σ_i Σ_D |C_i(D)|·Δ_i/(2M_i).
- (iii) The ratio of the per-cell weights in (ii) and (i) is
  (Δ/(2M))/(2Δ/(PM) + 1) ≥ min(P/8, Δ/(4M)).
- Combine, and divide by 2^{n−1} = P²/2. ∎

`test_theorem_L_toy_rule_mass` checks the extremal rule exactly on S_20.

**Corollaries.**
- **L1 (local inverse theorems).** Take 𝒬 = pairs uncovered at the band primes dividing ∏M_i. If those
  are o(Y/|V|) primes, Theorem N (dchworst §4) gives |𝒬| = (1 − o(1))|S_n|.
  - If Δ_i ≥ PM_i, the average mass is ≥ cP/n². So there are ≥ cN^{1/4}/n² Coppersmith containers on
    average.
  - Theorem-G boxes need worst-case cost ≥ √(max m) ≥ cN^{1/4}/n. A value < N has ≤ 2 prime factors
    in [P/2, P], and Σ√V_b ≥ √(ΣV_b).
  - Such rules are no better than brute-force Coppersmith over all N^{1/4} windows.
- **L2 (bit budget).** For a rule with Δ/M = P^a, the average container count is
  ≥ P^{a−1/2}/n^{O(1)}. Staying ≤ X = P^{12/31} forces a ≤ 55/62.
  - So each sub-rule must observe ≥ (69/62)·log₂P bits of N in AP-cell form, i.e. 69/124 of all n bits
    (`theory.max_forgetful_modulus_exponent` = 69/62).
  - Lehman/stationary-point containers observe only the top ≈ n/4 bits (Δ ≈ P^{3/2}). They are covered by
    L1, consistent with their known catch rate ≈ (#containers)·P^{−1/2}.
- **L3 (local triggers reduce to coverage).** Now let validity be required only on globally uncovered
  pairs, with mass ≤ XP^{1/2}.
  - Combine (i) with (ii) in the form Σ_D|C_i(D)| ≤ (2M_i/Δ_i)Σ_N m_i(N). The rule then catches
    #caught ≤ Σ_i (2M_i/Δ_i)(2Δ_i/(PM_i) + 1)·Σ_N m_i(N) ≤ (8/P)·(P²/2)·max_N m(N) = 4P·max_N m(N)
    pairs (when Δ_i ≥ PM_i/2).
  - With max m ≤ XP^{1/2}n^{O(1)} this is ≤ n^{O(1)}XP^{3/2} pairs, i.e. density n^{O(1)}N^{−7/124} in
    S_n.
  - So its validity is *equivalent* to coverage of every pair outside an explicit forgetful set of that
    density. Dichotomy designs ("covered or Lehman-easy") are a special case. Theorem U, B1 and Prop L of
    coverage38 are about certifying M ≥ 1 for the remaining pairs, and are unchanged.
- **Scope.** Theorem L is a volume argument of the same kind as below16 T1/T2, applied to containers
  rather than selectors. It is not a factoring lower bound.
  - It says nothing about rules that use N non-forgetfully. Such rules use the exact integer hyperbola,
    the only object that distinguishes p from the "virtual" members of U(N).
  - Those virtual members sit on the modular hyperbolas xy ≡ N (mod ℓ) for every ℓ, but on no common
    integer hyperbola.

## 5. Diffuseness and the Minkowski endpoint

**Prop D (exhaustive local covering; proved; quantitative form of codex's bounded-k statement).**
- For S ⊂ 𝒫, let A_S = {x ∈ I : x mod ℓ ∉ R_ℓ(N) ∀ℓ ∈ S}. Then
  |A_S| ≥ |I|(1 − Σ_{ℓ∈S}|V|/ℓ) − |S||V|.
- Any family of intervals of length w (or residue classes of modulus ≥ |I|/w) that covers A_S therefore
  has ≥ |A_S|/w members.
- Mass ≤ XP^{1/2} = |I|·N^{−7/124+o(1)} forces Σ_{ℓ∈S}ρ_ℓ ≥ 1 − o(1). That means
  |S| ≥ (1 − o(1))Y/|V| ≍ J/λ = J^{4/5} = N^{4/31+o(1)} primes, with modulus ∏_S ℓ ≥ exp(N^{4/31+o(1)}).
- In the independent model, conditioning on global avoidance moves the law of any bounded set of
  coordinates by total variation ≤ Σ_Sρ_ℓ = o(1) (codex; proved there).

**Cor. D2 (where inverse theorems do exist).** Let |S| = s be bounded. Then Prop D needs
ρ ≥ (1 − o(1))/s per prime, i.e. H² ≳ Y/s.
- This is Minkowski height. Example (dchworst P1, E4 P9b): with (1,1) excluded and H = ⌊√ℓ⌋,
  "uncovered at ℓ ⇒ t_ℓ = 1 ⇒ p ≡ ±√N (mod ℓ)".
- Two such primes give 4 classes mod ℓ₁ℓ₂ > P^{1/2}, which is Coppersmith-decodable.
- The design cost at H ≍ √Y is HP/Y² + sY. Its minimum is at Y = P^{2/5}, giving **N^{1/5}**: Theorem W
  at c = 0. Inverse theorems with boundedly many primes therefore live exactly at Harvey's exponent.

**Prop Q (planted-survivor query model; codex, re-checked).**
- Hide z uniformly in I. An algorithm whose public data are independent of z may ask Q "is x = z?"
  queries, then output containers holding B integers.
- Then Pr[success] ≤ (Q + B)/|I|.
- With Q = N^{1/5} and B = N^{1/5+1/4}, this is ≤ N^{−1/20+o(1)}.
- This is a model statement. Transferring it to arithmetic N would need exactly the pseudorandomness that
  §7 supports but nobody can prove.

## 6. Literature (Prop I; KNOWN; checked by codex and the literature critic, §9)

| result | hypothesis | met by U(N) or {p}? |
|---|---|---|
| Green–Harper 1311.6176 Thm 1.1/1.3 | A mod ℓ ⊂ S_ℓ, \|S_ℓ\| = (ℓ+1)/2 with energy ≥ (1/16+δ)ℓ³ for all ℓ ≤ T^{1/2}; or APs of length (1−ε)ℓ | no: U keeps (1−ρ)ℓ classes on one band |
| Helfgott–Venkatesh (2009) | \|S mod ℓ\| ≤ αℓ, S ⊂ [0,T]² | only vacuously ({p}, or the known hyperbola) |
| Walsh (Duke 2012; 1307.0259) | \|S mod ℓ\| ≤ αℓ^k, or ≪ ℓ^κ; prime mass ≫ log Q | no; for {p} the polynomial is x − p |
| Gallagher larger sieve; Croot–Elsholtz | \|A mod ℓ\| ≤ ν_ℓ with Σ log ℓ/ν_ℓ > log L | denominator is negative or trivial |
| Croot–Mao–Pohoata–Yip 2607.15311 Prop 5.1 | nonempty E, band product ≤ √T/2 | applies, but the gain is e^{O(√(ρ/ℓ))} = e^{o(1)}, non-constructive |

All of these concern sets confined to few classes, the large-sieve-extremal regime. Non-coverage is the
opposite regime: few classes are *forbidden*. A survivor of a sieve that removes a fraction ρ → 0 per prime
is not forced into any structure. Random sets survive at the same rate.

## 7. Experiments vs frozen predictions (`experiments/analysis/inversesieve_summary.json`)

**E1** (100 N per config; exact U(N); frozen nulls). The table gives P1 |U|/model; P2 the interval cover
ratio at L = N^{1/4} with its z-score; P3 Coppersmith-AP triples; P4 the χ² occupancy ratio; P5 the pair
ratio vs the local model; P5b the halves ratio; P6 integer energy; P8c the Lehman ratio.

| cfg | \|U\| mean | P1 | P2 (z) | P3 | P4 χ² | P5 | P5b | P6 | P8c |
|---|---|---|---|---|---|---|---|---|---|
| 24M | 193 | 0.997 | 1.000 (0.2) | 0.985 | **0.80** | 0.994 | 0.997 | 1.004 | 1.00 |
| 24S | 102 | 0.998 | 1.001 (0.3) | 0.987 | **0.87** | 0.993 | 1.000 | 0.992 | 1.03 |
| 28M | 242 | 1.006 | 1.007 (3.5) | 0.987 | **0.87** | 1.013 | 1.008 | 1.003 | 0.99 |
| 28S | 132 | 1.012 | 0.997 (−0.9) | 0.998 | 0.97 | 1.029 | 1.010 | 1.000 | 0.99 |
| 32M | 171 | 0.983 | 1.003 (0.8) | 0.995 | 0.96 | 0.989 | 0.986 | 0.999 | 0.98 |
| 32S | 210 | 1.002 | 0.999 (−0.3) | 1.001 | 0.99 | 1.009 | 1.002 | 0.997 | 0.98 |
| 36M | 115 | 1.003 | 0.999 (−0.2) | 0.990 | 0.97 | 0.983 | 1.003 | 1.000 | 1.04 |
| 36S | 218 | 0.994 | 1.001 (0.3) | 1.001 | 0.95 | 0.984 | 0.994 | 1.000 | 1.00 |
| 40M | 167 | 1.012 | 1.006 (1.5) | 0.988 | 0.99 | 1.056 | 1.011 | 1.000 | 0.98 |
| 40S | 240 | 0.999 | 0.996 (−1.1) | 0.997 | 0.97 | 1.038 | 0.999 | 1.000 | 1.00 |

- P1, P2, P3, P5 (gap KS p ≥ 0.077), P5b, P6 and P8c pass in all 10 configs.
- No U(N) was empty, and no missing classes mod r ≤ 23 occurred anywhere. The survivor trajectory tracks
  the model within 1.3% at every checkpoint.
- **P4 fails in 24M, 24S and 28M (χ² ratio 0.80–0.87 < 0.9).** U(N) is *more* uniform mod small primes
  than random: the dispersed direction.
  - Explanation (proved mechanism, magnitude heuristic): U = I minus a union of APs of difference ℓ. Each
    such AP is spread over the classes mod r as evenly as possible.
  - The effect vanishes as the removed fraction falls (0.95–0.99 at n ≥ 32).
- **P4b fails in 40S (CRT collision ratio 0.78; 34 vs 44 expected, −1.5σ).** The tolerance was too tight
  for ~44 counts. Dispersed direction.
- |U(N)| is over-dispersed across N (var/mean 1.5–3.5). This is a count fluctuation across inputs, not
  within-set clustering: the cover numbers are random-like. Mechanism UNVERIFIED; candidate: the ±x
  symmetry mod each ℓ.

**E2.** P7 passes. Median E(R)/E_sym = 0.989 and median E(F)/E_sym = 1.52 (40S: 4.46). E×(R) = E×(F) in
250/250 cases.

**E3** (1842 uncovered true pairs, 2000 controls). **P8 passes.**
- Lehman RR 1.02 [0.87, 1.20].
- NN-percentile KS p = 0.13.
- Near-count ratio 0.955.
- Held-out ratio 0.932.
- Unthresholded reports: L* ≤ X_s holds for 100% of both groups, uninformative at these n. Dirichlet
  quality KS p = 0.68. L* KS p = 0.008 and cf-jump KS p = 0.004.
- The last two prompted the decision-rule follow-up (exploratory, disclosed):
  - With **matched controls** (same p, q a uniform prime in the same range; `explore_matched.py`): cf-jump
    p = 0.006, L* p = 0.08, Lehman 13.6% vs 11.0%.
  - A **fresh-seed replication** (`explore_replicate.py`, 3627 uncovered pairs, matched controls):
    Lehman RR 1.05 [0.92, 1.19]; mean Δlog(cf jump) = 0.021 ± 0.030; dead weights (0 < |δ_v| < Y) 0.433
    vs 0.423.
  - **Not replicated**, so no structure. Even the largest point estimate is a ~20% enrichment, whereas an
    inverse theorem needs enrichment by a factor N^{7/124} with probability 1.
- **Held-out side finding (exploratory, `explore_held_prime.py`).** In S-shapes the partner survival on
  band half 2, given survival on half 1, is 0.90 / 0.86 / 0.81 of the independent prediction at
  n = 32 / 36 / 40. The same holds for prime partners (0.90 / 0.83 / 0.81), with full-band uncovered rates
  0.87 / 0.80 / 0.76 of the model.
  - Random classes of the same size give 1.00, and U(N) halves give 1.00. So it is arithmetic and
    pair-specific.
  - The same-weight exclusion (|δ_v| < Y² ⇒ at most one band factor) explains only ≈ 0.95. Co-coverage by
    distinct weights is 0.996 of independent (`explore_pair_explore.py`). The rest is UNVERIFIED.
  - Direction: uncovered pairs are *rarer* than Poisson. This is dchworst e1-P3's "better than random
    coverage". It is irrelevant for containers, and at the real scale (Y² = P^{20/31} < P) the exclusion
    mechanism is absent.

**E4** (dchworst worst pairs). **P9 and P9b pass.**
- All 172 records are uncovered with p ∈ U(N).
- For the 88 records with J ≥ 16: L*-percentile mean 0.507 (KS p = 0.31); Lehman 9.1% vs control
  13.1% / 9.0%.
- The J ≤ 2 censored signed records have p ≡ q (mod ℓ₁): the known (1,1)-exclusion artefact, the only
  *structured* worst cases. They are exactly the Minkowski-endpoint inverse theorem of Cor. D2.
- **U(N) = {p} exactly in all 88 J ≥ 16 records.** Other survivors: 0 observed against 0.57 predicted.
  For J < 16: Σ(|U| − 1) = 40 389 against Σ model = 40 366.
  - The frozen analysis's "U_over_model = 156" is selection bias: these N were chosen so that p survives.
    Corrected for that, the other members follow the model.
  - So the NN-percentile is undefined for these records, and the frozen code returned None. P9 is
    evaluated on the L* percentile alone.
  - **This is the cleanest illustration of the kill.** At genuinely extreme N the uncovered set *is* a
    perfect container: U(N) = {p}. The obstacle is purely computational. Locating the unique survivor of
    a sieve that removes only ρ → 0 per prime needs either the sieve over I (cost ≍ P), or the
    information of ≳ J/λ primes combined (Prop D). It cannot be done from a forgetful description of N
    (Theorem L).

**Decision (frozen rule).** All H1-direction predictions pass. The failures are dispersed or tolerance
failures. Verdict: **no conditional structure detected in any tested container family; kill in the tested
families.**

## 8. Candidate card

| field | content |
|---|---|
| mechanism | DCL design at λ = J^{1/5} + inverse theorem "uncovered ⇒ explicit containers" + Coppersmith |
| payoff if true | deterministic N^{6/31} for every balanced N |
| status | **refuted for forgetful/local containers (Theorem L, Prop D, proved)**; equivalent to explicit-exception coverage (Prop C); experimentally pseudorandom (E1–E4) |
| sharpest obstruction | information diffuseness: ρ = N^{−4/31} per prime, ≤ 3 primes per Coppersmith modulus; ≥ J^{4/5} primes must be combined (Prop D); bounded-prime inverse theorems sit at Minkowski height = N^{1/5} (Cor. D2) |
| what would revive it | an explicit, non-forgetful (≥ 56% of the bits) description of a density-N^{−7/124} set outside of which coverage at λ = J^{1/5} is provable. That is a new coverage theorem, not an inverse-sieve theorem |

## 9. Reviews

(filled in below after the GPT-6.1-sol referee and two Opus critics)

## 10. Commands

```
.venv/bin/python -m pytest -q tests/test_inversesieve.py
.venv/bin/python -m src.inversesieve.experiments e1 --workers 9   # likewise e2, e3, e4
.venv/bin/python -m src.inversesieve.analyze
.venv/bin/python src/inversesieve/explore_matched.py
.venv/bin/python src/inversesieve/explore_replicate.py
.venv/bin/python src/inversesieve/explore_held_prime.py 40 S 30   # also 36 S 40, 32 S 60
.venv/bin/python src/inversesieve/explore_pair_explore.py 40 S
```
