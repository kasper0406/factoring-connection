# ptesplit: synthesis on route 1's open escape (dense N-dependent split families / multigrade PTE)

Agent `ptesplit`, phase 16, 2026-10-20. This note is revised after codex ptesplit_2 (MAJOR REVISION, all points applied)
and two Opus critics (§7).

**Files.**
- Formulation, shared with ptebarrier and ptecons: `notes/ptesplit_formulation.md` (commit 641c3e2).
- Pre-registration: `notes/ptesplit_prereg.md` (commit 987a9c7, made before E1/E2).
- Code: `src/ptesplit/{census,nccheck,run_experiment,verify}.py`.
- Tests: `tests/test_ptesplit.py`.
- Config: `experiments/configs/ptesplit.json`.
- Raw results: `experiments/raw/ptesplit_E{1,2}.jsonl`, plus the post-hoc check `ptesplit_E1_verify_posthoc.jsonl`.
- Codex: blind round `codex/ptesplit_1.{md,out}`; adversarial round `codex/ptesplit_2.{md,out}`.

**Sources synthesised.**
- `notes/ptebarrier_results.md`: final revision 4ac68dc, after its codex round and its Opus critics.
- `notes/ptecons_results.md`: final, 6ce0300.
- `notes/resbarrier_results.md`: §2.6 and §4.2.
- `notes/systeq_results.md`, `notes/smoothenum_results.md`, `notes/heurnfs_results.md`.

**Labels.** PROVED (proof here or in the cited note) / KNOWN (cited) / HEURISTIC (assumption named) / measured / open /
assessment. Novelty is never claimed.

---

## 0. Verdict

**Plain language.** The escape asked for many split ("smooth by construction") polynomials that share one small
reduction modulo an N-dependent presentation. This is the integer version of how the finite-field algorithms get about
q³ relations from one congruence.

**What is proved: the faithful, systematic version gives no amplification.**
- This is the version where the reduction is an *exact* polynomial identity with small multipliers and a small remainder,
  the analogue of the BGJT identity.
- In the stated rigidity regime, one presentation yields one relation (ptebarrier Theorem W and Corollary P).
- The reason is archimedean. Over Z, norms and heights are integers that cannot be small for many independent pieces.
  Over F_q, every constant is a unit.

**What is not proved.** Any barrier for reductions that hold only modulo N, or only in value at an integer. These remain
open, under conditional barriers only:
- ptecons' no-cancellation heuristic NC;
- the window heuristic below;
- resbarrier's G_r.

No construction meets the target. ptecons examined seven families, and our census shows fibres collapsing by M ≈ 8.

**Precise verdict on route 1's open escape.**
1. **Exact compressions: refuted, in the regime.** By ptebarrier Theorem W, Corollary P and Prop S [PROVED, with their
   degree, height and regime hypotheses], an exact compressed family whose degrees and log-heights are
   o(x/log B) provides ≤ 1 relation per presentation. This covers every BGJT-shaped family.
   - Outside the regime it is open: large multiplier or quotient degree, and corner parameters (ptebarrier O3).
   - For small M, Lemma NC gives an NFS-like *description* of the members. It does **not** give an exclusion (codex #2).
2. **Seeded additive (PTE fibres): bounded amplification.**
   - A k = 1 fibre has ≤ 2c·B^{1/2}(log B)^κ + 4 members for M ≥ 4 (ptebarrier, from Binyamini–Cluckers–Novikov).
   - Its residuals are spaced by (M − 1)! (Kleiman).
   - So one seed buys at most B^{1/2+o(1)} relations.
   - No sub-L[1/3] seed source is known. That is an assessment: we prove no seed-cost lower bound, and cheap special seeds
     exist (codex #5).
   - Near-ideal fibres (k ≥ 2) are controlled only heuristically (NC; the window lemma).
3. **Non-exact and evaluated compressions: open (ptebarrier O1).** No barrier and no candidate. Any construction would
   have to live here: congruences mod N that are coherent across a family without being exact over Z.
   - The sufficient primitive is smoothenum's output-sensitive split or doubly-smooth listing.
4. **Relation to resbarrier §2.6.** In the residue view (no rational side), the exact two-order implementation of the
   escape is a bounded-coefficient rank-≤ k linked module between two nonlinear orders: an *instance* of resbarrier
   §2.6 (§2.4). That case is closed only under the generic-modulus heuristic G_r, with a conditional formal floor of
   L[1/4] for aggregate families.
   - The containment is scoped. It covers that implementation only; other presentations and extraction schemes are not
     shown to reduce to it (codex #7).
   - The converse fails: §2.6 needs no split structure.
   - **Joint statement:** in the residue view, what remains of route 1 that is not already O1 lies inside resbarrier's
     single open case. A rigorous rank-r ledger, or a proof of G_r, would close both at once. In the integer-evaluation
     view (rational side present), ptebarrier Prop S applies in its regime. resbarrier KL bounds heights only under its
     own hypotheses (fixed module, rank ≥ 2), which we have not verified here for every member shape.

**Assessment.** Probability < 1% that this escape yields sub-L[1/3]. The other estimates:
- ptebarrier: < 1%; ptecons: < 1%;
- codex: 1–3% (blind), ~2% (ptebarrier's round);
- systeq: < 1%.

**My own contributions.**
- The shared formulation: mechanisms A/B/C and the cost targets.
- The window lemma (conditional).
- Lemma NC, the same identity as ptecons LEDGER, read as an NFS-like description.
- Prop NF.
- The joint statement with resbarrier.
- The literature map.
- Two pre-registered experiments:
  - E1, an exhaustive PTE fibre census, M ≤ 10, 62 213 970 multisets: P1–P3 and P5–P6 held, P4 failed as stated;
  - E2, a norm-ledger check: P7–P8 held; P9 failed because the cells lie below the lemma's regime.

---

## 1. Exact formulation (summary of `ptesplit_formulation.md`)

**Objects.**
- A **relation** is (u, v) with u ≡ v (mod N) and u ≠ v. Here u = ∏_j π_j(m) is split by construction over a factor base
  𝒫 with |𝒫| = L[α], and v is L[α]-smooth.
- A **split member** is P = ∏_{j≤M}(b_jX − a_j), with (a_j, b_j) ∈ [−B, B]².
- A **presentation** is (F̃, φ), where φ: Z[X]/(F̃) → A is a map to a Z/N-algebra. View A: φ(X) = m ∈ Z, with N | F̃(m).
  View B: φ(X) = μ ∈ Z/N.
- The **reduction** is r_P = P mod F̃, of degree < k and height ≤ h.

**Mechanisms.**
- **(A) exact:** P ≡ r_P mod F̃ in Z[X].
- **(B) seeded additive:** one seed relation, and fibre members P with deg(P − P₁) < k.
- **(C) non-exact:** a congruence only mod (N, F̃).

**Targets** (L-cost model, LA quadratic in |𝒫|). Write b = log B ≍ x^α(log x)^{1−α} and S = x^{2α}(log x)^{1−2α}.
- **(T-v)** log|v| ≲ S. For α = 1/3, 1/4, 1/10 this is x^{2/3}, x^{1/2}, x^{1/5} up to logarithms. Codex: a v that is
  common to all members needs no smoothness, since one formal column suffices.
- **(T-u)** M·log(max piece) ≥ x − 1.
- **(T-N)** About L[α] useful relations at cost L[α].

**Window lemma [the inequality is PROVED; the sizes are HEURISTIC: generic Δe and random residual smoothness].**
- Setting: mechanism (B) with pieces m − a_j and λ = log m. Then
  (P − P₁)(m) = Σ_{l<k} Δe_{M−l} m^l, with |Δe_{M−l}| ≤ 2C(M, l)B^{M−l}.
- Put q = M − k + 1 and t = k − 1. If the Δe have generic size, (T-v) forces q·b ≲ S and t·λ ≲ S, while (T-u) forces
  M·λ ≳ x. Hence

  λ ≳ x^{1−α}(log x)^α,  k − 1 ≲ (x/log x)^{3α−1},  M ≲ S/b + t ≍ x^α(log x)^{−α}.

- So under the heuristic, k = 1 for fixed α < 1/3, and M is bounded *above* by about x^α. Smaller M with larger λ is
  allowed.
- Caveat: the escape is exactly about exceptional, moment-conditioned fibres, for which "generic Δe" is not justified.
  Smoothness by identity could also bypass (T-v). This is a heuristic guide, not a barrier (codex #4).

## 2. Status by mechanism

### 2.1 Mechanism A (exact presentation)

**ptebarrier (final, 4ac68dc) [PROVED, with hypotheses].**
- Lemma R: let P, P′ be monic and split over Z with roots in [−B, B]. Suppose QP − Q′P′ has degree < K and height ≤ h,
  with deg Q, deg Q′ ≤ j. If M ≥ K + 2⌊log₂(2^{K−1}√K·h·B^{K−1})⌋ + 4j + 3, then QP = Q′P′.
- Theorem W: in that regime, every exact family h_A·P_A = F·Q_A + r_A gives one common relation. F may have any degree and
  height.
- Corollary P: an exact PTE family with d ≥ k has ≤ 1 member.
- Prop S: in View A, under (ρ+1)H(|m|+Y)^ρ(d+1)Y^d < N, an exact PTE family has ≤ 1 member once M − k ≥ d + 2ρ − 1.
- Scope: (k + j)·log B + log(H·H_Q) = o(x/log B) in View B. Outside it (O3) the question is open.

**Lemma NC (norm ledger) [PROVED identity; the same identity as ptecons LEDGER].**
- *Setting.* F̃ is monic. G is a monic irreducible factor of F̃ of degree δ ≥ 2 with N | G(m). If no factor G_i has
  N | G_i(m), then some gcd(N, G_i(m)) is a proper factor of N. Suppose P = ∏(b_jX − a_j) ≡ r (mod F̃) in Z[X]; r is
  integral because F̃ is monic. Let θ be a root of G.
- *Statement.* P(θ) = r(θ) in Z[θ], and

  ∏_j |Ĝ(a_j, b_j)| = |N(r(θ))| ≤ ‖r‖₁^δ‖F̃‖₂^{deg r}.

- *Proof of the bound.* Use |r(θ_i)| ≤ ‖r‖₁·max(1, |θ_i|)^{deg r} and ∏_i max(1, |θ_i|) = M(G) ≤ M(F̃) ≤ ‖F̃‖₂
  (Landau).
- *Consequences.*
  - At most log₂(‖r‖₁^δ‖F̃‖₂^{deg r}) pieces are non-units. This counts piece *types*; unit pieces may repeat without limit.
  - Unit types bθ − a: ≤ 2δ monic ones; O(δ) primitive ones for δ ≥ 3 (Thue; the Bombieri–Schmidt / Schmidt counts);
    O(log B) in the box for real quadratic orders and O(1) for imaginary quadratic orders.
  - δ = 1 is the rational case: the identity is an exact integer identity, and the only N-information is m ≡ c (mod N).
  - Non-monic F̃ needs denominator accounting and is not covered. For 2X² − 1, N(θ) = −1/2 (codex #7).
- *Withdrawn (codex #1).* The draft said that members with the same non-unit part are redundant. That is **false**.
  - Counterexample: F = X² − 2, m = 11, N = 119. Both (X + 1)² ≡ 2X + 3 and X² − 1 ≡ 1 have an empty non-unit part.
  - Yet the first gives 12² ≡ 5² (mod 119) and factors N.
  - Unit-dominated members are informative through the rational map. Redundancy must be judged in the full relation
    matrix.
- *Reading.* An exact member is an NFS-like relation: a sieving element r(θ) of degree < k with an element factorisation
  on the algebraic side.
  - Box-linear elements need not lie in a bounded-norm prime-ideal factor base, so NFS smoothness probabilities do not
    transfer automatically.
  - **NC alone excludes nothing for small M.** Its force is the regime M ≫ log₂‖r‖₁^δ, where ptebarrier's Lemma R already
    applies.

**Function-field contrast [assessment].** In BGJT the identity lives in the finite ring F_{q²}[X]/(I). Every piece X − a
is a unit there and the norm carries no size. Over Z the ledger charges every root with |F(a)| ≥ 2 (ptecons LEDGER).

### 2.2 Mechanism C (non-exact) and N-free families

**Non-exact members.**
- ptebarrier Prop NX [PROVED]: in the PTE setting with d ≥ k, distinct members have distinct reduced N-defects.
- This does **not** imply independent costs. Members can inherit compatibility from one difference identity.
- Codex gives an outside-regime example: F = X² − (N+2), P = X², P′ = X² − 1, with equal defects.
- For k = 1 the inherited compatibility is itself a split relation (ptebarrier NX1).
- Counting the attainable defect cells is **open** (O1).

**Codex blind [PROVED].** If d ≥ k, a shared remainder mod (N, F) forces P = P′ whenever some p | N exceeds 2B.

**Prop NF [PROVED; elementary].**
- *Statement.* Let S be a set of integers of bit-size ≤ T, fixed independently of N, and put T* = ⌈log₂(2^T + V)⌉. The
  balanced semiprimes N = pq ∈ [X, 2X] (p, q > √X/2) for which some s ∈ S and |v| ≤ V with s ≠ v satisfy N | s − v
  number at most

  |S|(2V+1)(T*/log₂(√X/2))².

- *Proof.* s − v ≠ 0 has bit-size ≤ T*, so it has at most T*/log₂(√X/2) prime factors above √X/2.
- *Consequence.* Balanced semiprimes in the interval number e^{x+o(x)}. So for |S|, V, T = e^{o(x)} the successful
  fraction is **at most** e^{−x+o(x)}.
- *Scope.* This shows that families succeeding for almost all N under these resource bounds must adapt to N. It does not
  establish ptecons' uniform-rate heuristic AVG. Codex's counterexamples to the universal form (squares; lcm(range) + 1)
  have exponentially many elements or exponential bit length.

### 2.3 Mechanism B (seeded additive PTE fibres)

- **Injection bound [PROVED; ptesplit, codex, ptebarrier].** |𝓕_K(h)| ≤ C(2B + k, k).
- **Kleiman 1975 [KNOWN].**
  - For an ideal pair, (M − 1)! divides the constant C = P − P′.
  - So at most 1 + 2V/(M − 1)! members of one seeded k = 1 cluster have |v| ≤ V (ptecons §3.4).
  - An L[α]-size cluster inside (T-v) therefore needs M log M ≲ S. Kleiman alone restricts the shape (few, huge pieces)
    but does not kill it.
- **Determinant method, k = 1 [PROVED from KNOWN; ptebarrier final].**
  - For M ≥ 4 a fibre has ≤ 2c·B^{1/2}(log B)^κ + 4 members, with c and κ absolute. The squarefree levels are bounded by
    the Binyamini–Cluckers–Novikov uniform bound (arXiv 2302.04209, Thm 2: c·d²·B^{1/d}(log B)^κ). Rolle's theorem
    allows at most two non-squarefree levels.
  - M = 2 is genuinely exceptional: X² − a² is one fibre of size B + 1 (codex #6).
- **Reading.** Mechanism B is a relation amplifier, with gain ≤ B^{1/2+o(1)} = L[α, O(1)] per seed.
  - Total cost ≥ (#seeds)·(cost per seed).
  - So any exponent statement needs an explicit **seed-cost hypothesis**, and we prove none.
  - Cheap special seeds exist. Codex: A = {0, 3, 4, 7}, A′ = {1, 1, 6, 6}, m = 24, N = P_A(24) − 1 = 349·491 gives
    residuals 1 and 37.
  - For general N with a fixed Q, no seed method below relation search is known [assessment].
- **Near-ideal fibres (k ≥ 2).** Controlled only under ptecons' NC heuristic. This is ptecons' open item 1.

**Vinogradov's mean value theorem [KNOWN; the analysis is agreed by all three agents].**
- The bound J_{s,t}(X) ≪_{s,t,ε} X^{s+ε} + X^{2s−t(t+1)/2+ε} controls Σ fibre², i.e. collisions.
- Its pointwise consequence is max fibre ≪ B^{M/2+ε}. That is weaker than the injection bound whenever k < M/2.
- Its constants are not uniform in (M, t), and the loss B^ε = L[α, ε] sits at the scale in question.
- **VMVT is the wrong tool.** Where a barrier exists, it is archimedean: remainder heights and norms.

### 2.4 Relation to resbarrier §2.6 (orchestrator's question)

**Representation [PROVED, scoped].**
- *Setting.* View B, exact members, implemented with two orders. O = Z[θ] carries the split side; O′ = Z[θ′], with
  φ′(θ′) = μ, factors the same residual r; the maps to Z/N agree on the residual generators.
- *Claim.* The candidates {(r(θ), r(θ′)) : deg r < k, ‖r‖ ≤ h} form a linked module of effective rank ≤ k, with
  bounded coordinates.
- *Classification.* For k ≥ 3 with O′ nonlinear, this is resbarrier §2.6. Rank 2 is resbarrier's C2 case.
- *Role of the split requirement.* It only strengthens the O-side condition, via ∏|Ĝ| = |N(r(θ))|.
- *Why two orders.* A single order is the one-ring situation of resbarrier Theorem E and Lemma ND, which extracts nothing
  *for the extraction schemes that ND covers*.

**Scope [codex #7].**
- Theorem E and Lemma ND do not prove that *every* useful extraction has this two-order, same-r form. The geometric
  representation also does not automatically carry resbarrier's sampling and smoothness assumptions.
- Not covered:
  - multivariate presentations without a finite étale specialisation;
  - individually recentred residues;
  - general piece dictionaries;
  - non-monic presentations;
  - pieces of higher degree. Norm multiplicativity survives there, but the linear Thue count does not.
- G_r is a generic-modulus heuristic, not an all-N closure.

**Consequence.**
- The exact two-order residue-view implementation of route 1's escape is an instance of resbarrier's single open case.
  The converse fails.
- The two open problems should be stated **jointly**. A rigorous rank-r ledger, or a proof of G_r, closes both for that
  implementation. Neither problem has a candidate construction.

## 3. Constructions (ptecons, final; summary)

No candidate meets (a)–(c). ptecons examined seven families:
- Dickson/Chebyshev;
- PTE, including Prouhet–Thue–Morse, Frolov and Choudhry chains, which fix O(log M) moments;
- cyclotomic and Galois orbits;
- elliptic preimages;
- (C₂)^r group rings;
- lattice/LLL;
- aligned clusters.

Proved there: LIN-R (real rigidity, via de Gua), CUN, and LEDGER. STAR combines proved counting with the heuristic NC, so
its consequence is heuristic.

Open there:
- systematic near-ideal cancellation;
- non-star complex-direction families;
- infinite-unit-group fields;
- an output-sensitive enumerator.

ptecons E2 (toy scale; confounded by different root bounds and pool truncation): the generic quartic's median residual
exponent was 0.565 against COS's 0.574. That is within the registered 0.02 margin, so it is not a gain.

Codex blind agrees:
- D_n(X) − c has no integer-split fibre for n ≥ 7;
- dilation multiplies the ideal constant by a^M.

## 4. Experiments vs frozen predictions (pre-registration 987a9c7)

### E1: exhaustive PTE fibre census

**Setup.** All multisets of size M over {0..R} (62 213 970 in total), grouped by (p₁..p_s), s = M − k. Grouping used a double
64-bit hash (see the exact verification below); the tests cross-check against brute force. Raw data:
`experiments/raw/ptesplit_E1.jsonl`.

**Post-hoc exact verification (disclosed; codex #10).** `src/ptesplit/verify.py` recomputes the exact power-sum
vectors of every member of every nontrivial hash group: 12660528 member-checks over all cells and k, **0 mismatches**
(`experiments/raw/ptesplit_E1_verify_posthoc.jsonl`). The census is therefore exact rather than probabilistic. Peak RSS
was 2.52 GB, again above the target.

**Resources.** 176 s; peak RSS **2.35 GB**, above the 2 GB target. This is a disclosed deviation; nothing else was running.

| (M, R) | T | k = 1 pairs (box) | max | k = 2 pairs (box) | max | k = 3 pairs (box) | max |
|---|---|---|---|---|---|---|---|
| (4, 60) | 635 376 | 4 878 (0.067) | 4 | 1 928 985 (5.8e4) | 29 | 1 548 937 735 (8.4e8) | 6 786 |
| (5, 60) | 8 259 888 | 623 (9e−8) | 2 | 396 757 (5.8) | 7 | 229 368 115 (6.3e6) | 180 |
| (6, 40) | 9 366 819 | 123 (5e−15) | 2 | 6 157 (3e−6) | 3 | 1 675 461 (49) | 14 |
| (7, 30) | 10 295 472 | 0 | 1 | 1 023 (2e−13) | 2 | 26 973 (4e−5) | 4 |
| (8, 24) | 10 518 300 | 0 | 1 | 2 (2e−21) | 2 | 3 924 (3e−12) | 3 |
| (9, 20) | 10 015 005 | 0 | 1 | 0 | 1 | 0 | 1 |
| (10, 18) | 13 123 110 | 0 | 1 | 0 | 1 | 0 | 1 |

The minimum k = 1 constants are 36, 5040 and 14 400 at M = 4, 5, 6: multiples of 3!, 4! and 5!. The (4, 60) row reproduces
codex's independent counts exactly (4878 / 1 928 985 / 1 548 937 735; max 4 / 29 / 6786).

**Predictions.**
- **P1 (injection bound): held.** The largest ratio is 6786/39711 at (4, 60, k = 3).
- **P2 (Kleiman): held**, 0 violations.
- **P3: held.** k = 1 pairs exist for M = 4–6 against box predictions ≤ 0.07, and none exist for M ≥ 7.
- **P4: FAILED as stated.**
  - Part 1, observed ≥ 10 × box wherever box ≥ 10⁻³, fails at (4, 60, k = 3): the ratio is 1.85. It holds in every other
    eligible cell (ratios 33 to ~10⁵).
  - Part 2, nontrivial k = 3 fibres for M ≤ 8: held.
  - *Post hoc:* with s = 1 prescribed moment the box model is already near-correct; structure dominates only once s ≥ 2.
- **P5 (k = 1 max ≤ 2√(R+1) + 2): held.** The maxima are 4, 2, 2, 1, 1, 1, 1.
- **P6 (max fibre ≤ 0.1 × injection for M ≥ 6, k ∈ {2, 3}): held.** The largest ratio is 2/325, at (8, 24, k = 2).

**Reading [measured].**
- Structure (symmetric and affine families) beats the box model by up to 10¹⁵ at small M. But fibres collapse to size ≤ 4
  by M = 7 and to size 1 by M = 9 at these entry bounds.
- No k ≤ 3 fibre of size > 4 exists for M ≥ 7. That is the opposite of the dense L[α] fibres (M → ∞) the escape needs.
- These are finite sizes. The asymptotic statements rest on Kleiman, the determinant method and Theorem W, not on E1.

### E2: norm compression on exact presentations

**Setup.** F ∈ {X³ − 2, X³ + X + 1, X² − 2, X² + 1, X⁴ + X + 1}, (M, B) ∈ {(4, 8), (5, 8), (6, 6)}, h ∈ {20, 200}; 30 cells,
299–5128 members each. Raw data: `experiments/raw/ptesplit_E2.jsonl`.

**Predictions.**
- **P7 (ledger identity and bound): held.** In all 30 cells ∏|F(a_j)| = |Res(F, r)| exactly, within the bound; 0 violations.
- **P8 (#(r, non-unit part) ≤ 3·#r): held.** This measures residual and non-unit-part multiplicity, not the rank of useful
  relations. The largest ratio is 1.94 (X² − 2, M = 6, h = 20). X² − 2 and X² + 1 show
  genuine multiplicity, from unit pieces a = ±1 and from factorisation non-uniqueness. The cubic and quartic F have
  ratio 1.00.
- **P9 (d ≥ 3, M ≥ 5: > ½ of members have ≤ 2 non-unit pieces): FAILED** in all 12 eligible cells (0.02–0.48).
  - *Post-hoc diagnosis:* the cells sit **below the lemma's regime**. The allowed non-unit count log₂(‖r‖₁^d‖F‖₂^{d−1}) is
    ≈ 20–30 for d = 3, h ≥ 20, while M ≤ 6. So the ledger permits every piece to be a non-unit, and the maximum observed
    non-unit count equals M in every cell.
  - The prediction was badly calibrated. The ledger has force only for M ≫ d·log₂ h, while the identity itself holds for
    every M. E2 is a correctness check of the identity and does **not** test unit domination.

## 5. Literature (task 5; cross-checked against the codex, ptecons and ptebarrier lists)

- **PTE.**
  - Prouhet (1851); Tarry/Escott (1910s); Wright, Quart. J. Math. 6 (1935), the conjecture that ideal solutions exist for
    every size.
  - Hua, *Additive theory of prime numbers*, on Tarry's problem and mean values.
  - Borwein–Ingalls, Enseign. Math. 40 (1994); Borwein, *Computational Excursions in Analysis and Number Theory* (2002),
    Ch. 11.
  - Borwein–Lisoněk–Percival, Math. Comp. 72 (2003).
  - Ideal solutions are known for sizes ≤ 12 except 11 (size 12: Kuosa–Meyrignac–Chen 1999).
  - Coppersmith–Mossinghoff–Scheinerman–VanderKam, arXiv 2304.11254: searches to size 16, including Z[i] and Z[ω].
  - Alpers–Tijdeman, J. Number Theory 123 (2007): the two-dimensional PTE problem.
  - Caley: thesis, Waterloo 2012, on the constants C_n; Math. Comp. 82 (2013), on PTE over Z[i].
  - Kleiman, J. reine angew. Math. 278/279 (1975): (n − 1)! | C_n. Rees–Smyth (1990); Filaseta–Markovich.
- **PTE and VMVT.**
  - Wooley, "Some remarks on Vinogradov's mean value theorem and Tarry's problem", Monatsh. Math. 122 (1996). This is the
    direct PTE/VMVT link: it bounds the degree attainable with s variables, not fibre sizes.
  - BDG, Ann. Math. 184 (2016); Wooley, Adv. Math. 294 (2016) (cubic) and Proc. LMS 118 (2019) (all degrees).
- **Determinant method.**
  - Bombieri–Pila, Duke 59 (1989); Pila (1995);
  - Heath-Brown, Ann. Math. 155 (2002); Salberger (2007, 2023);
  - Castryck–Cluckers–Dittmann–Nguyen, Algebra Number Theory 14 (2020);
  - Binyamini–Cluckers–Novikov (uniform affine curve bound; cited via ptebarrier, UNVERIFIED by us).
- **Thue counts.** Bombieri–Schmidt, Invent. Math. 88 (1987).
- **PTE in cryptography/factoring.**
  - Costello–Meyer–Naehrig, EUROCRYPT 2021 (ePrint 2020/1283): twin-smooth integers via PTE, N-free parameter generation.
  - Follow-ups optimise PTE-based twin-smooth searches (codex cites Ahrens; UNVERIFIED).
  - Moment subset sum / Reed–Solomon decoding (Gandikota–Ghazi–Grigorescu, ECCC 2016/176).
  - **No use of PTE in integer factoring or in breaking a modulus-dependent problem was found** (by us, codex, ptecons or
    ptebarrier).
- **Prior informal statement of the obstruction.** Granger–Kleinjung–Zumbrägel, AMC 12 (2018) §7: prime fields need
  "fundamentally new ideas".

## 6. What remains open (joint list)

1. **Non-exact and evaluated compressions** (ptebarrier O1). These are congruences mod N that are coherent across a
   family without being exact over Z. There is no barrier and no candidate. A sufficient primitive is smoothenum's
   output-sensitive split or doubly-smooth listing.
2. **resbarrier §2.6/§4.2:** rank-r ≥ 3 linked modules with two nonlinear sides. They are closed only under G_r. The
   exact two-order residue-view implementation of route 1 is an instance of this case (§2.4).
3. **Near-ideal fibres with systematic cancellation at an N-aligned point** (ptecons open 1; heuristic NC). Also any
   seed source for mechanism B below relation search.
4. **ptebarrier O2–O4.** These are:
   - pencils of ≥ 2 presentations (the missing lemma is three-term rigidity);
   - exact families outside the regime (large multiplier degree or height; tiny m with a high-degree small side);
   - other root rings and multivariate presentations;
   - non-monic presentations, which also lie outside Lemma NC.

## 7. Reviews

### 7.1 Codex blind round (`ptesplit_1.out`)

Same verdict: no construction; 1–3%. Codex contributed:
- the shared-remainder obstruction for d ≥ k;
- the injection bound (found independently);
- the 2B-smoothness of ideal constants;
- the pointwise VMVT bound;
- the EX qualification (exact identity plus N-dependent evaluation is informative);
- the PTM, doubling and Dickson audits;
- the M = 4 fibre counts, which our E1 reproduces exactly.

### 7.2 Codex adversarial round (`ptesplit_2.out`): MAJOR REVISION; all points applied

| # | Severity | Issue | Fix |
|---|---|---|---|
| 1 | MAJOR | "same non-unit part ⇒ redundant" is false: X² − 2, m = 11, N = 119 factors N from a unit-only member | withdrawn (§2.1) |
| 2 | MAJOR | "element-factor-base NFS" is a description, not a lower-bound reduction; KL hypotheses | relabelled; "NC alone excludes nothing for small M" |
| 3 | MAJOR | W/P/NX hypotheses dropped; NX does not imply independent costs (codex example) | regimes attached; NX scoped; ptebarrier final folded in |
| 4 | MAJOR | window lemma: logs; M is bounded above, not ≈ x^α; "generic Δe" unjustified on exceptional fibres | rewritten as a conditional heuristic (§1) |
| 5 | MAJOR | the seed argument proves nothing without a seed-cost bound; cheap special seeds exist | "bounded amplification; seed-cost hypothesis needed"; codex seed quoted |
| 6 | MINOR | determinant bound: M = 2 exception, repeated roots | M ≥ 4 and the +4 term from ptebarrier final |
| 7 | MAJOR | §2.4 is a scoped representation, not universal containment; non-monic denominators | scoped; uncovered cases listed |
| 8 | MINOR | Prop NF: bit-size of s − v; "at most"; lcm example | T* and an upper-bound reading |
| 9 | MAJOR | the verdict is stronger than its sources; STAR is heuristic; ptecons E2 margin | verdict rewritten; STAR relabelled; margin stated |
| 10 | MINOR | numbers (62 213 970; 2/325; 299–5128); P8 meaning; hash exactness | corrected; exact post-hoc verification, 0 mismatches |

Codex independently re-scored P1–P9 from the raw files and agrees: P4 and P9 failed, the rest held.

### 7.3 Opus critics

See §7.4 (addendum).

## 8. Reproduction

```bash
.venv/bin/python -m pytest -q tests/test_ptesplit.py
.venv/bin/python -m src.ptesplit.run_experiment E1   # about 3 min, 2.35 GB peak
.venv/bin/python -m src.ptesplit.run_experiment E2   # about 1 s
.venv/bin/python -m src.ptesplit.verify              # post hoc, about 3 min, 2.52 GB peak
```
