# ptesplit: synthesis on route 1's open escape (dense N-dependent split families / multigrade PTE)

Agent `ptesplit`, phase 16, 2026-10-20.

**Files.**
- Formulation (shared with ptebarrier and ptecons): `notes/ptesplit_formulation.md` (commit 641c3e2).
- Pre-registration: `notes/ptesplit_prereg.md` (commit 987a9c7, made before E1/E2).
- Code: `src/ptesplit/{census,nccheck,run_experiment}.py`. Tests: `tests/test_ptesplit.py`.
- Config: `experiments/configs/ptesplit.json`. Raw results: `experiments/raw/ptesplit_E{1,2}.jsonl`.
- Codex: blind round `codex/ptesplit_1.{md,out}`; adversarial round `codex/ptesplit_2.{md,out}` (§7).

**Sources synthesised.**
- `notes/ptebarrier_results.md` (task 2): **draft** ab0be49, whose codex/Opus reviews were pending when this note was
  written. Its results are used as stated there and are marked [pb-draft].
- `notes/ptecons_results.md` (task 3, final 6ce0300).
- `notes/resbarrier_results.md` (§2.6, §4.2), `notes/systeq_results.md`, `notes/smoothenum_results.md`,
  `notes/heurnfs_results.md`.

**Labels.**
- **PROVED**: a proof is given here or in the cited note.
- **KNOWN**: published, with a citation.
- **HEURISTIC**: the assumption is named.
- **measured**.
- **open**.
- **assessment**.

Novelty is never claimed.

---

## 0. Verdict

**Plain language.** The last open door of route 1 is closed in every form in which it differs from ordinary relation
search. That door was: many split ("smooth by construction") polynomials sharing one small reduction modulo an
N-dependent presentation.

The finite-field trick gets roughly q³ relations out of *one* congruence. Over the integers one congruence gives at most
*one* relation. The reason is archimedean.
- An exact identity in an integral ring has norms in Z. A small residual therefore forces almost every split piece to be
  a unit.
- Two distinct split polynomials with bounded roots cannot differ by a small, low-degree polynomial.

What survives is not a compression mechanism. It is either:
- the original relation-finding problem, i.e. each member pays its own congruence mod N (that is smoothenum's
  doubly-smooth listing); or
- in the "residue" view with no rational side, an *instance* of resbarrier's open rank-r case (§2.6/§4.2). There it is
  closed only under the heuristic G_r, and the PTE structure adds nothing to that case.

**The escape and resbarrier's open case: contained in it, not equal to it.**
- Route 1's escape, with its split-specific amplification removed, is a special case of resbarrier's rank-k,
  two-nonlinear-sides linked module. The split side is an extreme form of norm-smoothness.
- The converse fails. The rank-r case is broader, and nothing in it needs split structure.
- So the two open problems should be stated jointly as **one** problem: resbarrier §2.6. Closing §2.6 (proving G_r, or a
  rigorous rank-r ledger) would also close what is left of route 1.
- With the rational side present ("View A", integer evaluation), the escape is outside §2.6. It is closed rigorously
  in the size model: resbarrier KL, ptebarrier Prop S, and Lemma NC here.

**Verdict on route 1's open escape.**
- **Refuted** in its systematic (BGJT) form, unconditionally. This rests on ptebarrier Theorem W / Corollary P [pb-draft],
  ptecons LEDGER/CUN/Kleiman, and Lemma NC below.
- **Refuted** in the integer-evaluation view, rigorously in the size model (Prop S [pb-draft]; Lemma NC; resbarrier KL).
- **No exponent gain** in seeded additive form for ideal fibres. Kleiman spacing and the B^{1/2} determinant bound
  (ptecons §3.4; ptebarrier §4) make it a relation amplifier with gain ≤ B^{1/2+o(1)} per seed, so it cannot beat its
  own seed source.
- **Heuristically refuted** for near-ideal fibres (ptecons NC, the window lemma below).
- **Open** only as:
  - (i) resbarrier's rank-r case (G_r);
  - (ii) relation search itself (non-exact members; smoothenum's doubly-smooth listing);
  - (iii) ptebarrier's corners O2–O4, which have no candidate.

Assessment: < 1% that this escape yields sub-L[1/3]. This agrees with ptebarrier (< 1%), ptecons (< 1%), codex blind
(1–3%) and systeq (< 1%).

**My own contributions.**
- The formulation: mechanisms A/B/C, targets, and the window lemma.
- Lemma NC (norm compression): the same identity as ptecons LEDGER, read here as a reduction to element-factor-base NFS.
- Prop NF (N-free families are useless for almost all balanced semiprimes).
- The joint statement with resbarrier.
- The literature map.
- Pre-registered experiments:
  - E1, an exhaustive PTE fibre census, M ≤ 10, 7.2·10⁷ multisets: 5 of 6 predictions held, 1 failed as stated;
  - E2, a norm-compression check: 2 of 3 held; P9 failed because the cells sit below the lemma's regime.

---

## 1. Exact formulation (summary of `ptesplit_formulation.md`)

A **relation** is (u, v) with u ≡ v (mod N), u ≠ v. Here u = ∏_j π_j(m) is split by construction over a factor base 𝒫,
|𝒫| = L[α], and v is L[α]-smooth.

A **split member** is P = ∏_{j≤M}(b_jX − a_j) with (a_j, b_j) ∈ [−B, B]².

A **presentation** is (F̃, φ), with φ: Z[X]/(F̃) → A a Z/N-algebra map.
- View A: φ(X) = m ∈ Z with N | F̃(m).
- View B: φ(X) = μ ∈ Z/N, a residue only.

The *reduction* is r_P = P mod F̃, of degree < k and height ≤ h.

Three mechanisms by which N enters:
- **(A) exact presentation:** P ≡ r_P mod F̃ in Z[X];
- **(B) seeded additive:** one seed relation, and fibre members P with P − P₁ of degree < k;
- **(C) non-exact:** congruence only mod (N, F̃).

**Targets.** In the L-cost model (LA quadratic in |𝒫|):
- (T-v) log|v| ≤ c·x^{2α}(log x)^{1−2α}. For α = 1/3, 1/4, 1/10 this is log|v| ≲ x^{2/3}, x^{1/2}, x^{1/5}.
- (T-u) M·log(max piece) ≥ x − 1.
- (T-N) The family must supply ≈ L[α] useful relations at cost L[α].

Codex adds that a *common* residual v need not be smooth: one formal column suffices. A common residual is, however,
exactly what Theorem W shows yields only one relation.

**Window lemma [PROVED inequality + HEURISTIC sizes].** In mechanism (B), with pieces m − a_j and λ = log m:
- The residual difference is (P − P₁)(m) = Σ_{l<k} Δe_{M−l} m^l, with |Δe_{M−l}| ≤ 2C(M,l)B^{M−l}.
- For generic Δe, (T-v) and (T-u) give λ ≥ x^{1−α} and k − 1 ≤ x^{3α−1}.
- So for α < 1/3 only **ideal-PTE fibres (k = 1)**, with M ≈ x^α, can serve. At α = 1/3, k = O(1).
- Beating L[1/3], L[1/4] or L[1/10] thus requires all of the following:
  - k = 1 fibres of size L[α, c > 0] with M ≈ x^α → ∞;
  - ~|𝒫|/|fibre| seeds of that special shape at total cost L[α];
  - smooth v-sides.

The ideal-PTE requirement is exactly where Kleiman and the determinant method bite (§2.3).

## 2. Barrier status, by mechanism

### 2.1 Mechanism A (exact presentation): refuted

**Lemma NC (norm compression) [PROVED; the same identity as ptecons LEDGER].**
- Let F̃ be monic, and G | F̃ monic irreducible of degree δ with N | G(m). If no such G exists, some gcd(N, G_i(m)) splits N.
- If P = ∏(b_jX − a_j) ≡ r (mod F̃) in Z[X], then P(θ) = r(θ) in Z[θ] (θ a root of G), and

  ∏_j |Ĝ(a_j, b_j)| = |N(r(θ))| ≤ ‖r‖₁^δ‖F̃‖₂^{d−1}.

- *Proof of the bound:* |r(θ_i)| ≤ ‖r‖₁·max(1, |θ_i|)^{deg r}, and ∏_i max(1, |θ_i|) = M(G) ≤ M(F̃) ≤ ‖F̃‖₂ (Landau).
- Hence at most log₂(‖r‖₁^δ‖F̃‖₂^{d−1}) pieces are non-units. The rest are units bθ − a:
  - ≤ 2δ monic ones;
  - O(δ) primitive ones for δ ≥ 3 (Bombieri–Schmidt 1987; *count cited, not re-checked*);
  - O(log B) in the box for δ = 2.
- Members with the same non-unit part differ by exact unit relations, which give zero information (systeq EX).

**Reading.** Every exact member is an NFS relation for the pair (G, X − m), with sieving element r(θ) of degree ≤ k − 1 and
an *element* (not ideal) factorisation on the algebraic side.
- In View A this is a rank-k linked module **with** a rational side. It falls under resbarrier Lemma KL (rigorous, size
  model) and heurnfs B2 (t-dimensional sieving, α ≥ 1/3 under G).
- In the function field the identity lives in the finite ring F_{q²}[X]/(I). There every piece is a unit and the norm
  carries no size. That is the precise point where Z differs (also ptecons LEDGER).

**ptebarrier [pb-draft] closes the regime with M large, uniformly in the presentation.**
- Lemma R (rigidity): two split polynomials with roots in [−B, B] whose (multiplier-weighted) difference is small and of
  low degree are equal, once M ≥ K + 2⌊log₂(2^{K−1}√K·h·B^{K−1})⌋ + 4j + 3.
- Theorem W: hence a single presentation of any degree and height yields **one** relation.
- Prop S (View A): the split side uses ≤ deg F fixed pieces.
- With pieces ≤ L[β], β < 1/3, the regime condition (★) holds for every BGJT-shaped family and for every View-A
  parameter compatible with the L[β] target, except the corner O3.
- Lemma NC is the complement for small M: there the members are just NFS relations, so no dense family appears either.

### 2.2 Mechanism C (non-exact): this is relation search

- ptebarrier Prop NX [pb-draft]: in the PTE setting with d ≥ k, distinct members have distinct N-defects, so nothing is
  inherited.
- Codex blind [PROVED, codex]: with d ≥ k, a shared remainder mod (N, F) forces P = P′ whenever some p | N exceeds 2B.
- What is left is "find split u with u ≡ small (mod N)", the original problem. A sufficient primitive would be an
  output-sensitive enumerator of {P split} ∩ {P(m) ≡ small mod N}. That is smoothenum's algebraic or doubly-smooth
  listing, worth (4/9)^{1/3} or polynomial time, and not known.

**Prop NF [PROVED; elementary].** Let S be a set of integers of bit-size ≤ T, fixed independently of N.
- Statement: the balanced semiprimes N = pq ∈ [X, 2X] with some (s, v), |v| ≤ V, s ≠ v, N | s − v number at most
  |S|(2V+1)(2T/log₂(√X/2))².
- Proof: s − v has at most T/log₂(√X/2) prime factors above √X/2.
- For |S|, V, T = e^{o(x)} this is a fraction e^{−x+o(x)}.
- This makes rigorous, for semiprimes, what ptecons' AVG states conditionally. Codex's counterexamples to the
  universal form (squares; lcm(range) + 1) need either exponentially many elements or N-dependence. N must enter the
  construction.

### 2.3 Mechanism B (seeded additive PTE): refuted for ideal fibres, heuristic for near-ideal

- **Injection bound [PROVED; ptesplit/codex/ptebarrier].** |𝓕_K(h)| ≤ C(2B + k, k).
- **Kleiman 1975 [KNOWN].** For an ideal pair the constant C = P − P′ is divisible by (M − 1)!.
  - So the residuals v₁ + C_i of one seeded k = 1 cluster are pairwise spaced by multiples of (M − 1)!, and at most
    1 + 2V/(M − 1)! members have |v| ≤ V (ptecons §3.4).
  - A cluster of size L[α] inside the (T-v) window therefore needs M log M ≲ x^{2α}.
  - (T-u) then forces pieces of log-size λ ≳ x^{1−2α}·log x, consistent with the window lemma's λ ≥ x^{1−α}.
  - So Kleiman alone does not kill the cluster. It restricts the shape: few, huge pieces.
- **Determinant method, k = 1 [PROVED from KNOWN; ptebarrier §4, pb-draft].**
  - A k = 1 fibre has ≤ C·B^{1/2}(log B)^κ + 2 members, uniformly in M. The constants come from the uniform Bombieri–Pila
    bound of Binyamini–Cluckers–Novikov; the citation is flagged UNVERIFIED in ptebarrier.
  - So one fibre cannot fill a factor base of size B, and ≳ B^{1/2} independent N-dependent seeds are needed.
  - Each seed is itself a relation of the target shape. With a fixed Q, finding m with Q(m) ≡ small (mod N) is
    root-finding mod N.
  - Mechanism B is therefore a *relation amplifier* with gain ≤ B^{1/2+o(1)} = L[α, O(1)] per seed.
  - Total cost ≥ (#seeds)·(cost per seed) ≥ cost of one seed. An amplifier changes L-constants at most, never the
    exponent of the seed source [PROVED, trivial].
  - Seeds are relations of the target shape. With a fixed Q, finding m with Q(m) ≡ small (mod N) is root-finding mod N.
    Lattice (Coppersmith) methods solve Q(m₀ + y) ≡ v only when Y^M·V < N, where generically no solution exists
    (HEURISTIC count YV/N).
  - **Assessment:** no seed source cheaper than relation search is known.
- **Near-ideal (k ≥ 2).** The residual is controlled only under ptecons' no-cancellation heuristic NC. This is ptecons'
  open item 1. Our E1 shows how thin such fibres are at M ≥ 8 (§4).

**VMVT [KNOWN; analysis agreed by all three agents].**
- J_{s,t}(X) ≪_{s,t,ε} X^{s+ε} + X^{2s − t(t+1)/2 + ε} (BDG 2016; Wooley 2016, 2019) bounds Σ fibre², i.e. collisions.
- Its best pointwise consequence is J_{M/2,t} ≪ B^{M/2+ε}, weaker than the injection bound whenever k < M/2.
- Constants are non-uniform in (M, t), and the X^ε loss is B^ε = L[α, ε], the very scale in question.
- **VMVT does not give the barrier.** The barrier is archimedean (the height of the remainder; norms), not cardinality.

### 2.4 Coincidence with resbarrier §2.6 (requested by the orchestrator)

**Claim [PROVED reduction + assessment].** In View B (no integer lift of X), an exact member gives P(θ) = r(θ) in one
order O = Z[θ].
- By resbarrier Theorem E / Lemma ND (and ptebarrier §5.5, "single rings"), relations comparing the images of *one*
  element of *one* ring extract nothing.
- Information requires a second order O′ = Z[θ′] with φ′(θ′) = μ, in which the same r is factored. The candidate set
  {(r(θ), r(θ′)) : deg r < k, ‖r‖ ≤ h} is then a **rank-k linked module between two nonlinear orders**. That is resbarrier
  §2.6's open case (rank r = k ≥ 3, no rational side).
- The split requirement constrains only the O-side: ∏|Ĝ(a_j, b_j)| = |N(r(θ))|. It is a stronger smoothness condition
  (factorisation into box-linear elements), not a new source of candidates.
- By Theorem W [pb-draft], one presentation still yields ≤ 1 relation in the regime (★). So the module's size, not the PTE
  fibre, sets the count, and that count is resbarrier's G_r question.

**Consequence.**
- Route 1's open escape ⊆ resbarrier §2.6 (with View A excluded: there the rational side makes it KL-rigorous).
- §2.6 is strictly larger: the aggregate family of §4.2 uses {0, ±1}-combinations of r ≪ d linked elements, which are not
  split.
- **Joint statement:** under G_r both are closed with α ≥ 1/3. Without G_r, resbarrier gives a conditional formal floor
  L[1/4] for aggregate families. Split families are not aggregate families, but they inherit no better rigorous bound
  than §2.6 in View B.
- *Assessment:* View B with split O-side is, if anything, harder to populate than §4.2. Kleiman / LEDGER force any small
  r(θ) with split factorisation to be unit-dominated once M exceeds log₂|N(r(θ))|.

## 3. Constructions (ptecons, final; summarised)

No candidate meets (a)–(c). Seven families were examined:
- Dickson/Chebyshev: Niven; D_n(X) − c splits over Z only for n | 6 or degenerate c;
- PTE (Prouhet–Thue–Morse, Frolov, Choudhry chains, which fix O(log M) moments);
- cyclotomic/Galois orbits: the finite subgroups of PGL₂(Q) have order ≤ 12;
- elliptic preimages;
- (C₂)^r group rings;
- lattice/LLL;
- aligned clusters.

LIN-R (real rigidity, de Gua), STAR and CUN are proved there.

The codex blind round agrees: Dickson D_n(X) − c has no integer-split fibre for n ≥ 7, and dilations multiply the
ideal constant by a^M.

E2 of ptecons: N-aligned higher-degree split families did not beat the degree-2 COS shape at 64/80 bits.

## 4. Experiments vs frozen predictions (pre-registration 987a9c7)

### E1: exhaustive PTE fibre census

**Setup.** All multisets of size M over {0..R} (7.2·10⁷ in total), grouped by (p₁..p_s), s = M − k. Grouping used a double
64-bit hash (false-merge probability ≲ 10⁻²³ per cell, disclosed); the tests cross-check against brute force. Raw data:
`experiments/raw/ptesplit_E1.jsonl`.

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
- **P6 (max fibre ≤ 0.1 × injection for M ≥ 6, k ∈ {2, 3}): held.** The largest ratio is 14/12341.

**Reading [measured].**
- Structure (symmetric and affine families) beats the box model by up to 10¹⁵ at small M. But fibres collapse to size ≤ 4
  by M = 7 and to size 1 by M = 9 at these entry bounds.
- No k ≤ 3 fibre of size > 4 exists for M ≥ 7. That is the opposite of the dense L[α] fibres (M → ∞) the escape needs.
- These are finite sizes. The asymptotic statements rest on Kleiman, the determinant method and Theorem W, not on E1.

### E2: norm compression on exact presentations

**Setup.** F ∈ {X³ − 2, X³ + X + 1, X² − 2, X² + 1, X⁴ + X + 1}, (M, B) ∈ {(4, 8), (5, 8), (6, 6)}, h ∈ {20, 200}; 30 cells,
353–5128 members each. Raw data: `experiments/raw/ptesplit_E2.jsonl`.

**Predictions.**
- **P7 (ledger identity and bound): held.** In all 30 cells ∏|F(a_j)| = |Res(F, r)| exactly, within the bound; 0 violations.
- **P8 (#(r, non-unit part) ≤ 3·#r): held.** The largest ratio is 1.94 (X² − 2, M = 6, h = 20). X² − 2 and X² + 1 show
  genuine multiplicity, from unit pieces a = ±1 and from factorisation non-uniqueness. The cubic and quartic F have
  ratio 1.00.
- **P9 (d ≥ 3, M ≥ 5: > ½ of members have ≤ 2 non-unit pieces): FAILED** in all 12 eligible cells (0.02–0.48).
  - *Post-hoc diagnosis:* the cells sit **below the lemma's regime**. The allowed non-unit count log₂(‖r‖₁^d‖F‖₂^{d−1}) is
    ≈ 20–30 for d = 3, h ≥ 20, while M ≤ 6. So the ledger permits every piece to be a non-unit, and the maximum observed
    non-unit count equals M in every cell.
  - The prediction was badly calibrated. The lemma's force needs M ≫ d·log₂ h. E2 is a correctness check of the identity
    only.

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

1. **resbarrier §2.6/§4.2:** rank-r ≥ 3 linked modules with two nonlinear sides. Closed only under G_r (conditional floor
   L[1/4] for aggregate families). Route 1's View-B residue is an instance of it.
2. **Relation search itself.** Non-exact members; an output-sensitive split enumerator in {P(m) ≡ small mod N}
   (smoothenum's listing primitive; ptecons open 3).
3. **Near-ideal fibres with systematic cancellation at an N-aligned point** (ptecons open 1, heuristic NC).
4. **ptebarrier O2–O4:** pencils of ≥ 2 presentations, which need three-term rigidity; the corner |m| tiny with k ≍ x^{2β};
   other root rings and multivariate presentations.

No candidate is known for any of them.

## 7. Reviews

- **Codex blind (`ptesplit_1.out`).**
  - Same verdict: no construction, 1–3%.
  - It contributed:
    - the shared-remainder obstruction for d ≥ k (§2.2);
    - the injection bound (independently);
    - the 2B-smoothness of ideal constants;
    - the pointwise VMVT bound;
    - the EX qualification: an exact identity plus an N-dependent evaluation is informative, which is our View A and
      Lemma NC;
    - the PTM/doubling and Dickson audits;
    - the M = 4 fibre counts, reproduced exactly by our E1.
- Codex adversarial (`ptesplit_2`), Opus proof critic, Opus literature critic: see §7.1 (addendum).

## 8. Reproduction

```bash
.venv/bin/python -m pytest -q tests/test_ptesplit.py
.venv/bin/python -m src.ptesplit.run_experiment E1   # ~3 min, 2.35 GB peak
.venv/bin/python -m src.ptesplit.run_experiment E2   # ~1 s
```
