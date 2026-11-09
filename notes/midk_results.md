# midk: the uncovered intermediate-k regime — results

Agent `midk`, phase 31 follow-up, 2026-11-09. **v2** — the Opus proof critic (ACCEPT WITH FIXES on
the ledger; ACCEPT WITH FIXES on SPL-INT's core, REJECT of its over-broad consequence) and codex
blind `midk_1` are applied; codex adversarial `midk_2` recorded in §9. **The central honesty fix
(critic): the L[1/4] figure is the typical-size *model* floor, not a rigorous floor — only the AM–GM
step is proved, and resbarrier §2.6 is explicit that "without these heuristics we have no positive
floor for this sub-case."** Lemma SPL-INT is narrowed to its correct form and its consequence demoted
to an assessment with named gaps.

Files:
- this note; pre-registration `notes/midk_prereg.md` (commit `df4a933`, before experiments).
- Code: `src/midk/ledger.py` (E1 exponents; reuses `src/resbarrier/ledger.py`),
  `src/midk/independence.py` (E2; Lemma SPL-INT), `src/midk/experiments.py`.
- Tests: `tests/test_midk.py` (10 pass).
- Data: `experiments/raw/midk_e1.jsonl`, `experiments/raw/midk_e2.jsonl`.
- Codex GPT-6.1-sol: `codex/midk_1.{md,out}` (blind), `codex/midk_2.{md,out}` (adversarial).

Labels: **proved** (elementary, here, or internal/unrefereed chain + cited published results),
**proved (model)** (exact consequence of a stated size model), **heuristic**, **open**, **measured**
(exact toy), **UNVERIFIED**. No novelty claimed beyond "no prior formulation found by us".

Notation: x = log N; L[α, c] = exp((c+o(1))(log N)^α(log log N)^{1−α}). "α" is the L-exponent.

---

## 0. Verdict

**Plain language. It is not possible to beat L[1/3] by any intermediate-k regime we can analyse —
with one genuine but unexploited hole, which we both bound and argue is factoring-equivalent to
enter.**

First job — *which* regime. Scanning the catalogue, **three** families (an assessment of completeness,
not a theorem — a formal algorithm class + exhaustive reduction would be needed, both critics) carry a
discrete parameter k whose intermediate values could in principle do better than the endpoints, and
whose block at intermediate k is a counting/independence hypothesis rather than a proved impossibility.
Only family A's block is purely counting; C's is deterministic (Coppersmith), listed for completeness:

| family | parameter k | counting hypothesis | endpoint(s) | intermediate k if hypothesis fails favourably |
|---|---|---|---|---|
| **A** rank-r linked modules (`resbarrier` §2.6/§4.2, `rankdegen` §5–6) | rank r, 3 ≤ r ≤ d | **G_r** (Schmidt point counts + independence of the r−1 linking congruences mod N); `rankdegen` LP_exc | r=2 → 1/3 (proved, Cor RB2) | **α → 1/4** in the typical-size model (no rigorous positive floor < 1/3) |
| **B** descent ladder (`multidescent`/`pc3`/`c12rat`/`c12alg`) | descent level j, 2 ≤ j < ∞ | **C₁/C₂** pair-correlation | j=1 → 2/5 (proved L[2/5]) | α_j = (j+1)/(3j+2) ∈ (1/3, 2/5); **never ≤ 1/3** |
| **C** floor-sum ladder (`floorhard`/`floorsum`/`topblock`) | degree j / box side | Coppersmith known-bits barrier | j=2 → N^{1/4} | N^{1/(j+2)}, a **deterministic** metric ≫ L[1/3] |

**Only Family A has a favourable-failure exponent below 1/3.** Family B's ladder is strictly above
1/3 at every *finite* intermediate j (its infimum 1/3 is the j→∞ endpoint, needing nested
descendability for which there is no tool — `recurse` §3.4). Family C lives in the deterministic
N^{1/(j+2)} metric, where every value is ≫ L[1/3]; it is a route to beat the deterministic N^{1/5},
not L[1/3]. So the **weakest-looking intermediate-k regime that could beat L[1/3] is the rank-r
aggregate linked-module family** (the "no-rational-side rank ≥ 3" case of `resbarrier`, the
frame-less case of `rankdegen`), whose typical-size-model balance is L[1/4].

Overriding question — *can any intermediate k beat L[1/3]?* For Family A:
1. **The exponent ledger.** The per-candidate log-norm is log X = κ·x/d + d·c, κ = height cost of the
   linking, c = log(E√r) > 0. AM–GM gives **log X ≥ 2√(κ c x)** (proved, exact, E1). The subtlety the
   critic insists on: the only *rigorous* term is the height term κx/d (2-plane resultant, Bean/Thunder
   area on a rank-2 subplane, κ ≥ 1 up to an O(d log r) loss) and it → 0 as d grows; the x^{1/2}
   exponent — hence the 1/4 — comes from the **degree term d·c, the typical-size/balanced-profile
   heuristic**. So:
   - **Fallback (G_r fails, κ=1): α = 1/4 in the typical-size model** (attained at d ≈ √(x/log x),
     r ≈ x^{1/4}(log x)^{3/4}, u ≈ x^{1/4}); within that model no (d, r) gives α < 1/4. **This is not a
     rigorous floor** — `resbarrier` §2.6 relabelled exactly this (critic A, M3): "without these
     heuristics we have no positive floor for this sub-case"; `rankdegen` §0: "no ledger" for the
     frame-less configurations. Codex `midk_1` sharpens it: a congruence penalty N^{−r^γ} (partial
     collapse, 0 ≤ γ ≤ 1) gives optimum **α = 1/(4−γ)**, interpolating 1/4 (γ=0) ↔ 1/3 (γ=1).
   - **Under G_r (κ = r−1): α = 1/3**, r = 2 optimal — from the KL ledger 2√((r−1)xC/r), increasing in
     r (gr_min/kl_min height floors also non-decreasing, E1). No intermediate r beats 1/3 under G_r.
2. **Attacking G_r both ways.**
   - *Construct a favourable failure (dependent linking at height N^{1/d}).* Dependence sources:
     (I) a global char-0 relation — ring-power/subfield, pays geometric height (`rankdegen` Lemma SUB,
     E2), or the Galois/automorphism-graph case (`rankdegen` type A, Lemma CC: a finite useless split
     set); (II) a dependence through a Z/N intertwiner whose rank differs across the prime factors of
     N — **Lemma SPL-INT (proved, narrowed): a rank-discrepant intertwiner splits N**. The explicit
     witness (codex `midk_1`): the cyclotomic automorphism-graph family links all r generators at
     **height 1** (the independence *count* is genuinely false, E2), but the two sides are
     Galois-conjugate, so norms and modular square roots coincide — no gcd: a counting-hypothesis
     failure that is **useless**. *Scope (critic):* SPL-INT does **not** show every exploitable failure
     splits N; it misses (1) a prime-uniform mod-N dependence (how base-m NFS links mod N without
     factoring), (2) a global *non*-ring-power frame-less dependence (Weil restrictions over a
     non-imaginary-quadratic subfield; `rankdegen`'s finite-volume frame-less example) — the actual
     open hole, and (3) nonlinear / ring-map dependence. So "exploitable failure ⇒ factoring" is an
     **assessment**, proved only for the rank-discrepant-matrix case.
   - *Prove G_r.* Not done (Schmidt-at-growing-(d,r) + genuine mod-N independence remain assumptions).
     The independence baseline holds at small scale (E2 P2a, |z| < 1.5).
3. **Second-order obstruction.** The minor-product invariant Δ = ∏_{|I|=r} M_I² ∈ Z∖{0} with
   N^{2·C(D−2,r−2)} | Δ (unramified linking) + Hadamard gives a **deterministic** height penalty
   h₂ ≥ (r−1)x/(D−1) (codex `midk_1`; = `rankdegen` Corollary FN / frame–Hölder), which **re-imposes
   α ≥ 1/3 across the whole NONSINGULAR range** (all maximal minors nonzero), with no independence
   heuristic. The hole is exactly the **singular** case (some M_I = 0 ⇔ frame-less / rank-deficient),
   where Δ = 0 and the penalty vanishes; LP_exc would close the large-box part but is unproved, and no
   area theorem is known for the frame-less corner. The decoupled fibre-product sub-corner is
   L[1/2]-type (worse than useless).

**Bottom line.** No intermediate-k regime we can analyse beats L[1/3]. The one place the barrier is
only a counting hypothesis (Family A, rank r ≥ 3, no rational side, **singular/frame-less**) has a
typical-size-model balance of L[1/4] — **not a rigorous floor; there is no proven positive floor
below 1/3 there** — and a genuine, uninhabited hole. Every favourable failure we could construct is
either nonsingular (then α ≥ 1/3 by the minor-product penalty) or useless (Galois-conjugate /
rank-discrepant ⇒ splits N). This is an assessment plus one narrow elementary lemma, not a lower
bound for factoring.

---

## 1. The three families and why only A matters for L[1/3]

The target is barrier arguments blocked at an *intermediate* k by a counting/independence hypothesis
(not a proved area/discriminant/resultant/Coppersmith wall). The catalogue has exactly three.

### 1.1 Family A — rank-r linked modules (the G_r regime)

`resbarrier` §2.6, §4.2 and `rankdegen` §5–6. A linked rank-r module is a family
(w_{1j}, w_{2j}) ∈ O₁ × O₂, j = 1..r, with φ₁(w_{1j}) = φ₂(w_{2j}) mod N; candidates are
Σ_j a_j(w_{1j}, w_{2j}). The two norm forms W₁, W₂ have degree d each. Blocks:
- **r = 2 (endpoint):** rigorous. Cor RB2 (`resbarrier`): N² | disc, Bean 1994 area theorem +
  Thunder/Davenport lattice counting ⇒ #{box points with small |F|} ≤ 16 N^{−k/(D(D−1))}X^{2/D}+…,
  giving α ≥ 1/3. Proved impossibility.
- **r ≥ 3 with a rational side:** rigorous (Lemma KL, Minkowski on the rank-(r−1) kernel lattice):
  min(h₁+h₂) = min(1, (r−1)/d₂), non-decreasing in r, α ≥ 1/3.
- **r ≥ 3, two nonlinear sides, no rational side (intermediate k):** **only G_r closes it.** This is
  the uncovered regime. G_r = Schmidt 1993 point counts + independence/equidistribution of the r−1
  linking congruences mod N. Under G_r, aggregate families are absent for typical N (expected count
  N^e, e ≤ 1 − (r−1)(1−1/d) < 0 for r ≥ 3). If G_r fails favourably, the formal balance is L[1/4].

### 1.2 Family B — the descent ladder

`multidescent` Theorem 2L, `pc3`/`c12rat`/`c12alg`. Parameter j = descent level. Cost in the size
model α_j = (j+1)/(3j+2), blocked by the pair-correlation hypotheses C₁/C₂ (`pc3` Prop DC). The
endpoints are j=1 (α=2/5, **proved** L[2/5], Theorem U) and j→∞ (α→1/3). **Key fact (E1):** α_j is
strictly decreasing with inf 1/3, and **α_j > 1/3 for every finite j**. So even a fully favourable
resolution of C₁/C₂ at intermediate j gives α_j ∈ (1/3, 2/5) — it *approaches* but never reaches or
beats 1/3. Reaching 1/3 needs j = ω(1) levels with geometrically contracting constants and a nested
descendability lemma for which **no tool exists** (`recurse` §3.4, blocked by endogenous children).
So Family B is not a route below 1/3; its conditional prize is L[3/8], L[4/11], … , all above 1/3.

### 1.3 Family C — the floor-sum ladder

`floorhard`/`floorsum`/`topblock`. Parameter j = Taylor degree / box side, intermediate regime
√b ≪ H ≪ b^{1−o(1)}, j ≥ 3. This is a **deterministic** algorithm with cost **N^{1/(j+2)}**
(partial-key-exposure: "factor p from its top (1/(j+2))·log₂N bits"). Endpoint j=2 is Coppersmith
N^{1/4}; j ≥ 3 would give N^{1/(j+2)} < N^{1/5}, blocked by the Coppersmith known-bits barrier. The
metric is N^α, and N^{1/(j+2)} ≫ L[1/3] for every fixed j. So Family C is a possible route to beat
the *deterministic* record N^{1/5}, not L[1/3]. (It is included because the task names "box-side /
degree j ≥ 3" explicitly; it is recorded, then set aside for the L[1/3] question.)

**Conclusion of the mapping.** Only Family A has a favourable-failure exponent below 1/3. It is the
regime the reviewer flagged. The rest of the note is about Family A.

## 2. The cost ledger of Family A, exactly (task item 1)

Work per side with degree d, total degree D = 2d, rank r, box side E (a_j ∈ {0, ±1} unless noted),
c := log(E√r) the per-degree box factor (c ≥ 1 whenever r ≥ 3 or E ≥ 1). The candidate's per-side
log-norm is
    log X = (family log-height) + (degree term) = κ·x/d + d·c,
with κ the *height cost of the linking*:
- **κ = 1 (G_r fails favourably):** the r−1 linking conditions are dependent, so the family height is
  pinned only by the rigorous 2-plane resultant, x/d, independent of r.
- **κ = r − 1 (G_r holds):** each of the r−1 independent congruences mod N costs one height unit, so
  the family height is (r−1)x/d (= `resbarrier` `gr_min`'s leading term for d ≥ r).

**Lemma LEDGER (AM–GM step exact/proved, E1; the a=1/2 exponent is the typical-size model).** min over
d>0 of κ·x/d + d·c is **2√(κ c x)** at d = √(κx/c). Hence, *within the typical-size model* (balanced
embedding profiles, the candidate taken at coefficient box E):
- **Fallback floor:** log X ≥ 2√(c x) = x^{1/2}·(2√c). With the aggregate box taken from the r sign
  patterns (a_j ∈ {0,±1}, E = O(1)), c = log(E√r) grows only logarithmically in the candidate count,
  so the box entropy enters log X with exponent **b = 0**, giving α = alpha_floor(1/2, 0) = **1/4**
  (`src/resbarrier/ledger.py::alpha_floor`). (For large E, c ≈ C/r and one recovers b = 1/2, the rank-2
  behaviour; the aggregate route is precisely the small-E/large-r corner.) **Within this model no
  (d, r) gives a < 1/2, so no α < 1/4.** The L[1/4] optimum is attained at d ≈ √(x/log x),
  r ≈ x^{1/4}(log x)^{3/4}, log X ≈ √(x log x), B = L[1/4], u ≈ x^{1/4} (E1 `attaining_scaling`;
  matches `resbarrier` §4.2 / codex `resbarrier_2` #7). **Not rigorous:** the x^{1/2} exponent rests on
  the degree term d·c (typical-size heuristic); the only rigorous term, κx/d (2-plane resultant, up to
  an O(d log r) loss), vanishes as d grows. resbarrier §2.6: "no positive floor for this sub-case"
  without the heuristic.
- **G_r floor:** log X ≥ 2√((r−1)c x); now r−1 ≍ the box entropy (candidates from the r sign
  patterns), so **b = 1/2**, α = alpha_floor(1/2, 1/2) = **1/3**, and log X grows like √(r−1), so the
  optimum is r = 2 (the standard rank-2 ledger). G_r and KL floors are non-decreasing in r
  (E1 `floor_monotone`): no intermediate r helps.

So the exact answer to "does any intermediate r give α < 1/3": **only if G_r fails, and then the
optimum in the typical-size model is α = 1/4** (α = 1/(4−γ) for a partial penalty N^{−r^γ}),
attained at r ≈ x^{1/4}(log x)^{3/4}. **Caveat (critic):** only the AM–GM step is proved; the x^{1/2}
exponent rests on the heuristic degree term d·c, so **this is the model floor, not a rigorous floor**
— `resbarrier` §2.6 ("no positive floor for this sub-case" without the typical-size heuristic).

## 3. Attacking G_r (task item 2)

### 3.1 Construct a favourable failure

To get κ = 1 at r ≥ 3 one needs the r−1 linking congruences dependent *without* paying the
geometric height that ring-power constructions pay. The possible sources of dependence split cleanly:

- **(I) Global (characteristic-0) dependence.** A Z-linear relation among the generators that holds
  before reduction. Then the module is a subfield/ring-power module: `rankdegen` Lemma SUB reduces it
  to a smaller non-degenerate module of the same candidates, and the top generator w = α^{r−1} has
  height H₀^{r−1} — **geometric in r** (E2 `P2b_ring_power_height`). This does not reach height
  N^{1/d}; it lands back on the r = 2 / KL floor, α ≥ 1/3.
- **(II) Dependence realised only mod N.** The two sides are intertwined by a Z/N-module map T
  (second-side reductions = T·first-side reductions) whose rank differs across the prime factors of N.

**Lemma SPL-INT (proved, elementary; narrowed after both critics).** Let N be composite, T ∈
Mat_m(Z/N), λ ∈ Z/N. Fraction-free Gaussian elimination of (T − λI) over Z/N has exactly three
outcomes: (0) **full rank** — λ is not a Z/N eigenvalue, no kernel, no factor; (i) a nonzero **kernel
vector** v with Tv ≡ λv (a *global* Z/N eigenvector), no factor; (ii) a nonzero **non-unit pivot** g,
and gcd(g, N) is a nontrivial factor of N. **The sufficient condition for outcome (ii) is
rank_p(T − λI) ≠ rank_q(T − λI) for some p, q | N** (unit pivots alone give the same rank mod every
prime, so a non-unit pivot must appear). Being non-scalar is **not** sufficient (outcome (i): T =
[[λ,1],[0,μ]]); nor is "λ an eigenvalue mod N" (outcome (0): N=21, T = [[0,−1],[1,0]], det(T−λI) =
λ²+1 a unit for all λ — codex `midk_2`).
*Proof.* A nonzero pivot in Z/N is a unit or has 1 < gcd(·, N) < N (a zero divisor ⇒ factor). If all
pivots are units the elimination equals the field computation mod every prime, giving equal rank mod
each p | N; so rank_p ≠ rank_q forces a non-unit pivot — outcome (ii). ∎

*Consequence (assessment, NOT a theorem).* A favourable G_r failure realised by a **rank-discrepant**
Z/N intertwiner is factoring-equivalent to exhibit. The explicit witness (codex `midk_1`): the
cyclotomic automorphism-graph family (K = Q(ζ_m), φ₂ = φ∘σ⁻¹, w_{2j} = σ(w_{1j})) links all r
generators at **height 1** — so the counting half of G_r is genuinely false — but the two sides are
Galois-conjugate, so norms and modular square roots coincide and no gcd results (`rankdegen` type A,
Lemma CC). **Gaps the lemma does not cover (critic):** (1) a prime-uniform mod-N dependence (base-m
NFS links mod N without factoring — case (i), not a split, and no height saving shown); (2) a global
*non*-ring-power dependence — the genuinely open frame-less configurations (`rankdegen`: Weil
restrictions over a non-imaginary-quadratic subfield, the finite-volume frame-less example); (3)
nonlinear / ring-map dependence. So the dichotomy is NOT exhaustive; it covers the matrix-intertwiner
case only, and formalises `resbarrier` §4.2's "cyclotomic needs roots of unity mod p" / "N-is-a-norm
is useless" for that case. E2: rank-discrepant intertwiners split N 18/18 (recovering p or q); global
scalars and the [[λ,1],[0,μ]] non-scalar do not; the automorphism-graph witness links at height 1 but
is diagonal/useless.

### 3.2 Prove G_r

Not achieved. G_r has two halves: (a) a Schmidt-type upper bound on the number of rank-r norm forms
of bounded height (Schmidt 1993, Monatsh. Math. 115; a fixed-dimension estimate applied at growing
d, r — an assumption), and (b) independence/equidistribution of the r−1 linking congruences mod N.
The independence half holds at small scale (E2 P2a: the rate that r uniform projective pairs are all
proportional mod p tracks p^{−(r−1)} within a 3σ Poisson band for r = 2, 3, 4 and p up to 10⁴).
Lemma SPL-INT shows the *failure* of (b) is factoring-hard to use, but it is not a proof of (b) as a
counting statement. So G_r stays a heuristic; what we add is that its favourable failure is
unexploitable without splitting N.

## 4. Second-order obstruction (task item 3)

Does an area theorem / union bound / "no intermediate prize" re-impose α ≥ 1/3 across the whole
intermediate range? Synthesis of `rankdegen`:
- **Frame-coverable configurations: closed.** `rankdegen` Theorem FH (frame–Hölder, a side-aware
  weighted form of Bean 1994's Hölder step, Bean–Thunder's area theorem per frame) gives Φ ≥ (r−1)/r,
  hence α ≥ 1/3, for general-position modules and for degenerate types C1/C2 with small loss ε′ — the
  whole part of the intermediate range that admits r+1 rows in general position. This *is* a
  second-order area obstruction and it covers most of the range.
- **Large boxes: closed under LP_exc only.** Points in excess of the volume bound occur at the L[1/3]
  scale (`rankdegen` §6.1), but wherever analysed they lie on lower-rank sublattices charged by a
  smaller ledger. LP_exc would make this a theorem; it is unproved.
- **Frame-less configurations: genuinely open.** Relative modules over a common subfield k ≠ Q not
  imaginary quadratic (infinite volume, no frames); `rankdegen`'s finite-volume frame-less example;
  decoupled fibre products. These are exactly where a κ = 1 aggregate family would have to live (its
  r−1 dependent conditions make it rank-deficient = frame-less). The decoupled fibre-product corner
  gives an **L[1/2]-type** cost (worse than useless). No frame-less family that beats rank 2 is known.
- **"No intermediate prize" analogue.** `cov38b` Prop 0 is the exact parallel in the deterministic
  world: inside its design family any coverage exponent c ≤ 3/8 gives exactly N^{1/5}; there is no
  intermediate c. Family A has the same flavour under G_r (the floor is monotone in r, so the
  endpoint r = 2 is optimal), but *without* G_r the monotonicity is exactly what fails, and no area
  theorem is known for the frame-less corner. So, unlike `cov38b`, Family A's intermediate range is
  **not** provably closed.

**So a genuine hole remains**, confined to frame-less (singular, Δ = 0) rank-r modules — floored at
L[1/4] only in the typical-size model (no proven positive floor below 1/3), experimentally
uninhabited, and partly (rank-discrepant case only, Lemma SPL-INT) factoring-equivalent to enter.

## 5. Experiments (pre-registered, `notes/midk_prereg.md`, commit `df4a933`)

**E1 (exact exponents, `experiments/raw/midk_e1.jsonl`).** All held.
- P1a: fallback α = alpha_floor(1/2,0) = 1/4; AM–GM floor log X = 2√(cx), a = 1/2, no α < 1/4. Held.
- P1b: under G_r α = 1/3; gr_min and kl_min non-decreasing in r at d ∈ {6,10,20}. Held.
- P1c: ladder α_j = (j+1)/(3j+2) strictly decreasing, all > 1/3, inf 1/3, j=1 → 2/5. Held.
- P1d: floor-sum N^{1/(j+2)} recorded as deterministic metric ≫ L[1/3]. Held.

**E2 (exact independence/construction, `experiments/raw/midk_e2.jsonl`).** All held.
- P2a independence baseline: empirical all-proportional rate vs p^{−(r−1)} prediction within |z| < 1.5
  for (p, r) ∈ {101,1009,10007}×{2,3,4}, 2·10⁵ trials each. Held.
- P2b(II)/Lemma SPL-INT: 9/9 **rank-discrepant** intertwiners (m = 2,3,4; three N = pq up to 41 bits;
  genuine non-unit pivot via the unit-preferring elimination, not an entrywise-zero matrix — the
  earlier tautological construction and a spurious-entry-gcd bug were fixed per the critics) split N,
  recovering p or q; the 9 global scalars, the non-scalar [[λ,1],[0,μ]], and the outcome-0 example
  (N=21, [[0,−1],[1,0]], all λ) do not split. Held. The cyclotomic automorphism-graph witness links
  all r generators at height 1 (counting failure) but is diagonal/useless (m=5,7,11). Held.
- P2b(I): ring-power family height is geometric in r (top power H₀^{r−1}); "N is a norm" is full-rank
  (no intermediate r), useless. Held.

(P2c, the excess-on-sub-rank classification, was folded into the §4 synthesis of `rankdegen` §6.1–6.2
rather than re-run; `rankdegen` E3 already measured it. Disclosed deviation.)

## 6. Answers to the task

1. **Weakest regime:** Family A, rank r ≥ 3, no rational side, singular/frame-less — blocked only by
   G_r; its favourable-failure exponent is L[1/4] (typical-size model), the only intermediate-k
   exponent below 1/3 in the catalogue.
2. **Ledger & minimiser:** log X = κx/d + dc; AM–GM (proved) ⇒ log X ≥ 2√(κcx); in the typical-size
   model α = 1/4 (fallback) at r ≈ x^{1/4}(log x)^{3/4}, d ≈ √(x/log x) (α = 1/(4−γ) for partial
   penalty N^{−r^γ}); under G_r α = 1/3 at r = 2. No intermediate r beats 1/3 unless G_r fails; the
   1/4 is a model floor, not rigorous (resbarrier §2.6: no positive floor without the typical-size
   heuristic).
3. **Both ways on G_r:** construction — a height-free dependence is either a rank-discrepant mod-N
   intertwiner, which
   splits N (Lemma SPL-INT, proved; E2 18/18); proof — G_r not proved (Schmidt-at-growing-d +
   independence remain assumptions), independence holds at small scale.
4. **Second order:** the minor-product height penalty (= frame–Hölder) closes the whole NONSINGULAR
   range (α ≥ 1/3, no heuristic); LP_exc would close large boxes; the singular/frame-less corner is a
   genuine hole, floored at L[1/4] only in the typical-size model. No construction inhabits it.

## 7. Reproduction

```bash
PYTHONPATH=. .venv/bin/python -m pytest -q tests/test_midk.py          # 12 tests
PYTHONPATH=. .venv/bin/python -m src.midk.experiments all              # writes midk_e1/e2.jsonl
```

## 8. Candidate card

- **Mechanism.** Rank-r aggregate linked module (two degree-d fields, r ≪ d linked generators,
  height N^{1/d}), no rational side, dependence among the r−1 linking conditions (singular/frame-less).
- **Model.** Single N, RAM, B² linear algebra, smoothness heuristic H; `resbarrier`/`rankdegen` box
  scope.
- **Result.** Typical-size-model floor L[1/4] (AM–GM proved; the x^{1/2} exponent is the heuristic
  degree term; no rigorous positive floor < 1/3, resbarrier §2.6). Under G_r, L[1/3], r = 2 optimal.
  Nonsingular configurations give α ≥ 1/3 (minor-product penalty). A rank-discrepant mod-N intertwiner
  splits N (Lemma SPL-INT, narrow). Frame-less corner open.
- **Failure scope.** No construction known; G_r unproven; LP_exc unproven; frame-less area theorem
  unknown; SPL-INT does not cover prime-uniform or global-non-power dependence. Not a candidate
  algorithm.
- **Next step.** Either a frame-less (singular, Δ=0) area/volume theorem (would close the hole to
  L[1/3]), or an explicit global-but-non-power dependence at height N^{1/d} (would open it toward
  L[1/4]). We found neither.

## 9. Reviews

- **Codex blind `midk_1` (`codex/midk_1.out`):** independent and corroborating. Reproduced the
  fallback optimum α = 1/4 at δ=1/2, ρ=1/4 and α ≥ 1/3 under the independence penalty; added the
  interpolation **α = 1/(4−γ)** for a partial penalty N^{−r^γ}; gave the explicit cyclotomic
  automorphism-graph construction (dependent links at height 1, useless — "the modular square roots
  agree … no nontrivial gcd"); and gave the deterministic minor-product height penalty
  h₂ ≥ (r−1)x/(D−1) with the clean "no intermediate prize" conclusion under nonsingularity (all
  maximal minors nonzero), isolating the singular/frame-less corner as the remaining hole. It declined
  to attribute a universal "G_r" to the literature; the hypothesis is stated here as resbarrier's.
- **Opus proof critic (read-only; Lemma LEDGER, Lemma SPL-INT):** ACCEPT WITH FIXES on both;
  all applied. (1) The L[1/4] figure is the typical-size *model* floor, not rigorous — only AM–GM is
  proved; relabelled in §0, §2, §6, §8. (2) SPL-INT's case-(ii) hypothesis as written was vacuous;
  the correct condition is **rank-discrepancy across the prime factors**, "non-scalar" is insufficient,
  and the E2 construction was tautological — replaced by a genuine rank-discrepant intertwiner; the
  "factoring-equivalent" consequence demoted to an assessment with gaps (1)–(3) named. (3) Confirmed
  Families B and C never beat L[1/3] at any finite intermediate value. Minor: c > 0 (not c ≥ 1);
  r=2-optimality cited from the KL ledger, not the monotonicity test.
- **Codex adversarial `midk_2` (`codex/midk_2.out`):** MAJOR REVISION; all points align with the Opus
  critic and are applied. (1) MAP's "exactly three" is an **assessment of completeness, not a theorem**
  — a formal algorithm class + exhaustive reduction would be needed; Coppersmith's multiple-polynomial
  NFS adds a polynomial-count parameter (improves constants only, no 4th sub-1/3 family); C's block is
  deterministic (Coppersmith), so strictly it fails MAP's "counting-only" criterion and is listed for
  completeness. For B/C, *growing* j can reach the L[1/3] scale but that is not an *intermediate* value
  and needs uniform cost bounds — consistent with "no finite intermediate value beats 1/3". (2) LEDGER:
  the 1/4 is valid for the displayed proxy only; a coefficient-height envelope bounds norms from
  **above**, so a pointwise lower bound and a count of useful exceptional candidates are missing — this
  is why the figure is a model floor (already relabelled). (3) SPL-INT: supplied the three-outcome
  correction (outcome 0 added; N=21 example) and confirmed outcome (i) gives no factor, so
  factoring-equivalence of a general exploitable dependence is **not** proved — matches the demotion to
  assessment. The dependence-source classification (I)/(II) "needs a theorem" (a linear relation need
  not be a ring power; a nonlinear correspondence need not linearise) — flagged in §3.1's gaps.
