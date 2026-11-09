# completeness_audit: a cheap completeness audit of the barrier catalogue (phase 5 meta-audit)

Agent `caudit`, 2026-11-09. Reviewer item #3 ("a cheap completeness audit").

**Question.** The project has dozens of barrier theorems, each ruling out a class of routes to sub-L[1/3] (or
sub-2/5, sub-3/8) factoring *in a restricted model*. Taken together, do they cover ALL routes, or are there
model-classes / paradigms / parameter corners that no theorem touches and that were quietly dropped rather than
closed? Overriding user question: **is it possible to beat L[1/3]?**

**Files.**
- Code: `src/caudit/coverage.py` (exact checks Q1–Q5; reuses `src/suffmap/principle.py` by import).
- Tests: `tests/test_caudit_coverage.py` (6 tests, exact arithmetic, no leakage).
- Raw: `experiments/raw/caudit_coverage_checks.json`.
- GPT-6.1-sol: blind round `codex/caudit_1.{md,out}` (name every route NOT closed); adversarial round
  `codex/caudit_2.{md,out}` (attack this coverage matrix). §8.
- Opus critic: read-only, on the single most credible uncovered route. §9.
- No experiment tests a prediction (the checks in §7 are exact formula checks), so there is no pre-registration,
  per the task rule and the `suffmap` precedent.

**Labels.** PROVED (in a stated model) · model-PROVED · HEURISTIC · OPEN · UNPROVEN EXTRAPOLATION · assessment.

---

## 0. Verdict (plain language)

**The catalogue is NOT literally complete, and it does not claim to be.** It is a patchwork of ~40 model
theorems plus roughly a thousand explicitly-flagged "open / not covered / outside the model" remarks in the notes
(1037 such lines by grep). Every barrier is a MODEL theorem, never a factoring lower bound, and the gap between
each model and a general algorithm is stated honestly in the notes themselves. So "the theorems jointly rule out
all routes" is false as a literal claim.

**But the gaps do not translate into a credible way to beat L[1/3].** Every uncovered corner we found is one of:
1. **endpoint-dominated** — even if the mechanism worked perfectly it only improves the L[1/3] *constant*
   (Prop SI → 1.2114; MNFS/factory → 1.6386), or it costs N^{Θ(1)} to cash in (every leak delivers O(1) bits of a
   factor, and n/4 bits are needed: exact §7 Q4); or
2. **secretly factoring- or DL-equivalent** — the doubly-smooth lister, Cayley-graph navigation, hidden-root
   short-quotient navigation (navigate Theorem DLB), #E(Z/N), class groups; or
3. **a genuine coverage gap with NO candidate algorithm** — the cross-family union (§5.1) and the PRIM nested-gate
   frontier (§5.4) are not closed by any single theorem, but nobody has a mechanism that lives there.

**The single most important structural weakness** is that the catalogue is a union of *single-family* barriers
proved in *disjoint* models, and it contains **no compositional, amortized, or adaptive-conditioning theorems**.
Concretely: (i) the leak barriers (PT/NE1/PT-J/AC0-ε) and the relation barriers (SS/heurnfs/listing) have no
common model, so no theorem bounds a hybrid (`leakmap` §8: "unions ... are not covered, and representatives make
any union containing them complete"); and (ii) every relation barrier is for *oblivious* candidate generation and
*completed* candidates — the cost of an *adaptive* search that builds the smoothness certificate jointly with the
candidate is unbounded by anything in the repo. All three independent reviewers (Opus critic §9, both GPT-6.1-sol
rounds §8) converged on this: "the catalogue lacks the conditional, amortized, and compositional theorems needed
to turn its individual barriers into the claimed universal prohibition" (caudit_1).

**Where the believer should push, corrected by the critic.** The cross-family *leak→relation* hybrid (our first
draft's top candidate) is **demoted**: the Opus critic showed the leak component is *vacuous* — a leak is just an
N-dependent candidate-selection rule (`leakmap` §1), and resbarrier's Corollary RB2 / Lemma C2(e) already bound
*every* rank-2 selection rule and every N, with no counting heuristic; "candidates smaller than the resultant
floor" is proved impossible at rank 2. So §5.1 reduces to the relation-method heuristics that resbarrier already
lists as uncovered — **H5** (sparse selection from boxes larger than B; not covered) and **G_r** (rank ≥ 3 with no
rational side; closed only under Schmidt's heuristic) — plus the **adaptive-constructor cost lemma** (§5.2). The
aggregate rank-r≥3 family is the one place where the size-entropy law *formally* permits sub-1/3 (it realizes the
a=½,b=0 law → α=¼, §7 Q1) and the barrier is explicitly heuristic, not a theorem — but it has no known
construction and is believed false for typical N.

**Answer to the user's question.** On the current evidence, **no** — we found no uncovered route that credibly
beats L[1/3]. The place with a real algorithm already outside a proved barrier is **MNFS / many-factor-base NFS**
(outside model S): provably (heuristically) exponent-1/3, constant-only (§7 Q2). The single most credible
*genuinely uncovered* corner is the **rank-r≥3, no-rational-side linked-module family (resbarrier G_r)** together
with **H5** — a relation-method heuristic, not a new paradigm and not a leak. Subjective probability that any item
here yields sub-L[1/3]: **< 1%** (the reviewers decline to put a number on it, noting the missing compositional
theorems; we keep <1% as an assessment, concentrated on G_r/H5 and the adaptive-constructor cost lemma).

---

## 1. Method

1. Read the barrier map in full (`BARRIERS.md`, `APPROACHES.md`), the three survey syntheses
   (`connections_results`, `native_results`, `leakmap_results`), the size–entropy note (`suffmap_results`), and
   `STATUS.md`.
2. Built the coverage matrix (§2): every distinct route family → the exact theorem that closes it → the MODEL it
   closes it in → the gap between that model and a general algorithm → status.
3. Census of open remainders (§3): grep for `open|OPEN|not covered|outside the model|does not imply|not
   excluded|extrapolat` across `notes/` (1037 lines), deduplicated into the recurring gap-types.
4. Looked specifically (§5) for routes impossible to express in any covered model: non-oblivious/adaptive;
   N-via-an-auxiliary-object-of-size->N; cross-family hybrids; super-constant-k / many-factor-base (what the
   size–entropy law does NOT constrain); and whether single-family barriers imply a cross-family barrier.
5. Exact sanity checks (§7) on the quantitative claims that decide plausibility.
6. GPT-6.1-sol blind + adversarial (§8); Opus critic on the top uncovered route (§9).

---

## 2. Coverage matrix

Columns: **route family** · **closing theorem(s)** · **model it closes it in** · **gap to a general algorithm** ·
**status for beating L[1/3]**. "Closes" always means *in the stated model*; none is a lower bound.

| # | Route family | Closing theorem(s) | Model | Gap to a general algorithm | Status |
|---|---|---|---|---|---|
| 1 | Single-family sampling / sieve (Dixon, QS, NFS, SNFS) | Size–entropy law SS (`suffmap`) | model S: one shared factor base, O(1) smoothness integers/candidate, tested in full, ω≤2 | **SS is a reduction, not a prohibition.** It says α=a/(2−b); sub-1/3 needs a family with a<1/2 (b=½) or b<½. It does NOT forbid such a family — §7 Q1, the stylised rank-r law formally gives 1/4 | endpoint: no sub-1/3 family realized |
| 2 | Multi-family / MNFS / factory (shared rational side, V fields, many factor bases) | heurnfs (MNFS 1.9019; factory 1.6386) | heuristic NFS model, validated on CADO | **Outside model S entirely** (SS does not cover it). Covered only by heurnfs, a heuristic; exponent-1/3 robustness rests on that | constant-only (exp. 1/3), §7 Q2 |
| 3 | Smooth-number *listing* (doubly / one-sided / Bohr-set) | Theorem LS; Prop DLT | listing model (cost per output) | doubly-lister ⇒ poly time but **factoring/DL-equivalent**; one-sided costs L[1−α]/output | equivalence-closed |
| 4 | Construct smooth numbers near kN (Prop SI) | Prop SI; sihard Theorem J; siconstruct Prop CV | oblivious constructor; promise-CVP at factor 1+o(1) | **adaptive / non-oblivious constructors NOT covered** (sihard, sinonobliv, smoothbuild all flag this); but navigate: hidden-root navigation is FACT∧DL | constant-only (1.2114) even if it worked |
| 5 | Coverage / local designs (deterministic, exp. regime) | T1/T2, multimod A/A′/B, cycloharvey LB, dchworst W, cov38b B | 𝔐⁺ / identity-model covering, labels charged | true-factor floors N^{1/6} (full coverage) / N^{1/5} (one modulus). **Average-case & arbitrary N-dependent placements not covered** (coverage38); worst-case <N^{1/5} needs 7/25 saving | exponential regime; not an L[1/3] route |
| 6 | Descent ladders (rigorous randomized) | Theorem U (2/5), 2L-final (3/8), nested NR | explicit-error degree-uniform estimates; union-bound multilevel | depth j → (j+1)/(3j+2)→1/3 but **depth ≥4 fails** (precision < covariance); L[1/3] needs ERH + all-target supply | reaches 3/8 rigorously; 1/3 open & conditional |
| 7 | Floor-sum / analytic / short-circuit | floorsum/topblock (=Coppersmith known-bits), analytic, slpcirc SP, periodfind T | oblivious floor-sum ladder; sparse circuits; Fourier observables | degree-j ladder = Coppersmith barrier; **top blocks OPEN**; general oracle NP/#P-hard | no mechanism; open corners exponential |
| 8 | Representation / self-reduction | factorbits (prime-exp addition ⇔ FACT), selfred SR1/SR2, Theorem PT | algebraic self-reduction classes; modulus-as-modulus | **no modulus-changing self-reduction in algebraic classes**; general poly-time maps not covered (complete) | equivalence-closed where covered |
| 9 | Period-finding / Fourier side | periodfind Theorem T | observable g with g(a^x) Fourier mass | only large-order multiplicative-character observables peak → back to smoothness | closed to smoothness |
| 10 | Leaks without smoothness (make p,q differ w/o splitting) | leakmap: PT-X, NE1, NE1-R, PT-J, AC0-ε; genusleak GB; newfamily | decision trees over fixed Z/N-algebras; N-exponent gates; AC⁰ on bits; Jacobi of pure values | three families (M/A/Q) only; **nested/retained gates (PRIM) OPEN; representatives are complete; unions not covered** | §5.1, §5.4 — open, but leaks are O(1)-bit |
| 11 | Native invariants as new views | nativeinv TRI, nativemod, nativeana LIN | quasi-poly local data; linear analytic estimators; powering in fixed algebra | **single a_f(N) of a fixed non-CM form OPEN** (no reduction, no barrier); nonlinear/adaptive estimators not covered | open view, native ≥ N^{1/2} conjectured |
| 12 | Unexpected connections (40 fields) | connections survey (7 mechanisms) | per-connection named failure step | distributional / aggregate / restatement — each closed at a step, none a theorem | no route |
| 13 | **Structured lattice / Diophantine constructors** (Schnorr–Adleman prime-log lattice; simultaneous approximation U/V≈N) | none complete | — | **The cost lemma is MISSING** (caudit_2 #3, Opus critic): no barrier bounds the amortized cost of producing enough independent useful smooth triples via lattice solving. Best endpoint = Prop SI 1.2114; schnorrcvp refutes Schnorr 2021 only at γ→1 | open construction problem, endpoint-dominated (1.2114), DL-adjacent |

---

## 3. The open-remainder census (what was flagged open and never closed)

Grep over `notes/` returns 1037 lines matching `open|not covered|outside the model|does not imply|not excluded`.
Deduplicated, they fall into **seven recurring gap-types**, each appearing across many agents:

- **G-ADAPT — adaptive / non-oblivious constructions.** Repeatedly dropped, never closed: adaptive CRT decoding
  far below the Johnson radius (`sinonobliv`, `crtdecode`); adaptation driven by many small residue queries whose
  CRT combination is uncontrolled (`below16` §3.3, line 131/665); adaptive target-dependent growing products
  (`smoothbuild` line 61/318); nonlinear & adaptive analytic estimators (`nativeana` line 55); adaptive selection
  among sub-rules using other data of N (`inversesieve` line 166); random thinning / adaptive small prime support
  (`siconstruct` line 178). **Every covered barrier is oblivious.**
- **G-CROSS — cross-family / multi-field hybrids.** `leakmap` §8 (unions not covered; representatives complete);
  `fruit` line 51/486 (multi-field hidden-root, batch-/field-level cross-field elimination not covered);
  `escapeR2` §2.3 (tied fields: parity transfer exists, usefulness open; directed recombination not covered).
- **G-SIZE — outside the size-entropy model S.** many factor bases; super-constant k (the listing boundary);
  non-sampling; MNFS/factory (`suffmap` §1.1 "What S excludes").
- **G-SUPN — N used through an auxiliary object of size > N.** arithmetic-surface / adaptive algorithms not
  covered (`resbarrier` line 383); rank-≥3 / no-rational-side corner closed only heuristically (`resbarrier`,
  `BARRIERS.md`); working mod N^t (`trackBD` W-4 covers the derivative case only).
- **G-PRIM — nested/retained N-exponent gates.** `leakmap` §3.6, §8 (the class-G frontier; deterministic depth
  ≥2 in `prim`).
- **G-NATIVE — a single modular coefficient.** `cuspcoef`/`nativemod`: a_f(N) of a fixed non-CM form is OPEN
  (no reduction, no impossibility, none in the literature).
- **G-ENDPT — explicitly constant-only openings.** sub-quadratic sparse GF(2) LA (MNFS→1.8595); percolation
  relation selection; near-regular NFS matroids. All flagged as constant-level and expected dead.

**Reading.** G-SIZE, G-SUPN and G-ENDPT are endpoint-dominated or exponent-1/3-robust (§7). G-NATIVE is a "new
view," conjectured native ≥ N^{1/2}, not a route. **G-CROSS collapses** under the Opus critic (§9): the leak part
of a leak→relation hybrid is a vacuous selection rule, and the residue is the relation-method heuristics H5/G_r,
which belong to G-SIZE/resbarrier, not to a new paradigm. The gap-type genuinely *not* reducible to a costed
argument is **G-ADAPT** (adaptive certified-relation construction): oblivious/promise constructors are bounded,
adaptive ones are not — but the best achievable endpoint is still the L[1/3] constant 1.2114 (Prop SI), and the
all-N solver is DL-equivalent (navigate DLB). So even G-ADAPT is endpoint-dominated *if it works at all*; what is
missing is the cost lemma, not a sub-1/3 exponent.

---

## 4. Does the union of single-family barriers imply a cross-family barrier?

**No (UNPROVEN EXTRAPOLATION).** This is the central completeness finding. Formally:
- The **leak** barriers live in the decision-tree / fixed-Z/N-algebra / AC⁰ models (`leakmap`). Their conclusion
  is "non-split advantage is negligible for target F."
- The **relation** barriers live in model S (`suffmap`) / the heuristic NFS model (`heurnfs`) / the listing model
  (`pilist`,`smoothenum`). Their conclusion is "cost ≥ exp(C*) / ≥ MNFS."
- There is **no single model containing both**, and no theorem of the form "if family-1 is barred in model M₁ and
  family-2 is barred in M₂ then their composition is barred." `leakmap` §8 states the union is uncovered and that
  integer representatives (which every real relation method uses) make any union containing them *computationally
  complete* — i.e. a classification there would be a classification of all factoring algorithms. So the
  cross-family barrier is not merely unproven; in the complete model it is *equivalent to the whole problem*.

This is exactly the "patchwork" weakness a reviewer would flag. The mitigation (not a proof) is §7 Q4: the only
known leaks deliver O(1) bits of a factor, which is endpoint-dominated.

**Partial closure found by caudit_2 (bounded advice).** If the leak stage communicates only **t = O(log log N)
bits in total** (not per query), enumerate all 2^t transcripts; the overhead is polynomial in the input length,
so *if every fixed-transcript continuation satisfies the relation barrier's hypotheses uniformly*, the hybrid is
closed by this standard bounded-advice argument. This genuinely closes the *total-O(log log N)-bit* hybrid. It
does **not** close O(1)-bits-*per-query* hybrids that accumulate or adapt to candidates, nor does it bound an
*efficient conditioned constructor* that produces biased candidates directly (rejection sampling pays for the
filter's rarity: Pr(S|E) ≤ Pr(S)/Pr(E), so a rare filter needs 1/Pr(E) work — caudit_2). And caudit_2's caveat
is fair: Coppersmith's n/4-bit threshold (§7 Q4) is a *sufficient* recovery condition, not a *necessary*
information threshold, so "O(1)-bit ⇒ endpoint-dominated" is an argument about *contiguous* factor bits, not a
theorem about all partial information.

---

## 5. The most credible uncovered routes, with sketches and exponents

Ordered most-credible-first. "Credible" = not obviously closed; it does NOT mean we believe it beats L[1/3].

### 5.1 Rank-r≥3 / H5 relation-method heuristics (was: cross-family leak hybrid) [DEMOTED by the Opus critic]
- **Original sketch (demoted).** Use a cheap leak to *bias* candidate generation so the effective size-law has
  a<½, producing candidates smaller than the resultant floor conditioned on the leak.
- **Why the LEAK part is VACUOUS (Opus critic §9; caudit_2 agrees).** A leak is just an N-dependent
  candidate-selection rule (`leakmap` §1: "'leak' is computational"). resbarrier **Corollary RB2** (§2.3) and
  **Lemma C2(e)** (§2.3, l.211ff) bound *all* box points for *every* rank-2 method and *every* N "with no counting
  heuristic for the polynomials", so no selection rule — leak-guided or not — produces more small-norm candidates
  than C2 allows; "a<½ at rank 2" is *impossible* under H+H5. The premise was also mislabelled: suffmap's H_S
  exclusion is about pre-filtering on *smoothness*, not about choosing *smaller* candidates; a size-biased
  sub-family is still covered by SS and C2(e). An ℓ-bit external leak just picks one of 2^ℓ leak-free generators
  (worth only advice; closed for ℓ=O(log n) by §4's bounded-advice argument).
- **What actually survives, with no leak (the real residue).** resbarrier's own uncovered hypotheses:
  - **H5** — sparse selection from boxes/lattices of side > B^{1+o(1)}: *explicitly not covered* (needs an
    X-dependent uniform-Thue error term, Thunder's 𝓗(F)^{D−2}) (`resbarrier` §2.5 l.231–235, §2.7 H5).
  - **G_r / H4** — rank r≥3 with no rational side: the **aggregate rank-r family** (`resbarrier` §4.2) formally
    realizes the a=½,b=0 law → **α=¼** (§7 Q1), and is closed *only* under Schmidt's heuristic G_r (expected count
    N^{e}, e≤1−(r−1)(1−1/d)<0). **No construction known; believed false for typical N; "open as a theorem".**
  - **H / H6** — smoothness-biased selection at fixed size: maps to *partial listing* (smoothenum T2–T7), constant
    only for θ<1, and sub-1/3 needs doubly listing (1−θ→0) which is DL-hard via Prop DLT.
- **Exponent if G_r/H5 could be beaten.** α=¼ formally (the stylised rank-r law). So this *is* the one corner
  where the size-entropy law permits sub-1/3 — but only conditional on constructing a family nobody can construct.
- **Verdict: a real open THEOREM (G_r, H5) but NOT a candidate algorithm, and NOT a leak phenomenon.** The
  single-most-credible-structural-hole billing from our first draft was wrong; the honest statement is that the
  relation-method lower bound is heuristic exactly at H5 and G_r, which the project already lists.

### 5.2 Adaptive / non-oblivious smooth-number construction near kN [G-ADAPT; OPEN, DL-equivalent]
- **Sketch.** An adaptive constructor that, given failed attempts in intervals [kN,kN+h], steers the next query
  (residue-guided beam search over divisor products, `crtdecode` line 397) to a smooth u, beating the
  oblivious-constructor L[2/3]/output floor.
- **Why uncovered.** sihard Theorem J and Prop PC bound *oblivious* and *promise* constructors; adaptive
  target-dependent products are explicitly not covered (smoothbuild, sinonobliv).
- **Exponent if it worked.** Prop SI: L[1/3, 1.2114]. **Constant-only even in the best case.**
- **Verdict: illusory for sub-1/3 (endpoint-dominated) and probably DL-equivalent.** navigate Theorem DLB:
  a solver for all N gives prime-field DL in L[1/3,1.80]. So it is both constant-only and an equivalence.

### 5.3 Many-factor-base / super-constant-k sampling [G-SIZE; covered by heurnfs, constant-only]
- **Sketch.** Leave model S by using M→∞ factor bases (MNFS) or k→∞ smoothness integers per candidate (the
  listing boundary), to lower the effective b.
- **Why "uncovered" by SS.** SS explicitly excludes M>1 and k≍u.
- **Verdict: not actually uncovered.** heurnfs optimizes over all M and finds exponent 1/3, constant 1.9019
  (factory 1.6386). §7 Q2 reproduces the constants. **The exponent is robust; this is the clearest "gap that
  isn't."**

### 5.4 PRIM nested/retained N-exponent gates [G-PRIM; OPEN, O(1)-bit]
- **Sketch.** A program computing e.g. (a^N+1)^N over a fixed Z/N-algebra, whose zero set in r=q mod(p−1) is
  large and not explained by order alignment or a short-interval pin (leakmap §3.6).
- **Why uncovered.** leakmap proves only the one-layer (NE1) and terminal-Jacobi parts of PRIM.
- **The displayed example is itself closed (caudit_2, exact).** For balanced odd primes gcd(N,(p−1)(q−1))=1
  (q∤p−1, and p|q−1 forces q=p+1), so x↦x^N permutes each residue field and gcd((a^N+1)^N, N)=gcd(a+1,N). Harness
  check `caudit` diagnostic: 0 mismatches in 12000 cases (p,q,a). Nested (x^E)^F collapses to x^{EF}; the *additions*
  are what could matter.
- **Verdict: illusory for sub-1/3.** Even a perfect PRIM leak is a *leak* (O(1) bits / a predicate), not a
  relation engine; cashing it in is §7 Q4 (N^{1/4}). A real open model-classification question, not a route.

### 5.5 A single modular coefficient a_f(N) as a native view [G-NATIVE; OPEN]
- **Sketch.** Compute one coefficient a_f(N) of a fixed non-CM cusp form at composite N below L[1/3] without
  separating σ(N).
- **Why uncovered.** cuspcoef/nativemod: no reduction (it is in ZPP^FACT, NP-hardness would give NP=coNP under
  GRH) and no impossibility; none in the literature.
- **Verdict: open but not a route.** It is a *new view*, not an algorithm; the native cost is conjectured
  ≥ N^{1/2}. It cannot by itself beat L[1/3]; at best it is a target for future reduction work.

---

### 5.6 Prioritized list (the deliverable)

**Real candidates for an uncovered gap (ranked), after the critic + both codex rounds:**
1. **G_r + H5 — rank-r≥3, no-rational-side linked-module family, and sparse selection from large boxes.** The one
   place the size-entropy law *formally* permits sub-1/3 (α=¼) and the barrier is explicitly heuristic, not a
   theorem (resbarrier §2.7, §4.2). *No construction known; believed false for typical N.* This is a missing
   THEOREM, not a candidate algorithm.
2. **G-ADAPT — adaptive construction of certified relations** (structured lattice/Diophantine constructors;
   matrix row 13). The cost lemma is genuinely missing (all three reviewers). *But* endpoint-dominated at 1.2114
   (Prop SI) and the all-N solver is DL-equivalent (navigate DLB) — so uncovered as a *bound*, not as a sub-1/3
   *route*.
3. **The missing compositional / amortized theorems themselves** (§4). The catalogue's completeness defect is
   real and reviewer-confirmed; closing it is "the whole relation-method lower bound" (critic), beyond reach.
4. **a_f(N) native view** (G-NATIVE) and **PRIM nested-additive gates** (G-PRIM): genuinely open
   model-classification questions, but neither is a sub-1/3 route (native ≥ N^{1/2}; leaks are O(1)-bit).

**Illusory (secretly equivalent or endpoint-dominated) — demoted:** the cross-family leak hybrid (leak vacuous,
§5.1/§9); MNFS/factory/many-bases (exponent-1/3, §7 Q2); doubly-smooth lister / S-unit / Cayley navigation / hidden-root
navigation (factoring/DL-equivalent); the displayed PRIM example (identity, §5.4); super-N auxiliary objects
(size enters λ, worsens the exponent, §8 caudit_1 #5); Pochhammer SLP (τ-conjecture, landmark not route).

## 6. Routes that are illusory (secretly equivalent or endpoint-dominated) — the demotions

| Candidate uncovered route | Why it is illusory |
|---|---|
| Adaptive smooth-number constructor (5.2) | constant-only (1.2114) AND DL-equivalent (navigate DLB) |
| Doubly-smooth lister / S-unit / Cayley navigation | poly-time but factoring/DL-equivalent (Prop DLT, E3) |
| Hidden-root short-quotient navigation | FACT∧DL (navigate Prop H), DL in L[1/3,1.80] |
| MNFS / factory / many factor bases (5.3) | exponent 1/3, constant-only (heurnfs; §7 Q2) |
| PRIM nested gates (5.4) | O(1)-bit leak; cashing in is N^{1/4} (§7 Q4) |
| Pochhammer/factorial SLP | P/poly but contradicts τ-conjecture; landmark either way, not a route |
| #E(Z/N), class groups, φ/σ/λ | factoring-equivalent (E7) |
| Super-N auxiliary object (mod N^t, big variety) | size enters λ; SS/resbarrier charge it; endpoint-dominated |

---

## 7. Exact checks (`src/caudit/coverage.py`; `experiments/raw/caudit_coverage_checks.json`)

All exact (sympy/fractions); floats only for display. 6 tests pass.

- **Q1 — the sub-1/3 wall.** α=a/(2−b); α<1/3 ⟺ 3a<2−b. Every *realized* construction (Dixon/QS at (1,0)→½;
  GNFS/SNFS/t=3 at (½,½)→⅓) sits at or above 1/3 and passes the balance identity C²=2λlogλ. **But the stylised
  aggregate rank-r law (½,0) formally gives ¼ and also passes** — so **Theorem SS does not forbid sub-1/3**; it
  reduces it to constructing a family with that law, whose absence is resbarrier's heuristic G_r, not a theorem.
  This is a genuine subtlety the matrix must state plainly.
- **Q2 — outside model S is constant-only.** GNFS (64/9)^{1/3}=1.9230; MNFS ((92+26√13)/27)^{1/3}=1.9019;
  factory 1.6386. All exponent 1/3. The one regime with a real algorithm outside a proved barrier buys a constant,
  never the exponent.
- **Q3 — Prop SI.** (16/9)^{1/3}=1.21141 (matches suffmap's `single_interval_program`), exponent 1/3, below MNFS.
  The cleanest open sufficient condition is still an L[1/3] constant.
- **Q4 — leak-exploitation floor.** Coppersmith: ⌈n/4⌉ known top bits of p factor N (tight worst case). A leak of
  O(1) or O(log n) bits per event needs N^{Θ(1)} events — exponential. This is why §5.1/5.4 are endpoint-dominated.
- **Q5 — cross-family union.** Logical check: leak model and relation model are disjoint; no theorem bounds their
  composition; representatives make the union complete. Status: UNPROVEN EXTRAPOLATION.

---

## 8. GPT-6.1-sol rounds (`codex/caudit_1.out`, `codex/caudit_2.out`)

**caudit_1 (blind — "name every route it does NOT close").** Converged on our verdict independently: the
catalogue's "surveyed without finding an algorithm", "factoring-equivalent", "contradicts τ-conjecture" are not
unconditional lower bounds, and a leak classification does not bound every algorithm using those leaks. Its three
best directions: (1) prove *fruitfulness* for a randomized NFS variant (rigorous L[1/3], beats the rigorous 2/5
and 3/8, not heuristic NFS — matches the project's own open problem #4); (2) **adaptive construction of certified
relations** (its strongest strictly-sub-1/3 direction, but "no established sub-1/3 estimate"; our SI endpoint caps
it at 1.2114); (3) depth-two retained PRIM programs. Exact points we fold in: the descent formula
(j+1)/(3j+2)=1/3+1/(9j+6) *never crosses 1/3*; and for an auxiliary object with log M=(log N)^{r+o(1)}, an
L_M[α] method has N-exponent rα, so merely enlarging N worsens the exponent (super-N objects need extra structure,
not size). Its closing line: the catalogue "lacks the conditional, amortized, and compositional theorems needed to
turn its individual barriers into the claimed universal prohibition." **No strict sub-1/3 candidate earned a
favorable exponent.**

**caudit_2 (adversarial — attack this matrix).** Two concrete contributions folded in above: (i) the
**bounded-advice partial closure** of the cross-family hole (§4): total t=O(log log N) leak bits ⇒ enumerate 2^t
transcripts; but O(1)-bits-per-query does not close, and an efficient conditioned constructor is exactly what is
not bounded. (ii) the **exact identity closure** of the PRIM example (§5.4). It also flagged: do NOT close the
a_f(N) route via Bach–Charles (that uses both a_f(N) and a_f(N²)); and it independently asked for matrix **row 13**
(structured lattice/Diophantine constructors), whose cost lemma is missing. It found no defensible sub-L[1/3]
algorithm and declined to certify "below 1%".

## 9. Opus critic on the first-draft top route (§5.1) — DECISIVE

A read-only Opus critic was tasked to show the cross-family leak→relation hybrid is secretly closed/equivalent. It
**succeeded in defeating the "leak" framing** (folded into §5.1 and §0):
- A leak is a computable N-dependent selection rule (`leakmap` §1); resbarrier RB2 (§2.3) + Lemma C2(e) bound
  *every* rank-2 selection rule and every N with no counting heuristic, so "candidates smaller than the resultant
  floor" is *proved impossible at rank 2*. The leak component "adds nothing."
- Even a perfect leak only loosens N²|disc (k=2) to p²≈N (k=1): constant (16k/9)^{1/3} moves 1.526→1.2114, α≥1/3
  for every fixed k. Many-bit leaks are already factoring (GAR 3n/4 bits; Coppersmith n/4).
- The H_S exclusion is about smoothness pre-filtering, not smaller candidates — the audit's premise was
  mislabelled.
- **What genuinely survives, none needing a leak:** H (smoothness-biased selection → partial listing → DL-hard via
  DLT), **H5** (sparse selection from boxes > B, not covered), **G_r** (rank≥3 no rational side, heuristic only),
  non-commutative rings. All already listed by resbarrier §2.7 / selfred SR3.
- **To close §5.1 a theorem would need:** for every family in the linked-order model, every rank, and every
  poly-time N-dependent selector, smooth yield ≤ count(Y)·ρ(u_Y)·L^{o(1)} — i.e. prove H for all efficient
  selectors AND remove H5 and G_r. "That is the whole relation-method lower bound, beyond current reach (SR3)."
- Recommendation (applied): reclassify §5.1 as "leak vacuous; reduces to H, H5, G_r"; fix the sketch wording.

The critic's verdict stands: the route does not survive as an independent hole; its residue is the already-known
H5/G_r relation-method heuristics. Per the task rule, no *new* uncovered route is reported as credible that failed
the critic — and this one did fail, so it is reported as demoted, not credible.

---

## 10. Reproduction
```
.venv/bin/python -m pytest -q tests/test_caudit_coverage.py
.venv/bin/python -m src.caudit.coverage > experiments/raw/caudit_coverage_checks.json
```
