# auditRA: is Theorem RA / Corollary P true?

Agent `auditRA`, phase 11, 2026-10-20. This note audits `notes/recurse_results.md` §1, which contains Theorem RA,
Lemmas OS, W, F and T_rec, and Corollary P. Those results are built on Theorem A (`notes/thmA_proof.md`,
`notes/thmA_audit.md`, `notes/lem85_results.md`).

Files:
- This note: `notes/auditRA_audit.md`.
- `notes/auditRA_blind_claims.md` is the blind excerpt given to every reviewer: §1 of recurse_results verbatim, with the review sections and verdicts
  removed.
- `src/auditRA/constants_check.py` writes the exact constants to `experiments/raw/auditRA_constants.json`. It was written from scratch and does not
  import the recurse module; it agrees with it.
  Command: `.venv/bin/python -m src.auditRA.constants_check experiments/raw/auditRA_constants.json`.
- Codex blind referee: `codex/auditRA_1.{md,out}`. The full report is copied to `notes/auditRA_codex1_report.md` from
  the git-ignored `codex/scratch/auditRA_1/auditRA_referee_report.md`.

## Verdict

**Theorem RA and Corollary P are correct consequences of Theorem A, conditional on Theorem A's proof being correct.
The two whole-claim verdicts are RA: ACCEPT WITH FIXES and P: ACCEPT WITH FIXES.** I reached them myself, and so did
all three blind reviewers (1 GPT-6.1-sol, 2 Opus). Nobody found a FATAL or MAJOR issue in the recursion layer.

In plain terms:
- **The recursion layer is elementary and holds up.** It consists of three parts:
  - Wald's identity for adaptive sub-calls.
  - An oracle-substitution lemma: any exact factoring oracle with canonical output gives the same outer law.
  - The recurrence F(x) ≤ L^{p}+L^{a}F(√(2x)), solved by an induction whose margin beats the constant at every level.
- **a really is free.** Every a-dependent constant (C₀ ≍ 1/a², C₁ ≍ 1/a, K = 6 + ⌈2ξ/a²⌉, the L^{θ/(4a)} shifts) costs only
  (log N)^{O_a(1)} = L^{o(1)}. Every upper limit those constants meet is a fixed power of N (Y = N^{2ν}, x^{η_D},
  H = N^h). Small a and large a both pass. The places where a = θ/√2 was used are cost remarks only, for example
  "B ≤ S²", which RA replaces by Õ(B + S²).
- **The portfolio is a legitimate single algorithm.** It is a restart schedule, and for each fixed j its expected time is
  O(2^j·E T_j·poly(s₀)). No uniformity in j is needed.
- **The fixes are all local.** None changes c(a), c\* or the limit o(1).
  - Make the parameter sequence explicit, with η₀_j → 0.
  - Remove the analytic threshold N₀, or prove it computable.
  - Weaken the wording on effectivity.
  - State "0.851 at d = 9" as an infimum.
  - Delete one false aside (§1.5, "|u₀v₀| ≍ N").

**What this does NOT establish.**
- RA and P add essentially no new mathematical risk. They do raise the stakes on Theorem A a great deal: a correct
  Theorem A would then imply a rigorous L[1/2, o(1)], not just a constant 0.92.
- All of the depth lies upstream, in the relation-supply estimates FM, Lemma 8.5′ and LD\*, and in LV6.
- Those are internal and unrefereed. Blind AI review, including the five-referee thmA audit, is not peer review.
- The single most likely point of failure in the whole chain is still the relation supply, especially FM's transfer
  to the moving-residue family. That is codex's weakest link, and I agree with it.
- An extraordinary claim that rests on unrefereed internal lemmas should not be stated publicly as a theorem until those
  lemmas are checked by humans.

**Sanity: there is no conflict with known heuristics or barriers (§4).**
- Heuristic fixed-degree NFS gives L[1/2, 2/√(d+1)]. RA's infimum is √(3(1+√2)/(d+1)), which is exactly √(3(1+√2)/4) ≈ 1.346
  times larger at every d. RA therefore claims less than the heuristic. L[1/2, o(1)] is also weaker than heuristic NFS's L[1/3].
- The "constant 1" folklore is not a theorem. It is the random-squares ledger max(a + ξ/(2a), 2a) at ξ = 1/2, sharp for
  random integers by Croot–Granville–Pemantle–Tetali and Balister–Bollobás–Morris.
- That ledger does not bind here. The √N halves are not required to be smooth: they must be products of k fixed-size
  primes, an event of probability L^{−o(1)}, and they are then completely factored. Only the descent values, of total size
  N^θ, enter the matrix.
- This evasion is legitimate *if* FM/LD\* hold.

---

## 1. Issue table

Severity: FATAL breaks the claim; MAJOR is a gap needing a new argument; MINOR has a local fix. All issues are in
`recurse_results.md` unless stated otherwise. I did not edit other agents' files; owners should apply the fixes.

| # | location | issue | severity | found by |
|---|---|---|---|---|
| 1 | §1.1–1.3; thmA_proof §3 l.117, U3 | **N₀ branch.** "N′ < N₀ ⇒ trial division" needs a computable N₀(d, k, μ, η₀, **a**). The notes only assert effectivity; lem85 notes Drappeau's constants are not stated to be effective. RA is fine for any fixed N₀, which can be hard-coded. P, however, needs the map j ↦ N₀(j) to be computable. A wrong N₀ can make E T_j = ∞, because SPLIT may loop forever on a small cofactor. **Fix:** a wrapper that races the core against trial division at every node, with no N₀ (lem85 §7). OS then applies to the uncapped SPLIT core, not to the wrapper's race. Wrapper cost ≤ 2·core cost, and the recurrence only gains a factor 2. | MINOR (MAJOR for P if left unrepaired) | me; all three reviewers |
| 2 | §1.6 | **The parameter sequence of P is under-specified.** "η₀ small and rational" must be η₀_j → 0, for two reasons. First, θ ≥ η₀, so a fixed η₀ gives c\*_j ≥ √((1+√2)η₀); η₀ = 10⁻³ gives 0.049. Second, μ ≥ ν + η₀ fails at fixed η₀ once 1/(d+1) < η₀. Admissibility was checked only for j ≤ 7. **Fix:** take η₀_j = 1/d_j². Admissibility then holds for all j, since the inequalities are monotone in j. My script checks j ≤ 60, critic 1 checks j ≤ 2000, and c\*(j = 60) = 0.239. Also specify a_j as a rational with \|a_j − a\*_j\| → 0. | MINOR | me; all three |
| 3 | verdict, §1.6 | "The rate is not effective" overstates. The argument *supplies* no effective rate, but noncomputability is not proved. **Fix:** say "no effective rate is supplied". | MINOR (wording) | me; codex; critic 1 |
| 4 | verdict table, task statement | "0.851 at d = 9" is inf over (k, η₀). It is not attained by one parameter set: at d = 9, k = 13 one gets c\* = 0.99521. A single algorithm with 0.851 + o(1) at d = 9 exists only as a portfolio over k, by the argument of P. The table is labelled correctly, but the prose is not. | MINOR (wording) | me; codex |
| 5 | §1.5 "The halves" | "Gauss reduction makes \|u₀v₀\| ≍ N" is false as a universal statement: R = 1 gives (u₀, v₀) = (1, 1). It is used nowhere, since the cost proof needs only \|u₀\|, \|v₀\| ≤ √(2N). **Fix:** delete it, or qualify it as typical. | MINOR | codex |
| 6 | §1.3 Lemma OS | Presentation is circular: OS assumes an exact, a.s.-terminating oracle, and for O = FACTOR_rec that is what is being proved. **Fix:** apply OS inside the strong induction on n, with inputs ≤ √(2n) < n. Make canonical output part of the algorithm's specification, not only of the proof. State that nested calls share no mutable state. | MINOR | critic 1; codex |
| 7 | §1.3 accounting | Define F(x) to include n = 1. The node count 2Ω(N) − 1 omits the perfect-power branch, which adds O(log N). | MINOR | codex; critic 2 |
| 8 | §1.2 citations; thmA_proof U3 | The a-freedom list cites hypD l.396 (the outline of Lemma 8.5′). The operative statements are pcest §2.2 and lem85 regime (R). It omits thmA_proof l.206 ("B ≤ S²") and lem85 l.110 (the interface remark "a = θ/√2"; regime (R) holds for any a). U3 should say that N₀ depends on a. | MINOR | critic 2 |
| 9 | inherited (thmA audit item 3) | LD\* and PZ are stated for the implicit S² = (log N)^{C₃}/(ρρ). RA needs the general-S form, γ ≥ 1/(P³ + P/m) without E N ≥ 1. This is already fixed in thmA_proof U5, but not in ldbound l.459. | MINOR (inherited) | codex |
| 10 | §1.3 optimisation | a\* is algebraic, so parameters should be computable rationals. The note offers this already. | MINOR | codex |
| 11 | everything | **Inherited dependence on Theorem A**: FM (incl. Lemma 8.5′ ← Drappeau Lemme 5 ← Harper's unpublished §3.3), LD\*, Kern, L2_fam and LV6. All are internal and unrefereed. RA/P add no new analytic input but amplify the consequence of Theorem A. | not an RA error; the real risk | all |

Checked and found correct, with no issue:
- Lemma W.
- Lemma F, including the margin η′ = η(1 − 1/√2)/4, the scaling L_{√(2x)}^b ≤ e·L_x^{b/√2}, and the extension
  [2, Y] → [2, Y²/2].
- Finiteness by strong induction.
- The count of two O-calls per DESCEND_ONCE, with E τ ≤ L^{a+o(1)}. This includes field selection, factor exits and
  the geometric number of LV6 rounds.
- Charging 2k·T_SQ to every attempt.
- The halves: they are nonzero, and inputs 1, primes, prime powers and tiny values are handled by FACTOR_rec's own
  branches.
- Inner failures and retries are inside T(m).
- T_rec. It holds even if inner transcripts *and running times* are exposed: they are functions of (R, inner coins), and
  r ↦ εr leaves them unchanged.
- The restart computation: phase cost 2^s − 1, Pr[reach s₀ + m] ≤ 2^{−m(m+1)/2}, and the total over 2^{s₀} is ≤ 7.3.
  Only finite E T_j is used, so there is no E[T log T] issue.

## 2. What I checked myself

1. **a-freedom.** I read the statements and the a-dependent lines:
   - FM: pcest l.304; C₁ > ξ_F/(2a) at l.352.
   - Theorem R: C₀ > max(1, θ/(2a²)) at l.496.
   - LD\*/L2_fam: R₀ = B^{C₀} ≥ 64S²L/(ϱϱ), R₀ < Y, at ldbound l.303–304; any fixed C₀ in DS at l.362.
   - lem85 regime (R): any a > 0; K = 6 + ⌈2ξ_max/a²⌉; Drappeau's range (log x)^{c_D} ≤ y; R₁ = y^K ≤ x^{η_D}.
   - The absorption of R₁^{1−α} = (log x)^{O(K)} by H(u)^{−δ}. This works because u/log²u ≫ log log x, with u ≍ a^{−1}(log N/log log N)^{1/2}.

   In every case the constant is fixed once a is fixed, its loss is polylogarithmic, and its upper limits are
   N-powers. Further points:
   - 2^ω ≤ exp(log N·log 2/(aΛ)) = L^{o(1)}.
   - The h-floors (Kern 2ν, level μ − (d−2)ν, window 1 − dμ + ν), μ ≥ ν + η₀ and the θ ≤ 1 size conditions do not
     involve a.
   - The only a = θ/√2 uses are cost remarks: "B ≤ S²" (thmA l.206, lv6 l.567, ldbound l.535) and "2a < (1+θ)/√2".
     At a\*, B ≤ S² holds anyway, since a\*² < θ/2.
   - Large a needs nothing like S² ≥ B. The S² ≥ ρ_r^{−2}L^c condition is used only by Theorem R, which is off the
     critical path.
2. **Recursion.** I re-derived Lemma F line by line, the Wald step, the node accounting (nodes are SPLIT calls, not FACTOR calls,
   so there is no circularity), and the fixed point c = max(p, a/(1 − 1/√2)) with p = max(2a, a + θ/(2a)).
3. **OS and (T).** There is no cap, timeout or race in the specified algorithm: U7, and SQ_SEARCH is a full enumeration of at most C·S²
   points. The one race appears only in the N₀ repair (issue 1), and it sits at the wrapper level, outside SPLIT. A
   single shared random tape is fine, because later bits are still independent and uniform.
4. **Corollary P.** I recomputed the Markov and phase sum exactly, and checked admissibility of the sequence for j ≤ 60 with
   η₀_j = 1/d², and with a fixed η₀ (issue 2). There is no d(N) diagonalisation: each fixed j gives its bound
   separately. So the fixed-d thresholds N₀(j) never need uniformity, only computability (issue 1).
5. **Sanity.** I derived the heuristic fixed-degree NFS constant √(2θ_h) with θ_h = 2/(d+1), compared it with c\*, and ran the
   literature search (§4).

## 3. Reviewer verdicts (all blind)

Every reviewer read `notes/auditRA_blind_claims.md` and was forbidden to open recurse_results.md, any review or audit
note, STATUS.md and codex/. Each took Theorem A at a = θ/√2 as given.

| reviewer | focus | RA | P | weakest link named |
|---|---|---|---|---|
| GPT-6.1-sol `auditRA_1` (xhigh) | everything (items 1–6) | **ACCEPT WITH FIXES** | **ACCEPT WITH FIXES** | Inherited: FM's transfer from restricted smooth-modulus BV to the weighted moving-residue family, where a wrong relative-error or induced-weight estimate would destroy the relation supply. The recursion and restart arguments do not depend on it. |
| Opus critic 1 | analysis and probability (OS, W, F, T_rec, P) | **ACCEPT WITH FIXES** | **ACCEPT WITH FIXES** | P needs N₀(j) computable in j, or the alternation repair. Beyond that, RA rests on Theorem A at every fixed a. |
| Opus critic 2 | parameter ranges and a-freedom | **ACCEPT WITH FIXES** | **ACCEPT WITH FIXES** | a-freedom holds, since all a-constants cost (log N)^{O_a(1)} against fixed N-power margins. So the real weak point is Theorem A itself, and inside RA/P it is the computability of N₀. |
| me | all | ACCEPT WITH FIXES | ACCEPT WITH FIXES | Within RA/P: issue 1 (N₀ computability, repaired by the wrapper). Overall: the relation supply FM/LD\* of Theorem A. |

All three reviewers independently found issues 1 and 2. Codex alone found issue 5. Nobody found anything FATAL or MAJOR.

## 4. Sanity and literature

Sources: my WebSearch, plus the literature sections of the codex report.

- **No rigorous c < 1 is known.** Lenstra–Pomerance (JAMS 5, 1992) remains the rigorous record at L[1/2, 1 + o(1)]. Their p. 509
  remark that re-tuning y cannot improve *their* Algorithm 10.1 is not a general barrier.
  - Lee–Venkatesan (arXiv 1805.08873, JNT 2018) prove rigorous L[1/3] congruences of squares for a randomised NFS,
    but not fruitful factoring.
  - de Boer–Pellet-Mary–Wesolowski (arXiv 2512.01588) is ERH-rigorous for class groups. It does not treat factoring.
  - Nothing found in 2025–2026 changes this.
- **The folklore obstruction.** For relations built from random-like integers of size N^ξ, the ledger max(a + ξ/(2a), 2a) gives
  L[1/2, √(2ξ)], which is 1 at ξ = 1/2.
  - Pomerance's surveys tie the L[1/2, 1] of QS, CFRAC and class-group methods to their √N-size smooth auxiliaries
    (Pomerance, "The role of smooth numbers…", §5).
  - Croot–Granville–Pemantle–Tetali (Annals 2012) and Balister–Bollobás–Morris (Annals 2018) make the random-squares
    threshold sharp, for independent uniform integers.
  - **None of these is a lower bound for algorithms that eliminate large primes by an external descent.** In Theorem A
    and RA, the √N halves must have a k-prime shape (probability L^{−o(1)} at fixed k), not B-smoothness.
  - Only descent values of total size N^θ are combined, so the ledger gives √(2θ) as a *heuristic floor*.
  - RA's √((1+√2)θ) lies above it, because the recursive shape test costs L^{a + c/√2}.
  - The evasion is legitimate, but only as strong as FM/LD\*. Calling a step "descent" does not prove the moments.
- **The heuristic fixed-degree NFS** (balanced non-monic, Buhler–Lenstra–Pomerance §12.9) gives c_h(d) = 2/√(d+1): 0.632 at d = 9, against
  RA's 0.851.
  - The ratio is √(3(1+√2)/4) ≈ 1.346 for every d.
  - Letting d → ∞ gives L[1/2, o(1)] heuristically too, and actual NFS gives L[1/3]. So P claims nothing beyond the heuristic picture.
  - The reason no-one "noticed" is not a barrier. A rigorous fixed-degree relation supply with fruitfulness has never been
    published, and that supply is exactly Theorem A's content.
- **Prior art for the recursion.** Neither codex nor I found a rigorous analysis that recursively factors √N-size auxiliaries with the
  algorithm itself. recurse §4 names Pomerance 1987 (REC applied once to p − 1) as the closest precedent. I did not
  re-verify it.

Sources: [Lenstra–Pomerance 1992](https://www.ams.org/jams/1992-05-03/S0894-0347-1992-1137100-0/);
[Lee–Venkatesan](https://arxiv.org/abs/1805.08873); [dBPW 2512.01588](https://arxiv.org/pdf/2512.01588);
[Balister–Bollobás–Morris](https://arxiv.org/pdf/1608.03857); [CGPT](https://annals.math.princeton.edu/2012/175-3/p10);
[Pomerance, role of smooth numbers](https://math.dartmouth.edu/~carlp/PDF/paper104.pdf);
[van Leeuwen 2007.02689](https://arxiv.org/abs/2007.02689).

## 5. Exact constants (src/auditRA/constants_check.py)

Symbolic results (sympy, exact):
- c(a) crosses at a\* = √((√2−1)θ/2), with residual 0.
- c\*² − (1+√2)θ = 0.
- a\*² − θ/2 = θ(√2 − 2)/2 < 0.
- On the grid θ = n/1000, n < 1000, the three equivalences (c\* < 1, c\* < (1+θ)/√2, θ < √2 − 1) agree, and
  2a\* < c\* holds throughout. There are 0 failures.

| d | θ_inf = 3/(d+1) | a\* | c\* (RA) | Theorem A inf |
|---|---|---|---|---|
| 7 | 3/8 | 0.27868 | 0.951488 | 0.972272 |
| 9 | 3/10 | 0.24926 | **0.851037** | 0.919239 |
| 11 | 1/4 | 0.22754 | 0.776887 | 0.883883 |
| 13 | 3/14 | 0.21067 | 0.719258 | 0.858630 |
| 19 | 3/20 | 0.17626 | 0.601774 | 0.813173 |
| 29 | 1/10 | 0.14391 | 0.491346 | 0.777817 |
| 49 | 3/50 | 0.11147 | 0.380595 | 0.749533 |
| 99 | 3/100 | 0.07882 | 0.269122 | 0.728320 |

Further checks:
- **d = 9, k = 13** (θ = 16/39, μ = 17/156): c\* = 0.995212 and Theorem A gives 0.997202. At k = 12, θ = 23/54, so c\* > 1.
- **Minimal k with c\* < 1:** d = 7: 21; 9: 13; 11: 13; 13: 14; 19: 17. These are the same as Theorem A's thresholds, as claimed.
- **Fixed-depth iteration at θ = 0.3,** by direct minimisation of max(2a, a + c/√2, a + θ/(2a)) rather than the closed form:
  1, 0.91924, 0.88077, 0.86365, 0.85631, 0.85323, … → 0.851037, with g′(c\*) = √2 − 1.
- **P sequence** d_j = 2j+7, k_j = ⌈d_j²/4⌉ + j, η₀ = 1/d²: admissible for all j ≤ 60. j = 7: θ = 0.1465, c\* = 0.5948;
  j = 60 (d = 127): c\* = 0.2394.
