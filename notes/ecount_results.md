# ecount: #E(Z/N) for one fixed curve, Dieulefait–Urroz, and partial counts (phase 23)

Agent `ecount`, 2026-10-22. Version **v4**. It follows the codex adversarial round (MAJOR REVISION), the Opus literature
critic (MAJOR PROBLEMS: prior art) and the Opus proof critic (ACCEPT WITH FIXES). Every required fix is applied (§8,
§9).

- Pre-registration: `notes/ecount_prereg.md`, commit `d856b0b`, made before any experiment ran.
- Code:
  - OPERATIONAL: `src/ecount/zn_curve.py`. It never sees p or q, never calls integer factorisation and never imports the
    harness. An AST test checks this, and codex inspected it by hand.
  - HARNESS: `src/ecount/harness.py` (PARI via cypari2 2.2.4, added to `.venv`) and `src/ecount/run_experiments.py`.
- Tests: `tests/test_ecount.py` (32 pass).
- Raw data: `experiments/raw/ecount_{e1,e2,e3,e4,posthoc,posthoc3,posthoc_xe}.jsonl` and `ecount_e1_summary.json`.
- GPT-6.1-sol:
  - blind round `codex/ecount_1.{md,out}`, cited as **cx1** (scratch scripts in `codex/scratch/ecount_1/`);
  - adversarial round `codex/ecount_2.{md,out}`, cited as **cx2**. The .out is copied from codex's scratch dir, because
    its sandbox could not write the repo copy.
- Labels:
  - **KNOWN** (with reference);
  - **PROVED** (proof here or in cx1/cx2, re-checked here);
  - **CONDITIONAL**, **HEURISTIC**, **OPEN**;
  - **experimentally supported**;
  - **REFUTED (prediction)**;
  - **UNVERIFIED**.
  The strongest novelty wording is "no prior statement found by us".

Notation:
- N = pq with p < q primes and gcd(N, 6Δ_E) = 1.
- E: y² = f(x) = x³ + ax + b is a FIXED curve over Z. E′_r denotes the quadratic twist of E/F_r by a non-residue mod r
  (all such twists are F_r-isomorphic). We write E′ when r is clear.
- A_r = #E(F_r) = r + 1 − a_r and A′_r = #E′(F_r) = r + 1 + a_r. λ_r and μ_r are the exponents of E(F_r) and E′(F_r).
- M = #E(Z/N) := A_pA_q, the projective count (CRT: E(Z/N) = E(F_p) × E(F_q)). M* = (A_p − 1)(A_q − 1) is the
  affine count.
- r_r = #{roots of f mod r}.
- Exc(E) = {N : μ_p | M and μ_q | M} for M = #E(Z/N). For a general multiple M, write Exc_M for the same condition
  with that M.
- X_E = Exc(E) ∩ {M odd} (cx1).

---

## 0. Verdict (plain language)

**The selfred open target ("can #E(Z/N) for one fixed ordinary curve be computed without factoring?") was in substance
already answered in the literature.** The project marked it OPEN because of a literature miss.

**What is known.**
- One count M of one given curve lets you factor N, by a randomized twist trick.
  - Pick x₀ at random. It is the abscissa of a point on E at one prime and on the twist E′ at the other.
  - Multiplying by M kills only the E-component, and a gcd exposes the factor.
- Sources:
  - Okamoto–Uchiyama (EUROCRYPT '98, §4 and Thm 9, under a large-prime-factor hypothesis).
  - Galbraith–McKee 2005 §2: "a single elliptic curve over Z/NZ with a known number of points M", with both the twist
    method and the 2-adic (Miller) method, presented heuristically.
  - Galbraith's book, Example 24.2.30.
  - Dąbrowski–Pomykała–Shparlinski Algorithm 3.1: exactly our point (δx₀, δ²), for a varying curve family.
  - Kunihiro–Koyama 1998 and Martín–Morillo–Villar 2001, for ℓ-adic and 2-adic descent.
- We rediscovered the twist method (Theorem TW). cx1 rediscovered it with 2-adic stripping (TW*).

**What this round adds.** No prior statement of these was found by us; priority is UNVERIFIED.
1. A rigorous, uniform statement for ANY fixed E/Q and EVERY N = pq of good reduction outside an explicit set X_E: the
   per-trial success probability is ≥ 1/8 − O(p^{−1/2}).
2. X_E has density zero among balanced prime pairs (Serre open image plus Chebotarev; for CM, via split primes; cx1).
3. **X_E = ∅ for every curve with a rational 2-torsion point.** For such a fixed curve, computing #E(Z/N) is
   randomized-poly-time equivalent to factoring on ALL good N.
4. Given the factorisation of M, a randomized reduction works for every N (cx1; repairs D–U).

**Still OPEN.**
- An everywhere reduction for an arbitrary fixed curve on X_E. X_E is nonempty for some fixed non-CM curves: cx1's
  contrived example, and one census member for 11a1.
- A deterministic one-curve reduction.

**Consequence for the project.**
- A poly-time algorithm computing #E(Z/N) for one fixed curve would factor every N outside the density-zero set X_E
  (all N, for curves with rational 2-torsion).
- So the selfred §6(a) "moved-divisor" prototype is not a route below factoring, except conceivably on X_E.

**Dieulefait–Urroz, verified in both versions** (arXiv v1 2019; *Mathematics* 8 (2020) 2126; §1):
- Their single-curve malleability result needs the **factorisation** of the count.
- It is a random-curve statement, and its divisor lemma (Prop 1) has an invalid step.
- Their projective two-count Theorem 1 rests on a circular lemma, which uses 4(p+1)(q+1) and so already factors N.
- Their ERH theorem queries twist counts: it is a multi-curve oracle.
- Their "small difference" circumradius formula is wrong, and its short-arc claim has a concrete counterexample (cx2,
  verified here).
- They did not cite Okamoto–Uchiyama, whose Reduction 1 already treats a single given curve.

**Equivalences (§3).**
- M ⇒ factoring off X_E.
- (a_N, a_{N²}) ⇒ factoring (Bach–Charles).
- **(a_p, a_q) alone, or a_N alone: OPEN.** Trace-only data at many ℓ give exactly (a_p, a_q), so FROB's "numbers, not
  operators" scope gap stays open for trace outputs.

**Partial information (§4).**
- M mod ℓ with ℓ an input ⇒ factoring off X_E, for balanced N; for general N, take W = O(N) and B ≈ (ln 2 + o(1))n.
  This is the fixed-curve endgame of Kunihiro–Koyama Lemma 2 / Thm 2.
- Fixed-ℓ counts are OPEN in general, though some are trivially constant: rational ℓ-torsion gives M ≡ 0.
- Their values leak fixed-base residuosity bits. Chebotarev advantages:
  - M mod ℓ: 9/128, 73/4608 and 241/41472 for ℓ = 3, 5, 7;
  - M mod 2: 2/9.
- **The root count of f mod N decides (Δ_f/p) exactly on (Δ_f/N) = 1** (cx2). So computing it is at least as hard as
  fixed-base residuosity for base Δ_f.
- For y² = x³ + 6x − 7 the root count is *equivalent* to the −3-residuosity bit (cx1).

---

## 1. Dieulefait–Urroz (verified from primary text)

### 1.1 Sources
- **arXiv:1911.11004 v1** (25 Nov 2019; the only arXiv version). We read the LaTeX source in full.
- **Published version.** L. V. Dieulefait, J. Urroz, "Factorization and Malleability of RSA Moduli, and Counting Points
  on Elliptic Curves Modulo N", *Mathematics* 8(12) (2020) 2126, doi:10.3390/math8122126. Received 29 Sep 2020,
  accepted 18 Nov 2020, published 27 Nov 2020.
  - We read the full text extracted from the Universitat de Barcelona repository PDF; cx1 located it.
  - The literature critic read it independently.
- **Numbering.** Published (arXiv): Thm 1 (1), Rem 1 (2), Lem 1 (3), Thm 2 (4), Rem 2 (–), Thm 3 (5), Thm 4 (6),
  Prop 1 (–), Thm 5 Coppersmith (7), Cor 1 (–), Rem 3–4 (8–9).
- Published Thm 1 is **projective-only**, since |E(Z/NZ)| is defined as the group. The affine two-count statement is
  arXiv Thm 1 / published Thm 2.
- D–U's bibliography cites Kunihiro–Koyama as "LNCS 1043". The correct citation is EUROCRYPT '98, LNCS 1403, 47–58,
  doi 10.1007/BFb0054116.

### 1.2 What each statement says, and its status

| Published (arXiv) | Statement | Status here |
|---|---|---|
| Thm 2 (4), affine | E′(Z/N) = (p − a_p)(q − a_q) and E′_d(Z/N), (d/p) = −1 ⇒ factor N | **Correct for (a_p, a_q) ≠ (0, 0).** Published Remark 2 itself excludes a_p = a_q = 0, where every affine twist count equals N (E3b: 60/60). arXiv omits that qualification. If exactly one trace vanishes, gcd(E′, N) already splits N. E3a: 180/180 + 180/180 |
| Lem 1 (3), projective | "Knowing two among Ẽ, Ê, Ē and E, we know the four" | **Proof circular in both cases.** Case 1 uses Ẽ + Ē = 4PQ − L, and Case 2 uses (M+1)Ê = 4PQ − L, with PQ = (p+1)(q+1) = N + p + q + 1. This is a gap in the proof, not a counterexample to the conclusion. Found by cx1, cx2, the literature critic and us |
| Thm 1 (1), projective | \|E(Z/N)\|, \|E_d(Z/N)\| and N pairwise distinct ⇒ factor N deterministically | **Unproved as published.** Repaired in §1.3 for (d/N) = −1, under balance. For (d/p) = (d/q) = −1 there is a randomized repair by TW* on N ∉ X_E. A deterministic repair of that case is OPEN |
| Thm 3 (5), ERH | Factoring ≡ "counting points, affine or projective, of elliptic curves E mod N" (published adds "non-supersingular for both primes") | The proof "appl[ies] Theorem 2 to the pair E, E_d for every d up to this bound", so it uses **twist counts** for d ≤ O(log²N): a multi-curve oracle. Projective inherits Lem 1's gap. The repair needs the least d with (d/N) ≠ 1, which is < 2 ln²N under ERH (Bach, Math. Comp. 55 (1990) 355–380); if gcd(d, N) > 1 that also factors N |
| Thm 4 (6) + Prop 1 | Random curve E mod N (isomorphism-class model); the **factorisation** of the projective or affine count gives the factors of N "with probability one" | **Core correct.** The divisor u = q + 1 − a_q of the count has \|u − q\| ≤ 2√q + 1, then Coppersmith. Deterministic when d(count) ≤ poly(n) (E3c 180/180). Three caveats. (i) Prop 1's step "E ∈ S_N ⇒ E ∈ S_p or E ∈ S_q" is invalid (cx1). Raising the global threshold to (log N)⁸ repairs that inclusion only. The displayed local estimates give O(1/log p), not the claimed o(1/log N), and nothing gives "probability one" (cx2). Local exponent 5 with global exponent 10 would supply the missing logarithm. (ii) It covers random curves, not a fixed curve against adversarial N. (iii) arXiv v1's justification was only the average order of d(n). E3c measured d(M) up to 16896 at n = 96 |
| "Small difference" section | Circumradius ⇒ ≤ 2 lattice points on arcs of xy = T of length < (N/32)^{1/3}, if \|p − q\| < c′N^{1/3} | **Wrong formula; false assertion.** The correct radius is R = √(∏_{i<j}((x_ix_j)² + T²))/(2Tx₀x₁x₂); the paper prints R² (cx1). Near the diagonal R ≍ √T. **Counterexample (cx2; re-verified exactly here):** p = 1080557, q = 1080559, E: y² = x³ + 583803335960x + 583803335961, supersingular at both (≡ y² = x³ + 1 mod p, p ≡ 2 (3); ≡ y² = x³ − x mod q, q ≡ 3 (4)). Then T = (p+1)(q+1) = 1167607752480, and xy = T has three lattice points (1078482, 1082640), (1079520, 1081599), (1080558, 1080560) on an arc of length < 2942 < (N/32)^{1/3} ≈ 3316.8 |
| Intro / Cor 1 | "not known how to factor N only with the number E_N as input" (arXiv); "no known polynomial time algorithm that can factor N using this information" (published) | **Inaccurate when written.** Okamoto–Uchiyama 1998 (Thm 9: given curve, under a large-prime hypothesis) and Galbraith–McKee 2005 (§2, heuristic) already address it. TW* gives the rigorous fixed-curve version off X_E. The affine one-curve count remains OPEN |

### 1.3 Repaired projective two-count reduction (PROVED, routine; balanced case)
- Let (d/N) = −1, say (d/p) = −1 and (d/q) = +1. Then E = (P − a)A_q and E_d = (P + a)A_q, with P = p + 1 and
  a = a_p ≠ 0 (published Thm 1's distinctness excludes a = 0).
- Put S = E + E_d = 2P·A_q and D = E_d − E = 2a·A_q. Then h := S/gcd(S, |D|) = P/k with k = gcd(P, a) ≤ 2√p, so
  p = kh − 1.
- **Balance assumption (cx2).** p ≥ N^{1/2 − O(1/log N)}, e.g. q/p bounded.
  - Then f(y) = y − h^{−1} mod N has the root k = O(N^{1/4}) modulo a factor of size N^{1/2 − o(1)}.
  - A univariate small-root theorem with a dimension schedule growing with n recovers k deterministically (May 2003,
    §3.2 Thm 7). Here gcd(h, N) = 1 always, because h | p + 1 and h < q.
  - This also works when the flipped prime is q, so the algorithm need not know which prime is flipped.
  - For an unbalanced smaller flipped prime r = N^β with β < ½, k = O(N^{β/2}) exceeds the N^{β²} range, so that case
    is not covered.
- **Implementation is a toy.**
  - `du_projective_twist_repaired` scans k ≤ 4096 by gcd, then falls back to a fixed-dimension (10) lattice.
  - That lattice has radius 2^{(2/9)n}, so covering the window needs ~2^{n/36} calls (cx2). The same holds for
    `du_malleability`.
  - E3d: 180/180, all via the gcd scan. The lattice branch passes a forced test (n = 64, 96; 6/6).
  - None of this is asymptotic evidence.
- **Case (d/p) = (d/q) = −1.** There are two equations in (s, a_p, a_q). TW* covers it off X_E.
- **cx1's randomized four-twist variant.** The counts of E, E_d, E_e and E_de, with independent sign vectors (probability
  3/8 per attempt), sum to 4(p+1)(q+1).

---

## 2. Theorems TW and TW*: one fixed curve, count ⇒ factoring

### 2.1 Lemma ZN: the group law over Z/N (PROVED; cx2: CORRECT, 10985 independent checks)
- Elements are O (both components O) or (x, y) affine in both components.
- `ec_add` proceeds as follows.
  - If gcd(x₂ − x₁, N) = 1, use the chord formula.
  - If that gcd is proper, output it.
  - If it is N, let h = gcd(y₁ + y₂, N). Then h = N gives O, h = 1 gives the tangent formula, and any other h is
    output.
- Claim: the result is the CRT-correct sum, or a proper factor.
- *Proof.* Work componentwise. With x₁ = x₂ at both primes, y₂ = ±y₁ at each, and y₁ + y₂ = 0 covers inverse pairs
  including y = 0. Mixed configurations make one gcd proper. ∎
- **Corollary.** `ec_mul(k, P)` (left-to-right double-and-add) evaluates [k′]P for every binary prefix k′ of k. If any
  of these is O in exactly one component, it outputs a proper factor.

### 2.2 Theorem TW (twist trick with Jacobi filter)
This is a rediscovery of Okamoto–Uchiyama's Reduction 1 mechanism, here without their hypotheses.

Let M be any integer with λ_p | M and λ_q | M.

**One trial.**
1. Pick x₀ uniform in Z/N and set δ = f(x₀). If gcd(δ, N) is proper, output it.
2. If (δ/N) = −1, compute [M]P for P = (δx₀, δ²) on E^(δ): Y² = X³ + aδ²X + bδ³.

**Claim.** One trial splits with probability ≥ ρ(p, q) (for the given M; the exceptional set below is Exc_M), where

  ρ(p, q) = [(A_p − 1 − r_p)(A′_q − |H′_q|) + (A′_p − |H′_p|)(A_q − 1 − r_q)] / (4pq), with H′_r = E′(F_r)[M].

Outside Exc(E), ρ ≥ 1/8 − O(p^{−1/2}). Bounded small p are handled by trial division.

*Proof (cx2: CORRECT).*
1. Write x₀ = (u, v) by CRT.
2. If f(u) is a nonzero square mod p, then P_p lies on a curve isomorphic to E/F_p, so M kills it.
3. Mod q, write δ = c·s² with c a fixed non-square. The map (X, Y) ↦ (X/δ, sY/δ²) is an F_q-isomorphism onto
   c·y² = f(x). It sends P_q to (v, s) and preserves M-torsion.
4. Non-square v correspond two-to-one to points of E′(F_q) with y ≠ 0. Points with y = 0 lie in H′_q, because then
   2 | λ_q | M. So there are (A′_q − |H′_q|)/2 good v.
5. For such (u, v), [M]P is O only at p, so Lemma ZN outputs p. The mirror case is symmetric.
6. If μ_q ∤ M, then |H′_q| ≤ A′_q/2. ∎

### 2.3 Theorem TW* (no filter, 2-adic descent; cx1's argument, transcription checked by cx2)
TW* is TW without the Jacobi filter. At each prime, P lies on E or on E′ independently.

**Claim.** If N ∉ X_E, one trial splits with probability ≥ 1/8 − O(p^{−1/2}).

The proof critic sharpened this. Off Exc, TW* has α_r(M) = (A_r + L_r − 2)/(2r) ∈ [½ − o(1), ¾ + o(1)] at a prime with
μ_r ∤ M, so the true bound is ≥ ¼ − O(p^{−1/2}) in both branches. That gives an expected ≤ 4 + o(1) trials.

*Proof.*
- **Twist part.** As in TW.
- **2-adic part.** Suppose N ∈ Exc and M = 2^s m with s ≥ 1.
  - Take r with A_r even, and j ≥ 1 minimal with 2^j m annihilating both E(F_r) and E′(F_r).
  - At k = 2^{j−1}m both kernels have index ≤ 4, and one is proper. So α_r(k) ∈ [¼ − O(r^{−1/2}), ¾ + O(r^{−1/2})].
  - The split probability is α_p(1 − α_q) + (1 − α_p)α_q ≥ min(α_r, 1 − α_r).
  - k is a prefix of M, so the chain visits it. ∎
- **Exact non-root kernel mass (cx2).** It is K_r(k) + L_r(k) − 2 − 2r_r·[2 | k], divided by 2r, or by 2(r − r_r)
  after conditioning on f(x₀) ≠ 0. This is an O(1/r) correction.

**Corollary 2T (PROVED, cx1; cx2 CORRECT).** If E(Q) has a point of order 2, every A_r is even, so X_E = ∅. For such
a fixed curve, #E(Z/N) ≡ factoring (randomized) on all N of good reduction. Example: C3: y² = x³ + 2x + 3 =
(x + 1)(x² − x + 3).

### 2.4 Corollary EQ, density, exceptional sets

**Corollary EQ (PROVED).** Fix E (CM or non-CM). On the promise N ∉ X_E, #E(Z/N) and factoring are equivalent under
randomized polynomial-time reductions, via **TW*** (TW alone covers only N ∉ Exc):
- one count and expected ≤ 4 + o(1) trials of Õ(n³) bit cost in one direction;
- Schoof/SEA in the other.

**Doubly supersingular N.** Here M = (p+1)(q+1), and SR2 factors N deterministically. Such N lie in Exc, and are also
covered by TW*, since M is even.

**Density zero of X_E (PROVED; cx1's argument with the CM restriction restored per cx2).** The pairs (p, q) are iid
uniform primes in [X, 2X], and density is taken as X → ∞.
- **Non-CM curves.**
  - Take a finite set S of large ℓ with full image in ∏_{ℓ∈S} GL₂(F_ℓ). This uses the adelic open image (Serre,
    Invent. Math. 15 (1972)) plus Goursat and a surjective determinant; per-ℓ surjectivity alone is not enough.
  - Consider the event "ρ_ℓ(Frob_p) has eigenvalues −1 and t ≠ ±1, and ρ_ℓ(Frob_q) has no eigenvalue 1". Its
    probability is ((ℓ−3)/(ℓ−1)²)·(1 − O(1/ℓ)).
  - The event forces ℓ | A′_p and ℓ ∤ M, so μ_p ∤ M.
  - Chebotarev (Serre, Publ. IHÉS 54 (1981)) bounds the density of Exc ⊇ X_E by ∏_{ℓ∈S}(1 − c/ℓ), which tends to 0.
- **CM curves.** Exc does NOT have density zero: pairs of inert primes have a_r = 0 and lie in Exc, with density 1/4
  (E2).
  - X_E does have density zero. M odd excludes inert primes, since there A_r = r + 1 is even.
  - So restrict to split p and q, and run the same argument with split auxiliary ℓ, where the image is the Cartan
    subgroup (F_ℓ^×)² (Serre 1972 §4.5).
  - Details added at the proof critic's request:
    - Pairs involving an inert prime (inert–inert or inert–split) have M even, so they are not in X_E. Ramified primes
      are finitely many.
    - Chebotarev is applied to degree-one primes of K.
    - The event is invariant under swapping the two Cartan coordinates (complex conjugation).
    - The open image over K gives the full product ∏(F_ℓ^×)² for large split ℓ.
    - So X_E has density zero among all balanced pairs.
- No rate is claimed.

**Census (E2, PH4).**
- Non-CM Exc fractions at b = 18: 3.9e−6, 7.8e−7 and 1.9e−6.
- X_E has exactly **one** member among ≈ 7.6·10⁷ pairs: 11a1 with (p, q) = (54443, 62819) and M odd. Both twist
  groups are cyclic, and μ_p = A′_p = 3·13·19·73 and μ_q = A′_q = 3²·5²·281 both divide M = 3³·5²·13·19·73·281.
- cx1's X_E example y² = x³ + 27331x + 2 at N = 27331 is contrived: the curve reduces to a CM curve mod N.

**X_E is not this code's failure set (cx2).**
- Intermediate non-divisor prefixes can split. TW* split the census X_E member, and cx2 found 4797 TW and 9153 TW*
  structural splits among all 27331 x₀ for cx1's example.
- The registered claim "TW fails on every trial iff N ∈ Exc(E)" (prereg) is **REFUTED**. Exc(E) is where the *proof*
  fails, not where the algorithm fails.

**Factored count ⇒ factoring for every N (PROVED, cx1; cx2 CORRECT).**
- Two odd local exponents cannot coincide. If both equal an odd L, then L | r + 1. Writing the groups as C_d × C_L and
  C_e × C_L with d, e | gcd(L, r − 1) forces d = e = 1 and L = r + 1, which is even.
- So stripping the primes of the supplied factorisation one at a time separates the kernels at some stage. This needs
  scalars M/ℓ^j, which we did not implement.

**Deterministic version: OPEN.**

### 2.5 Relation to the literature
Sources: the literature critic and cx2. OU, GM and DPS were read in primary text by the critics.
- **Okamoto–Uchiyama**, "Security of an identity-based cryptosystem and the related reductions", EUROCRYPT '98, LNCS
  1403, 546–560.
  - §4: anomalous E_n; a random x-coordinate is multiplied by n with x-only arithmetic. The mixed E/twist cases give p
    or q with probability ½.
  - Thm 9 (Reduction 1): a given curve, with "order" λ = lcm(#E_p, #E_q), and random x with (f(x)/n) = −1. It assumes
    #E_p and #E_q have large prime factors, and that p′ ∤ q + 1 + s "with overwhelming probability" (average-case).
  - Thm 10 covers random curves.
  - **TW is OU's Reduction 1 with λ replaced by any common multiple M, and with their hypotheses replaced by Exc(E).** An
    exact theorem-level comparison is UNVERIFIED, because the full proof was only partly accessible to cx2.
- **Galbraith–McKee**, "Pairings on elliptic curves over finite commutative rings", IMA Cryptography and Coding 2005,
  LNCS 3796, 392–409, §2.
  - Setting: "a single elliptic curve over Z/NZ with a known number of points, M".
  - Method (1): multiply a point on a random quadratic twist by M using x-coordinates (citing OU).
  - Method (2): the Miller-style odd part followed by doublings (citing MMV).
  - Presented heuristically ("for most curves we would not expect …"). **TW + TW* as a known heuristic algorithm.**
- **S. D. Galbraith**, *Mathematics of Public Key Cryptography* (CUP 2012), §24.2.4, Example 24.2.30: the fixed-curve
  twist example, with an unqualified success assertion (cx2).
- **Dąbrowski–Pomykała–Shparlinski**, "On oracle factoring of integers", J. Complexity 76 (2023) 101741,
  arXiv:1912.00345.
  - Algorithm 3.1: P = (αx, α²y) with α = f_b(x)/y². At y = 1 this is our point.
  - Thm 1.1: deterministic, for the family y² = x(x² − b) with b ≤ (log N)^γ, for all n ≤ N with few exceptions.
- **Kunihiro–Koyama**, EUROCRYPT '98, LNCS 1403, 47–58.
  - Thm 1: FCT ≡_RP counting #E_n(a, b) for arbitrary given (a, b). Random curves through a chosen point, S-adic
    stripping, and a distributional heuristic.
  - Remark (3): a one-curve variant.
  - Lemma 1: a random multiple of #E_n (the Miller analogue).
  - **Lemma 2 + Thm 2: #E_n mod small d, then CRT, then factoring. This is our Prop CRT-ℓ for arbitrary curves.**
- **Martín–Morillo–Villar**, Appl. Math. Lett. 14(3) (2001) 341–346: a point-order oracle; GM credit it for the 2-adic
  method. Full text paywalled; content UNVERIFIED.
- Further varying-curve work:
  - Dryło–Pomykała, Colloq. Math. 159 (2020) 259–284 (point-exponent oracle, error 53/80);
  - Jiménez Urroz–Pomykała, Ramanujan J. 64 (2024) 265–273, Thm 1 (GRH, random curve plus twist; verified by the
    literature critic).
- **No prior statement found by us** of:
  - the uniform fixed-E/Q, all-N bound outside the explicit X_E;
  - the density-zero proof;
  - the rational-2-torsion corollary (everywhere equivalence);
  - the randomized factored-count repair for every N.
- The "Miller analogue" framing (M/N = ∏ L_r(E, 1)^{−1}, §3) is presentation only. The analogy itself is KK Lemma 1 /
  GM §2.

---

## 3. Can #E(Z/N) be computed without factoring? Complexity and equivalences

**Answer.** Not on N ∉ X_E (all N for rational-2-torsion curves) unless factoring is easy there (Cor. EQ). On X_E it is
OPEN. The standard algorithm is to factor N (NFS, heuristic L_N[1/3, (64/9)^{1/3}]) and then run SEA on p and q. This
is not a proved optimality statement.

**Schoof over Z/N.**
- Producing Frobenius on (Z/N)[x]/ψ_ℓ factors N when (N/ℓ) = −1 (nonsmooth FROB(d); cx1).
- **Correction to v2 (cx2).** The claim "trace-only Schoof data is factoring-hard" is withdrawn.
  - Our argument used p mod ℓ, which is factor data: enough residues p mod ℓ give p directly.
  - Trace-only data (a_p mod ℓ, a_q mod ℓ) for many ℓ amount to the integers (a_p, a_q), whose factoring power is OPEN
    (below).
  - So FROB's "numbers, not operators" scope gap remains open for trace outputs.

**Hecke operators and modular symbols.**
- a_N = a_pa_q, whereas M = N + p + q + 1 − (p+1)a_q − (q+1)a_p + a_pa_q.
- T_N sums over σ₁(N) matrices (Stein, *Modular Forms: A Computational Approach*, GSM 79, Ch. 2, level 1).
- Modular-symbol Hecke algorithms (Stein Ch. 3; Ch. 8, §8.3) are not polynomial in log N.
- a_N = Σ_{x mod N}(f(x)/N) (cx1) has N easy terms, so evaluating it is exponential.

**Euler product (PROVED, one line).** M/N = ∏_{r|N} A_r/r = ∏_{r|N} L_r(E, 1)^{−1}.

**Bach–Charles** (arXiv:0708.1192 v2, Thm 2.2, Cor 2.7–2.8).
- A poly-time algorithm for a(n) at all n, for a fixed non-CM even-weight eigenform, factors a density-1 set of RSA
  moduli. The proof queries a(N) and a(N²).
- Weight-2 explicit form (cx1): S = a_N² + N − a_{N²} = qa_p² + pa_q², and both terms are roots of Z² − SZ + N·a_N².
- It fails iff a_p = a_q = 0. Then M = (p+1)(q+1) still factors N.

**Oracle relations given N** (R = randomized, D = deterministic).
- M ≡_R FACT, off X_E.
- FACT ⇒_D (a_p, a_q) ⇒_D a_N.
- (a_N, a_{N²}) ⇒_D FACT, unless a_p = a_q = 0.
- **(a_p, a_q) ⇒ FACT, and a_N ⇒ FACT: OPEN.** We found no computable multiple of A_p from (N, a_p, a_q). The formulas
  that need p are not impossibility proofs.
- #E(Z/N^k) = N^{k−1}M (cx1).

**Lattice view.**
- M = (p + B)(q + A) with |A|, |B| = O(N^{1/4}) under constant-ratio balance.
- Implicit factoring with shared MSBs needs t ≥ 2α + 3 bits (Faugère–Marinier–Renault, PKC 2010, LNCS 6056, 70–87).
  Here t ≈ n/4 and α = n/2.
- The general three-variable problem includes M = N, i.e. plain factoring.
- TW* does not make lattices moot on X_E. Lattices are one possible route for M* and X_E, not established as the only
  one.

**Edixhoven–Couveignes.**
- Thm 15.2.1 is at level one (T(1, k)) and takes n *with its factorisation*.
- At general level (Bruin; Javanpeykar, ANT 8 (2014), Thm 5.0.1) the factored index is also needed.
- An adaptation that output M at composite N off X_E would be a factoring algorithm.

---

## 4. Partial information

### 4.1 Prop CRT-ℓ (PROVED; the fixed-curve endgame of KK Lemma 2 / Thm 2)
- For balanced N (q/p < 4/3): √p + √q ≤ ((4/3)^{1/4} + (3/4)^{1/4})N^{1/4} < 2.01N^{1/4}. Hence |M − N| ≤ W with
  W ≤ 4.1N^{3/4} + O(√N).
- Primes ℓ with ∏ℓ > 2W + 1 determine M by CRT. The cutoff is B = (0.52 + o(1))n, i.e. asymptotically π(B) calls; in
  practice, choose B by the product condition.
- Then TW*. So "M mod ℓ with ℓ an input" ⇒ FACT holds on N ∉ X_E, and on all N for rational-2-torsion curves.

### 4.2 Fixed ℓ

**M mod ℓ.** M mod ℓ = det(1 − g_p)·det(1 − g_q), with g_r = ρ_{E,ℓ}(Frob_r).
- **Predicting β_ℓ = (p/ℓ) on (N/ℓ) = 1 from (N mod ℓ, M mod ℓ).** Under independent Haar classes in GL₂(F_ℓ) (full
  image), the Bayes advantage is:

  | ℓ | 3 | 5 | 7 |
  |---|---|---|---|
  | advantage | 9/128 = 0.0703 | 73/4608 = 0.0158 | 241/41472 = 0.0058 |

  - These are exact rationals, cross-checked by cx2.
  - They are model predictions, not finite-prime theorems.
- Some cases are easy. Rational ℓ-torsion gives M ≡ 0 (C1 at ℓ = 5; PH2 increment exactly 0), and C3's M is always
  even.
- v2's "strictly weaker than φ(N) mod ℓ" is withdrawn (cx2). M mod ℓ involves trace data, and no deterministic
  information ordering is established.

**ℓ = 2: three different observables (cx2 correction of v2's false "same information").**
- Setting: f irreducible with Gal(f) = S₃, conditioned on (Δ_f/N) = 1, and β = (Δ_f/p).
- S₃-Chebotarev joint masses:

  | observable | value | β = +1 | β = −1 |
  |---|---|---|---|
  | M mod 2 | odd | 2/9 | 0 |
  | M mod 2 | even | 5/18 | 1/2 |
  | r_p·r_q | 0 | 4/9 | 0 |
  | r_p·r_q | 9 | 1/18 | 0 |
  | r_p·r_q | 1 | 0 | 1/2 |
  | #E(Z/N)[2] | 1 | 2/9 | 0 |
  | #E(Z/N)[2] | 16 | 1/18 | 0 |
  | #E(Z/N)[2] | 4 | 2/9 | 1/2 |

- **M mod 2: advantage 2/9.**
  - "M odd ⇒ β = +1" is a one-sided certificate (0 exceptions in E4).
  - In a finite population the advantage is Pr[M odd] + max(Pr[β = −1], Pr[M even, β = +1]) − ½.
- **Root count r_N(f) = r_p·r_q: advantage ½, i.e. it decides β exactly (PROVED, cx2). This holds for any separable
  cubic, not only S₃, by Stickelberger: (Δ_f/r) = −1 iff r_r = 1.**
  - On (Δ_f/N) = 1 the two Frobenius signs agree.
  - A product of 1 means two transpositions, so β = −1. Products 0 or 9 mean even permutations, so β = +1.
  - So **counting roots of the fixed cubic f mod N decides fixed-base quadratic residuosity of Δ_f mod N**, given
    (Δ_f/N) = 1.
  - Tested in `test_cubic_root_count_decides_disc_bit`.
- **#E(Z/N)[2]: advantage 5/18.**
- The pre-registration's "equivalently" (E4 (ii)) was false. E4 measured only M mod 2.

**Exact equivalence (cx1, PROVED; cx2 CORRECT).**
- E*: y² = x³ + 6x − 7 = (x − 1)(x² + x + 7), with j = 2048/3 (non-CM), has r_r = 2 + (−3/r).
- So R_N ∈ {9, 3, 1} is equivalent, given N, to "is −3 a square mod N" on N ≡ 1 (mod 3). It is also equivalent to the
  information in φ(N) mod 3 and in |E*[2](Z/N)| ∈ {16, 8, 4}.
- Not claimed factoring-equivalent.

**Status.**
- M mod ℓ for fixed ℓ, M mod 2 and r_N(f) are OPEN, except the trivially constant cases.
- The root count is at least as hard as the fixed-base Δ_f residuosity bit. That bit is itself OPEN, like φ(N) mod ℓ
  (selfred §2.1).

**Consistency with earlier results.**
- **Theorem PT.** A pure tree with advantage δ must split with probability ≥ δ/2. TW-type algorithms use Jacobi
  symbols, uniform sampling and M-dependent schedules, all outside PT's model. No conflict.
- **RANK.** Linear divisor formulas give s mod ℓ. The elliptic residues carry trace data and are not ordered relative to
  s mod ℓ. Prop RANK's list of easy cases extends: rational ℓ-torsion gives a constant residue.

---

## 5. Experiments
Pre-registered: E1–E4. Post-hoc, labelled: PH1–PH4. Peak RSS was ≤ 450 MB. Heavy runs ran alone; the light post-hoc
runs (< 5 s) overlapped E3.

### E1. TW (Jacobi filter): 40 instances per (curve, n), 64 trials each, inputs (N, a, b, M) only

| curve | solved, n = 32, 64, 96, 128 | per-trial rate: empirical / exact ρ (z against ρ + bonus) | first split: mean / max trials |
|---|---|---|---|
| C1 = 11a1 | 40, 40, 40, 40 | .496/.498 (−0.2), .512/.500 (1.2), .498/.500 (−0.2), .501/.500 (0.1) | 1.9 / 7 |
| C2 = 37a1 | 40, 40, 40, 40 | .495/.493 (0.2), .481/.500 (−1.9), .514/.500 (1.4), .478/.500 (−2.3) | 1.8 / 7 |
| C3: y² = x³ + 2x + 3 | 40, 40, 40, 40 | .498/.493 (0.5), .505/.500 (0.5), .494/.500 (−0.6), .490/.500 (−1.0) | 2.0 / 7 |
| CM1: y² = x³ − x (10, 17, 10, 10 in Exc) | 40, 40, 40, 40 | .43/.24, .38/.19, .42/.24, .44/.24 (z = 21–24) | 2.4 / 16 |
| CM2: y² = x³ + 1 (10, 11, 9, 11 in Exc) | 40, 40, 40, 40 | .44/.27, .42/.26, .41/.24, .44/.26 (z = 18–21) | 2.5 / 15 |
| CM1, p ≡ q ≡ 3 (mod 4) | TW 20/20 per n; SR2 20/20 per n | — | — |

- **P1 HOLDS:** 480/480 non-CM instances split; none is in Exc.
- **P2.**
  - The ±3 SE clause holds in 12/12 non-CM cells.
  - The one-sided "≥ ρ" clause fails in **6** non-CM cells (C1 at 32 and 96, C2 at 64 and 128, C3 at 96 and 128) by
    sampling noise; it was badly posed. The SE uses the pooled mean; the exact heterogeneous variance (a sum of
    Bernoulli variances) is no larger.
  - CM cells: **REFUTED as registered** (z = 18–24), caused by prefix 2-adic splitting.
- **P3 (M + 2 null: 0 splits) is REFUTED. The control was mis-designed.** The explanation is corrected per cx2. There
  were 21711 structural splits:
  - 20480 with 4 | M: the chain visits [M/2]P.
  - 1225 with M odd: the last addition is [M+1]P + P, and [M+1]P = P at the annihilated component, which gives a
    component-mixed tangent case.
  - 6 with v₂(M) = 1.
  - The corrected null PH1 gave 0/15360 structural splits.
- **P4 is REFUTED, in the favourable direction.** TW split 80/80 doubly supersingular instances via 2-adic descent, and
  SR2 split 80/80.
- **Prereg claim "TW fails on every trial iff N ∈ Exc(E)": REFUTED** (§2.4, PH3).

### E2. Exact census of Exc(E) and ρ (all unordered pairs from S_b, b = 12–18)

| curve | Exc fraction b = 12 / 14 / 16 / 18 | b = 18: both-ss / one-ss / ordinary | min ρ outside Exc (all b) | median ρ |
|---|---|---|---|---|
| C1 | 9.0e−4 / 5.4e−5 / 1.3e−5 / 3.9e−6 | 45 / 5 / 5 | 0.124 | 0.500 |
| C2 | 0 / 0 / 1.5e−5 / 7.8e−7 | 6 / 2 / 3 | 0.125 (b = 16) | 0.500 |
| C3 | 0 / 4.3e−5 / 2.8e−5 / 1.9e−6 | 10 / 14 / 3 | 0.124 | 0.500 |
| CM1 | 0.271 / 0.247 / 0.250 / 0.251 | 3517878 / 1260 / 0 | 0.193 | 0.250 |
| CM2 | 0.241 / 0.252 / 0.249 / 0.256 | 3584503 / 1888 / 1 | 0.184 | 0.250 |

- **P5: partly REFUTED.** The bound ≤ 10⁻³ holds in 12/12 cells. Monotonicity fails for C2 and C3, whose counts start
  at zero.
- **P6 HOLDS.**
- **P7.** Every both-supersingular pair is in Exc (C(2653, 2) = 3517878 and C(2678, 2) = 3584503 exactly). The
  "rest obeys P5" clause fails at b = 12 and 14: **partly REFUTED**.
- cx2 independently reproduced the b = 12 vectorised ρ and Exc counts for all five curves.

### E3. Dieulefait–Urroz on toy examples (C1–C3, 20 seeds per cell)

| part | cells | result |
|---|---|---|
| (a) affine two-count, (d/N) = −1 and d QNR at both | n = 32, 64, 128 | 180/180 and 180/180 (**P8**) |
| (b) CM1, p ≡ q ≡ 3 (4): affine counts of E and E_d (d ∈ {2, 3, 5, 6, 7}) all equal N | n = 32, 64, 128 | 60/60 (**P9**) |
| (c) factorisation of M → divisors in [√N/4, 4√N] → toy Coppersmith (dim 10) | n = 48, 64, 96 | 180/180 (**P10**). d(M) median 192–3072, max 16896. In-window divisors ≤ 1316. Tried ≤ 644. ≈ 2.9·10⁵ LLL calls, ≈ 1707 s |
| (d) projective repair (the as-written Lem 1/3 cannot run without PQ) | n = 32, 64, 96 | 180/180 (**P11**); all via the k ≤ 4096 gcd scan |

### E4. Partial information: ordered pairs from S₁₆ (1491 primes), conditioned on (N/ℓ) = 1, resp. (Δ_f/N) = 1

| curve | ℓ = 3 (pred 0.070) | ℓ = 5 (pred 0.016) | ℓ = 7 (pred 0.006) | M mod 2 (pred 2/9) |
|---|---|---|---|---|
| C1 | 0.072 | 0.017 (M ≡ 0 mod 5) | 0.014 | 0.213; M odd ⇒ β = +1 always |
| C2 | 0.050 | 0.019 | 0.014 | 0.214; same |
| C3 | 0.074 | 0.035 | 0.014 | 0.013 (M always even; prediction inapplicable) |

- **P12: 8/10 applicable cells within 0.01.** It fails for C2 at ℓ = 3 and C3 at ℓ = 5: **partly REFUTED**. The
  post-hoc bands (PH2) contain both, but that does not retroactively satisfy the registered tolerance.
- **P13: REFUTED as registered.** The 0.017 is the finite-population N mod 5 baseline. The increment from M is exactly 0.
- **P14: REFUTED (prediction error).** 2/9 is correct subject to the majority condition (§4.2).
- E4 did not separately measure the root count or #E[2], although they were registered "equivalently". That
  equivalence was false (§4.2).

### Post-hoc (labelled)
- **PH1.** TW with M′ uniform in [N − 8N^{3/4}, N + 8N^{3/4}] gave 0/15360 structural splits. This is descriptive,
  not a proof that M is the only possible source of splits.
- **PH2.** Haar-GL₂ synthetic populations (1491 primes, 200 replicates); the bands are 5% / median / 95% of the Bayes
  advantage:

  | ℓ | 5% | median | 95% |
  |---|---|---|---|
  | 3 | .044 | .071 | .098 |
  | 5 | .009 | .024 | .054 |
  | 7 | .008 | .025 | .049 |

  - All observed values are inside the bands.
  - Increments over the N mod ℓ baseline: .070 / .048 / .071 at ℓ = 3; 0 / 0 / .0003 at ℓ = 7; 0 for C1 at ℓ = 5.
  - E4 cannot resolve a 0.006 effect with 1491 primes.
- **PH3.** On Exc members at b = 16 (all 14 / 17 / 31 for C1–C3; 300-pair samples for CM1 and CM2), TW split 13/14,
  17/17, 31/31, 300/300 and 300/300. TW* split all of them, including the single M-odd member.
- **PH4.** X_E over 5 curves × 4 sizes (≈ 7.6·10⁷ pairs) has one member, as above. The CM curves have none.

---

## 6. Most decisive open questions
1. **X_E for an arbitrary fixed curve.** Is there an everywhere randomized reduction from one count, or are there curves
   where the count on X_E carries less? X_E is nonempty but has density zero.
2. **The affine one-curve count M*.** It is not a group exponent, and it equals N on doubly supersingular N.
3. **(a_p, a_q) or a_N alone** (Galois-only data). Do they factor N? This is also FROB's open trace-output case.
4. **Fixed-base residuosity bits reached by partial counts.** These include r_N(f) for base Δ_f, the root counts of E*,
   and M mod ℓ. Their status is that of φ(N) mod ℓ: OPEN both ways.
5. A **deterministic** one-curve reduction, and a rate for the density of X_E.

**For the project.**
- selfred §6(a) is settled off X_E, and in substance was known (OU 1998, GM 2005).
- What remains of the moved-divisor question: selfred §6(b)–(c), plus the X_E / M* / (a_p, a_q) residue above.

## 7. Reproduction and deviations
```
.venv/bin/python -m pytest -q tests/test_ecount.py                 # 32 tests, ~10 s
.venv/bin/python -m src.ecount.run_experiments e1                 # ~60 s
.venv/bin/python -m src.ecount.run_experiments e2                 # ~8 s
.venv/bin/python -m src.ecount.run_experiments e4                 # <1 s
.venv/bin/python -m src.ecount.run_experiments e3                 # ~29 min (toy Coppersmith dominates)
.venv/bin/python -m src.ecount.run_experiments posthoc            # PH1, PH2 (~2 s)
.venv/bin/python -m src.ecount.run_experiments posthoc3           # PH3 (~1 s)
.venv/bin/python -m src.ecount.run_experiments posthoc_xe         # PH4 (~2 s)
```

Deviations (disclosed):
- **D1.** PARI's multithreaded `polmodular` deadlocked under cypari2 in the first E1 run. It was killed and rerun with
  `nbthreads = 1` and a 64 MB → 1 GB stack. The output was regenerated from scratch; seeds are deterministic.
- **D2.** The registered M + 2 null was invalid; PH1 replaces it.
- **D3.** P2's "≥", P13 and P14 were mis-formulated. They are reported as REFUTED, not re-registered.
- **D4.** E3d's lattice branch is never reached on the registered instances. It is forced in a unit test.
- **D5.** `du_affine_twist` follows D–U literally and does not first try gcd(E′, N).
- **D6.** TW* (`jacobi_filter=False`) was added after cx1 and is used only in PH3.
- **D7.** PH3 counts only structural splits. An earlier δ-gcd hit (probability ~r_p/p) could mask a later one.
- **D8.** In E3c, d(M) reached 16896, against "of the order of the logarithm" (arXiv), which would be ≈ 66.
- **D9 (cx2).** In v1/v2, `ec_add` did not increment `Counters.invs`, so the `invs` stored in E1 raw data undercount
  inversions. Fixed in the code; no conclusion is affected.
- **D10 (cx2).** `test_du_projective_lemma_needs_PQ` contained a tautological assertion. It now checks the
  reconstructed counts against the true mixed-twist counts.
- **D11.** The Coppersmith routines are toys with fixed lattice dimension (§1.3).

## 8. Codex adversarial round (cx2: MAJOR REVISION), applied in v3

| cx2 item | Response |
|---|---|
| Global "≡ FACT" / "closed" claims exceed the domain N ∉ X_E | Restricted everywhere (§0, §3, §4.1, §6); X_E is now open question 1 |
| Prior art missed: Galbraith Ex. 24.2.30, OU Thm 9, DPS Alg. 3.1 | §0 and §2.5 rewritten. The mechanism is KNOWN and we rediscovered it. Novelty is limited to the rigorous refinements, with priority UNVERIFIED |
| CM density transcription dropped the split-prime restriction; Exc for CM has density ¼ | Fixed (§2.4) |
| Parity, root count and #E[2] are not the same information | Fixed: advantages 2/9, ½ and 5/18; the root count decides β exactly; E4's omission disclosed |
| The PAR identity needs a majority condition | Fixed |
| The projective repair needs balance; a fixed-dimension lattice is not polynomial | Balance added; toy label (D11) |
| D–U Prop 1: the (log N)⁸ change repairs the inclusion only; small-difference counterexample | Both added; the counterexample was re-verified exactly |
| "Trace-only Schoof data is factoring-hard" overclaims | Withdrawn (§3) |
| "Strictly weaker than RANK" is unsupported; there are easy constant cases | Withdrawn; easy cases added |
| P2 fails in 6 cells, not 5; the M + 2 explanation is parity-dependent; the "fails iff" refutation was unreported | All fixed (§5) |
| `Counters.invs`; tautological test; exact non-root kernel mass; CRT-ℓ cutoff wording; O(N^{1/4}) shifts | Fixed (D9, D10, §2.3, §4.1, §3) |
| Items cx2 found CORRECT | Lemma ZN (10985 checks); TW and ρ; the TW* transcription; the rational-2-torsion corollary; the factored-count repair; SR2; non-CM density; Lem 1 circularity; the affine doubly supersingular obstruction; the ERH multi-curve reading; the Bach–Charles weight-2 form; the Euler product; the GL₂ constants (exact rationals); cx1's E* example; root-count/QR equivalence for quadratics; E2/E4 code (independently reproduced at b = 12); no operational leakage |

## 9. Opus critics

**Literature critic: MAJOR PROBLEMS, fixes applied.** It read D–U (both versions), KK, OU, GM, JU–P, DPS, Bach–Charles,
E–C and Serre in primary text. Fixes:
- Prior art (OU §4 and Thm 9, GM 2005 §2, KK): §0 and §2.5 rewritten.
- The OU description corrected.
- KK Lemma 2 / Thm 2 credited for CRT-ℓ.
- Published D–U Thm 1 is projective-only.
- D–U's LNCS 1043 misprint noted.
- Stein Ch. 3 and Ch. 8.
- E–C at level 1.
- The full Bach 1990 reference, with ln² and the gcd case.
- DP 2020, DPS 2023 and JU–P 2024 added.
- "Misreads" softened, as follows.

**On nonsmooth's table** (critic's wording, adopted):
- Its row "any single curve under ERH (their Thm 5)" repeats arXiv Thm 5's own wording ("any elliptic curve"), which the
  proof does not support because it uses twist counts.
- Its "D from a curve and one twist (Thm 1)" is unproved for projective counts, because of the Lem 1 circularity.
- We do not edit that file; this note and the STATUS line flag both.

**Proof critic: ACCEPT WITH FIXES, all applied.**
- It found that no PROVED claim is wrong.
- Checks it ran:
  - recomputed the GL₂ constants and the S₃ table exactly;
  - ran an exact TW* enumeration over all x₀ for 44 small pairs, 26 of them in Exc (every TW* rate was ≥ 0.57 except
    cx1's X_E point, at 0.335);
  - re-verified the PH4 member and the D–U short-arc counterexample with PARI;
  - cross-checked every reported number against the raw data.
- Fixes applied:
  1. Per-prime twist notation E′_r.
  2. Exc_M for general M.
  3. Cor. EQ is via TW*, and the domain reads "CM or non-CM".
  4. TW* bound sharpened to ¼ − o(1), so ≤ 4 + o(1) expected trials.
  5. CM density details: degree-one primes, swap invariance, inert and ramified pairs.
  6. E2 C2 min ρ corrected: 0.125, not 0.166.
  7. The X_E member sentence corrected (v3's "primes of each twist order divide the other count" was false).
  8. E3c time ≈ 1707 s.
  9. E1 z is computed against ρ + bonus.
  10. CRT-ℓ balance qualifier.
  11. This entry.
  12. Removed the vacuous gcd(h, N) clause, and noted that the root-count decision holds for any separable cubic.
