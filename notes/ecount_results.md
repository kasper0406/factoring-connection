# ecount: #E(Z/N) for one fixed curve, Dieulefait–Urroz, and partial counts (phase 23)

Agent `ecount`, 2026-10-22. Version **v2**: after the codex blind round, before the adversarial and critic rounds.
§8–§9 record those rounds.

- Pre-registration: `notes/ecount_prereg.md`, commit `d856b0b`, made before any experiment ran.
- Code:
  - OPERATIONAL: `src/ecount/zn_curve.py`. It never sees p or q, never calls integer factorisation and never imports the
    harness; an AST test checks this.
  - HARNESS: `src/ecount/harness.py` (PARI via cypari2 2.2.4, added to `.venv`) and `src/ecount/run_experiments.py`.
- Tests: `tests/test_ecount.py` (31 pass).
- Raw data: `experiments/raw/ecount_{e1,e2,e3,e4,posthoc,posthoc3,posthoc_xe}.jsonl` and `ecount_e1_summary.json`.
- GPT-6.1-sol:
  - blind round `codex/ecount_1.{md,out}`, abbreviated **cx1** (its scratch scripts are in `codex/scratch/ecount_1/`);
  - adversarial round `codex/ecount_2.{md,out}` (§8).
- Labels:
  - **KNOWN** (with reference);
  - **PROVED** (proof here or in cx1, re-checked here; "no prior formulation found by us" is the strongest novelty
    statement);
  - **CONDITIONAL**, **HEURISTIC**, **OPEN**;
  - **experimentally supported**;
  - **REFUTED (prediction)**;
  - **UNVERIFIED**.

Notation:
- N = pq with p < q primes, gcd(N, 6Δ_E) = 1.
- E: y² = f(x) = x³ + ax + b is a FIXED curve over Z, and E′ is its quadratic twist.
- A_r = #E(F_r) = r + 1 − a_r and A′_r = #E′(F_r) = r + 1 + a_r. λ_r and μ_r are the exponents of E(F_r) and E′(F_r).
- M = #E(Z/N) := A_pA_q, the projective count (E(Z/N) = E(F_p) × E(F_q) by CRT). M* = (A_p − 1)(A_q − 1) is the affine
  count.
- r_r = #{roots of f mod r}.
- Exceptional sets:
  - Exc(E) = {N : μ_p | M and μ_q | M};
  - X_E = Exc(E) ∩ {M odd} (cx1's set).

---

## 0. Verdict (plain language)

**The open target is closed, negatively.**
- For ONE fixed elliptic curve, the number #E(Z/N) alone factors N in randomized polynomial time. The inputs are
  ordinary or CM curves, with no factorisation of the count and no twist counts.
- Exceptions form a density-zero set X_E. X_E is **empty** for every curve with a rational 2-torsion point.
- So computing #E(Z/N) for a fixed curve is as hard as factoring (Theorems TW / TW*, §2). Supersingular-at-both-primes
  inputs are additionally handled by SR2.
- The selfred prototype "moved-divisor map" N ↦ #E(Z/N) (selfred §0, §6(a)) is therefore factoring-equivalent. It is
  not a route below factoring.

**The mechanism.**
- Take a random x₀, put δ = f(x₀), and use the point (δx₀, δ²) on the twist-by-δ curve over Z/N. No square roots are
  needed.
- At each prime this is a point on E or on its twist E′. Multiplying by M kills the E-components. A twist-component
  generically survives, and if it does not, the trailing doublings separate the 2-parts (Miller).
- A component mismatch shows up as a non-invertible denominator, i.e. a factor.
- This is the elliptic analogue of Miller's φ(N) ⇒ factoring. Note M/N = ∏_{r|N} L_r(E, 1)^{−1}.
- We (TW, with the Jacobi filter, §2.2) and cx1 (via division polynomials, with 2-stripping) found it independently.

**Dieulefait–Urroz, verified** (arXiv v1, 2019, and the published version, *Mathematics* 8 (2020) 2126; §1):
- Their single-curve result (Thm 4 published, Thm 6 arXiv) needs the **factorisation** of the count. It is a
  random-curve statement, and its divisor-count lemma (Prop. 1) has a gap.
- Both versions say no method is known to factor N from the count alone. TW answers that.
- Their two-count Theorem 1 is fine for affine counts, except the doubly supersingular case. For projective counts the
  proof is circular: Lemma 1/3 uses 4(p+1)(q+1), which already factors N. We give a deterministic repair when
  (d/N) = −1. When (d/p) = (d/q) = −1, TW gives a randomized repair.
- Their ERH theorem queries twist counts: it is a multi-curve oracle, not a fixed-curve result.
- The project's nonsmooth table misreads it as "any single curve under ERH".

**Equivalences (§3).**
- #E(Z/N) ≡ factoring.
- (a_N, a_{N²}) ⇒ factoring (Bach–Charles).
- Whether the ordered pair (a_p, a_q), or a_N alone, factors N is **OPEN**.

**Partial information (§4).**
- M mod ℓ with ℓ an input ≡ factoring.
- Fixed-ℓ counts are **OPEN**. Their values would leak fixed-base residuosity bits: Chebotarev advantage 9/128 for
  ℓ = 3, and 2/9 for M mod 2.
- For y² = x³ + 6x − 7 (non-CM), the root count of f mod N is *exactly equivalent* to the −3-residuosity bit (cx1).
- This is consistent with Theorem PT and refines Prop RANK.

---

## 1. Dieulefait–Urroz (verified from primary text)

### 1.1 Sources
- **arXiv:1911.11004 v1** (25 Nov 2019; the only arXiv version). We read the LaTeX source in full.
- **Published version:** L. V. Dieulefait, J. Urroz, "Factorization and Malleability of RSA Moduli, and Counting
  Points on Elliptic Curves Modulo N", *Mathematics* 8(12) (2020) 2126, doi:10.3390/math8122126. Received 29 Sep 2020,
  accepted 18 Nov 2020, published 27 Nov 2020.
  - cx1 located it. We read the PDF in full, as text extracted from the Universitat de Barcelona repository copy.
  - MDPI itself returned 403 to our fetcher.
- Numbering, published (arXiv): Thm 1 (1), Remark 1 (2), Lemma 1 (3), Thm 2 (4), Remark 2 (none), Thm 3 (5), Thm 4 (6),
  Prop. 1 (none), Thm 5 Coppersmith (7), Corollary 1 (none), Remarks 3–4 (8–9).
- Notation flips between sections. In §2, E′ or E* means affine. In §3 of arXiv, E_N means affine and E*_N projective.

### 1.2 What each statement says, and its status

| Published (arXiv) | Statement | Status here |
|---|---|---|
| Thm 2 (4), affine | E′(Z/N) = (p − a_p)(q − a_q) and E′_d(Z/N), (d/p) = −1 ⇒ factor N | **Correct** when (a_p, a_q) ≠ (0, 0). Published Remark 2 itself excludes a_p = a_q = 0: then every affine count of E and of its twists is N. E3a 180/180 + 180/180; E3b 60/60 counts = N. If a_q = 0 only, gcd(E′, N) = q directly |
| Lemma 1 (3), projective | "Knowing two among Ẽ, Ê, Ē and E, we know the four" | **Proof circular in both cases.** Case 1 uses Ẽ + Ē = 4PQ − L; Case 2 uses (M+1)Ê = 4PQ − L. With PQ = (p+1)(q+1) = N + p + q + 1 unknown, knowing it is already factoring. cx1 independently reports the same |
| Thm 1 (1) | Projective \|E(Z/N)\|, \|E_d(Z/N)\| and N pairwise distinct ⇒ factor N deterministically | **Unproved as published** (it rests on Lemma 1). Repaired in §1.3: deterministic when (d/N) = −1, randomized (TW) when (d/p) = (d/q) = −1. Deterministic repair of the latter: OPEN |
| Thm 3 (5), ERH | Factoring ≡ "counting points, affine or projective, of elliptic curves E mod N, non-supersingular for both primes" (published adds non-supersingular) | Proof: "apply Theorem 2 to the pair E, E_d for every d up to this bound". So it uses **twist counts** for d ≤ O(log²N): a multi-curve oracle. Projective inherits Lemma 1's gap; the repair needs a d with (d/N) = −1 (≤ 2log²N under ERH, Bach 1990) |
| Thm 4 (6) + Prop. 1 | Random curve E mod N; given the **factorisation** of the projective or affine count, factor N "with probability one" | The core is correct (cx1 agrees): a divisor u = q + 1 − a_q of the count with \|u − q\| ≤ 2√q + 1, then Coppersmith. It is deterministic and polynomial whenever d(count) ≤ poly(n). But: (i) the published Prop. 1 step "E ∈ S_N ⇒ E ∈ S_p or E ∈ S_q" is false as written (cx1). d(ab) > (log N)⁴ does not force d(a) > (log p)⁴ or d(b) > (log q)⁴. Raising the threshold to (log N)⁸ fixes this step. (ii) The probability is over a random curve in the isomorphism-class model, so for a fixed curve and adversarial N nothing is claimed. (iii) The arXiv v1 "probability one" rested only on the average order of d(n). E3c measured d(M) ≤ 9216 at n = 96 |
| "Small difference" section | Heron/circumradius: arcs of xy = T of length < (N/32)^{1/3} hold ≤ 2 lattice points | **Formula error (cx1, re-checked here).** The circumradius of (x_i, T/x_i) is R = √(∏_{i<j}((x_ix_j)² + T²)) / (2Tx₀x₁x₂). The paper prints the square of this, without the root. Near x = √T, R ≍ √T, not ≍ T, so the N^{1/3} conclusion is not established |
| Intro / Cor. 1 | "it is not known how to factor N only with the number E_N as input"; "there is no known polynomial time algorithm that can factor N using this information" | **Answered by TW / TW*** (projective count, randomized). The affine one-curve count remains OPEN |

### 1.3 Repaired projective two-count reduction (PROVED, routine)
- Let (d/N) = −1, say (d/p) = −1 and (d/q) = +1. Then E = (P − a)A_q and E_d = (P + a)A_q, with P = p + 1, a = a_p.
- S := E + E_d = 2P·A_q and D := E_d − E = 2a·A_q. Hence r := S/gcd(S, |D|) = P/k with k = gcd(P, a) ≤ 2√p
  (a ≠ 0; if a = 0, E = E_d, excluded by published Thm 1's distinctness).
- So p ≡ −1 (mod r), and Coppersmith/Howgrave-Graham on f(y) = r·y − 1 (made monic mod N), with |y| ≤ 2√p + 1 and
  O(1) sub-intervals, recovers k and hence p. Deterministic.
- `du_projective_twist_repaired`. E3d: 180/180 (n = 32, 64, 96). All were solved by the small-k gcd scan. The lattice
  branch is exercised by `test_du_projective_repaired_lattice_path` (6/6, n = 64, 96).
- Case (d/p) = (d/q) = −1: the data are E = A_pA_q and Ê = A′_pA′_q, i.e. two equations in s, a_p, a_q. We know no
  deterministic repair. TW* factors N from E alone.
- A randomized multi-twist variant (cx1, PROVED): counts of E, E_d, E_e, E_de with independent sign vectors sum to
  4(p+1)(q+1). The vectors are independent with probability 3/8 per attempt.

---

## 2. Theorems TW and TW*: one fixed curve, count ⇒ factoring (PROVED)

### 2.1 Lemma ZN: group law over Z/N (PROVED; standard in substance, cf. Lenstra's ECM)
- Elements are O (O in both components) or (x, y) affine in both components.
- `ec_add`:
  - gcd(x₂ − x₁, N) = 1: chord formula.
  - Proper gcd: output it.
  - gcd = N: put h = gcd(y₁ + y₂, N). If h = N, return O. If h = 1, use the tangent formula (2y₁ is a unit).
    Otherwise output h.
- The result is the CRT-correct sum, or a proper factor is output.
- *Proof.* Componentwise:
  - x₁ ≠ x₂ in both components: the chord formula is right in each.
  - x₁ = x₂ in both, so y₁ = ±y₂ in each. h = N means y₁ = −y₂ everywhere (sum O). h = 1 means y₁ = y₂ ≠ 0 everywhere
    (tangent).
  - Any other configuration differs between components, and the corresponding gcd is proper. ∎
- Corollary: `ec_mul(k, P)` (double-and-add) evaluates [k′]P for every binary prefix k′ of k. If any such [k′]P is O
  in exactly one component, it outputs a proper factor. Tested against separate mod-p and mod-q arithmetic.

### 2.2 Theorem TW (twist trick, with Jacobi filter; ours)
- Let M be any integer with λ_p | M and λ_q | M (e.g. M = #E(Z/N)).
- One trial: x₀ uniform in Z/N, δ = f(x₀). Output gcd(δ, N) if proper. If (δ/N) = −1, compute [M]P for
  P = (δx₀, δ²) ∈ E^(δ)(Z/N), where E^(δ): Y² = X³ + aδ²X + bδ³.
- Per-trial split probability:

  ρ(p, q) ≥ [(A_p − 1 − r_p)(A′_q − |H′_q|) + (A′_p − |H′_p|)(A_q − 1 − r_q)] / (4pq), with H′_r = E′(F_r)[M].

- Outside Exc(E): ρ ≥ 1/8 − O(p^{−1/2}).

*Proof.*
1. Write x₀ = (u, v) by CRT.
2. Suppose f(u) is a nonzero square mod p and f(v) a non-square mod q. Then P_p lies on a curve ≅ E/F_p, so
   [M]P_p = O.
3. Mod q, δ = c·s² with c a fixed non-square. The F_q-isomorphism to E′_c: c y² = f(x) sends P_q to (v, ±s). It
   preserves M-torsion.
4. The v with f(v) a non-square are exactly the abscissae of points of E′_c(F_q) with y ≠ 0, two points per v. If
   r_q > 0, A_q is even, so 2 | M and the y = 0 points lie in H′_q. Hence there are (A′_q − |H′_q|)/2 good v, and
   (A_p − 1 − r_p)/2 good u.
5. For such (u, good v), [M]P is O at p only, so Lemma ZN outputs p.
6. Add the mirror case.
7. If μ_q ∤ M, then |H′_q| ≤ A′_q/2. Hasse gives the bound. ∎

### 2.3 Theorem TW* (no filter, 2-adic descent; cx1's proof via division polynomials, transcribed to our algorithm and re-checked)
- TW* is TW without the Jacobi filter (`tw_trial(..., jacobi_filter=False)`). The point (δx₀, δ²) is used for every
  unit δ. At each prime it lies on E (δ a square) or on E′ (δ a non-square), independently, with the local law of
  cx1's formula (2).
- Claim: if N ∉ X_E, one trial splits with probability ≥ 1/8 − O(p^{−1/2}).

*Proof.*
- **Twist part.** If μ_r ∤ M for some r, the local annihilation probability at scalar M is
  α_r(M) = (A_r + |H′_r| − 2)/(2r) ∈ [1/2 − o(1), 3/4 + o(1)].
  - The split probability α_p(1 − α_q) + (1 − α_p)α_q is ≥ min(α_r, 1 − α_r) whatever the other prime does.
- **2-adic part.** If Exc holds but M = 2^s·m with s ≥ 1 and m odd:
  - Take r with A_r even and j ≥ 1 minimal with 2^j m annihilating both E(F_r) and E′(F_r).
  - At k = 2^{j−1}m, both kernels have index ≤ 4, because the cokernel is killed by 2 and E[2] has ≤ 4 points. At least
    one kernel is proper.
  - So α_r(k) ∈ [1/4 − o(1), 3/4 + o(1)], and the same min bound applies.
  - k is a binary prefix of M, so the chain computes [k]P (Lemma ZN). ∎
- **Corollary (rational 2-torsion; PROVED, cx1).** If E(Q) has a point of order 2, every A_r is even, so M is even and
  X_E = ∅. For such a fixed curve, #E(Z/N) and factoring are randomized-poly-time equivalent on all N of good
  reduction. Example: C3: y² = x³ + 2x + 3 = (x + 1)(x² − x + 3).

### 2.4 Corollary EQ, and the size of the exceptional sets
- **Corollary EQ (PROVED).** For a fixed E (ordinary or CM) and N of good reduction with N ∉ X_E, #E(Z/N) and factoring
  are equivalent under randomized polynomial-time reductions:
  - one oracle call and expected ≤ 8 + o(1) trials of Õ(n³) bit cost each;
  - Schoof/SEA in the other direction.
- **Doubly supersingular N.** M = (p+1)(q+1) and SR2 factors N deterministically. Such N are in Exc. TW* also handles
  them whenever M is even, which holds for r > 3.
- **Density zero (PROVED; cx1's argument, re-checked here).**
  - Let E be non-CM, and S a finite set of large primes ℓ for which the image of Gal in ∏_{ℓ∈S} GL₂(F_ℓ) is the full
    product (Serre, Invent. Math. 15 (1972); full product via Goursat and surjective det).
  - The event "ρ_ℓ(Frob_p) has eigenvalues −1 and t ≠ ±1, and ρ_ℓ(Frob_q) has no eigenvalue 1" has probability
    ((ℓ−3)/(ℓ−1)²)·(1 − O(1/ℓ)) ≍ 1/ℓ. It forces ℓ | A′_p and ℓ ∤ M, hence μ_p ∤ M and N ∉ Exc.
  - Chebotarev, applied to the iid pair (p, q) and to independent ℓ ∈ S, bounds the density of Exc ⊇ X_E by
    ∏_{ℓ∈S}(1 − c/ℓ). This tends to 0 as S grows.
  - CM curves: the same argument with split ℓ in the CM field (Serre 1972 §4.5).
  - No rate is claimed.
- **Census (E2, PH4).**
  - Exc fractions at b = 18 for the three non-CM curves: 3.9e−6, 7.8e−7, 1.9e−6. Most members are supersingular at one
    or both primes.
  - X_E has **one** member in the entire census: C1, (p, q) = (54443, 62819), b = 16, M odd.
    - A′_p = 3·13·19·73 and A′_q = 3²·5²·281. Each twist order's primes divide the other count.
    - TW* split it anyway, through a non-divisor binary prefix of M. This is outside the theory.
  - CM curves: X_E = ∅ throughout. Both CM test curves have a rational 2-torsion point, so M is always even.
  - cx1 exhibits a nonempty X_E for a contrived non-CM curve: y² = x³ + 27331x + 2 at N = 27331 = 151·181. Since
    27331 ≡ 0 mod N, this curve reduces mod N to the CM curve y² = x³ + 2.
- **Factored count ⇒ factoring for every N (PROVED, cx1).**
  - Strip each supplied prime ℓ of M (scalars M/ℓ^j).
  - Odd exponents λ_r ≠ μ_r always hold: if both equalled an odd L, then L | r + 1 forces both groups cyclic of
    order r + 1, which is even.
  - This repairs D–U Thm 4 (randomized, all N, fixed curve), without divisor enumeration.
- **Deterministic version: OPEN.** We need an x₀ with a given Jacobi pattern and an unkilled twist point. ERH bounds
  for (f(x)/N) along short intervals are not available.

### 2.5 Relation to the literature (to be checked by the literature critic, §9)
- Kunihiro–Koyama (EUROCRYPT '98, LNCS 1403) and Okamoto–Uchiyama (EUROCRYPT '98) reduce factoring to an oracle for
  #E(Z_n) of **arbitrary given** curves. Martín–Morillo–Villar (Appl. Math. Lett. 14 (2001)) use point orders.
- Dieulefait–Urroz need a twist count or the factorisation. Jiménez Urroz–Pomykała (Ramanujan J. 64 (2024); arXiv
  2210.04835) give a GRH randomized reduction with varying curves and twists (cx1; not re-read by us).
- The twist mechanism (x-only/Kummer arithmetic reaching E′) is standard in cryptography: twist security,
  invalid-curve attacks.
- The one-fixed-curve, count-only statement with the explicit exceptional set: no prior formulation found by us or by
  cx1.

---

## 3. Can #E(Z/N) be computed without factoring? Complexity and equivalences

**No, unless factoring is easy (Cor. EQ).** Any route that outputs M for one fixed curve is a factoring algorithm.
The cheapest known algorithm is to factor N (NFS, heuristic L_N[1/3, (64/9)^{1/3}]) and then run SEA on p and q.

- **Schoof over Z/N.**
  - Producing Frobenius on (Z/N)[x]/ψ_ℓ factors N when (N/ℓ) = −1 (nonsmooth FROB(d); cx1 re-derives it).
  - A method that output only the *numbers* (a_p mod ℓ, a_q mod ℓ, p mod ℓ) for enough ℓ would give M mod ℓ, hence M
    (Prop CRT-ℓ), hence the factors.
  - So trace-only Schoof output is also factoring-hard. This closes FROB's "numbers, not operators" scope gap for
    this output.
- **Hecke / modular symbols.**
  - a_N = a_p·a_q, whereas M = N + p + q + 1 − (p+1)a_q − (q+1)a_p + a_pa_q (cx1 (7)).
  - Hecke formulas sum over ≥ σ₁(N) matrices of determinant N.
  - Modular-symbol algorithms are not polynomial in log N (Stein, *Modular Forms: A Computational Approach*, GSM 79,
    Ch. 2, via cx1).
  - Also a_N = Σ_{x mod N} (f(x)/N) (cx1): N easy terms, which is exponential.
- **Euler product (PROVED, one line).** M/N = ∏_{r|N} A_r/r = ∏_{r|N} L_r(E, 1)^{−1}, with
  L_r(E, s) = (1 − a_r r^{−s} + r^{1−2s})^{−1}. So M plays the role of φ(N) for L(E, s), and TW is the analogue of
  Miller's reduction.
- **Bach–Charles** (arXiv:0708.1192 v2, Thm 2.2, Cor. 2.7–2.8; read by our fetch agent and cx1):
  - an algorithm for a(n) at every n, for a fixed non-CM eigenform, factors a density-1 set of RSA moduli;
  - the proof queries a(N) and a(N²);
  - weight-2 explicit form (cx1): S = a_N² + N − a_{N²} = qa_p² + pa_q², and qa_p², pa_q² are the roots of
    Z² − SZ + N·a_N²; take gcds;
  - it fails iff a_p = a_q = 0, and then M = (p+1)(q+1) still factors (SR2).
- **Oracle relations given N** (R = randomized, D = deterministic):
  - M ≡_R FACT, off X_E;
  - FACT ⇒_D (a_p, a_q) ⇒_D a_N;
  - (a_N, a_{N²}) ⇒_D FACT unless a_p = a_q = 0.
  - (a_p, a_q) ⇒ M or ⇒ a_{N²}: each needs p.
  - **(a_p, a_q) alone, or a_N alone, ⇒ FACT: OPEN.** We found no computable multiple of A_p = p + 1 − a_p from
    (N, a_p, a_q). These are "Galois-only" data with no algebraic occurrence of p.
  - #E(Z/N^k) = N^{k−1}M (cx1), so powers add nothing.
- **Lattice view.** M = (p + B)(q + A) with |A|, |B| ≤ 2N^{1/4} + 1.
  - Implicit factoring with shared MSBs needs t ≥ 2α + 3 bits (Faugère–Marinier–Renault, PKC 2010, Thm 1). Here
    t ≈ n/4 and α = n/2, so it is far out of range.
  - The general three-variable problem (cx1 (8)) includes M = N (A = B = 0), i.e. plain factoring, so no universal
    lattice theorem should be expected.
  - TW makes this moot for the projective count. For the affine count it is the remaining route (OPEN).
- **Edixhoven–Couveignes / p-adic methods.** These are prime-only. E–C Thm 15.2.1 needs the factored index (via cx1).
  An adaptation outputting M at composite N would be a factoring algorithm by TW*.

---

## 4. Partial information

### 4.1 Prop CRT-ℓ (M mod ℓ with ℓ an input ≡ FACT; PROVED)
- For balanced N (p < q < 4p/3, the generator's range), √p + √q ≤ 2.01N^{1/4}. Hasse then gives |M − N| ≤ W with
  W ≤ 4.1N^{3/4} + O(√N).
- Prime moduli ℓ with ∏ℓ > 2W + 1 give M by CRT. That is θ(B) ≈ 0.75n ln 2 + O(1), so about π(0.52n) calls.
- Then TW*. In the other direction, M mod ℓ ≤ FACT.

### 4.2 Fixed ℓ
- M mod ℓ = det(1 − g_p)·det(1 − g_q), where g_r = ρ_{E,ℓ}(Frob_r) and det g_r = r mod ℓ.
- **Predicting β_ℓ = (p/ℓ) = (ℓ*/p) on (N/ℓ) = 1 from (N mod ℓ, M mod ℓ).** Under Haar-independent g_p, g_q in
  GL₂(F_ℓ), the Bayes advantage is:

  | ℓ | 3 | 5 | 7 |
  |---|---|---|---|
  | advantage | 9/128 = 0.0703 | 0.0158 | 0.0058 |

  - Exact enumeration: `harness.chebotarev_adv_l`.
  - Exceptions: rational ℓ-torsion makes M ≡ 0 (C1 at ℓ = 5, increment 0 exactly, PH2). C3 has rational 2-torsion,
    so M is always even.
- **Prop PAR (ℓ = 2; PROVED).** Assume f is irreducible with Gal(f) = S₃.
  - M is odd ⇔ r_p = r_q = 0 ⇔ both Frobenii are 3-cycles. Then Δ_f is a square mod N: a one-sided certificate. In
    E4, 0 pairs had M odd and β = −1.
  - Bayes advantage on β_Δ given (Δ_f/N) = 1: Pr[M odd] + Pr[β = −1] − ½, which is 2/9 under S₃-Chebotarev.
  - #E(Z/N)[2] = (1 + r_p)(1 + r_q) carries the same information.
- **Exact equivalence with an open bit (cx1, PROVED; re-checked).**
  - Take E*: y² = x³ + 6x − 7 = (x − 1)(x² + x + 7) (j = 2048/3, so non-CM). Then r_r = 2 + (−3/r).
  - The root count R_N ∈ {9, 3, 1} is *equivalent* to the fixed-base bit "is −3 a square mod N" on N ≡ 1 (mod 3).
  - It is also equivalent to the information in φ(N) mod 3, and to |E*[2](Z/N)| ∈ {16, 8, 4}.
  - Also: root counting of x² − d mod N *is* the QR decision for d.
- **Status: OPEN** for fixed-ℓ M mod ℓ, M mod 2 and r_N(f). No sub-factoring algorithm is known, and no reduction
  from factoring is known. They are fixed-base residuosity-type bits, like φ(N) mod ℓ (selfred §2.1).
- **Consistency.**
  - Theorem PT: a pure tree with advantage δ on such a bit must split with probability ≥ δ/2. TW-type algorithms
    learn only by splitting. PT does not cover M-dependent schedules or Jacobi symbols, so there is no conflict.
  - RANK: linear divisor formulas give s mod ℓ, which decides β_ℓ. The GL₂ count M mod ℓ gives a strictly weaker,
    nonzero pair-class function, and for some curves (E*) an exact equivalent.

---

## 5. Experiments (pre-registered E1–E4; post-hoc PH1–PH4 labelled)

All runs stayed ≤ 450 MB RSS. Heavy runs ran alone. The light post-hoc runs (< 5 s, < 100 MB) overlapped E3.

### E1. TW (Jacobi filter) on random balanced semiprimes

Curves: C1 = 11a1, C2 = 37a1, C3: y² = x³ + 2x + 3 (non-CM); CM1: y² = x³ − x and CM2: y² = x³ + 1. Each (curve, n)
cell has 40 instances with n ∈ {32, 64, 96, 128}. Each instance gets 64 trials, and TW sees only (N, a, b, M).

| curve | solved per n (32, 64, 96, 128) | per-trial rate: empirical / exact ρ (z) | mean / max trials to first split |
|---|---|---|---|
| C1 | 40, 40, 40, 40 | .496/.498 (−0.2), .512/.500 (1.2), .498/.500 (−0.2), .501/.500 (0.1) | 1.9 / 7 |
| C2 | 40, 40, 40, 40 | .495/.493 (0.2), .481/.500 (−1.9), .514/.500 (1.4), .478/.500 (−2.3) | 1.8 / 7 |
| C3 | 40, 40, 40, 40 | .498/.493 (0.5), .505/.500 (0.5), .494/.500 (−0.6), .490/.500 (−1.0) | 2.0 / 7 |
| CM1 (10, 17, 10, 10 instances in Exc) | 40, 40, 40, 40 | .43/.24, .38/.19, .42/.24, .44/.24 (z = 21–24) | 2.4 / 16 |
| CM2 (10, 11, 9, 11 in Exc) | 40, 40, 40, 40 | .44/.27, .42/.26, .41/.24, .44/.26 (z = 18–21) | 2.5 / 15 |
| CM1, p ≡ q ≡ 3 (4) | TW 20/20 at each n; SR2 20/20 at each n | — | — |

- **P1: HOLDS.** 480/480 non-CM instances split, and none of them is in Exc.
- **P2.**
  - The ±3 SE clause holds in 12/12 non-CM cells (|z| ≤ 2.25). The one-sided "≥ ρ" clause fails in 5 non-CM cells by
    sampling noise; it was badly posed.
  - CM cells: **REFUTED as registered** (z = 18–24). The excess is the 2-adic prefix splitting of §2.3. Instrumented
    on CM1 (v₂(M) = 5): 75/100 Jacobi-(−1) trials split, at chain positions 4–5 from the end.
- **P3 (null M + 2 gives 0 splits): REFUTED. The control was mis-designed.**
  - M + 2 shares M's binary prefixes down to ⌊M/2⌋, so the chain evaluates [M/2]P and splits whenever λ_p | M/2.
    M + c with small public c is not information-free.
  - Corrected null PH1: 0 structural splits in 15360 trials.
  - Bonus (gcd(δ, N)) hits: 0 observed, ≤ 3·10⁻⁵ predicted. No ACCIDENTAL splits.
- **P4 (CM1 p ≡ q ≡ 3 gives 0 TW splits): REFUTED, in the favourable direction.** TW split 80/80 by 2-adic descent.
  SR2 split 80/80, as predicted.

### E2. Exact census of Exc(E) and ρ (all unordered pairs from S_b, b = 12–18; up to 1.4·10⁷ pairs per curve)

| curve | Exc fraction b = 12 / 14 / 16 / 18 | at b = 18: both-ss / one-ss / ordinary | min ρ outside Exc | median ρ |
|---|---|---|---|---|
| C1 | 9.0e−4 / 5.4e−5 / 1.3e−5 / 3.9e−6 | 45 / 5 / 5 | 0.124 (b = 16) | 0.500 |
| C2 | 0 / 0 / 1.5e−5 / 7.8e−7 | 6 / 2 / 3 | 0.166 | 0.500 |
| C3 | 0 / 4.3e−5 / 2.8e−5 / 1.9e−6 | 10 / 14 / 3 | 0.124 | 0.500 |
| CM1 | 0.271 / 0.247 / 0.250 / 0.251 | 3517878 / 1260 / 0 | 0.193 | 0.250 |
| CM2 | 0.241 / 0.252 / 0.249 / 0.256 | 3584503 / 1888 / 1 | 0.184 | 0.250 |

- **P5 (non-CM ≤ 10⁻³; non-increasing in b).**
  - The bound holds in all 12 cells, and C1 is monotone.
  - C2 and C3 rise from initial zeros (b = 12–14) before decreasing: **partly REFUTED**.
- **P6: HOLDS.** No pair outside Exc has ρ < 1/8 − 4/√p.
- **P7.**
  - Exc contains every both-supersingular pair: C(2653, 2) = 3517878 and C(2678, 2) = 3584503, exactly.
  - The "rest obeys P5" clause fails at b = 12 and 14 (one-ss members are 3.7e−3 and 1.5e−3 of the remaining pairs for
    CM1, 6.6e−3 and 1.5e−3 for CM2). It holds at b = 16 and 18: **partly REFUTED**.
- Ordinary–ordinary Exc members have smooth twist orders (e.g. A′ = 54720 = 2⁶·3²·5·19).

### E3. Dieulefait–Urroz on toy examples (C1–C3, 20 seeds per cell)

| part | cells | result |
|---|---|---|
| (a) affine two-count, case (d/N) = −1 and case d QNR at both | n = 32, 64, 128 | 180/180 and 180/180 (**P8 holds**) |
| (b) CM1, p ≡ q ≡ 3 (4): affine counts of E and E_d (d = 2, 3, 5, 6, 7) all equal N | n = 32, 64, 128 | 60/60 (**P9 holds**) |
| (c) Thm 4/6: factorisation of M → divisors in [√N/4, 4√N] → Coppersmith (flint LLL, dim 10) | n = 48, 64, 96 | FILL |
| (d) projective, repaired (Lemma 1/3 as written cannot run: it needs PQ) | n = 32, 64, 96 | 180/180 (**P11 holds**) |

### E4. Partial information, exact census, ordered pairs from S₁₆ (1491 primes; ≈ 1.1·10⁶ pairs per cell)

| curve | ℓ = 3 (pred 0.070) | ℓ = 5 (pred 0.016) | ℓ = 7 (pred 0.006) | ℓ = 2 (pred 2/9) |
|---|---|---|---|---|
| C1 | 0.072 | 0.017 (M ≡ 0 mod 5) | 0.014 | 0.213; M odd ⇒ β = +1 always |
| C2 | 0.050 | 0.019 | 0.014 | 0.214; same |
| C3 | 0.074 | 0.035 | 0.014 | 0.013 (M always even) |

- **P12 (within 0.01): 8/10 applicable cells.** It fails for C2 ℓ = 3 (−0.020) and C3 ℓ = 5 (+0.019): **partly
  REFUTED**. PH2 shows both lie inside the 5–95% finite-population null band.
- **P13 (C1 ℓ = 5 advantage 0): REFUTED as registered.** The 0.017 is the finite-population baseline from N mod 5 alone.
  The increment from M mod 5 is exactly 0 (PH2), which is the intended statement.
- **P14 (advantage = Pr[M odd]/2): REFUTED (prediction error).**
  - The Bayes identity is Pr[M odd] + Pr[β = −1] − ½, with Chebotarev value 2/9. Observed: 0.213 and 0.214.
  - The certificate direction is exact (0 counterexamples).

### Post-hoc (labelled; not pre-registered)
- **PH1 (corrected null).** TW with M′ uniform in [N − 8N^{3/4}, N + 8N^{3/4}] (non-CM, n = 64, 128, 40 seeds,
  64 trials): **0 structural splits in 15360 trials**. The information comes from M.
- **PH2 (E4 calibration).** Synthetic populations of 1491 primes with Haar GL₂(F_ℓ) classes, 200 replicates:

  | ℓ | null 5% / median / 95% |
  |---|---|
  | 3 | .044 / .071 / .098 |
  | 5 | .009 / .024 / .054 |
  | 7 | .008 / .025 / .049 |

  - Every observed advantage lies inside its band; null ranks 0.09–0.69.
  - The finite-population Bayes advantage at ℓ = 7 is dominated by the N mod ℓ baseline (the increment from M mod 7 is
    0 for C1 and C2). At this population size E4 cannot resolve a 0.006 effect.
  - At ℓ = 3 the increment over the baseline is 0.070, 0.048 and 0.071.
- **PH3 (TW vs TW* on Exc members, b = 16).**
  - Exc members tested: 14 (C1), 17 (C2), 31 (C3); 300-pair samples for CM1 and CM2.
  - TW split 13/14, 17/17, 31/31, 300/300 and 300/300. TW* split all, including the single M-odd member.
- **PH4 (size of X_E).** Over all 5 curves × 4 sizes (≈ 7.6·10⁷ pairs), X_E has exactly one member: C1, (54443, 62819).
  CM curves contribute 0.

---

## 6. Most decisive open questions
1. **Affine one-curve count.** Does M* = (A_p − 1)(A_q − 1) factor N? It is not a group exponent. It equals N on doubly
   supersingular N.
2. **(a_p, a_q) or a_N alone** (Galois-only data): do they factor N? Bach–Charles need a_{N²} as well.
3. **Fixed-ℓ partial counts** (M mod ℓ, r_N(f), #E(Z/N)[2]). These are fixed-base residuosity-type bits, with the same
   status as φ(N) mod ℓ: OPEN both ways. Some are exactly equivalent to the −3-residuosity bit (E*).
4. A **deterministic** one-curve reduction, and a quantitative rate for the density of X_E.

**For the project.**
- The selfred §6(a) prototype is settled: a poly-time one-curve count would itself factor N.
- What remains of the moved-divisor question is selfred §6(b) and (c), plus non-group "moved" quantities like M* and
  (a_p, a_q).

## 7. Reproduction and deviations
```
.venv/bin/python -m pytest -q tests/test_ecount.py                 # 31 tests, ~10 s
.venv/bin/python -m src.ecount.run_experiments e1                 # ~60 s
.venv/bin/python -m src.ecount.run_experiments e2                 # ~8 s
.venv/bin/python -m src.ecount.run_experiments e4                 # <1 s
.venv/bin/python -m src.ecount.run_experiments e3                 # FILL (Thm 6 Coppersmith dominates)
.venv/bin/python -m src.ecount.run_experiments posthoc            # PH1, PH2 (~2 s)
.venv/bin/python -m src.ecount.run_experiments posthoc3           # PH3 (~1 s)
.venv/bin/python -m src.ecount.run_experiments posthoc_xe         # PH4 (~2 s)
```

Deviations (disclosed):
- **D1.** PARI's multithreaded `polmodular` deadlocked under cypari2 in the first E1 run (C1, n = 128; 0% CPU, in
  `pthread_join`). We killed the run and reran it with `nbthreads = 1` and a 64 MB → 1 GB stack. Seeds are
  deterministic and the file was regenerated from scratch.
- **D2.** The registered M + 2 null is invalid (P3). PH1 replaces it.
- **D3.** P2's "≥", P13 and P14 were mis-formulated. They are reported REFUTED, not re-registered.
- **D4.** E3d's lattice branch is never reached on registered instances. A forced unit test covers it.
- **D5.** `du_affine_twist` follows D–U literally and does not first try gcd(E′, N).
- **D6.** TW* (`jacobi_filter=False`) was added after cx1. It is used only in PH3. The registered E1 used TW.
- **D7.** PH3's `tw_factor_h` counts only structural splits. If a δ-gcd hit came first (probability ~r_p/p per trial),
  a later structural split would be missed. This was not observed to matter.

(§8 codex adversarial round and §9 Opus critics: to be filled.)
