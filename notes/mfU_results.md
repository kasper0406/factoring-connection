# mfU: Theorem U-MF — a self-contained rigorous L[2/5] via the shared-rational-side multi-field design

Agent `mfU`, phase 30, 2026-10-28. The clean write-up for referees is **`notes/mfU_theorem.md`**; this note holds the
verdict, the audit, the dependency analysis and the reviews.

Files:
- `notes/mfU_theorem.md`: the theorem, the algorithm (pseudocode), every lemma with a proof or an exact citation, the
  dependency graph, and the cost ledger.
- `src/mfU/params.py`: exact leading-order margins (Lemma P, 18 checks including one negative control) and the cost
  constant c_MF with 2^{2/5} bracketed by rationals. `src/mfU/toy.py`: exact identity checks (base expansion and forced
  residue with composite w, the residue law of integer w, the norm-class identity of Lemma 8.5).
- `tests/test_mfU.py` (12 tests). `experiments/raw/mfU_cost.json` (cost-constant grid optima).
- Codex: `codex/mfU_1.{md,out}` (blind), `codex/mfU_2.{md,out}` (adversarial on the write-up), `codex/mfU_3.*` if run (§6).
- No experiments were run (only exact identity checks), so there is no pre-registration.

---

## 0. Verdict

**Theorem U-MF [proved here; unrefereed; AI-reviewed — see §6 for the final review status].** For fixed rationals
a, κ, d₀ > 0, FACTOR_MF is Las Vegas and completely factors every N ≥ 2 in expected time L_N[2/5, c_MF + o(1)], with an
explicit formula for c_MF in terms of one non-numerical constant C_γ.

In plain terms:
- **The multi-field route to L[2/5] is now written out in full, and it is simpler than the original Theorem U.** Its proof
  uses no algebraic second moment (no LD*, Kern-U, L2_fam, Lemma W) and **no Bombieri–Vinogradov-type estimate at all**.
- **The unpublished input is gone.** The earlier chain needed Drappeau's Lemme 5 (Compositio 2015), whose published proof is
  a pointer to §3.3 of Harper's unpublished arXiv:1208.5992. Theorem U-MF does not use it, nor any other smooth-number
  equidistribution in progressions. The replacement is elementary: draw the multiplier w uniformly from the *integers* of
  (W, 2W] instead of the primes (Lemma OP). Then the forced residue c_P·w of F(P₀) mod G(P₀) is uniformly distributed, and
  the smoothness probability becomes an interval density ≥ ρ_F/61. **GPT-6.1-sol found the same modification independently
  in its blind round** (codex/mfU_1 Lemmas 8–9, including the content normalisation and a height floor with constant 30).
- **The external inputs are all published numbered theorems**: Hildebrand–Tenenbaum 1986 (Thms 1, 2(i), 3), Canfield–
  Erdős–Pomerance, Rosser–Schoenfeld, Stark 1974 (Lemma 8), Landau (via Davenport), Lagarias–Odlyzko 1977 (Thm 1.3), AKS,
  Cantor–Zassenhaus, Kaltofen–Saunders, plus textbook algebra.
- **Simultaneous square certification is solved more cheaply than the triage expected.** After character filtering, the
  residual obstruction in each field is either trivial or the class of a *rational* number (Lemma 8.2, via Stark and
  Kummer theory in odd degree). The norm detects a rational class in odd degree. So π(B) + 1 shared "norm-parity" columns
  annihilate all J_f residual obstructions at once, and one kernel vector suffices (Lemma 8.5, Proposition PP). The
  triage's formulation also holds — at most J_f linear functionals (Remark in §8.6) — but that version needs J_f + 1 kernel
  vectors, at cost exponent 2a + 3σ_F instead of 2(a + σ_F).
- **The two predicted failure modes are handled explicitly.**
  - *Rejected fields absorbing the pointwise mass:* Lemma REJ bounds every rejection of a side, for each fixed rational
    side g, by e^{−(σ̄+1)Λ} ≪ ρ_F. This needs the MF height floor log H ≥ 30(σ̄ + 1)dΛ log log N (constant 30 instead of
    13, because certificate primes dividing a composite w are lost), and content is normalised, not rejected (its
    probability can be ≍ 2^{−d} ≫ ρ_F).
  - *Unaccounted certification cost:* it is in the ledger (§10.1 of the theorem). The rows grow to
    M_rows = L^{a+σ_F+o(1)}, because each of the J_f fields has its own prime-ideal columns, and linear algebra costs
    L^{2(a+σ_F)}.
- **The price is the constant, not the exponent.**
  c_MF = max(2(a + σ_F), a + σ_F + s + max(a, σ_r, σ_F), (a + σ_F + s)/(1 − 2^{−2/5})), with
  σ_r = 2/(5ad₀), σ_F = (2/(5a))(1/d₀ + d₀/(4κ)), s = 2κ(1 + C_γ).
  Since the recursion term always dominates, c_MF = (a + σ_F + s)/(1 − 2^{−2/5}) (codex/mfU_2). Grid optima: c_MF ≤ 9.45 / 10.85 / 13.52 / 17.41 for C_γ = 0 / 1 / 5 / 20, against Theorem U's 3.51 / 4.04 / 5.04 / 6.47.
  C_γ is not computed: it contains unspecified implied constants of Hildebrand–Tenenbaum. It is a different constant from
  Theorem U's C_γ.
- **Audit of the shared rational supply** (task 3, §3): no gap. The rational supply of MF is Proposition R1 (first moment),
  Proposition R2 (= [pp] Theorem R, the rational pair correlation) and Lemma DS. Three presentational gaps in [pp]'s proof of
  Theorem R are repaired in the write-up: squares straddling the filter boundary, the covering count, and q ≤ 4S² with the
  common box. All are MINOR. LD*, L2_fam and Lemma 8.5′ (used by Theorem U, not by U-MF) were re-read line by line, and
  no new gap was found. Their only external weak point is unchanged: Lemma 8.5′ rests on Drappeau's Lemme 5.
- **Weakest links** (my assessment; §6 has the reviewers'):
  1. *Everything is unrefereed.* It needs a human analytic number theorist, ideally on Lemma OP, Lemma REJ, Prop R2,
     Lemma 8.5 / Prop PP and Lemma 7.3.
  2. *Lemma 8.5 (the norm trick)* is new to this note. It relies on Lemma 8.2's conclusion that every exceptional class is
     rational. That follows from the same Stark/Landau argument the project already uses for |𝖤| ≤ 1, and the alternative
     J_f-functional certification is available as a fallback at a larger constant.
  3. *Asymptotics start astronomically late* (N^{η₀}-margins beat L[2/5] by D^{1/10}). This is not a proof gap, because the
     algorithm uses no threshold.
- **Not below 2/5.** The design changes the proof, not the one-level descent geometry (theorem §11.4).

**Prior work.** The closest is Lee–Venkatesan, JNT 187 (2018). Their c-randomisation of f and their §6 Stark–Landau analysis are used here (theorem §0, Lemma 8.2). Their Thm 2.1 gives rigorous L[1/3] for possibly trivial congruences; their Thm 2.3 gives fruitfulness only for p ≡ q ≡ 3 (mod 4) semiprimes under Conj. 7.1.

**Novelty.** We found no prior rigorous L[α] with α < 1/2 for factoring arbitrary N (search results of lv5unif §14 and auditU
§6, plus codex/mfU_1's literature table). This is a search result, not a novelty claim.

---

## 1. Task 1: the multi-field post-processing (summary; full proofs in the theorem, §8)

| Step | Statement | Where |
|---|---|---|
| row identity | r_i² ≡ 𝔮_i∏_ℓ f̃_{j_ℓ,d}^{−e_{iℓ}}ψ_{j_ℓ}(z_{P_ℓ})^{e_{iℓ}}, one shared m, a separate homomorphism ψ_j: Z[β_j] → Z/N per field | (8.1) |
| columns | shared rational sign + primes; **shared norm sign + norm parities (π(B) + 1)**; per field: parity, (r, s), t_χ^{(j)} characters | §8.2 |
| characters | t_χ^{(j)} = ⌈5d(ln|ℋ_j|_bd + ln(4J_f))⌉; all fields good simultaneously with probability ≥ 3/4; factorial-bit character primes cost e^{O(d log d)} = L^{o(1)} | Lemma 8.3 |
| obstruction count | on the good event, b_j(S) := [γ_{j,S} ∉ K_j^{×2}] is F₂-linear on the kernel: at most J_f functionals (the triage's lemma) | Remark §8.6 |
| **sharper** | each residual obstruction is a rational class [δ] (Lemma 8.2); Norm(δξ²) ≡ δ in odd degree, so δ_j(S) := sqfree(∏_{P∈j}F̃_j(P)^{e_P(S)}) satisfies δ_j(S)γ_{j,S} ∈ K_j^{×2}; the norm-parity columns make ∏_jδ_j(S) a square | Lemmas 8.2, 8.5 |
| certified roots | SQ_j on δ_j(S)g_j′(β_j)²∏z_P^{|e_P|}, exact verification, 2^{d−1} sign vectors | Lemma 8.4 |
| congruence | 𝖸 = √(𝔮_S/∏_jδ_j(S))·∏_j f̃_{j,d}^{−E_j/2}ψ_j(η_j)/(ψ_j(γ⁻_{j,S})ψ_j(g_j′(β_j))) | Prop PP |
| DD transcript invariance | every quantity is a function of (A, certificates, coins); character retries reuse rows (decisions are transcript-measurable); a trivial final gcd discards the rows | Lemma 7.3 |
| costs | M_rows ≤ J_f(1 + dπ(B) + O(d² log N)) + 2π(B) + 14 = L^{a+σ_F}; T_post = L^{2(a+σ_F)} (one KS kernel vector) | §8.2, Prop PP |
| MF height floor in the ledger | multiplies X_F by e^{30(σ̄+1)dΛλ}; changes log(1/ρ_F) by O(D^{1/5}λ²) = o(Λ): exponent +0 | theorem §10.1 |

*Why the norm trick works, in one line.* Stark puts an exceptional zero of ζ_{K(√𝔮)} in a quadratic subfield F. Since d
is odd, K(√𝔮) = KF, so 𝔮 ≡ δ_F mod K^{×2} (Kummer). The norm map sends δ_Fξ² to δ_F^d·□ ≡ δ_F. So the obstruction is
read off the norm — which is the product of the smooth values F̃_j(P), i.e. known rational data.

*Comparison with codex/mfU_1.* Codex's Lemmas 1–7 give the J_f-functional version (J + 1 kernel vectors, c_dense with
3(a + σ_F) or c_sparse with 2a + 3σ_F). It notes that the same exceptional quadratic field can recur across fields but
does not use the rationality of the classes. The norm-parity version is ours and is cheaper. Codex's ledger charges J·S_r²
algebraic evaluations per special-q, because it tests every rationally smooth candidate. FACTOR_MF tests only the first
one, P₀, which is all the proof uses, so its search costs Õ(B + S² + J_f).

## 2. Task 2: Theorem U-MF and its inputs

The statement, the algorithm (FACTOR_MF, SPLIT, SAMPLE_FIELD_MF, DESCEND_ONCE, SQ_SEARCH_MF, POSTPROCESS_MF, SQ_j), all
lemmas, the dependency graph (§11.1) and the list of unrefereed entry points (§11.2) are in `notes/mfU_theorem.md`.

### 2.1 Where an unrefereed or unpublished input enters

| Input | In Theorem U (lv5unif) | In Theorem U-MF |
|---|---|---|
| Drappeau 2015 Lemme 5 (proof = pointer to Harper arXiv:1208.5992 §3.3) | yes, via Lemma 8.5′ (bad moduli are rare) in FM-U | **no** |
| Iwaniec–Kowalski 7.13 / Montgomery–Vaughan large sieve | yes (Lemma 8.5′, L2_fam) | no |
| Brun–Titchmarsh | yes (Lemma W, prime w) | no |
| Tenenbaum's Ψ_q theorem | yes (good moduli) | no (all residues averaged) |
| de la Bretèche–Tenenbaum local bound | yes (Lemma 8.5′) | no |
| Lenstra–Pomerance 1992 | no (recursive shape test) | no |
| Hildebrand–Tenenbaum 1986 Thms 1, 2(i), 3; CEP; Rosser–Schoenfeld | yes | yes |
| Stark Lemma 8; Landau; Lagarias–Odlyzko Thm 1.3 | yes | yes (+ Kummer theory for Lemma 8.2) |
| AKS, Cantor–Zassenhaus, Kaltofen–Saunders, Stanley log-concavity, remainder trees | yes | yes |
| internal, unrefereed | [pp] §§3–9, lem85, recurse, lv5unif | `notes/mfU_theorem.md` only (self-contained) |

### 2.2 Can one avoid Harper's preprint inside the prime-w route? And what does a weaker input give?

- The integer multiplier makes the question moot for MF (and would also remove Lemma 8.5′ from Theorem U's first moment; its
  second moment LD* would then need Lemma W with all residues, which c12alg's DW analysis says goes through — not re-checked
  here, so recorded as plausible only).
- For the prime-w route, codex/mfU_1 checked the published candidates (statements read in primary texts or author copies;
  "verified" = statement inspected, not proof re-done). None replaces the input:
  - Fouvry–Tenenbaum, Proc. LMS 63 (1991), Thm 6 (UNVERIFIED): absolute error x(log x)^{−A}, useless against the sparse smooth main term.
  - Fouvry–Tenenbaum, Proc. LMS 72 (1996), Thm 1 / Cor 1: relative error, but only for y > exp((log x)^{2/3+ε}), while here
    y = L_x[1/2].
  - Granville, Acta Math. 170 (1993), Thm 2 (number UNVERIFIED): Q only x^{o(1)} at this y, and capped at √x/log^{B(A)}x.
  - Harper, JNT 132 (2012): individual moduli ≤ y^{4√e−δ}, not y^{Θ(u)}.
  - Harper, Compositio 152 (2016): exponential sums, upper bounds only — no lower bound in a prescribed residue class.
  - de la Bretèche–Drappeau, JEMS 22 (2020), Thm 2.1: its proof returns to Harper §3.3.
  - Pascadi, Compositio 161 (2025), Lemma 3.2: restates Lemme 5 with the same provenance.
- *Exponent from a weaker range* (codex/mfU_1 (25), re-derived here). An input valid only for log y ≳ (log x)^β, with
  d = D^δ, needs δ ≥ 1 − α/β, while the one-level geometry needs δ ≤ α/2. So α ≥ 2β/(2 + β). Fouvry–Tenenbaum's β = 2/3 + ε
  gives α > 1/2, so that route is no better than Lenstra–Pomerance. A precision loss τ instead gives (2 + τ)/(5 + 2τ),
  e.g. 3/7 at τ = 1.

## 3. Task 3: hostile audit of the shared rational supply

Read line by line: [pp] Lemma `lem:fibres`, Theorem `thm:FM` (rational part), Lemma `lem:DS`, Proposition `prop:R`
(Theorem R), Lemma `lem:reg`, Lemma `lem:pairs`, Theorem `thm:L2` (L2_fam), Theorem `thm:LD` (LD*), Lemma `lem:W`, and
Lemma 8.5′ in its (R′) form (lv5unif §5) together with Appendix A's E-lemmas. Verdicts:

| # | Sev. | Location | Finding | Status |
|---|---|---|---|---|
| A1 | MINOR | [pp] Prop R (Theorem R), covering | "Near the axes the filter removes everything" does not address squares that straddle the filter boundary \|n\| = X_r/128; DS needs its interval inside [XB^{−C₀}, X]. | Repaired in the theorem's Prop R2: each square's intervals are replaced by length-w intervals inside ±[X_r/256, X_r] (possible because w ≤ X_r/(8S²)). |
| A2 | MINOR | [pp] Prop R, count of squares | The O(qS²) covering count rests on an unproved claim about the parallelogram's side directions. | Replaced by area/w² + 4·perimeter/w + 4 with explicit side lengths. |
| A3 | MINOR | [pp] Prop R with a common box | "b > 2S² ≥ q" becomes q ≤ 4S² with one V for all ℓ; b is still a unit mod q because b > V/2 ≫ 4S² (P1). The case b = b′ is impossible when ℓ \| Δ. | Stated in Prop R2. |
| A4 | MINOR | [pp] Lemma fibres (ii) | The Möbius error 2^{ω(ℓ\|n\|)} must be compared with M/V, not N^{η₀/2} (auditU #2). | Done in Lemma FIB. |
| A5 | — | [pp] L2_fam classes (A)–(C) | Re-derived: the claim m′ \| \|n\|q ⟺ m₁ \| n holds, because gcd(m₁, q/m₂) = 1 (checked p-adically); (A) needs R₀ ≤ X_Fρ_F and R₀ ≥ 1/ρ_F (both true); (B) and (C) ratios as stated; (L2b) via Parseval with t ≥ MV/32. | No gap. |
| A6 | — | [pp] LD* Steps 0–7 | Step 1's lattice Γ_Φ and the count Σ_jA_jA_{κ^{d−1}j}; Step 3's injectivity z ↦ F; Step 4 needs only the marginal residue law of each point (AM–GM first); Step 5's fibre progression has difference q (G changes by −Δ per step); the exponent count V⁴/Y². | No gap (as auditU found). |
| A7 | — | Lemma 8.5′(R′) | Step C(a): Drappeau's sum restricted to primitive characters of modulus r ≤ R₁ is a sub-sum of nonnegative terms; E3/E8 exchange; E6 converges for β, β′ ≥ 3/4; (b) and Part 2 via LS′ with the density floor. | No new gap. Single external weak point: (D), Drappeau's Lemme 5 with Harper provenance (unchanged). |
| A8 | scope | Theorem U vs U-MF | The task's "shared rational supply (LD*, L2_fam, 8.5′)" is the supply of the *single-field* Theorem U. The MF route shares only the rational side, which is R1 + R2 + DS + FIB + CS. U-MF depends on A1–A4 (repaired) and on none of A5–A7. | — |

## 4. The constant (`src/mfU/params.py`, exact rationals)

`.venv/bin/python -m src.mfU.params experiments/raw/mfU_cost.json` (≈ 40 s): all 18 margin checks hold for (a, d₀, κ) =
(3/5, 1, 1/8), including the negative control: without the MF floor, the irreducibility exponent loses to log(1/ρ_F).
Upper brackets for the grid minimum of c_MF:

| C_γ | c_MF ≤ | (a, d₀, κ) | Theorem U (lv5unif, auditU grid) |
|---|---|---|---|
| 0 | 9.446 | (9/10, 1, 47/200) | 3.51 |
| 1 | 10.849 | (21/20, 3/4, 53/400) | 4.04 |
| 5 | 13.522 | (13/10, 1/2, 23/400) | 5.04 |
| 20 | 17.410 | (7/4, 1/4, 7/400) | 6.47 |

These illustrate the shape only. The C_γ of the two theorems are different, uncomputed constants. MF pays mainly through
M_rows = L^{a+σ_F} (J_f fields' ideal columns), which enters both the recursion term (a + σ_F + s)/(1 − 2^{−2/5}) and the
linear algebra 2(a + σ_F).

## 5. What was not done

- No toy implementation of FACTOR_MF (the identities it relies on are tested exactly; nothing asymptotic is testable).
- Theorem U's single-field LD* with integer w (§2.2) is not re-proved.
- No human referee.

## 6. Reviews

**Opus proof critic (read-only, blind to this note; on theorem v0): ACCEPT WITH FIXES, no FATAL or MAJOR.** It checked
every internal proof: Lemma OP (uniform |F̂_w| ≤ x₁; averaging over w commutes with the layer cake; all residues
reachable, no n–t coprimality needed; constants), Lemma REJ (both r ∤ g₁ and r | g₁ cases; brute-force irreducible
fractions over F₁₇, F₁₉ ≈ 1/d against the claimed 1/(4d); no forgotten rejection path; no loss in d), Prop MF, R1, R2,
CS, §8 (Stark + Kummer with 4 ∤ 2d; δ_j is the norm class; the sign and leading-coefficient bookkeeping of 𝖸²; units;
column sufficiency), Theorem K, S″, DD, Lemma 7.3, Prop select, the ledger and F_α. It found no exponential-in-d loss
inside γ, and all margins pass for three parameter triples. Its probability of correctness after the fixes: ≈ 0.8.
Weakest link: the exact general (non-normal) form of Stark's Lemma 8 behind Lemma 8.2.

| # | Sev. | Finding | Action |
|---|---|---|---|
| 1 | MINOR | Lemma 4.4: "g₁ ≤ M < ℓ" is false (M ≫ ℓ) | justified by gcd(g₀, g₁) = 1 |
| 2 | MINOR | Prop R2 needs Y > 32S², not 2S² < Y | (P1) strengthened; params check covers it |
| 3 | MINOR | Prop R2: "w ≤ X_r/(28S²)" wrong | w ≤ X_r(1 + o(1))/(8S²) ≤ X_r/256 |
| 4 | MINOR | Lemma DS needs Δ₀ ≥ 1 | stated (A/q ≥ XB^{−2C₀} ≥ 1) |
| 5 | MINOR | Landau needed with multiplicity | stated in §1.2 |
| 6 | MINOR | "one step each" alternation undefined with nested calls | one bit operation of the whole nested computation |
| 7 | MINOR | Lemma P constants implicit | explicit inequalities listed in (P2) |
| 8 | MINOR | exact forms of Stark Lemma 8, LO Thm 1.3, HT Thm 3 not checked against print | sent to the literature critic |
| 9 | MINOR | φ(N)/N should be φ(N)/(N − 1); loose 8 in D(x) | fixed / cosmetic |

**GPT-6.1-sol adversarial round `codex/mfU_2` (xhigh, on theorem v0.1, blind to this note): ACCEPT WITH MINOR FIXES; no
FATAL or MAJOR.** Probability correct after the fixes: 80%. Weakest link: Proposition R2 through Lemma DS (the rational
supply). It re-checked OP's layer averaging, REJ/MF conditioning, CS computability, Stark (p. 148; the factorial
threshold needs no normality), Landau, LO Thm 1.3 (primary text), Lemma 8.5's signs, negative exponents and leading
coefficients, PP, the size-biased selection, transcripts and recursion. It also observed that the third term of c_MF
always dominates, so **c_MF = (a + σ_F + s)/(1 − 2^{−2/5})** (now in the theorem; test added). Its findings, all MINOR and
all applied:

| # | Finding | Action |
|---|---|---|
| 1 | (P4) "d log S = o(ΛD^{1/5})" false (it is O) | corrected; the needed o(D^{4/5}λ) holds |
| 2 | Lemma 4.4 "g₁ ≤ M < ℓ" (same as proof critic #1) | fixed |
| 3 | Prop R2 w ≤ X_r/(28S²) (same as proof critic #3) | fixed |
| 4 | Lemma DS applies U1 slightly below the "in range" floor; U1's intermediate points can exceed N | DS requires X ≥ B^{D^{1/3}+2C₀+1}; range extended to N² |
| 5 | KS projection "coordinate at a fixed basis element" does not justify uniformity | F₂-linear π with π(1) = 1, solver sees only A^⊤x; argument written |
| 6 | standard algebraic inputs (units, trace dual, φ, ω bounds) not listed | listed in §1.2 |
| 7 | §11.4 norm exponent: algebraic part is 1/d + dν, total 2/d + dν | clarified |
| 8 | "expected running time is" should be "at most"; trial-division primes not AKS-certified | fixed |

Because round 2 was not MAJOR REVISION, no third codex round was run.

**Opus literature critic (read-only, on v0.1; relayed by the orchestrator): citations essentially sound, no proof breaks, one
important citation missing.** All fixes applied:

| # | Finding | Action |
|---|---|---|
| 1 | **Missing prior work: Lee–Venkatesan JNT 187 (2018).** Their (3.1)/§4 already randomise f by a large random multiple of (x − my), which is our SM1/OP c-randomisation. Their §6 (Lemma 6.10, Claim 6.25, Facts 6.22/6.24, Rem. 6.27) has the core of Lemma 8.2. | Cited prominently (theorem §0, Def. 4.1, Lemma 8.2). Lemma 8.2 is presented as an extension of LV §6. New: |𝖤| ≤ 1, and the rationality of the class (Kummer). LV instead bound the classes by (4/3) log log n and remove them by brute force. No prior use was found of the integer multiplier w or of norm columns detecting the rational exceptional class (BLP §5: a square norm is "necessary … but very far from sufficient"). |
| 2 | Record statement should contrast with LV Thm 2.1 (unconditional L[1/3, (64/9)^{1/3}], possibly trivial congruences) and Thm 2.3 (semiprimes p ≡ q ≡ 3 (mod 4), conditional on Conj. 7.1) | added to theorem §0; the claim "no rigorous L[α < 1/2] for arbitrary N" stands |
| 3 | CEP has no "Corollary 3.1": it is the unnumbered Corollary after Thm 3.1, with an explicit error term | label and statement corrected; §11.2 updated |
| 4 | Labels UNVERIFIED against print: LO "Thm 1.3" (content confirmed), HT Thm 2 "(i)", Stanley "Prop. 2", Granville Thm 2 (whose range also caps Q at √x/log^{B(A)}x), FT 1991 Thm 6 | marked UNVERIFIED |
| 5 | Kaltofen–Saunders Thm 3 (rank) is Monte Carlo | Las Vegas behaviour comes from our exact verification, and F₂-kernel uniformity from our own projection argument; stated |
| 6 | Landau: Davenport ch. 14 is fair; Montgomery–Vaughan MNT I ch. 11 states the multiplicity | added |

Confirmed: HT Thms 1, 3 (via the HT 1993 survey), Rosser–Schoenfeld Cor. 3, Stark Lemma 8 verbatim (including F ⊂ L in the
non-normal case), the Kummer and subfield arguments, AKS, Cantor–Zassenhaus, and the Harper/Drappeau dependency description.
Harper arXiv:1208.5992 is still unpublished as far as can be checked.
