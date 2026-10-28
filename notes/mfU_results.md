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
  Grid optima: c_MF ≤ 9.45 / 10.85 / 13.52 / 17.41 for C_γ = 0 / 1 / 5 / 20, against Theorem U's 3.51 / 4.04 / 5.04 / 6.47.
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
  - Fouvry–Tenenbaum, Proc. LMS 63 (1991), Thm 6: absolute error x(log x)^{−A}, useless against the sparse smooth main term.
  - Fouvry–Tenenbaum, Proc. LMS 72 (1996), Thm 1 / Cor 1: relative error, but only for y > exp((log x)^{2/3+ε}), while here
    y = L_x[1/2].
  - Granville, Acta Math. 170 (1993), Thm 2: Q only x^{o(1)} at this y.
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

(Filled in below as the reviews complete.)
