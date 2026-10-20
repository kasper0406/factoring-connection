# auditU: the two blind Opus critics (final reports, lightly condensed)

Both critics were read-only and blind. They read `notes/auditU_blind_claims.md`, `paper2/main.tex` and the internal
sources. They were told not to read STATUS.md, notes/lv5unif_*, codex/lv5unif_*, or any audit note. Their scratch scripts
are under the session scratchpad (`critic_analytic/`, `critic_algebraic/`) and are not in the repository.

---

## Critic 1: analytic number theory and uniformity — **ACCEPT WITH FIXES**

Scope: the analytic and uniformity layer only. LV6, Theorem DD and the recursion were not audited. The critic read [pp]
§§2–6, 8, 9, 11, [lem85] §§1–6 and [rec] §§3.2–3.3, and re-derived the key estimates. **Nothing in the relation-supply
chain loses exp(ω(log log N)). Every loss is d^{O(1)}, η₀^{−1}, (log N)^{O(1)} or Λ^{O(1)}.**

| Sev. | Location | Problem | Fix |
|---|---|---|---|
| MINOR | §4 Lemma FIB (Möbius error) | "2^{ω(ℓ\|n\|)} … negligible against M/V ≥ N^{η₀/2}" compares with the wrong quantity. 2^ω can reach N^{(log 2+o(1))(μ+ν)/λ}, i.e. log ≈ 0.7·D^{4/5}/d₀, while log N^{η₀/2} = D^{1/2}λ/2. So 2^ω ≫ N^{η₀/2}. The conclusion still holds, because the main term is N^{μ−ν−o(1)} and μ − ν ≍ D^{−1/5} ≫ μ/λ. `params.margins` already encodes the correct comparison. | Compare with M/V using μ − ν ≥ μ/2. Do not use the η₀ floor. |
| MINOR | Lemma W-U, (P2), (R′4), (P10) | Every N^{η₀} margin beats Λ only by the factor η₀ log N/Λ = D^{1/10}. Lemma W-U, for example, needs D^{1/10} ≥ 3σ + o(1). This is correct asymptotically, but N₀ is astronomical. The algorithm never uses N₀. | Say that each margin is D^{1/10}-thin. |
| MINOR | Lemma reg-U (common V) | "ℓ ≤ 2Y < V²/2" needs S² > 2, which fails at S = 1. The inequality actually used, ℓ < 2(1−1/d)V², holds for every S ≥ 1. | Reword. |
| MINOR | LD*-U Step 5 | Two signs times enlarged length give U = 4M/V·ρ_r(log N)^C, not 3M/V. | Constant only. |
| MINOR | §5 Lemma 8.5′(R′) | This is a list of substitutions into [lem85], not a rewritten proof. Each substitution checks. | Write SM(a), SM(b), Part 1, Part 2 and (i) out in full under (R′). |
| scope | §14, [lem85] | The analytic input of substance is still Drappeau's Lemme 5, whose proof points to unpublished Harper §3.3. It is used only for conductors ≤ y^K. Harper's range (log x)^K ≤ y covers y = L[2/5]. | None possible without re-proving Harper. |

**Weakest link: Lemma 8.5′(R′).** It is the one place where a published external estimate is pushed into a new regime
(y = L[2/5], u ≍ D^{2/5}). It rests on Drappeau Lemme 5, and its level gap has only a D^{1/10} margin. No error was found.

Key checks:
- **Window.** (bb′)^{d−2} ≥ e^{−2}V^{2d−4}, (1+1/d)^{d−2} ≤ e and (1−1/d)^{d−1} ≥ 1/e, checked exactly for 3 ≤ d < 400.
- **Envelope.** |F̂| ≤ 2^{2d+5}MHV^dN^{−η₀/2} ≤ MHV^dN^{−η₀/4}.
- **Kern-U.** True. Brute force on 162 random pairs (d = 3, 4, 5) gave 0 violations and a maximum constant of 1.176 < e.
- **LD*-U Step 7.** Re-derived with sympy: the ratio to S⁴ρ_r²ρ_F² is 4800e^{11/3}(2d+1)²/(d log(2S²Y)), times C_W(log N)^CΣ_q τ(q)/q.
- **Lemma CS.** Correct; the S^d trap is handled.
- **Tenenbaum's Ψ_q.** Its range is fine (ω(t) ≤ log N ≪ e^{c′D^{2/5}}).
- **CEP and HT.** Both are uniform in the range used.
- **Hidden fixed-d uses in [pp].** All of them are replaced: C_d in Kern, the Step 7 envelopes, c_k and 2^{−ω}, r₀(d) = 16, O_d in Lemma H and t_χ, and "O(1) primes of J divide N".
- **Lemma S″.** Brute-force illustrations on 4 small N passed.
- **C_γ.** Independent of N, d and k.

---

## Critic 2: algebraic, algorithmic and cost aspects — **ACCEPT WITH FIXES** (within scope)

Scope: the algebraic, algorithmic and cost layers. "Theorem U still stands or falls with the analytic relation-supply
chain … I did not re-prove that chain line by line, so this is not an ACCEPT of the whole theorem." No fatal or major error
was found.

| Sev. | Location | Problem | Fix |
|---|---|---|---|
| MINOR | §1 "Common V" | S_ℓ = V/√ℓ ∈ [S(1−o(1)), √2S], not [S/√2·(1−o(1)), √2S]. The stated range is weaker, so nothing breaks. Lemma reg's counts check: λ₂² ≤ (16/3)V² and ≤ 358S². | Correct the range. |
| MINOR | Theorem K-U | "2k dyadic windows": k windows suffice. The bound re-checks exactly. | Cosmetic. |
| MINOR | §10 cost table | The quoted optima are upper bounds, not grid optima. Finer grid: C_γ = 0 gives 3.5103, 1 gives 4.028, 5 gives 5.031, 20 gives 6.493. | Say "≤". |
| MINOR | Prop select-U(ii) | The displayed "(b/p)·2(…)" is garbled; the argument behind it is right. | Rewrite as an integral against the size-biased law. |
| MINOR | (10.1) | Charging 2k special-q batches to every attempt is wasteful; this is why a+s+max(a,σ) binds. | Optional amortisation. |
| MINOR | §9 / Fact LO | The Lagarias–Odlyzko theorem number is uncertain; the content used is right. | Pin the citation. |

**Weakest link: Theorem D-U** (γ ≥ (log N)^{−C_γ} uniformly in d). It needs LD*-U, Lemma 8.5′(R′) and GM-U, all outside
this critic's scope.

The critic's sanity point: with log x ≍ θ₀D^{4/5}λ, the smoothness bound B = L_x[1/2, a/√θ₀]. So the analytic regime is
the fixed-degree one, and the step from 1/2 to 2/5 adds only bookkeeping. The risk sits in [pp]'s second-moment machinery.

Key checks:
- **Stark's Lemma 8.** Applied correctly with c = 4(2d)!.
- **At most one exceptional class.** d is odd, so K has no quadratic subfield; L has at most one. Kummer injectivity and Landau finish it.
- **Bit size.** s = ⌈(2d)!(…)⁴⌉ gives ν₁s ≫ log s and log s ≈ 2d log d = o(Λ).
- **Characters.** t_χ = O(d² log N).
- **Square roots.** Lemma SQ costs e^{O(D^{1/5})}L^{a+o(1)}.
- **Lemma IR.** The bound 1/(4d) holds for r ≥ 16 and every d ≥ 3.
- **Fruitfulness.** Transcript invariance survives the recursive shape test.
- **Lemma S″.** Re-derived.
- **Recursion.** Lemma F_α holds, with η′ = η(1−2^{−α})/4.
- **Exponents.** θ₀, σ, the cost formula and the 2/5 floor were all re-derived.
- **Finiteness of c.** c is affine in C_γ and finite; it is not effective. The algorithm never uses C_γ.
- **Obstructions.** None known for this family.
- **Probe of Kern-U against the index-|Δ| image lattice.** The per-pair loss is ≈ S²/q, which sums over q to a log factor, consistent with V⁴/Y² ≤ 4S⁴.
