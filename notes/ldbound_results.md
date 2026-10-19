# ldbound: Hypothesis LD, via a randomized resultant multiplier

Agent `ldbound`, phase 9, 2026-10-19.
- Code: `src/ldbound/` (`params.py`: exponents, margins, constants, the μ-optimization; `checks.py`:
  exact identity checks). Tests: `tests/test_ldbound.py`.
- GPT-6.1-sol:
  - blind round `codex/ldbound_1.{md,out}`, with its full report copied to `notes/ldbound_codex1_report.md`;
  - adversarial round `codex/ldbound_2.{md,out}`;
  - re-check of this final version `codex/ldbound_3.{md,out}` (§8).
- Opus critics: one on proofs, one on literature (§8).
- No prediction-testing experiment was run, so there is no pre-registration.

## Verdict

**Theorem A's LD input is now proved, for a slightly modified design (a‴_k).** The modification adds one
random parameter to the field. Its infimum constant is the same as before, (1 + 3/(d+1))/√2, which is
0.9192 at d = 9. At d = 9 the minimal k for c < 1 is 13 (for pcest's design: 12).

So Theorem A becomes **Theorem A‴, conditional only on LV6**, plus the explicit caveats pcest already
listed: HECM, or Theorem A′ without it; the irreducibility transfer; and the internal hypD/lemmaM
results. Hypothesis LD *for pcest's design (a″_k) itself* remains open, but it is no longer needed.

Status of the proof:
- it passed an Opus proof critic (no fatal error; fixes applied);
- it passed an Opus literature critic (no prior statement found; credits fixed);
- it passed a GPT-6.1-sol adversarial round (CORRECT WITH FIXES; applied).
- GPT-6.1-sol found the same design and proof independently in its blind round.
- These are internal, unrefereed results.

In plain terms:
- **Why LD was hard.** In (a″_k), an algebraic value F_A(P) is forced to lie in one residue class
  modulo the large rational value G(P): r₀ ≡ N(b/g₁)^d (mod G) (pcest Lemma 1.2).
  - Given the rational data, this class is deterministic.
  - It moves along each rational fibre together with the partner point's rational-smoothness indicator.
  - So LD needed a correlation estimate between two things: smooth numbers in sparse, moving, structured
    progressions, and smoothness of an affine sequence. No such estimate is known.
- **The fix.** Build f with F(g₀, g₁) = N·w instead of N, where w is a random prime of size slightly
  above G. This is still an NFS pair, since g₁^d f(g₀/g₁) ≡ 0 (mod N).
  - The forced class becomes r₀ ≡ w·N(b/g₁)^d.
  - Because w is independent of the rational data, Brun–Titchmarsh makes the law of r₀ at most a
    constant times uniform over the reduced classes mod G, whatever the rational data are.
  - The structured-residue problem (the critic's piece (ii)) disappears. Only an all-residue statement
    (piece (i)) is left.
- **Piece (i) is elementary at precision L^{o(1)}.** LD only asks for an upper bound within L^{o(1)}. So
  the all-residue second moment of smooth numbers modulo the family ℓ²nq follows from two classical
  tools:
  - divisor switching (pcest Lemma DS) for conductors ≤ B^{C₀};
  - the large sieve for larger conductors, using the average over the special-q primes ℓ for
    conductors divisible by ℓ.
- **The cost.** w lives in the constant coefficient a₀ ≈ Nw/M^d. It stays invisible in the algebraic
  norm if h ≥ 1 − dμ + ν, where M = N^μ and ν = 1/(4k). The base-size exponent μ is now free.
  - The best μ (GPT-6.1-sol, ldbound_2) gives θ = 3s + (3(d−1)/(d+1) + 1)ν for large k, against
    3s + ν for (a″_k).
  - The infimum over k is unchanged.

Labelled results:
1. **Design (a‴_k; μ), Lemma GB_w [proved; elementary]** (§1). Re-checks of FM, Theorem R, Lemma Kern,
   Lemma 1.3 and Prop 4.2 [proved] (§2).
2. **Lemma W (residue spreading) [proved; Brun–Titchmarsh]** (§3).
3. **Theorem L2_fam (all-residue second-moment upper bound for smooth numbers modulo ℓ²nq) [proved;
   DS + large sieve + ℓ-average]** (§4). This is critic piece (i). It is BDH-type only in shape: an upper
   bound at the main-term scale, with no saving.
4. **Theorem LD* for (a‴_k) [proved]** (§5): E[N(N−1)] ≤ (log N)^{O(1)}·S⁴ρ_r²ρ_F². With FM,
   PC ≤ (log N)^{O(1)}, so D_k(γ) holds with γ ≥ L^{−o(1)}.
5. **Theorem A‴ [conditional on LV6 + listed caveats]** (§6). The time is L_N[1/2, (1 + θ)/√2 + O(η₀) + o(1)],
   with θ = min_μ θ_μ(d, k) from the table.
6. **Design (a″_k)** (§7):
   - LD stays open [conjectured];
   - the exact sufficient missing correlation (Prop 7.1);
   - what the available averaging does and does not give (Prop 7.2) [proved facts plus assessment];
   - two barriers [proved]: a coefficient-only repair at fixed resultant cannot move r₀ (Prop 7.4); any
     resultant multiplier must have support ≥ t^{1−o(1)} (Prop 7.3), which forces h ≥ 1 − dμ + ν − o(1)
     within this method.

### Constants (L[1/2, c], Lenstra–Pomerance shape test; `src/ldbound/params.py`)

All values are at zero slack (η₀ → 0).
- "min k" means c < 1 holds for all sufficiently small η₀ > 0. At d = 9, k = 13 this needs
  η₀ < √2 − 1 − 16/39 ≈ 0.0040.
- The infimum needs η₀ → 0 and k → ∞.

| d | (a″_k): min k | (a‴_k), μ = 1/(d+1): min k | (a‴_k), best μ: min k | best μ, c at k = 30 | k = 50 | k = 100 | k = 1000 | inf_k |
|---|---|---|---|---|---|---|---|---|
| 7 | 13 | 45 | 21 | 0.9914 | 0.9838 | 0.9780 | 0.9728 | 0.9723 |
| 9 | 12 | 20 | 13 | 0.9393 | 0.9313 | 0.9252 | 0.9198 | 0.9192 |
| 11 | 13 | 17 | 13 | 0.9053 | 0.8963 | 0.8901 | 0.8845 | 0.8839 |
| 13 | 13 | 17 | 14 | 0.8975 | 0.8713 | 0.8649 | 0.8593 | 0.8586 |
| 19 | 16 | 18 | 17 | 0.8988 | 0.8519 | 0.8197 | 0.8138 | 0.8132 |

- **Best μ.** Whenever μ* := (1 + (d−1)ν)/(d+1) ≥ dν (so that the level floor dominates the Kern floor;
  this holds exactly when k ≥ (d² + 1)/4; `params.best_mu` handles all k), the best μ is μ*. Then h = μ* − (d−2)ν and θ = 3μ* + ν. At d = 9, k = 30 this
  gives μ* = 8/75 and θ = 197/600.
- **Theorem A′ (Pollard–Strassen instead of HECM), best μ:**
  - d = 9: 0.9864 (k = 50), 0.9766 (k = 100), 0.9680 (k = 1000), limit 0.9670;
  - d = 11: 0.9230 (k = 100), limit 0.9142;
  - d = 19: limit 0.8226.
  - A fully HECM-free statement also needs LV6 to meet its cost without HECM (pcest §5, open).

---

## 1. Design (a‴_k; μ)

Notation is pcest's (§1 there):
- L = L_N[1/2, 1], Λ = (log N log log N)^{1/2}, B = L[1/2, a], ϱ(x) = Ψ(x, B)/x;
- s = 1/(d+1), ν := 1/(4k);
- ℓ uniform among the primes of J = (Y, 2Y], with Y ≍ N^{2ν}; π_J := #(primes in J) ≥ Y/(2 log Y);
- S = L^{σ+o(1)}, V = V_ℓ = ⌊√ℓ·S⌋, with envelopes V₋ := ⌊√Y·S⌋ ≤ V_ℓ ≤ V_max := ⌊√(2Y)·S⌋;
- 𝒫_ℓ = {(a, b) : b prime in (V/2, V], |a| ≤ V, gcd(a, b) = 1}, G_g(x, y) = g₁x − g₀y,
  Λ_ℓ(g) = {P : ℓ | G_g(P)};
- trial division to B as preprocessing.

"L^{o(1)}" means exp(o(Λ)). Fixed powers L^ε and N^{o(1)} are not L^{o(1)}.

**Parameters.** Fix:
- odd d ≥ 3, k > (d+1)/4, a > 0, a compact σ-range, η₀ > 0;
- a base exponent μ with ν + η₀ ≤ μ < 1 (μ = s is pcest's choice);
- h ≥ h_*(μ) + η₀, where
    h_*(μ) := max(2ν [Kern], μ − (d−2)ν [FM level ½], 1 − dμ + ν [window]);
- H := ⌊N^h⌋, M := ⌊N^μ⌋, and W := ⌈N^{μ+ν+η₀/2}⌉.

Then θ = θ_μ := (μ − ν) + (μ + h + dν) = 2μ + h + (d−1)ν. For μ = s, h_*(s) = s + ν and θ = 3s + dν.

**Field law 𝒜‴ (the moment law).** The variables ℓ, g, w, c are independent:
- g is uniform on 𝒢 = {(g₀, g₁) ∈ [M/2, M]² ∩ Z² : gcd(g₀, g₁) = 1};
- **w** is uniform among the primes in (W, 2W];
- c is uniform on C_H = [−H, H]^d ∩ Z^d.

Then F_A := F̂_{g,w} + G_g·R_c, where F̂_{g,w} is pcest Lemma GB applied to the integer N·w. So
F_A(g₀, g₁) = N·w.
- **Gcd checks.**
  - If gcd(ℓ, N), gcd(w, N), gcd(g₁, N) or gcd(G(P), N) is > 1, it is a proper factor, since all four
    numbers are < N, and the algorithm halts.
  - Since w is prime, gcd(w, N) > 1 means w | N.
  - ψ: x ↦ g₀g₁^{−1} (mod N) is well defined on Z[x]/(f), since g₁^d f(g₀/g₁) = Nw ≡ 0.
- **Candidates and counting** are exactly as in (a″_k): the special-q lattices, the filter
  |G(P)| ∈ [MV/32, 2MV], the capped search and the filtered count N_ℓ(A). We set N_ℓ(A) := 0 when ℓ | N,
  because the algorithm halts.
- **Order of operations (codex ldbound_2).**
  - The moment bounds (FM, LD*) are proved under the raw product law 𝒜‴.
  - Rejection of reducible f (and of non-unit leading coefficients, if LV6 requires that) is applied
    afterwards. It costs ≤ max N_ℓ·Pr(rejected) in the first moment and only decreases N(N−1), by
    positivity.
  - Theorem K's field selection is applied after that, as in pcest.
  - Conditioning on acceptance *before* Step 4 of §5 would destroy the product law, and we never do it.
- Sampling w uniformly among primes: draw uniform integers in (W, 2W] and test with AKS. The expected
  number of trials is O(log N).

**Lemma GB_w [proved].** Let gcd(g₀, g₁) = 1, g₁ ≥ 2, gcd(g₁, N) = 1. Lemma GB applied to T_d := Nw
gives:
- Σ a_i g₀^i g₁^{d−i} = Nw;
- |a_i| ≤ g₁/2 for 1 ≤ i ≤ d;
- |a₀| ≤ 2^d Nw/M^d + d2^{d−1}M;
- a_d ≡ Nw g₀^{−d} ≢ 0 (mod g₁).
*Proof.*
- The recursion and identity are GB's, verbatim.
- a₀ = (Nw − Σ_{i≥1} a_i g₀^i g₁^{d−i})/g₁^d, with g₁ ≥ M/2 and |a_i|(g₀/g₁)^i ≤ 2^{i−1}M.
- gcd(g₁, Nw) = 1, because w > M ≥ g₁ is prime. ∎
(Tested: `field_w`, `test_field_w_resultant_and_lemma12`.)

**Sizes [proved arithmetic].** For P ∈ 𝒫_ℓ:
- |F̂_{g,w}(P)| ≤ (dM/2 + |a₀|)V^d ≤ C_d(M + Nw/M^d)V^d. Since Nw/M^{d+1} ≤ 3·2^{d+1}N^{1−dμ+ν+η₀/2} ≤ H·N^{−η₀/2}
  by the window floor, we get |F̂_{g,w}(P)| ≤ C_d MHV^d·N^{−η₀/2}.
- |G(P)⟨c, v(P)⟩| ≤ 2dMHV^d.
- Put X* := C′_d MH V_max^d, with C′_d ≥ C_d + 2d. Then every F̂(P) + G(P)z with |z| ≤ T := dHV^{d−1} lies in
  [−X*, X*], uniformly in ℓ, g, w and P.
- X_F ≍ X* and X_r := 2MV/ℓ. So ρ_F = ϱ(X_F), ρ_r = ϱ(X_r) and θ = 2μ + h + (d−1)ν.

**Lemma 1.2_w [proved].** g₁^d F_A(P) ≡ Nw·b^d (mod G(P)) (homogeneity, as in pcest). (Tested.)

**Coprimality [proved; scope corrected after review].** Let ℓ ∤ N, P ∈ 𝒫_ℓ ∩ Λ_ℓ(g) with Y_P = 1 (that
is, n := G(P)/ℓ is B-smooth), and t := |G(P)|. Then:
- (b, t) = 1, because b is a prime > B, b ≠ ℓ, and n is smooth.
- (g₁, t) = (g₁, g₀b) = (g₁, b) = 1, because (g₀, g₁) = 1 and b ∤ g₁ (else b | t).
- (N, t) = 1, by preprocessing and ℓ ∤ N.
- (w, t) = 1, because w > W > 2MV_max ≥ t is prime.
- Hence the moving residue r₀ := F̂_{g,w}(P) mod t satisfies
    r₀ ≡ c_P·w (mod t),  c_P := N b^d g₁^{−d} ∈ (Z/t)^×,   (1.1)
  and **c_P depends only on (ℓ, P, g), not on w.**
- Without Y_P = 1 the conclusion can fail; codex gives an example with (g₁, t) = b. Such candidates
  never enter the moment sums.

## 2. Re-check of pcest's proved results for (a‴_k; μ)

- **Lemma 1.3, Lemma Q, Lemma DS, Theorem R, Lemma Kern.** These involve only the rational form, the
  candidate set and the linear map c ↦ (⟨c, v⟩, ⟨c, v′⟩). Their proofs use M only through
  M/V ≥ N^{μ−ν} → ∞ polynomially. **Unchanged** [proved there].
  - Lemma Kern needs H ≥ V², which the Kern floor gives.
  - Theorem R needs S² ≥ ρ_r^{−2}, which holds as in pcest.
- **Prop 4.2.** It contains F̂ through its windows and lift classes. Its proof is uniform under arbitrary
  translations of the windows, so it transplants verbatim. §5 uses only Kern(iii) and Kern(i) directly.
- **Theorem FM for (a‴_k; μ) [proved].** pcest's proof applies for each fixed admissible w, with
  constants uniform in w. Averaging over w then gives FM for 𝒜‴. Line by line:
  - *Step 1 (fibres, coprimality, Möbius).* Rational only. It needs M/V ≫ N^{o(1)}·2^{ω}, which follows
    from μ ≥ ν + η₀.
  - *Step 2(a).* Use (1.1).
  - *Steps 2(b), 2(c).* These concern c only. They need h > ν + η and H > dωV, both implied by the Kern
    floor.
  - *Step 2(d).*
    - The sign argument needs |F̂| ≤ ty/4 for y ≥ Y_ω = (2H+1)V′^{d−1}/ω. With t ≥ MV/32 and the size bound,
      |F̂|/(tY_ω) ≤ C_d ω N^{−η₀/2} → 0.
    - The additive error 3(|F̂|/t + 1)/(y + 1) is ≤ C_dωN^{−η₀/2}, which is polynomially small against
      ρ_F = L^{−O(1)}.
    - The grid 𝒳 uses X_F^max := X*. It is independent of t, P, ℓ, w and the residue.
  - *Step 3 (Lemma 8.5′).* The moduli t ∈ [MV/32, 2MV] have exponent q₀ = μ + ν, and x = X_F has exponent
    ξ = μ + h + dν.
    - q₀ > 2ν + η₁ holds, since μ > ν.
    - q₀ < ξ/2 ⟺ h > μ − (d−2)ν, which is the level floor.
    - Goodness concerns all reduced residues, so it is uniform in w.
  - *Step 4.* Unchanged.
  - Conclusion: E_{𝒜‴} N_ℓ(A) ≥ S²ρ_rρ_F/(log N)^{O(1)}.
  - No step of FM uses F(g₀, g₁) = N exactly. It uses only Lemma 1.2_w, the coprimality (1.1) and the
    size bound (codex ldbound_2 concurs).
- **Irreducibility.** The LV Lemma 5.1 transfer reads as in pcest §1, with the affine family
  f̂_{g,w} + (g₁x − g₀)Z[x]_{d−1} and the value condition ĝ(g₀, g₁) ≡ Nw (mod r) for r ∤ g₁Nw.
  - It is still checked only at the level of statements, now for a constant coefficient of size
    ≈ Nw/M^d ≤ H·N^{−η₀/2}.
  - Codex's blind round gives a self-contained finite-field (Turán-type) sieve proof for its version of
    the design: Pr(reducible or degree < d) ≤ H^{−1/4+o(1)} uniformly in the other parameters
    (`notes/ldbound_codex1_report.md` §4). We record it, not as a replacement for the assumption.
  - Reducible f can split as F₁F₂ = Nw without factoring N, so they must be rejected.
- **Relations.** ψ(a − bα) = G(P)/g₁ (mod N), unchanged. There is no w-parity column: w never occurs in a
  relation, and w > B, w ∉ J.

## 3. Lemma W (the forced residue becomes spread)

**Cited.**
- **(BT)** Montgomery–Vaughan, "The large sieve", Mathematika 20 (1973) 119–134, Theorem 2:
  π(x + y; q, a) − π(x; q, a) < 2y/(φ(q) log(y/q)) for q < y. We use the specialization
  π(x; q, a) ≤ 2x/(φ(q) log(x/q)).
- **(PNT)** π(2W) − π(W) ≥ W/(2 log W) for W ≥ W₀.

**Lemma W [proved].** For all N ≥ N₀(d, k, μ, σ, η₀), every t ∈ [1, 2MV_max] and every a ∈ Z,
  Pr_w[w ≡ a (mod t)] ≤ C_W/φ(t),  with C_W := 24/η₀.
*Proof.*
- If (a, t) > 1 the probability is 0.
- Otherwise the count is ≤ π(2W; t, a) ≤ 4W/(φ(t) log(2W/t)). Divide by W/(2 log W).
- log W ≤ log N, and log(2W/t) ≥ (η₀/2) log N − O(Λ) ≥ (η₀/3) log N. ∎
- The interval form of BT would give about 12/η₀ (literature critic).
- Choosing W/(2MV_max) and H/W to be suitable fixed powers of L, instead of N^{η₀/2}, would make the slack
  o(1) at a polylogarithmic cost in C_W (codex). This changes no exponent.
(Exact toy check: `bt_ratio`.)

**Corollary W [proved].** Let (ℓ, P, g) satisfy the hypotheses of (1.1), and let Φ: Z/t → [0, ∞) be fixed
independently of w. Then
  E_w Φ(r₀) ≤ (C_W/φ(t)) Σ_{c ∈ (Z/t)^×} Φ(c),
since a ↦ c_P·a is a bijection of (Z/t)^×. ∎

## 4. Theorem L2_fam (critic piece (i))

**Notation.**
- 𝒮* := 𝒮_B ∩ [−X*, X*] ∖ {0}, with #𝒮* ≤ 3Ψ(X*).
- S(α) := Σ_{F∈𝒮*} e(αF), and 𝔔(m) := #{(F, F′) ∈ 𝒮*² : F ≡ F′ (mod m)}.
- X′₊ := 2MV_max/Y. For ℓ ∈ J, 𝒩_ℓ := {n ∈ Z : P⁺(|n|) ≤ B, ℓ|n| ∈ [MV/32, 2MV]}. So |n| ≤ X′₊ and
  #𝒩_ℓ ≤ 3Ψ(X′₊).

**Additive Parseval [proved; exact].** Every h/m is a unique reduced fraction a/m′ with m′ | m, so
  𝔔(m) = (1/m) Σ_{h mod m} |S(h/m)|² = (1/m) Σ_{m′|m} f(m′),  f(m′) := Σ_{(a,m′)=1} |S(a/m′)|².
The term m′ = 1 equals (#𝒮*)²/m. (Exact integer form via Ramanujan sums: `pair_count_ramanujan`.)

**Cited.**
- **(LS)** The large sieve: for a_F supported on K consecutive integers,
  Σ_{m′≤Q} Σ_{(a,m′)=1} |Σ a_F e(aF/m′)|² ≤ (K + Q²)Σ|a_F|².
  - Montgomery–Vaughan, Mathematika 20 (1973), Thm 1 (in the form N + δ^{−1}).
  - The sharp K − 1 + Q² is due to Selberg and to Montgomery–Vaughan (with P. Cohen); see Montgomery,
    Bull. AMS 84 (1978) 547–567.
- **(DS)** pcest Lemma DS: every residue, non-reduced included; modulus ≤ B^{C₀}; intervals of length
  ≥ XB^{−C₀} inside [XB^{−C₀}, X].
- **(Hi)** ϱ(z₂) ≤ 2ϱ(z₁) for z₂ ≥ z₁. **(C8.16)** ϱ(x) = ϱ(xB^c)(log x)^{O(⌈c⌉)} for fixed c ≥ 0.

Fix C₀ such that R₀ := B^{C₀} ≥ 64S²·L/(ϱ(X′₊)ϱ(X*)) for large N. This is possible since all three factors
are L^{O(1)}. Note R₀ < Y.

**Theorem L2_fam [proved].** Uniformly for 1 ≤ q ≤ 2S²,
  (***)  E_ℓ Σ_{n∈𝒩_ℓ} Σ_{m′ | ℓ²|n|q, m′>1} f(m′) ≤ τ(q)(log N)^{O(1)}·Ψ(X′₊)Ψ(X*)².
Consequently, with the envelope V₋,
  E_ℓ Σ_{n∈𝒩_ℓ} 𝔔(ℓ²|n|q)/φ(ℓ|n|) ≤ τ(q)(log N)^{O(1)}·Ψ(X′₊)Ψ(X*)²/(Y q (MV₋)²).   (4.1)

*Proof.* Recall ℓ > max(B, 2S²) ≥ max(P⁺(n), q), so ℓ ∤ |n|q.
- **Rearrangement.** f(m′) does not depend on ℓ, because X* is common to all ℓ ∈ J. So the left side of
  (***) is Σ_{m′>1} f(m′)ν(m′), with ν(m′) := E_ℓ #{n ∈ 𝒩_ℓ : m′ | ℓ²|n|q}.
- **Partition (independent of ℓ).**
  - (A): m′ has no prime factor in J, and 1 < m′ ≤ R₀.
  - (B): m′ has no prime factor in J, and m′ > R₀.
  - (C): m′ has a prime factor in J.
  - Since R₀ < Y, the classes are disjoint and exhaustive.
  - For (A) and (B), m′ | ℓ²|n|q ⟺ m′ | |n|q, for every ℓ.
  - For (C), m′ | ℓ²|n|q forces ℓ | m′. Then ℓ is m′'s unique J-prime factor, with exponent i ∈ {1, 2}. An m′
    with two J-prime factors has ν = 0.
- **(A).**
  - *Bound on f.* f(m′) ≤ Σ_{a mod m′}|S(a/m′)|² = m′Σ_c N_c², with N_c := #{F ∈ 𝒮* : F ≡ c (mod m′)}. DS on
    ±[X*R₀^{−1}, X*], plus the trivial bound X*R₀^{−1}/m′ + 1 on |F| < X*R₀^{−1}, gives
    N_c ≤ (X*/m′)ϱ(X*)(log N)^{O(1)}. This uses R₀^{−1} ≤ ϱ(X*) and m′ ≤ X*ϱ(X*). Hence
    f(m′) ≤ m′·max_c N_c·#𝒮* ≤ Ψ(X*)²(log N)^{O(1)}.
  - *Bound on ν.* Let m₂ := gcd(m′, q) and m₁ := m′/m₂. Then m′ | |n|q ⟺ m₁ | n, since gcd(m₁, q/m₂) = 1.
    If m₁ is not B-smooth, ν = 0. Otherwise, by (Hi) and (C8.16),
    ν ≤ 2Ψ(X′₊/m₁) ≤ 4(X′₊/m₁)ϱ(X′₊R₀^{−1}) ≤ (Ψ(X′₊)/m₁)(log N)^{O(C₀)}.
  - The map m′ ↦ (m₁, m₂) is injective, with m₂ | q. So Σ_{m′≤R₀} 1/m₁ ≤ τ(q)(1 + log R₀).
  - Hence (A) ≤ τ(q)(log N)^{O(1)}Ψ(X′₊)Ψ(X*)².
- **(B).**
  - Here m′ ≤ X′₊q. On dyadic blocks (R, 2R], R = R₀2^i: ν(m′) ≤ #{1 ≤ |n| ≤ X′₊ : m₁ | n} ≤ 2X′₊/m₁ ≤ 2X′₊q/R.
  - By (LS) with K = 2⌊X*⌋ + 1 and Q = 2R: Σ_{m′∈(R,2R]} f(m′) ≤ (3X* + 4R²)·3Ψ(X*).
  - So (B) ≤ 36X′₊qX*Ψ(X*)/R₀ + 192X′₊²q²Ψ(X*).
  - First term: 36X′₊qX*Ψ(X*)/R₀ ≤ Ψ(X′₊)Ψ(X*)², because R₀ ≥ 36q/(ϱ(X′₊)ϱ(X*)).
  - Second term: 192X′₊q² ≤ ϱ(X′₊)ϱ(X*)X*, because X*/X′₊ ≥ c_dHV₋^{d−1}Y ≥ N^h and q²/(ϱϱ) = L^{O(1)}.
- **(C).**
  - Only ℓ = ℓ(m′) contributes, so ν(m′) = π_J^{−1}#{n ∈ 𝒩_ℓ : m″ | |n|q}, where m″ := m′/ℓ^i.
  - On a block (R, 2R] with Y/2 ≤ R ≤ 4Y²X′₊q (O(log N) blocks), two bounds hold: ν ≤ 3Ψ(X′₊)/π_J, and, since
    m″ > R/(4Y²), ν ≤ 8X′₊qY²/(π_JR).
  - Using the first bound against 3X* and the second against 4R²:
      (C) ≤ π_J^{−1}[27 log N·Ψ(X′₊)X*Ψ(X*) + 192X′₊qY²·8Y²X′₊q·Ψ(X*)].
  - First term: X*/π_J ≤ Ψ(X*)/log² N, since π_J = N^{2ν+o(1)} and 1/ϱ(X*) = L^{O(1)}.
  - Second term: it suffices that X′₊q²Y³ log Y ≤ c·ϱϱ·X*. Since X*/X′₊ ≥ c_dHV₋^{d−1}Y, this follows from
    h + (d−1)ν > 4ν, i.e. h > (5−d)ν. The Kern floor h ≥ 2ν + η₀ implies this with margin, for all d ≥ 3.
- **(4.1).**
  - For n ∈ 𝒩_ℓ: m := ℓ²|n|q ≥ ℓqMV₋/32 and φ(ℓ|n|) ≥ cMV₋/log log N.
  - So 𝔔(m)/φ(ℓ|n|) ≤ C log log N·[(#𝒮*)² + Σ_{m′|m, m′>1} f(m′)]/(Yq(MV₋)²).
  - Sum over n, average over ℓ, and use #𝒩_ℓ ≤ 3Ψ(X′₊) and (***). ∎

(Exponent audit by codex ldbound_2, ratios to the target Ψ(X′₊)Ψ(X*)²:
- (B) first term: q/(R₀ρ_rρ_F);
- (B) second term: N^{−h−(d+1)ν}L^{O(1)};
- (C) first term: N^{−2ν}L^{O(1)};
- (C) second term: N^{(5−d)ν−h}L^{O(1)}.
All are small. `params.margins` checks the exponents.)

*Remarks.*
- **What is used.**
  - No Siegel–Walfisz input and no equidistribution are used.
  - The pointwise DS bound covers conductors ≤ B^{C₀} for any fixed C₀.
  - The large sieve's diagonal X*Ψ(X*) is beaten in two ways: by the divisibility sparsity
    ν(m′) ≲ X′₊q/m′ in class (B), and by the average over N^{2ν} primes ℓ in class (C).
- **The ℓ-average is essential in (C).**
  - For one fixed ℓ, the frequencies with denominator ℓ² contribute ℓ²Σ_c N_c² over classes mod ℓ². Its
    trivial bound is X*Ψ(X*) ≫ Ψ(X*)².
  - A single-modulus bound would need smooth numbers to be L²-equidistributed mod ℓ² = B^{Θ(u)}.
  - That is unconditionally out of reach: the primes ≤ B could lie in a subgroup of (Z/ℓ²)^× of index up
    to 1/ϱ(ℓ², B) = L^{Θ(1)}, since Burgess-type bounds are far above B = ℓ^{o(1)} (cf. Soundararajan,
    CRM Proc. 46 (2008); Harper 2012).
- **Codex's blind proof is different.** Its Theorem AR (`notes/ldbound_codex1_report.md` §3) proves the
  same statement through multiplicative characters: induced-character multiplicities, Harper's Prop. 1
  and zero-density results for the ℓ, ℓ² conductors, and exact gcd layers for non-reduced classes. Its
  worst exponent is the same, h + (d−5)ν > 0. The additive route above avoids the induced-character
  bookkeeping.
- **Prop B5 is not contradicted.** L2_fam is an L² statement over all classes. It never passes through
  max-class ratios.
- **Attribution (literature critic).**
  - The route is the classical large-sieve route to Barban–Davenport–Halberstam upper bounds (Gallagher,
    Mathematika 14 (1967); Montgomery, Michigan Math. J. 17 (1970)).
  - The rearrangement Σ f(m′)·#{moduli divisible by m′} is in the spirit of Hooley's switching
    (J. reine angew. Math. 274/275 (1975)) and of Harper 2025's "non-concentration on multiples".
  - Smooth-number BDH-type bounds over *all* moduli: Harper, arXiv 1208.5992, Thm 2 (preprint), with
    saving H(u)^{−c}. It does not imply L2_fam: the family's main term is Ψ²·L^{−Θ(1)}, below Harper's
    error term Ψ²L^{−o(1)}.
  - Upper bounds of this type over all moduli are folklore: Granville, MSRI Publ. 44 (2008) §4.2,
    attributes one to Granville–Vaughan.
  - Harper, J. London Math. Soc. 112 (2025) e70293, Cor. 1 requires x^{0.51} ≤ Q ≤ x and does not apply.
  - We found no statement of L2_fam itself. What is specific here is the restricted family, the L^{o(1)}
    bookkeeping and the ℓ-average.

## 5. Theorem LD* for (a‴_k; μ)

**Theorem LD\* [proved].**
- **Hypotheses.** §1 and §2, with S² = (log N)^{C₃}/(ρ_rρ_F) as in pcest Theorem 4.3, and N_ℓ(A) the filtered
  count under the moment law 𝒜‴.
- **Conclusion.** E[N(N−1)] ≤ (log N)^{O(1)}·S⁴ρ_r²ρ_F².
- The same argument applied to Z (Step 1) also proves the bound in pcest's normalization:
    E_{ℓ,w}Σ_{P≠P′, ℓ|Δ}|𝒢|^{−1}Σ_g Y_PY_{P′}(S²/q)κ_alg ≤ (log N)^{O(1)}S⁴ρ_r².
  This uses (S²/q)κ_alg ≍ Z/(4d²H²V^{2d−4}ϱ(T|G|)ϱ(T|G′|)) and ϱ(T|G|) ≥ ρ_F/(log N)^{O(1)} under the filter.
- With FM: PC ≤ (log N)^{O(1)}, so D_k(γ) holds with γ ≥ L^{−o(1)} (Prop PZ).

*Proof.* Put T := dHV^{d−1}, t := |G_g(P)|, t′ := |G_g(P′)| and q := |Δ|/ℓ ∈ [1, 2S²]. Pair sums run over ordered
pairs P ≠ P′ in 𝒫_ℓ ∩ Λ_ℓ(g), with ℓ | Δ and both filtered.
- **Step 0 (law).** For ℓ | N, N_ℓ := 0. For ℓ ∤ N, under 𝒜‴,
    E[N(N−1)] = E_{ℓ,g,w} Σ_{P≠P′} Y_PY_{P′}·Pr_c[F_A(P), F_A(P′) ∈ 𝒮_B].
  - Pairs with ℓ ∤ Δ contribute nothing (pcest Lemma Q, step 0).
  - Only F̂_{g,w} depends on w.
- **Step 1 (kernel).**
  - Lemma Kern(iii), with ‖v ∧ v′‖ ≥ |Δ|(bb′)^{d−2} ≥ |Δ|(V/2)^{2d−4}, gives Pr_c[Φ(c) = z] ≤ C_d/(H²V^{2d−4}).
  - Φ(C_H) ⊂ [−T, T]², and Γ = {z₂ ≡ κ^{d−1}z₁ (mod Δ)} (Kern(i)), where P′ ≡ κP (mod Δ). Hence
      Pr_c[both smooth] ≤ (C_d/(H²V^{2d−4}))·Z,  Z := Σ_{j mod |Δ|} A_j(P)A_{κ^{d−1}j}(P′),
    with A_j(P) := #{|z| ≤ T : z ≡ j (mod |Δ|), F̂(P) + G(P)z ∈ 𝒮_B}.
- **Step 2 (AM–GM).**
  - κ is a unit mod Δ, since P′ is primitive.
  - So Z ≤ ½(Q(P) + Q(P′)), with Q(P) := Σ_j A_j(P)², which depends on (P, g, w, |Δ|) only.
  - This is not a pointwise comparison (Z may vanish while Q > 0), but both sides have the same
    main-term scale.
- **Step 3 (enlargement).**
  - z ↦ F = F̂(P) + G(P)z is injective, with F ≡ r₀ (mod t). Also z ≡ z′ (mod |Δ|) ⟺ F ≡ F′ (mod t|Δ|), for
    either sign of G. All such F lie in [−X*, X*] (§1). Hence
      Q(P) ≤ 𝔔₀(t, |Δ|; r₀) := #{(F, F′) ∈ 𝒮*² : F ≡ F′ ≡ r₀ (mod t), F ≡ F′ (mod t|Δ|)}.
  - 𝔔₀ is a fixed function of r₀ that does not depend on w.
  - Under the filter, tT ≍_d X*, so the enlargement has bounded scale.
- **Step 4 (average over w).**
  - Both candidates have Y = 1, so (1.1) and Corollary W apply. The classes c partition the counted
    pairs, so
      E_w 𝔔₀(t, |Δ|; r₀) ≤ (C_W/φ(t)) Σ_{c∈(Z/t)^×} 𝔔₀(t, |Δ|; c) ≤ (C_W/φ(t))·𝔔(t|Δ|).
  - The same holds for P′.
  - The sum is invariant under (P, P′) ↦ (P′, P), with Δ ↦ −Δ and symmetric filter and weights. Hence
      E[N(N−1)] ≤ (C_dC_W/(H²V₋^{2d−4}))·E_{ℓ,g} Σ_{P≠P′} Y_PY_{P′}·𝔔(tℓq)/φ(t).   (5.1)
  - No joint law of (r₀(P), r₀(P′)) under the shared w is needed.
- **Step 5 (fibres; the partner's rational weight).**
  - Fix ℓ, P = (a, b) and P′ with |Δ| = ℓq. The algebraic factor 𝔔(tℓq)/φ(t) depends on g only through
    n := G_g(P)/ℓ, and Y_P = 1 ⟺ n ∈ 𝒩_ℓ.
  - For fixed n, the g with G_g(P) = ℓn are g = g* + j(a, b) in (g₀, g₁)-coordinates. The direction is
    primitive with |(a, b)| > V/2, so there are ≤ √2M/V + 1 values of j.
  - Along the fibre, G_g(P′) = G_{g*}(P′) − jΔ. So n′ = G_g(P′)/ℓ runs through a progression of difference q,
    in an interval of length ≤ A := q(√2M/V + 1).
  - The filter puts |n′| ∈ [X_r/64, X_r]. For each sign, enlarge the n′-interval to length
    max(|I_±|, X_rB^{−C₀}) inside [X_rB^{−C₀}, X_r], which needs R₀ ≥ 64.
  - DS (modulus q ≤ 2S² ≤ R₀; X_rB^{−C₀}/q ≤ M/V) then gives
      #{g ∈ 𝒢 : G_g(P) = ℓn, Y_{P′} = 1, P′ filtered} ≤ U := (3M/V₋)ρ_r(log N)^{O(1)},
    uniformly.
  - Hence Σ_g Y_PY_{P′}𝔔(tℓq)/φ(t) ≤ U·Σ_{n∈𝒩_ℓ} 𝔔(ℓ²|n|q)/φ(ℓ|n|).
- **Step 6 (partners).** For fixed P and q, at most 12 points P′ ∈ 𝒫_ℓ have det(P, P′) = ±ℓq. They lie on two
  lines parallel to P; on each, consecutive integer points differ by P, with |P| ≥ V/2, inside a segment
  of length ≤ 2√2V.
- **Step 7 (assembly).**
  - Bound all ℓ-dependent geometric factors by their envelopes. These are #𝒫_ℓ ≤ CV_max²/log V₋,
    |𝒢| ≥ cM², U and V^{−(2d−4)} ≤ V₋^{−(2d−4)}; V_max/V₋ ≤ 2. Then apply (4.1):
      E[N(N−1)] ≤ C·(H²V₋^{2d−4})^{−1}M^{−2}·V_max²·(M/V₋)ρ_r·Σ_{q≤2S²} τ(q)(log N)^{O(1)}Ψ(X′₊)Ψ(X*)²/(Yq(MV₋)²).
  - Here Ψ(X′₊) ≤ (MV_max/Y)ρ_r(log N)^{O(1)} and Ψ(X*) ≤ C MHV_max^dρ_F(log N)^{O(1)}, by (Hi) and (C8.16).
  - Σ_{q≤2S²} τ(q)/q ≪ (1 + log S)².
  - The powers give V⁴/Y² ≍ S⁴, up to constants from the floors.
  - Hence E[N(N−1)] ≤ (log N)^{O(1)}S⁴ρ_r²ρ_F².
  - **PC.** By FM, (E N)² ≥ S⁴ρ_r²ρ_F²/(log N)^{O(1)} and E N ≥ 1, so PC ≤ (log N)^{O(1)}. ∎

*Remarks.*
- **No large-values input appears.** The residue is integrated out over w (Step 4) before any arithmetic
  estimate is applied. The partner's rational weight is then handled by divisor switching along fibres
  (Step 5). That is possible only because the algebraic factor no longer depends on the fibre index.
- **Theorem R is not used.** The rational pair correlation is re-derived inside Step 5.
- **Uniformity.** All constants depend only on (d, k, μ, a, the σ-range, η₀). Every "for large N" is a
  polynomial margin (`params.margins`) or an L^{O(1)} comparison.
- **LD versus LD\*.** pcest's literal LD has right side L^{o(1)} × (rational pair count). We prove the
  bound Theorem 4.3 consumes. The two agree if the rational pair count is ≥ S⁴ρ_r²L^{−o(1)}; that lower
  bound is not needed, and we do not prove it.
- **Independent confirmation.** Codex's blind round found the same order of operations:
  - a pointwise AM–GM split κ_alg(z) ≤ ½(K_P(z) + K_{P′}(z)) under one shared multiplier;
  - the multiplier average;
  - fibre divisor switching.
  Its multiplier is z uniform in [1, N^h] with F = zF̂ + GR, rather than a prime w. Non-unit z then needs
  exact gcd layers, and FM needs a new scalar-kernel lower bound, because the offset is as large as the
  window (`notes/ldbound_codex1_report.md` §§4–5). The prime w, kept polynomially below H/M^{d−1}, avoids
  both.

## 6. Theorem A‴

**Theorem A‴ [conditional].**
- **Hypotheses.** Fix odd d ≥ 3, k > (d+1)/4, μ ∈ [ν + η₀, 1) and η₀ > 0. Let h := h_*(μ) + η₀, and assume
  θ := 2μ + h + (d−1)ν ≤ 1. Put a := θ/√2 and B = L[1/2, a]. Assume:
  - **LV6** for the law 𝒜‴, and for the field-selected law derived from it. This means LV §6 adapted to:
    pure-descent rows at height V; the non-monic pair (f, g₁x − g₀) with Res = ±Nw; a leading coefficient
    a_d ≡ Nw g₀^{−d} (mod g₁) that depends on w; possibly a restriction to unit leading coefficients mod N
    (codex blind §4 shows this costs only a logarithm in FM); the content of f (it divides Nw); the
    g₁-parity column; the cost (rows)^{2+o(1)}; and transcript invariance (T). LV §6 also supplies
    characters, exceptional square classes and complexity bounds. None of this is established at these
    heights by the published theorem.
  - **HECM's running-time theorem**, as in pcest. Theorem A′ replaces HECM by Pollard–Strassen, at the
    constants of the table; this also needs LV6 without HECM.
  - **The LV Lemma 5.1 irreducibility transfer**, at the level of statements, now with constant
    coefficient ≈ Nw/M^d.
  - **The internal, unrefereed hypD/lemmaM results** used by pcest's Theorem A: Theorems K, S′ and DD.
- **Conclusion.** Every N is factored in expected time L_N[1/2, (1 + θ)/√2 + o(1)].
- Optimizing μ gives the table. The infimum is (1 + 3/(d+1))/√2: d = 9 gives 0.9192, d = 11 gives 0.8839,
  d = 19 gives 0.8132.
- c < 1 holds for d = 9 and d = 11 from k = 13 on, for sufficiently small η₀.

*Proof.* pcest Theorem A's proof (steps 1–5) applies verbatim to (a‴_k; μ):
- Theorem 4.3(ii) is replaced by Theorem LD* + FM + Prop PZ, applied in the order of §1 (moments,
  rejection, selection).
- T_SQ = S²L^{o(1)} = L^{θ/(2a)+o(1)}.
- The cost is max(a + 1/√2, a + θ/(2a), 2a) = (1+θ)/√2 at a = θ/√2.
- Sampling w costs (log N)^{O(1)} per field. ∎

## 7. Design (a″_k): what remains open

**Prop 7.1 (a sufficient missing estimate for (a″_k)) [proved reduction; the estimate itself is open].**
- Steps 0–3 of §5 hold for (a″_k) verbatim, with w = 1 and the old X*. Step 4 is the only step that
  averages over w.
- Without Step 4, the algebraic factor depends on g through the fibre index j (via r₀), not only through
  n. So Step 5's decoupling is unavailable.
- LD* for (a″_k) follows from this estimate:

      E_ℓ Σ_{P∈𝒫_ℓ} Σ_{q≤2S²} Σ_{P′: |det(P,P′)|=ℓq} Σ_{n∈𝒩_ℓ} Σ_{j ∈ 𝒥(P,P′,n): Y_{P′}(j)=1} 𝔔₀(ℓ|n|, ℓq; c_{P}(j))
          ≤ L^{o(1)}·E_ℓ Σ_{P} Σ_q Σ_{P′} Σ_{n∈𝒩_ℓ} U·𝔔(ℓ²|n|q)/φ(ℓ|n|),

  where:
  - 𝒥(P, P′, n) is the set of admissible fibre indices (g ∈ 𝒢, both filters);
  - c_P(j) = N b^d (g₁* + jb)^{−d} (mod ℓn);
  - n′(j) = n′₀ ∓ qj.
- It is a correlation, along each fibre, between the smoothness of an affine sequence and the lift-L² mass
  at a "Möbius progression" of d-th powers. It is sufficient, not necessary.

**Prop 7.2 (what the averaging in (a″_k) gives) [proved facts; the conclusion is an assessment].**
- (a) For fixed rational data (ℓ, P, P′, g), r₀ is determined; c does not move it. Averaging over g moves
  r₀ and the partner indicator *together*, by (2.1) of codex's report: both depend on the same j.
  - This is not a proof that no averaging argument can work (codex ldbound_2: a deterministic function of
    random data can be equidistributed even after conditioning).
  - It only locates the coupling.
- (b) Marginally, the structured residues are harmless (codex blind, §1).
  - For admissible candidates, (b, g₁) = (g₁, t) = 1. Distinct admissible fractions b/g₁ ≡ u (mod t) are
    separated by ≥ t/M². So each u has O(1) representations: there are no Farey spikes.
  - u ↦ Nu^d is ≤ 2d^{ω(t)}-to-1, and Σ_{n smooth} A^{ω(n)} ≤ Ψ·L^{o(1)}.
  - So the marginal structured-residue bound, without the partner weight, follows from L2_fam (weighted
    form, codex blind (3.13)).
  - The pcest §4.3 premise "representation count ≈ K on a set of measure ≈ 1/K" holds only without
    primitivity or the filter.
- (c) Discarding the partner weight costs 1/ρ_r = L^{Θ(1)} in this argument. Keeping it requires Prop 7.1.
  - Codex gives a model in which marginal all-residue control is exact but the partner-weighted energy is
    1/ρ_r times too large. So marginal estimates alone cannot imply Prop 7.1.
  - The model is not an arithmetic disproof of LD for (a″_k).
- We know no tool, published or here, that gives Prop 7.1 [assessment]. LD for (a″_k) remains
  **conjectured**.

**Prop 7.3 (multiplier barrier) [proved; scope: designs with M = N^μ fixed and a whole-box value bound].**
- **Hypothesis.** A design randomizes the forced residue by F(g₀, g₁) = Nw, with w ∼ μ_w independent of
  the rational data, and wants E_w Φ(c_Pw) ≤ (λ/φ(t))ΣΦ for every Φ ≥ 0 on Z/t.
- **Conclusions.**
  - μ_w(w ≡ a (mod t)) ≤ λ/φ(t) (take Φ an indicator), so #supp ≥ φ(t)/λ and max|w| ≥ φ(t)/λ − 1 = t^{1−o(1)}
    when λ = L^{o(1)}.
  - F(g₀, g₁) = Nw with g₀, g₁ ≤ M forces max|a_i| ≥ Nw/((d+1)M^d). By norm equivalence on binary forms over
    the candidate region, sup_P |F̂(P)| ≥ c_d(Nw/M^d)V^d.
  - The whole-box bound X_F = O(MHV^d) then forces h ≥ 1 − dμ + ν − o(1), the window floor.
- **Consequence.** Together with the Kern and level floors (inputs of the method), the optimized μ of the
  table is optimal *within this method*.
- The increase over (a″_k) is Δθ = min_μ θ_μ − θ_{a″}.
  - For μ = s it is min(s − ν, (d−1)ν) (codex ldbound_2; the first draft's (d−1)ν holds only for large k).
  - For the best μ and large k it is (3(d−1)/(d+1))ν.

**Prop 7.4 (no coefficient-only repair at fixed resultant) [proved; codex blind §4].**
- If two integral forms of degree d take the same value N at the primitive point (g₀, g₁), their
  difference is divisible in Z[x, y] by G = g₁x − g₀y (Gauss's lemma, G primitive).
- So the family F̂ + G·R already contains every integral perturbation with F(g₀, g₁) = N.
- By Lemma 1.2, every member has the same coarse class N(b/g₁)^d (mod G).
- So more coefficient randomness at fixed resultant can never move r₀. Changing the resultant, as
  (a‴_k) does, is necessary for this route.

**Status of the critic's route.**
- (i): proved (Theorem L2_fam).
- (ii): proved marginally (Prop 7.2(b)); the partner-weighted version is open for (a″_k), and made
  unnecessary by Lemma W.
- LV*: not needed and not proved.

## 8. Reviews

**Which results passed which reviews.**

| Result | Opus proof critic | Opus literature critic | codex ldbound_2 (adversarial) | codex ldbound_1 (blind) | codex ldbound_3 (final) |
|---|---|---|---|---|---|
| Lemma GB_w, sizes, Lemma 1.2_w, (1.1) | CORRECT (F1 fixed) | — | CORRECT; (1.1) scope fixed | same device, independently | see §8.5 |
| FM re-check | CORRECT | — | CORRECT WITH FIXES | own FM proof for its design | see §8.5 |
| Lemma W, Cor. W | CORRECT (F1) | standard; citation OK | CORRECT; scope fixed | analogous (uniform z) | see §8.5 |
| Theorem L2_fam | CORRECT WITH FIXES (F2, F4) | not in the literature; routine combination; credits added | CORRECT WITH FIXES (envelopes, citation) | independent proof (Theorem AR) | see §8.5 |
| Theorem LD* | CORRECT (F1) | — | CORRECT WITH FIXES | independently proved | see §8.5 |
| Theorem A‴ | correct conditional reduction | prior-art search: none | correct conditional reduction | same conclusion | see §8.5 |
| μ-optimized design | — (introduced after) | — | proposed by codex | — | see §8.5 |

### 8.1 Opus proof critic (all accepted)
- **Bottom line.** No fatal error. All losses are polylogarithmic, or O(1) from C_W; there is no L^ε or
  N^{o(1)} loss. In fact PC ≤ (log N)^{O(1)}. The constants table was recomputed independently.
- **F1.** The ℓ | N case: Corollary W needs c_P to be a unit. Fixed by setting N_ℓ := 0 and the scope of (1.1).
- **F2.** Classes in L2_fam were defined relative to ℓ after the ℓ-average. Fixed by an ℓ-independent
  partition by J-prime factors.
- **F3.** A constant in FM Step 2(d).
- **F4.** A typo, Ψ(X*) → Ψ(X*)², in (B).
- **F5.** The spurious η₀ < 1/(8k) was dropped; w ≈ G·N^{η₀/2}.
- **F6.** "Verbatim" in Prop 7.1 was overstated.
- **F7.** The A′ table line.
- **Added caveats.** LV6 with a w-dependent a_d; the irreducibility transfer with a large a₀; the
  internal hypD/lemmaM results.

### 8.2 Opus literature critic (all accepted)
- **L2_fam.** No published statement of it or of anything implying it was found.
  - Harper 1208.5992 Thm 2 (all moduli, saving H(u)^{−c}) does not imply it.
  - Harper JLMS 2025 Cor. 1 needs Q ≥ x^{0.51}.
  - LV Lemma 8.5 is L¹. LV's "naturally provide" sentence refers to Harper's *all-moduli* BDH, and LV prove
    no L² statement.
  - Granville 1993, Fouvry–Tenenbaum 1996, Drappeau 2015, de la Bretèche–Drappeau 2020, Pascadi 2025 and
    Baier–Zhao 2006 do not cover it.
  - The proof is a routine combination of classical tools (attribution paragraph added in §4).
  - Renamed from "BDH_fam" to "L2_fam", since there is no saving.
- **Lemma W, prior art for multipliers.**
  - Lenstra–Pomerance 1992 ("the main device … is the use of multipliers").
  - Schnorr–Lenstra, Math. Comp. 43 (1984); Seysen, Math. Comp. 48 (1987); A. K. Lenstra 1988.
  - Knuth–Schroeppel kN in CFRAC (Morrison–Brillhart 1975) and MPQS.
  - kN in NFS polynomial selection: Murphy, thesis, ANU 1999 §2.3.2.
  - **We claim only:** using a large random prime multiplier of the resultant, of the size of the
    rational modulus, to spread the forced residue in a rigorous second moment. No earlier instance was
    found.
  - Brun–Titchmarsh as a pointwise density bound is standard.
- **Citations fixed.** The large sieve is Montgomery–Vaughan, Mathematika 1973, Thm 1, with the sharp
  constant by Selberg and Montgomery–Vaughan–Cohen (Montgomery 1978). "J. London Math. Soc. 1974" is
  their Hilbert-inequality paper. BT is Montgomery–Vaughan 1973, Thm 2.
- **Prior art for c < 1.**
  - None found, 1992–2026, conditional or unconditional, for every N. The rigorous record is
    Lenstra–Pomerance 1992 (c = 1); under GRH, A. K. Lenstra 1988 (c = 1).
  - Lee–Venkatesan (JNT 187 (2018)) obtain unconditional congruences of squares in L(1/3). Their
    factoring is conditional on a character-decorrelation conjecture, and only for N = pq with
    p ≡ q ≡ 3 (mod 4).
  - LV write that they "use randomness to remove dependence on second moment bounds for which proofs known
    to us use the Riemann Hypothesis". Our second moment avoids RH for two reasons: we need it only to
    precision L^{o(1)}, and the w-randomization reduces it to an all-residue large-sieve bound.
  - This is a search result, not a certification.
- **Contradictions.** None found. Known obstructions concern individual moduli, and asymptotics exist
  only at Q ≥ x^{1/2+δ}.

### 8.3 GPT-6.1-sol adversarial round ldbound_2 (all accepted)
- **Verdict.** CORRECT WITH FIXES for LD*. Theorem A‴ is a correct conditional reduction.
- **Fixes applied.**
  1. (1.1) and Corollary W are restricted to Y_P = 1, ℓ ∤ N (with a counterexample otherwise).
  2. The moment law is defined before rejection and selection, with the transfer made explicit (§1).
  3. ℓ-independent envelopes V₋, V_max are used before the ℓ-average.
  4. The harmonic sum is ≪ (1 + log S)²; floors give ≍; the normalized form is "also proved"; Prop 4.2 is
     described accurately.
  5. Prop 7.2 is rewritten (it is not an impossibility theorem); Prop 7.3 is scoped, and its cost
     comparison corrected.
  6. Slack qualifications and the Pollard–Strassen values (d = 11: 0.9425, 0.9170, 0.9142).
  7. The LV6 sentence was too sweeping.
- **Improvement adopted.** Let M = N^μ vary. The window floor becomes 1 − dμ + ν, and at d = 9 the minimal k
  drops from 20 to 13 (§1, table). Re-checked by us (params, tests) and in ldbound_3.

### 8.4 GPT-6.1-sol blind round ldbound_1
- **Independent findings.**
  - The same multiplier device, with F = zF̂ + GR and z uniform in [1, N^h].
  - The same budget at μ = s: θ = 3/(d+1) + d/(4k) + η.
  - LD* proved for its design, with the same order of operations.
  - An independent character-sum proof of piece (i) (Theorem AR).
  - The admissible Farey multiplicity is O(1) (Prop 7.2(b)).
  - The partner-weight coupling is the residual obstruction for (a″_k) (Prop 7.2).
  - The Gauss-lemma barrier (Prop 7.4).
  - An FM proof for large offsets, and an irreducibility sieve.
- **Differences.** Codex uses uniform z (exact gcd layers needed) rather than prime w, and gets its
  piece (i) via characters and Harper's Prop. 1 rather than the additive large sieve.
- **Adopted.** Props 7.2(b) and 7.4, and the caveats on leading coefficient and content in LV6.

### 8.5 GPT-6.1-sol ldbound_3 (re-check of this final version)
*(see `codex/ldbound_3.out`; adjudication below)*

## 9. Candidate cards

```text
ID: KD-3 bounded-k Dixon descent, design (a'''_k; mu): (a''_k) with M = N^mu free and a random prime
  multiplier w ~ N^{mu + 1/(4k) + eta} of the resultant, F(g0, g1) = N w.
Target: every N (after trial division to B).
Proved: Lemma GB_w; FM, Theorem R, Lemma Kern re-checked; Lemma W (Brun-Titchmarsh spreading of the forced
  residue N w (b/g1)^d mod G); Theorem L2_fam (all-residue second-moment upper bound for smooth numbers
  modulo l^2 n q at precision (log N)^{O(1)}, via Lemma DS + large sieve + average over l); Theorem LD*.
Conditional on: LV6 for this law (Res = N w), HECM (or Pollard-Strassen, A'), irreducibility transfer,
  internal hypD/lemmaM results.
Cost: L[1/2, (1 + theta)/sqrt 2], theta = min_mu 2mu + max(2nu, mu - (d-2)nu, 1 - d mu + nu) + (d-1)nu;
  inf over k (1 + 3/(d+1))/sqrt 2 (d = 9: 0.919); c < 1 for d = 9 from k = 13.
Status: proved modulo the listed conditions; reviewed by two Opus critics and two GPT-6.1-sol rounds.

ID: B-mult (barrier) a resultant multiplier must have support >= t^{1-o(1)}, so h >= 1 - d mu + nu - o(1).
Status: proved (elementary). Scope: fixed-M designs randomizing the forced residue via Res(f, g) = N w,
  with a whole-box value bound.

ID: B-coef (barrier) at fixed resultant N, no coefficient randomization moves the coarse class N(b/g1)^d.
Status: proved (Gauss's lemma). Credit: codex ldbound_1.
```
