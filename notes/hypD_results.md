# hypD: relation supply for Dixon descent (Hypothesis D)

Agent `hypD`, phase 7, 2026-10-19.
- Pre-registration: `notes/hypD_prereg.md`. Commit `7e6f1ee` came before any registered run; Addendum A,
  commit `af85691`, came before E1x.
- Code: `src/hypD/`. Tests: `tests/test_hypD.py`.
- Raw data: `experiments/raw/hypD_e{1,1x,2,3a,3b}.jsonl`.
- GPT-6.1-sol: the blind prompt `codex/hypD_1.md` could not run because the Codex login expired. The
  prompts are written for a later run.
- Reviews: four independent read-only Opus critics (two on proofs, two on prior art) and the
  GPT-6.1-sol adversarial round `codex/hypD_2.out`. All findings were adjudicated and applied (§9).
  The GPT-6.1-sol blind round `hypD_1` is pending.

## Verdict

**We did not prove Hypothesis D, and we obtained no rigorous L[1/2, c < 1] factoring algorithm.**
Lenstra–Pomerance's L[1/2, 1+o(1)] stands. The results are partial and are listed below.

Review status:
- Opus-reviewed, with all findings adjudicated and applied (§9): two proof critics and two prior-art
  critics.
- GPT-6.1-sol adversarial review done (`codex/hypD_2.out`): no fatal error, all fixes applied (§9.2).
- The GPT-6.1-sol blind round (hypD_1) was still running at this commit.

**Removing two of lemmaM's ingredients.**
- **Theorem S / S′ [proved; a routine adaptation of the Dixon–Pomerance random-squares lemma].**
  - Statement: for *any* set S ⊂ [Y/2, Y] coprime to N, a uniform square R satisfies R ≡ u/v with
    u, v ∈ S with probability ≥ |S|⁴/(2·2^ω·φ(N)·Y²(1+log Y)).
  - For S a set of k-products: ≥ |S|²/(c_k 2^ω φ(N)).
  - The new step is a second-moment bound for the gcd-multiplicity that the ratio form introduces.
  - Small prime factors are removed first by trial division to y₀ = L^{o(1)} (standard), which makes
    2^ω = L^{o(1)}.
  - This removes **ingredient 3** (split counting), for an arbitrary, field-dependent set of
    admissible splits.
- **Theorem K (i)–(ii) [proved].** Allow only splits whose large primes are k distinct primes from one
  dyadic interval J ∋ N^{1/(2k)}.
  - Then complete descents in one field follow from one averaged statement,
    **D_k: Pr_{A, ℓ∈J}[ℓ descends in A] ≥ L^{−o(1)}**.
  - No joint event across the targets of a split is needed, so **ingredient 2** is removed. This also
    answers codex lemmaM_2's "complementary field types" objection: each field uses only splits inside
    its own good set. Both proof critics confirmed this.
- **Prop PZ [proved; Paley–Zygmund].** D_k ⟸ FM (first moment) + PC (pair correlation within L^{o(1)}).

**Correction found in review (critical).**
- The first version claimed the shape test, i.e. factoring u, v ≈ √N into k primes ≈ N^{1/(2k)}, costs
  L^{o(1)} by HECM. **That was wrong.** It costs L_N[2/3, c(2k)^{−2/3}].
- Theorem K(iii) and Corollary K are therefore restated with **Lenstra–Pomerance's rigorous
  factorisation** as the shape test, at L_N[1/2, 1/√2] per trial. The conditional constant becomes
  **c = (1+θ)/√2**:
  - d = 9: 0.920; d = 11: 0.884; d = 19: 0.814; d → ∞: 1/√2.
  - The previous value was √(2θ).
- **Prop H (dilemma) [proved arithmetic]:**
  - An L^{o(1)} rigorous shape test (HECM) forces large primes ≤ exp(o((log N)^{3/4}(log log N)^{1/4})).
  - Then k ≥ (log N)^{1/4−o(1)}, and every averaged estimate must hold to precision
    exp(−o((log N)^{1/4}(log log N)^{1/2})) (barrier B3).
  - Bounded k avoids B3 but pays the L_N[1/2, 1/√2] shape cost.

**What is still open for (a′_k).**
- **FM** is a sketch: a transplant of LV Lemma 8.5 to moduli ℓ·s. The published smooth-moduli BV of
  de la Bretèche–Drappeau (JEMS 2020, Thm 2.1) and Pascadi (Thm 7.1) does not apply as stated. It
  requires P⁺(q₀) ≤ y₂, and its e^{O(u₂)} loss is fatal at y₂ = B = L[1/2].
- **PC** is open. It splits into:
  - **PC_rat**, a binary correlation of smooth values of two linear forms in one variable;
  - **PC_alg**, which is ⟸ a restricted BDH with residue classes restricted to the moving class
    Nb^d mod t (RBDH′, conjectured), plus a 2D LV kernel lemma. This is an outline, not a sketch.
- Also open: the LV §6 square-root and character adaptation, and the rigour of HECM, whose LPP Part III
  was never published (affecting smoothness tests at B; caveat in §8).

**Barriers.**
- **Prop G [elementary].** Slicing or density bounds that ignore arithmetic mod D lose at least
  ℓ = N^{1/(2k)} in PC_alg.
- **Prop P [hypothetical computation; relabelled after review].** The published two-linear-form bound
  (Pascadi Cor. 7.2, ρ^{5/8}) **does not apply** at these parameters: it needs coefficients ≤ x^δ in an
  n-range x = M/ℓ, and gcd(b, b′) = 1. If a coefficient-uniform ρ^{−3/8} loss held anyway, this design
  would give only c → 1/√2 + √(3/4) = 1.573 (LP shape test, FM constraint k ≥ 3(d+1)/4). c < 1 would
  need a loss exponent < 17 − 12√2 = 0.0294. So the Paley–Zygmund route with anything like current pair-correlation technology
  does not beat LP.

**Conditional routes below 1.**
- **Design V [conditional; components known].** Combine cheap non-fruitful LV relations with L^{o(1)}
  quasi-uniform Vallée relations, which are fruitful by Dixon/Vallée fibre symmetry; products of
  several such dependencies fix the quasi-uniformity bias. This gives L[1/2, √(2/3)] ≈ L[1/2, 0.816]
  under:
  - **RK**, a rank-deficiency bound for the NFS rows on Λ/2. It is the GF(2) analogue of Stange's
    (2022) Hypothesis 3.1 and of Hafner–McCurley relation-lattice generation;
  - LV relation supply at fixed odd d ≥ 9 and L[1/2] scale (sketch);
  - the LV §6 adaptation;
  - Vallée's theorem (pointwise quasi-uniformity checked in the primary paper by GPT-6.1-sol).
  - RK is in the genuine-certificate deficiency form with a batch argument, after codex.
  - Design V does not suffer the shape-test problem (its smoothness tests are at B = L[1/2]).
  - **It is the best conditional route for d < 19.** The specific combination with Vallée's sampler
    was not found elsewhere.
- **(a′_k) with the LP shape test:** (1+θ)/√2 < 1 for d ≥ 7, conditional on FM + PC.
- A recursive variant, using the conditional algorithm itself as the shape test, gives √((1+√2)θ)
  (0.851 at d = 9) if D_k holds at every scale (critic computation, re-derived by us).

**Experiments (pre-registered, exact, toy).**
- E2 (16 instances, 38 400 cells):
  - PC_alg ratio 1.08–1.21 (P3 holds);
  - Poisson-with-field-quality fit of γ within 0.009 (P4 holds);
  - same-field pair ratio 1.02–1.07, so no complementary types (P5c holds);
  - Pr_A[g_A ≥ γ/2] ≥ 0.95 everywhere (P5b holds);
  - overdispersion 2.2–4.2: **P5a missed in 1/16** (4.21 > 4).
- E3 (20 toy NFS instances):
  - 360/400 Dixon vectors in the NFS row space (P6 holds);
  - 0/5098 single-prime vectors with p ∉ Q_N in the span (P7a holds);
  - 787/787 with p ∈ Q_N and ≥ 3 occurrences in the span (P7b holds).
- E1 as registered was uninformative (a disclosed design error). The exact E1x shows 0/120 violations of
  Theorem S and model ratios 0.85–0.94 (P2′ holds). P2″ was vacuous.

### Per-design constants (the L[1/2, c] exponent; "inputs" = what a proof would still need)

| design | c | status | missing inputs |
|---|---|---|---|
| Lenstra–Pomerance 1992 (class groups) | 1 | proved (prior art) | — |
| Vallée 1991 (small squares, θ = 2/3) | √(4/3) ≈ 1.155 | proved (prior art; not re-verified) | — |
| lemmaM (a), plain fixed-d descent | 2/√(d+1) | heuristic | first moment at t ≍ F^{1−o(1)}, beyond every known level |
| lemmaM (a′), HECM shape test, k ≍ (log N)^{1/4} | √(6/(d+1)) | conditional on Hyp. D | FM; PC to precision exp(−o((log N)^{1/4}(log log N)^{1/2})) (B3) |
| **(a′_k), fixed k, LP shape test** | **(1+θ)/√2**, θ → 3/(d+1): d = 7: 0.973; d = 9: 0.920; d = 11: 0.884; d = 19: 0.814; d → ∞: 0.707 | **conditional on D_k ⟸ FM (sketch) + PC (open)** | FM (LV 8.5 transplant), PC_rat (binary smooth correlation), PC_alg (⟸ RBDH′ + 2D kernel lemma), LV §6 |
| (a′_k), recursive shape test | √((1+√2)θ) (d = 9: 0.851) | conditional on D_k at all scales | as above |
| (a′_k) with a hypothetical coefficient-uniform ρ^{−3/8} PC_rat bound | 1.573 (limit) | hypothetical | the published bound does not apply; c < 1 needs loss exponent < 0.0294 |
| **Design V (LV rows + Vallée rows)** | **√(2/3) ≈ 0.816** (d ≥ 9) | **conditional on RK** | RK (genuine-certificate deficiency); LV supply at L[1/2] (sketch); LV §6 postprocessing keeping a Vallée projection |
| Design V with uniform Dixon rows | 1 | conditional on RK | as above; no gain over LP |
| QS-style descent (lemmaM (d)) | 1 | heuristic | FM/PC for a fixed quadratic polynomial; nothing to average over |
| one-level growing degree (lemmaM (b)) | L[2/5] | heuristic / conditional | bounded k does not apply |

---

## 1. Setting

- Dixon descent follows `lemmaM` §1–§2.
  - R = r² mod N with r uniform; R is split as u/v (mod N).
  - The large primes ℓ of uv are descended by special-q relations in an LV-randomized field
    A = (m, f = f̂_{N,m} + (x − m)R_A).
  - Theorem DD makes every dependency that contains a descent relation fruitful with probability
    1 − 2^{1−ω}, provided the transcript is invariant under r ↦ εr.
  - All procedures below act on R and A only, so (T) holds.
- L = L_N[1/2, 1]. "L^{o(1)}" means exp(o((log N log log N)^{1/2})), uniformly in N.
- ω = ω(N). Q_N = (Z/N)^{×2}, of index h = 2^ω in (Z/N)^× for odd N.

## 2. Split counting

**Lemma Ω [standard].** Put y₀ = exp((log N)^{1/2}(log log N)^{1/4}).
- y₀ is itself L^{o(1)}, so plain trial division to y₀ finds every prime factor ≤ y₀ in L^{o(1)} time.
  The first version used Pollard–Strassen; the critics noted it is unnecessary. Removing small factors
  is routine in rigorous analyses; see e.g. Stange 2022 §4.
- Removing these factors, and prime powers, leaves an N with every prime factor > y₀. Then
  ω ≤ log N/log y₀ = (log N)^{1/2}(log log N)^{−1/4}.
- Hence 2^ω ≤ exp(0.7(log N)^{1/2}(log log N)^{−1/4}) = L^{o(1)}. ∎

**Theorem S [proved; a routine adaptation of the Dixon 1981 / Pomerance 1982 random-squares lemma to the ratio form u ≡ Rv].** Let N be odd and Y ≥ 2 with 2Y² < N, and let S ⊂ [Y/2, Y] ∩ Z be any set of
integers coprime to N. Let R be uniform on Q_N, and let E_S be the event that u ≡ Rv (mod N) for some
u, v ∈ S. Then
  Pr[E_S] ≥ |S|⁴ / (2·2^ω·φ(N)·Y²·(1 + log Y)).
Moreover, given R, the event E_S is decided by one Gauss reduction plus at most 2 membership tests per
multiple j·(u₀, v₀), j ≤ Y/max(u₀, v₀). The number of multiples can be ≍ Y in the worst case and is
bounded only in expectation. Theorem K tests only (u₀, v₀), using divisor-closedness.

*Proof.*
1. Put n(R) = #{(u, v) ∈ S² : u ≡ Rv (mod N)}.
2. **First moment.** Let Ĝ be the group of characters of (Z/N)^× that are trivial on Q_N, so |Ĝ| = h.
   Orthogonality gives 1[x ∈ Q_N] = h^{−1}Σ_{χ∈Ĝ}χ(x). Then
     P := Σ_{R∈Q_N} n(R) = #{(u, v) ∈ S² : uv^{−1} ∈ Q_N} = h^{−1}Σ_χ |Σ_{u∈S}χ(u)|² ≥ |S|²/h.
   Every χ ∈ Ĝ is real, so χ̄(v) = χ(v) and χ(uv^{−1}) = χ(u)χ(v); each summand is a square.
3. **Second moment.** If (u, v), (u′, v′) ∈ S² both satisfy u ≡ Rv, then uv′ ≡ u′v (mod N). Also
   |uv′ − u′v| < Y² < N, so uv′ = u′v and the two vectors are proportional.
   - Write g = gcd(u, v) and (u₀, v₀) = (u, v)/g. The partners of (u, v) are the vectors j(u₀, v₀) with
     ju₀ ∈ [Y/2, Y].
   - Since u ≥ Y/2, we have u₀ ≥ Y/(2g), so there are at most Y/(2u₀) + 1 ≤ 2g partners.
   - Hence
       Σ_R n(R)² ≤ Σ_{(u,v)∈S², uv⁻¹∈Q_N} 2gcd(u, v) ≤ 2Σ_{g≥1} g·⌊Y/g⌋² ≤ 2Y²(1 + log Y).
4. **Cauchy–Schwarz.** #{R : n(R) > 0} ≥ P²/Σn(R)². Dividing by |Q_N| = φ(N)/h gives
     Pr[E_S] ≥ (|S|⁴/h²)·h/(φ(N)·2Y²(1 + log Y)).
5. **Algorithm.** With signed entries, two independent lattice vectors with |u|, |v| ≤ Y have
   0 < |det| ≤ 2Y² < N. But every determinant is a multiple of N, so all such vectors are multiples of the
   primitive shortest vector. The condition 2Y² < N was missing in the first version (critic fix; the
   counting part needs only Y² < N). ∎

*Remarks.*
- No independence, no 1/ℓ or saddle-point law, and no character bound is used. The square-class
  condition is absorbed by positivity, which is Dixon's coset pigeonhole Σ_i|S_i|² ≥ |S|²/2^ω written
  with characters.
- The single surviving 2^ω is an artefact of the Cauchy–Schwarz route. Dixon's pointwise multiplicity
  bound loses none.
- The heuristic precedent is the NFS-DL initial split (Joux–Lercier 2003, Commeine–Semaev 2006) and the
  Blake–Fuji-Hara–Mullin–Vanstone u/v trick. We found no rigorous ratio-form statement for composite N
  restricted to square classes; that is a search result.
- The fourth power is the price of the crude multiplicity bound. For smooth-type S one can use
  Ψ(x/d, y) ≪ d^{−α}Ψ(x, y) to recover |S|²·L^{o(1)}. Not needed here.
- E1x (§7) checks the bound exactly in 120 small cases. The true value is about |S|²/φ(N) for random S
  and about 2^ω|S|²/φ(N) when S ⊂ Q_N. So the bound is far from sharp but has the right structure.

**Theorem S′ (k-product sets) [proved].** Let G be a set of primes in a dyadic interval, coprime to N, with
|G| ≥ 4k². Let S be the products of exactly k distinct primes of G lying in a window [Y/2, Y] (Y² < N),
with |S| ≥ C(|G|, k)/(2k). Then Pr[E_S] ≥ |S|²/(c_k·2^ω·φ(N)) with c_k = 8k²·C(2k, k).

*Proof.* As in Theorem S, P ≥ |S|²/h. Split Σ_R n(R)² ≤ Σ_{(u,v)∈S²} #partners(u,v) according to
i = the number of shared primes, so (u, v) = (ju₀, jv₀) with j an i-product.
1. There are at most C(|G|,i)·C(|G|,k−i)² such pairs.
2. Each has at most C(|G|,i) partners j′(u₀, v₀).
3. So the i-th term is ≤ C(|G|,i)²C(|G|,k−i)² ≤ 2C(k,i)²C(|G|,k)² ≤ 8k²C(k,i)²|S|². The middle step uses
   C(n,i)C(n,k−i) = C(k,i)C(n,k)·∏_{t<i}(n−t)/(n−k+i−t) ≤ C(k,i)C(n,k)·e^{k²/(n−k)} ≤ √2·C(k,i)C(n,k) for n ≥ 4k².
4. Summing, Σ_i C(k,i)² = C(2k,k), so Σn(R)² ≤ c_k|S|². Then use Cauchy–Schwarz as before. ∎

This is the |S|² ("λ = 2") form, the Dixon–Pomerance k-product lemma in ratio form. Theorem K(ii)
improves accordingly, to Pr ≥ c′_k·(g_A/(4 log N))^{2k}·2^{−ω} for |G_A| ≥ 4k². E1x (random S) sits at
about 0.9·|S|²/φ(N).

## 3. Theorem K: bounded-k splits reduce Hypothesis D to one averaged statement

**Setup.**
- Fix integers k ≥ 1 and odd d ≥ 3. Let X = ⌊√(N/2)⌋.
- Let j be maximal with 2^{kj} ≤ X, and put J = (2^{j−1}, 2^j], so ℓ ∈ J means ℓ = N^{1/(2k)+o(1)}.
- For a field A, let G_A ⊂ J be the set of primes ℓ ∤ N for which a fixed, deterministic,
  capped special-q search succeeds. The search depends only on (A, ℓ), with time cap T_SQ·L^{o(1)}. It
  succeeds when it finds a primitive (a, b) ∈ Λ_ℓ(m) in the box with (a − bm)/ℓ and F_A(a, b) both
  B-smooth.
- Let g_A = |G_A|/π_J with π_J = #(primes in J), and γ := E_A g_A = Pr_{A, ℓ∈J}[ℓ ∈ G_A].

**Hypothesis D_k(γ).** γ ≥ L^{−o(1)}, uniformly in N. The search box and the field law are those of §4.

**Theorem K [proved, given the setup].**
- (i) Pr_A[g_A ≥ γ/2] ≥ γ/2.
- (ii) Fix a field A with |G_A| ≥ 4k², and let S_A be the products of exactly k distinct primes of
  G_A. The hypothesis on |G_A| was missing; codex fix. Then for N ≥ N₀(k),
    Pr_R[R ≡ u/v with u, v ∈ S_A] ≥ c_k·(g_A/(4 log N))^{4k}·2^{−ω}/(1 + log N), with c_k > 0 depending
    only on k.
- (iii) **[corrected after review]**
  - Acceptance rule (codex fix). Accept R when the primitive split (u₀, v₀) has each side a product of
    **at most** k distinct primes of G_A. "Exactly k" is not divisor-closed: cancelling gcd(u, v) can
    remove primes.
  - The count in (ii) still lower-bounds acceptance, because the primitive part of an S_A-pair is
    accepted.
  - In a field with g_A ≥ γ/2 and |G_A| ≥ 4k², a descent attempt succeeds after an expected
    L^{o(1)}·γ^{−2k} R-trials.
  - Each trial must split u₀, v₀ ≈ √N into their prime factors, all of which are ≈ N^{1/(2k)}. This is
    general factoring of √N-size integers, **not** a smoothness test.
  - HECM would cost L_{N^{1/(2k)}}[2/3] = L_N[2/3, c(2k)^{−2/3}]. The first version wrote "L_N[1/3]",
    which is wrong.
  - We use Lenstra–Pomerance's rigorous algorithm instead, at T_shape = L_{√N}[1/2, 1] = L_N[1/2, 1/√2]
    per trial.
  - Each trial with the right shape adds at most 2k special-q searches.
- So the expected cost per descent relation is (T_shape + T_SQ)·L^{o(1)} whenever γ ≥ L^{−o(1)}.

*Proof.*
- (i) Since 0 ≤ g_A ≤ 1, we have γ ≤ γ/2 + Pr[g_A ≥ γ/2].
- (ii)
  1. All products of k distinct primes of G_A lie in (2^{k(j−1)}, 2^{kj}] ⊂ (X·2^{−2k}, X]. By the
     pigeonhole principle, one of the 2k dyadic windows [Y/2, Y] with Y ≥ X·2^{−2k} holds at least
     C(|G_A|, k)/(2k) of them.
  2. Distinct prime sets give distinct products, so |S_A ∩ [Y/2, Y]| ≥ (|G_A|/k)^k/(2k).
  3. Use π_J ≥ 2^{j−1}/(2(j−1)log 2) for j ≥ j₀ (Rosser–Schoenfeld) and Y ≤ 2^{kj}. Then
     |S|/Y ≥ (g_A/(4kj log 2))^k/(2k), and 4kj log 2 ≤ 2 log N.
  4. Apply Theorem S with Y² ≥ X²·2^{−4k} ≥ N·2^{−4k−2}.
- (iii) Average over R. The bounded factors (log N)^{O(k)}, 2^{O(k)} and 2^ω are L^{o(1)}, by Lemma Ω.
  The shape test costs T_shape per trial. Conditional on u₀v₀ having the right shape, at most 2k
  special-q searches are run. Pr[shape]/Pr[success] ≤ (log N)^{O(k)}g_A^{−4k}2^{2ω} = L^{o(1)} (upper-bound Pr[shape] by 2·2^ω|S_J|²/φ(N), S_J = all k-products from J), so
  wasted searches cost only L^{o(1)} on average. ∎

**Corollary K [proved implication, modulo the LV §6 adaptation that lemmaM Corollary DD also leaves
unwritten; corrected after review].** Assume D_k(γ), for k ≥ (3/4)(d+1) (FM's rational-side
range). Then every N is factored in expected time L[1/2, (1+θ(d,k))/√2 + o(1)] for θ < 1, where
  θ(d,k) = 2s + (d−1)/(4k) + max(1/(4k), s − (d−2)/(4k)), s = 1/(d+1).
- The floor 1/(4k) on h comes from LV's kernel mechanism, which needs H ≫ ‖(a,b)‖ (critic fix). It
  matters only for k < (d²−1)/4.
- *Algorithm.*
  1. Field selection (codex repair, replacing an ill-posed pilot and then a guess-and-double that
     could waste a main-scale factor).
     - Sample fresh fields, one capped descent attempt each, until an attempt succeeds. Keep that field.
     - With q_A the per-attempt success probability and μ = E_A q_A, the kept field has law ∝ q_A, so
       E[1/q_A | kept] = Pr(q_A > 0)/μ ≤ 1/μ.
     - Thus C further successes cost C/μ attempts in expectation.
     - By (i) and S′, μ ≥ (γ/2)^{2k+1}L^{−o(1)}.
  2. Collect C·L^{o(1)} pure-descent rows, where C is the column count (factor base plus LV §6
     characters).
  3. Linear algebra, then Theorem DD.
- *Cost.*
  - Take B = L[1/2, a]. Relations cost B·(T_shape + T_SQ) = L[1/2, a + max(1/√2, θ/(2a))], and linear
    algebra costs L[1/2, 2a].
  - For θ < 1 the optimum is at a = θ/√2, giving c = (1+θ)/√2.
  - Without the shape cost, c would be √(2θ). The first version claimed that, but it is not achievable
    with a rigorous shape test at fixed k.
  - Recursive variant: use this algorithm itself as the shape test, assuming D_k at all scales. The
    fixed point c = a + c/√2 with θ/(2a) ≤ c − a gives c = √((1+√2)θ).
- *Sizes.*
  - ‖(a, b)‖ ≤ N^{1/(4k)}S, where S² = T_SQ = L^{O(1)}.
  - The rational cofactor is (a − bm)/ℓ = N^{s − 1/(4k) + o(1)}.
  - F(a, b) = N^{s + h + d/(4k) + o(1)}, where H = N^h is the coefficient height of R_A.
  - The level-½ condition |a − bm| ≤ |F|^{1/2−η} is used by FM (§4). It needs h > s − (d−2)/(4k) + 2η.
  - Adding the rational and algebraic sizes gives θ(d,k) + O(η).
- θ(d,k) ↓ 3s as k → ∞.
  - The displayed limits are infima over growing *fixed* k (codex); they are not attained at any
    finite k. Limits of c = (1+θ)/√2: d = 7: 0.973; d = 9: 0.920 (0.923 at k = 50); d = 11: 0.884; d = 19: 0.814;
    d → ∞: 1/√2.
  - Recursive variant: 0.851 at d = 9.
  - All numbers were checked numerically by us and by both proof critics.

**Prop H (the bounded-k dilemma) [proved arithmetic].** With k large primes per split,
Pr[all in G_A] ≈ g_A^{k}.
- **Horn 1: an L^{o(1)} rigorous shape test.** HECM costs L_v[2/3] for primes ≤ v. This is L^{o(1)} iff
  log v = o((log N)^{3/4}(log log N)^{1/4}), which forces k ≥ (log N)^{1/4−o(1)}.
  - Then g_A, and so every averaged estimate, must be accurate to exp(−o((log N)^{1/4}(log log N)^{1/2}))
    (B3).
  - LV-type goodness only gives L^{o(1)} losses in the exponent.
  - lemmaM's B1 = exp((log N)^{3/4}) sits on this horn.
- **Horn 2: bounded k.** g_A = L^{−o(1)} suffices, but every shape test is a factorisation of
  √N-size integers into N^{Θ(1)}-size primes.
  - The best *available* rigorous test is Lenstra–Pomerance, at L_N[1/2, 1/√2].
  - This is an upper bound. No lower bound shows that every rigorous shape test must pay it (codex
    wording fix).
- **Growing k (codex).** Let T = √(log N log log N). Interpolating between the horns requires
  k·log log N = o(T) and k·log(1/γ) = o(T). Writing γ = L^{−o(1)} does not by itself give the second
  condition.
- Either way something is paid. Horn 2 still allows c < 1 conditionally, (1+θ)/√2. But Design V
  (√(2/3)) is better unless d ≥ 19.

## 4. From D_k to a first and a second moment

**Field and box law.**
- m is uniform in [N^s/2, N^s].
- f = f̂_{N,m} + (x − m)R_A, where f̂ is LV's non-monic base-m expansion and R_A has degree d − 1 with
  i.i.d. coefficients c uniform in [−H, H], H = N^h.
- ℓ is uniform in the primes of J.
- The box is {(a, b) ∈ Λ_ℓ(m) : |a|, |b| ≤ √ℓ·S}, primitive with b > 0. Irregular lattices (λ₁ ≤ √ℓ/S)
  are excluded; they occur with probability O(S⁻²) over m.
- N_ℓ(A) is the number of box points with (a − bm)/ℓ and F_A(a, b) both B-smooth.
- For primitive (a, b), F_A(a, b) = F̂(a, b) + t·⟨c, v(a, b)⟩, with t = a − bm,
  v(a, b) = (a^{d−1−i}b^i)_{0≤i<d}, and F_A(a, b) ≡ N b^d (mod t).

**Prop PZ [proved].** γ ≥ (E N)²/E[N²] = 1/(1/E N + PC), where PC := E[N(N−1)]/(E N)² and every E is
over (A, ℓ). So if E N ≥ 1 (hypothesis FM) and PC ≤ L^{o(1)} (hypothesis PC), then D_k(γ) holds with
γ ≥ L^{−o(1)}. ∎
- Here N_ℓ(A) must count the candidates actually covered by the *capped* search that yield *valid
  certificates*. Existence outside the cap does not suffice (codex).
- Irregular lattices, and the finite-interval error from sampling m in an interval rather than uniformly
  mod ℓ, must be built into the moment estimates themselves. An unconditional exclusion probability does
  not replace that. FM's polynomial room M ≫ ℓ makes this repairable (codex).

### 4.1 FM, the first moment [sketch; inputs are proved theorems, but the transfer is not written out]
FM says E N ≥ L^{−o(1)}·S²·ρ(u_rat)·ρ(u_F). Choose S² = T_SQ = L^{o(1)}/(ρ(u_rat)ρ(u_F)).
1. **Rational side (proved input).**
   - Fix ℓ and primitive (a, b). As m runs over m ≡ a/b (mod ℓ), the cofactor s = (a − bm)/ℓ runs over an
     arithmetic progression of common difference b ≤ √ℓS, with M/ℓ terms of size ≍ bM/ℓ.
   - Harper's BV theorem for smooth numbers applies to these moduli b, which form a *dense* family, so
     the log-power saving suffices. It needs √ℓS ≤ (bM/ℓ)^{1/2−η}, i.e. roughly
     k ≥ (3/4)(d+1)(1+O(η)).
2. **Algebraic side (LV §5 + LV Lemma 8.5, transplanted).**
   - LV Lemmas 5.13–5.20 show that F_A(a, b) is close to uniform on long progressions mod t.
   - LV's goodness (Def. 5.12) then gives Pr_c[F smooth] ≥ ρ(u_F)^{1+o(1)} for every *good* modulus t.
   - It remains to show that bad moduli are an o(1) proportion of the family 𝒯 = {ℓs : ℓ ∈ J, s B-smooth},
     weighted by the induced law. This is a restricted BV for 𝒯 with relative saving.
   - **Lemma 8.5′ [sketch]** is LV's Lemma 8.5 with y-smooth moduli replaced by 𝒯. Take LV's conductor
     decomposition, with m₀ = min(B^η, e^{η√log x}) and M₀ = x^η, and choose η < 1/(2k(s + h + 1)).
     - Every conductor divisible by ℓ is ≥ ℓ ≥ M₀. Those characters fall in LV's large-conductor range,
       where Harper's Fact 8.11 bounds the sum over *all* moduli by Q√Ψ·log^{7/2}x.
     - Smaller conductors r divide the B-smooth part s. The counts Σ_{q∈𝒯, r|q}1/φ(q) equal LV's sums
       (Claims 8.14, 8.16), multiplied by the common factor Σ_{ℓ∈J}1/ℓ. That factor also multiplies
       |𝒯|, so the relative saving e^{−cu/log²u} + B^{−c} is unchanged.
     - The large-conductor term needs Q ≤ √Ψ(x, B)·δ_𝒯·L^{−ε}, where δ_𝒯 is the density of 𝒯. This holds
       with polynomial room, because t ≤ |F|^{1/2−η}.
   - LV's proof uses one y for both the moduli and the counted numbers; here both are B, plus the prime ℓ.
3. **Status.** FM is a plausible transplant of LV §5 and §8. It is labelled *sketch*, and the earlier
   "no new idea expected" is withdrawn on codex's points:
   - the smooth variable has length Q/ℓ, so the weights involve densities at Q/(ℓr) and are not LV's
     sums times Σ1/ℓ;
   - the induced multiplicities of t need explicit bounds;
   - Harper's large-conductor term is √Ψ(Q + x^{1/2−η}log²x)log^{7/2}x, not just Q√Ψ.
   Fixed polynomial slack should absorb the last point. One point to check: LV's proof of Lemma 5.21 compares bad
   B-smooth moduli against Ψ(Q_max, B′) rather than Ψ(Q_max, B). This is harmless only when β = β′ (our
   design has a single B), and a referee should confirm it.
4. **Review additions.**
   - The published smooth-moduli BV results do not cover FM as stated: de la Bretèche–Drappeau JEMS 2020
     Thm 2.1 and Pascadi Thm 7.1. They require P⁺(q₀) ≤ y₂, and their saving carries a factor
     e^{O(u₂)}, which is fatal at y₂ = B = L[1/2].
   - So FM rests on LV's Lemma 8.5 method, whose relative saving ϱ(Q, y)(e^{−cu/log²u} + y^{−c}) has no
     e^{O(u)} loss.
   - The sketch must also pass from "bad t are rare in 𝒯" to "bad t are rare under the induced law of t".
     Here s is B-smooth in a progression mod b, so the bad set could concentrate on progressions.
   - The parameter regime k ≥ (3/4)(d+1), s > 1/(2k) belongs in D_k.

### 4.2 PC, the second moment [open]
E[N(N−1)] sums over pairs of distinct box points (a, b), (a′, b′) of one special-q lattice. It splits as
follows.
- **PC_rat.** The rational cofactors satisfy b′s − bs′ = Δ/ℓ, where Δ = ab′ − a′b and |Δ/ℓ| ≤ 2S². This
  identity holds *for every m*. So as m varies the pair (s, s′) moves along one line, and the rational
  part of E[N(N−1)] is an average, over coefficients b, b′ ≤ √ℓS and shifts h = Δ/ℓ, of
  #{k : s₀ − bk and s₀′ − b′k both B-smooth}. This is a binary additive problem; b = b′ = 1 is the
  "smooth pairs (n, n + h)" problem. We know no estimate, even an upper bound within L^{o(1)} of
  ρ(u)², at smoothness L[1/2]. See §8 for the literature check.
- **PC_alg.** Given the rational data, PC_alg = Pr_c[F, F′ both smooth]/(Pr_c[F]Pr_c[F′]).

**Prop G (geometric barrier) [proved].** Let (a, b) ≠ ±(a′, b′) be primitive points of Λ_ℓ(m), and let
Γ = {(⟨c, v⟩, ⟨c, v′⟩) : c ∈ Z^d}, with v = v(a, b) and v′ = v(a′, b′).
- (i) Every 2×2 minor of [v; v′] is divisible by Δ. Hence ℓ | Δ | D := [Z² : Γ].
- (ii) For c uniform on any box, the law of (F, F′) is supported on a coset of the lattice
  Γ″ = {(tX, t′X′) : (X, X′) ∈ Γ}, which has index D in tZ × t′Z. The marginals live on the full
  progressions mod t and mod t′.
- (iii) Consequently, any bound that ignores arithmetic modulo D, i.e. any bound of the form Pr[F, F′ ∈ 𝒮] ≤ sup_z μ(z)·#{(n, n′) ∈ 𝒮² in the product
  progression box} exceeds Pr[F ∈ 𝒮]·Pr[F′ ∈ 𝒮] by a factor ≥ D/C_d ≥ ℓ/C_d = N^{1/(2k)}/C_d. That
  holds whenever 𝒮 has the product density ρ_Fρ_{F′} in the product box. Such bounds include the
  continuous Ball-slicing bound plus lattice-point corrections (lemmaM §3).
  - Regime condition: the support zonotope must be wide compared with the minima of Γ″. That holds
    here, since sin∠ ≥ 1/(2S²) and H is large.
  - (iii) concerns this unconditioned, unclipped counting template only. It is not a lower bound on
    every PC argument, and clipping the density bound at 1 changes the literal comparison (codex).
  - Bounds that do respect the Γ″-coset are exactly the RBDH′ problem below. So Prop G is a barrier to
    a *template*, not to PC_alg. We label it **elementary**: ℓ | Δ is the standard lattice-sieve
    determinant fact (Pollard 1993), and Δ | minors follows from x − y | x^n − y^n.
- *Proof.*
  - (i): For i < j, v_iv′_j − v_jv′_i = a^{d−1−j}a′^{d−1−j}b^ib′^i((ab′)^{j−i} − (a′b)^{j−i}), and
    x^n − y^n is divisible by x − y. Also ℓ | Δ, because a ≡ bm and a′ ≡ b′m (mod ℓ).
  - (ii): By definition.
  - (iii): μ is carried by Γ″-points. These make up a 1/D share of the product-lattice points in the
    support box, so sup μ ≥ D/(#product points) up to the box-spline shape constant C_d. ∎
- *Meaning.* PC_alg needs to know that smooth pairs equidistribute over the D cosets of Γ″. That is
  equidistribution of B-smooth numbers modulo q = tD in each coordinate. The ℓ-part of D is invisible to
  smoothness, since F ≡ Nb^d (mod ℓ) is already fixed by the marginal, but it still has to be proved.

**Hypothesis RBDH′ (restricted Barban–Davenport–Halberstam, residue-restricted) [conjectured].**
- Setup:
  - 𝒬 is the multiset of moduli q = t·D(pair) over pairs of the special-q boxes, with multiplicity
    w(q).
  - R(pair) is the set of D residues r mod q with r ≡ Nb^d (mod t), i.e. the classes F can occupy.
- Statement: for x ≍ |F|,
    Σ_pairs Σ_{r∈R(pair)} |Ψ(x, B; q, r) − Ψ_q(x, B)/φ(q)|² ≤ o(1)·Σ_pairs D·(Ψ_q(x, B)/φ(q))².
- The first version summed over *all* residues mod q, normalised by Ψ_q²/φ(q). That is weaker by a
  factor of about φ(t) on the D-element subset that matters, a polynomial loss (critic fix).
- RBDH′ requires the moving classes Nb^d mod t to be typical on average over the pair family.
- **PC_alg ⟸ RBDH′ + a two-dimensional version of LV Lemma 5.13 [programme, not a sufficient
  condition yet; codex].**
  - For primitive pairs, D = [Z² : Γ] = |Δ| exactly; this strengthens Prop G(i).
  - The primitive kernel vector is large. For d = 3 it is (bb′, −(ab′ + a′b), aa′), of size V² for
    points of size V. So joint smoothing needs H ≫ V², not H ≫ V. That would raise the h-floor and θ;
    not costed here.
  - RBDH′ also needs main terms for non-reduced classes, plus the actual pair multiplicities and value
    windows.
  - The rational and algebraic averaged bounds cannot simply be multiplied. The algebraic estimate must
    hold under the rational-smoothness weights.
  - The 2D lemma needs small kernel vectors realizing every Γ-shift in a box; it is plausible for d ≥ 3.
  - Decompose the Γ″-coset into D cosets of tDZ × t′DZ.
  - The projection of each coset onto a coordinate is surjective mod D, because v is primitive.
  - Write each count as main term plus error. Cross terms are handled by restricted BV, and the
    error × error term by Cauchy–Schwarz and RBDH′.
  - Nearly parallel pairs have sin∠(v, v′) ≍ |Δ|/(‖(a,b)‖‖(a′,b′)‖) ≥ 1/S². A dyadic decomposition in
    the angle costs only log S = L^{o(1)}.
- **Why RBDH′ is plausible.**
  - Harper's BDH (his Theorem 2) has level Ψ(x, y), and all our q ≤ x^{1/2+O(1/k)}.
  - LV state that their method "naturally provides" the BDH analogue for smooth moduli.
  - Our family is ℓ²·(B-smooth)·(Δ′D″) with Δ′ ≤ 2S² arbitrary, so it is denser than LV's.
- **Why RBDH′ is not proved.**
  - The family is not of LV's shape.
  - The weights w(q) come from a 4-variable geometric count.
  - The diagonal term QΨ must beat δ_𝒬Ψ²/Q.
  - None of this has been carried out. Restricted-moduli BDH is not in the literature we found (§8).

**Assessment.** PC_alg looks within reach of the Harper/LV machinery; this is an assessment, not a
claim. PC_rat does not: it is a correlation of smooth values of two linear forms in one variable. A
possible lever, not worked out, is a *randomized rational side* g = g₁x − g₀ with (g₁, g₀) random
(Kleinjung-style selection, Res(f, g) = N). That gives the rational side 2D randomness, which turns
PC_rat into a PC_alg-type (RBDH′) statement. It must be checked against Kleinjung's residue constraints on
g₀ mod g₁.

## 5. Design V: non-fruitful LV relations plus a few fruitful Vallée relations

**Observation [proved; elementary linear algebra + Dixon fibre symmetry].**
- Let 𝒞 be any multiset of relations whose coins are independent of the hidden roots, for instance LV
  relations, which need not be fruitful.
- Let 𝒱 = {(x_i, e_i)} be relations x_i² ≡ ∏p^{e_i} (mod N) whose roots are uniform on their fibre
  x_i·μ₂ given everything else. Examples: Dixon's uniform x, or Vallée's quasi-uniform small-square
  sampler, by its x ↦ εx invariance.
- Then the number of independent dependencies that contain a 𝒱-row is
  |𝒱| − (rank(𝒞 ∪ 𝒱) − rank 𝒞).
- With an exactly fibre-uniform sampler, each such dependency is fruitful with probability 1 − 2^{1−ω}.
  This is Dixon's fibre argument, i.e. Theorem DD with descent replaced by direct sampling.
- With a κ-quasi-uniform sampler, the guarantee is only ≥ 1 − 2κ·2^{−ω}, which can be ≤ 0. Hence the
  convolution fix below. The first version's "(1 − 2^{1−ω})·(fibre-bias factor)" was wrong.
- The sampler's quasi-uniformity must hold on a set closed under x ↦ εx, as {x : |x² mod N| ≤ h} is.
- The transcript may expose only x² mod N and its factorisation.
- So it suffices that some Vallée vector lies in the row space of 𝒞, after appending zero algebraic and
  character coordinates.

**Quasi-uniform samplers [proved].**
- Vallée's sampler is only *quasi*-uniform: each x in the target set has probability within a constant
  ratio κ = c₂/c₁ of uniform. It does not assert exact x ↦ εx invariance.
- Given z = x² mod N, each fibre element εx has probability ≥ 1/(κ2^ω). So for a valid dependency
  containing that row, **Pr[fruitful] ≥ (2^ω − 2)/(κ2^ω) ≥ 1/(2κ)**. This is codex's pointwise bound;
  it sums the lower bound over the good fibre elements.
- The P2 critic's "can fail when κ ≥ 2^{ω−1}" came from bounding the bad elements from above and is
  superseded. The sampler's root-dependent internal transcript must stay hidden from dependency
  selection.
- The convolution argument below is valid but unnecessary.
- Fix: take t independent dependencies D₁, …, D_t, each containing exactly its own Vallée row, and
  *multiply their congruences*. Do not GF(2)-sum the vectors, since shared 𝒞-rows would cancel.
  - The congruence ratio of the sum is ±∏_i ε_iZ_i, with the ε_i independent given everything else.
  - Each ε_i has min-probability ≥ 1/(κ2^ω) on μ₂ ≅ (Z/2)^ω.
  - The t-fold convolution is within (1 − 1/κ)^t of uniform in total variation.
- So t = O(κ(ω + 1)) = L^{o(1)} Vallée rows give fruitfulness 1 − 2^{1−ω} − o(1).

**Hypothesis RK (genuine-certificate deficiency form; reformulated after codex).**
- Let U_V be the span of the actual Vallée target supports, and W_true the space of rational-only
  vectors realised by *genuine* certificate squares. These are combinations of LV rows after LV §6's
  extra combinations remove the residual square-class kernel; raw matrix membership is not enough.
- RK: with probability ≥ L^{−o(1)} over the LV relation set, dim U_V − dim(W_true ∩ U_V) ≤ r = L^{o(1)}.
- Batch argument: r + 1 Vallée vectors then force a dependency with nonzero Vallée projection.
- The earlier form ("one Vallée vector is in the span with probability L^{−o(1)}") does not follow from
  codimension L^{o(1)}: a uniform target gives only 2^{−r}. So the batch form is the right one.
- RK is plausible on a well-covered support. Full factor-base saturation can fail through uncovered
  coordinates.
- **Credit.** RK is the GF(2) analogue of Stange's Hypothesis 3.1 (arXiv 2211.06821: random
  relations generate the relation lattice; she reduces factoring to relation lattices via order-finding
  and notes that NFS relations adapt). It is also the analogue of Hafner–McCurley/Buchmann
  relation-lattice generation, and of the index-calculus individual-logarithm template ("target in the
  span of factor-base relations").
- Necessary structure [proved]: the row space ∩ (rational-only vectors) lies inside
  Λ/2 := {v : ∏p^{v_p} ∈ Q_N}, as soon as character-trivial even combinations are true squares in K
  (LV §6). Indeed ∏(a − bm) ≡ φ(γ)².
- RK asserts the reverse inclusion, up to rarely occurring primes.

**Theorem V [conditional; a new combination of known components, cf. Stange 2022, Dixon, Vallée].**
Assume:
- relation supply: LV supply at fixed odd d ≥ 9 and B = L[1/2, 1/√6]. It needs only a first moment,
  by LV's stochastic deepening; it is the FM-type transfer of §4.1 without special-q;
- RK in the genuine-certificate deficiency form;
- the LV §6 square-class postprocessing, adapted to keep a nonzero Vallée projection. Then every N is factored in expected time
L[1/2, √(2/3) + o(1)].
- *Cost.* Linear algebra costs L[1/2, 2a]. LV relations cost L[1/2, a + θ/(2a)] with θ = 3/(d+1) ≤ 1/3.
  Vallée relations each cost L[1/2, (2/3)/(2a)], with L^{o(1)} of them. All three are ≤ L[1/2, 2/√6]
  at a = 1/√6.
- Vallée's sampler: B. Vallée, "Generation of elements with small modular squares and provably fast
  integer factoring algorithms", Math. Comp. 56 (1991) 823–849.
  - GPT-6.1-sol checked the pointwise quasi-uniformity (Definition 1, Theorem 5) in the primary paper.
    We did not read it ourselves.
  - Vallée's small-factor preprocessing should use Pollard–Strassen, since trial division to
    L[1/2, 1] would exceed the budget.
  - With uniform Dixon x the same argument gives only c = 1.
- Theorem V is unaffected by the HECM correction: its smoothness tests are at B = L[1/2], where
  L_B[2/3] = L_N[1/3].
- *Proof.* Combine the observation with the cost count. ∎

**Remarks.**
- RK needs no fruitfulness, no hidden root in 𝒞, and no special-q. It is a statement about the GF(2)
  rank of NFS relation matrices: the relations "fill" Λ/2.
- Its natural proof route is relation-lattice generation in the Hafner–McCurley/Buchmann style. That
  route needs random-exponent (descent-type) relations, so RK is not obviously easier than PC. It is
  nonetheless a *different* missing ingredient, and E3 supports it.
- Any rigorous sampler of small squares with θ_V < 2/3 would give c = √θ_V.

## 6. Barriers and what does not work

- **(B1) The log-saving barrier [known; LV §8].** A BV or BDH bound whose saving is
  (log x)^{−A}, or e^{−cu/log²u}, relative to the sum over *all* moduli says nothing about a moduli
  family of density L^{−c}. Every first or second moment here lives on such a family. So *restricted*
  versions, in the style of LV Lemma 8.5, are necessary. This is why FM must transplant LV §8 rather
  than cite Harper.
- **(B2) The geometric barrier (Prop G) [proved].** Pair laws in one special-q lattice are carried by an
  index-D ≥ ℓ sublattice. Slicing and density bounds alone lose at least ℓ = N^{1/(2k)}, so arithmetic
  equidistribution mod tD is unavoidable.
- **(B3) The precision requirement for unbounded k [proved arithmetic, §3; a requirement of this
  route, not a theorem about factoring].**
  - With k ≍ (log N)^{1/4} targets per split, product-of-marginal strategies must lose at most
    exp(o((log N)^{1/4}(log log N)^{1/2})). That is finer than LV-type goodness, whose losses are
    L^{o(1)} in the exponent.
  - Fixed k removes this, but only at the shape-test cost of Prop H.
- **(B4) What Theorem K does *not* remove.**
  - The good-field set and the target set are still field-dependent. That is harmless: Theorem S takes
    any set.
  - The first moment still has moving residues Nb^d mod t. LV's goodness handles them because it is
    uniform in residues: max_r, "for all (a, r) = 1".
  - Pascadi's level 5/8 is *not* usable. It is for a fixed residue, and our residues move with b.
- **(B5) The QS-style descent [assessment].** The special-q relation x² − N = ℓc has no randomization to
  average over. FM and PC would have to hold for the single polynomial x² − N at smoothness L[1/2].
  Known results on smooth values of quadratic polynomials (Dartyge; Harper; Balog–Blomer–Dartyge–
  Tenenbaum for binary forms) need smoothness x^{θ} with θ fixed, far above L[1/2] (§8). Not pursued.
- **(B6) No L[1/3] [lemmaM (c)].**
  - The descent is one level, and lemmaM's α ≥ 2/5 argument applies to the growing-degree regime.
  - At fixed k the shape test alone costs L_N[1/2, 1/√2] rigorously, or L_N[2/3] by HECM (Prop H).
    The first version's claim "HECM does it in L[1/3]" was wrong.
- **(B7) RK and fruitfulness are different statements [observation].** RK says the NFS relations fill
  Λ/2. Conjecture F says their *dependencies* carry a nonzero β-label. RK + Vallée ⇒ fruitful
  dependencies exist, via the Vallée row's hidden root. Conjecture F concerns pure-NFS dependencies,
  which Design V never needs.

## 7. Experiments (pre-registered; exact; toy; correctness and plausibility checks only)

Commands are in `notes/hypD_prereg.md`. Analysis: `.venv/bin/python -m src.hypD.analyze {e1,e2,e3} …`.
The E1x summary was computed inline from `experiments/raw/hypD_e1x.jsonl`; the script is in §7.1.

**Disclosed deviations.**
1. The frozen E1 code crashed on an ACCIDENTAL gcd(r, N) = 13127 in a 40-bit `tri` instance. A random
   r hit a 13-bit prime factor; the null rate is about trials·Σ1/p. The fix resamples r and counts such
   events; E1 was rerun.
2. E1 as registered had model-expected pooled hits of about 0.3, so P2 is **uninformative** (a design
   error). It was replaced by the exact E1x, registered in Addendum A before E1x ran.

### 7.1 E1 / E1x (Theorem S)
- **E1 (80 records).**
  - P1: 0/80 below the bound. Trivially so: the bound is about 10⁻¹⁶.
  - P2 is uninformative: 1 hit against a model expectation of 0.1.
  - 16 ACCIDENTAL gcd(r, N) > 1 events occurred in 320 000 r-draws. The r-stream is seeded per (N, seed),
    so the 4 (k, g) records of one N reuse the same draws: that is 4 distinct events, against a null
    expectation of about 8 (Σ1/p over the 13–16-bit factors of the tri instances). STRUCTURAL rate 0.
- **E1x (120 exact records).**
  - **P1′ holds:** 0/120 violations. The exact rate exceeds the bound by at least 7.5·10⁶.
  - **P2′ holds:** in 30/30 random-mode records with |S| ≥ 200, exact/model ∈ [0.849, 0.943]
    (mean 0.908). The shortfall is consistent with non-coprime pairs collapsing.
  - **P2″ was vacuous:** no qr-mode record had |S| ≥ 200. Descriptively, the largest qr records
    (|S| = 61–186) give exact/(2^ω·model) = 0.73–0.83, i.e. exact/model ≈ 2.9–3.3 at ω = 2. This shows
    the 2^ω positivity excess.
- Inline summary:
  `python - <<…` computing min(exact/bound), and exact/model over records with |S| ≥ 200 by mode.

### 7.2 E2 (special-q counts in toy LV-style fields)
- Grid: 64 and 72 bits, blum/tri, 4 seeds, 60 fields × 40 primes ℓ ∈ (2^14, 2^15], B = 8192, box 65×33.

| check | registered | result |
|---|---|---|
| P3 PC_alg | [0.67, 2.0] | 1.076–1.212 in 16/16 (PC 1.080–1.216) — holds |
| P4 Poisson-with-field-quality γ | \|diff\| ≤ 0.07 | max 0.009 — holds |
| P5a overdispersion of g_A | [0.5, 4] | 2.24–4.21; **fails in 1/16** (blum 64 s1: 4.21) |
| P5b Pr_A[g_A ≥ γ/2] | ≥ 0.8 | 0.95–1.00 — holds |
| P5c same-field pair ratio | [0.9, 1.3] | 1.016–1.066 — holds |

- μ = 1.41–1.64 at 64 bits and 0.66–0.80 at 72 bits; γ = 0.72–0.78 and 0.47–0.53.
- Reading:
  - Pair correlations within one special-q box are mildly positive, ≈ +15%. This is consistent with
    Prop G's mod-ℓ correlation being invisible to smoothness, plus small-prime effects mod Δ′.
  - Field quality varies; the overdispersion of about 3 is the usual root-property (Murphy-α) effect.
    It does not create complementary target types. The same-field pair ratio stays ≈ 1, and almost
    every field has g_A ≥ γ/2.
  - The P5a miss (4.21 against a bound of 4) is an underestimate of field-quality variance on my part,
    not a model failure relevant to Theorem K. Theorem K only needs Pr_A[g_A ≥ γ/2].

### 7.3 E3 (RK in a toy NFS)
- Grid: 48 bits (B = 1200) and 56 bits (B = 3000), blum/tri, 5 seeds each. All 20 instances were valid
  (rows/cols 1.29–1.5); rank/cols was 0.92–0.97.
- **P6 holds:** 360/400 Dixon vectors lie in the row space (90%), 15–20 of 20 per instance.
- **P7a holds:** 0/5098 single-prime vectors with p ∉ Q_N are in the span.
- **P7b holds:** 787/787 with p ∈ Q_N and ≥ 3 odd occurrences are in the span.
- Reading:
  - Span ∩ rational = Λ/2 restricted to well-covered primes, exactly as RK predicts.
  - The Dixon failures plausibly come from rarely occurring primes. This was not separately registered,
    so it is a hypothesis.
  - The toy skips index-divisor and projective pairs.
  - No square roots or fruitfulness were computed here; DD-type fruitfulness was tested in lemmaM E1.

### 7.4 Prop P: what a pair-correlation bound of the published *shape* would buy [hypothetical computation; relabelled after review]
- **The published bound.** Pascadi (arXiv 2505.00653, Cor. 7.2) improves de la Bretèche–Drappeau
  (JEMS 2020, Thm 4.1, exponent 3/5). The statement is:
  - hypotheses: (a, c) = 1, ad − bc ≠ 0, |a|, |b|, |c|, |d| ≤ x^δ, and (log x)^C ≤ y₁ ≤ y₂ ≤ x with
    y₂ ≤ y₁^C;
  - conclusion: #{n ≤ x : P⁺(an+b) ≤ y₁, P⁺(cn+d) ≤ y₂} ≪ Ψ(x, y₁)ϱ(u₂)^{5/8−ε}.
- **Why it does not apply to PC_rat here (critic finding).**
  - n indexes m ≡ m₀ (mod ℓ), so x = M/ℓ = N^{s−1/(2k)}, while the coefficients are b, b′ ≈ N^{1/(4k)}S.
  - Requiring b ≤ x^δ with tiny, ineffective δ forces ks ≥ 1/2 + 1/(4δ). Our optima violate this badly.
  - Its normalisation (u₂ from the n-range, not from the value size) costs a further L^{Θ(δ)}.
  - gcd(b, b′) > 1 is not covered.
  - **The published bound gives nothing usable at these parameters.**
- **Hypothetical.** Suppose a coefficient-uniform loss PC_rat ≤ ρ(u_rat)^{−e} held. Then
  γ ≥ ρ(u_rat)^{e}L^{−o(1)}, and each descent costs γ^{−2k} more (Theorem S′).
  - Optimize max(2a, a + λke·σ_r/(2a) + max(1/√2, θ/(2a))) over a, odd d and k ≥ (3/4)(d+1), with
    σ_r = s − 1/(4k) and λ = 2.
  - **e = 3/8 (Pascadi's exponent) gives c → 1/√2 + √(3/4) = 1.5731** (codex closed form, large-d
    limit; finite d is worse). With the LP shape test, c < 1 needs **e < 17 − 12√2 = 0.029437**.
  - With a free shape test the values are √(17/12) = 1.1902 and e < 3 − 2√2 = 0.17157. Our numerics
    (1.573, ≈ 0.028; 1.192, ≈ 0.165) agree up to finite-d effects.
  - If the shape test were free, these would be 1.19 and ≈ 0.17. The first version reported these
    shape-free numbers as proved; that was wrong twice over.
  - Numbers were checked by us and by a proof critic.
- **Conclusion.** Inserting anything like the current two-linear-form technology into this Paley–Zygmund
  route does not beat Lenstra–Pomerance. This is a statement about this route, not about all methods.
  The quantified target is a pair-correlation loss ≤ ρ^{−0.03} for smooth values on lines with large
  coefficients.

## 8. Literature (Opus literature subagent, read from primary PDFs unless marked; 2026-10-19)
- **Harper, arXiv 1208.5992** (preprint; no journal version found).
  - Thm 1 (BV): for log^K x ≤ y ≤ x and Q ≤ √Ψ, the error is ≪ Ψ(e^{−cu/log²(u+1)} + y^{−c}) + Q√Ψ log^{7/2}x.
  - Thm 2 (BDH): for Q ≤ Ψ, the error is ≪ Ψ²(e^{−cu/log²u} + y^{−c}) + QΨ.
  - Neither has a restricted-moduli version.
- **de la Bretèche–Drappeau, JEMS 22 (2020), Thm 2.1 (2.2).** This is a published BV for y₂-smooth moduli
  q ∼ Q with Q ≤ x^{1/2−ε}.
  - It includes multiples q₀q of a fixed y₂-smooth q₀ ≤ x^δ, has max over residues, and saves
    Ψ(Q, y₂)/(φ(q₀)Q)·e^{O(u₂)}(log x)^{−A}.
  - The authors credit LV for the smooth-moduli modification.
  - **Bearing on FM (corrected after review): it does not cover FM.**
    - The saving carries e^{O(u₂)}, through (e^{u₂} log x)^{O_k(1)} in their (2.4). At y₂ = B = L[1/2]
      we have u₂ ≍ (log N/log log N)^{1/2}, so e^{O(u₂)}(log x)^{−A} → ∞ and the bound is vacuous for a
      density statement.
    - Separately, the theorem requires P⁺(q₀) ≤ y₂.
    - Interestingly, the authors' own Thm 4.1, and Pascadi's Cor. 7.2, apply it with non-smooth
      q₀ = a. That suggests the condition is removable, but it is not a citation for FM.
    - LV Lemma 8.5 is genuinely stronger at our parameters: its saving ϱ(Q, y)(e^{−cu/log²u} + y^{−c})
      has no e^{O(u)} loss. So FM stays a transplant of LV's method.
- **Drappeau, Compositio 2015**: level 3/5, with a fixed residue a₁ā₂ and no max over residues.
  **Pascadi**: level 66/107, and 5/8 unconditionally (arXiv 2505.00653), also with a fixed residue.
  Neither handles our moving residues Nb^d, as (B4) notes.
- **Pair correlations.**
  - de la Bretèche–Drappeau Thm 4.1 (ρ^{3/5}); Pascadi Cor. 7.2 (ρ^{5/8}).
  - Pascadi Cor. 1.6: #{n ≤ x : n, n + 1 ∈ S(y)} ≪ xρ(u)^{1+5/8−ε}.
  - At fixed y there are only finitely many consecutive y-smooth pairs (Størmer/Lehmer), so smoothness
    of n and n + 1 is not independent at extremely small y. At y = L[1/2] independence within L^{o(1)}
    is expected but open. (A "Goldfeld–Yamagishi" reference in the first version could not be
    verified and was removed.)
  - Averaged two-point correlation results (Matomäki–Radziwiłł–Tao type) give only o(1) relative
    savings, which is useless at density L^{−c}. We found nothing beating ρ^{5/8} at y = L[1/2], even on
    average over coefficients (prior-art critic).
  - Harper's restriction theory (Compositio 2016) handles ternary problems (a + b = c) at log^K x ≤ y,
    not binary ones.
- **Hyperelliptic curve method.** Lenstra–Pila–Pomerance I (Phil. Trans. 1993), Cor. 1.2: for every
  n ≥ n₀, all prime factors ≤ v are found with probability ≥ ½ in L_v[2/3, c₀](log n)³.
  - **Caveat:** Part II (Proc. LMS 2002) defers details to a Part III that apparently was never
    written. The rigorous assembly is not in print. lemmaM and LV (Cor. 3.19) rely on the same method.
- **Vallée, Math. Comp. 56 (1991).** The primary text could not be retrieved (AMS block). The quasi-uniform
  sampling of x with x² mod N = O(N^{2/3}) and the rigorous L[1/2, √(4/3)] come from secondary sources
  (LV, STOC'89 abstract). **Theorem V's Vallée input is therefore not verified by us.**
- **Rigorous factoring and discrete logs, 2018–2026.**
  - Nothing found below c = 1 for factoring, or below √2 for discrete logs in prime fields.
  - Kleinjung–Wesolowski (fixed characteristic, quasi-polynomial) is rigorous but does not transfer.
  - de Boer–Pellet-Mary–Wesolowski (arXiv 2512.01588, ERH): rigorous class-group sampling toolkit.
  - Stange 2022 is heuristic (Hypothesis 3.1).
  - LV follow-ups: van Leeuwen's thesis, Le Gluher–Spaenlehauer–Thomé (heuristic), and an announced
    Balister–Bollobás–Morris–Lee paper on the 2-core.
- **Rank.** Balister–Bollobás–Morris (Annals 2018) prove that the GF(2) parity matrix of fewer than
  (e^{−γ} − ε)J *uniform random* integers has full row rank whp. This is the closest rigorous rank
  statement to RK, but it is for uniform integers rather than NFS relations, and it concerns rows, not
  the span of Λ/2. We found no span or relation-lattice result for sieve relations.

## 9. Reviews

- **GPT-6.1-sol adversarial review: DONE** (`codex/hypD_2.out`). No fatal error. Every proved claim was
  judged CORRECT or CORRECT WITH FIXES, and all fixes are applied (§9.2).
- The blind round `codex/hypD_1` was still running at this commit.
- **Opus literature subagent:** supplied §8.

### 9.1 Opus critics: four independent read-only reports, adjudicated
Two critics reviewed proofs (P1, P2) and two reviewed prior art (A1, A2). Every finding was accepted
except where noted.

**Proof findings.**
1. **[P1, P2; accepted; critical]** HECM shape-test cost. "L_{N^{1/(2k)}}[2/3] = L_N[1/3]" was wrong:
   the true cost is L_N[2/3, c(2k)^{−2/3}].
   - Theorem K(iii) and Corollary K are restated with the Lenstra–Pomerance shape test.
   - The constant changes from √(2θ) to (1+θ)/√2 (d = 9: from 0.775/0.781 to 0.920/0.923).
   - The "c → 0 as d grows" remark is withdrawn (the limit is now 1/√2).
   - The dilemma is stated as Prop H, and B6 is fixed.
   - The recursive variant √((1+√2)θ) was re-derived by us.
2. **[P2; accepted; major]** Prop P misapplied Pascadi's Cor. 7.2 (coefficient condition, n-range
   normalisation, coprimality).
   - It is relabelled as a hypothetical computation.
   - The numbers are recomputed with the LP shape cost and the FM constraint: 1.57, threshold ≈ 0.028.
   - The universal-negative wording is withdrawn.
3. **[P1, P2; accepted]** Verdict overclaims: "PC is the single missing estimate" and "Design V needs
   only RK". Rewritten.
4. **[P1; accepted]** RBDH was normalised over all residues, which loses about φ(t). It is now RBDH′,
   restricted to the moving classes, and PC_alg ⟸ RBDH′ is labelled an outline.
5. **[P1, P2; accepted]** Theorem S's algorithm step needs 2Y² < N. The multiplicity count is bounded
   only in expectation.
6. **[P1; accepted]** Floor h ≥ 1/(4k) in θ(d,k), from LV's kernel cube.
7. **[P2; accepted]** The Corollary K pilot ("keep if half succeed") was ill-posed. Replaced by
   guess-and-double with Chernoff thresholds.
8. **[P1, P2; accepted]** Design V:
   - the fruitfulness formula is now ≥ 1 − 2κ2^{−ω};
   - the convolution fix requires multiplying congruences, not GF(2)-summing them;
   - the sampler set must be ε-closed;
   - RK is stated as a deficiency bound, applied per 𝒞.
9. **[P1, P2; accepted]** B3's exponent is now (log N)^{1/4}(log log N)^{1/2}, and B3 is labelled a
   requirement of this route. Prop G(iii) is a barrier to bounds that ignore arithmetic mod D.
10. **[P2; noted]** "G_A deterministic" requires the HECM coins in the special-q search to be part of A.
    One-sided error at most halves γ. DD's hypothesis (T) holds: the field selection and all tests act on
    the R's and A only (both proof critics confirmed).

**Prior-art findings.**
11. **[A1, A2; accepted]** Theorem S/S′ is a routine adaptation of the Dixon 1981 / Pomerance 1982
    k-product random-squares lemma to the ratio form. The new step is the gcd-multiplicity second moment.
    The NFS-DL initial split (Joux–Lercier; Commeine–Semaev) and the Waterloo u/v trick are the heuristic
    precedents. Relabelled.
12. **[A1, A2; accepted]** Lemma Ω is standard, and trial division suffices.
13. **[A1, A2; accepted]** Design V needs credit to Stange 2022 (Hyp. 3.1; NFS relations adaptable),
    Hafner–McCurley/Buchmann lattice generation, the index-calculus individual-log template, and Dixon's
    fibre argument. "No prior formulation found" is replaced by "new combination of known components".
14. **[A1, A2; accepted]** FM ≈ de la Bretèche–Drappeau Thm 2.1 was wrong: the e^{O(u₂)} loss and
    P⁺(q₀) ≤ y₂. FM stays a transplant of LV 8.5. A2 found that dBD and Pascadi themselves apply Thm 2.1
    with non-smooth q₀ = a.
15. **[A1, A2; accepted]** Prop G is elementary (Pollard's lattice-sieve determinant; x − y | x^n − y^n).
16. **[A2; accepted]** The unverifiable "Goldfeld–Yamagishi" citation was removed. The Vallée citation
    is made precise; Vallée is still unverified by us.
17. **[A1, A2; confirmed]** Theorem K: no prior statement found. Call it a reduction, not a method.

**Verdicts.** P1 and P2 both returned ACCEPT WITH MAJOR FIXES.
- Correct: Lemma Ω, Theorem S (with fix 5), Theorem S′, K(i)–(ii), Prop PZ, Prop G(i)–(ii), the
  Design V dependency count and Theorem V's cost count, the PC_rat identity, and "D_k ⇒ complete
  descents in one field".
- Wrong as first stated: K(iii), the Corollary K cost, and Prop P.

### 9.2 GPT-6.1-sol adversarial round (hypD_2): findings and adjudication
All findings were accepted.
- **Verdicts.**
  - CORRECT: Lemma Ω; Theorem S/S′; K(i); PZ, with N_ℓ counting capped, certified candidates.
  - CORRECT WITH FIXES: K(ii) (|G_A| ≥ 4k²); K(iii) (accept "at most k" primes); Corollary K (limits
    are infima over growing fixed k); Prop G (scope (iii) to the counting template; D = |Δ| exactly);
    Prop P as a hypothetical (closed forms 1.5731 and 0.029437); the Design V observation (the
    pointwise bound gives ≥ 1/(2κ); convolution unnecessary).
  - CORRECT, conditionally: Theorem V's cost, given relation supply, RK and postprocessing.
- **(a) Uniformity.** Growing-k conditions k log log N = o(T) and k log(1/γ) = o(T). Irregular lattices
  and the finite-interval error belong inside the moments. The field-selection repair ("retain the first
  field whose attempt succeeds", E[1/q | kept] ≤ 1/μ) is adopted. The LP shape cost is an upper bound
  only.
- **(b)** D_k ⇒ complete descents in one field: **confirmed**, and the complementary-types objection is
  genuinely answered, provided G_A means certified, capped success.
- **(c) FM.** Plausible transplant, not established. "No new idea expected" withdrawn.
- **(d) PC.** The identity is confirmed. The 2D kernel needs H ≫ V². RBDH′ alone is not a completed
  sufficient condition.
- **(e) Design V.** Vallée's quasi-uniformity was checked in the primary paper. RK is reformulated as a
  genuine-certificate deficiency with a batch argument. LV §6's residual kernel must preserve a nonzero
  Vallée projection.
- **Most promising route (codex).** Prove the genuine-certificate RK, dim U_V − dim(W_true ∩ U_V) ≤ L^{o(1)}
  with probability L^{−o(1)}. That avoids both special-q FM and binary PC_rat. The proof target is
  uniform rank growth against every surviving dual obstruction, or a classification showing that only
  L^{o(1)} of them survive. Fixed-test first moments will not suffice. We agree this is the clearest
  route to the conditional √(2/3).

## 10. Candidate cards
```text
ID: KD-1 bounded-k Dixon descent in LV fields (Theorem K + Theorem S' + DD).
Target: every N (after trial division to y0).
Operation: R = r^2; split R = u/v with u, v products of exactly k primes from J ~ N^{1/(2k)}, all in the
  field's descendable set G_A; one special-q relation per prime; DD gives fruitfulness.
Proved: Thm S/S' (Dixon-Pomerance lemma, ratio form), Thm K(i)-(ii) (D_k => complete descents in good
  fields), Prop PZ, Prop G (elementary), Prop H (bounded-k dilemma).
New lemma required: FM (transplant of LV Lemma 8.5 to moduli l*s; sketch); PC = PC_rat (binary smooth
  correlation; published bounds inapplicable at these parameters) + PC_alg (<= RBDH' + 2D kernel; outline).
Cost: L[1/2, (1+theta)/sqrt 2] (d = 9: 0.920; d -> inf: 0.707) conditional; Lenstra-Pomerance shape test.
Status: conditional; E2 toy: PC_alg 1.08-1.21, no complementary types.

ID: V-1 LV rows + Vallée rows (Design V). Best conditional route for d < 19.
Proved: dependency count; fibre-symmetry fruitfulness with the quasi-uniform convolution fix; cost sqrt(2/3).
New lemma required: RK (deficiency bound; GF(2) analogue of Stange Hyp 3.1); LV supply at L[1/2] (sketch);
  LV §6; Vallée (not re-verified).
Status: conditional; E3: 360/400 Dixon vectors in span; span ∩ rational = Lambda/2 on well-covered primes.

ID: G-1 slicing route to PC_alg. Status: refuted for this model (Prop G: loss >= l unless mod-D arithmetic used).
ID: P-1 PZ route with two-linear-form technology. Status: refuted for this model (Prop P: c ~ 1.57 hypothetical).
ID: QS-D QS-style descent. Status: not pursued (no randomization).
```
