# cov38b: can Bourgain–Garaev-type bounds for inverses of band primes push the deterministic worst case below N^{1/5}?

Agent `cov38b`, phase 30, 2026-10-28 (conn_alg shortlist #2).
- Pre-registration: `notes/cov38b_prereg.md` (commit ed36d22; addendum A in ed46052, both before any run).
- Code: `src/cov38b/{design,fourier,experiments,analyze}.py` (harness). Tests: `tests/test_cov38b.py` (13 pass).
- Raw data: `experiments/raw/cov38b_e{1,2,3}.jsonl`; verdicts: `experiments/analysis/cov38b_summary.json`.
- GPT-6.1-sol: blind round `codex/cov38b_1.{md,out}`; adversarial round `codex/cov38b_2.{md,out}` (§8).
- Opus critics (read-only): one on proofs, one on literature (§8).
- Notation as in notes/coverage38_results.md and notes/dchworst_results.md. ‖x‖_q := |x mods q|.

## 0. Verdict (plain language)

**No. None of the Bourgain–Garaev / Kloosterman-sum estimates we checked reaches the needed exponent, and the worst case
stays at Harvey–Hittmeir's N^{1/5}. No new deterministic record.** For the standard way of using such estimates
(phase-discarding Erdős–Turán certificates) the failure is a theorem with an explicit exponent gap.

- **There is no intermediate prize.** Inside the design family of Theorem U (cost model with WIN), any coverage exponent c ≤ 3/8 gives
  exactly N^{1/5} (Theorem W: min(1/5, (1−c)/(7/2−c))), and even a logarithmic gain needs the full power statement
  λ < J^{1/4}. So "some c < 3/8 with an explicit δ" does not exist in this framework (Prop. 0).
- **The target (§1, Prop. T).** For every prime q ≍ P and every prime p ∈ (q/2, q), some prime ℓ ≡ 3 (mod 4) in
  (Y, 2Y], Y = q^{2/(7−2c)}, and some B ≤ H/4 satisfy ‖pBℓ̄‖_q ≤ qH/(4Y), with H = Y^{1/2}J^{−c}, c > 3/8. At
  c = 3/8: Y = q^{8/25}, H = q^{1/25}, the arcs have length q^{18/25}, and the expected number of hits is
  λ = J^{1/4} ≈ q^{0.08}. Uniform form: the J dilates ℓ·S₀ of one Farey-arc set S₀ of measure ≈ q^{−6/25} cover F_q^*
  (i.e. F_q^* ⊆ ∪ℓ·S₀).
- **Every phase-discarding (term-by-term) certificate fails, with a sharp threshold (Theorem B, proved).**
  Take a pointwise bound max_{h≠0}|Σ_ℓ e(hℓ̄/q)| ≤ J q^{−δ}, used through an Erdős–Turán-type certificate that bounds
  the error term by term in the frequency h. Then:
  - it certifies the uniform target for δ > δ*(c) + ε, and needs δ ≥ δ*(c) − o(1), with δ*(c) = (1+2c)/(7−2c): a
    sharp leading-exponent threshold. This holds for every certificate function that is ≤ 0 off the outer Farey set
    (any sign, scale or mass there). At c = 3/8, δ* = 7/25 = 0.28. The statement is asymptotic, for every fixed
    c > 0: the hidden loss is exp(O(log q/log log q)), and log-scale improvements over Thue are not excluded.
  - but Parseval caps every such bound at δ ≤ 1/(7−2c) = 4/25 = 0.16 (square-root cancellation).
  - **The gap is 3/25 in the exponent of q, beyond square-root cancellation.** Known results at length q^{0.32}
    give δ = 0 for primes (all prime-variable Kloosterman bounds start at q^{1/2+ε}) and only a (log q)^{1/2−o(1)} saving for integers (Bourgain–Garaev 2014 Thm 16).
  - New ingredient: the Fourier ℓ¹ norm of a Farey-arc set with H′ denominators is ≥ H′/polylog
    (Lemma F, via an exact Freiman isomorphism to rationals and an lcm-sum energy bound). Every nonzero Fourier
    coefficient of the set is ≤ K q^{o(1)} (Lemma F2, from GPT-6.1-sol's blind round). So the union over the H′ dilates
    buys nothing in any phase-discarding argument: one arc does as well.
  - Multiplicative characters (GRH included), joint phase discarding over (B, ℓ), and a-adaptive magnitudes with a
    fixed weight fail in the same way (Lemma M, Remark B2, Prop. B′).
  - **Not excluded:** estimates that keep the phases across h or B, e.g. a bilinear bound for
    Σ_h c_h Σ_B Σ_ℓ e(hBaℓ̄/q) with the prescribed Farey-arc coefficients (coverage38 §3.4), and weights depending
    on ℓ or on a. Phase-retaining bilinear estimates do exist at this length, but none checked here is strong enough.
- **Bourgain–Garaev specifically (§2).** Their exponential-sum bounds, used by phase discarding, fall under
  Theorem B whatever their strength. Their phase-retaining bilinear bounds (and Kerr–Shparlinski–Wu–Xi's), applied to
  the hit count, miss positivity by a power (error q^{7/25} against main term q^{2/25}).
  Their additive-energy input is already optimal here: whenever 16Y³ < q (true at c = 3/8), the inverses of band primes have the minimal
  possible additive energy 2J² − J (an exact identity). It still gives nothing, because the obstruction is pointwise
  in the dilate p, not a lack of energy control. The covering theorems for sums of reciprocals (Shparlinski 2002;
  Glibichuk 2006; not Bourgain–Garaev) concern additive sumsets, while the target is a covering by multiplicative
  dilates. Bourgain–Garaev 2014 contains no positivity statement for product sets or dilates.
- **Two-pronged algorithms (§4, Theorem TP, proved in the model).** "DCL coverage OR Lehman/Harvey windows" is itself a
  design of Theorem U's family. At budget N^{1/5−η}, a one-modulus Lehman prong covers at most a fraction ≍ N^{−5η/4} of pairs
  (E3 measures the predicted G^{−5/4} law: slopes −1.27 to −1.31). So the DCL prong must still cover every
  Lehman-hard pair at λ < J^{1/4}. Barrier B3 is genuinely escaped (linear families are Lehman-easy). In our assessment the residual
  statement is no easier: it has the same λ-threshold, and for coprime-pair proofs it still contains a sparse survivor (x, ⌊αx⌋).
- **Experiments (pre-registered; all 10 named predictions pass).**
  - E1: exhaustive worst case over all primes p ∈ (q/2, q) for 15 primes q up to 1.7·10⁹ (up to 40.6M pairs per q) at the c = 3/8
    scaling. The worst height h* lies in the independent-model band in 15/15 cases, and below the Thue height for
    every q ≥ 10⁷. The empirical exponent is c_emp = 0.11–0.14 at q ≈ 10⁸–10⁹, rising slowly, as the model predicts
    (the model's limit is 1/2).
  - **At these sizes the target itself is false at the prescribed height.** With H = ⌊√Y·J^{−3/8}⌋ = 4, 5, 6 at
    q ≈ 10⁷, 10⁸, 10⁹, there are 81 869, 299 426 and 2 584 494 uncovered prime pairs (GPT-6.1-sol's adversarial
    recount of our histograms). c > 3/8 is an asymptotic conjecture that toy sizes cannot test directly.
  - E2: the Fourier ℓ¹ norm of Farey-arc minorants grows like H′^{0.80–0.85} (Lemma F's H′/polylog), and
    phase-discarding certificates appear only at H ≈ 9–14·√Y, against a true worst case at or below √Y.
- **Bottom line.** A sub-1/5 proof by this route needs a positivity mechanism for unions over moduli that does not
  discard phases. The candidates are: a phase-retaining multilinear estimate; a Poisson-quality centred moment over
  the dilates of order ≥ 26; or some joint geometry-of-numbers or arithmetic argument. None strong enough is known,
  even for a log-factor gain over Thue/Minkowski. Our assessment is that it is out of reach with current technology. This is an assessment, not a
  theorem: it is not a lower bound for factoring.

## 1. The target statement (Task 1)

**Setting.** Fix c ∈ (3/8, 1/2]. For n-bit N let P = 2^{n/2}, let Y = Y(n) be n-only with
Y = N^{1/(7−2c)+o(1)} = P^{2/(7−2c)+o(1)}, let 𝒫 be the J ≍ Y/log Y primes ℓ ≡ 3 (mod 4) in (Y, 2Y], and let
H = ⌊Y^{1/2}J^{−c}⌋. Write q for the larger prime of a pair (p, q) ∈ S_n, so q ∈ (P/√2, √2P) and p ∈ [q/2, q).

- **Cov(c)** (the pair form; Theorem W's hypothesis, paper Cor. (b)). Every (p, q) ∈ S_n has ℓ ∈ 𝒫 and a primitive
  (A, B), 1 ≤ |A|, B ≤ H, with ℓ | Bp − Aq.
- **Pos(c)** (reciprocal form). For every pair: ∃ ℓ ∈ 𝒫, 1 ≤ B ≤ H/4 with ‖pBℓ̄‖_q ≤ k₀ := ⌊qH/(4Y)⌋.
- **Pos∀(c)** (uniform form). For every prime q in range and every a ∈ F_q^*: ∃ ℓ, B as above with ‖aBℓ̄‖_q ≤ k₀.
  Equivalently **F_q^* ⊆ ∪_{ℓ∈𝒫} ℓ·S₀**, where S₀ = S(H/4, k₀) and S(H′, k₀) := {B̄k : 1 ≤ B ≤ H′, |k| ≤ k₀} is the
  Farey-arc set. Equivalently F_q^* ⊆ 𝒫·[−k₀, k₀]·[1, H/4]^{−1}.

**Proposition T (proved).**
- (i) Pos(c) ⇒ Cov(c), and Pos∀(c) ⇒ Pos(c).
- (ii) Cov(c) at height H ⇒ the reciprocal statement with (H, 2qH/Y) in place of (H/4, qH/(4Y)). So Cov and Pos are
  equivalent up to a constant factor in H, which changes c by o(1).
- (iii) **Size.** |S(H′, k₀)| = (12/π² + o(1))H′k₀ when 2H′k₀ < q (injectivity), H′ → ∞ and k₀/H′ → ∞
  (for H′ = 1 it is 2k₀ + 1). The expected number of hits is
  λ = J|S₀|/q ≍ JH²/Y ≍ J^{1−2c}. At c = 3/8: Y = q^{8/25}, H = q^{1/25}, k₀ = q^{18/25}, the measure of S₀ is
  ≈ q^{−6/25}, and λ ≈ J^{1/4} ≈ q^{2/25}.
- (iv) Cov(c), c > 3/8 ⇒ deterministic worst case N^{(1−c)/(7/2−c)+o(1)} < N^{1/5} on S_n (Theorems R, W).

*Proof.*
- (i) Let k ≡ pBℓ̄ (mod q) with |k| ≤ k₀, and put A := (pB − ℓk)/q ∈ Z.
  - |A| ≤ (pB + ℓ|k|)/q < B + 2Y·k₀/q ≤ H/4 + H/2 < H.
  - A ≠ 0: otherwise pB = ℓk with ℓ ∤ B (B < ℓ) and ℓ ≠ p.
  - Then ℓ | Bp − Aq = ℓk. Divide (A, B) by g = gcd(A, B): ℓ ∤ g since g ≤ H < ℓ, so the primitive weight still
    works.
  - (`test_reciprocal_test_implies_coverage` checks this exactly.)
- (ii) Write Bp − Aq = ℓk. Then pB ≡ ℓk (mod q), so pBℓ̄ ≡ k, and |k| ≤ (Bp + |A|q)/ℓ ≤ 2Hq/Y.
- (iii) The rationals k/B are distinct mod q (Lemma F0 below). Count coprime pairs.
- (iv) Paper Corollary (b). ∎

**Proposition 0 (what weaker or averaged forms Theorem R can use; proved in the model, from Theorems U and W).**
- (a) Only prime a = p ∈ [q/2, q) is needed, and only q in the n-bit balanced range. Pos∀ is a convenient sufficient
  form, not a necessity.
- (b) **No intermediate exponent.** For every design of the coverage38 family (N-uniform groups, per-group charge,
  WIN), X ≥ cN^{1/5}(λ⁴/J)^{1/5} (Theorem U). Coverage of every pair at λ = J^θ with θ ≥ 1/4 therefore costs
  ≥ cN^{1/5}. Consistently, Theorem W's (achievable) exponent is min(1/5, (1−c)/(7/2−c)), which equals 1/5 for all c ≤ 3/8.
  - So a weaker positivity statement (c ≤ 3/8) gives **no** N^{1/5−δ} through Theorem R.
  - Log gains: λ = J^{1/4}(log N)^{−a} buys at most (log N)^{4a/5} (U4), and it is still a power-strength statement.
- (c) **Averaged forms help only if they make the exceptional set empty.** Theorem R guarantees success only on
  covered pairs. An almost-all statement (all but E pairs, E ≥ 1) needs a second algorithm for E.
  - If that algorithm is in the family, Theorem U applies to the union design (§4), and the union must cover every
    pair at combined λ < J^{1/4}.
  - Outside the family, nothing is known that recognises or factors the uncovered set. inversesieve shows it carries
    no usable local structure, and forgetful containers fail.
  - **An averaged estimate strong enough to force E = ∅ does suffice.** GPT-6.1-sol's blind round corrected our
    draft's "no averaged form helps".
    - Let C_q(a) = Σ_ℓ 1_{S₀}(aℓ̄), with exact nonzero-residue mean Λ_q = J(|S₀| − 1)/(q − 1).
    - Per q, a centred-moment bound Σ_{a≠0}|C_q(a) − Λ_q|^{2r} ≪ q^{1+o(1)}Λ_q^r over the dilates a would imply
      positivity for every a once Λ_q^r > q^{1+ε}. A single zero contributes Λ_q^{2r}.
    - At c = 3/8, Λ_q = q^{2/25}, so r = 13 (a **26th** moment) suffices. In general r > (7−2c)/(2(1−2c)), for c < 1/2.
    - Over the pair population instead (sieve-style), Chebyshev with |S_n| ≈ q² pairs needs an even moment order
      > 2(7−2c)/(1−2c) = 50, i.e. ≥ 52. (This corrects "25" in v1, after both reviewers.)
  - **Literature status: we know of no Poisson-quality per-q centred moment over a, even of order 4.**
    - Heuristic picture, for a smoothed interval: the 2r-th moment involves 2r points u_i = B_iℓ̄_i, and its error
      comes from relations Σm_iu_i ≡ 0 with small m_i (|m_i| ≲ q/k₀ = Y/H).
    - After clearing denominators these are congruences mod q between integers of size ≈ Y^{2r}. They exceed q already
      for 4 points, so a simple rational lift fails.
    - This shows the easy route fails. It does not show that no fourth-moment theorem exists.
  - Different moments are proved elsewhere in the project. They are moments over the pair population (D′'s large
    sieve, sieve2d's two-variable fourth moment), with tails λ^{−1}, λ^{−2}. The pair-box level (≤ |S_n| ≈ Y^{7−2c})
    supports such sieve moments of order ≤ 6 at best, and sieve2d's error analysis stops at 4–5. Truncated
    inclusion–exclusion needs depth ≳ λ.

## 2. Tool exponents (Task 2)

### 2.1 The phase-discarding (Erdős–Turán) certificate and its exact threshold

Fix q, the band 𝒫, J = |𝒫|, and K(h) = Σ_{ℓ∈𝒫} e(hℓ̄/q). A **certificate** is a real g on Z/q plus a number κ ≥ max_{h≠0}|K(h)|. Since
Σ_ℓ g(aℓ̄) = (J/q)ĝ(0) + (1/q)Σ_{h≠0} ĝ(h)K(ah), it certifies Σ_ℓ g(aℓ̄) > 0 for every a once

  J·Σg/q > κ·W(g),   W(g) := (1/q)Σ_{h≠0}|ĝ(h)|.     (★)

**Scope (corrected after GPT-6.1-sol's adversarial round).** Theorem B concerns certificates with g ≤ 0 off the
outer set S_out := S(H, 2qH/Y).
- Every residue x at which a prime input p ∈ [q/2, q) can be hit (x = pℓ̄) with the hit implying coverage via ℓ lies in
  S_out, by Prop. T(ii).
- A certificate that only has to work on prime inputs could add residues x such that xℓ mod q is never a prime in
  [q/2, q). Such certificates are **not** covered by Theorem B; nothing is known about that set.
- The draft's "valid for every a" definition was wrong. For a = ℓ, the hit aℓ̄ = 1 ∈ S₀ gives A = 0 and no coverage,
  so even S₀ is not "valid for every a". Pos∀ is a statement about residues, and it implies coverage only for prime a.

For c > 1/4 and Y large, 8H⁴ < Y, which is the Freiman condition below. Lemma F bounds
the full ℓ¹ norm W_full = W + Σg/q; since Σg/q ≤ |S_out|/q = o(1) and W_full ≥ sup g, we have W ≥ W_full/2.

**Lemma P (Parseval; proved).** For J distinct points, max_{h≠0}|K(h)| ≥ √(J(q−J)/(q−1)). So in (★)
κ ≥ √(J/2) for J ≤ q/2. Writing κ = J q^{−δ} with J = q^{α}, α = 2/(7−2c): **δ ≤ α/2 = 1/(7−2c)** for every valid κ.

**Lemma F0 (Freiman isomorphism; proved).** Let 4H′³k₀ < q and R = {k/B ∈ Q : |k| ≤ k₀, 1 ≤ B ≤ H′}. The map
k/B ↦ B̄k mod q is a bijection R → S(H′, k₀). For r₁, …, r₄ ∈ R, r₁ + r₂ ≡ r₃ + r₄ (mod q) iff r₁ + r₂ = r₃ + r₄ in Q.
- *Proof.* Clear denominators by D = B₁B₂B₃B₄. The numerator of r₁ + r₂ − r₃ − r₄ has absolute value
  ≤ 4k₀H′³ < q, and it is ≡ 0 (mod q) iff it is 0. Injectivity is the case of two terms. ∎
- In particular E(S) = E_Q(R). `test_freiman_energy_identity` checks this exactly.

**Lemma F1 (energy of Farey-arc sets; proved).** If moreover k₀ ≥ H′³, then

  E(S(H′, k₀)) ≤ 8(2k₀+1)³·Σ_{B₁..B₄≤H′} min_i B_i / lcm(B₁..B₄) ≤ C k₀³ H′ (log 2H′)^{16}.

Also |S| ≥ k₀H′/2 for H′ ≥ H₀.
- *Counting step.* Count over representations, which is an upper bound. For fixed denominators let L = lcm, c_i = L/B_i.
  The equation is Σ±c_ik_i = 0 with gcd(c) = 1.
  - Let i be the index with the largest c_i = L/B_min ≤ H′³ ≤ 2k₀ + 1.
  - The other three k's lie in a sublattice of Z³ of index c_i that contains c_iZ³. So a box of side n = 2k₀+1 holds
    ≤ (n + c_i)³/c_i of them, and k_i is then determined. That gives ≤ 8n³B_min/L.
- *Divisor step.* Use min ≤ B₁ and B₁/lcm(B₁, M) = gcd(B₁, M)/M, with M = lcm(B₂, B₃, B₄).
  - Σ_{B₁≤H′} gcd(B₁, M) ≤ H′τ(M).
  - #{(B₂, B₃, B₄) : lcm = M} ≤ τ(M)³.
  - So the sum is ≤ H′Σ_{M≤H′³}τ(M)⁴/M ≪ H′(log H′)^{16}.
- `test_lemma_F1_energy_bound` and `test_lcm_sum_divisor_bound` check both inequalities. ∎
- F1 is sharp up to the polylog. The sets A*_B = {k/B : gcd(k, B) = 1, |k| ≤ k₀} are disjoint, and Cauchy–Schwarz
  gives E(A*_B) ≥ |A*_B|⁴/(4k₀+1) ≍ k₀³(φ(B)/B)⁴. Summing over B gives E(S) ≥ ck₀³H′. (Corrected after the
  critic: the earlier "single-denominator progressions" argument over-counted the shared integers.)

**Lemma F (Fourier ℓ¹ norm of Farey-arc minorants; proved).** Let S = S(H′, k₀) with 4H′³k₀ < q and k₀ ≥ H′³.
- **(a) Nonnegative, θ-mass.** If 0 ≤ g ≤ 1_S and Σg ≥ θ|S|, then (1/q)Σ_h|ĝ(h)| ≥ cθ³H′/(log 2H′)^8.
- **(b) Any sign.** If g is real, g ≤ 0 off S, sup g = 1 and Σg = θ|S| > 0, then
  (1/q)Σ_h|ĝ(h)| ≥ max(1, cθ^{4/3}H′^{2/3}/(log 2H′)^{16/3}). Consequently Σg / ((1/q)Σ_h|ĝ|) ≤ C k₀√H′(log 2H′)^4.

*Proof.*
- (a) Hölder gives Σ|ĝ| ≥ (Σ|ĝ|²)^{3/2}/(Σ|ĝ|⁴)^{1/2} = q(Σg²)^{3/2}/E(g)^{1/2}. Then:
  - Σg² ≥ (Σg)²/|S| ≥ θ²|S|;
  - E(g) ≤ E(S), since 0 ≤ g ≤ 1_S;
  - apply F1 and |S| ≍ k₀H′.
- (b) The physical ℓ¹ norm is ‖g‖₁ ≤ 2|S|, so ‖ĝ‖_∞ ≤ 2|S|.
  - Since g ≤ 0 off S: θ|S| ≤ ⟨g, 1_S⟩ = (1/q)Σĝ·conj(1̂_S) ≤ (1/q)‖ĝ‖_{4/3}‖1̂_S‖₄.
  - Here ‖ĝ‖_{4/3} ≤ ‖ĝ‖_∞^{1/4}‖ĝ‖₁^{3/4} and ‖1̂_S‖₄ = (qE(S))^{1/4}.
  - Rearranging gives (1/q)‖ĝ‖₁ ≥ cθ^{4/3}|S|/E(S)^{1/3}. The bound ≥ 1 is ‖ĝ‖₁/q ≥ sup|g|.
  - The ratio bound follows by maximising θ|S|/max(1, cθ^{4/3}H′^{2/3}/L) over θ; the maximum is at θ ≍ H′^{−1/2}L^{3/4}.
- The signed Hölder step was first sketched in GPT-6.1-sol's working log (codex/cov38b_1.log). It is not in its final
  output, which uses the different Lemma F2 route. We re-derived it independently. ∎

The same k₀√H′ ratio bound holds for indicators g = 1_A of an arbitrary A ⊆ S, by a direct route:
- E(1_A) ≤ |A|·max_{x∈S} r₃(x), where r₃(x) = #{(r₂, r₃, r₄) ∈ R³ : r₃ + r₄ − r₂ = x}. The Freiman lift needs x ∈ S.
- For fixed denominators the solutions in (k₂, k₃, k₄) lie in a coset of a sublattice. With c_i the largest of
  c₂, c₃, c₄ (c_i ≤ H′² ≤ 2k₀ + 1), the density of the coset is B_i/lcm(B₂, B₃, B₄).
- So r₃ ≤ 4(2k₀+1)²Σ_{B₂,B₃,B₄} B_min/lcm(B₂, B₃, B₄) ≪ k₀²H′(log H′)^8. The critic verified this bound numerically.

**Lemma F2 (the largest nonzero Fourier coefficient of a Farey-arc set; proved).** This is GPT-6.1-sol's blind-round
argument, re-derived and checked here. Let 2hK < q and S = S(h, K). Then

  max_{t≠0} |1̂_S(t)| ≤ (K + h log q)·q^{o(1)}.

*Proof.*
- **Möbius decomposition.** S is the disjoint union over d ≤ h of S*_d = {k/d : gcd(k, d) = 1, |k| ≤ K}. The
  reduced fractions are distinct mod q since 2hK < q. By Möbius inversion,
  1_{S*_d} = Σ_{e|d} μ(e)·1[(d/e)^{−1}·[−K/e, K/e]].
- **Each piece.** A dilated interval m̄·[−L, L] has |transform at t| ≤ min(2L+1, q/(2|v_m|)), where v_m := t m̄ mods q.
  For fixed m = d/e the pieces satisfy Σ_{e≤h/m} min(K/e, A) ≤ min((h/m)A, K(1 + log h)).
- **Counting small v_m (divisor bound).** m·v_m ≡ t (mod q) and |m v_m| ≤ UT for m ≤ U, |v_m| ≤ T. So m divides one of
  O(1 + UT/q) nonzero integers t + nq of size ≤ hq. Hence #{m ≤ U : |v_m| ≤ T} ≤ (1 + UT/q)q^{o(1)}.
- **Zero term.** For e > K the interval is {0}, with transform 1. These and the "+1"s add O(h log h) in total.
- **Assembly (corrected after both reviewers).** Take a cell with m ∈ (U/2, U] and |v_m| ∈ (T/2, T].
  - Each m in it contributes ≤ min(2hq/(UT), K(1 + log h)).
  - The cell holds ≤ (1 + UT/q)D_q such m, with D_q = max_{n≤q²}τ(n).
  - So the cell contributes ≤ D_q[K(1 + log h) + (UT/q)·2hq/(UT)] = D_q(K log h + 2h).
  - Summing over the O(log h·log q) cells gives (K + h)·D_q·polylog(q). ∎
- **Size of the loss.** D_q = exp(O(log q/log log q)) is not a polylog, and this route cannot do better.
  - For t with many divisors d ≤ h and tK ≪ q, the arcs with d | t add in phase, so |1̂_S(t)| ≳ K·#{d ≤ h : d | t}.
  - The proof critic measured |1̂_S(t₀)|/K = 26–47 at primorial-type t₀ (q = 10¹⁵ + 37, h = 64), against ≤ 0.01 at random t.
  - With Wigert's constant, the proved loss at q ≈ 2^512 is ≈ q^{0.12}, about the size of the 3/25 gap. So the threshold
    statements below are asymptotic (fixed c > 0, q → ∞).
- Only 2hK < q is needed (no Freiman condition, no K ≥ h).
- Unregistered numerical proof check (q ≈ 10⁶, 10⁷; 4 ≤ h ≤ 64): max_{t≠0}|1̂_S|/(K + h ln q) = 2.7–11.

**Theorem B (sharp leading-exponent threshold for phase-discarding certificates; proved).** Let c ∈ (0, 1/2] be fixed,
Y ≥ Y₀(c) (so that 2HK′ < q, i.e. Y > 4H², J^{2c} > 4), and κ = Jq^{−δ}.
- *Uniform κ is no restriction for Pos∀* (proof critic): a bound κ_h valid for all a must satisfy
  κ_h ≥ max_a|K(ah)| = max_{x≠0}|K(x)|.
- **(i) Sufficiency.** If max_{h≠0}|K(h)| < J(k₀+1)/q with k₀ = ⌊qH/(4Y)⌋, then Pos∀(c) holds. In exponents, this
  holds as soon as δ > δ*(c) := (1+2c)/(7−2c) + o(1).
  - *Proof.* Take G = Σ_{B≤H/4} Fej_B, with Fej_B(B̄k) = 1 − |k|/(k₀+1) the Fejér triangle on the arc B̄[−k₀, k₀].
  - Ĝ ≥ 0, so W(G) ≤ G(0) = H/4. Also ΣG = (H/4)(k₀+1), and supp G ⊆ S₀. Apply (★) and Prop. T(i).
  - Note that **H′ cancels**: the single arc B = 1 does exactly as well. By (ii), no phase-discarding certificate
    can use the union over dilates to lower the threshold by a power of q.
- **(ii) Necessity, every certificate supported (positively) in S_out.** Let g be any real function with g ≤ 0 off
  S_out = S(H, K′), K′ = 2qH/Y. Any scaling, sign or mass is allowed on S_out. Then Σg·(1 − |S_out|/q) ≤ W(g)·max_{t≠0}|1̂_{S_out}(t)|.
  - So (★) forces κ < J·Σg/(qW(g)) ≤ J·K′q^{o(1)}/q = J(2H/Y)q^{o(1)}, i.e. **δ ≥ δ*(c) − o(1)**.
  - *Proof.* g ≤ 0 off S_out gives Σ_x g·1_{S_out} ≥ Σg. Write the left side as
    (1/q)Σ_t ĝ(t)·conj(1̂_{S_out}(t)). The t = 0 term is Σg·|S_out|/q, and the rest is ≤ W(g)·max_{t≠0}|1̂_{S_out}|.
    Then apply Lemma F2 with h = H, K = K′ ≥ H log q.
  - No Freiman condition is needed, only 2HK′ < q. So (ii) holds for every fixed c ∈ (0, 1/2] and q ≥ q₀(c).
  - This gives sufficiency for δ > δ* + ε and necessity δ ≥ δ* − o(1). Equality and log-scale boundary behaviour are
    undecided. In particular c = c(q) → 0 (log-factor improvements over Thue) is not excluded: the Parseval floor
    exceeds the required κ only by Y^c(log Y)^{1/2−c}.
  - (This supersedes the draft's weaker "any certificate ≥ (1/2+3c)/(7−2c) = 13/50", which came from the energy
    route of Lemma F(b). GPT-6.1-sol's eq. (17)–(18) closes the window.)
- **(iii) Impossibility.** δ*(c) = (1+2c)/(7−2c) exceeds the Parseval cap 1/(7−2c) (Lemma P) **for every fixed c > 0**; they
  coincide at c = 0, i.e. the Thue level. At the target c = 3/8: required δ* = 7/25 = 0.28, cap 4/25 = 0.16,
  gap 2c/(7−2c) = 3/25. **No pointwise bound of any quality, known or conjectural, certifies any fixed c > 0
  through (★), whatever certificate is used.** ∎
- Novelty: the literature critic found no prior bound on Wiener (Fourier ℓ¹) norms of Farey-arc or major-arc sets;
  the nearest are Littlewood-type polylog bounds (Green–Konyagin 2009; Sanders). "No prior formulation found by us."
- Lemma F (energy route) remains as an independent check. It shows the Fourier ℓ¹ norm of full-mass nonnegative
  minorants is ≥ H/polylog, with an explicit polylog. Lemma F2 gives the same conclusion up to q^{o(1)} for all
  certificates. GPT-6.1-sol (its eq. (15)) gives the matching upper bound A(1_S) ≪ h log h log q.

**Remark B2 (joint phase discarding over (B, ℓ); proved; critic's suggestion).** The JH′ points Bℓ̄ (B ≤ H′, ℓ ∈ 𝒫)
are distinct mod q: Bℓ′ ≡ B′ℓ with both sides < q forces Bℓ′ = B′ℓ, and then ℓ = ℓ′ since ℓ ∤ B.
- The hit count is Σ_{B,ℓ} 1_I(aBℓ̄), which is a single-arc certificate over JH′ points.
- Parseval forces max_m|Σ_{B,ℓ} e(maBℓ̄/q)| ≥ √(JH′/2). With a Fejér triangle on the single arc, W ≍ 1 (a sharp arc
  has W ≍ log).
- So phase discarding needs λ ≳ √(JH′), i.e. J^{2−4c} ≳ J^{3/2−c}, i.e. c ≤ 1/6. It fails at c > 3/8.

*Remarks.*
- coverage38 Prop. L(b) obtained "c ≤ 1/4" from W ≥ 1/2 − μ. Theorem B identifies the true price: the mass-to-norm
  ratio of any certificate is ≤ K′q^{o(1)}, so the required pointwise saving is exactly δ*(c).
- E2 measures W(g_tri)/H′ = 0.41–0.65 and slope 0.80–0.85 in H′ (P2b). Several Wiener-scan points (e.g. q = 10⁷,
  H′ = 32, k₀ = 38) violate k₀ ≥ H′³. Their Hölder statistics are still valid, but they do not test Lemma F's regime.

**Proposition B′ (adaptive phase discarding with exact magnitudes; proved).** Suppose 16Y³ < q. Then
1/ℓ₁ + 1/ℓ₂ = 1/ℓ₃ + 1/ℓ₄ (mod q) forces {ℓ₁, ℓ₂} = {ℓ₃, ℓ₄}: the numerators are < q, and for primes the reduced
fractions determine the pair. So the inverses of band primes have additive energy exactly 2J² − J.
- This gives Σ_x|K(x)|⁴ = q(2J² − J) and Σ_x|K(x)|² = qJ.
- Hölder then gives mean_{x≠0}|K(x)| ≥ (1 − o(1))√(J/2).
- For any fixed g, averaging over a gives E_a Σ_{h≠0}|ĝ(h)||K(ah)| = q·W(g)·mean|K|. Paley–Zygmund applied to
  Z(a) = Σ|ĝ(h)||K(ah)|, with E[Z²] ≤ (qW)²·mean|K|², gives: **for at least (1 − o(1))/8 of all a ∈ F_q^***, the
  a-specific phase-discarded error is ≥ (1/2)(1 − o(1))√(J/2)·W(g).
- So for a fixed weight g (independent of a), the worst-case phase-discarded error with exact magnitudes |K(ah)|
  cannot fall below a constant multiple of √J·W(g). By Lemma F2 that is already too large for every fixed c > 0
  (Theorem B(iii)). Weights g_a that depend on a are not covered.
- For the prime-only Pos, the averaging covers a positive proportion of residues but cannot by itself say where the
  primes p ∈ (q/2, q) fall.
- This is where Bourgain–Garaev-type energy information enters, and it is already optimal: the energy is exactly
  minimal. Precedent: Bourgain–Garaev 2014 Thm 6 bounds the 2k-fold energy of inverses of primes ≤ N by
  (2k)^k(N^{2k−1}/p + 1)N^k; for k = 2 this is J₄ < 16(N³/p + 1)N². Our exact identity is the elementary sharp form.

**Lemma M (multiplicative characters, GRH included; proved).** Expand the hit count with characters mod q:
N(a) = (1/(q−1))Σ_χ χ̄(a)S_𝒫(χ)F̂(χ), where F̂(χ) = Σ_{|k|≤k₀}Σ_{B≤H′}χ(k)χ̄(B).
- This concerns a common uniform bound over χ ≠ χ₀; character-dependent exact magnitudes are a different question.
- Parseval over characters gives max_{χ≠χ₀}|S_𝒫(χ)| ≥ √(J/2).
- (1/(q−1))Σ_{χ≠χ₀}|F̂(χ)| ≥ #{(k, B) : k = B} − O(λ/J) ≥ H′ − 1.
- Phase discarding therefore needs λ ≳ √J·H′, i.e. J^{1/2−c} ≳ Y^{1/2}. With J ≍ Y/log Y, this fails for every
  c ≥ 0.
- With arbitrary nonnegative weights F′ supported on the same product set in place of F, it still needs λ ≳ √J,
  i.e. c ≤ 1/4. GPT-6.1-sol's Cauchy–Schwarz version gives the threshold δ > 1/2 + αc (31/50 at c = 3/8).
- GRH gives |S_𝒫(χ)| ≪ √Y log²(qY), at the Parseval level, so GRH used this way gives c ≤ 1/4 at best.

**Lemma S (scale; proved in the model).** By Brun–Titchmarsh (coverage38 §1), a group with modulus M ≤ P^{1−ε} has
π_g ≤ C(log P)/M. For M near P, only π_g ≲ log P/M is available.
- A design of cost X has ≤ X groups.
- So the groups with M > Z carry λ_{>Z} ≤ C(log P)X/Z ≤ C/log P for Z = X(log P)².
- Hence all but a fraction C/log P of the pairs must be covered by moduli ≤ X·polylog ≤ q^{2/5+o(1)}. There every
  known prime-variable Kloosterman bound mod a prime is trivial.
- The leftover pairs could be exactly the hard ones, covered by large moduli. λ_{>Z} < 1 (a population average)
  rules out covering the whole population with the large-modulus prong alone. It does not, by itself, rule out large
  moduli covering a residual exceptional set.

### 2.2 The table (Task 2's question: "the exact exponent each tool gives")

δ is the pointwise saving in max_{h≠0}|K(h)| ≤ J q^{−δ} at the length the design needs (q^{0.32} at c = 3/8). "Certified
c" is the best coverage exponent the tool yields for Pos(c). Thue/Minkowski (c = 0, exponent 1/5) is always available
and is the baseline.

| tool (references in §6) | applies at length q^{0.32}? | δ it gives | certified c | reaches c > 3/8? |
|---|---|---|---|---|
| Thue / Dirichlet / Minkowski, one modulus (geometry, no analysis) | yes | — | **0** (exponent 1/5; Harvey's) | no |
| Prime-variable Kloosterman sums: Bourgain 2005 (prime q, q^{1/2+ε} ≤ N ≤ q, saving ≍ ε⁴); Garaev 2010 (prime q, nontrivial from q^{3/4+ε}); Fouvry–Shparlinski 2011 (any q, q^{3/4+ε} ≤ N ≤ q^{4/3−ε}); Baker 2012 (squarefull part ≤ q^{1/4}, q^{1/2+ε} ≤ N ≤ q); Korolev 2018 (elementary, q^{1/2+ε}) and 2019; Korolev's composite-modulus results (from q^{5/8+ε}) | **no**: nontrivial only from q^{1/2+ε}, and Lemma S puts the moduli below q^{2/5} | 0 | none (nothing beyond Thue) | no |
| Irving (averaged over moduli q) | no: one fixed q here | — | none | no |
| Integer-variable incomplete Kloosterman: Karatsuba; Bourgain–Garaev 2014 Thm 16 (prime q: N(log log q)³log q/(log N)^{3/2}, i.e. (log q)^{1/2−o(1)} at N = q^{0.32}); Bourgain–Garaev, Acta Arith. 164 (2014) Thm 5 (any m, saving (log m)^{1/2}/(log log m)^{O(1)}); Korolev 2016 (prime q, integer variable, short ranges) | yes | log saving | none (δ = o(1) < δ*) | no |
| Bourgain's multilinear Kloosterman sums; Bourgain–Garaev 2014 Thms 8, 11–13 and Cor. 2; prime-variable version Bourgain–Garaev, Acta Arith. 164 (2014) Cor. 3 (moduli = products of primes: two primes of size q^{4/25} give δ ≈ 1/800 at total length q^{8/25}; GPT-6.1-sol and the literature critic agree) | yes | small power, e.g. 1/800 | none (δ ≪ δ* = 0.28) | no |
| Heath-Brown identity + Weil/Deligne completion | Type I sums of length < q^{1/2}: completion costs √q·log q > trivial | 0 | none | no |
| Bilinear sums with phases e_q(am n̄): Kerr–Shparlinski–Wu–Xi, JLMS 108 (2023), Thm 2.4 (r = q case), nontrivial for all M, D ≥ q^ε (cited by GPT-6.1-sol; the literature critic reads that paper as mainly about complete sums, so this entry is UNVERIFIED by us); Bourgain–Garaev 2014 Thms 7–8 | yes (phases kept in the (m, n) form) | power saving | none. Applied to the hit count with m = tB (GPT-6.1-sol), the error is ≪ HY^{3/4}q^{o(1)} = q^{7/25} against the main term λ = q^{2/25}; it needs H > Y^{3/4} | no |
| Hooley R* / square-root cancellation (conjectural) | — | 4/25 (the Parseval cap) | **none** (Theorem B: needs 7/25) | no |
| GRH, multiplicative characters | yes | Parseval level | ≤ 1/4 with arbitrary weights, ≤ 0 with the natural one (Lemma M) | no |
| Bourgain–Garaev additive energy / sum-product for {ℓ̄} | yes | energy exactly minimal (Prop. B′) | none (B′) | no |
| Large sieve / k-th moments (D′, sieve2d); centred moments over a | yes | almost-all only (pair-population moments of order ≤ 4 proved; no Poisson-quality per-q moment over a known) | worst case needs order ≥ 26 over a (per q) or > 50 over pairs; nothing known beyond order 4 | no |

**The exact gap (Task 3, "isolate").**
- *Phase-discarding route.* Positivity needs δ > (1+2c)/(7−2c) = 7/25. The most any bound can give is 1/(7−2c) = 4/25.
  The deficit is 2c/(7−2c) = 3/25 in the exponent of q. Equivalently, the hit count needs an error o(λ) = o(J^{1/4})
  on a sum of J terms: a saving of J^{3/4} over trivial, J^{1/4} beyond square root.
- *Known.* The saving is 0 for a single band of primes at q^{0.32}. For integers it is (log q)^{1/2−o(1)}. For products of two
  band-scale primes (composite moduli) it is δ ≈ 10⁻³ (Bourgain–Garaev bilinear bounds). All of these are far below
  4/25.
- *What would suffice instead.* A positivity mechanism that does not discard phases term by term. Either a
  phase-retaining multilinear estimate (none is known at q^{0.32}) or a non-analytic one. The only one we know is
  single-modulus geometry (Thue/Minkowski), and B2 shows it stops at c = 0. No joint-geometry statement for unions
  over ℓ is known, even with a log-factor gain over Thue.

### 2.3 Bourgain–Garaev in detail (the shortlist's first question)

conn_alg asks: "Do Bourgain–Garaev-type multiplicative-energy bounds for {ℓ̄ : ℓ prime in a band} extend from length
q^{1/2+ε} to q^{0.32} for positivity in a union of arcs?" Three separate inputs answer it.
- **Exponential-sum bounds** (Bourgain 2005; Bourgain–Garaev 2014, multilinear Kloosterman with x_i ≤ q^ε). Every such
  bound is pointwise in the frequency, so used through a phase-discarding certificate Theorem B applies whatever its
  strength. Even square-root cancellation
  (δ = 4/25) falls short of δ* = 7/25.
- **Energy / sum-product bounds** for reciprocals. At our scale the energy of {ℓ̄} is already exactly minimal
  (Prop. B′), so better energy bounds cannot help. Prop. B′ also shows that energy information, used through
  averaging over a, controls only a positive proportion of dilates. The worst case asks for all of them.
- **Covering theorems** ("every residue is a sum of k reciprocals of integers ≤ q^ε": Shparlinski, Arch. Math. 78
  (2002), k ≪ ε^{−3}; Glibichuk, Math. Notes 79 (2006), k ≪ ε^{−2}). Bourgain–Garaev 2014 gives sumset *sizes*
  (Cor. 1) and counting lemmas for products (Lemmas 10 and 12), not covering. These cover F_q by additive k-fold sumsets of {x̄}. The target asks to cover F_q^* by the
  multiplicative dilates ℓ·S₀, i.e. by the product set 𝒫·[−k₀, k₀]·[1, H′]^{−1}.
  - For product sets the known all-classes ranges are the modular-hyperbola ones: xy ≡ a is known to be solvable for
    all a once XY ≫ q^{3/2+o(1)} (Weil). This is a sufficient range, not a barrier. The trivial character bound needs |A||B| > q².
  - The target's product has size λq = q^{1+2/25}, far below both.
  - Thue reaches product size q only because one factor, [1, H′]^{−1}·[−k₀, k₀] with H′k₀ ≥ q, is a full Dirichlet box.
- **Answer: no.** The obstruction is not a lack of distribution information. It is that positivity at
  λ < J^{1/4} is a statement below the fluctuation scale of every phase-discarding certificate (Theorem B), uniformly in
  the dilate. Phase-retaining multilinear estimates are not excluded, but none is known at this length.

### 2.4 Direct answer to the triage claim: does any tool give *positivity* (not an average) at length q^{0.32}?

GPT-6.1-sol's triage (codex/r29_triage.out, lead D, ≈ 1%) says that positivity in the tiny Farey-arc unions needs
arithmetic information beyond generic cancellation in sums of inverses. We tested this directly. **None of the
estimates checked here proves positivity at length q^{0.32}.** They fall into four kinds:

| kind of statement | tools | what it gives at q^{0.32} | positivity for every a (every pair)? |
|---|---|---|---|
| geometry, one modulus | Thue / Dirichlet / Minkowski | positivity for **every** a, but only at λ ≍ J (c = 0) | yes, at c = 0 only (exponent 1/5) |
| average over a (second or higher moments) | large sieve, Montgomery–Vaughan; sieve2d fourth moment; energy of {ℓ̄} (Bourgain–Garaev type, here exactly minimal) | uncovered fraction ≤ C/λ, C/λ² | **no**: average-type. A worst case needs a centred moment of order ≥ 26 over a, which is unknown (Prop. 0(c)) |
| pointwise exponential and character sums, used by phase discarding | prime/integer/multilinear Kloosterman, Heath-Brown + Weil, R*, GRH | positivity iff δ > δ*(c) = 7/25; any bound has δ ≤ 4/25 (Parseval); known δ = 0 for primes | **no** (Theorem B): the only tools that are pointwise in a are capped below the threshold |
| covering theorems for sets in F_q | sums of k reciprocals (Shparlinski 2002, Glibichuk 2006); xy ≡ a in boxes (Garaev 2006: |x|, |y| ≤ Cq^{3/4} for all a); Erdős–Odlyzko–Sárközy | positivity, but for additive k-fold sumsets, or for product sets of size ≥ q^{3/2} | **no**: the target is a multiplicative covering of size q^{1+2/25} |

- The triage claim is confirmed and made quantitative.
- Among the results checked, positivity for every residue appears in only two regimes: single-modulus geometry
  (c = 0), and harmonic analysis at product size ≥ q^{3/2} (or ≥ q² with characters).
- Phase-retaining bilinear bounds uniform in a do exist at this length (Bourgain–Garaev Thms 7–8; Kerr–Shparlinski–Wu–Xi).
  Applied to the hit count they leave an error ≈ q^{7/25} against λ = q^{2/25} (§2.2).
- The pointwise exponential- or character-sum bounds, used by phase discarding (Theorem B), need
  a saving δ > (1+2c)/(7−2c) = 7/25 in max_h|Σ_ℓ e(hℓ̄/q)|.
- No bound can exceed δ = 1/(7−2c) = 4/25, and the known value at this length is δ = 0.
- **The exact gap is 2c/(7−2c) = 3/25 in the exponent of q, beyond square-root cancellation.** In counting terms: an
  error o(J^{1/4}) is needed on a sum of J terms, against the J^{1/2} floor of every phase-discarding argument.
- Average-type tools fail for a different reason. They bound the number of uncovered a's by about qλ^{−r}, and this
  drops below 1 only at moment order 2r ≥ 26 (per q). Nothing beyond order 4 is known, and the 4th moment is already
  a Kloosterman-type problem at this length. A Poisson-quality 26th-moment estimate is the one averaged route not
  excluded (GPT-6.1-sol).

## 3. Prove or isolate (Task 3)

- **Nothing is proved toward c > 3/8**, nor toward any c > 0 for all pairs, nor toward a log-factor union gain over
  Thue.
- What is proved:
  - Prop. T, the precise target and its equivalent forms;
  - Prop. 0: no intermediate regime; averaged forms help only through moment bounds strong enough to empty the
    exceptional set (order ≥ 26 per q), none of which is known;
  - Theorem B, Remark B2, Lemmas F0, F1, F, P, M and S, and Prop. B′: the exact exponents and gap of all
    phase-discarding uses of the known tools.
- **Dependency list for Theorem B.**
  - Lemma F2: Möbius decomposition, Dirichlet-kernel bound, and the divisor bound τ(n) = n^{o(1)}.
  - Lemma P (Parseval).
  - Prop. T(i)–(ii) (coverage38 §3.1 with the corrected centre).
  - Fejér's kernel (Ĝ ≥ 0).
  - The size relations of Theorem W (Y = q^{2/(7−2c)+o(1)}, J ≍ Y/log Y by PNT in progressions mod 4).
  - Nothing else; in particular no unproved hypothesis.
  - Lemmas F0, F1 and Hölder support only Lemma F, the independent energy-route check.
- **Isolated gap.** See §2.2: 7/25 needed against at most 4/25 possible for phase-discarding arguments, and 0 known
  at q^{0.32}. A proof must retain phases (a multilinear estimate with the Farey coefficients), or be a non-analytic
  positivity argument for unions over moduli.

## 4. Two-pronged algorithms (Task 4)

**Theorem TP (proved in the coverage38 model; Lehman prong per Prop. D(a)).** Let 𝒜 run DCL-wc with a band design 𝒟
and a Lehman/Harvey prong ℒ (stationary or full windows at its own moduli). Each prong is N-uniform, both are charged
per group, and WIN holds.
- (i) **Success.** 𝒜 factors every pair covered by 𝒟 (Theorem R) or by ℒ (Prop. D(a), or Theorem R at ℒ's moduli).
- (ii) **Same threshold.** 𝒟 ∪ ℒ is one design of the family. Theorem U gives X ≥ cN^{1/5}(λ⁴/J)^{1/5} with
  λ = λ_𝒟 + λ_ℒ and J = J_𝒟 + J_ℒ, when the modulus sets are disjoint. Otherwise J counts distinct moduli, and
  merging duplicate (v, M) groups changes λ by at most the overlap. The estimates survive up to constants, with WIN
  checked for the merged windows. No two-pronged combination of prongs from the family has a weaker coverage
  threshold than λ < J^{1/4}.
- (iii) **The Lehman prong covers few pairs.** Theorem U applied to ℒ alone gives λ_ℒ ≤ C(X/P^{2/5})^{5/4}J_ℒ^{1/4}.
  At X = N^{1/5−η} with J_ℒ = 1 (ordinary Lehman), the covered fraction is ≤ CN^{−5η/4}. A stationary-window prong
  with J_ℒ moduli covers up to CN^{−5η/4}J_ℒ^{1/4}.
  - The direct count agrees: Pr[L* ≤ X] ≈ Σ_{h(v)≤√X} X√(P/h)/(hP) ≍ X^{5/4}/√P.
  - E3 measures slopes −1.31, −1.28, −1.27 against −5/4.
- (iv) **Residual statement.** 𝒜 handles every pair if **Res(η)** holds: every pair with L*(p, q) > N^{1/5−η} is
  covered by 𝒟. Here λ_𝒟 ≤ CJ_𝒟^{1/4}N^{−5η/4}, by Theorem U applied to 𝒟 itself. This is sufficient. It is
  necessary only for proofs through the coverage guarantee, since accidental collisions are not excluded.
- (v) **B3 is escaped.**
  - For p = ax + b, q = cx + d with fixed coefficients, δ_{(a,c)/g} is constant, so L* = O(1). Linear families never
    reach Res(η).
  - So Res(η) does not imply the interval statement U(x, H) via (x, x+1).
- (vi) **A survivor remains for coprime-pair proofs.** Let α ∈ (1, 2) be badly approximable with |α − B/A| ≥ c_α/A².
  - Then |Bx − A⌊αx⌋| ≥ c_αx/h − h, and minimising over h gives L*(x, ⌊αx⌋) ≥ c′_α P^{2/5}.
  - So these pairs are Lehman-hard at every budget N^{1/5−η}, and any proof of Res(η) valid for all coprime pairs
    implies band-prime coverage of the sparse sets {Bx − A⌊αx⌋ : |A|, B ≤ H}. These are ≍ H² integers per x, not an
    interval.
  - Prime-pair-only proofs are not obstructed this way.
- (vii) **Other second prongs.**
  - p − 1 (Ward–Maestro Prop. 35.3: a positive relative proportion of pairs in (x^{1−ε}, x], not S_n, and not a
    worst-case result), Hittmeir/Harvey, and class-group methods: none has a
    proved easy set containing the uncovered pairs.
  - Uncoveredness is a condition on (pq̄ mod ℓ)_ℓ that carries no local structure (inversesieve Prop. D, Theorem L).
    p − 1 smoothness and orders are unrelated data (coverage38 §4(f)).
  - We found no provable disjointness; that none exists is UNVERIFIED.

*Proof of (iii) and (vi).*
- (iii) Theorem U: X ≥ cP^{2/5}λ^{4/5}J^{−1/5}, solved for λ. The covered fraction is ≤ λ by the union bound.
- (vi) |Bx − A⌊αx⌋| = |x(B − Aα) + A{αx}| ≥ x·c_α/A − A.
  - With h = max(|A|, B) ≍ A (B/A ≈ α), L* ≥ min_h max(h², (c_αP/h − h)√(h/P)).
  - The two terms balance at h ≍ (c_α√P)^{2/5}, which gives L* ≍ c_α^{4/5}P^{2/5}. ∎

**Answer to Task 4.** A two-pronged worst-case algorithm handles every pair iff (for coverage-based proofs) Res(η)
holds. Res(η) has the same λ < J^{1/4} threshold, applies to all but an N^{−5η/4} fraction of pairs, and for coprime
pairs retains B3's sparse survivor. Our assessment, which Theorem U does not rank: the dichotomy changes the
quantifier, not the difficulty.

## 5. Experiments (pre-registered; `notes/cov38b_prereg.md`, addendum A)

All 10 named frozen predictions pass (P1a–d, P2a–d, P3a–b).
- The analysis script checks grid completeness. An exact E2 Wiener-grid assertion was added after the run, at
  GPT-6.1-sol's request; it is output-neutral, and the verdicts are identical. GPT-6.1-sol also recomputed the full
  summary from the 260 E2 records.

### 5.1 E1: worst-case coverage height at the c = 3/8 scaling

Exhaustive over every prime p ∈ (q/2, q). Band: primes ≡ 3 (mod 4) in (Y, 2Y], Y = ⌊q^{8/25}⌋. Signed weights plus (1,1).

| Q0 | q | Y | J | n_p | Thue ⌊√ℓ_min⌋ | h* | model [H₀.₀₅, H₉₅] | c_emp |
|---|---|---|---|---|---|---|---|---|
| 10⁵ | 100003 / 130003 / 170003 | 39 / 43 / 47 | 5 / 6 / 5 | 4.5k–7.2k | 6 / 6 / 7 | 6 / 6 / 7 | [6,6] / [6,6] / [7,7] | 0.02 / 0.05 / −0.01 |
| 10⁶ | 1000003 / 1300021 / 1700021 | 83 / 90 / 98 | 7 / 9 / 10 | 37k–61k | 10 | 9 / 9 / 9 | [9,9] / [8,9] / [8,9] | 0.01 / 0.02 / 0.04 |
| 10⁷ | 10000019 / 13000027 / 17000023 | 173 / 189 / 205 | 14 / 16 / 16 | 316k–520k | 13 / 13 / 14 | 11 / 11 / 12 | [11,11] / [11,11] / [11,12] | 0.07 / 0.08 / 0.06 |
| 10⁸ | 100000007 / 130000001 / 170000009 | 363 / 394 / 430 | 29 / 31 / 35 | 2.8M–4.6M | 19 / 20 / 20 | 13 / 13 / 13 | [13,14] ×3 | 0.11 / 0.12 / 0.13 |
| 10⁹ | 1000000007 / 1300000003 / 1700000009 | 758 / 825 / 898 | 53 / 61 / 63 | 24.5M–40.6M | 28 / 28 / 30 | 17 / 16 / 17 | [15,17] / [15,17] / [16,17] | 0.12 / 0.14 / 0.14 |

- **P1a** (h* ≤ Thue): pass.
- **P1b**: pass, 15/15 in band, 0 on the high side.
- **P1c** (h* strictly below Thue for Q0 ≥ 10⁷): pass.
- **P1d**: pass. For the 9 all-residue runs (every a ∈ [1, q−1]), h*_all − h* ∈ {0, 1}. Primality of a is invisible
  at the worst-case level.
- **Reading (experimentally supported).** At the exact Theorem-W scaling the worst case is the independent-model
  extreme. It already beats Thue by a factor 1.6–1.8 at q ≈ 10⁹. The finite-size exponent rises slowly (log
  corrections); the model's limit is c → 1/2.
- **Correction (GPT-6.1-sol adversarial round).** The draft said the data showed the target true "with room to
  spare". That is false.
  - At the prescribed height H = ⌊√Y·J^{−3/8}⌋ (4, 5, 6 for q ≈ 10⁷, 10⁸, 10⁹), the histograms contain 81 869,
    299 426 and 2 584 494 uncovered prime pairs. Every run has c_emp < 3/8.
  - The data fit the independent model and beat Thue. They say nothing about the asymptotic conjecture c > 3/8, which
    the model predicts only in the limit (c_emp → 1/2 very slowly), and nothing about provability.
  - The experimental Y = ⌊q^{8/25}⌋ depends on q. The algorithm needs an n-only Y(n), which changes constants only.

### 5.2 E2: Wiener norms and Erdős–Turán certificates

| q | J | κ/√(J ln q) | W(g_tri)/H′ at H′ = 2, 4, 8, 16, 32 | slope (H′ ≥ 4) | LB/H′ at H′ = 32 | H_ET farey / arc | ⌊√Y⌋ | h* (E1) |
|---|---|---|---|---|---|---|---|---|
| 100003 | 5 | 0.66 | 0.65, 0.56, 0.49, —, — | 0.81 | — | 52 / >76 | 6 | 6 |
| 1000003 | 7 | 0.71 | 0.65, 0.56, 0.51, 0.42, — | 0.80 | — | 96 / >164 | 9 | 9 |
| 10000019 | 14 | 0.85 | 0.65, 0.56, 0.51, 0.45, 0.41 | 0.85 | 0.137 | 176 / 332 | 13 | 11 |

- **P2a, P2b, P2c, P2d: all pass.**
- W(1_S) (sharp) is larger: 5–22, with log-k₀ tails.
- The Hölder lower bound LB/H′ decays from 0.42 to 0.14 over H′ = 2 … 32, consistent with a polylog loss in
  Lemma F. P2c's floor of 0.1 would eventually fail; this is disclosed.
- **Reading.** Certificates exist only at H ≈ Y (≈ 9–14 × √Y), far above the true worst case h* ≤ √Y. At these sizes
  κ ≈ 0.7–0.85·√(J ln q) is close to J itself, so the certificate needs almost all of F_q. This is the toy-size face
  of Theorem B.

### 5.3 E3: the Lehman prong

| n | frac(L* ≤ P^{2/5}/G), G = 2, 4, 8, 16, 32, 64 | slope (G = 8…64) | max L*/P^{2/5} |
|---|---|---|---|
| 48 | 0.725, 0.318, 0.127, 0.0484, 0.0191, 0.0084 | −1.31 | 0.90 |
| 64 | 0.741, 0.333, 0.139, 0.0579, 0.0234, 0.0098 | −1.28 | 0.95 |
| 80 | 0.752, 0.344, 0.145, 0.0600, 0.0246, 0.0103 | −1.27 | 0.96 |

- **P3a, P3b: pass.** The Lehman prong at budget P^{2/5}/G covers a fraction ≍ G^{−5/4}, as Theorem U predicts with
  J = 1.
- Asymptotically, at G = N^η, it covers N^{−5η/4} (Theorem TP(iii)).
- Dirichlet's bound max L* ≤ 2P^{2/5} holds (0.90–0.96).

## 6. Literature used (to be checked by the literature critic; §8)

- Bourgain, "More on the sum-product phenomenon in prime fields and its applications", IJNT 1 (2005): Kloosterman sums
  over primes, nontrivial for q^{1/2+ε}.
- Bourgain–Garaev, "Sumsets of reciprocals in prime fields and multilinear Kloosterman sums", Izv. Math. 78 (2014).
- Garaev, "Estimation of Kloosterman sums with primes and its application", Math. Notes 88 (2010).
- Fouvry–Shparlinski, "On a ternary quadratic form over primes", Acta Arith. 150 (2011): q^{3/4+ε} ≤ N ≤ q^{4/3}.
- Baker, "Kloosterman sums with prime variable", Acta Arith. 156 (2012).
- Korolev, Math. Notes 100 (2016) (integer variable, prime modulus); Math. Notes 103 (2018) (elementary, primes,
  q^{1/2+ε}); Math. Notes 106 (2019) (primes, q^{1/2+ε}).
- Bourgain–Garaev, "Kloosterman sums in residue rings", Acta Arith. 164 (2014) 43–64 (Thm 5; Cor. 3 for primes).
- Shparlinski, Arch. Math. 78 (2002); Glibichuk, Math. Notes 79 (2006): sums of reciprocals cover F_p.
- Garaev, Acta Arith. 124 (2006): xy ≡ a for all a with |x|, |y| ≤ Cp^{3/4}. Shparlinski's survey "Modular
  hyperbolas" (Jpn. J. Math. 2012) calls going below p^{3/2} out of reach.
- Irving, "Average bounds for Kloosterman sums over primes", Funct. Approx. Comment. Math. 51 (2014) 221–235,
  arXiv 1301.6372.
- Kerr–Shparlinski–Wu–Xi, J. London Math. Soc. 108 (2023) 578–621, Thm 2.4 (as cited by GPT-6.1-sol; UNVERIFIED by us).
- Matomäki–Teräväinen (ternary prime products covering residue classes): checked by GPT-6.1-sol; not applicable.
- Karatsuba's method. (Removed after the literature critic: Shkredov 2021 concerns complete Kloosterman sums, and
  Kerr–Shkredov–Shparlinski–Zaharescu concerns modular square roots. Neither is a Kloosterman-fraction result.)
- Heath-Brown's identity (Canad. J. Math. 34, 1982); Weil/Deligne.
- Modular hyperbolas: Shparlinski's survey; XY ≫ q^{3/2+o(1)} for all classes.
- Prior project notes: coverage38 (Theorem U, Prop. L, §3), dchworst (Theorems R, W; B1–B4), inversesieve, sieve2d,
  avgtime.

## 7. Candidate card

| field | content |
|---|---|
| mechanism | Pos(c): positivity of inverses of band primes in unions of Farey arcs, via Bourgain–Garaev / Kloosterman / character tools |
| payoff if true | deterministic worst case N^{(1−c)/(7/2−c)} < N^{1/5} for c > 3/8 |
| status | **refuted for every phase-discarding (term-by-term in h) certificate supported in the outer Farey set, at every fixed c > 0, with the sharp leading threshold δ*(c) = (1+2c)/(7−2c) against the Parseval cap 1/(7−2c) (Theorem B, proved)**; no intermediate exponent (Prop. 0); averaged forms need an unknown ≥ 26th moment; two-pronged designs face the same threshold (Theorem TP) |
| experimentally supported | worst-case heights follow the independent-model extreme and beat Thue, but at toy sizes they lie far above the prescribed c = 3/8 height (E1); Wiener norm ≍ H′ (E2); Lehman prong ≍ G^{−5/4} (E3) |
| sharpest obstruction | a positivity statement at λ = J^{1/4} is below the fluctuation scale √J of every phase-discarding certificate; for every phase-discarding certificate the union over dilates is invisible (H′ cancels in Theorem B(i)) |
| what would revive it | (a) a phase-retaining multilinear estimate with the prescribed Farey coefficients, error o(J^{1/4}); (b) a Poisson-quality centred moment over a of order ≥ 26; (c) a joint geometry-of-numbers positivity for unions of lattices Λ_ℓ = Z(p, q) + ℓZ². As a first test of any of these: a log-factor union gain over Thue |

## 8. Reviews

**Opus proof critic (read-only; scripts in the session scratchpad): ACCEPT WITH FIXES.**
- *Verified:*
  - Prop. T(i)–(ii), exhaustively on 6096 reciprocal hits and 79 497 coverage events;
  - Lemmas F0, F1 (index counting, divisor sums, exact energies for H′ ≤ 6), F(a), F(b) (with random signed tests);
  - Theorem B(i) (Ĝ ≥ 0 numerically);
  - all exponent arithmetic (7/25, 13/50, 4/25, 3/25), the Freiman condition 8H⁴ < Y for c > 1/4;
  - Prop. B′ (exact energy, Hölder, Paley–Zygmund), Lemma M;
  - Theorem TP(iii), (v), (vi) (L*/P^{2/5} ∈ [0.64, 0.77] for α = golden ratio up to x = 10²⁴); Prop. 0(b).
- *Fixes applied:*
  - **(major)** "every linear tool fails" rescoped to *phase-discarding* certificates, with phase-retaining estimates
    and ℓ- or a-dependent weights stated as not excluded;
  - the moment count corrected;
  - Theorem B(ii) normalisation and k₀′;
  - the "H′ cancels" scope (later made unconditional by Lemma F2);
  - B(iv) range;
  - the F1 sharpness argument (disjoint A*_B);
  - the indicator remark's proof;
  - Prop. T(iii) hypotheses;
  - Prop. B′ wording ((1 − o(1))/8; fixed g only; 16Y³ < q);
  - Lemma S (Brun–Titchmarsh source, C/log P, the leftover pairs);
  - the validity definition of S;
  - ⊆ instead of =;
  - the modular-hyperbola range stated as sufficient;
  - "consistently" for Theorem W;
  - assessments labelled as such.
- It suggested Remark B2, which was added.

**GPT-6.1-sol blind round (`codex/cov38b_1.out`).** It reached the same verdict independently.
- *Integrated, with attribution:*
  - δ* = (1+2c)/(7−2c) and the 7/25 vs 4/25 comparison, found independently;
  - the signed-minorant argument of Lemma F(b);
  - **Lemma F2** (max_{t≠0}|1̂_S| ≪ Kq^{o(1)}), which upgrades Theorem B(ii) to every certificate and every c > 0;
  - the upper bound A(1_S) ≪ h log h log q;
  - the correction that a sufficiently strong averaged (26th centred moment) estimate *would* give the worst case;
  - the Bourgain–Garaev Thm 8 computation (δ ≈ 1/800 for products of two primes);
  - the bilinear (m = tB) computation giving error HY^{3/4};
  - the multiplicative-character threshold δ > 1/2 + αc;
  - the J_L^{1/4} factor for multi-modulus Lehman prongs;
  - the table of failed non-linear attempts (sum-product growth, reciprocal sums with the height loss
    B ≪ Y^{1−1/r}, joint geometry, pigeonhole on shortest vectors, sieve over p, polynomial method).
- *Its corrections to our framing:*
  - Theorem B's obstruction concerns c > 0, not the recovery of Thue at c = 0 (S = F_q once hK ≥ q);
  - R* is an integer-variable hypothesis.

**Opus literature critic (read-only): no conclusion changes.**
- *Fixes applied:*
  - Shkredov 2021 (complete sums) and Kerr–Shkredov–Shparlinski–Zaharescu (modular roots) removed as Kloosterman-
    fraction sources;
  - the covering theorems for sums of reciprocals credited to Shparlinski 2002 and Glibichuk 2006, not to
    Bourgain–Garaev;
  - the integer saving corrected to (log q)^{1/2−o(1)};
  - Baker's range corrected (any modulus with squarefull part ≤ q^{1/4}, up to N ≤ q);
  - Korolev 2016/2018/2019 attributions;
  - Bourgain–Garaev Acta Arith. 164 (Thm 5, Cor. 3) added;
  - Bourgain–Garaev Izv. Thm 6 cited as the precedent for Prop. B′;
  - Garaev 2006 (xy ≡ a with |x|, |y| ≤ Cp^{3/4}) added;
  - the Ward–Maestro caveat added.
- *Not found by it:*
  - any all-residue positivity result for prime × small × inverse-small at total size q^{1+ε};
  - any all-a modular-hyperbola result below p^{3/2}, even with x prime;
  - any prior bound on Wiener norms of Farey-arc sets.
- *Unresolved:* the Kerr–Shparlinski–Wu–Xi Thm 2.4 reading (marked UNVERIFIED; it does not affect any conclusion).

**Opus proof critic, round 2 (Lemma F2 and the rewritten Theorem B): ACCEPT WITH FIXES.**
- *Verified:* the Möbius decomposition (exact formula against the FFT), the per-m bound, the divisor count,
  Theorem B(ii)–(iii), Remark B2, and Prop. 0(c) (26 per q, 50 over pairs).
- *Fixes applied:*
  - **(major)** the F2 assembly step (keep the min with K log h on the "1" part of the count);
  - the size of the q^{o(1)} loss (exp(O(log q/log log q)), sharp for this route; so the threshold is asymptotic);
  - the condition Y > 4H²;
  - the reference B(iv) → B(iii);
  - the dependency list;
  - the moment-type disambiguation;
  - averaged-form consistency in §3 and §7;
  - two §0 wordings;
  - the line "uniform κ is WLOG for Pos∀".

**GPT-6.1-sol adversarial round (`codex/cov38b_2.out`): ACCEPT WITH FIXES.** It recomputed the summary JSON from all 24
E1, 260 E2 and 3 E3 records, reproduced two complete E1 histograms, and found no deviation from the amended
pre-registration.
- *Fixes applied:*
  - **(major)** §5.1 claimed the data show the target true "with room to spare". At the prescribed c = 3/8 height,
    the toy sizes have 10⁵–10⁶ uncovered pairs. Corrected in §0 and §5.1.
  - **(major)** the "valid for every a" definition excluded S₀ itself (a = ℓ gives A = 0). Theorem B is now stated for
    certificates nonpositive off S_out, and prime-input-only certificates are flagged as not covered.
  - **(major)** the literature conclusions ("no tool", "every tool uniform in a is pointwise") narrowed to "none of
    the estimates checked here". Phase-retaining bilinear bounds exist but are too weak.
  - "iff" → sharp leading-exponent threshold, and only fixed c > 0 (log-scale improvements not excluded);
  - the F2 zero-term and assembly;
  - Irving's journal; Kerr–Shparlinski–Wu–Xi added to §6; 1/800 confirmed;
  - attribution of the signed Hölder step;
  - moment wording (2r points, ≥ 52 over pairs, Λ_q defined, literature-status label);
  - TP distinct-moduli accounting and the J_ℒ = 1 qualification;
  - B2 (Fejér), B′ (precise conclusion), M (common bound, support), S (modulus range, population vs residual);
  - reporting: 10 named predictions, 13 tests, slopes 0.80–0.85, Wiener points outside k₀ ≥ H′³, E2 grid assertion gap.
- *Its additional findings (adopted):*
  - the survivor family has infinitely many coprime members (lower convergents), with ≍ H² distinct values;
  - a scale-free covering statement would be false (a prime ℓ ≍ log q times [1, q−1] covers trivially), so the
    short-band constraint is essential;
  - Matomäki–Teräväinen's ternary prime products do not apply.

## 9. Reproduction

```
PYTHONPATH=. .venv/bin/python -m pytest -q tests/test_cov38b.py           # 13 tests, ~2 s
PYTHONPATH=. .venv/bin/python -m src.cov38b.experiments e1                 # ~2 min, < 1 GB
PYTHONPATH=. .venv/bin/python -m src.cov38b.experiments e2                 # ~25 min, < 1 GB (FFT length 1e7)
PYTHONPATH=. .venv/bin/python -m src.cov38b.experiments e3                 # ~3 min (8 processes)
PYTHONPATH=. .venv/bin/python -m src.cov38b.analyze
```
