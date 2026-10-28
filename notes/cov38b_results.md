# cov38b: can Bourgain–Garaev-type bounds for inverses of band primes push the deterministic worst case below N^{1/5}?

Agent `cov38b`, phase 30, 2026-10-28 (conn_alg shortlist #2).
- Pre-registration: `notes/cov38b_prereg.md` (commit ed36d22; addendum A in ed46052, both before any run).
- Code: `src/cov38b/{design,fourier,experiments,analyze}.py` (harness). Tests: `tests/test_cov38b.py` (11 pass).
- Raw data: `experiments/raw/cov38b_e{1,2,3}.jsonl`; verdicts: `experiments/analysis/cov38b_summary.json`.
- GPT-6.1-sol: blind round `codex/cov38b_1.{md,out}`; adversarial round `codex/cov38b_2.{md,out}` (§8).
- Opus critics (read-only): one on proofs, one on literature (§8).
- Notation as in notes/coverage38_results.md and notes/dchworst_results.md. ‖x‖_q := |x mods q|.

## 0. Verdict (plain language)

**No. Nothing in the Bourgain–Garaev / Kloosterman-sum toolbox reaches the needed exponent, and the worst case stays at
Harvey–Hittmeir's N^{1/5}. No new deterministic record.** The reasons can now be stated as exact exponents.

- **There is no intermediate prize.** Inside the design family of Theorem U, any coverage exponent c ≤ 3/8 gives
  exactly N^{1/5} (Theorem W: min(1/5, (1−c)/(7/2−c))), and even a logarithmic gain needs the full power statement
  λ < J^{1/4}. So "some c < 3/8 with an explicit δ" does not exist in this framework (Prop. 0).
- **The target (§1, Prop. T).** For every prime q ≍ P and every prime p ∈ (q/2, q), some prime ℓ ≡ 3 (mod 4) in
  (Y, 2Y], Y = q^{2/(7−2c)}, and some B ≤ H/4 satisfy ‖pBℓ̄‖_q ≤ qH/(4Y), with H = Y^{1/2}J^{−c}, c > 3/8. At
  c = 3/8: Y = q^{8/25}, H = q^{1/25}, the arcs have length q^{18/25}, and the expected number of hits is
  λ = J^{1/4} ≈ q^{0.08}. Uniform form: the J dilates ℓ·S₀ of one Farey-arc set S₀ of measure ≈ q^{−6/25} cover F_q^*.
- **Every "linear" tool fails, with an exact threshold (Theorem B, proved).** A pointwise bound
  max_{h≠0}|Σ_ℓ e(hℓ̄/q)| ≤ J q^{−δ}, used through any Erdős–Turán-type (phase-discarding) certificate:
  - certifies the target **iff** δ > δ*(c) = (1+2c)/(7−2c) (up to q^{o(1)}), for the natural certificates; at c = 3/8,
    δ* = 7/25 = 0.28. No certificate of any shape (signed, any support) works below (1/2+3c)/(7−2c) = 13/50;
  - but Parseval caps every such bound at δ ≤ 1/(7−2c) = 4/25 = 0.16 (square-root cancellation).
  - **The gap is 3/25 in the exponent of q, beyond square-root cancellation.** Known results at length q^{0.32}
    give δ = 0 for primes (all prime-variable Kloosterman bounds start at q^{1/2+ε}) and a log saving for integers.
  - New ingredient: the Fourier ℓ¹ norm of a Farey-arc set with H′ denominators is ≥ H′/polylog
    (Lemma F, via an exact Freiman isomorphism to rationals and an lcm-sum energy bound). This means the union over
    the H′ dilates buys nothing in a phase-discarding argument: one arc does as well as all of them.
  - Multiplicative characters (GRH included) and a-adaptive phase discarding fail the same way (Lemma M, Prop. B′).
- **Bourgain–Garaev specifically (§2).** Their exponential-sum bounds fall under Theorem B whatever their strength.
  Their additive-energy input is already optimal here: for Y large, the inverses of band primes have the minimal
  possible additive energy 2J² − J (an exact identity). It still gives nothing, because the obstruction is pointwise
  in the dilate p, not a lack of energy control. Their covering theorems (sums of reciprocals) concern sumsets, while
  the target is a covering by multiplicative dilates.
- **Two-pronged algorithms (§4, Theorem TP, proved in the model).** "DCL coverage OR Lehman/Harvey windows" is itself a
  design of Theorem U's family. At budget N^{1/5−η}, the Lehman prong covers at most a fraction ≍ N^{−5η/4} of pairs
  (E3 measures the predicted G^{−5/4} law: slopes −1.27 to −1.31). So the DCL prong must still cover every
  Lehman-hard pair at λ < J^{1/4}. Barrier B3 is genuinely escaped (linear families are Lehman-easy), but the residual
  statement is just as strong, and for coprime-pair proofs it still contains a sparse survivor (x, ⌊αx⌋).
- **Experiments (pre-registered; all 12 predictions pass).**
  - E1: exhaustive worst case over all primes p ∈ (q/2, q) for 15 primes q up to 1.7·10⁹ (up to 40.6M pairs per q) at the c = 3/8
    scaling. The worst height h* lies in the independent-model band in 15/15 cases, and below the Thue height for
    every q ≥ 10⁷. The empirical exponent is c_emp = 0.11–0.14 at q ≈ 10⁸–10⁹, rising slowly, as the model predicts
    (the model's limit is 1/2).
  - E2: the Fourier ℓ¹ norm of Farey-arc minorants grows like H′^{0.81–0.85} (Lemma F's H′/polylog), and
    phase-discarding certificates appear only at H ≈ 9–14·√Y, against a true worst case at or below √Y.
- **Bottom line.** A sub-1/5 proof by this route needs a positivity mechanism for unions over moduli that is not a linear
  statistic: some joint geometry-of-numbers or arithmetic argument. None is known, even for a log-factor gain over
  Thue/Minkowski. Our assessment is that it is out of reach with current technology. This is an assessment, not a
  theorem: it is not a lower bound for factoring.

## 1. The target statement (Task 1)

**Setting.** Fix c ∈ (3/8, 1/2]. For n-bit N let P = 2^{n/2}, let Y = Y(n) be n-only with
Y = N^{1/(7−2c)+o(1)} = P^{2/(7−2c)+o(1)}, let 𝒫 be the J ≍ Y/log Y primes ℓ ≡ 3 (mod 4) in (Y, 2Y], and let
H = ⌊Y^{1/2}J^{−c}⌋. Write q for the larger prime of a pair (p, q) ∈ S_n, so q ∈ (P/√2, √2P) and p ∈ [q/2, q).

- **Cov(c)** (the pair form; Theorem W's hypothesis, paper Cor. (b)). Every (p, q) ∈ S_n has ℓ ∈ 𝒫 and a primitive
  (A, B), 1 ≤ |A|, B ≤ H, with ℓ | Bp − Aq.
- **Pos(c)** (reciprocal form). For every pair: ∃ ℓ ∈ 𝒫, 1 ≤ B ≤ H/4 with ‖pBℓ̄‖_q ≤ k₀ := ⌊qH/(4Y)⌋.
- **Pos∀(c)** (uniform form). For every prime q in range and every a ∈ F_q^*: ∃ ℓ, B as above with ‖aBℓ̄‖_q ≤ k₀.
  Equivalently **F_q^* = ∪_{ℓ∈𝒫} ℓ·S₀**, where S₀ = S(H/4, k₀) and S(H′, k₀) := {B̄k : 1 ≤ B ≤ H′, |k| ≤ k₀} is the
  Farey-arc set. Equivalently F_q^* ⊆ 𝒫·[−k₀, k₀]·[1, H/4]^{−1}.

**Proposition T (proved).**
- (i) Pos(c) ⇒ Cov(c), and Pos∀(c) ⇒ Pos(c).
- (ii) Cov(c) at height H ⇒ the reciprocal statement with (H, 2qH/Y) in place of (H/4, qH/(4Y)). So Cov and Pos are
  equivalent up to a constant factor in H, which changes c by o(1).
- (iii) **Size.** |S(H′, k₀)| = (12/π² + o(1))H′k₀ when k₀ ≥ H′³. The expected number of hits is
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
  ≥ cN^{1/5}. Equivalently, Theorem W's exponent is min(1/5, (1−c)/(7/2−c)), which equals 1/5 for all c ≤ 3/8.
  - So a weaker positivity statement (c ≤ 3/8) gives **no** N^{1/5−δ} through Theorem R.
  - Log gains: λ = J^{1/4}(log N)^{−a} buys at most (log N)^{4a/5} (U4), and it is still a power-strength statement.
- (c) **No averaged form suffices by itself.** Theorem R guarantees success only on covered pairs. An almost-all
  statement (all but E pairs) needs a second algorithm for E.
  - If that algorithm is in the family, Theorem U applies to the union design (§4), and the union must cover every
    pair at combined λ < J^{1/4}.
  - Outside the family, nothing is known that recognises or factors the uncovered set. inversesieve shows it carries
    no usable local structure, and forgetful containers fail.
  - Moments: an exceptional count |S_n|·λ^{−j} < 1 needs j > (7−2c)/(1−2c), which is 25 at c = 3/8. The pair box
    supports level ≤ |S_n| ≈ Y^{7−2c}, i.e. moments of order ≤ 6 (sieve2d's error analysis stops at 4–5).

## 2. Tool exponents (Task 2)

### 2.1 The phase-discarding (Erdős–Turán) certificate and its exact threshold

Fix q, the band 𝒫, J = |𝒫|, and K(h) = Σ_{ℓ∈𝒫} e(hℓ̄/q). A **certificate** is a real g on Z/q with g ≤ 0 off a set S
all of whose hits imply coverage, plus a number κ ≥ max_{h≠0}|K(h)|. Since
Σ_ℓ g(aℓ̄) = (J/q)ĝ(0) + (1/q)Σ_{h≠0} ĝ(h)K(ah), it certifies coverage of every a once

  J·Σg/q > κ·W(g),   W(g) := (1/q)Σ_{h≠0}|ĝ(h)|.     (★)

By Prop. T(ii), any valid S lies in the outer set S_out := S(H, 2qH/Y). For c > 1/4 and Y large, 8H⁴ < Y, which
is the Freiman condition below.

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
- The lower bound E ≥ (2/3)H′(2k₀)³ (the H′ single-denominator progressions) shows that F1 is sharp up to the polylog.

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
- (The signed Hölder step is GPT-6.1-sol's, from the blind round; we re-derived it.) ∎

The same argument with g = 1_A for an arbitrary A ⊆ S, and the 3-term representation count
max_x #{r₃ + r₄ − r₂ = x} ≤ C k₀²H′(log)^8, gives the same k₀√H′ ratio bound for indicators.

**Theorem B (exact threshold for phase-discarding certificates; proved).** Let c ∈ (1/4, 1/2], Y ≥ Y₀, and κ = Jq^{−δ}.
- **(i) Sufficiency.** If max_{h≠0}|K(h)| < J(k₀+1)/q with k₀ = ⌊qH/(4Y)⌋, then Pos∀(c) holds. In exponents, this
  holds as soon as δ > δ*(c) := (1+2c)/(7−2c) + o(1).
  - *Proof.* Take G = Σ_{B≤H/4} Fej_B, with Fej_B(B̄k) = 1 − |k|/(k₀+1) the Fejér triangle on the arc B̄[−k₀, k₀].
  - Ĝ ≥ 0, so W(G) ≤ G(0) = H/4. Also ΣG = (H/4)(k₀+1), and supp G ⊆ S₀. Apply (★) and Prop. T(i).
  - Note that **H′ cancels**: the single arc B = 1 does exactly as well. Phase discarding cannot use the union over
    dilates.
- **(ii) Necessity, nonnegative certificates of positive mass.** If g ≥ 0 is supported in S_out with Σg ≥ θ|S_out|
  for a fixed θ > 0, then (★) forces δ ≥ δ*(c) − o(1).
  - *Proof.* By Lemma F(a), W(g) ≥ cθ³H/polylog − Σg/q. So (★) needs κ ≤ J|S_out|/(qW) ≲ J k₀/q·polylog.
- **(iii) Necessity, any certificate.** (★) with any real g ≤ 0 off S_out forces
  δ ≥ δ*(c) − (1−2c)/(2(7−2c)) − o(1) = (1/2 + 3c)/(7−2c) − o(1).
  - *Proof.* Lemma F(b) with H′ = H, k₀ = 2qH/Y.
- **(iv) Impossibility.** All of (ii) and (iii) exceed the Parseval cap 1/(7−2c) (Lemma P) as soon as c > 1/6. At the
  target c = 3/8 the numbers are: sufficiency/necessity δ* = 7/25 = 0.28; any certificate ≥ 13/50 = 0.26; Parseval
  cap 4/25 = 0.16. **No pointwise bound of any quality, known or conjectural, certifies the target through (★).** ∎

*Remarks.*
- coverage38 Prop. L(b) obtained "c ≤ 1/4" from W ≥ 1/2 − μ. Theorem B identifies the true price: W ≍ H′ for
  full-mass shapes, so the required pointwise saving is exactly δ*(c).
- E2 measures W(g_tri)/H′ = 0.41–0.65 and slope 0.81–0.85 in H′ (P2b).

**Proposition B′ (adaptive phase discarding with exact magnitudes; proved).** Suppose 16Y³ < q. Then
1/ℓ₁ + 1/ℓ₂ = 1/ℓ₃ + 1/ℓ₄ (mod q) forces {ℓ₁, ℓ₂} = {ℓ₃, ℓ₄}: the numerators are < q, and for primes the reduced
fractions determine the pair. So the inverses of band primes have additive energy exactly 2J² − J.
- This gives Σ_x|K(x)|⁴ = q(2J² − J) and Σ_x|K(x)|² = qJ.
- Hölder then gives mean_{x≠0}|K(x)| ≥ (1 − o(1))√(J/2).
- For any fixed g, averaging over a gives E_a Σ_{h≠0}|ĝ(h)||K(ah)| = q·W(g)·mean|K|. Paley–Zygmund applied to
  Z(a) = Σ|ĝ(h)||K(ah)|, with E[Z²] ≤ (qW)²·mean|K|², gives: **for at least 1/8 of all a ∈ F_q^***, the
  a-specific phase-discarded error is ≥ (1/2)(1 − o(1))√(J/2)·W(g).
- So even using the exact values |K(ah)| for each a, the uniform statement Pos∀ needs the same δ*(c).
- For the prime-only Pos, the averaging covers a positive proportion of residues but cannot by itself say where the
  primes p ∈ (q/2, q) fall.
- This is where Bourgain–Garaev-type energy information enters, and it is already optimal: the energy is exactly
  minimal.

**Lemma M (multiplicative characters, GRH included; proved).** Expand the hit count with characters mod q:
N(a) = (1/(q−1))Σ_χ χ̄(a)S_𝒫(χ)F̂(χ), where F̂(χ) = Σ_{|k|≤k₀}Σ_{B≤H′}χ(k)χ̄(B).
- Parseval over characters gives max_{χ≠χ₀}|S_𝒫(χ)| ≥ √(J/2).
- (1/(q−1))Σ_{χ≠χ₀}|F̂(χ)| ≥ #{(k, B) : k = B} − O(λ/J) ≥ H′ − 1.
- Phase discarding therefore needs λ ≳ √J·H′, i.e. J^{1/2−c} ≳ Y^{1/2}. With J ≍ Y/log Y, this fails for every
  c ≥ 0.
- With arbitrary nonnegative weights in place of F, it still needs λ ≳ √J, i.e. c ≤ 1/4.
- GRH gives |S_𝒫(χ)| ≪ √Y log²(qY), at the Parseval level, so GRH used this way gives c ≤ 1/4 at best.

**Lemma S (scale; proved in the model).** Under WIN, a group with modulus M has π_g ≤ C(log P)/M (coverage38 §1).
- A design of cost X has ≤ X groups.
- So the groups with M > Z carry λ_{>Z} ≤ C(log P)X/Z ≤ 1/log P for Z = X(log P)².
- Hence all but a fraction 1/log P of the pairs must be covered by moduli ≤ X·polylog ≤ q^{2/5+o(1)}.
- **The moduli that carry the coverage are always below q^{1/2}**, where every known prime-variable Kloosterman
  bound mod a prime is trivial.

### 2.2 The table (Task 2's question: "the exact exponent each tool gives")

δ is the pointwise saving in max_{h≠0}|K(h)| ≤ J q^{−δ} at the length the design needs (q^{0.32} at c = 3/8). "Certified
c" is the best coverage exponent the tool yields for Pos(c). Thue/Minkowski (c = 0, exponent 1/5) is always available
and is the baseline.

| tool (references in §6) | applies at length q^{0.32}? | δ it gives | certified c | reaches c > 3/8? |
|---|---|---|---|---|
| Thue / Dirichlet / Minkowski, one modulus (geometry, no analysis) | yes | — | **0** (exponent 1/5; Harvey's) | no |
| Prime-variable Kloosterman sums mod a prime: Bourgain 2005; Garaev 2010; Fouvry–Shparlinski 2011 (q^{3/4} for general moduli); Baker 2012; Korolev 2019 | **no**: nontrivial only from q^{1/2+ε}, and Lemma S puts the moduli below q^{2/5} | 0 | none (nothing beyond Thue) | no |
| Irving (averaged over moduli q) | no: one fixed q here | — | none | no |
| Integer-variable incomplete Kloosterman: Karatsuba; Bourgain–Garaev 2014 Thm 16; Korolev 2016 (composite moduli are allowed in Theorem R) | yes | log saving | none (δ = o(1) < δ*) | no |
| Bourgain's multilinear Kloosterman sums (moduli = products of k primes each ≥ q^ε) | yes | tiny power δ(ε, k) | none (δ ≪ δ* = 0.28) | no |
| Heath-Brown identity + Weil/Deligne completion | Type I sums of length < q^{1/2}: completion costs √q·log q > trivial | 0 | none | no |
| Bilinear Kloosterman fractions (Bourgain–Garaev; Shkredov; Kerr–Shparlinski et al.) as Type II input | possibly a small δ (UNVERIFIED ranges) | < Parseval cap 0.16 | none (Theorem B) | no |
| Hooley R* / square-root cancellation (conjectural) | — | 4/25 (the Parseval cap) | **none** (Theorem B: needs 7/25) | no |
| GRH, multiplicative characters | yes | Parseval level | ≤ 1/4 with arbitrary weights, ≤ 0 with the natural one (Lemma M) | no |
| Bourgain–Garaev additive energy / sum-product for {ℓ̄} | yes | energy exactly minimal (Prop. B′) | none (B′) | no |
| Large sieve / k-th moments (D′, sieve2d) | yes | almost-all only | worst case needs moment order > 25; box allows ≤ 6 | no |

**The exact gap (Task 3, "isolate").**
- *Linear-statistic route.* Positivity needs δ > (1+2c)/(7−2c) = 7/25. The most any bound can give is 1/(7−2c) = 4/25.
  The deficit is 2c/(7−2c) = 3/25 in the exponent of q. Equivalently, the hit count needs an error o(λ) = o(J^{1/4})
  on a sum of J terms: a saving of J^{3/4} over trivial, J^{1/4} beyond square root.
- *Known.* The saving is 0 for primes at q^{0.32}, and (log q)^{1−o(1)} for integers.
- *What would suffice instead.* A positivity mechanism that is not a linear statistic. The only one we know is
  single-modulus geometry (Thue/Minkowski), and B2 shows it stops at c = 0. No joint-geometry statement for unions
  over ℓ is known, even with a log-factor gain over Thue.

### 2.3 Bourgain–Garaev in detail (the shortlist's first question)

conn_alg asks: "Do Bourgain–Garaev-type multiplicative-energy bounds for {ℓ̄ : ℓ prime in a band} extend from length
q^{1/2+ε} to q^{0.32} for positivity in a union of arcs?" Three separate inputs answer it.
- **Exponential-sum bounds** (Bourgain 2005; Bourgain–Garaev 2014, multilinear Kloosterman with x_i ≤ q^ε). Every such
  bound is a pointwise linear statistic, so Theorem B applies whatever its strength. Even square-root cancellation
  (δ = 4/25) falls short of δ* = 7/25.
- **Energy / sum-product bounds** for reciprocals. At our scale the energy of {ℓ̄} is already exactly minimal
  (Prop. B′), so better energy bounds cannot help. Prop. B′ also shows that energy information, used through
  averaging over a, controls only a positive proportion of dilates. The worst case asks for all of them.
- **Covering theorems** ("every residue is a sum of k reciprocals of integers ≤ q^ε"; Bourgain–Garaev,
  Shparlinski). These cover F_q by additive k-fold sumsets of {x̄}. The target asks to cover F_q^* by the
  multiplicative dilates ℓ·S₀, i.e. by the product set 𝒫·[−k₀, k₀]·[1, H′]^{−1}.
  - For product sets the known all-classes ranges are the modular-hyperbola ones: xy ≡ a needs XY ≫ q^{3/2+o(1)}
    (Weil). The trivial character bound needs |A||B| > q².
  - The target's product has size λq = q^{1+2/25}, far below both.
  - Thue reaches product size q only because one factor, [1, H′]^{−1}·[−k₀, k₀] with H′k₀ ≥ q, is a full Dirichlet box.
- **Answer: no.** The obstruction is not a lack of distribution information. It is that positivity at
  λ < J^{1/4} is a statement below the fluctuation scale of any linear statistic (Theorem B), uniformly in the dilate.

## 3. Prove or isolate (Task 3)

- **Nothing is proved toward c > 3/8**, nor toward any c > 0 for all pairs, nor toward a log-factor union gain over
  Thue.
- What is proved:
  - Prop. T, the precise target and its equivalent forms;
  - Prop. 0, no intermediate regime and no usable averaged form in the model;
  - Theorem B, Lemmas F0, F1, F, P, M and S, and Prop. B′: the exact exponents and gap of all linear-statistic tools.
- **Dependency list for Theorem B.**
  - Lemma F0 (elementary).
  - Lemma F1 (elementary divisor sums: Σ_{M≤x}τ(M)⁴/M ≪ (log x)^{16}).
  - Hölder; Parseval.
  - Prop. T(i)–(ii) (coverage38 §3.1 with the corrected centre).
  - Fejér's kernel (Ĝ ≥ 0).
  - The size relations of Theorem W (Y = q^{2/(7−2c)+o(1)}, J ≍ Y/log Y by PNT in progressions mod 4).
  - Nothing else; in particular no unproved hypothesis.
- **Isolated gap.** See §2.2: 7/25 needed against at most 4/25 possible for linear statistics, and 0 known at
  q^{0.32}. A proof must be a non-linear positivity argument for unions over moduli.

## 4. Two-pronged algorithms (Task 4)

**Theorem TP (proved in the coverage38 model; Lehman prong per Prop. D(a)).** Let 𝒜 run DCL-wc with a band design 𝒟
and a Lehman/Harvey prong ℒ (stationary or full windows at its own moduli). Each prong is N-uniform, both are charged
per group, and WIN holds.
- (i) **Success.** 𝒜 factors every pair covered by 𝒟 (Theorem R) or by ℒ (Prop. D(a), or Theorem R at ℒ's moduli).
- (ii) **Same threshold.** 𝒟 ∪ ℒ is one design of the family. Theorem U gives X ≥ cN^{1/5}(λ⁴/J)^{1/5} with
  λ = λ_𝒟 + λ_ℒ and J = J_𝒟 + J_ℒ. No two-pronged combination of prongs from the family has a weaker coverage
  threshold than λ < J^{1/4}.
- (iii) **The Lehman prong covers few pairs.** Theorem U applied to ℒ alone gives λ_ℒ ≤ C(X/P^{2/5})^{5/4}J_ℒ^{1/4}.
  At X = N^{1/5−η} with J_ℒ = 1 (Lehman), the covered fraction is ≤ CN^{−5η/4}.
  - The direct count agrees: Pr[L* ≤ X] ≈ Σ_{h(v)≤√X} X√(P/h)/(hP) ≍ X^{5/4}/√P.
  - E3 measures slopes −1.31, −1.28, −1.27 against −5/4.
- (iv) **Residual statement.** 𝒜 handles every pair if **Res(η)** holds: every pair with L*(p, q) > N^{1/5−η} is
  covered by 𝒟. Here λ_𝒟 ≤ CJ_𝒟^{1/4}N^{−5η/4}, by Theorem U for the whole algorithm. This is sufficient. It is
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
  - p − 1 (Ward–Maestro's positive proportion below 1/6), Hittmeir/Harvey, and class-group methods: none has a
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
pairs retains B3's sparse survivor. The dichotomy changes the quantifier, not the difficulty.

## 5. Experiments (pre-registered; `notes/cov38b_prereg.md`, addendum A)

All 12 frozen predictions pass. The analysis script checks that every grid is complete.

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
- The data say the target is true with room to spare, consistent with K(θ) heuristics. They say nothing about
  provability.

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
- Korolev, Math. Notes 100 (2016) and 106 (2019).
- Irving, "Average bounds for Kloosterman sums over primes", Forum Math. 27 (2015), arXiv 1301.6372.
- Karatsuba's method; Shkredov, "Modular hyperbolas and bilinear forms of Kloosterman sums" (2019/2021);
  Kerr–Shkredov–Shparlinski–Zaharescu, "Energy bounds for modular roots and their applications" (2020s).
- Heath-Brown's identity (Canad. J. Math. 34, 1982); Weil/Deligne.
- Modular hyperbolas: Shparlinski's survey; XY ≫ q^{3/2+o(1)} for all classes.
- Prior project notes: coverage38 (Theorem U, Prop. L, §3), dchworst (Theorems R, W; B1–B4), inversesieve, sieve2d,
  avgtime.

## 7. Candidate card

| field | content |
|---|---|
| mechanism | Pos(c): positivity of inverses of band primes in unions of Farey arcs, via Bourgain–Garaev / Kloosterman / character tools |
| payoff if true | deterministic worst case N^{(1−c)/(7/2−c)} < N^{1/5} for c > 3/8 |
| status | **refuted for every linear-statistic (phase-discarding) certificate, with the exact threshold δ*(c) = (1+2c)/(7−2c) against the Parseval cap 1/(7−2c) (Theorem B, proved)**; no averaged/intermediate form helps (Prop. 0); two-pronged designs face the same threshold (Theorem TP) |
| experimentally supported | the target is true at toy sizes with the Poisson extreme (E1); Wiener norm ≍ H′ (E2); Lehman prong ≍ G^{−5/4} (E3) |
| sharpest obstruction | a positivity statement at λ = J^{1/4} is below the fluctuation scale √J of every linear statistic; the union over dilates is invisible to phase discarding (H′ cancels in Theorem B(i)) |
| what would revive it | a non-linear positivity mechanism for unions of lattices Λ_ℓ = Z(p, q) + ℓZ² (joint geometry of numbers), even a log-factor union gain over Thue as a first test |

## 8. Reviews

(to be filled)

## 9. Reproduction

```
PYTHONPATH=. .venv/bin/python -m pytest -q tests/test_cov38b.py           # 11 tests, ~2 s
PYTHONPATH=. .venv/bin/python -m src.cov38b.experiments e1                 # ~2 min, < 1 GB
PYTHONPATH=. .venv/bin/python -m src.cov38b.experiments e2                 # ~25 min, < 1 GB (FFT length 1e7)
PYTHONPATH=. .venv/bin/python -m src.cov38b.experiments e3                 # ~3 min (8 processes)
PYTHONPATH=. .venv/bin/python -m src.cov38b.analyze
```
