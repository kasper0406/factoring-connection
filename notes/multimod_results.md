# multimod: several moduli and multi-base searches beyond cycloharvey's N^{1/5} bound

Agent `multimod`, 2026-10-18. Pre-registration: `notes/multimod_prereg.md` (commit `81dd28e`,
before any pre-registered run). Code: `src/multimod/`. Tests: `tests/test_multimod.py`.
Raw data: `experiments/raw/multimod_*.jsonl`. Summary: `experiments/analysis/multimod_summary.json`.

## 0. Verdict (plain language)

**The gap is open, not closed. With several moduli, the model allows going below N^{1/5}.**
Unconditionally the worst-case exponent of the enlarged model lies in [1/6, 1/5]. It is
1/6 for almost all inputs (proved), and 1/6 in the worst case under a divisor-cover
conjecture.

- **Lower bound (proved).** With any number of moduli, one base and sumset baby/giant
  steps, every covering search costs ≥ c·N^{1/6} (Thm A). There is no M^{o(1)} loss. Two
  refinements:
  - Targets of degree ≥ 2 never beat N^{1/5}, however many moduli are used (Thm A′).
  - With J distinct moduli the bound is N^{1/5}·(J·N^{o(1)})^{−1/5}. Going below 1/5
    therefore needs N^{Ω(1)} moduli (Thm B).
- **Scheme (DCL, "divisor-cover Lehman").**
  - It uses Lehman/Harvey–Hittmeir targets (Bp + Aq − c)/M² with small fixed weights
    v = (A,B).
  - The moduli M range over all primes near N^{1/6}. No class is guessed: targets exist
    exactly for the roots of Br² ≡ AN (mod M).
  - The true p is hit iff some small v has Bp − Aq divisible by one of the primes.
  - One shared baby list serves every modulus via γ = β^{∏M²}.
- **Upper bounds for DCL (search model only).**
  - **Proved (Thm D):** cost N^{1/6+O(ε)} for all but an O(N^{−ε}) fraction of balanced
    semiprimes, via the arithmetic large sieve.
  - **Conjectured:** for *every* N (DCH, §3.6). Under DCH the enlarged-model exponent is
    exactly 1/6.
- **Guess-neutrality.** CRT does not break it.
  - **Proved (Thm C):** class-lattice indices multiply over coprime moduli, as φ does.
  - **Heuristic (cycloharvey's volume heuristic):** schemes covering all classes of their
    moduli therefore stay neutral. The rigorous single-composite-modulus statement is
    cycloharvey's LB, neutral up to M^{−o(1)}.
  - What escapes is *selective* covering. Each p needs only one modulus whose class
    holds a short weight vector. Across many coprime moduli these "lucky" classes cover
    almost everything (Thm D) and, conjecturally, everything (DCH).
- **Multi-base.**
  - **Proved for independent bases (Thm E):** they reduce coordinatewise to the 1-D bounds.
  - **Proved, trivially:** bases that are known powers of one generator are a 1-D search.
  - Mixed known-relation lattices: **not covered**.
- **Experiments (pre-registered; §6).**
  - Model counts: DCL at 95% success (R1) has slope 0.171 over 56–96 bits, against HH's
    0.196. It is cheaper than HH from 56 bits on: 0.45× at 96 bits.
  - Sampled success: 94–95% at 56–96 bits, against the predicted ≈95%.
  - The full-coverage regime (R2) has 0 uncovered candidates exhaustively at 32–48 bits
    and 100% sampled success at 56–96 bits. Its finite-size slope is 0.196, no better
    than HH's, because λ ∝ log N.
  - Operational runs: 64/64 factored in each regime, and the implementation list equals
    the model count in every list-building run.
  - **Caveat on the verdict.** The pre-registered verdict is **SUPPORTED** only after a
    classifier bug was fixed (disclosed in §6.4). With the original classifier, rule (d)
    failed at 56 and 64 bits.
- **Scope (important).**
  - Everything is in the *collision-search model* on the balanced range
    p ∈ [√(N/2), √N].
  - Unbalanced ratios N^{o(1)} < q/p < N^{1/3} are **not** handled at N^{1/6} by DCL
    (§3.9).
  - Two costs outside the search are unsolved: a base γ = β^L of adequate order (§3.8),
    and small-factor removal (Strassen, N^{1/6}).
  - **This is not a new deterministic factoring bound.** A deterministic worst-case
    N^{1/6} would need DCH, an order lemma for β^L, and the unbalanced range.

## 1. Models

- **𝔐₁ (cycloharvey §§2, 3.4).**
  - Targets are Laurent polynomials g ∈ Λ_I with I ⊆ [−d,d]∖{0}, or g ∈ L_{r,M} for
    classes modulo **one** modulus M.
  - Windows are integer intervals.
  - Coverage must be proved by the identity for every integer candidate p in
    [√(N/c), √N]: if p were a divisor, g(p) would be an integer in the window.
  - Cost is X ≥ #targets with Σ W_τ ≤ X² (shared or per-run baby lists).
- **𝔐⁺ (enlarged).** The same, except every target τ carries its own modulus M_τ ≥ 1,
  with gcd(M_τ, N) = 1 (Lemma K′ needs it), and its own class r_τ. Both are chosen freely:
  nested, coprime, any number. The cost hypotheses are unchanged:
  #targets ≤ X and Σ W_τ ≤ X².
  - These hold for per-modulus baby lists, a shared list, or any mix, because
    Σ W_τ ≤ Σ m_{j(τ)}⌈W_τ/m_{j(τ)}⌉ ≤ X·max m ≤ X².
  - The model is therefore *generous*: it does not charge the distinct-base problem that
    codex r10 warns about. DCL solves that problem anyway (§2.3).
- **𝔐^k (multi-base).** There are k bases α_i; a target is a handle ∏α_i^{F_i(p,q)} with a box
  window ∏[U_i, U_i + W_i). Baby set B ⊂ Z^k, giant sets G_τ, and coverage needs F(p) ∈
  B + G_τ, i.e. B + G_τ ⊇ box_τ, the whole box and not only the curve {F(p)}.
  Cost ≥ |B| + Σ|G_τ|.

## 2. The scheme: divisor-cover Lehman (DCL)

### 2.1 Class lock (proved; HH Prop. 3.3)

Let M be a prime with M ∤ 2ABN, and let r be a root of Br² ≡ AN (mod M). Write s = N·r^{−1}
(mod M²; the inverse is taken mod M², as codex points out). Then for every integer
p ≡ r (mod M) with q = N/p ∈ Z:

  Bp + Aq ≡ c := Br + As (mod M²).

*Proof.* Write p = r + Mx. Then q ≡ s − sMx·r^{−1} (mod M²), so
Bp + Aq ≡ Br + As + Mx(B − AN r^{−2}) = Br + As + Mx·(Br² − AN)/r² ≡ Br + As (mod M²). ∎

So g = (Bp + Aq − c)/M² ∈ Z. Over the balanced range, F = Bp + AN/p lies in [F_lo, F_hi]:
F is convex, with its minimum 2√(ABN) when A/B ∈ [1/2, 1] and an endpoint value otherwise.
The whole-range window is therefore g ∈ [⌈(F_lo−c)/M²⌉, ⌊(F_hi−c)/M²⌋]. Its width is
W ≤ (F_hi − F_lo)/M² + 1 = O(H√N/M² + 1), with H = max(A,B).

The exact check is `dcl.verify_windows`: rational arithmetic, convexity, both class
endpoints and the class members around the stationary point. It is run on every
experimental design and records 0 violations.

### 2.2 Coverage is a sieve

With whole-range windows, target (v, M, r) covers exactly the integers p ≡ r (mod M) in the
range. Since Bp² − AN = p(Bp − Aq), the condition p ≡ r for some root r is the same as
M | Bp² − AN, which for an actual divisor is M | Bp − Aq. So:

  **p is covered ⇔ some chosen (v, M) has M | Bp² − AN** (equivalently M | Bp − Aq).

The algorithm never guesses. Its target set is determined by N, and failure happens
only when no chosen (v, M) divides.

### 2.3 One baby list for all moduli (proved)

Let L = ∏_{M used} M² and γ = β^L. Since p ≡ 1 and q ≡ N (mod p − 1),

  (β^{L/M²})^{B + AN − c} ≡ β^{(L/M²)(Bp + Aq − c)} = β^{L g} = γ^{g} (mod p).

Codex verified the identity B + AN − (Bp + Aq) = (p − 1)(Aq − B).

- **Shared list.** Every target is a giant step (β^{L/M²})^{B+AN−c}·γ^{−g_lo−jm} against
  the single baby list γ^i.
- **Answer to the r10 warning.** The baby lists of different moduli coincide as group
  elements because the handle, not the base, is raised to L/M².
- **Cost of the powers.** All β^{L/M²} come from a product tree in the exponent
  (`operational.all_but_one_powers`) in O(log L·log #moduli) = Õ(Y) multiplications.
  Every handle then costs O(log N).

### 2.4 Cost and design

- **Cost.** m + Σ_τ ⌈W_τ/m⌉, cycloharvey's counting.
- **λ.** Each (v, M) has ρ_v(M) ∈ {0, 2} targets, and covers a fraction ρ/M of the
  candidates. Write λ = Σ ρ_v(M)/M for the expected cover multiplicity.
- **Efficiency.** With m given, a pair (v, M) with M ≳ (H√N/m)^{1/2} needs one giant step
  per target, and its efficiency ρ/M decreases in M.
- **The greedy design** (`dcl.design`) takes (v, M) pairs in decreasing order of
  ρ/(M⌈W/m⌉) until λ reaches its target. `best_design` then optimises m.
- **Heuristic size.** With M ≈ Y, cost ≈ m + λY·(1 + o(1)) and m ≈ √N/Y², so

  cost ≈ 3·(λ/2)^{2/3}·N^{1/6}·(polylog from weight heights).

### 2.5 Heuristic coverage (Poisson; conjectured)

Treat the events "M | Bp − Aq" as independent across (v, M), with probabilities
ρ_v(M)/M. Then a fraction ≈ e^{−λ} of the candidates is uncovered.
- **R1.** λ = 3, success ≈ 95%.
- **Full coverage.** Needs e^{−λ}·#candidates < 1, so λ ≈ ln √N.

The exhaustive check (§6, E-C) finds coverage somewhat *better* than Poisson.
- R1 leaves 1.3–4.3% uncovered, against Poisson's ≈5%.
- The smallest λ giving full coverage, λ*, grows like ln √N, as the Poisson tail predicts.
  It is 5.5–8 at n = 32–36 and 11–15.5 at n = 44–48, against ln √N + 2 = 13.1–18.6.

## 3. Theory

Notation: P = √N, and t = p/P ∈ J = [c^{−1/2}, 1]. #candidates ≍ P.

### 3.1 Per-target coverage (proved; from cycloharvey Lemmas C, K′)

Let τ be a target with modulus M (M = 1 when nothing is guessed), support I, |I| = K,
e = max|k| ≤ d, window width W.

- **(i) Degree e = 1, K = 2.** Lemma K′ (with K = 2) lets us write g = (Aq + Bp − c)/(C_d M²)
  with integers A, B over the **fixed** denominator C_d M². Define H := max(|A|,|B|) =
  C_d M²‖b‖_∞/P ≥ 1. This normalisation is also the one Thm B's count uses (critic,
  finding 4).
  Lemma C gives μ_τ ≤ C(W/‖b‖)^{1/2}, so

    #covered integers ≤ Pμ_τ/M + C_I ≤ C·(PW/H)^{1/2} + C.        (1)

  **The M-dependence cancels exactly.** That cancellation is the per-target form of
  guess-neutrality.
- **(ii) e = 1, K = 1** (support {1} or {−1}). Here ‖b‖ ≥ P/(C M), so the coverage is
  ≤ C·W + C ≤ C(PW)^{1/2} + C when W ≤ P. If some W > P, then X ≥ √W > N^{1/4} already.
- **(iii) e ≥ 2.** Lemma K′ gives ‖b‖ ≥ P^{e}/(C_d M^{K}), so

    #covered ≤ C·M·P^{−e/K}W^{1/K}·P/M + C = C·P^{1−e/K}W^{1/K} + C.   (2)

  Again M cancels.

Targets with the same (M, class, vector) are merged and Lemma C is applied to the union of
their windows. Lemma C only needs the measure, so the merged width is ≤ Σ W. "Distinct
targets ≤ X" therefore refers to distinct triples.
- The +C rounding term is paid once per *original* window, i.e. a preimage made of ≤ C_I
  intervals.
- The total is still ≤ CX, because every window costs at least one giant step.

### 3.2 Theorem A (proved): 𝔐⁺ costs ≥ c·N^{1/6}

**Statement.** Fix d and c > 1. Any 𝔐⁺ covering family of degree ≤ d has X ≥ c_{d,c}·N^{1/6}.

*Proof.*
- Coverage must reach c′P.
- By (1) and (ii), the e = 1 targets cover at most
  C P^{1/2} Σ_τ W_τ^{1/2} + CX ≤ C P^{1/2} X^{1/2}(Σ W)^{1/2} + CX ≤ C P^{1/2} X^{3/2} + CX
  (Cauchy–Schwarz).
- By (2) and Hölder over at most X targets of each of the finitely many (e, K) classes,
  the e ≥ 2 targets cover at most C Σ_{(e,K)} P^{1−e/K} X^{1+1/K} + CX.
- Put X = δP^{2/5} with δ ≤ 1. The (e,K) term is then
  ≤ C δ^{1+1/K} P^{1 − (5e − 2K − 2)/(5K)} ≤ C δ P, because 5e ≥ 2K + 2 when K ≤ 2e and e ≥ 2.
- So for X ≤ δP^{2/5} with δ small, the e ≥ 2 part covers ≤ c′P/2. The e = 1 part must then
  give C P^{1/2}X^{3/2} ≥ c′P/2 − CX, i.e. X ≥ c P^{1/3} = c N^{1/6}.
- If instead X > δP^{2/5}, then X ≥ δ N^{1/5} ≥ N^{1/6} for large N. ∎

**Theorem A′ (proved).** If the e = 1 targets of an 𝔐⁺ family cover fewer than c′P/2
candidates, then X ≥ c N^{1/5}, with any moduli. This is the second branch of the proof.

**Tightness.** Equality in (1) with H = O(1) and W ≈ P/M² means a target covers its *whole*
class across the range, P/M ≈ (PW)^{1/2} candidates. That is exactly what DCL's targets do
(M ≈ N^{1/6}, W ≈ N^{1/6}). So A is tight up to the polylog/overlap factor, if coverage can
be completed (Thm D, DCH).

### 3.3 Theorem B (proved): the bound as a function of the number of moduli

**Statement.** Suppose the e = 1 targets use J distinct moduli, all ≤ P (targets with M > P
cover ≤ 1 + C candidates each and are absorbed in the +CX term). Then

  X ≥ c·N^{1/5}·(J·N^{o(1)})^{−1/5}.

*Proof.*
- **Counting.** By cycloharvey §3.4 (e = 1 branch), for each modulus M and each (A,B), with
  H measured over the fixed denominator C_d M² as in §3.1(i), the number of classes on
  which (Aq + Bp − c)/(C_d M²) is integer-valued is ≤ O(1)^{ω(M)}·gcd(A,B,M).
  - The critic brute-forced composite M ≤ 105 with denominators M² and 2M². The count is
    always ≤ 4^{ω(M)}·gcd(A,B,M).
  - Measured over the *reduced* denominator the count fails: M = 25, D = 125. So the number of
  (class, vector) pairs with H ≤ h is ≤ M^{o(1)} h². Over J moduli it is ≤ J' h², with
  J' = J·N^{o(1)}.
- **Sum of 1/H.** For the ≤ X distinct triples, Σ_τ 1/H_τ ≤ C(J'X)^{1/2}: the extremal case
  takes the X smallest heights, so h_max ≍ (X/J')^{1/2} and Σ ≍ J' h_max.
- **Coverage.** Cauchy–Schwarz on (1) bounds the e = 1 coverage by
  C P^{1/2}(Σ W)^{1/2}(Σ 1/H)^{1/2} + CX ≤ C P^{1/2} X (J'X)^{1/4} + CX.
- **K = 1 targets.** Each covers ≤ CW + C (§3.1(ii)), so ≤ CX² + CX = o(P) in total when
  X ≤ δP^{2/5}.
- **Conclusion.** Theorem A's (e,K) computation handles the e ≥ 2 part. What remains is
  X^{5/4} ≥ c P^{1/2}J'^{−1/4}. ∎

**Corollary.**
- J = 1 recovers cycloharvey's LB.
- Exponent 1/5 − ε needs J ≥ N^{5ε−o(1)} moduli.
- J ≤ X gives Theorem A again, up to N^{o(1)}.

### 3.4 Theorem C (proved): CRT composes guess-neutrality multiplicatively

**Statement.** Let M = M₁M₂ with gcd(M₁, M₂) = 1, let r be a class mod M, and let
r_i = r mod M_i. Then

  L_{r,M} = L_{r₁,M₁} + L_{r₂,M₂}  and  [L_{r,M} : Λ_I] = [L_{r₁,M₁} : Λ_I]·[L_{r₂,M₂} : Λ_I].

*Proof.*
- **Decomposition.** Elements of L_{r,M} lie in M^{−∞}Λ_I = Λ_I ⊗ Z[1/M]. Partial fractions
  give M^{−∞}Λ/Λ = M₁^{−∞}Λ/Λ ⊕ M₂^{−∞}Λ/Λ, so g = g₁ + g₂ uniquely mod Λ.
- **Each part is integer-valued.** Suppose g(p) − c ∈ Z for all p ≡ r (mod M). The
  fractional part of g₁(p) lies in M₁^{−∞}Z/Z and depends only on p mod M₁^{a} for a large
  (with q = N p^{−1} mod M₁^{a}). By CRT, every class mod M₁^{a} lying over r₁ contains
  integers ≡ r₂ (mod M₂).
- On such integers, g₂(p) has M₂-power denominator. The fractional part of g(p) − c is 0, so
  the M₁-part g₁(p) − c₁ ≡ 0. Hence g₁ ∈ L_{r₁,M₁}, and likewise g₂ ∈ L_{r₂,M₂}.
- **Converse.** If g₁ ∈ L_{r₁,M₁} and g₂ ∈ L_{r₂,M₂}, then g₁ + g₂ is integer-valued up to a
  constant on the class r. So the sum is direct modulo Λ, and the indices multiply. ∎

Experimental check: `experiments crt` compares a direct HNF computation of
[L_{r,M₁M₂} : Λ] with the product of the two indices for I ∈ {(−1,1), (−2,−1,1)} and
(M₁,M₂) ∈ {(5,7), (3,7), (5,11)} (§6, E-CRT).

**Consequence (heuristic: it rests on cycloharvey's volume heuristic and the HNF-measured
κ = K(K+1)/2; the rigorous single-modulus statement is cycloharvey's LB, neutral up to
M^{−o(1)}).** Neutrality is
the identity "index M^{K(K+1)/2} against φ(M) runs". Both sides are multiplicative over
coprime parts, so within that heuristic a scheme covering **all** φ(M) classes of a composite M with a
shared baby list is exactly neutral, whatever its CRT structure.
- This includes CRT sign-splitting (signsplit's result) and the user's Gray-code bucketing:
  odd-m_i bucketing costs more than it saves, and only the 2-adic r₀ ↦ r₀ + 4 identity gives
  a constant.
- *CRT does not break neutrality; it composes multiplicatively and stays neutral.*

**What does break it.** DCL does not cover all classes of any modulus. For a fixed small v
only the ρ_v(M) ∈ {0, 2} classes with a root are used; that is a 2/M fraction, but each is
reached by a vector of height O(1) instead of the typical √M.
- These "lucky" classes from N^{1/6+o(1)} coprime moduli are independent by CRT. They union
  to almost everything (Thm D, proved), and conjecturally to everything (DCH).
- In Thm B's terms, many moduli supply J·h² cheap (class, vector) pairs, where one modulus
  supplies only h².

### 3.5 Theorem D (proved, in the search model): average-case N^{1/6+ε}

**Statement.** Fix ε ∈ (0, 1/10). For n-bit balanced semiprimes N = pq (primes p < q ≤ 2p),
take the following design. **H₀, Y and m depend on n only**, so the sifted classes for a
fixed q do not move with p (critic, finding 7).
- V = {(A,B) : gcd(A,B) = 1, A, B ≤ H₀} with H₀ = ⌈2^{nε/2}·n⌉;
- 𝒫 = the primes in (Y, 2Y] with Y = ⌈2^{n/6}H₀^{−1/3}⌉;
- every root r of Br² ≡ AN (mod ℓ), for each v ∈ V and ℓ ∈ 𝒫;
- m = ⌈C H₀ 2^{n/2}/Y²⌉.

Its model cost is O(H₀^{5/3}N^{1/6}) = N^{1/6+5ε/6+o(1)}. The fraction of balanced semiprimes
whose p is covered by no target is O(N^{−ε}).

*Proof.*
- **Cost.** At most 2|V||𝒫| targets with ⌈W/m⌉ = 1, by §2.1 with W ≤ C H₀√N/Y² + 1. So the
  cost is m + 2|V||𝒫| ≤ C H₀√N/Y² + 4H₀² Y/log Y = O(H₀^{5/3} N^{1/6}).
- **Reduction to a sieve.** Fix a prime q. The admissible p form an interval I_q of length
  ≍ √N (p ∈ [q/2, q] with pq of n bits). By §2.2, p is covered iff ∃(v, ℓ) with
  ℓ | Bp − Aq. For ℓ ∈ 𝒫 we have ℓ ∤ q (ℓ < q), so this is p ≡ Aq·B^{−1} (mod ℓ).
- **Distinct classes.** For ℓ > 2H₀² the slopes A/B are distinct mod ℓ, since
  |AB' − A'B| < H₀² < ℓ. So the uncovered p avoid ω(ℓ) = |V| residue classes modulo every
  ℓ ∈ 𝒫.
- **Large sieve.** The arithmetic large sieve (Montgomery 1968; Montgomery–Vaughan 1973)
  gives

    #{p ∈ I_q uncovered} ≤ (|I_q| + Q²)/J,  J ≥ Σ_{ℓ∈𝒫} |V|/(ℓ − |V|) ≥ |V|·Σ_{ℓ∈𝒫}1/ℓ ≫ H₀²/log N,

  with Q = 2Y. Here Q² = O(N^{1/3}) = o(|I_q|).
- **Normalising.** The primes in I_q number ≍ √N/log N, so the uncovered fraction for this
  q is O(log N·log N/H₀²) = O(N^{−ε}). For the few boundary q with short I_q, summing the
  large-sieve bound over q makes the Σ_q Q² terms O(N^{1/3}·√N) = o(#pairs). Summing over q
  gives the claim. ∎

The critic checked the cost and sieve steps:
- ⌈W/m⌉ = 1, since ℓ > Y gives span/ℓ² < m;
- ℓ ∤ 2ABNq;
- the classes are distinct;
- "p ≡ Aq/B (mod ℓ)" is exactly "the root p mod ℓ of Br² ≡ AN is in the design for N = pq".

**Remark (polylog version).** H₀ = log² N in the same proof gives cost N^{1/6}·polylog with
failure fraction O(1/log² N).

Remarks.
- Only the search is covered. Success also needs γ = β^L to have order mod p (or mod N) larger
  than the exponent range compared, about m + W, or one-sided detection otherwise (§3.8).
- For a *fixed* N, the same sieve in p (classes = roots of Br² ≡ AN) bounds the uncovered
  candidates by ≈ #cand/J_N. But J_N involves Σ_v (ABN/ℓ), and lower-bounding it needs
  quadratic-character cancellation over the short range ℓ ≈ N^{1/6}. That is UNVERIFIED in
  general; it would follow on average over ℓ from Heath-Brown's quadratic large sieve.
  Either way it gives "almost all candidates", not "all".

### 3.6 Conjecture DCH (divisor-cover hypothesis) and the worst case

**DCH(κ).** For all large N coprime to 6 and every integer p ∈ [√(N/2), √N], some coprime
A, B ≤ (log N)^κ and some prime ℓ ∈ (N^{1/6}, N^{1/6}(log N)^κ] satisfy ℓ | Bp² − AN.

- Under DCH, DCL is a deterministic covering in 𝔐⁺ with cost N^{1/6+o(1)}. With Thm A, the
  𝔐⁺ exponent is then **exactly 1/6**.
- **Status: conjectured.**
  - Experimentally supported at 32–48 bits: zero uncovered candidates among all ≈10⁵–5·10⁶
    integer candidates per instance at λ = ln√N + 2 (§6).
  - The large sieve proves only the "almost all" version.
- **Relation to Umans–Wang (codex check).** UW's Strong Prefactored Divisor Conjecture asks
  for structured S, T whose differences hit multiples of every d ≤ n. DCH asks that each
  *candidate* p make one of the forms Bp² − AN divisible by a prime from a window. Neither
  implication was found, and both are "every element of a long range is hit by a
  structured divisor family" statements.
  - He–Sahai's obstruction concerns UW's single-AP version of a different covering
    statement (every d ≤ n dividing a difference). DCH is a sieve statement: each (v, ℓ)
    removes ρ ≤ 2 N-dependent quadratic residue classes mod a prime ℓ from a window
    around N^{1/6}. We know no transfer of their argument to it, but we have not proved
    that none exists.

### 3.7 Theorem E: multi-base searches

**Setting.** In 𝔐^k, let Λ_known ⊂ Z^k be the relations ∏α_i^{x_i} ≡ 1 provable mod N for
all p, i.e. the identity-model relations. Coverage must hold in G = Z^k/Λ_known, because two
collisions mod p differing by an unknown relation are not identity-provable.

**Statement.**
- (a) **Proved.** The sumset bound |B|·|G_τ| ≥ |image of box_τ in G| holds in any abelian
  group. So Σ_τ |image| ≤ X² with a shared B, and #targets ≤ X.
- (b) **Proved, trivially.** If all bases are known powers of one β (G ≅ Z), as with DCL's γ
  and β^{L/M²}, the search is a 1-D search with target Σ e_i F_i ∈ L, and 𝔐⁺ applies.
- (c) **Proved for independent bases (Λ_known = 0).**
  - Each α_i^{F_i(p,q)} must be identity-valid on its own, so each non-constant F_i lies in
    Λ ∪ L_{r,M}, and Lemma K′ applies to it.
  - Every coordinate must lie in its window, so μ_τ ≤ min_i μ(F_i, W_i) ≤ μ(F_{i*}, vol_τ),
    using W_i ≤ ∏_j W_j = vol_τ.
  - With Σ vol_τ ≤ X² (box ⊆ B + G_τ, injective in Z^k), the coverage bounds (1)–(2) and
    Thms A–B hold verbatim. The same holds when Λ_known splits the bases into groups, each
    group a set of known powers of one generator: apply (b) inside each group and (c)
    across the groups.
- (d) **Not covered (GAP, critic finding 2).** General Λ_known ≠ 0 with free rank k′ ≥ 2 and
  non-split relations.
  - A single F_i need then only be valid modulo Λ_known, so Lemma K′ need not apply
    coordinatewise.
  - Projected windows are sums of APs rather than intervals, so Lemma C's per-interval
    constant is not controlled.
  - We know no mechanism that would exploit this, but the reduction is **not proved**
    there.

**Cross-difference.** A test is multipoint-evaluable with one product tree iff it is
separable in G: E(τ, j; i) = f_τ(j) − h(i). Baby sets mixing σ-powers or bases with
p-dependent coefficients are non-separable, as cycloharvey §3.5 showed. Splitting them
into separable pieces gives back case (c).
- **Verdict.** For independent or grouped bases, k-dimensional exponents do not help
  beyond what several moduli already give: **refuted for this model**. General
  known-relation lattices are open (d).

### 3.8 Costs outside the search (not solved here)

- **Order of γ = β^L.** ord_p(γ) = ord_p(β)/gcd(ord_p(β), L). It becomes small when p − 1 is
  divisible by two or more of the moduli, which happens with probability
  ≈ (Σ_{ℓ∈𝒫}1/ℓ)²/2 ≍ (log u/log Y)².
  - Then the baby list repeats mod p. That reveals p by gcd (class PM1) unless it is
    simultaneously small mod q.
  - **The deterministic order lemma for L-th powers is open (UNVERIFIED).**
    Oznovich–Volk give order ≥ D in D^{1/2+o(1)}, but for β, not β^L.
- **Small factors.** Removing factors ≤ N^{1/3} by Pollard–Strassen costs N^{1/6+o(1)}
  (r10). That is at the same exponent, so it is no obstruction to 1/6 but blocks going
  lower.

### 3.9 Unbalanced range (limitation; heuristic)

For q/p ≈ R ≫ 1, the weights need B/A ≈ R, and F's span grows to ≈ √(RN).
- With per-band baby lists the band costs ≈ (RN)^{1/6}. With the shared list it costs
  R^{1/4}N^{1/6}.
- Summing over dyadic R up to N^{1/3} gives N^{2/9} or N^{1/4}. That is worse than Harvey,
  whose Lehman pairs ab ≤ r cover all ratios at once.
- So DCL as specified is a balanced-range method. A uniform-ratio version is
  **not yet tested**.

## 4. The three orchestrator inputs

- **Gray-code / Stage-3 leaves (user).**
  - Agreed, and it is the content of Thm C. Any all-classes scheme keeps
    #leaves ≍ φ(M)·m₀ = N^{1/5}; only the generation cost changes.
  - DCL reduces the number of **distinct leaves** to N^{1/6+o(1)}.
  - Each DCL leaf covers a whole residue class across the range, P/M ≈ (PW)^{1/2}
    candidates, the maximum allowed by (1). An HH leaf covers about P^{1/3}W^{2/3}·N^{−o(1)}.
- **r10 warnings.**
  - (3) Separate bases: solved by γ = β^L (§2.3). The baby lists coincide as group elements.
  - (4) Base generation: §3.8, open for L-th powers.
  - (2) UW / He–Sahai: §3.6.
- **Chirp / quadratic-in-k leaves (user).** Leaves Y_k = Y₀γ₁^kγ₂^{k²} cost O(1) to
  generate, but Stage 3 is the product tree ∏_k(X − Y_k) evaluated at the baby points. Its
  cost is Õ(#leaves + #babies) however the leaves were made, so cheap generation saves only
  Stage 2.
  - Could chirp leaves cover more exponents per leaf? No. The covered exponent set is
    E = {e_k + j}, with |E| ≤ #leaves·#babies, so the sumset bound |B|·|G| ≥ W is unchanged.
  - Each window W_k is set by the unknown p inside a bucket, (√(ap) − √(bq))², not by
    centre truncation (orchestrator (ii)). So the number of leaves Σ_k⌈W_k/m⌉ is the same as
    without the chirp.
  - When leaves and babies share β, the double product over (leaves × babies) is literally
    the 1-D test over E. The chirp is a Stage-2 constant, as you expected.
  - **Refuted for this model** (as an exponent improvement).

## 5. Prior art (GPT-6.1-sol check, 2026-10-18; curated)

- **Ingredients (PROVED as prior art).**
  - Lehman 1974: multipliers.
  - McKee 1999: a modulus-squared restriction "if a prime m divides …", heuristic N^{1/4+ε}.
  - Hittmeir 2018: modular-hyperbola sieving over many small primes.
  - Harvey 2021: a shared collision search over multiplier families, and Alg. 4.1.
  - HH 2022 Prop. 3.3: the class lock, with a single modulus and lattice-chosen weights.
  - Hart 2012: one-line factoring, many multipliers, heuristic N^{1/3}.
- **Not found (UNVERIFIED that none exists).** No source with the union over many prime
  moduli near N^{1/6} with fixed small weights, nor with the cofactor-exponent shared baby
  list. We claim nothing beyond "no prior formulation found by us".
- **Context.** A heuristic N^{1/6} is not a general factoring milestone: ECM and QS are
  conjecturally subexponential. The interest is only within the deterministic
  collision-search line (Harvey/HH/UW).

## 6. Experiments (pre-registered; results)

All runs follow `notes/multimod_prereg.md`. The one deviation is the classifier correction
in §6.4. Summary: `experiments/analysis/multimod_summary.json`.

### 6.1 E-M: model counts (exact; N = 2^n + 2^{n/2} + 1)

Cost is list elements, m + Σ⌈W/m⌉.

| n | DCL R1 (λ=3, ≈95%) | DCL R2 (λ=ln√N+2, full) | HH (calibrated) | Harvey (calibrated) | R1/HH |
|---|---|---|---|---|---|
| 40 | 320 | 1 399 | 267 | 1 025 | 1.20 |
| 48 | 799 | 4 238 | 797 | 3 733 | 1.00 |
| 56 | 2 095 | 12 666 | 2 297 | 12 598 | 0.91 |
| 64 | 5 483 | 37 483 | 6 891 | 42 030 | 0.80 |
| 72 | 13 843 | 110 321 | 20 435 | 138 603 | 0.68 |
| 80 | 36 167 | 329 878 | 60 731 | 453 165 | 0.60 |
| 88 | 92 780 | 976 005 | 177 595 | 1 467 601 | 0.52 |
| 96 | 239 092 | 2 916 279 | 536 185 | 4 713 259 | 0.45 |

Slopes (least squares, n = 56–96): **R1 0.1707**, R2 0.1962, HH 0.1963, Harvey 0.2137.

The heuristic is cost ≈ 3(λ/2)^{2/3}N^{1/6}×(height factor). That gives R1 → 1/6
asymptotically. For R2, λ grows like ln√N, adding ≈ (2/3)/(n ln 2) ≈ 0.01–0.02 to the slope,
plus a growing weight height.

### 6.2 E-C: exhaustive coverage of all integer candidates (n = 32–48; 4 random + special N per size)

- **R2.** **0 uncovered** in all 25 instances, out of 1.6·10⁴ to 4.9·10⁶ candidates each.
  Window violations: 0.
- **R1.** Uncovered fraction 1.3–4.3%, below the Poisson prediction e^{−λ} in every
  instance. Examples: n = 48, 4.2–4.3% against 4.9%; n = 32, 1.3–2.4% against 4.9–5.4%.
- **λ\*** (smallest grid λ with full coverage), against ln√N + 2:

| n | λ* | ln√N + 2 |
|---|---|---|
| 32 | 5.5–7.5 | 13.1 |
| 36 | 7.0–8.0 | 14.5 |
| 40 | 9.0–9.5 | 15.9 |
| 44 | 11.0–12.0 | 17.2 |
| 48 | 12.0–15.5 | 18.6 |

  λ* grows by ≈ 0.4 per bit, close to the Poisson tail's ln(#candidates) slope of 0.35 per
  bit. Full coverage at λ* costs ≈3.5–4.5× the R1 cost at 48 bits.
- Experimentally supported: DCH holds for these instances with λ ≈ ln√N, at polylog extra
  cost.

### 6.3 E-S: sampled success of the true p (harness)

| n | instances | R1 success | R2 success | R1 median cost | R2 median cost |
|---|---|---|---|---|---|
| 56 | 2000 | 0.9535 | 1.000 | 1 997 | 11 770 |
| 64 | 2000 | 0.9465 | 1.000 | 5 149 | 35 193 |
| 80 | 1000 | 0.942 | 1.000 | 34 383 | 309 977 |
| 96 | 300 | 0.937 | 1.000 | 229 236 | 2 769 822 |

The slight downward drift of R1 (0.954 → 0.937) is within ≈2σ at 96 bits (σ ≈ 0.014). It is
still consistent with a mild size dependence of the overlap structure; **not yet tested**
beyond 96 bits.

### 6.4 E-O: operational runs (16 balanced semiprimes per size, full lists)

| n | R1 succ. | R1 classes (S/PM1/A) | R1 list = model | R2 succ. | R2 classes | R2 list = model | Harvey succ. | mean list cost R1 / R2 / Harvey | mean mults R1 / R2 / Harvey |
|---|---|---|---|---|---|---|---|---|---|
| 40 | 16/16 | 11/3/2 | 13/13 | 16/16 | 12/4/0 | 12/12 | 16/16 | 296 / 1 264 / 967 | 1.1e4 / 3.8e4 / 1.0e3 |
| 48 | 16/16 | 15/0/1 | 16/16 | 16/16 | 12/3/1 | 13/13 | 16/16 | 760 / 3 847 / 3 482 | 4.1e4 / 1.4e5 / 1.1e5 |
| 56 | 16/16 | 13/1/2 | 15/15 | 16/16 | 10/5/1 | 11/11 | 16/16 | 1 983 / 11 743 / 11 903 | 1.2e5 / 3.8e5 / 4.5e5 |
| 64 | 16/16 | 14/0/2 | 16/16 | 16/16 | 14/1/1 | 15/15 | 16/16 | 5 156 / 35 364 / 39 814 | 3.7e5 / 1.7e6 / 1.7e6 |

- **Factoring.** 128/128 DCL runs factored N.
- **List = model.** The implementation's giant list equals the model count in every run
  that built a list. Runs that stop at a PM1 gcd before the list build no list.
- **Wall time at 64 bits:** R1 0.07 s, R2 0.25 s, Harvey 0.34 s.
- **Multiplications.** These include the exponent product tree (§2.3) and one O(log N)
  powering per target. At 40 bits Harvey runs at r = 1 (plain BSGS), so its count is tiny.
- **Classifier correction (deviation, disclosed).** The pre-registered classifier checked
  the collision equation only for x = p.
  - DCL also finds q: target (A,B,M,r) collides mod q when q ≡ r (mod M) and
    (Bq + Ap − c)/M² falls in the window. That is the same mechanism with the factors'
    roles swapped, as in cycloharvey's classifier, which uses "the factor found".
  - With the original classifier, STRUCTURAL+PM1 was 12/16 at 56 bits and 13/16 at 64 bits
    for R1, so rule (d) **failed**.
  - With the corrected classifier it is 14/16 at both sizes.
  - Original records: `experiments/raw/multimod_oper_v1_origclassifier.jsonl`. The rerun used the
    same seeds, and its runs are identical apart from the class labels.
- **ACCIDENTAL hits are order-collapse hits** (post-hoc diagnostic, not pre-registered;
  `experiments diag` → `experiments/raw/multimod_oper_diag.jsonl`).
  - Every remaining ACCIDENTAL instance has a null expectation (m·#giants)/ord_p(γ) +
    (m·#giants)/ord_q(γ) of ≥ 1.5. For the four cleanest it is 1.5–160, computed from
    the exact ord_p(β^L).
  - The pre-registered null was not defined. These hits are collisions in a group whose
    order L has partly killed (§3.8), not DCL structure. They still factor N (one-sided).

### 6.5 E-CRT

Direct HNF indices [L_{r,M₁M₂} : Λ] equal the products [L_{r₁,M₁} : Λ][L_{r₂,M₂} : Λ] in
**20/20** cases: I = (−1,1) at (5,7), (3,7), (5,11), and I = (−2,−1,1) at (3,7), (5,7).
The independent critic also checked prime powers: (9,5), (4,7), (3,13).

### 6.6 Decision (pre-registered rule)

| rule | result |
|---|---|
| (a) R1 slope ≤ 0.185 | 0.1707 ✓ |
| (b) HH − R1 ≥ 0.015 | 0.0256 ✓ |
| (c) E-S R1 ≥ 0.93 everywhere | 0.937–0.954 ✓ |
| (d) list = model; STRUCTURAL+PM1 ≥ 14/16 | list ✓; class ✓ only after the classifier correction (✗ with the original code) |

**Verdict: SUPPORTED, with the disclosed classifier correction.** In the calibrated
list-element model, at fixed ≈95% success, DCL's exponent is below 1/5 (slope 0.171 against
HH 0.196 over 56–96 bits), and it is cheaper than HH from 56 bits on.

- **R2 (matched 100% success).** As predicted: full coverage everywhere tested, but the
  finite-size slope is 0.196, equal to HH. At ≤ 96 bits there is no advantage at matched
  100% success. R2 lists are 0.99× and 0.89× Harvey's at 56 and 64 bits, but 5.4× HH's model count
  at 96 bits.
- Matched-success comparisons against HH use HH's *model* count, because HH is not
  operationally implemented here.

## 7. Candidate card

- **Mechanism.** A Lehman/HH one-sided exponent collision. Weighted sums Bp + Aq are
  locked mod M² on the N-determined classes Br² ≡ AN (mod M), for **many prime moduli**
  M ≈ N^{1/6} and small fixed weights. One baby list γ^i (γ = β^{∏M²}) is shared by all
  moduli through cofactor-powered handles.
- **Representation.** (Z/N)^*. Targets (v, M, r). The design uses only N (Legendre data
  mod small public primes).
- **Factor feature.** g = (Bp + Aq − c)/M² in a whole-range window of size ≈ N^{1/6}.
  **Extraction.** Product tree / Bluestein gcd (Alg. 4.1), or exact match plus Lemma 3.1.
- **Outcome.**
  - Model lower bound N^{1/6} (proved).
  - DCL attains N^{1/6+o(1)} on average (proved, search model) and in the worst case under
    DCH (conjectured).
  - Experiments: slope 0.171 at 95% success; full coverage at λ ≈ ln√N exhaustively to 48
    bits.
- **Scoped obstruction.**
  - Below N^{1/6}, nothing in the one-base sumset model works (Thm A).
  - Below N^{1/5}, nothing works with N^{o(1)} moduli or with degree ≥ 2 targets
    (Thms A′, B).
- **What would escape / is missing.**
  - A proof of DCH, or a structured choice of (v, M) with provable full coverage.
  - An order lemma for β^L.
  - A ratio-uniform version (§3.9).
  - Going below 1/6 needs a non-sumset collision mechanism, and also faster small-factor
    removal than Strassen.

## 8. Reproduction

```
.venv/bin/python -m pytest -q tests/test_multimod.py
PYTHONPATH=. .venv/bin/python -m src.multimod.experiments model
PYTHONPATH=. .venv/bin/python -m src.multimod.experiments coverage
PYTHONPATH=. .venv/bin/python -m src.multimod.experiments mc        # needs model first
PYTHONPATH=. .venv/bin/python -m src.multimod.experiments oper
PYTHONPATH=. .venv/bin/python -m src.multimod.experiments crt
PYTHONPATH=. .venv/bin/python -m src.multimod.experiments diag      # post-hoc, after oper
PYTHONPATH=. .venv/bin/python -m src.multimod.experiments summarize
```
