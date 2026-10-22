# dsr: Hypothesis DSR (deterministic smooth representation): derandomised walks, certificates, compression

Agent `dsr`, phase 23, 2026-10-22. **Version 2.** It incorporates the GPT-6.1-sol adversarial round and both Opus
critics (§10). v1 is commit `3801973`.

**Files.**
- Pre-registration: `notes/dsr_prereg.md` (commit `409641a`, before any statistic was computed).
- Code:
  - `src/dsr/core.py`: operational; sees only N; deterministic.
  - `src/dsr/harness.py`: E1 diagnostic (discrete logs, float FFT); E2 seeded random controls.
  - `src/dsr/analyze.py`.
  - `src/dsr/posthoc_e1.py`, `src/dsr/posthoc_adv.py`: post hoc, labelled.
- Tests: `tests/test_dsr.py` (7 tests).
- Data: `experiments/raw/dsr_e1.jsonl`, `dsr_e2.jsonl`, `dsr_e1_posthoc.jsonl`, `dsr_posthoc_adv.jsonl`; summary
  `experiments/analysis/dsr_summary.json`.
- GPT-6.1-sol: blind round `codex/dsr_1.{md,out}`; adversarial round `codex/dsr_2.{md,out}`, full report
  `notes/dsr_codex2_report.md`.

**Labels.**
- **proved**; **proved (folklore-level / elementary)**; **conditional (GRH)** etc.; **open**;
  **observation** (a limitation of a proof template or instantiation, not a lower bound); **heuristic**;
  **experimentally supported**; **refuted**.
- Novelty: at most "no prior statement found by us"; the literature search is limited.

## Verdict (plain language)

**DSR was not proved. Nothing here brings it materially closer.** The work instead sharpens where it is stuck,
through four negative results and two reductions, all checked by three reviewers.

1. **Black-box derandomisation of the JMV walk does not give DSR [proved for the stated instantiations; not a lower
   bound for arithmetic-specific methods].**
   - Under GRH the random walk on Cay((Z/N)^×, primes ≤ B) hits the smooth set after (1 + o(1))ln N/ln B steps. That
     is ≈ log₂ N random bits, i.e. N^{1+o(1)} seeds rather than L^{O(1)} (Prop RW).
   - Any guarantee that holds for every density-δ accepting set (total-variation closeness of the endpoint law;
     Nisan/INW/NZ with arbitrary terminal predicates; AKS/IZ samplers) needs more than (1 − δ)φ(N) endpoints
     (Theorem SB(a)).
   - Character-wise generators (GKM Fourier shapes, small-bias spaces) need per-character error ≤ N^{−1/2+o(1)} for the
     canonical smooth set, because that set's Fourier-algebra norm is N^{1/2−o(1)} (Theorem FN). That costs Θ(log N)
     seed bits.
   - The single-accepting-state width-free results (HPV, BHPP) do not cover a δφ-element accepting set.
   - In short, GRH spectral bounds describe the *walk*. A proof of DSR must use arithmetic of the specific smooth set
     beyond density and Fourier bias (Theorem SB(b)).
2. **Smooth sets are Fourier-large [proved; folklore-level duality].** Theorem FN gives ‖1_{S_can}‖_A ≥ N^{1/2−o(1)}
   for every admissible N *unconditionally* (Pólya–Vinogradov), and ‖1_{S_Gauss}‖_A ≥ N^{1/4−o(1)} under GRH.
   - This settles zroute's speculative FC(c) conjecture for these two sets: Bach/Ankeny-type small-norm certificates
     are impossible for them and for their dense minorants.
   - It does not exclude other testable sets or non-norm arguments.
3. **Short multiplier boxes are provably non-uniform [proved, unconditional, elementary; no prior statement found by
   us].** Prop BOX-∞: for every b and every large h there is an admissible N < 2^h for which the DSR target 2^H mod N is
   a tiny non-smooth integer r⁸.
   - For that N, no positive multiplier v ≤ N^{1−o(1)} gives a smooth least residue.
   - No B-smooth rational multiplier of height ≤ N^{1/2−o(1)} gives a smooth Gauss split.
   - This **refutes** the natural missing lemma proposed in GPT-6.1-sol's blind round (k-fold prime-product
     correlation, uniform in N). It does not refute DSR: lifts, signed or large representatives, and large w escape.
4. **Almost-all versions are easy [conditional (GRH); folklore-level].** Prop AA: for every N, the explicit list of
   k-fold prime products works for all targets outside a set E_N with |E_N|/φ ≤ L^{−kb+1/(2b)+o(1)}. E_N contains every
   small non-smooth unit residue, at least (1 − o(1))φ/B^k of them. DSR needs every target q^H to avoid E_N, and by
   BOX-∞ they need not.
5. **DSR includes a prime-modulus problem [proved; folklore padding].** Lemma PT: DSR, even restricted to N whose
   prime factors all exceed L^C, yields the same smooth-representation constructor modulo every prime p. Run it on
   N = p·r with r a known small prime.
   - So DSR must solve the mod-p problem on inputs whose factorisation is supplied; knowing the factorisation does not
     trivialise it.
   - What Theorem DZ needs minimally (dominant rows, or PIN) is equivalent to deterministic factoring under ERH
     (Prop EQ-dom). DSR is a strictly more structured sufficient condition, and we know no reduction from factoring
     to DSR.
6. **Smooth compression (u ≤ N³) [proved reformulation].** Prop SC: it is a height-coupled translate-hitting problem
   with inverse multipliers d^{−1} (d ≤ N²) and B lifts. The N³ freedom drives GRH existence (GX); we see no way it
   helps a constructor.
7. **Conditional derandomisation [conditional; folklore-level].** GRH plus Impagliazzo–Wigderson circuit hardness
   implies DSR_{1/2}, hence ERH + IW ⇒ deterministic L[1/2] factoring (Prop IWD; also polyimpl DR).
8. **Closest analogue:** deterministic prime finding (Polymath4 / Tao–Croot–Helfgott 2012). It has the same shape: a
   dense, efficiently testable set, a trivial randomized algorithm, and a deterministic version open even under RH.

**What remains (§8).** A deterministic L^{O(1)} list of legal multiplier data per N, together with a proof that some
entry gives a smooth representation of q^H for *every* admissible N and every q ≤ B. BOX-∞ forces such a list, if it
uses canonical or Gauss tests, to contain wrapping integer multipliers (≥ N^{1−o(1)}) or inverses (canonical test), or
rationals of height ≥ N^{1/2−o(1)} (Gauss test); lifted and other tests are unconstrained. FN and SB show that
spectral and Fourier-norm arguments cannot supply the proof. By PT, the same is needed modulo primes. This is a
restricted sufficient route, not a characterisation of DSR.

**Probability (subjective, not calibrated by these results).** DSR provable from GRH with known techniques: ≲ 3%
(ours); ≈ 10% (codex blind).

---

## 1. Setting

- N is **admissible**: odd, not a prime power, no prime factor ≤ max(B, 3 ln² N). G = (Z/N)^×, φ = φ(N).
- L = L_N[1/2, 1] and B = B_N = ⌊L^b⌋. Here b is zroute's a, renamed so that a can denote the target.
  H = 8⌈log₂ N⌉. Targets a = q^H mod N, q ≤ B prime.
- **S_can** := {x ∈ G : least positive residue of x is B-smooth}. δ_can = Ψ(N − 1, B)/φ = L^{−1/(2b)+o(1)}.
- **S_Gauss**: x whose Gauss-reduced vector (u, v) (u ≡ xv) has both halves B-smooth. **S_bal** (experiments) is the
  subset with u² + v² < N; there the shortest vector is unique up to sign.
- Fourier: Ŝ(χ) = Σ_{s∈S}χ(s), 1_S = φ^{−1}Σ_χ Ŝ(χ)χ̄, bias β(S) = max_{χ≠χ₀}|Ŝ(χ)|, and
  ‖1_S‖′_A := φ^{−1}Σ_{χ≠χ₀}|Ŝ(χ)|. Plancherel: Σ_x f(x)\overline{g(x)} = φ^{−1}Σ_χ f̂(χ)\overline{ĝ(χ)}.
- **Translate hitting TH(S, T, W):** every a ∈ T has some w ∈ W with a·w ∈ S, i.e. T ⊆ S·W^{−1}.
- **What DD1 uses [proved, trivial].** For each q ≤ B, some row j·e_q + e(w) − e(u) + e(v) with
  j > Ω(u) + Ω(v) + Ω(w) and q^j·w·v ≡ u. So orbit targets q^{H+j} and multipliers are equally usable
  (`test_orbit_hit_gives_dominant_row`). Theorem DZ itself may need less (PIN; Prop EQ-dom).
- DSR is asymptotic ("N ≥ N₀"). Codex's example (N = 1009·1013, B = 2) shows outputs need not exist at tiny B.

## 2. Approach 1: deterministic search and derandomised walks

### 2.1 Explicit spectral data (GRH)

**Prop RW (the randomized constructor) [conditional (GRH); folklore-level: JMV 2009 + Fourier inversion].**
- **Spectral bound.** Let μ be uniform on the primes ≤ y, with y ≤ B (all are units), and λ := max_{χ≠χ₀}|μ̂(χ)|.
  - JMV's proof of Thm 1.1 (arXiv:0811.0647; checked by us and the literature critic) bounds the weighted prime
    character sums for every nonprincipal χ, induced ones included. Partial summation gives
    |Σ_{p≤y}χ(p)| ≪ √y·ln(Ny), so λ ≪ ln y·ln(Ny)/√y.
  - JMV's Cor 1.3 assumes y ≥ (log N)^{B₀} with B₀ > 2, which holds for y = L^b.
  - Codex derives a loose explicit form λ ≤ min(1, 10⁶√y ln³(Ny)/π(y)) from Ernvall-Hytönen–Palojärvi's
    character-level estimates. It is not JMV's constant, and we have not audited it.
- **Mixing.** For S of density δ and any start a,
  |Pr[a·p₁⋯p_t ∈ S] − δ| ≤ λ^t(φδ(1 − δ))^{1/2} (Parseval + Cauchy–Schwarz). This needs λ < 1.
- **Lengths** at y = B, ln(1/λ) ≥ (b/2 − o(1))√(ln N ln ln N):
  - Hitting with probability ≥ δ/2 needs t_hit = ⌈ln(4φ(1 − δ)/δ)/(2 ln(1/λ))⌉. For δ = N^{−o(1)} (true for S_can)
    this is (1/b + o(1))√(ln N/ln ln N).
  - Pointwise mixing (GX) needs t_∞ = ⌈ln(2φ)/ln(1/λ)⌉ ≈ 2t_hit.
- **Randomness and endpoints.** t_hit·log₂π(B) = (1 + o(1))log₂ N bits; 2 log₂ N at t_∞. Products are v ≤ N^{1+o(1)}
  at t_hit and N^{2+o(1)} at t_∞. Fixed-bit sampling uses a power-of-two prefix of the primes, which has the same
  spectral bound. An implementation uses certified bounds in place of the unknown φ and λ.
- **The DSR output.** u = least residue of a·v < N, v = p₁⋯p_t, w = 1, all B-smooth, Ω(v) = t ≪ log₂ N. Expected cost
  L^{b+1/(2b)+o(1)} per target (codex agrees).
- **Consequence.** Under GRH, DSR is a randomized L[1/2] search problem; the question is purely derandomisation. L^{c}
  seeds means seed length c√(ln N ln ln N)/ln 2 = o(log₂ N).

### 2.2 Theorem SB (spectral blindness)

**Theorem SB [proved; (a) trivial, found independently by codex].** Let G be finite abelian of order φ ≥ 16, and W ⊆ G
(or the support of a multiset).
- **(a)** If |W| ≤ φ − ⌈δφ⌉, some S ⊆ G ∖ W has |S| = ⌈δφ⌉. So no hitting guarantee can follow from hypotheses on
  (W, generators, spectrum) plus the density of S alone. The same holds for adaptive membership queries: answer "no"
  to every query.
- **(b)** If |W| ≤ φ/2, 0 < δ ≤ 1/4 and δφ ≥ 24 ln(8φ), there is S ⊆ G ∖ W with ||S| − δφ| ≤ (6δφ ln 8φ)^{1/2} and
  β(S) ≤ 2δ·β(1_W) + (12δφ ln 8φ)^{1/2}.

*Proof of (b).*
- Put δ′ := δφ/(φ − |W|) ≤ 2δ ≤ 1/2, and include each g ∉ W independently with probability δ′.
- For χ ≠ χ₀, E Ŝ(χ) = −δ′·1̂_W(χ).
- Re(Ŝ − EŜ), Im(Ŝ − EŜ) and |S| − δφ are sums of independent mean-zero terms of modulus ≤ 1, each with variance
  ≤ σ² := δφ.
- Bernstein with s := (6σ² ln 8φ)^{1/2} ≤ σ²/2 gives each failure probability ≤ 2(8φ)^{−2}.
- There are 2φ − 1 events, so the union bound is < 1. Hence some outcome has |Ŝ(χ) − EŜ(χ)| ≤ √2·s for all χ ≠ χ₀. ∎
  (All three reviewers re-derived this.)

*Corollary SB′ (Fourier sup-norm certificates).*
- Apply (b) with 2δ in place of δ, assuming δ ≤ 1/8 and 2δφ ≥ 24 ln 8φ. This gives |S| ≥ δφ, and for
  |W| ≤ (φ ln(8φ)/δ)^{1/2} also β(S) ≤ C(δφ ln 8φ)^{1/2}.
- Every set of that size has β ≥ ((|S|(φ − |S|))/(φ − 1))^{1/2} ≥ (δφ/2)^{1/2} (Parseval). So the constructed S is
  within O((ln φ)^{1/2}) of optimally Fourier-uniform.
- Hence no theorem "β(S) ≤ X, |S| ≥ δφ ⇒ W ∩ S ≠ ∅" with X ≥ C(δφ ln 8φ)^{1/2} holds for L^{O(1)}-size lists.
  (v1 claimed a factor 7, which was not justified at finite size. Corrected.)

*What SB does and does not say about derandomisers (rescoped after review).*
- **Covered by (a).** Instantiations whose correctness is stated for every accepting set of the width-φ walk program:
  - total-variation closeness of the endpoint law;
  - Nisan/INW/Nisan–Zuckerman with arbitrary terminal predicates;
  - AKS/IZ samplers from a uniform start.
  Any such guarantee needs > (1 − δ)φ distinct endpoints.
- **Character-wise generators (GKM Fourier shapes, small-bias spaces).** These are covered not by SB but by Theorem
  FN. The generic summation needs per-character error ε with ε‖1_S‖′_A < δ/2, and ‖1_{S_can}‖′_A ≥ N^{1/2−o(1)}
  forces log(1/ε) ≥ (1/2 − o(1))ln N.
- **HPV/BHPP** state results for one accepting state. Summing HPV's errors over δφ states needs ε = O(1/φ). BHPP's HSG
  does not transfer that way (codex).
- **KNP/De** have group-size dependence: KNP constants depending on the fixed group; De's seed is
  O((w⁸ + log 1/ε)log n) (literature critic, via Steinke).
- **Reingold/UES** are traversal results (length ≥ φ − 1).
- None of this is a lower bound against a method that uses the arithmetic of S_N.

### 2.3 The tools, one by one

| tool | what it saves | instantiation for DSR / failure point |
|---|---|---|
| AKS 1987; Impagliazzo–Zuckerman 1989 | randomness for *repeating* an experiment | still one uniform start of ≈ log₂ φ bits: N^{1+o(1)} seeds |
| Nisan 1992; INW 1994 | seed O(log ℓ·log(Wℓ/ε)) | W = φ, ℓ = O(log N): O(log N log log N) |
| Nisan–Zuckerman 1996 | poly coins in space S to O(S) | S ≈ log N: N^{O(1)} seeds |
| Cohen–Doron–Goldgraber 2026 (INW forward-backward, permutation programs) | O((log w + log 1/ε)log n) | log w = log φ: Θ(log N) |
| Hoza–Pyne–Vadhan 2021; Bogdanov–Hoza–Prakriya–Pyne 2022 | width-free, *one* accepting state | δφ accepting states: ε ≈ 1/φ for HPV; BHPP's width-w HSG keeps log w |
| Gopalan–Kane–Meka 2015 (Fourier shapes) | fools each endpoint character | needs ε < δ/(2‖1_S‖′_A) ≤ N^{−1/2+o(1)} by Theorem FN |
| KNP 2011; De 2011 (group products) | small bias for a fixed group | constants (KNP) or poly(w) (De) with w = |G| ≈ N |
| Rozenman–Vadhan 2005 | expansion at low degree | biases must reach ≈ φ^{−1/2}; support ≥ φ/(1 + φβ²) by Parseval |
| Reingold 2008; UES | deterministic traversal | length ≥ φ − 1; products far above N³ |
| Impagliazzo–Wigderson 1997 | seed O(log T) for size-T circuits | **works conditionally** (Prop IWD) |

**Prop IWD [conditional; folklore-level; codex agrees].**
- *Statement.* If some f ∈ E has circuit complexity ≥ 2^{εm} at all large lengths m, then under GRH DSR_{1/2}(b, c)
  holds for some c = c(b, ε).
- *Proof.*
  - Repeat the Prop RW experiment ⌈8/δ₀⌉ times, with the computable δ₀ := L^{−1/(2b)−ε′} ≤ δ for large N. The
    failure probability is ≤ e^{−4}, and the circuit size is T = L^{O(1)}, including prime-table generation and trial
    division.
  - The IW generator with seed O(log T) fools size-T circuits with error < 1/2.
  - Enumerate the T^{O(1)} seeds and verify each output.
- With Theorem DZ: ERH + IW ⇒ deterministic L[1/2, O(1)] factoring, consistent with polyimpl Theorem DR.
- We know of no result showing that a deterministic DSR algorithm would imply circuit lower bounds. (v1's stronger
  "it would not imply" is withdrawn.)

### 2.4 What would suffice (codex, checked), and what FN removes

- **(i) A Fourier-small minorant.** f ≤ 1_S with E f ≥ δ/2 and ‖f‖_A ≤ L^{O(1)}. Then k fixed prime steps suffice
  (cost L^{kb+O(1)} times the test).
  - For **S_can** and **S_Gauss** no such f exists. The proof of Theorem FN applies verbatim to such minorants: f is
    ≤ 0 on the disjoint translate set T, so ‖f‖′_A ≥ (δ/2)|T|/max|T̂|. This holds unconditionally for S_can and under
    GRH for S_Gauss.
  - Any useful f must be a minorant of a *different* testable set (e.g. tests that try lifts x̄ + jN).
- **(ii) A conditional-expectation estimator** for h_j(x) = Pr[x·p₁⋯p_j ∈ S] to additive error δ/(8t). GRH supplies
  it only for large j.
- **(iii)** IW-type hardness (Prop IWD).

### 2.5 Prop AA (almost all targets, every N)

**Prop AA [conditional (GRH); folklore-level: Parseval second moment].**
- *Statement.* Let W_k be the multiset of ordered k-tuples of primes ≤ y (size m = π(y)^k, products ≤ y^k), and S ⊆ G
  of density δ. Then E := {a ∈ G : aw ∉ S ∀w ∈ W_k} has |E| ≤ λ^{2k}φ(1 − δ)/δ.
- *At y = B and fixed k.* |E|/φ ≤ L^{−kb+1/(2b)+o(1)} for S_can. The enumeration has π(B)^k = L^{kb+o(1)} products,
  and trial-division testing adds a factor L^{b+o(1)}.
- *Proof.*
  - h(a) := #{w : aw ∈ S} = φ^{−1}Σ_χ χ̄(a)Ŝ(χ)\overline{Ŵ(χ)}, with Ŵ = (Σ_{p≤y}χ(p))^k; the χ₀ term is δm.
  - Σ_a|h(a) − δm|² = φ^{−1}Σ_{χ≠χ₀}|Ŝ|²|Ŵ|² ≤ (λ^k m)²φδ(1 − δ).
  - Each a ∈ E contributes (δm)². ∎
- *Comparison of bounds (not an asymptotic for |E|).* By BOX(i), E contains every unit a whose least residue
  r ≤ N/y^k is not B-smooth: (1 − o(1))φ/B^k elements. The upper bound exceeds this by at most L^{1/(2b)+o(1)}.
- *Contrast.* A random list of size m leaves φ(1 − δ)^m uncovered in expectation.
- *What DSR needs.* The L^{b+o(1)} targets q^H must *all* avoid E. Cardinality information certifies that only when
  the bound on |E| is < 1, i.e. for |W| ≳ φ/δ. (v1's §2.7 "exceptional-set" bullet confused one hit with all hits.
  Corrected.)

### 2.6 Prop BOX and BOX-∞

**Prop BOX [proved, trivial].** Let a ∈ G have least positive residue r, not B-smooth, and let ℓ_r > B be a prime
dividing r.
- **(i)** For every integer 1 ≤ w < N/r, a·w mod N = rw is not B-smooth.
- **(ii)** For B-smooth integers s, t ≥ 1 with (rs)² + t² < N, the Gauss-reduced vector of a·s·t^{−1} is ±(rs, t)/g,
  g = gcd(rs, t). Since g | t is B-smooth, ℓ_r ∤ g, so the first half rs/g is not B-smooth.
  - *Why the shortest vector is that.* The lattice {(u, v) : u ≡ a s t^{−1} v} contains (rs, t) and has determinant N.
    Two independent vectors of squared length < N would give 0 < |det| < N. So the primitive vector (rs, t)/g is the
    unique shortest vector up to sign.
- **(iii)** For the canonical test, inverse multipliers escape: r·2^{−1} ≡ (r + N)/2 ≈ N/2 (proof critic).

**Prop BOX-∞ (adversarial N in every size range) [proved, unconditional, elementary; no prior statement found by us;
the small-offset modulus N = 2^h − r is classical background (pseudo-Mersenne / Solinas)].**
- *Statement.* Fix b > 0, put B_* := ⌊L_{2^h}[1/2, b]⌋ ≥ B_N for N < 2^h, u_h := ⌈3 ln ln B_*⌉, and R := B_*^{u_h} = N^{o(1)}.
  For all large h there is an admissible N ∈ (2^{h−1}, 2^h) with 2^h ≡ r (mod N), r ∈ (B_*, R), r not B_*-smooth.
  Then H = 8⌈log₂ N⌉ = 8h and 2^H ≡ r⁸ (mod N), with r⁸ < N having a prime factor > B_N. Consequently:
  - for every integer 1 ≤ v ≤ N/r⁸ = N^{1−o(1)}, 2^H·v mod N = r⁸v ∉ S_can;
  - for B-smooth s, t with (r⁸s)² + t² < N (height ≤ N^{1/2−o(1)}), 2^H·s/t ∉ S_Gauss.
- *Proof.*
  - Let z := max(B_*, 3h² ln² 2) + 1, which exceeds max(B_N, 3 ln² N) for every N < 2^h. Sift
    A := (2^h − R, 2^h − B_*) ∩ Z by the primes < z.
  - **Fundamental lemma of the sieve** (Friedlander–Iwaniec, *Opera de Cribro*, Ch. 6; numbering not checked), with
    bounded weights at level D = R^{1/2}: the remainders are |r_d| ≤ 1 uniformly in the interval's position, and
    s = ln D/ln z = u_h/2 − o(1) → ∞. So S(A, z) = |A|V(z)(1 + o(1)) + O(R^{1/2}) ≥ (e^{−γ} − o(1))R/ln B_*.
  - **Remove primes.** By Brun–Titchmarsh (Montgomery–Vaughan 1973) there are ≤ 2R/ln R = o(R/ln B_*) of them.
  - **Remove prime powers.** There are O(h), at most one per exponent, since consecutive j-th powers near 2^h differ
    by ≫ 2^{h/2} > R.
  - **Remove N with r = 2^h − N B_*-smooth.** There are ≤ Ψ(R, B_*) ≪ R e^{−u_h/2} ≤ R(ln B_*)^{−3/2}
    (Tenenbaum III.5, Thm 1).
  - A remaining N is z-rough (odd; no prime factor ≤ max(B_N, 3 ln² N)), composite and not a prime power, so it is
    admissible. Also gcd(r, N) = gcd(2^h, N) = 1. Then apply Prop BOX. ∎
  - (All three reviewers checked the sieve steps.)
- *Consequences.*
  - **(1) Refutation [proved].** GPT-6.1-sol's proposed missing lemma (`codex/dsr_1.out` §6) is
    Σ_{p₁..p_k≤B}1_{S_can}(q^H p₁⋯p_k) ≥ η·π(B)^k·δ for every admissible N ≥ N₀ and every q ≤ B, with k = k(b) fixed.
    For q = 2 on BOX-∞'s N the sum is 0. A version restricted to odd q is not refuted by this construction. Shifts
    2^j·r⁸ stay small, so the mechanism is not specific to q = 2, but we have not constructed odd-q instances.
  - **(2) What is not refuted.** DSR itself allows signed u, |u|, |v| ≤ N³ (lifts r⁸v ± jN) and a separate w with
    Ω(w) ≤ log₂ N, of height up to B^{log₂N}. These escape BOX. So does deterministic factoring by other means (DD1's
    q = 2 row is one design, not a necessity).
  - **(3) Orbits.** Orbit lists q^{H+j} escape after O(log N) steps (r⁸2^j wraps). Wrapping is not a proof of a hit.
- *Experiments (post hoc, labelled; `src/dsr/posthoc_adv.py`).*
  - (i) In all 30 E2 instances, D_box's worst-case coverage satisfies m* ≥ Ψ(N/ℓ, B) + 1, with ℓ the least prime > B
    (ℓ ∤ N in all instances, since their factors exceed 2^9). The gap is 0–3 in the 20 uncensored instances. Codex
    recomputed ℓ's own first hit: it equals m* in all 20, so a = ℓ is the worst target there. For the 10 censored
    instances, Ψ + 1 > cap.
  - (ii) A seeded scan of 4000 balanced 24-bit semiprimes found 7 (N, q) pairs with q^H mod N = r ≤ N/2¹², r
    non-smooth. On those targets D_box first hits at Ψ(N/r, B) + 1..7 (2 685–8 038), against 3–76 for the random walk.
    This is toy-scale frequency only.

### 2.7 Second moments over the targets (template computations; limitations of these upper-bound methods only)

- **Over q.** Σ_{q≤B}|h(q^H) − δm|² ≤ λ^{2k}m²[π(B)|Ĝ[H]|δ(1 − δ) + C√B log²(NB)·‖1_S‖′_A²] (codex re-derived it).
  With ‖1_S‖′_A² ≤ δφ, certifying every q needs λ^{2k} ≲ δ/(√B φ), i.e. |W| ≥ φ^{1−o(1)}. q ↦ q^H need not be
  injective. FN shows ‖1_{S_can}‖′_A² ≥ N^{1−o(1)}, so this template cannot be rescued by a better norm bound.
- **Over an orbit q^{H+j}, j < J.** The pointwise route needs λ^k‖1_S‖′_A < δ, which is impossible for S_can by FN.
  A Parseval route needs J·m ≳ φ/δ over *distinct* residues.

## 3. Approach 2: Bach-style certificates and coset structure

- **Why Bach works [standard].**
  - A coset has ‖1_C‖_A = 1 and a subgroup complement 2(1 − 1/k), and GRH needs only the characters trivial on the
    subgroup.
  - Hitting a coset of index k needs primes up to ≈ ((k − 1)log N)² (Lamzouri–Li–Soundararajan 2015, Thm 1.4; Math.
    Comp. 84, 2391–2412, corrigendum 86 (2017) 2551–2554). So it is cheap only at small index.
- **Coset decomposition.**
  - The set of classes with *some* smooth representative ≤ N³ is all of G under GRH: one coset, with no
    representative supplied.
  - The testable sets are truncations, so they are not closed under multiplication. Codex's example: N = 1009·1013,
    B = 2, S_can = {2^j : j ≤ 19}, which contains no nontrivial coset.

**Theorem FN (smooth sets have large Fourier-algebra norm) [(i) proved unconditionally, folklore-level duality +
Pólya–Vinogradov; (ii) conditional (GRH). No prior statement of these lower bounds at x ≈ N found by us or the
literature critic. Closest work: Hardy–Xu 2026 (upper bounds on the same first moment for q ≥ x^{1+ε}) and Harper.]**
- **(i)** For every admissible N ≥ N₀(b): ‖1_{S_can}‖′_A ≥ N^{1/2−o(1)}. Since ‖1_S‖_A ≤ |S|^{1/2}, this is sharp up to
  N^{o(1)}.
- **(ii)** Under GRH for Dirichlet L-functions mod N: ‖1_{S_Gauss}‖′_A ≥ N^{1/4−o(1)}.
- Both bounds hold for every minorant f ≤ 1_S with E f ≥ δ/2, with δ replaced by δ/2.

*Proof.*
- **Duality.** If T ⊆ G is disjoint from S (or f ≤ 0 on T), Plancherel gives
  0 ≥ Σ_x f(x)1_T(x) = φ^{−1}[f̂(χ₀)|T| + Σ_{χ≠χ₀}f̂(χ)\overline{T̂(χ)}]. Hence
  ‖f‖′_A ≥ (E f)·|T|/max_{χ≠χ₀}|T̂(χ)|.
- **(i)** Choose a prime ℓ ∈ (B, 2B] with ℓ ∤ N. One exists for large N, since such primes outnumber ω(N); codex's
  N = 489733·18889…821 shows that the least prime > B can divide N.
  - Let T := {ℓm mod N : 1 ≤ m ≤ M, gcd(m, N) = 1}, with M := ⌊(N − 1)/ℓ⌋. These are distinct units, whose least
    residues ℓm are divisible by ℓ > B, so T ∩ S_can = ∅. |T| ≥ N^{1−o(1)}.
  - T̂(χ) = χ(ℓ)Σ_{m≤M}χ(m). For nonprincipal χ induced by χ* mod q* > 1, Möbius over d | N gives
    |Σ_{m≤M}χ(m)| ≤ τ(N)·max|Σ χ*| ≪ τ(N)√q* log q* = N^{1/2+o(1)} (Pólya–Vinogradov).
  - δ_can = N^{−o(1)}. ∎(i)
- **(ii)** Take μ uniform on primes ≤ B, λ ≤ B^{−1/2+o(1)} under GRH (Prop RW), ℓ as above, and
  k := ⌊ln(√N/(4B))/ln B⌋.
  - For every k-tuple, s = ℓp₁⋯p_k satisfies s² + 1 < N, so s·1 is its own Gauss split (BOX(ii)) and is not smooth.
  - Expanding E[f(ℓp₁⋯p_k)] = φ^{−1}Σ_χ f̂(χ)χ̄(ℓ)\overline{μ̂(χ)}^k ≤ 0 gives ‖f‖′_A ≥ (E f)λ^{−k}.
  - ln λ^{−k} ≥ (k ln B)/2 − O(k ln ln N) = (1/4)ln N − O_b(√(ln N ln ln N)), so λ^{−k} = N^{1/4−o(1)}. The loss is at
    L-scale, which is harmless here.
  - **δ_Gauss ≥ N^{−o(1)}** (codex's construction). Take disjoint prime sets I₁, I₂ ⊂ (B/2, B] of size ≍ B/ln B, and
    t := ⌊ln(√N/4)/ln B⌋. Let u, v range over products of t distinct primes from I₁ and from I₂ respectively.
    - Then u² + v² < N and gcd(u, v) = 1, so distinct pairs give distinct x = u/v (BOX(ii)).
    - There are ≥ ((|I₁|/t)^t)² = N^{1−o(1)} such pairs.
  - The same walk argument reproves (i) under GRH. ∎

*Meaning.*
- Bach/Ankeny-type certificates need ‖·‖_A = L^{O(1)} (zroute FC(c)). For S_can and its dense minorants they are
  impossible for every N, unconditionally; for S_Gauss, under GRH.
- In L¹ Fourier norm, S_can is as spread out as a random set of its size up to N^{o(1)}. The S_Gauss bound is not
  shown sharp.
- A bad translate set T with small character sums forces a large norm. This is the duality behind both FN and the BOX
  obstruction.
- FN does not exclude other testable sets, lifted tests, or arithmetic arguments that are not norm certificates.

**E1 (pre-registered; §7).**
- The L¹ Fourier norm of S_can and S_bal at toy size is *smaller* than for random sets: r_A = 0.55–0.63 and 0.32–0.35,
  so P1a ([0.7, 1.3]) **failed**. The maximum bias is 5–72× larger (P1c passed). The norm grows like |S|^{0.47}
  (P1b passed; random 0.51), consistent with FN.
- *Post hoc (`posthoc_e1.py`; suggestive only).* The largest coefficients sit at characters of *large* order, with
  mean Re χ(p) over the factor base 0.33–0.86. In the post hoc subset the factor base has 8–36 primes, with empirical
  λ = 0.57–0.99; the GRH bound is vacuous there.
- A non-expanding small factor base is a *possible* explanation for the sub-random constant. Codex points out that it
  is not established: λ = 0.57 is compatible with a good expander, and FN permits sub-random constants. We regard the
  P1a failure as **unexplained**. P1c's registered rationale (small-order characters) was wrong even though the
  direction held.

## 4. Approach 3: smooth compression

**Prop SC (height-coupled reformulation) [proved, trivial; corrected after codex's counterexample to v1].**
- *Statement.* For a ∈ G, a B-smooth u ≡ a (mod N) with 0 < u ≤ N³ exists iff there are a B-smooth d ≤ N² and
  0 ≤ j < B such that c := (least residue of a·d^{−1}) + jN is B-smooth and d·c ≤ N³.
- *Proof.*
  - (⇒) If u ≤ N², take d = u and c = 1. Otherwise multiply prime factors of u into d while d ≤ N². Then d > N²/B,
    c = u/d < NB, and c ≡ a d^{−1}.
  - (⇐) Take u = dc. ∎
- *Why the coupling matters.* v1 omitted d·c ≤ N³. Codex's counterexample: N = 1031·1033, B = 2, d = 2^{40},
  c = 2^{21} gives dc = 2^{61} > N³, and no 2-smooth u ≤ N³ is ≡ a.
- *Meaning.*
  - Smooth compression is translate hitting with inverse multipliers d^{−1} (≤ N² labels, at most φ distinct residues)
    and B lifts, coupled through the height constraint.
  - The N³ freedom is what GX's existence proof uses. We see no way it helps a constructor, but SC is a reformulation,
    not a proof that it cannot (codex).
- **Relation to Prop SI / Theorem J / sihard.**
  - Prop SI asks for L^{o(1)}-cost constructors at L[1/3] scale. Theorem J says single-integer Johnson-class
    certificate decoders yield on average no more than random tests. DSR at L[1/2] can afford random testing on
    average (δ = L^{−O(1)}), and is purely about worst-case determinism. Theorem J's average statement neither obstructs
    nor helps it.
  - Codex checked Boneh's CRT decoder with the full factor base: its agreement threshold ≈ √(2B ln N)/ln B exceeds the
    ≈ 3 ln N/ln B prime divisors from a dyadic interval that an integer ≤ N³ can have.
  - sihard's NP-hardness (Prop XS) concerns explicit prime lists, a different problem. Under GRH, DSR is a total search
    problem with poly-size verifiable witnesses, so NP-hardness under deterministic reductions would give NP = coNP
    (Megiddo–Papadimitriou; sihard Prop NC).

## 5. Structural results

**Prop EQ-dom [proved; ERH; folklore-level].**
- *Statement.* For fixed α, these are equivalent up to constants:
  - (a) deterministic L[α, O(1)] complete factoring;
  - (b) a deterministic L[α, O(1)] algorithm producing, for each prime q ≤ B, a certified relation row ρ_q
    (∏p^{ρ_{q,p}} ≡ 1) with |ρ_{q,q}| > Σ_{p≠q}|ρ_{q,p}| and L^{O(1)}-bit entries.
- *Proof.* (a)⇒(b): factor and output φ(N)e_q. (b)⇒(a): DD1 + Miller (Theorem DZ). ∎
- DSR ⇒ (b). PIN on the primes < 2 ln² N is weaker still (codex). So DSR is a more structured sufficient condition;
  its structure (smoothness) is what gives GRH existence (GX).

**Lemma PT (prime-modulus transfer) [proved; folklore-level padding argument].**
- *Statement.* Fix α ∈ (0, 1) and b, c, C. Suppose DSR_α(b, c) holds, even only for admissible N whose prime factors
  all exceed L_N[α, C]. Then for every large prime p, after a one-time preprocessing of cost L_p[α, c₀] with
  c₀ = O(max(b, C)) (finding r), there is a deterministic L_p[α, c + o(1)]-per-target algorithm with the following
  output. For every prime q ≤ B_N (B_N = L_p[α, b + o(1)]) it gives B_N-smooth u, v, w with |u|, |v| ≤ p^{3+o(1)},
  Ω(w) ≤ (1 + o(1))log₂ p and u ≡ q^{H′}wv (mod p), where H′ = 8⌈log₂(pr)⌉. All targets together cost
  L_p[α, b + c + o(1)].
- *Proof.*
  - Let r be the least prime > max(L_{p²}[α, max(b, C)], 3 ln²(p²)), found by sieving. Then r = p^{o(1)} < p.
  - N := pr is admissible: two distinct primes, both > max(B_N, 3 ln² N) and > L_N[α, C], since N < p².
  - Run DSR on N and reduce mod p. ln N = (1 + o(1))ln p. ∎
- *Meaning.*
  - DSR must solve the mod-p smooth-representation problem on inputs whose factorisation is supplied. Knowing the
    factorisation does not by itself trivialise DSR. This is not a separation: factoring-based methods could still
    contain useful ideas (codex).
  - The mod-p task is: produce *short prescribed-target* smooth witnesses q^{H′}wv ≡ u. Arbitrary full-rank relations
    mod p are trivial ((p − 1)e_q), and an unsaturated lattice does not give discrete logarithms.
  - We and the literature critic found no deterministic subexponential algorithm for prescribed-target short smooth
    witnesses mod p, or for a smooth number in a prescribed class, under GRH or otherwise.
  - Theorem DZ's recursion meets inputs like N = pr: preprocessing removes only factors ≤ max(B, 3 ln² N). Restricting
    DSR to N with all prime factors ≥ N^ε escapes PT, but DZ would then need deterministic removal of factors in
    (L^C, N^ε), which costs N^{ε/2+o(1)} by Pollard–Strassen (v1 said N^{ε/4}; corrected). That is exponential.

## 6. Literature (approach 4)

Verified by the literature critic unless marked (u) = unverifiable / not re-checked. Full list of fixes in §10.4.

- **Expanders under GRH.**
  - Jao–Miller–Venkatesan, JNT 129 (2009) 1491–1504: Thm 1.1 |λ| = O(n x^{1/2}log(xq)); Cor 1.3 for x ≥ (log q)^{B₀},
    B₀ > 2; the (Z/q)^× example with p^{±1}, p ≤ (log q)^{2+δ}.
  - Explicit GRH estimates: Ernvall-Hytönen–Palojärvi (Math. Comp. 2022; arXiv:2003.01925) (u: constants);
    Bach–Sorenson, Math. Comp. 65 (1996); Lamzouri–Li–Soundararajan 2015 + corrigendum 2017.
- **Derandomisation tools** (used only at the level of stated seed bounds, §2.3):
  - AKS STOC 1987; Impagliazzo–Zuckerman FOCS 1989; Nisan Combinatorica 12 (1992); INW STOC 1994; Nisan–Zuckerman
    JCSS 52 (1996);
  - Cohen–Doron–Goldgraber ECCC TR26-123 (2026); Hoza–Pyne–Vadhan ITCS 2021 (one accepting vertex); BHPP CCC 2022;
  - Gopalan–Kane–Meka FOCS 2015; Koucký–Nimbhorkar–Pudlák STOC 2011; De CCC 2011 (seed O((w⁸ + log 1/ε)log n), per
    Steinke ECCC TR12-083); Braverman–Rao–Raz–Yehudayoff FOCS 2010; Reingold–Steinke–Vadhan RANDOM 2013;
  - Rozenman–Vadhan RANDOM 2005; Reingold JACM 55 (2008); Reingold–Trevisan–Vadhan STOC 2006; Azar–Motwani–Naor
    Combinatorica 18 (1998); Impagliazzo–Wigderson STOC 1997; Bach JCSS 42 (1991) (random bits for GRH character
    hitting).
- **GRH derandomisation successes.** All exploit cosets or Galois structure (Fourier-small sets), exactly what FN
  excludes for the canonical smooth sets.
  - Ankeny 1952 and Bach 1990.
  - Shoup 1992, Math. Comp. 58 369–380: a small candidate set containing a primitive root; identifying one still needs
    the factorisation of p − 1.
  - Adleman–Lenstra STOC 1986: deterministic irreducible polynomials under GRH; unconditional for approximate degree.
  - Huang 1991; Rónyai 1992; Evdokimov 1994 (quasi-polynomial in the degree); Ivanyos–Karpinski–Saxena ISSAC 2009.
- **Closest analogue: deterministic prime finding.** Polymath4 (arXiv:1009.3956, "D.H.J. Polymath") /
  Tao–Croot–Helfgott, Math. Comp. 81 (2012) 1233–1246.
  - Randomized: log^{O(1)}N. Deterministic: N^{1/2+o(1)} (Odlyzko; also the RH bound).
  - The paper's partial result decides the parity of the number of primes in intervals of length ≤ N^{1/2+c} in
    N^{1/2−c+o(1)}. It is not a prime-finding algorithm below the square-root barrier.
  - Pseudodeterministic constructions (Oliveira–Santhanam STOC 2017; Chen–Lu–Oliveira–Ren–Santhanam FOCS 2023) are
    randomized and input-length-indexed, so they do not apply to the worst-case, input-dependent DSR.
- **Path-finding in expanders.** Charles–Goren–Lauter, J. Cryptology 22 (2009) 93–113; Petit–Lauter–Quisquater
  (SCN 2008); Sardari (LPS navigation). This is an analogy only: DSR is deterministic short-path finding in the GRH
  expander Cay(G, primes ≤ B).
- **Smooth numbers and characters.**
  - Hardy–Xu, "Character sums over smooth numbers" (arXiv:2607.00592, 2026) bounds the same first moment
    φ^{−1}Σ_χ|Σ_{n≤x, y-smooth}χ(n)| = o(√Ψ) for q ≥ x^{1+ε}. That is the complementary regime to FN (x ≈ q).
  - Hardy–Xu, "Helson's conjecture for smooth numbers" (arXiv:2511.03430); Harper (arXiv:2301.04390; 2607.01184,
    low moments with multiplicative twists; applicability to smooth indicators unverified); Konyagin–Shkredov
    (arXiv:1403.8129, Wiener norms in Z_p).
  - Soundararajan 2008 and Harper 2012 (moduli ≤ y^{4√e−ε}); Pascadi 2025 (level 5/8, averaged).
  - Bourgain–Glibichuk–Konyagin (subgroups ≥ p^ε; additive sums).
  - Harman 1999; Booker–Pomerance (Proc. AMS 2017); Munsch–Shparlinski (arXiv:1710.04705): existence of smooth
    representatives.
  - Bach–Sorenson (arXiv:2006.07445): random factored smooth integers.
  - Boneh JCSS 64 (2002); Zrałek, Math. Comp. 79 (2010) (deterministic p − 1).
  - None gives pointwise hits for L^{O(1)} explicit lists.
- **Factoring.** Lenstra–Pomerance JAMS 5 (1992) (rigorous randomized L[1/2]); Pomerance 1987 (rigorous randomized
  DL); Harvey, Math. Comp. 90 (2021) and Harvey–Hittmeir 2022 (deterministic N^{1/5}); Pollard–Strassen (factors ≤ M in
  M^{1/2+o(1)}).
- **Sieve inputs for BOX-∞.** Fundamental lemma (Friedlander–Iwaniec, *Opera de Cribro*, AMS Colloq. 57, Ch. 6;
  Halberstam–Richert; numbering (u)); Montgomery–Vaughan, Mathematika 20 (1973); Tenenbaum III.5 (Ψ ≪ x e^{−u/2});
  Mertens.

## 7. Experiments (pre-registered; toy sizes; calibration only)

Pre-registration `notes/dsr_prereg.md` (`409641a`). Summary `experiments/analysis/dsr_summary.json`.

**Disclosed deviations.**
- (1) The rule-1 exception (float64 FFT, E1) was declared in the prereg. Parseval holds to 6.7e−16.
- (2) A shared scratchpad log was overwritten by another agent. The JSONL files are authoritative.
- (3) D_slide turned out to have only |P|·t distinct multipliers (Sturmian windows), which I did not anticipate at
  registration. P2c/P2e therefore test that design defect, not explicit walks in general.

**E1 (Fourier-algebra norm; 18 moduli: primes 2^16–2^22, semiprimes 2^20, 2^22; factor bases of 5–39 primes).**

| prediction | result | verdict |
|---|---|---|
| P1a: r_A ∈ [0.7, 1.3] | S_can 0.549–0.632; S_bal 0.315–0.353 | **FAIL** (below; unexplained, §3) |
| P1b: slope ln A′ vs ln|S| ∈ [0.4, 0.6] | 0.474 (random 0.510) | pass |
| P1c: r_∞ > 1 in ≥ 75% | 18/18 (5.2–12.9; S_bal 13.9–72) | pass (registered rationale wrong: large-order characters) |

**E2 (worst-case coverage of all of G; 30 semiprimes at 20/22/24 bits; S_can at b = 0.8, δ = 0.06–0.10; S_bal at
b = 0.5, δ = 0.02–0.07).**

| prediction | result | verdict |
|---|---|---|
| P2a: ν(R_G) ∈ [0.8, 1.5] in ≥ 90% | 60/60, ν ∈ [0.80, 1.28] | pass |
| P2b: median m*(R_walk)/m*(R_G) ∈ [0.7, 1.5] | 1.02 | pass |
| P2c: D_box, D_slide cover within cap | > cap: D_box/S_can 10/30, D_slide/S_can 30/30, D_box/S_bal 29/30 | **FAIL** |
| P2d: median m*(D_box)/m*(R_G) ≥ 1.5 | 37.5 (censored values counted as > cap) | pass |
| P2e: median m*(D_slide)/m*(R_G) ∈ [0.7, 1.5] | > cap in all | **FAIL** (deviation 3) |
| P2f: median ν(D_box) at 24 bits ≤ 2× at 20 bits | 23.1 → 40.3 (22) → > cap (24) | **FAIL** (also on the available lower bounds) |

- **Diagnosis.**
  - Proved lower bound: m*(D_box) ≥ Ψ(N/ℓ, B) + 1 = N^{1−o(1)} via the target a = ℓ (Prop BOX).
  - Finite-data observation: ℓ is the worst target in all 20 uncensored instances (codex's recomputation).
  - The P2c/P2f failures are this exponential mechanism.
- **DSR targets only (descriptive).** Median T* = max_q first hit (20/22/24 bits):

  | list | T* | mean first hit |
  |---|---|---|
  | R_G | 34/50/59 | 10.2/12.8/15.4 |
  | R_walk | 41/54/71 | — |
  | D_orbit | 44/58/93 | — |
  | D_slide | 40/69/85 | — |
  | D_box (S_can) | 68/116/167 | — |
  | D_box (S_bal) | 1 663 / 3 304 / 9 894 | 553 / 1 229 / 3 180 |

  - The S_bal row also changes b and δ. It is consistent with zroute E4's inherited-vector mechanism, but without
    matched-parameter controls we do not attribute the whole factor to it.

## 8. What remains

A restricted sufficient route, not a characterisation of DSR. EQ-dom's dominant rows, or PIN, are weaker targets for
factoring.

> **(R)** A deterministic L^{c}-time map N ↦ W_N of legal multiplier data, together with a proof that for every
> admissible N ≥ N₀ and every prime q ≤ B some entry yields a legal DSR witness for q^H. Legal means: the
> combined multiplier, after reducing numerator and denominator, clears to integers u, v, w with |u|, |v| ≤ N³ and
> Ω(w) ≤ log₂ N.

Constraints on (R) that this note proves:
- **Format.** With the canonical test, W_N cannot consist only of positive integers < N^{1−ε}. With the Gauss test it
  cannot consist only of B-smooth rationals of reduced height < N^{1/2−ε} (BOX-∞). Lifted, signed or other tests are
  not constrained by BOX.
- **Proof method.** It cannot follow from walk spectra plus density, or from density plus near-optimal Fourier bias
  (SB). Through a Fourier-small minorant it works for neither S_can nor S_Gauss (FN). Second-moment templates certify
  all targets only with |W| ≳ φ/δ (AA, §2.7).
- **Scope.** By PT, it must also work modulo every prime with known small padding.

What could still work: an arithmetic input about how q^H interacts with smooth numbers modulo N that random sets lack;
a Fourier-small minorant of some other testable set; or a conditional-expectation estimator. We know no instance of
this pattern; deterministic prime finding is the analogous open problem.

## 9. Candidate cards

```text
ID: dsr-1 derandomised JMV walk (approach 1).
Mechanism: enumerate seeds of an explicit expander-walk / PRG / HSG; test canonical residue or Gauss split.
Cost: plain walk N^{1+o(1)} seeds (hitting) / N^{2+o(1)} (pointwise).
Failure point: black-box guarantees covering all density-delta accepting sets need > (1-delta)phi endpoints (SB(a));
  character-wise generators need eps <= N^{-1/2+o(1)} (Theorem FN); single-accept-state generators don't apply.
Status: black-box instantiations refuted (proved); arithmetic-specific HSG OPEN; IW-conditional (Prop IWD).

ID: dsr-2 explicit bounded-height multiplier boxes / k-fold prime products (codex blind round's missing lemma).
Status: REFUTED uniformly in N for canonical/Gauss tests (Prop BOX-infinity, unconditional; q = 2);
  almost-all-targets version holds (Prop AA, GRH).

ID: dsr-3 Bach/coset certificates.  Need a Fourier-small minorant.
Status: REFUTED for S_can (every N, unconditional) and S_Gauss (GRH) and their dense minorants (Theorem FN);
  other test sets OPEN.  E1: toy-scale norms below random (P1a failed; unexplained).

ID: dsr-4 smooth compression (u <= N^3).  Prop SC: height-coupled inverse-multiplier hitting; no constructive gain seen.

Structural: Prop EQ-dom (DD1 rows == deterministic factoring, ERH); Lemma PT (DSR => mod-p prescribed-target witnesses).
```

## 10. Reviews

### 10.1 GPT-6.1-sol blind round (`codex/dsr_1`, xhigh)
- Independently obtained Prop RW (with an explicit loose constant), the tool-by-tool failures, SB(a), the Fourier
  minorant and conditional-expectation sufficient conditions, Prop IWD, the Boneh threshold, the finite-N caveat, and
  "DSR not known to follow from factoring".
- Its proposed missing lemma is refuted by BOX-∞ (q = 2). Its probability: ≈ 10%.

### 10.2 GPT-6.1-sol adversarial round (`codex/dsr_2`, xhigh; `notes/dsr_codex2_report.md`)
- **Verdicts.**
  - SB: correct with fixes. Blanket derandomiser claim: overclaimed.
  - RW: correct with fixes. AA: correct, but "sharp" and §2.7 overclaimed.
  - BOX: correct. BOX-∞: correct with clarifications. Refutation of dsr_1 §6: correct.
  - FN: correct with fixes (ℓ | N counterexample; δ_Gauss proof; minorant inconsistency).
  - PT: correct with cost fixes. Pollard–Strassen exponent wrong. Separation rhetoric overclaimed.
  - EQ-dom: correct; "minimal" inaccurate.
  - **SC: WRONG as stated** (counterexample, height coupling).
  - IWD: correct, but the non-implication claim was overclaimed.
  - TH-wrap: overclaimed. E1 diagnosis: overclaimed. Literature: incomplete.
- **Its FATAL items, all fixed in v2:**
  - SC replaced by the height-coupled statement;
  - the blanket derandomiser impossibility withdrawn and rescoped (§2.2);
  - TH-wrap withdrawn as "exact weakest" and replaced by the restricted route (R) (§8).
- **MAJOR and MINOR items, all applied:** FN's ℓ ∤ N and δ_Gauss proof; SB constants and quantifiers; PT costs and
  Pollard–Strassen; AA wording and §2.7; first-moment prior art; the E1 diagnosis downgraded; E2 censoring wording;
  RW ranges, sampling and heights; B_* vs B_N.

### 10.3 Opus proof critic (read-only) — ACCEPT WITH FIXES, no FATAL
- It re-derived SB, RW, AA, BOX, BOX-∞, FN, PT, EQ-dom, SC (⇒) and IWD.
- It numerically checked FN's inequality on p ∈ {100003, 524287, 999983} (30 cases; probability exactly 0 whenever
  ℓB^k < p; vacuous at that size).
- **MAJOR:**
  - derandomiser overclaim (fixed, §2.2);
  - TH-wrap necessity and the inverse-multiplier escape (fixed, §8 and BOX(iii)).
- **MINOR, all applied:**
  - ℓ ∤ N; δ_Gauss count; "as random" for S_can only; SB constant;
  - AA labelled GRH; §2.7 one-versus-all;
  - SC bound BN³ (superseded by codex's exact fix);
  - PT cost of r and per-target time; Pollard–Strassen;
  - q = 2 scope; toy-scale qualifiers; IWD computable δ.

### 10.4 Opus literature critic (read-only) — OK with fixes
- **Corrections:**
  - Oliveira–Santhanam is STOC 2017;
  - Pollard–Strassen gives N^{ε/2};
  - De's seed is poly(w);
  - Shoup qualified;
  - TCH length restriction and Polymath credit;
  - the LLS corrigendum;
  - JMV's Cor 1.3 condition;
  - IKS venue.
- **Key finding:** the S_can half of FN is *unconditional* by Pólya–Vinogradov duality. This is now the proof of
  FN(i), with the S_can part downgraded to folklore-level.
- **Prior art added:** Hardy–Xu (two papers), Harper, Konyagin–Shkredov, Booker–Pomerance, Munsch–Shparlinski,
  Zrałek, Lenstra–Pomerance, Bach 1991, Bach–Sorenson, RTV/RSV/Steinke/BRRY/AMN, Petit–Lauter–Quisquater.
- PT and BOX-∞ relabelled folklore and elementary respectively.
- **No deterministic subexponential algorithm** for prescribed-target relation collection mod p, or for a smooth
  number in a prescribed class, was found.

### 10.5 Confirmation round on v2 (Opus proof critic) — pending.

## 11. Reproduction

```
PYTHONPATH=. .venv/bin/python -m src.dsr.harness e1 experiments/raw/dsr_e1.jsonl
PYTHONPATH=. .venv/bin/python -m src.dsr.harness e2 experiments/raw/dsr_e2.jsonl
PYTHONPATH=. .venv/bin/python -m src.dsr.analyze > experiments/analysis/dsr_summary.json
PYTHONPATH=. .venv/bin/python -m src.dsr.posthoc_e1 > experiments/raw/dsr_e1_posthoc.jsonl      # post hoc
PYTHONPATH=. .venv/bin/python -m src.dsr.posthoc_adv > experiments/raw/dsr_posthoc_adv.jsonl    # post hoc
.venv/bin/python -m pytest -q tests/test_dsr.py
```
E1 peaks at ≈ 1.24 GB RSS (22-bit 2-D FFTs); E2 takes ≈ 15 min.
