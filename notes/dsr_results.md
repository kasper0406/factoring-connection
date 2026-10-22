# dsr: Hypothesis DSR (deterministic smooth representation): derandomised walks, certificates, compression

Agent `dsr`, phase 23, 2026-10-22. Version 1 (before reviews; §10 records them).
- Pre-registration: `notes/dsr_prereg.md` (commit `409641a`, before any statistic was computed).
- Code: `src/dsr/core.py` (operational: sees only N, deterministic), `src/dsr/harness.py` (E1 diagnostic with
  discrete logs and float FFTs; E2 seeded random controls), `src/dsr/analyze.py`, `src/dsr/posthoc_e1.py` and
  `src/dsr/posthoc_adv.py` (post hoc, labelled). Tests: `tests/test_dsr.py` (7 tests).
- Data: `experiments/raw/dsr_e1.jsonl`, `dsr_e2.jsonl`, `dsr_e1_posthoc.jsonl`, `dsr_posthoc_adv.jsonl`;
  summary `experiments/analysis/dsr_summary.json`.
- GPT-6.1-sol: blind round `codex/dsr_1.{md,out}` (§10.1); adversarial round `codex/dsr_2.{md,out}` (§10.2).
- Labels: **proved**; **proved (folklore-level)**; **conditional** (on a named hypothesis); **open**; **observation**
  (a limitation of a proof template, not a lower bound); **heuristic**; **experimentally supported**; **refuted**.
  Novelty: at most "no prior statement found by us".

## Verdict (plain language)

**DSR was not proved, and nothing here comes close. The derandomised-walk idea fails for a reason that can be made
exact, and the simplest explicit constructors are provably non-uniform.**

1. **Derandomising the JMV walk with generic tools cannot work [proved].** Under GRH the random walk on
   Cay((Z/N)^×, primes ≤ B) hits the smooth set after (1+o(1)) ln N/ln B steps. That costs ≈ log₂ N random bits:
   N^{1+o(1)} seeds, not L^{O(1)} (Prop RW).
   - Every general-purpose derandomiser (AKS/IZ amplification, Nisan/INW/Nisan–Zuckerman, the width-independent
     HPV/BHPP generators, GKM Fourier shapes, KNP/De group products, derandomised squaring, Reingold) proves its
     guarantee for a class of tests that contains *every* density-δ subset of the group. Any family hitting all of
     those has more than (1−δ)φ(N) endpoints (Theorem SB(a)).
   - Even "S has density δ and near-optimal Fourier uniformity" does not force a hit for lists of size ≤ √(φ ln φ/δ)
     (Theorem SB(b)).
   - So the GRH spectral bounds, which are statements about the *walk*, cannot imply the needed hitting property. A
     proof must use arithmetic of the specific smooth set S_N beyond density and Fourier bias.
2. **The natural explicit lists are provably non-uniform [proved, unconditional; new as far as we found].**
   Prop BOX-∞: for every b, k and every large h there is an admissible N < 2^h for which the DSR target 2^H mod N is
   a tiny non-smooth integer r⁸. For that N, no multiplier v ≤ N^{1−o(1)} makes 2^H·v mod N smooth (canonical
   residue), and no v ≤ N^{1/2−o(1)} gives a smooth Gauss split.
   - This **refutes** the most natural missing lemma, GPT-6.1-sol's k-fold prime-product correlation estimate
     (§2.6), as stated.
   - A uniform DSR constructor must use wrapping multipliers of height ≥ N^{1−o(1)}, or non-canonical
     representations.
   - Experiments (E2, post hoc) show the same mechanism at 20–24 bits. Adversarial instances occur naturally: 7 of
     4000 random 24-bit N.
3. **Almost-all versions are easy, and the exceptional set is genuinely there [proved; folklore-level].** Prop AA:
   under GRH, for every N, the explicit list of k-fold prime products works for all targets outside an exceptional set
   E_N with |E_N|/φ ≤ L^{−kb+1/(2b)+o(1)}. BOX shows E_N contains every small non-smooth residue, so the bound is
   tight up to L^{1/(2b)+o(1)}. DSR asks whether the π(B) specific targets q^H avoid E_N. They need not.
4. **DSR is a prime-modulus problem in disguise [proved, trivial].** Lemma PT: DSR (even restricted to N whose prime
   factors all exceed L^C) implies the same smooth-representation constructor modulo every prime p. Apply it to
   N = p·r with r the least prime above L^C; the algorithm then knows the factorisation and still has to solve the
   mod-p problem.
   - So DSR packs in the derandomisation of index-calculus relation collection in F_p^×, which we found open.
     Factoring N cannot help, and DSR is *not* known to follow from deterministic factoring (Prop EQ-dom: what
     Theorem DZ minimally needs is equivalent to deterministic factoring itself).
5. **Approaches 2 and 3 reduce to the same hitting problem [proved reductions + observations].**
   - Bach/Ankeny certificates need a Fourier-algebra-small minorant of S (codex). **Theorem FN [proved under GRH]:**
     for every admissible N, ‖1_{S_can}‖_A ≥ N^{1/2−o(1)} (sharp) and ‖1_{S_Gauss}‖_A ≥ N^{1/4−o(1)}. One bad
     translate (target ℓ, the least prime > B) against GRH-biased k-fold prime products forces it. This settles
     zroute's speculative conjecture: Bach-type certificates are impossible for the canonical smooth sets, and for
     their dense minorants.
   - E1 found toy-scale norms *below* random (r_A ≈ 0.55 / 0.32; pre-registered P1a failed), growing like |S|^{0.47}.
     This was diagnosed post hoc as an artefact: with 8–36 primes the factor base is not yet an expander (empirical
     λ = 0.57–0.99; GRH bound vacuous).
   - Smooth compression with u ≤ N³ is exactly translate-hitting with the N^{2−o(1)} inverse multipliers d^{−1}
     (Prop SC). The height freedom helps existence proofs (GX), not constructors.
6. **Conditional derandomisation [proved (folklore-level); conditional].** Under GRH plus Impagliazzo–Wigderson
   circuit hardness, DSR_{1/2} holds (Prop IWD), so ERH + IW ⇒ deterministic L[1/2] factoring (also polyimpl DR).
   Nothing shown here makes GRH supply such a hitting set.
7. **Closest analogue:** deterministic prime finding (Polymath4; Tao–Croot–Helfgott 2012). It has the same shape:
   a dense, efficiently testable set, a trivial randomized algorithm, and an open deterministic version even under
   RH (best deterministic N^{1/2+o(1)}).

**The exact remaining obstacle (§8).** For every admissible N and every prime q ≤ B, deterministically produce an
L^{O(1)}-size explicit list of *wrapping* multipliers m (products of smooth numbers and their inverses, of height
≥ N^{1−o(1)}) and prove that q^H·m mod N lands in a dense testable smooth set for one of them. Equivalently, prove that
the specific targets q^H avoid the exceptional set E_N of such a list. The same must hold modulo every prime p
(Lemma PT).

**Probability assessment (subjective).** DSR provable from GRH with known techniques: ≲ 3% (ours; lowered from
zroute's 5% by Lemma PT and BOX-∞); codex blind ≈ 10%.

---

## 1. Setting

- N is admissible: odd, not a prime power, no prime factor ≤ max(B, 3 ln² N). G = (Z/N)^×, φ = φ(N).
  L = L_N[1/2, 1], B = L^b (b is zroute's a, renamed to free the letter a), H = 8⌈log₂ N⌉. Targets a = q^H mod N,
  q ≤ B prime.
- S_can := {x ∈ G : least positive residue of x is B-smooth}. δ_can = Ψ(N − 1, B)/φ = L^{−1/(2b)+o(1)} (all B-smooth
  integers below N are units).
- S_Gauss: x whose Gauss-reduced vector (u, v) (u ≡ xv) has both halves B-smooth. S_bal (experiments): the subset with
  u² + v² < N, where the shortest vector is unique up to sign (§7).
- Fourier: Ŝ(χ) = Σ_{s∈S} χ(s); bias β(S) = max_{χ≠χ₀}|Ŝ(χ)|; algebra norm ‖1_S‖_A = φ^{−1}Σ_χ|Ŝ(χ)|.
- **Translate-hitting TH(S, T, W):** for every target a ∈ T some w ∈ W has a·w ∈ S; equivalently T ⊆ S·W^{−1}. Every
  DSR variant here is an instance (§2.6, §4).
- **Orbit/relaxed form of what Theorem DZ needs [proved, trivial].** DD1 only needs, for each q ≤ B, *some* row
  ρ_q = j·e_q + e(w) − e(u) + e(v) with j > Ω(u) + Ω(v) + Ω(w) and q^j·w·v ≡ u. So targets q^{H+j} (orbit), q^H·w
  (multipliers) and any j ≤ L^{O(1)} are all admissible. `tests/test_dsr.py::test_orbit_hit_gives_dominant_row`.

## 2. Approach 1: deterministic search and derandomised walks

### 2.1 Explicit spectral data (GRH)

**Prop RW (the randomized constructor, explicit) [proved, folklore-level; JMV 2009 + Fourier inversion].**
- Assume GRH for Dirichlet L-functions mod N. Let μ be uniform on the primes ≤ y (y ≤ B; all are units) and
  λ := max_{χ≠χ₀}|μ̂(χ)|.
- JMV Thm 1.1 (specialised to Q; verified on arXiv:0811.0647) bounds the nontrivial adjacency eigenvalues of the
  undirected graph (generators p^{±1}) by O(√y log(Ny)). The same character-sum argument bounds |Σ_{p≤y}χ(p)|.
  Hence λ ≤ C·log y·log(Ny)/√y. Codex gives an explicit loose constant via Ernvall-Hytönen–Palojärvi:
  λ ≤ Λ(N, y) := min(1, 10⁶√y·ln³(Ny)/π(y)). We have not audited that constant (§10).
- For any S ⊆ G of density δ and any start a: |Pr[a·p₁⋯p_t ∈ S] − δ| ≤ λ^t·(φδ(1−δ))^{1/2}. This is Parseval and
  Cauchy–Schwarz on (1/φ)Σ_{χ≠χ₀}χ̄(a)μ̂(χ)^t \overline{Ŝ(χ)}.
- **Lengths at y = B = L^b.** ln(1/λ) ≥ (b/2 − o(1))√(ln N ln ln N).
  - Hitting S with probability ≥ δ/2 needs t_hit = ⌈ln(4φ/δ)/(2 ln(1/λ))⌉ = (1/b + o(1))√(ln N/ln ln N).
  - Pointwise mixing (every class, GX) needs t_∞ = ⌈ln(2φ)/ln(1/λ)⌉, i.e. twice that.
- **Randomness and endpoints.** A walk uses t_hit·log₂π(B) = (1 + o(1))log₂ N bits (pointwise: 2 log₂ N). The seed
  space has N^{1+o(1)} elements (N^{2+o(1)}). Endpoints a·v with v ≤ B^{t} = N^{1+o(1)}.
- **The DSR output.** u = least residue of a·v < N, v = p₁⋯p_t ≤ N^{1+o(1)}, w = 1, all B-smooth.
  Ω(v) = t ≪ log₂ N. Expected cost L^{b+1/(2b)+o(1)} per target (codex, same). Under GRH, DSR is therefore a
  *randomized* L[1/2] search problem; the question is purely derandomisation.
- **The budget.** L^{c} seeds means seed length c·√(ln N ln ln N)/ln 2 = o(log₂ N). The deficit against the plain walk
  is a factor N^{1−o(1)}.

### 2.2 Theorem SB (spectral blindness)

**Theorem SB [proved; (a) is trivial and was found independently by codex].** Let G be a finite abelian group of order
φ ≥ 16, and W ⊆ G (or the support of a multiset).
- (a) If |W| ≤ (1 − δ)φ, some S ⊆ G ∖ W has |S| ≥ δφ. In particular no hitting guarantee whose hypotheses concern only
  (W, the generators, their spectrum) and the density of S can hold. This applies to adaptive query algorithms too:
  answer "not in S" to every query.
- (b) Suppose |W| ≤ φ/2, 0 < δ ≤ 1/4 and δφ ≥ 24 ln(8φ). Then there is S ⊆ G ∖ W with
  ||S| − δφ| ≤ (6δφ ln 8φ)^{1/2} and β(S) ≤ 2δ·β(1_W) + (12δφ ln 8φ)^{1/2}.
  If moreover |W| ≤ (3φ ln(8φ)/δ)^{1/2}, then β(S) ≤ 2(12δφ ln 8φ)^{1/2}. That is within a factor 7(ln 8φ)^{1/2} of
  the Parseval minimum (δφ(1 − δ))^{1/2}(1 − o(1)) that every set of that size satisfies.

*Proof of (b).*
- Put δ′ := δφ/(φ − |W|) ≤ 2δ ≤ 1/2. Include each g ∈ G ∖ W independently with probability δ′.
- For χ ≠ χ₀, E Ŝ(χ) = δ′Σ_{g∉W}χ(g) = −δ′·1̂_W(χ), since Σ_G χ = 0. So |E Ŝ(χ)| ≤ 2δ·β(1_W).
- X := Re(Ŝ(χ) − EŜ(χ)) = Σ_{g∉W}(ξ_g − δ′)Re χ(g) is a sum of independent mean-zero terms of absolute value ≤ 1,
  with variance ≤ δ′(φ − |W|) = δφ =: σ². The same holds for Im and for |S| − δφ.
- Bernstein: Pr[|X| ≥ s] ≤ 2exp(−s²/(2σ² + 2s/3)). Take s := (6σ² ln 8φ)^{1/2}. Since σ² ≥ 24 ln 8φ, s ≤ σ²/2, so the
  exponent is ≥ 6σ² ln(8φ)/((7/3)σ²) ≥ 2 ln 8φ, and each probability is ≤ 2/(64φ²).
- There are 2(φ − 1) + 1 events, so the union bound is < 1. Some outcome satisfies all of them, and then
  |Ŝ(χ) − EŜ(χ)| ≤ √2·s.
- The last claim uses β(1_W) ≤ |W|. ∎

*Consequences.*
- **Generic derandomisers.** Each tool in §2.3 states its guarantee for a class of tests containing, for the walk's
  natural branching program (state x ∈ G, step x ↦ xp), every accepting set S′ ⊆ G of density δ: width-φ permutation
  programs with arbitrary accepting sets.
  - By (a), any generator correct for that class has ≥ (1 − δ)φ distinct endpoints, so seed length ≥ log₂((1−δ)φ).
  - This is a lower bound for those *guarantees*, not for algorithms that use S_N's arithmetic.
- **Fourier certificates.** By (b), no theorem of the form "β(S) ≤ X and |S| ≥ δφ ⇒ W ∩ S ≠ ∅" with
  X ≥ 2(12δφ ln 8φ)^{1/2} can hold for lists of size ≤ (3φ ln(8φ)/δ)^{1/2}. This includes every L^{O(1)} list.
  Combined with zroute's Prop FC, sup-norm, Parseval and Hölder certificates are all excluded. What remains is the
  algebra norm, i.e. structure (§3).
- **GRH as used here.** GRH enters through character sums over primes, i.e. through W's spectrum. Bounds on twisted
  smooth sums Ŝ_N(χ) are sup-norm information about S and fall under (b). Any successful proof must use a property of
  S_N that some (density, near-optimal bias) set lacks. This is not a claim that GRH has no further consequences; it
  locates what they would have to say.

### 2.3 The tools, one by one

The (1 − δ)φ bound of SB(a) explains every row's failure point. Codex `dsr_1` §1(b) supplied the HPV/BHPP, De, RV
and NZ entries; we checked the arithmetic, not every cited theorem (literature critic, §10).

| tool | what it saves | cost for DSR / failure point |
|---|---|---|
| AKS 1987; Impagliazzo–Zuckerman 1989 (expander-walk amplification) | randomness for *repeating* an experiment | still needs one uniform start, about log₂ φ bits: N^{1+o(1)} seeds |
| Nisan 1992; INW 1994 (width W, length ℓ ROBPs) | seed O(log ℓ·log(Wℓ/ε)) | W = φ, ℓ = O(log N): seed O(log N log log N), worse than the plain walk |
| Nisan–Zuckerman 1996 | poly random bits in space S to O(S) bits | S ≈ log N: N^{O(1)} seeds |
| Hoza–Pyne–Vadhan 2021; Bogdanov–Hoza–Prakriya–Pyne 2022 (width-independent, regular/permutation programs) | log W removed | only for *one* accepting state. δφ accepting states need ε ≈ 1/φ, which restores log φ |
| Gopalan–Kane–Meka 2015 (Fourier shapes) | fools each endpoint character | 1_S needs ε·‖1_S‖_A < δ. Generic ‖1_S‖_A ≤ √(δφ) forces log(1/ε) ≈ log N. With a hypothetical ‖1_S‖_A = L^{O(1)}, plain k-step walks already suffice (§3) |
| KNP 2011; De 2011 (group products) | small-bias for fixed groups | constants or terms in log|G| = Θ(log N) |
| Rozenman–Vadhan 2005 (derandomised squaring) | expansion at low degree | error must reach ≈ φ^{−1/2}; any distribution with all biases ≤ β has support ≥ φ/(1 + φβ²) (Parseval) |
| Reingold 2008; universal exploration sequences | deterministic traversal | length ≥ φ − 1, products far above N³ |
| Impagliazzo–Wigderson 1997 (circuit hardness) | seed O(log T) for size-T circuits | **works**, conditionally (Prop IWD) |

**Prop IWD [conditional; folklore-level; codex `dsr_1` §4 agrees].**
- *Statement.* If some f ∈ E has circuit complexity ≥ 2^{εm} at all large lengths m, then for every b > 0 under GRH,
  DSR_{1/2}(b, c) holds for some c = c(b, ε).
- *Proof.* Repeat the Prop RW experiment ⌈8/δ⌉ times with independent coins. Success probability ≥ 0.98, and the
  circuit size is T = L^{O(1)}. The IW generator with seed O(log T) fools size-T circuits with error < 1/2, so some
  seed succeeds. Enumerate all T^{O(1)} seeds and verify each output. Primes are sampled from coins uniformly on the
  first 2^{⌊log₂π(B)⌋} primes, to which the same GRH bound applies.
- With Theorem DZ: ERH + IW ⇒ deterministic L[1/2, O(1)] factoring. This is consistent with polyimpl Theorem DR, which
  derandomises Lenstra–Pomerance directly.
- **No complexity barrier is known in the other direction.** DSR is one explicit problem, so a deterministic DSR
  algorithm would be a hitting set for one test family only. It would not imply circuit lower bounds (compare AKS
  primality).

### 2.4 What *would* suffice (codex, checked)

- (i) A minorant f ≤ 1_S with E f ≥ δ/2 and ‖f‖_A ≤ L^{O(1)}. Then a fixed number k of plain prime steps suffices, by
  enumeration of k-fold products: cost L^{kb+O(1)}. But for S_can and S_Gauss no such minorant exists for any N under
  GRH. The proof of Theorem FN applies verbatim to f ≤ 1_S with E f ≥ δ/2, because the translate ℓ·W_k misses S.
  Any useful f must live on a different test set.
- (ii) An efficient estimator of h_j(x) = Pr[x·p₁⋯p_j ∈ S] to additive error δ/(8t): the method of conditional
  expectations. GRH supplies it only for large j, when the remaining walk mixes.
- (iii) IW-type hardness (Prop IWD).

### 2.5 Prop AA (almost all targets, every N)

**Prop AA [proved under GRH; folklore-level: Parseval / large-sieve second moment].**
- *Statement.* Let W_k be the multiset of ordered k-tuples of primes ≤ y, of size m = π(y)^k, with products ≤ y^k, and
  let S ⊆ G have density δ. The set E = {a ∈ G : a·w ∉ S for all w ∈ W_k} satisfies |E| ≤ λ^{2k}φ(1 − δ)/δ.
  At y = B = L^b and fixed k: |E|/φ ≤ L^{−kb+1/(2b)+o(1)} for S = S_can. The cost is one L^{kb+o(1)} enumeration per
  target.
- *Proof.*
  - h(a) := #{w ∈ W_k : aw ∈ S} = φ^{−1}Σ_χ χ̄(a)Ŝ(χ)\overline{Ŵ(χ)}, with Ŵ(χ) = (Σ_{p≤y}χ(p))^k. The χ₀ term
    is δm.
  - Parseval over a: Σ_a|h(a) − δm|² = φ^{−1}Σ_{χ≠χ₀}|Ŝ|²|Ŵ|² ≤ (λ^k m)²·|S|(1 − δ).
  - Each a ∈ E contributes (δm)². ∎
- *Contrast.* A random list of size m leaves φ(1 − δ)^m uncovered in expectation: exponentially small in m, not
  polynomially.
- *Sharpness.* By BOX (§2.6), E ⊇ {a : least residue r ≤ N/y^k, r not B-smooth}, which has size ≥ φ·y^{−k}(1 − o(1)).
  So AA is tight up to the factor λ^{2k}y^k/δ = L^{1/(2b)+o(1)}.
- *What DSR needs.* The π(B) = L^{b+o(1)} specific targets q^H (a φ^{−1+o(1)} fraction of G) must avoid E_N. AA says nothing about any particular
  set of φ^{−1+o(1)} density. Second moments taken over the targets themselves (over q, or over an orbit j) need
  |W|·#targets ≳ φ/δ by the same Parseval loss (computations in §2.7).

### 2.6 Prop BOX and BOX-∞ (short boxes of multipliers are not uniform)

**Prop BOX [proved, trivial].**
- (i) Let the least positive residue r of a be not B-smooth. Then for every integer w with 1 ≤ w < N/r, a·w mod N = rw
  is not B-smooth, so a·w ∉ S_can.
- (ii) If moreover r²w² + 1 < N, the lattice {(u, v) : u ≡ awv} has (rw, 1) as its unique (up to sign) shortest
  vector. Two independent vectors of squared length < N would have |det| < N. So the Gauss split is (rw, 1) and
  a·w ∉ S_Gauss.

**Prop BOX-∞ (adversarial N exist for every size) [proved; unconditional; no prior statement found by us].**
- *Statement.* Fix b > 0. For all sufficiently large h there is an admissible N with 2^{h−1} < N < 2^h (so
  H = 8⌈log₂N⌉ = 8h) such that 2^H ≡ r⁸ (mod N) with B < r < R := B^{u_h}, r not B-smooth. Here
  u_h := ⌈3 ln ln B⌉ and B := ⌊L_{2^h}[1/2, b]⌋, so R = L^{O(b ln ln N)} = N^{o(1)}. Consequently:
  - for every v ≤ N/r⁸ = N^{1−o(1)}, 2^H·v mod N = r⁸v is not B-smooth;
  - for every v ≤ (N/2)^{1/2}/r⁸ = N^{1/2−o(1)}, the Gauss split of 2^H·v is (r⁸v, 1), not smooth.
- *Proof.*
  - Let z := max(B, 3h² ln² 2) + 1, which bounds max(B(N), 3 ln² N) for N < 2^h. Sift the interval
    A := (2^h − R, 2^h − B) ∩ Z, of length R(1 − o(1)), by the primes < z.
  - **Fundamental lemma of the sieve** (e.g. Friedlander–Iwaniec, *Opera de Cribro*, Cor. 6.10; Halberstam–Richert
    Thm 2.5), at level D = R^{1/2}: |r_d| ≤ 1 for intervals, and s = ln D/ln z ≥ u_h/2 − o(1) → ∞. So
    S(A, z) = |A|·V(z)(1 + O(e^{−s})) + O(R^{1/2}). With V(z) = ∏_{p<z}(1 − 1/p) ≍ 1/ln B (Mertens),
    S(A, z) ≥ (e^{−γ} − o(1))R/ln B.
  - **Remove primes.** By Brun–Titchmarsh (Montgomery–Vaughan) there are ≤ 2R/ln R = 2R/(u_h ln B) = o(R/ln B) of them.
  - **Remove prime powers p^j (j ≥ 2).** For each j at most one lies in the interval (consecutive j-th powers near
    2^h differ by ≥ 2^{(h−1)/2} > R), so at most h in total.
  - **Remove N with r := 2^h − N B-smooth.** These number ≤ Ψ(R, B) ≪ R·e^{−u_h/2} ≤ R/(ln B)^{3/2} = o(R/ln B)
    (Rankin's bound Ψ(x, y) ≪ x e^{−u/2}, u = ln x/ln y; Tenenbaum III.5).
  - Some N ∈ A remains. It is z-rough (hence odd, with no prime factor ≤ max(B, 3 ln² N)), composite and not a prime
    power, so admissible. And r = 2^h − N ∈ (B, R) is not B-smooth.
  - Then 2^h ≡ r, so 2^H = 2^{8h} ≡ r⁸ (mod N), with r⁸ ≤ R⁸ = N^{o(1)}. Prop BOX applies with a = 2^H and least
    residue r⁸. Its prime factor > B divides r. ∎
- *Consequences.*
  - (1) **Refutation [proved].** Codex `dsr_1` §6's missing lemma asks for Σ_{p₁..p_k≤B} 1_{S_can}(q^H p₁⋯p_k)
    ≥ η·π(B)^k·δ for every admissible N and every q ≤ B. It is false for q = 2 on the N of BOX-∞, where the sum is 0.
    So is any uniform lemma pairing a multiplier box of height ≤ N^{1−ε} with the canonical test, or height
    ≤ N^{1/2−ε} with the Gauss test, and every Fourier-minorant route (§2.4(i)) for S_can or S_Gauss with
    bounded-height multipliers.
  - (2) Theorem DZ needs a row for q = 2. DSR's freedom (|u|, |v| ≤ N³; w with Ω ≤ log₂N, hence height up to
    B^{log₂N}) is not violated: wrapping multipliers, or lifts u = r⁸v + jN, escape BOX. A uniform constructor
    **must** use them. That is exactly where explicit control is lost: the walk endpoints of Prop RW have v ≈ N^{1+o(1)}.
  - (3) Orbit lists q^{H+j} escape after O(log N) steps (r⁸2^j wraps once 2^j > N/r⁸). BOX-∞ does not refute them. We
    have no adversarial construction for wrapping explicit lists.
- *Experiments (post hoc, labelled; `src/dsr/posthoc_adv.py`).*
  - (i) In all 30 E2 instances, D_box's worst-case coverage satisfies m* ≥ Ψ(N/ℓ, B) + 1, where ℓ is the least prime
    > B and the target a = ℓ is the witness. The gap is 0–3 in all 20 uncensored instances; in the 10 censored ones
    Ψ + 1 exceeds the cap.
  - (ii) Scanning 4000 random admissible 24-bit N found 7 (N, q) pairs with q^H mod N = r ≤ N/2¹², r non-smooth. On
    those targets D_box first hits at index Ψ(N/r, B) + 1..7 (2 685–8 038 tries), against 3–76 for the random walk.

### 2.7 Second moments over the targets (template computations, not lower bounds)

- **Over q.** Σ_{q≤B}|h(q^H) − δm|² splits into a diagonal (χ^H = χ′^H), bounded by π(B)·|Ĝ[H]|·(λ^km)²·δ, and an
  off-diagonal bounded via GRH by C√B log²(NB)·(λ^km)²·‖1_S‖′_A². With only ‖1_S‖_A² ≤ |S| = δφ, making the number of
  bad q < 1 needs λ^{2k} ≲ δ/(√B φ): |W| ≥ φ^{1−o(1)}.
- **Over an orbit q^{H+j}, j < J.** The pointwise bound needs λ^k‖1_S‖_A < δ (FC(c) again). Parseval-type averaging
  needs J·m ≳ φ/δ.
- **Exceptional-set version.** T = {q^H} can lie inside E only if |E| ≳ B/polylog. That again needs |W| ≳ φ/(Bδ).
- All three lose exactly where zroute's Prop FC does: only L² information on S is available.

## 3. Approach 2: Bach-style certificates and coset structure

- **Why Bach works [proved, standard].** A coset has ‖1_C‖_A = 1 and a subgroup complement 2(1 − 1/k), and GRH needs
  only the few characters trivial on the subgroup. Codex adds that arbitrary cosets of index k need primes up to
  ≈ k² log² N (Lamzouri–Li–Soundararajan Thm 1.4), so even coset hitting is cheap only at bounded index.
- **Decomposition into cosets [no positive result; one example].**
  - The set of classes with *some* smooth representative ≤ N³ is all of G under GRH (GX): one coset, with no
    representative supplied.
  - The efficiently testable sets (S_can, S_Gauss) are truncations, and truncation destroys semigroup closure. Codex's
    example: N = 1009·1013, B = 2, S_can = {2^j : j ≤ 19} contains no nontrivial coset.
  - A useful decomposition would be a minorant with ‖f‖_A = L^{O(1)} (§2.4(i)). Theorem FN excludes even
    ‖1_S‖_A ≤ N^{1/2−ε} for S_can (N^{1/4−ε} for S_Gauss), for every N under GRH. For minorants, and for other test
    sets, it is open.
**Theorem FN (smooth sets have near-maximal Fourier-algebra norm) [proved under GRH; no prior statement found by
us; settles zroute's speculative FC(c) conjecture for these sets].**
- *Statement.* Assume GRH for Dirichlet L-functions mod N, b > 0 fixed. For every admissible N ≥ N₀(b):
  ‖1_{S_can}‖′_A := φ^{−1}Σ_{χ≠χ₀}|Ŝ_can(χ)| ≥ N^{1/2−o(1)}, and ‖1_{S_Gauss}‖′_A ≥ N^{1/4−o(1)}.
  Since ‖1_S‖_A ≤ |S|^{1/2} ≤ N^{1/2} always, the first bound is sharp up to N^{o(1)}.
- *Proof.*
  - Let ℓ be the least prime > B (so ℓ ≤ 2B), μ uniform on the primes ≤ B, and λ := max_{χ≠χ₀}|μ̂(χ)|. Under GRH
    λ ≤ C log B·log(NB)/√B = B^{−1/2+o(1)} (Prop RW).
  - Expand E[1_S(ℓ·p₁⋯p_k)] = φ^{−1}Σ_χ Ŝ(χ)χ̄(ℓ)\overline{μ̂(χ)}^k. The χ₀ term is δ, and the rest has modulus
    ≤ λ^k‖1_S‖′_A. Hence ‖1_S‖′_A ≥ (δ − Pr[ℓp₁⋯p_k ∈ S])·λ^{−k}.
  - *S_can.* Take k := ⌊ln(N/(2B))/ln B⌋. Then ℓ·p₁⋯p_k ≤ 2B^{k+1} < N, so ℓp₁⋯p_k is its own least residue and is
    divisible by ℓ > B: the probability is 0 (Prop BOX).
    - λ^{−k} = B^{k/2}·exp(−O(k ln ln N)) = N^{1/2−o(1)}, because k = O(√(ln N/ln ln N)).
    - δ_can = L^{−1/(2b)+o(1)} = N^{−o(1)}.
  - *S_Gauss.* Take k := ⌊ln(√N/(4B))/ln B⌋. Then s := ℓp₁⋯p_k satisfies s² + 1 < N, so (s, 1) is the unique shortest
    vector (Prop BOX(ii)) and the split is not smooth. λ^{−k} = N^{1/4−o(1)}. δ_Gauss ≥ N^{−o(1)}: the coprime B-smooth pairs (u, v),
    v > 0, u² + v² < N give pairwise distinct x = u/v (Prop BOX(ii) uniqueness), and there are
    ≫ N·(Ψ(√N/2, B)/√N)² = N·L^{−O(1/b)} of them. ∎
- *Meaning.*
  - Bach/Ankeny-type certificates, which need ‖1_S‖_A = L^{O(1)} (zroute FC(c)), are impossible for S_can and S_Gauss
    for **every** N under GRH. In L¹ Fourier norm the canonical smooth sets are as spread out as random sets, up to
    N^{o(1)}.
  - A single bad translate (the target ℓ) for a small-bias family forces a large norm. Conversely, a small-norm
    minorant would rule out bad translates for every short biased family (§2.4(i)). The two obstructions are dual.
  - This does **not** exclude Fourier-small minorants of *other* testable sets (e.g. tests that try lifts x̄ + jN).
- **E1 (pre-registered; §7).** The L¹ Fourier norm of S_can and S_bal is *smaller* than for random sets of equal size
  (r_A = 0.55–0.63 and 0.32–0.35; P1a [0.7, 1.3] **failed**), and the maximum bias is 5–72× larger (P1c passed). The
  norm still grows like |S|^{0.47} (P1b passed; random 0.51).
  - *Post hoc diagnosis* (`posthoc_e1.py`): the excess bias is at characters of *large* order, nearly trivial on the
    whole factor base (mean Re χ(p) over p ≤ B = 0.33–0.86). With 8–36 primes, Dirichlet approximation forces such
    characters. The factor base is far from an expander at these sizes: empirical λ = 0.57–0.99, while the GRH bound
    log²(pB)/√B = 31–67 is vacuous.
  - Under GRH at L-scale no character can be near 1 on all of the factor base, so this mechanism should disappear.
    E1's constant (A′ ≈ 0.47|S|^{1/2} for S_can) is consistent with Theorem FN's N^{1/2−o(1)}, which settles the
    L-scale question rigorously under GRH.

## 4. Approach 3: smooth compression

**Prop SC [proved, trivial].**
- *Statement.* A B-smooth u ≡ a (mod N), 0 < u ≤ N³, exists iff there are a B-smooth d ≤ N² and 0 ≤ j < B such that
  (least residue of a·d^{−1}) + jN is B-smooth.
- *Proof.*
  - (⇒) If u ≤ N², take d = u and c = 1. Otherwise multiply prime factors of u into d while d ≤ N². Then d > N²/B,
    c := u/d ≤ NB, c ≡ a·d^{−1}, so c is one of the first B lifts.
  - (⇐) u := d·c. ∎
- *Meaning.* Smooth compression is TH with the N^{2−o(1)} multipliers d^{−1} and a test of density ≤ B·δ. The N³
  freedom enlarges the multiplier family; it does not raise the per-test success probability above L^{O(1)}·δ. It helps
  existence proofs (GX needs N^{2+ε}) and does not help constructors that test candidates.
- **Relation to Prop SI / Theorem J / sihard.**
  - Prop SI asks for L^{o(1)}-cost constructors at L[1/3] scale. Theorem J says single-integer Johnson-class
    certificate decoders yield on average no more than random tests.
  - DSR is on a different axis. At L[1/2], random testing is affordable on average (δ = L^{−O(1)}), and the whole
    problem is worst-case determinism. Theorem J (an average statement) is therefore neither an obstruction to DSR nor
    a help.
  - Codex checked Boneh's CRT decoder with the full factor base: its agreement threshold ≈ √(2B ln N)/ln B exceeds the
    ≈ 3 ln N/ln B dyadic-interval prime divisors an integer ≤ N³ can have.
  - sihard's NP-hardness (Prop XS) is for explicit prime lists. DSR uses all primes ≤ B, its solutions are dense, and
    under GRH it is total. As in sihard Prop NC, NP-hardness would then give NP = coNP (Megiddo–Papadimitriou).

## 5. Structural results

**Prop EQ-dom [proved; ERH; folklore-level].**
- *Statement.* For fixed α, the following are equivalent up to constants:
  - (a) deterministic L[α, O(1)] complete factoring;
  - (b) a deterministic L[α, O(1)] algorithm giving, for each prime q ≤ B, a certified relation row ρ_q
    (∏ p^{ρ_{q,p}} ≡ 1 mod N) with |ρ_{q,q}| > Σ_{p≠q}|ρ_{q,p}| and L^{O(1)}-bit entries.
- *Proof.* (a)⇒(b): factor N, compute λ(N), output λ(N)e_q. (b)⇒(a): DD1 + Miller (Theorem DZ). ∎
- DSR ⇒ (b) by Theorem DZ. Codex adds that PIN on the primes < 2 ln² N is weaker still.
- So the target Theorem DZ minimally needs is deterministic factoring itself. DSR is a *sufficient* strengthening
  whose extra structure (smoothness) is what makes GRH existence (GX) available.

**Lemma PT (prime-modulus transfer) [proved, trivial; no prior statement found by us].**
- *Statement.* Fix α, b, c, C. Suppose DSR_α(b, c) holds, even only for admissible N all of whose prime factors exceed
  L_N[α, C]. Then for every large prime p there is a deterministic L_p[α, c + o(1)] algorithm that outputs, for every
  prime q ≤ B_p := L_{p}[α, b + o(1)], B_p-smooth u, v, w with |u|, |v| ≤ p^{3+o(1)}, Ω(w) ≤ (1 + o(1))log₂ p and
  u ≡ q^{H′}·w·v (mod p), for an explicit H′ ∈ [8 log₂ p, 9 log₂ p].
- *Proof.*
  - Let r be the least prime > max(L_{p²}[α, max(b, C)], 3 ln²(p²)). It is found deterministically by trial division in
    L^{O(1)}, r < 2·max(…) by Bertrand, and r < p for large p.
  - N := p·r is admissible: odd, two distinct primes, both > max(B(N), 3 ln² N), both > L_N[α, C].
  - Run DSR on N with targets q^{H_N} mod N, H′ := H_N = 8⌈log₂ N⌉. Reduce the output congruence mod p.
  - Since ln N = (1 + o(1))ln p, the bounds L_N = L_p^{1+o(1)} and N³ = p^{3+o(1)} follow. ∎
- *Meaning.*
  - The algorithm knows N's factorisation (r is tiny) and still has to solve the mod-p problem. Composite structure,
    and any bootstrapping through the factoring algorithm being built, cannot be what proves DSR.
  - DSR implies a deterministic, GRH-free-in-its-conclusion construction of smooth relations q^{H′}w v ≡ u in F_p^×
    for every p and every q ≤ B_p: the relation-collection step of index calculus with prescribed targets. Via DD1 this
    gives a full-rank relation sublattice for F_p^× deterministically. That alone does not give discrete logarithms:
    a multiple of p − 1 is trivially known, and saturation is missing.
  - We found no deterministic subexponential algorithm for this relation collection, under GRH or otherwise
    (literature §6).
  - Theorem DZ's recursion meets such N: preprocessing removes only factors ≤ max(B, 3 ln² N). Restricting DSR to N
    with all factors ≥ N^ε would escape PT, but DZ would then need deterministic removal of factors in (L^C, N^ε).
    The best known deterministic method costs N^{ε/4+o(1)} (Pollard–Strassen), which is exponential.
- **Is DSR implied by deterministic factoring?** Not known (codex agrees). PT shows DSR implies the prime-modulus
  problem, which factoring does not obviously touch. As far as we know, DSR is strictly stronger than what it is used
  for.

**Finite-N caveat (codex).** Read literally, with B = exp(b√(ln N ln ln N)) at small N, DSR outputs can fail to exist
(N = 1009·1013 with B = 2). DSR is an asymptotic statement ("for N ≥ N₀"). Theorem DZ only needs that.

## 6. Literature (approach 4)

Verified by me where marked (v); (c) reported by codex `dsr_1` with links, not re-verified by me.
- **Expanders under GRH.**
  - Jao–Miller–Venkatesan, JNT 129 (2009) 1491–1504; arXiv:0811.0647 (v): Thm 1.1, |λ| = O(n x^{1/2}log(xq)) on
    narrow ray class groups. Cor 1.3: walks of length ≥ C log|G|/log log q from any vertex land in any S with
    probability ≥ |S|/(2|G|).
  - The (Z/q)^× example uses primes ≤ (log q)^{2+δ} and p^{±1}.
  - Explicit GRH character-sum constants: Ernvall-Hytönen–Palojärvi (c); Lamzouri–Li–Soundararajan, Math. Comp. 84
    (2015) (c).
- **Derandomisation tools.** AKS (STOC 1987); Impagliazzo–Zuckerman (FOCS 1989); Nisan (Combinatorica 1992); INW
  (STOC 1994); Nisan–Zuckerman (JCSS 1996); Hoza–Pyne–Vadhan (ITCS 2021) (c); Bogdanov–Hoza–Prakriya–Pyne (CCC 2022)
  (c); Gopalan–Kane–Meka (FOCS 2015) (c); Koucký–Nimbhorkar–Pudlák (STOC 2011) (c); De (CCC 2011) (c);
  Rozenman–Vadhan (RANDOM 2005); Reingold (JACM 2008); Impagliazzo–Wigderson (STOC 1997). All are used only at the
  level of their stated seed bounds, in §2.3.
- **GRH derandomisation successes.** All rely on cosets or algebraic structure (Fourier-small sets or Galois
  structure), the case FC(c) / §3 isolates.
  - Ankeny 1952 and Bach 1990: least nonresidue, generation.
  - Shoup 1992: primitive roots.
  - Adleman–Lenstra STOC 1986 (v, via the Leiden record and secondary statements): deterministic polynomial-time
    construction of irreducible polynomials under GRH; unconditional for approximate degree.
  - Huang 1991; Rónyai 1992; Evdokimov 1994 (quasi-polynomial deterministic polynomial factoring under GRH);
    Ivanyos–Karpinski–Saxena 2009.
- **Closest analogue: deterministic prime finding.** Polymath4 / Tao–Croot–Helfgott, Math. Comp. 81 (2012) 1233–1246
  (v, abstract): randomized log^{O(1)} N; the best deterministic is N^{1/2+o(1)} (Odlyzko), also the RH-conditional
  bound. Their partial result decides the parity of π on intervals in N^{1/2−c}.
  - Same shape as DSR: dense, testable, trivially randomized, deterministically open even under RH.
  - Pseudodeterministic constructions (Oliveira–Santhanam FOCS 2017; Chen–Lu–Oliveira–Ren–Santhanam FOCS 2023, primes
    for infinitely many lengths) are randomized and input-length-indexed, so they do not apply to the input-dependent,
    worst-case DSR.
- **Path-finding in expanders.** Cayley and isogeny hash functions (Charles–Goren–Lauter, J. Cryptology 22 (2009))
  rely on rapid mixing coexisting with hard path-finding. DSR is deterministic short-path finding in the GRH expander
  Cay(G, primes ≤ B), with "path" meaning smooth representation. This is an analogy, not a reduction.
- **Smooth numbers.**
  - Soundararajan 2008 and Harper 2012 (moduli ≤ y^{4√e−ε}); Pascadi 2025 (level 5/8, averaged).
  - Bourgain–Glibichuk–Konyagin (subgroups ≥ p^ε): additive sums over subgroups, not smoothness correlations, and the
    sizes needed are p^{o(1)}.
  - Harman 1999 (unconditional smooth representatives); Boneh JCSS 2002 (CRT decoding).
  - None gives pointwise hits for L^{O(1)} explicit lists (codex `dsr_1` §1(d), consistent with zroute Obs. LV).
- **Sieve inputs for BOX-∞.** The fundamental lemma (Friedlander–Iwaniec, *Opera de Cribro*, Ch. 6;
  Halberstam–Richert); Brun–Titchmarsh (Montgomery–Vaughan 1973); Hildebrand's Ψ(x, y) asymptotics; Mertens.
- **Relation collection in F_p^×.** Pomerance 1987 (rigorous, randomized). We found no deterministic subexponential
  version, under GRH or otherwise.

## 7. Experiments (pre-registered; toy sizes; calibration only)

Pre-registration `notes/dsr_prereg.md` (`409641a`). Summary `experiments/analysis/dsr_summary.json`.

**Disclosed deviations.**
- (1) The rule-1 exception (float64 FFT in E1) was declared in the prereg. Parseval holds to 6.7e−16.
- (2) A shared scratchpad log was overwritten by another agent. Results were read from the JSONL files, which are
  authoritative.
- (3) D_slide has only |P|·t distinct multipliers (Sturmian windows: consecutive Beatty differences are 1 or 2). I did
  not anticipate this when registering. P2c/P2e therefore test that design defect, not "explicit walks" in general.
- No other deviations.

**E1 (Fourier-algebra norm; 18 moduli: primes 2^16–2^22, semiprimes 2^20, 2^22).**

| prediction | result | verdict |
|---|---|---|
| P1a: r_A ∈ [0.7, 1.3] | S_can 0.549–0.632; S_bal 0.315–0.353 | **FAIL** (below: smooth sets have smaller L¹ norm) |
| P1b: slope ln A′ vs ln|S| ∈ [0.4, 0.6] | 0.474 (random controls 0.510) | pass |
| P1c: r_∞ > 1 in ≥ 75% | 18/18 (5.2–12.9; S_bal 13.9–72) | pass |

- Post hoc: §3. The mechanism is factor-base characters near 1, with a factor base far from an expander at these
  sizes. It is not extrapolable.

**E2 (worst-case coverage of all of G; 30 semiprimes, 20/22/24 bits; S_can at b = 0.8, δ = 0.06–0.10; S_bal at
b = 0.5, δ = 0.02–0.07).**

| prediction | result | verdict |
|---|---|---|
| P2a: ν(R_G) ∈ [0.8, 1.5] in ≥ 90% | 60/60, ν ∈ [0.80, 1.28] | pass |
| P2b: median m*(R_walk)/m*(R_G) ∈ [0.7, 1.5] | 1.02 | pass |
| P2c: D_box, D_slide cover within cap | censored: D_box/S_can 10/30, D_slide/S_can 30/30, D_box/S_bal 29/30 | **FAIL** |
| P2d: median m*(D_box)/m*(R_G) ≥ 1.5 | 37.5 | pass |
| P2e: median m*(D_slide)/m*(R_G) ∈ [0.7, 1.5] | ∞ (never covers) | **FAIL** (deviation 3) |
| P2f: median ν(D_box) at 24 bits ≤ 2× at 20 bits | 23.1 → 40.3 (22 bits) → ∞ (24 bits) | **FAIL** |

- **Diagnosis [proved, post hoc].** D_box's worst case is the target a = ℓ (Prop BOX): m* ≥ Ψ(N/ℓ, B) + 1 = N^{1−o(1)}.
  The P2c/P2f failures are this exponential mechanism, not noise.
- **DSR targets only (descriptive).** Median T* = max_q first hit (20/22/24 bits):

  | list | T* | mean first hit |
  |---|---|---|
  | R_G | 34/50/59 | 10.2/12.8/15.4 |
  | R_walk | 41/54/71 | — |
  | D_orbit | 44/58/93 | — |
  | D_slide | 40/69/85 | — |
  | D_box (S_can) | 68/116/167 | 21/25/35 |
  | D_box (S_bal) | 1 663 / 3 304 / 9 894 | 553 / 1 229 / 3 180 |

  - The random targets q^H rarely hit the BOX set at these sizes, but adversarial N exist (§2.6(ii)).
  - D_box with the Gauss test is ≈ 100× worse even on q^H. This is zroute E4's inherited-lattice-vector mechanism:
    multiplying by small w mostly rescales the target's own short vector.

## 8. The exact remaining obstacle

Stated as the weakest form we can isolate that would give Theorem DZ, given Prop EQ-dom, BOX-∞ and PT:

> **(TH-wrap)** For fixed b and some fixed c, there is a deterministic L^{c}-time map N ↦ W_N (a list of pairs (m, test),
> where m is an explicit product of B-smooth integers and inverses with Ω ≤ log₂ N, and *wrapping*: height
> ≥ N^{1−o(1)}) such that for every admissible N ≥ N₀ and every prime q ≤ B, some (m, test) ∈ W_N certifies a smooth
> representation of q^H·m. By Lemma PT the same then holds modulo every prime p.

- *Why each word is needed.*
  - Bounded height is impossible by BOX-∞.
  - "Every N" rather than almost all is needed by Theorem DZ.
  - "Every q" rather than almost all targets is where Prop AA stops.
  - Spectral or GRH-walk information alone is insufficient by SB.
- *What could prove it.* A property of the explicit smooth sets that Theorem SB's random sets lack: a Fourier-small
  minorant on some testable set other than S_can / S_Gauss (Theorem FN), a conditional-expectation estimator, or a new arithmetic input on how q^H interacts
  with smooth numbers modulo N.
- *Known instances of this pattern:* none (deterministic prime finding is the analogous open problem).

## 9. Candidate cards

```text
ID: dsr-1 derandomised JMV walk (approach 1).
Mechanism: enumerate seeds of an explicit expander-walk / PRG / HSG; test canonical residue or Gauss split.
Cost: plain walk N^{1+o(1)} seeds (hitting) / N^{2+o(1)} (pointwise). Generic generators: seed >= log2((1-delta)phi).
Failure point: Theorem SB -- guarantees for classes containing all density-delta sets need (1-delta)phi endpoints;
  near-optimal Fourier uniformity of S does not help (SB(b)).
Status: refuted for generic derandomisers (proved); arithmetic-specific HSG OPEN; IW-conditional (Prop IWD).

ID: dsr-2 explicit bounded-height multiplier boxes / k-fold prime products (codex's missing lemma).
Status: REFUTED uniformly in N (Prop BOX-infinity, unconditional); almost-all-targets version proved (Prop AA, GRH),
  tight up to L^{1/(2b)}.

ID: dsr-3 Bach/coset certificates.  Needs a Fourier-small minorant.  REFUTED for S_can / S_Gauss and their dense
  minorants, every N, under GRH (Theorem FN: norm >= N^{1/2-o(1)} / N^{1/4-o(1)}).  Other test sets OPEN.
  E1: toy-scale norms below random (artefact of a non-expanding factor base).

ID: dsr-4 smooth compression (u <= N^3).  Prop SC: same hitting problem with N^{2-o(1)} inverse multipliers; no gain.

Structural: Prop EQ-dom (minimal DZ target == deterministic factoring); Lemma PT (DSR => mod-p relation collection).
```

## 10. Reviews

### 10.1 GPT-6.1-sol blind round (`codex/dsr_1`, xhigh)
- Independently obtained:
  - Prop RW (with an explicit loose GRH constant);
  - the tool-by-tool failures, adding HPV/BHPP, De, NZ and Cohen–Doron–Goldgraber 2026;
  - SB(a) (the complete-graph counterexample);
  - the Fourier-minorant and conditional-expectation sufficient conditions;
  - Prop IWD;
  - the Boneh threshold calculation;
  - the finite-N caveat;
  - "DSR not known to follow from factoring".
- Its proposed missing lemma (k-fold prime-product correlation) is refuted here by BOX-∞. Its probability: ≈ 10%.

### 10.2 GPT-6.1-sol adversarial round (`codex/dsr_2`) — pending.
### 10.3 Opus proof critic — pending.
### 10.4 Opus literature critic — pending.

## 11. Reproduction

```
PYTHONPATH=. .venv/bin/python -m src.dsr.harness e1 experiments/raw/dsr_e1.jsonl
PYTHONPATH=. .venv/bin/python -m src.dsr.harness e2 experiments/raw/dsr_e2.jsonl
PYTHONPATH=. .venv/bin/python -m src.dsr.analyze > experiments/analysis/dsr_summary.json
PYTHONPATH=. .venv/bin/python -m src.dsr.posthoc_e1 > experiments/raw/dsr_e1_posthoc.jsonl      # post hoc
PYTHONPATH=. .venv/bin/python -m src.dsr.posthoc_adv > experiments/raw/dsr_posthoc_adv.jsonl    # post hoc
.venv/bin/python -m pytest -q tests/test_dsr.py
```
E1 peaks at ≈ 1.24 GB RSS (22-bit 2-D FFTs); E2 ≈ 15 min.
