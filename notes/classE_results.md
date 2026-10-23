# classE: continued-fraction functionals at structured arguments, lattice-cone sums, reduced-form coefficients (phase 25)

Agent `classE`, 2026-10-23. Version v1, written before the adversarial round and the critics.

- Pre-registration: `notes/classE_prereg.md` (commit `2bf707a`), frozen before any labelled run.
- Code:
  - OPERATIONAL: `src/classE/{structure,fd,battery}.py`.
  - HARNESS/DIAGNOSTIC: `src/classE/{harness,toy3,posthoc,diag}.py`.
  - Tests: `tests/test_classE.py`.
- Raw data: `experiments/raw/classE_main_n{48,64,128}.jsonl.gz` and `classE_toy3_n{20,22}.jsonl`. Summaries are in
  `experiments/analysis/classE_*.json`.
- GPT-6.1-sol:
  - Blind round: `codex/classE_1.{md,out}`; the full report is copied to `codex/classE_1_report.md`. Several results
    below are codex's and are credited "(cx1 §…)".
  - Adversarial round: `codex/classE_2.{md,out}` (§7).
- Labels:
  - **PROVED**: the proof is here or in the cited cx section and was re-checked by us numerically. "No prior
    formulation found by us" is the strongest novelty claim.
  - **KNOWN**: with a reference.
  - Also used: **CONDITIONAL**, **HEURISTIC**, **experimentally supported**, **OPEN**, **UNVERIFIED**.

Notation:
- N = pq with p < q < 2p, as in the battery law D_n (both primes have the top two bits set).
- a/N = [0; a₁, …, a_r] with gcd(a, N) = 1.
- Remainders: r₀ = N, r₁ = a, r_{i+1} = r_{i−1} mod r_i.
- Convergents h_k/s_k (s₀ = 1, s_r = N).
- ε = (−1/p) on N ≡ 1 (mod 4); tb = [8p ≥ 7q].

---

## 0. Verdict (plain language)

**No leak was found. We also have no full floor theorem.** What we can now say:

1. **Where the factors could enter a CF path (PROVED).** The Euclid path of a/N is the lattice-reduction path of
   L_a = {(s, r) : a s ≡ ±r (mod N)}. Reduced mod p, this lattice is the line r ≡ ±a_p s. A computable statistic of
   the path can touch p in only two exact ways.
   - **Divisibility (Prop VIS; the exact iff is cx1 Thm 1).** A proper convergent denominator is divisible by p
     iff a = uq ± 1 together with a cylinder inequality. So the denominator then equals p, and no proper denominator
     is ever divisible by q.
     - For the Euler and Fermat residues this event is an Euler/Fermat-liar event at q, i.e. group order (M).
     - For ⌊√N⌋ it is impossible.
     - For a uniform random a it has probability 2(q−p−1)/((p−1)(q−1)) = O(N^{−1/2}).
     - The CF path is a two-dimensional lattice. It sees "a ≡ small mod q" only at height O(q/p) = O(1).
       Coppersmith's lattices reach N^{1/4}.
   - **Quadratic characters (Prop SPLIT).** Every path product satisfies
     ((−1)^k s_k r_{k+1} | p) = (a | p). That is one hidden bit, and it is known when a is a square with known root
     (2^{N−1}; 3^{(N−1)/2} when N ≡ 1 (mod 4); Fermat residues).
     - The products are ≥ N/(a_{k+1}+2). So the path generates "relations" of size about N. The CF of √N generates
       relations of size about 2√N.
     - HEURISTIC reading: the class-E genus channel behaves like a weaker CFRAC (families Q and M via smoothness).
       Products of size N^{1−c} need a_{k+1} ≥ N^c, which for random a has probability ≍ N^{−c} log N. We have not
       proved that statistics such as (s_k|p) alone are unbiased with respect to ε.
     - Reason: the form s(as − Nh) has *square* discriminant N², so it has no genus theory beyond (a|p). The
       principal form x² − Ny² of √N has discriminant 4N and carries the full genus.
2. **Floor theorems (PROVED, partial).**
   - **Random argument.** One end of the path, up to denominator N^{1/2−δ}, has a factor-free law within
     O(N^{−2δ}) in total variation (cx1 Thm 2; elementary). Both ends together, up to N^{1/8−δ}, are within
     O(N^{−4δ+o(1)}) (Weil). Beyond about N^{1/4} the two ends provably cannot be modelled as independent
     (cx1 §2.3). The full-path floor is OPEN.
   - **Structured archimedean arguments: Theorem F-SQ (new; no prior formulation found by us).** For a = ⌊√N⌋, every
     function of ⌊√N⌋ mod m₀, of t, and of the CF of (N mod ⌊√N⌋)/⌊√N⌋ up to denominator Q has the same law in
     both ε classes, up to O(m₀Q²N^{−1/6} log N). For tb the bound is m₀Q²N^{−1/12+o(1)}.
     - The proof uses Bombieri–Iwaniec's double large sieve on Σ χ(p)χ'(q) e(h√(pq)/m₀).
     - The path of ⌊√N⌋/N falls into three zones (corrected after the proof critic):
       - (a) convergent denominators ≤ N^{7/12−δ}: covered by F-SQ, for digit/matrix-state functionals;
       - (b) N^{7/12} to ≈ N^{3/4}: the rest of the forward half of CF(s′/r). This is OPEN, but it is an
         *archimedean* problem: equidistribution of √(pq) mod 1 in an ε class down to scale ≈ N^{−1/2};
       - (c) above ≈ N^{3/4}: the reversed CF of ±(N⁻¹ mod r)/r, which is a Kloosterman-type problem. OPEN.
   - **Exponent residues.** Only a conditional theorem: target-conditional rectangle equidistribution of
     (a/N, ā/N) at the readout's resolution suffices. At full depth this hypothesis is essentially the conclusion.
3. **Lattice-cone sums, d = 2 and 3.**
   - Exact evaluation (PROVED): poly-time for d = 2 via floor sums. For d = 3 we have an O(N) toy evaluator. A
     poly-time route via Barvinok exists (cx1 §3.5) but is not implemented.
   - Every sum has a public "Rademacher congruence": 12N·F₂ and 24N·F₃ mod N are explicit Laurent polynomials in
     the arguments and the shift (cx1 Thm 4; verified on 300 moduli).
   - A DSF analogue holds: the unweighted local part F_d(p) + F_d(q) factors N whenever C_d is a unit.
   - **The unshifted d = 3 sum is exactly three classical Dedekind sums**:
     F₃(0; a, b, c; N) = (N−1)/(8N) − ½ Σ s(a_i a_j⁻¹, N). Nothing new there.
   - The genuinely three-dimensional content is the shift-dependent part. In a toy exhaustive-population test
     (n = 20, 22) it showed no association.
4. **Reduced-form coefficients.** Reduction paths are quadratic-irrational CF paths with unbounded state (class H,
   classical). Composition histories and full coefficient access are computationally complete, so no unrestricted
   floor exists (cx1 §4.2; cf. leakmap REP-C).
   - **Random-class floor (KNOWN-based, references to be verified).** By Duke's theorem and its genus/subconvexity
     refinements (Michel 2004; Harcos–Michel 2006), the shape statistics of a uniformly random reduced form in a
     fixed genus coset equidistribute independently of the coset. A uniformly random class therefore leaks no genus
     bit beyond o(1).
   - The known Q-leaks are cusp excursions: a small leading coefficient is a Heegner point high in the cusp.
   - Short walks f^j are not covered by any theorem.
5. **Battery (pre-registered, 130 statistics, 2·10⁵ instances at each of n = 48, 64, 128): no hit for either target.**
   - The largest |z| among non-genus statistics is 3.15 (ε, n = 48) and 3.03 (tb, n = 48). It is ≤ 2.8 at
     n ≥ 64.
   - Genus-expected statistics: |z| ≤ 2.15.
   - Every positive control fires with 0 violations (z = 24–388 at n ≤ 64).
   - MDE |ρ| = 0.0153 (ε) and 0.0188 (tb, effective n).
   - One disclosed deviation: the registered V1 expectation formula is wrong for the size-predictable target tb. It
     passes with the corrected formula (§5).

**Most decisive open questions.**
- (1) Zone (b): target-conditional equidistribution of √(pq) mod 1 down to scale N^{−1/2+η}. The double large
  sieve stops at N^{−1/6}; this is an archimedean bilinear problem.
- (2) The *second half* of structured paths: the reversed CF of ā(N)/N, e.g. N⁻¹ mod ⌊√N⌋ over ⌊√N⌋, or the
  inverse of 2^{(N−1)/2} mod N, jointly with ε. This is a bilinear Kloosterman-type problem over primes with a
  moving modulus.
- Both together for ⌊√N⌋ would give the first full-depth floor for a structured class-E argument.

---

## 1. Q1: the structure of a/N at structured arguments

### 1.1 Identities and the two halves (KNOWN; cx1 §1.1, §1.9; tested)
- a s_k − N h_k = (−1)^k r_{k+1}, so a s_k ≡ ±r_{k+1} (mod N).
- N = s_k r_k + s_{k−1} r_{k+1} for k ≥ 1.
- **Mirror (Lemma TH).** The reversed word [0; a_r, …, a₁] equals s_{r−1}/N, with
  s_{r−1} ≡ (−1)^{r−1} a⁻¹ (mod N) and s_{r−1} = min(ā, N − ā). This is continuant symmetry.
- **Split.** Let m be maximal with s_m ≤ √N. Then N = s_m w_{r−m} + s_{m−1} w_{r−m−1}, where w are the reversed
  continuants. The path is:
  - a forward prefix of denominator ≤ √N, which is a function of the Farey cell of a/N;
  - one junction digit a_{m+1};
  - the reversal of a prefix of the CF of ā*/N of denominator < √N.
- So every class-E functional is A(a/N) + B(ā*/N) + junction, after transporting the matrix residues mod m₀
  through the terminal matrix (cx1 §1.9 warns that the bare reversal is wrong for state-dependent F).
- Test: `test_cf_identities_and_two_halves`.

### 1.2 Prop VIS (factor-bearing convergents). PROVED (necessity: ours; exact iff: cx1 Thm 1)
For p < q < 2p and gcd(a, N) = 1:
- p divides a proper convergent denominator s_k iff there are σ = ±1 and 1 ≤ u < p with a = qu + σ, such that
  w := (−σu⁻¹ mod p) < q − p.
- The denominator is then exactly p, with remainder p.
- q never divides a proper denominator.
- gcd(s_k, N) = gcd(r_{k+1}, N), so remainder exposure is the same event.

*Necessity.*
1. p | s_k forces p | r_{k+1}.
2. s_k r_{k+1} < N then gives (s_k/p)(r_{k+1}/p) < q/p < 2.
3. More generally, for unbalanced N the event needs a ≡ ±r′/s′ (mod q) with r′s′ < q/p.

*Sufficiency* is the cylinder inequality (cx1 (12)–(13)). The congruence alone is insufficient: cx1 gives
899 = 29·31, a = 63.

Tested: `test_prop_vis_necessity_and_planted` (3 000 cases, 0 necessity violations; planted a ≡ ±1 (mod q) hits
in ≈ 10%).

**Consequences.**
- (i) ⌊√N⌋ (unit case): never; the nonunit case is exactly q = p + 2 (cx1 §1.2).
- (ii) 2^{N−1}, 3^{(N−1)/2} (mod N): the event a ≡ ±1 (mod q) means ord_q(2) | p − 1 or the analogue. That is an
  Euler/Fermat-liar event at one prime, the group-order deviation of leakmap NE1 (family M). Battery: 0 path hits
  for these arguments at every size.
- (iii) Uniform a: exact probability 2(q−p−1)/((p−1)(q−1)) (cx1 (14)).
- (iv) **Comparison with Coppersmith.** Along the CF, p appears only when a mod q is a fraction of height
  < q/p = O(1). The unknown-divisor small-root method (Coppersmith; Howgrave-Graham) factors whenever a mod q lies
  within N^{1/4−η} of a known value. The 2-dimensional Euclid lattice is the weakest member of that hierarchy.
- (v) For a unit a, p/q itself is never a convergent: |a/N − p/q| < 1/q² forces a = p². A convergent h/s with
  |h/s − p/q| ≪ 1/s² needs a/N within ≈ 1/s² of p/q. This is the (A) covering problem
  (leakmap §4.3), not a CF effect.
- (vi) If a computable a(N) produced such an event for a positive fraction of N, scanning O(log N) gcds would
  factor. Poly-time constructors are therefore factoring algorithms.

### 1.3 Prop SPLIT (the genus content of class E). PROVED for the identity and bound; the "Meaning" is HEURISTIC/structural
- For every k with gcd(s_k r_{k+1}, N) = 1 and ℓ ∈ {p, q}: ((−1)^k s_k r_{k+1} | ℓ) = (a | ℓ).
- Moreover s_k r_{k+1} ≥ N/(a_{k+1} + 2).

*Proof.*
1. Multiply a s_k ≡ (−1)^k r_{k+1} by s_k: a s_k² ≡ (−1)^k s_k r_{k+1} (mod N).
2. N = s_{k+1}r_{k+1} + s_k r_{k+2} ≤ (a_{k+1}+1)s_k r_{k+1} + s_k r_{k+1}. ∎

Tested: `test_prop_split_genus_constraint_and_product_size`; 9 342 path products for 2^{N−1} and 3^{(N−1)/2}, 0
violations.

**Meaning.**
- F_a(s, h) = s(as − Nh) has discriminant N², a square. The values of a split form carry no character constraint
  except the one multiplicative constant (a|p).
- Compare the principal cycle of √(kN): A² − kNB² = (−1)^i Q_i, discriminant 4kN, gives (Q_i | p) = 1
  *individually* at every step. This is the genus source of the eps/leakmap n = 32 hits.
- For a known square a = x² (all our exponent residues have even exponent and known root), the class-E path is a
  **relation generator**: (xs_k)² ≡ (−1)^k s_k r_{k+1} (mod N), with right-hand sides of size N/(a_{k+1}+2). These
  are CFRAC relations of size ≈ N instead of ≈ 2√N.
  - A product reveals a fixed-conductor character only if it factors over small primes. That costs smoothness of
    numbers of size N^{1−o(1)}, unless a partial quotient of size N^{1/2+} occurs, which for random-like a happens
    with probability O(N^{−1/2} log N).
- For a with (a|N) = 1 and unknown root, all products share the one unknown bit (a|p) = (a|q). Their pairwise
  ratios are again relations.
- So the quadratic channel of class E is a weaker copy of the (Q)/(M) relation mechanism. This explains the
  pre-registered prediction P1 and the null `menu` battery at the square arguments.

### 1.4 Archimedean arguments. PROVED (structure)
- **a = ⌊√N⌋ = r (Prop SQ-CF).**
  - N = r(r+t) + s′ with 0 < s′ < r and t ∈ {0, 1} in the unit case (t = 2 ⟺ s′ = 0 ⟺ q = p + 2, a nonunit).
    So CF(r/N) = [0; r+t, w(s′/r)], where w(·) denotes the digit word of a rational in (0,1).
  - s′ = N mod r ≡ −uv (mod r), where u = r − p and v = q − r.
  - With δ = {√N}: s′/r = {2δ + δ²/r}. So outside a boundary layer of width 1/r, the CF of s′/r agrees with the CF
    of {√(4N)} = {2δ} through every cylinder whose endpoints are more than δ²/r (< 1/r) from {2δ}, with no wrap at 1.
    Typically this reaches denominators ≍ N^{1/4}. (Measured median
    agreement depth N^{0.234}.)
  - There is **no** uniform agreement: cx1's 3127 = 53·59 shows the boundary effect, and the recursion for √N
    differs.
  - The rest of the path is the reversed CF of (±N⁻¹ mod r)/r (Lemma TH applied to s′/r).
  - Also r² ≡ −(N − r²) (mod N) gives the single public quadratic datum ((r² − N | p) = 1) (cx1 §1.2).
  - Tests: `test_sqrt_argument_structure`.
- **a = ⌊N^{2/3}⌋.** CF = [0; C, CF(A/E)] with C ∈ {⌊N^{1/3}⌋, ⌊N^{1/3}⌋+1} (cx1 (7)). The forward half is the CF
  of the real number N^{−1/3} up to cylinder boundaries, i.e. archimedean (A).
- **a = ⌊N/k⌋ ± c, small k, c (Prop LIN).** One large quotient, then a bounded-modulus tail. Every bounded-state
  statistic is a function of N mod O_{k,c}(1) together with the large quotient, which is an affine floor of N
  (cx1 §1.4). Battery: the `div3` statistics are degenerate or size-only, which matches.
- **a = ⌊ωN⌋, fixed ω.** The forward half is the CF of ω (public; for bounded-type ω up to denominator c_A√N,
  cx1 §1.5). All N-dependence sits in the reversed half, ā*/N.
- **General principle.** For a = ⌊Nω(N)⌋ with ω a smooth real function, the forward half is archimedean and the
  backward half is the CF of an inverse residue ā/N. Hence structured archimedean arguments split into (A) data
  plus one "pseudo-random" inverse residue.

### 1.5 Exponent residues and GAR
- CRT components (cx1 (10)–(11)): 2^{(N−1)/2} ≡ (2/p)·2^{(q−1)/2} (mod p) and 2^{N−1} ≡ 2^{q−1} (mod p).
- By Lemma GAR: rep(x)/N ≡ (x_p − x_q)θ + x_q/N. A CF readout of rep(x)/N is a function of this phase. Holding d
  and x_q fixed, a depth-V cylinder resolves θ only to 1/(|d|V²).
- Reading θ to precision N^{−3/4} would need x_q to be known to within about N^{1/4}. That is exactly
  Coppersmith's regime on x − x_q ≡ 0 (mod q). Knowing x_q exactly is already gcd(x − x_q, N) = q.
- So GAR adds no channel beyond Prop VIS. The phase is visible only through knowledge of the CRT components that
  would factor anyway (cx1 §1.10). OBSERVATION; not a theorem that all readouts must go through θ.

---

## 2. Q2: floor theorems

### 2.1 Random arguments. PROVED (cx1 Thm 2, (22)–(23), §2.3; we re-derived (i) and (ii))
Let a be a uniform unit mod N. Let P_T(x) be the CF word stopped at the last denominator ≤ T, without the
overshooting digit.
- **(i) One end.** d_TV(P_T(a/N), P_T(U)) ≪ τ(N)T²/φ(N), where U is uniform on (0,1). The reference law does not
  depend on N.
  - Proof: P_T is constant on O(T²) intervals with Farey endpoints of order 2T, and
    #{units in I} = φ(N)|I| + O(τ(N)) by Möbius.
  - With fixed residue labels a mod m the same holds with error O(m·τ(N)T²/φ(N)), so remainder residues are covered
    too. Inversion permutes
    units, so the same holds for the reversed end, against the folded reference min(U, 1 − U).
- **(ii) Both ends.** Weil gives |S(h,k;N)| ≤ 4√(N·gcd(h,k,N)). The gcd terms need p | h and p | k (so max(|h|,|k|) ≥ p) and so drop out
  below the cutoff H = c√N < p. Hence the rectangle discrepancy of (a/N, ā/N) is ≪ N^{−1/2} log² N, and
  d_TV(joint T-prefixes, independent model) ≪ T⁴N^{−1/2} log² N, which is useful for T ≤ N^{1/8−η}.
- **(iii) Obstruction (cx1 §2.3).** For T ≥ N^{1/4+η}, the joint law is at TV distance → 1 from the independent
  model: N atoms cannot fill T⁴ cells. No endpoint-equidistribution argument can give a full-path floor. A
  full-path floor needs the dynamics: the transfer operator of the Gauss map with congruence twists, at a *fixed*
  modulus.
  - Baladi–Vallée (2005) and Bettin–Drappeau (2022) average over moduli and are not that.
- **Literature at fixed modulus** (cx1 §2.4; checked by cx against the papers; our literature critic to re-check):
  - Porter: the mean length has a uniform O(N^{−1/6+ε}) error.
  - Bykovskii–Frolenkov (Mat. Sb. 2017) improve this to a normalised error O(N^{−6/29+ε}) for balanced N.
  - Aistleitner–Borda–Hauke (arXiv 2210.14095): tails of Σ a_i and of the maxima at fixed N.
  - David–Kim–Mor–Shapira (arXiv 2401.15586) Thm 1.2: polynomial large deviations for the length and for fixed-word
    frequencies, with *ineffective* exponents.
  - None of these compares the full law of an arbitrary class-E functional.
- **Consequence for targets (cx1 (29)).** A reference law within ρ_N of the actual law gives
  |E[Y g] − E[Y B(C)]| ≤ E|Y|ρ_N. With a fresh uniform a, every bounded predictor reading only P_T(a/N) for
  T = N^{1/2−δ} has target advantage ≤ O(N^{−2δ}) beyond nothing at all: the reference law is N-independent.

### 2.2 Theorem F-SQ (structured archimedean argument ⌊√N⌋). PROVED (no prior formulation found by us)
**Setting.**
- p₁, p₂ are independent and uniform over the primes in I = [3P/4, P), with P = 2^{n/2}. This is the D_n
  sampler; conditioning on N ≡ 1 mod 4 and on p₁ ≠ p₂ changes only constants.
- A target class c is defined by product weights α(p₁)β(p₂): ε = ±1 means p₁ ≡ p₂ ≡ u (mod 4).
- Write y = √N/m₀ mod 1, r = ⌊√N⌋, δ = {√N}, t = ⌊(N − r²)/r⌋, s′ = N − r² − tr.

**Statement.** Let m₀ ≤ N^{1/6} and 2 ≤ Q ≤ √N. Let G be any [−1,1]-valued function of (r mod m₀, t,
P_Q(s′/r)). This includes the digit and matrix-state (γ_k mod m₀) contributions of every step of the path of
⌊√N⌋/N whose full convergent denominator is ≤ Q(r+t). Remainder-residue statistics also need N mod m₀. They are
covered by adding N mod m₀ to G's data and writing 1[p₁p₂ ≡ ν (mod m₀)] as ≤ φ(m₀) product weights, at a cost
polynomial in m₀ (critic finding 2; not written out). Then for the two ε classes c, c′:

  |E[G | c] − E[G | c′]| ≪ m₀ Q² N^{−1/6} log N.

For the archimedean target tb the bound holds with N^{−1/12+o(1)} in place of N^{−1/6}.

*Proof.*
1. **Reduction to a step function of y.**
   - s′/r = {ξ} with ξ = 2δ + δ²/r and |ξ − 2δ| ≤ 1/r.
   - The data (r mod m₀, t, P_Q({2δ})) form a step function g(y) with K ≤ m₀(2|F_{2Q}| + O(1)) ≪ m₀Q² pieces
     (δ ↦ {2δ} covers [0,1) twice): P_Q is constant between consecutive Farey fractions of order 2Q (cx1 Thm 2
     proof), and t jumps at 2δ = 1. In the unit case t ∈ {0, 1}.
   - The real data agree with g(y) except when y lies within O(1/(m₀ r)) of one of the K breakpoints.
2. **Discrepancy bound.** So the deviation is ≤ 2K·(D_c(y) + O(1/(m₀√N))), where D_c is the discrepancy of y in
   class c.
3. **Erdős–Turán.** D_c ≪ 1/H + Σ_{h≤H} h⁻¹ |E_c e(h√N/m₀)|.
4. **Double large sieve** (Bombieri–Iwaniec 1986, Lemma 2.4; Iwaniec–Kowalski §7.? — exact numbering UNVERIFIED).
   - Use x_p = √p − c₀ (centred, |x_p| ≤ X ≍ √P) and y_q = (h/m₀)√q; the factor e(c₀ y_q) goes into β_q.
     Then (1+XY) ≍ hP/m₀.
   - Spacing: |x_p − x_{p′}| ≥ |p−p′|/(2√P). Within 1/Y there are ≤ 1 + 4m₀/h primes, and within 1/X there are
     ≤ 1 + O(m₀/h).
   - Hence |Σ_{p,q} α_p β_q e(h√(pq)/m₀)| ≪ (hP/m₀)^{1/2} π(I)(1 + m₀/h).
   - Divide by π(I)² ≍ (P/log P)², with √P = N^{1/4}: |E_c e(h√N/m₀)| ≪ (log N) N^{−1/4}(h^{1/2}m₀^{−1/2} + m₀^{1/2}h^{−1/2}).
5. **Summing over h ≤ H = N^{1/6}:** D_c ≪ N^{−1/6} log N + m₀^{1/2}N^{−1/4} log N.
6. **Target tb.** The region 8p ≥ 7q is not a product.
   - Split G = G⁺ − G⁻.
   - Sandwich 1[log p₁ − log p₂ ∈ J] between Beurling–Selberg majorant and minorant trigonometric polynomials of
     degree T in log p₁ − log p₂. Each term is a product weight.
   - The sandwich gap has prime-pair mass O(1/T). This needs primes in intervals of length P/T; Huxley's
     short-interval prime number theorem covers lengths P^{7/12+ε}, so T ≤ P^{5/12−ε}.
   - With T = N^{1/12} this gives D_tb ≪ N^{−1/12+o(1)}. (The 1/k decay of the coefficients allows
     T = N^{1/6}, i.e. N^{−1/6+o(1)}; not written out.) ∎

**Numerical check (proof critic).** An exact computation over all prime pairs at n = 40, m₀ = 1 gives
|E_c e(h√N)| = 7.2·10⁻³, 1.7·10⁻³ and 3.6·10⁻⁴ for h = 1, 16, 1024, against the scale N^{−1/4} = 2.0·10⁻³. The
observed decay is ≈ N^{−1/4}h^{−1/2}, and the two mod-4 classes agree to within noise.

**Remarks.**
- F-SQ covers zone (a), the path of ⌊√N⌋/N up to full denominators ≈ N^{7/12−δ}.
  - The CF of ⌊√N⌋/N starts with the digit r + t ≈ √N.
  - The CF of s′/r then has its own forward half up to inner denominators ≈ N^{1/4}, i.e. full denominators
    ≈ N^{3/4}.
  - Zone (b), from N^{7/12} to N^{3/4}, is archimedean but beyond the double large sieve's N^{−1/6} resolution.
- The same method handles ⌊√(kN)⌋, and ⌊N^{2/3}⌋ via phases h(pq)^{1/3}: for the latter the spacing gives
  D ≪ N^{−1/6+o(1)} (sketch, not written out).
- Zone (c), the reversed half, is **not** covered: the CF of (±N⁻¹ mod r)/r, i.e. Kloosterman fractions with the
  moving modulus ⌊√(pq)⌋. Equivalently, it is the reversed CF of ±r̄/N up to denominator N^{5/12}. These are the
  open problems of §0.

### 2.3 Conditional theorem for exponent residues (cx1 Thm 3). CONDITIONAL
- **Hypothesis H(Δ).** In each target class, the law of (a(N)/N, ā(N)/N) and of finitely many residue labels
  differs from one common reference by ≤ Δ on every rectangle. (Use the correct reference: for ⌊ωN⌋ the first
  coordinate is degenerate.)
- Then for every event that is a union of ≤ R rectangles (with sandwich error η), the target-class difference is
  ≤ 2(RΔ + η).
- For two T-prefixes, R = O(T⁴).
- For full-depth class-E functionals R is ≈ N cells, so H must hold at resolution 1/N with Δ ≪ N^{−1−}. That is
  essentially "a(N) is indistinguishable from a random unit by cylinder tests", i.e. the conclusion itself.
- For 2^{(N−1)/2} we know no target-conditional equidistribution even at coarse resolution. It would need sums
  Σ_{p,q} e(h(2^{(q−1)/2} q̄/p ± 2^{(p−1)/2} p̄/q)): exponential functions in the exponent mixed with inverses.

---

## 3. Q3: lattice-cone (Fourier–Dedekind) sums in dimension d = 2, 3

F_d(n; a; M) = (1/M) Σ_{k=1}^{M−1} e(kn/M) / ∏_j (1 − e(k a_j/M)), with gcd(a_j, M) = 1.

### 3.1 Exact evaluation. PROVED
- **Formula.** Expand 1/(1−x) = −(1/M) Σ_{r<M} r x^r for x^M = 1, x ≠ 1. This gives
  F_d = (−1)^d M^{−d−1}[M·C − (M(M−1)/2)^d], where C = Σ_{r ∈ [0,M)^d, n + Σa_j r_j ≡ 0} ∏ r_j.
  This is `fd.py`, checked against the brute-force sum and floating point (`test_fourier_dedekind_exact_vs_brute_and_float`).
- **d = 2.** Normalising a₂ = 1 gives C = Σ_r r·((c₀ + c₁ r) mod M). This is a weighted floor sum, computed by the
  Euclid-type `floor_sums` recursion in O(log M) (poly-time).
- **d = 3.**
  - Our evaluator is O(M): G(c) = Σ_{r₂} r₂((c − βr₂) mod M) obeys
    G(c+1) = G(c) + M(M−1)/2 − M(β⁻¹(c+1) mod M). It takes 0.34 s at M ≈ 2²².
  - cx1 (44): F₃ = −M⁻² Σ_{x,y<M} (x−h)(y−h)⌊(γ+αx+βy)/M⌋ with h = (M−1)/2. This is a weighted lattice count in a
    3-dimensional polytope, so it is polynomial-time by Barvinok (1994) with polynomial weights (De Loera et al.).
    It was not implemented here.

### 3.2 Prop FD3-0: the unshifted d = 3 sum is classical. PROVED
- F₃(0; a, b, c; M) = (M−1)/(8M) − ½[s(ab⁻¹, M) + s(ac⁻¹, M) + s(bc⁻¹, M)].
- *Proof.* (1 − e(x))⁻¹ = ½ + (i/2)cot(πx). Expand the product. Terms with an odd number of cotangents cancel
  under k ↦ −k. The pair terms are (1/4M) Σ_k cot cot = s(a_i a_j⁻¹, M). ∎ (Also cx1 (37).)
- Tested on 150 random moduli (`test_d3_unshifted_is_three_dedekind_sums`).
- Corollaries:
  - The unshifted d = 3 sum is class E at the three ratio arguments, and newfamily's DSF applies to it.
  - The first genuinely higher unshifted object is Zagier's 4-cotangent sum (d = 4).
  - For d = 3 the new content is in the shift: F₃(n) − F₃(n + c) = F₂(n) (cx1 (39)), and the n-dependence carries
    the odd-cotangent terms Σ_k sin(2πkn/M) ∏_{j≤3} cot(πka_j/M).

### 3.3 Rademacher-type congruences. PROVED (cx1 Thm 4; we verified on 300 random moduli, all M including even and 3 | M)
- 12M·F₂(n; a, b; M) ≡ (−a² − b² − 3ab + 6n(a+b) − 6n²)/(ab) (mod M).
- 24M·F₃(n; a, b, c; M) ≡ (4n³ − 6An² + 2(A²+B)n − AB)/(abc) (mod M), with A = Σa and B = Σ_{i<j} a_i a_j.
- Our independent fit of the d = 3 hyperplane sum, C/M² ≡ (1/24)Σ_{i≠j}a_iā_j − (n/12)Σ_i a_i∏_{j≠i}ā_j −
  (n/4)Σā_i + (n²/4)Σ_{i<j}ā_iā_j − (n³/6)ā₁ā₂ā₃, is the same statement (`test_aggregate_congruences_d2_d3_shifted`).
- **Meaning.**
  - N·F_d mod N is an aggregate: a Laurent polynomial in the public a_j, ā_j and n mod N.
  - The non-aggregate content is the integer J_d = (D_d N F_d − rep(C_d))/N, with D₂ = 12 and D₃ = 24.
  - For d = 2 and n = 0, J₂ is the CF (Barkan–Hickerson–Knuth) functional, i.e. Rademacher's Φ. For a moving
    shift n(N), the Euclid recursion carries an O(N)-size shift, and membership in the narrow class E is not
    established (cx1 §3.4).
  - For shifted d = 3 it is the Barvinok/lattice-count functional (44).

### 3.4 DSF analogue. PROVED (cx1 §3.3)
- **Unweighted local part.** For L_d = F_d(p) + F_d(q) and X = D_d N L_d ∈ Z:
  X ≡ C_d(p + q) (mod N), so p + q = X·C_d⁻¹ mod N whenever C_d is a unit.
- **Root-order local part.** The local part of the frequency split, F_d(p)/q + F_d(q)/p, is a different oracle.
  At fixed small arguments it factors trivially; for example F₂(0; 1, 1; M) = −(M−1)(M−5)/(12M) gives a quadratic
  in p + q.
- **Degenerate case.** F₃(2; 1, 1, 2; M) ≡ 0, so some local parts are identically zero.
- **Spatial vs frequency splits.** For Dedekind sums the sawtooth split s(a,N) = s(a,p) + s(a,q) + U (newfamily
  DS1) and the frequency split s(a,N) = s(a,p)/q + s(a,q)/p + U′ are different. "The local part" is not canonical;
  only the first factors via the Rademacher congruence. (This is cx1 §3.1, including the Dedekind special
  case is ours.)

### 3.5 Toy test of the shifted d = 3 sums (pre-registered §6; exhaustive populations)
- **Populations.** Exhaustive: n = 20 has 324 semiprimes (10-bit primes, N ≡ 1 mod 4) and n = 22 has 1 089.
- **Cost.** The exact O(N) evaluator takes 0.36 s per instance (4 sums) at n = 20 and 1.5 s at n = 22. The whole
  run took 9.8 min on 3 workers.
- **Results.** 16 statistics × 2 targets per size; one statistic is degenerate in each cell.
  - max |z| = 1.64 (ε) and 1.80 (tb) at n = 20.
  - max |z| = 1.66 (ε) and 1.86 (tb) at n = 22.
  - Mean z² is 0.68–1.07. There is no hit at the registered |z| ≥ 4.5.
- **Status.** A sanity check only, as registered. The MDE is |ρ| ≈ 0.35 at n = 20 and 0.19 at n = 22, primes are
  heavily reused, and the populations are tiny.
- It shows only that the shifted d = 3 sums carry no *strong* toy-size signal. Within the operational evaluator
  this is consistent with the exact structure: an aggregate congruence (§3.3) plus the integer J₃.

---

## 4. Q4: reduced-form coefficients (genusleak's open class)

### 4.1 What is proved
- **Reduction is class H** (KNOWN; cx1 §4.1). Indefinite reduction of a form of discriminant 4kN is the CF of a
  quadratic irrational of that discriminant, with unbounded state (P_i, Q_i). The principal cycle satisfies
  A² − kNB² = (−1)^i Q_i, which gives the genus constraints.
- **Positive definite (−4kN).** Gauss reduction is reduction in the upper half-plane, not a rational a/N path.
- **Completeness** (cx1 §4.2; cf. leakmap REP-C, genusleak §6.2). Coefficients together with the discriminant
  retain N, so unrestricted post-processing computes anything. A floor needs a restricted statistic grammar. The
  natural one is bounded-state additive statistics along a polynomial-length walk, which is what the battery tests.
- **Split vs non-split (Prop SPLIT, §1.3).** The structural difference between class E and class H is the
  discriminant of the form whose values the path visits: N² (split, one character bit) against 4kN (non-split,
  the full genus).
  - Q-leaks need small values. In class H, Q_i ≤ 2√(kN) with smooth values at rate about u^{−u}. In class E the
    values are ≥ N/(a+2).
  - HEURISTIC/structural: this is the sense in which "class E plus Q" covers the CF/form paths. It covers their *character*
    content. It does not exhaust non-multiplicative content; the battery tests that empirically.

### 4.2 Random classes: Heegner equidistribution. KNOWN-based (references UNVERIFIED in detail)
- A reduced form (a, b, c) of discriminant D = −4kN corresponds to the Heegner point
  τ = (−b + √D)/(2a) ∈ SL₂(Z)\H.
- Duke (Invent. Math. 92, 1988) proved that these points equidistribute as |D| → ∞, for fundamental D.
- Restricted to a genus coset, equidistribution follows from Weyl sums twisted by genus characters. By
  Gross–Zagier/Waldspurger-type formulas (Zhang), these factor as L(½, φ⊗χ_{d₁})L(½, φ⊗χ_{d₂}) with d₁d₂ = D.
  Subconvexity for quadratic twists (Blomer–Harcos; Michel 2004 Annals; Harcos–Michel 2006 Invent.) then bounds
  them. Level structure (a, b mod m₀) is the same on Γ₀(m₀) or Γ(m₀).
- **Consequence.** For a uniformly random class in a fixed genus coset, the law of every bounded "shape" statistic
  (b/2a, a/√|D|, a and b mod m₀) is the same for all cosets up to |D|^{−η}.
  - So shape statistics of random classes leak no genus character beyond o(1).
  - The Q-leaks of eps/leakmap are *cusp excursions*: a small first coefficient a, i.e. Im τ = √|D|/(2a) large.
    This is a region of measure ≍ a/√|D|.
- **Not covered.** Short orbits {f^j : j ≤ poly(n)}. These are Linnik-type "short orbit" equidistribution problems
  and are OPEN. Battery FW tests j ≤ 32.

---

## 5. Q5: battery (pre-registered; `notes/classE_prereg.md`)
**Setup.**
- 2·10⁵ instances of D_n at each n ∈ {48, 64, 128}.
- Seeds: main = 25.52·10¹² + n·10⁹.
- 130 statistics: 113 `menu`, 5 `menuG`, 5 NC, 3 PC, plus the 4 degenerate NC.
- Stratified permutation z. Strata: N mod 2⁸ for ε; the top 12 bits of N for tb. A hit is |z| ≥ 6.
- Cost: 1.9–3.9 ms per instance. Wall time 140 / 180 / 280 s on 3 workers.

| n | target | P(y=1) | V1 planted z (registered exp. / corrected exp.) | V2 PCs: fires / violations / z | V3 shuffle mean z² | menu mean z² | max \|z\| menu | max \|z\| menuG | hits | MDE \|ρ\| |
|---|---|---|---|---|---|---|---|---|---|---|
| 48 | ε | 0.498 | 21.2 (22.4 / 22.4) | pc2 25 800/0/173; pc3 85 300/0/388 | 1.28 | 1.10 | 3.15 (`cf3_sumbits`) | 2.15 | 0 | 0.0153 |
| 48 | tb | 0.715 | 16.1 (22.4 / 18.2) | fermat 31 734/0/105 | 0.93 | 1.27 | 3.03 (`isqrt2_ms_mod3`) | 0.94 | 0 | 0.0188 |
| 64 | ε | 0.499 | 21.2 (22.4 / 22.4) | pc2 1 858/0/43.3; pc3 765/0/27.8 | 0.85 | 0.99 | 2.27 | 1.84 | 0 | 0.0153 |
| 64 | tb | 0.712 | 16.8 (22.4 / 18.2) | fermat 2 052/0/24.3 | 1.03 | 0.97 | 2.78 (`cbrt2_n_a1`) | 1.21 | 0 | 0.0188 |
| 128 | ε | 0.500 | 21.2 (22.4 / 22.4) | never fire (P3) | 1.14 | 0.95 | 2.62 | 0.62 | 0 | 0.0153 |
| 128 | tb | 0.713 | 18.3 (22.4 / 18.2) | never fire (P3) | 0.68 | 0.96 | 2.75 (`div3_*`) | 1.87 | 0 | 0.0188 |

**Validity.**
- V2: every PC has 0 violations, and at n ≤ 64 every one fires with z ≥ 24.
- V3: the shuffle mean z² lies in [0.68, 1.28] in every cell.
- V1:
  - ε passes, with the same single flip pattern as newfamily, hence 21.2.
  - tb fails the *registered* expectation at n = 48 and n = 64 and passes the corrected one in all three cells
    (deviation D1, §9).
  - The cause is that the top bits of N explain part of tb's variance. The plumbing is correct.
- The tb MDE uses the effective n = 4ΣS_yy ≈ 1.32·10⁵.

**Conclusion (registered rule).**
- **Refuted for this menu** at n ∈ {48, 64, 128} for ε (|ρ| ≥ 0.0153) and for tb (|ρ| ≥ 0.0188, with the D1
  qualification).
- This is a power statement for *linear* association beyond the floor. Categorical codes (mod 7, mod 5, mod 3)
  are tested linearly only.
- P1–P4 were all met:
  - P2: no menuG hit, even at n = 48. The genus statistics of √(3N) at 256 steps are too weak at this size.
  - P4: 0 nonunit arguments and 0 Euclid-path gcd hits at every size.
    - The only splits are 3 ambiguous forms on the n = 48 form walks, f^j with ord f | 2j ≤ 64, a small-order
      (family M) event. They are STRUCTURAL, with rate 1.5·10⁻⁵, and there are none at n ≥ 64.

**Post hoc (`posthoc.py`; not registered as tests).**
- tb under finer floors (top 16 bits; top 8 bits × N mod 16): max |z| 3.35 / 3.06 (n = 48), ≤ 2.74 (n ≥ 64).
- **Secondary view, tb against N mod 2⁸ strata only.** Many size statistics "hit": `*_dr1_bits`, the `div3`
  digits, up to |z| = 80.7. This is the registered size confounding (prereg §1) and is not a finding: those
  statistics are functions of the fine size of N, which predicts tb (cf. leakmap's top-bits leak). Under the size
  floor they are null.
- Top |z| per cell is in `experiments/analysis/classE_posthoc.json`. No statistic exceeds |z| = 3.35 in any
  primary or finer-floor cell.

---

## 6. Candidate cards (rule 9, condensed)

| Mechanism | Model / scope | Obstruction (precise) | Status |
|---|---|---|---|
| Factor-bearing convergents of a/N | any a(N) | iff a = uq ± 1 with a cylinder inequality (VIS). Exponent residues: Euler/Fermat liar at q (M). ⌊√N⌋: impossible | refuted for this model (PROVED) |
| θ via rep(x)/N (GAR) | CF of CRT residues | phase visible only with CRT-component knowledge ⇒ Coppersmith/gcd | observation |
| Quadratic characters on the class-E path | any a | split discriminant N² ⇒ one bit (a\|p); relations of size ≥ N/(a+2) | Q/M via smoothness (PROVED structure) |
| Random-argument class E | prefixes ≤ N^{1/2−δ}; joint ≤ N^{1/8} | TV bounds (cx1 Thm 2, Weil) | floor PROVED (truncated); full path OPEN |
| ⌊√N⌋, first N^{7/12} of path | digit/γ-state functionals, zone (a) | double large sieve (F-SQ) | floor PROVED |
| ⌊√N⌋ zone (b) | √(pq) mod 1 at scale N^{−1/2} | double large sieve stops at N^{−1/6} | OPEN |
| ⌊√N⌋ / exponent residues, reversed half | Kloosterman fractions with a moving modulus | none known | OPEN |
| d = 2, 3 Fourier–Dedekind | exact | aggregate congruence + J_d; unshifted d = 3 = Dedekind; DSF analogue | no leak; d = 3 shifted toy-null |
| Reduced-form walks | bounded-state statistics, j ≤ 32 | random-class floor (Duke/Michel); short walks open | battery null |

## 7. GPT-6.1-sol rounds
CODEX_PLACEHOLDER

## 8. Opus critics
CRITIC_PLACEHOLDER

## 9. Deviations from the pre-registration
- **D1.** The registered V1 expectation ρ√n assumes within-stratum target variance ¼. For tb under the top-12-bit
  floor, the size of N explains part of the variance, so the planted z is 2ρ√(Σ_h S_yy,h) (`posthoc.py`; derivation
  in its docstring).
  - By the registered formula, V1 fails for tb at n = 48 (16.1 vs 22.4 ± 3) and n = 64 (16.8).
  - By the corrected formula (18.2) it passes at every size (16.1, 16.8, 18.3).
  - The tb MDE is reported with the effective n = 4ΣS_yy.
- **D2.** `harness.py` was changed to stream records to disk before the main run (commit after the prereg; no
  change to statistics or analysis).
- The post-hoc finer floors (top 16 bits; top 8 bits × N mod 16) were run for tb although there was no hit
  (registered only as an audit step).

## 10. Reproduction (project root)
```
.venv/bin/python -m pytest -q tests/test_classE.py
PYTHONPATH=. .venv/bin/python -m src.classE.harness --n 48 64 128 --range main --count 200000 --workers 3
PYTHONPATH=. .venv/bin/python -m src.classE.posthoc 48 64 128
PYTHONPATH=. .venv/bin/python -m src.classE.toy3 --n 20 22 --workers 3
```
