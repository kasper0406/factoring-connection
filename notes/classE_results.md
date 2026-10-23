# classE: continued-fraction functionals at structured arguments, lattice-cone sums, reduced-form coefficients (phase 25)

Agent `classE`, 2026-10-23. Version v2: revised after the codex adversarial round (MAJOR REVISION, all findings applied), the Opus proof critic (ACCEPT WITH FIXES) and the Opus literature critic (OK WITH FIXES).

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

**No leak was found, and no full floor theorem was proved.** The registered battery had no main-menu hit. What we
can say precisely:

1. **Exact algebraic contacts between a CF path and the factors (PROVED identities; their exhaustiveness is NOT
   claimed).** The Euclid path of a/N reduces the lattice {(s, r) : a s ≡ r (mod N)}, with signed remainders.
   Reduced mod p, this lattice is the line r ≡ a_p s. Two families of exact consequences follow.
   - **Divisibility (Prop VIS; the exact iff is cx1 Thm 1).** A proper convergent denominator is divisible by p
     iff a = uq ± 1 together with a cylinder inequality. The denominator is then exactly p, and q never divides a
     proper denominator.
     - For the exponent residues the event is an Euler/Fermat-liar event at q (group order, M).
     - For ⌊√N⌋ it is impossible.
     - For a uniform a its probability is exactly 2(q−p−1)/((p−1)(q−1)) = O(N^{−1/2}). At n = 48 this predicts
       ≈ 0.02 events in the whole battery; 0 were observed.
     - The two-dimensional Euclid lattice sees "a ≡ small (mod q)" only at height < q/p = O(1). Coppersmith-type
       lattices reach N^{1/4}.
   - **Power-residue certificates (Prop SPLIT and its order-d form, the latter due to codex_2).**
     W_{k,d} := (−1)^k s_k^{d−1} r_{k+1} ≡ a·s_k^d (mod N). So χ(W_{k,d}) = χ(a) for every character of order
     dividing d at p or q.
     - For d = 2: ((−1)^k s_k r_{k+1} | p) = (a | p), one hidden bit. It is known when a has a known square root
       (2^{N−1}; 3^{(N−1)/2} on N ≡ 1 (mod 4); Fermat residues). 2^{N−1} even has a known 4th root.
     - The products satisfy N/(a_{k+1}+2) ≤ s_k r_{k+1} < N. The CF of √(kN) gives (−1)^i Q_i ≤ 2√(kN) instead.
     - HEURISTIC reading (not a theorem): read as relations, class-E paths at square arguments are CFRAC-like with
       values of size N/(a_{k+1}+2). Such values are small only after a huge partial quotient. The form s(as − Nh)
       has square discriminant N², so ordinary genus theory does not apply. That is a structural remark, not a
       lower bound.
2. **Floor theorems (PROVED, partial).**
   - **Random argument.**
     - One end of the path, up to denominator T = N^{1/2−δ}: TV distance O(N^{−2δ}) from a reference that depends
       only on T (cx1 Thm 2).
     - Both ends together, up to T = N^{1/8−δ}: TV distance O(N^{−4δ+o(1)}) (Weil; cx1 (23)).
     - Beyond T ≈ N^{1/4} no independent-endpoints model can work (cx1 §2.3).
     - The full-path floor is OPEN.
   - **Structured archimedean argument, Theorem F-SQ (no prior formulation found by us).** For a = ⌊√N⌋:
     - Any bounded function of (⌊√N⌋ mod m₀, t, the CF word of (N mod ⌊√N⌋)/⌊√N⌋ up to denominator Q) has the
       same law in both ε classes up to O(m₀Q²N^{−1/6} log N). The proof uses the Bombieri–Iwaniec double large
       sieve on Σ α(p)β(q) e(h√(pq)/m₀).
     - This controls the tail word and matrix states up to full denominators N^{7/12−δ}. It does *not* control
       arbitrary functions of the huge first digit ⌊√N⌋ + t.
     - Not covered (OPEN): zone (b), the archimedean range N^{7/12}–N^{3/4}; and zone (c), the reversed CF of
       (±N⁻¹ mod ⌊√N⌋)/⌊√N⌋.
   - **Exponent residues.** Only a conditional theorem (cx1 Thm 3). At full depth its hypothesis is essentially the
     conclusion.
3. **Lattice-cone sums, d = 2 and 3.**
   - Exact evaluation (PROVED): poly-time for d = 2 via floor sums. For d = 3 we have an O(N) toy evaluator.
     Barvinok gives a poly-time route (cx1 §3.5), not implemented.
     - Remark: F₃(n) − F₃(n + c) = F₂(n) (cx1 (39)) plus a short vector of the lattice
       {(j₁, j₂, j₃) : Σ j_i a_i ≡ n} would give an O(N^{1/3+o(1)}) evaluator. This is not implemented.
   - Public "Rademacher congruences" (cx1 Thm 4; verified by us on all moduli): 12N·F₂ and 24N·F₃ mod N are
     explicit Laurent polynomials. A DSF analogue holds: the unweighted local part F_d(p) + F_d(q) factors N
     whenever C_d is a unit.
   - **The unshifted d = 3 sum is three classical Dedekind sums** (found independently here and in cx1 (37)).
     Telescoping makes every shift n ≡ j·a_i with small j classical too.
   - **Our registered toy test of shifted d = 3 sums used shift 1 with third argument 1. It is therefore a test of
     classical Dedekind combinations, not of 3-dimensional content (codex_2).**
     - A disclosed post-hoc rerun with shift (N−1)/2, which does not telescope in few steps, is in §3.5.
     - Both are null but tiny (MDE |ρ| ≈ 0.2–0.35).
4. **Reduced-form coefficients.**
   - Reduction paths are quadratic-irrational CF paths with unbounded state (class H, classical).
   - Coefficient access plus post-processing retains N, so no unrestricted floor exists.
   - Duke-type equidistribution (Duke 1988; Harcos–Michel 2006 for genus cosets of fundamental D) gives only weak
     convergence against fixed regular test sets. It is **not** a floor for arbitrary bounded statistics: at fixed
     D distinct genus cosets have disjoint supports (codex_2 found this FATAL in v1; withdrawn).
   - Short walks f^j are uncovered. Our battery's walk statistics turned out to be largely uninformative (§5): for
     small j, f^j = (ℓ^j, ·, ·) is public.
5. **Battery (pre-registered; 130 statistics, of which 117 `menu` and 5 `menuG`; 2·10⁵ instances at each of
   n = 48, 64, 128).**
   - No |z| ≥ 6. The largest is 3.15 (ε) and 3.03 (tb); at n ≥ 64 it is ≤ 2.8.
   - Every positive control fires with 0 violations (z = 24–388 at n ≤ 64).
   - V1 fails for tb at n = 48 and 64 under the *registered* expectation formula, which is wrong for a
     size-predictable target. It passes under the corrected first-order formula (deviation D1).
   - MDE |ρ| ≈ 0.0153 for linear association beyond the floor. For tb the feature-specific value lies between
     0.0154 and 0.0188.
   - Several null features are uninformative by construction: `fw_*` at n ≥ 64, `fw_minbits`, `fw_sq`, `div3_*`,
     and `dr1_*` (classical Dedekind).

**Most decisive open questions.**
- (1) Zone (b) of ⌊√N⌋: target-conditional equidistribution of √(pq) mod 1 at scale N^{−1/2+η}.
- (2) Second halves of structured paths: the reversed CF of ā(N)/N (N⁻¹ mod ⌊√N⌋, or the inverse of
  2^{(N−1)/2} mod N), jointly with ε. This is a bilinear Kloosterman-type problem over primes with a moving modulus.
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
  < q/p = O(1). The unknown-divisor small-root method (Coppersmith; Howgrave-Graham) recovers a small root x₀ with
  a − x₀ ≡ 0 (mod q) and |x₀| < N^{1/4−η}; it factors when gcd(a − x₀, N) is proper (it can be N, e.g. for a known
  global integer root). The 2-dimensional Euclid lattice is the weakest member of that hierarchy.
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

Tested: `test_prop_split_genus_constraint_and_product_size` (fermres and a random argument), plus a scratch check of 9 342 path products for 2^{N−1} and 3^{(N−1)/2} with 0
violations.

- **Order-d form (codex_2).** W_{k,d} = (−1)^k s_k^{d−1} r_{k+1} ≡ a s_k^d (mod N). So χ(W_{k,d}) = χ(a) for every
  local character χ of order dividing d. For example 2^{N−1} = (2^{(N−1)/4})² is a known 4th power on N ≡ 1 (mod 4).
  Evaluating a higher character at an unknown prime is a separate computational problem, so this is a certificate,
  not a leak.
- **Argument-specific identities.** For a = ⌊√N⌋ and D = N − r²: r_{k+1}² + D s_k² = N((s_k − r h_k)² + D h_k²).
  This restates (−D | p) = 1 (cx1 §1.2).
- **The identities and the bound are all that is proved.** We do *not* claim that these are the only ways a
  computable path statistic can depend on p (withdrawn after codex_2).

**Meaning (HEURISTIC/structural).**
- F_a(s, h) = s(as − Nh) has discriminant N², a square, so the ordinary genus theory of non-square discriminants
  does not apply. The values of the path satisfy the one multiplicative constraint above.
- Compare the principal cycle of √(kN): A² − kNB² = (−1)^i Q_i with discriminant 4kN gives ((−1)^i Q_i | p) = 1
  at every step, i.e. (Q_i | p) = (−1|p)^i. This is the genus source of the eps/leakmap n = 32 hits. (v1 dropped the
  sign; codex_2's example: N = 103·107, Q₁ = 205, (205 | 103) = −1.)
- For a known square a = x², the path yields congruences (xs_k)² ≡ (−1)^k s_k r_{k+1} (mod N) whose right-hand
  sides lie in [N/(a_{k+1}+2), N). They are small only after a large partial quotient.
  - For a *uniform random* a, a partial quotient ≥ N^c has probability ≍ N^{−c} log N (Aistleitner–Borda–Hauke).
    For structured arguments this tail estimate does not transfer: for a = ⌊√N⌋ the k = 0 product is a ≈ √N.
  - Magnitude alone does not determine the relation yield (square factors, the separate factors s_k and r_{k+1},
    correlations). We make no complexity comparison with CFRAC.
- This is consistent with prediction P1 and with the null battery at the square arguments. It is not a proof that
  they cannot leak.

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
- **(ii) Both ends (KNOWN in substance: equidistribution on the modular hyperbola, Weil/Estermann; e.g. Shparlinski,
  "Modular hyperbolas", Japan. J. Math. 7 (2012)).** Weil gives |S(h,k;N)| ≤ 4√(N·gcd(h,k,N)). The gcd terms need p | h and p | k (so max(|h|,|k|) ≥ p) and so drop out
  below the cutoff H = c√N < p. Hence the rectangle discrepancy of (a/N, ā/N) is ≪ N^{−1/2} log² N, and
  d_TV(joint T-prefixes, independent model) ≪ T⁴N^{−1/2} log² N, which is useful for T ≤ N^{1/8−η}.
- **(iii) Obstruction (cx1 §2.3).** For T ≥ N^{1/4+η}, the joint law is at TV distance → 1 from the independent
  model: N atoms cannot fill T⁴ cells. No endpoint-equidistribution argument can give a full-path floor. A
  full-path floor needs the dynamics: the transfer operator of the Gauss map with congruence twists, at a *fixed*
  modulus.
  - Baladi–Vallée (2005) and Bettin–Drappeau (2022) average over moduli and are not that.
- **Literature at fixed modulus** (cx1 §2.4; checked by cx against the papers; our literature critic to re-check):
  - Porter, "On a theorem of Heilbronn", Mathematika 22 (1975) 20–28: the mean length has a uniform O(N^{−1/6+ε}) error.
  - Bykovskii–Frolenkov, Mat. Sb. 208:5 (2017), improve this to a normalised error O(N^{−6/29+ε}) for balanced N.
    The exponents were cx1's reading; the literature critic could not access the paper (UNVERIFIED).
  - Aistleitner–Borda–Hauke (arXiv 2210.14095): tails of Σ a_i and of the maxima at fixed N.
  - David–Kim–Mor–Shapira, "On the rate of convergence of continued fraction statistics of random rationals"
    (arXiv 2401.15586), Thm 1.2: polynomial large deviations for the length and for fixed-word frequencies, with
    exponents that are not quantified.
  - None of these compares the full law of an arbitrary class-E functional.
- **Consequence for targets (cx1 (29)).** A reference law within ρ_N of the actual law gives
  |E[Y g] − E[Y B(C)]| ≤ E|Y|ρ_N.
  - With a fresh uniform a, every bounded predictor reading only P_T(a/N) for T = N^{1/2−δ} has target advantage
    ≤ O(N^{−2δ}) beyond what the scale T (and any residue labels used) already gives.
  - The reference depends on N only through T = T(N) and those labels.

### 2.2 Theorem F-SQ (structured archimedean argument ⌊√N⌋). PROVED for ε; the tb step is abbreviated (no prior formulation found by us; a new application of the standard type-II method used for {√p} mod 1, cf. Balog 1983, Harman 1983)
**Setting.**
- p₁, p₂ are independent and uniform over the primes in I = [3P/4, P), with P = 2^{n/2}. This is the D_n
  sampler; conditioning on N ≡ 1 mod 4 and on p₁ ≠ p₂ changes only constants.
- A target class c is defined by product weights α(p₁)β(p₂): ε = ±1 means p₁ ≡ p₂ ≡ u (mod 4).
- Write y = √N/m₀ mod 1, r = ⌊√N⌋, δ = {√N}, t = ⌊(N − r²)/r⌋, s′ = N − r² − tr.

**Statement.** Let m₀ ≤ N^{1/6} and 2 ≤ Q ≤ √N. Let G be any [−1,1]-valued function of (r mod m₀, t,
P_Q(s′/r)), where P_Q is the stopped word without the overshooting digit.
- This covers the tail word and the convergent matrices mod m₀ of the path of ⌊√N⌋/N up to full denominator
  ≤ Q(r+t). Those are reconstructed from A mod m₀, with A = r + t.
- It does **not** cover arbitrary functions of the huge first digit A itself. codex_2's collision:
  N = 193·229 and N = 197·233 share t and the stopped word but have first digits 210 and 214. Remainder-residue statistics also need N mod m₀. They are
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
4. **Double large sieve** (Bombieri–Iwaniec 1986, Lemma 2.4, as located by codex_2):
   |Σ_m Σ_n α_m β_n e(x_m y_n)|² ≪ (1 + XY)·Σ_{|x_m−x_{m′}|≤1/Y} |α_mα_{m′}|·Σ_{|y_n−y_{n′}|≤1/X} |β_nβ_{n′}|.
   - Use x_p = √p − c₀ (centred, |x_p| ≤ X ≍ √P) and y_q = (h/m₀)√q; the factor e(c₀ y_q) goes into β_q.
     Then (1+XY) ≍ hP/m₀.
   - Spacing: |x_p − x_{p′}| ≥ |p−p′|/(2√P). Within 1/Y there are ≤ 1 + 4m₀/h primes, and within 1/X there are
     ≤ 1 + O(m₀/h).
   - Hence |Σ_{p,q} α_p β_q e(h√(pq)/m₀)| ≪ (hP/m₀)^{1/2} π(I)(1 + m₀/h).
   - Divide by π(I)² ≍ (P/log P)², with √P = N^{1/4}: |E_c e(h√N/m₀)| ≪ (log N) N^{−1/4}(h^{1/2}m₀^{−1/2} + m₀^{1/2}h^{−1/2}).
5. **Summing over h ≤ H = N^{1/6}:** D_c ≪ N^{−1/6} log N + m₀^{1/2}N^{−1/4} log N.
6. **Target tb.** For unsorted primes the region is |log p₁ − log p₂| ≤ log(8/7), which is not a product.
   - Split G = G⁺ − G⁻.
   - Sandwich the indicator between Beurling–Selberg majorant and minorant trigonometric polynomials of degree T
     in log p₁ − log p₂ (period > 2 log(4/3)). Each term is a product weight, and the coefficient ℓ¹ norm is
     O(log T).
   - The discrete error is the prime-pair mass of a boundary strip of width η = 1/T. By Brun–Titchmarsh it is
     ≤ O(π(I)(Pη + 1)) pairs per boundary, i.e. a normalized mass O(log P·(η + 1/P)) (codex_2).
   - With T = P^{1/6} this gives D_tb ≪ N^{−1/12+o(1)}, and hence the bound m₀Q²N^{−1/12+o(1)}.
   - The corresponding coverage is full denominators ≈ N^{13/24−δ}. ∎ (Detailed write-up of this step is
     abbreviated; the ε case is the main statement.)

**Numerical check (proof critic).** An exact computation over all prime pairs at n = 40, m₀ = 1 gives
|E_c e(h√N)| = 7.2·10⁻³, 1.7·10⁻³ and 3.6·10⁻⁴ for h = 1, 16, 1024, against the scale N^{−1/4} = 2.0·10⁻³. The
observed decay is ≈ N^{−1/4}h^{−1/2}, and the two mod-4 classes agree to within noise.

**Remarks.**
- F-SQ covers zone (a), the path of ⌊√N⌋/N up to full denominators ≈ N^{7/12−δ}.
  - The CF of ⌊√N⌋/N starts with the digit r + t ≈ √N.
  - The CF of s′/r then has its own forward half up to inner denominators ≈ N^{1/4}, i.e. full denominators
    ≈ N^{3/4}.
  - Zone (b), from N^{7/12} to N^{3/4}, is archimedean but beyond the double large sieve's N^{−1/6} resolution.
- The same method should handle ⌊√(kN)⌋, and ⌊N^{2/3}⌋ via phases h(pq)^{1/3}. For the latter the spacing
  suggests D ≪ N^{−1/6+o(1)}. This is a SKETCH, not proved.
- **Normalisation details (codex_2).**
  - Normalise by M_u(M_u − 1), where M_u ≈ π(I)/2 counts the primes ≡ u (mod 4). Removing the diagonal p₁ = p₂
    costs O(log P/P).
  - Erdős–Turán with H = P^{2/3} and m₀ ≤ P^{1/3} gives D_c ≪ X^{−1/6} log X (X = P² ≍ N).
  - The reference is Bombieri–Iwaniec, *On the order of ζ(1/2+it)*, Ann. Sc. Norm. Sup. Pisa (4) 13 (1986),
    449–472, Lemma 2.4 (pp. 452–453).
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

These are the Fourier–Dedekind sums of Beck–Diaz–Robins (J. Number Theory 96 (2002) 1–21). Their σ_t has
(λ^{a_k} − 1) in the denominator, so F_d = (−1)^d σ_n.

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
    3-dimensional polytope, so it is polynomial-time by Barvinok (Math. Oper. Res. 19 (1994) 769–779), with
    polynomial weights via differential operators (De Loera–Hemmecke–Köppe–Weismantel, Math. Oper. Res. 31 (2006);
    lemma number UNVERIFIED). It was not implemented here.
  - **KNOWN.** Polynomial-time computation of Fourier–Dedekind sums is announced by Xin–Xu, arXiv:2303.01185,
    which answers a Beck–Robins open problem via a modular Barvinok decomposition.

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

### 3.3 Rademacher-type congruences. PROVED (cx1 Thm 4; verified by us on 300 random moduli, all M including even and 3 | M)
- Shift n = 0: KNOWN. For d = 2 this is Rademacher's congruence (Rademacher–Grosswald, *Dedekind Sums*, 1972)
  together with F₂(0) = (M−1)/(4M) − s(ab⁻¹, M) (Beck–Diaz–Robins §4). For d = 3 it follows from FD3-0. Compare
  Chae–Jun–Lee (arXiv:1406.3563, Thms 1.2/5.1), who give such congruences for unshifted generalised Todd/Bernoulli
  Dedekind sums, including Zagier's.
- The shift-dependent (polynomial in n) form: no prior formulation found by us (literature critic).
- 12M·F₂(n; a, b; M) ≡ (−a² − b² − 3ab + 6n(a+b) − 6n²)/(ab) (mod M).
- 24M·F₃(n; a, b, c; M) ≡ (4n³ − 6An² + 2(A²+B)n − AB)/(abc) (mod M), with A = Σa and B = Σ_{i<j} a_i a_j.
- Our independent fit of the d = 3 hyperplane sum, C/M² ≡ (1/24)Σ_{i≠j}a_iā_j − (n/12)Σ_i a_i∏_{j≠i}ā_j −
  (n/4)Σā_i + (n²/4)Σ_{i<j}ā_iā_j − (n³/6)ā₁ā₂ā₃, is the same statement (`test_aggregate_congruences_d2_d3_shifted`, which uses gcd(M, 6) = 1; cx1 Thm 4 was checked separately on all M, including even M and 3 | M).
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
  DS1) and the frequency split s(a,N) = s(a,p)/q + s(a,q)/p + U′ are different, so "the local part" is not
  canonical (cx1 §3.1).
  - The unweighted modular decoder does not automatically decode the weighted (frequency) oracle.
  - That oracle can still factor. At a = 1, X = 12N(s(1,p)/q + s(1,q)/p) = H² − 2N − 3H + 4 with H = p + q
    (codex_2). v1's "only the first factors" is withdrawn.

### 3.5 Toy test of the shifted d = 3 sums (pre-registered §6; exhaustive populations)

**Correction (codex_2 finding 9).** The registered toy sums all collapse to classical Dedekind sums:
- F₃(1; a, b, 1; N) = F₃(0; a, b, 1; N) − F₂(0; a, b; N), by telescoping F₃(n) − F₃(n+c) = F₂(n);
- this equals −(N−1)/(8N) + ½[s(ab⁻¹, N) − s(a, N) − s(b, N)].

Likewise F₂(1; a, 1; N) = −(N−1)/(4N) − s(a, N), so the main-battery `dr1_*` features are classical Dedekind
sums. (We checked the telescoping identity exactly; codex_2 checked 356 cases.)

So the registered toy test below is a test of *classical Dedekind combinations* at the paired arguments, not of
three-dimensional content.
- **Populations.** Exhaustive: n = 20 has 324 semiprimes (10-bit primes, N ≡ 1 mod 4) and n = 22 has 1 089.
- **Cost.** The exact O(N) evaluator takes 0.36 s per instance (4 sums) at n = 20 and 1.5 s at n = 22. The whole
  run took 9.8 min on 3 workers.
- **Results.** 16 statistics × 2 targets per size; one statistic is degenerate in each cell.
  - max |z| = 1.64 (ε) and 1.80 (tb) at n = 20.
  - max |z| = 1.66 (ε) and 1.86 (tb) at n = 22.
  - Mean z² is 0.68–1.07. There is no hit at the registered |z| ≥ 4.5.
- **Status.** A sanity check only, as registered. The MDE is |ρ| ≈ 0.35 at n = 20 and 0.19 at n = 22, primes are
  heavily reused, and the populations are tiny.
- **Nonunits (rule 3).** Four argument nonunits were zero-coded, all ACCIDENTAL:
  - n = 20: `cbrt2` at 674413 = 769·877 (⌊N^{2/3}⌋ = 10·769); `hash2` at 727393 = 773·941.
  - n = 22: `cbrt2` at 3019553 and 3380477.
  - The null expectation Σ(1/p + 1/q) per argument and instance gives ≈ 1.4 (n = 20) and ≈ 2.4 (n = 22) per
    argument.
- **Inference.** Exhaustive-population permutation scores are descriptive. Prime reuse makes their sampling
  interpretation unclear (codex_2).

**Post-hoc rerun (NOT registered) with shift n = (N−1)/2** (`toy3.py --shift half`). This shift does not telescope
in a few steps: it would need j = n·a_i⁻¹ ≡ −a_i⁻¹/2 (mod N) steps along any single argument.
- n = 20 (324 instances): max |z| = 1.63 (ε) and 1.76 (tb).
- n = 22 (1 089 instances): max |z| = 1.71 (ε) and 2.64 (tb); mean z² = 0.69 and 1.46.
- No |z| ≥ 4.5. This is null at toy power. It is the only test of non-collapsing shifted d = 3 content, and it is
  very weak.

---

## 4. Q4: reduced-form coefficients (genusleak's open class)

### 4.1 What is proved
- **Reduction is class H** (KNOWN; cx1 §4.1). Indefinite reduction of a form of discriminant 4kN is the CF of a
  quadratic irrational of that discriminant, with unbounded state (P_i, Q_i). The principal cycle satisfies
  A² − kNB² = (−1)^i Q_i, which gives the genus constraints.
- **Positive definite (−4kN).** Gauss reduction is reduction in the upper half-plane, not a rational a/N path.
- **Retention of N** (cx1 §4.2; cf. leakmap REP-C, genusleak §6.2). Coefficients together with the discriminant
  retain N, so unrestricted post-processing can compute any function of N. This is not a polynomial-time
  completeness theorem. A floor needs a restricted statistic grammar, e.g. bounded-state additive statistics along a
  polynomial-length walk.
- **Split vs non-split (structural, HEURISTIC; §1.3).** The form whose values the path visits has discriminant N²
  (split) for class E and 4kN (non-split) for class H.
  - In class H the visited values are |Q_i| ≤ 2√(kN), each with a quadratic constraint ((−1)^i Q_i | p) = 1.
  - In class E the visited values are ≥ N/(a_{k+1}+2), with the constraint χ(W_{k,d}) = χ(a).
  - Neither statement exhausts non-multiplicative content. There is no classification theorem "class E plus Q";
    genusleak's open class remains open.

### 4.2 Random classes: Heegner equidistribution. KNOWN, weak form only (v1's "floor" claim WITHDRAWN as FATAL by codex_2)
References: Duke, Invent. Math. 92 (1988) 73–90; Harcos–Michel, Invent. Math. 163 (2006) 581–655, Thm 2, with
error ε_g([G_K:G]D^{−δ}) for continuous g. Indicator statistics need smoothing, which degrades the rate. Michel,
Ann. Math. 160 (2004) 185–236, supplies the subconvexity input. Level-m₀ coefficient residues: no citable theorem
located (UNVERIFIED).
- A reduced form (a, b, c) of negative discriminant D corresponds to the Heegner point τ = (−b + √D)/(2a) in
  SL₂(Z)\H.
- **The theorems.**
  - Duke (Invent. Math. 92 (1988), Thm 1(i)): for fundamental D < 0, these points equidistribute on fixed convex
    regions with piecewise-smooth boundary, with ineffective constants.
  - Harcos–Michel (Invent. Math. 163 (2006), Thm 2, (1.4)): the same for cosets of subgroups of index
    o(|D|^{1/23042}). That includes genus cosets (the square subgroup), against fixed continuous tests.
  - Our D = −4N (N ≡ 1 mod 4) is fundamental. Non-maximal orders need the extra local estimates of
    Einsiedler–Lindenstrauss–Michel–Venkatesh §2.3. Level structure needs separate statements (Harcos–Michel
    Thm 6, §6.4).
- **What this does and does not give.**
  - It gives **weak convergence**: for any *fixed* regular test set, the genus-coset proportions agree up to o(1).
  - It does **not** give total-variation closeness. At fixed D distinct genus cosets have disjoint supports, and the
    coefficients (a, b) determine the class.
  - Coefficient residues can carry genus characters. For example χ₋₄([a,b,c]) = (−1 | a) for odd represented a.
    For our walk this character is public: χ₋₄(f_ℓ^j) = (−1|ℓ)^j.
  - **So there is no floor for arbitrary bounded coefficient statistics, even for random classes.** v1's claim is
    withdrawn.
- **The cusp remark.** The genus leaks of eps/leakmap come from the *positive*-discriminant principal cycle of
  √(kN) (closed geodesics), not from these Heegner points. "Small leading coefficient ↔ high in the cusp" is only a
  qualitative geometric analogy (for definite forms, Im τ = √|D|/(2a)). Squareness and smoothness are arithmetic
  events, and fixed-region equidistribution says nothing about shrinking cusp regions.
- **Not covered.** Short orbits {f^j : j ≤ poly(n)} are OPEN. Moreover, for j ≤ log(√(|D|/3))/log ℓ the reduced
  form f^j is just (ℓ^j, ·, ·), a public object. Our battery's FW statistics (j ≤ 32) are therefore largely
  uninformative (§5).

## 5. Q5: battery (pre-registered; `notes/classE_prereg.md`)
**Setup.**
- 2·10⁵ instances of D_n at each n ∈ {48, 64, 128}.
- Seeds: main = 25.52·10¹² + n·10⁹.
- 130 statistics: 117 `menu`, 5 `menuG`, 5 NC, 3 PC.
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
- V1, reported both ways:
  - ε passes, with the same single flip pattern as newfamily, hence 21.2.
  - **Registered rule:** tb FAILS at n = 48 (16.10 vs 22.36 ± 3) and n = 64 (16.75), and passes at n = 128.
  - **Amended rule (D1):** with the leading-order expectation 2ρ√(ΣS_yy) = 18.18 / 18.21 / 18.20, tb passes all
    three cells.
    - The amendment is a mathematical correction: the top bits of N explain part of tb's variance (codex_2 confirms
      the first-order formula).
    - It changes a registered validity rule post hoc, so the tb cells at n = 48 and 64 are valid only under the
      amended protocol.
- **MDE.** For ε, ρ_MDE = 0.0153 for linear residual correlation. For tb the effective size is feature-specific
  (S_xx S_yy / V). It ranges from ≈ 4ΣS_yy = 1.32·10⁵ (MDE 0.0188, the planted-flip calibration) to ≈ 2·10⁵ for
  features uncorrelated with size (e.g. `isqrt2_ms_mod3`: 1.98·10⁵, MDE 0.0154; codex_2).

**Conclusion.**
- **ε:** "refuted for this menu" by the registered rule at n ∈ {48, 64, 128}, with 80% power at |ρ| = 0.0153.
- **tb:** the same, at n = 128 under the registered rule, and at n = 48 and 64 under the amended V1 rule (D1).
- The defensible plain statement is **no registered threshold hit**. This is a power statement for *linear*
  association beyond the stated coarse floors, not an exclusion. Categorical codes are tested linearly only, and the
  top-12-bit floor does not remove all archimedean information. A within-stratum normalized-size check by codex_2
  gives |z| ≤ 1.12.

**Uninformative or redundant features (found after the run; codex_2 finding 15).**
- **FW walk.** For j ≤ log(√(|D|/3))/log ℓ the reduced form f^j is (ℓ^j, ·, ·), a public function of ℓ and j.
  With ℓ = 3 that is j ≤ 15 at n = 48 and j ≤ 20 at n = 64; at n = 128 it covers all 32 steps.
  - So `fw_*` are informative only on the later steps at n = 48 and 64, and not at all at n = 128.
  - `fw_minbits` equals the bit length of the starting prime in 199 998 of 200 000 records at n = 48 and in all
    records at n ≥ 64.
  - `fw_sq` counts the even j.
  - The FW group is therefore not a test of reduced-form walks. Q4's walk statistics remain untested at scale.
- **`div3_*`.** Bounded-modulus tail after one large quotient (Prop LIN): public residue functions plus size.
- **`*_dr1_*`.** Classical Dedekind sums (§3.5). The genuinely shifted d = 2 features are `fd2_*`
  (F₂(1; a, b) with b ≠ 1).
- P1–P4 were all met:
  - P2: no menuG hit, even at n = 48. The genus statistics of √(3N) at 256 steps are too weak at this size.
  - P4: 0 nonunit arguments and 0 Euclid-path gcd hits at every size.
    - The registered rationale "Σ(1/p+1/q) per remainder" was wrong: it predicts ≈ 8 events at n = 48.
    - Prop VIS gives the right rate for random-like arguments: ≈ 2(q−p)/N per argument, i.e. ≈ 0.02 events in the
      whole n = 48 run. The observed 0 agrees.
    - The only splits are 3 ambiguous forms on the n = 48 form walks, f^j with ord f | 2j ≤ 64, a small-order
      (family M) event. They are STRUCTURAL, with rate 1.5·10⁻⁵, and there are none at n ≥ 64.

**Post hoc (`posthoc.py`; not registered as tests).**
- tb under finer floors: max |z| 3.35 with top 16 bits (n = 48) and ≤ 2.74 at n ≥ 64.
  - The joint (top 8 bits, N mod 16) floor (3.06) is not a refinement of top-12; it is a different coarse floor.
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
**Blind round (classE_1, xhigh).** The report is `codex/classE_1_report.md` (914 lines). It reached the same
picture independently, and several of its results are used above with credit:
- Thm 1: the exact VIS iff with the cylinder inequality, and the exact random-a probability.
- Thm 2: one-end TV. (22)–(23): Weil joint discrepancy. §2.3: the TV → 1 obstruction for two fine ends.
- Thm 3: the conditional theorem.
- Thm 4: the d = 2, 3 congruences, with proof via the shift-difference identities (38)–(39).
- The DSF decoder (41) and its exceptions (42); the Barvinok reduction (44).
- The counterexamples to over-strong prefix-agreement claims.

Its top-ranked open mechanism is "partial local-value separation in two- and three-factor sums", for which the
decoder (41) is a precise target. Next comes junction/reverse data for exponent residues.

**Adversarial round (classE_2, xhigh): MAJOR REVISION.** The report is `codex/classE_2_report.md`, with audit
scripts in `codex/scratch/classE_2/`. All findings were applied in v2:
- **FATAL #6.** The Heegner "floor for every bounded shape statistic" is false (weak vs TV convergence; disjoint
  coset supports; χ₋₄ visible mod 4). Withdrawn (§4.2).
- **MAJOR.**
  - #1, #2: the "only two exact ways" exhaustiveness and the CFRAC domination are withdrawn. The order-d
    certificate W_{k,d} is added. The sign ((−1)^iQ_i|p) = 1 is restored.
  - #3: F-SQ does not control the first digit A. The statement is restricted and the collision example added.
  - #5: the tb step is repaired (Beurling–Selberg plus Brun–Titchmarsh strip).
  - #9: **the registered toy d = 3 sums, and `dr1_*`, collapse to classical Dedekind sums by telescoping.**
    Relabelled, plus a post-hoc non-collapsing rerun.
  - #10: the frequency local part can also factor (a = 1).
  - #11: toy nonunits disclosed; n = 22 completed.
  - #12, #13: V1 is reported under both the registered and the amended rule, and the MDE is feature-specific.
  - #15: FW, `div3` and `dr1` features are uninformative.
  - #16: the verdict is narrowed; the N^{2/3} extension is now a sketch; the Coppersmith gcd caveat is added; FD3-0
    is credited to cx1 (37) as well.
- **MINOR #4, #14, #17.** The double-large-sieve reference was located (Bombieri–Iwaniec 1986, Lemma 2.4), the
  normalisation was written out, the floor wording fixed, the counts corrected (117 menu), and the P4 rationale
  corrected via VIS.
- **Verified correct by codex_2:**
  - the VIS cylinder iff and the √N gcd obstruction;
  - the CF identities;
  - the SPLIT identity and bound;
  - the ε F-SQ analytic core (restricted statement);
  - FD3-0, Thm 4 and the DSF decoder (240 + 24 checks);
  - the evaluators;
  - our 11 tests;
  - the battery maxima, recomputed from raw data (all six match), and 9/9 operational vectors recomputed from N.

## 8. Opus critics
**Proof critic: ACCEPT WITH FIXES.** All fixes are applied.
- **MAJOR.** F-SQ's coverage was misdescribed. There are three zones: (a) proved; (b) archimedean, open; (c)
  Kloosterman-type, open.
- **MINOR.**
  - Remainder residues need N mod m₀ (via product weights).
  - The tb step is a sketch: Beurling–Selberg sandwich plus a short-interval prime bound.
  - A numerical re-check of the F-SQ Weyl sums was added: exact over all prime pairs at n = 40, decay
    ≈ N^{−1/4}h^{−1/2}, and the classes agree.
  - t ∈ {0, 1} in the unit case. The agreement wording is fixed.
  - The "Meaning" paragraphs are relabelled heuristic.
  - J₂ is Φ only at n = 0.
  - Attribution to cx1 §3.1 is fixed.
- **NIT.** Residue-label factor m; folded reference; the gcd condition p | h, k; the piece count; "p/q is never a
  convergent".
- **Verified by the critic:**
  - VIS and cx1 Thm 1, exhaustively, on 2 936 960 cases, including the exact count 2(q−p−1) per N;
  - SPLIT, on 32 725 products;
  - SQ-CF, on 1 200 instances at n ≤ 128;
  - §2.1 (i)–(iii);
  - the F-SQ steps;
  - all of §3 numerically (FD formula, floor sums, FD3-0, (39), (44), Thm 4 on all M, the DSF decoder on 200
    semiprimes, (42));
  - no factor leakage in the operational code.

**Literature critic: OK WITH FIXES.** All fixes are applied.
- §4.2 references are corrected:
  - Duke is for fundamental D, which our D = −4N is.
  - Harcos–Michel Thm 2 gives genus cosets with a test-dependent rate.
  - Michel 2004 is the subconvexity input only.
  - The level-m₀ claim is UNVERIFIED.
- n = 0 congruences relabelled KNOWN (Rademacher; Beck–Diaz–Robins; Chae–Jun–Lee). The shift-dependent form has
  no prior formulation found.
- Poly-time Fourier–Dedekind evaluation is KNOWN (Xin–Xu 2303.01185).
- §2.1(ii) relabelled KNOWN (modular hyperbola; Estermann's bound for composite moduli).
- Bombieri–Iwaniec Lemma 2.4 confirmed; Brun–Titchmarsh suffices for tb.
- F-SQ is a new application of the standard type-II method (Balog, Harman on {√p} mod 1). No prior formulation was
  found for F-SQ, Prop VIS or SPLIT.
- Full citations were supplied:
  - Coppersmith EUROCRYPT '96 / J. Cryptology 10 (1997);
  - Howgrave-Graham CaLC 2001;
  - Barkan 1977; Hickerson, Crelle 290 (1977); Knuth, Acta Arith. 33 (1977);
  - Zagier, Math. Ann. 202 (1973); Mordell 1951; Pommersheim, Math. Ann. 295 (1993), as context for FD3-0;
  - Baladi–Vallée, JNT 110 (2005); Bettin–Drappeau, PLMS 2022.
  - Page numbers given from memory are UNVERIFIED.
- The critic judged "Q-leaks as cusp excursions" a fair *description*. codex_2 judged it an analogy only; we keep
  the weaker wording.

## 9. Deviations from the pre-registration
- **D1.** The registered V1 expectation ρ√n assumes within-stratum target variance ¼. For tb under the top-12-bit
  floor, the size of N explains part of the variance, so the planted z is 2ρ√(Σ_h S_yy,h) (`posthoc.py`; derivation
  in its docstring).
  - By the registered formula, V1 fails for tb at n = 48 (16.1 vs 22.4 ± 3) and n = 64 (16.8).
  - By the corrected formula (18.2) it passes at every size (16.1, 16.8, 18.3).
  - The tb MDE is reported with the effective n = 4ΣS_yy.
- **D2.** `harness.py` was changed to stream records to disk before the main run (commit after the prereg; no
  change to statistics or analysis).
- **D3.** `toy3.py` gained a `--shift half` option for the post-hoc non-collapsing rerun
  (`experiments/raw/classE_toy3_half_n{20,22}.jsonl`). The registered default behaviour is unchanged.
- The registered toy test and the `dr1_*` features turned out to be classical Dedekind combinations (§3.5). This is
  a design error, disclosed, not a deviation in execution.
- The post-hoc finer floors (top 16 bits; top 8 bits × N mod 16) were run for tb although there was no hit
  (registered only as an audit step).

## 10. Reproduction (project root)
```
.venv/bin/python -m pytest -q tests/test_classE.py
PYTHONPATH=. .venv/bin/python -m src.classE.harness --n 48 64 128 --range main --count 200000 --workers 3
PYTHONPATH=. .venv/bin/python -m src.classE.posthoc 48 64 128
PYTHONPATH=. .venv/bin/python -m src.classE.toy3 --n 20 22 --workers 3
PYTHONPATH=. .venv/bin/python -m src.classE.toy3 --n 20 22 --workers 3 --shift half   # post hoc
```
