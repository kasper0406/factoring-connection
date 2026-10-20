# floorhard: hardness side of the Theorem L_j floor-sum oracle (phase 17)

Agent `floorhard`, 2026-10-20.

**Files.**
- Pre-registration: `notes/floorhard_prereg.md`, commit `6b3a37e`, made before any code existed.
- Code: `src/floorhard/{core,run_experiments}.py`. Tests: `tests/test_floorhard.py` (9 pass).
- Raw data: `experiments/raw/floorhard_e{1,2,3}.jsonl`. Summary: `experiments/analysis/floorhard_summary.json`.
- GPT-6.1-sol: blind round `codex/floorhard_1.{md,out}`; adversarial round `codex/floorhard_2.{md,out}` (§9).
- Opus critics (proof, literature): §9.

**Labels.** KNOWN (with reference), PROVED (proof here; all elementary; at most "no prior formulation
found by us"), HEURISTIC, CONJECTURED, OPEN, experimentally supported, UNVERIFIED.

**Consistency with `floorsum`.** floorsum (`notes/floorsum_prereg.md`) works on the algorithm side. Its
claims C1 (the convexity scale summed over blocks is N^{1/3} for every j) and C5 (Coppersmith-type
lattices already detect divisors in j = 2 blocks) agree with §5 here. Nothing below depends on
floorsum's experimental results.

## 0. Verdict (plain language)

**As Theorem L_j states it, the oracle hypothesis is far too strong to be plausible, for every j ≥ 2.**
The hypothesis is exact Σ⌊P(h)⌋ in polylog time for all rational P of degree ≤ j, with arbitrary
shifts and subintervals. It is #P-hard and NP-hard. It is hard even when restricted to the ladder's
own Taylor polynomials of N/(x+h), provided the oracle must work for every N, or for every prime N.

A fast uniform oracle would not merely factor in N^{1/(j+2)}. It would factor any m directly, with
O(log m) calls on intervals of length ≤ m. So the ladder is pointless for uniform oracles. It can
only matter for an algorithm that exploits the special arithmetic of the blocks of one given N.

**Where the hardness lives.** Every hard instance we can build sits in blocks well below √N, at
x ≈ N^{1/3}–N^{3/8}, at specially structured N. In these blocks the hyperbola locally reproduces a
quadratic-residue pattern h² mod b with many wraps. A balanced semiprime only needs the top blocks
x ∈ [√N/κ, √N]. There:
- the j = 2 arcs wrap only O(κ²) times, and Coppersmith–May already detects divisors in them;
- for j ≥ 3 the arcs are short relative to their arithmetic modulus, which is the lattice
  (Coppersmith) regime, not the NP regime.

**We found no reduction from any believed-hard problem other than factoring itself to the top-block
oracle.** Its only known consequence is the ladder. That amounts to "factoring with (1/(j+2))·log₂N
known high bits of p": beyond Coppersmith's n/4 for j ≥ 3, and approaching full factoring as j → ∞.

**Plausibility (HEURISTIC).** A top-block oracle with T_j = N^{o(1)} for some fixed j ≥ 4 would beat
the deterministic record N^{1/5}. We see no mechanism for it. We put it at about 2–5%; codex
estimates about 10% for a uniform H^{o(1)} oracle. The NP-hardness does **not** contribute to this
number. The hard families are themselves solvable in H^{o(1)} (§4.3), so hardness does not exclude
any oracle cost the ladder needs.

**Sharp barrier statement (§7).** For every fixed j ≥ 2:
- a uniform exact Taylor-class oracle in time poly(log N) gives P = NP and FP = #P;
- in time 2^{o(√log N)}, it refutes ETH;
- in deterministic time T(H), it factors every m deterministically in O(log² m · T(m)), which beats
  Theorem L_j itself whenever T(H) ≤ H^{2/(j+4)};
- for the top-block oracle (the only part a balanced-semiprime ladder needs), no implication beyond
  Theorem L_j is known: OPEN.

## 1. Basic identities

**Lemma W (window identity; PROVED).** Let f be an integer polynomial, D ≥ 1, and 0 ≤ u < v ≤ D. Then
#{h ∈ [h₁, h₂) : f(h) mod D ∈ [u, v)} = F(f − u) − F(f − v), where F(g) = Σ_{h₁≤h<h₂}⌊g(h)/D⌋.
- Proof: ⌊(f−u)/D⌋ − ⌊(f−v)/D⌋ is 1 exactly when f mod D ∈ [u, v), and 0 otherwise.
- For real P, the same argument gives #{h : {P(h)} ∈ [α, β)} = Σ(⌊P − α⌋ − ⌊P − β⌋).
- Consequences:
  - with constant shifts, a floor-sum oracle **is** a residue-distribution oracle;
  - Σ⌊f/D⌋ = (Σf − Σ(f mod D))/D, so floor sums and residue sums are interchangeable;
  - the window [0, 1) counts the solutions of f(h) ≡ 0 (mod D) in the interval.

## 2. Reductions from fast degree-2 floor sums

### 2.1 Arbitrary intervals: NP-, #P- and ETH-hardness

**Theorem SS (PROVED; parsimonious).**
- Given positive integers w₁, …, w_k and T, one can construct in polynomial time a, b and at most two
  intervals I such that #{x ∈ I : x² ≡ a (mod b)} = #{δ ∈ {0,1}^k : Σδ_i w_i = T}.
- The count is 2 (or 4) exact floor sums of x² − a and x² − a − 1 over D = b (Lemma W).
- The bit length is L = O(k log(kS)), where S = Σw_i.

*Proof.*
1. **Parameters.** Let M = 4S + 4. Choose k distinct primes p_i > max(kM, 4M) (deterministic search
   with a primality test; the primes in (P, 2P] suffice). Put b = ∏p_i.
2. **The θ_i.** Let u_i = round(p_i w_i/M); then 1 ≤ u_i < p_i/4. Put θ_i = (b/p_i)u_i, and
   s_i = θ_i mod p_i ≠ 0.
3. **The square roots.** Let a = (Σθ_i)² mod b, so a ≡ s_i² (mod p_i). Since each p_i is odd and
   s_i ≢ −s_i, the roots of x² ≡ a (mod b) are exactly x_ε = Σε_iθ_i mod b, ε ∈ {±1}^k, all distinct.
4. **Their positions.** x_ε/b ≡ Σε_iu_i/p_i ≡ Σε_iw_i/M + η (mod 1), with
   |η| ≤ Σ1/(2p_i) < 1/(2M).
5. **Separating the target.** Put ε_i = 2δ_i − 1. Then Σε_iw_i = 2Σδ_iw_i − S. All such values have the
   parity of S and lie in [−S, S], so they occupy distinct points of (−¼, ¼) spaced ≥ 2/M apart. The
   window x/b ∈ [(2T′−1)/(2M), (2T′+1)/(2M)) mod 1, with T′ = 2T − S, therefore contains exactly the
   roots with Σδ_iw_i = T.
6. **No boundary hits.** Every root lies at distance > 0 from the window ends. This gives the integer
   endpoints used in `core.subset_sum_to_qc`. ∎

**Corollaries.**
- **(i) #P- and NP-hardness (PROVED given KNOWN inputs).** Exact degree-2 floor sums with L-bit
  coefficients over intervals of length < 2^L are #P-hard. #SubsetSum is #P-complete: the textbook
  3SAT → SubsetSum reduction, with clause target 4 and slacks {1, 2}, is parsimonious (3 = 1 + 2,
  2 = 2, 1 = 1 uniquely), and #3SAT is #P-complete.
  - Hence polynomial time (in L) would give FP = #P and P = NP.
  - NP-hardness alone is KNOWN: Manders–Adleman, JCSS 16 (1978) 168–184. Bounded quadratic congruence
    "∃x < c, x² ≡ a (mod b)" is NP-complete, even with the factorisation of b given (Garey–Johnson
    [AN1]).
  - The codex blind round derived (3) from Manders–Adleman with two calls.
- **(ii) ETH (PROVED given ETH ⇒ Subset Sum needs 2^{Ω(k)} with weights ≤ 2^{O(k)}; folklore via
  sparsification).**
  - Here L = O(k²). So exact degree-2 floor sums are not computable in 2^{o(√L)} time unless ETH
    fails.
  - This is the strongest fine-grained statement we have. It excludes nothing between 2^{√L} and
    2^{εL}, i.e. nothing between N^{o(1)} and N^{ε} (§4.3).

**Prop BIS (deterministic factoring by bisection; PROVED; codex found it independently).**
- Let m be odd with ω(m) ≥ 2. Then x² ≡ 1 (mod m) has 2^{ω(m)} roots, and roots other than ±1 exist.
- Count the roots in [lo, hi) with two floor sums (Lemma W, f = x² − 1, D = m). Bisect [2, m − 1) to
  isolate a root r. Then 1 < gcd(r − 1, m) < m.
- This takes O(log m) calls on intervals of length ≤ m. Complete factorisation takes O(log² m) calls,
  after removing perfect powers.
- Square roots of any a mod m are located the same way.

### 2.2 Full periods: class numbers, residuosity, character sums

**Prop FP (PROVED; Dirichlet for prime m).** For odd squarefree m and gcd(a, m) = 1,

  Σ_{h=0}^{m−1} ⌊ah²/m⌋ = a(m−1)(2m−1)/6 − (m−1)/2 + Σ_{d | m, d ≡ 3 (mod 4)} (a/d)·h*(−d),

where h*(−3) = 1/3 and h*(−d) = h(−d) for d > 3 (−d is a fundamental discriminant).

*Proof.*
1. **Root counts.** The number of h mod m with ah² ≡ r is Σ_{d|m}(ā r/d), by CRT
   (∏_{p|m}(1 + (·/p)) expanded).
2. **Residue sum.** So Σ_h (ah² mod m) = m(m−1)/2 + Σ_{d|m, d>1}(a/d)(m/d)Σ_{s<d}s(s/d).
3. **Character sums.** Σ_{s<d}s(s/d) is 0 for even characters (d ≡ 1 mod 4). It is −d·h*(−d) for
   odd ones (Dirichlet's class number formula).
4. Divide by m and subtract from aΣh². ∎

- Checked exactly: 6 520 (m, a) pairs, 0 failures (E1a). Codex states the same identity as its (6).
- For prime m it is the barrier Q of `notes/analytic_results.md`.

**Corollaries.**
- **(FP1) Class numbers.** One full-period call gives h(−p) for prime p ≡ 3 (mod 4).
  - Shifted full-period calls give the incomplete sums Σ_{t≤X}(t/p) (codex (4)), hence h(−4p) via
    Σ_{t<p/4}.
  - Best known: L[½] under GRH (Hafner–McCurley; Biasse–Erukulangara). No polylog algorithm is known.
  - That h(−p) needs super-polylog time is CONJECTURED; nothing is proved.
- **(FP2) Factoring from full periods (PROVED; randomized; experimentally supported).** Let m = pq,
  3 ∤ m.
  - Put M = m if m ≡ 3 (mod 4). Then the class part is ±h(−q) + (a/m)h(−m), where q is the factor
    ≡ 3 (mod 4).
  - Otherwise put M = 3m. Then the class part is (a/3)/3 ± (h(−p) ± h(−q)) + (a/3m)h(−3m).
  - Restrict to a with (a/M) = 1 (and (a/3) = 1). The unknown sign takes both values with probability
    ½ each, and the mean of the two observed values is h(−M).
  - Then h(−4M) = h(−M)(2 − (−M/2)), and ambiguous forms (Shanks; Schoof 1982) split M.
  - E1b: h(−M) decoded correctly in 30/30, median 3 oracle calls (max 6); m factored in 29/30. The one
    miss was M = 3m, where the ambiguous-form step returned only the factor 3 in 8 tries.
  - Codex left the composite case UNVERIFIED for unbounded ω(m). The randomized decoder here covers
    ω(m) = 2; general ω needs ≈ 2^{ω} sign patterns.
- **(FP3) Residuosity.** The window [0, 1) over a full period counts the square roots of c mod m,
  ∏(1 + (c/p)). This decides quadratic residuosity mod a composite (the Goldwasser–Micali assumption).
  - Windows [0, Y) give Σ_{d|m}(ā/d)·(incomplete Jacobi sums mod d).
  - So full-period quadratic floor sums select local quadratic characters (the U1 of
    `notes/nonsmooth_results.md`).

**Which ranges each reduction needs.** Each reduction needs the quadratic coefficient to wrap many
times over the interval, i.e. curvature K = (a/m)·|I|² ≫ 1:

| Consequence | Interval | Wraps K = a·\|I\|²/m |
|---|---|---|
| NP / #P (Thm SS) | length ≈ b/M ≈ b^{1−1/k} | ≈ b^{1−2/k} |
| factoring (BIS) | adaptive subintervals of [0, m) | up to m |
| h(−p), QR, incomplete sums (FP) | full period | m |

Short-interval *exact-congruence* counting with |I| ≤ √(b/a) is trivial (x² = a + by with y ≤ 1). The
regime |I| < b^{1/2} is the Coppersmith small-root range.

## 3. The ladder's restricted class

**Theorem TE (Taylor embedding; PROVED).**
- Fix j ≥ 2, A ≥ 1, b ≥ 2, c′ ≥ 1. Let t be an integer with t ≥ 2c′/b and t > 8Ac′³/b. Put N = Ab²t³
  and x = bt.
- Then the degree-j Taylor polynomial P_j of N/(x+h) at x satisfies, for 0 ≤ h < c′,

  P_j(h) = Abt² − Ath + Ah²/b + R(h),  |R(h)| < 1/(4b).

- Hence #{h ∈ I : Ah² ≡ r (mod b)} = #{h ∈ I : {P_j(h) + 1/(4b)} ∈ [r/b, r/b + 1/(2b))}, for every
  subinterval I ⊆ [0, c′). This is two oracle calls with rational shifts (Lemma W).
- Choosing t ≥ (const)·A c′⁴/b² also makes the ladder block [x, x + ⌊cx/B⌋) contain [x, x + c′).

*Proof.*
1. P_j(h) = Σ_{i≤j}(−1)^i N h^i/x^{i+1} = Σ_i (−1)^i A h^i t^{2−i}/b^{i−1}.
2. The terms i ≤ 1 are integers; the term i = 2 is Ah²/b.
3. The terms i ≥ 3 total at most (Ac′³/(b²t))·Σ_k(c′/(bt))^k ≤ 2Ac′³/(b²t) < 1/(4b).
4. The block condition is H = ⌊cx/B⌋ ≈ c·b^{1−2/(j+2)}t^{1−3/(j+2)}A^{−1/(j+2)}. This is
   ≥ c b^{1/2}t^{1/4}A^{−1/4} for j ≥ 2. ∎

**Variants.**
- **Prime N (PROVED, Las Vegas).** N = Ab²t³ + e with 0 ≤ e < bt²/(16c′) and t > 16Ac′³/b adds a
  perturbation with variation < 1/(8b).
  - The admissible e form an interval of length ≍ N^{2/3}·(const), which exceeds N^{0.525}. By
    Baker–Harman–Pintz (Proc. LMS 83 (2001)) it holds ≫ N^{0.525}/log N primes once t is large.
  - Random e plus AKS find one in expected O(log N) trials.
  - So an oracle that is correct only for prime N already gives NP ⊆ ZPP if polylog. Balanced-semiprime
    N would need balanced semiprimes in intervals of length N^{2/3}: HEURISTIC, not proved here.
- **Codex's one-call embedding (PROVED in codex/floorhard_1.out (7)–(11); re-derived here).**
  - Take x = 16m⁴ and N = (ax³ − bx²)/m + 1. Then gcd(N, x) = 1 and x² ≤ N ≤ x³.
  - Any (ah² + bh + c)/m floor sum on [0, m) equals one Taylor-class query minus an explicit integer.
- E3a: 60/60 exact window checks (j = 2, 3, 4). E3b: 20/20 with prime N. E3c: 40/40 composites m ≤ 20 000
  factored using only embedded Taylor-class window counts (≤ 32 calls, ≤ 2.2·log₂m).

**Corollary DOM (uniform oracles make the ladder pointless; PROVED).**
- Suppose an oracle evaluates the Taylor class for all N exactly, on intervals of length ≤ H, in
  deterministic time T(H) (polynomial in the bit length absorbed). Then, by Prop BIS through
  Theorem TE (A = 1, b = m, c′ = m, N ≈ m⁸), every m is completely factored in O(log² m·T(m)).
- Theorem L_j gives N^{1/(j+2)}·T(N^{1/2−1/(j+2)}).
- For T(H) = H^δ, the direct route is better iff δ < 2/(j+4). In particular:
  - polylog T: polylog factoring, and P = NP by §2.1;
  - H^{o(1)}: m^{o(1)} deterministic factoring;
  - H^{1/5−ε}: beats Harvey's deterministic N^{1/5}.
- So L_j is informative only for oracles that are **not** uniform over N: fast on the blocks of the
  given N, slow on the embedded structured blocks.

**Degree ≥ 3.** The general degree-≤j class contains degree 2. Theorem TE holds for every j ≥ 2 (the
i ≥ 3 terms are absorbed in R). Codex also handles exactly-degree-j numerators. So every hardness
statement above holds verbatim for all j ≥ 2. Degree adds nothing to the uniform barrier.

## 4. Where the hardness lives

### 4.1 Location of the hard blocks (PROVED)

- In Theorem TE the hard block sits at x = bt = (Nb/A)^{1/3}. Hard instances need c′ ≥ b^{1−o(1)}
  (§2.2), so t ≥ 8Ac′³/b ≥ b^{2−o(1)} and N ≥ b^{8−o(1)}. The minimal choice gives x = N^{3/8+o(1)};
  larger t moves x down towards N^{1/3}.
- In codex's embedding x = N^{4/11+o(1)} (x = 16m⁴, N ≈ 4096 a m^{11}).
- Both are polynomially below √N. Both need N to be close to a value with x² dividing it approximately:
  N/x² ≈ t ∈ Z and N/x³ ≈ A/b.

### 4.2 The top blocks (all a balanced-semiprime ladder needs)

For N = pq with p ∈ [√N/κ, √N], the ladder only scans x ∈ [√N/κ, √N]. There:

- **(T1) j = 2: bounded wraps (PROVED).** N/x³ ≤ κ²/x, and B = N^{1/4} ≥ √x gives H ≤ c√x. So the
  quadratic part N h²/x³ moves by ≤ c²κ² over the block.
  - ⌊P₂⌋ is therefore a linear floor plus a count of {ℓ(h)} under O(κ²) parabolic thresholds.
  - Any vertical-frame congruence encoding (N/x³ ≈ u/b) has b ≳ x/κ² ≥ H²/(c²κ²): the
    "trivial/Coppersmith" range |I| ≲ √b of §2.2.
- **(T2) j = 2: detection is already solved (KNOWN).** Coppersmith–Howgrave-Graham–May find a divisor
  d ≈ N^β in an interval of length N^{β²}. Since β² ≥ β − ¼ (equality only at β = ½), every j = 2
  ladder block H = xN^{−1/4} lies within that range (floorsum C5). The j = 2 oracle is not needed for
  the ladder's purpose.
- **(T3) j ≥ 3: short relative to the modulus (HEURISTIC obstruction).**
  - The quadratic coefficient is N/x³ ≈ 1/x (real), so any rational model u/b has b ≳ x/κ². The block
    covers a fraction H/b ≲ κ²cN^{−1/(j+2)} of one period of h² mod b.
  - CRT encodings (Theorem SS, Manders–Adleman) need a period fraction ≳ k·b^{−1/k} with k → ∞ primes.
    Here at most k ≈ (j+2)/2 primes fit, so there are only 2^{O(j)} roots: no NP-hardness.
  - The non-vertical (continued-fraction) frames give the same exponent balance up to constants at
    j = 2. We could not make them embed a CRT instance, and we have no proof that none exists.
- **(T4) What the top-block oracle is equivalent to (PROVED one way).**
  - Window counts of top blocks are lattice counts in thin strips N/d − β < m ≤ N/d − α. Equivalently
    they are Σ_{d∈block}[N mod d ∈ [αd, βd)], i.e. divisor counts of N − k summed over a window of k.
  - The width-0 limit is "does N have a divisor in [x, x + N^{1/2−1/(j+2)})". That is factoring with
    (1/(j+2))·log₂N known top bits of p, beyond the Coppersmith bound (¼) for j ≥ 3.
  - The converse (divisor detection ⇒ exact strip counts) is not known.

### 4.3 NP-hardness does not touch the ladder's regime (PROVED)

- The #SubsetSum family of Theorem SS has b = ∏p_i with p_i = 2^{O(k)} = b^{O(1/k)}. Trial division
  finds all p_i in b^{O(1/k)} = b^{o(1)} time. The 2^k CRT roots can then be enumerated in
  2^k = 2^{O(√log b)} = b^{o(1)} time.
- So the hard family is solvable in H^{o(1)}. Codex makes the same point for Manders–Adleman with known
  factorisation: 2^{ω(m)} = m^{o(1)}.
- Every fixed-j consequence of Theorem L_j needs only T_j = N^{ε}. **The NP/#P/ETH hardness is
  therefore irrelevant to the question that matters.** It refutes only the literal polylog hypothesis.

## 5. Degree ≥ 3, cubic sums, cubic residuosity

- **(C1) Pure complete cubic floor sums are elementary (PROVED; codex (13)).** For squarefree m and
  gcd(a, m) = 1, Σ_{h<m}⌊ah³/m⌋ = am(m−1)²/4 − (m−1)/2. This follows from the odd symmetry
  h → −h; in character language, all cubic characters are even, so Σr·ψ(r) = 0.
  - So the cubic analogue of Prop FP carries **no** class-number information, unlike degree 2.
- **(C2) Shifted full periods (PROVED).** Σ_h⌊(ah³ + c)/m⌋ reduces to incomplete cubic character sums
  −mΣ_{s∈[m−c, m)}ψ(s), over products ψ of cubic characters of the prime factors.
  - The window [0, 1) counts cube roots, ∏_{p|m, p≡1(3)}(1 + ψ_p(c) + ψ̄_p(c)) ∈ {0, 9} for
    m = pq: a cubic-residuosity decision.
- **(C3) Relation to analytic's C3 (PROVED as a statement about what is NOT implied).**
  - The complete cubic Weyl sum S(a, m) = Σ_r r₃(r)e(ar/m) is a Fourier coefficient of the
    cube-root-count function. Floor sums give its interval sums.
  - O(1)-error recovery of one Fourier coefficient from interval sums needs ≍ m intervals in general
    (codex: "O(m)-query reconstruction").
  - No polylog reduction in either direction is known (OPEN). The general cubic class factors anyway
    via its degree-2 subclass.
- **(C4) Hecke / cubic fields.** Incomplete sums Σ_{n<m/k}ψ(n) of even cubic characters are tied to
  L(1, ψχ) for odd characters χ mod k (sextic CM-type data): "hard at primes", class (b) of nonsmooth.
  - We prove no reduction to a specific cubic-field class number.
  - Q(√−3) has class number 1, so the fixed field carries no information; only varying conductors
    could.
- **Answer to Task 2.** Polylog exact degree-j floor sums for j ≥ 3 imply everything degree 2 implies
  (P = NP, factoring) and nothing new. For the restricted class, higher degree only lengthens the top
  blocks (T3).

## 6. Literature: exact lattice-point counting under curves

| Problem | Best known | Lower bound / equivalence |
|---|---|---|
| D(x) = Σ_{n≤x}τ(n) exactly | O(x^{1/3} log x): Polymath / Tao–Croot–Helfgott (Math. Comp. 81 (2012), Prop. 2.2); Sladkey arXiv 1206.3369 (O(log) space). Li, arXiv 2609.29949 (2026 preprint; UNVERIFIED): N^{9/28+o(1)} expected, via Voronoi, for **isolated full** values only | none. Li cites "circumstantial evidence" that 1/3 is characteristic of local-linear methods only. Hull vertex count Θ(N^{1/3} log N) (Alcántara–Blanco–Criado–Santos arXiv 2501.19193; Balog–Bárány) bounds enumerating methods only |
| truncated Σ_{d∈[u,v)}⌊N/d⌋ | same N^{1/3} | interval access factors N by bisection (folklore; codex) |
| Gauss circle Σ_{n≤X} r₂(n) exactly | O(X^{1/3} polylog): Har-Peled's discrete convex hull (CGTA 10 (1998)) plus Pick; or Jacobi's 4Σ(⌊X/(4k+1)⌋ − ⌊X/(4k+3)⌋) with hyperbola blocking | none; conic-sector counts factor N (p ≡ q ≡ 1 mod 4) |
| π(x) | x^{2/3} (LMO 1985), x^{1/2+ε} analytic (Lagarias–Odlyzko 1987), Õ(√x) elementary (Hirsch–Kessler–Mendlovic, arXiv 2212.09857, Math. Comp.) | none. Tao–Croot–Helfgott decide the parity of π on short intervals in x^{1/2−c}; π(x) in x^{1/2−c} is their open hope |
| Mertens M(x) | x^{3/5} (Helfgott–Thompson 2023); Õ(√x) (HKM) | none |
| polygon counts, fixed dimension | polynomial (Barvinok 1994) | — |
| quadratic boundary (parabola in a box) | — | NP-hard (Manders–Adleman), #P-hard (Thm SS) |

**Answer to Task 3.**
- No conditional lower bound ("D(x) in x^{1/3−ε} ⇒ …") is known to us or to codex.
- No equivalence between exact divisor-sum complexity and a standard hard problem is known. The only
  hardness for lattice counting under algebraic curves is NP-hardness for quadratic boundaries with
  succinct coefficients (Manders–Adleman), and its #P form above.
- These are worst-case statements over coefficient classes. They say nothing about the hyperbola
  xy = N or about the top blocks.

## 7. The sharp barrier statement

For fixed j ≥ 2, let 𝒪_j be the Taylor-class oracle of Theorem L_j: Taylor polynomials of N/(x+h) at
admissible x, rational shifts, subintervals of the block.

1. **(Uniform, PROVED.)** If 𝒪_j is computed for **all** N (or all prime N, Las Vegas):
   - in time poly(log N): FP = #P and P = NP (NP ⊆ ZPP for prime-only);
   - in time 2^{o(√log N)}: ETH fails;
   - in deterministic time T(H): every m is factored deterministically in O(log² m·T(m)). This beats
     the ladder iff T(H) ≤ H^{2/(j+4)} (power-law T).
2. **(Hard instances are not in the ladder's path, PROVED/HEURISTIC.)**
   - All constructed hard blocks (c′ ≥ b^{1−o(1)}) lie at x ≤ N^{3/8+o(1)} (PROVED for both
     embeddings: TE gives N^{3/8+o(1)} at minimal t, codex gives N^{4/11+o(1)}).
   - The hard families are solvable in H^{o(1)} (PROVED). So no statement of type "T_j ≤ H^δ ⇒ X"
     with X believed hard follows from them for any δ > 0.
3. **(Top blocks, OPEN.)** For 𝒪_j restricted to x ∈ [√N/κ, √N] of one N:
   - its consequences: Theorem L_j, and factoring with (1/(j+2))·log₂N known high bits of p;
   - its known hardness: none, beyond being factoring-type.
   - j = 2 top blocks are already inside the Coppersmith–May detection range (T2).

There is **no** reduction from a known hard problem to the restricted oracle the ladder actually uses.
The plausibility of a fast one rests on mechanism, not on complexity theory:
- the hull-walk cost per block is the affine length N^{1/3}/B, so the ladder only redistributes the
  N^{1/3} of Hiary / Tao–Croot–Helfgott (floorsum C1);
- Fourier methods need ≍ H frequencies (analytic barrier F; floorsum C4);
- one Euclid step turns a parabolic arc into a tilted conic. Conic-sector counting with full affine
  length factors via sums of two squares.

None of these is a lower bound.

## 8. Pre-registration outcomes

| Prediction | Result |
|---|---|
| P1.1: 0 failures of Prop FP (E1a) | ✓ 6 520 pairs, 0 failures |
| P1.2: decoded h(−M) correct 30/30 | ✓ 30/30 |
| P1.3: factored ≥ 29/30 | ✓ 29/30 (miss: M = 3m, ambiguous-form step returned only 3) |
| P1.4: median calls ≤ 4, max ≤ 12 | ✓ median 3, max 6 |
| P2.1 (E2a, k = 3, real floor sums) | ✓ 40/40 (32 with solutions) |
| P2.2 (E2b, k = 4–10, CRT core) | ✓ 200/200 (129 with solutions; b up to 181 bits) |
| P3.1 (E3a, j = 2, 3, 4) | ✓ 60/60 |
| P3.2 (E3b, prime N) | ✓ 20/20 |
| P3.3 (E3c, bisection factoring) | ✓ 40/40, max 32 calls ≤ bound |

- **Disclosed deviation.** The prereg said E1b used 16–22-bit m. The generator (two primes of
  ⌊bits/2⌋ and ⌈bits/2⌉ bits) produced m of 15–22 bits: three instances have 15 bits.
- **Disclosed peek.** One development smoke run before the experiments, on m < 200 and 3 subset-sum
  instances, all consistent. The pytest suite was also run before the experiments.
- All oracles are brute force. The experiments check identities and reductions, never costs. Toy
  factoring is a correctness check.

Reproduction:
```bash
.venv/bin/python -m pytest -q tests/test_floorhard.py
PYTHONPATH=. .venv/bin/python -m src.floorhard.run_experiments e1   # ~1 min
PYTHONPATH=. .venv/bin/python -m src.floorhard.run_experiments e2   # ~1 min
PYTHONPATH=. .venv/bin/python -m src.floorhard.run_experiments e3   # ~1 min
PYTHONPATH=. .venv/bin/python -m src.floorhard.run_experiments summary
```

## 9. Reviews

**GPT-6.1-sol blind round (floorhard_1).** Found independently:
- the NP-hardness via Manders–Adleman with two calls;
- deterministic bisection factoring;
- the restricted Taylor embedding (a different, one-call construction);
- Prop FP (its (6));
- the cubic elementary identity (13);
- "NP-hardness does not rule out H^{o(1)}";
- the literature table entries (Har-Peled, Li, HKM, Helfgott–Thompson as Mertens).

It did not give a #P or ETH statement, the domination corollary, or the top-block analysis. It flagged
the top-block/narrower-oracle question as OPEN, as we do.

(Adversarial round and critics: pending at this draft.)
