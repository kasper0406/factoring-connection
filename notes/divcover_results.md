# divcover: higher-rank divisor covers for the Umans–Wang route to N^{1/6}

Agent `divcover`, 2026-10-18. Pre-registration: `notes/divcover_prereg.md` (commit `77d8b67`,
made before any experiment-driver run). Code: `src/divcover/`. Tests: `tests/test_divcover.py`.
Tables: `experiments/analysis/divcover_tables.md`. GPT-6.1-sol rounds: `codex/divcover_1.{md,out}`
(blind attempt) and `codex/divcover_2.{md,out}` (adversarial review).

## Verdict

**We found no construction. No improvement over N^{1/5} is claimed.** Labels follow
AGENT_RULES: proved / experimentally supported / conjectured / refuted for this model /
UNVERIFIED (literature or search-completeness statements we did not check ourselves).

1. **The integer-factoring corollary needs much less than Conj. 5.1 states (proved, R1/R3).**
   - UW's proof of Thm 5.2 uses its covering hypothesis only for *primes* p ≤ n. With a
     Pollard–Strassen pre-step (R3), only primes in (n^{2δ}, n] with δ = max(α, β) are needed.
   - He–Sahai's refutation of the AP version of Conj. 5.1 runs entirely through semiprimes
     pq | (one term). It therefore does not apply to the weakened prime-only hypothesis (1′),
     which is all the integer application needs, not even at rank 1. HS do refute the AP
     version of Conj. 5.1 as literally stated.
   - The prefactored condition can be replaced by a combinatorial "separation" condition
     (Prop. R2; proved in the repaired form after both reviews).
   - For integer factoring the problem therefore reduces to: **make the nonzero differences
     of two low-rank GAPs of small height cover all large primes.**

2. **For that reduced problem we found no construction with β < 1/2 and any fixed α < 1,
   either in the literature or ourselves.** The literature search was done by GPT-6.1-sol
   through 2026-10-18; its completeness is UNVERIFIED.
   - The best low-rank constructions we know (subset sums of block primorials; CRT/Minkowski)
     have log H = Θ(n/log n), i.e. α = 1 − o(1) (proved).
   - UW's α = 1 − β product construction and the binomial-coefficient construction (§3) are
     not low-rank GAPs.
   - Reaching α = β = 0.39 (which would beat N^{1/5}) needs ≫ n^{0.61} elements of a
     rank-n^{o(1)} GAP, each with a divisor ≥ exp(c·n^{0.22}) composed of primes in
     (n^{0.78}, n] (Lemma M, proved).
   - We know of no mechanism producing that (conjectured to be impossible; §4).

3. **Random-model barrier (conjectured) and exact small-n data (experimentally supported).**
   - In a random-residue model, a rank-r GAP behaves like a CRT design:
     log H ≈ Σ_{p>L} log(p/L)/(r+1). With r = n^{o(1)} that forces α → 1 whenever β < 1/2.
   - The exact exhaustive searches (n ≤ 100, ranks 1–3, prime and full covers, 198
     pre-registered cells) found **no structure beyond the random model**. The median
     ρ = log H*/log H_pred is 1.04 for rank 1, and 1.17 and 1.35 for ranks 2 and 3
     (details in §7). The pre-registered
     structure signal (ρ = log H*/log H_pred ≤ 0.5) never fired.
   - Structured families at n = 10^3–10^4 also showed no signal: CRT/LLL designs give
     α_eff = 0.76–0.92, rising with n; cyclotomic/polynomial-value GAPs cover at most 0.34× the
     naive random expectation.

4. **He–Sahai at rank r (proved: an extension, plus a scoped explanation of its failure).**
   - Their linear-space lemma generalises to v ≤ 1 + Δ(μΔ − 1), where μ is the largest number
     of band primes shared by two covering elements.
   - For a rank-1 AP of length L ≤ a²n, μ = 1.
   - For rank ≥ 2, the a-priori bounds on μ (up to n^{α+o(1)}) and on Δ (up to L/L_max for
     balanced boxes) are too weak: plugged in, the lemma gives nothing beyond the trivial mass
     bound for every α < 1/2.
   - This shows the *method with these bounds* fails. It does not show that sharper bounds on
     Δ or μ for small-height covers are impossible.
   - A conditional rank-2 theorem holds: μ = 1 and short side o(n^{1/4}/√log n) give HS's bound.
   - If Δ ≲ L/√n and log H = o(√n), any rank-r n-divisor GAP at (1/3, 1/3) must contain two
     elements whose difference is divisible by ≳ n^{1/6}/log n distinct primes from [2√n/3, √n].

Bottom line:
- The Umans–Wang route remains **open but unsupported**.
- Its integer-factoring version is a prime-cover problem, which He–Sahai's argument does not
  touch.
- We conjecture that problem is as hard as strong "large-prime parts in short GAPs"
  statements.
- We have no evidence at any scale for exponents below 2/5.

---

## 1. Exact statements (from the primary sources)

Sources were read in full from the arXiv PDFs (text extraction): UW = Umans–Wang,
arXiv 2511.10851v1; HS = He–Sahai, arXiv 2608.06681v1.

- **Def. (UW 3.1).** A set A of positive integers has the *n-divisor property* if every
  i ∈ {1, …, n} divides some a ∈ A.
- **Conj. 3.2 ((α,β)-Divisor).** For infinitely many n there are S, T ⊆ Z_{>0} such that
  A = {s − t} has the n-divisor property, all elements of S, T are ≤ exp(n^α), and
  |S|, |T| ≤ n^β.
- **Conj. 3.3 (Strong).** In addition, S = S_1 + … + S_c and T = T_1 + … + T_{c'} with
  c, c' ≤ n^{o(1)}, where the S_i, T_j are arithmetic progressions (generalized APs, "GAPs").
- **Conj. 5.1 (Strong Prefactored).** Conj. 3.3 holds, and for every s ∈ S, t ∈ T the prime
  factorization of s − t can be output deterministically in Õ(n^α).
- **Mass constraint.** The product of all differences is ≥ lcm-type mass e^{(1+o(1))n}, which
  gives α ≥ 1 − 2β (UW §3). UW: "we do not even know if α < 1 − β is achievable" (for 3.2).
- **Thm 5.2 (framework).** Let n = ⌊N^{1/2}⌋ and A_1, …, A_k integers with:
  - (1) every integer ≤ n divides Π A_i;
  - (2) interval products Π_{a≤i≤b} A_i mod d (d ≤ N) cost Õ(log d · n^γ);
  - (3) the prime factors of each A_i cost Õ(n^γ).

  Then N can be factored deterministically in Õ(N^{γ/2} log k).
- **Thm 5.5.** Conj. 5.1 ⇒ deterministic factoring in Õ(N^{max(α,β)/2+o(1)}). It uses
  A_{j+ℓ|S|} = s_j − t_ℓ, Lemma 5.3 (S, T mod d in Õ(log d (n^{α+o(1)} + n^β))) and Lemma 5.4
  (multipoint evaluation).
- **AP version (UW §6, Prop. 3.4).** A single AP A of length n^{2β}, height exp(n^α), with the
  n-divisor property, implies Conj. 3.3.
- **HS Thm 1.1.** An n-divisor AP of height H with log H = o(√n) has length
  L ≥ (√(8/27) − o(1)) n^{3/4}/√log n.

  Hence the AP version fails for α < 1/2, β < 3/8 (Cor. 1.2), in particular at (1/3, 1/3).
  The proof (HS §3): primes P in [2√n/3, √n] not dividing c; blocks B_i = {p : p | u + ic};
  every pair pq ≤ n divides a common term; distinct blocks meet in ≤ 1 point because
  pq | (i − j)c is impossible for 0 < |i − j| < L ≤ pq. Then the linear-space Lemma 2.1 gives
  v ≤ Δ(Δ − 1) + 1 with Δ ≤ 1 + L/(a√n).
- **HS Remark 4.2.** "This theorem concerns the one-dimensional arithmetic-progression version
  only … For a rank-two progression, the corresponding difference contains two independent
  coefficients, and the at-most-one-intersection argument does not follow."

**Exponent bookkeeping.** Any GAP Q of box size L splits as S − T with |S|·|T| = L: split the
coordinates, and split the longest side as in Prop. 3.4. So we work with the single difference
GAP Q and β = log L/(2 log n). Beating Harvey's N^{1/5} requires max(α,β) < 2/5, i.e.
L < n^{4/5} and log H < n^{2/5}.

## 2. Reformulations (what the integer application really needs)

**Prop. R1 (prime-only cover). PROVED.** In Thm 5.2, property (1) can be replaced by
(1′): every prime p ≤ n divides some nonzero A_i.

*Proof.* The proof of Thm 5.2 uses (1) only in the sentence "for every prime p ≤ n that
divides N, it occurs as a factor of some d_i, by Property 1". That holds under (1′): p | N_0 =
gcd(Π A_i, N), and the recursive splitting keeps Π d_i = N_0. Multiplicities are recovered by
trial division (UW: "we trial divide each with N for at most O(⌈log N⌉) rounds"). The cofactor
is 1 or the unique prime > n. Zero A_i must be excluded (0 has no factorization), e.g. by
S ∩ T = ∅. ∎

*Checked against the text by GPT-6.1-sol (divcover_1 §1: PROVED).*

**Corollary (proved).** He–Sahai's theorem does not apply to the integer-factoring use of UW, even at
rank 1. Its only use of the n-divisor property is "the n-divisor property supplies a progression
term divisible by pq" (HS §3). A prime-only cover never supplies that.

**Prop. R3 (only large primes). PROVED.** Primes ≤ n^{2δ} (δ = max(α,β)) can be removed
first by Pollard–Strassen in Õ(n^δ) ring operations: take F(X) = Π_{j≤m}(X + j), m = ⌈n^δ⌉,
evaluate at multiples of m, and split by gcds. So only primes in (n^{2δ}, n] need covering.

**Prop. R2 (prefactoring can be replaced by separation). Proved, in the repaired form
below.** The first version (support-distinctness only among primes in (n^{2γ}, n]) was refuted
by the GPT-6.1-sol review (codex/divcover_2.out, C2). Its counterexample is N = 77, n = 8,
(A_i) = (2,3,5,77): the prime 11 > n has the same support as 7, so the coprime base returns 77.
Let Z_q = {i : q | A_i, A_i ≠ 0} for any prime q. Replace (3) in Thm 5.2 by:

- (3″) for every prime p ∈ (n^{2γ}, n] and every prime q ≠ p dividing Π A_i (of any size,
  including q > n), Z_p ≠ Z_q; and
  M := max_{p ∈ (n^{2γ},n]} |Z_p| ≤ n^{o(1)}.

Then (1′) + (2) + (3″) still give deterministic factoring in Õ(N^{γ/2+o(1)}·M·log k).
We also assume that no A_i is 0 (e.g. S ∩ T = ∅). A zero index would make N_0 = N and stay live
in every subtree (critic finding 2).

*Proof.*
1. Remove prime factors ≤ n^{2γ} of N (R3).
2. Compute N_0 = gcd(Π A_i, N) via (2).
3. Run UW's recursive splitting, but send d to *both* halves:
   d_lo = gcd(Π_{a..c} A_i, d) and d_hi = gcd(Π_{c+1..b} A_i, d).
   - Discard d = 1.
   - **Stop** at any d that is a prime power (perfect-power test plus deterministic
     primality test).
   - Successive gcds keep every leaf equal to gcd(A_i, N_0).
4. Cost. A non-stopped node's d is not a prime power, so it has ≥ 2 distinct prime factors.
   At most one of them exceeds n, so it contains a controlled prime p ∈ (n^{2γ}, n], and it
   lies on a root-to-leaf path of an index in Z_p.
   - Hence each level has ≤ ω(N_0)·M non-stopped nodes, plus at most twice as many stopped
     children.
   - Each node costs Õ(log N · n^γ) by (2), and there are ⌈log₂ k⌉ levels.
   - Total: Õ(n^γ · M · log² N · log k).
5. Collect the leaf values and the stopped prime powers, and compute a coprime base
   (Bernstein; quasi-linear).
6. Two distinct primes lie in the same base element only if they have identical supports
   among the inputs. That is ruled out by (3″) whenever one of them is a controlled prime.
   So every base element is a prime power, identified by a perfect-power test.
7. Multiplicities follow by trial division; the remaining cofactor is 1 or the unique prime
   > n. ∎

*Remarks.*
- (3″) is a property of the cover alone, independent of N, but it constrains *all* primes
  dividing the differences, including uncontrolled large ones.
- Two covers whose "block partitions" of the controlled primes are orthogonal (each
  controlled prime divides exactly one element of each cover, and no two controlled primes
  share an element in both) give M = 2. They satisfy (3″) except for coincidences where some
  other prime q divides both designated elements of a controlled p. Excluding those remains
  a hypothesis.
- So, for integer factoring, UW's conjecture can be used **without** the prefactored
  condition, at the price of a combinatorial separation hypothesis.
- This does not change any exponent we can reach. In every construction we know,
  prefactoring is not the binding constraint; the GAP structure is (conjectured to hold in
  general).

**Lemma M (mass, threshold form). PROVED.** Let D be the distinct nonzero differences and
R(a) = Σ_{p | a, n^{2β} < p ≤ n} log p. Prime coverage gives Σ_{a∈D} R(a) ≥ ϑ(n) − ϑ(n^{2β}) = (1 − o(1))n for fixed β < 1/2. Since R(a) ≤ n^α and
|D| ≤ n^{2β}, we get α + 2β ≥ 1. Moreover ≫ n^{1−α} differences have R(a) ≫ n^{1−2β}
(differences below a small multiple of that threshold contribute ≤ n/2). At α = β = 0.39:
≫ n^{0.61} elements of the GAP must each have a divisor composed of primes in (n^{0.78}, n] of size ≥ exp(c n^{0.22}).

**Uniformity caveat. UNVERIFIED as a formal requirement; it follows from the algorithmic use.**
Thm 5.5 needs, for the given N, a cover for some n′ ≥ ⌊√N⌋ with n′ ≤ N^{1/2+o(1)}. So "infinitely
many n" must be dense, and the GAP description (b, c_k, L_k) must be computable in
Õ(n^{max(α,β)}). In particular a construction may not enumerate the primes up to n (that costs
n^{1−o(1)} = N^{1/2−o(1)}). Every per-prime CRT design is therefore algorithmically useless even
before its height is considered.

## 3. Constructions tried, with exact (α, β)

C = prime cover of (n^{2δ}, n]; G = low-rank GAP with |S|,|T| ≤ n^β; H = height ≤ exp(n^α);
F = prefactored (or separated, R2); U = uniformly constructible in Õ(n^δ).

| Family | (α, β) achieved | C | G | H | F | U | status |
|---|---|---|---|---|---|---|---|
| Interval / Pollard–Strassen: S = {1+im}, T = {1+j}, S−T = [1,n] | (height α→0, but factoring the differences ≤ n by Pollard–Strassen costs n^{1/4}; β = 1/2) | ✓ | ✓ | ✓ | ✓ | ✓ | **proved**; gives N^{1/4} |
| UW products: partition [n] into n^β parts, S = products, T = {0} | (1−β, β) | ✓ | ✗ (rank n^β) | ✓ | ✓ | ✗ (needs all primes ≤ n) | proved (UW §3) |
| Binomials: S = {C(jk,k) : j ≤ ⌊n/k⌋+1}, k = n^α/log n, T = {0} | (α, 1−α) | ✓ (primes in (k, n]) | ✗ | ✓ | ✓ (Kummer; segmented sieve, Õ(k + √n)) | ✓ | **proved** (codex-1 + critic); not a GAP; vacuous for the reduced problem since δ ≥ 1/2 |
| Subset sums of block primorials: Q = {Σ_{k≤r} ε_k M_k}, r = 2β log₂ n, blocks of mass n/r | (1 − o(1): log H = n/(2β log₂ n), any β) | ✓ | ✓ (rank r = O(log n)) | ✗ for α<1 | ✗ (mixed sums); separation: Z_p ≠ Z_q holds with 2 orthogonal partitions but the bound M on accidental hits is not shown, and ε = 0 gives a zero element (UNVERIFIED) | ✗ | **proved** (C, G, H) |
| CRT/Minkowski: assign primes to box points, short vector of the (r+1)-dim lattice | log H ≤ ϑ(n)/(r+1) + O(log L) | ✓* | ✓ | ✗ | ✗ | ✗ | **proved** (Minkowski bound; *the short vector may make an assigned element zero, so nonvanishing needs an extra counting argument not given here); E2(a): α_eff 0.76–0.92 at n ≤ 10^4 |
| Polynomial values / cyclotomic: generators t^i, box of small-coefficient f, S−T ⊇ {f(t)} (contains t^k−1, Φ-type values) | β = 1/2 by pigeonhole (box > p forces a collision mod p); special t: none found | — | ✓ | ✓ | ✗ (Φ_m(t) not factored) | ✓ | pigeonhole part **proved** (codex-1); E2(b): coverage 0.26–0.34× naive expectation (no signal, explained by ±f and X·f redundancy) |
| Selected cyclotomic values (n^{o(1)} of them) | impossible for α < 1 | ✗ | | | | | **proved** (mass) |
| UW Lemma 6.1 / selector combination of k APs | β grows by (k log 2)/2/log n; k = O(log n) costs n^{Θ(1)} | — | ✓ | O(k n^a) | ✗ (mixed sums) | | **proved** (codex-1) |
| Number-field norms a² − Db² | inert p ∤ D: p | norm ⇒ p | a, b; boxes below p miss all inert primes | ✗ | ✗ (not linear) | | | | **proved** limitation |
| Smooth-number identities (Balog–Wooley strings) | reach smoothness X^{1/u} with u fixed or ~log_4 X/log_5 X, far from y = (log X)^{1/α} | ✗ | | | | | UNVERIFIED (literature, codex-1) |

**Exponent (order) variant (reduction proved; no gain found, conjectured).** One can replace s − t by a^s − a^t mod N.
Then p | N is caught iff ord_p(a) | s − t, and S, T mod d cost Õ(n^α + n^β) multiplications by
repeated squaring. So the *polynomial-style* full divisor property of the exponent sets S − T
(UW's Conj. 3.3, which needs no prefactoring at the integer level) would suffice. The orders,
however, must still be covered, and these include large prime divisors of p − 1. This is
order-finding BSGS / Pollard p − 1 with a structured exponent set (cf. `notes/pm1bsgs_analysis.md`).
Splitting leaves with equal orders needs the factorization of s − t again. Nothing improves.

**The core obstruction for every family (conjectured; a heuristic summary, not a theorem).** Integer factoring needs ≫ n^{1−α} GAP elements
with designed large-prime parts. Designing a part costs its full mass in height (CRT), unless
the GAP coefficients carry the primes. The coefficients are only r + 1 = n^{o(1)} numbers, so
they carry at most n^{α+o(1)} of the mass n. Algebraic identities (Fermat, Wilson, Kummer,
cyclotomic) produce prime-rich numbers, but not *linearly compressible* families of them.
Any number ≡ g(p) (mod p) with g a polynomial is ≡ g(0); Wilson-type residues m! ≡ ±1/(p−m−1)!
grow super-polynomially in p − m.

## 4. The random-model barrier

**Random model (conjectured; heuristic).** For a rank-r GAP with box L and generic coefficients, a prime
p > L divides a nonzero element with probability ≈ |Q mod p|/p ≈ L/p, independently across
primes. Then:
- Pr(cover (L, n]) = exp(−(1 − 2β + o(1))n) [ϑ-weighted Σ log(p/L)];
- the number of rank-r GAPs of height ≤ H is ≈ H^{r+1}/((r+1)! Π(L_k − 1));
- so log H_pred ≈ Σ_{p>L} log(p/L)/(r+1).

With r = n^{o(1)} this forces α ≥ 1 − o(1) unless β = 1/2. It matches the CRT/Minkowski
construction up to the assignment entropy.

**EXPERIMENTALLY SUPPORTED (small n).** See E1 in §7. The exact minimal heights follow
log H* ≈ log H_pred for rank 1, and lie slightly above it for ranks 2–3, over all 198
pre-registered cells with n ≤ 100. No cell comes close to the pre-registered signal ρ ≤ 0.5.
The (r+1)·log H*_r column is roughly rank-independent, as the model predicts. So small-n
optima show no hidden structure: higher rank buys exactly the (r+1)-th root of the CRT mass,
and nothing more.

**Literature analogues (codex-2 C5; literature statements UNVERIFIED by us).**
- The rank-1, c = 1 prime-cover problem is: find b ≥ 1 with b mod p ∈ {0, −1, …, −(L−1)} for
  every prime p ∈ (L, n]. The CRT density of such b is exactly Π_{L<p≤n} L/p =
  exp(−(1 − 2β + o(1))n) (proved).
- Placing the *least* solution near the inverse density is the heuristic step.
- The only unconditional lower bound we have is the elementary mass bound
  b + L ≥ exp((1 + o(1))n/L), which is superpolynomial but far from exp(cn).
- Pseudosquares are the closest analogue (least x-pseudosquare Q(x)):
  - unconditional lower bound Q(x) ≫ x^{4√e − ε} (Schinzel via Burgess);
  - under ERH, Q(p) > exp(√(p/2)) (Lukes–Patterson–Williams);
  - the random-model size exp(cx) is conjectural.
- Pseudosquare arguments exploit multiplicative subgroups via characters, while our residue
  sets are additive intervals. No located theorem gives a stronger obstruction here.
- So an unconditional refutation of the prime-cover version, even at rank 1, looks out of
  reach of current methods.

## 5. He–Sahai at rank r

**Prop. HS-r. PROVED.** Let P be a set of v ≥ 2 primes and E a finite set of nonzero integers
such that every pair {p, q} ⊆ P divides a common e ∈ E, and no e is divisible by all of P.
Let B_e = {p ∈ P : p | e}, Δ = max_p #{e : p | e}, and μ = max_{e≠e′} |B_e ∩ B_{e′}|
(≤ max ω_P(e − e′)). Then

  v ≤ 1 + Δ(μΔ − 1).

*Proof.*
- Fix e and choose w ∈ P ∖ B_e. Each q ∈ B_e co-occurs with w in some e′, and e′ ≠ e
  because w ∉ B_e. Hence B_e ⊆ ∪_{e′∋w, e′≠e}(B_e ∩ B_{e′}) and |B_e| ≤ μΔ.
- Fix u ∈ P. Then P ∖ {u} ⊆ ∪_{e∋u}(B_e ∖ {u}), so v − 1 ≤ Δ(μΔ − 1). ∎

With μ = 1 this is HS Lemma 2.1. (Property-tested on 3000 random block systems in
`tests/test_divcover.py`.)

**Application to a GAP.** Let Q = {b + Σ x_k c_k} have box L, longest side L_max, height H,
and band P ⊆ [a√n, √n] with p ∤ c_{k_max}. Then:
- Δ ≤ (L/L_max)·⌈L_max/p⌉ ≤ L/L_max + L/(a√n);
- μ ≤ log(2H)/log(a√n), since band primes shared by two elements divide their difference,
  which is an element of Q − Q of size ≤ 2H.

**Why the method fails at rank ≥ 2 (proved for the a-priori bounds only).**

(a) **μ > 1 is possible.** Example (codex-1): the proper box b = pq, c = 1, e = pq − 1,
m = k = 2 has elements {pq, pq+1, 2pq−1, 2pq}, containing two multiples of pq at height 2pq.
In general a rank-2 GAP's difference set Q − Q is itself a rank-2 GAP of height 2H, and its
elements may hold up to log(2H)/log(a√n) band primes. For rank 1, Q − Q = {jc}, and with
|j| < L ≤ a²n at most one band prime divides j: μ = 1 automatically. This is exactly HS's
"(i − j)c" step.

(b) **The degree bound collapses for balanced boxes.** A prime may divide one element per
line in the longest direction, so Δ ≲ L/L_max. For a balanced rank-2 box L_max = √L ≪ √n, and
even with μ = 1 the lemma gives only √L ≳ √v, i.e. L ≳ √n/log n (weaker than trivial). For
balanced rank r the bound is n^{r/(4(r−1))−o(1)} ≤ n^{1/2} (codex-2).
This shows only that the a-priori bound is weak. We have not shown that Δ ≈ L/L_max is
attained by a small-height n-divisor GAP. One hit per line for most band primes forces the
ratio of the generators mod p into a short window for most p, which is a CRT constraint with
its own mass cost, so a sharper Δ bound may exist (critic finding 3).

Qualification (codex-2 C3): a *bounded* μ alone would not hurt. Any constant μ keeps the HS
exponent 3/4 whenever Δ is AP-like (≈ L/√n). The real losses are μ growing polynomially (up to
n^{α+o(1)}) and unfavourable box geometry. Rank ≥ 2 does not prevent stronger bounds for
specially controlled GAPs, e.g. the conditional theorem below.

**Exponent consequence (proved).** With the a-priori bounds (taking even the favourable
Δ ≈ L/√n; the worst case L/L_max is weaker still) the inequality gives, for every
α < 1/2 and rank ≥ 2, nothing beyond the trivial mass bound L ≥ n^{1−α−o(1)}:
- with Δ ≈ L/√n and μ ≈ n^α/log n, it gives L ≳ n^{3/4−α/2};
- 3/4 − α/2 ≤ 1 − α ⟺ α ≤ 1/2.

**Conditional rank-2 theorem (proved; codex-1 §4, checked here and by the Opus critic).** Let
A = {b + ic + je : i < m, j < k} ⊂ [1, H] be proper, with full n-divisor property,
log H = o(√n), m ≥ k. Suppose:
1. no nonzero element of A − A is divisible by two distinct primes in [2√n/3, √n] (μ = 1);
2. k = o(n^{1/4}/√log n).

Then L = mk ≥ (√(8/27) − o(1)) n^{3/4}/√log n.

*Proof:* Δ ≤ k + 3L/(2√n), then HS-r with μ = 1.

**Necessary structure at (1/3, 1/3) (corollary of HS-r; proved).** Let a rank-r n-divisor GAP
have log H = o(√n), L ≤ n^{2β}, 2β < 3/4, and degree Δ ≤ C·L/√n on the band. Then it has two distinct
elements whose difference is divisible by ≥ v/Δ² ≳ n^{3/2−4β}/log n band primes. At β = 1/3
that is n^{1/6}/log n. The height bound permits up to 2n^{1/3}/log n, so this is no
contradiction. Under the stated Δ hypothesis (and log H = o(√n), so that blocks are proper), this is the
feature a rank-r construction must have and HS's argument
cannot exclude.

## 6. Best guess and open question

**CONJECTURED (low confidence).**
- The uniformly constructible prime-cover version (hence the Strong Prefactored conjecture) is
  false for every max(α, β) < 2/5. The random model predicts no GAP below α ≈ 1 with β < 1/2,
  and nothing known beats that.
- A proof of falsity at rank 2 would require a sum-product-type statement: a proper rank-2 GAP
  of size n^{2β} and height exp(n^α) cannot have total (n^{2β}, n]-radical mass ≥ (1 − o(1))n.
  We believe this is open even at rank 1 (pseudosquare-type difficulty; conjectured).

**Most promising next step** (shared with codex-1): bound the weighted target-prime radical
Σ_{a∈Q} R(a) for proper rank-2 GAPs, first at α = β = 0.39 and without prefactoring. An
o(n) bound would be the first genuine prime-cover obstruction. A matching construction would
be a breakthrough. The exact search code here can test any proposed rank-2 family at small n.

## 7. Computations

### E1 — exact minimal heights (pre-registered)

Full tables: `experiments/analysis/divcover_tables.md`. Raw data:
`experiments/raw/divcover_exact.jsonl` (198 cells). Compute: 6.1 CPU-hours on 24 cores.

**Censoring.**
- Prime mode: 99 cells; 6 had a shape time out, 5 have no value.
- Full mode: 99 cells; 39 timed out, 40 have no value, and 2 of those are genuinely
  "none ≤ Hmax = 200·H_pred".

Censored cells are never imputed.

**Pre-registered structure signal** (prime, n ≥ 40, L ≤ n/4, ρ = log H*/log H_pred ≤ 0.5):
**NONE.** Within L ≤ n/4 the smallest ρ is 0.87 (n = 100, L = 24, rank 1: H* = 3622 against
H_pred = 12925). Below 0.87 there are only L > n/4 cells: ρ = 0.71 at n = 100, L = 32 with
H* = 194 ≈ 2n. That is the trivial interval-type regime, which is why the criterion excludes
it.

**Experimentally supported (n ≤ 100):** exact optimal prime-cover heights follow the
random-model prediction:

| rank | uncensored cells (n ≥ 40) | min ρ | median ρ | max ρ |
|---|---|---|---|---|
| 1 | 33 | 0.71 | 1.04 | 1.21 |
| 2 | 25 | 0.92 | 1.17 | 1.38 |
| 3 | 18 | 1.06 | 1.35 | 1.54 |

- Rank 1 sits on the model.
- Ranks 2–3 are slightly *above* it. The model probably overcounts distinct rank-2/3 GAPs
  (non-proper boxes, symmetric parameters, residue collisions); this explanation is post hoc.
- Higher rank lowers heights roughly as the (r+1)-th root (table "Rank comparison"), and
  never beyond.
- Full covers need markedly larger heights than prime covers at the same (n, L): ρ_full
  measured against the prime model is ≈ 1.5–2.5. So the composite targets that He–Sahai
  exploit are already the binding constraint at small n.

**Exponent fits (descriptive only; pre-asymptotic, as pre-stated).**
- At α = 1/3 the cap exp(n^{1/3}) ≤ 104 for n ≤ 100. No grid cell reaches it in either mode,
  so L*(n; 1/3) is undefined on the whole range.
- At α = 1/2, prime mode, L*(n) gives fitted slopes 0.950 / 0.966 / 0.966 (ranks ≤ 1/2/3),
  i.e. β_fit ≈ 0.48.
- At α = 1/2, full mode, the fits give β_fit ≈ 0.43, but on censored and partial data.
- None of these is evidence about asymptotics.

**Calibration remark.** He–Sahai's n^{3/4}/√log n is an asymptotic lower bound, with
hypothesis log H = o(√n), and it is not computable at n ≤ 100. Our rank-1 data are consistent
with it, but do not test it. At these sizes the height caps exp(n^α) are polynomial in n.

### E2 — structured families at n = 10^3, 10^4 (pre-registered)

Raw: `experiments/raw/divcover_structured.jsonl`.

**(a) CRT/LLL designs.**
- Ranks 1–4, L ≈ n^{0.5..0.8}, n = 10^3 and 10^4. log H agrees with mass/(r+1) to within 0.3%.
- α_eff ranges over 0.760–0.890 at n = 10^3 and 0.824–0.923 at n = 10^4. It *increases*
  with n, as the theory predicts (α_eff → 1).
- 31 of 32 designs were verified as covers. The exception (n = 10^4, r = 1, L = 631) left 2
  small primes uncovered after the pre-set 6 repair rounds; it is recorded as failed.
- **Pre-registered E2 success: NONE.**

**(b) Polynomial-value / cyclotomic GAPs.**
- Covered fraction / naive random expectation is 0.26–0.34 for every non-vacuous (t, r, K).
  The t = 3 rows are 0: their elements are smaller than every target prime.
- 25 of the 40 cells are vacuous because the box exceeds n.
- **Pre-registered signal (≥ 2× at n = 10^4): NONE.**
- Post hoc: the shortfall below 1× comes from the family's redundancy. f and −f, and f and
  X·f, have the same prime divisors, which cuts the effective number of distinct values by
  roughly 3×.

## 8. Reviews

1. **GPT-6.1-sol blind attempt** (`codex/divcover_1.out`). It independently:
   - confirmed R1;
   - found no construction below 2/5;
   - corrected our first mass statement (Lemma M threshold form);
   - supplied the binomial construction, the norm limitation, the selector cost and the
     conditional rank-2 theorem.
   - Its best guess: the prime-only conjecture is false below 2/5 (low confidence).
2. **GPT-6.1-sol adversarial review** (`codex/divcover_2.out`).
   - C1 CORRECT.
   - C2 (our first R2) GAP: the prime > n can share a support (counterexample N = 77), and
     can blow up the recursion. **Repaired** as (3″) with prime-power stops (§2).
   - C3 CORRECT with qualifications, now incorporated: a bounded μ alone does not hurt; for
     balanced boxes the bound is n^{r/(4(r−1))}.
   - C4: construction calculation correct. "None known" was downgraded to "none found in the
     literature search".
   - C5: density calculation correct; the least-solution prediction is a heuristic; pseudosquare
     bounds were added.
3. **Independent Opus 5.5 critic** (read-only subagent). It reviewed the pre-repair draft:
   - It independently found the same R2 gap (prime > n), plus zero elements (now an explicit
     hypothesis).
   - It found overclaims in "why rank ≥ 2 escapes". Now scoped: these are a-priori bounds
     only, and a sharper Δ bound may exist.
   - It found a wrong prefactoring cost for binomials (Õ(k + √n), not Õ(k)) and an off-by-one
     in the index range.
   - It flagged label problems; the labels are now mapped to AGENT_RULES.
   - It verified R1, R3, Lemma M, HS-r, the GAP Δ/μ bounds, the conditional rank-2 theorem
     and the §1 transcriptions as CORRECT.

   Its remaining open points:
   - the subset-sum row's separation needs an M bound (marked UNVERIFIED);
   - the CRT row's nonvanishing of assigned elements (footnoted).

No claim of improvement over N^{1/5} is made, so the mandatory improvement-claim review does
not apply. The reviews above were done anyway.

## 9. Reproduction

```
.venv/bin/python -m pytest -q tests/test_divcover.py          # 39 tests
.venv/bin/python -m src.divcover.run_exact --jobs 24           # E1 (≈1–2 h on 24 cores)
.venv/bin/python -m src.divcover.run_structured                # E2 (≈1 min)
.venv/bin/python -m src.divcover.analyze                       # tables -> experiments/analysis/divcover_tables.md
```
The exact search is `src/divcover/gapsearch.c`, compiled on demand by `src/divcover/core.py`.
Usage: `gapsearch n mode(0=prime,1=full) Hmax L1 [L2 [L3]]`.
