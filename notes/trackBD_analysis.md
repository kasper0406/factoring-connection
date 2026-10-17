# Track B + D: quotient sketches, lifts, annihilators — analysis

> **r3 audit applied 2026-10-17 (agent `repair`).** Corrections are marked "Corrected 2026-10-17 after r3 audit". The adjudication is in `notes/repair_adjudication.md`. Reproductions: `src/repair/audit_trackBD.py`.

Agent `sketchBD`, 2026-10-17. Code: `src/trackBD/*`, `src/quotient_sketches.py`; tests: `tests/test_trackBD.py`
(23 tests, all pass). Raw data: `experiments/raw/trackBD_*.jsonl`; configs: `experiments/configs/trackBD_*.json`;
aggregated tables: `experiments/analysis/trackBD_tables.txt` (produced by `src/trackBD/analyze.py`).

Status labels follow AGENT_RULES. "We" = this agent; no novelty is claimed.
Review: one independent Opus critic pass (read-only) on §§2–5 and the code; its 15 findings (Euclid "iff" too strong,
squarefree hypothesis in B3, Model-M overclaim, non-equidistribution of q mod (p−1) for the balanced generator,
BM step and Stirling-window details, D1 hypotheses) are incorporated below.

---------------------------------------------------------------------------------------------------

## 0. Summary of results

| # | Statement | Status |
|---|---|---|
| B1 | **Exact single-sketch criterion.** For any monic f in Z[X], any N = p^k m with p not dividing m: S := F_N mod (f, N) satisfies S ≡ F_m(X)^{p^k} (mod p, f). Hence S(α)=0 at a root α of f mod p iff (1+α)^m = 1+α^m. | proved |
| B2 | Readout corollaries: content gcd, resultant, multiplication rank, generic Hankel rank and the *degree of* gcd(f, S) reduce to the *vanishing pattern* {α : (1+α)^m = 1+α^m} on the roots of f mod each prime. Res(f,S) ≡ Res(f,F_q) (mod p). (Euclid can also fail on non-generic remainder sequences: one-way implication only.) | proved |
| B3 | For squarefree N, column-gcd elimination over Z/N returns a factor iff the pivot profiles mod the primes differ (up to a 'stuck' branch of probability ≤ (1 − φ(N)/N)^64. That is ≤ 2^-64 only when φ(N)/N ≥ 1/2, e.g. odd semiprimes, and **not** for every squarefree N: N = 105 gives ≥ 1.0·10^-17. Corrected 2026-10-17 after r3 audit). **Profile split ≠ rank split in general:** [3 5] over R_15 has rank 1 mod 3 and mod 5, but pivots {1} vs {0}, so elimination returns 3. Rank split ⇔ profile split holds for multiplication-by-S in the power basis (Krylov), checked on 3 000 random cases. Euclid over Z/N returns a factor iff the remainder-degree sequences differ; in particular whenever deg gcd(f,S) differs (needs the constant-remainder check). | proved; 0 violations in 99 000 matrices |
| B4 | The user's hypothesis, made precise: the vanishing condition depends on q only through q mod L_α, L_α = lcm(ord α, ord(1+α)) dividing p^d − 1. "Frobenius coincidence" q ∈ ⟨p⟩ (mod L_α) is **sufficient but not necessary**. | proved |
| B5 | The non-Frobenius ("sporadic") roots **outside U** form orbits of the anharmonic group {x, 1/x, −1−x, …} of order 6 (on F_p ∪ {∞}). Every sporadic count observed (2 960 prime-side censuses in §6.2 and 22 178 (p, e) pairs in §6.3) is ≡ 0 (mod 6). **Corrected 2026-10-17 after r3 audit:** with B4's literal definition (exclude only 0, −1 and FROB) the count is *not* ≡ 0 mod 6. For p = 7, q = 11 the roots 2, 4 (the cube roots of unity) are non-FROB and give count 2. This fails in 1 039 of 5 260 (p, e) pairs with p < 400. SPORADIC must exclude U = {0, −1, ω, ω²}. The six-orbit structure is classical (Cauchy–Mirimanoff; Boyd–Martin–Thom 2015, Lemma 5.2). | proved (p > 3, U excluded); the census is new data, the symmetry is not |
| B6 | On balanced random semiprimes the Frobenius part is a heavy-tailed, Pollard/Bach–Shallit-type term (mean ≈ 2–6 roots, median 0, P(>0) ≈ 0.2–19 %); the sporadic part has mean ≈ 1.0 and P(>0) ≈ 15 %. Single-sketch structural success per random weight is exactly (#non-public F_p-roots of F_q)/p (+ the q-side term); its observed average is Θ(1/p), the same order as the random-gcd null. | exact formula proved; Θ(1/p) and the means are experimental (16–40 bits) |
| B7 | Rank-gap probability of multiplication-by-sketch on F_l[X]/(X^r − c), uniform admissible c: exactly predicted by the F_l-roots of F_m *on the measured instances*; r-independent on those instances (r = 2..16); zero for 65–86 % of balanced instances. **Corrected 2026-10-17 after r3 audit: r-independence is NOT a theorem.** For N = 209 = 11·19, c = −1: gcd(X^6+1, F_19) = X^4 − X^2 + 1 mod 11, gcd(X^6+1, F_11) = 1 mod 19, and P_gap = 0 at r = 2 but 1/9 at r = 6, due to extension-field roots (own enumeration and `exp_rankdist.side` agree). Codex's r = 2 zeros (101·113, 107·193, 1009·1481) are the special case "no non-universal F_l-roots on either side", and persist for every r ≤ 16. | exact enumeration; F_p-part proved, extension-field part heuristic (observed ≤ 2 weights per (size, r)) |
| B8 | Lucas-support gaps (coefficient ≡ 0 mod p because no support exponent folds there) occur when q − p ≲ r/2 (balanced, q < 2p; threshold approximate): a weak special case of Fermat's method. | mechanism proved, threshold heuristic; all 24 observed instances at ≥ 24 bits have q − p ≤ 162 |
| D1 | The divided quantities we construct from the N²-lift (F_N'/N, Hasse D^[k]F_N/N for k below the least prime, ∂_c of the sketch for r, c units) are unit multiples of G_k = (1+X)^{N−k} − X^{N−k} or are determined by S itself. G_k vanishing at α is a pure multiplicative-order condition (Fermat test of β = 1+1/α). | proved |
| D2 | Codex claim (c): the scalar G_1 split probability (d−1)(p+q−2(d−1))/(pq), public root −1/2, nothing left when d = 2. | verified exactly (census: #roots = g−1 on all 1 120+360 instances; measured vs predicted hits agree) |
| D3 | BM over Z/N on c ↦ F_N(c): for e_p := q mod (p−1) ≥ 3 the linear complexity mod p is exactly e_p, mod q it is p (for e_p = 1 the p-side complexity is 0). For the balanced generator (q < 4p/3) e_p = q − p + 1, so this is a Fermat-gap method with cost Θ(q − p). Window predicate: Δ^K F_N(0) = (K+1)!·S(N,K+1) (K ≥ 1) is ≡ 0 mod p for all K ≥ e_p, and is a proper nonunit for e_p ≤ K ≤ p−1 up to accidental q-side zeros. | LC and window proved; **new (repair, proved): a nonunit discrepancy occurs at index ≤ p + e_p − 1, and any nonunit at index ≥ 2e_p reveals p** (60/60 runs). "≈ 2e_p" remains experimental (58 at 2e_p, 1 at 2e_p − 1; 31/90 runs hit an accidental nonunit first); LC_p = e_p in 90/90. **Corrected 2026-10-17 after r3 audit:** the Δ^K zero predicate is *not* the prefix predicate. It is non-monotone (N = 391: gcds 391, 17, 391 at K = 12, 13, 14) |
| — | Codex claim (a) (cyclotomic factors of F_N are exactly X+1 for odd N and Φ_3 with multiplicity 2 / 1 for N ≡ 1 / 5 mod 6) | verified for all 3 ≤ N ≤ 160, Φ_m with m ≤ 40 |
| — | Codex claim (d) (rank gap ⇒ nonunit pivot; F_{pq} mod p = F_q^p) | agrees with B1/B3 (we prove the sharper "iff" for column-gcd elimination with N squarefree) |

**Bottom line (refuted for model M, §5; scope corrected 2026-10-17 after r3 audit: Model M as written has NO proved lower bound. Arbitrary public R_N-polynomial readouts include Pollard p−1, e.g. N = 77 and the constant 1×1 matrix [2^6 − 1] = [63] returns 7. "Refuted" refers to the specific readouts B1–B3 analysed, as measured.)** In the model "compute sketches F_N mod (f, N) for public f of degree
r = poly(n), then any R_N-linear-algebra readout that sees only the vanishing pattern", the structural
success probability per sketch on balanced random semiprimes is (number of non-public roots of F_q in the
roots of f mod p)/p + (same with p, q swapped), which is Θ(1/p) on average and exactly 0 on a majority of
instances for the binomial family X^r − c. For every construction we analysed (M1–M5, Hankel, BM over the
weight, lifts) the split condition is a coincidence of the power map x ↦ x^q on F_{p^d} at public points
(proved). On the p-side the power map is x ↦ x^q = Frob(x)·x^{q−p}, with δ = q − p a "random" integer in
(0, p/3) for this generator; empirically the resulting coincidences behave like those of a random condition
(rate ≈ 1/p per condition). That no other polynomial readout of F_N-sketches beats Θ(1/p) is a conjecture
(§5 (iii)), not a theorem.

---------------------------------------------------------------------------------------------------

## 1. Notation

N = p^k m with p prime, p ∤ m. For a positive integer e put F_e(X) = (1+X)^e − 1 − X^e ∈ Z[X] (so F_N is the
brief's polynomial). For N = pq (p < q primes) the relevant cofactor is m = q (for the prime p) and m = p
(for the prime q). For a monic f ∈ Z[X] of degree r, the **sketch** is S = F_N mod (f, N) ∈ R_N[X]/(f),
computed by repeated squaring with O(n) multiplications in R_N[X]/(f) (`sketch.sketch`).
Universal (public) points: U = {0, −1, ω, ω²}, ω a primitive cube root of unity. W = X(X+1)(X²+X+1).

## 2. Track B: the exact criterion

### Theorem B1 (proved)
In F_p[X], F_N(X) = F_m(X)^{p^k}. Consequently S ≡ F_m^{p^k} (mod p, f), and for every root α ∈ F̄_p of f,
S(α) = F_m(α)^{p^k}; in particular S(α) = 0 ⇔ (1+α)^m = 1 + α^m.

*Proof.* In characteristic p the map u ↦ u^{p^k} is a ring endomorphism of F_p[X], so
F_m(X)^{p^k} = ((1+X)^m)^{p^k} − 1 − (X^m)^{p^k} = (1+X)^N − 1 − X^N. Reduce mod f and evaluate. ∎

For a prime power N = p^k (m = 1) F_1 = 0, so every sketch has all coefficients divisible by p
(the PRIMEPOW class; `test_prime_power_always_vanishes_mod_p`). Checked exactly in
`test_theorem_B1_frobenius_identity` (semiprimes, 3- and 4-prime N, prime powers, p²·m) and
for every irreducible component of every sketch of the tiny run (§6.1: 0 violations).

### Corollary B2 (readouts; proved)
Let N be squarefree and f squarefree mod p.
1. *Content gcd:* p | gcd(N, all coefficients of S) ⇔ f | F_m in F_p[X] ⇔ every root of f mod p is a root of F_m.
2. *Resultant:* Res(f,S) ≡ Res(f, F_m)^{p} ≡ Res(f,F_m) (mod p) (Res(f, g^p) = Res(f,g)^p; Fermat on the integer). For N = pq, Res(f, F_q) = q^r Res(f, H_q) with H_q = F_q/q ∈ Z[X], so p | Res(f,S) ⇔ p | Res(f,H_q) ⇔ f and H_q share a root mod p. (Checked: `test_resultant_corollary`.)
3. *Rank:* over F_p, multiplication by S on F_p[X]/(f) has rank r − deg gcd(f, F_m^p) (kernel = (f/gcd)·A). Hankel matrices [ℓ(X^{i+j}S)] (size L ≥ rank) have rank ≤ this, with equality unless ℓ lies in a proper subvariety (probability O(r/p) for random ℓ); codex §4.3 gives the exact pairing for ℓ = top coefficient. This item holds for any f.
4. *Single coefficient gcds (classification, not a theorem):* we call a zero coefficient mod p structural if it is explained by (1) or by a Lucas-support gap (B8); other zero coefficients are labelled accidental, and their frequency matches the uniform-residue null (§6.7) — experimentally supported.

### Proposition B3 (proved; N squarefree)
Column-gcd elimination (`linalg.elim_profile`: at column j compute g = gcd(N, remaining entries); g = N skip,
1 < g < N output, g = 1 choose a unit pivot or a unit row combination) returns a proper factor iff the
reduced-row-echelon pivot sets of M mod the different primes of N differ.
*Proof.* While all pivots are units the row operations are invertible over Z/N, hence over each F_p, and
the processed columns are pivot columns for every p simultaneously. At column j the remaining entries are,
mod each p, the coordinates of column j outside the span of earlier pivot columns; they all vanish mod p
iff j is a non-pivot column mod p. So g ∉ {1, N} exactly when j is a pivot column for some primes and not
others; and if no such j exists the profiles coincide. ∎ (Squarefreeness is needed: N = 175, M = diag(35, 1)
has equal profiles mod 5 and 7 but returns 35. ~~The random-combination branch fails with probability ≤ 2^-64
and then returns no factor.~~ **Corrected 2026-10-17 after r3 audit:** the proof is conditional on finding unit
pivots. A uniform random combination of a column whose entries generate the unit ideal is uniform on R_N, so it is a unit
with probability exactly φ(N)/N (Monte Carlo 0.3837 vs 0.3836 at N = 15015). The 64 trials therefore fail with probability
≥ (1 − φ(N)/N)^64, which exceeds 2^-64 whenever φ(N)/N < 1/2. Examples: N = 105 (≥ 1.0·10^-17), N = 15015 (≥ 3.6·10^-14). The 2^-64 bound
holds for odd semiprimes, where φ(N)/N ≥ 8/15. In the stuck case the routine returns no factor.)
**Rank split vs profile split (r3 audit).** "Factor ⇔ pivot-profile split" (this proposition) is not "factor ⇔ rank split".
[3 5] over R_15 has rank 1 mod both primes, pivots {1} mod 3 and {0} mod 5, and elimination returns 3. The two coincide for
Krylov-type matrices such as multiplication-by-S in the power basis (M1 below).
Euclid over Z/N (`sketch.euclid_readout`) is weaker: it returns a factor iff the *remainder-degree sequences*
mod the primes differ (a nonunit leading coefficient is a degree drop mod some primes only). A difference in
deg gcd(f, S) forces this (provided the *final constant remainder* is gcd-checked: our first version lacked that
check and missed every r = 2 split; fixed and regression-tested), but abnormal remainder sequences also occur
with equal gcd degrees (e.g. N = 101·113, f = X^5 + 1768X^4 + 8330X^3 + 10875X^2 + 7191X + 1416 returns 113 at a
degree-1 remainder) — these are our Euclid "ACC" hits.
Experimental check: 99 000 matrices (§6.4), "factor returned" ⇔ "profiles differ" with 0 exceptions.

### Proposition B4 — where the "Frobenius coincidence" sits (proved)
Let α ∈ F_{p^d} \ {0,−1} be a root of f mod p, L_α = lcm(ord α, ord(1+α)) (a divisor of p^d − 1). Then
(1+α)^e = 1 + α^e depends on e only through e mod L_α, and the solution set E_α ⊂ Z/L_α contains the
Frobenius subgroup ⟨p⟩: for e ≡ p^i (mod L_α), (1+α)^e = (1+α)^{p^i} = 1 + α^{p^i} = 1 + α^e. Hence
  S(α) = 0  ⇔  m mod L_α ∈ E_α  ⇐  m ∈ ⟨p⟩ (mod L_α)   ("Frobenius coincidence", FROB),
and the converse fails exactly on E_α \ ⟨p⟩ (SPORADIC). For d = 1 the FROB condition reads
L_α | gcd(p−1, q−1) = g, i.e. α and 1+α both lie in μ_g — a Pollard/Bach–Shallit-type order condition
on the *other* prime; its count over α ∈ F_p is the cyclotomic number #{x : x, 1+x ∈ μ_g} = g²/p + O(√p)
(informative only for g ≫ p^{3/4}).

So the user's hypothesis is **confirmed as a sufficient condition and refuted as an iff**: the exact
statement is "S(α) = 0 ⇔ q mod L_α ∈ E_α", and on balanced semiprimes the observed successes at ≥ 28 bits are
all SPORADIC (§6.2–6.3). Note that "F_N ≡ 0 iff (1+X^p)^q ≡ 1+X^{pq} in R_p" as literally stated is a
tautology (that is F_N mod p); the content is the reduction to F_q(α) and to q mod L_α.

### Proposition B5 — anharmonic symmetry (proved)
For odd e, Z_e = {x ∈ F̄_p : (1+x)^e = 1 + x^e} is stable under x ↦ 1/x (for x ≠ 0; on the punctured or projective
domain, r3 audit) and x ↦ −1−x, hence under the order-6 anharmonic group. *Proof.* If x ≠ 0 and (1+x)^e = 1+x^e, divide by x^e: (1+1/x)^e = 1 + (1/x)^e.
And (1 + (−1−x))^e − 1 − (−1−x)^e = (−x)^e − 1 + (1+x)^e = −x^e − 1 + (1+x)^e = 0. ∎
*Corollary (proved, p > 3, e = q mod (p−1)):* the SPOR count on F_p is ≡ 0 (mod 6), **with SPOR defined as Z_e ∩ F_p minus U minus FROB**.
(Corrected 2026-10-17 after r3 audit: B4 literally defines SPORADIC as E_α \ ⟨p⟩ for α ∉ {0, −1}, which keeps the cyclotomic orbit
{ω, ω²}. For p = 7, q = 11, e = 5 the non-FROB roots are exactly 2, 4 = ω, ω², giving count 2. The proof below already treats {ω, ω²}
as universal. The census excluded U, so the recorded data stand.) *Proof.* The FROB
condition x^{e−1} = 1 = (1+x)^{e−1} is orbit-invariant (e−1 is even: (1/x)^{e−1} = 1, (1+1/x)^{e−1} =
((1+x)/x)^{e−1} = 1, (−1−x)^{e−1} = (1+x)^{e−1}, (1+(−1−x))^{e−1} = (−x)^{e−1} = x^{e−1}), so SPOR is a union of
orbits. An orbit has fewer than 6 points iff it contains a fixed point of a non-identity element: the
involutions fix {1, −1}, {−1/2, ∞}, {0, −2}, the 3-cycles fix the roots of X²+X+1. The resulting small
orbits are {0, −1, ∞} and {ω, ω²} (universal) and {1, −2, −1/2}, which lies in Z_e iff 2^{e−1} = 1 and is then
FROB. Hence every SPOR orbit has exactly 6 points. ∎ Observed for every prime in every census (§6.2, §6.3).
(These are the root symmetries of the Cauchy–Mirimanoff polynomials. **Verified 2026-10-17 (r3 audit lead, primary source read):**
D. W. Boyd, G. Martin, M. Thom, *Squarefree values of trinomial discriminants*, LMS J. Comput. Math. 18 (2015) 148–169, §5. Lemma 5.2
gives the same six-element action on Z/pZ ∪ {∞} with exceptional orbits {0, −1, ∞} (trivial), {1, −2, −1/2} (Wieferich) and
{ζ, ζ^{−1}} (cyclotomic). §5 traces the S_3 action to Helou's Cauchy–Mirimanoff polynomials and to Cauchy (1839). Their polynomial
is f_p(x) = ((x+1)^p − x^p − 1)/p, not our cross-prime F_e mod p, but the group action is identical. The symmetry is classical. Our census is new data only.)
The count #Z_e ∩ F_p is the entry δ(1,1) of the differential spectrum of the power permutation x ↦ x^e of F_p
(e = q mod (p−1) is a unit since q > p is prime) — a connection to the cryptographic literature on
differential uniformity of power maps (not verified against a primary source in this session).

### Corollary B6 — what the weight family X^r − c can see (proved, modulo extension terms)
For f = X^r − c with c uniform mod N (admissible: Res(X^r − c, W) a unit), the p-side rank of
multiplication by S is deficient iff some r-th root of c is a non-universal root of F_q. Splitting roots
by field of definition: the F_p-roots contribute exactly the weights c ∈ {x^r : x ∈ A_p},
A_p := (Z_q ∩ F_p) \ U, which has size ≤ |A_p| *independent of r*; roots in F_{p^d}, d ≥ 2, contribute weights
c = α^r ∈ F_p only for α with α^{r(p−1)} = 1, and were observed for ≤ 2 weights per instance (B7 table).
Hence P(exactly one side deficient) = a_p + a_q − 2a_p a_q with a_l = |{x^r : x ∈ A_l} ∩ adm|/#adm + (ext); the
profile-split probability is this plus P(both deficient with different gcd degrees) (a lower bound; the two
coincided in all our data). The F_p-part satisfies a_l ≤ |A_l|/#adm (the bound is r-independent: proved). The **image**
{x^r : x ∈ A_l}, the admissible set and the extension part all depend on r (heuristically O(r/p) weights; observed ≤ 2 per (size, r) on our instances).
So "vary the degree" (codex's repair) cannot raise the F_p-part beyond |A_l|/#adm, and it did not raise P_gap on the measured instances (B7).
**Corrected 2026-10-17 after r3 audit:** degree independence of the success probability is **false as a general statement**.
For N = 209 = 11·19 and c = −1, the extension roots X^4 − X^2 + 1 | gcd(X^6 + 1, F_19) mod 11 give P_gap = 1/9 at r = 6 versus 0 at r = 2
(`src/repair/audit_trackBD.py: n209`).

## 3. Track B question 3: weighted/paired sketches as linear algebra over R_N

All matrices below are built from N and public data only; we measured rank-profile splits (§6.4) and
derived their criteria.

* **M1 multiplication by S** (rank r − deg gcd(f, F_q^p) mod p; Cor. B2.3). For f squarefree mod p and q:
  split ⇔ deg gcd(f, F_q) mod p ≠ deg gcd(f, F_p) mod q (the pivot profile of a multiplication matrix is
  {0,…,rank−1}, so profile split = rank split); public-root factors present on one side only (TD) included.
  Prediction checked exactly: M1 rank mod p and mod q equal the B1 prediction in 19 800/19 800 cases.
* **M2 Hankel/Toeplitz of projected coefficients** u_k = ℓ(X^k S mod f): rank = degree of the minimal
  polynomial of the sequence ≤ rank M1; split ⇔ M1-split or ℓ accidentally degenerate mod one prime
  (prob ≈ r/p). Measured split rate ≈ 1.5× M1's at every size.
* **M3 Berlekamp-type Q_N: g(X) ↦ g(X^N) mod f** (columns X^{iN} mod f). In evaluation coordinates at the roots
  α_i of f mod p, Q_N becomes the Vandermonde at the points α_i^N = (α_i^q)^p, so rank Q_N mod p =
  #{distinct α_i^q} (proved for f squarefree mod p: a Vandermonde with r columns at s distinct points has rank s).
  A collision α_i^q = α_j^q needs ord(α_i/α_j) = q, i.e. q | p^D − 1 for D ≤ deg of the splitting field —
  an order coincidence of heuristic probability O(D/q) over random p. **M3b = Q_N − I**: ~~a kernel jump needs α^N = α (ord α | pq − 1,
  a Fermat condition) or a coincidence among the 2r points {α_i, α_i^N}.~~ **Corrected 2026-10-17 after r3 audit: that necessity
  is false.** A kernel jump needs a polynomial relation g(α_i^N) = g(α_i) for all roots, not a point coincidence. Example: N = 247 = 13·19,
  f = X³ + 5. Mod 13 the roots 2, 5, 6 have N-th powers 11, 8, 7 (six distinct points, no fixed point), yet Q_N = diag(1, −1, 1), so X² is an
  extra kernel vector of Q_N − I (α^N = −α for every root). Mod 19, Q_N = diag(1, 5, 6): rank split, factor found. Measured: 0 splits at ≥ 28 bits.
  For X^r − c with c a unit and gcd(r, N) = 1, X^{iN} ≡ c^{⌊iN/r⌋} X^{iN mod r}: Q_N is a monomial matrix with unit
  entries, full rank on both sides; for Φ_m with gcd(m, N) = 1 it is the Galois automorphism ζ ↦ ζ^N (unimodular).
  Q_N never splits for these families (Q_N − I can: diagonal entries c^k − 1 with public k — a Pollard-type test).
* **M4 Krylov/Hankel over the weight c** (sequence c ↦ [X^j](F_N mod X^r − c)): see D3 — the linear
  complexity mod p equals the degree+1 of the reduced function c ↦ P_j(c) on F_p, which is e_p = q mod (p−1)
  (= q − p + 1 for this generator) for r = 1 and ≈ p for r = 2 (measured LC_p = p in 30/30 cases): rank split only for Hankel size ≳ e_p.
* **M5 multiplication by G_1** (Track D comparison): split ⇔ some root α with β = 1+1/α of order dividing
  pq − 1 on one side only (Fermat in the extension). Rate ≈ (τ(p−1)−2)/p for d = 1 components.

**Why no rank split for generic balanced N ("both components generic").** Every matrix above has entries that
are R_N-polynomials in the sketches; mod p its rank is a function of (f mod p, the values F_q(α) at the roots).
The only q-dependence is through the power map x ↦ x^q on the finite groups F_{p^d}^*, i.e. through
q mod (p^d − 1). (Caveat: this residue is *not* equidistributed for our instances — q < 4p/3, so q mod (p−1) =
q − p + 1 and q mod (p^d − 1) = q for d ≥ 2; x^q = Frob(x)·x^{δ}, δ = q − p. Our rate statements therefore rest on
per-instance exact censuses, not on a Dirichlet argument.) A split therefore needs
a public polynomial relation that holds for the specific power map x ↦ x^{q mod (p^d−1)} on the p-side and not
for x ↦ x^{p mod (q^d−1)} on the q-side. Relations that hold *structurally* are exactly the Frobenius
coincidences (e ∈ ⟨p⟩: the power map is F_p-linear on the relevant subgroup); all others are "sporadic"
coincidences of a fixed algebraic condition, empirically of probability ≈ (number of conditions)/|F_{p^d}|.
(Status: the reduction to q mod (p^d − 1) is proved; "sporadic ≈ 1/|F_{p^d}| per condition" is
conjectured / experimentally supported: mean sporadic count ≈ 1.0 per prime in the balanced censuses, §6.2.)

**What a rank-split construction needs (structural property).** An R_N-algebra A built from N in poly time and an
operation Ψ on A (here Ψ = N-th power) such that, with inverse-polynomial probability over the hidden primes,
Ψ is linear (a Frobenius power, or zero, or the identity) on a nontrivial piece of A_p and not on A_q. For
Ψ = x ↦ x^N on a unit group G this means |G(F_p)| (or the order of a known element) divides q − p^i — a
group-order coincidence between the two primes. Classical methods obtain the analogous coincidence by
*choosing the exponent* (Pollard p−1, Williams p+1, Bach–Shallit cyclotomic: exponent = smooth number,
success when |G(F_p)| is smooth); with the exponent pinned to N, no choice is left and the coincidence
probability is Θ(1/p)·(heavy-tailed constant). A polynomial-time variant would need a group (or algebra)
with |G(F_p)| | E(N) for an N-computable E and a 1/poly fraction of p — not available from these sketches.

## 4. Track D: lifts, Hasse derivatives, divided quantities

### Proposition D1 (proved)
1. In Z[X], F_N' = N·G_1 with G_k := (1+X)^{N−k} − X^{N−k}. Computing (1+X)^{N−1}, X^{N−1} mod (N², f),
   forming N(…) and dividing by N returns G_1 mod (N, f) (`lifted_derivative_divided`,
   `test_lift_derivative_equals_G`). G_1 is directly computable mod N without the lift.
2. Hasse derivatives: D^[k]F_N = C(N,k)·G_k (k ≥ 1). For gcd(k, N) = 1, C(N,k)/N = C(N−1,k−1)/k, so
   D^[k]F_N / N ≡ (C(N−1,k−1)/k)·G_k (mod N); for k below the least prime of N the scalar is (−1)^{k−1}/k, a unit
   (for larger k it can be a nonunit, e.g. N = 35, k = 11: C(34,10) ≡ 0 mod 5 — a trial-division-sized event)
   (`test_hasse_divided_is_unit_times_G`).
3. Quotient terms do not commute with differentiation, but they carry nothing new mod N: writing
   F_N = Q·(X^r − c) + S over Z, Q mod (X^r − c) = ∂S/∂c (fold of t·c^{t−1}), and F_N' = N G_1 gives
   ∂_X S = −r X^{r−1} ∂_c S (mod N, X^r − c) (`test_quotient_identity_dX_S_equals_minus_r_X_dc_S`, computed
   with dual numbers c+ε). For r and c units mod N this gives c·∂_c S_j = −(j/r) S_j (j ≥ 1) and r·∂_c S_0 ≡ 0:
   the quotient fold is determined by S. Mod N² it is determined by S mod N² and G_1 mod N.
4. *Remark (not a theorem):* "(·)^N − (·) divided by N" is only defined when the difference is ≡ 0 mod N, which for composite N and
   generic arguments it is not (it is exactly the event the gcd tests detect); no further divided quantity
   arises.
So every divided quantity *among those constructed here* is a unit multiple of G_k or is determined by S.

### Proposition D2 (proved; confirms codex (c))
For α ∉ {0, −1} in F_{p^d}: G_k(α) = 0 ⇔ β^{N−k} = 1, β = 1 + 1/α. For d = 1 (scalar evaluation) the roots are
β ∈ μ_{gcd(N−k, p−1)} \ {1}, i.e. gcd(m − k, p−1) − 1 roots; for k = 1, g − 1 roots with g = gcd(p−1,q−1),
including the public root α = −1/2 (β = −1). Uniform scalar sampling splits N with probability
(g−1)(p+q−2(g−1))/(pq). Census: #roots = g − 1 for every prime of all 1 120 balanced + 360 diagnostic instances;
measured hit counts agree with this formula at every size (§6.2).
So G_k is a *pure multiplicative-order* (Fermat-type) condition, vanishing "more often" than F_N's two-element
condition: mean nonpublic roots τ(p−1) − 2 when e is uniform over units (Menon's identity; B2c: e.g. p = 1093: 22) versus F_N's ≈ 1 sporadic + heavy-tailed FROB.
Still Θ(τ(p−1)/p): exponential. This is the classical gcd(a^{N−1} − 1, N) Fermat-witness phenomenon.
For Carmichael N, F_N(c) ≡ 0 modulo every prime for every c (e_p = 1 ⇒ F_1 = 0): the scalar F-sketch never
splits (0 hits in 66 000 evaluations); G_1 splits only through c ≡ 0, −1 modulo some prime (public TD).

### Proposition D3 — Berlekamp–Massey and the finite-difference window
Let N = pq, p < q odd primes, e = e_p := q mod (p−1) (odd, 1 ≤ e ≤ p−2). For the balanced generator
(q < 4p/3) e_p = q − p + 1 exactly.
*(proved)* The function c ↦ F_N(c) mod p equals the polynomial (1+c)^e − 1 − c^e (exponent reduction for
c ≢ 0, −1; at c = 0, −1 both sides are 0 since q, e are odd). For e ≥ 3 it has degree e − 1 and leading
coefficient e ≢ 0; for e = 1 (i.e. p−1 | q−1, e.g. p = 3 or Carmichael-type pairs) it is identically 0 and any
F_N(c) with c ≢ 0, −1 (mod q) reveals p. Mod q the polynomial is F_p, of degree p − 1 < q. Hence for e ≥ 3 the
sequence s_c = F_N(c), c = 0, 1, 2, … has linear complexity exactly e mod p and p mod q, and BM over Z/N
(`linalg.bm_mod_N`, gcd on every nonzero discrepancy) ~~must meet a nonunit by step ≈ 2e~~. **Corrected 2026-10-17 after r3 audit**
(the two infinite-sequence complexities do not imply the step). **Proposition (repair, proved).** For e ≥ 3, BM over Z/N on s_c meets a
nonunit discrepancy at some index t ≤ p + e − 1, and every nonunit at index t ≥ 2e is ≡ 0 mod p and ≢ 0 mod q, so it reveals p.
*Proof.* Before any nonunit, every discrepancy is a unit or 0 mod N, so the runs mod p and mod q make identical decisions. They are genuine BM
runs over F_p and F_q with equal lengths L_t. The p-side has LC e, so after 2e terms its LFSR is final and all later p-discrepancies vanish. From index 2e on,
the q-side LFSR has length e and generates the prefix. Let t be its first failure, where the q-discrepancy is ≠ 0 and the p-discrepancy is 0, so the nonunit reveals p.
By Massey's lemma, any LFSR generating s_0..s_t has length ≥ t + 1 − e. The minimal q-side LFSR (length p) generates everything, so t ≤ p + e − 1.
Such t exists because LC_q = p > e. ∎ Checked: 60/60 balanced 16/20-bit runs, max step/(p+e−1) = 0.47 (`audit_trackBD.bm_step_bound`).
"≈ 2e" holds when the q-side has no accidental zero discrepancy near 2e. That is heuristic (≈ 1/q per step) and experimental.
*(experimental)* The first structural nonunit occurred at step 2e (58 runs) or 2e − 1 (1 run; s_0 = 0 lets the
profile converge early); 31/90 runs met an accidental nonunit first.
*Window (proved, K ≥ 1).* Δ^K F_N(0) = Σ_{i≤K} (−1)^{K−i} C(K,i)[(1+i)^N − i^N] = (K+1)!·S(N, K+1)
(Stirling numbers of the second kind; from Δ^K x^N|_{x=1} = K!S(N,K) + (K+1)!S(N,K+1)). Since Δ^K kills
polynomials of degree < K: Δ^K F_N(0) ≡ 0 (mod p) for all K ≥ e, ≢ 0 for K = e − 1 (e ≥ 3); mod q it equals
Δ^K F_p(0) = (K+1)!·S(p, K+1), which vanishes for K ≥ p and for K ≤ p − 1 only accidentally. So
gcd(Δ^K F_N(0), N) = p for e ≤ K ≤ p − 1 (up to q-side zeros S(p, K+1) ≡ 0 mod q, whose number we do not bound; calling them
"accidental" is a label, not a bound) and = N for K ≥ p. ~~The predicate "Δ^K F_N(0) ≡ 0 mod N" is the brief's prefix predicate Z(N, ·) in disguise;~~
**Corrected 2026-10-17 after r3 audit: WRONG.** The Δ^K zero predicate is non-monotone. For N = 391 = 17·23 (e = 7) the gcds at K = 12, 13, 14 are
391, 17, 391, whereas Z(N, ·) is monotone. It is a sibling of the prefix predicate, not the same predicate; the nonunit window is [q − p + 1, p − 1]
for this generator (length ≈ 0.9p on average here: e_p/p ≤ 0.33, mean ≈ 0.11 over 1 120 instances).
This is also a sketch: Δ^K F_N(0)/K! is the top coefficient of F_N mod ∏_{i=0}^{K}(X − i) — a *dense* modulus of
degree K + 1 > e. Evaluating Δ^K F_N(0) for one K ≈ √N faster than Θ(K) is an open primitive (we know no
algorithm even matching the factorial's Õ(√K)). The BM route itself costs Θ(q − p) terms here, i.e. it is a
(much weaker) Fermat-gap method: Fermat's method needs ≈ (q−p)²/(8√N) steps. Tests:
`test_bm_linear_complexity_over_c`, `test_stirling_window`.

## 5. The model and what it excludes

**Model M (single-sketch / sketch-linear-algebra).** Input N. The algorithm chooses public monic f_1, …, f_s ∈ Z[X]
of degree ≤ r = poly(n) and coefficients ≤ 2^{poly(n)} (possibly adaptively, from N and previous outputs),
computes S_i = F_N mod (f_i, N), and outputs gcds of N with: coefficients, resultants, Euclid / inversion
leading coefficients, pivots of any Gaussian elimination of matrices whose entries are fixed R_N-polynomials
in the S_i coefficients (multiplication, Hankel, Krylov, Berlekamp-type), and BM discrepancies of sequences of
such entries.
**Scope (corrected 2026-10-17 after r3 audit).** As written, Model M has **no proved factoring lower bound**. "Fixed R_N-polynomials" includes
constants and polynomials of exponential degree computed by repeated squaring, so Pollard p−1 lives inside M. Example: N = 77 and the allowed constant
matrix [2^6 − 1] = [63] gives elimination factor 7. A meaningful lower-bound model must restrict readout degree and constants, e.g. a fixed degree and
coefficients independent of N, and even then nothing is proved. Everything below is about the specific readouts named.
**What is proved in M.** (0) Every M-readout mod p is a function of f mod p and the values F_q(α)^p at the
roots α of the f_i (B1); the q-dependence is through q mod (p^d − 1) only. Readouts are *not* restricted to the
vanishing pattern: e.g. N = 1009·1481, f = X − 60722 has F_N(c) ≡ 2 (mod 1009) and the 1×1 matrix [S − 2]
reveals 1009 — a level-set test whose success rate (≈ 1/p per tested level) we classify as accidental.
(i) For the readouts content gcd, resultant, deg gcd(f, S) (and hence Euclid's forced failures), multiplication
rank (M1) and generic-ℓ Hankel rank (M2), a split for the prime p happens iff the vanishing pattern
{roots α of f_i mod p : (1+α)^q = 1+α^q} differs between the primes (Lucas-support gaps, B8, affect single coefficients only); each
α-event depends only on q mod L_α and is implied by, but not equivalent to, a Frobenius coincidence (B1–B4).
(ii) For Q_N, Q_N − I (M3) and BM over the weight (D3) the exact conditions are the order/degree statements of
§3 and D3. (iii) For an *arbitrary* matrix whose entries are fixed R_N-polynomials in several sketches, mod p the
rank is a function of the values F_q(α) at the public algebraic points, so a split requires a public algebraic
relation among these values that holds on one side only; that such relations are, apart from public roots and
Frobenius coincidences, satisfied with probability O(#conditions/|F_{p^d}|) is **conjectured** (supported by
every experiment here, not proved).
**Measured rate (experimental):** per random weight, (|A_p| + |A_q|)/p-ish with E|A_p| ≈ 1 (sporadic) + heavy-tailed FROB term;
per-sketch split rate ≈ 2.3–2.6/p (B7) — i.e. the accidental-gcd null 1/p + 1/q up to a constant ≈ 1.2–1.3.
**Excluded by M (intended, not enforced by the written definition; see Scope above):** exponents other than N (Pollard/Williams/ECM-type smooth exponents E — these change the
coincidence condition from q mod L_α to "L_α | E"); quotient degrees ≥ e_p or ≈ p (D3 needs them);
non-polynomial-size coefficient growth; algebraic groups other than unit groups of R_N[X]/(f) (elliptic curves,
tori with N-dependent twists); interaction with lattice methods (Track E) using partial information; any
adaptive use of a *non-vanishing-pattern* invariant (e.g. exact values of Res(f,S) mod p as data for
Coppersmith) — we did not find such a use, but M does not exclude it.

## 6. Experiments (exact arithmetic throughout; commands reproduce the raw files)

```
.venv/bin/python -m pytest -q tests/test_trackBD.py
.venv/bin/python -m src.trackBD.exp_tiny --limit 5000 --R 200 --M 100 --procs 12 --out experiments/raw/trackBD_tiny.jsonl
.venv/bin/python -m src.trackBD.exp_scalar --bits 16 20 24 28 32 --inst 200 --K 4000 --procs 12 --out experiments/raw/trackBD_scalar.jsonl
.venv/bin/python -m src.trackBD.exp_scalar --bits 36 40 --inst 60 --K 4000 --procs 12 --out experiments/raw/trackBD_scalar.jsonl
.venv/bin/python -m src.trackBD.exp_scalar --family diag_bigg --bits 24 32 --s 1 2 4 8 16 32 --inst 30 --K 4000 --out experiments/raw/trackBD_scalar_diag.jsonl
.venv/bin/python -m src.trackBD.exp_scalar --family carmichael --inst 33 --K 2000 --out experiments/raw/trackBD_scalar_carm.jsonl
.venv/bin/python -m src.trackBD.exp_expect --bits 8 10 11 12 --count 12
.venv/bin/python -m src.trackBD.exp_rankdist --procs 8
.venv/bin/python -m src.trackBD.exp_rank --bits 16 20 24 28 32 40 --inst 150 --procs 8
.venv/bin/python -m src.trackBD.exp_bm --bits 16 20 24 --inst 30 --procs 8
.venv/bin/python -m src.trackBD.exp_sketch_bal --bits 16 20 24 28 32 --inst 100 --K 100 --procs 10
.venv/bin/python -m src.trackBD.exp_sketch_bal --bits 36 40 --inst 30 --K 100 --procs 10
.venv/bin/python -m src.trackBD.analyze > experiments/analysis/trackBD_tables.txt
```
(`trackBD_tiny.jsonl` and `trackBD_rank.jsonl` were gzipped after the run; `analyze.py` reads either form.)
Instances: `src.common.instances.balanced_semiprime(bits, seed)` (seeds 0..inst−1; at 16 bits there are only
~66 distinct semiprimes of this shape, so 16-bit rows contain repeats). Weights c come from
`random.Random(f"…|{N}|{seed}…")` (public). Total CPU ≈ 2–3 core-hours.

### 6.1 Exhaustive tiny composites (B1 run), all 4 ≤ N < 5000 (4 329 composites)
Moduli: X^r − c (r = 1..200, c ∈ {1, −1, 2, 3, 5}), Φ_m (m = 1..100), 20 fixed random monic squarefree f of
degree 2..8. Readouts: coefficient gcds, content gcd, Euclid. Regime "cmp" = deg f < smallest prime p
(compressed), "unc" otherwise (the sketch then holds most of F_N; for N < 5000 most (N, r) are uncompressed).
Theorem B1 was checked on every irreducible component of every sketch: **0 violations**.
Classes per hit (diagnostic): PRIMEPOW; TD (p divides a public integer of f: f(0), f(−1), disc f, Res(f, Φ_3));
FULL_x / SPLIT_x (component vanishing, labelled FROB / SPORADIC); SUPPORT (Lucas gap, BD-6); ACC (none of these).

Odd semiprimes, compressed regime (hits per class; one sketch can count in several readouts):

| family (sketches) | readout | hits | TD (of which via Res(f,Φ3)) | FULL/SPLIT FROB | FULL/SPLIT SPORADIC | SUPPORT | ACC |
|---|---|---|---|---|---|---|---|
| X^r − c (57 540) | coef | 28 924 | 3 033 (1 338) | 2 513 | 15 | 12 077 | 11 286 |
| | content | 5 080 | 2 552 (857) | 2 513 | 15 | – | – |
| | Euclid | 32 364 | 3 869 (2 180) | 5 057 | 771 | – | 22 671 |
| Φ_m (22 547) | content | 2 570 | 1 684 (259) | 876 | 10 | – | – |
| | Euclid | 11 719 | 2 199 (259) | 1 043 | 218 | – | 8 259 |
| random f (11 304) | content | 944 | 784 (434) | 160 | 0 | – | – |
| | Euclid | 5 308 | 2 301 (1 190) | 952 | 104 | – | 1 951 |

Readings. (1) For tiny p the FROB class dominates the structural hits (5 057 vs 771 for Euclid) because
g = gcd(p−1, q−1) is comparable to p; at ≥ 28 bits it disappears (§6.2). (2) Most Euclid and coefficient hits
in the tiny range are ACC or SUPPORT (p ≤ 70: accidental zeros have probability ≥ 1/70 per coefficient).
(3) Codex claim (a), "most early successes are W/Res(X^r−c, Φ_3) baseline discoveries", is **not** what we
see on the tiny exhaustive set: public-integer (TD) hits are 10–20 % of Euclid/coefficient hits and Res(·,Φ_3)
is about half of TD; for *scalar* sketches on balanced N, however, W-roots (c ≡ 0, −1, ω mod p) do explain
≈ 70 % of hits (§6.2) — of which only ≈ 1/3 come from Φ_3. Both statements are consistent with codex's
filtering advice; the cyclotomic-factor statement itself is verified (§0).
(4) Prime powers: every hit is PRIMEPOW (F_1 = 0); 97–100 % of sketches split them. Odd N with ≥ 3 primes
behave like semiprimes with more chances. Carmichael N < 5000 (561, 1105, 1729, 2465, 2821): scalar sketches
(r = 1) vanish modulo every prime (Korselt ⇒ e_p = 1) and never split; higher-degree components do split
(SPLIT_FROB 24 of 110 compressed binomial sketches), via the FROB condition in F_{p^d}, d ≥ 2.
(5) Even N: p = 2 makes 0 ≡ −1 and almost every sketch splits trivially (TD / SPLIT_TRIVIAL).
Full tables: `experiments/analysis/trackBD_tables.txt`.

### 6.2 Scalar sketches F_N(c), G_1(c), G_2(c) (B2/D2), K = 4000 random c per instance
Classes (diagnostic, for the prime revealed): TD_0 / TD_m1 / TD_phi3 = p divides the public integer c, c+1,
c²+c+1 (public W-roots); FROB / SPOR as in B4; G1: ORDER (β^{N−1} = 1) / TD_half (p | 2c+1). "pred" = exact
prediction from the per-prime root census, combined by CRT (a + b − 2ab).

| bits | inst | mean p | F hits: TD (0,−1,Φ3) | FROB | SPOR | all (pred) | G1 ORDER | G1 TD_half | G1 all (pred) | NULL hits |
|---|---|---|---|---|---|---|---|---|---|---|
| 16 | 200 | 213 | 6891+7099+7858 | 1835 | 5392 | 29075 (29220) | 27208 | 7003 | 34211 (34349) | 7023 |
| 20 | 200 | 853 | 1859+1772+1703 | 21 | 1399 | 6754 (6786) | 5553 | 1678 | 7231 (7320) | 1789 |
| 24 | 200 | 3409 | 466+436+415 | 186 | 375 | 1878 (1965) | 3997 | 474 | 4471 (4402) | 478 |
| 28 | 200 | 13622 | 95+104+119 | 0 | 127 | 445 (467) | 710 | 112 | 822 (895) | 119 |
| 32 | 200 | 54575 | 17+26+43 | 0 | 33 | 119 (114) | 272 | 30 | 302 (315) | 27 |
| 36 | 60 | 218322 | 3+0+3 | 0 | 2 | 8 (8.5) | 27 | 3 | 30 (28.6) | 0 |
| 40 | 60 | 879733 | 1+0+1 | 0 | 0 | 2 (1.9) | 6 | 0 | 6 (3.7) | 0 |

Readings: (i) the exact criterion predicts the hit counts at every size (within Poisson error);
(ii) ≈ 70 % of F hits are public (W-roots: c ≡ 0, −1, ω mod p); (iii) the non-public F rate ≈ the NULL rate
(e.g. 32 bits: 33 SPOR vs 27 NULL; 28 bits: 127 vs 119); (iv) at ≥ 28 bits no FROB hit occurred: the
Frobenius coincidences are rare-instance events; (v) G_1 non-public rate ≈ (g−2)/p ≈ 10× NULL (mean g ≈ 12).
Every sporadic count in the census is ≡ 0 (mod 6) (400/400 primes per size).

Diagnostic family `diag_bigg` (q = 1 + t(p−1)/s, so g ≥ (p−1)/s): FROB dominates and the rate follows
the cyclotomic count ≈ g²/p² per c, e.g. 32 bits: s = 1: 90 030 hits (pred 89 999) of 120 000; s = 4: 12 570
(12 573); s = 16: 896 (915). This validates the FROB term of B4 across three orders of magnitude.

### 6.3 Exact expectation over the hidden cofactor (B2c, `exp_expect`)
Caveat: this averages over *all* units e mod p−1, which is not the distribution of e_p = q − p + 1 for our
balanced generator (e_p ≤ 0.33p); the balanced censuses of §6.2 are the instance-relevant numbers (they agree:
mean non-public F_p-root count ≈ 0.7–1.4).
For 43 primes p ∈ [137, 4027] and every unit e mod p−1:
FROB count mean 2.0–6.3 (median 0, P(>0) 0.2–19 %), SPOR count mean 0.0–1.75 (≈ 1.0 for p > 1000; P(>0) ≈ 15 %,
always ≡ 0 mod 6), G_1 non-public count mean τ(p−1) − 2. The F_N single-weight structural success rate averaged
over q is thus (E[FROB] + E[SPOR])/p ≈ (3–7)/p, dominated by rare q with large gcd(p−1, q−1).

### 6.4 Rank splits (B3) — 99 000 matrices
Moduli per instance: X^r − c (r ∈ {2,3,4,6,8,12,16,24}, small public c ∈ [2,50]), Φ_m (m ∈ {5,7,9,12,15,21}),
8 random f; 150 instances per size; no public-root filtering (so TD splits are included).

| matrix | 16 bits | 20 | 24 | 28 | 32 | 40 |
|---|---|---|---|---|---|---|
| M1 mult. by S | 57 | 5 | 4 | 0 | 1 | 0 |
| M2 Hankel ℓ(X^{i+j}S) | 91 | 12 | 6 | 0 | 2 | 0 |
| M3 Q_N | 12 | 2 | 1 | 0 | 0 | 0 |
| M3b Q_N − I | 100 | 18 | 12 | 2 | 0 | 0 |
| M5 mult. by G_1 | 87 | 12 | 24 | 1 | 0 | 0 |

(splits out of 3 300 matrices each). Prop. B3 ("factor ⇔ profile split") held in all 99 000 matrices; M1 ranks
mod p and mod q equalled the Theorem B1 prediction in 19 800/19 800 cases. Rates fall like 1/p (≈ ×4 per +4 bits)
— no construction holds a constant split rate.

### 6.5 Exact rank distributions over the weight c (B3b/B7)
For each N and r, all admissible c mod p and mod q enumerated (`exp_rankdist`); P_gap exact under uniform c mod N.

* Codex's N = 101·113, 107·193, 1009·1481: **P_gap = 0 for every r ∈ {2,…,8,12,16}**, because neither F_p
  contains a non-universal root of F_q nor vice versa (A_p = A_q = ∅) and no extension root occurs. This is the
  r-independent special case of B6/B7.
* N = 47959 = 199·241 (codex's repaired example): A_199 is a single 6-orbit, A_241 = ∅; P_gap = 1422/46215 ≈ 3.08 %
  for every r (6 bad weights mod 199 out of 195 admissible), as predicted.
* Balanced 16/20/24-bit (28/37/40 distinct instances): mean P_gap is r-independent (24 bits: 0.00071–0.00078 for
  r = 2..16), equals the F_l-root prediction to 3 digits, extension-field roots contributed ≤ 2 weights in total per
  (size, r); the fraction of instances with P_gap > 0 is 14 % / 35 % / 30 %; mean P_gap·p ≈ 0.8–1.8 / 2.2–2.3 / 2.3–2.6.

### 6.6 Berlekamp–Massey (D3/D5)
30 instances per size; T = 768 / 2560 / 4096 terms at 16 / 20 / 24 bits.

| sequence | bits | STRUCT | ACC | none | notes |
|---|---|---|---|---|---|
| S1: F_N(c), c = 0.. | 16 / 20 / 24 | 18 / 21 / 20 | 12 / 9 / 10 | 0 | LC_p = e_p in 90/90; STRUCT step = 2e_p (58) or 2e_p − 1 (1) |
| S2: [X^1](F_N mod X²−c) | 16 / 20 / 24 | 0 / 2 / 0 | 30 / 28 / 25 | 0 / 0 / 5 | LC_p = p (no structure) |
| S3: [X^0](F_N mod X^m−1) | 16 / 20 / 24 | 0 | 28 / 30 / 27 | 2 / 0 / 3 | maximal LC; ACC at ≈ p/2–p/4 steps |

STRUCT = the connection polynomial at the event annihilates the whole window mod the revealed prime. Every
STRUCT event needs ≈ 2 e_p = 2(q − p + 1) terms (Θ(p) on average for this generator); ACC events (a discrepancy vanishing mod p by chance, ≈ 2/p per step)
arrive after Θ(p) steps as well. Both are exponential in n.

### 6.7 Degree-r sketches on balanced semiprimes (B2b)
100 instances per size at 16–32 bits, 30 at 36/40; per instance 100 random weights for each
r ∈ {1,2,3,4,6,8,12,16,32,64,128,200}, plus Φ_3..Φ_60 and 20 random f. "pred" = per-sketch split probability
from the F_l-roots (Cor. B6); "null" = r·(1 − φ(N)/N − 1/N).

| bits (mean p) | r | sketches | Euclid SPLIT (pred) | coef ACC (null) | TD | SUPPORT |
|---|---|---|---|---|---|---|
| 20 (854) | 2 | 10 000 | 20 (23) | 41 (45) | 58 | 0 |
| 20 | 16 | 10 000 | 26 (23) | 324 (359) | 58 | 495 |
| 20 | 200 | 10 000 | 24 (23) | 478 (4 490*) | 73 | 8 559 |
| 24 (3 395) | 2 | 10 000 | 13 (10) | 12 (11) | 17 | 0 |
| 24 | 8 | 10 000 | 7 (9) | 47 (45) | 17 | 0 |
| 24 | 200 | 10 000 | 9 (9) | 763 (1 130) | 14 | 1 799 |
| 28 (13 602) | 4 | 10 000 | 5 (1.7) | 7 (5.6) | 4 | 0 |
| 32 (55 164) | 8 | 10 000 | 0 (0.3) | 5 (2.8) | 1 | 0 |
| 32 | 200 | 10 000 | 1 (0.3) | 60 (70) | 2 | 200 |
| 36 (218 253) | all r | 3 000 each | 0 (≈0.03) | ≤ 8 | 0 | 0 |
| 40 (892 300) | all r | 3 000 each | 0 (≈0.004) | ≤ 1 | 0 | 100 (one instance, q − p = 54) |

(*for r = 200 at 20 bits the null exceeds the observed ACC count because SUPPORT instances absorb hits.)
SPLIT counts track the prediction; coefficient ACC counts track the uniform-residue null; SUPPORT appears in
whole-instance blocks exactly when q − p is small (BD-6). Φ_m and random f behave the same way (e.g. 24 bits:
Φ_m Euclid SPLIT 10 / ACC 91 of 5 800; random f SPLIT 3 / ACC 6 of 2 000).

## 7. Cost audit
Per sketch: O(n) multiplications in R_N[X]/(f) (q_mul ≤ 2n, each O(r²) R_N-ops schoolbook; flint uses
fast multiplication), plus O(r) gcds for coefficient readouts, O(r) Euclid steps with one gcd each, O(r³)
for elimination. Everything is polynomial *per attempt*; the success probability per attempt is Θ(1/p)
(B6/B7), so expected work is Θ(p·poly(n)) = exponential in n — no better than trial division in exponent,
and worse than Pollard rho (Θ(√p)). The BM/finite-difference route (D3) needs Θ(e_p) = Θ(p) terms.
Nothing in this track grows like √p.

## 8. Most promising open question (for the orchestrator)
The criterion shows the information in F_N-sketches is exactly "does x ↦ x^q behave additively at α".
A genuinely new route must change *the operation*, not the readout: find an R_N-algebra A (poly-size,
built from N) and an N-computable operation Ψ that is forced to act F_p-linearly on A_p for a 1/poly fraction
of primes p independently of q. Within this track the sharpest concrete primitive left is D3's:
evaluate Δ^K F_N(0) = Σ_{i≤K} (−1)^{K−i} C(K,i)[(1+i)^N − i^N] = (K+1)!·S(N, K+1) mod N at a single
K ≈ √N/2 in time o(K^{1/2}). Its "nonunit" window [e_p, p − 1] has constant density in [1, √N], like the
factorial's window [p, q) for gcd(K!, N) on balanced N, so it is a sibling of the Track C predicate, not an
improvement; we do not know an algorithm for it even matching the factorial's Õ(√K). Belongs with Track C;
we rate it low priority.
