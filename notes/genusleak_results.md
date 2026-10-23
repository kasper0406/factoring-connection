# genusleak: family (Q) without smoothness (phase 24)

Agent `genusleak`, 2026-10-23. Pre-registration: `notes/genusleak_prereg.md` (commit fdef700). Code:
`src/genusleak/` (operational: `forms.py`, `cf.py`, `cayley.py`, `sqrtred.py`, `genus.py`; harness: `harness.py`).
Tests: `tests/test_genusleak.py`. Raw data: `experiments/raw/genusleak_E{1..5}.jsonl`.

Throughout, N = pq with p < q odd primes, k ≥ 1 squarefree with gcd(k, N) = 1, D = −4kN (imaginary case) or
D = 4kN (real case), and h = h(D). "Hidden 2-torsion" means a class of Cl(D)[2] outside the subgroup A_pub generated
by the classes that k alone determines (Lemma 0).

## 0. Verdict (plain language)

**Negative, as a search outcome. Among the routes and methods we reviewed, we found none that reaches a hidden
ambiguous class (or a useful ambiguous/square form) in subexponential time without a group order, a relation, or a
smoothness event.** We did not prove that no such route exists. The barriers below hold only in their stated models.
Arithmetic uses of form coefficients remain unbounded by them.

What is proved or checked, route by route:
- **Genus characters (Prop GEN).** Deciding whether a class lies in the principal genus (= Cl²) is equivalent to
  quadratic residuosity mod N, by randomized reductions. One direction uses GRH, with a small factored multiplier k.
  QR_N is not known to be easy, nor known to be as hard as factoring. Even full genus data does not locate the
  involutions.
- **Square roots (Prop SQ).** Square roots in Cl(D) yield factors, given near-uniform class sampling (GRH, via Bach's
  generators). Known; Hühnlein–Meyer–Takagi 1998.
- **Cayley-graph walks.**
  - Every non-formal closed walk is a smooth value x² + kNy² ≥ kN (Lemma GIRTH, D < 0). So return and collision
    statistics are relation statistics.
  - For independent fixed-time return bits, the hidden genus term needs ≳ |H| walks (Prop WALK, a counterfactual
    bound). Its eigenvalue is the mean of (ℓ/p) over the generators, a QR-type quantity.
  - Endpoint *coefficients* are outside these statements.
- **Infrastructure.**
  - Any certified multiple of R (hR suffices; h need not be separated) leads in polynomial time to the midpoint
    candidate at R/2 (Prop HR, N(ε) = +1). That candidate gives a factor when it is useful. Usefulness is an
    observation: 72% at k = 1; 97% for some k ≤ 11 (E4).
  - An approximation of hR to precision Δ is finished by an O(√Δ) search. With Euler products (GRH; Shanks,
    Lenstra, Schoof) or random summation (unconditional, randomized; Srinivasan 1998), the specified architecture
    costs |D|^{1/5}. That is an upper bound for this architecture, not a floor.
  - SQUFOF's proper squares have density ≍ D^{−1/4} per cycle position. Giant steps to random positions do no
    better (E2).
- **CM / volcanoes.** Every known way to construct a curve with Frobenius discriminant −4kN·f² costs |D|^{1/2+o(1)}
  or more. We prove no lower bound for that construction.
- **Ring class groups and multiplier pencils.** Under the operations we examined, they add more hidden-order groups.
- **Barrier (Theorem GB, generic model).** In a generic hidden-order group whose odd part has a P-component with P
  uniform-like in [X, 2X], a q-query algorithm reaches hidden 2-torsion with probability O(q³/X).
  - The shape is the same as Sutherland's 2007 unknown-order bound and Hhan 2024 Thm 4.1. We claim no novelty.
  - Any generic success outputs a nonzero relation.
  - Applying the model to real class groups needs distributional hypotheses that we do not prove (§6.1).
- **What escapes the model (§6.2).** Each known escape we found uses one of these features:
  - smooth norms (relations: Seysen, Hafner–McCurley, Lenstra–Pomerance, SQUFOF2), i.e. (M);
  - order smoothness (SPAR, Schnorr–Lenstra), i.e. (M);
  - square norms (SQUFOF, N^{1/4});
  - collisions (rho in Cl, N^{1/4});
  - analytic order approximation (N^{1/5});
  - rational approximation / lattices / small roots (Lehman, Hart, McKee, Hittmeir, Harvey–Hittmeir, Coppersmith-type
    attacks), i.e. (A).

  This is a literature survey, not a proof of exhaustiveness. The open class is non-multiplicative, non-lattice uses
  of the integers of reduced forms.
- **Structural remark (OBSERVATION, not a reduction).** h(−4kN) is an Euler product of public local factors (D/ℓ).
  The hidden characters (ℓ/p) enter the group *structure*, not the local factors of h. In every method we found, the
  hidden 2-torsion is reached through the order (or a relation), not through the characters. This does not show that
  h is efficiently available, nor that every use of h is relational.

Nothing here is claimed as novel. Labels: PROVED (proof here), KNOWN (reference), OBSERVATION, CONJECTURED,
HEURISTIC.

## 1. Background facts used (KNOWN)

- **(F1) Ambiguous reduced forms** (elementary; Buell, *Binary Quadratic Forms* (1989) Ch. 3; Cox, *Primes of the
  form x² + ny²*, §3).
  - A reduced form f = (a, b, c) has [f]² = 1 iff f ~ (a, −b, c). By uniqueness of reduced forms this happens iff
    b = 0, a = b or a = c.
  - So the reduced form of each class in Cl(D)[2] is ambiguous; there is exactly one per class.
  - The reduced ambiguous forms are enumerable from the divisors of D (tested: `genus.ambiguous_reduced_forms` against
    PARI's 2-rank).
- **(F2) Gauss genus theory** (Cox, *Primes of the form x² + ny²*, 1st ed. 1989: Thm 3.15, Thm 6.1 — theorem
  numbers UNVERIFIED across editions; Gauss DA §§229–287).
  - Cl/Cl² ≅ (Z/2)^{μ−1}, cut out by the assigned characters.
  - For D = −4kN these are (·/p), (·/q), (·/ℓ) for odd ℓ | k, and the 2-adic characters among χ₋₄, χ₈, χ₋₈ (which
    ones depends on kN mod 8).
  - The principal genus equals Cl² (Gauss's duplication theorem).
- **(F3) The class number formula.**
  - h(D) = (w/2π)√|D|·L(1, χ_D) for D < 0. This holds for non-fundamental D too, with χ_D the Kronecker symbol.
  - For D > 0: hR = ½√D·L(1, χ_D).
- **(F4) Generation** (Bach, Math. Comp. 55 (1990) 355–380; GRH). The classes of primes ≤ 6 log²|D| generate Cl(D)
  for fundamental D.

## 2. Route 1: 2-Sylow access

### 2.1 Lemma 0: what an ambiguous form gives (PROVED; elementary; unit-tested)
Let N > 12k. Every reduced ambiguous form of D = −4kN is one of two kinds:
- it has first coefficient a | 4k (a *public* form; its class lies in A_pub := the subgroup generated by the classes
  of the primitive public ambiguous forms, e.g. (d, 0, kN/d) for d | k, together with (2, 2, ·) when primitive); or
- one of a, b, 2a ± b has gcd with N equal to p or q.

The class of (p, 0, kq) is not in A_pub. So A_pub ≠ Cl[2] and at least half of Cl[2] is useful.

(Correction after the proofs critic. v1 defined A_pub through "prime ideals above 2k". In the non-maximal order
Z[√−kN] (kN ≡ 3 mod 4) the prime above 2 is not invertible and (2d, 2d, c) is imprimitive, so the definition is now
by primitive public forms. The critic checked 4 493 triples (p, q, k), k ≤ 105, for both maximal and non-maximal
orders. In every case the public set is a subgroup, the class of (p, 0, kq) is never in it, and useful/|Cl[2]| = ½
exactly.)

*Proof.*
1. **Forms (a, 0, c).** Here ac = kN and a ≤ c, so a ≤ √(kN) < N. Hence gcd(a, N) ∈ {1, p, q}, and gcd = 1 means
   a | k.
2. **Forms (a, a, c).** Here a(4c − a) = 4kN and 4c − a ≥ 3a, so a ≤ √(4kN/3) < N. Again gcd(a, N) ∈ {1, p, q}, and
   gcd = 1 means a | 4k.
3. **Forms (a, b, a), 0 < b < a.** Put d₁ = 2a − b and d₂ = 2a + b. Then d₁d₂ = 4kN and d₁ ≤ d₂ ≤ 3d₁.
   - If gcd(d₁, N) = 1, then d₁ | 4k, so d₂ ≥ N > 12k ≥ 3d₁, which is impossible.
   - If N | d₁, then d₁d₂ ≥ N² > 4kN, also impossible.
   - So gcd(d₁, N) ∈ {p, q}.
4. **[𝔭] ∉ A_pub.**
   - *Non-maximal case* (kN ≡ 3 mod 4, so k is odd and the conductor is 2).
     - The primitive public reduced ambiguous forms are exactly (d, 0, kN/d) with d | k:
       - (2, 2, (kN+1)/2) and (2d, 2d, (kN/d + d)/2) have even c, so they are imprimitive;
       - a = 4d gives a non-integer c;
       - odd a with b = a forces 4 | a.
     - The ideals 𝔡 = (d, √−kN) for odd d | k are invertible, since odd norm is prime to the conductor, and they
       satisfy 𝔡² = (d).
     - So A_pub = {[𝔡] : d | k}. Then [𝔭] = [𝔡] makes the primitive invertible ideal 𝔭𝔡 principal, of norm pd. That
       gives x² + kNy² = pd, which is impossible as below with 2k replaced by k.
   - *Maximal case.* A relation [𝔭] = ∏[𝔩_i] over distinct ramified 𝔩_i above primes dividing 2k (these classes
   have order ≤ 2) would make the primitive ideal 𝔭∏𝔩_i principal, with norm p·d, d | 2k. Its generator x + y√−kN
   satisfies x² + kNy² = pd.
   - If y ≠ 0 this forces kN ≤ 2kp, i.e. q ≤ 2.
   - If y = 0 it forces x² = pd with p ∤ d.

   Both are impossible.
5. **Classes outside A_pub give p** (both cases). By F1 the reduced form of a class outside A_pub is ambiguous. If it were public
   (a | 4k), the class would be that of an ideal of norm a, a product of ramified primes above 2k, hence in A_pub.

∎ (Checked for k ∈ {1, 2, 3, 5, 6, 7, 10, 11} on six semiprimes; the count of reduced ambiguous forms matches
2^{2-rank} from PARI.)

### 2.2 Prop SQ: square roots in Cl(D) ⇔ factoring (KNOWN; precise version)
- **FACT ≤ SQRT** (randomized; conditional on GRH through near-uniform sampling). First trial-divide N by the primes
  ≤ 6 log²|D| (a hit is a factor).
  - Take D = −4kN fundamental: k = 1 if N ≡ 1 (mod 4), k = 2 otherwise. Let S be all primes ≤ 6 log²|D| not
    dividing N: the split ones and the ramified ones above 2k. The ramified classes are public, and Bach's theorem
    (F4) covers all prime ideals, not only split ones (critic fix).
  - Sample x = ∏ f_ℓ^{e_ℓ} with e_ℓ uniform in [0, N²). Then x is within statistical distance |S|·h/N² of uniform on
    ⟨S⟩ = Cl (F4).
  - Give y = x² to an arbitrary root oracle, which returns some z with z² = y. Then t = z·x⁻¹ ∈ Cl[2].
  - Given y, x is uniform on a coset of Cl[2], and z depends on y only. So t is uniform on Cl[2], and
    Pr[t ∉ A_pub] ≥ ½ − |S|h/N².
  - Lemma 0 then gives p.
  - Experiment E3 used a deterministic oracle, blind to x: the lexicographically least root. Success rates were
    0.485 and 0.4875 per query (n = 24, 32; predicted ½), with 389/389 correct factors.
- **SQRT ≤ FACT.** With the full factorisation of D known (p, q and the factored k), membership in Cl² is decided by genus characters. A square root is then
  computed in randomized polynomial time by Gauss's duplication/ternary-form method:
  - Hühnlein–Meyer–Takagi, ICISC 1998, §3.1 (algorithm) and §4 (factoring reduction) — VERIFIED by codex_1;
  - Lagarias, J. Algorithms 1 (1980) 142–186, Thm 6.10 (complexity) — theorem text UNVERIFIED.
- The equivalence is folklore in imaginary-quadratic cryptography (Buchmann–Hamdy survey, 2001; Biehl–Buchmann–
  Hamdy–Meyer, DCC 25 (2002)) — UNVERIFIED wording. **So "take a square root to land in Cl[2]" is factoring, not a
  route.**

### 2.3 Prop GEN: genus characters ≡ quadratic residuosity mod N (PROVED; randomized reductions; one direction uses
GRH for primes in progressions)
Let GEN-SQ be: given a form f of discriminant D = −4kN, decide whether [f] ∈ Cl² (the principal genus, F2). Let
QR_N be: given a with (a/N) = 1, decide whether a is a square mod N.

- **GEN-SQ ≤ QR_N.**
  - Evaluate f at random small (x, y) to get a represented m coprime to 2kN (a nontrivial gcd gives a factor).
  - The public assigned characters of m (2-adic characters and (m/ℓ) for ℓ | k) and (m/N) are computable.
  - [f] ∈ Cl² iff all of these are 1 and QR_N(m) holds, because (m/p) = (m/q) = 1 ⇔ (m/N) = 1 ∧ (m/p) = 1.
- **QR_N ≤ GEN-SQ** (`genus.qr_via_principal_genus`; unit-tested against a brute-force "is a square in Cl" oracle).
  1. Given a, use CRT to pick a prime m ≡ a (mod N) with m ≡ 1 (mod 8k), m < 8kN³ (the code draws m = m₀ + 8kN·s,
     s < N²).
  2. Then (D/m) = (−1/m)(k/m)(N/m) = 1·1·(m/N) = (a/N) = 1, by quadratic and Jacobi reciprocity using m ≡ 1 (mod 8)
     and m ≡ 1 (mod ℓ) for ℓ | k. So the prime form f_m exists and is computable (Tonelli–Shanks mod the prime m).
  3. All public characters of m equal their value at 1, namely 1. The hidden ones are (m/p) = (a/p) and
     (m/q) = (a/q) = (a/p).
  4. Hence [f_m] ∈ Cl² ⇔ (a/p) = 1 ⇔ a is a square mod N.
  5. Under GRH a random m in that progression below 8kN³ is prime with probability ≍ 1/log N, since
     8kN³ ≫ (8kN)²log⁴N.

∎

- **Hypotheses (codex_2).** k is small and given factored, which is needed to evaluate the public characters. The
  GRH range 8kN³ ≫ (8kN)²log⁴N needs k = N^{o(1)}; otherwise enlarge the range.
- **2-adic assigned characters**, with M = kN squarefree (codex_2 table, checked on 627 queries over all six
  admissible residues mod 8):

  | M mod 8 | 1, 5 | 3, 7 | 2 | 6 |
  |---|---|---|---|---|
  | extra character | χ₋₄ | none | χ₋₈ | χ₈ |

  All of them equal 1 at m ≡ 1 (mod 8).
- **Scope.** The equivalence is for *membership* (GEN-SQ). Evaluating a *labelled* character χ_p rather than χ_q on
  Jacobi-(−1) inputs needs an extra anchor choice. That is a stronger specification, not covered here.
- **Prior art.** The direction GEN-SQ ≤ QR_N is immediate from genus theory. Castagnos–Laguillaumie–Tucker (ASIACRYPT
  2022, ePrint 2022/1143) state that QR on Jacobi-(+1) elements is the same as deciding the principal genus for an
  RSA-type discriminant; the literature critic saw this only in a snippet. We regard Prop GEN as folklore, made
  explicit.
- **Consequence.** "Genus characters are hidden" is exactly the QR_N assumption. It is not circular in the strict
  sense: QR_N is not known to be factoring-hard (Adleman–McCurley 1994, open problem O11; see
  `notes/codex_r4_sos_classgroup.md` §1).
- But no subexponential non-relational algorithm for QR_N is known either.
- Even a GEN-SQ oracle is not known to yield Cl[2] elements. Testing squareness does not produce a root.
- Complete genus information does not even locate the involutions (codex_1). Cl(−56) ≅ C₄, and its nonidentity
  involution [(2, 0, 7)] = [(3, 2, 5)]² is a square, so every genus character is trivial on it.
- Codex_1 independently proved the same QR_N reduction for D = −4N, using primes r ≡ a (mod N), r ≡ 1 (mod 4).

### 2.4 Lemma GIRTH and Prop WALK: Cayley-graph statistics are relation statistics (PROVED)
Let S be a set of split primes ℓ ∤ D with prime forms f_ℓ (ideals 𝔩, conjugates 𝔩̄ = f_ℓ⁻¹). Let H = ⟨S⟩ and
Λ_S = {e ∈ Z^S : ∏ f_ℓ^{e_ℓ} = 1}.

**Lemma GIRTH** (D < 0 only; for D > 0 the norm bound fails through Pell-type cancellation). "Girth" means girth
after the formal (commutation and f·f⁻¹) relations are cancelled. The plain Cayley graph has 4-cycles.
- (i) For e ≠ 0 in Λ_S, ∏ℓ^{|e_ℓ|} = x² + kNy² with gcd(x, y) = 1 and y ≠ 0. So ∏ℓ^{|e_ℓ|} ≥ kN and
  ‖e‖₁ ≥ g₀ := ⌈log(kN)/log ℓ_max⌉.
- (ii) Conversely, every primitive S-smooth value of x² + kNy² with y ≠ 0 gives an element of Λ_S.

*Proof.*
- 𝔞_e = ∏_{e>0}𝔩^{e}∏_{e<0}𝔩̄^{|e|} is primitive (for each ℓ only one of 𝔩, 𝔩̄ occurs) and has norm ∏ℓ^{|e|}.
- [𝔞_e] = 1 iff 𝔞_e = (α) with α = x + y√−kN ∈ Z[√−kN].
- If y = 0, then (x) is primitive only for |x| = 1, i.e. e = 0.
- For (ii): (x + y√−kN) is primitive of S-smooth norm, so it factors into 𝔩's and 𝔩̄'s.

∎

So the closed-walk count W_t = Σ_{λ∈Λ_S} #{walks in Z^S of length t ending at λ} is a function of Λ_S ∩ {‖·‖₁ ≤ t}.
In particular W_t = W_t^free for t < g₀. **Every return or collision observed by a walk statistic is a smooth value
of the principal form. That is the Hafner–McCurley/Seysen relation source, i.e. family (M).** E1 confirms
W_t = W_t^free below g₀ in 60/60 instances. The first excess appears at g₀ + {0..4}.

**Prop WALK (return-frequency statistics; PROVED).**
- The spectral expansion is P_t = W_t/(2|S|)^t = (1/|H|)Σ_{ψ∈Ĥ} λ_ψ^t, with λ_ψ = |S|⁻¹Σ_ℓ Re ψ(f_ℓ).
- The hidden genus characters restrict to χ_p(f_ℓ) = (ℓ/p), so λ_{χ_p} = |S|⁻¹Σ_{ℓ∈S}(ℓ/p). This is a QR-type
  aggregate.
- Let G_t be the contribution of the hidden characters that are distinct on H from the public ones (deviation D1,
  §7). There are at most c_μ := 2^{μ−2} of them (c_μ = 2 for k = 1).
- **Hypothesis.** g₀ > 2 (no relation of ℓ₁-length ≤ 2), so that P_t ≤ ½ for t ≥ 2. Without it the bound fails.
  Codex_2's example is N = 15, S = {17}, P₂ = 1. E1 satisfies the hypothesis (g₀ ≥ 4).
- For even t, every term is ≥ 0, so b := P_t − G_t ≥ 1/|H| and G_t ≤ c_μ/|H|.
- The squared Hellinger distance between Bernoulli(P_t) and Bernoulli(b) is
  ½[(P_t − b)²/(√P_t + √b)² + (P_t − b)²/(√(1−P_t) + √(1−b))²]. For t ≥ 2 we have P_t ≤ ½, and this is
  ≤ (P_t − b)²/(P_t + b) ≤ c_μ²/(2|H|). (v1 omitted the second term. Critic fix; the bound is unchanged up to the
  constant.)
- M independent walks therefore distinguish "with" from "without" the hidden term with advantage ≤ c_μ√(M/|H|).
  Advantage δ needs M ≥ δ²|H|/c_μ² ≈ δ²·N^{1/2}/c_μ².
- **Scope.** Each walk contributes one Bernoulli at a single fixed t. Whole trajectories, several t, adaptive
  strategies, and the distribution of reduced-form coefficients at walk endpoints are **not** covered. The "without"
  hypothesis is counterfactual: deleting spectral terms does not produce another valid class-group instance. So
  this bounds an artificial testing problem, not factoring.
- The general statement is GIRTH: *equality events* are relations. Endpoint coefficients are non-generic features
  (§6.2), and they can leak at toy size. Codex_2's example is N = 21, S = {17}: the one-step endpoint (5, 4, 5) is
  itself a hidden ambiguous form, and gcd(2a − b, N) = 3. For large N that happens with probability ≍ |Cl[2]|/h.
- For collision counts (pairs of walks), the signal-to-noise heuristic gives M ≳ √(|H|/2). This is the birthday
  bound, consistent with Theorem GB.
- The direct collision method (codex_1). For uniform x, y, Pr[x² = y², xy⁻¹ ∉ A] = (|Cl[2]| − |A|)/h. So birthday
  collisions among squares reach hidden 2-torsion in Θ(√(h/|Cl[2]|)) = N^{1/4+o(1)} samples. This is a generic
  N^{1/4} route, not a subexponential one.
- E1, for the variance criterion M = P_t/G_t²: the minimum over t of M/|H| was ≥ 1.0 (as-registered) and ≥ 2.0
  (D1) in every instance. Medians were 1.4·10³, 4.1·10³, 1.3·10⁴ at n = 20, 24, 28. The predicted upper bound of 100
  (P2) was **not met**: in the median instance, detection is far *harder* than the proved floor. (The critic notes
  that the provable floor is ½ for the registered G_t and ¾ for D1, so the registered ¼ was loose.)

### 2.5 4-rank via Rédei
- Rédei matrices are indexed by the prime-discriminant factors (·/p) and (·/q). Building them needs the factors.
- Without the factors, the 4-rank still decides a fixed-element QR bit (codex_r4 F2). For N ≡ 1 (mod 8):
  [2 is a square mod N] ⇔ r₄(Cl(−8N)) > 0 if p ≡ q ≡ 1 (mod 4), and ⇔ r₄(Cl(−4N)) > 0 if p ≡ q ≡ 3 (mod 4).
  - We re-verified this independently: 92/92 random pairs p, q < 400, using PARI `quadclassunit`
    (`scratchpad/genusleak/r4check.py`).
- So computing r₄ is at least as hard as that QR instance. It is also not a route to Cl[2]: r₄ describes Cl[4]/Cl[2],
  and Cl[2]'s *size* is already public.

### 2.6 Route-1 verdict
Every way into Cl[2] that we found is one of three things:
- a root (factoring-equivalent);
- a group-order computation: a multiple of the exponent kills the odd part, then 2-Sylow linear algebra (Shanks;
  Schnorr–Lenstra; see the §3.3 correction);
- a relation. Walk returns and collisions are relations by GIRTH.

Genus characters are QR_N-hard and do not by themselves produce Cl[2] elements.

## 3. Route 2: infrastructure (D = 4kN)

### 3.1 Midpoint (KNOWN; E4)
- When N(ε) = +1, the principal cycle is symmetric under conjugation. The reduced principal ideal at distance R/2 is
  ambiguous, with a | 4kN, and a ∤ 4k gives p.
- In the continued-fraction picture this is the symmetry point of the even period of √(kN) (Legendre; Lehmer–Powers
  1931).
- E4: at k = 1 the midpoint splits N in 72% (n = 28) and 64% (n = 32) of instances; some k ≤ 11 splits 97% / 98%.

### 3.2 Prop HR: any multiple of R suffices (PROVED, given standard infrastructure arithmetic)
- **Statement** (a midpoint *search*; it yields a factor when the midpoint is useful). Assume N(ε) = +1. This is automatic if a prime ≡ 3 (mod 4) divides kN. If N(ε) = −1 the period is
  odd, and the symmetry point is a form (a, b, −a) with 4kN = b² + 4a², not a divisor (critic fix). Suppose X = M·R
  is known to absolute precision 0.1, for an unknown integer M ≥ 1. Then the ambiguous ideal at R/2 is reached in
  poly(log X, log D) operations (the neighbour check uses the bounded gap
  between consecutive reduced principal ideals near a given distance).
- **Proof.**
  - For j = 1, 2, …, ⌊log₂ X⌋ + 1, compute the reduced principal ideal nearest distance X/2^j. Shanks–Lenstra
    giant steps track distances to O(log X) bits: Lenstra 1982; Buchmann–Williams; Jacobson–Williams, *Solving the
    Pell Equation* (2009), Ch. 10–11.
  - Inspect its O(1) neighbours for an ambiguous form.
  - At j = v₂(M) + 1 the distance X/2^j = (odd)·R/2 ≡ R/2 (mod R).

  ∎
- **Consequence.** Separating h from R is unnecessary: hR itself (F3) is such an X. That the midpoint ideal splits
  N with constant probability per k is an **OBSERVATION** (E4), not proved.
  - It can fail. Codex_1's example is N = 731 = 17·43, with ε = 730 + 27√731 of norm +1. The midpoint reflects the
    public factor 2: N(27 + √731) = −2.
  - Mind the conventions: h⁺R⁺ = 2hR.
  - An approximation to a large multiple is not by itself a period-membership oracle. Certified precision, e.g.
    (f, p)-representations, is assumed.
- This is the classical "regulator ⇒ factoring" reduction (Lenstra 1982; Schoof 1982; used by Hallgren, JACM 54
  (2007)). Murru–Salvatori (arXiv 2409.03486, §5) explicitly discuss regulator multiples and repeated halving (codex_2,
  VERIFIED there). No novelty claimed.

### 3.3 What a precision-Δ approximation of hR buys (PROVED as an upper bound; KNOWN in substance)
- Given |A − hR| ≤ Δ, find a distance in [A − Δ, A + Δ] at which the principal ideal recurs. Any such distance is a
  multiple of R, which is enough by Prop HR.
- BSGS (baby steps of total distance w around 0, giant steps of length w) finds one in O(w + min(Δ, R)/w) steps,
  i.e. O(√min(Δ, R)).
- **Imaginary case (corrected; v1 was WRONG).** v1 said: "BSGS for an element's order, take an odd multiple, square
  one random class down". The critic's counterexample is N = 2311·2671, D = −4N, h = 1280 = 2⁸·5. Squaring down
  g^{h_odd} for a single class g gives a factor for only 1/128 of all g, because the 2-Sylow has a long cyclic factor
  whose involution is public. Also, BSGS returns a multiple of ord(g), not of h. The correct finish:
  1. Run BSGS in the window for a few random g, and take the lcm of the annihilators. This gives a multiple M of the
     exponent with high probability.
  2. Sample the 2-Sylow uniformly as g^{M_odd}.
  3. Compute Cl[2] from a basis of the 2-Sylow by 2-adic linear algebra / Pohlig–Hellman. This costs
     poly(2^{μ}, log h), since the 2-rank is ≤ μ − 1.

  The window cost √Δ is unchanged.
- **Unconditional randomized estimator** (codex_1; reasoning checked by us). This is in substance Srinivasan's
  "Random Summation Technique" (Math. Comp. 67 (1998) 1285–1308), so credit is hers.
  - The tail is bounded by |L(1,χ) − Σ_{n≤T}χ(n)/n| ≤ 2|D|/T.
  - Choose a dyadic block of [1, T] uniformly, then n uniformly in it. The estimator Z = χ(n)/(n·p_n) is unbiased
    and satisfies |Z| ≤ K = O(log T).
  - Hoeffding then gives absolute precision ε from Õ(ε⁻²) samples. So precision Δ for h (or hR) costs Õ(|D|/Δ²), and
    the total Õ(|D|/Δ² + √Δ) gives |D|^{1/5} with no GRH. The result is randomized and uncertified; the BSGS finish
    verifies it.
  - The imaginary finish also needs a sampling hypothesis (GRH, §2.2). A verified annihilator certifies only itself,
    not the analytic estimate.
- **Euler products.**
  - Under GRH the truncated product to x has |log L(1,χ) − log L_x| ≪ x^{−1/2}log²(x|D|) (Bach, "Improved
    approximations for Euler products", Number Theory (Halifax 1994), CMS Conf. Proc. 15 (1995) 13–28 — exact constant
    UNVERIFIED). So Δ(x) ≈ √|D|·x^{−1/2}.
  - The total cost x + |D|^{1/4}x^{−1/4} is minimised at x = |D|^{1/5}. That is the Lenstra/Schoof exponent (also
    `notes/analytic_results.md`, `notes/below16_results.md` §5). The deterministic GRH N^{1/5} class-group algorithm is
    credited to Shanks 1971 (Proc. Symp. Pure Math. 20, 415–440) by Harvey, Math. Comp. 90 (2021) §1.
  - E5 measures the fitted slope d log Δ/d log x ∈ {−0.52, −0.53, −0.56} and log C_min/log|D| ∈ {0.20, 0.18, 0.18}
    at n = 32, 40, 48.
- **Exponent 1/5 is the optimised cost of this architecture, not a proved floor** (codex_2). Better estimators,
  arithmetic cancellation, or other finishes are not excluded.
  - A subexponential algorithm along this route needs Δ = |D|^{o(1)} in subexponential time, i.e. relative precision
    |D|^{−1/2+o(1)} for L(1, χ).
  - Euler products need x ≈ |D|^{1−o(1)} for that: tail variance ≍ 1/(x log x) (below16 §5.3; heuristic).
  - We know no other non-relational subexponential approximation of L(1, χ_D) to that precision.

### 3.4 Square forms and giant steps (OBSERVATION + E2)
- A proper square Q_i = r² in the principal cycle with a *primitive* root form is 𝔯², with N𝔯 = r ≲ D^{1/4} and
  [𝔯] ∈ Cl⁺[2] ∖ {1}.
  - Correction (codex_2). When r is even, the root form can be imprimitive. For N = 21, Q₂ = 4 has root form
    (2, 2, −10). It belongs to the maximal order of discriminant N, not to Cl⁺(4N).
  - 36/90 E2 instances contain a proper square with even r (10, 12, 14 at n = 28, 32, 36). The factor computations
    stand; the class-group reading applies to discriminant D/4 there.
- So SQUFOF finds 2-torsion through a *multiplicative self-coincidence* of a reduced norm (a = r·r). The density of
  such points per cycle position is ≍ |Cl⁺[2]|·D^{−1/4} (heuristic; Gower–Wagstaff 2008).
- Giant steps reach chosen distances. The distance of 𝔯² is 2·dist(𝔯), where dist(𝔯) is measured in the cycle of a
  hidden class, so giant steps cannot aim at squares.
- In the generic infrastructure model below, jumping to random distances costs 1/ρ forms, exactly like walking.
- **E2.**
  - Median i₁·ρ (first proper square × density): 0.45, 0.54, 1.00 at n = 28, 32, 36. At n = 28 one instance has no
    proper square in its period, so that median is over 29 instances. P4 (∈ [0.5, 2]) is **met at
    32 and 36, not met at 28**: sequential walking is, if anything, slightly better than random jumps.
  - Fitted slope of log₂ median i₁: 0.38 (P5, ∈ [0.18, 0.32], **not met**).
  - Post hoc: the slope of log₂ median(1/ρ) is 0.24, consistent with N^{1/4}. The excess in i₁ at small n comes from
    i₁·ρ drifting upward with n.
- CFRAC is the relational completion of SQUFOF: products of Q_i that form squares.

### 3.5 Are Buchmann/Vollmer subexponential regulator algorithms relational? (KNOWN; references UNVERIFIED in detail)
Yes. They all collect relations among prime ideals of smooth norm, then do lattice reduction / HNF (Hafner–McCurley
style):
- Buchmann, "A subexponential algorithm for the determination of class groups and regulators of algebraic number
  fields", Sém. Théorie des Nombres Paris 1988–89 (1990);
- Vollmer, ANTS-IV (LNCS 1838, 2000) and ANTS-V (LNCS 2369, 2002);
- Biasse–Jacobson, ANTS-IX (2010).

Smooth norms put them in (M).

### 3.6 Prop GI: generic infrastructure barrier (HEURISTIC ABSTRACTION / proof proposal; downgraded after codex_2)

The oracle semantics (reduction corrections, finite-precision distances) and the transcript simulation are not
specified. Only the fixed-difference measure computation below is proved.
- **Model.** The algorithm holds handles of reduced principal ideals with known distances, tracked exactly (e.g.
  via (f, p)-representations), so errors do not accumulate. Each handle sits within c of its nominal distance. The
  algorithm may compose handles (distances add) and compare them. The regulator R is hidden, with density ≤ 1/X on
  [X, 2X], or ≤ 1/Δ on an analytic window of width Δ.
- **Collisions.** A collision means some known difference d = δ₁ − δ₂ with |d| > c is within c of mR, for some
  integer m ≥ 1. Since R ∈ [X, 2X], m lies in [(|d| − c)/2X, (|d| + c)/X], a range of ratio about 2. The set of such
  R has measure ≤ Σ_m 2c/m = O(c), independent of |d|.
- **Bound (conditional on a transcript simulation that is not supplied).** As in Theorem GB, Pr[success] ≤ O((q + 2)²·c/X), so q = Ω(√(X/c)). BSGS achieves this. (v1 had a
  spurious factor (1 + L ln 2) and concluded Ω(X^{1/3}); the critic's sharpening. The q³ in GB is genuine because P
  is discrete; here R is continuous.)
- **Reading.** Without relations and without square coincidences, the infrastructure is a hidden-period problem. R
  is the order of R/RZ.

## 4. Route 3: isogeny volcanoes / CM (failure points)

### 4.1 Construction
- An ordinary E/F_ℓ with Z[π] of discriminant t² − 4ℓ = −4kN f² needs ℓ ≥ kNf² (4ℓ = t² + 4kNf²).
- The curves with trace t number H(4ℓ − t²) ≈ f√(kN)·polylog (Deuring; Kronecker–Hurwitz class number), out of
  ≈ 2ℓ curves in total. So random curves plus SEA need ≈ 2ℓ/H ≥ √(kN)/polylog trials, which is worse than trial
  division.
- Fixing ℓ does not enlarge the target. At most 2 traces (±t) satisfy (t/2)² + kNf² = ℓ.
- In H(4kNf²) ≈ f√(kN)·polylog, the polylog holds under GRH; unconditionally it is |D|^{o(1)}.
- The CM method needs the class polynomial H_O mod ℓ.
  - deg H_O = h(O) ≥ h(O_K)/3 = |d_K|^{1/2+o(1)} by Siegel (ineffective), with K = Q(√−kN).
  - Computing it costs Õ(|D|): Enge, Math. Comp. 78 (2009); Sutherland, Math. Comp. 80 (2011).
  - Sutherland, "Accelerating the CM method", LMS J. Comput. Math. 15 (2012), finds one root without writing out
    H_O (codex_1, VERIFIED there). It still costs |D|^{1+o(1)} under GRH and uses class-group enumeration at
    auxiliary primes.
- We know of no method that produces any curve whose endomorphism algebra is K in time |d_K|^{o(1)}. Producing one
  is the well-known "CM bottleneck"; see e.g. Bröker–Stevenhagen, "Constructing elliptic curves of prime order", Contemp. Math. 463 (2008) 17–28, where |D|
  must be small.
- **Failure point:** every *known* construction costs |D|^{1/2+o(1)} or more. A large class-polynomial degree does
  not force an algorithm producing one curve to compute that polynomial, so this is not a lower bound (codex_2).

### 4.2 Even given the curve
- **Crater walks.** Small split ℓ′-isogenies realise the class-group *action*: a hidden-order group action, weaker
  than a black-box group. The crater length is ord([𝔩′]), an order computation.
- **The ambiguous classes.** On the crater (End ⊇ Z[√−kN], p ∤ f), [𝔭] acts by an isogeny of degree p with kernel
  E[p] ∩ ker(π − t/2), since (π − t/2)² = −kNf². Computing that representative needs p. The same class could act
  through another ideal representative, e.g. a chain of small isogenies. Finding one is the class-group problem
  again, so this is not a general obstruction by itself.
- **Choosing t ≡ 2 (mod 2N).** This makes ker(√−kN) ∩ E[N] F_ℓ-rational: E(F_ℓ) contains Z/N. But a point of order N
  does not separate p from q.
- **Group structure.** With O = End_{F_ℓ}(E) exactly, E(F_ℓ) ≅ O/(π − 1) as O-modules (H. W. Lenstra Jr., J. Number Theory 56 (1996)). Its order
  #E = (t/2 − 1)² + kN f² is public.
- **No leak found.**

## 5. Route 4: non-maximal orders and multiplier pencils

### 5.1 Ring class groups
- 1 → K_m → Cl(Dm²) → Cl(D) → 1, with |K_m| = m∏_{ℓ|m}(1 − (D/ℓ)/ℓ)/[O*: O_m*] (Cox Thm 7.24). The order of K_m is
  public.
- A known multiple of |K_m| kills K_m, not the hidden odd part of Cl(D).
- In Theorem GB the Z/P-component survives unless P | m∏(…), i.e. unless the algorithm "guessed" a multiple of P.
  That is an L-bit integer divisible by P, which has probability ≤ Lρ/log₂X.
- If K_m is odd and x already maps to an involution, then x^{K_m} is an involution upstairs with the same image
  (codex_1). This lifts a known target; it does not discover one.
- No gain *for the powering operation examined*. Explicit maps, transfer, norms, or coefficient operations are
  outside GB and are not excluded (codex_2).

### 5.2 Multiplier pencil
- The hidden data common to all Cl(−4kN) is the character (·/p), not a group element. The orders h(−4kN) are
  *public* Euler products in k (F3), and they vary independently enough to be separate lottery tickets: Schnorr–
  Lenstra's multiplier trick, i.e. (M).
- **The compositum L = Q(√−kN, √−k′N) ⊃ Q(√kk′).**
  - Since p ∤ kk′, p is unramified in Q(√kk′). So 𝔭_k O_L = 𝔭_{k′}O_L (both are the product of the primes of L above
    p), and [𝔭_k], [𝔭_{k′}] extend to the same class.
  - Finding such a pair is a cross-group collision in Cl(L). There, equality testing is a principal-ideal problem in a
    quartic field (subexponential only relationally), over a group of order ≈ h_k h_{k′}h(kk′) (Kuroda's class
    number formula; constant UNVERIFIED).
- Many discriminants give more (M)-type chances, not a non-smooth genus character, as far as we found. "Vary
  independently enough" is a heuristic. Applying GB across correlated groups would need a joint oracle model, which
  we do not supply.
- Extension to L can principalise classes. A larger class number alone does not prove a harder collision problem.

## 6. Route 5: the barrier, and the non-generic escapes

### 6.1 Theorem GB (generic hidden-order barrier; PROVED in the model)
**Model.**
- G = G′ ⊕ Z/P is finite abelian. P is an odd prime drawn from a law π with min supp π ≥ X and max point mass ρ.
  P ∤ |G′|.
- A ≤ G′[2] is public, given with a basis.
- Generators g_i = (g_i′, u_i), i ≤ s. The u_i are uniform and independent in Z/P, and g_i′ is arbitrary — the
  algorithm may even know G′ completely.
- Random injective encoding σ: G → Σ. The oracle returns σ(x − y). The inputs are σ(0), σ(g_i), and σ(a_j) for
  j ≤ r.
- Success: output σ(t) with t ∈ G[2] ∖ A.

**Theorem GB.** Count sampling calls (if any) among the q calls; each sample gets a fresh formal variable. Then any
algorithm making q oracle calls succeeds with probability

  ≤ ((q + s + r + 1)² + 1)·(1/X + ρ(q + 1)/log₂X) + (q + 1)|G|/|Σ|.

For P uniform on the primes of [X, 2X] this is O(q³/X). So q = Ω(X^{1/3}).

*Proof (Shoup-style simulation; the order is unknown, so divisibility by P is bounded by counting).*
1. **Formal vectors.** Track, for every handle, a formal vector (e, c) ∈ Z^s × F₂^r. Inputs are unit vectors and
   oracle answers are differences, and each subtraction at most doubles the largest entry, so |e| ≤ 2^q. Here
   s, r = O(q).
2. **The simulation.** Simulate with fresh random strings, identifying two handles only when their formal vectors
   are equal. Since A has a known basis, e = e′ and c ≠ c′ really do give different elements.
3. **When simulation and reality differ.** They can diverge only if some pair with d = e − e′ ≠ 0 has d·u ≡ 0
   (mod P).
   - If d ≢ 0 (mod P), this has probability 1/P ≤ 1/X over u.
   - If d ≡ 0 (mod P), then P divides a nonzero integer of absolute value ≤ 2^{q+1}. That integer has
     ≤ (q + 1)/log₂X prime factors ≥ X, so the probability over P is ≤ ρ(q + 1)/log₂X.
   - In the simulation the pairs are independent of (u, P), so a union bound over pairs applies.
4. **Success in the simulation.** The output is either a fresh string (probability ≤ |G|/|Σ|) or a handle (e, c).
   For t ∈ G[2] ∖ A we need e ≠ 0 (otherwise t = c·a ∈ A) and e·u ≡ 0 (mod P), because P is odd and t has order
   ≤ 2. That has the same probability bound.

∎

**Corollary GB-X (extraction; any G, no distribution).**
- Suppose a generic algorithm outputs a known combination t = e·g + c·a ∈ G[2] ∖ A. Then e·π_odd(g) = 0 and
  e·g ≠ c·a, so 2e is a nonzero relation of G.
- For s = 1, e is a nonzero multiple of the odd order of g.
- Conversely, an odd multiple M′ of ord_odd(g) gives g^{M′} ∈ G₂, and squaring down reaches the involution of ⟨g⟩.
  That involution is hidden only if the involution of ⟨g⟩ lies outside A. Codex_2's counterexample: G = C₄ ⊕ C₂,
  A = ⟨(2, 0)⟩, g = (1, 1).

So in rank 1, *under that reachability hypothesis*, "hidden 2-torsion in ⟨g⟩" and "a multiple of ord_odd(g)" are
interchangeable. In general a relation is a necessary output certificate, not an equivalence with finding a useful
relation. In general rank, success
hands over a relation. **In this model, family (Q) is a group-order / relation problem.** Several multipliers k give
several groups; the union bound adds their success probabilities. A small P⁺(h_odd) for some k is exactly the
Schnorr–Lenstra smoothness event, i.e. (M).

**Sampling oracles.** Charged sampling calls do not break the theorem: each sample gets a fresh formal variable
(codex_1 Theorem C). Unlimited free samples could reveal |G|. The single-generator simulation a·g^u needs G = A⟨g⟩
(codex_2).

**Scope.**
- The bound O(q³/X) is weaker than the √X of BSGS. The q³ comes from exponents divisible by many primes in [X, 2X].
  - This is the same shape as Sutherland's Ω(P^{1/3}) unknown-order lower bound (PhD thesis, MIT 2007, §§2.2 and 4;
    VERIFIED by codex_1).
  - Sutherland's O(√(P/log log P)) primorial-steps algorithm shows that a literal Ω(√P) is false in the unrestricted
    model.
  - With all exponents < X, the q³ term disappears and the bound becomes Ω(√X).
  - Hhan, arXiv 2402.11269 (2024), §4, reproves unknown-prime-order bounds.
- Related work:
  - Damgård–Koprowski, EUROCRYPT 2002: root extraction; needs order-distribution hypotheses.
  - Rotem–Segev–Shahaf, EUROCRYPT 2020: delay functions; a different target.
  - The low-order assumption: Boneh–Bünz–Fisch, "A survey of two verifiable delay functions", ePrint 2018/712,
    Def. 1.
  - Belabas–Kleinjung–Sanso–Wesolowski, ePrint 2020/1310: odd low order at special discriminants, not 2-torsion.
  - Shoup 1997 and Maurer 2005 are known-order discrete-log bounds and do not apply directly.
  - Exact statements are checked in §10.
- **Application hypotheses (not proved for Cl(±4kN); codex_2):**
  - the transcript and side information are independent of (P, u);
  - small prime forms have independent uniform projections to the P-component;
  - the public D carries no extra information;
  - Cl has a direct P-component G′ ⊕ Z/P;
  - s, r = O(q).

  Real class groups are not drawn from π. The model is a Cohen–Lenstra-type idealisation, and it is broken by the
  non-generic escapes below.
- **Prior art.** GB is the same as Sutherland's thesis Thm 2.3 / Cor 2.4 and Hhan's Thm 4.1, plus the extraction
  step GB-X. Treat it as known in substance. Closest analogue: Seres–Burcsi, ePrint 2020/402 (low order ≡ factoring
  in RSA groups).

### 6.2 Which non-generic features escape, and to where (literature survey; not a proof of exhaustiveness)

Categories, defined operationally:
- **norm relations**: smooth norms factored over a base;
- **order smoothness**: exponent = a product of small primes;
- **collisions**: equality of group elements;
- **square coincidence**: a norm equals r²;
- **analytic order approximation**;
- **approximation / lattice / small-root reconstruction**: rational approximation, short vectors, Coppersmith-type.

| Feature or method | In the generic model? | Mechanism | Lands in |
|---|---|---|---|
| Canonical reduced forms (equality) | yes (encoding) | — | generic |
| Norm a of reduced forms: smoothness | no | norm relations: Seysen, Math. Comp. 48 (1987); Hafner–McCurley, JAMS 2 (1989); Lenstra–Pomerance, JAMS 5 (1992) (rigorous); Buchmann 1990; Vollmer 2000/02; Biasse–Jacobson 2010; SQUFOF2 (Bradford–Wagstaff, Bull. Pol. Acad. Sci. Math. 70 (2022), QS-style sieving builds the square form) | (M) |
| Exponent smoothness | partly | SPAR (Shanks–Pollard–Atkin–Rickert), Schnorr–Lenstra, Math. Comp. 43 (1984) | (M) |
| Collisions in Cl | yes | rho in Cl, BSGS, primorial steps (Sutherland 2007) | generic, N^{1/4} |
| Norm a: squareness | no | SQUFOF and racing multipliers (Gower–Wagstaff 2008) | (Q), N^{1/4} (heuristic) |
| gcd(a, N) on random forms | no | hits with probability ≍ 1/p (or ≍ \|Cl[2]\|/h for ambiguous endpoints) | accidental |
| Explicit forms of prescribed prime norm | no | sampling ideals; QR_N ≤ GEN-SQ (Prop GEN) | QR_N |
| Genus characters | no | ≡ QR_N (membership) | QR_N (open) |
| Gauss descent (square roots) | no | roots without h, *given* the factorisation of D | circular |
| Distances in the infrastructure | partly (GI, heuristic) | period R; multiple of R ⇒ midpoint (HR); Elia; Murru–Salvatori 2024 | order of R/RZ |
| Analytic class number formula | no | hR to precision Δ ⇒ √Δ search; \|D\|^{1/5} (Shanks; Lenstra; Schoof; Srinivasan) | order approximation |
| Public approximation of φ(N) | n/a (unit group) | Hittmeir; Harvey–Hittmeir N^{1/5} (deterministic): rational approximation + lattices | (A) |
| Representations of N by forms | no | Fermat, Lehman, Hart (N^{1/3}), McKee (Euler's two representations, N^{1/4+ε}) | (A) / collisions of representations |
| Coppersmith / lattices on form coefficients | no | Castagnos–Joux–Laguillaumie–Nguyen, ASIACRYPT 2009: factors pq² given a NICE kernel hint. Needs a hidden conductor q and auxiliary input; no squarefree-pq analogue known | (A), needs a promise |
| CM curves / isogenies | as a group action | known constructions cost \|D\|^{1/2+o(1)} or more | — |
| Class-field values (traces of class invariants generate the genus field ∋ √p*) | no | enumerates Cl², exponential | — |
| Conductors / multipliers / composita | yes (more groups), under the examined operations | more lottery tickets | (M) |

**Open class.** Algorithms that use the integer coefficients of reduced forms (or CF partial quotients, or
distances) in ways that are none of the above. That excludes multiplicative coincidences (smooth, square, equal),
approximation/lattice reconstruction, and order approximations. The Coppersmith-type pq² attacks show that
non-multiplicative coefficient arithmetic can leak, given a structural promise; we know of no such promise for
squarefree pq.
- A pure-integer statistic of the (a_i) along the principal cycle is an example. eps_results §5 found such statistics
  correlating with ε, but only through ((−1)^iQ_i | p) admissibility, i.e. the QR_N-type genus data of Prop GEN, not
  through Cl[2] elements.
- We have no barrier for this class. A barrier there would be a factoring lower bound.

## 7. Experiments (pre-registered; outcomes)

| P | Statement | Outcome |
|---|---|---|
| P1 | W_t = W_t^free for t < g₀ | **met** (60/60) |
| P2 | min M_ret/\|H\| ≥ ¼ (proved) and median ≤ 100 | lower bound **met** (min 1.0); upper **not met** (medians 1.4e3, 4.1e3, 1.3e4) |
| P3 | min M_coll/√\|H\| ≥ 1/√2 (proved) and median ≤ 50 | lower bound **met** (min 1.41); upper **not met** (53, 90, 160) |
| P4 | median i₁ρ ∈ [0.5, 2] | met at n = 32 (0.54) and 36 (1.00); **not met** at 28 (0.45) |
| P5 | slope of log₂ median i₁ ∈ [0.18, 0.32] | **not met** (0.38); post hoc, the slope of 1/ρ is 0.24 |
| P6 | sqrt-oracle success 0.5 ± 0.08 | **met** (0.485, 0.4875; H = Cl in 20/20) |
| P7 | claimed factors correct | **met** (389/389) |
| P8 | midpoint splits at k = 1 in [0.2, 0.8] | **met** (0.72, 0.64) |
| P9 | some k ≤ 11 splits ≥ 0.9 | **met** (0.97, 0.98) |
| P6′ | Euler-product error slope ∈ [−0.75, −0.25] | **met** (−0.52, −0.53, −0.56) |
| P7′ | log C_min / log\|D\| ∈ [0.15, 0.30] | **met** (0.20, 0.18, 0.18) |

**Deviation D1 (disclosed; wording corrected after the critic).** The registered G_t double-counts when a hidden
character coincides on S with a public one or with the other hidden one. There are two cases:
- χ_p|_S or χ_q|_S is trivial. This is the same event as the other one equalling χ₋₄|_S. It occurs in 4 instances.
- χ_p|_S = χ_q|_S, i.e. χ₋₄ is trivial on S (all ℓ ≡ 1 mod 4). This occurs at n = 28, seeds 0 and 42.

H = h/2 occurs in 6 instances (1, 1, 4 at n = 20, 24, 28): 4 of the first kind and 2 of the second. The code
deduplicates characters as functions on S. We report both versions. The proved lower bounds hold for both, and the D1
medians also miss the P2/P3 upper predictions, so no status changes. Instances with no distinct hidden character on
H: 1, 1, 2.

**Deviation log (codex_2 finding 16; no status changes).**
- **E3** rows omit the promised skipped-seed records (seeds failing N ≡ 1 mod 4). They are recoverable from
  `instances()`.
- **E2**, n = 28, seed 29: no proper square in the period. It is excluded from the medians; this failure rule was not
  pre-registered. Failure-inclusive (as +∞) median i₁ = 110, versus 98 excluded. P5 stays not met either way.
- **Floating point.** E1's derived sample sizes and E2's densities are computed in floating point from exact counts.
  E1 JSON contains `Infinity` for undefined ratios.
- **E4** also checks gcd(N, P_i). Codex_2 re-ran the registered Q_{L/2}-only rule and got identical outcomes.
- **E5.** C_min uses the true, harness-known error and the minimum over a retrospective grid. At n = 32, 17/20
  minima sit at the smallest grid point x = 2⁸. So P7′ is a retrospective score, not evidence of an architectural
  floor.
- **Code fix after the runs.** `forms.public_ambiguous` listed imprimitive (2, 2, c) (codex_2 finding 11). It is not
  used by any experiment. It now filters for primitivity, with a test.

**Reproduction.** `.venv/bin/python -m src.genusleak.harness E1` (likewise E2–E5); `.venv/bin/python -m pytest -q
tests/test_genusleak.py`. Wall times: E1 6 s, E2 18 s, E3 9 s, E4 2 s, E5 10 s.

## 8. Candidate cards (rule 9)
- **GL-1 Cayley spectral detection.**
  - Target: λ_{χ_p} = mean (ℓ/p).
  - Missing lemma: detecting it below |H| samples.
  - Failure: Lemma GIRTH (returns are relations); Prop WALK (≳ δ²|H| walks for fixed-time return bits, counterfactual).
  - Scope: refuted for this model only. Endpoint-coefficient statistics are not covered.
- **GL-2 Infrastructure jump to R/2.**
  - Works given any multiple of R (Prop HR).
  - Failure: obtaining that multiple. The known ways are relations, an analytic window (|D|^{1/5}), or BSGS; GI is a
    heuristic abstraction of the generic case. The midpoint is useful only with observed constant frequency.
- **GL-3 Square roots / genus characters.** Square roots give factors under GRH sampling (SQ). Principal-genus
  membership is QR_N-equivalent (GEN). No route.
- **GL-4 CM volcano.**
  - Failure: curve construction. The known methods cost √(kN) trials or |D|^{1+o(1)} CM; no lower bound is claimed.
  - Also: the ambiguous ideal representative has the unknown degree p. Other representatives are again a
    class-group problem.
- **GL-5 Multiplier pencil / compositum.** Only more hidden-order groups (GB union bound). Cross-group collisions
  live in Cl(L), which is larger.

## 9. Most decisive open question
Is there a polynomial-time computable statistic of the reduced forms in the principal cycle of √(kN) (or of the
reduced forms reached by composition), *not* of the form "multiplicative coincidence among norms", whose correlation
with a hidden 2-torsion class (rather than with the QR_N-type characters (ℓ/p)) is non-negligible?

By Prop GEN, correlation with (ℓ/p) alone buys at most QR_N information. By Theorem GB, any generic construction of
the class itself is a relation.

## 10. Reviews

**GPT-6.1-sol blind round (`codex/genusleak_1.md` → `.out`).**
- Independent conclusions: the same negative verdict, with no verified qualifying method.
- It proved independently:
  - the generic extraction (success ⇒ relation) and the Ω(P^{1/3}) uniform-prime bound;
  - the rank-1 equivalence with odd order multiples;
  - the QR_N ↔ hidden-genus reduction;
  - the √(h/|Cl[2]|) collision method;
  - girth after cancelling formal relations;
  - the R/2 midpoint with a 2-adic issue.
- Merged from it:
  - the N = 731 midpoint counterexample;
  - the Cl(−56) ≅ C₄ example (genus data does not locate involutions);
  - the unconditional random-summation estimator (credited to Srinivasan 1998);
  - HMT 1998 for SQRT ≤ FACT;
  - Sutherland's thesis and Hhan 2024;
  - the sampling-oracle treatment;
  - Sutherland's one-root CM method.
- Its open question: locate a proper square form by coefficient arithmetic faster than density D^{−1/4}. It matches
  ours.

**Opus proofs critic (read-only).**
- WRONG, fixed: the v1 imaginary finish in §3.3 ("square down one class"). Counterexample N = 2311·2671, 1/128
  success; we re-ran the critic's script and confirmed it.
- GAPs, fixed:
  - N(ε) = +1 in HR;
  - GI's spurious L factor (now √X; GI was later downgraded to heuristic by codex_2);
  - Bach generation needs the ramified primes in S;
  - Lemma 0 in the non-maximal order (A_pub redefined; proof now written);
  - the Hellinger second term, the 2^{μ−2} scaling, and the WALK scope;
  - the full D1 disclosure;
  - the GB exponent constant;
  - §4 nits (≤ 2 traces, crater hypothesis, GRH polylog);
  - the E2 n = 28 median is over 29 instances.
- Confirmed OK: Prop GEN (including even k and all 2-adic characters), GIRTH, GB, §5. All experiment numbers
  reproduce.

**GPT-6.1-sol adversarial referee (`codex/genusleak_2.md` → `.out`): MAJOR REVISION; all findings applied.**
- The headline is restated as a search outcome.
- SQUFOF even-r squares are imprimitive (36/90 E2 instances), and §3.4 is corrected.
- GB-X's converse needs reachability (C₄ ⊕ C₂ counterexample).
- WALK needs g₀ > 2 (N = 15 counterexample). Endpoint coefficients can leak at toy size (N = 21).
- GB's application hypotheses are listed. Sampling must be charged; the a·g^u simulation needs G = A⟨g⟩.
- GI is downgraded to a heuristic abstraction.
- HR is a midpoint *search*.
- 1/5 is an architecture cost, not a floor.
- GEN: the 2-adic table and the factored-small-k hypothesis are added; membership is distinguished from labelled
  characters.
- Lemma 0: the code bug (imprimitive public form) is fixed and the non-maximal proof completed.
- SQ: the GRH qualification is in the headline.
- CM: known costs, not lower bounds; other ideal representatives.
- Ring/pencil: "no gain" is restricted to the examined operations.
- The §6.2 survey is expanded (SQUFOF2, Lenstra–Pomerance, SPAR, rho, Hart, McKee, Hittmeir, Harvey–Hittmeir,
  Coppersmith-type, Murru–Salvatori).
- The deviation log is added. It independently regenerated all 380 experiment rows.

**Opus literature critic (read-only).**
- Corrected:
  - Bach 1995 is CMS Conf. Proc. 15, not Math. Comp.;
  - Bröker–Stevenhagen is Contemp. Math. 463 (2008);
  - the low-order assumption is BBF ePrint 2018/712;
  - BKSW concerns odd low order;
  - the 1/5 algorithm is credited to Shanks 1971;
  - the random-summation estimator is credited to Srinivasan;
  - Cox theorem numbers depend on the edition (marked UNVERIFIED).
- Added: Lenstra–Pomerance; the CJLN Coppersmith-type attack; McKee's representation collisions; class-field values;
  Pollard/Hittmeir public order approximations.
- VERIFIED: Lenstra 1982; Hallgren 2007; Srinivasan 1998; Gower–Wagstaff 2008; Schnorr–Lenstra; Seysen;
  Hafner–McCurley; Buchmann; Vollmer; Biasse–Jacobson; Bach 1990 (maximal orders); BBHM 2002; Adleman–McCurley 1994;
  Enge 2009; Sutherland 2011/2012; Lenstra 1996; Sutherland thesis Thm 2.3; Hhan Thm 4.1; Damgård–Koprowski;
  Rotem–Segev–Shahaf; Lemmermeyer 1994.
- UNVERIFIABLE: Lagarias Thm 6.10; HMT §3.1/§4 content (codex_1/codex_2 verified it in substance); the "O11" label;
  Schoof 1982 internals.
- Prior art:
  - Lemma 0 is classical (Shanks 1971; Schnorr–Lenstra; Cohen GTM 138 §8.6).
  - Prop SQ is folklore.
  - Prop GEN is folklore made explicit (CLT 2022 snippet).
  - GIRTH is the Seysen/Hafner–McCurley relation source.
  - HR is in substance Lenstra/Schoof/Hallgren; Murru–Salvatori discuss halving.
  - GB is Sutherland/Hhan plus GB-X.
  - No prior generic infrastructure bound was found, but GI is only heuristic here.
