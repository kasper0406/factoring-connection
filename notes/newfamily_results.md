# newfamily: searching for a fourth family of factor leakage (phase 24)

Agent `newfamily`, 2026-10-23.
- Pre-registration: `notes/newfamily_prereg.md` (commit `67fb78a`, before any battery instance).
- Code: `src/newfamily/{recip,modsym,battery,harness,diag}.py`. Tests: `tests/test_newfamily.py` (11 pass).
- Raw data: `experiments/raw/newfamily_main_n{32,48,64,128}.jsonl.gz`. Summaries:
  `experiments/analysis/newfamily_main_n{32,48,64,128}.json`, `experiments/analysis/newfamily_salie_sweep.json`.
- GPT-6.1-sol: blind round `codex/newfamily_1.md`, report `codex/newfamily_1_report.md` (25 candidates; copied from its
  scratch dir); adversarial round `codex/newfamily_2.{md,out}` (§7).
- Two read-only Opus critics, proofs and literature (§8).
- Labels: **KNOWN** (with reference), **PROVED** (proof here; routine unless stated; "no prior formulation found by
  us" is the strongest novelty claim), **HEURISTIC**, **CONDITIONAL**, **experimentally supported**, **OPEN**,
  **UNVERIFIED**. Negative literature statements ("no algorithm known") come from our and codex's searches and are
  UNVERIFIED as negatives.

Notation: N = pq balanced, ε = (−1/p) on N ≡ 1 (mod 4), QR₂ = [2 square mod N] on N ≡ 1 (mod 8).
- *Aggregate*: a function of N mod m for a fixed m (the fiber floor of Theorem PT-m, `notes/factorbits_results.md`
  §3.2). For symbols this means a product χ(p)χ(q) of one fixed character χ mod m.
- *FACT*: randomized polynomial-time equivalent to factoring N.
- *Class E*: Euclid-path functionals of a/N (§2.5).

---

## 0. Verdict (plain language)

**We found no fourth family.** Every candidate fell into one of four bins:
- an aggregate (it carries nothing beyond N mod m);
- factoring-equivalent;
- not computable in polynomial time even with the factors in hand;
- a continued-fraction statistic of a rational a/N ("class E").

Class E is the one bin that is computable, is not an aggregate, and is not known to be hard. A pre-registered battery on
2·10⁵ semiprimes per size found no signal there.

What is new here (no prior formulation found by us; the mathematics is routine):

1. **Theorem RA (PROVED).** Take any number field K and any abelian extension L/K, and evaluate the Artin symbol of
   the *rational* ideal N·O_K. The result depends only on N mod m, where m is fixed by the conductor.
   - This covers Jacobi symbols, all power-residue symbols computed by de Boer–Pagano at rational N, and Kronecker
     symbols in number fields, including extensions that are non-abelian over Q.
   - Corollary: pure trees (Theorem PT) that may also query such symbols at N-independent arguments are PT-m trees.
     Without splitting N they learn nothing beyond the fiber floor.
   - So the Jager–Schwenk escape (Jacobi is easy by reciprocity although generically hard) never opens a new channel.
     Reciprocity computes exactly the aggregates.
2. **The conjecture "every reciprocity-computable invariant is an aggregate" is false as stated. The counterexample
   shows where the factors hide.**
   - The Dedekind sum s(a, N) is computable by reciprocity in O(log N) steps. It splits exactly as
     s(a, pq) = s(a, p) + s(a, q) + U(a, N), where U is a sum over units only.
   - **Prop DSF (PROVED):** the local part s(a, p) + s(a, q) for any single a factors N deterministically.
     Equivalently, 6N·U(a, N) ≡ (a²+1)·φ(N)/(2a) (mod N).
   - What remains beyond the aggregate a + a⁻¹ mod N, (a/N) and N mod 24 is the alternating sum of the partial
     quotients of a/N, plus the parity of their number (Barkan–Hickerson–Knuth).
   - So the one non-aggregate reciprocity object we know is a continued-fraction statistic, and separating its local
     part is factoring.
3. **Root-type exponential sums are factoring-equivalent (Prop SAL, PROVED; decoder success rate measured).**
   - Salié sums S(a, b; N), and Kloosterman sums at square moduli K(a, b; N²), have closed forms over the square
     roots of ab.
   - An oracle for their normalised value factors N. A dyadic decoder needs only *constant* absolute precision
     (10⁻² was enough in ≥ 19/20 trials up to n = 128) and log₂N + 3 queries.
   - Quadratic Gauss sums are aggregates. Kloosterman sums at squarefree N have no known polylog algorithm even at a
     prime modulus.
4. **Battery (pre-registered, 84 menu statistics, n = 32, 48, 64, 128, 2·10⁵ instances each): no hit.**
   - The statistics are Jacobi symbols, Dedekind sums, continued-fraction and Euclid statistics, and level-11 modular
     symbols at six N-dependent arguments a(N).
   - All validity checks pass. The largest |z| is 3.25 at n = 48–128 and 3.76 at n = 32.
   - Effects with |ρ| ≥ 0.0153 (ε) and ≥ 0.0217 (QR₂) are excluded for this menu at n ∈ {48, 64, 128}.

**Most promising open candidate (§4):** class-E functionals at structured N-dependent arguments, in particular the
Rademacher integer Φ(a(N)/N) and fixed-level modular symbols. Prop DSF shows that the local data sit inside s(a, N)
in exactly computable form, and that separating them is factoring. Neither a mechanism nor a barrier for this class is
known.

---

## 1. Candidate table

The "Ref" column names the source (KNOWN) or the place where we prove the claim. cx = GPT-6.1-sol blind report
(`codex/newfamily_1_report.md`, entry number); nonsmooth/analytic/ecount/factorbits = earlier project notes.

**Verdict key**
- **AGG**: aggregate.
- **FACT**: factoring-equivalent.
- **NC**: no polynomial algorithm known.
- **E**: class E (computable, not an aggregate, OPEN).

### Reciprocity / Euclid-computable

| # | Invariant | Poly-time from N? | Verdict | Ref |
|---|---|---|---|---|
| 1 | Jacobi/Kronecker (c/N), fixed c; Zolotarev sign | yes (reciprocity) | AGG | KNOWN; Thm RA (K = Q) |
| 2 | Artin symbol of N·O_K in any abelian L/K | yes for fixed (K, L) | AGG | **Thm RA** |
| 3 | Power-residue symbol (α/(N))_n, fixed α ∈ Z[ζ_n] (de Boer–Pagano ISSAC 2017) | yes (DP heuristic-expected-poly; table lookup for fixed α) | AGG; depends on N mod Norm(α)·… | Thm RA; tested: cubic α = 1+3ω, 2+3ω, 4+3ω (Norm 7, 7, 13) |
| 4 | Rational quartic residuosity ρ₄(t²; N) := ρ₄(t²; p)ρ₄(t²; q) | not known | NC; an oracle decides ε (with (t/N) = −1) | cx 04 (PROVED there). Shows that "product of local data" ≠ aggregate |
| 5 | Dedekind sum s(a, N), fixed a | yes | AGG + size: 12s(a,N) − a/N − N/a − 1/(aN) + 3 = −12s(N mod a, a) | Prop RA2; cx 05 |
| 6 | s(a, N), N-dependent a; Rademacher integer Φ = (a + ā)/N − 12s(a, N) | yes | **E** (not AGG: Prop DS) | Prop DS, DSF; battery |
| 7 | Local part s(a,p) + s(a,q); unit part U(a, N) | — | **FACT** (one value, deterministic) | **Prop DSF** |
| 8 | Shifted Dedekind–Rademacher sums (Hall–Wilson–Zagier) | yes | E (heuristic: shifts evolve by the CF matrices mod a fixed denominator) | cx 07; not run |
| 9 | Higher Dedekind/Zagier cotangent sums; Bernoulli D_{1,3}; Fourier–Dedekind sums | yes in fixed dimension (Barvinok; Beck–Robins; Xin–Xu 2023) | OPEN (lattice-cone class, adjacent to (A)) | cx 08–10; not run |
| 10 | Dedekind eta multiplier; theta multiplier | yes | AGG (Jacobi × explicit 24th/8th root of unity) | Petersson/Knopp; cx 11–12 |
| 11 | 12N·s(a, N) mod 24N, gcd(N, 6) = 1 | yes | AGG in (a + ā mod N, (a/N), N mod 24) | Rademacher congruences; tested |
| 12 | Modular symbols {∞, a/N} of a fixed level M (e.g. 11a) | yes (Manin trick, O(log N)) | fixed a: AGG; N-dependent a: **E** | Stein, *Modular Forms* ch. 3; battery |
| 13 | Torus-knot signatures σ(T(a, N)), lens-space Casson–Walker / d-invariants | yes (Euclid recursions) | E (Dedekind-sum expressions) | Hirzebruch–Zagier; not run separately |
| 14 | Euclid-remainder residue statistics of (N, a(N)) | yes | E | battery (`rem1mod4`, `rem0mod3`) |
| 15 | Jacobi symbols of N-dependent elements (a(N)/N) | yes | local product χ_{a(N)}(p)χ_{a(N)}(q) whose conductor moves with N; not covered by RA; battery null | battery (`*_jac`) |
| 16 | Jacobi vector ((x/N))_{x ≤ H} | yes | AGG (modulus lcm(4x)); o(1) advantage for H ≤ c·log N under GRH (Thm L + Parseval) | RA-PT + selfred Thm L |
| 17 | Root numbers of N-twisted elliptic families (e.g. y² = x(x−1)(x−N): w_p = −(−1/p)) | when computable | AGG (local root numbers at tame primes are fixed-conductor Kronecker characters; Rohrlich 1993; Hilbert reciprocity) | HEURISTIC in general; example PROVED |
| 18 | Rédei symbols [a, b, N] | needs an ideal of norm N, i.e. √a mod N | NC (Artin symbol of N·O_K itself is RA-trivial) | Stevenhagen 2022 |

### Exponential sums

| # | Invariant | Poly-time from N? | Verdict | Ref |
|---|---|---|---|---|
| 19 | Quadratic Gauss G(a; N), G(a, b; N); Landsberg–Schaar | yes | AGG: ε_N (a/N)√N · e(−(4a)⁻¹b²/N) | KNOWN; cx 13 |
| 20 | Fixed-rank Gaussian / abelian Chern–Simons Z_A(N) | yes | AGG: ε_N^d (det A/N) | cx 14 |
| 21 | Salié S(a, b; N) | iff factoring | **FACT**, constant precision | **Prop SAL**; cx 16 (precision N⁻⁴, one query) |
| 22 | Kloosterman K(a, b; N²) | iff factoring | **FACT**, constant precision | **Prop SAL** (square case) |
| 23 | Kloosterman K(a, b; N), N squarefree | no known algorithm, even for prime N | NC; K(a,b;pq) = K(a q̄, b q̄; p)K(a p̄, b p̄; q) | Katz 1988 (rank-2 sheaf, no bounded-genus curve); cx 15; F_{2^m}, F_{3^m} are the elliptic exceptions (Lachaud–Wolfmann 1990) |
| 24 | Cubic Weyl sums Σ e(ax³/N); cubic Gauss/Jacobi sums mod N | no | FACT (analytic C3); zero-test decides (−3/p) | analytic §0 item 3; cx 17 |
| 25 | Non-quadratic characters of (Z/N)^* (and their Gauss/Jacobi sums) | no (needs a prime ideal above p) | NC; rational arguments make (·/(N))_n trivial | Thm RA remark |

### Half-integral weight / theta

| # | Invariant | Poly-time from N? | Verdict | Ref |
|---|---|---|---|---|
| 26 | r₂(N) = 16[ε = 1] | no | NC (open, = ε) | factorbits §2.3; cx 19 |
| 27 | r₃(N) = 12h(−4N) etc.; Hurwitz H(N); ternary representation numbers | no (class numbers, L[1/2]) | FACT-hard direction (h(−4N) ⇒ ambiguous forms); family Q/M | Gauss; Shanks; cx 20 |
| 28 | r₄(N) = 8(N+1+p+q) | iff factoring | FACT | Jacobi; cx 21; nonsmooth |
| 29 | Half-integral weight cusp coefficients c(N) (Waldspurger c(N)² ∝ L(½, f⊗χ_N)) | no (≈ √N terms) | NC; not local at p | Waldspurger 1981 |
| 30 | Cohen–Eisenstein H(r, N), r ≥ 2 (L(1−r, χ_D)) | no known | NC; no Euler splitting over p, q | Cohen 1975 |

### Quantum / topological

| # | Invariant | Poly-time from N? | Verdict | Ref |
|---|---|---|---|---|
| 31 | WRT of lens spaces L(N, a); Chern–Simons of lens spaces | yes | AGG + E (Gauss sums + Dedekind sums) | Jeffrey 1992 |
| 32 | Weil representation traces of SL₂(Z/N) | yes | AGG (Gauss sums); \|tr\|² = #fixed vectors is visible only through gcd(det(γ−1), N) | — |
| 33 | Trefoil colored Jones / Kontsevich F(e(a/N)) (quantum modular) | approximations via reciprocity | E (cocycle evaluated along the CF of a/N); exact value has φ(N) coordinates | Zagier 2010 |
| 34 | Kashaev invariant of 4₁ / hyperbolic knots at e(1/N) | no (O(N) recurrences; asymptotics only) | NC | Garoufalidis–Zagier; cx 22 |
| 35 | Habiro-ring values at ζ_N: f(ζ_N) ≡ f(ζ_q^b) mod (1 − ζ_p) | exact values exponential | NC (local structure exists algebraically, invisible archimedeanly) | Habiro 2004 |

### Dynamical / analytic

| # | Invariant | Poly-time from N? | Verdict | Ref |
|---|---|---|---|---|
| 36 | Period of x ↦ x²+c mod N | no | NC; lcm of local periods (rho, family M) | cx 23 |
| 37 | ord_N(a), selectable a; λ(N) | no | FACT (Miller 1976; Shor post-processing). Fixed base 2 factors a constant fraction of N heuristically | cx 24 |
| 38 | CF of N^{1/3} | yes (prefix) | relations x³ − Ny³ = m (NFS-type, family M) | — |
| 39 | CF of e^{1/N} = [1; N−1, 1, 1, 3N−1, …] | yes | trivial function of N | Euler |
| 40 | Farey / modular-group word of a/N | yes | E | = #6 |
| 41 | Fixed transfer matrix tr(A^N) over F_ℓ; Lehmer–Pierce Res(x^N − 1, P) mod ℓ | yes | AGG (period of the sequence) | cx 25 |
| 42 | Res(x^N − 1, P) mod N ≡ Res(x^q − 1, P) (mod p) | yes | N-dependent-exponent algebra (leakmap's class) | — |

### Wild / cross-domain

| # | Invariant | Poly-time from N? | Verdict | Ref |
|---|---|---|---|---|
| 43 | Duadic-code existence at length N (all p ≡ ±1 mod 8) | no | NC; it is the target itself | — |
| 44 | Linear complexity / autocorrelation of the Jacobi sequence (x/N) | O(N) | NC; divisor detection (U3) | Ding et al.; nonsmooth |
| 45 | Necklace / Lyndon counts of length N; cyclic codes of length N | no | NC (divisor sums, ord_d(2)) | — |
| 46 | Spanning trees of circulant graphs C_N(1, a) | exact: exponential size; mod ℓ: AGG | AGG / NC | Lehmer–Pierce divisibility |
| 47 | QR-graph subgraph counts on Z/N (tensor product of local graphs) | no | NC (class (a) of nonsmooth) | — |
| 48 | Bost–Connes polylogs Li_k(e(a/N)) | yes | trivial (Bernoulli polynomial of a/N) | — |
| 49 | Barvinok counts of N-dilated fixed-normal polytopes | yes | AGG (Ehrhart quasi-polynomial in N) | Ehrhart; Barvinok 1994 |
| 50 | Same, with variable normals | yes | E / lattice; divisor detection needs a hull with Θ(N^{1/3}) vertices | analytic §0 |
| 51 | Rademacher symbol / ζ_C(0) on forms of discriminant t² − 4 with N \| t² − 4 (Richaud–Degert) | yes | family Q (genus characters via relations); handed to `genusleak` | Zagier 1975 |
| 52 | Partition function p(N) | N^{1/2} | NC; no multiplicative structure | Rademacher |

**Reading the table.**
- Everything computable splits into AGG (rows 1–3, 5, 10, 11, 16, 17, 19, 20, 31, 32, 41, 49) or class E (rows 6,
  8, 12–14, 33, 40, 50).
- Rows 9 and 15 are poly-time and OPEN but fit neither pattern cleanly: 9 is a lattice-cone class, 15 is a symbol at a
  moving argument.
- Everything whose closed form contains a square root of a nonlinear congruence (rows 7, 21, 22, 24) is FACT.
- Everything else (rows 4, 18, 23, 25–30, 34–37, 43–47, 52) has no polynomial algorithm known even with p, q.
- Row 42 belongs to leakmap's N-dependent-exponent class, and row 51 to genusleak's family Q.

---

## 2. Theorems and propositions

### 2.1 Theorem RA (reciprocity aggregate). PROVED (standard class field theory; no prior statement in this form found by us)

**Statement.** Let K be a number field and L/K a finite abelian extension of conductor 𝔣 = 𝔣₀𝔣_∞. Let m > 0 generate
𝔣₀ ∩ Z. Let ψ: I_K^{𝔣} → Gal(L/K) be the Artin map. For positive integers N coprime to m, the value ψ(N·O_K) depends
only on N mod m. Hence N ↦ ψ(N·O_K) is a homomorphism (Z/m)^* → Gal(L/K), and ψ(pq·O_K) = χ(p)χ(q) for this fixed
χ.

*Proof.*
1. N·O_K is coprime to 𝔣₀: a prime 𝔭 | 𝔣₀ contains m, so 𝔭 ∩ Z = ℓZ with ℓ | m.
2. Let N ≡ N′ (mod m). Pick a positive integer M with NM ≡ 1 (mod m).
3. Then NM − 1 ∈ mZ ⊂ 𝔣₀ and NM is a unit at every 𝔭 | 𝔣₀. So NM ≡* 1 mod 𝔣₀, and NM > 0 at every real place.
4. Hence (NM) ∈ P_{K,𝔣}, and ψ(N·O_K)ψ(M·O_K) = 1. The same holds for N′.
5. So ψ(N·O_K) = ψ(M·O_K)^{−1} = ψ(N′·O_K). Multiplicativity of ψ gives the homomorphism. ∎

**Instances.**
- K = Q: Jacobi/Kronecker symbols. This is Kronecker–Weber, as in nonsmooth §5.4.
- K ∋ ζ_n, L = K(α^{1/n}): the power-residue symbol (α/(N))_n, which is what de Boer–Pagano compute for a rational
  modulus.
- Quadratic, quartic and cubic symbols in Z[i] or Z[ω] with a rational modulus.
- L/Q need not be abelian. Only L/K must be.
- **Test (exact):** for α ∈ {1+3ω, 2+3ω, 4+3ω}, (α/(N))₃ was computed from the factorisation for 1 500 pairs of primes
  < 1 500 per α. It depends only on N mod Norm(α), and it takes all three values. For α = 1+3ω, cubic reciprocity
  gives the explicit form (α/(N))₃ = χ₇(N), with χ₇ the cubic character mod 7.

**Corollary RA-PT (pure trees with symbol nodes). PROVED.**
- Model: Theorem PT's pure trees, extended by nodes that branch on a symbol ψ_v(c_v, N). Here c_v is an algebraic
  integer computed by an N-independent straight-line program over O_{K_v}, and ψ_v is any of the symbols above at the
  rational ideal (N).
- Let m = lcm of the conductors of all symbol nodes; the tree is finite.
- Given N mod m, every symbol answer is fixed. The tree is therefore a PT-m tree: the root reads N mod m, and the
  symbol nodes become hard-wired branches.
- So Theorem PT-m applies. Without splitting, the tree learns at most the fiber floor at m. For m ≤ X^{1−η} this is
  o(1) advantage on ε under GRH (selfred Theorem L).
- For conductors ≫ X the floor statement is vacuous. The answer is still a product χ(p)χ(q), and its correlation with
  ε is |E[εχ(p)]|·|E[χ(q)]|-type. Bounding that needs prime character sums of large conductor (Heath-Brown's
  quadratic large sieve is the natural tool); not done.

**Scope (what RA does not cover).**
- (i) Symbols at ideals 𝔫 | N not extended from Z, e.g. (N, ζ − r). Building them needs a root of Φ_n mod N, which is
  Frobenius/factoring-type (nonsmooth Prop FROB).
- (ii) Symbols at N-dependent arguments α(N) (row 15). The value is still χ_{α(N)}(p)χ_{α(N)}(q), but χ moves with N,
  so there is no fixed floor. This is the general problem in disguise: a single (a(N)/N) could be any bit.
- (iii) "Products of local data" in general. cx 04's quartic test ρ₄(t²;p) is local at p but is not a fixed-conductor
  character in p, and its product decides ε.

### 2.2 Prop RA2 (fixed-argument Euclid recursions). PROVED (trivial)

**Statement.**
- Hypotheses: f(a, c) = f(a mod c, c), and f(a, c) = R(a, c) ⊙ f(c, a). Here R depends only on (a mod m, c mod m)
  and on a, c as real numbers.
- Conclusion: for fixed a, f(a, N) = R(a, N) ⊙ f(N mod a, a) is a function of N mod lcm(a, m) and the real number N.

**Example.** The Dedekind case is 12s(a,N) − a/N − N/a − 1/(aN) + 3 = −12s(N mod a, a) (cx 05).

**Meaning.** Any reciprocity-type recursion at a *fixed* argument is floor plus archimedean size. Only N-dependent
arguments can go beyond it.

### 2.3 Prop DS (Dedekind decomposition) and Prop DSF (the local part factors). PROVED; tested

**DS1. Decomposition.** For N = pq with p ≠ q, s(a, N) = s(a, p) + s(a, q) + U(a, N), where
U(a, N) = Σ_{k ∈ (Z/N)^*} ((k/N))((ak/N)).
- *Proof.* Inclusion–exclusion over k ≡ 0 mod p and k ≡ 0 mod q. For k = pj, ((pj/N)) = ((j/q)) and
  ((apj/N)) = ((aj/q)), so that block sums to s(a, q). Symmetrically for q. ∎
- Tested exactly for 4 small N and all a.

**DS2. Only a continued-fraction statistic is non-aggregate.**
- Write a/N = [0; a₁, …, a_r], 0 < a < N, and ā = a⁻¹ mod N in (0, N). Then
  12 s(a, N) = Σ_{i=1}^{r} (−1)^{i+1} a_i + (a + ā)/N − 2 + (−1)^r (Barkan–Hickerson–Knuth; verified on 3 000 random
  pairs).
- With the classical congruences 12N s ≡ a + ā (mod N) and 12N s ≡ N + 1 − 2(a/N) (mod 8) (N odd), and
  12N s ≡ 0 (mod 3) when 3 ∤ N (all tested), the value 12N s mod 24N is an aggregate of a.
- Equivalently, the Rademacher integer Φ = (a + ā)/N − 12s(a, N) = −Σ(−1)^{i+1}a_i + 2 − (−1)^r carries everything
  else.

**DSF. Separating the local part is factoring.**
- Statement: for any a with gcd(2a(a²+1), N) = 1, the exact value L = s(a, p) + s(a, q) determines p + q, and hence
  the factors (one value; deterministic). Conversely the factors give L by reciprocity at p and q.
- *Proof.*
  1. k_p := 6p·s(a, p) ∈ Z, and 2k_p ≡ a + a⁻¹ (mod p) (the Rademacher congruence at the prime p > 3). Same for q.
  2. X := 6N·L = k_p q + k_q p ∈ Z.
  3. Then 2aX ≡ q(a² + 1) (mod p) and 2aX ≡ p(a² + 1) (mod q).
  4. Hence 2aX ≡ (p + q)(a² + 1) (mod N), and p + q = 2aX(a² + 1)⁻¹ mod N, because p + q < N. ∎
- Equivalently, 6N·U(a, N) ≡ (a² + 1)φ(N)/(2a) (mod N) (verified on 500 random cases). The unit part of a Dedekind sum
  hides φ(N).
- Test: `test_prop_DSF_local_part_of_dedekind_sum_factors`, 30 instances, n up to 256, all factored.

**Meaning.**
- The reciprocity-computable s(a, N) is "aggregate + entangled". It contains the local values s(a, p), s(a, q)
  exactly, mixed with U.
- Pulling them apart is factoring, just as recovering (a/p) from (a/N) is QR-hard. Here the hardness is FACT, not
  merely QR.
- The conjecture "every reciprocity-computable invariant is an aggregate" therefore holds for the multiplicative
  (symbol) class (Theorem RA; cx's Theorem R is the same statement for a finite-group axiomatisation).
- It fails for the additive (quasimorphism) class. That class is exactly the continued-fraction statistics of class E.

### 2.4 Prop SAL (Salié at N and Kloosterman at N²). PROVED; decoder success measured

**Statement.**
- Let N = pq be odd, k ∈ {1, 2}, M = N^k, and c a unit mod M. Let e be the square root of 1 mod M with e ≡ 1 (mod p^k)
  and e ≡ −1 (mod q^k). Define T_M(c) = Σ_{v² ≡ c² (M)} e(2v/M).
- Then:
  - (i) S(c², 1; N) = ε_N √N·T_N(c). This is the Salié closed form, valid for composite odd N (cx 16 gives a Fourier
    proof; Iwaniec, *Topics in Classical Automorphic Forms* Lemma 4.10 for the prime case). Tested numerically.
  - (ii) K(c², 1; N²) = N·T_{N²}(c). This is the evaluation of Kloosterman sums at even prime powers (Iwaniec–Kowalski
    §12.3) combined with twisted multiplicativity. Tested numerically for N = 15, 21, 35.
  - (iii) T_M(c) = 2cos(4πc/M) + 2cos(4πce/M).
- **Decoder.** Query T_M(2^k c) with absolute error δ for k = 0, …, K = bitlen(M) + 2.
  1. Subtracting the known first cosine leaves d_k = ‖x_k‖ with x_k = 2^{k+1}ce/M mod 1.
  2. Reconstruct x₀ backwards: x_k ∈ {x_{k+1}/2, (x_{k+1}+1)/2}, chosen by ‖x_k‖ = d_k. The error halves at each
     step.
  3. Round to a multiple of 1/M to get ±2ce mod M.
  4. Then gcd(ce − c, N) = p, since ce ≡ c (mod p) and ce ≡ −c (mod q).
- **Correctness when it succeeds:** PROVED (the output is gcd-verified).
- **Success probability:** HEURISTIC. A step fails only when x_{k+1} lies within O(δ) of ½. Measured (all with
  20/20 trials unless noted):

| n | δ = 10⁻² | 3·10⁻³ | 10⁻³ | 10⁻⁴ | 10⁻⁶ |
|---|---|---|---|---|---|
| 64 | 20/20 | 20 | 20 | 20 | 20 |
| 128 | 19/20 | 20 | 20 | 20 | 20 |
| 256 | 15/20 | 18 | 18 | 19 | 20 |

  (Queries: 67, 131, 259; `newfamily_salie_sweep.json`. Kloosterman-at-N² version tested at n = 40.)
- cx 16 gives the complementary single-query decoder at precision N⁻⁴.
- **Converse:** the factors give S and K(·; N²) to any precision in polynomial time (square roots mod p, q; Hensel;
  CRT).

**Dichotomy (HEURISTIC, as an organising principle).** For complete exponential sums mod N:
- a quadratic phase gives an aggregate (Gauss);
- a phase whose stationary points are roots of a nonlinear congruence gives FACT (Salié: v² ≡ ab; K at N²; cubic Weyl
  through Jacobi sums in Z[ω], analytic C3);
- the remaining sums (Kloosterman at squarefree N, higher degree) have no known algorithm even at a prime modulus.

This matches the degree-2/3 line of Cai–Chen–Lipton–Lu noted in analytic §0.

### 2.5 Class E and what would be needed for a floor theorem. HEURISTIC / CONDITIONAL

**Definition of class E.** Class E consists of the additive functionals Σ_k F(a_k, γ_k mod m₀) along the CF path of
a/N, where γ_k are the convergent matrices and m₀ is fixed. It contains:
- Φ and the Dedekind sums (DS2);
- level-M modular symbols (Manin symbols are (q_{k−1} : q_k) mod M);
- shifted Dedekind–Rademacher sums (shifts transported mod their denominator; heuristic);
- Euclid-remainder residue counts;
- quantum-modular values computed through their cocycles.

**Random arguments.** For a public random a (independent of p, q given N), every class-E statistic g has advantage
equal to that of G(N) = E_a[g]. G(N) is an average over a ∈ (Z/N)^*.
- An effective equidistribution theorem for CF functionals over (Z/N)^*, uniform in the factorisation of N with error
  N^{−δ}, would give advantage O(N^{−δ}).
- Bettin–Drappeau-type limit laws are the candidate tools. Whether their error terms are uniform in the factorisation
  of an individual modulus is UNVERIFIED.

**Structured arguments.** For structured a(N) (⌊√N⌋, ⌊N^{2/3}⌋, Fermat residues, 2^{(N−1)/2}) we have no theorem.
The CF of a/N gives congruences a·q_k ≡ ±r_{k+1} (mod N) with q_k r_{k+1} ≍ N/a_{k+1}. These are never small, so the
relation mechanism that made the CF of √N leak at n = 32 (eps §5) is absent. This is the registered prediction P1,
and the battery confirms it.

---

## 3. Battery results (pre-registered; `notes/newfamily_prereg.md`)

**Setup.** 2·10⁵ instances of D_n per size (QR₂: the ≈ 10⁵ with N ≡ 1 mod 8). 98 statistics. Stratified permutation z,
strata N mod 2⁸; a hit needs |z| ≥ 6. Cost was 1.4–2.7 ms per instance; the run took about 5 min wall on 6 workers.

| n | target | instances | V1 planted z (exp.) | V2 PCs: fires / violations / z | V3 shuffle mean z² | menu mean z² | max \|z\| (menu+NC) | hits | ρ_MDE |
|---|---|---|---|---|---|---|---|---|---|
| 32 | ε | 200 000 | 21.2 (22.4) | pc2 98 099/0/440; pc3 99 233/0/445 | 0.89 | 2.18* | 5.69 (`nc_ded5_bits`)†; menu 3.40 | 0 | 0.0153 |
| 32 | QR₂ | 99 580 | 15.6 (15.8) | pcq2 49 083/0/299; pcq3 51 432/0/314 | 1.02 | 1.70* | 3.76 | 0 | 0.0217 |
| 48 | ε | 200 000 | 21.2 (22.4) | 26 122/0/174; 85 227/0/387 | 0.97 | 1.20 | 2.80 | 0 | 0.0153 |
| 48 | QR₂ | 99 976 | 15.7 (15.8) | 9 584/0/103; 45 708/0/290 | 0.86 | 1.13 | 2.57 | 0 | 0.0217 |
| 64 | ε | 200 000 | 21.2 (22.4) | 1 891/0/43.7; 747/0/27.4 | 1.22 | 1.51 | 3.25 | 0 | 0.0153 |
| 64 | QR₂ | 99 724 | 15.6 (15.8) | 665/0/25.9; 392/0/19.8 | 0.82 | 0.95 | 3.19 | 0 | 0.0217 |
| 128 | ε | 200 000 | 21.2 (22.4) | never fire (as predicted) | 1.16 | 0.48 | 1.53 | 0 | 0.0153 |
| 128 | QR₂ | 100 211 | 15.7 (15.8) | never fire | 0.91 | 0.83 | 2.48 | 0 | 0.0216 |

**Validity.**
- All size × target cells are valid by the registered rules. V1: every |z − expected| ≤ 1.2. V2: 0 violations in
  every cell; at n ≤ 64 the PCs fire with z ≥ 19.8. V3: shuffle mean z² ∈ [0.82, 1.22].
- `nc_N_mod8` and `nc_N_mod256` are exactly degenerate everywhere (P2).
- Also degenerate, as anticipated in the prereg: `ferm_jac`, `euler_jac`, and the size statistics of `sqrt`/`cbrt2`
  whose first partial quotient fixes them. `nc_ded5_sign`, `ferm_dsign`, `ferm_altsign` are degenerate too (constant
  sign). Degenerate statistics detect nothing.

**Conclusion (registered rule).**
- **Refuted for this menu** at n ∈ {48, 64, 128} for both targets. No statistic correlates with ε at |ρ| ≥ 0.0153, or
  with QR₂ at |ρ| ≥ 0.0217, beyond the N mod 2⁸ floor (80% power at |z| = 6).
- n = 32 also has no hit, but see the two footnotes.
- P1, P2 and P3 were all met.

\* **n = 32 inflation.** Only 1 491 primes qualify per side, so each prime appears in ≈ 270 instances, and the
permutation z are anti-conservative (eps §5.3). We did not rerun the cluster-robust variance. No |z| reached 6, so the
conclusion does not depend on it.

† **`nc_ded5_bits` (z = −5.69) and `nc_lead16` (−4.97) at n = 32, ε.** These are size statistics, NC by design. They
reflect a finite-set artefact: among the 1 491 sixteen-bit primes with the top two bits set, the share ≡ 1 (mod 4) is
0.504/0.503/0.515/0.477 in the four quarters by size (corr(size, [p ≡ 1 mod 4]) = −0.0073). Combined with prime reuse,
this gives an expected |z| ≈ 4–5. It is the same kind of effect as the QR₂ finite-set offset in eps §3. Post-hoc
explanation; |z| < 6, so not a hit. Both are null at n ≥ 48 (|z| ≤ 2.2).

**Nonunits (rule 3).**
- All 12 instances with gcd(a(N), N) > 1 are at n = 32: 6 from `cbrt2` and 6 from `hash`. The ACCIDENTAL null
  expectation Σ(1/p + 1/q) is 7.0 per argument, so both match.
- Four arguments are structurally coprime:
  - `sqrt`: ⌊√N⌋ ∈ (p, q) for balanced N;
  - `ferm`: p | (r+1)² would need r + 1 ∈ {p, q};
  - `euler`: a unit;
  - `gold`: ⌊αN⌋ ≡ ⌊p{αq}⌋ (mod p), which vanishes only if {αq} < 1/p. Since α = (√5−1)/2 is badly approximable,
    {αq} > 1/(√5q), so the event is rare. Observed 0.
- No STRUCTURAL splits; n ≥ 48 has none.

---

## 4. The most promising open candidate, and why

**Candidate.** Class-E functionals at structured N-dependent arguments, sharpened by Prop DSF. Concretely: the
Rademacher integer Φ(a(N)/N), or the level-11 modular symbol, at a(N) ∈ {⌊√N⌋, ⌊N^{2/3}⌋, (⌊√N⌋+1)² − N, …}.

**Why this one.**
- (1) It is the only class of poly-time, reciprocity-computed invariants that is provably **not** an aggregate
  (Theorem RA covers all symbols; this class escapes it).
- (2) The factors are *inside* it in exact form. s(a, N) = s(a, p) + s(a, q) + U(a, N), and U hides φ(N) mod N
  (DSF). That is more than can be said for any other computable object in the table.
- (3) No barrier is known. PT does not apply (Euclid on integer representatives is not generic-ring). RA does not
  apply (additive, not multiplicative). The relation mechanism is absent (§2.5).

**Why it is probably not a fourth family.**
- The battery is null at |ρ| ≥ 0.015.
- For random arguments a floor would follow from an effective, factorisation-uniform equidistribution theorem for CF
  statistics in (Z/N)^* (§2.5).
- DSF shows the local part is exactly as hard as factoring. A leak would have to come from a *partial* separation
  that the CF path performs by itself, and we see no mechanism for that.

**Sharpest question.** Is there a polynomial-time functional of the CF path of a/N (a = a(N) public) whose correlation
with the "local" Dedekind value s(a, p) + s(a, q) exceeds that of the trivial predictor s(a, N)? Equivalently, does
the unit part U(a, N) have a computable component correlated with φ(N) mod small moduli?

**The next experiment.**
- Codex's lattice-cone sums (row 9: Zagier cotangent, Fourier–Dedekind, via Barvinok/Xin–Xu) are the untested
  members.
- They need an optimised cone implementation and were not run.

---

## 5. Candidate cards (rule 9, condensed)

| Mechanism | Model / scope | Obstruction (precise) | Status |
|---|---|---|---|
| Symbols at (N) | any abelian L/K, N-independent argument | value factors through (Z/m)^* (Thm RA) ⇒ PT-m floor | refuted for this model (PROVED) |
| Dedekind / Rademacher at fixed a | Euclid recursion | N mod a + size (RA2) | refuted (PROVED) |
| Dedekind local part | exact L(a,N) | ≡ FACT (DSF) | FACT (PROVED) |
| Salié / Kloosterman at N² | approximate-value oracle | ≡ FACT at constant precision (SAL) | FACT (PROVED; success measured) |
| Class E at N-dependent a | 84 statistics, 4 sizes | no signal at \|ρ\| ≥ 0.0153 | experimentally refuted for this menu |
| Lattice-cone sums | fixed-dimension Barvinok | none known | not yet tested |

---

## 6. Deviations from the pre-registration

- `src/newfamily/diag.py` (DIAGNOSTIC, not imported by the battery) was extended after the freeze:
  - the Kloosterman-at-N² oracle;
  - the Prop DSF decoder;
  - a fix of the cubic-symbol embedding (k₂ = 2j₂). This came before the RA test was used in any conclusion.
- `tests/test_newfamily.py` gained the corresponding tests.
- The frozen battery modules (`recip`, `modsym`, `battery`, `harness`) were not changed.
- The n = 32 cluster-robust variance named in prereg §5 was not computed. It was not needed, since no |z| ≥ 6.
- Run order was 48, 64, 128, 32. This has no effect on the results.

## 7. GPT-6.1-sol rounds

**Blind round (newfamily_1, xhigh).** Independent list of 25 candidates (`codex/newfamily_1_report.md`). It agrees
with every overlapping verdict here. Its contributions:
- quartic residuosity as a non-aggregate local product that decides ε (row 4);
- the single-query Salié decoder;
- the lattice-cone Dedekind generalisations (rows 8–9);
- Theorem R, a finite-group version of RA;
- a descriptive 5 × 10⁵-instance check of Φ mod 5/7, asinh Φ and two Manin coordinates (all |corr| ≤ 0.0115, not
  replicated).

Its top three open candidates are Φ, fixed-level modular symbols, and shifted Dedekind–Rademacher sums. These are
class E, which matches §4.

**Adversarial round (newfamily_2).** PENDING.

## 8. Opus critics

PENDING.

## 9. Reproduction (project root)

```
.venv/bin/python -m pytest -q tests/test_newfamily.py
PYTHONPATH=. .venv/bin/python -m src.newfamily.harness --n 48 64 128 32 --range main --count 200000 --workers 6
PYTHONPATH=. .venv/bin/python -m src.newfamily.sweep          # Prop SAL precision sweep (~3 min)
```
