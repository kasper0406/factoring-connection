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

**We found no fourth family.** This is the outcome of a scoped search, not a classification theorem. Most candidates
fell into one of four bins:
- an aggregate (it carries nothing beyond N mod m);
- factoring-equivalent or factoring-hard;
- no polynomial algorithm known from N alone;
- a continued-fraction statistic of a rational a/N ("class E").

Class E is the principal computable, non-aggregate, not-known-hard bin, and the one we tested. Other computable bins
also remain OPEN (§1, "Reading the table"):
- fixed-dimension lattice-cone sums;
- Jacobi symbols at moving arguments;
- N-dependent-exponent algebra;
- automata on the bits of N.

A pre-registered battery on 2·10⁵ semiprimes per size found no signal in class E.

Results (items 2–3: no prior formulation found by us, mathematics routine; item 1 is standard class field theory, and
the contribution there is the corollary for pure trees):

1. **Theorem RA (PROVED; standard, cf. Lemmermeyer, *Reciprocity Laws* 2000).** Take any number field K and any abelian extension L/K, and evaluate the Artin symbol of
   the *rational* ideal N·O_K. The result depends only on N mod m, where m is fixed by the conductor.
   - This covers Jacobi symbols, the fixed-numerator power-residue symbols of de Boer–Pagano at the full rational
     ideal (N), and Kronecker symbols in number fields. L/Q need not be abelian.
   - Corollary RA-PT: pure trees (Theorem PT) that may also query such symbols at N-independent arguments are PT-m
     trees. Without splitting they are bound by PT-m's fiber-floor inequality.
     - The fiber floor refers to the *hidden* point p mod m. It can contain ε when 4 | m.
     - The symbol answers themselves are functions of N mod m. By selfred Theorem L (GRH) they carry o(1) advantage on
       ε when m ≤ X^{1−η}.
   - RA does not cover symbols at N-dependent arguments (Jager–Schwenk's own setting, row 15), nor the additive class
     of item 2.
2. **The conjecture "every reciprocity-computable invariant is an aggregate" is false as stated. The counterexample
   shows where the factors hide.**
   - The Dedekind sum s(a, N) is computable by reciprocity in O(log N) steps. It splits exactly as
     s(a, pq) = s(a, p) + s(a, q) + U(a, N), where U is a sum over units only.
   - **Prop DSF (PROVED):** the local part s(a, p) + s(a, q), for any single a with gcd(2a(a²+1), N) = 1, factors N deterministically.
     Equivalently, 6N·U(a, N) ≡ (a²+1)·φ(N)/(2a) (mod N).
   - s(a, N) is determined by the easily computed data a + a⁻¹ mod N, (a/N) and N mod 24, together with the
     alternating sum of the partial quotients of a/N and the parity of their number (Barkan–Hickerson–Knuth).
   - So the principal non-aggregate reciprocity object we know is a continued-fraction statistic, and separating
     this local part is factoring.
3. **Root-type exponential sums are factoring-equivalent (Prop SAL, PROVED; decoder success rate measured).**
   - Salié sums S(a, b; N), and Kloosterman sums at square moduli K(a, b; N²), have closed forms over the square
     roots of ab.
   - An oracle for their normalised value factors N. PROVED at precision N⁻⁴ with one query (cx 16), or with an
     exact oracle. A dyadic decoder (HEURISTIC; success rates measured) needs precision δ = O(1/log N), i.e. about
     log log N bits, and log₂N + 3 queries. Its success rate falls at fixed δ: at δ = 10⁻² it is 20/20 at n = 64,
     15/20 at n = 256, and 9/20 at n = 1024 (critic run).
   - Quadratic Gauss sums are aggregates. Kloosterman sums at squarefree N have no known classical polylog algorithm even at a
     prime modulus.
4. **Battery (pre-registered, 84 menu statistics, n = 32, 48, 64, 128, 2·10⁵ instances each): no hit.**
   - The statistics are Jacobi symbols, Dedekind sums, continued-fraction and Euclid statistics, and level-11 modular
     symbols at six N-dependent arguments a(N).
   - All validity checks pass.
   - Largest menu |z|: 3.25 at n = 48–128, and 3.40 (ε) / 3.76 (QR₂) at n = 32.
   - Overall maximum: 5.69, for the NC size statistic `nc_ded5_bits` at n = 32. This is a finite-prime-set size
     artefact (§3).
   - At n ∈ {48, 64, 128} the battery had 80% power to detect a *linear* association of a coded statistic with ε at
     |ρ| ≥ 0.0153 (QR₂: 0.0217). None was found.

**Most promising open candidate (§4):** class-E functionals at structured N-dependent arguments, in particular the
Rademacher integer Φ(a(N)/N) and fixed-level modular symbols. Prop DSF shows that one natural local/global split of
s(a, N) is exactly as hard as factoring. Neither a mechanism nor a barrier for this class is
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
| 18 | Rédei symbols [a, b, N] | needs an ideal of norm N (√a mod N), or, via Rédei reciprocity [a,b,N] = [a,N,b], a solution of x² − ay² = Nz² (a conic with coefficient N; needs the factorisation in general: Simon; Cremona–Rusin) | NC (the Artin symbol of N·O_K itself is RA-trivial) | Stevenhagen, MPCPS 172 (2022) |

### Exponential sums

| # | Invariant | Poly-time from N? | Verdict | Ref |
|---|---|---|---|---|
| 19 | Quadratic Gauss G(a; N), G(a, b; N); Landsberg–Schaar | yes | AGG: ε_N (a/N)√N · e(−(4a)⁻¹b²/N) | KNOWN; cx 13 |
| 20 | Fixed-rank Gaussian / abelian Chern–Simons Z_A(N) | yes | AGG: ε_N^d (det A/N) | cx 14 |
| 21 | Salié S(a, b; N) | iff factoring | **FACT** (precision N⁻⁴ PROVED; δ = O(1/log N) heuristic) | **Prop SAL**; cx 16 (precision N⁻⁴, one query) |
| 22 | Kloosterman K(a, b; N²) | iff factoring | **FACT** (as row 21) | **Prop SAL** (square case) |
| 23 | Kloosterman K(a, b; N), N squarefree | no known algorithm, even for prime N | NC; K(a,b;pq) = K(a q̄, b q̄; p)K(a p̄, b p̄; q) | rank-2 Kloosterman sheaf (Deligne, SGA 4½); "no bounded-genus curve computes it" is HEURISTIC; cx 15; elliptic exceptions F_{2^m} (Lachaud–Wolfmann 1990) and F_{3^m} (Katz–Livné 1989); classical only: Bruin (arXiv 1807.03600) has a quantum state-preparation algorithm, while a single classical value costs ≈ √q |
| 24 | Cubic Weyl sums Σ e(ax³/N); cubic Gauss/Jacobi sums mod N | no | FACT-hard on the promise p ≡ q ≡ 1 (mod 3) with selectable queries (analytic C3); NC in general; zero-test decides (−3/p) | analytic §0 item 3; cx 17 |
| 25 | Non-quadratic characters of (Z/N)^* (and their Gauss/Jacobi sums) | no (needs a prime ideal above p) | NC; rational arguments make (·/(N))_n trivial | Thm RA remark |

### Half-integral weight / theta

| # | Invariant | Poly-time from N? | Verdict | Ref |
|---|---|---|---|---|
| 26 | r₂(N) = 16[ε = 1] | no | NC (open, = ε) | factorbits §2.3; cx 19 |
| 27 | r₃(N) = 12h(−4N) etc.; Hurwitz H(N); ternary representation numbers | no (class numbers, L[1/2]) | FACT-hard direction (h(−4N) ⇒ ambiguous forms); family Q/M | Gauss; Shanks; cx 20 |
| 28 | r₄(N) = 8(N+1+p+q) | iff factoring | FACT | Jacobi; cx 21; nonsmooth |
| 29 | Half-integral weight cusp coefficients c(N) (Waldspurger c(N)² ∝ L(½, f⊗χ_N)) | no (conductor ≍ N², so ≈ N terms) | NC; not local at p | Waldspurger 1981; Kohnen–Zagier 1981 |
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
| 48 | Bost–Connes polylogs Li_k(e(a/N)) | yes (numerically) | only the symmetrisation is elementary: Li_k(e(x)) + (−1)^k Li_k(e(−x)) = −(2πi)^k B_k(x)/k!; a function of a/N, no factor structure | — |
| 49 | Barvinok counts of N-dilated fixed-normal polytopes | yes | AGG + height (Ehrhart quasi-polynomial: polynomial in N with coefficients periodic in N) | Ehrhart; Barvinok 1994 |
| 50 | Same, with variable normals | yes | E / lattice; divisor detection needs a hull with Θ(N^{1/3}) vertices | analytic §0 |
| 51 | Rademacher symbol / ζ_C(0) on forms of discriminant t² − 4 with N \| t² − 4 (Richaud–Degert) | yes | family Q (genus characters via relations); handed to `genusleak` | Zagier 1975 |
| 52 | Partition function p(N) | Õ(N^{1/2}) | NC; no multiplicative structure | Rademacher; Johansson 2012 |
| 53 | Fixed finite automata / AC⁰ on the bits of N (cx_2) | yes | OPEN as a class (factorbits §3.3 F; leakmap's) | selfred §0 item 5 |
| 54 | Central binomial residue binom(2N, N) mod N (cx_2) | no known fast evaluator | FACT-hard on q < 3p/2: gcd = q (Lucas-type congruences; checked by codex on 74 pairs) | codex_2 #17 |

**Reading the table.**
- **Computable, AGG (normalised; possibly plus explicit height terms):** rows 1–3, 5, 10, 11, 16, 19, 20, 31, 32, 41,
  49. Row 17 is AGG in the proved example and heuristic in general.
- **Computable, class E:** rows 6, 8, 12–14, 33, 40.
- **Computable, OPEN, outside E (not tested here):**
  - row 9 (lattice-cone sums);
  - row 15 (symbols at moving arguments);
  - row 42 (N-dependent exponents; leakmap);
  - row 50 (variable-normal polytopes);
  - row 53 (automata on bits; leakmap).
- **FACT-equivalent:** rows 7, 21, 22 (root-type closed forms), 28 and 37. Row 24 and row 54 are FACT-hard on a
  promise.
- **No polynomial algorithm known from N alone:**
  - easy given p, q: rows 4, 18, 26, 43, 45;
  - not known to be easy even given p, q (they need class numbers, L-values, p−1 factorisations or long sums): rows
    23, 27, 29, 30, 34–36, 44, 46, 47, 52.
- Row 51 belongs to genusleak's family Q.

---

## 2. Theorems and propositions

### 2.1 Theorem RA (reciprocity aggregate). PROVED; a standard consequence of class field theory

The Artin map is compatible with the transfer: ψ_{L/K}(𝔞O_K) = Ver(ψ_{L′/Q}(𝔞)) (Neukirch, *Algebraic Number
Theory* VI, functoriality of the reciprocity map). We claim no novelty. The point is only its consequence for the
fiber floor.

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
- So Theorem PT-m applies, including its split bounds and the PT2 cross-population qualification.
- Caveat (codex_2). PT-m's fiber floor is the best prediction from the *hidden* fiber point p mod m, not from N mod m.
  - With 4 | m the fiber point determines ε.
  - Example: a pure node testing C = ∏_{r ≡ 1 (4), r in the window} r decides ε perfectly on N ≡ 1 (mod 4) without
    splitting, and splits the mixed population (PT2). It has exponential cost.
  - So Theorem L applies only to the vector of *symbol answers*. Those are functions of N mod m and carry o(1)
    advantage under GRH for m ≤ X^{1−η}.
  - The zero-test part is governed by PT-m/PT2 as before.
- A short straight-line program can denote an enormous c_v, and hence an enormous conductor. Compilation into PT-m is
  not a small-conductor statement.
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
  12N s ≡ 0 (mod 3) when 3 ∤ N (all tested), the value 12N s mod 24N is a function of (a + ā mod N, (a/N), N mod 24),
  all easily computed.
- Equivalently, the Rademacher integer Φ = (a + ā)/N − 12s(a, N) = −Σ(−1)^{i+1}a_i + 2 − (−1)^r carries everything
  else.

**DSF. Separating the local part is factoring.**
- Statement: for any a with gcd(2a(a²+1), N) = 1, the exact value L = s(a, p) + s(a, q) determines p + q, and hence
  the factors (one value; deterministic). Conversely the factors give L by reciprocity at p and q.
- *Proof.*
  1. k_p := 6p·s(a, p) ∈ Z, and 2k_p ≡ a + a⁻¹ (mod p) (the Rademacher congruence at the prime p; it holds at p = 3 too, as the critic checked). Same for q.
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
    proof; it is Iwaniec–Kowalski Lemma 12.4 for any odd modulus coprime to a). Tested numerically.
  - (ii) K(c², 1; N²) = N·T_{N²}(c). This is the evaluation of Kloosterman sums at even prime powers (Iwaniec–Kowalski
    §12.3) combined with twisted multiplicativity. Tested numerically for N = 15, 21, 35.
  - (iii) T_M(c) = 2cos(4πc/M) + 2cos(4πce/M).
  - *Direct proof of (ii)* (codex_2):
    1. Write x = y + Nt with t mod N. The sum over t vanishes unless a − b y⁻² ≡ 0 (mod N), and otherwise
       contributes N.
    2. Each stationary class has a unique lift with ax² ≡ b (mod N²), where the phase is 2ax.
    3. Put v = ax. No Gauss factor arises, and none is needed for p ≡ 3 (mod 4).
- **Decoder.** Query T_M(2^k c) with absolute error δ for k = 0, …, K = bitlen(M) + 2.
  1. Subtracting the known first cosine leaves d_k = ‖x_k‖ with x_k = 2^{k+1}ce/M mod 1.
  2. Reconstruct x₀ backwards: x_k ∈ {x_{k+1}/2, (x_{k+1}+1)/2}, chosen by ‖x_k‖ = d_k. The error halves at each
     step.
  3. Round to a multiple of 1/M to get ±2ce mod M.
  4. Then gcd(ce − c, N) = p, since ce ≡ c (mod p) and ce ≡ −c (mod q).
- **Correctness when it succeeds:** PROVED (the output is gcd-verified).
- **Success probability:** HEURISTIC.
  - A step fails only when x_{k+1} lies within O(δ) of ½. There are about log₂M steps, so success is roughly
    exp(−cδ·log N), and δ = O(1/log N) is needed. It is not a constant.
  - Measured (successes out of 20; the critic's runs at δ = 10⁻² gave 22/30 at n = 512 and 9/20 at n = 1024, and at
    δ = 10⁻³ gave 18/20 at n = 1024):

| n | δ = 10⁻² | 3·10⁻³ | 10⁻³ | 10⁻⁴ | 10⁻⁶ |
|---|---|---|---|---|---|
| 64 | 20/20 | 20 | 20 | 20 | 20 |
| 128 | 19/20 | 20 | 20 | 20 | 20 |
| 256 | 15/20 | 18 | 18 | 19 | 20 |

  (Queries: 67, 131, 259; `newfamily_salie_sweep.json`. Kloosterman-at-N² version tested at n = 40.)
- cx 16 gives the complementary single-query decoder at precision N⁻⁴. That decoder is the basis of the PROVED FACT
  claim: one query at precision N⁻⁴ recovers min(v, N − v) by rounding.
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

**Random arguments.** Take a bounded prediction rule g(N, a), with a a fresh uniform unit independent of p, q given N.
Then E[ε·g(N, a)] = E[ε·G(N)] with G(N) = E_a[g(N, a) | N]. This preserves the linear prediction objective of that
rule, not Pearson correlation of a raw feature.
- A conditional floor needs a uniform bound |G(N) − B(N)| ≤ CN^{−δ}, for a baseline B with an independently
  established factor-bit bound (e.g. a function of N's size and N mod m).
- A limiting distribution alone does not suffice: residue-class postprocessing (parities, `dmod7`) can distinguish
  distributions with the same continuous limit.
- The battery's `hash` argument is deterministic, not a fresh unit, so this reduction is a heuristic model for it.
- Bettin–Drappeau (arXiv 1903.00457) average over all denominators q ≤ Q, so they are not the tool for a fixed
  modulus. Fixed-modulus results exist:
  - Porter (Mathematika 22, 1975), improved by Ustinov (2008) and Bykovskii–Frolenkov (Mat. Sb. 2017): the mean
    Euclid length over a ∈ (Z/N)^* is (12 ln 2/π²) ln N + C + O(N^{−1/6+ε}), uniformly in N. So the
    factorisation-uniform statement holds for the *mean* of `len` (KNOWN).
  - Aistleitner–Borda–Hauke (arXiv 2210.14095): concentration and tail bounds for partial-quotient sums and Dedekind
    sums at fixed N.
  - Mean squares of s(a, N) over units (Walum 1982; Conrey–Fransen–Klein–Scott, JNT 1996; Zhang, JTNB 1996). Their
    main terms carry explicit ∏_{p|N} factors, i.e. the factorisation enters at relative size ≈ 1/p (exact form
    UNVERIFIED by us). This parallels DSF and still gives only o(1) advantage.
- David–Kim–Mor–Shapira (2025, Thm 1.2) give polynomial-rate large-deviation bounds for normalised CF length and
  fixed-word frequencies for every sufficiently large denominator.
- Status: the random-argument floor is KNOWN-in-substance for the length mean. It is OPEN for general class-E
  characteristic functions.

**Structured arguments.** For structured a(N) (⌊√N⌋, ⌊N^{2/3}⌋, Fermat residues, 2^{(N−1)/2}) we have no theorem.
The CF of a/N gives congruences a·q_k ≡ ±r_{k+1} (mod N) with q_k r_{k+1} ≍ N/a_{k+1}.
- These are small only at very large partial quotients, and no quadratic relation (a square times a known small
  number) is produced.
- So we see no analogue of the relation mechanism that made the CF of √N leak at n = 32 (eps §5). This is a heuristic
  observation, not a theorem. This is the registered prediction P1,
The battery null is consistent with it, though a null cannot confirm a mechanism.

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
- All size × target cells are valid by the registered rules. V1: every |z − expected| ≤ 1.2. (V1 uses one fixed flip pattern per
  target, which is why it reads 21.2 at every n. It checks the plumbing, not sampling variability.) V2: 0 violations in
  every cell; at n ≤ 64 the PCs fire with z ≥ 19.8. V3: shuffle mean z² ∈ [0.82, 1.22].
- `nc_N_mod8` and `nc_N_mod256` are exactly degenerate everywhere (P2).
- Also degenerate, as anticipated in the prereg: `ferm_jac`, `euler_jac`, and the size statistics of `sqrt`/`cbrt2`
  whose first partial quotient fixes them. `nc_ded5_sign`, `ferm_dsign`, `ferm_altsign` are degenerate too (constant
  sign). Degenerate statistics detect nothing.

**Conclusion (registered rule).**
- **Refuted for this menu** at n ∈ {48, 64, 128} for both targets.
  - No confirmed hit. The design had 80% power at |z| = 6 to detect a *linear* association of the integer-coded
    statistic, beyond the N mod 2⁸ floor, at |ρ| ≥ 0.0153 (ε) or 0.0217 (QR₂).
  - This is a power statement, not an exclusion.
  - A non-monotone dependence of the categorical codes (`dmod7`, `dmod5`, `ms_mod3`, `nc_cubic7`, `rem*`) would be
    missed by this linear test.
- n = 32 also has no hit, but see the two footnotes.
- P1, P2 and P3 were all met.

\* **n = 32 inflation.** Only 1 491 primes qualify per side, so each prime appears in ≈ 270 instances, and the
permutation z are anti-conservative (eps §5.3). We did not rerun the cluster-robust variance. No |z| reached 6, so the
conclusion does not depend on it.

† **`nc_ded5_bits` (z = −5.69) and `nc_lead16` (−4.97) at n = 32, ε.** These are size statistics, NC by design. They
reflect a real association under the finite D₃₂ law. Among the 1 491 sixteen-bit primes with the top two bits set,
the share ≡ 1 (mod 4) is 0.504/0.503/0.515/0.477 in the four quarters by size (corr(size, [p ≡ 1 mod 4]) = −0.0073).
- codex_2 enumerated all 555 025 admissible pairs. The population-implied z at 2·10⁵ draws are −4.62 (`nc_lead16`)
  and −5.47 (`nc_ded5_bits`), against observed −4.97 and −5.69.
- Prime reuse does not create dependence between independent draws from a fixed population. It only makes this
  genuine finite-population association measurable.
- "Chebyshev-type" is an analogy, not a proof of the precise imbalance. It is the same kind of effect as the QR₂ finite-set offset in eps §3. Post-hoc
explanation; |z| < 6, so not a hit. Both are null at n ≥ 48 (|z| ≤ 2.2).

**Nonunits (rule 3).**
- All 12 instances with gcd(a(N), N) > 1 are at n = 32: 6 from `cbrt2` and 6 from `hash`. The ACCIDENTAL null
  expectation Σ(1/p + 1/q) is 7.0 per argument, so both match.
- Three arguments are structurally coprime on D_n, and one is coprime by observation:
  - `sqrt`: ⌊√N⌋ ∈ (p, q) for balanced N. Twin primes give ⌊√N⌋ = p, but then N ≡ 3 (mod 4), which D_n excludes;
  - `ferm`: p | (r+1)² would need r + 1 ∈ {p, q};
  - `euler`: a unit;
  - `gold` (not structural): ⌊αN⌋ ≡ ⌊p{αq}⌋ (mod p), which vanishes only if {αq} < 1/p.
    - Exact enumeration of the whole D₃₂ population (555 025 pairs) gives 0 such pairs, against a uniform-residue
      benchmark of 7.0 per 2·10⁵. This is a finite-window effect of the golden ratio's three-distance structure.
    - It is not a theorem: codex_2's counterexample 193·233 = 44 969 has ⌊αN⌋ = 193·144.
- No STRUCTURAL splits; n ≥ 48 has none.

---

## 4. The most promising open candidate, and why

**Candidate.** Class-E functionals at structured N-dependent arguments, sharpened by Prop DSF. Concretely: the
Rademacher integer Φ(a(N)/N), or the level-11 modular symbol, at a(N) ∈ {⌊√N⌋, ⌊N^{2/3}⌋, (⌊√N⌋+1)² − N, …}.

**Why this one.**
- (1) It is the main class of poly-time, reciprocity-computed invariants that Theorem RA does not reduce to N mod m
  plus size (row 9's lattice-cone sums are a higher-dimensional relative). RA is multiplicative and this class is additive.
  - We have *not* proved that it is not such a function at the structured arguments. "Not an aggregate" is literally
    true only because of the size dependence.
- (2) DSF shows that one natural split, s(a, N) = s(a, p) + s(a, q) + U(a, N), is exactly as hard as factoring, with
  U hiding φ(N) mod N. That is a statement about hardness of separation, not evidence that s(a, N) is closer to
  leaking than (a/N) = (a/p)(a/q).
- (3) No barrier is known. PT does not apply (Euclid on integer representatives is not generic-ring). RA does not
  apply (additive, not multiplicative). The relation mechanism is absent (§2.5).

**Why it is probably not a fourth family.**
- The battery is null at |ρ| ≥ 0.015.
- For random arguments a floor would follow from an effective, factorisation-uniform equidistribution theorem for CF
  statistics in (Z/N)^* (§2.5).
- DSF shows the local part is exactly as hard as factoring. A leak would have to come from a *partial* separation
  that the CF path performs by itself, and we see no mechanism for that.

**Sharpest question.** Is there a polynomial-time functional of the CF path of a/N (a = a(N) public) whose correlation
with the "local" Dedekind value s(a, p) + s(a, q) exceeds that of the trivial predictor s(a, N)? A related, not proved equivalent, question: does the unit part
U(a, N) have a computable component correlated with φ(N) mod small moduli?

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
| Salié / Kloosterman at N² | approximate-value oracle | ≡ FACT at precision N⁻⁴; dyadic decoder at δ = O(1/log N) (SAL) | FACT (PROVED at N⁻⁴; decoder heuristic, measured) |
| Class E at N-dependent a | 84 statistics, 4 sizes | no |z| ≥ 6; nominal 80% power for linear effects at \|ρ\| ≥ 0.0153 | refuted for this menu (registered wording; power, not exclusion) |
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

**Adversarial round (newfamily_2, xhigh): MAJOR REVISION; all items applied.**

Confirmed correct:
- RA (no extra discriminant condition needed);
- DS1, DS2, DSF (28 664 checks including p = 3) and the U formula;
- the SAL closed forms and the scalable-precision reduction;
- battery execution and registered validity (40 raw records recomputed; frozen modules unchanged).

Fixed:
- **MAJOR.**
  - (a) RA-PT wrongly chained Theorem L (an observation floor) onto PT-m's hidden fiber floor. The C = ∏ r
    counterexample is noted; L now applies only to the symbol answers.
  - (b) The constant-precision decoder is HEURISTIC throughout. The PROVED claim rests on precision N⁻⁴.
  - (c) The random-argument reduction is restated for a bounded predictor, with a baseline B(N) and the limit-law
    caveat. `hash` is noted to be deterministic.
  - (d) MDE language is restated as power, not exclusion.
  - (e) The table and verdict no longer claim an exhaustive AGG/E split. Rows 9, 15, 42, 50 and 53 are separate open
    bins, and row 24 is FACT-hard only on a promise.
- **MINOR.**
  - The class-E definition now requires polynomial output size; "never small" is heuristic.
  - David–Kim–Mor–Shapira 2025 is cited.
  - The n = 32 size association is a real finite-population effect (population-implied z −4.62/−5.47); clustering is
    not the explanation.
  - The golden-ratio "structural coprimality" was wrong (193·233). The D₃₂ count of 0 is a finite-window fact,
    verified by exact enumeration.
  - Row 48 is corrected to the symmetrised Bernoulli identity.
  - Row 49 now carries height terms.
  - Rows 53 and 54 are added.
  - The DSF hypothesis is kept in all summaries, and the decoder handles gcd = N.

## 8. Opus critics

**Proof critic: ACCEPT WITH FIXES.** All fixes are applied above.

Verified independently (scripts in `/private/tmp/claude-501/newfamily_critic/`):
- DS1, DS2 and DSF: 2 036/2 036 cases, including p = 3.
- The U-residue formula.
- The Salié closed form at N, including p ≡ q ≡ 3 (mod 4), and the Kloosterman closed form at N².
- The RA proof and the cubic test.
- The battery implementation against the prereg, and operational isolation.
- The n = 32 size artefact: the exact population correlation of log N with ε predicts z ≈ −4.1.

Fixes applied:
- **MAJOR 1.** "Constant precision" was false for the dyadic decoder, which needs δ = O(1/log N); the critic
  measured 9/20 at n = 1024, δ = 10⁻². The PROVED claim is now based on cx 16's N⁻⁴ single query.
- **MAJOR 2.** The §0 Jager–Schwenk sentence overclaimed. It is now restricted to N-independent arguments with
  conductor ≤ X^{1−η}.
- **MINOR.**
  - §4's "provably not an aggregate" and "more than any other object" are softened.
  - The n = 32 max |z| is now reported (5.69, NC).
  - "Excluded" is replaced by a power statement, noting the linear-test limitation for categorical codes.
  - RA is cited as a standard transfer consequence.
  - V1 is noted to be a single flip pattern.
- **NIT.**
  - The golden-ratio inequality is corrected.
  - The twin-prime remark is added.
  - `factor_from_local_part` no longer crashes when N | a² + 1.
  - The DSF congruence also holds at p = 3.
  - The aggregate wording for a + ā is fixed.
  - "Confirms" is removed.

**Literature critic: ACCEPT WITH FIXES.** All fixes are applied.

- **MAJOR.** §2.5 cited Bettin–Drappeau, which treat averages over all denominators. It now cites the fixed-modulus
  literature: Porter, Ustinov and Bykovskii–Frolenkov for the uniform Euclid-length mean; Aistleitner–Borda–Hauke;
  and the mean-square results for Dedekind sums.
- **MINOR.**
  - Salié (i) is now cited as Iwaniec–Kowalski Lemma 12.4.
  - SAL is a routine composition of the closed forms with Rabin 1979 (a nontrivial square root of 1 factors N). No
    prior formulation was found.
  - Row 23: the sheaf is Deligne's; Katz–Livné is added for F_{3^m}; "classical" is specified and Bruin's quantum
    algorithm noted.
  - RA is now presented as standard (cf. Lemmermeyer). The contribution is Corollary RA-PT.
  - Row 18 now gives both Rédei obstructions.
- **NIT.**
  - Waldspurger term count is corrected (≈ N terms).
  - References: Rademacher–Grosswald, *Dedekind Sums* (Carus 1972) ch. 3 for the congruences; Barkan 1977,
    Hickerson 1977 (Crelle 290) and Knuth 1977 (Acta Arith. 33) for BHK; Kirby–Melvin 1994 and Ghys 2007 for Φ;
    Garoufalidis–Zagier arXiv 2111.06645; Xin–Xu arXiv 2303.01185 (Exp. Math. 2024); Johansson 2012 for p(N);
    Miller 1976 is deterministic under ERH and randomized unconditionally.
- **Not verified by either critic:** Garoufalidis–Zagier Prop. 10.1 (cited by cx 22).

## 9. Reproduction (project root)

```
.venv/bin/python -m pytest -q tests/test_newfamily.py
PYTHONPATH=. .venv/bin/python -m src.newfamily.harness --n 48 64 128 32 --range main --count 200000 --workers 6
PYTHONPATH=. .venv/bin/python -m src.newfamily.sweep          # Prop SAL precision sweep (~3 min)
```
