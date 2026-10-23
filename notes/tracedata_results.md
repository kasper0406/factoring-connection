# tracedata: trace data a_p, a_q, a_N, fixed-ℓ partial counts, X_E and determinism (phase 25)

Agent `tracedata`, 2026-10-23. Version **v1** (before reviews; §9 will record them).

- Pre-registration: `notes/tracedata_prereg.md`, commit `f2df008`, made before any experiment.
- Code: OPERATIONAL `src/tracedata/{reductions,decode,xe_descent,galois}.py` (inputs: N and the stated data only; an
  AST test forbids harness/PARI/factorisation imports and calls); HARNESS `src/tracedata/{harness,run_experiments}.py`.
- Tests: `tests/test_tracedata.py` (20 pass, < 1 s).
- Raw data: `experiments/raw/tracedata_{e1,e2,e3a,e3b,e3c,e3d,e4,e5,posthoc}.jsonl`.
- GPT-6.1-sol blind round: `codex/tracedata_1.{md,out}`, cited **cx1**.
- Labels: KNOWN (with reference), PROVED, MODEL (exact computation or proof inside the Chebotarev/Haar model),
  HEURISTIC, experimentally supported, REFUTED (prediction), OPEN, UNVERIFIED. Strongest novelty wording: "no prior
  statement found by us".

Notation: N = pq, p < q, gcd(N, 6Δ_E) = 1; E/Q fixed; a_r = r + 1 − #E(F_r); A = a_N = a_pa_q; T = a_p + a_q;
s = p + q; M = #E(Z/N) = (p+1−a_p)(q+1−a_q); S = qa_p² + pa_q²; g_r = ρ_{E,ℓ}(Frob_r) ∈ GL₂(F_ℓ); t_r = tr g_r,
d_r = det g_r = r mod ℓ. "X ⇒ FACT" means a polynomial-time algorithm factoring N from N and X.

---

## 0. Verdict (plain language)

**Neither the trace pair (a_p, a_q) of one fixed curve nor its Hecke eigenvalue a_N is known to factor N. We found no
reduction, and we found four precise reasons why the natural routes fail.** Bach–Charles's extra coefficient
a_{N²} works because it carries the determinant of Frobenius, i.e. p itself.
1. **Squarefree Hecke data are determinant-free.** a_N, and (a_p, a_q), are polynomials in traces. For non-CM E no
   algebraic relation ties them to p on a positive-density set (Prop DEC). Any one datum that brings in the
   determinant closes the algebra: a_{N²} = a_N(Sym²), a_{N³}, a weight-3 CM partner, or M. Each then factors N
   deterministically (§2.2–2.3). DEC excludes only algebraic recovery. M alone is also algebraically free of p, and
   it factors N through group structure (TW*), so DEC is not hardness evidence.
2. **One curve's congruences carry about 0.3 bits.** In the Chebotarev/Haar model, all congruences a_p, a_q mod ℓ^k
   together tell about 0.3 bits about the residues of p (Prop GC). CRT needs n/2 bits.
3. **Trace products see only quadratic residuosity (Theorem SC, MODEL; no prior statement found by us).** For ANY
   number of curves, the mod-m Galois content of (N, a_N(E_1), …, a_N(E_K)) is invariant under
   (p, q) ↦ (λ²p, λ⁻²q). So congruence methods applied to Hecke eigenvalues at N can reach only fixed-base QR bits,
   never factors. E3b shows exactly this plateau, and E3d fails 0/30.
4. **Promise obstruction (cx1).** For y² = x³ − x on Blum moduli, all trace data are 0 or powers of N, computable
   without factoring. So an everywhere trace-only reduction for this curve would factor Blum integers outright.
   Doubly supersingular N play the same role for every curve.

**Positive (HEURISTIC, strongly supported).** Ordered trace NUMBERS of polynomially many fixed curves do factor N:
- In the Haar model, (a_p(E_i))_{i<K} determines p mod ℓ once K ≈ 1.4ℓ²log₂ℓ (one-trace information ≈ 0.72/ℓ² bits).
  CRT then gives p.
- E3a: 40/40 cells match the model. E3c: 30/30 toy N factored from 40000 traces.
- Unordered per-curve pairs and multi-curve M mod ℓ identify s mod ℓ instead, so they also factor heuristically.
- So FROB's "numbers, not operators" gap closes heuristically for trace pairs at many curves. It does not close for
  a_N (Theorem SC), and it does not close for one curve.

**Fixed-ℓ counts** (M mod ℓ, a_N mod ℓ, root counts of a fixed polynomial):
- Each has O(1) bits. Under N-local reductions none is factoring-hard unless factoring is in BPP (Prop IL).
- They sit above fixed-base QR: r_N(x² − Δ) ≡ QR_N(Δ); QR_N(Δ_f) ≤ r_N(f) for S₃ cubics (ecount); a_N for a CM curve
  decides QR_N(−D) exactly.
- They carry non-abelian information that no finite set of quadratic symbols determines (cx1, PROVED).
- In Theorem PT's pure-tree model, computing any of them exactly forces splitting at rate ≥ 0.16–0.36.

**X_E.**
- Everywhere randomized reductions now cover j = 0 (cx1: sextic twists, PROVED), j = 1728 and rational 2-torsion
  (X_E = ∅).
- For other j it is OPEN. Every twist group over every quadratic étale Z/N-algebra is annihilated by M on X_E (cx1),
  so twist-type methods are blind there.
- A deterministic one-curve reduction is OPEN. Bach's least-non-residue bound does not apply, because x ↦ (f(x)/r) is
  not multiplicative.

---

## 1. Equivalence table (oracle relations given N; D deterministic, R randomized, ≤ poly-time Turing)

| Data (fixed curve(s)) | Relation to FACT | Status / source |
|---|---|---|
| M = #E(Z/N) | ≡_R FACT off X_E; everywhere if E has rational 2-torsion, j = 1728, or **j = 0** | ecount TW*; Theorem J0 (cx1) §2.9 |
| (a_N, a_{N²}) | ⇒_D FACT unless a_p = a_q = 0 | KNOWN, Bach–Charles; E1a 450/450 as predicted |
| (a_N, a_{N³}) | ⇒_D FACT iff a_N ≠ 0 | PROVED §2.2; E1b 450/450 as predicted |
| (a_N(f), a_N(g)), a_p(g) = c·a_p(f)² − e·p | ⇒_D FACT (Sym²: c = e = 1; CM weight 3: c = 1, e = 2) | PROVED §2.2; E1d 90/90 |
| ({a_p, a_q}, M) or (T, A, M) | ⇒_D FACT (both orders, quadratic); a = b = 1 via TW* | PROVED §2.3 (cx1 for a = b = 1); E1c 450/450 |
| {a_p, a_q} + any one datum polynomial in p with s-dependence | ⇒_D FACT | PROVED §2.3 |
| (a_p, a_q), ordered or unordered (≡ up to ×2), one curve | **OPEN**; everywhere impossible for CM unless Blum factoring is easy | §2.1, §2.4 |
| a_N, one curve | **OPEN**; ≤ (a_p, a_q) + factoring a_N (n/2 bits); ⇒ s mod m if E(Q)[m] ≠ 0; ⇒ QR_N(−D) if CM | §2.5 |
| (a_N(E_i))_{i ≤ K}, any K | Galois content = quadratic residuosity bits only | **Theorem SC** (MODEL) §2.6; E3b, E3d |
| (a_p(E_i))_{i ≤ K} ordered, K = O(n² log n) | ⇒ FACT | HEURISTIC §2.8; E3a, E3c 30/30 |
| ({a_p(E_i), a_q(E_i)})_i or (M_i mod ℓ)_{i, ℓ ≤ L} | ⇒ s mod ℓ ⇒ FACT (K ~ ℓ² resp. ℓ⁴ per ℓ) | HEURISTIC (model identifiability, PH1/PH2; untested end-to-end) |
| M mod ℓ, a_N mod ℓ, r_N(f) for fixed ℓ, f | O(1) bits: not FACT-hard under N-local reductions unless FACT ∈ BPP | Prop IL §2.7 |
| r_N(x² − Δ) | ≡_D QR_N(Δ) | PROVED (cx1) |
| r_N(f), f an S₃ cubic | QR_N(Δ_f) ≤_D r_N(f) ≤_D FACT; not a function of finitely many quadratic symbols | ecount (cx2); Prop NA (cx1) |
| M mod ℓ (one curve) vs φ(N) mod ℓ | neither determines the other at the Frobenius level | PROVED (cx1) |

---

## 2. Statements and proofs

### 2.1 Prop DEC (determinant-free data are algebraically independent of p). PROVED (non-CM)

**Claim.** Let E/Q be non-CM, and let p, q be iid uniform primes in [X, 2X] (the project's instance law). For every
nonzero G ∈ R[P, Q, U, V],

  Pr[G(p, q, a_p, a_q) = 0] → 0  (X → ∞).

In particular, for any F(N, U, V, Z) with F(PQ, U, V, P) ≢ 0 (i.e. the relation is not a formal consequence of N =
PQ), Pr[F(N, a_p, a_q, p) = 0] → 0. So p, q, s are not algebraic functions of (N, a_p, a_q) on any event of
non-vanishing probability.

*Proof.*
- Write p = uX, q = vX, a_p = 2x√p and a_q = 2y√q, with (u, v, x, y) ∈ [1, 2]² × [−1, 1]².
- Give P and Q weight 2 and U and V weight 1, and let w be the top weight among the monomials of G. Then

  G = X^{w/2}·(Φ(u, v, x, y) + O_G(X^{−1/2})),

  where Φ = Σ_{top weight} g_{ijkl} u^{i+k/2} v^{j+l/2} (2x)^k (2y)^l.
- Φ is a nonzero real-analytic function: distinct monomials give distinct functions. So its zero set has Lebesgue
  measure 0.
- G = 0 forces |Φ| ≤ C_G X^{−1/2}.
- By Sato–Tate for non-CM E/Q (Barnet-Lamb–Geraghty–Harris–Taylor, Publ. RIMS 47 (2011)), applied in the windows
  [uX, u′X], the law of (u, x) converges weakly to (uniform on [1, 2], by the prime number theorem) ×
  μ_ST. μ_ST is absolutely continuous. The same holds for (v, y), and the two are independent.
- Portmanteau on the closed set {|Φ| ≤ η} gives limsup_X Pr[G = 0] ≤ μ(|Φ| ≤ η), which tends to 0 as η → 0. ∎

**CM curves.** The claim fails: G = U vanishes on the probability-½ event "p inert", because the CM
Sato–Tate law has an atom at x = 0. Conditioned on split primes, the angles are equidistributed (Hecke) and the proof
goes through.

**Scope (important).**
- DEC excludes algebraic recovery: straight-line programs over Q, radicals and root selection, Coppersmith-type
  polynomial relations ("a_pa_q = A" contains no p, so there is no small-root equation at all).
- It does not exclude arithmetic: gcds, group structure, Galois congruences, divisibility, or factoring a_N.
- M = (p+1−a_p)(q+1−a_q) is ALSO algebraically independent of p over Q(N, M), yet M factors N through the group law
  (TW*). So DEC is not evidence of hardness. It locates the missing ingredient: the determinant.
- **Polynomial annihilators.** No nonzero K ∈ Z[N, X, Y] has (P + 1 − X) | K(PQ, X, Y) in Z[P, Q, X, Y]. Substituting
  X = P + 1 gives K(PQ, P+1, Y) = 0, and (PQ, P+1, Y) are algebraically independent.
- So any multiple of #E(F_p) that the traces could supply cannot come from an identity. It must be smoothness-type
  (ECM constants) or arithmetic.

### 2.2 Bach–Charles and what a_{N²} supplies. KNOWN + PROVED variants

- Hecke at p²: a_{p²} = a_p² − p. The determinant det ρ(Frob_p) = p enters here and nowhere in squarefree a_N.
- a_{N²} = a_N² − S + N with S = qa_p² + pa_q², and the roots of Z² − SZ + N·A² are qa_p² and pa_q².
  Their gcds with N give q and p when a_pa_q ≠ 0. If exactly one trace vanishes, gcd(S, N) splits N. If both vanish,
  the data are (0, N) and carry nothing (KNOWN form: Bach–Charles, arXiv:0708.1192 v2, Thm 2.2; weight-2 form ecount
  §3).
- Interpretation: a_{N²} = a_N(Sym²f), the Hecke eigenvalue at the squarefree N of the GL₃ form Sym²f, since the trace
  of Sym²ρ(Frob_p) is a_p² − p. Bach–Charles use two "trace numbers at N", one of which contains det.
- **Prop BC3 (PROVED).** a_{p³} = a_p³ − 2pa_p, so a_{N³} = A(A² − 2S + 4N). For A ≠ 0 this gives S, hence the
  factors. A = 0 kills all odd powers.
- **Prop LINK (PROVED).** If a_p(g) = c·a_p(f)² − e·p at both primes (c, e ≠ 0 fixed), then
  a_N(g) = c²A² − ceS + e²N, which gives S.
  - Instance: f the weight-2 CM form of 32a2 (y² = x³ − x), g = η(4z)⁶ ∈ S₃(Γ₀(16), χ₋₄), with
    a_p(g) = π² + π̄² = a_p(f)² − 2p.
  - So the eigenvalues at N of two FIXED CM forms factor N deterministically on split–split N.
  - E1d: 90/90. Checked against the η-product q-expansion in the tests.
  - No prior statement found by us. It is an immediate variant of Bach–Charles.
- **Single a(N^k):** OPEN. cx1: on Blum N for 32a2, a_{N^{2m}} = N^m and a_{N^{2m+1}} = 0, so no single higher
  coefficient suffices everywhere for CM curves.

### 2.3 Prop MIX (traces + one coupling datum). PROVED

- Given {a_p, a_q} (equivalently (T, A): roots of Z² − TZ + A) and any datum Φ(p, q, a_p, a_q) such that
  P^k·Φ(P, N/P, a, b) is a non-constant polynomial in P, the integer roots for both assignments give p.
  This is deterministic and needs no balance.
- Instance M: with u = 1 − a_p, v = 1 − a_q,

  v·p² + (N + uv − M)·p + uN = 0,

  equivalently C = M − N − 1 + T − A = (1 − a_q)p + (1 − a_p)q (cx1).
- Degenerate a_p = a_q = 1: then M = N. cx1 shows N ∉ X_E, since twist orders p+2 and q+2 cannot both divide pq. So TW*
  applies.
- E1c: 450/450.
- **Minimal symmetric data.** With traces only, (T, A) is complete symmetric trace data, and its power is OPEN. ONE
  further datum carrying p closes it. There is no canonical "minimum" among incomparable oracles (cx1).

### 2.4 Prop PROMISE (cx1). PROVED

- For E: y² = x³ − x and p, q ≡ 3 (mod 4): a_p = a_q = 0 (f is odd and (−1/r) = −1), so a_{N^k} ∈ {0, N^{k/2}}.
- These values are computable from N without factoring. So an everywhere reduction "(a_p, a_q) ⇒ FACT", or
  "(a_N, a_{N^k}) ⇒ FACT", for this curve would factor every Blum integer with no oracle.
- For any E (CM or not), the same holds on the promise "both primes supersingular". That set is density zero for non-CM
  E, so the obstruction is conditional on factoring being hard there.
- Consequence: trace-only reductions can at best be almost-everywhere. In contrast, M on doubly supersingular N equals
  (p+1)(q+1), which factors N (SR2).

### 2.5 a_N alone: what it gives. PROVED pieces; main question OPEN

- **Hasse box (cx1).** Given A ≠ 0, the integer trace pairs explaining it number ≤ 2τ(|A|). This bounds candidates,
  but identifying primes having those traces is the whole problem. a_N ⇒ (a_p, a_q) up to factoring A (n/2 bits;
  subexponential, not polynomial) and τ(A) candidates.
- **Lemma TOR.** If E(Q) has a point of order m, then m | #E(F_r) at good r ∤ m. So a_r ≡ r + 1 and
  a_N ≡ (p+1)(q+1) = N + 1 + s (mod m), which gives s mod m (= the φ(N) mod m information).
  - Ordered (a_p, a_q) give p mod m.
  - Mazur bounds m: lcm reachable = 5040 via Z/2×Z/8, Z/9, Z/5, Z/7. This is fixed-modulus information only.
  - E1e: 270/270.
- **Lemma CM-QR.** For CM by O_K, K = Q(√−D) of class number one, on (−D/N) = 1, a_N ≠ 0 iff both primes split.
  At split r > D, a_r = π + π̄ ≠ 0, since π = −π̄ would force r = D·(square)/4.
  So a_N decides QR_N(−D) exactly. For 32a2 this decides ε (whether N is a sum of two squares). E1f: 45/45.
- **ℓ = 2.** a_N is odd iff M is odd iff r_p = r_q = 0 (ecount's 2/9 one-sided certificate).
- **Congruence content (Prop GC, next).** About 0.3 bits.

### 2.6 Theorem SC (scalar-twist symmetry of trace products). MODEL (Chebotarev); no prior statement found by us

**Setting.**
- E_1, …, E_K/Q and a modulus m. G ⊆ ∏_i GL₂(Z/m) is the image of Galois under ⊕_i ρ_{E_i, m}.
- Assume G ⊇ {(λI, …, λI) : λ ∈ Λ} for a subgroup Λ ⊆ (Z/m)^×. This holds with Λ = (Z/m)^× whenever G is the full
  fibre product {(g_i) : det g_i all equal}, e.g. for independent Serre curves.
- (p, q) range over pairs of primes ∤ mΔ, with the product Chebotarev density.
- O(p, q) := (N mod m, (a_p(E_i)·a_q(E_i) mod m)_{i ≤ K}).

**Claim.** For every λ ∈ Λ and every set B of residue pairs, the density of {O = o, (p, q) mod m ∈ B} equals the
density of {O = o, (λ⁻²p, λ²q) mod m ∈ B}. Hence:
- the posterior of p mod m given O is invariant under multiplication by Λ²;
- O can identify p mod m at most up to (Z/m)^×/Λ²;
- for m = ℓ prime and Λ = F_ℓ^×, O identifies at most the bit (p/ℓ) when (N/ℓ) = 1, and nothing at all (beyond N)
  when (N/ℓ) = −1, because then p ↔ q swaps the classes.

*Proof.*
- By Chebotarev applied to the compositum (independently at p and q), the pair (Frob_p, Frob_q) is equidistributed in
  Haar measure on G × G, on conjugation-stable sets.
- The map τ_λ(g, h) = (λg, λ⁻¹h) is a bijection of G × G, because (λI)_i ∈ G and (λ⁻¹I)_i ∈ G. It commutes with
  conjugation.
- It preserves the observation: tr(λg_i)·tr(λ⁻¹h_i) = tr g_i·tr h_i, and det(λg)·det(λ⁻¹h) = det g·det h.
- It moves the hidden residues: det(λg) = λ²·det g, i.e. p ↦ λ²p and q ↦ λ⁻²q.
- The relevant sets are conjugation-stable unions of classes, and τ_λ maps one onto the other with equal Haar mass. ∎

**Remarks.**
- The symmetry is broken by integrality and the Hasse bound: λ·a_p is not a trace when λ ≠ ±1.
- So any method extracting more than QR bits from Hecke eigenvalues at N must use their archimedean/integer structure,
  e.g. factoring a_N into a_p·a_q. Congruence decoding cannot.
- **Contrast.**
  - Ordered traces t_p (or (t_p, t_q)) are not τ-invariant. They identify p mod ℓ completely.
  - Unordered pairs {t_p, t_q} and M mod ℓ = det(1 − g_p)det(1 − g_q) identify exactly {p, q} mod ℓ, i.e. s mod ℓ
    (PH1: exact identifiability classes for ℓ = 5, 7, 11, 13, all N mod ℓ; unit-tested).
- **Positive corollary (MODEL).** With many curves, the products do identify (p/ℓ) on (N/ℓ) = 1, i.e. QR_N(ℓ*). So the
  Galois content of "Hecke eigenvalues at N" is exactly fixed-base quadratic residuosity.
- **Experiment.** E3b: the product decoder's accuracy for s mod ℓ is flat in K (ℓ = 5: 0.64, 0.64, 0.64 at
  K = 10³, 10⁴, 4·10⁴), exactly the square-class ceiling, and matches the model in 15/15 cells.
  E3d (end-to-end from products) 0/30, model 0.03. The pre-registration wrongly expected s mod ℓ to be decodable
  from products; the plateau exposed the symmetry (deviation D2).

### 2.7 Fixed-ℓ partial counts

**Prop IL (instance-local bits). PROVED (trivial).**
- Definitions:
  - Let g be a function on moduli whose value on N lies in a set of size c independent of N, e.g. M mod ℓ, a_N mod ℓ,
    r_N(f) ∈ {0, …, (deg f)²}.
  - Call a reduction N-local if each oracle answer it requests is computable in polynomial time from g(N) and the
    query, e.g. queries N^j·k with k's local data computable, using multiplicativity:
    #E(Z/N^jk) = N^{j−1}M·#E(Z/k), and r_{N^j}(f) = r_N(f) for separable f.
- Claim: an N-local reduction FACT ≤ g with success σ yields a randomized factoring algorithm with success ≥ σ/c, by
  guessing g(N).
- So no fixed-ℓ count is factoring-hard under N-local reductions unless factoring is in BPP. Fixed-base QR_N(Δ) has
  the same status.
- Not covered: reductions that query g on moduli coprime to N with unknown factorisations.

**Placement.**
- QR_N(Δ) ≡_D r_N(x² − Δ) (cx1).
- QR_N(Δ_f) ≤_D r_N(f) ≤_D FACT for S₃ cubics (ecount/cx2). The root count of a fixed polynomial is computable after
  factoring via deg gcd(f, x^r − x).
- **Prop NA (cx1, PROVED).** r_N(f) for an S₃ cubic, and M mod ℓ for full image, are not functions of N together with
  finitely many quadratic (or abelian) symbols of p and q.
  - The S₃ (resp. GL₂) field meets any multiquadratic or cyclotomic field inside its maximal abelian subfield. So
    Chebotarev supplies primes with equal abelian data but Frobenius identity vs 3-cycle (resp. I vs −I:
    det(1 − I) = 0, det(1 + I) = 4).
  - This is non-determination, not a computational separation from QR oracles.
- **M mod ℓ vs φ(N) mod ℓ (cx1, PROVED).** Neither determines the other at the Frobenius level.
  - Fixed determinants allow varying M mod ℓ.
  - diag(1, d) has det(1 − g) = 0 for all d, so determinant pairs (1, 1) and (−1, −1) give the same (N, M) mod ℓ and
    different φ(N) mod ℓ.
- **Many curves (MODEL, PH1/PH2).** (M_i mod ℓ)_i over many independent curves identifies exactly s mod ℓ, i.e. the
  φ(N) mod ℓ information, with per-curve information ≈ 0.74/ℓ⁴ bits. So "fixed-ℓ counts of all curves" carry φ(N) mod ℓ
  and, over ℓ ≤ L, heuristically factor N. That is consistent with Kunihiro–Koyama's CRT reduction for arbitrary input
  curves.
- **Chebotarev predictors (MODEL, E2; census E4).** For (p/3) on N ≡ 1 (mod 3), the Bayes advantage is:

  | Predictor | Bayes advantage | Source |
  |---|---|---|
  | M mod 3 | 9/128 ≈ 0.070 | ecount |
  | a_N mod 3 | 0.156 | here |
  | r_N(ψ₃) | 0.226 | here |

  - For ψ₃, the 4 roots correspond to the 4 lines of F₃² and r_p = number of Frob-stable lines.
  - So r_p(ψ₃) = 2 ⇒ p ≡ 2 (mod 3), and r_p ∈ {1, 4} ⇒ p ≡ 1. On N, r_N ∈ {1, 16} ⇒ p ≡ 1, but r_N = 4 is
    ambiguous (2·2 vs 1·4).
  - E4: 0 exceptions to the certificate rules.
- **Theorem PT corollaries (PROVED from factorbits PT; constants MODEL).**
  - A pure tree that computes F exactly on iid primes splits with probability s ≥ (1 − max_b Pr[F = b])/2:

    | F | s ≥ |
    |---|---|
    | M mod 3 | 0.158 |
    | M mod 5 | 0.289 |
    | M mod 7 | 0.350 |
    | a_N mod 3 | 0.195 |
    | a_N mod 5 | 0.313 |
    | a_N mod 7 | 0.365 |
    | r_N(ψ₃) | 0.195 |
    | r_N(S₃ cubic) | 0.222 |

  - With the Chebotarev class law in place of μ, error o(1).
  - PT-J transfers with −1 replaced by Δ, since (Δc|x) = (Δ|x)(c|x). So Jacobi symbols of small pure values give
    negligible advantage on (Δ/p) (Landau–Page, as in leakmap PT-J1).
  - Therefore any polynomial algorithm for these counts must be non-pure: sampling, integer representatives, N as
    exponent, etc.
- **Status: OPEN** whether any of them is computable in polynomial time. A poly-time r_N(f) for an S₃ cubic would
  decide QR_N(Δ_f).

### 2.8 Multi-curve Galois decoding. MODEL constants + HEURISTIC transfer + experimentally supported

**Model.** The Haar likelihood is P(t | d) = (ℓ + χ(t² − 4d))/(ℓ² − 1). Exact information per curve in bits
(PH2, E2):

| observation | ℓ²·I at ℓ = 11 / 31 (→ limit) | what K → ∞ identifies |
|---|---|---|
| t_p (one prime, ordered) | 0.658 / 0.698 (→ 1/(2 ln 2) ≈ 0.72) | p mod ℓ |
| (t_p, t_q) ordered | 1.314 / 1.397 (→ 1/ln 2) | p mod ℓ |
| {t_p, t_q} unordered | 0.715 / 0.721 | s mod ℓ |
| t_pt_q (a_N) | ℓ³I → 0.73 | (p/ℓ) on (N/ℓ) = 1 only (Theorem SC) |
| det(1−g_p)det(1−g_q) (M) | ℓ⁴I → 0.74 | s mod ℓ |

**Prop GC (one curve; MODEL).**
- I(p mod ℓ; (t_p, t_q) | N mod ℓ) summed over odd ℓ ≤ 47 is 0.202 bits.
- The tail beyond 47 is ≈ 0.01. Prime powers add +0.013 (ℓ = 3 → 9) and +0.002 (5 → 25); the 2-adic part to mod 8
  is 0.072 (PH3).
- Under independence of the ℓ-adic images, informations add. So one curve's full congruence data carry ≈ 0.3 bits
  about p's residues, while CRT needs n/2 bits.
- Galois congruences cannot be the route from one curve's (a_p, a_q) to p.
- Scope: congruence-only use, in the Haar model. Integrality is not covered.

**HEURISTIC MC.**
- Assumption: for a fixed large prime p, the Frobenius classes of p in the ℓ-division fields of a fixed family of
  curves behave like independent Haar draws conditioned on det = p mod ℓ. This is a short-family version of
  Castryck–Hubrechts-type equidistribution, which is proved for random curves over F_p (citation UNVERIFIED).
- Under it, ML decoding recovers p mod ℓ from K ≈ 1.4ℓ²·log₂ℓ ordered traces.
- With ℓ ≤ L ≈ 0.35n (∏ℓ > p), or L ≈ 0.17n with Coppersmith on p mod W ≥ N^{1/4}, K = O(n² log n) ordered traces
  factor N.
- Unordered pairs (s mod ℓ) need K = O(n² log n); M mod ℓ needs K = O(n⁴ log n). Neither is tested end-to-end.
- Effective Chebotarev cannot justify this: the field degrees exceed p. It stays a heuristic.
- Evidence:
  - E3a: 200 primes of 30 bits, curves y² = x³ + Ax + B. Real accuracy matches the Haar simulation within 0.10 in
    40/40 cells. Mean real − sim = −0.020 over the 27 unsaturated cells, with 18/27 negative (PH4): a slight, consistent
    shortfall, plausibly from curves with non-full image.
  - The null control is at chance in 40/40 cells.
  - E3c: 30/30 toy N (28-bit primes) factored from 40000 ordered traces of one prime. The truth was top-ranked at every
    ℓ ≤ 29 in every instance.
- These are toy correctness checks of the model, not evidence of a factoring algorithm: the harness supplies the traces.

**Meaning for FROB's scope gap.**
- Trace NUMBERS with a consistent ordering, or unordered per-curve pairs, at poly many fixed curves carry factoring
  power (heuristic).
- Products a_N do not, through congruences (Theorem SC).
- One curve's pair does not, through congruences (GC).

### 2.9 X_E and determinism

- **Lemma SEP (cx1 in ecount; re-proved in cx1 here).** If λ_r and μ_r both divide an odd M, then λ_r ≠ μ_r. Some
  odd ℓ | M therefore separates curve and twist at r. Given that ℓ, descent with M/ℓ^j splits with constant
  probability.
- **Prop QA (cx1, PROVED).** On X_E, every finite étale quadratic Z/N-algebra R has E(R) annihilated by M.
  - Non-split component: P ± πP lie in E(F_r) and in the twist part, so [2M]P = O. And E(F_{r²})[2] = 0, because
    Frobenius is a 3-cycle on E[2] \ 0 when M is odd.
  - Split component: E(F_r)².
  - The same holds for the Weil restriction. So final-scalar twist methods are blind on X_E. Intermediate scalars
    (descent) are not covered.
- **Theorem J0 (cx1, PROVED; checked here).** Every fixed j = 0 curve E_b: y² = x³ + b has an everywhere randomized
  reduction M ⇒ FACT.
  - If M is even, use TW*.
  - If M is odd, both primes are ≡ 1 (mod 3). Draw uniform (x, y), set c = y² − x³, and compute [M](x, y) on
    E_c: y² = x³ + c.
  - Sextic classes: c ∈ b·(F_r^×)⁶ gives E_c ≅ E_b, all killed (mass ≥ 1/6 − O(r^{−1/2})).
  - The two classes with −c a cube give rational 2-torsion. Odd M has a kernel of index ≥ 2 there, so mass
    ≥ 1/6 − O(r^{−1/2}) is not killed. b is not in these classes, because E_b has no 2-torsion.
  - CRT then gives mixed probability ≥ 1/18 − o(1), and Lemma ZN turns a mixed infinity into a split.
  - We re-derived the class counts. cx1's exact enumeration at N = 27331 gives terminal mixed probability 0.5035.
  - No prior statement found by us or cx1.
  - Together with rational 2-torsion and j = 1728 (always 2-torsion), X_E is now empty-or-handled for every curve with
    extra automorphisms or 2-torsion.
- **General j: OPEN.** Only quadratic twists exist, and Prop QA kills them. Cubic extensions have unknown orders
  B_r = r² + r(a_r − 1) + a_r² + a_r + 1 for the trace-zero part. Isogenous curves have the same M.
- **Small-prime descent (E5).** TW* with descent over ℓ ≤ 1000 dividing M split both known X_E members at the first
  trial. These are the ecount census member for 11a1 and cx1's contrived N = 27331. ecount had already split them with
  TW*.
- **Determinism: OPEN.**
  - A deterministic version needs a small x₀ with (f(x₀)/p) ≠ (f(x₀)/q), plus a non-kernel condition.
  - x ↦ (f(x)/r) is not multiplicative (cx1: f = x³ − x − 1 mod 13, f(2), f(4) both non-squares). So Bach's ERH
    least-non-residue bound does not apply.
  - We know no GRH-type bound for least witnesses of non-multiplicative character sums of this kind.
  - Rational points on fixed twists E^{(d)}(Q) remove sampling. But they still need ord(G mod q) ∤ M for a fixed G,
    an elliptic Artin-type condition with no uniform guarantee.

---

## 3. Experiments (pre-registered in notes/tracedata_prereg.md; commit f2df008)

Peak RSS ≤ 245 MB; jobs ran one at a time.

**E1 (exact reductions; n = 32, 64, 128; 30 seeds; no instance skipped).**

| curve | BC2 ok (pred) | BC3 ok (pred) | MIX | TOR | CM-QR | both a = 0 |
|---|---|---|---|---|---|---|
| C1 = 11a1 (Z/5) | 90 (90) | 90 (90) | 90 | 90/90 | — | 0 |
| C2 = 37a1 | 90 (90) | 90 (90) | 90 | — | — | 0 |
| C3: y² = x³+2x+3 (Z/2) | 90 (90) | 90 (90) | 90 | 90/90 | — | 0 |
| CM1 = 32a2 | 62 (90 match) | 17 (90 match) | 90 | — | 45/45 | 28 |
| T7 = 26b1 (Z/7) | 90 (90) | 90 (90) | 90 | 90/90 | — | 0 |

- "pred" counts instances where success matched the prediction; every BC2 failure is exactly a_p = a_q = 0.
- CM2 (E1d): 30/30 at each n.
- **P1a–P1f all hold.**

**E2 (model constants).**
- P2a holds: ℓ²I_ord ∈ [1.31, 1.41] for ℓ ≥ 11, and 1.381–1.412 for ℓ ≥ 23, against 1/ln 2 = 1.443.
- P2b holds: 0.202 bits.
- P2c holds: ℓ³I_prod ∈ [0.730, 0.769].
- P2d holds: min floor 0.316 (M mod 3).

**E3a (ordered decoding, 200 primes × 10⁴ curves).**
- P3a holds: 40/40 cells within 0.10, max |diff| 0.0945 (ℓ = 11, K = 250).
- P3b holds: 40/40 cells with |null z| < 2.
- Selected accuracies (real / sim):

  | ℓ | K = 1000 | K = 4000 | K = 10000 |
  |---|---|---|---|
  | 13 | — | 0.995 / 0.990 | 1.0 / 1.0 |
  | 23 | 0.29 / 0.31 | 0.71 / 0.73 | 0.96 / 0.98 |
  | 31 | — | — | 0.785 / 0.826 |

**E3b (products).**
- P3c holds: 15/15. P3d holds: 15/15.
- The accuracies are flat in K, per Theorem SC: ℓ = 3: 1.0; ℓ = 5: 0.64; ℓ = 7: 0.44–0.46; ℓ = 11: 0.27–0.30;
  ℓ = 13: 0.20–0.26.

**E3c.** 30/30, P3e holds. Every ℓ had the truth ranked first.

**E3d.** 0/30 against sim 0.03, so P3f holds formally. It is uninformative: the registered design rested on a false
premise (D2).

**E4 (census, 1491 primes in [3·2^14, 2^16), 1.11·10⁶ ordered pairs, N ≡ 1 mod 3).**

| curve | adv(a_N mod 3) | adv(r_N(ψ₃)) | mod-3 proxy surjective |
|---|---|---|---|
| model | 0.156 | 0.226 | — |
| 37a1 | 0.153 | 0.223 | yes |
| 92: y² = x³ − x + 1 | **0.187** | 0.242 | yes |

- **P4a: REFUTED as registered in 1 of 4 cells** (curve 92, a_N: |0.187 − 0.156| = 0.031 > 0.03). The
  finite-population Bayes advantage is upward-biased and noisy; ecount PH2's band at ℓ = 3 has width ≈ ±0.027. This
  explanation is post hoc.
- P4b holds: 0 prime and 0 pair exceptions.

**E5.** Both X_E members split at trial 1. P5 holds.

**Post-hoc (labelled).**
- PH1: identifiability classes (§2.6).
- PH2: information per kind (§2.8).
- PH3: prime-power moduli (§2.8).
- PH4: E3a real − sim summary.

---

## 4. Literature (cx1 checked primary texts where stated; the literature critic's verdicts go in §9)
- **Bach–Charles**, "The hardness of computing an eigenform", arXiv:0708.1192 v2, §2, Thm 2.2, Rem 2.5: oracle values
  a(N), a(N²). This is the BC2 mechanism. A single value is not covered.
- **Kunihiro–Koyama**, EUROCRYPT '98, LNCS 1403, 47–58. The CRT construction (counts mod small primes ⇒ count) is
  verified in Kunihiro's thesis (U. Tokyo 2001), Ch. 4, Lemma 4.2, Thm 4.2 (cx1). Conference numbering UNVERIFIED.
- **Okamoto–Uchiyama**, EUROCRYPT '98, LNCS 1403, 546–560, §4, Thms 9–10. **Galbraith–McKee**, IMA C&C 2005, LNCS 3796,
  §2. **Martín–Morillo–Villar**, Appl. Math. Lett. 14 (2001) 341–346 (point orders; full text UNVERIFIED).
  **Dąbrowski–Pomykała–Shparlinski**, J. Complexity 76 (2023) 101741, arXiv:1912.00345 v4: the deterministic theorem is
  **Thm 1.2** in v4 (cx1; ecount wrote 1.1). **Jiménez Urroz–Pomykała**, Ramanujan J. 64 (2024), Thm 1 (GRH, varying
  curves and twists).
- **Edixhoven–Couveignes** (Ann. Math. Studies 176, 2011), Thm 15.2.1 / Cor 15.2.2: the input includes the
  factorisation of n.
- **Sato–Tate:** Barnet-Lamb–Geraghty–Harris–Taylor, Publ. RIMS 47 (2011) 29–98 (non-CM over totally real fields,
  including Q).
- **Castryck–Hubrechts**, distribution of #E(F_p) mod m over random curves (Ramanujan J., c. 2013): UNVERIFIED, cited as
  motivation for the MC heuristic only.
- **Negative findings (cx1 and us):** no prior statement of a_N-alone or (a_p, a_q) ⇒ FACT, of a scalar-twist
  obstruction for trace products, of a multi-curve trace decoding, or of the j = 0 sextic everywhere reduction.

## 5. Candidate cards (rule 9)

| Mechanism | Verdict | Scope |
|---|---|---|
| Coppersmith/lattice on a_pa_q = A | refuted for this model | No equation contains p (DEC) |
| Galois congruences of one curve's traces → p mod ℓ → CRT | refuted (MODEL) | ≈ 0.3 bits total (GC) |
| Congruences of a_N at many curves → s mod ℓ | refuted (MODEL) | Theorem SC: QR bits only; E3b/E3d |
| Ordered/unordered traces at many curves → CRT | works HEURISTICALLY | Needs a source of such trace numbers; E3c 30/30 |
| Twist methods on X_E | blind (Prop QA) | Except descent and j = 0 sextic twists |
| Polynomial-identity annihilators from traces | impossible (DEC) | Smoothness or arithmetic annihilators not excluded |

## 6. Most decisive open questions
1. **a_N or (a_p, a_q) of one non-CM curve.** Does any arithmetic (non-algebraic, non-congruence) use of the integer
   a_N factor N almost everywhere? Factoring a_N (n/2 bits) is the obvious first step.
2. **A source of trace numbers.** A poly-time method producing a_N(E) at composite N for many E would only give QR bits
   (SC). Ordered traces or M mod ℓ would factor N (heuristically).
3. X_E for j ∉ {0, 1728} without 2-torsion. Determinism.
4. Fixed-ℓ counts: computational relation to QR_N (no reduction either way beyond §2.7).

## 7. Reproduction and deviations
```
.venv/bin/python -m pytest -q tests/test_tracedata.py           # 20 tests, <1 s
.venv/bin/python -m src.tracedata.run_experiments e1            # 17 s
.venv/bin/python -m src.tracedata.run_experiments e2            # <1 s
.venv/bin/python -m src.tracedata.run_experiments e3a           # 67 s
.venv/bin/python -m src.tracedata.run_experiments e3b           # 5 min
.venv/bin/python -m src.tracedata.run_experiments e3c           # 36 s
.venv/bin/python -m src.tracedata.run_experiments e3d           # 73 s
.venv/bin/python -m src.tracedata.run_experiments e4            # <1 s
.venv/bin/python -m src.tracedata.run_experiments e5            # <1 s
.venv/bin/python -m src.tracedata.run_experiments posthoc       # PH1-PH4, ~1 min
```
- **D1.** The end-to-end decoders test up to 4 CRT lifts (`max_lifts`). This is not in the prereg and is harmless:
  W > p already, and every candidate is verified exactly.
- **D2.** E3b/E3d were designed on the false premise that products identify s mod ℓ (they identify only the square
  class; Theorem SC, found from E3b's plateau). P3c/P3d/P3f pass as registered, but E3d is uninformative.
- **D3.** The first E4 run crashed on a JSON serialisation of numpy booleans. It was rerun unchanged otherwise.
- **D4.** E4's curve C5 = y² = x³ − x + 1 was fixed in the prereg text before any run (commit f2df008).

## 8. Codex rounds
- cx1 (blind) supplied: Prop PROMISE, Prop NA, the M/φ incomparability, Prop QA, Theorem J0, the (T, A, M) a = b = 1
  case, the Hasse-box count, the Bach-bound counterexample, and the literature corrections (DPS 1.2, Kunihiro thesis).
- Adversarial round: pending (cx2).

## 9. Reviews
Pending: codex adversarial (cx2), Opus proof critic, Opus literature critic.
