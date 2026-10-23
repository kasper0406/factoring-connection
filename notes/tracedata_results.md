# tracedata: trace data a_p, a_q, a_N, fixed-ℓ partial counts, X_E and determinism (phase 25)

Agent `tracedata`, 2026-10-23. Version **v2.1**: after the codex adversarial round (cx2: MAJOR REVISION) and the Opus
proof critic (MAJOR PROBLEMS, all text-fixable); every required fix applied (§9). Opus literature critic (OK WITH FIXES) applied in v2.1; proof-critic round 2 on
Theorem EXT and Prop CM-PAIR: §9.

- Pre-registration: `notes/tracedata_prereg.md`, commit `f2df008`, made before any experiment.
- Code: OPERATIONAL `src/tracedata/{reductions,decode,xe_descent,ext_reduction,galois}.py` (inputs: N and the stated
  data only; an AST test forbids harness/PARI/factorisation imports and calls; cx2 also inspected by hand);
  HARNESS `src/tracedata/{harness,run_experiments}.py`.
- Tests: `tests/test_tracedata.py` (24 pass, < 1 s).
- Raw data: `experiments/raw/tracedata_{e1,e2,e3a,e3b,e3c,e3d,e4,e5,posthoc,posthoc5,posthoc6}.jsonl`.
- GPT-6.1-sol: blind round `codex/tracedata_1.{md,out}` (**cx1**), adversarial round `codex/tracedata_2.{md,out}`
  (**cx2**). Scratch scripts in `codex/scratch/tracedata_{1,2}/`.
- Labels: KNOWN (with reference), PROVED, MODEL (exact computation or proof inside the Chebotarev/Haar model),
  HEURISTIC, experimentally supported, REFUTED (prediction), OPEN, UNVERIFIED. Strongest novelty wording: "no prior
  statement found by us".

Notation: N = pq, p < q, gcd(N, 6Δ_E) = 1; E/Q fixed; a_r = r + 1 − #E(F_r); A = a_N = a_pa_q; T = a_p + a_q;
s = p + q; M = #E(Z/N) = (p+1−a_p)(q+1−a_q); S = qa_p² + pa_q²; g_r = ρ_{E,ℓ}(Frob_r) ∈ GL₂(F_ℓ); t_r = tr g_r,
d_r = det g_r = r mod ℓ. "X ⇒ FACT" means a polynomial-time algorithm factoring N from N and X.

---

## 0. Verdict (plain language)

**Main positive result (cx2, re-implemented and checked here): ecount's exceptional set X_E is gone. For EVERY fixed
elliptic curve E/Q, the single count M = #E(Z/N) factors every good N = pq in randomized polynomial time.**
- The method computes [M]P for a twisted point over a random degree-5 algebra (Z/N)[z]/(F).
- With probability → 1/20, F has a linear factor mod p and is irreducible mod q. Then M kills the F_p component
  about half the time but can never annihilate a curve over F_{q⁵}, whose exponent exceeds q^{5/2} > M.
- PH6: 240/240 instances (four curves, n = 48–96) and both known X_E members split, in 1.2–1.9 trials on average.
  The null control (wrong M) gave 0 splits in 2400 trials.
- Galbraith–McKee 2005 §3 assert in one sentence, without analysis, that their §2 methods extend to (Z/N)[x]/(f) by
  twisting with random elements. The degree-5 exponent argument and the every-curve/every-N claim are not there.
  Dieulefait–Urroz (2020) state explicitly that factoring from one count was not known. No prior statement of the
  uniform result was found by cx2, the literature critic or us; priority UNVERIFIED.
- **A deterministic one-curve reduction remains OPEN.**

**One curve's trace pair (a_p, a_q), or its Hecke eigenvalue a_N alone, is still not known to factor N.** Four
precise obstacles:
1. **Squarefree Hecke data contain no determinant.** For non-CM E, no FIXED algebraic relation ties (N, a_p, a_q) to p
   on events of positive probability (Prop DEC). Bach–Charles's a_{N²} = a_N(Sym²f) works because it carries
   det Frob_p = p. Any single datum carrying p closes the algebra: a_{N³}, a CM weight-3 partner, or M (§2.2–2.3).
2. **One curve's congruences carry O(1) bits.** In the full-Haar model, all congruences of (a_p, a_q) carry ≈ 0.3 bits
   about p's residues (Prop GC). Real curves add entanglement terms: +≈0.7 bits from Q(√Δ) ⊂ Q(E[2]), and exact p mod 5
   for 11a1. CRT needs n/2 bits.
3. **Congruences of trace products see only quadratic residuosity, given scalars in the image (Theorem SC, MODEL; no
   prior statement found by us).** If the joint mod-m image of E_1, …, E_K contains the scalars, the Galois content of
   (N, a_N(E_1), …, a_N(E_K)) is invariant under (p, q) ↦ (λ²p, λ⁻²q).
   - Curves with rational torsion or isogenies at ℓ escape. For them a_N gives s mod ℓ, but only at finitely many ℓ.
   - E3b shows the plateau; E3d fails 0/30.
   - SC is not a factoring barrier. It concerns congruences only, and the integers themselves can factor N: next item.
4. **Promise obstruction (cx1).** For y² = x³ − x on Blum moduli, all trace data are 0 or powers of N, computable
   without factoring. So an everywhere trace-only reduction for this curve would factor Blum integers outright.

**Several fixed curves (proof critic, PROVED, re-checked here).**
- For p, q ≡ 1 (mod 4), a_N(y² = x³ − x) together with a_N(y² = x³ − Dx), for any D that is a non-residue at p or at
  q, factor N deterministically through two representations as sums of two squares.
- PH5: every one of 217 such N (n = 64, 128) factored with D ≤ 59. Under ERH, D ≤ 2 ln²N suffices.
- The critic reports that sextic twists of y² = x³ + D do the same for p, q ≡ 1 (mod 3) (not re-verified).
- So Hecke eigenvalues at N of a handful of fixed CM curves factor a positive density (7/16) of N, even though their
  congruences see only QR bits (SC).
- Ordered traces of two CM curves give p = (a_p(E_1)² + a_p(E_D)²)/4 directly, on density 15/16.
- For non-CM families, ordered trace numbers of O(n² log n) fixed curves factor N under a Chebotarev-independence
  heuristic (decode p mod ℓ, then CRT): E3a matches the model in 40/40 cells; E3c 30/30 toy N.
- So FROB's "numbers, not operators" gap closes for trace numbers of several curves: rigorously for CM pairs,
  heuristically for non-CM families. It stays open for one curve.

**Fixed-ℓ counts** (M mod ℓ, a_N mod ℓ, root counts of a fixed polynomial):
- Each has O(1) bits, so none is factoring-hard under N-local reductions unless factoring has a randomized
  polynomial-time algorithm (Prop IL).
- They sit above fixed-base QR: r_N(x² − Δ) ≡ QR_N(Δ); QR_N(Δ_f) ≤ r_N(f) for S₃ cubics; a CM curve's a_N decides
  QR_N(−D).
- They carry non-abelian information that no finite set of abelian symbols determines (Prop NA, cx1).
- Theorem PT: a pure algorithm computing any of them exactly is itself a splitting algorithm at rate ≥ 0.16–0.36.
  This does not exclude efficient pure algorithms.
- Whether any is computable in polynomial time: OPEN.

---

## 1. Equivalence table (oracle relations given N; D deterministic, R randomized)

| Data (fixed curve(s)) | Relation to FACT | Status / source |
|---|---|---|
| M = #E(Z/N), any fixed E | ≡_R FACT on **every** good N | Theorem EXT (cx2) §2.9; PH6 240/240; TW* (ecount) off X_E |
| (a_N, a_{N²}) | ⇒_D FACT unless a_p = a_q = 0 | KNOWN, Bach–Charles; E1a: 450/450 matched the prediction |
| (a_N, a_{N³}) | ⇒_D FACT iff a_N ≠ 0 | PROVED §2.2; E1b: 450/450 matched the prediction |
| (a_N(f), a_N(g)) with a_p(g) = c·a_p(f)² − e·p at both primes | ⇒_D FACT unless a_p = a_q = 0 | PROVED §2.2 (Sym²: c = e = 1; CM weight 3 at p, q ≡ 1 mod 4: c = 1, e = 2); E1d 90/90 |
| a_N(y² = x³ − x), a_N(y² = x³ − Dx), D ≤ B | ⇒_D FACT if p ≡ q ≡ 1 (mod 4) and some D ≤ B is a non-residue at p or q (ERH: B = 2 ln²N) | PROVED §2.5 (proof critic); PH5 217/217 |
| ({a_p, a_q}, M) or (T, A, M) | ⇒_D FACT except a_p = a_q = 1 (then M = N; ⇒_R via TW*) | PROVED §2.3 |
| (a_p, a_q), ordered or unordered (≡ up to ×2), one curve | **OPEN**; an everywhere reduction for a CM curve would factor Blum integers | §2.1, §2.4 |
| a_N, one curve | **OPEN**; ≤ (a_p, a_q) + factoring a_N (n/2 bits); ⇒ s mod |T| if E(Q)_tors = T; ⇒ QR_N(−D) if CM | §2.5 |
| (a_N(E_i))_{i ≤ K}, scalars in the joint image | congruence content ⊆ quadratic classes of p | Theorem SC (MODEL) §2.6; E3b, E3d |
| (a_p(E_i))_{i ≤ K} ordered, non-CM family, K = O(n² log n) | ⇒ FACT | HEURISTIC §2.8; E3a, E3c 30/30 |
| ({a_p(E_i), a_q(E_i)})_i | identifies s mod ℓ ⇒ FACT, K = O(ℓ² log ℓ) per ℓ | HEURISTIC; identifiability MODEL (PH1); untested end-to-end |
| (M_i mod ℓ)_i, many curves | identifies s mod ℓ; generic per-curve information Θ(ℓ⁻⁵) (cx2) | MODEL; sample bound K = Ω(ℓ⁵ log ℓ); sufficiency UNVERIFIED |
| M mod ℓ, a_N mod ℓ, r_N(f), fixed ℓ, f | O(1) bits: not FACT-hard under N-local reductions unless factoring is in randomized poly time | Prop IL §2.7 |
| r_N(x² − Δ) | ≡_D QR_N(Δ) | PROVED (cx1) |
| r_N(f), f an S₃ cubic | QR_N(Δ_f) ≤_D r_N(f) ≤_D FACT; not determined by finitely many abelian Frobenius data | ecount (cx2); Prop NA (cx1) |
| M mod ℓ (one curve) vs φ(N) mod ℓ, odd ℓ | neither determines the other at the Frobenius level | PROVED (cx1) |

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
- DEC excludes FIXED (n-independent degree and coefficients) polynomial equalities and fixed algebraic formulas for p
  in terms of (N, a_p, a_q). It does not exclude growing-degree equations, modular small-root relations with unknown
  modulus p (Coppersmith), auxiliary variables, arithmetic branching, gcds, group structure, Galois congruences, or
  factoring a_N (proof critic #4, cx2 #11). "a_pa_q = A contains no p" disposes of that one direct equation only.
- Likewise p is not algebraic over the formal field Q(N, M) on positive-probability events, yet M factors N through the
  group law (TW*, and §2.9 everywhere). So DEC is not evidence of hardness. It locates the missing ingredient: the
  determinant.
- **Polynomial annihilators.** No nonzero K ∈ Z[N, X, Y] has (P + 1 − X) | K(PQ, X, Y) in Z[P, Q, X, Y]. Substituting
  X = P + 1 gives K(PQ, P+1, Y) = 0, and (PQ, P+1, Y) are algebraically independent.
- So any multiple of #E(F_p) that the traces could supply cannot come from an identity. It must be smoothness-type
  (ECM constants) or arithmetic.

### 2.2 Bach–Charles and what a_{N²} supplies. KNOWN + PROVED variants

- Hecke at p²: a_{p²} = a_p² − p. The determinant det ρ(Frob_p) = p enters here and nowhere in squarefree a_N.
- a_{N²} = a_N² − S + N with S = qa_p² + pa_q², and the roots of Z² − SZ + N·A² are qa_p² and pa_q².
  - Their gcds with N give q and p when a_pa_q ≠ 0.
  - If exactly one trace vanishes, gcd(S, N) splits N.
  - If both vanish, the data are (0, N) and carry nothing.
  - KNOWN form: Bach–Charles, arXiv:0708.1192 v2, §2, Thm 2.2 and Rem 2.5; weight-2 form in ecount §3.
- Interpretation: a_{N²} = a_N(Sym²f), the Hecke eigenvalue at the squarefree N of the GL₃ form Sym²f, since
  tr Sym²ρ(Frob_p) = a_p² − p. Bach–Charles use two "trace numbers at N", one of which contains det.
- **Prop BC3 (PROVED).** a_{p³} = a_p³ − 2pa_p, so a_{N³} = A(A² − 2S + 4N). For A ≠ 0 this gives S, hence the
  factors. A = 0 kills all odd powers.
- **Prop LINK (PROVED).** If a_r(g) = c·a_r(f)² − e·r at r = p and r = q (c, e ≠ 0 fixed), then
  a_N(g) = c²A² − ceS + e²N, which gives S. Exception: A = 0 = S, e.g. Sym² of 32a2 at Blum N (cx2).
  - Instance: f is the weight-2 CM form of 32a2 (y² = x³ − x), and g = η(4z)⁶ ∈ S₃(Γ₀(16), χ₋₄).
  - At split r ≡ 1 (mod 4), a_r(g) = π² + π̄² = a_r(f)² − 2r. At r ≡ 3 (mod 4), a_r(g) = 0 ≠ −2r (proof critic #7).
  - So the eigenvalues at N of two FIXED CM forms factor N deterministically on split–split N.
  - E1d: 90/90. Checked against the η-product q-expansion in the tests.
  - No prior statement found by us. It is an immediate variant of Bach–Charles.
- **Single a(N^k): OPEN.** cx1: on Blum N for 32a2, a_{N^{2m}} = N^m and a_{N^{2m+1}} = 0, so no single higher
  coefficient suffices everywhere for CM curves.

### 2.3 Prop MIX (traces + one coupling datum). PROVED

- Given {a_p, a_q} (equivalently (T, A): roots of Z² − TZ + A) and an observed datum φ_obs = Φ(p, q, a_p, a_q), with Φ
  a fixed rational function.
- For each assignment of the traces, clear denominators in Φ(P, N/P, a, b) − φ_obs. If the resulting polynomial in P is
  not identically zero for the ACTUAL (a, b, φ_obs), its integer roots include p.
- This is deterministic and needs no balance. The condition depends on the values, not only on Φ: Φ = T is useless
  (proof critic #8, cx2 #5).
- Instance M: with u = 1 − a_p, v = 1 − a_q,

  v·p² + (N + uv − M)·p + uN = 0,

  equivalently C = M − N − 1 + T − A = (1 − a_q)p + (1 − a_p)q (cx1).
  - The polynomial vanishes identically exactly when a_p = a_q = 1, which forces M = N.
  - cx2's example: y² = x³ − 7x − 2, p = 5, q = 7, with `mix_factor(35, [1, 1], 35)` = None (unit test).
  - In that case use TW*. cx1's argument: if both twist EXPONENTS divided pq = M, then since p ∤ p + 2 the twist group
    at p would have order q or q². q² is too large, so q = p + 2; symmetrically q + 2 = p², which is impossible. So
    N ∉ X_E. (Theorem EXT now covers it anyway.)
  - So: ⇒_D except at a_p = a_q = 1, and ⇒_R everywhere.
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

### 2.5 a_N alone, and a_N of several CM curves

- **Hasse box (cx1).** Given A ≠ 0, the integer trace pairs explaining it number ≤ 2τ(|A|). This bounds candidates,
  but identifying primes having those traces is the whole problem. a_N ⇒ (a_p, a_q) up to factoring A (n/2 bits;
  subexponential, not polynomial) and τ(A) candidates.
- **Lemma TOR (PROVED).** If E(Q)_tors has order m, then m | #E(F_r) at good r ∤ m. So a_r ≡ r + 1 and
  a_N ≡ (p+1)(q+1) = N + 1 + s (mod m), which gives s mod m (the φ(N) mod m information).
  - Ordered (a_p, a_q) give p mod m.
  - Mazur bounds the torsion subgroup. lcm over curves = 5040, using |T| = 16 (Z/2×Z/8), 9, 5, 7. This is
    fixed-modulus information only.
  - E1e: 270/270.
- **Lemma CM-QR (PROVED).** Let E have CM by O_K, K = Q(√−D) of class number one, and suppose (−D/N) = 1.
  - At split r > D, a_r = π + π̄ ≠ 0, since π = −π̄ would force r = D·(square)/4. So a_N ≠ 0 iff both primes split.
  - So a_N decides QR_N(−D) exactly. For 32a2 this decides ε (whether N is a sum of two squares). E1f: 45/45.
- **ℓ = 2.** a_N is odd iff M is odd iff r_p = r_q = 0 (ecount's 2/9 one-sided certificate).
- **Prop CM-PAIR (proof critic #1; PROVED; re-checked here).** Let r ≡ 1 (mod 4) and write r = α_r² + β_r², α_r odd.
  - For E_1: y² = x³ − x, a_r(E_1) = ±2α_r.
  - For E_D: y² = x³ − Dx, a_r(E_D) = ±2α_r if D is a square mod r, and ±2β_r if not (quartic residue symbol of D;
    Gauss/Jacobsthal).
  - Given A_1 = a_N(E_1) and A_D = a_N(E_D), put X = A_1/4 and Y = A_D/4.
  - If D is a non-residue at both primes, then X = ±α_pα_q and Y = ±β_pβ_q. The numbers u = α_pβ_q and v = β_pα_q
    satisfy u² + v² = N − X² − Y² and u²v² = X²Y². So u², v² are the integer roots of a known quadratic, and
    gcd(N, X² + v²) = gcd(N, α_q²p) = p.
  - If D is a non-residue at exactly one prime, then Y = ±α_pβ_q (say), and gcd(N, X² + Y²) = gcd(N, α_p²q) = q.
  - If D is a residue at both primes, there is no information.
  - `cm_quartic_pair_factor`, `cm_quartic_family_factor`.
  - PH5 (post hoc): with D ∈ {2, …, 59}, all 217 N with p ≡ q ≡ 1 (mod 4) among 800 random N factored (n = 64 and 128). With D = 2 alone,
    exactly the classes (p, q) mod 8 ∈ {(1, 5), (5, 5)} factored, as the proof predicts. No other class factored: there
    a_N = 0.
  - Under ERH the least D that is a non-residue at p is ≤ 2 ln²p (Bach 1990). So the first ⌈2 ln²N⌉ curves of the fixed
    sequence y² = x³ − Dx suffice for every N with p ≡ q ≡ 1 (mod 4).
  - The proof critic reports the sextic analogue (y² = x³ + D, p ≡ q ≡ 1 mod 3; 2966/3000 recoveries). We have not
    re-verified it. Together these give density 7/16.
  - Ordered traces give p = (a_p(E_1)² + a_p(E_D)²)/4 directly whenever p splits in Q(i) (or Q(ζ₃)): density 15/16.
  - Consistent with SC: CM images contain the scalars, so the congruence content is QR-only. The factoring power comes
    from the archimedean identity α² + β² = r.
  - No prior statement found by us or the literature critic. Bach–Charles Rem 2.6 only notes that CM forms fail when
    a(p) = a(q) = 0.
  - Closest antecedent: DPS 2023 use the same family y² = x(x² − b) with b up to a least-non-residue bound. Their
    mechanism is different (2-adic orders from count multiples, p, q ≡ 3 mod 8), so the results are complementary.
  - DU's Thm 5 (one curve plus quadratic twists d = O(log²N), ERH) has the same shape. Quadratic twists change a_N only
    by (d/N), which is why quartic twists are needed.
  - Ingredients are classical, with α odd and β even:
    - Jacobsthal, J. Reine Angew. Math. 132 (1907) 238–245;
    - Berndt–Evans–Williams, *Gauss and Jacobi Sums* (1998), Ch. 6;
    - Ireland–Rosen GTM 84, Ch. 18 §4;
    - two representations as sums of two squares factor N.

### 2.6 Theorem SC (scalar-twist symmetry of trace products). MODEL (Chebotarev); no prior statement found by us

**Setting.**
- E_1, …, E_K/Q and a modulus m. G ⊆ ∏_i GL₂(Z/m) is the image of Galois under ⊕_i ρ_{E_i, m}.
- Assume G ⊇ {(λI, …, λI) : λ ∈ Λ} for a subgroup Λ ⊆ (Z/m)^×. This holds with Λ = (Z/m)^× for the full fibre
  product {(g_i) : det g_i all equal}, at every prime ℓ where all images are surjective, and for CM images (Cartan
  normalisers).
  - It FAILS at the finitely many ℓ where some E_i has a rational ℓ-isogeny or ℓ-torsion. For 11a1 mod 5, Λ = {1}.
    There a_N ≡ χ(N)(1+p)(1+q) gives s mod ℓ (Lemma TOR). The headline "QR bits only" therefore needs this hypothesis
    (proof critic #3, cx2 #2).
  - (Serre-curve entanglement makes the full fibre product false at composite m divisible by the entanglement level.
    Scalars still lie in G there.)
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
- The finite-image transformation τ_λ has no automatic lift to the integer prime/trace instances (cx2's wording).
- So any method extracting more than QR bits from Hecke eigenvalues at N (under the hypothesis) must use their
  archimedean/integer structure, and some do:
  - Prop CM-PAIR (sums of two squares);
  - Bach–Charles's pair (a_N, a_N(Sym²)). SC holds for it too, since tr Sym^j(λg) = λ^j tr Sym^j g (cx2), yet the
    integers factor N.
  So SC is a statement about fixed-modulus congruence decoding, NOT a factoring barrier. Fixed-(m, K) Chebotarev cannot
  be used uniformly when m and K grow with N.
- **Contrast.**
  - Ordered traces t_p (or (t_p, t_q)) are not τ-invariant. They identify p mod ℓ completely.
  - Unordered pairs {t_p, t_q} and M mod ℓ = det(1 − g_p)det(1 − g_q) identify exactly {p, q} mod ℓ, i.e. s mod ℓ
    (PH1: exact identifiability classes for ℓ = 5, 7, 11, 13, all N mod ℓ; unit-tested).
- **Attainability (MODEL; cx2, needs independent full-image copies, not just scalars).** With many curves the products
  do identify (p/ℓ) on (N/ℓ) = 1, i.e. QR_N(ℓ*). Already the zero frequencies differ: P(t = 0 | d) = (ℓ + χ(−d))/(ℓ² − 1).
  SC is the upper bound; this is a separate attainability statement. No computational reduction to QR oracles is
  claimed.
- **Experiment.** E3b: the product decoder's accuracy for s mod ℓ is flat in K (ℓ = 5: 0.64, 0.64, 0.64 at
  K = 10³, 10⁴, 4·10⁴). It matches the model in 15/15 cells.
  - This is the plateau of OUR decoder (`rank_s` maximises over determinant preimages). It is not exactly the Bayes
    ceiling. The limiting Bayes ceilings are 1, 0.625, 0.500, 0.300, 0.250 for ℓ = 3–13; our tie convention gives
    1, 0.625, 0.444, 0.280, 0.236 (cx2 #8).
  E3d (end-to-end from products) 0/30, model 0.03. The pre-registration wrongly expected s mod ℓ to be decodable
  from products; the plateau exposed the symmetry (deviation D2).

### 2.7 Fixed-ℓ partial counts

**Prop IL (instance-local bits). PROVED (trivial).**
- Definitions:
  - Let g be a function on moduli whose value on N lies in a set of size c independent of N, e.g. M mod ℓ, a_N mod ℓ,
    r_N(f) ∈ {0, …, (deg f)²}.
  - Call a reduction N-local if each oracle answer it requests is computable in polynomial time from g(N) and the
    query, e.g. queries N^j·k with k's local data computable, using multiplicativity:
    #E(Z/N^jk) = N^{j−1}M·#E(Z/k), and r_{N^j}(f) = r_N(f) for p, q ∤ disc(f)·lc(f).
- Claim: an N-local reduction that factors N with probability σ yields, by running it with each of the c possible
  answers and verifying the outputs, an oracle-free algorithm with success ≥ σ and c-fold cost.
- Quantifiers (cx2 #12): a per-input (worst-case) reduction gives worst-case randomized polynomial-time factoring, and a
  deterministic reduction gives deterministic factoring. Distributional success gives only distributional factoring.
- So no fixed-ℓ count is factoring-hard under N-local reductions unless factoring is easy in the same sense.
  Fixed-base QR_N(Δ) has the same status.
- Not covered: reductions that query g on moduli coprime to N with unknown factorisations.

**Placement.**
- QR_N(Δ) ≡_D r_N(x² − Δ) (cx1).
- QR_N(Δ_f) ≤_D r_N(f) ≤_D FACT for S₃ cubics (ecount/cx2). The root count of a fixed polynomial is computable after
  factoring via deg gcd(f, x^r − x).
- **Prop NA (cx1, PROVED; wording per proof critic #10, cx2 #12).** r_N(f) for an S₃ cubic, and M mod ℓ for full image
  (odd ℓ), are not determined by N mod m together with finitely many abelian Frobenius data (quadratic or other
  abelian symbols) of p and q, for any fixed m. (N itself determines everything mathematically.)
  - The S₃ (resp. GL₂) field meets any multiquadratic or cyclotomic field inside its maximal abelian subfield. So
    Chebotarev supplies primes with equal abelian data but Frobenius identity vs 3-cycle (resp. I vs −I:
    det(1 − I) = 0, det(1 + I) = 4).
  - This is non-determination, not a computational separation from QR oracles.
- **M mod ℓ vs φ(N) mod ℓ (cx1, PROVED, odd ℓ).** Neither determines the other at the Frobenius level. At ℓ = 2,
  φ(N) mod 2 is constant and I = −I; use the S₃ identity/3-cycle witness instead.
  - Fixed determinants allow varying M mod ℓ.
  - diag(1, d) has det(1 − g) = 0 for all d, so determinant pairs (1, 1) and (−1, −1) give the same (N, M) mod ℓ and
    different φ(N) mod ℓ.
- **Many curves (MODEL, PH1/PH2, cx2 #7).** (M_i mod ℓ)_i over many independent full-image curves identifies exactly
  s mod ℓ, i.e. the φ(N) mod ℓ information.
  - The AVERAGE per-curve information is ≈ 0.736/ℓ⁴ bits at ℓ = 31 and still decreasing (0.726 at ℓ = 101).
  - The average is dominated by exceptional determinant strata. On generic strata (c ≠ 1; d, e ≠ 1; d ≠ e), cx2 proves
    the likelihoods differ by O(ℓ^{−7/2}) per nonzero output, so the information is O(ℓ⁻⁵) and Fano gives
    K = Ω(ℓ⁵ log ℓ). A matching sufficient bound is UNVERIFIED. v1's "K ~ ℓ⁴, O(n⁴ log n)" is WITHDRAWN.
  - This is consistent with Kunihiro–Koyama's CRT reduction for arbitrary input curves.
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

  - These constants are exact under the full-image / S₃ laws (cx2 and the proof critic recomputed them, e.g. M mod 3:
    81/512). They are not universal over torsion curves: rational ℓ-torsion makes M mod ℓ constant.
  - With the Chebotarev class law in place of μ, error o(1).
  - Meaning (corrected, cx2 #10): an efficient pure algorithm computing one of these counts IS an efficient splitting
    algorithm on that population at the stated rate. PT does not exclude efficient pure algorithms.
  - PT-J transfers to a fixed non-square Δ, since (Δc|x) = (Δ|x)(c|x): m_±(S) = β(c_S) ± β(Δc_S), and at most one of
    c_S, Δc_S is a square. Landau–Page goes through with the fixed conductor factor. The restrictions of PT-J1 (heights,
    number of tests, equal priors, units, non-adaptive observations) carry over unchanged.
- **Status: OPEN** whether any of them is computable in polynomial time. A poly-time r_N(f) for an S₃ cubic would
  decide QR_N(Δ_f).

### 2.8 Multi-curve Galois decoding. MODEL constants + HEURISTIC transfer + toy agreement

**Model.** The Haar likelihood is P(t | d) = (ℓ + χ(t² − 4d))/(ℓ² − 1) (KNOWN: Gekeler 2003; Castryck–Hubrechts Thm 1(ii)). Exact information per curve in bits
(PH2, E2; v1's unordered law double-counted diagonal pairs, fixed per cx2 #4, numbers regenerated; identifiability
partitions are unaffected):

| observation | scaled information at ℓ = 11 / 31 (finite-ℓ values, not proved limits) | what K → ∞ identifies |
|---|---|---|
| t_p (one prime, ordered) | ℓ²I = 0.658 / 0.698 | p mod ℓ |
| (t_p, t_q) ordered | ℓ²I = 1.314 / 1.397 | p mod ℓ |
| {t_p, t_q} unordered | ℓ²I = 0.660 / 0.699 | s mod ℓ |
| t_pt_q (a_N) | ℓ³I = 0.769 / 0.735 (0.725 at ℓ = 101) | (p/ℓ) on (N/ℓ) = 1 only (SC) |
| det(1−g_p)det(1−g_q) (M) | ℓ⁴I = 0.753 / 0.736 (average; generic strata Θ(ℓ⁻⁵)) | s mod ℓ |

**Prop GC (one curve; full-Haar MODEL, cx2 #3, proof critic #2).**
- I(p mod ℓ; (t_p, t_q) | N mod ℓ) summed over odd ℓ ≤ 47 is 0.202 bits. Adding 3 → 9 (+0.013), 5 → 25 (+0.002) and
  the mod-8 part (0.072) gives 0.288. Further tails are small but not bounded here (8 → 16: +0.003; 9 → 27: +0.001).
- "≈ 0.3 bits" summarises the independent full-Haar product model. It is not a proved bound.
- Real curves are entangled. Q(√Δ) ⊂ Q(E[2]) ∩ Q^ab, and a_p odd ⇒ (Δ/p) = +1 (0 violations in 20000 primes for
  37a1 and 92). This adds I((Δ/p); a_p, a_q mod 2 | (Δ/N)) ≈ 0.70 bits for curves with S₃ mod-2 image (proof critic).
  Rational ℓ-torsion gives p mod ℓ exactly (11a1: two bits).
- With finite-index adelic image, the total is plausibly O(1) bits for every fixed curve (HEURISTIC: not proved). CRT
  needs n/2 bits.
- So, in the restricted statistical model where (a_p, a_q) are used only through their residues as Frobenius data, one
  curve cannot reach p. This is not a theorem about arbitrary congruence uses of the integer traces: H(p | N) = 0.

**HEURISTIC MC (non-CM families).**
- Assumption: for a fixed large prime p, the Frobenius classes of p in the ℓ-division fields of a fixed family of
  curves behave like independent Haar draws conditioned on det = p mod ℓ. This is a short-family version of
  Castryck–Hubrechts equidistribution (Ramanujan J. 30 (2013) 223–242, final arXiv v3 Thm 2: random curves over F_q,
  error C·m²·ln ln m/√q for m coprime to q). That theorem does NOT cover the box 1 ≤ A, B ≤ 200.
- Sample complexity (cx2 #7, proof critic #5).
  - Distinct single-trace laws have squared Hellinger separation Ω(ℓ⁻²) (quartic character sums, Hasse). So
    K = C·ℓ²·log(ℓ/ε) suffices for top-1 decoding with error ε.
  - Numerically, reliable decoding needs K ≈ 2.8ℓ²log₂(ℓ/δ). At 1.4ℓ²log₂ℓ, the mutual-information count, ML accuracy
    is only ≈ 0.68.
- With ℓ ≤ L ≈ ½(ln 2)n ≈ 0.347n (∏ℓ > p), K = O(n² log n) ordered traces factor N. Alternatively L ≈ 0.173n plus
  Coppersmith's p mod W ≥ N^{1/4} bound with a margin.
- The polynomial-time claim needs reliable top-1 residues. Our toy list decoder searches 2^{π(L)} combinations, which
  is exp(Θ(n/log n)) for L = Θ(n) (cx2 #7).
- Unordered pairs (s mod ℓ) are identifiable. A separation-based sample bound is not given here.
- Effective Chebotarev cannot justify this: the field degrees exceed p. It stays a heuristic.
- Evidence (toy agreement with this decoder and model; cx2 #6):
  - E3a: 200 primes of 30 bits, curves y² = x³ + Ax + B. Real accuracy matches the Haar simulation within 0.10 in
    40/40 cells. The cells reuse the same 200 primes across nested K, so they are correlated checks.
  - Mean real − sim = −0.020 over the 27 unsaturated cells, with 18/27 negative (PH4).
  - Possible causes, untested:
    - the family is not independent: 38588 distinct j among 40000 curves, quadratic-twist chains like (1,1), (4,8),
      (9,27), …, and 401 members with rational 2-torsion (cx2);
    - non-full images.
  - The null control is at chance in 40/40 cells.
  - E3c: 30/30 toy N (28-bit primes) factored from 40000 ordered traces of one prime, with the truth top-ranked at every
    ℓ ≤ 29 in all 270 (N, ℓ) cells. The binomial 95% lower bound on the success rate is ≈ 0.905; this is one toy size
    and one family.
- These are toy correctness checks of the model, not evidence of a factoring algorithm: the harness supplies the traces.
- For CM pairs, the stronger rigorous Prop CM-PAIR applies instead.

**Meaning for FROB's scope gap.**
- Trace NUMBERS of several curves can carry factoring power: rigorously for CM pairs (CM-PAIR, via a_N or ordered
  traces), heuristically for ordered/unordered traces of non-CM families.
- Congruence decoding of products a_N cannot (SC, under its hypothesis).
- For one curve's pair, congruences carry O(1) bits (GC, model). Integrality remains OPEN.

### 2.9 X_E and determinism

**Theorem EXT (cx2 #1; PROVED; re-implemented independently in `src/tracedata/ext_reduction.py` and checked).**
For every fixed E/Q: y² = x³ + ax + b and every good N = pq with p, q ≥ 101, M = #E(Z/N) gives an expected-polynomial
factoring algorithm. No factorisation of M, no second count, and no extra automorphisms are needed. Primes < 101 are
found by trial division.

*One trial.*
- Pick a uniform monic quintic F over Z/N and set R = (Z/N)[z]/(F).
- Pick a uniform x ∈ R and set δ = x³ + ax + b. Take P = (δx, δ²) on Y² = X³ + aδ²X + bδ³ over R.
- Compute [M]P with affine arithmetic over R. Every non-unit met is examined: a proper gcd of its norm (5×5 determinant)
  or of its coefficients with N is output. Restart on an obstruction that reveals nothing.

*Proof.*
- **The event.** F ≡ (linear)(irreducible quartic) mod p and F irreducible mod q. By CRT its probability is
  (1 − p⁻²)/4 · (1 − q⁻⁴)/5 → 1/20, and on it R ≅ F_p × F_{p⁴} × F_{q⁵}.
- **Factor extraction on the event.** Every nonzero zero divisor w reveals a factor.
  - If w vanishes in F_{q⁵}, all its coefficients are divisible by q, so their gcd with N is q.
  - Otherwise w vanishes in some p-component, so gcd(norm w, N) = p.
  - Hence affine arithmetic either computes the componentwise sum or splits N. A terminal point that is infinity in
    one component and affine in another forces a split.
- **The F_p component is killed with constant probability.** In F_p, δ is a nonzero square with probability
  α_p = (#E(F_p) − 1 − r_p)/(2p) ≥ (p − 2√p − 3)/(2p). Then the twist is ≅ E, and M kills that component.
- **The F_{q⁵} component is never annihilated.** Put Q = q⁵. Any elliptic group over F_Q, including the quadratic twist,
  has exponent ≥ √#E(F_Q) ≥ √Q − 1 = q^{5/2} − 1. But M ≤ (√p + 1)²(√q + 1)² < q^{5/2} − 1 for q ≥ 101.
  - So M is not a multiple of either exponent. The kernel has index ≥ 2, and the non-killed probability is
    ≥ 1/2 − 5/(2Q).
- The two component conditions are independent, so a trial succeeds with probability ≥ 1/80 − o(1). Each trial costs
  O(log M) operations in a 5-dimensional algebra. ∎
- Not needed: recognising the splitting type, and balance (only M < q^{5/2} for the larger prime q, which holds always,
  because p < q).

*Checks.*
- PH6 (post hoc): 240/240 balanced N for 11a1, 37a1, y² = x³ + 2x + 3 and y² = x³ − x + 1 at n = 48, 64, 96, in
  1.2–1.9 trials on average. Other splitting types also split, so success is far above 1/80.
- Both known X_E members split (11a1 census member in 2 trials; cx1's N = 27331 in 1).
- Null control: M′ uniform in a ±2^{3n/4} window gave 0 splits in 2400 trials.
- cx2's own implementation split the same 11a1 member and a 60-bit N.

*Literature.*
- Galbraith–McKee (LNCS 3796, 2005, §3) is one paragraph. It defines R_{N,f} = (Z/N)[x]/(f), notes that f may split
  differently mod different primes, and says the §2 results "readily extend to this setting, twisting by random elements
  of R". There is no proof, no choice of degree, and no claim about exceptional N or curves (literature critic). That is
  the antecedent.
- Dieulefait–Urroz (arXiv:1911.11004 v1; Mathematics 8 (2020) 2126) wrote that "it is not known how to factor N only
  with the number E_N as input". JUP 2024 and DPS 2023 use varying curves.
- Bach–Shallit (Math. Comp. 52 (1989)) is a methodological antecedent for extension rings of Z/N.
- The uniform "every curve, every N" statement: no prior statement found by cx2 or us; priority UNVERIFIED.

*Consequences.*
- ecount's open question 1 (X_E for arbitrary fixed curves) is CLOSED.
- "#E(Z/N) for one fixed curve" is randomized-polynomial-time equivalent to factoring on all good N.
- The selfred §6(a) moved-divisor prototype is not a route below factoring anywhere.

**Earlier partial results (kept; superseded for the randomized question).**
- **Lemma SEP.** If λ_r and μ_r both divide an odd M, then λ_r ≠ μ_r, so some odd ℓ | M separates curve and twist at r.
  It gives no small-prime search bound.
- **Prop QA (cx1).** On X_E, for every finite étale QUADRATIC Z/N-algebra R, M annihilates E(R): from [2M]P = O and the
  fixed-point-free action of Frob² on nonzero 2-torsion, doubling is invertible. Twist methods over quadratic algebras
  are therefore blind on X_E. Degree 5 escapes because an F_{q⁵} group cannot be annihilated by M.
- **Theorem J0 (cx1).** For j = 0 curves, sextic twists c = y² − x³ give mixed probability ≥ 1/18 − o(1). cx1 does
  trial division below 37 and gets a uniform bound of 1/72. Exact check at N = 27331: 0.5035.
- **E5.** TW* with descent over ℓ ≤ 1000 split both X_E members at trial 1. That is not evidence of a universal
  small-prime descent.

**Determinism: OPEN.**
- A deterministic version needs, for instance, a small x₀ with (f(x₀)/p) ≠ (f(x₀)/q) plus a non-kernel condition.
- x ↦ (f(x)/r) is not multiplicative (cx1: f = x³ − x − 1 mod 13: f(2) = 5 and f(4) = 7 are both non-squares, so the
  symbol at 4 = 2² is not the square of the symbol at 2). So Bach's ERH least-non-residue bound does not apply.
- We know no GRH bound for least witnesses of such non-multiplicative character sums.
- EXT's random quintic F would also need derandomising: a deterministic choice of F with a prescribed splitting type
  mod both unknown primes.
- Rational points on fixed twists E^{(d)}(Q) remove the sampling of x. But they still need ord(G mod q) ∤ M for a fixed
  G, an elliptic Artin-type condition with no uniform guarantee.

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

**E3d.** 0/30 against sim 6/200 = 0.03, so P3f holds formally.
- It is uninformative: the registered design rested on a false premise (D2).
- The "sim" checks independent per-ℓ top-3 inclusion. It is not a joint integer-trace model.

**E4 (census, 1491 primes in [3·2^14, 2^16), 1.11·10⁶ ordered pairs, N ≡ 1 mod 3).**

| curve | adv(a_N mod 3) | adv(r_N(ψ₃)) | mod-3 proxy surjective |
|---|---|---|---|
| model | 0.156 | 0.226 | — |
| 37a1 | 0.153 | 0.223 | yes |
| 92: y² = x³ − x + 1 | **0.187** | 0.242 | yes |

- **P4a: REFUTED as registered in 1 of 4 cells** (curve 92, a_N: 69155/370018 = 0.1869 vs 0.15625, off by 0.0306).
- Post-hoc finite-prime variation, not a demonstrated cause:
  - The independent units are 1491 primes (strata 747/744), not 1.11·10⁶ pairs.
  - cx2 simulated 100000 Haar populations of that size: mean 0.1559, SD 0.0150, upper-tail frequency 0.019 at the
    observed value.
  - The proof critic found 0.1555 for curve 92 (and 0.1587 for 37a1) on 30000 primes near 2^22.
  - v1's "upward bias" and the borrowed ecount band are withdrawn (cx2 #9).
- P4b holds: 0 prime and 0 pair exceptions.

**E5.** Both X_E members split at trial 1. P5 holds.

**Post-hoc (labelled).**
- PH1: identifiability classes (§2.6).
- PH2: information per kind (§2.8).
- PH3: prime-power moduli (§2.8).
- PH4: E3a real − sim summary.
- PH5: CM quartic pairs (§2.5).
- PH6: Theorem EXT and its null control (§2.9).
- These were added after the reviews, are labelled post hoc, and are reproducible via `posthoc5` and `posthoc6`.

---

## 4. Literature (cx1/cx2 checked primary texts where stated; the literature critic's verdict is in §9)
- **Bach–Charles**, "The hardness of computing an eigenform", arXiv:0708.1192 v2; Contemp. Math. 463 (AMS 2008) 9–15
  (published numbering UNVERIFIED).
  - §2, Thm 2.2 (non-CM, density-1 set of RSA moduli), Rem 2.5 (balanced), Rem 2.6 (CM forms fail when
    a(p) = a(q) = 0: the antecedent of Prop PROMISE), Cor 2.8 (τ).
  - Oracle values a(N), a(N²). Single values, several forms and Sym² are not discussed.
  - Restated in Edixhoven–Couveignes Rem 2.2.4.
- **Kunihiro–Koyama**, EUROCRYPT '98, LNCS 1403, 47–58. The CRT construction (counts mod small primes ⇒ count) is
  verified in Kunihiro's thesis (U. Tokyo 2001), Ch. 4, Lemma 4.2, Thm 4.2 (cx1). Conference numbering UNVERIFIED.
- **Okamoto–Uchiyama**, EUROCRYPT '98, LNCS 1403, 546–560 (theorem numbers UNVERIFIED, paywall). Per secondary sources
  their attack factors N from a curve of order N itself: the a_p = a_q = 1, M = N case of §2.3 (also Kunihiro's thesis
  §5.5).
- **Galbraith–McKee**, IMA C&C 2005, LNCS 3796, 392–409: §2 (one curve with known M; twist and 2-adic methods), §3
  (extension algebras with differing splitting types: the antecedent of Theorem EXT).
- **Martín–Morillo–Villar**, Appl. Math. Lett. 14 (2001) 341–346 (point orders; full text UNVERIFIED).
- **Dąbrowski–Pomykała–Shparlinski**, J. Complexity 76 (2023) 101741, arXiv:1912.00345 v4: the deterministic theorem is
  **Thm 1.2** in v4 (cx1; ecount wrote 1.1).
- **Jiménez Urroz–Pomykała**, Ramanujan J. 64 (2024), Thm 1 (GRH, varying curves and twists).
- **Edixhoven–Couveignes** (Ann. Math. Studies 176, 2011), Thm 15.2.1 / Cor 15.2.2: the input includes the
  factorisation of n.
- **Sato–Tate:** Barnet-Lamb–Geraghty–Harris–Taylor, "A family of Calabi–Yau varieties and potential automorphy II",
  Publ. RIMS 47 (2011) 29–98, Cor 8.9 (non-CM E over totally real F); equivalently Thm B together with modularity
  (Breuil–Conrad–Diamond–Taylor 2001).
- **Castryck–Hubrechts**, "The distribution of the number of points modulo an integer on elliptic curves over finite
  fields", Ramanujan J. 30 (2013) 223–242, doi 10.1007/s11139-012-9444-0; final arXiv v3 (0902.4332), "completely
  rewritten because of an error".
  - Thm 2 error: C·p^e·N′²·ln ln N′/√q, over F_q-isomorphism classes.
  - Their Thm 1(ii) with n = 1, and Gekeler (IMRN 2003, no. 37, 1999–2018), give the §2.8 likelihood
    (ℓ + χ(t² − 4d))/(ℓ² − 1): KNOWN.
- **Bach**, Math. Comp. 55 (1990) 355–380: under ERH the least element outside a proper subgroup of (Z/m)^× is
  < 2 ln²m, so the least QNR mod p is < 2 ln²p (QNR corollary label UNVERIFIED). Sharper: Lamzouri–Li–Soundararajan,
  Math. Comp. 84 (2015) 2391–2412.
- **Kunihiro thesis** (U. Tokyo 2001) Thm 4.1 (random curves; heuristic uniformity), Lemma 4.2 (#E from #E mod d,
  CRT over ℓ ≤ log²(n+1)), Thm 4.2: confirmed by the literature critic from scans. This is consistent with Prop IL.
- **JUP**: Ramanujan J. 64 (2024) 265–273, doi 10.1007/s11139-023-00822-1.
- **GM Thm 2 and Thm 4** tie count comparisons and pairings to QR; relevant to §2.7.
- **Jacobsthal sums / Gauss:** a_p of y² = x³ − Dx in terms of p = α² + β² (classical).
- **Negative findings (cx1, cx2 and us):** no prior statement found of a_N-alone or (a_p, a_q) ⇒ FACT, of a scalar-twist
  obstruction for trace products, of multi-curve trace decoding, of CM-PAIR, of the j = 0 sextic reduction, or of
  Theorem EXT's uniform form.

## 5. Candidate cards (rule 9)

| Mechanism | Verdict | Scope |
|---|---|---|
| Direct small-root equation from a_pa_q = A | no such equation (it contains no p) | DEC excludes fixed algebraic relations only; growing-degree or modular Coppersmith constructions are NOT excluded |
| Galois congruences of one curve's traces → p mod ℓ → CRT | fails in the model | O(1) bits (GC; ≈0.3 full-Haar, ≈1 with entanglement) |
| Congruences of a_N at many curves → s mod ℓ | fails in the model | SC: QR bits only, given scalars in the image; E3b/E3d |
| a_N of CM quartic/sextic twist pairs | **works (PROVED)** on p, q split in the CM field | CM-PAIR; PH5 217/217 |
| Ordered/unordered traces at many non-CM curves → CRT | works HEURISTICALLY | Needs a source of such trace numbers; E3c 30/30 |
| Twist methods over quadratic algebras on X_E | blind (Prop QA) | — |
| Twist over a random quintic algebra | **works (PROVED)** for every curve | Theorem EXT; PH6 240/240 |
| Polynomial-identity annihilators from traces | impossible (DEC remark) | Smoothness or arithmetic annihilators not excluded |

## 6. Most decisive open questions
1. **a_N or (a_p, a_q) of ONE non-CM curve.** Does any arithmetic (non-fixed-algebraic, non-congruence) use of these
   integers factor N almost everywhere? Factoring a_N (n/2 bits) is the obvious first step. CM-PAIR shows that the
   archimedean structure of two such integers can suffice.
2. **A source of trace numbers.**
   - A polynomial-time method producing a_N(E) at composite N would factor N if it worked for the CM curves
     y² = x³ − Dx (CM-PAIR).
   - For non-CM families, ordered traces would factor N heuristically. Congruence decoding of a_N alone gives only QR
     bits (SC).
3. **A deterministic one-curve reduction** (EXT, TW and J0 are randomized). The obstacles are least witnesses for
   non-multiplicative character sums and the choice of F.
4. **Fixed-ℓ counts:** their computational relation to QR_N. There is no reduction either way beyond §2.7.

## 7. Reproduction and deviations
```
.venv/bin/python -m pytest -q tests/test_tracedata.py           # 24 tests, <1 s
.venv/bin/python -m src.tracedata.run_experiments e1            # 17 s
.venv/bin/python -m src.tracedata.run_experiments e2            # <1 s
.venv/bin/python -m src.tracedata.run_experiments e3a           # 67 s
.venv/bin/python -m src.tracedata.run_experiments e3b           # 5 min
.venv/bin/python -m src.tracedata.run_experiments e3c           # 36 s
.venv/bin/python -m src.tracedata.run_experiments e3d           # 73 s
.venv/bin/python -m src.tracedata.run_experiments e4            # <1 s
.venv/bin/python -m src.tracedata.run_experiments e5            # <1 s
.venv/bin/python -m src.tracedata.run_experiments posthoc       # PH1-PH4, ~1 min (regenerated after the unordered fix)
.venv/bin/python -m src.tracedata.run_experiments posthoc5      # PH5, CM pairs
.venv/bin/python -m src.tracedata.run_experiments posthoc6      # PH6, Theorem EXT + null
```
- **D1.** The end-to-end decoders test up to 4 CRT lifts (`max_lifts`). This is not in the prereg and is harmless:
  W = 6469693230 > p already, and every candidate is verified exactly.
- **D2.** E3b/E3d were designed on the false premise that products identify s mod ℓ. They identify only the square
  class (Theorem SC, found from E3b's plateau). P3c/P3d/P3f pass as registered, but E3d is uninformative.
- **D3.** The first E4 run crashed on a JSON serialisation of numpy booleans. It was rerun unchanged otherwise.
- **D4.** E4's curve C5 = y² = x³ − x + 1 was fixed in the prereg text before any run (commit f2df008).
- **D5.** `galois.law(..., 'unordered')` double-counted diagonal pairs in v1 (cx2 #4). The fix leaves identifiability
  unchanged and changes PH2's unordered row only.
- **D6.** `harness.fam_traces` had a vacuous assert, now a comment. Good reduction holds because
  |4A³ + 27B²| < 2^25 < p for every prime used. A direct count found 0 bad-reduction pairs among the 3.6·10⁶ (prime,
  curve) pairs of E3c/E3d.

## 8. Codex rounds
- **cx1 (blind)** supplied: Prop PROMISE, Prop NA, the M/φ incomparability, Prop QA, Theorem J0, the a = b = 1 case,
  the Hasse-box count, the Bach-bound counterexample, and literature corrections (DPS 1.2, Kunihiro thesis).
- **cx2 (adversarial; MAJOR REVISION)**:

| cx2 item | Response |
|---|---|
| 1 Degree-5 algebra: everywhere reduction from M for every curve | Adopted as Theorem EXT. Re-implemented and checked (PH6, null). §0, §1, §2.9, §5, §6 rewritten |
| 2 SC headline omits the scalar hypothesis; Sym² also invariant; attainability needs full images | Hypothesis stated everywhere; Sym²/BC remark added; attainability separated |
| 3 GC is a full-Haar estimate; entanglement; "cannot be the route" too strong | Relabelled; entanglement terms added; restricted-model wording |
| 4 Unordered law bug | Fixed (D5); PH2 regenerated; normalisation test added |
| 5 MIX condition; a = b = 1 degeneracy; LINK double-zero; M is randomized | Fixed (§2.3, §2.2, §1) |
| 6 MC: family dependence, correlated cells, "strongly supported" | Disclosed; reworded to "toy agreement"; Castryck–Hubrechts cited exactly |
| 7 Sample complexity vs exponential list search; M needs Ω(ℓ⁵ log ℓ) | Adopted; v1's ℓ⁴ claim withdrawn |
| 8 E3b plateau vs Bayes ceiling | Fixed |
| 9 E4 explanation | Replaced with cx2's simulation and the critic's 30000-prime check; refutation retained |
| 10 PT does not force non-pure algorithms | Corrected |
| 11 DEC scope | Narrowed (§2.1, §5) |
| 12 IL quantifiers; NA wording; odd ℓ; TOR torsion order | Fixed |
| 13 J0, QA, SEP scope | Fixed |
| 14 Priority wording | Kept "no prior statement found by us"; GM §3 credited for EXT |

## 9. Reviews

**Opus proof critic (v1): MAJOR PROBLEMS, all text-fixable; all applied in v2.** It confirmed DEC, BC3, LINK (split
primes), MIX, PROMISE, TOR, CM-QR, SC (under its hypothesis), IL, NA, QA, J0 and every recomputed number.
1. CRITICAL: "a_N does not close" was false. Two fixed CM curves factor N. → Prop CM-PAIR, re-verified (PH5).
2. GC ignores entanglement. → §2.8.
3. SC headline drops its hypothesis. → §0, §1, §2.6.
4. DEC scope. → §2.1, §5.
5. MC constant. → §2.8.
6. Limits are finite-ℓ values. → §2.8 table.
7. LINK at split primes. → §2.2.
8. MIX hypothesis. → §2.3.
9. a = b = 1 paraphrase. → §2.3.
10. NA wording. → §2.7.
11. Minor wording. → TOR, IL, randomized poly time, DEC/M phrasing, PT-J needs non-square Δ.
12. Vacuous assert; twist duplicates. → D6, §2.8.

**Opus literature critic (v2): OK WITH FIXES; all applied in v2.1.**
- It found no prior statement of EXT's uniform form, CM-PAIR, SC, or the multi-curve decoding.
- Fixes:
  - GM §3 described exactly;
  - Dieulefait–Urroz added as the explicit "open" statement;
  - DPS's y² = x(x² − b) family added as the CM-PAIR antecedent;
  - Jacobsthal / Berndt–Evans–Williams / Ireland–Rosen attribution;
  - BC published version, and Rem 2.6 as the PROMISE antecedent;
  - OU for the M = N case;
  - BLGHT Cor 8.9;
  - Castryck–Hubrechts exact error term;
  - the likelihood formula labelled KNOWN (CH Thm 1(ii), Gekeler 2003);
  - Bach 1990 precise;
  - Kunihiro thesis confirmed.
