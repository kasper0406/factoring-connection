# nativemod: coefficients and Hecke traces of FIXED-level modular forms at a composite index (phase 29)

Agent `nativemod`, 2026-10-27. **Version v2.1** (second Opus proof critic applied, §8). v2: after the codex adversarial round (cx2, MAJOR REVISION) and the Opus
proof critic (ACCEPT WITH FIXES), with every required fix applied (§8). The Opus literature critic's fixes are in
v2.1 (§8).
- Pre-registration: `notes/nativemod_prereg.md` (commit `9acb3cc`, before any experiment).
- Code: OPERATIONAL `src/nativemod/{reductions,ecalg}.py`. These take N and the stated datum only; an AST test
  forbids PARI, harness and factorisation imports and calls. HARNESS `src/nativemod/{harness,run_experiments}.py`.
- Tests: `tests/test_nativemod.py` (28 pass, < 1 s).
- Raw data: `experiments/raw/nativemod_{e1,e2,e3,e3c,posthoc}.jsonl`.
- GPT-6.1-sol: blind round `codex/nativemod_1.{md,out}` (**cx1**); adversarial round `codex/nativemod_2.{md,out}`
  (**cx2**).
- Labels:
  - KNOWN (with reference);
  - PROVED (proof here; routine unless said);
  - HEURISTIC, OPEN, experimentally supported;
  - REFUTED (prediction), UNVERIFIED.
  - The strongest novelty wording used is "no prior statement found by us".

Notation:
- N = pq, with p < q odd primes coprime to the fixed level M and to the fixed small primes of the form; s = p + q.
- f is a fixed newform of weight k and nebentypus χ, with Galois representations ρ_{f,λ}.
- g_r = ρ̄_{f,λ}(Frob_r) and a_r = a_f(r).

---

## 0. Verdict (plain language)

**We found no way in which fixing the level helps.**
- Every member of the family whose value we can show determines the factors is factoring-equivalent. The best
  native algorithms we know for those members cost N^{c}.
- Some members carry only partial information and are not known to be factoring-equivalent:
  - a weight-1 Eisenstein coefficient gives one quadratic-residue bit;
  - 11a1's coefficient gives s mod 5;
  - one lattice count gives one linear relation.
- **The member whose status is genuinely open is one coefficient a_f(N) of one non-CM cusp form.** For it we found no
  algorithm, no hardness proof, and no model barrier that covers every route. Each natural route is ruled out only
  within a stated model.

1. **Edixhoven–Couveignes at a composite index (Q1).**
   - Every step except one transfers unchanged. The torsion algebra of J₁(ℓ)[ℓ] (or of E[ℓ]) is a FIXED algebra and
     reduces mod N without difficulty. The step that breaks is "Frobenius".
   - "The Frobenius trace at N" has five readings, and they separate cleanly:
     - **The componentwise Frobenius**, as an operator, as its pair of conjugacy classes, or as its trace on
       H¹(E, O) over Z/N: each **factors N**.
       - Prop FROB (from `nonsmooth`; applies when (det-character/ℓ) can be made −1).
       - Prop CLS: the determinants of the two classes.
       - Prop CRT-T: the CRT trace gives a_r mod r for the larger prime r and Coppersmith finishes, under the
         balance hypothesis 2√p + 2√q < min(p, q). Under that hypothesis it fails exactly when a_p = a_q; without
         it the exception set is {s | a_r − a_s} (v2.1, second proof critic).
     - **The literal power map x ↦ x^N** is not a ring map and need not permute roots. Its read-outs are **PRIM
       programs** in extension components.
       - Element tests fall under leakmap's NE1.
       - Trace and Jacobi read-outs lie in PRIM's OPEN fragment.
       - Our pre-registered null search there, on three specific read-outs (E3; 4 × 20 000 runs, 12 targets), found
         nothing: max |z| = 2.37, also after full gcd checking (E3c).
     - **a(N) itself** is neither. It is N^{(k−1)/2}·2∫_{SU(2)} tr(θ_p·kθ_qk⁻¹)dk, the trace of a *convolution* of
       two Frobenius classes (Prop CONV). There is no "Frobenius of N", so E–C's output (one conjugacy class) has no
       composite analogue.
   - No (c)-type candidate was found.

2. **Eisenstein series (Q2).**
   - For weight k ≥ 2, one coefficient of the eigen-Eisenstein series E_k(χ, χ′) at N factors N deterministically
     (Prop EIS-k; 720/720). Some linear combinations of them lose this information.
   - Native algorithms exist, all exponential:
     - Õ(N^{1/2}) by the hyperbola method;
     - N^{0.49+o(1)} (Polymath);
     - plausibly Õ(N^{1/3}) by the convex-hull method.
   - The q-expansion, the Hecke action via Merel's matrix families, Fourier extraction and the Ramanujan-sum series
     cost about N or need moduli ≥ p.
   - The coefficient sequence of a modular form is not P-recursive.
   - Weight-1 Eisenstein coefficients are O(1)-bit QR-type data (r₂(N)).

3. **Fixed-level Hecke action (Q3).**
   - **No "fast exponentiation" exists in the literal sense (Prop LRS, PROVED).** Powering a fixed element of a fixed
     finite-rank algebra over a number field gives a sequence that is eventually periodic modulo every prime λ ∤ S, i.e. an
     aggregate. σ(N) mod ℓ and a_f(N) mod λ are not aggregates.
   - The integer Hecke matrix of T_N on a space with an Eisenstein part has σ_{k−1}(N) as an eigenvalue, so it is
     factoring-equivalent.
   - **In every integral-weight reduction we know, the factor information enters through a cross term** (Prop CROSS):
     a(N) of ρ₁ ⊕ ρ₂ contains a_p(ρ₁)a_q(ρ₂) + a_p(ρ₂)a_q(ρ₁). This one mechanism covers:
     - σ(N): 1 ⊕ ε^{k−1};
     - Bach–Charles: ρ⊗ρ = Sym²ρ ⊕ det;
     - tracedata's CM pairs: ψ ⊕ ψ̄ over K;
     - Brandt row sums;
     - the Eichler–Selberg divisor term.
   - A single non-CM coefficient has no infinite-order abelian summand on any open subgroup.
   - This is a taxonomy, not a theorem. cx2 notes that half-integral weight escapes it: one value of r₃(N) is a class
     number.

4. **Brandt matrices and quaternary lattices (Q4).**
   - Each quaternary theta coefficient is (Eisenstein: a σ-type divisor sum) + (cusp part).
   - At level 11 we proved r_L(n) = (12/5)σ(n) + c_L·a_{11a1}(n) (11 ∤ n), with c_L = 8/5, −12/5, 18/5 for the three
     lattices L1, L2, L3. The proof uses the Sturm bound on M₂(Γ₀(11)); we also checked it directly for n ≤ 3000.
   - So **two of these lattice counts at N factor N, and the single non-CM coefficient is the theta *difference*
     a_{11a1}(N) = (r_{L1}(N) − r_{L2}(N))/4**.
   - One count alone gives only a linear relation. That relation factors very unbalanced N (cx2), but not balanced N.
   - Counting algorithms:
     - Õ(N) for any fixed quaternary form (cx1);
     - Õ(N^{1/2}) for class-number-one forms, via divisor summatory functions;
     - Barvinok's polytope algorithms do not apply: the integer hull of a 4-ball of radius √N has ≍ N^{6/5} vertices.

5. **Shortlist (Q5).** Is "the N-th coefficient of a fixed modular form in poly(log N)" possible?
   - **(b) Equivalent to factoring** for:
     - eigen-Eisenstein series of weight ≥ 2, θ⁴ and class-number-one quaternary lattices;
     - Brandt rows, and integer Hecke matrices on spaces with Eisenstein part;
     - enough lattice counts of a genus to separate the Eisenstein part (two at level 11);
     - the Eichler–Selberg class-number moments of the dimension-zero weights.
     - Pairs of CM coefficients factor N on a positive-density set: p ≡ q ≡ 1 (mod 4) with D a non-residue
       (tracedata).
   - **(c) Genuinely open** for one coefficient of one non-CM cusp form. No reduction from factoring is known, and no
     native algorithm faster than factoring is known. The routes ruled out, each in its own model:
     - fast exponentiation (Prop LRS);
     - P-recursive recurrences;
     - aggregates;
     - the E–C transfer (FROB, or PRIM with an OPEN fragment that is empirically null on the read-outs tested);
     - lattice-by-lattice counting, where two counts already give σ(N);
     - analytic extraction (sister agent `nativeana`).
   - **First question.** For a fixed non-CM curve written as a fixed combination of theta series of a fixed genus,
     compute that combination at an unfactored balanced N = pq in time below L_N[1/3].
     - At level 11 the combination is the difference (r_{L1}(N) − r_{L2}(N))/4. 11a1 leaks s mod 5 through its
       rational 5-torsion.
     - 37a1 is the torsion-free analogue, using the level-37 Brandt genus.
     - The constraint is to compute the cuspidal combination without separately computing enough counts to recover
       σ(N). This is a target, not a theorem.
   - **Hardness of the generalisation.** With level and weight in unary, the problem is polynomial-time given a
     factoring oracle (E–C Thm 15.2.1, level 1, deterministic; Bruin and Javanpeykar for general level, GRH, Las
     Vegas). So NP-hardness would imply NP ⊆ ZPP^{FACT} (under GRH). That is not believed, but it is not proved
     impossible either. The lattice-as-input version is #P-complete (Charles 2007), but it leaves the fixed-level
     family.

Subjective: the probability that this family yields sub-L[1/3] factoring is < 1%. The open item is a hardness
question about a(N), not a speed route.

---

## 1. Q1 — Edixhoven–Couveignes at a composite index

### 1.1 Step by step

**How E–C works at a prime** (Ann. Math. Studies 176, 2011).
- They compute a(p) for prime p as follows. For each small prime ℓ ≤ O(log p), compute the mod-ℓ representation ρ̄_ℓ
  as a polynomial P_ℓ ∈ Q[x]. Its roots correspond to the nonzero points of V_ℓ ⊂ J₁(ℓ)[ℓ]; in weight 2 they
  correspond to E[ℓ], which is Schoof's ψ_ℓ.
- Read the conjugacy class of Frob_p from F_p[x]/(P_ℓ), then CRT the traces.
- Bruin, Zeng, Tian, Derickx–van Hoeij–Zeng and Mascot change how P_ℓ is computed, not the last step.

| step | at a prime p | at N = pq | status |
|---|---|---|---|
| P_ℓ over Q (fixed ℓ, level, weight) | poly(ℓ), independent of p | identical | fine (fixed data) |
| torsion algebra A = Z[1/S][x]/(P_ℓ) ⊗ Z/N | étale F_p-algebra | ≅ A_p × A_q, with components F_{p^f}, F_{q^{f′}}; a proper gcd with a denominator or disc is a factor | fine |
| "Frobenius" on A | x ↦ x^p, a ring automorphism | componentwise (Frob_p, Frob_q): producing it factors N when FROB(b) applies; literal x ↦ x^N: §1.2 | **breaks** |
| read the class / trace | class of Frob_p, trace t_p | class *pair* {C_p, C_q}: factors N (Prop CLS); a(N) mod λ = t_pt_q | — |
| CRT over ℓ ≤ O(log N) | a(p) | from traces: the unordered pair {a_p, a_q}, or a(N) | tracedata OPEN |

**FROB transfer for E–C's algebra (PROVED; nonsmooth Prop FROB(b), (d)).**
- The sign of g ∈ GL₂(F_ℓ) acting on F_ℓ² ∖ {0} is Legendre(det g). Reason: diag(g, 1) has ℓ cycles of length ℓ − 1,
  and SL₂(F_ℓ) is perfect for ℓ ≥ 5. The proof critic brute-forced all of GL₂(F_ℓ) for ℓ = 3, 5, 7: 0 mismatches.
- So (disc P_ℓ / N) = (χ(N)N^{k−1}/ℓ), up to a fixed sign, when the residue field is F_ℓ. For a larger residue field
  use the quadratic character of F_λ^×.
- When some ℓ makes this −1, the componentwise Frobenius on (Z/N)[x]/(P_ℓ) factors N. This always happens for even
  k with χ = 1, i.e. (N/ℓ) = −1.
- It never happens for odd k with χ(N) = 1. There FROB(b) does not apply, and Prop CLS is the route.

### 1.2 Prop POW: what the literal power map is (PROVED; cx1, cx2)

**The map.**
- On a component F_{p^f} of A_p, Pow_N = Frob_p ∘ Pow_q. It equals some Frob_p^r iff pq ≡ p^r (mod p^f − 1).
- It is not additive in general: in F_5, (1 + 1)^{35} = 3 ≠ 2.
- So x ↦ x^N is not an algebra endomorphism and need not permute the roots of P_ℓ.

**Its natural read-outs:**
- **Element tests**, e.g. P_ℓ(x^N) = 0, or x^N = h(x) for a fixed h. These are one-layer exponent-gate element tests
  in a fixed algebra. leakmap Lemma NE-C / Theorem NE1 cover them also in extension components (leakmap §3.5).
  - That bound is conditional on the small-order events (cx2): it is a deviation bound, not unconditional
    negligibility for every base.
  - E3: S4 = [ψ_ℓ(x^N) = 0] never fired in 80 000 runs.
- **Coordinate tests.** Tr_{A/(Z/N)}(x^N), the trace of multiplication by x^N, is ≡ s_q(P_ℓ) (mod p) and
  ≡ s_p(P_ℓ) (mod q). Here s_m is the m-th power sum of the roots.
  - PROVED: Σα^{pq} = (Σα^q)^p, and s_q is p-integral. P_ℓ/ℓ need not be integral, e.g. s_7(ψ₃/3) = −91/3 for 37a1
    (cx2), but ℓ is a unit at p. Tested.
  - So the "trace of x ↦ x^N" is the Perrin-type value s_N mod N, whose p-component is a power sum at the
    *partner's* exponent.
  - Coordinate tests, and Jacobi symbols of such values, are the OPEN fragment of PRIM: extension components, and
    Jacobi symbols of deterministic gate values (leakmap §3.5, prim §6–§7).
- **The literal Hasse invariant** [x^{N−1}] F^{(N−1)/2} mod N, for y² = F(x) with F cubic.
  - It is ≡ H_p·c_q (mod p), with H_p ≡ a_p and c_q = [x^{q−1}]F^{(q−1)/2} ∈ Z unreduced.
  - PROVED: F^{(N−1)/2} = F^{(p−1)/2}·(F^{(q−1)/2})^p mod p, and the degree bound leaves only i = p − 1, j = q − 1.
    Checked at 4 pairs here and 5 by the critic.
  - It mixes the partner's integer into the p-component. It is not a trace of Frobenius.

### 1.3 Five meanings of "the trace of Frobenius at N"

| # | quantity | what it is | factors N? |
|---|---|---|---|
| T1 | a(N) = a_pa_q | trace of g_p ⊗ g_q on ρ ⊠ ρ (cx1); = N^{(k−1)/2}·2∫_{SU(2)} tr(θ_p kθ_q k⁻¹)dk (Prop CONV) | OPEN (tracedata) |
| T2 | {tr g_p, tr g_q} mod λ, all λ up to O(log N) | gives the integers {a_p, a_q} by CRT of symmetric functions | OPEN (tracedata) |
| T3 | class pair {C_p, C_q}, all λ up to O(log N) | contains det C_r = χ(r)r^{k−1} | **yes** for k ≥ 2, deterministic (Prop CLS) |
| T4 | t_N = CRT(a_p mod p, a_q mod q) ∈ Z/N | componentwise Frobenius on H¹(E_{Z/N}, O) (Hasse–Witt over Z/N) | **yes** for balanced N unless a_p = a_q; unbalanced: unless s divides a_r − a_s (Prop CRT-T) |
| T5 | Tr_A(x^N) = s_N(P_ℓ) mod N | §1.2 | no mechanism known; E3 null |

**Prop CONV (PROVED; Schur).**
- Statement: let ρ: G → GL_d(C) be an irreducible continuous representation of a compact group, with character χ and
  Haar probability dk. Then ∫_G χ(x·kyk⁻¹) dk = χ(x)χ(y)/d.
- *Proof:* ∫ρ(kyk⁻¹)dk commutes with ρ(G), so it is the scalar χ(y)/d by Schur. Multiply by ρ(x) and take traces. ∎
- It is applied to complex compact groups only. An ℓ-adic image carries no usable Haar integral, and v1's "fix
  Q̄_ℓ ≅ C" is withdrawn (cx2, critic).
  - *Sato–Tate form* (p ≠ q, p, q ∤ M). For trivial χ write a_r = r^{(k−1)/2}·tr θ_r with θ_r ∈ SU(2) a conjugacy
    class (Deligne). With nebentypus, θ_r ∈ U(2) with det θ_r = χ(r), and the U(2) form below applies. Then
    **a(N) = N^{(k−1)/2}·2∫_{SU(2)} tr(θ_p·kθ_qk⁻¹) dk.** Equivalently, ∫_{U(2)} tr(XkYk⁻¹)dk = tr X·tr Y/2 (cx2).
  - *Artin form.* For an irreducible Artin representation of a finite group, the same identity holds with the average
    over G.
- Meaning: the product set C_p·C_q is a union of many classes, with different traces in general. a(N) is the trace of
  the *convolution* of two Frobenius class measures, not of an element.
- For a 1-dimensional ρ the product set is the single class Art((N)). That is the abelian aggregate case (nonsmooth
  §5.4).

**Prop CLS (PROVED, elementary; repaired after cx2 and the critic, using cx2's version).**
- Setting: weight k ≥ 2, χ of order h, e := h(k − 1) > 0; p, q ∤ level, and p, q larger than every ℓ used (small
  factors fall to trial division). From the class pair at λ, take e₁ := det(C_p)^h + det(C_q)^h ≡ p^e + q^e (mod λ).
- Decoding:
  - e₁ = p^e + q^e is a rational integer < 2N^e, so its residues at several λ over the same ℓ are redundant. CRT over
    DISTINCT rational primes ℓ (one λ | ℓ each) with ∏ ℓ > 2N^e + 1 recovers it exactly: O(e log N) primes of size
    O(e log N), by the prime number theorem.
  - p^e and q^e are the roots of X² − e₁X + N^e. Exact e-th roots give p.
- So a composite-index E–C that output Frobenius classes would be a polynomial-time factoring algorithm.
- Weight 1 (e = 0) is excluded: there det C_p = χ(p) carries only character values.
- v1 decoded s mod ℓ from the determinants one ℓ at a time. That fails for nontrivial χ (cx2: k = 2, χ = (·/5), ℓ = 7,
  (3, 17) vs (11, 151)) and for odd k.
- The code (`class_pair_factor`) covers χ = 1, k = 2, the case tested.

**Prop CRT-T (PROVED; standard Coppersmith; corrected exception set after cx2 and the critic).**
- Setting: E/Q fixed, N = pq of good reduction, t ≡ a_p (mod p), t ≡ a_q (mod q). Write s < r for the smaller and
  larger prime factor (s here is not p + q). Use the centred representative of t: if a_p < 0 then t = N + a_p.
- Coppersmith is run modulo the LARGER factor r, with root x₀ = a_r, |x₀| ≤ 2√r. A split then needs
  gcd(t − a_r, N) = r, i.e. s ∤ (a_r − a_s).
- **Balance hypothesis:** 2√s + 2√r < s (roughly s > (4N)^{1/3}). Then |a_r − a_s| < s, so the exception set is
  exactly **{a_p = a_q}**: there t is the small integer a_p and carries only the common value.
  - This includes the CM Blum case a_p = a_q = 0 (tracedata PROMISE).
  - For non-CM E it has density zero (heuristically O(p^{−1/2})).
- **Without balance** the exception set is {s | a_r − a_s}. Counterexample (second proof critic, 37a1): s = 1009,
  a_s = −47, r = 254899, a_r = 962, a_r − a_s = 1009 = s, N = 257193091, t = 962 and gcd(t − a_r, N) = N. v2's
  "fails exactly when a_p = a_q" was false for unbalanced N.
- Under balance and a_p ≠ a_q:
  - Apply May's theorem (May 2010, Thm 6: divisor b ≥ N^β, monic degree δ, |x₀| ≤ c_N N^{β²/δ}) with b = r.
  - β is an input to that theorem. Try β₀ = log₂R/log₂N over the O(n) dyadic bins R ≤ r < 2R (cx2). This gives
    N^{β₀²} ≥ R^{1/2}/O(1), so O(1) windows per bin.
  - Lattice dimension Θ(log N) for the polynomial-time statement.
- E2b:
  - 37a1: 80/80; 11a1: 80/80.
  - 32a2 (y² = x³ − x): 53/80. All 27 failures have a_p = a_q = 0, and none of the successes do.
  - Null control t′: 0/80.
  - cx2 and the critic produced nonzero collisions: 37a1 (1013, 1097), a = 36; 11a1 (1009, 1019), a = −10; three more
    for 37a1 at 20 bits. All give None, as now stated. Unit test added.
- The toy uses fixed d = 9, m = 4 and the balance q/p < 4/3, with ≤ 38 windows at n = 128 (rule 8: toy).
- No earlier exact statement found by cx2 or us.

### 1.4 Placement (the question as posed)

The natural Z/N analogue of the E–C mechanism is:
- **(a) a PRIM program** if "Frobenius" means x ↦ x^N: one exponent gate in a fixed algebra with components F_{p^f},
  f ≥ 2.
  - Its element-test read-outs are NE1.
  - Its trace and Jacobi read-outs lie in PRIM's unclassified extension-component fragment. E3 tested three of them,
    not the fragment.
- **(b) a FROB instance** if it means the componentwise Frobenius, its classes or its CRT trace. Each factors N, with
  the stated exceptions.
- **(c) nothing found.**

For growing ℓ the algebras form a family rather than one fixed algebra (cx1).

---

## 2. Q2 — Eisenstein coefficients

**Prop EIS-k (PROVED; routine generalisation of Bach–Miller–Shallit for σ and nonsmooth Prop MF′).**
- Setting: E_k(χ, χ′) has a(n) = Σ_{d|n} χ(n/d)χ′(d)d^{k−1}, with k ≥ 2 and real characters of conductor dividing c.
  N = pq is coprime to c. For k = 2, χ = χ′ = 1, use the holomorphic level-2 modification, since E₂ is quasimodular
  (cx2).
- Then a(N) = χ(N) + χ′(N)N^{k−1} + u·p^{k−1} + v·q^{k−1}, with u = χ(q)χ′(p) and v = χ(p)χ′(q), both in {±1}.
- For each of the ≤ 4 sign patterns, consider v X² − R X + u N^{k−1} = 0 with R = a(N) − χ(N) − χ′(N)N^{k−1}. For the
  true pattern it has the root X = q^{k−1}. Take an exact (k−1)-th root and test divisibility.
- This is deterministic and never degenerate, because v ≠ 0.
- Complex characters: enumerate the finitely many root-of-unity patterns (cx2).
- **Scope (corrected).** The statement is for the eigen-series E_k(χ, χ′). A fixed linear combination
  const + U p^{k−1} + V q^{k−1} can have U = V = 0 on some residue patterns. Example (critic, cx2):
  E₃(χ₋₄, 1) − E₃(1, χ₋₄) vanishes at every N ≡ 1 (mod 4).
- Weight 1, real characters: a(N) ∈ {0, ±4}. For (1, χ₋₄) it is r₂(N)/4, which decides only ε (both split vs both
  inert).
- E2a: 720/720 (six families × {64, 128, 256} bits × 40).

**Native algorithms for σ(N) (all exponential; costs are best-known upper bounds, not lower bounds).**

| route | cost | status |
|---|---|---|
| σ(N) = S(N) − S(N−1), S(X) = Σ_{d≤X} d⌊X/d⌋, hyperbola blocks | Õ(N^{1/2}) | PROVED (cx1) |
| same, Polymath affine-progression decomposition (Math. Comp. 81 (2012), Prop 2.2) | N^{0.49+o(1)} | KNOWN (cx1) |
| same, convex-hull walk (Sladkey, arXiv 1206.3369, for Σd(n)) plus polynomial floor sums | Õ(N^{1/3}) expected | UNVERIFIED (extension to the weighted sum not checked) |
| r₄(N) by enumeration; r₄ = Σ_t r₃(N − t²) | Θ(N); ≈ 2√N class numbers | KNOWN |
| q-expansion of θ⁴ or E₂ up to index N | Θ(N) | output size |
| Hecke action, Merel's family S_n (Merel 1994, §3.1 Prop 18; the set used by Stein GSM 79 §8.3.3) | ≍ (12 log 2/π²)σ(n)log n matrices | KNOWN (cx1, cx2) |
| Hecke action, Merel's X_n = {ad − bc = n, a > b ≥ 0, d > c ≥ 0} (Merel Prop 20) | (3/π²)p log²p + O(p log p) at primes (cx2's proof; measured 0.40 → 0.378 for p ≤ 10⁵) | PROVED (cx2) / experimentally supported |
| Fourier extraction a(N) = ∫ f(x + iy)e(−N(x+iy))dx | Ω(N) samples for the trapezoid rule (aliasing), or precision e^{−2πNy} | → `nativeana` |
| Ramanujan series σ_{k−1}(N) = N^{k−1}ζ(k)Σ_c c_c(N)/c^k | the terms with c < p coincide with those of a number with no small prime factors: moduli ≥ p are needed | PROVED for k = 2 (nonsmooth §4.1); same argument for all k |
| fixed linear recurrence in n for a(n) | impossible: not P-recursive | PROVED below |

**Not P-recursive (PROVED from KNOWN facts; argument after the critic).**
- Let f ≢ 0 be a holomorphic modular form of weight k > 0 for a finite-index Γ. The function φ = y^{k/2}|f| is
  Γ-invariant.
- Suppose f extended holomorphically across a real interval I. Then φ → 0 on approach to I.
- The Γ-orbit of every z ∈ H accumulates at every real point, because the limit set is P¹(R). Hence φ ≡ 0 and f ≡ 0.
- So |q| = 1 is a natural boundary. This includes cusp forms.
- A D-finite power series has only finitely many singularities, among the roots of its leading coefficient (Cauchy
  existence theorem; Stanley 1980).
- Hence (a(n))_n satisfies no linear recurrence with polynomial coefficients. Bostan–Gaudry–Schost √N-step
  evaluation does not apply.
- Scope (cx2): this excludes only P-recursive recurrences. It does not exclude nonlinear identities or other
  isolated-coefficient algorithms.

---

## 3. Q3 — fixed-level Hecke action

### 3.1 Known routes to T_N on S_k(Γ₀(M)) or M_k(Γ₀(M))

| route | uses the factorisation? | best known cost |
|---|---|---|
| definition: ψ(N)-many coset matrices (a b; 0 d), ad = N | yes (divisors) | — |
| q-expansion: a_{T_N f}(n) = Σ_{d|(n,N)} d^{k−1}a(nN/d²) | yes, and coefficients up to ≍ N | Θ(N) |
| Merel's matrix families S_N / X_N on Manin symbols | **no** | ≍ σ(N)log N / ≍ N log² N matrices |
| Hecke multiplicativity T_N = T_pT_q, with E–C/Bruin for T_p | yes | poly(log N) |
| trace formula (Eichler–Selberg) | yes (divisor term), or ≈ 2√N class numbers | ≈ √N·L[1/2] |
| Brandt matrices (Pizer) | no (theta counts) | Õ(N) in general; Õ(N^{1/2}) for class-number-one orders (§4) |
| PARI `mfheckemat(n)` | apparently evaluates coefficients up to ≍ n (inferred from memory use) | E1b: stack overflow near n ≈ 7·10⁷ at 512 MB |

So "every known route needs the factorisation" is literally false (cx1). Merel sums, Brandt counts and q-expansions
are native. The accurate statement: **every known poly(log N) route goes through T_p and T_q, and every native route
costs N^{Ω(1)}.** On M₂(SL₂Z) = 0 there is nothing to compute; the Eisenstein argument needs a nonzero Eisenstein
component.

### 3.2 Prop LRS: no fast exponentiation (PROVED; (iii) strengthened following cx2)

**Statement.** Let K be a number field, S a finite set of primes, and u(n) = c^T Aⁿ b with A ∈ M_r(O_{K,S}) and
b, c ∈ O_{K,S}^r. Equivalently, u is a linear recurrence sequence. This covers any fixed polynomial in the entries of
gⁿ for g in a fixed finite-rank O_{K,S}-algebra, since sums and products of LRS are LRS. Then:
- (i) For every prime λ ∤ S of O_K, u(n) mod λ is eventually periodic.
- (ii) On semiprimes, σ(N) mod ℓ (ℓ ≥ 3) is not eventually a function of N mod m, for any m.
- (iii) Let f be a cuspidal newform and λ a prime with ℓ odd such that ρ̄_{f,λ} is absolutely irreducible. This holds
  for all but finitely many λ (Ribet), and also for weight 1. Then on semiprimes a_f(N) mod λ is not eventually a
  function of N mod m.

Hence no LRS satisfies u(N) = σ(N), or u(N) = a_f(N), for all sufficiently large semiprimes. (For a_f, u lives in
O_{K,S} and a_f in K_f: pass to the compositum K·K_f and choose a prime Λ | λ with Λ ∤ S; v2.1.)

*Proof.*
- **(i)** A mod λ lies in a finite set, so its powers are eventually periodic.
- **(ii)**
  - Put L = lcm(m, ℓ). Choose λ₀ ∈ (Z/L)^× with λ₀ ≢ 1 (mod ℓ), and choose p₀, q₀ with p₀ ≢ λ₀⁻¹q₀ (mod ℓ).
  - Dirichlet gives arbitrarily large primes p ≡ p₀, q ≡ q₀, p′ ≡ λ₀p₀, q′ ≡ λ₀⁻¹q₀ (mod L).
  - Then pq ≡ p′q′ (mod m), while (p′ + q′) − (p + q) ≡ (λ₀ − 1)(p₀ − λ₀⁻¹q₀) ≢ 0 (mod ℓ).
- **(iii)** (cx2's argument.)
  - Let G be the image and G′ = [G, G]. Not every element of G′ has trace 2. Otherwise, since det = 1 on G′, every
    element of G′ would be unipotent. Then G′ would be a normal ℓ-group, which acts trivially on an absolutely
    irreducible module in characteristic ℓ, so G′ = 1 and G would be abelian. That is a contradiction, since an
    abelian group has no absolutely irreducible 2-dimensional representation.
  - Pick h ∈ G′ with tr h ≠ 2 = tr 1. Every abelian quotient of G, in particular the one cut out by
    K_λ ∩ Q(ζ_m), kills G′.
  - So Chebotarev in K_λ(ζ_m) gives distinct primes p, p′, q, all coprime to Mℓ, with p ≡ p′ (mod m), g_p = h and
    g_{p′} = 1 up to conjugacy, and g_q = 1. Then tr g_q = 2 ≠ 0, because ℓ is odd; this is where "ℓ odd" is used.
  - Then N = pq ≡ p′q = N′ (mod m), and a(N) − a(N′) ≡ (tr h − 2)·2 ≢ 0 (mod λ). ∎

**Scope.**
- It rules out "N as an exponent in a fixed global algebra": square-and-multiply on a fixed matrix, exponential
  polynomials in N, and generalised power sums.
- It does not rule out digit-recursive (k-regular) schemes. Mod ℓ those reduce to finite automata reading the digits
  of N and predicting s mod ℓ on semiprimes, which is OPEN in general (leakmap §6 covers only AC⁰ and
  Thue–Morse/Rudin–Shapiro).
- It does not rule out arithmetic mod N.
  - Powering in T ⊗ Z/N is PRIM.
  - Corrected after the critic (v1 misapplied PT-X, which covers fixed-output trees, not value outputs): the
    Eisenstein eigenvalue mod N is σ_{k−1}(N) ≡ 1 + s^{k−1} (mod N), since the cross terms of (p + q)^{k−1} are
    multiples of pq. Unit test.
  - For k = 2 this is 1 + s and factors N.
  - For k ≥ 4, recovering s ≈ √N from s^{k−1} mod N is an RSA-type root extraction beyond Coppersmith's
    N^{1/(k−1)}. It is not known to factor N.
  - So "integer Hecke matrix ≡ FACT" does not transfer to "Hecke matrix mod N" for k ≥ 4.
- It does not rule out other computations of the N-dependent correspondence (cx2), e.g. a Brandt module mod an
  auxiliary prime projected to its cuspidal part.

### 3.3 Prop CROSS: the ingredient that makes coefficients factor N

**Identity (PROVED).** For N = pq and compatible systems ρ₁, ρ₂,
a_N(ρ₁ ⊕ ρ₂) = a_N(ρ₁) + a_N(ρ₂) + X, with X = a_p(ρ₁)a_q(ρ₂) + a_p(ρ₂)a_q(ρ₁).

| data | decomposition | cross term X |
|---|---|---|
| σ_{k−1}(N) (Eisenstein) | 1 ⊕ ε^{k−1} | p^{k−1} + q^{k−1} |
| Bach–Charles a(N²) = a_N(Sym²f), χ = 1 | ρ ⊗ ρ = Sym²ρ ⊕ det | q^{k−1}(a_p² − p^{k−1}) + p^{k−1}(a_q² − q^{k−1}); so a(N²) = a(N)² + N^{k−1} − S, S = q^{k−1}a_p² + p^{k−1}a_q² (E1b(iii): 100/100) |
| CM pairs (tracedata CM-PAIR) | ρ_f restricted to G_K = ψ ⊕ ψ̄, ψ(𝔭)ψ̄(𝔭) = p | ψ(𝔭)ψ̄(𝔮) + ψ̄(𝔭)ψ(𝔮), via sums of two squares |
| Brandt row sums, genus averages (§4) | Brandt module = Eisenstein ⊕ cusp | Eisenstein eigenvalue σ(N) |
| Eichler–Selberg divisor term | hyperbolic (Eisenstein-type) term | −½Σ_{dd′=N} min(d, d′)^{k−1} = −(1 + p^{k−1}) |
| residual: 11a1 mod 5 | ρ̄₅^{ss} = 1 ⊕ ε | a(N) ≡ (1 + p)(1 + q) (mod 5), giving s mod 5 (tracedata Lemma TOR; unit test) |

**Observation (HEURISTIC; a taxonomy, not a theorem; scope narrowed after cx2 and the critic).**
- Claim: among reductions from integer coefficient data of integral-weight compatible systems at a squarefree index,
  every one we know gets the factors from a cross term with a one-dimensional summand of infinite order. That summand
  is a power of the cyclotomic character, or a Hecke character of an imaginary quadratic field, possibly after
  restriction to an open subgroup. Residual (mod λ) Eisenstein congruences give the finite-information version.
- Outside this scope:
  - CRT-T uses archimedean size and Coppersmith on a value in Z/N.
  - Bach–Charles uses the index N², which is not squarefree.
  - Half-integral weight escapes the taxonomy (cx2). One coefficient r₃(N) = 24h(−N) (N ≡ 3 mod 8, squarefree) gives
    a class number, and class-group powering then finds an ambiguous form and factors N, under GRH.
- A single non-CM coefficient has no such summand on any open subgroup. tracedata Prop DEC shows that (N, a_p, a_q)
  satisfy no fixed algebraic relation with p. This is where the open case sits.

**Comparison with Bach–Charles (KNOWN, arXiv:0708.1192 §2, Thm 2.2).**
- Ingredients:
  - the two oracle values a(N) and a(N²);
  - the Hecke relation at p², which is where det = χ(p)p^{k−1} enters (the Sym² ⊕ det cross term above);
  - non-vanishing of a(p) for density-one primes (Serre; Lehmer for τ);
  - Deligne's bound, for the integrality step (cx1).
- An oracle for Tr T_{rN} and Tr T_{rN²} at several fixed r on a fixed nonzero S_k(1) also suffices. Invert a fixed
  system to get a_i(N) and a_i(N²), then apply BC (cx1). This is a multi-index oracle.
- One value a(N), or one trace Tr T_N: OPEN.

### 3.4 Eichler–Selberg at level 1 (KNOWN identity; E1a verified)

**The identity.**
- C_k(N) := Σ_{t²<4N} P_k(t, N)H(4N − t²), with P_k = [z^{k−2}](1 − tz + Nz²)⁻¹ (Zagier, *Modular forms of one
  variable*, §6.3; Cohen 1975).
- Then Tr(T_N | S_k(1)) = −½C_k(N) − ½Σ_{dd′=N} min(d, d′)^{k−1} for k ≥ 4. For N = pq the divisor term is
  −(1 + p^{k−1}).

**Dimension-zero weights** k ∈ {4, 6, 8, 10, 14}.
- C_k(N) = −2(1 + p^{k−1}), and C₂(N) = 2(N + q) (Kronecker–Hurwitz). Each moment factors N deterministically.
- E1a: 0 mismatches over all 785 odd semiprimes ≤ 4000, for every k. The critic recomputed with independent class
  numbers.
- Best known native cost: ≈ 2√N Hurwitz class numbers.

**Nonzero spaces** k ∈ {12, 16, 18, 20, 22, 26}.
- C_k = −2Tr T_N − 2(1 + p^{k−1}). E1a: 0 mismatches over 293 semiprimes ≤ 1500.
- The moment carries p^{k−1}, but it is masked by Tr T_N, which is of the same order for balanced N.

**Charles 2006** (Ramanujan J. 11, 221–224).
- It computes τ(N) through this formula after factoring N.
- Bounds: N^{1/2+ε} randomized under GRH, and N^{3/4+ε} unconditional (cx1).
- It is not a native composite-index algorithm.

---

## 4. Q4 — Brandt matrices and quaternary theta series

**Structure (KNOWN: Siegel–Weil, Eichler's basis problem, Pizer's Brandt matrices).**
- For a positive quaternary lattice L of level M and N coprime to M, r_L(N) = E_L(N) + C_L(N). Here E_L is a
  combination of σ-type divisor sums, and C_L is a cusp form.
- The genus average, with weights 1/|O(L)|, is E_L.
- For an Eichler order the Brandt matrix B(N) has the constant vector as an eigenvector with eigenvalue σ(N) for
  (N, M) = 1. These are its row sums: the p + 1 lines of F_p² at each prime (cx1).

**Level 11 made explicit (E1c, PH; critic and cx2 re-verified for n ≤ 500).**
- The lattices:
  - L1 = (x² + xy + 3y²) ⊕ (z² + zw + 3w²), |O(L1)| = 32: the norm form of the maximal order with 4 units.
  - L2 = [[4,−1,−1,2],[−1,4,2,−1],[−1,2,4,0],[2,−1,0,4]] (Gram of 2Q), |O| = 24: the connecting-ideal lattice.
    r_{L2}(1) = 0, so it is not an order (cx2).
  - L3 = [[2,−1,−1,−1],[−1,2,0,0],[−1,0,8,−3],[−1,0,−3,8]], |O| = 72: the norm form of the order with 6 units.
- They are pairwise non-isometric. All three have the same anisotropic discriminant form on (Z/11)², so they lie in
  one genus (Nikulin; critic). We do not claim the genus is complete.
- Identities, for 11 ∤ n where needed:
  - (i) a_{11a1}(n) = (r_{L1}(n) − r_{L2}(n))/4;
  - (ii) the Siegel average over {L1, L2, L3} equals (12/5)σ(n);
  - (iii) r_{L1} = (12/5)σ + (8/5)a, r_{L2} = (12/5)σ − (12/5)a, r_{L3} = (12/5)σ + (18/5)a;
  - (iv) the Brandt rows 3r_{L1} + 2r_{L2} = 12σ(n) and 3r_{L2} + 2r_{L3} = 12σ(n) (critic).
- (iii) and (iv) are **PROVED for all n**: θ_L lies in M₂(Γ₀(11)), which is 2-dimensional with Sturm bound 2, and
  each identity has been checked to n = 3000.
- Registered E1c(ii) used only L1 and L2, with weights 1/|O(L)|: **REFUTED as registered** (2312/2727 failures). L3
  was found post hoc by a wider search.
- Consequences (PROVED from (iii)):
  - Two counts at N give σ(N), hence the factors, because the c_L are distinct. This is because the cusp space is
    1-dimensional. In general, (dim of cusp part + 1) counts with independent cusp vectors are needed.
  - One count is one linear relation:
    s = (5/12)r_L(N) − (N + 1) − (5c_L/12)·a_pa_q, with |(5c_L/12)a_pa_q| ≤ 6√N.
    - It approximates s, and hence q, within O(√N). By Coppersmith that factors unbalanced N with
      q ≥ N^{1/√2+ε} (cx2).
    - For balanced N it is not known to factor N (OPEN, like tracedata's mixed data). For L1 it also gives a(N) mod 3,
      since 5r_{L1} ≡ 8a (mod 3).
  - The non-CM coefficient is the theta difference.
- Class-number-one lattices: r_{Z⁴} = 8σ (odd n), r_{D4} = 24σ (odd n; the Hurwitz order, level 2) and
  r_{A2⊕A2} = 12σ (3 ∤ n). E1c(iii): 0 failures for n ≤ 3000. Each single count is FACT.

**Counting algorithms (best known upper bounds).**
- Any fixed positive quaternary form: Õ(N), by meet-in-the-middle on two-coordinate sums (cx1).
- Class-number-one forms: via σ, Õ(N^{1/2}) or N^{0.49} (§2).
- Slicing to ternary: r_L(N) = Σ_t r_{L′}(N − Q(t)) is ≈ 2√N ternary counts. Each is a class number only when the
  ternary genus has no cusp part (cx2).
- Gauss-circle-type exact counting is Õ(R^{2/3}) for planar curves (Hensley, Pacific J. Math. 166 (1994), per cx1).
  It does not lift to dimension 4: slices cost ≈ N^{4/3} (cx1).
- **Barvinok** (fixed dimension, polynomial in the facet description) applies to polytopes.
  - The natural polytope with the same lattice points, the integer hull of the 4-ball of radius √N, has
    ≍ R^{d(d−1)/(d+1)} = N^{6/5} vertices (Bárány–Larman, Math. Ann. 312 (1998); Balog–Bárány 1991 for d = 2).
  - This obstructs the polytope route only, not every curved-body algorithm (cx2).
  - Barvinok flags the 4-ball as the factoring obstruction ("Something curvilinear?", IAS/Park City notes 2007; per
    cx1).
- Analytic (theta-function / circle-method) evaluation: the precision wall, sister agent `nativeana`.

---

## 5. Q5 — verdict, first question, hardness of the generalisation

| member | status | best known native cost |
|---|---|---|
| E_k(χ, χ′) coefficient, k ≥ 2; r₄; σ_{k−1} | ≡ FACT (Prop EIS-k) | Õ(N^{1/2}), N^{0.49} |
| weight-1 Eisenstein (r₂/4) | QR-type, O(1) bits | one QR bit (OPEN) |
| integer Hecke matrix on a space with Eisenstein part; Brandt row; genus average | ≡ FACT (eigenvalue σ_{k−1}(N)) | Õ(N^{1/2}) at best |
| Hecke matrix mod N, k ≥ 4 | gives 1 + s^{k−1} mod N: RSA-type, not known FACT | — |
| enough lattice counts of a genus (two at level 11) | ≡ FACT | Õ(N) |
| Hurwitz moments C_k(N), dim S_k(1) = 0 | ≡ FACT | ≈ 2√N class numbers |
| pairs of CM coefficients (y² = x³ − Dx family) | ⇒ FACT on p ≡ q ≡ 1 (4) with D a non-residue (tracedata CM-PAIR); a_N = 0 on Blum N | — |
| (a(N), a(N²)), non-CM | ⇒ FACT, density one (Bach–Charles) | — |
| **a_f(N), one non-CM newform; equivalently one cuspidal theta combination** | **OPEN both ways** (11a1 leaks s mod 5) | factor + E–C; natively Õ(N) |
| one Tr T_N on a nonzero S_k(1) | OPEN (single-coefficient type) | — |
| one lattice count, multi-class genus | OPEN for balanced N; FACT for q ≥ N^{1/√2+ε} | Õ(N) |
| weight-1 non-abelian (η(z)η(23z)) | O(1) bits; Prop IL: not FACT-hard under N-local reductions unless factoring is easy | OPEN |
| half-integral weight r₃(N) | a class number (CLASS tier) | L[1/2] (GRH) |

**(a), (b) or (c)?**
- **(b)** for the FACT rows: poly(log N) evaluation ⟺ polynomial-time factoring, randomized or deterministic as
  stated.
- **(c)** for a single non-CM coefficient. No barrier covers all algorithms. In their models, the following are ruled
  out:
  - fast exponentiation in a fixed global algebra (Prop LRS);
  - P-recursive recurrences (§2);
  - aggregates N mod m (Prop LRS(iii));
  - the E–C transfer: FROB (operator, classes, CRT trace) or PRIM, whose element tests are NE1-bounded and whose
    coordinate/Jacobi fragment is OPEN but null on the three read-outs E3 tested;
  - lattice-by-lattice counting, where enough counts give σ(N);
  - Fourier/analytic extraction (precision wall; `nativeana`).

**First question.**
- For a fixed non-CM curve E written as a fixed linear combination Σ c_L θ_L of theta series of one genus, with
  Σ c_L·(genus weight) = 0: is there an algorithm computing Σ c_L r_L(N) for unfactored balanced N = pq in time below
  L_N[1/3], or polynomial?
  - 11a1: the difference (r_{L1} − r_{L2})/4, which leaks s mod 5.
  - 37a1: the torsion-free case, the level-37 Brandt genus. Explicit lattices are not computed here.
- The constraint is to avoid computing enough individual counts to separate σ(N). This is a target, not a theorem;
  computing one count natively at balanced N would itself be new.
- This is cx1's question (a native a_N(E)) in lattice form. The companion hardness question is tracedata §6.1: does
  a_N alone factor N?

**Natural generalisation (KNOWN + corollary; literature per cx2).**
- With the level and weight in unary and the form given by Hecke data, coefficients are computable in polynomial time
  given the factorisation of the index:
  - E–C Thm 15.2.1: level 1, deterministic; unconditional for fixed weight, and uniform in the weight under GRH.
  - Bruin, Publ. Math. Besançon 2011, Thm 1.1: levels n₀n₁ with n₁ squarefree; GRH, Las Vegas.
  - Javanpeykar, ANT 8 (2014), Thm 5.0.1: general levels.
- So the problem lies in ZPP^{FACT} (under GRH). NP-hardness would imply NP ⊆ ZPP^{FACT}. That is not believed, but
  it is not proved impossible.
- The lattice-as-input variant is #P-complete under polynomial-time Turing reductions (D. X. Charles, "Counting
  lattice vectors", JCSS 73 (2007) 962–972, Thm 3.1). Its level grows with the input, so it is outside this family.

---

## 6. Experiments (pre-registered; `notes/nativemod_prereg.md`, commit 9acb3cc)

Peak RSS ≤ 256 MB. Jobs ran one at a time: E1 14 s, E2 17 s, E3 47 s, E3c 52 s.

| exp | registered prediction | result |
|---|---|---|
| E1a | 0 mismatches | 0 over (785 × 6 + 293 × 6) checks; **holds** |
| E1b | 0 mismatches | 4 levels × 25 N × 3 identities: 0 failures; **holds** (D1: N ≤ 5000, not 10⁴) |
| E1c(i) | 0 | 0 / 3000; **holds** |
| E1c(ii) | 0 (two lattices, weights 1/\|O\|, constant 12/5) | 2312 failures; **REFUTED as registered** (genus has a third class; PH) |
| E1c(iii) | 0 | 0 for Z⁴, D4, A2⊕A2; **holds** |
| E1d | descriptive | #X_n/((12 log 2/π²)σ(n) log n) ranges over 1.77–3.37 for n ∈ [300, 600] (wrong reference set; see §2) |
| E2a | "240/240 per family" | 120/120 per family (720/720). The registered count was a miscount of the design (D2); **holds** |
| E2b | 80/80 per curve; null 0 | 37a1 80/80; 11a1 80/80; **32a2 53/80: REFUTED as registered** (all 27 failures have a_p = a_q = 0; should have been predicted); null 0/80 |
| E2c | 240/240 | not run as an experiment; only a 10-pair unit test exists (D5) |
| E3 P3a | all \|z\| < 4 | max \|z\| = 2.37 over 4 cells × 3 statistics × (ℓ + 3) targets (cluster-robust: 2.36); **holds** |
| E3 P3b | S4 ≤ 0.1% | 0 / 80 000; **holds** |
| E3 P3c | ≤ 3 accidental splits per cell, every value gcd-checked | v1 code checked only the traces and S4 (cx2), so the registered gcd coverage was **not met**. Corrected rerun E3c: 113, 1, 454, 3 splits against null expectations ≤ 167, 0.9, 562, 3.0. **REFUTED as registered** for n = 32: the prediction ignored intermediate coefficients. All ACCIDENTAL. E3c max \|z\| = 2.36 |
| E3 P3d | positive control \|z\| > 10 | 141.4; **holds** but degenerate: the control statistic equals the target, so z = √n. It checks wiring, not power |

**Post hoc (PH; `run_experiments posthoc`).**
- Three-lattice level-11 genus: the Siegel average, the Brandt 3:2 relation and the per-lattice decompositions all
  have 0 failures; the three lattices are pairwise non-isometric.
- Merel's X_p at primes p ≤ 8009: #X_p/(p log² p) decreases 0.51 → 0.40, and cx2 measured 0.378 at p = 100003.
  This matches cx2's proof of (3/π²)p log²p + O(p log p).

**Reading of E3.** Three specific deterministic read-outs (Tr y, Tr xy, Tr y² Jacobi) of y = x^N in Schoof's algebra
for 37a1, ℓ ∈ {3, 5}, show no association with a_N mod ℓ, ε or (p/ℓ), beyond N mod 4ℓ. The scope is limited (cx2):
- it is not the whole fragment;
- ℓ does not grow;
- the two ℓ cells reuse moduli;
- at n = 32 there are 1491 distinct primes behind 19 802 distinct N.

## 7. Candidate cards (rule 9)

| mechanism | model | outcome |
|---|---|---|
| E–C with x ↦ x^N on the torsion algebra | PRIM, extension components | element tests NE1; three trace/Jacobi read-outs null in E3/E3c; OPEN fragment |
| E–C with componentwise Frobenius / class pair / CRT trace | FROB | factors N (FROB(b) when applicable; CLS for k ≥ 2; CRT-T unless a_p = a_q): a reduction, not a native route |
| fast exponentiation of a fixed Hecke element | global LRS | **refuted for this model** (Prop LRS) |
| fixed recurrence for a(n) | P-recursive | **refuted** (natural boundary) |
| Merel sums for T_N | native, explicit | ≍ σ(N)log N or N log² N matrices: exponential |
| θ⁴ / σ via summatory functions | native | Õ(N^{1/2}), N^{0.49}: exponential |
| lattice counts in a genus | native | enough counts give σ(N); one count is OPEN for balanced N; best Õ(N) |
| Barvinok for 4-ellipsoids | polytope methods | integer hull has N^{6/5} vertices: route not applicable |
| Eichler–Selberg moments | class numbers | ≈ 2√N class numbers; FACT |

## 8. Reviews

**cx2 (GPT-6.1-sol adversarial round, MAJOR REVISION; 13 items). All applied in v2.**

| cx2 item | response |
|---|---|
| 1 POW (s_q integrality, "need not", NE1 conditional, FROB per algebra) | applied (§1.2, §1.1) |
| 2 CONV application to ℓ-adic images unjustified | withdrawn; Sato–Tate / U(2) and Artin forms (§1.3) |
| 3 CLS wrong for χ ≠ 1 (counterexample) and odd k | replaced by cx2's p^e + q^e repair (§1.3) |
| 4 CRT-T: equal nonzero traces (37a1 (1013, 1097); 11a1 (1009, 1019)); unknown β; toy parameters | exception set {a_p = a_q}; dyadic β bins; unit test (§1.3) |
| 5 EIS-k scope (combinations; E₂ quasimodular; weight-1 values) | applied (§2) |
| 6 LRS (iii) gap for larger residue fields and inner twists | replaced by cx2's absolute-irreducibility argument (§3.2) |
| 7 CROSS taxonomy: r₃(N) counterexample in half-integral weight | scope narrowed; listed (§3.3) |
| 8 natural-boundary sketch incomplete for cusp forms | critic's invariant-function proof (§2) |
| 9 lattice consequences: one count ≠ σ; unbalanced partial reduction; L2 not an order; ternary slices; Barvinok scope | applied (§4) |
| 10 "not NP-hard" unjustified; E–C/Bruin/Javanpeykar status; M₂(SL₂Z) = 0 | applied (§0, §3.1, §5) |
| 11 E3 gcd coverage (111 / 417 unchecked coefficient splits); positive control degenerate; scope | E3c rerun with full gcd checks; P3c reported as REFUTED; scope stated (§6) |
| 12 Merel: X_n vs S_n (Prop 18 vs Prop 20); (3/π²)p log²p proved | §2 and §3.1 corrected |
| 13 missing native routes (Jacobi forms, theta lifts, Shimura, CM values, Borcherds, p-adic, Brandt mod ℓ) | §8.1 |

**Second Opus proof critic (run by the orchestrator on v2; CLS, LRS(iii), CONV, CRT-T). All applied in v2.1.**
- Verdicts: CLS ACCEPT WITH FIXES; LRS(iii) ACCEPT WITH FIXES; CONV ACCEPT WITH FIXES; CRT-T MAJOR as stated. It
  found no counterexample to CLS, LRS(iii) or the CONV identity.
- (1) MAJOR, CRT-T: the exception set {a_p = a_q} is complete only under balance (2√s + 2√r < s), with Coppersmith
  run modulo the larger factor. Unbalanced counterexample on 37a1, s = 1009, r = 254899 → §1.3 rewritten, §0 updated.
- (2) CLS: CRT over distinct rational ℓ, not ∏ Nλ. (3) CLS: p, q ∤ level and > every ℓ used.
- (4) LRS: compositum K·K_f with Λ ∤ S; "λ ∤ S" in §0. (5) LRS(iii): g_q = 1; p, p′, q distinct and coprime to Mℓ.
- (6) CONV: SU(2) only for trivial χ; U(2) with det θ_r = χ(r) otherwise; p ≠ q, p, q ∤ M.
- (7) CRT-T: centred representative of t.
- The code's toy decoder assumes q/p < 4/3, which lies inside the balance hypothesis; E2b is unaffected.

**Opus proof critic (ACCEPT WITH FIXES; 15 findings). All applied in v2.**
- It independently verified POW, the Hasse identity, the GL₂ sign (ℓ = 3, 5, 7), EIS-k, LRS (i)–(ii), CROSS, the
  Eichler–Selberg moments, all level-11 identities, the genus via discriminant forms, and every number in §6.
- It reproduced E3 cell ℓ = 3, n = 32 exactly with independent code.
- Majors, all applied: (1) the CRT-T equal-trace exception; (2) CLS; (3) the CONV application; (4) the LRS/PT-X
  scope (mod N gives 1 + s^{k−1}).
- Minors, all applied:
  - (5) LRS (iii);
  - (6) β;
  - (7) FROB sign with χ and odd k;
  - (8) EIS-k for combinations;
  - (9) CM pairs and "any two lattices";
  - (10) Observation scope and the 11a1 mod-5 leak;
  - (11) N² in the multi-index route;
  - (12) the lattice wording, (iii) proved via Sturm, genus justification, second Brandt row, no completeness claim;
  - (13) costs as upper bounds;
  - (14) the E2c and gcd-coverage deviations and the E1d wording;
  - (15) the degenerate positive control, the natural-boundary proof, and the headline.

### 8.1 Further native routes (cx2 item 13; assessment, not proofs)

- **Shimura / theta lifts / Jacobi forms.**
  - Shimura coefficients relate c(tn²) to Σ_{d|n} μ(d)χ_t(d)d^{λ−1}a(n/d), so they again need divisors.
  - Waldspurger's c(N)² ∝ L(f ⊗ χ_N, 1) is a central value of conductor ≍ N², not a(N) (`nativeana`).
  - Jacobi-form theta decompositions repackage the same coefficients.
- **CM values, Borcherds products, Rankin–Cohen.**
  - CM-value traces need class data of growing discriminant.
  - Borcherds and Rankin–Cohen produce divisor transforms or additive convolutions.
  - No isolated-coefficient shortcut found.
- **p-adic / overconvergent methods at a fixed auxiliary prime.** They control U_ℓ and ℓ-power indices. No recurrence
  for N coprime to ℓ follows.
- **Brandt module mod an auxiliary prime ℓ.**
  - B(N) on the fixed supersingular set is still an N-dependent correspondence. It counts N-isogenies, i.e. theta
    coefficients, and its row sums are σ(N).
  - Projecting to the cuspidal part discards the Eisenstein row sum. A direct evaluation of the projected operator,
    without the row sum, is exactly the first question; LRS does not exclude it.

## 9. Reproduction and deviations
```
.venv/bin/python -m pytest -q tests/test_nativemod.py                      # 28 tests, < 1 s
PYTHONPATH=. .venv/bin/python -m src.nativemod.run_experiments e1         # 14 s
PYTHONPATH=. .venv/bin/python -m src.nativemod.run_experiments e2         # 17 s
PYTHONPATH=. .venv/bin/python -m src.nativemod.run_experiments e3         # 47 s (as registered code, v1 gcd coverage)
PYTHONPATH=. .venv/bin/python -m src.nativemod.run_experiments e3c        # 52 s (full gcd checks, v2)
PYTHONPATH=. .venv/bin/python -m src.nativemod.run_experiments posthoc    # ~10 s
```
- **D1.** E1b drew its semiprimes from N ≤ 5000 instead of 10⁴. PARI's `mfheckemat(N²)` overflows a 512 MB stack for
  N ≳ 8000.
- **D2.** E2a's registered "240/240 per family" was a miscount. The registered design (3 sizes × 40) has 120 per
  family.
- **D3.** E3's registered z is the stratified CMH z. A prime-cluster-robust z was added post hoc. cx2 notes it
  overcounts duplicated dyads; it changes nothing here.
- **D4.** Smoke runs (`--quick`) were made before the frozen runs and were overwritten by them.
- **D5.** E2c (registered 240/240) was implemented only as a unit test (10 pairs, all factored).
- **D6.** The registered "every computed Z/N value is gcd-checked" was not implemented in v1 (`ecalg` checked the
  traces and S4 only). v2 adds `ecalg.Checker` and the rerun E3c. The v1 E3 file is kept.
