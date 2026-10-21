# abcsunit: factoring relations as S-unit equations — abc, Evertse, Lagarias–Soundararajan

Agent `abcsunit`, phase 19, 2026-10-21. Revised (v3) after codex round 2 (MAJOR REVISION), the Opus proof critic (ACCEPT WITH FIXES) and the Opus
literature critic (OK WITH FIXES); §7 records the reviews.

**Files.**
- Code: `src/abcsunit/`.
  - `core.py`: exact smooth-number enumeration, residue buckets, relations, F₂ algebra, square extraction,
    N-free primitive S-unit triples, radicals; DIAGNOSTIC |⟨S⟩| in (Z/N)^* (uses p, q).
  - `model.py`: the counting models M0/M1 (shift singular series; floating point for model evaluation only).
  - `run_experiment.py` (E1–E4), `score.py` (P1–P9), `posthoc.py` (PH1, PH2; disclosed post hoc).
- Tests: `tests/test_abcsunit.py` (11).
- Pre-registration: `notes/abcsunit_prereg.md`, commit 604ab31, made before any registered run.
- Config `experiments/configs/abcsunit.json`; raw `experiments/raw/abcsunit_{E123,E4,posthoc,design_smoke}.jsonl`.
- Codex (GPT-6.1-sol): blind round `codex/abcsunit_1.{md,out}`; adversarial round `codex/abcsunit_2.{md,out}`.

**Labels.** PROVED (proof here) · PROVED from KNOWN (cited theorem as black box) · KNOWN (cited) · model-PROVED
(exact consequence of a named model) · HEURISTIC (assumption named) · measured · open · assessment. Novelty is never
claimed; "no prior statement found by us" is the most we say.

**Notation.** x = log N; L[α, c] as in README. N is squarefree with least prime factor P⁻(N) > y (the case of
interest: N = pq balanced; a factor-base prime dividing N is found by a gcd and ends the computation).
S = S_y = {ℓ ≤ y}, s = |S| = π(y), Γ_S = {±∏ℓ^{e_ℓ}}, θ(y) = Σ_{ℓ≤y} log ℓ, R_S = e^{θ(y)}.
e(u) ∈ Z^s is the exponent vector of an S-unit u. φ_N : Z^s → (Z/N)^*, e ↦ ∏ℓ^{e_ℓ} mod N; H = im φ_N;
**Λ_N = ker φ_N** (the relation lattice; full rank, index |H|).

---

## 0. Verdict

**Plain language.**
- **The cited abc and S-unit bounds establish no general barrier for factoring.** The reason is one line: every
  relation contains N in its radical, and N is squarefree, so abc never excludes a primitive relation of height
  ≤ C·N^{1+ε}, and above that it can only exclude relations whose multiplier kN or smooth side is N^ε-powerful
  (Prop ABCV). Such powerful families are genuinely constrained (e.g. 2^a + N = 3^b forces 3^b ≪ (6N)^{1+ε}; an
  exact high-quality example is 2¹⁹ + 7153 = 3¹², quality > 1.2, codex), but no known method needs them. Counting bounds (Evertse) apply only to the triply-smooth
  subfamily and allow exponentially more relations (3·7^{2s+3}) than the ≈ s a method needs; the relations factoring
  actually uses (k free) are the points of a full-rank lattice and are infinite in number (Prop EV).
- **"Solve S-unit equations with an N-congruence" is not a new target.** Listing smooth pairs u ≡ v (mod N) is,
  verbatim, the doubly-split product-listing problem Π of pilist / nonexact, i.e. listing short vectors of Λ_N
  (Prop SU). It is SUFFICIENT for heuristic polynomial-time factoring with a polylog factor base (smoothenum T8),
  with one extra condition that we state explicitly (Lemma FR): the dependencies must reach a square root of unity other than
  ±1 — guaranteed (probability exactly 1/2) when the found relations generate Λ_N up to finite odd index; the exact
  criterion is Lemma FR(2), and maximal parity rank alone is not enough.
- **Known S-unit algorithms have no N-dependent step that helps.** de Weger's p-adic reduction linearises the
  *higher-level* conditions ℓ^m | x − y by ℓ-adic logarithms; the level-one condition N | x − y that relations
  need is exactly the lattice Λ_N, and computing that lattice (or its determinant) factors N under GRH (Prop DW).
  This is hardness of the complete kernel, not of producing selected short relations. The real-place reduction with target log N
  is Schnorr's CVP method (agent `schnorrcvp`). Modularity methods need level ≥ N.
- **Lagarias–Soundararajan (LS) does not transfer usefully, and the reason is not where one might expect.**
  - The literal "a, b, c all smooth with c = kN" is empty (c has the prime factors of N).
  - The two correct N-analogues are *doubly smooth* (u ≡ v, k free) and *triply smooth* (u, v, k smooth).
  - For the doubly smooth one, which is all factoring needs, **existence is free**: pigeonhole gives
    unconditional abundance at any polylog smoothness y = (log N)^κ (κ > 1) for heights N^{κ/(κ−1)+ε}
    (Prop EX); a GRH asymptotic follows from LS's own Euler-product bound (codex, eq. 4.3). So LS's hard
    analytic content (minor arcs for a + b = c) is not needed, and adds nothing.
  - What fails is **constructivity**: LS is a circle-method count with no search procedure in it, and no
    polynomial-time construction of polylog-smooth abc triples is known even without N. A constructive
    N-version with matching heights, output count and prefix cost would be exactly Π-doubly-split, i.e. the open
    primitive itself. LS thus re-describes the target; it
    is not a route to it.
- **Experiments (pre-registered; 2 of 11 registered predictions failed, both diagnosed post hoc).** The
  N-congruence count follows the LS-type singular-series model M1 within ±10% for κ ≥ 1.5 or K = 3/2 (median ratio
  0.96 over all 51 runs), but M1 **overestimates by 10–45% at the smallest smoothness κ = 1.25** (median −29% at
  K = 9/8, ≈ 5.6σ under the clustering noise); the one out-of-band run (R1 = 0.56) is a point of that bias. Local
  factors are essential (the plain equidistribution model is off by up to ×8). F₂-rank is exactly π(y) − 2
  (quadratic-character obstruction) in every run with enough relations. Four runs (all with ≤ 32 relations) had
  Z-rank < s and hence relations not of odd index; three of them had no fruitful dependency — as Lemma FR predicts.
  None of 5.9·10⁶ primitive relations has abc quality > 1 (max 0.83).

**Main results.**

| Result | Status |
|---|---|
| Translation: doubly smooth = Λ_N points; triply smooth = fixed-coefficient S-unit equation; one-sided (Dixon/QS) = abc-type but not S-unit; NFS = mixed S-unit congruence in K × Q (§1) | PROVED (definitions) |
| Prop ABCV: abc can only exclude relations with kN or u N^ε-powerful; fraction ≪ N^{−ε/2} + Σ_{p\|N}1/p among multipliers; abc never excludes a primitive relation of height ≤ C·N^{1+ε} | PROVED (conditional only in reading abc) |
| Cor ABC-T: under abc, primitive triply-smooth relations have height ≤ C_ε(N·R_S)^{1+ε} (N-relative LS 2011 Thm 1.1/2.1 = LS 2012 Prop 1.1) | PROVED assuming abc |
| Prop EV: #primitive triply-smooth relations ≤ 3·7^{2s+3} (any height, any N); never below the ≈ s needed; doubly smooth: infinitely many | PROVED from KNOWN (Evertse 1984) |
| Lemma SC: no set of s primes has more S-units ≤ X than the first s primes | PROVED |
| Prop SU: S-unit congruence listing = Π doubly split = short vectors of Λ_N | PROVED (definitional) |
| Lemma FR: square roots of unity reachable = image of (Λ₀ ∩ 2Z^s)/2 in H[2]; full when [Λ_N : Λ₀] is finite and odd; z-dim = ρ − r | PROVED |
| Prop DW: det Λ_N (S ⊇ primes < 2 log² N) = φ(N) under GRH, so the level-one step of de Weger's method factors N | PROVED from KNOWN (Bach 1990) |
| Prop EX: unconditional abundance of doubly (log N)^κ-smooth relations at height N^{κ/(κ−1)+ε} | PROVED from KNOWN (Ψ(x,(log x)^κ) = x^{1−1/κ+o(1)}) |
| Codex eq. (4.3): GRH asymptotic for collisions at height N^A, A > κ/(κ−2) | PROVED by codex from LS Prop 5.1; checked by us (§4.2) |
| Closure lemma (codex): the map (a, b) ↦ c = a(b+N)/(b−a) (and any map whose output vector lies in the Z-span of its inputs) cannot raise rank | PROVED |
| LS constructivity: no polynomial-time construction known; constructive N-version ≡ Π | assessment / KNOWN (absence in LS) |

**Assessment (subjective).** Probability that the abc/S-unit route yields L[α < 1/3]: ≤ 1% (ours; codex
blind: ≈ 1% over twenty years). The route contributes a clean reformulation and shows that the cited Diophantine
bounds give no general barrier; it contributes no mechanism.

---

## 1. The exact translation (Task 1)

### 1.1 Doubly smooth relations are lattice points [PROVED]

Let u, v be positive S-units, u ≠ v, u ≡ v (mod N); write v − u = kN (k ≠ 0, WLOG v > u, k ≥ 1). Put d = gcd(u, v).
Since d | kN and gcd(d, N) = 1 (P⁻(N) > y), d | k. With u′ = u/d, v′ = v/d, k′ = k/d:

  u′ + k′N = v′,  gcd(u′, v′) = 1,  and gcd(k′, u′v′) = 1

(a prime dividing k′ and u′ divides v′). The map (u, v) ↦ e(v) − e(u) depends only on (u′, v′) and is a
bijection between primitive pairs of *positive* S-units and nonzero vectors of Λ_N up to sign (relations
u ≡ −v, e.g. QS's negative values or Schnorr's w < 0, need a sign coordinate: Λ_N ⊂ Z/2 × Z^s), with u′ = ∏ℓ^{e⁻}, v′ = ∏ℓ^{e⁺} and
Σ_ℓ |e_ℓ| log ℓ = log(u′v′). So "relations with |u|, |v| ≤ X" are the vectors of Λ_N with both halves of weighted
ℓ¹-norm ≤ log X (Schnorr's prime-number-lattice norm).

Three formulations must be separated (codex §3a, independently):
1. **Doubly smooth, k free** (u, v smooth): an S-unit *congruence*, not an S-unit equation. Λ_N has full rank, so
   there are infinitely many; e.g. ℓ^{j·ord_N(ℓ)} ≡ 1 (codex). The support of the abc triple (u′, k′N, v′) is
   S ∪ P(N) ∪ P(k′) and varies with the relation.
2. **Triply smooth** (u, v, k smooth): dividing by k gives v′/k′ − u′/k′ = N, i.e. the **fixed-coefficient unit
   equation** N⁻¹X + N⁻¹Y = 1 in X, Y ∈ Γ_S (codex), equivalently an S′-unit equation for S′ = S ∪ P(N). It does
   not require knowing P(N).
3. **"a, b, c all smooth with N | c"**: empty when P⁻(N) > y (codex). This is the literal reading of "LS with
   c = kN", and it must be replaced by 1 or 2.

### 1.2 Which relations fit

| Method | Relation | Fits |
|---|---|---|
| Doubly smooth congruence | u − v = kN, u, v S-units | 1.1(1): Λ_N point |
| Schnorr 1991 | u − vN = w with u, v smooth, \|w\| small and smooth | 1.1(2): triply smooth (k = v) |
| Dixon | t² − r = kN, r smooth, t random | abc-type triple (r, kN, t²); **not** an S-unit relation. In lattice form: adjoin a free generator t, relation 2ε_t − e(r) in ker(Z^{s+1} → (Z/N)^*) |
| QS | (⌊√N⌋ + t)² − N = r | as Dixon, with k = 1 |
| NFS rational side alone | a − bm smooth | an abc identity r + bm = a only; a, b, m not smooth |
| NFS (two-sided) | φ(a − bθ) = a − bm in Z[θ]/(N, θ − m) ≅ Z/N | **not a rational S-unit relation.** It is a *mixed congruence*: α = a − bθ is an S_K-unit of O_K (its principal ideal is S_K-supported; nonmonic f adds denominator qualifications), but writing its exponent vector on fixed generators requires class group and unit data, which NFS sidesteps with ideals and quadratic characters/Schirokauer maps. α − (a − bm) = −b(θ − m), and for monic f, Norm(m − θ) = f(m), in general a *multiple* of N (codex: f = x² + 1, m = 98, N = 17·113 gives 5N), so adjoining only the degree-one primes above N does not make θ − m an S-unit. With b smooth the clean form is the fixed-coefficient unit equation (α/b)/(m − θ) − ((a − bm)/b)/(m − θ) = 1 in S_K-units of K, with m − θ kept as a coefficient (codex); Evertse's bound then applies with all places counted (3·7^{d+2t}, t = #finite places + r₁ + r₂) after removing common scaling |

**Uniform statement [PROVED, definitional].** A relation is an element of ker(Γ₁ × Γ₂ → (Z/N)^*) for finitely
generated groups Γ_i with known reductions (Γ_S; Γ_S × ⟨t⟩; Γ_S × O_{K,S_K}^*), and "size" is a height. Factoring
needs ≥ rank-many short kernel elements with fruitful dependencies (§3.2).

---

## 2. What abc and S-unit theory rule out (Task 2)

### 2.1 abc is vacuous for relations [PROVED]

**Prop ABCV.** Let (u′, k′N, v′) be a primitive relation as in §1.1 (N squarefree, P⁻(N) > y), and suppose
abc fails for it at (ε, C): v′ > C·rad(u′v′k′N)^{1+ε}, with 0 < ε ≤ 1, C ≥ 1. Put g = gcd(rad k′, N). Then
  (i) k′g/rad(k′) > C·N^ε/2, or (ii) u′/rad(u′) > C·N/2.

*Proof.* The three terms are pairwise coprime, so rad(u′v′k′N) = rad(u′)rad(v′)rad(k′N), and since N is squarefree
rad(k′N) = rad(k′)·N/g ≥ N. Also v′ = u′ + k′N ≤ 2·max(u′, k′N).
- If k′N ≥ u′: 2k′N ≥ v′ > C·rad(k′N)^{1+ε} ≥ C·rad(k′N)·N^ε, so k′N/rad(k′N) > CN^ε/2, and
  k′N/rad(k′N) = k′g/rad(k′).
- If u′ > k′N: 2u′ > v′ > C(rad(u′)·N)^{1+ε} ≥ C·rad(u′)·N, so u′/rad(u′) > CN/2. ∎

**Corollary (sparsity of the exceptions) [PROVED].** For K ≥ 1 and T ≥ 1,
#{k ≤ K : k/rad(k) ≥ T} ≤ 2ζ(3/2)·K/√T. (k/rad(k) ≤ P(k), the powerful part of k, since p^{e−1} ≤ p^e for e ≥ 2
and p^0 for e = 1; the number of powerful m ≤ t is ≤ Σ_b √(t/b³) ≤ ζ(3/2)√t, as every powerful m is a²b³; partial
summation gives Σ_{m≥T, powerful} 1/m ≤ 2ζ(3/2)/√T.) Multipliers with g > 1 (p | k′) are a fraction ≤ Σ_{p|N} 1/p.
So under (i) abc can exclude at most a fraction ≪ N^{−ε/2} + Σ_{p|N}1/p of multipliers k′ ≤ K; condition (ii)
requires the smooth side itself to have powerful part > N/2.

**Reading.** abc constrains only relations in which kN or u is N^ε-powerful. A fully general form, covering also
Dixon/QS triples (r, kN, t²) whose terms are not S-units: since every primitive triple containing the term kN has
radical ≥ N (N squarefree, coprime to the other terms), **abc never excludes a primitive relation of height
≤ C·N^{1+ε}**. Doubly smooth relations with generic multipliers are therefore unconstrained; Schnorr-type
(triply smooth) relations receive only the upper height bound of Cor ABC-T, vacuous at polylog y. In particular
abc cannot give a statement "no relation family with candidates of log-size < f(N)": the only lower bound on
relation heights is the elementary one, max(u, v) > N for distinct positive congruent integers (codex; the
"aggregate bound" of systeq §4), and abc gives only *upper* height bounds.

**Comparison with the smoothness-size principle.** That principle (candidates of log-size ≍ x^{2α}(log x)^{1−2α}
for L[α]) is a statement about the probability that a random integer is smooth (Canfield–Erdős–Pomerance, J. Number Theory 17 (1983)), plus the
assumption that a method samples candidates "randomly". It has no Diophantine content. abc neither implies nor
strengthens it; a product lister would by design violate its sampling assumption (codex §2a).

### 2.2 The N-relative Lagarias–Soundararajan bound [PROVED assuming abc]

**Cor ABC-T.** Assume abc. For every ε > 0 there is C_ε such that every primitive triply-y-smooth relation
(u′, k′, v′ all y-smooth) satisfies v′ ≤ C_ε(N·R_S)^{1+ε}, i.e. log v′ ≤ (1+ε)(log N + (1+o(1))y) + O_ε(1).
*Proof.* rad(u′v′k′N) ≤ R_S·N. ∎ (codex eq. 2.2.) This is LS 2011 Thm 1.1/2.1 = LS 2012 Prop 1.1 (abc ⇒ κ₀ ≥ 1) with N added
to the radical.

**Reading.** For y = (log N)^κ, κ > 1, the bound allows heights up to exp((log N)^κ), far beyond any polynomial
search space; it binds only for y = O(log N), where triply-smooth relations of height N^A, A > (1+ε)(1 + y/log N),
are excluded. That is an *upper* height bound. It can constrain a proposed tall triply-smooth family (codex), but
it gives no lower bound on relation sizes and no obstruction for families of height N^{O(1)} at polylog y.

**Effective versions [KNOWN].** Stewart–Yu (Duke Math. J. 108 (2001) 169–181, Thm 1) give
log c ≪ R^{1/3}(log R)³ with R = rad(abc) (the earlier Stewart–Tijdeman 1986 and Stewart–Yu 1991 bounds are
weaker). With R ≥ N this allows
c up to exp(N^{1/3+o(1)}): vacuous. Baker-type lower bounds for |u − v|, u ≠ v S-units, carry constants exponential
in s and are weaker than the trivial |u − v| ≥ N.

### 2.3 Counting bounds do not bite [PROVED from KNOWN]

**Prop EV.**
(a) *Doubly smooth.* Λ_N is a full-rank lattice, so there are infinitely many primitive doubly smooth relations;
no S-unit counting theorem applies. Under a height bound X the number of all congruent pairs is Θ(Ψ_S(X)²/|H|) in the models of §5
and ≥ ½(Ψ²/|H| − Ψ) unconditionally; the primitive ones number ≥ (Ψ/2)/Ψ(X/N) when Ψ ≥ 2|H| (Prop EX).
(b) *Triply smooth.* The map (u′, k′, v′) ↦ (X, Y) = (v′/k′, −u′/k′) is injective into the solutions of
N⁻¹X + N⁻¹Y = 1 in Γ_S (u′/k′ is in lowest terms since gcd(k′, u′) = 1). Evertse (Invent. Math. 75 (1984)
561–584, Thm 1: ax + by = 1 over a number field of degree d has ≤ 3·7^{d+2t} solutions in S-units, t = number of
places in S including the infinite ones) gives, for Q with s finite primes (t = s + 1):

  #{primitive triply-smooth relations} ≤ 3·7^{2s+3},  for every N and every height.

Evertse's theorem allows arbitrary fixed coefficients a, b ∈ K^* (literature critic: checked). (Beukers–Schlickewei,
"The equation x + y = 1 in finitely generated groups", Acta Arith. 78 (1996) 189–199: 2^{8(r+1)} for subgroups
of (K^*)² of rank r — weaker here.)
(c) *Does it bite?* A method needs about s + O(1) relations (columns). 3·7^{2s+3} > s for every s ≥ 1, so this
global ceiling never forbids a factor base of size s from supplying enough triply-smooth relations, in any regime
of s (L[β], polylog, O(log N), o(log N)). It says nothing about how many lie below the height a method can use, or
about their span (codex); and it provides no existence. Erdős–Stewart–Tijdeman (1988) conjectured an upper bound ≪ exp(s^{2/3+ε}) for
**every** set S of s primes (as stated in Konyagin–Soundararajan 2007; literature critic). Even this conjectural
ceiling, applied to S′ = S ∪ {p, q}, is ≫ s.
(d) *Measured.* In every registered run the triply-smooth primitive count exceeded π(y) (min ratio 1.12; P8).

**Lower bounds [KNOWN; relevance nil].** Erdős–Stewart–Tijdeman (Compositio 66 (1988) 37–56):
exp((4 − ε)(s/log s)^{1/2}) solutions for suitable S; Konyagin–Soundararajan (J. Number Theory 124 (2007)
193–199): exp(s^{2−√2−ε}); LS 2011 Thm 1.4 under GRH: exp(s^{1/8−ε}) for the first s primes. These concern chosen
or initial prime sets and N-free equations. They defeat a hypothetical "too few unit solutions" barrier but do not
produce N-relations.

**Lemma SC [PROVED].** For any set of s primes ℓ₁ < ⋯ < ℓ_s and the first s primes p₁ < ⋯ < p_s,
Ψ_{ℓ}(X) ≤ Ψ_{p}(X) for all X: ∏ℓ_i^{e_i} ↦ ∏p_i^{e_i} is injective and does not increase values. This controls
only the number of S-units; the injection preserves neither congruences mod N nor |H_S|, so it gives **no**
optimality statement for N-relations (codex: for N = 101·103, ord_N(3) = 1700 < ord_N(2) = 5100, so at height
3^{1700} the base {3} has a relation and {2} has none). What remains is the heuristic reading [HEURISTIC]: in the
equidistribution regime the relation count is ≈ Ψ_S²/(2|H_S|), and |H_S| can be lowered only by choosing S inside a
subgroup, which beyond Jacobi-symbol conditions requires the factorisation. The ES-T/KS gains come from N-free
arithmetic coincidences, which carry no information on N.

### 2.4 Is this known in substance?

- The abc → smoothness mechanism is LS 2011 (Thm 1.1, 2.1; unconditional Thm 2.2 from Stewart–Yu). Cor ABC-T is
  its immediate N-relative variant.
- That abc is vacuous for factoring relations because rad ⊇ N (Prop ABCV, with the explicit exception set): no
  prior statement found by us; it is elementary and probably folklore.
- That Evertse's bound applies to the fixed-coefficient equation N⁻¹X + N⁻¹Y = 1: immediate from Evertse; no
  prior use in the factoring context found by us.

---

## 3. What algorithmic S-unit solving would give (Task 3)

### 3.1 Same target as Π [PROVED, definitional]

**Prop SU.** Let SU(N, y, X, R) ask for R distinct primitive pairs (u, v) of y-smooth integers ≤ X with
u ≡ v (mod N), u ≠ v. Then SU is literally the doubly-split instance of Π (pilist §1.1, first rows; nonexact §6):
the form bN + v on the triangular region {b ≥ 1, 1 ≤ v, bN + v ≤ X} with base = primes ≤ y, v split; the
bijection is b = (u − v)/N. (On the rectangular box [1, X/N] × [1, X] heights up to 2X occur; codex.)
Equivalently: list R short vectors of Λ_N in the norm of §1.1. Requiring k smooth (triply) restricts the output set
to the fixed-coefficient unit equation; requiring c = kN smooth empties it.

**Sufficiency [PROVED conditional on the stated hypotheses; the existence parenthetical is HEURISTIC; = smoothenum
T8, restated].** Let P⁻(N) > y (else a gcd factors N) and H ⊇ (Z/N)^*[2] (e.g. H = (Z/N)^*; under GRH this
follows from Bach for y ≥ 2 log² N, i.e. κ > 2; unproved for 1 < κ ≤ 2). Suppose SU can be solved with
polynomial delay for y = (log N)^κ, X = N^K, with K > 1 and 2K(1 − 1/κ) > 1 + δ (so that, heuristically,
N^{δ+o(1)} ≥ π(y) relations exist), and its outputs generate Λ_N up to finite odd index (the sufficient condition
of Lemma FR). Then N factors in polynomial time:
s = π(y) = (log N)^{κ+o(1)} columns, R = s + O(1) relations, dense F₂ elimination poly(s), and a dependency gives
Z² ≡ 1 with Z uniform in H[2] (Lemma FR), splitting with probability 1/2 for N = pq (1 − 2^{1−ω} for
squarefree N with ω prime factors). Total time
T_SU + poly(log N). With y = L[β] the same reduction preserves the exponent β. Listing *all* solutions is
exponential (the counts are N^{Θ(1)}), so only a useful prefix can be asked for (smoothenum T8, codex §4d).

### 3.2 Lemma FR: fruitfulness is governed by Λ₀ ∩ 2Z^s, not by parity rank [PROVED]

Let Λ = Λ_N ⊂ Z^s, H = Z^s/Λ, Λ′ = {x : 2x ∈ Λ}, so Λ′/Λ ≅ H[2]. Let e₁, …, e_R ∈ Λ be relation vectors,
Λ₀ their Z-span, r the F₂-rank of their reductions in Z^s/2Z^s, ρ the F₂-rank of their classes in Λ/2Λ, and
D = {c ∈ F₂^R : Σc_i e_i ∈ 2Z^s} the dependency space. Define z : D → Λ′/Λ, c ↦ [½Σc_i e_i] (lift c to {0,1}).
The square root of unity of a dependency is Z = ∏ℓ^{−½Σc_i e_i} (mod N) (§core.square_from_dependency; Z = W/V).

1. z is a well-defined homomorphism (changing lifts changes Σc_i e_i by 2Λ₀ ⊂ 2Λ).
2. im z = ((Λ₀ ∩ 2Z^s)/2 + Λ)/Λ. (Write y ∈ Λ₀ ∩ 2Z^s as Σa_i e_i; c = a mod 2 lies in D and
   Σc_i e_i ≡ y mod 2Λ₀.)
3. **dim im z = ρ − r.** The kernel of Λ/2Λ → Z^s/2Z^s is (Λ ∩ 2Z^s)/2Λ = 2Λ′/2Λ ≅ H[2]; the image V of Λ₀ in Λ/2Λ
   has dim ρ and maps onto a space of dim r, and V ∩ kernel corresponds exactly to im z.
4. If [Λ : Λ₀] is finite and odd, then ρ = s and im z = H[2] (conversely ρ = s iff the index is finite and odd). (For x ∈ Λ′ and m = [Λ : Λ₀] odd, 2mx ∈ Λ₀ ∩ 2Z^s and
   [mx] = [x] in H[2].)
5. r ≤ s − dim(H/2H) (rows lie in (Λ + 2Z^s)/2Z^s, of codimension dim H/2H). For N = pq and H = (Z/N)^* this is
   **r ≤ s − 2**: the parity rows are orthogonal to c_p, c_q (indicators of ℓ with (ℓ/p) = −1, resp. q), because
   (uv/p) = (u/p)² = 1.

So with H[2] = {±1, ±ε} ≅ F₂² (N = pq, H = (Z/N)^*), a uniformly random dependency is fruitful with
probability exactly 1/2 if ρ = r + 2 (z uniform on H[2]); if ρ = r + 1 the probability is 0 or 1/2 according as
im z = {1, −1} or one of {1, ε}, {1, −ε}; if ρ = r it is 0. ρ = r + 2 holds when [Λ_N : Λ₀] is finite and odd.
If the parity rank is maximal (r = s − 2) and rank_Z Λ₀ < s, then ρ − r ≤ rank_Z Λ₀ − r < 2 and fruitfulness is
lost or halved. (Without r = s − 2 a Z-rank deficit need not hurt: the critic's example N = 1147 = 31·37,
S = {2, 3, 5, 7}, relations (225², 256²), (175², 972²) has Z-rank 2, r = 0, ρ = 2 and fruitful dependencies giving
31 and 37; codex's example N = 77, S = {2, 3, 5}, relations (30, 0, 0), (0, 30, 0) has Z-rank 2, r = 0 and full
image {1, 34, 43, 76}, and adding (0, 0, 30) gives a full-rank Λ₀ of *even* index 450 with full image. So odd
index is sufficient, not necessary; the exact criterion is item 2: im z ⊄ ⟨−1⟩.)
Reaching the maximal parity rank r = s − 2 is not sufficient. Codex found the extreme case independently: integer
combinations of parity-independent kernel vectors have only unfruitful dependencies (codex §3a). For one-sided
(QS/Dixon) relations each relation brings its own free generator t; this changes the sampling argument (Dixon's
fruitfulness proof uses the randomness of t) but does not remove the square-root-image requirement — e.g. the
relations (jN + 1)² ≡ 1 all give the trivial root (codex).
*Prior statement:* this is the standard structure of index-calculus square roots; the nearest prior art is the
rigorous relation-lattice analyses of Schnorr–Lenstra (1984), Seysen (1987), Hafner–McCurley (1989) and
Lenstra–Pomerance (1992) (literature critic). We state the odd-index / square-root-image criterion explicitly for
purely S-unit relations; no prior statement in exactly this form found by us. **Measured:** ρ − r = 2 ⟺ rank_Z = s in all 51 runs
(§6, PH2).

### 3.3 Known S-unit algorithms and the N-condition

**The algorithms [KNOWN].** de Weger (CWI Tract 65, 1989; also J. Number Theory 26 (1987)); Tzanakis–de Weger;
Smart (1995, 1998); the Sage solver of Alvarado–Koutsianas–Malmskog–Rasmussen–Vincent–West (2019/2020); von
Känel–Matschke (arXiv 1605.06079; Memoirs AMS 2023). Pattern: (1) height bounds — Baker/Matveev and Yu linear forms in logarithms for de Weger, Smart and
Alvarado et al.; Faltings/modularity-based bounds in terms of N_S = ∏_{ℓ∈S} ℓ for von Känel–Matschke — exponential
in s; (2) LLL reduction at the real place (log ℓ_i) and at each ℓ ∈ S (ℓ-adic logarithms); (3) sieving/enumeration
of the remaining box. None has complexity polynomial in s; even the output is exp(Ω(√(s/log s))) for some S (ES-T), exp(s^{2−√2−ε}) (KS), and
grows like 2^s empirically for the first s primes (VKM 2016 Thm A, which covers n ≤ 16: |Σ(S)| = 1, 4, 17, 63, 190,
545, 1433, …, 839638 for n = 1, …, 15, counted modulo the S₃ symmetry, total 6|Σ| − 3; our E4 reproduces n ≤ 5).

**Prop DW (where the N-condition sits) [PROVED from KNOWN].**
1. For ℓ ∈ S the ℓ-adic step encodes ℓ^m | X − Y, m large, through ℓ-adic logarithms of the ℓ-units involved
   (normalised to principal units; ℓ itself does not enter the logarithm, ℓ = 2 is handled separately, and
   inhomogeneous coefficient logarithms appear in modern versions, cf. Alvarado et al. §5.1), after the
   level-one residue/torsion condition mod ℓ is handled through |F_ℓ^*| = ℓ − 1 (known, small).
2. A factoring relation needs only the **level-one** condition N | v − u, i.e. e ∈ Λ_N. There is nothing to
   linearise: the level-one lattice is the target itself.
3. Computing a basis of Λ_N — even only det Λ_N = |H| — factors N when S ⊇ {primes < 2 log² N}: under GRH these
   primes generate (Z/N)^* (Bach, Math. Comp. 55 (1990): the least integer outside a proper subgroup is
   < 2 log² m, hence also the least prime; Lamzouri–Li–Soundararajan give < 3 log² m with a prime p ∤ m; with
   P⁻(N) > y the factor-base primes are coprime to N; literature critic),
   so |H| = φ(N) = N − p − q + 1. Conversely, computing Λ_N from p, q takes 2s discrete logarithms (of each ℓ ∈ S mod p and mod q).
   (That a multiple of φ(N), or an order computation in (Z/N)^*, factors N is classical: Miller 1976; Bach,
   "Discrete logarithms and factoring", 1984.)
   So the step de Weger handles trivially for small ℓ is, for the composite N, at least as hard as factoring —
   in the sense of computing the complete kernel (or its determinant). This does not show that every way of
   obtaining *selected* short vectors of Λ_N is factoring-hard (that is the open primitive itself; codex).
   (Unconditionally, |H| is a multiple of the order of every element of H; Miller-type splitting from it works
   when the 2-parts of the orders differ mod p and q, which fails e.g. if |H| is odd. We do not need this.)
4. The real-place step with an inhomogeneous target log N is Schnorr's lattice (EUROCRYPT 1991; Adleman 1995
   manuscript). Schnorr's 2021 claims of fast factoring by this route (IACR ePrint 2021/232, 2021/933) were not
   borne out experimentally (Ducas, github.com/lducas/SchnorrGate). Find e with
   Σe_iℓ log ℓ_i ≈ log N, producing u − vN = w small, a triply-smooth relation if w is smooth (§1.2). That is the
   remit of agent `schnorrcvp`; we do not duplicate it.
5. Modularity methods (Frey curve y² = x(x − u′)(x + k′N)) need modular forms of level ∝ rad(u′v′k′N) ≥ N, and
   von Känel–Matschke's height bounds are in terms of N_S, which here would include N: the known enumeration
   algorithms then cost N^{Ω(1)} (an explanation of their expense, not a universal lower bound).
6. Codex adds: N is a known ℓ-adic unit for ℓ ∈ S, so ℓ-adic logarithms of N are free; substituting composite N
   for a prime place is invalid; whether *every* composite-modulus logarithmic method is factoring-hard is open.

**Closure lemma (codex §4c) [PROVED].** From pairs (a, a + N), (b, b + N), c = a(b + N)/(b − a) satisfies
c + N = b(a + N)/(b − a); when integral both are smooth, but e(c) − e(c + N) = [e(a) − e(a + N)] − [e(b) − e(b + N)].
So this closure (and any construction whose output vector lies in the Z-span of the input vectors) never
enlarges Λ₀; it cannot raise ρ or r, and by Lemma FR it cannot create fruitfulness. Other rational constructions
can introduce non-cancelling exponent vectors and are not covered.

### 3.4 Consistency with `suffmap`

suffmap (notes/suffmap_results.md, E3 and row H8) files abc/S-unit counting as **unrelated** and algorithmic
S-unit solving as **sufficient** (= its doubly-smooth lister E3), citing this note; that agrees with §§2–3. Two
remarks on the interface:
- suffmap's one-line reason for abc vacuity ("rad(uv) ≤ e^{B(1+o(1))}") covers the triply-smooth case with
  B ≫ log N (our Cor ABC-T reading). For relations with k free, the operative reason is the *lower* bound
  rad(u′v′k′N) ≥ N·rad(k′)/g of Prop ABCV.
- suffmap E3(iii) proves that constructive *navigation* in Cay((Z/N)^*, S) to hidden random squares a_i² factors N
  with no heuristic. That construction attaches a fresh hidden generator to each relation, which is exactly why
  Lemma FR's odd-index condition does not arise there; for purely S-unit relations (Prop SU) it must be imposed.
  suffmap Prop DLT (a modulus-generic lister also yields prime-field discrete logs) applies verbatim to SU.

---

## 4. Lagarias–Soundararajan and smooth abc solutions (Task 4)

*Notation in §4:* H denotes a height (LS's notation); the subgroup ⟨S⟩ ⊂ (Z/N)^* of §§1–3 is written H_S here.

### 4.1 What LS prove [KNOWN; checked against the arXiv text]

- *JTNB* 23 (2011) 209–234 ("the xyz conjecture"): primitive A + B + C = 0 with all prime factors
  ≤ (log H)^κ, H = max(|A|, |B|, |C|); κ₀ := the threshold exponent.
  - abc ⇒ κ₀ ≥ 1 (Thm 1.1; more generally weak abc with exponent κ₁ ⇒ κ₀ ≥ κ₁);
  - unconditionally (Stewart–Yu) finitely many with smoothness ≤ (3 − ε) log log H (Thm 2.2);
  - GRH ⇒ κ₀ ≤ 8 (Thm 1.2); strong xyz conjecture κ₀ = 3/2 (heuristic Ψ³/H = H^{3(1−1/κ)−1+o(1)} → ∞ ⇔ 3(1 − 1/κ) > 1 ⇔ κ > 3/2);
  - Thm 1.4 (GRH): ≫ exp(s^{1/8−ε}) primitive S-unit solutions for the first s primes.
- *Proc. LMS* 104 (2012) 770–798 ("Counting smooth solutions to A + B = C"), Thms 2.1–2.3 under GRH, for
  (log H)^{8+δ} ≤ y ≤ exp((log H)^{1/2−δ}): weighted counts
  N(H, y; Φ) = 𝔖_∞(1 − 1/κ, Φ)·𝔖_f(1 − 1/κ, y)·Ψ³/H·(1 + O(log log y/log y)), the primitive version with 𝔖_f*
  and error O(Ψ³/(H(log y)^{1/4})),
  and the primitive density 1/ζ(2 − 3/κ). 𝔖_f(c, y) = Π_{p≤y}(1 + (p−1)/(p(p^{3c−1} − 1))·((p − p^c)/(p−1))³)
  Π_{p>y}(1 − 1/(p−1)²). Method: circle method over exponential sums of smooth numbers, Hildebrand–Tenenbaum
  saddle point; GRH bounds partial Euler products (their Prop 5.1) on the minor arcs. The exponent 8 is an artefact
  of the error term H^{3/4+ε}Ψ against Ψ³/H (codex: 2 − 3/κ > 7/4 − 1/κ ⇔ κ > 8).
- Later: Drappeau (Math. Proc. Cambridge 154 (2013), GRH, unweighted asymptotics; codex); Harper (Compositio 152
  (2016) 1121–1158) with Drappeau's major arcs: the asymptotic **unconditionally** for (log H)^C ≤ y ≤ H, C large
  (Harper, abstract, VERIFIED by us; codex confirms; the major-arc input is Drappeau, "Sommes friables d'exponentielles
  et applications").
- LS give no algorithm. Their explicit examples come from de Weger's exhaustive lists, Gross–Zagier singular
  moduli, Poonen–Schaefer–Stoll and the Monster.

### 4.2 Making it N-dependent

**Prop EX (doubly smooth existence is free) [PROVED from KNOWN].** Let M = Ψ(X, y) (all y-smooth n ≤ X are
units mod N). By Cauchy–Schwarz over the ≤ |H_S| occupied classes,
  #{(u, v) : u ≠ v ≤ X smooth, u ≡ v (mod N)} (ordered) ≥ M²/|H_S| − M.
For y = (log N)^κ, κ > 1, Ψ(X, (log X)^{κ′}) = X^{1−1/κ′+o(1)} (de Bruijn 1966; codex gives the binomial/Rankin
proof), so with X = N^{κ/(κ−1)+ε} one gets M ≥ 2|H_S| and at least M/2 = N^{1 + ε(1−1/κ) + o(1)} unordered relations.
A primitive pair (u′, v′) has at most Ψ(X/v′) ≤ Ψ(X/N) multiples (v′ > N), so there are at least
(M/2)/Ψ(X/N) = N^{1 − 1/κ + o(1)} distinct primitive ones. No GRH, no circle method. (This is nonexact
Prop AB with |H_S| in place of N.)
- *GRH asymptotic (codex eq. 4.3), checked.* With smooth weights and P⁻(N) > y, LS Prop 5.1 (GRH; uniform in the
  modulus q, literature critic) bounds the truncated Euler products, and codex's Perron/Mellin contour shift bounds the
  twisted smooth sums by H^{1/2+η}N^η for every nonprincipal χ mod N (an imprimitive χ has the same truncated Euler
  product as its primitive inducer, since P⁻(N) > y); character orthogonality then gives
  C_{N,w} = M_w²/φ(N) − D_w + O(H^{1+ε}N^ε), an asymptotic once H = N^A with A > κ/(κ − 2), κ > 2. (Our check:
  the error term is φ(N)⁻¹·Σ_χ |T_χ|² ≤ (H^{1/2+η}N^η)²; it is o(M²/φ(N)) iff 2A(1 − 1/κ) − 1 > A, i.e.
  A > κ/(κ−2).) This is the N-dependent analogue of LS 2012, and it is *easier*: no minor arcs, only characters.
- *Triply smooth.* The count of u′ + k′N = v′ with all three smooth is the LS ternary problem with one variable
  restricted to multiples of N. Heuristically ≈ N^{2K − 1 − (3K − 1)/κ} at height N^K (pairs N^{2K(1−1/κ)−1} times
  Pr[k smooth] = N^{−(K−1)/κ}), positive iff κ > (3K − 1)/(2K − 1) → 3/2 as K → ∞ — LS's threshold reappears
  [HEURISTIC]. A proof would need LS/Drappeau/Harper uniform in the auxiliary modulus N: open (codex §4b). It is
  not needed: factoring does not need k smooth.
- *Triples with a ≡ b (mod N) and a + b = c smooth* (codex): local model T*/N; open.

### 4.3 Where exactly it fails

1. **Not existence.** The N-analogue factoring needs (doubly smooth) is unconditional (Prop EX) and, under GRH,
   asymptotic (4.3). LS's actual achievement — controlling a + b = c with all three smooth, through minor arcs —
   concerns an N-free problem that factoring does not need.
2. **Not N-dependence.** The congruence enters through characters mod N, which are easier than LS's additive
   problem.
3. **Constructivity — this is the failure.** The circle method computes ∫S(α)²S(−α)dα; nothing in it locates a
   solution. Counting is not extraction (codex §4c): Prop EX establishes N^{1−o(1)} relations without identifying
   one. To exploit the count one must *find* collisions of the map u ↦ u mod N on smooth integers:
   - generic (birthday) search: T ≍ √(R|H_S|) = N^{1/2+o(1)} smooth numbers for R relations [model-PROVED in the
     independent-residue model; codex §4d];
   - Dixon-type search (test the residue for smoothness): N/Ψ(N, y) = N^{1/κ+o(1)} per relation, L[1/2, √2]
     after optimising y;
   - a constructive LS-with-congruence must beat MNFS L[1/3, 1.902], i.e. at polylog y beat N^{1/κ} by an
     exponential factor. By Prop SU, such a constructor *is* the doubly-split Π lister.
4. **Even without N no efficient construction is known.** For N-free polylog-smooth abc triples at prescribed large
   height we (and codex) know no polynomial-time algorithm. Sampling a, b and testing a + b costs H^{1/κ}; 3SUM-type
   bucketing costs H^{Ω(1)}. Polynomial identities f + g = h with split f, g, h of bounded degree d give smoothness
   y ≍ H^{1/d}, not polylog; unbounded d would need ideal Prouhet–Tarry–Escott-type split identities (known only for
   small sizes; systeq §6), and then the N-condition would require the identity to be non-exact modulo N at a small
   point, which systeq Lemma EX turns into an N-sized height (exact identities carry no information on N).
   Pascal-type and rational-closure constructions collapse after gcd normalisation or do not add rank (codex;
   closure lemma §3.3). Practical constructions are N-free: Lehmer's Pell method is exhaustive and provably complete for a fixed S
   but exponential (2^s − 1 Pell equations); Conrey–Holmstrom–McLaughlin's closure method is heuristic.
5. **The fruitfulness condition is an extra, independent requirement.** Even a constructive count would have to
   produce relations whose dependencies reach a nontrivial square root of unity (Lemma FR; odd index suffices);
   counting arguments say nothing about rank or index (codex §4c).

**Conclusion of §4.** LS is an existence theorem for an N-free ternary problem. Its N-analogue for factoring is a
binary congruence problem whose existence part is elementary, and whose constructive part is exactly the open
primitive Π (when the binary target, heights, output count and prefix cost match; constructing one triple, or a
triply-smooth subfamily, is a different guarantee). LS/Drappeau/Harper provide explicit major-arc structure and
counting estimates that could inform a search, but no extraction procedure follows from them; and the polylog
smoothness that looks so striking is already available, non-constructively, by pigeonhole.

---

## 5. Experiments (Task 5; pre-registered, commit 604ab31)

Grid: N = balanced_semiprime(n, seed), n ∈ {24, 28, 32, 36}, seeds 0–2; y = round((n ln 2)^κ), κ ∈ {1.25, 1.5, 2};
X = ⌊N^K⌋, K ∈ {9/8, 5/4, 3/2}; 17 cells, 51 runs, none capped. Peak RSS < 1.5 GB, one job at a time. Counts, ranks,
gcds and the q > 1 test (v > R as integers) are exact; model values and the q > 1.1 / ABCV predicates use floating
logarithms (no classification is near a threshold: max q = 0.83; codex). M1's numerical uncertainty (subsample,
left window near X) is ≈ 1–2% in the design cell (codex: 123.0 → 121.9 without subsampling, 120.8 with extended
windows), far below the observed deviations.

**Scorecard (all 51 runs; the same verdicts with the 4 design runs excluded).**

| P | Prediction | Result |
|---|---|---|
| P1 | pairs/M1 ∈ [0.6, 1.5] for M1 ≥ 100; median ∈ [0.8, 1.2] | **FAILED** in 1 of 51 runs: (36, 1.25, 9/8, seed 1) R1 = 0.557 (M1 = 283, observed 158). Median 0.960; range 0.557–1.073 |
| P1′ | equidistribution model M0 overestimates (R0 < 0.85) for K = 9/8 | held (24/24; R0 = 0.12–0.77) |
| P1b | quadratic characters mod N negligible (\|R1b − R1\| ≤ 0.05) | held (max 0.012) |
| P2 | \|H\| = φ(N) | held (51/51) |
| P3 | Σc_r² ≥ Ψ²/\|H\| (identity) | held |
| P4 | rank ≤ π(y) − 2; = π(y) − 2 with ≥ R relations | held (36/36 with ≥ R relations; no run above the bound) |
| P5 | at predicted rank: pooled fruitful ∈ [0.4, 0.6]; z-dim 2 and ≥ 1 fruitful in every run | **FAILED**: pooled 0.462 held, but 2 of 49 runs had z-dim < 2 and 0 fruitful: (24, 1.25, 9/8, 2) z-dim 1; (36, 1.25, 9/8, 0) z-dim 0 |
| P6 | triply-smooth count / M1 restricted to smooth k ∈ [0.6, 1.5] (runs with X/N > y) | held (21/21; 0.81–1.16) |
| P7 | abc quality > 1.1 in ≤ 10⁻³ of primitive relations; 0 Prop ABCV violations | held: **0** of 5,861,386 primitive relations have q > 1 at all; max q = 0.830 |
| P8 | triply ≤ Evertse bound; triply ≥ π(y) | held (min triply/π(y) = 1.12) |
| P9 | N-free primitive S-unit triples saturate at the VKM totals (1, 4, 17, 63, 190) by H = 2⁴⁰ | held (largest c: 2, 9, 128, 4375, 18225) |

**Reading.**
- E1: the N-congruence count of doubly smooth pairs is described by the N-analogue of the LS singular series
  (M1) within ±10% for κ ≥ 1.5 or K = 3/2; at κ = 1.25 M1 overestimates by 10–45% (R1 = 0.56–0.82 in 11 of the 12
  runs at K = 9/8, median 0.71; 0.81–0.99 at K = 5/4); 20 of 51 runs are outside ±10%, all low. Without local
  factors the error reaches ×8 at X/N small. The multiplicative quadratic correction 1 + Σ(T_χ/Ψ)² is negligible (P1b), but codex
  points out that this is not the exact quadratic-character term of the orthogonality identity, Σ_χ quad T_χ²/(2φ(N)),
  which reaches 9.6% of the observed pairs in the sparse cell (28, 1.25, 9/8, 0) and is offset there by the other
  characters (the finite-size effect). So "no role" is not established; the quadratic characters are simply not
  separately modelled. We read the
  κ = 1.25 deficit as model imperfection at very small y (one saddle point α(X, y) ≈ 0.35 for numbers of very
  different sizes; Haar-uniform unit residues assumed), not as N-structure: it is the same across seeds and n.
- E2: rank exactly π(y) − 2, the quadratic-character obstruction of Lemma FR(5), in every run with enough
  relations. Fruitfulness 0.46 pooled (the sampler draws uniformly from D∖{0}, so with full square-root image the null expectation is
2^{d−1}/(2^d − 1), d = dim D, slightly above 1/2; e.g. 8/15 for d = 4).
  Including the design run, three runs had 0/40 fruitful; a fourth had z-dim 1 but 22/40 fruitful (§6 PH2).
- E3: triply-smooth relations exist in quantity even at toy size; abc is irrelevant in practice (no primitive
  relation even reaches quality 1).
- E4: the N-free S-unit equation is saturated by tiny solutions; it has no N-content.

---

## 6. Post-hoc analysis (disclosed; not pre-registered; `src/abcsunit/posthoc.py`, `abcsunit_posthoc.jsonl`)

**PH1 (P1 failure).** Pairs cluster: each primitive pair (u′, v′) carries its Ψ(X/v′) multiples, so the count is a
sum over few primitive events with multiplicities. Under a Poisson model for primitive pairs, sd = √Σm² (the observed Σm² is a plug-in estimate; independence and
calibration are assumed, so z is a plausibility measure, not a calibrated test — codex). The failing
run has 18 primitive pairs and z = (158 − 283)/48.4 = −2.6. It is not an isolated fluctuation: the κ = 1.25,
K = 9/8 family as a whole has R1 median 0.71 and a combined (Stouffer) z ≈ −5.6 (proof critic). Across all runs
z ranges −4.4 … +3.7; M1 overestimates at small κ and small K and slightly underestimates at K = 3/2. The registered band was set on
pair counts rather than primitive counts — a design flaw; the failure stands as registered.

**PH2 (P5 failure).** The Z-rank of the E2 relation set equals s in every run with z-dim 2, and is < s in exactly
the runs with z-dim < 2: (24, 1.25, 9/8, 0) [design; also short of R] Z-rank 8/11, z-dim 0; (24, 1.25, 9/8, 2)
10/11, z-dim 1; (36, 1.25, 9/8, 0) 14/16, z-dim 0; (36, 1.25, 9/8, 2) 14/16, z-dim 1 (that run still had 22/40
fruitful, since its one nontrivial square root of unity was ±ε). This is exactly Lemma FR(3): z-dim = ρ − r ≤
rank_Z − r. The registered condition (r = s − 2) was the wrong one; the right one is ρ = s (odd index). All four
runs had fewer than R = 3π(y) + 30 relations (17–32). The Z-rank deficiency comes from multiplicative coincidences
among few relations; it disappears once the relation count exceeds ≈ 2s.

**Scoring notes (proof critic).** The pre-registration writes M0 with prefactor 1/|H|; the code (correctly, see
model.py) uses N/|H| on the density sum — the registered *code* was scored. P3 is the Cauchy–Schwarz inequality,
not an identity. The "triply ≤ Evertse" half of P8 was not scored in the first version of score.py (codex); it is now scored
against both the registered looser bound 3·7^{2(s+2)+3} (the raw field `evertse_bound_log10`) and the bound
3·7^{2s+3} proved in §2.3; every run passes both.

**Design-stage note.** The first E2 design (all chain edges, including non-primitive scaled copies (gu, gv)) gave
0/15 fruitful dependencies for the same reason: scaled copies have identical e-vectors, so their dependencies are
Z-relations. This led to Lemma FR before the registered run (prereg, design disclosure).

---

## 7. Reviews

- **Codex blind round (`codex/abcsunit_1.out`).** Independently: the same translation (incl. the empty "c = kN
  smooth" reading and the fixed-coefficient unit equation), abc vacuity for unrestricted k, Evertse never bites,
  Π-equivalence, the near-full-rank-is-not-fruitfulness phenomenon, the closure lemma, the GRH collision
  asymptotic (4.3), "eight is an artefact", birthday benchmarks; probability ≈ 1%. We adopted (4.3), the closure
  lemma, the fixed-coefficient formulation and the Drappeau reference, with attribution. Codex did not give Prop
  ABCV's explicit exception set, Prop DW, Lemma SC, or the odd-index criterion; it did not mention Harper 2016.
- **Codex adversarial round (`codex/abcsunit_2.out`): MAJOR REVISION.** All applied:
  - verdict overclaim ("cannot") → "the cited bounds establish no general barrier", with the powerful-family
    caveat and the exact example 2¹⁹ + 7153 = 3¹² (#3, the "single most important correction");
  - NFS row rewritten (S_K-unit membership is exact; Norm(m − θ) = f(m) is a multiple of N; fixed-coefficient
    form) (#2); sign coordinate (#1);
  - primitive vs all-pairs count in Prop EV(a) (#4; also proof critic #5);
  - Cor ABC-T reading and "never bites" qualified (#5); Lemma SC's optimality claim deleted (#6, counterexample
    ord_N(3) < ord_N(2) for N = 101·103);
  - Lemma FR: odd index sufficient, not necessary; necessity only given r = s − 2 (#7, example N = 77; also proof
    critic #2–3); sampler null 2^{d−1}/(2^d − 1) (#8); one-sided relations do not escape the image requirement (#9);
  - Prop DW limited to the complete kernel/determinant; ℓ-adic normalisation; modular level is not a universal
    lower bound (#10); SU region made triangular (#11); LS primitive error term, Harper confirmed (#12);
    constructivity wording (#13);
  - experiments: M1 accuracy restated (#14; also proof critic #1); quadratic-character "no role" withdrawn (#15);
    M1 numerical uncertainty reported (#16); PH1 z is a plug-in plausibility measure (#17); P8 upper half now
    scored (#18); floating-point predicates disclosed (#19).
- **Opus proof critic: ACCEPT WITH FIXES.** Re-derived every proof; confirmed Lemma FR(3) on 300 random relation
  sets for N = 1147 (0 mismatches) and reproduced the scorecard. Findings applied: M1 accuracy misreported
  (MAJOR; 20 of 51 runs outside ±10%, κ = 1.25 family median −29%, ≈ 5.6σ); false Z-rank sentence in Lemma FR
  (counterexample N = 1147, now a regression test); "must have odd index" → sufficient; H ⊇ (Z/N)^*[2] hypothesis
  and relabelling of §3.1; primitive count; ABCV reading (Dixon/QS covered by the height ≤ C·N^{1+ε} form; Schnorr
  relations get Cor ABC-T); Σ1/p term in the table; closure lemma scope; Lemma SC label; run counts in §0; NITs
  (xyz inequality typo, Evertse place count in K, sign coordinate, finite odd index, 2s discrete logs, Miller
  criterion, prereg M0 prefactor, P3 is an inequality, sampler bias).
- **Opus literature critic (v2): OK WITH FIXES**, no major attribution errors (report relayed by the
  orchestrator). It confirmed the key black boxes: Evertse with arbitrary coefficients; LS theorem numbers, ranges
  and Prop 5.1 uniformity in q (with our imprimitive-character extension); Harper's unconditional range (also
  checked by us against arXiv 1408.1662); Bach's constant under P⁻(N) > y. Applied: ES-T lower bound
  exp(Ω(√(s/log s))) and KS (#1); ES-T 1988 conjecture ≪ exp(s^{2/3+ε}) for every S, still ≫ s for S ∪ {p, q} (#2);
  VKM uses Faltings/modularity height bounds in terms of N_S (#3, strengthens Prop DW 5); VKM Thm A covers n ≤ 16,
  counts modulo S₃ (#4); LS 2011 Thm 1.1/2.1 = LS 2012 Prop 1.1 (#5); Perron/Mellin step attributed to codex, H
  notation clash resolved (#7); Bach < 2 log² m, Lamzouri–Li–Soundararajan < 3 log² m, Miller 1976 and Bach 1984 added
  (#13); Lehmer exhaustive but exponential (#15); Schnorr 2021 / Ducas / Adleman 1995 (#17); prior art for Lemma FR
  and softened §0 wording (#18); dangling §7.1 reference removed (#20). With these, the earlier UNVERIFIED labels
  are cleared.
- **Bibliographic completions.** Stewart–Tijdeman, Monatsh. Math. 102 (1986) 251–257; Stewart–Yu, Math. Ann. 291
  (1991) 225–230; Stewart–Yu, Duke Math. J. 108 (2001) 169–181; Drappeau, Canad. J. Math. (2015) (as cited by the literature critic) and
  "Sommes friables d'exponentielles et applications" (Harper's major-arc input); Canfield–Erdős–Pomerance, J. Number Theory 17 (1983) 1–28; Smart's S-unit papers (mid-1990s) and "The algorithmic resolution of Diophantine
  equations", LMS Student Texts 41 (1998); Beukers–Schlickewei 1996 (subgroups of (K^*)²); de Weger, J. Number
  Theory 26 (1987) 325–367; Erdős–Stewart–Tijdeman, Compositio Math. 66 (1988) 37–56.

---

## 8. Reproduction

```
memory_pressure | tail -1
.venv/bin/python -m pytest -q tests/test_abcsunit.py
.venv/bin/python -m src.abcsunit.run_experiment E123   # experiments/raw/abcsunit_E123.jsonl (≈ 6 min)
.venv/bin/python -m src.abcsunit.run_experiment E4     # experiments/raw/abcsunit_E4.jsonl
.venv/bin/python -m src.abcsunit.score                 # P1–P9
.venv/bin/python -m src.abcsunit.posthoc               # PH1, PH2 (post hoc)
```
