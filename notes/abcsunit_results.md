# abcsunit: factoring relations as S-unit equations — abc, Evertse, Lagarias–Soundararajan

Agent `abcsunit`, phase 19, 2026-10-21. DRAFT v1 (before codex round 2 and the Opus critics; §7 records the reviews).

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
- **abc and S-unit theory give no barrier for factoring, and they cannot.** The reason is one line: every
  relation contains N in its radical, and N is squarefree, so abc's inequality is automatically satisfied unless
  the relation is unusually "powerful" (Prop ABCV). Counting bounds (Evertse) apply only to the triply-smooth
  subfamily and allow exponentially more relations (3·7^{2s+3}) than the ≈ s a method needs; the relations factoring
  actually uses (k free) are the points of a full-rank lattice and are infinite in number (Prop EV).
- **"Solve S-unit equations with an N-congruence" is not a new target.** Listing smooth pairs u ≡ v (mod N) is,
  verbatim, the doubly-split product-listing problem Π of pilist / nonexact, i.e. listing short vectors of Λ_N
  (Prop SU). It is SUFFICIENT for heuristic polynomial-time factoring with a polylog factor base (smoothenum T8),
  with one precise extra condition we isolate: the found relations must generate Λ_N up to odd index
  (Lemma FR) — near-full parity rank is not enough.
- **Known S-unit algorithms have no N-dependent step that helps.** de Weger's p-adic reduction linearises the
  *higher-level* conditions ℓ^m | x − y by ℓ-adic logarithms; the level-one condition N | x − y that relations
  need is exactly the lattice Λ_N, and computing it factors N (Prop DW). The real-place reduction with target log N
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
    N-version would be exactly Π-doubly-split, i.e. the open primitive itself. LS thus re-describes the target; it
    is not a route to it.
- **Experiments (pre-registered; 2 of 11 registered predictions failed, both explained post hoc).** The
  N-congruence count follows the LS-type singular-series model M1 to within ≈ ±10% (median ratio 0.96 over 51 runs;
  one low-count run at 0.56, z = −2.6 under the clustering noise); local factors are essential (the
  equidistribution model is off by up to ×8); F₂-rank is exactly π(y) − 2 (quadratic-character obstruction) in
  every run with enough relations; fruitfulness is 1/2 except in two runs where the relations did not generate
  Λ_N up to odd index — exactly Lemma FR. None of 5.9·10⁶ primitive relations has abc quality > 1 (max 0.83).

**Main results.**

| Result | Status |
|---|---|
| Translation: doubly smooth = Λ_N points; triply smooth = fixed-coefficient S-unit equation; one-sided (Dixon/QS) = abc-type but not S-unit; NFS = mixed S-unit congruence in K × Q (§1) | PROVED (definitions) |
| Prop ABCV: abc can only exclude relations with kN or u N^ε-powerful; fraction ≪ N^{−ε/2} among multipliers | PROVED (conditional only in reading abc) |
| Cor ABC-T: under abc, primitive triply-smooth relations have height ≤ C_ε(N·R_S)^{1+ε} (N-relative LS Prop 1.1) | PROVED assuming abc |
| Prop EV: #primitive triply-smooth relations ≤ 3·7^{2s+3} (any height, any N); never below the ≈ s needed; doubly smooth: infinitely many | PROVED from KNOWN (Evertse 1984) |
| Lemma SC: no set of s primes has more S-units ≤ X than the first s primes | PROVED |
| Prop SU: S-unit congruence listing = Π doubly split = short vectors of Λ_N | PROVED (definitional) |
| Lemma FR: fruitful dependencies ⟺ image of (Λ₀ ∩ 2Z^s)/2 in H[2]; full when [Λ_N : Λ₀] is odd; z-dim = ρ − r | PROVED |
| Prop DW: det Λ_N (S ⊇ primes < 2 log² N) = φ(N) under GRH, so the level-one step of de Weger's method factors N | PROVED from KNOWN (Bach 1990) |
| Prop EX: unconditional abundance of doubly (log N)^κ-smooth relations at height N^{κ/(κ−1)+ε} | PROVED from KNOWN (Ψ(x,(log x)^κ) = x^{1−1/κ+o(1)}) |
| Codex eq. (4.3): GRH asymptotic for collisions at height N^A, A > κ/(κ−2) | PROVED by codex from LS Prop 5.1; checked by us (§4.2) |
| Closure lemma (codex): N-dependent rational closure of seed pairs cannot raise rank | PROVED |
| LS constructivity: no polynomial-time construction known; constructive N-version ≡ Π | assessment / KNOWN (absence in LS) |

**Assessment (subjective).** Probability that the abc/S-unit route yields L[α < 1/3]: ≤ 1% (ours; codex
blind: ≈ 1% over twenty years). The route contributes a clean reformulation and a no-go for "Diophantine
barriers", not a mechanism.

---

## 1. The exact translation (Task 1)

### 1.1 Doubly smooth relations are lattice points [PROVED]

Let u, v be positive S-units, u ≠ v, u ≡ v (mod N); write v − u = kN (k ≠ 0, WLOG v > u, k ≥ 1). Put d = gcd(u, v).
Since d | kN and gcd(d, N) = 1 (P⁻(N) > y), d | k. With u′ = u/d, v′ = v/d, k′ = k/d:

  u′ + k′N = v′,  gcd(u′, v′) = 1,  and gcd(k′, u′v′) = 1

(a prime dividing k′ and u′ divides v′). The map (u, v) ↦ e(v) − e(u) depends only on (u′, v′) and is a
bijection between primitive pairs and nonzero vectors of Λ_N up to sign, with u′ = ∏ℓ^{e⁻}, v′ = ∏ℓ^{e⁺} and
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
| NFS (two-sided) | φ(a − bθ) = a − bm in Z[θ]/(N, θ − m) ≅ Z/N | **not a rational S-unit relation.** It is a *mixed congruence*: α = a − bθ is an S_K-unit of O_K only up to class group and units (which NFS sidesteps with ideals and quadratic characters); α − (a − bm) = −b(θ − m), and (θ − m) has norm N. With b smooth and S_K′ ∋ the degree-one primes above N it becomes a genuine S-unit equation in K of degree d (Evertse bound 3·7^{d+2s′}) |

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

**Reading.** abc constrains only relations in which kN or u is N^ε-powerful. Dixon, QS, NFS (rational side) and
Schnorr-type relations have generic multipliers; abc implies nothing about their sizes or numbers. In particular
abc cannot give a statement "no relation family with candidates of log-size < f(N)": the only lower bound on
relation heights is the elementary one, max(u, v) > N for distinct positive congruent integers (codex; the
"aggregate bound" of systeq §4), and abc gives only *upper* height bounds.

**Comparison with the smoothness-size principle.** That principle (candidates of log-size ≍ x^{2α}(log x)^{1−2α}
for L[α]) is a statement about the probability that a random integer is smooth (Canfield–Erdős–Pomerance), plus the
assumption that a method samples candidates "randomly". It has no Diophantine content. abc neither implies nor
strengthens it; a product lister would by design violate its sampling assumption (codex §2a).

### 2.2 The N-relative Lagarias–Soundararajan bound [PROVED assuming abc]

**Cor ABC-T.** Assume abc. For every ε > 0 there is C_ε such that every primitive triply-y-smooth relation
(u′, k′, v′ all y-smooth) satisfies v′ ≤ C_ε(N·R_S)^{1+ε}, i.e. log v′ ≤ (1+ε)(log N + (1+o(1))y) + O_ε(1).
*Proof.* rad(u′v′k′N) ≤ R_S·N. ∎ (codex eq. 2.2.) This is LS 2011 Prop./Thm 1.1 (abc ⇒ κ₀ ≥ 1) with N added
to the radical.

**Reading.** For y = (log N)^κ, κ > 1, the bound allows heights up to exp((log N)^κ), far beyond any polynomial
search space; it binds only for y = O(log N), where triply-smooth relations of height N^A, A > (1+ε)(1 + y/log N),
are excluded. That is an *upper* height bound (it confines the triply-smooth solution set, which is helpful, not
obstructive) and concerns a factor base too small for any method anyway.

**Effective versions [KNOWN].** Stewart–Yu (Duke Math. J. 108 (2001) 169–181, Thm 1) give
log c ≪ R^{1/3}(log R)³ with R = rad(abc) (the earlier Stewart–Tijdeman 1986 and Stewart–Yu 1991 bounds are
weaker). With R ≥ N this allows
c up to exp(N^{1/3+o(1)}): vacuous. Baker-type lower bounds for |u − v|, u ≠ v S-units, carry constants exponential
in s and are weaker than the trivial |u − v| ≥ N.

### 2.3 Counting bounds do not bite [PROVED from KNOWN]

**Prop EV.**
(a) *Doubly smooth.* Λ_N is a full-rank lattice, so there are infinitely many primitive doubly smooth relations;
no S-unit counting theorem applies. Under a height bound X their number is Θ(Ψ_S(X)²/|H|) in the models of §5
and ≥ ½(Ψ²/|H| − Ψ) unconditionally (Prop EX).
(b) *Triply smooth.* The map (u′, k′, v′) ↦ (X, Y) = (v′/k′, −u′/k′) is injective into the solutions of
N⁻¹X + N⁻¹Y = 1 in Γ_S (u′/k′ is in lowest terms since gcd(k′, u′) = 1). Evertse (Invent. Math. 75 (1984)
561–584, Thm 1: ax + by = 1 over a number field of degree d has ≤ 3·7^{d+2t} solutions in S-units, t = number of
places in S including the infinite ones) gives, for Q with s finite primes (t = s + 1):

  #{primitive triply-smooth relations} ≤ 3·7^{2s+3},  for every N and every height.

(Beukers–Schlickewei, Acta Arith. 78 (1996) 189–199: 2^{8(r+1)} for rank r subgroups — weaker here.)
(c) *Does it bite?* A method needs about s + O(1) relations (columns). 3·7^{2s+3} > s for every s ≥ 1, so the
bound never prevents a factor base of size s from supplying enough relations, in any regime of s (L[β], polylog,
O(log N), o(log N)). It provides no existence either. The conjectured truth for the first s primes,
exp(s^{2/3+o(1)}) (Erdős–Stewart–Tijdeman's conjectural exponent, as recalled by LS 2011 §1.5), also exceeds s (the precise form of the ES-T conjecture is as recalled by LS; we did not check ES-T's
original wording).
(d) *Measured.* In every registered run the triply-smooth primitive count exceeded π(y) (min ratio 1.12; P8).

**Lower bounds [KNOWN; relevance nil].** Erdős–Stewart–Tijdeman (Compositio 66 (1988) 37–56):
exp((4 − ε)(s/log s)^{1/2}) solutions for suitable S; Konyagin–Soundararajan (J. Number Theory 124 (2007)
193–199): exp(s^{2−√2−ε}); LS 2011 Thm 1.4 under GRH: exp(s^{1/8−ε}) for the first s primes. These concern chosen
or initial prime sets and N-free equations. They defeat a hypothetical "too few unit solutions" barrier but do not
produce N-relations.

**Lemma SC [PROVED].** For any set of s primes ℓ₁ < ⋯ < ℓ_s and the first s primes p₁ < ⋯ < p_s,
Ψ_{ℓ}(X) ≤ Ψ_{p}(X) for all X: ∏ℓ_i^{e_i} ↦ ∏p_i^{e_i} is injective and does not increase values. So, under the
equidistribution heuristic (pairs ≈ Ψ_S²/(2|H_S|)), no clever choice of S (Erdős–Stewart–Tijdeman style) beats the
first s primes for N-relations, up to computable index-2^{O(1)} effects such as restricting to ℓ with a fixed
Jacobi symbol (ℓ/N). The ES-T/KS gains come from N-free arithmetic coincidences, which carry no information on N.

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
the form bN + v on the box [1, X/N] × [1, X] with base = primes ≤ y, v split; the bijection is b = (u − v)/N.
Equivalently: list R short vectors of Λ_N in the norm of §1.1. Requiring k smooth (triply) restricts the output set
to the fixed-coefficient unit equation; requiring c = kN smooth empties it.

**Sufficiency [model-PROVED; = smoothenum T8, restated].** Suppose SU can be solved with polynomial delay for
y = (log N)^κ, X = N^K, with K > 1 and 2K(1 − 1/κ) > 1 + δ (so that, heuristically, N^{δ+o(1)} ≥ π(y) relations
exist), and its outputs satisfy the generation condition of Lemma FR. Then N factors in polynomial time:
s = π(y) = (log N)^{κ+o(1)} columns, R = s + O(1) relations, dense F₂ elimination poly(s), and a dependency gives
Z² ≡ 1 with Z uniform in H[2] (Lemma FR), splitting with probability 1/2 (codex §3a). Total time
T_SU + poly(log N). With y = L[β] the same reduction preserves the exponent β. Listing *all* solutions is
exponential (the counts are N^{Θ(1)}), so only a useful prefix can be asked for (smoothenum T8, codex §4d).

### 3.2 Lemma FR: fruitfulness needs odd index, not just parity rank [PROVED]

Let Λ = Λ_N ⊂ Z^s, H = Z^s/Λ, Λ′ = {x : 2x ∈ Λ}, so Λ′/Λ ≅ H[2]. Let e₁, …, e_R ∈ Λ be relation vectors,
Λ₀ their Z-span, r the F₂-rank of their reductions in Z^s/2Z^s, ρ the F₂-rank of their classes in Λ/2Λ, and
D = {c ∈ F₂^R : Σc_i e_i ∈ 2Z^s} the dependency space. Define z : D → Λ′/Λ, c ↦ [½Σc_i e_i] (lift c to {0,1}).
The square root of unity of a dependency is Z = ∏ℓ^{−½Σc_i e_i} (mod N) (§core.square_from_dependency; Z = W/V).

1. z is a well-defined homomorphism (changing lifts changes Σc_i e_i by 2Λ₀ ⊂ 2Λ).
2. im z = ((Λ₀ ∩ 2Z^s)/2 + Λ)/Λ. (Write y ∈ Λ₀ ∩ 2Z^s as Σa_i e_i; c = a mod 2 lies in D and
   Σc_i e_i ≡ y mod 2Λ₀.)
3. **dim im z = ρ − r.** The kernel of Λ/2Λ → Z^s/2Z^s is (Λ ∩ 2Z^s)/2Λ = 2Λ′/2Λ ≅ H[2]; the image V of Λ₀ in Λ/2Λ
   has dim ρ and maps onto a space of dim r, and V ∩ kernel corresponds exactly to im z.
4. If [Λ : Λ₀] is odd, then ρ = s and im z = H[2]. (For x ∈ Λ′ and m = [Λ : Λ₀] odd, 2mx ∈ Λ₀ ∩ 2Z^s and
   [mx] = [x] in H[2].)
5. r ≤ s − dim(H/2H) (rows lie in (Λ + 2Z^s)/2Z^s, of codimension dim H/2H). For N = pq and H = (Z/N)^* this is
   **r ≤ s − 2**: the parity rows are orthogonal to c_p, c_q (indicators of ℓ with (ℓ/p) = −1, resp. q), because
   (uv/p) = (u/p)² = 1.

So with H[2] = {±1, ±ε} ≅ F₂² (N = pq, H = (Z/N)^*), a uniformly random dependency is fruitful with
probability exactly 1/2 if ρ = r + 2 (z uniform on H[2]); if ρ = r + 1 the probability is 0 or 1/2 according as
im z = {±1} or {1, ±ε}; if ρ = r it is 0. ρ = r + 2 holds when [Λ_N : Λ₀] is odd, and **fails whenever the
relations have Z-rank < s** (ρ ≤ rank_Z Λ₀).
Reaching the maximal parity rank r = s − 2 is not sufficient. Codex found the extreme case independently: integer
combinations of parity-independent kernel vectors have only unfruitful dependencies (codex §3a). For one-sided
(QS/Dixon) relations each relation brings its own free generator t, which is why the issue is invisible there.
*Prior statement:* this is the standard structure of index-calculus square roots; the explicit odd-index criterion
for purely S-unit relations: no prior statement found by us. **Measured:** ρ − r = 2 ⟺ rank_Z = s in all 51 runs
(§6, PH2).

### 3.3 Known S-unit algorithms and the N-condition

**The algorithms [KNOWN].** de Weger (CWI Tract 65, 1989; also J. Number Theory 26 (1987)); Tzanakis–de Weger;
Smart (1995, 1998); the Sage solver of Alvarado–Koutsianas–Malmskog–Rasmussen–Vincent–West (2019/2020); von
Känel–Matschke (arXiv 1605.06079; Memoirs AMS 2023). Pattern: (1) Baker/Matveev and Yu height bounds, exponential
in s; (2) LLL reduction at the real place (log ℓ_i) and at each ℓ ∈ S (ℓ-adic logarithms); (3) sieving/enumeration
of the remaining box. None has complexity polynomial in s; even the output is exp(Ω(√s)) for some S (ES-T) and
grows like 2^s empirically for the first s primes (VKM 2016 Thm A: 1, 4, 17, 63, 190, 545, 1433, …, 839638 for
n = 1, …, 15; our E4 reproduces n ≤ 5).

**Prop DW (where the N-condition sits) [PROVED from KNOWN].**
1. For ℓ ∈ S the ℓ-adic step encodes ℓ^m | X − Y, m large, through e ↦ Σe_i log_ℓ(ℓ_i) mod ℓ^m after
   killing the level-one torsion by the power ℓ − 1 = |F_ℓ^*| (known, small).
2. A factoring relation needs only the **level-one** condition N | v − u, i.e. e ∈ Λ_N. There is nothing to
   linearise: the level-one lattice is the target itself.
3. Computing a basis of Λ_N — even only det Λ_N = |H| — factors N when S ⊇ {primes < 2 log² N}: under GRH these
   primes generate (Z/N)^* (Bach, Math. Comp. 55 (1990), the least prime outside a proper subgroup is < 2 log² m),
   so |H| = φ(N) = N − p − q + 1. Conversely, computing Λ_N from p, q is two discrete-log problems (mod p, mod q).
   So the step de Weger handles trivially for small ℓ is, for the composite N, at least as hard as factoring.
   (Unconditionally, |H| is a multiple of the order of every element of H, which suffices for Miller-type
   splitting when H is large; we do not need this.)
4. The real-place step with an inhomogeneous target log N is Schnorr's lattice (EUROCRYPT 1991): find e with
   Σe_iℓ log ℓ_i ≈ log N, producing u − vN = w small, a triply-smooth relation if w is smooth (§1.2). That is the
   remit of agent `schnorrcvp`; we do not duplicate it.
5. Modularity methods (Frey curve y² = x(x − u′)(x + k′N)) need modular forms of level ∝ rad(u′v′k′N) ≥ N: cost
   N^{Ω(1)}.
6. Codex adds: N is a known ℓ-adic unit for ℓ ∈ S, so ℓ-adic logarithms of N are free; substituting composite N
   for a prime place is invalid; whether *every* composite-modulus logarithmic method is factoring-hard is open.

**Closure lemma (codex §4c) [PROVED].** From pairs (a, a + N), (b, b + N), c = a(b + N)/(b − a) satisfies
c + N = b(a + N)/(b − a); when integral both are smooth, but e(c) − e(c + N) = [e(a) − e(a + N)] − [e(b) − e(b + N)].
Rational closure of seeds never enlarges Λ₀, so it cannot raise ρ (or r), and by Lemma FR it cannot create
fruitfulness either.

---

## 4. Lagarias–Soundararajan and smooth abc solutions (Task 4)

### 4.1 What LS prove [KNOWN; checked against the arXiv text]

- *JTNB* 23 (2011) 209–234 ("the xyz conjecture"): primitive A + B + C = 0 with all prime factors
  ≤ (log H)^κ, H = max(|A|, |B|, |C|); κ₀ := the threshold exponent.
  - abc ⇒ κ₀ ≥ 1 (Thm 1.1; more generally weak abc with exponent κ₁ ⇒ κ₀ ≥ κ₁);
  - unconditionally (Stewart–Yu) finitely many with smoothness ≤ (3 − ε) log log H (Thm 2.2);
  - GRH ⇒ κ₀ ≤ 8 (Thm 1.2); strong xyz conjecture κ₀ = 3/2 (heuristic Ψ³/H ≥ 1 ⇔ 3(1 − 1/κ) > 2);
  - Thm 1.4 (GRH): ≫ exp(s^{1/8−ε}) primitive S-unit solutions for the first s primes.
- *Proc. LMS* 104 (2012) 770–798 ("Counting smooth solutions to A + B = C"), Thms 2.1–2.3 under GRH, for
  (log H)^{8+δ} ≤ y ≤ exp((log H)^{1/2−δ}): weighted counts
  N(H, y; Φ) = 𝔖_∞(1 − 1/κ, Φ)·𝔖_f(1 − 1/κ, y)·Ψ³/H·(1 + O(log log y/log y)), the primitive version with 𝔖_f*,
  and the primitive density 1/ζ(2 − 3/κ). 𝔖_f(c, y) = Π_{p≤y}(1 + (p−1)/(p(p^{3c−1} − 1))·((p − p^c)/(p−1))³)
  Π_{p>y}(1 − 1/(p−1)²). Method: circle method over exponential sums of smooth numbers, Hildebrand–Tenenbaum
  saddle point; GRH bounds partial Euler products (their Prop 5.1) on the minor arcs. The exponent 8 is an artefact
  of the error term H^{3/4+ε}Ψ against Ψ³/H (codex: 2 − 3/κ > 7/4 − 1/κ ⇔ κ > 8).
- Later: Drappeau (Math. Proc. Cambridge 154 (2013), GRH, unweighted asymptotics; codex); Harper (Compositio 152
  (2016) 1121–1158) with Drappeau's major arcs: the asymptotic **unconditionally** for (log H)^C ≤ y ≤ H, C large
  (abstract; to be checked by the literature critic).
- LS give no algorithm. Their explicit examples come from de Weger's exhaustive lists, Gross–Zagier singular
  moduli, Poonen–Schaefer–Stoll and the Monster.

### 4.2 Making it N-dependent

**Prop EX (doubly smooth existence is free) [PROVED from KNOWN].** Let M = Ψ(X, y) (all y-smooth n ≤ X are
units mod N). By Cauchy–Schwarz over the ≤ |H| occupied classes,
  #{(u, v) : u ≠ v ≤ X smooth, u ≡ v (mod N)} (ordered) ≥ M²/|H| − M.
For y = (log N)^κ, κ > 1, Ψ(X, (log X)^{κ′}) = X^{1−1/κ′+o(1)} (de Bruijn 1966; codex gives the binomial/Rankin
proof), so with X = N^{κ/(κ−1)+ε} one gets M ≥ 2|H| and at least M/2 = N^{1 + ε(1−1/κ) + o(1)} unordered relations.
A primitive pair (u′, v′) has at most Ψ(X/v′) ≤ Ψ(X/N) multiples (v′ > N), so there are at least
(M/2)/Ψ(X/N) = N^{1 − 1/κ + o(1)} distinct primitive ones. No GRH, no circle method. (This is nonexact
Prop AB with |H| in place of N.)
- *GRH asymptotic (codex eq. 4.3), checked.* With smooth weights and P⁻(N) > y, LS Prop 5.1 (GRH) bounds the
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
   - generic (birthday) search: T ≍ √(R|H|) = N^{1/2+o(1)} smooth numbers for R relations [model-PROVED in the
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
   closure lemma §3.3). Practical constructions (Lehmer's Pell method; Conrey–Holmstrom–McLaughlin) are N-free and
   empirical.
5. **The fruitfulness condition is an extra, independent requirement.** Even a constructive count would have to
   produce relations generating Λ_N up to odd index (Lemma FR); counting arguments say nothing about rank or index
   (codex §4c).

**Conclusion of §4.** LS is an existence theorem for an N-free ternary problem. Its N-analogue for factoring is a
binary congruence problem whose existence part is elementary, and whose constructive part is exactly the open
primitive Π. Nothing in LS's method is constructive, and the polylog smoothness that looks so striking is already
available, non-constructively, by pigeonhole.

---

## 5. Experiments (Task 5; pre-registered, commit 604ab31)

Grid: N = balanced_semiprime(n, seed), n ∈ {24, 28, 32, 36}, seeds 0–2; y = round((n ln 2)^κ), κ ∈ {1.25, 1.5, 2};
X = ⌊N^K⌋, K ∈ {9/8, 5/4, 3/2}; 17 cells, 51 runs, none capped. Peak RSS < 1.5 GB, one job at a time. Every number
below is exact except model values (floating point).

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
  (M1) to within ≈ ±10%; without local factors the error reaches ×8 at X/N small. Quadratic characters mod N play
  no role at these sizes. M1 has a systematic bias (≈ −5 to −10% at K = 9/8, +3 to +7% at K = 3/2; §6 PH1), which
  is model imperfection, not structure.
- E2: rank exactly π(y) − 2, the quadratic-character obstruction of Lemma FR(5), in every run with enough
  relations. Fruitfulness 0.46 pooled.
- E3: triply-smooth relations exist in quantity even at toy size; abc is irrelevant in practice (no primitive
  relation even reaches quality 1).
- E4: the N-free S-unit equation is saturated by tiny solutions; it has no N-content.

---

## 6. Post-hoc analysis (disclosed; not pre-registered; `src/abcsunit/posthoc.py`, `abcsunit_posthoc.jsonl`)

**PH1 (P1 failure).** Pairs cluster: each primitive pair (u′, v′) carries its Ψ(X/v′) multiples, so the count is a
sum over few primitive events with multiplicities. Under a Poisson model for primitive pairs, sd = √Σm². The failing
run has 18 primitive pairs and z = (158 − 283)/48.4 = −2.6. Across all runs z ranges −4.4 … +3.7, with sign tied to
K (negative at 9/8, positive at 3/2): M1 is biased by up to ≈ 10% at these sizes. The registered band was set on
pair counts rather than primitive counts — a design flaw; the failure stands as registered.

**PH2 (P5 failure).** The Z-rank of the E2 relation set equals s in every run with z-dim 2, and is < s in exactly
the runs with z-dim < 2: (24, 1.25, 9/8, 0) [design; also short of R] Z-rank 8/11, z-dim 0; (24, 1.25, 9/8, 2)
10/11, z-dim 1; (36, 1.25, 9/8, 0) 14/16, z-dim 0; (36, 1.25, 9/8, 2) 14/16, z-dim 1 (that run still had 22/40
fruitful, since its one nontrivial square root of unity was ±ε). This is exactly Lemma FR(3): z-dim = ρ − r ≤
rank_Z − r. The registered condition (r = s − 2) was the wrong one; the right one is ρ = s (odd index). All four
runs had fewer than R = 3π(y) + 30 relations (17–32). The Z-rank deficiency comes from multiplicative coincidences
among few relations; it disappears once the relation count exceeds ≈ 2s.

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
- Codex adversarial round, Opus proof critic, Opus literature critic: pending (§7 to be completed).

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
