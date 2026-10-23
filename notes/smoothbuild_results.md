# smoothbuild: constructing smooth numbers instead of searching for them

Agent `smoothbuild`, phase 24, 2026-10-23. Version 2 (after the codex blind round and the Opus proof critic; §9–§10 record the reviews).

**Files.**
- Pre-registration: `notes/smoothbuild_prereg.md` (commit `85d70fe`, before any run).
- Code: `src/smoothbuild/`.
  - `core.py`: operational helpers (factor-base trial division, Gauss split, the hidden-root reduction of Prop EQ1)
    and exact checks (Burnside ranks, roots mod p, Lemma AT counts, GF(2) kernels).
  - `e1_recurrence.py`, `e2_class.py`, `e3_cgpt.py`: the three pre-registered experiments (harness).
  - `posthoc_e2.py`: post-hoc diagnostics for E2's two failed predictions (labelled post hoc).
  - `analyze.py`: the prediction table.
- Tests: `tests/test_smoothbuild.py` (10 tests).
- Data: `experiments/raw/smoothbuild_E{1,2,3}.jsonl`, `smoothbuild_E2_posthoc.jsonl`; summary
  `experiments/analysis/smoothbuild_summary.json`.
- GPT-6.1-sol: blind round `codex/smoothbuild_1.{md,out}`; adversarial round `codex/smoothbuild_2.{md,out}`.

**Labels.** KNOWN (published; cited) · PROVED (proof here; elementary unless said) · HEURISTIC (named assumption) ·
MODEL (exact consequence of a named model) · MEASURED · POST HOC · OPEN · assessment. "No prior statement found by us"
is the strongest novelty claim made.

**Notation.** x = ln N; L[α, c] = exp((c + o(1))x^α(ln x)^{1−α}); y = L[1/3, κ] the smoothness bound; h = L[2/3, ν]
the window of Prop SI (suffmap §1.4). A *free* side of a relation is one whose factorisation the construction
supplies; a *tested* side is one whose smoothness must be checked. Theorem SS (suffmap §1.2): a sampling method whose
tested sides have total log-size ≍ x^σ costs L[σ/2]; in particular a tested side of log-size x (a full residue mod N)
means L[1/2].

---

## 0. Verdict

**Plain language.**
- **No construction found here beats sampling at L-scale.** For each of the five angles we can say exactly where the
  construction stops: either the tested side of every relation stays the size that Theorem SS charges for, or the
  construction step is itself as hard as the original problem.
- **"Known factorisation" is not a weaker currency than smoothness.** A congruence of squares needs only atoms that
  recur, not small primes. That is true, but:
  - an atom that occurs in only one relation can never be in a dependency (Lemma KI). A large known prime helps only
    if it recurs, and each non-recurring one costs at most one dimension of the kernel;
  - making a given large prime P recur is the same problem as writing P's residue class as smooth × square
    (Lemma LP), i.e. the residue-class construction problem itself, up to random self-reduction;
  - on a side that is uniformly random (the residue, under H_S), no set of K multiplicatively independent
    known-factored atoms beats the K smallest primes (Lemma AT, PROVED; the consequence is MODEL). If atoms may be
    chosen after seeing the residues (all primes allowed),
    Croot–Granville–Pemantle–Tetali 2012 prove that the gain over the best single smoothness bound is at most a
    constant factor. E3 measures it: allowing every prime reaches the first square dependency after 0.49–0.51·J₀
    integers, against 1.30–1.37 times as many for the J₀-optimal bound y* fixed in advance;
  - structured families make large primes recur more often, but by a bounded factor. For values of an irreducible
    degree-d polynomial the pair-recurrence rate is (r_G − 1) times that of random integers, where r_G is the rank of
    the Galois group on the roots. So the gain is ≤ d − 1, reached exactly by regular groups such as cyclotomic Φ_n
    (Prop GR, PROVED for the Chebotarev average). E1 confirms the law for eight families (observed/predicted
    0.94–1.01; the two lowest are explained by ± root pairs, a toy-range effect).
- **Pre-factored families cannot be steered.** Cyclotomic and Aurifeuillian values, Lucas/Chebyshev values,
  binomials, products of consecutive integers and norms × units either:
  - have no tie to N at all (their factorisations are exact identities over Z, which carry no information mod N);
  - or are tied to N only through a residue side that is a full-size residue (Dixon, L[1/2]);
  - or must land in a short window near kN, which sparse families cannot do (sinonobliv Lemma PL / Prop FF).
  Steering kN onto a Cunningham-type value needs a root of a fixed polynomial mod N. For x² − c with c random that is
  Rabin's problem (as hard as factoring); for a fixed c no equivalence is known; in general it is polynomial selection
  (NFS). Among one-parameter families the square gives the smallest residue (QS).
- **Sum-product gives existence, not construction.** Bourgain–Konyagin, Bourgain–Glibichuk–Konyagin, Garaev and
  Shkredov bound how products of structured sets *distribute*. That is evidence *for* the flatness under which
  meet-in-the-middle does not amortise; E2 measures smooth × smooth products in a window mod N at 0.99–1.01 of the
  flat prediction. The one fast constructive primitive in this area is the 2-dimensional lattice (continued
  fractions). It produces *ratio* splits t ≡ u/v with |u|, |v| ≈ √N in polynomial time, and these dominate *product*
  splits xy ≡ t, which have no known fast search (Prop RP). A ratio split leaves two √N-size tested pieces: total
  tested log-size x, i.e. CFRAC/Dixon level L[1/2].
- **Residue-class construction is index calculus (Theorem EQ).**
  - An average-case constructor of smooth representations u ≡ t·v for random squares t (roots hidden from it)
    factors N, rigorously (EQ1, Dixon's argument).
  - With unbounded exponents the representation problem is randomized-equivalent to "factoring plus discrete logs
    modulo each prime factor", given that the small primes generate (Z/N)^× (EQ2; the DL direction uses an
    exchangeability argument from the Opus proof critic).
  - A modulus-generic constructor also does prime-field index calculus, descent included (EQ3).
  - The *weakest* named smooth construction that still factors is a hidden-root relation sampler, and that is
    randomized-equivalent to factoring itself (EQ4). The targets form a partial order: HRS ≡ FACT, below
    CLASS_∞ ≡ FACT ∧ DL, below bounded CLASS. A construction is useful only if it is *cheaper* than sampling, not
    merely possible.
- **Descent does not transfer, and the missing ingredient is not size-halving.** In Z/N and in F_q[X]/(I) alike,
  every non-identity relation has one side at least as large as the modulus (Lemma DSC); F_q escapes that floor only
  because the big side is split for free. What BGJT and
  Kleinjung–Wesolowski have is a family of identities in which that big side splits *by construction* into atoms
  from a fixed polynomial-size set (the linear polynomials), which therefore recur in every relation. Z has no such
  family: additive Frobenius is impossible (systeq Prop F), and our Lemmas KI/LP/AT show that the integer substitute,
  known-factored atoms, must recur, which no fixed family achieves. In the one-sided Z/N model every node of a
  descent tree is a full-size instance of the root problem.
- **Strongest remaining construction target:** unchanged. It is Prop SI (a y-smooth integer in [kN, kN + h] at
  L[1/3, τ], τ < 0.757), or its Bohr-set form for exponent < 1/3. Replacing "y-smooth" by "known factorisation"
  in Prop SI gives a target worth no more than Dixon's L[1/2, √2] unless the outputs happen to be y-smooth
  (Prop KS, HEURISTIC).
- **Assessment.** Probability that any route in this note beats L[1/3, 1.902]: ≤ 1% (ours).

**Main results.**

| Result | Statement (short) | Status |
|---|---|---|
| Lemma KI | an atom occurring in one relation lies in no dependency; each non-recurring large atom costs ≤ 1 kernel dimension | PROVED (pure linear algebra; folklore: singleton removal) |
| Lemma LP | a second relation containing a large atom P ⟺ a representation of P's class as (square) × (smooth); = CLASS up to random self-reduction | PROVED (identity) |
| Lemma AT | K multiplicatively independent atoms compose ≤ Ψ(X, p_K) integers ≤ X | PROVED, elementary; no prior statement found by us |
| CGPT | uniform random integers: first square product at T ∈ [(π/4)e^{−γ}, e^{−γ}]·J₀, any primes allowed | KNOWN (Annals 2012) |
| Prop GR | irreducible F of degree d: Chebotarev average of ρ_F(p)(ρ_F(p) − 1) = r_G − 1 ∈ [1, d − 1]; = d − 1 iff G regular | PROVED (Chebotarev + Burnside); pair-count law HEURISTIC, E1 0.94–1.01 |
| Prop KS | Prop SI with "known factorisation" in place of "y-smooth" needs ≳ J₀(N) = L[1/2, √2] outputs unless the outputs are y-smooth | HEURISTIC (H_S + CGPT-type counting on the y-rough parts) |
| Prop RP | in every box, ratio splits (lattice, poly-time enumeration) are at least as numerous as product splits (modular hyperbola) and cheaper to find | PROVED (counting + enumeration); no fast product search known (assessment) |
| Theorem EQ | EQ1 hidden-root class constructor ⇒ factoring (also found by codex); EQ2 CLASS_∞ ≡_R FACT ∧ DL(p_i) (Bach 1984 in factor-base form); EQ3 modulus-generic ⇒ F_p index calculus; EQ4 hidden-root sampler ≡_R FACT | PROVED (EQ1, EQ4); EQ2 KNOWN in substance; EQ2/EQ3 given generation of (Z/N)^× by small primes (e.g. GRH), rank step by exchangeability |
| Lemma DSC | u ≡ t·v (mod N), (u, v) ∦ (t, 1) ⇒ \|u\| + \|t\|·\|v\| ≥ N; same with degrees in F_q[X]/(I) | PROVED (trivial) |

**Experiments** (pre-registered; §8).
- **E1 (Prop GR):** all 8 families and the random control inside their bands. Observed/predicted 0.94–1.01.
  Cyclotomic Φ₅/S4 = 3.07 (law: 3), Φ₇/S6 = 5.09 and Φ₉/S6 = 5.03 (law: 5), D4/S4 = 1.91 (law: 2), A4/S4 = 0.99 (law: 1).
- **E2 (constructed vs tested side):** 3 of 5 predictions held at every size.
  - P2b failed at 40 bits (0.725, band edge 0.75): Gauss-split pieces were 0.73–0.83 as smooth as independent
    same-size integers at every size. Post hoc, a *coprime* same-size null explains the systematic deficit completely
    (1.00–1.01): the two pieces of a split are coprime. The 40-bit failure itself is noise at the band edge (0.752 on
    the post-hoc sample).
  - P2d failed at every size (0.78–0.79). Post hoc: 36–41% of multiplier-steered smooth hits repeat a relation
    already found for the same target (17% of steered splits are exact rescalings of the unsteered split). A Poisson
    model on the number of *distinct* relations per target reproduces the per-target success to within 3%.
  - Flatness (P2e): smooth × smooth products in a window mod N, O/E 0.99–1.01.
- **E3 (CGPT):** T/Ĵ₀ median 0.49, 0.51, 0.51 at 32/40/48 bits (P3a held).
  - P3b failed at 32 bits: one seed's 97th integer was a perfect square. The threshold was miscalibrated: perfect
    squares alone give ≈ 18% a-priori failure probability across 12 seeds at 32 bits.
  - P3c held, but trivially: T_best ≡ T by construction, a design flaw disclosed in §8.
  - Post hoc, the fixed-in-advance bound y* needs 1.30–1.37 × T integers.

---

## 1. Setting and the degrees-of-freedom ledger

**What a relation is.** A congruence A² ≡ w (mod N), with A a known residue and w a nonzero rational with known
factorisation over a set of pairwise coprime *atoms* (primes, or any integers whose divisibility we can track). This
merged rational form covers Dixon, CFRAC, QS/MPQS/SIQS and integer relations u ≡ v. It does not cover NFS-type
relations (prime-ideal columns, square root in a number ring); Lemma KI is pure linear algebra and holds for any
column set, but the free/tested ledger below is stated for the merged rational form. An
F₂-dependency among the exponent vectors of the w_i gives X² ≡ Y² (mod N). A relation u ≡ v (mod N) between two
integers is the case w = v/u, A = 1, or w = v when u = A² is a square.

**The ledger.** For each construction we record:
- the *free* degrees of freedom (parameters chosen, whose contribution is factored for free);
- the *tested* side and its log-size, which Theorem SS prices;
- the *tie to N*: where the congruence mod N enters.

Two facts constrain every row of the ledger.
- **Barrier 1** (systeq §4, codex; PROVED): if U ≡ V (mod N) and U ≠ V, then max(|U|, |V|) ≥ N/2. A relation therefore
  always has a side of log-size ≥ x − ln 2. Either that side is free (a square, or smooth by construction) or it is
  tested at full size.
- **Theorem SS** (suffmap): a tested side of log-size x gives L[1/2]; reaching L[1/3] needs tested log-size ≍ x^{2/3}.

So a construction can help only by making the *tested* side small while keeping the *free* side's atoms recurrent.
§2–§6 check each angle against this. Previously excluded routes are not redone:
- siconstruct: Lemma SR (scale invariance), Lemma RC (placed factor, tested cofactor), Prop CV (promise CVP);
- sihard: Theorem J (oblivious certificate constructors), Prop MM (meet-in-the-middle), Prop DLI;
- sinonobliv: Lemma PL, Prop FF (sparse families cannot hit windows), Prop NV (norms with units = integer knapsack),
  Prop LC;
- crtdecode: CRT list decoding;
- dsr: Theorem SB, FN, BOX-∞, Lemma PT, Prop EQ-dom;
- systeq: Prop F, F′, Lemma CF, DES, EX;
- ptebarrier/ptesplit: Lemma R, Theorem W;
- schnorrcvp: the log-lattice route.

---

## 2. Angle 1: pre-factored algebraic families

### 2.1 Catalogue

| Family | What is "known" | Tie to N | Ledger | Verdict |
|---|---|---|---|---|
| Cyclotomic Φ_n(a), a^n − 1 = ∏_{d\|n}Φ_d(a) | pieces Φ_d(a) (not their primes); primes ≡ 1 (mod d) or \| d | none, except N \| Φ_n(m) (special N, systeq M13: no gain) | free side = pieces; residue a^n mod N is full-size | no tie; Dixon if residue tested |
| Aurifeuillian Φ_n(x) = C² − nxD² at x = n·m² | two half-size pieces | none (special N: constant only, snfsrig) | splitting a candidate into k = O(1) known pieces changes ρ by k^u = L^{o(1)} (suffmap §1.1) | constant at most |
| Lucas/Lehmer U_n, Chebyshev T_n − T_m | cyclotomic factorisation over a quadratic / real cyclotomic field; primitive divisors ≡ ±1 (mod n) | mod N only as group orders (p ± 1 methods) | not a relation family | group-order family (M), not relations |
| Binomials C(n, k), consecutive products | complete factorisation (Kummer/Legendre): n-smooth | none | C(p, 1) = p, so binomial products are exactly the smooth numbers; large-k binomials are sparse | nothing new (systeq M3); sparse ⇒ Prop FF |
| Norms of ε∏π_i^{e_i} in fields with large unit groups | norm = smooth by construction | through φ: Z[θ] → Z/N | norms are an integer knapsack (Prop NV); φ(γ) of a long product is a full residue mod N; small φ(γ) needs small coordinates = NFS | Dixon, or NFS (Prop LC) |
| Cunningham-type kN = F(a) | factorisation of F(a) over Z[x] | needs a root of F mod N | for F = x² − c, square roots mod N (Rabin's problem when c is random); in general polynomial selection | NFS |

**Why exact identities carry no information** (systeq Lemma EX, codex Barrier 1). An identity among integers below N
that holds over Z is a trivial congruence. Information about N enters only through a residue, which is a
full-size integer unless the construction places the free side within h of a multiple of N.

### 2.2 The known-factorisation relaxation

**Lemma KI (kernel inclusion) [PROVED; folklore: singleton removal in every large-prime variation].**
- Let relations have exponent matrix M = [M_s | M_ℓ] over F₂, where M_ℓ holds the columns of a set 𝒜_ℓ of "large"
  atoms. Then the left kernel of M equals ker(M_s) ∩ ker(M_ℓ) ⊆ ker(M_s).
- (a) A relation whose row has an odd entry in a column of 𝒜_ℓ in which no other row is odd lies in no dependency.
- (b) For a fixed set of rows, each large atom that occurs (oddly) in the rows lowers the kernel dimension by at most
  one relative to ker(M_s).
- *Proof.* Immediate from the definitions: a left-kernel vector of M annihilates both blocks; (a) is the column
  with one odd entry; (b) adding one column lowers the rank-deficiency by at most one. ∎
- *Reading.* A large atom contributes to a dependency only if it recurs; each one that does not costs at most one
  relation. Rows *with* large atoms are rows one would not otherwise have (E3: allowing every prime reaches the first
  dependency ≈ 25% sooner), so known large primes can help, but only through recurrence, and Lemma LP and CGPT bound
  how much. KI does not say that a square free side is optimal; §2.3's preference for squares rests on residue size.
- *Merged form.* A free-side atom P can also be matched by a tested-side P in another relation; that second relation
  is again one containing P, which is what Lemma LP describes.

**Lemma LP (recurrence = class representation) [PROVED; an identity].**
- Let P be an atom, P > y, gcd(P, N) = 1. A relation A² ≡ P^{±1}·w (mod N) with w a y-smooth rational is the same
  thing as a representation of the class of P^{∓1} as (square) × (y-smooth rational). Hence a *second* occurrence of
  P (the only way P can join a dependency, by KI(a)) is a solution of the residue-class problem of §4 at target P,
  up to random self-reduction: CLASS applied to P·q·s² (s a random unit, q a factor-base prime making Pq a square
  class) yields such a relation, and any such relation is a relation containing P. (§4's CLASS is guaranteed only on
  uniformly random squares, so the equivalence is via that reduction, not literal.)
- *Chance recurrence.* If large atoms behave like independent primes with Pr[P] ≍ c/P (HEURISTIC, the large-prime
  model), R relations give ≈ R²Σ_{P∈(y,P_max]}c²/P² ≍ R²/(y ln y) coincidences, dominated by P near y. For atoms of
  size ≥ L[α] with α > 1/3 (for example a prime near kN, or a cofactor of Φ_n(a)) chance recurrence needs
  R ≳ L[α]^{1/2} relations: hopeless. So the known-prime relaxation is the classical large-prime variation (a
  constant, or L^{o(1)}) for atoms near y, and the residue-class problem itself for larger atoms.

**Lemma AT (atom optimality) [PROVED, elementary; no prior statement found by us].**
- *Statement.* Let a₁ < … < a_K be multiplicatively independent integers > 1 (e.g. pairwise coprime), and p_K the
  K-th prime. Then the number of integers ≤ X that are products of the a_i (with repetition) is at most Ψ(X, p_K).
- *Proof.*
  - The prime-exponent vectors of a₁, …, a_i have rank i, so together they involve at least i distinct primes, all
    ≤ a_i. Hence a_i ≥ p_i. (For pairwise coprime atoms: their least prime factors are distinct.)
  - By independence, distinct exponent vectors e give distinct products.
  - The map e ↦ ∏p_i^{e_i} is injective and ∏p_i^{e_i} ≤ ∏a_i^{e_i} ≤ X. ∎
- *Dependent atoms.* Linear algebra works over a coprime base (Bernstein's refinement), and the bound applies with K
  the number of columns. The Opus critic's brute force (K = 2, 3, atoms in [2, 40), X = 10⁴, 10⁵) found no
  violation of Ψ(X, p_K) even for dependent atom sets with K = the number of atoms.
- *Consequence (MODEL: H_S on the residue).* If the tested side is uniform in [1, V], then for any K atoms fixed in
  advance Pr[tested side composed of them] ≤ Ψ(V, p_K)/V. Known-factored atoms (cyclotomic pieces, Cunningham
  primes, large primes) never beat the K smallest primes as columns for a random side.
- Checked exhaustively on random coprime atom sets in `test_lemma_AT_atom_optimality`.

**Adaptive atoms: CGPT [KNOWN].**
- Croot, Granville, Pemantle, Tetali, "On sharp transitions in making squares", Ann. of Math. 175 (2012) 1507–1550.
  For a_i uniform in [1, x], let T be the first index at which some subsequence product is a square. *Any* primes may
  be used. Then T ∈ [(π/4)(e^{−γ} − o(1)), (e^{−γ} + o(1))]·J₀(x) with probability 1 − o(1), where
  J₀(x) = min_y π(y)x/Ψ(x, y) is Pomerance's count for the best single smoothness bound. They conjecture
  T ∼ e^{−γ}J₀.
- So for random-like integers, allowing every prime with its known factorisation beats a fixed smoothness bound by at
  most a constant factor. E3 (§8) measures T/Ĵ₀ ≈ 0.50, and 1.30–1.37 × T for the strategy that fixes y in advance.

**Prop GR (Galois recurrence) [Chebotarev average PROVED; pair-count law HEURISTIC; E1].**
- *Statement.* Let F ∈ Z[x] be irreducible of degree d ≥ 2, G its Galois group acting on the d roots, and
  r_G = number of G-orbits on ordered pairs of roots (the rank). For p ∤ disc(F)·lc(F) let ρ_F(p) = #roots mod p.
  - (i) The natural-density average of ρ_F(p)(ρ_F(p) − 1) over primes equals r_G − 1.
  - (ii) 1 ≤ r_G − 1 ≤ d − 1. The lower bound is attained iff G is 2-transitive, the upper bound iff G acts
    regularly (|G| = d; e.g. Φ_n, every abelian G).
  - (iii) (HEURISTIC: roots equidistributed in short intervals mod p.) For a ∈ an interval [a₀, a₀ + A) with
    A ≤ P₁, the number of pairs a ≠ a′ with F(a), F(a′) sharing a prime p ∈ (P₁, P₂] is ≈ (r_G − 1)·C(A, 2)·Σ_p p⁻².
    For independent uniform integers it is ≈ C(A, 2)·Σ_p p⁻². The heuristic needs a₀ ≫ P₂, or no pairs of roots
    related by an affine map with small coefficients: if r′ = −r, then a′ ≡ −a forces p | a + a′ < 2(a₀ + A), so such
    pairs never fall in the window once p > 2(a₀ + A) (Opus critic; this is E1's deficit for D4 and V4, §8.1).
- *Proof of (i)–(ii).*
  - For unramified p, ρ_F(p) = fix(Frob_p). By Chebotarev the density average of the class function fix(g)² − fix(g)
    is (1/|G|)Σ_g (fix(g)² − fix(g)).
  - Burnside on X × X gives (1/|G|)Σ fix(g)² = r_G, and (1/|G|)Σ fix(g) = 1 by transitivity.
  - r_G − 1 counts orbits on off-diagonal pairs. By transitivity these correspond to the orbits of a point stabiliser
    G_x on X ∖ {x}. There are at least 1 and at most d − 1 of them, with d − 1 iff G_x fixes every point, i.e.
    G_x = 1 (faithful action), i.e. |G| = d. ∎
- *Reading.*
  - Structured pre-factored families make large primes recur up to d − 1 times more often than random integers.
    Cyclotomic families attain the maximum.
  - For Φ_n the two effects are explicit: an *eligible* prime (p ≡ 1 mod n) has φ(n) roots and a pair rate
    φ(n)(φ(n) − 1) times that of a random prime, but only a 1/φ(n) share of primes is eligible. The net factor is
    φ(n) − 1.
  - For families of degree d ≤ (ln N)^{O(1)} this is an L^{o(1)} factor.
  - The same families are *less* smooth at small primes (only p ≡ 1 mod n divide Φ_n(a) beyond p | n): systeq
    measured R = 0.61 and 0.25.
  - Families with unbounded recurrence, such as {a^n − 1 : n ≤ n₀} sharing Φ_d(a), are exact identities on the free
    side. Lemma KI says their best case is to behave like a square.
- E1 (§8) confirms (iii) for S4, A4, D4, C4, V4, S6 and two C6 families: observed/predicted 0.94–1.01.

### 2.3 Steering

- **Near kN.** sinonobliv's Lemma PL / Prop FF (PROVED, averaged): a family fixed independently of the window must
  contain a p/((h + 1)ρ) fraction of all smooth integers near X. Every family above with O(X^{1−ε}) members near X
  fails by X^{−ε}: cyclotomic, Aurifeuillian, binomial with k ≥ 2, Lucas, and one-parameter polynomial values.
- **Residue floor of one-parameter families [MODEL; KNOWN in substance].** If F has degree d ≥ 2 and the free side is
  F(x) near kN, the tested residue F(x) − kN is ≈ |F′(x₀)|·θ ≈ d(kN)^{1−1/d}·θ, with θ the distance of the real root
  x₀ of F(x) = kN to the nearest integer. Over k ≤ K the best residues are ≈ d(kN)^{1−1/d}/K. d = 2 with F = x² is
  the quadratic sieve, and the residue size d(kN)^{1−1/d} is smallest at d = 2. So among one-parameter pre-factored
  families the square gives the smallest residue. Multi-parameter families with small residues are
  polynomial selection (resbarrier, heurnfs).
- **Into a residue class.** "Steer a product of family members to ≡ □ (mod N)" with the square's root known is the
  residue-class problem of §4 restricted to the family's semigroup. Lemma AT says the family's atoms are no better than
  primes.
- **Cunningham-type multiples.** kN = F(a) for a fixed F needs a root of F mod N. For F = x² − c this is a square
  root of c mod N. For c ranging over random squares that is Rabin's problem, equivalent to factoring. For a single
  fixed c no equivalence is known to us, but no algorithm either. For general F it is the polynomial-selection step of
  NFS, whose height/norm floors are resbarrier's.

### 2.4 Verdict on angle 1

Refuted as a route to faster factoring, in the following sense. Every family either has no tie to N (exact
identities), or ties to N through a full-size tested residue (Theorem SS: L[1/2]), or must hit windows (Prop FF:
exponential for sparse families). "Known factorisation" adds at most a constant (CGPT) or a factor d − 1 (Prop GR) in
recurrence, and for large atoms it *is* the residue-class problem (Lemma LP).

**Prop KS (known-factorisation Prop SI) [HEURISTIC].**
- *Setting.* A constructor outputs, for each k, some u_k ∈ [kN, kN + h] together with its full factorisation, at
  cost L^{o(1)}. It does not promise smoothness. The residue v_k = u_k − kN ≤ h is tested for y-smoothness.
- *Claim.* Assume H_S for the u_k: their prime factors above y are distributed like those of uniform integers of size
  kN. Then every dependency needs the y-rough parts of the u_k in it to multiply to a square, because primes > y
  cannot occur in a y-smooth v. By CGPT-type counting at size N, that needs ≳ J₀(N)^{1−o(1)} = L[1/2, √2]
  outputs with smooth residues.
- *Reading.* "Known factorisation" in place of "y-smooth" turns Prop SI's 1.2114 into Dixon's L[1/2, √2]. The
  relaxed oracle is worth what Prop SI is worth only on its y-smooth outputs, i.e. when it is a Prop SI constructor.
  A constructor that *controls* the recurrence of large primes is solving Lemma LP's class problem.
- *Status.* HEURISTIC. CGPT prove the counting for uniform integers. Applying it to the y-rough parts of
  window-restricted integers is the H_S step.

---

## 3. Angle 2: multiplicative closure and sum-product

**What the theorems say.**
- Bourgain–Konyagin, Bourgain–Glibichuk–Konyagin: exponential sums over multiplicative subgroups of size ≥ p^δ have
  power savings. Small subgroups are additively equidistributed.
- Garaev, Garaev–Karatsuba, Shparlinski: the congruence xy ≡ λ (mod p) with x, y in short intervals has the expected
  number of solutions, for every λ once the box is large enough (Weil/Kloosterman range), for almost all λ in smaller
  boxes.
- Shkredov, Erdős–Szemerédi type: lower bounds on |A·A| + |A + A|, upper bounds on multiplicative energy.

All of these are *non-concentration* statements: products of structured sets spread out. They are existence and
equidistribution results with no algorithm.

**How they relate to the barriers already proved.**
- Flatness, i.e. Pr[a·b mod N ∈ W] ≈ |W|/N for a, b from structured sets, is exactly the hypothesis H_flat under which
  meet-in-the-middle does not amortise (siconstruct Lemma MA). sihard's Prop MM proves an averaged-over-N form
  unconditionally.
- So sum-product theorems are *evidence for* the barrier, not tools against it. Character-sum proofs of flatness
  for a fixed window need |A||B| ≥ N³/h² (Pólya–Vinogradov error √N log N per character against main term
  |A||B|h/N), which is useless at L-scale. The averaged form is trivial (each product lies in ≤ h + 1 windows), and
  that is what Prop FF and Prop MM use.
- E2-M4 (§8) measures smooth × smooth products landing in t·[1, W] mod N: O/E 0.99–1.01 at 40–64 bits.
- **Closure adds nothing.** The relations form a lattice (exponent vectors of products ≡ 1). Multiplying or dividing
  relations adds lattice vectors and creates no new F₂-information. A product of two Prop SI outputs is a relation
  with residue v₁v₂ ≤ h², not a new short one.

**Prop RP (ratio splits dominate product splits) [PROVED for the counting and enumeration; the hardness of product
search is an assessment].**
- *Ratio splits.* For t ∈ (Z/N)^×, Λ_t = {(u, v) ∈ Z² : u ≡ tv (mod N)} is a lattice of determinant N. Extended Euclid
  stopped at the first remainder below √N finds (u, v) ∈ Λ_t with |u| < √N and |v| ≤ N/r_{i−1} ≤ √N in O(x²) bit
  operations (`gauss_split`). This is the classical rational-reconstruction step of index calculus.
- *Product splits.* H_t = {(x, y) : xy ≡ t (mod N)} is the modular hyperbola, which is not a lattice. Points with
  x, y ≤ H exist for H ≥ N^{1/2+ε} for almost all t, and for all t once H ≥ N^{3/4+ε} (Weil/Kloosterman; for
  squarefree composite N via Estermann's bound). We know no search faster than enumerating x ≤ H and testing
  tx⁻¹ mod N ≤ H, which costs ≈ N/H per point.
- *Domination, box by box.* In the box [−H, H]², Λ_t has 4H²/N + O(H/λ₁ + 1) points, each enumerable in polynomial
  time from a reduced basis, while H_t has ≈ H²/N points found at cost ≈ N/H each. A product split xy ≡ t with x, y
  y-smooth gives the relation t ≡ xy; a ratio split u ≡ tv with u, v y-smooth gives t ≡ u/v. So for every box the
  ratio splits are at least as many, carry relations of the same shape, and are found in polynomial time per point.
  (Multiples j·b₁ of one lattice vector give the same relation, which is E2's duplicate effect.) ∎
- *Reading.* The sum-product world (products) is a strictly worse construction than the lattice world (ratios), and
  the lattice world gives two tested pieces of total log-size x (Theorem SS: L[1/2], the constant of Dixon). Going
  below needs more pieces, i.e. a higher-dimensional lattice: Schnorr's route, refuted for this model (schnorrcvp).

**What is new relative to sihard/siconstruct.**
- They treated integer products in intervals near kN (the hyperbola over Z) and showed that meet-in-the-middle is
  exponential.
- Here, residue-class products mod N (the modular hyperbola) are dominated by ratio splits. Those are polynomial-time
  but leave a full-size tested ledger. Placing one free factor by the 2-D lattice is the same operation as steering
  by a free multiplier (E2-M2), and it changes the target, not the size of the tested pieces (P2c: 0.98–1.01).
- E2 also shows a correlation: free smooth multipliers often reproduce the *same* relation. In 17% of steered splits
  the split of t·s is the rescaled split (s·u, v), which is the identical relation. Only 59–64% of steered hits are
  distinct (post hoc).

**Verdict on angle 2:** no constructive content at L-scale. The one constructive primitive (2-D lattice) is Lemma RC /
rational reconstruction, which is already priced.

---

## 4. Angle 3: residue-class construction and index calculus

**Definitions.** Fix y and A > 0. A *class constructor* CLASS_A(y; T, δ) takes (N, t) with t ∈ (Z/N)^×. In time T it
outputs exponent vectors (e, f) over {−1} ∪ {primes ≤ y} with ∏p^{e_p} ≡ t·∏p^{f_p} (mod N) and
Σ(|e_p| + |f_p|) ln p ≤ A·x, succeeding with probability ≥ δ when t is a uniformly random square. CLASS_∞ drops the
size bound (exponents of poly(x) bits). A *hidden-root sampler* HRS(y; T) outputs pairs (X_i, w_i), w_i y-smooth,
with X_i² ≡ w_i (mod N), such that conditioned on all w_j and on X_j (j ≠ i), X_i is uniform among the square roots
of w_i. Each pair costs expected time T.

**Theorem EQ.** Let N be odd with ≥ 2 distinct prime factors, none ≤ y, and K = π(y) + 1.
- **EQ1 [PROVED; Dixon's argument].** CLASS_A(y; T, δ) ⇒ factoring N in expected time
  O((K + 1)(T + x²)/δ + K^ω x + K·A·x³) with success probability ≥ 1/2 (ω the matrix-multiplication exponent; the
  last term computes Y from exponents of size A·x; w(N) below is the number of distinct prime factors).
  - *Proof.* Draw a uniform unit r, put t = r² mod N, and call the constructor on (N, t). Each success gives
    r² ≡ ∏p^{e_p − f_p}. The constructor sees only t, so r is uniform among t's square roots given its view. After
    K + 1 successes an F₂-dependency gives X = ∏r_i and Y = ∏p^{(Σ(e−f))/2} with X² ≡ Y². Multiplying one r_i by a
    nontrivial square root of 1 does not change the view, so Pr[X ≡ ±Y] ≤ 2^{1−w(N)} ≤ 1/2.
    `hidden_root_class_factor` implements it with the Gauss split as the constructor: 32/40/48-bit N factor
    STRUCTURALLY in tests. ∎
- **EQ2 [KNOWN in substance: Bach 1984 (DL mod N ≡ FACT ∧ DL mod primes), here in factor-base form; PROVED given
  that the primes ≤ y generate (Z/N)^×, e.g. under GRH by Bach's bound when y ≥ 2 ln² N; the rank step of (⇒ DL) by
  the exchangeability argument below (Opus proof critic), or codex's primitive-kernel-vector argument (§9.1)].**
  CLASS_∞ ≡_R FACT ∧ DL(p) for every p | N.
  - (⇐) Factor N and compute discrete logs of the factor-base primes and of t modulo each p^k ∥ N. A solution e of
    Σe_p log p ≡ log t modulo each group exponent exists by generation and is found by linear algebra over Z.
  - (⇒) Factoring is EQ1. For DL of h modulo p | N (work mod p^k ∥ N if p² | N), let g generate (Z/p)^× and
    replace h by g·h if h is a non-square. Make m queries, each t_i = g^{2k_i} mod p (k_i uniform) lifted by CRT with
    a uniform square mod N/p, so t_i is a uniform square mod N. At one uniformly random position put the target
    h·g^{2k} instead, which has the same distribution. Each success gives a vector v_i = e_i − f_i ∈ Z^K with
    v_i·(log ℓ)_ℓ ≡ 2k_i (mod p − 1), plus log h at the target position.
  - *Rank step (exchangeability).* Over Z/(p − 1), the indices i with v_i outside the span of the other vectors
    give a strict chain of submodules of (Z/(p − 1))^K, so there are at most K·log₂ p of them. The oracle's view does
    not reveal the target's position, so the target is one of them with probability ≤ K·log₂ p/s, s the number of
    successes. Otherwise v_T = Σc_i v_i (coefficients from a Howell/Hermite form, without factoring p − 1), and
    log h ≡ Σc_i·2k_i − 2k (mod p − 1). Cost O(K·x·(T + x²)/(εδ)) for error ε. ∎
  - This replaces the rank heuristic of v1; the appeal to Pomerance 1987 is withdrawn (an adversarial oracle's
    relations are not uniform).
  - The bounded version CLASS_A asks in addition for a *short* vector of the relation lattice in a coset. That is a
    Schnorr-type CVP and is not implied by FACT ∧ DL (OPEN).
- **EQ3 [PROVED; padding as in dsr Lemma PT, rank step as in EQ2].** A constructor that works for all admissible N, or only
  for N = p·r with r a known prime of size p^{o(1)} above y, gives CLASS modulo every large prime p at the same cost.
  So it performs prime-field index calculus with descent at cost K·T/δ. At y = (ln N)^{O(1)} and polynomial T that is
  polynomial-time factoring *and* polynomial-time prime-field DL. At L[1/3, τ] it would undercut the 1.902 of
  prime-field DL (Matyukhin 2003) once τ is small (sihard Prop DLI).
- **EQ4 [PROVED].** HRS(y; T) ⇒ factoring in O(K(T + x²) + K^ω x), as in EQ1. Conversely, given the factorisation,
  an HRS with T = poly(x) exists. Draw w = w₁w₂ with w₁, w₂ i.i.d. products of random factor-base primes. Then
  Pr[w is a square mod every p | N] = Σ_g Pr[χ(w₁) = g]² ≥ 2^{−w(N)} by Cauchy–Schwarz, χ the vector of Legendre
  symbols (a single product can fail with probability 1, e.g. if every factor-base prime is a non-residue mod p₁).
  Take square roots mod each p by Tonelli–Shanks, lift by Hensel (p odd, w a unit), and combine by CRT with independent
  random signs: X is uniform among the 2^{w(N)} roots. So HRS ≡_R FACT. ∎ (≡_R is equivalence of problems, not of
  costs: HRS at cost T gives factoring at ≈ K·T.)

**The hierarchy** (a partial order, not a chain; from stronger to weaker where an implication is shown):
- CLASS_A, modulus-generic: this implies FACT and prime-field DL with descent;
- CLASS_A for one N;
- CLASS_∞ ≡ FACT ∧ DL;
- HRS ≡ FACT.
- Off this chain: the doubly-smooth lister (suffmap E3) implies FACT under heuristic R, and is an HRS if its outputs
  have hidden roots. No implication between it and CLASS_∞ is shown.

Prop SI (one-sided, intervals) is not on this chain. It supplies half-relations whose residues are tested, and it
needs H_S and R (sicheck shows that a deterministic constructor can make R fail).

**Weakest construction that still factors:** an HRS, the weakest named target, which is no harder than factoring.
CLASS_∞ ≡ FACT ∧ DL sits strictly between it and bounded CLASS unless DL reduces to factoring (open). This is why the project's
sufficient conditions take the form "construct at *cost* below sampling", not "construct at all". Dixon's algorithm
*is* an HRS at its own cost. QS and NFS are HRS-like under heuristic R: their square roots are determined by the
relation, and fruitfulness is heuristic.

**Where Prop EQ-dom fits** (dsr). Dominant relation rows ≡ deterministic factoring under ERH. EQ4 is the randomized,
smooth-square analogue. Both say that the minimal relation-supply target is factoring itself.

---

## 5. Angle 4: smooth numbers by recursion (descent)

**Lemma DSC (size conservation) [PROVED; trivial].**
- Let 0 < |t| < N and u ≡ t·v (mod N) with u ≠ t·v. Then N | u − tv ≠ 0, so |u| + |t||v| ≥ N and
  max(|u|, |tv|) ≥ N/2.
- The same holds in F_q[X]/(I) with deg I = k, degrees in place of logs: max(deg u, deg t + deg v) ≥ k.
- *Reading.* A "special-q" relation for a target q of size Q has a side of size ≥ N/2 or a cofactor v of size ≥ N/(2Q).
  A descent step from q never produces a relation of total size |u|·Q·|v| below N/2. Smaller targets are not easier in the
  one-sided Z/N model: every node of a descent tree is a full-size instance of the class problem of §4 (and CLASS is
  random self-reducible by multiplying t by random squares). In F_{q^k} the same floor holds, so **size-halving is not
  the missing ingredient**.

**What BGJT / KW actually escape with.**
- In F_{q^k} the side of degree ≥ k is the systematic product ∏_{a∈F_q}(X − a) transported by homographies. It is
  *split by construction into atoms from a fixed set of q² + 1 linear polynomials*, and these atoms recur in every
  relation.
- The other side has degree ≤ 2D, independent of k. The axioms are systeq's E0–E5; KW: "applying Frobenius is
  equivalent to a small degree rational function".
- So the transferable requirement, in this note's language, is a family of N-tied identities whose ≥ N/2 side
  splits into atoms from a set of size L^{o(1)}·(#relations) that recur, while the other side is small.

**Why integers do not have it.**
- Additive Frobenius on a faithful Z/N-algebra forces λ(N) | s − 1, so exhibiting one already factors N (systeq
  Prop F). The restricted BGJT compatibility is impossible in rank ≥ 4 when all primes of N exceed s (Prop F′).
- Systematic products vanishing on all scalars have degree ≥ p (Lemma CF). Sparse split families do not exist over R
  (Lemma DES).
- Integer products of shifted pieces ∏(x − a) have pieces that are *new integers for each x*, so they do not recur.
  By Lemma KI they are dead weight unless they recur, and by Lemma AT recurring atoms are no better than primes. The
  integer analogue of "the linear polynomials" is the set of small primes, and a product of small primes ≥ N/2 tied
  to N is a Prop SI / CLASS output, the very thing being constructed.
- Polynomial selection does give an N-tied presentation, but its norm floor is N^{1/(d+1)} (RES/G), which depends on
  N: the classical integer descent (NFS-DL) is L[1/3].
- The open escape remains systeq §6: dense split families with an N-dependent compression (a multigrade
  Prouhet–Tarry–Escott problem). systeq §6(a) already asks for pieces ≤ B^{O(1)}, which implies recurrence; this note
  makes the recurrence requirement explicit (Lemmas KI/LP/AT), which rules out
  shifted-integer pieces.

**Prop SI by recursion [no gain from SR; not a lower bound].** To construct u ∈ [kN, kN + h] as m·u′ with m smooth, u′
must lie in [kN/m, (kN + h)/m], a window of the same relative width h/(kN) = N^{−1+o(1)}. At size X′ = kN/m the window
has length h/m. Recursion does not relax the relative precision h/X (Lemma SR is an identity, not a lower bound
against recursive algorithms). For m ≲ L[2/3] the sub-instances exist but have a smaller ν; once m > L[2/3] the
expected number of smooth integers in a sub-window is below one. So recursion buys nothing that we can see.

**Verdict on angle 4:** the quasi-polynomial machinery does not transfer. The precise missing ingredient is an
identity family whose big side is split into *recurring* atoms by construction (Frobenius additivity supplies it in
F_q). The DSC floor holds in both settings; F_q escapes it only because the big side is split for free.

---

## 6. Angle 5: wild ideas (quick kill or keep)

| # | Idea | Kill / keep | Reason |
|---|---|---|---|
| W1 | Orbits with built-in factorisations: Sylvester s_{n+1} − 1 = ∏s_i; x ↦ x(x + 1); Euclid–Mullin | kill | exact identities over Z (no tie to N); reduced mod N they are full-size residues (Dixon) |
| W2 | Dynamics mod N preserving smoothness (x ↦ s·x, s smooth) | kill | a smooth canonical residue stays smooth only until a wrap; each wrap is one Dixon sample of a full-size residue |
| W3 | Tropical / valuation constructions | kill | smoothness is a valuation condition and the tie to N is archimedean; the tropicalised problem is the log-lattice CVP (schnorrcvp) |
| W4 | Derandomisable random walks (Cayley graph on small primes) | kill | dsr Theorem SB / FN: density and Fourier guarantees cannot certify a hit |
| W5 | Prescribe residues of u mod many small moduli (CRT), hope for u mod N | kill | residues mod small moduli say nothing mod N unless their product exceeds u: crtdecode |
| W6 | Use smooth polynomials over F_ℓ (free in F_ℓ[X]) | kill | no ring map Z/N → F_ℓ-algebra for ℓ ∤ N |
| W7 | Steer units ε in a field with large unit rank so that φ(γε) is small mod N | kill | a multiplicative knapsack in (Z/N)^× with full-size outputs; Prop NV / LC |
| W8 | Binomial / factorial products near kN | kill | binomial products = smooth numbers (C(p, 1) = p); large-k ones sparse (Prop FF); factorials: Lemma CF, N^{1/4} |
| W9 | Pre-split candidates into many known pieces (k ≍ u^θ pieces) | keep only as known | gains k^u = L-constant changes, but this is the listing regime (suffmap §1.1, Theorem LS); placement is still Prop SI |
| W10 | Large-prime recurrence engineered by sharing structure (a^n − 1 over Φ_d(a)) | kill | exact on the free side; KI ⇒ at best a square |

Codex's blind ideas are adjudicated in §9.1 (endpoint-conditioned walks / conditional sampling over exponent
vectors: KEEP as the open navigation problem; z ↦ z(z + 1) orbits: killed, each step introduces a new prime; finite
field lifts: killed, the defect ℓR(M) survives evaluation).

---

## 7. Strongest remaining construction target

- **Prop SI** (suffmap §1.4), unchanged: a y-smooth integer in [kN, kN + h], y = L[1/3, κ], h = L[2/3, ν],
  ν > 2/(3κ), at L[1/3, τ] per output with τ < 0.757. Its Bohr-set form (doubly smooth, windows ≤ y) would give
  exponent < 1/3.
- This note narrows it in three ways.
  - "Known factorisation" does not relax it (Prop KS, Lemmas KI/LP/AT, CGPT).
  - Pre-factored families cannot supply it (§2).
  - Residue-class versions of it are index calculus with descent (Theorem EQ). A modulus-generic residue-class
    constructor would break prime-field DL outright (EQ3); a modulus-generic Prop SI constructor would also improve
    prime-field DL, to L_p[1/3, 1.232] (sihard Prop DLI). Neither form escapes that implication.
- **Codex's alternative target** (blind round §6), stated precisely: given odd N and a uniformly random square a
  with its root withheld, output y-smooth u, v with u ≡ a·v (mod N), ln u + ln v ≤ 3 ln N, on an inverse-polynomial
  fraction of inputs, at y = L[1/3, 0.60] and cost L[1/3, 1.20] per call. By EQ1 that factors N in L[1/3, 1.80]
  *rigorously* (no H_S, no heuristic R), and under GRH such quotients exist (JMV). Both sides are constructed, so
  Theorem SS does not apply. The search problem is an endpoint-conditioned walk of length O(x/ln y) on the small
  primes. We list it as the cleanest *interface*; by EQ3 a modulus-generic solution would also give prime-field DL,
  which we take as evidence against it (assessment).
- For descent the target is systeq §6 plus recurrence: a dense, N-tied family of identities whose big side splits into
  atoms from a recurring set.
- **Assessment:** ≤ 1% that any route here beats L[1/3, 1.902].

---

## 8. Experiments (pre-registered in `notes/smoothbuild_prereg.md`, commit `85d70fe`)

### 8.1 E1: Galois recurrence (Prop GR)

A = 2¹⁵ consecutive a, five disjoint ranges from 2¹⁶; primes p ∈ (2¹⁶, 2²²]; S = Σ_p C(c_p, 2).

| Family | G | r_G − 1 | mean S | S / E (Chebotarev) | S / E (arrangement) | P1a / P1b |
|---|---|---|---|---|---|---|
| x⁴ + x + 1 | S4 | 1 | 662.0 | 0.983 | 0.979 | held |
| x⁴ + 8x + 12 | A4 | 1 | 652.8 | 0.970 | 0.973 | held |
| x⁴ − 2 | D4 | 2 | 1264.0 | 0.939 | 0.945 | held |
| Φ₅ | C4 | 3 | 2031.0 | 1.005 | 1.005 | held |
| Φ₈ = x⁴ + 1 | V4 | 3 | 1937.8 | 0.959 | 0.963 | held |
| x⁶ + x + 1 | S6 | 1 | 665.6 | 0.989 | 1.017 | held |
| Φ₇ | C6 | 5 | 3388.0 | 1.006 | 1.005 | held |
| Φ₉ | C6 | 5 | 3349.4 | 0.995 | 0.994 | held |
| uniform [2⁶³, 2⁶⁴) | — | (1) | 669.6 | 0.994 | — | held |

P1c held: Φ₅/S4 = 3.07, Φ₈/S4 = 2.93, D4/S4 = 1.91, Φ₇/S6 = 5.09, Φ₉/S6 = 5.03, A4/S4 = 0.99.
- *Post hoc (Opus critic).* The two lowest ratios are explained by root pairs r′ = −r, which can never both fall in
  the window once p > 2(a₀ + A) (Prop GR (iii)). Splitting S: D4 ± pairs 594.8 observed vs 669.7 heuristic (0.888),
  other pairs 669.2 vs 668.5 (1.001); V4 ± pairs 606.6/670.9 (0.904), other pairs 1331.2/1341.9 (0.992). Φ₅, Φ₇, Φ₉,
  S4 and S6 have no ± root pairs. A toy-range artefact (a₀ ≪ P₂), irrelevant at L-scale.
- *Reading.* The recurrence boost matches the Galois rank minus one, up to that artefact. A4 (2-transitive, non-symmetric) behaves like
  S4, so the law is about the rank, not about the group being small or abelian. Root sieving was cross-checked
  against flint factorisations (`test_root_sieve_matches_factorisation`).

### 8.2 E2: residue-class construction vs sampling

N balanced semiprimes, 6 per size; y = 2^{bits/4}; targets t = r² mod N.

| bits | M0 rate | M0 O/E | M1 rate | M1 O/E (indep. null) | post hoc sample: indep. / coprime null | M2/M1 per split | M2 target / indep. pred. | M4 O/E |
|---|---|---|---|---|---|---|---|---|
| 40 | 0.0068 | 0.970 | 0.0958 | **0.725** | 0.752 / 1.014 | 0.994 | **0.775** | 0.992 |
| 48 | 0.0061 | 0.972 | 0.0957 | 0.775 | 0.773 / 1.012 | 0.979 | **0.784** | 0.998 |
| 56 | 0.0059 | 0.961 | 0.0951 | 0.805 | 0.804 / 1.006 | 1.010 | **0.794** | 0.991 |
| 64 | 0.0061 | 1.041 | 0.0949 | 0.825 | 0.828 / 1.000 | 0.994 | **0.792** | 1.006 |

- P2a held, P2c held, P2e held.
- **P2b failed at 40 bits** (0.725 < 0.75); it held at 48–64. *Post hoc:* the two pieces of a Gauss split are coprime,
  because gcd(u, v) divides N. On a fresh post-hoc sample the independent null gives 0.752–0.828 (0.752 would have
  passed at 40 bits, so the failure itself is noise at the band edge), while a null drawn as coprime same-length pairs
  gives 1.000–1.014 at every size. Coprimality explains the systematic 17–27% deficit completely. It is sicheck's
  Lemma CP effect, and it shrinks with size.
- **P2d failed at every size** (0.78–0.79). *Post hoc:* 16.6–17.1% of the 16 steered splits of t·s are exact
  rescalings (s·u, v) of t's own split, which yield the *identical* relation, and 36–41% of steered smooth hits repeat
  a relation already found for that target. The rescalings alone would predict ≈ 0.95, not 0.78. What explains P2d
  quantitatively is the number λ of *distinct* relations per target (0.909, 0.952, 0.954, 0.970 at 40–64 bits):
  Poisson 1 − e^{−λ} = 0.597, 0.614, 0.615, 0.621 against observed 0.614, 0.639, 0.627, 0.634 (Opus critic's check).
  Independence of re-randomisation, assumed by P2d, is false for small smooth multipliers.
- *Relations per smoothness test* (M0 vs M1): 0.006 vs 0.048 at every size. The ratio-split constructor beats plain
  sampling by ≈ 8 at u = 4. That is ρ(u/2)²/ρ(u) ≈ 2^u-type, an L^{o(1)} factor (Theorem SS: same total tested
  log-size x).
- *Reading.* Free factors (multipliers, placed factors) change the target, not the tested ledger, and often not even
  the relation. Smooth × smooth products are flat in a window mod N.

### 8.3 E3: all primes allowed vs a fixed smoothness bound (CGPT)

12 seeds per size; Ĵ₀ from 200 000 independent integers (grid y = 2^{k/2}).

| bits | Ĵ₀ | y* | median T/Ĵ₀ | min T/Ĵ₀ | median T_{y*}/T (post hoc) | median (largest prime in dependency)/y* | median dependency size |
|---|---|---|---|---|---|---|---|
| 32 | 4 338 | 1 448 | 0.489 | **0.022** | 1.370 | 183 | 46 |
| 40 | 19 807 | 4 096 | 0.511 | 0.491 | 1.303 | 327 | 118 |
| 48 | 85 486 | 11 585 | 0.508 | 0.469 | 1.318 | 896 | 266 |

- P3a held: the medians lie inside CGPT's proved interval [0.441, 0.561]. The shortfall from the conjectured
  e^{−γ} = 0.561 is unexplained (finite size). It is not a grid artefact: the Opus critic re-estimated Ĵ₀ on a
  2^{1/32}-step grid and found the √2 grid overestimates J₀ by only 0.3–0.6%.
- **P3b failed at 32 bits.** Seed 3's 97th integer was 1 075 315 264 = (8·4099)², a single-element dependency.
  The failure event (T/Ĵ₀ < 0.25, i.e. a perfect square among the first ≈ 1084 draws) has probability ≈ 1084·2⁻¹⁶
  ≈ 1.65% per seed, ≈ 18% across 12 seeds. The pre-registered threshold was miscalibrated at 32 bits.
- **P3c held trivially; design flaw disclosed.** T_y ≥ T for every y, and T_y = T once y exceeds the dependency's
  largest prime, so T_best ≡ T by construction. The informative statistic is T at the bound y* fixed in advance:
  1.30–1.37 × T (post hoc). (T_y = T held within the grid cap 2^{bits−2} in all 36 runs.)
- *Reading.* The first dependency uses primes ≈ 200–900 × y* (the large-prime variation happens by itself), but
  it arrives only ≈ 25% sooner than at the best fixed bound. That is the constant-factor gain CGPT proves.

**Deviations.**
- (1) E3 records the largest prime (descriptive, registered) and the size of the first dependency (descriptive, not
  registered).
- (2) The post-hoc statistics above are labelled as such.
- (3) Ĵ₀'s minimum was taken over grid points y ≤ 2²⁶. Larger y have π(y)/f̂ > 3.9·10⁶ > Ĵ₀ trivially.

---

## 9. GPT-6.1-sol rounds

### 9.1 Blind round (`codex/smoothbuild_1.out`, launched before our analysis; two earlier attempts died on
"model at capacity" and were relaunched unchanged)

**Agreement.** Codex reached the same verdict independently on every angle: pre-factored families supply identities,
not factorisations (Φ₁₁(2) = 23·89) and no tie to N; sum-product results are distribution statements without
navigation; integer descent needs a nonzero modular defect tA − B = kN, which forces a side of log-size x − O(1)
(our Lemma DSC); direct recursion on smaller moduli consumes the window budget (our §5, sharper: N′ = ⌈N/D⌉ forces
D ≤ h). Codex's research prior: ≈ 5% for beating 1.902 by a constructive route, ≈ 0.5% for exponent < 1/3.

**Codex results we adopt.**
- *Hidden-root reduction* (§3.1): identical to our EQ1, found independently, including the explicit statement that
  the constructor may be adversarial and know the factors. Its cost bookkeeping (max(κ + τ + d, 2κ) with sparse
  algebra) is used in §7.
- *Prime-only Lemma AT* (§1.4): the same injection for a set of K primes. Our Lemma AT extends it to any
  multiplicatively independent atoms.
- *Exact cyclotomic incidence* (§1.3): Pr[p | Φ_n(a)] = φ(n)/(p − 1) for n | p − 1, else 0; pair collisions enhanced by
  ≈ φ(n). This is the C_n case of our Prop GR, which gives the exact factor r_G − 1 = φ(n) − 1 for p > A.
- *Shared-base recurrence* (§1.3): p occurs in a^n − 1 exactly for ord_p(a) | n; much of it is the inherited divisor
  a^d − 1 (v_p(a^{dt} − 1) = v_p(a^d − 1) + v_p(t)). Consistent with our W10.
- *Bach 1984* (§3.2, abstract verified by codex): DL modulo composite N is polynomial-time equivalent to factoring
  together with DL modulo primes. Our EQ2 is therefore **KNOWN in substance** (the factor-base-representation version
  of Bach's theorem); we keep only the representation bookkeeping.
- *Rank-free DL from representations* (§3.3): represent g^{r_i}h^{s_i}; a primitive kernel vector c makes Σc_is_i
  uniform, so z = log h is recovered with probability φ(M)/M. An alternative to our exchangeability argument in EQ2;
  both avoid the rank heuristic.
- *Covering count* (§2.2): products of r integers from [1, H] meeting every circular interval of length h + 1 need
  r ln H ≥ x − ln(h + 1); for H ≤ y that is r ≳ x/ln y factors. Existence of enough entropy, not an algorithm.
- *Counterexample to "smooth u with smooth residue suffices"* (§3.5): for N = 2^k − 1, 2^{2k} ≡ 1 is a valid but
  useless relation, so fruitfulness needs hidden roots (or heuristic R): the reason EQ4 is phrased with hidden roots.

**Where codex differs.** Codex names the bounded-height hidden-root quotient constructor as the strongest target
(§7 above); we keep Prop SI as the reference target and list codex's as the cleanest interface. Codex does not use
the Galois-rank formulation or CGPT; we do not use its formal quasi-polynomial descent axiom (§4.1), which we record
here: a polynomial-cost elimination with ≤ x^C children each of size ≤ θz (θ < 1), valid for every target modulo N,
would give quasi-polynomial factoring.

### 9.2 Adversarial round (`codex/smoothbuild_2.out`)

(pending)

## 10. Opus critics

(pending)

## 11. Reproduction

```
.venv/bin/python -m src.smoothbuild.e1_recurrence      # ~35 s
.venv/bin/python -m src.smoothbuild.e3_cgpt            # ~2 min
.venv/bin/python -m src.smoothbuild.e2_class           # ~20 s
.venv/bin/python -m src.smoothbuild.posthoc_e2         # ~1 min (post hoc)
.venv/bin/python -m src.smoothbuild.analyze
.venv/bin/python -m pytest -q tests/test_smoothbuild.py
```

## References (as used; verification status in §10)

- E. Croot, A. Granville, R. Pemantle, P. Tetali, On sharp transitions in making squares, Ann. of Math. 175 (2012).
- C. Pomerance, The role of smooth numbers in number theoretic algorithms, ICM 1994; Fast, rigorous factorization and
  discrete logarithm algorithms (1987).
- J. D. Dixon, Asymptotically fast factorization of integers, Math. Comp. 36 (1981).
- E. Bach, Discrete logarithms and factoring (1984); Explicit bounds for primality testing and related problems,
  Math. Comp. 55 (1990).
- M. O. Rabin, Digitalized signatures and public-key functions as intractable as factorization (1979).
- J. Bourgain, S. Konyagin, Estimates for the number of sums and products and for exponential sums over subgroups in
  fields of prime order (2003); J. Bourgain, A. Glibichuk, S. Konyagin (2006).
- M. Z. Garaev; I. Shparlinski; I. Shkredov (sum-product and modular hyperbola results).
- R. Barbulescu, P. Gaudry, A. Joux, E. Thomé, EUROCRYPT 2014; T. Kleinjung, B. Wesolowski, J. AMS 35 (2022).
- A. K. Lenstra, M. S. Manasse, Factoring with two large primes (1994).
- D. Coppersmith, A. Odlyzko, R. Schroeppel, Discrete logarithms in GF(p), Algorithmica 1 (1986).
- Chebotarev density theorem; Burnside's lemma (standard).
