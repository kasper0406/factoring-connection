# smoothbuild: constructing smooth numbers instead of searching for them

Agent `smoothbuild`, phase 24, 2026-10-23. Version 3 (after the codex blind and adversarial rounds and both Opus critics; §9–§10 record the reviews; v1
commit `7e7037c`, v2 `bd4f4b1`).

**Files.**
- Pre-registration: `notes/smoothbuild_prereg.md` (commit `85d70fe`, before any run).
- Code: `src/smoothbuild/`.
  - `core.py`: operational helpers (factor-base trial division, Gauss split, the hidden-root reduction of Prop EQ1)
    and exact checks (Burnside ranks, roots mod p, Lemma AT counts, GF(2) kernels).
  - `e1_recurrence.py`, `e2_class.py`, `e3_cgpt.py`: the three pre-registered experiments (harness).
  - `posthoc_e2.py`: post-hoc diagnostics for E2's two failed predictions (labelled post hoc).
  - `posthoc_m4.py`: post-hoc fixed-window flatness test replacing the uninformative registered M4.
  - `analyze.py`: the prediction table.
- Tests: `tests/test_smoothbuild.py` (10 tests).
- Data: `experiments/raw/smoothbuild_E{1,2,3}.jsonl`, `smoothbuild_E2_posthoc.jsonl`, `smoothbuild_E2_posthoc_m4.jsonl`; summary
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
- **No construction found here beats sampling at L-scale, and none of the standard pre-factored, sum-product or
  descent ideas supplies one.** Each claim below is marked as one of three kinds: *standard method fails* (a concrete
  method, analysed), *excluded within a stated model* (a theorem or model statement), or *open* (not excluded). The
  exclusions are narrower than "construction is impossible": one construction interface, endpoint-conditioned
  navigation in the small-prime Cayley graph, is open and is not priced by Theorem SS at all (§7).
- **"Known factorisation" is not a cheaper currency than smoothness, in the models we can analyse.** A congruence of
  squares needs only atoms that recur, not small primes. That is true, but:
  - an atom that occurs to an odd power in exactly one relation can never be in a dependency (Lemma KI, pure linear
    algebra). A large known prime helps only if it recurs, alone or in a block with other large atoms;
  - an *isolated* large atom P (all other factors smooth) recurs exactly when P's class is written as
    square × smooth (Lemma LP): a prescribed-target representation problem;
  - for a uniformly random tested side, no K multiplicatively independent atoms beat the K smallest primes
    (Lemma AT, PROVED; the consequence is MODEL). If atoms may be chosen after seeing uniform random integers (all
    primes allowed), the first square dependency comes at e^{−γ}J₀ (Croot–Granville–Pemantle–Tetali 2012;
    Balister–Bollobás–Morris 2018), at most a factor e^{γ} ≈ 1.78 sooner than with the J₀-optimal fixed bound.
    E3 measures 0.49–0.51·J₀, and 1.30–1.37 times as many integers for the bound y* fixed in advance;
  - polynomial-value families make large primes recur more often by a bounded factor: the second moment of the root
    count is the Galois rank r_G (Prop GR, KNOWN in substance), so under a root-equidistribution heuristic the
    pair-recurrence rate is (r_G − 1) ≤ d − 1 times that of random integers, the maximum reached by regular groups
    such as cyclotomic Φ_n. E1 matches the heuristic for eight toy families (0.94–1.01; the two lowest are explained
    by ± root pairs).
- **Pre-factored families: the standard ones fail.** Cyclotomic, Aurifeuillian, Lucas/Chebyshev values, binomials,
  products of consecutive integers and norms × units supply identities, not factorisations, and either have no tie to
  N (exact identities over Z), or tie to N through a full-size tested residue (Dixon, L[1/2]), or must land in a
  window near kN, which sparse fixed families cannot (sinonobliv Lemma PL / Prop FF). Steering kN onto a
  Cunningham-type value needs a root of a fixed polynomial mod N (Rabin's problem for x² − c with c random;
  polynomial selection in general). *Open:* adaptive, target-dependent growing products (codex), which no
  fixed-family count excludes.
- **Sum-product gives distribution statements, not navigation.** Bourgain–Konyagin, Bourgain–Glibichuk–Konyagin,
  Garaev and Shkredov bound how products of structured sets spread; flatness is the hypothesis under which
  meet-in-the-middle does not amortise, and a post-hoc fixed-window test finds smooth × smooth products landing near
  multiples of N at 1.00–1.02 of the flat rate. The fast constructive primitive is the 2-dimensional lattice (rational
  reconstruction; cf. Vallée 1991): ratio splits t ≡ u/v with |u|, |v| ≤ √N. Product splits xy ≡ t (x, y ≤ H) are
  exactly divisor pairs of t + kN, i.e. Dixon samples (Prop RP). Both leave total tested log-size ≈ x: L[1/2].
- **Residue-class construction is index calculus (Theorem EQ).**
  - An average-case constructor of smooth representations u ≡ t·v for random squares t, with the roots hidden from
    it, factors N with a proved splitting probability (EQ1; Dixon's argument in oracle form, found also by codex).
  - With unbounded exponents, representation is randomized-equivalent to "factoring plus discrete logs modulo each
    prime-power factor", given that the small primes generate (Z/N)^× (EQ2; Bach 1984 in factor-base form, with a
    rank-free DL step).
  - A modulus-generic constructor also does prime-field index calculus, descent included (EQ3).
  - The weakest named target, a hidden-root relation sampler, is randomized-equivalent to factoring (EQ4). So the
    content of any construction target is its *cost*, not its possibility.
- **Descent: the quasi-polynomial machinery does not transfer by any route we can see.** Every non-identity relation
  modulo N (or modulo I in F_q[X]) has a side at least as large as the modulus (Lemma DSC). F_q escapes that floor
  because a target-conditioned elimination splits the big side for free (systematic equation / Frobenius), with a
  small other side, a decreasing degree, and controlled branching. Z has no such elimination: additive Frobenius on a
  Z/N-algebra already factors N (systeq Prop F), and the integer substitutes analysed here (known-factored atoms,
  shifted products, rational reconstruction of small targets) give no decreasing measure. This sharpens
  Joux–Odlyzko–Pierrot 2014 §6; it is not a lower bound.
- **Strongest remaining construction targets.** (1) Prop SI (a y-smooth integer in [kN, kN + h] at L[1/3, τ],
  τ < 0.757), or its Bohr-set form for exponent < 1/3. (2) Codex's hidden-root short-quotient constructor
  (endpoint-conditioned navigation), which would give rigorous L[1/3, 1.80] factoring by EQ1 and is outside Theorem SS.
  Both, if modulus-generic, would also improve prime-field DL (sihard Prop DLI; EQ3), which we read as evidence
  against them.
- **Assessment.** Probability that a route within this note's scope beats L[1/3, 1.902]: ≤ 1% (ours); ≈ 5% (codex,
  both rounds). The difference is in how much weight to give the open navigation interface.

**Main results.**

| Result | Statement (short) | Status |
|---|---|---|
| Lemma KI | an atom occurring oddly in exactly one row lies in no dependency; each such large atom costs ≤ 1 kernel dimension for fixed rows | PROVED (pure linear algebra; folklore: singleton removal) |
| Lemma LP | an isolated large atom P (smooth complement) recurs ⟺ P's class is written as (square) × (smooth) | PROVED (identity); block recurrence of several large atoms is the large-prime hypergraph, not covered |
| Lemma AT | K multiplicatively independent atoms compose ≤ Ψ(X, p_K) integers ≤ X | PROVED, elementary (prime-only case also found by codex); no precise prior statement found by us or the critics |
| CGPT / BBM | uniform random integers: first square product at T ∼ e^{−γ}J₀, any primes allowed | KNOWN (CGPT Annals 2012 interval; Balister–Bollobás–Morris Annals 2018 limit) |
| Prop GR | irreducible F of degree d: second moment of ρ_F(p) = Galois rank r_G; 1 ≤ r_G − 1 ≤ d − 1, = d − 1 iff G regular; pair-recurrence law under a root-equidistribution heuristic | (i)–(ii) KNOWN in substance (Frobenius–Burnside + Chebotarev; Akbary–Wong 2021 Thm 1.7); (iii) HEURISTIC, E1 0.94–1.01 |
| Prop KS | Prop SI with "known factorisation" in place of "y-smooth": under a joint random-incidence model for the odd-exponent rough parts, ≳ J₀(N) = L[1/2, √2] outputs are needed unless the outputs' odd-exponent support is ≤ y | MODEL (named assumption JI) |
| Prop RP | product splits xy ≡ t (x, y ≤ H) are divisor pairs of t + kN; ratio splits test two √N pieces; both total tested log-size ≈ x | PROVED (identity); v1's domination claim withdrawn |
| Theorem EQ | EQ1 hidden-root class constructor ⇒ factoring; EQ2 CLASS_∞ ≡_R FACT ∧ DL(p^a ∥ N); EQ3 modulus-generic ⇒ F_p index calculus with descent; EQ4 hidden-root sampler ≡_R FACT | EQ1, EQ4 KNOWN in substance (Dixon 1981 in oracle form), proofs here; EQ2 KNOWN in substance (Bach 1984), factor-base form PROVED given generation of (Z/N)^×; EQ3 PROVED (qualitative), cf. GKZ §2, KW Thm 3.1 |
| Lemma DSC | u ≡ t·v (mod N), (u, v) ∦ (t, 1) ⇒ \|u\| + \|t\|·\|v\| ≥ N; same with degrees in F_q[X]/(I) | PROVED (trivial) |

**Experiments** (pre-registered; §8).
- **E1 (Prop GR):** all 8 families and the random control inside their bands. Observed/predicted 0.94–1.01.
  Cyclotomic Φ₅/S4 = 3.07 (law: 3), Φ₇/S6 = 5.09 and Φ₉/S6 = 5.03 (law: 5), D4/S4 = 1.91 (law: 2), A4/S4 = 0.99 (law: 1).
- **E2 (constructed vs tested side):** P2a, P2c held at every size; P2b failed at 40 bits; P2d failed at every size;
  P2e held but is uninformative (design flaw found by codex).
  - Gauss-split pieces were 0.73–0.83 as smooth as independent same-size integers at every size (P2b's registered
    band failed at 40 bits). Post hoc, a coprime same-size null gives 1.00–1.01: the pieces of a split are coprime
    except when a nonunit exposes a factor.
  - Multiplier steering (P2d): per-target success 0.78–0.79 of the independence prediction. Post hoc, 36–41% of
    steered smooth hits repeat a relation already found for that target (17% of steered splits are exact rescalings
    of the unsteered split).
  - P2e (M4) measured only random-target dilation: its lists' products never wrap mod N. A post-hoc fixed-window
    replacement (random smooth lists near N^{0.55}, window [0, W] near multiples of N) gives O/E 1.00–1.02.
- **E3 (CGPT/BBM):** T/Ĵ₀ median 0.49, 0.51, 0.51 at 32/40/48 bits (P3a held).
  - P3b failed at 32 bits (one seed hit a perfect square at its 97th integer); the threshold was miscalibrated (≈ 18%
    a-priori failure probability across 12 seeds from perfect squares alone).
  - P3c held trivially: T_best ≡ T by construction, a design flaw disclosed in §8.
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

So a construction can help in one of two ways: by making the *tested* side small while keeping the *free* side's
atoms recurrent (Theorem SS prices the tested side), or by constructing *every* side, so that nothing is tested
(navigation; Theorem SS does not apply, and only the cost of the construction itself matters; §7). §2–§6 check each
angle against this. Previously excluded routes are not redone:
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
| Aurifeuillian (Brent 1993: Φ_n((−1)^{(n−1)/2}x) = C_n² − nxD_n² for odd squarefree n; Φ_{n/2}(−x²) = C_n² − nxD_n² for even), at x = n·m² | two half-size pieces | none (special N: constant only, snfsrig) | splitting a candidate into k = O(1) known pieces changes ρ by k^u = L^{o(1)} (suffmap §1.1) | constant at most |
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

**Lemma LP (isolated recurrence = class representation) [PROVED; an identity].**
- Let P be an atom, P > y, gcd(P, N) = 1. A relation A² ≡ P^{±1}·w (mod N) with w a y-smooth rational is the same
  thing as a representation of the class of P^{∓1} as (square) × (y-smooth rational). So a second relation in which P
  is the *only* large atom is a solution of a prescribed-target representation problem at P.
- *Scope (codex adversarial).* A second relation containing P may also contain other large atoms; recurrence of a
  *block* is the large-prime graph/hypergraph mechanism, which LP does not cover. Example: N = 77, y = 5, P = 13,
  Q = 17: 12² ≡ 23² ≡ PQ ≡ 67, and the two rows give gcd(12·23 − 13·17, 77) = 11, although neither row isolates P.
- *No reduction to §4's CLASS is claimed.* CLASS is guaranteed only on uniformly random squares. Turning target P into
  such a square needs a factor-base product q with Pq a square class mod every p | N, which cannot be recognised
  without the factors (codex's N = 187 example: no single factor-base prime works, the product 6 does). v1's
  "up to random self-reduction" is withdrawn.
- *Chance recurrence.* If large atoms behave like independent primes with Pr[P] ≍ c/P (HEURISTIC, the large-prime
  model), R relations give ≈ R²Σ_{P∈(y,P_max]}c²/P² ≍ R²/(y ln y) coincidences, dominated by P near y. For atoms of
  size ≥ L[α] with α > 1/3 (for example a prime near kN, or a cofactor of Φ_n(a)) chance recurrence needs
  R ≳ L[α]^{1/2} relations: hopeless. So the known-prime relaxation is the classical large-prime variation (a
  constant, or L^{o(1)}) for atoms near y, and the residue-class problem itself for larger atoms.

**Lemma AT (atom optimality) [PROVED, elementary, presumably folklore (the injection is the standard Beurling-integer
comparison); the prime-only case was found independently by codex].**
- *Statement.* Let a₁ < … < a_K be multiplicatively independent integers > 1 (e.g. pairwise coprime), and p_K the
  K-th prime. Then the number of integers ≤ X that are products of the a_i (with repetition) is at most Ψ(X, p_K).
- *Proof.*
  - The prime-exponent vectors of a₁, …, a_i have rank i, so together they involve at least i distinct primes, all
    ≤ a_i. Hence a_i ≥ p_i. (For pairwise coprime atoms: their least prime factors are distinct.)
  - By independence, distinct exponent vectors e give distinct products.
  - The map e ↦ ∏p_i^{e_i} is injective and ∏p_i^{e_i} ≤ ∏a_i^{e_i} ≤ X. ∎
- *Dependent atoms.* Linear algebra works over a coprime base (Bernstein, "Factoring into coprimes in essentially linear time", J. Algorithms 54 (2005)), and the bound applies with K
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
  J₀(x) = min_y π(y)x/Ψ(x, y) is the count for the best single smoothness bound (this is the min form used by
  Balister–Bollobás–Morris; CGPT define J₀ via a maximiser y₀ of Ψ(x, y)/y). CGPT conjectured T ∼ e^{−γ}J₀;
  Balister, Bollobás, Morris, "The sharp threshold for making squares", Ann. of Math. 188 (2018) 49–143, proved it.
  Attribution: the upper bound (1 + ε)J₀ for a fixed bound is Schroeppel's; the J₀^{1−o(1)} lower bound is Pomerance,
  "Multiplicative independence for random integers" (1996). CGPT Thm 1.3(c): all elements of the first square product
  are y₀²·exp((2 + ε)√(ln y₀ ln ln y₀))-smooth, which explains E3's largest primes ≈ 200–900·y*. The single-large-prime
  limit ≈ 0.7499·J₀ is computed in Croot–Granville–Pemantle–Tetali, ANTS VIII, LNCS 5011 (2008) 1–36: direct prior art
  for Lemma LP's "large-prime variation = a constant".
- So for random-like integers, allowing every prime with its known factorisation beats the J₀-optimal fixed bound by at
  most the factor e^{γ} ≈ 1.78 asymptotically (BBM against Schroeppel's (1 + o(1))J₀). E3 (§8) measures T/Ĵ₀ ≈ 0.50, and 1.30–1.37 × T for the strategy that fixes y in advance.

**Prop GR (Galois recurrence) [(i)–(ii) KNOWN in substance: the k-th moment of the root count is the number of
G-orbits on X^k (Frobenius–Burnside + Chebotarev; e.g. Akbary–Wong, IJNT 17 (2021) Thm 1.7, cited by codex, not
checked by us); proof recalled here; (iii) HEURISTIC; E1].**
- *Statement.* Let F ∈ Z[x] be irreducible of degree d ≥ 2, G its Galois group acting on the d roots, and
  r_G = number of G-orbits on ordered pairs of roots (the rank). For p ∤ disc(F)·lc(F) let ρ_F(p) = #roots mod p.
  - (i) The natural-density average of ρ_F(p)(ρ_F(p) − 1) over primes equals r_G − 1.
  - (ii) 1 ≤ r_G − 1 ≤ d − 1. The lower bound is attained iff G is 2-transitive, the upper bound iff G acts
    regularly (|G| = d; e.g. Φ_n, every abelian G).
  - (iii) (HEURISTIC RE: for p ∈ (P₁, P₂], the positions of the distinct roots of F mod p, reduced modulo p relative
    to the interval [a₀, a₀ + A), behave like independent uniform points of Z/p.) For A ≤ P₁, the number of pairs
    a ≠ a′ in the interval with F(a), F(a′) sharing a prime p ∈ (P₁, P₂] is then ≈ (r_G − 1)·C(A, 2)·Σ_p p⁻²; for
    independent uniform integers it is ≈ C(A, 2)·Σ_p p⁻².
  - RE is a genuine hypothesis, not implied by Chebotarev. It fails for root pairs r′ = −r when the interval sits near
    0 modulo every p in the band: then p | a + a′ < 2(a₀ + A) is impossible for p > 2(a₀ + A) (Opus critic: E1's deficit
    for D4 and V4, §8.1). Taking a₀ ≫ P₂ does not rescue it (codex: for F = x² + 1 and a₀ the product of all primes in
    the band, the interval is ≡ [0, A) mod every such p and the pair count is 0). Uniformity in F (degree or field
    growing with N) is a further assumption. Equidistribution of roots: Hooley 1964 (all moduli, on average);
    Duke–Friedlander–Iwaniec 1995, Tóth 2000 (prime moduli, quadratic F only).
  - Collisions among *all* family values are not the collisions among *retained* partial relations; the conditional
    distribution is a further modelling step.
- *Proof of (i)–(ii).*
  - For unramified p, ρ_F(p) = fix(Frob_p). By Chebotarev the density average of the class function fix(g)² − fix(g)
    is (1/|G|)Σ_g (fix(g)² − fix(g)).
  - Burnside on X × X gives (1/|G|)Σ fix(g)² = r_G, and (1/|G|)Σ fix(g) = 1 by transitivity.
  - r_G − 1 counts orbits on off-diagonal pairs. By transitivity these correspond to the orbits of a point stabiliser
    G_x on X ∖ {x}. There are at least 1 and at most d − 1 of them, with d − 1 iff G_x fixes every point, i.e.
    G_x = 1 (faithful action), i.e. |G| = d. ∎
- *Reading (under RE).*
  - Polynomial-value families make large primes recur up to d − 1 times more often than random integers.
    Cyclotomic families attain the maximum. This is a statement about eight toy-tested fixed families and the RE
    model, not a universal recurrence barrier.
  - For Φ_n the two effects are explicit: an *eligible* prime (p ≡ 1 mod n) has φ(n) roots and a pair rate
    φ(n)(φ(n) − 1) times that of a random prime, but only a 1/φ(n) share of primes is eligible. The net factor is
    φ(n) − 1.
  - For families of degree d ≤ (ln N)^{O(1)} this is an L^{o(1)} factor.
  - The same families are *less* smooth at small primes (only p ≡ 1 mod n divide Φ_n(a) beyond p | n): systeq
    measured R = 0.61 and 0.25.
  - Families with unbounded recurrence, such as {a^n − 1 : n ≤ n₀} sharing Φ_d(a), are exact identities on the free
    side; their recurrence is free but carries no information about N.
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

*Standard methods fail:* every family analysed either has no tie to N (exact identities), or ties to N through a
full-size tested residue (Theorem SS: L[1/2]), or must hit windows (Prop FF: exponential for sparse fixed families).
*Within stated models:* "known factorisation" adds at most a constant for uniform random integers (CGPT/BBM), a factor
≤ d − 1 in pair recurrence for polynomial values under RE (Prop GR), and isolated recurrence of a large atom is a
prescribed-target representation problem (Lemma LP). *Open:* adaptive, target-dependent growing products with
controlled odd-exponent support (codex); nothing here excludes them, and they are the navigation problem of §7.

**Prop KS (known-factorisation Prop SI) [MODEL under assumption JI].**
- *Setting.* A constructor outputs, for each k, some u_k ∈ [kN, kN + h] together with its full factorisation, at
  cost L^{o(1)}. The residue v_k = u_k − kN ≤ h is tested for y-smoothness. Every gcd(u_k, N) is checked first (the
  relaxed interface allows u = kN, whose factorisation factors N; an honest model must exclude such outputs).
- *What matters is the odd-exponent support,* not smoothness: u = s·a² with s y-smooth contributes no large column
  however large the primes of a (codex: 147 = 3·7² = 143 + 4 at y = 5). The relevant relaxation of Prop SI is
  "outputs whose squarefree kernel is y-smooth", which is Prop SI up to a square factor.
- *Assumption JI (joint random incidence).* Conditioned on v_k being y-smooth, the odd-exponent parts of the u_k
  above y are independent across k and distributed like those of independent uniform integers of size kN.
- *Claim (MODEL).* Under JI, every dependency needs the y-rough odd parts of its u_k to multiply to a square, since
  primes > y cannot occur in a y-smooth v. The CGPT/BBM threshold for that projected matrix is then ≳ J₀(N)^{1−o(1)}
  = L[1/2, √2] outputs with smooth residues. (Applying the uniform-integer threshold to the projection onto primes
  > y is itself part of the model, not a corollary of CGPT.)
- *Reading.* Under JI, "known factorisation" in place of "y-smooth" turns Prop SI's 1.2114 into Dixon's L[1/2, √2].
  JI is exactly what an adversarial or structured constructor need not satisfy; a constructor that controls the
  odd-exponent support of its outputs is a Prop SI constructor (up to squares) or is solving a navigation problem.

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
- A post-hoc fixed-window test (§8.2; it replaces the registered M4, which codex showed measured only target
  dilation) counts pairs of random y-smooth integers near N^{0.55} whose product mod N lies in [0, W], the Prop SI
  window near multiples of N: O/E 1.00–1.02 at 40–64 bits.
- **Closure adds nothing.** The relations form a lattice (exponent vectors of products ≡ 1). Multiplying or dividing
  relations adds lattice vectors and creates no new F₂-information. A product of two Prop SI outputs is a relation
  with residue v₁v₂ ≤ h², not a new short one.

**Prop RP (ratio splits vs product splits) [PROVED identity; domination factor MODEL; v1 corrected after the Opus
literature critic].**
- *Ratio splits.* For t ∈ (Z/N)^×, Λ_t = {(u, v) ∈ Z² : u ≡ tv (mod N)} is a lattice of determinant N. Extended Euclid
  stopped at the first remainder below √N finds (u, v) ∈ Λ_t with |u| < √N and |v| ≤ N/r_{i−1} ≤ √N in O(x²) bit
  operations (`gauss_split`). This is rational reconstruction (Thue's lemma for existence; used in index calculus by
  Blake–Fuji-Hara–Mullin–Vanstone 1984 and Coppersmith–Odlyzko–Schroeppel 1986).
- *Product splits.* H_t = {(x, y) : xy ≡ t (mod N)} is the modular hyperbola. Points with x, y ≤ H exist for
  H ≥ N^{1/2+ε} for almost all t, and for all t once H ≫ N^{3/4+ε} (Weil/Kloosterman; Estermann's bound for every
  modulus). With 0 ≤ t < N, a point with 1 ≤ x, y ≤ H is exactly a divisor pair of t + kN for some 0 ≤ k ≤ H²/N. So
  product splits are found by enumerating k and factoring t + kN: polynomial or subexponential time per candidate.
  (v1 claimed "no search faster than ≈ N/H per point"; that was wrong.) Lattice small-root methods add nothing:
  xy − t is a polynomial in z = xy, so Coppersmith-type bounds reduce to |z| < N.
- *Comparison.* A smooth product split is a y-smooth value of t + kN, i.e. one tested integer of log-size ≥ x (Dixon
  sampling); a smooth ratio split tests two pieces of log-size ≈ x/2 each. Both have total tested log-size ≈ x, so
  both are L[1/2] (Theorem SS); the ratio split wins by the L^{o(1)} factor ρ(u/2)²/ρ(u), which E2 measures as ≈ 8
  at u = 4 (M1 vs M0). Box-by-box point counts do *not* order the two for individual t (literature critic, N =
  1 000 003, H = 1000: t = 720720 has 24 product points vs 6 lattice points; t = 123457 has 0 vs 24); counted as
  distinct relations, all points with one k give the same relation. ∎
- *Reading.* The sum-product world (products) collapses to Dixon sampling of t + kN, and the lattice world gives two
  tested pieces of total log-size x (Theorem SS: L[1/2], the constant of Dixon). Known 2-D constructions do better
  only by shrinking a *single* tested residue (CFRAC: one residue < 2√N; Vallée 1991: x² mod N ≤ N^{2/3}), still
  L[1/2]. Going below needs more pieces, i.e. a higher-dimensional lattice: Schnorr's route, refuted for that model
  (schnorrcvp).

**What is new relative to sihard/siconstruct.**
- They treated integer products in intervals near kN (the hyperbola over Z) and showed that meet-in-the-middle is
  exponential.
- Here, residue-class products mod N (the modular hyperbola) are divisor pairs of t + kN, and ratio splits are
  polynomial-time; both leave a full-size tested ledger. Placing one free factor by the 2-D lattice is the same operation as steering
  by a free multiplier (E2-M2), and it changes the target, not the size of the tested pieces (P2c: 0.98–1.01).
- E2 also shows a correlation: free smooth multipliers often reproduce the *same* relation. In 17% of steered splits
  the split of t·s is the rescaled split (s·u, v), which is the identical relation. Only 59–64% of steered hits are
  distinct (post hoc).

**Verdict on angle 2:** the cited results supply no constructor at L-scale (standard methods fail). The one
constructive primitive found (2-D lattice) is Lemma RC / rational reconstruction, already priced. Covering theorems
show only that enough entropy exists (codex: r ln H ≥ x − ln(h + 1) factors); turning existence into navigation is
the open problem of §7, not excluded here.

---

## 4. Angle 3: residue-class construction and index calculus

**Definitions.** Fix y and A > 0. A *class constructor* CLASS_A(y; T, δ) takes (N, t) with t ∈ (Z/N)^×. In time T it
outputs exponent vectors (e, f) over {−1} ∪ {primes ≤ y} with ∏p^{e_p} ≡ t·∏p^{f_p} (mod N) and
Σ(|e_p| + |f_p|) ln p ≤ A·x, succeeding with probability ≥ δ when t is a uniformly random square. CLASS_∞ drops the
size bound (exponents of poly(x) bits). A *hidden-root sampler* HRS(y; T) outputs pairs (X_i, w_i), w_i y-smooth,
with X_i² ≡ w_i (mod N), such that conditioned on all w_j and on X_j (j ≠ i), X_i is uniform among the square roots
of w_i. Each pair costs expected time T.

**Theorem EQ.** Let N be odd with ≥ 2 distinct prime factors, none ≤ y, and K = π(y) + 1.
- **EQ1 [KNOWN in substance (Dixon 1981); proof here].** CLASS_A(y; T, δ) ⇒ factoring N in expected time
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
  CLASS_∞ ≡_R FACT ∧ DL(p^a) for every p^a ∥ N, with overhead polynomial in K, x and the output length.
  - (⇐) Factor N and compute discrete logs of the factor-base primes and of t modulo each p^k ∥ N. A solution e of
    Σe_p log p ≡ log t modulo each group exponent exists by generation and is found by linear algebra over Z.
  - (⇒) Factoring is EQ1. For DL of h modulo q = p^a ∥ N (p odd, so (Z/q)^× is cyclic of order M = φ(q)), let g
    generate (Z/q)^×; if h is a non-square replace it by g·h and subtract 1 from the recovered logarithm at the end.
    Make m queries, each t_i ≡ g^{2k_i} (mod q) (k_i uniform mod M) combined by CRT with a uniform square modulo N/q,
    so t_i is a uniform square mod N. At one uniformly random position put the target h·g^{2k} instead, which has the
    same distribution. Each success gives v_i = e_i − f_i ∈ Z^K with v_i·(log ℓ)_ℓ ≡ 2k_i (mod M), plus log h at the
    target position.
  - *Rank step (exchangeability; Opus proof critic).* Over Z/M, the indices i with v_i outside the span of the
    other vectors give a strict chain of submodules of (Z/M)^K, so there are at most K·log₂ M of them. Conditioned on
    the successful queries, the oracle's view does not reveal the target's position, so the target is one of them
    with probability ≤ K·log₂ M/s, s the number of successes. Otherwise v_T = Σc_i v_i (coefficients from a Howell
    form, without factoring M), and log h ≡ Σc_i·2k_i − 2k (mod M). Codex's alternative (§9.1): query g^{r_i}h^{s_i}
    and use a primitive kernel vector, which makes Σc_is_i uniform.
  - *Cost.* O(s/δ) oracle calls with s = O(K·x/ε) for error ε, plus module linear algebra on an s × K matrix over
    Z/M (polynomial in K and x; cubic in K by elimination, cf. KW Thm 3.1). Polynomial in K is not polynomial in x
    when K = L[1/3, κ]. ∎
  - This replaces the rank heuristic of v1; the appeal to Pomerance 1987 is withdrawn (an adversarial oracle's
    relations are not uniform).
  - The bounded version CLASS_A asks in addition for a *short* vector of the relation lattice in a coset. That is a
    Schnorr-type CVP and is not implied by FACT ∧ DL (OPEN).
- **EQ3 [PROVED, qualitative; padding as in dsr Lemma PT, rank step as in EQ2; descent-to-DL without a rank
  hypothesis is explicit in Granger–Kleinjung–Zumbrägel §2 and Kleinjung–Wesolowski Thm 3.1].** A constructor that
  works for all admissible N, or only for N = p·r with r a known prime of size p^{o(1)} above y, gives CLASS modulo
  every large prime p at the same cost per call. So it performs prime-field index calculus with descent at
  O(K·x·T/δ) oracle cost plus the linear algebra of EQ2. At y = (ln N)^{O(1)} and polynomial T that is
  polynomial-time factoring *and* polynomial-time prime-field DL. At L[1/3, τ] it would undercut the 1.902 of
  prime-field DL (Matyukhin 2003) once τ is small (sihard Prop DLI).
- **EQ4 [KNOWN in substance; proof here].** HRS(y; T) ⇒ factoring in O(K(T + x²) + K^ω x), as in EQ1. Conversely,
  given the factorisation, an HRS with T = poly(x) exists: output w = 1 and X a uniformly random square root of 1
  (independent CRT signs mod each p^a ∥ N). Then X is uniform among the 2^{w(N)} roots of w, independently across
  calls. (v1 sampled smooth w that are squares mod every p and took roots; the acceptance probability ≥ 2^{−w(N)}
  is not polynomial when w(N) is large, so that construction was replaced, after codex.) So HRS ≡_R FACT. ∎
  (≡_R is equivalence of problems, not of costs: HRS at cost T gives factoring at ≈ K·T.)
- EQ1 and EQ4 are Dixon's 1981 random-squares fruitfulness argument in oracle form; codex found EQ1 independently.

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
CLASS_∞ ≡ FACT ∧ DL lies between it and bounded CLASS; whether strictly is open (it is, unless DL reduces to
factoring). This is why the project's
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
- *Reading.* A "special-q" relation for a target q of size Q has a side of size ≥ N/2 or a cofactor v of size
  ≥ N/(2Q): the relation's total size |u|·Q·|v| is ≥ N/2. Rational reconstruction halves a generic full-size target
  into two √N pieces, but applied to a small target q it may simply return (q, 1), an exact identity with no descent.
  In F_{q^k} the same floor holds. So **the floor itself is not what separates the two settings**; what separates
  them is whether something cheap splits the big side. (v1 also claimed that every node of an integer descent tree
  is a full-size instance of the root problem; that does not follow from DSC and is withdrawn, codex.)

**What BGJT / KW actually have** (codex adversarial; Joux–Odlyzko–Pierrot 2014; GKZ 2018; KW 2022).
- *Relation collection:* the big side is the systematic product ∏_{a∈F_q}(X − a) transported by homographies, split
  by construction into the fixed factor base of linear polynomials, which recur in every relation.
- *Descent:* for a degree-D target P, the systematic side consists of translates P(X) − a of degree D, a
  *target-dependent* set (BGJT §4, equation (E_m), Prop. 2), followed by elimination within that set; the other side
  has degree ≤ (1 + δ)D ≤ 3D (δ = max deg h_i ≤ 2), independent of k. GKZ add extension fields, quadratic
  elimination, norm-back steps and trap avoidance; KW use an elliptic presentation (Frobenius = translation).
- So the transferable package is: a target-conditioned elimination; a decreasing complexity measure; controlled
  branching and exponent sizes; an efficiently available compressed presentation; control of exceptional targets
  throughout. Recurrence is needed only of the terminal generators, not of every intermediate piece. This is
  essentially Joux–Odlyzko–Pierrot's diagnosis (§6: the new methods "deeply rely on specific properties of
  polynomials which do not seem readily adaptable to numbers"); our DSC adds only that the size floor is common to
  both settings.

**Why integers do not have it (standard ingredients fail; no lower bound).**
- Additive Frobenius on a faithful Z/N-algebra forces λ(N) | s − 1, so exhibiting one already factors N (systeq
  Prop F). The restricted BGJT compatibility is impossible in rank ≥ 4 when all primes of N exceed s (Prop F′).
- Systematic products vanishing on all scalars have degree ≥ p (Lemma CF). Sparse split families do not exist over R
  (Lemma DES).
- Polynomial selection does give an N-tied presentation, but its norm floor is N^{1/(d+1)} (RES/G), which depends on
  N: the classical integer descent (NFS-DL) is L[1/3].
- Shifted-integer products ∏(x − a): for one x the pieces are distinct new integers, but across a growing family they
  overlap (x − a = (x + 1) − (a + 1)), so recurrence is not the obstruction (codex). The unresolved challenge is a
  useful N-dependent compression of the big side together with a decreasing measure.
- The open escape remains systeq §6: dense split families with an N-dependent compression (a multigrade
  Prouhet–Tarry–Escott problem), and more generally any target-conditioned integer elimination meeting the package
  above (codex's formal axiom, §9.1, would give quasi-polynomial factoring).

**Prop SI by recursion [no gain from SR; not a lower bound].** To construct u ∈ [kN, kN + h] as m·u′ with m smooth, u′
must lie in [kN/m, (kN + h)/m], a window of the same relative width h/(kN) = N^{−1+o(1)}. At size X′ = kN/m the window
has length h/m. Recursion does not relax the relative precision h/X (Lemma SR is an identity, not a lower bound
against recursive algorithms). For m ≲ L[2/3] the sub-instances exist but have a smaller ν; once m > L[2/3] the
expected number of smooth integers in a sub-window is below one. Codex (blind §4.5) makes this exact for smaller
moduli: with N′ = ⌈N/D⌉, a guarantee for all outputs forces D ≤ h, so ln N′ ≥ x − ln h − O(1).

**Verdict on angle 4:** the quasi-polynomial machinery does not transfer through any ingredient we could identify
(standard methods fail; no lower bound). The missing ingredient is a target-conditioned integer elimination with a
decreasing measure and a small other side, i.e. a cheap split of the ≥ N/2 side; the DSC floor holds in both settings
and F_q escapes it only because the big side is split for free.

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
| W10 | Large-prime recurrence engineered by sharing structure (a^n − 1 over Φ_d(a)) | kill | exact identities on the free side: recurrence is free but carries no information about N |
| W11 | Endpoint-conditioned navigation: choose small-prime word letters by counting/sampling completions that reach the target (codex) | **keep (open)** | entropy suffices (r ln(2π(y)) ≈ (2 + ε)x against one x-bit constraint); no tested side, so Theorem SS does not price it; no algorithm known; a modulus-generic solution would also give prime-field DL (EQ3) |
| W12 | Adaptive products with controlled squarefree kernel (u = s·a², s smooth) | keep (open) | not excluded by KI, AT or fixed-family counts (codex); it is Prop SI up to squares |

Codex's blind ideas are adjudicated in §9.1 (endpoint-conditioned walks / conditional sampling over exponent
vectors: KEEP as the open navigation problem; z ↦ z(z + 1) orbits: killed, each step introduces a new prime; finite
field lifts: killed, the defect ℓR(M) survives evaluation).

---

## 7. Strongest remaining construction target

- **Prop SI** (suffmap §1.4), unchanged: a y-smooth integer in [kN, kN + h], y = L[1/3, κ], h = L[2/3, ν],
  ν > 2/(3κ), at L[1/3, τ] per output with τ < 0.757. Its Bohr-set form (doubly smooth, windows ≤ y) would give
  exponent < 1/3.
- This note narrows it in three ways.
  - "Known factorisation" does not relax it in the models analysed (Prop KS under JI; Lemmas KI/LP/AT; CGPT/BBM),
    except by a square factor (outputs with y-smooth squarefree kernel).
  - Pre-factored families cannot supply it (§2).
  - Residue-class versions of it are index calculus with descent (Theorem EQ). A modulus-generic residue-class
    constructor would break prime-field DL outright (EQ3); a modulus-generic Prop SI constructor would also improve
    prime-field DL, to L_p[1/3, 1.232] (sihard Prop DLI). Neither form escapes that implication.
- **Codex's alternative target** (blind round §6), stated precisely: given odd N and a uniformly random square a
  with its root withheld, output y-smooth u, v with u ≡ a·v (mod N), ln u + ln v ≤ 3 ln N, on an inverse-polynomial
  fraction of inputs, at y = L[1/3, 0.60] and cost L[1/3, 1.20] per call. By EQ1 that factors N in L[1/3, 1.80]
  *rigorously* (no H_S, no heuristic R), and under GRH such quotients exist with some O(1) in place of 3 (JMV's
  expansion; the constant 3 is not derived here). The relation shape u ≡ a·v is Schnorr's (EUROCRYPT '91). Both sides are constructed, so
  Theorem SS does not apply. The search problem is an endpoint-conditioned walk of length O(x/ln y) on the small
  primes. We list it as the cleanest *interface*; by EQ3 a modulus-generic solution would also give prime-field DL,
  which we take as evidence against it (assessment).
- For descent the target is a target-conditioned integer elimination with a decreasing measure and a small other
  side (systeq §6; codex's axiom in §9.1).
- **Assessment:** ≤ 1% (ours) that any route here beats L[1/3, 1.902]; codex ≈ 5% in both rounds, ≈ 0.5% for an
  exponent below 1/3. The gap is how much weight the open navigation interface (W11) deserves.

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

| bits | M0 rate | M0 O/E | M1 rate | M1 O/E (indep. null) | post hoc sample: indep. / coprime null | M2/M1 per split | M2 target / indep. pred. | M4 O/E (registered; uninformative) | fixed-window O/E (post hoc) |
|---|---|---|---|---|---|---|---|---|---|
| 40 | 0.0068 | 0.970 | 0.0958 | **0.725** | 0.752 / 1.014 | 0.994 | **0.775** | 0.992 | 1.019 |
| 48 | 0.0061 | 0.972 | 0.0957 | 0.775 | 0.773 / 1.012 | 0.979 | **0.784** | 0.998 | 1.000 |
| 56 | 0.0059 | 0.961 | 0.0951 | 0.805 | 0.804 / 1.006 | 1.010 | **0.794** | 0.991 | 1.017 |
| 64 | 0.0061 | 1.041 | 0.0949 | 0.825 | 0.828 / 1.000 | 0.994 | **0.792** | 1.006 | 0.999 |

- P2a held, P2c held.
- **P2e held but is uninformative (design flaw, found by codex).** M4's lists were the 1 024 smallest y-smooth
  integers above 2^{bits/2−2}, so every product a·a′ lies in a band of width ≈ 4.5·10⁻⁷·N with no modular wrap
  (64-bit seed 0: all products in [0.0799922, 0.0799926]·N; at t = 1 the registered window has 0 hits against a
  nominal 400). Multiplying by a random t⁻¹ makes the averaged count agree with the flat expectation whatever the
  list's structure, so M4 measured target dilation, not flatness. *Post-hoc replacement*
  (`posthoc_m4.py`): 1 500 random y-smooth integers in [N^{0.55}, 2N^{0.55}] per N (products wrap ≈ N^{0.1} times),
  counting pairs with a·b mod N in the *fixed* window [0, W] near multiples of N (expected 400 per N, 2 400 per
  size): O/E 1.019, 1.000, 1.017, 0.999 (Poisson s.d. ≈ 0.02).
- **P2b failed at 40 bits** (0.725 < 0.75); it held at 48–64. *Post hoc:* gcd(u, v) divides N, so the two pieces of a
  Gauss split are coprime except when a nonunit exposes a factor (codex: gauss_split(67, 77) = (7, 7)). On a fresh
  post-hoc sample the independent null gives 0.752–0.828 and a null drawn as coprime same-length pairs gives
  1.000–1.014 at every size. Coprimality explains the systematic 17–27% deficit. The registered failure stands; the
  fresh sample only shows that 40 bits sits at the band edge. It is sicheck's Lemma CP effect, and it shrinks with
  size.
- **P2d failed at every size** (0.78–0.79). *Post hoc:* 16.6–17.1% of the 16 steered splits of t·s are exact
  rescalings (s·u, v) of t's own split, which yield the *identical* relation, and 36–41% of steered smooth hits repeat
  a relation already found for that target. The rescalings alone would predict ≈ 0.95, not 0.78. What explains P2d
  better is the number λ of *distinct* relations per target (0.909, 0.952, 0.954, 0.970 at 40–64 bits): Poisson
  1 − e^{−λ} = 0.597, 0.614, 0.615, 0.621 against observed 0.614, 0.639, 0.627, 0.634, i.e. it underpredicts by
  1.2–2.5 percentage points at every size. λ is measured on the same data, so this is a descriptive fit, not a
  validated mechanism (codex). Independence of re-randomisation, assumed by P2d, is false for small smooth
  multipliers.
- *Relations per smoothness test* (M0 vs M1): 0.006 vs 0.048 at every size. The ratio-split constructor beats plain
  sampling by ≈ 8 at u = 4. That is ρ(u/2)²/ρ(u) ≈ 2^u-type, an L^{o(1)} factor (Theorem SS: same total tested
  log-size x).
- *Reading.* Free factors (multipliers, placed factors) change the target, not the tested ledger, and often not even
  the relation. In the post-hoc fixed-window test, smooth × smooth products land near multiples of N at the flat
  rate.

### 8.3 E3: all primes allowed vs a fixed smoothness bound (CGPT)

12 seeds per size; Ĵ₀ from 200 000 independent integers (grid y = 2^{k/2}).

| bits | Ĵ₀ | y* | median T/Ĵ₀ | min T/Ĵ₀ | median T_{y*}/T (post hoc) | median (largest prime in dependency)/y* | median dependency size |
|---|---|---|---|---|---|---|---|
| 32 | 4 338 | 1 448 | 0.489 | **0.022** | 1.370 | 183 | 46 |
| 40 | 19 807 | 4 096 | 0.511 | 0.491 | 1.303 | 327 | 118 |
| 48 | 85 486 | 11 585 | 0.508 | 0.469 | 1.318 | 896 | 266 |

- P3a held: the medians lie inside CGPT's proved interval [0.441, 0.561]. The shortfall from the limit e^{−γ} = 0.561
  (proved by Balister–Bollobás–Morris 2018) is a finite-size deviation. It is not a grid artefact: the Opus critic re-estimated Ĵ₀ on a
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
- (4) The registered M2 null (one same-bit-length pair per steered split) was not implemented; P2c compares M2 with
  M1 directly, which is what it registered, but the M2 O/E against a null is missing (found by codex).
- (5) v1 of `e3_cgpt.py` appended to its raw file, so re-running duplicated rows; it now overwrites. The committed raw
  file has no duplicate keys (codex checked).
- (6) `hidden_root_class_factor` sampled r from [2, N − 2]; it now samples the full range [1, N − 1] so that r is
  exactly uniform among the roots (codex). It is used only in tests.

---

## 9. GPT-6.1-sol rounds

### 9.1 Blind round (`codex/smoothbuild_1.out`, launched before our analysis; two earlier attempts died on
"model at capacity" and were relaunched unchanged)

**Agreement and disagreement.** Codex agreed that none of the examined standard methods supplies a fast
constructor: pre-factored families supply identities,
not factorisations (Φ₁₁(2) = 23·89) and no tie to N; sum-product results are distribution statements without
navigation; integer descent needs a nonzero modular defect tA − B = kN, which forces a side of log-size x − O(1)
(our Lemma DSC); direct recursion on smaller moduli consumes the window budget (our §5, sharper: N′ = ⌈N/D⌉ forces
D ≤ h). But codex declined the broader exclusions: it kept adaptive growing families, endpoint-conditioned
navigation and growing-piece elimination explicitly open, and warned that parameter counts are not computational
lower bounds. v2 of this note overstated the agreement; v3 adopts codex's scoping (standard method fails / excluded
in a model / open). Codex's research prior: ≈ 5% for beating 1.902 by a constructive route, ≈ 0.5% for exponent
< 1/3.

**Codex results we adopt.**
- *Hidden-root reduction* (§3.1): identical to our EQ1, found independently, including the explicit statement that
  the constructor may be adversarial and know the factors. Its cost bookkeeping (max(κ + τ + d, 2κ) with sparse
  algebra) is used in §7.
- *Prime-only Lemma AT* (§1.4): the same injection for a set of K primes. Our Lemma AT extends it to any
  multiplicatively independent atoms.
- *Exact cyclotomic incidence* (§1.3): Pr[p | Φ_n(a)] = φ(n)/(p − 1) for a uniform over units and n | p − 1 (φ(n)/p for a uniform mod p), else 0; pair collisions enhanced by
  ≈ φ(n). This is the C_n case of our Prop GR, which gives the exact factor r_G − 1 = φ(n) − 1 for p > A.
- *Shared-base recurrence* (§1.3): p occurs in a^n − 1 exactly for ord_p(a) | n; much of it is the inherited divisor
  a^d − 1 (v_p(a^{dt} − 1) = v_p(a^d − 1) + v_p(t), lifting the exponent, p odd, p | a^d − 1). Consistent with our W10.
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

### 9.2 Adversarial round (`codex/smoothbuild_2.out`, on v2): MAJOR REVISION, all points adjudicated

| # | Codex point | Severity | Adjudication |
|---|---|---|---|
| 1 | KI: say "occurs oddly in exactly one row"; KI compares kernels for fixed rows, not costs; no "at best a square" | MINOR | accepted (§2.2, W10) |
| 2 | LP: a second relation may contain a *block* of large atoms (N = 77 example); the reduction to random-square CLASS needs an unrecognisable square-class adjustment (N = 187 example) | MAJOR | accepted: LP restricted to isolated atoms, reduction claim withdrawn |
| 3 | GR (iii): "a₀ ≫ P₂" is not sufficient (x² + 1 with a₀ = product of the band's primes gives 0 pairs); uniformity in F; retained-relation conditioning; cite Akbary–Wong for the moments | MAJOR | accepted: hypothesis RE stated explicitly, reading restricted to RE and the tested families |
| 4 | RP: the "every box" domination is false (N = 77, t = 6); factoring t + kN is a faster product search | FATAL | accepted (also raised by the literature critic); RP restated as an identity, domination withdrawn |
| 5 | KS: CGPT + H_S does not give the projected threshold; outputs s·a² have no large odd columns; u = kN would factor N | MAJOR / FATAL to the unqualified verdict | accepted: KS is MODEL under the named assumption JI; relaxation phrased via the squarefree kernel; gcd check added to the interface |
| 6 | EQ1 implementation excluded r = ±1; EQ2 prime powers, CRT with p^a, non-square adjustment, linear-algebra cost, K vs log N; EQ3 cost; EQ4 converse not polynomial (2^{−w(N)}): use w = 1; "strictly" → "with strictness open" | MINOR / MAJOR | all accepted (code fixed; EQ2–EQ4 rewritten) |
| 7 | DSC does not show equal difficulty of descendants; BGJT's descent side is target-dependent translates P(X) − a, not the linear factor base; shifted pieces do overlap across families; missing literature (Joux–Pierrot, JOP) | MAJOR | accepted: §5 rewritten around the elimination package; v2's "recurring atoms" description and "every node is full-size" withdrawn |
| 8 | Verdict language exceeds the qualifications; adaptive growing products, navigation and growing-piece elimination remain open; separate "fails / excluded in a model / open" | MAJOR | accepted: §0, §2.4, §6 (W11, W12), §7 |
| 9 | CGPT's conjecture was proved by Balister–Bollobás–Morris 2018; Shparlinski Thms 13, 16 for the hyperbola ranges; Garaev's results are fixed-power | MAJOR / MINOR | accepted |
| 10 | M4 measured target dilation, not flatness; coprimality claim has exceptions; Poisson fit overstated (pp, descriptive); M2 null missing; E3 appended to its raw file | MAJOR / MINOR | accepted: P2e declared uninformative and replaced by a post-hoc fixed-window test (O/E 1.00–1.02); deviations (4)–(6) |
| 11 | "Same verdict on every angle" overstated the blind round's agreement | MAJOR | accepted (§9.1) |

Codex found no constructor meeting τ < 0.757 and none below exponent 1/3; its points invalidated exclusions, not the
"no construction found" verdict. Its probability: ≈ 5% (beats 1.902), ≈ 0.5% (exponent < 1/3).

## 10. Opus critics

**Proof critic (read-only, on v1): ACCEPT WITH FIXES; all applied in v2.** Main points: KI's reading overstated (it
compares fixed rows); EQ2's DL direction was labelled proved but rested on a rank heuristic, now replaced by the
critic's exchangeability argument; EQ4's converse sampler had a gap (single products can never be squares mod p₁);
LP's link to CLASS loose; AT extends to multiplicatively independent atoms; GR (iii) misses ± root pairs, which the
critic verified explain E1's D4/V4 deficits exactly; RP's original proof incomplete; §5 overstated SR and
size-halving; E2's "mostly re-finds" contradicted by the data (59–64% distinct), with the Poisson check supplied;
E2 table mixed samples; E3's P3a grid explanation false (the √2 grid is within 0.3–0.6%); P3b's probability
misquoted (≈ 18% across seeds, not 0.15%); hierarchy is a partial order.

**Literature critic (read-only, on v2): OK WITH FIXES; all applied in v3.** Main points: Prop RP's search claim wrong
(product splits are divisor pairs of t + kN); CGPT's conjecture proved by Balister–Bollobás–Morris (Ann. of Math. 188
(2018)); J₀ definitions (CGPT's maximiser form vs BBM's min form); attributions to Schroeppel, Pomerance 1996,
CGPT ANTS 2008 (single-large-prime limit ≈ 0.7499·J₀, prior art for LP's reading); §5 is a sharpening of
Joux–Odlyzko–Pierrot §6; BGJT's other side is ≤ (1 + δ)D ≤ 3D, not 2D; Vallée 1991 (2-D lattice construction of small
squares) is prior art for "the lattice is the constructive primitive"; "CFRAC/Dixon" conflation (CFRAC tests one
residue < 2√N, a better ledger than a ratio split); rational reconstruction citations; JMV's constant 3 not derived;
Prop GR (i)–(ii) standard; Brent's Aurifeuillian forms; EQ1/EQ4 are Dixon's argument; Estermann's bound holds for
every modulus. The critic found no prior statement of Lemma AT and suggests "presumably folklore".

## 11. Reproduction

```
.venv/bin/python -m src.smoothbuild.e1_recurrence      # ~35 s
.venv/bin/python -m src.smoothbuild.e3_cgpt            # ~2 min
.venv/bin/python -m src.smoothbuild.e2_class           # ~20 s
.venv/bin/python -m src.smoothbuild.posthoc_e2         # ~1 min (post hoc)
.venv/bin/python -m src.smoothbuild.posthoc_m4         # ~1 min (post hoc, fixed-window flatness)
.venv/bin/python -m src.smoothbuild.analyze
.venv/bin/python -m pytest -q tests/test_smoothbuild.py
```

## References (verification status: "critic" = checked in the primary source by the Opus literature critic or codex; otherwise from memory)

- E. Croot, A. Granville, R. Pemantle, P. Tetali, On sharp transitions in making squares, Ann. of Math. 175 (2012)
  1507–1550, Thm 1.2, Thm 1.3(c) (critic); Running time predictions for factoring algorithms, ANTS VIII, LNCS 5011
  (2008) 1–36.
- P. Balister, B. Bollobás, R. Morris, The sharp threshold for making squares, Ann. of Math. 188 (2018) 49–143,
  Thm 1.1 (critic).
- C. Pomerance, Multiplicative independence for random integers, Analytic Number Theory vol. 2, Progr. Math. 139
  (1996) 703–711; Fast, rigorous factorization and discrete logarithm algorithms (1987).
- J. D. Dixon, Asymptotically fast factorization of integers, Math. Comp. 36 (1981) 255–260.
- E. Bach, Discrete logarithms and factoring, UC Berkeley Tech. Rep. UCB/CSD-84-186 (1984) (abstract, critic and
  codex); Explicit bounds for primality testing and related problems, Math. Comp. 55 (1990) 355–380.
- B. Vallée, Generation of elements with small modular squares and provably fast integer factoring algorithms,
  Math. Comp. 56 (1991) 823–849 (abstract, critic).
- M. A. Morrison, J. Brillhart, A method of factoring and the factorization of F₇, Math. Comp. 29 (1975).
- P. S. Wang, M. J. T. Guy, J. H. Davenport, P-adic reconstruction of rational numbers, SIGSAM Bull. 16(2) (1982);
  I. F. Blake, R. Fuji-Hara, R. C. Mullin, S. A. Vanstone, SIAM J. Algebraic Discrete Methods 5 (1984);
  D. Coppersmith, A. Odlyzko, R. Schroeppel, Discrete logarithms in GF(p), Algorithmica 1 (1986) 1–15.
- M. O. Rabin, Digitalized signatures and public-key functions as intractable as factorization, MIT/LCS/TR-212 (1979).
- C. P. Schnorr, Factoring integers and computing discrete logarithms via Diophantine approximation, EUROCRYPT '91.
- D. Jao, S. D. Miller, R. Venkatesan, Expander graphs based on GRH with an application to elliptic curve
  cryptography, J. Number Theory 129 (2009) 1491–1504 (abstract, critic; Thm 1.1, codex).
- J. Bourgain, N. Katz, T. Tao, GAFA 14 (2004) 27–57; J. Bourgain, S. Konyagin, C. R. Math. Acad. Sci. Paris 337
  (2003) 75–80 (codex); J. Bourgain, A. Glibichuk, S. Konyagin, J. London Math. Soc. (2) 73 (2006) 380–398.
- M. Z. Garaev, A. A. Karatsuba, arXiv:math/0503164 (critic); M. Z. Garaev, arXiv:1701.07119 (codex);
  J. Bourgain, M. Z. Garaev, S. Konyagin, I. Shparlinski, arXiv:1203.0017 (codex); I. Shparlinski, Modular
  hyperbolas, Japan. J. Math. 7 (2012) 235–294, Thms 13, 16 (codex); T. Estermann, On Kloosterman's sum,
  Mathematika 8 (1961) 83–86; A. Ayyad, T. Cochrane, Z. Zheng, J. Number Theory 59 (1996) 398–413.
- A. Akbary, P.-J. Wong, On the moments of torsion points modulo primes and their applications, IJNT 17 (2021)
  1029–1045, Thm 1.7 (cited by codex, not checked by us); C. Hooley, Mathematika 11 (1964); W. Duke, J. Friedlander,
  H. Iwaniec, Ann. of Math. 141 (1995); Á. Tóth, Duke Math. J. 104 (2000).
- R. P. Brent, On computing factors of cyclotomic polynomials, Math. Comp. 61 (1993) 131–149 (abstract, critic);
  K. Zsigmondy, Monatsh. Math. Phys. 3 (1892); Yu. Bilu, G. Hanrot, P. Voutier, J. reine angew. Math. 539 (2001)
  75–122 (codex).
- A. K. Lenstra, M. S. Manasse, Factoring with two large primes, EUROCRYPT '90, LNCS 473 (1991); Math. Comp. 63 (1994)
  785–798; P. Leyland, A. K. Lenstra, B. Dodson, A. Muffett, S. Wagstaff, MPQS with three large primes, ANTS-V,
  LNCS 2369 (2002).
- D. J. Bernstein, Factoring into coprimes in essentially linear time, J. Algorithms 54 (2005) 1–30.
- R. Barbulescu, P. Gaudry, A. Joux, E. Thomé, A heuristic quasi-polynomial algorithm for discrete logarithm in
  finite fields of small characteristic, EUROCRYPT 2014, arXiv:1306.4244 (critic, codex); T. Kleinjung,
  B. Wesolowski, J. AMS 35 (2022), arXiv:1906.10668 (critic, codex); R. Granger, T. Kleinjung, J. Zumbrägel,
  Trans. AMS 370 (2018) 3129–3145, arXiv:1507.01495; A. Joux, C. Pierrot, ASIACRYPT 2014; A. Joux, A. Odlyzko,
  C. Pierrot, The past, evolving present and future of discrete logarithm, Springer 2014, §§4, 6 (critic);
  A. Joux, SAC 2013; F. Göloğlu, R. Granger, G. McGuire, J. Zumbrägel, CRYPTO 2013.
- A. Matyukhin, Discrete Math. Appl. 13 (2003).
