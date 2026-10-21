# suffmap: the map of SUFFICIENT conditions for faster factoring

Agent `suffmap`, phase 19, 2026-10-21.

**Files.**
- Code: `src/suffmap/principle.py` (Theorem SS, size–entropy laws, the σ–α table, Prop SI, listing costs; all exact
  sympy/fractions), `src/suffmap/listing.py` (Lemma SL, Prop PC, toy exact checks; harness only).
- Tests: `src/suffmap/test_suffmap.py` (run explicitly, see §8).
- GPT-6.1-sol: blind round `codex/suffmap_1.{md,out}` (§6.1); adversarial round `codex/suffmap_2.{md,out}` (§6.2).
- Reviews: one read-only Opus proof critic and one read-only Opus literature critic (§7).
- No experiment tests a prediction, so there is no pre-registration (task rule). The checks in §5 are exact
  formula checks.
- Cross-references: `schnorrcvp` (lattice / approximate-CVP route) and `abcsunit` (S-unit / abc route) are parallel
  phase-19 agents.
  - `notes/abcsunit_results.md` appeared during this work. Its verdict agrees with E3/H8 here: Prop SU, Lemma FR,
    Prop EX, Prop ABCV and Prop EV.
  - schnorrcvp had only a pre-registration at the time of writing.
  - §3 classifies their objects only at the level needed here and defers the details to them.

**Labels.** KNOWN (published; cited) · PROVED (proof here; elementary unless said) · model-PROVED (exact consequence of a
named model/heuristic) · HEURISTIC (named assumption) · OPEN · assessment (subjective). "No prior statement found by us"
is the strongest novelty claim made.

**Notation.** x = log N, F_α = x^α (log x)^{1−α}, L[α, c] = exp((c + o(1))F_α). Note F_α² = F_{2α}·log x.

---

## 0. Verdict

**Plain language.**
- We found **no sufficient condition that is both plausible and gives an exponent below 1/3.** Every sufficient
  condition on the map is one of three kinds:
  1. **as strong as polynomial-time factoring in disguise**, and usually also a breakthrough elsewhere (prime-field
     discrete logarithms, algebraic circuit complexity, quantum simulation);
  2. **only an L[1/3] constant**, some of them tied to concrete, attackable problems;
  3. **already refuted for typical N**, or a model in which the target is absent under the standard counting
     heuristics.
- **The organising theorem is the size–entropy law (Theorem SS).** For a method that *samples* candidates and waits for
  them to be smooth over *one shared factor base*, with every examined candidate tested in full (model S, §1.1),
  everything is decided by one function λ(C): the smallest total log-size of the smoothness integers when e^C
  candidates are examined. The cost is e^{C*}, where C* is the least C with C² ≥ 2λ(C)·log λ(C). For a power
  law λ = A·x^a(log x)^{1−a−b}C^b this gives **α = a/(2 − b)** and **c^{2−b} = 4αA**. The formula reproduces Dixon √2,
  QS 1, GNFS (64/9)^{1/3}, SNFS (32/9)^{1/3} and t-dimensional sieving exactly (§5). An L[α] sampling method needs
  candidates of log-size σ = 2α in the exponent: **σ = 1/2 for α = 1/4, 2/5 for 1/5, 1/5 for 1/10**, and
  λ = O((log x)²/log log x) for polynomial time. In this model faster linear algebra (any ω ≤ 2) does not move the
  cost.
- **Scope of model S.** MNFS (shared rational filtering, then V fields on the survivors; L[β_a + v] columns), unions
  of many polynomial pairs, and the factory lie *outside* S. That is why MNFS reaches 1.9019 < 1.9230 and feels ω.
  Theorem SS is a statement about S, not about all relation methods.
- **Listing escapes the size law** (Theorem LS). If smooth objects are *listed* rather than waited for, size enters only
  through existence (box sizes), and the cost is the per-output listing cost. A complete search oracle gives a lister
  with ≤ 2K + 1 calls (Lemma SL, PROVED). A *promise* constructor (one that only works on intervals guaranteed to contain
  a smooth number) is worth exactly one output per tile of its guaranteed length (Prop PC, PROVED).
- **The cleanest new sufficient condition** is about constructing smooth numbers in a given short interval (Prop SI).
  Suppose a y-smooth integer can be constructed in polynomial time in every interval [kN, kN + h], with
  h = L[2/3, ν], y = L[1/3, κ] and ν > 2/(3κ). Then N factors in **L[1/3, (16/9)^{1/3} + ε] = L[1/3, 1.2114 + ε]**
  for every ε > 0 (HEURISTIC: the residue behaves randomly and the dependencies are fruitful). That is below MNFS 1.9019 and SNFS 1.526, with no number field at all. Single intervals can
  never beat exponent 1/3 (model-PROVED: max(α, 1 − 2α)). The sparse-set version, smooth numbers in the union of
  windows around many multiples of N (a Bohr set), is a *different* hypothesis (an oracle on unions of windows). With
  windows of width ≤ y it is the doubly lister E3 and gives heuristic polynomial time, even from a promise constructor.
  Prop SI itself is a one-sided form of Schnorr's target. No prior statement of the 1.2114 reduction was found by us, and codex found the
  same constant independently.

**The map in one table** (details §2, high-profile verdicts §3, ranking §4).

| Tier | Sufficient condition (short) | Gives | Believed |
|---|---|---|---|
| E | sampling family with C* = O(F_α), α < 1/3 (rank-r aggregate, toric) | L[α], formally 1/4 | false for typical N (G_r) |
| E | one-sided Bohr-set lister (smooth u, u mod N ≤ L[2α]) at L[α, O(1)] per output | L[α] (with u mod N ≤ y it is the doubly lister: poly time) | false |
| E | doubly-smooth lister / S-unit enumeration / relation-lattice vectors / Cayley-graph navigation | heuristic poly time (navigation to arbitrary targets: rigorous randomized poly time) | false; in coset form it would also give prime-field DL |
| E | short SLPs for Pochhammer polynomials (⇐ permanent in constant-free VP) | factoring ∈ P/poly, **rigorous** | false (τ-conjecture) |
| E | P = NP ∩ coNP; BQP ⊆ BPP | P; BPP | false |
| E | algebraic-group families whose annihilation separates CRT components with probability L_p[α]^{−O(1)} | L_p[α] | no structure known |
| C | smooth numbers constructed in single short intervals near kN (Prop SI) | **1.2114** | open; a reformulation of one-sided Schnorr CVP (no reduction either way) |
| C | complete search / sampler for smooth numbers in NFS rational rows | 1.6386 | open |
| C | the same on the algebraic side | inf 0.7631 | open |
| C | beyond-Johnson output-sensitive CRT list decoding; H_sel violations | c(θ) | open, implausible |
| C | sparse F₂ linear algebra with ω < 2 | 1.8595 (MNFS only) | open, plausible-ish |
| C | small-height polynomial pairs for general N | (32κ/9)^{1/3} | **false for almost all N** (GAE) |
| C | Kirchner's ramified fields | 1.2114 (special N), 1.5263 (most N) | special: false for a.a. N; relaxed: unconstructed |
| X | Coppersmith from a fraction f < 1/2 of p's bits | N^{f/2} (exponential) | open |

**High-profile verdicts (§3).**
- Matomäki–Radziwiłł / Harper-type existence of smooth numbers in short intervals: **irrelevant to speed** (no
  algorithm, used by none).
  They become the *existence half* of the promise-constructor route (Prop PC, Prop SI) only at the scale
  y = L_X[1/3], h = L_X[2/3], far beyond current results.
- Deterministic construction of a smooth number in a given interval: **sufficient for L[1/3, 1.2114]** at the Prop SI
  scale, and **for heuristic polynomial time** in the Bohr-set form.
- Polymath4 (primes): **unrelated**.
- Sum-product: **unrelated**.
- Expansion: existence of short paths is **unrelated**; constructive navigation is **sufficient** (= doubly listing).
- Small-discriminant fields with a root mod N: **constant only**, and the height form is refuted for almost all N.
- abc and S-unit counting: **unrelated**; algorithmic S-unit solving is **sufficient**.
- Approximate CVP in prime lattices: a *method* for the listing targets; whether any approximation factor suffices is
  schnorrcvp's question.
- Coppersmith: **exponential regime only**.
- CRT list decoding beyond Johnson: **constant only**, unless doubly.
- Sparse linear algebra: **constant only**.
- Factorial / Pochhammer straight-line programs: **sufficient (nonuniform P), and a landmark either way**.

**Ranking (§4).**
- *Most worth attacking (assessment; codex prefers C2 or C5 for clearer interfaces):* the Prop SI problem. Construct, or prove hard in a model, L[1/3]-smooth integers in intervals
  of length L[2/3] near multiples of N. It is concrete, no barrier beyond method classes is known, it is
  modulus-generic, and a solution would cut the factoring constant from 1.90 to 1.21.
- *Landmark either way:* the straight-line complexity of factorials / Pochhammer polynomials. A proof of easiness gives
  nonuniform polynomial-time factoring. A proof of hardness, quantified over power-of-two multiples, gives a
  constant-free permanent lower bound (Bürgisser 2009).
- *Highest-profile direct exponent target:* the doubly-smooth / Bohr-set lister (Schnorr's target). Its consequence for
  prime-field DL (Prop DLT) makes it implausible.

---

## 1. Two principles

### 1.1 Models

**Sampling model S.**
- From N, the algorithm generates candidates z from a family 𝓕 and *examines* them, paying ≥ 1 each.
- Each z carries k = O(1) nonzero smoothness integers n₁(z), …, n_k(z), of total log-size λ(z) = Σ log|n_i(z)|.
- Each z also carries a multiplicative link to Z/N: congruences of products, or ring maps as in resbarrier's LOM.
  F₂-dependencies among B-smooth candidates therefore give X² ≡ Y² (mod N).
- All candidates share **one factor base** of π(B)(1 + o(1)) columns, and each examined candidate is tested on *all*
  its k integers, at cost ≥ 1. There is no cheaper pre-filtering on one component.
- Linear algebra on K columns costs K^ω.

**What S excludes (proof critic M1).**
- MNFS: one shared rational side is sieved over E² points, and only its survivors are tested in V fields. The cost is
  E²(1 + V·L[−r]), not E²·V, and there are L[β_a + v] columns.
- Unions of M polynomial pairs: entropy 2 log E + log M with M separate factor bases. Their law is
  2x/D + D(C − log M)/2 + O(log M), better than 2√(xC) for one pair.
- The factory, which uses precomputation.
- Listing (§1.3).
Theorem SS therefore does not cover the best known algorithm (MNFS). Its role is to say exactly what a *single-family*
sampling method must achieve.

**Heuristic H_S.**
- An examined z is B-smooth with probability ρ(λ(z)/log B)^{1+o(1)}, where −log ρ(u) = (1 + o(1))·u log u.
- This is per examined candidate, independent of how z was selected. That clause is a **strengthening** of heurnfs H,
  which only lets joint events multiply and so allows conditioning. The clause fails for candidates pre-filtered on one
  component, which is why such filtering is excluded from S.
- Smoothness tests and candidate generation cost L^{o(1)}.

**Heuristic R.** π(B)·L^{±o(1)} relations are needed and suffice, and the dependencies are fruitful. For random-like
relations, "needed" is Croot–Granville–Pemantle–Tetali 2012 (up to a constant). For structured families it is an
assumption.

**Size–entropy profile.**
- C_𝓕(λ) := log #{z ∈ 𝓕 : λ(z) ≤ λ}.
- The *law* λ(C) is its inverse: the smallest λ such that e^C candidates of size ≤ λ exist.
- Example: two-polynomial NFS under G has λ(C) = 2√(xC) (heurnfs B1 with C = 2 log E). SNFS has √(2xC).

**Why k = O(1) matters.** If a candidate's smoothness integers are many small pieces (k ≍ u), then ∏ρ(λ_i/μ) is no
longer ρ(λ/μ). The extreme is "smooth by construction": a product of factor-base elements. That is the *listing* regime
of §1.3, not a counterexample to §1.2. systeq §6, ptesplit and nonexact live there.

### 1.2 Theorem SS (sampling) [model-PROVED in S + H_S + R; KNOWN in substance]

Write T(λ) := √(2λ log λ).

**(a) Fixed size [KNOWN].** If every examined candidate has size λ, the cost is exp((1 + o(1))T(λ)). In Pomerance's
notation this is L(X)^{√2+o(1)} with X = e^λ, L(X) = exp(√(log X log log X)) (Pomerance 1982; 2008 survey "Smooth numbers
and the quadratic sieve", MSRI Publ. 44). The optimum is attained at B = exp(T(λ)/2), with e^{T(λ)} candidates.

**(b) L-units.** If λ = s·F_{2α}, the cost is L[α, 2√(αs)], with factor base L[α, √(αs)] and L[α, 2√(αs)] candidates.
Equivalently: **L[α, c] needs candidates of total log-size ≤ (c²/(4α))·F_{2α}**, i.e. σ = 2α.

**(c) Size–entropy law.** Let C*_δ := inf{C : C² ≥ (1 − δ)·2λ(C) log λ(C)} and C* := C*_0. Two regularity
assumptions are used:
- *Levels:* λ_max ≤ exp(x^{O(1)}), so that (1+ε)-adic size levels number L^{o(1)}.
- *Stability:* C*_δ = (1 + o(1))C* as δ → 0, and λ((1 + o(1))C*) = (1 + o(1))λ(C*). This holds for every regularly
  varying law, in particular every law in the catalogue.
- Without stability the theorem fails: a step profile with λ = λ₀ up to just below T(λ₀), then huge, has C* huge while
  cost ≈ e^{T(λ₀)} cannot be excluded under ρ^{1+o(1)} (proof critic m1).
- (lower bound) Every algorithm in model S costs ≥ exp((1 − o(1))C*). This needs **no** assumption on linear algebra,
  but it does use heuristic R (≈ π(B) relations are *needed*). H_S alone does not imply that: a generator could supply
  relations with extra structure (codex). For *raw independent uniform integers* there are theorems on the
  first-square-product threshold: Croot–Granville–Pemantle–Tetali 2012; Balister–Bollobás–Morris (arXiv 1608.03857,
  Thm 1.1) give T(X) ∼ e^{−γ}·min_y π(y)X/Ψ(X, y). Their model allows dependencies through primes outside a prescribed
  base, so they support R only as an analogy. For structured families R is an assumption.
- (upper bound) If 𝓕 is enumerable at L^{o(1)} per candidate, H_S holds as an equality, and ω ≤ 2, the cost is
  ≤ exp((1 + o(1))C*). This needs the population to extend a factor e^{o(C*)} beyond the threshold, which stability
  guarantees, and yield concentration, which follows from H_S's independence by Chernoff. Without stability the
  general statement is unjustified (codex suffmap_2 #4). The clean restricted theorem is for regularly varying laws,
  in particular the power laws of (d).
- So for single-family sampling, **ω < 2 does not help**. This agrees with heurnfs (GNFS stays at 1.9230 at ω = 1.5).
  MNFS gains from ω < 2 only because its V fields have different column structures.

**(d) Power laws.** For λ(C) = A·x^a(log x)^{1−a−b}·C^b with a > 0, 0 ≤ b < 2 and a < 2 − b, so that 0 < α < 1:

  C* = (c + o(1))F_α,   α = a/(2 − b),   c^{2−b} = 4αA.

The log-power 1 − a − b is the only one for which the law is L-homogeneous. Other log-powers move the constant to 0 or
∞ at the same α. If b ≥ 2 there is no subexponential solution.

**(e) Polynomial and quasi-polynomial time.** Here generation and tests must cost poly(x), not L^{o(1)}, and the
level argument loses an additive O(log x) in C (proof critic m2). With these changes:
- C* = O(log x) iff λ(C*) = O((log x)²/log log x).
- λ = (log x)^k gives exp(O((log x)^{k/2}(log log x)^{1/2})), which is quasi-polynomial.

**Proof.**
1. *Optimal B for fixed λ.* With μ = log B and u = λ/μ, the relations needed are B, and each examined candidate succeeds
   with probability ρ(u). So the number of candidates is ≥ exp(μ + u log u)·(1 − o(1)).
   g(μ) := μ + (λ/μ)log(λ/μ) is minimised at μ₀ = √(λ log λ/2)(1 + o(1)). There log u = (½ + o(1))log λ, so
   u log u = μ₀(1 + o(1)) and g(μ₀) = 2μ₀ = T(λ)(1 + o(1)).
   Linear algebra costs B^ω = e^{ωμ₀} ≤ e^{2μ₀} when ω ≤ 2. This proves (a), and the ω-independence for fixed λ.
2. *(b).* Substitute λ = sF_{2α}, μ = bF_α. Then u = (s/b)x^α(log x)^{−α}, log u = (α + o(1))log x, and
   u log u = (αs/b)F_α. Minimising b + αs/b gives b = √(αs) and the value 2√(αs) = 2b.
   This is `optimum_fixed_size`; exact for the tested rationals.
3. *(c), lower bound.* The yield is Σ_{examined z} ρ(λ(z)/μ)^{1+o(1)} ≥ π(B)(1 − o(1)). Split the examined z into
   (1+ε)-adic size levels, with ε → 0 slowly. Some level λ carries a 1/L^{o(1)} share of the yield.
   - It contains ≤ e^{C_𝓕((1+ε)λ)} candidates.
   - It requires ≥ exp(g(μ))·L^{−o(1)} ≥ exp((1 − o(1))T(λ)) examined candidates, by step 1, since sizes vary by a
     factor ≤ 1 + ε.
   - Hence C_𝓕((1+ε)λ) ≥ (1 − o(1))T(λ), and cost ≥ exp((1 − o(1))T(λ)).
   - Minimising over the admissible λ: put C₁ = (1 − o(1))T(λ). Since λ(·) is nondecreasing, C_𝓕((1+ε)λ) ≥ C₁ gives
     λ(C₁) ≤ (1+ε)λ, so C₁² ≥ (1 − o(1))·2λ(C₁) log λ(C₁). Hence C₁ ≥ C*_δ for some δ = o(1), and stability gives
     cost ≥ exp((1 − o(1))C*). Here T((1+ε)λ) = (1 + O(ε))T(λ) was used.
4. *(c), upper bound.* Take λ* = λ(C*), B = e^{T(λ*)/2}, and examine the e^{C*} candidates of size ≤ λ*. The expected
   yield is ≥ e^{C*}ρ(λ*/μ) ≥ B·L^{−o(1)} (one more L^{o(1)} slack in C absorbs it). Linear algebra costs
   B^ω ≤ B² = e^{T(λ*)} ≤ e^{C*}.
5. *(d).* Substitute C = cF_α. Then λ = Ac^b x^{a+αb}(log x)^{1−a−αb}. With a + αb = 2α this is Ac^b F_{2α}, and
   log λ = (2α + o(1))log x. So 2λ log λ = 4αAc^b F_α², while C² = c²F_α². Equality holds iff c^{2−b} = 4αA.
   `verify_law_balance` checks the ratio is exactly 1 for every catalogued law.
6. *(e).* C* = K log x iff 2λ log λ = K²(log x)², and log λ = (2 + o(1))log log x. ∎

**Status.**
- (a) is KNOWN (Pomerance's L(X)^{√2}). The size principle itself, and hence σ = 2α in L-units, is stated explicitly in
  Pomerance 2008 (MSRI 44, eq. (8) and §"A general principle, a moral, and three bullets", p. 78): the analysis
  "can serve, in an almost intact manner, for many factoring algorithms. What may change is the number X."
- The matching exponent-level lower bound for random squares is Pomerance, ICM 1994 (Birkhäuser 1995, pp. 411–422),
  sharpened to constants by CGPT 2012. The lemma number in the ICM paper is not checked.
- The exponent bookkeeping α = (1+α)/2-type is heurnfs B5.
- The fixed-point form c^{2−b} = 4αA, the ω-independence for single families, and (e) are our packaging. No prior
  statement was found by us. They are routine.

**What known methods look like under (d)** (`principle.LAWS`; every value exact):

| Method | law λ(C) | (A, a, b) | α | c |
|---|---|---|---|---|
| Dixon | x | (1, 1, 0) | 1/2 | √2 |
| QS | x/2 | (½, 1, 0) | 1/2 | 1 |
| GNFS, any 2-poly selection under G | 2√(xC) | (2, ½, ½) | 1/3 | (64/9)^{1/3} = 1.9230 |
| SNFS (κ = 1) | √(2xC) | (√2, ½, ½) | 1/3 | (32/9)^{1/3} = 1.5263 |
| sieving dimension t = 3 | 2√(4xC/3) | (4/√3, ½, ½) | 1/3 | (256/27)^{1/3} = 2.1165 |
| aggregate rank-r (resbarrier §4.2, stylised) | ≍ √(x log x) | (A, ½, 0) | 1/4 | √A |

**Reading.**
- Exponent 1/3 is the line a/(2 − b) = 1/3. NFS realises (a, b) = (½, ½): candidate size grows like the square root
  of the entropy, times √x. (2/3, 0) is another point on the line: a fixed size x^{2/3}.
- Sub-1/3 needs either a smaller a, i.e. the candidate's size depends less on N, or a smaller b, i.e. it depends less on
  how many candidates there are.
- **α = 1/4, 1/5, 1/10 at b = 0 need laws λ ≍ x^{1/2}, x^{2/5}, x^{1/5}** (times logs): candidates whose size does not
  grow with the family. At b = ½ they need a = 3/8, 3/10, 3/20.
- Pigeonhole gives only λ(C) ≥ C − O(k log λ), when the candidates' smoothness tuples are distinct. That constrains
  nothing.
- For integer congruences u ≡ v (mod N), u ≠ v, with *both* sides random, max(|u|, |v|) ≥ N/2 (PROVED;
  `congruence_max_lower_bound`). So a = 1 and α ≥ 1/2. Below 1/2 one needs ring maps (NFS) or one side smooth by
  construction (QS: a square; listing: a product).

### 1.3 Theorem LS (listing) [model-PROVED; pilist Table 1.1(A), smoothenum T2–T8 restated]

A *lister* outputs candidates already known to be smooth on one side (one-sided) or on all sides (doubly). Cost is
charged per output, L[α, τ] each. Factor base L[α, κ].

- **(a) Doubly.** The cost is L[α, max(κ + τ, 2κ)], for every α ∈ (0, 1), provided the region contains ≥ L[α, κ]
  relations and setup fits the bound. A fixed-α lister does **not** give polynomial time by shrinking B: L[α, τ] is
  superpolynomial for every fixed α > 0. Polynomial time needs a *separate* hypothesis (smoothenum T8): with
  y = (log N)^K, K > 1, and U = N^A, A > K/(K − 1), there must be polynomial setup and polynomial cost for every
  requested output prefix, including factorisations, duplicate removal and extraction.
- **(b) One-sided.** Suppose the residue side has log-size νF_{2α}. The cost is L[α, max(κ + τ + αν/κ, 2κ)].
  At τ = 0 this is Theorem SS(b) with s = ν: only the residue side counts.
- **(c) Existence fixes box sizes, not cost.**
  - Congruence instances need box area (U/N)·V ≥ L[1 − α]^{Θ(1)} (pilist §1.2). This reduces to U/N ≥ L[1 − α] when
    log V = o(F_{1−α}), i.e. V ≤ y, or V = L[2α] with α < 1/3. It does **not** reduce at α = 1/3 with V = L[2/3]
    (Prop SI), where U/N = L[1/3] suffices.
  - NFS-shaped doubly listing at factor base L[α] needs a box log-size x^γ with γ ≥ 1 − 2α.
  - In both cases the box is never paid.
- **(d) NFS-shaped partial listing.** Listing at E²·L[−θ·penalty] gives the curves c_rat(θ), c_alg(θ), c_both(θ) of
  smoothenum T2/T6/T7 and pilist Table 1.1(B).

**Lemma SL (search ⇒ listing) [PROVED; `listing.list_by_search`].** Let O be a *complete* search oracle for a set S: on
[a, b] it returns some element of S ∩ [a, b], or correctly reports none. Then S ∩ [a, b] (K elements) can be listed with
≤ 2K + 1 calls. A union of J intervals needs ≤ 2K + J calls.
- *Proof.* After each hit z, recurse on [a, z − 1] and [z + 1, b], skipping empty intervals. Each element is returned
  exactly once, because the sub-intervals are disjoint and exclude z. Every call is either the root or one of the ≤ 2
  children of a hit, so there are ≤ 1 + 2K calls. ∎
- *Lexicographic version.* The same holds for any totally ordered region, e.g. pairs (k, v) in the Bohr set
  {kN + v : k ≤ K, |v| ≤ V}, provided the oracle accepts lexicographic intervals. Such an interval is ≤ 3 rectangles.
- *Sampler version.* A sampler returning near-uniform elements of S ∩ [a, b] lists a constant fraction in O(K) calls
  (coupon collector).

**Prop PC (promise constructors) [PROVED, elementary; `AdversarialPromiseOracle`, `tiling_outputs`].** Let O be
guaranteed only on intervals of length ≥ h₀ (each such interval meets S). Let the range have length H, and let every gap
of S be ≤ h₀/4.
- (a) Against an adversarial choice of outputs, any algorithm obtains ≤ 4H/h₀ + 1 distinct elements.
  - The adversary's strategy: return an already-revealed element if the query contains one. Otherwise reveal the
    element of the query nearest the midpoint of the revealed-free gap containing it.
  - Each new point is then ≥ h₀/2 − h₀/4 from every earlier revealed point. If the gap midpoint lies inside the query,
    the new point is within one S-gap of the midpoint. If it lies outside, the new point is within one S-gap of the
    query end nearest the midpoint, and the query has length ≥ h₀.
- (b) Tiling by disjoint intervals of length h₀ achieves ≥ ⌊H/h₀⌋ distinct elements.
- So a promise constructor is worth **one output per h₀-tile**. Since the promise needs h₀ ≥ 4·(max gap), its value
  is governed by the ratio of h₀ to the mean gap 1/δ only under the HEURISTIC that max gap ≤ mean gap·L^{o(1)}
  (Poisson-like gaps; unproved for smooth numbers).
  - If h₀ = δ^{−1}·L^{o(1)}, it is as good as full listing at L-scale (θ = 1).
  - If h₀ = δ^{−1}·L[1/3, ηr], it *thins* the output set to one per tile (rational count penalty (1 + η)r), at no
    sieving cost. This is not the per-output-cost parameter θ of pilist.
  - In a doubly Bohr-set setting the region can be enlarged at polynomial cost, while ranges stay polynomially encoded
    (u ≤ N^A). There the useful number of tiles is what matters, not h₀/δ^{−1}. This requires a *blockwise* promise
    (every queried block of consecutive k contains a valid u), which does not follow from existence in the whole set
    (codex suffmap_2 #8).
  - Prop PC bounds distinct *oracle answers*. The θ readings are count-model substitutions: they need cheap calls and
    the same downstream randomness and fruitfulness assumptions as full listing.
- codex's blind statement "a constructor guaranteed only for sufficiently long intervals does not give that reduction"
  is correct when h₀ ≫ δ^{−1}. Prop PC quantifies it.

### 1.4 Prop SI (single short intervals near multiples of N) [model-PROVED under H_S on the residue + R; `single_interval_program`]

**Setting.**
- Factor base y = L[1/3, κ], intervals I_k = [kN, kN + h] with h = L[2/3, ν], k = 1, 2, ….
- An oracle constructs a y-smooth u_k ∈ I_k in L^{o(1)} time. The oracle's domain is the intervals I_k only. A promise
  constructor suffices, and complete search gives nothing more at L-scale (proof step 7).
- Each u_k gives the relation u_k ≡ v_k (mod N), v_k = u_k − kN ∈ [0, h].

**Claim.**
- (i) The *mean* count of y-smooth integers per interval is L[2/3, ν − 2/(3κ)], which tends to infinity iff
  ν > 2/(3κ). Mean density is not uniform existence. For fixed N, y, h, existence for *every* k is false: CRT can place
  a large prime factor at every position of some interval. So the promise is restricted to the k ≤ L[1/3, C] actually
  used, and it is a hypothesis.
- For a constructor at fixed (κ, ν) with ν > 2/(3κ), the cost is L[1/3, max(κ + ν/(3κ), 2κ)]. The constant
  1.2114 + ε needs constructors with (κ, ν) near κ* = (2/9)^{1/3}, ν* = 2/(3κ*).
- (ii) The cost is L[1/3, max(κ + ν/(3κ), 2κ)]. Its infimum is approached as ν ↓ 2/(3κ) with κ³ = 2/9, giving
  **c = (16/9)^{1/3} + ε = 1.2114 + ε** for every ε > 0. This needs Heuristic R: fruitful, full-rank dependencies.
  A deterministic constructor whose output is a function of kN could make v non-random or the dependencies trivial,
  so the random-like-output assumption is explicit.
- (iii) With smoothness L[α], single intervals cost L[max(α, 1 − 2α)], so **α = 1/3 is optimal and exponent < 1/3 is
  impossible with single intervals**.
- (iv) *A separate, stronger hypothesis, not a consequence of (i)–(ii).* Take an oracle on Bohr-set unions
  ∪_{k∈[k₀, k₀+t₀]} [kN − V, kN + V] with V ≤ y, so the relation is doubly smooth automatically and the region can be
  enlarged for free (Prop PC). This is an instance of the doubly lister E3 and gives heuristic polynomial time, even
  with a promise constructor.

**Proof.**
1. u ≍ N is y-smooth with probability L[2/3, −2/(3κ)]. This is CEP with u = x/(κF_{1/3}) and log u = (2/3)log x.
2. v ≤ h = L[2/3, ν] is y-smooth with probability L[1/3, −ν/(3κ)], by the same computation one level down. This is
   HEURISTIC: the oracle's choice must not bias v.
3. So L[1/3, κ] relations need L[1/3, κ + ν/(3κ)] calls, and linear algebra costs L[1/3, 2κ].
4. At ν = 2/(3κ) the trial term κ + 2/(9κ²) is decreasing for κ < (4/9)^{1/3} = 0.763. Its intersection with 2κ is at
   κ = (2/9)^{1/3} = 0.606 < 0.763, so the max is minimised there: c = 2κ = (16/9)^{1/3}.
5. (iii): an interval must have log-length ≥ (1 − α)F_{1−α}/κ to contain a smooth u. The residue then has log-size of
   order F_{1−α} and is smooth with probability L[1 − 2α]^{−Θ(1)}.
6. (iv): relations are B-smooth on both sides by construction. One oracle call per window of t₀ consecutive k gives one
   relation, whatever t₀ is. Existence needs smooth u ≤ N·k_max in residue windows of width 2V, i.e. equidistribution
   of very smooth numbers modulo N (HEURISTIC, far beyond the Harper/Drappeau range). Fruitfulness and rank need the
   oracle's outputs to be random-like (HEURISTIC). An adversarial oracle that knows p, q could return only useless
   relations.
7. *Complete search gives nothing more.*
   - Intervals below the threshold (ν < 2/(3κ)) are almost all empty, which wastes L[2/3, 2/(3κ) − ν] calls per hit.
   - Listing *all* hits in one interval leaves each output's residue penalty L[1/3, −ν/(3κ)] unchanged.
   - Searching for the smooth number nearest kN shrinks v below the gap g only with probability w/g per interval,
     which again costs an L[2/3] factor. ∎

**Robustness to oracle cost.** If each oracle call costs L[1/3, τ], the trial term becomes κ + τ + 2/(9κ²). At ω = 2
this beats MNFS 1.9019 iff τ < 1.9019 − (3/2)(4/9)^{1/3} = 0.757 (at κ = (4/9)^{1/3} = 0.763; proof critic m4(d)).
So the comparison with MNFS does not depend on the oracle being exactly L^{o(1)}.

**Remarks.**
- The constant (16/9)^{1/3} coincides with Kirchner's special-n constant and with resbarrier's k = 1 rank-2 floor. The
  arithmetic reason: the relation's only random side has log-size (2/(3κ))F_{2/3}, which is half the base-m norm ledger.
- Prop SI is QS with the square replaced by a smooth number. Squares near kN have gaps ≍ √N, which forces v ≍ √N and
  L[1/2]. L[1/3]-smooth numbers have gaps L[2/3], which gives v ≍ L[2/3] and L[1/3].
- Codex's blind round derived the same optimum independently: b = 0.6057, s = 1.1006, c = 1.2114 (§6.1).
- Relation to Schnorr: u = ∏p^{e_p} ∈ I_k is a sparse inhomogeneous approximation Σe_p log p ∈ [log kN, log(kN + h)],
  with ‖e‖₁ ≤ x/log 2. That is a one-sided form of Schnorr's CVP target (schnorrcvp).
- smoothenum Prop W shows that balanced, height-bounded k-trees cannot do it (k ≤ 2, meet-in-the-middle, exponential).
  The pilist Johnson-radius class certifies one point per call in this congruence scaling, which is outside its NFS
  scope.

### 1.5 The σ–α table

`principle.sigma_alpha_table()`. Sampling rows use Theorem SS. Listing rows use Theorem LS.

| target | σ (λ ≍ x^σ) | λ allowed for constant c | candidates / factor base | size law needed (b = 0 / b = ½) | one-sided listing | doubly listing |
|---|---|---|---|---|---|---|
| α = 1/3 (reference) | 2/3 | (3c²/4)F_{2/3} | L[1/3, c] / L[1/3, c/2] | a = 2/3 / a = 1/2 (NFS) | Prop SI: 1.2114 from single intervals; Bohr sets any c | any c |
| α = 1/4 | 1/2 | c²F_{1/2} | L[1/4, c] / L[1/4, c/2] | λ ≍ √(x log x) / a = 3/8 | per-output L[1/4, τ], window L[1/2, ν], τ + √ν < c | box U/N = L[3/4] |
| α = 1/5 | 2/5 | (5c²/4)F_{2/5} | L[1/5, c] / … | x^{2/5}(log x)^{3/5} / a = 3/10 | τ + 2√(ν/5) < c | U/N = L[4/5] |
| α = 1/10 | 1/5 | (5c²/2)F_{1/5} | L[1/10, c] / … | x^{1/5}(log x)^{4/5} / a = 3/20 | τ + 2√(ν/10) < c | U/N = L[9/10] |
| quasi-polynomial | → 0 | (log x)^{O(1)} | exp(polylog x) | — | — | B = exp(polylog) |
| polynomial | → 0 | O((log x)²/log log x) | poly(x) | — | — | B = (log N)^K, U = N^{K′}, K′ > K/(K−1) (smoothenum T8) |

The one-sided listing column assumes the residue-window constraint √(αν) ≤ c/2 (pilist τ_max).

**What exists.**
- Only the α = 1/3 sampling row is realised (NFS).
- No family is known for α = 1/4: the aggregate rank-r family is absent under G_r.
- No family is known for α ≤ 1/5. At b = 0 it would need candidates whose total log-size ≍ x^{2/5} does not grow with
  the family. In ring-map terms that means norms far below the resultant floor for any degree that the linking allows
  (resbarrier C2, KL).

---

## 2. Sufficient conditions, one by one

Columns: exact statement · what it gives · extra assumptions · relation to a recognised problem
(⇐ implied by, ⇔ equivalent, ⇒ implies) · believed · status of the implication.

### Tier E: an exponent below 1/3, or polynomial time

**E1. A sampling family whose law gives C* = O(F_α), α < 1/3.**
- *Statement.* An L^{o(1)}-enumerable family whose law has C* ≤ cF_α, α < 1/3, with O(1) smoothness integers per
  candidate and a ring-map or congruence link.
- *Gives* L[α, c]. *Extra:* H_S, R.
- *Relation.*
  - ⇔ (within model S) any sub-1/3 sampling method (Theorem SS).
  - Instances:
    - E1a, the rank-r aggregate linked module: two degree-d ≍ √x fields, r ≍ x^{1/4} linked elements, family height
      N^{1/d}; formally 1/4 (resbarrier §4.2);
    - E1b, toric / entropy–height separation (heurnfs §5(i)).
  - Recognised neighbours: Schmidt's counts of algebraic points of bounded height, used as the heuristic G_r.
- *Believed* false for typical N. G_r predicts N^{e} families with e ≤ 1 − (r−1)(1 − 1/d) < 0, and none is
  constructed.
- *Implication:* model-PROVED.

**E2. One-sided Bohr-set lister.**
- *Statement.* Lists y-smooth u ≤ U, y = L[α, κ], with |u mod N| ≤ L[2α, ν], at L[α, τ] per output, τ + 2√(αν) < c.
- *Gives* L[α, c].
- *Extra:* H on the residue, R, and modulus-generic random-like outputs.
- *Relation.*
  - ⇔ Π-nonexact one-sided (pilist).
  - Contains the dense-PTE / systematic-equation escapes (systeq §6, ptesplit O1, ptecons). Their exact versions are
    refuted (ptebarrier Theorem W).
  - Smooth numbers in Bohr sets; one-sided Schnorr CVP.
  - ⇒ Prop DLT (below).
- *Believed* false: the box sieve costs L[1 − α] per output.
- *Implication:* model-PROVED.

**E3. Doubly-smooth lister.**
- *Statement.* Lists u ≡ v (mod N), both y-smooth (or both NFS norms smooth), at L[α, τ] per output.
- *Gives* any α (Theorem LS(a)), and heuristic polynomial time.
- *Extra:* R; existence (equidistribution of very smooth numbers mod N, HEURISTIC; or, unconditionally for some
  catalogues, nonexact Prop AB).
- *Equivalent forms* (PROVED, elementary):
  - (i) S-unit pairs u − v = kN with u, v ∈ Γ_S (abcsunit Prop SU). abcsunit Lemma FR sharpens the fruitfulness
    requirement: the relations must generate Λ_N up to odd index; near-full parity rank is not enough. abcsunit Prop EX
    proves existence (abundance) unconditionally at y = (log N)^κ. This is a correspondence of *witnesses*, not of
    algorithms:
    standard fixed-coefficient S-unit solvers do not enumerate over the free multiplier k at the required cost;
  - (ii) nonzero vectors of the relation lattice Λ_N = {e ∈ Z^{π(B)} : ∏p^{e_p} ≡ 1 (mod N)}, of polynomial ℓ₁-norm
    (Schnorr 1991; Adleman 1995). Any vector of Λ_N is a doubly-smooth relation, whatever its length. But the target is
    a *fruitful, full-rank* set of vectors of polynomial ℓ₁-norm: φ(N)·e_i are in Λ_N and useless. Approximate CVP in Schnorr's lattice is a *method* for producing them, not the target
    itself;
  - Precedent for (ii): Stange, arXiv 2211.06821, Thm 3.2. Random small multiplicative relations mod n over a factor
    base give the order, hence a factor (Monte Carlo, poly(b, log n)).
  - (iii) constructive navigation in Cay((Z/N)^*, {p ≤ B}). Words equal to 1 are E3 directly. The stronger *target*
    version is **PROVED sufficient with no heuristic** (codex suffmap_1 §3.3).
    - *Hypothesis:* known primes p ≤ B generate (Z/N)^*, and any target can be expressed as a word with
      polynomial-bit exponents in polynomial time.
    - *Proof:* navigate a_i² for hidden uniform units a_i and take a parity dependency among π(B) + 1 exponent
      vectors. This gives b with b² = ∏a_i². The hidden ∏a_i is uniform among the square roots of ∏a_i², so
      gcd(b − ∏a_i, N) splits N with probability ≥ 1/2.
    - *Conclusion:* randomized polynomial-time factoring. Navigation to arbitrary targets is DL mod N relative to the
      factor base, so this is KNOWN in substance (Bach 1984) and belongs with the one-way oracles of E7.
  - ⇒ Prop DLT.
- *Believed* false.
  - Generic lower bounds forbid computing orders or relations in black-box groups of unknown order cheaply (Babai 1999;
    Sutherland, MIT thesis 2007, §2.2: Ω(N^{1/3}) generic order computation; Damgård–Koprowski 2002 for root
    extraction; Shoup 1997 is for known order). So any method must use integer structure.
  - The only known integer structure is sieving.
- *Implication:* model-PROVED (smoothenum T7/T8).

**Prop DLT (listing transfers to prime-field DL) [PROVED reduction, modulo the standard index-calculus heuristics].**
- *Statement.* Let a lister for E2 or E3 use the modulus M only through arithmetic mod M (modulus-generic). Suppose it
  also handles cosets u ≡ t·v (mod M) for a given t.
- *Consequence.* Running it with M = p prime collects factor-base relations and individual-log relations for F_p^* at
  the same L-*exponent*, which gives L_p[α] discrete logarithms. The constant is not claimed: it depends on the linear
  algebra over Z/(p−1), the rank modulo every prime power ℓ^e ∥ p − 1, generator normalisation, and the cost of
  factoring p − 1 (codex suffmap_2 #2). Without the coset capability the individual-log step is not covered.
- *Proof.* Relations ∏p_i^{e_i} ≡ t^{ε}∏p_j^{f_j} (mod p) are linear equations in the logs over Z/(p−1). The linear
  algebra is polynomial in the factor-base size, which preserves α. Rank mod each ℓ | p − 1 is the usual heuristic. ∎
- *Reading.* Prime-field DL has been stuck at 1.902 (Matyukhin 2003; Commeine–Semaev 2006) (heurnfs §5(ii)). So any
  modulus-generic E2/E3 lister would break two problems believed hard. This parallels prepread's Corollary O1-N. A
  lister that exploits compositeness escapes Prop DLT, but no such mechanism is known.

**E4. Short straight-line programs for Pochhammer polynomials.**
- *Statement.* For n = 2^j ≤ √N, the polynomials P_n(X) = 2^{e(n)}∏_{a=1}^{n}(X − a) have constant-free SLPs of size
  (log n)^{O(1)}.
- *Gives* factoring in **P/poly**, rigorously and with no heuristic.
  - Store SLPs for P_{2^j}, j ≤ log₂√N, as poly(log N) advice.
  - Compute k! mod N for any k from the binary blocks, as ∏_j P_{2^j}(s_j) up to signs and powers of 2, which are units
    mod odd N.
  - Binary search on k with gcd(k! mod N, N).
  - If the SLPs are constructible in polynomial time from n, factoring is in P.
  - A *uniform* variant needs no SLPs at all: computing k! mod N in time poly(log N + log k) gives deterministic
    polynomial-time factoring by binary search (codex suffmap_1 §2.8). The integer-SLP question (τ) and the uniform
    modular question are distinct.
  - The weakest known sufficient statement for N-oblivious SLP methods is slpcirc's root-density condition W_U(s, μ)
    (slpcirc §1.3). It is necessary up to a factor 2 within that model, and it generalises E4.
- *Relation.*
  - ⇐ the permanent has polynomial-size constant-free circuits. That is Bürgisser 2009 (Comput. Complexity 18,
    81–103; ECCC TR06-113 Theorem 1.1, stated contrapositively; the hypothesis τ(Per_n) = n^{O(1)} is constant-free
    *and division-free*). The Pochhammer form is Thm 4.1(2): τ(2^{e(n)}∏(X − k)) = (log n)^{O(1)}, division-free and
    constant-free. "Allowing divisions" in the abstract refers only to the Taylor polynomials (Thm 4.1(3)). So the
    division-free requirement of E4 is met; this was checked by the literature critic in the ECCC text, with the journal
    numbering not checked. The n! form
    (τ(n!) polylog for all n) needs the uniformity or advice caveat below. The permanent hypothesis implies VP = VNP.
  - Contradicted by the Shub–Smale τ-conjecture. P_n has n integer roots, so the conjecture forces τ(P_n) ≥ n^{Ω(1)}.
  - Strassen 1976 and Lipton 1994 for the factoring link.
- *Believed* false, on the evidence of the τ-conjecture. VP ≠ VNP is *not* evidence: the implication runs Per easy ⇒
  Pochhammer easy, and VP ≠ VNP is consistent with short Pochhammer SLPs.
- *Implication:* KNOWN (folklore, slpcirc Prop BB).
- **Landmark either way:** a hardness proof would be a constant-free permanent lower bound (Bürgisser), and an easiness
  proof would be nonuniform polynomial-time factoring.

**E5. Complexity collapses.**
- *Statements and consequences.*
  - "N has a factor ≤ k" is in UP ∩ coUP (sorted factorisation + AKS certificate; polyimpl). So P = UP ∩ coUP, or
    P = NP ∩ coNP, puts factoring in P.
  - BQP ⊆ BPP puts it in BPP (Shor 1997). It is enough to simulate classically the specific order-finding circuits of
    Shor's algorithm, which is weaker than BQP ⊆ BPP. Classical order finding for a ↦ a^x mod N is itself
    randomized-equivalent to factoring (Miller 1976).
  - Jeřábek 2016 (JCSS 82, Thm 3.7) gives randomized reductions of factoring to PPA and to PWPP; Thm 3.8
    derandomises them under GRH_q (GRH for quadratic characters). So polynomial algorithms for PPA or PWPP give randomized polynomial-time factoring.
- *Believed* false. *Implication:* KNOWN, trivial. No structure.

**E6. Group-order families with biased smoothness.**
- *Statement (corrected after codex suffmap_2 #12).* An efficiently represented, efficiently samplable family of
  groups and points over Z/N, with L_p[α, β]-cost annihilation stages. For each p | N, with probability
  ≥ L_p[α, γ]^{−1} the annihilator kills the point's reduction mod p but not mod some other prime factor, which is
  detectable by a gcd.
- Smooth orders alone are not enough: the trivial group has order 1 everywhere and reveals nothing.
- *Gives* L_p[α] (ECM-style scalar multiplication), so L_N[α] for balanced N.
- *Barrier side.*
  - nonsmooth Prop OA.
  - Uniform random curves cannot be biased towards any set of orders by more than polylog factors (Lenstra 1987, Annals
    126, Prop. 1.9: class-number bounds for the number of curves with order in a given subset of the Hasse interval).
  - Torsion forced over Q is ≤ 16 (Mazur), and over degree-d fields it is bounded by Merel. Forcing torsion of size m
    needs points on X₁(m) over Z/N, of genus ≍ m².
  - CM families give ≤ 32 orders per p (ecmavg).
- *Gives* L_p[α, β + γ] by repetition.
- *Believed* false; no structure is known. *Implication:* PROVED, standard (Lenstra 1987 §2; codex suffmap_1 §2.10).

**E7. Circular equivalences and one-way oracles.**
- *Factoring-equivalent:* φ(N), λ(N), σ(N), square roots mod N, Frobenius operators on étale Z/N-algebras, additive
  power maps (nonsmooth §1; Prop FROB; systeq Prop F). These are also *necessary* in the trivial sense that fast
  factoring computes them.
- Also factoring-equivalent or circular: class-group or class-number computation for Q(√−N) (ambiguous forms;
  Schnorr–Lenstra 1984), #E(Z/N) point counting (Kunihiro–Koyama 1998; Dieulefait–Urroz), and modular square roots
  (Rabin 1979).
- *One-way sufficient only:* DL mod N and order finding (Bach 1984; Miller 1976). Factoring alone leaves prime-field DL,
  so these are not known to be equivalent (codex suffmap_2 #13).
- Either way they have no independent structural content.

**E8. Meta-complexity and learning (pointers from the literature critic; not checked by us).** These conditions mostly
concern Blum integers or the average case, and several give exponents far below 1/3:
- a natural property useful against TC⁰ ⇒ Blum integers factor in 2^{n^ε} for every ε (Razborov–Rudich 1997 with the
  Naor–Reingold PRFs, JACM 2004);
- MCSP ∈ P ⇒ fast factoring of Blum integers (Kabanets–Cai 2000); Factoring ∈ ZPP^{MCSP} (Allender et al., SICOMP 2006);
- efficient learning of TC⁰ / AC⁰-type classes under the uniform distribution with membership queries
  (Kharitonov 1993; Kearns–Valiant 1994);
- feasible interpolation for TC⁰-Frege (Bonet–Pitassi–Raz 2000), and Extended Frege versus RSA (Krajíček–Pudlák 1998);
- non-existence of one-way functions, or average-case easiness of K^t (Liu–Pass 2020) ⇒ random semiprimes factor on
  average.
All are believed false and give no structure; they belong with E5.

### Tier C: an L[1/3] constant only

**C1. Prop SI (single intervals near kN).**
- *Gives* 1.2114.
- *Statement and assumptions:* §1.4. The existence half is "every interval [X, X + L_X[2/3, ν]] contains an
  L_X[1/3, κ]-smooth integer for ν > 2/(3κ)". That is a smooth-numbers-in-all-short-intervals statement at a scale far
  beyond current results. All-interval results: Hildebrand 1986; Friedlander–Granville 1993 (y ≥ exp((log x)^{5/6+o(1)}));
  Croot 2007; Matomäki 2016. Almost-all results: Hildebrand–Tenenbaum 1993; Jain, arXiv 2502.10530
  (y ≥ exp((log X)^{2/3+ε})). Nothing reaches y = L[1/3].
- *Recognised neighbour:* deterministic construction of smooth numbers in a given interval. Boneh (STOC 2000 / JCSS 64
  (2002), Thm 3.1) works only in the polylog-y regime, log x < y < (log² x)/5 per Granville 2008 §5.7, on windows of
  length x^{(log x)/(4y)}. At Prop SI's scale those windows are shorter than 1, so the method is vacuous there.
  Boneh's §4 sketches a QS use, states no reduction and no constant, and notes the bounds are insufficient for QS.
- *Novelty (literature critic).* No prior statement of the reduction or of 1.2114 was found. It is, however, an
  immediate instance of Pomerance's principle with X = h = L[2/3] under the constraint that intervals must contain
  smooth numbers. Its closest structural analogues are Schroeppel's linear sieve and the cubic sieve
  (Coppersmith–Odlyzko–Schroeppel 1986), where one side is constructed and the residue is N^{1/2} or N^{1/3}.
- *Believed:* open; no method below exhaustive sieving per output is known. Prop W and the pilist Johnson-radius class
  fail here.
- *Implication:* model-PROVED. No prior statement found by us.

**C2. Complete search or sampler for smooth values on NFS rational rows.**
- *Parameters:* X = L[2/3], h = L[1/3], y = L[1/3]; in X-units h, y = L_X[1/2].
- *Gives* 1.6386 (smoothenum T3, KNOWN as the factory per-N program), via Lemma SL. A promise constructor gives the same
  if h₀ = gap·L^{o(1)} (Prop PC). A partial version gives c(θ) with slope −0.21.
- *Dominated by C1 in constant*, but it has different parameters.
- Practical precedent: Hittmeir, "Smooth Subsum Search" (arXiv 2301.10529), which constructs candidates divisible by
  many factor-base primes; a constant-factor heuristic.

**C3. Algebraic-side rows (smooth values of F(·, b)).**
- *Gives* inf 0.7631 (smoothenum T6; 0.8246 at the e ≤ 3 cap).
- Doubly row search gives inf 0 at L[1/3]. Growing boxes then give E3.
- **Caveat (codex suffmap_2 #3).** Row-by-row search pays the 2K + J startup term: E·V row calls, i.e. L[1/3, e + v].
  The T6 infimum is approached with e → ∞, so a row-wise implementation does not reach it. Codex reports 1.16477 for
  the program with the e + v term added; we have not checked this. Reaching 0.7631 or 0 needs a whole-box
  two-dimensional oracle that can skip empty rows. C2 is unaffected, because E ≪ E²L[−r] at its optimum.

**C4. Beyond-Johnson, output-sensitive CRT list decoding of the factor-base code; H_sel violations.**
- *Gives* c(θ′) via pilist Remark SEL (c). Within the Johnson radius the gain is L^{o(1)} (Prop JR).
- *Relation:* CRT list decoding (Goldreich–Ron–Sudan; Boneh; Guruswami–Sahai–Sudan). Plain-RS capacity decoding is a
  2026 preprint (Brakensiek et al.), with no CRT analogue at vanishing rate.
- *Believed* implausible: pilist puts it below 1%.

**C5. Sparse F₂ linear algebra in K^{ω+o(1)}, ω < 2, for K^{1+o(1)} nonzeros.**
- *Gives* MNFS 1.8595 at ω ≲ 1.67 (heurnfs B4). By Theorem SS(c) no single-family method moves.
- *Relation:* "solve sparse linear systems over finite fields faster than Wiedemann's n·nnz" (Wiedemann 1986; block
  Wiedemann, Coppersmith 1994; Eberly–Giesbrecht–Giorgi–Storjohann–Villard 2006/07). Peng–Vempala (SODA 2021) give
  O(n^{2.33}) for poly-conditioned real systems with O(n) nonzeros. That beats matrix multiplication but is *worse*
  than n·nnz ≈ n², so it says nothing about ω < 2.
- *Believed* open; the least implausible item on the map.

**C6. Small-height polynomial pairs for general N.**
- *Statement:* heights with κ < 2 in heurnfs's ledger.
- *Gives* (32κ/9)^{1/3}, never an exponent. resbarrier C2/RB2 gives α ≥ 1/3 for every rank-2 pair and every N *within
  its model*: §2.2-form extraction, D ≥ 3, box or lattice side ≤ B^{1−o(1)} (H5), heuristic H, B² linear algebra, single
  N.
- *Status:* **PROVED false for almost all N** at the height level (resbarrier GAE). The discriminant version is open.
- **Small-discriminant number fields with a root mod N** are this item. N² | disc is forced for two-map extraction and
  N | disc for ramified extraction (Lemma ND). Bounded root discriminant (Golod–Shafarevich towers) gives small-norm
  ideals but no linked rank-2 family.
- Small discriminant plus a ring map is not even a size certificate on its own (codex). Q has discriminant 1, and its
  ring map Z → Z/N gives candidates a − bm as large as m. What matters is the height or discriminant of the *linked
  pair*.

**C7. Kirchner 2018 fields.**
- *Gives* 1.2114 for very special n (k = 1). Excluded by GAE-ram for almost all N.
- The relaxed version gives 1.5263 "for most n". It is not excluded by GAE-ram, but needs a Galois-subfield structure and
  a generator that nobody has.
- Kirchner's abstract: (16/9)^{1/3} for special n; for general numbers it "may reach" (32/9)^{1/3} "if the method can
  be adapted"; "there is no evidence it can be done". The full-text phrases quoted in resbarrier §4.1 were not
  re-checked by our literature critic.
- *Implication:* conditional claim, not verified by us (resbarrier §4.1).

**C8. Non-uniform precomputation (the factorization factory).** 1.6386 per N after L[1/3, 2.0069], KNOWN (Coppersmith
1993; Bernstein–Lange 2014). Not single-N RAM.

### Tier X: the exponential regime

**X1. Coppersmith / partial key exposure.**
- *Statement:* poly-time factoring from a fraction f < 1/2 of p's top bits.
- *Gives* N^{f/2}. It beats the deterministic N^{1/5} record (Harvey, Math. Comp. 90, 2021; log-log speedup by
  Harvey–Hittmeir, Math. Comp. 91, 2022) iff f < 2/5. The topblock ladder needs f = 2/(j+2).
- For L[α] one would need f·x/2 = O(F_α). That is factoring itself up to guessing: circular (topblock).

**X2. Coverage.** dchworst Theorem R: c > 3/8 gives a worst case below N^{1/5}; DCH_pair gives N^{1/6}.

**X3. Polylog evaluation of degree-2 floor sums.** nonsmooth §6, floorsum: exponential-regime ladders only.

---

## 3. High-profile candidates: sufficient, necessary, or unrelated

| # | candidate | verdict | precise reason |
|---|---|---|---|
| H1 | Matomäki–Radziwiłł 2016 (Cor. 1: X^ε-smooth in [X, X + C√X]; Cor. 6: almost all intervals); Hildebrand, Friedlander–Granville, Croot, Matomäki, Jain 2025 (small y). Harper's papers are on APs and minor arcs, not short intervals | **irrelevant to speed**: yields no algorithm and is used by none; relevant to *rigour* | Existence theorems at y = X^ε carry no algorithm. Factoring needs y = X^{o(1)}. Their L-scale analogues are the *existence half* of C1/C2's promise routes (Prop PC) and inputs to rigorous relation supply (the project's Theorem U uses Drappeau/Harper for smooth numbers in APs; Lee–Venkatesan 2018 for rigorous NFS) |
| H2 | "every [x, x + x^ε] contains an x^ε-smooth number" | **unrelated** (stated in Granville 2008 §1.2/§4.4, who shows it implies Vinogradov's least-nonresidue conjecture for p ≡ 3 mod 4) | wrong scale: positive density ρ(1/ε), y = x^ε; no construction |
| H3 | deterministic construction of a y-smooth integer in a given interval | **sufficient**: 1.2114 at (X ≈ kN, h = L_X[2/3, ν], y = L_X[1/3, κ], ν > 2/(3κ)) (Prop SI); 1.6386 at NFS rational rows if complete or with h₀ ≈ gap (Lemma SL, Prop PC); **heuristic poly time** in the Bohr-set form. **Insufficient** for exponent < 1/3 in single-interval form | §1.3–1.4 |
| H4 | Polymath4 (deterministic primes in [x, 2x]; Tao–Croot–Helfgott 2012) | **unrelated** in the sufficient direction | primes are dense (1/log x) and need no smoothness; the smooth analogue is H3 |
| H5 | sum-product (Erdős–Szemerédi; Bourgain–Katz–Tao; Chang) | **unrelated** | These are *upper* bounds on additive structure of multiplicative sets. They are consistent with H and too weak to bound box counts below H's prediction. Inverse theorems (PFR) concern small doubling, which box values lack |
| H6 | expansion of Cay((Z/N)^*, small primes) (GRH) | existence of short paths, diameter and mixing **unrelated**, since they yield no algorithm. (Abundance of collisions, nonexact Prop AB, is a different statement: collisions can lie in a proper subgroup.) **Constructive navigation sufficient**. It is DL mod N relative to a factor base, so it is KNOWN in substance (Bach 1984: DL mod N ⇒ factoring); GRH expansion is Jao–Miller–Venkatesan 2009: randomized polynomial time, PROVED with no heuristic (E3(iii), codex) | generic-group bounds make navigation without integer structure exponential |
| H7 | number fields of small discriminant with a prescribed root mod N (Kirchner; resbarrier) | **sufficient for constants only** (≤ 1.5263 two-map, 1.2114 ramified); never an exponent (C2/RB2) | heights refuted for almost all N (GAE, GAE-ram); discriminant version open |
| H8 | abc; S-unit counts (Evertse 1984; Erdős–Stewart–Tijdeman 1988; Lagarias–Soundararajan's xyz conjecture) | **unrelated** | With S = primes ≤ B, rad(uv) ≤ e^{B(1+o(1))}, so abc is vacuous for u − v = kN. Evertse's 3·7^{3+2s} (s finite primes; his count includes the infinite place) is exp(O(B)). Existence is not the bottleneck (S²/N − S pairs, nonexact Prop AB). *Algorithmic* S-unit solving is **sufficient** (E3). abcsunit Props ABCV/EV agree: abc is automatically satisfied unless the relation is unusually powerful, and Evertse counts only the triply-smooth subfamily |
| H9 | approximate CVP in Schnorr–Adleman prime lattices | a **method** for E2/E3; sufficiency for any approximation factor **not established** | resbarrier (h): |uw| ≥ N − 1 and the lattice-shortness step fails; Ducas 2021 experiments. Details: schnorrcvp |
| H10 | Coppersmith / partial-key bounds | **sufficient in the exponential regime only** (X1) | L[α] needs f = O(F_α/x), circular |
| H11 | CRT list decoding beyond Johnson | **sufficient for constants** if output-sensitive (C4); a doubly version would be E3 | pilist Prop JR, Remark SEL |
| H12 | sub-quadratic sparse linear algebra | **sufficient for a constant** (C5) | Theorem SS(c): single-family sampling ignores ω ≤ 2 |
| H13 | factorial / Pochhammer SLPs; τ-conjecture; VP vs VNP | short Pochhammer SLPs (suitable 2-power multiples) are **sufficient** (E4), nonuniform and rigorous. Constant-free, division-free poly-size permanent circuits imply them. The τ-conjecture, if true, *refutes* this route. VP ≠ VNP says nothing about it | Strassen 1976, Lipton 1994, Bürgisser 2009 |
| H14 | GRH, RH | **unrelated to speed** | They buy determinism and rigour (Miller; detGRH), not exponents (codex r11: [QRH] gives no gain) |
| H15 | P = NP ∩ coNP; BQP ⊆ BPP | **sufficient**, trivially (E5) | |
| H16 | RSA inversion vs factoring (Boneh–Venkatesan 1998: no algebraic reduction from factoring to low-degree RSA inversion) | **out of scope** | breaking RSA might be easier than factoring; nothing on this map needs or gives that |

**What is *necessary*?**
- No *independent structural hypothesis* on the map is known to be necessary for faster factoring. (The
  factoring-equivalent invariants of E7 are necessary trivially.) For example, factoring in P is not known to give a
  doubly lister in the useful sense: a *fruitful, full-rank set of polynomial-ℓ₁ relations*. Given p and q the trivial
  vectors φ(N)·e_i lie in Λ_N but are useless, and polynomial-size ones seem to need DL mod p.
- The necessary statements available are model-internal, and they are the barrier literature: Theorem SS(c), resbarrier
  C2/ND/GAE, LQ/JR/TAIL.
- Theorem SS is an *equivalence within model S*: an L[α] sampling method exists iff a family with C* ≍ F_α exists. MNFS,
  unions of pairs and the factory are outside S (§1.1).

---

## 4. Ranking: directness × plausibility × profile

Scores are 1–5 (assessment). Directness: 5 = immediate exponent or polynomial time with no heuristic; 1 = exponential
regime only. Plausibility: our credence that the hypothesis is true and attackable. Profile: recognition outside this
project.

| target | D | P | R | product | note |
|---|---|---|---|---|---|
| C1 Prop SI (smooth numbers in L[2/3]-intervals near kN) | 3 | 2 | 3 | 18 | big constant (1.90 → 1.21); concrete. An oracle on *all* intervals [X, X + h] would also give prime-field DL below 1.902 (relations as in Prop SI with N → p; individual logs by a descent through intervals [kp/q, (kp + h)/q], about L[1/3, 0.83] per proof critic m5; not checked in detail). The I_k-only oracle of §1.4 does not support Prop DLT's coset hypothesis |
| C5 sparse F₂ LA, ω < 2 | 2 | 3 | 4 | 24 | small constant; a mainstream TCS problem |
| E4 Pochhammer / factorial SLPs | 5 | 1 | 5 | 25 | **landmark either way** |
| E3 doubly / Bohr-set lister (Schnorr target, S-units, navigation) | 4 | 1 | 4 | 16 | would also break prime-field DL (Prop DLT) |
| C2 NFS rational-row search | 2 | 2 | 2 | 8 | dominated by C1 |
| C4 beyond-Johnson CRT decoding | 2 | 1 | 3 | 6 | pilist < 1% |
| E1a aggregate rank-r | 3 | 1 | 1 | 3 | absent under G_r |
| C7 Kirchner fields | 2 | 1 | 2 | 4 | special case refuted for a.a. N |
| X1 Coppersmith below ½ of the bits | 1 | 2 | 4 | 8 | exponential only |
| E5 collapses | 5 | 1 | 5 | 25 | no structure; not a research target |
| C6 small-height pairs, general N | 2 | 0 | 2 | 0 | PROVED false for a.a. N |

**Most worth attacking: C1 (Prop SI).**
- It is a single, clean question: construct L[1/3]-smooth numbers in intervals of length L[2/3] near multiples of N.
- It sits at an untested scale (between Boneh's CRT decoding and Schnorr's CVP).
- It has no barrier beyond method classes (Prop W, LQ, JR in other scalings).
- Its payoff would be the largest known constant drop, 1.90 → 1.21.
- A natural first test: is it equivalent to the NFS-scale smooth-interval problem C2, or to one-sided Schnorr CVP at
  dimension L[1/3]?
- E4 and E5 score higher on the raw product, but they are not attackable research targets in this project.
- This is an assessment. The payoff rests on extra oracle-distribution assumptions (random residues, fruitfulness).
  Codex's adversarial round prefers C2 and C5 for their clearer interfaces (§6.2).

**Landmark either way: E4.** Both outcomes are major theorems: nonuniform polynomial-time factoring, or a constant-free
permanent lower bound. The hardness must be quantified over the allowed power-of-two multiples.
- This is not literally unique. Any unconditional superpolynomial lower bound for an NP search problem, e.g. for the
  Prop SI constructor, would imply P ≠ NP (codex suffmap_2 #17).
- E4 is distinguished because its negative side is a natural, actively studied lower-bound problem (Bürgisser's
  route to VP⁰ ≠ VNP⁰), rather than P ≠ NP itself.

**Highest-profile exponent target: E3**, together with its high-profile faces: S-unit enumeration, Schnorr's lattice,
Cayley-graph navigation. Prop DLT makes it implausible.

---

## 5. Feasibility checks (exact; `src/suffmap/test_suffmap.py`, 12 tests)

1. **Theorem SS(b).** `optimum_fixed_size` returns b* = √(αs) and c = 2√(αs) exactly, with trials = LA at the optimum,
   for (α, s) ∈ {(½, 1), (½, ½), (¼, 1), (1/5, 5/4), (1/10, 3)}. So Dixon √2 and QS 1 hold, and α = 1/4 with s = 1
   gives c = 1.
2. **Theorem SS(d).** The law table gives exactly:
   - Dixon α = ½, c = √2; QS c = 1;
   - GNFS c³ = 64/9; SNFS c³ = 32/9; t = 3 c³ = 256/27 (= heurnfs's (128(t−1)/(9t))^{1/3});
   - aggregate α = ¼.
   The symbolic ratio C²/(2λ log λ) is exactly 1 for every law. That check is partly circular, because c is computed
   from the same identity, so it confirms the substitution algebra (log powers and exponents) only. The independent check
   is the agreement with the published constants.
3. **σ table.** σ = 2/3, 1/2, 2/5, 1/5 for α = 1/3, 1/4, 1/5, 1/10. The exponent rule α = a/(2−b) is checked on
   (½, ½) → 1/3, (2/5, 0) → 1/5, (3/10, ½) → 1/5 and (1/5, 0) → 1/10.
4. **Prop SI.** κ³ = 2/9 and c³ = 16/9 exactly. The intersection lies left of the trial-term minimiser. The
   single-interval exponent max(α, 1 − 2α) is minimised at α = 1/3 on a grid of 59 rationals.
5. **Lemma SL.**
   - `python -m src.suffmap.listing`: exact smoothness oracles with three different hit-choice rules list all 33 smooth
     integers of [10⁶, 10⁶ + 2000] (y = 50) in 66 ≤ 67 calls. On toy NFS rational rows (m = 10⁵ + 3, E = 200, y = 60,
     20 rows) the 234 outputs take 485 ≤ 488 calls, against 8020 box points.
   - The tests use [10⁶, 10⁶ + 1500] with y = 40, the sparse case (y = 7 near 2^31), and rows with E = 150, y = 50,
     12 rows.
6. **Prop PC.** With the adversarial promise oracle, tiling plus sliding queries never exceed 4H/h₀ + 1 distinct
   outputs, and revealed points are pairwise ≥ h₀/4 apart. An honest constructor tiled gives ≥ ⌊H/h₀⌋.
7. **Integer congruences.** max(|u|, |v|) ≥ N/2 for u ≡ v, u ≠ v, checked exhaustively for N ∈ {15, 21, 35, 77}.

These are formula checks at toy scale. They confirm algebra, not asymptotic behaviour.

---

## 6. GPT-6.1-sol rounds

### 6.1 Blind round (`codex/suffmap_1.out`, run before our analysis)

**Independent agreements.**
- The fixed-size optimum exp((√2 + o(1))√(λ log λ)), and c = 2√(αs) in L-units.
- The size–entropy equation D² ∼ 2λ(D) log λ(D), with α = a/(2−b), t^{2−b} = 4ακ, and the log-power condition
  c = 1 − a − b. It reproduces GNFS 64/9 and SNFS 32/9 from the same ledgers.
- σ = 2α, with the warning that λ ≍ x^σ without the log factor gives exp(Θ(x^{σ/2}√log x)).
- The polynomial-time boundary λ log λ = O((log x)²).
- ω ≤ 2 does not move fixed-ledger sampling.
- The search ⇒ listing bisection, and the limitation of long-interval promises.
- The single-interval constant **(16/9)^{1/3} = 1.211414**. Codex lists all smooth u in one interval
  [N + 1, N + H] instead of constructing one u per interval [kN, kN + h]. Both give the same program.
- 1.6386 for a perfect rational lister.
- Height profiles (32(1 + η)/9)^{1/3}, the almost-all counting floor, and Kirchner's claims (conditional, unconstructed).
- Bürgisser: τ(Per) poly ⇒ τ(n!) polylog, "landmark either way".
- MR and Harper: existence or average statements only.
- Polymath4: no factoring implication.
- RS capacity: no CRT transfer.
- Peng–Vempala: not an F₂ sparse solver.

**Adopted from codex.**
- The PROVED randomized reduction from constructive navigation (§2, E3(iii)). The literature critic notes it is KNOWN in
  substance (Bach 1984).
- The remark that small discriminant plus a ring map is not a size certificate.
- Jeřábek's PPA/PWPP reductions.
- The Balister–Bollobás–Morris citation, and the point that H alone does not give the relation-count lower bound
  (assumption R made explicit).
- The impossibility of an every-interval constructor at α < 1/3 (density e^{−Θ(F_{1−α})} vs e^{O(F_{2α})}), which is
  Prop SI(iii) stated as a counting fact.
- The multi-parameter ledger t³ = 64q/(9k), with k = x^θ giving exponent (1 − θ)/3.

**Ours.**
- Prop PC, the quantified promise-constructor value. It shows that codex's 1.2114 needs only a *promise* constructor on
  many intervals near kN, provided all such intervals contain smooth numbers.
- The Bohr-set polynomial-time version.
- Prop DLT.
- The Pochhammer form of E4, which needs only poly-size advice. A nonuniform τ(k!) bound for each k does not by itself
  give poly-size advice; codex notes the uniformity gap.
- Tier assignments, the necessity column, and the ranking.

**Disagreement.**
- Codex ranks the NFS rational-side lister (1.6386) as the best *initial* target and the interval search (1.2114)
  second.
- We rank Prop SI first, because it is the same kind of problem at a different scale with a larger payoff, and it needs
  no number field.
- Both rankings are assessments.

### 6.2 Adversarial round (`codex/suffmap_2.out`): MAJOR REVISION, all points adjudicated

Codex confirmed the fixed-size optimum, the ω-claim, the power-law algebra, the 1.2114 computation (no
double-counting), Lemma SL, Prop PC's packing and tiling arguments, the constructive-navigation proof and E4's core
implication. Its findings, many of which overlap the proof critic's:

| # | Sev. | Finding | Adjudication |
|---|---|---|---|
| 1 | MAJOR | Prop SI quantifiers: fixed (κ, ν) gives its own cost; R needed | accepted: tunable-constructor wording, "1.2114 + ε" as the optimised conditional constant, R explicit |
| 2 | FATAL (to SI ⇒ DL) | SI supplies no cosets; DLT preserves α, not the constant | accepted: the C1 DL remark is now conditional on an all-interval oracle plus an unchecked descent; DLT claims the exponent only |
| 3 | MAJOR | C3: row-wise search pays L[e + v] startup; the 0.7631 infimum needs e → ∞ | accepted: caveat added, 2-D oracle required; codex's 1.16477 not checked by us |
| 4 | MAJOR | SS(c) generality not justified without stability and population beyond threshold | accepted: stability hypothesis (C*_δ), population and concentration; restricted theorem for regularly varying laws |
| 5 | MAJOR | BBM/CGPT do not establish R | accepted: cited for their actual model; R an assumption for structured families |
| 6 | MINOR | domain 0 < α < 1 (a < 2 − b); polynomial endpoint costs | accepted |
| 7 | MAJOR | fixed-α listing does not give polynomial time; T8 needs prefix cost | accepted: LS(a) rewritten |
| 8 | MAJOR | "enlarged for free; any h₀ works" needs encoding bounds and a blockwise promise | accepted: Prop PC and SI(iv) |
| 9 | MINOR | PC bounds oracle answers; θ readings are count substitutions | accepted |
| 10 | MINOR | SI existence: mean density ≠ uniform existence; every k is false (CRT) | accepted: k ≤ L[1/3, C]; promise is a hypothesis |
| 11 | MAJOR (tables) | algebraic-complexity arrows; τ-conjecture is not "sufficient"; 2-power scalings | accepted: H13 and E4 rewritten |
| 12 | FATAL (to E6) | smooth orders do not give CRT separation (trivial group) | accepted: E6 restated with a separation event (codex suffmap_1 §2.10 form) |
| 13 | MAJOR | DL mod N is one-way, not equivalent; necessity paragraph | accepted: E7 split; "no independent structural hypothesis is necessary" |
| 14 | MAJOR | C6/H7 dropped RB2's scope | accepted: scope restored |
| 15 | MAJOR | abundance ≠ diameter; S-unit equivalence is witness-level | accepted: H6 and E3(i) |
| 16 | MINOR | Boneh already links interval search to factoring; Jain (arXiv 2502.10530) for small-y short intervals | accepted; novelty of the exact 1.2114 optimisation remains "not located" (§7.2) |
| 17 | MAJOR | the ranking is not supported as unique; missing items (uniform modular factorial, slpcirc W_U, Shor-specific simulation, RSA vs factoring) | partly accepted: the ranking is labelled an assessment, codex's preference for C2/C5 is recorded, the "landmark" claim is qualified, and the missing items are added (E4, E5, §3 H16). We keep Prop SI first: same interface type as C2, larger payoff, no number field. |

## 7. Opus critics

### 7.1 Proof critic: ACCEPT WITH FIXES, all applied

The critic ran the tests and both `__main__` blocks and confirmed every constant: √2, 1, 64/9, 32/9, 256/27, the σ table,
κ³ = 2/9 and 1.2114. Its findings:

| # | Sev. | Finding | Fix |
|---|---|---|---|
| M1 | MAJOR | MNFS is a sampling method outside model S (shared filtering, more columns); H_S's selection clause strengthens heurnfs H; unions of pairs have a better law | §1.1 defines S with one factor base and full testing; lists what S excludes; §0 scope bullet; §3 equivalence scoped |
| m1 | MINOR | the C* lower bound needs a stability/regularity hypothesis (step-profile counterexample) | C*_δ and stability added to §1.2(c); monotonicity step written out |
| m2 | MINOR | the polynomial regime needs poly(x) costs and loses O(log x) | (e) qualified; `polytime_size_bound` string fixed |
| m3 | MINOR | the congruence existence condition is on the box area (U/N)·V | §1.3(c) corrected; the α = 1/3, V = L[2/3] exception noted |
| m4 | MINOR | Prop SI: R needed in (ii); "complete search gives nothing more" unproved; 1.2114 is an infimum; oracle-cost robustness; (iv) is a separate Bohr-set hypothesis (= E3) | step 7 added; "+ε"; robustness threshold τ < 0.757; (iv) relabelled; §0 and table fixed |
| m5 | MINOR | the ranking's DL claim for C1 overreached Prop DLT (coset hypothesis; oracle domain) | weakened; all-interval oracle plus descent (critic's sketch, ≈ L[1/3, 0.83]) stated as unchecked |
| m6 | MINOR | E4: division-free form overstated; VP ≠ VNP is not evidence against short Pochhammer SLPs | UNVERIFIED tag; evidence is the τ-conjecture only |
| m7 | MINOR | φ(N)·e_i ∈ Λ_N contradicts "factoring in P gives no doubly lister" | restricted to fruitful, full-rank, polynomial-ℓ₁ relation sets |
| m8 | MINOR | Prop PC's mean-gap reading needs a max-gap heuristic; "(1+η)r" is thinning, not θ | labelled HEURISTIC; wording fixed |
| n1–n11 | NIT | Lemma SL call count; R in the lower bound; the a/(2−b) = 1/3 line; distinct tuples for pigeonhole; "beating a law"; H1 wording; §5 sources; circularity of `verify_law_balance`; §6.1 | applied |

### 7.2 Literature critic: fixes applied

The critic read the primary texts: Bürgisser ECCC TR06-113, Jeřábek, the Lenstra 1987 scan, Boneh JCSS, the Granville
and Pomerance chapters of MSRI 44, Polymath4, Stange, Sutherland's thesis, and MR.

**Confirmed.**
- Pomerance 1982/2008.
- MR Corollaries 1 and 6.
- Bürgisser Thm 1.1 and Thm 4.1(2), division-free and constant-free. This resolves E4's UNVERIFIED tag; the
  journal-version numbering is still not checked.
- Strassen, BCSS, Lipton and Shub–Smale as cited by Bürgisser; the τ-conjecture is Smale's 4th problem.
- Lenstra 1987 **Prop. (1.9)**, both the upper bound used here and the lower bound for the inner half of the Hasse
  interval. Codex's "Prop. 2.6" is the single-curve success condition.
- Jeřábek Thm 3.7 and 3.8 (GRH_q).
- Evertse, Erdős–Stewart–Tijdeman, Lagarias–Soundararajan, Evertse–Schlickewei–Schmidt.
- BBM, CGPT, Coppersmith 1993, Kirchner's abstract, Brakensiek et al.
- Polymath4: it uses no factoring oracle, so the verdict "unrelated" stands.

**Fixes applied.**
- Pomerance 2008 states the size principle explicitly, so Theorem SS(a) and σ = 2α are KNOWN; the ICM 1994 lower bound
  is added.
- Boneh is vacuous at Prop SI's scale. The novelty of Prop SI is modest: it is an instance of Pomerance's principle,
  with the linear and cubic sieves as analogues.
- Harper is misattributed for short intervals; the correct short-interval literature is cited instead.
- Granville 2008 is cited for H2.
- The generic-group citations are corrected (Babai, Sutherland).
- Peng–Vempala is reframed: n^{2.33} is worse than n·nnz.
- Kirchner quotes are taken from the abstract only.
- The N^{1/5} record is attributed to Harvey 2021.
- The Evertse convention is stated.

**Missed prior art, added.**
- Stange 2022 (E3(ii)).
- Bach 1984 for navigation (E3(iii) is KNOWN in substance).
- Hittmeir's Smooth Subsum Search (C2).
- Lee–Venkatesan (rigour).
- Jao–Miller–Venkatesan (H6).
- Cubic and linear sieves (Prop SI).

**Missing items, added.**
- The meta-complexity and learning conditions (E8, as pointers, unchecked).
- Class groups, #E(Z/N) and square roots in E7.

**Not added.**
- A separate "no Frobenius over Z" entry: systeq Prop F and nonsmooth Prop FROB already cover it, and E7 cites them.
- Fellows–Koblitz: not cited, since the AKS certificate argument suffices.

## 8. Reproduction

```bash
.venv/bin/python -m pytest -q src/suffmap/test_suffmap.py     # 12 exact checks, < 1 s
.venv/bin/python -m src.suffmap.principle                      # law table and sigma-alpha table
.venv/bin/python -m src.suffmap.listing                        # Lemma SL / toy rational rows
```
