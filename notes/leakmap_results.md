# leakmap: where "information costs splitting" stops (phase 24)

Agent `leakmap`, 2026-10-23. Version v1 (before the codex adversarial round and the Opus critics; see §11).
- Pre-registration: `notes/leakmap_prereg.md` (commit `d050db1`, before any code or experiment).
- Code: `src/leakmap/{netrees,e1,e1_p5,cfgenus,e2_posthoc,jacobi}.py`. Tests: `tests/test_leakmap.py` (20 pass).
- Raw data: `experiments/raw/leakmap_e1*.jsonl`, `leakmap_e2_{main,fresh}_n{32,40}.jsonl.gz`, `leakmap_e3.jsonl`.
  Summaries: `experiments/analysis/leakmap_*.json`.
- GPT-6.1-sol blind round: `codex/leakmap_1.{md,out}` (several results below are codex's and are credited as such).
- Labels: **PROVED** (proof here; "no exact prior formulation located" is the strongest novelty claim),
  **KNOWN** (reference), **HEURISTIC**, **CONJECTURED**, **OPEN**, **experimentally supported**, **refuted for this
  model**, **UNVERIFIED**.

Notation. N = pq, p < q balanced primes; for theorems, p, q iid with law μ on the primes of an interval I = (X, 2X]
(the virtual-diagonal convention of factorbits §3.2 handles p = q). π(I) = number of primes in I. A *tree* is a
finite adaptive decision tree whose nodes compute an element of (Z/N)-algebra and test it for zero, reporting a
split when a component-separating zero gives a proper gcd. s = Pr[split]. B = max_b Pr[F = b] for the target F.

---

## 0. Verdict (plain language)

**No fourth leak family was found, and for four of the five natural extensions of the pure-tree model we can now
say exactly what leaks.** The map is:

1. **N as an exponent (class G) leaks only through group-order coincidences, which degenerate into archimedean
   pins.** For trees that may raise pure values to exponents E(N) = αN + β or (N − c)/2^j and test the result
   (Fermat, Miller–Rabin, Solovay–Strassen's power, Lucas, Frobenius, terminal AKS-type element tests):
   - every bit of information obtained without splitting comes from a *deviation event*: the gate output
     y^{E(N)} hits a root of the node polynomial at p although the test is formally nonzero (Theorem PT-X);
   - such an event needs either a pure base of small order at p (a pure group-order event, family M), or the
     *cross coincidence* ord_p(y) | E(N) − t, i.e. gcd(E(N) − t, |G_p|) ≥ ord_p(y) (family M, the conjectured
     gcd(N−1, p−1) structure), whose extreme case exp(G_p) | E(N) − t is a Diophantine pin such as q = p + 2 or
     q = p + k − 1 (family A) (Theorem NE1);
   - quantitatively, for iid balanced primes the cross-coincidence term is ≤ h·D·(|I|/L + 1)/π(I) per tree of height h,
     i.e. N^{−1/4+o(1)} with L = √X; with fresh uniform bases it is exactly the product
     gcd(E(N), p−1)gcd(E(N), q−1)/φ(N) form (codex; order-pattern lemma).
   - E1 checked every inequality on 1.9·10⁶ instrumented runs (0 violations) and found the archimedean pins
     exactly where predicted (e.g. twin primes q = p + 2 with (5|p) = −1, (5|q) = +1 are non-split Lucas liars).
   - This needs the exponent menu to be fixed: codex shows that an arbitrary computable exponent E(N) = 1 + h(N) on
     the base −1 reads any polynomial-time predicate h without splitting.
2. **Integer representatives are computationally complete, so no PT-type theorem can exist there.** Reading
   rep(−1) = N − 1, or even only the low bits of 2^j mod N (codex), gives N, and then any function of N. All three
   families live in this model (QS/NFS relation tests, Lehman, SQUFOF). What survives is a *structure* statement:
   every representative of a residue is rep(x)/N ≡ (x_p − x_q)θ + x_q/N (mod 1) with one hidden phase
   θ = (q^{−1} mod p)/p ≈ j/q, so all cross-prime information in one-shot readouts is Diophantine information about
   θ; 3n/4 bits of θ factor N (Howgrave-Graham's partial approximate common divisor), the same N^{1/4} barrier as
   Coppersmith's known bits.
3. **Continued fractions (class H): genus or splitting, experimentally.** The principal-cycle constraint
   ((−1)^iQ_i | p) = 1 (PROVED, classical) makes every quadratic symbol (ℓ|p) of small ℓ leak at n = 32–40. The
   pre-registered decisive test E2 found **no leak of non-quadratic characters of p** (cubic mod 7, quartic mod 5:
   max |z| 3.2, MDE |ρ| ≈ 0.03–0.04), while all four quadratic targets leaked (|z| up to 16). An archimedean target
   (q/p < 33/32) also leaked strongly (|z| = 21–35), but post hoc the whole effect sits on instances where a square
   form at an even index *splits N* (Lehman's identity inside the CF); on the non-split remainder |z| ≤ 2.5.
   So the CF prefix leaks (Q) genus symbols without splitting and (A)/(Q) square-form information only by splitting.
   That every CF statistic obeys this is CONJECTURED (H-GEN); it is not proved.
4. **AC⁰ on the bits of N cannot predict ε (PROVED, codex's argument, checked here).** The proof is unconditional:
   a character-pair trick reduces the needed Walsh–type-II estimate to bilinear prime sums Σχ₄(p)e(αpq), which are
   small on minor arcs by the elementary bilinear bound and on major arcs because χ and χχ₄ cannot both be the
   Landau–Page exceptional character. Thue–Morse and Rudin–Shapiro automata follow from Mauduit–Rivat/Müllner type-II
   bounds; general finite automata stay OPEN.
5. **Jacobi symbols of pure values (Q).** A Jacobi test (c|N) = (c|p)(c|q) on N ≡ 1 (mod 4) carries ε-information
   iff β(c)β(−c) ≠ 0, β(c) = E_x(c|x) (Theorem PT-J). That is impossible (negligible) for log|c| ≤ c₀√log X
   unconditionally (same character-pair trick) and for log|c| ≤ X^{1/4−δ} under GRH, possible for CRT-built c
   (codex), and OPEN only for SLP-compressed doubly-exponential c. E3 searched 345 such values (Fermat numbers,
   p−1 constants, Lucas values, factorials, ∏(x²+1), random SLPs) over 5·10⁵ primes: 0 hits.

**Classification theorem (§8).** In each of the models (i) pure trees with side information N mod m, (ii) one-layer
N-exponent trees with affine/dyadic exponents, (iii) non-adaptive Jacobi layers on pure values with
log|c| ≤ X^{1/4−δ} (GRH), (iv) AC⁰ circuits on the bits of N, every non-split leak is charged to (M) group-order
events, (A) Diophantine pins, or is negligible; (Q) appears in (iii) only through compressed c. These are separate
theorems for separate models, not one theorem for their union.

**The strongest natural model where the classification is open** (codex reached the same model independently):
ring operations in fixed Z/N-algebras, **retained and nested** N-th powers (x ↦ x^N whose output is reused,
e.g. (a^N + 1)^N, or several gates added together), and gcd/zero tests — the home of every primality test — plus
Jacobi symbols of compressed pure values. A breakthrough of the requested kind would have to be a deviation event
there that is neither a subgroup coset nor a short-interval pin in the exponent variable r = q mod (p − 1): an
exponential-polynomial equation mod p with ≫ p^{1−o(1)}/poly solutions r that are not explained by order alignment.
Nothing we know points to one; bounds of that type are the open "exponential congruence" problems (§3.6).

Subjective: probability that any route in this map yields sub-L[1/3] factoring < 1%.

---

## 1. The three families, operationally

In a tree with zero tests, everything is a function of N; "leak" is computational. Following factorbits §3, the
diagonal transcript τ_p of a prime p | N is the run answered by [p | c_v], where c_v is computed from N as the tree
computes it. For pure trees τ_p depends on p only (the class). Outside pure trees τ_p may depend on q, and we call
the dependence a *deviation*. The families:
- **(M)** a deviation or split explained by the order of a pure element at p (smoothness, small order) or by a
  common divisor of a group order at p and a fixed polynomial in q (gcd(N − 1, p − 1)).
- **(A)** a deviation explained by an exact small-height relation between p and q (q = p + 2, q = (j(p−1) + t − β)/α,
  |aq − bp| small).
- **(Q)** a quadratic character of p (genus symbol) or an ambiguous/square form.

## 2. Theorem PT-X (decoupling). PROVED

**Setting.** T is any deterministic tree whose node v computes an element c_v(N) ∈ A/NA (A a fixed free Z-algebra)
by *any* procedure from N, tests c_v = 0 and reports a split when the gcd of its coordinates with N is proper. For a
pair (p, q) let τ_p(p, q) be p's diagonal transcript. Let κ be any map from primes to transcripts and
D = {τ_p ≠ κ(p)} ∪ {τ_q ≠ κ(q)}.

**Lemma (diagonal lemma).** For p ≠ q the run on N = pq splits iff τ_p ≠ τ_q, and otherwise ends with transcript
τ_p. *Proof:* the gcd of the coordinates with pq is divisible by p iff c_v ≡ 0 (mod p); it is proper iff exactly one
of p, q divides it. Induct along the path. ∎ (This is factorbits' partition lemma with the class replaced by τ_p.)

**Theorem PT-X.** With p, q iid,

  Pr[T outputs F(p, q) or splits] ≤ B + 2s + 2·Pr[D ∧ no split] ≤ B + 2s + 2Pr[D].

*Proof.* Let π be the law of κ(p), ℓ* its heaviest class, η = 1 − π_{ℓ*}, s₀ = Pr[κ(p) ≠ κ(q)] = 1 − Σπ_ℓ² ≥ η.
On ¬D the lemma applies to κ: the run splits iff κ(p) ≠ κ(q), and otherwise outputs o(κ(p)). Hence
Pr[correct ∧ no split] ≤ Pr[D ∧ no split] + Pr[κ(p) = κ(q) = ℓ*, F = o(ℓ*)] + Σ_{ℓ≠ℓ*}π_ℓ²
≤ Pr[D ∧ no split] + B + η². Further s₀ ≤ s + Pr[κ(p) ≠ κ(q) ∧ no split] ≤ s + Pr[D ∧ no split], because
no split forces τ_p = τ_q, so κ(p) ≠ κ(q) forces a deviation. With η² ≤ η ≤ s₀ the claim follows. ∎

For uniform *distinct* ordered pairs from a finite prime set the same proof holds, since
c(c−1)/(n(n−1)) ≤ (c/n)² (used in E1). PT is the case D = ∅ (κ = class). **Meaning:** every bit of non-split
information is paid for by a deviation of one diagonal transcript from a one-prime reference. All results below
bound Pr[D] in specific models.

## 3. Class G: N as an exponent

### 3.1 Model NE1 (one-layer N-exponent trees)
- A fixed free Z-algebra A of rank r with A ⊗ F_x étale for x ∈ I (e.g. Z; Z[x]/(x² − bx − c) with bc-discriminant
  prime to I).
- *Pure values*: N-independent straight-line programs over A (any depth, coins fixed).
- *Gate node*: element P_v(Z) = Σ_{k≤D} c_k Z^k with pure coefficients c_k ∈ A, evaluated at the gate output
  Z = y^{E(N)} of a pure base y; element test P_v = 0 in A/NA. Pure nodes (no gate) are PT nodes.
- *Exponent menu*: E(N) = αN + β (α ≠ 0) or E(N) = (N − c)/2^j on N ≡ c (mod 2^j), or ⌊(N − 1)/2^j⌋ (fixed
  α, β, c, j). Gate outputs are not reused (one layer); several gate nodes may share a base.
- Covers: Fermat, Miller–Rabin ladders, the power in Solovay–Strassen, Lucas/strong Lucas, quadratic Frobenius,
  Baillie–PSW's two parts, and the AKS identity tested as an element (not coordinate-wise).

**The menu restriction is necessary (codex_1 §1, PROVED).** With an arbitrary computable exponent E(N) = 1 + h(N)
and the pure base −1, the test (−1)^{E(N)} + 1 = 0 reads any polynomial-time predicate h without splitting. In NE1
this is excluded: for affine/dyadic E and a base of public order o, E(N) mod o is a function of N mod o·2^j (PT-m
side information).

### 3.2 Counting lemma and Theorem NE1. PROVED
Fix a prime x ∈ I and the coins. Reference class κ(x): the formal transcript, in which a gate node answers "zero"
iff P_v ≡ 0 as a polynomial in Z in every component of A ⊗ F_x (i.e. all c_k ≡ 0 at x; a component where the base is
0 is substituted Z_i = 0). κ(x) depends on x only.

**Lemma NE-C.** Let v be a gate node whose formal test at x is nonzero, and let i be a component of A ⊗ F_x
(a finite field F_{x^d}) in which P_{v,i} ≢ 0 and y_i ≠ 0; let o_i = ord(y_i). Then the actual answer at x is
"zero" only if E(N) ≡ t (mod o_i) for one of at most D residues t. Hence

  Pr_q[node v deviates at x | x] ≤ D · μ_x(o_i),   μ_x(o) := max_t Pr_q[E(xq) ≡ t (mod o) | x].

*Proof.* "Zero" in A/xA needs P_{v,i}(y_i^{E(N)}) = 0 in F_{x^d}; P_{v,i} has ≤ D roots; each root in ⟨y_i⟩ fixes
E(N) mod o_i. ∎

**Theorem NE1.** For an NE1 tree of height h and any x, Pr_q[τ_x ≠ κ(x) | x] ≤ Σ_{v ∈ path κ(x)} D_v μ_x(o*_v(x)),
where o*_v(x) is the largest component order available in the lemma. With p, q iid, Theorem PT-X then gives

  Pr[correct or split] ≤ B + 2s + 4·E_x[ Σ_{v∈path κ(x)} D_v μ_x(o*_v(x)) ].

*Proof.* The first deviation along τ_x happens at a node of κ(x)'s path; union bound with Lemma NE-C; then
Pr[D] ≤ 2 E_x[...]. ∎

**Corollary NE1-bal (iid balanced primes).** For E(N) = αN + β, μ_x(o) ≤ |α|(|I|/o + 1)/π(I) (the congruence
αxq ≡ t − β mod o has ≤ gcd(α, o) solution classes, and gcd(x, o) = 1). For (N − c)/2^j and ⌊(N−1)/2^j⌋,
μ_x(o) ≤ (|I|/o + 2^j)/π(I). Splitting at a threshold L,

  Pr[D] ≤ 2 Pr_x[SO_L(x)] + 2h·D·|α|·(|I|/L + 2^j)/π(I),

where SO_L(x) = "some gate base on κ(x)'s path has o*_v(x) ≤ L" is a pure, one-prime event. With L = √X and bases
whose small-order mass is negligible (coins; fixed small integers, for which #{x : ord_x(y) ≤ L} ≤ L²), the
non-split leak is ≤ h·D·N^{−1/4+o(1)} + Pr[SO]: a polynomial-size NE1 tree learns nothing non-negligible without
splitting unless its pure part has already put a base into a small subgroup at p.

### 3.3 What the deviations are (classification). PROVED (definitions + Lemma NE-C)
A deviation at x of gate node v (base y, root w = y^t when the target is a known power of the base) is
- **SMALL** if o*_v(x) ≤ L. Pure (M): the pure part has annihilated all but L of the group order at x; such x are
  detected by the pure value ∏_{k≤L}(y^k − 1), a PT class event. (Example: y = 3^{lcm(1..64)}, the p−1 method.)
- **SUBGROUP** if o > L and the hit is the cross coincidence o | E(N) − t: gcd(E(N) − t, |G_x|) ≥ o. This is (M):
  a common factor of a group order at x and a fixed polynomial in the partner. For Fermat, o | gcd(N − 1, x − 1) =
  gcd(q − 1, p − 1).
- **UNIVERSAL** if exp(G_x) | E(N) − t, so every base is fooled. For |G_x| > |I| this pins the partner to O(1)
  integers: an exact small-height relation (A). Examples (all observed in E1, §9): ARCH node 2^N = 2^k at the smaller
  prime ⇔ q = p + k − 1; Lucas α^{N+1} = 1 in the torus at both primes ⇔ |T_p| = p+1 and |T_q| = q−1 divide
  (p+1)² ⇔ q = p + 2 with (5|p) = −1, (5|q) = +1 (a two-sided, non-split, archimedean leak; it is of course useless:
  q − p = 2 factors N by Fermat).

So **in NE1 the conjecture holds in the form: every non-split leak is (M), and its degenerate end is (A)**; (Q) does
not occur. The counting is one-sided (it charges a deviation at either prime). The two-sided product form:

**Theorem NE1-R (codex_1 §1, "terminal probes"; PROVED there, re-checked here).** If every gate node is a probe
y^{E} = 1 with a *fresh uniform* base y in the local unit groups (independent CRT components), the gate value is
discarded after the test, and the probe is gcd-checked, then

  Pr[correct or split] ≤ B + 2s + 2R + 2/M,   R = E Σ_{probes reached} κ_p(E)κ_q(E),
  κ_r(E) = ∏_i gcd(E, d_{r,i})/d_{r,i} for G_r ≅ ∏ C_{d_{r,i}}.

For scalar Fermat probes κ_pκ_q = gcd(N−1, p−1)gcd(N−1, q−1)/φ(N), the rate seen by the eps battery (PC4/MR)
and derived exactly in primality §1 (order-pattern lemma). *Proof sketch:* couple with the pure shadow tree in which
every probe answers "false"; the first disagreement is a split or a both-pass event, whose probability given the
history is κ_pκ_q by freshness. ∎ Our Theorem NE1 is the deterministic-coin, arbitrary-target version; NE1-R is
the sharp random-base version.

### 3.4 Several gates in one node (F_p only). PROVED from Theorem SP / Kelley 2016
Let a node be P(y_1^{E_1(N)}, …, y_k^{E_k(N)}) with affine E_j and all y_j ∈ F_x^*. Since x ≡ 1 (mod x−1),
E_j(xq) ≡ α_j q + β_j (mod x − 1), so as a function of r = q mod (x−1) the node is a Laurent polynomial
f(g^r) = Σ_m γ_m g^{e_m r} with T ≤ C(D + k, k) terms (g a generator). By Theorem SP (slpcirc §4.2; Kelley,
LMS J. Comput. Math. 19A (2016), Thms 2.2–2.3), zero density ≥ λ on F_x^* forces every prime ℓ > (2/λ)^{T−1}
dividing x − 1 to divide some exponent difference e_m − e_{m'}, i.e. some monomial ratio ∏ y_j^{α_j(m_j−m'_j)} has
order prime to ℓ: an (M) alignment. Otherwise the zero set has < λ(x−1) residues and
Pr_q[deviation] ≤ λ(x − 1)(|I|/(x−1) + 1)/π(I) ≤ 3λ ln X. With λ = (ln X)^{−2} the non-aligned contribution is
O(1/ln X) per node: weaker than NE1 but the same classification.

### 3.5 Exponents larger than the partner's range
In F_{x^d}, d ≥ 2, the exponent residue E(N) mod o can range over o > |I|; Lemma NE-C still applies to element
tests (μ ≤ 2|α|/π(I)). For **coordinate** tests (one coefficient of an element of Z/N[X]/(X^r − 1), as in AKS) the
node is a sparse exponential sum over Galois conjugates evaluated along the short progression {xq}; incomplete
sums of this type have no general bound. primality §3 measured exactly such a coefficient leak, confined to the
Kummer regime q − p ≲ r: an (A) leak, consistent with the classification but not covered by a theorem.

### 3.6 Nested and retained gates: OPEN
For (a^N + 1)^N, or a node adding two gate outputs of different bases with retained values, the p-side transcript
is a function of r = q mod (x−1) through exponential-polynomial towers such as r ↦ (g^r + 1)^r. A deviation set of
size ≫ |I|^{1−o(1)}·poly among the r ∈ I mod (x−1), not explained by order alignment or a short-interval pin, would
be a non-(M), non-(A) leak. Bounding such sets is the problem of counting solutions of exponential congruences
(of the type g^r ≡ r, r^r ≡ c mod p; known results give p^{1−δ}-type bounds for some single equations; literature
status UNVERIFIED here, see §11). **This is the frontier for class G.**

## 4. Integer representatives

### 4.1 Prop REP-C (completeness). PROVED (codex_1 §2)
- A tree that can compare rep(x) ∈ [0, N) with constants reads rep(−1) = N − 1 in n comparisons and can then label
  its leaves by any function of N. PT is size-free, so **no size-free PT-type statement exists with representatives.**
- Even the least significant bits b_j = LSB(2^j mod N) = ⌊2^j/N⌋ mod 2, j ≤ 2n + 2, of pure units determine N
  (they are the binary digits of 1/N).
- Poly-size uniform trees with representatives and Boolean gates realised as x ∧ y = xy on {0,1}-residues simulate
  every polynomial-time algorithm. So a classification of leaks in this model is a classification of all factoring
  algorithms (cf. selfred SR3: unconditional statements would separate P from NP-type classes).
- All three families live here: QS/NFS (smoothness of the integer representative), Lehman/Harvey (integer
  arithmetic on kN), SQUFOF (forms), Coppersmith (lattices). The (M) family *without* representatives is only the
  group-order part (ECM is pure); relation (M) needs representatives.

### 4.2 Lemma GAR (one hidden phase). PROVED
For x ∈ Z/N with components x_p ∈ [0, p), x_q ∈ [0, q), and θ := (q^{−1} mod p)/p,

  rep(x)/N ≡ (x_p − x_q)·θ + x_q/N  (mod 1),   and θ = j/q + 1/N with qθp = 1 + jp, 0 < j < q.

*Proof.* With q̄ = q^{−1} mod p, p̄ = p^{−1} mod q, the integer qq̄ + pp̄ is ≡ 1 mod p and mod q and lies in (1, 2N),
so it equals N + 1; thus q̄/p + p̄/q = 1 + 1/N. CRT gives rep(x) ≡ x_p qq̄ + x_q pp̄ (mod N); divide by N and
substitute p̄/q = 1 + 1/N − θ. The second claim is qq̄ = 1 + jp divided by N. ∎

**Consequences.**
- One-shot readouts of representatives of pure values see the one-prime classes (x_p, x_q) and *one* real number θ,
  a rational of denominator q up to 1/N. All cross-prime information is Diophantine information about θ: family (A).
- θ exactly gives gcd(θN − 1, N) = p. θ to k top bits gives θN = 1 + jp up to N·2^{−k}: an approximate multiple of
  p, so k ≥ 3n/4 + O(1) bits factor N by Howgrave-Graham's partial approximate common divisor algorithm
  (Coppersmith's method; Howgrave-Graham, CaLC 2001) — **KNOWN method**, the same N^{1/4}-error barrier as n/4
  known bits of p (topblock §2.3).
- Readouts that feed representatives back into ring arithmetic are complete (4.1); Lemma GAR says nothing about them.

### 4.3 The archimedean algorithms as covering events (map, not a theorem)
- Lehman/Harvey/our DCL use the identity (aq + bp)² − 4abN = (aq − bp)² (codex): a small |aq − bp| puts a
  factor-dependent integer in a short searchable interval. Success = "q/p is covered by a rational b/a of small
  height", a Diophantine event. coverage38 Theorem U bounds the cost of any covering design in that family by
  N^{1/5}(λ⁴/J)^{1/5}; dchworst Theorem R is the recovery step.
- Coppersmith: success = "p lies within N^{1/4} of an enumerated approximation"; topblock Prop T shows window
  lattices add nothing to the univariate lattice.
- In E2 the archimedean target q/p < 33/32 is predicted by CF square forms only on instances that those square forms
  split (§5.2): (A) inside the CF enters through Lehman's identity, i.e. by splitting.
- **OPEN (and unformalisable as stated):** "every leak in the representative model is Diophantine".

## 5. Class H: continued fraction of √N and infrastructure

### 5.1 What is proved (classical; codex_1 §3 gives the proofs)
- A_{i−1}² − NB_{i−1}² = (−1)^iQ_i, hence ((−1)^iQ_i | p) = ((−1)^iQ_i | q) = 1 when gcd(Q_i, N) = 1 (genus).
  Odd i with Q_i = k² certifies ε = +1; even i with Q_i = ℓk² certifies (ℓ|p) = +1; odd i with Q_i = ℓk² certifies
  (−ℓ|p) = +1 (E2 controls, 0 violations in 8·10⁵ instances).
- Products of Q_i that are squares give z² = ±1 mod N: a √−1 (ε certificate) or a congruence of squares (split,
  CFRAC). Square forms give SQUFOF splits.
- Random-walk models for the principal cycle (Gower–Wagstaff, Math. Comp. 77 (2008), §4) are heuristic; a model
  "subject only to genus constraints" is literally false (period, reversal symmetry; codex).

### 5.2 E2: the CF prefix leaks genus symbols, not other characters (experimentally supported)
Pre-registered (§9). 16 CF statistics × 7 targets, 2·10⁵ instances per size, fresh-seed replication.

| target | type | n = 32 hits (confirmed) | max |z| n=32 | n = 40 hits | max |z| n=40 | MDE |ρ| |
|---|---|---|---|---|---|---|
| G1 ε = (−1|p) | genus | 9 (9) | 15.9 | 4 (4) | 8.7 | 0.019 |
| G2 (−3|p) | genus | 3 (3) | 13.8 | 1 (1) | 5.8 | 0.019 |
| G3 (5|p) | genus | 1 (1) | 9.0 | 0 | 2.7 | 0.026 |
| G4 (−7|p) | genus | 1 (1) | 6.0 | 0 | 3.0 | 0.032 |
| NG1 cubic χ mod 7 | non-genus | 0 | 3.2 | 0 | 2.6 | 0.032 |
| NG2 quartic χ mod 5 | non-genus | 0 | 1.7 | 0 | 2.3 | 0.037 |
| A1 q/p < 33/32 | archimedean | 5 (5) | 21.5 | 4 (4) | 34.3 | 0.013 |

(n = 32: two-way cluster-robust z; n = 40: permutation z. Strata N mod 2^8·m; for A1 also the top 6 bits of N.)

**Post hoc (not pre-registered; `e2_posthoc.py`):** square forms Q_i = r² at even i ≤ 64 occur in 28.7% (n=32) /
8.1% (n=40) of instances; 56–58% of them split N via gcd(A_{i−1} − r, N). P(A1) is 0.42 on the split instances,
0.15–0.19 otherwise. On the instances with no splitting square form the A1 correlation of cf64_n_sq_even is
z = −2.4 (n=32) and −1.1 (n=40). The archimedean leak is entirely a splitting leak (Lehman via the CF).

**Conjecture H-GEN (refined).** For poly(n)-length CF/infrastructure prefixes of √(kN), k fixed, every
polynomial-time statistic predicts a function of the factors beyond the fiber floor and beyond the events where the
prefix's square relations split N only through the quadratic symbols (ℓ|p), ℓ ∈ {−1} ∪ {primes ≤ poly}, constrained
by ((−1)^iQ_i | p) = 1. Status: CONJECTURED; supported by E2 (non-quadratic characters null at MDE 0.03–0.04;
genus targets hit). Codex notes that without a postprocessing restriction the statement is vacuous (Q_1 = N − P_1²
already encodes N); "statistic" must mean a symmetric low-complexity functional of the (a_i, Q_i) multiset, which we
have not formalised.

### 5.3 Why this is (Q) and not new
The leaked symbols (ℓ|p) for small ℓ are fixed-base quadratic residuosity bits (nonsmooth Prop RANK; Seurin's
2-Φ/4-hiding for ℓ = −1). Their strength decays with N because small Q_i become rare (eps §6). Getting one symbol
with constant advantage costs SQUFOF-type N^{1/4} or CFRAC-type L[1/2] work (eps §6, PC2/PC3).

## 6. Class F: AC⁰ and automata

### 6.1 Theorem AC0-ε (codex_1 §4; PROVED there; proof re-derived and checked here, see §11 for critic status)
For iid uniform primes p, q ∈ (X, 2X] conditioned on N ≡ 1 (mod 4), every family of depth-d, size-s Boolean
circuits on the bits of N with s = poly(n) satisfies Pr[A(N) = ε(p)] = ½ + o(1) (indeed + O_C(n^{−C}) for every C).

*Proof outline.*
1. Character pair: χ₄(p)·[pq ≡ 1 (4)] = (χ₄(p) + χ₄(q))/2, so Σ_{pq≡1(4)} χ₄(p)K(pq) = Σ_{p,q} χ₄(p)K(pq).
2. Walsh functions W_S of b-bit integers have additive Fourier ℓ¹-norm ≤ (Cb)^{|S|} (Bourgain's Lemma 1).
3. Σ_{p,q} χ₄(p)e(αpq) ≪ X²e^{−c(log X)^{1/4}} uniformly in α: minor arcs (r ≥ Q = e^{L^{1/4}}) by the bilinear
   bound X⁴/r + X³log r + X²r log r with arbitrary bounded coefficients; major arcs by expanding e(at/r) in
   characters χ mod lcm(4, r): each term factors as (Σ_p χχ₄(p))(Σ_q χ(q)), and χ, χχ₄ cannot both be principal or
   the Landau–Page exceptional character ξ (ξ ≠ χ₄ since L(s, χ₄) > 0 on (0, 1)).
4. Hence β_t = max_{|S|≤t}|E[ε W_S(N) | N ≡ 1 (4)]| ≤ exp(−c b^{1/4} + O(t log b)).
5. LMN: the degree > t Fourier tail of A is ≤ 2s·2^{−t^{1/d}/20}; the semiprime law has density O(b²) w.r.t. the
   uniform cube, so |E[εA]| ≤ √(Σ_{j≤t}C(b,j))·β_t + O(b·tail^{1/2}). Take t = polylog. ∎

**Scope.** This is the target ε on this iid law. It extends verbatim to any target χ(p) with χ a fixed real
character such that the pair trick applies; other targets (bits of p) need other arguments. Finite automata:
Thue–Morse and Rudin–Shapiro follow from Mauduit–Rivat (JEMS 2015, Prop. 2) and Müllner (Duke 2017, Prop. 5.5),
which allow arbitrary bounded coefficients on balanced ranges (codex; hypotheses re-checked only at the level of
statements). General automata and AC⁰[⊕]: OPEN. The selfred §3.2 "OPEN" entry for AC⁰ is thereby closed (pending
the critic, §11).

## 7. Jacobi symbols of pure values (Q)

**Theorem PT-J. PROVED.** Population: e ∈ {±1} uniform, p, q iid from μ_e (the primes of I with ε = e), so
N ≡ 1 (mod 4). Observations: J_i = (c_i | N) = λ_i(p)λ_i(q), i ≤ k, for N-independent integers c_i that are units
at every prime of I (otherwise a gcd splits). For S ⊆ [k] put c_S = ∏_{i∈S} c_i and m_e(S) = E_{μ_e}[(c_S | x)]. Then

  Pr[correct] ≤ ½ + ¼ Σ_{∅≠S⊆[k]} |m_+(S)² − m_−(S)²|,

and when μ gives ε = ±1 equal mass, m_±(S) = β(c_S) ± β(−c_S) with β(c) := E_μ[(c|x)], so the bound is
½ + Σ_{S≠∅} |β(c_S)β(−c_S)|. A single test leaks iff β(c)β(−c) ≠ 0.

*Proof.* The law of u = λ(p)⊙λ(q) ∈ {±1}^k under e has Fourier coefficients E_e[χ_S(λ(p))]² = m_e(S)².
Bayes success is ½ + ¼Σ_u|P_+(u) − P_−(u)| ≤ ½ + ¼Σ_S|P̂_+(S) − P̂_−(S)|. With equal masses,
m_e(S) = β(c_S) + eβ(−c_S) since (−c|x) = ε(x)(c|x), and m_+² − m_−² = 4β(c_S)β(−c_S). ∎ (Unit-tested by brute force.)

**Corollary PT-J1.** (a) For log|c_S| ≤ c₀√(log X), unconditionally |β(c_S)β(−c_S)| ≤ e^{−c√log X}: (c_S|·) and
(−c_S|·) = χ₄(c_S|·) are characters of conductor ≤ 4|c_S| whose ratio is χ₄, so at most one is principal or
Landau–Page exceptional (the §6 pair trick). (b) Under GRH, for every non-square c with log|c| ≤ X^{1/4−δ},
|β(c)| ≪ X^{−1/2}log³(|c|X), so k = O(log n) tests leak negligibly. (c) Without restriction PT-J's criterion is met:
codex's CRT-constructed c (QNR exactly on half of the ε = +1 primes) gives success → ¾ — but c has no small SLP.

**OPEN:** SLP-compressed c with log|c| ≥ X^{1/4} (e.g. 2^{2^k} + 1 with 2^k > X^{1/4}). Whether a poly-size SLP can
make both (c|x) and (−c|x) biased over primes is a clean question about Legendre symbols of compressed integers.
**E3 (pre-registered):** 459 menu values, 114 excluded as structural squares (V_{2m}(3) + 2 = V_m², and 88 random SLPs
ending in a squaring), 345 tested over all 513 708 primes in [2^23, 2^24): **0 hits** (max min(|z₁|,|z₂|) = 2.27,
max |z₁z₂| = 5.3). No compressed value in the menu has detectable Legendre bias on both c and −c.

## 8. Synthesis: classification and the open frontier

**Theorem CLASS (union of §§2–7; PROVED in each model separately).** Let X be one of:
1. pure trees with side information N mod m (factorbits PT, PT-m, PT2);
2. NE1 trees (affine/dyadic exponents, one gate per node, element tests), any pure depth, adaptive (Theorem NE1;
   NE1-R for fresh bases; §3.4 for several F_p gates);
3. a non-adaptive layer of k = O(log n) Jacobi tests of pure integers with log|c| ≤ c₀√log X (unconditional) or
   ≤ X^{1/4−δ} (GRH) (Theorem PT-J);
4. AC⁰ circuits on the bits of N, target ε (Theorem AC0-ε);
5. functions of N mod m, m ≤ X^{1−η} (GRH) or low bits 2^j ≤ N^c (selfred Theorem L).
Then every non-split prediction advantage beyond the fiber floor is bounded by (M) terms (small-order pure events;
cross gcd coincidences gcd(E(N) − t, |G_p|) ≥ ord), (A) terms (Diophantine pins exp(G_p) | E(N) − t), or is
negligible (o(1), mostly N^{−Ω(1)} or e^{−c(log N)^{1/4}}). (Q) is unreachable in 1–5.

These are separate models; the union (e.g. Jacobi tests after an NE1 tree, or AC⁰ on the bits of pure residues) is
not covered, and representatives (§4) make any union containing them complete.

**Strongest natural open model.** Polynomial-size uniform programs with ring operations in fixed Z/N-algebras,
retained and nested gates x ↦ x^N (and x^{N±1}, x^{(N−1)/2}), gcd/zero tests, and Jacobi symbols of computed values
whose integer representatives are compressed (SLPs). Call it PRIM: it contains every primality test, Lucas
sequences, AKS, Frobenius tests and Solovay–Strassen with its Jacobi symbol. Our theorems cover its one-layer part
(NE1), its terminal Jacobi part with moderate conductors (PT-J), and nothing else. Codex_1 independently names the
same model ("+, −, ×, x ↦ x^N, gcd/zero tests, no representatives").

**Where a breakthrough would have to live.** By PT-X, a non-splitting leak in PRIM is a deviation of one prime's
transcript from its one-prime reference. In NE1 those are cosets of subgroups (M) or pins (A). A new family would
be a retained-gate computation whose zero set in r = q mod (p − 1) is large (≥ |I|^{1−o(1)}/poly solutions in the
relevant progression) and unexplained by order alignment or short intervals — or a compressed c with
β(c)β(−c) ≠ 0. Neither is excluded by anything we can prove; E1 (retained gates absent) and E3 (compressed c) give
no hint of one.

## 9. Experiments (pre-registered in `notes/leakmap_prereg.md`, commit d050db1)

### E1: instrumented N-exponent trees
Populations: all pairs of b-bit primes with top two bits set, b = 12, 14 (iid), 2·10⁵ sampled pairs at b = 16, and
the Fermat-regime pairs q − p ≤ 2^{b/2}. 7 trees (FERMAT, MR2 with ⌊(N−1)/2^j⌋ exponents, LUCAS, QFROB, ARCH,
REPAIR, PM1). 42 cells, 1.9·10⁶ pair-runs (×3 transcripts each).
- **P1 (diagonal lemma): 0 mismatches** in every cell. Every split is STRUCTURAL (τ_p ≠ τ_q); ACCIDENTAL: 0.
- **P2 (PT-X): holds in all 126 (cell × target) evaluations** (Bayes-optimal leaf labels). Largest non-split gain:
  pure p−1 trees (≈ 0.02–0.03, PT's known smoothness/residue-class correlation), NE trees ≤ 0.004 except MR2 on the
  b = 12 Fermat population (0.0077).
- **P3 (NE1 counting bound): 0 violations** over all primes in all cells; worst ratio freq/bound 0.74.
- **P4 (random bases):** b = 14, 2·10⁴ pairs × 20 bases: Fermat liars 7 vs 9.7 expected, MR strong liars 4 vs 3.9
  (both z < 1; underpowered as registered). Disclosed added cell b = 10, all 666 pairs × 2000 bases: 359 vs 344.4
  (z = 0.79) and 171 vs 159.0 (z = 0.95). Passed.
- **P5:** the registered claim "UNIVERSAL deviations occur only for ARCH in iid populations" is **refuted as
  registered** (a prediction error, not a theorem failure): LUCAS has UNIVERSAL deviations at twin primes (q = p + 2,
  (5|p) = −1, (5|q) = +1; two-sided, non-split) and MR2's floor exponents ⌊(N−1)/2^j⌋ give pins as well. The exact
  ARCH pin count holds in all 6 cells (pins = UNIVERSAL pairs: 432, 1307, 682, 432, 1307, 3870). The directional
  prediction holds: the UNIVERSAL share is larger in the Fermat population for every tree and size where UNIVERSAL
  occurs (e.g. LUCAS 0.21–0.50 → 0.95–0.96; ARCH 0.33–0.50 → 0.96–0.98).

### E2: CF genus exclusivity — table in §5.2
- Q1 (G1 replicates at n = 32): **yes**, 9/9 confirmed. Q2 (G2–G4 reach): **yes**, 5 confirmed hits.
- **Q3 (decisive): NG1, NG2 0 hits at both sizes, max |z| 3.25 < 4 (fresh 3.20). H-GEN not refuted.**
- Q4 (decay): genus hits 14 at n = 32, 5 at n = 40. Holds.
- Controls PC{3,5,7}_{even,odd}: 0 violations in 8·10⁵ instances (fires 13 604–48 868 per size and control).
- A1 (exploratory): 5 and 4 confirmed hits; explained post hoc by splitting square forms (§5.2).

### E3: Legendre bias of compressed pure values — §7. R1 (0 hits): **held.**

### Deviations
- D1: P4 b = 14 cell was underpowered (expected counts < 10); an extra b = 10 cell was added after the run and is
  labelled as such.
- D2: E2 post-hoc split analysis of A1 (not registered).
- D3: P5's first sentence was a mis-specified prediction; reported as refuted.
- No registered statistic or cell was dropped.

## 10. Candidate cards (rule 9)

| Mechanism | Model | Outcome |
|---|---|---|
| Fermat/MR/Lucas/Frobenius exponent N ± 1 as a factor leak | NE1 | **refuted for this model** beyond (M)/(A): Theorem NE1, NE1-R; E1 |
| Archimedean exponent probes 2^N = 2^k | NE1 | (A) pins q − p = k − 1, exactly as predicted; cost linear in the window |
| Arbitrary exponent selectors E(N) | NE with free E | predicate-complete (codex); not a mechanism |
| Retained/nested N-th powers | PRIM | **OPEN** (frontier) |
| Integer-representative phase θ | one-shot readouts | (A); needs 3n/4 bits of θ (HG ACD) |
| CF prefix statistics | H | genus symbols leak, non-quadratic characters do not (E2); (A) only via splits |
| AC⁰ on bits of N | F | **refuted** (Theorem AC0-ε, unconditional) |
| Jacobi of pure values | PT-J | refuted for log|c| ≤ X^{1/4−δ} (GRH); open for compressed c; E3 null |

## 11. Reviews
(To be filled after the codex adversarial round `codex/leakmap_2` and the two Opus critics.)

## 12. Reproduction (project root)
```
.venv/bin/python -m pytest -q tests/test_leakmap.py                     # 20 tests, ~1 s
.venv/bin/python -m src.leakmap.e1                                      # E1 -> experiments/raw/leakmap_e1*.jsonl (~4 min)
.venv/bin/python -m src.leakmap.e1_p5                                   # P5 -> experiments/analysis/leakmap_e1_p5.json
.venv/bin/python -m src.leakmap.cfgenus run --n 32 --range main         # E2 (~1 min on 6 workers), same for n 40, --range fresh
.venv/bin/python -m src.leakmap.cfgenus analyse experiments/raw/leakmap_e2_main_n32.jsonl.gz --cluster --out experiments/analysis/leakmap_e2_main_n32.json
.venv/bin/python -m src.leakmap.e2_posthoc experiments/raw/leakmap_e2_main_n32.jsonl.gz experiments/raw/leakmap_e2_main_n40.jsonl.gz
.venv/bin/python -m src.leakmap.jacobi --bits 24                        # E3 (~1 min)
```
Peak RSS < 400 MB for every job.
