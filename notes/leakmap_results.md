# leakmap: where "information costs splitting" stops (phase 24)

Agent `leakmap`, 2026-10-23. **Version v3**: after the codex adversarial round (`codex/leakmap_2.out`, MAJOR REVISION),
the Opus proof critic (ACCEPT WITH FIXES) and the Opus literature critic (OK WITH FIXES); all findings applied (§11).
- Pre-registration: `notes/leakmap_prereg.md` (commit `d050db1`, before any code or experiment).
- Code: `src/leakmap/{netrees,e1,e1_p5,e1_mrl,cfgenus,e2_summary,e2_posthoc,jacobi}.py`.
  Tests: `tests/test_leakmap.py` (20 pass).
- Raw data: `experiments/raw/leakmap_e1*.jsonl`, `leakmap_e2_{main,fresh}_n{32,40}.jsonl.gz`, `leakmap_e3.jsonl`.
  Summaries: `experiments/analysis/leakmap_*.json`.
- GPT-6.1-sol: blind round `codex/leakmap_1.{md,out}` (several results below are codex's and are credited as such);
  adversarial round `codex/leakmap_2.{md,out}`.
- Labels: **PROVED** (proof here or in the cited codex file, re-checked by the proof critic; "no exact prior
  formulation located" is the strongest novelty claim), **KNOWN** (reference), **HEURISTIC**, **CONJECTURED**,
  **OPEN**, **experimentally supported**, **refuted for this model**, **UNVERIFIED**.

Notation. N = pq, p < q; for theorems p, q are independent primes from laws μ₁, μ₂ on the primes of an interval
I = (X, 2X] (iid when μ₁ = μ₂). π(I) = number of primes in I. A *tree* is a finite adaptive decision tree whose
nodes compute an element of a fixed (Z/N)-algebra, test it for zero, and report a split when the gcd of its
coordinates with N is proper. s = Pr[split]; B = max_b Pr[F = b] for the target F.

---

## 0. Verdict (plain language)

**No fourth leak family was found.** What we can now say precisely, model by model:

1. **One general principle (Theorem PT-X, PROVED).** In *any* tree with zero tests mod N, every bit learned about
   the factors without splitting is paid for by a *deviation*: one prime's diagonal transcript differing from a
   reference transcript that depends on that prime alone. PT (factorbits) is the case with no deviations.
2. **N as an exponent (class G), one layer (Theorem NE1, PROVED).** Trees that test pure values raised to fixed
   exponents E(N) = αN + β or ⌊(N − c)/2^j⌋ (Fermat, Miller–Rabin, the power in Solovay–Strassen, Lucas, quadratic
   Frobenius, AKS tested as an element) deviate only when the gate output hits one of ≤ D roots, i.e. when
   E(N) ≡ t mod ord_p(y). For iid balanced primes this costs ≤ h·D·g·(|I|/L + 1)/π(I) per tree beyond a
   *pure small-order event* (a gate base of order ≤ L at p): N^{−1/4+o(1)} for random bases, N^{−1/6+o(1)} for
   fixed integer bases. With fresh random bases the both-primes rate is exactly κ_pκ_q, e.g.
   gcd(N−1, p−1)gcd(N−1, q−1)/φ(N) for Fermat (Monier, TCS 12 (1980); Baillie–Wagstaff, Math. Comp. 35 (1980) for
   Lucas; the tree-level bound NE1-R is codex's). So the conjectured statement holds for
   this model in the form of a quantitative bound. Interpreting the events: small order and common divisors of a
   group order with E(N) − t are group-order (M) events, and their extreme case exp(G_p) | E(N) − t pins the partner
   exactly (q = p + 2 for Lucas with (5|p) = −1, (5|q) = +1; q = p + k − 1 for 2^N = 2^k), an archimedean (A) event.
   This labelling is an interpretation of a case split, not a theorem that M and A are disjoint mechanisms; quadratic
   characters do enter Euler-type gates (a^{(N−1)/2} = (a|p)a^{(q−1)/2} at p), they just need no separate term.
   The exponent menu must be fixed: an arbitrary exponent E(N) = 1 + h(N) on the base −1 reads any poly-time
   predicate h without splitting (codex). Side information N mod m and the target ε need an NE1 analogue of PT2,
   which is OPEN (§3.3).
3. **Integer representatives are computationally complete**, so no PT-type theorem exists there (codex: even the
   low bits of 2^j mod N give N). What survives is an identity: every representative satisfies
   rep(x)/N ≡ (x_p − x_q)θ + x_q/N (mod 1) with one hidden phase θ ≈ j/q, and (3/4 + ε)n bits of θ suffice to
   factor (Howgrave-Graham; May), the same N^{1/4} wall as Coppersmith. This is not a theorem that leaks must go
   through θ.
4. **Continued fractions (class H): genus or splitting, experimentally.** The principal-cycle constraint
   ((−1)^iQ_i | p) = 1 (PROVED, classical) makes quadratic symbols leak: all four registered quadratic targets
   hit at n = 32 (|z| up to 14 with a conservative prime-sharing variance) and two at n = 40. The two registered
   non-quadratic targets (a cubic character mod 7, a quartic one mod 5) were null at both sizes (max |z| 3.3; detection
   limit |ρ| ≈ 0.15 at n = 32 after the prime-sharing design effect, ≈ 0.05 at n = 40). An archimedean target
   (q/p < 33/32) leaked strongly (|z| up to 34), but post hoc every one of the 16 statistics is null (max |z| 2.9) on
   the instances where no even-index square form splits N: here (A) entered the CF only through Lehman-type splits.
   "Every CF statistic leaks only genus symbols without splitting" is CONJECTURED (H-GEN), not proved.
5. **AC⁰ on the bits of N cannot predict ε (PROVED; codex's argument, checked by codex's own adversarial round and
   by the Opus proof critic).** Unconditional, with advantage O_C(n^{−C}): a character-pair identity reduces it to
   bilinear prime sums Σχ₄(p)e(αpq), small on minor arcs by the elementary bilinear bound and on major arcs because
   χ and χχ₄ cannot both be the Landau–Page exceptional character. This is for the target ε; AC⁰ *does* leak
   factor sizes (codex: the top bits of N predict min(p, q) > 3X/2 with constant advantage), an elementary (A) leak.
6. **Jacobi symbols of pure values (Q) (Theorem PT-J, PROVED).** (c|N) = (c|p)(c|q) on N ≡ 1 (mod 4) carries
   ε-information iff β(c)β(−c) ≠ 0, β(c) = E_x(c|x). This is negligible for log|c| ≤ (log X)^{1−δ} unconditionally
   and for log|c| ≤ X^{1/4−δ} under GRH, fails for CRT-built c (codex), and is OPEN for SLP-compressed c of
   doubly exponential size. E3 tested 345 such values over 5·10⁵ primes: 0 hits.

**The strongest natural model where the classification is open** (codex reached it independently): uniform
programs over fixed Z/N-algebras with N-independent shape whose only N-dependence is exponent gates x ↦ x^{E(N)}
from a fixed menu, **with gate outputs retained and fed into further gates or sums** (e.g. (a^N + 1)^N), plus
gcd/zero tests and Jacobi symbols of computed values. Call it PRIM ("generic ring + N-th power map"). We prove the
one-layer part (NE1) and terminal Jacobi tests with moderate conductors (PT-J) and nothing more. A new family would
have to appear there as a deviation whose zero set in r = q mod (p − 1) is neither a coset of an order subgroup nor
a short-interval pin — a counting question for exponential congruences such as (g^r + 1)^r ≡ c (mod p) — or as an
SLP-compressed c with Legendre bias on both c and −c. Nothing we found points to either.

Subjective: probability that any route in this map yields sub-L[1/3] factoring < 1%.

---

## 1. The three families, operationally

Everything a tree computes is a function of N, so "leak" is computational. The *diagonal transcript* τ_p of a
prime p | N is the run answered by [p | c_v], where c_v is computed from N exactly as the tree computes it. For pure
trees τ_p depends on p alone (factorbits' class); otherwise it may depend on q, and a *deviation* is a difference
from a chosen one-prime reference. Families, as used here (descriptive labels, not formal definitions):
- **(M)** events explained by the order of an element at p: a pure base of small order (smoothness), or a large
  common divisor of a local group order and E(N) − t (gcd(N − 1, p − 1)).
- **(A)** events explained by an exact small-height relation between p and q (q = p + 2,
  (p − 1) | αq + β − t with small quotient, |aq − bp| small), or by factor sizes.
- **(Q)** quadratic characters of p (genus symbols), ambiguous or square forms.

## 2. Theorem PT-X (decoupling). PROVED

**Setting.** T is a deterministic tree with fixed leaf outputs; node v computes c_v(N) ∈ A/NA (A a fixed free
Z-algebra) by any procedure, the same value being used for the zero test and for the gcd. For p ≠ q let τ_p(p, q)
be p's diagonal transcript. Let κ map primes to transcripts and D = {τ_p ≠ κ(p)} ∪ {τ_q ≠ κ(q)}.

**Lemma (diagonal lemma).** For p ≠ q the run on N = pq splits iff τ_p ≠ τ_q, and otherwise ends with transcript
τ_p. For N = p² a non-split run follows τ_p. *Proof:* p divides the coordinate gcd iff c_v ≡ 0 (mod p), and the gcd
with pq is proper iff exactly one of p, q divides it; for p², gcd = p is a split and otherwise [p² | c] = [p | c].
Induct along the path. ∎ (factorbits' partition lemma with the class replaced by τ_p.)

**Theorem PT-X.** For independent p ~ μ₁, q ~ μ₂ (p = q allowed), with α, β the laws of κ(p), κ(q) and
ℓ* = argmax α:

  Pr[T outputs F(p, q) or splits] ≤ B + 2s + 2·Pr[D ∧ no split].

*Proof.* Put s₀ = Pr[κ(p) ≠ κ(q)] = 1 − Σ_ℓ α_ℓβ_ℓ. On ¬D the lemma applies to κ: no split forces κ(p) = κ(q),
and then the output is o(κ(p)). So Pr[correct ∧ no split] ≤ Pr[D ∧ no split] + Pr[κ(p) = κ(q) = ℓ*, F = o(ℓ*)]
+ Σ_{ℓ≠ℓ*}α_ℓβ_ℓ ≤ Pr[D ∧ no split] + B + (1 − α_{ℓ*}). Next, 1 − α_{ℓ*} ≤ 1 − Σα_ℓβ_ℓ = s₀, and
s₀ ≤ s + Pr[κ(p) ≠ κ(q) ∧ no split] ≤ s + Pr[D ∧ no split] (no split gives τ_p = τ_q, so differing references force
a deviation). Adding Pr[split] = s gives the claim. ∎

- For iid laws the middle term can be sharpened to η² as in PT; the stated form suffices here.
- Uniform *distinct* ordered pairs from a finite prime set (used in E1) also satisfy it: the same-class
  probabilities c(c−1)/(n(n−1)) are ≤ (c/n)².
- **Meaning.** Every bit of non-split information is charged to deviations of one diagonal transcript from a
  one-prime reference. §3 bounds Pr[D] for N-exponent trees. With independent laws, PT-X also gives the fiber version
  (condition on (p mod m, q mod m) = (u, c/u)), as PT-m does.

## 3. Class G: N as an exponent

### 3.1 Model NE1 (one-layer N-exponent trees)
- A fixed free Z-algebra A of rank r with A ⊗ F_x étale for x ∈ I (e.g. Z; Z[x]/(x² − bx − c) with discriminant
  prime to the primes of I).
- *Pure values*: N-independent straight-line programs over A (any depth; coins fixed).
- *Gate node*: element P_v(Z) = Σ_{k≤D} c_k Z^k with pure coefficients c_k ∈ A, evaluated at the gate output
  Z = y^{E(N)} of a pure base y; element test P_v = 0 in A/NA (all coordinates, the CRT criterion of primality §0).
  Pure nodes (no gate) are PT nodes.
- *Exponent menu*: E(N) = αN + β with α ≥ 1, or E(N) = ⌊(N − c)/2^j⌋; fixed α, β, c, j; E(N) ≥ 1 on the support.
  Gate outputs are not reused (one layer); several gate nodes may share a base.
- Covers, with fixed parameters: Fermat, MR-type ladders, the power in Solovay–Strassen, Lucas/strong Lucas with a
  fixed D, quadratic Frobenius with a fixed algebra, the AKS identity tested as an element. Versions that choose
  parameters from N mod m (Selfridge's D via (D|N), the exact ladder (N−1)/2^j with j ≤ v₂(N−1), the comparison
  with (a|N)) branch on side information and need §3.3's fiber version.

**The menu restriction is necessary (codex_1 §1, PROVED).** With an arbitrary computable exponent E(N) = 1 + h(N),
the pure base −1 and the test (−1)^{E(N)} + 1 = 0, a tree reads any polynomial-time predicate h without splitting.
Within NE1 such small public-order bases are charged to the small-order term below.

### 3.2 Counting lemma and Theorem NE1. PROVED
Fix a prime x ∈ I and the coins. **Reference κ(x)** (formal transcript): a gate node answers "zero" iff, in every
component F_{x^d} of A ⊗ F_x, the component polynomial P_{v,i} becomes zero after substituting Z_i := 0 in components
where y_i = 0 (legitimate since E ≥ 1) and leaving Z_i formal elsewhere. κ(x) depends on x only. Call a component
*eligible* for v if y_i ≠ 0 and P_{v,i} ≢ 0.

**Lemma NE-C.** If v's formal answer at x is "zero", the actual answer is "zero". If it is "nonzero", there is a
component i with P_{v,i} ≢ 0; if no such component is eligible (all have y_i = 0) the actual answer equals the formal
one; otherwise, for every eligible i with o_i = ord(y_i), the actual answer is "zero" only if E(N) ≡ t (mod o_i)
for one of ≤ D residues t. Hence

  Pr_q[v deviates at x | x] ≤ D · min_{i eligible} μ_x(o_i),   μ_x(o) := max_t Pr_q[E(xq) ≡ t (mod o) | x].

*Proof.* "Zero" needs P_{v,i}(y_i^{E(N)}) = 0 in F_{x^d}; P_{v,i} has ≤ D roots; each root lying in ⟨y_i⟩ fixes
E(N) mod o_i. ∎

**Theorem NE1.** Write ν_v(x) := D_v·min_{i eligible} μ_x(o_i) (0 if no eligible component). For an NE1 tree,
Pr_q[τ_x ≠ κ(x) | x] ≤ Σ_{v ∈ path κ(x)} ν_v(x), and for independent p, q

  Pr[correct or split] ≤ B + 2s + 2(E_p[Σ_{v∈path κ(p)} ν_v(p)] + E_q[Σ_{v∈path κ(q)} ν_v(q)])

(= B + 2s + 4E_x[...] in the iid case).
*Proof.* Each c_v(N) is a fixed function of N, so the first index where τ_x and κ(x) differ is a node on κ(x)'s
path at which the answers differ; union bound with Lemma NE-C; then PT-X. ∎

*Prior art.* "Output a gcd when a test is inconsistent between components" is standard in the generic ring model
(Leander–Rupp, ASIACRYPT 2006; Aggarwal–Maurer, EUROCRYPT 2009, journal version IEEE Trans. IT 62(11) (2016)
6251–6259); PT-X/NE1 have no counterpart there. The closest analogue of Lemma NE-C is the hidden-order generic group
model (Damgård–Koprowski, EUROCRYPT 2002, LNCS 2332, 256–271; Rotem–Segev, CRYPTO 2020, LNCS 12172, 481–509), where
exponent arithmetic modulo an unknown order is analysed behind a hidden input; here N itself is the exponent and
there is no hidden input. No exact prior formulation of NE1 located.

**Corollary NE1-bal (μ uniform on the primes of I).** For E = αN + β and g = gcd(α, o): q is confined to ≤ g classes
mod o, so μ_x(o) ≤ (g|I|/o + g)/π(I) ≤ α(|I|/o + 1)/π(I). For E = ⌊(N − c)/2^j⌋: q is confined to ≤ 2^j classes
mod 2^j·o, so μ_x(o) ≤ (|I|/o + 2^j)/π(I). With SO_L(x) := "some gate node on κ(x)'s path has all eligible
component orders ≤ L" (an event of x alone):

  Pr[D] ≤ 2Pr_x[SO_L(x)] + 2h·D·α·(|I|/L + 1)/π(I)   (affine), with α(|I|/L + 1) replaced by |I|/L + 2^j (floor).

Small-order mass: for a uniform random base, Pr[ord_x(y) ≤ L] ≤ Σ_{d|x−1, d≤L}φ(d)/(x−1) ≤ Lτ(x−1)/(x−1), and
L = √X gives Pr[D] ≤ h·D·α·N^{−1/4+o(1)}; for a fixed integer base |y| ≥ 2 with log|y| = O(1), the primes with
ord_x(y) = k divide y^k − 1, so #{x ∈ I : ord_x(y) ≤ L} ≤ L² log|y|/log X, and L = X^{1/3} gives
h·D·α·N^{−1/6+o(1)}. (y ∈ {0, ±1}: order ≤ 2 at every x; such bases are entirely charged to SO.) The decay rates
assume α, 2^j, h, D bounded by N^{o(1)}.
**In words:** a polynomial-size NE1 tree learns essentially nothing without splitting, unless its pure part has
already placed a gate base in a small subgroup at p (a one-prime, pure group-order event).

### 3.3 What the deviations are, and what is not covered
Case split for a deviation of gate node v at x (o = largest eligible component order; t the residue of
Lemma NE-C). These labels are descriptive; for general P_v the residue t_x is a discrete logarithm of a root and may
depend on x arbitrarily (codex_2 #6), so the M/A reading is sharp only for targets that are fixed powers of the base
(w = y^t with fixed t, e.g. w = 1).
- **SMALL** (o ≤ L): a pure small-order event, family (M) (e.g. y = 3^{lcm(1..64)}, the p−1 method); it includes
  public-order bases such as −1 or i, which read N mod m.
- **SUBGROUP** (o > L): o | E(N) − t, i.e. gcd(E(N) − t, |G_x|) ≥ o. For Fermat at x = p this is
  ord_p(y) | gcd(N − 1, p − 1) = gcd(q − 1, p − 1): the conjectured gcd structure, family (M).
- **UNIVERSAL** (a subcase): exp(G_x) | E(N) − t, so every base is fooled. When exp(G_x) ≥ |I| this pins the
  partner to O(α) integers: family (A). Examples (all observed in E1, §9):
  - ARCH: 2^N = 2^k at the smaller prime p iff ord_p(2) | q − k (SUBGROUP in general; codex's p = 17, q = 23, k = 15);
    UNIVERSAL iff (p − 1) | q − k iff q = p + k − 1 in the balanced range.
  - Lucas α^{N+1} = 1 in the norm-one torus at both primes: (p+1) | N+1 ⇔ (p+1) | q − 1 and (q−1) | N+1 ⇔
    (q−1) | p + 1; together q = p + 2 (then N + 1 = (p+1)²), with (5|p) = −1, (5|q) = +1. A two-sided, non-split,
    archimedean liar (F₃₂₄ ≡ 0 mod 323 is the classical instance). KNOWN: twin-prime Lucas pseudoprimes are in Somer,
    Fibonacci Quart. 46/47 (2008/09) no. 3, Thm 2.1(i); Rotkiewicz, Ann. Math. Sil. 17 (2003) 17–39 (Fibonacci case,
    p ≡ 7 mod 10); Arnault, Math. Comp. 66 (1997) 869–881. (For the Fibonacci α of norm −1, α^{p+1} = −1 at an inert
    p; the conclusion still holds since (p+1)² is even. Our E1 algebra x² = 3x − 1 has norm 1.) useless, since q − p = 2 factors N by Fermat.
- **(Q) in NE1.** At p, a^{(N−1)/2} = (a|p)·a^{(q−1)/2}, so quadratic characters do enter Euler-type gates, and an
  (M)-type event can carry quadratic information (MR exhibiting √−1 certifies ε = +1, eps PC4). The theorem needs no
  separate term for them; this is not a proof that (M) and (Q) are disjoint mechanisms.

**Theorem NE1-R (codex_1 §1, "terminal probes"; PROVED there; re-checked by both reviewers).** Suppose every gate
node is a probe y^{E} = 1 whose base is fresh and uniform in the local unit groups (independent CRT components,
uniform given the history), the gate value is discarded after the test, the probe is gcd-checked, and E is specified
before the base is drawn. Then

  Pr[correct or split] ≤ B + 2s + 2R + 2/π(I),   R = E Σ_{probes reached} κ_p(E)κ_q(E)·1[p ≠ q],
  κ_r(E) = ∏_i gcd(E, d_{r,i})/d_{r,i} for G_r ≅ ∏_i C_{d_{r,i}};

the 2/π(I) term charges p = q and disappears under the virtual-diagonal convention. *Proof:* couple with the pure
shadow tree in which every probe answers "false"; the first disagreement is a split (pass at one prime) or a
both-pass event of conditional probability κ_pκ_q by freshness. ∎ Per probe the both-pass probability is exactly
κ_pκ_q. For Fermat this is gcd(N−1, p−1)gcd(N−1, q−1)/φ(N), Monier's liar count (Monier, "Evaluation and comparison
of two efficient probabilistic primality testing algorithms", TCS 12 (1980) 97–108); Lucas analogues are in
Baillie–Wagstaff, "Lucas pseudoprimes", Math. Comp. 35 (1980) 1391–1417 (theorem numbers UNVERIFIED); averages over N
in Erdős–Pomerance, Math. Comp. 46 (1986) 259–279. It is the eps PC4/MR rate and primality §1's order-pattern lemma;
R is an upper bound for the whole tree. NE1 is the deterministic-coin, arbitrary-target version.

**Fiber version and ε: OPEN as stated.** By PT-X for independent laws, conditioning on N ≡ c (mod m) and
(p mod m, q mod m) = (u, c/u) gives Pr[correct or split | c] ≤ E_{u|c}[max_b Pr[F = b | u] + 2s + 2Pr[D | fiber]];
public-order gates (orders dividing m) are then evaluated exactly in the reference, and μ is computed in
progressions (losing ≈ φ(m)). For ε on N ≡ 1 (mod 4) the per-fiber floor is already 1, so an NE1 analogue of PT2
(a cross-population split rate) would be needed. We have not proved it. So: **NE1 bounds non-split leakage of any
target under independent laws without side information; the ε question for NE1 trees that branch on N mod 4 is
OPEN.**

### 3.4 Several gates in one node (F_p only). PROVED (per node) from Kelley 2016 / Theorem SP
Let a node be P(y_1^{E_1(N)}, …, y_k^{E_k(N)}) with affine E_j and all y_j ∈ F_x^*. Since x ≡ 1 (mod x − 1),
E_j(xq) ≡ α_jq + β_j (mod x − 1), so as a function of r = q mod (x − 1) the node is f(g^r) = Σ_m γ_m g^{e_m r}, a
Laurent polynomial with T ≤ C(D + k, k) terms after merging (g a generator; identically-zero f is a formal zero).
Kelley ("Roots of sparse polynomials over a finite field", LMS J. Comput. Math. 19(A) (2016) 196–204,
Thms 2.2–2.3; = slpcirc Theorem SP) gives, for zero density ≥ λ on F_x^*:
every prime ℓ > (2/λ)^{T−1} dividing x − 1 satisfies ℓ^{v_ℓ(x−1)} | e_m − e_{m'} for a partner pair, so some monomial
ratio ∏y_j^{α_j(m_j − m'_j)} has order prime to ℓ — an (M) alignment. Three cases per node:
(i) zero density < λ: Pr_q[deviation] ≤ λ(x − 1)(|I|/(x−1) + 1)/π(I) ≤ (3 + o(1))λ ln X;
(ii) density ≥ λ with alignment (M). (v3, literature critic: Kelley's Thm 2.2 gives every exponent a partner
congruent modulo the large zero-coset size C ≥ (λ/2)^{T−1}(x−1), whether or not x − 1 is smooth; v2's separate
"smooth" case (iii) was too pessimistic and is dropped: the alignment then holds modulo C rather than prime by prime.) This is a per-node classification; O(1/ln X) per node does
**not** make a polynomial-height tree negligible, so §3.4 is not part of the quantitative statement in §8.

### 3.5 Exponents larger than the partner's range
In F_{x^d}, d ≥ 2, E(N) mod o can exceed |I|; Lemma NE-C still covers element tests (μ ≤ 2α/π(I)). **Coordinate**
tests (one coefficient of an element of Z/N[X]/(X^r − 1), as in AKS) are sparse exponential sums over Galois
conjugates along the short progression {xq}: incomplete sums with no general bound. primality §3 measured such a
coefficient leak confined to the Kummer regime q − p ≲ r, an (A) leak, consistent with the picture but not covered.

### 3.6 Nested and retained gates: OPEN (the class-G frontier)
For (a^N + 1)^N, or nodes adding retained gate outputs of different bases, p's transcript depends on
r = q mod (p − 1) through exponential-polynomial towers r ↦ (g^r + 1)^r. A non-(M), non-(A) leak would be such a
node whose zero set in r is large (≥ |I|^{1−o(1)}/poly among the residues the partner can take) without order
alignment or a short-interval pin. Bounding such sets is the problem of counting solutions of exponential
congruences (g^r ≡ r, r^r ≡ c mod p, …). Single equations of this type have nontrivial upper bounds
(Glebsky–Shparlinski, Des. Codes Cryptogr. 56 (2010), Thm 2; Balog–Broughan–Shparlinski, Acta Arith. 148 (2011);
Cilleruelo–Garaev, Proc. AMS 144 (2016), Thm 1 / Lemma 5); towers such as (g^r + 1)^r are not covered.

## 4. Integer representatives

### 4.1 Prop REP-C (completeness). PROVED (codex_1 §2)
- A tree that can compare rep(x) ∈ [0, N) with constants reads rep(−1) = N − 1 in n comparisons and can then
  label its leaves by any function of N. PT is size-free, so **no size-free PT-type statement exists with
  representatives**.
- Even b_j = LSB(2^j mod N) = ⌊2^j/N⌋ mod 2, j ≤ 2n + 2, determine N (binary digits of 1/N; all values are units).
- Poly-size uniform trees with representatives simulate every polynomial-time algorithm (Boolean gates as
  x ∧ y = xy on {0,1}-residues). A classification of leaks here is a classification of all factoring algorithms
  (cf. selfred SR3).
- All three families live here: QS/NFS (smoothness of integer representatives), Lehman/Harvey (integer
  arithmetic on kN), SQUFOF (forms), Coppersmith (lattices). Without representatives (M) is only its group-order part
  (ECM is pure); relation (M) needs representatives.

### 4.2 Lemma GAR (one hidden phase). PROVED (identity)
For x ∈ Z/N with components x_p ∈ [0, p), x_q ∈ [0, q), and θ := (q^{−1} mod p)/p,

  rep(x)/N ≡ (x_p − x_q)·θ + x_q/N  (mod 1),   θ = j/q + 1/N, where q·(q^{−1} mod p) = 1 + jp, 0 < j < q.

*Proof.* With q̄ = q^{−1} mod p, p̄ = p^{−1} mod q, the integer qq̄ + pp̄ is ≡ 1 mod p and mod q and lies in (1, 2N),
so it equals N + 1 and q̄/p + p̄/q = 1 + 1/N. CRT gives rep(x) ≡ x_p qq̄ + x_q pp̄ (mod N); divide by N and substitute
p̄/q = 1 + 1/N − θ. The second claim is qq̄ = 1 + jp divided by N. ∎ (Checked on 1000 random exact cases by the critic.)

Consequences (sufficient conditions only):
- θ exactly gives gcd(θN − 1, N) = p.
- θ to k top bits gives θN = 1 + jp up to O(N2^{−k} + 1): an approximate multiple of p next to the exact multiple N,
  so 3n/4 + O(1) bits factor N by the partial approximate common divisor method (Howgrave-Graham, "Approximate
  integer common divisors", CaLC 2001, LNCS 2146, 51–66; the δ = 1 form with error N^{β²} in May, "Using
  LLL-reduction for solving RSA and factorization problems", in The LLL Algorithm, Springer 2010; theorem number
  UNVERIFIED). Settled by the literature critic (v2 had marked the endpoint UNVERIFIED). This is the same N^{1/4} error barrier as n/4 known bits of p (topblock §2.3), and it is
  not a lower bound.
- GAR is an exact parametrization by factor-dependent values. It is **not** a proof that useful readouts must reveal
  an approximation to θ (codex_2 #11); readouts fed back into ring arithmetic are complete (4.1).

### 4.3 The archimedean algorithms as covering events (map, not a theorem)
- Lehman/Harvey/our DCL use (aq + bp)² − 4abN = (aq − bp)² (codex): a small |aq − bp| puts a factor-dependent
  integer in a short searchable interval; success = "q/p is covered by a rational of small height". coverage38
  Theorem U bounds any covering design in that family by N^{1/5}(λ⁴/J)^{1/5}; dchworst Theorem R is the recovery
  step.
- Coppersmith: success = "p within N^{1/4} of an enumerated approximation"; topblock Prop T: window lattices add
  nothing to the univariate lattice.
- Elementary size leakage is also (A): the top bits of N predict min(p, q) > 3X/2 with constant advantage, without
  splitting (codex_1 §1). "Family A is always exponential" is true only for targets that are not functions of the
  factor sizes.
- In E2 the target q/p < 33/32 is predicted by CF square forms, and post hoc only on instances those square forms
  split (§5.2).
- **OPEN and unformalisable as stated:** "every leak in the representative model is Diophantine".

## 5. Class H: continued fraction of √N and infrastructure

### 5.1 What is proved (classical; codex_1 §3 gives the proofs)
- A_{i−1}² − NB_{i−1}² = (−1)^iQ_i, hence ((−1)^iQ_i | p) = ((−1)^iQ_i | q) = 1 when gcd(Q_i, N) = 1 (genus).
  Odd i with Q_i = k² certifies ε = +1; even i with Q_i = ℓk² certifies (ℓ|p) = +1; odd i with Q_i = ℓk² certifies
  (−ℓ|p) = +1 (E2 controls: 0 violations in 8·10⁵ instances).
- Products of Q_i that are squares give z² = ±1 mod N: a √−1 (ε certificate) or a congruence of squares (split,
  CFRAC); square forms give SQUFOF and Lehman-type splits.
- Statistical models of the principal cycle are heuristic: Gower–Wagstaff, "Square form factorization",
  Math. Comp. 77 (2008) 551–588, use equidistribution assumptions (§4, Assumptions 4.8–4.23); Assumption 4.13 is the
  genus-independence heuristic closest to H-GEN; "subject only to genus constraints" is literally false (period, reversal symmetry; codex).

### 5.2 E2 (pre-registered): genus targets leak, non-genus targets null; archimedean leak via splits
16 CF statistics (src/eps/cf.py menu) × 7 targets, 2·10⁵ instances per size, fresh-seed replication. **v2
statistics:** the permutation z now includes the fixed contribution of singleton strata (a v1 bug found by codex_2:
it inflated some n = 40 A1 scores, e.g. cf1024_n_Qeven 7.35 → 0.69), and the cluster z is a *dyadic-robust* z
(observations dependent when they share a prime in either role; v1 clustered smaller and larger factors separately).
Primary z: n = 32 dyadic-robust (prime reuse), n = 40 permutation (both reported in the json).

| target | type | n=32 hits (confirmed) | max |z| n=32 | n=40 hits (conf.) | max |z| n=40 | MDE |ρ| n=32 / n=40 (design-effect adjusted) |
|---|---|---|---|---|---|---|
| G1 ε = (−1|p) | genus | 9 (9) | 14.4 | 4 (4) | 8.7 | 0.15 / 0.037 |
| G2 (−3|p) | genus | 3 (3) | 12.2 | 1 (1) | 5.8 | 0.15 / 0.037 |
| G3 (5|p) | genus | 1 (1) | 8.4 | 0 | 2.7 | 0.16 / 0.041 |
| G4 (−7|p) | genus | 1 (1) | 5.8 | 0 | 3.0 | 0.14 / 0.045 |
| NG1 [p ≡ ±1 mod 7] (cubic) | non-genus | 0 | 3.3 | 0 | 2.6 | 0.15 / 0.045 |
| NG2 [p ≡ 1 mod 5 \| p ≡ ±1] (quartic) | non-genus | 0 | 1.7 | 0 | 2.3 | 0.22 / 0.057 |
| A1 q/p < 33/32 | archimedean | 6 (6) | 16.8 | 2 (2) | 34.1 | 0.033 / 0.015 |

The design effect at n = 32 is large (up to 65) because 16-bit primes are reused ~270 times each; the n = 32 nulls
are therefore weak (MDE ≈ 0.15), the n = 40 nulls moderate (≈ 0.05). G4 and NG1 share a population (N ≡ 1 mod 7):
same rows, the quadratic target hits and the cubic one does not.

**Post hoc (not pre-registered; `e2_posthoc.py`).** Square forms Q_i = r² at even i ≤ 64 occur in 28.7% (n=32) /
8.1% (n=40) of instances; 56–58% of them split N via gcd(A_{i−1} − r, N); P(A1) is 0.42 on split instances and
0.15–0.19 otherwise. On the instances with no splitting square form, all 16 statistics have |z| ≤ 2.9 (n=32) and
≤ 1.6 (n=40) against A1 (dyadic-robust: ≤ 2.7, ≤ 1.6). Qualified reading: at these sizes and detection limits, the
archimedean signal in the CF prefix is carried by Lehman-type square-form splits; a residual leak below the
detection limit is not excluded.

**Conjecture H-GEN (refined).** For poly(n)-length CF/infrastructure prefixes of √(kN), k fixed, every
"low-complexity symmetric statistic" of the (a_i, Q_i) predicts a function of the factors, beyond the fiber floor
and beyond events where the prefix's square relations split N, only through quadratic symbols (ℓ|p) with small ℓ
(incl. −1, 2), via ((−1)^iQ_i | p) = 1. Status: CONJECTURED; E2 is consistent with it (registered non-quadratic
targets null at the MDEs above; quadratic targets hit). Without a restriction on post-processing the statement is
vacuous (Q₁ = N − P₁² encodes N; codex); the statistic class is not formalised. E2 tests two non-quadratic
characters only, on their registered populations; it does not show independence from all non-genus information.

### 5.3 Why this is (Q) and not new
The leaked symbols are fixed-base quadratic residuosity bits (nonsmooth Prop RANK; Seurin's 2-Φ/4-hiding for −1).
Their strength decays with N because small Q_i become rare (eps §6); a constant-advantage symbol costs SQUFOF-type
N^{1/4} or CFRAC-type L[1/2] work (eps §6, PC2/PC3).

## 6. Class F: AC⁰ and automata

### 6.1 Theorem AC0-ε (codex_1 §4). PROVED; checked in codex_2 #12 and by the Opus proof critic
For iid uniform primes p, q ∈ (X, 2X] conditioned on N ≡ 1 (mod 4), every family of depth-d, size-n^{O(1)} circuits
on the bits of N has Pr[A(N) = ε(p)] = ½ + O_{C,d}(n^{−C}) for every C.

*Proof* (codex_1 §4, steps restored):
1. *Character pair.* χ₄(p)·1[pq ≡ 1 (4)] = (χ₄(p) + χ₄(q))/2, so Σ_{pq≡1(4)}χ₄(p)K(pq) = Σ_{p,q}χ₄(p)K(pq) for any K.
2. *Walsh to additive characters.* W_S(n) = (−1)^{Σ_{i∈S}bit_i(n)} on b-bit integers has additive Fourier ℓ¹-norm
   ≤ (Cb)^{|S|} (each bit's square wave has ℓ¹-norm O(b); products convolve; Bourgain, "Möbius–Walsh correlation bounds and an
   estimate of Mauduit and Rivat", J. Anal. Math. 119 (2013) 147–163, Lemma 1).
3. *Bilinear prime sum.* Uniformly in α, Σ_{p,q∈(X,2X]}χ₄(p)e(αpq) ≪ X²e^{−c(log X)^{1/4}}. Put L = log X,
   Q = e^{L^{1/4}}, choose a/r with r ≤ X²/Q, |α − a/r| ≤ Q/(rX²).
   - *Minor arcs* (Q ≤ r ≤ X²/Q): for arbitrary bounded coefficients,
     |Σ_{m,n}a_mb_ne(αmn)|² ≪ X⁴/r + X³log(2r) + X²r log(2r) (Cauchy–Schwarz, then Σ_{h≤X}min(X, ‖αh‖^{−1})).
   - *Major arcs* (r < Q): with D = lcm(4, r) and p, q > D, e(at/r) = Σ_{χ mod D}c_χχ(t), |c_χ| ≤ 1, so each
     rectangular sum Σ_{p≤u, q≤v} is a sum of c_χ(Σ_pχχ₄(p))(Σ_qχ(q)). By Landau–Page (Ford–Green–Konyagin–Maynard–Tao,
     "Long gaps between primes", JAMS 31 (2018) 65–105, Lemma 7.1; the prime-character estimate in Watkins,
     "Distribution of the 2-Selmer rank under twisting", Publ. Math. Besançon 2022, 59–133, §3.4.1; textbook form in
     Davenport, Multiplicative Number Theory, ch. 20), all
     non-principal characters of conductor ≤ 4Q except possibly one real ξ have prime sums ≪ Xe^{−c√L} uniformly in
     the interval. χ and χχ₄ cannot both lie in {principal, ξ}: that would force ξ = χ₄, but L(s, χ₄) > 0 on (0, 1).
     Two-variable partial summation with H(u, v) = e(βuv), H_{uv} = (2πiβ + (2πiβ)²uv)H, costs
     O(1 + |β|X² + |β|²X⁴) = O(Q²).
4. *Normalise* by #pairs ≍ X²/L²: β_t := max_{|S|≤t}|E[εW_S(N) | N ≡ 1 (4)]| ≤ exp(−cb^{1/4} + O(t log b)).
5. *LMN* (Linial–Mansour–Nisan, JACM 40 (1993)): the Fourier tail of A above degree t is ≤ 2s·2^{−t^{1/d}/20}. The
   conditioned semiprime law has density O(b²) w.r.t. the uniform cube (each N has ≤ 2 representations), so
   |E[εA]| ≤ √(Σ_{j≤t}C(b,j))·β_t + O(b·tail^{1/2}). Take t = (C′log b)^d. ∎

**Scope.**
- Target ε. For another fixed real character χ(p) the same proof works for large X because a fixed character is not
  exceptional; the L(s, χ₄) > 0 step is specific to χ₄.
- For general targets the statement is false (factor-size leakage from the top bits of N, §4.3).
- Finite automata: Thue–Morse and Rudin–Shapiro follow from Mauduit–Rivat (JEMS 17 (2015), Prop. 2; number UNVERIFIED) and Müllner
  ("Automatic sequences fulfill the Sarnak conjecture", Duke 166 (2017), Prop. 5.5, verified), which allow arbitrary bounded coefficients on balanced ranges (codex; hypotheses
  checked only at the level of statements). General automata, AC⁰[⊕] and AC⁰ on bits of pure residues: OPEN.
- Prior art: Green (CPC 21 (2012)) and Bourgain (2013) treat Möbius–Walsh/AC⁰; AC⁰ lower bounds for other
  number-theoretic targets: Allender–Saks–Shparlinski, JCSS 62 (2001); Bernasconi–Damm–Shparlinski, Inf. Comput. 168
  (2001); Shparlinski, ECCC TR98-054; Shparlinski's book (1999) and Kopparty
  (STOC 2011) treat Boolean complexity of residuosity/powering in other settings. No exact predecessor for the hidden
  common sign of balanced semiprime factors was located (UNVERIFIED as a negative). This closes the AC⁰ entry of
  selfred §3.2 for ε.

## 7. Jacobi symbols of pure values (Q)

**Theorem PT-J. PROVED.** Population: e ∈ {±1} uniform, p, q iid from μ_e (primes of I with ε = e), so N ≡ 1 (mod 4).
Observations: J_i = (c_i | N) = λ_i(p)λ_i(q), i ≤ k, for N-independent integers c_i that are units at every prime of
I (otherwise a gcd splits). For S ⊆ [k] put c_S = ∏_{i∈S}c_i, m_e(S) = E_{μ_e}[(c_S | x)]. Then

  Pr[correct] ≤ ½ + ¼ Σ_{∅≠S⊆[k]} |m_+(S)² − m_−(S)²|,

with equality for k = 1. If μ gives ε = ±1 equal mass (an artificial equal-prior population; for the uniform prime
law the masses are only approximately equal and the general form applies), m_±(S) = β(c_S) ± β(−c_S),
β(c) := E_μ[(c|x)], and the bound is ½ + Σ_{S≠∅}|β(c_S)β(−c_S)|. A single test leaks iff β(c)β(−c) ≠ 0.
*Proof.* The law of u = λ(p)⊙λ(q) ∈ {±1}^k under e has Fourier coefficients m_e(S)² (independence); Bayes success is
½ + ¼Σ_u|P_+ − P_−| ≤ ½ + ¼Σ_S|P̂_+ − P̂_−|; and (−c|x) = ε(x)(c|x). ∎ (Unit-tested by brute force.)

**Corollary PT-J1.**
- (a) Unconditional. If Σ_i log|c_i| ≤ (log X)^{1−δ} (so every c_S qualifies), then
  |β(c_S)β(−c_S)| ≤ e^{−c(log X)^{δ'}}: (c_S|·) and (−c_S|·) = χ₄(c_S|·) are characters of conductor ≤ 4|c_S| whose
  ratio is χ₄, and Landau–Page with T = exp((log X)^{1−δ}) leaves at most one exceptional character, which cannot be
  both (§6 step 3). At most one of ±c_S is a square.
- (b) GRH. With the standard bound ψ(x, χ) ≪ √x·log x·log(fx) (Iwaniec–Kowalski Thm 5.15 / Montgomery–Vaughan
  (13.21); numbers UNVERIFIED), |β(c)| ≪ X^{−1/2}log X·log(|c|X) for non-principal (c|·); so for Σ_i log|c_i| ≤ X^{1/4−δ} the leak of
  k = O(log n) tests is ≤ 2^k X^{−1/2−2δ}·polylog (the range extends towards X^{1/2−δ}; we state the safe range).
  (v1 printed a log³ bound that did not give its own range; corrected.)
- (c) Without restriction the criterion is met: codex's CRT-constructed c (a non-residue exactly on half of the
  ε = +1 primes) gives success → ¾; but such c has no small SLP.

**OPEN:** SLP-compressed c with log|c| ≫ X^{1/4} (e.g. 2^{2^k} + 1 with 2^k > X^{1/4}): can a poly-size SLP make both
(c|x) and (−c|x) biased over primes? **E3 (pre-registered):** 459 menu values; 114 excluded as structural squares
(V_{2m}(3) + 2 = V_m², and 88 random SLPs ending in a squaring); 345 tested over all 513 708 primes in [2^23, 2^24):
**0 hits** (max min(|z₁|, |z₂|) = 2.27, max |z₁z₂| = 5.3).

## 8. Synthesis: what is classified and what is open

**Summary theorem (a list of separate quantitative statements; items 1 and 5 are imported).**
1. Pure trees, side information N mod m: any target, iid/independent laws (factorbits PT, PT-m, PT2).
2. NE1 trees (fixed affine/floor exponents, one gate per node, element tests, any pure depth, adaptive, no side
   information): any target F under independent laws; non-split advantage ≤ 2Pr[D] ≤ 4Pr[SO_L] +
   O(h·D·α(|I|/L + 1 + 2^j)/π(I)) (Theorem NE1, Cor NE1-bal); with fresh random probes ≤ 2R + 2/π(I) (NE1-R).
3. A non-adaptive layer of k Jacobi tests of pure integers: target ε; negligible for Σ log|c_i| ≤ (log X)^{1−δ}
   (unconditional) or ≤ X^{1/4−δ} (GRH) (Theorem PT-J, Cor PT-J1).
4. AC⁰ on the bits of N: target ε (or a fixed real character of p); advantage O(n^{−C}) (Theorem AC0-ε).
5. Functions of N mod m: ε or f(p mod m); m ≤ X^{1−η} under GRH, or 2^j ≤ N^c unconditionally via Bourgain's
   prescribed digits (selfred §3.1).
In 2 the leak is charged to pure small-order events and to the residue coincidence E(N) ≡ t mod ord, whose
interpretation is (M), degenerating to (A) pins; in 1, 3–5 non-split leaks are negligible for the stated targets.
These are separate models; unions (Jacobi tests after an NE1 tree, AC⁰ on bits of pure residues, NE1 with N mod 4)
are not covered, and representatives (§4) make any union containing them complete.

**Strongest natural open model: PRIM.** Uniform programs of N-independent shape over fixed Z/N-algebras whose only
N-dependent primitives are exponent gates x ↦ x^{E(N)} from a fixed polynomial menu (E(N) ∈ {N, N ± 1, (N−1)/2, …})
**with retained, nested and summed gate outputs**, plus gcd/zero tests and Jacobi symbols of computed values. (The
gates must be primitives: computing x^N by square-and-multiply would be an N-dependent SLP, which is vacuous,
slpcirc Prop R.) PRIM contains the fixed-parameter forms of the standard primality tests. Our results cover its
one-layer part and its terminal Jacobi part with moderate conductors. Relevant generic-model antecedents, none of
which classifies PRIM: Aggarwal–Maurer (EUROCRYPT 2009; IEEE Trans. IT 62(11) (2016) 6251–6259; generic RSA with E
an input), Jager–Schwenk (J. Cryptology 26 (2013); generic Jacobi/QR hardness — the reverse direction of PT-J, not a
predecessor), Boneh–Venkatesan (EUROCRYPT 1998), Maurer–Wolf (EUROCRYPT 1998). Whether some existing
theorem settles PRIM is UNVERIFIED (we found none).

**Where a breakthrough would have to live.** By PT-X, a non-splitting leak in PRIM is a deviation of one prime's
transcript from a one-prime reference. In NE1 these are residue coincidences with an element order (cosets of
subgroups of Z/(p−1)) or exact pins. A new family would be a retained-gate node whose zero set in r = q mod (p − 1)
is large and unexplained by order alignment or short intervals (§3.6), or a compressed c with β(c)β(−c) ≠ 0 (§7).
Neither is excluded by anything we can prove; E1 (no retained gates) and E3 (compressed c) give no hint of one.

## 9. Experiments (pre-registered in `notes/leakmap_prereg.md`, commit d050db1)

### E1: instrumented N-exponent trees (checks of proved statements)
Populations: all pairs of b-bit primes with the top two bits set, b = 12, 14 (iid); 2·10⁵ sampled pairs at b = 16;
the Fermat-regime pairs q − p ≤ 2^{b/2}. 7 trees. 42 cells, 1.9·10⁶ pair-runs (× 3 transcripts each).
- **P1 (diagonal lemma): 0 mismatches** in every cell; every split is STRUCTURAL (τ_p ≠ τ_q), ACCIDENTAL 0.
  The registered MR ladder (exact exponents (N−1)/2^j, j ≤ v₂(N−1), a per-N tree) was added in v2
  (`e1_mrl.py`): 0 mismatches in all 6 populations.
- **P2 (PT-X): holds in all 126 (cell × target) evaluations** (Bayes-optimal leaf labels). Largest non-split gain:
  pure p−1 trees (≈ 0.02–0.03, PT's smoothness/residue-class correlation); NE trees ≤ 0.004 except MR2 on the b = 12
  Fermat population (0.0077).
- **P3 (NE1 counting bound, registered form with Σ over components, which is weaker than the min form): 0
  violations**; worst ratio freq/bound 0.74.
- **P4 (random bases):** b = 14, 2·10⁴ pairs × 20 bases: Fermat liars 7 vs 9.7 expected, MR strong liars 4 vs 3.9
  (underpowered as registered); added b = 10, all 666 pairs × 2000 bases: 359 vs 344.4 (z = 0.79), 171 vs 159.0
  (z = 0.95). Non-unit bases (probability ≈ 2/√N) are skipped and not corrected for. Passed.
- **P5:** "UNIVERSAL deviations occur only for ARCH in iid populations" **refuted as registered** (a prediction
  error): LUCAS has two-sided UNIVERSAL deviations at twin primes and MR2's floor exponents give pins too. The exact
  ARCH pin count holds in all 6 cells (432, 1307, 682, 432, 1307, 3870). Directional prediction holds: the UNIVERSAL
  share is larger in the Fermat population for every tree and size where UNIVERSAL occurs (LUCAS 0.21–0.50 →
  0.95–0.96; ARCH 0.33–0.50 → 0.96–0.98; MR2 0.11–0.19 → 0.26–0.39).
- Deviation classes (v2 uses the largest eligible order for SMALL, as in the theorem): e.g. n = 32 iid ARCH:
  SMALL 584, SUBGROUP 804, UNIVERSAL 680; LUCAS non-split deviations are all UNIVERSAL twin-prime pins.

### E2 — §5.2. Q1 (G1 replicates): yes, 9/9. Q2 (G2–G4 reach): yes, 5 confirmed. **Q3 (decisive): NG1, NG2 0 hits
at both sizes, max |z| 3.3 < 4 (fresh 3.55).** Q4 (decay): 14 genus hits at n = 32, 5 at n = 40. Controls PC{3,5,7}:
0 violations.

### E3 — §7. R1 (0 hits): held.

### Deviations
- D1: P4's b = 14 cell was underpowered; an extra b = 10 cell was added after the run (labelled).
- D2: E2 post-hoc split analysis of A1 (not registered).
- D3: P5's first sentence was a mis-specified prediction; reported as refuted.
- D4 (codex_2 #19): E1's `MR2` is a fixed battery of floor exponents ⌊(N−1)/2^j⌋, j ≤ 3, not the registered
  v₂(N−1)-ladder (v2 adds a P1 run of the registered ladder; P2/P3 were not run for it because its node set depends on
  N mod 2^k, i.e. it needs the fiber version). E1's targets ε and [p ≡ 1 mod 3] were evaluated on all pairs (the iid
  law PT-X needs), not on the registered sub-populations N ≡ 1 (mod 4) / (mod 3).
- D5 (codex_2 #16–17): E2 statistics corrected in v2 (singleton strata; dyadic-robust variance); all E2 numbers
  above are v2. Changes from v1: A1 n = 40 hits 4 → 2; n = 32 genus maxima ≈ 10% lower; conclusions unchanged.
- D6: v2 deviation classification uses the largest eligible component order (v1: smallest).
- No registered statistic or cell was dropped.

## 10. Candidate cards (rule 9)

| Mechanism | Model | Outcome |
|---|---|---|
| Fermat/MR/Lucas/Frobenius exponents N ± 1 as a factor leak | NE1 | **refuted for this model** beyond small-order and residue-coincidence events (Theorem NE1, NE1-R; E1) |
| Archimedean exponent probes 2^N = 2^k | NE1 | (A) pins q − p = k − 1 exactly as predicted; cost linear in the window |
| Arbitrary exponent selectors E(N) | NE with free E | predicate-complete (codex); not a mechanism |
| Retained/nested N-th powers | PRIM | **OPEN** (frontier) |
| Integer-representative phase θ | one-shot readouts | identity only; (3/4 + ε)n bits of θ suffice (HG) |
| CF prefix statistics | H | quadratic symbols leak, two non-quadratic characters null (E2); (A) only via splits (post hoc) |
| AC⁰ on bits of N, target ε | F | **refuted** (Theorem AC0-ε, unconditional) |
| Jacobi of pure values | PT-J | refuted for Σ log|c_i| ≤ X^{1/4−δ} (GRH); open for compressed c; E3 null |

## 11. Reviews

**Codex adversarial round (`codex/leakmap_2.out`, MAJOR REVISION; 20 findings), applied in v2:**
- OK: PT-X (#1), NE-C core (#3), NE1-R (#8), GAR identity and the 3/4 exponent (#10), AC0-ε (#12), PT-J (#14).
- WRONG, fixed: formal reference with zero bases and E = 0 (#2: E ≥ 1 required, eligible components defined);
  fixed-base small-order estimate for y ∈ {0, ±1} (#5); the ARCH equivalence (#6: now ord_p(2) | q − k, UNIVERSAL iff
  q = p + k − 1); E2 singleton-strata bug (#16, rerun); E1 implementation vs registration (#19, disclosed as D4,
  registered ladder P1 added).
- GAP/OVERCLAIM, fixed by narrowing: NE1-bal qualifications and the dyadic exact-division case (#4); M/A labels as
  interpretation, not a mechanism theorem (#6); "(Q) does not occur" withdrawn (#7); §3.4 is per-node only (#9);
  GAR is an identity, not a classification (#11); prior art for AC⁰ (#13); PT-J1(b) GRH range (#15); cluster
  variance and power (#17: dyadic-robust z and design-effect MDE); E2 non-quadratic inference narrowed (#18);
  CLASS rewritten as a list with targets and laws, PRIM defined with primitive gates, generic-model literature added
  (#20).

**Opus proof critic (ACCEPT WITH FIXES; 16 findings), applied in v2:** PT-X for independent laws and p = q (#1);
min over eligible components, α ≥ 1, no-eligible-component case (#2); NE1-bal hypotheses (μ uniform on primes,
|y| ≥ 2, floor form, |α| in §0) (#3); the N mod m gap — fiber version stated, ε under NE1 marked OPEN (#4); Lucas
wording, exp(G_x) ≥ |I|, "(Q) only as fiber information" (#5); NE1-R wording (#6); §3.4 constant, ℓ^{v} statement,
smooth case, removed from the quantitative summary (#7, #14); GAR "sufficient" (#8); AC⁰ scope for other characters,
TM/RS label (#9); equal-mass caveat in PT-J (#10); PT-J1(b) range (WRONG in v1) corrected (#11); PT-J1(a) on c_S
(#12); CLASS targets restricted and the size-leak counterexample added (#13); imported items labelled (#15); SMALL
wording (#16).

**Opus literature critic (OK WITH FIXES), applied in v3.** No miscitation changes a theorem.
- MISSING PRIOR ART, now credited: Monier (TCS 12 (1980)), Baillie–Wagstaff (Math. Comp. 35 (1980)), Erdős–Pomerance
  (Math. Comp. 46 (1986)) for the NE1-R per-probe rate; hidden-order generic groups (Damgård–Koprowski 2002;
  Rotem–Segev 2020) as the closest analogue of NE-C (NE1's N-as-exponent setting has no hidden input); twin-prime
  Lucas liars (Somer 2008/09; Rotkiewicz 2003; Arnault 1997) with the norm −1 remark; AC⁰ background
  (Allender–Saks–Shparlinski; Bernasconi–Damm–Shparlinski; Shparlinski ECCC TR98-054).
- SETTLED: the 3n/4 + O(1) endpoint in §4.2 (May 2010, δ = 1 form; Howgrave-Graham 2001); §3.6 exponential
  congruences (Glebsky–Shparlinski 2010; Balog–Broughan–Shparlinski 2011; Cilleruelo–Garaev 2016).
- CORRECTIONS: §3.4 case (iii) dropped (Kelley's partner congruences hold modulo C regardless of smoothness);
  Aggarwal–Maurer journal version; Jager–Schwenk is the reverse direction of PT-J, not a predecessor; full Watkins and
  FGKMT references, Davenport ch. 20; Bourgain Lemma 1 and Müllner Prop. 5.5 verified; Gower–Wagstaff described as
  equidistribution assumptions (Assumptions 4.8–4.23; 4.13 for genus independence).
- Still UNVERIFIED (numbering only): Montgomery–Vaughan (13.21), Iwaniec–Kowalski Thm 5.15, Mauduit–Rivat JEMS
  Prop. 2, Baillie–Wagstaff theorem numbers, May (2010) theorem number.

## 12. Reproduction (project root)
```
.venv/bin/python -m pytest -q tests/test_leakmap.py                     # 20 tests, ~1 s
.venv/bin/python -m src.leakmap.e1                                      # E1 -> experiments/raw/leakmap_e1*.jsonl (~4 min)
.venv/bin/python -m src.leakmap.e1_p5                                   # P5 -> experiments/analysis/leakmap_e1_p5.json
.venv/bin/python -m src.leakmap.e1_mrl                                  # registered MR ladder P1 (v2)
.venv/bin/python -m src.leakmap.cfgenus run --n 32 --range main         # E2 data (~1 min, 6 workers); also n 40, --range fresh
.venv/bin/python -m src.leakmap.cfgenus analyse experiments/raw/leakmap_e2_main_n32.jsonl.gz --cluster --out experiments/analysis/leakmap_e2_main_n32.json
.venv/bin/python -m src.leakmap.e2_summary                              # confirmed-hit table
.venv/bin/python -m src.leakmap.e2_posthoc experiments/raw/leakmap_e2_main_n32.jsonl.gz experiments/raw/leakmap_e2_main_n40.jsonl.gz
.venv/bin/python -m src.leakmap.jacobi --bits 24                        # E3 (~1 min)
```
Peak RSS < 400 MB for every job.
