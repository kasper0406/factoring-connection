# prim: classifying non-split leaks in PRIM (phase 25)

Agent `prim`, 2026-10-23. **Version v1** (before the adversarial and critic rounds; §11 will record them).
- Pre-registration: `notes/prim_prereg.md` (commit `f8952d0`, before any code). Second-stage predictions for the
  post-hoc classes: `notes/prim_posthoc_prereg.md` (commit `59906fb`, before the b = 16 post-hoc run).
- Code: `src/prim/{prog,zerosets,e1,e1_posthoc,theory}.py`; tests `tests/test_prim.py`.
- Data: `experiments/raw/prim_e1A_b{12,14,16}.jsonl`; summaries `experiments/analysis/prim_e1_b{12,14,16}.json`,
  `prim_e1_posthoc_b{14,16}.json`.
- GPT-6.1-sol blind round: `codex/prim_1.{md,out}` (several results below are codex's and are credited as such).
- Labels: **PROVED** (proof here or in the cited codex file), **KNOWN** (reference), **HEURISTIC**, **CONJECTURED**,
  **OPEN**, **experimentally supported**, **refuted for this model**, **UNVERIFIED**. "No prior formulation located
  by us" is the strongest novelty statement made.

Notation as in leakmap: N = pq, p ≠ q primes of I = (X, 2X] (for theorems; experiments use b-bit primes with the top
two bits set), iid laws unless stated, n = n_p = p − 1, s = Pr[split], B = max_b Pr[F = b]. A node's *answer* at p is
[node ≡ 0 mod p]; τ_p is p's diagonal transcript. PT-X (leakmap §2): Pr[correct or split] ≤ B + 2s + 2Pr[D ∧ no split]
for any one-prime reference κ.

---

## 0. Verdict (plain language)

**No leak was found, and a large part of PRIM is now classified. Full PRIM is not.**

1. **Everything reduces to one prime and a hidden exponent (Theorem U, PROVED).** For affine gates, at p every gate
   x^{E(N)} equals x^{E evaluated at r}, r = q mod (p − 1) = 1 + (q − p). Each residue r is hit by at most two partner
   primes. So the chance that a node misbehaves at p is at most (2/π(I)) times the number of "bad" exponents r. That
   count, the *uniform-r zero set*, depends on p alone. Non-split leakage is therefore a one-prime counting question
   about programs over F_p with a uniformly random hidden exponent. The reduction fails in extension components F_{p^k},
   k ≥ 2, where the partner is a short interval.
2. **With fresh random bases, PRIM is classified at every depth for chains (Theorems R + CB, PROVED).**
   - A new decoupling theorem R says: for any program with fresh coins, Pr[correct or split] ≤ B + 4s plus "modal
     deviations". Information costs splitting, up to the event that a node's typical answer at p changes.
   - For chains of gates and shifts, a node is typically zero only if ∏ gcd(E_i(N), p − 1) ≥ (p − 1)/2h. That is an
     (M) alignment, and it is rare. The proof weakens with depth k (loss ≈ X^{−1/k}) and becomes vacuous near
     k ≈ log X / log log X.
   - Chains of pure N-th-power gates are permutations of F_p at any depth (codex).
3. **Sums of one-layer gates are (M) (Theorem PS, PROVED via a Kelley-type structure lemma, codex's compression form).**
   Zeros lie either in zero cosets of large subgroups, where partner monomials are order-aligned, or in a residual of
   size ≤ 2n^{1−(1−θ)/(t−1)}. The bound is useful for t up to about log X / log log X terms.
4. **Deterministic two-layer gates are mostly (M) (Theorem D2, PROVED, with a gap).** For nodes A·(y^{ar+b} + c)^{αr+β} = B
   (A, B monomials), a divisor-and-Jacobi-sum argument gives
   |zeros| ≤ τ(d)·gcd(α,d)·(n/d + d(n/o)(√p + 2)) + n/o, for every divisor d of p − 1, where o = ord_p(y^a).
   That is a power saving at primes where p − 1 has a divisor in [p^δ, p^{1/4}] and y has large order.
   - The excluded primes have rough p − 1 or small ord_p(y), both (M)-type properties of p.
   - Under GRH, the excluded primes have density ≤ Cδ + o(1). So per node the deviation rate is o(1), **not negligible**.
   - Codex proved a clean √p-type bound for the target w = 1. General targets are a "Zech-hyperbola" counting
     problem (ind u · ind(u + c) ≡ t mod p − 1) with no known theorem.
5. **Depth ≥ 3 with two additive shifts, deterministic bases: OPEN.** There the needed equidistribution is "doubly
   Zech" and no character-sum structure is left. This is the first depth where our methods break.
6. **Jacobi of computed values.**
   - Fresh-base gate polynomials are classified (Theorem J1, PROVED; the Weil core is codex's). The only biases are
     (M) gcd alignments, or a degenerate square shape that reduces to PT-J for a pure integer, i.e. (Q).
   - Deterministic bases: OPEN. For example, (2^N − 1 | N) = (M_q | p)(M_p | q) couples Legendre symbols of
     Mersenne numbers across the two primes.
   - Pure SLP-compressed c: PT-J's criterion β(c)β(−c) ≠ 0 stays OPEN. Counting shows most sign patterns need SLPs
     of size Ω(π(I)/log); that says nothing about structured ones.
7. **Experiments (pre-registered, b = 12, 14, 16; 1.1·10⁶ pairs at b = 16 via Theorem U).**
   - All proved statements pass: diagonal lemma 0 mismatches; PT-X, Theorem U and KS 0 violations.
   - The registered "fresh-random model" for nested zero sets **failed as registered**. The failures trace to exact
     integer identities of the tower at small exponents (FROZEN, e.g. 2 + 3 = 5) and to "closure" identities (e.g.
     a^r = −(a+1) ⇒ (a^r + 1)^r = a + 1 for odd r). Both reduce the node to a lower-depth (M/A) event.
   - After removing them, residual zero counts match the model within 0.88–1.04 on held-out b = 16, as predicted in
     a second-stage registration. Non-split gains of the nested batteries are ≤ 0.0015.
   - All 11 registered Jacobi statistics per size are null (|z| < 4); the positive control is found (z = 52 at b = 12).

So the strongest statement we can make: within PRIM, every non-split leak we can analyse is (M), (A) or (Q). The
unclassified remainder is (i) deterministic towers with ≥ 2 additive shifts, (ii) deterministic two-layer nodes at
the o(1)-density exceptional primes, (iii) exponential-degree pure post-processing of gate outputs (pure root
densities, slpcirc's τ-type territory), (iv) extension components, and (v) Jacobi of deterministic gate values and
of compressed pure integers. Nothing we computed points to a fourth family.

Subjective: probability that PRIM contains a sub-L[1/3] factoring route ≤ 1%.

---

## 1. The model PRIM

A PRIM program is a finite adaptive decision tree of N-independent shape. Each node computes an element of A/NA,
for a fixed free Z-algebra A (rank 1: Z; rank 2: Z[x]/(f), disc f prime to the primes of I). Allowed operations:
- ring operations with N-independent integer constants;
- pure powers x^M with fixed M (any size given by an SLP);
- exponent gates x ↦ x^{E(N)} with E from a fixed menu: affine αN + β (α ≥ 1, E ≥ 1 on the support) or floor
  ⌊(N − c)/2^j⌋. Gate outputs may be retained, nested (fed to later gates) and summed;
- coins: elements drawn uniformly from A/NA (only in §3);
- zero tests with gcd (a split is reported when the coordinate gcd with N is proper), and Jacobi symbols (c|N) of
  computed c ∈ Z/N.

The shape restriction is essential: an N-dependent SLP is vacuous (slpcirc Prop R), and an arbitrary exponent
E(N) = 1 + h(N) on the base −1 reads any predicate (leakmap §3.1). As codex notes, standard generic-ring models
already allow forward powering by public N (by square-and-multiply), but their theorems concern RSA challenges, not
challenge-free factor-bit leaks (§8).

**Depth.** The nesting depth of a node is the largest number of gates on a path through its expression in which
each gate's base contains the previous gate. A *shift* is an addition of a nonzero element that is not a monomial in
the same gate outputs. The *shift depth* is the number of shifts between gates on such a path. For example,
(2^N + 1)^N has depth 2 and shift depth 1.

## 2. Reduction to one prime: Theorem U and the permutation lemma

**Theorem U (Δ-form and uniform-r domination). PROVED.** Let p ≠ q ∈ I and n = p − 1. Suppose every component of
A ⊗ F_p is F_p (e.g. A = Z).
- (a) *Affine gates.* For x ∈ F_p, x^{αN+β} = x^{α r + β} with r := N mod n = q mod n, where 0^{(·)} := 0. Hence
  every node evaluates at p to Φ_v(p; r), a function of p and r only (by induction over the expression). For p < q,
  r = 1 + (q − p) = 1 + Δ. At q, N mod (q − 1) = p ≡ 1 − Δ.
- (b) *Floor gates.* ⌊(N − c)/2^j⌋ mod n depends only on N mod 2^j n, which is determined by q mod 2^j n.
- (c) *Counting.* An interval of length X contains at most 2 integers of a class mod n (n ≥ X), and at most 1 of a
  class mod 2^j n (j ≥ 1). So for any Z ⊆ Z/n, Pr_{q∼μ}[q mod n ∈ Z] ≤ 2|Z|·max_q μ(q). For μ uniform on the
  primes of I this is ≤ 2|Z|/π(I).
- (d) *Domination.* Let κ_v(p) be any per-node reference answer and Z^dev_v(p) = {r : answer of v at (p; r) ≠ κ_v(p)}.
  For a tree T and its reference transcript κ(p):

  Pr_q[τ_p ≠ κ(p) | p] ≤ (2/π(I)) Σ_{v ∈ path κ(p)} |Z^dev_v(p)|,

  and with PT-X, for iid uniform primes:

  Pr[correct or split] ≤ B + 2s + (8/π(I))·E_p[Σ_{v∈path κ(p)} |Z^dev_v(p)|].

*Proof.* (a) Fermat: for a unit x, x^E depends on E mod n, and N = pq ≡ q (mod n). For p < q ≤ 2X < 2p − 1 we have
q − n ∈ (1, n), so q mod n = q − p + 1. (b) N mod 2^j n = (p·(q mod 2^j n)) mod 2^j n, and
⌊(N − c)/2^j⌋ = k·n + ⌊(s − c)/2^j⌋ when N = k·2^j n + s. (c) as stated. (d) The first index where τ_p and κ(p) differ is
a node on κ(p)'s path whose answer at p deviates (leakmap NE1 proof). Union bound, then (c), then PT-X with
Pr[D] ≤ Pr[dev at p] + Pr[dev at q]. ∎

**Meaning.** ρ_v(p) := |Z^dev_v(p)|/n is the deviation probability when the hidden exponent is uniformly random.
It is a property of p alone, computed by running the program over F_p with exponent αr + β. Non-split leakage is
at most ≈ 16 ln X · E_p[Σ ρ_v(p)]. A leak needs some node with ρ_v(p) ≳ 1/(ln X·poly) on a non-negligible set of p.
The bound is lossy by the prime density (ln X), which is harmless for negligibility. The *signed* averages needed for
Jacobi symbols are **not** controlled by (c) (codex 3.3): domination bounds rare events, not cancellation.

*Extension components.* In a component F_{p^k}, k ≥ 2, x^N = Frob(x)^q depends on q mod (p^k − 1) ≫ |I|. Then (c)
loses a factor p^{k−1}, unless every gate base lies in a subgroup of exponent O(X) (for instance the norm-one torus,
order p + 1, for Lucas). **OPEN** in general; this is leakmap §3.5's short-interval problem.

**Lemma P (N-gates are permutations; codex prim_1 §4.2). PROVED.** For distinct p, q ∈ I, gcd(N, p − 1) =
gcd(N, q − 1) = 1. *Proof:* gcd(q, p − 1) = 1 since q > p − 1 ≥ q/2 (and q ≠ p − 1); if p | q − 1 then q ≥ 2p + 1 > 2X.
∎ So x ↦ x^N permutes F_p^* and F_q^*. Every prime partner gives a unit r for E = N. A node (u)^N = w with a unit
exponent has the unique solution u = w^{1/r} (e.g. (x^N + c)^N = 1 at p ⇔ x^N + c = 1; codex (4)).

## 3. Fresh coins: Theorem R, the inverse-fiber bound, and all-depth chains

**Theorem R (modal decoupling). PROVED.** Let T be a finite adaptive tree with zero/gcd tests whose nodes may use
coins drawn independently and uniformly from Z/N (or A/NA). Write a_ℓ(p|q) := Pr_coins[τ_p = ℓ | p, q]. By CRT the
coin components mod p and mod q are independent, so τ_p and τ_q are conditionally independent given (p, q).
Let m_p := argmax_ℓ a_ℓ(p|q) (ties broken by a fixed order), and let κ₀ be any map from primes to transcripts. Put
D_mod := {m_p ≠ κ₀(p)} ∪ {m_q ≠ κ₀(q)}. For independent p ∼ μ₁, q ∼ μ₂ (distinct, or uniform distinct ordered pairs
of a finite set):

  Pr[T outputs F(p, q) or splits] ≤ B + 4s + 2·Pr[D_mod ∧ m_p = m_q].

*Proof.* Fix (p, q). Write α_x = max_ℓ a_ℓ(x|·) and s(p, q) = 1 − Σ_ℓ a_ℓ(p)a_ℓ(q). The latter is the split
probability, by the diagonal lemma applied for each fixed coin value.
1. Σ_ℓ a_ℓ(p)a_ℓ(q) ≤ α_p, so Σ_{ℓ≠m_p} a_ℓ(p)a_ℓ(q) ≤ 1 − α_p ≤ s(p, q).
2. If m_p ≠ m_q then s(p, q) ≥ ½. If α_q ≤ ½, Σ ≤ α_q ≤ ½. Symmetrically if α_p ≤ ½. Otherwise a_{m_p}(q) ≤ 1 − α_q
   and Σ ≤ α_p(1 − α_q) + (1 − α_p)α_q = ½ − 2(α_p − ½)(α_q − ½) ≤ ½.
3. Pr[correct ∧ no split | p, q] = Σ_ℓ a_ℓ(p)a_ℓ(q)·1[F = o(ℓ)]. If m_p = m_q this is ≤ 1[F = o(m_p)] + s(p, q)
   by step 1. If m_p ≠ m_q it is ≤ 1 − s(p, q) ≤ ½ ≤ s(p, q). So Pr[correct ∧ no split] ≤ Pr[m_p = m_q,
   F = o(m_p)] + s.
4. On ¬D_mod, m = κ₀. As in PT-X, Pr[κ₀(p) = κ₀(q), F = o(κ₀)] ≤ B + (1 − α_{ℓ*}) ≤ B + Pr[κ₀(p) ≠ κ₀(q)], with
   α the law of κ₀(p). Also Pr[κ₀(p) ≠ κ₀(q)] ≤ Pr[m_p ≠ m_q] + Pr[D_mod ∧ m_p = m_q] ≤ 2s + Pr[D_mod ∧ m_p = m_q]
   by step 2.
5. Add: s + (s + Pr[D_mod ∧ m_p = m_q] + B + 2s + Pr[D_mod ∧ m_p = m_q]). ∎

(Exact enumeration check: `tests/test_prim.py::test_theorem_R_exact`, three fresh-base chain probes on all ordered
pairs of primes in (100, 200), three targets; the inequality holds with slack in every case.)

*Meaning.* For deterministic programs (no coins) R is PT-X. With coins, every bit of non-split information beyond the
modal transcript is paid for by splitting, with constant 4. Rare both-pass events such as Monier liars are absorbed
into s, since κ_pκ_q ≤ s when κ ≤ ½. What remains is the *modal* deviation: a node whose typical answer at p, given q,
differs from a p-only reference. No prior formulation located by us. It is the randomized analogue of PT-X.

**Lemma CB (inverse-fiber bound; codex prim_1 (23)). PROVED.** For X uniform on F_p^* and a chain
((X^{e₁} + c₁)^{e₂} + c₂ ⋯)^{e_k}, with positive exponents and arbitrary shifts, Pr[chain = w] ≤ ∏_i gcd(e_i, p − 1)/(p − 1).
Interspersed nonconstant polynomial maps of degree D_i multiply the numerator by ∏ D_i. (Shifts are injective; a
power map has fibres of size ≤ gcd(e, p − 1).) Tested: `test_collision_bound`.

**Corollary R-chain. PROVED.** Consider a tree of height h whose nodes are pure nodes and *fresh-base chain probes*
(A·chain − B with a fresh coin at the bottom, affine menu exponents E_i, ≤ k gates, polynomial degree ≤ D). Take κ₀(p)
to be the path on which every probe answers its formal value ("nonzero") and pure nodes answer as at p. Then:
- m_p ≠ κ₀(p) forces Pr_coins[τ_p ≠ κ₀(p)] ≥ ½, hence some probe on κ₀(p)'s path has local zero probability ≥ 1/(2h);
- by CB, that probe has G_v(p, q) := ∏_i gcd(E_i(N), p − 1) ≥ (p − 1)/(2hD);
- Pr_q[G_v ≥ T] ≤ k·Σ_{d|n, d ≥ T^{1/k}} Pr_q[d | α_i q + β_i] ≤ kα τ(n)(|I|T^{−1/k} + 1)/π(I).

Hence

  Pr[correct or split] ≤ B + 4s + O(h·kα·τ(p−1)·(2hD)^{1/k}·X^{−1/k}·ln X)

(iid uniform; the error is summed over the reference path). For fixed depth k this is X^{−1/k+o(1)} per probe.
It becomes vacuous at k ≈ log X / log log X. For E = N gates every gcd is 1 (Lemma P), so pure-N chains have exact
uniform laws at every depth (codex). The loss with depth is a feature of CB, not known to be real.

**Theorem W2 (precise two-layer rates; codex prim_1 (7), (8), (10)). PROVED.**
- With Y = X^{e₁} uniform on H₁ (|H₁| = n/d₁), |Pr[(Y + c)^{e₂} = w] − d₂/n| ≤ √p/|H₁| for attainable w ≠ 0.
- The retained form U^A(U + c)^B has the analogous rate g₀/n with g₀ = gcd(n, B, Ak).
- The core is the subgroup-intersection estimate ||H ∩ (aK − c)| − |H||K|/n| ≤ √p (Jacobi sums).

These give the exact (M) both-pass rates κ_pκ_q with explicit Weil error. Beyond two layers a direct Weil estimate
becomes trivial because the composed degree grows (codex 4.3). The classification does not need it: R + CB suffice.

**What R does not cover.** Nodes that are polynomials in several chain outputs of one fresh base (sums of nested
chains) need "local zero probability ≥ 1/(2h)" to be classified. That is OPEN beyond one-layer sparse sums (§4) and
bounded-degree Weil. The ε fiber problem of leakmap §3.3 (side information N mod 4) is not addressed by R.

## 4. Summed one-layer gates (deterministic bases): Theorem PS

At p with A = Z, a node Σ_m γ_m ∏_j y_j^{E_j(N)} (pure y_j, affine menu, t monomials after merging) is f(g^r) for a
t-sparse Laurent polynomial f with exponents e_m = Σ α_j·ind_g(y_j) (p-dependent). Here g is a generator, and
r ↦ g^r is a bijection Z/n → F_p^*.

**Lemma KS (sparse structure; codex prim_1 (12), degree compression after Kelley 2016 / CFKLLS 2000). PROVED.**
Let G ≤ F_Q^* be cyclic of order n and f a t-sparse Laurent polynomial (t ≥ 2) with nonzero coefficients. For every
1 ≤ L ≤ n, the roots not lying in a zero coset of a subgroup of order ≥ L number at most

  |R_L| ≤ 2n/(⌊n/L⌋ + 1)^{1/(t−1)} < 2n^{1−1/(t−1)}L^{1/(t−1)}.

The removed zero cosets can be taken to be ≤ ⌊n/L⌋ disjoint cosets of a single subgroup of order ≥ L.
*Proof sketch* (codex's): Dirichlet gives 1 ≤ e < ⌊n/L⌋ + 1 with all |e·a_i mod n| ≤ M := n/(⌊n/L⌋+1)^{1/(t−1)}.
The map x ↦ x^e sends G onto H of order n/gcd(e, n) ≥ L. On each coset bH, f(bx^e) is a Laurent polynomial with
exponents in [−M, M]: either zero (a zero coset) or with ≤ 2M roots in x, i.e. ≤ 2M/gcd(e, n) on bH. ∎
- *Partner congruences* (Kelley Thm 2.2, verified from the arXiv source): if f vanishes on bH_C then grouped
  coefficients cancel mod C, so every exponent has a partner ≡ mod C. The ratio of the partnered gate monomials has
  order dividing n/C: an order alignment, (M).
- (v0 had a weaker recursive form, 2⌊log₂t⌋n^{1−1/(e(1+⌊log₂ t⌋)(t−1))}, still valid; superseded.)
- Tested by brute force: `test_lemma_KS_bruteforce`, `test_lemma_KS_compression_bruteforce`.

**Theorem PS. PROVED.** For such a node, reference "formally nonzero", and L = n^θ (0 < θ < 1):

  Pr_q[v deviates at p] ≤ (2 + o(1))θ^{−1}·(k/φ(k))·dens(V_L(p)) + 4n^{1−(1−θ)/(t−1)}/π(I),

where V_L(p) is the union of the removed zero cosets, i.e. q-classes r ≡ r₀ (mod k), k | p − 1, k ≤ n^{1−θ}.
*Proof.* The residual uses Theorem U(c). For V_L use Brun–Titchmarsh (Montgomery–Vaughan) for each class mod
k ≤ X^{1−θ}: Pr_q[q ≡ r₀ (k)] ≤ 2|I|/(φ(k) log(|I|/k)π(I)). There are k·dens(V_L) classes. ∎
- The first term is the (M) part. It is the rate a uniformly random exponent would see (up to k/φ(k) ≪ log log X),
  on cosets where partnered monomials are order-aligned. It contains leakmap's SMALL and SUBGROUP events, e.g.
  2^{N−1} + 3^{N−1} − 2 vanishes on r ≡ 1 mod lcm(ord 2, ord 3).
- The residual is negligible (exp(−(log X)^ε)) for t ≤ (log X)^{1−2ε}, and X^{−Ω(1)} for t = O(1). It is vacuous for
  t ≳ log X / log log X (codex (16)).
- This is a per-node statement; a tree pays a union over its path.

**Exponential-degree pure post-processing. OPEN.** A gate followed by a pure circuit P of exponential degree (e.g.
(g^N + 1)^M − 1, M = lcm(1..B)) has, by Theorem U, deviation ≤ 2·#{roots of P on g's coset}/π(I): a pure root
density at p.
- Codex: (g^r + 1)^M − 1 has exactly gcd(M, n) − 1 zeros, which form a translated subgroup kernel, not cosets in r.
  At p = 241 it has 119 roots and its largest zero coset has order 5. This is an (M) event (smoothness of p − 1), but
  "M = cosets in r" is false.
- Classifying root densities of pure univariate circuits is slpcirc's open τ-type problem. The sparse case is
  Theorem SP there.

## 5. Deterministic two-layer gates

**Codex's bounds (prim_1 (3), (5)). PROVED.** For primitive g and c ≠ 0, #{r mod n : (g^r + c)^r = 1} ≤ τ(n)(1 + √p).
For a unit outer exponent and a target w of order h, the unit solutions number ≤ φ(h)·n·gcd(a₁, m)/m. Both are
order-alignment statements for small-order targets.

**Theorem D2 (one shift, general target). PROVED.**
- *Node at p:* F(r) = A(r)·u(r)^{αr+β} − B(r), with u(r) = y^{ar+b} + c (y pure, c ≠ 0, a, α ≥ 1).
- *Monomial condition:* A, B are nonzero monomials in gate outputs of pure bases, so B/A = γ·g^{λr} with γ ∈ F_p^*,
  λ ∈ Z/n, g a generator. The inner value may also be u = M₁(M₂ + c) with monomials M_i; the proof is unchanged.
- *Bound:* with o := ord_p(y^a), for every divisor d of n = p − 1,

  |Z(p)| := #{r ∈ Z/n : F(r) = 0} ≤ gcd(α, d)·τ(d)·(n/d + d·(n/o)·(√p + 2)) + n/o.

*Proof.*
1. If u(r) = 0, F(r) = −B(r) ≠ 0 unless B = 0. The set {u = 0} has ≤ n/o elements (u depends on r mod o).
2. Otherwise F(r) = 0 ⟺ (αr + β)·L(r) ≡ Γ + λr (mod n), with L(r) = ind_g u(r) and Γ = ind_g γ.
   Reducing mod d, the pair (r mod d, L(r) mod d) lies in S_d = {(x, ℓ) : (αx + β)ℓ ≡ Γ + λx (mod d)}.
3. For each x there are ≤ gcd(αx + β, d) values ℓ. Hence |S_d| ≤ gcd(α, d)·Σ_{z mod d} gcd(z, d) ≤ gcd(α, d)·d·τ(d).
4. *Cell count.* For (x, ℓ) ∈ (Z/d)², let N(x, ℓ) := #{r : u(r) ≠ 0, r ≡ x, L(r) ≡ ℓ (mod d)}. Expand with additive
   characters e(s(r − x)/d) and multiplicative characters χ of order dividing d (χ(0) = 0).
   - The principal pair gives (n − #{u = 0})/d².
   - Pairs (s ≠ 0, χ = 1) give Σ_r e(sr/d)·1[u ≠ 0] = −Σ_{u(r)=0} e(sr/d), of size ≤ n/o.
   - For χ ≠ 1: T(s, χ) = Σ_{r∈Z/n} e(sr/d)χ(y^{ar+b} + c). Write r = z + oj (o | n). The j-sum vanishes unless
     d | so; then e(sz/d) = ψ(y^{az}) for a character ψ of K = ⟨y^a⟩. So T = (n/o)Σ_{h∈K} ψ(h)χ(y^b h + c).
     With k = n/o and an extension Ψ of ψ, this is (n/(ok))Σ_{w∈F_p^*} Ψ(w^k)χ(y^b w^k + c). The polynomial
     w^{·}(y^b w^k + c)^{·} has ≤ k + 1 distinct roots, the roots of y^b w^k + c being simple and χ ≠ 1, so it is not
     an n-th power. Weil (Lidl–Niederreiter Thm 5.41) gives |Σ_{w∈F_p}| ≤ k√p. Removing w = 0 costs 1. So
     |T| ≤ (n/o)(√p + 1).
   - Summing, |N(x, ℓ) − n/d²| ≤ (n/o)(√p + 2).
5. |Z| ≤ n/o + |S_d|·max N(x, ℓ). ∎

(The numerical companion `d2_bound` checks it at all 12-bit primes for D2a_w and D2d: `test_theorem_D2_bound_holds`.
At these sizes the bound is above n, so the check has no content beyond consistency.)

**Corollaries.**
- (i) *Per prime.* If o ≥ p^{1−η} and n has a divisor d ∈ [p^δ, p^{1/4−η/2}], then |Z(p)| ≤ p^{1−δ+o(1)} + p^{3/4+η/2+o(1)}.
  The deviation probability is then ≤ 2|Z|/π(I) = p^{−min(δ, 1/4−η/2)+o(1)}.
- (ii) *Exceptional primes.* The primes failing (i) have p − 1 without a prime factor in (p^δ, p^{1/4−η/2}], or have
  ord_p(y^a) < p^{1−η}. Both are multiplicative, (M)-type properties of p.
  - The first set has relative density ≤ C·δ/(1/4 − η/2) by an upper-bound sieve for shifted primes (Brun/Selberg;
    standard; exact reference UNVERIFIED).
  - Under GRH the second has density o(1) (Erdős–Murty-type results on ord_p(a); reference and form UNVERIFIED).
  - So under GRH, for any fixed δ, all but a proportion ≤ Cδ + o(1) of p satisfy (i). Unconditionally (i) holds at
    each prime that meets its hypotheses.
- (iii) *Not a negligible bound.* The per-node deviation rate in the iid population is ≤ Cδ + o(1) for every δ, i.e.
  o(1). It is not n^{−ω(1)}, so D2 does not give a poly-height tree theorem. It does say that a *large* zero set
  (density ≥ p^{−δ}) at a deterministic two-layer node can occur only at primes whose p − 1 lacks medium divisors,
  or where ord_p(y) is small.

**Reformulation (codex (6)). PROVED.** For primitive g and w ≠ 0, (g^r + c)^r = w ⟺ ind(u)·ind(u + c) ≡ ind(w)
(mod p − 1) with u = g^r: the "Zech hyperbola". D2 bounds it through its reductions mod divisors d of p − 1.
- No theorem bounds it uniformly. Closest known results concern g^x = x (Glebsky–Shparlinski 2010; Bourgain–Konyagin–
  Shparlinski 2008) and x^x = c (Balog–Broughan–Shparlinski 2011; Cilleruelo–Garaev 2016). None of them covers the
  product-of-indices equation (codex, verified table in prim_1 §1.3).
- Exceptional-prime-free version: **OPEN**. Equivalently, power-saving bounds for
  Σ_{s∈(Z/n)^*} e_p(λg^{1/s} + μw^s), the exponential map composed with modular inversion.

**Lemma FZ (frozen exponents). PROVED.** For an integer t at which all gate exponents are positive, let m_p(t) be the
lcm of the orders at p of every gate-base value in the evaluation at t. For r ≡ t (mod m_p(t)), every gate output
equals its value at t, by induction along the expression, so v(r) = T_v(t) mod p. Here T_v(t) ∈ Z[1/·] is the integer
tower. So frozen zeros are:
- identities (T_v(t) = 0, at every p), e.g. 2 + 3 − 5 = 0 for S1 at t = 1, (3/2)(4/3) − 2 = 0 for S5 at t = −1,
  (2 + 1)^1 − 3 = 0 for D2a_3 at t = 1;
- or low-height coincidences p | T_v(t).

In the q-world, r_q ≡ t (mod m_p(t)) is an (M) alignment. When m_p(t) = p − 1 it is the pin q = p + t − 1, an (A)
event. leakmap's ARCH (2^N = 2^k) and the twin-prime Lucas liars are instances.

**Lemma CL (closure). PROVED.** If at r an intermediate gate input equals ±∏ y_i^{k_i} (pure constants y_i), the gate
output is ±∏ (y_i^{E})^{k_i}, a monomial in lower-depth gate outputs. On the set where this happens, the node agrees
with a node of smaller depth, and the set itself is the zero set of a lower-depth node. Example: a^r = −(a + 1) and
r odd give (a^r + 1)^r = (−a)^r = a + 1. At depth 2 both pieces are one-layer nodes, hence (M) by NE1/PS.
- *Control (scratch, b = 14, random bases a, unit r, frozen classes removed):* solutions of (a^r + 1)^r = w number
  2.03× the fresh-random model for w = a + 1, and 0.87–1.13× for w ∈ {a + 2, random, (a + 1)², (a + 1)^{−1}}.
- This is exactly the excess seen in E1 (§9).

**Depth ≥ 3 (shift depth ≥ 2), deterministic. OPEN.** D2's mechanism needs the joint law of (r mod d, ind(u) mod d)
for an inner value u that is a binomial in gate outputs. At shift depth 2 the inner value is (y^{ar} + c₁)^{e(r)} + c₂.
Its multiplicative characters are not characters of an algebraic function of y^{ar} of bounded degree, so no Weil or
Jacobi-sum estimate applies. This is the first place where our methods break for deterministic bases. Experiments on
D3 (shift depth 2) show the same FROZEN/CLOSURE + fresh-random behaviour as depth 2 (§9).

## 6. Jacobi symbols of computed values

(c|N) = χ_p(c_p)χ_q(c_q) never splits unless c ≡ 0 at a prime, and then gcd splits. PT-X and R do not apply.

**Theorem J1 (fresh-base gate polynomials). PROVED** (core estimate codex prim_1 (21)–(22)).
- *Setting:* population N ≡ 1 (mod 4), e = ±1 equally likely, p, q iid μ_e. There are k Jacobi tests
  J_i = (P_i(x^{E(N)}) | N), with x a fresh uniform coin (shared), P_i ∈ Z[T] of degree ≤ D, and E an affine menu
  exponent. For S ⊆ [k] put P_S = ∏_{i∈S} P_i and g_x = gcd(E(N), x − 1). Then:

  Pr[correct or split] ≤ ½ + (2^{k−1} + 1)s + Σ_{S≠∅, degenerate} |β(a_S)β(−a_S)|
                         + ½Σ_{S≠∅} E[min(1, (kD√p+1)g_p/(p−1))·min(1, (kD√q+1)g_q/(q−1)) + 2kD(g_p/(p−1) + g_q/(q−1))].

- *Degenerate:* S is degenerate iff P_S = a_S T^j Q(T)², a_S a pure integer.
- *Proof outline:*
  - Fourier over {±1}^k: Pr[correct ∧ all nonzero] ≤ ½(1 − s) + ¼Σ_S |E_+[J_S] − E_−[J_S]| + 2^{k−2}(s_+ + s_−).
  - E_e[J_S | p, q] = m_{S,p}m_{S,q} with m_{S,x} = E_{h∈H_x} χ_x(P_S(h)), H_x = image of the gate at x, |H_x| = n_x/g_x.
  - Non-degenerate S: |m| ≤ (kD√p + 1)/|H| (Weil over the subgroup; codex (21)).
  - Degenerate S: m_{S,x} = χ_x(a_S)·ν + O(kD/|H_x|), with ν = 1 iff j is even or g_x is even (H_x ⊆ squares).
    Since g_x ≡ E(N) (mod 2), ν is public, and the term is PT-J for the pure integer a_S.
- *Classification:* a fresh-base gate-Jacobi leak of ε is either an (M) alignment (g_p g_q ≳ √(pq)/poly; for affine
  menus this has probability X^{−1/2+o(1)}) or (Q) via PT-J for a low-height pure a_S. The latter is negligible by
  PT-J1 (log|a_S| = O(kD·height)).
- Codex's correction is built in: P(T) = T on a subgroup of squares has bias 1. That is the j-odd, g-even case,
  i.e. public parity.

**Deterministic bases. OPEN.**
- Theorem U does not control signed averages. At p, (2^N − 1 | p) = (M_q | p) with M_q = 2^q − 1, so
  (2^N − 1 | N) = (M_q | p)(M_p | q).
- Bounds for Σ_{q≤T prime} χ_p(M_q) at a single modulus exist (Banks–Conflitti–Friedlander–Shparlinski 2004; exact
  statement UNVERIFIED here). The two-prime coupled product has none that we know of.
- E1's six deterministic gate-Jacobi statistics are null (§9).

**Compressed pure c. OPEN.** PT-J: (c|N) leaks ε iff β(c)β(−c) ≠ 0.
- Codex's counting (prim_1 (20)): SLPs of length s give ≤ (s + 2)[3(s+2)²]^s sign patterns on π(I) primes, so most
  CRT patterns need s log s = Ω(π(I)). This does not exclude a structured short c.
- The natural families (Fermat numbers, a^M − 1, Lucas values, ∏(x² + 1) over poly(n) ranges) give no positive
  answer. For the last, PT-J1(b) applies under GRH since log c = poly(n) ≤ X^{1/4}.
- Biases created by x | c are splits (M), not Jacobi leaks.
- E3 (leakmap: 345 values) and J7/J8 here are null.
- Prime averages of a real character of doubly-exponential conductor and small SLP complexity: no result known to us.

## 7. Depth hierarchy (Q4)

| fragment | status | loss / bound | where the proof breaks |
|---|---|---|---|
| pure (depth 0) | PT (imported) | 2s | — |
| fresh coins, chains of any depth | **PROVED** (R + CB) | 4s + X^{−1/k+o(1)} per probe | k ≈ log X/log log X (CB is lossy) |
| fresh coins, polynomials in several chain outputs | R reduces to local zero-prob ≥ 1/(2h) | — | OPEN beyond sparse one-layer / bounded-degree Weil |
| fresh coins, precise two-layer rates | **PROVED** (W2, codex) | √p/\|H₁\| | depth 3: composed degree ~ p, Weil trivial |
| deterministic, 1 gate, low degree | NE1 (imported) | D·μ(o) | — |
| deterministic, 1 layer, t-sparse sums | **PROVED** (PS) | (M)-part + 4n^{1−(1−θ)/(t−1)}/π(I) | t ≈ log X/log log X |
| deterministic, 1 layer + exponential-degree pure circuit | reduces to pure root density (U) | — | OPEN (slpcirc τ-type) |
| deterministic, depth 2, one shift | **PROVED** at good primes (D2); w = 1 √p bound (codex) | p^{−δ} at good p; o(1)-density exceptions (GRH) | exceptional primes (rough p − 1, small ord) |
| deterministic, shift depth ≥ 2 | OPEN | — | no character-sum structure ("doubly Zech") |
| Jacobi, fresh base | **PROVED** (J1) | (M) + PT-J(a_S) | — |
| Jacobi, deterministic gates | OPEN | — | signed averages over prime partners |
| Jacobi, pure compressed c | PT-J criterion; OPEN | — | character sums of huge conductor |
| extension components F_{p^k} | OPEN | — | short intervals (already at depth 1) |

So for **deterministic** bases the first break is at shift depth 2 (depth 3). Depth 2 is proved only off an
exceptional set of primes, with an o(1) loss. For **fresh** bases nothing breaks at a fixed depth for chains.
"Depth" is not the right parameter for post-processing: exponential-degree pure circuits already escape at gate
depth 1 (codex 4.4).

## 8. Literature (Q5)

No result located classifies PRIM or any challenge-free leak model with N-exponent gates (UNVERIFIED as a negative).
Codex verified the following; details are in prim_1 §5.
- **Damgård–Koprowski** (EUROCRYPT 2002, LNCS 2332, 256–271): generic root extraction in hidden-order groups, with
  no ring addition.
- **Leander–Rupp** (ASIACRYPT 2006, LNCS 4284, 241–251): generic low-exponent RSA ≡ factoring. Their §3.1 explains
  why exponential circuit degree defeats root counting, which is the same obstruction as our exponential-degree
  post-processing.
- **Aggarwal–Maurer** (EUROCRYPT 2009, LNCS 5479, 36–53; IEEE T-IT 62(11) 2016): generic RSA ≡ factoring.
- **Jager–Schwenk** (J. Cryptology 26 (2013) 225–245): generic Jacobi/QR hardness. This is the reverse direction:
  PRIM supplies Jacobi as a primitive.
- **Rotem–Segev** (CRYPTO 2020, LNCS 12172, 481–509): generic-ring delay functions; speeding up repeated squaring
  is equivalent to factoring. This is the closest in spirit, as a generic-ring analysis of N-dependent powering, but
  it is about computing x^{2^T} faster, with a challenge.
- **Dachman-Soled–Loss–O'Neill** (ITC 2024, LIPIcs 304, 8:1–8:24): generic RSA with preprocessing.
- **Maurer–Wolf** (EUROCRYPT 1998) and **Boneh–Venkatesan** (EUROCRYPT 1998).

All of these require inverting a random challenge, while a PRIM leak computes a factor bit from N alone. Their
"component-separating zero ⇒ gcd" engine is the diagonal lemma behind PT-X/R. Exponential congruences:
Glebsky–Shparlinski (DCC 56 (2010) Thm 2), Balog–Broughan–Shparlinski (Acta Arith. 148 (2011)), Cilleruelo–Garaev
(Proc. AMS 144 (2016) Thms 1–2), Kurlberg–Luca–Shparlinski (MRL 22 (2015)), Bourgain–Konyagin–Shparlinski (IMRN 2008).
Sparse roots: Kelley (LMS JCM 19A (2016) Thms 2.2–2.3; checked against the arXiv source), CFKLLS (Israel J. Math.
120 (2000) Lemma 7).

## 9. Experiments (pre-registered; `notes/prim_prereg.md`, commit f8952d0)

Populations: all pairs at b = 12 (7 750) and b = 14 (92 665); 200 000 sampled pairs at b = 16. Full uniform-r zero
sets were computed at every prime of each range (125, 431, 1 491 primes). The b = 16 post-hoc q-world uses all
1.1·10⁶ pairs via Theorem U. 32 affine nodes, 1 floor node, 1 Lucas node, 11 Jacobi statistics per size plus the J_CRT control at b ≤ 14. Peak RSS ≤ 0.93 GB.

**Registered predictions.**

| | b = 12 | b = 14 | b = 16 | verdict |
|---|---|---|---|---|
| P1 diagonal lemma | 0 mismatches | 0 | 0 (200 000 pairs × 36 nodes) | **pass** |
| P2 PT-X (all cells) | holds | holds | holds | **pass** |
| P3 Theorem U | 0 violations | 0 | 0 | **pass** |
| P4 KS (750/2 586/8 946 node-primes) | 0 | 0 | 0 | **pass** |
| P5 fresh-random model, full-r | (not registered; 4/8 would fail) | 6/8 fail | 6/8 fail | **refuted as registered** |
| P6 model, q-world ratio | — | D2a_1, D2a_3, D3 fail | D2a_1, D2a_3, D2b, D2d, D3 fail | **refuted as registered** (and underpowered) |
| P6 two-sided total ≤ 3·model + 5 | 1 vs 0.7 | 8 vs 2.2 | 7 vs 1.1 | **pass** |
| P7 summed residual | pass | S1, S5 fail | S1, S5 fail | **refuted as registered** for S1, S5 |
| P8 Jacobi, \|z\| < 4 | all | all (J7_{14,15,16}: dyadic variance ≤ 0, iid z ≤ 0.3) | all (max \|z\| 3.05, J7_19) | **pass** (with the noted undefined estimator) |
| P8 positive control J_CRT z > 10 | 51.8 | dyadic undefined; iid 151.6 | — | **pass at b = 12; estimator failure at b = 14** |
| ACCIDENTAL splits | 0 | 0 | 0 | **pass** |

Non-split Bayes gains (Pr[correct ∧ no split] − (1 − s)B, all pairs):
- NEST battery: |gain| ≤ 0.0022 (b = 12), ≤ 0.0015 (b = 14), ≤ 0.0006 (b = 16). s = 0.024, 0.008, 0.002.
- SUM battery: ≤ 0.0014. LUC: ≤ 0.0004.
- CTRL/ALL: 0.02–0.07 on ε and [p ≡ 1 mod 3]. This comes from the pure p − 1 controls (C_PM1, C_P), i.e. PT's
  known smooth-class correlation (factorbits §3.2), with s = 0.24–0.53.

**Why P5/P6/P7 failed (diagnosed on b ≤ 14, then tested on held-out b = 16).**
1. *FROZEN identities* (Lemma FZ). These are 38–57% of D2a_1/D2a_3/D2b/D3 zeros and 43–61% of S1/S2/S5 zeros.
   Examples: (2 + 1)^1 = 3; 2 + 3 = 5; (3/2)(4/3) = 2; and r ≡ 0 (mod ord 2) ⇒ (2^r + 1)^r = 1. They form cosets of
   all sizes, including the sub-L cosets that broke P7.
2. *Self-consistent subgroup cosets.* For example, 2^r = 1 makes the D2a_1 tower constant on a coset of size
   n/ord_p(2). These are the "unflagged primes with V_L" in P5.
3. *CLOSURE identities* (Lemma CL), e.g. a^r = −(a + 1). These are the 2× excess of D2a_3 after FROZEN.
4. *Permutation lemma.* For E = N and a prime partner (r a unit), D2a_1 has no zeros at all: u^r = 1 forces
   u = 1, impossible for u = 2^r + 1. The registered model counted u = 1, so P6's 0-vs-13 is exact, not noise. D3
   vanishes for unit r only if 2^q ≡ −1 (mod p), when the inner value is 0; this is a FROZEN-type (M) event, and none
   occurred.
5. *P6 power.* The expected one-sided counts were 7–25 per node, so the registered [0.6, 1.6] window was too narrow
   for the D2d failure (4 vs 7).

**Second-stage predictions (`notes/prim_posthoc_prereg.md`, fixed before the b = 16 post-hoc run).** Residual =
zeros outside COSET ∪ FROZEN ∪ F1 ∪ CLOSURE, compared with the fresh-random model with shift exclusion on the same
residues.

| node | b = 14 (development) full-r | b = 16 full-r (Q1: [0.8, 1.25]) | b = 16 q-world (Q2: [0.6, 1.6]) |
|---|---|---|---|
| D2a_1 | 0.884 | 0.955 | 0 vs 0 (exact, Lemma P) |
| D2a_2 | 0.866 | 1.036 | 1.254 (48 vs 38) |
| D2a_3 | 0.932 | 0.980 | 0.810 |
| D2b | 0.884 | 0.955 | 0.982 (541 vs 551) |
| D2c | 1.052 | 0.986 | 0.968 |
| D2d | 0.960 | 0.974 | 0.680 |
| D2e | 0.877 | 0.882 | 0.863 |
| D3 | 0.980 | 1.006 | 0 vs 0 (only 2^q ≡ −1 events possible) |

- Q1 and Q2 pass.
- **Q3:** two-sided OTHER/OTHER zeros over all 1.1·10⁶ pairs = 8 (all D2b) vs model 5.2. Pass (≤ 20.6).
- D2b's 41 two-sided zeros are mostly COSET/FROZEN: base-2 Fermat-type alignments, (M).
- Summed nodes: zero two-sided events at b = 16.

**Reading.** These are **experimentally supported** at b ≤ 16, i.e. primes ≤ 2^16, and are weak evidence
asymptotically. Deterministic nested zero sets were checked exhaustively at every prime. They are explained exactly
by FROZEN/CLOSURE identities and small-order cosets, which are (M/A) via Lemmas FZ, CL and PS. What remains
behaves like a fresh-random inner value, whose rates are (M) by W2. No unexplained structure was seen.

Deviations from the registration: (D1) the post-hoc classes FROZEN, F1 and CLOSURE and the shift-excluded model were
added after seeing b ≤ 14 data; only b = 16 is held out. (D2) For ε and [p ≡ 1 mod 3] on their registered
sub-populations, P2's PT-X bound is the fiber form and is trivially 1. The meaningful check is on all pairs, which was
also run (leakmap D4 had the same issue). (D3) The dyadic-robust variance was centred (as in leakmap) after a b = 12
debugging run that had used uncentred residuals. (D4) The b = 12 summaries were regenerated by later code versions;
the registered fields are unchanged.

## 10. Candidate cards

| Mechanism | Model | Outcome |
|---|---|---|
| Nested N-th powers (g^N + c)^N, deterministic | PRIM depth 2 | (M) at good primes (D2); experiments: FROZEN/CLOSURE + random |
| Nested chains with fresh bases, any depth | PRIM-R | **refuted for this model** beyond 4s + modal (R + CB) |
| Sparse sums of gates | PRIM depth 1 | **refuted for this model** beyond (M) cosets (PS) |
| Gate + smooth pure power (g^N + 1)^M | depth 1 + pure | (M) via smoothness; root-density classification OPEN |
| Jacobi of fresh-base gate polynomials | PRIM-R + Jacobi | **refuted for this model** beyond (M) + PT-J (J1) |
| Jacobi of deterministic gate values (Mersenne-type) | PRIM + Jacobi | OPEN; 6 statistics null |
| Compressed pure Jacobi (Fermat numbers, 2^M − 3) | PT-J | OPEN; null |
| Three-layer towers ((2^N+1)^N+1)^N | PRIM shift depth 2 | OPEN; experiments as depth 2 |

## 11. Reviews

(to be filled: codex adversarial round prim_2; Opus proof critic; Opus literature critic.)

## 12. Reproduction (project root)
```
.venv/bin/python -m pytest -q tests/test_prim.py                 # 14 tests, ~2 s
.venv/bin/python -m src.prim.e1 --b 12                           # ~5 s
.venv/bin/python -m src.prim.e1 --b 14                           # ~1 min
.venv/bin/python -m src.prim.e1 --b 16 --pairs 200000            # ~13 min, peak RSS 0.93 GB
.venv/bin/python -m src.prim.e1_posthoc --b 14                   # ~1 min
.venv/bin/python -m src.prim.e1_posthoc --b 16                   # ~15 min
```
