# prim: classifying non-split leaks in PRIM (phase 25)

Agent `prim`, 2026-10-23. **Version v2**: after the codex adversarial round (`codex/prim_2.out`, MAJOR REVISION), the
Opus proof critic (ACCEPT WITH FIXES) and the Opus literature critic (OK WITH FIXES); all findings applied (§11).
- Pre-registration: `notes/prim_prereg.md` (commit `f8952d0`, before any code). Second-stage predictions for the
  post-hoc classes: `notes/prim_posthoc_prereg.md` (commit `59906fb`, before the b = 16 post-hoc run; addendum commit
  `090997f`, before the b = 15 confirmation run on an untouched range).
- Code: `src/prim/{prog,zerosets,e1,e1_posthoc,theory}.py`; tests `tests/test_prim.py` (15 pass).
- Data: `experiments/raw/prim_e1A_b{12,14,16}.jsonl`; summaries `experiments/analysis/prim_e1_b{12,14,16}.json`,
  `prim_e1_posthoc_b{14,15,16}.json`.
- GPT-6.1-sol: blind round `codex/prim_1.{md,out}` (several results below are codex's and are credited as such);
  adversarial round `codex/prim_2.{md,out}`.
- Labels: **PROVED** (proof here or in the cited codex file), **KNOWN** (reference), **HEURISTIC**, **CONJECTURED**,
  **OPEN**, **experimentally supported**, **refuted for this model**, **UNVERIFIED**. "No prior formulation located
  by us" is the strongest novelty statement made.

Notation as in leakmap: N = pq, p ≠ q primes of I = (X, 2X] with X an integer (for theorems; experiments use b-bit
primes with the top two bits set), iid laws on distinct pairs unless stated, n = n_p = p − 1, s = Pr[split],
B = max_b Pr[F = b]. A node's *answer* at p is [node ≡ 0 mod p]; τ_p is p's diagonal transcript. PT-X (leakmap §2):
Pr[correct or split] ≤ B + 2s + 2Pr[D ∧ no split] for any one-prime reference κ.

---

## 0. Verdict (plain language)

**No leak was found. PRIM is not classified as a whole; a list of fragments now has proved quantitative bounds.**
(M), (A), (Q) remain descriptive labels; the theorems are bounds, not a taxonomy theorem.

Proved fragment bounds (zero/gcd trees unless stated):
1. **Reduction to one prime and a hidden exponent (Theorem U).** For affine gates at a prime p (all components F_p),
   every gate x^{E(N)} equals x^{αr+β} with r = q mod (p − 1) ≡ 1 + (q − p). Each residue r is hit by at most two
   partner primes. So a node's deviation probability at p is ≤ (2/π(I))·|{bad r}|, a count depending on p alone (the
   *uniform-r zero set*). Floor gates need q mod 2^j(p − 1). The reduction bounds rare events only; it gives nothing
   for signed (Jacobi) averages, and nothing in extension components F_{p^k}, k ≥ 2 (short intervals).
2. **Fresh coins (Theorem R, a challenge-free reformulation of the Aggarwal–Maurer typical-path / Leander–Rupp
   simulation argument).** For any program with coins uniform in Z/N: Pr[correct or split] ≤ B + 4s plus "modal
   deviations" (the typical answer at p, given q, differs from a p-only reference).
   - With the inverse-fiber bound (codex), fresh-base **chains** of fixed depth k are classified: a modal deviation
     needs ∏ gcd(E_i(N), p − 1) ≳ (p − 1)/(2hD), an (M) alignment of probability X^{−1/k+o(1)}.
   - The bound is vacuous for k ≳ log X/log log X. Chains of pure N-th-power gates are permutations of F_p at every
     depth (codex).
3. **Sums of one-layer gates (Theorem PS, on a degree-compression lemma after Bi–Cheng–Rojas and Kelley).** Zeros lie
   in zero cosets of one large subgroup, where partnered monomials are order-aligned (M), or in a residual of size
   ≤ 2n^{1−(1−θ)/(t−1)}. The bound is useful for t up to ≈ log X/log log X terms.
4. **Deterministic two-layer nodes with a binomial inner value (Theorem D2).** For A·(y^{ar+b} + c)^{αr+β} = B with
   A, B monomials, for every divisor d of p − 1:
   |zeros| ≤ gcd(α,d)·τ(d)·(n/d + d(n/o)(√p + 2)), o = ord_p(y^a).
   - This is a power saving at primes where p − 1 has a divisor in [p^δ, p^{1/4}] and y has order ≥ p^{1−η}.
   - For a fixed non-torsion base under GRH, the exceptional primes have density ≤ Cδ + o(1). This does **not**
     establish negligibility.
   - For SLP-sized bases (e.g. 2^{lcm(1..256)}) most primes can be exceptional.
   - With a unit exponent and a small-order target (e.g. w = 1), codex's √p-type bound applies. On prime partners
     this also reduces the registered three-layer D3 exactly to a depth-2 condition.
5. **Jacobi of fresh-base gate polynomials (Theorem J1, repaired after a codex counterexample).** For bounded-height
   polynomials, the only biases are (M) gcd alignments, splits/double zeros, or a degenerate square shape that
   reduces to PT-J for a pure integer, i.e. (Q). With SLP-sized coefficients the polynomial can degenerate modulo
   every sampled prime, and codex's example then leaks ε with success 0.75. That is an exponential-height CRT
   construction, the analogue of PT-J's CRT c.

Where our methods stop (OPEN):
- deterministic towers with ≥ 2 additive shifts and a target that is not of small order;
- depth-2 nodes with ≥ 3-term inner values or sums of nested gates;
- D2's exceptional primes, and large-height bases;
- exponential-degree pure post-processing (pure root densities, slpcirc's τ-type problem);
- extension components;
- Jacobi of deterministic gate values, e.g. (2^N − 1 | N) = (M_q | p)(M_p | q) with Mersenne numbers;
- Jacobi of compressed pure integers (PT-J's β(c)β(−c) ≠ 0);
- the ε fiber problem (side information N mod 4) for all zero-test reductions.

Experiments (pre-registered; b = 12, 14, 16 and a confirmation range b = 15):
- **Proved statements pass:** diagonal lemma 0 mismatches; PT-X, Theorem U and the sparse bounds 0 violations.
- **The registered fresh-random model for nested zero sets failed as registered.** The failures trace to exact
  identities: tower identities at small exponents (FROZEN, e.g. 2 + 3 = 5) and closure identities (a^r = −(a + 1)
  ⇒ (a^r + 1)^r = a + 1 for odd r). These remove a gate and leave lower-depth (M/A) conditions.
- **Residual after removing them.** On the remaining zeros, counts match the model:
  - 0.88–1.04 on b = 16 (a prospectively registered reanalysis of partially exposed data);
  - 0.85–0.99 (full-r) and 0.68–1.16 (q-world) on the untouched b = 15 range, as predicted in an addendum
    committed before any b = 15 computation.
  - This is empirical agreement only, not a fitted stochastic test.
- **Gains and Jacobi.** Non-split Bayes gains of the nested batteries are ≤ 0.0022. All defined registered Jacobi
  statistics are null (|z| < 4); four b = 14 cells and the b = 14 positive control have an undefined dyadic variance
  and are "not established".
- **A further attack failed.** Codex tried genuine two-shift towers such as ((2^N+1)^N+2)^N − 3 on all b = 12 pairs:
  no two-sided zero, and a non-split gain of exactly s/2.

Subjective (not a mathematical conclusion): probability that PRIM contains a sub-L[1/3] factoring route ≤ 1%.

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
- (c) *Counting* (X an integer). The X integers of I contain at most one integer of each class mod n (n ≥ X), and
  at most 1 of each class mod 2^j n (j ≥ 1); we keep the conservative constant 2 below. So for any Z ⊆ Z/n, Pr_{q∼μ}[q mod n ∈ Z] ≤ 2|Z|·max_q μ(q). For μ uniform on the
  primes of I this is ≤ 2|Z|/π(I).
- (d) *Domination.* Let κ_v(p) be any per-node reference answer and Z^dev_v(p) = {r : answer of v at (p; r) ≠ κ_v(p)}.
  For a tree T and its reference transcript κ(p):

  Pr_q[τ_p ≠ κ(p) | p] ≤ (2/π(I)) Σ_{v ∈ path κ(p)} |Z^dev_v(p)|,

  and with PT-X, for iid uniform primes:

  Pr[correct or split] ≤ B + 2s + (8/π(I))·E_p[Σ_{v∈path κ(p)} |Z^dev_v(p)|].

  For floor gates the deviation sets live in Z/2^j n and the factor 2 becomes 1.

*Proof.* (a) Fermat: for a unit x, x^E depends on E mod n, and N = pq ≡ q (mod n). For p < q ≤ 2X < 2p − 1 (X an integer) we have
q − n ∈ (1, n), so q mod n = q − p + 1. (b) N mod 2^j n = (p·(q mod 2^j n)) mod 2^j n, and
⌊(N − c)/2^j⌋ = k·n + ⌊(s − c)/2^j⌋ when N = k·2^j n + s. (c) as stated. (d) The first index where τ_p and κ(p) differ is
a node on κ(p)'s path whose answer at p deviates (leakmap NE1 proof). Union bound, then (c), then PT-X with
Pr[D] ≤ Pr[dev at p] + Pr[dev at q]. ∎

**Scope.** U concerns zero/gcd trees with affine gates (floor gates: the modulus 2^j n; mixed menus: a common
enlarged modulus). The reference κ(p) is a function of p alone; contrast R's partner-dependent mode m_p(p, q). Under
uniform distinct pairs from M primes the partner mass is 1/(M − 1); the factor 8/π(I) stays conservative.

**Meaning.** ρ_v(p) := |Z^dev_v(p)|/n is the deviation probability when the hidden exponent is uniformly random.
It is a property of p alone, computed by running the program over F_p with exponent αr + β. Non-split leakage is
at most ≈ 16 ln X · E_p[Σ ρ_v(p)]. A leak needs some node with ρ_v(p) ≳ 1/(ln X·poly) on a non-negligible set of p.
The bound is lossy by the prime density (ln X), which is harmless for negligibility. The *signed* averages needed for
Jacobi symbols are **not** controlled by (c) (codex 3.3): domination bounds rare events, not cancellation.

*Extension components.* In a component F_{p^k}, k ≥ 2, x^N = Frob(x)^q depends on q mod (p^k − 1) ≫ |I|. Then (c)
loses a factor p^{k−1}, unless every gate base lies in a subgroup of exponent O(X) (for instance the norm-one torus,
order p + 1, for Lucas). **OPEN** in general; this is leakmap §3.5's short-interval problem.

**Lemma P (N-gates are permutations; codex prim_1 §4.2). PROVED.** For distinct p, q ∈ I, gcd(N, p − 1) =
gcd(N, q − 1) = 1. *Proof:* WLOG p < q. gcd(N, p − 1) = gcd(q, p − 1) = 1 since q > p − 1. If p | q − 1 then
q − 1 ≥ 2p (q − 1 even), so q > 2X. ∎ So x ↦ x^N permutes F_p and F_q (zero included). A uniform field input stays uniform on the field; a uniform
unit input becomes uniform on the field with one point omitted. Every prime partner gives a unit r for E = N. A node (u)^N = w with a unit
exponent has the unique solution u = w^{1/r} (e.g. (x^N + c)^N = 1 at p ⇔ x^N + c = 1; codex (4)).

## 3. Fresh coins: Theorem R, the inverse-fiber bound, and all-depth chains

**Theorem R (modal decoupling). PROVED.** Let T be a finite adaptive tree with zero/gcd tests whose nodes may use
coins drawn independently and exactly uniformly from Z/N (or A/NA), and whose leaf outputs are fixed functions of
the transcript. (Treat the coins as one tape uniform in (Z/N)^m; which tape entry a node uses is determined by the
transcript prefix, so sharing and reuse across nodes are allowed.) Write a_ℓ(p|q) := Pr_coins[τ_p = ℓ | p, q]. By CRT the
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

*Meaning.* For deterministic programs (no coins) R reduces to PT-X, with 4s in place of PT-X's 2s. With coins, every bit of non-split information beyond the
modal transcript is paid for by splitting, with constant 4. Rare both-pass events such as Monier liars (Monier, TCS 12 (1980) 97–108; Rabin, J. Number Theory 12 (1980)
128–138) are absorbed
into s, since κ_pκ_q ≤ s when κ ≤ ½. What remains is the *modal* deviation: a node whose typical answer at p, given q,
differs from a p-only reference. **Prior art (literature critic):** the mechanism is known. Aggarwal–Maurer
(EUROCRYPT 2009, Lemmas 3, 7, 8) show that a non-extreme answer probability forces splits (η ≥ δ²) and follow the
typical-answer path; Leander–Rupp (ASIACRYPT 2006, Lemma 1) bound simulation failure by the split probability. R is a
challenge-free, Bayes/factor-bit reformulation with an explicit constant and a per-prime reference; we claim no more.

**Lemma CB (inverse-fiber bound; codex prim_1 (23)). PROVED.** For X uniform on F_p^* and a chain
((X^{e₁} + c₁)^{e₂} + c₂ ⋯)^{e_k}, with positive exponents and arbitrary shifts, Pr[chain = w] ≤ ∏_i gcd(e_i, p − 1)/(p − 1).
For a coin uniform on F_p (a coin from Z/N), add 1/p.
Interspersed nonconstant polynomial maps of degree D_i multiply the numerator by ∏ D_i. (Shifts are injective; a
power map has fibres of size ≤ gcd(e, p − 1).) Tested: `test_collision_bound`.

**Corollary R-chain. PROVED.** Consider a tree of height h whose nodes are pure nodes and *fresh-base chain probes*
(A·chain − B with a fresh coin at the bottom; A ≢ 0 at p; A, B, the shifts and the coefficients fixed before the
coin; interspersed polynomial maps nonconstant mod p; affine menu exponents E_i; ≤ k gates; D the product of the
polynomial degrees). Take κ₀(p)
to be the path on which every probe answers its *formal value at p* ("nonzero" unless A ≡ B ≡ 0 or the probe is
constant at p, e.g. A ≡ 0 (mod p), in which case its deterministic answer) and pure nodes answer as at p. Then:
- m_p ≠ κ₀(p) forces Pr_coins[τ_p ≠ κ₀(p)] ≥ ½, hence some probe on κ₀(p)'s path has local zero probability ≥ 1/(2h);
- by CB, that probe has G_v(p, q) := ∏_i gcd(E_i(N), p − 1) ≥ (p − 1)(1/(2h) − 1/p)/D;
- Pr_q[G_v ≥ T] ≤ Σ_{i ≤ k} Σ_{d|n, d ≥ T^{1/k}} Pr_q[d | α_i q + β_i] ≤ kα τ(n)(|I|T^{−1/k} + 1)/π(I).

Hence

  Pr[correct or split] ≤ B + 4s + O(h·kα·E_p[τ(p−1)]·(2hD)^{1/k}·X^{−1/k}·ln X)

(iid uniform on distinct pairs, or add 1/π(I) for p = q; the error is summed over the reference path). Here E_p[τ(p − 1)] = O(log X) (Titchmarsh divisor problem / Brun–Titchmarsh summation). For fixed depth k
this is X^{−1/k+o(1)} per probe. Separately, as a statement about usefulness rather than negligibility, it becomes
vacuous at k ≈ log X / log log X. For E = N gates every gcd is 1 (Lemma P), so pure-N chains have exact
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

**Lemma KS (sparse structure; codex prim_1 (12); the degree-compression method of Bi–Cheng–Rojas, ISSAC 2013 /
SIAM J. Comput. 45(4) (2016) 1433–1447, Thm 1.1, and Kelley 2016 Thm 2.2; CFKLLS 2000). PROVED.**
Let G ≤ F_Q^* be cyclic of order n and f a t-sparse Laurent polynomial (t ≥ 2) with nonzero coefficients and exponents distinct mod n. For every
1 ≤ L ≤ n, the roots not lying in a zero coset of a subgroup of order ≥ L number at most

  |R_L| ≤ 2n/(⌊n/L⌋ + 1)^{1/(t−1)} < 2n^{1−1/(t−1)}L^{1/(t−1)}.

The removed zero cosets can be taken to be ≤ ⌊n/L⌋ disjoint cosets of a single subgroup of order ≥ L (the
*selected* union V*_L). This differs from the *maximal* union V_L of all zero cosets of order ≥ L, computed by
`coset_union`, which need not consist of cosets of one subgroup: (T⁶ − 1)(T⁵ − 1) over F_31 with L = 5 (codex).
*Proof sketch* (codex's; normalise a₁ = 0 by a monomial; if δ = gcd(n, a₂, …, a_t) ≥ L every root is removed):
Dirichlet gives 1 ≤ e < ⌊n/L⌋ + 1 with all |e·a_i mod n| ≤ M := n/(⌊n/L⌋+1)^{1/(t−1)}.
The map x ↦ x^e sends G onto H of order n/gcd(e, n) ≥ L. On each coset bH, f(bx^e) is a Laurent polynomial with
exponents in [−M, M]: either zero (a zero coset) or with ≤ 2M roots in x, i.e. ≤ 2M/gcd(e, n) on bH. ∎
- *Partner congruences* (Kelley Thm 2.2, verified from the arXiv source): if f vanishes on bH_C then grouped
  coefficients cancel mod C, so every exponent has a partner ≡ mod C. The ratio of the partnered gate monomials has
  order dividing n/C: an order alignment, (M). Alignment is necessary; the coefficient cancellations are additionally required.
- (v0 had a weaker recursive form, 2⌊log₂t⌋n^{1−1/(e(1+⌊log₂ t⌋)(t−1))}, still valid; superseded.)
- Tested by brute force: `test_lemma_KS_bruteforce`, `test_lemma_KS_compression_bruteforce`.

**Theorem PS. PROVED.** For such a node, reference "formally nonzero", and L = n^θ (0 < θ < 1):

  Pr_q[v deviates at p] ≤ (2 + o(1))θ^{−1}·(k/φ(k))·dens(V_L(p)) + 4n^{1−(1−θ)/(t−1)}/π(I),

where V_L(p) here denotes the selected union V*_L of the single subgroup supplied by KS, i.e. q-classes
r ≡ r₀ (mod k) for one k | p − 1, k ≤ n^{1−θ}. The reference is "formally nonzero at p" (coefficients reduced and
exponents merged mod p − 1).
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
- *Node at p:* F(r) = A(r)·u(r)^{αr+β} − B(r), with u(r) = y^{ar+b} + c (y pure, y ≢ 0 and c ≢ 0 mod p, a, α ≥ 1).
- *Monomial condition:* A, B are nonzero monomials in gate outputs of pure bases, so B/A = γ·g^{λr} with γ ∈ F_p^*,
  λ ∈ Z/n, g a generator. The inner value may also be u = M₁(M₂ + c) with monomials M_i, M₂ = g^{μ₂r+·}, o := ord(g^{μ₂}). The proof changes
  slightly: ind M₁ = γ₁ + λ₁r moves (αr + β)(γ₁ + λ₁r) to the right side, which becomes quadratic in r. For each fixed
  x = r mod d the number of admissible ℓ is still ≤ gcd(αx + β, d), so |S_d| and the bound are unchanged.
- Not covered: inner values with ≥ 3 terms (e.g. y^r + z^r + c) and sums of several nested gates.
- *Bound:* with o := ord_p(y^a), for every divisor d of n = p − 1,

  |Z(p)| := #{r ∈ Z/n : F(r) = 0} ≤ gcd(α, d)·τ(d)·(n/d + d·(n/o)·(√p + 2))

  (local hypotheses: y, c, and the coefficients of A, B are units mod p). A non-vacuous instance (codex): p = 65537,
  y = 3, d = 8 gives the bound 41 025 < 65 536; the observed zero counts for w = 1, 2, 3, 5, 17 are 13, 6, 2, 0, 0
  (`test_theorem_D2_nonvacuous`).

*Proof.*
1. If u(r) = 0, F(r) = −B(r) ≠ 0 (B is a nonzero monomial), so such r are never zeros; the final "+ n/o" is
   therefore unnecessary but harmless. The set {u = 0} has ≤ n/o elements (u depends on r mod o).
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
5. |Z| ≤ |S_d|·max N(x, ℓ) (u = 0 never gives a zero). ∎ (Proof critic: 2 664 numerical checks of the cell estimate
   in step 4 at primes 1 009–10 007, 0 violations, worst ratio 0.44.)

(The numerical companion `d2_bound` checks it at all 12-bit primes for D2a_w and D2d: `test_theorem_D2_bound_holds`.
At these sizes the bound is above n, so the check has no content beyond consistency.)

**Corollaries.**
- (i) *Per prime* (α = p^{o(1)}, 0 < δ < 1/4 − η/2). If o ≥ p^{1−η} and n has a divisor d ∈ [p^δ, p^{1/4−η/2}], then |Z(p)| ≤ p^{1−δ+o(1)} + p^{3/4+η/2+o(1)}.
  The deviation probability is then ≤ 2|Z|/π(I) = p^{−min(δ, 1/4−η/2)+o(1)}.
- (ii) *Exceptional primes.* The primes failing (i) have p − 1 without a prime factor in (p^δ, p^{1/4−η/2}], or have
  ord_p(y^a) < p^{1−η}. Both are multiplicative, (M)-type properties of p.
  - The first set has relative density ≤ C·δ/(1/4 − η/2) by an upper-bound sieve for shifted primes (Selberg sieve
    of dimension 1 with Bombieri–Vinogradov level; Halberstam–Richert, *Sieve Methods* (1974); Koukoulopoulos,
    arXiv 0905.0163, Thm 1.8 / Lemma 2.5 per codex; theorem numbers UNVERIFIED).
  - For a base y of bounded height (log|y| = O(1)), under GRH the second has density o(1) (Erdős–Murty-type results
    on ord_p(a): P. Erdős and M. R. Murty, "On the order of a (mod p)", CRM Proc. Lecture Notes 19 (1999) 87–97: under
    GRH ord_p(a) ≥ p/f(p) for almost all p and any f → ∞; confirmed via secondary sources only). It requires a FIXED
    non-torsion base: y = ±1 and bases growing with X are not covered. For large-height pure bases this fails (proof critic: for
    y = 2^{lcm(1..256)}, 79% of the primes in [3·2^14, 2^16) have ord_p(y) < p^{3/4}).
  - So under GRH and for bounded-height bases, for any fixed δ, all but a proportion ≤ Cδ + o(1) of p satisfy (i). Unconditionally (i) holds at
    each prime that meets its hypotheses.
- (iii) *Negligibility not established* (fixed bases, GRH). For each fixed δ the per-node deviation rate in the iid
  population is ≤ Cδ + o(1); letting δ → 0 slowly gives o(1). It is not shown to be n^{−ω(1)}, so D2 does not give a
  poly-height tree theorem. It does say that a *large* zero set
  (density ≥ p^{−δ}) at a deterministic two-layer node can occur only at primes whose p − 1 lacks medium divisors,
  or where ord_p(y) is small.

**Reformulation (codex (6)). PROVED.** For primitive g and w ≠ 0, (g^r + c)^r = w ⟺ ind(u)·ind(u + c) ≡ ind(w)
(mod p − 1) with u = g^r: the "Zech hyperbola". D2 bounds it through its reductions mod divisors d of p − 1.
- No theorem bounds it uniformly. Closest known results:
  - g^x = x: Glebsky–Shparlinski, DCC 56 (2010) Thm 2, per prime. Bourgain–Konyagin–Shparlinski (IMRN 2008,
    rnn090) treat g^h ≡ h on average over p.
  - x^x = c: Balog–Broughan–Shparlinski, Acta Arith. 148 (2011) 93–103; Cilleruelo–Garaev 2016. These also reduce
    over divisors of p − 1.
  - None covers the product-of-indices equation (codex prim_1 §1.3; literature critic).
- Step 4 of D2 is a subgroup-restricted cyclotomic-number estimate (Storer, *Cyclotomy and Difference Sets* (1967);
  Berndt–Evans–Williams, *Gauss and Jacobi Sums* (1998)).
- Exceptional-prime-free version: **OPEN**.
  - For the unit-exponent sector (s = r^{−1}), one *sufficient possible approach* (not an equivalence; codex prim_2
    #15) is a power-saving bound for Σ_{s∈(Z/n)^*} e_p(λg^{1/s} + μw^{s}).
  - That is an "exponential analogue of Kloosterman sums" and an explicitly open problem: Shparlinski, *Open problems
    on exponential and character sums* (2016), Problem 8.
  - Known partial results: Shparlinski, Rocky Mountain J. Math. 34 (2004) 1497–1502 (on average over g of order
    t ≥ p^{3/4+δ}); single-term sums Σ e_p(a g^{1/x}) by Bourgain–Shparlinski, Acta Arith. 134 (2008) 83–91.
  - None gives individual bounds for a fixed pure base.

**Lemma FZ (frozen exponents). PROVED.** Fix a bound T (independent of p; T = 64 in the code). For an integer t,
|t| ≤ T, at which all gate exponents are positive and every gate base value is a unit at p, let m_p(t) be the
lcm of the orders at p of every gate-base value in the evaluation at t. (The exact freezing modulus is
lcm_i ord(b_i(t))/gcd(α_i, ord(b_i(t))); the lcm of orders is sufficient. The code uses the implemented rule
"all gate outputs equal those at t", with t ≤ 0 read through reduced exponents. For negative t the rational tower
needs invertible denominators, and zero bases need separate treatment.) For r ≡ t (mod m_p(t)), every gate output
equals its value at t, by induction along the expression, so v(r) = T_v(t) mod p. Here T_v(t) ∈ Z[1/·] is the integer
tower. So frozen zeros are:
- identities (T_v(t) = 0, at every p), e.g. 2 + 3 − 5 = 0 for S1 at t = 1, (3/2)(4/3) − 2 = 0 for S5 at t = −1,
  (2 + 1)^1 − 3 = 0 for D2a_3 at t = 1;
- or low-height coincidences p | T_v(t) (T_v(t) has height exp(O(T^{depth})) for bounded-height bases, so at most
  that many / log X primes of I).

Without the bound |t| ≤ T the lemma is vacuous: every zero r is "frozen at t = r" (proof critic).

In the q-world, r_q ≡ t (mod m_p(t)) is an (M) alignment occurring at rate ~1/m_p(t) at each p; it is rare in the
population only through small-order counting, which again needs bounded-height bases. When m_p(t) = p − 1 it is the
pin q = p + t − 1, an (A) event. leakmap's ARCH (2^N = 2^k) and the twin-prime Lucas liars are instances.

**Lemma CL (closure). PROVED.** Fix an exponent vector (k_i) independent of r and p (|k_i| ≤ 2 in the code). If at r
an intermediate gate input equals ±∏ y_i^{k_i} (pure constants y_i), the gate
output is (±1)^E ∏ (y_i^{E})^{k_i} (negative k_i: units), a monomial in lower-depth gate outputs. On the set where
this happens the selected gate is removed. For chains the depth strictly drops, and the set itself is the zero set of
a lower-depth node. In branched expressions the maximum depth need not drop. At depth 3 closure can leave a depth-2
condition with D2's exceptions. Example: a^r = −(a + 1) and
r odd give (a^r + 1)^r = (−a)^r = a + 1. At depth 2 both pieces are one-layer nodes, hence (M) by NE1/PS: e.g.
{a^r + 1 + a = 0} is one class mod ord_p(a). (For r-dependent exponents the lemma is vacuous: for primitive y every
value is ±y^{k(r)}.)
- *Control (scratch, b = 14, random bases a, unit r, frozen classes removed):* solutions of (a^r + 1)^r = w number
  2.03× the fresh-random model for w = a + 1, and 0.87–1.13× for w ∈ {a + 2, random, (a + 1)², (a + 1)^{−1}}.
- This is exactly the excess seen in E1 (§9).

**Shift depth ≥ 2, deterministic.**
- *Reducible cases (codex prim_2 #18; PROVED).* If the outer exponent is a unit mod p − 1 (on prime partners, E = N
  gives r a unit by Lemma P) and the target w has small order h, peeling gives
  ((g^r + c₁)^r + c₂)^r = w ⟺ (g^r + c₁)^r + c₂ ∈ {w^s : s ∈ (Z/h)^*}, a finite union of depth-2 conditions.
  For w = 1: (g^r + c₁)^r = 1 − c₂.
  - The registered D3 (c₁ = c₂ = 1, w = 1) reduces on prime partners to g^r = −1, i.e. a one-layer condition.
  - Codex checked the equivalence in 2 144 instances.
- *General case: OPEN.* For a large-order target, D2's mechanism needs the joint law of (r mod d, ind(u) mod d) with
  inner value u = (y^{ar} + c₁)^{e(r)} + c₂. Its multiplicative characters are not characters of an algebraic function
  of y^{ar} of bounded degree, so no Weil or Jacobi-sum estimate applies. Shift depth 2 with a large-order target is
  the boundary of the present general method. Multi-term inner values at depth 2 also lie beyond it.
- Experiments on D3 show the FROZEN/CLOSURE + fresh-random behaviour (§9). Codex's genuine two-shift test
  ((2^N+1)^N+2)^N − 3, b = 12: 1 split, no two-sided zero.

## 6. Jacobi symbols of computed values

(c|N) = χ_p(c_p)χ_q(c_q) never splits unless c ≡ 0 at a prime, and then gcd splits. PT-X and R do not apply.

**Theorem J1 (fresh-base gate polynomials). PROVED after repair** (core estimate codex prim_1 (21)–(22); v1 was
WRONG as stated: codex prim_2 #1 gave a counterexample, and the proof critic flagged the same gap).

*Counterexample to v1 (codex).* P(T) = a T² + b with b the product of all sampled primes and a chosen by CRT with
prescribed Legendre signs. P is non-degenerate over Q but ≡ aT² modulo every sampled prime. (P(x^N)|N) then leaks ε
with success 0.750 (v1's right side: 0.509), also on all primes of (1000, 2000]. Its coefficients have height
≈ ∏ primes, an exponential CRT construction like PT-J's CRT c.
- *Setting:* population N ≡ 1 (mod 4), e = ±1 equally likely, p ≠ q iid μ_e (distinct pairs); β is taken w.r.t. the equal mixture
  ½(μ₊ + μ₋). There are k **non-adaptive** Jacobi tests J_i = (P_i(x^{E(N)}) | N) (in an adaptive tree, let k count all
  Jacobi nodes), with x a fresh coin uniform on Z/N (shared), P_i ∈ Z[T] of degree ≤ D, and E an affine menu
  exponent. For S ⊆ [k] put P_S = ∏_{i∈S} P_i and g_x = gcd(E(N), x − 1).
- *Degeneracy is defined mod x:* S is degenerate at x iff P_S ≡ a T^j Q(T)² over F̄_x. Let a_S be the pure integer
  for S degenerate over Q, and let Bad be the event that p or q is a prime at which some S is degenerate mod the prime
  but not over Q. Such primes divide a fixed nonzero integer (the leading/constant coefficients and the discriminant
  of the squarefree part of P_S) of log-size O(kD·(log H + log(kD))), where H bounds the coefficients of the P_i. So
  Pr[Bad] ≤ 2^{k+1}·O(kD(log H + log kD)/log X)/π(I), which is negligible for log H ≤ X^{1−ε}. **The height
  hypothesis is necessary**, by the counterexample above (log H ≈ π(I) log X).
- *Double zeros:* let dz = Pr[some J_i = 0 with no split] (P_i(x^E) ≡ 0 at both primes). Independence of the coin's
  components gives dz ≤ k·E[(Dg_p/(p−1) + 1/p)(Dg_q/(q−1) + 1/q)].
- *Bound:*

  Pr[correct or split] ≤ ½ + (2^{k−1} + 1)(s + dz) + Pr[Bad] + Σ_{S≠∅, degenerate over Q} |β(a_S)β(−a_S)|
                         + ½Σ_{S≠∅} E[min(1, (kD√p+1)g_p/(p−1))·min(1, (kD√q+1)g_q/(q−1)) + 2kD(g_p/(p−1) + g_q/(q−1)) + 2/p + 2/q].

- *Proof outline:*
  - Fourier over {±1}^k on the event that no J_i is 0. The event "some J_i = 0" has probability z_e = s_e + dz_e
    under e. This gives Pr[correct ∧ all nonzero] ≤ ½(1 − z) + ¼Σ_{S≠∅}|E_+[J_S] − E_−[J_S]| + 2^{k−2}(z_+ + z_−).
  - E_e[J_S | p, q] = m_{S,p}m_{S,q} with m_{S,x} = E_{coin} χ_x(P_S(h)), where h = x^E is uniform on the gate image
    H_x (|H_x| = n_x/g_x) except with probability 1/x (coin ≡ 0).
  - Non-degenerate at x: |m| ≤ (kD√x + 1)/|H_x| + 1/x (Weil over the subgroup; codex (21)).
  - Degenerate: m_{S,x} = χ_x(a_S)·ν + O(kD/|H_x| + 1/x), with ν = 1 iff j is even or g_x is even (H_x ⊆ squares).
    Since g_x ≡ E(N) ≡ α + β (mod 2) at both primes, ν is a public constant, and the term is ν·m_e(a_S)². PT-J
    gives ¼|m₊² − m₋²| = |β(a_S)β(−a_S)| for the equal mixture.
- *Classification:* a fresh-base gate-Jacobi leak of ε is either an (M) alignment (g_p g_q ≳ √(pq)/poly; for affine
  menus this has probability X^{−1/2+o(1)}), a split/double-zero event, or (Q) via PT-J for a low-height pure a_S.
  The last is negligible by PT-J1 (log|a_S| = O(kD·height)).
- Codex's correction is built in: P(T) = T on a subgroup of squares has bias 1. That is the j-odd, g-even case,
  i.e. public parity.

**Deterministic bases. OPEN.**
- Theorem U does not control signed averages. At p, (2^N − 1 | p) = (M_q | p) with M_q = 2^q − 1, so
  (2^N − 1 | N) = (M_q | p)(M_p | q).
- At a single modulus, sums Σ_{n≤T} Λ(n)χ(g^n + a) are bounded by Banks–Friedlander–Garaev–Shparlinski, "Exponential
  and character sums with Mersenne numbers", J. Aust. Math. Soc. 92 (2012) 1–13, Thm 3.3. (v1 miscited the additive
  paper Banks–Conflitti–Friedlander–Shparlinski, Compositio 140 (2004) 15–30.) Its useful range needs summation
  lengths well beyond the balanced partner range (codex), so the single-modulus input is already insufficient; the
  two-prime coupled product has no result we know of.
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
| fresh coins, chains of fixed depth k | **PROVED** (R + CB) | 4s + X^{−1/k+o(1)} per probe | k ≈ log X/log log X (CB is lossy) |
| fresh coins, polynomials in several chain outputs | R reduces to local zero-prob ≥ 1/(2h) | — | OPEN beyond sparse one-layer / bounded-degree Weil |
| fresh coins, precise two-layer rates | **PROVED** (W2, codex) | √p/\|H₁\| | depth 3: composed degree ~ p, Weil trivial |
| deterministic, 1 gate, low degree | NE1 (imported) | D·μ(o) | — |
| deterministic, 1 layer, t-sparse sums | **PROVED** (PS) | (M)-part + 4n^{1−(1−θ)/(t−1)}/π(I) | t ≈ log X/log log X |
| deterministic, 1 layer + exponential-degree pure circuit | reduces to pure root density (U) | — | OPEN (slpcirc τ-type) |
| deterministic, depth 2, binomial inner value, monomial A, B | **PROVED** at good primes (D2); w = 1 √p bound (codex) | p^{−δ} at good p; Cδ + o(1)-density exceptions (GRH, fixed non-torsion base) | exceptional primes (rough p − 1, small ord); large-height bases |
| deterministic, depth 2, ≥ 3-term inner value or sums of nested gates | OPEN | — | Weil degree ~ p |
| deterministic, shift depth ≥ 2, unit exponent, small-order target | reduces to depth 2 (peeling, codex) | as depth 2 | — |
| deterministic, shift depth ≥ 2, large-order target | OPEN | — | no character-sum structure ("doubly Zech") |
| Jacobi, fresh base, coefficient height log H ≤ X^{1−ε} | **PROVED** (J1, repaired) | (M) + PT-J(a_S) + Pr[Bad] | SLP-sized coefficients (codex counterexample) |
| Jacobi, deterministic gates | OPEN | — | signed averages over prime partners |
| Jacobi, pure compressed c | PT-J criterion; OPEN | — | character sums of huge conductor |
| extension components F_{p^k} | OPEN | — | short intervals (already at depth 1) |

So for **deterministic** bases the present general method stops at depth 2 for multi-term inner values or sums of
nested gates, and at shift depth 2 with a large-order target. Small-order targets with unit exponents peel down a level.
Binomial depth-2 nodes are bounded only off an exceptional set of primes, without established negligibility, and only
for fixed bases. For **fresh** bases nothing breaks at a fixed depth for chains.
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
- **Maurer–Wolf** (EUROCRYPT 1998, LNCS 1403, 72–84) and **Boneh–Venkatesan** (EUROCRYPT 1998, LNCS 1403, 59–71).

The tasks differ: the RSA papers (Leander–Rupp, Aggarwal–Maurer, Dachman-Soled–Loss–O'Neill, Damgård–Koprowski)
concern inverting a challenge; Jager–Schwenk concerns computing Jacobi/QR generically; Rotem–Segev concerns
accelerating forward evaluation; Boneh–Venkatesan concern algebraic (straight-line) reductions from factoring to
low-exponent RSA. None computes a factor bit from N alone. Further SLP-themed references: D. R. L. Brown, "Breaking RSA
may be as difficult as factoring", J. Cryptology 29 (2016) 220–241; U. Maurer, "Abstract models of computation in
cryptography", IMA C&C 2005, LNCS 3796, 1–12. The mechanism of Theorem R is Aggarwal–Maurer's (Lemmas 3, 7, 8) and
Leander–Rupp's (Lemma 1). Their
"component-separating zero ⇒ gcd" engine is the diagonal lemma behind PT-X/R. Exponential congruences:
Glebsky–Shparlinski (DCC 56 (2010) Thm 2), Balog–Broughan–Shparlinski (Acta Arith. 148 (2011)), Cilleruelo–Garaev
(Proc. AMS 144 (2016) Thms 1–2), Kurlberg–Luca–Shparlinski (MRL 22 (2015)), Bourgain–Konyagin–Shparlinski (IMRN 2008).
Sparse roots: Kelley (LMS JCM 19A (2016) 196–204, Thms 2.2–2.3; checked against the arXiv source), Bi–Cheng–Rojas
(ISSAC 2013; SIAM J. Comput. 45(4) (2016) 1433–1447, Thm 1.1), CFKLLS (Israel J. Math. 120 (2000) 23–46, Lemma 7).
Modular hyperbolas: Shparlinski, Jpn. J. Math. 7 (2012) 235–294.

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
| P8 Jacobi, \|z\| < 4 | 9 defined cells pass; J7_14, J7_15 undefined | 8 defined cells pass; J7_14, J7_15, J7_16 undefined | all 11 defined, pass (max \|z\| 3.05, J7_19) | **pass at b = 16; not established for the undefined cells** |
| P8 positive control J_CRT z > 10 | 51.8 | dyadic estimate ≤ 0: undefined | — | **pass at b = 12; calibration not established at b = 14** |
| ACCIDENTAL splits | 0 | 0 | 0 | **pass** |

Non-split Bayes gains (Pr[correct ∧ no split] − (1 − s)B, all pairs):
- NEST battery: |gain| ≤ 0.0022 (b = 12), ≤ 0.0015 (b = 14), ≤ 0.0006 (b = 16). s = 0.024, 0.008, 0.002.
- SUM battery: ≤ 0.0014. LUC: ≤ 0.0004.
- CTRL/ALL: 0.02–0.07 on ε and [p ≡ 1 mod 3]. This comes from the pure p − 1 controls (C_PM1, C_P), i.e. PT's
  known smooth-class correlation (factorbits §3.2), with s = 0.24–0.53.

**Why P5/P6/P7 failed (diagnosed on b ≤ 14).**
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

**Second-stage predictions (`notes/prim_posthoc_prereg.md`).**
- Q1–Q3 were fixed before the b = 16 post-hoc run. b = 16 is a *prospectively registered reanalysis of partially
  exposed data* (codex prim_2 #21): the registered E1 summaries at b = 16 had been seen.
- Q1′–Q3′ were fixed (commit 090997f) before any computation on the **untouched** range b = 15 (787 primes, all
  309 291 pairs).
- Residual = zeros outside COSET ∪ FROZEN ∪ F1 ∪ CLOSURE, compared with the fresh-random model with shift exclusion
  on the same residues. In v2 the q-world model also masks PIN residues, as the numerator does (codex prim_2 #20);
  the b = 14/16 numbers below are from the corrected code.

| node | b = 14 (development) full-r | b = 16 full-r (Q1: [0.8, 1.25]) | b = 16 q-world (Q2: [0.6, 1.6]) | **b = 15 full-r (Q1′)** | **b = 15 q-world (Q2′)** |
|---|---|---|---|---|---|
| D2a_1 | 0.884 | 0.955 | 0 vs 0 (exact, Lemma P) | 0.920 | 0 vs 0 |
| D2a_2 | 0.866 | 1.036 | 1.284 (48 vs 37.4) | 0.877 | 1.162 (24 vs 20.7) |
| D2a_3 | 0.932 | 0.980 | 0.830 | 0.933 | 0.678 (14 vs 20.6) |
| D2b | 0.884 | 0.955 | 0.996 (541 vs 543) | 0.920 | 0.853 |
| D2c | 1.052 | 0.986 | 0.992 | 0.953 | 1.102 |
| D2d | 0.960 | 0.974 | 0.696 | 0.992 | 0.727 |
| D2e | 0.877 | 0.882 | 0.876 | 0.848 | 0.862 |
| D3 | 0.980 | 1.006 | 0 vs 0 (reduces to 2^r = −1, §5) | 0.988 | 0 vs 0 |

- Q1, Q2 pass. Q1′, Q2′ pass on the untouched range.
- **Q3:** two-sided OTHER/OTHER zeros over all 1.1·10⁶ b = 16 pairs = 8 (all D2b) vs model 5.24. Pass.
- **Q3′:** 0 vs 2.63 at b = 15. Pass.
- D2b's 41 (b = 16) and 13 (b = 15) two-sided zeros are mostly COSET/FROZEN: base-2 Fermat-type alignments, (M).
  Summed nodes have no two-sided events at b = 15 or 16.
- The full-r ratios sit mostly just below 1 (0.85–1.04). The fresh-random model slightly overcounts the residual,
  which is the safe direction for "no hidden structure".

**Reading** (narrowed after codex prim_2 #22).
- FROZEN and CLOSURE are *exact* explanations of the zeros they remove (Lemmas FZ, CL). They leave lower-depth
  conditions.
- Agreement of the residual with the fresh-random model is *empirical evidence only*:
  - the ratio windows are descriptive, not calibrated tests;
  - the full-r totals are exhaustive finite-population counts with no sampling error;
  - a stochastic null would need the coset/shared-prime covariance structure, which we have not modelled.
- The scale is b ≤ 16 (primes < 2^16): weak evidence asymptotically.
- No unexplained structure was seen.

Deviations from the registration:
- (D1) The post-hoc classes FROZEN, F1 and CLOSURE and the shift-excluded model were added after seeing b ≤ 14 data.
  b = 16 is a partially exposed reanalysis; b = 15 is the clean confirmation.
- (D2) For ε and [p ≡ 1 mod 3] on their registered sub-populations, P2's PT-X bound is the fiber form and is
  trivially 1. The meaningful check is on all pairs, which was also run (leakmap D4 had the same issue).
- (D3) The dyadic-robust variance was centred (as in leakmap) after a b = 12 debugging run that had used uncentred
  residuals.
- (D4) The b = 12 summaries were regenerated by later code versions; the registered fields are unchanged.
- (D5) P4 in the production code checks the older recursive KS bound. The compression bound (Lemma KS as now stated)
  was rechecked by codex on the stored residuals (0 violations in 750/2 586/8 946 cases) and by the brute-force tests.
- (D6) The q-world PIN-mask bug in v1 changed the b = 16 Q2 ratios by ≤ 0.03; Q2 passes either way.

## 10. Candidate cards

| Mechanism | Model | Outcome |
|---|---|---|
| Nested N-th powers (g^N + c)^N, deterministic | PRIM depth 2 | (M) at good primes (D2); experiments: FROZEN/CLOSURE + random |
| Nested chains with fresh bases, any depth | PRIM-R | **refuted for this model** beyond 4s + modal (R + CB) |
| Sparse sums of gates | PRIM depth 1 | **refuted for this model** beyond (M) cosets (PS) |
| Gate + smooth pure power (g^N + 1)^M | depth 1 + pure | (M) via smoothness; root-density classification OPEN |
| Jacobi of fresh-base gate polynomials | PRIM-R + Jacobi | **refuted for this model** beyond (M) + PT-J (J1), for coefficient height log H ≤ X^{1−ε}; SLP-sized coefficients leak (codex counterexample, exponential CRT) |
| Jacobi of deterministic gate values (Mersenne-type) | PRIM + Jacobi | OPEN; 6 statistics null |
| Compressed pure Jacobi (Fermat numbers, 2^M − 3) | PT-J | OPEN; null |
| Three-layer towers ((2^N+1)^N+1)^N | PRIM shift depth 2 | target 1 with unit exponent reduces to depth 2 (here to 2^q ≡ −1); large-order targets OPEN |

## 11. Reviews

**Codex adversarial round (`codex/prim_2.out`, MAJOR REVISION; 25 findings), applied in v2.**
- OK: Theorem U's affine reduction, adaptive reference and 8/π(I) (#3); Lemma P (#6); Theorem R with shared coins
  and adaptivity (#7); KS/PS's compression and counting argument (#10); D2's finite-field core and |S_d| (#12) and
  its extension to M₁(M₂ + c) (#13).
- WRONG, fixed:
  - J1 false as stated: counterexample with coefficients degenerate modulo every sampled prime. J1 is now restricted
    to coefficient height log H ≤ X^{1−ε} with an explicit Bad-prime term; the hypothesis is shown necessary (#1).
  - "Equivalently" for the exponential-sum reformulation is now "a sufficient approach to the unit sector" (#15).
  - P8 misreported undefined dyadic cells; now "not established" (#19).
  - The q-world model did not mask PIN residues; code fixed and rerun (#20).
  - Mersenne citation corrected to Banks–Friedlander–Garaev–Shparlinski 2012, and the task distinctions in §8
    rewritten (#24).
- GAP/OVERCLAIM, fixed by narrowing:
  - J1's distinct pairs, coin model and double zeros (#2).
  - U's scope: affine zero/gcd trees, floor modulus, integral X (#4, #5).
  - R-chain hypotheses (A ≢ 0, coefficients fixed before the coin, nonconstant maps) and E[τ(p − 1)] (#8, #9).
  - Selected vs maximal coset union (#11).
  - D2 density corollaries need fixed non-torsion bases and do not establish negligibility (#14).
  - FZ's admissible t and denominators (#16); CL's sign and depth claim (#17).
  - Depth hierarchy rewritten around exact peeling, under which D3 reduces to depth 2 (#18).
  - b = 16 relabelled a partially exposed reanalysis, and a clean b = 15 confirmation added (#21).
  - Statistical reading narrowed (#22); verdict's gain number corrected and a non-vacuous D2 test added (#23).
  - Headline rewritten as a list of proved fragment bounds (#25).
- Codex's further leak attempt on two-shift towers found nothing (#25).

**Opus proof critic (ACCEPT WITH FIXES; 18 findings), applied in v2.**
- R, U, P, CB, W2, KS, PS and D2 verified; D2's cell estimate also numerically (2 664 checks, worst ratio 0.44).
- Fixes:
  - large-height bases defeat D2's GRH density (y = 2^{lcm(1..256)}: 79% of the 16-bit primes exceptional) (#11);
  - FZ/CL vacuous without fixed t and k_i (#12, #13);
  - J1 bookkeeping (double zeros, degeneracy mod p, non-adaptive k, equal-mixture β) (#14);
  - §0/§7 overclaims on depth and multi-term inner values (#15, #16);
  - minor: Lemma P proof, floor-gate factor, CB's 1/p, R's constant vs PT-X, D2 local unit hypotheses, the
    unnecessary + n/o.

**Opus literature critic (OK WITH FIXES), applied in v2.**
- Theorem R's novelty claim dropped; mechanism credited to Aggarwal–Maurer Lemmas 3/7/8 and Leander–Rupp Lemma 1.
- Bi–Cheng–Rojas credited for the compression method behind KS.
- Shparlinski's Open Problem 8 and the Rocky Mountain J. Math. 2004 average bound cited for the exponential
  Kloosterman analogue; BKS described as an average-over-p result.
- Citations added or fixed: the Erdős–Murty reference; the Mersenne citation; cyclotomy references for D2 step 4;
  Monier/Rabin; Maurer–Wolf and Boneh–Venkatesan pages; Brown 2016; Maurer 2005.
- Nothing found that supersedes D2 or the OPEN claims.

## 12. Reproduction (project root)
```
.venv/bin/python -m pytest -q tests/test_prim.py                 # 15 tests, ~2 s
.venv/bin/python -m src.prim.e1 --b 12                           # ~5 s
.venv/bin/python -m src.prim.e1 --b 14                           # ~1 min
.venv/bin/python -m src.prim.e1 --b 16 --pairs 200000            # ~13 min, peak RSS 0.93 GB
.venv/bin/python -m src.prim.e1_posthoc --b 14                   # ~1 min
.venv/bin/python -m src.prim.e1_posthoc --b 15                   # ~3 min (clean confirmation range)
.venv/bin/python -m src.prim.e1_posthoc --b 16                   # ~13 min
```
