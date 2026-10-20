# smoothenum: output-sensitive enumeration of doubly-smooth NFS pairs

Agent `smoothenum`, phase 15, 2026-10-20.
- Code: `src/smoothenum/` (`target.py` cost programs, `barrier.py` barrier helpers, `nfsdata.py`/`mech.py`
  mechanisms, `run_proto.py` the pre-registered runs).
- Tests: `tests/test_smoothenum.py` (15 pass).
- Pre-registration: `notes/smoothenum_prereg.md`, commit 1d4dde2, made before any registered run.
- Raw data: `experiments/raw/smoothenum/`.
- GPT-6.1-sol: blind round `codex/smoothenum_1.out`, adversarial review `codex/smoothenum_2.out` (§6).
- Reviews: Opus cost-model and literature critics (§7).

Labels:
- **KNOWN** means published, with a citation.
- **model-PROVED** means an exact consequence of the stated heuristic model.
- **HEURISTIC** means it rests on a named assumption.
- **HYPOTHETICAL** means it is conditional on an algorithm nobody has.
- **measured** means one of our runs.

## Verdict

**No algorithmic speedup was found. The single-N constant stays at MNFS 1.9019.**

**What the lever is worth (model-PROVED in heuristic H).**
- *Listing all the rational-smooth pairs at output cost* would give 1.6386, i.e.
  c³ = (20+8√6)/9. This is exactly the factory's per-N program.
- *Partial listing.* If listing costs E²·L[−θr], every fixed θ > 0 helps, and the gain is linear at the
  start, about −0.21 per unit θ. A 0.01 gain needs θ ≈ 0.047. An *absolute* sieve saving L[σ] is worth
  ≈ −0.84 per unit σ.
- *Algebraic side only.* Output-sensitive listing on the algebraic side alone would be worth far more:
  infimum (4/9)^{1/3} = 0.7631.
- *Both sides.* Listing the doubly-smooth pairs output-sensitively drives the constant to 0. For two linear
  forms it even gives heuristic polynomial time. That makes the hypothesis about as strong as polynomial-time
  factoring.

**Why every mechanism we tried fails (model-PROVED in explicit query models, plus HEURISTIC cofactor
randomness).**
- *What the mechanisms do.* They select candidate pairs through a known divisor or a known factorization:
  - (a) bottom-up products;
  - (b) multi-special-q on one or both sides;
  - (c) smooth-modulus lattices, Coppersmith/Boneh CRT decoding and meet-in-the-middle;
  - (d) factory transfer to one N;
  - (e) k-tree/birthday products and the relation-lattice reformulation;
  - (f) adaptive special-q cycle closing.
- *Proposition LQ (§3.2).* Each unit of cost yields at most P(y)·(C·u log u)^{v_box} smooth pairs, where
  v_box = log(area)/log B = O(1). That factor is L^{o(1)}, so θ = 0. The finite-size gain behind this is
  real (§3.2): about 14× at C61 (measured) and about 70× at RSA-1024 scale (computed). That explains why
  lattice sieving pays in practice, but it is not an asymptotic gain.
- *Bottom-up generation (§3.3).* Residue-blind generation loses a factor m/E = L[2/3]. Generation with a
  chosen factor reduces exactly to LQ.
- *k-tree (§3.4).* The size constraint forces k ≤ 2, which needs lists of L[2/3] size.
- *Boneh/CRT batching (§3.5).* The decoding bound caps the batched prime product at log P ≲ (log S)²/log E.
  The batch then gains only L^{o(1)}.
- *Factory transfer (§3.6).* For one N this *is* MNFS, because kN with a shared m gives MNFS's fields.
  Precomputation 2.0069 > 1.9019 is unavoidable.
- *Adaptive special-q (§3.7).* With a bounded number of tail primes per relation, it is MNFS with the head
  bound as factor base: 1.9019.

**Prototypes (CADO c61 polynomial, pre-registered).**
- The per-point heuristic that Proposition LQ rests on holds at this size.
  - Lattice cofactor yields: observed/predicted 0.82–0.99 in all 10 buckets.
  - The yield per unit cost peaks at M ≈ area (2^28 for an area of 2^31), as predicted.
  - Beyond the peak it decays with slope −0.87 bit/bit (predicted −0.90).
- The best mechanism is two-sided special-(q₁,q₂). At equal output (4436 relations) it uses 10.6× fewer
  examined points per relation than las sieves. But it is 69× slower in wall time (Python, exact tests per
  point, against las's C sieve at about 23 ns/point).
- This is a u^{O(1)} effect, consistent with the barrier, not an L-power.
- One registered band failed. The P2 point-count ratio in bucket t = 28 was 0.925, outside [0.95, 1.05].
  That is a z = −1.4 Poisson fluctuation on 367 expected points; the band was set too tight.

**What would be needed.** A way to find the *nonempty* lattices among L[2/3]-size moduli, or to correlate
residue windows mod m with multiplicative structure, at cost below one unit per candidate modulus. No
known arithmetic structure supplies either. There is also no unconditional lower bound: sieving is not
proved optimal (§3.8).

---

## 1. The target (Task 1)

**Program.** The model is heuristic H of notes/heurnfs_results.md §1. Variables:
- e: box E = L[e];
- β, β_a: factor bases;
- v: V = L[v] fields;
- δ: degree.

Derived quantities:
- r = 1/(3βδ), a = (δe + 1/δ)/(3β_a), s = max(β, β_a + v);
- relations 2e + v − r − a ≥ s;
- LA cost 2s.

Only the enumeration and test terms change (`target.py::costs`):

| Mode | Enumeration cost | Test cost |
|---|---|---|
| box sieve | 2e − σ | 2e − r + v |
| rational listing | 2e − θr | 2e − r + v |
| algebraic listing | 2e + v − a | 2e + v − a (rational tests) |
| doubly listing | 2e + v − r − a (= output) | none |

**Table 1: constants (L_N[1/3, ·], one N, RAM).**

| # | Enumeration hypothesis | Constant | Status | Optimum / note |
|---|---|---|---|---|
| T1 | box sieve (MNFS) | **1.9019** | KNOWN (Coppersmith 1993) | e = β = 0.951, β_a = 0.826, v = 0.125, δ = 1.402 |
| T2 | rational listing, cost E²·L[−θr] | c(θ): 1.9019, 1.8912, 1.8804, 1.8467, 1.7860, 1.7177, 1.6386 at θ = 0, .05, .1, .25, .5, .75, 1 | model-PROVED (HYPOTHETICAL oracle) | closed reduced form `c_rat_theta_closed` (codex eq. 5), checked against the 5-variable program |
| T3 | rational listing, θ = 1 | **1.6386** = 2((2+√6)/6)^{2/3}, c³ = (20+8√6)/9 | model-PROVED; KNOWN as the factory per-N constant (Coppersmith 1993; Bernstein–Lange 2014/921 Table 1.4) | v* = 0, β = β_a = 0.8193, δ = β^{−1/2} = 1.105, 2e = 2.0069 (= the factory precomputation) |
| T4 | gain thresholds | θ = 0.0047 / 0.047 / 0.2275 / 0.437 for gains 0.001 / 0.01 / 0.05 / 0.10 | model-PROVED (numeric) | no smallest useful θ: every fixed θ > 0 helps |
| T5 | absolute saving: box cost 2e − σ | 1.8935, 1.8852, 1.8598, 1.8179 at σ = .01, .02, .05, .1 | model-PROVED | slope ≈ −0.84 per unit σ |
| T6 | algebraic listing only | **inf (4/9)^{1/3} = 0.7631**, not attained | model-PROVED | needs δ ↑ 6β and e → ∞, with v* = 0; a run capped at e ≤ 10 gives 0.779 |
| T7 | doubly listing | **inf 0** (L[1/3, ε] for every ε) | model-PROVED | codex's family: β = β_a = h, δ = 3h, e = h + 2/(9h²), cost 2h |
| T8 | doubly listing, two linear forms of height √N | heuristic polynomial time | model-PROVED from CEP | B = (log N)^κ, box N^K, K > 1/(2(κ−1)); about N^{2K−(2K+1)/κ} doubly-smooth pairs |

**Re-derivation of 1.6386 (model-PROVED).**
- In T3 the binding set is 2e − r = 2β = c, a = β, β_a = β, v = 0.
- So e = β + 1/(6βδ), and a ≤ β becomes δβ + 1/(6β) + 1/δ ≤ 3β².
- The left side is minimal at δ = β^{−1/2}, which gives 2√β + 1/(6β) = 3β².
- With y = β^{3/2} this is 18y² − 12y − 1 = 0, so y = (2+√6)/6 and c = 2y^{2/3}.
- The precomputation is 2e = 2β + 1/(3√β) = 2.006850. This matches Bernstein–Lange's 2.006… and
  heurnfs §1. It agrees with codex to 9 digits.

**Interpretation.**
- *Asymmetry.* Output-sensitivity is worth much more on the *algebraic* side (T6) than on the rational side
  (T3). A free algebraic list lets the box grow, and the degree can then rise until the rational norm is
  cheap.
- *Doubly listing.* Output-sensitive doubly-smooth listing (T7, T8) is an implausibly strong primitive: it
  makes factoring heuristically polynomial.
- *What the rational lever needs.* The rational lever needs only θ > 0 in the exponent. Constant-factor or
  u^{O(1)} savings are θ = 0.

## 2. Mechanisms (Task 2)

**Table 2: mechanisms.** "Constant" is the single-N L[1/3] constant the mechanism can reach in the model.

| Mech. | Idea | Heuristic cost of listing the rational-smooth pairs | Constant | Barrier (§3) | Prototype (§4) |
|---|---|---|---|---|---|
| (a) bottom-up | generate smooth s, solve G(a,b) = s | residue-blind: ≥ (#found)·m/(2E), i.e. L[2/3] waste; chosen-factor variant = LQ | ≥ 1.9019 (blind: no L[1/3] at all) | Prop GEN | E3: hit rate 0.982 × predicted 2A/m (toy); c61 rate 1.1·10⁻⁸ |
| (b) multi-special-q (q₁ on G, q₂ on F, nested) | lattice of index M = ∏q | E²·L^{o(1)}: the forced divisor is ≤ area, gain (u log u)^{O(1)} | **1.9019** | Prop LQ | E2 (one side), E4 (two sides): 14.6× box per-point yield; 10.6× fewer points than las; 69× slower |
| (c1) huge smooth modulus M ≫ E² (Coppersmith-type) | CVP in Λ_M | one setup per M, expected box points E²/M; yield per cost (E²/M)·P(y/M) decays like 1/M | **1.9019** | Prop LQ (M > area branch) | E2 t = 34..40: slope −0.87 (pred −0.90) |
| (c2) Boneh CRT list decoding (batch all M dividing P) | one decode finds every n in a row with gcd(n,P) ≥ S | decoding bound log P ≲ (log S)²/log E caps P's primes at ≈ log N; hits per call E·S^{−c}, c ≥ 1/2 | **1.9019** | Prop BC | — (analysis) |
| (c3) meet-in-the-middle over CRT residues | per b, intersect q₁- and q₂-classes | per row ≥ \|Q₁\| + \|Q₂\| = sieving the row: E·B ≈ E² | **1.9019** | §3.5 | — |
| (d) factory transfer (kN, many m, shifted m) | one rational table for many related numbers | kN with shared m is exactly MNFS's V fields; shift m→m+t is a coordinate change; distinct m need distinct tables | **1.9019**; factory per-N 1.6386 only after L[2.0069] precomputation | Prop FT | — |
| (e1) k-tree / generalized birthday | product of k listed smooth numbers in the window | needs (1+lg k)/k ≥ log(m/E)/log(Em) → k ≤ 2; k = 2 lists L[2/3] | no L[1/3] | Prop W | — |
| (e2) relation lattice in (Z/m)^* | ∏p^{e_p} ≡ 1 mod m with small weight | equivalent reformulation; generic short-vector enumeration is exponential in dim = π(B) | no gain | §3.4 | — |
| (e3) sum-product / additive structure | box values form a 2-D generalized arithmetic progression (GAP); the smooth numbers are a multiplicative set | no algorithmic bridge known | — | SPECULATIVE | — |
| (f) adaptive special-q cycle closing (r21) | generate relations that close labelled cycles | bounded tail ⇒ yield = MNFS with factor base = head bound | **1.9019** | §3.7 | — |

## 3. Barriers (Task 3)

### 3.1 Black-box floor [model-PROVED; also codex §3]

Suppose the algorithm learns smoothness only by testing points, and each point is smooth independently
with probability p. Then q tests find pq smooth points in expectation, so a constant fraction needs Ω(U)
tests.

This is a floor in the black-box model only. Arithmetic access (divisibility, CRT, products) escapes it,
so the models below add that access explicitly. There is no information-theoretic Ω(E²) bound: the output
list has only O(K log U) bits.

### 3.2 Proposition LQ (lattice-query model) [model-PROVED under H_LQ + H_sel]

**Model.** The algorithm may:
- set up Λ_{M,R} for any M with known factorization into factor-base primes or ideals, with R the root it
  induces (cost ≥ 1 per setup, even if the lattice is empty);
- enumerate its box points (cost 1 each);
- test them (cost L^{o(1)}).

Adaptivity is allowed.

**Assumptions.**
- *H_LQ.* For points of Λ_M, G/M behaves like a random integer of its size.
- *H_sel.* For an M selected by the algorithm, the expected number of box points is U/M·(1+o(1)) + O(1).
  The selection is blind to the box hits; the batched selection of §3.5 is treated separately.

**Claim.** Expected smooth pairs found, counted with multiplicity, satisfy ≤ cost × P(y)·g*, with
g* ≤ (C·u log u)^{v_box}, v_box = log U/log B and u = log y/log B. Here g* = L^{o(1)}.

**Proof.**
- *Gain per point.* A point of Λ_M with M = B^v has success probability P(y/M) ≈ ρ(u − v). By de Bruijn,
  ρ(u−1)/ρ(u) ~ u·ξ(u) with ξ ~ log u. So ρ(u−v)/ρ(u) ≤ (C·u log u)^v.
- *Yield per unit cost.* For v ≤ v_box the yield per unit cost is ≤ P(y)(C·u log u)^{v}. For v > v_box it is
  ≤ (U/M)·P(y/M) ≤ P(y)·(C·u log u)^{v_box}·((C·u log u)/B)^{v−v_box}, which is smaller because
  B ≫ u log u.
- *Two sides.* For two-sided lattices (M = M_r·M_a), the same argument applies with v_r log B_r + v_a log B_a
  ≤ log U.
- *Size of g*.* At the MNFS optimum, v_box = 2e/β = 2 and u ~ (x/log x)^{1/3}. So g* = (log N)^{O(1)} =
  L^{o(1)}.
- *Consequence.* Any LQ algorithm needs cost ≥ (#rational-smooth pairs needed)/(P(y)·L^{o(1)}). That is the
  box cost for the same output. In Table 1 this means σ = θ = 0 and c = 1.9019. ∎

**Finite size.** `barrier.lq_gain_table` evaluates g* at the MNFS optimum.

| bits of N | 512 | 768 | 1024 | 2048 | 4096 | 16384 | 2^16 | 2^20 |
|---|---|---|---|---|---|---|---|---|
| u_rat | 3.49 | 3.84 | 4.11 | 4.88 | 5.84 | 8.51 | 12.6 | 28.8 |
| g* (upper bound, with multiplicity) | 36 | 53 | 69 | 126 | 228 | 727 | 2246 | 19672 |
| σ_eff = log g*/L-unit | 0.156 | 0.144 | 0.135 | 0.115 | 0.096 | 0.066 | 0.045 | 0.019 |
| fraction of the rational penalty | 0.88 | 0.81 | 0.76 | 0.64 | 0.52 | 0.35 | 0.22 | 0.09 |

Reading the table:
- σ_eff → 0, but slowly: like (log x)/x^{1/3}.
- At practical sizes the o(1) is not small: at 1024 bits, g* would offset 76% of log(1/ρ(u_rat)).
- *But g* is an upper bound counted with multiplicity.*
  - A smooth n is counted once per divisor M ≈ U that it has.
  - Visiting every useful modulus costs about U·ρ(v_box)·(1+1) units.
  - So the gain in *distinct* yield is far smaller. Choosing M = q₁q₂ with q_i ∈ [B/2, B] gives a gain
    of about polylog(B).
- This finite-size regime is why lattice sieving and special-q work so well in practice. It is not an
  L-power, and it is already exploited by las.

### 3.3 Proposition GEN (bottom-up) [PROVED: injectivity and density; HEURISTIC: blind equidistribution]

**Unique preimage [PROVED].**
- The solutions of Y₁a + Y₀b = s are (a₀ + kY₀′, b₀ − kY₁′).
- If |Y₀| > 2A, at most one solution has |a| ≤ A. G is therefore injective on the box, and the box values
  have density U/(2·max|G|) ≈ E/m.
- A generator of smooth s that is blind to s mod Y₀ hits the box with that probability. Its cost per
  rational-smooth pair found is m/E = L[2/3, 1/δ]. This is worse than the box by a factor L[2/3].
- Generating all smooth s ≤ Em costs Ψ(Em, B) = Em·L[−r] ≫ E².

**Reduction to LQ [PROVED].**
- Fix a smooth factor s₁ and let s₂ range over the values that put s₁s₂ in the box.
- The admissible s₂ are in bijection with the box points of Λ_{s₁}: s₁ | G(a,b) ⇔ (a, b) ∈ Λ_{s₁}.
- So chosen-factor bottom-up generation *is* LQ.

**Measured (E3, toy N = 64 bits, d = 3, m = 60311).** The hit rate was 0.016675, against 2A/m = 0.016979
predicted (ratio 0.982). At c61 the rate is 2A/|Y₀| = 1.09·10⁻⁸.

### 3.4 Proposition W (k-tree / birthday) [model-PROVED size accounting; HEURISTIC algorithm class]

**Size accounting.**
- Wagner's k-tree needs bottom lists of size ≥ (m/E)^{1/(1+lg k)} to hit a window of E out of m residues.
- A product of k factors that is ≤ Em allows at most (Em)^{1/k} values per list.
- Feasibility therefore needs (1 + lg k)/k ≥ log(m/E)/log(Em) = 1 − O(log E/log m).
- At L[1/3] against L[2/3], only k ≤ 2 survives (`wagner_max_k`: k = 8 at 512 bits, k = 2 from about
  2^16 bits).
- k = 2 is plain meet-in-the-middle, with lists of size ≥ (m/E)^{1/2} = L[2/3] ≫ E².

**Two further obstacles.**
- The window |s mod m| ≤ E is not compatible with multiplicative merging. One would have to pass to
  discrete logs in (Z/m)^*, where the window is an unstructured set.
- The "relation lattice" reformulation, {e : ∏p^{e_p} ≡ 1 mod m} with small weighted ℓ¹ norm, is
  equivalent to the original problem. It has dimension π(B) = L[β], and generic enumeration in it is
  exponential in the dimension.

### 3.5 Proposition BC (Boneh CRT batching) [HEURISTIC; uses Boneh 2000's decoding bound]

**Boneh's bound.** Boneh's CRT list decoding finds, in one polynomial-time call, every n in an interval of
length H whose gcd with a fixed P is ≥ S, provided log H ≲ (log S)²/(D log P).

**What this allows here.**
- Use one row (D = 1, H = 2E). Then log P ≲ (log S)²/log E. Since log S ≤ log y ~ x^{2/3}, this gives
  log P ≲ x^{1+o(1)}. So P can contain only primes up to y_P ≈ log N.
- An n with a y_P-smooth part ≥ S has probability ≈ ρ(log S/log y_P) = S^{−c} with c ≥ 1/2 − o(1) in every
  regime (c ≈ 2/3 when S is L[2/3]-sized).
- So a call returns about E·S^{−c} candidates, each with cofactor y/S.
- The yield per unit cost is then ≤ P(y)·(C·u log u)^{log S/log B}·min(1, E·S^{−c}).
- That is maximised at S ≈ E^{1/c} = L[1/3], i.e. at v = O(1), so the gain is L^{o(1)}.

**What would break it.** Making c → 0 needs P to contain the factor base, i.e. log P ≈ B. Boneh's bound
forbids that. Codex reached the same conclusion (`smoothenum_1` §2(c), eq. 7).

**Meet-in-the-middle over CRT residues (c3).**
- For each b and q₁ > 2E there is at most one admissible a. Finding the q₁ that have one in row b is
  sieving that row.
- Batch remainder trees cost (E + |Q|)·polylog per row.
- Total E·(E + B) ≈ E² at e = β.

### 3.6 Proposition FT (factory transfer to one N) [PROVED]

1. *Precomputation cannot be avoided.* Every algorithm that uses a rational table built by box sieving
   pays 2e ≥ 2.0069 at the factory optimum. Amortising it needs ≥ L[0.368] targets, and the original N
   still pays the whole table.
2. *kN is MNFS.* The base-m polynomials f_k of kN for many k with the same m all have the root m mod N.
   Relations from f_k combine mod N exactly as MNFS's fields do: the rational side is shared and the
   algebraic columns are pooled. The program is the MNFS program, so c = 1.9019.
3. *Shifted m is a coordinate change.* m → m+t with f(X) → f(X − t) and a → a + tb preserves both norms
   (codex). Distinct m with distinct f need distinct rational tables. Their sheared boxes overlap only in
   |b| ≤ 2E/|t − t′|, so there is no sharing.
4. *A universal m (say m = 2^t) is the factory itself.* One m serves every N in [m^d, m^{d+1}]. It gives
   1.6386 per N only after a once-for-all L[2.0069] precomputation, which is a non-uniform algorithm
   (Coppersmith 1993, Bernstein–Lange 2014).

### 3.7 (f) Adaptive special-q collection for a smaller dependency core (r21) [model-PROVED under H with bounded tails; HEURISTIC otherwise]

**The idea.** Collect relations whose primes lie in a head of size L[η], plus at most k = O(1) tail primes
in (L[η], L[β]]. Then choose special-q among the tail primes so that cycles close.

**Bounded tails.**
- By the large-prime estimate (heurnfs row 14), the relation probability is L[−r(η) − a(η)]·u^{O(k)}. That
  is the probability for factor base L[η], up to L^{o(1)}.
- A closing relation with q a tail prime is a q-lattice point. By Proposition LQ its yield per unit cost
  is the same L^{o(1)}-corrected probability.
- So the collection program is the MNFS program with β := η, and LA costs at least 2η. The constant is
  1.9019, however the cycles are scheduled.

**Unbounded tails (k ~ u).** The relations are then just L[β]-smooth. The tail is a hypergraph, not a
graph, and the core is not provably smaller. The toy bound of codex r21 applies: c ≤ h independent tail
cycles with independent labels contain a dependency with probability ≤ 2^{c−h}. Bender–Canfield's model
predicts a linear-size core.

**Verdict.** No asymptotic lever. A provable arithmetic bias in the labels would be new; we found none.

### 3.8 What is *not* proved

- No lower bound excludes an output-sensitive rational enumerator. Rational listing is smooth-number
  enumeration in ⋃_{|b|≤E}[bm − E, bm + E], i.e. E short intervals on an arithmetic progression.
- We know no reduction showing that this is hard. The hypotheses that sit on its *algebraic* side and on
  both sides are calibrated as very strong (T6–T8).
- Propositions LQ, GEN, W and BC are barriers for explicitly modelled access, not theorems about all
  algorithms. Pomerance's warning applies (heurnfs §4.5).

## 4. Prototypes against las (Task 4; pre-registered)

**Instance.**
- N = `balanced_semiprime(200, 1501)` (C61).
- Polynomial: CADO polyselect + ropt, degree 4, skew 9594 (`experiments/raw/smoothenum/c61.poly`).
- Parameters: CADO params.c60.
- Box: |a| < 2^21, 1 ≤ b ≤ 2^9, U = 2^31.

**Memory and pressure.** Every run streamed rows or lattices. Peak RSS: calib 47 MB, E1 57 MB, E2 49 MB,
E3 45 MB, E4 50 MB, las 31 MB. Memory pressure was normal (48% free) before every run.

**Registered outcomes.**

| Pred. | Statement | Result | Verdict |
|---|---|---|---|
| P1 | box rational yield / Σ P_s(\|G\|) ∈ [0.80, 1.25] | 36774 / 40727 = 0.903 | **pass** |
| P1′ | sieve recall ≥ 0.95 | 6904/7117 = 0.970 | **pass** |
| P2 | lattice point counts / (U/M) ∈ [0.95, 1.05] for t ≤ 28 | 1.0000, 1.0000, 1.0000, 1.0002, 1.0011, **0.925** (t = 28) | **fail at t = 28** (340 vs 367.5 expected over 60 lattices; z = −1.43, Poisson; band too tight) |
| P2 | same ratio ∈ [0.90, 1.10] for t ≥ 31 | 0.997, 1.073, 1.081, 0.925 | **pass** |
| P3 | observed / Σ P_s(\|G\|/M) ∈ [0.75, 1.33] per bucket | 0.903, 0.870, 0.916, 0.937, 0.955, 0.818, 0.938, 0.901, 0.982, 0.990 (t = 8..40) | **pass** |
| P4(i) | argmax_t Y_t ∈ {28, 31, 34} | t = 28 (Y = 0.2425) | **pass** |
| P4(ii) | max Y_t/Y_box ≤ 1.5 × predicted max | 13.8 ≤ 1.5 × 18.1 | **pass** |
| P4(iii) | Y₃₇ < Y₃₄, Y₄₀ < Y₃₇; slope within ±0.3 of the prediction | 0.0329 > 0.0057 > 0.0009; slope −0.866 vs −0.904 | **pass** |
| P5 | bottom-up hit rate / (2A/m) ∈ [0.90, 1.10] | 0.982 | **pass** |
| P6a | ≥ 99.5% of las relations pass our exact tests | 4436/4436 | **pass** |
| P6b | las ≥ 20× faster per relation | 1.65 s vs 113.4 s wall (69×); 2.2 vs 72.4 CPU-s (33×) | **pass** |
| P6c | our two-sided special-(q₁,q₂) yield per unit cost ∈ (1, 10⁴) × box two-sided yield per point | 4.88·10⁻⁴ / 3.34·10⁻⁵ = 14.6 | **pass** |

**E2: yield per unit cost of smooth-modulus lattices (rational side).**
- Definitions: Y_t is rational relations per (points + setups). Y_box = 0.01754 per box point.
- The predicted gain uses the coprime fraction 6/π². The E1 rows (b odd) have a higher coprime fraction
  (0.818), which deflates the measured gain relative to the prediction.

| t = log₂M | 8 | 12 | 16 | 20 | 24 | 28 | 31 | 34 | 37 | 40 |
|---|---|---|---|---|---|---|---|---|---|---|
| Y_t/Y_box (measured) | 2.0 | 3.3 | 5.0 | 7.8 | 11.6 | **13.8** | 8.1 | 1.9 | 0.33 | 0.05 |
| Ŷ_t/Ŷ_box (predicted, geometry + P_s only) | 2.9 | 4.8 | 7.3 | 10.6 | 15.3 | 18.1 | 10.8 | 2.5 | 0.41 | 0.06 |

The shape is the Proposition LQ shape:
- the gain rises like P_s(y/M)/P_s(y) up to M ≈ area;
- beyond area it falls like U/M (slope −0.87 measured, −0.90 predicted, per bit).

**E4: best mechanism vs las at equal output (K = 4436 distinct two-sided relations).**

| | las (C, 2 threads) | ours: special-(q₁,q₂), Python, 1 thread |
|---|---|---|
| relations | 4436 (192 special-q, sqside 1, q ∈ [61961, 63961]) | 4436 distinct |
| examined points / relation | 9.65·10⁷/4436 = 21 750 | (5.99·10⁶ points + 3.10·10⁶ setups)/4436 = 2049 |
| wall / relation | 0.37 ms | 25.6 ms |
| CPU per examined unit | ≈ 23 ns | ≈ 8 µs |

The two-sided double special-q examines 10.6× fewer units per relation. But each unit is a lattice
reduction or an exact two-sided test, about 350× dearer than one las sieve point. A C implementation with
batch cofactorization would need under about 230 ns per unit to tie las. The gain is the finite u^{O(1)}
gain of §3.2: two forced divisors on norms of 56 and 100 bits. It does not scale to an L-power. (Composite
and hybrid special-q are KNOWN: Boudot et al. 2020 §4; codex §4.)

**Disclosed deviations.**
- las needed `-dup-qmin 61961,61961 -dup-qmax 63961,63961` together with `-dup`, and `-sqside 1` was made
  explicit (it was the default). No other change.
- In E4 the `tests` counter in `e4_ours_summary.json` was overwritten by a key collision (it reads 0). The
  number of exact tests equals the number of coprime points, about 0.61 × 5.99·10⁶. The code was fixed
  after the run.
- Only 979 of the 4436 las relations lie inside our box, because las's q-lattices reach |a| ≈ 5·10⁶. The
  point-per-relation comparison is therefore between two different regions of the same size class.

## 5. Status of the r20 §5 connection

- r20 #5 (search entropy against norm height) needs a family where many search coordinates share one
  height cost.
- Output-sensitive enumeration is the *enumeration-side* analogue: the cost of a search would have to
  track the number of hits instead of the entropy.
- Both run into the same wall. The available arithmetic structure (divisibility lattices, products, CRT)
  only reduces the norm by known factors of total size ≤ the search area, so the gain is
  ρ(u − O(1))/ρ(u) = L^{o(1)}.

## 6. GPT-6.1-sol

**Blind round (`codex/smoothenum_1.out`).** It was run before our analysis. Agreements:
- every constant: 1.9019, 1.638586234 (closed form (20+8√6)/9 cubed-root), the algebraic infimum
  (4/9)^{1/3}, the doubly infimum 0, and the θ thresholds 0.047 / 0.2275 / 0.437;
- the ambient m/E waste for bottom-up generation;
- the coverage argument for special-q;
- the Boneh decoding-bound failure;
- factory transfer = MNFS;
- the black-box floor.

Its additions, adopted:
- the same prime on both sides forces p | f(m) = N;
- the reduced closed form for c(θ) (verified in `test_theta_closed_form_endpoints_and_numeric`);
- the m → m+t coordinate-change remark.

Our additions:
- the explicit LQ proposition with the finite-size g* table;
- the k-tree size inequality;
- the Boneh S^{−c} hit-rate accounting;
- (f);
- all measurements.

**Adversarial review (`codex/smoothenum_2.out`).** See §7.3.

## 7. Reviews

(To be filled after the critics: §7.1 cost-model critic, §7.2 literature critic, §7.3 codex adversarial.)

## 8. Reproduction

```bash
.venv/bin/python -m src.smoothenum.target          # Table 1 (~4.5 min)
.venv/bin/python -m src.smoothenum.barrier         # LQ g* table, k-tree table
.venv/bin/python -m pytest -q tests/test_smoothenum.py
# instance (CADO build of cadobench, read-only):
B=~/src/cado-nfs/build/kmini; cd experiments/raw/smoothenum
$B/polyselect/polyselect -N <N> -degree 4 -P 420 -admax 10000 -incr 60 -nq 64 -keep 10 -t 2 > c61_raw.txt
$B/polyselect/polyselect_ropt -inputpolys c61_sizeopt.poly -ropteffort 0.1 -I 10 -Bf 111342 -Bg 78682 -t 2
$B/sieve/makefb -poly c61.poly -lim 111342 -maxbits 10 -out c61.roots1.gz -t 2
$B/sieve/las -poly c61.poly -fb1 c61.roots1.gz -lim0 78682 -lim1 111342 -lpb0 18 -lpb1 19 -mfb0 17 -mfb1 38 \
   -ncurves0 2 -ncurves1 2 -I 10 -q0 61961 -q1 63961 -sqside 1 -t 2 -dup -dup-qmin 61961,61961 \
   -dup-qmax 63961,63961 -out c61_las.rels
cd ../../..
for c in calib e1 e2 e3; do /usr/bin/time -l .venv/bin/python -m src.smoothenum.run_proto $c; done
.venv/bin/python -m src.smoothenum.run_proto geo
.venv/bin/python -m src.smoothenum.run_proto e4check experiments/raw/smoothenum/c61_las.rels
/usr/bin/time -l .venv/bin/python -m src.smoothenum.run_proto e4ours 4436 3600
```
