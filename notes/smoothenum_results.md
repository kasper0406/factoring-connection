# smoothenum: output-sensitive enumeration of doubly-smooth NFS pairs

Agent `smoothenum`, phase 15, 2026-10-20.
- Code: `src/smoothenum/` (`target.py` cost programs, `barrier.py` barrier helpers, `nfsdata.py`/`mech.py`
  mechanisms, `run_proto.py` the pre-registered runs).
- Tests: `tests/test_smoothenum.py` (15 pass).
- Pre-registration: `notes/smoothenum_prereg.md`, commit 1d4dde2, made before any registered run.
- Raw data: `experiments/raw/smoothenum/`.
- GPT-6.1-sol: blind round `codex/smoothenum_1.out`; adversarial review `codex/smoothenum_2.out` (§7.3).
- Two independent Opus critics, one on the cost model and one on the literature (§7.1, §7.2).
- This note is the revised version: every fix from the three reviews has been applied.

Labels:
- **KNOWN** means published, with a citation.
- **model-PROVED** means an exact consequence of the stated heuristic model and the named assumptions.
- **HEURISTIC** means it rests on a named assumption.
- **HYPOTHETICAL** means it is conditional on an algorithm nobody has.
- **measured** means one of our runs.
- "No prior statement found by us" is the strongest novelty claim made anywhere in this note.

## Verdict

**No algorithmic speedup was found. The single-N constant stays at Coppersmith's MNFS 1.9019.**

**What the lever would be worth (model-PROVED in heuristic H).**
- *Listing all the rational-smooth pairs at output cost* would give 1.6386, with c³ = (20+8√6)/9.
  - This is exactly the factorization factory's per-N program. KNOWN: Coppersmith 1993; Bernstein–Lange
    2013 §5.1; Bernstein–Lange 2014/921 §2.5.
- *Partial listing.* If listing costs E²·L[−θr], every fixed θ > 0 helps, starting at about −0.21 per
  unit θ. A gain of 0.01 needs θ ≈ 0.047.
- *An absolute sieve saving L[σ]* is worth about −0.845 per unit σ.
- *Algebraic side only.* An output-sensitive algebraic-side list would be worth far more: infimum
  (4/9)^{1/3} = 0.7631, which is half the SNFS constant.
- *Doubly-smooth listing.* An output-sensitive doubly-smooth list drives the constant to 0. With early stopping
  it gives heuristic polynomial time for two linear forms; this is a quantified form of the Schnorr–Adleman
  lattice observation. That primitive is as strong as polynomial-time factoring.

**Why the mechanisms we tried fail.**

The mechanisms:
- (a) bottom-up products;
- (b) multi- or two-sided special-q (KNOWN as a technique);
- (c) smooth-modulus lattices, Boneh/Coppersmith CRT decoding, and meet-in-the-middle;
- (d) factory transfer to one N;
- (e) k-tree products, and the relation lattice of Schnorr/Adleman;
- (f) adaptive special-q cycle closing.

The reasons:
- **Divisibility gives only an L^{o(1)} gain.** These mechanisms select candidates through a known divisor.
  In an explicit lattice-query model with conditional cofactor randomness, a known divisor no larger than the
  search area buys at most (C·u log u)^{O(1)} = L^{o(1)} (Proposition LQ, §3.2). So θ = 0.
- **The model excludes cheap identification by assumption.** The assumption that does the work is that
  nonempty lattices with M > area cannot be identified for less than one unit each. The model *assumes* this
  rather than proving it arithmetically.
  - Boneh's CRT decoding is the one known batch identifier. Within its proven decoding range it still gains
    only L^{o(1)} (§3.5).
  - Whole-box multivariate decoders are *not* excluded.
- **Practical gains are real but do not scale.**
  - We measured about 14× per unit cost on C61, which is why lattice sieving and special-q pay off in
    practice.
  - The rational-side upper bound at 1024 bits (illustrative, counted with multiplicity) is about 69×.
  - These are u^{O(1)} effects, not L-powers.
- **Bottom-up and k-tree constructions waste factors of size L[2/3]** (§3.3, §3.4).
- **Factory transfer to a single N is MNFS**: 1.9019 (§3.6).
- **Adaptive special-q gives no gain when the tails are short.** With o(u) tail primes per relation it is
  MNFS with the head bound as factor base (§3.7). With Θ(u) tail primes, cheap elimination is *not*
  excluded; it remains open.

**Prototypes.** These used a CADO c61 polynomial and were pre-registered.
- *The per-point heuristic behind Proposition LQ holds at this size.*
  - Lattice cofactor yields came out at 0.82–0.99 of prediction in all 10 buckets.
  - The yield per unit cost peaks at M ≈ area.
  - Beyond the peak it decays with slope −0.87 per bit (predicted −0.90).
- *Two-sided special-(q₁,q₂) against las at equal output (4473 relations, with las configured correctly).*
  - It examines 10.5× fewer units per relation than las sieves.
  - It yields 13.8× the exact box yield per unit.
  - It is 24× slower in wall time and 26× slower in CPU time: Python with exact tests against las's C sieve.
- *One registered band failed.* The P2 point-count ratio at t = 28 was 0.925, outside [0.95, 1.05]. The
  failure stands.

**What would be needed.**
- One of the following:
  - an identifier of the nonempty lattices among moduli larger than L[1/3], working at less than one unit
    per candidate modulus and on a scale Boneh's bound does not allow;
  - an arithmetic correlation between residue windows mod m and multiplicative structure.
- No known structure supplies either.
- Sieving is *not* proved optimal; no lower bound exists (§3.8).

---

## 1. The target (Task 1)

**Program.** The model is heuristic H of notes/heurnfs_results.md §1. Variables:
- e: box E = L[e];
- β, β_a: factor bases;
- v: V = L[v] fields;
- δ: degree.

Derived quantities:
- rational penalty r = 1/(3βδ);
- algebraic penalty a = (δe + 1/δ)/(3β_a);
- columns s = max(β, β_a + v).

Constraint and cost:
- relations 2e + v − r − a ≥ s;
- LA 2s, with ω = 2.

Only the enumeration and test terms vary (`target.py::costs`).

| Mode | Enumeration cost | Test cost |
|---|---|---|
| box sieve | 2e − σ | 2e − r + v (each rational-smooth pair tested in V fields) |
| rational listing | 2e − θr | 2e − r + v |
| algebraic listing | 2e + v − a (per field) | 2e + v − a (rational test) |
| doubly listing | 2e + v − r − a (= output) | — |

**Table 1: constants (L_N[1/3, ·], one N, RAM).**

| # | Enumeration hypothesis | Constant | Status | Optimum / note |
|---|---|---|---|---|
| T1 | box sieve (MNFS) | **1.9019** | KNOWN (Coppersmith 1993; closed-form optimum in Lee–Venkatesan 2018 Rem. 4.4, Barbulescu–Pierrot 2014) | e = β = 0.951, β_a = 0.826, v = 0.125, δ = 1.402 (the optimum is flat in β_a and δ) |
| T2 | rational listing, cost E²·L[−θr] | c(θ): 1.9019, 1.8912, 1.8804, 1.8467, 1.7860, 1.7177, 1.6386 at θ = 0, .05, .1, .25, .5, .75, 1 | model-PROVED (HYPOTHETICAL oracle); no prior statement found by us | reduced closed form `c_rat_theta_closed` (codex eq. 5), checked against the 5-variable program |
| T3 | rational listing, θ = 1 | **1.6386** = 2((2+√6)/6)^{2/3}, c³ = (20+8√6)/9 | KNOWN as the factory per-N program (Coppersmith 1993; Bernstein–Lange 2013 §5.1; 2014/921 §2.5 and Table 1.4); re-derived here | v* = 0, β = β_a = 0.8193, δ = β^{−1/2} = 1.105 (BL's μ = 1/δ = 0.905), 2e = 2.0069 |
| T4 | gain thresholds | θ = 0.0047 / 0.047 / 0.2275 / 0.437 for gains 0.001 / 0.01 / 0.05 / 0.10 | model-PROVED (numeric) | no smallest useful θ: every fixed θ > 0 helps; initial slope −0.213 per unit θ |
| T5 | absolute saving: box cost 2e − σ | 1.8934, 1.8850, 1.8598, 1.8179 at σ = .01, .02, .05, .1 | model-PROVED | initial slope −0.845 per unit σ |
| T6 | algebraic listing only | **inf (4/9)^{1/3} = 0.7631** = ½·(32/9)^{1/3}; not attained | model-PROVED; no prior statement found by us | δ < 6β_a is *necessary*; the infimum is approached as δ ↑ 6β_a, e → ∞, v = 0, β = β_a ↓ 18^{−1/3}. It is an infimum over *fixed* parameter choices: each tolerance needs a fixed large e, then N → ∞. e ≤ 10 gives 0.7786 |
| T7 | doubly listing | **inf 0** | model-PROVED (same quantifier convention) | codex family: β = β_a = h, δ = 3h, e = h + 2/(9h²), cost 2h, for each fixed h > 0 |
| T8 | doubly listing, two linear forms of height √N | heuristic polynomial time, *if* the primitive has polynomial delay / prefix-output cost, and joint smoothness and fruitfulness hold | HYPOTHETICAL; a quantified form of Schnorr 1991 / Adleman 1995 | B = (log N)^κ with κ > 1, box N^K with K > 1/(2(κ−1)): N^{2K−(2K+1)/κ+o(1)} doubly-smooth pairs. Listing *all* of them would be exponential, so the hypothesis must allow stopping after B^{1+o(1)} outputs |

**Re-derivation of 1.6386 (model-PROVED).**
- At the T3 optimum the binding set is 2e − r = 2β = c, a = β, β_a = β, v = 0.
- So e = β + 1/(6βδ), and a ≤ β becomes δβ + 1/(6β) + 1/δ ≤ 3β².
- The left side is minimised at δ = β^{−1/2}, giving 2√β + 1/(6β) = 3β².
- With y = β^{3/2} this is 18y² − 12y − 1 = 0, so y = (2+√6)/6.
- The precomputation is 2e = 2β + 1/(3√β) = 2.006850. This matches Bernstein–Lange's 2.006… (Table 1.4),
  heurnfs §1, and codex to 9 digits.

**Interpretation.**
- *The algebraic side is worth far more (T6 vs T3).* A free algebraic list lets the box grow and the degree
  rise until the rational norm is cheap.
- *Doubly listing (T7, T8) is implausibly strong.* With early stopping it would make factoring heuristically
  polynomial.
- *The rational lever needs θ > 0 in the exponent.* Constant-factor or u^{O(1)} savings give θ = 0.

## 2. Mechanisms (Task 2)

**Table 2: mechanisms.** "Constant" is the single-N L[1/3] constant reachable within the stated model.

| Mech. | Idea | Heuristic cost of listing rational-smooth pairs | Constant | Barrier (§3) / status | Prototype (§4) |
|---|---|---|---|---|---|
| (a) bottom-up | generate smooth s, solve G(a,b) = s | residue-blind: cost ≥ (#found)·m/(2E), an L[2/3] waste; the chosen-factor variant that enumerates lattice points is LQ | ≥ 1.9019 (blind: no L[1/3]) | Prop GEN | E3: hit rate 0.982 × predicted 2A/m (toy); c61 rate 1.1·10⁻⁸ |
| (b) multi-/two-sided special-q (q₁ on G, q₂ on F, nested) | lattice of index M = ∏q | E²·L^{o(1)}: the forced divisor is ≤ area, gain (u log u)^{O(1)} | **1.9019** | Prop LQ. Technique KNOWN: Pollard 1993; Franke–Kleinjung 2005; composite / hybrid special-q in Kleinjung 2006, Aoki 2008, Boudot et al. 2020 §4.1 and fn. 8 | E2 (one side); E4 (two sides): 13.8× the exact box per-point yield, 10.5× fewer units than las, 24× slower |
| (c1) huge smooth modulus M ≫ E² (Coppersmith-type) | CVP in Λ_M | one setup per M, E²/M expected box points; yield per cost (E²/M)·P(y/M) falls like 1/M | **1.9019** | Prop LQ, M > area branch | E2 t = 34..40: slope −0.87 (pred −0.90) |
| (c2) Boneh CRT list decoding (batch every M dividing P) | one decode finds every n in a row with gcd(n,P) ≥ S | decoding bound log P ≲ (log S)²/log E; candidates per row ≤ E·S^{−1/2+o(1)} (Rankin) | **1.9019** within the decoding range | Prop BC; Boneh 2002 §4 already notes the bounds do not help QS | — |
| (c3) meet-in-the-middle over CRT residues | per b, intersect q₁- and q₂-classes | per row ≥ \|Q₁\| + \|Q₂\|, i.e. sieving the row: E·B ≈ E² | **1.9019** | §3.5 | — |
| (d) factory transfer (kN, shifted m, universal m) | one rational table for many related numbers | kN with a shared m = MNFS's V fields; shifted m with its own base-m field shares ≤ 75% of an adjacent table, which is dominated by MNFS's 100% sharing; a universal m is the factory | **1.9019**; factory per-N 1.6386 only after an L[2.0069] precomputation | Prop FT | — |
| (e1) k-tree / generalized birthday | product of k listed smooth numbers in the window | for balanced, height-bounded k-trees: (1+lg k)/k ≥ log(m/E)/log(Em), so k ≤ 2; k = 2 needs L[2/3] lists | no L[1/3] | Prop W (narrowed) | — |
| (e2) relation lattice | ∏p^{e_p} ≡ 1 mod m (or the doubly version) with small weight | equivalent reformulation; generic short-vector search is exponential in dim π(B) | no gain | KNOWN: Schnorr 1991, Adleman 1995; experiments Ritter–Rössner 1997, Ducas 2021; Mulder–Sterner–van Woerden 2025 | — |
| (e3) sum-product | box values form a 2-D generalized arithmetic progression; smooth numbers are multiplicative | no algorithmic bridge known | — | SPECULATIVE | — |
| (f) adaptive special-q cycle closing (r21) | generate relations that close labelled cycles | tails of o(u) primes: yield = MNFS with the head bound as factor base | **1.9019** (tails o(u)); open (tails Θ(u)) | §3.7 | — |

## 3. Barriers (Task 3)

### 3.1 Black-box floor [model-PROVED; folklore; also codex]

**The model.** The only access is point tests, and each point is independently smooth with probability p.

**The floor.** q tests find pq successes in expectation, so a constant fraction of the smooth points needs
Ω(U) tests.

**Scope.** There is no information-theoretic Ω(E²) bound: the output has O(K log U) bits. Arithmetic
access escapes the black-box model, so the models below add it explicitly.

### 3.2 Proposition LQ (lattice-query model) [conditional formalisation; model-PROVED under H_LQ + H_sel; no prior statement found by us, the folklore being "special-q gains u^{O(1)}"]

**Model.** The algorithm keeps a history 𝓕. At each step it may do one of three things.
- *Set up* Λ_{M,R} = {(a,b) : a ≡ Rb mod M}.
  - M has known factorisation into factor-base primes or ideals; R is the root it induces.
  - Cost ≥ 1, even when the lattice is empty.
- *Enumerate* the lattice's box points, at cost 1 each.
- *Test* points, at cost L^{o(1)}.

Adaptivity is allowed. Only *new* points, not revealed before, count toward the output.

Special-q above B, or M with prime factors in (B, U], adds new matrix columns. Those lattices fall under
bounded-large-prime accounting (heurnfs row 14: gain u^{O(1)}). Equivalently, apply LQ with B replaced by U.

**Assumptions.**
- *H_LQ (conditional cofactor randomness).* Conditional on 𝓕, an unrevealed point of Λ_M has its cofactor
  G/M distributed like a random integer of that size. The two-sided version assumes conditional *joint*
  randomness of G/M_r and F/M_a.
- *H_sel (blind selection), needed only for M > U.* For every M selected with M > U,
  E[# new box points of Λ_M | 𝓕] ≤ C·U/M, with no additive constant.
  - For M ≤ U no assumption is needed: cost ≥ points, and each point succeeds with probability ≤ P(y/M).
  - H_sel is the crux. It excludes cheap identification of nonempty large-modulus lattices *by assumption*.
    §3.5 checks the one known batch identifier (Boneh) against it.

**Claim.**
- The expected number of new smooth pairs is ≤ cost × P(y)·g*, where g* ≤ (C·u log u)^{v_box},
  v_box = log U/log B and u = log y/log B.
- So g* = L^{o(1)} whenever v_box = O(1). At the MNFS optimum v_box = 2e/β = 2.

**Proof.**
1. *Dickman step.* −(log ρ)′(t) = ρ(t−1)/(tρ(t)) = ξ(t) + O(1/t), with ξ(t) ~ log(t log t) (de Bruijn
   1951; Hildebrand–Tenenbaum 1993). Hence, uniformly for 0 ≤ v ≤ u − 1,
   log(ρ(u−v)/ρ(u)) = ∫_{u−v}^{u} −(log ρ)′ ≤ v·log(C·u log u).
2. *M ≤ U.* The yield per unit cost is ≤ P(y/M) ≤ P(y)·(C·u log u)^{v}.
3. *M > U.* Under H_sel the yield per unit cost is ≤ (C·U/M)·P(y/M)
   ≤ P(y)·(C·u log u)^{v_box}·((C·u log u)/B)^{v−v_box}, and the last factor is < 1 because B ≫ u log u.
4. *Two sides.* For M = M_r·M_a, the gain is (C·u_r log u_r)^{v_r}·(C·u_a log u_a)^{v_a}, with
   v_r log B + v_a log B_a ≤ log U.
5. *MNFS.*
   - A rational-side lattice helps all V fields. An algebraic lattice helps only its own field.
   - So the yield per unit cost is ≤ V·P_r·P_a·L^{o(1)}.
   - Testing a rational-smooth pair in V fields still costs L[v].
   - The binding program is therefore exactly the box program: σ = θ = 0 and c = 1.9019. ∎

**Finite size (`barrier.lq_gain_table`).**
- The table gives the upper bound g* *with multiplicity*.
- It is evaluated at the asymptotic MNFS point with the o(1) terms dropped and continuous d. It is
  illustrative, not tuned parameters.

| bits of N | 512 | 768 | 1024 | 2048 | 4096 | 16384 | 2^16 | 2^20 |
|---|---|---|---|---|---|---|---|---|
| u_rat | 3.49 | 3.84 | 4.11 | 4.88 | 5.84 | 8.51 | 12.6 | 28.8 |
| g*_rat | 36.3 | 53.4 | 68.5 | 126 | 229 | 729 | 2246 | 19764 |
| σ_eff,rat = log g*/λ | 0.156 | 0.144 | 0.135 | 0.115 | 0.096 | 0.066 | 0.045 | 0.019 |
| log g*_rat / log(1/ρ(u_rat)) | 0.88 | 0.81 | 0.76 | 0.64 | 0.53 | 0.35 | 0.22 | 0.09 |
| u_alg (v_box,alg = 2.30) | 9.2 | 10.4 | 11.3 | 13.8 | 17.0 | 25.9 | 39.5 | 93.0 |
| g*_alg | 2.4·10³ | 3.6·10³ | 4.8·10³ | 9.3·10³ | 1.8·10⁴ | 6.3·10⁴ | 2.2·10⁵ | 2.4·10⁶ |
| σ_eff,alg | 0.338 | 0.297 | 0.271 | 0.217 | 0.174 | 0.111 | 0.071 | 0.029 |

Reading the table:
- σ_eff → 0 like (log x)^{1/3}/x^{1/3}, slowly.
- At practical sizes the o(1) is large, especially on the algebraic side. That is why las puts special-q on
  the algebraic side.
- These are upper bounds counted with multiplicity: a smooth n counts once per suitable divisor. We did
  *not* compute the distinct-yield gain of any particular covering scheme. The earlier draft's
  distinct-yield paragraph was unsupported and has been withdrawn.
- Published superpolynomial-but-L^{o(1)} gains of a different kind exist. Example: early-abort ECM's
  exp((c+o(1))(log N)^{1/6}(log log N)^{5/6}) test-cost saving (Bernstein–Lange 2014 §3).

### 3.3 Proposition GEN (bottom-up) [PROVED: injectivity and density; HEURISTIC: blind equidistribution]

**Unique preimage [PROVED].**
- Let g = gcd(Y₀, Y₁). The solutions of Y₁a + Y₀b = s are (a₀ + kY₀/g, b₀ − kY₁/g).
- If |Y₀|/g > 2A, at most one solution has |a| ≤ A. G is then injective on the box, and the box values
  have density ≈ U/(2·max|G|) ≈ E/m. For c61, g = 1.
- A generator of smooth s that is blind to s mod Y₀ hits the box with that probability. Its cost per
  rational-smooth pair found is m/E = L[2/3, 1/δ].
- Generating every smooth s ≤ Em costs Ψ(Em, B) ≫ E². Smooth-number enumeration and generation are
  KNOWN: Bernstein's thesis chapter "Enumerating and counting smooth integers"; Bach 1988; Kalai 2003;
  Bach–Sorenson. None of them imposes the residue window.

**Chosen factor.** If s₁ is fixed, the admissible s₂ are in bijection with the box points of Λ_{s₁}. So a
generator that *enumerates those lattice points and tests the cofactors* is an LQ algorithm.

**Scope.** This is an equality of candidate sets, not an algorithmic equivalence. A generator that produced
admissible smooth cofactors without visiting the lattice points would fall outside LQ. None is known.

**Measured (E3).** Toy instance: N of 64 bits, d = 3, m = 60311. Hit rate 0.016675 against 2A/m = 0.016979
(ratio 0.982). The c61 rate is 2A/|Y₀| = 1.09·10⁻⁸.

### 3.4 Proposition W (balanced, height-bounded k-trees) [model-PROVED size accounting for that class only]

**The inequality.**
- Wagner's k-tree (CRYPTO 2002) needs bottom lists of size ≥ (m/E)^{1/(1+lg k)} to hit a window of E out
  of m residues.
- A product of k factors that is ≤ Em allows at most (Em)^{1/k} values per list.
- Feasibility therefore needs (1 + lg k)/k ≥ log(m/E)/log(Em) = 1 − O(log E/log m).
- At L[1/3] against L[2/3] only k ≤ 2 survives (`wagner_max_k`: k = 8 at 512 bits, 4 at 4096, 2 from about
  2^16).
- k = 2 is meet-in-the-middle, whose cost per solution is ≥ (m/(KE))^{1/2} = L[2/3]. The small-list
  extended k-tree (Minder–Sinclair, SODA 2009 / J. Cryptology 2012) does not escape this.

**What it does not cover.** It is not a lower bound for arbitrary product search: unequal height strata,
rejection sampling, ratios and other merges need their own analysis.

**A further obstacle.** Wagner's merging is additive, while the window |s mod m| ≤ E is not compatible with
multiplicative merging.

**The relation-lattice reformulation (e2).** It is KNOWN: Schnorr 1991; Adleman 1995. Its known algorithms
are exponential in the dimension π(B); see Ritter–Rössner 1997, Ducas 2021, and Mulder–Sterner–van Woerden
2025, who find smooth twins through SVP.

### 3.5 Proposition BC (Boneh CRT batching) [HEURISTIC; random-integer model; within Boneh's proven decoding range]

**The decoding theorem.**
- Boneh, JCSS 64 (2002), STOC 2000: Thm 3.1 for intervals, Thm 4.1 for monic f of degree D.
- One polynomial-time call finds every t with |t| < H and gcd(f(t), P) ≥ S, provided
  log H ≲ (log S)²/(D log P).
- This is the Coppersmith / Howgrave-Graham divisor-modulus bound (Cohn–Heninger, Thm 1.1).
- Use it per row b: D = 1, H = 2E, f made monic by multiplying by Y₁⁻¹ mod P, which needs gcd(Y₁, P) = 1.
  This gives log P ≲ (log S)²/log E.

**The hit rate (Rankin).** For a random integer n,

  Pr[gcd(n,P) ≥ S] ≤ S^{−1/2}·E[gcd(n,P)^{1/2}] ≤ S^{−1/2}·exp(Σ_{p|P} p^{−1/2}).

- Given Σ_{p|P} log p ≤ Λ, the sum is largest when P consists of the smallest primes, where it is
  O(√Λ/log Λ) = O(log S/(√(log E)·log Λ)) = o(log S).
- So a call returns ≤ 2E·S^{−1/2+o(1)} candidates.
- This needs no prime cutoff. The earlier "y_P ≈ log N" argument assumed a primorial-type P and has been
  replaced.
- Stratify by the actual gcd S′ ≥ S. A candidate's cofactor gain is (C·u log u)^{log S′/log B} = S′^{o(1)}.
- So the yield per unit cost is ≤ P(y)·(C·u log u)^{v}·min(1, E·S^{−1/2+o(1)}). This is maximised at
  S ≤ E^{2+o(1)}, i.e. v = O(1), giving L^{o(1)}.

**Scope.**
- Boneh's threshold is a *sufficient* guarantee, not an impossibility result. Boneh §4 already observes that
  the bounds "are insufficient for improving the basic quadratic sieve" (see also Guruswami's thesis §12.3.2).
- Bivariate, whole-box Coppersmith/Boneh constructions are *not* covered by this row-by-row accounting. They
  would need their own lattice, extraction and yield analysis. We do not exclude them.

**Meet-in-the-middle over CRT residues (c3).**
- For each b and each q₁ > 2E there is at most one admissible a.
- Finding the q₁ that have one in row b is sieving that row. Batch remainder trees cost (E + |Q|)·polylog
  per row, so the total is E·(E + B) ≈ E² at e = β.

### 3.6 Proposition FT (factory transfer to one N) [model-PROVED; items 3 and 4 KNOWN or folklore]

**1. The table is chargeable.** For one N the rational table is charged to that N. max(2e, per-N terms) is
exactly the box/MNFS program, so the result is 1.9019.
- 2e = 2.0069 is the table size *at the 1.6386 factory point*, not a universal minimum.
- Amortising it down to 1.6386 needs about L[0.368] targets.
- Amortising it to just below 1.9019 needs about L[0.105] targets (codex).
- One N has one target.

**2. kN is MNFS.** The base-m polynomials f_k of kN for many k, with one m, all have the root m mod N.
- Relations combine mod N as in MNFS: the rational side is shared and the algebraic columns are pooled.
- The coefficient growth from V = L[0.125] is of lower order.
- The program is MNFS, so 1.9019. No prior statement found by us of the kN phrasing; the construction is
  Coppersmith's.

**3. A shifted m with the same field is a coordinate change.** m → m+t, f(X) → f(X − t), a → a + tb
preserves both norms (folklore: Murphy 1999; Kleinjung 2006).
- A shifted m with its *own* base-m polynomial is a new field. Its rational table is a shear of the old
  one.
- Adjacent shifts overlap in 2E − b of row b, about 75% of the table, so J shifts cost about E²(2 + J/2)
  instead of 2JE². This corrects the earlier draft's false "no sharing" claim.
- MNFS with a single m shares 100% across its V fields, so shifted-m families are dominated by MNFS.

**4. A universal m (e.g. 2^t) is the factory itself.** It is non-uniform, with a once-for-all L[2.0069]
precomputation (Coppersmith 1993; Bernstein–Lange 2013 §5.1: "conventional NFS chooses m as a function of
N, while this attack does not").
- Practical precedent: the Mersenne factory (Kleinjung–Bos–Lenstra, ASIACRYPT 2014) shares *algebraic*
  sieving across 17 SNFS numbers.

**Aside (codex; standard, Buhler–Lenstra–Pomerance 1993).** A prime dividing both norms of a primitive pair
divides Res(f, g) = ±N (for c61, f(m)·Y₁⁴ ≡ 0 mod N). A two-sided special-q with q₁ = q₂ is therefore
useless or reveals a factor.

### 3.7 (f) Adaptive special-q collection for a smaller dependency core (r21)

**The idea.** Collect relations whose primes lie in a head of size L[η], plus at most k tail primes in
(L[η], L[β]]. Choose special-q among the tail primes so that cycles close.

**k = o(u) [model-PROVED under H + H_LQ + H_sel + the generic-dependency heuristic].** The latter means that
about L[η]^{1−o(1)} relations are needed (Pomerance 1996; Croot–Granville–Pemantle–Tetali 2012).
- The relation probability is L[−r(η) − a(η)]·u^{O(k)}, and u^{O(k)} = L^{o(1)} while k log u = o(λ).
- A closing relation for tail prime q is a q-lattice point. Its yield per unit cost is the same, by LQ.
- So the program is MNFS with β := η, with LA ≥ 2η: 1.9019.

**k = Θ(u) [OPEN].** Here relations are just L[β]-smooth and the tail is a hypergraph.
- The r21 toy bound (c ≤ h independent tail cycles with independent uniform head labels contain a
  dependency with probability ≤ 2^{c−h}) assumes independence. Adaptive selection could violate it.
- Bender–Canfield's random-core model predicts a linear core, but it is not a proof for adaptively selected
  relations.
- So cheap elimination of polylog-weight hyperedges is **not excluded**. The earlier draft's exclusion is
  withdrawn.
- Related literature: large-prime cycles (Lenstra–Manasse 1990; Dodson–Lenstra 1995; Leyland et al. 2002);
  Cavallar's filtering; special-q descent in DL individual logs.

### 3.8 What is *not* proved

- No lower bound excludes an output-sensitive rational enumerator. Rational listing is smooth-number
  enumeration in ⋃_{|b|≤E}[bm − E, bm + E], i.e. E short intervals along an arithmetic progression.
- We found no hardness reduction, and the literature critic found none either (cf. Pomerance's ICM 1994
  survey on smooth numbers in algorithms).
- Propositions LQ, GEN, W and BC are barriers for explicitly modelled access. H_sel assumes away exactly the
  capability an improvement would need.
- Pomerance's 1982 L[1/2] barrier was broken by NFS; the same could happen here.

## 4. Prototypes against las (Task 4; pre-registered)

**Instance.**
- N = `balanced_semiprime(200, 1501)` (C61).
- Polynomial: CADO polyselect + ropt, degree 4, skew 9594 (`experiments/raw/smoothenum/c61.poly`).
- Parameters: CADO params.c60.
- Box: |a| < 2^21, 1 ≤ b ≤ 2^9, U = 2^31.

**Memory.** Every run streamed rows or lattices. Peak RSS was ≤ 57 MB for our runs and 37 MB for las.
Memory pressure was normal (≥ 47% free) before each run; cadobench was running CADO jobs concurrently,
which explains the wall-time variance below.

**Registered outcomes.**

| Pred. | Statement | Result | Verdict |
|---|---|---|---|
| P1 | box rational yield / Σ P_s(\|G\|) ∈ [0.80, 1.25] | 36774 / 40727 = 0.903 | **pass** |
| P1′ | sieve recall ≥ 0.95 | 6904/7117 = 0.970 | **pass** |
| P2 | lattice point counts / (U/M) ∈ [0.95, 1.05] for t ≤ 28 | 1.0000, 1.0000, 1.0000, 1.0002, 1.0011, **0.925** (t = 28: 340 vs 367.5 expected over 60 lattices) | **FAIL at t = 28** |
| P2 | same ratio ∈ [0.90, 1.10] for t ≥ 31 | 0.997, 1.073, 1.081, 0.925 | **pass** |
| P3 | observed / Σ_{tested pts} P_s(\|G\|/M) ∈ [0.75, 1.33] per bucket (registered raw estimator) | 0.900, 0.870, 0.916, 0.937, 0.955, 0.818 (t = 8..28); 0.938, 0.901, 0.982, 0.990 (t = 31..40, all points tested) | **pass** |
| P4(i) | argmax_t Y_t ∈ {28, 31, 34} | t = 28 (Y = 0.2425) | **pass** |
| P4(ii) | max Y_t/Y_box ≤ 1.5 × predicted max | 13.8 ≤ 1.5 × 18.1 | **pass** |
| P4(iii) | Y₃₇ < Y₃₄, Y₄₀ < Y₃₇; slope within ±0.3 of the prediction | 0.0329 > 0.0057 > 0.0009; slope −0.866 vs −0.904 | **pass** |
| P5 | bottom-up hit rate / (2A/m) ∈ [0.90, 1.10] | 0.982 | **pass** |
| P6a | ≥ 99.5% of las relations pass our exact tests | registered run 4436/4436; corrected run 4473/4473 | **pass** |
| P6b | las ≥ 20× faster per relation | registered: 69× wall (1.65 s vs 113.4 s), 33× CPU. Corrected K = 4473: 24× wall (2.55 s vs 61.0 s), 26× CPU (2.27 s vs 59.4 s) | **pass** (margin small in the corrected run) |
| P6c | our two-sided yield per unit cost ∈ (1, 10⁴) × the box two-sided yield per point | 4.88·10⁻⁴ / 3.34·10⁻⁵ = 14.6 against the sieved baseline; 13.8 against the exact baseline (74 relations, not 70; `e1_exact.json`) | **pass** |

**About the P2 failure.**
- The t = 28 shortfall corresponds to a Poisson z of −1.43. A Poisson law for these counts is not
  established, and this explanation is post hoc.
- The registered failure stands.
- It concerns the geometry check only. The yield checks P3 and P4 in the same bucket pass.

**E2: yield per unit cost of smooth-modulus lattices (rational side).**
- Definitions: Y_t is rational relations per (points + setups). Y_box = 0.01754 per box point (sieved).
- The prediction uses the coprime fraction 6/π². The E1 rows (b odd) have a higher coprime fraction (0.818),
  so the measured gain is deflated relative to the prediction.

| t = log₂M | 8 | 12 | 16 | 20 | 24 | 28 | 31 | 34 | 37 | 40 |
|---|---|---|---|---|---|---|---|---|---|---|
| Y_t/Y_box (measured) | 2.0 | 3.3 | 5.0 | 7.8 | 11.6 | **13.8** | 8.1 | 1.9 | 0.33 | 0.05 |
| Ŷ_t/Ŷ_box (predicted from geometry + P_s only) | 2.9 | 4.8 | 7.3 | 10.6 | 15.3 | 18.1 | 10.8 | 2.5 | 0.41 | 0.06 |

**E4: best mechanism against las at equal output.**
- *Correction to the las configuration.* The registered las command needed `-dup-qmin/-dup-qmax`. The first
  run used `61961,61961`, which suppressed duplicates against side-0 jobs that were never run; codex caught
  this.
- The corrected run (`-dup-qmin 0,61961 -dup-qmax 0,63961`) returns 4473 relations: 37 more, none lost.
  We reran our mechanism to K = 4473. Both runs are reported.

| | las (C, 2 threads), corrected | ours: special-(q₁,q₂), Python, 1 thread (K = 4473) |
|---|---|---|
| relations | 4473 (192 special-q, sqside 1, q ∈ [61961, 63961]) | 4473 distinct |
| examined units / relation | 9.65·10⁷ sieved points / 4473 = 21 570 | (6.05·10⁶ points + 3.14·10⁶ setups) / 4473 = 2054 |
| wall / relation | 0.57 ms | 13.6 ms |
| CPU per examined unit | ≈ 23.5 ns | ≈ 6.5 µs (3.69·10⁶ exact two-sided tests) |

**Reading E4.**
- The two-sided double special-q examines about 10.5× fewer units per relation. Each unit is a lattice
  reduction or an exact two-sided test, about 275× dearer than one las sieve point.
- A C implementation would need under about 250 ns per unit to tie las.
- *The gain is the finite u^{O(1)} gain of §3.2* (two forced divisors on 56- and 100-bit norms). One
  instance is consistent with that explanation, but does not demonstrate it.
- *The comparison is descriptive.*
  - The regions differ: las relations reach |a| = 2.5·10⁷ and b ≤ 1958, and only 985 of the 4473 lie in
    our box.
  - The cost units are heterogeneous.
  - Equal raw relation counts do not guarantee equal dependency yield.

**Disclosed deviations.**
- las needed `-dup-qmin/-dup-qmax` with `-dup`. The registered run used the wrong side value; it was
  corrected and rerun, as above. `-sqside 1` was made explicit (it was the default).
- In the K = 4436 run, the `tests` counter was overwritten by a key collision; this was fixed before the
  K = 4473 run.
- Exact-baseline and corrected-las runs were added after review. They are not registered predictions.
- Wall times vary by up to about 2× between runs because of concurrent CADO jobs.

## 5. Relation to r20 §5

- r20 #5 (search entropy against norm height) needs many search coordinates to share one height cost.
- Output-sensitive enumeration is the enumeration-side analogue: the search cost would track the number of
  hits instead of the entropy.
- Inside the modelled access (divisibility lattices, products, CRT within Boneh's range), both hit the same
  wall. Known factors of total size ≤ area reduce u only by O(1), and the gain is ρ(u − O(1))/ρ(u) =
  L^{o(1)}.

## 6. GPT-6.1-sol blind round (`codex/smoothenum_1.out`)

It was run before our analysis.

**Agreements.**
- Every constant in Table 1.
- The ambient m/E waste of residue-blind generation.
- The coverage argument for special-q.
- The failure of Boneh's decoding bound.
- Factory transfer = MNFS.
- The black-box floor.

**Adopted from codex.**
- The reduced closed form for c(θ) (tested).
- The common-prime remark.
- The coordinate-change remark.

**Ours.**
- The LQ formalisation and the g* tables.
- The k-tree size inequality.
- The Rankin hit-rate accounting.
- Mechanism (f).
- All the measurements.

## 7. Reviews

### 7.1 Opus cost-model critic: ACCEPT WITH FIXES

The critic independently recomputed T1–T7 and the g* table (all reproduced) and found no missed improvement.

| # | Sev. | Finding | Fix applied |
|---|---|---|---|
| 1 | MAJOR | H_sel's additive O(1) destroys the M > U branch | H_sel restated: E[new points \| 𝓕] ≤ C·U/M for M > U, no additive term; not needed for M ≤ U |
| 2 | MINOR | the ρ-ratio step was justified only for v = 1 | log-derivative integral, uniform in v |
| 3 | MINOR | special-q above B was not in the model | column / large-prime accounting line added |
| 4 | MINOR | the MNFS case was not spelled out | step 5 of the proof |
| 5 | MINOR | the g* table was rational-only and over-read; the distinct-yield paragraph was unsupported | algebraic row added; "illustrative, with multiplicity"; paragraph withdrawn; verdict reconciled |
| 6 | NIT | decay rate | now (log x)^{1/3}/x^{1/3} |
| 7 | MINOR | T8 listing all outputs is exponential | prefix-output hypothesis |
| 8 | MINOR | "2.0069 unavoidable" was overstated | FT(1) rewritten; L[0.105] noted; relabelled |
| 9 | MINOR | BC used y_P ≈ log N and ρ at polylog smoothness | replaced by a Rankin bound with no cutoff |
| 10 | MINOR | §3.7 status | labels and assumptions added |
| 11–12 | NIT | gcd(Y₀,Y₁), Minder–Sinclair, T5 1.8850, δ ↑ 6β_a | applied |

### 7.2 Opus literature critic

Fixes applied:
- Schnorr 1991 / Adleman 1995 for (e2) and T8, plus Ritter–Rössner, Ducas, and Mulder–Sterner–van Woerden.
- Boneh: Thm 3.1/4.1 citations; P built from prime powers; the monic requirement; Boneh §4 already notes the
  insufficiency.
- Bernstein–Lange 2013 §5.1 for T3 and FT(4).
- The Mersenne factory.
- Composite and hybrid special-q: Boudot §4.1 and fn. 8, Kleinjung 2006, Aoki 2008.
- Pollard 1993 and Franke–Kleinjung; de Bruijn and Hildebrand–Tenenbaum.
- The early-abort L^{o(1)} example.
- Wagner and Minder–Sinclair.
- Cycle literature.
- The shift and common-prime remarks marked folklore.
- Bernstein / Bach / Kalai for smooth generation; Pomerance's ICM 1994 survey.
- The T6 = ½·SNFS identity.

Answer to the critic's Q4: no published L-power gain was found for any mechanism here.

### 7.3 GPT-6.1-sol adversarial review (`codex/smoothenum_2.out`): "major revision", all items applied

1. **LQ selection gap.** Same as 7.1 #1. Additions: the history 𝓕, new-points-only counting, conditional
   joint randomness, and the statement that H_sel excludes cheap identification *by assumption*.
2. **T8.** Prefix/delay, joint smoothness and fruitfulness hypotheses added.
3. **Exclusions withdrawn.** Unbounded-tail elimination (§3.7, now OPEN, with the extension to k = o(u))
   and whole-box multivariate decoding (§3.5, not excluded).
4. **las.** Duplicate-suppression misconfiguration found; rerun (4473); E4 rerun.
5. **FT.** The false no-sharing claim corrected to 75% adjacent overlap, dominated by MNFS; the L[0.105]
   threshold added.
6. **GEN and W narrowed; BC hit rate repaired.**
   - GEN: the reduction is restricted to generators that enumerate lattice points.
   - W: restricted to balanced, height-bounded k-trees.
   - BC: the prime cutoff is replaced by the Rankin bound.
7. **Numbers and literature.**
   - T5 corrected to 1.8850.
   - The g* grid now includes v_box exactly.
   - Decay rate corrected.
   - Exact baseline (74, not 70): gain 13.8×.
   - P2 failure kept without excuse.
   - P3 uses the registered raw estimator.
   - Guruswami thesis and Cohn–Heninger cited.

After the fixes, the verdict "no speedup; 1.9019" stands. All three reviewers agree with it.

## 8. Reproduction

```bash
.venv/bin/python -m src.smoothenum.target          # Table 1 (~4.5 min)
.venv/bin/python -m src.smoothenum.barrier         # (g* tables via lq_gain_table(), k-tree table)
.venv/bin/python -m pytest -q tests/test_smoothenum.py
# instance (CADO build of cadobench, read-only):
B=~/src/cado-nfs/build/kmini; cd experiments/raw/smoothenum
$B/polyselect/polyselect -N <N> -degree 4 -P 420 -admax 10000 -incr 60 -nq 64 -keep 10 -t 2 > c61_raw.txt
$B/polyselect/polyselect_ropt -inputpolys c61_sizeopt.poly -ropteffort 0.1 -I 10 -Bf 111342 -Bg 78682 -t 2
$B/sieve/makefb -poly c61.poly -lim 111342 -maxbits 10 -out c61.roots1.gz -t 2
$B/sieve/las -poly c61.poly -fb1 c61.roots1.gz -lim0 78682 -lim1 111342 -lpb0 18 -lpb1 19 -mfb0 17 -mfb1 38 \
   -ncurves0 2 -ncurves1 2 -I 10 -q0 61961 -q1 63961 -sqside 1 -t 2 -dup -dup-qmin 0,61961 \
   -dup-qmax 0,63961 -out c61_las_fixed.rels          # registered run used 61961,61961 -> c61_las.rels
cd ../../..
for c in calib e1 e1exact e2 e3; do /usr/bin/time -l .venv/bin/python -m src.smoothenum.run_proto $c; done
.venv/bin/python -m src.smoothenum.run_proto geo
.venv/bin/python -m src.smoothenum.run_proto e4check experiments/raw/smoothenum/c61_las_fixed.rels
/usr/bin/time -l .venv/bin/python -m src.smoothenum.run_proto e4ours 4473 3600   # (4436 = registered run)
```
