# sicheck: validating Prop SI (smooth numbers in single intervals near kN ⇒ L[1/3, 1.2114])

Agent `sicheck`, phase 20, 2026-10-21.
- Pre-registration: `notes/sicheck_prereg.md`, commit `2bba07e`, made before every E1/E3/E4 run.
- Code: `src/sicheck/core.py` (operational: interval sieve, constructors, relations, F₂ algebra, square roots,
  `factor_si`, the Blum bit β), `model.py` (exact sympy re-derivation), `diag.py` (harness: nulls, Blum
  classification), `run_experiment.py`, `score.py`, `posthoc.py` (disclosed post-hoc). Tests: `tests/test_sicheck.py`.
- Raw data: `experiments/raw/sicheck_E{1,3,4}.jsonl`, scores `sicheck_E{1,3,4}_score.txt`, post-hoc
  `sicheck_E1_posthoc.{jsonl,txt}`.
- GPT-6.1-sol: blind round `codex/sicheck_1.{md,out}` (§6.1); adversarial round `codex/sicheck_2.{md,out}` (§6.2).
- Reviews: one read-only Opus proof/cost critic and one read-only Opus literature critic (§7).
- Siblings: `siconstruct` and `sihard` had only pre-registrations when this note was written; nothing here depends on
  them.

**Labels.** KNOWN (published; cited) · PROVED (proof here; elementary) · model-PROVED (exact consequence of a named
model or heuristic) · HEURISTIC (named assumption) · experimentally supported / refuted for this model (toy scale) ·
OPEN · UNVERIFIED (citation not checked). "No prior statement found by us" is the strongest novelty claim made.

**Notation.** x = log N, F_α = x^α(log x)^{1−α}, L[α, c] = exp((c + o(1))F_α). y = L[1/3, κ], h = L[2/3, ν],
I_k = [kN, kN + h], u_k ∈ I_k y-smooth (the oracle's output), v_k = u_k − kN. N = pq with p, q > y.

---

## 0. Verdict

**Plain language.** The arithmetic of Prop SI is right, and so is the constant: if a constructor really delivers
y-smooth numbers in these short intervals at negligible cost, and its outputs behave like the standard heuristics
assume, then N factors in L[1/3, (16/9)^{1/3} + ε] ≈ L[1/3, 1.2114]. Nothing in the accounting is given away. The free
lunch is the oracle itself: it replaces a search that costs L[2/3] by every implementation we analysed (this is a
benchmark, not a lower bound) with one that costs L^{o(1)}. But the two heuristics are **not automatic** for a
deterministic constructor. What we establish is the *conditional* cost optimisation and exact degeneracy results;
constructor availability, fresh-relation yield and fruitfulness remain **separate, explicit hypotheses** (§5):

1. **H_S fails at toy scale in its naive form, for a structural reason. A refined form agrees with a post-hoc
   benchmark, and if it holds, its rate loss is asymptotically negligible.**
   - gcd(u_k, v_k) = gcd(u_k, k) exactly (Lemma I, PROVED).
   - So (i) most smooth residues are *inherited*: they are scaled copies (u, v) = g·(u′, v′) of a relation already
     obtained at k/g, with the identical relation vector. (ii) A *fresh* relation's residue is coprime to every prime
     of u, which is a correlation against smoothness.
   - Measured on 54 000 intervals (40–80 bits, four constructors):
     - 49–71% of smooth residues are duplicates.
     - On primitive outputs, the smooth rate is only 0.20–0.86 of the size-matched null (pooled 0.66–0.68).
     - A post-hoc coprimality-corrected benchmark (conditioning on rad(u/g)) agrees within 3% pooled, on all
       first-seen reduced pairs (O/E = 1.008–1.026). Conditioning also on k/g gives 1.14–1.17, so the agreement
       depends on the conditioning choice. It is a benchmark, not a calibrated validation (codex #4).
     - Net effect: fresh relations arrive at **10–69% of the naive H_S rate**, with the deficit shrinking from 40 to 80 bits
       (monotonically for `next`; cells also differ in y, h and K, so the trend is confounded).
   - **Refined H_S [HEURISTIC]:** v/g behaves like a random integer of its size coprime to u/g (and to k/g). It is
     supported at 40–80 bits only. *If* it holds, the rate loss relative to naive H_S is ≥ 1 − o(1) from coprimality
     (Lemma CP; KNOWN in cumulative form, Tenenbaum 1993 / Xuan 1995) and a constant from duplicates (Lemma D, under a
     joint divisibility model; or L^{o(1)} model-free by restricting to k without prime factors ≤ y). So the
     refined heuristic is consistent with the L-scale rate ρ(log h/log y)·L^{−o(1)} that Prop SI needs. It is not
     derived from the naive one.
2. **R is a real requirement: a deterministic constructor that sees only kN and N can make every dependency
   trivial.**
   - For Blum N, fruitfulness of a Prop SI dependency is a linear functional computable from N alone (Prop B, PROVED):
     β_i = ((Ω_T(u_i) − Ω_T(v_i))/2) mod 2, with T = {ℓ : J_N(ℓ) = −1}.
   - Hence (Cor. B, Blum N) a J_N = +1 factor base (the shape of a QS factor base), or a constructor that avoids β = 1
     outputs, gives only trivial square roots.
   - Experiments (E3, E3b):
     - The four natural constructors gave fruitful dependency spaces in 144/144 runs on 36 moduli (49.6% of sampled
       dependencies fruitful, which is algebra once Φ ≢ 0).
     - The J+ base on Blum N gave 13/13 degenerate runs, against 16/16 fruitful on the same N with the full base.
     - beta0 gave 18/18 degenerate runs, all at the maximal parity rank π(y) − 2.
   - No hidden-root construction or Dixon-style guarantee is supplied for the I_k interface (§3.4; an argument, not an
     impossibility theorem). Randomising k or the choice in I_k
     reduces R to a static condition on the pool (Prop R′), at L^{o(1)} cost if that condition holds. For Blum N, R is
     checkable and can be imposed as one extra row.
3. **End-to-end toy.** The operational algorithm with a brute-force (sieving) constructor factored 30/30 N of 40–80
   bits, all through dependencies, with 0.51 of dependencies fruitful. Its fresh-relation yield was 14–60% of the
   naive model, and its cost is L[2/3], as it must be without the oracle.
4. **Claim (iii)** (single intervals cannot beat exponent 1/3) is correct in the typical-interval / random-residue
   model and only there (codex). **The domain of k matters:** k ≍ N^θ gives (16(1+θ)/9)^{1/3}, so θ = 1 is SNFS's
   1.526 (codex, re-derived exactly).

**Scorecard of registered predictions** (details §2.3, §3.5, §4).

| # | Prediction | Result |
|---|---|---|
| P1 | O/E (size-matched null) ∈ [0.85, 1.15] per cell, [0.92, 1.08] pooled | **FAILED**: 10 of 24 cell×constructor cells outside, from 0.76 to 2.78. The pooled values 0.98–1.01 are a cancellation, not a confirmation (§2.3) |
| P2a | odd-k ratio < 1 < even-k ratio | held (next 0.887 / 1.135; rand 0.903 / 1.096) |
| P2b | parity model within ±15% per stratum | held (0.952 / 1.063; 0.971 / 1.024) |
| P3a | 1 ≤ P(empty)/e^{−μ̂} ≤ 2.5 | **FAILED** in B10 (μ̂ ≈ 8–9: ratio 27–33); held in the other 10 cell-seeds (1.32–2.30) |
| P3b | KS(next), KS(last) ≤ 0.10 vs truncated exponential | **FAILED** in 3 of 12 cell-seeds (B10 next 0.126 / 0.157; A s0 last 0.108) |
| P3c | KS(rand v/h), KS(cofactor v/u₀) ≤ 0.05 | **FAILED**: cofactor 0.82–0.93 in every cell (the constructor coincides with `next`, §2.3); rand failed once (A s0 0.071) |
| P4 | duplicates: next > rand, next ∈ [0.05, 0.5]; non-primitive ≥ 0.2 | **FAILED**: next 0.708 > 0.5 (next > rand 0.560 and non-primitive ≥ 0.71 held) |
| R1 | honest constructors fruitful | held: 144/144 Φ ≢ 0; 0.496 of random dependencies fruitful |
| R2 | J+ base: Blum degenerate, non-Blum not | held but thin in E3 (1 + 1 runs); **held** in E3b: 13/13, 16/16 controls, 12/12 |
| R3 | beta0 degenerate | held: 18/18 |
| R4 | Prop B exact on Blum | held: 0 mismatches in 7825 dependencies |
| R5 | prefix fruitless rate ≈ 2^{−d} | held: 0.576, 0.285, 0.146, 0.062 (the factor-2 tolerance is weak) |
| R6 | parity rank fb − 2 | held: 142/144 |
| F1 | factor found in ≥ 95% of runs, all STRUCTURAL | held: 30/30, 0 ACCIDENTAL |
| F2 | smooth outputs / null ∈ [0.85, 1.15] | held (0.966), but the count includes duplicates; the fresh yield is 0.14–0.60 of the null (§4) |

---

## 1. Re-derivation (Task 1)

### 1.1 The relation, and why v must be smooth
- u − v = kN, so u ≡ v (mod N). With e = exp(u) − exp(v) ∈ Z^{π(y)}, ∏ℓ^{e_ℓ} ≡ 1. So e lies in the relation lattice
  Λ_N = {e : ∏ℓ^{e_ℓ} ≡ 1 (mod N)}.
- A dependency c (Σc_ie_i ≡ 0 mod 2) gives Z = ∏ℓ^{(Σc_ie_i)/2} with Z² ≡ 1. It splits N iff Z ∉ {±1}. The code
  computes Z = W/V, with W² = ∏u_iv_i and V = ∏v_i.
- **v must be smooth** because e needs v's full factorisation over the factor base. A non-smooth v adds one column per
  large prime. Large-prime variants need their own bounded analysis (codex #15); the program analysed here is
  strictly smooth. Smoothness tests cost L^{o(1)}: ECM (Lenstra, Ann. Math. 126 (1987)) at L_y[1/2] = L[1/6], or
  Bernstein's batch smooth-parts algorithm (2004 preprint, "How to find smooth parts of integers").

### 1.2 Cost [model-PROVED; `model.si_program`, exact sympy]
CEP in L-units: an integer of log-size s·F_b is L[a, k]-smooth with probability L[b − a, −(b − a)s/k].
- u ≈ kN has log-size x = F₁ (if log k = o(x)), so it is y-smooth with probability L[2/3, −2/(3κ)]. Mean count per
  interval L[2/3, ν − 2/(3κ)]; existence threshold ν₀ = 2/(3κ).
- v ≤ h has log-size νF_{2/3}, so it is smooth with probability L[1/3, −ν/(3κ)] (under H_S).
- Relations needed π(y) = L[1/3, κ]. Oracle calls L[1/3, κ + ν/(3κ) + τ], where τ is the oracle's per-call cost.
  LA L[1/3, 2κ]. Smoothness tests of v (ECM or batch) cost L^{o(1)}.
- ∂(trials)/∂ν = 1/(3κ) > 0, so ν = ν₀. The trial term is then κ + 2/(9κ²).
  - It meets 2κ at κ*³ = 2/9.
  - Its minimiser is κ_min³ = 4/9 > 2/9, so it is decreasing at the crossing, and the max is minimised there.
  - **c = 2κ* = (16/9)^{1/3} = 1.211414**, κ* = 0.605707, ν* = 6^{2/3}/3 = 1.100642.
  - All exact; test `test_model_exact_constants`.
- This is an infimum. With fixed slack ν > ν₀ the constant is max(κ + ν/(3κ), 2κ) > 1.2114 (suffmap §1.4, codex).

### 1.3 Claim (iii) [model-PROVED in the typical-interval / random-residue model]
At smoothness level y = L[α, κ]:
- Existence needs h ≥ L[1 − α, (1 − α)/κ].
- The residue is then smooth with probability L[1 − 2α, −(1 − α)(1 − 2α)/κ²] (`single_interval_exponent_program`).
- The cost exponent is max(α, 1 − 2α) ≥ 1/3, with equality only at α = 1/3.
- Shorter intervals (ν < ν₀) waste L[2/3, ν₀ − ν] calls on empty intervals and gain only L[1/3, (ν₀ − ν)/(3κ)] on
  the residue (`below_threshold_never_helps`).
- **Scope** (codex sicheck_1 §1). This is not a lower bound for arbitrary single-interval constructors. A constructor
  that biases v towards small or smooth values leaves the model, and in the limit it is a doubly-smooth lister
  (suffmap E3).

### 1.4 What kind of relation it is
- **Not a Dixon/QS relation.**
  - In Dixon, x² ≡ s with x freely sampled. The hidden choice of x among the roots of x² makes a dependency fruitful
    with probability ≥ 1/2.
  - Here both sides are factor-base products. The relation is a vector of Λ_N (an S-unit pair, abcsunit Prop SU), and
    Z is a deterministic function of the relation set.
  - QS also has no *hidden* root (its x is deterministic), so R is a heuristic there too. But a QS relation is not an
    S-unit relation: its root x + m is not a factor-base product, so J_N(Z) = ∏J_N(x_i + m) carries an extra coin per
    relation (critic, finding 1). In QS one side is a square, which removes one smoothness requirement; Prop SI pays
    for that side with the oracle.
- **Exact classification.** It is the degree-1 two-form sieve of nfsfruit Cor. J (J2), with L₁(k, v) = kN + v and
  L₂(k, v) = v. Then d₁L₁ ≡ c₁L₂ with c₁ = d₁ = 1, so no |S|-parity constraint is imposed. Fruitfulness is therefore
  the Lemma FR / Lemma H question (abcsunit, nfsfruit), not Dixon's.
- Classical analogues with one side constructed: Schroeppel's linear sieve for factoring (analysed in Pomerance,
  Math. Centre Tracts 154 (1982)) and the cubic sieve (Coppersmith–Odlyzko–Schroeppel, Algorithmica 1 (1986) 1–15,
  a discrete-log paper), with residues N^{1/2} and N^{1/3}. Their relations over the extended base {H + c} are also
  S-unit relations, so Prop B / Cor. B apply to them.
- **Closest prior form (literature critic): Schnorr's fac-relations** (Schnorr, EUROCRYPT '91, LNCS 547, 281–293;
  Adleman 1995; Ducas 2021 experiments): u and |u − kN| both smooth, with k smooth too and u built by lattice CVP.
  A fac-relation is a Prop SI relation with smooth k, and Prop B / Cor. B apply verbatim. Schnorr's u and v have
  disjoint support, so every relation is primitive by construction (the "gcd(u, k) = 1" device of Lemma D).

### 1.5 Consistency with known heuristics; what the oracle gives away
- **Size principle** (suffmap Theorem SS(b); Pomerance 2008). A random part of log-size s·F_{2/3} gives c = 2√(s/3)
  at base L[1/3, √(s/3)].
  - Prop SI's only random part is v, with s = ν* = 1.1006.
  - √(ν*/3) = κ* exactly (`size_principle_consistency`).
  - So 1.2114 is what the principle predicts once the N-sized side is free.
- **Random-size budgets** s = 3c²/4 (exact, `comparison_table`, except the factory row, which is numeric). The
  "consistency" check is the same equation κ³ = 2/9 seen from Theorem SS(b), not independent evidence.

  | Method | c | s |
  |---|---|---|
  | Prop SI | 1.2114 | 1.1006 |
  | SNFS | 1.5263 | 1.7472 |
  | factory, per N | 1.6386 | 2.0137 |
  | GNFS | 1.9230 | 2.7735 |

- **Why it beats SNFS and the factory without special form.**
  - SNFS uses special N to make *both* norms small.
  - The factory lists the L[2/3]-sized rational side "for free" after precomputation, but still pays the random
    algebraic norm.
  - Prop SI makes the **N-sized** side free and pays only a residue of size L[2/3]. The coincidence with Kirchner's
    special-n constant is the same arithmetic (suffmap §1.4).
- **What is given away, exactly:** a y-smooth integer of size N in a window of length ≈ 1/δ = L[2/3, 2/(3κ)].
  - Every implementation we know of pays ≈ 1/δ per call (a benchmark, not a lower bound): sieving the window, scanning with ECM, or the cofactor scan
    (codex; our toy constructors sieve the whole window).
  - With such an implementation the algorithm is L[2/3], which is our toy (§4).
  - The oracle is therefore a collapse of one search problem from exponent 2/3 to o(1). The rest of the accounting
    (relations, LA, residue penalty) is honest.
  - Whether such a constructor can exist is siconstruct/sihard's question.
- **Domain of k** (codex sicheck_1 §2; `si_constant_with_k_exponent`, exact). If k ≍ N^θ, the threshold becomes
  2(1+θ)/(3κ) and the constant c³ = 16(1 + θ)/9. θ = 1 gives (32/9)^{1/3} = 1.5263. So log k = o(log N) is part of
  the hypothesis. The program uses k ≤ L[1/3, c], for which log k = O(F_{1/3}) = o(x).

---

## 2. Heuristic H_S for deterministic constructors (Task 2)

### 2.1 Structure (Lemma I PROVED; Lemma CP KNOWN in cumulative form; Lemma D HEURISTIC)

**Lemma I (inheritance).** Let gcd(N, ∏_{ℓ≤y} ℓ) = 1, let u ∈ I_k be y-smooth and v = u − kN > 0.
- (a) gcd(u, v) = gcd(u, k) =: g.
- (b) (u′, v′) := (u/g, v/g) satisfies u′ − v′ = (k/g)N, u′ ∈ I_{k/g}, v′ ≤ h/g, and exp(u) − exp(v) =
  exp(u′) − exp(v′). It is the same relation vector.
- (c) With k′ = k/g: gcd(u′, k′) = 1 because g = gcd(u, k); hence gcd(u′, v′) = gcd(u′, k′N) = 1 and
  gcd(v′, k′) = gcd(u′ − k′N, k′) = 1. So v′ is coprime to every prime of u′ and of k′. In particular, for every prime ℓ ≤ y that divides u but not k, ℓ ∤ v.
- *Remark (critic).* The PH2 null p_cop conditions on rad(u) only. That is the right null *averaged* over consecutive
  or uniform k: the extra coprimality with k cancels against the enrichment P(ℓ | v) = 1/(ℓ − 1) at primes ℓ ∤ uk.
  It is not a per-k null, and for prime k there is no such cancellation.

*Proof.* gcd(u, v) = gcd(u, u − v) = gcd(u, kN) = gcd(u, k), since gcd(u, N) = 1. Divide by g for (b). (c) holds by
construction, and for ℓ | u, ℓ ∤ kN we have v ≡ −kN ≢ 0 (mod ℓ). ∎ (Codex states (a) and (c) independently.)
*Status:* standard in substance. It is the NFS rule of using only coprime (a, b) on the forms (kN + v, v). The
coprime-pair density 6/π² appears in NFS yield analyses (Buhler–Lenstra–Pomerance, LNM 1554, 1993), and the
coprimality correction is the analogue of Murphy's α (Murphy, PhD thesis, ANU 1999). The k-parity effect is the
k-analogue of the Knuth–Schroeppel multiplier function (TAOCP vol. 2 §4.5.4; Morrison–Brillhart, Math. Comp. 29
(1975)). We found no prior discussion of inheritance across multiples kN as such (literature critic).

**Consequences.**
- Distinct relations correspond to distinct *primitive pairs* (u′, v′). A constructor that, at k = g·k′, returns
  g·u′_{k′} (this happens whenever g·v′ ≤ h and no other choice wins) re-emits an old relation.
- `next` re-emits it whenever no smooth number lies in [kN, kN + g·v′). `rand` re-emits it with probability
  1/(number of smooth integers in I_k) when g·v′ ≤ h.
- Every fresh relation has a residue coprime to rad(u′). For y-smooth u′ that is a strong condition at toy scale,
  because smooth numbers are rich in small primes.

**Local uniformity (the reason the naive H_S looks plausible).** Suppose k runs over consecutive integers (or uniform
random k), M ≪ h, and the constructor's choice of u mod M is independent of kN mod M. Then v = u − kN mod M is
exactly uniform. So the *marginal* divisibility of v by every small modulus is that of a random integer. The parity
strata of P2 show the local correlation directly: for odd k, v is even iff u is odd. Since π₀ = P(u even) ≈ 0.6–0.8 at
toy scale, odd k depresses P(v smooth) and even k raises it, and they cancel on average. Restricting k to odd numbers
or primes > y (codex's duplicate-free domain) removes the cancellation. This is an L^{o(1)} effect at the L-scale
(Lemma CP), but a constant factor at toy scale.

**Lemma CP (coprimality penalty) [KNOWN in cumulative form; local-window version HEURISTIC].**
- *Model.* For a set S of primes ≤ y and z in the range of v, the density of y-smooth integers near z among integers
  coprime to S is the unrestricted smooth density times ∏_{p∈S} f_p, with f_p = (1 − p^{−α})/(1 − p^{−1}) and
  α = α(z, y) the Hildebrand–Tenenbaum saddle point.
  - For cumulative counts this is a theorem covering our range (literature critic, checked in Granville 2008 §4.2,
    MSRI 44): Xuan, Acta Arith. 69 (1995) 303–327, gives Ψ_q(x, y) ~ ∏_{p|q, p≤y}(1 − p^{−α})Ψ(x, y) when q has
    ≤ y^{1/2} prime factors ≤ y, for (log x)^{1+o(1)} < y < x^{o(1)}; Tenenbaum, Phil. Trans. R. Soc. A 345 (1993)
    377–384, gives Ψ_q ~ (φ(q)/q)Ψ when ω_y(q) ≤ y^{o(1/log u)}. Here |S| ≤ 2x/log 2 ≪ y^{1/2}, so the claim below
    is in substance Tenenbaum 1993. See also de la Bretèche–Tenenbaum, Ramanujan J. 9 (2005) 139–202, Thm 2.1
    (codex). These are asymptotic; at 40–80 bits the exact-count null is the right tool.
  - The step from counts to *constructor-selected* residues (refined H_S) remains HEURISTIC.
  - At toy scale it is the null that fits E1 (§2.3).
- *Claim.* At the Prop SI scale (z ≤ h = L[2/3, ν], y = L[1/3, κ]), uniformly over S = the prime factors of any
  y-smooth u ≤ N² (or of u·k, k ≤ N), ∏_{p∈S} f_p ≥ 1 − o(1), for u_z → ∞ (for z below about y·e^{−0.58}, α ≥ 1 and f_p ≥ 1). The proof is for z ≍ h; since α(z, y) ≥ α(h, y) for z ≤ h and f_p
  increases with α, the product at any z ≤ h is at least its value at h, which is the lower bound needed (codex #2).
- *Proof.*
  1. Write 1 − α = δ. Then δ log y = ξ(u_z) + O(1), with ξ(u) = log(u log u) + o(1) (Hildebrand–Tenenbaum, Trans. AMS 296
     (1986) 265–290, Lemma 2; Hildebrand, J. Number Theory 22 (1986), Lemma 4; KNOWN), and
     u_z = log z/log y ≤ (ν/κ)x^{1/3}(log x)^{−1/3}. So y^δ ≤ e^{O(1)}·u_z log u_z = x^{1/3+o(1)} and
     δ = O((log x)^{1/3}x^{−1/3}).
  2. −log f_p = (1 + o(1))(p^δ − 1)/(p − 1) uniformly, since p^{−α} = p^{δ−1} ≤ y^δ/p.
  3. *Small primes*, p ≤ x²: p^δ − 1 ≤ 2δ·2 log x, so their total contribution is
     ≤ 4δ log x·Σ_{p≤x²} 1/(p−1) = O(δ log x·log log x) = o(1).
  4. *Large primes*, p > x² in S: there are ≤ log(N²)/log(x²) ≤ x of them, each contributing ≤ 2y^δ/p ≤
     x^{1/3+o(1)}/x². The total is ≤ x^{−2/3+o(1)} = o(1).
  5. Hence −log ∏f_p = o(1). ∎
- *Toy scale.* α ≈ 0.4–0.6 and f₂f₃f₅f₇ ≈ 0.1–0.3, which matches the measured primitive deficits 0.2–0.4 at 40–48
  bits. A rough illustration at N ≈ 2^1024 with κ*, ν*: y ≈ 2^27, u_z ≈ 8.6, α ≈ 0.85, and the product over the typical
  small primes of u is ≈ 0.6–0.66 (critic: α ≈ 0.82), a constant. The error bound δ·log x·log log x of step 3
  exceeds 1 at 1024 bits, so the o(1) is not small at any practical size (consistent with C7).

**Lemma D (duplicates cost a constant) [HEURISTIC: needs a joint law; codex #1].** Let k be uniform on [1, K].
*Model assumption (joint, not marginal):* for the chosen u, the events {ℓ | u} are independent of k and of one another,
with P(ℓ | u) = ℓ^{−α_u}(1 + o(1)), α_u = α(N, y) → 1. Marginals alone do not give the product (a coupling of [ℓ | u]
with [ℓ | k] keeps the marginals and kills primitivity), and no constructor's selection law is shown to satisfy it.
So the "0.6ρ" bound is **not proved**; the model-free route below is. Then
- P(output primitive) = P(gcd(u, k) = 1) ≥ ∏_ℓ (1 − ℓ^{−1−α_u})(1 − o(1)). For α_u < 1 this is below 1/ζ(2) ≈ 0.61,
  and it tends to 1/ζ(2) only in the limit, so "≥ 0.6" holds eventually;
- every primitive output with smooth v is a fresh relation (its primitive pair determines (k, u), by Lemma I).
- So the fresh rate is ≥ 0.6·P(v smooth | primitive)·(1 − o(1)), and with Lemma CP this is ≥ 0.6·ρ(log h/log y)·
  L^{−o(1)}. ∎
- *Model-free alternative (critic).* If k has no prime factor ≤ y, gcd(u, k) = 1 automatically. Such k have density
  ≍ 1/log y = L^{−o(1)}, so restricting to them costs only L^{o(1)} without any divisibility model.
- Alternatively, prime k > y makes every output primitive (codex). The fresh rate is then P(v smooth | coprime to
  rad u) directly, with the k-parity bias of P2 uncancelled (L^{o(1)}).

**Size distribution.** The constructors differ only in the law of v's size:
- `next`: a forward gap, ≈ exponential with rate δ, truncated at h;
- `last`: h minus a backward gap;
- `rand`: ≈ uniform on [0, h];
- `cofactor`: **coincides with `next`** in 97–99.7% of intervals. The least u₀ > v dividing a smooth u is just above
  v (median v/u₀ = 0.99), because smooth numbers have dense divisors. So the picture "v uniform below u₀" is false.

A smaller v only helps, so at the L-scale every one of these laws gives P(v smooth) ≥ ρ(log h/log y)·L^{−o(1)}. With
slack ν > ν₀, `next` even does better than ρ(log h/log y) (codex). At the infimum the slack vanishes.

**What H_S cannot cover.** A constructor that sees the candidate v's (it sees kN, so it can) may prefer non-smooth or
smooth residues. Preferring smooth ones makes it a doubly lister (suffmap E3); preferring non-smooth ones defeats
Prop SI. So H_S must be stated as a **requirement**: the choice among smooth u ∈ I_k is not anti-correlated with the
smoothness of (u − kN)/g, beyond the forced coprimality of Lemma I.

### 2.2 Experiment E1 (registered; §2 of the pre-registration)
Six cells, 40–80 bits (A: 40/y = 100; B, B10: 48/200; C: 56/600; D: 64/1500; E: 80/6000), two N each, K = 2000–6000
consecutive k, 54 000 intervals in all. Each interval is sieved exactly, and all four constructors are applied to the
same smooth list. Null p_null(v): exact smooth density in [v ± max(200, v/20)]. Peak RSS 1.5 GB (cell E, from the
enumeration of smooth integers ≤ 1.06h for the null), one job at a time.

### 2.3 E1 results

**Registered predictions (`experiments/raw/sicheck_E1_score.txt`).**
- **P1 FAILED.** O/E against the size-matched null, over all outputs: 2.78 (A next) down to 0.76 (B next). 10 of 24
  cells are outside ±15% (next 3, last 2, rand 2, cofactor 3; corrected from "14" after codex sicheck_2 #7). The pooled values (next 1.011, last 0.980, rand 0.999, cofactor 1.013) sit inside
  [0.92, 1.08] only because cell A's large excess cancels the deficits of B–D. They do not confirm H_S.
- **P2 held.**
  - Odd k: naive O/E 0.887 (next) and 0.903 (rand). Even k: 1.135 and 1.096.
  - The parity model, with π₀ = P(u even) measured, gives 0.952 / 1.063 (next) and 0.971 / 1.024 (rand).
  - So the k-dependence is the predicted local correlation (v even iff u odd, for odd k). The remaining even-k
    excess (1.063) is plausibly inheritance raising P(u even | k even), since π₀ was pooled over both parities
    (critic).
- **P3 partly failed.**
  - P(empty)/e^{−μ̂} is 1.3–2.3 in 10 of 12 cell-seeds; smooth integers are over-dispersed by inheritance. In B10
    (μ̂ ≈ 8–9) e^{−μ̂} ≈ 10^{−4} and the ratio is 27–33: a mis-specified prediction (the mean count drifts four-fold
    across the k-range, column "mean count by k-quarter").
  - KS for next/last ≤ 0.10 except B10 next (0.13, 0.16) and A s0 last (0.108); rand ≤ 0.05 except A s0 (0.071).
  - **cofactor 0.82–0.93 everywhere:** it coincides with `next` (§2.1). The pre-registered "v uniform below u₀"
    model was wrong.
- **P4 FAILED.** Duplicate fraction among smooth-residue outputs: next 0.708, cofactor 0.709, rand 0.560, last 0.492
  (pooled). Non-primitive ≥ 0.71 for all constructors. next > rand held; next ≤ 0.5 failed. I under-estimated
  inheritance by a factor of about 2.

**Disclosed post-hoc analysis** (PH1–PH3, `src/sicheck/posthoc.py`; motivated by the P1/P4 failures, not registered).
Primitive outputs are those with gcd(u, v) = 1. The coprime null p_cop(v; u) is the exact smooth density among
integers coprime to rad(u) in the same window.

| cell (bits) | constructor | O/E all (P1) | dup. | primitive O/E, naive (PH1) | primitive O/E, coprime null (PH2) | fresh / naive E (PH3) |
|---|---|---|---|---|---|---|
| A (40) | next | 2.778 | 0.96 | 0.331 | 1.394 | 0.100 |
| B (48) | next | 0.761 | 0.78 | 0.414 | 1.057 | 0.165 |
| B10 (48) | next | 0.838 | 0.77 | 0.469 | 1.195 | 0.189 |
| C (56) | next | 0.862 | 0.72 | 0.586 | 0.960 | 0.244 |
| D (64) | next | 0.904 | 0.63 | 0.744 | 0.998 | 0.337 |
| E (80) | next | 1.030 | 0.61 | 0.848 | 1.010 | 0.407 |
| A | rand | 2.236 | 0.93 | 0.257 | 1.225 | 0.150 |
| B | rand | 0.788 | 0.68 | 0.390 | 1.054 | 0.250 |
| B10 | rand | 0.982 | 0.55 | 0.416 | 1.215 | 0.438 |
| C | rand | 0.907 | 0.63 | 0.601 | 1.005 | 0.337 |
| D | rand | 0.917 | 0.49 | 0.710 | 0.981 | 0.471 |
| E | rand | 1.045 | 0.41 | 0.863 | 1.042 | 0.621 |
| A | last | 1.664 | 0.89 | 0.199 | 1.029 | 0.183 |
| E | last | 1.060 | 0.35 | 0.804 | 0.983 | 0.690 |
| A–E | cofactor | ≈ next | ≈ next | ≈ next | ≈ next | ≈ next |

Pooled over all cells, the primitive outputs give O/E = 0.657–0.680 against the naive null and **1.001–1.016**
against the coprime null (z between 0.0 and +0.8; 1650–2317 successes per constructor).

**PH2b (after codex sicheck_2 #4–5; `posthoc.py ph2b`, `experiments/raw/sicheck_E1_posthoc_ph2b.txt`).** PH2 dropped
the fresh relations that first appear through a non-primitive output (45% of `last`'s, 35% of `rand`'s). PH2b uses
every output whose reduced pair (u/g, v/g) has not occurred before among *all* earlier outputs, smooth or not. That
selection does not depend on the outcome, and it contains every fresh relation. Pooled O/E against the coprime null
over rad(u/g): next 1.008, last 1.012, rand 1.026, cofactor 1.014. Per cell: 0.96–1.39 (the excesses are in the
40–48-bit cells, with 42–147 successes). As a sensitivity check, the null over rad(u/g·k/g), which also imposes the
proved coprimality of v′ with k′, gives 1.14–1.17 pooled (1.09–1.89 per cell). The rad(u)-only null is the right one
*on average over consecutive k* (the k-coprimality cancels against enrichment at the other primes; proof critic), but
neither null is a generative model of the constructor's selection, and the Bernoulli z-scores are not calibrated for
these dependent, selected observations.

**Reading [experimentally supported at 40–80 bits, as a post-hoc benchmark].**
- The naive H_S ("v is a random integer of size ≤ h") is wrong at toy scale, in both directions. Inherited residues
  are smooth more often (they are g·v′ with v′ smooth); fresh residues are smooth less often (coprime to rad u). The
  P1 cells mix the two.
- The refined H_S ("v/g is a random integer of its size coprime to u/g") agrees with the pooled data to ≤ 3% under
  the rad(u/g) benchmark. The per-cell excesses (A, B10: 1.2–1.4) rest on 15–147 successes. A pre-registered joint
  selection model tested on fresh instances would be needed to call it validated.
- **The quantity that the cost model needs, fresh relations per call, is 10–69% of the naive H_S rate**, and the
  deficit shrinks with size (next, monotone: 0.10, 0.17, 0.24, 0.34, 0.41 at 40, 48, 56, 64, 80 bits; rand is not
  monotone, B10 0.44 > C 0.34; y, h and K also change across cells, so this is a confounded trend). This is
  what Lemma CP predicts (and Lemma D under its joint model): both factors are constants or 1 − o(1) at the L-scale,
  but they are large at toy scale. No asymptotic inference is drawn from the trend. `last` and `rand` lose less than `next`, because a larger v makes inheritance less likely.

---

## 3. Heuristic R: fruitfulness and rank (Task 3)

### 3.1 The criterion [PROVED; KNOWN in substance]
- For relation vectors e_i ∈ Λ_N and the dependency space D = {c ∈ F₂^R : Σc_ie_i ≡ 0 (mod 2)}, the map
  Φ : D → μ₂(Z/N)/{±1}, c ↦ [∏ℓ^{(Σc_ie_i)/2}] is a homomorphism (abcsunit Lemma FR(1); nfsfruit Lemma H; codex
  sicheck_1 §3). For N = pq the target is F₂.
- So either every dependency is trivial (Φ ≡ 0), or exactly half of D is fruitful.
- Lemma FR(2)–(4): im Φ ≠ 0 iff (Λ₀ ∩ 2Z^s)/2 reaches a mixed square root of unity, where Λ₀ is the span of the
  relations. Generating Λ_N up to finite odd index makes im Φ the 2-torsion of the subgroup generated by the factor
  base, so it suffices *provided that subgroup contains a mixed root of unity* (it need not: 2 has odd order 33 mod
  161; codex #8). Maximal parity rank alone is not enough (Lemma FR(5), codex's N = 77 example).
- **Rank.** The parity rows of uv are orthogonal to the indicator vectors of {ℓ : (ℓ/p) = −1} and {ℓ : (ℓ/q) = −1},
  since (uv/p) = (u/p)² = 1. So r ≤ π(y) − rank_{F₂}{c_p, c_q}: this is π(y) − 2 only when the two vectors are
  distinct and nonzero (counterexample N = 161, y = 3: rank 1 > 0; codex #8); on a J+ base they coincide and the
  bound is only fb − 1 (so R6's "all runs" clause was mis-specified for qsfb; it held).
  Literal full rank is impossible, and that is harmless.

### 3.2 Prop B: for Blum N, fruitfulness is computable from N [PROVED]
- *Setting.* N = pq, p ≡ q ≡ 3 (mod 4), p, q > y. Relations (u_i, v_i) with u_i ≡ v_i (mod N), both y-smooth.
  T = {ℓ ≤ y : J_N(ℓ) = −1}, Ω_T(n) = Σ_{ℓ∈T} v_ℓ(n), a_i = Ω_T(u_i), b_i = Ω_T(v_i).
- *Claim.* (1) a_i ≡ b_i (mod 2). (2) A dependency S is fruitful iff Σ_{i∈S} β_i ≡ 1 (mod 2), where
  **β_i = ((a_i − b_i)/2) mod 2**.
- *Proof.*
  1. J_N(u_i) = J_N(v_i) because u_i ≡ v_i (mod N), and J_N(n) = (−1)^{Ω_T(n)} by multiplicativity.
  2. Write Z ≡ (s_p, s_q) ∈ {±1}² by CRT. Since (−1/p) = (−1/q) = −1, J_N(Z) = s_ps_q. So Z ∉ {±1} iff
     J_N(Z) = −1.
  3. Z = W/V with W = ∏ℓ^{E_ℓ/2}, E_ℓ = Σ_S (v_ℓ(u_i) + v_ℓ(v_i)), and V = ∏_S v_i. Then J_N(W) =
     (−1)^{Σ_S (a_i + b_i)/2}, each term an integer by (1), and J_N(V) = (−1)^{Σ_S b_i}.
  4. J_N(Z) = J_N(W)J_N(V) = (−1)^{Σ_S [(a_i + b_i)/2 + b_i]}. Finally (a + b)/2 + b ≡ (a − b)/2 (mod 2). ∎
- *Status.* This is nfsfruit's Corollary J (J2) specialised to the forms (kN + v, v), and it was derived blind by codex
  (sicheck_1 §3, "t_i = j^Te_i/2"). We claim no novelty. The ingredients are classical: a square root of 1
  with J_N = −1 is nontrivial for Blum N (Rabin 1979; Blum 1982), and the proper-subgroup mechanism is in Seysen,
  Math. Comp. 48 (1987) 757–780, p. 777. The literature critic found no precise prior statement of β or of Cor. B
  for S-unit relations. Exhaustively checked: test
  `test_beta_formula_and_lemma_H_blum`, and E3 R4.
- *General N* (nfsfruit Prop G specialised; PROVED there). Let ψ = χ_pχ_q with χ_p, χ_q characters of 2-power order
  2^{v₂(p−1)}, 2^{v₂(q−1)}, c = max of the two valuations, ψ(ℓ) = ζ^{τ_ℓ}, ζ = e^{2πi/2^c}. Then
  2^c | Σ_ℓ τ_ℓ(v_ℓ(u_i) − v_ℓ(v_i)), and Φ(S) = Σ_S β_i with β_i = Σ_ℓ τ_ℓ(v_ℓ(u_i) − v_ℓ(v_i))/2^c mod 2. For
  non-Blum N, χ_p has order ≥ 4 and β is hidden.

### 3.3 Degenerate constructors [PROVED]
**Cor. B.** For Blum N:
- (a) if every factor-base prime has J_N(ℓ) = +1, every dependency is trivial;
- (b) if every relation used has β_i = 0, every dependency is trivial.
- (a) is the shape of a QS factor base ({ℓ : (N/ℓ) = 1} = {ℓ : J_N(ℓ) = 1} for odd ℓ, since N ≡ 1 mod 4). It does
  not bite in QS because a QS relation is not an S-unit relation: its root x + m lies outside the factor base and
  contributes J_N(x + m) per relation. (The proof critic checked this: QS on 6 Blum N ≈ 2^35 with a J+ base containing
  2 had fruitful dependencies in 6/6; our E3b J+ base also contains 2 and is fully degenerate. An earlier draft wrongly
  attributed the difference to ℓ = 2.)
- (b) is realised by the constructor `beta0`. It sees only kN and N and returns the least smooth u ∈ I_k that does not
  produce a relation with β = 1.
  - The degeneracy is PROVED. Its yield is HEURISTIC: if β behaves like a fair coin, beta0 keeps a relation with
    probability ≈ 1/2, so its relation yield is ρ·L^{o(1)} (codex's marked-candidate model). Our beta0 also sieves
    the whole interval, so it is not a fast constructor; a wrapper that rejects a fast oracle's β = 1 relations is the
    fast version, with the same conditional yield (codex #10).
  - **So the L-scale smoothness rate does not imply R.**
  - Totality caveat (codex): an interval whose smooth candidates all give β = 1 relations has no admissible output.
    beta0 then outputs nothing, which costs only a fraction of intervals.
- Hence "the constructor sees only kN" does not protect against degeneracy. The constructor needs a hidden-character
  neutrality requirement (§5, C5).

### 3.4 Randomisations and what they buy
- **Hidden roots (Dixon/Rabin) are not available through the I_k interface.** What is PROVED is tautological: I_k
  relations contain no element outside the factor base. The cost statements below are HEURISTIC estimates for
  specific constructions, not lower bounds (critic finding 7). Codex agrees.
  - Every Prop SI relation is φ(e) = 1, and no freely sampled root enters it.
  - Choosing k as a function of a random square does not insert its root.
  - A single hidden-root relation r² ≡ a/b with a, b smooth would make every dependency containing it fruitful with
    probability 1/2 (lemmaM Theorem DD). Producing one generically means a, b ≈ √N jointly smooth (rational
    reconstruction gives ab ≈ N), at cost L[2/3, 2/(3κ)] [HEURISTIC]: the same exponent as one brute-force oracle
    call. A QS-type relation (x² − kN smooth, root outside the base but not hidden) costs L[2/3, 1/(3κ)] and would
    heuristically restore a fruitfulness coin. Every such route found is L[2/3].
  - It would need a *coset* query ("smooth u ≡ t·smooth"), which is exactly the capability Prop DLT asks for and the
    I_k oracle lacks (suffmap §1.4, §2 Prop DLT).
- **Prop R′: random k or random choice in I_k** reduce R to a static property of a pool [model-PROVED; nfsfruit
  Prop R transcribed].
  - Let X be the set of relations obtainable from all k ≤ K′, with square-class vectors x and Bockstein bits β.
  - Let δ_X = min_w Pr[β ≠ ⟨w, x⟩] under the *sampling law*, i.e. the pushforward of the distribution of k (and of
    the in-interval choice). Inherited relations appear with multiplicity, so this is not uniform on X.
  - M relations drawn i.i.d. from X are degenerate with probability ≤ 2^{π(y)}(1 − δ_X)^M. So
    M ≥ (π(y) ln 2 + ln(1/ε))/δ_X suffices: an L^{o(1)} enlargement iff δ_X ≥ L^{−o(1)}.
  - The bound survives sampling without replacement (nfsfruit), which is what deduplication does.
  - X, its weights after acceptance (only smooth residues kept) and the deduplication rule must be specified, and
    δ_X computed under exactly that law (codex #12). With i.i.d. sampling with replacement no lower bound on |X| is
    needed; without replacement, M ≤ |X| distinct relations are required.
  - The default C1 domain (consecutive k, deterministic constructor) is a fixed relation set: Prop R′ does not apply
    to it, and R is then a yes/no property of N.
  - Against the adversarial beta0 rule this changes nothing (δ_X = 0 on its pool). Randomisation helps against
    *accidental* degeneracy only.
- **Random multipliers m** (query I_{km}): the same as restricting k. Multiplying both sides by a known smooth integer
  leaves e unchanged (codex).
- **Blum N: check and impose.** By Prop B, Φ ≢ 0 is decidable from the relations in time poly(π(y)). If it holds,
  solving Σc_ix_i = 0, Σc_iβ_i = 1 gives a dependency that is fruitful with certainty (one extra dense row; nfsfruit
  Cor. J). This *detects* fruitfulness but cannot create it.

### 3.5 Experiment E3 (registered) and E3b (Addendum A1)
**E3.** 36 N (40, 48, 56 bits; 6 Blum + 6 non-Blum each); constructors next, last, rand, cofactor, plus J+ base
(when J_N(2) = +1) and beta0 (Blum). Each collects |base| + 20 distinct primitive relations, or stops at k = 15 000.
`experiments/raw/sicheck_E3_score.txt`.

| # | Prediction | Result |
|---|---|---|
| R1 | honest: Φ ≢ 0 in ≥ 95%; random-dependency fruitful rate ∈ [0.42, 0.58] | **held**: 144/144 runs Φ ≢ 0; 4567/9216 = 0.496 (per constructor 0.48–0.51) |
| R2 | J+ base: Blum Φ ≡ 0 always; non-Blum Φ ≢ 0 in ≥ 50% (≥ 5 deps) | held but **thin**: one Blum run with dependencies (56 bits, 22 deps, 0/64 fruitful); one non-Blum run (40 bits, p ≡ q ≡ 1 mod 4: 30/64 fruitful). The J+ base yields almost no relations at 40–56 bits; hence E3b |
| R3 | beta0 (Blum): Φ ≡ 0 always | **held**: 18/18 runs (17 reached the relation target), all dependencies trivial |
| R4 | Prop B: 0 mismatches on Blum; mismatches on non-Blum | **held**: 0 / 7825 Blum dependencies; non-Blum mismatch rate 3130/6280 = 0.498 (the Blum formula is uninformative there, as expected) |
| R5 | P(Φ ≡ 0 on first d deps) ∈ [0.5, 2]·2^{−d} | **held**: 0.576, 0.285, 0.146, 0.062 for d = 1..4 (n = 144) |
| R6 | r ≤ fb − 2 always; = fb − 2 in ≥ 90% | **held**: 142/144 honest runs at fb − 2 exactly; beta0 18/18 at fb − 2 (degenerate *at maximal parity rank*) |

**E3b** (Addendum A1; 32 and 36 bits; 8 Blum and 8 non-Blum N per size, all with J_N(2) = +1; constructor `next`;
`experiments/raw/sicheck_E3b_score.txt`).

| # | Prediction | Result |
|---|---|---|
| R2b | Blum, J+ base: Φ ≡ 0 in every run with a dependency (≥ 10 runs) | **held**: 13/13 runs (7–12 dependencies each), 0 fruitful of 416 sampled |
| R2c | Blum, full base, same N: Φ ≢ 0 in ≥ 90% | **held**: 16/16; 218/512 = 0.43 fruitful |
| R2d | non-Blum, J+ base: Φ ≢ 0 in ≥ 50% (≥ 5 deps) | **held**: 12/12 (p, q mod 4 = (1,1): 2, (1,3): 6, (3,1): 4); 209/384 = 0.54 |

So Cor. B(a) is a Blum-specific obstruction, exactly as proved. The same N with the full base factor every time.

**Reading.**
- For the four natural deterministic constructors, Φ ≢ 0 in every run. That is the evidential content. Once Φ ≢ 0,
  "half of D is fruitful" is algebra, and the sampled fraction of nonzero dependencies is 2^{d−1}/(2^d − 1); so 0.496
  checks the sampler, not neutrality (codex #14). The prefix statistic (R5) is the one test of the fair-functional
  model, and its tolerance is weak. The 144 runs share 36 moduli and overlapping candidate pools; they are not 144
  independent tests. This matches nfsfruit's toy NFS/QS results. **Experimentally supported at 40–56 bits; heuristic
  R is plausible for oblivious constructors.**
- E3's raw `k_last` is the shared loop counter, not each constructor's stopping time (codex #17); the code now also
  records `k_done`. Only the "mean k_last" display in the score file is affected.
- Degeneracy is easy to produce on purpose without knowing p (beta0), and it appears at the maximal attainable parity
  rank fb − 2. So rank statistics cannot certify R; only Φ (or, for Blum N, β) can.

## 4. End-to-end toy (Task 4)

`core.factor_si` is operational. It sees N, y, h and a seed. For k = 1, 2, … it sieves I_k exactly (the brute-force
constructor), applies `rand` or `next`, keeps distinct primitive relations until it has π(y) + 20, then runs F₂
elimination, square roots and gcd. The constructor receives only X = kN. The harness checks factors and computes the
null; it never passes p, q to the algorithm.

**Registered results** (`experiments/raw/sicheck_E4_score.txt`; 30 runs, 40–80 bits, seeds 0–2).
- **F1 held.** A nontrivial factor was found in 30/30 runs, every one STRUCTURAL (from a dependency), with 0
  ACCIDENTAL.
  - Random dependencies were fruitful in 980/1920 = 0.51.
  - Every run had 22–29 dependencies, as expected from π(y) + 20 relations at rank π(y) − 2.
- **F2 held, but only in the same cancelling sense as P1.** Pooled, smooth-residue outputs / Σp_null = 16 087 /
  16 648.5 = 0.966. This count includes duplicates.
- **F3 (descriptive), pooled per size:**

| bits | constructor | intervals | duplicates among smooth outputs | **fresh relations / Σp_null** | fruitful deps |
|---|---|---|---|---|---|
| 40 | next / rand | 10 138 / 6 673 | 0.85 / 0.72 | **0.14 / 0.26** | 0.49 / 0.46 |
| 48 | next / rand | 10 483 / 7 181 | 0.81 / 0.68 | 0.19 / 0.28 | 0.49 / 0.54 |
| 56 | next / rand | 7 951 / 6 011 | 0.70 / 0.55 | 0.25 / 0.38 | 0.49 / 0.54 |
| 64 | next / rand | 6 442 / 4 920 | 0.63 / 0.48 | 0.35 / 0.52 | 0.54 / 0.58 |
| 72 | rand | 8 100 | 0.52 | 0.49 | 0.48 |
| 80 | rand | 9 878 | 0.42 | 0.60 | 0.49 |

**Reading.**
- The algorithm works end to end, and its dependencies behave as R predicts.
- Its relation yield is set by *fresh* relations. These arrive at 14–60% of the naive H_S rate, rising with size. This
  is the same deficit as E1 PH3, now inside the real algorithm.
- The cost is that of the brute-force constructor: h sieved positions per interval, i.e. L[2/3]. Toy successes are
  correctness checks only (AGENT_RULES 8). Dickman ρ(log h/log y) under-predicts the size-matched null by 5–50% at
  these sizes (expected at small u; cadobench calibrated ρ only at real NFS sizes).

## 5. Verdict and the precise constructor requirements (Task 5)

**Verdict.** Prop SI's L[1/3, (16/9)^{1/3} + ε] is a **correct heuristic implication** (model-PROVED in the stated
model), *provided* the constructor satisfies C1–C7 below. C1–C3 are the content of the hypothesis. C4–C6 are what
H_S and R amount to for a deterministic constructor. For natural oblivious constructors, C4 and C5 held at toy scale
(C4 in the refined form). For C4 the refinement costs only 1 − o(1) times a constant at the L-scale (Lemma CP; duplicates by
Lemma D under its joint model, or model-free on k without small prime factors).
C5 stays a heuristic, reduced only by Prop R′.
Prop B / Cor. B show that C5 is a real condition and cannot be dropped.

Let T* = L[1/3, κ + ν/(3κ)] be the number of calls the program needs, with κ ↑ κ* = (2/9)^{1/3} and
ν ↓ 2/(3κ) slowly.

- **C1 (domain of k).** The queried multipliers form a set 𝒦 with |𝒦| ≥ T*·L^{o(1)} and log k = o(log N) for
  k ∈ 𝒦 (k ≍ N^θ gives (16(1+θ)/9)^{1/3}).
  - 𝒦 = an interval of consecutive integers, or uniform random k in [1, K]: inherited duplicates occur (Lemma I) but
    cost a constant factor under Lemma D's joint model (L^{o(1)} model-free if k has no prime factor ≤ y), and the local k-correlations of v cancel (P2).
  - 𝒦 = primes in [K, 2K] with K > y (codex): no duplicates, but the k-parity bias is uncancelled (an L^{o(1)}
    factor).
  - Either way, every probabilistic statement in C4–C5 must hold on that exact 𝒦.
- **C2 (availability and cost).** On at least an L^{−o(1)}-fraction of k ∈ 𝒦, the constructor returns a y-smooth
  u ∈ [kN, kN + h] with its factorisation, at total cost T*·L[1/3, τ] with τ = o(1). (Each τ > 0 adds τ to the trial
  term. MNFS is beaten iff τ < 0.757 *after re-optimising κ* (the optimum moves to κ³ = 4/9 for τ > 0.381); at fixed
  κ = κ* the threshold is ≈ 0.69 (codex #16; suffmap).) A promise on intervals of length h with ν − 2/(3κ) → 0⁺ is needed
  for the limit constant; with fixed slack, the constant is max(κ + ν/(3κ), 2κ).
- **C3 (no leakage).** The constructor sees only kN, h, y (and possibly N and randomness): no p, q, φ(N). Our toy
  constructors satisfy this; `test_operational_code_has_no_leakage` checks the code.
- **C4 (refined H_S: fresh-relation rate).** Among the first T outputs, the number of **distinct primitive** pairs
  (u/g, v/g) with v smooth is ≥ T·ρ(log h/log y)·L^{−o(1)}.
  - *Sufficient:* refined H_S [HEURISTIC]: v/g is distributed like a random integer of its size coprime to u/g
    (and k/g). Then the rate follows [via Lemma CP, and Lemma D or the k-coprime domain]: the coprimality factor is ≥ 1 − o(1) and
    duplicates cost a constant. A choice rule that ignores v's factorisation does **not** by itself make v random
    (for `next`, v is a gap between smooth numbers); that is exactly the heuristic content.
  - *Measured (post-hoc benchmark):* the refined null agrees at 40–80 bits (pooled O/E 1.008–1.026 on all first-seen
    reduced pairs; 1.14–1.17 if also conditioned on k). The fresh rate is 10–69% of the naive
    H_S rate, improving with size.
  - *Excluded:* any rule that prefers residues with large prime factors. Rules that prefer smooth residues are
    allowed: they only help, and in the limit they are doubly listers.
- **C5 (R: hidden-character neutrality).** On the collected relation set, the Bockstein functional β (2-adic
  characters χ_pχ_q, §3.2) is not in the F₂-span of the parity rows, i.e. Φ ≢ 0. With M ≥ (π(y) ln 2 + ln(1/ε))/δ
  relations, failure probability ≤ ε.
  - *Sufficient:* outputs drawn (by random k or by random choice in I_k) from a pool X with δ_X ≥ L^{−o(1)} (Prop R′).
  - *Necessary in practice:* for Blum N, the factor base must not be restricted to J_N(ℓ) = +1 primes (Cor. B(a);
    E3b R2d shows this base is fine for non-Blum N), and the selection must not discriminate *against* β = 1
    (Cor. B(b)). Using Jacobi data is not in itself harmful: imposing Σβ = 1 uses it and helps.
  - For Blum N, C5 is **checkable** from N (Prop B) and can be imposed with one extra row.
  - No hidden-root (Dixon) construction is supplied for the I_k interface (§3.4), so C5 stays a heuristic, as for QS
    and NFS.
  - *Measured:* 144/144 natural runs fruitful; 0.496 of dependencies.
- **C6 (rank / count).** At least rank + 1 distinct relations. The parity rank is ≤ π(y) − rank{c_p, c_q}
  (= π(y) − 2 generically), with equality in 142/144 runs. Maximal parity rank is not evidence for C5 (beta0 is degenerate at
  rank π(y) − 2).
- **C7 (scale).** The constant is an infimum: approaching it needs ν − 2/(3κ) → 0 slowly enough that the promise C2
  still holds, with uniform control of the o(1) terms (codex). At toy and practical sizes the o(1) terms are large:
  constant factors 2–10 in relation yield at 40–80 bits.
- **C8 (joint, quantified control; codex #13).** C2, C4 and C5 must hold jointly, under one sampling law, uniformly
  over the intended N and the queried multipliers. The domain 𝒦 must be efficiently enumerable or samplable, and the
  availability fraction, the cost of failed calls and the failure response must be specified. Global smooth-number
  density does not give availability on the thin union of intervals near kN; the finite-domain promise is an explicit
  assumption.

**Answer to "is something being given away for free?"** Only the oracle's search, which is the hypothesis. No
relation, residue or linear-algebra cost is missing. The two heuristics are not free: H_S needs the refined form C4,
which fails naively at toy scale and whose rate loss is asymptotically negligible (the refined form itself is a
heuristic); R needs C5, which a p-oblivious constructor can violate.

## 6. GPT-6.1-sol rounds

### 6.1 Blind round (`codex/sicheck_1.out`, run before our analysis)
**Independent agreements.** The cost program and the infimum (16/9)^{1/3} (κ* = 0.6057, ν* = 1.1006); claim (iii) in
the typical-interval model only; v must be smooth for this factor-base representation; the relation is not a Dixon
relation; gcd(u, v) = gcd(u, k) and the duplicate mechanism; v and u cannot share a prime ∤ kN; the Blum functional
t_i = j^Te_i/2 (= our β, Prop B); an explicit p-oblivious degenerate selection rule (= beta0); "no hidden root from
this interface"; uniform-lattice-lift model for R.
**Adopted from codex.** The k ≍ N^θ constant (16(1+θ)/9)^{1/3}; prime k > y as a duplicate-free domain; the N = 77
degenerate example at maximal parity rank; the totality caveat for beta0; the marked-candidate yield model; the
requirement list's availability/uniformity items.
**Ours.** The measurements (E1, E3, E3b, E4), the inheritance/coprimality explanation of the toy H_S failure and its
size trend, Lemma CP at the Prop SI scale, the J+ base degeneracy (Cor. B(a)) and its experimental test, Prop R′.

### 6.2 Adversarial round (`codex/sicheck_2.out`): MAJOR REVISION, all points adjudicated
Codex confirmed the cost optimisation, Lemma I, Prop B (including ℓ = 2, signs, lifts), the sieve threshold, the
dedup, Z = W/V, the elimination (250 matrices against exhaustive kernels), the cofactor choice, the KS and P2
computations, and PH3's arithmetic (0.0981–0.6901).

| # | Sev. | Finding | Adjudication |
|---|---|---|---|
| 1 | MAJOR | Lemma D needs a joint law, not marginals | accepted: relabelled HEURISTIC; model-free route (k with no prime ≤ y) added; "0.6ρ" no longer claimed |
| 2 | MINOR | CP's range overstated for small z | accepted: proof at z ≍ h, monotone lower bound for z ≤ h |
| 3 | MINOR | citation for Ψ_q | accepted: de la Bretèche–Tenenbaum 2005 Thm 2.1; literature critic adds Xuan 1995, Tenenbaum 1993 |
| 4 | MAJOR | PH2 is a benchmark, not a calibrated null; k-conditioning changes it (A: 1.39 → 1.89) | accepted: PH2b reports both nulls; wording "post-hoc benchmark" |
| 5 | MAJOR | PH2 omits fresh relations from g > 1 outputs | accepted: PH2b on all first-seen reduced pairs (pooled 1.008–1.026) |
| 6 | MAJOR | size trend confounded by y, h, K | accepted: trend stated as descriptive, `next` only monotone, no asymptotic inference |
| 7 | MINOR | P1 failure count is 10, not 14 | accepted |
| 8 | MAJOR | rank bound and odd-index sufficiency need subgroup hypotheses (N = 161) | accepted: r ≤ s − rank{c_p, c_q}; subgroup condition stated; C6 fixed |
| 9 | MAJOR | QS explanation wrong | accepted (also the proof critic): QS's root x + m lies outside the base |
| 10 | MINOR | beta0's yield/cost conditional | accepted |
| 11 | MAJOR | hidden-root impossibility is an argument | accepted: "no construction supplied"; costs HEURISTIC |
| 12 | MAJOR | Prop R′ needs the retained sampling law | accepted |
| 13 | MAJOR | C1–C7 need joint, quantified uniformity | accepted: C8 added; C2 wording fixed |
| 14 | MAJOR | E3's 1/2 is algebra once Φ ≢ 0; 36 moduli not 144 tests | accepted |
| 15 | MINOR | large-prime remark overbroad; cite batch testing | accepted |
| 16 | MINOR | τ < 0.757 needs re-optimisation (≈ 0.69 at fixed κ*) | accepted |
| 17 | MINOR | E3 `k_last` is the shared counter | accepted: code records `k_done`; display-only effect noted |
| 18 | MAJOR | §0 overclaims | accepted: §0 restricted to the conditional implication plus exact results |

## 7. Opus critics

### 7.1 Proof/cost critic: ACCEPT WITH FIXES, all applied
It confirmed the cost program, the τ threshold (and its move to κ³ = 4/9 for τ > 0.381), claim (iii), the θ-constant,
Lemma I, Lemma CP's proof (α ≈ 0.82 at 2^1024), Prop B, the general-N specialisation, Prop R′'s inequality, and found
no code bug that changes a number. Fixes: (1) the false "ℓ = 2" explanation of why Cor. B(a) spares QS (its scratch
QS test on 6 Blum N had 6/6 fruitful with a J+ base containing 2); (2) refined H_S relabelled HEURISTIC with the rate
consequence separated; (3) the §5 sentence limited to C4; (4) Lemma I(c) extended to coprimality with k′, with the
averaging remark for PH2; (5) Lemma D's model assumption and the model-free k-coprime route; (6) CP stated as a lower
bound, practical-size caveat; (7) §3.4 relabelled; (8) Prop R′ sampling-law remarks; (9) C2 "at least"; (10) C5
"necessary" points reworded; (11) rank-bound hypotheses; (12) factory row numeric, consistency check not independent;
(13) monotonicity only for `next`; (14) P2's pooled π₀ and R5's weak tolerance noted.

### 7.2 Literature critic: OK WITH FIXES, all applied
It replaced the UNVERIFIED Ψ_q citation with Xuan (Acta Arith. 69, 1995) and Tenenbaum (Phil. Trans. A 345, 1993),
both checked via Granville 2008 §4.2 and both covering ω(q) ≈ x primes ≤ y = L[1/3], so Lemma CP is KNOWN in
cumulative form; completed HT 1986 (Trans. AMS 296, Lemma 2), COS 1986 (Algorithmica 1), Pomerance 1982 for the linear
sieve, Dixon 1981, Kirchner 2018, Coppersmith 1993 for the factory, Bernstein's batch test and Lenstra's ECM; added
**Schnorr's 1991 fac-relations** as the closest prior form of a Prop SI relation (Prop B / Cor. B apply verbatim),
the NFS coprime-pair rule, BLP 1993, Murphy's α and Knuth–Schroeppel as prior art for Lemma I's phenomena; confirmed
that no precise prior statement of Prop B / Cor. B was found; and noted that Boneh 2002 / Granville 2008 §5.7 posed the
"remove the window search" question (polylog-y regime only) and that Jain 2026 (y ≥ exp((log X)^{2/3+ε})) gives no
existence at y = L[1/3]. One nit: α < 1 already for z slightly below y; wording adjusted.

## 8. Reproduction
```bash
.venv/bin/python -m pytest -q tests/test_sicheck.py          # exact model checks, sieve, Prop B, Cor. B, factor_si
.venv/bin/python -m src.sicheck.model                         # the exact program
.venv/bin/python -m src.sicheck.run_experiment E1             # ~25 min, peak 1.5 GB
.venv/bin/python -m src.sicheck.run_experiment E3
.venv/bin/python -m src.sicheck.run_experiment E4
.venv/bin/python -m src.sicheck.score                         # registered predictions
.venv/bin/python -m src.sicheck.run_experiment E3b            # Addendum A1
.venv/bin/python -m src.sicheck.posthoc                       # disclosed post-hoc (PH1-PH3)
.venv/bin/python -m src.sicheck.posthoc ph2b                  # disclosed post-hoc PH2b (after codex round 2)
```
