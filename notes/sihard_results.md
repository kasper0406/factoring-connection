# sihard: the hardness side of Prop SI (smooth numbers in short intervals near kN)

Agent `sihard`, phase 20, 2026-10-21. Version 2: after the codex blind and adversarial rounds and both Opus critics (§10, §11).

> **Erratum (2026-10-27): ETH.** OpenAI, "The Exponential Time Hypothesis is False" (preprint, 17 Oct 2026, rev. 25 Oct 2026), https://github.com/openai/math/blob/main/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/paper.pdf gives a deterministic #k-SAT algorithm in time 2^{O(n log log n/√log n)}. This refutes ETH, randomized ETH, #ETH and Gap-ETH (unrefereed; Theorem 1.1 there). Every ETH-conditional statement below is therefore restated under **ETH⁻**: *for every ε > 0, 3-SAT on n variables with O(n) clauses has no 2^{O(n^{1−ε})}-time algorithm* (randomized ETH⁻ for randomized reductions). On such instances the paper's bound is 2^{O(n/√log n)} = 2^{n^{1−o(1)}}, so it is consistent with ETH⁻. ETH⁻ is stated for O(n)-clause instances because sparsification costs 2^{εn}, which no power-saving hypothesis can absorb. Conclusion here: the blow-up exponent 3/2 is unchanged. A reduction with c < 3/2 would give 3SAT in 2^{O(m^{2c/3}(log m)^{1/3})} = 2^{O(m^{1−δ})}, which contradicts ETH⁻ as well as ETH.


**Files.**
- Code: `src/sihard/`.
  - `core.py`: exact uint64 smoothness sieve, Dickman ρ (design only), Johnson-certificate counts, two-list merge
    counts.
  - `reduction.py`: the X3C → S-smooth-in-interval reduction (Prop XS). `build_instance` is the registered toy
    generator; `build_instance_theory` uses the provable parameters.
  - `consequences.py`: Prop SI and Prop DLI constants, the collapse table.
  - `run_experiments.py`: E1–E4 and the analysis.
- Tests: `tests/test_sihard.py` (20 pass).
- Pre-registration: `notes/sihard_prereg.md`, commit `7f9e398`, plus addendum A1 (`bb7fbd0`), both made before any run.
- Raw data: `experiments/raw/sihard_E{1,2,3,4}.jsonl`. Summary: `experiments/analysis/sihard_summary.json`.
- GPT-6.1-sol: blind round `codex/sihard_1.{md,out}` (§10.1); adversarial round `codex/sihard_2.{md,out}` (§10.2).
- Parallel agents: `siconstruct` (constructions; its note appeared during this work and is cross-referenced),
  `sicheck` (validation of Prop SI; pre-registration only at the time of writing).

**Labels.** KNOWN (cited; checked in the source where stated) · PROVED (proof here) · model-PROVED (exact consequence
of a named model) · HEURISTIC (named assumption) · CONJECTURED · OPEN · MEASURED · POST-HOC · assessment.
"No prior statement found by us" is the strongest novelty claim made.

**Notation.** n = log X (natural log) for the size X ≈ kN of the integers sought; F_α = n^α(log n)^{1−α};
L[α, c] = exp((c + o(1))F_α). Factor base y = L[1/3, κ], window h = L[2/3, ν], u₀ = n/log y (so ρ(u₀) =
L[2/3, −2/(3κ)]), ρ the Dickman function, Ψ(x, y) the count of y-smooth integers ≤ x. Prop SI's optimum is
κ* = (2/9)^{1/3} = 0.6057, ν ↓ 2/(3κ*) = 1.1006, constant (16/9)^{1/3} = 1.2114. A constructor costing L[1/3, τ] per
output still beats MNFS iff τ < 0.757.

---

## 0. Verdict

**Plain language.**
- **No barrier is known for Prop SI in general, and none can be proved with current techniques.** An unconditional
  superpolynomial lower bound would separate P from NP.
- What we can say is narrower. Every *oblivious* construction we could formalise fails completely: it costs at least
  L[2/3, 2/(3κ) − o(1)] per output, the leading cost of sieving the interval. In those models the 1.2114 target
  collapses to L[2/3] per relation, nowhere near MNFS's 1.902. The models are:
  - point-query sampling with independent labels;
  - Johnson-radius certificate decoders whose parameters do not depend on the window's position (Theorem J);
  - divisor completion;
  - fixed product tables (meet-in-the-middle; exponential).
- Two further routes look hopeless but are not proved so. Wagner-style k-trees are argued only in a simplified list
  model. The Schnorr log-lattice route costs exp(L[1/3]) for the generic algorithms we know of; that is an estimate,
  not a lower bound.
- **Not excluded:** adaptive choice of moduli or lattices from the input; correlated polynomial identities of the
  Friedlander–Lagarias type; multi-sample common-divisor decoders (Cohn–Heninger); smooth values of polynomials or
  norms near (kN)^{1/d}; uncertified decoder outputs. We have no evidence that any of them works.
- **The main new result is Theorem J (PROVED, averaged over where the window sits).** Take any single-integer
  Coppersmith / Howgrave-Graham / Boneh / Guruswami–Sahai–Sudan certificate decoder, with any modulus built from
  factor-base prime powers, any weights and any lattice dimension, applied to one window or batched over many
  multipliers k.
  - With its parameters fixed in advance, one call certifies on average at most L[2/3, −2/(3κ) + o(1)] smooth
    integers. That is the yield of L[2/3]^{o(1)} random tests.
  - The trivial divisor trick z = d·⌈x/d⌉ attains this at leading order. So at leading L[2/3] order the whole
    single-integer Johnson class is no better than the trick.
  - It is a quantitative, averaged form of Boneh's 2000 remark that "a random interval … is unlikely to contain
    integers with sufficiently large smooth factors".
  - Codex proved the fixed-modulus case independently, by a different method.
  - The averaging replaces any random-divisibility heuristic. Only certified outputs are counted, and only call
    schedules fixed in advance are covered, or call types drawn from a repertoire of size L[2/3]^{o(1)}.
- **Consequences beyond factoring (Prop DLI).** A constructor that works for every modulus would also give
  prime-field discrete logarithms in L_p[1/3, 1.232] (heuristic).
  - Collecting relations costs 1.2114.
  - The constructor also performs the descent of medium-sized primes, at L[1/3, ≈ 1], when queried with the modulus
    ⌊jp/q⌋. This is Granville's t_j = [jp/q] device.
  - The classical "booting" step dominates: 1.232 with Barbulescu's early abort (as stated in
    Fried–Gaudry–Heninger–Thomé 2017, App. A.1; checked), or Commeine–Semaev's 3^{1/3} = 1.442 with plain ECM.
  - For general primes, prime-field DL has stood at 1.902 (Matyukhin 2003). So a modulus-generic constructor would
    break two problems believed hard.
  - No implication was found for class groups or for other groups.
- **The existence half alone is deep (Prop QNR, codex; Granville's argument).** If every interval at Prop SI's
  scale contains a smooth integer, then for primes p ≡ 3 (mod 4) the least quadratic nonresidue is
  ≤ L_p[2/3, ν + ε] = p^{o(1)}. That is Vinogradov's conjecture for those primes. So only *promise* constructors are
  realistic targets, and their promise cannot be proved with current methods.
- **Factoring does not obviously help.** In an independent-label point-query model, which is a surrogate and not a
  real factoring-oracle separation, finding one smooth integer needs ≥ 1/δ tests, δ being the density (Lemma BB).
  Free factoring speeds up testing, not selection.
- **NP-hardness does not reach Prop SI's scale** (§6).
  - "Smooth over an explicit list of primes, in an interval" is NP-complete under randomized reductions. That follows
    from the known hardness of DIVISOR (Alon–Kilian; Shor's 2011 proof). Our Prop XS extends it to unbounded
    exponents and makes it parsimonious.
  - For the real problem (all primes ≤ y) no NP-hardness is known in any range.
  - It is provably not NP-hard unless NP = coNP in three cases:
    - when h ≤ poly(log x);
    - where every interval provably contains a smooth integer;
    - at Prop SI's scale, if its existence half (itself Vinogradov-hard) holds.
  - Unconditionally the sieve solves Prop SI's scale in L[2/3]. So a deterministic reduction from 3SAT must blow up
    instance size by a power ≥ 3/2 unless ETH⁻ fails (stated under ETH in v2; ETH itself is refuted, see the erratum).
- **Assessment.** Probability that an L[1/3, τ] constructor with τ < 0.757 exists: about 1% (ours), about 5%
  (codex).
  - The oblivious models fail for a shared reason (siconstruct's Lemma SR): with the other factors fixed, every factor
    of u inherits the relative precision h/x = N^{−1+o(1)}.
  - Correlated constructions, such as (m − a)(m + a) = m² − a², are not covered by that observation.

**Main results.**

| Result | Status |
|---|---|
| Lemma EQ: a modulus-uniform Prop SI constructor = an arbitrary-interval constructor at the same scale | PROVED (trivial) |
| Prop DLI: modulus-generic constructor ⇒ F_p DL in L_p[1/3, max(1.2114, c_boot)], c_boot = 1.232 (early abort) or 3^{1/3} (ECM); descent 0.99–1.05 | HEURISTIC (H_S, H-affine, rank mod ℓ \| p−1) |
| Prop QNR: all-interval existence at Prop SI scale ⇒ least QNR ≤ p^{o(1)} for p ≡ 3 (4) | PROVED (codex; Granville's argument) |
| Remark FO: P^{FACT} ⊆ NP ∩ coNP | PROVED (codex) |
| Lemma BB: independent-label point-query model needs ≥ (1−o(1))/δ queries per output (a surrogate, not a factoring-oracle separation) | model-PROVED |
| Theorem J: single-integer Johnson-class certificates (weighted, any P, batched over k ≤ h) average ≤ L[2/3, −2/(3κ)+o(1)] smooth certified outputs per call; sharp at leading order | PROVED (unweighted, averaged over x or N, parameters fixed in advance); weighted PROVED (sketch); codex independent proof for fixed P |
| Prop S: Prop SI-type relations in sampling model S cost ≥ L[1/2, √2] | model-PROVED (suffmap Theorem SS(a); assumes H_S, R and no prefiltering, so it excludes MNFS-type filtering) |
| Prop MM: fixed two-list product tables, interval and Bohr (all k at once): ≥ √(X/3)/(h+1) per output, resp. ≥ √(X/K)/(h(1+ln K)+K) (exponential for K = X^{o(1)}) | PROVED (averaged); siconstruct Lemma MA under H_eq |
| Wagner/k-tree bit budget: cost ≥ √X/h | HEURISTIC sketch (simplified list model without representation multiplicities) |
| Prop NC / TF / ETH: not NP-hard (h ≤ poly, proven-existence ranges, Prop SI scale under Hyp E); deterministic-reduction blow-up ≥ m^{3/2} under ETH⁻ (was ETH; erratum) | PROVED (Megiddo–Papadimitriou plus the transcript argument for Cook reductions) |
| Prop XS: explicit-S version NP-complete (randomized reductions), parsimonious; extends the known DIVISOR hardness (Alon–Kilian; Shor 2011) to unbounded exponents | PROVED; E4, E4b 300/300 |
| Lemma FP: h ≥ 2x/(y−1) is in FP | PROVED (trivial) |

**Experiments** (pre-registered, §8).
- E2 (P4, exact), E3 (P5 and P6) and E4 (P7, 300/300) passed.
- E1:
  - **P1 failed in 1 of 6 cells**: the near-kN/control mean ratio was 0.86 at y = 2^10, μ = 2.
  - **P2 failed in 4 of 6 cells**: overdispersion. Diagnosed post hoc as size heterogeneity across k; the binned
    dispersion index is 1.00–1.02.
  - P3 passed.
  - Post hoc: intervals placed along kN show a modulus-dependent bias of either sign, about ±10% at 48 bits.

---

## 1. The problem, and what counts as a barrier

- **SMOOTH-INT(x, h, y).** Given x, h, y, output a y-smooth integer in [x, x + h], or (decision) say whether one
  exists. Prop SI asks for it at scale (h, y) = (L_x[2/3, ν], L_x[1/3, κ]) with ν > 2/(3κ), at x = kN, in
  L^{o(1)} time, or in L[1/3, τ] time with τ < 0.757.
- **Lemma EQ (PROVED).** If a constructor works for every modulus M and every k ≤ L[1/3, C], then calling it with
  M := x, k := 1 solves SMOOTH-INT(x, L_x[2/3, ν], L_x[1/3, κ]) for every x. Conversely an arbitrary-interval
  constructor serves the intervals I_k = [kN, kN + h]. ∎
  - The equivalence needs uniformity in the modulus. A constructor that uses the shape of N (for example that N is a
    semiprime) is not covered. No such mechanism is known.
- **The trivial costs.**
  - Sieving [x, x + z] costs ≪ z log log y + u·y (Granville 2008, §5.7; KNOWN, checked in the text). One output
    needs z ≈ 1/ρ(u₀) = L[2/3, 2/(3κ)].
  - Granville: any such algorithm likely needs time ≥ π(y), "since the primes up to y are part of the input". He then
    writes that "one might guess that it is not necessary to have z in the running time", and presents Boneh's
    (y log x)^{O(1)} method. Prop SI asks for exactly that, at its scale. The
    π(y) term does not bound the cost per call, because the factor base is shared across calls.
- **What a barrier could be.** An unconditional superpolynomial lower bound for an NP search problem implies P ≠ NP
  (suffmap §4). So barriers here are of three kinds:
  - (i) model lower bounds (§5);
  - (ii) hardness reductions, with their limits (§6);
  - (iii) consequences that would break other problems believed hard (§2).

## 2. Consequences of a constructor (Task 1a)

### 2.1 Factoring
L[1/3, 1.2114 + ε] (suffmap Prop SI; heuristics H_S on the residue and R). With per-call cost L[1/3, τ] the constant
is min_κ max(κ + τ + 2/(9κ²), 2κ), below MNFS iff τ < 0.757 (`consequences.prop_si_optimum`; siconstruct §1 table).

### 2.2 Prime-field discrete logarithms (Prop DLI) [HEURISTIC; no prior statement found by us]

**Hypothesis.** 𝒪 is modulus-generic: for every integer M ≥ M₀ and every k ≤ L_M[1/3, C] it outputs a y-smooth
u ∈ [kM, kM + h] in L^{o(1)} time, with y = L_M[1/3, κ] and h = L_M[2/3, ν]. It may be a promise constructor, and
the promise must hold on the queried intervals.

**Claim.** Assume:
- H_S: residues are smooth like random integers of their size;
- H-affine: as j and k vary, the affine residuals qu − kjp of step 2 are smooth at that same rate, despite their
  congruence conditions modulo q (codex's name).
  - H_S and H-affine are hypotheses on the integers *the oracle chooses*. An adversarial 𝒪 picks u, and with it
    u − kM, so these are assumptions about 𝒪, not facts about random integers;
- the usual index-calculus rank heuristic modulo every ℓ^e ∥ p − 1, with relations normalised by expressing a
  generator g in the factor base;
- ECM's heuristic running time.

Then discrete logarithms in F_p^* cost L_p[1/3, max((16/9)^{1/3} + ε, c_ind)].
- With the classical ECM booting, c_ind = 3^{1/3} = 1.4422.
- With early-abort booting, c_ind = 1.232. *Checked:* Fried–Gaudry–Heninger–Thomé (EUROCRYPT 2017; arXiv 1610.02874,
  §5 and App. A.1) state B_init = L_p(2/3, 0.811) and T_init = L_p(1/3, 1.232). These come from Barbulescu's thesis,
  Ch. 4: early abort with 6 ECM filters, numerically optimised. FGHT's "Fact" box misprints the bound as L_p(1/3, 0.811).
  FGHT also state the plain-ECM value 3^{1/3} at c = 1/∛3. The initialization does not use the NFS polynomials, so it
  transfers.
- So the total is **≤ L_p[1/3, 1.232]**, under the heuristics above plus Barbulescu's early-abort analysis, and
  ≤ L_p[1/3, 1.4422] with plain ECM booting.

**Proof sketch.**
1. *Relations.* Run Prop SI with N := p. Each u ≡ v (mod p), u, v y-smooth, is a linear relation among factor-base
   logarithms modulo p − 1. The cost is 1.2114 + ε, and the linear algebra mod p − 1 costs L[1/3, 2κ].
2. *Descent of a medium prime q = L_p[2/3, c].*
   - Pick j ≥ 1, set M := ⌊jp/q⌋, and query 𝒪(M, k) for k = 1, 2, …, getting y-smooth u ∈ [kM, kM + h].
   - Then v := qu − kjp = q(u − kM) − k(jp mod q) satisfies |v| ≤ q(h + k) = L_p[2/3, c + ν], and v ≡ qu (mod p).
   - If v is y-smooth (probability L[1/3, −(c + ν)/(3κ)] under H_S), then log q = log v − log u.
   - M ≈ p/q has L_M = L_p to leading order, so 𝒪 is used at its own scale.
   - Cost per medium prime: L[1/3, (c + ν)/(3κ)]. Testing v for y-smoothness by ECM costs L^{o(1)}.
3. *Booting.* Write the target as t·g^r mod p and wait until it is L_p[2/3, c]-smooth.
   - The probability is L[1/3, −1/(3c)]. Finding the medium factors by ECM costs L_{q}[1/2, √2] = L_p[1/3, √(4c/3)]
     per candidate.
   - The cost min_c 1/(3c) + √(4c/3) = 3^{1/3} is at c = 3^{−1/3}. KNOWN: Commeine–Semaev, PKC 2006 (LNCS 3958),
     as credited in FGHT App. A.1; `consequences.dl_booting_ecm`.
   - At that c the descent costs (0.693 + 1.1006)/(3·0.6057) = 0.987 < 1.442. So booting dominates:
     c_ind = 3^{1/3}.
4. *Early abort.* FGHT's booting (B = L[2/3, 0.811], T = L[1/3, 1.232]) does not use special NFS polynomials, so it
   transfers. The descent at c = 0.811 then costs (0.811 + 1.1006)/(3·0.6057) = 1.052 < 1.232. So
   c_DL = max(1.2114, 1.232) = 1.232. Codex computed the same 1.052. Neither 1.442 nor 1.232 is a proved bottleneck.
5. *Can 𝒪 also boot?* Booting needs t·s ≡ w (mod p) with s, w both smooth. The interval oracle offers only additive
   control, u ∈ [t + jp, t + jp + h], so it gives u ≡ t + v, not a multiplicative relation. Multiplicative control
   would be a Bohr-set (coset) lister, suffmap E2 with Prop DLT. We found no way to make 𝒪 boot.

**Remarks.**
- This sharpens suffmap §4's unchecked remark: the proof critic's "≈ L[1/3, 0.83]" is the descent without the ECM
  cost of booting (`consequences.dl_descent_optimum`: 0.827).
- The I_k-only constructor also suffices, *provided* it is uniform over moduli: the descent queries modulus
  ⌊jp/q⌋, not p. It does **not** need suffmap's coset hypothesis.
- Codex's blind round obtained the same construction (with k = 1 and M = ⌈jp/q⌉) and the same constants
  independently: boot 1.4422, descent 0.987, and 1.052 with the early-abort boot (§10.1).
- *Linear algebra and normalisation (codex).*
  - Express the generator g as well as the target: boot and descend g symbolically.
  - Signs are irrelevant (u, v > 0).
  - Solve the sparse system modulo the prime powers dividing p − 1. Factor p − 1 with Prop SI's own factoring
    algorithm (1.2114, under its heuristics). Alternatively, run Lanczos/Wiedemann modulo the unfactored part and
    split it whenever a non-unit pivot appears. Small prime factors of p − 1 go to Pohlig–Hellman.
  - Enough trials come from varying j and k.
- *Forerunners.* Granville's proof (§4.4) that smooth numbers in short intervals imply Vinogradov uses the same
  device: t_j = [jp/q] and q(t_j + i) − jp. Schnorr (EUROCRYPT 1991, LNCS 547) already used relations with u and
  |u − vN| smooth both for factoring and for discrete logarithms in F_p. The specific statement "interval oracle ⇒ DL
  at ≤ L[1/3, 1.232]" was not found by us or the literature critic.
- *Reading.* For general primes (no special form), prime-field DL has resisted any improvement on L[1/3, 1.902]
  (Matyukhin 2003, the constant ((92 + 26√13)/27)^{1/3} per Barbulescu–Pierrot 2014; see also Commeine–Semaev 2006). A
  modulus-generic Prop SI constructor would bring it to ≤ 1.232 (≤ 1.442 with plain ECM booting) and
  factoring to 1.2114. That is evidence against
  such constructors, of the same kind as suffmap's Prop DLT (assessment).

### 2.3 Class groups, other groups [no implication found; assessment]
- *Class groups of imaginary quadratic fields.* Relations need smooth *norms of elements of a given ideal*: values of
  a positive binary quadratic form f(a, b) on a 2-D lattice. An integer constructed in an interval is not represented
  by f, except with probability ≈ (number of representable integers)/(interval). A useful hypothesis would be
  "smooth values of f in boxes", which is the NFS-row type C2/C3 of suffmap, not Prop SI. Hafner–McCurley's
  L[1/2] is untouched.
- *Elliptic curves, Jacobians, F_{p^n} in small or medium characteristic.* Relations there do not pass through
  integers in intervals. No implication.
- *DL modulo composite N.* That follows from factoring plus prime-field DL (Bach 1984), so it inherits 2.1 and 2.2.

### 2.4 Smooth numbers in intervals at other scales [PROVED, elementary]
- *Monotonicity.* A constructor at (κ, ν) gives one at every (y′, h′) ≥ (y, h) at the same x: query a sub-interval.
  It gives nothing at smaller y or shorter h.
- *Consequence.* A deterministic constructor that is correct on *every* x (non-promise, via Lemma EQ; this entails
  Hyp E of §6) would give, deterministically in time x^{o(1)}, x^ε-smooth integers in
  [x, x + C√x]. That is Matomäki–Radziwiłł's existence range, and also Harman's and Balog's ranges [x, x + x^{1/2+ε}].
  - Randomized constructions there are easy. Matomäki–Radziwiłł Cor. 1 gives at least √x(log x)^{−4} of them in
    [x, x + C√x], density ≥ 1/polylog, so sample and factor with ECM or NFS in subexponential time. Global density
    alone would not suffice; their count is what guarantees the local density. The open part is determinism. Deterministic factoring is exponential (N^{1/5}, Harvey 2021), so
    the trivial method does not derandomize.
  - The constructive results known to us are much weaker:
    - sieving (exponential);
    - Boneh (y ≤ polylog x, time (y log x)^{O(1)}; §4);
    - Friedlander–Lagarias 1987, as reported by Granville §4.1: x^{1/2}-smooth integers in every interval of
      length x^{1/4+ε}, by explicit construction (UNVERIFIED in the primary source);
    - the trivial Lemma FP below.
  - We found no deterministic subexponential construction in the MR range.
- **Prop QNR (PROVED; codex blind §1(a)).** Suppose Hyp E(κ, ν) holds: every [x, x + L_x[2/3, ν]] with x ≥ x₀ contains
  an L_x[1/3, κ]-smooth integer. Then for every prime p ≡ 3 (mod 4) large enough, the least quadratic nonresidue is
  ≤ 2H with H = L_p[2/3, ν + ε].
  - *Proof.* Take a smooth u ∈ [p − 2H, p − H]. If every prime ≤ 2H were a residue, then u and p − u ≤ 2H would
    both be residues. But u ≡ −(p − u) and −1 is a nonresidue. ∎
  - Granville (2008, §1.2, §4.4) gives the x^ε version and calls Vinogradov's conjecture "intrinsically difficult".
    So Hyp E is out of reach of current methods; only *promise* constructors are meaningful targets.
  - The argument uses all-interval existence, not the constructor.
- **Lemma FP (PROVED).** If h ≥ 2x/(y − 1), then z = 2^a·⌈x/2^a⌉ with 2^a ∈ [x/(y−1), 2x/(y−1)] is a y-smooth integer
  in [x, x + h], computed in polynomial time.
  - *Proof.* z − x < 2^a ≤ h, and ⌈x/2^a⌉ ≤ y − 1 + 1 = y. ∎
  - Below h ≈ x/y, divisor completion stalls: the relative precision h/x is inherited by every factor (siconstruct
    Lemma SR).
- *Primes.* Smooth numbers are not primes. The constructor gives nothing for Polymath4 (deterministic primes), whose
  "Strategy 3" needs *non*-smooth numbers in intervals (polyimpl §3).

### 2.5 Map

| Target | Implied by a modulus-generic Prop SI constructor? | Constant / range | Status |
|---|---|---|---|
| factoring | yes | L[1/3, 1.2114] | HEURISTIC (H_S, R) |
| F_p discrete log | yes | ≤ L[1/3, 1.232] (early-abort boot; 1.442 with plain ECM); relations 1.2114; descent 0.99–1.05 | HEURISTIC (Prop DLI) |
| least quadratic nonresidue, p ≡ 3 (mod 4) | by the *existence* half (Hyp E), not the constructor | ≤ L_p[2/3, ν + ε] | PROVED (Prop QNR) |
| DL mod composite N | yes | via Bach 1984 | HEURISTIC |
| smooth-in-interval at larger (y, h) | yes, deterministic if 𝒪 is | e.g. x^ε-smooth in [x, x+√x] deterministically in time x^{o(1)} | PROVED (monotone) |
| class groups (imaginary quadratic), units | no implication found | — | assessment |
| EC / genus-g / small-char DL | no | — | assessment |
| deterministic primes | no | — | assessment |
| rigorous factoring | no (residue randomness is heuristic) | — | assessment |

## 3. The converse: is Prop SI implied by factoring, or by anything known? (Task 1b)

**Lemma BB (model-PROVED; the independent-label point-query model).** Model: the algorithm learns whether z is y-smooth, or its full factorisation, only by
querying z, and it may do arbitrary other computation. The smooth set is replaced by a random set in which each z
in [x, x + h] is included independently with probability ≤ δ. Then T queries find a member of [x, x + h] with
probability ≤ (T + 1)δ. The +1 covers outputting an unqueried z.
- *Proof.* Union bound over the queried points and the output point; unqueried points are independent of the
  transcript. ∎
- *Reading.* With δ = ρ(u₀)(1 + o(1)), one output needs ≥ L[2/3, 2/(3κ)]/2 queries for success probability 1/2.
- *Limits (codex).* This is a surrogate, not a factoring-oracle separation.
  - A real factoring oracle is not an independent-label oracle: its answers obey multiplicative identities, and
    queries outside the interval can reveal information about many candidates.
  - Hildebrand's density estimate is global; it does not justify independent labels or a uniform local density.
  - What the lemma shows is that speeding up *point tests* (factoring included) cannot by itself reach Prop SI.
- The interval sieve already tests h integers in h·O(log log y) operations. So even free factoring improves
  the cost per examined integer by at most O(log log y) (compare polyimpl Lemma PF, where a remainder tree replaces
  factoring).
- The model ignores arithmetic: real algorithms can *construct* smooth numbers (products). §5 treats those.
- **Remark FO (PROVED; codex).** A full-factorisation oracle answer can be certified: the prime factorisation, with
  AKS. So every language decidable in polynomial time with a factoring oracle lies in NP ∩ coNP.
  - So if the *decision* version of SMOOTH-INT at some scale were NP-hard, no polynomial-time algorithm with a factoring
    oracle could decide it unless NP = coNP.
  - At Prop SI's scale NP-hardness is not known (§6), so this is a conditional separation only.

**Known constructive results do not reach Prop SI's scale.**
- Lemma FP needs h ≥ x/y.
- Friedlander–Lagarias type identities reach h = x^{1/4+ε}, still exponential.
- Boneh's decoder lists, and it lists nothing in random intervals at this scale (§4, Theorem J).
- Existence theorems (Hildebrand; Friedlander–Granville; Balog; Harman; Croot; Matomäki; Matomäki–Radziwiłł;
  Jain 2025, arXiv 2502.10530) are non-constructive, and none reaches y = L_x[1/3] in all intervals of length
  L_x[2/3]. Granville §4.1: the asymptotic (4.1) is known for exp((log x)^{5/6+o(1)}) ≤ y and
  √x·y²·exp((log x)^{1/6}) ≤ z ≤ x (Friedlander–Granville 1993).
- **Conclusion.** Prop SI is implied by nothing known. Polynomial-time factoring is not known to imply it. In the
  point-query surrogate it does not, and no real factoring-oracle separation is known. Conversely Prop SI gives factoring only at L[1/3, 1.21], not
  polynomial time. The two are incomparable as far as is known.

## 4. Literature status (Task 2)

- **Sieving.** Granville 2008 §5.7 (MSRI 44, checked): running time ≪ z log log y + uy for [x, x + z].
- **Boneh 2000/2002** (STOC 2000; JCSS 64 (2002) 768–784; Theorem 3.1 checked in the author's PostScript).
  - Boneh's logarithms are base 2, which matters for the constants 4|I| and ¼ below.
  - Setting: S = ∏_{q<s} q^{⌊log s/log q⌋}. Every *strongly* s-smooth integer (p^m ∤ N whenever p^m > s) divides S.
  - For any d ≥ 1 + √(log S/log|I|), it finds in poly(d, size) time all N ∈ I with gcd(N, S) > T, whenever
    T > S^ε with ε = √(log(4|I|)/log S) + 5/(4d). Approximately |I| < ¼·T^{(log T/log S) − 2.5/d}; a (5/(4d))² term is dropped.
  - If moreover V < 2T, every returned integer divides S, so it is strongly s-smooth. *Completeness* (every strongly
    s-smooth integer of I = [U, V] is returned) also needs U > T, i.e. T < U ≤ V < 2T (codex). The cited statement
    omits this; the gcd-qualified guarantee is what fits Theorem J's class.
  - Granville's summary: for (1/5)(log x)² > y > log x, there are ≪ y/log x integers in [x, x + x^{(log x)/(4y)}]
    dividing lcm[1..y]. Boneh: "a random interval I satisfying the bounds of Theorem 3.1 is unlikely to contain
    integers with sufficiently large smooth factors". §4, about the *generalized* CRT problem (smooth values f(x) of a
    polynomial, Theorem 4.1), says the bounds "are insufficient for improving the basic quadratic sieve".
  - **At Prop SI's scale** (log x ≈ n, log S ≈ y = L[1/3]): log|I| ≲ n²/y ≪ 1, so |I| < 1 and it is vacuous for full
    smoothness. Partial smoothness is Theorem J: no gain beyond L[2/3]^{o(1)}.
  - Boneh also dismisses the subset-sum route: "the dimension of the lattice … is exponential in the size of the
    input (log s)".
- **Coppersmith / Howgrave-Graham / GSS weighted CRT decoding.**
  - These give the same Johnson-type guarantee: (log T)² ≥ log V·log P, writing T for the gcd threshold to avoid a
    clash with Boneh's S.
  - Weighted: Σw·agreement ≥ √(Q log V) (GSS FOCS 2000, Thm 4, with an extra tolerance term).
  - Herrmann–May bivariate linear: ∏X_i < P^{g(μ)} with g(μ) ≤ μ² (pilist §3.1), as a certification region under
    their Assumption 1.
  - These single-integer certificates are in Theorem J's class.
  - **Not in it:** Cohn–Heninger's multi-sample approximate common divisors (ANTS X, 2013; errors up to
    P^{β^{(m+1)/m}}), and Sarkar–Maitra's Jochemsz–May variant for smooth integers in short intervals (ePrint
    2010/146, not checked). See §5.2.
- **Schnorr's log-lattice.** siconstruct Prop CV (PROVED there, with slack factors): a squarefree Prop SI output is
  a CVP_{1+O(1/n)} answer in the π(y) = L[1/3]-dimensional prime lattice, and conversely up to an O(√n) enlargement of
  the window and O(1) changes of multiplier.
- **Hardness literature.**
  - For the real problem (all primes ≤ y) we found no NP-hardness result, and no statement that its complexity is
    open.
  - For an explicit prime list it is known. DIVISOR (a divisor of N in [L, U]) is NP-hard under randomized
    reductions (Alon–Kilian, credited in Sudan 2010; proof by Shor, cstheory.SE Q. 4769, 2011), with primes whose sizes
    encode subset-sum weights. That is the smooth-over-an-explicit-list problem with exponents ≤ 1 (§6, Prop XS).
  - Granville §5.7 suggests the interval length need not enter the running time, and presents Boneh.
    Naccache–Shparlinski 2009 (§6.2, checked by the literature critic) call constructing smooth numbers in a given
    interval "a challenging problem". They tie it to signature forgery, where smooth values a·μ(m) − b·N near
    multiples of N are sought (Coron–Naccache–Stern 1999; Coron–Naccache–Tibouchi–Weinmann 2009). Those are Prop SI's
    relation shape.
  - Related NP-complete problems (KNOWN):
  - Subset Product (Garey–Johnson 1979, SP14, credited to Yao, by reduction from X3C; X3C is SP2, Karp 1972);
  - Manders–Adleman 1978 (bounded quadratic congruences; used in floorhard Theorem SS);
  - subset sum.
  None of them has the all-primes-≤-y structure.
- **More existence results** (from codex's blind round; statements as reported there, not re-checked by us):
  - Hildebrand 1986 (JNT 22, Thm 3): the asymptotic for h ≥ x/y^{5/12};
  - Croot 2007 (IJNT 3): x^{47/(190√e)+ε}-smooth integers in [x, x + C√x];
  - Matomäki 2016 (IJNT 12, Thms 1.1, 1.3): y ≥ exp((log x)^{2/3}(log log x)^{4/3+ε}), all intervals of length
    ≈ √x·exp(…);
  - Younis 2024 (arXiv 2409.05761): h ≥ x^θ, θ > 17/30, y ≥ exp((log x)^{2/3+ε});
  - Jain, "Existence of smooth numbers in short intervals" (QJM 2026, doi 10.1093/qmath/haag010; arXiv 2502.10530):
    almost all intervals at y ≥ exp(C(log X)^{2/3}(log log X)^{4/3}) (Thm 1.1), and an all-interval version with
    h ≥ √x·exp(…) (Thm 1.2);
  - Matomäki 2010 (IJNT 6): x^{1/(5√e)+ε}-smooth integers in intervals of length C√x (stronger than Croot).
  None reaches y = L[1/3].
- **The computational problem is recognised.** Naccache–Shparlinski, "Divisibility, smoothness and cryptographic
  applications" (2009, §6.2; arXiv 0810.2067), describe construction under an interval constraint as challenging,
  with Boneh's restricted parameters (codex; UNVERIFIED by us).
- **Existence.** See §3. Granville's "Challenge problem 2000" (x^α-smooth in [x, x + √x] for all α > 0) is
  answered up to a constant factor in the interval length by Matomäki–Radziwiłł 2016, Corollary 1 (as
  cross-checked by suffmap's literature critic).

## 5. Model lower bounds (Task 3)

Scale throughout: x ∈ [X, 2X), y = L_X[1/3, κ], window length V ≤ L_X[2/3, ν_max]. "Per output" means expected
cost per y-smooth integer found in the target window(s).

### 5.1 Model S (sampling / sieving) [model-PROVED]
- **Lemma BB** (§3): ≥ (1 − o(1))/ρ(u₀) = L[2/3, 2/(3κ)] tests per output.
- **Prop S** (assumes suffmap's H_S, its heuristic R, and model S's ban on prefiltering; it therefore excludes
  MNFS-type shared filtering). In suffmap's model S, a Prop SI-type relation has smoothness integers u ≍ kN and v ≤ h. So
  λ = log u + log v ≥ n, and Theorem SS(a) gives total cost ≥ exp((√2 − o(1))√(n log n)) = L[1/2, √2]. This is
  siconstruct row 0.
- Amortisation over k is irrelevant here: each examined integer is smooth with the same probability.

### 5.2 Model J (single-integer Johnson-radius certificate decoders): Theorem J [PROVED, averaged]

**The class 𝒥.**
- A call is (P, w, R).
  - P = ∏_{q≤y} q^{e_q} ≥ 2, any exponents.
  - Weights w_q ≥ 0, not all zero.
  - R = [x, x + V] is the region.
- Score W(z) = Σ_q min(v_q(z), e_q)·w_q log q, with Q = Σ_q e_q w_q² log q.
- **Certified set** C(R) = {z ∈ R : W(z) ≥ √(Q log V)}.
- *Unweighted* (w ≡ 1): certified ⇔ (log gcd(z, P))² ≥ log V·log P. This is the idealised Coppersmith–Howgrave-Graham
  bound.
- Boneh's Theorem 3.1 has a strictly larger threshold, so his guaranteed set is contained in ours. GSS's weighted
  theorem has an extra tolerance and coprime moduli; the prime-power score is a GSS-type extrapolation, as in
  pilist §2.2.
- We grant the decoder its entire certified set at unit cost. Uncertified roots it may also return are not counted
  (pilist's convention). Dimension, shifts and multiplicities are arbitrary; only the guarantee enters.
- *Obliviousness.* (P, w, V) are fixed independently of x, or chosen from a family 𝒫 fixed in advance.
- *Larger regions do not help.* The threshold √(Q log V) increases with V. So for a window W ⊂ R,
  C(R) ∩ W ⊆ C(W), and it suffices to treat regions no longer than the target window. Smooth outputs outside the
  windows are useless to Prop SI. A region meeting several windows yields at most (number of windows)·(the bound);
  with ≤ L[1/3, C] windows this is absorbed.
- *Primes of P.* WLOG P uses only primes ≤ y. A prime > y never divides a smooth z, and only raises π or Q, i.e. the
  threshold. So "any modulus" is covered.
- *Adaptive call sequences.* Earlier outputs, including certified non-smooth ones, carry information about x. They
  are covered only through (c), by a union over the possible transcripts. A general adaptive version is not proved.

**Theorem J.** Fix κ, ν_max > 0. Uniformly over all oblivious (P, w) and all V ≤ L_X[2/3, ν_max]:

  (1/X)·Σ_{x=X}^{2X−1} #{z ∈ C([x, x + V]) : z y-smooth} ≤ exp(−(2/(3κ) − o(1))·F_{2/3}) = ρ(u₀)^{1−o(1)}.

- (a) The same holds for the I_k-shaped regions [kN, kN + h] averaged over N ∈ [X, 2X) (any fixed k ≤ h).
- (b) It also holds for batched regions R_N = {kN + t : k₀ ≤ k < k₀ + K, 0 ≤ t ≤ h} with k₀ + K − 1 ≤ h (every k in
  the batch must be ≤ h; for k > h one would need smooth numbers in arithmetic progressions), averaged over N, with
  V := K(h + 1) in the threshold. This is the Herrmann–May-type bivariate guarantee.
- (c) For an algorithm choosing P from a family of size |𝒫|, multiply the bound by |𝒫|.
- (d) **Sharpness.** The divisor trick attains this bound. Take d = 2^{⌊log₂h⌋}, P = d, w = 1, and query a
  sub-window of length V := d, so certified ⇔ d | z. (With d < V nothing is certified, since
  √(log V·log d) > log d; proof critic.) The averaged yield is (1 + o(1))ρ(u₀ − log d/log y) =
  ρ(u₀)·L[1/3, 2ν/(3κ)] = ρ(u₀)^{1−o(1)}, using log d ∼ log h (codex).

**Proof (unweighted).** Write t = log V, π = log P and σ = √(tπ). Certified means log gcd(z, P) ≥ σ.
1. *Averaging (exact).* Each z lies in at most V + 1 of the regions [x, x + V]. So the left side is
   ≤ ((V+1)/X)·Z, where Z = #{z ∈ [X, 2X + V] : log gcd(z, P) ≥ σ, z y-smooth}.
2. *Divisor decomposition.* Write z = d·m with d = gcd(z, P). Then m is y-smooth, so
   Z ≤ Σ_{d|P, log d ≥ σ, d ≤ 3X} Ψ(3X/d, y) (only d ≤ 3X can occur). The sum is empty unless σ ≤ n + 2, and
   nonempty only if π ≥ σ, i.e. π ≥ t. Put c := π/t ≥ 1, so σ = t√c. Dickman arguments below are clamped at 0
   (ρ := 1 on (−∞, 1]). Bounded P, with π = O(1), is trivial: then σ = O(1), so t = O(1), and step 5 applies.
3. *Ψ bound.* Ψ(Y, y) ≤ C₀·Y·ρ(log Y/log y) for 1 ≤ Y ≤ 3X. This is Hildebrand 1986, uniform for
   y ≥ exp((log log Y)^{5/3+ε}), which y = L[1/3] satisfies; ρ := 1 on [0, 1]. Group the d by
   j = ⌊log d⌋ ≥ ⌊σ⌋ and let N_j = #{d | P : j ≤ log d < j + 1}. Then
   LHS ≤ 3C₀(V + 1)·Σ_{j≥⌊σ⌋} N_j e^{−j}·ρ((n − j − 2)/log y).
4. *Three counts.*
   - (a) N_j ≤ τ(P) ≤ exp(1.538·log 2·π/log π) (Nicolas–Robin 1983).
   - (b) N_j ≤ #{m ≤ e^{j+1} composed of the first ω(P) primes} ≤ e^{0.3(j+1)}·∏_{i≤ω(P)}(1 − p_i^{−0.3})^{−1}
     ≤ e^{0.3(j+1) + C₂π^{0.7}}. The first step maps the i-th prime of P to the i-th prime, an injection that does
     not increase the value. The second is Rankin's trick. The last uses ω(P) ≤ π/log 2.
   - (c) N_j ≤ e^{j+1}.
5. *Small regions.* If t ≤ t₀ := F_{2/3}/log n, use Z ≤ Ψ(3X, y): LHS ≤ 3C₀(V+1)ρ(u₀ − O(1)) = ρ(u₀)^{1−o(1)}.
6. *Case A: 1 ≤ c ≤ (log t)²/4.* Use (a). The exponent of (V+1)N_j e^{−j} is
   ≤ t(1 − √c + 1.07c/log t) + O(1). As a function of a = √c ∈ [1, (log t)/2] the bracket 1 − a + 1.07a²/log t is convex, so it is at
   most its larger endpoint value: 1.07/log t at a = 1, and 1 − 0.23 log t < 0 at the other end. So the factor is
   ≤ exp(1.07t/log t + O(1)) = exp(o(F_{2/3})).
   - Every d here has log d ≤ π ≤ t(log t)² = o(n/log n), so ρ((n − j − 2)/log y) ≤ ρ(u₀(1 − 1/log n)) = ρ(u₀)^{1−o(1)}.
   - The number of j is ≤ π + 1 = poly(n).
7. *Case B: c > (log t)²/4.* Use (b). The exponent is ≤ t − 0.7j + C₂π^{0.7} + log ρ(·) + O(1), with
   j ≥ σ − 1 = t√c − 1 ≥ t(log t)/2 − 1.
   - Also π ≤ (n + 2)²/t (from σ ≤ n + 2) and t ≥ t₀ give π^{0.7}/σ ≤ c^{0.2}t^{−0.3} ≤ (n+2)^{0.4}t^{−0.7} → 0.
   - If j ≤ n/(log n)², the ρ-term is ρ(u₀)^{1−o(1)}, and t − (0.7 − o(1))j ≤ t(1 − 0.34 log t) < 0. This also
     absorbs C₂π^{0.7} = o(j).
   - If j > n/(log n)², then −0.7j alone is ≤ −F_{2/3}·ω(1).
8. Summing over poly(n) values of j gives the claim. ∎

**Weighted case (PROVED (sketch); the key inequalities are given).** Normalise max w = 1, so Q ≥ log 2. Apply step 5
of the unweighted proof first, so t ≥ t₀; the tail estimate needs this.
1. *Classes.* Split the primes of P into classes C_i = {q : w_q ∈ ((1+η)^{−i−1}, (1+η)^{−i}]}, i < I, plus a tail,
   with η = 1/log log n and I = O(log n/η).
2. *Tail.* Take ε = 1/log n and I = ⌈log((n + 2)/(ε√(t₀ log 2)))/log(1 + η)⌉ = O(log n·log log n). Then the tail
   contributes ≤ (1+η)^{−I}(n + 2) ≤ ε√(Q log V) to W. Here Q > 0, because some weight is positive.
3. *Cauchy–Schwarz.* Let P_i = ∏_{q∈C_i} q^{e_q}, s_i = log gcd(z, P_i) and π_i = log P_i. Then
   (Σ_i (1+η)^{−i}s_i)² ≤ (1+η)²Q·Σ_i s_i²/π_i. So certified ⇒ Σ_i s_i²/π_i ≥ t(1 − ε)²/(1 + η)² = t(1 − o(1)).
4. *Good and bad classes.* For each certified z, call class i *good* if r_i := (⌊s_i⌋ + 1)²/π_i ≥ R₀ := C·n^{4/7}·log n.
   - The bad classes carry r-mass ≤ I·R₀ = o(t), because t ≥ t₀ ≫ n^{4/7+o(1)}.
   - Condition only on the divisors of the good classes (dropping a class only enlarges the count). There are
     ≤ 2^I·(n+3)^I = exp(o(F_{2/3})) choices of good set and buckets.
5. *Per good class.* The unweighted Case A/B inequalities give −j_i + log N_i(j_i) ≤ −r_i(1 − o(1)). In Case B one
   needs π_i^{0.7} = o(s_i). With π_i = (⌊s_i⌋+1)²/r_i and s_i ≤ n + 2: π_i^{0.7}/s_i ≲ s_i^{0.4}/r_i^{0.7}
   ≤ (n+3)^{0.4}/R₀^{0.7} ≤ C^{−0.7}(log n)^{−0.7} → 0. The extra log n in R₀ is essential (codex): r_i ≥ C s_i^{4/7} alone
   gives only O(C^{−0.7}).
   - *Large Σj (proof critic).* In Case A, π_i = (j_i + 1)√c_i ≤ (j_i + 1)(log r_i)/2 and log π_i ≥ log r_i − O(1). So
     1.07π_i/log π_i ≤ 0.54(j_i + 1), and −j_i + log N_i ≤ −0.46j_i + O(1).
   - In Case B, −j_i + log N_i ≤ −0.69j_i.
   - So every good class contributes −Θ(j_i). If Σ_good j_i > n/(log n)², the extra negativity is ω(F_{2/3}).
6. *Total.* The exponent is ≤ t − t(1 − o(1)) + o(F_{2/3}) + log ρ(u₀(1 − o(1))). ∎

**Extensions (a), (b).**
- *(a) Average over N.* Fix k ≤ h. For z ∈ [kX, 2kX + h], #{N : 0 ≤ z − kN ≤ h} ≤ h/k + 1 ≤ 2h/k. So the averaged
  count is ≤ (2h/(kX))·Z_k, with Z_k = #{z ≤ 3kX : certified, smooth} ≤ Σ_{d|P, d≥T} Ψ(3kX/d, y) ≤ 3kX·C₀Σ_d d^{−1}ρ(…).
  The prefactor is 6C₀h, and steps 3–8 apply with V := h and size bound n + log k + O(1) in place of n + 2.
- *(b) Batched regions.* Sum (a) over the K values k₀ ≤ k < k₀ + K, all ≤ h. The prefactors add up to
  6C₀K(h + 1) = 6C₀|R|, matching the threshold √(Q log|R|). Steps 3–8 apply with t := log|R| and
  n + log(k₀ + K) + O(1). The batch parameters must not depend on N.
- Herrmann–May's linear bivariate result supplies a *certification region* of this form under their Assumption 1
  (pilist: g_n(μ) ≤ μ², PROVED there). It is a recovery guarantee, not a ceiling on multivariate constructions.

**Independent proof (codex blind round, §3(b)).** Codex proved the same bound for fixed P without the dyadic
case analysis.
- *Tilting.* Give each smooth m weight m^{−a}/Z_y(a), with a = 1 − log u₀/log y. Then X^{a−1}Z_y(a) = ρ(u₀)^{1−o(1)},
  and prime valuations are independent geometric variables.
- *Markov.* Apply Markov's inequality to gcd(m, P)^s with s ≈ log σ/log π − (1 − a).
- *Weights.* For GSS's coprime-modulus weighted certificate it uses the Kearns–Saul sub-Gaussian bound for
  Bernoulli variables (Berend–Kontorovich 2013; eq. (3) of arXiv 1210.3248).
- Two independent proofs of Theorem J (unweighted; and weighted with coprime moduli) therefore exist. Ours also
  covers prime-power scores J_q = min(v_q, e_q).

**Reading.**
- *Runtime corollary (scope stated after codex's review).* Consider an algorithm whose useful outputs are certified
  outputs of 𝒥-calls, and whose calls are fixed in advance or drawn from a repertoire of 2^b complete call types
  (P, w, region shape), fixed in advance. Over T calls its expected number of smooth outputs is
  ≤ T·2^b·ρ(u₀)^{1−o(1)}. So it needs b + log₂T ≥ (2/(3κ) − o(1))F_{2/3}/log 2.
  - With b = o(F_{2/3}) the Prop SI constructor in this sense costs ≥ L[2/3, 2/(3κ) − o(1)] calls per output, and the
    program costs L[2/3] per relation.
  - Equivalently (proof critic): the distinct outputs lie in the union of the repertoire's certified sets, so their
    expected number is ≤ 2^b·ρ(u₀)^{1−o(1)} whatever T is. One output needs a repertoire of size
    ≥ ρ(u₀)^{−1+o(1)}. A fixed schedule of T calls is the case 2^b = T.
  - This is a cardinality statement about fixed repertoires. It is not an information-theoretic or computational
    lower bound for choosing calls adaptively from x. The input has n ≫ F_{2/3} bits, so such a choice is not
    excluded, and we do not claim it is as hard as the original problem.
  - Uncertified roots of a real decoder, and arithmetic between calls, are outside the theorem.
- The o(1) in Theorem J hides a possible gain of exp(o(F_{2/3})). The unweighted proof gives O(F_{2/3}/log n); the
  weighted sketch with η = 1/log log n gives only O(F_{2/3}/log log n). That is larger than any L[1/3], so Theorem J
  does not decide whether the class gains L[1/3, Θ(1)]. The divisor trick gains L[1/3, 2ν/(3κ)] (siconstruct Lemma
  RC(iii)). Either way the cost stays L[2/3].
- *Scope: single-integer certificates.* Cohn–Heninger (ANTS X, 2013) decode m samples sharing one hidden divisor,
  with errors up to P^{β^{(m+1)/m}}, beyond the β² envelope. That is a joint certificate over several integers, for
  example in windows for m different k. It is not of the form W(z) ≥ T, so Theorem J does not cover it. The analogous
  count Σ_{d} ∏_i (V/d)·ρ (with V < d) suggests the same conclusion, but it is not proved here. Sarkar–Maitra (ePrint
  2010/146) revisit smooth integers in short intervals with Jochemsz–May lattices; we have not checked whether their
  guarantee lies in 𝒥.
- *Prior art.* Boneh's Theorem 3.2 bounds the number of N ∈ I with gcd(N, S) > T in the worst case. His remark that
  "a random interval I satisfying the bounds of Theorem 3.1 is unlikely to contain integers with sufficiently large
  smooth factors" is the qualitative content of Theorem J. Coppersmith–Howgrave-Graham–Nagaraj (Math. Comp. 2008)
  give lattice-derived counts of divisors in residue classes. We found no averaged statement of this kind; "no prior
  statement found by us".
- Compared with pilist's Prop JR (NFS scaling, congruence scaling excluded): here the regions are L[2/3]-large, and
  the quantity bounded is the absolute expected number of outputs, not the gain over the region. The proof replaces
  JR's random-divisibility model H_div by exact averaging over x, plus divisor counting (Nicolas–Robin and Rankin).

### 5.3 Model D (divisor completion) [PROVED; siconstruct Lemma RC]
- Constructions z = D·⌈x/D⌉ with D chosen obliviously (or z = D·m with m tested) succeed, averaged over x, with
  probability equal to the Ψ-density at x/D, i.e. (1 + o(1))ρ(u₀ − log D/log y).
- For D ≤ h this is at most ρ(u₀ − log h/log y) = ρ(u₀)·L[1/3, 2ν/(3κ)]. Equality holds for D = 2^{⌊log₂h⌋}, and
  small D gains less (D = 1 gains nothing). This is Theorem J's extremal case. siconstruct's E3 measured exactly this.

### 5.4 Model M2 (oblivious products, meet-in-the-middle, Wagner) [PROVED averaged; Wagner model-PROVED]

**Prop MM.** Let A, B be finite sets of positive integers fixed independently of x (resp. N). Assume the merge reads
both lists, so its cost is ≥ |A| + |B|.
- (i) *Interval.* E_{x∈[X,2X)}#{(a, b) : ab ∈ [x, x + h]} ≤ |A||B|(h + 1)/X. So one expected output needs
  max(|A|, |B|) ≥ √(X/(h + 1)).
  - *Per output* (corrected after the proof critic, who gave a counterexample to the draft's 2√(X/(h + 1))). A useful
    pair has ab ≤ 2X + h ≤ 3X, so min(a, b) ≤ √(3X). The useful pairs number ≤ √(3X)(|A| + |B|). Each lands in the
    window with averaged probability ≤ (h + 1)/X.
  - So expected outputs ≤ √(3X)(|A| + |B|)(h + 1)/X, and the cost per output is ≥ √(X/3)/(h + 1) = e^{(1/2 − o(1))n}.
    This is tight up to constants: take A = B = {1, …, ⌈√(3X)⌉}.
- (ii) *Bohr, all k at once.* E_{N∈[X,2X)}#{(a, b) : ab mod N ∈ [0, h], N ≤ ab ≤ KX} ≤ |A||B|·(h(1 + ln K) + K)/X.
  Charge each useful pair to its smaller factor, which is ≤ √(KX) since ab ≤ KX. The useful pairs then number
  ≤ √(KX)(|A| + |B|), so the outputs per unit cost are ≤ √(KX)(h(1 + ln K) + K)/X. Hence the cost per output is
  ≥ √(X/K)/(h(1 + ln K) + K).
  - For K = X^{o(1)} (Prop SI uses K = L[1/3, O(1)]) this is e^{(1/2 − o(1))n}.
  - When the K term dominates the bound is ≍ √X/K^{3/2}, constant already at K ≍ X^{1/3}. The draft's "unless
    K ≳ X/h²" was wrong (codex).
- *Proof.* (i) For a fixed product w, #{x : x ≤ w ≤ x + h} ≤ h + 1. (ii) For each pair, Pr_N[ab ∈ ∪_k [kN, kN + h]]
  ≤ Σ_{k≤K}(h/k + 1)/X. ∎
- siconstruct Lemma MA states (i)–(ii) under its equidistribution heuristic H_eq. The averaged form removes H_eq.
- E3 checks (i) exactly (ratio 1.000) and (ii) by Monte Carlo (measured ≈ 0.06 of the bound: the bound sums over all
  k, while products in [2^36, 2^37) meet only ln 2 worth of 1/k).

**Wagner / k-tree (HEURISTIC sketch; downgraded after codex's review).**
- *Simplified list model.* A node merges two lists and keeps the pairs whose partial log-sum (fixed point,
  n′ = log₂(X/h) bits) has a prescribed block of b bits.
  - The node's size is |L||R|/2^b.
  - A leaf list has at most as many entries as its value range.
  - There are no repeated representations: each value occurs once.
- *Bound.* Let the root's children A, B have value ranges with X_A X_B ≍ X. Then
  log₂|A| + log₂|B| ≤ log₂X − f_A − f_B, where f counts the bits filtered in the subtree. The root inherits at most
  min(path bits) pre-fixed bits and needs ≥ 1 expected survivor. Hence the path bits are ≤ log₂h + O(1), and
  max(|A|, |B|) ≥ 2^{(n′ − log₂h)/2} ≈ √X/h.
- *Caveats (codex).*
  - Lists of ordered factor pairs can contain about M log M representations of values ≤ M.
  - Deduplication changes the multiply-and-divide recurrence.
  - The "inherited bits" of unbalanced trees need a precise definition.
  So this is not a proof for Wagner's algorithm on integer lists. Codex's balanced-template analysis gives the weaker,
  but cleanly modelled, log T = Ω(√n). Both exceed L[1/3] by far.
- *Bohr variant.* The target {ab : ab mod N small} is not a coset of a subgroup with a known bit filtration, and none
  that can be computed efficiently is supplied. We do not claim one would require discrete logarithms.

### 5.5 Lattice reduction on Schnorr's log-lattice [estimates of known approaches, not a lower bound]
- siconstruct Prop CV (PROVED there, with its qualifications: squarefree witnesses, windows enlarged by O(√n),
  possible changes of multiplier and denominator, two-sided residuals): Prop SI is near-exact CVP_{1+O(1/n)} in
  dimension π(y) = L[1/3], up to those slacks.
- **Dimension lemma (PROVED, averaged, for a prime set fixed in advance).** If only the first d primes are allowed
  (or any fixed set of d primes, by the injection lemma), the expected number of admissible integers in [x, x + h] is
  ≤ (h + 1)·Ψ(3X, p_d)/X. This is ≥ 1 only if p_d ≥ L[1/3, 2/(3ν)(1 − o(1))]. At ν ≈ 2/(3κ) a fixed prime set must
  have ≥ L[1/3, κ(1 − o(1))] elements, essentially π(y). A prime subset selected from the input is not covered.
- *Estimates.* With the BKZ/GSA success heuristic of Alkim–Ducas–Pöppelmann–Schwabe (the "2016 estimate"), recovering
  a vector at distance ≈ λ₁/√2 needs blocksize β = Θ(d). Exact or sieving CVP costs 2^{Θ(d)}. In dimension
  d = L[1/3] this is exp(L[1/3]).
  - These are the costs of the generic approaches we know of. They do not prove that this structured family needs
    that time. The family's geometry (λ₁ an L[1/3] factor below the Gaussian heuristic, siconstruct Prop CV(e))
    makes generic estimates unreliable in either direction.
  - Background: Schnorr (EUROCRYPT 1991; DIMACS 1993), Adleman 1995, Ducas 2021 (experiments), schnorrcvp.

### 5.6 The collapse table

| Model | Per-output constructor cost | Prop SI program becomes | Status |
|---|---|---|---|
| point queries, independent labels (Lemma BB) | ≥ L[2/3, 2/(3κ)] queries | L[2/3] per relation | model-PROVED (surrogate model) |
| S (suffmap's sampling model: H_S, R, no prefiltering) | — | L[1/2, √2] overall | model-PROVED |
| J (single-integer Johnson certificates, call types from a fixed repertoire of size L[2/3]^{o(1)}) | ≥ L[2/3, 2/(3κ) − o(1)] calls | L[2/3] | PROVED (averaged) |
| D (divisor completion, D fixed in advance) | ≥ L[2/3, 2/(3κ)]·L[1/3, −2ν/(3κ)] (upper bound on success ⇒ lower bound on cost) | L[2/3] | PROVED (averaged) |
| M2 (fixed two-list product tables, all k at once, K = X^{o(1)}) | ≥ e^{(1/2−o(1))n} | exponential | PROVED (averaged) |
| Wagner k-tree (simplified list model) | ≥ √X/h | exponential | HEURISTIC sketch |
| generic lattice algorithms on the log-lattice | exp(L[1/3]) (estimate) | worse than trial division | estimate, not a bound |

None of these reaches L[1/3, anything], so within them "collapses back to ≥ 1.902" is an understatement: they never
reach exponent 1/3. The table says nothing about the mechanisms listed as not excluded in §0 and §7.

## 6. NP-hardness, TFNP, and where Prop SI's scale sits (Task 4)

Problem SMOOTH-INT(x, h, y) (all primes ≤ y allowed; exponents unbounded). Membership in NP: the witness is z with
its factorisation into primes ≤ y. No primality test is needed: it suffices that the factors are ≤ y and multiply
to z. In NP for every range.

**Prop NC (PROVED).**
- If h ≤ poly(log x), the decision problem is in NP ∩ coNP. The coNP witness gives, for each of the h + 1 integers,
  a prime factor > y with a primality certificate.
- So the search problem "output a smooth z or certify none" is in TFNP.
- *Conventions.* Megiddo–Papadimitriou (TCS 81, 1991, Thm 2.1) state the result for FNP reductions. It extends to
  Cook reductions that succeed for **every valid oracle answer**: an NP machine guesses the whole polynomial-length
  transcript of answers and verifies each one. Arbitrary promise search problems are not covered. With that
  convention, NP-hardness ⇒ NP = coNP.
- h = 0 is "is x y-smooth?". It is in NP ∩ coNP, and with y varying it is polynomial-time Turing equivalent to
  factoring: binary-search the largest prime factor, divide it out, repeat (codex).

**Prop TF (PROVED from KNOWN).** If a theorem guarantees that **every** interval [x, x + h(x)] with x ≥ x₀ contains
a y(x)-smooth integer, the *search* problem is total, so it is in TFNP and not NP-hard unless NP = coNP (same
conventions). The *decision* problem is then trivially YES and lies in P. Known
all-interval ranges (KNOWN, as summarised by Granville 2008 §4.1 and suffmap §3):
- Balog 1987: h = x^β, β > 1/2, y = x^α (any α > 0);
- Harman 1991: h = x^β, β > 1/2, y ≥ exp((log x)^{2/3+o(1)});
- Matomäki–Radziwiłł 2016, Cor. 1: h = C(ε)√x, y = x^ε, with at least √x(log x)^{−4} such integers (not positive
  local density);
- Jain (QJM 2026, Thm 1.2): an all-interval version with h ≥ √x·exp(…), y ≥ exp(C(log x)^{2/3}(log log x)^{4/3});
- Friedlander–Lagarias 1987: y = x^{1/2}, h = x^{1/4+ε}; y = x^{1/3}, h = x^{5/12+ε}.
- Lemma FP: h ≥ 2x/(y − 1) is in FP.
- Also in FP (codex; PROVED, elementary):
  - fixed y, by enumerating the O(log x)^{π(y)} exponent vectors;
  - h and y both poly(log x), by trial division of every integer in the window.

**Prop SI scale.**
- *Under Hyp E(κ, ν)*: every [x, x + L_x[2/3, ν]] with x ≥ x₀ contains an L_x[1/3, κ]-smooth integer. This is the
  existence half of Prop SI, far beyond the known results. Under it the problem is in TFNP, so NP-hardness ⇒
  NP = coNP.
- *Unconditionally*: the sieve decides and searches in time L_x[2/3, ν]·poly(n) = 2^{O(n^{2/3}(log n)^{1/3})}.
  Suppose a deterministic Karp or Cook reduction maps 3SAT instances with m variables to instances with
  n = O(m^c) bits. Then 3SAT runs in 2^{O(m^{2c/3}(log m)^{1/3})}, so ETH⁻ forces c ≥ 3/2, since c < 3/2 makes the exponent O(m^{1−δ}) (PROVED, routine; c is a
  constant). Randomized reductions need randomized ETH⁻. (v2 said ETH; see the erratum.)
- The NP-hardness of the real problem is OPEN in every range (to our knowledge). For the ranges just listed, and
  for h ≤ poly, it is implausible.

**Prop XS (explicit factor base; PROVED; parsimonious; `src/sihard/reduction.py`).** The problem S-SMOOTH-INT: given
a list S of primes and x, h, is some S-smooth z in [x, x + h]? It is NP-complete under randomized (zero-error)
polynomial-time reductions.
- **Prior art (literature critic, codex).** The core technique is known.
  - DIVISOR ("does N have a divisor in [L, U]?") is NP-hard under randomized reductions. The result is attributed to
    Alon–Kilian (Sudan, "The P vs. NP problem", 2010, credits Kilian).
  - Shor gave a proof (cstheory.SE, Q. 4769, 2011): it reduces from subset sum, chooses random primes p_i in short
    intervals so that log p_i encodes the weights, and sets N = ∏p_i.
  - A DIVISOR instance built that way is an S-SMOOTH-INT instance with exponents ≤ 1.
  - What Prop XS adds is modest: unbounded exponents, excluded by base-b digits with b > q and a fixed product count;
    parsimony, from X3C; and an explicit prime-supply parameter. We do not claim more.
- *Reduction from X3C* (exact cover by 3-sets, NP-complete; Karp 1972). Take universe {0, …, 3q−1} and sets
  T_1 … T_m.
  - b = q + 1, w_j = Σ_{e∈T_j} b^e, W = Σ_{e<3q} b^e, A_max = (b^{3q} + 1)R.
  - Q ≥ 2^{q+2}A_max²/R, R′ = R/(4q).
  - **Choice of R** (proof critic): R = K² with K = 2^{q+2}(b^{3q} + 1)². Then Q ≍ K³ and R′ ≍ Q^{2/3}/(4q). By
    Huxley's prime number theorem in short intervals (lengths ≥ Q^{7/12+ε}), [Q′, Q′ + R′) has density
    ~1/log Q of primes. So a random point is prime with probability ≫ 1/log Q, and sampling until AKS certifies a
    prime makes the step Las Vegas. All bit lengths are polynomial in q.
  - *Code.* `build_instance_theory` implements this choice with unbounded retries. It uses sympy's BPSW test in place
    of AKS, so the code is not itself a certified zero-error reduction. `build_instance` with R = 2^20 is the
    registered toy generator of E4.
  - p_j is a random prime in [Q + w_jR, Q + w_jR + R′).
  - E4 used the registered toy value R = 2^20 (prime supply heuristic there). E4b reruns with R = K², a disclosed
    post-review rerun (§8).
  - x = Q^q + Q^{q−1}RW, h = Q^{q−1}R/2.
- *Correctness.*
  - A product of r ≠ q primes from S lies outside [x, x + h]: (2Q)^{q−1} < x, and Q^{q+1} > x + h.
  - For r = q, the product is Q^q + Q^{q−1}(RΣw + Σε) + E with 0 ≤ Σε < qR′ ≤ R/4 and
    0 ≤ E ≤ Q^{q−1}·2^qA_max²/Q ≤ Q^{q−1}R/4.
  - So z ∈ [x, x + h] ⇔ Σ_{multiset} w = W.
  - Base b > q means no carries, so this holds ⇔ every digit is exactly 1 ⇔ the multiset is an exact cover with no
    repeated set.
  - z determines the multiset by unique factorisation. So the reduction is parsimonious. ∎
- *Parameters.* y ≈ Q, x ≈ y^q, h/x ≈ R/(2Q).
- *Counting.* #X3C is #P-complete under parsimonious reductions. This is Hunt–Marathe–Radhakrishnan–Stearns (SIAM J.
  Comput. 27, 1998), as cited by arXiv 2606.29157; we have not checked the original. Our reduction is randomized, so
  #S-SMOOTH-INT is #P-hard under randomized parsimonious reductions.
- *Codex's blind round* considered only the composite-generator route (Subset Product with generators ∏_{i∈C_j}p_i).
  It noted that splitting the generators into primes destroys the constraint, and that the exact-target version with
  an explicit prime set is polynomial. Both are correct. Prop XS avoids them: it encodes the weights in the *sizes* of
  the primes and uses an interval target, i.e. an approximate log-knapsack.
- *Why it does not transfer to all primes ≤ y.* At these parameters the full y-smooth count in the window is about
  h·ρ(q) ≫ 1. Indeed Lemma FP applies, since h ≥ 2x/(y − 1). Hardness needs the explicit sparse S to forbid the
  generic smooth integers. With all primes allowed, smooth integers are abundant exactly where the explicit
  instances live.
- E4 (R = 2^20) and E4b (R = K²) check the parsimony exactly: 300/300 instances each, q ≤ 4.

**Range table.**

| Range | Complexity status | Label |
|---|---|---|
| h ≥ 2x/(y−1); fixed y; h, y ≤ poly(log x) | FP | PROVED (Lemma FP; codex) |
| h = x^{β}, β > 1/2, y ≥ exp((log x)^{2/3+o(1)}); h = C√x, y = x^ε | search in TFNP (existence theorems), so NP-hard ⇒ NP = coNP; decision trivial | PROVED from KNOWN |
| h ≤ poly(log x), any y | NP ∩ coNP; NP-hard ⇒ NP = coNP | PROVED |
| Prop SI: h = L[2/3, ν], y = L[1/3, κ], ν > 2/(3κ) | in DTIME(L[2/3]); TFNP under Hyp E (which is Vinogradov-hard); ETH⁻ ⇒ deterministic reductions blow up ≥ m^{3/2} (was ETH; erratum) | PROVED (conditional parts as stated) |
| Boneh's range y = polylog x, h = x^{o(1)} | open; solutions rare (typically none in a random interval) | OPEN |
| explicit S (sparse prime list) | NP-complete (randomized reductions): DIVISOR (Alon–Kilian; Shor), extended in Prop XS | KNOWN / PROVED |

So **Prop SI's scale lies outside every range where hardness is provable, and inside the ranges where NP-hardness is
excluded conditionally.** The condition, Hyp E, is itself at least as hard as Vinogradov's least-nonresidue
conjecture for p ≡ 3 (mod 4) (Prop QNR). That is the opposite of floorhard's situation, where the hard instances sat in a different
region (x ≤ N^{3/8}) of the same problem.

## 7. Barrier statement (Task 5)

**Sharpest statement we can support.**
- **No barrier is known for the Prop SI constructor, and none is provable with current techniques**: a superpolynomial
  lower bound would imply P ≠ NP.
- **Model bounds.** In each of these models the constructor costs at least L[2/3, 2/(3κ) − o(1)] per output, or
  exponentially much. That is the leading cost of sieving, so the 1.2114 target collapses to L[2/3] per relation,
  and in suffmap's model S to L[1/2, √2] overall:
  - (i) independent-label point queries (a surrogate for sampling with free tests, factoring included);
  - (ii) single-integer Johnson-radius certificate decoders (Coppersmith, Howgrave-Graham, Boneh, GSS-weighted;
    univariate or Herrmann–May-batched over k ≤ h). Their certified outputs are counted, with call types from a
    repertoire fixed in advance of size L[2/3]^{o(1)} (Theorem J, PROVED, averaged);
  - (iii) divisor completion with the divisor fixed in advance;
  - (iv) fixed two-list product tables, including amortisation over all k ≤ K = X^{o(1)} (Prop MM, exponential).
  Wagner k-trees (a simplified list model) and generic lattice algorithms on the log-lattice (exp(L[1/3])) are
  estimates, not proofs.
- **Not excluded.** These are open mechanisms, with no positive evidence for any of them:
  - adaptive selection of moduli, lattices or prime subsets from the input;
  - correlated algebraic constructions in the style of Friedlander–Lagarias, where factors vary jointly, e.g.
    (m − a)(m + a) = m² − a²; scale invariance (siconstruct Lemma SR) does not address these;
  - multi-sample common-divisor decoders (Cohn–Heninger);
  - smooth values of polynomials or norm forms near (kN)^{1/d};
  - uncertified decoder outputs.
- **NP-hardness cannot be the barrier.**
  - At Prop SI's scale it is excluded under the existence hypothesis (TFNP), and limited by ETH⁻ (was ETH; erratum) for deterministic
    reductions.
  - The only provably hard variant has an explicit sparse prime list (DIVISOR / Prop XS).
- **The strongest evidence of implausibility is Prop DLI.**
  - A modulus-generic constructor would also give prime-field DL at L[1/3, 1.232] (1.442 with plain ECM booting).
  - For general primes, prime-field DL has stood at 1.902 since Matyukhin 2003.
- **The existence half is itself deep.** All-interval existence at Prop SI's scale implies Vinogradov's
  least-nonresidue conjecture for p ≡ 3 (mod 4) (Prop QNR). So only promise constructors are realistic targets, and
  the promise cannot be proved with current methods.

**Plausibility (assessment).** About 1% for an L[1/3, τ < 0.757] constructor, uniform or not. We put a
modulus-generic one lower, because of Prop DLI. Codex (blind): about 5% for a uniform promise constructor.

## 8. Experiments (pre-registered: commit 7f9e398 and addendum A1 bb7fbd0)

Commands: `.venv/bin/python -m src.sihard.run_experiments e1|e2|e3|e4|analyze`. Total runtime < 3 min; RSS < 300 MB.

**E1. Smooth counts near kN versus random positions** (48-bit N, 10 seeds, 200 intervals each).

| y | μ | h | mean near | mean control | ratio | ID near / ctrl | empty near / ctrl | ctrl/Dickman | P1 | P2 | P3 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2^10 | 2 | ≈1.8e4 | 1.877 | 2.176 | **0.863** | 1.51 / 1.50 | 0.190 / 0.147 | 1.09 | **FAIL** | pass | pass |
| 2^10 | 8 | ≈7.3e4 | 8.26 | 8.63 | 0.957 | **3.77 / 3.36** | 0 / 0.001 | 1.08 | pass | **FAIL** | pass |
| 2^10 | 32 | ≈2.9e5 | 34.49 | 34.79 | 0.991 | **12.5 / 11.0** | 0 / 0 | 1.09 | pass | **FAIL** | — |
| 2^12 | 2 | ≈1.5e3 | 1.760 | 1.943 | 0.906 | 1.06 / 1.16 | 0.163 / 0.149 | 0.97 | pass | pass | pass |
| 2^12 | 8 | ≈6.0e3 | 7.32 | 7.87 | 0.930 | **2.06 / 1.96** | 0.002 / 0.001 | 0.98 | pass | **FAIL** | pass |
| 2^12 | 32 | ≈2.4e4 | 30.71 | 31.75 | 0.967 | **5.78 / 5.02** | 0 / 0 | 0.99 | pass | **FAIL** | — |

- **P1 failed in 1 cell, P2 in 4 cells, P3 passed.**
- *Diagnosis of P2 (POST-HOC, `src/sihard/diag_e1.py` part 1).* The registered design pools k = 1…200,
  across which the density varies about fivefold, together with the control's matching size range. That inflates the
  dispersion index for both samples alike: the ID near and control agree to within 0.15–1.5, and both grow with μ.
  Within bins of 100 consecutive k (1500 k, two seeds) the ID is 1.00–1.02 near kN and 0.99–1.00 for size-matched
  random offsets. So the counts are Poisson-like once size is controlled. The registered P2 band was a design error.
- *Diagnosis of P1 (POST-HOC, `src/sihard/diag_e1.py` part 2).* With 3000 multipliers and size-matched controls (x = kN + random
  offset) at y = 2^12, h = 1509, the near/control ratio varies by modulus: 0.97 and 0.96 (two semiprimes, z = −1.0
  and −1.5), 1.11, 0.93 and 0.87 (three random odd moduli, z = +3.6, −2.6, −4.9), and 1.01 (an even modulus).
  - So intervals placed along one progression kN *appear* to carry a modulus-dependent bias of either sign, about
    ±10% at this size.
  - Caveat (codex): the z-scores use independent-Poisson standard errors. Counts along one fixed progression are
    deterministic and arithmetically correlated, so "beyond noise" is not established. Neither is the
    binned-ID "Poisson-like" reading, which rests on two moduli.
  - The registered pooled deficit (all six cells below 1) is consistent with that spread plus chance.
  - A constant-factor effect of this kind would not affect L-scale heuristics. HYPOTHESIS, not established: "kN is a
    random position" holds only on average over N at toy sizes, perhaps through the alignment of kN modulo large
    smooth divisors d > h.
  - A proper test would compare common shifts of whole progressions, over many more moduli, with block-wise
    uncertainty estimates. Not done.
- The control density matches the Dickman prediction within 0.97–1.09 (exploratory).

**E2. Johnson certificates, exact averages** (X = 2^40; E2: y = 2^8, M = 2^16; E2b: y = 2^10, M = 2^20).
- **P4 (A ≤ B ≤ divisor-sum form): PASS** in all 84 cells.
- E2 has only 25–33 smooth integers in range, so A is mostly 0, as anticipated in A1.
- E2b (exploratory): the certified smooth yield per call, as A/dens, i.e. in units of one random test:

| V | divisor trick | best of family | best P | random squarefree P (c = 4) |
|---|---|---|---|---|
| 2^8 | 5.5 | 7.1 | lcm, c = 2 | 0.17–0.28 |
| 2^10 | 7.5 | 11.3 | lcm, c = 2 | 0–0.23 |
| 2^12 | 13.6 | 15.4 | primorial, c = 2 | 0 |

- At toy size a call is worth 5–15 random tests: this is the L[1/3, 2ν/(3κ)]-type gain of the divisor trick.
- The best P beats the divisor trick by ≤ 1.5×. Random squarefree P (large primes) are far worse.
- Consistent with Theorem J's sharpness at leading order. The toy cannot see asymptotics.

**E3. Oblivious two-list merging** (A, B: 2000 random 2^8-smooth integers in [2^18, 2^18√2)).
- Interval (exact): measured/bound = 1.000 at h = 2^6 and 2^10 (**P5, P6 PASS**).
- Bohr (4096 N): measured/bound = 0.065 and 0.059 (**P5 PASS**).

**E4. X3C reduction.** 300 instances (q ≤ 4; 166 with covers; max 2 covers). #S-smooth in window = #exact covers in
**300/300 (P7 PASS)**.
- **E4b** (disclosed post-review rerun, `run_experiments e4b`): the same design with the provable choice R = K² (x up to
  741 bits). 300/300 again (154 with covers, at most 3 covers).

**Deviations.** Addendum A1 (E2b) was added before any run. E4b was added after the proof critic asked for a
provable prime-supply parameter. It is a rerun of E4's design with one parameter changed, labelled POST-HOC. The E1 diagnostics are post hoc and unregistered; their
scripts are reproduced in `src/sihard/diag_e1.py`.

## 9. Relation to the parallel agents
- **siconstruct** (constructions). Its scale invariance (Lemma SR), random cofactor (Lemma RC), amortised two-list
  bound (Lemma MA, under H_eq) and Prop CV (CVP) are the construction-side counterparts of §5.
  - Theorem J makes its row 4 (Boneh/Coppersmith, "MODEL") rigorous on average, and extends it to weighted and
    batched decoders.
  - Prop MM removes H_eq from Lemma MA by averaging.
  - Its §4.5 remark on DL ("not checked in detail") is Prop DLI.
- **sicheck**: no results at the time of writing.

## 10. GPT-6.1-sol rounds
### 10.1 Blind round (`codex/sihard_1.out`; run before our analysis was shown to it)

**Independent agreements.**
- Lemma EQ.
- The DL descent through another modulus (M = ⌈jp/q⌉, k = 1, no descent tree), with constants 0.987 (ECM boot),
  1.4422 and 1.052 (early-abort boot).
- The class-group and generic-group negatives (codex cites Shoup 1997 for the latter).
- The point-query model (its "permuted-candidate model"; = Lemma BB).
- **Theorem J for fixed P**, unweighted and GSS-weighted, with a different proof (§5.2).
- The smooth-multiple trick as the comparison point.
- Fixed product tables need A + B ≳ √(X/h) (= Prop MM(i)).
- Amortisation over L[1/3] outputs subtracts only O(F_{1/3}).
- NP membership, NP ∩ coNP for h ≤ poly, the sieve's 2^{O(n^{2/3}(log n)^{1/3})}, the ETH⁻ exponent 3/2 (was ETH; erratum).
- Proven-existence ranges are TFNP (Megiddo–Papadimitriou Thm 2.1).
- SI itself is not TFNP without Hyp E.
- τ < 0.757.

**Adopted from codex.**
- The early-abort booting constant 1.232 (FGHT 2017) and hence DL ≤ 1.232.
- Hypothesis H-affine and the normalisation remark in Prop DLI.
- **Prop QNR.**
- Remark FO.
- The FP cases (fixed y; polynomial h and y).
- The extra existence citations (Croot 2007, Matomäki 2016, Younis 2024, Jain 2026) and Naccache–Shparlinski 2009.
- The explicit caveat that a decoder may find useful *uncertified* outputs.

**Ours, beyond codex.**
- Prop XS. Codex's blind round found no hardness for explicit prime sets; its adversarial round and the literature
  critic then pointed to the known DIVISOR hardness, which Prop XS only extends.
- Theorem J for prime-power scores, I_k- and batched (Herrmann–May-type) regions, and families of moduli.
- The Wagner bound √X/h. Codex's balanced template gives only log T = Ω(√n); ours uses the fact that leaf lists
  cannot exceed their value range.
- The Bohr-set two-list bound.
- The lattice dimension lemma.
- Lemma FP.
- The E1–E4 experiments.

**Differences.**
- Codex labels the adaptive extension "model-PROVED in a fresh-divisibility model". We state only the union bound (c).
- Codex's plausibility is 5%; ours is 1%.
### 10.2 Adversarial round (`codex/sihard_2.out`): MAJOR REVISION, all points adjudicated
Codex reviewed draft v1 and revision ded0b41, and ran the tests (17 pass). It confirmed:
- the unweighted proof of Theorem J (after small endpoint repairs), the weighted bookkeeping, and the batched
  extension once it is summed explicitly;
- Prop DLI's algebra and constants, including FGHT's 1.232;
- Prop NC, Prop XS's inequalities, no-carry argument and parsimony;
- the E1 measurements.

| # | Sev. | Finding | Adjudication |
|---|---|---|---|
| 1 | MAJOR | Theorem J promoted to an algorithmic lower bound beyond its scope | accepted: runtime corollary restricted to certified outputs and fixed repertoires; broader reading marked unproved (§5.2 Reading, §0, §7) |
| 2 | MAJOR | "bits of information about x" does not follow from the union bound; selection ≠ hardness | accepted: replaced by "log of the size of a fixed repertoire of complete call types"; the equivalence claim deleted |
| 3 | MINOR | d ≤ 3X, clamped ρ, small P, uniform Ψ citation | accepted (step 2) |
| 4 | MINOR | weighted: the log n in R₀ is essential; η loses F_{2/3}/log log n | accepted: calculation displayed; remainder claim weakened to exp(o(F_{2/3})); tail cutoff explicit; Q > 0 |
| 5 | MINOR | batched extension needs its own summation; n + log k; HM as certification region | accepted (Extensions (a), (b)) |
| 6 | MINOR | sharpness / Model D overstate the gain for arbitrary d ≤ h | accepted: d = 2^{⌊log₂h⌋}; Model D as an upper bound |
| 7 | MINOR | Boneh's completeness needs T < U | accepted (§4) |
| 8 | MINOR | DL normalisation and linear algebra over Z/(p−1) | accepted (§2.2, linear-algebra bullet) |
| 9 | MINOR | booting constants are heuristic upper bounds; "general primes" | accepted |
| 10 | MAJOR | Lemma BB is not a factoring-oracle separation; Prop S's assumptions hidden | accepted: renamed "independent-label point-query model", limits stated; Prop S lists H_S, R and no prefiltering |
| 11 | MINOR | MM(ii)'s transition "K ≳ X/h²" is wrong; charge pairs to the smaller factor | accepted |
| 12 | MAJOR | Wagner bit budget lacks a premise (representation multiplicities, deduplication) | accepted: downgraded to a HEURISTIC sketch in a stated simplified model; DL remark softened |
| 13 | MAJOR | the lattice row is not a lower bound; CVP correspondence has slacks; dimension lemma only for fixed prime sets | accepted (§5.5 rewritten; collapse table relabelled "estimate") |
| 14 | MINOR | Megiddo–Papadimitriou conventions for Cook reductions; decision vs search; randomized ETH | accepted (§6) |
| 15 | MAJOR | the code is not a ZPP reduction (default R, retry cap, BPSW) | accepted: `build_instance_theory` (R = K², unbounded retries); BPSW vs AKS stated; "bit lengths polynomial in q" |
| 16 | MINOR | Alon–Kilian / Shor DIVISOR antecedent; randomized counting reductions | accepted (also found by the literature critic) |
| 17 | MINOR | Friedlander–Granville range, Younis condition, Jain Thm 1.2, MR count ≠ positive density | accepted |
| 18 | MAJOR | scale invariance does not exclude correlated algebraic constructions (Friedlander–Lagarias) | accepted: "shared reason" claims replaced by a list of unresolved mechanisms (§0, §7) |
| 19 | MINOR | E1's statistical and causal conclusions too strong | accepted: labelled hypotheses; Poisson-SE caveat |
| 20 | MINOR | `johnson_S` not exact; `lemma_fp` rejects y = 2 | fixed (mpmath, labelled numerical; y = 2 and x ≤ y branches). E2 rerun: identical |

## 11. Opus critics

### 11.1 Proof critic: ACCEPT WITH FIXES, all applied
The critic checked Theorem J step by step (averaging, divisor decomposition, the uniformity of the Hildebrand Ψ
bound, Nicolas–Robin, injection and Rankin, Case A convexity, Case B, uniformity of the o(1)). It confirmed that
Boneh's Theorem 3.1 and Herrmann–May lie in 𝒥. It verified Prop DLI's algebra and constants (3^{1/3}, 0.987 < 1.442,
1.052 < 1.232), Prop QNR, FO, FP, EQ, NC, TF, the ETH argument, the Wagner bit budget, the dimension lemma, and every
inequality of Prop XS against the code. 17 tests pass.

| # | Sev. | Finding | Fix |
|---|---|---|---|
| 13 | MAJOR | Prop MM(i) "cost per output ≥ 2√(X/(h+1))" is false (counterexample A = B = {1..⌈√(3X)⌉}) | replaced by the useful-pair bound √(X/3)/(h+1), tight up to constants; conclusion unchanged |
| 1 | MINOR | Case B "j ≥ 1.4t ⇒ t − 0.7j ≤ 0" is wrong (0.98 < 1) | uses j ≥ t(log t)/2 − 1 |
| 2 | MINOR | weighted large-Σj case covered only Case-B classes | Case A gives −0.46j_i; added |
| 4 | MINOR | sharpness example with d ∈ [V/2, V] certifies nothing | V := d |
| 5 | MINOR | (b) needs every k in the batch ≤ h | k₀ + K − 1 ≤ h |
| 9 | MINOR | 1.232 rested on an unchecked citation | now checked in FGHT App. A.1 (Barbulescu's analysis), with the misprint noted |
| 11 | MINOR | the deterministic x^ε-in-√x consequence needs a non-promise constructor (Hyp E) | qualified in §0, §2.4 |
| 18 | MINOR | R unspecified in Prop XS; prime supply heuristic | R = K², Huxley; AKS; E4b rerun |
| 19–21 | MINOR | §0/§7 overclaims: "no better than the trivial trick", lattice row unlabelled, "only a mechanism that selects" | "at leading L[2/3] order"; certified-only caveat; HEURISTIC labels; the selection statement PROVED only inside 𝒥 |
| 3, 6–8, 10, 12, 14–17, 22–24 | NIT | weighted label "PROVED (sketch)"; t ≥ t₀ for the tail; many-window remark; primes > y WLOG; H-affine is about 𝒪's choices; Lemma BB (T+1)δ; Remark FO wording; MM(ii) constants unified; ETH "c ≥ 3/2"; NP membership needs no AKS; Hyp E ≥ Vinogradov; Boneh nits; test count | applied |

### 11.1b Proof critic, confirmation round on v2 (7a19fac): ACCEPT
- Every changed step was re-checked: Theorem J steps 2 and 7, weighted steps 2 and 5, Extensions (a)/(b), the
  sharpness, the runtime corollary, Prop MM(i)/(ii), and Prop XS with `build_instance_theory`.
- Five NITs, all applied: the π_i definition in weighted step 5; the explicit constant in (log n)^{−0.7}; the MM(ii)
  constant in §0; the K = X^{o(1)} qualifier in §7; a stale 1.442 phrase. The critic also suggested the
  repertoire-size form of the runtime corollary, now added. The reduction.py docstring was aligned.

### 11.2 Literature critic: OK WITH FIXES, all applied
The critic read Boneh, Granville, FGHT, GSS, Barbulescu–Pierrot, Matomäki–Radziwiłł, Jain, Croot, Naccache–Shparlinski,
Berend–Kontorovich, Nevo, and the cstheory.SE thread.

| # | Sev. | Finding | Fix |
|---|---|---|---|
| 1 | MAJOR | Prop XS's core (random primes whose sizes encode weights; interval target) is the known DIVISOR hardness (Alon–Kilian via Sudan 2010; Shor, cstheory.SE 4769) | Prop XS reframed as an extension (unbounded exponents, parsimony); §0, §4, §6 and the table rewritten |
| 2 | MAJOR | Theorem J's "any decoder" claim misses multi-sample Cohn–Heninger decoders (beyond β²) | scope narrowed to single-integer certificates; Cohn–Heninger and Sarkar–Maitra listed as not covered |
| 3 | MINOR | Friedlander–Granville range is √x·y², not √(xy) | fixed |
| 4 | MINOR | attributions: 1.902 to Matyukhin (via Barbulescu–Pierrot), 3^{1/3} to Commeine–Semaev | fixed |
| 5 | MINOR | Prop DLI forerunners: Granville's t_j = [jp/q] device; Schnorr EUROCRYPT 1991 | cited |
| 6 | MINOR | Theorem J forerunners: Boneh Thm 3.2 and his random-interval remark; CHN 2008 | cited; Theorem J presented as their averaged quantitative form |
| 7 | MINOR | Naccache–Shparlinski confirmed; signature-forgery line of work relevant | UNVERIFIED tag removed; CNS 1999, CNTW 2009 cited |
| 8 | MINOR | #X3C reference: Hunt–Marathe–Radhakrishnan–Stearns 1998 (via a citing paper) | cited as such |
| 9 | MINOR | lattice section lacks primary citations | Schnorr, Adleman, Ducas, ADPS 2016 added |
| 10–15 | NIT | Granville wording (≪, "≥ π(y)"); Boneh base-2 logs and the QS quote's context; Jain's title and doi; Matomäki 2010; Younis's y condition; MR count; GJ SP14 credited to Yao; Harvey 2021; Berend–Kontorovich numbering; gcd threshold renamed T | applied |

Confirmed by the critic:
- FGHT's 1.232 and 0.811, the misprint, and that the boot does not use the NFS polynomials;
- Granville §1.2/§4.4 (Prop QNR attribution) and §4.1 (Balog, Harman, Friedlander–Lagarias, Challenge Problem 2000);
- Boneh's Theorem 3.1 details;
- GSS Theorem 4;
- MR Corollary 1.

No prior statement was found for Theorem J's averaged bound, for Prop DLI, for the 1.2114 constant, or for any
NP-hardness of the all-primes problem.

## 12. Reproduction
```bash
.venv/bin/python -m pytest -q tests/test_sihard.py
.venv/bin/python -m src.sihard.run_experiments e1      # ~2 min
.venv/bin/python -m src.sihard.run_experiments e2      # < 1 s
.venv/bin/python -m src.sihard.run_experiments e3      # ~40 s
.venv/bin/python -m src.sihard.run_experiments e4      # < 1 s
.venv/bin/python -m src.sihard.run_experiments analyze
.venv/bin/python -m src.sihard.diag_e1                 # POST-HOC diagnostics (~3 min)
.venv/bin/python -m src.sihard.consequences            # constants
```
