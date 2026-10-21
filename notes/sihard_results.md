# sihard: the hardness side of Prop SI (smooth numbers in short intervals near kN)

Agent `sihard`, phase 20, 2026-10-21. DRAFT v1 (before the codex adversarial round and the Opus critics).

**Files.**
- Code: `src/sihard/`.
  - `core.py`: exact uint64 smoothness sieve, Dickman ρ (design only), Johnson-certificate counts, two-list merge
    counts.
  - `reduction.py`: the X3C → S-smooth-in-interval reduction (Prop XS).
  - `consequences.py`: Prop SI and Prop DLI constants, the collapse table.
  - `run_experiments.py`: E1–E4 and the analysis.
- Tests: `tests/test_sihard.py` (16 pass).
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
  superpolynomial lower bound would separate P from NP. But **every natural way of building the constructor fails,
  and fails completely.** In each model below the constructor costs at least L[2/3, 2/(3κ) − o(1)] per output. That
  is the same leading cost as sieving the interval, so the 1.2114 target collapses to L[2/3] per relation.
  - In the sampling/sieving model the whole relation program then costs L[1/2, √2].
  - Meet-in-the-middle costs exponentially much, and lattice reduction costs exp(L[1/3]).
  - Nothing ends up near the 1.902 of MNFS, let alone below it.
- **The cleanest new result is Theorem J (PROVED, averaged over the interval position).** Take any Coppersmith /
  Howgrave-Graham / Boneh / Guruswami–Sahai–Sudan CRT list decoder:
  - any modulus P built from factor-base prime powers, any weights, any lattice dimension;
  - univariate, or bivariate batched over many multipliers k.
  Averaged over where the window sits, one call returns at most L[2/3, −2/(3κ) + o(1)] smooth certified integers.
  That is no more than testing L[2/3]^{o(1)} random integers. The trivial "divisor trick" z = d·⌈x/d⌉ attains this
  bound at leading order, so the whole Johnson-radius class is no better than the trivial trick.
  - The bound needs no heuristic: averaging over x replaces the random-divisibility assumption. It holds whenever
    the modulus is chosen without looking at x (or from a family fixed in advance of size L[2/3]^{o(1)}).
  - It does not cover moduli selected from x. Making that selection is the problem itself.
- **Consequences beyond factoring (Prop DLI).** A constructor that works for every modulus would also give
  prime-field discrete logarithms in L_p[1/3, ≤ 3^{1/3} = 1.442] (heuristic). Collecting relations costs 1.2114.
  The constructor also handles the descent of medium-sized primes, at L[1/3, 0.99] when queried with the modulus
  ⌊jp/q⌋. The classical ECM "booting" step, at 3^{1/3}, then dominates. Prime-field DL has been stuck at 1.902, so a
  modulus-generic constructor would break **two** problems believed hard.
  - No implication was found for class groups or for other groups.
  - It would also give the first subexponential-time constructions of x^ε-smooth integers in [x, x + √x]. That is
    the constructive counterpart of the Matomäki–Radziwiłł existence theorem; none is known to us.
- **The converse fails in the black-box model.** A factoring oracle does not help: in that model finding a smooth
  integer in a window needs ≥ 1/δ smoothness tests, where δ is the local density (Lemma BB). The constructor's
  difficulty is *selection*, not *testing*.
- **NP-hardness does not reach Prop SI's scale** (§6).
  - "Smooth over an explicitly given prime list S, in an interval" is NP-complete under randomized reductions
    (Prop XS; parsimonious, from exact cover by 3-sets).
  - For the real problem (all primes ≤ y) no NP-hardness is known in any range.
  - It is provably not NP-hard unless NP = coNP (via NP ∩ coNP or TFNP) in three cases:
    - when h ≤ poly(log x);
    - in every range where all intervals provably contain a smooth integer;
    - at Prop SI's own scale, if the existence half of Prop SI holds.
  - Unconditionally, the interval sieve solves Prop SI's scale in L[2/3]. So a reduction from 3SAT would need
    output size ≥ m^{3/2−o(1)} unless ETH fails.
- **Assessment.** Probability that an L[1/3, τ] constructor with τ < 0.757 exists: about 1% (ours). The models fail
  for one shared reason, siconstruct's scale invariance: every factor of u faces the same relative precision
  h/x ≈ N^{−1+o(1)}. Only a mechanism that selects from x with ≈ (2/(3κ))F_{2/3} bits of x-dependent information
  could escape.

**Main results.**

| Result | Status |
|---|---|
| Lemma EQ: a modulus-uniform Prop SI constructor = an arbitrary-interval constructor at the same scale | PROVED (trivial) |
| Prop DLI: modulus-generic constructor ⇒ F_p DL in L_p[1/3, max(1.2114, c_boot)], c_boot ≤ 3^{1/3}; descent 0.99 | HEURISTIC (H_S on residues, rank mod ℓ \| p−1) |
| Lemma BB: black-box/factoring-oracle model needs ≥ (1−o(1))/δ tests per output | model-PROVED |
| Theorem J: Johnson-class decoders (weighted, any P, batched) average ≤ L[2/3, −2/(3κ)+o(1)] smooth certified outputs per call; sharp | PROVED (averaged over x or N; oblivious modulus) |
| Prop S: Prop SI-type relations in sampling model S cost ≥ L[1/2, √2] | model-PROVED (suffmap Theorem SS(a)) |
| Prop MM: oblivious two-list merging, interval and Bohr (all k at once): ≥ √(X/(h+1)) per output, resp. ≥ √(N/K)/(2(h+1)(1+ln K)) | PROVED (averaged); siconstruct Lemma MA under H_eq |
| Wagner/k-tree bit budget: cost ≥ √X/h | model-PROVED (Wagner's filter semantics) |
| Prop NC / TF / ETH: not NP-hard (h ≤ poly, proven-existence ranges, Prop SI scale under Hyp E); ETH blow-up ≥ m^{3/2} | PROVED (from KNOWN Megiddo–Papadimitriou) |
| Prop XS: explicit-S version NP-complete (ZPP reductions), parsimonious | PROVED; E4 300/300 |
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
  - Sieving [x, x + z] costs ≍ z log log y + u·y (Granville 2008, §5.7; KNOWN, checked in the text). One output
    needs z ≈ 1/ρ(u₀) = L[2/3, 2/(3κ)].
  - Granville notes that the y-term is a likely lower bound when the primes ≤ y are part of the input, and that
    "one might guess that it is not necessary to have z in the running time". Prop SI asks exactly for that. The
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

**Claim.** Under H_S (residues are smooth like random integers of their size), the usual index-calculus rank
heuristic modulo every ℓ^e ∥ p − 1, and ECM's heuristic running time, discrete logarithms in F_p^* cost
L_p[1/3, max((16/9)^{1/3} + ε, c_ind)], where c_ind ≤ 3^{1/3} = 1.4422. So the total is ≤ L_p[1/3, 1.4422].

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
   - The cost min_c 1/(3c) + √(4c/3) = 3^{1/3} is at c = 3^{−1/3} (KNOWN classical booting constant;
     `consequences.dl_booting_ecm`).
   - At that c the descent costs (0.693 + 1.1006)/(3·0.6057) = 0.987 < 1.442. So booting dominates:
     c_ind = 3^{1/3}.
4. Smaller booting constants from early-abort strategies would apply unchanged. Values near 1.23 are reported for
   NFS-DL booting (Barbulescu's thesis; UNVERIFIED). With them c_DL = max(1.2114, c_boot).
5. *Can 𝒪 also boot?* Booting needs t·s ≡ w (mod p) with s, w both smooth. The interval oracle offers only additive
   control, u ∈ [t + jp, t + jp + h], so it gives u ≡ t + v, not a multiplicative relation. Multiplicative control
   would be a Bohr-set (coset) lister, suffmap E2 with Prop DLT. We found no way to make 𝒪 boot.

**Remarks.**
- This sharpens suffmap §4's unchecked remark: the proof critic's "≈ L[1/3, 0.83]" is the descent without the ECM
  cost of booting (`consequences.dl_descent_optimum`: 0.827).
- The I_k-only constructor also suffices, *provided* it is uniform over moduli: the descent queries modulus
  ⌊jp/q⌋, not p. It does **not** need suffmap's coset hypothesis.
- Codex's blind round obtained the same constants independently (boot 1.4422, descent 0.987; §10.1).
- *Reading.* Prime-field DL has resisted any improvement on L[1/3, 1.902] (Matyukhin 2003; Commeine–Semaev 2006). A
  modulus-generic Prop SI constructor would bring it to ≤ 1.442 and factoring to 1.2114. That is evidence against
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
- *Consequence.* It would give, in time x^{o(1)} (subexponential in log x), x^ε-smooth integers in [x, x + C√x]. That
  is Matomäki–Radziwiłł's existence range, and also Harman's and Balog's ranges [x, x + x^{1/2+ε}].
  - The constructive results known to us are much weaker:
    - sieving (exponential);
    - Boneh (polylog y, §3);
    - Friedlander–Lagarias 1987, as reported by Granville §4.1: x^{1/2}-smooth integers in every interval of
      length x^{1/4+ε}, by explicit construction (UNVERIFIED in the primary source);
    - the trivial Lemma FP below.
  - We found no subexponential construction in the MR range.
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
| F_p discrete log | yes | ≤ L[1/3, 1.442] (booting); relations 1.2114; descent 0.99 | HEURISTIC (Prop DLI) |
| DL mod composite N | yes | via Bach 1984 | HEURISTIC |
| smooth-in-interval at larger (y, h) | yes, deterministic if 𝒪 is | e.g. x^ε-smooth in [x, x+√x] in time x^{o(1)} | PROVED (monotone) |
| class groups (imaginary quadratic), units | no implication found | — | assessment |
| EC / genus-g / small-char DL | no | — | assessment |
| deterministic primes | no | — | assessment |
| rigorous factoring | no (residue randomness is heuristic) | — | assessment |

## 3. The converse: is Prop SI implied by factoring, or by anything known? (Task 1b)

**Lemma BB (model-PROVED).** Model: the algorithm learns whether z is y-smooth, or its full factorisation, only by
querying z, and it may do arbitrary other computation. The smooth set is replaced by a random set in which each z
in [x, x + h] is included independently with probability ≤ δ. Then T queries find a member of [x, x + h] with
probability ≤ Tδ.
- *Proof.* Union bound over the queried points; unqueried points are independent of the transcript. ∎
- *Reading.* With δ = ρ(u₀)(1 + o(1)) (Hildebrand), one output needs ≥ L[2/3, 2/(3κ)]/2 queries for success
  probability 1/2. A factoring oracle is such a query oracle, so it changes the cost of a test, not the number of
  tests.
- The interval sieve already tests h integers in h·O(log log y) operations. So even free factoring improves
  the cost per examined integer by at most O(log log y) (compare polyimpl Lemma PF, where a remainder tree replaces
  factoring).
- The model ignores arithmetic: real algorithms can *construct* smooth numbers (products). §5 treats those.

**Known constructive results do not reach Prop SI's scale.**
- Lemma FP needs h ≥ x/y.
- Friedlander–Lagarias type identities reach h = x^{1/4+ε}, still exponential.
- Boneh's decoder lists, and it lists nothing in random intervals at this scale (§4, Theorem J).
- Existence theorems (Hildebrand; Friedlander–Granville; Balog; Harman; Croot; Matomäki; Matomäki–Radziwiłł;
  Jain 2025, arXiv 2502.10530) are non-constructive, and none reaches y = L_x[1/3] in all intervals of length
  L_x[2/3]. Granville §4.1: (4.1) is known for exp((log x)^{5/6+o(1)}) ≤ y and z ≥ √(xy)·exp((log x)^{1/6}).
- **Conclusion.** Prop SI is implied by nothing known. Polynomial-time factoring is not known to imply it, and in
  the black-box model it provably does not. Conversely Prop SI gives factoring only at L[1/3, 1.21], not
  polynomial time. The two are incomparable as far as is known.

## 4. Literature status (Task 2)

- **Sieving.** Granville 2008 §5.7 (MSRI 44, checked): running time ≍ z log log y + uy for [x, x + z].
- **Boneh 2000/2002** (STOC 2000; JCSS 64 (2002) 768–784; Theorem 3.1 checked in the author's PostScript).
  - Setting: S = ∏_{q<s} q^{⌊log s/log q⌋}. Every *strongly* s-smooth integer (p^m ∤ N whenever p^m > s) divides S.
  - For any d ≥ 1 + √(log S/log|I|), it finds in poly(d, size) time all N ∈ I with gcd(N, S) > T, whenever
    T > S^ε with ε = √(log(4|I|)/log S) + 5/(4d). Equivalently |I| < ¼·T^{(log T/log S) − 2.5/d}.
  - If moreover V < 2T, it outputs all strongly s-smooth integers in I.
  - Granville's summary: for (1/5)(log x)² > y > log x, there are ≪ y/log x integers in [x, x + x^{(log x)/(4y)}]
    dividing lcm[1..y]. Boneh: "a random interval I satisfying the bounds of Theorem 3.1 is unlikely to contain
    integers with sufficiently large smooth factors". §4 says the bounds "are insufficient for improving the basic
    quadratic sieve".
  - **At Prop SI's scale** (log x ≈ n, log S ≈ y = L[1/3]): log|I| ≲ n²/y < 0, so it is vacuous for full
    smoothness. Partial smoothness is Theorem J: no gain beyond L[2/3]^{o(1)}.
  - Boneh also dismisses the subset-sum route: "the dimension of the lattice … is exponential in the size of the
    input (log s)".
- **Coppersmith / Howgrave-Graham / GSS weighted CRT decoding.** The same Johnson-type guarantee
  (log S)² ≥ log V·log P, weighted: Σw·agreement ≥ √(Q log V) (GSS FOCS 2000). Herrmann–May bivariate linear:
  ∏X_i < P^{g(μ)} with g(μ) ≤ μ² (pilist §3.1). All are in Theorem J's class.
- **Schnorr's log-lattice.** siconstruct Prop CV (PROVED there): a squarefree Prop SI output is a CVP_{1+O(1/n)}
  answer in the π(y) = L[1/3]-dimensional prime lattice, and conversely.
- **Hardness literature.** We found no NP-hardness result for finding smooth integers in intervals, and no statement
  that the question is open in complexity terms. Granville §5.7 frames finding smooth numbers without paying the
  interval length as open beyond Boneh. Related NP-complete problems (KNOWN):
  - Subset Product (Garey–Johnson 1979, SP14, by reduction from X3C; full proof in Moret's textbook, per the
    literature search; UNVERIFIED by us);
  - Manders–Adleman 1978 (bounded quadratic congruences; used in floorhard Theorem SS);
  - subset sum.
  None of them has the all-primes-≤-y structure.
- **Existence.** See §3. Granville's "Challenge problem 2000" (x^α-smooth in [x, x + √x] for all α > 0) is
  answered up to a constant factor in the interval length by Matomäki–Radziwiłł 2016, Corollary 1 (as
  cross-checked by suffmap's literature critic).

## 5. Model lower bounds (Task 3)

Scale throughout: x ∈ [X, 2X), y = L_X[1/3, κ], window length V ≤ L_X[2/3, ν_max]. "Per output" means expected
cost per y-smooth integer found in the target window(s).

### 5.1 Model S (sampling / sieving) [model-PROVED]
- **Lemma BB** (§3): ≥ (1 − o(1))/ρ(u₀) = L[2/3, 2/(3κ)] tests per output.
- **Prop S.** In suffmap's model S, a Prop SI-type relation has smoothness integers u ≍ kN and v ≤ h. So
  λ = log u + log v ≥ n, and Theorem SS(a) gives total cost ≥ exp((√2 − o(1))√(n log n)) = L[1/2, √2]. This is
  siconstruct row 0.
- Amortisation over k is irrelevant here: each examined integer is smooth with the same probability.

### 5.2 Model J (Johnson-radius certificate decoders): Theorem J [PROVED, averaged]

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

**Theorem J.** Fix κ, ν_max > 0. Uniformly over all oblivious (P, w) and all V ≤ L_X[2/3, ν_max]:

  (1/X)·Σ_{x=X}^{2X−1} #{z ∈ C([x, x + V]) : z y-smooth} ≤ exp(−(2/(3κ) − o(1))·F_{2/3}) = ρ(u₀)^{1−o(1)}.

- (a) The same holds for the I_k-shaped regions [kN, kN + h] averaged over N ∈ [X, 2X) (any fixed k ≤ h).
- (b) It also holds for batched regions R_N = {kN + t : k₀ ≤ k < k₀ + K, 0 ≤ t ≤ h}, K ≤ h, averaged over N, with
  V := K(h + 1) in the threshold. This is the Herrmann–May-type bivariate guarantee.
- (c) For an algorithm choosing P from a family of size |𝒫|, multiply the bound by |𝒫|.
- (d) **Sharpness.** The divisor trick, P = d ∈ [V/2, V] y-smooth with w = 1 (so certified ⇔ d | z), attains
  ρ(u₀)·L[1/3, 2ν/(3κ)] = ρ(u₀)^{1−o(1)}.

**Proof (unweighted).** Write t = log V, π = log P and σ = √(tπ). Certified means log gcd(z, P) ≥ σ.
1. *Averaging (exact).* Each z lies in at most V + 1 of the regions [x, x + V]. So the left side is
   ≤ ((V+1)/X)·Z, where Z = #{z ∈ [X, 2X + V] : log gcd(z, P) ≥ σ, z y-smooth}.
2. *Divisor decomposition.* Write z = d·m with d = gcd(z, P). Then m is y-smooth, so
   Z ≤ Σ_{d|P, log d ≥ σ} Ψ(3X/d, y). The sum is empty unless σ ≤ n + 2, and nonempty only if π ≥ σ, i.e. π ≥ t. Put
   c := π/t ≥ 1, so σ = t√c.
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
   ≤ t(1 − √c + 1.07c/log t) + O(1). As a function of a = √c ∈ [1, (log t)/2] the bracket is concave, so it is at
   most its larger endpoint value: 1.07/log t at a = 1, and 1 − 0.23 log t < 0 at the other end. So the factor is
   ≤ exp(1.07t/log t + O(1)) = exp(o(F_{2/3})).
   - Every d here has log d ≤ π ≤ t(log t)² = o(n/log n), so ρ((n − j − 2)/log y) ≤ ρ(u₀(1 − 1/log n)) = ρ(u₀)^{1−o(1)}.
   - The number of j is ≤ π + 1 = poly(n).
7. *Case B: c > (log t)²/4.* Use (b). The exponent is ≤ t − 0.7j + C₂π^{0.7} + log ρ(·) + O(1), with
   j ≥ σ − 1 = t√c − 1 ≥ 1.4t.
   - Also π ≤ (n + 2)²/t (from σ ≤ n + 2) and t ≥ t₀ give π^{0.7}/σ ≤ c^{0.2}t^{−0.3} ≤ (n+2)^{0.4}t^{−0.7} → 0.
   - If j ≤ n/(log n)², the ρ-term is ρ(u₀)^{1−o(1)} and t − 0.7j ≤ 0.
   - If j > n/(log n)², then −0.7j alone is ≤ −F_{2/3}·ω(1).
8. Summing over poly(n) values of j gives the claim. ∎

**Weighted case (PROVED; sketch with the key inequalities).** Normalise max w = 1, so Q ≥ log 2.
1. *Classes.* Split the primes of P into classes C_i = {q : w_q ∈ ((1+η)^{−i−1}, (1+η)^{−i}]}, i < I, plus a tail,
   with η = 1/log log n and I = O(log n/η).
2. *Tail.* The tail contributes ≤ (1+η)^{−I}(n + 2) ≤ ε√(Q log V) to W.
3. *Cauchy–Schwarz.* Let P_i = ∏_{q∈C_i} q^{e_q}, s_i = log gcd(z, P_i) and π_i = log P_i. Then
   (Σ_i (1+η)^{−i}s_i)² ≤ (1+η)²Q·Σ_i s_i²/π_i. So certified ⇒ Σ_i s_i²/π_i ≥ t(1 − ε)²/(1 + η)² = t(1 − o(1)).
4. *Good and bad classes.* For each certified z, call class i *good* if r_i := (⌊s_i⌋ + 1)²/π_i ≥ R₀ := C·n^{4/7}·log n.
   - The bad classes carry r-mass ≤ I·R₀ = o(t), because t ≥ t₀ ≫ n^{4/7+o(1)}.
   - Condition only on the divisors of the good classes (dropping a class only enlarges the count). There are
     ≤ 2^I·(n+3)^I = exp(o(F_{2/3})) choices of good set and buckets.
5. *Per good class.* The unweighted Case A/B inequalities give −j_i + log N_i(j_i) ≤ −r_i(1 − o(1)). In Case B one
   needs π_i^{0.7} = o(s_i), which follows from r_i ≥ C s_i^{4/7}. If Σ_good j_i > n/(log n)², the extra negativity
   is ≫ F_{2/3}, exactly as in step 7.
6. *Total.* The exponent is ≤ t − t(1 − o(1)) + o(F_{2/3}) + log ρ(u₀(1 − o(1))). ∎

**Extensions (a), (b).**
- Average over N instead of x. For z ∈ [kX, 2kX + h], #{N : 0 ≤ z − kN ≤ h} ≤ h/k + 1 ≤ 2h/k. Step 1 becomes
  LHS ≤ (2h/(kX))·#{z ≤ 3kX : certified, smooth}, and steps 2–7 are unchanged with V replaced by h or K(h + 1).
- Herrmann–May's guarantee for linear bivariate polynomials (under their Assumption 1) satisfies the same inequality
  with V = |R| (pilist: g_n(μ) ≤ μ², PROVED there).

**Reading.**
- One Johnson-class call is worth ρ(u₀)^{1−o(1)} outputs, i.e. L[2/3]^{o(1)} random tests. So the Prop SI
  constructor in 𝒥 costs ≥ L[2/3, 2/(3κ) − o(1)] calls per output, and the program costs L[2/3] per relation.
- The o(1) in Theorem J hides a possible gain exp(O(F_{2/3}/log n)). That is larger than any L[1/3], so Theorem J
  does not decide whether the class gains L[1/3, Θ(1)]. The divisor trick gains exactly L[1/3, 2ν/(3κ)] (siconstruct
  Lemma RC(iii)). Either way the cost stays L[2/3]: the collapse is complete.
- Compared with pilist's Prop JR, which covers NFS scaling and excludes the congruence scaling: here the regions are
  L[2/3]-large, and the quantity bounded is the absolute expected number of outputs, not the gain over the region.
  The proof replaces JR's random-divisibility model H_div by exact averaging over x, plus divisor counting (Wigert
  and Rankin).
- *What escapes.* A modulus or weight vector chosen *from x* with ≥ (2/(3κ) − o(1))F_{2/3}/log 2 bits of
  x-dependence (by (c)). Making such a choice is itself the selection problem.

### 5.3 Model D (divisor completion) [PROVED; siconstruct Lemma RC]
- Constructions z = D·⌈x/D⌉ with D chosen obliviously (or z = D·m with m tested) succeed, averaged over x, with
  probability Ψ-density at x/D, i.e. ρ(u₀ − log D/log y)(1 + o(1)).
- For D ≤ h this is ρ(u₀)·L[1/3, 2ν/(3κ)]. This is Theorem J's extremal case. siconstruct's E3 measured exactly this.

### 5.4 Model M2 (oblivious products, meet-in-the-middle, Wagner) [PROVED averaged; Wagner model-PROVED]

**Prop MM.** Let A, B be finite sets of positive integers fixed independently of x (resp. N). Assume the merge reads
both lists, so its cost is ≥ |A| + |B|.
- (i) *Interval.* E_{x∈[X,2X)}#{(a, b) : ab ∈ [x, x + h]} ≤ |A||B|(h + 1)/X. So one expected output needs
  max(|A|, |B|) ≥ √(X/(h + 1)), and the cost per output is ≥ 2√(X/(h + 1)) = e^{(1/2 − o(1))n}.
- (ii) *Bohr, all k at once.* E_{N∈[X,2X)}#{(a, b) : ab mod N ∈ [0, h], N ≤ ab ≤ KX} ≤ |A||B|·(h(1 + ln K) + K)/X.
  Taking the side with elements ≤ √(2KX + h) as the counted factor, outputs per unit cost ≤ √(2KX)(h(1 + ln K) + K)/X.
  So the cost per output is ≥ √(X/(2K))/(h(1 + ln K) + K): exponential unless K ≳ X/h².
- *Proof.* (i) For a fixed product w, #{x : x ≤ w ≤ x + h} ≤ h + 1. (ii) For each pair, Pr_N[ab ∈ ∪_k [kN, kN + h]]
  ≤ Σ_{k≤K}(h/k + 1)/X. ∎
- siconstruct Lemma MA states (i)–(ii) under its equidistribution heuristic H_eq. The averaged form removes H_eq.
- E3 checks (i) exactly (ratio 1.000) and (ii) by Monte Carlo (measured ≈ 0.06 of the bound: the bound sums over all
  k, while products in [2^36, 2^37) meet only ln 2 worth of 1/k).

**Wagner / k-tree (model-PROVED in Wagner's filter semantics).**
- In Wagner's semantics a node merges two lists and keeps the pairs whose partial log-sum (fixed point, n′ = log₂(X/h)
  bits of precision) has a prescribed block of b bits. List sizes multiply and then divide by 2^b.
- Let the root's children A, B have value ranges X_A X_B ≍ X. Then log₂|A| + log₂|B| ≤ log₂X − f_A − f_B, where f
  counts the bits filtered inside the subtree, which is at least the bits filtered along any root-to-leaf path. The
  root must filter n′ minus the path bits, and must have ≥ 1 expected survivor.
- Hence the path bits are ≤ log₂h + O(1), and max(|A|, |B|) ≥ 2^{(n′ − log₂h)/2} ≈ √X/h.
- So filtering below the root buys at most a factor h; this is smoothenum Prop W extended to unbalanced trees.
- *Bohr variant.* The target {ab : ab mod N small} is not a coset of a subgroup with a bit filtration. Wagner's group
  structure would need discrete logarithms mod N. So the k-tree does not apply at all.

### 5.5 Model L (lattice reduction on Schnorr's log-lattice) [HEURISTIC, cost model]
- siconstruct Prop CV (PROVED there): Prop SI ⇔ CVP_{1+O(1/n)} in dimension π(y) = L[1/3].
- **Dimension lemma (PROVED, averaged).** If only the first d primes are allowed, the expected number of admissible
  integers in [x, x + h] is ≤ (h + 1)·Ψ(3X, p_d)/X. This is ≥ 1 only if p_d ≥ L[1/3, 2/(3ν)(1 − o(1))]. At ν ≈ 2/(3κ)
  the dimension must be ≥ L[1/3, κ(1 − o(1))], essentially π(y).
- In the BKZ/GSA model (2016 estimate), recovering a vector at distance ≈ λ₁/√2 needs blocksize β = Θ(d), so the
  cost is 2^{Θ(d)} = exp(L[1/3]). That is worse than trial division. siconstruct row 1 agrees.

### 5.6 The collapse table

| Model | Per-output constructor cost | Prop SI program becomes | vs MNFS 1.902 |
|---|---|---|---|
| S (sampling/sieving, black box) | L[2/3, 2/(3κ)] | L[2/3] per relation; best use in S: L[1/2, √2] overall | collapses |
| J (Johnson class, oblivious modulus) | ≥ L[2/3, 2/(3κ) − o(1)] calls | L[2/3] | collapses |
| D (divisor completion) | L[2/3, 2/(3κ)]·L[1/3, −2ν/(3κ)] | L[2/3] | collapses |
| M2 (oblivious MITM, all k at once) | ≥ e^{(1/2−o(1))n} | exponential | collapses |
| Wagner k-tree | ≥ √X/h | exponential | collapses |
| L (BKZ on the log-lattice) | exp(L[1/3]) | worse than trial division | collapses |

In none of these models does the constructor cost even L[1/3, anything]. "Collapses back to ≥ 1.902" is
an understatement: they never reach exponent 1/3.

## 6. NP-hardness, TFNP, and where Prop SI's scale sits (Task 4)

Problem SMOOTH-INT(x, h, y) (all primes ≤ y allowed; exponents unbounded). Membership in NP: the witness is z with
its factorisation, checked with AKS. In NP for every range.

**Prop NC (PROVED).**
- If h ≤ poly(log x), the decision problem is in NP ∩ coNP. The coNP witness gives, for each of the h + 1 integers,
  a prime factor > y with a primality certificate.
- So the search problem "output a smooth z or certify none" is in TFNP. If it were NP-hard under Cook reductions,
  then NP = coNP (Megiddo–Papadimitriou 1991, TCS 81).
- h = 0 is "is x y-smooth?", which is in NP ∩ coNP and close to factoring.

**Prop TF (PROVED from KNOWN).** If a theorem guarantees that **every** interval [x, x + h(x)] with x ≥ x₀ contains
a y(x)-smooth integer, the search problem is total, so it is in TFNP and not NP-hard unless NP = coNP. Known
all-interval ranges (KNOWN, as summarised by Granville 2008 §4.1 and suffmap §3):
- Balog 1987: h = x^β, β > 1/2, y = x^α (any α > 0);
- Harman 1991: h = x^β, β > 1/2, y ≥ exp((log x)^{2/3+o(1)});
- Matomäki–Radziwiłł 2016, Cor. 1: h = C(ε)√x, y = x^ε;
- Friedlander–Lagarias 1987: y = x^{1/2}, h = x^{1/4+ε}; y = x^{1/3}, h = x^{5/12+ε}.
- Lemma FP: h ≥ 2x/(y − 1) is in FP.

**Prop SI scale.**
- *Under Hyp E(κ, ν)*: every [x, x + L_x[2/3, ν]] with x ≥ x₀ contains an L_x[1/3, κ]-smooth integer. This is the
  existence half of Prop SI, far beyond the known results. Under it the problem is in TFNP, so NP-hardness ⇒
  NP = coNP.
- *Unconditionally*: the sieve decides and searches in time L_x[2/3, ν]·poly(n) = 2^{O(n^{2/3}(log n)^{1/3})}.
  Suppose a Karp or Cook reduction maps 3SAT instances with m variables to instances with n = O(m^c) bits. Then 3SAT
  runs in 2^{O(m^{2c/3}(log m)^{1/3})}, so ETH forces c ≥ 3/2 − o(1) (PROVED, routine).
- The NP-hardness of the real problem is OPEN in every range (to our knowledge). For the ranges just listed, and
  for h ≤ poly, it is implausible.

**Prop XS (explicit factor base; PROVED; parsimonious; `src/sihard/reduction.py`).** The problem S-SMOOTH-INT: given
a list S of primes and x, h, is some S-smooth z in [x, x + h]? It is NP-complete under randomized (ZPP)
polynomial-time reductions.
- *Reduction from X3C* (exact cover by 3-sets, NP-complete; Karp 1972). Take universe {0, …, 3q−1} and sets
  T_1 … T_m.
  - b = q + 1, w_j = Σ_{e∈T_j} b^e, W = Σ_{e<3q} b^e, A_max = (b^{3q} + 1)R.
  - Q ≥ 2^{q+2}A_max²/R, R′ = R/(4q).
  - p_j is a random prime in [Q + w_jR, Q + w_jR + R′); expected O(log Q) trials, AKS or Miller–Rabin.
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
- *Counting.* #S-SMOOTH-INT is #P-hard if #X3C is #P-complete under parsimonious reductions (believed standard;
  reference UNVERIFIED).
- *Why it does not transfer to all primes ≤ y.* At these parameters the full y-smooth count in the window is about
  h·ρ(q) ≫ 1. Indeed Lemma FP applies, since h ≥ 2x/(y − 1). Hardness needs the explicit sparse S to forbid the
  generic smooth integers. With all primes allowed, smooth integers are abundant exactly where the explicit
  instances live.
- E4 checks the parsimony exactly: 300/300 instances, q ≤ 4.

**Range table.**

| Range | Complexity status | Label |
|---|---|---|
| h ≥ 2x/(y−1) | FP | PROVED (Lemma FP) |
| h = x^{β}, β > 1/2, y ≥ exp((log x)^{2/3+o(1)}); h = C√x, y = x^ε | TFNP (existence theorems); NP-hard ⇒ NP = coNP | PROVED from KNOWN |
| h ≤ poly(log x), any y | NP ∩ coNP; NP-hard ⇒ NP = coNP | PROVED |
| Prop SI: h = L[2/3, ν], y = L[1/3, κ], ν > 2/(3κ) | in DTIME(L[2/3]); TFNP under Hyp E; ETH ⇒ any reduction blows up ≥ m^{3/2} | PROVED (conditional parts as stated) |
| Boneh's range y = polylog x, h = x^{o(1)} | open; solutions rare (typically none in a random interval) | OPEN |
| explicit S (sparse prime list), any shape | NP-complete (ZPP reductions), parsimonious | PROVED |

So **Prop SI's scale lies outside every range where hardness is provable, and inside the ranges where NP-hardness is
excluded conditionally.** That is the opposite of floorhard's situation, where the hard instances sat in a different
region (x ≤ N^{3/8}) of the same problem.

## 7. Barrier statement (Task 5)

**Sharpest statement we can support.**
- In five natural models the Prop SI constructor costs at least L[2/3, 2/(3κ) − o(1)] per output, or more:
  - sampling/sieving with black-box tests;
  - Johnson-radius CRT/Coppersmith certificate decoders with oblivious moduli and weights (Theorem J: any P, any
    dimension, univariate or batched over k);
  - divisor completion;
  - oblivious meet-in-the-middle and Wagner k-trees (exponential);
  - lattice reduction on the log-lattice (exp(L[1/3])).
- That is the leading cost of sieving. So the 1.2114 target collapses completely in each model: L[2/3] per relation,
  and L[1/2, √2] for the best relation program in model S. That is far above 1.902.
- **Outside these models no barrier is known, and none is provable with current techniques**: a superpolynomial
  lower bound would imply P ≠ NP.
- **NP-hardness cannot be the barrier.** At Prop SI's scale it is excluded under the existence hypothesis (TFNP) and
  constrained by ETH. The only provably hard variant has an explicit sparse prime list.
- **The strongest evidence of implausibility is Prop DLI.** A modulus-generic constructor would also give prime-field
  discrete logs at ≤ L[1/3, 1.442], while that problem has resisted any improvement on 1.902.
- **Where an escape would have to live.** Every model fails for one shared reason. The target window has relative
  width h/x = N^{−1+o(1)}, and every factor of u inherits that precision (siconstruct Lemma SR). Oblivious objects (a
  modulus, lists, lattice data not chosen from x) hit it only at the random rate. An escape must select, using
  ≈ (2/(3κ))F_{2/3} bits of information about x, which is exactly the selection problem.

**Plausibility (assessment).** About 1% for an L[1/3, τ < 0.757] constructor, uniform or not. We put a
modulus-generic one lower, because of Prop DLI.

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
  - So intervals placed along one progression kN carry a modulus-dependent bias of either sign, about ±10% at this
    size, beyond Poisson noise.
  - The registered pooled deficit (all six cells below 1) is consistent with that spread plus chance.
  - A constant-factor effect of this kind does not affect L-scale heuristics. But it shows that "kN is a random
    position" (H_S) holds only on average over N at toy sizes. One suspected source is the alignment of kN modulo
    large smooth divisors d > h. Not investigated further.
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

**Deviations.** Addendum A1 (E2b) was added before any run. The E1 diagnostics are post hoc and unregistered; their
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
### 10.1 Blind round (`codex/sihard_1.out`)
(pending)
### 10.2 Adversarial round (`codex/sihard_2.out`)
(pending)

## 11. Opus critics
(pending)

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
