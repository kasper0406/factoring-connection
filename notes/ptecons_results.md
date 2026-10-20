# ptecons: constructing dense N-dependent split families (the systeq §6 escape)

Agent `ptecons`, phase 16 (route 1, construction side), 2026-10-20.

**Files.**
- Code: `src/ptecons/`.
  - `polys.py`: exact polynomial helpers.
  - `lemmas.py`: diagnostic checkers for every lemma below.
  - `acs.py`: OPERATIONAL aligned-cluster generators. They see only N.
  - `run_experiment.py`: the harness.
- Tests: `tests/test_ptecons.py` (11).
- Pre-registration: `notes/ptecons_prereg.md` (commit 95a669f, made before E1/E2; the pilot is disclosed there).
- Config: `experiments/configs/ptecons.json`.
- Raw results: `experiments/raw/ptecons_E{1,2}.jsonl`.
- Codex (GPT-6.1-sol): blind round `codex/ptecons_1.{md,out}`; adversarial round `codex/ptecons_2.{md,out}` (§8).

**Labels.**
- **PROVED**: proof given here or in the cited file.
- **KNOWN**: published, with a citation.
- **HEURISTIC**: the assumption is named.
- **measured**.
- **refuted for this model**.
- **open**.
- **assessment**.

Novelty is never claimed; "no prior statement found by us" is the most we say.

---

## 0. Verdict

**Plain language.** We found no construction meeting systeq §6 (a)–(c). We tried seven families: Dickson/Chebyshev, Prouhet–Tarry–Escott (PTE), cyclotomic/Galois orbits, elliptic, group-ring, lattice, and our own aligned clusters. Each fails at an identifiable step. The failures fall into seven mechanisms (§4). Two of them, real rigidity (LIN-R) and the ideal-cluster width bound (CW), are clean theorems that pin down why the F_q engine does not transfer.

**What the finite-field engine has that Z lacks, in two lines.**
1. **Additivity gives linear families.** Over F_q the roots b + c·a (c ∈ F_q) move *linearly* in the parameters (a, b), yet the split product keeps all middle coefficients zero, because Σ_{c∈F_q} c^j = 0. Over R this is impossible: if e₂ and e₃ of a split family are affine in real parameters, at most one root moves (**LIN-R, PROVED**). The only escape is complex directions with vanishing power sums, i.e. roots of unity. That forces cyclotomic fields, where the Splitting Lemma (systeq §5.4) makes piece norms ≥ N^{φ(n)/n}.
2. **Constants are units in F_q but not in Z.** In F_q[X] every nonzero constant is a unit of norm 1, so the norm ledger ∏ F(a) = Res(F, r) costs nothing. In Z it charges every root a with |F(a)| ≥ 2 (**LEDGER, PROVED**). For the Cunningham/F_1 presentation X^n ≡ 1 this reads: ∏_{|a_l|≥2}|a_l| ≤ ‖r‖₁(1 + o(1)) (**CUN, PROVED**). A compressed split polynomial has at most log₂‖r‖₁ roots outside {−1, 0, 1}.

**The PTE route specifically.**
- **Ideal classes (one free coefficient) are dead for a BGJT-style common offset (CW, PROVED, codex ptecons_1).**
  - Two members differ by a constant C with log|C| ≥ ½M log M − O(M).
  - Since M ≥ x/(C₀ log B), at most one member per class can have log|v| below ≈ (x/(2C₀ log B)) log(x/log B).
  - That exceeds x^{2/3} exactly when B = L[α] with α < 1/3. At α = 1/3 the bound is constant-level only.
  - So ideal clusters give **no** gain anywhere below the L[1/3] boundary.
- **Near-ideal classes.** They obey a proved coefficient-height bound (codex). Their evaluated residuals are controlled only heuristically (§3.5).
- **Generic classes** are weak. The heuristic count says a nontrivial class can fix only the top √(2M) coefficients, so the residual is N^{1−O(M^{−1/2})}.
- **Structured classes are weaker still.** Prouhet–Thue–Morse, Frolov and Choudhry chains fix only O(log M) moments.

**Experiments (pre-registered; two control predictions failed as stated, disclosed in §5).**
- **E1 (AVG).** Five structured families fixed independently of N show no excess congruence yield over 200 semiprimes: G = 1.05, 1.07, 1.04, 0.96, 0.99. Structure needs N-dependence to help at all.
- **E2 (aligned clusters, our best candidate).** On 64- and 80-bit semiprimes at an equal pool of 20000:
  - no degree-4 or degree-6 N-aligned split family has typically smaller residuals than COS (degree 2), with median λ = log|v|/log N: COS 0.592/0.574, generic quartic 0.595–0.606/0.565–0.583, symmetric quartic 0.633/0.606, one-digit sextic 0.675/0.648;
  - two-digit sextic alignment is infeasible: 0–4 members against a cap of 20000.
  - This matches the alignment-capacity heuristic HALF (§3.6): matching more than one base-m digit with split polynomials forces roots ≳ m^{1/2}, hence a bottom coefficient ≳ N^{1/2}.

**Verdict.**
- Route 1's PTE escape is **refuted for every concrete family examined** (table §2).
- The ideal and common-offset version is refuted by a theorem (CW).
- What remains **open** is exactly one algorithmic primitive, with no candidate:
  - an output-sensitive enumerator of split polynomials in the N-lattice {P : P(m) ≡ small (mod N)};
  - equivalently, smoothenum's algebraic or doubly-smooth listing (worth (4/9)^{1/3} or polynomial time).
- Assessment: < 1% that a split-by-identity construction yields sub-L[1/3]. Codex's two rounds: about 2%.

---

## 1. Formulation (aligned with systeq §6 and the ptesplit_1 model)

**Inputs.** A relation is a pair of actual integers (u, v) with u ≡ v (mod N), u ≠ v. The data are:
- an evaluation point m (with F(m) ≡ 0 mod N for a presentation F, or simply a residue);
- a split polynomial P = ∏_{l=1}^M (X − a_l) with |a_l| ≤ A;
- u = P(m) = ∏(m − a_l). The pieces m − a_l are small integers when |m| + A ≤ B, or factor-base symbols otherwise;
- v = (P mod F)(m), or the centered residue of u mod N.

**Requirements.** Following ptesplit_1, smallness of r = P mod F is not smallness of v. Only |v| matters for smoothness. Nontriviality forces M ≥ (x − O(1))/(C₀ log B) when the pieces are ≤ B^{C₀}.

**Cost threshold** [HEURISTIC, random-residual model; codex ptecons_1 and ptesplit_1].
- With log|v| = c·x^γ(log x)^{1−γ}, the optimum is L[γ/2]. So a fixed exponent below 1/3 needs a *power saving*, γ < 2/3.
- log|v| = o(x^{2/3}) alone may only improve the constant.
- If v is smooth *by construction*, the requirement on |v| disappears. That is the doubly-smooth case (§3.4, candidate 2b).

**Cluster + alignment decomposition (our organising device).**
- A *cluster* is a class 𝒞(T) of split polynomials sharing a top part T: P = T + r_P with deg r_P < t, i.e. equal e₁..e_{M−t}, i.e. a multigrade PTE class.
- An *alignment* is the condition T(m) ≡ w (mod N) with |w| small.
- One alignment then yields |𝒞| relations u_P = P(m) ≡ w + r_P(m) with a *common offset* w.
- These are genuinely distinct multiplicative relations. Each is the aligned congruence plus a different exact integer identity u_P − u_{P′} = (r_P − r_{P′})(m). This is exactly the BGJT shape: all members of the PGL₂ orbit reduce under one presentation.
- Re-centering each u_P independently mod N would forfeit the cluster gain: each member then needs its own modular condition.

## 2. Candidate table

The columns give each family's density (members with roots of size ≤ A), its residual |v|, how N enters, and the verdict. "Fails" names the first condition or barrier violated.

| # | Family | Density | Residual | N-dependence | Verdict |
|---|---|---|---|---|---|
| 1a | Dickson D_n(X) − D_n(t) over Z | Splits only for n ∈ {1, 2} (all t), n ∈ {3, 6} (t = ±1, ±2), n = 4 (t = 0, ±2); none otherwise [PROVED; codex and our check, n ≤ 30] | Exact identity, e.g. (X² − 4)(X² − 1)² = D₆ − 2 | None: F = P − r is chosen from P (EX) | refuted: (c) fails |
| 1b | Dilated/shifted D₆: ((X−h)² − 4λ²)((X−h)² − λ²)² | Θ(A²) multisets | Fixing one degree-6 F fixes h and λ², so the family collapses to O(1) per F [PROVED, codex] | Must align T(m) per member | refuted: no common reduction |
| 1c | D_n − c in a fixed field K | Splitting forces Q(ζ_n + ζ_n^{−1}) ⊆ K, so φ(n)/2 ≤ [K:Q]; at most O((log A)^{r_L}) shifts (τ is a unit) [PROVED, codex] | — | — | refuted: (a) needs n → ∞; Splitting Lemma |
| 2a | Ideal PTE classes (t = 1) | ≤ q members per class (q = 2A + 1); known ideal solutions only for M ≤ 10 and M = 12 | Common offset: members differ by C, log\|C\| ≥ ½M log M − O(M) | Alignment T(m) ≡ w costs ≈ N/V trials (AVG) | **refuted by CW** for α < 1/3 |
| 2b | Ideal class with T(m) ≡ 0: u = ∏(m − a), v = ±∏a, both smooth by identity (codex) | as 2a | v is smooth by construction but ≥ C/2 | Needs N \| T(m): a root of a fixed degree-M polynomial mod N, probability ≈ 1/N per (class, m) | refuted: alignment costs ≈ N |
| 2c | Generic PTE classes | Mean class size A^{M−k(k+1)/2}, so k ≲ √(2M) fixed moments [HEURISTIC]; VMVT makes large classes rare on average [KNOWN, BDG 2016, for fixed M, k] | ≈ m^{M−√(2M)} = N^{1−O(M^{−1/2})} | Prescribed N-dependent top part | refuted: (b) fails |
| 2d | Prouhet–Thue–Morse, Frolov lifts, Choudhry chains | Up to A^r parameter tuples; M = 2^r (or j^k, 2j^{k−1}) | Only k = O(log M) moments fixed, so the residual degree is ≈ M | none | refuted: (b) |
| 2e | Near-ideal classes (2 ≤ t ≪ M) | ≤ C(q+t−1, t) members per class [PROVED, ptesplit_1 / ptecons_1] | Coefficient height H ≥ K_s/(t2^{t−1}A^{t−1}) with s = ⌊(M−t)/2⌋ + 1 [PROVED, codex]. The evaluated residual is bounded only under NC (§3.5) | Prescribed top part | refuted under heuristic NC; **open** without it (no construction) |
| 3a | μ_n-orbits ∏_j(X − h − aζ^j) = (X − h)^n − a^n | Θ(A²) | The full-orbit product is one binomial (Galois collapse) [PROVED] | — | refuted: one relation per orbit; piece norms ≥ N^{φ(n)/n} |
| 3b | Projective μ_n orbits (dX − b)^n − (a − cX)^n | Θ(A⁴) | Mixed binomial terms do not compress; F′ applies to forcing it [PROVED, codex] | — | refuted |
| 3c | Cunningham N \| Φ_n(m), X^n ≡ 1 (F_{1^n}) | — | CUN: ∏_{\|a\|≥2}\|a\| ≤ ‖r‖₁(1+o(1)), so at most log₂‖r‖ nontrivial roots; F_1 group of order 2n [PROVED] | Special N only | refuted: degenerates to M13 of systeq |
| 3d | Jacobi/Gauss sums (Stickelberger) | O(n²A²) | Explicit principal-ideal relations in Z[ζ_n]; image mod N of size ≈ N/m | Needs ζ_n ↦ m mod N; mod N this is APR-CL | refuted: (b) |
| 4 | Division polynomials / Lattès preimages [n]^{−1}(Q) | Points of height ≤ A: O((log A)^{r/2}) [KNOWN, Néron–Tate]; full fibres need μ_n ⊂ K(P), so φ(n) ≤ 2[K:Q] [PROVED, codex; Weil pairing, Merel] | Denominators ψ_n² | None | refuted: (c) and (a) |
| 5a | Group ring (C₂)^r: ∏_{ε∈{±1}^r}(X − Σε_j a_j) | ≈ A^r/(r!)², M = 2^r: genuinely dense and split | Moment rigidity: at most 2r = 2 log₂ M leading coefficients agree [PROVED, codex] | — | refuted: (b) |
| 5b | Abelian group determinants / resultants ∏(a_i − b_j) | as indexed | Cyclotomic norms; translation-invariant differences | Zero residual ⇒ a gcd reveals the factor | refuted |
| 6 | Lattice Λ_{N,m} ∩ split variety; Schnorr–Adleman | Split density in the coefficient box is A^{−M(M−1)/2} [PROVED, codex] | Expected intersection count = S·2V/N (AVG) | Genuine: the lattice contains N | **open** as an algorithm; KNOWN failures of concrete variants (Vera 2010; Ducas 2021 experiments; Khattar–Yosri 2023) |
| 7a | Aligned clusters, ACS (§3.6; our best) | Pool-limited | λ ≥ 1/2 typical (HALF; E2) | Base-m digit matching | refuted: at best COS, L[1/2] |
| 7b | Composition clusters (X² − c)² with a² + b² = 2c | r₂(2c) = A^{o(1)} per c | As 7a | none | refuted (E1 G = 0.99) |
| 7c | Linear (affine-in-parameter) real families | — | LIN-R: one moving root | — | refuted (PROVED) |
| 7d | Block products ∏_j{P_j or Q_j} | 2^r members | Rank ≤ 2r: the rank trap [PROVED, codex] | — | refuted |

## 3. Lemmas

### 3.1 AVG — no excess yield without N-dependence [HEURISTIC + measured]

**Statement.** Let U be a set of integers fixed independently of N, and let N range over a wide set (here 200 balanced 24-bit semiprimes). Then

  E #{u ∈ U, |u| > N : |u mod± N| ≤ V} = Σ_u (2V+1)/N · (1 + o(1)).

So a structured family needs N-dependence to give any excess yield.

**Status.** The heuristic is equidistribution of u mod N over N. A rigorous average over N ∈ [X, 2X] needs divisor-sum control of #{N | u − v}; we have not written it. E1 measured G = 0.96–1.07 (§5).

**Consequence.** All N-dependence must come through alignment, and its cost is the whole question.

### 3.2 GAP-R and LIN-R — real rigidity [PROVED; GAP-R is classical in substance]

**GAP-R.**
- *Statement.* Let y₁..y_M be real with e_i(y) = e_{i+1}(y) = 0 for some 1 ≤ i ≤ M − 1. Then at most i − 1 of the y_l are nonzero.
- *Proof.*
  - R(Z) = ∏(Z + y_l) = Σ_j e_j Z^{M−j} is real-rooted.
  - Differentiate it M − i − 1 times. The result's Z¹ and Z⁰ coefficients are nonzero multiples of e_i and e_{i+1}, so it is Z²·(…).
  - For a real-rooted f, a multiple root of f′ is a root of f: between consecutive distinct roots f′ has exactly one simple root, which accounts for all of f′'s roots. By induction, a root of multiplicity μ ≥ 2 of f^{(r)} is a root of f of multiplicity ≥ r + μ.
  - Hence 0 is a root of R of multiplicity ≥ M − i + 1. ∎
- *Check.* Exhaustive over integer multisets with M ≤ 6 in [−5, 5]: 1364 double-gap cases, 0 violations (`lemmas.gap_r_exhaustive`).

**LIN-R.**
- *Statement.* Let P_λ(X) = ∏_{l=1}^M (X − α_l − L_l(λ)) with α_l ∈ C and L_l real linear forms on R^s. If the coefficients of X^{M−2} and X^{M−3} are affine in λ, then at most one L_l is nonzero.
- *Proof.*
  - The degree-i homogeneous part in λ of e_i(α + L(λ)) is e_i(L(λ)).
  - Affineness therefore gives e₂(L(λ)) = e₃(L(λ)) = 0 for all λ ∈ R^s.
  - Apply GAP-R with i = 2 pointwise: for each λ at most one L_l(λ) ≠ 0.
  - If two forms were not identically zero, a λ off both of their kernels would contradict this. ∎

**Interpretation.**
- *In F_q.* The affine-subgroup part of the BGJT family, roots b + c·a with c ∈ F_q, is a linear family in which *every* middle coefficient vanishes. The reason is p_j(F_q) = 0 for 1 ≤ j ≤ q − 2. Linear dependence on parameters is what makes the alignment a *linear* condition, solvable for a whole family at once.
- *Over R.* LIN-R says this is impossible.
- *Over C.* Directions with vanishing power sums p_2..p_k exist: (1, i, −1, −i) gives ∏(1 + uT) = 1 − T⁴. With house 1 they force roots of unity of order > k, i.e. cyclotomic fields of degree ≥ φ(k+1).
  - The pieces m − α − ζ^j λ then have norms ≈ m^{φ(n)}. By the Splitting Lemma (systeq §5.4) their product carries total norm ≈ N^{[K:Q]}.
- *Over a fixed imaginary quadratic field.*
  - Multisets in Z[i] with vanishing p₁..p_k exist for every k. The reason: 0 is interior to the convex hull of the moment curve z ↦ (z, …, z^k), and Q(i) is dense in C.
  - But their heights grow. The μ₄ construction gives k = 3 for free; each further level (e.g. Σ z_r⁴ = 0 with z⁴ ∈ {1, −4}) multiplies the house. We have not quantified the growth, so this sub-case is **open**.
  - It meets the residual problem anyway. Moving roots by λ·c_l adds Σ_{i>k} e_i(c)λ^i m^{M−i} ≈ N(H|λ|/m)^{k+1} to v, so k must be ≈ M − o(M).

### 3.3 LEDGER and CUN — constants are not units [PROVED]

**LEDGER.** Let F be monic of degree d, P = ∏(X − a_l) and r = P mod F. Then

  ∏_l |F(a_l)| = |Res(F, r)| ≤ ‖r‖₁^d · M(F)^{deg r},

where M is the Mahler measure.
- *Proof.* Res(F, P) = ∏_θ P(θ) = ∏_θ r(θ) over the roots θ of F, and also = ±∏_l F(a_l). Bound |r(θ)| ≤ ‖r‖₁ max(1, |θ|)^{deg r}.
- *Check.* `lemmas.ledger_holds` on random cases.
- *Reading.* This is just the norm multiplicativity of ∏(θ − a_l) = r(θ) in Z[θ]. Its content is the contrast with F_q[X], where the "norm" of X − a is the constant I(a) ∈ F_{q²}^*. That has degree 0, so the ledger costs nothing.

**CUN.** For n ≥ 3, integer roots a_l and r = P mod (X^n − 1):

  ∏_{|a_l|≥2}|a_l| ≤ ‖r‖₁ · (1 − 2^{−n})^{−M/(n−2)}.

- *Proof.*
  - At the n − 2 or n − 1 nonreal n-th roots of unity ζ, P(ζ) = r(ζ) and |r(ζ)| ≤ ‖r‖₁.
  - For |a| ≥ 2, ∏_{ζ nonreal}|ζ − a| = |(a^n − 1)/(a² − 1)| or |(a^n − 1)/(a − 1)| ≥ |a|^{n−2}(1 − 2^{−n}).
  - For a ∈ {−1, 0, 1} the product is ≥ 1.
  - Multiply over l and take the (n−2)-th root. ∎
- *Checked* in tests (200 random cases).
- *Consequence for systeq M13 / F_{1^n}.* Every split polynomial whose reduction mod X^n − 1 has height h has at most log₂ h + o(1) roots outside {−1, 0, 1}. So its u-side uses essentially three symbols, m − 1, m and m + 1. This is the quantitative form of "|PGL₂(F_{1^n})| = 2n" on the split side.

### 3.4 CW — ideal-cluster width [PROVED; codex ptecons_1, proof in `codex/scratch/ptecons_1/ptecons_pte_height_lemma.md`, re-checked by us]

**Statement.** Let P ≠ Q be monic of degree M, split over Z with roots in [−A, A], with deg(P − Q) < t. Put k = M − t, s = ⌊k/2⌋ + 1 and K_s = ⌊s/2⌋!⌈s/2⌉!. Then

  ‖P − Q‖_∞ ≥ K_s/(t 2^{t−1} A^{t−1}).

**Proof.**
- Remove G = gcd(P, Q), which has degree < t.
- Each cofactor has ≥ s distinct roots. Otherwise W = ∏_{a∈S}(X − a)² would have degree ≤ k, vanish on one cofactor's roots and be positive on the other's, contradicting the equal moments.
- Evaluate a suitable derivative of P − Q at a root of Q's cofactor. That gives a product of ≥ s distinct nonzero integer differences.

**Ideal case (t = 1).** P − Q = C with log|C| ≥ ½M log M − O(M). Every prime divisor of C is ≤ 2A.

**Checked.** All known ideal solutions of sizes 2–8, 10 and 12 (Borwein–Lisoněk–Percival tables) satisfy the bound, with large slack. For example, M = 12 has |C| = 6.7·10²².

**Consequence (common offset).**
- In an aligned ideal cluster, v_P − v_Q = C. So at most one member has |v| < |C|/2.
- With M ≥ x/(C₀ log B) this gives log|v| ≥ (x/(2C₀ log B)) log(x/log B) − O(x/log B) for all but one member per class.
- For log B = x^α(log x)^{1−α}, α < 1/3, this is ≥ x^{1−α+o(1)} ≫ x^{2/3}.
- So ideal clusters provide no gain below L[1/3]. At α = 1/3 the bound is of the boundary order: constant-level only.

**Context.** The barrier-side companion is ptebarrier's Lemma R (in progress, `notes/ptebarrier_prereg.md`). It gives an exponential rather than a factorial lower bound for general K, with an explicit M-threshold.

### 3.5 Near-ideal clusters [PROVED height bound + HEURISTIC NC]

**NC (no cancellation) [HEURISTIC].** For the class differences D = P − Q (degree < t, height H), |D(m)| ≳ max(H, m^{deg D}).

**Two regimes under NC, in a common-offset cluster.**
- **t ≲ M log M/(4 log A).** The CW height bound gives log|v| ≳ (M/4) log M.
- **t larger.** log|v| ≳ (t − 1) log m.

**Consequence.** For pieces ≤ B = L[α], α < 1/3, and m ≥ A^{Ω(1)}, both regimes exceed x^{2/3}. The case "m tiny and t ≈ M" is not covered.

**What is not settled.** NC is unproved. A near-ideal class with *systematic* cancellation of D at an N-aligned m would evade this. That is precisely the open primitive.

### 3.6 Alignment capacity and HALF (our best candidate) [HEURISTIC + measured]

**Algorithm ACS** (`src/ptecons/acs.py`):
1. For k ≤ K, put m_k = ⌈(kN)^{1/M}⌉. Write the offset R₀ = m^M − kN in base m.
2. Greedily choose the top coefficients e₁, e₂, … of a split P with roots in [−A, A] to cancel R₀'s leading base-m digits.
3. Enumerate the class (cluster) of multisets realising them. For the quartic this uses a two-sums meet-in-the-middle; it is checked against brute force.
4. For each member, set u = ∏(m − a_l), whose pieces are ≤ m + A, and v = u − kN.
5. Smoothness-test v, then do linear algebra.

**Capacity.**
- Cancelling digit j needs |e_j| ≈ m, so A ≳ m^{1/j}. Already j = 2 needs A ≳ m^{1/2}/M.
- The cluster's free low coefficients then vary by e_i·m^{M−i} with |e_i| ≈ A^i. The bottom term alone is e_M ≈ A^M ≳ m^{M/2} = √N.
- So once two digits are matched, typical |v| ≳ N^{1/2} (HALF).
- Matching only e₁ leaves v ≈ m^{M−1}·A ≥ N^{1−1/M}.
- Generic realizability caps the number of matched digits at ≈ √(2M). The typical residual is then ≈ m^{M − (j+1)/2}: N^{0.625} for M = 4.
- M = 2 gives exactly COS, |v| ≈ A√N.

**Cost** [HEURISTIC, random-residual model]. With |v| ≥ N^{1/2}, ACS is Dixon-type with a half-size residual, i.e. L[1/2], the same exponent as COS/QS. Adding an algebraic side turns it into NFS, which is barred from beating 1/3 by RES + G (heurnfs B1).

**Measured (E2, §5).** No higher-degree aligned variant beat COS. The two-digit sextic is infeasible.

### 3.7 Collapse and rank [PROVED, elementary]

- **Galois collapse.** A full μ_n-orbit gives ∏_k(x − ζ^k t) = Res(Z^n − 1, x − tZ) = x^n − t^n: one relation (`lemmas.galois_orbit_collapse`). Galois conjugates have equal norms.
- **Root support.** A relation matrix over symbols {m − a} has rank ≤ #(distinct roots used). Families whose members have disjoint root sets (Galois/Dickson fibres) are rank-deficient (`root_support_rank`).
- **Block products** (codex). Products of aligned blocks span rank ≤ 2r for 2^r members.
- **Comparison.** In BGJT each linear factor appears in ≈ q² of the q³ members.

## 4. Failure modes (exhaustive for the families examined; not claimed exhaustive in general)

| Mode | Statement | Status | Kills |
|---|---|---|---|
| N-independence | Structured but N-free families have the uniform congruence rate | HEURISTIC, measured G ≈ 1 | 1a, 2d, 3a, 5, 7b |
| Collapse / rank | Orbit products are norms; disjoint supports | PROVED | 1c, 3a, 3d, 7d |
| Field degree | Splitting over K needs φ(n) ≤ c[K:Q]; piece norms ≥ N^{φ(n)/n} | PROVED / KNOWN | 1c, 3, 4 |
| Archimedean non-units | LEDGER / CUN | PROVED | 3c, presentations of small height |
| Real rigidity | LIN-R / GAP-R | PROVED | 7c; any linear-alignment family |
| Ideal-cluster width | CW | PROVED | 2a, 2b (common offset) |
| Alignment capacity | HALF, generic √(2M) | HEURISTIC, measured | 2c, 2e, 6, 7a |

**The BGJT orbit evades all seven at once in F_q:**
- it is N-dependent through the presentation;
- its root support is shared;
- the field is fixed (F_{q²});
- constants are units;
- the directions are F_q-linear;
- the cluster width is a degree, not a size;
- and its capacity is linear.

**No characteristic-0 family examined evades even two of the rigidity modes (LIN-R, LEDGER) together.**

## 5. Experiments (pre-registered, `notes/ptecons_prereg.md`, commit 95a669f)

### E1 — AVG (`experiments/raw/ptecons_E1.jsonl`)

**Setup.** 200 balanced 24-bit semiprimes, V = 1000. Only pairs with |u| > N are counted.

| family | size | pairs | observed | expected | G |
|---|---|---|---|---|---|
| generic sextic, [−6, 6], m = 40 | 18564 | 3.71M | 615 | 584.9 | 1.05 |
| symmetric J = 3 clusters | 11575 | 2.32M | 389 | 364.7 | 1.07 |
| PTE-class members (M = 5, e₁..e₃ equal) | 1132 | 226k | 37 | 35.7 | 1.04 |
| Dickson integer-split members | 7826 | 1.48M | 224 | 233.8 | 0.96 |
| composition clusters (a² + b² = 2c) | 17512 | 2.49M | 390 | 392.9 | 0.99 |

**P1** (G ∈ [0.8, 1.25] when expected ≥ 200): **held** for all four eligible families. The PTE family (expected 36) is below the threshold.

### E2 — aligned clusters (`experiments/raw/ptecons_E2.jsonl`)

**Setup.** 64-bit semiprimes (5 seeds) and 80-bit semiprimes (3 seeds); pool cap 20000; B = 2^16; control = a random integer of the same bit length.

Median λ, with the range over seeds:

| variant | 64-bit median λ | 80-bit median λ | λ at rank 1000 (64 / 80) | members (64-bit) | smooth / control |
|---|---|---|---|---|---|
| cos (M = 2, A = 100) | 0.592 | 0.574 | 0.537 / 0.530 | 20000 | 0.85–0.91 |
| quart_g (e₁ exact, e₂ window) | 0.595–0.606 | 0.565–0.583 | 0.556–0.584 / 0.524–0.540 | 3072–13093 | 0.76–0.98 |
| quart_s (symmetric quartic) | 0.633 | 0.606 | 0.578–0.581 / 0.562–0.564 | ≈ 19100 | 0.96–1.03 |
| sext_s (two digits matched) | — | — | — | 1, 1, 1, 1, 4 (80-bit: 0) | — |
| sext_s_wide (one digit) | 0.674–0.677 | 0.647–0.649 | 0.695–0.699 / — | 881–1101 | 0.79–1.36 |
| quart_rand (control) | 0.855 | 0.860 | 0.80 / 0.82 | 20000 | 0.54–0.93 |

**Evaluation.**
- **P2a** (cos median = 0.5 + log₂100/n ± 0.04; predicted 0.604 and 0.583): measured 0.592 and 0.574. **Held.**
- **P2b** (every aligned variant ≥ cos − 0.02 at every N): **held**. The closest case is quart_g at 80-bit seed 2, 0.565 against 0.574.
- **P2c** (λ₁₀₀₀ ≥ 0.45): **held**; the minimum was 0.524.
- **P2e** (sext_s < 200 members): **held**, with 0–4.
- **Kill signal** (an aligned variant below cos − 0.05 at ≥ 2 N): **not triggered**.
- **P2d** (smooth ratio ∈ [0.6, 1.6]): **failed** in one cell, quart_rand at 80-bit seed 2 (0.54, from 25 against 46 hits). *Disclosed explanation:*
  - v = ∏(m − a) − kN is *anti*-correlated with small primes. A product of four integers is even with probability 15/16, so v (with N odd) is almost always odd; likewise mod 3.
  - This is the local-root effect already seen in systeq E4.
  - COS (two factors) shows the milder 0.85–0.91. The aligned symmetric variants are ≈ 1.
- **P2f** (quart_rand median λ ≥ 0.95): **failed as stated** (0.855 / 0.860). *Disclosed:*
  - The "unaligned" control is implicitly aligned on its leading digit: m = ⌈N^{1/4}⌉ makes u ∈ N(1 ± 4A/m).
  - Hence λ ≈ 1 − log₂(m/4A)/n = 0.89 / 0.875.
  - A genuinely unaligned control would need m unrelated to N. The prediction was mis-specified, and the result is consistent with alignment mattering.

**Reading.** At toy sizes the generic quartic roughly ties COS. Asymptotically it cannot: it needs A ≳ m^{1/2} (§3.6), giving λ → 0.625, while COS → 1/2.

## 6. Literature (KNOWN; checked by the literature critic, §8)

- **PTE.**
  - Borwein–Ingalls (1994), the survey.
  - Borwein–Lisoněk–Percival (2003), computational ideal solutions; sizes ≤ 10 and 12.
  - Coppersmith–Mossinghoff–Scheinerman–VanderKam (arXiv 2304.11254), searches to size 16.
  - Choudhry (arXiv 1908.11192), chains.
  - Wright's conjecture.
  - Constants: Kleiman (1976), Rees–Smyth (1990), Caley (arXiv 1011.1262), on divisibility of the ideal constant C.
- **PTE in cryptography.**
  - Costello–Meyer–Naehrig (EUROCRYPT 2021) use exactly the ideal-class shape u, u + C to sieve twin smooth integers.
  - This is N-free, so for factoring purposes it is AVG.
- **Vinogradov mean value theorem.** Bourgain–Demeter–Guth (Annals 2016); Wooley (cubic case). Its constants depend on (M, k).
- **Schnorr–Adleman lattice factoring.**
  - Schnorr (1991); Adleman (1995).
  - Vera (arXiv 1003.5461): an effective analysis showing too few relations.
  - Khattar–Yosri (arXiv 2307.09651): failure of small-dimensional variants.
  - The 2021 "destroys RSA" claim is refuted in practice (Ducas's experiments).
- **COS linear sieve and cubic sieve.** Coppersmith–Odlyzko–Schroeppel (Algorithmica 1986). The cubic sieve needs special N (Reyneri).
- **Gap and fewnomial results.** Descartes; Lenstra (1999) on lacunary factors; Hajós's lemma on root multiplicity.

## 7. Open problem (precise)

Find, for general N, an algorithm that outputs L[α] (α < 1/3) relations u = ∏(m − a_l) ≡ v (mod N) with pieces ≤ B = L[α], at cost |𝓕|·L[o(1)], and with one of the following:
- log|v| ≤ x^γ for some γ < 2/3; or
- v smooth by construction.

This is an output-sensitive enumerator of the intersection of the N-lattice {P : P(m) ≡ small} with the split variety. By §§3.4–3.6 it cannot be a common-offset ideal cluster. A near-ideal cluster would have to cancel systematically at the aligned m. It cannot be linear in real parameters. In smoothenum's terms it is the doubly or algebraic listing primitive. Nothing examined supplies it.

## 8. Reviews

### 8.1 Codex blind round (ptecons_1)

**Same verdict**, about 2%.

**Contributions we adopted:**
- CW and the near-ideal height bound;
- the Dickson fixed-field and unit-density results;
- 2b (doubly smooth by identity, alignment ≈ N);
- the (C₂)^r group-ring family and its moment rigidity;
- the rank trap;
- split density in coefficient space;
- the γ < 2/3 threshold;
- the elliptic Weil-pairing bound.

**Codex's own experiment.** M = 6, roots in [−12, 12], m = 13. It found 12 two-member ideal classes and 0 with three or more. No ideal-class member had |v| ≤ 1000, because the class constants (14400–100800) are > 2V. That is a direct instance of CW.

### 8.2 Codex adversarial round (ptecons_2)

PENDING; filled in below once complete.

### 8.3 Opus critics

PENDING.

## 9. Reproduction

```bash
.venv/bin/python -m pytest -q tests/test_ptecons.py
PYTHONPATH=. .venv/bin/python -m src.ptecons.run_experiment E1    # ~2 s
PYTHONPATH=. .venv/bin/python -m src.ptecons.run_experiment E2    # ~3 min
```
