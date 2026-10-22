# sinonobliv: the non-oblivious routes to Prop SI that sihard left open

Agent `sinonobliv`, phase 22, 2026-10-22. Version 3 (after the codex adversarial round and both Opus critics).

**Files.**
- Code: `src/sinonobliv/` — `core.py` (exact smoothness masks, density curves, identity families, divisors of
  lcm(1..B), certificate thresholds), `msample.py` (multi-sample certificates, Lemma IH scan), `uncert.py`
  (Howgrave-Graham lattice, uncertified roots), `theory.py` (cost formulas), `run_exp.py` (E1–E3), `posthoc.py`
  (disclosed post-hoc analyses).
- Tests: `tests/test_sinonobliv.py` (12 pass).
- GPT-6.1-sol adversarial round: `codex/sinonobliv_2.{md,out}` (§9.2).
- Pre-registration: `notes/sinonobliv_prereg.md`, commit 424ee86 (E1, E2), addendum A1 commit 01eb40d (E3), each made
  before the runs it registers.
- Raw data: `experiments/raw/sinonobliv_E{1,2A,2B,3}.jsonl`; post-hoc `sinonobliv_E1_posthoc.jsonl`,
  `sinonobliv_E2B_posthoc.jsonl`.
- GPT-6.1-sol: blind round `codex/sinonobliv_1.{md,out}` (§9.1).

**Labels.** KNOWN (cited) · PROVED (proof here) · PROVED in MODEL (exact consequence of a stated model) ·
HEURISTIC (named assumption) · MEASURED · POST-HOC · OPEN · UNVERIFIED · assessment. "No prior statement found by
us" is the strongest novelty claim made.

**Notation** (as in sihard). X ≈ kN, x = ln X, F_α = x^α(ln x)^{1−α}, L[α, c] = exp((c + o(1))F_α). Factor base
y = L[1/3, κ], window h = L[2/3, ν], t := ln(h + 1), u₀ = x/ln y, ρ the Dickman function, so
ρ(u₀) = L[2/3, −2/(3κ)]. Prop SI's optimum κ* = (2/9)^{1/3}, ν* = 2/(3κ*). P denotes a y-smooth modulus, π := ln P.
A constructor costing L[1/3, τ] per output beats MNFS (1.9019) iff τ < 0.75717.

---

## 0. Verdict

**Plain language.**
- **We found no constructor below L[2/3] per output on any of the five routes.** So no finite τ was obtained,
  let alone τ < 0.757. The best constructor in every route is still the known divisor trick z = d·⌈x/d⌉ (or the
  interval sieve), at L[2/3, 2/(3κ)] per output up to an L[1/3] factor.
- **Two routes now have barriers in stated models (narrow ones), one reduces to Prop SI itself, and one has no
  barrier.** None of this closes a route for arbitrary algorithms.
  - *Multi-sample (Cohn–Heninger) decoders* — barrier. With independently placed windows, a joint certificate
    yields at most ρ(u₀)^{1−o(1)} smooth outputs per call, for every modulus and every 2 ≤ m ≤ ln h, even when Cohn–Heninger's
    enlarged decoding region is granted for free (Prop MS, PROVED, averaged). With the real windows near k_iN of one
    modulus, every joint certificate with a divisor d > Kh is a disguised single window: the offsets are forced to be
    proportional and all m outputs give the same relation (Lemma IH, PROVED; codex found it too); its yield is
    Theorem J's once one grants that dependent samples carry one sample's information (Prop IH-J, MODEL). Below d ≤ Kh the
    untested cofactor is ≥ N/(Kh) and the yield per certified tuple is again ρ(u₀)^{1−o(1)}, for decoders that
    examine all certified tuples (Prop MS′).
  - *Correlated polynomial identities* (t(t+1), products of consecutive integers, Aurifeuillian and Friedlander–
    Lagarias products, powers, PTE clusters) — barrier for families fixed in advance. Their smoothness gain is real
    (measured 8× to 1100× at toy sizes; r^{u₀} asymptotically), but a usable fixed family must contain at least a
    p/((h+1)ρ) fraction of *all* y-smooth integers near X, p being its per-window success rate (Prop FF, PROVED,
    averaged). Every correlated family we know (one parameter, or balanced-range forms of degree above their
    dimension) has O(X^{1−ε}) distinct values near X and lands in a window with probability ≤ (h+1)X^{−ε}:
    exponentially small. A dense family is a large share of the smooth integers of some shape, and listing it in a
    window is Prop SI restricted to that shape (assessment).
  - *Smooth norm and polynomial values near (kN)^{1/d}* — no gain; the route reduces to Prop SI (a reduction, not a barrier). Norms of constructed elements form a sub-knapsack
    of the integer one (Prop NV, PROVED); one-parameter values t^d − a are exponentially sparse (Lemma PL); a
    degree-d analogue of the Prop SI relation must use an element of height ≥ N^{1/d} (Prop LC; in Q(i) its norm is
    ≥ N^{1−o(1)}; the general norm scale is open). The NFS
    sampling model reproduces SNFS 1.526, GNFS 1.923 and MNFS 1.902 (codex), but those are relation generators, not
    Prop SI constructors.
  - *Adaptive choice of moduli/lattices* — no barrier, and we explain why none of the information-theoretic kind can
    exist. Fixed chains and fixed adaptive trees obey Theorem J (Prop JC, Prop AT). But every bit of information
    needed is cheap: the residues of x modulo all prime powers ≤ y cost L[1/3, κ] and determine every smooth position
    in the window (Lemma IC). What remains is a computational problem: list-decoding a Chinese-remainder code far
    below its Johnson radius (Prop RD; this is Boneh's generalized CRT decoding problem at full smoothness). A decoder of cost π(y)^a per window would beat MNFS iff a < 1.1416. This is
    the one place a constructor could hide, and it is the original problem in coding language.
- **Uncertified lattice outputs** (sihard's fifth item; E3, pre-registered): the extra integer roots of a
  Howgrave-Graham divisor lattice are "near-certified" points whose gcd with P sits just below the provable threshold
  (median e^{18.1} vs e^{18.9}); most are even inside Theorem J's ideal Johnson class. The 42 outside it (in 5000
  calls) are smooth at 0.60 of what their known divisor predicts. No evidence of a hidden bias, at toy size and with
  small counts.
- **Assessment.** Probability that any of these routes yields τ < 0.757: about 1% (ours); 5% (codex, blind),
  mostly on target-dependent adaptive constructions. We agree that route 1 is where any remaining probability sits.

**Per-route table.**

| Route | Best construction found | Per-output cost (τ) | Barrier / model | Status |
|---|---|---|---|---|
| 1. Adaptive moduli/lattices | divisor trick; fixed chains = Theorem J | L[2/3, 2/(3κ)]·L[1/3, −2ν/(3κ)]; **no finite τ** | Prop JC (fixed chains), Prop AT (fixed trees, M nodes: ≤ M·ρ^{1−o(1)}); no information barrier possible (Lemma IC); residual problem = CRT decoding below Johnson (Prop RD) | barrier PROVED in MODEL; unrestricted OPEN |
| 2. Correlated identities | t(t+1), t(t+1)(t+2), 4t⁴+1, FL products, t^d | sparse families: ≥ X^{ε}/h per output (exponential); dense families: = Prop SI | Lemma PL, Prop FF (fixed families, averaged) | PROVED (placement), HEURISTIC (gain); target-dependent growing families OPEN (→ route 1) |
| 3. Multi-sample decoders | Cohn–Heninger PACD/ACD on m windows | ≥ ρ^{−1+o(1)} calls per output (PACD region, any P, 2 ≤ m ≤ ln h); ≥ L[2/3 + 1/(3m)] (general ACD, codex) | Prop MS (independent windows), Lemma IH + Prop MS′ (same N) | PROVED (averaged, region granted); adaptive helper windows OPEN |
| 4. Norm / polynomial values | t^d − a; Norm(α); degree-d coset analogue | exponential (placement) or = integer knapsack; no finite τ obtained | Lemma PL (PROVED); Prop NV, Prop LC (reductions / height bound, not barriers) | no gain; target-dependent fields OPEN (→ route 1) |
| 5. Softer target | — | none of the constructions considered reaches any L[1/3, τ] (a finding, not a theorem) | threshold τ < 0.75717; decoder π(y)^a needs a < 1.1416 | threshold PROVED (algebra of the cost model) |
| (uncertified outputs) | Howgrave-Graham extra roots | as the divisor trick | E3: smoothness explained by the known divisor | MEASURED (toy) |

---

## 1. Setting

- **Prop SI** (suffmap §1.4; sicheck §1.2): if y-smooth u ∈ I_k = [kN, kN + h] can be constructed (with its
  factorisation) at L^{o(1)} per output, N factors in L[1/3, 1.2114] under H_S on the residue and R. With per-output
  cost L[1/3, τ] the constant is c(τ) = min_κ max(κ + τ + 2/(9κ²), 2κ).
- **What sihard proved** (sihard §5): every *oblivious* model collapses to L[2/3] per output. Theorem J bounds
  single-integer Johnson-class certificates with parameters fixed in advance by ρ(u₀)^{1−o(1)} certified smooth
  outputs per call, averaged over the window position; the divisor trick attains it. sihard's list of mechanisms
  *not* excluded is this note's task: adaptive moduli, correlated identities, multi-sample decoders, norm values,
  uncertified outputs.
- **Constraints inherited.** siconstruct's Lemma SR (every factor of u ∈ [X, X + h] faces relative precision h/X)
  and Lemma RC (an independently placed factor s ≲ Kh leaves a tested cofactor ≥ N/(Kh)); sicheck's C1–C8 (in
  particular log k = o(log N), the fresh-relation rate C4 and fruitfulness C5); ptebarrier's Lemma R / Theorem W
  (exact split compressions give ≤ 1 relation per presentation).

## 2. Route 1: adaptive choice of moduli and lattices

### 2.1 Fixed chains obey Theorem J

**Prop JC (PROVED; reduction to Theorem J).** A *chain call* fixes a y-smooth D ≥ 1 and a Johnson-class certificate
(P, w) for the cofactor: on the window [x, x + V] it certifies z with D | z and W(z/D) ≥ √(Q·ln(V/D)) (the
Coppersmith call on the progression D·Z has window V/D in the cofactor variable). For V ≤ L[2/3, ν_max], uniformly
over D and (P, w):
  E_{x∈[X,2X)} #{certified y-smooth z} ≤ ρ(u₀)^{1−o(1)}.
Depth-r chains (D = d₁⋯d_{r−1} fixed, the last call a certificate) are the case D = d₁⋯d_{r−1}.

*Proof.* z = Dm is y-smooth iff m is (D is y-smooth).
- *Case D ≤ V.* Each z lies in ≤ V + 1 windows, so the average is ≤ ((V + 1)/X)·#{m ≤ 3X/D : m certified at window
  V/D, m y-smooth}. With X′ = X/D and V′ = V/D ≥ 1, (V + 1)/X ≤ (V′ + 1)/X′, and the right side is exactly the
  quantity bounded in steps 2–8 of Theorem J's proof at scale (X′, V′). Since X′ ≥ X/V = X^{1−o(1)}, ln X′ ∼ x and
  those steps hold with the same y, giving ρ(ln X′/ln y)^{1−o(1)} = ρ(u₀ − ln D/ln y)^{1−o(1)}. Finally
  ρ(u₀ − δ)/ρ(u₀) ≤ exp(δ(ln(u₀ ln u₀) + O(1))) with δ = ln V/ln y = O(F_{2/3}/F_{1/3}), an exp(O(F_{1/3})) factor.
- *Case D > V.* The region has at most one point, D·⌈x/D⌉; no lattice is needed. The average is
  ≤ ((V + 1)/X)·Ψ(3X/D, y) ≤ 3C₀((V + 1)/D)·ρ(u₀ − ln(D/3)/ln y), and ρ(u₀ − ln D/ln y)/ρ(u₀) = D^{o(1)}
  (same ratio bound), so the average is ≤ ρ(u₀)·(V/D)^{1−o(1)}·e^{o(F_{2/3})}. ∎

Three caveats (proof critic). The cofactor score W(z/D) models a real Coppersmith call on f(m) = Dm + x only when
P is coprime to D (D must be invertible mod P). The weighted case inherits sihard's "PROVED (sketch)". The uniformity
needed is in n = ln X′ ∈ [x − νF_{2/3}, x], and Theorem J's o(1) depends only on n.

**Prop AT (adaptive trees with a fixed repertoire; PROVED; codex blind §1, = sihard J(c)).** Let the call at each
node of a decision tree be a chain call whose parameters are fixed in advance, with previous answers selecting the
child. If the tree has M distinct nodes, the number of smooth certified outputs on input x is
≤ Σ_{nodes θ} S_θ(x), so its average is ≤ M·ρ(u₀)^{1−o(1)}. A tree of depth ≤ B with ≤ J leaves has ≤ 1 + BJ nodes;
constant success with B = L[1/3, τ] calls needs ln J ≥ (2/(3κ) − o(1))F_{2/3} − τF_{1/3}. ∎

**Why this is not a running-time bound (codex blind §1).** A binary tree of polynomial depth has exponentially many
leaves, and averages of fixed-parameter calls cannot control adaptive ones: if A is uniform on {0, …, M − 1} and
certificate θ succeeds iff A = θ, every fixed certificate succeeds with probability 1/M, but ⌈log₂M⌉ binary
queries find A. Codex also gives a conditional barrier: if every reachable history H satisfies
Pr(A ∈ E | H) ≤ e^I Pr(A ∈ E) with I = o(F_{2/3}), then calls per output are ≥ e^{−I}/ρ(u₀)^{1−o(1)} (PROVED in that
MODEL). The next lemma shows that no real algorithm needs to respect that hypothesis.

### 2.2 Information is cheap

**Lemma IC (PROVED, elementary; this is the sieve — or Bernstein's batch smooth-parts computation — seen as
information).** For x ≤ 3X let e_q = ⌊log_q(3X + h)⌋ and R(x) := (x mod q^{e_q})_{q ≤ y}.
- (a) R(x) determines the set S(x) = {j ∈ [0, h] : x + j is y-smooth}: j ∈ S(x) iff
  Σ_{q ≤ y} ln q·#{1 ≤ e ≤ e_q : j ≡ −x (mod q^e)} = ln(x + j).
- (b) R(x) costs π(y) reductions of an x-bit integer, i.e. L[1/3, κ + o(1)] bit operations. For the windows I_k of
  one N, R(kN) = k·R(N) mod q^{e_q} costs O(π(y)·x) bit operations per k after one such precomputation (the moduli
  q^{e_q} have ≈ x bits).
- *Proof.* (a) v_q(x + j) = #{e ≥ 1 : q^e | x + j}, all such q^e ≤ 3X + h, and x + j is y-smooth iff
  Σ_q v_q(x + j) ln q = ln(x + j). (b) is direct. ∎
- *Consequence (a statement about query models, not an algorithm).* In a model that charges per residue query (or
  per bit learned about x) and allows unrestricted computation on the answers, the query cost of Prop SI is at most
  L[1/3, κ + o(1)] per window. This does not output anything efficiently, and it assumes a smooth integer exists in
  the window. Since x is public anyway, the point is negative: a query restriction cannot be what makes Prop SI
  hard. In cost-model terms an information cost τ = κ (a = 1 in Prop RD) gives
  c = 6^{1/3} = 1.8171 < 1.9019, after re-optimising κ. Codex's bounded-likelihood hypothesis
  fails for real algorithms: after R(x) the posterior on S(x) is a point mass. **So every barrier for route 1 must be
  computational, not informational.** Theorem J and Props JC/AT are barriers only because they restrict which
  computations are allowed (a fixed repertoire).

### 2.3 The residual problem is CRT list decoding below the Johnson radius

**Prop RD (PROVED, elementary; a specialisation of Boneh's generalized CRT decoding).** Let C_h = {c(j) = (j mod q^e)_{q ≤ y, e ≤ e_q} : 0 ≤ j ≤ h}, a
Chinese-remainder (CRT/Mandelbaum-type) code with coordinate weights ln q, and let the received word be
w(x) = (−x mod q^e). The weighted agreement of c(j) with w(x) is A(j) = ln(y-smooth part of x + j). Hence
**SMOOTH-INT(x, h, y) ⇔ find j ∈ [0, h] with A(j) ≥ ln(x + j)**: a weighted list-decoding problem.
- *Prior art (literature critic).* This is Boneh's "generalized CRT list decoding" problem (finding integers in an
  interval with a large smooth part; Boneh, STOC 2000, 265–272; JCSS 64 (2002) 768–784, §§3–4), with the threshold
  set to full smoothness; the far-below-Johnson reading is already in pilist §2.3. The coordinates q, q², … are not
  pairwise coprime: the code is the CRT code modulo ∏q^{e_q} with agreement ln gcd(x + j, P), Boneh's gcd form.
  In coding terms the regime is far below the Johnson radius but far above capacity (x ≫ ln h).
- Total weight W = Σ_q e_q ln q ≍ π(y)·x; message length ln(h + 1) = νF_{2/3}; required agreement ≈ x.
- The Johnson radius of known decoders (Goldreich–Ron–Sudan, IEEE Trans. IT 46 (2000), a √2 factor short; Boneh,
  who removes the √2; Guruswami–Sahai–Sudan, FOCS 2000, weighted) is
  agreement ≥ √(ln(h + 1)·ln P) for the modulus P used (Theorem J's threshold). With P = ∏ q^{e_q}, √(νF_{2/3}·W) ≫ x
  by an L[1/3] factor: far below Johnson.
- The list is short *on average*: averaged over x, its expected size is (h + 1)ρ(u₀)^{1+o(1)} = L[2/3, ν − 2/(3κ)],
  which is O(1) if h is chosen as c/ρ(u₀). This is an average, not a worst-case or promise bound; with leading
  constants cancelling, only exp(o(F_{2/3})) follows in general, which may exceed every power of π(y) (codex round 2).
  So on average the combinatorial obstruction to beyond-Johnson decoding (super-polynomial lists for full-length
  Reed–Solomon codes over F_{2^m}; Ben-Sasson–Kopparty–Radhakrishnan, IEEE Trans. IT 56 (2010) 113–120) does not apply; for individual received words it is not
  excluded.
- Capacity results. Folded/multiplicity codes (Parvaresh–Vardy; Guruswami–Rudra; the polynomial-ideal framework of
  Bhandari–Harsha–Kumar–Sudan, arXiv 2103.07930) change the code or work over polynomial rings. Plain prime-field
  Reed–Solomon codes now have capacity list decoders at constant rate (Brakensiek–Chen–Putterman–Zhang–Zheng, arXiv
  2609.08005, preprint; also Jeronimo, arXiv 2609.05870; cited in pilist §2.3), but their cost exponent blows up as
  the agreement fraction tends to 0, and here it is ≍ x/W ≈ 1/π(y). Guruswami's number-field codes (IEEE Trans. IT 49
  (2003)) are list-decoded only up to Johnson (Biasse–Quintin, arXiv 1107.2321). Interleaved-CRT/power decoding
  (Li–Sidorenko–Nielsen, ISIT 2013) goes beyond unique decoding under probabilistic assumptions, not to this regime.
  Cohn–Heninger call multi-input ACD the integer analogue of Parvaresh–Vardy/Guruswami–Rudra decoding; that "folding"
  is route 3, treated in §4. No CRT or number-field decoder for this regime is known to us; OPEN.
- *Softer target in decoder terms (PROVED, algebra).* A decoder costing L[1/3, aκ] per window (polynomial of degree
  a in π(y)) gives c = (3/2)(4/9)^{1/3}(1 + a)^{2/3} for a ≥ 1/3. It beats MNFS iff **a < 1.1416**; a = 1 (linear in
  the information size) gives 6^{1/3} = 1.8171 (siconstruct §1).
- The explicit-prime-list version of this decoding problem is NP-hard under randomized reductions (DIVISOR,
  credited to Alon–Kilian in Sudan 2010, proof by Shor, cstheory.SE Q. 4769, 2011; sihard Prop XS). With all primes ≤ y, nothing is known.

### 2.4 Natural adaptive schemes, costed
- *Prime-by-prime accumulation* (choose d₁ | z, recurse in the progression): while D ≤ h the cofactor is
  ≥ X/h and its smoothness costs L[2/3, 2/(3κ)] (Lemma RC; codex §1); once D > h there is ≤ 1 candidate, and choosing
  D so that it is smooth is the placement problem itself. No finite τ.
- *Residue-driven selection of P* (keep primes whose residue class hits a chosen sub-window): selecting j hit by
  many primes is sieving the window, cost ≈ h. No finite τ.
- *Adaptive helper windows for multi-sample decoding*: reduces to choosing shifts as functions of x; covered only by
  Prop AT. OPEN.

## 3. Route 2: correlated (Friedlander–Lagarias-type) identities

### 3.1 The gain is real

**Density gain (HEURISTIC H_piece; KNOWN in substance).** If a member z ≈ X factors by algebra into r pieces of
sizes X^{c_i} (Σc_i = 1) that are smooth independently at their sizes, then
Pr[z smooth] ≈ ∏ρ(c_iu₀) = ρ(u₀)·exp(u₀H(c) + o(u₀)), H(c) = −Σc_i ln c_i ≤ ln r (codex §2). Related literature: smooth polynomial values (Martin, JNT 93 (2002)); consecutive smooth
numbers (Teräväinen, Forum Math. Sigma 2018; Conrey–Holmstrom–McLaughlin, Exp. Math. 22 (2013)); PTE
parametrisations raising smoothness (Costello–Meyer–Naehrig, EUROCRYPT 2021); smooth twins from lattice vectors
(Mulder–Sterner–van Woerden, arXiv 2509.17699). For fixed r this gain
r^{u₀} = exp(Θ(x^{2/3}(ln x)^{−2/3})) is larger than every L[1/3] but is L[2/3, o(1)]: the cost stays
L[2/3, 2/(3κ)]. With r = x^δ pieces the cost is L[2/3, (2/3 − δ)/κ]; the gain only reaches L[1/3] scale when the
pieces have size y^{O(1)}, i.e. when the member is "smooth by construction" (codex §2).

*Measured (E1; §8.1):* forced pieces are smooth independently up to a local-residue constant (R = 0.46–1.02 raw;
0.84–1.00 after the post-hoc local correction), and the gains over random integers of the same size are 8–22×
(4t⁴ + 1), 16–170× (t(t + 1)) and 220–1100× (t(t + 1)(t + 2); two cells with unreliable controls omitted) at 40–60 bits. The irreducible control t² + 1 has
gain 0.1–0.5 (its values avoid primes ≡ 3 mod 4).

### 3.2 But correlated families do not land in windows

**Lemma PL (placement; PROVED, averaged; codex §2, §4).**
- (a) For any multiset 𝓕 of positive integers fixed independently of x:
  E_{x∈[X,2X)} #(𝓕 ∩ [x, x + h]) ≤ (h + 1)·#(𝓕 ∩ [X, 2X + h])/X.
- (b) If F ∈ Z[t] has degree d ≥ 2, positive leading coefficient and F′ ≥ c·X^{1−1/d} on the branch where
  F ≍ X, then each window contains ≤ 1 + h/(cX^{1−1/d}) members (≤ 1 eventually), and the probability that a random
  window contains one is ≤ C(h + 1)X^{1/d − 1}. For M presentations fixed in advance: ≤ CM(h + 1)X^{1/d−1}, so a
  constant hit rate needs ln M ≥ (1 − 1/d)x − νF_{2/3}: exponential.
- *Proof.* (a) Each member lies in ≤ h + 1 windows. (b) The mean value theorem bounds the count; (a) gives the
  probability, with #(F(Z) ∩ [X, 3X]) = O(X^{1/d}). ∎
- *Instances:* t(t + 1)⋯(t + r − 1) (hit rate ≈ (h + 1)/(rX^{1−1/r})), 4t⁴ + 1 and other Aurifeuillian
  presentations, t^d − a with |a| ≤ A ((A + h)X^{1/d−1}, codex §4), t^d (smoothness ρ(u₀/d) but hit rate hX^{1/d−1}).
  E1's P2 confirms the hit law exactly at toy size (ratios 1.000, 1.001).
- **Friedlander–Lagarias** (JNT 25 (1987), Thm 4; checked by codex in the primary text): the products
  (m − a₁)⋯(m − a_{r−1})(m + Σa_j) with |a_j| ≪ m^{1/2^j} place X^{1/r+ε}-smooth integers in every interval
  [X, X + X^{β_r}], β_r = 1 − 2/r + 1/(r2^{r−1}) = 1/4, 5/12, 17/32 for r = 2, 3, 4. The exponent β_r *increases*
  with r: the method trades smoothness for interval length in the wrong direction for Prop SI, and no fixed r reaches
  y = X^{o(1)}. FL credit the case r = 2 independently to Balog and Sárközy, and note that their set S_r has
  O(X^{2/r − 1/(r2^{r−1})}) elements below X ("well-spaced"): (#members)·(window) ≈ X, i.e. FL's family is an
  equality case of Lemma PL(a), the counting precedent for Lemma PL / Prop FF. The identity 4t⁴ + 1 is Sophie
  Germain's (a⁴ + 4b⁴ with a = 1).

**Prop FF (fixed families must be dense in the smooth numbers; PROVED, averaged).** Let a constructor output, for the
window [x, x + h], only y-smooth members of a multiset 𝓕 fixed independently of x (𝓕 may depend on the size X;
members carry their factorisations). If its success probability per window is ≥ p, then
  #{distinct y-smooth members of 𝓕 in [X, 2X + h]} ≥ p·X/(h + 1) = p·Ψ(2X, y)/((h + 1)ρ(u₀)^{1+o(1)}).
- *Proof.* The probability that [x, x + h] meets a set S of integers is ≤ (h + 1)·#(S ∩ [X, 2X + h])/X (each element
  lies in ≤ h + 1 windows); repeated presentations do not raise it (codex round 2). Ψ(2X, y) = Xρ(u₀)^{1+o(1)}. ∎
- *Reading.* A constructor that spends ≥ 1 unit per window it queries and costs L[1/3, τ] per output must have
  p ≥ L[1/3, −τ]. So a fixed family must contain at least a p/((h + 1)ρ(u₀)^{1+o(1)}) fraction of *all* y-smooth
  integers near X. Only this exact form is claimed: (h + 1)ρ(u₀) can be exp(√(F_{1/3}F_{2/3})) or larger while
  h = L[2/3, 2/(3κ)] still holds, so the fraction is not an L[1/3]^{−1} statement (codex round 2).
  Every family with O(X^{1−ε}) distinct values near X fails by a factor X^{−ε}: exponential. This covers
  one-parameter polynomials of degree ≥ 2, Aurifeuillian and cyclotomic values, PTE fibres of bounded size, and
  coercive (balanced-range) m-parameter forms of degree d > m such as st(s + t) on a box. Degree alone does not
  suffice: st(t + 1) with t = 1 contains every even integer (codex), i.e. unbalanced ranges are the divisor trick. A family that passes contains, up to multiplicity, a large share of all
  smooth integers of some factorisation shape (e.g. all products of r balanced factors, which covers almost all
  smooth integers because their divisors are dense), and listing *it* in a window is Prop SI restricted to that shape
  (assessment).
- *Scope.* Families chosen from x (target-dependent, growing) are not covered. They are adaptive constructions
  (route 1).
- *Exact split compressions* (one presentation, many members) are limited to ≤ 1 relation per presentation by
  ptebarrier's Lemma R / Theorem W; Prop FF is the complementary statement for *smoothness-based* (non-exact)
  families: their members must be numerous, hence uncorrelated.

### 3.3 Toy economics versus asymptotics (D1, MODEL)
At toy sizes h ≥ √X, so t(t + 1) lands about 21 times per window at 48 bits (E1 P2) and is a far better constructor
than the sieve. In the o(1)-dropped model at ν*, κ* the expected smooth members per window of t(t + 1) is
e^{−143} at 512 bits and e^{−301} at 1024 bits; its per-output cost crosses the sieve's near 120 bits
(`theory.consecutive_product_log_hits`). This is the same toy artefact siconstruct found for MITM.

## 4. Route 3: multi-sample (Cohn–Heninger) decoders

**The decoders.** Cohn–Heninger (ANTS X, 2013; Thms 1–2, checked by codex): *partial* ACD with one exact multiple P
of the hidden divisor d = P^β and m approximate multiples with errors ≤ R recovers d when ln R < π·β^{(m+1)/m}, i.e.
ln d ≥ σ_m := t^{m/(m+1)}π^{1/(m+1)} (t = ln R), under an algebraic-independence assumption and the side condition
β²π → ∞ (their "N^{β²} ≫ 1"; for m > 1 it means "p cannot be any smaller than when m = 1"); the bound applies to the
geometric mean of the error bounds, up to (1 + o(1)); m = 1 is Howgrave-Graham/Johnson. The *general* ACD version (no exact multiple) has radius X^{C_mβ^{m/(m−1)}}. Herrmann–May
(ASIACRYPT 2008, LNCS 5350, 406–424, Thm 4) and Sarkar–Maitra (ePrint 2010/146, Thm 2: an 8× constant-factor improvement of Boneh's radius, checked by codex in
SM's text, not against Boneh's) are single-
sample or constant-factor variants.

**Prop MS (independent windows; PROVED, averaged; the decoding region granted at unit cost).** Fix 2 ≤ m ≤ t and any
y-smooth P ≥ 2 (for P = 1 the threshold is vacuous and the statement fails). Place windows W_i = [x_i, x_i + h] with x_1, …, x_m independent and uniform in [X, 2X). A *certified tuple*
is (d, z_1, …, z_m) with d | P, ln d ≥ σ_m, d | z_i ∈ W_i. Then
  E #{certified tuples with z_1 y-smooth} ≤ C·ρ(u₀)^{1−o(1)},
with C absolute and the o(1) uniform in P and m. (Summing over coordinates: ≤ C·m·ρ(u₀)^{1−o(1)} smooth outputs per
call.)
- *Proof.* Write ℓ = ln d. Every certified d has ℓ ≤ π and ℓ ≥ σ_m = t^{m/(m+1)}π^{1/(m+1)} ≥ t^{m/(m+1)}ℓ^{1/(m+1)},
  hence ℓ ≥ t and d ≥ h + 1; so each window contains ≤ 1 multiple of d, with E #multiples ≤ ((h + 1)/d)(1 + (h + d)/X)
  (for d ≤ X/m the m − 1 boundary factors multiply to ≤ e²; for d > X/m, ((h + 1)/d)^m ≤ (m(h + 1)/X)^m swamps them),
  and E #y-smooth multiples in W_1 ≤ ((h + 1)/X)·Ψ(3X/d, y) ≤ 3C₀((h + 1)/d)ρ(u₀ − ln(d/3)/ln y) (Hildebrand's
  uniform bound, as in Theorem J step 3). By independence,
  E ≤ 3C₀ Σ_{d | P, ℓ ≥ σ_m} ((h + 1)/d)^m ρ(u₀ − ln(d/3)/ln y).
  Put j₁ = mt/(m − 1) + 1.
  - *Case ℓ ≥ j₁.* Group d by j = ⌊ℓ⌋ ≥ j₁ − 1 and use #{d : ⌊ln d⌋ = j} ≤ e^{j+1}. The group contributes
    ≤ e^{1 + mt − (m − 1)j}·ρ(u₀ − (j + 1)/ln y) ≤ e^{1 − (m − 1)(j − j₁ + 1)}·ρ(u₀ − (j + 1)/ln y). Since
    ρ(v − δ)/ρ(v) ≤ exp(δ(ln(v ln v) + O(1))) and ln(u₀ ln u₀)/ln y → 0, the ρ factor grows by e^{o(1)} per unit of j
    while the exponential falls by e^{−(m−1)}. The sum is O(ρ(u₀ − (j₁ + 1)/ln y)) = ρ(u₀)^{1−o(1)}, because
    j₁ ≤ 2t + 1 = O(F_{2/3}).
  - *Case σ_m ≤ ℓ < j₁.* Then π/t ≤ (ℓ/t)^{m+1} < (m/(m − 1) + 1/t)^{m+1} ≤ 8e² (using m ≤ t), so P ≤ e^{60t},
    τ(P) ≤ exp(O(t/ln t)) = e^{o(F_{2/3})} (Nicolas–Robin), and each term is ≤ ρ(u₀ − 60t/ln y) = ρ(u₀)^{1−o(1)}. ∎
- *Reading.* The factor ((h + 1)/d)^{m−1} — the price of finding the *same* d in every window — wipes out the
  enlarged region. For m = 1 the argument fails (Johnson's threshold is needed; that is Theorem J). MEASURED (E2B):
  the smooth yield per call (distinct smooth z₁, symmetric normalisation) of m = 2, 3 is 8–92× below the m = 1
  Johnson cell, and the tuple counts match
  Σ_d((h + 1)/d)^m within 2% (C2).
- *Codex (blind §3; PROVED in the same MODEL).* Granting the general-ACD radius instead (no modulus P), the averaged
  cost is ≥ L[2/3 + 1/(3m), (m − 1)(ν/C_m)^{(m−1)/m}] per output, e.g. L[5/6, √(8ν/3)] for m = 2: worse than sieving.
  Also, the full modulus P = lcm(1..y) violates Cohn–Heninger's side condition β²π → ∞ ((ln d)²/y → 0), so it gets no
  guarantee at all.

**Lemma IH (same modulus: proportional offsets; PROVED, elementary — the standard two-sample ACD elimination
(Howgrave-Graham 2001); codex blind §3 found the same, with an exact count).** Let
k_i, k_j ≤ K and d | k_iN + r_i, d | k_jN + r_j with 0 ≤ r_i, r_j ≤ h. If d > Kh then k_jr_i = k_ir_j.
- *Proof.* d | k_j(k_iN + r_i) − k_i(k_jN + r_j) = k_jr_i − k_ir_j, and |k_jr_i − k_ir_j| ≤ Kh < d. ∎
- *Consequence (codex).* With g = gcd(k_1..k_m), k_i = gλ_i: r_i = λ_i s and k_iN + r_i = λ_i(gN + s) for one
  integer s ≤ h/λ_max. All m certified integers are multiples of the single number gN + s ∈ [gN, gN + h/λ_max], and
  they give the same relation (sicheck Lemma I). The joint certificate is a single-window certificate in disguise.
  Codex's exact count: Q_d = (G/d)(⌊T/G⌋ + 1), G = gcd(g, d), T = ⌊h/λ_max⌋, instead of (h/d)^m.
- **Prop IH-J (yield for d > Kh; PROVED in MODEL "dependent samples carry one sample's information"; proof
  critic finding 1).** Lemma IH bounds *relations per tuple*, not *certified tuples per call*: if Cohn–Heninger's
  enlarged region σ_m were granted to this degenerate configuration, the call would certify essentially every
  smooth gN + s in [gN, gN + T] — Prop SI for free. It cannot be granted, by Cohn–Heninger's own restriction (Open Book Ser. 1 (2013) 271–293, §1.2.1: inputs made by
  "taking a_i to be a small multiple of a_1 plus an additional perturbation" are excluded, and "if a_1, …, a_m are
  related in trivial ways, then the algorithm will simply recover the corresponding relations"). Proportional samples k_iN put the m − 1
  linear forms λ_jX_i − λ_iX_j, with tiny coefficients, into the Cohn–Heninger lattice, so its algebraic-independence
  hypothesis fails; the samples carry the information of the single sample gN with error ≤ T. *Model:* such a call is
  a single-sample Johnson-class call on the window [gN, gN + ⌊h/λ_max⌋] (g ≤ K ≤ h). Theorem J(a) then gives
  ≤ ρ(u₀)^{1−o(1)} certified smooth outputs per call, averaged over N. A single-window certificate *beyond* Johnson
  would be exactly the beyond-Johnson CRT decoder of Prop RD, i.e. route 1's open problem, not a route-3 gain.
- MEASURED (E2A): 0 violations over 20 moduli × 28 pairs × 513² offsets; every pair with gcd > Kh was
  proportional (3516 per modulus), and no pair with h < gcd ≤ Kh was.

**Prop MS′ (same modulus, d ≤ Kh; PROVED in the examine-all-certified-tuples model, averaged over N).** If d ≤ Kh, the cofactor of the
first coordinate is ≥ k_1N/(Kh) and, conditional on the certificate (an event depending on N mod d only), z_1/d runs
through an arithmetic progression of difference k_1 ≤ K over a range of length ≍ X/d. So each certified tuple has
z_1 y-smooth with probability ≤ C·K·ρ(u₀ − ln(Kh)/ln y) = ρ(u₀)^{1−o(1)}, and a decoder that pays ≥ 1 per certified
tuple it examines pays ≥ ρ(u₀)^{−1+o(1)} per smooth output. *Scope (codex round 2).* The conditioning is on a fixed
certificate event (fixed d, multipliers, offsets); it is valid for decoders that examine all certified tuples (as a
lattice decoder that enumerates its roots does) or select among them using residues mod d only. A selector that uses
the quotient or the full input can emit only smooth tuples (codex's example: emit only when N + 1 is a power of two);
that is a route-1 algorithm and is not covered. (The factor K bounds smooth numbers in a progression by all smooth
numbers; it is exp(O(F_{1/3})).)

**Conclusion for route 3.** Independent windows: no gain (Prop MS). Windows of one modulus: large common divisors
force a single window, where only a Johnson-class (Theorem J) or a beyond-Johnson (route 1) certificate is
available (Lemma IH, Prop IH-J, MODEL); small ones leave the divisor-trick cofactor (Prop MS′). No finite τ. Adaptively chosen
helper windows are route 1.

## 5. Route 4: smooth norm and polynomial values near (kN)^{1/d}

- **Powers and one-parameter values** (u = t^d, t^d − a, F(t)): Lemma PL(b), exponential placement. The smoothness
  gain of t^d (ρ(u₀/d) instead of ρ(u₀)) does not help.
- **Prop NV (norm values; PROVED).** Let K be a number field and α ∈ O_K be constructed as ε∏π_i^{e_i} with π_i
  generating prime ideals 𝔭_i of norm p_i^{f_i} ≤ y and ε a unit. Then |Norm(α)| = ∏p_i^{f_ie_i}: the attainable
  norms are y-smooth integers, and Norm(α) ∈ [x, x + h] is the integer knapsack Σe_if_i ln p_i ∈ [ln x, ln(x + h)]
  over the sub-factor base of primes with a principal prime of small norm (non-principal primes only enter through
  class-group relations, a further restriction). Units do not move the norm. So **norm constructors are integer
  constructors with a restricted factor base**. Under the NFS heuristic H, a non-constructed α of norm ≍ X is smooth
  with probability ρ(u₀)^{1+o(1)}: no gain either way. E1's control t² + 1 = Norm(t + i) is even 2–10× *less* often
  smooth than random integers of its size (only primes ≡ 1 mod 4 and 2 divide it).
- **Prop LC (a degree-d analogue of the Prop SI relation; height bound PROVED, norm scale OPEN).** Take f with
  f(m) ≡ 0 (mod N), θ a root, φ: Z[θ] → Z/N. A degree-d Prop-SI relation is a constructed y-smooth α with
  φ(α) ≡ v (mod N), 0 ≤ v ≤ h (tested). Then α − v lies in the kernel ideal, of norm N, so N | Norm(α − v); if every
  embedding of α has absolute value ≤ A and α ≠ v, then N ≤ (A + h)^d (codex round 2). So the constructed element has
  height ≥ N^{1/d} − h in some embedding. This is a height constraint, not a lower bound on |Norm(α)|; the v1 claim
  that the smallest constructed side is ≍ N/(h c_K) is **withdrawn** (unit orbits and the arithmetic of generator
  residues are not uniform). Example (codex): in Q(i) with m² + 1 ≡ 0 (mod N), every α = a + bi ≠ v with
  a + bm ≡ v has |α − v| ≥ √N and |Norm α| ≥ (√N − h)² = N^{1−o(1)}: no smaller than Prop SI's u ≈ kN. In real fields
  units of norm 1 leave the norm scale open. We found no field in which the constructed side is smaller than
  N^{1−o(1)}, and no construction exploiting one (assessment).
- **NFS constants (codex blind §4; MODEL).** In the size–entropy sampling model, two-sided degree-d relations give
  SNFS (32/9)^{1/3} = 1.5263, GNFS (64/9)^{1/3} = 1.9230 and, reusing rational smoothness across fields,
  MNFS ((92 + 26√13)/27)^{1/3} = 1.9019 (re-derived exactly by codex; KNOWN: Coppersmith, J. Cryptology 6 (1993) 169–180; Lenstra–Lenstra (eds.),
  LNM 1554 (1993)). Their amortised per-relation exponents
  (0.763, 0.961, 0.951) are **not** Prop SI τ's: those relations do not output a factored smooth u ∈ [kN, kN + h].
  The SNFS per-relation exponent 0.7631 = (4/9)^{1/3} sits just above the 0.7572 threshold, a numerical coincidence.
- **Verdict.** Route 4 reproduces the known NFS constants when it abandons the Prop SI shape, and gives nothing when
  it keeps it. Target-dependent fields (K chosen from x) are route 1.

## 6. Uncertified decoder outputs (sihard's fifth item)

A Coppersmith/Howgrave-Graham lattice (Howgrave-Graham, IMA C&C 1997, LNCS 1355; CaLC 2001, LNCS 2146)
returns all integer roots of its reduced polynomials, not only the certified
ones. **E3 (pre-registered, Addendum A1)** ran 5000 calls in each of two configurations (§8.4). Uncertified roots
occur 0.04 and 0.001 times per call. In U1 (198 of them) their gcd with P sits just below the certificate threshold
(median ln gcd 18.1 vs ln T_cert 18.9; 10th percentile 16.6) and their smooth fraction (37%) is 0.70 of what their
known divisor predicts. The registered threshold is the lattice's own conservative certificate; against the *ideal*
Johnson threshold √(ln V·ln P) used by Theorem J, 156 of U1's 198 (and all 6 of U2's) are already inside Theorem J's
class (codex round 2). Only 42 roots in 5000 calls (0.008 per call) fall outside it; 10 of them are smooth against
16.8 predicted by the same cofactor model (ratio 0.60). So we see no evidence of a smoothness bias beyond the known
divisor, at toy size and with small counts; this is weaker than a demonstration of "no bias", and it says nothing
about decoders of a different kind. (Conditioning on an exact gcd also imposes divisibility exclusions on the
cofactor that the unconditional density curve ignores.)

## 7. Route 5: the softer target

- c(τ) = min_κ max(κ + τ + 2/(9κ²), 2κ) < 1.9019 iff τ < 0.757169 (PROVED; codex gives the closed form: c(τ) ≥
  τ + (3/2)(4/9)^{1/3} with equality at κ = (4/9)^{1/3} once τ ≥ 0.381).
- Per-route τ: **no finite τ obtained on any route** (no L[1/3] constructor found), see the table in §0.
- The only quantities in this note at L[1/3] scale are *upper bounds on information* (Lemma IC: τ = κ, a = 1) and the
  decoder threshold a < 1.1416 (Prop RD). They describe what a constructor would have to achieve computationally.

## 8. Experiments (pre-registered: 424ee86; Addendum A1: 01eb40d)

Commands: `.venv/bin/python -m src.sinonobliv.run_exp E1|E2A|E2B|E3`; post-hoc `python -m src.sinonobliv.posthoc`
and `posthoc e2b`. Runtimes: E1 22 s (RSS 0.69 GB), E2A 5 s, E2B 67 s, E3 12 s.

**Scorecard.**

| # | Prediction | Result |
|---|---|---|
| P1 | R = joint/independence ∈ [0.67, 1.5] (cp2, cp3, aur; pred ≥ 30) | **FAILED in 4 of 13 cells** (cp2 56b y=2^8: 0.57; cp3 48b, 56b, 60b y=2^8: 0.61, 0.54, 0.46); passed in 9 |
| P2 | hit law in I_k: pooled ratio ∈ [0.95, 1.05] (pred ≥ 200) | **held**: cp2 0.9998 (40 b), 1.0008 (48 b); cp3 cells below 200 (1.02, 0.90; reported only) |
| P3 | Lemma IH: 0 violations | **held**: 0 in 20 moduli |
| P4 | tuple count ratio ∈ [0.85, 1.15] (pred ≥ 50) | **FAILED in 1 of 10 cells** (C1, m = 1 Johnson: 1.23); others 0.99–1.10 |
| P5 | z₁-smooth ratio ∈ [0.6, 1.6] (≥ 30 tuples) | **FAILED in 1 of 8 eligible cells** (C1 Johnson: 1.63); others 0.74–1.11 |
| P6 | CH (m = 2, 3) per-window smooth yield below Johnson (m = 1) | **held** in both configs |
| P7 | uncertified roots per call ≤ 0.5 | **held**: 0.040 (U1), 0.0012 (U2) |
| P8 | uncertified smooth / known-divisor model ∈ [0.33, 3] (≥ 30 roots) | **held** (U1: 0.70); U2 not evaluated (6 roots) |
| P9 | certificate: 0 violations | **held** (101 and 11 certified points, all roots) |

### 8.1 E1 (identity families)

| family | bits | y | members | smooth members | indep. pred | R | gain vs random | local factor (post-hoc) | R corrected |
|---|---|---|---|---|---|---|---|---|---|
| t(t+1) | 40 | 2^8 | 65536 | 624 | 820.0 | 0.761 | 29.2 | 0.781 | 0.975 |
| t(t+1) | 40 | 2^10 | 65536 | 4254 | 4712.3 | 0.903 | 15.6 | 0.902 | 1.001 |
| t(t+1) | 48 | 2^8 | 1048576 | 1237 | 1823.8 | 0.678 | 61.8 | 0.694 | 0.977 |
| t(t+1) | 48 | 2^10 | 1048576 | 16782 | 19351.9 | 0.867 | 30.8 | 0.870 | 0.996 |
| t(t+1) | 56 | 2^8 | 2097152 | 245 | 427.3 | **0.573** | 170.6 | 0.618 | 0.928 |
| t(t+1) | 56 | 2^10 | 2097152 | 6906 | 8486.3 | 0.814 | 53.8 | 0.827 | 0.984 |
| t(t+1)(t+2) | 48 | 2^8 | 65536 | 543 | 887.2 | **0.612** | 654 | 0.675 | 0.907 |
| t(t+1)(t+2) | 48 | 2^10 | 65536 | 5238 | 6011.1 | 0.871 | 221 | 0.884 | 0.986 |
| t(t+1)(t+2) | 56 | 2^8 | 416127 | 554 | 1031.0 | **0.537** | (4193)* | 0.591 | 0.909 |
| t(t+1)(t+2) | 56 | 2^10 | 416127 | 9830 | 12238.3 | 0.803 | 683 | 0.823 | 0.976 |
| t(t+1)(t+2) | 60 | 2^8 | 1048576 | 477 | 1044.4 | **0.457** | (1739)* | 0.547 | 0.835 |
| t(t+1)(t+2) | 60 | 2^10 | 1048576 | 11974 | 15967.8 | 0.750 | 1063 | 0.784 | 0.957 |
| 4t⁴+1 | 40 | 2^10 | 2175 | 25 | 30.3 | 0.824 | 7.9 | — | — |
| 4t⁴+1 | 46 | 2^10 | 6144 | 23 | 27.8 | (0.83) | 12.7 | — | — |
| 4t⁴+1 | 52 | 2^10 | 17379 | 20 | 19.6 | (1.02) | 21.7 | — | — |
| t²+1 | 40/48/56 | 2^10 | — | 123/226/43 | — | (control) | 0.45/0.41/0.32 | — | — |

(*) control density from a curve with empty cells: unreliable. Parenthesised R: pred < 30, not registered. 4t⁴ + 1
at y = 2^8 had 0–2 expected smooth members (not evaluated).

- **P1 failure, diagnosed (POST-HOC, `posthoc.py`).** Consecutive pieces cannot share a prime, while smooth integers
  are rich in small primes. A local model (a smooth integer is ≡ 0 mod p with the measured frequency a_p, uniformly
  in the other classes; independent across p) predicts the joint/independent ratio
  ∏_p p^{r−1}Σ_c∏_jq_p(c + j). It accounts for most of the deficit: R/local factor = 0.84–1.00 in all 12 cp cells (residual 0–17%, largest at
  y = 2^8). This is
  the r-piece analogue of sicheck's coprimality penalty (Lemma CP). At the Prop SI scale a_p → p^{−α} with α → 1 and
  the factor is a constant (L^{o(1)}); it does not affect the asymptotic gain.
- The gain grows with size at fixed y, as r^{u} predicts; its size is irrelevant to the barrier, which is placement.

### 8.2 E2A (Lemma IH): 20 moduli (32 and 40 bits), K = 8, h = 2^9. 0 violations; 3516 proportional pairs with
gcd > Kh per modulus (the inherited events, the same count for every N, as Lemma IH predicts); 5455–6039 pairs with
h < gcd ≤ Kh per modulus, none proportional.

### 8.3 E2B (independent windows)

| cfg | cell | σ (d ≥) | calls | tuples / pred | z₁-smooth tuples / model | distinct smooth z₁ per window |
|---|---|---|---|---|---|---|
| C1 (2^36, h 2^10, lcm(1..32), y 2^8) | J, m = 1 | 3.4e6 | 20000 | 11417 / 9277 = **1.23** | 5930 / 3634 = **1.63** | 0.0321 |
| C1 | CH m = 2 | 1.1e5 | 20000 | 1012 / 919 = 1.10 | 134 / 121 = 1.11 | 0.00205 |
| C1 | CH m = 3 | 2.7e4 | 20000 | 287 / 265 = 1.08 | 25 / 21.6 (n/a) | 0.00030 |
| C1 | bare m = 2 | 2050 | 2000 | 86496 / 86793 = 1.00 | 3080 / 3122 = 0.99 | 0.192 |
| C1 | bare m = 3 | 2050 | 2000 | 27186 / 26811 = 1.01 | 970 / 872 = 1.11 | 0.077 |
| C1 | divisor trick d = 2^10 | — | 20000 | — | — | 0.0215 |
| C2 (2^44, h 2^12, lcm(1..41), y 2^10) | J, m = 1 | 8.2e7 | 40000 | 10637 / 10641 = 1.00 | 4543 / 4317 = 1.05 | 0.0207 |
| C2 | CH m = 2 | 1.2e6 | 40000 | 736 / 723 = 1.02 | 80 / 108 = 0.74 | 0.00064 |
| C2 | CH m = 3 | 2.2e5 | 40000 | 188 / 185 = 1.02 | 9 / 16.5 (n/a) | 0.000075 |
| C2 | bare m = 2 | 8194 | 2000 | 205439 / 206918 = 0.99 | 7348 / 7618 = 0.96 | 0.494 |
| C2 | bare m = 3 | 8194 | 2000 | 64112 / 64279 = 1.00 | 2243 / 2184 = 1.03 | 0.197 |
| C2 | divisor trick d = 2^12 | — | 20000 | — | — | 0.0240 |

- **P4/P5 failures, diagnosed (POST-HOC, `posthoc.py e2b`; six fresh seeds of the C1 Johnson cell).** The m = 1
  counts are heavy-tailed: a single z with a huge gcd contributes all its divisors above the threshold. Fresh seeds
  gave tuple ratios 0.94–1.13 and z₁-smooth ratios 0.96–1.40. Codex (round 2) computed the exact finite-window
  expectations: 9277.37 tuples and 4005.03 smooth tuples per 20000 calls; the cofactor density model (3633.94)
  underestimates the latter by 10.2%. Against the exact expectation the registered smooth ratio is 1.48 and the six
  fresh seeds pool to 1.07; a single z = 128701918800 contributed 853 smooth tuples to the registered run. So the
  failures are a heavy tail plus a 10% model deficit.
- *Normalisation (proof critic).* The registered column "distinct smooth z₁ per window" divides coordinate-1 counts by
  all m windows, understating the m ≥ 2 cells by a factor m. With the symmetric normalisation (distinct smooth z₁
  per call) the CH-vs-Johnson gaps are 7.8× and 36× (C1, m = 2, 3) and 16× and 92× (C2), instead of the registered
  16–276×; the bare cells become 0.38, 0.23 (C1) and 0.99, 0.59 (C2). P6 holds either way.
- "Bare" cells grant every common divisor d > 2(h + 1) of P with no lattice guarantee: an idealised oracle that lists
  13.6–103 tuples per call. Per certified tuple its smooth rate (3.6%) matches the cofactor model; per window it is the
  only cell above the divisor trick, because it is not a decoder.
- The Johnson cell (0.021–0.032) and the divisor trick (0.022–0.024) agree, as Theorem J's sharpness predicts.

### 8.4 E3 (uncertified roots)

| cfg | calls | T_cert (median ln) | certified roots | uncertified roots | per call | uncert. smooth / model | certified smooth / model |
|---|---|---|---|---|---|---|---|
| U1 (x ~ 2^40, V 2^12, lcm(1..40), ℓ = s = 3, y 2^10) | 5000 | 18.9 | 101 | 198 | 0.040 | 73 / 104.0 = 0.70 | 47 / 72.6 = 0.65 |
| U2 (x ~ 2^48, V 2^10, lcm(1..60), ℓ = s = 4, y 2^12) | 5000 | 23.2 | 11 | 6 | 0.0012 | 3 / 4.0 (n/a) | 9 / 8.4 = 1.07 |

### 8.5 Deviations and post-hoc items
- `posthoc.py` (E1 local factor; E2B seed spread) was written after the failures it diagnoses.
- E3 was registered as Addendum A1 after E1, E2A and E2B had run, motivated by codex's blind round.
- The DensityCurve was patched before any run to keep windows local at small sizes (window [2^b, 2^b + min(W, 2^b))),
  after a design check showed the unpatched curve overestimated small-cofactor densities; no run used the old code.

## 9. GPT-6.1-sol

### 9.1 Blind round (`codex/sinonobliv_1.out`, xhigh, run before our analysis was written)
**Independent agreements.** No constructor; no finite τ obtained on any route; the fixed-tree union bound (= Prop AT); fixed
one-parameter families exponentially sparse (= Lemma PL); the r^{u} density gain and its L[2/3, o(1)] size;
independent-window multi-sample counting barrier; proportional offsets for windows of one modulus (= Lemma IH); NFS
constants not SI τ's; τ < 0.757169. Plausibility 5% (ours 1%).
**Adopted from codex.** The toy counterexample showing averages do not bound adaptive selection; the
bounded-likelihood barrier (which Lemma IC then shows real algorithms violate); Friedlander–Lagarias Thm 4 exponents
β_r, checked in the primary text; Cohn–Heninger's side condition β²π → ∞ and the general-ACD bound
L[2/3 + 1/(3m)]; the exact count Q_d; Sarkar–Maitra's 8× constant; the exact re-derivation of SNFS/GNFS/MNFS.
**Ours, beyond codex.** Prop JC (chains); the consequence of Lemma IC that no query/information barrier exists
(the lemma itself is the sieve, elementary); the decoder threshold a < 1.1416 attached to Prop RD (Prop RD itself is
Boneh's generalized CRT decoding at full smoothness); Prop FF; Prop MS for the partial-ACD region with an arbitrary
modulus P (codex treated general ACD); Prop MS′; Prop NV; experiments E1–E3. For FF, MS, MS′ and JC: no prior
statement found by us (literature critic). Lemma IH is the standard two-sample ACD elimination.

### 9.2 Adversarial round (`codex/sinonobliv_2.{md,out}`, on v1 04a1f52): MAJOR REVISION, all points adjudicated
Codex reproduced all six raw/post-hoc files exactly, ran the 12 tests, re-scored P1–P9 (same results as ours after
the P5 count fix), and confirmed Prop JC (with weights), Prop MS's inequalities, Lemma IC, the CRT agreement identity
and a* = 1.1415601395.

| # | Sev. | Finding | Adjudication |
|---|---|---|---|
| 1 | MAJOR | Prop LC's N/(h c_K) norm scale unjustified; Q(i) counterexample (norm ≥ (√N − h)²) | accepted: scale withdrawn; LC is now a PROVED height bound N ≤ (A + h)^d plus the Q(i) example; norm scale OPEN |
| 2 | MAJOR | Prop MS′ conditioning invalid for decoders that select outputs with more information | accepted: restricted to decoders examining all certified tuples or selecting by residues mod d; codex's selector example cited |
| 3 | MAJOR | Prop FF reading: L[1/3] fraction does not follow from leading constants; C1 gives no window bound | accepted: exact fraction p/((h+1)ρ) only; justification of p via per-window cost ≥ 1 |
| 4 | MAJOR | degree > dimension does not imply sparsity (st(t+1)); max multiplicity claim false | accepted: count distinct values (multiplicity removed from the proof); sparsity stated for balanced/coercive ranges |
| 5 | MAJOR | RD list size is an average, not worst-case; "purely algorithmic" overclaims | accepted: average vs worst-case separated |
| 6 | MINOR | MS: boundary factors raised to m − 1; P = 1 | accepted: boundary split (also proof critic #4); P ≥ 2 |
| 7 | MINOR | IC is a query-model statement | accepted: reworded |
| 8 | MINOR | missing CRT literature; BKR overstated | accepted: Li–Sidorenko–Nielsen, Bhandari–Harsha–Kumar–Sudan added; BKR "superpolynomial in specified regimes" |
| 9 | MAJOR | route-wide "τ = ∞" and barrier labels exceed the models | accepted: "no finite τ obtained"; NV and LC labelled reductions/height bound; §0 says no route is closed for arbitrary algorithms |
| 10 | MINOR | most "uncertified" E3 roots are inside Theorem J's ideal Johnson class | accepted: both classifications reported (42 outside, 10 smooth vs 16.8 predicted); "no hidden bias" weakened |
| 11 | MINOR | P5 has 8 eligible cells; P2 0.9998 | accepted (also proof critic) |
| 12 | MINOR | exact E2B expectation: model deficit 10.2%; registered smooth ratio vs exact 1.48, fresh seeds 1.07 | accepted: §8.3 rewritten |

## 10. Opus critics

### 10.1 Proof/cost critic (read-only): ACCEPT WITH FIXES, all applied
The critic re-checked Prop JC (both cases), AT, IC, RD, the constant a* = 1.141560 and τ* = 0.7571694, Lemma PL,
Prop FF's inequality, every inequality of Prop MS (case-2 cap 15.625·t at m = 2, ≤ 59.1t), Lemma IH (the E2A count 3516
= Σ_pairs(⌊h/λ_max⌋ + 1)), Prop MS′'s conditioning, posthoc.py's local model, and uncert.py's certificate; it re-scored
P1–P9 and found no code bug changing a number.

| # | Sev. | Finding | Fix |
|---|---|---|---|
| 1 | MAJOR | same-modulus d > Kh: Lemma IH bounds relations per tuple, not certified tuples per call; granting the CH region there would give Prop SI for free | Prop IH-J added: dependent samples carry one sample's information (MODEL), so Theorem J(a) applies; a beyond-Johnson single-window certificate is route 1 |
| 2 | MINOR | max multiplicity can be exp(Θ(x/ln x)) | multiplicity removed (count distinct) |
| 3 | MINOR | "L[1/3]^{−1} fraction", "list tends to L^{o(1)}" overclaim | exact forms stated |
| 4 | MINOR | MS boundary factor (1 + O(d/X))^{m−1} | split at d = X/m |
| 5 | MINOR | P6 column divides coordinate-1 counts by m | disclosed; symmetric gaps 7.8–92× |
| 6 | MINOR | P5: 8 eligible cells | fixed |
| 7 | MINOR | 1.63 is above all six fresh seeds | fixed (with codex's exact expectation) |
| 8 | MINOR | local factor explains most, not all, of the deficit | "most; residual 0–17%" |
| 9 | MINOR | gains quoted from unreliable control cells | 8–1100× |
| 10 | MINOR | IC(b) word-operation count | O(π(y)·x) bit operations per k |
| 11 | MINOR | "τ_info ≤ κ < 0.757" compares constants at different optima | restated via a = 1 ⇒ 1.8171 |
| 12 | MINOR | NV is a reduction, LC "no easier" is an assessment | relabelled |
| 13 | MINOR | "none reaches L[1/3, τ]" labelled PROVED | relabelled as a finding |
| 14–16 | NIT | JC needs gcd(P, D) = 1 and inherits the weighted "sketch"; theory.ms_case_bound used a different split; P2 rounding, bare tuple range, stray (1 + o(1)) | applied; ms_case_bound and its test aligned with the note |

### 10.2 Literature critic (read-only): OK WITH FIXES, all applied
The critic checked Friedlander–Lagarias Thm 4 and Cohn–Heninger Thms 1–2 in the primary texts, the Sophie Germain
identity, all NFS constants and thresholds, and uncert.py's use of Howgrave-Graham's lemma. Herrmann–May and
Sarkar–Maitra primary texts could not be fetched (UNVERIFIED by the critic; codex had checked them).

| # | Sev. | Finding | Fix |
|---|---|---|---|
| 1 | MAJOR (attribution) | Prop RD is Boneh's generalized CRT decoding (STOC 2000; JCSS 2002 §§3–4) at full smoothness; far-below-Johnson reading already in pilist §2.3 | credited; removed from "ours"; gcd form and "above capacity" stated |
| 2 | MAJOR (outdated) | plain RS codes now have capacity list decoders (Brakensiek et al., arXiv 2609.08005; Jeronimo 2609.05870) | bullet replaced; cost blow-up as agreement fraction ≈ 1/π(y) → 0 noted |
| 3 | MEDIUM | Prop IH-J's model is Cohn–Heninger's own restriction (§1.2.1) | quoted and cited |
| 4–6 | MINOR | BKR precise statement; GRS is √2 short of Johnson; number-field codes only to Johnson (Biasse–Quintin); CH call multi-input ACD the integer Parvaresh–Vardy | applied |
| 7 | MINOR | NFS constants credited only to codex | Coppersmith 1993; LNM 1554 |
| 8 | MINOR | smooth consecutive/polynomial-value literature missing | Martin 2002; Teräväinen 2018; Conrey–Holmstrom–McLaughlin 2013; Costello–Meyer–Naehrig 2021; Mulder–Sterner–van Woerden |
| 9 | MINOR | FL: r = 2 also Balog–Sárközy; FL's well-spacing remark makes their family an equality case of Lemma PL(a) | added |
| 10–12 | MINOR | CH side condition wording, geometric mean, (1 + o(1)); HM/SM venues and provenance of SM's 8×; Howgrave-Graham citations | applied |
| 13 | MINOR | novelty: IH is the standard ACD elimination, IC is the sieve; "no prior statement found by us" only for FF, MS, MS′, JC | applied |
| 14 | MINOR | DIVISOR attribution | sihard's form used |

## 11. Reproduction
```bash
.venv/bin/python -m pytest -q tests/test_sinonobliv.py
.venv/bin/python -m src.sinonobliv.run_exp E1     # ~25 s, 0.7 GB
.venv/bin/python -m src.sinonobliv.run_exp E2A    # ~5 s
.venv/bin/python -m src.sinonobliv.run_exp E2B    # ~70 s
.venv/bin/python -m src.sinonobliv.run_exp E3     # ~12 s
.venv/bin/python -m src.sinonobliv.posthoc        # POST-HOC E1 local factor
.venv/bin/python -m src.sinonobliv.posthoc e2b    # POST-HOC E2B seed spread
.venv/bin/python -m src.sinonobliv.theory         # constants
```
