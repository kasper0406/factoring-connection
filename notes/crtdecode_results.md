# crtdecode: adaptive CRT decoding far below the Johnson radius (sinonobliv route 1)

Agent `crtdecode`, phase 23, 2026-10-22. Version 1 (draft before the adversarial round and the critics).

**Files.**
- Code: `src/crtdecode/`.
  - `core.py`: exact window profiles, prime powers, smooth ranges, the exact idealised Johnson certificate.
  - `decoders.py`: the operational toy decoders (DT, DT8, PH-θ, RT) and the information decoder of Prop DT.
  - `theory.py`: c(a), a*, the regime dictionary.
  - `run_exp.py`, `analyze.py`: E1, E2 and their scoring.
- Tests: `tests/test_crtdecode.py` (9 pass).
- Pre-registration: `notes/crtdecode_prereg.md`, commit 7d07797, made before any code existed.
- Raw data: `experiments/raw/crtdecode_E{1,2}.jsonl`.
- GPT-6.1-sol: blind round `codex/crtdecode_1.{md,out}` (§5.1).

**Labels.** KNOWN (cited) · PROVED (proof here) · PROVED in MODEL (exact consequence of a stated model) ·
HEURISTIC (named assumption) · MEASURED · POST-HOC · OPEN · UNVERIFIED · assessment. "No prior statement found by
us" is the strongest novelty claim made.

**Notation** (as in sihard and sinonobliv). X ≈ kN, x = ln X, F_α = x^α(ln x)^{1−α}, L[α, c] = exp((c + o(1))F_α).
y = L[1/3, κ], h = L[2/3, ν], V = h + 1, u₀ = x/ln y, ρ the Dickman function, ρ(u₀) = L[2/3, −2/(3κ)].
n := π(y). κ* = (2/9)^{1/3}, ν* = 2/(3κ*). A decoder of cost n^a per window beats MNFS iff a < a* = 1.14156.

---

## 0. Verdict

**Plain language.**
- **No decoder with a < 1.1416 was found, and no barrier for arbitrary algorithms was proved.** The best adaptive
  decoder we could build or analyse is still the oblivious divisor trick, at L[2/3] per output.
- **Where the problem sits.** Prop SI is list decoding of a Chinese-remainder code exactly *at* that code's
  capacity and far *below* its Johnson radius.
  - The code has n = π(y) prime coordinates. Its "dimension" is k = ln h/ln y, and Prop SI needs agreement in
    t = u₀ coordinates. At Prop SI's optimum, t = k²/2 and n ≈ k^k.
  - A Johnson-radius decoder needs agreement ≈ √(kn) = n^{1/2+o(1)} coordinates, not k²/2.
  - The average list size is (h + 1)ρ(u₀) = L[2/3, ν − 2/(3κ)]. Prop SI's optimum ν ↓ 2/(3κ) is exactly where it
    crosses 1. That is the code-specific capacity. The textbook capacity formulas (large alphabet, or a uniform q-ary
    alphabet) put the threshold elsewhere, because the primes have very different sizes.
  - Brute force over messages (the sieve) costs h = n^{k(1+o(1))}. The target is n^{1.14}.
  - No CRT, number-field or Reed–Solomon decoder we know of works in this regime. That includes the 2026
    capacity-achieving RS decoders: in their published bounds, the derivative order needed here exceeds the message
    dimension.
- **The natural adaptive model cannot carry a barrier.** Consider a decision tree whose queries are decoder calls,
  chosen adaptively, each at unit cost. It solves Prop SI with depth ≤ 4π(y) + 1 on average, i.e. a ≤ 1 < a*
  (Prop DT, PROVED; E2 confirms it exactly).
  - The tree's leaf computation is the sieve.
  - A comparison tree does it in depth O(x).
  - So depth and query counts are useless as cost measures here. Only the size of the call repertoire, or real
    computation, can carry a lower bound.
- **What adaptivity provably cannot do, in restricted classes:**
  - Zooming into a divisor certified by an earlier call never finds anything that call did not already certify
    (Lemma GC).
  - A fixed pointer schedule with adaptively chosen tests yields at most ρ(u₀)^{1−o(1)} smooth outputs per pointer
    (Prop PT).
  - A Johnson call that certifies full smoothness by agreement alone, without a cofactor test, finds a smooth
    integer in a window with probability at most (h + 1)X^{−3/4+o(1)} (Prop SG). Pure decoding is therefore
    exponentially worse than the divisor trick.
  - Fixed repertoires of any shape need ρ(u₀)^{−1+o(1)} calls (Prop AT, sinonobliv).
- **A barrier for arbitrary algorithms is out of reach.** An unconditional lower bound of n^{a} for any fixed
  a > 0 would imply P ≠ NP (Prop PNP).
  - The remaining question is purely computational. Given the residue table, which costs n^{1+o(1)}, find a point
    of [0, h] covered by weight ≥ x of the progressions j ≡ −X (mod q^e). This is the problem HEAVY-POINT.
  - With an explicit prime list instead of all primes ≤ y it is NP-hard (DIVISOR).
- **Experiments (pre-registered, all five predictions held).**
  - E1: the natural adaptive aggregation rule tests first the positions hit by several smooth divisors. It
    concentrates smoothness 4.4–5.5× among tested positions. But it costs 3.5–3.9× *more* per first smooth output
    than the oblivious divisor trick, because building collisions costs pointers.
  - E2: the information decoder (Prop DT) reconstructs every window exactly with 2.25–2.27·π(y) calls. Its leaf
    computation equals the sieve, incidence for incidence.
- **Assessment.** We put the probability that a CRT decoder with a < 1.1416 exists at about 1%; codex's blind
  estimate is 10% for an average-case version and 3% for all promised windows. In decoding terms the problem is a
  "planted, at-capacity, sub-Johnson" regime. The RS analogue of this regime (split polynomials in a coset) is tied to
  discrete logarithms by Cheng–Wan, mirroring sihard's Prop DLI.

**Main results.**

| Result | Status |
|---|---|
| Prop F1: exact formulation (C_y, valuation agreement, threshold exact for h < X) | PROVED (= Prop RD/Boneh; codex) |
| Prop F2: regime: agreement fraction ≈ 1/n; Johnson ratio L[1/3, −κ/2]; average list L[2/3, ν − 2/(3κ)], so the code-specific capacity threshold is ν = 2/(3κ); the large-alphabet (α = R) and uniform q-ary (ν = 1/κ) formulas misplace it | PROVED (averaged); average list via Hildebrand–Tenenbaum |
| Prop F3: dictionary (n, k, t): t/k² → κ/ν², ln n/(k ln k) → 3κ²/ν; at (κ*, ν*): t = k²/2, n = k^{k(1+o(1))}; sieve = n^{k(1+o(1))} | PROVED (algebra) |
| Prop F4: smooth ⇔ minimum-weight CRT codeword; Prop SI = minimum-weight codeword in the coset X + [0, h]; RS analogue = split polynomial in g + F_q[T]_{<k} (Cheng–Wan's BDD regime) | PROVED (F4); analogue KNOWN in substance, UNVERIFIED details |
| Prop DT: decision trees of idealised decoder calls: depth ≤ 2π(y) + 1 + Σ_q m_q, E ≤ 4π(y) + 1 + o(π(y)), worst case π(y)^{1+o(1)}; leaf = sieve | PROVED; E2 |
| Lemma NU: a comparison tree of depth ≤ log₂ Ψ + 1 = O(x) solves Prop SI nonuniformly | PROVED (trivial) |
| Lemma GC: zoom chains stay inside the first call's certified set | PROVED (trivial) |
| Prop PT: oblivious pointer schedules (d ≥ V/2), adaptive tests: ≤ T·ρ(u₀)^{1−o(1)} smooth outputs | PROVED (averaged; via Prop JC's estimate) |
| Prop SG: full-agreement Johnson certification: ≤ (h + 1)·X^{−3/4+o(1)} smooth outputs per fixed call | PROVED (averaged) |
| Prop PNP: an unconditional n^{a} lower bound (a > 0 fixed) implies P ≠ NP | PROVED (standard) |
| Lemma SC: exact interval counting of smooth numbers in time C gives Prop SI in C·log₂(h + 2) + poly | PROVED (trivial) |
| No decoder with a < 1.1416 among the mechanisms in §3.2 | finding, not a theorem |

---

## 1. Formalisation and regime map (Question 1)

### 1.1 The code and the problem

**Definition.** For X ∈ [X₀, 2X₀) let e_q = ⌊ln(3X₀)/ln q⌋ and P = ∏_{q≤y} q^{e_q}. The code C_y has
- messages j ∈ [0, h] (an integer interval, size V);
- for each prime q ≤ y, one coordinate holding j mod q^{e_q};
- the *valuation agreement* at q between two words a, b given by min(v_q(a − b), e_q)·ln q (Boneh's gcd form);
- received word w = −X mod P.

The agreement of message j is A(j) = ln gcd(X + j, P).

**Prop F1 (PROVED; Prop RD of sinonobliv, Boneh's generalized CRT problem at full smoothness; codex §1).**
X + j is y-smooth iff A(j) = ln(X + j). If h < X this is equivalent to A(j) ≥ ln X: a proper divisor of
X + j < 2X is smaller than X. The factorisation is read off the agreeing prime-power prefixes. Prop SI asks for
*one* element of the list 𝓛_X = {j : A(j) = ln(X + j)}, with its factorisation.

*Codex's metric remark (PROVED).* The coordinates q, q², … of one prime are nested. Treating them as independent
symbols gives the wrong rate, and replacing a tower by one equality test loses the partial-valuation information. We
therefore measure rates in nats with the compressed coordinates.

### 1.2 Rate, agreement, Johnson radius, capacity

**Prop F2 (PROVED; averaged statements over X ∈ [X₀, 2X₀)).** Write W = ln P = (1 + o(1))·n·x, b = ln V = νF_{2/3}.
- (a) *Rate and agreement.* The rate is R = b/W and the required agreement fraction is α = x/W ∼ 1/n. So
  α/R ∼ x/b = (x/ln x)^{1/3}/ν → ∞, but only polylogarithmically, while α/√R ∼ √(x/(nb)) = L[1/3, −κ/2] → 0.
  - The Johnson threshold for C_y (two distinct messages agree in weight ≤ ln h) is A_J ≈ √(bW), and
    A_J/x = L[1/3, κ/2] (codex; Theorem J's threshold).
  - The squarefree code (one coordinate mod q per prime, W_sf = θ(y) ∼ y) shows the same separation
    R_sf ≪ α_sf ≪ √R_sf.
- (b) *Average list size (codex §1; Hildebrand–Tenenbaum).*
  E|𝓛_X| = ((h + 1)/X₀)·(Ψ(2X₀, y) − Ψ(X₀, y)) + O(h²/X₀) = (h + 1)ρ(u₀)^{1+o(1)} = L[2/3, ν − 2/(3κ)].
  - For ν < 2/(3κ), most windows have an empty list.
  - For ν > 2/(3κ), the mean list size L[2/3, ν − 2/(3κ)] exceeds every power of n. So listing all solutions at
    cost n^{a} per window is impossible on average (codex's output-size obstruction).
  - **The threshold ν = 2/(3κ) where the average list crosses 1 is the code-specific (average-case) list-decoding
    capacity of C_y, and Prop SI's optimum ν ↓ 2/(3κ) sits exactly on it.** By averaging, some received word has a
    list at least the mean, so C_y's worst-case list-decoding radius for lists of size L[2/3]^{o(1)} is at most this
    threshold. Whether worst-case lists stay small *at* the threshold is OPEN (BKR-type phenomena are not excluded,
    sinonobliv §2.3).
- (c) *The textbook formulas misplace the threshold.*
  - Large-alphabet capacity α = R would place Prop SI far inside the decodable region, by the polylog factor in (a).
  - The uniform q-ary formula R = 1 − H_q(1 − α), with alphabet q ≈ y on n coordinates, gives expected list
    V·C(n, t)·y^{−t} = V·exp(−u₀(ln u₀ + ln ln y − 1 + o(1))). Its threshold is ν = 1/κ = (2/(3κ) + 1/(3κ)), not
    2/(3κ).
  - The difference u₀(ln ln y − ln ln u₀) = (1/(3κ) + o(1))F_{2/3} comes from the small primes. Their coordinates
    have small alphabets, so agreement there is cheap, and the smooth-number count Ψ, not a q-ary entropy, governs
    the list.
  - Numerically (`theory.dictionary`, o(1) dropped), R/R_cap(q-ary) = 0.80, 0.79, 0.77, 0.76 at 512, 1024, 2048
    and 4096 bits, drifting towards νκ = 2/3.
- *Proof of (c).* Stirling and de Bruijn's expansion ln ρ(u) = −u(ln u + ln ln u − 1 + o(1)), with
  ln ln y = (1/3 + o(1)) ln x and ln ln u₀ = ln ln x + O(1). ∎

### 1.3 The Reed–Solomon dictionary

**Prop F3 (PROVED; algebra, o(1) terms explicit).** For the squarefree prime code (n = π(y) coordinates of alphabet
≈ y), put k := ln V/ln y (the number of "symbols" in a message) and t := x/ln y = u₀ (the number of agreeing primes a
squarefree smooth X + j needs). Then:
- t/k² = κ/ν²·(1 + o(1)) and ln n/(k ln k) = (3κ²/ν)(1 + o(1));
- at (κ*, ν*): t/k² → 9κ*³/4 = 1/2 and ln n/(k ln k) → 9κ*³/2 = 1;
- brute force over messages costs V = y^k = n^{k(1+o(1))}; the sieve is this brute force;
- the Johnson threshold needs √(kn) = n^{1/2+o(1)} agreeing primes.

*Proof.* k = (ν/κ)(x/ln x)^{1/3}, t = (1/κ)(x/ln x)^{2/3}, ln n = ln y − ln ln y = κF_{1/3}(1 + o(1)), and
ln k = (1/3 + o(1)) ln x. ∎

Numerically (o(1) dropped): at 1024 bits n ≈ 8.9·10⁶, k ≈ 8.7, t ≈ 37.5, √(kn) ≈ 8800. At 2048 bits
n ≈ 4.7·10⁹, k ≈ 10.5, t ≈ 55.6, √(kn) ≈ 2.2·10⁵.

**Reading.** Prop SI asks to decode, in time n^{1.14}, a code with dimension k ≈ ln n/ln ln n, from agreement in
t ≈ k²/2 of n ≈ k^k coordinates (polylogarithmic agreement). Brute force over messages, or over k-subsets of
coordinates, costs n^{Θ(k)}: quasi-polynomial in n, and that is the sieve. In RS language this is polynomial
reconstruction far below √(kn) and above capacity, in the large-alphabet sense. That regime is the subject of
cryptographic hardness assumptions: Kiayias–Yung's polynomial reconstruction assumption (ICALP 2002) and Naor–Pinkas'
noisy polynomial interpolation. Bleichenbacher–Nguyen (EUROCRYPT 2000) broke some parameter sets of the latter with
lattices, and also treat noisy CRT. (Citations from memory; UNVERIFIED pending the literature critic.)

### 1.4 Minimum-weight codewords and the function-field analogue

**Prop F4 (PROVED).**
- For 1 ≤ m ≤ 3X₀, the weighted weight of c(m) is ln P − ln gcd(m, P) ≥ ln P − ln m. Equality holds iff m | P, iff m
  is y-smooth (v_q(m) ≤ e_q holds automatically for m ≤ 3X₀).
- So **Prop SI = find a minimum-weight word in the affine coset {c(X + j) : 0 ≤ j ≤ h}**: a Singleton-type equality
  case.
- RS analogue. With evaluation set A ⊆ F_q, the evaluation vector of f has weight ≥ |A| − deg f, with equality iff f
  splits into distinct linear factors over A. The analogue is: find f in the coset g + F_q[T]_{<k} (deg g = t) that
  splits over A. Equivalently, bounded-distance decode the received word −g|_A to RS[|A|, k] at agreement t.
- The dictionary of §1.3 is the classical smooth-integer ↔ split-polynomial analogy, which we use only as an
  analogy.

**Literature on the analogue (UNVERIFIED pending the literature critic).**
- Cheng–Wan (FOCS 2004; SIAM J. Comput. 2007) study exactly this regime. Received words are evaluations of a
  polynomial of degree slightly above k. They show that bounded-distance decoding of RS codes there is at least as
  hard as discrete logarithms in an extension field.
- Gandikota–Ghazi–Grigorescu (codex) extend NP-hardness of RS decoding to radii N − K − d with d ≲ log N/log log N,
  at dimension about N/2 and exponential field size. That is a different regime from k ≈ ln n/ln ln n, t ≈ k²/2.
- The DLP link mirrors sihard's Prop DLI: a modulus-generic Prop SI constructor gives prime-field discrete logarithms
  in L[1/3, 1.232].
- In small characteristic, Frobenius identities (x^q − x = ∏(x − a)) do produce split polynomials in structured
  cosets. That is the engine of the quasi-polynomial DLP descent (Barbulescu–Gaudry–Joux–Thomé 2014). We know of no
  integer analogue of that identity (assessment).

### 1.5 Literature comparison

| Result | What it decodes | Here |
|---|---|---|
| Goldreich–Ron–Sudan (IEEE IT 2000) | CRT codes, agreement ≳ √(2k_c n log p_n/log p_1) | Johnson-type, √2 short; L[1/3, κ/2] too high |
| Boneh (STOC 2000; JCSS 2002) | removes the √2; smooth integers in intervals (gcd form) | Johnson; vacuous at full smoothness (sihard §4) |
| Guruswami–Sahai–Sudan (FOCS 2000) | weighted ("soft") CRT, √(Q·ln V) | Johnson, weighted; Theorem J covers it |
| Coppersmith/Howgrave-Graham; CHHS (ASIACRYPT 2016) | univariate auxiliary polynomials; CHHS prove optimality of the method's bound (capacity theory) | the certificate family stops at Johnson; arbitrary algorithms not covered |
| CRT/number-field capacity (Guruswami 2003 constructions; Cohn–Heninger multi-sample ACD) | change the code or need several correlated samples | C_y is fixed; one sample (sinonobliv route 3) |
| Brakensiek–Chen–Putterman–Zhang–Zheng, arXiv 2609.08005 (Thm 4.1 as read by codex; UNVERIFIED by us) | plain RS over prime fields to capacity: agreement ε, slack θ, derivative order d = ⌈ε^{−3/θ}⌉, needs k_RS > d, time Q^{O(ε^{−12/θ})} | here ε ≈ L[1/3, −κ], so d = L[1/3, 3κ/θ] ≫ k_RS = poly(x): hypothesis fails; no transfer theorem to CRT either |
| Jeronimo, arXiv 2609.05870 (UNVERIFIED) | Q^{O_γ(1)} at fixed additive slack γ | the additive slack here → 0 |
| Ben-Sasson–Kopparty–Radhakrishnan (IEEE IT 2010) | superpolynomial RS lists beyond Johnson for special parameters | worst-case lists; not finding one element |
| Guruswami–Vardy (2005); GGG | NP-hardness of ML / BDD RS decoding for constructed instances | different parameters, no CRT transfer |
| Cheng–Wan (2004/2007) | RS BDD at agreement slightly above k ≥ DLP | the function-field analogue of Prop SI (§1.4) |

*Where the 2026 RS decoder fails (answering the brief).* It is a polynomial-ring decoder with hidden Hasse
derivatives. C_y has no global derivation, and the nested prime-power coordinates are not derivatives. Even granting
a transfer, its derivative order exceeds the message dimension in this regime, so its theorem does not apply (codex
§1).

## 2. Adaptive decoding (Question 2)

### 2.1 Decision trees of decoder calls: the natural model is information-trivial

The model, as the brief suggests: each query is one decoder call, chosen as any function of earlier answers, at unit
cost. X enters only through answers. The calls are idealised unweighted Johnson certificates (Theorem J's class
𝒥): modulus M (y-smooth), region R = [X + ℓ, X + ℓ + |R| − 1], certified set
{z ∈ R : (ln gcd(z, M))² ≥ ln|R|·ln M}. The answer lists each certified z with gcd(z, M).

**Prop DT (PROVED).** Let X ∈ [2^{b−1}, 2^b) with b public, and h < 2^{b−1}. The following tree outputs a y-smooth
z ∈ [X, X + h] with its factorisation iff one exists.
1. For each q ≤ y, with a_q = ⌊log_q V⌋ ≥ 1, call (q^{a_q}, [X, X + q^{a_q} − 1]). The region has length equal to
   the modulus, so exactly the unique multiple is certified. This reveals X mod q^{a_q}.
2. For each q and e = a_q + 1, a_q + 2, …: call (q^e, window). Stop at the first e at which no point has gcd q^e.
   Since q^e > V, every point with v_q ≥ e is certified.
3. Leaf: all valuations in the window are now known. X + j is smooth iff its smooth part is ≥ 2^{b−1}: a non-trivial
   cofactor is ≥ y + 1 ≥ 3 and X + j < 3·2^{b−1}. Certify the least smooth z with one call (M = z, window).

Depth D(X) = 2π(y) + 1 + Σ_q m_q, where m_q is the number of levels e > a_q present. Always D ≤ π(y)(2 + log₂(3X)) + 1
= π(y)^{1+o(1)}, and E_X D ≤ 4π(y) + 1 + O(π(y)·b·V/2^{b−1}).

*Proof.* Correctness is by construction. For the expectation:
Pr[some j ∈ [0, h] has q^e | X + j] ≤ V/q^e + V/2^{b−1}, and Σ_{e>a_q} V/q^e < q/(q − 1) ≤ 2 because
q^{a_q+1} > V. ∎

- **So in depth the threshold is beaten informationally: a ≤ 1.** With regions not restricted to the window, one
  call (M = P, |R| = P) reveals X mod P, and the depth is 2.
- The leaf computation is the sieve on known residues. In E2 it equals the sieve's incidence count exactly
  (≈ 2.9·h).
- **Lemma NU (PROVED).** For the family X ∈ [X₀, 2X₀), binary search among the S smooth integers of
  [X₀, 2X₀ + h] is a comparison tree ("X ≥ c?") of depth ⌈log₂(S + 1)⌉ + 1 ≤ x/ln 2 + 3. It outputs a smooth
  element of the window whenever one exists. So every nonuniform depth measure is O(x).
- **Consequence.** A lower bound strong enough for Prop SI must charge the *size* of the call repertoire (Prop AT;
  codex's portfolio bound) or real computation. Codex makes the same point. A fixed b₀-ary tree of depth T has up to
  b₀^T distinct calls, so Theorem J only gives T ≥ ln(1/ρ)/ln b₀ = poly(x) (codex §3).

### 2.2 What adaptivity cannot do in restricted classes

**Lemma SR (siconstruct, restated).** Dividing a known smooth divisor D out of the window maps SMOOTH-INT(X, h, y)
to SMOOTH-INT(X/D, h/D, y). The relative width h/X is unchanged and u₀ drops by ln D/ln y. Every nested divisor or CRT
pointer step is such a self-reduction. The divisor trick is the case D ≈ h, which gains ρ(u₀ − k)/ρ(u₀) =
L[1/3, 2ν/(3κ)].

**Lemma GC (PROVED, trivial; tested).** Let a call certify z iff gcd(z, P₁) ≥ T (any monotone threshold, e.g.
Johnson's). If z₁ is certified with d = gcd(z₁, P₁), then every z in the window with d | z is certified by the same
call. *Proof.* gcd(z, P₁) ≥ d ≥ T. ∎
- *Consequence.* Boneh-style iteration that re-runs a decoder on the progression of a divisor found by an earlier
  call (zoom-in) never leaves the first call's certified set. Its smooth yield is bounded by Theorem J for that call.
  Useful adaptivity must use *sub-threshold* information.

**Prop PT (oblivious pointer schedules, adaptive tests; PROVED, averaged over X).** Fix a sequence d₁, d₂, … of
integers d_t ≥ h/2, independently of X. Consider any decoder that issues pointers in this order (paying ≥ 1 each),
may test only positions j with d_s | X + j for an issued d_s, and chooses its tests adaptively in any way. Then the
expected number of distinct y-smooth outputs after T pointers is ≤ T·ρ(u₀)^{1−o(1)}, uniformly in the schedule.
- *Proof.* Outputs lie in ∪_{t≤T}{j : d_t | X + j, X + j smooth}. Each smooth z with d | z lies in ≤ V windows, so
  E_X #{j : d | X + j, X + j smooth} ≤ (V/X₀)·Ψ(3X₀/d, y) ≤ 3C₀·(V/d)·ρ((ln(3X₀) − ln d)/ln y), by Hildebrand's
  uniform bound as in Theorem J. Summing over t ≤ T gives ≤ T·sup_{d ≥ h/2}(that).
  - For h/2 ≤ d ≤ V: V/d ≤ 3 and ρ(u₀ − ln d/ln y) = ρ(u₀)·e^{o(F_{2/3})}, by ρ(v − δ)/ρ(v) ≤ exp(δ(ln(v ln v) + O(1)))
    with δ ≤ ln V/ln y.
  - For d > V: the term is ≤ ρ(u₀)(V/d)^{1−o(1)}e^{o(F_{2/3})}, by Prop JC's case D > V. ∎
- *Scope.* This covers E1's DT, DT8 and PH-θ, whose pointer orders are fixed and whose tests are adaptive. Pointer
  orders chosen adaptively from a fixed repertoire R are covered by Prop AT (≤ |R|·ρ^{1−o(1)}). Pointer moduli
  chosen from residues of X are not covered: that is Prop DT's regime. Small moduli d ≪ h (many multiples per
  pointer) turn the schedule into a partial sieve and are not covered by the per-pointer statement.

**Prop SG (pure decoding: full agreement certificates; PROVED, averaged over X, fixed modulus).** Suppose a call
(P, window) is to certify full smoothness by agreement alone, i.e. certify z with z | P and no cofactor test. Then
Johnson's condition (ln z)² ≥ ln V·ln P forces ln P ≤ Λ_J := (ln(3X₀))²/ln V = x^{4/3−o(1)}. Every such z is supported
on S = supp(P), and θ(S) := Σ_{q∈S} ln q ≤ Λ_J. Hence
  E_X #{certifiable smooth z in the window} ≤ (V/X₀)·Ψ(3X₀, z_J) ≤ V·X^{−3/4+o(1)},
where z_J is the largest prime with θ(z_J) ≤ Λ_J, so z_J = (1 + o(1))Λ_J = (ln X)^{4/3+o(1)}.
- *Proof.* Map the i-th prime of S to the i-th prime overall. This injection on S-smooth integers does not increase
  values (Theorem J step 4(b)), so #{m ≤ 3X₀ : supp m ⊆ S} ≤ Ψ(3X₀, p_{|S|}). We have p_{|S|} ≤ z_J because the |S|
  smallest primes have θ ≤ θ(S) ≤ Λ_J. Then use Ψ(X, (ln X)^A) = X^{1−1/A+o(1)} for fixed A > 1 (KNOWN; de Bruijn),
  with A = ln z_J/ln x → 4/3. ∎
- *Reading.* A Johnson-class call reaches full agreement only with a modulus so small (ln P ≤ x^{4/3}, against
  ln P ≈ n·x for the full code) that the integers it can certify are exponentially sparse. Every Johnson-class
  success at Prop SI's scale (Theorem J's ρ^{1−o(1)}) therefore comes from partial agreement plus a cofactor test,
  i.e. from the divisor trick. "Decoding" in the strict sense contributes nothing. For a repertoire of K moduli
  multiply by K. At practical sizes A ≈ 1.22, so the exponent is ≈ −0.82, not −3/4 (`theory`, o(1) dropped).

**Prop AT (sinonobliv; codex's portfolio form).** Any adaptive use of a fixed portfolio of K certificate calls has
success probability ≤ K·ρ(u₀)^{1−o(1)}. The bound is on the portfolio size, not the running time.

**The ladder.**

| Adaptivity allowed | Bound | Status |
|---|---|---|
| Fixed schedule (Theorem J, Prop JC) | ρ^{1−o(1)} per call | PROVED |
| Fixed schedule of pointers, adaptive tests (Prop PT) | ρ^{1−o(1)} per pointer | PROVED |
| Zoom into certified divisors (Lemma GC) | inside the first call's set | PROVED |
| Adaptive selection from a fixed repertoire R (Prop AT) | ≤ \|R\|·ρ^{1−o(1)} | PROVED |
| Full agreement only (Prop SG) | V·X^{−3/4+o(1)} per call | PROVED |
| Moduli computed from X (decision tree, Prop DT) | depth π(y)^{1+o(1)}, a ≤ 1 | PROVED: no barrier possible |
| Uniform computation | OPEN; an n^{a} lower bound ⇒ P ≠ NP (Prop PNP) | — |

### 2.3 The residual computational problem

**HEAVY-POINT_y.** Input: h and the residue table r_{q^e} = −X mod q^e (q ≤ y, e ≤ e_q), of size n^{1+o(1)} words.
Output: j ∈ [0, h] with Σ_{q^e : j ≡ r_{q^e}} ln q ≥ ln(X + j).
- The table costs n^{1+o(1)} (Lemma IC). Given the table, Prop SI and HEAVY-POINT_y are the same problem.
- The sieve solves it in O(h ln ln y).
- With an explicit list of primes in place of all primes ≤ y it is NP-hard under randomized reductions (DIVISOR,
  Alon–Kilian; sihard Prop XS).
- Prop SI needs n^{1.14}, i.e. sublinear in h by a factor n^{Θ(k)}.

**Prop PNP (PROVED; standard).** Promise-Prop SI is an NP search problem. Its witness is j and the factorisation,
with ≤ log₂(2X) prime factors, each checkable in polynomial time. If P = NP it is solvable in time poly(log X), by
binary search on j through the NP predicate "[X + a, X + b] contains a y-smooth integer". A lower bound of n^{a} =
L[1/3, aκ] for any fixed a > 0, which is superpolynomial in log X, therefore implies P ≠ NP.
- Conversely, P ≠ NP gives no lower bound: at this scale the problem is not NP-hard unless NP = coNP, under sihard's
  Prop NC hypotheses.
- So barriers must be either model-restricted (§2.2) or conditional on a stronger structured assumption.

**Lemma SC (PROVED, trivial).** If Ψ(b, y) − Ψ(a − 1, y) can be computed exactly in time C for a, b ≤ 3X, then Prop
SI costs ≤ C·⌈log₂(h + 2)⌉ plus one trial division by the primes ≤ y.
- Exact interval counting is the lattice-point count of a thin slab {ln a ≤ Σe_i ln p_i ≤ ln b} in dimension
  π(y) (Prop CV's geometry).
- Approximate counting does not help. Saddle-point methods (Hildebrand–Tenenbaum) and FPRAS-type knapsack
  counting (Dyer) give *relative* error on Ψ(b, y) ≈ bρ. The window count is a difference of two such numbers and is
  of order 1.
- We know of no exact method below L[2/3] at this scale (assessment).

## 3. Lower bound or construction (Question 3)

### 3.1 Answer

- **No construction with a < 1.1416.** The divisor trick remains best in every adaptive class we analysed.
- **No lower bound in an unrestricted model.** That is impossible without proving P ≠ NP (Prop PNP).
- **The natural adaptive model has the opposite answer: a ≤ 1 in depth** (Prop DT). It shows that Lemma IC's cheap
  information is in fact exploitable by an adaptive decoder, and that everything left is the leaf computation.
- **Barriers in restricted adaptive classes** are Lemma GC, Prop PT, Prop SG and Prop AT (§2.2).

### 3.2 Mechanisms examined and costed

| Mechanism | Cost per output | Why |
|---|---|---|
| Pointer heavy-hitters (smooth positions collect many divisor hits) | ≥ ρ^{−1+o(1)} (Prop PT); E1: 3.5–3.9× worse than DT | Heuristic count: enhancement ≈ (ω ln y)^{k} = L[1/3, Θ(1)] (the divisor-count gain), but one hit on a smooth position costs ≈ V/E pointers = L[2/3] |
| Iterated Boneh / zoom-in | ≤ first call (Lemma GC) | certified progressions are already certified |
| Support guessing + Johnson (random S of Johnson-admissible size) | exponential | Prop SG: any fixed S supports X^{1/4+o(1)} integers ≤ 3X; biased S does not help, since a smooth integer's prime mass sits near y |
| High-density subset-sum methods (BKW/Lyubashevsky, Flaxman–Przydatek, Wagner k-tree) | not applicable | combinations ≈ (X/h)·(hρ): density L[2/3, ν − 2/(3κ)] only through the window slack; the size budget forces sparse supports (≤ u₀ primes) and lists of pieces X^{1/m} too small for the k-tree (sihard's Wagner budget) |
| Exact counting + binary search (Lemma SC) | as costly as exact slab counting | no exact method known below L[2/3] |
| Hidden-derivative / multiplicity decoding | not applicable | derivative order ≫ dimension; prime-power coordinates are not derivatives (§1.5) |
| Lattices (centred CVP, Prop CV) | exact promise-CVP in dimension n in time n^{1.14} | would need near-linear exact CVP for one structured family; codex's "most promising" mechanism, but no algorithm is known |
| Generic additive model (only Σe_i ln p_i comparisons) | n^{⌈u₀/2⌉} = X^{1/2} for r-linear decision trees with generic weights (Erickson 1999, KNOWN; transfer = assessment) | the sieve escapes it by using residues; unrestricted linear decision trees are nonuniformly polynomial (Kane–Lovett–Moran 2018), as Lemma NU predicts |

### 3.3 Conditional evidence
- A CRT decoder with a < 1.1416 that works for arbitrary X is a modulus-generic Prop SI constructor (Lemma EQ). By
  sihard's Prop DLI it would also give prime-field discrete logarithms in L[1/3, max(c(a), 1.232)], with
  c(a) < 1.902 (HEURISTIC). It would therefore beat two long-standing records.
- In the function-field analogue, Cheng–Wan tie the corresponding RS decoding problem to discrete logarithms
  (§1.4).

## 4. Experiments (pre-registered: 7d07797)

Commands: `.venv/bin/python -m src.crtdecode.run_exp E1` (14 s, 59 MB RSS), `… E2` (45 s, 42 MB), scoring
`python -m src.crtdecode.analyze`.

**Scorecard.**

| # | Prediction | Result |
|---|---|---|
| P1 | capped mean(PH) ≥ 0.8·capped mean(DT), both θ, both configs | **held**: PH/DT = 3.51, 3.76 (C1); 3.80, 3.95 (C2) |
| P2 | per-test smooth rate ≥ 2× for positions with ≥ 2 hits vs 1 hit | **held**: 4.39× (C1), 5.48× (C2) |
| P3 | DT < RT < SV (capped means) | **held** in both configs |
| P4 | information decoder exact in 100% of windows | **held**: 100/100 (C1), 100/100 (C2) |
| P5 | depth ≤ 3π(y) + 1 | **held**: max 231 ≤ 292 (C1), max 400 ≤ 517 (C2) |

### 4.1 E1 (300 windows per config; all windows contain ≥ 1 smooth integer)

| cfg | decoder | capped mean cost to first | success ≤ B | pointers / window | tests / window |
|---|---|---|---|---|---|
| C1 (2^44, y 2^9, h 2^14, B 3000; 8.6 smooth per window) | DT | 111.4 | 100% | 48.1 | 63.3 |
| C1 | DT8 | 140.3 | 100% | 101.0 | 39.2 |
| C1 | PH-1.25 | 390.7 | 100% | 380.6 | 10.1 |
| C1 | PH-1.5 | 419.0 | 100% | 411.2 | 7.8 |
| C1 | RT | 1523.0 | 77.7% | — | 1523 |
| C1 | SV (incidences to first) | 5612.9 | — | — | — |
| C2 (2^52, y 2^10, h 2^16, B 6000; 13.4 smooth per window) | DT | 205.2 | 100% | 86.6 | 118.5 |
| C2 | DT8 | 251.0 | 100% | 179.2 | 71.7 |
| C2 | PH-1.25 | 779.2 | 99.7% | 763.9 | 15.4 |
| C2 | PH-1.5 | 809.9 | 99.7% | 799.2 | 10.7 |
| C2 | RT | 3460.7 | 71.7% | — | 3461 |
| C2 | SV | 15459.8 | — | — | — |

- Per-test smooth rate by hits at test time (P2). C1: 600/5363 = 11.2% (≥ 2 hits) vs 300/11773 = 2.5% (1 hit). C2:
  598/7823 = 7.6% vs 300/21522 = 1.4%.
- The "1 hit" group comes from DT8, which tests at first landing. The "≥ 2" group comes from PH. Each decoder stops
  at its first success, so each window contributes ≤ 1 success per decoder. These are per-test rates (disclosed
  reading of the registered "pooled").
- **Reading.** Aggregation works as designed: positions hit twice are 4–5× likelier to be smooth, which is the
  divisor-count gain that makes smooth numbers heavy hitters. But PH spends ≈ 98% of its budget on pointers to create
  those collisions, and loses 3.5–3.9× to the oblivious divisor trick.
- Toy sizes cannot show asymptotics. Prop PT says the gap cannot reverse beyond an L^{o(1)} factor for any oblivious
  pointer schedule.
- DT8 loses 1.2× to DT, because moduli above h often have no multiple in the window.

### 4.2 E2 (100 windows per config)

| cfg | π(y) | calls (mean / max) | step-1 calls | exact | leaf incidences / h | leaf = sieve |
|---|---|---|---|---|---|---|
| C1 | 97 | 220.1 / 231 | 97 | 100/100 | 2.87 | yes (all windows) |
| C2 | 172 | 386.8 / 400 | 172 | 100/100 | 2.98 | yes (all windows) |

- Depth ≈ 2.25–2.27·π(y): step 1 contributes π(y) calls (one per prime: no two of the prime powers q^{a_q} fit
  under V together), step 2 contributes π(y) + Σ m_q, and one final certification call.
- The tree certified the least smooth integer in every window.
- The leaf's incidence count equals the full window sieve's count exactly. Prop DT's informational a = 1 is bought by
  a leaf that *is* the sieve.

### 4.3 Deviations
- The tested-position log is stored aggregated by hit count (`by_hits`), decided before the first run, for file size.
  It is not a registered quantity.
- P2's "pooled" is read as per-test rates by hits at test time, as described in §4.1.
- No other deviations. No post-hoc analysis was used for any registered claim.

## 5. GPT-6.1-sol

### 5.1 Blind round (`codex/crtdecode_1.out`, xhigh; prompt written before our analysis)

**Independent agreements.**
- The exact formulation (= F1).
- The Johnson ratio L[1/3, κ/2], with α between R and √R.
- The average list size L[2/3, ν − 2/(3κ)].
- Query counting collapses under free computation (= the content of Prop DT and Lemma IC).
- Theorem J gives only portfolio/size bounds (= Prop AT), and depth T ≥ ln(1/ρ)/ln b₀ is poly(x).
- No transfer from BKR, GR, GV, GGG or subset product.
- a* = 1.14156 (codex rounded MNFS to 1.9019, giving 1.1415878; with the exact MNFS it is 1.1415601).
- No construction, no computational lower bound.
- Codex's probability: 10% (average-case formulation) and 3% (all promised windows). Ours is 1%.

**Adopted from codex.**
- The nested-coordinate metric remark.
- The output-size obstruction: with fixed slack, full listing exceeds every n^{a}.
- BCPZZ Thm 4.1's parameters (d = ⌈ε^{−3/θ}⌉ > k_RS) and Jeronimo's slack condition (UNVERIFIED by us).
- The GGG regime.
- CHHS's scope.
- The promise caveat: conditioning on non-empty windows divides by Pr(|𝓛| > 0).
- Its experiment proposals (beam search over divisor products; adaptive Boneh lattices) are recorded for future
  work. We ran the pointer/heavy-hitter form instead.

**Ours, beyond codex.**
- The code-specific capacity statement and the misplacement by the textbook formulas (F2(b, c)).
- The RS dictionary t = k²/2, n ≈ k^k (F3).
- The minimum-weight reading and the Cheng–Wan analogue (F4).
- The explicit depth bound of Prop DT and its experiment.
- Lemma NU, Lemma GC, Prop PT, Prop SG, Prop PNP, Lemma SC.
- Experiments E1–E2.

### 5.2 Adversarial round
Pending.

## 6. Opus critics
Pending.

## 7. Reproduction
```bash
.venv/bin/python -m pytest -q tests/test_crtdecode.py
.venv/bin/python -m src.crtdecode.run_exp E1     # ~15 s, 60 MB
.venv/bin/python -m src.crtdecode.run_exp E2     # ~45 s, 45 MB
.venv/bin/python -m src.crtdecode.analyze        # P1-P5
.venv/bin/python -m src.crtdecode.theory         # a*, c(a), dictionary
```
