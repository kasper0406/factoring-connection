# crtdecode: adaptive CRT decoding far below the Johnson radius (sinonobliv route 1)

Agent `crtdecode`, phase 23, 2026-10-22. Version 2: after the codex adversarial round (MAJOR REVISION) and both Opus critics (ACCEPT WITH FIXES; OK WITH FIXES), all adjudicated in §5.2 and §6.

**Files.**
- Code: `src/crtdecode/`.
  - `core.py`: exact window profiles, prime powers, smooth ranges, the exact idealised Johnson certificate.
  - `decoders.py`: the operational toy decoders (DT, DT8, PH-θ, RT) and the information decoder of Prop DT.
  - `theory.py`: c(a), a*, the regime dictionary.
  - `run_exp.py`, `analyze.py`: E1, E2 and their scoring.
- Tests: `tests/test_crtdecode.py` (11 pass).
- Pre-registration: `notes/crtdecode_prereg.md`, commit 7d07797, made before any code existed.
- Raw data: `experiments/raw/crtdecode_E{1,2}.jsonl`.
- GPT-6.1-sol: blind round `codex/crtdecode_1.{md,out}` (§5.1); adversarial round `codex/crtdecode_2.{md,out}` (§5.2).
- Post-hoc: `experiments/raw/crtdecode_E2_posthoc_crt_tree.json` (Prop DT′).

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
- **Where the problem sits.** Prop SI is list decoding of a Chinese-remainder code far *below* its Johnson radius,
  with parameters that approach the code's leading-order average-list threshold.
  - In the squarefree surrogate the code has n = π(y) prime coordinates and "dimension" k = ln h/ln y. Prop SI needs
    agreement in at least t = u₀ coordinates. At Prop SI's optimum, t ≈ k²/2 and n ≈ k^{k(1+o(1))}.
  - A Johnson-radius decoder needs agreement ≈ √(kn) = n^{1/2+o(1)} coordinates, not k²/2.
  - The average list size is (h + 1)ρ(u₀)^{1+o(1)} = L[2/3, ν − 2/(3κ)]. Its leading exponent changes sign at
    ν = 2/(3κ), which Prop SI's optimum approaches. At equality the lower-order terms decide (codex round 2), so we
    claim only this leading-order statement. A uniform-alphabet Hamming benchmark has a different threshold
    (ν = 1/κ), because the primes have very different sizes.
  - Brute force over messages (the sieve) costs h = n^{k(1+o(1))}. The target is n^{1.14}.
  - No CRT, number-field or Reed–Solomon decoder we know of works in this regime. That includes the 2026
    capacity-achieving RS decoders. BCPZZ's derivative order would exceed the message dimension and its running time
    the sieve's; Jeronimo's guarantee needs fixed additive slack, which tends to 0 here.
- **The decision-tree model of decoder calls cannot carry a depth barrier.** In this model the queries are
  idealised Johnson calls, chosen adaptively from earlier answers, each at unit cost, and the leaf computation is
  free.
  - O(x/ln h) + 1 = n^{o(1)} calls suffice (5 at our toy sizes). Querying q^{⌊log_q(h+1)⌋} for successive small primes reconstructs X by CRT, and one
    certificate follows (Prop DT′, codex round 2; post hoc: exactly 5 calls in all 200
    E2 windows).
  - The registered tree (Prop DT) uses ≤ 4π(y) + 1 + o(π(y)) calls on average; E2 confirms it.
  - Either way, no lower bound on the number of calls can reach a* in this model. The leaf computation is the
    sieve.
  - A comparison tree with precomputed constants does it in depth O(x) (Lemma NU).
  - Only the size of the call repertoire, or charged computation, can carry a lower bound.
- **What adaptivity provably cannot do, in restricted classes:**
  - Zooming into a divisor certified by an earlier call never finds anything that call did not already certify
    (Lemma GC).
  - A fixed pointer schedule with adaptively chosen tests yields at most ρ(u₀)^{1−o(1)} smooth outputs per pointer
    (Prop PT).
  - A fixed *unweighted* Johnson call on the whole window that certifies a smooth z only when z divides its
    modulus (no cofactor test) succeeds with probability at most (h + 1)X^{−3/4+o(1)}. On a region of length
    |R| ≥ 2 the bound is |R|·X^{−1/2+o(1)} (Prop SG). The weighted analogue is FALSE: codex's weights-on-2
    counterexample reaches ρ(u₀)^{1−o(1)}.
  - Fixed repertoires of any shape need ρ(u₀)^{−1+o(1)} *distinct calls in the repertoire* (portfolio size, not
    calls executed; Prop AT, sinonobliv).
- **A barrier for arbitrary algorithms is out of reach.** An unconditional lower bound of n^{a} for any fixed
  a > 0 would imply P ≠ NP (Prop PNP).
  - The remaining question is purely computational. Given the residue table, which costs n^{1+o(1)}, find a point
    of [0, h] covered by weight ≥ x of the progressions j ≡ −X (mod q^e). This is the problem HEAVY-POINT.
  - With an explicit prime list instead of all primes ≤ y it is NP-hard (DIVISOR).
- **Experiments (pre-registered, all five predictions held).**
  - E1: the natural adaptive aggregation rule tests first the positions hit by several smooth divisors. It needs
    4.4–5.5× fewer tests per success than the same pointers tested at first landing. But it costs 3.5–3.9× *more* per first smooth output
    than the oblivious divisor trick, because building collisions costs pointers.
  - E2: the information decoder (Prop DT) reconstructs every window exactly with 2.25–2.27·π(y) calls (post hoc,
    the shallow tree DT′ needs 5). Its leaf computation *is* the sieve, by construction: it costs ≈ 2.9·h incidences.
- **Assessment.** We put the probability that a CRT decoder with a < 1.1416 exists at about 1%; codex's blind
  estimate is 10% for an average-case version and 3% for all promised windows. The RS analogue, split polynomials in
  a translate of a low-degree space, is a deep-hole / Moments-Subset-Sum problem. No hardness or algorithm is known
  in its regime. Cheng–Wan's DLP link concerns the residue-class problem (Prop DLI), not the interval.
  - (Outside scope.) Quantum search over offsets with reversible trial division costs n^{1+o(1)}√h = L[2/3, ν/2] on
    promised windows (BBHT; codex). That halves the classical exponent but is still L[2/3]. Everything here is
    classical.

**Main results.**

| Result | Status |
|---|---|
| Prop F1: exact formulation (C_y, valuation agreement, threshold exact for h < X) | PROVED (= Prop RD/Boneh; codex) |
| Prop F2: regime: agreement fraction ≈ 1/n; Johnson ratio L[1/3, −κ/2]; average list L[2/3, ν − 2/(3κ)], whose leading exponent changes sign at ν = 2/(3κ) (the optimum approaches it); the uniform-alphabet Hamming benchmark has threshold ν = 1/κ | PROVED (averaged; leading order); average list via Hildebrand–Tenenbaum |
| Prop F3: squarefree surrogate dictionary (n, k, t): t/k² → κ/ν², ln n/(k ln k) → 3κ²/ν; at (κ*, ν*): t ≈ k²/2, n = k^{k(1+o(1))}; sieve = n^{k(1+o(1))} | PROVED (algebra, leading order) |
| Prop F4: smooth ⇔ pointwise equality wt(c(m)) = ln P − ln m; Prop SI = find an equality case in the translate X + [0, h] (the minimum-weight word is the largest one); RS analogue = deep hole / Moments Subset Sum over D = F_q at d ≈ k²/2 | PROVED (F4); analogue: deep-hole/MSS literature |
| Prop DT / DT′: decision trees of idealised decoder calls, free leaf: depth ≤ 2π(y) + 1 + Σ_q m_q (E ≤ 4π(y) + 1 + o(π(y))), and O(x/ln h) + 1 = n^{o(1)} for the CRT tree DT′; leaf = sieve | PROVED in the unit-cost oracle model; E2 (DT), post hoc (DT′) |
| Lemma NU: a comparison tree with precomputed constants, depth ≤ log₂ Ψ + 1 = O(x), solves Prop SI nonuniformly | PROVED (trivial) |
| Lemma GC: zoom chains stay inside the first call's certified set | PROVED (trivial) |
| Prop PT: oblivious pointer schedules (d ≥ h/2), adaptive tests: ≤ T·ρ(u₀)^{1−o(1)} smooth outputs | PROVED (averaged; via Prop JC's estimate) |
| Prop SG: fixed unweighted full-divisor certificates (z \| P) on the window: ≤ (h + 1)·X^{−3/4+o(1)} smooth outputs per call; weighted analogue false | PROVED (averaged) |
| Prop PNP: an unconditional n^{a} lower bound (a > 0 fixed) implies P ≠ NP | PROVED (standard) |
| Lemma SC: exact interval counting of smooth numbers in time C gives Prop SI in C·log₂(h + 2) + poly | PROVED (trivial) |
| No decoder with a < 1.1416 among the mechanisms in §3.2 | finding, not a theorem |

---

## 1. Formalisation and regime map (Question 1)

### 1.1 The code and the problem

**Definition.** For X ∈ [X₀, 2X₀) (with X + h ≤ 3X₀, as holds at Prop SI's scale) let e_q = ⌊ln(3X₀)/ln q⌋ and P = ∏_{q≤y} q^{e_q}. The code C_y has
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
    R_sf ≪ α_sf ≪ √R_sf, with α_sf := x/W_sf. Only squarefree smooth z reach full agreement there.
- (b) *Average list size (codex §1; Hildebrand–Tenenbaum).*
  E|𝓛_X| = ((h + 1)/X₀)·(Ψ(2X₀, y) − Ψ(X₀, y)) + O(h²/X₀) = (h + 1)ρ(u₀)^{1+o(1)} = L[2/3, ν − 2/(3κ)].
  - For ν < 2/(3κ), most windows have an empty list.
  - For ν > 2/(3κ), the mean list size L[2/3, ν − 2/(3κ)] exceeds every power of n. So listing all solutions at
    cost n^{a} per window is impossible on average (codex's output-size obstruction).
  - **The leading exponent of the mean list size changes sign at ν = 2/(3κ), and Prop SI's optimised parameters
    approach this boundary** (ν ↓ 2/(3κ)). This is a leading-order average-list threshold for the special received
    words −X, X in a dyadic interval. It is neither a channel-capacity theorem nor a worst-case list-decoding theorem
    (codex round 2, finding 1).
    - At equality, lower-order terms decide. For example, h = ρ(u₀)^{−1}e^{±g} with ln n ≪ g ≪ F_{2/3} has the same
      ν, but the mean list vanishes for one sign and exceeds every power of n for the other.
    - To claim a mean list of about 1, h must be set from the actual friable density.
    - By averaging, some received word has a list at least the mean. Whether worst-case lists stay small near the
      threshold is OPEN.
    - RS precedent (literature critic): the Justesen–Høholdt counting bound (IEEE IT 47(4), 2001) and Cheng–Wan's
      threshold ĝ(n, k, q), which they show lies far below √(nk). For CRT codes we found no prior statement of
      F2(b, c); "no prior statement found by us".
- (c) *Benchmarks with other metrics have other thresholds* (codex round 2, finding 2: these are different metrics,
  not failures of standard capacity theory on its own model).
  - The large-alphabet benchmark α = R (in the compressed normalisation of (a)) is met with a polylog margin.
  - The uniform-alphabet Hamming surrogate has rate R_H = ln V/(n ln y) and agreement α_H = t/n, with q ≈ y on n
    coordinates. Its entropy formula R_H = 1 − H_q(1 − α_H) gives expected list
    V·C(n, t)·y^{−t} = V·exp(−u₀(ln u₀ + ln ln y − 1 + o(1))). Its threshold is ν = 1/κ = (2/(3κ) + 1/(3κ)), not
    2/(3κ).
  - The difference u₀(ln ln y − ln ln u₀) = (1/(3κ) + o(1))F_{2/3} comes from the small primes. Their coordinates
    have small alphabets, so agreement there is cheap, and the smooth-number count Ψ, not a q-ary entropy, governs
    the list.
  - Numerically (`theory.dictionary`, o(1) dropped), R_H/R_cap = 0.80, 0.79, 0.77, 0.76 at 512, 1024, 2048
    and 4096 bits, drifting towards νκ = 2/3.
- Status of (c): exact algebra on a stated q-ary model (PROVED in MODEL); the comparison is the point.
- *Proof of (c).* Stirling and de Bruijn's expansion ln ρ(u) = −u(ln u + ln ln u − 1 + o(1)), with
  ln ln y = (1/3 + o(1)) ln x and ln ln u₀ = ln ln x + O(1). ∎

### 1.3 The Reed–Solomon dictionary

**Prop F3 (PROVED; algebra, leading order; a surrogate).** For the squarefree prime code (n = π(y) coordinates of
alphabet ≈ y), put k := ln V/ln y (the number of "symbols" in a message) and t := x/ln y = u₀. t is a *lower bound* on
the number of agreeing primes a squarefree smooth X + j needs. A large power of 2 has a single agreeing prime, so the
valuation code of §1.1 is the actual problem and this dictionary is a surrogate (codex round 2). Then:
- t/k² = κ/ν²·(1 + o(1)) and ln n/(k ln k) = (3κ²/ν)(1 + o(1));
- at (κ*, ν*): t/k² → 9κ*³/4 = 1/2 and ln n/(k ln k) → 9κ*³/2 = 1;
- brute force over messages costs V = y^k = n^{k(1+o(1))}; the sieve is this brute force;
- the Johnson threshold needs √(kn) = n^{1/2+o(1)} agreeing primes.

*Proof.* k = (ν/κ)(x/ln x)^{1/3}, t = (1/κ)(x/ln x)^{2/3}, ln n = ln y − ln ln y = κF_{1/3}(1 + o(1)), and
ln k = (1/3 + o(1)) ln x. ∎

Convergence is slow: ln n/(k ln k) = 0.81, 0.86, 0.90, 0.93 at 512–4096 bits (`theory`), i.e. n ≈ k^{0.86k} at
1024 bits. Numerically (o(1) dropped): at 1024 bits n ≈ 8.9·10⁶, k ≈ 8.7, t ≈ 37.5, √(kn) ≈ 8800. At 2048 bits
n ≈ 4.7·10⁹, k ≈ 10.5, t ≈ 55.6, √(kn) ≈ 2.2·10⁵.

**Reading.** Prop SI asks to decode, in time n^{1.14}, a code with dimension k ≈ ln n/ln ln n, from agreement in
t ≈ k²/2 of n ≈ k^k coordinates (polylogarithmic agreement). Brute force over messages, or over k-subsets of
coordinates, costs n^{Θ(k)}: quasi-polynomial in n, and that is the sieve. In RS language this is polynomial
reconstruction far below √(kn), with agreement far above the large-alphabet threshold (t ≫ k), so it is
information-theoretically decodable.
- Related RS hardness assumptions: Kiayias–Yung, "Cryptographic hardness based on the decoding of Reed–Solomon codes"
  (ICALP 2002, LNCS 2380, 232–243), a decisional polynomial-reconstruction assumption for randomly corrupted
  instances.
- Noisy polynomial *interpolation* (several candidate values per point; Naor–Pinkas, STOC 1999) is a different and
  easier problem. Bleichenbacher–Nguyen (EUROCRYPT 2000, LNCS 1807, 53–69) reduce it to SVP and discuss noisy CRT by
  analogy.
- None of these is a hardness transfer to the coherent received words −X here (codex, literature critic).

### 1.4 Minimum-weight codewords and the function-field analogue

**Prop F4 (PROVED).**
- For 1 ≤ m ≤ 3X₀, the weighted weight of c(m) is ln P − ln gcd(m, P) ≥ ln P − ln m. Equality holds iff m | P, iff m
  is y-smooth (v_q(m) ≤ e_q holds automatically for m ≤ 3X₀).
- Pointwise: wt(c(m)) − (ln P − ln m) = ln(m/gcd(m, P)), which vanishes iff m is smooth. So **Prop SI = find an
  equality case of this Singleton-type bound in the translate {c(X + j) : 0 ≤ j ≤ h}** (an interval translate, not a
  group coset).
- Smoothness does not mean minimum raw weight. Codex's example: X = 8, h = 2, y = 5, P = 720; 8, 9 and 10 are all
  smooth, but only 10 minimises the weight. On a non-empty promised window, minimising the raw weight returns the
  *largest* smooth integer, since ln(y + 1) > ln((X + h)/X). Finding any smooth integer is a weaker task.
- RS analogue. With evaluation set A ⊆ F_q, the evaluation vector of f has weight ≥ |A| − deg f, with equality iff f
  splits into distinct linear factors over A. The analogue is: find f in the translate g + F_q[T]_{<k} (deg g = t)
  that splits over A. Equivalently, bounded-distance decode the received word −g|_A to RS[|A|, k] at agreement t.
- The dictionary of §1.3 is the classical smooth-integer ↔ split-polynomial analogy, which we use only as an
  analogy.

**Literature on the analogue (checked by the literature critic in the primary texts unless noted).**
- **The short-interval analogue is the deep-hole / Moments-Subset-Sum setting.** A received word that is the
  evaluation of a degree-(k + d) polynomial g has a codeword at agreement k + d iff g − c splits over the evaluation
  set D. Equivalently, a split polynomial has its top d coefficients prescribed. References:
  - Cheng–Murray, "On deciding deep holes of RS codes" (TAMC 2007, LNCS 4484);
  - Li–Wan, "On the subset sum problem over finite fields" (FFA 2008; second-hand);
  - Gandikota–Ghazi–Grigorescu, Moments Subset Sum (arXiv 1611.03069; SIAM J. Comput.). NP-hard for
    d ≤ c·log N/log log N at K = Θ(N), and up to c·log N under quasi-polynomial reductions; constructed evaluation
    sets, field size exp(n);
  - Wan–Zhang (arXiv 2609.29120, Sept 2026; abstract only). Deterministic NP-completeness for d = ⌊n^α⌋, α < 1/2,
    at rate ≈ 1/2, constructed evaluation sets.
- Prop SI's analogue has d = t − k ≈ k²/2, D = all of F_q, and q ≈ k^k. The structured D plays the role of "all
  primes ≤ y", mirroring the DIVISOR remark of §2.3. **No hardness result or algorithm known to us covers this
  regime.**
- **Cheng–Wan is the analogue of the residue-class problem, not of the short interval** (literature critic, MAJOR;
  codex round 2, finding 10).
  - Cheng–Wan: "On the list and bounded distance decodability of Reed–Solomon codes", FOCS 2004; SIAM J. Comput.
    37(1) (2007) 195–209.
  - They seek t·h + f, splitting over S, with h irreducible. That is a smooth element in a residue class modulo an
    irreducible, the function-field form of index-calculus relation finding, i.e. of sihard's Prop DLI (a
    modulus-generic Prop SI constructor). Their received word is a rational function on S.
  - Their Thm 1: list decoding at the counting threshold ĝ ⇒ DLP. Their Thm 2: BDD of RS[q, 3h + 4] at agreement
    4h + 4 ≈ 4k/3 ⇒ DLP in F_{q^h}. Both reductions are randomized.
  - Under T ↦ 1/T, the note's translate g + F_q[T]_{<k} becomes a residue class mod T^{t−k+1}. The DLP in
    (F_q[T]/T^m)^* is easy, so Cheng–Wan's DLP hardness does **not** transfer to Prop SI's analogue. The v1 claim that
    Cheng–Wan "study exactly this regime" is withdrawn.
  - Cheng–Wan's ĝ (with the Justesen–Høholdt counting bound, IEEE IT 2001) is the RS precedent of F2(b)'s
    average-list threshold lying far below √(nk). Their later "Complexity of decoding positive-rate RS codes" (arXiv
    0802.1220) extends the reduction to all positive rates for non-prime q.
- The DLP link that does exist is sihard's Prop DLI: a modulus-generic Prop SI constructor would give prime-field
  discrete logarithms in L[1/3, 1.232] (HEURISTIC, inherited).
- *Small characteristic.* Frobenius identities (homographic images of x^q − x = ∏(x − a)) produce split elements in
  structured quotients. They are the engine of the heuristic quasi-polynomial DLP of Barbulescu–Gaudry–Joux–Thomé
  (EUROCRYPT 2014), after Joux (SAC 2013) and Göloğlu–Granger–McGuire–Zumbrägel (CRYPTO 2013); a rigorous version is
  Kleinjung–Wesolowski (J. AMS 2022). We know of no integer analogue of the identity (assessment).

### 1.5 Literature comparison (bibliography and descriptions checked by the literature critic and codex)

| Result | What it decodes | Here |
|---|---|---|
| Goldreich–Ron–Sudan, IEEE IT 46(4) (2000) 1330–1338 | CRT codes, agreement ≳ √(2k_c n log p_n/log p_1) | Johnson-type, √2 short; L[1/3, κ/2] too high |
| Boneh, STOC 2000, 265–272; JCSS 64 (2002) 768–784 | removes the √2; smooth integers in intervals (gcd form) | Johnson; vacuous at full smoothness (sihard §4) |
| Guruswami–Sahai–Sudan, FOCS 2000 | weighted ("soft") CRT, relatively prime moduli, equality scores; unweighted t ≥ √(k(n + ε)) | Johnson. Its threshold is √(W·ln V) in nats. The nested valuation score is a specialisation in Theorem J's idealised class, not GSS's published theorem |
| Cohn–Heninger, "Ideal forms of Coppersmith's theorem and Guruswami–Sudan list decoding" (ICS 2011) | unifies Coppersmith, Howgrave-Graham, Boneh-CRT and GS | the natural reference for "Johnson-class certificates" |
| Chinburg–Hemenway–Heninger–Scherr, ASIACRYPT 2016 | no Coppersmith-type auxiliary polynomial finds roots of size N^{1/d+ε} **mod N** | modulus-N case only. The unknown-divisor (Boneh/CRT) case is listed as future work there, so CHHS does not cover this family. Theorem J is the only Johnson-limit statement for it that we know |
| Beyond-Johnson CRT variants: Cohn–Heninger multi-sample ACD (ANTS 2012); interleaved CRT (Li–Sidorenko–Nielsen, ISIT 2013) | need several correlated samples, or a different code | one sample here (sinonobliv route 3) |
| CRT capacity | none known to us. BHKS's ideal-theoretic framework (RANDOM 2021 / IEEE IT 2024) covers polynomial ideals only; Guruswami 2003 constructs number-field codes, without a capacity result | — |
| Brakensiek–Chen–Putterman–Zhang–Zheng, arXiv 2609.08005, Thm 4.1 (statement checked; proof not audited) | plain RS to capacity over prime q: d = ⌈ε^{−3/θ}⌉, needs k > d, rate ≤ (1 − θ)ε, q ≥ max(n, 4ε^{1−9/θ}n/k), time q^{O(ε^{−12/θ})} | ε ≈ L[1/3, −κ]: d ≫ k fails, the field-size condition fails, and the time q^{n^{Ω(1)}} is worse than the sieve; no CRT transfer either |
| Jeronimo, arXiv 2609.05870 (statement checked) | radius 1 − R − γ, prime q ≥ n, q^{C(γ,L)} time; uniform in k but not as γ → 0 | the additive slack γ = (t − k)/n → 0 |
| Ben-Sasson–Kopparty–Radhakrishnan, IEEE IT 56(1) (2010) 113–120; Guruswami–Rudra, IEEE IT 52(8) (2006) | large RS lists beyond Johnson (subspace polynomials); limits to list recovery | worst-case lists, not finding one element |
| Guruswami–Vardy, IEEE IT 51(7) (2005) 2249–2256; GGG; Wan–Zhang 2026 | NP-hardness of ML / BDD RS decoding, constructed instances | different regimes (§1.4); no CRT transfer |
| Cheng–Wan, FOCS 2004 / SICOMP 2007 | RS list decoding at ĝ, and BDD at agreement ≈ 4k/3, ⇒ DLP | analogue of the residue-class (Prop DLI) problem, not of Prop SI's interval (§1.4) |

*Where the 2026 RS decoders fail (answering the brief).*
- BCPZZ: derivative order ≫ dimension, the field-size condition fails, and the time is worse than the sieve.
- Jeronimo: the additive slack tends to 0.
- Both are polynomial-ring decoders. C_y has no global derivation, and its nested prime-power coordinates are not
  derivatives. So even favourable parameters would leave an unproved transfer theorem.

## 2. Adaptive decoding (Question 2)

### 2.1 Decision trees of decoder calls: information-trivial in depth

The model, as the brief suggests: each query is one decoder call, chosen as any function of earlier answers, at unit
cost. X enters only through answers. The calls are idealised unweighted Johnson certificates (Theorem J's class
𝒥): modulus M (y-smooth), region R = [X + ℓ, X + ℓ + |R| − 1], certified set
{z ∈ R : (ln gcd(z, M))² ≥ ln|R|·ln M}. The answer lists each certified z with gcd(z, M).

**Prop DT (PROVED in the unit-cost oracle model).** Let X ∈ [2^{b−1}, 2^b) with b public, h < 2^{b−1} and y ≤ h + 1 (so a_q ≥ 1). The following tree outputs a y-smooth
z ∈ [X, X + h] with its factorisation iff one exists.
1. For each q ≤ y, with a_q = ⌊log_q V⌋ ≥ 1, call (q^{a_q}, [X, X + q^{a_q} − 1]). The region has length equal to
   the modulus, so exactly the unique multiple is certified. This reveals X mod q^{a_q}.
2. For each q and e = a_q + 1, a_q + 2, …: call (q^e, window). Stop at the first e at which no point has gcd q^e.
   Since q^e > V, every point with v_q ≥ e is certified.
3. Leaf: all valuations in the window are now known. X + j is smooth iff its smooth part is ≥ 2^{b−1}: a non-trivial
   cofactor is ≥ y + 1 ≥ 3 and X + j < 3·2^{b−1}. Certify the least smooth z with one call (M = z, window).

Depth D(X) ≤ 2π(y) + 1 + Σ_q m_q, with equality unless a level cap ends a prime's step 2 early or the window has no
smooth integer (then there is no final call); here m_q is the number of levels e > a_q present. Always D ≤ π(y)(2 + log₂(3X)) + 1
= π(y)^{1+o(1)}, and E_X D ≤ 4π(y) + 1 + O(π(y)·b·V/2^{b−1}).

*Proof.* Correctness is by construction. For the expectation:
Pr[some j ∈ [0, h] has q^e | X + j] ≤ (V/q^e)(1 + h/2^{b−1}) + V/2^{b−1}, and Σ_{e>a_q} V/q^e < q/(q − 1) ≤ 2
because q^{a_q+1} > V. The extra terms add ≤ 2h/2^{b−1} + e_q·V/2^{b−1} per prime. Every step-2 answer has ≤ 1 certified
point with b² ≥ e·a_q (proof critic), so unit cost per call is honest. ∎

- **So in depth the threshold is beaten informationally: a ≤ 1.** What is proved is that no *depth* (call-count)
  lower bound above a = 1 exists in this model, not that no barrier exists. With regions not restricted to the
  window, one call (M = P, |R| = P) reveals X mod P, and the depth is 2.

**Prop DT′ (the shallow tree; codex round 2, finding 20; PROVED in the same model; POST-HOC check).** Query
(q^{a_q}, [X, X + q^{a_q} − 1]) for q = 2, 3, 5, … until ∏ q^{a_q} > 2^b + h. The answers give X mod a modulus
exceeding X, so X itself by CRT; the free leaf sieves and one call certifies.
- Since q^{a_q} > V/q, the number of calls is ≤ b ln 2/ln(V/y) + 2 = O(x/ln h) + 1 = n^{o(1)}.
- POST-HOC (`decoders.crt_query_tree`, `experiments/raw/crtdecode_E2_posthoc_crt_tree.json`): on all 200 E2 windows,
  exactly 5 calls (4 residue queries + 1 certificate), with exact reconstruction and the same output as the sieve.
- This strengthens Prop DT, as Lemma IC already suggested: the decision-tree model with free leaves does not just
  fail to carry a barrier at a*, it makes the problem cost n^{o(1)} calls. Everything is in the leaf.
- The leaf computation is the sieve on known residues. Its incidence count equals the sieve's by construction
  (checked in E2); it is ≈ 2.9·h.
- **Lemma NU (PROVED).** For the family X ∈ [X₀, 2X₀), binary search among the S smooth integers of
  [X₀, 2X₀ + h] is a comparison tree ("X ≥ c?") of depth ⌈log₂(S + 1)⌉ + 1 ≤ x/ln 2 + 3. It outputs a smooth
  element of the window whenever one exists. So depth in comparison trees with arbitrary precomputed constants is
  O(x). This is not a statement about every nonuniform measure: the permitted queries matter (codex).
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
call. *Proof.* gcd(z, P₁) ≥ d ≥ T. ∎ The same holds for weighted scores, which depend on z only through
gcd(z, P₁) and are monotone in it (proof critic).
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
- *Further scope (codex).* The bound needs a deterministic pointer cap and an X-independent schedule. It preserves
  only the leading L[2/3] exponent, so it does not exclude L[1/3]-factor improvements. Translating it to promised
  (non-empty) windows divides by Pr(|𝓛_X| > 0).

**Prop SG (fixed, unweighted, full-divisor certificates on the window; PROVED, averaged over X).** Suppose a
fixed unweighted call (P, window) certifies a smooth z only when z | P, with no cofactor test. (We call this a
"full-divisor" certificate; it is an imposed restriction, not standard decoding terminology. In decoding,
"full agreement" would mean P | z.) Then
Johnson's condition (ln z)² ≥ ln V·ln P forces ln P ≤ Λ_J := (ln(3X₀))²/ln V = x^{4/3−o(1)}. Every such z is supported
on S = supp(P), and θ(S) := Σ_{q∈S} ln q ≤ Λ_J. Hence
  E_X #{certifiable smooth z in the window} ≤ (V/X₀)·Ψ(3X₀, z_J) ≤ V·X^{−3/4+o(1)},
where z_J is the largest prime with θ(z_J) ≤ Λ_J, so z_J = (1 + o(1))Λ_J = (ln X)^{4/3+o(1)}.
- *Proof.* Map the i-th prime of S to the i-th prime overall. This injection on S-smooth integers does not increase
  values (Theorem J step 4(b)), so #{m ≤ 3X₀ : supp m ⊆ S} ≤ Ψ(3X₀, p_{|S|}). We have p_{|S|} ≤ z_J because the |S|
  smallest primes have θ ≤ θ(S) ≤ Λ_J. Then use Ψ(X, (ln X)^A) = X^{1−1/A+o(1)} for fixed A > 1 (KNOWN; de Bruijn
  1966; see Granville's survey arXiv 1212.1579),
  with A = ln z_J/ln x → 4/3. Since A varies with X, sandwich it between the fixed exponents 4/3 ± η and use
  monotonicity of Ψ in y (codex). ∎
- *Reading.* An unweighted Johnson-class call reaches z | P only with a modulus so small (ln P ≤ x^{4/3}, against
  ln P ≈ n·x for the full code) that the integers it can certify that way are exponentially sparse. For a
  repertoire of K moduli multiply by K. (The v1 deduction that *every* Johnson-class success must use a cofactor test
  is withdrawn: it fails for weighted calls, see below.) At practical sizes A ≈ 1.22, so the exponent is ≈ −0.82, not −3/4 (`theory`, o(1) dropped).
- *Other regions (proof critic).* For a region of length 2 ≤ |R| ≤ V the same proof gives Λ = (ln 3X₀)²/ln|R|,
  A → 2 and ≤ |R|·X^{−1/2+o(1)}: still exponentially small. For |R| = 1 the threshold is 0 and the "call" is just a
  smoothness test of one point (yield ρ(u₀)). The bound is loose: a modulus carrying full prime powers supports only
  ≈ x^{1/3} primes.
- *Weighted calls: the analogue is FALSE (codex round 2, finding 5).* Take d = 2^a ≥ V and
  P = d·∏_{3≤q≤y} q^{E_q}, with weight 1 on 2 and 0 elsewhere. Every multiple of d in the window passes the weighted
  Johnson condition. Smooth multiples with odd cofactor divide P, and their averaged yield is ρ(u₀)^{1−o(1)} ≫
  V·X^{−3/4}. (This is the divisor trick in weighted dress; it is still within Theorem J.)

**Prop AT (sinonobliv; codex's portfolio form).** Any adaptive use of a fixed portfolio of K certificate calls has
success probability ≤ K·ρ(u₀)^{1−o(1)}. The bound is on the portfolio size, not the running time.

**The ladder.**

| Adaptivity allowed | Bound | Status |
|---|---|---|
| Fixed schedule (Theorem J, Prop JC) | ρ^{1−o(1)} per call | PROVED |
| Fixed schedule of pointers, adaptive tests (Prop PT) | ρ^{1−o(1)} per pointer | PROVED |
| Zoom into certified divisors (Lemma GC) | inside the first call's set | PROVED |
| Adaptive selection from a fixed repertoire R (Prop AT) | ≤ \|R\|·ρ^{1−o(1)} | PROVED |
| Fixed unweighted full-divisor certificates (Prop SG) | V·X^{−3/4+o(1)} per call on the window; \|R\|·X^{−1/2+o(1)} for 2 ≤ \|R\| ≤ V | PROVED (weighted analogue false) |
| Moduli computed from answers (decision tree, free leaf; Props DT, DT′) | depth n^{o(1)} (DT′) | PROVED: no depth barrier at all |
| Sub-threshold arithmetic search (residue-guided beam search over divisor products, reusing failed candidates; codex) | none known | OPEN: not covered by GC/PT/AT; no gain shown |
| Uniform classical computation | OPEN; an n^{a} lower bound ⇒ P ≠ NP (Prop PNP) | — |

### 2.3 The residual computational problem

**HEAVY-POINT_y.** Input: X, h and the residue table r_{q^e} = −X mod q^e (q ≤ y, e ≤ e_q), of size n^{1+o(1)} words.
Output: j ∈ [0, h] with Σ_{q^e : j ≡ r_{q^e}} ln q ≥ ln(X + j).
- The table costs n^{1+o(1)} (Lemma IC). Given the table, Prop SI and HEAVY-POINT_y are the same problem.
- The sieve solves it in O(h ln ln y).
- With an explicit list of primes in place of all primes ≤ y it is NP-hard under randomized reductions (DIVISOR;
  sihard Prop XS; Shor's and Barak's cstheory.SE answers, by randomized reductions or deterministically under a
  prime-gap hypothesis. The Alon–Kilian credit is sihard's and was not confirmed by the literature critic).
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
- (Assessment, narrowed after codex.) Relative approximations to the two *cumulative* counts Ψ(b, y) ≈ bρ
  (saddle point; Hildebrand–Tenenbaum) cannot resolve their difference, which is of order 1. A relative
  approximation directly to the *interval* count would suffice, but we know of no such oracle. Dyer's FPRAS
  counts one-sided knapsack solutions and does not supply this thin-slab count.
- We know of no exact method below L[2/3] at this scale (assessment).

## 3. Lower bound or construction (Question 3)

### 3.1 Answer

- **No construction with a < 1.1416.** The divisor trick remains best in every adaptive class we analysed.
- **No lower bound in an unrestricted model.** That is impossible without proving P ≠ NP (Prop PNP).
- **The decision-tree model of decoder calls has the opposite answer: n^{o(1)} calls** (Props DT, DT′). Lemma
  IC's cheap information is exploitable by an adaptive decoder, and everything left is the leaf computation.
- **Barriers in restricted adaptive classes** are Lemma GC, Prop PT, Prop SG and Prop AT (§2.2).

### 3.2 Mechanisms examined and costed

| Mechanism | Cost per output | Why |
|---|---|---|
| Pointer heavy-hitters (smooth positions collect many divisor hits) | ≥ ρ^{−1+o(1)} (Prop PT); E1: 3.5–3.9× worse than DT | Heuristic count: enhancement ≈ (ω ln y)^{k} = L[1/3, Θ(1)] (the divisor-count gain), but one hit on a smooth position costs ≈ V/E pointers = L[2/3] |
| Iterated Boneh / zoom-in | ≤ first call (Lemma GC) | certified progressions are already certified |
| Support guessing + Johnson (random S of Johnson-admissible size) | exponential | Prop SG (unweighted): any fixed S supports X^{1/4+o(1)} integers ≤ 3X; biased S does not help, since a smooth integer's prime mass sits near y |
| High-density subset-sum methods (BKW/Lyubashevsky, Flaxman–Przydatek, Wagner k-tree) | no route found (assessment) | combinations ≈ (X/h)·(hρ): density L[2/3, ν − 2/(3κ)] only through the window slack. A solution has ≥ u₀ prime factors but ≤ log₂(2X) ≪ n nonzero exponents (sparse relative to n, not ≤ u₀ — v1's "≤ u₀" was reversed, codex). Fixed product tables, MITM included, need lists of pieces too large for their domain (sihard Prop MM; Wagner budget HEURISTIC). These methods are not proved inapplicable. |
| Exact counting + binary search (Lemma SC) | as costly as exact slab counting | no exact method known below L[2/3] |
| Hidden-derivative / multiplicity decoding | not applicable | derivative order ≫ dimension; prime-power coordinates are not derivatives (§1.5) |
| Lattices (centred CVP, Prop CV) | exact promise-CVP in dimension n in time n^{1.14} | would need near-linear exact CVP for one structured family; codex's "most promising" mechanism, but no algorithm is known |
| Generic additive models (comparisons of Σe_i ln p_i only) | no transferable bound | Erickson's Ω(n^{⌈r/2⌉}) for r-linear decision trees concerns arbitrary real inputs. It gives no bound for the fixed weights ln p with a moving interval target, so the v1 transfer is withdrawn (codex). Three gaps: r = u₀ grows, and a u₀^{−Θ(u₀)} loss in the constants is L[2/3]-sized; the inputs are fixed; the target is an interval. General linear decision trees solve knapsack in polynomial depth (Meyer auf der Heide, JACM 1984), and Kane–Lovett–Moran give O(rn log² n)-depth comparison trees for r-SUM, consistent with Lemma NU. |

### 3.3 Conditional evidence
- A CRT decoder with a < 1.1416 that works for arbitrary X is a modulus-generic Prop SI constructor (Lemma EQ). By
  sihard's Prop DLI it would also give prime-field discrete logarithms in L[1/3, max(c(a), 1.232)], with
  c(a) < 1.902. This inherits sihard's heuristics (H_S, H-affine and the booting step). It would therefore beat two
  long-standing records.
- Cheng–Wan's RS ⇒ DLP reductions are the function-field analogue of Prop DLI (the residue-class problem), not of
  Prop SI's interval. They are not evidence about Prop SI itself (§1.4).

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
- The groups are decoder-pure (proof critic). "1 hit" is exactly DT8's tests (300 successes); "≥ 2 hits" is exactly
  PH's (300 + 300, resp. 300 + 298); PH never reached its post-exhaustion phase. Each decoder stops at its first
  success, so each rate is (#successes)/(#tests). P2 therefore measures "PH needs 4–5× fewer tests per success than
  DT8". It is not a within-population enrichment, and the ratio estimator is biased under sequential stopping. The
  registered prediction held as worded.
- **Reading.** Aggregation works as designed: tests on doubly hit positions succeed 4–5× more often, consistent
  with the divisor-count gain that makes smooth numbers heavy hitters. But PH spends ≈ 98% of its budget on pointers to create
  those collisions, and loses 3.5–3.9× to the oblivious divisor trick.
- Toy sizes cannot show asymptotics. Prop PT says that every oblivious pointer schedule, PH's included, needs
  L[2/3, 2/(3κ) − o(1)] pointers per output. Its ρ^{1−o(1)} hides e^{O(F_{1/3})} factors, so it does not rank PH
  against DT at L[1/3] scale; that is irrelevant against an n^{a} = L[1/3] target.
- (POST-HOC) DT8 loses 1.2× to DT, plausibly because moduli above h often have no multiple in the window.
- P3's second inequality (RT < SV) is nearly automatic: SV counts incidences, ≈ 2.9× the first smooth offset.

### 4.2 E2 (100 windows per config)

| cfg | π(y) | calls (mean / max) | step-1 calls | exact | leaf incidences / h | leaf = sieve |
|---|---|---|---|---|---|---|
| C1 | 97 | 220.1 / 231 | 97 | 100/100 | 2.87 | yes (by construction; checked) |
| C2 | 172 | 386.8 / 400 | 172 | 100/100 | 2.98 | yes (by construction; checked) |

- Depth ≈ 2.25–2.27·π(y): step 1 contributes π(y) calls (one per prime: no two of the prime powers q^{a_q} fit
  under V together), step 2 contributes π(y) + Σ m_q, and one final certification call.
- The tree certified the least smooth integer in every window.
- The leaf's incidence count equals the full window sieve's count exactly. Prop DT's informational a = 1 is bought by
  a leaf that *is* the sieve.

### 4.3 Deviations
- Costs are *registered abstract costs* (pointers, tests, incidences). E1 excludes the lcm bookkeeping and sorting
  by design.
- E2's code caps step-2 exponents with the public bound 2^b + h rather than the pre-registered e_q = ⌊log_q(X + h)⌋
  (which would read X). This is benign and enforces oracle-only access (codex round 2).
- After E2 (codex finding 19), `JohnsonOracle.call` now raises on unsupported general-modulus calls. The old code
  silently omitted partial-gcd points there (e.g. X = 15, M = 30, |R| = 2). The exercised calls (|R| = M, prime
  powers, and M = X + j* in its own window) are unchanged, and E1/E2 results are unaffected (tests re-run).
- POST-HOC: Prop DT′'s shallow tree (`crt_query_tree`) on the 200 E2 windows.
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
- BCPZZ Thm 4.1's parameters (d = ⌈ε^{−3/θ}⌉ > k_RS) and Jeronimo's slack condition. Statements checked by codex and
  the literature critic; proofs not audited.
- The GGG regime.
- CHHS's scope.
- The promise caveat: conditioning on non-empty windows divides by Pr(|𝓛| > 0).
- Its experiment proposals (beam search over divisor products; adaptive Boneh lattices) are recorded for future
  work. We ran the pointer/heavy-hitter form instead.

**Ours, beyond codex's blind round.** Not novelty claims. F1, NU, GC, PNP and SC are standard or elementary.
- The leading-order average-list threshold and the uniform-Hamming comparison (F2(b, c)).
- The surrogate RS dictionary t ≈ k²/2, n ≈ k^{k(1+o(1))} (F3).
- The Singleton-equality reading (F4).
- The explicit tree of Prop DT and its experiment.
- Props PT and SG.
- Experiments E1–E2.
- Prop DT′ is codex's, from round 2.

### 5.2 Adversarial round (`codex/crtdecode_2.{md,out}`, on v1 5cb24e0): MAJOR REVISION, all points adjudicated
Codex reproduced E1, E2, the scoring and the tests in scratch: 800 records match after removing timing fields, and
P1–P5 and both tables check. It found no leakage or accounting bug.

| # | Sev. | Finding | Adjudication |
|---|---|---|---|
| 1 | MAJOR | "exactly at capacity" not established; lower-order terms decide at ν = 2/(3κ) | accepted: "leading-order average-list threshold", which the optimum approaches; e^{±g} example quoted; "planted" removed |
| 2 | MINOR | uniform-alphabet benchmark is a different metric, not a failure of capacity theory | accepted: surrogate symbols R_H, α_H; reworded |
| 3 | MINOR | t = x/ln y is a lower bound; the squarefree dictionary is a surrogate; X + h ≤ 3X₀ | accepted |
| 4 | MAJOR | smooth ≠ minimum raw weight (X = 8, h = 2, y = 5); the interval is not a coset | accepted: pointwise Singleton-equality form; example quoted (the proof critic flagged the same as a NIT) |
| 5 | MAJOR | Prop SG false for weighted calls (weights on 2); "every Johnson success uses a cofactor test" unjustified | accepted: SG restricted to fixed unweighted full-divisor certificates; weighted analogue stated FALSE; deduction withdrawn |
| 6 | MINOR | A(X) → 4/3 needs a fixed-exponent sandwich | added |
| 7 | MINOR | D = 2π + 1 + Σm_q not always exact; label | "≤", equality conditions; "PROVED in the unit-cost oracle model" |
| 8 | MINOR | NU only for comparison trees with constants; distinguish set containment / query / runtime | accepted |
| 9 | MINOR | PT preserves only the L[2/3] exponent; deterministic cap; promise conditioning | scope bullet added |
| 10 | MAJOR | Cheng–Wan regime mismatch (agreement O(k); list-decoding vs one-output) | accepted, and corrected further by the literature critic (§6.2 #1) |
| 11 | MAJOR | CHHS is full-modulus only; GSS uses coprime moduli and equality scores | accepted: §1.5 rows rewritten |
| 12 | MAJOR | "sparse supports ≤ u₀" reversed; Erickson transfer unjustified | accepted: corrected (≥ u₀ factors, ≤ log₂ 2X nonzero coordinates); transfer withdrawn; subset-sum row now "no route found (assessment)" |
| 13 | MINOR | BCPZZ/Jeronimo statements check; add field-size hypotheses | accepted |
| 14 | MINOR | KY/NP/BN titles and noise models | accepted (§1.3) |
| 15 | NIT | GV/GGG/BKR precise statements | accepted (§1.5) |
| 16 | MINOR | BGJT and DLI are heuristic | labelled |
| 17 | MINOR | "approximate counting does not help" too broad | narrowed to cumulative counts |
| 18 | MINOR | abstract costs; E2 exponent-cap change undisclosed | disclosed (§4.3) |
| 19 | MINOR | `JohnsonOracle.call` general-modulus contract false (X = 15, M = 30, \|R\| = 2) | fixed: unsupported calls raise; exercised calls unchanged; test added |
| 20 | MINOR | a 5-call query tree is missing | accepted: Prop DT′; verified post hoc (5 calls in all 200 windows) |
| 21 | MINOR | quantum search gives L[2/3, ν/2]; residue-guided sub-threshold beam search not covered | "classical" stated; beam search added to the ladder as OPEN |

## 6. Opus critics

### 6.1 Proof/cost critic (read-only): ACCEPT WITH FIXES, all applied
The critic checked every PROVED item step by step. That covered F2's ratios and the q-ary computation, F3's algebra,
F4, Prop DT, Lemmas NU/GC/SC, Prop PT, Prop SG and Prop PNP. It found that each step-2 answer has ≤ 1 certified point
and that no two step-1 prime powers fit under V. It reproduced every §4 number from the JSONL, re-ran the decoders
bit-identically on 8 windows, checked the accounting (budget checks, cost = tests + pointers, the exact `_int_pow_ceil`
threshold, no retests) and found no leakage. Nothing FATAL or MAJOR.

| # | Sev. | Finding | Fix |
|---|---|---|---|
| 1 | MINOR | Prop SG proved only for unweighted calls on the whole window; \|R\| = 1 is a counterexample to the unqualified form; weighted case not covered | qualified in §0, the ladder and §3.2; general-\|R\| form \|R\|·X^{−1/2+o(1)} added; weighted case OPEN |
| 2 | MINOR | "Fixed repertoires need ρ^{−1} calls" contradicts Prop DT; Prop AT bounds the portfolio size | reworded |
| 3 | MINOR | "no barrier possible" overstates Prop DT (it is a depth statement); probability bound short a factor (1 + h/2^{b−1}); hypothesis y ≤ h + 1 missing | all fixed |
| 4 | MINOR | leaf = sieve holds by construction, not as a measurement | relabelled |
| 5 | MINOR | P2 groups are decoder-pure; P2 = ratio of tests per success; biased under sequential stopping | disclosed in §4.1; §0 reworded |
| 6 | MINOR | Prop PT's ρ^{1−o(1)} cannot rank PH vs DT at L[1/3] scale | stated |
| 7 | NIT | minimum-weight word is the largest smooth one; Prop SI needs any equality case | F4 reworded |
| 8 | NIT | F3's o(1) terms not explicit; slow convergence | "leading order"; convergence numbers added |
| 9 | NIT | HEAVY-POINT input must include X | fixed |
| 10 | NIT | assessment inside Lemma SC | labelled |
| 11 | NIT | P3's RT < SV nearly automatic | noted |
| 12 | NIT | α_sf undefined | defined |
| 13 | NIT | DT8 explanation unlabelled | labelled POST-HOC |

Also adopted: Lemma GC holds for weighted scores; F2(c) labelled PROVED in MODEL.

### 6.2 Literature critic (read-only): OK WITH FIXES, all applied
Checked in primary texts:
- Cheng–Wan (arXiv math/0405082 and 0802.1220);
- GGG; GSS;
- BCPZZ and Jeronimo (arXiv HTML);
- the abstracts and listings of Wan–Zhang, CHHS, BHKS, Erickson, KLM, Kiayias–Yung, Bleichenbacher–Nguyen, Boneh,
  Guruswami 2003 and Cohn–Heninger.

Other items are second-hand, as marked.

| # | Sev. | Finding | Fix |
|---|---|---|---|
| 1 | MAJOR | Cheng–Wan's split polynomials lie in a residue class mod an irreducible (Prop DLI's analogue). Their received word is rational; Thm 2's agreement is ≈ 4k/3. Under T ↦ 1/T the note's translate is a residue class mod T^m, where DLP is easy. | §0, §1.4, §1.5, §3.3, F4 and §5.1 rewritten; CW removed as evidence for Prop SI |
| 2 | MEDIUM | the short-interval analogue lives in the deep-hole / Moments-Subset-Sum literature (Cheng–Murray, Li–Wan, GGG) | cited; "no known hardness or algorithm at d ≈ k², D = F_q" |
| 3 | MEDIUM | GGG regime; Wan–Zhang 2026 missing | updated |
| 4 | MEDIUM | CHHS is modulus-N only | fixed |
| 5 | MEDIUM | Guruswami 2003 has no capacity result; no integer CRT capacity decoder is known; BHKS covers polynomial ideals only | row rewritten |
| 6 | MEDIUM | RS precedent for F2(b): Justesen–Høholdt and CW's ĝ | cited; "no prior statement found by us" kept for the CRT version |
| 7–8 | MINOR | BCPZZ and Jeronimo verified; separate failure modes; BCPZZ's time is worse than the sieve | applied |
| 9 | MINOR | GRS/Boneh/GSS bibliography; GSS's threshold in nats | applied |
| 10 | MINOR | KY/NP/BN venues; noisy interpolation ≠ reconstruction; "above capacity" ambiguous | applied |
| 11–12 | MINOR | BKR, Guruswami–Rudra, Guruswami–Vardy venues | applied |
| 13 | MINOR | Erickson transfer gaps; Meyer auf der Heide | applied (transfer withdrawn anyway) |
| 14–15 | MINOR | Dyer, Hildebrand–Tenenbaum, de Bruijn attributions | de Bruijn/Granville added |
| 16 | MINOR | BGJT is heuristic; Joux, GGMZ, Kleinjung–Wesolowski | applied |
| 17 | MINOR | the Alon–Kilian DIVISOR credit is unconfirmed; cite the Shor/Barak answers | applied |
| 18 | MINOR | Cohn–Heninger "Ideal forms" (ICS 2011) is the natural Johnson-class reference | added to §1.5 |

## 7. Reproduction
```bash
.venv/bin/python -m pytest -q tests/test_crtdecode.py
.venv/bin/python -m src.crtdecode.run_exp E1     # ~15 s, 60 MB
.venv/bin/python -m src.crtdecode.run_exp E2     # ~45 s, 45 MB
.venv/bin/python -m src.crtdecode.analyze        # P1-P5
# Prop DT' (post hoc): src.crtdecode.decoders.crt_query_tree on the E2 windows (see §2.1)
.venv/bin/python -m src.crtdecode.theory         # a*, c(a), dictionary
```
