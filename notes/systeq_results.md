# systeq: is there an integer analogue of the systematic equation x^q − x = ∏(x − a)?

Agent `systeq`, phase 16 (route 1), 2026-10-20. This note is revised after codex systeq_2 and two Opus critics (§8).

**Files.**
- Code: `src/systeq/`.
  - `ffield.py`: table arithmetic for F_Q.
  - `bgjt.py`: the positive control, the post-hoc decomposition and the distinct-relation count.
  - `intanalog.py`: integer analogues and exact smoothness.
  - `barriers.py`: diagnostic checkers, including Prop F′ and the Descartes bound.
  - `run_experiment.py`.
- Tests: `tests/test_systeq.py` (16).
- Pre-registration: `notes/systeq_prereg.md` (commit 29d9ae8, made before the runs).
- Config: `experiments/configs/systeq.json`.
- Raw results:
  - pre-registered: `experiments/raw/systeq_E{1,2,3,4}.jsonl`;
  - post-hoc, disclosed: `systeq_E1_posthoc.jsonl`, `systeq_E1_distinct.jsonl`, `systeq_E2_nondegenerate.jsonl`.
- Codex: blind round `codex/systeq_1.{md,out}`; adversarial round `codex/systeq_2.{md,out}`.

**Labels.**
- **PROVED**: proof given here.
- **KNOWN**: published, with a citation.
- **HEURISTIC**: the assumption is named.
- **measured**.
- **refuted for this model**.
- **open**.
- **assessment**.

---

## 0. Verdict

**Plain language.** We did not find an integer analogue of the systematic equation that could push factoring below L[1/3]. For every concrete analogue we could formulate, we can say exactly which ingredient fails. The small-characteristic engine rests on one fact: *the q-th power map is additive, so a two-term polynomial X^q − X has q roots, and one fixed congruence X^q ≡ h₀/h₁ compresses a whole q³-element family of fully split polynomials.*

**The integer obstructions.**
1. **Prop F / F′ [PROVED].**
   - An additive power map z ↦ z^s on a faithful Z/N-algebra forces N to be squarefree and λ(N) | s − 1. Exhibiting such an s already factors N.
   - The restricted version is what BGJT actually uses: (X + b)^s ∈ span(1, X^s) for s + 1 constants b with unit differences, in a target of rank ≥ 4. It is impossible when every prime of N exceeds s.
   - Escapes exist only in rank ≤ 3 or with ≤ s constants (codex's examples, §5.1).
2. **Lemma CF [PROVED].** A systematic product that vanishes on all scalars (for example "∏(z − c) = σ(z) − z" for an automorphism σ) has degree ≥ p. The Wilson/factorial identity is this case. That is the obstruction behind factorial methods. Their N^{1/4+o(1)} cost comes from baby-step giant-step on top of it. Sparse products over a small μ_n *do* exist; they are handled by 3.
3. **Lemma DES [PROVED].** In characteristic 0, fewnomials have few real roots. The four-term BGJT/GKZ/KW space V_n = span(X^{n+1}, X^n, X, 1) has no nondegenerate member that splits over R once n ≥ 4. Over F_{q²} it has q³ + q.
4. **The auxiliary-algebra escape (codex) [PROVED part + HEURISTIC part].**
   - A free presentation A = (Z/N)[X]/(F) gives X^d ≡ h for free. But every compression identity whose reduction has height < N is an *exact* identity over Z[X]/(F̃), so it carries no information about N (Lemma EX).
   - Making it non-exact costs d conditions mod N on F. Under counting heuristic G_aux this pushes factor-base norms to N^{d/(d+1)}, worse than polynomial selection with a root (one condition, N^{1/(d+1)}).
   - Dense split families with an N-dependent compression are **open**: a Prouhet–Tarry–Escott-type problem (§6).

**Experiment (pre-registered; one prediction failed as stated, disclosed).**
- **Positive control, exact over F_{q²}, q = 4..16.** The PGL₂-orbit relations beat the size-bound prediction by G = 0.4, 2.5, 58, 275, 2.0·10³, 1.5·10⁵, 2.2·10⁷, **4.3·10¹⁰**. Distinct informative relations number 16, 37, 88, 103, 143, 245, 462, 752, against a factor base of Q + 1. The pre-registered headline G_FF (up to 2.3·10¹²) also counted the degenerate Frobenius-power family, which yields about 1/2 per form but only ≈ Q/2 distinct relations (§3.2).
- **The closest integer analogue.** This is the "F_{1^n}" V-family on Cunningham-type N = Φ_n(m): μ_n = ⟨m⟩ ⊂ Z/N is explicit, and m^n ≡ 1 plays X^q ≡ h₀/h₁. It gives **G = 1.03–1.06**, and **1.003–1.020** once the degenerate candidates are removed. That is no gain: the systematic side is an unstructured integer of size > N.
- Cyclotomic pieces come out *less* smooth than the per-piece baseline (R = 0.61 and 0.25; local-root effect).
- Codex's binomial collision test gave 10 collisions against 8.6 expected (UNVERIFIED by us).

**Verdict on route 1: refuted for models M1–M13 of §2**, under the barriers of §5. Open, and with no candidate known to us:
- dense split families with an N-dependent compression (§6);
- an elimination operation whose output size is independent of N.

Assessment: probability < 1% that route 1 yields sub-L[1/3].

**Against smoothenum's targets.**
- An output-sensitive algebraic-side generator tied to N would be worth (4/9)^{1/3}.
- A doubly-smooth generator would be worth polynomial time. The F_q systematic family *is* such a generator (one side split by construction, the other of bounded degree), which is consistent with quasi-polynomial DL.
- No integer candidate here is either. The V-family supplies a smooth-by-construction *small* side (|r_f| ≤ B always), but its other side is N-sized. So in smoothenum's terms θ = 0, and the family is Dixon-type rather than NFS-type.
- Prop LQ does not apply directly, since these are identities rather than divisibility selections. Our barriers (F′, DES, EX) play that role.

---

## 1. The descent engine, axiomatised (Step 1)

**Sources.**
- Joux 2013 (SAC 2013; ePrint 2013/095, L[1/4]).
- Göloğlu–Granger–McGuire–Zumbrägel (CRYPTO 2013; the X^{q+1} + aX + b family).
- BGJT 2014 (Eurocrypt 2014; arXiv 1306.4244).
- Granger–Kleinjung–Zumbrägel, "On the powers of 2" (ePrint 2014/300; degree-2 elimination).
- GKZ, "On the discrete logarithm problem in finite fields of fixed characteristic" (Trans. AMS 370 (2018); arXiv 1507.01495). Rigorous for infinitely many explicit fields with q > 61, q not a power of 4, and k ≥ 18. Covers all fields only under a representation conjecture.
- Kleinjung–Wesolowski (J. AMS 35 (2022), doi 10.1090/jams/985; arXiv 1906.10668; §§2–3.4 read by us).
- Lido (arXiv 2206.10327; small, not fixed, characteristic).
- Traps: Cheng–Wan–Zhuang (arXiv 1310.5124).
- Split counts: Bluher (FFA 10 (2004)).
- The name "systematic equation" goes back to Blake–Fuji-Hara–Mullin–Vanstone 1984. BGJT note that their usage differs.

KW put the core in one sentence: "All we need is a property similar to Equation (1): applying Frobenius is equivalent to a small degree rational function."

| Axiom | Content in F_{q^{2k}} |
|---|---|
| **E0 Size = degree** | Non-archimedean size with no carries. The q² linear polynomials form a factor base polynomial in log Q |
| **E1 Small constant field** | k = F_q with q = (log Q)^{O(1)}. T^q − T = ∏_{c∈k}(T − c) is a 2-term polynomial that splits completely |
| **E2 Additive Frobenius** | z ↦ z^q is a ring endomorphism fixing k: (aX + b)^q = a^σX^q + b^σ |
| **E3 Presentation freedom** | Target ≅ k′[X]/(I) with X^q ≡ h₀/h₁, deg h_i ≤ 2. *Existence* of such a model is heuristic for BGJT/GKZ in general (rigorous for GKZ's explicit families). KW make it rigorous with an elliptic model (Frobenius = translation by Q ∈ E(F_q)) after a degree-O(log n) extension. The freedom exists because the target depends only on its cardinality |
| **E4 Homographic family + compression** | V = span(X^{q+1}, X^q, X, 1). The PGL₂(k′)-orbit of X^q − X (q³ + q cosets mod PGL₂(F_q)) consists of squarefree split forms. Every f ∈ V reduces mod I to degree ≤ 3 |
| **E5 Size-halving elimination at cost poly(q)** | Degree 4 → 3 → 2 elimination over F_{q^{2^e}}, avoiding traps. The splitting probability in V is 1/O(q³) *independent of the extension*. Norms go down the tower (zigzag); depth O(log n) gives (qn)^{O(log n)} |

E0–E4 *together with classical descent* give L[1/4] (Joux). Quasi-polynomial time needs E5. Codex's axiomatisation (`systeq_1.out`, E1–E6) matches this.

**Measured split structure of P(V)(F_{q²}) [measured; agrees with Bluher-type counts].** The split members are exactly:
- the squarefree PGL₂ orbit (q³ + q forms);
- all forms ℓ₁^q·ℓ₂ with ℓ₁, ℓ₂ ∈ P¹(F_{q²}), i.e. the Frobenius-power family (q² + 1)².

Their total is (q² + 1)² + q³ + q for every q in E1 (Opus critic; also 357, 806, …, 70161 in our tables).

## 2. Integer analogues, model by model (Step 2)

| # | Model / candidate | Axioms supplied | Smooth supply | Tie to N | Verdict |
|---|---|---|---|---|---|
| M1 | Cyclotomic X^n − 1 = ∏ Φ_d(X) | E1 only over Z[ζ_n] (μ_n, "F_{1^n}") | Pieces Φ_d(t) of log-size φ(d) log t; gain ≤ n/φ(n) (§5.4). Zsigmondy: not B-smooth for n ≥ B. E4: R = 0.61 and 0.25 (local-root effect) | Only for N \| Φ_n(m); then the systematic product is the SNFS algebraic norm Φ_n(t, s) | refuted for general N. For special N we know nothing better than SNFS (assessment) |
| M2 | Aurifeuillian | 2 pieces | Factor 2 in log-size | Special N | constant only (snfsrig §3) |
| M3 | Consecutive products / binomials | Smooth by construction (C(n, k) is n-smooth) | Abundant, but products of factor-base primes are already free | None; exact identities (codex Barrier 1). Codex: 10 vs 8.64 collisions (UNVERIFIED) | refuted (no excess yield) |
| M4 | Factorials mod N | Lemma CF: ∏_{c<p}(z − c) = σ(z) − z mod p | — | Degree ≥ p | N^{1/4+o(1)} (Pollard 1974, Strassen 1976, BGS 2007). Deterministic records are now N^{2/9} (Hittmeir 2018) and N^{1/5} (Harvey 2021) |
| M5 | Chebyshev/Dickson composition | Multiplicative iterate; cyclotomic in disguise | As M1 | p + 1 torus method | group-order, not relations |
| M6 | EDS / division polynomials | E1 over Q(E[n]) | Primitive part ≳ (6/π²) log\|W_n\| asymptotically (Everest–McLaren–Ward 2006). Primitive divisors for all but finitely many n (Silverman 1988, ineffective) | ECM, L_p[1/2, √2] (Lenstra 1987) | group-order |
| M7 | Lucas sequences | Cyclotomic over a quadratic field | Primitive divisors for n > 30 (Bilu–Hanrot–Voutier 2001) | p ± 1; Φ_k(p) methods (Bach–Shallit 1989) | group-order |
| M8 | Cyclotomic norms, Gauss/Jacobi sums | Automorphisms σ_a of (Z/N)[ζ] over every N. Jacobi sums give smooth norms ℓ^{d/2} | Explicit small ideals | Lemma CF if identities must hold on all scalars. A map to Z/N needs a root of Φ_m mod N | NFS in cyclotomic fields (RES/G). Mod N, Jacobi sums are the APR-CL primality test, which constrains the prime divisors r ≡ N^i (mod s) (Adleman–Pomerance–Rumely 1983; Cohen–Lenstra 1984; Lenstra 1985) |
| M9 | p-adic / Witt lift | Teichmüller: X^p − X has p roots in Z_p | — | Z/p^s → Z/N exists iff N \| p^s | refuted (r20 #6) |
| M10 | Carmichael field equation x^{λ(N)+1} ≡ x | The genuine field equation of Z/N | — | Degree λ(N); knowing it factors N | equivalent to the problem |
| M11 | Λ-rings / Adams ψ_q (substitution Frobenius) | ψ_q(T) = T^q on Z[T] | — | Lemma SUB: σ_q(f)g − fσ_q(g) = det·(T^q − T): the homography trick gives one relation (tested) | refuted |
| M12 | Integer systematic equations, q ≤ 3 | T³ − T = (T − 1)T(T + 1); holds over Z only for q ≤ 3 (tested q ≤ 8) | 4 pieces | No additivity | E4 fails |
| M13 | **The "F_{1^n}" V-family on Cunningham N = Φ_n(m)** (§3) | E1 (μ_n = ⟨m⟩ ⊂ Z/N), E3 (m^n ≡ 1), the *small side* of E4 (f(m) ≡ (A + C)m + (B + D)) | The split side of E4 fails (Lemma DES); f(m) is an unstructured integer > N | Special N | measured G ≈ 1.00–1.02; refuted |
| M14 | Auxiliary algebra (Z/N)[X]/(F) with F free (codex) | E3 for free (X^d ≡ h) | Exact lifts: smooth but uninformative (Lemma EX). N-informative lifts: norms ≥ N^{d/(d+1)} (G_aux) | — | refuted under G_aux. The dense/PTE version is **open** (§6) |

## 3. The best candidate (M13), developed and tested (Step 3)

### 3.1 Construction and cost

**Construction.**
- Let N | Φ_n(m). Then m has order n mod N, and m^n ≡ 1 is the analogue of X^q ≡ h₀/h₁.
- For f = A x^{n+1} + B x^n + C x + D:

  f(m) ≡ r_f := (A + C)m + (B + D) (mod N).

  This is the compression half of E4: |r_f| ≤ 2H(m + 1) < B in every E2 instance, so this side is smooth by construction.

**What fails.** By Lemma DES, for n ≥ 4 the only members of V_n that split over R are x^n(Ax + B). Their "relation" only restates m^n ≡ 1. So f(m) is unstructured, with f(m) ≈ H·m^{n+1−φ(n)}·N > N.

**Cost [HEURISTIC: smoothness of linear-form values as for random integers].** The relations are Dixon-type, with one side of size N·m^{n+1−φ(n)}·H. That gives cost ≥ L_N[1/2, √2]. This assumes fast smoothness tests and sparse linear algebra (Pomerance 1982; Croot–Granville–Pemantle–Tetali 2012). Dixon's 1981 analysis gives L[1/2, 2√2].

The same structure used as an NFS is SNFS, L[1/3, (32/9)^{1/3} = 1.526]. For general N, small-height m, h with m^n ≡ h (mod N) violate heuristic G.

### 3.2 Experiments vs pre-registration

**E1 (positive control; exact enumeration of P(V)(F_{q²})).**
- Size prediction: P_rand(q+1)·P_rand(3), the exact fraction of split binary forms.
- Every relation was verified as h₁f ≡ r_f (mod I).

| q | \|P(V)\| | split | relations (raw) | G_FF raw (prereg) | G_FF* | orbit: rel / form | G_orbit | distinct informative rel. | Q + 1 |
|---|---|---|---|---|---|---|---|---|---|
| 4 | 4369 | 357 | 160 | 9.1 | — | 7/68 | 0.40 | 16 | 17 |
| 5 | 16276 | 806 | 362 | 38 | — | 24/130 | 2.5 | 37 | 26 |
| 7 | 120100 | 2850 | 1362 | 1.27e3 | — | 62/350 | 58 | 88 | 50 |
| 8 | 266305 | 4745 | 2215 | 8.7e3 | 2.7e3 | 70/520 | 275 | 103 | 65 |
| 9 | 538084 | 7462 | 3464 | 7.0e4 | 2.0e4 | 102/738 | 2.0e3 | 143 | 82 |
| 11 | 1786324 | 16226 | 7626 | 6.4e6 | 1.5e6 | 184/1342 | 1.5e5 | 245 | 122 |
| 13 | 4855540 | 31110 | 14996 | 8.6e8 | 1.7e8 | 376/2210 | 2.2e7 | 462 | 170 |
| 16 | 16843009 | 70161 | 33776 | 2.3e12 | 3.8e11 | 623/4112 | 4.3e10 | 752 | 257 |

**Evaluation against the pre-registration.**
- (P1a) |log₁₀ G_FF − log₁₀ G_FF*| ≤ 1: held as stated (deviations 0.50–0.79 decades). *Disclosed:* it held because two errors offset. The 8/q³ model overstates the orbit about 8-fold, and the degenerate family (X + b)^q(X + c) adds about q-fold. That family yields ≈ 1/2 per form but cancels to ≈ Q/2 distinct Frobenius relations. G_orbit lies about one decade below G_FF*, at the tolerance edge.
- (P1b) Monotone: held, both raw and orbit-only.
- (P1c) Relations ≥ Q + 1: held for raw counts. For *distinct* relations (post-hoc) it holds for every predicted q (8–16) but fails at q = 4 (16 < 17; q = 4 was not predicted).
- (P1d) All relations verified: held.
- Rank was not computed.

**E2 (integer V-family; 200 000 candidates; B = 2^16, H = 32).**

| (n, m) | N bits | G (prereg) | degenerate removed | G (non-degenerate) |
|---|---|---|---|---|
| (5, 520) | 37 | 1.047 | 64 | 1.018 |
| (7, 130) | 43 | 1.034 | 52 | 1.003 |
| (11, 19) | 43 | 1.040 | 251 | 1.003 |
| (13, 11) | 42 | 1.059 | 433 | 1.020 |

- (P2a) G ∈ [0.5, 2]: held.
- (P2b) The kill signal G ≥ 4 was not triggered. But since r_f is always smooth, E2 tests only the f(m) side.
- (P2c) "f-side log-size ≥ log N": **failed as stated**. The degenerate candidates Am + B = 0 give f(m) = r_f (codex), and Cm + D = 0 gives f(m) = m^n r_f (Opus critic). These are all automatically "smooth", and they caused the +2–6% excess. The note's earlier "log-uniform within bins" explanation was wrong. After removing them, P2c holds (min |f(m)| ≥ 2^45 > N) and G = 1.003–1.020.
- The first N has the small factor 11. That does not affect yield.

**E3 (Descartes; exhaustive with H = 4).** Nondegenerate split members of V_n(Z):

| n | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|
| count | 208 | 16 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

(P3: 0 for n ≥ 7) held. Lemma DES proves 0 for n ≥ 4.

**E4 (cyclotomic calibration).**
- R_piece = 113/184.8 = 0.61 for n = 6, and 8/31.9 = 0.25 for n = 12. (P4: [0.1, 1.5]) held.
- *Interpretation (corrected):* this is the local-root effect. For example, Φ₁₂(m) is never divisible by 2 or 3, and eligible primes p ≡ 1 (mod d) occur with multiplicity φ(d). It is not Zsigmondy, which bites only for n ≥ B.
- The baseline multiplies marginal rates of pieces at the same m, which assumes independence. n = 12 rests on 8 successes. Treat both rows as finite calibrations, not as a universal bound.

## 4. Why "smooth by construction" is not the bottleneck in Z

Smooth integers by construction are free. What BGJT supplies is a *congruence* between a free-smooth object and a small one, valid for a large family at once.

- Codex Barrier 1 [PROVED]: products of small integers that are both below N and congruent mod N are equal.
- Codex's aggregate bound [PROVED]: if U ≡ V (mod N) and U ≠ V, then max(|U|, |V|) ≥ N/2. It still allows one side to be assembled from many small known factors, and that is exactly the BGJT shape the barriers of §5 address.

## 5. Barriers (Step 4)

### 5.1 Prop F and Prop F′ [PROVED; elementary. Cf. the generalized Korselt criterion (Korselt 1899, Carmichael 1912), AKS Lemma 2.1, Howe "Higher-order Carmichael numbers" (2000), Lenstra 1985. The use as a BGJT-transfer obstruction: to our knowledge not stated before]

**Prop F.** Let A be a faithful commutative Z/N-algebra and s ≥ 2 with (x + y)^s = x^s + y^s on A.
- The map z ↦ z^s is then a unital ring endomorphism, so it fixes every integer scalar: a^s ≡ a (mod N) for all a ∈ Z.
- By generalized Korselt this holds iff N is squarefree and λ(N) | s − 1. Directly: p² | N would give p^s ≡ p (mod p²), and on F_p one gets (p − 1) | (s − 1).
- A multiple of λ(N) factors N in randomized polynomial time (folklore built on Miller–Rabin; e.g. Crandall–Pomerance §5, Bach–Miller–Shallit 1986). Miller 1976 gives the deterministic ERH version.
- Tests: on (Z/15)[t]/(t² + 1) the additive s ≤ 40 are exactly the set predicted by the component Frobenius congruences, which is stronger than Prop F's necessary condition. On Z/21 they are s ≡ 1 (mod 6).

**Prop F′ (restricted compatibility; this is what E2/E4 of BGJT use).**
- Let k′ ⊇ Z/N be a nonzero commutative ring and A = k′[X]/(F) with F monic of degree d ≥ 4.
- Let s ≥ 3 with s, C(s, 2), C(s, 3) units in k′. This holds if every prime of N exceeds s.
- Suppose (X + b_j)^s ∈ k′ + k′X^s for s + 1 constants b_j ∈ k′ with pairwise differences units.
- Then we have a contradiction.

*Proof.*
- Q(T) := (X + T)^s = Σ_i C(s, i)X^i T^{s−i} is a degree-s polynomial in T with coefficients in A.
- The Vandermonde matrix of the b_j is invertible over k′. Inverting it shows each coefficient C(s, i)X^i lies in the k′-module M = k′ + k′X^s.
- Hence 1, X, X², X³ ∈ M. Choosing preimages in k′² of these four elements gives an injective map k′⁴ → k′² (1, …, X³ are free in A), because composing with k′² → A recovers the injective map k′⁴ → A.
- This contradicts McCoy's theorem (an injective k′-linear map k′^m → k′^n forces m ≤ n). ∎

**Scope.**
- Cross-reference: `nonsmooth` Prop FROB (phase 16) proves the companion statement that a Frobenius *operator* on an étale Z/N-algebra already factors N.
- Codex's escapes satisfy F′'s exclusions:
  - R[X]/(X² − X − 1), where every (aX + b)² ∈ span(1, X²), has rank 2;
  - R[X]/(C_ℓ) with C_ℓ = ((T + 1)^ℓ − T^ℓ − 1)/ℓ, where (aX + b)^ℓ = a^ℓX^ℓ + b^ℓ for a, b ∈ {0, 1}, uses ≤ 2 constants.
- Both are tested (`test_prop_F_prime_restricted`; rank 4–6 random F give only b = 0).
- With ≤ s constants the homographic family has ≤ s^{O(1)} members. Whether that could suffice for a descent with polylog s is not excluded (**open**). But E5 also needs compatibility for every later target, over growing extensions.
- If s ≥ p_min, then gcd(s!, N) already reveals a factor, and any systematic product has degree ≥ p (Lemma CF).

### 5.2 Lemma CF [PROVED; elementary root covering]

**Statement.** Let A be a Z/N-algebra, S ⊂ A finite and P(T) = ∏_{c∈S}(T − c). If P vanishes on Z·1_A, then the image of S in each residue field of characteristic p contains F_p, so |S| ≥ p.

**Application.** This covers automorphism-difference identities P(z) = σ(z) − z, and more generally any bounded-degree expression E(z, σz) that holds on all scalars (codex). A unital σ fixes the scalars, so P vanishes there.

**Instance (tested).** ∏_{c=0}^{p−1}(z − c) = σ(z) − z in F_p[√D]. This is the Wilson/factorial shadow of X^p − X.

**Not covered.** Restricted identities, e.g. T^{n+1} − T = T∏(T − ζ^j) in R[ζ_n] evaluated at z = uζ^j, giving P(z) = (u^n − 1)z (codex). Those are sparse fewnomial identities, handled by §5.3 and §5.5.

**KNOWN context.** Strassen 1976 and Blum–Cucker–Shub–Smale 1998 (as cited by Bürgisser, Comput. Complexity 18 (2009)): short straight-line programs for n! would give *nonuniform* polynomial-time factoring. Lipton 1994 (ANTS-I) relates average-case hardness of factoring to a weak τ-conjecture. The function-field counterpart is that X^{q^k} − X has an O(k log q) SLP, which is distinct-degree factorisation.

### 5.3 Lemma DES [PROVED from KNOWN]

**Statement.**
- By Descartes' rule, applied to f(x) and f(−x), f = A x^{n+1} + B x^n + C x + D with all terms nonzero has at most 3 (n even) or 4 (n odd) nonzero real roots (`sign_variation_bound`).
- If D = 0 ≠ C, then f = x·g with g a trinomial, which has at most 2 or 3 nonzero real roots.
- So for **n ≥ 4** every member of V_n(R) with A ≠ 0 that splits over R is x^n(Ax + B). This strengthens the original n ≥ 6, following codex.

**Over number fields.** Lenstra 1999 ("On the factorization of lacunary polynomials", Number Theory in Progress I, de Gruyter, 277–291) bounds the total degree, with multiplicity, of the irreducible factors of degree ≤ d of a lacunary f with f(0) ≠ 0 over K by a function of [K : Q], the number of terms, and d. Roots are the d = 1 case.

**Prime-field analogue (KNOWN).** Bi–Cheng–Rojas (ISSAC 2013) and Cheng–Gao–Rojas–Wan (arXiv 1411.6346) note that sparse polynomials with many roots are easy over F_{p^k} but no explicit ones are known over prime fields F_p. Prime fields thus share part of the characteristic-0 obstruction, consistent with prime-field DL being stuck at L[1/3].

**Over Q_p.** X^p − X does have p roots (Teichmüller), but only at the hidden prime (M9).

### 5.4 Splitting Lemma (integer values only) [PROVED / KNOWN]

**Statement.** For F ∈ Z[x], t ∈ Z and α a root of F in a number field K,

N_{K/Q}(t − α) = ±minpoly_α(t)^{[K:Q(α)]}, up to leading-coefficient factors.

So on *integer values*, a factorisation of F over any K yields nothing finer than the Q-irreducible factorisation. The *ideal* factorisation of t − α in O_K is finer, and that is what NFS uses.

**Consequences.**
- *Cyclotomic:* the gain in log-size is ≤ n/φ(n) ≤ e^γ log log n + O(1/log log n) (Rosser–Schoenfeld 1962, Thm 15). Since n ≤ log N here, this is O(log log log N) inside the L-constant: harmless, though not zero in principle.
- *EDS:* the gain is ≤ π²/6, asymptotically (Everest–McLaren–Ward 2006).
- *Zsigmondy (Bang 1886; Zsigmondy 1892):* m^n − 1 has a primitive prime ≡ 1 (mod n) except for (2, 1), (2, 6), and n = 2 with m + 1 a power of 2. So m^n − 1 is not B-smooth for n ≥ B ≥ 7 (tested).

### 5.5 Auxiliary algebras and rigidity

**F_q side [KNOWN].** The BGJT/KW freedom is the freedom to choose I, because all fields of one size are isomorphic. Model *existence* is heuristic except in KW and GKZ's explicit families.

**Lemma RES [PROVED, snfsrig §1.1].** It governs presentations Z[x]/(f) → Z/N through a root mod N. Heuristic G then gives heights N^{1/(d+1)} (heurnfs B1).

**Auxiliary algebras A = (Z/N)[X]/(F), F chosen freely (codex systeq_2 §4).**
- *Lemma EX (exactness) [PROVED].* Let F̃ ∈ Z[X] be a monic lift of F and u, v ∈ Z[X] with u ≡ v mod (N, F̃). Write R₀ = (u − v) mod F̃ ∈ Z[X]. Then N | R₀ coefficientwise. So either R₀ = 0, in which case u = v in Z[X]/(F̃) and the norm valuation rows coincide (zero information), or ‖R₀‖_∞ ≥ N.
  - Codex's examples illustrate both branches. P_M − X − 1 gives an exact identity with smooth, uninformative norms (the factorization of M!). P₄ − X − 1 + N gives an informative one with N-sized norms.
- *G_aux [HEURISTIC; counting in the style of G].* Impose a Frobenius-type compression X^q ≡ h(X) on A that is *not* exact over Z, i.e. R₁ = (X^q mod F̃) − h ≡ 0 (mod N) with R₁ ≠ 0. That is d conditions mod N on the d + 1 coefficients of F̃ (height H_F) and on h (height H_h). Generic solutions need H_F^{d+1}H_h^{d+1} ≳ N^d, so the factor-base norms F̃(c) ≈ H_F c^d are ≳ N^{d/(d+1)}. Polynomial selection through a root (one condition) gives N^{1/(d+1)}. Free presentations therefore trade one condition for d and are worse for every d ≥ 2. For d = 1 they *are* polynomial selection.
- *What G_aux does not cover:* families in which the compression is not imposed by one fixed congruence. See §6.
- Splitting A itself need not factor N: (1 + X)/2 is an idempotent of R[X]/(X² − 1) for every odd N (codex). A reduction must extract a zero divisor of the base ring.

### 5.6 E5 has no counterpart we can construct [HEURISTIC]

**Function fields.** The advantage of E5 is the splitting probability 1/O(q³) in V, *independent of the extension* F_{q^{2^e}}, together with a small side of degree (1 + δ)D independent of the target degree k. Base change by itself lowers no absolute size in either setting: a degree-D irreducible becomes D linear factors of norm q^D each. The earlier "[PROVED] base change" remark is withdrawn.

**Integers.** Under RES/G every candidate tied to N has log-norm product ≥ 2√(2·log N·log E) (heurnfs B1). The bottom level of any descent then costs L[1/3] (heurnfs B5, under H + HD + T).

### 5.7 The F_1 dictionary [PROVED computation; interpretation labelled]

- Over "F_{1^n}" = μ_n ∪ {0}, GL₂ = μ_n ≀ S₂ = μ_n² ⋊ S₂ (Kapranov–Smirnov, unpublished preprint c. 1995; Le Bruyn arXiv 1304.6532 §2). So |PGL₂| = 2n.
- This group sends T^n − 1 only to unit multiples of T^n − 1 and 1 − T^n: O(1) relations, against q³ + q.
- *Interpretation:* additivity is what enlarges the monomial group to PGL₂. This concerns the F_1 group only, not every homography over the real coefficient ring.

## 6. What it would take (the open escape)

A route-1 success needs a family 𝓕 of pairs (u, v) with u ≡ v (mod N), u ≠ v, such that:
- (a) u is a product of ≳ log N/log B pieces, each ≤ B^{O(1)}, by an identity;
- (b) v has log-size o(x^{2/3}), or at least independent of log N;
- (c) |𝓕| ≥ L[α] with α < 1/3, computable without the factors.

**Excluded:**
- Frobenius-type power compatibility on ≥ s + 1 constants in rank ≥ 4 (F′);
- systematic products valid on all scalars (CF);
- sparse four-term spaces (DES);
- one fixed non-exact compression in an auxiliary algebra (EX + G_aux);
- polynomial selection (RES + G).

**Not excluded (open).** *Dense* split families, many split polynomials ∏_j(X − a_j) with small a_j, sharing a common reduction to a small side modulo some presentation (N, F). Over Z this asks for many integer multisets of size M with prescribed elementary symmetric functions e₁, …, e_{M−4}. That is a multigrade Prouhet–Tarry–Escott problem, whose ideal solutions are known only for M ≤ 12, and it would additionally need the N-dependence. Codex's Dickson example (X² − 4)(X² − 1)² ≡ X + 3 is one such split polynomial, exact over Z.

We have no construction and no barrier here (assessment: implausible, since F_q gets the family from the PGL₂ symmetry of a 2-nomial, which has no characteristic-0 counterpart by §5.7). It is the precise remaining target for route 1.

## 7. Codex rounds

### 7.1 Blind round (`systeq_1.out`, before our analysis)

**Same verdict**, < 1% for a direct transfer.

**Codex contributions:**
- E1–E6;
- the coset count q(q² + 1);
- Z/p^s → Z/N iff N | p^s;
- Λ-ring equivariance density ≤ p²/N;
- T^p − T = ∏(T − a) + pC_p(T);
- Jacobi-sum norms ℓ^{d/2};
- a⁴ + 4b⁴;
- the binomial collision kill test: 10 collisions against 8.635 expected (random control 8). The script is in git-ignored scratch and was not re-run by us, so it is UNVERIFIED.

### 7.2 Adversarial round (`systeq_2.out`)

**Major points, all applied:**
- the barriers were stated too universally;
- partial power compatibility escapes exist, motivating F′ and its scope;
- the CF classification was overstated;
- DES strengthens to n ≥ 4;
- the auxiliary-algebra escape with norm accounting, motivating Lemma EX, G_aux and §6;
- E1 counts redundant Frobenius-family relations;
- E2's P2c fails because of the Am + B = 0 degenerates, and r_f is always smooth;
- E4 relies on an independence assumption and its exceptions;
- Prop F context in Howe 2000;
- the μ_n-set convention in Le Bruyn.

## 8. Opus critics

### 8.1 Proof / cost critic: ACCEPT WITH FIXES (applied)

| Issue | Severity | Fix applied |
|---|---|---|
| Prop F's global hypothesis does not match the BGJT use | MAJOR | Prop F′ added; verdict scoped to M1–M14 |
| §0 statement of CF was false (μ_n products exist) | MAJOR | restated |
| The E1 headline counted redundant and tautological relations | MAJOR | orbit-only G and distinct counts; P1a/P1c disclosure |
| The E2 bias explanation was wrong (degenerate Cm + D = 0) | MAJOR | measured, excluded; r_f always smooth noted |
| §5.6 "PROVED" said nothing that separates the settings | MAJOR | withdrawn; function-field sentence added |
| Assessments were labelled PROVED (SNFS-best, anti-smooth) | MAJOR | relabelled; local-root effect separated from Zsigmondy |
| Faithful ⇒ squarefree should be derived, not assumed | MINOR | derived |
| Splitting Lemma should be stated for values only; (log log)^{O(1)} remark | MINOR | stated |
| Dixon constant should be "≥ √2" | MINOR | stated |
| E3 caveat; L[1/4] needs classical descent; exact split count | MINOR | stated |
| P2c and P3 bookkeeping; UNVERIFIED codex test; Lenstra citation | MINOR | done |

### 8.2 Literature critic (applied)

- Corrected the GKZ titles and scope; added GGMZ, Bluher, Cheng–Wan–Zhuang, Blake et al., Lido.
- Gave KW's J. AMS venue.
- Corrected Lipton / Shub–Smale: Strassen 1976 and BCSS; nonuniform.
- Added the factorial records (Hittmeir, Harvey).
- Corrected the Miller attribution.
- Gave the Lenstra d = 1 case; Zsigmondy exceptions and Bang; Silverman ineffective; Everest–McLaren–Ward; APR-CL wording; Bach–Shallit.
- Reclassified Prop F as elementary/known in substance (Korselt, AKS, Lenstra 1985, Howe).
- Prior informal statement of the mechanism: Granger–Kleinjung–Zumbrägel, "Indiscreet logarithms in finite fields of small characteristic" (AMC 12 (2018), arXiv 1604.03837, §7): the new algorithms "exploit … the splitting of X^q − X"; for prime fields "fundamentally new ideas will be required".
- Added Bi–Cheng–Rojas and Cheng–Gao–Rojas–Wan.
- No published integer or F_1 analogue of BGJT descent was found.

No codex confirmation round was run on this revision.

## 9. Reproduction

```bash
.venv/bin/python -m src.systeq.run_experiment all       # E1, E1_posthoc, E2, E3, E4 (~5 min)
.venv/bin/python -m pytest -q tests/test_systeq.py
```

The post-hoc E1 distinct counts come from `bgjt.distinct_relations` over the E1 fields. The non-degenerate E2 run is `intanalog.vfamily_yield(..., exclude_degenerate=True)` with the E2 config and seed. Both are `run_experiment.py E1_distinct` and `run_experiment.py E2_nondeg`, which write the same files.
