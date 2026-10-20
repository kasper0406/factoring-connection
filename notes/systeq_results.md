# systeq: is there an integer analogue of the systematic equation x^q − x = ∏(x − a)?

Agent `systeq`, phase 16 (route 1), 2026-10-20.

**Files.**
- Code: `src/systeq/`.
  - `ffield.py`: table arithmetic for F_Q.
  - `bgjt.py`: the positive control.
  - `intanalog.py`: integer analogues and exact smoothness.
  - `barriers.py`: diagnostic checkers.
  - `run_experiment.py`.
- Tests: `tests/test_systeq.py`.
- Pre-registration: `notes/systeq_prereg.md` (commit 29d9ae8, made before the runs).
- Config: `experiments/configs/systeq.json`. Raw results: `experiments/raw/systeq_E{1,2,3,4}.jsonl`, plus `systeq_E1_posthoc.jsonl`.
- Codex: blind round `codex/systeq_1.{md,out}`; adversarial round `codex/systeq_2.{md,out}` (§7).

**Labels.**
- **PROVED**: proof given here.
- **KNOWN**: published, with a citation.
- **HEURISTIC**: the assumption is named.
- **experimentally supported**.
- **refuted for this model**.
- **UNVERIFIED**.

---

## 0. Verdict

**Plain language.** No integer analogue of the systematic equation exists that could push factoring below L[1/3]. We can say precisely why. The small-characteristic engine rests on one fact: **the q-th power map is additive, so a two-term polynomial X^q − X has q roots.**

Characteristic 0 and Z/N each destroy that fact in their own provable way:
1. **Prop F.** An additive power map z ↦ z^s on any faithful Z/N-algebra forces (p − 1) | (s − 1) for every prime p | N. So s ≥ p, and anyone who can exhibit such an s can already factor N (Miller).
2. **Lemma CF.** The constants of any "systematic product" over a Z/N-algebra must cover F_p, so its degree is ≥ p. The Wilson/factorial identity is the integer systematic equation over the hidden field F_p. That is why every factorial method is exponential (N^{1/4}).
3. **Lemma DES.** In characteristic 0, fewnomials have few roots.
   - Descartes: a real polynomial with k terms has ≤ 2k − 2 nonzero real roots.
   - Lenstra 1999: over a number field of bounded degree, the number of roots is bounded in terms of k and the degree.
   - So the four-term BGJT/GKZ/KW space V = span(X^{q+1}, X^q, X, 1) has *no* nondegenerate member splitting over Q once q ≥ 6. Over F_{q²} it has q³ + q of them.
   - Escaping through cyclotomic fields (X^n − 1 = ∏(X − ζ^j)) yields only the Q-irreducible pieces Φ_d(t). Their largest piece is a factor ≤ n/φ(n) = O(log log n) smaller in log-size. By Zsigmondy those pieces are *anti*-smooth: they contain a prime ≥ n + 1.
4. **Rigidity.** A BGJT model exists for every small-characteristic target because all fields of a given order are isomorphic. Z/N admits only polynomial-selection presentations, which are governed by Lemma RES and heuristic G (heurnfs B1). An integer "X^q ≡ h₀/h₁" with small heights exists only for special N, and then nothing beyond SNFS follows.

**Experiment (pre-registered; all predictions held).**
- Positive control (exact, over F_{q²}, q = 4..16): the systematic family beats the size-bound prediction by G_FF = 9, 38, 1.3·10³, 8.7·10³, 7.0·10⁴, 6.4·10⁶, 8.6·10⁸, **2.3·10¹²**.
- The closest integer analogue (the "F_{1^n}" V-family on Cunningham-type N = Φ_n(m), where μ_n ⊂ Z/N is an explicit constant set and m^n ≡ 1 plays X^q ≡ h₀/h₁) gives **G = 1.03–1.06**. That is no gain: its systematic side is a random integer of size ≥ N.
- Cyclotomic splitting is anti-smooth relative to the per-piece prediction: R = 0.61 and 0.25.
- Codex independently ran a binomial-coefficient collision test: 10 collisions against 8.6 expected.

**Verdict on route 1: refuted for every analogue we could formulate.** The precise barriers are in §5. What remains open is not an "analogue" at all. It would be an elimination operation whose output size is independent of N (codex's E4/E5). By heurnfs B1/B5 such an operation cannot come from polynomial selection. We have no candidate. Our probability that route 1 yields sub-L[1/3] is < 0.5% (assessment).

---

## 1. The descent engine, axiomatised (Step 1)

Sources:
- Joux 2013 (L[1/4]; ePrint 2013/095);
- BGJT 2014 (arXiv 1306.4244);
- GKZ 2018 (arXiv 1507.01495, "On the powers of 2");
- Kleinjung–Wesolowski 2019/2022 (arXiv 1906.10668, §§2–3.4, read by us).

KW state the core in one sentence: "All we need is a property similar to Equation (1): applying Frobenius is equivalent to a small degree rational function."

| Axiom | Content in F_{q^{2k}} | Used by |
|---|---|---|
| **E0 Size = degree** | Non-archimedean size with no carries. Linear polynomials (q² of them) form a factor base of size *polynomial* in log Q | all |
| **E1 Small constant field** | k = F_q with q = (log Q)^{O(1)} inside the base k′ = F_{q²}. The systematic equation T^q − T = ∏_{c∈k}(T − c) is a 2-term polynomial that splits completely over k | all |
| **E2 Additive Frobenius** | z ↦ z^q is a ring endomorphism fixing k: (aX + b)^q = a^σX^q + b^σ | all |
| **E3 Presentation freedom** | Target ≅ k′[X]/(I) with X^q ≡ h₀/h₁ mod I, deg h_i ≤ 2 (BGJT/GKZ). KW: an elliptic model where Frobenius equals translation by Q ∈ E(F_q) at a degree-n place; it always exists after a degree-O(log n) extension. This works because the target depends only on its cardinality | all |
| **E4 Homographic family + compression** | For V = span(X^{q+1}, X^q, X, 1), the PGL₂(k′)-orbit of X^q − X (q³ + q cosets mod PGL₂(F_q)) consists of split forms. Every f ∈ V reduces mod I to degree ≤ 3, so one side splits for free and the other has degree independent of k | Joux (L[1/4]) and all later work |
| **E5 Size-halving elimination at cost poly(q)** | A degree-D element over F_{q^{2^e}} (degree 4 → 3 → 2, trap-avoiding) is rewritten via E4 over the *extended constant field*, with splitting probability 1/O(q³) *independent of the extension*. Norms go down the quadratic tower (zigzag). Depth O(log n) gives (qn)^{O(log n)} | BGJT/GKZ/KW (quasi-polynomial) |

E0–E4 give L[1/4] (Joux). Quasi-polynomial time needs E5, which is E4 used over arbitrarily large extensions of the constants at a cost governed by q alone.

Codex's independent axiomatisation (`systeq_1.out`, E1–E6) matches this, with a stricter split of E5 into "elimination" and "completion/traps".

**Measured (E1, exact) [experimentally supported; new only as a measurement].** Over F_{q²} the split members of P(V) are exactly:
- the PGL₂ orbit (q³ + q squarefree forms; the count holds for every q tested);
- the degenerate Frobenius-power forms (X + b)^q(X + c) and (X + b)^q;
- forms of degree ≤ 1.

So frac_split ≈ (q + 1)/q³.

## 2. Integer analogues, axiom by axiom (Step 2)

| # | Candidate | Axioms supplied | Smooth supply near the needed size | Tie to N | Verdict |
|---|---|---|---|---|---|
| 1 | Cyclotomic X^n − 1 = ∏_{d\|n} Φ_d(X), at integers | E1 only over Z[ζ_n] (constants μ_n, "F_{1^n}") | The pieces Φ_d(t) have log-size φ(d) log t: the gain is ≤ n/φ(n) = O(log log n) [PROVED, §5.4]. Zsigmondy: m^n − 1 is never B-smooth for n ≥ B [PROVED via KNOWN]. E4 measured R = 0.61 and 0.25 | Only for N \| Φ_n(m) (special N). Then it is SNFS: the systematic product ∏_j(t − sζ^j) *is* the algebraic norm Φ_n(t, s) | refuted for this model (general N); special N gives only SNFS, L[1/3, 1.526] |
| 2 | Aurifeuillian factorisations | Splits Φ_n(m) into 2 pieces for special m | Factor 2 in log-size | Special N | constant factor only (snfsrig §3) |
| 3 | Consecutive products / binomials (Kummer, Legendre) | "Smooth by construction": C(n, k) is n-smooth | Abundant (O(B²) values ≤ 2^B). But products of factor-base primes are already free: supply was never the bottleneck | None. Exact integer identities carry no information mod N (codex Barrier 1). Codex test: 10 collisions vs 8.64 expected | refuted (no excess yield) |
| 4 | Factorials mod N (Pollard–Strassen) | Lemma CF: ∏_{c<p}(z − c) ≡ σ(z) − z mod p is the systematic equation over F_p | N/A | Needs degree ≥ p | exponential, N^{1/4} (KNOWN, BGS 2007). The exponent is forced by Lemma CF. A polylog SLP for n! would factor N (Lipton 1994; Shub–Smale τ) |
| 5 | Chebyshev/Dickson composition D_m ∘ D_n = D_{mn} | A cheap iterate (E2-like, but multiplicative only). It is cyclotomic in disguise (x = u + 1/u) | As in row 1 | Mod N: the p + 1 (torus) group-order method | group-order method, not relations |
| 6 | Elliptic divisibility sequences / division polynomials ([m] as "Frobenius") | E1 over Q(E[n]) (degree ~n⁴) | Primitive part ≥ (6/π²)·log\|W_n\|: gain ≤ π²/6. Primitive divisors ℓ ≥ (√n − 1)² (Silverman) | Mod N: ECM, L_p[1/2] | group-order method |
| 7 | Lucas sequences | Cyclotomic over a quadratic field | As in row 1 (Bilu–Hanrot–Voutier primitive divisors, n > 30) | p ± 1 methods | group-order method |
| 8 | Cyclotomic norm forms, Gauss/Jacobi sums, Stickelberger | Galois automorphisms σ_a of (Z/N)[ζ] exist over every N (E2 as an automorphism, not as a power). Jacobi sums give explicit smooth norms ℓ^{d/2} | Explicit generators of small ideals | Fixed ring ⊇ Z/N, so Lemma CF applies. A map to Z/N needs a root of Φ_m mod N (special N) | NFS in cyclotomic fields: covered by RES/G. Jacobi sums mod N are the APR-CL *primality* test (a Frobenius exists iff N is prime) |
| 9 | p-adic / Witt lift of X^q − X | Teichmüller: X^p − X = ∏(X − [a]) in Z_p[X] (a 2-term polynomial with p roots in Q_p) | — | A unital map Z/p^s → Z/N exists iff N \| p^s (codex). With ℓ small and known, the ℓ-adic constants never reach Z/N | refuted (r20 #6; heurnfs §5(ii)) |
| 10 | Group rings, x^{λ(N)+1} ≡ x (Carmichael) | The *genuine* field equation of Z/N (squarefree N): T^{λ+1} − T vanishes on Z/N | — | Degree λ(N) ≥ p − 1. Knowing it factors N (Miller) | equivalent to the problem |
| 11 | Λ-rings / Adams operations / Borger's F_1 | ψ_q(T) = T^q is a Frobenius lift on Z[T] by substitution | — | Lemma SUB: σ_q(f)g − fσ_q(g) = det·(T^q − T) for linear f, g. The homography trick yields **one** relation, not q³ [PROVED, tested] | refuted |
| 12 | Integer systematic equations for q ≤ 3: T² − T = T(T − 1), T³ − T = (T − 1)T(T + 1) | E1 holds exactly over Z for q = 2, 3 only (tested q ≤ 8) | 4 pieces | No additivity: (aT + b)³ is not linear in T³ | no compression (E4 fails) |
| 13 | **Best: the "F_{1^n}" V-family on Cunningham N = Φ_n(m)** (§3) | E1 (μ_n = ⟨m⟩ ⊂ Z/N explicit), E3 (m^n ≡ 1, the analogue of X^q ≡ h₀/h₁), E4 *small side* (f(m) ≡ (A + C)m + (B + D)) | E4 *split side* fails (Lemma DES). f(m) is a random integer of size ≥ N | Special N only | measured G = 1.03–1.06; refuted |

## 3. The best candidate, developed (Step 3)

**Construction (the closest instance of E1 + E3 + half of E4).**
- Let N | Φ_n(m). Then m has order n mod N, μ_n = ⟨m⟩ ⊂ (Z/N)^× is an explicit "constant field F_{1^n}", and m^n ≡ 1 is the analogue of X^q ≡ h₀/h₁.
- For f = A x^{n+1} + B x^n + C x + D ∈ V_n:

  f(m) ≡ r_f := (A + C)m + (B + D) (mod N).

  The small side has size ≈ 2Hm: the compression half of E4 holds.
- A BGJT relation needs f to split into linear factors with small values.

**What fails: Lemma DES.**
- For n ≥ 6, the only members of V_n(R) that split over R are x^n(Ax + B). Their relation reads m^n(Am + B) ≡ Am + B, i.e. just m^n ≡ 1.
- So f(m) is an unstructured integer of size H·m^{n+1} ≥ H·m^{n+1−φ(n)}·N, i.e. log-size ≥ log N.

**Heuristic cost [HEURISTIC: CEP smoothness of linear-form values].**
- Relations need f(m) B-smooth at size ≥ N. This is the Dixon/random-squares regime, with optimum L_N[1/2, √2], worse than SNFS.
- The same structure used as a number field sieve (algebraic side t − sζ, rational side t − sm) is SNFS, L[1/3, 1.526]. There the "systematic product" ∏_j(t − sζ^j) appears as the irreducible norm form Φ_n(t, s), whose values are generic.
- For general N, a pair (m, n) with m^n ≡ small (mod N) at small height violates heuristic G.

**Experiment E2 (pre-registered).** N of 37–43 bits, B = 2^16, H = 32, 200 000 candidates each. Baseline: the exact empirical smooth rate of 200 000 uniform integers per bit length.

| (n, m) | N bits | observed | predicted (size bound) | G |
|---|---|---|---|---|
| (5, 520) | 37 | 2325 | 2221.6 | 1.047 |
| (7, 130) | 43 | 1684 | 1628.1 | 1.034 |
| (11, 19) | 43 | 3820 | 3674.0 | 1.040 |
| (13, 11) | 42 | 5574 | 5263.9 | 1.059 |

- (P2a) G ∈ [0.5, 2]: **held**.
- (P2b) Kill signal G ≥ 4: not triggered.
- The consistent +4–6% is HEURISTIC: candidate values are log-uniform inside each dyadic bit-length bin, while the baseline is uniform, so candidates sit lower in the bin and are slightly smoother. This was not tested further.

**Positive control E1 (pre-registered; exact enumeration of P(V)(F_{q²})).**
- G_FF = (relations/|P(V)|) / (P_rand(q+1)·P_rand(3)), where P_rand(D) is the exact fraction of split binary forms of degree D.
- Every relation was verified as h₁f ≡ r_f mod I.

| q | k = deg I | \|P(V)\| | split in V | relations | rel/(Q+1) | G_FF | G_FF* (prereg model) |
|---|---|---|---|---|---|---|---|
| 4 | 6 | 4369 | 357 | 160 | 9.4 | 9.1 | (not predicted) |
| 5 | 6 | 16276 | 806 | 362 | 13.9 | 38 | (not predicted) |
| 7 | 8 | 120100 | 2850 | 1362 | 27.2 | 1.27e3 | (not predicted) |
| 8 | 8 | 266305 | 4745 | 2215 | 34.1 | 8.7e3 | 2.7e3 |
| 9 | 9 | 538084 | 7462 | 3464 | 42.2 | 7.0e4 | 2.0e4 |
| 11 | 13 | 1786324 | 16226 | 7626 | 62.5 | 6.4e6 | 1.5e6 |
| 13 | 14 | 4855540 | 31110 | 14996 | 88.2 | 8.6e8 | 1.7e8 |
| 16 | 15 | 16843009 | 70161 | 33776 | 131.4 | 2.3e12 | 3.8e11 |

- (P1a) Within one decade: **held** (deviations 0.50–0.79 decades). The model's 8/q³ underestimates frac_split ≈ (q + 1)/q³.
- (P1b) Monotone: **held**.
- (P1c) Relations exceed the factor base: **held**.
- (P1d) All relations verified: **held**.

**Post-hoc decomposition (disclosed, not pre-registered; `systeq_E1_posthoc.jsonl`).**
- The PGL₂ orbit (q³ + q forms) yields relations at 0.10–0.18 per form, ≈ 1/6. G_orbit runs from 0.4 (q = 4) to 4.3·10¹⁰ (q = 16).
- The Frobenius-power family (X + b)^q(X + c) yields ≈ 1/2 (the relation holds when the quadratic h₀ + b^q h₁ splits). G reaches 2.3·10¹² at q = 16.
- That family is the exact counterpart of the integer degenerate x^n(Ax + B). In F_q it carries the nontrivial Frobenius relation (X + b)^q ≡ (h₀ + b^q h₁)/h₁. In Z/N its counterpart says only m^n ≡ 1.

**E3 (pre-registered; exhaustive, H = 4).** Nondegenerate members of V_n(Z) splitting over Q:

| n | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|
| count | 208 | 16 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

- (P3: 0 for n ≥ 7) **held**. Lemma DES in fact proves 0 for n ≥ 6.

**E4 (pre-registered).**
- n = 6, m ∈ [2^20, 2^20 + 30000): R_piece = 113/184.8 = 0.61.
- n = 12, m ∈ [4096, 16384): R_piece = 8/31.9 = 0.25.
- (P4: [0.1, 1.5]) **held**. Cyclotomic pieces are *less* smooth than random integers of their size, because their primes are confined to ≡ 1 (mod d).

## 4. Why "smooth by construction" is not the bottleneck in Z

Smooth integers by construction are free: any product of factor-base primes, binomials, factorials. What BGJT supplies is a *congruence* between a free-smooth object and a *small* object, and that congruence comes from E2 + E3.

- Codex Barrier 1 [PROVED]: two products of small integers that are both below N and congruent mod N are equal.
- So the integer question is entirely "where does the small side come from?". By Lemma RES / heuristic G (heurnfs B1), every polynomial presentation pays N^{2/D}E^D.
- In F_q, E3 + E4 make the small side's degree independent of the target degree k. Integers have no presentation with that property for general N (§5.5).

## 5. Barriers (Step 4)

### 5.1 Prop F (additive power maps exist only where factoring is already done) [PROVED]

**Statement.** Let A ≠ 0 be a commutative Z/N-algebra and s ≥ 2 with (x + y)^s = x^s + y^s for all x, y ∈ A. Then (p − 1) | (s − 1) for every prime p | char A. If A is faithful (char A = N) and N is squarefree, then λ(N) | s − 1. Given s, N then factors in randomised time poly(log N, log s) (Miller 1976; any multiple of λ(N) suffices).

**Proof.**
- π(z) = z^s is additive and multiplicative with π(1) = 1, so it is a ring endomorphism, and π(𝔪) ⊆ 𝔪 for every ideal 𝔪.
- For p | char A, A/pA ≠ 0 (CRT on Z/N), so A has a maximal ideal 𝔪 of residue characteristic p. π induces z ↦ z^s on κ = A/𝔪, which is a ring endomorphism.
- On the prime field F_p ⊆ κ, a ring endomorphism is the identity. So c^s = c for all c ∈ F_p, and with c a generator of F_p^× this gives (p − 1) | (s − 1).
- If char A = N is squarefree, lcm_{p|N}(p − 1) = λ(N). ∎

**Consequences.**
- s ≥ p for each p | N. The BGJT exponent s = |k| therefore has size ≥ the hidden prime.
- Tests: `test_prop_F_additive_exponents`. On (Z/15)[t]/(t² + 1) ≅ F₉ × F₅ × F₅ the additive s ≤ 40 are exactly those predicted; on Z/21 they are s ≡ 1 mod 6.
- Codex adds that on Z/N itself the additive power maps are the identity map (e.g. u ↦ u^{561} on Z/561), and that a nonunital endomorphism is z ↦ ez with e idempotent, which already yields a factor.

### 5.2 Lemma CF (the constant field is the hidden prime field) [PROVED]

**Statement.** Let A be a Z/N-algebra, S ⊂ A finite and P(T) = ∏_{c∈S}(T − c). If P vanishes on the image of Z in A, then for each p | char A the image of S in a residue field of characteristic p contains F_p, so |S| ≥ p. This holds in particular when S is the fixed ring of a ring endomorphism, since every fixed ring contains Z·1.

**Proof.** In κ[T] the product P̄ = ∏(T − c̄) has every c ∈ F_p as a root, so c = c̄ for some c ∈ S. ∎

**Instance (tested).** In F_p[√D] with D a non-residue, ∏_{c=0}^{p−1}(z − c) = σ(z) − z, where σ is conjugation, which equals Frobenius. This is the Wilson/factorial shadow of X^p − X. Pollard–Strassen and Bostan–Gaudry–Schost evaluate it in ≈ √p, which gives N^{1/4}. That is the precise sense in which "factorial methods are exponential": the only computable systematic equation over Z/N has the hidden prime field as its constant field.
- KNOWN context: Lipton 1994 and Shub–Smale 1995 show that short straight-line programs for n! would factor fast. This is the integer counterpart of the function-field fact that X^{q^k} − X has an O(k log q) SLP, which is distinct-degree factorisation.

### 5.3 Lemma DES (fewnomials have few roots in characteristic 0) [PROVED from KNOWN]

**Statement.**
- Descartes' rule: f ∈ R[x] with k nonzero terms has at most k − 1 positive and k − 1 negative roots, with multiplicity.
- Hence a k-term f that splits over R has deg f ≤ 2(k − 1) + ord₀ f.
- For V_n with A ≠ 0:
  - if D ≠ 0, splitting forces n + 1 ≤ 6;
  - if D = 0 ≠ C, then f = x(Ax^n + Bx^{n−1} + C) forces n ≤ 4;
  - if C = D = 0, f is the degenerate x^n(Ax + B).
- So for n ≥ 6 there are **no nondegenerate split members**.
- Over a number field K, Lenstra (1999, "On the factorization of lacunary polynomials") bounds the number of roots in K of a k-nomial, apart from 0, by c(k, [K : Q]), independent of the degree.

**Consequence.**
- In F_q the 2-nomial X^q − X has q roots, which is E1 + E2. No number field of bounded degree has an analogue.
- Unbounded-degree fields (Q(ζ_n), Q(E[n])) give splitting only at the price of the Splitting Lemma (§5.4).
- Over Q_p the 2-nomial X^p − X does have p roots (Teichmüller). That is the only characteristic-0 escape, and it lives at the hidden prime (row 9).

### 5.4 Splitting Lemma and anti-smoothness [PROVED / KNOWN]

**Splitting Lemma.** For F ∈ Z[x], t ∈ Z and any number field K ∋ α with F(α) = 0:

N_{K/Q}(t − α) = ±(minpoly_α(t))^{[K:Q(α)]}, up to leading-coefficient factors.

So any factorisation of F over any K induces, on the integer F(t), nothing finer than the Q-irreducible factorisation F(t) = ∏ F_j(t). Hence:
- **Cyclotomic:** the largest piece of t^n − 1 is Φ_n(t), with log-size φ(n) log t. The gain in log-size is ≤ n/φ(n) ≤ e^γ log log n + O(1/log log n) (Rosser–Schoenfeld 1962, Thm 15). This changes L[α, c] by at most a factor (log log)^{O(1)} inside the constant, never the exponent.
- **EDS:** the primitive part has log-size ≥ (6/π²)·log |W_n|, so the gain is ≤ π²/6.

**Anti-smoothness [KNOWN ⇒ PROVED corollary].**
- Zsigmondy (1892): m^n − 1 (m ≥ 2) has a primitive prime ℓ ≡ 1 (mod n), apart from (2, 6) and n = 2.
- So m^n − 1 is not B-smooth whenever n ≥ B. Test `test_zsigmondy_antismooth`, and E4 (R = 0.25 at n = 12).
- Lucas sequences: Bilu–Hanrot–Voutier 2001, n > 30. EDS: Silverman 1988, with ℓ ≥ (√n − 1)² by Hasse.

### 5.5 Rigidity of Z/N (no presentation freedom) [PROVED part + HEURISTIC part]

**F_q side [KNOWN].** BGJT/GKZ need some irreducible I of degree k dividing h₁X^q − h₀. Any such I presents *the* field F_{Q^k}. KW prove that an elliptic model exists after a degree-O(log n) extension. In both cases the freedom comes from all targets of one size being isomorphic.

**Z/N side.**
- [PROVED] Every ring presentation Z[x]/(f) → Z/N with x ↦ m satisfies Lemma RES (snfsrig §1.1).
- [HEURISTIC G] An integer counterpart of X^q ≡ h₀/h₁, namely N | h₁(m)m^q − h₀(m) with small h_i and small m, has ≈ (Hm²)²·M parameter choices against probability 1/N. For general N this needs Hm² ≥ N^{1/3}-ish.
- [PROVED] Even when such a presentation exists (special N), Prop F and Lemma DES leave the f-side unstructured. The best use of the structure is SNFS.
- [PROVED] Extensions of the target ((Z/N)[x]/(F)) do not change the prime field, so Lemma CF still applies.

### 5.6 E5 has no integer counterpart [PROVED for the size-invariance part]

**Quasi-polynomial descent.** It halves the absolute degree at cost poly(q), independent of the constant field F_{q^{2^e}} in which it runs: the splitting probability in V is 1/O(q³) over every extension. The size of the other side is (1 + δ)D, with no dependence on the target degree k.

**Integer side.**
- Under RES/G every candidate tied to N has log-norm-product ≥ 2√(2·log N·log E) (heurnfs B1). The other side's size therefore depends on N, and the descent's bottom level costs L[1/3] (heurnfs B5, §3.1 bottom bound).
- [PROVED] In a number field a rational prime ℓ has prime ideals of norm ℓ^f ≥ ℓ, so base change never lowers absolute size. There is no "extend the constants" step.

### 5.7 The F_1 dictionary [PROVED, elementary]

- Over "F_{1^n}" = μ_n ∪ {0}, GL₂ is the monomial group μ_n² ⋊ S₂ (Kapranov–Smirnov convention), so PGL₂(F_{1^n}) has order 2n.
- It sends T^n − 1 only to unit multiples of T^n − 1 and 1 − T^n. That is O(1) relations, against q³ + q for PGL₂(F_{q²})/PGL₂(F_q).
- Additivity is exactly what turns the monomial group into PGL₂. Its absence (Prop F) collapses E4.

## 6. What it would take (precise missing step)

A route-1 success would need a family 𝓕 of pairs (u, v) with u ≡ v (mod N), u ≠ v, satisfying:
- (a) u is a product of ≳ log N/log B pieces, each ≤ B^{O(1)}, by an identity;
- (b) v has log-size independent of log N (or o(x^{2/3}));
- (c) |𝓕| ≥ L[α] for α < 1/3, computable without the factors.

By Prop F, Lemma CF and Lemma DES, (a) cannot come from a Frobenius, from a systematic product over Z/N, or from a fewnomial identity over a number field of bounded degree. By RES/G, (b) cannot come from polynomial selection for general N. Neither we nor codex found any other source.

## 7. Reviews

### 7.1 Codex blind round (`codex/systeq_1.out`, launched before our analysis)

**Same verdict.** It gives < 1% for a straightforward transfer.

**Independent contributions:**
- the E1–E6 axiomatisation;
- the correct coset count q(q² + 1);
- unital Z/p^s → Z/N exists iff N | p^s;
- the Λ-ring equivariance failure, i.e. a^p ≡ a mod N has density ≤ p²/N;
- the explicit carry T^p − T = ∏(T − a) + pC_p(T);
- the Jacobi-sum smooth norms ℓ^{d/2};
- the a⁴ + 4b⁴ Aurifeuillian example;
- the binomial collision kill test: N = 1 000 003 · 1 000 033, B = 4096, M = 4 155 871 binomials of ≥ 80 bits, **10 collisions observed against 8.635 expected for uniform units** (random control: 8). The script is `codex/scratch/systeq_1/` (git-ignored; not re-run by us).

**What we add.** Prop F (in its quantitative "s − 1 ≡ 0 mod λ(N)" form), Lemma CF, Lemma DES applied to V, the measured positive-control contrast and the Cunningham V-family. For each of these: no prior statement found by us.

### 7.2 Codex adversarial round and Opus critics

See §7.3 (filled in after the reviews).

## 8. Reproduction

```bash
.venv/bin/python -m src.systeq.run_experiment all      # E1-E4, ~3 min
.venv/bin/python -m pytest -q tests/test_systeq.py
```

The post-hoc decomposition is `bgjt.posthoc_by_type(p, e0)` over the E1 fields. It is also `run_experiment.py E1_posthoc` (same output file).
