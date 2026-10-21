# schnorrcvp: factoring via CVP/SVP in prime-number lattices (Schnorr 1991/2021, Adleman 1995)

Agent `schnorrcvp`, phase 19, 2026-10-21. Version 2: after codex adversarial review and two Opus critics (§9).

- Code: `src/schnorrcvp/`:
  - `core.py`: exact relation objects, Λ_N, GF(2) square roots;
  - `psi.py`: exact Ψ brackets and the counting models;
  - `lattice.py`: fpylll lattices and list-CVP;
  - `enum_rel.py`: exact box enumeration;
  - `run_e0*.py`, `run_e1.py`, `run_e2.py`, `check_list0.py`, `analyze.py`.
- Tests: `tests/test_schnorrcvp.py`.
- Pre-registration: `notes/schnorrcvp_prereg.md`, commit `ecd3c25`, made before E1/E2. Raw data:
  `experiments/raw/schnorrcvp_E*.jsonl`.
- GPT-6.1-sol: blind round `codex/schnorrcvp_1.{md,out}`; adversarial round `codex/schnorrcvp_2.{md,out}`.

Labels:
- **PROVED**: proof here, or a theorem checked in the primary source.
- **KNOWN**: published, with citation.
- **HEURISTIC (name)**: rests on a named assumption.
- **MODEL**: an exact consequence of named heuristics, computed in code.
- **SCENARIO**: a cost model with stipulated algorithmic assumptions; not a lower bound.
- **MEASURED**.
- **POST-HOC**: chosen after seeing data.
- **CONJECTURED**.
- **OPEN**.

"No prior statement found by us" is the strongest novelty claim made here. It reports a search, not a proof of novelty.

---

## 0. Verdict

**Plain language.**

- **Schnorr 2021 is wrong on its own terms.**
  - His key Lemma 3.1 drops a factor of about N from the residual u − vN. Daniel S pointed out the N vs N′ confusion
    publicly in March 2021.
  - His 800-bit norm estimate is impossible.
- **Worse, at his parameters the relations he needs almost certainly do not exist at all.** The total expected
  number of fac-relations of every size is about 10^{−52} for 400 bits with n = 47, and 10^{−115} for 800 bits with
  n = 95. He needs 48 and 96. (Heuristic: smooth-number densities, validated at toy sizes.)
- **Required dimension.** Our model's crossover, where enough relations exist, is n ≈ 275 at 400 bits and ≈ 775 at
  800 bits. It grows roughly like (ln N)^{3/2}.
- **Prior work.** No published count of this kind was found by us or the literature critic.
- **Factoring does reduce to approximate CVP in these lattices, but heuristically and only for γ → 1:**
  - γ = 1 + n^{−ε} gives polynomial time;
  - γ = 1 + (log n)^{−c} gives L[1/(2+c)];
  - γ = 1 + η gives L[1/2, √(γ²−1)], which beats QS only for γ < √2.

  These figures assume a fast oracle on this lattice family. The ℓ1 version of this reduction is essentially in
  Schnorr 1991; we redo it in ℓ2 with explicit constants and rigorous distinctness.
- **At large factors the certificate gives nothing.** For the approximation factors lattice algorithms achieve and
  lattice cryptography uses (poly(n), 2^{n^ε}, 2^{n/log n}), the guaranteed residual bound becomes useless. The
  residual's size is exponentially sensitive to γ.
- **What this route cannot show.**
  - Whether a cleverer reduction could work at large γ is open. Our tightness argument is only a conjecture.
  - Rigorously, factoring already reduces to GapCVP at much larger factors (n^{c/log log n}) on general lattices,
    simply because factoring is in NP (folklore).
  - **No link** between factoring and the lattice problems of post-quantum cryptography is known in either
    direction, and this route does not supply one.
- **Cost with real lattice algorithms (SCENARIO).**
  - CVP at γ = 1 + O(1) costs 2^{Θ(n)} with every known algorithm, and the dimension must exceed n_min. So the method
    costs exp((ln N)^{3/2+o(1)}) in this scenario, asymptotically slower than trial division.
  - Concretely: 2^{187} at 400 bits and 2^{467} at 800 bits.

**Precisely.**

1. **Theorem R [PROVED].** For the ℓ2 lattice with diagonal √ℓ_i (Adleman's weights; ℓ_i = ln p_i):
   - ‖b(e) − t‖² = Q(e) + C²Δ²;
   - w = u − vN = 2e^{(A+L)/2} sinh(Δ/2);
   - A = ln(uv) ≤ Q;
   - gcd(uv, w) = 1;
   - |uw| ≥ N − 1. This last bound is resbarrier's.

   Schnorr 1991's ℓ1 lattice (weights ℓ_i) has A = ‖·‖₁ exactly. Vera 2010 Thms 1–2 give ℓ_p versions.
2. **Theorem P (promise reduction) [PROVED in the stated parameter regime].**
   - If a squarefree witness relation with |w₀| ≤ ω exists, every γ-CVP answer at scale C = e^{(L+A₀)/2}/ω has
     |w| ≤ ω·γ√(A₀+1)·exp(((γ²−1)A₀+γ²)/2 + R/(2C)).
   - Corollary P1: an exact oracle with a scan over C returns relations.
   - This is the ℓ2/approximate form of Schnorr 1991 Thm 2 and Adleman 1995 §2.
3. **Lemma LA [PROVED].** Testing a GF(2)-kernel *basis* suffices. It factors N deterministically if (a) the relations
   generate an odd-index sublattice of Λ_N (a strong sufficient condition) and (b) ⟨P⟩ = (Z/N)^*. (b) holds under
   ERH for B ≥ 2 ln²N when N has no prime factor ≤ B (Bach 1990 Thm 2).
4. **Theorem H (heuristic reduction) [HEURISTIC H_S per target; H_LA; ERH; N not a prime power].**
   - Setting: fixed κ > 2, B = (ln N)^κ, n = π(B), targets kN for primes k ∈ (B, 4B].
   - Then n + 2 pairwise-distinct relations come from O(n·poly L) calls to an exact or γ-CVP oracle, with
     γ − 1 ≤ c_κ ln B/L.
   - Distinctness is PROVED.
   - So factoring ≤ CVP_{1+n^{−ε}} on these lattices, for any fixed κ > max(2, 1/ε).
   - The witness count is KNOWN (Schnorr 1991 Thm 5, under his "Hypothesis").
5. **Theorem G (cost vs γ) [MODEL, benign oracle, unit-cost calls].**
   - γ = 1 + τL^{2α−1}(ln L)^{1−2α} gives L[α, 2√(ατ)].
   - With oracle cost n^c this becomes L[α, 2√((1+c)ατ)].
   - **Conjecture A** (tightness for adversarial oracles; downgraded from "Prop A" after codex found a gap) would make
     the regime of item 4 the exact one for worst-case oracles.
6. **Prop NP [PROVED from known theorems; folklore].**
   - Factoring ≤_T GapCVP_γ for γ ≤ n^{c/lglg n} (DKRS 2003).
   - Factoring ≤_R GapSVP_γ for every constant γ (Khot 2005; Haviv–Regev 2007).
   - So Adleman's 1995 goal was superseded by Ajtai 1998 (SchnorrGate README).
   - An unrefereed 2026 manuscript claims deterministic NP-hardness of GapCVP_{n^{1/400}}.
   - Unknown window: (n^{c/lglg n}, √(n/log n)), or (n^{1/400}, √(n/log n)) if that claim is correct.
7. **Schnorr 2021 [PROVED errors; MODEL counts].**
   - Lemma 3.1 is missing a factor N/N′ = N^{n/(n+1)}. KNOWN since March 2021; codex and we make the factor
     explicit.
   - "‖b₁‖ < 0.917" is impossible.
   - (3.3) is inconsistent with the diagonal f(i).
   - E_tot = 10^{−51.9} at (400, 47) and 10^{−114.8} at (800, 95). His own 800-bit box holds 10^{−176}.
8. **Experiments (pre-registered; every pre-registered prediction passed).**
   - E1, exact box counts: H_S′ ratio 1.63 pooled. The post-hoc parity/CRT model gives 0.92.
   - E2, 78,963 CVP outputs:
     - pre-registered H_res′ ratio 0.87;
     - after codex found probabilities > 1, the normalised predictor gives 1.29, with the residual excess due to
       parity;
     - Schnorr-2021 scaling leaves |w| ≈ N^{0.92–0.96};
     - all N of 20–36 bits were factored end-to-end (correctness check).

**Best open question.** Is there a lattice family, norm or embedding in which a *non-NP-hard* approximation factor
(γ ≥ √(n/log n)) controls the residual of a factoring relation? Theorem P's exponential sensitivity is what makes the
Schnorr/Adleman strategy useless there. A positive answer needs a residual whose log-size is not linear in the
lattice norm.

---

## 1. Setting and exact identities

Notation:
- N is odd, not a prime power, with no prime factor ≤ B (trial-divide first).
- P = {p₁, …, p_n} are the primes ≤ B; ℓ_i = ln p_i; L = ln N.
- e ∈ Z^n gives coprime u = ∏_{e_i>0}p_i^{e_i}, v = ∏_{e_i<0}p_i^{−e_i}, and w = u − vN.
- A *fac-relation* has |w| P-smooth. It gives u ≡ w (mod N).

**The lattices** (literature critic):
- **Schnorr 1991**: diagonal ℓ_i, last row N^c ℓ_i, target (0, …, 0, N^c L), **ℓ1 norm**. B = (ln N)^α, c > 1.
  Rigorous Thm 2: ‖z − N‖₁ ≤ (2c−1)L + 2δ ln B ⇒ |u − vN| ≤ B^{1/α+δ+o(1)}.
- **Adleman 1995** (draft of 16 March 1995, SchnorrGate `LATFAC.pdf`, read in full): diagonal √ℓ_j, ℓ2, exact SVP,
  first column M·ln p_j with M = N⁴, row V₀ = (M(ln N + γ), 0, …), B = ln^{20}(2N^{1.25}).
- **Vera 2010** (arXiv 1003.5461): ℓ_p versions with weights ℓ_i^{1/p}. Remark 9 asks for an ℓ2 analysis.
- **Schnorr 2021**: diagonal f(i), a random permutation, C = N′ = N^{1/(n+1)}, SVP on the (n+1)-dimensional
  embedding.
- **We** use Adleman's ℓ2 weights in CVP form: b(e) = (√ℓ₁e₁, …, √ℓ_ne_n, C⟨e,ℓ⟩), t = (0, …, 0, CL).
  Q(e) = Σℓ_ie_i², A(e) = Σ|e_i|ℓ_i = ln(uv), Δ = ⟨e,ℓ⟩ − L.

**Theorem R [PROVED].**
- (i) ‖b(e) − t‖² = Q + C²Δ².
- (ii) w = vN(e^Δ − 1) = 2e^{(A+L)/2} sinh(Δ/2).
- (iii) A ≤ Q, with equality iff uv is squarefree.
- (iv) gcd(uv, w) = 1.
- (v) |uw| ≥ N − 1 whenever u − vN = w with u, v, w ≠ 0.

*Proof.*
- (ii): vN = e^{(A+L−Δ)/2}.
- (iii): |e_i| ≤ e_i².
- (iv): p | u, w ⇒ p | vN; p | v, w ⇒ p | u.
- (v): N ≤ |v|N ≤ |u| + |w| ≤ |uw| + 1. ∎

Tests: `test_residual_identity_exact`, `test_uw_lower_bound`, `test_cvp_output_obeys_certificate`.

**Relation lattice.**
- Λ_N = {(a, s) : (−1)^s∏p^a ≡ 1 (mod N)}, with relation vectors exps(u) − exps(|w|).
- The seminorm Σ|a_i|ℓ_i vanishes on (0, 2). So work in the projection π(Λ_N) = {a : ∏p^a ≡ ±1}, of index ≤ φ(N).
- There (v) says every nonzero vector has Σ|a_i|ℓ_i = ln(u|w|) ≥ ln(N − 1).
- Minkowski in R^n gives λ₁(π(Λ_N)) ≤ (n!∏ℓ_i φ(N))^{1/n} ≈ L(ln L + lnln L − 1 + o(1)) for n ≈ L [PROVED; the
  critic checked it numerically].
- So relations u ≡ ±w with u|w| ≤ N^{(1+o(1))lnln N} always exist. But v = (u ∓ w)/N is arbitrary, not smooth, so
  they are relations of Adleman's multiplier form, not vectors of Schnorr's lattice.

**Gaussian heuristic warning [HEURISTIC].**
- For n ≫ L the short vectors are mostly in {0, ±1}^n.
- Ball counts follow smooth-pair counts, which grow like e^{T(1−1/κ)} in ln(uv) ≤ T, not like the volume T^n.
- Ducas (RISC slides 2021-05-20, slide 18): "Gaussian Heuristic seems invalid for certain C".

## 2. Reduction theorems

### 2.1 Theorem P (promise reduction) [PROVED]

*Witness.* e⁰ ∈ {0, ±1}^n with A₀ = ln(u₀v₀) and 0 < |w₀| ≤ ω. Set C = e^{(L+A₀)/2}/ω.

*Claim.* For every γ ≥ 1 and every b(e) with ‖b(e) − t‖ ≤ γ·dist(t, Λ_C):

  |w(e)| ≤ ω·R·exp((R² − A₀)/2 + R/(2C)),   R = γ√(A₀+1),  R² − A₀ = (γ²−1)A₀ + γ².

*Proof.*
- |w₀| ≥ e^{(A₀+L)/2}|Δ₀|, so C|Δ₀| ≤ 1 and dist² ≤ A₀ + 1.
- For the answer, A ≤ Q ≤ R² and |Δ| ≤ R/C.
- sinh y ≤ y·e^y for y ≥ 0. ∎

Every step was checked by the proof critic and by codex.

**Corollary P1 (exact oracle, A₀ unknown).**
- Scan j ∈ [⌊L⌋ − 1, A_max], with A_max = poly(L) part of the promise.
- For each j, put ω_j = B·e^{−0.7}/√(j+2) and C_j = e^{(L+j)/2}/ω_j.
- At j = ⌊A₀⌋ (δ = A₀ − j), C_j²Δ₀² ≤ e^{−δ}, and every exact answer has |w| ≤ ω_j√(A₀+1)e^{(δ+e^{−δ})/2}
  ≤ 0.99·B. So it is a fac-relation, since every integer ≤ B is B-smooth, and w ≠ 0 since N has no factor ≤ B.
- The promise is a witness with |w₀| ≤ ω_{⌊A₀⌋}.
- These are the critic's constants. Our draft's e^{−0.51} with "extra e^{1/2}" did not close.

**Rounding [PROVED in the regime A₀, A_max = O(L), γ = O(1)].**
- Use integer entries round(S√ℓ_i), round(SCℓ_i), round(SCL).
- Q changes by a factor 1 + O(1/S). The last coordinate changes by at most (‖e‖₁+1)/(2S) ≤ (Q/ln 2 + 1)/(2S) in
  units of CΔ.
- S ≥ L² suffices in this regime. For general heights or γ take S ≫ R³ (codex).
- The logarithms must be certified to that precision. `lattice.py` uses 80-digit mpmath, which covers our toy sizes
  and is not a general certified construction.

### 2.2 Lemma LA [PROVED]

- Relations r₁, …, r_m ∈ Λ_N generate Λ′. K is the GF(2)-kernel of the parity map.
- c ↦ x(c) = (−1)^{y_s}∏p^{y(c)} with y(c) = ½Σc_ir_i is a homomorphism K → μ₂(Z/N). Carries change y by elements
  of Λ′.
- So a basis suffices to find any x ≠ ±1 in the image.
- If [Λ_N : Λ′] = d is odd, the image is all of μ₂ ∩ ⟨−1, P⟩. If x² = 1 with x = χ(z), then 2dz ∈ Λ′, and the
  corresponding c gives x^d = x.
- Under ERH, Bach 1990 Thm 2 (primary text, via the literature critic) says a proper subgroup containing all
  n < x (all units) has x < 2 ln²N. So with N free of prime factors ≤ B and B ≥ 2 ln²N, ⟨P⟩ = (Z/N)^*, which has
  ≥ 4 square roots of 1.
- Without trial division, use 3 ln²N (Bach Thm 3).
- **Caveat (codex).** Odd index is a strong *sufficient* condition, which we call H_LA. Congruence-of-squares methods
  can succeed without it, and ERH does not imply it.
- **Code fix.** `squares_dependency_factor` truncated the kernel basis at 64 dependencies by default. Codex
  constructed a counterexample (N = 35). The default now tests all, and a regression test was added. The E2-P6 runs
  used the old default and succeeded anyway.

### 2.3 Distinctness [PROVED] and Adleman's distinctness step [PROVED gap; no prior statement found by us]

- **Prime multipliers.**
  - Targets ln(kN) for primes k ∈ (B, 4B] give relations u − vkN = w ≡ u (mod N).
  - k ≠ k′ giving the same (u, w) would force vk = v′k′, which is impossible for P-smooth v, v′ and distinct primes
    k, k′ > B.
  - This proves distinct relations, not independent ones.
- **Adleman 1995 §2.**
  - To get different congruences he replaces "V₀ by V₀ + V_i + V_j and rerun[s]", intending to force p_ip_j | P.
  - Since V₀ is a basis vector, this is a unimodular change. The lattice, and hence the set of valid SVP answers, is
    unchanged (`test_adleman_rerun_lattice_unchanged`). So the literal argument has a gap.
  - A repair (delete V_i, V_j and use target ln(p_ip_jN)) is a *different construction*. It needs its own witness and
    shortest-vector argument. Its distinctness remains heuristic, because p_ip_j ∈ ⟨P⟩.
  - Searched without finding this observation: Vera 2010, Ducas (README, RISC and PKC slides), crypto.SE 88582,
    Mulder–Sterner–van Woerden.
  - Not checked (paywalled): Kan–Shen, "On Schnorr–Adleman lattice", PDCAT 2003.

### 2.4 Theorem H [HEURISTIC]

**Setting.** Fixed κ > 2, B = (ln N)^κ, n = π(B), targets kN for primes k ∈ (B, 4B] (≥ n + 2 of them).

**Hypotheses.**
- **H_S, for every target** (uniformity is part of the hypothesis): a squarefree witness exists with
  v ≤ N^a, A₀ ≤ (1 + 2a + o(1))L and |w₀| ≤ ω.
  - Expected witness count under H_S: N^{a(1−2/κ) − 1/κ − o(1)} → ∞ for a > 1/(κ−2).
  - This is Schnorr 1991 Thm 5's count N^{(c−1) − (2c−1)/α}, with his α = κ and c − 1 = a [KNOWN; literature critic].
- **H_LA.**
- **ERH.**

**Conclusion.**
- Corollary P1 gives n + 2 distinct relations from O(n·L) exact-CVP calls, and Lemma LA factors N.
- With a γ-oracle, Theorem P allows (γ² − 1)A₀ ≤ 2 ln B − 3 ln L − O(1), i.e.
  γ − 1 ≤ c_κ ln B/L with c_κ = (1 − 3/(2κ))/(1 + 2a) + o(1).
- Since ln B/L = n^{−1/κ+o(1)}:

> **factoring ≤ CVP_{1+n^{−ε}} on Schnorr–Adleman prime lattices of dimension n = (ln N)^κ, for any fixed
> κ > max(2, 1/ε), with O(n log N) oracle calls. Heuristic in H_S (per target) and H_LA; ERH; N not a prime power.**

**What is new versus Schnorr 1991 (literature critic).** His Thms 2 + 5 already give the ℓ1 analogue: an ℓ1-CVP with
slack 1 + O(ln B/L) yields relations under his Hypothesis. New here, as far as we found:
- the ℓ2 (Adleman-weight) certificate, which answers Vera's Remark 9;
- the explicit constant c_κ and approximation-factor form;
- the rigorous distinctness;
- the kernel-basis Lemma LA.

### 2.5 Theorem G (cost as a function of γ) [MODEL] and Conjecture A

**Model.**
- A *benign* γ-oracle returns a vector whose residual behaves like a random integer of the size allowed by Theorem P,
  ln|w| ≈ θL + ln B with θ = (γ²−1)(1+2a)/2 → (γ²−1)/2 as a → 0, under the forced local conditions (H_res′).
- With B = L[α, b], each call succeeds with probability exp(−u ln u), u = θL/ln B. We need n relations and n²
  linear algebra.

| oracle factor γ(n) | cost, unit-cost oracle (×√(1+c) on each constant for oracle cost n^c) |
|---|---|
| 1 + n^{−ε} | poly(log N) (Theorem H) |
| 1 + τL^{2α−1}(ln L)^{1−2α} ≈ 1 + (log n)^{2−1/α+o(1)} | L[α, 2√(ατ)] |
| 1 + (log n)^{−c} | L[1/(2+c)]; α < 1/3 ⇔ c > 1 |
| γ fixed > 1 | L[1/2, √(γ²−1)]; beats QS iff γ < √2; γ = 2 gives L[1/2, √3] |
| γ → ∞ with n: poly(n), 2^{n^ε}, 2^{n/log n} | the guaranteed bound becomes N^{Θ(γ²)}: the certificate gives nothing, and for γ = poly(n), n = L^{Ω(1)}, the benign model is worse than trial division |

- Derivation (§2.6). In the draft the γ ≥ 2 row read "no gain". Codex correctly noted that γ = 2 still gives
  L[1/2, √3] in this model.
- No oracle as fast as this table assumes is known. All known algorithms for γ = 1 + O(1) cost 2^{Θ(n)}.

**Conjecture A (tightness) [CONJECTURED].** If (γ² − 1)A₀ ≫ ln B, the γ-ball around t contains non-relations, so an
adversarial oracle defeats this query strategy.
- **Why "Prop A" was withdrawn (codex, FATAL to the claim).** The draft used radius γ√(A₀+1), but the true distance d
  may be much smaller than the witness bound. Also, A ≈ γ²A₀ does not give Q ≈ A without squarefreeness, and the
  §3 shell counts do not show *non-smooth* vectors inside the admissible ball jointly over the C-scan.
- A proof needs actual-distance bounds and a joint count inside the ball.
- Even if true, it would say nothing about other lattices or embeddings.

### 2.6 Derivation [MODEL; checked by the proof critic and codex]

- ln B = bL^α(ln L)^{1−α}, and u ln u = sL^α(ln L)^{1−α}. So u = (s/α)L^α(ln L)^{−α} and
  θ = (sb/α)L^{2α−1}(ln L)^{1−2α}.
- With unit-cost calls, the calls cost L[α, b + s] and linear algebra L[α, 2b]. The optimum s = b gives L[α, 2b], so
  θ = τ(…) gives L[α, 2√(ατ)].
- With cost n^c, the calls cost L[α, (1+c)b + s], giving 2√((1+c)ατ).
- α = 1/2: θ = (γ²−1)/2 gives L[1/2, √(γ²−1)].
- α = 1/3: beating GNFS 1.923 needs τ < 2.77/(1+c).

### 2.7 Prop NP [PROVED from known theorems; folklore]

- **CVP.** {(N, k) : N has a divisor in [2, k]} ∈ NP. GapCVP_γ is NP-hard for γ = n^{c/lglg n} under deterministic
  Karp reductions (Dinur–Kindler–Raz–Safra, Combinatorica 23 (2003)). Adaptive binary search therefore gives a
  deterministic polynomial-time *Turing* reduction of the factoring function to GapCVP_γ on lattices of dimension
  poly(log N). The queries satisfy the gap promise.
- **SVP.** GapSVP_γ is NP-hard under randomized reductions for every constant γ (Ajtai 1998; Micciancio 2001;
  Khot 2005; Haviv–Regev 2007/ToC 2012 Thm 3.1). So factoring ≤_R GapSVP_γ.
- **Folklore, not ours.** The SchnorrGate README says Ajtai "instead found a proof of a much stronger statement:
  namely that SVP is NP-hard", while trying to complete Adleman's goal. Micciancio (2012 talk) says the prime lattice
  "played a fundamental role in establishing the NP-hardness of SVP". Ducas's RISC slides p. 22 say the same.
  - Adleman's own framing is in §1 of his draft: "the first substantial evidence that Lattice is not in P".
- **Newer claims (not refereed, not checked by us).**
  - "Ten Advances in Mathematics and Theoretical Computer Science" (OpenAI, updated 6 Aug 2026), Ch. 7, claims a
    deterministic Karp reduction 3SAT → GapCVP^{(2)}_{n^{1/400}}. If correct, Prop NP extends to γ ≤ n^{1/400}.
  - Hair–Sahai (ePrint 2025/2181) claim deterministic SVP hardness in ℓ_p, p > 2. That does not cover ℓ2.
- **Upper side.** GapCVP_γ ∈ NP ∩ coAM for γ = O(√(n/log n)) (Goldreich–Goldwasser 2000), and ∈ NP ∩ coNP for
  γ = C√n (Aharonov–Regev 2005, Thm 1.1). LWE-based cryptography uses γ = Õ(n)–poly(n) (Regev, STOC 2005 /
  JACM 56(6) 2009).
- **Consequences.**
  - (i) Theorem H's γ = 1 + n^{−ε} lies inside the NP-hard range. Its content is the restriction to an explicit
    family, which is incomparable with the general-lattice reduction.
  - (ii) No reduction from factoring to GapCVP_γ is known for γ ≥ √(n/log n). Such a reduction would contradict
    nothing, since factoring ∈ NP ∩ coNP. This route's certificate cannot supply one; whether any route can is open.
  - (iii) Conversely, faster factoring is not known to help any lattice problem. Prime-lattice CVP asks for smooth u
    near smooth multiples of N, and N's factors do not help. Ducas–Pierrot (ePrint 2018/146) decode discrete-log
    lattices *by* factoring small integers; that is a different family and the trivial direction.
  - (iv) A poly-time CVP_{1+n^{−ε}} algorithm on prime lattices alone would heuristically factor. We know of no
    evidence that CVP on this family is easier than worst-case CVP.

## 3. Counting relations [HEURISTIC H_S; MODEL]

**Shell estimate** (codex: these are *shell* exponents, not cumulative counts).
- For B = (ln N)^κ, Ψ(X, B) = X^{1−1/κ+o(1)}.
- The expected number of relations in the shell v ≍ N^a, |w| ≍ N^b (u ≍ vN, fixed a, b) is
  N^{e(a,b)+o(1)} with e(a, b) = a(1−2/κ) + b(1−1/κ) − 1/κ.
- Close-vector relations have b < 1 + a. Maximising over admissible shells:

| κ | max shell exponent over b ≤ 1 + a | small residual (b → 0) |
|---|---|---|
| κ < 3/2 | negative for every shell (e(a, 1+a) = a(2 − 3/κ) + 1 − 2/κ < 0) | negative |
| 3/2 < κ < 2 | positive only for a > (2−κ)/(2κ−3) | negative |
| κ > 2 | positive already at a = 0 (b > 1/(κ−1)) | positive for a > 1/(κ−2) (Theorem H) |

- These are heuristic shell thresholds, not proofs of non-existence. The boundary cases and the tail over
  unbounded heights are not treated asymptotically. For fixed B, Ψ grows only polylogarithmically, so the shell
  contributions decay; E0 checks the tail numerically.
- With multiplier targets (v free, Adleman's form), shells are positive for every κ > 1, but need
  N^{(2−κ)/(2(κ−1))} targets.
- For 3/2 < κ < 2, the closest vector at moderate scale C has |w| ≈ N^{1/κ}, smooth with probability N^{−1/κ²}.

**Finite sizes (E0) [MODEL, naive H_S].**
- Primitive and non-primitive pairs (v, w) are summed alike. ln v is truncated at 1.5L + 50, and widening it to 3L
  changed log₁₀E_tot by 0.000.
- Densities are midpoints of the proved cumulative brackets. Codex notes that the *density* lower brackets are 0 at
  these sizes, so the densities are model values, not certified. The cumulative counts are tight, e.g.
  10^{62.20} ≤ Ψ(2^{800}, 499) ≤ 10^{62.35}.

| bits | n | E_tot |
|---|---|---|
| 40 | 47 | 1.9·10^{14} |
| 200 | 47 | 6.6·10^{−11} |
| 400 | 47 (Schnorr) | **1.4·10^{−52}** |
| 400 | 95 | 3.0·10^{−31} |
| 800 | 95 (Schnorr) | **1.5·10^{−115}** |

Resolution sensitivity is 0.03 in log₁₀.

**Model crossover n_min** (least n on our grid with E_tot ≥ n + 1). This is a crossover of this model, not a proved
necessary dimension; n + 1 is the generic dependency count.

| bits | 64 | 128 | 200 | 256 | 400 | 512 | 800 | 1024 |
|---|---|---|---|---|---|---|---|---|
| n_min | 20 | 50 | 90 | 130 | **275** | 400 | **775** | 1125 |
| κ_eff = ln p_{n_min}/ln L | 1.13 | 1.21 | 1.20 | 1.27 | 1.33 | 1.35 | 1.37 | 1.39 |

- The local log-slope of n_min against L is 1.56 (256–800 bits). That is consistent with, not proof of,
  (ln N)^{3/2+o(1)}.
- E1 shows naive H_S ≥ observed *primitive* counts in all 9 toy cells, by factors 6–28. That supports, but does not
  establish, reading E_tot as an over-estimate of what a lattice can find (lattice vectors are primitive). Local
  arithmetic can move predictions either way (codex). The margin at Schnorr's points is ≥ 10^{50}.

## 4. Schnorr 2021, dissected

Text read: "corrected" ePrint 2021/933, PDF dated 16.06.2021 (ePrint page: received 2021-07-09), via web.archive.org.

1. **The claim.**
   - n = 47 / 95 for 400 / 800 bits.
   - One fac-relation per random permutation f.
   - Block size 24; 4.2·10⁹ and 8.4·10^{10} operations.
   - "This destroys the RSA cryptosystem."
   - At 800 bits he derives |u − vN| ≲ 522.9 ≈ p_n^{1.0075}, assuming v ≲ p_n, and writes "ρ(0.007544) = 0.991". In
     the standard normalisation that is ρ(1.0075) ≈ 0.99.
2. **Algebraic failures [PROVED; primary text checked]. These alone invalidate the proof.**
   - Lemma 3.1's proof concludes "u − vN′ < v|ẑ|/(1 − ε/2)". The correct conclusion is
     |u − vN| < (N/N′)·v|ẑ|/(1 − ε/2), so the factor N^{n/(n+1)} is missing.
     - **KNOWN:** Daniel S, crypto.SE answer 88601 (update of 16 March 2021): the bounds are on |u − vN′|, "of very
       different size" from the |u − vN| that must be smooth.
     - Codex (schnorrcvp_1) and we make the missing factor explicit (Theorem R(ii)).
   - "‖b₁‖ < 0.917" is impossible: every nonzero vector of L(R_{n,f}) has a coordinate |f(i)e_i| ≥ 1, or norm
     ≥ N′ ln N.
     - Daniel S (same answer) recomputed the bound as ≈ 5.67·10^{10}.
     - The "norm ≥ 1" argument is codex's; no prior statement was found by the critic.
   - (3.3) "‖b‖² ≥ ln uv + ẑ²" holds for diagonal √ln p_i, not for the displayed f(i). A coordinate with f(i) = 1 and
     p_i > e contradicts it.
   - Context: kodlu and Daniel S on a determinant error in the March 3 version; Couteau (88590) on reliance on flawed
     pruning estimates.
3. **What exists [MODEL, H_S].**
   - His 800-bit box (v ≤ 499, |w| ≤ 523) holds an expected **3.7·10^{−176}** relations. The analogous 400-bit box
     (v ≤ 211, |w| ≤ 250) is our construction, not his: 9.2·10^{−85}.
   - Over all sizes the totals are 1.5·10^{−115} and 1.4·10^{−52}.
   - The mechanism is |uw| ≥ N − 1: with |w| ≤ 523, u is a 499-smooth integer ≥ N/523. Their density near 2^{800} is
     3.6·10^{−180} (model), and there are 10^{62.2–62.35} of them below 2^{800} (bracket).
   - This is heuristic. It says that, under the standard smoothness model, the required relations should not exist
     at these parameters. Combined with item 2, the claim has neither a valid proof nor a plausible mechanism.
4. **SVP embedding.**
   - Ducas tested the **March 3 version** (scale N, not N′): 0 relations in 1000 trials at 400 bits with n = 47.
     Exponent vectors were "typically very sparse, say a single 1 and a single −1".
   - Only last coefficient ±1 is useful, and the paper never enforces it: ccorn (crypto.SE 88647) and Vera's
     Remark 6.
   - Short vectors with last coefficient 0 correspond to smooth near-twins (Mulder–Sterner–van Woerden 2025 §1:
     "the shortest vectors in the prime number lattice often correspond to smooth twins").
   - "Shortest vectors always have last coefficient 0" is not a theorem (codex). It is the typical behaviour at
     Schnorr's parameters.
5. **Cost with lattice algorithms [SCENARIO].**
   - **Asymptotic.** If CVP with γ = 1 + O(1) costs 2^{Θ(n)} and n ≥ n_min ≈ (ln N)^{3/2+o(1)}, the cost is
     exp((ln N)^{3/2+o(1)}). With multiplier targets it is exp(Θ(L lnL/lnlnL)). This is a scenario under stated
     assumptions, not a lower bound for specialised algorithms.
   - **β-dependence.** BKZ-β has root-Hermite factor δ_β ≈ (β/2πe)^{1/(2(β−1))}. Under a GSA-type scenario,
     Babai/enumeration on a BKZ-β basis gives residuals that blow up like N^{Θ(γ_β²)} unless β approaches n.
     - Codex warns that a root-Hermite factor is not a CVP ratio relative to the actual target distance. So this is a
       sensitivity scenario only.
     - At toy n ≤ 40 (E2), LLL-Babai outputs had median |w| ≈ N^{0.64}. Ten LLL-Babai outputs had u or v with more
       than 4000 bits.
   - **Concrete (E0b) [SCENARIO].**
     - log₂T = 0.292n + max(0, log₂V − 0.2075n), where V(n) is the number of lattice vectors an ideal unit-cost
       enumerator (increasing |w|) must output to see n + 1 relations.
     - V uses naive H_S and counts non-coprime (u, v) too, so its direction of bias is not established (codex).
     - The sieve's 2^{0.2075n} close vectors per call are an idealisation.

     | bits | best n (coarse grid) | log₂T | log₂V | crude log₂ QS / GNFS |
     |---|---|---|---|---|
     | 128 | 100 | 42.8 | 34.4 | 28.8 / 33.7 |
     | 256 | 300 | 104.2 | 78.8 | 43.7 / 46.7 |
     | 400 | 600 | 186.7 | 136.0 | 57.0 / 57.2 |
     | 512 | 800 | 264.1 | 196.5 | 65.9 / 63.9 |
     | 800 | 1600 | 467.4 | 332.2 | 85.4 / 77.9 |
     | 1024 | 1600 | 699.1 | 563.9 | 98.5 / 86.8 |

     - These are scenario numbers. The QS/GNFS values are uncalibrated L-expressions with the o(1) dropped.
6. **Other literature on the required dimension.**
   - Grebnev et al. (IEEE Access 11 (2023), arXiv 2303.04656) state that no estimate of the required n as a function
     of N had been derived "neither theoretically, nor empirically". §3/E0 supplies a model estimate.
   - Empirical scaling, all finding exponential behaviour:
     - Aboumrad–Widdows–Kaushik (arXiv 2308.07804, Figs 1, 3, Table 1);
     - Sato–Auzemery–Katayama–Yasuda (IWSEC 2024);
     - Khattar–Yosri (arXiv 2307.09651; ≤ 70 bits);
     - Ritter–Rössner (ePrint 1997/008; 60 bits).
   - Tesoro et al. (Phys. Rev. A 113, 032418 (2026)) claim numerical polynomial scaling up to 100 bits. Their
     Yan-type variant tests smoothness of w over an enlarged base (B₂ ≈ n² primes), which changes the count. Our model
     with B₂ ≫ B would raise E_tot, and has not been run. It is not a contradiction at ≤ 100 bits, where E0 already
     shows modest n suffices.
   - Ducas (RISC slide 21): "A solid average-case complexity analysis is still missing". §3 is a step toward it, not a
     proof.

## 5. Rigour map (Q2)

| Step | Status |
|---|---|
| Theorem R; λ₁(π(Λ_N)) bounds | PROVED |
| Theorem P + Corollary P1 (+ rounding, in the stated regime) | PROVED |
| Distinctness via primes k > B | PROVED |
| Lemma LA (kernel basis; odd index + ERH ⇒ factor) | PROVED (conditional) |
| Witness for every target (smooth u within ≈ B/√L of smooth multiples of kN) | HEURISTIC H_S; a smooth-numbers-in-short-intervals statement far beyond known results; Schnorr 1991's "Hypothesis" |
| Odd index | HEURISTIC H_LA (strong sufficient condition) |
| Benign-oracle residual smoothness | HEURISTIC H_res′ (E2: ratio 0.87 pre-registered, 1.29 normalised) |
| Theorem G | MODEL |
| Conjecture A (tightness) | CONJECTURED |
| E0 counts / n_min | MODEL (naive H_S) |
| BKZ / E0b costs | SCENARIO |
| Prop NP | PROVED from cited theorems (folklore) |

**Cleanest provable statement.** Let N not be a prime power, with no prime factor ≤ B, B ≥ 2 ln²N, and assume ERH.
Assume (i) for each of n + 2 primes k ∈ (B, 4B] a squarefree witness exists for target kN, with A₀ ≤ poly(L) and
|w₀| ≤ ω_{⌊A₀⌋}; and (ii) the relations generate an odd-index sublattice of Λ_N. Then an exact-CVP oracle on the prime
lattices of dimension n = π(B) factors N in deterministic polynomial time with O(n·poly L) calls. Under H_S, (i) needs
B = (ln N)^κ with κ > 2. B ≥ 2 ln²N is only the linear-algebra floor.

## 6. Variants (Q5)

| Variant | Effect | Cost / status |
|---|---|---|
| Multiplier targets ln(kN), v free | shells positive for κ > 1; rigorous distinctness for k > B | N^{Ω(1)} targets when κ < 3/2; SCENARIO exp(Θ(L lnL/lnlnL)) with 2^{Θ(n)} CVP |
| Adleman SVP form u − kN^γ = w | same relations | its distinctness step has a gap (§2.3); same cost |
| ℓ1 / ℓ_p / other weights | log\|w\| stays linear in the norm of e | Theorem P-type upper bound PROVED (Vera's ℓ_p versions); tightness CONJECTURED |
| Enlarged residual base B₂ (Yan, Grebnev) | raises the yield per output | not modelled here; changes constants and dimension, not the 2^{Θ(n)} premise (assessment) |
| NFS rational side via prime lattices | a poly-delay enumerator of dimension L[1/3] would give the factory constant 1.6386 (smoothenum T3) | HYPOTHETICAL; known enumeration is exp(L[1/3]) |
| Algebraic side (S-unit lattice CVP) | toward smoothenum T6/T7 | SPECULATIVE: needs S-unit/class-group data plus principality, dimension L[1/3] |
| Two lattices (u-side and w-side) | w = u − vN is not multiplicative in e, so w-smoothness cannot be imposed linearly | assessment |

No variant is plausibly better.
- Only a *polynomial-time* CVP_{1+n^{−ε}} algorithm on this family gives polynomial-time factoring (Theorem H).
- A 2^{n^c} algorithm at n = L^κ gives exp(L^{κc}).
- No converse reduction is known.

## 7. Experiments

### 7.1 E1: exact relation counts (pre-registered)

Setup: 9 cells × 10 semiprimes; v ≤ 10³, |w| ≤ 10⁵ (2·10⁴ for n = 25).

| (bits, n) | all | primitive | pred H_S′ | ratio | naive | post-hoc H_S″ ratio |
|---|---|---|---|---|---|---|
| (24,10) | 2207 | 74 | 36.0 | 2.05 | 2064 | 0.84 |
| (24,15) | 22186 | 1304 | 729 | 1.79 | 20226 | 0.86 |
| (24,25) | 76176 | 10403 | 6482 | 1.61 | 71604 | 0.93 |
| (28,10) | 99 | 10 | 4.1 | 2.45 | 283 | 0.96 |
| (32,10) | 403 | 2 | 0.43 | 4.64 | 36 | 1.75 |
| (32,15) | 1182 | 27 | 13.5 | 2.01 | 528 | 0.88 |
| (32,25) | 3520 | 356 | 213 | 1.67 | 3038 | 0.89 |
| (40,15) | 4 | 1 | 0.20 | 4.88 | 10.8 | 1.98 |
| (40,25) | 159 | 16 | 5.6 | 2.85 | 100 | 1.41 |

- **E1-P1 passed** (pooled 1.63; every cell with prediction ≥ 10 within [1/3, 3]).
- **E1-P2 passed** (naive ≥ primitive in every cell).
- The ≈ ×2 excess is parity: exactly one of u, v, w is even. POST-HOC H_S″ (`psi.local_factor`, with u uniform on the
  p − 1 admissible residues, so g₂ = 2f₂) gives 0.92.
- E1 does not test squarefree witnesses or windows as short as Theorem H needs (codex). It tests the density model
  only.

### 7.2 E2: CVP-output yield (pre-registered)

Setup:
- 3600 lattices: 6 sizes × 3 N × n ∈ {12, 20, 30, 40} × 5 variants × 10 seeds.
- Outputs per lattice: LLL-Babai, BKZ-10-Babai, and the ≤ 20 closest vectors found. 78,963 outputs in all.

| variant | relations | obs/H_res′ (pre-reg.) | obs/naive | obs/normalised (post-hoc) | post-hoc H_S″ |
|---|---|---|---|---|---|
| S21perm (f(i), N′) | 296 | 1.75 | 1.34 | 1.78 | 1.03 |
| S21sqrt | 379 | 1.49 | 1.09 | 1.53 | 0.91 |
| A0.5 (√ln p_σ(i), C = N^{1/2}) | 2506 | 1.24 | 0.90 | 1.40 | 0.94 |
| A1.0 | 3494 | 0.82 | 0.74 | 1.24 | 0.73 |
| A1.25 | 3113 | 0.68 | 0.73 | 1.21 | 0.62 |
| pooled | 9788 | **0.87** | 0.80 | **1.29** | 0.75 |

- **E2-P1 to E2-P4 passed** as pre-registered:
  - pooled 0.87 ∈ [1/2, 2];
  - refined fits better than naive;
  - yield at 40 bits < yield at 20 bits in every cell;
  - median log₂|w|/bits for S21 at ≥ 32 bits is 0.96 / 0.92.
- **Predictor defect (codex #18; disclosed).**
  - The pre-registered predictor gave 2066 probabilities > 1 (up to 12.7), and 66 relations had p = 0. It multiplied
    real-interval densities by local factors without discrete normalisation.
  - The POST-HOC normalised predictor (`analyze.corrected_e2`) uses exact discrete probabilities for |w| ≤ 2·10⁵:
    smooth and coprime-to-supp(uv) among the window's integers coprime to supp(uv); 1 for |w| ≤ B. Larger |w| is
    clipped. It gives **1.29 pooled** (1.21–1.78).
  - The remaining excess matches E1's parity effect, which the coprime-only predictor omits.
  - Conclusion: residuals of CVP outputs are smooth at roughly the rate of random integers of their size with the
    forced local conditions, within ×1.8. E2 tests smoothness given the observed residual. It does not test whether
    approximate answers saturate Theorem P's bound.
- **Yield per list output.**
  - n = 40, A1.0: 0.82 / 0.66 / 0.50 / 0.40 / 0.16 / 0.085 at 20–40 bits.
  - n = 40, S21sqrt: 0.18 → 0.000.
  - Exact-CVP (list0) residuals for A-variants have median |w| ≈ N^{0.475}.
- **E2-P5.** 700 list0 relations, of which 416 distinct (59%). Random weights do not guarantee distinctness (as
  Aboumrad–Widdows–Kaushik observed). Prime multipliers do (§2.3).
- **E2-P6 (correctness only).**
  - All 15 N of 20–36 bits were factored from pooled relations. Factors are STRUCTURAL: N's prime factors exceed
    p₄₀ = 173.
  - At 40 bits, only 17–22 distinct relations were found (< 42).
- **Certification of list0 (codex #19).**
  - list0 is "closest found", not certified exact. Targets enter fplll's enumeration as binary64.
  - Spot check (`check_list0.py`, 100 random lattices) against fplll `CVP.closest_vector(method="proved")`:
    - 91 agree;
    - in 8, list0 is strictly *closer*, so fplll's "proved" mode is not exact on these large-entry lattices;
    - in 1 (40 bits, n = 30, A1.25), list0 missed the closest vector (distance ratio 2.53).
- **Deviation.** The first E2 run crashed on a Babai output with u > 10^{4300}. We patched serialization and reran
  deterministically; the crashed log is kept.

### 7.3 Commands

- `.venv/bin/python -m src.schnorrcvp.run_e0`
- `run_e0b`
- `run_e0b_ext`
- `run_e1`
- `run_e2`
- `check_list0`
- `analyze e1|e2|posthoc_e1|posthoc_e2|corrected_e2`

E0b, its extension, `corrected_e2`, `check_list0` and the post-hoc models were added after the pre-registration.

## 8. Coordination

**suffmap** (H9: "sufficiency for any approximation factor not established; details: schnorrcvp"). Our answer:
- sufficient for polynomial time at γ = 1 + n^{−ε} (Theorem H; heuristic);
- sufficient for L[α < 1/3] at γ = 1 + (log n)^{−c}, c > 1 (benign oracle, Theorem G);
- for γ fixed > 1, L[1/2, √(γ²−1)] (benign);
- nothing guaranteed at growing γ.
- Near-exact prime-lattice CVP is an instance of suffmap's E3 (Bohr-set lister over all smooth multipliers v), so it
  inherits Prop DLT (it would also give prime-field DL, as Schnorr 1991's title suggests).
- **Consistency with suffmap's Prop SI (1.2114).** Prop SI fixes v = 1, so its residual must span the smooth-number
  gap 1/d(kN), which sets the constant. Schnorr's lattice ranges over N^{o(1)} smooth v, so witnesses with |w| ≤ B
  exist and the floor is paid by the γ = 1 + o(1) oracle instead.

**abcsunit.** Relations are S-unit solutions of u − vN = w. Evertse 1984 bounds λx + μy = 1 by 3·7^{d+2s} (s = |S|
including ∞), i.e. 3·7^{2s+3} over Q with s finite primes. That is 2^{O(n)}, far above n + 1, so finiteness is no
obstruction. The obstruction is the small-height count of §3.

## 9. Reviews

### 9.1 Codex blind round (schnorrcvp_1)

Codex independently found:
- the identity and certificate;
- the N = 35 counterexample to "exact CVP always returns a relation";
- Schnorr 1991's conditions;
- the explicit missing factor in Lemma 3.1 and the 0.917 impossibility;
- |uw| ≥ N − 1.

### 9.2 Codex adversarial round (schnorrcvp_2)

**FATAL (to the specific claims).**
- (a) Prop A's tightness argument. Withdrawn; now Conjecture A.
- (b) The λ₁(Λ_N) lower bound with the sign coordinate. Fixed via the projection π(Λ_N).

**MAJOR, all addressed.**
- Corollary P1's constants (fixed with the proof critic).
- The rounding regime (now stated).
- Lemma LA: odd index presented as standard (now "strong sufficient"), and Bach's constants (fixed: Thm 2 needs trial
  division).
- Kernel truncation at 64 in code (fixed; regression test).
- The Adleman repair overclaimed (now "different construction").
- Theorem H's per-target uniformity and the limiting dimension (now fixed κ; uniformity is in the hypothesis).
- Theorem G's θ and oracle cost (fixed).
- The false discontinuity at γ = 2 (fixed: L[1/2, √3]).
- Shell vs cumulative counts (fixed).
- "All sizes", the upper-bound reading and n_min as necessity (now: model, crossover).
- Ψ-bracket stability vs validity (stated).
- Over-reach in the Schnorr critique (algebraic errors now lead; statistics qualified).
- The BKZ lower bound and E0b's bias direction (now SCENARIO).
- E2 probabilities > 1 (normalised re-analysis, 1.29).
- E2 used to validate untested oracle assumptions (restricted).
- Universal barrier claims (removed).

**MINOR.**
- Prop NP notation (Turing reduction of a function; Aharonov–Regev C√n; coAM at O(√(n/log n))).
- Certified logarithms.
- list0 certification (spot-checked: 91/100 agree, 8 better than fplll, 1 miss).
- Literature positioning.

### 9.3 Opus critics

**Proof critic.** No FATAL issues. The core mathematics was checked line by line. Its two MAJOR points (overclaims;
Theorem G's oracle-cost constants) and twelve MINOR points are applied (version 1.1, commit 3f6fbd7). Its λ₁ seminorm
point coincides with codex's FATAL item (b).

**Literature critic (sources read in primary form are listed in its report).**
- The Lemma 3.1 / 0.917 errors were public in March 2021 (Daniel S, crypto.SE 88601). Now credited.
- Schnorr 1991 Thm 5 already contains the witness count, and Thms 2 + 5 give the ℓ1 analogue of Theorem H. Now
  credited, with the new content restated.
- Schnorr 1991 is an ℓ1 lattice with weights ℓ_i; the √ℓ_i weights are Adleman's. Fixed.
- Ducas tested the March 3 version, and his observation is only sparsity. Fixed; ccorn, Vera and MSW credited.
- Prop NP is folklore (SchnorrGate README; Micciancio). Fixed.
- Adleman quote location (§1) and Ducas slide number (22). Fixed.
- Missing empirical literature (Grebnev et al.; AWK; Sato et al.; Tesoro et al.; Khattar–Yosri). Added (§4.6).
- The 400-bit box is ours. Fixed.
- Reference hygiene. Fixed.
- Newer claims (OpenAI 2026 manuscript Ch. 7; Hair–Sahai). Added as unrefereed.
- Bach Thm 2 vs Thm 3. Fixed.
- Found no prior statement of: the Adleman distinctness gap; the E_tot / n_min counts or a (ln N)^{3/2} threshold;
  the "norm ≥ 1" argument. Not checked: Kan–Shen 2003.

## 10. References

- C. P. Schnorr, *Factoring integers and computing discrete logarithms via diophantine approximation*, EUROCRYPT '91,
  LNCS 547, pp. 281–293; DIMACS 13 (1993), pp. 171–181. Read in primary form by the literature critic.
- L. M. Adleman, *Factoring and lattice reduction*, draft, 16 March 1995 (SchnorrGate `LATFAC.pdf`). Read.
- C. P. Schnorr, *Fast factoring integers by SVP algorithms, corrected*, ePrint 2021/933 (PDF 16.06.2021; replaces the
  withdrawn 2021/232; Eurocrypt 2009 rump-session precursor). Read.
- Daniel S et al., crypto.stackexchange.com/questions/88582 (answers 88590, 88601, 88609, 88647), March 2021.
- L. Ducas, SchnorrGate (GitHub: README, purpose.md, fac.sage); RISC slides, 20 May 2021. Read.
- B. Vera, *A note on integer factorization using lattices*, arXiv 1003.5461 (2010).
- Aboumrad, Widdows, Kaushik, arXiv 2308.07804; Grebnev et al., IEEE Access 11 (2023) 134760 / arXiv 2303.04656;
  Sato, Auzemery, Katayama, Yasuda, IWSEC 2024, LNCS pp. 142–157; Tesoro et al., Phys. Rev. A 113, 032418 (2026) /
  arXiv 2410.16355; Yan et al., arXiv 2212.12372; Khattar–Yosri, arXiv 2307.09651; Ritter–Rössner, ePrint 1997/008;
  Mulder–Sterner–van Woerden, arXiv 2509.17699; Ducas–Pierrot, ePrint 2018/146; Kan–Shen, PDCAT 2003 (not checked).
- van Emde Boas 1981; Dinur–Kindler–Raz–Safra, Combinatorica 23 (2003); Ajtai, STOC 1998; Micciancio, SICOMP 2001;
  Khot, JACM 2005; Haviv–Regev, STOC 2007 / ToC 8 (2012); Goldreich–Goldwasser, JCSS 2000; Aharonov–Regev, JACM 2005;
  Regev, STOC 2005 / JACM 56(6) 2009; Bach, Math. Comp. 55 (1990) 355–380 (Thms 2, 3); Evertse, Invent. Math. 75
  (1984). Complexity statements were checked against the critic's knowledge and the cited abstracts, not all against
  full texts.
- OpenAI, "Ten Advances in Mathematics and Theoretical Computer Science" (updated 6 Aug 2026), Ch. 7: claimed, not
  refereed, not checked. Hair–Sahai, ePrint 2025/2181: claimed.
