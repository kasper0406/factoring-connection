# detGRH: deterministic subexponential factoring under GRH? Multiplicative square roots, the Z route, and what is left

Agent `detGRH`, phase 13, 2026-10-20.
- Pre-registration: `notes/detGRH_prereg.md` (commit `7cab676`, before any code or run).
- Code:
  - `src/detGRH/core.py` (operational: sees only N);
  - `src/detGRH/diag.py` (harness/diagnostic: uses the factorisation);
  - `src/detGRH/harness_e1.py`, `src/detGRH/analyze.py`, `src/detGRH/posthoc_e1.py` (post hoc, labelled).
- Tests: `tests/test_detGRH.py`.
- Data: `experiments/raw/detGRH_e1.jsonl`, `experiments/analysis/detGRH_e1_summary.json`; config
  `experiments/configs/detGRH_e1.json`.
- Reviews:
  - GPT-6.1-sol blind round: `codex/detGRH_1.{md,out}`, full report copied to `notes/detGRH_codex1_report.md`.
  - GPT-6.1-sol adversarial referee: `codex/detGRH_2.{md,out}`, full report copied to
    `notes/detGRH_codex2_report.md`.
  - Two independent read-only Opus critics, one on proofs and one on literature (§9).
  - All fixes are applied in this version.
- Labels:
  - **proved**; **proved (folklore-level)** for elementary combinations of classical results;
  - **conditional**; **open**;
  - **refuted** (for a stated hypothesis);
  - **observation**: a limitation of a specific proof method, not a lower bound;
  - **experimentally supported**.

## Verdict (plain language)

**No. We did not obtain a deterministic subexponential algorithm for complete factorisation under GRH, for all N or
for almost all N.**

The orchestrator's idea does not survive as stated. Analysing it does locate the difficulty: **under GRH,
fruitfulness is not the real obstruction. Deterministic relation supply is.**

1. **The multiplicativity idea [proved].**
   - The implication is valid. If the algorithm's square root y(R) were multiplicative (even modulo ±1) on all
     squares mod N, GRH would give a fruitful seed below 2 log² N.
   - Its hypothesis can hold only in a degenerate case.
     - **Lemma MR.** A multiplicative square-root function on (Z/N)^{×2}, even modulo ±1, exists iff every p | N is
       ≡ 3 (mod 4). It is then unique: the principal root. Computing it factors N deterministically under GRH.
     - **Proposition CD.** The square roots that index calculus or NFS produce are multiplicative on *certificates*,
       not on residues. All admissible roots of one residue agree modulo ±1 iff G_B ∩ μ₂ = {±1}.
     - **Theorem M(G).** Suppose G_B ∩ μ₂ ⊋ {±1}, and the root function is honest on small prime squares
       (y(r²) = r). Then y cannot be multiplicative. Some pair (R, R′) gives y(RR′) ≢ ±y(R)y(R′) mod N, and that
       pair is itself a factor.
     - These statements are about abstract total root functions. Theorem A's actual descent is partial: it rejects
       small seeds outright (Lemma SR). Nothing locates a fruitful pair in a short list.
2. **Fruitlessness over the whole factor base [proved].**
   - Every relation set over a factor base is fruitless iff G_B ∩ μ₂(Z/N) = {±1}, where G_B is the subgroup the
     factor base generates.
   - lemmaM's barrier B1 is an instance of this, within its degree-1 two-form scope.
   - Under GRH this case is impossible once the factor base contains the primes below 2 log² N.
   - What remains is a property of the relation *family*. Proposition FL gives a sufficient criterion: fruitful if
     the 2-rank of Λ/Λ_c is < ω − 1. Proposition FL(ii) gives the exact one.
3. **The Z route [proved, folklore-level: Miller 1976 + Bach 1990].**
   - Do the linear algebra over Z instead of F₂. Suppose the kernel relations ∏ q^{e_q} ≡ 1 (mod N) have full Q-rank
     over a factor base that contains the primes < 2 log² N, with no prime factor of N there. Then any nonzero
     maximal minor is a multiple of φ(N), and Miller's ERH test factors N deterministically.
   - No F₂ square roots, no Φ, no quadratic characters and no hidden roots are needed.
   - The orchestrator's "multiplicative χ plus small seed" is realised correctly by Miller's *power map*, not by an
     F₂ square root.
   - Conversely, a factoring algorithm yields such relations trivially (λ(N)e_q). So under ERH, deterministic
     factoring is equivalent to deterministically producing an *arbitrary* full-rank kernel-relation set
     (Proposition EQ).
   - It is **not** shown to be equivalent to producing such relations from a prescribed NFS source (Hypothesis
     DS-NFS). Only DS-NFS ⇒ factoring is proved, and that needs an extraction lemma (XK-Z) that we state but have
     not proved.
4. **Derandomizing Theorem A [assessment, with proved items; §5].**
   - Most random choices can be replaced deterministically, unconditionally or under GRH:
     - kernel sampling: deterministic F₂ elimination at L^{3a} or L^{ωa};
     - character primes: GRH for the Kummer extensions;
     - Hensel prime and finite-field routines: deterministic at fixed degree;
     - the separate shape-test subroutine: recursion, *if* the outer algorithm is deterministic;
     - the trivial-gcd restarts: removed by the Z route.
   - The seed's *fruitfulness* role is replaced by the Z route, conditional on XK-Z and on full projected rank. The
     deterministic exact linear algebra costs L^{(ω+1)a+o(1)}. That fits Theorem RA's slack iff ω < 1 + √2 (true for
     ω = 2.371) and Theorem A's iff θ < 1/ω.
   - **Unresolved** (no known derandomization; not a lower bound): the field law (g₀, g₁, w, c) and the seed law, in
     their **relation-supply** roles.
5. **Field and seed derandomization [observations; §6].**
   - The existing first and second moments (FM, LD*) work by averaging over coefficient boxes and multiplier sets of
     size N^{Ω(1)}. Their analytic input (smooth numbers in full-length progressions at level ≤ 1/2) does not apply
     to short boxes.
   - This is a limitation of that proof, not a lower bound for every explicit sub-family.
   - Even granting a square-root cancellation bound for the target set (not known, even under GRH), the natural
     Fourier certificate that an explicit seed list hits the descendable set only works for lists of N^{1+o(1)}
     distinct seeds.
   - Checkability of success does not help without a new pointwise or short-family smooth-value theorem. We know of
     none, conditional on GRH or not.
6. **Literature [§8].**
   - We found no deterministic subexponential complete-factoring algorithm under GRH, and no impossibility result.
   - The deterministic GRH record is N^{1/5+ε}: Shanks 1971 and Schoof 1982 §6. It is now unconditional (Harvey 2021;
     Harvey–Hittmeir 2022).
   - Deterministic subexponential results exist only with oracles (Żrałek 2010: L[1/3] given iterated φ) or under
     circuit-hardness hypotheses (Impagliazzo–Wigderson, polyimpl Theorem DR).
   - The Z route's ingredients are classical: Miller 1976; Bach 1984 §3 (deterministic under ERH); Bach 1990;
     Shoup §§10.4–10.5. The closest prior formulation is Stange 2022 (heuristic, randomized). Proposition EQ is a
     folklore-level reformulation.
   - Finding *one* factor of almost all integers is trivial (trial division to log N). The open "almost all" question
     concerns complete factorisation, or a hard ensemble such as balanced semiprimes.
7. **Experiments (pre-registered E1; toy sizes, correctness only) [§7].**
   - 240 instances at 28–40 bits, plus an exhaustive check over N < 3000.
   - Every exact prediction held: Lemma MR (0/1499 violations), Proposition FL, the cocycle identity, and the Z route.
     The Z route factored 237/237 full-rank instances, including the doubled control in which every F₂ dependency is
     fruitless by construction.
   - The orchestrator's multiplicativity hypothesis is refuted: re-certified products of squares disagree with
     products of roots at rate 0.500 for semiprimes and 0.746 for ω = 3, and every disagreement is a factor.
   - One pre-registered heuristic window failed (P-CY3: non-diagonal seed rate 0.24). Post hoc, about half of the
     seeds' canonical completions are their own honest relation.

**The most promising open statement** is Hypothesis DS-NFS + Lemma XK-Z (§5): an explicit, deterministic list of
L^{O(1)} NFS fields and seeds whose descents, after exact elimination of the algebraic side, give rational kernel
relations of full Q-rank. Under ERH this implies deterministic L[1/2, O(1)] factoring.
- The new idea it needs is a pointwise or short-family smooth-values theorem together with a spanning statement.
- F₂ fruitfulness is no longer part of the problem.
- It is not easier than factoring in any proved sense: Proposition EQ shows the unrestricted version is equivalent.

---

## 1. Setting and notation

- **Input.** N is odd and not a prime power, ω = ω(N) ≥ 2, and N has no prime factor ≤ B. Preprocessing removes 2,
  primes, perfect powers and every prime factor ≤ max(B, 3 log² N).
- **Factor base.** FB = (−1, p₁, …, p_{n−1}), with p_i the primes ≤ B.
- **Groups.**
  - Γ := Z^n, with −1 a free generator. ev: Γ → U := (Z/N)^×, ev(e) = ∏ q^{e_q}.
  - Λ := ker ev ∋ 2e_{−1}. G_B := ev(Γ) ≅ Γ/Λ. μ₂ := U[2] ≅ {±1}^ω.
  - Λ_c ⊆ Λ is the lattice of found relations. It always contains 2e_{−1} (the sign control).
- **Generality.**
  - The F₂ statements hold for any finitely generated certificate group with a homomorphism ev; replace {±1} by
    ev(Γ_tors).
  - An NFS certificate group Q^×_S × K^×_S enters the Z route only after exact elimination of the algebraic side
    (§5, Lemma XK-Z).
- **Square dependencies and Φ.**
  - A *square dependency* of Λ_c is v ∈ Λ_c ∩ 2Γ, and Φ(v) := ev(v/2) ∈ μ₂.
  - The reduced map Φ̄ := Φ mod {±1} takes values in μ₂/{±1} ≅ F₂^{ω−1}.
  - v is *fruitful* iff Φ̄(v) ≠ 0, iff gcd(Φ(v) − 1, N) is proper.
  - Φ is never trivial: Φ(2e_{−1}) = −1. Fruitfulness statements are therefore about Φ̄.
  - In this language, Lemma H (nfsfruit) says Φ̄ is a homomorphism on Λ_c ∩ 2Γ that vanishes on 2Λ_c.
- **Bach's bounds (ERH for Dirichlet L-functions mod N).**
  - Every proper subgroup of U omits a positive integer < 2 log² N, which may divide N.
  - If the omitted element must be a prime ℓ ∤ N, the bound is 3 log² N.
  - Source: Bach, Math. Comp. 55 (1990) 355–380. The exact theorem numbers are unverified by us; the constants are
    as summarised by Lamzouri–Li–Soundararajan, Math. Comp. 84 (2015), §1.1, and Shoup, *A Computational
    Introduction*, §10.5.
  - LLS Thm 1.1(2) sharpens this to (log q)² under extra conditions.
  - **Consequence.** After the preprocessing above, if FB ⊇ {primes < 2 log² N}, then G_B = U and G_B ∩ μ₂ = μ₂.

## 2. Multiplicative square roots and the orchestrator's χ

Let Q := U², and σ: U/{±1} → Q, [u] ↦ u².

**Lemma MR [proved].**
- (a) A homomorphism y: Q → U with y(s)² = s for all s exists iff every p | N is ≡ 3 (mod 4). It is then unique,
  and it is the principal root (the root lying in Q).
- (b) For ω ≥ 2, the same holds for homomorphisms y: Q → U/{±1} with σ(y(s)) = s. For ω = 1, σ is an isomorphism.

*Proof.* Write U = ∏_i C_{n_i} by CRT, where n_i = p_i^{e_i−1}(p_i − 1) is even, and n_i/2 is odd iff p_i ≡ 3 (mod 4).
- **(a), only if.** If some n_i/2 is even, C_{n_i} has an element a of order 4. Then q := a² ∈ Q has order 2, so a
  homomorphism gives y(q)² = y(q²) = 1 ≠ q. Contradiction.
- **(a), if.** If all n_i/2 are odd, Q = ∏ C_{n_i/2} has odd order. Squaring is then an automorphism of Q, and its
  inverse (followed by the inclusion into U) is a section.
- **(a), uniqueness.** Two sections differ by a homomorphism Q → μ₂, which is trivial because |Q| is odd.
- **(b), the map χ.** Put χ(r) := y(r²)·[r]^{−1} ∈ μ₂/{±1}. Then χ is a homomorphism on U.
- **(b), its restriction to μ₂.** Because (εr)² = r², we get χ(εr) = [ε]χ(r). So χ restricted to μ₂ is the quotient
  map μ₂ → μ₂/{±1}.
- **(b), the contradiction.** Suppose some n_i/2 is even, and let ε^{(i)} be −1 at factor i and 1 elsewhere. It is a
  square (ε^{(i)} = (g_i^{n_i/4})²), and every homomorphism into an exponent-2 group kills squares. So χ(ε^{(i)}) = 0.
  But [ε^{(i)}] ≠ 0 when ω ≥ 2.
- **(b), uniqueness.** As in (a), since |Q| is odd. ∎

*Checks.*
- E1a: all 1 499 odd N < 3000, plus an independent backtracking count for N < 160 (test).
- The Opus critic counted 30 further N, including prime powers.
- Codex detGRH_1 checked 150 498 identities for N ≤ 250.
- We found no published statement of MR. It is elementary: Rabin–Williams use the principal root on Blum integers.

**Corollary MR [proved; Dirichlet GRH].** Assume y as in MR(b), with ω ≥ 2. Then:
- all p | N are ≡ 3 (mod 4);
- y is the principal root modulo ±1;
- χ(r) = ((r/p₁), …, (r/p_ω)) modulo the diagonal;
- some integer 2 ≤ r < 2 log² N either has gcd(r, N) > 1 or satisfies χ(r) ≠ 0, i.e. gcd(y(r²) − r, N) is proper.

*Proof.* χ: U → μ₂/{±1} is a nonzero homomorphism (MR(b)). Its kernel is proper, and Bach gives a small element
outside it. ∎

So the orchestrator's implication is **valid**. Its hypothesis forces the principal root, which factors N
deterministically under GRH. This is Rabin–Williams; cf. also Jeřábek 2016 (§8). The question is whether
index-calculus roots are multiplicative (§3).

## 3. Certificate-level multiplicativity; the full-lattice functional

**Certificates and the canonical convention.**
- Put C := 2Γ + Λ_c. A *certificate* of R is c ∈ C with ev(c) = R.
- A *decomposition* is c = 2g + λ with λ ∈ Λ_c. The *admissible roots* of c are the values ev(g).
- **Canonical convention** (codex detGRH_2 §3; this is what `core.canonical_root` computes).
  - Fix b_j ∈ Λ_c whose reductions form a basis of the image of Λ_c in Γ/2Γ.
  - Write c ≡ Σ a_j b_j (mod 2), with a_j ∈ {0, 1}, and set Y(c) := ev((c − Σ a_j b_j)/2).
  - Y is a homomorphism C → U: adding parity coefficients produces only carries 2b_j, whose halved values evaluate
    to 1.
  - Also Y(2g) = ev(g) and Y(c + 2h) = Y(c)·ev(h).

**Proposition FL (full lattice) [proved].**
- (i) Φ_full: (Λ ∩ 2Γ)/2Λ → G_B ∩ μ₂, v ↦ ev(v/2), is an isomorphism, and [2e_{−1}] ↦ −1. Write
  W′ := (Λ ∩ 2Γ)/2Λ; dim W′ = rk₂(G_B ∩ μ₂).
- (ii) Let V_c be the image of Λ_c in Λ/2Λ. Then Φ(Λ_c ∩ 2Γ) = Φ_full(V_c ∩ W′). Hence Λ_c has a fruitful dependency
  iff dim(V_c ∩ W′) ≥ 2. **This is the exact criterion.**
- (iii) If Λ_c has full rank and c := dim Λ/(Λ_c + 2Λ) ≤ dim W′ − 2, then Λ_c has a fruitful dependency.
  **Sufficient, not necessary.**

*Proof.*
- (i) Well-defined: v ∈ 2Λ gives ev(v/2) = 1. Injective: ev(v/2) = 1 gives v/2 ∈ Λ. Surjective: for
  ε = ev(g) ∈ G_B ∩ μ₂, the vector 2g lies in Λ ∩ 2Γ.
- (ii), the inclusion ⊆: a v ∈ Λ_c ∩ 2Γ has a class in V_c ∩ W′ with the same Φ-value.
- (ii), the inclusion ⊇: a class in V_c ∩ W′ has representatives v ∈ Λ_c and w ∈ Λ ∩ 2Γ with v − w ∈ 2Λ ⊆ 2Γ. So
  v ∈ Λ_c ∩ 2Γ and Φ(v) = Φ_full(w), since (v − w)/2 ∈ Λ.
- (ii), the criterion: [2e_{−1}] ∈ V_c ∩ W′ is nonzero and maps to −1.
- (iii) dim(V_c ∩ W′) ≥ (n − c) + dim W′ − n. ∎

*Checks.* The Opus critic tested 600 random Λ_c for N with ω = 2, 3: 0 mismatches with brute force. In E1, (iii)'s
premise held in 42 tri instances, with 0 violations.

**Proposition CD (certificate dependence) [proved].**
- (a) The admissible roots of a fixed c are exactly ev(g₀)·Φ(Λ_c ∩ 2Γ), for any one admissible root ev(g₀). If Λ_c
  is fruitless (Φ̄ ≡ 0), the root class modulo ±1 is intrinsic, and it is a homomorphism C → U/{±1}.
- (b) For certificates c, c′ of the same R, with the canonical convention, Y(c)/Y(c′) = Φ_full(·) of a vector in
  Λ ∩ 2Γ. As c′ ranges over the certificates of R, this ratio ranges over all of G_B ∩ μ₂.
- Hence: **all admissible roots of R, over all certificates, agree modulo ±1 iff G_B ∩ μ₂ ⊆ {±1}.** The same holds
  for the canonical convention Y.

*Proof.*
- (a) Two decompositions differ by λ′ − λ = 2(g − g′) ∈ Λ_c ∩ 2Γ, and every element of Λ_c ∩ 2Γ arises this way.
- (b) c − c′ ∈ Λ. Taking c′ := c + 2h with ev(h) ∈ G_B ∩ μ₂ realises every value, since Y(c + 2h) = Y(c)ev(h). ∎

*Caveat (codex detGRH_2 §3).* "Iff" fails for an arbitrary *selected* root function. For N = 77 and FB = (−1, 2, 3),
selecting the principal root for every certificate is admissible (the certificate 2e₂ admits both 2 and 9), and that
function does descend to residues.

**Theorem M (dichotomy) [proved; abstract].** Exactly one of the following holds.

**(S) G_B ∩ μ₂ = {±1}.** Then:
- every square dependency of every Λ_c ⊆ Λ is fruitless;
- every admissible root class is intrinsic and multiplicative on G_B², modulo ±1;
- y(r²)/r ∈ {±1} for every r ∈ G_B;
- Miller's test with unit bases in G_B and any multiple of exp G_B never splits N.
  - *Reason (corrected):* the Miller sequence lies in G_B, but need not lie in μ₂. For N = 65, FB = (−1, 2) and base
    2, the sequence 8, 64, 1 starts with an element of order 4. Its last non-identity element is a square root of 1
    in G_B, hence in G_B ∩ μ₂ = {±1}.
- **lemmaM's B1 is an instance**, for a Blum semiprime N = pq (distinct p, q ≡ 3 mod 4) in the even-weight
  two-form degree-1 scope.
  - All FB primes have (q/N) = +1, so G_B ⊆ ker J_N.
  - J_N(ε_p, ε_q) = ε_pε_q, and −1 ∈ G_B, so G_B ∩ μ₂ = {±1}.
  - Checked by the critic for N = 2641, 6169, 18001.
  - It does not extend to all N whose primes are ≡ 3 (mod 4) with multiplicities, nor to QS or hidden-root
    relations.
- Under ERH and the preprocessing of §1, (S) is impossible once FB ⊇ {primes < 2 log² N}.

**(G) G_B ∩ μ₂ ⊋ {±1}.**
- Fruitfulness of Λ_c is then decided by Proposition FL(ii).
- Suppose moreover G_B = U (ERH, as above), and F is any total root function on Q with:
  - **root-value honesty** on small prime squares: F(r²) ≡ r (mod ±1) for every prime r < 2 log² N (for example,
    the canonical Y applied to c(r²) = 2e_r);
- Then F is not multiplicative modulo ±1. So some R, R′ have δ := F(RR′)/(F(R)F(R′)) ∉ {±1}, and
  gcd(F(RR′) − F(R)F(R′), N) is proper.

*Proof of (G).*
- If F were multiplicative, Corollary MR would give a prime r < 2 log² N with χ(r) ≠ 0, which contradicts honesty.
- For non-Blum N, MR alone already excludes multiplicativity. ∎
- Certificate honesty (c(r²) = 2e_r) gives root honesty only under the canonical convention: the certificate 2e₂ mod
  77 also admits the root 9.

*Non-vacuity (codex detGRH_2 §4).*
- Take N = 1019·1031, where −1 and 2 generate U.
- Let F be r on small prime squares and the principal root elsewhere.
- Then F(4) = 2, F(9) = 3 and F(36) = 1 049 564, with gcd(F(36) − 6, N) = 1031.

*Application to Theorem A: not available [open].* Theorem A's DESCEND_ONCE is partial. It does not certify every
square, and it rejects small seeds (Lemma SR). Theorem M(G) is a statement about abstract total root functions.
Extending the descent to all squares at controlled cost is not supplied.

**Lemma SR (small seeds are rejected) [proved; codex detGRH_1 §2.3].** If r ≥ 2 and r⁴ + 1 < N, the shortest
nonzero vector of {(u, v) : u ≡ r²v (mod N)} is ±(r², 1).
- *Proof.* If u − r²v = kN ≠ 0, Cauchy–Schwarz gives ‖(u, v)‖ ≥ N/√(1 + r⁴) > √(1 + r⁴). If k = 0 the vector is a
  multiple of (r², 1).
- So Theorem A's shape test rejects every seed r ≤ L^{o(1)}, because r² is neither squarefree nor composed of
  J-primes ≈ N^{2ν}.
- Accepting the trivial certificate instead gives the tautological congruence: χ(r) = 1.
- The orchestrator's literal small seeds produce no rows. Useful seeds need r ≳ N^{1/4}.

**Lemma SC [proved, trivial].** If R, R′ and RR′ are B-smooth integers below N, certified by their own
factorisations, then e(RR′) − e(R) − e(R′) = 0, and under the canonical convention the defect δ is 1.
- So with own-factorisation certificates a fruitful pair needs wrap-around modulo N.
- *Corrected (codex detGRH_2 §4; Opus critic):* the earlier claim that "cocycles generate only over all of G_B²" is
  withdrawn. The defects of an arbitrary section need not generate the kernel: for Z → Z/3 with lifts 0, 4, 2 they
  generate 6Z. Conversely, a few well-chosen ones can generate it.
- What is true is weaker: **no known method locates a short list of pairs containing a fruitful one.**

**Answer to question 2 (structural fruitlessness and Φ).**
- (a) If Φ̄(Λ_c ∩ 2Γ) ≠ 0, testing a basis finds a fruitful dependency deterministically (Lemma H). No seed and no
  GRH are needed.
- (b) Fruitlessness for *every* relation family over the factor base is exactly case (S). There, multiplicativity
  holds and forces failure, in the F₂ route and for Miller bases in G_B. GRH removes (S).
- (c) Fruitlessness of a *given* family in case (G) is a property of that family. It can be total: codex
  detGRH_1 §3.3 notes that the family c_i = r_i², with stored roots r_i, has Φ̄ ≡ 0 for every N.
  - Its exact form is nfsfruit's rowspace criterion: Φ̄ ≡ 0 iff the factor-detecting rows lie in the row space of
    the square-class matrix.
  - Putting a nonresidue into the factor base removes (S). It does not force a given family to use it in a
    fruitful pattern.
  - So the reduction the orchestrator hoped for ("Φ ≢ 0 + GRH + small seed") has this form: GRH removes only (S),
    and Φ̄ ≢ 0 for the actual family remains a hypothesis. §4 replaces it by a Q-rank hypothesis.

## 4. The Z route: fruitfulness under ERH from full Q-rank

**Proposition Z [proved; folklore-level: Miller 1976, Lemma 5 and Theorems 3–4; Bach 1990].**
- *Hypotheses.* Assume ERH and the preprocessing of §1. Let FB ⊇ {−1} ∪ {primes < 2 log² N}, and let Λ_c ⊆ Λ have
  full rank n.
- *Claim.* D := |any nonzero n × n minor of a generating matrix| is a multiple of [Γ : Λ_c], hence of
  |G_B| = φ(N), hence of λ(N).
- *Algorithm.* Miller's test with M = D: for a = 2, …, ⌈2 log² N⌉, first check gcd(a, N); then write D = 2^s t, and
  inspect the last x ≠ 1 in a^t, a^{2t}, …. It returns a proper factor.
- *Cost.* Deterministic, Õ(log² N·(log D + s)) multiplications mod N. D is computed in deterministic polynomial
  time in the matrix's bit size (fraction-free elimination or HNF).

*Proof.*
- **D is a multiple of λ(N).** Γ/Λ_c (of order [Γ : Λ_c], which divides D) surjects onto Γ/Λ ≅ G_B = U.
- **The subgroup H.** Let k be maximal with some b^{2^k t} = −1. Then k ≥ 0, since (−1)^t = −1, and k < s. Put
  H := {a : a^{2^k t} ∈ {±1}}, a subgroup.
- **H is proper.** By CRT take c ≡ b at one prime power and c ≡ 1 at the others; then c^{2^k t} ∉ {±1}.
- **Every a ∉ H splits N.** The last x ≠ 1 in the chain from a^{2^k t} cannot be −1 (by maximality of k, or because
  a ∉ H). So x ∈ μ₂ ∖ {±1}.
- **A small a ∉ H exists.** Bach gives one below 2 log² N; preprocessing makes it a unit. ∎

*Checks.* The Opus critic ran 2 995 random composites with D = λ·random, all split. E1 factored 237/237 instances.

**Proposition EQ [proved; folklore-level reformulation of Miller's ERH equivalence].**
- *Statement.* Assume ERH. Consider all odd N that are not prime powers and have no prime factor < 3 log² N. Then the
  following are deterministically equivalent, with T(N) ↦ O(log N)·T*(N) + poly, where T* is the monotone envelope.
  - (1) splitting N;
  - (2) producing integer vectors of full Q-rank on FB = {−1} ∪ {primes < 2 log² N} that satisfy
    ∏ q^{e_q} ≡ 1 (mod N), encoded in binary.
- *Proof.*
  - (2) ⇒ (1): Proposition Z.
  - (1) ⇒ (2): split recursively, handling powers of 2, primes (AKS), perfect powers (integer roots) and small
    factors (trial division) along the way. This gives the factorisation, hence λ(N). Output λ(N)e_q for each q, and
    2e_{−1}. ∎
- *Scope.* The equivalence concerns **arbitrary** relation sets. An almost-all splitting guarantee would not survive
  the recursion.

**Consequences.**
- **(Z1) [proved].** Given such relations, no F₂ fruitfulness, hidden roots or quadratic characters are needed. This
  is a sufficient architecture. It does not show that any proved relation source supplies the relations.
- **(Z2) [proved; with one interpretive remark].**
  - The orchestrator's "homomorphism χ + small seed lemma" is exactly Miller's argument. The homomorphism is the power
    map a ↦ a^{2^k t} (k need not be known: the algorithm tests every level); the small seed is Miller's base.
  - Z-linear algebra keeps the 2-adic information that F₂ elimination discards.
  - *Remark (interpretation, not a theorem):* this is the "second-order/Bockstein" information of nfsfruit.
  - The doubled control (every relation, including 2e_{−1}, doubled) has Φ̄ ≡ 0, yet the Z route still factors
    (E1 P-Z2).
- **(Z3) [proved for the unrestricted lattice; open for any source].**
  - If w ∈ Q^n vanishes on Λ then w = 0, because Γ/Λ is finite. So there is no rank obstruction of Bockstein type for
    the whole lattice.
  - A restricted source can be rank-deficient for structural reasons. Codex's example, for N = 77: ⟨(2,0,0),
    (0,30,30)⟩ ⊂ Λ uses every column but has rank 2.
  - A missing prime (as in E1) is one detectable defect. Fixing it requires enlarging or changing the relation search
    *and proving* the new source spans.
  - Proving full rank for an explicit source is part of the open hypothesis.
- **(Z4) [proved with qualification].**
  - Usual QS rows 2e_x − e_{Q(x)} need x as a new generator. Exact elimination of these columns kills every such row,
    so the usual QS supply does not feed the Z route.
  - The exception is a row whose x is itself B-smooth (a second smoothness event). For example, 9² − 77 = 4 gives
    3⁴·2^{−2} ≡ 1 (mod 77).
  - The Z route therefore does not rescue detqs's QS fruitfulness (H2).
  - NFS rows can feed it only after exact elimination of the algebraic side (Lemma XK-Z, §5), which is not proved.
- **(Z5) [conditional; UNVERIFIED].** The same idea could replace LV-type fruitfulness conjectures by
  "full projected Q-rank + cost of exact elimination" under GRH. It does not discharge Conjecture 7.1.

## 5. Every source of randomness in Theorem A (question 4)

Theorem A is as in `notes/thmA_proof.md`, with the recursive shape test of `recurse` (Theorem RA). "Unresolved"
means that no deterministic replacement is known. It is never a lower bound.

| # | randomness | role | deterministic replacement | status |
|---|---|---|---|---|
| 1 | (g₀, g₁) uniform in [M/2, M]² | FM rational fibre count; Lemma DS in LD* | none known (§6) | **unresolved (supply)** |
| 2 | w uniform prime in (W, 2W] | LD* via Lemma W (Brun–Titchmarsh spreading of c_P·w mod t) | none known. Finding *a* prime in (W, 2W] deterministically is easy, but LD* uses the average over N^{Θ(1)} values of w | **unresolved (supply)** |
| 3 | c uniform in [−H, H]^d | FM progression F̂ + G·⟨c, v(P)⟩; Kern in LD* | none known (§6) | **unresolved (supply)** |
| 4a | Dixon seed r (uniform), supply role | Theorem S′/K: R = r² uniform on squares hits accepted splits w.p. L^{−o(1)} | none known for explicit seeds. Literal small seeds are rejected (Lemma SR). Replacing them by factor-base products changes the input law, and S′/K do not transfer | **unresolved (supply)** |
| 4b | Dixon seed r, fruitfulness role (Theorem DD) | the hidden root on its μ₂-fibre | **under ERH, replaced by the Z route**: seeds x_e ∈ ⟨FB⟩, descents, exact elimination (Lemma XK-Z), Proposition Z | **conditional** on XK-Z and on full projected rank |
| 5 | special-q search | — | already deterministic in (A, ℓ) | none needed |
| 6 | Lenstra–Pomerance shape test | factor \|u₀\|, \|v₀\| ≤ √(2N) | recursion is deterministic *if* a deterministic outer algorithm with uniform cost and supply guarantees exists (strong induction; small inputs by trial division). Non-recursive deterministic: Harvey N^{1/10+o(1)} per half | conditional closure property only |
| 7 | LV6 kernel vectors (Wiedemann, random) | F₂ dependencies | deterministic F₂ elimination: L^{3a+o(1)} (cubic) or L^{ωa+o(1)}. Testing a basis suffices (Lemma H) | **removable** for the F₂ route |
| 7′ | exact Z-linear algebra (Z route) | rational rank, minors, Dixon lifting | dlog Prop XK step 2 picks a **random** prime q (Las Vegas). Deterministic replacement: a rank profile and a nonzero minor modulo each of L^{a+o(1)} word-size primes whose product exceeds the Hadamard bound, then deterministic p-adic lifting / fraction-free elimination: L^{(ω+1)a+o(1)} (Storjohann-type bounds) | removable at cost L^{(ω+1)a}; see the slack note |
| 8 | quadratic-character primes (LV6 Lemma X) | cut K^×/K^{×2} | under GRH for the Kummer extensions K(√h)/K, effective Chebotarev (Lagarias–Montgomery–Odlyzko 1979; Bach 1990) gives separating first-degree primes. Since log disc K(√h) = L^{a+o(1)}, about L^{2a+o(1)} primes are needed. They must also avoid the supports of the atoms, so the exclusion set enters the bound; character evaluation costs ≈ L^{3a}. Or avoid characters entirely via the Z route | removable under GRH (with this bookkeeping still to be written out) |
| 9 | Hensel prime, square root in Z[β] | exact square root | an inert prime need not exist (the Galois group may lack a d-cycle). Use the least prime splitting f completely: polylogarithmic under GRH-Chebotarev for the fixed-degree splitting field. Exhaustive roots, deterministic square roots in F_p (nonresidue search costs O(p) at p ≍ B), 2^d sign patterns verified over Z | removable at fixed degree |
| 10 | roots of f mod r ≤ B (columns) | — | exhaustive O(rd), total L^{2a} | removable |
| 11 | restarts after a trivial gcd | geometric number of LV6 rounds | the Z route has no fruitfulness retries (Miller is deterministic under ERH). Retries for relation supply remain | removable under ERH |
| 12 | Lemma IR certificate, gcd factor exits, content checks | — | deterministic | none needed |
| 13 | computable S, parameters, AKS, Gauss reduction | — | deterministic | none needed |

**Slack for the deterministic exact linear algebra (Opus critic M1; codex detGRH_2 §8) [proved arithmetic].**
- **In Theorem RA.** a\*/c\* = 1/(2 + √2) = 1 − 1/√2 = 0.2929.
  - A deterministic cost L^{Ca} fits iff C < 2 + √2.
  - With C = ω + 1 this needs ω < 1 + √2 ≈ 2.414. With ω = 2.3714 the ratio is 0.987, so it fits.
  - With cubic elimination (C = 4) the ratio is 1.17, so it does **not** fit.
- **In Theorem A's own design.** a = θ/√2, and (ω+1)a < (1+θ)/√2 iff θ < 1/ω ≈ 0.4217.
  - This holds in the regime c < 1, θ < √2 − 1.
  - It fails for the larger θ < 1 that Theorem A allows. The earlier "which holds" was wrong for that range.
- **The F₂ route** costs L^{3a} for cubic elimination; it is slack in RA (0.879) and needs θ < 1/2 in Theorem A.

**Lemma XK-Z (needed, NOT proved) [open].**
- *Statement.* Given M NFS descent rows over a fixed field (rational columns, ideal columns, the g₁ column, and
  seed-exponent columns), deterministic exact elimination of the ideal and unit parts yields pure rational kernel
  relations ∏_{q∈FB} q^{e_q} ≡ ±1 (mod N). Their projection to Z^{FB} has rank equal to the rank of the rows'
  rational parts modulo the algebraic dependencies, in time L^{(ω+1)a+o(1)}, with output bit sizes L^{a+o(1)}.
- *What exists.* dlog Prop XK eliminates *all* columns and outputs one primitive torsion dependency, in the
  discrete-log setting, using a random prime.
- *What is missing.*
  - Deterministic multimodular rank (as in row 7′).
  - The unit-elimination step (Lemma XK4) applied to a basis of the algebraic-side kernel, not to d vectors.
  - A bit-size and cost ledger for many output relations.
- We expect it to go through at fixed d, but have not written it out.

**Hypothesis DS-NFS(c₀) [open].** There is a deterministic procedure that, given N, runs in total time L_N[1/2, c₀].
Total time includes listing the (field, seed) pairs, executing all descents and shape tests, and the exact
elimination of Lemma XK-Z. Its output is rational kernel relations of full Q-rank on FB ⊇ {primes < 2 log² N}.

**Corollary DZ [the deduction is proved].** ERH + DS-NFS(c₀) ⇒ deterministic splitting in time
L_N[1/2, c₀ + o(1)] + poly(output). It gives complete factorisation of every N with ≤ log₂ N calls on cofactors,
assuming DS-NFS for all N.
- *Proof.* Proposition Z on the output; its cost is polynomial in the output's bit size, which DS-NFS's time bound
  covers.
- *Status.* We do **not** claim that DS-NFS is equivalent to factoring. Proposition EQ's converse produces λ(N)e_q,
  which are not NFS descents (Opus critic M2; codex detGRH_2 §6). DS-NFS is a sufficient, explicitly checkable
  target.

## 6. Derandomizing the field (question 3) and the seeds

**Checkability alone gives nothing [proved, trivial].** For any explicit list E of m fields in a family 𝒜, and any
γ ≤ 1 − m/|𝒜|, some subset of density γ avoids E. So global moments are compatible with total failure on E: codex
detGRH_1 §5.1 gives an exact example with (E Z)²/E Z² = γ. Any proof that E contains a good field must use the
arithmetic of good fields.

**Observation FS (a limitation of the present full-progression argument; not a lower bound) [observation].**
- **How FM uses c.** For fixed P, as c varies, F_A(P) runs through F̂(P) + G(P)·j with |j| ≤ dHV^{d−1} and modulus
  t = |G(P)| = N^{μ+ν+o(1)}.
- **The analytic input.** FM's input (LV Lemmas 5.13–5.20 plus Lemma 8.5′) handles smooth numbers in *full-length*
  progressions, averaged over the moduli t at level t ≤ X_F^{1/2}. That uses X_F/t ≥ t^{1+ε} terms.
- **Short boxes.** A short coefficient box realises only J terms. If J ≪ X_F^{1/2}/t, the needed statement is about
  smooth numbers in short progressions (intervals of length tJ ≤ X_F^{1/2}), which no known theorem provides at
  B = L^{O(1)}.
  - The earlier "J ≪ t ⇒ tJ ≪ X_F^{1/2+o(1)}" was arithmetically wrong (codex detGRH_2 §9). The correct condition is
    J ≪ X_F^{1/2}/t. A subexponential J satisfies it.
- **Lemma W and Lemma DS** use N^{Θ(1)}-size averages over w and (g₀, g₁) in the same way.
- **What this does not exclude.**
  - A sparse explicit sub-family spread over the full range, which is not a short box.
  - Another proof technique.
- **Individual-progression theorems.** The best known results for smooth numbers in a *single* progression are
  unconditional: Granville 1993, Soundararajan 2008, Harper 2012. They need q ≤ y^{4√e−δ} and log x/log q → ∞, and
  both fail here, since q = N^{Θ(1)} ≫ y^{A} for every A when y = L^{O(1)}.
  - We found no GRH-conditional count of smooth numbers in an individual progression reaching q = N^{Θ(1)} (Opus
    literature critic). The GRH results we found concern character sums over smooth numbers.

**Observation CS (Fourier certification of explicit seeds) [observation; conditional computation].**
- **Setup.**
  - E is a family of m **distinct** residues x_e. For a multiset, replace m by the collision energy, i.e. m^{−1/2}
    by √(𝒞_E)/m.
  - S ⊂ [1, √N] is the target set of descendable halves, with s := |S| = N^{1/2}L^{−o(1)}.
  - Z := #{(e, u, v) : x_e ≡ u/v} = φ(N)^{−1} Σ_χ Ê(χ)|Ŝ(χ)|², with main term ms²/φ(N).
- **Inputs.**
  - Parseval: Σ_χ|Ê(χ)|² = φ(N)·m (equality, for distinct residues).
  - U_S := max_{χ≠χ₀}|Ŝ(χ)|.
- **Hölder.** With conjugate exponents p ≤ 2 ≤ q,
  error/main ≤ (φ(N)/√m)·s^{−2+1/q}·U_S^{2−2/q}.
  - (i) Under the **assumption** U_S ≤ s·N^{−1/4+o(1)} (square-root cancellation), this is ≤ m^{−1/2}N^{1/2+o(1)} for
    every q ≥ 2. The bound certifies a hit only if m ≥ N^{1+o(1)}.
  - (ii) With only U_S ≤ s, the best choice q = 2 gives m^{−1/2}N^{3/4}, so it needs m ≥ N^{3/2}.
  - (iii) p > 2 is worse.
  - Both thresholds exceed the trivial covering size φ(N).
- **What changed after review** (Opus critic M3; codex detGRH_2 §9).
  - The square-root cancellation in (i) is **not** a GRH consequence for S. GRH bounds character sums over *all*
    primes of J, not over a field-dependent good subset or its products.
  - A character-biased S can have Ŝ(χ) = s, e.g. S = {a ≤ √N : χ(a) = 1}.
  - The conclusion is only that *this absolute-value certificate* cannot certify small lists. It is not a necessary
    seed count. E = {1} trivially meets S/S (u = v), though uselessly.
- **Expander walks.** Under GRH the Cayley graph of U with generators the primes ≤ (log N)^{2+δ} is an expander
  (Jao–Miller–Venkatesan, JNT 129 (2009), Thm 1.1). Random walks of length O(log N/log log N) mix (Cor. 1.3). That is
  a statement about random sampling; enumerating all such walks costs N^{2+δ+o(1)}. This is no lower bound for
  cleverer deterministic uses.

**Q3 verdict [assessment].**
- No known method, conditional on GRH or not, proves that an explicit L^{O(1)} list of fields or seeds contains a
  good one.
- The short-box variant of FM/LD* is outside the present analytic input (FS). A second moment over an explicit
  initial segment is a short box only for consecutive enumerations; for other segments it needs new list-specific
  estimates.
- The only known derandomization route is the Impagliazzo–Wigderson hypothesis (polyimpl Theorem DR: deterministic
  L[1/2, O(1)] via Lenstra–Pomerance). It is not known to follow from GRH.

## 7. Experiment E1 (pre-registered, `7cab676`; toy sizes, correctness only)

All predictions were frozen in `notes/detGRH_prereg.md` before any code existed.

**Disclosed deviations.**
- (1) The tri class cannot have all three factors > B(28) = 768 at 28 bits, since the smallest factor is
  ≤ N^{1/3} < 2^{9.4}. tri/28 was therefore replaced by tri/40, in the config before the run. Every tri prime has
  ⌈n/3⌉ bits.
- (2) The blind prompt pointed to Miller/Bach in question 2(c) (disclosed in the prereg).
- No other deviations.

**Results (240 instances + exhaustive E1a; `experiments/analysis/detGRH_e1_summary.json`).**

| prediction | result | verdict |
|---|---|---|
| P-MR1: #sections = [all p ≡ 3 mod 4], every odd N < 3000 | 1499 N, 0 violations; an independent backtracking count agrees for N < 160 (test) | **pass** |
| P-MR2: the principal root is multiplicative (Blum-type N < 400) | 85 N, 0 violations | **pass** |
| primes < 2 log² N generate (Z/N)^× (diagnostic) | 240/240 | as predicted |
| full Q-rank within the caps (descriptive) | 237/240; the 3 failures are tri/40 at the cap K = 600 | — |
| P-Z1: full rank ⇒ Miller factors N | **237/237**; base 2 in 165 cases, 3 in 48, 5 in 14, at most 13 | **pass** |
| P-Z2: control 2Λ_c, F₂ route dead and Z route alive | F₂ fruitless 240/240; Z route factors 237/237 | **pass** |
| P-FL1: dim W = ω − 1 | 240/240 | **pass** |
| P-FL2: c < dim W ⇒ fruitful | premise true in 42 instances (all tri, c = 1), 0 violations | **pass** |
| P-CY1: δ ∈ μ₂; δ ∉ {±1} ⇒ proper gcd | 0 violations; 61 734 nontrivial defects, all giving proper factors | **pass** |
| P-CY2: nontrivial-defect rate, semiprimes [0.30, 0.70] | 52 401/104 800 = **0.500** | **pass** |
| P-CY2: tri [0.55, 0.95] | 9 333/12 506 = **0.746** | **pass** |
| P-CY3: non-diagonal χ(r_i), semiprimes [0.30, 0.70] | 6 227/25 992 = **0.240** | **FAIL** |
| P-CY3: tri [0.55, 0.95] | 1 344/5 747 = **0.234** | **FAIL** |

**What E1 shows.**
- **Exact algebra.** MR, FL, Z and CY are confirmed with no exceptions. These are correctness checks of proved
  statements.
- **The orchestrator's multiplicativity hypothesis is refuted on every class.** Products of certified squares
  re-certified independently disagree with the product of the roots in exactly the fair-coin proportion,
  1 − 2^{1−ω}. Every disagreement is a factor, as Theorem M(G) predicts.
- **The Z route works with no F₂ input**, including the doubled control, where every F₂ dependency is fruitless by
  construction.

**Post-hoc analysis (labelled; `src/detGRH/posthoc_e1.py`; 16 instances at 32 bits).**
- **Why P-CY3 failed.** For 1 211/2 388 seeds, the canonical completion of e(R_i) is the seed's own relation
  2e(r_i) − e(R_i). Then y(R_i) = r_i exactly and χ(r_i) = 1 by construction: an honest certificate, as in Lemma SC.
  Among the other 1 177 seeds the non-diagonal rate is 0.538.
- The pre-registered window wrongly assumed that the completion is a uniform dependency. The failure is in our
  heuristic prediction, not in a claim.
- It does illustrate Theorem M(G): χ at the seeds is mostly trivial, and the fruitful information sits in the
  cocycles.
- **Why c ≥ 1 in every instance (c_hist: 1 → 149, 2 → 49, 3 → 21, 4 → 12, 5 → 4, 6 → 1, 21 → 1).**
  - Both relation sources have positive values, so no found relation has an odd −1 coordinate (0/11 499 rows).
    Hence Λ_c lies in the index-2 sublattice {e ∈ Λ : e_{−1} even}, and Λ does contain elements with e_{−1} odd,
    since −1 ∈ ⟨primes⟩.
  - So [Λ : Λ_c] is always even. Proposition FL(iii)'s sufficient condition then never applies to semiprimes, yet the
    F₂ route was fruitful in 240/240 instances (195/195 with c ≥ dim W).
  - The criterion is sufficient, not necessary, as stated.
- **The rank shortfall.** It is a pinning failure: for example one prime (971, in a 32-bit trial run with a small
  cap) occurred in no relation. This is the detectable degeneracy of Z3.
- **No complexity claim.** Toy sizes, rational relation sources with Dixon-size values.


## 8. Literature (question 5)

Sources: my own searches, codex detGRH_1 §7, codex detGRH_2 §10 and the Opus literature critic. Primary sources
were checked where reachable. Exact theorem numbers for Bach 1990 were not verified.

- **Deterministic factoring.**
  - Under GRH: N^{1/5+ε}, by Shanks, "Class number, a theory of factorization, and genera", Proc. Sympos. Pure
    Math. 20 (1971) 415–440, and Schoof, "Quadratic fields and factorization", Computational Methods in Number
    Theory, MC Tract 154/155 (1982) 235–286, §6 (Algorithms 6.4, 6.5; class groups and infrastructure). Verified in
    the scan by the critic.
  - Unconditional: Harvey, Math. Comp. 90 (2021) 2937–2950; Harvey–Hittmeir, Math. Comp. 91 (2022) 1367–1379.
  - Newer preprints: Gao–Feng–Hu–Pan, arXiv 2512.19076; Harvey–Hittmeir, arXiv 2601.11131. Umans–Wang (arXiv
    2511.10851) still call 1/5 the current exponent.
  - Classical: Pollard–Strassen, N^{1/4}.
- **No deterministic subexponential complete factorisation under GRH was found**, for all N or for a density-one set.
  - Oracle results are not counterexamples. Żrałek, Math. Comp. 79 (2010) 513–533, Thm 8.1/Cor 8.2, gives
    deterministic L[1/3] *given the iterates of φ*. See also Morain–Renault–Smith (arXiv 1802.08444) and
    Dąbrowski–Pomykała–Shparlinski (oracle, almost all n).
  - Finding one factor of almost all integers is trivial (trial division to log N). The open question concerns
    complete factorisation or a hard ensemble.
- **Rigorous subexponential methods are all randomized.**
  - Dixon 1981.
  - Schnorr 1982 (J. Algorithms 3), heuristic and randomized.
  - Schnorr–Lenstra 1984, Monte Carlo.
  - Pomerance 1987: L[1/2, √2], randomized, also for discrete logs in F_p.
  - Seysen 1987 (Math. Comp. 48): probabilistic, ERH.
  - A. K. Lenstra 1988 (Indag. Math. 50): probabilistic, GRH.
  - Hafner–McCurley 1989 (JAMS 2): Las Vegas class groups, ERH.
  - Vallée 1991; Lenstra–Pomerance 1992 (JAMS 5).
  - Buchmann 1990: heuristic.
  - de Boer–Pellet-Mary–Wesolowski (arXiv 2512.01588): rigorous under ERH, but probabilistic.
  - Lee–Venkatesan 2018 (JNT 187): a randomized NFS; fruitfulness assumes their Conjecture 7.1, which is not a
    consequence of GRH.
  - We found no deterministic subexponential class-group or F_p discrete-log algorithm under GRH. The deterministic
    GRH bound for class groups is |D|^{1/5+ε} (Lenstra 1982; Schoof 1982).
- **Miller and Bach.**
  - Miller, JCSS 13 (1976) 300–317, Lemma 5 and Theorems 3–4: under ERH, factoring reduces deterministically to
    suitable exponent functions (φ, λ, or a multiple of λ). See Żrałek Thm 6.2: given φ(n), deterministic O(log⁶ n)
    under ERH. See also Shoup §§10.4–10.5.
  - Bach, *Discrete Logarithms and Factoring*, UCB/CSD-84-186 (1984), §3: factoring from exponents obtained via
    discrete logs. It notes the reduction becomes deterministic under ERH. (Codex detGRH_2 cites §3; the literature
    critic judged from the abstract that it does not contain the relation-lattice statement. We cite it only for the
    exponent-to-factor step.)
  - Bach, Math. Comp. 55 (1990): the constants 2 log² m (any element) and 3 log² m (prime ∤ m). Unconditional
    deterministic factoring from φ is open (Adleman–McCurley problem 23).
- **Index calculus mod N.**
  - McCurley, J. Cryptology 1 (1988) 95–105, §3.1: order multiples from exponent collisions, then Miller.
  - McCurley, "The discrete logarithm problem", Proc. Symp. Appl. Math. 42 (1990).
  - **Stange, "Factoring using multiplicative relations modulo n" (arXiv 2211.06821), §2 and Thm 3.2.** The closest
    prior formulation of the Z route: multiplicative relations mod n, rational linear algebra, an order multiple.
    Her version is heuristic in relation supply and randomized.
  - The class-group analogue (relation lattice of small prime forms → h → ambiguous forms) is Seysen, Hafner–McCurley
    and Lenstra–Pomerance.
  - Proposition EQ in its stated form was not found by any reviewer. It is a folklore-level reformulation, and we
    make no novelty claim.
- **Square roots.**
  - Rabin 1979: randomized reduction.
  - Jeřábek, "Integer factoring and modular square roots", JCSS 82 (2016) 380–394, Thm 3.8: under quadratic GRH,
    factoring ≡ FacRoot (a root or a factor for every a with (a/N) = 1) by deterministic reductions.
  - For an oracle defined only on squares and answering adversarially, no deterministic reduction is known.
    Jeřábek notes that a "devious" oracle can return the root the caller already knows.
  - Lemma MR: no published statement found. It is elementary.
- **Expanders.** Jao–Miller–Venkatesan, JNT 129 (2009) 1491–1504, Thm 1.1 and Cor 1.3.
- **Smooth numbers in a single progression.** All unconditional: Granville 1993 (q^{1+ε} ≤ y); Soundararajan 2008
  (q ≤ y^{4√e−δ}); Harper 2012 (JNT 132). Averaged results (Fouvry–Tenenbaum 1996; Drappeau; Pascadi, arXiv
  2304.11696) cannot single out a modulus. Nothing reaches q = N^{Θ(1)} with y = L^{O(1)}.
- **Derandomization.**
  - Impagliazzo–Wigderson 1997 (and Umans 2003) give deterministic L[1/2, O(1)] via Lenstra–Pomerance (polyimpl
    Theorem DR).
  - Pseudodeterministic factoring is trivial, since the output is unique (Gat–Goldwasser).
  - No factoring-specific barrier to derandomization was found. Kabanets–Impagliazzo 2004 is a general obstacle to
    unconditional derandomization.

## 9. Reviews

**GPT-6.1-sol blind round (`codex/detGRH_1`; report `notes/detGRH_codex1_report.md`).**
- It reached, independently, the same main conclusions:
  - multiplicative selectors are impossible at p ≡ 1 (mod 4) and forced to be the principal root otherwise;
  - the Z route via full-rank relations and Miller–Bach;
  - field moments do not give a short list;
  - the same final hypothesis.
- Disclosure: question 2(c) of its prompt pointed to Miller and Bach.
- Merged from it: Lemma SR, the order-4 proof of MR, family-level degeneracy (c_i = r_i²), Jeřábek 2016, the
  "almost all" triviality, the rewording of "irreducible", and the deterministic character-separation argument.

**GPT-6.1-sol adversarial referee (`codex/detGRH_2`; report `notes/detGRH_codex2_report.md`).**
- Verdicts:
  - CORRECT: FL.
  - CORRECT WITH FIXES: MR, CD, Z, EQ.
  - GAP: M as applied to Theorem A; Z1–Z5; the §5 table; the overall assessment.
  - WRONG as stated: FS and CS.
- Every fix is applied:
  - the canonical convention and root-value honesty in CD and M(G), with the N = 77 counterexample and the
    non-vacuity example N = 1019·1031;
  - Theorem M is not applicable to Theorem A's partial descent;
  - the Miller-in-G_B explanation (N = 65);
  - the Schreier claim withdrawn (Z → Z/3 example);
  - Φ̄ in place of Φ;
  - the Bach unit/nonunit constants and preprocessing;
  - EQ restricted to arbitrary relations, with DS-NFS ⇒ factoring only;
  - Z3 scoped (rank-2 example), Z4 qualified (9² − 77);
  - XK's random prime, the deterministic L^{(ω+1)a} cost and the θ range;
  - row 8 bookkeeping and row 9 corrected;
  - FS relabelled and its J-arithmetic fixed;
  - CS made conditional, with collisions handled;
  - Stange 2022 and Shoup added.

**Opus proof critic (read-only).**
- No FATAL issue. Three MAJOR issues, all applied:
  - (M1) exact linear algebra is not deterministic as cited; the deterministic cost is L^{(ω+1)a}, which needs
    ω < 1 + √2 in RA and θ < 1/ω in Theorem A;
  - (M2) the DS-NFS equivalence was an overclaim;
  - (M3) CS's sup bound is not a GRH consequence, and the Parseval direction was wrong.
- Minor fixes applied:
  - the MR(a) automorphism wording;
  - CD's "homomorphism of F₂-type data" sentence (replaced by the canonical-convention homomorphism);
  - the Miller explanation in (S);
  - "B1 is an instance", not "exactly";
  - FL "exact" versus "sufficient";
  - row 8 count (L^{2a} primes);
  - row 9 (inert prime need not exist);
  - "known only under IW".
- Numerical checks: MR on 30 further N, including prime powers; FL on 600 random lattices; Z on 2 995 composites;
  B1 on 3 moduli.

**Opus literature critic (read-only).**
- Applied:
  - credit Shanks 1971 with Schoof 1982 §6 (verified);
  - Żrałek 2010 as an oracle result;
  - Bach's constants (2 log² any element; 3 log² prime ∤ m) and LLS 2015;
  - the citation for index calculus mod N (Bach 1984 kept only for the exponent step);
  - individual-progression smooth-number results are unconditional, with q ≤ y^{4√e−δ};
  - GRH √Y cancellation is for primes;
  - "not known to follow from GRH";
  - JMV mixing time;
  - Lagarias–Montgomery–Odlyzko for effective Chebotarev;
  - Jeřábek; Burthe/Bellotti (below).
- **Noted, not used.** For almost all N, the ERH hypothesis in Proposition Z can be replaced by unconditional bounds
  on the least generating set:
  - Burthe, Acta Arith. 80 (1997): G(n) ≤ (log n)^{97+ε} for almost all n;
  - Bellotti, arXiv 2610.10177 (Oct 2026, unrefereed): (log n)².
- This does not help, because supply is the open part.

**Surviving claims.**
- [proved]: Lemma MR, Corollary MR (GRH), Proposition FL, Proposition CD (with the canonical convention), Theorem M
  (abstract), Lemma SR, Lemma SC.
- [proved, folklore-level]: Propositions Z and EQ (ERH), and the deduction ERH + DS-NFS ⇒ factoring.
- [observations]: FS, CS (conditional).
- [open]: Lemma XK-Z, DS-NFS.
- [experimentally supported]: E1, toy correctness only.

## 10. Candidate cards

```text
ID: detGRH-1 multiplicative square-root selector + Ankeny/Bach small seed (orchestrator's idea).
Target: every N, deterministic, GRH.
Operation: chi(r) = y(r^2)/r homomorphic => small r with non-diagonal chi.
Proved: implication valid (Cor MR); hypothesis possible only for all p = 3 mod 4, and then y = principal root
  (Lemma MR); certificate roots are not functions of residues when G_B cap mu_2 > {+-1} (Prop CD);
  honest total root functions are non-multiplicative, so some pair is fruitful (Thm M(G)), but nothing locates it;
  Theorem A rejects small seeds (Lemma SR).
Status: refuted for this model (as a route to deterministic fruitfulness).

ID: detGRH-2 Z route (full-rank kernel relations + Miller/Bach).
Target: every N, deterministic, ERH.
Proved: full Q-rank on FB containing primes < 2 log^2 N => factor (Prop Z); equivalent to factoring for arbitrary
  relation sets (Prop EQ); no F2 fruitfulness, characters or hidden roots needed.
Missing: deterministic relation supply with full rank (DS-NFS) + exact elimination lemma XK-Z.
E1: 237/237 full-rank toy instances factored, including the F2-dead doubled control.
Status: proved reduction (folklore-level); supply OPEN.

ID: detGRH-3 short explicit field list via FM/LD* / second moment / GRH tools.
Status: open; the present full-progression argument needs N^{Omega(1)} families (Observation FS, a limitation of
  that proof only); Fourier certification of seed lists needs N^{1+o(1)} distinct seeds even under a square-root
  cancellation assumption not known under GRH (Observation CS).
```

## 11. Reproduction

```
PYTHONPATH=. .venv/bin/python -m src.detGRH.harness_e1 experiments/configs/detGRH_e1.json experiments/raw/detGRH_e1.jsonl
PYTHONPATH=. .venv/bin/python -m src.detGRH.analyze experiments/raw/detGRH_e1.jsonl experiments/analysis/detGRH_e1_summary.json
PYTHONPATH=. .venv/bin/python src/detGRH/posthoc_e1.py      # post hoc, labelled
.venv/bin/python -m pytest -q tests/test_detGRH.py
```
The E1 run takes about 13 minutes on the development machine.
