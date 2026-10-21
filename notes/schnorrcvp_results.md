# schnorrcvp: factoring via CVP/SVP in prime-number lattices (Schnorr 1991/2021, Adleman 1995)

Agent `schnorrcvp`, phase 19, 2026-10-21.
- Code: `src/schnorrcvp/`:
  - `core.py`: exact relation objects, Λ_N, GF(2) square roots;
  - `psi.py`: exact Ψ brackets and the counting models;
  - `lattice.py`: fpylll lattices and exact list-CVP;
  - `enum_rel.py`: exact box enumeration;
  - `run_e0*.py`, `run_e1.py`, `run_e2.py`, `analyze.py`.
- Tests: `tests/test_schnorrcvp.py` (11 pass).
- Pre-registration: `notes/schnorrcvp_prereg.md`, commit `ecd3c25`, made before E1/E2. Raw data:
  `experiments/raw/schnorrcvp_E*.jsonl`.
- GPT-6.1-sol: blind round `codex/schnorrcvp_1.{md,out}` (§9.1); adversarial round `codex/schnorrcvp_2.{md,out}`
  (§9.2).
- Reviews: one read-only Opus proof critic and one read-only Opus literature critic (§9.3).

Labels:
- **PROVED**: proof here, or a cited theorem checked in the primary source.
- **KNOWN**: published, with a citation.
- **HEURISTIC (name)**: rests on the named assumption.
- **MODEL**: an exact consequence of named heuristics, computed in code.
- **MEASURED**: our runs.
- **POST-HOC**: chosen after seeing data.
- **OPEN**.

"No prior statement found by us" is the strongest novelty claim made in this note.

---

## 0. Verdict

**Plain language.**
- Factoring does reduce to CVP in prime-number lattices, but only for approximation factors **γ → 1**:
  - γ = 1 + n^{−ε} gives heuristic polynomial time;
  - γ = 1 + (log n)^{−c} gives L[1/(2+c)];
  - γ = 1 + Θ(1) gives L[1/2].
  - These assume a polynomial-time oracle on this lattice family (heuristic; §2.4–2.5).
- Every factor that lattice algorithms or lattice cryptography care about (poly(n), 2^{n^ε}, 2^{n/log n}) gives *nothing* through this route. A γ-approximate answer can have a residual up to N^{Θ(γ²−1)} times larger (heuristically tight, Prop A).
- In γ, these statements are weaker than what NP-hardness of CVP already gives rigorously for *general* lattices. The only new content is that the oracle is needed only on an explicit lattice family.
- So **no high-profile link** between factoring and the lattice problems underlying post-quantum cryptography is known in either direction, and this route heuristically cannot supply one.
- Schnorr 2021 fails before any lattice algorithm is involved. At his parameters (n = 47 at 400 bits, n = 95 at 800 bits), heuristically **no fac-relations exist at all**: the expected total number of all sizes is 10^{−52} and 10^{−115}, and 48 or 96 are needed.
- The dimension at which enough relations exist is n ≈ 275 at 400 bits and ≈ 775 at 800 bits. It grows like (ln N)^{3/2+o(1)}.
- With any CVP algorithm costing 2^{Θ(n)}, the approach is asymptotically **slower than trial division**. Our concrete cost model gives ≈ 2^{187} at 400 bits and ≈ 2^{467} at 800 bits, against roughly 2^{57} and 2^{80} for QS and GNFS. At 800 bits that is above trial division (2^{400}).

**Precisely.**

1. **Theorem R (distance–residual identity) [PROVED].** For the lattice with diagonal √ℓ_i (ℓ_i = ln p_i) and
   scale C:
   - ‖b(e) − t‖² = Q(e) + C²Δ²;
   - u − vN = 2e^{(A+L)/2} sinh(Δ/2);
   - A = ln(uv) ≤ Q(e).

   Moreover gcd(uv, w) = 1, and |uw| ≥ N − 1. The last is resbarrier's bound, and it is exactly λ₁(Λ_N) ≥ ln(N − 1) for the
   relation lattice Λ_N in the norm Σ|a_i|ℓ_i.
2. **Theorem P (promise reduction) [PROVED; a quantitative CVP form of Adleman 1995 §2 and Schnorr 1991 §3].**
   Suppose some squarefree relation (witness) e⁰ has |w₀| ≤ ω. Then *every* γ-CVP answer, at scale C = e^{(L+A₀)/2}/ω, satisfies
   |w| ≤ ω·γ√(A₀+1)·exp(((γ²−1)A₀ + γ²)/2 + o(1)). For γ = 1 and ω = B·e^{−0.7}/√(A₀+1), every exact-CVP answer is a
   fac-relation, even when A₀ is unknown and C is scanned (Corollary P1).
3. **Lemma LA [PROVED].** Suppose the relations generate a sublattice of odd index in Λ_N and ⟨P⟩ mod N contains a
   square root of 1 other than ±1. That holds under GRH when B ≥ 2 ln²N (Bach), for N not a prime power. Then
   testing a GF(2)-kernel *basis* factors N deterministically.
4. **Theorem H (heuristic reduction) [HEURISTIC: H_S witness existence, H_LA odd index; GRH; N not a prime power].**
   - Setting: B = (ln N)^κ with κ > 2, n = π(B), targets k·N for primes k ∈ (B, 4B].
   - Then n + 2 relations come from O(nL) calls to an exact-CVP oracle, or to a γ-CVP oracle with
     γ − 1 ≤ c_κ·ln B/L, where c_κ = (1 − 3/(2κ))/(1 + 2a) + o(1) and v ≈ N^a.
   - Prop A says that larger constants fail.
   - Relations from distinct primes k are automatically distinct [PROVED], and they factor N.
   - In terms of n: factoring ≤ CVP_{1+n^{−ε}} on prime lattices of dimension (log N)^{κ}, for any fixed
     κ > max(2, 1/ε).
5. **Theorem G (cost as a function of γ) [MODEL: H_res′ for a "benign" oracle].**
   - γ = 1 + τ·L^{2α−1}(ln L)^{1−2α} gives L_N[α, 2√(ατ)] for 0 < α ≤ 1/2, with n = L[α] and a *unit-cost*
     (n^{o(1)}) oracle.
   - With oracle cost n^{c} the constant is 2√((1+c)ατ). The exponents α are unchanged.
   - Hence:
     - γ = 1 + n^{−ε}: polynomial;
     - γ = 1 + (log n)^{−c}: L[1/(2+c)], so α < 1/3 needs c > 1;
     - γ − 1 = Θ(1): L[1/2, √(γ²−1)] (unit-cost oracle), which beats QS iff γ < √2;
     - γ ≥ 2, poly(n), 2^{n^ε}, 2^{n/log n}: no gain at all.
   - **Prop A [HEURISTIC H_S].** The bound in Theorem P is tight up to O(ln L) in the exponent. So an adversarial
     oracle defeats the reduction once (γ² − 1)A₀ ≫ ln B.
6. **Prop NP (subsumption) [PROVED from known theorems].**
   - Factoring ≤_T GapCVP_γ for every γ ≤ n^{c/log log n}, on *general* lattices of dimension poly(log N),
     because factoring ∈ NP and GapCVP is NP-hard there (Dinur–Kindler–Raz–Safra 2003).
   - Similarly factoring ≤_R GapSVP_γ for every constant γ (Khot 2005; Haviv–Regev 2007).
   - So Adleman's 1995 heuristic "factoring ≤_R SVP" has been a theorem since Ajtai 1998. Ajtai's proof itself used the
     Schnorr–Adleman lattice.
   - The prime-lattice route reaches only γ = 1 + o(1), deep inside the NP-hard regime.
   - No reduction is known for γ ≥ √(n/log n), where GapCVP_γ ∈ NP ∩ coAM (Goldreich–Goldwasser 2000).
     Heuristically (Prop A), this route's query strategy cannot supply one.
   - Faster factoring is not known to help any lattice problem.
7. **Schnorr 2021 dissection.**
   - **(a) No relations exist [MODEL, H_S].** E_tot(400, 47) = 10^{−51.9}; E_tot(800, 95) = 10^{−114.8}. His
     claimed box (v ≲ p_n, |w| ≲ 523) holds 10^{−176} expected relations at 800 bits. The 499-smooth integers below
     2^{800} number 1.9·10^{62}, a density of 10^{−180}.
   - **(b) Lemma 3.1 drops a factor N^{n/(n+1)} [PROVED, primary text checked].** Codex found this first.
   - **(c)** The 800-bit "‖b₁‖ < 0.917" is impossible, since every nonzero lattice vector has norm ≥ 1 [PROVED].
   - **(d) Real cost [MODEL].**
     - Relations exist only for n ≥ n_min(N) = 275 / 400 / 775 / 1125 at 400 / 512 / 800 / 1024 bits, growing as
       (ln N)^{3/2+o(1)} (the local slope measured is 1.56).
     - With BKZ-β, γ_β = δ_β^{Θ(n)} ≫ 2 unless β ≈ n. So the cost is ≥ 2^{0.292 n_min}, which is
       exp((ln N)^{3/2+o(1)}).
8. **Experiments (pre-registered; all pre-registered predictions passed).**
   - E1: exact box counts of primitive relations agree with H_S′ within ×1.63 pooled. The POST-HOC parity/CRT
     model brings this to 0.92.
   - E2: ≈ 79,000 CVP outputs (78,963) at 20–40 bits. Relation yield vs H_res′: pooled 0.87. Schnorr-2021 scaling leaves
     |w| ≈ N^{0.92–0.96}. Every N of 20–36 bits was factored end-to-end as a correctness check.

**Best open question.** Is there any family of lattices, or any norm, in which a *non-NP-hard* approximation factor
(γ ≥ √n) controls the residual of a factoring relation?
- Heuristically, the exponential sensitivity (Theorem P with Prop A) defeats the Schnorr/Adleman query strategy, in
  which the residual's log-size is linear in the lattice norm.
- It says nothing about reductions that query other lattices or embeddings.
- A positive answer would need a residual that depends only polynomially on the norm.

---

## 1. Setting and exact identities

Notation:
- N is odd and has no prime factor ≤ B (trial-divide first).
- P = {p₁, …, p_n} are the primes ≤ B.
- ℓ_i = ln p_i and L = ln N.
- e ∈ Z^n gives u = ∏_{e_i>0} p_i^{e_i}, v = ∏_{e_i<0} p_i^{−e_i} (coprime), and w = u − vN.
- A *fac-relation* is a triple with |w| P-smooth (w ≠ 0 automatically). It gives the congruence u ≡ w (mod N)
  between P-smooth integers.

**Lattices.**
- Schnorr 1991 / Adleman 1995 weights: Λ_C = {b(e) = (√ℓ₁e₁, …, √ℓ_n e_n, C⟨e, ℓ⟩)}, target t = (0, …, 0, CL).
- Schnorr 2021 uses diagonal f(i), a random permutation of 1..n, and C = N′ = N^{1/(n+1)}. It uses SVP on the
  (n+1)-dimensional embedding.
- Adleman 1995 uses SVP with an extra row M(ln N + γ) and M = N⁴. We checked his manuscript, which is in Ducas's
  SchnorrGate repository as `LATFAC.pdf`, March 1995 draft.
- Write Q(e) = Σℓ_i e_i², A(e) = Σ|e_i|ℓ_i = ln(uv), and Δ(e) = ⟨e, ℓ⟩ − L = ln(u/(vN)).

**Theorem R [PROVED].**
- (i) ‖b(e) − t‖² = Q(e) + C²Δ².
- (ii) w = vN(e^Δ − 1) = 2e^{(A+L)/2} sinh(Δ/2).
- (iii) A ≤ Q, with equality iff uv is squarefree.
- (iv) gcd(uv, w) = 1.
- (v) |uw| ≥ N − 1 for any integers u − vN = w with v ≠ 0 and u, w ≠ 0.

*Proof.*
- (i) is by definition.
- (ii): u/(vN) = e^Δ and uv = e^A. So vN = e^{(A+L−Δ)/2}, and vN(e^Δ − 1) = e^{(A+L)/2}(e^{Δ/2} − e^{−Δ/2}).
- (iii): |e_i| ≤ e_i² for integers.
- (iv): a prime dividing u and w divides vN, which is impossible because gcd(u, v) = 1 and p ∤ N. Symmetrically, a
  prime dividing v and w divides u.
- (v) (codex): N ≤ |v|N = |u − w| ≤ |u| + |w| ≤ |uw| + 1. ∎

Tested in `test_residual_identity_exact` and `test_uw_lower_bound`.

**Relation lattice.**
- Λ_N = {(a, s) ∈ Z^n × Z : (−1)^s ∏p_i^{a_i} ≡ 1 (mod N)}. Its index is |⟨−1, P⟩ mod N| ≤ φ(N).
- A fac-relation is the Λ_N-vector exps(u) − exps(|w|) (with sign). Its weighted ℓ1 norm is ‖a‖_{1,ℓ} = ln(u|w|).
- So (v) says **λ₁ ≥ ln(N − 1)** in this seminorm, over vectors with a ≠ 0. ((0, 2) ∈ Λ_N has seminorm 0, so the
  bound is stated on the projection π(Λ_N) = {a : ∏p^a ≡ ±1}.) Every relation, found by any method, has u|w| ≥ N − 1.
- Minkowski in R^n, applied to π(Λ_N), whose index is ≤ φ(N), with the convex body {Σ|a_i|ℓ_i ≤ R} of volume
  2^n R^n/(n!∏ℓ_i), gives the upper bound λ₁ ≤ (n!·∏ℓ_i·φ(N))^{1/n} [PROVED].
  - For n ≈ L this is ≈ L(ln L + lnln L − 1 + o(1)). (Critic's numerics: λ₁/L ≤ 6.20, 7.02, 8.10 at 400, 800 and
    2000 bits.) So a relation u ≡ ±w with u|w| ≤ N^{(1+o(1)) lnln N} always exists.
  - That relation has v = (u − w)/N *arbitrary*, not smooth. So it is a relation for Adleman's form
    (u − kN^γ = w), not necessarily a vector of Schnorr's lattice.
  - These are rigorous existence statements. They do not guarantee small relations, and they do not count them.

**Gaussian heuristic warning [HEURISTIC + MODEL].**
- For n ≫ L the short lattice vectors are mostly {0, ±1}-vectors.
- The number of lattice vectors with ln(uv) ≤ T is a count of smooth pairs, which grows like e^{T(1−1/κ)+o(T)} for
  B = (ln N)^κ.
- The volume (Gaussian-heuristic) prediction grows like T^n. These differ greatly when n ≫ T (§3).
- Ducas (RISC slides, May 2021, slide 18) makes the same warning ("Gaussian Heuristic seems invalid for certain C").

## 2. The reduction theorems

### 2.1 Theorem P (promise reduction) [PROVED]

Let e⁰ ∈ {−1, 0, 1}^n be a *witness*: a squarefree relation with A₀ = ln(u₀v₀) and 0 < |w₀| ≤ ω. Put
C := e^{(L+A₀)/2}/ω. Then for every γ ≥ 1, every lattice vector b(e) with ‖b(e) − t‖ ≤ γ·dist(t, Λ_C) satisfies

  |w(e)| ≤ ω·R·exp((R² − A₀)/2 + R/(2C)),   with R = γ√(A₀+1),  R² − A₀ = (γ²−1)A₀ + γ².

*Proof.*
- |w₀| = 2e^{(A₀+L)/2}|sinh(Δ₀/2)| ≥ e^{(A₀+L)/2}|Δ₀|. So C|Δ₀| ≤ 1, and dist(t, Λ_C)² ≤ Q(e⁰) + C²Δ₀² ≤ A₀ + 1.
- For the answer e: A(e) ≤ Q(e) ≤ R² and |Δ(e)| ≤ R/C.
- Since sinh y ≤ y·e^y for y ≥ 0, we get |w| ≤ e^{(A+L)/2}|Δ|e^{|Δ|/2} ≤ e^{(R²+L)/2}(R/C)e^{R/(2C)}.
- Substitute C. ∎

**Corollary P1 (exact oracle).**
- A₀ is unknown, so scan j ∈ [⌊L⌋ − 1, A_max] (A_max = poly(L) is part of the promise).
- For each j, put ω_j = B·e^{−0.7}/√(j+2) and C_j = e^{(L+j)/2}/ω_j.
- At j = ⌊A₀⌋, write δ = A₀ − j ∈ [0, 1). Then C_j²Δ₀² ≤ e^{−δ}, and every exact-CVP answer has
  |w| ≤ ω_j√(A₀+1)·e^{(δ+e^{−δ})/2} ≤ ω_j√(A₀+1)e^{0.684} < B. This is the critic's sharper accounting, replacing
  our draft's "extra e^{1/2}", which was wrong.
- So under the promise that a witness with |w₀| ≤ ω_{⌊A₀⌋} exists, every exact-CVP answer is a fac-relation, since
  every integer of absolute value ≤ B is B-smooth.
- w ≠ 0 because N has no prime factor ≤ B.
- The scan costs O(A_max) calls.

**Rounding [PROVED].**
- Use integer entries round(S√ℓ_i), round(S·C·ℓ_i) and round(S·C·L).
- Q changes by a factor 1 + O(1/S).
- The last coordinate changes by at most (‖e‖₁ + 1)/(2S) ≤ (Q/ln 2 + 1)/(2S) in units of CΔ.
- With S ≥ L² both perturbations are o(1) for every vector with Q = O(L), so the bounds hold with o(1) changes.
- Entries have O(L) bits.

This is Adleman's §2 argument, generalised to γ and stated as a promise reduction. Adleman: B = ln^{20}(2N^{1.25}),
residual ≤ 72 ln²(N²), exact SVP. Schnorr 1991's §3 conditions are similar. Codex gives them as c > 1 and Schnorr's a > (2c−1)/(c−1), where Schnorr's
a is our κ.
Codex's equation (2) is the same certificate.

### 2.2 Lemma LA (linear algebra) [PROVED]

Let r₁, …, r_m ∈ Λ_N be relation vectors generating Λ′. Let K ⊆ GF(2)^m be the kernel of the parity map
c ↦ Σc_i r_i mod 2.
- For c ∈ K, lift to {0,1}^m and put y(c) = ½Σc_i r_i and x(c) = (−1)^{y_s}∏p^{y(c)} mod N. Then x(c)² ≡ 1.
- **c ↦ x(c) is a group homomorphism K → μ₂(Z/N).** Lifts of c + c′ differ from c + c′ by 2(c ∧ c′), which changes y
  by an element of Λ′, which maps to 1.
- So if some x(c) ∉ {±1}, a kernel *basis* already contains such an element.
- If [Λ_N : Λ′] is odd, every z with 2z ∈ Λ_N has an odd multiple in ½Λ′. So the image of K is all of
  {x ∈ ⟨−1, P⟩ : x² = 1}.
- Under GRH, ⟨P⟩ = (Z/N)^* once B ≥ 2 ln² N (Bach 1990). This group has ≥ 4 square roots of 1 when N has two odd
  prime factors.
- So under GRH plus odd index, N is factored deterministically from a kernel basis. ∎

`squares_dependency_factor` tests all kernel basis vectors. Tests: `test_lambda_vectors_and_factoring_toy` and
`test_gf2_kernel_basis_suffices`.

The odd-index hypothesis is **H_LA**. It is unproved for relations found by lattices, and it is the standard
"random relations" assumption of every congruence-of-squares method.

### 2.3 Distinctness [PROVED] and a flaw in Adleman's distinctness step [PROVED; no prior statement found by us]

**Our multipliers.**
- Use targets ln(kN) for primes k ∈ (B, 4B]. A relation u − v k N = w is still u ≡ w (mod N).
- Two targets k ≠ k′ giving the same (u, w) would force vk = v′k′. That is impossible, since v, v′ are P-smooth and
  k, k′ > B are distinct primes.
- So the m relations are pairwise distinct, without any heuristic.

**Adleman 1995 §2.** He obtains different congruences by "replacing V₀ by V₀ + V_i + V_j and rerunning".
- That does not change the lattice ⟨V₀, …, V_z⟩_Z: it is a unimodular basis change. So the shortest vectors, and any
  SVP oracle's output set, are unchanged.
- `test_adleman_rerun_lattice_unchanged` checks this.
- The repair is to delete V_i, V_j and use target ln(p_i p_j N). This restores his *heuristic* argument.
  - Distinctness then stays heuristic, because p_i p_j ∈ ⟨P⟩. A squarefree relation with p_i p_j p_{i′} p_{j′} | v can
    be returned for both pairs.
  - Rigorous distinctness needs multipliers outside P, as above.

### 2.4 Theorem H (the heuristic reduction) [HEURISTIC: H_S for witnesses, H_LA; GRH]

- Take B = (ln N)^κ with κ > 2, n = π(B), and targets kN for primes k ∈ (B, 4B]. There are ≥ n + 2 of them.
- **H_S** (u = vkN + w is smooth with the local density of smooth numbers, as tested in E1) predicts, for v ≤ N^a,
  N^{a(1−2/κ) − 1/κ − o(1)} squarefree witnesses with |w₀| ≤ ω. This → ∞ for a > 1/(κ−2).
  This matches Schnorr 1991's condition κ > (2c−1)/(c−1) (his "a" is our κ) with v = N^{c−1}, i.e. our
  a = c − 1.
- So every target has a witness with A₀ ≤ (1 + 2a + o(1))L.
- Corollary P1 (exact CVP, O(L) scales each) then gives n + 2 distinct relations, and Lemma LA factors N.
- Theorem P allows the same with any γ satisfying (γ² − 1)A₀ ≤ 2 ln B − 3 ln L − O(1).
  - With A₀ = (1+2a)L this is γ − 1 ≤ c_κ ln B/L, with c_κ = (1 − 3/(2κ))/(1+2a) + o(1).
  - The −3 ln L can be weakened to −ln L − O(1).
  Take ω := B·e^{−(γ²−1)A₀/2}/(γ√(A₀+1)e^{γ²/2+o(1)}) ≥ 1. The witness count is linear in ω, so witnesses with
  |w₀| ≤ ω still exist under H_S.
- In terms of the dimension, ln B/L = κ ln L/L = n^{−1/κ+o(1)}. So:

> **factoring ≤ CVP_{1+n^{−ε}} on Schnorr–Adleman prime lattices of dimension n = (log N)^{max(2, 1/ε)+o(1)},
> for any fixed κ > max(2, 1/ε), with O(n log N) oracle calls. Heuristic in H_S and H_LA; uses GRH (Bach) and
> N not a prime power; deterministic otherwise.**

### 2.5 Theorem G (cost as a function of γ) [MODEL under H_res′] and Prop A

**H_res′.** A "benign" γ-oracle returns a vector whose residual behaves like a random integer of the size allowed by
Theorem P, subject to the forced local conditions (coprime to uv). E2 tests this at toy sizes: pooled ratio 0.87.

**Cost.** The size is ln|w| ≈ θL + ln B with θ = (γ² − 1)(1 + 2a)/2, and a → 0 once B ≥ L^{ω(1)}. The method is then
"smooth u for free, w of size N^θ". With B = L[α, b] (so n = L[α, b]):
- smoothness costs exp(u ln u) with u = θL/ln B;
- n relations and n² sparse linear algebra are needed.

Balancing gives (derivation in §2.6):

| oracle factor γ(n) | factoring cost (unit-cost oracle; for oracle cost n^c multiply each L-constant by √(1+c)) |
|---|---|
| 1 + n^{−ε} | poly(log N) (Theorem H, even adversarial) |
| 1 + τ·L^{2α−1}(ln L)^{1−2α} ≈ 1 + (log n)^{2−1/α+o(1)} | L[α, 2√(ατ)] |
| 1 + (log n)^{−c} | L[1/(2+c)]; α < 1/3 ⇔ c > 1 |
| 1 + η, η = Θ(1) | L[1/2, √(2θ)] = L[1/2, √(γ²−1)], since θ → (γ²−1)/2: beats QS iff γ < √2 |
| γ ≥ 2; poly(n); 2^{n^ε}; 2^{n/log n} | residual up to N^{Θ(γ²)}: **no factoring gain**; for γ = poly(n) and n = L^{Ω(1)}, worse than trial division |

*The table assumes an oracle costing n^{o(1)} per call on dimensions n = L[α]. With cost n^c, the oracle work is
L[α, (1+c)b + s], and the optimum becomes L[α, 2√((1+c)ατ)] (proof critic). The exponents α are unchanged. No such
oracle is known: all known algorithms for γ = 1 + O(1) cost 2^{Θ(n)}.*

**Prop A (tightness / adversarial oracle) [HEURISTIC H_S].**
- Apply the counting of §3 to vectors with A ≈ γ²A₀ and C|Δ| ≤ 1. These have v ≈ e^{(γ²A₀−L)/2}.
- H_S predicts ≥ 1 such lattice vector with |w| ≥ ω·e^{(γ²−1)A₀/2 − O(ln L)}, and many once (γ² − 1)A₀ ≫ ln B.
- They lie in the γ-ball, and their residuals are smooth with probability ρ(·) < 1. So an adversarial γ-oracle can
  return a non-relation.
- Hence, heuristically, this query strategy with worst-case oracles works exactly in the regime
  γ − 1 = O(ln B/L) of Theorem H. Prop A says nothing about strategies that query other lattices or embeddings.

### 2.6 Derivation of the L[α] row [MODEL]

- Write ln B = bL^α(ln L)^{1−α}. The yield per call is exp(−u ln u) with u = θL/ln B.
- Setting u ln u = sL^α(ln L)^{1−α} gives u = (s/α)L^α(ln L)^{−α}, hence θ = (sb/α)L^{2α−1}(ln L)^{1−2α}.
- Oracle calls cost L[α, b + s] and linear algebra L[α, 2b]. The optimum is s = b, giving L[α, 2b] with θ = (b²/α)·(…).
- So θ = τL^{2α−1}(ln L)^{1−2α} gives L[α, 2√(ατ)]. α = 1/2 gives √(2τ), α = 1/3 gives 2√(τ/3).
- To beat GNFS 1.923 at α = 1/3 needs τ < 2.77. That assumes a unit-cost oracle; with cost n^c, the requirement
  is τ < 2.77/(1+c).

### 2.7 Prop NP: what is already rigorous [PROVED from known theorems]

- The language {(N, k) : N has a divisor in [2, k]} is in NP.
- Exact CVP is NP-hard (van Emde Boas 1981), and GapCVP_γ is NP-hard for γ = n^{c/log log n} under deterministic Karp
  reductions (Dinur–Kindler–Raz–Safra, Combinatorica 23 (2003)).
- Binary search therefore gives **factoring ≤_T GapCVP_{n^{c/lglg n}}** on lattices of dimension poly(log N),
  rigorously and deterministically.
- SVP: NP-hard under randomized reductions (Ajtai 1998; Micciancio 2001 for γ < √2; Khot 2005 and Haviv–Regev 2007 for
  every constant). So factoring ≤_R GapSVP_γ for every constant γ.
  - Adleman's 1995 abstract presents "factoring ≤_R Lattice (SVP)" as "the first substantial evidence that Lattice is
    not in P". Ajtai's 1998 theorem superseded it.
  - Ducas's slides (slide 21) note that Ajtai's and Micciancio's NP-hardness proofs reuse the Schnorr–Adleman lattice.
- **Consequences for the high-profile question.**
  - (i) Theorem H's γ = 1 + n^{−ε} is far inside the NP-hard range, so it adds nothing in γ. Its only new content is
    that the oracle is needed only on an explicit lattice family, of dimension (log N)^{κ} with κ > max(2, 1/ε).
    That is incomparable with the general-lattice NP-hardness reduction.
  - (ii) For γ ≥ √(n/log n), GapCVP_γ ∈ NP ∩ coAM (Goldreich–Goldwasser 2000), and ∈ NP ∩ coNP for γ ≥ √n
    (Aharonov–Regev 2005). LWE-based cryptography lives at γ = Õ(n) to poly(n) (Regev 2005).
  - No reduction from factoring is known there. Such a reduction would not contradict anything, since factoring is
    in NP ∩ coNP.
  - The prime-lattice route cannot supply one (Theorem G / Prop A): the residual is *exponentially* sensitive to
    the approximation factor, because log|w| is linear in the lattice norm.
  - (iii) Conversely, nothing is known to follow for lattice problems from faster factoring.
    - Prime-lattice CVP asks for smooth u near smooth multiples of N, and N's factors do not help.
    - Ducas–Pierrot (DCC 2019) decode *discrete-log* lattices by factoring small integers (rational reconstruction).
      That is the trivial direction, and it is a different family.
- **Worst case vs these structured lattices.**
  - Theorem H needs an oracle only on the prime-lattice family. A poly-time CVP_{1+n^{−ε}} algorithm on that family
    alone would (heuristically) factor; NP-hardness says nothing about this family.
  - We know of no evidence either way that CVP on prime lattices is easier than worst-case CVP.

**Clean statement (answer to Q4).** Factoring ≤ γ-CVP holds rigorously for γ ≤ n^{c/lglg n} on general lattices.
On prime lattices it holds heuristically for γ = 1 + n^{−ε}, and for no γ that is bounded away from 1 in the worst
case. No reduction, clean or heuristic, exists through this lattice family at the approximation factors relevant to
post-quantum cryptography. The reason is quantitative: Theorem P is tight (Prop A).

## 3. Counting relations: when do they exist at all? [HEURISTIC H_S; MODEL]

**Asymptotics.**
- For B = (ln N)^κ, Ψ(X, B) = X^{1−1/κ+o(1)} for X = N^{O(1)}.
- The expected number of relations with v ≤ N^a and |w| ≤ N^b is N^{a(1−2/κ) + b(1−1/κ) − 1/κ + o(1)} (H_S).
- Lattice (close-vector) relations have |w| < vN, i.e. b ≤ 1 + a. The exponent at b = 1 + a is a(2 − 3/κ) + 1 − 2/κ.

| κ | relations with \|w\| < vN | small residual \|w\| ≤ B |
|---|---|---|
| κ < 3/2 | N^{−Ω(1)}: **none** | none |
| 3/2 < κ < 2 | exist, all with v ≥ N^{(2−κ)/(2κ−3)} | none |
| κ > 2 | exist with v = N^{O(1)} | exist for v ≥ N^{1/(κ−2)} (Theorem H regime) |

- So the dimension needed for *any* Schnorr-type relation to exist is n ≥ (ln N)^{3/2+o(1)}.
- Exact CVP returns relations with polynomial yield only once n ≥ (ln N)^{2+o(1)}. For 3/2 < κ < 2 the closest vector
  has |w| ≈ N^{1/κ}, which is smooth with probability N^{−1/κ²}.
- With multiplier targets (v free), relations exist for every κ > 1 but need N^{(2−κ)/(2(κ−1))} targets.

**Finite sizes (E0, MODEL; H_S naive).** The total count includes every v and every |w| < vN, both signs.

| bits | n | E_tot |
|---|---|---|
| 40 | 47 | 1.9·10^{14} |
| 200 | 47 | 6.6·10^{−11} |
| 400 | 47 (Schnorr) | **1.4·10^{−52}** |
| 400 | 95 | 3.0·10^{−31} |
| 800 | 95 (Schnorr) | **1.5·10^{−115}** |

Sensitivity:
- Resolution K = 50 → 100 changes log₁₀ E_tot by 0.03.
- Widening the v-range changes it by 0.

n_min (least n on our grid with E_tot ≥ n + 1):

| bits | 64 | 128 | 200 | 256 | 400 | 512 | 800 | 1024 |
|---|---|---|---|---|---|---|---|---|
| n_min | 20 | 50 | 90 | 130 | **275** | 400 | **775** | 1125 |
| p_{n_min} | 71 | 229 | 463 | 733 | 1777 | 2741 | 5881 | 9049 |
| κ_eff = ln p/ln L | 1.13 | 1.21 | 1.20 | 1.27 | 1.33 | 1.35 | 1.37 | 1.39 |

- The local slope of ln n_min against ln L between 256 and 800 bits is 1.56, close to the asymptotic 3/2 + o(1).
- **Naive H_S overestimates primitive relations.** E1 shows naive ≥ observed primitive in all 9 cells, by factors
  6–28. So the E_tot values are heuristic *upper* estimates for what a lattice can find, since lattice vectors are
  always primitive.

## 4. Schnorr 2021, dissected

We checked the primary text: the "corrected" 2021/933 version of 16 June 2021, via a web.archive.org copy, since the
ePrint page blocked us.

1. **The claim.**
   - n = 47 for 400 bits, n = 95 for 800 bits.
   - One fac-relation per random permutation f.
   - Block size 24; 4.2·10⁹ and 8.4·10^{10} operations.
   - "This destroys the RSA cryptosystem."
   - At 800 bits he concludes |u − vN| ≲ 522.9 ≈ p_n^{1.0075}, assuming v ≲ p_n. So each lattice yields a relation
     with probability ≈ ρ(1.0075) ≈ 0.99.
2. **What exists (MODEL, H_S).**
   - In his own box (v ≤ 499, |w| ≤ 523) the expected number of relations is **3.7·10^{−176}** (800 bits), and
     9.2·10^{−85} at 400 bits for the box (v ≤ 211, |w| ≤ 250).
   - Over all sizes it is 1.5·10^{−115} and 1.4·10^{−52} (§3). He needs 96 and 48.
   - The reason is |uw| ≥ N − 1: with |w| ≤ 523, u is a 499-smooth integer ≥ N/523. Smooth integers near 2^{800} with
     B = 499 have density 3.6·10^{−180}, and there are only 1.9·10^{62} of them below 2^{800}.
   - This failure does not depend on the lattice algorithm: **no algorithm can find relations that heuristically do
     not exist.**
   - The required dimensions are n ≥ 275 (400 bits) and n ≥ 775 (800 bits), even before asking a lattice to find them.
3. **Where the argument breaks [PROVED].**
   - Lemma 3.1's proof writes "u − vN′ < v|ẑ|/(1 − ε/2)" with N′ in place of N. The correct bound is
     |u − vN| < (N/N′)·v|ẑ|/(1 − ε/2): a missing factor N^{n/(n+1)}, as in Theorem R(ii). Codex found this
     (schnorrcvp_1 §3); we confirmed it in the primary text.
   - The 800-bit estimate "‖b₁‖ < 0.917" is impossible. Every nonzero vector of L(R_{n,f}) has a coordinate f(i)e_i
     with |f(i)e_i| ≥ 1, or is supported on the last basis vector with norm ≥ N′ ln N (codex).
   - The paper's (3.3) "‖b‖² ≥ ln uv + ẑ²" holds for diagonal √ln p_i, not for the displayed diagonal f(i). This is an
     internal inconsistency.
4. **SVP vs CVP.**
   - [MEASURED by Ducas; heuristic in general] In the homogeneous (n+1)-dimensional SVP formulation the shortest
     vectors have last coefficient 0. They are pairs
     u ≈ v ("smooth near-twins": Ducas observed "a single 1 and a single −1").
   - The target enters only with coefficient ±1, which needs a BDD-type promise. That promise is the existence of a
     witness, which does not hold (item 2).
   - The smooth-twin lattice of Mulder–Sterner–van Woerden (arXiv 2509.17699) uses exactly these vectors, for a
     different purpose.
5. **What the approach actually costs (corrected exponent) [MODEL].**
   - Asymptotically, with any CVP/BKZ method of cost ≥ 2^{c·n} for approximation γ = 1 + O(1), it costs
     ≥ exp((ln N)^{3/2+o(1)}), because n ≥ n_min. With multiplier targets the floor is
     exp(Θ(ln N·lnln N/lnlnln N)). Both are super-polynomial in N, i.e. worse than trial division.
   - **As a function of β [MODEL: GSA + H_res′; the weakest step, since prime lattices are not GSA-typical].**
     BKZ-β reaches only γ_β ≈ δ_β^{Θ(n)} with δ_β ≈ (β/2πe)^{1/(2(β−1))}. Then θ_β ≈ γ_β² and |w| ≈ N^{θ_β}, which
     is useless unless β = n − o(n). So T(β) ≥ max(2^{0.292β}, ρ(θ_β L/ln B)^{−1}) is minimised only at
     β ≈ n ≥ n_min. At toy dimensions (n ≤ 40, E2) LLL-Babai outputs still had median |w| ≈ N^{0.64}, because
     δ^n is small there.
   - **Concrete (E0b, MODEL, optimistic for the method).**
     - V(n) is the number of lattice vectors an *ideal* unit-cost enumerator (increasing |w|) must output to see
       n + 1 relations. It assumes naive H_S, which overestimates relations.
     - Lattice cost per sieve-CVP call: 2^{0.292n}, yielding up to 2^{0.2075n} close vectors.
     - So log₂T = 0.292n + max(0, log₂V − 0.2075n).

     | bits | best n | log₂T | log₂V there | crude log₂ QS / GNFS |
     |---|---|---|---|---|
     | 128 | 100 | 42.8 | 34.4 | 28.8 / 33.7 |
     | 256 | 300 | 104.2 | 78.8 | 43.7 / 46.7 |
     | 400 | 600 | 186.7 | 136.0 | 57.0 / 57.2 |
     | 512 | 800 | 264.1 | 196.5 | 65.9 / 63.9 |
     | 800 | 1600 | 467.4 | 332.2 | 85.4 / 77.9 |
     | 1024 | 1600 | 699.1 | 563.9 | 98.5 / 86.8 |

   - Even the enumeration count V alone (with free lattice work) is 2^{136} at 400 bits.
   - The 1024-bit row and n > 1200 at 800 bits come from `run_e0b_ext` (§7.3).
   - The n-grid is coarse (…, 1200, 1600, 2000, 2400, …), so these are grid minima, accurate to a few tens in log₂T.
   - The "crude" L-values drop o(1). They are reference points only.
   - At 800 and 1024 bits the modelled cost exceeds trial division (2^{400}, 2^{512}).
6. **What Ducas's analysis adds.** We read his RISC slides (2021-05-20).
   - Pitfalls: projection, ℓ1 vs ℓ2, invalid Gaussian heuristic, and "p_i | u ∨ p_i | v ⇒ p_i ∤ S" (our Theorem
     R(iv); E2's H_res′).
   - His experiment on log S/log N against n shows that beating QS's non-lattice cost needs n ≥ 50 already at
     100 bits.
   - His verdict: "a solid average-case complexity analysis is still missing". Our §3/E0 counting is a step
     toward it. It does not prove anything. It turns the question into Ψ-counts that can be computed.

## 5. Rigour map (Q2)

| Step | Status |
|---|---|
| Distance–residual identity, coprimality, \|uw\| ≥ N−1, λ₁(Λ_N) bounds | PROVED |
| Promise reduction (witness ⇒ every γ-CVP answer has \|w\| ≤ bound); exact oracle ⇒ relation | PROVED (incl. rounding) |
| Distinctness via prime multipliers k > B | PROVED |
| Linear algebra: kernel basis suffices; deterministic factor under odd index + GRH (Bach) | PROVED (conditional) |
| Existence of witnesses for each target (smooth u within B/√L of smooth multiples vkN) | HEURISTIC H_S; unproved, since it is a smooth-numbers-in-short-intervals statement far beyond known results; E1 supports the density model within ×1.6 (×0.92 post-hoc) |
| Odd index / relations generate Λ_N up to odd index | HEURISTIC H_LA (standard) |
| Residual of an approximate answer is smooth like a random integer (benign oracle) | HEURISTIC H_res′; E2 within ×0.68–1.75 per variant |
| Cost tables (Theorem G, E0b) | MODEL |
| NP-hardness subsumption | PROVED from DKRS 2003 / Khot 2005 / Haviv–Regev 2007 |

**Cleanest provable statement.** Let N not be a prime power and have no prime factor ≤ B, with B ≥ 2 ln²N. Assume:
- (i) for each of n + 2 primes k ∈ (B, 4B] a squarefree witness exists for target kN, with A₀ ≤ poly(L) and
  |w₀| ≤ ω_{⌊A₀⌋} (§2.1);
- (ii) the resulting relations generate an odd-index sublattice of Λ_N.

Then an exact-CVP oracle on the prime lattices of dimension n = π(B) factors N in deterministic polynomial time
(O(n·poly(L)) oracle calls), under GRH. Under H_S, (i) needs B = (ln N)^{κ} with κ > 2. So B = 2 ln²N is the
rigorous floor for the linear algebra, not for the witnesses. (i) is the whole difficulty. It is a statement about smooth integers in intervals of
length B/O(√L) around smooth multiples of N, which is not known for any B = (log N)^{O(1)}.

## 6. Variants (Q5)

| Variant | Effect | Cost | Status |
|---|---|---|---|
| Multiplier targets ln(kN) (v free) | existence threshold κ > 1 instead of 3/2; rigorous distinctness | still exp(Θ(L lnL/lnlnL)) or more with 2^{Θ(n)} CVP | MODEL |
| Adleman SVP form (u − kN^γ = w) | same relations | same; its distinctness step is flawed but repairable (§2.3) | PROVED / MODEL |
| ℓ1 or other norms | log\|w\| is linear in any norm of e, so the exponential sensitivity is unchanged | none | Theorem P upper bound PROVED; tightness HEURISTIC (Prop A) |
| NFS rational side via prime lattices (list smooth a − bm near bm) | turns a poly-delay enumerator of dimension L[1/3] into the factory constant 1.6386 (smoothenum T3) | needs poly-delay listing in dim L[1/3]; known enumeration is 2^{Θ(n)} = exp(L[1/3]) | MODEL / HYPOTHETICAL |
| Algebraic side (S-unit lattice CVP for norms of a − bα) | would approach smoothenum T6/T7 (0.763 / 0) | needs class group / S-unit computations plus principality, and dimension L[1/3] | SPECULATIVE |
| Two lattices (u-side and w-side) | w-smoothness cannot be imposed linearly: w = u − vN is not multiplicative in e | — | assessment |
| The project's devices (Prop LQ, Prop JR, Lemma C2) | the prime lattice is a "list u-smooth points" device, so resbarrier's D = 2 hyperbola applies to the *density* (\|uw\| ≥ N − 1), not to the cost | — | assessment |

No variant is plausibly better. Any gain requires a fast near-exact CVP algorithm specific to prime lattices.
- Only a *polynomial-time* CVP_{1+n^{−ε}} algorithm on this family gives polynomial-time factoring (Theorem H).
- A 2^{n^c} algorithm at n = L^{κ} gives exp(L^{κc}) instead.
- No converse reduction (CVP on prime lattices ≤ factoring) is known.

## 7. Experiments

### 7.1 E1: exact relation counts (pre-registered, MEASURED)

Setup: 9 cells × 10 balanced semiprimes, exact enumeration of every relation with v ≤ 10³ and |w| ≤ 10⁵ (2·10⁴ for
n = 25).

| (bits, n) | all | primitive | pred H_S′ | ratio | pred naive | post-hoc H_S″ ratio |
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

- **E1-P1 passed:** pooled 1.63 ∈ [1/2, 2], and every cell with prediction ≥ 10 is within [1/3, 3].
- **E1-P2 passed:** naive ≥ observed primitive in every cell.
- The systematic factor ≈ 2 has a simple cause: exactly one of u, v, w is even, so parity is *forced*.
  - The POST-HOC model H_S″ (`psi.local_factor`) puts u uniform on the p − 1 residues left by coprimality, so
    g₂ = 2f₂. It gives 0.92 pooled.
  - Most relations in a box are non-primitive (gcd(v, w) > 1), i.e. multiples of smaller ones. Lattices only see
    primitive ones.

### 7.2 E2: relation yield of CVP outputs (pre-registered, MEASURED)

Setup: 6 sizes (20–40 bits) × 3 N × n ∈ {12, 20, 30, 40} × 5 variants × 10 permutation seeds = 3600 lattices.
Outputs per lattice: LLL-Babai, BKZ-10-Babai, and the ≤ 20 closest vectors (exact enumeration, mpfr GSO, exact
re-sort). That is 78,963 outputs.

| variant | outputs | relations | Σ p_ref | obs/H_res′ | obs/naive | post-hoc H_S″ |
|---|---|---|---|---|---|---|
| S21perm (diag f(i), C = N^{1/(n+1)}) | 15840 | 296 | 169 | 1.75 | 1.34 | 1.03 |
| S21sqrt | 15840 | 379 | 255 | 1.49 | 1.09 | 0.91 |
| A0.5 (√ln p_σ(i), C = N^{1/2}) | 15840 | 2506 | 2027 | 1.24 | 0.90 | 0.94 |
| A1.0 | 15840 | 3494 | 4277 | 0.82 | 0.74 | 0.73 |
| A1.25 | 15603 | 3113 | 4555 | 0.68 | 0.73 | 0.62 |

- **E2-P1 passed** (pooled 0.87; every variant in [1/3, 3]).
- **E2-P2 passed** (refined 0.87 vs naive 0.80).
- **E2-P3 passed** (yield at 40 bits < yield at 20 bits in every cell).
- **E2-P4 passed:** median log₂|w|/bits for Schnorr-2021 scaling at ≥ 32 bits is 0.96 (perm) and 0.92 (sqrt).
  N′ = N^{1/(n+1)} cannot target small residuals: |w| ≈ vN/N′·(…).
- The post-hoc local model fits S21 and A0.5. It over-predicts A1.0/A1.25 by ×1.4–1.6, probably because the 20 closest
  vectors of one lattice are correlated (e* + short z). Not pursued.

**Yield per list output (selected).**
- n = 40, A1.0: 0.82 / 0.66 / 0.50 / 0.40 / 0.16 / 0.085 at 20 / 24 / 28 / 32 / 36 / 40 bits.
- n = 40, S21sqrt: 0.18 → 0.000.
- Exact CVP (list0) residuals for A-variants have median |w| ≈ N^{0.475}, so these toy dimensions are far below
  the Theorem H regime.
- LLL-Babai produced 10 outputs with u or v of more than 4000 bits, an image of the exponential blow-up of
  approximate answers.

**E2-P5 (duplication).** Across the 10 permutation seeds, the 700 exact-CVP relations in cells with ≥ 3 are only 416
distinct (59%). Random weights do not reliably give distinct relations, as Aboumrad–Widdows–Kaushik also observed.
§2.3's prime multipliers do, provably.

**E2-P6 (correctness, not evidence of speed).**
- All 15 N of 20–36 bits were factored by the GF(2) square-root step from the pooled relations. Factors are
  STRUCTURAL: N's prime factors exceed p₄₀ = 173.
- At 40 bits only 17–22 distinct relations were found (< 42 needed).

**Deviation (disclosed).** The first E2 run crashed at 40 bits on a Babai output with u > 10^{4300} (Python's
int-to-str limit). We patched serialization (huge integers stored as bit lengths only) and reran everything; the
run is deterministic. The crashed log is kept as `schnorrcvp_E2_crashed_run1.log`.

### 7.3 E0 / E0b (MODEL)

- E0: pre-registered, deterministic (§3).
- E0b (ideal enumeration cost + sieve-CVP cost) and its extension were added after the pre-registration, and are
  disclosed as such.

Commands:
- `.venv/bin/python -m src.schnorrcvp.run_e0`
- `run_e0b`
- `run_e0b_ext`
- `run_e1`
- `run_e2`
- `analyze e1|e2|posthoc_e1|posthoc_e2`

## 8. Coordination

- `suffmap` (notes/suffmap_results.md, read after our draft) lists this route as H9: "sufficiency for any
  approximation factor not established; details: schnorrcvp". Our answer for its map:
  - Approximate CVP on prime lattices is *sufficient* for polynomial time at γ = 1 + n^{−ε} (Theorem H, even for
    adversarial oracles).
  - It is sufficient for L[α < 1/3] at γ = 1 + (log n)^{−c} with c > 1 (benign oracle, cost poly(n)).
  - It is not sufficient at any γ bounded away from 1 for worst-case oracles, and it gives nothing at crypto-relevant
    γ.
  - Exact or near-exact prime-lattice CVP is a concrete instance of suffmap's E3. The CVP ranges over all smooth
    multipliers v ≤ N^{a} at once, which is a union of windows (Bohr set). So it inherits suffmap's Prop DLT: such an
    oracle would also give prime-field discrete logarithms, as Schnorr 1991's title already suggests.
  - **Consistency with suffmap Prop SI (1.2114).** Prop SI's constructor is one-sided with v = 1. There the residual
    must cover the gap 1/d(kN) between smooth numbers, and that existence constraint sets the constant. In Schnorr's
    lattice, v ranges over N^{o(1)} smooth values, so witnesses with |w| ≤ B exist (§2.4, a = o(1)). Theorem G's
    L[α, 2√(ατ)] carries no such floor, because the floor is paid by the γ = 1 + o(1) oracle. No contradiction.
- `abcsunit` (notes/abcsunit_results.md): relations here are S-unit solutions of u − vN = w. Evertse's bound
  3·7^{2s+3}, as cited there, caps their number by 2^{O(n)} for every N, far above n + 1. So finiteness gives no
  obstruction. §3 counts how many have small height, and the obstruction is that count.

## 9. Reviews

### 9.1 Codex blind round (schnorrcvp_1)

Independently found:
- the identity and certificate (Theorem R, P);
- the exact-CVP failure example N = 35;
- the Schnorr 1991 conditions;
- the Lemma 3.1 error and the 0.917 impossibility;
- |uw| ≥ N − 1.

It did not do the existence count, the dimension threshold or the NP subsumption. It labelled "universal threshold
depending only on γ(n)" UNVERIFIED. We supply such a threshold only for the family B = (ln N)^κ with weights ln p_i
(Theorem H/G), not in general.

### 9.2 Codex adversarial round (schnorrcvp_2)

[pending]

### 9.3 Opus critics

**Proof critic (read-only Opus).** No FATAL issues. It checked:
- Theorem R, Theorem P and the rounding argument;
- Lemma LA, including the odd-index image argument;
- k > B distinctness;
- the Adleman flaw, after reading his §2 (it found no reading that saves it);
- the Theorem H exponent;
- the §2.6 derivation;
- the §3 thresholds and floors;
- the Schnorr 2021 errors.

Two MAJOR issues, both fixed:
1. Overclaims in §0, §2.7 and §6:
   - "provably cannot" became "heuristically, this query strategy cannot";
   - "no link" became "no link known";
   - "weaker than NP-hardness for any lattice" became "weaker in γ; the new content is the explicit family";
   - "precisely as hard as the target" was removed.
2. Theorem G's constants assumed a unit-cost oracle. With oracle cost n^c they gain a factor √(1+c). The exponents
   are unchanged.

MINOR fixes applied:
- the γ = 1 + Θ(1) constant is L[1/2, √(γ²−1)], which beats QS iff γ < √2 (not η < 1/2);
- Corollary P1's scan constants: ω_j = Be^{−0.7}/√(j+2), extra factor e^{0.184}, range from ⌊L⌋ − 1;
- GRH added to Theorem H's labels;
- the λ₁ seminorm and projection, and the lower-order term of the Minkowski bound;
- the Adleman repair's distinctness is only heuristic;
- the explicit constant c_κ;
- the hypotheses of the cleanest provable statement;
- the dimension statement "any fixed κ > max(2, 1/ε)";
- the notation clash with Schnorr's "a";
- three mislabels;
- the certificate test now checks the bound;
- 0.92–0.96.

**Literature critic (read-only Opus).** [pending]

## 10. References (checked in primary form unless noted)

- C. P. Schnorr, *Factoring integers and computing discrete logarithms via diophantine approximation*, EUROCRYPT '91,
  LNCS 547; DIMACS 13 (1993). Via codex.
- L. M. Adleman, *Factoring and lattice reduction*, draft, 16 March 1995 (SchnorrGate `LATFAC.pdf`). Read in full.
- C. P. Schnorr, *Fast factoring integers by SVP algorithms, corrected*, ePrint 2021/933, version 16.06.2021. Read.
- L. Ducas, SchnorrGate (GitHub) and RISC slides, 20 May 2021. Read.
- Aboumrad, Widdows, Kaushik, arXiv 2308.07804, Table 1. Read via fetch.
- Vera, arXiv 1003.5461 (Theorems 1–2). Read via fetch.
- Ducas–Pierrot, DCC 87 (2019), ePrint 2018/146.
- Mulder–Sterner–van Woerden, arXiv 2509.17699.
- Yan et al. arXiv 2212.12372; Khattar–Yosri arXiv 2307.09651; Ritter–Rössner 1997 (ePrint 1997/008). Via codex,
  not re-checked.
- Dinur–Kindler–Raz–Safra, Combinatorica 23 (2003); van Emde Boas 1981; Ajtai STOC 1998; Micciancio SICOMP 2001; Khot
  JACM 2005; Haviv–Regev STOC 2007; Goldreich–Goldwasser JCSS 2000; Aharonov–Regev JACM 2005; Regev JACM 2009; Bach
  Math. Comp. 55 (1990). Standard. Statements are to be checked by the literature critic.
- Canfield–Erdős–Pomerance 1983. Codex checked the Ψ asymptotics.
