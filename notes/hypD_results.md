# hypD: relation supply for Dixon descent (Hypothesis D)

Agent `hypD`, phase 7, 2026-10-19.
- Pre-registration: `notes/hypD_prereg.md`. Commit `7e6f1ee` came before any registered run; Addendum A,
  commit `af85691`, came before E1x.
- Code: `src/hypD/`. Tests: `tests/test_hypD.py`.
- Raw data: `experiments/raw/hypD_e{1,1x,2,3a,3b}.jsonl`.
- GPT-6.1-sol: the blind prompt `codex/hypD_1.md` could not run because the Codex login expired. The
  prompts are written for a later run.
- Reviews: two independent read-only Opus critics, one on proofs and one on literature (§9).
  **GPT-6.1-sol review is pending.**

## Verdict

**We did not prove Hypothesis D, and we obtained no rigorous L[1/2, c < 1] factoring algorithm.**
Lenstra–Pomerance's L[1/2, 1+o(1)] stands. The main gain is a reduction: of the three missing ingredients
of `lemmaM`, two are removed by proved arguments, and what remains is one averaged second-moment statement.
That statement must hold only at LV's own precision, L^{o(1)}.

- **Theorem S [proved].** A uniform square R splits as u/v with u, v ∈ S with probability at least
  |S|⁴/(2·2^ω·φ(N)·Y²·(1+log Y)), for *every* set S ⊂ [Y/2, Y] coprime to N (Y² < N).
  - The proof is the Dixon/Pomerance positivity Σ_χ|Σ_S χ|² ≥ |S|², plus a gcd-multiplicity count.
  - Lemma Ω: deterministic Pollard–Strassen in L^{o(1)} time makes 2^ω = L^{o(1)}.
  - This removes **ingredient 3** (split counting under the saddle law with the square-class condition)
    for an arbitrary, field-dependent set of admissible splits.
- **Theorem K [proved reduction].** Use splits whose large primes are exactly k primes from one dyadic
  interval J ∋ N^{1/(2k)}, with k a fixed integer. Then Hypothesis D (all targets of a split descend in
  one field) follows from one averaged statement:
  - **D_k:** Pr_{A, ℓ∈J}[special-q search for ℓ succeeds in field A] ≥ L^{−o(1)}.
  - D_k needs no complete-target-list joint event. A field only needs a good set G_A of density
    L^{−o(1)}, and Theorem S then finds splits inside G_A.
  - This removes **ingredient 2**. It also answers codex lemmaM_2's objection ("two field types may descend
    complementary prime types"): such fields are fine, because each field only uses splits that lie in
    its own G_A.
  - With bounded k the density need only be L^{−o(1)}, which is LV's precision. With lemmaM's
    k ≍ (log N)^{1/4} it would have to be exp(−o((log N)^{1/4})).
  - Cost: θ(d, k) = 2s + (d−1)/(4k) + max(0, s − (d−2)/(4k)) with s = 1/(d+1). This tends to lemmaM's
    3/(d+1) as k → ∞.
- **Prop PZ [proved].** D_k ⟸ FM + PC, by Paley–Zygmund over the joint law of (A, ℓ):
  - (FM) a first moment E N_ℓ(A) ≥ 1 at box size T_SQ·L^{o(1)};
  - (PC) a pair-correlation bound E[N(N−1)] ≤ L^{o(1)}(E N)².
- **What is left.**
  - **FM** reduces to Harper's BV over dense moduli (rational side, a proved input) and to LV §5 plus LV
    Lemma 8.5 transplanted to the moduli family {ℓ·s : ℓ ∈ J, s B-smooth} (algebraic side). The latter is
    a **[proof sketch, not written out line by line]**. Because ℓ ≥ x^η, every conductor divisible by ℓ
    falls in LV's large-conductor range, and the rest of LV's proof goes through with an extra common
    factor Σ_{ℓ∈J}1/ℓ.
  - **PC is open** and is now the single missing estimate. It splits into two parts:
    - **PC_alg:** an algebraic-side pair correlation in one random field;
    - **PC_rat:** the variance over (ℓ, m) of rational-side smooth counts. The rational values of two
      special-q pairs lie on the fixed line b′s − bs′ = Δ/ℓ, so this is a binary problem: two linear forms
      in one variable, averaged over coefficients. No estimate is known at L[1/2] smoothness.
  - **Prop P [proved implication; §7.4].** The best published bound for PC_rat is Pascadi's
    ρ^{5/8−ε} for two linear forms, a loss of ρ^{−3/8} against independence. Even granting FM and
    PC_alg, it yields only c ≈ 1.17, using Theorem S′. Reaching c < 1 needs a loss exponent below ≈ 0.17. **So no combination of known theorems beats
    Lenstra–Pomerance through this design.**
  - **Prop G [proved; barrier].** The joint image lattice of two special-q pairs has index D ≥ ℓ. Therefore
    the geometric route of lemmaM §3 (Ball slicing plus lattice corrections) loses at least
    ℓ = N^{1/(2k)} and cannot prove PC_alg.
  - PC_alg needs smooth numbers to equidistribute modulo q = t·D on average over the pair family, i.e.
    a **restricted Barban–Davenport–Halberstam** estimate (Hypothesis RBDH; conjectured).
- **Design V [conditional theorem; new route found by us; no prior formulation found by us].** Take cheap
  LV relations, which are not fruitful, and add L^{o(1)} Vallée relations x² ≡ ∏p^e (mod N), which are
  fruitful by fibre symmetry. Then every N is factored in L[1/2, √(2/3) + o(1)] ≈ L[1/2, 0.816].
  - This needs only **RK**: the row space of the NFS relations contains the Vallée vectors with
    probability L^{−o(1)}.
  - It also needs LV relation supply at fixed odd degree d ≥ 9 at L[1/2] scale (first moment only, via LV
    stochastic deepening; a sketch) and the LV §6 square-root adaptation.
  - RK is a rank statement about non-fruitful relations. It needs no fruitfulness and no special-q.
- **Experiments (pre-registered, exact, toy).**
  - E2 (special-q counts in toy LV fields; 16 instances, 38 400 cells):
    - PC_alg 1.08–1.21 (P3 holds);
    - the Poisson-with-field-quality model fits γ to within 0.009 (P4 holds);
    - same-field pair ratio 1.02–1.07, so there is no complementary-type effect (P5c holds);
    - between-field overdispersion 2.2–4.2: P5a **failed in 1/16** (4.21 > 4);
    - Pr_A[g_A ≥ γ/2] ≥ 0.95 everywhere (P5b holds).
  - E3 (RK in a toy NFS; 20 instances):
    - 360/400 Dixon vectors lie in the NFS row space (P6 holds);
    - 0/5098 single-prime vectors e_p with p ∉ Q_N do (P7a);
    - 787/787 with p ∈ Q_N and ≥ 3 occurrences do (P7b).
  - E1 as registered was uninformative, a disclosed design error. Its exact replacement E1x shows the
    Theorem S bound never fails (0/120), random-S model ratios of 0.85–0.94 (P2′ holds), and the predicted
    2^ω positivity excess for QR-biased S. P2″ was vacuous at the registered |S| ≥ 200 threshold.

### Per-design constants (the L[1/2, c] exponent; "inputs" = what a proof would still need)

| design | c | status | missing inputs |
|---|---|---|---|
| Lenstra–Pomerance 1992 (class groups) | 1 | proved (prior art) | — |
| Vallée 1991 (small squares, θ = 2/3) | √(4/3) ≈ 1.155 | proved (prior art) | — |
| lemmaM (a), plain fixed-d descent | 2/√(d+1) (0.816 at d = 5) | heuristic | first moment at t ≍ F^{1−o(1)}: beyond every known level |
| lemmaM (a′), k ≍ (log N)^{1/4} | √(6/(d+1)) | conditional on Hyp. D (3 ingredients) | first moment, complete-target-list average, split counting; PC with loss exp(o((log N)^{1/4})) |
| **this note, (a′_k), fixed k** | √(2θ(d,k)) → √(6/(d+1)) (d = 7: 0.866; d = 9: 0.775; d = 11: 0.707) | **conditional on D_k ⟸ FM (sketch; ≈ de la Bretèche–Drappeau 2020 Thm 2.1) + PC (open)** | PC = PC_rat (binary smooth correlation) + PC_alg (⟸ RBDH, conjectured) |
| (a′_k) using only the *proved* PC_rat bound (Pascadi ρ^{5/8}) | ≈ 1.17 (best over d, k) | conditional on FM + PC_alg | — (the proved bound is too weak; threshold loss exponent ≈ 0.17) |
| **Design V (LV rows + Vallée rows)** | √(2/3) ≈ 0.816 (needs d ≥ 9) | **conditional on RK** | RK; LV supply at L[1/2] (sketch); LV §6 adaptation |
| Design V with uniform Dixon rows | 1 | conditional on RK | as above; no gain over LP |
| QS-style descent (lemmaM (d)) | 1 | heuristic | FM and PC for a *fixed* quadratic polynomial; no randomization to average over |
| one-level growing degree (lemmaM (b)) | L[2/5] | heuristic / conditional | bounded-k trick does not apply: d → ∞ forces k ≫ d |

Under PC + FM, c → 0 as the fixed odd d grows, because (a′_k) is an L[1/2] analysis at fixed d. This is
consistent with the L[2/5] and L[1/3] heuristics; it is not a separate claim.

---

## 1. Setting

- Dixon descent follows `lemmaM` §1–§2.
  - R = r² mod N with r uniform; R is split as u/v (mod N).
  - The large primes ℓ of uv are descended by special-q relations in an LV-randomized field
    A = (m, f = f̂_{N,m} + (x − m)R_A).
  - Theorem DD makes every dependency that contains a descent relation fruitful with probability
    1 − 2^{1−ω}, provided the transcript is invariant under r ↦ εr.
  - All procedures below act on R and A only, so (T) holds.
- L = L_N[1/2, 1]. "L^{o(1)}" means exp(o((log N log log N)^{1/2})), uniformly in N.
- ω = ω(N). Q_N = (Z/N)^{×2}, of index h = 2^ω in (Z/N)^× for odd N.

## 2. Split counting

**Lemma Ω [proved; elementary].** Put y₀ = exp((log N)^{1/2}(log log N)^{1/4}).
- Deterministic Pollard–Strassen finds every prime factor ≤ y₀ of N in time y₀^{1/2+o(1)} = L^{o(1)}.
- Removing these factors, and prime powers, leaves an N with every prime factor > y₀. Then
  ω ≤ log N/log y₀ = (log N)^{1/2}(log log N)^{−1/4}.
- Hence 2^ω ≤ exp(0.7(log N)^{1/2}(log log N)^{−1/4}) = L^{o(1)}. ∎

**Theorem S [proved].** Let N be odd and Y ≥ 2 with Y² < N, and let S ⊂ [Y/2, Y] ∩ Z be any set of
integers coprime to N. Let R be uniform on Q_N, and let E_S be the event that u ≡ Rv (mod N) for some
u, v ∈ S. Then
  Pr[E_S] ≥ |S|⁴ / (2·2^ω·φ(N)·Y²·(1 + log Y)).
Moreover, given R, the event E_S is decided by one Gauss reduction plus at most 2 membership tests per
multiple j·(u₀, v₀), j ≤ Y/max(u₀, v₀).

*Proof.*
1. Put n(R) = #{(u, v) ∈ S² : u ≡ Rv (mod N)}.
2. **First moment.** Let Ĝ be the group of characters of (Z/N)^× that are trivial on Q_N, so |Ĝ| = h.
   Orthogonality gives 1[x ∈ Q_N] = h^{−1}Σ_{χ∈Ĝ}χ(x). Then
     P := Σ_{R∈Q_N} n(R) = #{(u, v) ∈ S² : uv^{−1} ∈ Q_N} = h^{−1}Σ_χ |Σ_{u∈S}χ(u)|² ≥ |S|²/h.
   Every χ ∈ Ĝ is real, so χ̄(v) = χ(v) and χ(uv^{−1}) = χ(u)χ(v); each summand is a square.
3. **Second moment.** If (u, v), (u′, v′) ∈ S² both satisfy u ≡ Rv, then uv′ ≡ u′v (mod N). Also
   |uv′ − u′v| < Y² < N, so uv′ = u′v and the two vectors are proportional.
   - Write g = gcd(u, v) and (u₀, v₀) = (u, v)/g. The partners of (u, v) are the vectors j(u₀, v₀) with
     ju₀ ∈ [Y/2, Y].
   - Since u ≥ Y/2, we have u₀ ≥ Y/(2g), so there are at most Y/(2u₀) + 1 ≤ 2g partners.
   - Hence
       Σ_R n(R)² ≤ Σ_{(u,v)∈S², uv⁻¹∈Q_N} 2gcd(u, v) ≤ 2Σ_{g≥1} g·⌊Y/g⌋² ≤ 2Y²(1 + log Y).
4. **Cauchy–Schwarz.** #{R : n(R) > 0} ≥ P²/Σn(R)². Dividing by |Q_N| = φ(N)/h gives
     Pr[E_S] ≥ (|S|⁴/h²)·h/(φ(N)·2Y²(1 + log Y)).
5. **Algorithm.** Every lattice vector of {u ≡ Rv} with |u|, |v| ≤ Y < √N is a multiple of the primitive
   shortest vector: two independent ones would have 0 < |det| < N, but every determinant is a multiple of
   N. ∎

*Remarks.*
- No independence, no 1/ℓ or saddle-point law, and no character bound is used. The square-class
  condition is absorbed by positivity, and 2^ω cancels except for one factor.
- The fourth power is the price of the crude multiplicity bound. For smooth-type S one can use
  Ψ(x/d, y) ≪ d^{−α}Ψ(x, y) to recover |S|²·L^{o(1)}. Not needed here.
- E1x (§7) checks the bound exactly in 120 small cases. The true value is about |S|²/φ(N) for random S
  and about 2^ω|S|²/φ(N) when S ⊂ Q_N. So the bound is far from sharp but has the right structure.

**Theorem S′ (k-product sets) [proved].** Let G be a set of primes in a dyadic interval, coprime to N, with
|G| ≥ 4k². Let S be the products of exactly k distinct primes of G lying in a window [Y/2, Y] (Y² < N),
with |S| ≥ C(|G|, k)/(2k). Then Pr[E_S] ≥ |S|²/(c_k·2^ω·φ(N)) with c_k = 8k²·C(2k, k).

*Proof.* As in Theorem S, P ≥ |S|²/h. Split Σ_R n(R)² ≤ Σ_{(u,v)∈S²} #partners(u,v) according to
i = the number of shared primes, so (u, v) = (ju₀, jv₀) with j an i-product.
1. There are at most C(|G|,i)·C(|G|,k−i)² such pairs.
2. Each has at most C(|G|,i) partners j′(u₀, v₀).
3. So the i-th term is ≤ C(|G|,i)²C(|G|,k−i)² ≤ 2C(k,i)²C(|G|,k)² ≤ 8k²C(k,i)²|S|². The middle step uses
   C(n,i)C(n,k−i) = C(k,i)C(n,k)·∏_{t<i}(n−t)/(n−k+i−t) ≤ C(k,i)C(n,k)·e^{k²/(n−k)} ≤ √2·C(k,i)C(n,k) for n ≥ 4k².
4. Summing, Σ_i C(k,i)² = C(2k,k), so Σn(R)² ≤ c_k|S|². Then use Cauchy–Schwarz as before. ∎

This is the |S|² ("λ = 2") form. Theorem K(ii) improves accordingly, to
Pr ≥ c′_k·(g_A/(4 log N))^{2k}·2^{−ω}. The λ = 2 entries of Prop P (§7.4) are therefore *proved* bounds
for this design. E1x (random S) sits at about 0.9·|S|²/φ(N).

## 3. Theorem K: bounded-k splits reduce Hypothesis D to one averaged statement

**Setup.**
- Fix integers k ≥ 1 and odd d ≥ 3. Let X = ⌊√(N/2)⌋.
- Let j be maximal with 2^{kj} ≤ X, and put J = (2^{j−1}, 2^j], so ℓ ∈ J means ℓ = N^{1/(2k)+o(1)}.
- For a field A, let G_A ⊂ J be the set of primes ℓ ∤ N for which a fixed, deterministic,
  capped special-q search succeeds. The search depends only on (A, ℓ), with time cap T_SQ·L^{o(1)}. It
  succeeds when it finds a primitive (a, b) ∈ Λ_ℓ(m) in the box with (a − bm)/ℓ and F_A(a, b) both
  B-smooth.
- Let g_A = |G_A|/π_J with π_J = #(primes in J), and γ := E_A g_A = Pr_{A, ℓ∈J}[ℓ ∈ G_A].

**Hypothesis D_k(γ).** γ ≥ L^{−o(1)}, uniformly in N. The search box and the field law are those of §4.

**Theorem K [proved, given the setup].**
- (i) Pr_A[g_A ≥ γ/2] ≥ γ/2.
- (ii) Fix a field A, and let S_A be the products of exactly k distinct primes of G_A. Then for
  N ≥ N₀(k),
    Pr_R[R ≡ u/v with u, v ∈ S_A] ≥ c_k·(g_A/(4 log N))^{4k}·2^{−ω}/(1 + log N), with c_k > 0 depending
    only on k.
- (iii) Hence, in a field with g_A ≥ γ/2, a descent attempt succeeds after an expected L^{o(1)}·γ^{−4k}
  R-trials. Each trial costs L^{o(1)}: rigorous HECM (Lenstra–Pila–Pomerance) finds prime factors
  ≤ 2^j = N^{1/(2k)} in L_{2^j}[2/3] = L_N[1/3] time. Each trial with a correct shape adds at most 2k
  special-q searches.
- So the expected cost per descent relation is T_SQ·L^{o(1)} whenever γ ≥ L^{−o(1)}.

*Proof.*
- (i) Since 0 ≤ g_A ≤ 1, we have γ ≤ γ/2 + Pr[g_A ≥ γ/2].
- (ii)
  1. All products of k distinct primes of G_A lie in (2^{k(j−1)}, 2^{kj}] ⊂ (X·2^{−2k}, X]. By the
     pigeonhole principle, one of the 2k dyadic windows [Y/2, Y] with Y ≥ X·2^{−2k} holds at least
     C(|G_A|, k)/(2k) of them.
  2. Distinct prime sets give distinct products, so |S_A ∩ [Y/2, Y]| ≥ (|G_A|/k)^k/(2k).
  3. Use π_J ≥ 2^{j−1}/(2(j−1)log 2) for j ≥ j₀ (Rosser–Schoenfeld) and Y ≤ 2^{kj}. Then
     |S|/Y ≥ (g_A/(4kj log 2))^k/(2k), and 4kj log 2 ≤ 2 log N.
  4. Apply Theorem S with Y² ≥ X²·2^{−4k} ≥ N·2^{−4k−2}.
- (iii) Average over R. The bounded factors (log N)^{O(k)}, 2^{O(k)} and 2^ω are L^{o(1)}, by Lemma Ω.
  The shape test costs L_N[1/3]·L^{o(1)}. Conditional on u₀v₀ having the right shape, at most 2k
  special-q searches are run. Pr[shape]/Pr[success] ≤ (log N)^{O(k)}g_A^{−4k}2^{2ω} = L^{o(1)} (upper-bound Pr[shape] by 2·2^ω|S_J|²/φ(N), S_J = all k-products from J), so
  wasted searches cost only L^{o(1)} on average. ∎

**Corollary K [proved implication, modulo the LV §6 adaptation that lemmaM Corollary DD also leaves
unwritten].** Assume D_k(γ). Then every N is factored in expected time
L[1/2, √(2θ(d,k)) + o(1)], where
  θ(d,k) = 2s + (d−1)/(4k) + max(0, s − (d−2)/(4k)), s = 1/(d+1).
- *Algorithm.*
  1. Field search with restarts: try fields; in each, run a pilot of L^{o(1)} capped descents and keep
     the field if at least half succeed. By (i) and Markov, a kept field arrives after L^{o(1)} tries.
  2. Collect C·L^{o(1)} pure-descent rows, where C is the column count (factor base plus LV §6
     characters).
  3. Linear algebra, then Theorem DD.
- *Cost.* Take B = L[1/2, a]. The relation cost is B·T_SQ = L[1/2, a + θ/(2a)] and linear algebra is
  L[1/2, 2a]. At a = √(θ/2) both equal L[1/2, √(2θ)].
- *Sizes.*
  - ‖(a, b)‖ ≤ N^{1/(4k)}S, where S² = T_SQ = L^{O(1)}.
  - The rational cofactor is (a − bm)/ℓ = N^{s − 1/(4k) + o(1)}.
  - F(a, b) = N^{s + h + d/(4k) + o(1)}, where H = N^h is the coefficient height of R_A.
  - The level-½ condition |a − bm| ≤ |F|^{1/2−η} is used by FM (§4). It needs h > s − (d−2)/(4k) + 2η.
  - Adding the rational and algebraic sizes gives θ(d,k) + O(η).
- θ(d,k) ↓ 3s as k → ∞. For instance, d = 9: k = 50 gives c = 0.781; k → ∞ gives 0.775.

**Why bounded k matters [proved; arithmetic].** With k large primes per split,
Pr[all in G_A] ≈ g_A^{k}.
- lemmaM uses B1 = exp((log N)^{3/4}), so k ≍ (log N)^{1/4}. Then g_A must be exp(−o((log N)^{1/4})),
  so every second-moment loss would have to be that small. That is finer than anything LV's methods
  give: their "goodness" carries an L^{o(1)} exponent loss.
- With k fixed, g_A = L^{−o(1)} suffices.
- The price is that the special-q primes are N^{1/(2k)} rather than N^{o(1)}. This costs only the
  O(1/k) term in θ(d,k), and the shape test is still L_N[1/3] by HECM.

## 4. From D_k to a first and a second moment

**Field and box law.**
- m is uniform in [N^s/2, N^s].
- f = f̂_{N,m} + (x − m)R_A, where f̂ is LV's non-monic base-m expansion and R_A has degree d − 1 with
  i.i.d. coefficients c uniform in [−H, H], H = N^h.
- ℓ is uniform in the primes of J.
- The box is {(a, b) ∈ Λ_ℓ(m) : |a|, |b| ≤ √ℓ·S}, primitive with b > 0. Irregular lattices (λ₁ ≤ √ℓ/S)
  are excluded; they occur with probability O(S⁻²) over m.
- N_ℓ(A) is the number of box points with (a − bm)/ℓ and F_A(a, b) both B-smooth.
- For primitive (a, b), F_A(a, b) = F̂(a, b) + t·⟨c, v(a, b)⟩, with t = a − bm,
  v(a, b) = (a^{d−1−i}b^i)_{0≤i<d}, and F_A(a, b) ≡ N b^d (mod t).

**Prop PZ [proved].** γ ≥ (E N)²/E[N²] = 1/(1/E N + PC), where PC := E[N(N−1)]/(E N)² and every E is
over (A, ℓ). So if E N ≥ 1 (hypothesis FM) and PC ≤ L^{o(1)} (hypothesis PC), then D_k(γ) holds with
γ ≥ L^{−o(1)}. ∎

### 4.1 FM, the first moment [sketch; inputs are proved theorems, but the transfer is not written out]
FM says E N ≥ L^{−o(1)}·S²·ρ(u_rat)·ρ(u_F). Choose S² = T_SQ = L^{o(1)}/(ρ(u_rat)ρ(u_F)).
1. **Rational side (proved input).**
   - Fix ℓ and primitive (a, b). As m runs over m ≡ a/b (mod ℓ), the cofactor s = (a − bm)/ℓ runs over an
     arithmetic progression of common difference b ≤ √ℓS, with M/ℓ terms of size ≍ bM/ℓ.
   - Harper's BV theorem for smooth numbers applies to these moduli b, which form a *dense* family, so
     the log-power saving suffices. It needs √ℓS ≤ (bM/ℓ)^{1/2−η}, i.e. roughly
     k ≥ (3/4)(d+1)(1+O(η)).
2. **Algebraic side (LV §5 + LV Lemma 8.5, transplanted).**
   - LV Lemmas 5.13–5.20 show that F_A(a, b) is close to uniform on long progressions mod t.
   - LV's goodness (Def. 5.12) then gives Pr_c[F smooth] ≥ ρ(u_F)^{1+o(1)} for every *good* modulus t.
   - It remains to show that bad moduli are an o(1) proportion of the family 𝒯 = {ℓs : ℓ ∈ J, s B-smooth},
     weighted by the induced law. This is a restricted BV for 𝒯 with relative saving.
   - **Lemma 8.5′ [sketch]** is LV's Lemma 8.5 with y-smooth moduli replaced by 𝒯. Take LV's conductor
     decomposition, with m₀ = min(B^η, e^{η√log x}) and M₀ = x^η, and choose η < 1/(2k(s + h + 1)).
     - Every conductor divisible by ℓ is ≥ ℓ ≥ M₀. Those characters fall in LV's large-conductor range,
       where Harper's Fact 8.11 bounds the sum over *all* moduli by Q√Ψ·log^{7/2}x.
     - Smaller conductors r divide the B-smooth part s. The counts Σ_{q∈𝒯, r|q}1/φ(q) equal LV's sums
       (Claims 8.14, 8.16), multiplied by the common factor Σ_{ℓ∈J}1/ℓ. That factor also multiplies
       |𝒯|, so the relative saving e^{−cu/log²u} + B^{−c} is unchanged.
     - The large-conductor term needs Q ≤ √Ψ(x, B)·δ_𝒯·L^{−ε}, where δ_𝒯 is the density of 𝒯. This holds
       with polynomial room, because t ≤ |F|^{1/2−η}.
   - LV's proof uses one y for both the moduli and the counted numbers; here both are B, plus the prime ℓ.
3. **Status.** FM is a parameter transfer of LV §5 and §8; we expect no new idea to be needed. It is
   labelled *sketch* until written out. One point to check: LV's proof of Lemma 5.21 compares bad
   B-smooth moduli against Ψ(Q_max, B′) rather than Ψ(Q_max, B). This is harmless only when β = β′ (our
   design has a single B), and a referee should confirm it.

### 4.2 PC, the second moment [open]
E[N(N−1)] sums over pairs of distinct box points (a, b), (a′, b′) of one special-q lattice. It splits as
follows.
- **PC_rat.** The rational cofactors satisfy b′s − bs′ = Δ/ℓ, where Δ = ab′ − a′b and |Δ/ℓ| ≤ 2S². This
  identity holds *for every m*. So as m varies the pair (s, s′) moves along one line, and the rational
  part of E[N(N−1)] is an average, over coefficients b, b′ ≤ √ℓS and shifts h = Δ/ℓ, of
  #{k : s₀ − bk and s₀′ − b′k both B-smooth}. This is a binary additive problem; b = b′ = 1 is the
  "smooth pairs (n, n + h)" problem. We know no estimate, even an upper bound within L^{o(1)} of
  ρ(u)², at smoothness L[1/2]. See §8 for the literature check.
- **PC_alg.** Given the rational data, PC_alg = Pr_c[F, F′ both smooth]/(Pr_c[F]Pr_c[F′]).

**Prop G (geometric barrier) [proved].** Let (a, b) ≠ ±(a′, b′) be primitive points of Λ_ℓ(m), and let
Γ = {(⟨c, v⟩, ⟨c, v′⟩) : c ∈ Z^d}, with v = v(a, b) and v′ = v(a′, b′).
- (i) Every 2×2 minor of [v; v′] is divisible by Δ. Hence ℓ | Δ | D := [Z² : Γ].
- (ii) For c uniform on any box, the law of (F, F′) is supported on a coset of the lattice
  Γ″ = {(tX, t′X′) : (X, X′) ∈ Γ}, which has index D in tZ × t′Z. The marginals live on the full
  progressions mod t and mod t′.
- (iii) Consequently, any bound of the form Pr[F, F′ ∈ 𝒮] ≤ sup_z μ(z)·#{(n, n′) ∈ 𝒮² in the product
  progression box} exceeds Pr[F ∈ 𝒮]·Pr[F′ ∈ 𝒮] by a factor ≥ D/C_d ≥ ℓ/C_d = N^{1/(2k)}/C_d. That
  holds whenever 𝒮 has the product density ρ_Fρ_{F′} in the product box. Such bounds include the
  continuous Ball-slicing bound plus lattice corrections (lemmaM §3).
- *Proof.*
  - (i): For i < j, v_iv′_j − v_jv′_i = a^{d−1−j}a′^{d−1−j}b^ib′^i((ab′)^{j−i} − (a′b)^{j−i}), and
    x^n − y^n is divisible by x − y. Also ℓ | Δ, because a ≡ bm and a′ ≡ b′m (mod ℓ).
  - (ii): By definition.
  - (iii): μ is carried by Γ″-points. These make up a 1/D share of the product-lattice points in the
    support box, so sup μ ≥ D/(#product points) up to the box-spline shape constant C_d. ∎
- *Meaning.* PC_alg needs to know that smooth pairs equidistribute over the D cosets of Γ″. That is
  equidistribution of B-smooth numbers modulo q = tD in each coordinate. The ℓ-part of D is invisible to
  smoothness, since F ≡ Nb^d (mod ℓ) is already fixed by the marginal, but it still has to be proved.

**Hypothesis RBDH (restricted Barban–Davenport–Halberstam) [conjectured].** Let 𝒬 be the multiset of
moduli q = t·D(pair) over pairs of the special-q boxes, weighted by multiplicity w(q). Then for x ≍ |F|,
  Σ_q w(q) Σ_{(r,q)=1} |Ψ(x, B; q, r) − Ψ_q(x, B)/φ(q)|² ≤ o(1)·Σ_q w(q)Ψ_q(x, B)²/φ(q).
- **PC_alg ⟸ RBDH + a two-dimensional version of LV Lemma 5.13 [sketch].**
  - The 2D lemma needs small kernel vectors realizing every Γ-shift in a box; it is plausible for d ≥ 3.
  - Decompose the Γ″-coset into D cosets of tDZ × t′DZ.
  - The projection of each coset onto a coordinate is surjective mod D, because v is primitive.
  - Write each count as main term plus error. Cross terms are handled by restricted BV, and the
    error × error term by Cauchy–Schwarz and RBDH.
  - Nearly parallel pairs have sin∠(v, v′) ≍ |Δ|/(‖(a,b)‖‖(a′,b′)‖) ≥ 1/S². A dyadic decomposition in
    the angle costs only log S = L^{o(1)}.
- **Why RBDH is plausible.**
  - Harper's BDH (his Theorem 2) has level Ψ(x, y), and all our q ≤ x^{1/2+O(1/k)}.
  - LV state that their method "naturally provides" the BDH analogue for smooth moduli.
  - Our family is ℓ²·(B-smooth)·(Δ′D″) with Δ′ ≤ 2S² arbitrary, so it is denser than LV's.
- **Why RBDH is not proved.**
  - The family is not of LV's shape.
  - The weights w(q) come from a 4-variable geometric count.
  - The diagonal term QΨ must beat δ_𝒬Ψ²/Q.
  - None of this has been carried out. Restricted-moduli BDH is not in the literature we found (§8).

**Assessment.** PC_alg looks within reach of the Harper/LV machinery; this is an assessment, not a
claim. PC_rat does not: it is a correlation of smooth values of two linear forms in one variable. A
possible lever, not worked out, is a *randomized rational side* g = g₁x − g₀ with (g₁, g₀) random
(Kleinjung-style selection, Res(f, g) = N). That gives the rational side 2D randomness, which turns
PC_rat into a PC_alg-type (RBDH) statement. It must be checked against Kleinjung's residue constraints on
g₀ mod g₁.

## 5. Design V: non-fruitful LV relations plus a few fruitful Vallée relations

**Observation [proved; elementary linear algebra + Dixon fibre symmetry].**
- Let 𝒞 be any multiset of relations whose coins are independent of the hidden roots, for instance LV
  relations, which need not be fruitful.
- Let 𝒱 = {(x_i, e_i)} be relations x_i² ≡ ∏p^{e_i} (mod N) whose roots are uniform on their fibre
  x_i·μ₂ given everything else. Examples: Dixon's uniform x, or Vallée's quasi-uniform small-square
  sampler, by its x ↦ εx invariance.
- Then the number of independent dependencies that contain a 𝒱-row is
  |𝒱| − (rank(𝒞 ∪ 𝒱) − rank 𝒞).
- Each such dependency is fruitful with probability ≥ (1 − 2^{1−ω})·(fibre-bias factor). This is
  Theorem DD with descent replaced by direct sampling.
- So it suffices that some Vallée vector lies in the row space of 𝒞, after appending zero algebraic and
  character coordinates.

**Quasi-uniform samplers [proved; standard convolution argument].**
- Vallée's sampler is only *quasi*-uniform: each x in the target set has probability within a constant
  ratio κ = c₂/c₁ of uniform. Given z = x² mod N, each fibre element εx then has probability
  ≥ 1/(κ2^ω). With a single Vallée row this can fail to give fruitfulness when κ ≥ 2^{ω−1}.
- Fix: sum t independent dependencies D₁, …, D_t, each containing its own Vallée row.
  - The congruence ratio of the sum is ±∏_i ε_iZ_i, with the ε_i independent given everything else.
  - Each ε_i has min-probability ≥ 1/(κ2^ω) on μ₂ ≅ (Z/2)^ω.
  - The t-fold convolution is within (1 − 1/κ)^t of uniform in total variation.
- So t = O(κ(ω + 1)) = L^{o(1)} Vallée rows give fruitfulness 1 − 2^{1−ω} − o(1).

**Hypothesis RK.** For LV fields at fixed odd d and B = L[1/2, a], with C·L^{o(1)} relations, a Vallée
vector lies in the GF(2) row space with probability ≥ L^{−o(1)}.
- Necessary structure [proved]: the row space ∩ (rational-only vectors) lies inside
  Λ/2 := {v : ∏p^{v_p} ∈ Q_N}, as soon as character-trivial even combinations are true squares in K
  (LV §6). Indeed ∏(a − bm) ≡ φ(γ)².
- RK asserts the reverse inclusion, up to rarely occurring primes.

**Theorem V [conditional].** Assume RK, and LV relation supply at fixed odd d ≥ 9 and B = L[1/2, 1/√6].
The supply needs only a first moment, by LV's stochastic deepening; it is the FM-type transfer of §4.1
without special-q. Assume also the LV §6 adaptation. Then every N is factored in expected time
L[1/2, √(2/3) + o(1)].
- *Cost.* Linear algebra costs L[1/2, 2a]. LV relations cost L[1/2, a + θ/(2a)] with θ = 3/(d+1) ≤ 1/3.
  Vallée relations each cost L[1/2, (2/3)/(2a)], with L^{o(1)} of them. All three are ≤ L[1/2, 2/√6]
  at a = 1/√6.
- Vallée's sampler is rigorous: Math. Comp. 56 (1991). With uniform Dixon x the same argument gives
  only c = 1.
- *Proof.* Combine the observation with the cost count. ∎

**Remarks.**
- RK needs no fruitfulness, no hidden root in 𝒞, and no special-q. It is a statement about the GF(2)
  rank of NFS relation matrices: the relations "fill" Λ/2.
- Its natural proof route is relation-lattice generation in the Hafner–McCurley/Buchmann style. That
  route needs random-exponent (descent-type) relations, so RK is not obviously easier than PC. It is
  nonetheless a *different* missing ingredient, and E3 supports it.
- Any rigorous sampler of small squares with θ_V < 2/3 would give c = √θ_V.

## 6. Barriers and what does not work

- **(B1) The log-saving barrier [known; LV §8].** A BV or BDH bound whose saving is
  (log x)^{−A}, or e^{−cu/log²u}, relative to the sum over *all* moduli says nothing about a moduli
  family of density L^{−c}. Every first or second moment here lives on such a family. So *restricted*
  versions, in the style of LV Lemma 8.5, are necessary. This is why FM must transplant LV §8 rather
  than cite Harper.
- **(B2) The geometric barrier (Prop G) [proved].** Pair laws in one special-q lattice are carried by an
  index-D ≥ ℓ sublattice. Slicing and density bounds alone lose at least ℓ = N^{1/(2k)}, so arithmetic
  equidistribution mod tD is unavoidable.
- **(B3) The precision barrier for unbounded k [proved arithmetic, §3].** With k ≍ (log N)^{1/4}
  targets per split, every averaged estimate must lose at most exp(o((log N)^{1/4})). That is finer than
  LV-type goodness, whose losses are L^{o(1)} in the exponent. Fixed k removes this.
- **(B4) What Theorem K does *not* remove.**
  - The good-field set and the target set are still field-dependent. That is harmless: Theorem S takes
    any set.
  - The first moment still has moving residues Nb^d mod t. LV's goodness handles them because it is
    uniform in residues: max_r, "for all (a, r) = 1".
  - Pascadi's level 5/8 is *not* usable. It is for a fixed residue, and our residues move with b.
- **(B5) The QS-style descent [assessment].** The special-q relation x² − N = ℓc has no randomization to
  average over. FM and PC would have to hold for the single polynomial x² − N at smoothness L[1/2].
  Known results on smooth values of quadratic polynomials (Dartyge; Harper; Balog–Blomer–Dartyge–
  Tenenbaum for binary forms) need smoothness x^{θ} with θ fixed, far above L[1/2] (§8). Not pursued.
- **(B6) No L[1/3] [lemmaM (c), unchanged].** Theorem K needs a shape test for N^{1/(2k)}-smooth √N-size
  integers, which HECM does in L[1/3]. The descent itself is one level, and lemmaM's α ≥ 2/5 argument
  applies to the growing-degree regime.
- **(B7) RK and fruitfulness are different statements [observation].** RK says the NFS relations fill
  Λ/2. Conjecture F says their *dependencies* carry a nonzero β-label. RK + Vallée ⇒ fruitful
  dependencies exist, via the Vallée row's hidden root. Conjecture F concerns pure-NFS dependencies,
  which Design V never needs.

## 7. Experiments (pre-registered; exact; toy; correctness and plausibility checks only)

Commands are in `notes/hypD_prereg.md`. Analysis: `.venv/bin/python -m src.hypD.analyze {e1,e2,e3} …`.
The E1x summary was computed inline from `experiments/raw/hypD_e1x.jsonl`; the script is in §7.1.

**Disclosed deviations.**
1. The frozen E1 code crashed on an ACCIDENTAL gcd(r, N) = 13127 in a 40-bit `tri` instance. A random
   r hit a 13-bit prime factor; the null rate is about trials·Σ1/p. The fix resamples r and counts such
   events; E1 was rerun.
2. E1 as registered had model-expected pooled hits of about 0.3, so P2 is **uninformative** (a design
   error). It was replaced by the exact E1x, registered in Addendum A before E1x ran.

### 7.1 E1 / E1x (Theorem S)
- **E1 (80 records).**
  - P1: 0/80 below the bound. Trivially so: the bound is about 10⁻¹⁶.
  - P2 is uninformative: 1 hit against a model expectation of 0.1.
  - 16 ACCIDENTAL gcd(r, N) > 1 events occurred in 320 000 r-draws. The r-stream is seeded per (N, seed),
    so the 4 (k, g) records of one N reuse the same draws: that is 4 distinct events, against a null
    expectation of about 8 (Σ1/p over the 13–16-bit factors of the tri instances). STRUCTURAL rate 0.
- **E1x (120 exact records).**
  - **P1′ holds:** 0/120 violations. The exact rate exceeds the bound by at least 7.5·10⁶.
  - **P2′ holds:** in 30/30 random-mode records with |S| ≥ 200, exact/model ∈ [0.849, 0.943]
    (mean 0.908). The shortfall is consistent with non-coprime pairs collapsing.
  - **P2″ was vacuous:** no qr-mode record had |S| ≥ 200. Descriptively, the largest qr records
    (|S| = 61–186) give exact/(2^ω·model) = 0.73–0.83, i.e. exact/model ≈ 2.9–3.3 at ω = 2. This shows
    the 2^ω positivity excess.
- Inline summary:
  `python - <<…` computing min(exact/bound), and exact/model over records with |S| ≥ 200 by mode.

### 7.2 E2 (special-q counts in toy LV-style fields)
- Grid: 64 and 72 bits, blum/tri, 4 seeds, 60 fields × 40 primes ℓ ∈ (2^14, 2^15], B = 8192, box 65×33.

| check | registered | result |
|---|---|---|
| P3 PC_alg | [0.67, 2.0] | 1.076–1.212 in 16/16 (PC 1.080–1.216) — holds |
| P4 Poisson-with-field-quality γ | \|diff\| ≤ 0.07 | max 0.009 — holds |
| P5a overdispersion of g_A | [0.5, 4] | 2.24–4.21; **fails in 1/16** (blum 64 s1: 4.21) |
| P5b Pr_A[g_A ≥ γ/2] | ≥ 0.8 | 0.95–1.00 — holds |
| P5c same-field pair ratio | [0.9, 1.3] | 1.016–1.066 — holds |

- μ = 1.41–1.64 at 64 bits and 0.66–0.80 at 72 bits; γ = 0.72–0.78 and 0.47–0.53.
- Reading:
  - Pair correlations within one special-q box are mildly positive, ≈ +15%. This is consistent with
    Prop G's mod-ℓ correlation being invisible to smoothness, plus small-prime effects mod Δ′.
  - Field quality varies; the overdispersion of about 3 is the usual root-property (Murphy-α) effect.
    It does not create complementary target types. The same-field pair ratio stays ≈ 1, and almost
    every field has g_A ≥ γ/2.
  - The P5a miss (4.21 against a bound of 4) is an underestimate of field-quality variance on my part,
    not a model failure relevant to Theorem K. Theorem K only needs Pr_A[g_A ≥ γ/2].

### 7.3 E3 (RK in a toy NFS)
- Grid: 48 bits (B = 1200) and 56 bits (B = 3000), blum/tri, 5 seeds each. All 20 instances were valid
  (rows/cols 1.29–1.5); rank/cols was 0.92–0.97.
- **P6 holds:** 360/400 Dixon vectors lie in the row space (90%), 15–20 of 20 per instance.
- **P7a holds:** 0/5098 single-prime vectors with p ∉ Q_N are in the span.
- **P7b holds:** 787/787 with p ∈ Q_N and ≥ 3 odd occurrences are in the span.
- Reading:
  - Span ∩ rational = Λ/2 restricted to well-covered primes, exactly as RK predicts.
  - The Dixon failures plausibly come from rarely occurring primes. This was not separately registered,
    so it is a hypothesis.
  - The toy skips index-divisor and projective pairs.
  - No square roots or fruitfulness were computed here; DD-type fruitfulness was tested in lemmaM E1.

### 7.4 What the best *known* pair-correlation bound buys (Prop P) [proved implication; numerics exact to 3 digits]
- **Input.** The sharpest published upper bound for two linear forms is Pascadi (arXiv 2505.00653,
  Cor. 7.2), improving de la Bretèche–Drappeau (JEMS 2020, Thm 4.1, exponent 3/5):
  #{n ∈ [x, 2x] : an + b ∈ S(y₁), cn + d ∈ S(y₂)} ≪ Ψ(x, y₁)·ρ(u₂)^{5/8−ε}.
  - It requires |a|, |b|, |c|, |d| ≤ x^δ and (log x)^C ≤ y ≤ x^{1/C}.
  - This applies to PC_rat: the coefficients are b, b′ ≤ N^{1/(4k)}S and the shifts are ≤ N^{O(1/k)}. It
    gives PC_rat ≤ ρ(u_rat)^{−3/8−ε}, *uniformly*, so no averaging is needed. Independence would need
    ρ^{−o(1)}.
- **Consequence.** Assume FM and PC_alg (the RBDH part) but use only this proved bound for PC_rat. Then
  γ ≥ ρ(u_rat)^{3/8+ε}/L^{o(1)}, and each descent costs γ^{−4k} more.
- Optimizing over odd d and k (`θ_eff = θ(d,k) + 4k·(3/8)·max(0, s − 1/(4k))`, then c = √(2θ_eff))
  gives c ≥ 1.30 with Theorem S's fourth power (λ = 4) and **c ≥ 1.17 with Theorem S′ (λ = 2, proved)**. As d → ∞, θ_eff → √z − z/4 with z = λ·(loss exponent).
  - **Threshold.** c < 1 needs √z − z/4 < 1/2, i.e. z < 0.343. That means a pair-correlation loss
  exponent below ≈ 0.086 (λ = 4) or ≈ 0.17 (λ = 2), against Pascadi's 0.375 and the conjectured 0.
- So **known tools, even granting FM and PC_alg, do not beat Lenstra–Pomerance through this design.** The
  missing input is quantified: smooth pairs on lines must be within ρ^{−0.17} of independence.

## 8. Literature (Opus literature subagent, read from primary PDFs unless marked; 2026-10-19)
- **Harper, arXiv 1208.5992** (preprint; no journal version found).
  - Thm 1 (BV): for log^K x ≤ y ≤ x and Q ≤ √Ψ, the error is ≪ Ψ(e^{−cu/log²(u+1)} + y^{−c}) + Q√Ψ log^{7/2}x.
  - Thm 2 (BDH): for Q ≤ Ψ, the error is ≪ Ψ²(e^{−cu/log²u} + y^{−c}) + QΨ.
  - Neither has a restricted-moduli version.
- **de la Bretèche–Drappeau, JEMS 22 (2020), Thm 2.1 (2.2).** This is a published BV for y₂-smooth moduli
  q ∼ Q with Q ≤ x^{1/2−ε}.
  - It includes multiples q₀q of a fixed y₂-smooth q₀ ≤ x^δ, has max over residues, and saves
    Ψ(Q, y₂)/(φ(q₀)Q)·e^{O(u₂)}(log x)^{−A}.
  - The authors credit LV for the smooth-moduli modification.
  - **Bearing on FM.** This is LV Lemma 8.5 in clean form, and *per fixed q₀*. With q₀ = ℓ it would give
    FM per target ℓ, not only on average, except that the theorem requires P⁺(q₀) ≤ y₂. With y₂ = B < ℓ
    that fails. Whether the condition can be dropped for a single prime q₀ = ℓ ≤ x^δ is the one point to
    check; this is an assessment. If it can, FM becomes a citation. The ℓ ≤ x^δ condition holds for
    k ≥ k₀(δ).
- **Drappeau, Compositio 2015**: level 3/5, with a fixed residue a₁ā₂ and no max over residues.
  **Pascadi**: level 66/107, and 5/8 unconditionally (arXiv 2505.00653), also with a fixed residue.
  Neither handles our moving residues Nb^d, as (B4) notes.
- **Pair correlations.**
  - de la Bretèche–Drappeau Thm 4.1 (ρ^{3/5}); Pascadi Cor. 7.2 (ρ^{5/8}).
  - Pascadi Cor. 1.6: #{n ≤ x : n, n + 1 ∈ S(y)} ≪ xρ(u)^{1+5/8−ε}.
  - Goldfeld–Yamagishi-type results show that at *extremely* small y the smoothness of n and n + 1 is
    genuinely not independent. So PC_rat at precision L^{o(1)} is not even formally guaranteed; at
    y = L[1/2] it is expected but open.
  - Harper's restriction theory (Compositio 2016) handles ternary problems (a + b = c) at log^K x ≤ y,
    not binary ones.
- **Hyperelliptic curve method.** Lenstra–Pila–Pomerance I (Phil. Trans. 1993), Cor. 1.2: for every
  n ≥ n₀, all prime factors ≤ v are found with probability ≥ ½ in L_v[2/3, c₀](log n)³.
  - **Caveat:** Part II (Proc. LMS 2002) defers details to a Part III that apparently was never
    written. The rigorous assembly is not in print. lemmaM and LV (Cor. 3.19) rely on the same method.
- **Vallée, Math. Comp. 56 (1991).** The primary text could not be retrieved (AMS block). The quasi-uniform
  sampling of x with x² mod N = O(N^{2/3}) and the rigorous L[1/2, √(4/3)] come from secondary sources
  (LV, STOC'89 abstract). **Theorem V's Vallée input is therefore not verified by us.**
- **Rigorous factoring and discrete logs, 2018–2026.**
  - Nothing found below c = 1 for factoring, or below √2 for discrete logs in prime fields.
  - Kleinjung–Wesolowski (fixed characteristic, quasi-polynomial) is rigorous but does not transfer.
  - de Boer–Pellet-Mary–Wesolowski (arXiv 2512.01588, ERH): rigorous class-group sampling toolkit.
  - Stange 2022 is heuristic (Hypothesis 3.1).
  - LV follow-ups: van Leeuwen's thesis, Le Gluher–Spaenlehauer–Thomé (heuristic), and an announced
    Balister–Bollobás–Morris–Lee paper on the 2-core.
- **Rank.** Balister–Bollobás–Morris (Annals 2018) prove that the GF(2) parity matrix of fewer than
  (e^{−γ} − ε)J *uniform random* integers has full row rank whp. This is the closest rigorous rank
  statement to RK, but it is for uniform integers rather than NFS relations, and it concerns rows, not
  the span of Λ/2. We found no span or relation-lattice result for sieve relations.

## 9. Reviews
- **GPT-6.1-sol: PENDING.** Both the blind round `codex/hypD_1.md` and the adversarial round
  `codex/hypD_2.md` are written. They could not run because the Codex login expired
  (refresh_token_expired, 401, 2026-10-19).
- **Opus literature subagent.** It supplied §8: the exact Harper, de la Bretèche–Drappeau, Pascadi, LPP
  and BBM statements, and the Vallée caveat.
- **Independent Opus critics (proofs; prior art): launched, NOT YET RECEIVED at hand-back.** Every
  'proved' label in this note is therefore *self-checked only* and pending both critics and GPT-6.1-sol.

## 10. Candidate cards
```text
ID: KD-1 bounded-k Dixon descent in LV fields (Theorem K + Theorem S + DD).
Target: every N (after Lemma Omega).
Operation: R = r^2; split R = u/v with u, v products of exactly k primes from J ~ N^{1/(2k)}, all in the
  field's descendable set G_A; one special-q relation per prime; DD gives fruitfulness.
Proved: Theorem S (any S), Theorem K (D_k => complete descents in good fields at cost T_SQ L^{o(1)}),
  Prop PZ (D_k <= FM + PC), Prop G (geometric barrier), Prop P (known PC_rat bound gives only c ~ 1.17).
New lemma required: PC (pair correlation, L^{o(1)} precision) = PC_rat (binary smooth correlation on
  lines; best known loss rho^{-3/8}) + PC_alg (<= RBDH, conjectured); FM (sketch; ~ dlB-D 2020 Thm 2.1).
Cost: L[1/2, sqrt(2 theta(d,k))] -> sqrt(6/(d+1)) (d = 9: 0.775) conditional; ~1.17 with the known PC_rat bound.
Status: conditional; experimentally supported in toys (E2: PC_alg 1.08-1.21, no complementary types).

ID: V-1 LV rows + Vallée rows (Design V).
Proved: dependency count and fibre-symmetry fruitfulness (observation); cost sqrt(2/3).
New lemma required: RK (Vallée vectors lie in the NFS row space with prob L^{-o(1)}); LV supply at L[1/2].
Status: conditional; E3: 360/400 Dixon vectors in span; span ∩ rational = Lambda/2 on well-covered primes.
Vallée's theorem not re-verified from the primary source.

ID: G-1 slicing route to PC_alg. Status: refuted for this model (Prop G: loss >= l = N^{1/(2k)}).
ID: QS-D QS-style descent. Status: not pursued (no randomization; fixed-polynomial smooth values open).
```
