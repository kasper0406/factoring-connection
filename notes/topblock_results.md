# topblock: do the ladder's structured windows beat known-top-bits factoring?

Agent `topblock`, phase 18, 2026-10-20.

**Files.**
- Pre-registration: `notes/topblock_prereg.md`, commit `d661588`, made before any experiment.
- Code: `src/topblock/{core,lattices,run_experiments}.py`. Tests: `tests/test_topblock.py`.
- Raw data: `experiments/raw/topblock_x{1,2,2post,3,4,5,6}.jsonl`. Summary:
  `experiments/analysis/topblock_summary.json`. Logs: `experiments/logs/topblock_*.log`.
- GPT-6.1-sol: blind round `codex/topblock_1.{md,out}`; adversarial round `codex/topblock_2.{md,out}`.

**Labels.** KNOWN (with reference), PROVED (proof here), experimentally supported, HEURISTIC,
CONJECTURED, OPEN, UNVERIFIED. Novelty is claimed at most as "no prior formulation found by us".

## 0. Verdict (plain language)

**The structured windows give nothing beyond plain known-top-bits factoring.** The floor-sum
ladder at degree j is, exactly, "guess the top (1/(j+2))·log₂N bits of p, then solve a
partial-key-exposure instance". No lattice built from its window conditions does better than
Coppersmith's univariate lattice. Its limit is the Coppersmith barrier. With every block solver known
to us, the ladder costs N^{1/4+o(1)} for every j ≥ 2.

What is established:

1. **The equivalence (§1, PROVED).**
   - *Windows are near-divisor sets.* A top-block window holds exactly the d = x + h for which N − k
     is a multiple of d with |k| = O(B), in an exactly computable k-range (identity (W)).
   - *The Taylor structure is the trivial relation.* For a divisor d with cofactor q, the window
     residue is r = −q(−h)^{j+1}. It lies on the curve F_j = (x+h)r + N(−h)^{j+1} = 0, which is the
     identity N = q(x+h) multiplied out (T2). Every degree, block and shift gives such a curve.
   - *Oracle ⇒ known-bits search:* bisection with shrinking windows uses a proved number of calls
     (Prop LOC). Codex's degree-(j+1) query removes all false candidates, so a polylog
     degree-(j+1) oracle gives O(log H) calls.
   - *The converse needs more:* answering window counts takes O(window width·x) calls of a general
     divisors-in-interval finder on the integers N − k (Prop IV). A known-bits solver for N alone does
     not suffice. The full oracle is therefore a **stronger** primitive than the ladder needs, never
     a weaker one.
   - *Hence (Theorem EQ):* any ladder implementation at exponent 1/(j+2) + e gives known-bits
     factoring at the same per-block cost (up to the proved candidate factor). Conversely,
     known-bits factoring plugged into the ladder's outer loop gives the same exponent. The windows
     carry no extra information about p. The extra content of a window count is the near-divisor
     multiset, which concerns the other integers N − k.
2. **Multi-window / multi-block lattices (§3).**
   - **Prop T (transfer theorem, PROVED).** Take any lattice built, Coppersmith-style, from the
     divisor congruence x + h ≡ 0 (mod p) and any number of window congruences, from any blocks,
     Taylor degrees and centres. The substitution r_b ↦ −N(−h_b)^{j_b+1}/(x+h) maps it into
     ∏D_b^{m2_b}·(N, x+h)^{m1+T} ⊂ Z[h], the ordinary univariate Coppersmith ideal. Its kernel
     consists of multiples of the curve relations F_b, which vanish on the whole curve and locate
     nothing.
   - A successful mixed lattice therefore yields a univariate polynomial in the standard ideal whose
     norm is within an explicit 2^{O(d)} factor of the threshold: an o(1) loss in the exponent.
   - Windows add no range beyond what the univariate ideal itself supports.
3. **Why naive counts mislead (PROVED + measured).** F·N^{m1−1}D^{m2−1} lies in the mixed lattice and
   is Howgrave-Graham-short up to X ≈ x·(x/N)^{m1/(j+1)}, i.e. N^{1/3} for j = 2 and m1 = 1. So mixed
   lattices show "good" short vectors far beyond N^{1/4}, but every one of them is trivial.
   - E.g. at 128 bits the good vectors reach log₂X = 41.5, while the actual range is 24.5.
4. **Experiments (pre-registered).**
   - **X3:** in all 54 cells, the mixed range equals the transfer target U_img to within 0.5 bit. It
     is 1.0–4.5 bits *below* a univariate lattice of the same dimension. The trivial-vector extent
     matches the F prediction exactly (difference 0.0 bits in all 18 m1 = 1 cells).
   - **X4 (two windows: adjacent blocks, or two degrees):** never above U_img (max gap 0.0 bits, 24 cells), 1.5–3.0 bits below matched-dimension univariate.
   - **X5 (Herrmann–May split into block index and offset):** 2.0–3.5 bits below matched-dimension univariate in all 12 cells, and all 2 769 of its short vectors are functions of z = Hu + v.
   - **X6 (exact SVP vs LLL):** exact SVP gives exactly the LLL range in 30/30 cells (dimensions 8–16): LLL is not the bottleneck, the lattice is.
   - **X1:** 3960/3960 window sets equal the near-divisor sets.
   - **X2:** p was found in 60/60 blocks, with median 1.93 calls per ⌈log₂H⌉.
   - **Failure:** P2.2 FAILED as registered. The registered call bound was mis-stated: it missed a
     ceiling and assumed one divisor per block, but p and q can share a block. The corrected bound
     (Prop LOC) holds 60/60, post hoc.
5. **Partial-key-exposure frontier (§2).**
   - Polynomial time is known from half of p's bits (Coppersmith 1996/97; May 2010 Thm 1). Below
     that, the cost is linear in the missing range: N^{θ−1/4} for an interval of length N^θ.
   - No special setting known to us applies to one balanced N with contiguous top bits. That covers
     Herrmann–May (which needs *more* bits as chunks multiply), implicit factoring (needs several
     moduli), p^r q, and Heninger–Shacham (needs bits of both factors).
   - Lower bounds exist only for the method. The determinant certificate is optimal for the basic
     lattice, which equals the full truncated ideal lattice (Lemma MI). AASW give optimality under a
     Gaussian-heuristic assumption. Chinburg–Hemenway–Heninger–Scherr exclude auxiliary polynomials
     beyond M^{1/d} for a known modulus. List-size arguments do not apply to semiprimes.
   - **Nothing proves that factoring from fewer than half of p's bits is impossible in polynomial
     time.** It would separate partial key exposure from factoring, since p's full half is factoring
     itself.

**Barrier statement (where the limit bites).**
- Let e(θ) be the exponent of the best block solver for "find p | N in an interval of length N^θ near
  √N". The degree-j ladder then costs N^{1/2 − θ_j + e(θ_j)} with θ_j = 1/2 − 1/(j+2).
- The best known solvers give e(θ) = max(0, θ − 1/4), capped at 1/5 by Harvey–Hittmeir. That makes
  the ladder N^{1/4} for every j ≥ 2, whatever window structure is used.
- Reaching its nominal N^{1/(j+2)} needs e(θ_j) = o(1): polynomial-time factoring from a fraction
  2/(j+2) of p's top bits, i.e. 2/5 of them for j = 3 and 1/3 for j = 4.
- That lies strictly beyond the Coppersmith barrier. This is believed out of reach for lattice
  methods (heuristically, and provably for restricted classes, §2.3); no complexity-theoretic
  evidence either way is known.

## 1. The equivalence

### 1.1 Notation and identities (PROVED)

- Block start x, offset h, degree j, D = x^{j+1}.
- S_j(h) = Σ_{i≤j}(−h)^i x^{j−i}, num_j(h) = N·S_j(h), P_j(h) = num_j(h)/D.
- R(h) = N(−h)^{j+1}/(x^{j+1}(x+h)), so that N/(x+h) = P_j(h) + R(h).
- ρ = NH^{j+1}/x^{j+2} and K = ρx. At the top (x ≍ √N), K ≍ c^{j+1}B.

**(T1)** (x+h)S_j(h) = x^{j+1} − (−h)^{j+1}. This is the geometric series.

**(T2)** For every integer q, num_j(h) − qD = S_j(h)(N − q(x+h)) − q(−h)^{j+1}.

*Proof.* Expand S_j(h)N − qS_j(h)(x+h) − q(−h)^{j+1} with (T1). ∎

**(W) Window = near-divisor set.** Let 0 ≤ a ≤ b < 1 and d = x + h. Then {P_j(h)} ∈ [a, b] iff some
q ∈ Z gives k := N − qd with

  a·d + N(−h)^{j+1}/D ≤ k ≤ b·d + N(−h)^{j+1}/D.

*Proof.* P_j(h) − q = (N − qd)/d − R(h) = k/d − R(h), and R(h)·d = N(−h)^{j+1}/D. So P_j − q ∈ [a, b]
iff k lies in the stated range. {P_j} ∈ [a, b] iff P_j − q ∈ [a, b] for some q. ∎
- Divisors (k = 0) lie in the window [|R(h)|, |R(h)|] (even j) or [1 − R, 1 − R] (odd j).
- The ladder's window [0, ρ] (even j) catches exactly the d with |k| ≲ ρd ≤ K(1 + H/x).
- X1: 3960/3960 (block, window) pairs agree (20–28 bits; j = 2, 3, 4; divisor windows and random
  windows).

**(F) The curve relation.** Put w_j(h, r) = r − num_j(h) and r_j(h) := −N(−h)^{j+1}/(x+h). Then

  (x+h)·w_j(h, r_j(h)) = −N·D,   F_j := (x+h)r + N(−h)^{j+1} = (x+h)w_j + N·D.

For a divisor d = x + h₀ with cofactor q₀, the window residue r₀ := num_j(h₀) − q₀D = −q₀(−h₀)^{j+1}
(by T2) equals r_j(h₀). So **the divisor's window data lie on the curve F_j = 0, i.e. N = q(x+h)
multiplied out.** For a near-divisor the same computation gives F_j(h, r) = D·k (codex topblock_1 has
this identity). Checked symbolically for j ≤ 5 and in tests.

**Lemma TS (Taylor structure adds nothing, PROVED).** For integers h ≥ 0 and q, with x > h:
"num_j(h) − qD = −q(−h)^{j+1}" ⟺ N = q(x+h).

*Proof.* By (T2), the left side is S_j(h)(N − q(x+h)) = 0, and S_j(h) > 0 for 0 ≤ h < x by (T1). ∎

So the degree-j Taylor "window equation" for exact divisors is the factoring equation N = q(x+h),
times a nonzero factor, for every j at once. Different degrees, blocks or centres only change the
cofactor S_j. With q = q₀ + μ this is Coppersmith's bivariate integer equation (EUROCRYPT '96), whose
range is again H < N^{1/4}.

### 1.2 Reductions

**Prop KB (promised search = partial key exposure; PROVED; codex topblock_1).** "Find p given
p ∈ [x, x+H), H = Θ(N^{1/2 − 1/(j+2)})" and "find p given its top (1/(j+2))log₂N + O(1) bits" reduce to
each other with O(1) calls, by covering intervals with O(1) dyadic prefixes and vice versa.

**Prop LOC (window oracle ⇒ promised search; PROVED; corrected post hoc).**
- *Algorithm.* Bisect [1, H). On a subinterval I = [h₁, h₂), query the exact-divisor window
  W_I = hull{|R(h)| : h ∈ I} (even j; for odd j use 1 − R). Recurse only into intervals with a
  nonzero count. Test the block start separately.
- *Correctness.* Every divisor lies in its subinterval's window (W).
- *Cost.* Write ν₀ for the number of exact divisors in the block (≥ 1; it can be 2, see below).
  Write Dep = ⌈log₂(H−1)⌉ and K_ℓ = (j+1)ρ(x+H)(⌈(H−1)/2^ℓ⌉ − 1)/H. Then

    calls ≤ 1 + 2·Σ_{ℓ<Dep} (ν₀ + #{(h, k) : 0 < |k| ≤ K_ℓ, x + h | N − k}).

*Proof.*
1. R is increasing on [0, H) with R′ ≤ (j+1)N h^j/x^{j+2} ≤ (j+1)ρ/H. So |W_I| ≤ (j+1)ρ(|I|−1)/H.
2. By (W), an h ∈ I counted in W_I has |k| ≤ (x+H)|W_I| ≤ K_ℓ at depth ℓ.
3. A node at depth ℓ has a nonzero count only if it contains a divisor or such a near-divisor. Nodes
   at one depth are disjoint.
4. Each nonzero non-leaf node spawns two queries. ∎

*Average size.* Σ over the ladder's blocks of the near-divisor counts is ≪ B log N (Shiu, as in
Theorem L_j). So the bound is O(log N) per block on average, but it has no worst-case polylog form.

*The registered bound was wrong.* It used 2(log₂H + 1) and ν₀ = 1. It failed in 7/60 X2 blocks:
- 5 failures are the missing ceiling plus root term;
- 2 are blocks containing both p and q (|q − p| < H; e.g. 44 bits, seed 2, j = 4: q − p = 2526 <
  H = 12287).

The corrected bound holds in 60/60 (post hoc; `topblock_x2post.jsonl`).

**Removing the candidate factor (PROVED; codex topblock_1).** A degree-(j+1) query on the degree-j
block has window width ρ′ = NH^{j+2}/x^{j+3}, with ρ′(x+H) ≤ c^{j+2}(1 + c/B) < 1. So every counted h
has |k| < 1, i.e. k = 0. One count plus bisection then finds the divisors in O(log H) calls, with no
false candidates. At degree j alone, re-centring on pieces of length xN^{−1/(j+1)} does the same at
O(N^{1/((j+1)(j+2))} + log H) calls.

**Prop IV (divisor finder ⇒ window counts; PROVED).** By (W), #{h ∈ I : {P_j(h)} ∈ [a, b]} is
determined by the divisors in [x + h₁, x + h₂) of the integers N − k, for k in an explicit range of
length (b − a)(x + H) + O(ρ(x+H)) + 1. One call to a general "all divisors of M in an interval" solver
per k, plus exact checks, answers the query.
- X2(iv): 400/400 queries agree.
- *What this does not give.* Those M = N − k are not balanced semiprimes. A known-bits solver for N
  does not answer them, and with factoring allowed this is the circular bypass of floorsum §6.5. We
  know of no reduction from window counting to known-bits search for N alone (codex agrees:
  UNVERIFIED).

**Theorem EQ (PROVED, as the composition of the above).** Fix j ≥ 2 and balanced N. Let
θ_j = 1/2 − 1/(j+2).
- (a) A top-block window oracle with per-query cost T gives promised search on blocks of length
  N^{θ_j} in O(C·T) operations, with C as in Prop LOC. A degree-(j+1) oracle gives O(T log N).
- (b) A promised-search solver with cost T′ factors N in O(B log N)·T′. The ladder with the window
  oracle factors N in O(B log² N)·T.
- (c) Hence the ladder's exponent with any oracle implementation is at least that of
  "guess-and-decode" with the induced block solver, up to the candidate factor. Its target
  N^{1/(j+2)+o(1)} is reachable through *some* block primitive iff promised search at θ_j costs
  N^{o(1)}.
- (d) The window oracle is never easier than promised search. It also counts near-divisors, which
  are facts about the integers N − k (k ≠ 0), so it can only be harder.

**What the "extra information" is, quantitatively.**
- A window count over a block is #exact divisors + #near-divisors in the window's k-range. The
  second term is O(c^{j+2}) on average per block (heuristically ≈ ρH; ≪ log N on average by Shiu).
  It concerns N − k, not p.
- *Shifts* rotate the window (codex). They select other k-ranges.
- *Several degrees* give curves F_j = 0 that all reduce to N = q(x+h) (Lemma TS).
- *Consecutive blocks* are disjunctive: p lies in exactly one block (codex), so their conditions
  cannot be multiplied as if all held. Their aggregate over many blocks is
  Σ_{d∈[u,v)}(⌊(N+K)/d⌋ − ⌊(N−K−1)/d⌋), a truncated hyperbola count. Its best known exact algorithm is
  N^{1/3} (floorsum, floorhard §6), and interval access to it factors N by bisection (folklore).
- So the ladder's structure supplies no conjunctive information about p. Every lattice-usable
  relation it supplies lies in the ideal of the single relation N = q(x+h) (§3).

## 2. The partial-key-exposure frontier

### 2.1 Known positive results (references checked by codex topblock_1 and the literature critic, §7)

| Result | Statement | Relevance to one balanced N, contiguous top bits of p |
|---|---|---|
| Coppersmith, J. Cryptology 10 (1997) 233–260, Thm 4 (EUROCRYPT '96) | top (1/4)log₂N bits of p give p in poly time | the j = 2 ladder (half of p's bits) |
| Howgrave-Graham, CaLC 2001 (LNCS 2146); May 2010 survey Thm 1 | divisor b ≥ N^β, monic deg-δ f: all roots \|x₀\| ≤ c·N^{β²/δ} in O(cδ⁵log⁹N) | β = ½, δ = 1: N^{1/4}; **linear in c**, so extending the range by c costs c (subdivision) |
| fewer bits | brute-force the missing bits down to N^{1/4}: cost N^{θ−1/4} for an interval N^θ | the j ≥ 3 ladder; gives N^{1/4} in total |
| Herrmann–May, ASIACRYPT 2008, Thms 4–5 (Assumption 1) | k unknown chunks, ∏X_i ≤ N^{γ_k}, γ₁ = 1/4, γ₂ ≈ 0.2071, γ_k → (1−ln2)/2 ≈ 0.1534 | more chunks need **more** known bits (→ ln 2 ≈ 69% of p); adjacent chunks merge (§3.3) |
| May–Ritzenhofen, PKC 2009 (implicit factoring) | several moduli whose p_i share bits | needs several moduli; one N gives none |
| Boneh–Durfee–Howgrave-Graham, CRYPTO 1999 | N = p^r q: about 1/(r+1) of p's bits | r = 1 gives back ½ |
| Heninger–Shacham, CRYPTO 2009 | random bits of p **and** q (δ > 2 − √2 for p, q alone) | top bits of q are implied by those of p; not new information |
| Harvey–Hittmeir, Math. Comp. 91 (2022) | deterministic N^{1/5} with log savings, ignoring hints | caps e(θ) ≤ 1/5 |

### 2.2 The composition (PROVED arithmetic; the e(θ) is the best known to us)

- Let the block solver cost N^{e(θ)} on intervals of length N^θ, θ ∈ [1/4, 1/2]. Then the ladder costs
  N^{1/2−θ+e(θ)}.
- With e(θ) = θ − 1/4 (Coppersmith subdivision) this is N^{1/4} for every θ, i.e. every j ≥ 2.
- With e = 1/5 (Harvey–Hittmeir) it is worse except at θ = 1/2, where it is plain Harvey–Hittmeir.
- Beating N^{1/4} through the ladder needs e(θ) < θ − 1/4 for some θ. Beating the deterministic
  N^{1/5} needs e(θ) < θ − 3/10.
  - For j = 3 (θ = 3/10) that means polynomial-time recovery of p from 2/5 of its top bits.
  - For j = 4 (θ = 1/3) it means e < 1/30.

### 2.3 Lower-bound arguments for Coppersmith-type methods (precise scope)

1. **Determinant certificate (PROVED, classical).** For the lattice N^{m−i}f^i, h^k f^m (dimension w,
   f = x + h), det = N^{m(m+1)/2}X^{w(w−1)/2}. The certificate det^{1/w} < p^m/√w (times the LLL
   factor) needs ξ < 2βθ − θ² ≤ β², with θ = m/w, X = N^ξ, p ≥ N^β. So X < N^{1/4} for β = ½. This
   limits the certificate, not λ₁.
2. **Lemma MI (PROVED).** For monic linear f, the basic lattice *is* the full lattice
   (N, f)^m ∩ Z[h]_{<w}.
   *Proof.* Substitute z = x + h. The ideal becomes the monomial ideal (N, z)^m, whose elements are
   exactly the polynomials whose z^i-coefficient is divisible by N^{max(m−i, 0)}. The triangular rows
   N^{m−i}z^i and z^k·z^m span exactly that set in degree < w. ∎
   So within ideal-power lattices with a degree bound, no choice of shift polynomials improves the
   determinant. The only freedom left is λ₁ versus det^{1/w}.
3. **AASW** (Aono–Agrawal–Satoh–Watanabe, ACISP 2012; AAECC 29 (2018), Thm 2, per codex). For
   canonical constructions, failure at U ≥ 16N^{β²/d} is proved under their Assumption 1, a
   Gaussian-type lower bound λ₁ ≥ GH/2. This is conditional and specific to their lattice class.
   - *Our observation (X6, descriptive).* For the univariate ideal lattice beyond the threshold,
     λ₁ is far **below** GH (log₂(λ₁/GH) ∈ [−5.5, +0.6] over all 30 cells), because the lattice is highly structured (the
     constant row N^m stays short as X grows).
   - Yet exact SVP still fails the Howgrave-Graham test (P6.1). Whether this conflicts with AASW's
     Assumption 1 depends on their exact lattice and regime, which we could not retrieve (ePrint
     blocked scripted access). UNVERIFIED; see §7.
4. **Capacity theory** (Chinburg–Hemenway–Heninger–Scherr, ASIACRYPT 2016, Thm 2, per codex). For a
   *known* modulus M and monic degree-d f, no nonzero h = Σa_{uv}z^u(f/M)^v is uniformly < 1 on the
   disc of radius M^{1/d+ε}. So Coppersmith's known-modulus auxiliary-polynomial method is optimal
   for all degrees. The unknown-divisor (β < 1) case is not covered by that theorem as cited (§7).
5. **List size (PROVED, scope).** f = z^d mod s^d has Θ(M^ε) roots of size M^{1/d+ε}, so complete
   listing fails for some moduli. A balanced semiprime has at most two divisors in any short
   interval, so this gives **no** barrier for our problem.
6. **Johnson-bound analogy** (Cohn–Heninger, Adv. Math. Commun. 9 (2015), Thm 1.1). N^{β²} is the
   Johnson radius of the CRT code. Decoding beyond Johnson is not known to be impossible in general,
   and for a promised unique answer list size is no obstacle.

**Status.** "Polynomial-time factoring from (1/4 − ε)log₂N contiguous top bits of one balanced factor"
is **open**. No impossibility result exists for general algorithms, and none could exist without
proving something about factoring itself (θ = 1/2 is factoring). For lattice methods the barrier is:
- proved for the determinant certificate and for ideal-power lattices with degree bounds (Lemma MI);
- conditional (AASW) for canonical constructions with λ₁ assumptions;
- proved (CHHS) only in the known-modulus analogue.

## 3. Multi-window and multi-block lattices

### 3.1 The transfer theorem

**Setting.**
- Unknown h₀ = p − x ∈ [0, X), q₀ = N/p. Windows b = 1..k: offsets Δ_b, degrees j_b, centres
  x_b = x − Δ_b, D_b = x_b^{j_b+1}, h_b = h + Δ_b.
- w_b(h, r_b) = r_b − g_b(h), where g_b ≡ num_{j_b}(h_b) coefficientwise mod D_b. The root is
  r_{b0} = −q₀(−h_{b0})^{j_b+1}, and w_b(h₀, r_{b0}) ≡ 0 (mod D_b) by (T2), since x_b + h_{b0} = p.
- 𝓛 is any lattice spanned by finitely many polynomials of the form
  (monomial in h, r) × N^{α}(x+h)^{s} × ∏_b w_b^{t_b}D_b^{m2_b − t_b}, with α + s ≥ m1.
  This covers every shift strategy (basic, Jochemsz–May, Herrmann–May).
  - All of them vanish mod p^{m1}∏D_b^{m2_b} at the root.
  - Let T be the largest total r-degree among the spanning polynomials, and d the largest value of
    (h-degree + Σ_b j_b·r_b-degree) + T.
- Define σ: r_b ↦ ρ_b(h) := −N(−h_b)^{j_b+1}/(x+h). Note ρ_b(h₀) = r_{b0}, since x + h₀ = p.

**Prop T (PROVED; no prior formulation found by us).** For every g ∈ 𝓛:
- (i) G := (x+h)^T σ(g) ∈ ∏_b D_b^{m2_b}·(N, x+h)^{m1+T} ⊂ Z[h], with deg G ≤ d.
- (ii) G = 0 iff g vanishes identically on the curves r_b = ρ_b(h), i.e. g lies in the ideal of the
  F_b = (x+h)r_b + N(−h_b)^{j_b+1} over Q(h).
- (iii) Let Y_b ≥ max_{0≤h≤X}|ρ_b(h)| (e.g. Y_b = ⌈N/x⌉(X + |Δ_b|)^{j_b+1}). If
  ‖g(hX, r·Y)‖₁ < p^{m1}∏D_b^{m2_b} and G ≠ 0, then G′ := G/∏D_b^{m2_b} is a nonzero element of the
  univariate lattice L_{M,d} = (N, x+h)^M ∩ Z[h]_{≤d}, M = m1 + T, with

    ‖G′(hX)‖₁ < (3 + 3√2)^d·(1 + X/x)^T·p^M,

  and G′(h₀) = 0.

*Proof.*
1. Each spanning polynomial is a product of factors. Under (x+h)·σ:
   - (x+h)σ(r_b) = −N(−h_b)^{j_b+1} ∈ (N);
   - (x+h)σ(w_b) = (x+h)(ρ_b − num_b) + (x+h)(num_b − g_b) = −N·D_b + D_b(x+h)c_b(h) ∈ D_b·(N, x+h),
     by (F) applied at centre x_b. Here x_b + h_b = x + h, and num_b − g_b = D_b c_b with c_b ∈ Z[h];
   - D_b, N and (x+h) are fixed by σ.
2. Distribute (x+h)^T over the r-degree of each factor and pad with (x+h)^{T − deg_r}. Each spanning
   polynomial then maps into ∏D_b^{m2_b}·(N, x+h)^{m1 + T}: the N/(x+h)-part gives m1, each r or w
   gives one more power, and the padding supplies the rest. The degree count is direct.
3. σ is a ring homomorphism, so integer combinations stay in the ideal. This proves (i).
4. For (ii): σ(g) = 0 iff g ∈ ker σ, which over Q(h)[r] is generated by r_b − ρ_b, i.e. by F_b/(x+h).
5. For (iii):
   - On [0, X], |σ(g)(h)| ≤ ‖g(hX, rY)‖₁ because |ρ_b| ≤ Y_b. So sup|G| < (x+X)^T p^{m1}∏D_b^{m2_b},
     and (x+X) ≤ p(1 + X/x) since p ≥ x.
   - A real polynomial of degree d bounded by 1 on [0, X] has ‖·(hX)‖₁ ≤ 3^d(1+√2)^d: Chebyshev on
     [−1, 1], then s = 2u − 1.
   - G′(h₀) = 0 because g(h₀, r₀) = 0 over Z (Howgrave-Graham). ∎

**Consequences.**
- *(PROVED.)* Mixed success at (X, m1, T, d) requires a nonzero vector of L_{M,d}(X) of norm
  < c^d p^M with c = 3 + 3√2. Univariate success at the same parameters is certified by norm
  < p^M/√(d+1). In exponent terms (X = N^ξ, log N → ∞) the two thresholds differ by O(d log c /
  (M log N)) → 0.
- So **any** window-augmented lattice reaches at most the λ₁-range of the univariate ideal lattices
  L_{M,d}, plus o(1) in the exponent. Its determinant-certified range is N^{1/4} (Lemma MI, §2.3).
- *(HEURISTIC/CONDITIONAL.)* If, as believed (AASW-type assumptions), the univariate ideal lattices
  have no vector below p^M beyond N^{1/4+o(1)}, then no window-augmented lattice does either.
- The same substitution with q ↦ N/(x+h) maps Coppersmith's integer-equation lattices for
  N = (x+h)(q₀+μ) into the same univariate ideal: W·h^aμ^b ↦ W·h^a(N − q₀(x+h))^b/(x+h)^b. This
  explains why the bivariate integer formulation also stops at N^{1/4}. That is a remark; the known
  equivalence is Coron's.

### 3.2 The F-artifact (PROVED + X3)

- With m1, m2 ≥ 1 the polynomial F·N^{m1−1}D^{m2−1} = (row (t=1, e=1)) + (row (0, 0)) lies in the
  mixed lattice.
- Its scaled norm is ≍ N^{m1}X^{j+1}D^{m2−1}, against the threshold x^{m1}D^{m2}. So it is
  Howgrave-Graham-short iff X^{j+1} ≲ x^{j+1}(x/N)^{m1}, i.e. X ≲ N^{1/2 − m1/(2(j+1))} at x ≍ √N:
  - N^{1/3} for j = 2, m1 = 1;
  - N^{3/8} for j = 3, m1 = 1.
- With m1 = 0 (windows only) the boundary is xN^{−1/(j+1)}.
- **Naive determinant counts therefore suggest mixed ranges up to N^{1/3}+, but every such vector is
  in ker σ.**
- *A tempting heuristic, refuted.* In units of log N, a window row of r-degree t contributes
  t((j+1)ξ − j/2) beyond the D-threshold, which is negative for ξ < j/(2(j+1)). This suggests that
  window rows "help". The negative contribution is exactly the trivial relation F.
- X3 (m1 = 1 cells): the largest X with any Howgrave-Graham-good vector equals the largest X at which
  the F row is short, in all 18 cells (difference 0.0 in all 18). Zero nontrivial good vectors appear
  beyond U_img + 3 bits.

### 3.3 Splitting the unknown (Herrmann–May; PROVED + X5)

- Writing p = x₀ + Hu + v (block index u, offset v) is an affine change of variables of the single
  contiguous unknown z = Hu + v.
- The HM lattice u^a f^b N^{(m−b)+}, a + b ≤ E, is the full truncated ideal (N, v′)^m in total degree
  ≤ E, with v′ = f. It contains the univariate lattice in z (Lemma MI in the variables (u, v′)).
- Its extra u-shifts cannot exploit the box U × H beyond the segment z < UH: the two have the same
  lattice points.
- The generic HM bound γ₂ ≈ 0.2071 < 1/4 is for non-adjacent chunks, whose coefficients are
  independent. For adjacent chunks the univariate sublattice dominates.
- X5: X5: in all 12 cells the HM range is 2.0–3.5 bits *below* the best univariate lattice of no larger dimension (15 or 28). All 2 769 Howgrave-Graham-good HM vectors have the linear factor b(H_s u + v) + c, i.e. they are functions of z, and recovery used that factor. Splitting therefore only wastes dimension.

## 4. Experiments

All pre-registered (`d661588`). Toy sizes are correctness and range checks, not costs.

| Exp | What | Outcome |
|---|---|---|
| X1 | (W) on 3960 (block, window) pairs, 20–28 bits, j = 2–4 | **P1 ✓** 0 mismatches |
| X2 | Prop LOC on 60 top blocks (40–48 bits); Prop IV on 400 queries (32 bits) | **P2.1 ✓** 60/60; **P2.2 ✗** 53/60 (bound mis-stated; corrected bound 60/60 post hoc); **P2.3** ✓ median 1.93 (max 3.64); **P2.4 ✓** 400/400 |
| X3 | single-window MIX, 54 cells (64/96/128 bits, j = 2, 3, three configs) | **P3.1 ✓** (max gap 0.0); **P3.2 ✓** (1.0–4.5 bits below); **P3.3 ✓** (18/18; 0 nontrivial); **P3.4 ✓** (0.0) |
| X4 | two windows (adjacent blocks; degrees 2 + 3), 24 cells | **P4.1 ✓** (max gap 0.0); **P4.2 ✓** (1.5–3.0 below) |
| X5 | HM split, 12 cells | **P5.1 ✓** (2.0–3.5 below); **P5.2 ✓** (100%) |
| X6 | exact SVP vs LLL, univariate, 30 cells | **P6.1 ✓** (gap 0.0 in 30/30) |

**X3 table** (medians over seeds 0–2; log₂ values).

| bits | j | (m1, m2, E) | dim | X_max MIX | X_max U_img (M, n) | best U with dim ≤ MIX | largest log₂X with an HG-good vector | F row short up to |
|---|---|---|---|---|---|---|---|---|
| 64 | 2 | (1,1,8) | 16 | 12.0 | 12.0 (2,10) | 14.5 | 20.0 | 20.0 |
| 64 | 2 | (2,1,10) | 20 | 13.0 | 13.0 (3,12) | 14.5 | 14.0 | 10.0 |
| 64 | 2 | (2,2,10) | 27 | 13.5 | 13.5 (4,13) | 14.5 | 14.5 | 10.0 |
| 64 | 3 | (1,1,9) | 17 | 12.0 | 12.0 (2,11) | 14.5 | 23.0 | 23.0 |
| 64 | 3 | (2,1,11) | 21 | 13.0 | 13.0 (3,13) | 14.5 | 14.0 | 14.0 |
| 64 | 3 | (2,2,12) | 30 | 13.5 | 13.5 (4,15) | 15.0 | 14.5 | 14.5 |
| 96 | 2 | (1,1,8) | 16 | 18.5 | 18.5 (2,10) | 21.5 | 30.5 | 30.5 |
| 96 | 2 | (2,1,10) | 20 | 19.5 | 19.5 (3,12) | 22.0 | 20.5 | 15.0 |
| 96 | 2 | (2,2,10) | 27 | 20.5 | 20.5 (4,13) | 22.5 | 21.5 | 15.0 |
| 96 | 3 | (1,1,9) | 17 | 18.5 | 18.5 (2,11) | 21.5 | 34.5 | 34.5 |
| 96 | 3 | (2,1,11) | 21 | 19.5 | 19.5 (3,13) | 22.0 | 20.5 | 20.5 |
| 96 | 3 | (2,2,12) | 30 | 20.5 | 20.5 (4,15) | 22.5 | 21.5 | 21.5 |
| 128 | 2 | (1,1,8) | 16 | 24.5 | 25.0 (2,10) | 29.0 | 41.5 | 41.5 |
| 128 | 2 | (2,1,10) | 20 | 26.5 | 26.5 (3,12) | 29.5 | 27.5 | 20.5 |
| 128 | 2 | (2,2,10) | 27 | 27.5 | 27.5 (4,13) | 30.0 | 28.5 | 20.5 |
| 128 | 3 | (1,1,9) | 17 | 24.5 | 25.0 (2,11) | 29.0 | 46.5 | 46.5 |
| 128 | 3 | (2,1,11) | 21 | 26.5 | 26.5 (3,13) | 29.5 | 27.5 | 27.5 |
| 128 | 3 | (2,2,12) | 30 | 27.5 | 27.5 (4,15) | 30.0 | 28.5 | 28.5 |

Reading the table:
- The mixed range depends only on the transfer target (M, E + T), not on j. For example, j = 2 with
  (2,1,10) and j = 3 with (2,1,11) give the same range.
- In the m1 = 1 rows, the "good" extent (≈ N^{1/2 − 1/(2(j+1))}, e.g. 41.5 of 64 half-bits at
  128 bits, j = 2) is pure F-artifact.
- In m1 = 2 rows the F row is longer, and good vectors stop about 1 bit above X_max.

## 5. Pre-registration outcomes

| Prediction | Result |
|---|---|
| P1: 0 mismatches | ✓ 0/3960 |
| P2.1: p found 100% | ✓ 60/60 |
| P2.2: calls ≤ registered bound, 100% | **✗ 53/60** (registered bound mis-stated; corrected Prop LOC bound 60/60, post hoc) |
| P2.3: median calls/⌈log₂H⌉ ≤ 3.0 | ✓ 1.93 (max 3.64) |
| P2.4: (iv) 100% agreement | ✓ 400/400 |
| P3.1: MIX ≤ U_img + 3, every cell | ✓ max gap 0.0 (equality in 51/54; −0.5 in 3) |
| P3.2: MIX ≤ U_dim + 3, every cell | ✓ MIX is 1.0–4.5 bits **below** |
| P3.3: good vectors beyond U_img + 4 in ≥ 80% of m1 = 1 cells; 0 nontrivial beyond U_img + 3 | ✓ 18/18; 0 nontrivial |
| P3.4: good extent − F prediction ∈ [−3, 3] | ✓ 0.0 in all 18 |
| P4.1: MIX2 ≤ U_img + 3 | ✓ max gap 0.0 (equality 20/24) |
| P4.2: MIX2 ≤ U_dim + 3 | ✓ 1.5–3.0 bits below |
| P5.1: HM ≤ U_dim + 1 | ✓ 2.0–3.5 bits below |
| P5.2: ≥ 80% of good HM vectors line-type | ✓ 2 769/2 769 (100%) |
| P6.1: 0 ≤ SVP − LLL ≤ 1.0 | ✓ 0.0 in 30/30 |

Decision rules:
- H-MIX is **rejected as a mechanism**: windows add no range.
- H-HM is supported.
- H-SVP is supported.

**Disclosures.**
- P2.2's bound was wrong as registered. The corrected statement (Prop LOC) was derived after
  seeing the 7 failures, and its 60/60 check is post hoc.
- The development peeks listed in the pre-registration were made before it was committed.

## 6. Answers to the four questions

1. **Equivalence.**
   - Promised divisor search in a top block is *exactly* partial key exposure with
     (1/(j+2))log₂N top bits (Prop KB).
   - The window oracle implies it (Prop LOC; candidate-free at degree j+1) and is otherwise strictly
     stronger (Prop IV).
   - The window structure carries no extra information about p:
     - the Taylor structure is the identity N = q(x+h) (Lemma TS);
     - consecutive blocks are disjunctive;
     - multiple windows and degrees lie on trivial curves;
     - the only extra content is near-divisor data about N − k.
   - So the ladder reduces exactly to partial key exposure, and its limit is the Coppersmith barrier.
2. **Frontier.** §2. Polynomial time from ½ of p's bits; below that the cost is linear in the missing
   range. No structured setting applies to one N with contiguous top bits. Lower bounds are
   method-specific: the determinant certificate, Lemma MI, conditional AASW, and CHHS for known
   moduli. No general impossibility.
3. **Multi-window lattices.**
   - They cannot beat the univariate ideal lattice (Prop T), and naive counts are inflated by the
     trivial F-vectors (§3.2).
   - Splitting into chunks only reparametrises (§3.3).
   - Measured: X3–X5.
4. **Conclusion.** No mechanism. The ladder at degree j is exactly as hard as factoring from a fraction
   2/(j+2) of p's top bits. With known solvers it costs N^{1/4} for every j ≥ 2, and the barrier sits
   at Coppersmith's half-the-bits.

## 7. Reviews

__REVIEWS__

## 8. Reproduction

```bash
.venv/bin/python -m pytest -q tests/test_topblock.py
for e in x1 x2 x3 x4 x5 x6 summary; do PYTHONPATH=. .venv/bin/python -m src.topblock.run_experiments $e; done
# x3 ~25 min, x4/x5/x6 a few minutes each; fpylll needs cysignals in .venv
```

