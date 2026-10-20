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

**We found nothing in the structured windows beyond plain known-top-bits factoring, and we can
prove a precise version of that for lattice methods.**
- The ladder's per-block task at degree j is exactly "find p from its top (1/(j+2))·log₂N bits"
  (partial key exposure). Its window oracle implies that task, and nothing shows it to be easier.
- Every Coppersmith-type lattice built from the window congruences (any blocks, degrees or centres,
  within the product ideal, with the Howgrave-Graham test) is dominated by the plain univariate
  Coppersmith lattice (exact SVP). The loss is a poly(d) factor in the range (Prop T): o(1) in the
  exponent whenever the degree d is N^{o(1)}.
  This is a comparison theorem, not a proof that N^{1/4} is a universal barrier.
- With the hint-using block solvers known to us (Coppersmith plus subdivision), the ladder costs
  N^{1/4+o(1)} for every j ≥ 2. Its nominal N^{1/(j+2)} needs factoring from a 2/(j+2) fraction of
  p's top bits in N^{o(1)} time. That is beyond the Coppersmith barrier for j ≥ 3, and open.

What is established:

1. **The equivalence (§1, PROVED).**
   - *Windows are near-divisor sets.* A top-block window holds exactly the d = x + h for which N − k
     is a multiple of d, with k in an exactly computable range (identity (W)). For the ladder's
     width-ρ windows that range has |k| = O(B).
   - *The Taylor structure is the trivial relation.* For a divisor d with cofactor q, the window
     residue is r = −q(−h)^{j+1}. It lies on the curve F_j = (x+h)r + N(−h)^{j+1} = 0, which is the
     identity N = q(x+h) multiplied out (T2). Every degree, block and centre gives such a curve. (This
     is an algebraic statement. It does not by itself show that the integer points on F_j carry no
     usable structure.)
   - *Oracle ⇒ known-bits search:* bisection with shrinking windows uses a proved number of calls
     (Prop LOC). Codex's degree-(j+1) query removes all false candidates, so a polylog
     degree-(j+1) oracle gives O(log H) calls.
   - *The converse is not known:* answering window counts takes O((b − a)x + |I| + ρx) calls of a
     general divisors-in-interval finder on the integers N − k (Prop IV). A known-bits solver for N
     alone does not obviously suffice. So the oracle implies the ladder's real need and is not known
     to be implied by it.
   - *Hence (Theorem EQ):*
     - a window oracle with worst-case per-query cost T gives promised search in O(C·T), or
       O(T log N) at degree j+1;
     - promised search in N^{o(1)} gives the ladder's exponent through its outer loop.
     - Ladder implementations that amortise work across blocks are not covered.
     - The extra content of a window count is the near-divisor multiset, which concerns the other
       integers N − k. It says nothing more about p.
2. **Multi-window / multi-block lattices (§3).**
   - **Prop T (transfer theorem, PROVED).** Take any lattice inside the product ideal of the divisor
     congruence x + h ≡ 0 (mod p) and any number of window congruences (from any blocks, Taylor
     degrees and centres), with the Howgrave-Graham certificate. The substitution r_b ↦ −N(−h_b)^{j_b+1}/(x+h) maps it into
     ∏D_b^{m2_b}·(N, x+h)^{m1+T} ⊂ Z[h], the ordinary univariate Coppersmith ideal. Its kernel
     consists of Q(h)-multiples of the curve relations F_b, which vanish on the whole curve and locate
     nothing, even jointly.
   - Keeping the vanishing at h₀, a successful mixed lattice at box X yields a Howgrave-Graham-
     certified vector of the univariate lattice L_{M,d} at box X/O(d^{3.5}) (Prop T(iv)). So
     mixed range ≤ O(d^{3.5}) × (exact-SVP univariate range).
   - Windows add no range beyond what the univariate ideal itself supports, up to that factor.
3. **Why naive counts mislead (PROVED + measured).** Multiples of F lie in the mixed lattice. The
   vector (x+h)^{m1−1}F·D^{m2−1} is Howgrave-Graham-short up to X ≈ N^{j/(2(j+1))} for every m1:
   N^{1/3} for j = 2, N^{3/8} for j = 3. So mixed lattices show "good" short vectors far beyond
   N^{1/4}, but all of them lie in ker σ.
   - E.g. at 128 bits the good vectors reach log₂X = 41.5, while the actual range is 24.5.
4. **Experiments (pre-registered).**
   - **X3:** in all 54 cells, the mixed range equals the transfer target U_img to within 0.5 bit. It
     is 1.0–4.5 bits *below* the best tested univariate lattice of no larger dimension. The trivial-vector extent
     matches the F prediction exactly (difference 0.0 bits in all 18 m1 = 1 cells).
   - **X4 (two windows: adjacent blocks, or two degrees):** never above U_img (max gap 0.0 bits, 24 cells), 1.5–3.0 bits below matched-dimension univariate.
   - **X5 (Herrmann–May split into block index and offset):** 2.0–3.5 bits below matched-dimension
     univariate in all 12 cells. All 2 769 of its short vectors have a linear factor b(Hu + v) + c (not
     all are functions of z).
   - **X6 (exact SVP vs LLL):** exact SVP gives exactly the LLL range in 30/30 cells (dimensions 8–16,
     with the HG certificate). For these parameters LLL is not the bottleneck; the lattice is.
   - All experimental statements are scoped to the tested parameters and dimensions.
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
   - Lower bounds exist only for the method.
     - The determinant certificate cannot pass N^{1/4}, and by Lemma MI no shift choice improves it
       within ideal-power lattices with a degree bound. That is not a proof that those lattices fail.
     - AASW prove failure of the certificate under a Gaussian-heuristic assumption. On this family the
       posited inequality is false for every instance once U ≥ N^{1/4}, so it cannot be applied here
       (§2.3).
     - Chinburg–Hemenway–Heninger–Scherr exclude auxiliary polynomials beyond M^{1/d}, but only for a
       known modulus.
     - List-size arguments do not apply to semiprimes.
   - **Nothing proves that factoring from fewer than half of p's bits is impossible in polynomial
     time.** It would separate partial key exposure from factoring, since p's full half is factoring
     itself.

**Barrier statement (where the limit bites).**
- Let e(θ) be the exponent of the best block solver for "find p | N in an interval of length N^θ near
  √N". The degree-j ladder then costs N^{1/2 − θ_j + e(θ_j)} with θ_j = 1/2 − 1/(j+2).
- The only hint-using block solver known to us is Coppersmith plus subdivision, with
  e(θ) = max(0, θ − 1/4). It makes the ladder N^{1/4+o(1)} for every j ≥ 2.
- Hint-ignoring solvers (Harvey–Hittmeir, e = 1/5) would formally give N^{1/5 + 1/(j+2)}. That is
  always worse than running Harvey–Hittmeir once on N, so they are irrelevant to the ladder.
- Reaching the nominal N^{1/(j+2)} needs e(θ_j) = o(1): factoring in N^{o(1)} from a fraction 2/(j+2)
  of p's top bits, i.e. 2/5 of them for j = 3 and 1/3 for j = 4. Polynomial time would be stronger.
- That lies beyond the Coppersmith barrier. For lattice methods the barrier is supported only by:
  - determinant arguments;
  - the conditional AASW result;
  - the known-modulus analogue (CHHS);
  - our comparison theorem Prop T.

  No complexity-theoretic evidence either way is known.

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
- Divisors (k = 0) with h ≥ 1 lie in the window [|R(h)|, |R(h)|] (even j) or [1 − R, 1 − R] (odd j).
  At h = 0, {P(0)} = 0.
- *Odd j.* The integer k of (W) is ≈ d, because q = ⌊P⌋ = q₀ − 1. "Near-divisor with small |k|" then
  means the representative k − d, which gives the same divisibility.
- The ladder's window [0, ρ] (even j) catches a subset of the d with |k| ≤ ρd ≤ K(1 + H/x), exactly
  described by (W): k ∈ [−|R|d, (ρ − |R|)d].
- X1: 3960/3960 (block, window) pairs agree (20–28 bits; j = 2, 3, 4; divisor windows and random
  windows).

**(F) The curve relation.** Put w_j(h, r) = r − num_j(h) and r_j(h) := −N(−h)^{j+1}/(x+h). Then

  (x+h)·w_j(h, r_j(h)) = −N·D,   F_j := (x+h)r + N(−h)^{j+1} = (x+h)w_j + N·D.

With the coefficientwise-reduced g ≡ num_j (mod D) of the lattices, the same identity reads
(x+h)(r − g) + ND = F + D(x+h)c(h), with c ∈ Z[h] (codex topblock_2).

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
*method guarantee* is again H < N^{1/4}.

*Scope (codex).* This rewrites the problem; it does not prove that the rewritten structure is
computationally useless. In particular, integer points on F_j = 0 satisfy x + h | N h^{j+1}. When
gcd(x+h, x) = 1 that forces x + h | N, so bounded-integer-point algorithms on that curve are exactly
divisor search again.

### 1.2 Reductions

**Prop KB (promised search = partial key exposure; PROVED; codex topblock_1).** "Find p given
p ∈ [x, x+H), H = Θ(N^{1/2 − 1/(j+2)})" and "find p given its top (1/(j+2))log₂N + O(1) bits" reduce to
each other with O(1) calls, by covering intervals with O(1) dyadic prefixes and vice versa.

*Prior art for the rungs (literature critic).*
- In residue-class (low-bits) form, Lenstra 1984 is the j = 1 rung: S³ > N gives N^{1/3+o(1)}.
- Coppersmith–Howgrave-Graham–Nagaraj 2008 is the j = 2 rung (N^{1/4}).
- Rivest–Shamir (EUROCRYPT '85) need (1/3)log₂N top bits, i.e. 2/3 of p's bits.
- Khomovsky (arXiv 1802.03658) already uses Taylor approximations of n/x with floored coefficients
  (degrees 1, 2, 4), but makes no link to known bits.
- We found no explicit general statement "degree-j block ⇔ (1/(j+2))log₂N top bits".

**Prop LOC (window oracle ⇒ promised search; PROVED; corrected post hoc).**
- *Algorithm.* Bisect [1, H). On a subinterval I = [h₁, h₂), query the exact-divisor window
  W_I = hull{|R(h)| : h ∈ I} (even j; for odd j use 1 − R). Recurse only into intervals with a
  nonzero count. Test the block start separately.
- *Correctness.* Every divisor lies in its subinterval's window (W).
- *Cost.* Write ν₀ = max(1, #{exact divisors x + h of N : 1 ≤ h < H}). It can be 2, see below.
  Write Dep = ⌈log₂(H−1)⌉ and K_ℓ = (j+1)ρ(x+H)(⌈(H−1)/2^ℓ⌉ − 1)/H. Then

    calls ≤ 1 + 2·Σ_{ℓ<Dep} (ν₀ + #{(h, k) : 0 < |k| ≤ K_ℓ, x + h | N − k}).

*Proof.*
1. |R| is increasing on [0, H) (R itself decreases for even j), with
   |R|′ ≤ (j+1)N h^j/x^{j+2} ≤ (j+1)ρ/H. So |W_I| ≤ (j+1)ρ(|I|−1)/H.
2. By (W), an h ∈ I counted in W_I has |k| ≤ (x+H)|W_I| ≤ K_ℓ at depth ℓ. For odd j use the
   representative k − d.
   - The split mid = ⌊(h₁+h₂)/2⌋ gives maximal length ⌈(H−1)/2^ℓ⌉ at depth ℓ.
   - Non-leaf nodes have length ≥ 2, hence depth < Dep.
3. A node at depth ℓ has a nonzero count only if it contains a divisor or such a near-divisor. Nodes
   at one depth are disjoint.
4. Each nonzero non-leaf node spawns two queries. ∎

*Average size.* Σ over all the ladder's blocks of the near-divisor counts is ≪ B log N (Shiu, as in
Theorem L_j). So the bound is O(log N) per block averaged over all O(B log N) blocks, and O(log² N)
for the top dyadic range alone. It has no worst-case polylog form. The proof critic re-checked it on
394 random cases, odd and even j, 0 violations.

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
length ≤ (b − a)(x + H) + min(a, 1 − b)|I| + O(ρ(x+H)) + 1, choosing the representative q or q + 1.
The draft omitted the window-position term min(a, 1 − b)|I| (both critics); the code's range was
correct. One call to a general "all divisors of M in an interval" solver
per k, plus exact checks, answers the query.
- X2(iv): 400/400 queries agree.
- *What this does not give.* Those M = N − k are not balanced semiprimes. A known-bits solver for N
  does not answer them, and with factoring allowed this is the circular bypass of floorsum §6.5. We
  know of no reduction from window counting to known-bits search for N alone (codex agrees:
  UNVERIFIED).

**Theorem EQ (PROVED, as the composition of the above; restated after review).** Fix j ≥ 2 and balanced N. Let
θ_j = 1/2 − 1/(j+2).
- (a) A top-block window oracle with per-query cost T gives promised search on blocks of length
  N^{θ_j} in O(C·T) operations, with C as in Prop LOC. A degree-(j+1) oracle gives O(T log N).
- (b) A promised-search solver with cost T′ factors N in O(B log N)·T′. The ladder with the window
  oracle factors N in O(B log² N)·T.
- (c) Suppose the window oracle has worst-case per-query cost T. Then promised search at θ_j costs
  O(C·T), or O(T log N) with a degree-(j+1) oracle. Hence a ladder whose per-block work is uniformly
  N^{o(1)} yields promised search in N^{o(1)}·C.
  - Conversely, promised search in N^{o(1)} gives an N^{1/(j+2)+o(1)} factoring algorithm through the
    ladder's outer loop.
  - Implementations that amortise or share work across blocks are **not covered**.
  - N^{o(1)} is not polynomial time.
- (d) The window oracle implies promised search, up to C (or exactly, at degree j+1). No reduction in
  the other direction is known. Its extra output, near-divisor counts for the integers N − k (k ≠ 0),
  is not known to be easier or harder.
- (d′) Context: general oracles need far fewer bits. Maurer ("Factoring with an oracle", EUROCRYPT '92,
  LNCS 658) factors heuristically in polynomial time with ε·log N answers to arbitrary oracle
  questions, using smooth numbers in short intervals. The window oracle is a restricted oracle.

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
| Herrmann–May, ASIACRYPT 2008, LNCS 5350, 406–424, Thms 3–5 (Assumption 1); time poly(log N, (e/ε)^k) | k unknown chunks, ∏X_i ≤ N^{γ_k}, γ₁ = 1/4, γ₂ ≈ 0.2071, γ_k → (1−ln2)/2 ≈ 0.1534 | more chunks need **more** known bits (→ ln 2 ≈ 69% of p); adjacent chunks merge (§3.3) |
| May–Ritzenhofen, PKC 2009 (implicit factoring) | several moduli whose p_i share bits | needs several moduli; one N gives none |
| Boneh–Durfee–Howgrave-Graham, CRYPTO 1999 | N = p^r q: about 1/(r+1) of p's bits | r = 1 gives back ½ |
| Heninger–Shacham, CRYPTO 2009 | random bits of p **and** q (δ > 2 − √2 for p, q alone) | top bits of q are implied by those of p; not new information |
| Harvey, Math. Comp. 90 (2021); Harvey–Hittmeir, Math. Comp. 91 (2022) | deterministic N^{1/5} (Harvey), improved by a log-log factor to O(N^{1/5}log^{16/5}N/(log log N)^{3/5}); ignores hints | caps e(θ) ≤ 1/5 |

### 2.2 The composition (PROVED arithmetic; the e(θ) is the best known to us)

- Let the block solver cost N^{e(θ)} on intervals of length N^θ, θ ∈ [1/4, 1/2]. Then the ladder costs
  N^{1/2−θ+e(θ)}.
- With e(θ) = θ − 1/4 (Coppersmith subdivision) this is N^{1/4} for every θ, i.e. every j ≥ 2.
- With e = 1/5 (Harvey–Hittmeir, which ignores the hint) the formula gives 1/5 + 1/(j+2). That is
  below 1/4 for j ≥ 19, but always above running Harvey–Hittmeir once on N. So hint-ignoring solvers
  never make the ladder worthwhile.
- Beating N^{1/4} through the ladder needs e(θ) < θ − 1/4 for some θ. Beating the deterministic
  N^{1/5} needs e(θ) < θ − 3/10.
  - For j = 3 (θ = 3/10) even e = 0 only *matches* N^{1/5}.
  - For j = 4 (θ = 1/3) it needs e < 1/30.
  - Large j needs e(θ_j) < 1/5 − 1/(j+2).

### 2.3 Lower-bound arguments for Coppersmith-type methods (precise scope)

1. **Determinant certificate (PROVED, classical).** For the lattice N^{m−i}f^i, h^k f^m (dimension w,
   f = x + h), det = N^{m(m+1)/2}X^{w(w−1)/2}. The certificate det^{1/w} < p^m/√w (times the LLL
   factor) needs ξ < 2βθ − θ² ≤ β², with θ = m/w, X = N^ξ, p ≥ N^β. So X < N^{1/4} for β = ½. This
   limits the certificate, not λ₁.
2. **Lemma MI (KNOWN: AASW 2018 §4.1–4.2, Def. 1, for general degree; proof here for completeness).**
   For monic linear f, the basic lattice *is* the full lattice
   (N, f)^m ∩ Z[h]_{<w}.
   *Proof.* Substitute z = x + h. The ideal becomes the monomial ideal (N, z)^m, whose elements are
   exactly the polynomials whose z^i-coefficient is divisible by N^{max(m−i, 0)}. The triangular rows
   N^{m−i}z^i and z^k·z^m span exactly that set in degree < w. ∎
   So within ideal-power lattices with a degree bound, no choice of shift polynomials improves the
   determinant. The only freedom left is λ₁ versus det^{1/w}.
3. **AASW** (Aono–Agrawal–Satoh–Watanabe, ACISP 2012; AAECC 29 (2018), Thm 2; §4.2 read by codex
   topblock_2).
   - *Their lattice.* For degree δ, L_d is the span of h^{j_i}N^{m−i_i}F^{i_i} with
     i_i = min(m, ⌊i/δ⌋) and j_i = i − δi_i, for 0 ≤ i ≤ d, scaled by U. For δ = 1 this is exactly
     our univariate lattice.
   - *Their results.*
     - Thm 2: under Assumption 1 (λ₁(L_d) ≥ ½GH(L_d) with high probability over (F, N, U) and d),
       the technique fails for U ≥ 4²N^{β²/δ}, in an overwhelming proportion of (F, N, U).
     - Thm 1 (unconditional) gives only U ≥ 4N^{β/δ}, i.e. N^{1/2} at β = ½, which says nothing here.
     - Their "failure" is certificate failure: every canonical h has sup_{[−U,U]}|h| > W. This matches
       item 1.
   - **On this factoring family the posited inequality is false (PROVED; codex, literature critic).**
     - For m = d = 1, F = h + x with x ≍ √N and 1 ≪ U ≪ √N: λ₁ ≤ √(x² + U²) ≍ √N while GH ≍ √(NU),
       so λ₁/GH → 0.
     - In general ‖F^m‖ ≈ N^{m/2}, while GH ≈ N^{ξd/2 + m(m+1)/(2(d+1))}. So F^m < GH whenever
       ξ > m(d−m)/(d(d+1)), a threshold always < 1/4. Hence the inequality fails on every instance once
       U ≥ N^{1/4}, exactly AASW's claimed-failure regime.
     - Assumption 1 is stated for random (F, N, U, d), and random F = x + a has a ≍ N, while here
       a ≍ √N. So the assumption is not refuted, but it cannot be applied.
     - F^m also blocks AASW's proof route (their Claim 1 / condition (6)). F^m still satisfies their
       failure definition (sup|F^m| = (x+U)^m > p^m), so this is **not** evidence that their
       conclusion is false.
   - *Measured (X6, descriptive).* log₂(λ₁/GH) ∈ [−5.5, +0.6] beyond the fixed-parameter LLL endpoint
     (every X6 point lies below nominal N^{1/4}). The draft attributed the short vector to the
     constant row N^m; codex showed it is F^m.
   - So AASW's conditional optimality cannot simply be imported for partial key exposure. Exact SVP
     still fails the HG test at those X (P6.1).
   - May (2010, right after Thm 1) already remarks that better lattice-reduction approximations do not
     improve the bound. X6 is consistent with that.
4. **Capacity theory** (Chinburg–Hemenway–Heninger–Scherr, ASIACRYPT 2016, Thm 2, per codex). For a
   *known* modulus M and monic degree-d f, no nonzero h = Σa_{uv}z^u(f/M)^v is uniformly < 1 on the
   disc of radius M^{1/d+ε}. So Coppersmith's known-modulus auxiliary-polynomial method is optimal
   for all degrees.
   - The unknown-divisor (β < 1) case is not covered. CHHS §2.3.1 call it future research, and the
     literature critic found no follow-up (that negative is UNVERIFIED).
   - CHHS Thm 3 shows that integer-valued (binomial) auxiliary polynomials beating 1/d do exist. They
     are excluded only via the small-prime-factor argument (Thm 4). This is relevant to the "not
     covered" list of Prop T.
5. **List size (scope).**
   - *Known modulus (PROVED, elementary).* f = z^d mod s^d has Θ(M^ε) roots of size M^{1/d+ε}, so
     complete listing fails for some moduli.
   - *The relevant unknown-divisor analogue is the known-low-bits twin of Prop KB.*
     - Lenstra, Math. Comp. 42 (1984) 331–340: polynomial time for divisors ≡ r mod s with s > N^{1/3},
       and a count O((α − 1/4)^{−2}) for every α > 1/4 (non-constructive).
     - Coppersmith–Howgrave-Graham–Nagaraj, Math. Comp. 77 (2008) 531–545: constructive for α > 1/4,
       with count O((α − 1/4)^{−3/2}).
     - Hales (arXiv 2410.05030, §4.2) surveys lower-bound constructions: only 6–7 divisors near
       α ≈ 1/3, and no large list size below 1/4 is known.
   - A balanced semiprime has at most two divisors in any short interval, so list size gives **no**
     barrier for our problem.
   - NP-completeness of "is there a root in an interval" (Manders–Adleman 1978, quadratic congruences)
     does not apply to linear f modulo a factor of a semiprime.
6. **Johnson-bound analogy.** N^{β²} is the analogue of the Johnson radius.
   - Cohn–Heninger (Adv. Math. Commun. 9 (2015) 311–339; ICS 2011), Thm 1.2 and §4.2, show that the
     polynomial analogue is equivalent to Guruswami–Sudan list decoding. Their Thm 1.1 is the integer
     Coppersmith / Howgrave-Graham / May statement.
   - For CRT codes see Guruswami–Sahai–Sudan (FOCS 2000) and Boneh (STOC 2000 / JCSS 2002); page
     numbers UNVERIFIED.
   - Ben-Sasson–Kopparty–Radhakrishnan (IEEE Trans. IT 56 (2010)) give superpolynomial Reed–Solomon
     list sizes just beyond Johnson. That is an analogue only (UNVERIFIED as cited).
   - Decoding beyond Johnson is not known to be impossible in general, and for a promised unique
     answer list size is no obstacle.

**Status.** "Polynomial-time factoring from (1/4 − ε)log₂N contiguous top bits of one balanced factor"
is **open**. No impossibility result exists for general algorithms, and none could exist without
proving something about factoring itself (θ = 1/2 is factoring). For lattice methods the barrier is:
- proved only for the determinant certificate, which by Lemma MI is optimal among ideal-power
  lattices with a degree bound;
- actual failure is proved only conditionally (AASW, whose assumption cannot be applied to this
  family) or in the known-modulus analogue (CHHS).

Herrmann–May's γ_k are distribution-independent *sufficient* bounds, not lower bounds.

## 3. Multi-window and multi-block lattices

### 3.1 The transfer theorem

**Setting.**
- Unknown h₀ = p − x ∈ [0, X), q₀ = N/p. Windows b = 1..k: offsets Δ_b, degrees j_b, centres
  x_b = x − Δ_b, D_b = x_b^{j_b+1}, h_b = h + Δ_b.
- w_b(h, r_b) = r_b − g_b(h), where g_b ≡ num_{j_b}(h_b) coefficientwise mod D_b. The root is
  r_{b0} = −q₀(−h_{b0})^{j_b+1}, and w_b(h₀, r_{b0}) ≡ 0 (mod D_b) by (T2), since x_b + h_{b0} = p.
- 𝓛 is any lattice inside the product ideal 𝓘 = (N, x+h)^{m1}·∏_b(D_b, w_b)^{m2_b} ⊂ Z[h, r],
  plus arbitrary Z[h, r]-multiples of the F_b (which σ kills).
  - Concretely, 𝓛 may be spanned by finitely many polynomials of the form
    (monomial in h, r) × N^{α}(x+h)^{s} × ∏_b w_b^{t_b}D_b^{m2_b − t_b}, with α + s ≥ m1 and
    0 ≤ t_b ≤ m2_b.
  - This covers the usual shift strategies built inside 𝓘 (basic; ideal-power variants), and window-only
    constructions (m1 = 0). It also covers any integer representatives g_b ≡ num_b (mod D_b).
  - All of them vanish mod p^{m1}∏D_b^{m2_b} at the root.
  - Let T be the largest total r-degree among the spanning polynomials, and
    d = T + max over their monomials h^a∏r_b^{u_b} of (a + Σ_b j_b u_b).
- Box: h ∈ [0, X] with X ≤ x/max(T, 1), so that (1 + X/x)^T ≤ e.
- Define σ: r_b ↦ ρ_b(h) := −N(−h_b)^{j_b+1}/(x+h). Note ρ_b(h₀) = r_{b0}, since x + h₀ = p.

**Prop T (PROVED).**

*Prior art (literature critic).*
- AASW 2018 §6.2, Thm 4 uses the same architecture: it maps a bivariate canonical polynomial (small
  inverse / Boneh–Durfee) to a univariate canonical polynomial for an unknown divisor, by substitution,
  then applies univariate optimality.
- Herrmann–May's unravelled linearization (ASIACRYPT 2009; PKC 2010) substitutes along the defining
  relation of a linearization variable, but to build lattices, not to bound them.
- AASW Lemma 3 (a Chebyshev sup-to-coefficient bound) is the prior analogue of step (iii).
- No prior formulation found by us of the specific statement: the product ideal with window moduli
  D_b, the rational-curve substitution with the curve ideal as kernel, and the vanishing-preserving
  range comparison (iv).

For every g ∈ 𝓛: For every g ∈ 𝓛:
- (i) G := (x+h)^T σ(g) ∈ ∏_b D_b^{m2_b}·(N, x+h)^{m1+T} ⊂ Z[h], with deg G ≤ d.
- (ii) G = 0 iff g vanishes identically on the curves r_b = ρ_b(h), i.e. g lies in the ideal of the
  F_b = (x+h)r_b + N(−h_b)^{j_b+1} over Q(h).
- (iii) Let Y_b ≥ max_{0≤h≤X}|ρ_b(h)| (e.g. Y_b = ⌈N/x⌉(X + |Δ_b|)^{j_b+1}). Suppose
  ‖g(hX, r·Y)‖₁ < p^{m1}∏D_b^{m2_b} and G ≠ 0. Then G′ := G/∏D_b^{m2_b} is a nonzero element of the
  univariate lattice L_{M,d} = (N, x+h)^M ∩ Z[h]_{≤d}, M = m1 + T, with G′(h₀) = 0,
  sup_{[0,X]}|G′| < e·p^M, and, for λ = 4d(d+1),

    ‖G′(hX/λ)‖₁ < e²(d+1)·p^M.

- (iv) **(Range comparison.)** Let μ = 2e²(d+1)^{3/2} and X″ = X/(λμ) = X/O(d^{3.5}). Suppose
  h₀ ∈ [0, X″]. Then ‖G′(hX″)‖₂ < p^M/√(d+1). So G′ passes the univariate Howgrave-Graham test for
  L_{M,d} at box X″, and exact-SVP univariate Coppersmith on L_{M,d} recovers h₀.
  - Hence, suppose a mixed lattice recovers every root in its box [0, X], via some HG-small g ∈ 𝓛
    with G ≠ 0. Then the exact-SVP univariate lattice L_{M,d} recovers every root in
    [0, X/O(d^{3.5})].
  - Kernel vectors cannot help even jointly (resultants, Gröbner bases). Their common zero set contains
    the curve r_b = ρ_b(h), which projects onto every h, so the elimination ideal in Q[h] is 0.

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
   - *Coefficient bound (Bernstein ellipse + Cauchy).* Map [0, X] to [−1, 1] by s = 2h/X − 1.
     - For a degree-d polynomial with |P| ≤ 1 on [−1, 1], Bernstein's lemma gives |P| ≤ ρ^d ≤ e
       inside the ellipse E_ρ with foci ±1 and ρ = 1 + 1/d.
     - The point of E_ρ closest to −1 is the vertex. For (a cos φ, b sin φ) with u = cos φ and
       a² − b² = 1, the squared distance to −1 is u² + 2au + 1 + b², increasing in u ∈ [−1, 1]. So E_ρ
       contains the disc of radius a − 1 = (ρ−1)²/(2ρ) = 1/(2d(d+1)) about −1.
     - In h this is the disc of radius X/λ about 0. Cauchy's estimates give |a_i|(X/λ)^i ≤ e·sup|P|,
       so ‖P(hX/λ)‖₁ ≤ e(d+1)·sup_{[0,X]}|P|.
     - This replaces the draft's valid but loose (3+3√2)^d.
   - G′(h₀) = 0 because g(h₀, r₀) = 0 over Z (Howgrave-Graham).
6. For (iv): write G′ = Σa_i h^i. From G′(h₀) = 0 and 0 ≤ h₀ ≤ X″ we get
   |a₀| ≤ Σ_{i≥1}|a_i|h₀^i ≤ Σ_{i≥1}|a_i|X″^i. Then
   ‖G′(hX″)‖₁ ≤ 2Σ_{i≥1}|a_i|(X/λ)^iμ^{−i} ≤ (2/μ)‖G′(hX/λ)‖₁ < 2e²(d+1)p^M/μ = p^M/√(d+1), and
   ‖·‖₂ ≤ ‖·‖₁. ∎

*Why (iv) is needed (codex topblock_2).* A bound on the norm alone would be vacuous: f^M ∈ L_{M,d}
already has norm (x+X)^M ≤ 2^M p^M but does not vanish at h₀. The vanishing at h₀ is what turns the
transfer into a range comparison.

**Consequences.**
- *(PROVED.)* Let ξ be the range exponent of a window-augmented lattice (X = N^ξ, recovering all roots
  in its box). Let ξ_U be that of exact-SVP univariate Coppersmith on L_{M,d}. Then
  ξ ≤ ξ_U + O(log d/log N) for ξ ≤ 1/2 − log T/log N (box condition X ≤ x/T, x ≍ √N).
  - This is o(1) whenever d = N^{o(1)}. That holds automatically for polynomial-time constructions,
    since a monomial h^a costs a·log₂X bits in the basis. Lattice dimension alone does not bound d.
  - With LLL in place of exact SVP, the univariate side loses LLL's 2^{O(d)} approximation factor. On
    the X6 grid (U(m, n), n ≤ 16) that loss was nil; X6 did not test the MIX transfer targets.
- *(OPEN.)* Whether ξ_U ≤ 1/4 + o(1) is the "optimality of Coppersmith" question for partial key
  exposure:
  - the determinant certificate says yes;
  - AASW say yes only under an assumption whose inequality is false on this family;
  - no unconditional proof is known.

  So Prop T is a **comparison theorem**: windows cannot beat the univariate lattice. It is not a
  proof that the univariate lattice stops at N^{1/4}.
- **Not covered** (both critics):
  1. lattices outside 𝓘 + (F_b), e.g. other moduli or ideals, saturations or divisions leaving 𝓘, and
     constructions with extra unknowns (q or μ, the near-divisor k, a block index u);
  2. certificates other than the Howgrave-Graham ℓ1/ℓ2 test, e.g. sublattice or Gröbner arguments on
     vectors that are not individually HG-small;
  3. boxes other than [0, X] (a centred box needs x − X in Y_b and changes the constants);
  4. near-divisor (shifted) windows, and rational centres before clearing denominators;
  5. non-lattice algorithms, such as integer-point enumeration on F_b = 0 or batched counting across
     blocks;
  6. the comparison is with exact SVP on L_{M,d}, not with LLL on it.
- *Adjacent-block windows* in X4 are centre shifts evaluated outside their own block. They are
  legitimate congruences, but not ladder queries.
- The same substitution with q ↦ N/(x+h) maps Coppersmith's integer-equation lattices for
  N = (x+h)(q₀+μ) into the same univariate ideal: W·h^aμ^b ↦ W·h^a(N − q₀(x+h))^b/(x+h)^b. This
  is consistent with the bivariate integer formulation's method guarantee also being N^{1/4}. That is a
HEURISTIC remark, not covered by Prop T. The known reduction from a univariate polynomial modulo an
unknown divisor to a bivariate integer polynomial is Blömer–May (EUROCRYPT 2005, LNCS 3494,
Thm 13). Coron (EUROCRYPT 2004; CRYPTO 2007) simplified the integer method but did not show this
domination.

### 3.2 The F-artifact (PROVED + X3)

- *Membership.* With m1, m2 ≥ 1, the polynomial F·N^{m1−1}D^{m2−1} lies in the mixed lattice. It
  equals row(t=1, e=1) + row(0, 0) plus the r-free correction N^{m1−1}D^{m2}(x+h)c(h) (the reduced-g
  identity of §1.1). That correction is in the lattice when E ≥ j+1.
- *Its length.* Its scaled norm is ≍ N^{m1}X^{j+1}D^{m2−1}, against the threshold x^{m1}D^{m2}. So it
  is HG-short iff X ≲ x(x/N)^{m1/(j+1)}.
- **For m1 ≥ 2 a shorter trivial vector dominates (Opus proof critic).** (x+h)^{m1−1}F·D^{m2−1} has
  norm ≍ N x^{m1−1}X^{j+1}D^{m2−1}. It is HG-short iff X ≲ (x^{j+2}/N)^{1/(j+1)} ≍ N^{j/(2(j+1))},
  independently of m1:
  - N^{1/3} for j = 2;
  - N^{3/8} for j = 3.
  - The critic verified this at 128 bits, j = 2, (2,1,10): short at log₂X = 28–40, where LLL returns
    only trivial good vectors.
- With m1 = 0 (windows only) the boundary of FD^{m2−1} is xN^{−1/(j+1)}.
- **Naive determinant counts therefore suggest mixed ranges up to about N^{1/3}, but every such vector
  is in ker σ.**
- *A tempting heuristic.* In units of log N, a window row of r-degree t contributes
  t((j+1)ξ − j/2) beyond the D-threshold, which is negative for ξ < j/(2(j+1)). This suggests that
  window rows "help".
  - That the gain is the trivial relation F is HEURISTIC: kernel vectors explain the observed
    misleading shortness.
  - Prop T(iv), not this heuristic, is what rules out a real gain.
- *X3, m1 = 1 cells.*
  - The largest X with any HG-good vector equals the largest X at which the F row is short (sparse
    three-term test) in all 18 cells, on the half-bit grid.
  - Zero nontrivial good vectors appear beyond U_img + 3 bits.
- *X3, m1 = 2 cells.* The pre-registered "F row" column tracks the wrong (longer) vector. The scans
  stopped 1 bit past X_max, with trivial good vectors still present, so the trivial extent was not
  measured.

### 3.3 Splitting the unknown (Herrmann–May; partly PROVED, partly measured)

- *Search problems.* Writing p = x₀ + Hu + v (block index u, offset v) turns the single contiguous
  unknown z = Hu + v into two. The integer box maps bijectively onto an interval, so the two search
  problems are equivalent.
- *(PROVED.)* With v′ = v + Hu + x₀ (unimodular, degree-preserving), the HM lattice u^a f^b N^{(m−b)+},
  a + b ≤ E, is the full truncated ideal (N, v′)^m in total degree ≤ E. It contains the univariate
  lattice in z.
- *(NOT proved.)* That its extra u-shifts never help. Equal point sets do not bound λ₁ of a strictly
  larger bivariate lattice, and the inverse map uses quotient and remainder, not a polynomial
  substitution.
- The generic HM bound γ₂ ≈ 0.2071 < 1/4 is a sufficient bound for general chunk positions.
- *X5 (measured).*
  - In all 12 cells the HM range is 2.0–3.5 bits *below* the best tested univariate lattice of no
    larger dimension (15 or 28).
  - All 2 769 HG-good HM vectors have a linear factor b(H_s u + v) + c, and recovery used it.
  - They are not all functions of z: codex found rows with (∂_u − H_s∂_v)g ≠ 0, e.g. u·(Hu + v − c).
  - In the tested cells, splitting only costs dimension.

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
- The mixed range tracks M = m1 + T and is insensitive to j and to n. For example, j = 2 with
  (2,1,10) and j = 3 with (2,1,11) give the same range.
- In the m1 = 1 rows, the "good" extent (≈ N^{1/2 − 1/(2(j+1))}, e.g. log₂X = 41.5 at 128 bits, j = 2)
  is pure kernel artefact.
- In the m1 = 2 rows, the last two columns are censored by the stopping rule (the scan ends 1 bit
  past X_max). The F-row column tracks a longer vector than the dominant trivial one, (x+h)F (§3.2).
- Each scan point redraws h₀ and hence x. Endpoints are protocol measurements on a half-bit grid,
  not exact monotone boundaries.

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

Decision rules apply to the tested parameters and dimensions only.

**Disclosures.**
- P2.2's bound was wrong as registered. The corrected statement (Prop LOC) was derived after
  seeing the 7 failures, and its 60/60 check is post hoc.
- The development peeks listed in the pre-registration were made before it was committed.
- *Not independent or not informative* (critics):
  - P3.2 and P4.2 follow from P3.1 and P4.1, because the U_img comparator is in the U_dim table with
    admissible dimension.
  - P1, P2.1 and P2.4 are implementation checks of proved identities.
  - P3.3's first half and P3.4 compare exact diagnostics of the forced F vector.
  - P6.1's lower side holds automatically.
  - P5.2's threshold and the linear-factor recovery were set after a disclosed peek.

## 6. Answers to the four questions

1. **Equivalence.**
   - Promised divisor search in a top block is *exactly* partial key exposure with
     (1/(j+2))log₂N top bits (Prop KB).
   - The window oracle implies it (Prop LOC; candidate-free at degree j+1). It is not known to be
     implied by it (Prop IV needs divisors of the N − k).
   - The window structure carries no extra information about p:
     - the Taylor structure is the identity N = q(x+h) (Lemma TS);
     - consecutive blocks are disjunctive;
     - multiple windows and degrees lie on trivial curves;
     - the only extra content is near-divisor data about N − k.
   - So the ladder's per-block task is exactly partial key exposure, and its oracle is at least as
     strong. For per-block (non-amortised) implementations its limit is the partial-key-exposure
     limit.
2. **Frontier.** §2. Polynomial time from ½ of p's bits; below that the cost is linear in the missing
   range. No structured setting applies to one N with contiguous top bits. Lower bounds are
   method-specific: the determinant certificate, Lemma MI, conditional AASW, and CHHS for known
   moduli. No general impossibility.
3. **Multi-window lattices.**
   - Within the product ideal and the HG certificate, they cannot beat the exact-SVP univariate ideal
     lattice by more than a poly(d) factor in range (Prop T(iv)).
   - Naive counts are inflated by trivial kernel vectors (§3.2).
   - Splitting into chunks is a reparametrisation of the search problem. Its lattice contains the
     univariate one, and measured, it is never better (§3.3).
   - Measured: X3–X5.
4. **Conclusion.** No mechanism.
   - The ladder's per-block task at degree j is exactly factoring from a fraction 2/(j+2) of p's top
     bits, and its oracle is at least as strong.
   - With known hint-using solvers the ladder costs N^{1/4+o(1)} for every j ≥ 2. The barrier sits
     at Coppersmith's half-the-bits.
   - That barrier is proved for determinant certificates and, via Prop T, transfers to every
     window-augmented lattice of the stated kind. It is not proved for general algorithms.
   - **Not excluded** (codex):
     - amortised or batched computation across blocks;
     - non-lattice integer-point methods on the curves F_b = 0;
     - near-divisor cancellation tricks outside the ideal.

     Any per-block success would itself improve partial key exposure (Prop KB).

## 7. Reviews

**GPT-6.1-sol blind round (`codex/topblock_1`).** Same core verdict, reached independently:
- promised search = partial key exposure;
- the window oracle implies it, with a candidate loss;
- the reverse reduction is UNVERIFIED;
- consecutive blocks are disjunctive.

Contributed:
- the degree-(j+1) candidate-free query;
- the identity d·r = (−1)^j N h^{j+1} + Dk;
- the checked frontier references (Coppersmith Thm 4, May Thm 1, Herrmann–May γ_k, May–Ritzenhofen,
  BDH, Heninger–Shacham, AASW, CHHS, Cohn–Heninger, Harvey–Hittmeir).

**GPT-6.1-sol adversarial round (`codex/topblock_2`): OVERCLAIMED.** All fixes applied:
- Prop T's norm-only transfer was vacuous (f^M), so (iv), the vanishing-preserving range comparison,
  was added;
- the scope of "every shift strategy" was narrowed;
- the reduced-g F identity;
- the AASW lattice and Assumption 1 identified, and the posited inequality shown false on this family;
- the short vector is F^m, not N^m;
- "functions of z" withdrawn (HM);
- the LOC parity details;
- the IV window-position term;
- EQ(c)/(d) restricted;
- the hint-ignoring-solver arithmetic (j ≥ 19);
- P3.2/P4.2 dependence disclosed;
- scan-protocol caveats;
- unexcluded mechanisms listed.

**Opus proof critic: no FATAL; 4 MAJOR (all fixed).**
- EQ(c)/(d) overclaims;
- the exponent consequence of Prop T, replaced by Bernstein-ellipse + Cauchy + the a₀ argument, which
  the critic re-checked as correct;
- the m1 ≥ 2 trivial vector (x+h)^{m1−1}F D^{m2−1} with extent N^{j/(2(j+1))};
- the §2.3 status / "provably" overclaim.

Also the HM "PROVED" scope, and minor items (h = 0, odd-j representative, ν₀, |R|, Prop IV length, d
definition, table reading, tautological predictions).

Its independent checks:
- Prop LOC: 394 cases, 0 violations;
- degree-(j+1): 294 blocks, 0 false candidates;
- Prop T(i): exact ideal membership of every basis row in 12 configurations;
- all tables recomputed.

Our check of T(iii) on 675 nontrivial HG-good mixed vectors: 0 violations, max ratio 0.003.

**Opus literature critic (primary texts checked): fixes applied.**
- Lemma MI relabelled KNOWN (AASW §4.1–4.2);
- AASW Thm 1 / failure definition, and interpretation reworded ("cannot be applied", not "refuted");
- Prop T prior art: AASW §6.2 Thm 4 (same architecture), Herrmann–May unravelled linearization,
  AASW Lemma 3;
- Blömer–May 2005 Thm 13 replaces "Coron's equivalence";
- Cohn–Heninger Thm 1.2/§4.2 for the Johnson analogy;
- Lenstra 1984 / CHN 2008 as the unknown-divisor list-size results and the low-bit j = 1, 2 rungs;
- Maurer's oracle factoring (context for EQ);
- CHHS Thm 3 (integer-valued polynomials);
- Herrmann–May Thms 3–5 and running time;
- Harvey 2021 / Harvey–Hittmeir wording;
- Khomovsky as prior art for Taylor windows.

**Outside this note (not edited; for the orchestrator).** `notes/analytic_results.md:46` says Hiary
needs "1/3 of the bits of p". It should read "(1/3)·log₂N bits, i.e. 2/3 of p's bits" (literature
critic).

## 8. Reproduction

```bash
.venv/bin/python -m pytest -q tests/test_topblock.py
for e in x1 x2 x3 x4 x5 x6 summary; do PYTHONPATH=. .venv/bin/python -m src.topblock.run_experiments $e; done
# x3 ~25 min, x4/x5/x6 a few minutes each; fpylll needs cysignals in .venv
```

