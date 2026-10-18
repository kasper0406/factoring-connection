# dcltheorem: a deterministic N^{1/6+ε} algorithm for almost all balanced semiprimes

Agent `dcltheorem`, 2026-10-18. Pre-registration: `notes/dcltheorem_prereg.md` (commit
`3550e24`, before any pre-registered run). Code: `src/dcltheorem/`. Tests:
`tests/test_dcltheorem.py`. Raw data: `experiments/raw/dcltheorem_*.jsonl`. Summary:
`experiments/analysis/dcltheorem_summary.json`. Inputs: notes/multimod_results.md (Thm D) and
the GPT-6.1-sol audit notes/codex_r12_dcl_audit.md.

## 0. Verdict (plain language)

**Thm D can be completed into a real deterministic factoring algorithm for almost all balanced
semiprimes, in the ordinary bit-complexity model.** The previous "base order" gap closes, and it
needs no order lemma for β^L at all.

- **Theorem D′ (proved; three independent reviews, all ACCEPT WITH FIXES, fixes applied:
  GPT-6.1-sol dclt_1, an Opus critic, and the fresh GPT-6.1-sol r13 referee).** For each
  ε ∈ (0, 1/60] and all n ≥ n₀(ε) there is a deterministic algorithm with these properties:
  - it sees only N and costs N^{1/6+5ε/6}·(log N)^{O(1)} bit operations on **every** n-bit
    input;
  - it never outputs anything but a proper divisor (1 < d < N, d | N) or FAIL. Every candidate
    is validated before output; this was added after r13 found improper outputs for tiny
    prime N;
  - it factors all but an O(N^{−ε}) fraction of the **bounded-ratio, distinct-prime**
    semiprimes S_n = {pq : p < q ≤ 2p primes, n bits}, with exceptions **counted by
    pairs**.
- **Corollary D″.** N^{1/6}(log N)^{O(1)} cost with an O((log N)^{−A}) exceptional fraction,
  for any fixed A.
- **How the order gap closes.**
  - Mod p the true collision is an exact identity, so no order is needed.
  - Mod q one needs only that the true giant does *not* collide. Running four disjoint prime
    bands guarantees a band coprime to q − 1. On that band the collision condition becomes
    ord_q(2) | (A−B)(p−1) − ℓ²δ.
  - That condition is rare on average: an elementary count of q with small ord_q(2), plus a CRT
    count over p.
- **What is not proved.**
  - The worst case is open. It needs DCH coverage, large order mod the specific q, and the
    exclusion of all-FULL pairs.
  - Ratios q/p > 2 are not covered.
  - A Poisson-quality exceptional rate is open. The large sieve's n²/H² is why N^{−ε} costs
    N^{5ε/6}.
  - The statement is about the true factor averaged over S_n, not full-candidate coverage for
    every N.
- **Experiment (pre-registered, theorem's own design, end to end).**
  - The algorithm and all lemmas held: 0 wrong outputs and 0 lemma violations in 656 DCL runs (240 lemma + 416 operational; corrected 2026-10-18 from "896", which double-counted the lemma runs), with
    3002 clean covering targets checked.
  - Measured arithmetic-cost slope over 48–80 bits: **0.203 against calibrated operational HH's
    0.226**.
  - **But HH is still 2.3–4.9× cheaper at 40–80 bits.** No practical gain at these sizes, and
    none at heuristically full coverage (C2 slope 0.227).
- **Novelty.** Every ingredient is prior art (Lehman, McKee, HH Prop. 3.3, Harvey Alg. 4.1, the
  large sieve).
  - **McKee–Pinch 1996 §5 is a direct precursor.** Its divisibility event is λp ≡ ±μq (mod k)
    with small multipliers, and it gives a *heuristic* N^{1/4+ε} balanced method over primes
    k ~ N^{1/4} (§10).
  - D′ adds: the band near N^{1/6}; one shared cofactor-powered baby list; the prime-pair
    large-sieve exceptional bound; the clean-band one-sidedness argument; rigour.
  - We found no prior *proved* sub-1/5 almost-all result. Carella (arXiv 1308.2891) *claims*
    deterministic 1/6; the r13 referee judges its proof incomplete (§10).

## 1. Population, target and statement

Throughout, n is the bit length, X = 2^n and P = 2^{n/2}.

**Population.**

  S_n = {(p, q) : p < q ≤ 2p primes, 2^{n−1} ≤ pq < 2^n}.

It is counted by pairs, equivalently by semiprimes N = pq. Then p ∈ (P/2, P) and
q ∈ (P/√2, √2 P). Counting with the PNT in [0.85X′, 0.9X′] × [1.2X′, 1.25X′],
X′ = 2^{(n−1)/2}, gives |S_n| ≥ c₀ P²/n² with an absolute c₀ > 0 for n ≥ n₀.

**Target.** A deterministic algorithm 𝒜_{n,ε} with the following properties.
- It receives N (n bits) and outputs either a nontrivial divisor of N or FAIL.
- Its worst-case cost on every n-bit input is N^{1/6+O(ε)} bit operations.
- It outputs FAIL on at most O(N^{−ε})·|S_n| pairs of S_n.

**Theorem D′ (proved; reviews in §7).** Fix ε ∈ (0, 1/60]. There are n₀(ε) and absolute
constants C, c such that for n ≥ n₀ the algorithm DCL-det_{n,ε} of §2 has the following
properties.
- (i) **Correctness.** On every input it outputs a proper divisor d (1 < d < N, d | N) or
  FAIL. Every candidate is validated before output (`algorithm._proper`).
  - Inputs N < 4 return FAIL; even N ≥ 4 return 2; a band prime ℓ is returned only if
    1 < ℓ < N.
  - (Correction, r13: the pre-fix code returned 2 for N = 2 and ℓ for N = ℓ prime.)
  - An exhaustive regression test covers all 2 ≤ N ≤ 10⁴ under two designs.
- (ii) **Cost.** On every n-bit input it uses O(n^c·2^{(1/6+5ε/6)n}) bit operations, i.e.
  N^{1/6+5ε/6}(log N)^{O(1)}.
- (iii) **Success.** It outputs FAIL on at most C·2^{−εn}·|S_n| pairs (p,q) ∈ S_n.

Equivalently: for every δ ∈ (0, 1/72], cost N^{1/6+δ+o(1)} with an exceptional fraction
O(N^{−6δ/5}).

**Corollary D″ (proved, same proof; polylog regime).** For every fixed a > 0, the same
algorithm with H = ⌈n^{1+a}⌉ (and T = 2^{n/4}n^{−1−a} in the analysis) has cost
N^{1/6}(log N)^{O(1+a)} and fails on O((log N)^{−2a}) of S_n.

**Quantifier (audit item 5).**
- Theorem D′ is a statement about the true factor of an input drawn from S_n, counted over
  pairs.
- It does **not** say that for almost every N all integer candidates p are covered, so it
  does not show that the full-candidate covering model has exponent 1/6 for almost
  every N.
- It does not cover the worst case (§6), nor ratios q/p > 2.

## 2. The algorithm DCL-det_{n,ε} (operational; `src/dcltheorem/algorithm.py`)

### 2.1 Parameters

All parameters are functions of n and ε only; nothing depends on N.

- H = ⌈n·2^{εn/2}⌉. The code (`params.theorem_design`) uses n·⌊2^{εn/2}⌋ and ⌊H^{1/3}⌋ in Y. These are constant-factor conventions, which the proof allows.
- V = {(A,B) ∈ [1,H]² : gcd(A,B) = 1, A ≠ B}. Then H²/2 ≤ |V| ≤ H² for large H.
- Y = ⌈2^{n/6}H^{−1/3}⌉.
- For k = 0, 1, 2, 3: 𝒫_k = {ℓ prime : ℓ ≡ 3 (mod 4), 2^kY < ℓ ≤ 2^{k+1}Y} and
  L_k = ∏_{ℓ∈𝒫_k} ℓ².
- m̄_k = ⌊span̄/(2^kY)²⌋ + 2, with span̄ = ⌊(5/2)·H·(2^{n/2}+1)⌋ + 2H + 3. This is an upper
  bound for every window width in block k (§2.2).
- The base is β = 2.

### 2.2 Steps

Input N (n bits). If N is even, output 2. Set a₀ = ⌈√(N/2)⌉ and b₀ = ⌊√N⌋. Every (p,q) ∈ S_n
has p ∈ [a₀, b₀].

For v = (A,B) ∈ V put:
- Lo_v = ⌊2√(ABN)⌋. By AM–GM, Lo_v ≤ F_v(x) := Bx + AN/x for all x > 0.
- Hi_v = max_{x∈{a₀,b₀}} ⌈F_v(x)⌉. By convexity, F_v ≤ Hi_v on [a₀, b₀].
- Then Hi_v − Lo_v ≤ (B + √2A)√N + 2H + 2 ≤ span̄.

For each block k = 0, 1, 2, 3:

1. **Exponent tree.**
   - Build the balanced integer product tree of {ℓ² : ℓ ∈ 𝒫_k}.
   - Push 2 down the tree: a node holding 2^{L/∏node} hands 2^{L/∏node·∏sibling} to each
     child. This yields u_ℓ = 2^{L_k/ℓ²} mod N for all ℓ.
   - Set γ = u_{ℓ₀}^{ℓ₀²} = 2^{L_k} mod N.
   - If ℓ | N for some ℓ, output ℓ.
2. **Targets.** For ℓ ∈ 𝒫_k and v ∈ V:
   - a := A N B^{−1} mod ℓ. If a^{(ℓ−1)/2} ≢ 1 (mod ℓ), skip.
   - Otherwise r₀ = a^{(ℓ+1)/4} mod ℓ, and the roots of Br² ≡ AN (mod ℓ) are {r₀, ℓ−r₀}.
   - For each root r:
     - s = N r^{−1} mod ℓ² and c = (Br + As) mod ℓ²;
     - g_lo = ⌈(Lo_v − c)/ℓ²⌉ ≥ 0 and g_hi = ⌊(Hi_v − c)/ℓ²⌋;
     - W_τ = g_hi − g_lo + 1 ≤ span̄/ℓ² + 1 ≤ m̄_k.
   - The target is τ = (v, ℓ, r). Write T_k for the number of targets: T_k ≤ 2|V||𝒫_k|.
3. **Baby length.** Either m = m̄_k (the theorem's default), or the implementation's
   *endpoint-competitive* rule: minimise m + Σ_τ⌈W_τ/m⌉ over an integer geometric grid
   that ends at W_max = max_τ W_τ ≤ m̄_k.
   - The success proof needs only m ≤ m̄_k.
   - The **cost** proof needs one of these two rules. Under either, m + T′ ≤ m̄_k + T_k,
     with T′ = Σ_τ⌈W_τ/m⌉.
   - *Correction (GPT-6.1-sol dclt_1).* An earlier draft said "the proof needs only
     m ≤ m̄_k". For cost that is false: m = 1 on prime-square inputs gives
     T′ ≫ H²P/(Y log Y).
4. **Giants.**
   - Precompute w_ℓ = u_ℓ^N for each ℓ.
   - For each target and 0 ≤ j < ⌈W_τ/m⌉:

       G_{τ,j} = w_ℓ^A · u_ℓ^{B−c+ℓ²} · γ^{−(g_lo+1)} · γ^{−jm}  ≡  u_ℓ^{B+AN−c}·γ^{−(g_lo+jm)}  (mod N).

     The congruence holds because u_ℓ^{ℓ²} = γ. γ^{−1} exists because gcd(2,N) = 1.
5. **Collision detection.**
   - f(X) = ∏_{0≤i<m}(X − γ^i) ∈ (Z/N)[X], built by a product tree.
   - f(G_{τ,j}) for all giants, by a remainder tree. Only monic divisors are used, so this
     is valid over Z/N.
   - For each giant, d = gcd(f(G), N). If 1 < d < N, output d.

Output FAIL if no block yields a divisor.

**Remarks.**
- (a) No order is computed, assumed or certified. No randomness is used.
- (b) The symmetric weight (1,1) is excluded. For A = B the true giant always collides mod
  q as well (Lemma 4 with A − B = 0), so it can never be one-sided.
- (c) The four bands are needed only for Lemma 6, which guarantees a "clean" band.

## 3. Lemmas

Fix (p,q) ∈ S_n, N = pq, a block k, and write L = L_k, u = u_ℓ and d := ord_q(2).

**Lemma 1 (class lock; HH 2022 Prop. 3.3).** Let ℓ ∈ 𝒫_k, (A,B) ∈ V and r a root, and
suppose p ≡ r (mod ℓ). Then g := (Bp + Aq − c)/ℓ² ∈ Z.

*Proof.*
- ℓ > Y > H ≥ A, B, and ℓ < p, so ℓ ∤ ABN. Write p = r + ℓx.
- Mod ℓ²: p^{−1} ≡ r^{−1}(1 − ℓxr^{−1}), so q = Np^{−1} ≡ s − ℓx·s r^{−1}.
- Hence Bp + Aq ≡ c + ℓx(B − A s r^{−1}).
- s r^{−1} ≡ N r^{−2} (mod ℓ) and Br² ≡ AN, so B − ANr^{−2} ≡ 0 (mod ℓ). The last term is
  therefore ≡ 0 (mod ℓ²). ∎

**Coverage criterion.** p mod ℓ is one of the computed roots ⇔ Bp² ≡ AN (mod ℓ) ⇔
ℓ | Bp − Aq (cancel p, since ℓ ∤ p). We say τ = (v, ℓ, p mod ℓ) *covers* p.

**Lemma 2 (window).** If τ covers p, then g_lo ≤ g ≤ g_hi. Put i_full = g − g_lo,
j* = ⌊i_full/m⌋ and i* = i_full − j*m ∈ [0, m).

*Proof.* p ∈ [a₀, b₀], so Lo_v ≤ F_v(p) = Bp + Aq ≤ Hi_v. ∎

**Lemma 3 (mod p).** If τ covers p, then G_{τ,j*} ≡ γ^{i*} (mod p). No hypothesis on any
order is needed.

*Proof.*
- B + AN − (Bp + Aq) = (p − 1)(Aq − B), and 2^{p−1} ≡ 1 (mod p).
- So u^{B+AN−c} ≡ 2^{(L/ℓ²)(Bp+Aq−c)} = 2^{Lg} = γ^g (mod p).
- Multiply by γ^{−(g_lo + j*m)}. ∎

**Lemma 4 (mod q).** If τ covers p, then G_{τ,j*} ≡ γ^{i*}·u^{(A−B)(p−1)} (mod q).

*Proof.*
- Mod q − 1: N ≡ p, so B + AN − c ≡ B + Ap − c.
- Also ℓ²g = Bp + Aq − c ≡ Bp + A − c.
- The difference is (A − B)(p − 1). Multiplying by L/ℓ² preserves the congruence mod q − 1,
  and 2^{q−1} ≡ 1 (mod q). ∎

**Lemma 5 (one-sided extraction).** Call the true giant G = G_{τ,j*} *FULL* if
G ≡ γ^{i′} (mod q) for some i′ ∈ [0, m). If τ covers p and G is not FULL, then
gcd(f(G), N) = p.

*Proof.*
- f(G) = ∏_i (G − γ^i).
- The factor i = i* vanishes mod p (Lemma 3).
- Mod q, no factor vanishes, by definition, and Z/q is a field.
- So p | f(G) and q ∤ f(G). ∎

**Lemma 6 (clean band).** Call block k *clean for q* if no ℓ ∈ 𝒫_k divides q − 1. If
Y⁴ > √2P, at least one of the four blocks is clean for q.

*Proof.* Every band prime exceeds Y. Four distinct band primes dividing q − 1 would give
q − 1 > Y⁴ > √2P > q. So at most 3 blocks contain a prime factor of q − 1. ∎

Asymptotically: Y⁴ ≥ 2^{2n/3}H^{−4/3}, and 2^{n/6} > √2 H^{4/3} = √2 n^{4/3}2^{2εn/3}
for ε < 1/4 and n large.

**Lemma 7 (FULL in a clean block).** If block k is clean for q, then the true giant is FULL
iff d | (A−B)(p−1) − ℓ²δ for some δ ∈ [−i*, m−1−i*] ⊂ (−m, m).

*Proof.*
- By Lemma 4, FULL ⇔ u^{(A−B)(p−1)} ≡ γ^{δ} (mod q) with δ = i′ − i*.
- Equivalently 2^{(L/ℓ²)[(A−B)(p−1) − ℓ²δ]} ≡ 1 (mod q).
- Since d | q − 1 and gcd(L, q − 1) = 1, gcd(d, L/ℓ²) = 1. ∎

**Lemma 8 (small orders are rare; elementary).** For Q₀ ≥ 2 and T ≥ 1,

  #{q > Q₀ prime : ord_q(2) ≤ T} ≤ T(T+1)/(2 log₂ Q₀).

*Proof.* Such a q divides 2^k − 1 for some k ≤ T. Since 2^k − 1 < 2^k, it has fewer than
k/log₂Q₀ prime factors > Q₀. Sum over k ≤ T. ∎

**Lemma 9 (coverage: the arithmetic large sieve, in the audit's repaired form).** Fix q and
k. Let U_k(q) be the set of integers p ∈ [P/2, P] that no target of block k covers, i.e.
ℓ ∤ Bp − Aq for all (v, ℓ) ∈ V × 𝒫_k. For n ≥ n₀,

  |U_k(q)| ≤ (P/2 + 1 + 4^{k+1}Y²)/J_k,  J_k ≥ |V|·Σ_{ℓ∈𝒫_k} 1/ℓ ≥ c₁H²/n.

*Proof.*
- For ℓ ∈ 𝒫_k, ℓ ∤ q, so ℓ | Bp − Aq ⇔ p ≡ AqB^{−1} (mod ℓ).
- These |V| classes are distinct: AqB^{−1} ≡ A′qB′^{−1} implies ℓ | AB′ − A′B, and
  0 < |AB′ − A′B| < H² < Y < ℓ. This uses Y > H², i.e. 2^{n/6} > n^{7/3}2^{7εn/6}, true for
  ε < 1/7 and n large.
- So U_k(q) avoids ω(ℓ) = |V| < ℓ classes mod every ℓ ∈ 𝒫_k.
- The arithmetic large sieve (Montgomery 1968; Montgomery–Vaughan 1973, constant
  N₀ + Q² for an interval of N₀ integers) with Q = 2^{k+1}Y gives
  |U| ≤ (N₀ + Q²)/J, J = Σ_{d≤Q} μ²(d)∏_{ℓ|d} ω(ℓ)/(ℓ − ω(ℓ)).
- Only the prime terms are needed, and they follow directly. For each ℓ ∈ 𝒫_k,
  Parseval and Cauchy–Schwarz give Σ_{a=1}^{ℓ−1}|S(a/ℓ)|² ≥ ω/(ℓ−ω)·|U|². The fractions
  a/ℓ are Q^{−2}-separated, so the additive large sieve applies (referee dclt_1).
- Keeping only d = 1 and d = ℓ ∈ 𝒫_k: J ≥ |V| Σ_{ℓ∈𝒫_k} 1/ℓ.
- By the PNT in the progression 3 mod 4 (effective for the fixed modulus 4),
  Σ_{ℓ∈𝒫_k} 1/ℓ ≥ (1/2 − o(1))·ln 2/ln(2^{k+1}Y) ≥ c′/n.
- |V| ≥ H²/2 for large H, since the coprime pairs number (6/π² + o(1))H². ∎

**Lemma 10 (FULL is rare).** Fix q with d = ord_q(2) > T and a block k clean for q. Let
F_k(q) be the set of integers p ∈ [P/2, P] such that some (v, ℓ) ∈ V × 𝒫_k has ℓ | Bp − Aq
and d | (A−B)(p−1) − ℓ²δ for some |δ| < m̄_k. Then

  |F_k(q)| ≤ Σ_{v,ℓ,|δ|<m̄_k} (P·H/(2·2^kY·T) + 1)
           ≤ 2|V||𝒫_k| m̄_k·(PH/(2^{k+1}YT) + 1).

*Proof.*
- Fix (v, ℓ, δ) and let g′ = gcd(A − B, d) ≤ H (A ≠ B, |A − B| < H).
- The congruence (A−B)(p−1) ≡ ℓ²δ (mod d) is either unsolvable, or holds exactly on one
  class of p mod d/g′.
- ℓ ∤ d, because clean means ℓ ∤ q − 1. With p ≡ AqB^{−1} (mod ℓ), CRT gives one class mod
  ℓd/g′ ≥ ℓd/H.
- An interval of length P/2 contains ≤ P H/(2ℓd) + 1 such integers. ∎

## 4. Proof of Theorem D′

**(i) Correctness.** Only gcds with N in (1, N), or band primes dividing N, are output.

**(ii) Cost (bit operations; M(b) = O(b log b) multiplication; per block, with T_k
targets and m ≤ m̄_k).**

| step | operations | bit cost |
|---|---|---|
| prime sieve to 16Y, V, Lo_v/Hi_v (isqrt of O(n)-bit integers) | O(Y log log Y) + O(H²) isqrt | Õ(Y + H²n) |
| integer product tree of ℓ² (Σ bits = log₂L_k ≤ 3·2^kY) | balanced tree | O(M(Y)·log Y) |
| exponent tree (per level Σ exponent bits ≤ log₂L_k + #nodes; ⌈log₂|𝒫_k|⌉ levels) + γ, γ^{−1} | O(Y log Y) mults mod N + 1 inversion | O(Y log Y·M(n)) |
| roots (Euler test, powering mod ℓ), s, c, windows | O(\|V\|\|𝒫_k\| log Y) mults mod ℓ, ℓ² | Õ(H²Y) |
| w_ℓ = u_ℓ^N; giants (exponents of O(n) bits) | O(n\|𝒫_k\| + nT_k + T′) mults mod N | O((n(T_k + Y) + T′)·M(n)) |
| choice of m (endpoint-competitive grid: O(n) grid points × T_k targets); inversions B^{−1} mod ℓ and r^{−1} mod ℓ²; reductions N mod ℓ, ℓ² | O(nT_k) word/small-integer ops | Õ(T_k n^{O(1)}) |
| babies γ^i | m mults | O(m M(n)) |
| f by product tree; remainder tree at the T′ = Σ_τ⌈W_τ/m⌉ ≤ T_k + Σ_τ W_τ/m giants | O((M(m) + M(T′)) log(m+T′)) ring ops; Kronecker substitution makes a degree-D product O(M(D·(2n + log D))) bits | Õ((m + T′)n) |
| gcds | T′ gcds of n-bit numbers | O(T′ M(n) log n) |

With m = m̄_k (or any m that minimises the list count) the list is
m + T′ ≤ m̄_k + T_k + 1 ≤ 2.5HP/Y² + 4 + 2H²·2^{k+1}Y. Over the four blocks:

  total = O(n^{O(1)}·(HP/Y² + H²Y)) = O(n^{O(1)} H^{5/3} P^{1/3}),

because Y ≍ P^{1/3}H^{−1/3}. With H = n·2^{εn/2} this is O(n^{O(1)}·2^{(1/6+5ε/6)n}).

**Arithmetic model (referee dclt_1, item 4).**
- The bound assumes quasi-linear integer multiplication M(b) = O(b log b) (Harvey–van der
  Hoeven; Schönhage–Strassen gives an extra log log). This applies to the exponent product
  tree, the Kronecker-substituted polynomial products and the gcds.
- Kronecker substitution over composite N: coefficients lie in [0, N), and the slots are
  wider than 2n + ⌈log₂(D+1)⌉.
- Monic remainders: reversal plus Newton inversion of a power series with constant term 1.
  No inverse of a nonunit is ever needed.
- The experimental metric E (§8) is a *proxy*. It is not the proof of (ii). The Python
  implementation uses GMP integer products (gmpy2) and flint polynomial arithmetic. Its
  `powN` charges binary square-and-multiply counts.

**(iii) Exceptional set.**
- Let T = ⌊2^{(1/4−ε)n}⌋. For (p,q) ∈ S_n, let k(q) be the least block clean for q
  (Lemma 6). If
  - (a) d = ord_q(2) > T,
  - (b) p ∉ U_{k(q)}(q), and
  - (c) p ∉ F_{k(q)}(q),

  then DCL-det succeeds.
- Indeed by (b) some target τ of block k(q) covers p (Lemma 1, coverage criterion).
- By (c) and Lemma 7, its true giant is not FULL: the exact δ-range lies inside |δ| < m̄_k.
- By Lemma 5 that giant's gcd is p. The algorithm outputs a nontrivial divisor at block k(q)
  at the latest. It may output q, or output earlier from another block; both are fine.
- So the FAIL pairs lie in E_a ∪ E_b ∪ E_c, with E_a = {ord_q(2) ≤ T},
  E_b = {d > T, p ∈ U_{k(q)}(q)} and E_c = {d > T, p ∈ F_{k(q)}(q)}.

The counts sum over prime q ∈ (P/√2, √2P), of which there are ≤ 3P/n for n large:

- **|E_a|.** By Lemma 8 with Q₀ = 2^{(n−1)/2}, each bad q has ≤ P/2 + 1 partners, so
  |E_a| ≤ (T(T+1)/(n−1))·(P/2+1) ≤ 2^{(1−2ε)n}/n·(1 + o(1)). Then
  |E_a|/|S_n| ≤ (n/c₀)·2^{−2εn}(1 + o(1)).
  - *Correction (dclt_1):* T(T+1) replaces T². The estimate is unchanged.
- **|E_b|.** By Lemma 9, Σ_q |U_{k(q)}(q)| ≤ (3P/n)·(P/2 + 1 + 256Y²)·n/(c₁H²) ≤ C P²/H².
  So |E_b|/|S_n| ≤ C n²/H² ≤ C 2^{−εn}. This is the audit's repaired normalisation: one
  fixed interval, integers counted, no short-interval prime counts.
- **|E_c|.** By Lemma 10 with |V| ≤ H², |𝒫_k| ≤ 2^kY and m̄_k ≤ 3HP/(4^kY²):
  |F_k(q)| ≤ 6H⁴P²/(Y²T) + 6H³P/Y. So

    |E_c|/|S_n| ≤ (3P/n)(6H⁴P²/(Y²T) + 6H³P/Y)·n²/(c₀P²) ≤ C n (H^{14/3}P^{1/3}/T + H^{10/3}/P^{1/3}),

  using Y² ≥ P^{2/3}H^{−2/3}. Here H^{14/3}P^{1/3}/T ≤ (2n)^{14/3}2^{(−1/12 + 7ε/3 + ε)n} and
  H^{10/3}/P^{1/3} ≤ n^{4}2^{(5ε/3 − 1/6)n}.
- **Conclusion.** For ε ≤ 1/60, −1/12 + 10ε/3 ≤ −1/36 < −ε − (small). So |E_c|/|S_n| =
  O(2^{−εn}) for n ≥ n₀(ε), and the total is O(2^{−εn}) = O(N^{−ε}). ∎

**Corollary D″.**
- Take H = ⌈n^{1+a}⌉ and T = ⌊2^{n/4} n^{−1−a}⌋.
- (a) gives O(n^{−2a−1}) and (b) gives O(n^{−2a}).
- (c) gives n^{O(1)}·2^{−n/12}·n^{1+a}·n^{O(1)} = o(n^{−2a}).
- The cost is n^{O(1)}H^{5/3}2^{n/6}. ∎

**Where the constants come from.**
- ε ≤ 1/60 is used only in (c) (the requirement is −1/12 + 13ε/3 < 0, i.e. ε < 1/52) and in
  the scale separations Y > H² (ε < 1/7) and Y⁴ > √2P (ε < 1/4).
- The exceptional rate N^{−ε} comes from the large sieve (b), which is the binding term.
- The Poisson heuristic predicts e^{−λ} with λ ≍ H²/n. That would give N^{−ε} already at
  H ≍ n√ε, i.e. polylog cost. The large sieve is the lossy step.

## 5. Answers to the audit (notes/codex_r12_dcl_audit.md) and to the work items

| audit item | resolution here | status |
|---|---|---|
| (1) normalisation of Thm D's exceptional set (short I_q, prime counts in I_q) | Lemma 9 + §4(iii) E_b: one fixed interval [P/2, P], all integer survivors counted, summed over prime q, divided by \|S_n\| ≍ P²/n² | proved |
| (2) deterministic roots | band primes ℓ ≡ 3 (mod 4), root a^{(ℓ+1)/4}; sieve loses a constant (PNT mod 4) | proved |
| (3) bit complexity: balanced integer product tree; uncharged Bluestein powers | balanced tree (GMP products) in `all_but_one`; full cost table §4(ii). The HH baseline's Bluestein chirp powers are now charged (`stage3.harvey_find_collisions`) | proved (table); implemented |
| (4) **base order of γ = β^L** (the substantive gap; 1333 example) | **Closed without any order lemma for β^L.** No order is needed mod p (Lemma 3 is an identity). Mod q, the four disjoint bands guarantee a clean band (Lemma 6), so ord_q(γ) = ord_q(2) (Lemma 7). Two-sided true giants are counted directly (Lemma 10), and small ord_q(2) is rare by an elementary count (Lemma 8). Extra failures are O(2^{−εn}), not o(1). The 1333 example (L built from 3, 5, 7 \| p−1, q−1) cannot occur in a clean band, and an L-power image of size 2 mod **p** is harmless here. | proved (+ GPT-6.1-sol referee dclt_1: "closes the previous order/extraction gap") |
| (5) quantifier mismatch | §1: D′ concerns the true factor counted over S_n. It does not give full-candidate coverage for almost every N | stated |
| (6) prior art | §10 | McKee–Pinch §5 read (r13 scans): direct precursor (heuristic N^{1/4}); the combination is claimed only as "no prior formulation found by us"; Carella 1308.2891 claim noted with the r13 assessment |
| (7) experiments used greedy Hmax=12, not Thm D's design | §8 runs the theorem's algorithm (4 bands, ℓ ≡ 3 mod 4, base 2, fixed V(H), n-only Y) with an n-only H-schedule (λ₀ ≈ 3, i.e. H ≍ √n). The bridge: D′'s algorithm is defined for any n-only (H, Y), and the schedule H ≍ √n lies *below* Corollary D″'s range (a > 0, H = n^{1+a}). There the proved bound n²/H² is vacuous, so the experiment tests the algorithm, the lemmas and the cost, not the proved rate. The heuristic rate e^{−λ} is what the data reflect. | disclosed |

Work items not covered above:
- **Small-factor removal (item 3).** Unnecessary. On S_n both factors exceed P/2 > 16Y ≥ every
  band prime, and the target is only S_n. Off the population the algorithm may FAIL (allowed),
  and band primes dividing N are output.
- **One-sided collisions (item 2, last bullet).** These are Lemmas 4, 5, 7 and 10. The weight
  (1,1) is excluded because it is always two-sided.
- **Restricting to M ∤ p−1 (item 2).** Not needed. Cleanliness is needed only for q − 1, and is
  handled by trying 4 bands rather than by an exceptional set. Order mod p is irrelevant.
- **Harvey–Hittmeir 2026 / Oznovich–Volk (item 2, first bullet).** Not used. The audit notes they
  produce an arbitrary large-order element, not one in the L-power image. The clean-band argument
  makes the L-power harmless mod q and irrelevant mod p.

## 6. Where the result stops (precise remaining gaps)

1. **Worst case (every N in S_n) is NOT proved.** Three independent obstacles remain:
   - (a) **coverage** of every p, which is multimod's DCH (conjectured; the large sieve gives
     only "almost all");
   - (b) **q with small ord_q(2)** (≤ N^{1/4−ε}). The proof simply discards them (E_a). A
     worst-case version would need a base with large order mod the specific prime q *on a
     clean band*. Harvey–Hittmeir 2026 gives a large-order element mod N, but large order
     mod N does not give large order mod q;
   - (c) pairs where **every** covering true giant is FULL (E_c). These are rare on average,
     but not excluded for a fixed N.
2. **The exceptional rate is limited by the large sieve.** The proved fraction is n²/H², so
   N^{−ε} costs H = N^{ε/2}, i.e. N^{5ε/6} in the cost. The Poisson heuristic, which the
   multimod and §8 data support, predicts e^{−λ} with λ ≍ H²/n. That would give N^{−ε} at
   polylog cost, but it is **conjectured**.
3. **Ratios q/p > 2** are not covered (multimod §3.9).
4. **Finite sizes.** At 40–80 bits HH is 2.3–4.9× cheaper in E (§8). The measured slope is
   smaller (0.203 vs 0.226), and an extrapolated crossover lies near 130 bits. That is **not a
   claim**: the slope difference partly reflects the H-schedule.
5. **Below 1/6.** multimod Thm A (enlarged sumset model) is unchanged. Nothing here goes
   below 1/6.

## 7. Review

- **GPT-6.1-sol referee dclt_1** (`codex/dclt_1.out`, read-only, own exact scripts).
  - Verdict: "the clean-band argument closes the previous order/extraction gap … Lemmas 1–10
    mathematically sound". It reproduced the checks independently: 131/131 clean non-FULL
    covering targets gave gcd = p, and 62/62 certified pairs succeeded.
  - Required fixes, all applied:
    - (1) "m ≤ m̄_k suffices" was false for cost. m = m̄_k or the endpoint-competitive rule is now
      mandatory (§2.2 step 3).
    - (2) Missing T′ in the giant-generation cost, plus the m-grid and small inversions, are now
      charged (§4 table).
    - (3) T(T+1) in E_a.
    - (4) Fast-arithmetic assumptions stated; E declared a proxy (§4).
    - (5) Float-free parameters, a GMP exponent product tree, and a finite-size guard. These are
      post-run hygiene; recorded runs are reproduced exactly (§9).
    - (6) McKee–Pinch §5 priority left unresolved (since resolved by r13, see below).
  - It also gave an exact example where the enlarged δ-range |δ| < m would be wrong for the
    *equivalence*. The note uses the enlargement only as a counting superset.
- **Independent Opus critic** (read-only; own implementation from the §2 spec; adversarial
  inputs p², primes, 3·x, 2^35+1, 2^36−1; 380 pairs including 180 with ord_q(2) ≤ 3000 to force
  FULL events).
  - Found: 0 lemma violations; 812 clean covering targets with 108 FULL, closed form = direct
    check every time; certified ⇒ success in 301/301.
  - It confirmed pre-registration compliance: freeze commit 13:45:02, runs from 13:45:07;
    post-run edits are provably inert.
  - **VERDICT: ACCEPT WITH FIXES.** All fixes were write-up level and are applied:
    - §0/§7/§8.3 filled; R-e decided by its frozen predictions;
    - HH exclusions 16 → 17;
    - E_block0 slopes added;
    - R-d range stated;
    - T′ bound written as Σ W_τ/m;
    - theorem_design H/Y conventions aligned;
    - n^5 → (2n)^{14/3} in E_c.
- **Fresh GPT-6.1-sol r13 referee** (notes/codex_r13_final_referee.md).
  - Verdict: ACCEPT WITH FIXES. Every proof step was CORRECT. The O(N^{−ε}) rate and the bit
    complexity were confirmed, with no hidden dependence on p or q.
  - Required fixes, all applied 2026-10-18:
    - (1) **Output validation.** The pre-fix code output 2 for N = 2, and the band prime ℓ
      for N = ℓ ∈ {7, 11, 19, 23, 31, 43, 47} (exhaustive 2 ≤ N ≤ 256). Now N < 4 returns
      FAIL; a band prime is output only if 1 < ℓ < N; every candidate passes `_proper`.
      There is a new exhaustive test over 2 ≤ N ≤ 10⁴ with two designs. The S_n runs are
      unaffected: 4 recorded C1 instances re-run with identical E_full and outcome.
    - (2) McKee–Pinch §5 comparison rewritten from the scans (§10).
    - (3) The §11 sieve-barrier sentence is now labelled a heuristic.
    - (4) Scope made explicit in §0/§1: n ≥ n₀(ε), bounded-ratio distinct-prime semiprimes,
      pair-counted exceptions.
    - (5) Carella 1308.2891 added, with the referee's assessment labelled as such.

## 8. Experiments (pre-registered; notes/dcltheorem_prereg.md, commit 3550e24)

The population is uniform pairs from S_n (`harness.sample_pair`). The algorithm sees N only.
Classes are assigned post hoc by the harness.

### 8.1 lemma (n = 28–40, 60 pairs each, design C1, all 4 blocks recorded)

- 240/240 correct and factored.
- **0 violations** of Lemmas 1–4 on every covering target.
- **0 violations** of the Lemma 5/7 equivalence (gcd = p ⇔ not FULL, exact δ-range) on 3002
  clean covering targets.
- 240/240 certified, and certified ⇒ success held in all 240.

### 8.2 oper (C1, C3: 32 pairs per n; C2: 8 pairs per n ≤ 64; all 4 blocks always run)

| n | C1 H/Y | C1 success (block 0 alone) | med E_first C1 | med E_full C1 | med E_first C3 | HH med E | C1/HH | list C1 block 0 | list HH | first-hit class C1 (P/Q/ACC) |
|---|---|---|---|---|---|---|---|---|---|---|
| 40 | 9/81 | 32/32 (32) | 1.03e5 | 1.07e6 | 1.27e5 | 2.12e4 | 4.87 | 679 | 252 | 16/10/6 |
| 48 | 9/206 | 32/32 (30) | 3.26e5 | 2.91e6 | 3.35e5 | 7.57e4 | 4.31 | 1650 | 742 | 14/15/3 |
| 56 | 10/503 | 32/32 (31) | 1.01e6 | 9.01e6 | 1.05e6 | 2.96e5 | 3.40 | 4266 | 2223 | 19/8/5 |
| 64 | 11/1461 | 32/32 (31) | 4.08e6 | 3.75e7 | 2.85e6 | 9.43e5 | 4.33 | 13223 | 6570 | 17/14/1 |
| 72 | 11/3097 | 32/32 (31) | 1.06e7 | 8.84e7 | 9.63e6 | 3.43e6 | 3.10 | 30610 | 19414 | 18/10/4 |
| 80 | 11/7804 | 32/32 (31) | 2.76e7 | 2.36e8 | 2.54e7 | 1.18e7 | 2.35 | 71004 | 58133 | 16/16/0 |

- C2 (heuristic full coverage, H = 19–27): 32/32 succeed. Median E_first is
  4.4e5 / 1.7e6 / 5.8e6 / 2.1e7 at n = 40 / 48 / 56 / 64.
- HH: 96/96 succeed. 17 runs are excluded from E as pre-registered (collision-only rule): 16
  exact-match early exits and 1 gcd_baby_order exit. *(Corrected from "16" after the Opus
  critic; the slope is 0.226 either way.)*

**Slopes of log₂(median) over n = 48–80:**

| quantity | C1 | C3 | C2 (48–64) | HH |
|---|---|---|---|---|
| E_first | **0.203** | 0.196 | 0.227 | **0.226** (E, collision runs) |
| E_block0 | 0.203 | 0.196 | 0.227 | — |
| E_full (4 blocks) | 0.200 | 0.202 | 0.237 | — |
| list elements (block 0 / HH list) | 0.171 | 0.163 | 0.189 | 0.196 |

**Decisions (pre-registered rules):**
- **R-a holds.** 0 incorrect outputs out of 656 DCL runs. 0 lemma violations. Every certified
  instance factored.
- **R-b: "DCL-det has a smaller measured exponent".** slope_C1 − slope_HH = −0.023 ≤ −0.01.
  The predicted ranges held: C1 0.203 ∈ [0.16, 0.21] and HH 0.226 ∈ [0.18, 0.23].
- **R-c holds.** C1/HH > 1 at every size (4.9 → 2.35). A naive extrapolation of the slope gap
  puts the crossover near 130 bits; this is not a claim.
- **R-d holds.** C2's slope is over 48–64 only (C2 was run only to 64 bits): 0.227. It is
  ≥ HH − 0.005 both against HH over 48–80 (0.226) and over 48–64 (0.227). So there is no
  advantage at heuristically full coverage.
- **R-e holds.** All three frozen predictions held (§8.3).

**Classes and null.**
- In 768 C1 blocks, 83 first hits were ACCIDENTAL, against Σ min(1, null) = 164.
  - The null is the expected number of coincidences, #giants·m·(1/ord_p γ + 1/ord_q γ).
  - Hits are therefore below the null bound, consistent with random gcd coincidences in a
    γ whose order mod p is reduced. Such hits are allowed: any gcd ∈ (1, N) is a factor.
- C3 has 137 accidental first hits against a null of 221; C2 has 24 against 107.
- STRUCTURAL_Q (q found through the swapped-role lock) is common. The theorem uses only
  STRUCTURAL_P.

**Interpretation.**
- The measured E slopes exceed the list slopes by ≈ 0.03 for both algorithms. That is the
  polylog (log d per poly coefficient, n-bit handles) drifting across 48–80 bits.
- C1's list slope (0.171) matches multimod's R1 model slope (0.171), now with the theorem's
  design, end to end.
- C1's E_full/E_first ≈ 8.5. The 4-band worst-case guarantee costs a constant ≈ 8–10 at these
  sizes: higher bands have more primes, hence more targets.

### 8.3 exc (C1 certificate components, 1000 pairs per n; no stage 3)

| n | certified | block 0 clean (pred. 1 − Σ_{P₀}1/(ℓ−1)) | block 0 uncovered (e^{−λ₀}) | block 0 certified | all 4 blocks uncovered | FULL / covering targets | median bits ord_q(2) / q |
|---|---|---|---|---|---|---|---|
| 48 | 1.000 | 0.939 (0.944) | 0.040 (0.049) | 0.900 | 0 | 5 / 11130 | 24 / 25 |
| 64 | 1.000 | 0.956 (0.954) | 0.021 (0.024) | 0.935 | 0 | 0 / 13321 | 32 / 33 |
| 80 | 1.000 | 0.976 (0.963) | 0.033 (0.046) | 0.943 | 0 | 0 / 11273 | 40 / 41 |
| 96 | 1.000 | 0.969 (0.969) | 0.028 (0.029) | 0.942 | 0 | 0 / ≈1.1·10⁴ | 48 / 49 |

**R-e (frozen predictions): all three hold.**
- Certified ≥ 0.99 at every n: observed 1.000.
- Block-0 clean within ±0.02 of 1 − Σ1/(ℓ−1): the largest deviation is +0.013 at 80 bits.
- FULL ≤ 1% of covering targets: observed ≤ 0.045%.

Uncovered fractions sit at or below the Poisson value e^{−λ₀}, as in multimod. ord_q(2) is
typically ≈ q/2, far above T. With four blocks, no sampled pair was uncovered in all blocks.
At these sizes the three exceptional sets E_a, E_b and E_c of §4(iii) are empirically
negligible (heuristic e^{−λ} regime; the proved rate needs H = n^{1+a}).

## 9. Reproduction

```
.venv/bin/python -m pytest -q tests/test_dcltheorem.py
PYTHONPATH=. .venv/bin/python -m src.dcltheorem.experiments config
PYTHONPATH=. .venv/bin/python -m src.dcltheorem.experiments lemma
PYTHONPATH=. .venv/bin/python -m src.dcltheorem.experiments hh
PYTHONPATH=. .venv/bin/python -m src.dcltheorem.experiments oper
PYTHONPATH=. .venv/bin/python -m src.dcltheorem.experiments exc
PYTHONPATH=. .venv/bin/python -m src.dcltheorem.experiments summarize
```

**Post-run code changes** (referee dclt_1 hygiene). These are:
- an integer-only m-grid, `(119*m)//100+1` in place of `int(m*1.19)+1`;
- an integer cube root in `theorem_design`;
- GMP products in the exponent tree;
- a finite-size guard ℓ | AB.

None changes any recorded run. Re-running 12 recorded C1/C3 instances (n = 40, 56; seeds
100–102) reproduced E_full and the first-success block exactly.

## 10. Prior art (GPT-6.1-sol r12 + dclt_1 + own WebSearch; curated)

- **Ingredients (all prior art).**
  - Weighted Fermat / Lehman searches: Lawrence 1895, Lehman 1974.
  - Modulus-squared lock: McKee 1999 (heuristic N^{1/4}); HH 2022 Prop. 3.3.
  - Shared one-sided collision search with multipoint evaluation: Harvey 2021,
    Alg. 4.1–4.2.
  - Arithmetic large sieve: Montgomery 1968; Montgomery–Vaughan 1973.
  - Elementary counting of primes with small ord_q(2): divisors of 2^k − 1.
- **Exponent below 1/5 for almost all bounded-ratio prime pairs (unconditional, deterministic):
  no prior result found** by GPT-6.1-sol (two independent sweeps) or by our searches.
  - The known bounds are 1/5 (Harvey 2021; HH 2022; Gao et al. 2512.19076 for balanced N,
    log factors only).
  - The conditional route is 1/6 under Umans–Wang's conjecture; He–Sahai 2026 refute its
    rank-1 form.
  - Almost-all unrestricted semiprimes are easy by trial division (Crişan–Erban counting), but
    that does not apply to S_n.
- **The clean-band trick** (disjoint prime bands so that one is coprime to q − 1, plus large
  ord_q(2) for almost all q, plus a direct count of two-sided true giants) was not found in
  the factoring papers checked.
- **McKee–Pinch 1996 (ANTS II, LNCS 1122) §5, pp. 221–223: a direct precursor.** The r13
  referee obtained it, and we read the scans ourselves: references/mckee_pinch_1996/
  mckee-03.jpg and mckee-04.jpg; source https://cr.yp.to/bib/1996/mckee-deterministic.html.
  - **Lehmer's speed-up.** If p ≡ q ≡ 1 (mod k), Fermat runs in O(q/k²).
  - **The authors' extension.** For p ≡ ±q (mod k), Fermat runs in O(q/k^{2−ε}), via
    n = u² − v² with u or v determined mod k² from square roots mod k. That is the
    modulus-squared lock.
  - **Multipliers.** If λp ≡ ±μq (mod k), the same method splits λμn in
    O(max(λp, μq)/k²·n^ε).
  - **A deterministic O(n^{1/3+ε}) algorithm.** Take k = 2⌊n^{1/3}⌋ + 1, and a pigeonhole
    lemma that guarantees λμ ≤ k with λp ≡ ±μq (mod k).
  - **A heuristic balanced method, "not however an algorithm".** Take primes k ∈ (n^{1/4},
    2n^{1/4}) and small λ, μ, and assume λp − μq ≡ 0 (mod k) for some of them. Each (λμ, k)
    costs O(n^ε). The expected running time O(n^{1/4+ε}) rests on the heuristic that q − p
    has a prime factor in the band with probability ≫ 1/log n.
  - **So the core event** (a band prime divides a small-weight combination Bp − Aq, with
    N-determined residues and a modulus-squared reduction) **is theirs.**
  - **What D′ adds:**
    - (i) the band at N^{1/6}H^{−1/3} instead of N^{1/4}. This is made possible by
    - (ii) one baby list shared by all band primes via the cofactor-powered handles
      u_ℓ^{B+AN−c}, so each target costs one giant step against a window of size ≈ N^{1/6};
    - (iii) a proof that the event happens for all but O(N^{−ε}) of the pairs: the
      prime-pair large sieve with growing weight height H, replacing their probability
      heuristic;
    - (iv) one-sided extraction with no order hypothesis: the clean band, ord_q(2), and the
      Lemma 10 count;
    - (v) a full bit-complexity accounting.
  - Their per-k cost is O(n^ε) because q/k² = O(1). Ours trades a larger window for fewer,
    smaller primes through the shared baby list.
- **Carella, arXiv 1308.2891.** It explicitly claims a deterministic exponent-1/6 factoring
  algorithm.
  - *The r13 referee's assessment, not independently re-derived by us:* the displayed proof
    loses a necessary primitive-polynomial normalisation. Every coefficient of the correctly
    expanded polynomial is divisible by B. Removing that content lowers the scaled height
    from ≈ N to N^{5/6}, and the cited small-root condition then no longer supports
    X = Y = N^{1/3}.
  - The referee does not regard that proof as establishing the bound. We list it as a claimed
    but (per that assessment) unproved prior statement of the same exponent.
- **Claim limited to:** "no prior formulation found by us" of the *combination*: (i)–(v) above,
  on top of McKee–Pinch §5's divisibility event and HH 2022 Prop. 3.3's lock. Novelty
  judgments are left to the orchestrator.

## 11. Candidate card

- **Mechanism.** DCL with one-sided extraction certified on a clean band (§2).
  **Representation.** (Z/N)^*, base 2.
- **Factor feature.** g = (Bp + Aq − c)/ℓ² in a whole-range window; the true giant equals
  γ^{i*} mod p exactly.
- **Outcome.**
  - Proved: Theorem D′ and Corollary D″.
  - Experimentally supported: correctness, every lemma, measured E slope 0.203 vs HH 0.226 at
    48–80 bits, with HH 2.3–4.9× cheaper in absolute terms.
  - Refuted for this model: none.
  - Conjectured: e^{−λ} exceptional rate (Poisson), DCH (worst case).
- **Scoped obstruction.**
  - The large-sieve exceptional rate costs N^{5ε/6} for N^{−ε}.
  - Worst case needs DCH plus large order mod the specific prime q on a clean band.
  - Below 1/6 is blocked in the sumset model (multimod Thm A).
- **Most promising open question.** Improve the exceptional-rate/cost trade-off beyond the
  large sieve's n²/H².
  - **Heuristic (not a theorem; qualified after r13).** For a fixed q, the p-interval has
    length P ≈ Y³, so a sieve whose level is limited to ≈ P only sees moduli built from
    ≤ 3 band primes. We *expect* that standard one-variable upper-bound sieves restricted to
    such levels (Selberg, the large sieve) give only polynomial-in-1/λ savings, and this
    explains the loss in our proof. We have **not** proved any barrier: other sieve weights,
    bilinear or dispersion methods, or the averaging over q below might do better.
  - Averaging over q gives a two-variable sieve over the box (classes Bp ≡ Aq mod ℓ, density
    |V|/ℓ) with level up to P². That allows products of up to ~6 band primes and plausibly a
    polynomial gain O((n/H²)^{c}) with c > 1 (**conjectured, not yet tested**).
  - A Poisson-type e^{−cλ} bound would give N^{1/6}·polylog cost with N^{−ε} exceptions
    (**conjectured**).
